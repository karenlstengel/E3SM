module kessler_perf_log

   !-----------------------------------------------------------------------
   ! Fortran counterpart to kessler_perf_log.py (this same directory) --
   ! same design, same CSV schema, so Fortran (CPU and GPU/OpenACC) timings
   ! can be compared directly against the JAX-side summary. Intended to be
   ! copied alongside the Fortran kessler.F90 sources at:
   !   atmospheric_physics/schemes/kessler/
   !   GPU_ports/atmospheric_physics/schemes/kessler/
   !
   ! Set KESSLER_PERF_LOG_PATH to choose the CSV path; falls back to
   ! ./kessler_perf_log.csv in the run directory otherwise.
   !
   ! log_call() only ever touches in-process memory (a small table of
   ! running per-label totals) -- there is no file I/O per call. An earlier
   ! version opened the CSV and took a cross-process flock() on every call,
   ! which serializes every MPI rank against every other rank on every
   ! timestep and measurably slows runs down (this is exactly what was
   ! fixed on the Python side; mirrored here even though this module isn't
   ! wired into any kessler_run yet, so it doesn't get reintroduced later).
   !
   ! flush_log() does the file I/O instead, and is meant to be called
   ! exactly once per rank, at simulation finalize -- NOT YET WIRED UP;
   ! whoever wires this module into a kessler_run should also add a call
   ! to flush_log() at that scheme's finalize/timestep_final entry point,
   ! mirroring kessler.py's finalize() -> KesslerMicrophysics::finalize_impl().
   ! If a rank never reaches finalize (crash, MPI_Abort), that rank's
   ! timings are simply lost -- an accepted tradeoff for a dev/perf tool,
   ! not a science output.
   !
   ! flush_log() appends this rank's per-label totals to a small raw
   ! intermediate CSV (KESSLER_PERF_LOG_PATH with a ".raw" suffix) under
   ! one flock(), then recomputes the GPTL-style summary (count/walltotal/
   ! wallmax/wallmin per label, across all ranks that have flushed so far)
   ! and rewrites KESSLER_PERF_LOG_PATH itself with that summary -- so that
   ! path always holds the final summary directly, no separate
   ! post-processing step needed. Because every rank's append +
   ! summary-rewrite happens inside one flock() critical section, whichever
   ! rank's critical section runs last (in real time) is guaranteed to see
   ! every other rank's already-committed row, so the summary is exactly
   ! correct once every rank has finalized, regardless of what order ranks
   ! arrive in. Raw POSIX open()/write() (via ISO_C_BINDING) with O_APPEND
   ! is used for that append, rather than Fortran's own OPEN(POSITION=
   ! 'APPEND'): most Fortran runtimes implement POSITION='APPEND' as a
   ! single seek-to-end done at OPEN time, not as the kernel-level O_APPEND
   ! flag -- so two ranks opening around the same time can seek to the same
   ! offset and one write can clobber the other. A real O_APPEND file
   ! descriptor makes every write() atomically land at the current end of
   ! file regardless of that race. The summary is then read back and
   ! rewritten using plain Fortran formatted I/O on separate units -- safe
   ! to do while still holding the flock, since flock() only serializes
   ! against other flock() callers, not against a plain read() from the
   ! same process.
   !-----------------------------------------------------------------------

   use, intrinsic :: iso_c_binding
   use ccpp_kinds, only: kind_phys

   implicit none
   private

   public :: log_call, flush_log

   integer, parameter :: PATH_LEN  = 512
   integer, parameter :: LABEL_LEN = 128

   character(len=PATH_LEN), save :: log_path = ''
   logical,                 save :: path_resolved = .false.
   character(len=PATH_LEN), save :: raw_path = ''
   logical,                 save :: raw_path_resolved = .false.
   integer,                 save :: rank = 0
   logical,                 save :: rank_resolved = .false.

   ! In-memory per-label running totals (this rank only). Mutated by
   ! log_call(), read/reset only by flush_log().
   integer, parameter :: MAX_LABELS = 16
   character(len=LABEL_LEN), save :: tot_label(MAX_LABELS)     = ''
   integer,                  save :: tot_calls(MAX_LABELS)     = 0
   real(kind_phys),          save :: tot_walltotal(MAX_LABELS) = 0.0_kind_phys
   real(kind_phys),          save :: tot_callmax(MAX_LABELS)   = 0.0_kind_phys
   real(kind_phys),          save :: tot_callmin(MAX_LABELS)   = 0.0_kind_phys
   integer,                  save :: n_labels = 0

   ! Linux x86-64 syscall constants (see fcntl.h / sys/file.h). This
   ! module is POSIX/Linux-only, same as the Python fcntl.flock() side.
   integer(c_int), parameter :: O_WRONLY  = int(o'0000001', c_int)
   integer(c_int), parameter :: O_CREAT   = int(o'0000100', c_int)
   integer(c_int), parameter :: O_APPEND  = int(o'0002000', c_int)
   integer(c_int), parameter :: MODE_0644 = int(o'0000644', c_int)
   integer(c_int), parameter :: LOCK_EX   = 2_c_int
   integer(c_int), parameter :: LOCK_UN   = 8_c_int
   integer(c_int), parameter :: SEEK_END  = 2_c_int

   interface
      function c_open(path, flags, mode) bind(C, name="open")
         import :: c_int, c_char
         character(kind=c_char), intent(in) :: path(*)
         integer(c_int), value :: flags, mode
         integer(c_int) :: c_open
      end function c_open

      function c_close(fd) bind(C, name="close")
         import :: c_int
         integer(c_int), value :: fd
         integer(c_int) :: c_close
      end function c_close

      function c_flock(fd, operation) bind(C, name="flock")
         import :: c_int
         integer(c_int), value :: fd, operation
         integer(c_int) :: c_flock
      end function c_flock

      function c_lseek(fd, offset, whence) bind(C, name="lseek")
         import :: c_int, c_long
         integer(c_int), value :: fd
         integer(c_long), value :: offset
         integer(c_int), value :: whence
         integer(c_long) :: c_lseek
      end function c_lseek

      function c_write(fd, buf, count) bind(C, name="write")
         import :: c_int, c_char, c_long
         integer(c_int), value :: fd
         character(kind=c_char), intent(in) :: buf(*)
         integer(c_long), value :: count
         integer(c_long) :: c_write
      end function c_write
   end interface

contains

   subroutine resolve_log_path()
      integer :: vlen, stat
      character(len=PATH_LEN) :: val

      if (path_resolved) return

      call get_environment_variable('KESSLER_PERF_LOG_PATH', value=val, &
           length=vlen, status=stat)
      if (stat == 0 .and. vlen > 0) then
         log_path = val(1:vlen)
      else
         log_path = 'kessler_perf_log.csv'
      end if
      path_resolved = .true.
   end subroutine resolve_log_path

   subroutine resolve_raw_path()
      integer :: dot_pos

      if (raw_path_resolved) return
      call resolve_log_path()

      dot_pos = index(trim(log_path), '.', back=.true.)
      if (dot_pos > 0) then
         raw_path = log_path(1:dot_pos-1) // '.raw' // trim(log_path(dot_pos:))
      else
         raw_path = trim(log_path) // '.raw'
      end if
      raw_path_resolved = .true.
   end subroutine resolve_raw_path

   subroutine resolve_rank()
      integer, parameter :: n_vars = 5
      character(len=32) :: env_vars(n_vars)
      character(len=32) :: val
      integer :: i, vlen, stat, ios, r

      if (rank_resolved) return

      ! Common MPI launchers: OpenMPI, MPICH/Intel MPI, Slurm srun.
      env_vars(1) = 'OMPI_COMM_WORLD_RANK'
      env_vars(2) = 'PMI_RANK'
      env_vars(3) = 'SLURM_PROCID'
      env_vars(4) = 'MPI_LOCALRANKID'
      env_vars(5) = 'PALS_RANKID'

      rank = 0
      do i = 1, n_vars
         call get_environment_variable(trim(env_vars(i)), value=val, &
              length=vlen, status=stat)
         if (stat == 0 .and. vlen > 0) then
            read(val(1:vlen), *, iostat=ios) r
            if (ios == 0) then
               rank = r
               exit
            end if
         end if
      end do
      rank_resolved = .true.
   end subroutine resolve_rank

   ! Find full_label in the running-totals table, adding a new entry if
   ! it isn't there yet (silently dropping the sample if the table is
   ! already full -- best-effort, same spirit as the rest of this module).
   function find_or_add_label(full_label) result(idx)
      character(len=*), intent(in) :: full_label
      integer :: idx, i

      do i = 1, n_labels
         if (trim(tot_label(i)) == trim(full_label)) then
            idx = i
            return
         end if
      end do

      if (n_labels < MAX_LABELS) then
         n_labels = n_labels + 1
         idx = n_labels
         tot_label(idx)     = full_label
         tot_calls(idx)     = 0
         tot_walltotal(idx) = 0.0_kind_phys
         tot_callmax(idx)   = -huge(1.0_kind_phys)
         tot_callmin(idx)   = huge(1.0_kind_phys)
      else
         idx = -1
      end if
   end function find_or_add_label

   subroutine log_call(label, ncol, nz, dt, elapsed)
      character(len=*), intent(in) :: label
      integer,          intent(in) :: ncol, nz
      real(kind_phys),  intent(in) :: dt, elapsed

      character(len=LABEL_LEN) :: full_label
      integer :: idx

      full_label = 'a:EAMxx::kessler::run::F90_run::' // trim(label)

      !$omp critical (kessler_perf_log_totals)
      idx = find_or_add_label(full_label)
      if (idx > 0) then
         tot_calls(idx)     = tot_calls(idx) + 1
         tot_walltotal(idx) = tot_walltotal(idx) + elapsed
         if (elapsed > tot_callmax(idx)) tot_callmax(idx) = elapsed
         if (elapsed < tot_callmin(idx)) tot_callmin(idx) = elapsed
      end if
      !$omp end critical (kessler_perf_log_totals)
   end subroutine log_call

   subroutine flush_log()
      integer, parameter :: MAX_SUM = 32

      ! Snapshot of this rank's totals, copied out under the same critical
      ! section log_call() uses, mirroring the Python side's
      ! `with _lock: totals = dict(...)`.
      character(len=LABEL_LEN) :: my_label(MAX_LABELS)
      integer                  :: my_calls(MAX_LABELS)
      real(kind_phys)          :: my_walltotal(MAX_LABELS), my_callmax(MAX_LABELS), my_callmin(MAX_LABELS)
      integer                  :: my_n

      ! Aggregate across every rank's rows seen in the raw CSV so far.
      character(len=LABEL_LEN) :: sum_label(MAX_SUM)
      integer                  :: sum_count(MAX_SUM)
      real(kind_phys)          :: sum_walltotal(MAX_SUM)
      real(kind_phys)          :: sum_wallmax(MAX_SUM)
      integer                  :: sum_wallmax_rank(MAX_SUM)
      real(kind_phys)          :: sum_wallmin(MAX_SUM)
      integer                  :: sum_wallmin_rank(MAX_SUM)
      real(kind_phys)          :: sum_callmax(MAX_SUM)
      integer                  :: sum_callmax_rank(MAX_SUM)
      integer                  :: n_sum

      character(len=PATH_LEN+1) :: cpath
      character(len=400) :: row
      character(len=64)  :: header
      integer(c_int)  :: fd, rc
      integer(c_long) :: fsize, nbytes, nwritten

      character(len=512) :: line
      character(len=LABEL_LEN) :: row_label
      integer :: row_rank, row_count
      real(kind_phys) :: row_walltotal, row_callmax, row_callmin
      integer :: read_unit, sum_unit, iostat_val, ios
      integer :: i, j
      logical :: found

      call resolve_log_path()
      call resolve_raw_path()
      call resolve_rank()

      !$omp critical (kessler_perf_log_totals)
      my_n = n_labels
      do i = 1, my_n
         my_label(i)     = tot_label(i)
         my_calls(i)     = tot_calls(i)
         my_walltotal(i) = tot_walltotal(i)
         my_callmax(i)   = tot_callmax(i)
         my_callmin(i)   = tot_callmin(i)
      end do
      !$omp end critical (kessler_perf_log_totals)

      if (my_n == 0) return  ! nothing logged on this rank -- nothing to flush

      cpath = trim(raw_path) // c_null_char
      fd = c_open(cpath, ior(ior(O_WRONLY, O_CREAT), O_APPEND), MODE_0644)
      if (fd < 0) return  ! best-effort: never abort the model run over a logging failure

      rc = c_flock(fd, LOCK_EX)  ! one-time cost: only taken once per rank, at finalize

      fsize = c_lseek(fd, 0_c_long, SEEK_END)
      if (fsize <= 0) then
         header = 'label,rank,count,walltotal_s,callmax_s,callmin_s' // char(10)
         nbytes = len_trim(header)
         nwritten = c_write(fd, header, nbytes)
      end if

      do i = 1, my_n
         write(row, '(A,",",I0,",",I0,",",G0,",",G0,",",G0,A)') &
              trim(my_label(i)), rank, my_calls(i), my_walltotal(i), my_callmax(i), my_callmin(i), char(10)
         nbytes = len_trim(row)
         nwritten = c_write(fd, row, nbytes)
      end do

      ! Read the whole raw file back (a plain Fortran unit, separate from
      ! the C-level fd above) and recompute the aggregate summary. Safe to
      ! do while still holding the flock: flock() only serializes against
      ! other flock() callers, and the write()s above already landed in
      ! the OS page cache, so this read sees them immediately.
      n_sum = 0
      open(newunit=read_unit, file=trim(raw_path), status='old', action='read', &
           form='formatted', iostat=iostat_val)
      if (iostat_val == 0) then
         read(read_unit, '(A)', iostat=ios) line  ! header row, discarded
         do
            read(read_unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) == 0) cycle

            read(line, *, iostat=ios) row_label, row_rank, row_count, row_walltotal, row_callmax, row_callmin
            if (ios /= 0) cycle

            found = .false.
            do j = 1, n_sum
               if (trim(sum_label(j)) == trim(row_label)) then
                  found = .true.
                  exit
               end if
            end do
            if (.not. found) then
               if (n_sum >= MAX_SUM) cycle
               n_sum = n_sum + 1
               j = n_sum
               sum_label(j)        = row_label
               sum_count(j)        = 0
               sum_walltotal(j)    = 0.0_kind_phys
               sum_wallmax(j)      = -huge(1.0_kind_phys)
               sum_wallmax_rank(j) = -1
               sum_wallmin(j)      = huge(1.0_kind_phys)
               sum_wallmin_rank(j) = -1
               sum_callmax(j)      = -huge(1.0_kind_phys)
               sum_callmax_rank(j) = -1
            end if

            sum_count(j)     = sum_count(j) + row_count
            sum_walltotal(j) = sum_walltotal(j) + row_walltotal
            if (row_walltotal > sum_wallmax(j)) then
               sum_wallmax(j) = row_walltotal; sum_wallmax_rank(j) = row_rank
            end if
            if (row_walltotal < sum_wallmin(j)) then
               sum_wallmin(j) = row_walltotal; sum_wallmin_rank(j) = row_rank
            end if
            if (row_callmax > sum_callmax(j)) then
               sum_callmax(j) = row_callmax; sum_callmax_rank(j) = row_rank
            end if
         end do
         close(read_unit)
      end if

      ! Overwritten in full each time, same as the Python side, so
      ! log_path always holds the current summary across every rank that
      ! has flushed so far.
      open(newunit=sum_unit, file=trim(log_path), status='replace', action='write', &
           form='formatted', iostat=iostat_val)
      if (iostat_val == 0) then
         write(sum_unit, '(A)') 'label,count,walltotal_s,wallmax_s,wallmax_rank,wallmin_s,wallmin_rank,callmax_s,callmax_rank'
         do j = 1, n_sum
            write(sum_unit, '(A,",",I0,",",G0,",",G0,",",I0,",",G0,",",I0,",",G0,",",I0)') &
                 trim(sum_label(j)), sum_count(j), sum_walltotal(j), &
                 sum_wallmax(j), sum_wallmax_rank(j), sum_wallmin(j), sum_wallmin_rank(j), &
                 sum_callmax(j), sum_callmax_rank(j)
         end do
         close(sum_unit)
      end if

      rc = c_flock(fd, LOCK_UN)
      rc = c_close(fd)

   end subroutine flush_log

end module kessler_perf_log

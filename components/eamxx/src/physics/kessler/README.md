# EAMxx Kessler Bridge 

## Background: 

See `components/eam/src/physics/crm/pam/external/physics/micro/kessler/Microphysics.h` for background information on EAM implementation. 

Non-CAM Fortran version of code (with C bindings) found in: `components/eam/src/physics/crm/pam/external/physics/micro/kessler/kessler.f90`

CAM Fortran routines are in: [this repo](https://github.com/ESCOMP/atmospheric_physics)

We are using the Kessler provided in: [atmospheric_physics](git@github.com:ESCOMP/atmospheric_physics.git). 

## Design & Code Layout

The Kessler suite is in `atmospheric_physics/schemes/kessler`. 

List of parameters: 

Kessler.F90:
kessler_init(real lv_in, real pref_in, real rhoqr_in, char* errmsg, int errflg) {r,r,r,w,w}  
kessler_run(int ncol, int nz, real dt, int lyr_surf, int lyr_toa, real cpair, real rair, real rho, real z, &
        real pk, real theta, real qv, real qc, real qr, real precl, real relhum, char* scheme_name, char* errmsg, int errflg) {r,r,r,r,r,r,r,r,r,r,u,u,u,u,w,w,w,w,w}  

Kessler_update.F90:

kessler_update_init(real gravit_in, char* errmsg, int errflg) {r, w, w}  
kessler_update_timestep_init(real temp, real temp_prev, real ttend_t, char* errmsg, int errflg) {r,w,w,w,w}  
kessler_update_run(int nz, int ncol, real dt, real theta, real exner, real temp_prev, real ttend_t, char* errmsg, int errflg) {r,r,r,r,r,r,u,w,w}  
kessler_update_timestep_final(int nz, real cpair, real temp, real zm, real phis, real st_energy, char* errmsg, int errflg) {r,r,r,r,r,w,w,w}  


## GPU stuff

### Implementation Notes

1. For the Fortran link to work, you must add both the C++ side and Fortran side field variables to the `ATMBufferManager` with `Packagename::init_buffers(const ATMBufferManager &buffer_manager)` and `Packagename::requested_buffer_size_in_bytes()` (see the main interface C++ file).
2. It seems like a good idea to create structs for passing input fields, output fields, and (if needed) parameters/etc. Note that this requires each member of the struct that is a field to be defined for both the C++ side (usually with `Spack` and Kokkos views) and the Fortran side (_unmanaged_ view of type `Real` in the same dimensions as the C++ version UNLESS WE ARE USING OPENACC THEN WE USE MANAGED VIEWS WITH `Real`???). 
3. You should create a transpose function for converting between the Fortran and C++ fields. See `transpose()` in the `packagename_functions.hpp` file. 
4. I _think_ you can only pass to Fortran the `get_field_out("fieldname").get_view<Spack**>().data()` object into the C to Fortran binding. 

### GPU vs. CPU performance 

## Misc notes

Example OpenACC Fortran loop. Recall that it is more efficient to loop over `k` first in fortran 

```fortran
    !$acc parallel deviceptr(p_mid) reduction(max:p_mid_max)
    !$acc loop gang vector collapse(2)
    do k = 1,pver
        do i = 1,ncol
            p_mid(i,k) = p_mid(i,k) * 2.0
        end do
    end do
    p_mid_max = MAXVAL(p_mid)
    !$acc end parallel
```
# EAMxx Kessler Bridge 

## Background: 

See `components/eam/src/physics/crm/pam/external/physics/micro/kessler/Microphysics.h` for background information on EAM implementation. 

Non-CAM Fortran version of code (with C bindings) found in: `components/eam/src/physics/crm/pam/external/physics/micro/kessler/kessler.f90`

CAM Fortran routines are in: [this repo](https://github.com/ESCOMP/atmospheric_physics)

We are using the Kessler provided in: [atmospheric_physics](git@github.com:ESCOMP/atmospheric_physics.git). 

Note that additional initial conditions files might need to be downloaded from https://web.lcrc.anl.gov/public/e3sm/inputdata/atm/scream/init/. 

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
2. It seems like a good idea to create structs for passing input fields, output fields, and (if needed) parameters/etc. Note that this requires each member of the struct that is a field to be defined for both the C++ side (usually with `Pack` and Kokkos views) and the Fortran side (_unmanaged_ view of type `Real` in the same dimensions as the C++ version UNLESS WE ARE USING OPENACC THEN WE USE MANAGED VIEWS WITH `Real`???). 
3. You should create a transpose function for converting between the Fortran and C++ fields. See `transpose()` in the `packagename_functions.hpp` file. 
4. I _think_ you can only pass to Fortran the `get_field_out("fieldname").get_view<Pack**>().data()` object into the C to Fortran binding. 

### GPU vs. CPU performance 

## Misc notes

Example OpenACC Fortran loop. Recall that it is more efficient to loop over `k` first in fortran 

```fortran
    !$acc parallel deviceptr(p_mid) 
    !$acc loop gang vector collapse(2) reduction(max:p_mid_max)
    do k = 1,pver
        do i = 1,ncol
            p_mid(i,k) = p_mid(i,k) * 2.0
        end do
    end do
    p_mid_max = MAXVAL(p_mid)
    !$acc end parallel
```


# Analytic conditions & test setup

StormSPEED uses https://github.com/ESCOMP/CAM/blob/cam_development/src/dynamics/tests/initial_conditions/ic_baroclinic.F90 c/o Jesse N. This sets: 

requires: vcoord,latvals, lonvals, z_int, 
- U: wind at all levels, u
- V: wind at all levels, v
- T: temperature at all levels (T_mid?)
- PS: surface pressure
- PHIS: phis (surface geopotential)
- Q: [qv, qr, qc, qi] need to check on order here
- can use to save z_mid and p_mid too I think

The other option is: https://github.com/ESCOMP/CAM/tree/cam_development/src/dynamics/tests c/o Jesse N.

Need to implement the baroclinic conditions in EAMxx. Question - what do we do about the aquaplanet? Do we use aquaplanet in StormSPEED? 
Below I added how to change to 58 levels to match StormSPEED but we should be able to swap the vertical levels files for 72 or 58 between the two models directly. 

## TO DO

1. In the namelists.xml we need to add something like `initial_conditions_analytic` to use instead of `initial_conditions_filename`. should contain a string to specify name of analytical function we are using ex:
    ```xml
        <initial_conditions>
            <analytic type="array(string)">UNSET</analytic>
            <analytic>"ic_baroclinic"</analytic>
        </initial_conditions>
    ```
For Kessler, we probably want to make a COMPSET with `nlev=58`. Also need to set the vertical coordinates to match StormSPEED with: 
    ```xml
    <!-- Grids manager specs -->
    <grids_manager>
        <type>homme</type>
        <physics_grid_type>gll</physics_grid_type>
        <physics_grid_type hgrid=".*pg2">pg2</physics_grid_type>
        <physics_grid_rebalance>none</physics_grid_rebalance>
        <dynamics_namelist_file_name>./data/namelist.nl</dynamics_namelist_file_name>
        <vertical_coordinate_filename type="file">UNSET</vertical_coordinate_filename>
        <vertical_coordinate_filename nlev="58">/glade/campaign/cesm/cesmdata/inputdata/atm/cam/inic/cam_vcoords_L58_c250227.nc</vertical_coordinate_filename>
        ...
    </grids_manager>
    ```
and in case script we probably want to have `./create_newcase --case ${CASE_NAME} .... --user-mods-dir ${CCSMROOT}/components/eamxx//cime_config/testdefs/testmods_dirs/eamxx/output/preset/2 ${CCSMROOT}/components/eamxx//cime_config/testdefs/testmods_dirs/eamxx/L58`

where `${CCSMROOT}/components/eamxx//cime_config/testdefs/testmods_dirs/eamxx/L58/shell_comands` contains 
```./xmlchange SCREAM_CMAKE_OPTIONS="`./xmlquery -value SCREAM_CMAKE_OPTIONS | sed 's/SCREAM_NUM_VERTICAL_LEV [0-9][0-9]*/SCREAM_NUM_VERTICAL_LEV 72/'`"```

2. in `/glade/derecho/scratch/kstengel/E3SM/E3SM/components/eamxx/src/control/atmosphere_driver.cpp::create_grids()` line 281 add something like:
    ```cpp
    else if (ic_pl.isParameter("analytic")) {
        // Initial run, if ICs are an analytic function, pass the name.
        auto ic_analytic = ic_pl.get<std::string>("analytic");
        gm_params.set("ic_analytic", ic_analytic);
        m_atm_params.sublist("provenance").set("initial_conditions_analytic",ic_analytic);
    }
    ```
    and line 1255 something like:
    ```cpp
    // If analytic is specified, compute initial conditions on all grids
    if (ic_pl.isParameter("analytic")) {
        // Now loop over all grids, and load from file the needed fields on each grid (if any).
        const auto& ic_analytic = ic_pl.get<std::string>("analytic");
        m_atm_logger->info("    [EAMxx] IC analytic type: " + ic_analytic);

        for (const auto& it : m_grids_manager->get_repo()) {
        const auto& grid = it.second;
        const auto& grid_name = grid->name();

        if (ic_fields_names[grid_name].size()==0) continue;

        std::vector<Field> ic_fields;
        for (const auto& fn : ic_fields_names[grid_name]) {
            ic_fields.push_back(m_field_mgr->get_field(fn,grid_name));
        }
        
        // Call the analytic function. I think this should go to a switch statement to call the correct analytic function but should check with EAMxx team about how they would want this setup
        set_fields_analytic(ic_fields,grid,ic_analytic);

        // I don't think we need to worry about IOP enabled options since we should just do the computation at the lat/lon column in the grid
        
        for (auto& f : ic_fields) {
            f.get_header().get_tracking().update_time_stamp(m_current_ts);
        }
        }
    }
    ```

3. compute the analytical initial conditions. see `/glade/derecho/scratch/kstengel/E3SM/E3SM/components/eamxx/src/share/io/scorpio_input.cpp`. 

## Field perturbation

### StormSPEED: 
???

### EAMxx:

Should be able to do automatically with the `perturbed_fields` and corresponding fields in the namelists.xml options. 
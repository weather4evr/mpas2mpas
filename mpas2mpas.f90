program mpas2mpas

use netcdf

use kinds, only : i_kind, r_double, r_kind
use namelist_mod, only : read_namelist, source_grid_template_file, destination_grid_template_file, &
                         file_to_interpolate, output_file, weight_file, &
                         netcdf_output_type, print_before_and_after_interp_values, &
                         variables_to_copy_from_destination_file, do_vertical_interpolation
use mpas_netcdf_interface, only : mpas_grid, init_mpas_grid, deallocate_mpas_grid, &
                                  open_netcdf, close_netcdf, get_netcdf_var, &
                                  get_netcdf_dims, get_netcdf_info, define_output_file_from_template, & 
                                  output_variable, source_grid, &
                                  destination_grid, sanity_checks
use interpolation_stuff, only : interpolate_mpas, nearest_cell, project_to_edges, calc_horiz_interp_weights, &
                                check_actual_type, print_min_max, calc_vert_interp_weights, &
                                vert_interp_weights_edges, vert_interp_weights, vert_interp_levels_edges, vert_interp_levels

implicit none

character(len=500) :: varname, data_char
character(len=NF90_MAX_NAME) :: dim_names(4)
integer(i_kind) :: i, j, k, v
integer(i_kind) :: ncfileid, ncfileid_dest, ncidout, rcode, ncvarid
integer(i_kind) :: xtype, ndims_var, natts, dimids(10), dims(4), char_len
integer(i_kind) :: start_cell, nz, ncells, ndims, nvars
integer(i_kind) :: iunit_weight = 42
integer(i_kind), allocatable, dimension(:)   :: nearest_cell_to_output_point
integer(i_kind), allocatable, dimension(:,:) :: interp_cells
real(r_double), allocatable, dimension(:,:)  :: interp_weights
!real(r_kind), allocatable, dimension(:,:,:)  :: vert_interp_weights_edges, vert_interp_weights
!integer(i_kind), allocatable, dimension(:,:,:) :: vert_interp_levels_edges, vert_interp_levels

real(r_kind), allocatable, dimension(:,:)  :: input_data_real, output_data_real
real(r_kind), allocatable, dimension(:,:)  :: input_u, input_v, output_u, output_v, edge_normal_wind, u
real(r_kind), allocatable, dimension(:,:)  :: zgrid_interp
!real(r_kind), allocatable, dimension(:,:)  :: input_pp, input_pb
integer(i_kind), allocatable, dimension(:,:) :: input_data_int, output_data_int
logical :: fexist, is_integer, do_winds, mismatch

!type(mpas_grid) :: source_grid, destination_grid

!!!!!

! Read the namelist
call read_namelist

! Read a template file on the source grid and get its dimensions and grid info
!source_grid accessible everywhere from mpas_netcdf_interface
call init_mpas_grid(source_grid_template_file,source_grid)

! Read a template file on the destination/target grid and get its dimensions and grid info
!destination_grid accessible everywhere from mpas_netcdf_interface
call init_mpas_grid(destination_grid_template_file,destination_grid)

! Do some sanity checks to make sure number of vertical levels and soil levels match
call sanity_checks(source_grid, destination_grid, mismatch) ! output is mismatch
if ( mismatch ) then
   call cleanup
   stop
endif

! Find the nearest MPAS cells to each point on the output grid
start_cell = 1
allocate(nearest_cell_to_output_point(destination_grid%ncells))
do j=1,destination_grid%ncells
   nearest_cell_to_output_point(j) = nearest_cell( destination_grid%latCell(j), destination_grid%lonCell(j), &
                                                   start_cell, source_grid%ncells, source_grid%maxEdges, &
                                                   source_grid%nEdgesOnCell, source_grid%cellsOnCell , &
                                                   source_grid%latCell, source_grid%lonCell)
   start_cell = nearest_cell_to_output_point(j)
end do

! compute interpolation weights for all the destination/output grid cell centers
allocate(interp_weights(3,destination_grid%ncells))
allocate(interp_cells(3,destination_grid%ncells))
call calc_horiz_interp_weights(source_grid,destination_grid%ncells,destination_grid%latCell, &
                               destination_grid%lonCell, weight_file, interp_weights, interp_cells)

if ( do_vertical_interpolation ) then
   write(*,*)''
   write(*,*)'Calculating vertical interpolation weights for later use'
   write(*,*)''

   nz = destination_grid%nvertlevelsp1
   allocate(zgrid_interp(nz,destination_grid%ncells))
   allocate(vert_interp_weights_edges(2,nz,destination_grid%ncells))
   allocate(vert_interp_levels_edges(2,nz,destination_grid%ncells))

   ! First, interpolate height to the output grid, then find the weights needed to interpolate
   !  to the destination grid's height levels
   call interpolate_mpas(interp_weights, interp_cells, source_grid%ncells, nz, destination_grid%ncells, &
      source_grid%zgrid, zgrid_interp, 'zgrid', .false., .false. )
   call calc_vert_interp_weights(nz, destination_grid%ncells, zgrid_interp, &
        destination_grid%zgrid, vert_interp_weights_edges, vert_interp_levels_edges)
   deallocate(zgrid_interp)

   ! Now do the same but for unstaggered zgrid
   nz = destination_grid%nvertlevels
   allocate(zgrid_interp(nz,destination_grid%ncells))
   allocate(vert_interp_weights(2,nz,destination_grid%ncells))
   allocate(vert_interp_levels(2,nz,destination_grid%ncells))
   call interpolate_mpas(interp_weights, interp_cells, source_grid%ncells, nz, destination_grid%ncells, &
      source_grid%zgrid_unstaggered, zgrid_interp, 'zgrid_unstaggered', .false., .false. )
   call calc_vert_interp_weights(nz, destination_grid%ncells, zgrid_interp, &
        destination_grid%zgrid_unstaggered, vert_interp_weights, vert_interp_levels)
   deallocate(zgrid_interp)

endif

! Open the file with variables we want to interpolate
! Make sure the number of cells in the file we want to interpolate is identical
!  to the number of cells in the template file for the source grid
call open_netcdf(file_to_interpolate, ncfileid)
call get_netcdf_dims(ncfileid,'nCells',ncells)
if ( source_grid%ncells .ne. ncells ) then
   write(*,*)'dimension mismatch between file_to_interpolate and source_grid_template_file.'
   call close_netcdf(file_to_interpolate, ncfileid)
   call cleanup
   stop
endif

! Open the destination file.  There are some fields that we will want to copy
!  from this file to the output, controlled by namelist variable variables_to_copy_from_destination_file
call open_netcdf(destination_grid_template_file, ncfileid_dest)

! Create output file. Define all dimensions and variables, but don't yet fill with data
! This opens a file for writing; referenced by ncidout, which is intent(inout) and can be used
call define_output_file_from_template(ncfileid,output_file,destination_grid,netcdf_output_type,ncidout)

! Get the number of variables in the file we want to interpolate and loop over them
! Interpolate each variable, or, in special cases, just copy the variable to the output file
call get_netcdf_info(ncfileid,ndims,nvars) ! ndims not currently used

do_winds = .true.
do v = 1,nvars

   dims(:) = 1 ! reset each time through loop
   dim_names(:) = ''
   rcode = nf90_Inquire_Variable(ncfileid, v, varname, xtype, ndims_var, dimids, natts) ! Output is varname, xtype,ndims_var, dimids

   write(*,fmt='(a)') ''
   write(*,fmt='(a)') 'Processing '//trim(adjustl(varname))
   do j = 1,ndims_var
      rcode = nf90_inquire_dimension( ncfileid, dimids(j), name=dim_names(j), len=dims(j) )
      dim_names(j) = trim(adjustl(dim_names(j))) ! can't do trim on an array, so do it here
      write(*,fmt='(a)') trim(adjustl(dim_names(j)))
   enddo
   write(*,fmt='(a15,i2,4i8)')'ndims, dims = ',ndims_var,dims(1:ndims_var)
   if ( ndims_var .ge. 4 ) then
      write(*,*)'ndims_var = ',ndims_var,' not allowed. cycle'
      cycle
   endif

   ! Error checking. Everything should be on cell centers (%ncells) except for u
   if ( any(dim_names == 'nEdges') ) then
      if ( trim(adjustl(varname)).ne."u") then
         write(*,*)'Variable '//trim(adjustl(varname))//' is defined on edges but is not u.'
         write(*,*)'This is curious, so stopping.'
         call close_netcdf(file_to_interpolate, ncfileid)
         call close_netcdf(destination_grid_template_file, ncfileid_dest)
         call close_netcdf(output_file, ncfileid)
         call cleanup
         stop
      endif
   else if ( any(dim_names == 'nCells') ) then
      continue
   else
      if ( ndims_var > 0 ) then
         if ( xtype .ne. nf90_char ) then
            write(*,*)'Integer or float variable '//trim(adjustl(varname))//' is not defined on edges or cell centers.'
            write(*,*)'This is curious, so stopping.'
            call close_netcdf(file_to_interpolate, ncfileid)
            call close_netcdf(destination_grid_template_file, ncfileid_dest)
            call close_netcdf(output_file, ncfileid)
            call cleanup
            stop
         endif
      endif
   endif

   ! Figure out vertical dimension for this variable
   ! Same for both source_grid and destination_grid, per earlier check
   if ( any(dim_names == 'nSoilLevels') ) then
      nz = source_grid%nsoillevels
   else if ( any(dim_names == 'nVertLevels') ) then
      nz = source_grid%nvertlevels
   else if ( any(dim_names == 'nVertLevelsP1') ) then
      nz = source_grid%nvertlevelsp1
   else
      nz = 1
   endif

   ! Interpolate the current variable.  What we do depends on the variable type
   if ( xtype == nf90_float .or. xtype == nf90_double ) then
     
      ! Handle all the wind variables at once
      if ( any( trim(adjustl(varname)).eq.(/"u","uReconstructZonal","uReconstructMeridional"/)) ) then

         if ( do_winds ) then
            allocate(input_u(nz,source_grid%ncells), input_v(nz,source_grid%ncells))
            allocate(edge_normal_wind(nz,destination_grid%nedges))
            call get_netcdf_var(ncfileid,'uReconstructZonal',      nz, source_grid%ncells, input_u)
            call get_netcdf_var(ncfileid,'uReconstructMeridional', nz, source_grid%ncells, input_v)

            ! Interpolate zonal and meridional winds to the new grid, then derive edge wind from interpolated U/V
            allocate(output_u(nz,destination_grid%ncells), output_v(nz,destination_grid%ncells))
            call interpolate_mpas(interp_weights, interp_cells, source_grid%ncells, nz, destination_grid%ncells, &
                        input_u, output_u, 'uReconstructZonal', print_before_and_after_interp_values, do_vertical_interpolation)
            call interpolate_mpas(interp_weights, interp_cells, source_grid%ncells, nz, destination_grid%ncells, &
                        input_v, output_v, 'uReconstructMeridional', print_before_and_after_interp_values, do_vertical_interpolation)
            call project_to_edges(destination_grid,output_u,output_v,edge_normal_wind)

            call output_variable(ncidout,destination_grid%ncells,nz,output_u,'uReconstructZonal') ! output interpolated zonal wind
            call output_variable(ncidout,destination_grid%ncells,nz,output_v,'uReconstructMeridional') ! output interpolated meridional wind
            call output_variable(ncidout,destination_grid%nedges,nz,edge_normal_wind,'u') ! output interpolated cell edge wind

            if ( print_before_and_after_interp_values ) then
               allocate(u(nz,source_grid%nedges))
               call get_netcdf_var(ncfileid,'u', nz, source_grid%nedges, u) ! just need for printing
               call print_min_max('u',nz, source_grid%nedges, u, destination_grid%nedges, edge_normal_wind)
               deallocate(u)
            endif

            do_winds = .false. ! u,uReconstructZonal, and uReconstructMeridional are all done

            deallocate(input_u,input_v,output_u,output_v,edge_normal_wind)
         else 
            write(*,fmt='(a)') 'Already processed '//trim(adjustl(varname))
         endif ! do_winds 

      else ! not wind

         allocate(input_data_real(nz,source_grid%ncells))
         allocate(output_data_real(nz,destination_grid%ncells))

         ! Just copy any requested variables, without doing any interpolation
         if ( any( trim(adjustl(varname)).eq.variables_to_copy_from_destination_file) ) then
            write(*,fmt='(a)')'Copying '//trim(adjustl(varname))//' from '//trim(adjustl(destination_grid_template_file))//' to output'
            if ( ndims_var == 1 .or. ndims_var == 2 ) then ! ndims_var == 1 could happen for 'ter'
               call get_netcdf_var(ncfileid_dest,varname,destination_grid%ncells,output_data_real(1,:))
            endif
            if ( ndims_var == 3 ) then
               call get_netcdf_var(ncfileid_dest,varname, nz, destination_grid%ncells, output_data_real)
            endif

         else ! We are doing barycentric interpolation

            if ( ndims_var == 1 .or. ndims_var == 2 ) then ! ndims_var == 1 could happen for 'ter'
               call get_netcdf_var(ncfileid,varname,source_grid%ncells,input_data_real(1,:))
            endif
            if ( ndims_var == 3 ) then
              call get_netcdf_var(ncfileid,varname, nz, source_grid%ncells, input_data_real)
            endif

           ! Some variables, are encoded as floats, but really are integers
           !   We want to know if they are actually integers, so we can do
           !   nearest neighbor interpolation on them, rather than barycentric interpolation
            call check_actual_type(source_grid%ncells, nz, input_data_real, is_integer)

            ! Nearest neighbor interpolation for integers
            if ( is_integer ) then ! nearest neighbor interpolation for integers
               write(*,*)'Variable '//trim(adjustl(varname))//' is actually an integer'
               do j = 1,destination_grid%ncells
                  output_data_real(:,j) = input_data_real(:,nearest_cell_to_output_point(j))
               enddo
               if ( print_before_and_after_interp_values ) & 
                   call print_min_max(varname,nz, source_grid%ncells, input_data_real, destination_grid%ncells, output_data_real)
            else
               ! barcyentric interpolation for floats/doubles
               call interpolate_mpas(interp_weights, interp_cells, source_grid%ncells, nz, destination_grid%ncells, &
                        input_data_real, output_data_real, varname, print_before_and_after_interp_values, do_vertical_interpolation)
            endif
         endif

         call output_variable(ncidout,destination_grid%ncells,nz,output_data_real,varname)

         deallocate(input_data_real,output_data_real)
      endif

   else if ( xtype == nf90_int .or. xtype == nf90_int64 ) then
      allocate(input_data_int(nz,source_grid%ncells))
      allocate(output_data_int(nz,destination_grid%ncells))

      ! Just copy any requested variables, without doing any interpolation
      if ( any( trim(adjustl(varname)).eq.variables_to_copy_from_destination_file) ) then
         write(*,fmt='(a)')'Copying '//trim(adjustl(varname))//' from '//trim(adjustl(destination_grid_template_file))//' to output'
         if ( ndims_var == 1 .or. ndims_var == 2 ) then ! ndims_var == 1 could happen for 'ter'
            call get_netcdf_var(ncfileid_dest,varname,destination_grid%ncells,output_data_int(1,:))
         endif
         if ( ndims_var == 3 ) then
            call get_netcdf_var(ncfileid_dest,varname, nz, destination_grid%ncells, output_data_int)
         endif

         call output_variable(ncidout,destination_grid%ncells,nz,output_data_int,varname)

      else if ( ndims_var == 0 ) then ! scalar integer variables; just copy them over to the new file

         call get_netcdf_var(ncfileid,varname, input_data_int(1,1)) ! basically a call to get_netcdf_var_scalar_integer
         call output_variable(ncidout,input_data_int(1,1),varname) ! basically a call to output_variable_scalar

      else ! Nearest neighbor interpolation for integer variable

         if ( ndims_var == 1 .or. ndims_var == 2 ) then
            call get_netcdf_var(ncfileid,varname,source_grid%ncells,input_data_int(1,:))
         endif
         if ( ndims_var == 3 ) then
            call get_netcdf_var(ncfileid,varname, nz, source_grid%ncells, input_data_int)
         endif

         ! nearest neighbor interpolation
         do j = 1,destination_grid%ncells
            output_data_int(:,j) = input_data_int(:,nearest_cell_to_output_point(j))
         enddo
         if ( print_before_and_after_interp_values ) &
             call print_min_max(varname,nz, source_grid%ncells, real(input_data_int), destination_grid%ncells, real(output_data_int))

         call output_variable(ncidout,destination_grid%ncells,nz,output_data_int,varname)
      endif

      deallocate(input_data_int,output_data_int)

   else if ( xtype == nf90_char ) then
      ! just copy character variables over to the new file
      char_len = dims(1) ! could be 2d, (/char_len,Time/), but time doesn't matter
      call get_netcdf_var(ncfileid,varname,char_len,data_char)
      call output_variable(ncidout,char_len,trim(adjustl(data_char)),varname) ! basically a call to output_variable_scalar
   else
      write(*,*)'Unknown variable type for variable '//trim(adjustl(varname))
      call close_netcdf(file_to_interpolate, ncfileid)
      call close_netcdf(destination_grid_template_file, ncfileid_dest)
      call close_netcdf(output_file, ncfileid)
      call cleanup
      stop
   endif

enddo ! end loop over variables (loop over v)

! clean up
if ( do_vertical_interpolation ) then
   deallocate(vert_interp_weights_edges, vert_interp_levels_edges )
   deallocate(vert_interp_weights, vert_interp_levels )
endif
call close_netcdf(file_to_interpolate, ncfileid) ! close the open netCDF files
call close_netcdf(destination_grid_template_file, ncfileid_dest)
call close_netcdf(output_file, ncidout)          
call cleanup

write(*,*)'All done'

stop

contains

subroutine cleanup
   if (allocated(nearest_cell_to_output_point)) deallocate(nearest_cell_to_output_point)
   if (allocated(interp_weights)) deallocate(interp_weights)
   if (allocated(interp_cells)) deallocate(interp_cells)
   call deallocate_mpas_grid(source_grid)
   call deallocate_mpas_grid(destination_grid)
end subroutine cleanup

end program mpas2mpas

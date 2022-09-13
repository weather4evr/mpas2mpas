module namelist_mod

implicit none

! Everything is public
!private

! namelist variables
character(len=500)  :: source_grid_template_file = 'dum'
character(len=500)  :: destination_grid_template_file = 'dum'
character(len=500)  :: file_to_interpolate = 'dum' 
character(len=500)  :: output_file = 'dum'
character(len=500)  :: weight_file = 'dum'
character(len=500)  :: netcdf_output_type = 'netcdf4' ! options are 'classic', 'netcdf4', 'cdf5'
character(len=500), dimension(100) :: variables_to_copy_from_destination_file = ''
logical :: do_vertical_interpolation = .false.
logical :: print_before_and_after_interp_values = .false.

namelist /share/ source_grid_template_file, destination_grid_template_file, &
                 file_to_interpolate, output_file, weight_file, &
                 netcdf_output_type, print_before_and_after_interp_values, &
                 do_vertical_interpolation, variables_to_copy_from_destination_file

contains

subroutine read_namelist

   integer :: iunit = 10

   ! need to initialize this here until coming up with a better way of doing things
   variables_to_copy_from_destination_file(1:4) = (/ 'rho_base','theta_base','ter','xland' /)

   open(file='namelist.input',unit=iunit, status = 'old')
   read(iunit,nml=share)
   close(iunit)

end subroutine read_namelist

end module namelist_mod

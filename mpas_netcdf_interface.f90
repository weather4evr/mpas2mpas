module mpas_netcdf_interface

use kinds, only : i_kind, r_double, r_kind
use netcdf
use pnetcdf, only : nf_64bit_data ! needed to output CDF-5 type...unclear if we can get from netcdf

implicit none

private

public :: mpas_grid, init_mpas_grid, deallocate_mpas_grid
public :: open_netcdf, close_netcdf, define_output_file_from_template
public :: get_netcdf_var, get_netcdf_dims, get_netcdf_info
public :: output_variable, sanity_checks
public :: source_grid, destination_grid

interface get_netcdf_var
   module procedure get_netcdf_var_char
   module procedure get_netcdf_var_scalar_integer
   module procedure get_netcdf_var_1d_integer
   module procedure get_netcdf_var_1d_real
   module procedure get_netcdf_var_1d_double
   module procedure get_netcdf_var_2d_integer
   module procedure get_netcdf_var_2d_real
   module procedure get_netcdf_var_2d_double
end interface

interface output_variable
   module procedure output_variable_scalar
   module procedure output_variable_real
   module procedure output_variable_int
   module procedure output_variable_char
end interface

type mpas_grid 
    !character(len=19)                              :: config_start_time
     integer                                        :: ncells
     integer                                        :: nedges
     integer                                        :: nvertices
     integer                                        :: nvertlevels
     integer                                        :: nvertlevelsp1
     integer                                        :: nsoillevels
     integer                                        :: maxEdges  
     integer                                        :: vertexDegree
     integer                                        :: isice_lu
     integer                                        :: iswater_lu
     integer,        dimension(:),      allocatable :: nEdgesOnCell
     integer,        dimension(:,:),    allocatable :: cellsOnCell
     integer,        dimension(:,:),    allocatable :: verticesOnCell
     integer,        dimension(:,:),    allocatable :: cellsOnVertex
     integer,        dimension(:,:),    allocatable :: edgesOnCell
     real(r_double), dimension(:),      allocatable :: latEdge
     real(r_double), dimension(:),      allocatable :: lonEdge
     real(r_double), dimension(:),      allocatable :: latCell
     real(r_double), dimension(:),      allocatable :: lonCell
     real(r_double), dimension(:),      allocatable :: latVertex
     real(r_double), dimension(:),      allocatable :: lonVertex
     real(r_double), dimension(:,:),    allocatable :: edgeNormalVectors
     real(r_kind),   dimension(:,:),    allocatable :: zgrid
     real(r_kind),   dimension(:,:),    allocatable :: zgrid_unstaggered
     real(r_kind),   dimension(:,:),    allocatable :: zs
     real(r_kind),   dimension(:,:),    allocatable :: dzs
end type mpas_grid 

! public variables
type(mpas_grid) :: source_grid, destination_grid

! variables visible to this module only
character(len=100) :: DIMSNAME,VARSNAME,ATTSNAME
integer(i_kind) :: ncstatus
integer(i_kind) :: ncidin,ncidout,ndims,nvars,ngatts,unlimdimid,idims,dimsval,ivars,idims2,ivars2
integer(i_kind) :: varstype,varsndims,varsdimids(4),varsnatts, ivarsnatts,igatts
logical         :: fexist

contains

subroutine open_netcdf(fname,ncfileid)
   character(len=*), intent(in) :: fname
   integer(i_kind), intent(out) :: ncfileid

   inquire(file=trim(adjustl(fname)),exist=fexist)
   if ( .not. fexist ) then
      write(*,fmt='(a)') trim(adjustl(fname))//' does not exist'
      stop
   endif

   write(*,fmt='(a)') 'Opening '//trim(adjustl(fname))

   ncstatus = nf90_open(path=trim(adjustl(fname)),mode=nf90_nowrite,ncid=ncfileid)  ! open file
   if ( ncstatus .ne. 0 ) then
      write(*,fmt='(a)') 'error reading '//trim(adjustl(fname))
      stop
   endif
end subroutine open_netcdf

subroutine close_netcdf(fname,ncfileid)
   character(len=*), intent(in) :: fname
   integer(i_kind), intent(in) :: ncfileid
   ncstatus = nf90_close(ncfileid) ! close file
   if ( ncstatus .ne. 0 ) then
      write(*,fmt='(a)') 'error closing '//trim(adjustl(fname))
      stop
   endif
   write(*,fmt='(a)') 'Closing '//trim(adjustl(fname))
end subroutine close_netcdf

subroutine init_mpas_grid(fname,grid)

   character(len=*), intent(in)   :: fname
   type(mpas_grid), intent(inout) :: grid

   integer :: ncfileid, ncdimid, ncvarid
   integer :: dim1, dim2, k

   ! Read MPAS file for grid info
   call open_netcdf(fname,ncfileid)
   if (ncstatus /= 0) then
      write(*,*) ' error opening mpas file '//trim(fname)
      stop
   end if

   ! Get dimensions (integer values)
   call get_netcdf_dims(ncfileid,'nCells',grid%ncells)
   call get_netcdf_dims(ncfileid,'nEdges',grid%nedges)
   call get_netcdf_dims(ncfileid,'nVertices',grid%nvertices)
   call get_netcdf_dims(ncfileid,'nVertLevels',grid%nvertlevels)
   call get_netcdf_dims(ncfileid,'nVertLevelsP1',grid%nvertlevelsp1)
   call get_netcdf_dims(ncfileid,'nSoilLevels',grid%nsoillevels)
   call get_netcdf_dims(ncfileid,'maxEdges',grid%maxEdges)
   call get_netcdf_dims(ncfileid,'vertexDegree',grid%vertexDegree)

   ! Get some attributes
  !ncstatus = nf90_get_att(ncfileid,nf90_global,'config_start_time', grid%config_start_time)

   ! Get MPAS latitudes and longitudes
   dim1 = grid%nedges
   allocate( grid%latEdge(dim1), grid%lonEdge(dim1))
   call get_netcdf_var(ncfileid,'latEdge',dim1,grid%latEdge)
   call get_netcdf_var(ncfileid,'lonEdge',dim1,grid%lonEdge)

   dim1 = grid%ncells
   allocate( grid%latCell( dim1 ), grid%lonCell( dim1 ) )
   call get_netcdf_var(ncfileid,'latCell',dim1,grid%latCell)
   call get_netcdf_var(ncfileid,'lonCell',dim1,grid%lonCell)

   dim1 = grid%nvertices
   allocate( grid%latVertex( dim1 ), grid%lonVertex( dim1 ) )
   call get_netcdf_var(ncfileid,'latVertex',dim1,grid%latVertex)
   call get_netcdf_var(ncfileid,'lonVertex',dim1,grid%lonVertex)

   dim1 = grid%maxEdges
   dim2 = grid%ncells
   allocate( grid%verticesOnCell( dim1, dim2))
   call get_netcdf_var(ncfileid,'verticesOnCell',dim1, dim2, grid%verticesOnCell)

   dim1 = grid%vertexDegree ! should be 3
   dim2 = grid%nvertices
   allocate( grid%cellsOnVertex( dim1, dim2) )
   call get_netcdf_var(ncfileid,'cellsOnVertex',dim1, dim2, grid%cellsOnVertex)

   dim1 = grid%ncells
   allocate( grid%nEdgesOnCell( dim1 ))
   call get_netcdf_var(ncfileid,'nEdgesOnCell',dim1,grid%nEdgesOnCell)

   dim1 = grid%maxEdges
   dim2 = grid%ncells
   allocate( grid%cellsOnCell( dim1, dim2))
   call get_netcdf_var(ncfileid,'cellsOnCell',dim1, dim2, grid%cellsOnCell)

   dim1 = grid%maxEdges
   dim2 = grid%ncells
   allocate( grid%edgesOnCell( dim1, dim2))
   call get_netcdf_var(ncfileid,'edgesOnCell',dim1, dim2, grid%edgesOnCell)

   dim1 = 3
   dim2 = grid%nedges
   allocate( grid%edgeNormalVectors( dim1, dim2))
   call get_netcdf_var(ncfileid,'edgeNormalVectors',dim1, dim2, grid%edgeNormalVectors)

   dim1 = grid%nvertlevelsp1
   dim2 = grid%ncells
   allocate(grid%zgrid( dim1, dim2))
   call get_netcdf_var(ncfileid,'zgrid',dim1, dim2, grid%zgrid)

   dim1 = grid%nvertlevels
   dim2 = grid%ncells
   allocate(grid%zgrid_unstaggered( dim1, dim2))
   do k = 1,grid%nvertlevels
      grid%zgrid_unstaggered(k,:) = 0.5 * (grid%zgrid(k,:) + grid%zgrid(k+1,:))
   enddo

   dim1 = grid%nsoillevels
   dim2 = grid%ncells
   allocate( grid%zs( dim1, dim2 ), grid%dzs( dim1, dim2) )
   call get_netcdf_var(ncfileid,'zs',dim1, dim2, grid%zs)
   call get_netcdf_var(ncfileid,'dzs',dim1, dim2, grid%dzs)

   ! these are scalar integers
   call get_netcdf_var(ncfileid,'isice_lu', grid%isice_lu) 
   call get_netcdf_var(ncfileid,'iswater_lu', grid%iswater_lu)

   call close_netcdf(fname,ncfileid)

end subroutine init_mpas_grid

subroutine get_netcdf_dims(fileid,variable,output)
   integer, intent(in) :: fileid
   character(len=*), intent(in) :: variable
   integer, intent(out) :: output

   integer :: ncdimid

   ncstatus = nf90_inq_dimid(fileid,trim(adjustl(variable)),ncdimid)
   if ( ncstatus /= 0 ) write(*,*) 'Error inq_dimid '//trim(adjustl(variable))
   ncstatus = nf90_inquire_dimension(fileid,ncdimid, len=output)
   if ( ncstatus /= 0 ) write(*,*) 'Error inquire_dimension '//trim(adjustl(variable))
end subroutine get_netcdf_dims

subroutine get_netcdf_info(ncfileid,ndims,nvars)
   integer, intent(in) :: ncfileid
   integer, intent(inout) :: ndims,nvars
   ncstatus=nf90_inquire(ncfileid,ndims,nvars,ngatts,unlimdimid)
   if ( ncstatus /= 0 ) then
      write(*,*) 'Error in get_netcdf_info'
      write(*,*) 'ncstatus = ',ncstatus
      stop
   endif
end subroutine get_netcdf_info

subroutine get_netcdf_var_char(fileid,variable,dim1,output)

   integer, intent(in) :: fileid, dim1
   character(len=*), intent(in) :: variable
   character(len=dim1), intent(inout) :: output

   integer :: ncvarid, istatus

   istatus = 0
   ncstatus = nf90_inq_varid(fileid,trim(adjustl(variable)),ncvarid) 
   if ( ncstatus /= 0 ) write(*,*) 'Error inq_varid '//trim(adjustl(variable)) ; istatus = istatus + ncstatus
   ncstatus = nf90_get_var(fileid,ncvarid,output) 
   if ( ncstatus /= 0 ) write(*,*) 'Error inq_get_var '//trim(adjustl(variable)) ; istatus = istatus + ncstatus

   if ( istatus /= 0 ) then
      write(*,*)'Error in get_netcdf_var_char'
      stop
   endif
end subroutine get_netcdf_var_char

subroutine get_netcdf_var_1d_real(fileid,variable,dim1,output)

   integer, intent(in) :: fileid, dim1
   character(len=*), intent(in) :: variable
   real(r_kind), intent(inout), dimension(dim1)  :: output

   integer :: ncvarid, istatus

   istatus = 0
   ncstatus = nf90_inq_varid(fileid,trim(adjustl(variable)),ncvarid) 
   if ( ncstatus /= 0 ) write(*,*) 'Error inq_varid '//trim(adjustl(variable)) ; istatus = istatus + ncstatus
   ncstatus = nf90_get_var(fileid,ncvarid,output) 
   if ( ncstatus /= 0 ) write(*,*) 'Error inq_get_var '//trim(adjustl(variable)) ; istatus = istatus + ncstatus

   if ( istatus /= 0 ) then
      write(*,*)'Error in get_netcdf_var_1d_real'
      stop
   endif
end subroutine get_netcdf_var_1d_real

subroutine get_netcdf_var_1d_double(fileid,variable,dim1,output)

   integer, intent(in) :: fileid, dim1
   character(len=*), intent(in) :: variable
   real(r_double), intent(inout), dimension(dim1)  :: output

   integer :: ncvarid, istatus

   istatus = 0
   ncstatus = nf90_inq_varid(fileid,trim(adjustl(variable)),ncvarid) 
   if ( ncstatus /= 0 ) write(*,*) 'Error inq_varid '//trim(adjustl(variable)) ; istatus = istatus + ncstatus
   ncstatus = nf90_get_var(fileid,ncvarid,output) 
   if ( ncstatus /= 0 ) write(*,*) 'Error inq_get_var '//trim(adjustl(variable)) ; istatus = istatus + ncstatus

   if ( istatus /= 0 ) then
      write(*,*)'Error in get_netcdf_var_1d_double'
      stop
   endif
end subroutine get_netcdf_var_1d_double

subroutine get_netcdf_var_1d_integer(fileid,variable,dim1,output)

   integer, intent(in) :: fileid, dim1
   character(len=*), intent(in) :: variable
   integer, intent(inout), dimension(dim1)  :: output

   integer :: ncvarid, istatus

   istatus = 0
   ncstatus = nf90_inq_varid(fileid,trim(adjustl(variable)),ncvarid) 
   if ( ncstatus /= 0 ) write(*,*) 'Error inq_varid '//trim(adjustl(variable)) ; istatus = istatus + ncstatus
   ncstatus = nf90_get_var(fileid,ncvarid,output) 
   if ( ncstatus /= 0 ) write(*,*) 'Error inq_get_var '//trim(adjustl(variable)) ; istatus = istatus + ncstatus

   if ( istatus /= 0 ) then
      write(*,*)'Error in get_netcdf_var_1d_integer'
      stop
   endif
end subroutine get_netcdf_var_1d_integer

subroutine get_netcdf_var_2d_real(fileid,variable,dim1,dim2,output)

   integer, intent(in) :: fileid, dim1,dim2
   character(len=*), intent(in) :: variable
   real(r_kind), intent(inout), dimension(dim1,dim2)  :: output

   integer :: ncvarid, istatus

   istatus = 0
   ncstatus = nf90_inq_varid(fileid,trim(adjustl(variable)),ncvarid) 
   if ( ncstatus /= 0 ) write(*,*) 'Error inq_varid '//trim(adjustl(variable)) ; istatus = istatus + ncstatus
   ncstatus = nf90_get_var(fileid,ncvarid,output) 
   if ( ncstatus /= 0 ) write(*,*) 'Error inq_get_var '//trim(adjustl(variable)) ; istatus = istatus + ncstatus

   if ( istatus /= 0 ) then
      write(*,*)'Error in get_netcdf_var_2d_real'
      stop
   endif
end subroutine get_netcdf_var_2d_real

subroutine get_netcdf_var_2d_double(fileid,variable,dim1,dim2,output)

   integer, intent(in) :: fileid, dim1,dim2
   character(len=*), intent(in) :: variable
   real(r_double), intent(inout), dimension(dim1,dim2)  :: output

   integer :: ncvarid, istatus

   istatus = 0
   ncstatus = nf90_inq_varid(fileid,trim(adjustl(variable)),ncvarid) 
   if ( ncstatus /= 0 ) write(*,*) 'Error inq_varid '//trim(adjustl(variable)) ; istatus = istatus + ncstatus
   ncstatus = nf90_get_var(fileid,ncvarid,output) 
   if ( ncstatus /= 0 ) write(*,*) 'Error inq_get_var '//trim(adjustl(variable)) ; istatus = istatus + ncstatus

   if ( istatus /= 0 ) then
      write(*,*)'Error in get_netcdf_var_2d_double'
      stop
   endif
end subroutine get_netcdf_var_2d_double

subroutine get_netcdf_var_2d_integer(fileid,variable,dim1,dim2,output)

   integer, intent(in) :: fileid, dim1,dim2
   character(len=*), intent(in) :: variable
   integer, intent(inout), dimension(dim1,dim2)  :: output

   integer :: ncvarid, istatus

   istatus = 0
   ncstatus = nf90_inq_varid(fileid,trim(adjustl(variable)),ncvarid) 
   if ( ncstatus /= 0 ) write(*,*) 'Error inq_varid '//trim(adjustl(variable)); istatus = istatus + ncstatus
   ncstatus = nf90_get_var(fileid,ncvarid,output) 
   if ( ncstatus /= 0 ) write(*,*) 'Error inq_get_var '//trim(adjustl(variable)) ; istatus = istatus + ncstatus

   if ( istatus /= 0 ) then
      write(*,*)'Error in get_netcdf_var_2d_integer'
      stop
   endif
end subroutine get_netcdf_var_2d_integer

subroutine get_netcdf_var_scalar_integer(fileid,variable,output)

   integer, intent(in) :: fileid
   character(len=*), intent(in) :: variable
   integer, intent(inout) :: output

   integer :: ncvarid, istatus

   istatus = 0
   ncstatus = nf90_inq_varid(fileid,trim(adjustl(variable)),ncvarid) 
   if ( ncstatus /= 0 ) write(*,*) 'Error inq_varid '//trim(adjustl(variable)); istatus = istatus + ncstatus

   ncstatus = nf90_get_var(fileid,ncvarid,output) 
   if ( ncstatus /= 0 ) write(*,*) 'Error inq_get_var '//trim(adjustl(variable)) ; istatus = istatus + ncstatus

   if ( istatus /= 0 ) then
      write(*,*)'Error in get_netcdf_var_scalar_integer'
      stop
   endif
end subroutine get_netcdf_var_scalar_integer

subroutine output_variable_scalar(fileid,output_data,variable)
   integer, intent(in) :: fileid
   integer, intent(in) :: output_data
   character(len=*), intent(in) :: variable

   integer :: ncvarid

   ncstatus = nf90_inq_varid(fileid,trim(adjustl(variable)),ncvarid)
   if ( ncstatus /= 0 ) write(*,*) 'Error inq_varid '//trim(adjustl(variable)) 

   ncstatus = nf90_put_var(fileid,ncvarid,output_data)
   if ( ncstatus /= 0 ) write(*,*) 'Error put_var '//trim(adjustl(variable)) 
end subroutine output_variable_scalar

subroutine output_variable_real(fileid,npoints,nlevs,output_data,variable)
   integer, intent(in) :: fileid, npoints, nlevs
   real(r_kind), intent(in), dimension(nlevs,npoints) :: output_data
   character(len=*), intent(in) :: variable

   integer :: ncvarid

   ncstatus = nf90_inq_varid(fileid,trim(adjustl(variable)),ncvarid)
   if ( ncstatus /= 0 ) write(*,*) 'Error inq_varid '//trim(adjustl(variable)) 

   if ( nlevs.gt.1 ) then  ! 3-d variable
      ncstatus = nf90_put_var(fileid,ncvarid,output_data)
   else
      ncstatus = nf90_put_var(fileid,ncvarid,output_data(1,:))
   endif
   if ( ncstatus /= 0 ) write(*,*) 'Error put_var '//trim(adjustl(variable)) 
end subroutine output_variable_real

subroutine output_variable_int(fileid,npoints,nlevs,output_data,variable)
   integer, intent(in) :: fileid, npoints, nlevs
   integer, intent(in), dimension(nlevs,npoints) :: output_data
   character(len=*), intent(in) :: variable

   integer :: ncvarid

   ncstatus = nf90_inq_varid(fileid,trim(adjustl(variable)),ncvarid)
   if ( ncstatus /= 0 ) write(*,*) 'Error inq_varid '//trim(adjustl(variable)) 

   if ( nlevs.gt.1 ) then  ! 3-d variable
      ncstatus = nf90_put_var(fileid,ncvarid,output_data)
   else
      ncstatus = nf90_put_var(fileid,ncvarid,output_data(1,:))
   endif
   if ( ncstatus /= 0 ) write(*,*) 'Error put_var '//trim(adjustl(variable)) 
end subroutine output_variable_int

subroutine output_variable_char(fileid,char_len,output_data,variable)
   integer, intent(in) :: fileid, char_len
   character(len=char_len), intent(in) :: output_data
   character(len=*), intent(in) :: variable

   integer :: ncvarid

   ncstatus = nf90_inq_varid(fileid,trim(adjustl(variable)),ncvarid)
   if ( ncstatus /= 0 ) write(*,*) 'Error inq_varid '//trim(adjustl(variable)) 
   ncstatus = nf90_put_var(fileid,ncvarid,output_data)
   if ( ncstatus /= 0 ) write(*,*) 'Error put_var '//trim(adjustl(variable)) 
end subroutine output_variable_char

subroutine define_output_file_from_template(ncidin,fout,output_grid,netcdf_output_type,ncidout)
   integer(i_kind), intent(in)    :: ncidin
   character(len=*), intent(in) :: fout, netcdf_output_type
   type(mpas_grid), intent(in)    :: output_grid
   integer(i_kind), intent(inout) :: ncidout

   integer(i_kind) :: netcdf_file_type, group_ncid, i
   integer(i_kind) :: varids(1000)
   character(len=500) :: group_name

   if ( trim(adjustl(netcdf_output_type)) == 'netcdf4' ) then
      netcdf_file_type = NF90_NETCDF4  ! NETCDF4 file
   else if ( trim(adjustl(netcdf_output_type)) == 'cdf5' ) then
      netcdf_file_type = nf_64bit_data ! CDF-5 format
   else if ( trim(adjustl(netcdf_output_type)) == 'classic' ) then
      netcdf_file_type = NF90_CLOBBER
   else
      write(*,*)trim(adjustl(netcdf_output_type))//' is not allowed.'
      write(*,*)'Can be either netcdf4, cdf5, or classic.'
      write(*,*)'Try again.'
      stop
   endif

!  ncstatus=nf90_open(path=trim(adjustl(fin)),mode=nf90_nowrite,ncid=ncidin)
!  ncstatus=nf90_create(path=trim(adjustl(fout)),cmode=nf90_clobber,ncid=ncidout)
   ncstatus=nf90_create(path=trim(adjustl(fout)),cmode=netcdf_file_type,ncid=ncidout)

   if ( ncstatus == 0 ) then
      write(*,*)'Successfully opened '//trim(adjustl(fout))//' for writing.'
      write(*,*)'Output format is '//trim(adjustl(netcdf_output_type))
   endif

   ! dimensions
   ncstatus=nf90_inquire(ncidin,ndims,nvars,ngatts,unlimdimid)
   do idims=1,ndims
      ncstatus=nf90_inquire_dimension(ncidin,idims,DIMSNAME,dimsval)
      if ( trim(adjustl(DIMSNAME)) == 'nCells') then
         ncstatus=nf90_def_dim(ncidout,trim(adjustl(DIMSNAME)),output_grid%ncells,idims2)
      else if ( trim(adjustl(DIMSNAME)) == 'nEdges') then
         ncstatus=nf90_def_dim(ncidout,trim(adjustl(DIMSNAME)),output_grid%nedges,idims2)
      else if ( trim(adjustl(DIMSNAME)) == 'nVertices') then
         ncstatus=nf90_def_dim(ncidout,trim(adjustl(DIMSNAME)),output_grid%nvertices,idims2)
      else if ( trim(adjustl(DIMSNAME)) == 'Time') then
         ncstatus=nf90_def_dim(ncidout,trim(adjustl(DIMSNAME)),NF90_UNLIMITED,idims2) ! Time gets unlimited dimension
      else
         ncstatus=nf90_def_dim(ncidout,trim(adjustl(DIMSNAME)),dimsval,idims2)
      endif
   end do

   ! variables
   do ivars=1,nvars
      ncstatus=nf90_inquire_variable(ncidin,ivars,VARSNAME,varstype,varsndims,varsdimids,varsnatts)
      ncstatus=nf90_def_var(ncidout,trim(adjustl(VARSNAME)),varstype,varsdimids(1:varsndims),ivars2)
     !ncstatus=nf90_def_var(ncidout,VARSNAME,varstype,varsndims,varsdimids,ivars)
      do ivarsnatts=1,varsnatts
         ncstatus=nf90_inq_attname(ncidin,ivars,ivarsnatts,ATTSNAME)
         ncstatus=nf90_copy_att(ncidin,ivars,ATTSNAME,ncidout,ivars2)
      end do
   end do

   ! global attributes
   do igatts=1,ngatts
      ncstatus=nf90_inq_attname(ncidin,nf90_global,igatts,ATTSNAME)
      ncstatus=nf90_copy_att(ncidin,nf90_global,ATTSNAME,ncidout,nf90_global)
   end do
   ncstatus=nf90_enddef(ncidout)
end subroutine define_output_file_from_template

subroutine sanity_checks(grid1, grid2, mismatch)
   type(mpas_grid), intent(in) :: grid1, grid2
   logical, intent(inout) :: mismatch

   integer :: k
   mismatch = .false. ! assume everything is all good

   ! Sanity check to make sure number of vertical levels and soil levels match
   if ( (grid1%nvertlevels   .ne. grid2%nvertlevels) .or. &
        (grid1%nvertlevelsp1 .ne. grid2%nvertlevelsp1) .or. &
        (grid1%nsoillevels   .ne. grid2%nsoillevels) ) then
      write(*,*)'dimension mismatch'
      write(*,*)'nvertlevels src/dst = ',  grid1%nvertlevels,   grid2%nvertlevels
      write(*,*)'nvertlevelsp1 src/dst = ',grid1%nvertlevelsp1, grid2%nvertlevelsp1
      write(*,*)'nsoillevels src/dst = ',  grid1%nsoillevels,   grid2%nsoillevels
      mismatch = .true.
      return
   endif

   if ( (grid1%isice_lu    .ne. grid2%isice_lu) .or. &
        (grid1%iswater_lu  .ne. grid2%iswater_lu) ) then
      write(*,*)'ice/water lookup mismatch'
      write(*,*)'isice_lu src/dst = ',   grid1%isice_lu,   grid2%isice_lu
      write(*,*)'iswater_lu src/dst = ', grid1%iswater_lu, grid2%iswater_lu
      mismatch = .true.
      return
   endif
   
   if ( any(grid1%zs(:,1) - grid2%zs(:,1) .ne. 0) ) then
      write(*,*)'grid%zs mismatch'
      write(*,*)'grid1%zs(:,1)/grid2%zs(:,1) = ',grid1%zs(:,1), grid2%zs(:,1)
      mismatch = .true.
      return
   endif

   if ( any(grid1%dzs(:,1) - grid2%dzs(:,1) .ne. 0) ) then
      write(*,*)'grid%dzs mismatch'
      write(*,*)'grid1%dzs(:,1)/grid2%dzs(:,1) = ',grid1%dzs(:,1), grid2%dzs(:,1)
      mismatch = .true.
      return
   endif

end subroutine sanity_checks

subroutine deallocate_mpas_grid(grid)

 type(mpas_grid), intent(inout) :: grid

 if(allocated(grid%latCell)) deallocate(grid%latCell)
 if(allocated(grid%lonCell)) deallocate(grid%lonCell)
 if(allocated(grid%latVertex)) deallocate(grid%latVertex)
 if(allocated(grid%lonVertex)) deallocate(grid%lonVertex)
 if(allocated(grid%latEdge)) deallocate(grid%latEdge)
 if(allocated(grid%lonEdge)) deallocate(grid%lonEdge)
 if(allocated(grid%nEdgesOnCell)) deallocate(grid%nEdgesOnCell)
 if(allocated(grid%cellsOnCell)) deallocate(grid%cellsOnCell)
 if(allocated(grid%verticesOnCell)) deallocate(grid%verticesOnCell)
 if(allocated(grid%cellsOnVertex)) deallocate(grid%cellsOnVertex)
 if(allocated(grid%edgesOnCell)) deallocate(grid%edgesOnCell)
 if(allocated(grid%edgeNormalVectors)) deallocate(grid%edgeNormalVectors)
 if(allocated(grid%zgrid)) deallocate(grid%zgrid)
 if(allocated(grid%zgrid_unstaggered)) deallocate(grid%zgrid_unstaggered)
 if(allocated(grid%zs)) deallocate(grid%zs)
 if(allocated(grid%dzs)) deallocate(grid%dzs)

end subroutine deallocate_mpas_grid

end module mpas_netcdf_interface

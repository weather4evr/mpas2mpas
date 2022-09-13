module interpolation_stuff

use kinds, only : r_single, r_double, r_kind
use mpas_netcdf_interface, only : mpas_grid, source_grid, destination_grid

implicit none

private

! public subroutines
public :: interpolate_mpas, nearest_cell
public :: project_to_edges, check_actual_type, print_min_max
public :: calc_horiz_interp_weights, calc_vert_interp_weights

! public variables
public :: vert_interp_weights_edges, vert_interp_weights, vert_interp_levels_edges, vert_interp_levels

real(r_kind), allocatable, dimension(:,:,:)  :: vert_interp_weights_edges, vert_interp_weights
integer,      allocatable, dimension(:,:,:)  :: vert_interp_levels_edges, vert_interp_levels

!-----------------------------------------------
contains

subroutine interpolate_mpas(interp_weights,interp_cells,nCells,nVertLevels,ncells_target,input_data,output_data,varname, &
                                 do_print, do_vertical_interpolation)

   implicit none

   ! variables passed to this routine
   real (r_double), dimension(3, ncells_target), intent(in) :: interp_weights
   integer, dimension(3, ncells_target), intent(in)         :: interp_cells
   integer, intent(in)                                   :: nCells, nVertLevels, ncells_target
   character(len=*), intent(in)                            :: varname
   real(r_kind), dimension(nVertLevels, nCells), intent(inout) :: input_data ! data to interpolate; inout because vertical interpolation could change it
   logical, intent(in)                         :: do_print, do_vertical_interpolation
   real(r_kind), dimension(nVertLevels,ncells_target), intent(inout) :: output_data  ! data that has been interpolated

   ! local variables
   integer :: j, k
   logical :: vertical_interpolation_flag
   real(r_kind), dimension(nVertLevels, ncells_target) :: height_array ! destination grid heights
   real(r_kind), dimension(nVertLevels, ncells_target) :: tmp_val
   real(r_kind), dimension(2,nVertLevels) :: data_array
   real(r_kind), dimension(nVertLevels, nCells) :: source_height_array !source grid heights
   real(r_kind), dimension(nVertLevels, ncells_target) :: source_height ! source grid heights, interpolated onto destination grid
   real(r_kind) :: target_z

   ! local flag for vertical interpolation
   vertical_interpolation_flag = do_vertical_interpolation ! from namelist

   if ( nVertLevels == destination_grid%nvertlevelsp1 ) then
      height_array = destination_grid%zgrid
      source_height_array = source_grid%zgrid
   else if ( nVertLevels == destination_grid%nvertlevels ) then
      height_array = destination_grid%zgrid_unstaggered
      source_height_array = source_grid%zgrid_unstaggered
   else if ( nVertLevels == destination_grid%nsoillevels ) then
      vertical_interpolation_flag = .false. ! no vertical interpolation for soil levels
   else if ( nVertLevels == 1 ) then
      vertical_interpolation_flag = .false. ! no vertical interpolation for 1-d fields
   endif

   call interp_delaunay( interp_weights, interp_cells,  &
                         input_data, output_data,       &
                         nVertLevels, nCells, ncells_target )

        ! Now interpolate vertically
        ! This is our baseline that we think is correct
        ! But it's slow to do this for each variable, so try to precalculate weights instead
   if ( vertical_interpolation_flag ) then
      write(*,*)'Doing vertical interpolation for '//trim(adjustl(varname))
      do j=1,ncells_target
         source_height(:,j) = interp_weights(1,j)*source_height_array(:,interp_cells(1,j))  &
                            + interp_weights(2,j)*source_height_array(:,interp_cells(2,j))  &
                            + interp_weights(3,j)*source_height_array(:,interp_cells(3,j))
         data_array(1,:) = source_height(:,j) ! interpolated value of the source grid height; effective height after interpolation
         data_array(2,:) = output_data(:,j) ! data on destination grid
         do k=1,nVertLevels
            target_z = height_array(k,j)
            output_data(k,j) = vertical_interp(target_z, nVertLevels, data_array, extrap=0)
         enddo
      enddo
   endif

!  This way uses the weights, and it's much faster
!  if ( vertical_interpolation_flag ) then
!     write(*,*)'Doing vertical interpolation for '//trim(adjustl(varname))
!     do j=1,ncells_target
!        do k=1,nVertLevels
!           if ( nVertLevels == destination_grid%nvertlevelsp1 ) then
!              !Use the weights we found earlier
!              tmp_val(k,j) = vert_interp_weights_edges(1,k,j)*output_data(vert_interp_levels_edges(1,k,j),j) + &
!                             vert_interp_weights_edges(2,k,j)*output_data(vert_interp_levels_edges(2,k,j),j)
!           else if ( nVertLevels == destination_grid%nvertlevels ) then
!              tmp_val(k,j) = vert_interp_weights(1,k,j)*output_data(vert_interp_levels(1,k,j),j) + &
!                             vert_interp_weights(2,k,j)*output_data(vert_interp_levels(2,k,j),j)
!           endif
!                           !!! value_interp    = wm*zf(2,lm) + wp*zf(2,lp)
!        end do
!     end do ! loop over k
!     output_data = tmp_val
!  endif

   if ( do_print ) then
      call print_min_max(varname,nVertLevels,nCells,input_data,ncells_target,output_data)
   endif

end subroutine interpolate_mpas

subroutine interp_delaunay( interp_weights, interp_cells,    &
                            value_cell, value_interp,        &
                            nVertLevels, nCells, ncells_target)

   implicit none

   integer, intent(in) :: nCells, nVertLevels, ncells_target
   real (r_kind), dimension(nVertLevels, nCells), intent(inout) :: value_cell ! inout because interpolation could change it
   real (r_kind), dimension(nVertLevels, ncells_target), intent(out) :: value_interp
   real (r_double), dimension(3, ncells_target), intent(in) :: interp_weights
   integer, dimension(3, ncells_target), intent(in) :: interp_cells

   integer :: j, k

   do j=1,ncells_target
   do k=1,nVertLevels
      value_interp(k,j) = interp_weights(1,j)*value_cell(k,interp_cells(1,j))  &
                             + interp_weights(2,j)*value_cell(k,interp_cells(2,j))  &
                             + interp_weights(3,j)*value_cell(k,interp_cells(3,j))
   end do
   end do

end subroutine interp_delaunay

subroutine calc_horiz_interp_weights(grid,npoints_target,target_lats,target_lons,weight_file,interp_weights,interp_cells)
   type(mpas_grid), intent(in) :: grid
   integer, intent(in) :: npoints_target
   real(r_double), dimension(npoints_target), intent(in) :: target_lats, target_lons
   character(len=*), intent(in) :: weight_file
   real(r_double), dimension(3,npoints_target), intent(inout)  :: interp_weights
   integer, dimension(3,npoints_target), intent(inout)         :: interp_cells

   logical :: fexist
   integer :: iunit_weight = 42
   integer :: start_vertex, rcode, j
   integer, allocatable, dimension(:) :: nearest_vertex_to_output_points

   inquire(file=trim(adjustl(weight_file)), exist=fexist)
   if ( fexist ) then
      open(unit=iunit_weight, file=trim(adjustl(weight_file)), status='old', form='unformatted', iostat=rcode)
      if ( rcode .ne. 0 ) write(*,*) 'Error opening weights file '//trim(adjustl(weight_file))
      read(iunit_weight) interp_weights, interp_cells
      close(iunit_weight)
      write(*,fmt='(a)') 'Reading interpolation weights from '//trim(adjustl(weight_file))
   else
      ! Find the vertex on our source grid that is nearest to each of our desired output points
      write(*,*) 'computing interpolation weights '
      allocate(nearest_vertex_to_output_points(npoints_target)) ! number of target points
      start_vertex = 1
      do j=1,npoints_target
         nearest_vertex_to_output_points(j) = nearest_vertex( target_lats(j), target_lons(j), start_vertex, grid%ncells, &
                                              grid%nvertices, grid%maxEdges, &
                                              grid%nEdgesOnCell, grid%verticesOnCell, grid%cellsOnVertex, &
                                              grid%latCell, grid%lonCell, grid%latVertex, grid%lonVertex )
         start_vertex = nearest_vertex_to_output_points(j)
      end do

      call compute_interp_weights_all ( target_lats, target_lons, nearest_vertex_to_output_points, grid%ncells, &
                                     grid%nvertices, grid%cellsOnVertex, grid%latCell, grid%lonCell, &
                                     interp_weights, interp_cells, npoints_target )
      write(*,*) ' interpolation weights complete '
      open(unit=iunit_weight, file=trim(adjustl(weight_file)), status='new', form='unformatted', iostat=rcode)
      write(iunit_weight) interp_weights, interp_cells
      close(iunit_weight)
      deallocate(nearest_vertex_to_output_points)
   endif
end subroutine calc_horiz_interp_weights

   integer function nearest_vertex( target_lat, target_lon, &
                                     start_vertex, &
                                     nCells, nVertices, maxEdges, &
                                     nEdgesOnCell, verticesOnCell, &
                                     cellsOnVertex, latCell, lonCell, &
                                     latVertex, lonVertex )

        implicit none

        real(r_double), intent(in) :: target_lat, target_lon
        integer, intent(in) :: start_vertex
        integer, intent(in) :: nCells, nVertices, maxEdges
        integer, dimension(nCells), intent(in) :: nEdgesOnCell
        integer, dimension(maxEdges,nCells), intent(in) :: verticesOnCell
        integer, dimension(3,nVertices), intent(in) :: cellsOnVertex
        real(r_double), dimension(nCells), intent(in) :: latCell, lonCell
        real(r_double), dimension(nVertices), intent(in) :: latVertex, lonVertex


        integer :: i, cell1, cell2, cell3, iCell
        integer :: iVtx
        integer :: current_vertex
        real(r_double) :: cell1_dist, cell2_dist, cell3_dist
        real(r_double) :: current_distance, d
        real(r_double) :: nearest_distance

        nearest_vertex = start_vertex
        current_vertex = -1

        do while (nearest_vertex /= current_vertex)
            current_vertex = nearest_vertex
            current_distance = sphere_distance(latVertex(current_vertex), lonVertex(current_vertex), &
                                               target_lat, target_lon,                1.0_r_double)
            nearest_vertex = current_vertex
            nearest_distance = current_distance
            cell1 = cellsOnVertex(1,current_vertex)
            cell2 = cellsOnVertex(2,current_vertex)
            cell3 = cellsOnVertex(3,current_vertex)
            cell1_dist = sphere_distance(latCell(cell1), lonCell(cell1), target_lat, target_lon, 1.0_r_double)
            cell2_dist = sphere_distance(latCell(cell2), lonCell(cell2), target_lat, target_lon, 1.0_r_double)
            cell3_dist = sphere_distance(latCell(cell3), lonCell(cell3), target_lat, target_lon, 1.0_r_double)
            if (cell1_dist < cell2_dist) then
                if (cell1_dist < cell3_dist) then
                    iCell = cell1
                else
                    iCell = cell3
                end if
            else
                if (cell2_dist < cell3_dist) then
                    iCell = cell2
                else
                    iCell = cell3
                end if
            end if
            do i = 1, nEdgesOnCell(iCell)
                iVtx = verticesOnCell(i,iCell)
                d = sphere_distance(latVertex(iVtx), lonVertex(iVtx), target_lat, target_lon, 1.0_r_double)
                if (d < nearest_distance) then
                    nearest_vertex = iVtx
                    nearest_distance = d
                end if
            end do
        end do

    end function nearest_vertex

    real (r_double) function sphere_distance(lat1, lon1, lat2, lon2, radius)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Compute the great-circle distance between (lat1, lon1) and (lat2, lon2) on
    ! a sphere with given radius.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       implicit none

        real (r_double), intent(in) :: lat1, lon1, lat2, lon2, radius

        real (r_double) :: arg1

        arg1 = sqrt( sin(0.5_r_double*(lat2-lat1))**2 +  &
                     cos(lat1)*cos(lat2)*sin(0.5_r_double*(lon2-lon1))**2 )
        sphere_distance = 2.*radius*asin(arg1)

    end function sphere_distance


    real (r_double) function triangle_area(p1_lat,p1_lon,p2_lat,p2_lon,p3_lat,p3_lon,radius)

        implicit none

        real (r_double) p1_lat,p1_lon,p2_lat,p2_lon,p3_lat,p3_lon,radius
        real (r_double) a,b,c,s,e,tanqe

        a = sphere_distance(p1_lat,p1_lon,p2_lat,p2_lon,radius)
        b = sphere_distance(p2_lat,p2_lon,p3_lat,p3_lon,radius)
        c = sphere_distance(p3_lat,p3_lon,p1_lat,p1_lon,radius)
        s = 0.5*(a+b+c)

        tanqe = sqrt(tan(0.5_r_double*s)*tan(0.5_r_double*(s-a))*tan(0.5_r_double*(s-b))*tan(0.5_r_double*(s-c)))
        e = 4.*atan(tanqe)
        triangle_area = radius*radius*e

    end function triangle_area


    subroutine compute_interp_weights_all( target_lat, target_lon, nearest_vtx, &
                                        nCells, nVertices, cellsOnVertex, &
                                        lat_cell, lon_cell, &
                                        interp_weights, interp_cells, ncells_target)

        implicit none

        integer, intent(in) :: ncells_target
        integer, dimension(ncells_target), intent(in) :: nearest_vtx
        real (r_double), dimension(ncells_target), intent(in) :: target_lat, target_lon
        integer, intent(in) :: nCells, nVertices
        integer, dimension(3,nVertices), intent(in) :: cellsOnVertex
        real (r_double), dimension(nCells), intent(in) :: lat_cell, lon_cell
        real (r_double), intent(out), dimension(3,ncells_target) :: interp_weights
        integer, intent(out), dimension(3,ncells_target) :: interp_cells

        real (r_double), dimension(3) :: weights

        integer :: i, j, k
        real (r_double), parameter :: eps = 1.e-010_r_double
        real (r_double), parameter :: invpower = 1.0_r_double

        real (r_double) :: lat0, lon0, lat1, lon1, lat2, lon2, lat3, lon3
        real (r_double) :: area123, area012, area023, area013, sum_weights
        real (r_double) :: pii

        pii = acos(-1.0)     !  Value of PI From Henry W's code, Defined at top of module
       !pii = 2.0 * asin(1.0)   ! From Mike Duda's Gives same value as above

        do j=1,ncells_target

        lat0 = target_lat(j)
        lon0 = target_lon(j)
        lat1 = lat_cell(cellsOnVertex(1,nearest_vtx(j)))
        lon1 = lon_cell(cellsOnVertex(1,nearest_vtx(j)))
        lat2 = lat_cell(cellsOnVertex(2,nearest_vtx(j)))
        lon2 = lon_cell(cellsOnVertex(2,nearest_vtx(j)))
        lat3 = lat_cell(cellsOnVertex(3,nearest_vtx(j)))
        lon3 = lon_cell(cellsOnVertex(3,nearest_vtx(j)))

        if (lon0 < 0) lon0 = lon0+2.*pii
        if (lon1 < 0) lon1 = lon1+2.*pii
        if (lon2 < 0) lon2 = lon2+2.*pii
        if (lon3 < 0) lon3 = lon3+2.*pii

        area123 = triangle_area(lat1,lon1,lat2,lon2,lat3,lon3,1.0_r_double)
        area012 = triangle_area(lat0,lon0,lat1,lon1,lat2,lon2,1.0_r_double)
        area023 = triangle_area(lat0,lon0,lat2,lon2,lat3,lon3,1.0_r_double)
        area013 = triangle_area(lat0,lon0,lat1,lon1,lat3,lon3,1.0_r_double)

        !
        !  check areas
        !
        weights(1) = area023/area123
        weights(2) = area013/area123
        weights(3) = area012/area123

        sum_weights = weights(1)+weights(2)+weights(3)

        weights(1) = weights(1)/sum_weights
        weights(2) = weights(2)/sum_weights
        weights(3) = weights(3)/sum_weights

        do k=1,3
            interp_cells(k,j) = cellsOnVertex(k,nearest_vtx(j))
            interp_weights(k,j) = weights(k)
        end do

        end do

    end subroutine compute_interp_weights_all

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Finds the MPAS grid cell nearest to (target_lat, target_lon)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    integer function nearest_cell(target_lat, target_lon, start_cell, nCells, maxEdges, &
                              nEdgesOnCell, cellsOnCell, latCell, lonCell)

    implicit none

    real (r_double), intent(in) :: target_lat, target_lon
    integer, intent(in) :: start_cell
    integer, intent(in) :: nCells, maxEdges
    integer, dimension(nCells), intent(in) :: nEdgesOnCell
    integer, dimension(maxEdges,nCells), intent(in) :: cellsOnCell
    real (r_double), dimension(nCells), intent(in) :: latCell, lonCell

    integer :: i
    integer :: iCell
    integer :: current_cell
    real (r_double) :: current_distance, d
    real (r_double) :: nearest_distance

    nearest_cell = start_cell
    current_cell = -1

    do while (nearest_cell /= current_cell)
        current_cell = nearest_cell
        current_distance = sphere_distance(latCell(current_cell), lonCell(current_cell), target_lat, &
                                           target_lon, 1.0_r_double)
        nearest_cell = current_cell
        nearest_distance = current_distance
        do i = 1, nEdgesOnCell(current_cell)
            iCell = cellsOnCell(i,current_cell)
            if (iCell <= nCells) then
                d = sphere_distance(latCell(iCell), lonCell(iCell), target_lat, target_lon, 1.0_r_double)
                if (d < nearest_distance) then
                    nearest_cell = iCell
                    nearest_distance = d
                end if
            end if
        end do
    end do

    end function nearest_cell

 ! Modified from mpas_init_atm_cases.F
 ! Not currently used, but helpful to see how the weights are constructed
 real (r_kind) function vertical_interp(target_z, nz, zf, extrap)

      real(r_kind), intent(in) :: target_z
      integer, intent(in) :: nz
      real(r_kind), dimension(2,nz), intent(in) :: zf  ! zf(1,:) is column of vertical coordinate values, zf(2,:) is column of field values
      integer, intent(in), optional :: extrap   ! can take values 0 = constant, 1 = linear (default), 2 = lapse-rate
      integer :: k, lm, lp
      real(r_kind) :: wm, wp
      real(r_kind) :: slope

      integer :: extrap_type

      if (present(extrap)) then
         extrap_type = extrap
      else
         extrap_type = 1
      end if

      ! Extrapolation required
      if (target_z < zf(1,1)) then
         if (extrap_type == 0) then
            vertical_interp = zf(2,1)
         else if (extrap_type == 1) then
            slope = (zf(2,2) - zf(2,1)) / (zf(1,2) - zf(1,1))
            vertical_interp = zf(2,1) + slope * (target_z - zf(1,1))
         else if (extrap_type == 2) then
            vertical_interp = zf(2,1) - (target_z - zf(1,1))*0.0065
         end if
         return
      end if
      if (target_z >= zf(1,nz)) then
         if (extrap_type == 0) then
            vertical_interp = zf(2,nz)
         else if (extrap_type == 1) then
            slope = (zf(2,nz) - zf(2,nz-1)) / (zf(1,nz) - zf(1,nz-1))
            vertical_interp = zf(2,nz) + slope * (target_z - zf(1,nz))
         else if (extrap_type == 2) then
             write(*,*)'extrap_type == 2 not implemented for target_z >= zf(1,nz)'
         end if
         return
      end if

      ! No extrapolation required
      do k=1,nz-1
         if (target_z >= zf(1,k) .and. target_z < zf(1,k+1)) then
            lm = k
            lp = k+1
            wm = (zf(1,k+1) - target_z) / (zf(1,k+1) - zf(1,k))
            wp = (target_z - zf(1,k))   / (zf(1,k+1) - zf(1,k))
            exit
         end if
      end do

      vertical_interp = wm*zf(2,lm) + wp*zf(2,lp)

      return

   end function vertical_interp

   ! This subroutine calculates vertical interpolation weights and saves them for future use
   !  This is just to save time, rather than calling vertical_interp for each variable.
   !  By saving the weights, they can be calculated once, and then used for all variables with
   !  the same vertical grid.  
   ! Note that the weights need to be computed separately for staggered and unstaggered vertical
   !  levels.
   ! For logic about how the weights are ultimately applied, see 
   !     vertical_interp, which is  modified from mpas_init_atm_cases.F
   subroutine calc_vert_interp_weights(nz,npoints,zgrid_interpolated,zgrid_destination, vert_interp_weights, vert_interp_levels)

      integer, intent(in) :: nz, npoints ! npoints is on the destination grid
      real(r_kind), dimension(nz,npoints), intent(in) :: zgrid_interpolated, zgrid_destination
      real(r_kind), dimension(2, nz, npoints), intent(inout) :: vert_interp_weights !(2,nz, destination_grid%ncells)
      integer,      dimension(2, nz, npoints), intent(inout) :: vert_interp_levels  !(2,nz, destination_grid%ncells)
      integer :: i, j, k, kk, lm, lp
      real(r_kind) :: wm, wp, target_z
      real(r_kind), dimension(nz) :: zf  ! zf(:) is column of vertical coordinate values

      do j = 1,npoints
         zf(:) = zgrid_interpolated(:,j) ! interpolated value of the source grid height at each level for the jth point
         do k = 1,nz
            target_z = zgrid_destination(k,j)

            if (target_z < zf(1)) then
              !vertical_interp = zf(2,1) ! value at lowest model level
               vert_interp_weights(1,k,j) = 1.0 ! wm
               vert_interp_weights(2,k,j) = 0.0 ! wp
               vert_interp_levels(1,k,j) = 1  ! lm 
               vert_interp_levels(2,k,j) = -999 ! lp ... in this case, a dummy
               cycle
            end if
            if (target_z >= zf(nz)) then
              !vertical_interp = zf(2,nz) ! value at highest model level
               vert_interp_weights(1,k,j) = 1.0 ! wm
               vert_interp_weights(2,k,j) = 0.0 ! wp
               vert_interp_levels(1,k,j) = nz  ! lm 
               vert_interp_levels(2,k,j) = -999 ! lp ... in this case, a dummy
               cycle
            end if

            do kk=1,nz-1
               if (target_z >= zf(kk) .and. target_z < zf(kk+1)) then
                  lm = kk
                  lp = kk+1
                  wm = (zf(kk+1) - target_z) / (zf(kk+1) - zf(kk))
                  wp = (target_z - zf(kk))   / (zf(kk+1) - zf(kk))
                  vert_interp_weights(1,k,j) = wm
                  vert_interp_weights(2,k,j) = wp
                  vert_interp_levels(1,k,j) = lm
                  vert_interp_levels(2,k,j) = lp
                  exit
               end if
            end do

         end do ! loop over k
      end do ! loop over j

     !Ultimate value
     !vertical_interp = wm*zf(2,lm) + wp*zf(2,lp)

   end subroutine calc_vert_interp_weights

subroutine print_min_max(varname,nz,nhoriz1,data1,nhoriz2,data2)
   character(len=*), intent(in) :: varname
   integer, intent(in) :: nz, nhoriz1,nhoriz2
   real(r_kind), dimension(nz,nhoriz1), intent(in) :: data1
   real(r_kind), dimension(nz,nhoriz2), intent(in) :: data2

   integer :: k

   write(*,fmt='(a)') ' Information for '//trim(adjustl(varname))
   do k = 1,nz
      write(*,fmt='(a,i2,2(2x,f16.5,2x))')'min before/after interp k= ',k,minval(data1(k,:)),minval(data2(k,:))
      write(*,fmt='(a,i2,2(2x,f16.5,2x))')'max before/after interp k= ',k,maxval(data1(k,:)),maxval(data2(k,:))
   enddo
   write(*,*) '  '
end subroutine print_min_max

subroutine check_actual_type(nCells, nVertLevels, input_data, is_integer)
   integer, intent(in)  :: nCells, nVertLevels
   real(r_kind), dimension(nVertLevels,nCells), intent(in) :: input_data
   logical, intent(out) :: is_integer

   integer, dimension(nVertLevels,nCells) :: idata
   integer :: i,k

   ! turn input into an integer
   do k = 1,nVertLevels
      do i = 1,nCells
         idata(k,i) = int(input_data(k,i))
      enddo
   enddo

   ! see if the difference between the float and integer is 0 everywhere
   if (all( (idata-input_data) == 0)) then
      is_integer = .true.
   else
      is_integer = .false.
   endif
end subroutine check_actual_type

subroutine project_to_edges(grid,zonal,meridional,edge_normal_wind)

   type(mpas_grid), intent(in) :: grid
   real(r_kind),  dimension(grid%nvertlevels,grid%ncells), intent(in) :: zonal, meridional
   real(r_kind),  dimension(grid%nvertlevels,grid%nedges), intent(inout) :: edge_normal_wind

   integer :: jEdge, iEdge, i, j, k
   real(r_double),  dimension(:,:), allocatable :: east, north

   ! This code taken from DART models/mpas_atm/model_mod, in turn taken from MPAS

   allocate( east(3,grid%ncells))
   allocate( north(3,grid%ncells))

   do i = 1, grid%ncells
      east(1,i) = -sin(grid%lonCell(i))
      east(2,i) =  cos(grid%lonCell(i))
      east(3,i) =  dble(0.0)
      call r3_normalize(east(1,i), east(2,i), east(3,i))

      north(1,i) = -cos(grid%lonCell(i))*sin(grid%latCell(i))
      north(2,i) = -sin(grid%lonCell(i))*sin(grid%latCell(i))
      north(3,i) =  cos(grid%latCell(i))
      call r3_normalize(north(1,i), north(2,i), north(3,i))
   enddo

 ! Project data from the cell centers to the edges
   edge_normal_wind = 0.0 ! dble(0.0)
   do i = 1, grid%ncells
      do jEdge = 1, grid%nEdgesOnCell(i)
         iEdge = grid%edgesOnCell(jEdge, i)
         do k = 1, grid%nvertlevels
           edge_normal_wind(k,iEdge) = edge_normal_wind(k,iEdge) + 0.5 * zonal(k,i)   &
                     * (grid%edgeNormalVectors(1,iEdge) * east(1,i)  &
                     +  grid%edgeNormalVectors(2,iEdge) * east(2,i)  &
                     +  grid%edgeNormalVectors(3,iEdge) * east(3,i)) &
                     + 0.5 * meridional(k,i)            &
                     * (grid%edgeNormalVectors(1,iEdge) * north(1,i) &
                     +  grid%edgeNormalVectors(2,iEdge) * north(2,i) &
                     +  grid%edgeNormalVectors(3,iEdge) * north(3,i))
         enddo
      enddo
   enddo
   deallocate(east,north)

end subroutine project_to_edges

subroutine r3_normalize(ax, ay, az) !CSS added from DART models/mpas_atm/model_mod, in turn taken from MPAS

   !normalizes the vector (ax, ay, az)

   real(r_double), intent(inout) :: ax, ay, az
   real(r_double) :: mi

    mi = 1.0 / sqrt(ax**2 + ay**2 + az**2)
    ax = ax * mi
    ay = ay * mi
    az = az * mi

end subroutine r3_normalize

end module interpolation_stuff

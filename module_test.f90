!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!
!>      VUA and IPSL/LSCE by the iLOVECLIM / LUDUS coding group / Within the LUDUS code environement
!
!       LICENSING TERMS:
!>      \copyright
!!      This file is part of [insert sub-component name here, in following Foobar]
!!      Foobar is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
!!      as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!!
!!      Foobar is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!!      of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!!
!!      You should have received a copy of the GNU General Public License along with Foobar.
!!      If not, see <http://www.gnu.org/licenses/>.
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

! dmr DGC is the interpolation and distance over earth simple module
! --- 1 is using the test code pertaining to this module
! --- 0 is not using it
! ---   NOTA: does not make sense to use it with TOMS_760
! dmr

#define DGC_USE 0

! dmr TOMS_760 is the bicubic interpolation module from TOMS journal
! --- 1 is using the test code pertaining to this module
! --- 0 is not using it
! ---   NOTA: does not make sense to use it with DGC_USE
! dmr

#define TOMS_760 1

! dmr NCIO_USE proposes netCDF output of the test module
! --- 1 netCDF output
! --- 0 no netCDF output

#define NCIO_USE 1

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      MODULE: module_test
!
!>     @author  Didier M. Roche (dmr)
!
!>     @brief This module module test is handling the main testing of the grid_lens library
!
!>     @date Creation date: July, 30th, 2015
!>     @date Last modification: SLastChangedDate$
!>     @author Last modified by : dmr
!
!>     @version This is svn version: $LastChangedRevision$
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module module_test
      
       implicit none
       private

       public :: main

      ! NOTE_avoid_public_variables_if_possible

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: main
!
!>     @brief main call of the test module, used as a global wrapper
!
!      DESCRIPTION:
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      function main() result(returnValue)

       use grid_class, only: surface_grid, surface_grid_init
       use sub_grid_class, only: sub_grid, sub_grid_init
      
#if (DGC_USE == 1)
       use interpolate_mod, only: interpolate_init, interpolate
#endif

#if (NCIO_USE == 1)
       use ncio, only: nc_create, nc_write_attr, nc_write_dim, nc_write
#endif

#if ( TOMS_760 == 1 )
       use toms760_wrapper, only: bicubic_interpol
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       logical             :: returnValue
       
       type(surface_grid) :: lres_land_grid, hres_land_grid
       type(sub_grid)     :: zoom_grid
       logical            :: succeed = .false.
       character(len=256) :: input_file, input_file_l
       character(len=256), parameter     :: directory = "inputdata/"
       character(len=256) :: varname

#if (DGC_USE == 1)
       real(kind=8), allocatable, dimension(:,:) :: YLAT, XLONG
       real(kind=8), allocatable, dimension(:,:,:,:) :: tab_dat
       real(kind=8), allocatable, dimension(:,:,:) :: interpolatable, interpolated
       integer, parameter :: nw = 3, nz = 4, ex = 2
       integer :: i,j
       logical :: results
#endif
#if (NCIO_USE == 1)
       character(len=256) :: filename
#endif

#if ( DGC_USE == 1 || TOMS_760 > 0 )
       real(kind=8), allocatable, dimension(:,:,:) :: interpolated
       integer, parameter :: nbmois = 1
#endif

#if ( TOMS_760 == 1 )
       integer :: resulting_val, val_indx
       real(kind=8), allocatable, dimension(:,:) :: interpol_one
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!>    @bug Description of the stupid sticky bug that we know exist there but is not corrected yet!

#define LRES 1
#define HRES 1

#if ( LRES == 1 )
      input_file_l=""//trim(directory)//"regridedETOPO5-T21_out.nc"
#endif

#if ( HRES == 1 )
      input_file=""//trim(directory)//"regridedField_out-60min.nc"
#elif ( HRES == 2 )
      input_file=""//trim(directory)//"regridedField_out-30min.nc"
#elif ( HRES == 3 )
      input_file=""//trim(directory)//"regridedField_out-10min.nc"
#endif

#if ( LRES == 1 ) /*  ECBilt grid has different names for the axes */

      varname = "topo"
      succeed = lres_land_grid%set_nvars(1)
      succeed = lres_land_grid%surf_grid_vars(lres_land_grid%indx_elevation)%var_init(.true.,varname)

      varname = "topo"
      succeed = hres_land_grid%set_nvars(1)
      succeed = hres_land_grid%surf_grid_vars(lres_land_grid%indx_elevation)%var_init(.true.,varname)

#endif

      hres_land_grid%is_subgrid = .true.

! Setting up the low resolution grid from input_file_l
      succeed = surface_grid_init(lres_land_grid,input_file_l)
      write(*,*) "Setting Low resolution grid: ", succeed
! Setting up the high resolution grid from input_file_l
      succeed = surface_grid_init(hres_land_grid,input_file)
      write(*,*) "Setting High resolution grid: ", succeed

!      succeed = &
!      sub_grid_init(hres_land_grid,lres_land_grid,subg=zoom_grid,lat_min=33.23d0,lat_max=77.5d0,lon_min=346.5d0,lon_max=30.9d0)

      succeed = &
      !sub_grid_init(hres_land_grid,lres_land_grid,subg=zoom_grid,lat_min=35.0d0,lat_max=80.0d0,lon_min=340.0d0,lon_max=40.0d0)
      sub_grid_init(hres_land_grid,lres_land_grid,subg=zoom_grid,lat_min=-89.0d0,lat_max=89.0d0,lon_min=0.0d0,lon_max=360.0d0)

#if (DGC_USE == 1)

      allocate(YLAT(zoom_grid%n_lon,zoom_grid%n_lat))
      allocate(XLONG(zoom_grid%n_lon,zoom_grid%n_lat))
      allocate(tab_dat(zoom_grid%n_lon,zoom_grid%n_lat,nw,nz))

      do i=1,zoom_grid%n_lon
         YLAT(i,:) = zoom_grid%p_lats(:)
      enddo

      do j=1,zoom_grid%n_lat
         XLONG(:,j) = zoom_grid%p_lons(:)
      enddo

      results = interpolate_init(zoom_grid%n_lon,zoom_grid%n_lat,lres_land_grid%n_lat,lres_land_grid%n_lon,lres_land_grid%p_lats &
               ,lres_land_grid%p_lons,YLAT,XLONG,nw,nz,ex,tab_dat)

      allocate(interpolatable(lres_land_grid%n_lon,lres_land_grid%n_lat,nbmois))

      interpolatable(:,:,1) = lres_land_grid%surf_grid_vars(lres_land_grid%indx_elevation)%var_data(:,:)

      results = interpolate(tab_dat,interpolatable,interpolated,zoom_grid%n_lon,zoom_grid%n_lat,nw,nz,lres_land_grid%n_lon &
                           ,lres_land_grid%n_lat,nbmois)
#endif

#if ( DGC_USE == 1 || TOMS_760 > 0 )
      allocate(interpolated(zoom_grid%n_lon,zoom_grid%n_lat,nbmois))
#endif      


      write(*,*) "Globality of grids: ", lres_land_grid%is_global, hres_land_grid%is_global
      write(*,*) "Setup of sub_grid = ", succeed
#if (DGC_USE == 1)
      write(*,*) "Interpolation results", results
#endif

#if ( TOMS_760 == 1 )
      val_indx = lres_land_grid%indx_elevation
      resulting_val = bicubic_interpol(lres_land_grid, zoom_grid, interpol_one, val_indx)
      interpolated(:,:,1) = interpol_one(:,:)
      deallocate(interpol_one)
#endif

!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2----------|
! cdmr --- try writing out the results in a netCDF file ...
!-----|--1---------2---------3---------4---------5---------6---------7---------8---------9---------0---------1---------2----------|

#if (NCIO_USE == 1)

      filename = "outputdata/out_highres.nc"
      
      ! Create the netcdf file, write global attributes
      call nc_create(filename)
      call nc_write_attr(filename,"title",                              &
        "highres output from dgc interpolate test")
      call nc_write_attr(filename,"institution",                        &
        "CNRS/Laboratoire des Sciences du Climat et de l'Environnement | Vrije Universiteit Amsterdam")

      call nc_write_attr(filename,"creation date",             "2015-07-29")

      ! Write the dimensions (x, y), defined inline
      call nc_write_dim(filename,"x",x=zoom_grid%p_lons(:),units="degrees_east",axis="X")
      call nc_write_dim(filename,"y",x=zoom_grid%p_lats(:),units="degrees_north",axis="Y")

      call nc_write(filename,"lres_topo",interpolated(:,:,1),dim1="x",dim2="y")
      call nc_write(filename,"hres_topo",zoom_grid%surf_grid_vars(zoom_grid%indx_elevation)%var_data(:,:),dim1="x",dim2="y")

#endif

       returnValue = succeed

      end function main

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the module here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

end module module_test

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

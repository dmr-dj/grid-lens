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

!dmr -- Adding the choice of components through the pre-processing options
!#include 'choixcomposantes.h'
!dmr -- Adding the choice of components through the pre-processing options

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      MODULE: toms760-wrapper
!
!>     @author  Didier M. Roche (dmr)
!
!>     @brief Module toms760-wrapper is preparing matrices for the toms760 bicubic interpolation for surf_grid type variables
!
!>     @date Creation date: October, 2nd, 2015
!>     @date Last modification: SLastChangedDate$
!>     @author Last modified by : dmr
!
!>     @version This is svn version: $LastChangedRevision$
!
!>     Module long description:
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module toms760_wrapper

       use grid_class,         only: surface_grid ! The surface grid variable type that will be interpolated
       use sub_grid_class,     only: sub_grid     ! The subgrid on which we are interpolating
       use Grid_Interpolation, only: rgbi3p       ! The bicubic interpolation function from TOMS760 (see toms760.f90)

       implicit none
       private

       public :: bicubic_interpol

      ! NOTE_avoid_public_variables_if_possible

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      ROUTINE: bicubic_interpol
!
!>     @brief This subroutine is performing a bicubic interpolation for a given surface_variable following TOMS760
!
!      DESCRIPTION:
!
!>     Here add the long_description of the module ...
!!      more blablas ...
!!      more blablas ...
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      function bicubic_interpol(lres_land_grid, zoom_grid, interpolated, indx_var) result(ierror_code)

!       use AnotherModule_mod, only: some_variable
!       use AnotherModule_mod, only: some_otherfunction

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  By reference variables ...
!
!>    @param[in] lres_land_grid The low resolution grid input that will be interpolated
!>    @param[in] zoom_grid The high resolution grid onto which interpolation is done. Though inout, not modified in the call i think.
!>    @param[out] interpolated The values interpolated on high res grid with bicubic toms760 method
!>    @param[out] ierror_code The return code of the function. Non zero is anomalous.
!
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       type(surface_grid)                         , intent(in)       :: lres_land_grid
       integer                                    , intent(in)       :: indx_var
       type(sub_grid)                             , intent(inout)    :: zoom_grid ! need the inout attribute in call rgbi3p. Weird.
       
       real(kind=8), allocatable, dimension(:,:)  , intent(out)      :: interpolated
       
       integer                                                       :: ierror_code

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       integer                                   :: md, iyi, ixi, ier, expa = 2
       real(kind=8)                              :: dlon
       real(kind=8), allocatable, dimension(:)   :: geo_lonex
       real(kind=8), allocatable, dimension(:,:) :: geo_expanded

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!       Main code of the function starts here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!>    @bug Description of the stupid sticky bug that we know exist there but is not corrected yet!

!dmr need to provide to TOMS_760 an expanded tab in lat, lon to help finding additional data at borders
       

       allocate(geo_lonex(lres_land_grid%n_lon+expa*2))

       geo_lonex(3:lres_land_grid%n_lon+2) = lres_land_grid%p_lons(:)
       dlon = lres_land_grid%p_lons(lres_land_grid%n_lon) - lres_land_grid%p_lons(lres_land_grid%n_lon-1)

       allocate(geo_expanded(lres_land_grid%n_lon+expa*2, lres_land_grid%n_lat))

       ! geo_expanded(3:lres_land_grid%n_lon+2,:) = lres_land_grid%s_elevation(:,:)
       geo_expanded(3:lres_land_grid%n_lon+2,:) = lres_land_grid%surf_grid_vars(indx_var)%var_data(:,:)

       geo_lonex(2) = lres_land_grid%p_lons(1)-dlon
       geo_lonex(1) = lres_land_grid%p_lons(1)-2*dlon
       geo_lonex(lres_land_grid%n_lon+3) = lres_land_grid%p_lons(lres_land_grid%n_lon)+1*dlon
       geo_lonex(lres_land_grid%n_lon+4) = lres_land_grid%p_lons(lres_land_grid%n_lon)+2*dlon

       ! geo_expanded(2,:) = lres_land_grid%s_elevation(lres_land_grid%n_lon,:)
       geo_expanded(2,:) = lres_land_grid%surf_grid_vars(indx_var)%var_data(lres_land_grid%n_lon,:)
       ! geo_expanded(1,:) = lres_land_grid%s_elevation(lres_land_grid%n_lon-1,:)
       geo_expanded(1,:) = lres_land_grid%surf_grid_vars(indx_var)%var_data(lres_land_grid%n_lon-1,:)
       !geo_expanded(lres_land_grid%n_lon+3,:) = lres_land_grid%s_elevation(1,:)
       geo_expanded(lres_land_grid%n_lon+3,:) = lres_land_grid%surf_grid_vars(indx_var)%var_data(1,:)
       !geo_expanded(lres_land_grid%n_lon+4,:) = lres_land_grid%s_elevation(2,:)
       geo_expanded(lres_land_grid%n_lon+4,:) = lres_land_grid%surf_grid_vars(indx_var)%var_data(2,:)
       write(*,*) geo_lonex
       ! read(*,*)
       
       allocate(interpolated(zoom_grid%n_lon,zoom_grid%n_lat))
       
       DO  iyi = 1, zoom_grid%n_lat
         DO  ixi = 1, zoom_grid%n_lon
           IF (ixi == 1.AND.iyi == 1) THEN
             md = 1
           ELSE
             md = 2
           END IF

           CALL rgbi3p(md,lres_land_grid%n_lon+2*expa,lres_land_grid%n_lat, geo_lonex ,lres_land_grid%p_lats &
                      ,geo_expanded, 1, zoom_grid%p_lons(ixi), zoom_grid%p_lats(iyi), interpolated(ixi,iyi), ier)
           ! IF (ier > 0) STOP
!           dzi(ixi,iyi) = zi(ixi,iyi) - zie(ixi,iyi)

           ierror_code = max(ier,ierror_code)
         END DO
       END DO       

      end function bicubic_interpol

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      End of the module here
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

end module toms760_wrapper

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      The End of All Things (op. cit.)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

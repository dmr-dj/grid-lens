       module interpolate_mod

!-----|--1--------2---------3---------4---------5---------6---------7-|
!          This program shows and example usage of the dgc module
!
!     This is version 1.1, created June, 3rd, 2013
!     Updated version 1.2, created August, 1st, 2015
!     Last modification: August, 1st, 2015
!
! Copyright 2013,2015 Didier M. Roche a.k.a. dmr
! Didier M. Roche: Didier.roche@lsce.ipsl.fr
! Release under GPLv3 license.
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!-----|--1--------2---------3---------4---------5---------6---------7-|

      IMPLICIT NONE

      CONTAINS

      logical function interpolate_init(nx,ny,nlat,nlon,latEcb,lonEcb   &
                                       ,YLAT,XLONG,nw,nz,ex,tab_dat)

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       distance_great_circle module: main routines for computation
!-----|--1--------2---------3---------4---------5---------6---------7-|
       USE dgc, ONLY: compute_bounds, find_closest_EC_cell              &
                    , which_corner, create_interp_data, expand_tab      &
                    , circular_coord

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Tweaking FLAGS to change the program's behaviour
!-----|--1--------2---------3---------4---------5---------6---------7-|

!       DEBUG >  0 triggers additional printing for debug
#define DEBUG 0
#define USE_SUBG_OBJ 0

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       GENERAL DEFINITION VARIABLES
!       nx, ny    :: highres grid
!       nlat,nlon :: lowres grid
!       nw        :: ranks for the distance
!                    1: is lat index i
!                    2: is lon index j
!                    3: is distance
!       nz        :: Number of neighbouring cells to interpolate with
!                    It has currently been tested with 4, 9, 25 and 49
!       ex        :: array expansion number
!-----|--1--------2---------3---------4---------5---------6---------7-|

       INTEGER, INTENT(IN) :: nx, ny, nlat, nlon, nw, nz, ex


!-----|--1--------2---------3---------4---------5---------6---------7-|
!       XLONG and YLAT are coordinates of the highres grid
!-----|--1--------2---------3---------4---------5---------6---------7-|

       REAL(KIND=8), DIMENSION(nx,ny) :: XLONG, YLAT

!-----|--1--------2---------3---------4---------5---------6---------7-|
!      Output of this routine == tab of weights for the computation
!-----|--1--------2---------3---------4---------5---------6---------7-|
       REAL(KIND=8), DIMENSION(nx,ny,nw,nz),INTENT(out) :: tab_dat

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       latEcb, lonEcb are coordinates of the atmospheric grid
!-----|--1--------2---------3---------4---------5---------6---------7-|
       REAL(KIND=8), DIMENSION(nlat), INTENT(in) :: latEcb
       REAL(KIND=8), DIMENSION(nlon), INTENT(in) :: lonEcb

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Variables subscripted with "ex" are extended grids
!        latexEcb, lonexEcb are the coordinates of the atm. model
!        latex_bEcb, lonex_bEcb are the boundaries of the atm. grid
!-----|--1--------2---------3---------4---------5---------6---------7-|
       REAL(KIND=8), DIMENSION(-ex+1:nlat+ex) :: latexEcb, latex_bEcb
       REAL(KIND=8), DIMENSION(-ex+1:nlon+ex) :: lonexEcb, lonex_bEcb

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Useful integers for loops mainly + file name placeholder
!-----|--1--------2---------3---------4---------5---------6---------7-|
       INTEGER :: k, i_close, j_close, corner, ii, jj, ll

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Preparing arrays for computations
!-----|--1--------2---------3---------4---------5---------6---------7-|

!       2.1: expansion of arrays (longitude IS circular, latitude NOT)

       CALL expand_tab(latEcb,nlat,ex,latexEcb,0)
       CALL expand_tab(lonEcb,nlon,ex,lonexEcb,1)

!dmr   We need to fix manually the northest and southest bound ; if not
!       the last point will be the last grid point center. 

       latexEcb(nlat+1:nlat+ex) = 90.0d0
       latexEcb(-ex+1:0) = -90.0d0

#if ( DEBUG >= 2 )
       write(*,*) latEcb
       write(*,*) latexEcb
       read(*,*)
#endif

#if ( USE_SUBG_OBJ == 0 )
!       2.2: computing cells bounds on the basis of the cells center values

       CALL compute_bounds(latexEcb,nlat+2*ex,latex_bEcb)
       CALL compute_bounds(lonexEcb,nlon+2*ex,lonex_bEcb)

!       2.3: quick fixing of the extrema: one cannot expect the last value to be
!            correct since there is initially no +90°N,S bound

       latex_bEcb(-ex+1:-ex+2) = -90.0d0
       latex_bEcb(nlat+ex-2:nlat+ex) = 90.0d0
       lonex_bEcb(nlon+ex) = lonex_bEcb(circular_coord(nlon+ex,1,nlon)) &
                           +360.0d0

#if ( DEBUG >= 2 )
       DO i=nlat+ex,-ex+1,-1
          WRITE(*,*) latexEcb(i), latex_bEcb(i)
       ENDDO

       DO i=nlon+ex,-ex+1,-1
          WRITE(*,*) lonexEcb(i), lonex_bEcb(i)
       ENDDO

       READ(*,*)
#endif


#else /* when USE_SUBG_OBJ, the bounds area read in the outer file */

!       2.2 Initialize the bound table from the one existing in sub_grid
!           NOTA: in subgrid the bounds table have dimension (nx,2)
!                 whereas here they are only nx+ex

!       2.3 Compute the expanded table(s) that are still absent

#endif
!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Main loop over the cells of the ISM
!-----|--1--------2---------3---------4---------5---------6---------7-|

       DO ii = 1, nx
         DO jj = 1, ny

       k = 0

#if ( USE_SUBG_OBJ == 0 )
!       3.1: lookup the closest atmospheric cell of ISM(ii,jj)

       CALL find_closest_EC_cell(XLONG(ii,jj), YLAT(ii,jj),latexEcb     &
                    ,latex_bEcb,nlat,nlon, lonexEcb,lonex_bEcb, i_close &
                    ,j_close,ex)

!       3.2: in case we are out the range [1,nlon] of coordinates, reset

        j_close = circular_coord(j_close,1,nlon)

#if ( DEBUG >= 0 )
       if (ii.EQ.(nx-1)) then
!       WRITE(*,*) "Obtained coords:", i_close, j_close
!       WRITE(*,*) latexEcb(i_close-2:i_close+2)
!       WRITE(*,*) "Obtained : LAT", YLAT(ii,jj), latex_bEcb(i_close-1)  &
!                    , latexEcb(i_close), latex_bEcb(i_close)
!
!       WRITE(*,*) "Obtained : LONG", XLONG(ii,jj),lonex_bEcb(j_close-1) &
!                    , lonexEcb(j_close), lonex_bEcb(j_close)
!       READ(*,*)
       endif
#endif

#else

!       3.1 the closest EC_cell has been already found in sub_grid aggreg
!           use the clustered points table to this end ...

#endif /* on USE_SUBG_OBJ */

!       3.3: need to find in which corner of the atmospheric cell we are
!            to make sure that we interpolate correctly with neighbouting
!            cells

       call which_corner(XLONG(ii,jj), YLAT(ii,jj),latEcb               &
                    ,nlat,nlon, lonEcb,i_close,j_close, corner)

!       3.4: create interpolation data with distance ...
       call create_interp_data(XLONG(ii,jj), YLAT(ii,jj), latEcb, lonEcb&
                    ,nlat,nlon,i_close,j_close, corner,ii,jj,nx,ny,nw,nz&
                    ,tab_dat)

!       3.5: spanning the grid to check whether each original point
!             is taken into account in the cornered table version

        DO ll = 1, nz
          IF ((NINT(tab_dat(ii,jj,1,ll)).EQ.i_close).AND.               &
          (NINT(tab_dat(ii,jj,2,ll)).EQ.j_close)) THEN
            k = 1
          ENDIF
        ENDDO

        interpolate_init = .true.
         IF ( k.EQ.0) THEN
             WRITE(*,*) "PROBLEM !!!! ", ii,jj
             READ(*,*)
             interpolate_init = .false.
         ENDIF
        ENDDO
      ENDDO

      return
      end function interpolate_init

      logical function interpolate(tab_dat,sxsnowG,pfGi,nx,ny,nw,nz,nlon&
                                  ,nlat,nbmois)
!       INT_MODEL defines the type of interpolation model you wish
!       INT_MODEL 0 is no distance weight (equal contribution of all cells)
!       INT_MODEL 1 depends in cell distance quadratrically
!       INT_MODEL 2 depends in cell distance with exponential model
#define INT_MODEL 2

      INTEGER, INTENT(in) :: nx,ny,nw,nz,nlon,nlat,nbmois

      REAL(KIND=8), DIMENSION(nx,ny,nw,nz),INTENT(in) :: tab_dat
!-----|--1--------2---------3---------4---------5---------6---------7-|
!       sxsnowG is the reading in climatological variable
!-----|--1--------2---------3---------4---------5---------6---------7-|
      REAL(KIND=8), DIMENSION(nlon,nlat,nbmois), INTENT(in) :: sxsnowG

      REAL(KIND=8), DIMENSION(nx,ny, nbmois), INTENT(out) :: pfGi

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       sumw is the sum of weights from the neighbouring points
!       valmax is the maximum distance over the neighbouring cells
!-----|--1--------2---------3---------4---------5---------6---------7-|
       REAL(KIND=8) :: sumw

#if ( INT_MODEL == 1 || INT_MODEL == 2 )
       REAL(KIND=8) :: valmax
#endif

       REAL(KIND=8), DIMENSION(nx,ny,nz) :: weights

       INTEGER :: i,j,ii,jj,k,ll

!       4.2: Actual interpolation

       DO k = 1, nbmois
       DO ii = 1, nx
         DO jj = 1, ny

           pfGi(ii,jj,k) = 0.0d0
           sumw = 0.0d0
           weights(ii,jj,:)=0.0

#if ( INT_MODEL == 1 )
!dmr    valmax is the maximum distance from the given cells around
           valmax  = MAXVAL(tab_dat(ii,jj,nw,:))*1.1
#elif ( INT_MODEL == 2 )
!dmr    in model 2, valmax should be the e-fold  distance (a valmin in fact)
           ! valmax  = MINVAL(tab_dat(ii,jj,nw,:))
           valmax  = 400000.0 ! in meters
#endif

!dmr    nz is the number of neighbouring cells you want to interpolate with
           DO ll = 1, nz

             i = NINT(tab_dat(ii,jj,1,ll))
             j = NINT(tab_dat(ii,jj,2,ll))

             weights(ii,jj,ll) =                                        &
#if ( INT_MODEL == 1 )
           (1-tab_dat(ii,jj,nw,ll)**2/valmax**2)
#elif ( INT_MODEL == 2 )
           EXP(1-tab_dat(ii,jj,nw,ll)/valmax)
#elif ( INT_MODEL == 0 )
           1.0d0
#endif

             pfGi(ii,jj,k) = pfGi(ii,jj,k) + sxsnowG(j,i,k)             &
             * weights(ii,jj,ll)

             sumw = sumw + weights(ii,jj,ll)
           ENDDO

           pfGi(ii,jj,k) = pfGi(ii,jj,k) / sumw

#if ( DEBUG > 1 )
      DO i=1,NINT(SQRT(REAL(nz,KIND=4)))
        WRITE(*,'(3F18.8)') (tab_dat(ii,jj,nw,j)/1000.0,j=(i-1)*3+1,i*3)
      ENDDO
      DO i=1,NINT(SQRT(REAL(nz,KIND=4)))
        WRITE(*,'(3F10.5)') (weights(ii,jj,j)/sumw,j=(i-1)*3+1,i*3)
      ENDDO
      READ(*,*)
#endif

         ENDDO
      ENDDO
      ENDDO

      interpolate = .true.

      end function interpolate

      end module interpolate_mod

!-dmr The End of All Things (op. cit.)

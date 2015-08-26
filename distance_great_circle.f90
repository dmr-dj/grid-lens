       MODULE dgc
!-----|--1--------2---------3---------4---------5---------6---------7-|
!            dgc means initially "Distange on a GReat Circle"
!     The module is intended to compute the distance between two
!       points on a sphere and includes some additional routines for
!       easy interpolation
!
!     This is version 1.0, created June, 3rd, 2013
!     Last modification: June, 13th, 2013
!
! Copyright 2013, Didier M. Roche a.k.a. dmr
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

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Some generic definitions
!-----|--1--------2---------3---------4---------5---------6---------7-|


!       0.1: Earth's radius is fixed here. It is a mean radius as 
!            generally taken, in meters

           REAL(KIND=8), PARAMETER :: earth_r = 6371.0d3

!       0.2: PI, the number!!

           REAL(KIND=8), PARAMETER :: PI = DACOS(-1.0d0)

          CONTAINS

           SUBROUTINE dgc_main(lat_un,lat_deux,lon_un,lon_deux,dist)
!-----|--1--------2---------3---------4---------5---------6---------7-|
!       dgc_main performs the actual computation of the distance 
!        between two lat,lon points on a sphere
!
!       input variables: lat_un, lat_deux and lon_un,lon_deux are the 
!             coordinates of points un et deux
!
!       output variable: dist, distance in meters between un et deux
!-----|--1--------2---------3---------4---------5---------6---------7-|

           IMPLICIT NONE

           REAL(KIND=8), INTENT(IN) :: lat_un,lat_deux,lon_un,lon_deux
           REAL(KIND=8), INTENT(OUT) :: dist

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Local variables
!-----|--1--------2---------3---------4---------5---------6---------7-|

           REAL(KIND=8) :: rad_lat_un,rad_lat_deux,rad_lon_un           &
                         , rad_lon_deux, dx, dy, dz, Ch, d_sigma


!-----|--1--------2---------3---------4---------5---------6---------7-|
!       1.1: implementation of distance on a sphere 
!-----|--1--------2---------3---------4---------5---------6---------7-|

           rad_lat_un = deg2rad(lat_un)
           rad_lat_deux = deg2rad(lat_deux)
           rad_lon_un = deg2rad(lon_un)
           rad_lon_deux = deg2rad(lon_deux)

           dx = DCOS(rad_lat_deux)*DCOS(rad_lon_deux)                   &
               -DCOS(rad_lat_un)*DCOS(rad_lon_un)

           dy = DCOS(rad_lat_deux)*DSIN(rad_lon_deux)                   &
               -DCOS(rad_lat_un)*DSIN(rad_lon_un)

           dz = DSIN(rad_lat_deux) - DSIN(rad_lat_un)

           Ch = SQRT(dx**2+dy**2+dz**2)

!       d_sigma is radius at the center of the earth between the two points
           d_sigma = 2*ASIN(Ch/2)

           dist = d_sigma*earth_r

           END SUBROUTINE dgc_main

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Conversion of an angle from degrees to radians
!-----|--1--------2---------3---------4---------5---------6---------7-|

           REAL(KIND=8) FUNCTION deg2rad(angle)

           IMPLICIT NONE

           REAL(KIND=8), INTENT(IN) :: angle

           deg2rad = 1.0d0/180.0d0 * PI * angle

           END FUNCTION deg2rad

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Expand a array from lon to -Nex+1:lon+Nex, using circular 
!        assumption
!-----|--1--------2---------3---------4---------5---------6---------7-|
           SUBROUTINE expand_tab(tab,N,Nex,tab_ex,flag)

           IMPLICIT NONE

           INTEGER, INTENT(IN) :: N, Nex, flag
           REAL(KIND=8), DIMENSION(N), INTENT(IN) :: tab
           REAL(KIND=8), DIMENSION(-Nex+1:N+Nex), INTENT(OUT) :: tab_ex

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Local variables
!-----|--1--------2---------3---------4---------5---------6---------7-|

           INTEGER :: i


           IF (flag.eq.0) THEN ! non circular (latitude)
             DO i=-Nex+1,0
               tab_ex(i) = tab(1)
             ENDDO
             DO i=N+1,N+Nex
               tab_ex(i) = tab(N)
             ENDDO
             tab_ex(1:N) = tab(1:N)
           ELSE ! circular (longitude)
             DO i=-Nex+1,0
               tab_ex(i) = tab(circular_coord(i,1,N))-360.0d0
             ENDDO
             DO i=N+1,N+Nex
               tab_ex(i) = tab(circular_coord(i,1,N))+360.0d0
             ENDDO
             tab_ex(1:N) = tab(1:N)
           ENDIF

           END SUBROUTINE expand_tab

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       From a ordered suite of numbers, compute the geometric bound
!        between two consecutive numbers. 
!-----|--1--------2---------3---------4---------5---------6---------7-|

           SUBROUTINE compute_bounds(tab,N,tab_b)

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Input: tab, ordered suite of numbers ; N size of tab
!       Output: tab_b, bounds of tab
!
!       Nota Bene: tab(i) is between tab_b(i-1) and tab_b(i)
!-----|--1--------2---------3---------4---------5---------6---------7-|
           IMPLICIT NONE

           INTEGER, INTENT(IN) :: N
           REAL(KIND=8), DIMENSION(N), INTENT(IN) :: tab
           REAL(KIND=8), DIMENSION(N), INTENT(OUT) :: tab_b

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Local variables
!-----|--1--------2---------3---------4---------5---------6---------7-|
           INTEGER :: i
           REAL(KIND=8) :: delta

           tab_b = 0.0d0

           DO i=1,N-1
! --- I assume that tab(i+1) > tab
             delta = tab(i+1)-tab(i)
             tab_b(i) = tab(i) + delta/2.0d0 
           ENDDO

           END SUBROUTINE compute_bounds

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       For a circular coordinate between lon_l and lon_u, return the
!        value valu re-ordered between the bounds
!-----|--1--------2---------3---------4---------5---------6---------7-|

           INTEGER FUNCTION circular_coord(valu,lon_l,lon_u)

           IMPLICIT NONE

           INTEGER, INTENT(IN) :: lon_l,lon_u
           INTEGER, INTENT(IN) :: valu

             IF (valu.LT.lon_l) THEN
               circular_coord = lon_u + (valu-(lon_l-1))
             ELSE IF (valu.GT.lon_u) THEN
               circular_coord = valu - lon_u + (lon_l-1)
             ELSE
               circular_coord = valu
             ENDIF


           END FUNCTION circular_coord


!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Return a longitude between [0;360] in degrees
!-----|--1--------2---------3---------4---------5---------6---------7-|
           REAL(KIND=8) FUNCTION lon_360(long)

           IMPLICIT NONE

           REAL(KIND=8), INTENT(IN) :: long


           REAL(KIND=8), PARAMETER :: zero = 0.0d0
           REAL(KIND=8), PARAMETER :: tcs = 360.0d0
           IF (long.LT.zero) THEN  
             lon_360 = tcs + long
           ELSE IF (long.GT.tcs) THEN
             lon_360 = long - tcs
           ELSE
             lon_360 = long
           ENDIF
           END FUNCTION lon_360

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Test whether an integer is within a lower and upper limits
!-----|--1--------2---------3---------4---------5---------6---------7-|

           LOGICAL FUNCTION withini(valeur,l_limit,u_limit)

             IMPLICIT NONE 

             INTEGER, INTENT(IN) :: valeur, l_limit, u_limit

             withini = .FALSE.

             IF ((valeur.GE.l_limit).AND.(valeur.LE.u_limit)) THEN
               withini=.TRUE.
             ENDIF

           END FUNCTION withini

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Find the closest cell for given long, lat in tab_lat; tab_lon
!-----|--1--------2---------3---------4---------5---------6---------7-|

           SUBROUTINE find_closest_EC_cell(long,lat,tab_lat, tab_lat_b  &
                 ,nlat, nlon, tab_lon, tab_lon_b, i_close,j_close,ext)

           IMPLICIT NONE

           REAL(KIND=8), INTENT(IN) :: lat, long
           INTEGER, INTENT(IN) :: nlat, nlon, ext
           REAL(KIND=8), DIMENSION(-ext+1:nlat+ext),INTENT(IN) ::tab_lat&
                     , tab_lat_b 
           REAL(KIND=8), DIMENSION(-ext+1:nlon+ext),INTENT(IN) ::tab_lon&
                     , tab_lon_b
           INTEGER, INTENT(OUT) :: i_close, j_close

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Local variables
!-----|--1--------2---------3---------4---------5---------6---------7-|

           REAL(KIND=8), PARAMETER :: south_border = -90.0d0,           &
                                       west_border = 360.0d0,           &
                                              zero = 0.0d0

           REAL(KIND=8) :: step_lat, step_lon 


!           REAL(KIND=8) :: long_360
           INTEGER :: j_close_m1, j_close_P1


           step_lat = 180.0d0 / nlat
           step_lon = 360.0d0 / nlon


           !-dmr Initial Guess latitude / longitude
           i_close = NINT((lat-south_border)/step_lat)


           IF (.NOT.withini(i_close,LBOUND(tab_lat,DIM=1)               &
                                   ,UBOUND(tab_lat,DIM=1)-1)) THEN
            WRITE(*,*) "PROBLEM INITIAL GUESS !!!"
            CALL FLUSH()
            READ(*,*)
           ENDIF


           j_close = NINT((lon_360(long))/step_lon + 1)


           IF (.NOT.withini(j_close,LBOUND(tab_lon,DIM=1)               &
                                   ,UBOUND(tab_lon,DIM=1)-1)) THEN
            WRITE(*,*) "PROBLEM INITIAL GUESS !!!", j_close, long
            CALL FLUSH()
            READ(*,*)
           ENDIF

           DO WHILE(lat.GT.(tab_lat(i_close+1)))
             i_close = i_close + 1
           END DO

            IF ((lat.GT.tab_lat_b(i_close)).AND.                        &
                (lat.LT.tab_lat(i_close+1))) THEN
              i_close = i_close + 1
            ELSEIF ((lat.LT.tab_lat_b(i_close-1)).AND.                  &
                    (lat.GT.tab_lat(i_close-1))) THEN
              i_close = i_close - 1
            ENDIF

           j_close_m1 = j_close - 1
           j_close_p1 = j_close + 1

           IF (ABS(long-tab_lon_b(j_close+1)).GT.2.0*step_lon) THEN
              WRITE(*,*) "PROBLEM LONGITUDE IN find_closest", j_close
              WRITE(*,*) long, tab_lon_b(j_close_p1)
              READ(*,*)
           ENDIF

           IF ( long.GE.tab_lon_b(j_close).AND.                         &
                long.LT.tab_lon(j_close_p1)    ) THEN
             j_close = j_close_p1
           ELSEIF ( (long.LT.tab_lon_b(j_close_m1)).AND.                &
                    (long.GT.tab_lon(j_close_m1))       ) THEN
             j_close = j_close_m1
           ENDIF

           END SUBROUTINE find_closest_EC_cell

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Find the corner in which long, lat stands with respect to 
!        i_close, j_close
!-----|--1--------2---------3---------4---------5---------6---------7-|

           SUBROUTINE which_corner(long,lat,tab_lat,nlat, nlon, tab_lon &
                                  ,i_close,j_close,c)

           IMPLICIT NONE

           REAL(KIND=8), INTENT(IN) :: lat, long
           INTEGER, INTENT(IN) :: nlat, nlon
           REAL(KIND=8), DIMENSION(nlat), INTENT(IN) :: tab_lat
           REAL(KIND=8), DIMENSION(nlon), INTENT(IN) :: tab_lon
           INTEGER, INTENT(IN) :: i_close, j_close
           INTEGER, INTENT(OUT) :: c

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Local variables
!-----|--1--------2---------3---------4---------5---------6---------7-|


           REAL(KIND=8) :: lon_1, lon_2

           c = 0

           lon_1 = long
           lon_2 = tab_lon(j_close)

           IF (ABS(lon_1-lon_2).GT.180.0d0) THEN
               IF (lon_1.LT.120.0d0) THEN
                 lon_1 = lon_1+360.0d0
               ELSE
                 lon_2 = lon_2+360.0d0
               ENDIF
           ENDIF

           IF (lon_1.GT.lon_2) THEN ! right part of cell
             IF (lat.GT.tab_lat(i_close)) THEN ! upper right
               c = 1
             ELSE ! lower right
               c = 2
             ENDIF
           ELSE ! left part
             IF (lat.LT.tab_lat(i_close)) THEN ! lower left
               c=3
             ELSE ! upper left
               c=4
             ENDIF
           ENDIF

           END SUBROUTINE which_corner

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Now that we now the location of cells and corners, compute the
!        actual distance between the two locations
!-----|--1--------2---------3---------4---------5---------6---------7-|

           SUBROUTINE create_interp_data(lon,lat,tab_lat,tab_lon,nlat   &
                                  ,nlon,i_close,j_close,c,i_ism         &
                                  ,j_ism,nx,ny,nw,nz,data_int)

           IMPLICIT NONE

           REAL(KIND=8), INTENT(IN) :: lat, lon
           INTEGER, INTENT(IN) :: nlat, nlon
           REAL(KIND=8), DIMENSION(nlat), INTENT(IN) :: tab_lat
           REAL(KIND=8), DIMENSION(nlon), INTENT(IN) :: tab_lon
           INTEGER, INTENT(IN) :: i_close, j_close, i_ism, j_ism, nx, ny&
                                 ,nw,nz,c
           REAL(KIND=8), DIMENSION(nx,ny,nw,nz), INTENT(OUT) :: data_int

!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Local variables
!-----|--1--------2---------3---------4---------5---------6---------7-|

           INTEGER :: i,i_l,i_u, j,j_u, j_l, c_l, i_t, j_t
           REAL(KIND=8) :: distance, fix_lon
           INTEGER :: add_num

           if (nz.EQ.4) then
             add_num = 1

           SELECT CASE(c)
              CASE(1) !upper right
                i_l = i_close
                i_u = i_close + add_num
                j_l = j_close
                j_u = j_close + add_num
              CASE(2) !lower right
                i_l = i_close - add_num
                i_u = i_close
                j_l = j_close
                j_u = j_close + add_num
              CASE(3) !lower left
                i_l = i_close - add_num
                i_u = i_close
                j_l = j_close - add_num
                j_u = j_close
              CASE(4) !upper left
                i_l = i_close
                i_u = i_close + add_num
                j_l = j_close - add_num
                j_u = j_close
           END SELECT

           else if ((nz.eq.9).OR.(nz.EQ.25).OR.(nz.EQ.49)) then
                add_num = (NINT(SQRT(REAL(nz,KIND=4)))-1)/2 ! 9 => 1 ; 25 => 2
                i_l = i_close - add_num
                i_u = i_close + add_num

                j_l = j_close - add_num
                j_u = j_close + add_num
           endif

           c_l = 1


           DO i=i_l,i_u
             DO j=j_l,j_u


              i_t = i
              j_t = j


              IF (i_t.LT.1) THEN
               i_t = 1
              ENDIF

              IF (i_t.GT.nlat) THEN
               i_t = nlat
              ENDIF

               j_t = circular_coord(j_t,1,nlon)

               fix_lon = tab_lon(j_t)

               if (lon-tab_lon(j_t).GT.180.0) then
                 if (tab_lon(j_t).EQ.0.0d0) then
                    fix_lon= tab_lon(j_t)+360.0d0
                 endif
               endif

               call dgc_main(lat,tab_lat(i_t),lon,fix_lon,distance)

               data_int(i_ism,j_ism, 1, c_l) = i_t
               data_int(i_ism,j_ism, 2, c_l) = j_t

               data_int(i_ism,j_ism, 3, c_l) = distance

               c_l = c_l + 1

             ENDDO
           ENDDO

           END SUBROUTINE create_interp_data

       END MODULE dgc

!-dmr  The End of All Things (op. cit.)

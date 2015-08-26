#define DEBUG 0
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Class defining the grid types and properties
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      module sub_grid_class
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      use grid_class, only: surface_grid

      implicit none

      type, extends(surface_grid) :: sub_grid

        type(surface_grid), pointer :: my_parent

        integer, dimension(:,:), allocatable :: inparent_lat, inparent_lon

        integer, dimension(:,:), allocatable :: nb_child_points
        real(kind=8), dimension(:,:), allocatable :: min_elevation
        real(kind=8), dimension(:,:), allocatable :: max_elevation

      end type

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Module methods to work on the sub_grids
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Function to drive initialization a flat sub_grid
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      logical function sub_grid_init(h_grid,l_grid,subg,lat_min,lat_max,lon_min,lon_max)

      use grid_utils, only: find_indx_1d
      use grid_class, only: flat_grid_init_s

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  By reference variables
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      type(surface_grid)        , intent(inout) :: h_grid
      type(surface_grid), target, intent(in)    :: l_grid
      type(sub_grid)            , intent(out), optional   :: subg
      real(kind=8)              , intent(in),  optional   :: lat_min,lat_max,lon_min,lon_max

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Local variables
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       logical      :: lats_ok_d, lats_ok_u, lons_ok_l, lons_ok_r, coord_present
       integer      :: plat_min, plat_max, plon_min, plon_max
       integer      :: clat_min, clat_max, clon_min, clon_max
       real(kind=8) :: llon_min, ad_lonmin, ad_lonmax, ad_latmin, ad_latmax
       integer      :: i,j, s_lat, s_lon
       integer      :: prec_i, prec_j
       logical      :: succeed

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  First needs to performs checks on the spatial coverage of the grid:
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  First rapid check on the globality of the grid
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      if (l_grid%is_global) then


          h_grid%is_subgrid = .true.

      else

         lats_ok_d = (minval(l_grid%b_lats(:,:)).le.minval(h_grid%b_lats(:,:)))
         lats_ok_u = (maxval(l_grid%b_lats(2,:)).ge.maxval(h_grid%b_lats(2,:)))

#if ( DEBUG == 1 )
          write(*,*) "min b_lats l_grid / h_grid: ", minval(l_grid%b_lats(:,:)), minval(h_grid%b_lats(:,:))
#endif

! Have to deal here with the fact that the longitudes might have a discrepancy in being 0 360 or -180 180
! I assume for the moment that they HAVE TO be 0 : 360 ... if not exit with failure ...
! Or should this be assumed at the time of the grid setup to ensure conformity ???

         if (.not. (    (maxval(l_grid%p_lons(:)).le.180d0)             &
                    .or.(minval(h_grid%p_lons(:)).lt.0.d0)              &
                    .or.(minval(l_grid%p_lons(:)).lt.0.d0) )) then

         lons_ok_l = (minval(l_grid%b_lons(:,:)).le.minval(h_grid%b_lons(:,:)))
         lons_ok_r = (maxval(l_grid%b_lons(:,:)).ge.maxval(h_grid%b_lons(:,:)))

#if ( DEBUG == 1 )
         write(*,*) minval(l_grid%b_lons(:,:)), maxval(l_grid%b_lons(:,:))
         write(*,*) minval(h_grid%b_lons(:,:)), maxval(h_grid%b_lons(:,:))
#endif

         if ( lats_ok_d.and.lats_ok_u.and.lons_ok_l.and.lons_ok_r ) then

#if ( DEBUG == 1 )
           write(*,*) "Continue with setting up sub_grid"
#endif
           sub_grid_init = .true.

         else ! sub_grid is not contained in the grid, no zoom possible

           sub_grid_init = .false.

         endif

!   1| if the sub-grid is used entirely, check coherence of the outer bounds
!   2| if grid is NOT used entirely, define the sub-part of the main grid being "zoomed" in the parent
!        but also define the sub-portion of the grid (rectangular necessarily) that is to be used, from the parent
!        grid.
!   NOTA: in this setup, there is an issue related to the adequation of the coverage of both grids.
!     a) if the sub-area of the hires grid is defined from the parent grid, then there is no points from the parent that
!       have a partial child infilling. Easy to deal with dynamically.
!     b) if the area is defined in the hires grid from another consideration (e.g. ice-sheet grid) then there are potential
!       specific cases at the edges.

! For this first version, I will only deal with the a) possibility and not with b). The latter needs some thinking ...

! Need second a function that ties the subgrid, or portion of sub-grid to its parent

        else ! error on the grid being -180:180

          sub_grid_init = .false.

        endif

      endif ! on the globality of the grids ...

      if ( h_grid%is_subgrid ) then

        coord_present = present(lat_min).and.present(lat_max).and.present(lon_min).and.present(lon_max)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  lookup for the zoom bounds from the coordinate given in both the parent grid system and the zoomed grid system
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        if (present(subg).and.coord_present) then

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Here we have set the zoom grid part on the parent grid, we have to do the same with the child or zoomed grid
! dmr  Beware, to avoid limits of grids issues, lat and lon bounds should be set to the parents' limit and not to the actual limits
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if ( DEBUG == 1 )
          write(*,*) "On child grid coordinates ..."
#endif

           if ( lon_min.lt.minval(h_grid%b_lons)) then
             llon_min = 360-abs(lon_min)
           else
             llon_min = lon_min
           endif

           clat_min = find_indx_1d(h_grid%p_lats,real(lat_min,kind=8))
#if ( DEBUG == 1 )
           if (real(lat_min,kind=8).GE.h_grid%b_lats(1,clat_min).and.real(lat_min,kind=8).LE.h_grid%b_lats(2,clat_min)) then
             write(*,*) "Correctly set the lat_min @ index = ", clat_min
           else
             write(*,*) "Problem with lat_min : ", lat_min, h_grid%b_lats(1,clat_min), h_grid%b_lats(2,clat_min) &
                       , h_grid%p_lats(clat_min)
           endif
#endif

           clat_max = find_indx_1d(h_grid%p_lats,real(lat_max,kind=8),strtpt=clat_min)
#if ( DEBUG == 1 )
           if (lat_max.GE.h_grid%b_lats(1,clat_max).and.lat_max.LE.h_grid%b_lats(2,clat_max)) then
             write(*,*) "Correctly set the lat_max @ index = ", clat_max
           else
             write(*,*) "Problem with lat_max : ", lat_max, h_grid%b_lats(1,clat_max), h_grid%b_lats(2,clat_max) &
                       , h_grid%p_lats(clat_max)
           endif
#endif
           clon_min = find_indx_1d(h_grid%p_lons,real(llon_min,kind=8))

#if ( DEBUG == 1 )
           if (llon_min.GE.h_grid%b_lons(1,clon_min).and.llon_min.LE.h_grid%b_lons(2,clon_min)) &
           then
             write(*,*) "Correctly set the lon_min @ index = ", clon_min, llon_min
           else
             write(*,*) "Problem with llon_min : ", llon_min, h_grid%b_lons(2,clon_min), h_grid%b_lons(1,clon_min) &
                       , h_grid%p_lons(clon_min)
           endif
#endif
           clon_max = find_indx_1d(h_grid%p_lons,real(lon_max,kind=8))
#if ( DEBUG == 1 )
           if (lon_max.GE.h_grid%b_lons(1,clon_max).and.lon_max.LE.h_grid%b_lons(2,clon_max)) then
             write(*,*) "Correctly set the lon_max @ index = ", clon_max
           else
             write(*,*) "Problem with lon_max : ", lon_max, h_grid%b_lons(1,clon_max), h_grid%b_lons(2,clon_max) &
                       , h_grid%p_lons(clon_max)
           endif
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Start by looking up the parent grid for limits of the zoomed area
!      Update 01.06.15: I have changed the order of child / parent grids.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

! dmr update 01.06.15: lon_min => ad_lonmin etc.
           ad_lonmin = h_grid%p_lons(clon_min)
           ad_lonmax = h_grid%p_lons(clon_max)
           ad_latmin = h_grid%p_lats(clat_min)
           ad_latmax = h_grid%p_lats(clat_max)

           if (ad_lonmin.lt.minval(l_grid%b_lons)) then
             llon_min = 360-abs(ad_lonmin)
           else
             llon_min = ad_lonmin
           endif

#if ( DEBUG == 1 )
          write(*,*) "setting up the zoom in coordinates lat = ", lat_min, lat_max, "and lon = ", llon_min, lon_max
          write(*,*)
          write(*,*) "On parent grid coordinates ..."
#endif

           plat_min = find_indx_1d(l_grid%p_lats,ad_latmin)

#if ( DEBUG == 1 )

           if (ad_latmin.ge.l_grid%b_lats(1,plat_min).and.ad_latmin.le.l_grid%b_lats(2,plat_min)) then
             write(*,*) "Correctly set the lat_min @ index = ", plat_min, l_grid%b_lats(:,plat_min), ad_latmin
           else
             write(*,*) "Problem with lat_min : ", ad_latmin, l_grid%b_lats(1,plat_min), l_grid%b_lats(2,plat_min) &
                       , l_grid%p_lats(plat_min)
           endif
#endif

           plat_max = find_indx_1d(l_grid%p_lats,ad_latmax,strtpt=plat_min)
#if ( DEBUG == 1 )
           if (ad_latmax.le.l_grid%b_lats(2,plat_max).and.ad_latmax.ge.l_grid%b_lats(1,plat_max)) then
             write(*,*) "Correctly set the lat_max @ index = ", plat_max, l_grid%b_lats(:,plat_max), ad_latmax
           else
             write(*,*) "Problem with lat_max : ", ad_latmax, l_grid%b_lats(1,plat_max), l_grid%b_lats(2,plat_max) &
                       , l_grid%p_lats(plat_max)
           endif
#endif

           plon_min = find_indx_1d(l_grid%p_lons,llon_min)

#if ( DEBUG == 1 )
           if (llon_min.GE.l_grid%b_lons(1,plon_min).and.llon_min.LE.l_grid%b_lons(2,plon_min)) &
           then
             write(*,*) "Correctly set the lon_min @ index = ", plon_min, l_grid%b_lons(:,plon_min), llon_min
           else
             write(*,*) "Problem with llon_min : ", llon_min, l_grid%b_lons(2,plon_min), l_grid%b_lons(1,plon_min) &
                       , l_grid%p_lons(plon_min)
           endif
#endif
           plon_max = find_indx_1d(l_grid%p_lons,ad_lonmax)
#if ( DEBUG == 1 )
           if (ad_lonmax.GE.l_grid%b_lons(1,plon_max).and.ad_lonmax.LE.l_grid%b_lons(2,plon_max)) then
             write(*,*) "Correctly set the lon_max @ index = ", plon_max, l_grid%b_lons(:,plon_max), ad_lonmax
           else
             write(*,*) "Problem with lon_max : ", ad_lonmax, l_grid%b_lons(1,plon_max), l_grid%b_lons(2,plon_max) &
                       , l_grid%p_lons(plon_max)
           endif
#endif


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   NOTA: There is no guarantee that the outside border of the h_grid is within the borders of the l_grid. The only constraint
!        here is that the center point of the h_grid is within the l_grid cell
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Now we have to initialize the sub_grid
!      1 -> Compute its size
! dmr  slat = extent of subgridding in the child grid
!      slon = same as slat for the longitude
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        s_lat = clat_max - clat_min + 1

        if (clon_min.gt.clon_max) then ! case where the minimum is before zero line ...
          s_lon = (h_grid%n_lon-clon_min)+1+clon_max
        else
          s_lon = clon_max-clon_min + 1
        endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      2 -> Initialize the flat_grid of the sub_grid
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        succeed = flat_grid_init_s(subg%surface_grid%flat_grid,(s_lat),(s_lon))


        subg%p_lats(:) = h_grid%p_lats(clat_min:clat_max)
        subg%b_lats(:,:) = h_grid%b_lats(:,clat_min:clat_max)


        if (clon_min.gt.clon_max) then ! case where the minimum is before zero line ...

          subg%p_lons(1:(h_grid%n_lon-clon_min+1)) = h_grid%p_lons(clon_min:h_grid%n_lon)
          subg%p_lons((h_grid%n_lon-clon_min+2):subg%n_lon) = h_grid%p_lons(1:clon_max)

          subg%b_lons(:,1:(h_grid%n_lon-clon_min+1)) = h_grid%b_lons(:,clon_min:h_grid%n_lon)
          subg%b_lons(:,(h_grid%n_lon-clon_min+2):subg%n_lon) = h_grid%b_lons(:,1:clon_max)

        else

          subg%p_lons(:) = h_grid%p_lons(clon_min:clon_max)
          subg%b_lons(:,:) = h_grid%b_lons(:,clon_min:clon_max)

        endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      3 -> Initialize the surface_grid of the sub_grid (i.e. copy the given variables of h_grid on its portion subg
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        if (h_grid%ini_elevation) then

          allocate(subg%s_elevation(subg%n_lon,subg%n_lat))

          if (clon_min.gt.clon_max) then ! case where the minimum is before zero line ...

            subg%s_elevation(1:(h_grid%n_lon-clon_min+1),:)             &
              = h_grid%s_elevation(clon_min:h_grid%n_lon,clat_min:clat_max)

            subg%s_elevation((h_grid%n_lon-clon_min+2):subg%n_lon,:)    &
                                  = h_grid%s_elevation(1:clon_max,clat_min:clat_max)

          else

            subg%s_elevation(:,:) = h_grid%s_elevation(clon_min:clon_max,clat_min:clat_max)

          endif

        endif ! ini_elevation

        elseif (present(subg)) then ! the sub_grid limits are used
          write(*,*) "Setting up the sub_grid from its coordinate"
          write(*,*) "Not implemented yet!!"
          read(*,*)
        else ! no subg provided, then the "subgrid" is global as is the parent grid
          write(*,*) "Setting up the sub_grid globally"
          write(*,*) "Not implemented yet!!"
          read(*,*)
        endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      4 -> Cluster the subg points within the grid points of its parent
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       subg%my_parent => l_grid

       allocate(subg%inparent_lat(subg%n_lon,subg%n_lat))
       allocate(subg%inparent_lon(subg%n_lon,subg%n_lat))

       allocate(subg%nb_child_points(subg%my_parent%n_lon,subg%my_parent%n_lat))
       allocate(subg%min_elevation(subg%my_parent%n_lon,subg%my_parent%n_lat))
       allocate(subg%max_elevation(subg%my_parent%n_lon,subg%my_parent%n_lat))

       subg%nb_child_points(:,:) = 0

       subg%min_elevation(:,:) =  20000d0
       subg%max_elevation(:,:) = -20000.0d0


!       prec_i = 1
       do i=1,subg%n_lon
!         prec_j = 1
         do j=1,subg%n_lat

           subg%inparent_lat(i,j) = find_indx_1d(subg%my_parent%p_lats,real(subg%p_lats(j),kind=8)) ! ,strtpt=prec_j)
           subg%inparent_lon(i,j) = find_indx_1d(subg%my_parent%p_lons,real(subg%p_lons(i),kind=8)) ! ,strtpt=prec_i)

           prec_j = subg%inparent_lat(i,j)
           prec_i = subg%inparent_lon(i,j)

           if ((prec_i.eq.subg%my_parent%n_lon).and.(subg%p_lons(i).GT.maxval(subg%my_parent%b_lons))) then
             prec_i = 1
             subg%inparent_lon(i,j) = 1
           endif

           subg%nb_child_points(prec_i,prec_j) = subg%nb_child_points(prec_i,prec_j) + 1
           subg%min_elevation(prec_i,prec_j) = min(subg%min_elevation(prec_i,prec_j),subg%s_elevation(i,j))
           subg%max_elevation(prec_i,prec_j) = max(subg%max_elevation(prec_i,prec_j),subg%s_elevation(i,j))
!           write(*,*) subg%min_elevation(prec_i,prec_j), subg%max_elevation(prec_i,prec_j)

         enddo
       enddo

#if ( DEBUG == 1 )
      do j=plat_max+1, plat_min-1, -1
        write(*,*) "NB points", j , (subg%nb_child_points(i,j),i=plon_min-1,subg%my_parent%n_lon)   &
                              , (subg%nb_child_points(i,j),i=1,plon_max+1)
      enddo

      write(*,*) "Clustering of zoomed grid achieved: ", maxval(subg%nb_child_points), minval(subg%nb_child_points)
#endif

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!      5 -> Initialize interpolation method, using expanded grid
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      endif ! on h_grid is sub_grid
      end function sub_grid_init

      end module sub_grid_class
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

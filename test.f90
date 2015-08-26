      program test

#define DGC_USE 1
#define NCIO_USE 1

      use grid_class, only: surface_grid, surface_grid_init
      use sub_grid_class, only: sub_grid, sub_grid_init

#if (DGC_USE == 1)
      use interpolate_mod, only: interpolate_init, interpolate
#endif

#if (NCIO_USE == 1)
       USE ncio, only: nc_create, nc_write_attr, nc_write_dim, nc_write
#endif

      implicit none

      type(surface_grid), target :: lres_land_grid, hres_land_grid
      type(sub_grid)     :: zoom_grid
      logical :: succeed
      character(len=256) :: input_file, input_file_l
      character(len=256), parameter     :: directory = "inputdata/"


#if (DGC_USE == 1)
      real(kind=8), allocatable, dimension(:,:) :: YLAT, XLONG
      real(kind=8), allocatable, dimension(:,:,:,:) :: tab_dat
      real(kind=8), allocatable, dimension(:,:,:) :: interpolatable, interpolated
      integer, parameter :: nw = 3, nz = 9, ex = 2, nbmois = 1
      integer :: i,j
      logical :: results

#if (NCIO_USE == 1)
!-----|--1--------2---------3---------4---------5---------6---------7-|
!       Added for output with NCIO
!-----|--1--------2---------3---------4---------5---------6---------7-|
       character(len=256) :: filename
#endif

#endif

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
!      lres_land_grid%lat_name = "lat"
!      lres_land_grid%latbnd_name = "bounds_lat"
!      lres_land_grid%lon_name = "lon"
!      lres_land_grid%lonbnd_name = "bounds_lon"
      lres_land_grid%elevation_name = "topo"
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
      sub_grid_init(hres_land_grid,lres_land_grid,subg=zoom_grid,lat_min=35.0d0,lat_max=80.0d0,lon_min=340.0d0,lon_max=40.0d0)

!       YLAT = DIMENSION(zoom_grid%nlon,zoom_grid%nlat)
!       XLONG = DIMENSION(zoom_grid%nlon,zoom_grid%nlat)
!       INTEGER, PARAMETER :: nw = 3, nz = 9, ex = 2
!       tab_dat == REAL(KIND=8), DIMENSION(zoom_grid%nlon,zoom_grid%nlat,nw,nz) :: tab_dat

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
      allocate(interpolated(zoom_grid%n_lon,zoom_grid%n_lat,nbmois))

      interpolatable(:,:,1) = lres_land_grid%s_elevation

      results = interpolate(tab_dat,interpolatable,interpolated,zoom_grid%n_lon,zoom_grid%n_lat,nw,nz,lres_land_grid%n_lon &
                           ,lres_land_grid%n_lat,nbmois)

#endif

      write(*,*) "Globality of grids: ", lres_land_grid%is_global, hres_land_grid%is_global
      write(*,*) "Setup of sub_grid = ", succeed
      write(*,*) "Interpolation results", results

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
      call nc_write(filename,"hres_topo",zoom_grid%s_elevation(:,:),dim1="x",dim2="y")

#endif
      end program test

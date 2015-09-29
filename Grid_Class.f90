#define DEBUG 0
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Class defining the grid types and properties
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      module grid_class
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      implicit none

      integer, parameter, public                   :: str_len =256

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Type definitions for the grids
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      type flat_grid

         integer                                   :: n_lat        ! size in latitude
         integer                                   :: n_lon        ! size in longitude
         integer                                   :: n_bounds = 2 ! assuming two points for each bound (N and S)

         real(kind=8), dimension(:), allocatable   :: p_lats       ! center point latitude
         real(kind=8), dimension(:,:), allocatable :: b_lats       ! latitude bounds

         real(kind=8), dimension(:), allocatable   :: p_lons       ! center point longitude
         real(kind=8), dimension(:,:), allocatable :: b_lons       ! longitude bounds

         real(kind=8), dimension(:,:), allocatable :: g_surf       ! surface of the grid cells

         logical                           :: is_global = .false.  ! global grid ?

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr    These name definitions are the default for the .nc file defining the grid
! dmr      they are to be modified before the call to the init functions if
! dmr      a non std name exists in the .nc file
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

         character(len=str_len) ::                                       &
                                          lat_name   ="latitude",        &
                                          lon_name   ="longitude",       &
                                          latbnd_name="bounds_latitude", &
                                          lonbnd_name="bounds_longitude",&
                                          g_surf_name="gridbox_area"

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   has the grid been initialized?
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

         logical                                   :: is_defined=.false.

      end type flat_grid

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Variable on a surface_grid
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      type surface_var

         logical                  :: ini_var = .false.
         character(len=str_len)   :: var_name = " ", f_name = "  "
         real(kind=8), dimension(:,:), allocatable :: var_data

         contains

         procedure :: var_init => set_fields_surf_var

      end type surface_var

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Surface_grid is an extension of the flat_grid type, adding the sub_grid capability and the variables
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      type, extends (flat_grid) :: surface_grid

         logical                :: is_subgrid = .false.

         integer                :: nb_vars = 1

         type(surface_var), dimension(1) :: surf_grid_vars

! dmr For now, I need this conventional indexing of the different variables, though not limiting.
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Topographic elevation of the surface grid, mandatory
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
         integer                 :: indx_elevation = 1
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   If not a subgrid, then you may have fractional land cover in the grid possibly
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
         integer                 :: indx_landfrac  = 2
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   On land grids, you certainly may want to drain water at the surface
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
         integer                 :: indx_drainage  = 3

!        logical                :: ini_elevation = .true.
!         character(len=str_len) :: elevation_name = "topo"
!         real(kind=8), dimension(:,:), allocatable :: s_elevation
!
!         logical                :: ini_landfrac = .false.
!         character(len=str_len) :: landfrac_name = "landfrac"
!         real(kind=8), dimension(:,:), allocatable :: s_landfrac
!
!        logical                :: ini_drainage = .false.
!        character(len=str_len) :: drainage_name = "drainage"
!        real(kind=8), dimension(:,:), allocatable :: s_drainage

        contains

          procedure :: set_nvars => set_number_vars

      end type surface_grid

      contains

      function set_number_vars(this,int_val) result(set_nvars)

        class(surface_grid), intent(inout) :: this
        integer            , intent(in)    :: int_val

        logical                            :: set_nvars

        this%nb_vars = int_val
!        allocate(this%surf_grid_vars(1:this%nb_vars))
        set_nvars = .true.

      end function set_number_vars

      function set_fields_surf_var(this,ini,namee) result(var_init)

        class(surface_var), intent(inout) :: this
        character(len=str_len), intent(in) :: namee
        logical,                intent(in) :: ini

        logical                            :: var_init

        this%ini_var  = ini
        this%var_name = namee
        var_init      = .true.

      end function set_fields_surf_var


!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Module methods to work on the grids
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!      logical function surface_grid_init(f_grid,f_name)

!      logical function flat_grid_init(f_grid,f_name)

!      logical function dims_init(f_grid,f_name)

!      logical function read_metrics(grid,f_name)

!      logical function read_OneVar(var_name,var_array,size_x,size_y,f_name)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Function to initialize a flat grid
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      logical function surface_grid_init(s_grid,fs_name,fg_name)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        type(surface_grid), intent(inout)            :: s_grid  ! the surface grid variable
        character(len=str_len), intent(in)           :: fs_name ! netCDF file name to read in the surface data from
        character(len=str_len), intent(in), optional :: fg_name ! netCDF file name to read in the grid data from, if different from fs

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        logical :: base_grid_ok = .false. , var_ok = .false.

        integer :: i   !,j

        base_grid_ok = s_grid%is_defined

        if ( .not. base_grid_ok ) then ! basic grid is not defined ... please init!
          if ( present(fg_name) ) then ! use fg_name for the metrics
            base_grid_ok = flat_grid_init(s_grid,fg_name)
          else ! assume that the metrics are in the same file as the other characteristics
            base_grid_ok = flat_grid_init(s_grid,fs_name)
          endif
        endif

        if ( base_grid_ok ) then

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Allocate the arrays of type surface_variable
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!        allocate(s_grid%surf_grid_vars(s_grid%nb_vars))
! dmr 2015-09-29 -> moved to the type-bound function of surface_var

         do i=1, s_grid%nb_vars
            if (s_grid%surf_grid_vars(i)%ini_var) then ! need to initialize that variable on the grid

              allocate(s_grid%surf_grid_vars(i)%var_data(1:s_grid%n_lon,1:s_grid%n_lat))

              var_ok = read_OneVar(s_grid%surf_grid_vars(i)%var_name,s_grid%surf_grid_vars(i)%var_data,                         &
                                   s_grid%n_lon,s_grid%n_lat,fs_name)
            endif
         enddo

!         if ( s_grid%ini_elevation ) then ! do something for that variable
!
!           allocate(s_grid%s_elevation(s_grid%n_lon,s_grid%n_lat))
!
!           var_ok = read_OneVar(s_grid%elevation_name,s_grid%s_elevation&
!                 ,s_grid%n_lon,s_grid%n_lat,fs_name)
!
!         endif

! -dmr ### This piece of code need to be updated to take into account the new classes
! -dmr ### for now, it is commented before doing better

!         if ( s_grid%is_subgrid ) then ! this is a subgrid, no landfrac variable to be read
!
!           s_grid%ini_landfrac = .true.
!
!           allocate(s_grid%s_landfrac(s_grid%n_lon,s_grid%n_lat))
!
!           do j=LBOUND(s_grid%s_elevation,1), UBOUND(s_grid%s_elevation,1)
!             do i=LBOUND(s_grid%s_elevation,2), UBOUND(s_grid%s_elevation,2)
!               if (s_grid%s_elevation(j,i).GT.0.0d0) then
!                 s_grid%s_landfrac(j,i) = 1.0d0
!               else
!                 s_grid%s_landfrac(j,i) = 0.0d0
!               endif
!             enddo
!           enddo
!
!         else ! not a subgrid
!
!           if ( s_grid%ini_landfrac ) then ! do something for that variable
!
!             allocate(s_grid%s_landfrac(s_grid%n_lon,s_grid%n_lat))
!
!             var_ok = read_OneVar(s_grid%landfrac_name,s_grid%s_landfrac  &
!                 ,s_grid%n_lon,s_grid%n_lat,fs_name)
!
!           endif
!         endif
#if ( DEBUG == 1 )
         else ! base_grid_ok false = could not initialize base grid metrics
          write(*,*) "Cannot initialize the base grid in s_grid ... !!"
#endif

        endif ! base_grid_ok

        surface_grid_init = var_ok .and. base_grid_ok

      end function surface_grid_init

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Function to initialize a flat grid
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      logical function flat_grid_init(f_grid,f_name)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        class(flat_grid), intent(inout)     :: f_grid ! the flat grid variable
        character(len=str_len), intent(in)  :: f_name ! netCDF file name to read in the grid from

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
        logical :: file_exists = .false., dims_ok = .false., read_dims = .false., read_g_surf = .false.

        inquire(file=trim(f_name), exist=file_exists)

        if (file_exists) then

          dims_ok = dims_init(f_grid,f_name)

         if (dims_ok) then

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Allocation of flat grid specifics (metrics and points)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

          allocate(f_grid%p_lats(f_grid%n_lat))
          allocate(f_grid%b_lats(f_grid%n_bounds,f_grid%n_lat))

          allocate(f_grid%p_lons(f_grid%n_lon))
          allocate(f_grid%b_lons(f_grid%n_bounds,f_grid%n_lon))

          allocate(f_grid%g_surf(f_grid%n_lon,f_grid%n_lat))

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Reading input coordinates
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

          read_dims = read_metrics(f_grid,f_name)

         if (read_dims) then ! I have finished setting up the metrics of the grid ...
           f_grid%is_defined = .true.

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Now need to input the surface area of the boxes into the surface grid
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
           read_g_surf = read_OneVar(f_grid%g_surf_name,f_grid%g_surf,f_grid%n_lon,f_grid%n_lat,f_name)

#if ( DEBUG == 1 )
           write(*,*) "For the surface area variable success = ", read_g_surf
           write(*,*) "Test value gives = ", f_grid%g_surf(15,15)
#endif
#if ( DEBUG == 1 )
         else ! read_dims F
          write(*,*) "Cannot read metrics for the flat grid"
#endif
         endif
#if ( DEBUG == 1 )
         else ! dims_ok F
          write(*,*) "Dimensions of flat_grid read to zero ... "
          write(*,*) "Cannot allocate void array !!"
#endif
         endif
#if ( DEBUG == 1 )
        else ! file_exist F
          write(*,*) "Cannot open file = ", trim(f_name)
#endif
        endif ! file_exist

        flat_grid_init = file_exists .and. dims_ok .and. read_dims

        if ( flat_grid_init ) then
          call set_global_flag(f_grid)
        endif

      end function flat_grid_init

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Function to initialize a flat grid without a file name
!      This init assumes that the p_lats, p_lons, b_lats, b_lons
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      logical function flat_grid_init_s(f_grid,size_lat,size_lon)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   By reference variables ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        type(flat_grid), intent(inout) :: f_grid ! the flat grid variable
        integer        , intent(in)    :: size_lat, size_lon

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

       f_grid%n_lat = size_lat
       f_grid%n_lon = size_lon

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Allocation of flat grid specifics (metrics and points)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

        allocate(f_grid%p_lats(f_grid%n_lat))
        allocate(f_grid%b_lats(f_grid%n_bounds,f_grid%n_lat))

        allocate(f_grid%p_lons(f_grid%n_lon))
        allocate(f_grid%b_lons(f_grid%n_bounds,f_grid%n_lon))

        allocate(f_grid%g_surf(f_grid%n_lon,f_grid%n_lat))

        f_grid%is_defined = .true.

        flat_grid_init_s = f_grid%is_defined

        if ( flat_grid_init_s ) then
          call set_global_flag(f_grid)
        endif

      end function flat_grid_init_s

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Compute the grid extension, is it global?
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      subroutine set_global_flag(fgrid)

      type(flat_grid), intent(inout) :: fgrid

      real(kind=8) :: petit, bord_S, bord_N, bord_E, bord_O, max_lon, min_lon
      real(kind=8), parameter :: pi = 180.0d0

      petit = epsilon(fgrid%b_lats(1,1))

      bord_S = pi / (-2.0d0)
      bord_N = pi / 2.0d0
      bord_O = 0.0d0
      bord_E = 2.0d0 * pi

#if ( DEBUG == 1 )
! workaround to force a global grid temporarily in latitude
!      fgrid%b_lats(1,1) = 90.0d0
!      fgrid%b_lats(2,fgrid%n_lat) = -90.0d0
#endif


      if ( (minval(fgrid%b_lats(:,:)).LE.(bord_S+petit)) .and. (maxval(fgrid%b_lats(:,:)).GE. (bord_N-petit)) ) then

#if ( DEBUG == 1 )
        write(*,*) "given grid IS global in latitude"
#endif
        min_lon = minval(fgrid%b_lons(:,:))
        max_lon = maxval(fgrid%b_lons(:,:))

        if ( ( min_lon.LE.(bord_O+petit)) .and. (max_lon.GE.(bord_E-petit)) ) then
        ! bounds of the given grid are 0:360, easy case

          fgrid%is_global = .true.

#if ( DEBUG == 1 )
          write(*,*) "given grid IS global in longitude"
#endif
        elseif (abs(min_lon+2.0d0*pi-max_lon).LE.petit) then
        ! works only if the lower_bound = upper_bound [2*pi]

#if ( DEBUG == 1 )
          write(*,*) "given grid IS global in longitude but not 0:360"
#endif

          fgrid%is_global = .true.

        else

#if ( DEBUG == 1 )
          write(*,*) "given grid is NOT global in longitude"
#endif
        fgrid%is_global = .false.

        endif

      else

#if ( DEBUG == 1 )
       write(*,*) "given grid is NOT global in latitude"
#endif

       fgrid%is_global = .false.

      endif


      end subroutine set_global_flag
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Find dimensions of the grid in the netCDF read for definition
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
      logical function dims_init(f_grid,f_name)
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      use ncio, only: nc_size

      implicit none

      character(len=str_len), intent(in)    :: f_name
      type(flat_grid)       , intent(inout) :: f_grid

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#if ( DEBUG == 1 )
      write(*,*) "File used: ", trim(f_name)
      write(*,*) "Variable sought: ", trim(f_grid%lat_name)
      write(*,*) "Variable sought: ", trim(f_grid%lon_name)
#endif

      f_grid%n_lat = nc_size(f_name,f_grid%lat_name)
      f_grid%n_lon = nc_size(f_name,f_grid%lon_name)

#if ( DEBUG == 1 )
      write(*,*) "Size of array: ", f_grid%n_lat, f_grid%n_lon
#endif

      if ((f_grid%n_lat.GT.0) .and. (f_grid%n_lon.GT.0)) then
! dmr Beware, error handling only partially done for this function !!
        dims_init = .true.
      else
        dims_init = .false.
      endif

      end function dims_init

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   read dimension metrics for flat Grid ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      logical function read_metrics(grid,f_name)

      use ncio, only: nc_read

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  By reference variables
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      type(flat_grid)       , intent(inout) :: grid
      character(len=str_len), intent(in) :: f_name

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      call nc_read(f_name,grid%lat_name   ,grid%p_lats)
      call nc_read(f_name,grid%lon_name   ,grid%p_lons)
      call nc_read(f_name,grid%latbnd_name,grid%b_lats)
      call nc_read(f_name,grid%lonbnd_name,grid%b_lons)

! dmr Beware that the error handling is not set in this function. Potential problems ahead!
      read_metrics = .true.

      end function read_metrics

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   read given variable in the ad hoc file ...
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      logical function read_OneVar(var_name,var_array,size_x,size_y,f_name)

      use ncio, only: nc_read

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  By reference variables
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      character(len=str_len)                , intent(in) :: var_name
      integer                               , intent(in) :: size_x, size_y
      real(kind=8), dimension(size_x,size_y), intent(out):: var_array
      character(len=str_len)                , intent(in) :: f_name

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      call nc_read(f_name,var_name,var_array)

! dmr  Error handling not done
      read_OneVar = .true.

      end function read_OneVar

      end module grid_class
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

#define DEBUG 0
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Module defining a series of useful procedures for Grid Handling
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      module grid_utils

      implicit none

      contains

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr   Module method (function) to find the closest indx of given value in array
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      integer function find_indx_1d(tableau,valeur,strtpt)

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  By reference variable
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      real(kind=8),                 intent(in) :: valeur
      real(kind=8), dimension(:),   intent(in) :: tableau
      integer, optional,            intent(in) :: strtpt ! starting point in case we do not want to scan the whole array

!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|
! dmr  Local variables
!-----|--1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3-|

      integer      :: i, im
      real(kind=8) :: left, right

      if (present(strtpt)) then
        i = strtpt
        if (i <= lbound(tableau,1)+1 ) then
          i = lbound(tableau,1)+1
        endif
      else
        i = lbound(tableau,1)+1
      endif

      im = i -1

      do while ((i.ne.im).and.(i.le.ubound(tableau,1)))

       left = tableau(im)
       right = tableau(i)

       im = i

       i = i + abs(transfer(abs(valeur-left).GT.abs(valeur-right),i))

      enddo

      find_indx_1d = i-1

#if ( DEBUG == 1 )
       left = tableau(i-1)
       right = tableau(i)
       write(*,*) "Found the right index: ", i,im,valeur,"within : ", left,right
       write(*,*) "Bounds ...", lbound(tableau,1), ubound(tableau,1)
       read(*,*)
#endif

      return
      end function find_indx_1d

      end module grid_utils

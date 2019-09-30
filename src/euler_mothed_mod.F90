module euler_method_mod

  !! Euler method:
  !!   solution: / yn+1 = yn + f(xn + yn) * h       (n=1,2,3....n)
  !!             \  y0  = y(a)
  !!   error   : O(h2)

  implicit none

  private

  public :: euler_method

contains

  function euler_method(init_y, init_x, end_x, nstep, f) result(y)

    real     , intent(in)  :: init_y
    real     , intent(in)  :: init_x
    real     , intent(in)  :: end_x
    integer  , intent(in)  :: nstep

    procedure, pointer     :: f

    real                   :: y(nstep+1)

    integer                :: n
    real                   :: h
    real     , allocatable :: x(nstep+1)

    interface
      function f(x, y)
        real :: x
        real :: y
      end function
    end interface

    h = (init_y - end_x)/nstep

    x(0) = init_x
    y(0) = init_y

    do n = 0, nstep-1
      y(n+1) = y(n) + f(x(n), y(n)) * h
      x(n+1) = x(n) + h
    end do

  end function euler_method

end module euler_method_mod
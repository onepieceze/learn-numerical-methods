module improved_euler_method_mod

  !! Improved Euler method:
  !!
  !!   equation:  ____
  !!            / yn+1 = yn + f(xn, yn)h
  !!            |                                     ____
  !!            \ yn+1 = yn + h/2[f(xn, yn) + f(xn+1, yn+1)]
  !!
  !!   error   : O(h3)

  implicit none

  private

  public :: improved_euler_method

contains

  function improved_euler_method(init_y, init_x, end_x, nstep, f) result(y)

    real   , intent(in) :: init_y
    real   , intent(in) :: init_x
    real   , intent(in) :: end_x
    integer, intent(in) :: nstep
    real   , external   :: f

    integer             :: n
    real                :: h
    real                :: y(0:nstep)
    real                :: x(0:nstep)
    real                :: y_ba(0:nstep)

    h = (end_x - init_x)/nstep

    x(0) = init_x
    y(0) = init_y

    do n=0, nstep-1
      x(n+1) = x(n) + h
      y_ba(n+1) = y(n) + h * f(x(n), y(n))
      y(n+1) = y(n) + h / 2 * ((f(x(n), y(n)) + f(x(n+1), y_ba(n+1))))
    end do
      
  end function improved_euler_method

end module improved_euler_method_mod

module runge_kutta_O2_mod

  !! Second order Runge Kutta:
  !!
  !!   equation: / yn+1 = yn + h(w1k1 + w2k2)
  !!             |
  !!             | k1 = f(xn, yn) = y'(xn)
  !!             |
  !!             \ k2 = f(xn + ph, yn + phk1)
  !!
  !!     where: / w1 + w2 = 1
  !!            |
  !!            \ w2 * p = 1/2
  !!
  !!   erroe: O(h3)

  implicit none

  private

  public :: runge_kutta_O2

contains

  function runge_kutta_O2(init_y, init_x, end_x, nstep, f) result(y)

    real   , intent(in) :: init_y
    real   , intent(in) :: init_x
    real   , intent(in) :: end_x
    integer, intent(in) :: nstep
    real   , external   :: f

    integer             :: n
    real                :: h
    real                :: k1, k2
    real                :: w1, w2, p
    real                :: x(0:nstep)
    real                :: y(0:nstep)

    h = (end_x - init_x) / nstep

    w1 = 0.5
    w2 = 0.5
    p  = 1.0

    do n=0, nstep-1
      k1 = f(x(n), y(n))
      k2 = f(x(n)+p*h, y(n) + p*h*k1)
      y(n+1) = y(n) + h*(w1*k1 + w2*k2)
      x(n+1) = x(n) + h
    end do

  end function runge_kutta_O2

end module runge_kutta_O2_mod

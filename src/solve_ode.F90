program solve_ode

  use euler_method_mod

  implicit none

  real           :: x(11)
  real           :: y(11)
  real, external :: func

  x = ([0 + i * 0.1], i=0, i=10)

  y = euler_method(init_y=1, init_x=0, end_x=1, nstep=10, func)

  write(6, '("Euler method:")')
  write(6, '("x = ", 11F9.5)') x
  write(6, '("y = ", 11F9.5)') y

end program solve_ode

pure function func(x, y) result(res)

  real, intent(in) :: x
  real, intent(in) :: y

  real             :: re

  res = y - (x / y)

end function


  
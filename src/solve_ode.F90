program solve_ode

  use euler_method_mod
  use improved_euler_method_mod

  implicit none

  integer        :: i
  real           :: x(11)
  real           :: y(11)
  real, external :: func

  x = [ (i*0.1, i=0, 10) ]

  write(6, '("x = ", 11F9.5)') x

  y = euler_method(init_y=1., init_x=0., end_x=1., nstep=10, f=func)

  write(6, '("Euler method:")')
  write(6, '("y = ", 11F9.5)') y

  y = improved_euler_method(init_y=1., init_x=0., end_x=1., nstep=10, f=func)

  write(6, '("Improved Euler method:")')
  write(6, '("y = ", 11F9.5)') y


end program solve_ode

pure function func(x, y) result(res)

  real, intent(in) :: x
  real, intent(in) :: y

  real             :: re

  res = y - (2 * x / y)

end function


  

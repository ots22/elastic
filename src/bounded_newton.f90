module m_bounded_newton
  integer, parameter :: nr_newton_step    = -2
  integer, parameter :: nr_bisection_step = -1
  integer, parameter :: nr_converged      =  0
  integer, parameter :: nr_maxit          =  1
  integer, parameter :: nr_failed         =  2

  
contains
  subroutine bounded_newton(f,df,a,b,x,xtol_arg,ftol_arg,maxit_arg,err)
    interface
       function f(x) result(result)
         real, intent(in) :: x
         real result
       end function f
       function df(x) result(result)
         real, intent(in) :: x
         real result
       end function df
    end interface
    real, intent(inout) :: a,b,x
    real, intent(in), optional :: xtol_arg, ftol_arg
    integer, intent(in), optional :: maxit_arg
    integer, intent(out) :: err
    real xtol, ftol
    integer maxit
    real fa, fb, gsign
    integer i

    xtol = 1.0E-5
    ftol = 1.0E-5
    maxit = 1000
    if (present(xtol_arg))  xtol = xtol_arg
    if (present(ftol_arg))  ftol = ftol_arg
    if (present(maxit_arg)) maxit = maxit_arg

    fa = f(a)
    fb = f(b)
    gsign = fb - fa

    if ((fa.gt.0.and.fb.gt.0).or.(fa.lt.0.and.fb.lt.0)) then
       err = NR_FAILED
       return
    end if
    
    if (x.lt.a.or.x.gt.b) x = 0.5*(a+b)

    do i=1,maxit
       call bounded_newton_step
       if (err.eq.NR_CONVERGED) return
    end do
    err = NR_MAXIT

  contains
    subroutine bounded_newton_step
      real f0, df0, dx, x1
      gsign = merge(1.0,-1.0,gsign.ge.0.0)
      f0 = f(x)
      df0 = df(x)
      dx = -f0/df0
      x1 = x + dx

      if (abs(dx).lt.xtol.and.abs(df0*dx).lt.ftol) then
         x = x1
         err = NR_CONVERGED
      else if (x1.lt.b.and.x1.gt.a) then
         x = x1
         err = NR_NEWTON_STEP
      else
         if (gsign*f0.lt.0) then 
            a = x
         else
            b = x
         end if
         x = 0.5 * (a+b)
         err = NR_BISECTION_STEP
      end if
    end subroutine bounded_newton_step
  end subroutine bounded_newton
end module m_bounded_newton

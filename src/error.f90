module m_error
contains
  subroutine warn(msg)
    character(len=*), intent(in) :: msg
    write (0,*) "warning: ", msg
  end subroutine warn

  subroutine panic(msg)
    character(len=*), intent(in) :: msg
    write (0,*) msg
    error stop
  end subroutine panic

  subroutine assert(expr,msg)
    logical, intent(in) :: expr
    character(len=*), intent(in) :: msg
    if (.not.expr) call panic("Assert failure: " // msg)
    continue
  end subroutine assert
end module m_error

subroutine no_op(a,b,c,d,e,f,g,h)
  integer,optional,intent(in) :: a,b,c,d,e,f,g,h
  if (.false.) then
     if (present(a)) continue
     if (present(b)) continue
     if (present(c)) continue
     if (present(d)) continue
     if (present(e)) continue
     if (present(f)) continue
     if (present(g)) continue
     if (present(h)) continue
  end if
end subroutine no_op

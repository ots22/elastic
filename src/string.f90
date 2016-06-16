module m_string
contains
  function to_upper(in) result(out)
    character(*), intent(in) :: in
    character(len(in)) :: out
    integer i,c
    out = in
    do i = 1,len(in)
       c = iachar(in(i:i))
       select case (c)
       case (iachar('a'):iachar('z'))
          out(i:i) = achar(c-32)
       end select
    end do
  end function to_upper

  function to_lower(in) result(out)
    character(*), intent(in) :: in
    character(len(in)) :: out
    integer i,c
    out = in
    do i = 1,len(in)
       c = iachar(in(i:i))
       select case (c)
       case (iachar('A'):iachar('Z'))
          out(i:i) = achar(c+32)
       end select
    end do
  end function to_lower
end module m_string

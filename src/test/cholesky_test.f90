program cholesky_test
  use m_matutil
! example taken from wikipedia
  
  real A(3,3), U(3,3), expected(3,3)

  A(1,:) = [4, 12, -16]
  A(2,:) = [12, 37, -43]
  A(3,:) = [-16, -43, 98]

  U = cholesky(A)

  expected(1,:) = [2, 6, -8]
  expected(2,:) = [0, 1, 5]
  expected(3,:) = [0, 0, 3]

  if (.not.all(aprx_eq(U,expected))) STOP 1  
end program cholesky_test

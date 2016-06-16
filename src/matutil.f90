module m_matutil
  interface operator(.outer.)
     procedure outer_prod
  end interface

  real, parameter :: pi=4*atan(1.0)

contains
  pure function identity(n) result(id)
    integer, intent(in) :: n
    real id(n,n)
    integer i
    id = 0
    do i=1,n
       id(i,i) = 1
    end do
  end function identity

! matrix describing a rotation about the z-axis, through an angle `a'
  pure function rot_z(a) result(R)
    real, intent(in) :: a
    real R(3,3)
    R = reshape([cos(a),-sin(a), 0.0, &
         &       sin(a), cos(a), 0.0, &
         &       0.0,      0.0,      1.0],    [3,3])
  end function rot_z

  elemental function sgn(a)
    real, intent(in) :: a
    real sgn
    sgn = sign(1.0,a)
    if (a.eq.0.0) sgn = 0
  end function sgn

  elemental function minmod(x,y) result(mm)
    real, intent(in) ::  x,y
    real mm
    mm = (sign(0.5,x) + sign(0.5,y)) * min(abs(x), abs(y))
  end function minmod

  elemental function median(x,y,z)
    real, intent(in) :: x,y,z
    real median
    median = max(min(x,y), min(max(x,y),z));
  end function median
  
  pure function outer_prod(a,b) result(c)
    real, intent(in) :: a(:), b(:)
    real c(size(a),size(b))
    integer i,j
    do j=1,size(b)
       do i=1,size(a)
          c(i,j) = a(i) * b(j)
       end do
    end do
  end function outer_prod

  elemental function aprx_eq(a,b)
    logical aprx_eq
    real, parameter :: eps=1e-9
    real, intent(in) :: a,b
    real tol
    tol = eps
    aprx_eq = abs(b-a) < tol
  end function aprx_eq

  pure function is_symmetric(M) result(r)
    real, intent(in) :: M(:,:)
    logical r
    r = all(M.eq.transpose(M)) 
  end function is_symmetric

  pure function is_approx_symmetric(M) result(r)
    real, intent(in) :: M(:,:)
    logical r
    r = all(aprx_eq(M,transpose(M)))
  end function is_approx_symmetric

  pure function matrix_to_voigt(M) result(V)
    real, intent(in) :: M(3,3)
    real V(6)
    V(1) = M(1,1)
    V(2) = M(2,2)
    V(3) = M(3,3)
    V(4) = M(2,3)
    V(5) = M(1,3)
    V(6) = M(1,2)
  end function matrix_to_voigt

  pure function voigt_to_matrix(V) result(M)
    real, intent(in) :: V(6)
    real M(3,3)
    M(1,1) = V(1)
    M(2,2) = V(2)
    M(3,3) = V(3)
    M(2,3) = V(4);    M(3,2) = V(4)
    M(1,3) = V(5);    M(3,1) = V(5)
    M(1,2) = V(6);    M(2,1) = V(6)
  end function voigt_to_matrix

  pure function det3(M) result(det)
    real, intent(in) :: M(3,3)
    real det
    det =   M(1,1) * (M(2,2) * M(3,3) - M(2,3) * M(3,2)) &
         &- M(1,2) * (M(2,1) * M(3,3) - M(2,3) * M(3,1)) &
         &+ M(1,3) * (M(2,1) * M(3,2) - M(2,2) * M(3,1));
  end function det3
  
  pure function inv3(A) result(inv)
    real, intent(in) :: A(3,3)
    real inv(3,3), det, invdet
    
    det = det3(A)
    invdet = 1/det

    inv(1,1) =  (A(2,2)*A(3,3) - A(3,2)*A(2,3)) * invdet;
    inv(1,2) = -(A(1,2)*A(3,3) - A(1,3)*A(3,2)) * invdet;
    inv(1,3) =  (A(1,2)*A(2,3) - A(1,3)*A(2,2)) * invdet;
    inv(2,1) = -(A(2,1)*A(3,3) - A(2,3)*A(3,1)) * invdet;
    inv(2,2) =  (A(1,1)*A(3,3) - A(1,3)*A(3,1)) * invdet;
    inv(2,3) = -(A(1,1)*A(2,3) - A(2,1)*A(1,3)) * invdet;
    inv(3,1) =  (A(2,1)*A(3,2) - A(3,1)*A(2,2)) * invdet;
    inv(3,2) = -(A(1,1)*A(3,2) - A(3,1)*A(1,2)) * invdet;
    inv(3,3) =  (A(1,1)*A(2,2) - A(2,1)*A(1,2)) * invdet;
  end function inv3

  pure function curl_z(u)
    real, intent(in) :: u(:,-1:,-1:)
    real curl_z
    curl_z = 0.125 * ((u(2, 1, 1) + 2*u(2, 1, 0) + u(2, 1,-1)) &
         &          - (u(2,-1, 1) + 2*u(2,-1, 0) + u(2,-1,-1)) &
         &          + (u(1, 1,-1) + 2*u(1, 0,-1) + u(1,-1,-1)) &
         &          - (u(1, 1, 1) + 2*u(1, 0, 1) + u(1,-1, 1)))
  end function curl_z

  pure function div_2d(u)
    real, intent(in) :: u(:,-1:,-1:)
    real div_2d
       div_2d = 0.125 * ((u(1, 1, 1) + 2*u(1, 1, 0) + u(1, 1,-1)) &
            &          - (u(1,-1, 1) + 2*u(1,-1, 0) + u(1,-1,-1)) &
            &          + (u(2, 1, 1) + 2*u(2, 0, 1) + u(2,-1, 1)) &
            &          - (u(2, 1,-1) + 2*u(2, 0,-1) + u(2,-1,-1)))

    !div_2d = 0.5*(u(1,1,0) - u(1,-1,0)) + (u(2,0,1) - u(2,1,-1))
  end function div_2d

  pure function grad_2d(u)
    real, intent(in) :: u(-1:,-1:)
    real grad_2d(2)
    grad_2d(1) = 0.5*(u(1,0) - u(-1,0))
    grad_2d(2) = 0.5*(u(0,1) - u(0,-1))
  end function grad_2d

end module m_matutil

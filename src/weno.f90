module m_weno
  integer, parameter :: M=1
contains
  pure function gM(w,c)
    real gM
    real, intent(in) :: w,c
    gM = w*(c + c*c - 3.0*c*w + w*w)/(c*c + w*(1.0-2.0*c));
  end function gM

  subroutine wenom_reconstruct(w, wL, wR)
    use m_error, only: assert

!   w should have extent [-2,2], but use assumed shape array to avoid a copy
    real, dimension(-2:), intent(in) :: w
    real, intent(out) :: wL, wR
    real, dimension(3), parameter :: dL = [0.1,0.6,0.3]
    real, dimension(3), parameter :: dR = [0.3,0.6,0.1]
    real, dimension(3) :: b, aL, aR, oL, oR
    integer :: j

    real, parameter :: eps = 1E-8

!   call assert(size(w).eq.5,'')

    b(1) = (13.0/12.0)*(    w( 0) - 2.0*w(+1) +     w(+2))**2 &
          + (1.0/ 4.0)*(3.0*w( 0) - 4.0*w(+1) +     w(+2))**2
    b(2) = (13.0/12.0)*(    w(-1) - 2.0*w( 0) +     w(+1))**2 &
          + (1.0/ 4.0)*(    w(-1)             -     w(+1))**2
    b(3) = (13.0/12.0)*(    w(-2) - 2.0*w(-1) +     w( 0))**2 &
          + (1.0/ 4.0)*(    w(-2) - 4.0*w(-1) + 3.0*w( 0))**2

    forall (j=1:3)
       aL(j) = dL(j)/(eps + b(j)**2)
       aR(j) = dR(j)/(eps + b(j)**2)
    end forall

    oL(:) = aL(:) / sum(aL(:))
    oR(:) = aR(:) / sum(aR(:))

    forall (j=1:3)
       aL(j) = gM(oL(j),dL(j))
       aR(j) = gM(oR(j),dR(j))
    end forall

    oL(:) = aL(:) / sum(aL(:))
    oR(:) = aR(:) / sum(aR(:))

    wL = (1.0/6.0) * &
         ( oL(1) * ( 2.0*w(+2) - 7.0*w(+1) + 11.0*w( 0)) &
         + oL(2) * (    -w(+1) + 5.0*w( 0) +  2.0*w(-1)) &
         + oL(3) * ( 2.0*w( 0) + 5.0*w(-1) -      w(-2)))
    wR = (1.0/6.0) * &                          
         ( oR(1) * (    -w(+2) + 5.0*w(+1) +  2.0*w( 0)) &
         + oR(2) * ( 2.0*w(+1) + 5.0*w( 0) -      w(-1)) &
         + oR(3) * (11.0*w( 0) - 7.0*w(-1) +  2.0*w(-2)))
  end subroutine wenom_reconstruct

! at the gaussian quadrature points
  subroutine wenom_reconstruct_gqp(w, wL, wR)
!   w should have extent [-2,2], but use assumed shape array to avoid a copy
    real, dimension(-2:), intent(in) :: w
    real, intent(out) :: wL, wR
    real, dimension(3), parameter :: dL = [(210.0-sqrt(3.0))/1080.0, &
         &                                 11.0/18.0, &
         &                                 (210.0+sqrt(3.0))/1080.0]

    real, dimension(3), parameter :: dR = [(210.0+sqrt(3.0))/1080.0, &
         &                                 11.0/18.0, &
         &                                 (210.0-sqrt(3.0))/1080.0]

    real, dimension(3) :: b, aL, aR, oL, oR
    integer :: j
    real, parameter :: k = sqrt(3.0)/12.0

    real, parameter :: eps = 1E-8

    b(1) = (13.0/12.0)*(    w( 0) - 2.0*w(+1) +     w(+2))**2 &
          + (1.0/ 4.0)*(3.0*w( 0) - 4.0*w(+1) +     w(+2))**2
    b(2) = (13.0/12.0)*(    w(-1) - 2.0*w( 0) +     w(+1))**2 &
          + (1.0/ 4.0)*(    w(-1)             -     w(+1))**2
    b(3) = (13.0/12.0)*(    w(-2) - 2.0*w(-1) +     w( 0))**2 &
          + (1.0/ 4.0)*(    w(-2) - 4.0*w(-1) + 3.0*w( 0))**2

    forall (j=1:3)
       aL(j) = dL(j)/(eps + b(j)**2)
       aR(j) = dR(j)/(eps + b(j)**2)
    end forall

    oL(:) = aL(:) / sum(aL(:))
    oR(:) = aR(:) / sum(aR(:))

    forall (j=1:3)
       aL(j) = gM(oL(j),dL(j))
       aR(j) = gM(oR(j),dR(j))       
    end forall

    oL(:) = aL(:) / sum(aL(:))
    oR(:) = aR(:) / sum(aR(:))

    wL =   oL(1) * ( w(0) + k*( 3.0*w(+2) - 4.0*w(+1) +      w( 0))) &
         + oL(2) * ( w(0) + k*(    -w(+1)             +      w(-1))) &
         + oL(3) * ( w(0) + k*(-3.0*w( 0) + 4.0*w(-1) -      w(-2)))

    wR =   oR(1) * ( w(0) - k*( 3.0*w(+2) - 4.0*w(+1) +      w( 0))) &
         + oR(2) * ( w(0) - k*(    -w(+1)             +      w(-1))) &
         + oR(3) * ( w(0) - k*(-3.0*w( 0) + 4.0*w(-1) -      w(-2)))
  end subroutine wenom_reconstruct_gqp


  function median_state(d,w,wR) result(wRout)
    use m_matutil, only: median, minmod
    real, intent(in) :: d(:,-1:), w(:,-1:), wR(:)
    real, dimension(size(wR)) :: wRout, d_MD, d_LC, wR_UL, wR_MD, wR_LC, wR_min, wR_max
    
    d_MD = minmod(d(:, 0),d(:,1))    
    d_LC = minmod(d(:,-1),d(:,0))

    wR_UL = w(:,0) + 2.0 * (w(:,0) - w(:,-1))
    wR_MD = 0.5 * (w(:,0) + w(:,1) - d_MD)
    wR_LC = w(:,0) + 0.5 * (w(:,0) - w(:,-1)) + (4.0/3.0) * d_LC

    wR_min = max(min(w(:,0),w(:,1),wR_MD), min(w(:,0),wR_UL,wR_LC))
    wR_max = min(max(w(:,0),w(:,1),wR_MD), max(w(:,0),wR_UL,wR_LC))

    wRout = median(wR,wR_min,wR_max)
  end function median_state

end module m_weno

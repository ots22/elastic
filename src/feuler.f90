module m_feuler
  use m_error, only: assert
contains
  subroutine feuler(u,du,L)
    real, dimension(:,:), intent(inout) :: u, du
    interface
       subroutine L(u,du)
         real, intent(in) :: u(:,:)
         real, intent(out) :: du(:,:)
       end subroutine L
    end interface
    
    call L(u,du)
    u = u + du
  end subroutine feuler
end module m_feuler

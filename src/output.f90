module m_output
  private
  public visit_output
contains
!   subroutine make_output(u,it,t,eq,sol,psol)
!     use m_matutil, only: pi, matrix_to_voigt, rot_z, identity, curl_z, div_2d, inv3
!     use m_state
!     use m_eos
!     use m_domain

!     integer, intent(in) :: u, it
!     class(eos), intent(in) :: eq
!     real, intent(in) :: t, sol(:,:,:), psol(:,:,:)

!     real rho, sig(3,3), v(3), F(3,3), Fv(3,3,-1:1,-1:1), &
!          & rhoFv(3,3,-1:1,-1:1), S, gv(3,3,-1:1,-1:1)
!     integer ix,iy,i,j
!     real R(3,3)
!     real xy(3)
!     real curlg(3)
!     real divF(3)
!     real div_rhoF(3)

! !   R = rot_z(-pi/6.0)
!     R = identity(3)

!     write (u,'(A)') "x y rho v1 v2 v3 s11 s22 s33 s23 s13 s12 S F11 &
!             &F12 F13 F21 F22 F23 F31 F32 F33 curlg1 curlg2 curlg3 &
!             &div_rhoF1 div_rhoF2 div_rhoF3 it t"
    

!     do iy=1,ny
!        do ix=1,nx
!           xy = [ix-nx/2.0, iy-ny/2.0, 1.0] * (domain_length/nx)
!           xy = matmul(R,xy)
! !          ps = cons_to_prim(eq, sol(:,ix,iy))
!           v = matmul(R,prim_get_v(psol(:,ix,iy)))
!           F = matmul(R,prim_get_F(psol(:,ix,iy)))
!           S = prim_get_S(psol(:,ix,iy))
!           rho = F_density(eq%rho0, F)
!           sig = eq%stress(S, F)
!           if (ix.gt.1.and.ix.lt.nx.and.iy.gt.1.and.iy.lt.ny) then
!              do i=-1,1; do j=-1,1
! !                psv(:,i,j) = cons_to_prim(eq, sol(:,ix+i,iy+j))
!                 Fv(:,:,i,j) = matmul(R,prim_get_F(psol(:,ix+i,iy+j)))
!                 gv(:,:,i,j) = inv3(Fv(:,:,i,j))
!                 rhoFv(:,:,i,j) = matmul(R,cons_get_rhoF(sol(:,ix+i,iy+j)))
!              enddo; enddo
!              forall (j=1:3) 
!                 curlg(j) = curl_z(gv(j,:,:,:))/dx
!                 div_rhoF(j)  = div_2d(rhoFv(:,j,:,:))/dx
!              end forall
!           else
!              curlg(:) = 0
!              div_rhoF(:) = 0
!           end if

!           write (u,*) xy(1:2), rho, v, matrix_to_voigt(sig), S, F, curlg, div_rhoF, it, t
!        end do
!        write (u,*)
!     end do
!     write (u,*)
!     flush(u)
!   end subroutine make_output

  subroutine visit_output(u,it,t,eq,sol,psol)
    use m_matutil, only: curl_z, div_2d, inv3
    use m_state
    use m_eos
    use m_domain

    integer, intent(in) :: u, it
    class(eos), intent(in) :: eq
    real, intent(in) :: t, sol(:,:,:), psol(:,:,:)

    real rho, F(3,3), Fv(3,3,-1:1,-1:1), &
         & rhoFv(3,3,-1:1,-1:1), gv(3,3,-1:1,-1:1)
    integer ix,iy,i,j
!    real R(3,3)
!    real xy(3)
    real curlg(3)
!    real divF(3)
    real div_rhoF(3)

    write (u,'(A)') '# vtk DataFile Version 3.0'
    write (u,'(A)') 'vtk output'
    write (u,'(A)') 'ASCII'
    write (u,'(A)') 'DATASET RECTILINEAR_GRID'
    write (u,'(A)') 'FIELD FieldData 2'
    write (u,'(A)') 'TIME 1 1 double'
    write (u,'(E16.7E3)') t
    write (u,'(A)') 'CYCLE 1 1 int'
    write (u,'(I3)') it
    write (u,'(A,I7,I7,I3)') 'DIMENSIONS ', nx+1, ny+1, 1

    write (u,'(A,I7,1X,A)') 'X_COORDINATES', nx+1, 'FLOAT'
    write (u,'(I7)') 0
    do ix=1,nx
       write (u,'(I7)') ix
    end do

    write (u,'(A,I7,1X,A)') 'Y_COORDINATES', ny+1, 'FLOAT'
    write (u,'(I7)') 0
    do iy=1,ny
       write (u,'(I7)') iy
    end do

    write (u,'(A,I7,1X,A)') 'Z_COORDINATES', 1, 'FLOAT'
    write (u,'(I7)') 0

    write (u,'(A,I15)') 'CELL_DATA', nx*ny

    call write_scalars_header('rho0detF')
    do iy=1,ny; do ix=1,nx
       F = prim_get_F(psol(:,ix,iy))
       rho = F_density(eq%rho0, F)
       write (u,'(E16.7E3)') rho
    end do; end do

    call write_scalars_header('rho')
    do iy=1,ny; do ix=1,nx
       rho = prim_get_rho(psol(:,ix,iy))
       write (u,'(E16.7E3)') rho
    end do; end do

    write (u,'(A)') 'VECTORS velocity FLOAT'
    write (u,'(3E16.7E3)') psol(prim_v,:,:)

    write (u,'(A)') 'TENSORS stress FLOAT'
    do iy=1,ny; do ix=1,nx
    write (u,'(9E16.7E3)') eq%stress(prim_get_S(psol(:,ix,iy)),prim_get_F(psol(:,ix,iy)))
    end do; end do

    call write_scalars_header('S')
    write (u,'(E16.7E3)') psol(prim_S,:,:)

!   F written as many 'SCALARS' because 'TENSORS' is only for symmetric tensors
    call write_scalars_header('F11')
    write (u,'(E16.7E3)') psol(prim_F(1),:,:)
    call write_scalars_header('F21')
    write (u,'(E16.7E3)') psol(prim_F(2),:,:)
    call write_scalars_header('F31')
    write (u,'(E16.7E3)') psol(prim_F(3),:,:)
    call write_scalars_header('F12')
    write (u,'(E16.7E3)') psol(prim_F(4),:,:)
    call write_scalars_header('F22')
    write (u,'(E16.7E3)') psol(prim_F(5),:,:)
    call write_scalars_header('F32')
    write (u,'(E16.7E3)') psol(prim_F(6),:,:)
    call write_scalars_header('F13')
    write (u,'(E16.7E3)') psol(prim_F(7),:,:)
    call write_scalars_header('F23')
    write (u,'(E16.7E3)') psol(prim_F(8),:,:)
    call write_scalars_header('F33')
    write (u,'(E16.7E3)') psol(prim_F(9),:,:)

    write (u,'(A)') 'VECTORS div_rhoF FLOAT'
    do iy=1,ny; do ix=1,nx
       if (ix.gt.1.and.ix.lt.nx.and.iy.gt.1.and.iy.lt.ny) then
          rhoFv(:,:,:,:) = reshape(sol(CONS_RHOF,ix-1:ix+1,iy-1:iy+1),[3,3,3,3])
          forall (j=1:3) div_rhoF(j) = div_2d(rhoFv(:,j,:,:))/dx
          write (u,'(3E16.7E3)') div_rhoF
       else
          write (u,'(3E16.7E3)') 0.0,0.0,0.0
       end if
    end do; end do

    write (u,'(A)') 'VECTORS curl_g FLOAT'
    do iy=1,ny; do ix=1,nx
       if (ix.gt.1.and.ix.lt.nx.and.iy.gt.1.and.iy.lt.ny) then
          Fv(:,:,:,:) = reshape(psol(PRIM_F,ix-1:ix+1,iy-1:iy+1),[3,3,3,3])
          do i=-1,1; do j=-1,1
             gv(:,:,i,j) = inv3(Fv(:,:,i,j))
          end do; end do
          forall (j=1:3) curlg(j) = curl_z(gv(j,:,:,:))/dx
          write (u,'(3E16.7E3)') curlg
       else
          write (u,'(3E16.7E3)') 0.0,0.0,0.0
       end if
    end do; end do
       
  contains
    subroutine write_scalars_header(name)
      character(*), intent(in) :: name
      write (u,'(A,I3)') 'SCALARS ' // trim(name) // ' FLOAT', 1
      write (u,'(A)') 'LOOKUP_TABLE default'
    end subroutine write_scalars_header
  end subroutine visit_output
end module m_output

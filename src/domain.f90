module m_domain
! length of the domain
  real domain_length ! 40.0
! number of cells in the domain
  integer nx,ny
! number of boundary cells
  integer, parameter :: nb=4
! grid resolution
  real dx
! number of quadrature points for the reconstruction (probably a better place to put this)
  integer, parameter :: nquad=2
  real, parameter :: quad_weight(nquad) = [1.0, 1.0]

  namelist/DOMAIN/domain_length,nx,ny
contains
  subroutine init_domain(u)
    use m_error
    integer, intent(in) :: u
    nx = 0; ny = 0
    domain_length = 0
    read(u,nml=DOMAIN)
    call assert(nx.ge.2*nb.and.ny.ge.2*nb, "domain (nx,ny) must be large enough to contain boundary cells (dx and dy both >8)")
    call assert(domain_length.gt.0, "domain_length must be >0")
    dx = domain_length/nx
  end subroutine init_domain
end module m_domain

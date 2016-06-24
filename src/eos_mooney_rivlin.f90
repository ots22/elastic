module m_eos_mooney_rivlin
  use m_eos

  private

  type, public, extends(eos) :: eos_mooney_rivlin
     real lambda0, mu0, theta0, theta1
   contains
     procedure E
     procedure S
     procedure stress
     procedure init_from_config   
  end type eos_mooney_rivlin

end module m_eos_mooney_rivlin

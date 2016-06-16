module m_eos_simple
  use m_eos

  private

  type, public, extends(eos) :: eos_simple
  contains
    procedure E
    procedure S
    procedure stress
 end type eos_simple
 
contains
  function E(this, S, F)
    class(eos_simple) :: this
    real :: E, S, F(3,3)
    intent(in) :: S, F
    E = 0
  end function E

  function S(this, E, F)
    class(eos_simple) :: this
    real :: S, E, F(3,3)
    intent(in) :: E, F
    S = 0
  end function S

  function stress(this, S, F)
    class(eos_simple) :: this
    real :: stress(3,3), S, F(3,3)
    intent(in) :: S, F
    stress =  0
  end function stress
end module m_eos_simple

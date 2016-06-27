module m_plastic_mises_huber
  use m_plastic
  
  private
  
  type, public, extends(plastic_model) :: plastic_model_mises_huber
     real yield_stress
   contains
     procedure yield_f
     procedure df_dstress
     procedure init_from_config
  end type plastic_model_mises_huber
  
contains
  function yield_f(this, stress)
    use m_matutil, only: dev3
    class(plastic_model_mises_huber) this
    intent(in) stress
    real stress(3,3), yield_f
    yield_f = norm2(dev3(stress)) - sqrt(2.0/3.0)*this%yield_stress
  end function yield_f

  function df_dstress(this, stress)
    use m_matutil, only: dev3
    class(plastic_model_mises_huber) this
    intent(in) stress
    real stress(3,3), df_dstress(3,3)
    df_dstress = dev3(stress)/norm2(dev3(stress))
  end function df_dstress

  subroutine init_from_config(this,u)
    use m_error
    class(plastic_model_mises_huber) this
    integer, intent(in) :: u
    real yield_stress
    real, parameter :: sentinal=huge(1.0)
    namelist/PLASTIC_MISES_HUBER/yield_stress
    yield_stress = sentinal
    read(u,nml=PLASTIC_MISES_HUBER)
    this%yield_stress = yield_stress
  end subroutine init_from_config
end module m_plastic_mises_huber

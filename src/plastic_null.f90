! implements a `null' plastic model: one that never yields
module m_plastic_null
  use m_plastic

  private
  
  type, public, extends(plastic_model) :: plastic_model_null
   contains
     procedure yield_f
     procedure df_dstress
  end type plastic_model_null

contains
  function yield_f(this, stress)
    class(plastic_model_null) this
    intent(in) stress
    real stress(3,3), yield_f
    yield_f = -1.0 * huge(1.0)
  end function yield_f

  function df_dstress(this, stress)
    class(plastic_model_null) this
    intent(in) stress
    real stress(3,3), df_dstress(3,3)
    df_dstress = 0.0
  end function df_dstress

end module m_plastic_null

! implements a `null' plastic model: one that never yields
module m_plastic_null
  use m_plastic

  private
  
  type, public, extends(plastic_model) :: plastic_model_null
   contains
     procedure yield_f
     procedure df_dstress
     procedure df_dhardening
  end type plastic_model_null

contains
  function yield_f(this, stress, hardening)
    class(plastic_model_null) this
    intent(in) stress, hardening
    real stress(3,3), yield_f, hardening
    yield_f = -1.0 * huge(1.0)
  end function yield_f

  function df_dstress(this, stress, hardening)
    class(plastic_model_null) this
    intent(in) stress, hardening
    real stress(3,3), df_dstress(3,3), hardening
    df_dstress = 0.0
  end function df_dstress

  function df_dhardening(this, stress, hardening)
    class(plastic_model_null) this
    intent(in) stress, hardening
    real stress(3,3), df_dhardening, hardening
    df_dhardening = 0.0
  end function df_dhardening
end module m_plastic_null

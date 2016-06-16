module m_configurable
  type configurable
   contains
     procedure init_from_config
  end type configurable
contains
  subroutine init_from_config(this, u)
    use m_error
    integer, intent(in) :: u
    class(configurable) this
    call warn("configuration does nothing by default")

    ! suppress unused argument warning
    if (.false.) call no_op(this,u)
  end subroutine init_from_config
end module m_configurable

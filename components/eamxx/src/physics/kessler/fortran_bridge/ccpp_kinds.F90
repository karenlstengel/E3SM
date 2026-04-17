module ccpp_kinds

  #include "eamxx_config.f"
  #ifdef SCREAM_DOUBLE_PRECISION
  # define kind_phys c_double
  #else
  # define kind_phys c_float
  #endif

  implicit none
  private

  public kind_phys

  end module ccpp_kinds
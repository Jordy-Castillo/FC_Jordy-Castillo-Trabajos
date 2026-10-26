program calculo_fision
  use modulo_nucleo, only: dp, M_masa
  implicit none
  real(dp) :: m_u, m_pd, q

  m_u = M_masa(92, 236)
  m_pd = M_masa(46, 118)
  q = m_u - 2.0_dp * m_pd

  print *, '--- CÁLCULOS OBTENIDOS DE TU MÓDULO ---'
  print '(A, F24.8, A)', ' Masa de U-236 (Z=92, A=236):   ', m_u, ' MeV'
  print '(A, F24.8, A)', ' Masa de Pd-118 (Z=46, A=118):  ', m_pd, ' MeV'
  print '(A, F24.8, A)', ' Energía liberada Q calculada:  ', q, ' MeV'
end program calculo_fision

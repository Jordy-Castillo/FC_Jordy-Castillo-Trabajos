
!!! LLAMAMOS NUESTRO MODULO Y DEFINIMOS EL PROGRAMA !!!
program calculo_pregunta4
  use modulo_nucleo, only: dp, M_masa     
  implicit none
  real(dp) :: m_u, m_pd, q, porcentaje             !masa de U-236 y Pd-118 y energia liberada

  m_u = M_masa(92, 236)                            !masa de U-236
  m_pd = M_masa(46, 118)                           !masa de Pd-118
  q = m_u - 2.0_dp * m_pd                          !Energia Liberada
  porcentaje = (q/200.0_dp)*100                    !Porcentaje de energia liberada (regla de 3)

  print '(A, F24.8, A)', ' Masa de U-236 (Z=92, A=236):   ', m_u, ' MeV'
  print '(A, F24.8, A)', ' Masa de Pd-118 (Z=46, A=118):  ', m_pd, ' MeV'
  print '(A, F24.8, A)', ' Energía liberada Q calculada:  ', q, ' MeV'
  print '(A, F24.9, A)', ' Porcentaje de energia liberada: ', porcentaje, '%'

end program calculo_pregunta4

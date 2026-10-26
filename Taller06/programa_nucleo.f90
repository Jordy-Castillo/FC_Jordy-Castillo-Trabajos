!!!!!!!!!! Programa principal !!!!!!!!!! 

!!! DECLARAMOS EL PROGRAMA Y LAS VARIABLES !!!

program programa_nucleo
  use modulo_nucleo, only: dp, termino_vol, termino_sup, termino_coul, termino_asim, E_B, M_masa, z_estable
  implicit none

  integer :: i, a, z                                       !contadores
  integer :: u_nuc, u_valle, u_max                         !unidades logicas de los archivos a usar
  integer :: z_opt, a_global_max, z_global_max             !resultados de la funcion z_estable
  real(dp) :: EB_opt, EB_global_min                        !minimos de la energia de enlace


!!! DEFINIMOS LOS 8 NUCLIDOS DE REFERENCIA EN LISTAS !!! 
!! El orden dado es Pu, Ge, Hf, Ar, Fe, Ca, 235U, 238U
  integer, dimension(8), parameter :: list_z = [94, 32, 72, 18, 26, 20, 92, 92]
  integer, dimension(8), parameter :: list_a = [244, 70, 179, 38, 56, 40, 235, 238]





!!! GENERACIÓN DE EL ARCHIVO "NUCLEOS.DAT" !!!
open(newunit=u_nuc, file='nucleos.dat', status='replace', action='write')

!! escribimos la fila de nucleos.dat !!
  do i = 1, 8
    z= list_z(i)
    a= list_a(i)
    write(u_nuc, "(2I6, 5F18.8, F18.4)") z, a, &
      termino_vol(a), termino_sup(a), termino_coul(z, a), termino_asim(z, a), E_B(z, a), M_masa(z, a)
  end do
close(u_nuc)


!!! GENERACIÓN DE EL ARCHIVO "VALLE.DAT" !!!
open(newunit=u_valle, file='valle.dat', status='replace', action='write')

!! escribimos la fila de valle.dat !!
EB_global_min = huge(1.0_dp) 
 do a = 4, 260
   call z_estable(a, z_opt, EB_opt)                              !llamamos la subrutina para obtener z_opt y EB_opt
   write(u_valle, "(2I6, F18.8)") a, z_opt, EB_opt               !escribimos la fila de valle.dat

   if (EB_opt < EB_global_min) then                              !si EB_opt es menor que EB_global_min entonces se cumple;
     EB_global_min = EB_opt
     a_global_max = a
     z_global_max = z_opt
   endif
 end do
close (u_valle)

!!! GENERACIÓN DE EL ARCHIVO "MAXIMO.DAT" !!!
open(newunit=u_max, file='maximo.dat', status='replace', action='write')

!! escribimos la fila de maximo.dat !!
  write(u_max, "(2I6, F18.8)") a_global_max, z_global_max, EB_global_min
close(u_max)

end program programa_nucleo

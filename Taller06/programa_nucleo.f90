!!!!!!!!!! Programa principal !!!!!!!!!!

!!! DECLARAMOS EL PROGRAMA Y LAS VARIABLES !!!

program programa_nucleo
  use modulo_nucleo, only: dp, termino_vol, termino_sup, termino_coul, termino_asim, E_B, M_masa, z_estable
  implicit none

  integer :: a                                            !contador
  integer :: u_nuc, u_valle, u_max                         !unidades logicas de los archivos a usar
  integer :: z_opt, a_global_max, z_global_max             !resultados de la funcion z_estable
  real(dp) :: EB_opt, EB_global_min                        !minimos de la energia de enlace


!!! DEFINIMOS LOS 8 NUCLIDOS DE REFERENCIA EN LISTAS !!!
!! El orden dado es Pu, Ge, Hf, Ar, Fe, Ca, 235U, 238U

open(newunit=u_nuc, file='nucleos.dat', status='replace', action='write')

!Llamadas a la subrutina auxiliar para los 8 núclidos de referencia (no se pueden usar arreglos)
!La subrutina la pondremos al final para mayor comodidad

call escribir_nucleo(94, 244) ! Pu-244
call escribir_nucleo(32, 70)  ! Ge-70
call escribir_nucleo(72, 179) ! Hf-179
call escribir_nucleo(18, 38)  ! Ar-38
call escribir_nucleo(26, 56)  ! Fe-56
call escribir_nucleo(20, 40)  ! Ca-40
call escribir_nucleo(92, 235) ! U-235
call escribir_nucleo(92, 238) ! U-238
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

!!CREAMOS OTRA SUBRUTINA PARA ESCRIBIR LOS NUCLEOS
!!Subrutina auxiliar para escribir un solo núcleo en nucleos.dat!!

contains
  subroutine escribir_nucleo(z_val, a_val)
    implicit none
    integer, intent(in) :: z_val, a_val

    write(u_nuc, '(2I6, 5F18.8, F18.4)') &
      z_val, a_val, termino_vol(a_val), termino_sup(a_val), &
      termino_coul(z_val, a_val), termino_asim(z_val, a_val), &
      E_B(z_val, a_val), M_masa(z_val, a_val)
  end subroutine escribir_nucleo

end program programa_nucleo

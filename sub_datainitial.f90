  subroutine data_initial(freq)
      common/data_input/  rk0, wl0, eta0, distmin, cnste
      real    rk0, freq, wl0, eta0
      real    distmin
      real    pi
      complex ci
      complex cnste
      pi = 3.14159265     
      ci = (0.00000000, 1.00000000)
      wl0 =	0.3/freq     !²¨³¤  
      rk0 =	2.0*pi/wl0   !²¨Êý
      eta0 =120.0*pi     !²¨×è¿¹Z
      cnste = -30.0*ci/rk0  !-30j/k
      distmin = wl0*.2
      return
  end 
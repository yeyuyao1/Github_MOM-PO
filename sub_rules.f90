!计算判断两个面是否太近
      subroutine rules2(i1,i2, nearpat, irulef,irules,r1,r2,&
      &	   irulef_near,irules_near, distmin )
       implicit none
       integer irulef,irules,nearpat,irulef_near,irules_near,i1,i2
       real d,r1(3),r2(3), distmin,distmin_use
       nearpat = 0
    !      if(i1.eq.i2) then
     d=0.0
    !      else
     d=sqrt( (r1(1)-r2(1))**2+&
      &		 (r1(2)-r2(2))**2+&
      &		 (r1(3)-r2(3))**2  )
    !      end if
    
    distmin_use=distmin!*0.000001
    
    !	print*,'distmin=',distmin,distmin_use
    !	pause'---------------------------'
       if(d.le.distmin_use) then
    !      if(d.le.distmin) then
     nearpat = 1
     irulef = irulef_near
     irules = irules_near
       end if
       return
  end 
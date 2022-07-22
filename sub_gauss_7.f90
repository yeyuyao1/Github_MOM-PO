  !高斯七点积分
    subroutine getr_7(r,rho,r1,r2,r3,w_near)
      implicit none
      integer i,j
      real r(3,7),rho(3,7),r1(3),r2(3),r3(3),w_near(7)
      real v1,v2,v3,v4,v5,wa,wb,wc,vt1(7),vt2(7),vt3(7),st
      st=sqrt(15.0)
      v1=1.0/3.0
      v2=(6.0-st)/21.0
      v3=(9.0+2.0*st)/21.0
      v4=(6.0+st)/21.0
      v5=(9.0-2.0*st)/21.0
      wa=9.0/80.0
      wb=(155.0-st)/2400.0
      wc=(155.0+st)/2400.0
      vt2(1)=v1
      vt3(1)=v1
      vt2(2)=v2
      vt3(2)=v2
      vt2(3)=v2
      vt3(3)=v3
      vt2(4)=v3
      vt3(4)=v2
      vt2(5)=v4
      vt3(5)=v4
      vt2(6)=v4
      vt3(6)=v5
      vt2(7)=v5
      vt3(7)=v4
      do i=1,7
      vt1(i)=1-vt2(i)-vt3(i)
      enddo
      w_near(1)=wa
      do i=2,4
      w_near(i)=wb
      enddo
      do i=5,7
      w_near(i)=wc
      enddo
      do j=1,7
      do i=1,3
      r(i,j)=vt1(j)*r1(i)+vt2(j)*r2(i)+vt3(j)*r3(i)
      rho(i,j)=r(i,j)-r3(i)
      enddo
      enddo
  end 
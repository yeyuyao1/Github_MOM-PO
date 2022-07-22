  subroutine getr_9(r,rho,r1,r2,r3,w_near)
      implicit none
      integer i,j
      real r(3,9),rho(3,9),r1(3),r2(3),r3(3),w_near(9)
      real v1,v2,v3,v4,vt1(9),vt2(9),vt3(9)
      v1=1.0/18.0
      v2=2.0/9.0
      v3=7.0/18.0
      v4=13.0/18.0
      vt2(1)=v1
      vt3(1)=v2
      vt2(2)=v2
      vt3(2)=v1
      vt2(3)=v1
      vt3(3)=v4
      vt2(4)=v4
      vt3(4)=v1
      vt2(5)=v2
      vt3(5)=v3
      vt2(6)=v3
      vt3(6)=v2
      vt2(7)=v3
      vt3(7)=v3
      vt2(8)=v2
      vt3(8)=v4
      vt2(9)=v4
      vt3(9)=v2
      do i=1,9
      vt1(i)=1-vt2(i)-vt3(i)
      w_near(i)=1.0/18.0
      enddo
      do j=1,9
      do i=1,3
      r(i,j)=vt1(j)*r1(i)+vt2(j)*r2(i)+vt3(j)*r3(i)
      rho(i,j)=r(i,j)-r3(i)
      enddo
      enddo
  end 
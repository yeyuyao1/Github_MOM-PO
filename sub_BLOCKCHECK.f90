subroutine blockcheck(ii,sigma,node,maxnode,maxpatch,maxedge,patpnt,norm,rk,rc,iedge)
  implicit none 
  integer ii,maxnode,maxpatch,maxedge,patpnt(3,maxpatch),iedge(4,maxedge)
  real sigma,node(3,maxnode),norm(3),rk(3),rc(3)

  real PQ(3),b(3,3),ki(3)
  real A1(3,3),B1(3)
  integer i ,ia ,isw !signl是1说明是判断是否被PO遮挡，是2说明判断是否被MOM遮挡
  real dotmul ,eps,det 
  !P 是 rc Q 是 rk,rk应当是面心 ，而不是PO公共边中点  
  call vctadd(rk,rc,ki,-1)
  if(dotmul(ki,norm).ge.0.0) then 
      sigma = 0.0 
  else
      do i = 1,maxpatch
              if(i==iedge(3,ii).or.i==iedge(4,ii)) cycle 
          call getnode(b(1,1),node(1,patpnt(1,i)))
          call getnode(b(1,2),node(1,patpnt(2,i)))
          call getnode(b(1,3),node(1,patpnt(3,i)))

          do ia = 1,3
              A1(ia,1) = b(ia,1) - b(ia,3) 
              A1(ia,2) = b(ia,2) - b(ia,3) 
              A1(ia,3) = ki (ia) 
              B1(ia) = rk(ia) - b(ia,3) 
          end do 
 
          call gauss_real(3,A1,B1,eps,isw) 
          det=A1(1,1)*A1(2,2)*A1(3,3)+A1(1,2)*A1(2,3)*A1(3,1)+A1(1,3)*A1(2,1)*A1(3,2)
          det=det-A1(1,3)*A1(2,2)*A1(3,1)-A1(1,1)*A1(2,3)*A1(3,2)-A1(1,2)*A1(2,1)*A1(3,3)

          if(det==0.0) then 
              sigma = 1.0 
              elseif((0.0<=B1(1).and.B1(1)<=1.0).and.(0.0<=B1(2).and.B1(2)<=1.0).and.(0.0<=1.0&
              &-B1(1)-B1(2).and.1.0-B1(1)-B1(2)<=1.0).and.(0.0<=B1(3).and.B1(3)<=1.0))then 
                  sigma = 0.0
              else 
                  sigma = 1.0 
          end if

              if((sigma.ne.0.0).and.(sigma.ne.1.0)) then 
                  write(*,*) "error" 
              end if 
              if(sigma.eq.0.0) exit 
      end do 
  end if
end subroutine
!平面波遮挡判断 
    subroutine blockcheck_planewave(i1,sigma,ki,PO_node,maxnode_PO,maxpatch_PO,PO_ipatpnt,xyznorm,PO_ctr)
    implicit none
    integer maxnode_PO,maxpatch_PO,maxpatch,maxnode 
    integer PO_ipatpnt(3,maxpatch_PO)
    real:: b(3,3),ki(3)
    real PO_node(3,maxnode_PO),PO_ctr(3,maxpatch_PO)
    real:: xyzctr(3),xyznorm(3),c 
    real::sigma,dotmul,rvm(3),dotm1,dotm2,alpha
    integer i1,i2,j,isw,num,bin1,bin2,bin3,ia 
    real :: beta(3),ri(3),ri_1(3),r_(3),r_1(3)
    real A1(3,3),B1(3)
    real det,eps 

    !******* 
    real X(3)
    integer L 
    real MS(3) 

    xyzctr(1) = PO_ctr(1,i1)
    xyzctr(2) = PO_ctr(2,i1)
    xyzctr(3) = PO_ctr(3,i1) 
   !自遮挡判断 
    eps=0
    num = 0.0
    c = dotmul(ki,xyznorm)
    !sigma = 0.5
    !write(*,*) "sigma is :",sigma
    if(c.gt.-1.0d-7) then 
        sigma = 0.0        !被自遮挡了 
    else
    ! ! !     !互遮挡
        do j = 1,maxpatch_PO   
            if (j.eq.i1) cycle  !判断面元遍历 ,跳过自身

            if (j.ne.i1) then 
                call getnode(b(1,1) , PO_node(1,PO_ipatpnt(1,j)))
                call getnode(b(1,2) , PO_node(1,PO_ipatpnt(2,j)))
                call getnode(b(1,3) , PO_node(1,PO_ipatpnt(3,j)))

            end if
                ! !     !填充系数矩阵 
                    do ia = 1,3
                        A1(ia,1) = b(ia,1) - b(ia,3)
                        A1(ia,2) = b(ia,2) - b(ia,3)
                        A1(ia,3) = ki(ia) 
                        B1(ia) = xyzctr(ia) - b(ia,3)
                    end do

                    call gauss_real(3,A1,B1,eps,isw)
                    !call AGAUS(A1,B1,3,X,L,MS)
                    det=A1(1,1)*A1(2,2)*A1(3,3)+A1(1,2)*A1(2,3)*A1(3,1)+A1(1,3)*A1(2,1)*A1(3,2)
                    det=det-A1(1,3)*A1(2,2)*A1(3,1)-A1(1,1)*A1(2,3)*A1(3,2)-A1(1,2)*A1(2,1)*A1(3,3)
                    
                    if(-1.0d-7<det.and.det<1.0d-7)then 
                        sigma = 1.0

                    elseif(B1(3)<0) then 
                        sigma = 1.0

                    elseif(((0<B1(1).and.B1(1)<1.0).and.(0<B1(2).and.B1(2)<1.0)).and.(0<1.0&
                        &-B1(1)-B1(2).and.1.0-B1(1)-B1(2)<1.0))then
                        sigma = 0.0

                        else
                          sigma=1.0
                    end if
                    if((sigma.ne.0) .and.(sigma.ne.1)) then
                      write(*,*) "error" 
                    end if
                    if (sigma.eq.0.0) exit
                        
          end do

        end if
end subroutine
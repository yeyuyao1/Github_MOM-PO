subroutine makeBk1(ii,ki,pol,ceresult,PO_iedge,PO_xyznorm,maxedge_PO,maxpatch_PO,&
  &PO_node,maxnode_PO,PO_ipatpnt,maxnode,maxpatch,ipatpnt,xyznode,PO_ctr,sigma_PO)
  implicit none 
  integer ii 
  integer maxedge_PO,maxpatch_PO,maxnode_PO,maxpatch,maxnode
  real ki(3),pol(3)
  complex ceresult 
  integer PO_iedge(4,maxedge_PO),PO_ipatpnt(3,maxpatch_PO),ipatpnt(3,maxpatch)
  integer n1fld,n2fld,n3fld,ipfld
  integer igetn3  ,i ,number,number1 
  real PO_xyznorm(3,maxpatch_PO),PO_node(3,maxnode_PO),xyznode(3,maxnode),PO_ctr(3,maxpatch_PO)
  real n(3,2)
  real H(3),nXH(3)
  real rk(3),r1(3),r2(3),r3(3) !公共边中点  
  real t(3),t1(3),t2(3) 
  real dotmul 
  real pi, rk0,wl0,eta0,distmin
  complex cnste
  complex ci
  real sigma_PO(maxpatch_PO)
  common/aim_input/   pi,ci
  common/data_input/  rk0, wl0, eta0, distmin, cnste
 

        !计算n
        call getnode(n(1,1),PO_xyznorm(1,PO_iedge(3,ii)))
        call getnode(n(1,2),PO_xyznorm(1,PO_iedge(4,ii))) 
        !计算  H=ki X ei
        call xmul(ki,pol,H) 
        H = H / (120.0*pi)

        !计算rk(3)
        call getnode(r1,PO_node(1,PO_iedge(1,ii)))
        call getnode(r2,PO_node(1,PO_iedge(2,ii)))
        rk=0.5*(r1+r2)
        !得到H 
        !H = H * exp(-ci*rk0*dotmul(ki,rk))
        !计算t+t-
        n1fld = PO_iedge(1,ii)
        n2fld = PO_iedge(2,ii)
        ipfld = PO_iedge(3,ii) 
        n3fld = igetn3(n1fld,n2fld,PO_ipatpnt(1,ipfld))
        call getnode(r3,PO_node(1,n3fld))
        call get_t(r1,r2,r3,t1,1.0)
        ipfld = PO_iedge(4,ii)
        n3fld = igetn3(n1fld,n2fld,PO_ipatpnt(1,ipfld))
        call getnode(r3,PO_node(1,n3fld))
        call get_t(r1,r2,r3,t2,-1.0)
        t = t1 + t2
        ceresult = 0.0 
        do i=1,2
          call xmul(n(1,i),H,nXH) 
          ceresult = ceresult+sigma_PO(PO_iedge(2+i,ii))*dotmul(t,nXH)*exp(-ci*rk0*dotmul(ki,PO_ctr(1,PO_iedge(2+i,ii))))
          write(34,*) sigma_PO(PO_iedge(2+i,ii))
        end do
        !测试一下
        ceresult = ceresult/2.0
end subroutine
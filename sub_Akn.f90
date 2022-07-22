!被fill A引用 填充AKn 
    subroutine makeAKn(ii,jj,ceresult,maxnode,maxpatch,maxedge,xyznode,ipatpnt,iedge,xyzctr,xyznorm,edge,paera&
  &,maxnode_PO,maxpatch_PO,maxedge_PO,PO_node,PO_ipatpnt,PO_edge,PO_xyznorm,PO_ctr,PO_area,PO_iedge,sigma0,sigma1)
  implicit none
  integer ii,jj,maxnode,maxpatch,maxedge,maxnode_PO,maxpatch_PO,maxedge_PO
  complex ceresult,ceresult1
  real xyznode(3,maxnode),xyzctr(3,maxpatch),xyznorm(3,maxnode),edge(maxedge),paera(maxpatch)
  real PO_node(3,maxnode_PO),PO_edge(maxedge_PO),PO_xyznorm(3,maxpatch_PO),PO_ctr(3,maxpatch_PO),PO_area(maxpatch_PO)
  integer ipatpnt(3,maxpatch),iedge(4,maxedge),PO_ipatpnt(3,maxpatch_PO),PO_iedge(4,maxedge_PO)
  integer sigma0,sigma1,i,j 

  real pi,rk0, wl0, eta0, distmin,rirj,sigma(2) 
  complex ci ,cnste,cgreen,csum 
  real r1(3),r2(3),rc(3),rk(3),normfld(3,2),rsrc(3),rhosrc(3),rr(3)
  integer n1fld,n2fld,ipfld,n3fld
  integer lsrc,lfld,n1src,n2src,n3src,ipsrc 
  integer igetn3
  real rx1(3),rx2(3),r3(3),t1(3),t2(3),t(3)
  real rsrc_g(3,7),rhosrc_g(3,7),wsrc_g(7),dotmul 
  common/aim_input/   pi,ci
  common/data_input/  rk0, wl0, eta0,   distmin,  cnste
 
  
  ceresult = 0.0 
  
  !计算MOM公共边中点 
  call getnode(r1,xyznode(1,iedge(1,jj)))
  call getnode(r2,xyznode(1,iedge(2,jj)))
  rc = 0.5*(r1+r2)
  !计算PO区域公共边中点 
  call getnode(r1,PO_node(1,PO_iedge(1,ii)))
  call getnode(r2,PO_node(1,PO_iedge(2,ii)))
  rk=0.5*(r1+r2)
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

  do lfld = 1,2 
      ipfld = PO_iedge(2+lfld,ii)
      call getnode(normfld(1,lfld) , PO_xyznorm(1,ipfld))     
      call blockcheck(ii,sigma(lfld),PO_node,maxnode_PO,maxpatch_PO,maxedge_PO,PO_ipatpnt,normfld(1,lfld),&
      &PO_ctr(1,ipfld),rc,PO_iedge)
      if(sigma(lfld)==1.0) then 
          call blockcheck(jj,sigma(lfld),xyznode,maxnode,maxpatch,maxedge,ipatpnt,normfld(1,lfld),&
          &PO_ctr(1,ipfld),rc,iedge)
      end if
      !write(*,*) " sigma(fld) is : ",sigma(lfld)
  end do

  do lfld = 1,2 

      ceresult1 = 0.0 
      if(sigma(lfld)==0.0) cycle
      do lsrc = 1,2 

          n1src = iedge(1,jj)
          n2src = iedge(2,jj)
          ipsrc = iedge(2+lsrc,jj)
          n3src = igetn3(n1src,n2src,ipatpnt(1,ipsrc))
          call getnode(r1,xyznode(1,n1src))
          call getnode(r2,xyznode(1,n2src))
          call getnode(r3,xyznode(1,n3src)) 
          call getr_7(rsrc_g,rhosrc_g,r1,r2,r3,wsrc_g) 

          csum = 0.0 

          do i = 1,7
              do j = 1,3 
                  rsrc(j) = rsrc_g(j,i)
                  rhosrc(j) = rhosrc_g(j,i)
              end do

              call vctadd(rk,rsrc,rr,-1)           
              rirj = sqrt(dotmul(rr,rr))
              cgreen = exp(-ci*rk0*rirj)
              cgreen = cgreen/(4.0*pi*(rirj**3))
              cgreen = wsrc_g(i)*(1.0+ci*rk0*rirj)*cgreen 
              call xmul(rr,rhosrc,rx1)
              call xmul(normfld(1,lfld),rx1,rx2)
              csum = csum + dotmul(t,rx2)*cgreen 
          end do 

          ceresult1 = ceresult1 + float(2*lsrc-3)*edge(jj)*csum 
      end do

      ceresult = ceresult + sigma(lfld)*ceresult1

  end do

end subroutine
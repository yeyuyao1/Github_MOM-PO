  !被fill relate Z12 引用 
    subroutine makeZ12(ii,jj,ceh1,maxnode,maxpatch,maxedge,xyznode,ipatpnt,iedge,xyzctr,xyznorm,edge,paera&
    &,maxnode_PO,maxpatch_PO,maxedge_PO,PO_node,PO_ipatpnt,PO_edge,PO_xyznorm,PO_ctr,PO_area,PO_iedge)
    implicit none 
    integer ii,jj   !ii mom  jj po 
    integer maxnode, maxpatch, maxedge   !这是修改过的 
    integer iedge(4,maxedge),ipatpnt(3,maxpatch)
    real xyznode(3,maxnode), xyzctr(3,maxpatch), xyznorm(3,maxpatch)
    real edge(maxedge),paera(maxpatch)   

    integer maxnode_PO,maxpatch_PO,maxedge_PO  
    real PO_node(3,maxnode_PO),PO_xyznorm(3,maxpatch_PO) 
    real PO_area(maxpatch_PO),PO_edge(maxedge_PO),PO_ctr(3,maxpatch_PO) 
    integer PO_iedge(4,maxedge_PO),PO_ipatpnt(3,maxpatch_PO)  

    real rk0,wl0,eta0,distmin
    complex cnste 
    real pi
    complex ci 

    common/aim_input/   pi,ci
    common/data_input/  rk0, wl0, eta0,   distmin,  cnste

    integer i, j,     igetn3,k11
    integer lfld, n1fld,     n2fld,	n3fld
    integer lsrc, n1src,     n2src,	n3src 
    real signfld ,signsrc 
    real rfld(3), rsrc(3),   rhofld(3), rhosrc(3), rr(3)
    real rirj 
    real r1(3),   r2(3),     r3(3)
    real dotmul 
    real rfld_g(3,7),rhofld_g(3,7),wfld_g(7),rsrc_g(3,7),rhosrc_g(3,7),wsrc_g(7)
    complex cefie,csume
    complex cgreen ,ceh1,csume1 
    integer ipfld ,ipsrc 
    n1fld = iedge(1,ii)    
    n2fld = iedge(2,ii) 
    n1src = PO_iedge(1,jj)
    n2src = PO_iedge(2,jj) 
    call getnode(r1,PO_node(1,n1src))
    call getnode(r2,PO_node(1,n2src)) 
    
    ceh1 = (0.0,0.0) 
    do lfld = 1,2
      signfld = float(3-2*lfld)
      ipfld =  iedge(2+lfld,ii)
      n3fld = igetn3(n1fld,n2fld,ipatpnt(1,ipfld))  
      cefie = 0.0 
      do lsrc = 1,2
            signsrc = float(3-2*lsrc)
            ipsrc = PO_iedge(2+lsrc,jj) 
            n3src = igetn3(n1src,n2src,PO_ipatpnt(1,ipsrc))
            call getnode(r3,PO_node(1,n3src)) 
            csume = 0.0
            call getr_7(rfld_g,rhofld_g,xyznode(1,n1fld),xyznode(1,n2fld),xyznode(1,n3fld),wfld_g)
            call getr_7(rsrc_g,rhosrc_g,PO_node(1,n1src),PO_node(1,n2src),PO_node(1,n3src),wsrc_g)

              do j = 1,7
                do k11 = 1,3
                    rhofld(k11) = rhofld_g(k11,j) 
                    rfld(k11) = rfld_g(k11,j)
                end do

                do i = 1,7
                    do k11 = 1,3
                    rhosrc(k11) = rhosrc_g(k11,i)
                    rsrc(k11) = rsrc_g(k11,i)
                    end do

                    csume1 = dotmul(rhofld,rhosrc) - 4.0/(rk0*rk0) 
                    call vctadd(rfld,rsrc,rr,-1)
                    rirj = sqrt(dotmul(rr,rr))
                    cgreen = exp(-ci*rk0*rirj)/rirj 
                    csume = csume+wfld_g(j)*wsrc_g(i)*csume1*cgreen
                end do
              end do
                cefie = cefie + signsrc*csume
              
      end do
      csume = cefie 
      ceh1 = ceh1+signfld*csume 
   end do
    ceh1 = ceh1*edge(ii)*PO_edge(jj)
    return 
end subroutine
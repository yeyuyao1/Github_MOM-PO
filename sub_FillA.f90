subroutine fill_A(Akn,maxnode,maxpatch,maxedge,xyznode,ipatpnt,iedge,xyzctr,xyznorm,edge,paera,&
  &maxnode_PO,maxpatch_PO,maxedge_PO,PO_node,PO_ipatpnt,PO_iedge,PO_xyznorm,PO_edge,PO_ctr)
implicit none
integer maxnode, maxpatch, maxedge
integer iedge(4,maxedge),ipatpnt(3,maxpatch)
real xyznode(3,maxnode), xyzctr(3,maxpatch), xyznorm(3,maxpatch)
real edge(maxedge),paera(maxpatch) 
real    rk0, wl0, eta0
real    distmin
complex cnste,ceresult
complex Akn(maxedge_PO,maxedge)

integer maxnode_PO,maxpatch_PO,maxedge_PO  
real PO_node(3,maxnode_PO),PO_xyznorm(3,maxpatch_PO) 
real PO_area(maxpatch_PO),PO_edge(maxedge_PO),PO_ctr(3,maxpatch_PO) 
integer PO_iedge(4,maxedge_PO),PO_ipatpnt(3,maxpatch_PO)
real    pi
complex ci
integer ii,jj
integer sigma0,sigma1
common/data_input/  rk0, wl0, eta0, distmin, cnste
common/aim_input/   pi,ci
  sigma0 = 0
  sigma1 = 0
do ii = 1,maxedge_PO
  do jj = 1,maxedge  
    call makeAkn(ii,jj,ceresult,maxnode,maxpatch,maxedge,xyznode,ipatpnt,iedge,xyzctr,xyznorm,edge,paera&
    &,maxnode_PO,maxpatch_PO,maxedge_PO,PO_node,PO_ipatpnt,PO_edge,PO_xyznorm,PO_ctr,PO_area,PO_iedge,sigma0,sigma1)
    Akn(ii,jj) = ceresult/2.0
  end do
end do

!write(*,*) "sigma0 is ",sigma0
!write(*,*) "sigma1 is ",sigma1 
return 
end subroutine 
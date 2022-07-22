  subroutine fill_relateMatrix(zzPM,maxnode,maxpatch,maxedge,xyznode,ipatpnt,iedge,xyzctr,xyznorm,edge,paera&
            &,maxnode_PO,maxpatch_PO,maxedge_PO,PO_node,PO_ipatpnt,PO_iedge,PO_xyznorm,PO_edge,PO_ctr,PO_area)
    implicit none 
    integer maxnode, maxpatch, maxedge
    integer iedge(4,maxedge),ipatpnt(3,maxpatch)
    real xyznode(3,maxnode), xyzctr(3,maxpatch), xyznorm(3,maxpatch)
    real edge(maxedge),paera(maxpatch) 
    integer irule_fld, irule_src, irulef_near, irules_near
    integer igetn3 
    real    rk0, wl0, eta0
    real    distmin, rnear2
    complex cnste, ceh1
    complex zzPM(maxedge,maxedge_PO)

    integer maxnode_PO,maxpatch_PO,maxedge_PO  
    real PO_node(3,maxnode_PO),PO_xyznorm(3,maxpatch_PO) 
    real PO_area(maxpatch_PO),PO_edge(maxedge_PO),PO_ctr(3,maxpatch_PO) 
    integer PO_iedge(4,maxedge_PO),PO_ipatpnt(3,maxpatch_PO)

    common/aim_input2/ irule_fld, irule_src, irulef_near, irules_near
    common/data_input/  rk0, wl0, eta0, distmin, cnste
    common/aim_input/   pi,ci
    common/near/rnear2
    real    pi
    complex ci
    integer ii,jj
 
    open(2,file = " Z_Matirx_12.txt")
    do ii = 1,maxedge 
      do jj = 1,maxedge_PO 
        call makeZ12(ii,jj,ceh1,maxnode,maxpatch,maxedge,xyznode,ipatpnt,iedge,xyzctr,xyznorm,edge,paera&
        &,maxnode_PO,maxpatch_PO,maxedge_PO,PO_node,PO_ipatpnt,PO_edge,PO_xyznorm,PO_ctr,PO_area,PO_iedge)
        zzPM(ii,jj) = ceh1*30.0*rk0*ci
        write(2,*) zzPM(ii,jj) 
      end do 
    end do 
    close (2) 
    return 
  end subroutine
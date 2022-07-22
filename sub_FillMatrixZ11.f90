  !¾ØÁ¿·¨×Ô×è¿¹¾ØÕó
    subroutine fill_matrix( zz, maxnode,maxpatch,maxedge,&
      xyznode,ipatpnt,iedge, xyzctr, xyznorm, edge,paera)


   implicit none

   integer maxnode, maxpatch, maxedge
   integer iedge(4,maxedge),ipatpnt(3,maxpatch)
   real xyznode(3,maxnode), xyzctr(3,maxpatch), xyznorm(3,maxpatch)
   real edge(maxedge),paera(maxpatch)

   integer irule_fld, irule_src, irulef_near, irules_near
   real    rk0, wl0, eta0
   real    distmin, rnear2
   complex cnste, ceh1

   complex ZZ(maxedge,maxedge)

   common/aim_input2/ irule_fld, irule_src, irulef_near, irules_near
   common/data_input/  rk0, wl0, eta0, distmin, cnste
   common/aim_input/   pi,ci
   common/near/rnear2

   real    pi
   complex ci
   integer ii,jj
 
   !open(56,file="Z_MATRIX_11.txt")
   if(maxedge.eq.0) then
    zz = 0.0
   elseif(maxedge.ne.0) then
      do ii=1,maxedge
            do jj=ii,maxedge

              call	 efie( ii, jj, ceh1,maxnode,maxpatch,maxedge,&
      
              &   xyznode,ipatpnt,iedge,&
            &   xyzctr, xyznorm, edge,paera)
            zz(ii,jj)=ceh1*(30.0*rk0*ci)

            zz(jj,ii)=zz(ii,jj)
            !write(56,*) zz(ii,jj)
            enddo
        enddo
    end if
       !close(56)

    return
  end
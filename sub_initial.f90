  !三角形面元基本信息计算 
    subroutine initial( maxnode, maxpatch, maxedge,&
      &xyznode, ipatpnt, iedge, xyzctr, xyznorm, edge, paera)
     implicit none
     integer maxnode, maxpatch, maxedge
     integer iedge(4,maxedge),ipatpnt(3,maxpatch)
     real xyznode(3,maxnode),xyzctr(3,maxpatch)
     real edge(maxedge),paera(maxpatch)
     real xyznorm(3,maxpatch)
     real rnear2, edgemax
     common/near/rnear2
     real e1,e2,e3,  dx,dy,dz,  s,tolaera
     real r1(3),   r2(3),    r3(3),   rc(3)
     real dotmul
     real r12(3), r23(3), r31(3), tmp(3)
     integer ia,i,n1,n2,n3
     !	compute the length of edge--edge     
     !	-------------------------------------------------------------------
     edgemax=0.0
     do 10 ia=1,maxedge
     n1=iedge(1,ia)
     n2=iedge(2,ia)
     dx=xyznode(1,n1)-xyznode(1,n2)
     dy=xyznode(2,n1)-xyznode(2,n2)
     dz=xyznode(3,n1)-xyznode(3,n2)
     edge(ia)=sqrt(dx*dx+dy*dy+dz*dz)
     if(edge(ia).gt.edgemax) edgemax=edge(ia)
   10 continue
     rnear2=1.02*edgemax*edgemax
     !	compute the area of triangle--parea
     !	--------------------------------------------------------------------
     do 20 ia=1,maxpatch
     n1=ipatpnt(1,ia)
     n2=ipatpnt(2,ia)
     n3=ipatpnt(3,ia)
     do i=1,3
     r1(i)=xyznode(i,n1)
     r2(i)=xyznode(i,n2)
     r3(i)=xyznode(i,n3)
     end do
     call vctadd(r2,r1,r12,-1)
     call vctadd(r3,r2,r23,-1)
     call vctadd(r1,r3,r31,-1)
     e1=sqrt( dotmul(r12,r12) )
     e2=sqrt( dotmul(r23,r23) )
     e3=sqrt( dotmul(r31,r31) )
     s=(e1+e2+e3)*0.5
     paera(ia)=sqrt( abs(s*(s-e1)*(s-e2)*(s-e3) ) )
     !	compute the center of triangle--xyzctr
     !	and the normal directions of triangle--xyznorm
     !	-----------------------------------------------------------------------------------------------------------------------
     call vctadd(r1,r2,tmp,1)
     call vctadd(r3,tmp,rc,1)
     call xmul(r12,r23,tmp)
     s=sqrt(dotmul(tmp,tmp) )
     do i=1,3
     xyzctr( i,ia)=rc(i)/3.0
     xyznorm(i,ia)=tmp(i)/s
     end do
   20 continue
     close(1)
     tolaera=0.0
     do ia=1,maxpatch
     tolaera=tolaera+paera(ia)
     end do
     write(*,*) tolaera,'tol'
     return
  end 
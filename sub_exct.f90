  !É¢ÉäµÄ¼¤Àø¾ØÕó 
    subroutine exct( maxnode,maxpatch,maxedge,&
      &	 xyznode,ipatpnt, iedge, xyzctr, xyznorm, edge, paera &
      &	 ,thetai,phii,cexcit,ipol)
 
       implicit none
       integer maxnode, maxpatch, maxedge
 
       integer iedge(4,maxedge),ipatpnt(3,maxpatch)
       real xyznode(3,maxnode),xyzctr(3,maxpatch)
       real edge(maxedge),paera(maxpatch)
       real xyznorm(3,maxpatch)
 
       real wt(7),pol(3)
       real rk0, wl0, eta0
       integer  ipol
 
       complex ci
       complex cexcit(maxedge),cexcite
       integer i, n1,n2,n3,i3,ia, signl,m,n
       real  distmin, polrho
       real anorm(3) ,  dotmul
   !      real    ectr(3),  pctr(3)
       real   rfld(3,7) ,rhofld(3,7)
       complex	ct , cnste
       real(kind=8) thetai, phii
       real(kind=8) thesin, thecos, phisin, phicos 	!New
       real rk(3) ,pi						!New
     common/data_input/  rk0, wl0, eta0,  distmin,  cnste
       common/aim_input/    pi,ci
      
       !ipol=1
       thesin=sin(thetai)
       thecos=cos(thetai)
       phisin=sin(phii)
       phicos=cos(phii)
       rk(1)=rk0*thesin*phicos
       rk(2)=rk0*thesin*phisin
       rk(3)=rk0*thecos
 
     if(ipol.eq.1) then
       pol(1)=thecos*phicos
       pol(2)=thecos*phisin
       pol(3)=-thesin
     endif
     if(ipol.eq.2) then
       pol(1)=-phisin
       pol(2)= phicos
       pol(3)= 0.0
     endif
     open(2,file = "V_POMOM.txt")
     do 100 ia=1, maxedge
 
       n1 = iedge( 1, ia )
       n2 = iedge( 2, ia )
 
     cexcit(ia) = (0.0,0.0)
     cexcite    = (0.0,0.0)
 
     do 10 i=1,2
        i3=iedge(2+i,ia)
        n3= ipatpnt(1,i3)+ipatpnt(2,i3)+ipatpnt(3,i3)-n1-n2
 
        signl=float(3-2*i)
 
     call getr_7(rfld,rhofld,&
      &	    xyznode(1,n1),xyznode(1,n2),xyznode(1,n3), wt)
 
      do 20 m=1,7
          ct=cexp(ci*dotmul(rk,rfld(1,m)))
          polrho=dotmul(rhofld(1,m),pol)
          cexcite = cexcite+ct*polrho*wt(m)*signl
 
   20	   continue
   10	   continue
 
      cexcit(ia)=cexcite*edge(ia)
      write(2,*) cexcit(ia)
   100	  continue
    close(2) 
       return
end
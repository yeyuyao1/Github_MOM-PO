!±»Z11ÒýÓÃ
      subroutine efie( ii, jj, ceh1,maxnode,maxpatch,maxedge,&
      &	 xyznode,ipatpnt,iedge, xyzctr, xyznorm, edge,paera)
    
       implicit none
    
       integer maxnode, maxpatch, maxedge
       integer iedge(4,maxedge),ipatpnt(3,maxpatch)
       real xyznode(3,maxnode),xyzctr(3,maxpatch)
       real edge(maxedge),paera(maxpatch)
       real xyznorm(3,maxpatch)
       integer irule_fld, irule_src, irulef_near, irules_near
       real    rk0, wl0, eta0
       real    distmin, rnear2
       complex cnste
    
       common/aim_input2/ irule_fld, irule_src, irulef_near, irules_near
       common/data_input/  rk0, wl0, eta0,   distmin,  cnste
       common/aim_input/   pi,ci
       common/near/rnear2
       real    pi
       complex ci
    
    
       integer ii,jj, i, j,     igetn3,	irulef0, irules0
       integer lfld, n1fld,     n2fld,	n3fld,	 irulef
       integer lsrc, n1src,     n2src,	n3src,	 irules
       integer k11  ! New
    
       real signfld, signsrc
       real rfld(3), rsrc(3),   rhofld(3), rhosrc(3), rr(3)
       real rirj,    aera
       real r1(3),   r2(3),     r3(3)
       real anorm(3),bnorm(3),  dotmul
       real ectr(3),   pctr(3),	 rhot(3)
       real rf,	    rs,       d      !New
    real rfld_near(3,9),	rhofld_near(3,9)		! New
    real rsrc_near(3,7),	rhosrc_near(3,7)		! New
    real rsrc_near1(3,9),	rhosrc_near1(3,9)		! New
    real wfld_near(9),  wsrc_near(7),   wsrc_near1(9)	! New
       real rfld_far(3,7),    rhofld_far(3,7)		      ! New
    real rsrc_far(3,7),    rhosrc_far(3,7)			! New
    real wfld_far(7),  wsrc_far(7)				! New
       complex cefie,   csume
       complex	 csumh
       complex csume1,	  cgreen	       !New
       complex ceh1
       integer ipfld, ipsrc, nearpat
    
    !     parameter (pi	  = 3.14159265)
    !     parameter (ci	  = (0.00000000, 1.00000000))
    
    !	----------------------------------------------------
    
       n1fld = iedge( 1, ii )
       n2fld = iedge( 2, ii )
       n1src = iedge( 1, jj )
       n2src = iedge( 2, jj )
    
       call getnode(r1,xyznode(1,n1src))
       call getnode(r2,xyznode(1,n2src))
    
       irulef0 = irule_fld
       irules0 = irule_src
       d = 0.0
       do i=1,3
    rf = xyznode(i,n1fld)+xyznode(i,n2fld)
    rs = xyznode(i,n1src)+xyznode(i,n2src)
    d = d + abs((rf-rs)*(rf-rs))
       end do
       if( 0.25*d .le. rnear2 ) then
     irulef0 = irulef_near
     irules0 = irules_near
       end if
    ! ===================
    
       call getedgectr(xyznode(1,n1fld),xyznode(1,n2fld),ectr)
    
       ceh1=(0.0,0.0)
   
   do 600 lfld = 1, 2
     signfld = float(3-2*lfld)
     ipfld = iedge(2+lfld,ii)
     n3fld = igetn3( n1fld, n2fld, ipatpnt(1,ipfld) )
     call getnode( anorm, xyznorm(1,ipfld) )
     call getpatctr( xyzctr(1,ipfld), pctr )
     call vctadd( ectr, pctr, rhot, -1 )
    
     cefie = (0.0,0.0)
     do 400 lsrc = 1, 2
        signsrc = float( 3 - 2 * lsrc )
        ipsrc = iedge( 2+lsrc, jj )
        n3src = igetn3( n1src, n2src,ipatpnt(1,ipsrc) )
        call getnode(r3,xyznode(1,n3src))
    
    !		    ! check the distance of the two patches. if this
    !		    ! distance is small, special treatment is
    !		    ! required to calculate the matrix elements.
    
        irulef = irulef0
        irules = irules0
        call rules2(ipfld,ipsrc,nearpat,irulef,irules,&
      &	      xyzctr(1,ipfld),xyzctr(1,ipsrc),&
      &	      irulef_near,    irules_near, distmin )
        call getnode(bnorm,xyznorm(1,ipfld))
        aera = paera(ipsrc)
        csume = (0.0,0.0)
        csumh = (0.0,0.0)
    !***************************************************************************
    ! The following is for the near term that means ii=jj and source-patch
    !  is same as the field-patch cases
    !***************************************************************************
    if(nearpat.eq.1) then
    call getr_9(rfld_near,rhofld_near,&
      &		    xyznode(1,n1fld),xyznode(1,n2fld),xyznode(1,n3fld),&
      &		    wfld_near)
    call getr_7(rsrc_near,rhosrc_near,&
      &		    xyznode(1,n1src),xyznode(1,n2src),xyznode(1,n3src),&
      &		    wsrc_near)
    
    !    the following is for EFIE part 
    do 100 j=1,9
      do k11=1,3
         rhofld(k11)=rhofld_near(k11,j)
         rfld(k11)=rfld_near(k11,j)
      enddo
    do 100 i=1,7
      do k11=1,3
         rhosrc(k11)=rhosrc_near(k11,i)
         rsrc(k11)=rsrc_near(k11,i)
      enddo
    csume1=dotmul(rhofld,rhosrc)-4.0/(rk0*rk0)
    !	csume1=-dotmul(rhofld,rhosrc)*(rk0*rk0) +4.0
    call vctadd( rfld,rsrc,rr, -1 )
    rirj=sqrt(dotmul(rr,rr))
    cgreen=exp(-ci*rk0*rirj)/rirj
    csume=csume+wfld_near(j)*wsrc_near(i)*csume1*cgreen
    100	continue
        cefie = cefie + signsrc * csume
    else
    !***************************************************************************
    !  The following is for the far term that means ii is not equal to jj
    !***************************************************************************
    
    call getr_7(rfld_far,rhofld_far,&
      &		    xyznode(1,n1fld),xyznode(1,n2fld),xyznode(1,n3fld),&
      &		    wfld_far)
    call getr_7(rsrc_far,rhosrc_far,&
      &		    xyznode(1,n1src),xyznode(1,n2src),xyznode(1,n3src),&
      &		    wsrc_far)
    
    do 200 j=1,7
      do k11=1,3
         rhofld(k11)=rhofld_far(k11,j)
         rfld(k11)=rfld_far(k11,j)
      enddo
    do 200 i=1,7
      do k11=1,3
         rhosrc(k11)=rhosrc_far(k11,i)
         rsrc(k11)=rsrc_far(k11,i)
      enddo
    !	csume1=-dotmul(rhofld,rhosrc)*(rk0*rk0) +4.0
     csume1=dotmul(rhofld,rhosrc)-4.0/(rk0*rk0)
    call vctadd( rfld,rsrc,rr, -1 )
    rirj=sqrt(dotmul(rr,rr))
    cgreen=exp(-ci*rk0*rirj)/rirj
    csume=csume+wfld_far(j)*wsrc_far(i)*csume1*cgreen
    200	continue
        cefie = cefie + signsrc * csume
    endif
    400	 continue
    
      csume = cefie    !*cnste
    !	  csume = cefie    *cnste
      ceh1=ceh1+signfld*csume
    600   continue
    ceh1=ceh1*edge(jj)
    ceh1=ceh1*edge(ii)
       return
  end
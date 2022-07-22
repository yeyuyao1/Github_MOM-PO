!计算远场 
    subroutine farfld( maxnode,maxedge,maxpatch,&
      &	 xyznode, ipatpnt,iedge, edge,paera&
      &	 ,thetai,phii, cfld,crhs)

      implicit none
      integer maxnode, maxedge,maxpatch

      integer iedge(4,maxedge), ipatpnt(3,maxpatch)
      integer n1,n2,n3,i,l,i3, signl,m
      real xyznode(3,maxnode), edge(maxedge)

      real wt(7)
      real    rk0, wl0, eta0
      real pi
      complex ci, ct
      complex cfld(2),ctmp, csmu1, csmu2,crhs(maxedge)
  complex cnste
      real  distmin,paera(maxpatch)

      real dotmul,dot_tmp
      real(kind=8):: thetai, phii, thesin, thecos, phisin, phicos     !New
      real rk(3) ,thet(3), phi(3), rfld(3,7), rhofld(3,7),s(3)
  common/data_input/  rk0, wl0, eta0,  distmin,  cnste

      common/aim_input/    pi,ci


      thesin=sin(thetai)
      thecos=cos(thetai)
      phisin=sin(phii)
      phicos=cos(phii)
      !散射方向
      rk(1)=rk0*thesin*phicos
      rk(2)=rk0*thesin*phisin
      rk(3)=rk0*thecos
      s(1) = thesin*phicos
      s(2) = thesin*phisin
      s(3) = thecos
      !theta极化
  thet(1)=thecos*phicos
  thet(2)=thecos*phisin
  thet(3)=-thesin
  !phi极化
  phi(1)=-phisin
  phi(2)= phicos
  phi(3)= 0.0

  cfld(1)=0.0
  cfld(2)=0.0

  do 100 i=1,maxedge

      n1 = iedge( 1, i )
      n2 = iedge( 2, i )

  csmu1=0.0
  csmu2=0.0

  do 10 l=1,2
      i3=iedge(2+l,i)
      n3= ipatpnt(1,i3)+ipatpnt(2,i3)+ipatpnt(3,i3)-n1-n2

      signl=float(3-2*l)

  call getr_7(rfld,rhofld,&
      &	    xyznode(1,n1),xyznode(1,n2),xyznode(1,n3), wt)

    do 20 m=1,7
        ct=cexp(ci*dotmul(rk,rfld(1,m)))
        ctmp = ct*signl*wt(m)
        !以下是聂在平书中计算方法也是一般的论文的计算方法，实际采用的是经过推导化简之后的（考虑是单位向量且有特殊关系）
      ! dot_tmp = dotmul(thet,s)
      ! dot_tmp = dot_tmp*dotmul(s,rhofld(1,m))
      ! csmu1 =csmu1 + (dot_tmp-dotmul(rhofld(1,m),thet))*ctmp
      ! dot_tmp = dotmul(phi,s)
      ! dot_tmp = dot_tmp*dotmul(s,rhofld(1,m))
      ! csmu2 =csmu2 + (dot_tmp-dotmul(rhofld(1,m),phi))*ctmp
        csmu1=csmu1+dotmul(rhofld(1,m),thet)*ctmp
        csmu2=csmu2+dotmul(rhofld(1,m),phi)*ctmp
  20	   continue
  10	continue

  cfld(1)=cfld(1)+edge(i)*crhs(i)*csmu1
  cfld(2)=cfld(2)+edge(i)*crhs(i)*csmu2

  100	continue
  cfld(1)=cfld(1)*ci*rk0*30.0
  cfld(2)=cfld(2)*ci*rk0*30.0

  return
end
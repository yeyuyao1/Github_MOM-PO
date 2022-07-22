subroutine fill_B(Bk1,maxnode_PO,maxpatch_PO,maxedge_PO,PO_node,PO_ipatpnt,PO_iedge,PO_xyznorm,PO_edge,&
  &thetai,phii,maxnode,maxpatch,ipatpnt,xyznode,PO_ctr,ipol)
  implicit none
  integer maxnode_PO,maxpatch_PO,maxedge_PO,maxnode,maxpatch
  real PO_node(3,maxnode_PO),xyznode(3,maxnode)
  integer PO_ipatpnt(3,maxpatch_PO),PO_iedge(4,maxedge_PO),ipatpnt(3,maxpatch)
  real PO_xyznorm(3,maxpatch_PO),PO_edge(maxedge_PO),PO_ctr(3,maxpatch_PO)
  real(kind=8):: thetai,phii
  complex Bk1(maxedge_PO),ceresult
  real sigma_PO(maxpatch_PO)   !储存遮挡系数！！！！
  integer ii 
  integer ipol 
  real ki(3),pol(3) 
  real(kind = 8):: thesin,thecos,phisin,phicos
  real xyznorm(3)
  integer number,number1 
  number = 0 
  number1 = 0 
  sigma_PO = 0.0                                  
  !ipol=1       !theta极化！！！！
  thesin=sin(thetai)
  thecos=cos(thetai)
  phisin=sin(phii)
  phicos=cos(phii)

  ki(1)=-thesin*phicos      !入射方向应当加个-号吧  其他地方的入射方向也应该排除一下
  ki(2)=-thesin*phisin
  ki(3)=-thecos

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

  !开始计算sigma_PO() 
  do ii = 1 ,maxpatch_PO
    xyznorm(1) = PO_xyznorm(1,ii)
    xyznorm(2) = PO_xyznorm(2,ii)
    xyznorm(3) = PO_xyznorm(3,ii)
    !先判断PO区域 
    call blockcheck_planewave(ii,sigma_PO(ii),ki,PO_node,maxnode_PO,maxpatch_PO,PO_ipatpnt,xyznorm,PO_ctr)
      if(sigma_PO(ii).eq.1)  then
          !再矩量法区域判断，无需判断自身故noself
        call blockcheck_planewave_noself(ii,sigma_PO(ii),ki,xyznode,maxnode,maxpatch,ipatpnt,xyznorm,maxpatch_PO,PO_ctr)
      end if

      IF(sigma_PO(II).eq.1) number = number + 1 
      if(sigma_PO(ii).eq.0) number1 = number1 + 1 
  end do
  !sigma_PO() 计算结束
  open(34,file="sigma.txt")
  do  ii=1,maxedge_PO
    call makeBk1(ii,ki,pol,ceresult,PO_iedge,PO_xyznorm,maxedge_PO,maxpatch_PO,&
    &PO_node,maxnode_PO,PO_ipatpnt,maxnode,maxpatch,ipatpnt,xyznode,PO_ctr,sigma_PO)
    Bk1(ii) = ceresult
  end do
  close(34)
  write(*,*) " PO region incident light number is  : " ,number 
  write(*,*)  "PO region incident dark number1 is : ",number1 

end subroutine  
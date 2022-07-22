!****************·øÉäµÄ¼¤Àø¾ØÕó***************     
subroutine exct_radio( maxnode,maxpatch,maxedge,&
    &	xyznode,ipatpnt, iedge, xyzctr, xyznorm, edge, paera &
    &	,cexcit,num_v,sc,sc_num)
    implicit none
    integer maxnode, maxpatch, maxedge,sc_num
    integer iedge(4,maxedge),ipatpnt(3,maxpatch)
    real xyznode(3,maxnode),xyzctr(3,maxpatch)
    real edge(maxedge),paera(maxpatch)
    real xyznorm(3,maxpatch)
    real wt(7),pol(3)
    real rk0, wl0, eta0
    integer  ipol,num_v!À¡Ô´±ßÊý
    integer sc(sc_num),sc_od(sc_num)
    complex ci
    complex cexcit(maxedge),cexcite
    integer i, n1,n2,n3,i3,ia, signl,m,n
    real  distmin, polrho
    real anorm(3) ,  dotmul
    real   rfld(3,7) ,rhofld(3,7),h,hmin
    complex	ct , cnste
    real rk(3) ,pi	 !New
    common/data_input/  rk0, wl0, eta0,  distmin,  cnste
    common/aim_input/    pi,ci

    hmin=10000.0
    open(57,file="V2.txt")
!¶Ôsc ÅÅÐò ÒÔÃâ´íÂ© 
    sc_od = sc 
    do ia = 1,sc_num-1
        do i = ia+1,sc_num
            if(sc_od(ia).gt.sc_od(i)) then 
                m = sc_od(ia)
                sc_od(ia) = sc_od(i)
                sc_od(i) = m 
            end if
        end do 
    end do
    
      num_v = 1 
      do ia=1,maxedge
        cexcit(ia)=(0.0,0.0)
        if(num_V.le.sc_num) then 
          if(ia==sc_od(num_v))then
          cexcit(ia)=edge(ia)
          num_v = num_v + 1 
          end if
        end if
        write(57,*) cexcit(ia)
      end do
      close(57)
    return
end
!Stack reserve :2120000000
!Enable Large Addresses  Support...2GB
!Terminal Server
!可以计算最大矩阵15000*15000！
program PO_MOM
  implicit none
  integer scatteradio    !用来控制程序是算散射还是辐射
  integer i,j,ia,ij
  integer num_v   !馈源边号  
  integer sc_num !源边数量 
  integer,allocatable::sc(:)
  integer maxnode_PO,maxpatch_PO,maxedge_PO   !PO区域的
  integer maxnode,maxpatch,maxedge !矩量法区域的
  integer,allocatable:: iedge(:,:),ipatpnt(:,:),PO_iedge(:,:),PO_ipatpnt(:,:)
  real,allocatable:: xyznode(:,:),xyzctr(:,:), xyznorm(:,:)
  real,allocatable::edge(:),paera(:)
  real  distmin
  complex,allocatable:: crhs(:) ,zz(:,:) ,zzPM(:,:),Akn(:,:),Bk1(:) , Y(:,:),current(:),exct_1(:),tmp(:)!exct_1是算散射参数的电压矩阵
  complex cfld(2)
  real rnear2
  real(kind = 8)::THETAI,phii                          !决定统一以theta,phi来表述
  real freq,rk0,wl0,eta0, eps,	fld1, fld2
  complex::ci=(0.0,1.0)  !虚数单位
  real :: pi = 3.1415926
  complex cnste,Z0,Y0
  integer  ix ,ie,isw,ii,jj
  integer irule_fld, irule_src, irulef_near, irules_near ,ipol
  real,allocatable:: PO_node(:,:),PO_xyznorm(:,:)
  real,allocatable::PO_area(:),PO_edge(:),PO_ctr(:,:) 
  real pol(3),ki(3)   
  real maxfld!记录E的最大振幅，用于计算方向图
  real t_start,t_end
  real gainfld(360),sumgain 
  complex Es(3),cfld1(2),cfld2(2),cfld3(2)  !cfld1(2) 存放MOM区域的E面H面ES    clfd2(2)存放PO区域
  complex delta_Y
  
  
  real::finish,start 
  complex,allocatable:: Z_mat(:,:),current1(:),V_mat(:),current2(:),current3(:)
  common/aim_input2/ irule_fld, irule_src, irulef_near, irules_near
  common/data_input/  rk0, wl0, eta0, distmin, cnste
  common/aim_input/   pi,ci
  common/near/rnear2  
  !!!
  !PO公共区里要的是pi,ci,rk0,记得到时候把包含这三个的公共区替换进去！
  call cpu_time(start) 
  
  write(*,*) "Start producing data !"
  call Dataprocess_MoM()
  !call Dataprocess_PO()
  write(*,*) "Data process ending !"
  
  write(*,*) "please input the parameter ,1 for scattering ,2 for radiation "
  read(*,*) scatteradio 
  !读取保存PO数据
  open(1,file='edge_PO.txt')
  read(1,*) maxedge_PO
  allocate(PO_iedge(4,maxedge_PO),PO_edge(maxedge_PO),current2(maxedge_PO),current3(maxedge_PO),Bk1(maxedge_PO))
  do ia=1,maxedge_PO
   read(1,*) PO_iedge(1,ia),PO_iedge(2,ia),PO_iedge(3,ia),PO_iedge(4,ia)
  end do
  close(1)

  open(1,file="patch_PO.txt")
  read(1,*) maxpatch_PO
  allocate(PO_ipatpnt(3,maxpatch_PO),PO_xyznorm(3,maxpatch_PO),PO_area(maxpatch_PO),PO_ctr(3,maxpatch_PO))
  do ia = 1,maxpatch_PO
      read(1,*) i,PO_ipatpnt(1,ia),PO_ipatpnt(2,ia),PO_ipatpnt(3,ia)
  end do
  close(1)
  
  open(1,file="node_PO.txt")
  read(1,*) maxnode_PO
  allocate(PO_node(3,maxnode_PO))
  do ia = 1,maxnode_PO
      read(1,*) i,PO_node(1,ia),PO_node(2,ia),PO_node(3,ia)
  end do
  close(1)
  !读取保存矩量法数据 
  open(2,file='edge.txt')
  read(2,*) maxedge 
  if(scatteradio.eq.2) then
      read(2,*) sc_num !有几条加源边 
      allocate(sc(sc_num),current(sc_num),Y(sc_num,sc_num)) !加源边号数组，源边电流数组
  end if
  
  allocate(iedge(4,maxedge),edge(maxedge),crhs(maxedge),zz(maxedge,maxedge)&
      &,zzPM(maxedge,maxedge_PO),Akn(maxedge_PO,maxedge),current1(maxedge),V_mat(maxedge),Z_mat(maxedge,maxedge),exct_1(maxedge),tmp(maxedge))
  !read(2,*) sc_num 
  if(scatteradio.eq.2) then
      do ia = 1 ,sc_num   !加源边的边号存入数组 
        read(2,*) sc(ia)    
      end do 
  end if
  
do ia=1,maxedge
 read(2,*) iedge(1,ia),iedge(2,ia),iedge(3,ia),iedge(4,ia)
end do
close(2)

!	write(*,*) '  xyznode file name'
open(2,file='node.txt')
read(2,*) maxnode 
allocate(xyznode(3,maxnode))
do ia=1,maxnode
  read(2,*) ij,xyznode(1,ia),xyznode(2,ia),xyznode(3,ia)
end do
close(2)

!	write(*,*) '  patch file name'
open(2,file='patch.txt')
read(2,*) maxpatch
allocate(ipatpnt(3,maxpatch),xyzctr(3,maxpatch), xyznorm(3,maxpatch),paera(maxpatch))
do ia=1,maxpatch
  read(2,*) ij,ipatpnt(1,ia),ipatpnt(2,ia),ipatpnt(3,ia)
end do
  close(2)
  
  write(*,*) "Plz input frequence (GHZ):"
  read(*,*)  freq 
  if(scatteradio.eq.1) then 
    write(*,*) " Plz input theta and phi(Degree) : "
    read(*,*) thetai , phii 
    write(*,*) "Plz input Polarise angle : 1 for theta ,2 for phi "
    read(*,*) ipol
    if((ipol.ne.1).and.(ipol.ne.2)) stop
  end if 

  irule_fld=0
  irule_src=0
  irulef_near=0
  irules_near=0
  
  thetai=thetai*pi/180.0
    phii=phii*pi/180.0

    ki(1) = -sin(thetai)*cos(phii)
    ki(2) = -sin(thetai)*sin(phii)
    ki(3) = -cos(thetai) 

  call data_initial(freq)

  !PO区域处理   !先算PO区域的是为了防止MOM区域需要的rnear2被改变 
  call initial (maxnode_PO,maxpatch_PO,maxedge_PO,PO_node,PO_ipatpnt,&
  &PO_iedge,PO_ctr,PO_xyznorm,PO_edge,PO_area)

  !矩量法区域处理
  call initial ( maxnode, maxpatch, maxedge,&
  &xyznode, ipatpnt, iedge, xyzctr, xyznorm, edge, paera)

  !激励阵V1
  if (scatteradio.eq.1) then 
      call exct( maxnode,maxpatch,maxedge,&
      &xyznode,ipatpnt, iedge, xyzctr, xyznorm, edge, paera &
      &,thetai,phii,crhs,ipol)
  elseif (scatteradio.eq.2) then                                            
    call exct_radio( maxnode,maxpatch,maxedge,&
    &	xyznode,ipatpnt, iedge, xyzctr, xyznorm, edge, paera &
    &	,crhs,num_v,sc,sc_num)
  end if
      
  !自阻抗矩阵 Z11
  call CPU_TIME(t_start)
  call fill_matrix( zz, maxnode,maxpatch,maxedge,&
  & xyznode,ipatpnt,iedge, xyzctr, xyznorm, edge, paera)
  call CPU_TIME(t_end)
  write(*,*) "Z11 time consume is : ",t_end-t_start 

  !互阻抗矩阵 Z12
  call CPU_TIME(t_start)
  call fill_relateMatrix(zzPM,maxnode,maxpatch,maxedge,xyznode,ipatpnt,iedge,xyzctr,xyznorm,edge,paera,&
  &maxnode_PO,maxpatch_PO,maxedge_PO,PO_node,PO_ipatpnt,PO_iedge,PO_xyznorm,PO_edge,PO_ctr,PO_area)
  call CPU_TIME(t_end)
  write(*,*) "Z12 time consume is : ",t_end-t_start 
  !耦合矩阵 A 
  call CPU_TIME(t_start)
  call fill_A(Akn,maxnode,maxpatch,maxedge,xyznode,ipatpnt,iedge,xyzctr,xyznorm,edge,paera,&
  maxnode_PO,maxpatch_PO,maxedge_PO,PO_node,PO_ipatpnt,PO_iedge,PO_xyznorm,PO_edge,PO_ctr)
  call CPU_TIME(t_end)
  write(*,*) "A time consume is : ",t_end-t_start
  !A的里面遮挡系数判断，未加自遮挡  不知道是怎么想的 
  !修正系数矢量B
  ! call fill_B(Bk1,maxnode_PO,maxpatch_PO,maxedge_PO,PO_node,PO_ipatpnt,PO_iedge,PO_xyznorm,PO_edge,
  ! &thetai,phii)
  if(scatteradio.eq.1) then 
        call CPU_TIME(t_start)
        call fill_B(Bk1,maxnode_PO,maxpatch_PO,maxedge_PO,PO_node,PO_ipatpnt,PO_iedge,PO_xyznorm,&
          &PO_edge,thetai,phii,maxnode,maxpatch,ipatpnt,xyznode,PO_ctr,ipol)
        call CPU_TIME(t_end)
        write(*,*) "B time consume is : ",t_end-t_start
  elseif(scatteradio.eq.2) then 
        Bk1 = 0.0 
  end if
  !计算Z_mat
  Z_mat = zz + matmul(zzPM,Akn)
  V_mat = crhs - matmul(zzPM,Bk1)
  tmp = matmul(zzPM,Bk1)

  call CPU_TIME(t_start)
  call GAUSS(maxedge,Z_mat,V_mat,EPS,ISW)   !
  call CPU_TIME(t_end)
  write(*,*) "Gauss solve time consume is : ",t_end-t_start


  current1 = V_mat !I1    
  current2 = Bk1+matmul(Akn,current1)
  
  !OPEN(66,file = "CURRENT.txt")
  !do i = 1,maxedge 
  !  write(66,*) i,current1(sc(i))*edge(sc(i))
  !end do
  !close(66)


  !计算RCS部分  

  open(1,file = "Esq.txt")
  open(2,file = "Hsq.txt")
  open(3,file = "Efiled.txt")
  phii=00.0 
  phii = phii*pi/180
  maxfld = 0.0 
  do ii = 1,361
    cfld1 = 0.0
    cfld2 = 0.0
    thetai = float(ii-1)*pi/180.0
    call farfld(maxnode,maxedge,maxpatch,xyznode,ipatpnt,iedge,&
    &edge,paera,thetai,phii,cfld1,current1)
    call farfld(maxnode_PO,maxedge_PO,maxpatch_PO,PO_node,PO_ipatpnt,PO_iedge,PO_edge,PO_area&
    &,thetai,phii,cfld2,current2)
    cfld(1) = cfld1(1) + cfld2(1)         !总的E面（已点乘极化得到分量）
    cfld(2) = cfld1(2) + cfld2(2)         !总的H面
    fld1=cabs(cfld(1))*cabs(cfld(1))*4.0*pi
    fld2=cabs(cfld(2))*cabs(cfld(2))*4.0*pi
    ! if(fld1.lt.1.e-5) fld1=1.e-5
    ! if(fld2.lt.1.e-5) fld2=1.e-5
  
    write(1,*) ii-1,10.*(alog10(fld1)) 
    write(2,*) ii-1,10.*(alog10(fld2)) 
    write(3,*) ii-1,cabs(cfld(1))

  end do
  ! deallocate(iedge,ipatpnt,PO_iedge,PO_ipatpnt,xyznode,xyzctr, xyznorm,edge,paera&
  ! &,crhs,zz,zzPM,Akn,Bk1,PO_node,PO_xyznorm,PO_area,PO_edge,PO_ctr,Z_mat,&
  ! &current1,V_mat,current2)
  ! close(1) 
  ! close(2)
  ! close(3) 
  
  !******************************************************************************
  !！！！！计算散射参数！！！！！
  !*******************************************************************************
  !S11   lenth,number_sc,exct_1,maxedge
  if (scatteradio.eq.2) then 
          if(sc_num.gt.1) then 
              do ii = 1,sc_num
                  exct_1 = 0.0                    !&*&*&*&*暂时新添加
                  exct_1(sc(ii)) = edge(sc(ii))     !&*&*&*&*暂时新添加
                    V_mat = exct_1 - tmp            !&*&*&*&*暂时新添加
                    Z_mat = zz + matmul(zzPM,Akn)!&*&*&*&*暂时新添加
                    call GAUSS(maxedge,Z_mat,V_mat,EPS,ISW)     !&*&*&*&*暂时新添加
                    current1 = V_mat !矩量法区域的电流 
                    write(*,*) "源边 电流 为 ： ",current1(sc(ii))*edge(sc(ii))
                do ia = 1,sc_num 
            
                    !call S_calcu(maxedge,sc(ii),exct_1,edge,Z_mat,V_mat,EPS,ISW,tmp)   !只有该口加电压，其他口短路 !&*原有
                    Y(ia,ii) = current1(sc(ia))*edge(sc(ia))/1.0
                    write(*,*) "Y",ia,ii," is : " ,Y(ia,ii)
                end do       
              end do 

    
            !计算散射参数 
              WRITE(*,*) "输入特征阻抗Z0："
              read(*,*) Z0
              Y0 = 1.0/Z0
            delta_Y =(Y(1,1)+Y0)*(Y(2,2)+Y0)-Y(1,2)*Y(2,1)
    
            write(*,*) "S11 is : " ,((Y0-Y(1,1))*(Y0+Y(2,2))+Y(1,2)*Y(2,1))/delta_Y
            write(*,*) "S12 is : " ,-2.0*Y(1,2)*Y0/delta_Y
            write(*,*) "S21 is : " ,-2.0*Y(2,1)*Y0/delta_Y
            write(*,*) "S22 is : " ,((Y0+Y(1,1))*(Y0-Y(2,2))+Y(1,2)*Y(2,1))/delta_Y
          end if
  end if
  
  call cpu_time(finish) 
  write(*,*) " time consuming is : ",finish - start 
  write(*,*) " Compute complete !" 
  read(*,*) 
end program PO_MOM


    


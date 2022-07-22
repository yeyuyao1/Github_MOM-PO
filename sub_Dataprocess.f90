subroutine Dataprocess_MoM()
        implicit none
    character::GRID
    integer::i,maxnode,maxpatch,maxedge,n
    integer pm !决定是mom文件还是PO文件 
    character(len = 20) str 
    integer iedge(4,1000000),ipatpnt(3,3000000)
    integer  maxpatch3,ipat3(3,2000000)
    integer  ii,ij,jj,mm,ia,i0,nn,k 
    real buff,dotmul 
    real,allocatable::node(:,:)
    integer,allocatable::patch(:,:)
    integer:: sp_num = 0 ,a 
    integer,allocatable:: sp_edge(:,:),signl(:)
    real r1(3),r2(3),r3(3) 
    integer n1,n2,n3,p1,p2 
    real norm(3,2)
     
    !write(*,*) "QAQ please input maxnode maxpatch QAQ"
    !write(*,*) "please chose file ,1 for MoM,2 for PO :" 
    pm=1
    !*********************************
   !读取一下文件获取maxnode maxpatch 
    maxnode = 0
    maxpatch = 0
    open(2,file = "feko_MoM.nas")
    do ii = 1,9
        read(2,*) str 
    end do
    do 
          read(2,*,iostat = k ) str
          if(k/=0) exit
          if(str .eq. "GRID*") then
          maxnode = maxnode + 1 
          read(2,*) str,n,buff 
          elseif(str.eq."CTRIA3") then
          maxpatch = maxpatch + 1 
          end if
    end do
      close(2) 
  !***********************************
  !****************正式开始处理数据生成点面边
    open(2,file="feko_MoM.nas")
    do ii = 1,9
        read(2,*) str 
    end do
    maxpatch3=maxpatch*3
    allocate(node(3,maxnode),patch(4,maxpatch))
    do i=1,maxnode
        read(2,*) grid,n,node(1,i),node(2,i),n
        read(2,*) grid,n,node(3,i)
    end do
   ! backspace(2)
    do i=1,maxpatch
        read(2,*) grid,patch(1,i),n,patch(2,i),patch(3,i),patch(4,i)  !序号和点编号
        !write(*,*) i,grid,patch(1,i),n,patch(2,i),patch(3,i),patch(4,i)
    end do
    close(2)
  
    if((pm.ne.1).and.(pm.ne.2)) then 
        write(*,*) " pm error " 
    end if 
  
    if(pm .eq. 1) then 
        open(2,file="node.txt")
        open(3,file="patch.txt")
    elseif(pm.eq.2) then
        open(2,file="node_PO.txt")
        open(3,file="patch_PO.txt")
    end if
    write(2,*) maxnode 
    do i=1,maxnode
        write(2,*) i,node(1,i),node(2,i),node(3,i)
    end do
    close(2)
    write(3,*) maxpatch 
    do i=1,maxpatch
        write(3,*) patch(1,i),patch(2,i),patch(3,i),patch(4,i)
    end do
    close(3)
  
    if(pm.eq.1) then
        open(21,file='patch.txt')
        open(22,file='edge.txt')
    elseif(pm.eq.2) then 
        open(21,file='patch_PO.txt')
        open(22,file='edge_PO.txt')
    end if
         read(21,*) ij 
    do ii=1,maxpatch
        read(21,*) ij,ipatpnt(1,ii),ipatpnt(2,ii),ipatpnt(3,ii)
       
            i0=3*(ii-1)+1
            ipat3(1,i0)=ipatpnt(1,ii)
            ipat3(2,i0)=ipatpnt(2,ii)
            ipat3(3,i0)=ii
            ipat3(1,i0+1)=ipatpnt(2,ii)
            ipat3(2,i0+1)=ipatpnt(3,ii)
            ipat3(3,i0+1)=ii
            ipat3(1,i0+2)=ipatpnt(3,ii)
            ipat3(2,i0+2)=ipatpnt(1,ii)
            ipat3(3,i0+2)=ii
    end do
       
        mm=0
        nn=2*maxpatch
        open(1,file = "sp_num.txt")
    do ii=1,maxpatch3
         do jj=ii+1,maxpatch3
          if (ipat3(1,ii).ne.0) then
           if ((ipat3(1,ii).eq.ipat3(2,jj)).&
              &and.(ipat3(2,ii).eq.ipat3(1,jj))) then
            mm=mm+1
              do ia=1,3
               iedge(ia,mm) =ipat3(ia,ii)
              end do
              iedge(4,mm) =ipat3(3,jj)
            !  ipat3(1,jj)=0
            !  ipat3(2,jj)=0
            maxedge=mm
           end if
           if((ipat3(1,ii).eq.ipat3(1,jj)).&
           &and.(ipat3(2,ii).eq.ipat3(2,jj))) then 
            sp_num = sp_num + 1 
            mm=mm+1
            write(1,*) mm 
              do ia=1,3
               iedge(ia,mm) =ipat3(ia,ii)
              end do
              iedge(4,mm) =ipat3(3,jj)
            !  ipat3(1,jj)=0
            !  ipat3(2,jj)=0
            maxedge=mm
            end if
          end if
         end do
    end do
    close(1)
    write(*,*) " sp_num is : " ,sp_num 
    allocate(sp_edge(3,sp_num))
    open(1,file = "sp_num.txt")

    do ii = 1 ,sp_num
        read(1,*) sp_edge (1,ii) 
    end do
    close(1)

      !  write(22,*) maxedge
    allocate(signl(maxedge))
    signl = 0 
    !     do ii=1,maxedge
    !     write(22,110) iedge(1,ii),iedge(2,ii),iedge(3,ii),iedge(4,ii)
    !     end do
    !    110 format(3x,6I8)

       
       !*&&&&&**********
       !***********开始找特殊边形成最后的边矩阵
       !sp_edge 存放的是同样一条公共边对应的三个编号





       do ii = 1,sp_num
            a = 2 
            do i = 1,maxedge
                if(i == sp_edge(1,ii)) cycle 
                if((iedge(1,i)==iedge(1,sp_edge(1,ii)).and.iedge(2,i)==iedge(2,sp_edge(1,ii)))&
                &.or.(iedge(1,i)==iedge(2,sp_edge(1,ii)).and.iedge(2,i)==iedge(1,sp_edge(1,ii)))) then
                sp_edge(a,ii) = i 
                a = a + 1 
                end if 
                if(a .gt.4) exit 
            end do 
        end do

        !开始做删改操作，操作对象为sp_edge
        !写一个新数组存放删改后的东西，加一个控制符来控制
        !仅适用于【严格垂直】结构，即用于模拟天线的带状垂直于水平面
        signl = 0 
        do i = 1,sp_num
            do ii = 1,3
                n1 = iedge(1,sp_edge(ii,i)) 
                n2 = iedge(2,sp_edge(ii,i)) 
                p1 = iedge(3,sp_edge(ii,i))  !两个面号 
                p2 = iedge(4,sp_edge(ii,i))
                do a = 1,2 
                    n3 = patch(2,iedge(2+a,sp_edge(ii,i)) )+patch(3,iedge(2+a,sp_edge(ii,i)) )+patch(4,iedge(2+a,sp_edge(ii,i)) )
                    n3 = n3 - n1 - n2 
                        call getnode(r1,node(1,n1)) 
                        call getnode(r2,node(1,n2))
                        call getnode(r3,node(1,n3))

                        call getnorm(norm(1,a),r1,r2,r3)
                end do 
                if(dotmul(norm(1,1),norm(1,2)).ne.0.0)  then 
                    signl(sp_edge(ii,i)) = 1     !说明这条公共边是需要删改的公共边  
                end if
                    if (signl(sp_edge(ii,i)) == 1 ) exit 
            end do
                if (signl(sp_edge(ii,i)) == 1 ) then 
                    select case(ii)
                        case(1) 
                            if((iedge(3,sp_edge(2,i)).ne.p1).and.(iedge(3,sp_edge(2,i)).ne.p2)) then 
                                call swap1(iedge(3,sp_edge(2,i)),iedge(4,sp_edge(2,i)))
                            end if
                            if((iedge(3,sp_edge(3,i)).ne.p1).and.(iedge(3,sp_edge(3,i)).ne.p2)) then 
                                call swap1(iedge(3,sp_edge(3,i)),iedge(4,sp_edge(3,i)))
                            end if
                        case(2)
                            if((iedge(3,sp_edge(1,i)).ne.p1).and.(iedge(3,sp_edge(1,i)).ne.p2)) then 
                                call swap1(iedge(3,sp_edge(1,i)),iedge(4,sp_edge(1,i)))
                            end if
                            if((iedge(3,sp_edge(3,i)).ne.p1).and.(iedge(3,sp_edge(3,i)).ne.p2)) then 
                                call swap1(iedge(3,sp_edge(3,i)),iedge(4,sp_edge(3,i)))
                            end if
                        case(3)
                            if((iedge(3,sp_edge(1,i)).ne.p1).and.(iedge(3,sp_edge(1,i)).ne.p2)) then 
                                call swap1(iedge(3,sp_edge(1,i)),iedge(4,sp_edge(1,i)))
                            end if
                            if((iedge(3,sp_edge(2,i)).ne.p1).and.(iedge(3,sp_edge(2,i)).ne.p2)) then 
                                call swap1(iedge(3,sp_edge(2,i)),iedge(4,sp_edge(2,i)))
                            end if
                        case default
                    end select
                end if 
        end do  

        write(22,*) maxedge-sp_num
        do ii=1,maxedge
            if(signl(ii)==1) cycle 
        write(22,110) iedge(1,ii),iedge(2,ii),iedge(3,ii),iedge(4,ii)
        end do
    110 format(3x,6I8)

            
        
        close(21)
        close(22)
  
        deallocate(node,patch,sp_edge,signl)
  write(*,*) " MoM complete!! ",sp_num 
  end 

  

    






subroutine swap1(a,b)
    implicit none 
    integer a,b,c 
    c = a 
    a = b 
    b = c 
end subroutine


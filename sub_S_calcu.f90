subroutine S_calcu(maxedge,num,exct_1,edge,Z_mat,V_mat,EPS,ISW,tmp)
    implicit none 
    
    integer maxedge,num,ISW
    real eps,edge(maxedge)
    complex exct_1(maxedge),Z_mat(maxedge,maxedge),V_mat(maxedge),tmp(maxedge)
    
      CALL V4S_prmt(edge(num),num,exct_1,maxedge)
      V_mat = exct_1 - tmp
      call GAUSS(maxedge,Z_mat,V_mat,EPS,ISW) 
      
end subroutine  
    
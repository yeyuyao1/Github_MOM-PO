!填充二端口计算散射参数时候的激励矩阵（根据情况只在某条边加源其他不加）
    subroutine V4S_prmt(lenth,number_sc,exct_1,maxedge)
    implicit none 
    integer:: maxedge,number_sc
    real lenth 
    complex exct_1(maxedge)
    
    exct_1 = 0.0 
    exct_1(number_sc) = lenth 
    end subroutine  
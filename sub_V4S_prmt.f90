!�����˿ڼ���ɢ�����ʱ��ļ������󣨸������ֻ��ĳ���߼�Դ�������ӣ�
    subroutine V4S_prmt(lenth,number_sc,exct_1,maxedge)
    implicit none 
    integer:: maxedge,number_sc
    real lenth 
    complex exct_1(maxedge)
    
    exct_1 = 0.0 
    exct_1(number_sc) = lenth 
    end subroutine  
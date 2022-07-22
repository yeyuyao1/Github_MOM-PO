!获取垂直于公共边的单位向量  
    subroutine get_t(r1,r2,r3,t,sig)
  implicit none 
  real r1(3),r2(3),r3(3),t(3)
  real r13(3),unit(3)
  real lenth,sig,dotmul 
  call vctadd(r3,r1,r13,-1)
  call vctadd(r2,r1,unit,-1)
  !单位化 
  lenth = sqrt(unit(1)**2+unit(2)**2+unit(3)**2)
  unit(1) = unit(1)/lenth 
  unit(2) = unit(2)/lenth 
  unit(3) = unit(3)/lenth 
  lenth = dotmul(r13,unit)
  unit = unit*lenth 
  call vctadd(unit,r13,t,-1)
  t = t*sig
  lenth  = sqrt(dotmul(t,t)) 
  t = t/lenth 
end subroutine
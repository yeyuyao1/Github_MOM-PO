 !把一个外法向量换给另外一个点 
    subroutine getnorm( an, r1,r2,r3 )
      integer j
      real an(3), r1(3),r2(3),r3(3)
      real r12(3),r13(3)
      real dotmul
      call vctadd(r2,r1,r12,-1)
      call vctadd(r3,r1,r13,-1)
      call xmul(r12,r13,an)
      an = an/sqrt(dotmul(an,an))
      return
  end
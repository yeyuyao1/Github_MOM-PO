  !把面心赋给另外一个点 
    subroutine getpatctr( rc, r )
      integer j
      real r(3),rc(3)
      do j=1,3
    r(j) = rc(j)
      end do
      return
  end
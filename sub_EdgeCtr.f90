  !±ßÖĞµã 
    subroutine getedgectr( r1,r2,ectr )
      integer j
      real ectr(3), r1(3), r2(3)
      do j=1,3
    ectr(j) = 0.5*( r1(j) + r2(j) )
      end do
      return
  end
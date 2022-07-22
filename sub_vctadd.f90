!vector add 
    subroutine vctadd(x,y,xy,job)
      implicit none
      integer job
      real x(3),y(3),xy(3)
      if(job.eq. 1) then
      xy(1)=x(1)+y(1)
      xy(2)=x(2)+y(2)
      xy(3)=x(3)+y(3)
      end if
      if(job.eq.-1) then
      xy(1)=x(1)-y(1)
      xy(2)=x(2)-y(2)
      xy(3)=x(3)-y(3)
      end if
      return
  end 
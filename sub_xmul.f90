  !²æ³Ë
    subroutine xmul(x,y,xy)
      !	  ! cross multiplication of two vectors.
      implicit none
      real x(3),y(3),xy(3)
      xy(1)=x(2)*y(3)-x(3)*y(2)
      xy(2)=x(3)*y(1)-x(1)*y(3)
      xy(3)=x(1)*y(2)-x(2)*y(1)
      return
  end  
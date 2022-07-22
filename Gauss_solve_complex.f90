SUBROUTINE GAUSS(N,A,X,EPS,ISW) 
      COMPLEX A(N,N),X(N),C,W,T
    !	do ii=1,n
    !	write(1,*) ii,real(a(ii,5)),aimag(a(ii,5)),cabs(a(ii,5))
    !	end do
    !	close(1)
     
      NM1=N-1
      DO 10 K=1,NM1
      C=(0.0,0.0)
      DO 2 I=k,N
      IF(CABS(A(I,K)).LE.CABS(C)) GO TO 2
      C=A(I,K)
      I0=I
    2	CONTINUE
      IF(CABS(C).GE.EPS) GO TO 25
      ISW=0
      RETURN
    25	IF(I0.EQ.K) GO TO 50
      DO 60 J=K,N
      W=A(I0,J)
      A(I0,J)=A(K,J)
    60	A(K,J)=W
      T=X(K)
      X(K)=X(I0)
      X(I0)=T
    50	KP1=K+1
      C=1.0/C
      X(K)=X(K)*C
      DO 10 J=KP1,N
      A(K,J)=A(K,J)*C
      DO 20 I=KP1,N
    20	A(I,J)=A(I,J)-A(I,K)*A(K,J)
    10	X(J)=X(J)-A(J,K)*X(K)
      X(N)=X(N)/A(N,N)
      DO 40 I=1,NM1
      I1=N-I
      M=I1+1
      W=(0.0,0.0)
      DO 70 J=M,N
    70	W=W+A(I1,J)*X(J)
      X(I1)=X(I1)-W
    40	ISW=1
      RETURN
END
module precision
 implicit none
 integer, parameter :: dp=kind(1d0)
end module
module io_constants
 integer, parameter :: IW=6
end module
module eigen
 use iso_c_binding, only: c_int64_t
 contains
 function eigen_blas_scope_enter(m) result(n)
 integer,intent(in)::m
 integer(c_int64_t)::n
 n=0
 end function
 subroutine eigen_blas_scope_exit(n)
 integer(c_int64_t),intent(in)::n
 end subroutine
end module
module oqp_linalg
 ! Small dense kernels used by the one-dimensional Davidson regression.
 contains
 subroutine dgemm(ta,tb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
 character,intent(in)::ta,tb
 integer,intent(in)::m,n,k,lda,ldb,ldc
 double precision::alpha,beta,a(*),b(*),c(*)
 integer::i,j,l
 double precision::value
 if(ta/='T'.or.tb/='N')error stop 'unsupported dgemm mode'
 do j=1,n
  do i=1,m
   value=0
   do l=1,k
    value=value+a(l+(i-1)*lda)*b(l+(j-1)*ldb)
   enddo
   c(i+(j-1)*ldc)=alpha*value+beta*c(i+(j-1)*ldc)
  enddo
 enddo
 end subroutine
 subroutine dgemv(t,m,n,alpha,a,lda,x,ix,beta,y,iy)
 character,intent(in)::t
 integer,intent(in)::m,n,lda,ix,iy
 double precision,intent(in)::alpha,beta,a(*),x(*)
 double precision,intent(inout)::y(*)
 integer::i,j
 double precision::value
 if(t/='N')error stop 'unsupported dgemv mode'
 do i=1,m
  value=0
  do j=1,n
   value=value+a(i+(j-1)*lda)*x(1+(j-1)*ix)
  enddo
  y(1+(i-1)*iy)=alpha*value+beta*y(1+(i-1)*iy)
 enddo
 end subroutine
 subroutine dsyev(j,u,n,a,lda,w,work,lwork,info)
 character,intent(in)::j,u
 integer,intent(in)::n,lda,lwork
 integer,intent(out)::info
 double precision,intent(inout)::a(*),work(*)
 double precision,intent(out)::w(*)
 double precision::a11,a12,a22,disc,v1,v2,vn
 if(n==1)then
  w(1)=a(1);a(1)=1;info=0;return
 endif
 if(n/=2.or.j/='V'.or.u/='U')error stop 'unsupported dsyev problem'
 a11=a(1);a12=a(1+lda);a22=a(2+lda)
 disc=sqrt((a11-a22)**2+4*a12*a12)
 w(1)=0.5d0*(a11+a22-disc);w(2)=0.5d0*(a11+a22+disc)
 if(abs(a12)>1d-15)then
  v1=a12;v2=w(1)-a11
 else if(a11<=a22)then
  v1=1;v2=0
 else
  v1=0;v2=1
 endif
 vn=sqrt(v1*v1+v2*v2);v1=v1/vn;v2=v2/vn
 ! Fix a negative augmented-head component so the production sign correction
 ! is exercised deterministically.
 if(v1>0)then
  v1=-v1;v2=-v2
 endif
 a(1)=v1;a(2)=v2;a(1+lda)=-v2;a(2+lda)=v1
 info=0
 end subroutine
end module

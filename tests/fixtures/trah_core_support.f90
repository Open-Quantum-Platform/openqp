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
 if(j/='V'.or.u/='U')error stop 'unsupported dsyev problem'
 if(n>2)then
  call jacobi_eig(n,a,lda,w)
  info=0;return
 endif
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
 ! Cyclic Jacobi for the larger Rayleigh matrices of the augmented-Hessian
 ! Davidson regression; eigenvalues ascending, eigenvectors in the columns.
 subroutine jacobi_eig(n,a,lda,w)
 integer,intent(in)::n,lda
 double precision,intent(inout)::a(lda,*)
 double precision,intent(out)::w(*)
 double precision::s(n,n),v(n,n),off,th,t,c,sn,tmp(n)
 integer::p,q,k,sweep,imin
 s=a(1:n,1:n)
 do q=1,n
  do p=1,q-1
   s(q,p)=s(p,q)
  enddo
 enddo
 v=0
 do p=1,n
  v(p,p)=1
 enddo
 do sweep=1,100
  off=0
  do q=2,n
   do p=1,q-1
    off=off+s(p,q)**2
   enddo
  enddo
  if(off<1d-30*max(1d0,sum([(s(k,k)**2,k=1,n)])))exit
  do q=2,n
   do p=1,q-1
    if(abs(s(p,q))<1d-300)cycle
    th=(s(q,q)-s(p,p))/(2*s(p,q))
    t=sign(1d0,th)/(abs(th)+sqrt(th*th+1))
    c=1/sqrt(t*t+1);sn=t*c
    do k=1,n
     tmp(k)=s(k,p)
     s(k,p)=c*tmp(k)-sn*s(k,q)
     s(k,q)=sn*tmp(k)+c*s(k,q)
    enddo
    do k=1,n
     tmp(k)=s(p,k)
     s(p,k)=c*tmp(k)-sn*s(q,k)
     s(q,k)=sn*tmp(k)+c*s(q,k)
    enddo
    do k=1,n
     tmp(k)=v(k,p)
     v(k,p)=c*tmp(k)-sn*v(k,q)
     v(k,q)=sn*tmp(k)+c*v(k,q)
    enddo
   enddo
  enddo
 enddo
 do k=1,n
  w(k)=s(k,k)
 enddo
 do p=1,n-1
  imin=p
  do q=p+1,n
   if(w(q)<w(imin))imin=q
  enddo
  if(imin/=p)then
   t=w(p);w(p)=w(imin);w(imin)=t
   tmp=v(:,p);v(:,p)=v(:,imin);v(:,imin)=tmp
  endif
 enddo
 a(1:n,1:n)=v
 end subroutine
end module

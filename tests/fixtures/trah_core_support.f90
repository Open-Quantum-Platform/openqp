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
 ! Unused dense-solver entry points deliberately fail: these tests use CG.
 contains
 subroutine dgemm(ta,tb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
 character,intent(in)::ta,tb
 integer,intent(in)::m,n,k,lda,ldb,ldc
 double precision::alpha,beta,a(*),b(*),c(*)
 error stop 'unexpected dense solver'
 end subroutine
 subroutine dgemv(t,m,n,alpha,a,lda,x,ix,beta,y,iy)
 character::t
 integer::m,n,lda,ix,iy
 double precision::alpha,beta,a(*),x(*),y(*)
 error stop 'unexpected dense solver'
 end subroutine
 subroutine dsyev(j,u,n,a,lda,w,work,lwork,info)
 character::j,u
 integer::n,lda,lwork,info
 double precision::a(*),w(*),work(*)
 error stop 'unexpected dense solver'
 end subroutine
end module

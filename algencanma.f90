! *****************************************************************
! *****************************************************************

! This is a main program that calls Algencan to solve a simple
! problem. It is intended to be used as an (incomplete) example of
! usage. Algencan applies to problems of the form
!
! Minimize f(x)
!
! subject to
!
!   heq(x) = 0    heq : R^n \to R^m represents the equality constraints
!   hin(x) <= 0   hin : R^n \to R^p represents the inequality constraints
!   l <= x <= u   l and u \in R^n are the bound constraints.

! *****************************************************************
! *****************************************************************

program algencama

  use bmalgencan, only: algencan
  use iso_c_binding, only: c_ptr, c_loc,c_f_pointer

  implicit none

  ! Re-define this type (pdata_type) anyway you want. Algencan
  ! receives a pointer to a 'structure of this type'. Algencan has no
  ! access to the structure. It simple passes the pointer back to the
  ! user defined subroutines evalf, evalg, evalc, evalj, and
  ! evalhl. So, this is a trade safe way (no common blocks) of passing
  ! to the user-provided routines any information related to the
  ! problem. In this example, it is only being used for the user to
  ! count by itself the number of calls to each routine.
  type :: pdata_type
     integer :: counters(5) = 0
  end type pdata_type

  ! LOCAL SCALARS
  logical :: corrin,extallowed,rhoauto,scale
  integer :: allocerr,hlnnzmax,ierr,inform,istop,jnnzmax,m,maxoutit,n,nwcalls,nwtotit,outiter,p,totiter
  real(kind=8) :: bdsvio,csupn,epsfeas,epscompl,epsopt,f,finish,nlpsupn,rhoini,ssupn,start
  type(pdata_type), target :: pdata
  
  ! LOCAL ARRAYS
  logical, allocatable :: lind(:),uind(:)
  real(kind=8), allocatable :: c(:),lbnd(:),ubnd(:),lambda(:),x(:)

  ! Number of variables

  n = 2
  
  allocate(x(n),lind(n),lbnd(n),uind(n),ubnd(n),stat=allocerr)

  if ( allocerr .ne. 0 ) then
     write(*,*) 'Allocation error.'
     stop
  end if

  ! Initial guess and bound constraints
  
  x(1:n) = 0.0d0

  lind(1:n) = .false.

  uind(1:n) = .false.

  ! Number equality (m) and inequality (p) constraints.
  
  m = 0
  p = 3

  allocate(lambda(m+p),c(m+p),stat=allocerr)

  if ( allocerr .ne. 0 ) then
     write(*,*) 'Allocation error.'
     stop
  end if

  ! Initial guess for the Lagrange multipliers
  
  lambda(1:m+p) = 0.0d0

  ! Number of entries in the Jacobian of the constraints
  
  jnnzmax = n

  ! This should be the number of entries in the Hessian of the
  ! Lagrangian. But, in fact, some extra space is need (to store the
  ! Hessian of the Augmented Lagrangian, whose size is hard to
  ! predict, and/or to store the Jacobian of the KKT system). Thus,
  ! declare it as large as possible.
  
  hlnnzmax = huge( 1 ) / 3

  ! Feasibility, complementarity, and optimality tolerances
  
  epsfeas  = 1.0d-08
  epscompl = 1.0d-08
  epsopt   = 1.0d-08

  ! Maximum number of outer iterations
  
  maxoutit = 50

  ! rhoauto means that Algencan will automatically set the initial
  ! value of the penalty parameter. If you set rhoauto = .false. then
  ! you must set rhoini below with a meaningful value.
  rhoauto = .true.

  if ( .not. rhoauto ) then
     rhoini = 1.0d-08
  end if

  ! scale = .true. means that you allow Algencan to automatically
  ! scale the constraints. In any case, the feasibility tolerance
  ! (epsfeas) will be always satisfied by the UNSCALED original
  ! constraints.
  scale = .false.

  ! extallowed = .true. means that you allow Gencan (the active-set
  ! method used by Algencan to solve the bound-constrained
  ! subproblems) to perform extrapolations. This strategy may use
  ! extra evaluations of the objective function and the constraints
  ! per iterations; but it uses to provide overal savings. You should
  ! test both choices for the problem at hand.
  extallowed = .true.

  ! corrin = .true. means that you allow the inertia of the
  ! Jacobian of the KKT system to be corrected during the acceleration
  ! process. You should test both choices for the problem at hand.
  corrin = .false.
  
  print *, "PARAMETERS VALUES"
  print *, "jnnzmax = ", jnnzmax
  print *, "hlnnzmax = ", hlnnzmax
  print *, "n = ", n
  print *, "x = ", x
  print *, "lind = ", lind
  print *, "lbnd = ", lbnd
  print *, "uind = ", uind
  print *, "ubnd = ", ubnd
  print *, "m = ", m
  print *, "p = ", p
  print *, "lambda = ", lambda
  print *, "epsfeas = ", epsfeas
  print *, "epscompl = ", epscompl
  print *, "epsopt = ", epsopt
  print *, "maxoutit = ", maxoutit
  print *, "scale = ", scale
  print *, "rhoauto = ", rhoauto
  print *, "rhoini = ", rhoini
  print *, "extallowed = ", extallowed
  print *, "corrin = ", corrin
  print *, "inform = ", inform

 call cpu_time(start)
 
 !call algencan(user_evalf,user_evalg,user_evalc,user_evalj,user_evalhl,jnnzmax,hlnnzmax, &
 !    n,x,lind,lbnd,uind,ubnd,m,p,lambda,epsfeas,epscompl,epsopt,maxoutit, &
 !    scale,rhoauto,rhoini,extallowed,corrin,f,csupn,ssupn,nlpsupn,bdsvio, &
 !    outiter,totiter,nwcalls,nwtotit,ierr,istop,c_loc(pdata))
  
  call cpu_time(finish)
  
  write(*,*)
  write(*,*) 'Number of variables                                   = ',n
  write(*,*) 'Number of equality constraints                        = ',m
  write(*,*) 'Number of inequality constraints                      = ',p
  
  write(*,*)
  write(*,*) 'x values = ', x(:)
  write(*,*) '(REPORTED BY SOLVER) istop                            = ',istop
  write(*,*) '(REPORTED BY SOLVER) ierr                             = ',ierr
  write(*,*) '(REPORTED BY SOLVER) f                                = ',f
  write(*,*) '(REPORTED BY SOLVER) csupn                            = ',csupn
  write(*,*) '(REPORTED BY SOLVER) ssupn                            = ',ssupn
  write(*,*) '(REPORTED BY SOLVER) nlpsupn                          = ',nlpsupn
  write(*,*) '(REPORTED BY SOLVER) bounds violation                 = ',bdsvio
  write(*,*) '(REPORTED BY SOLVER) Number of outer iterations       = ',outiter
  write(*,*) '(REPORTED BY SOLVER) Number of inner iterations       = ',totiter
  write(*,*) '(REPORTED BY SOLVER) Number of Newton-KKT trials      = ',nwcalls
  write(*,*) '(REPORTED BY SOLVER) Number of Newton-KKT iterations  = ',nwtotit
  
  write(*,*)
  write(*,*) '(COMPUTED BY CALLER) Number of calls to evalf         = ',pdata%counters(1)
  write(*,*) '(COMPUTED BY CALLER) Number of calls to evalg         = ',pdata%counters(2)
  write(*,*) '(COMPUTED BY CALLER) Number of calls to evalc         = ',pdata%counters(3)
  write(*,*) '(COMPUTED BY CALLER) Number of calls to evalj         = ',pdata%counters(4)
  write(*,*) '(COMPUTED BY CALLER) Number of calls to evalhl        = ',pdata%counters(5)
  write(*,*) '(COMPUTED BY CALLER) CPU time in seconds              = ',finish - start

  
  ! *****************************************************************
  ! *****************************************************************
  ! Just checking ...

  inform = 0
  
  call evalf(n,x,f,inform,c_loc(pdata))
  
  if ( inform .ne. 0 ) then
   write(*,*) 'error when calling evalf in the main file. '
   stop
end if

call evalc(n,x,m,p,c,inform,c_loc(pdata))

  if ( inform .ne. 0 ) then
     write(*,*) 'error when calling evalc in the main file. '
     stop
  end if

  csupn = max( 0.0d0, max( maxval( abs( c(1:m) ) ), maxval( c(m+1:m+p) ) ) )

  bdsvio = max( 0.0d0, max( maxval( lbnd(1:n) - x(1:n), lind(1:n) ), maxval( x(1:n) - ubnd(1:n), uind(1:n) ) ) )

  write(*,*)
  write(*,*) '(COMPUTED BY CALLER) f                                = ',f
  write(*,*) '(COMPUTED BY CALLER) csupn                            = ',csupn
  write(*,*) '(COMPUTED BY CALLER) bounds violation                 = ',bdsvio

  write(*,*)
  write(*,*) 'When a quantity appears as computed by solver and computed by caller, they must coincide.'
  write(*,*) '(In case they do not coincide, please report it as a bug.)'
  ! *****************************************************************
  ! *****************************************************************
  
  deallocate(lind,lbnd,uind,ubnd,x,lambda,c,stat=allocerr)
  
  if ( allocerr .ne. 0 ) then
     write(*,*) 'Deallocation error.'
     stop
  end if
  
  stop

contains
  
  ! *****************************************************************
  ! *****************************************************************

  subroutine evalf(n,x,f,inform,pdataptr)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: n
    integer, intent(inout) :: inform
    real(kind=8), intent(out) :: f
    type(c_ptr), optional, intent(in) :: pdataptr

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)

    ! This routine must compute the objective function.
    
    ! LOCAL SCALARS
    type(pdata_type), pointer :: pdata
    real(kind=8) :: d1, d2, d3
    
    call c_f_pointer(pdataptr,pdata)
    pdata%counters(1) = pdata%counters(1) + 1
   
    d1 = sqrt((-5.0d0 - x(1))**2.0d0 + (10.0d0 - x(2))**2.0d0)
    d2 = sqrt((2.0d0 - x(1))**2.0d0 + (1.0d0 - x(2))**2.0d0)
    d3 = sqrt((10.0d0 - x(1))**2.0d0 + (5.0d0 - x(2))**2.0d0)

    f = d1 + d2 + d3
    
  end subroutine evalf

  ! *****************************************************************
  ! *****************************************************************

  subroutine evalg(n,x,g,inform,pdataptr)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: n
    integer, intent(inout) :: inform
    type(c_ptr), optional, intent(in) :: pdataptr
  
    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)
    real(kind=8), intent(out) :: g(n)
    
    ! This routine must compute the gradient of the objective
    ! function.
    
    ! LOCAL SCALARS
    type(pdata_type), pointer :: pdata
    real(kind=8) :: dfdx1, dfdx2, dfdx3, dfdy1, dfdy2, dfdy3
    
    call c_f_pointer(pdataptr,pdata)
    pdata%counters(2) = pdata%counters(2) + 1

    dfdx1 = (x(1) + 5.0d0)/sqrt(x(1)**2.0d0 + 10.0d0*x(1) + x(2)**2.0d0 + 125.0d0 - 20.0d0*x(2))
    dfdx2 = (x(1) - 2.0d0)/sqrt(x(1)**2.0d0 - 4.0d0*x(1) + x(2)**2.0d0 + 5.0d0 - 2.0d0*x(2))
    dfdx3 = (x(1) - 10.0d0)/sqrt(x(1)**2.0d0 - 20.0d0*x(1) + x(2)**2.0d0 + 125.0d0 - 10.0d0*x(2))

    dfdy1 = (x(2) - 10.0d0)/sqrt(x(2)**2.0d0 - 20.0d0*x(2) + x(1)**2.0d0 + 10.0d0*x(1) + 125.0d0)
    dfdy2 = (x(2) - 1.0d0)/sqrt(x(2)**2.0d0 - 2.0d0*x(2) + x(1)**2.0d0 + 5.0d0 - 4.0d0*x(1))
    dfdy3 = (x(2) - 5.0d0)/sqrt(x(2)**2.0d0 - 10.0d0*x(2) + x(1)**2.0d0 + 125.0d0 - 20.0d0*x(1))

    g(1) = (dfdx1 + dfdx2 + dfdx3)
    g(2) = (dfdy1 + dfdy2 + dfdy3)

  end subroutine evalg

  ! *****************************************************************
  ! *****************************************************************

  subroutine evalc(n,x,m,p,c,inform,pdataptr)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: m,n,p
    integer, intent(inout) :: inform
    type(c_ptr), optional, intent(in) :: pdataptr

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)
    real(kind=8), intent(out) :: c(m+p)
    real(kind=8) :: g1, g2, g3

    ! This routine must compute all the m+p constraints.
    
    ! LOCAL SCALARS
    type(pdata_type), pointer :: pdata
    
    call c_f_pointer(pdataptr,pdata)
    pdata%counters(3) = pdata%counters(3) + 1

    g1 = sqrt((-5.0d0 - x(1))**2.0d0 + (10.0d0 - x(2))**2.0d0) - 10
    g2 = sqrt((2.0d0 - x(1))**2.0d0 + (1.0d0 - x(2))**2.0d0) - 10
    g3 = sqrt((10.0d0 - x(1))**2.0d0 + (5.0d0 - x(2))**2.0d0) - 10
    
    c(1) = g1
    c(2) = g2
    c(3) = g3

  end subroutine evalc

  ! *****************************************************************
  ! *****************************************************************

  subroutine evalj(n,x,m,p,ind,sorted,jsta,jlen,lim,jvar,jval,inform,pdataptr)

    implicit none
    
    ! SCALAR ARGUMENTS
    integer, intent(in) :: lim,m,n,p
    integer, intent(inout) :: inform
    type(c_ptr), optional, intent(in) :: pdataptr

    ! ARRAY ARGUMENTS
    logical, intent(in) :: ind(m+p)
    real(kind=8), intent(in) :: x(n)
    logical, intent(out) :: sorted(m+p)
    integer, intent(out) :: jsta(m+p),jlen(m+p),jvar(lim)
    real(kind=8), intent(out) :: jval(lim)
    
    ! This routine must compute the Jacobian of the constraints. In
    ! fact, only gradients of constraints j such that ind(j) =
    ! .true. need to be computed.
    
    ! LOCAL SCALARS
    integer :: i
    type(pdata_type), pointer :: pdata
    real(kind=8) :: dg1dx, dg1dy, dg2dx, dg2dy, dg3dx, dg3dy
    
    call c_f_pointer(pdataptr,pdata)
    pdata%counters(4) = pdata%counters(4) + 1

    dg1dx = (x(1) + 5.0d0)/sqrt(x(1)**2.0d0 + 10.0d0*x(1) + x(2)**2.0d0 + 125.0d0 - 20.0d0*x(2))
    dg1dy = (x(2) - 10.0d0)/sqrt(x(2)**2.0d0 - 20.0d0*x(2) + x(1)**2.0d0 + 10.0d0*x(1) + 125.0d0)
    
    dg2dx = (x(1) - 2.0d0)/sqrt(x(1)**2.0d0 - 4.0d0*x(1) + x(2)**2.0d0 + 5.0d0 - 2.0d0*x(2))
    dg2dy = (x(2) - 1.0d0)/sqrt(x(2)**2.0d0 - 2.0d0*x(2) + x(1)**2.0d0 + 5.0d0 - 4.0d0*x(1))
    
    dg3dx = (x(1) - 10.0d0)/sqrt(x(1)**2.0d0 - 20.0d0*x(1) + x(2)**2.0d0 + 125.0d0 - 10.0d0*x(2))
    dg3dy = (x(2) - 5.0d0)/sqrt(x(2)**2.0d0 - 10.0d0*x(2) + x(1)**2.0d0 + 125.0d0 - 20.0d0*x(1))

    ! Only gradients of constraints j such that ind(j) = .true. need
    ! to be computed.
    
    if ( ind(1) ) then
       if ( lim .lt. n ) then
          inform = -94
          return
       end if
       
       jval(1) = dg1dx
       jval(2) = dg1dy

       jsta(1) = 1
       jlen(1) = n
       
       jvar(1) = 1
       jvar(2) = 2

       ! Says whether the variables' indices in jvar (related to this
       ! constraint) are in increasing order. In case they are,
       ! Algencan takes advantage of this. Implement sorted gradients
       ! of constraints if you can do this in a natural (cheap)
       ! way. Under no circumnstance use a sorting algorithm. (It is
       ! better to set sorted(1) = .false. in this case.)
       
       sorted(1) = .true.
    end if

    if (ind(2)) then
      if ( lim .lt. n ) then
         inform = -94
         return
      end if

      jval(3) = dg2dx
      jval(4) = dg2dy

      jsta(2) = 3
      jlen(2) = n

      jvar(3) = 1
      jvar(4) = 2

      sorted(2) = .true.
    end if

    if (ind(3)) then
      if ( lim .lt. n ) then
         inform = -94
         return
      end if

      jval(5) = dg3dx
      jval(6) = dg3dy
       
      jsta(3) = 5
      jlen(3) = n

      jvar(5) = 1
      jvar(6) = 2

      sorted(3) = .true.
    end if 
    
  end subroutine evalj

  ! *****************************************************************
  ! *********************       sorted(1) = .true.

  subroutine evalhl(n,x,m,p,lambda,lim,inclf,hlnnz,hlrow,hlcol,hlval,inform,pdataptr)

    implicit none
    
    ! SCALAR ARGUMENTS
    logical, intent(in) :: inclf
    integer, intent(in) :: m,n,lim,p
    integer, intent(out) :: hlnnz
    integer, intent(inout) :: inform
    type(c_ptr), optional, intent(in) :: pdataptr

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: lambda(m+p),x(n)
    integer, intent(out) :: hlrow(lim),hlcol(lim)
    real(kind=8), intent(out) :: hlval(lim)

    ! This routine must compute the Hessian of the Lagrangian. The
    ! Hessian of the objective function must NOT be included if inclf
    ! = .false.
    
    ! LOCAL SCALARS
    type(pdata_type), pointer :: pdata

    ! second order gradient of objective function
    real(kind=8) :: d2fdxx_term1, d2fdxx_term2, d2fdxx_term3, d2fdyy_term1, d2fdyy_term2, d2fdyy_term3
    real(kind=8) :: d2fdxx_num1, d2fdxx_num2, d2fdxx_num3, d2fdyy_num1, d2fdyy_num2, d2fdyy_num3
    real(kind=8) :: d2fdxy_num1, d2fdxy_num2, d2fdxy_num3, d2fdxy_den1, d2fdxy_den2, d2fdxy_den3

    ! second order gradient of g1 restriction
    real(kind=8) :: d2g1dxx_num1, d2g1dxx_term1
    real(kind=8) ::  d2g1dyy_num1, d2g1dyy_term1
    real(kind=8) :: d2g1dxy_num1, d2g1dxy_den1

    ! second order gradient of g2 restriction
    real(kind=8) :: d2g2dxx_num1, d2g2dxx_term1
    real(kind=8) ::  d2g2dyy_num1, d2g2dyy_term1
    real(kind=8) :: d2g2dxy_num1, d2g2dxy_den1

    ! second order gradient of g3 restriction
    real(kind=8) :: d2g3dxx_num1, d2g3dxx_term1
    real(kind=8) ::  d2g3dyy_num1, d2g3dyy_term1
    real(kind=8) :: d2g3dxy_num1, d2g3dxy_den1

    print *, "x_wrap = ", x(:)
    
    call c_f_pointer(pdataptr,pdata)
    pdata%counters(5) = pdata%counters(5) + 1

    hlnnz = 0

    ! If .not. inclf then the Hessian of the objective function must not be included
    
    if ( inclf ) then
       if ( hlnnz + 2 .gt. lim ) then
          inform = -95
          return
       end if
    
       hlnnz = hlnnz + 1
       
       hlrow(hlnnz) = 1
       hlcol(hlnnz) = 1

       d2fdxx_num1 = (x(2) - 20.0d0*x(2) + 100.0d0)
       d2fdxx_num2 = (x(2) - 2.0d0*x(2) + 1.0d0)
       d2fdxx_num3 = (x(2) - 10.0d0*x(2) + 25.0d0)

       d2fdxx_term1 = (x(1)**2.0d0 + 10.0d0*x(1) + x(2)**2.0d0 + 125.0d0 - 20.0d0*x(2))
       d2fdxx_term2 = (x(2)**2.0d0 - 4.0d0*x(1) + x(2)**2.0d0 + 5.0d0 - 2.0d0*x(2))
       d2fdxx_term3 = (x(1)**2.0d0 - 20.0d0*x(1) + x(2)**2.0d0 + 125.0d0 - 10.0d0 * x(2))
       hlval(hlnnz) = d2fdxx_num1/d2fdxx_term1*sqrt(d2fdxx_term1) &
                        + d2fdxx_num2/d2fdxx_term2*sqrt(d2fdxx_term2) &
                        + d2fdxx_num3/d2fdxx_term3*sqrt(d2fdxx_term3)
    
       hlnnz = hlnnz + 1

       hlrow(hlnnz) = 1
       hlcol(hlnnz) = 1

       d2fdxy_num1 = (x(2) - 10.0d0)*(x(1) + 5.0d0)
       d2fdxy_num2 = (x(2) - 1.0d0)*(x(1) - 2.0d0)
       d2fdxy_num3 = (x(2) - 5.0d0)*(x(1) - 10.0d0)

       d2fdxy_den1 = (x(2)**2.0d0 - 20.0d0*x(2) + x(1)**2.0d0 + 10.0d0*x(1) + 125.0d0)**3
       d2fdxy_den2 = (x(2)- 20.0d0*x(2) + x(1)**2.0d0 + 10.0d0 * x(1) + 125.0d0)**3
       d2fdxy_den3 = (x(2)**2.0d0 - 10.0d0*x(2) + x(1)**2.0d0 + 125.0d0 - 20.0d0 * x(1))**3

       hlval(hlnnz) = - 1.0d0 * d2fdxy_num1/sqrt(d2fdxy_den1) &
                     - 1.0d0 * d2fdxy_num2/sqrt(d2fdxy_den2) &
                     - 1.0d0 * d2fdxy_num3/sqrt(d2fdxy_den3)

       hlnnz = hlnnz + 1

       hlrow(hlnnz) = 2
       hlcol(hlnnz) = 2
       
       d2fdyy_num1 = (x(1)**2.0d0 + 10.0d0*x(1) + 25.0d0)
       d2fdyy_num2 = (x(1)**2.0d0 - 4.0d0*x(1) + 4.0d0)
       d2fdyy_num3 = (x(1)**2.0d0 - 20.0d0*x(1) + 100.0d0)

       d2fdyy_term1 = (x(2)**2.0d0 - 20.0d0*x(2) + x(1)**2.0d0 + 10.0d0*x(1) + 125.0d0)
       d2fdyy_term2 = (x(2)**2.0d0 - 2.0d0*x(2) + x(1)**2.0d0 + 5.0d0 - 4.0d0*x(1))
       d2fdyy_term3 = (x(2)**2.0d0 - 10.0d0*x(2) + x(1)**2.0d0 + 125.0d0 - 20.0d0*x(1))
       hlval(hlnnz) = d2fdyy_num1/d2fdyy_term1*sqrt(d2fdyy_term1) &
         + d2fdyy_num2/d2fdyy_term2*sqrt(d2fdyy_term2) &
         + d2fdyy_num3/d2fdyy_term3*sqrt(d2fdyy_term3)
    end if

    ! Note that entries of the Hessian of the Lagrangian can be
    ! repeated. If this is case, them sum of repeated entrances is
    ! considered. This feature simplifies the construction of the
    ! Hessian of the Lagrangian.
    
    if ( hlnnz + 1 .gt. lim ) then
       inform = -95
       return
    end if

   ! second order gradient of g1 restriction
    
   hlnnz = hlnnz + 1
       
   hlrow(hlnnz) = 1
   hlcol(hlnnz) = 1

   d2g1dxx_num1 = (x(2) - 20.0d0*x(2) + 100.0d0)

   d2g1dxx_term1 = (x(1)**2.0d0 + 10.0d0*x(1) + x(2)**2.0d0 + 125.0d0 - 20.0d0*x(2))
   hlval(hlnnz) = lambda(1) * (d2g1dxx_num1/d2g1dxx_term1*sqrt(d2g1dxx_term1))

   hlnnz = hlnnz + 1

   hlrow(hlnnz) = 2
   hlcol(hlnnz) = 1

   d2g1dxy_num1 = (x(2) - 10.0d0)*(x(1) + 5.0d0)

   d2g1dxy_den1 = (x(1)**2.0d0 + 10.0d0*x(1) + x(2)**2.0d0 + 125.0d0 - 20.0d0*x(2))**3

   hlval(hlnnz) = lambda(1) * (- 1.0d0 * d2g1dxy_num1/sqrt(d2g1dxy_den1))

   hlnnz = hlnnz + 1

   hlrow(hlnnz) = 2
   hlcol(hlnnz) = 2
   
   d2g1dyy_num1 = (x(1)**2.0d0 + 10.0d0*x(1) + 25.0d0)

   d2g1dyy_term1 = (x(1)**2.0d0 + 10.0d0*x(1) + x(2)**2.0d0 + 125.0d0 - 20.0d0*x(2))
   hlval(hlnnz) = lambda(1) * (d2g1dyy_num1/d2g1dyy_term1*sqrt(d2g1dyy_term1))

   ! second order gradient of g2 restriction

   hlnnz = hlnnz + 1
       
   hlrow(hlnnz) = 1
   hlcol(hlnnz) = 1

   d2g2dxx_num1 = (x(2) - 2.0d0*x(2) + 1.0d0)

   d2g2dxx_term1 = (x(1)**2.0d0 - 4.0d0*x(1) + x(2)**2.0d0 + 5.0d0 - 2.0d0*x(2))
   hlval(hlnnz) = lambda(2) * (d2g2dxx_num1/d2g2dxx_term1*sqrt(d2g2dxx_term1))

   hlnnz = hlnnz + 1

   hlrow(hlnnz) = 2
   hlcol(hlnnz) = 1

   d2g2dxy_num1 = (x(2) - 1.0d0)*(x(1) - 2.0d0)

   d2g2dxy_den1 = (x(1)**2 - 4.0d0*x(1) + x(2)**2.0d0 + 5.0d0 - 2.0d0*x(2))**3

   hlval(hlnnz) = lambda(2) * (- 1.0d0 * d2g2dxy_num1/sqrt(d2g2dxy_den1))

   hlnnz = hlnnz + 1

   hlrow(hlnnz) = 2
   hlcol(hlnnz) = 2
   
   d2g2dyy_num1 = (x(1)**2.0d0 - 4.0d0*x(1) + 4.0d0)

   d2g2dyy_term1 = (x(1)**2.0d0 + 4.0d0*x(1) + x(2)**2.0d0 + 5.0d0 - 2.0d0*x(2))
   hlval(hlnnz) = lambda(2) * (d2g2dyy_num1/d2g2dyy_term1*sqrt(d2g2dyy_term1))

   ! second order gradient of g3 restriction

   hlnnz = hlnnz + 1
       
   hlrow(hlnnz) = 1
   hlcol(hlnnz) = 1

   d2g3dxx_num1 = (x(2) - 10.0d0*x(2) + 25.0d0)

   d2g3dxx_term1 = (x(1)**2.0d0 - 20.0d0*x(1) + x(2)**2.0d0 + 125.0d0 - 10.0d0*x(2))
   hlval(hlnnz) = lambda(3) * (d2g3dxx_num1/d2g3dxx_term1*sqrt(d2g3dxx_term1))

   hlnnz = hlnnz + 1

   hlrow(hlnnz) = 2
   hlcol(hlnnz) = 1

   d2g3dxy_num1 = (x(2) - 5.0d0)*(x(1) - 10.0d0)

   d2g3dxy_den1 = (x(1)**2.0d0 - 20.0d0*x(1) + x(2)**2.0d0 + 125.0d0 - 10.0d0*x(2))**3

   hlval(hlnnz) = lambda(3) * (- 1.0d0 * d2g3dxy_num1/sqrt(d2g3dxy_den1))

   hlnnz = hlnnz + 1

   hlrow(hlnnz) = 2
   hlcol(hlnnz) = 2
   
   d2g3dyy_num1 = (x(1)**2.0d0 - 20.0d0*x(1) + 100.0d0)

   d2g3dyy_term1 = (x(1)**2.0d0 - 20.0d0*x(1) + x(2)**2.0d0 + 125.0d0 - 10.0d0*x(2))
   hlval(hlnnz) = lambda(3) * (d2g3dyy_num1/d2g3dyy_term1*sqrt(d2g3dyy_term1))
    
  end subroutine evalhl

end program algencama

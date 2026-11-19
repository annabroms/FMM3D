! Check st3ddirectstokg_opt against st3ddirectstokg for several densities,
! output flags, thresholds and one/eight-thread and nested execution.
! Built and run by 'make test' unless OMP=OFF. For bounds checking, run
! from the repository root:
! gfortran -O2 -fopenmp -fcheck=all src/Stokes/stokkernels.f \
!   test/Stokes/test_stokkernels_opt.f90 \
!   -o test/Stokes/int2-test-stokkernels-opt
! test/Stokes/int2-test-stokkernels-opt
program test_stokkernels_opt
  use omp_lib
  implicit none
  integer, parameter :: counts(6)=[1,2,3,4,8,16]
  integer :: k
  character(len=16) :: mode
  call get_command_argument(1,mode)
  call omp_set_dynamic(.false.)
  do k=1,size(counts)
    call check_case(int(counts(k),8),mode/='--serial')
  enddo
  print *, 'PASS: opt nd=1,2,3,4,8,16; flags, prefilled outputs, thresholds'
  print *, 'PASS: empty sources, excluded pairs, zero force and 1/8 threads'
  if (mode/='--serial') print *, 'PASS: nested calls'
contains
  subroutine check_case(nd,nested)
    integer(kind=8), intent(in) :: nd
    logical, intent(in) :: nested
    integer(kind=8), parameter :: ns=97,nt=83
    real(kind=8) :: source(3,ns),target(3,nt),force(nd,3,ns)
    real(kind=8) :: u(nd,3,nt),p(nd,nt),g(nd,3,3,nt)
    real(kind=8) :: v(nd,3,nt),q(nd,nt),h(nd,3,3,nt)
    real(kind=8) :: baseline(nd,3,nt), threshold,zero_force(nd,3,ns)
    integer(kind=8) :: i,j,k,flag
    integer :: threads,batch,failures,cutoff
    do j=1,ns
      source(:,j)=[sin(real(j,8)),cos(real(2*j,8)),sin(real(3*j,8))]
      do i=1,3
        do k=1,nd
          force(k,i,j)=cos(real(j+2*i+3*k,8))
        enddo
      enddo
    enddo
    do j=1,nt
      target(:,j)=[cos(real(j,8)),sin(real(2*j,8)),cos(real(3*j,8))]
    enddo
    target(:,1)=source(:,1)
    target(:,2)=source(:,2)+[1d-5,0d0,0d0]
    source(1,3)=0d0
    target(:,3)=source(:,3)
    target(1,3)=1d0/1024d0
    do cutoff=1,3
      threshold=1d-16
      if (cutoff==2) threshold=1d-3
      if (cutoff==3) threshold=1d0/1024d0
      do threads=1,8,7
        call omp_set_num_threads(threads)
        do flag=0,3
          u=0.25d0;p=0.75d0;g=-0.5d0
          v=u;q=p;h=g
          call st3ddirectstokg(nd,source,force,ns,target,nt,flag,u,p,g,threshold)
          call st3ddirectstokg_opt(nd,source,force,ns,target,nt,flag,v,q,h,threshold)
          if (flag/=1) then
            if (any(u/=v).or.any(p/=q).or.any(g/=h)) error stop 'Fallback differs'
          else
            if (maxval(abs(u-v))>3d-12*max(1d0,maxval(abs(u)))) error stop 'Velocity differs'
            if (any(q/=0.75d0).or.any(h/=-0.5d0)) error stop 'Disabled fields changed'
          endif
        enddo
      enddo
    enddo
    ! Empty sources and all-excluded inputs must leave every field intact.
    u=0.25d0;p=0.75d0;g=-0.5d0
    call st3ddirectstokg_opt(nd,source,force,0_8,target,nt,1_8,u,p,g,1d-16)
    if (any(u/=0.25d0).or.any(p/=0.75d0).or.any(g/=-0.5d0)) error stop 'Empty inputs'
    call st3ddirectstokg_opt(nd,source,force,ns,target,nt,1_8,u,p,g,100d0)
    if (any(u/=0.25d0).or.any(p/=0.75d0).or.any(g/=-0.5d0)) error stop 'Excluded pairs'
    zero_force=0d0
    call st3ddirectstokg_opt(nd,source,zero_force,ns,target,nt,1_8,u,p,g,1d-16)
    if (any(u/=0.25d0).or.any(p/=0.75d0).or.any(g/=-0.5d0)) error stop 'Zero force'
    if (nested) then
      call omp_set_num_threads(1)
      baseline=0.25d0;p=0.75d0;g=-0.5d0
      call st3ddirectstokg_opt(nd,source,force,ns,target,nt,1_8,baseline,p,g,1d-16)
      call omp_set_num_threads(8)
      call omp_set_max_active_levels(2)
      failures=0
      !$omp parallel do private(u,p,g) reduction(+:failures)
      do batch=1,8
        u=0.25d0;p=0.75d0;g=-0.5d0
        call st3ddirectstokg_opt(nd,source,force,ns,target,nt,1_8,u,p,g,1d-16)
        if (any(u/=baseline)) failures=failures+1
      enddo
      !$omp end parallel do
      if (failures/=0) error stop 'Nested outputs differ'
    endif
  end subroutine check_case
end program test_stokkernels_opt

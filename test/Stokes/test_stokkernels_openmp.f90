! Check direct Stokes target parallelism, output flags and nested calls.
! Built and run by 'make test' unless OMP=OFF. For bounds checking, run
! from the repository root:
! gfortran -O2 -fopenmp -fcheck=all src/Stokes/stokkernels.f \
!   test/Stokes/test_stokkernels_openmp.f90 \
!   -o test/Stokes/int2-test-stokkernels-openmp
! test/Stokes/int2-test-stokkernels-openmp
program test_stokkernels_openmp
  use omp_lib
  implicit none
  integer(kind=8), parameter :: nd=2, ns=79, nt=83
  real(kind=8) :: source(3,ns), target(3,nt), strength(nd,3,ns,7)
  real(kind=8) :: u(nd,3,nt), p(nd,nt), g(nd,3,3,nt)
  real(kind=8) :: u1(nd,3,nt), p1(nd,nt), g1(nd,3,3,nt)
  integer(kind=8) :: i,j,k,l,istress,irotlet,idoublet,flag
  integer :: failures,batch
  character(len=16) :: mode
  call get_command_argument(1,mode)
  call omp_set_dynamic(.false.)
  do j=1,ns
    do i=1,3
      source(i,j)=sin(real(3*j+i,8))
      do k=1,nd
        do l=1,7
          strength(k,i,j,l)=cos(real(j+2*i+3*k+5*l,8))
        enddo
      enddo
    enddo
  enddo
  do j=1,nt
    target(:,j)=[cos(real(j,8)),sin(real(2*j,8)),cos(real(3*j,8))]
  enddo
  target(:,1)=source(:,1) ! Exercise threshold exclusion as well.
  do istress=0,1
    do irotlet=0,1
      do idoublet=0,1
        do flag=0,3
          call omp_set_num_threads(1)
          call evaluate(u1,p1,g1,istress,irotlet,idoublet,flag)
          call omp_set_num_threads(8)
          call evaluate(u,p,g,istress,irotlet,idoublet,flag)
          if (any(u/=u1) .or. any(p/=p1) .or. any(g/=g1)) &
            error stop 'One/eight-thread outputs differ'
          if (flag<1 .and. any(u/=0)) error stop 'Disabled velocity evaluated'
          if (flag<2 .and. any(p/=0)) error stop 'Disabled pressure evaluated'
          if (flag<3 .and. any(g/=0)) error stop 'Disabled gradient evaluated'
        enddo
      enddo
    enddo
  enddo
  if (mode/='--serial') then
  call omp_set_num_threads(1)
  call evaluate(u1,p1,g1,1_8,1_8,1_8,3_8)
  call omp_set_num_threads(8)
  call omp_set_max_active_levels(2)
  failures=0
  !$omp parallel do private(u,p,g) reduction(+:failures)
  do batch=1,8
    call evaluate(u,p,g,1_8,1_8,1_8,3_8)
    if (any(u/=u1) .or. any(p/=p1) .or. any(g/=g1)) failures=failures+1
  enddo
  !$omp end parallel do
  if (failures/=0) error stop 'Direct outputs differ inside parallel region'
  endif
  print *, 'PASS: 32 kernel/flag cases, two densities and threshold'
  if (mode/='--serial') print *, 'PASS: nested calls with active levels=2'
contains
  subroutine evaluate(velocity,pressure,gradient,stress,rot,doub,outputs)
    real(kind=8), intent(out) :: velocity(nd,3,nt),pressure(nd,nt),gradient(nd,3,3,nt)
    integer(kind=8), intent(in) :: stress,rot,doub,outputs
    velocity=0
    pressure=0
    gradient=0
    call st3ddirectstokstrsrotdoubg(nd,source,strength(:,:,:,1), &
      stress,strength(:,:,:,2),strength(:,:,:,3), &
      rot,strength(:,:,:,4),strength(:,:,:,5), &
      doub,strength(:,:,:,6),strength(:,:,:,7), &
      ns,target,nt,outputs,velocity,pressure,gradient,1d-16)
  end subroutine evaluate
end program test_stokkernels_openmp

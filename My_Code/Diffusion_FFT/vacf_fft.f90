!-------------------------------------------------------
!
! VACF_FFT: calculate Flux Autocorrelation Function via FFT
!
! By Yanlei Wang,@THU,Mar 2013
!-------------------------------------------------------
program vacf_fft
  use, intrinsic :: iso_c_binding
  include 'fftw3.f03'
  
  INTEGER :: ty,nsteps,l,m,k,initstep,n, allstep, natoms, ions
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: heat,heat_cent,z,mask,adj,acorr0,acorr,nor_acorr
  DOUBLE PRECISION, DIMENSION(:,:,:),ALLOCATABLE :: pos,vel,acov
  DOUBLE PRECISION, DIMENSION(:),ALLOCATABLE :: c,uc
  DOUBLE COMPLEX,DIMENSION(:),ALLOCATABLE :: out
  DOUBLE PRECISION :: mean,x,y,norm,time,tag
  DOUBLE PRECISION :: scaler
  TYPE(C_PTR) :: plan_f,plan_inv,plan_m,datain1,datain2,dataout1,dataout2,datain_m,dataout_m
 ! REAL(C_DOUBLE),POINTER :: in_m(:)
  COMPLEX(C_DOUBLE_COMPLEX),POINTER :: in1(:),out1(:),in2(:),in_m(:),out_m(:),out2(:)
     
  !--------------------------------------------
  ! preparing read,and obtain nsteps&natoms
  !--------------------------------------------
  OPEN (10, FILE='run.info', STATUS='OLD')
  read(10, *) natoms
  read(10, *) allstep
  read(10, *) ndump
  read(10 ,*) ions  ! 0 is cation, 1 is anion 
  nsteps=allstep/ndump

  OPEN (11,FILE='../mol.vel.xyz',STATUS='OLD')
  OPEN (12,FILE='0vacf.data',STATUS='UNKNOWN')
  OPEN (13,FILE='1nor_facf.data',STATUS='UNKNOWN')
  !find the nextpow2

  m=0
  l = 2*nsteps-1
  k=1
  do while( k .LT. l)
    k=2*k
    m=m+1
  end do
  l=1
  do i=1,m
   l=2*l
  end do
  write(*, *)  "Original step = ", nsteps
  write(*, *)  "2^n step = ", l 
  !---------------------------------------------
  ! allocate data from dump file
  !---------------------------------------------
  ALLOCATE(c(nsteps),uc(nsteps),acov(nsteps,natoms,3))
  allocate(vel(nsteps,natoms,3))
  allocate(nor_acorr(nsteps,3))
  datain1=fftw_alloc_complex(INT(l,C_SIZE_T))
  dataout1=fftw_alloc_complex(INT(l,C_SIZE_T))
!  datain2=fftw_alloc_complex(INT(l,C_SIZE_T))
!  dataout2=fftw_alloc_complex(INT(l,C_SIZE_T))
    
  CALL c_f_pointer(datain1,in1,[l])
  CALL c_f_pointer(dataout1,out1,[l])
!  CALL c_f_pointer(datain2,in2,[l])
!  CALL c_f_pointer(dataout2,out2,[l])
  !--------------------------------------------
  ! read data from dump file
  !--------------------------------------------
  REWIND(11)
  READ (11,*)
  READ (11,*)
  READ (11,*)

  DO i = 1, nsteps
    READ (11,*)
    if (ions .eq. 0) then  
      do j=1,natoms
        read(11,*) ty,vel(i,j,:)
      enddo
      do j=1,natoms
        read(11,*)
      enddo
    else 
      do j=1,natoms
        read(11,*)
      enddo
      do j=1,natoms
        read(11,*) ty,vel(i,j,:)
      enddo
    endif
  end do
 
  !--------------------------------------------
  ! calculate Flux Autocorrelation Function
  !--------------------------------------------
  !================================
  !calculate the mean of flux heat
  write(*,*) 1
  !creat new array

do n =1,natoms
  do j=1,3
    do i=1,nsteps
      in1(i) = vel(i,n,j) !The nth atom's vel at time i along direc. j
    end do
    do i=nsteps+1,l
      in1(i) = 0
    end do
    !forward FFT
  ! plan_f = fftw_plan_dft_r2c_1d(l,in1,out1, FFTW_ESTIMATE)
  ! call fftw_execute_dft_r2c(plan_f,in1,out1)  
    plan_f= fftw_plan_dft_1d(l,in1,out1,-1,FFTW_ESTIMATE)
    call fftw_execute_dft(plan_f,in1,out1)
    do i=1,l
      in1(i) = out1(i)
    end do  
    do i=1,l
      x=real(in1(i))
      y=aimag(in1(i))
      norm=x*x+y*y
      in1(i) = norm
    end do
    !Inverse FFT
  !  plan_inv = fftw_plan_dft_c2r_1d(l,in2,out2, FFTW_ESTIMATE)
  !  call fftw_execute_dft_c2r(plan_inv,in2,out2)
    plan_inv = fftw_plan_dft_1d(l,in1,out1,+1,FFTW_ESTIMATE)
    call fftw_execute_dft(plan_inv,in1,out1)
    do i=1,nsteps
      acov(i,n,j) = real(out1(i))/l
    end do
  end do
end do

do i = 1,nsteps
  sumcorr = 0.d0
  do j=1,3
    do n=1,natoms
     sumcorr = sumcorr+ acov(i,n,j)
    end do
  end do
  c(i) = sumcorr/natoms
end do

uc(1) = c(1)/DBLE(nsteps)
 do i=2,nsteps
   uc(i) = (c(i)/DBLE(nsteps+1-i))/uc(1)
 end do
 uc(1) =1.d0
  !--------------------------------------------
  ! export
  !--------------------------------------------
write(*,*) 8 
  !export the autocorrelation function
  write(12,*) "nsteps"
  write(12,*) nsteps
  
  !!! unit to ns
  scaler=1.0E12
  
  do i=1,nsteps
    write(12,'(I20,E20.10)') i,c(i)/DBLE(nsteps+1-i)*scaler
  end do
  close(12)
 
write(13,*) nsteps 
  do i=1,nsteps
    write(13,'(I20,E20.10)') i,uc(i)
  end do
  close(13) 
 end program vacf_fft

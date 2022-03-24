!--------------------------------------------------------------------------
! ozymandias:cooling_module.f90
!--------------------------------------------------------------------------
!
! MODULE: cooling_module
!
!> @author F. Rodriguez Montero
!
!> @brief 
!> types and routines useful for the reading of the cooling curve of
!> RAMSES.
!
!> @details  
!> define hydroID, amr_info, sim_info and level derived types, as well
!> as functions for the initial setup of AMR and HYDRO data.
! 
!
!> @date 22/3/2022   0.1 fortran routines for ozymandias
!--------------------------------------------------------------------------
module cooling_module
    use local
    use constants

    type cooling_table
        integer::n1
        integer::n2
        real(kind=8),dimension(:)    ,allocatable::nH
        real(kind=8),dimension(:)    ,allocatable::T2
        real(kind=8),dimension(:,:)  ,allocatable::cool
        real(kind=8),dimension(:,:)  ,allocatable::heat
        real(kind=8),dimension(:,:)  ,allocatable::cool_com
        real(kind=8),dimension(:,:)  ,allocatable::heat_com
        real(kind=8),dimension(:,:)  ,allocatable::metal
        real(kind=8),dimension(:,:)  ,allocatable::cool_prime
        real(kind=8),dimension(:,:)  ,allocatable::heat_prime
        real(kind=8),dimension(:,:)  ,allocatable::cool_com_prime
        real(kind=8),dimension(:,:)  ,allocatable::heat_com_prime
        real(kind=8),dimension(:,:)  ,allocatable::metal_prime
        real(kind=8),dimension(:,:)  ,allocatable::mu
        real(kind=8),dimension(:,:,:),allocatable::n_spec
    end type cooling_table

    type(cooling_table)::ctable
    ! Utilisation de table%n_spec si necessaire
    logical, parameter :: if_species_abundances=.true.
    logical, parameter :: self_shielding=.true.
    ! Les parametres de la table par defaut
    integer,parameter     :: nbin_T_fix=101
    integer,parameter     :: nbin_n_fix=161
    real(kind=8),parameter:: nH_min_fix=1.d-10
    real(kind=8),parameter:: nH_max_fix=1.d+6
    real(kind=8),parameter:: T2_min_fix=1.d-2
    real(kind=8),parameter:: T2_max_fix=1.d+9

    real(kind=8):: logT2max,dlog_nH,dlog_T2,h,h2,h3,precoeff
    contains
  
    subroutine read_cool(filename)
        implicit none
        character(LEN=80)::filename
        open(unit=10,file=filename,form='unformatted')
        ! Get size of table
        read(10)ctable%n1,ctable%n2

        ! Allocate cooling table
        if(.not.allocated(ctable%cool)) then
            allocate(ctable%cool(ctable%n1,ctable%n2))
            allocate(ctable%heat(ctable%n1,ctable%n2))
            allocate(ctable%cool_com(ctable%n1,ctable%n2))
            allocate(ctable%heat_com(ctable%n1,ctable%n2))
            allocate(ctable%metal(ctable%n1,ctable%n2))
            allocate(ctable%cool_prime(ctable%n1,ctable%n2))
            allocate(ctable%heat_prime(ctable%n1,ctable%n2))
            allocate(ctable%cool_com_prime(ctable%n1,ctable%n2))
            allocate(ctable%heat_com_prime(ctable%n1,ctable%n2))
            allocate(ctable%metal_prime(ctable%n1,ctable%n2))
            allocate(ctable%mu  (ctable%n1,ctable%n2))
            allocate(ctable%nH  (1:ctable%n1))
            allocate(ctable%T2  (1:ctable%n2))
        endif

        ! And read values
        read(10)ctable%nH
        read(10)ctable%T2
        read(10)ctable%cool
        read(10)ctable%heat
        read(10)ctable%cool_com
        read(10)ctable%heat_com
        read(10)ctable%metal
        read(10)ctable%cool_prime
        read(10)ctable%heat_prime
        read(10)ctable%cool_com_prime
        read(10)ctable%heat_com_prime
        read(10)ctable%metal_prime
        read(10)ctable%mu
        if (if_species_abundances) read(10)ctable%n_spec
        close(10)

        ! Some initialisations
        logT2max=log10(T2_max_fix)
        dlog_nH=dble(ctable%n1-1)/(ctable%nH(ctable%n1)-ctable%nH(1))
        dlog_T2=dble(ctable%n2-1)/(ctable%T2(ctable%n2)-ctable%T2(1))
        h=1d0/dlog_T2
        h2=h*h
        h3=h2*h
        precoeff=2d0*XH/(3d0*kBoltzmann)
    end subroutine read_cool

    subroutine solve_cooling(nH,T2,zsolar,lambda,lambda_prime)
        ! nH [H/cc], T2 [T/mu in Kelvin], Zsolar [metallicity in Zsun]
        implicit none
        real(kind=8),intent(in)::nH,T2,zsolar
        real(kind=8),intent(inout)::lambda,lambda_prime

        integer::i_nH,i_T2
        real(kind=8)::boost
        real(kind=8)::facT,dlog_nH,dlog_T2
        real(kind=8)::metal,cool,heat,cool_com,heat_com,w1T,w2T,w11,w12,w21,w22,err,yy,yy2,yy3
        real(kind=8)::metal_prime,cool_prime,heat_prime,cool_com_prime,heat_com_prime,wcool
        real(kind=8)::fa,fb,fprimea,fprimeb,alpha,beta,gamma
        real(kind=8)::rgt,lft,tau
        real(kind=8)::facH,zzz
        real(kind=8)::w1H,w2H

        ! Compute radiation boost factor
        if(self_shielding)then
            boost=MAX(exp(-nH/0.01),1.0D-20)
        else
            boost=1.0
        endif

        ! Get the necessary values for the cell
        zzz = zsolar
        facH = MIN(MAX(log10(nH/boost),ctable%nH(1)),ctable%nH(ctable%n1))
        i_nH = MIN(MAX(int((facH-ctable%nH(1))*dlog_nH)+1,1),ctable%n1-1)
        w1H = (ctable%nH(i_nH+1)-facH)*dlog_nH
        w2H = (facH-ctable%nH(i_nH))*dlog_nH
        tau = T2

        facT=log10(T2)

        if(facT.le.logT2max)then

           i_T2 = MIN(MAX(int((facT-ctable%T2(1))*dlog_T2)+1,1),ctable%n2-1)
           yy=facT-ctable%T2(i_T2)
           yy2=yy*yy
           yy3=yy2*yy

           ! Cooling
           fa=ctable%cool(i_nH,i_T2  )*w1H+ctable%cool(i_nH+1,i_T2  )*w2H
           fb=ctable%cool(i_nH,i_T2+1)*w1H+ctable%cool(i_nH+1,i_T2+1)*w2H
           fprimea=ctable%cool_prime(i_nH,i_T2  )*w1H+ctable%cool_prime(i_nH+1,i_T2  )*w2H
           fprimeb=ctable%cool_prime(i_nH,i_T2+1)*w1H+ctable%cool_prime(i_nH+1,i_T2+1)*w2H
           alpha=fprimea
           beta=3d0*(fb-fa)/h2-(2d0*fprimea+fprimeb)/h
           gamma=(fprimea+fprimeb)/h2-2d0*(fb-fa)/h3
           cool=10d0**(fa+alpha*yy+beta*yy2+gamma*yy3)
           cool_prime=cool/tau*(alpha+2d0*beta*yy+3d0*gamma*yy2)

           ! Heating
           fa=ctable%heat(i_nH,i_T2  )*w1H+ctable%heat(i_nH+1,i_T2  )*w2H
           fb=ctable%heat(i_nH,i_T2+1)*w1H+ctable%heat(i_nH+1,i_T2+1)*w2H
           fprimea=ctable%heat_prime(i_nH,i_T2  )*w1H+ctable%heat_prime(i_nH+1,i_T2  )*w2H
           fprimeb=ctable%heat_prime(i_nH,i_T2+1)*w1H+ctable%heat_prime(i_nH+1,i_T2+1)*w2H
           alpha=fprimea
           beta=3d0*(fb-fa)/h2-(2d0*fprimea+fprimeb)/h
           gamma=(fprimea+fprimeb)/h2-2d0*(fb-fa)/h3
           heat=10d0**(fa+alpha*yy+beta*yy2+gamma*yy3)
           heat_prime=heat/tau*(alpha+2d0*beta*yy+3d0*gamma*yy2)

           ! Compton cooling
           fa=ctable%cool_com(i_nH,i_T2  )*w1H+ctable%cool_com(i_nH+1,i_T2  )*w2H
           fb=ctable%cool_com(i_nH,i_T2+1)*w1H+ctable%cool_com(i_nH+1,i_T2+1)*w2H
           fprimea=ctable%cool_com_prime(i_nH,i_T2  )*w1H+ctable%cool_com_prime(i_nH+1,i_T2  )*w2H
           fprimeb=ctable%cool_com_prime(i_nH,i_T2+1)*w1H+ctable%cool_com_prime(i_nH+1,i_T2+1)*w2H
           alpha=fprimea
           beta=3d0*(fb-fa)/h2-(2d0*fprimea+fprimeb)/h
           gamma=(fprimea+fprimeb)/h2-2d0*(fb-fa)/h3
           cool_com=10d0**(fa+alpha*yy+beta*yy2+gamma*yy3)
           cool_com_prime=cool_com/tau*(alpha+2d0*beta*yy+3d0*gamma*yy2)

           ! Compton heating
           fa=ctable%heat_com(i_nH,i_T2  )*w1H+ctable%heat_com(i_nH+1,i_T2  )*w2H
           fb=ctable%heat_com(i_nH,i_T2+1)*w1H+ctable%heat_com(i_nH+1,i_T2+1)*w2H
           fprimea=ctable%heat_com_prime(i_nH,i_T2  )*w1H+ctable%heat_com_prime(i_nH+1,i_T2  )*w2H
           fprimeb=ctable%heat_com_prime(i_nH,i_T2+1)*w1H+ctable%heat_com_prime(i_nH+1,i_T2+1)*w2H
           alpha=fprimea
           beta=3d0*(fb-fa)/h2-(2d0*fprimea+fprimeb)/h
           gamma=(fprimea+fprimeb)/h2-2d0*(fb-fa)/h3
           heat_com=10d0**(fa+alpha*yy+beta*yy2+gamma*yy3)
           heat_com_prime=heat_com/tau*(alpha+2d0*beta*yy+3d0*gamma*yy2)

           ! Metal cooling
           fa=ctable%metal(i_nH,i_T2  )*w1H+ctable%metal(i_nH+1,i_T2  )*w2H
           fb=ctable%metal(i_nH,i_T2+1)*w1H+ctable%metal(i_nH+1,i_T2+1)*w2H
           fprimea=ctable%metal_prime(i_nH,i_T2  )*w1H+ctable%metal_prime(i_nH+1,i_T2  )*w2H
           fprimeb=ctable%metal_prime(i_nH,i_T2+1)*w1H+ctable%metal_prime(i_nH+1,i_T2+1)*w2H
           alpha=fprimea
           beta=3d0*(fb-fa)/h2-(2d0*fprimea+fprimeb)/h
           gamma=(fprimea+fprimeb)/h2-2d0*(fb-fa)/h3
           metal=10d0**(fa+alpha*yy+beta*yy2+gamma*yy3)
           metal_prime=metal/tau*(alpha+2d0*beta*yy+3d0*gamma*yy2)

           ! Total net cooling
           lambda=cool+zzz*metal-heat+(cool_com-heat_com)/nH
           lambda_prime=cool_prime+zzz*metal_prime-heat_prime+(cool_com_prime-heat_com_prime)/nH

        else

           lambda = 1.42*1d-27*sqrt(tau)*1.1
           lambda_prime = lambda/2./tau

        endif
    end subroutine solve_cooling 

end module cooling_module
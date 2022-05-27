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
    subroutine retrieve_table(repository,mytable)
        implicit none
        character(128), intent(in) :: repository
        type(cooling_table),intent(inout) :: mytable

        call read_cool(repository)
        mytable = ctable
    end subroutine retrieve_table
  
    subroutine read_cool(filename)
        implicit none
        character(LEN=128)::filename
        open(unit=15,file=TRIM(filename),form='unformatted')
        ! Get size of table
        read(15)ctable%n1,ctable%n2

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
            if (if_species_abundances) allocate(ctable%n_spec(ctable%n1,ctable%n2,1:6))
        endif

        ! And read values
        read(15)ctable%nH
        read(15)ctable%T2
        read(15)ctable%cool
        read(15)ctable%heat
        read(15)ctable%cool_com
        read(15)ctable%heat_com
        read(15)ctable%metal
        read(15)ctable%cool_prime
        read(15)ctable%heat_prime
        read(15)ctable%cool_com_prime
        read(15)ctable%heat_com_prime
        read(15)ctable%metal_prime
        read(15)ctable%mu
        if (if_species_abundances) read(15)ctable%n_spec
        close(15)

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
        real(kind=8)::metal,cool,heat,cool_com,heat_com
        real(kind=8)::metal_prime,cool_prime,heat_prime,cool_com_prime,heat_com_prime,wcool
        real(kind=8)::fa,fb,fprimea,fprimeb,alpha,beta,gamma
        real(kind=8)::rgt,lft,tau
        real(kind=8)::facH,zzz

        ! Compute radiation boost factor
        if(self_shielding)then
            boost=MAX(exp(-nH/0.01),1.0D-20)
        else
            boost=1.0
        endif

        ! Get the necessary values for the cell
        dlog_nH = dble(ctable%n1-1)/(ctable%nH(ctable%n1)-ctable%nH(1))
        dlog_T2 = dble(ctable%n2-1)/(ctable%T2(ctable%n2)-ctable%T2(1))
        zzz = zsolar
        facH = MIN(MAX(log10(nH/boost),ctable%nH(1)),ctable%nH(ctable%n1))
        i_nH = MIN(MAX(int((facH-ctable%nH(1))*dlog_nH)+1,1),ctable%n1-1)
        facT=log10(T2)
        i_T2 = MIN(MAX(int((facT-ctable%T2(1))*dlog_T2)+1,1),ctable%n2-1)
        tau = T2

        if(facT.le.logT2max)then
           ! Cooling
           cool=10d0**(ctable%cool(i_nH,i_T2  ))
           cool_prime=10d0**(ctable%cool_prime(i_nH,i_T2  ))

           ! Heating
           heat=10d0**(ctable%heat(i_nH,i_T2  ))
           heat_prime=10d0**(ctable%heat_prime(i_nH,i_T2  ))

           ! Compton cooling
           cool_com=10d0**(ctable%cool_com(i_nH,i_T2  ))
           cool_com_prime=10d0**(ctable%cool_com(i_nH,i_T2  ))

           ! Compton heating
           heat_com=10d0**(ctable%heat_com(i_nH,i_T2  ))
           heat_com_prime=10d0**(ctable%heat_com_prime(i_nH,i_T2  ))

           ! Metal cooling
           metal=10d0**(ctable%metal(i_nH,i_T2  ))
           metal_prime=10d0**(ctable%metal_prime(i_nH,i_T2  ))

           ! Total net cooling
           lambda=cool+zzz*metal-heat+(cool_com-heat_com)/nH
           lambda_prime=cool_prime+zzz*metal_prime-heat_prime+(cool_com_prime-heat_com_prime)/nH

        else

           lambda = 1.42*1d-27*sqrt(tau)*1.1
           lambda_prime = lambda/2./tau


        endif
    end subroutine solve_cooling 

end module cooling_module
!--------------------------------------------------------------------------
! ozymandias:sunset.f90
!--------------------------------------------------------------------------
!
! MODULE: sunset
!
!> @author R. Teyssier, 01/11/08
!
!> @brief 
!> types and routines required to construct SUNSET images
!
!> @details  
!> Sunset is a program that computes realistic galaxy images using spectra
! of a single stellar population (SSP) of your favorite model (GALAXEV, STARDUST, PEGASE).
!>
!> The images are computed for various filters (U,V,B,R,I,J...) and take into account,
!> if desired, dust attenuation. In order to work, you need to provide the following data:
!> - a particle file containing stars, with position, velocity, birth date and metallicity
!> - a SSP spectral 3D cube (age, metallicity, lambda)
!> _ the spectral response datebase for each filer
!> If you need dust absorption, you also need:
!> - a dust opacity model (kappa g/cm2 versus lambda) 
!> - a dust mass density 3D array covering the same bounding box than the image
!> The output data are pretty standard stuff.
! 
!
!> @date Created v1.0: sunset, S. Martin-Alvarez ??/??/??
!> @date Updated v1.1 (10/08/17): allowing dust for off-axes projections
!> @date Upgrade v2.0 (16/01/20): module definitions in new files
!> @date Upgrade v2.1 (19/07/22): adapting sunset for the Ozymandias integration
!--------------------------------------------------------------------------

module sunset
    use local
    use constants
    use vectors
    use rotations
    use io_ramses
    use geometrical_regions
    use obs_instruments
    
    ! Emission parameters
    integer,parameter::NDUSTMAX=10000
    integer,parameter::nmaxlambda=2800
    integer,parameter::NMAXFILTERS=100
    integer,parameter::nmaxlines=100
    integer,parameter::nmaxtimes=100 

    ! Legacy variables
    integer :: nfilters
    integer :: ifilter
    integer::ilambda_min,ilambda_max
    real(dbl),dimension(1:NDUSTMAX)::ldust,kdust
    real(dble),dimension(1:NMAXFILTERS,1:1000)::lambdafilter,transfilter
    integer,dimension(1:NMAXFILTERS)::nlambdafilter,typetrans,typecalib
    character(LEN=80)::ordering,str
    character(LEN=128)::sspfile='/mnt/zfsusers/currodri/Codes/ozymandias/ozy/visualisation/DataFiles/bruzual_charlot_2003.dat'
    character(LEN=128)::filterfile='/mnt/zfsusers/currodri/Codes/ozymandias/ozy/visualisation/DataFiles/filters.dat'
    character(LEN=128)::dustfile='/mnt/zfsusers/currodri/Codes/ozymandias/ozy/visualisation/DataFiles/kext_albedo_WD_MW_3.1_60.txt'
    character(LEN=10),dimension(1:NMAXFILTERS)::fname

    ! Working parameters
    logical :: shiftFilter
    logical :: do_dust
    integer::maxcosmol
    logical::cosmo_levels=.false.
    logical :: maxram_levels=.true.
    integer :: nhist=10000
    integer :: nlambda,ntime,nmetal
    character(LEN=128)::interp_mode='nearest'

    ! Reading files variables
    character(len=128)::repository,repository1,repository2
    integer :: outputID

    ! Common arrays
    real(dbl),dimension(1:3,1:3)::trans_matrix
    real(dbl),dimension(:),allocatable :: lssp,tssp,mssp,sed,filter_contrib,kappa_dust
    real(dbl),dimension(:),allocatable :: time_hist,metal_hist
    real(dbl),dimension(:,:),allocatable::intssp
    real(dbl),dimension(:,:,:),allocatable::ssp,column
    real(dbl),dimension(:,:,:),allocatable::rho

    ! Camera/Observer variables
    type(region) :: map_box

    contains

    !---------------------------------------------------------------
    ! Subroutine: READ TRANSMISSION FILTERS
    !
    ! Read filter transmission from PEGASE file.
    !---------------------------------------------------------------
    subroutine read_transmission_filters(filtername,zeta)
        implicit none
        character(LEN=10) :: filtername
        real(dbl) :: zeta
        integer :: i,j
        write(*,'(": Read filters for the following bands")')
        open(30,status='old',file=filterfile)
        read(30,*)nfilters
        ifilter = 1
        do i=1, nfilters
            read(30,*)nlambdafilter(i),typetrans(i),typecalib(i),fname(i)
            if(TRIM(fname(i))==TRIM(filtername)) ifilter=i
            do j=1, nlambdafilter(i)
                read(30,*)lambdafilter(i,j),transfilter(i,j)
            end do
        end do
        close(30)
        ! Shift the read filter to the redshift of the source
        if (shiftFilter) then       
            do j=1,nlambdafilter(ifilter)
                lambdafilter(ifilter,j)=lambdafilter(ifilter,j)/(1.0d0+zeta)
            end do
        else
            write(*,'(": Filter redshift displacement deactivated")')
        end if
        write(*,'(": Working band ID: " ,I3,". Name: ",a7)') ifilter, TRIM(filtername)
    end subroutine read_transmission_filters

    !---------------------------------------------------------------
    ! Subroutine: READ DUST OPPACITIES
    !
    ! Read dust opacity file from Draine
    ! Units are micron and cm2/g
    !---------------------------------------------------------------
    subroutine read_dust_opacities
        implicit none
        logical :: err
        integer :: n,ndust
        real(dbl) :: lambda,dum,kappa

        open(1,file=dustfile,form='formatted',status='old')
        err=.true.
        do while(err)
            read(1,'(A80)',END=81)str
            if(TRIM(str)==' lambda   albedo    g     C_ext/H    K_abs')err=.false.
        end do
        read(1,*)
        read(1,*)
        err=.true.
        n=0
        do while(err)
            read(1,*,END=81)lambda,dum,dum,dum,kappa
            n=n+1
            ldust(n)=lambda
            kdust(n)=kappa
        end do
81      continue
        close(1)
        ndust=n
        ! Convert microns to Angstroms
        ldust(1:ndust)=ldust(1:ndust)*1d4
    end subroutine read_dust_opacities

    !---------------------------------------------------------------
    ! Subroutine: READ SSP FILE
    !
    ! Read SSP file (STARDUST 99 format)
    ! The mass of the SSP is 10^6 Msol.
    ! Fluxes are in erg/s/A.
    ! lambda are in Angstroms, 
    ! ages in years, 
    ! metallicity in solar units
    !---------------------------------------------------------------
    subroutine read_ssp_file
        implicit none

        write(*,'(": Reading spectra from Bruzual-Charlot 2003")')
        open(1,file=sspfile,form='unformatted',status='old')
        read(1)nlambda,ntime,nmetal
        allocate(lssp(1:nlambda))
        allocate(sed (1:nlambda))
        allocate(filter_contrib(1:nlambda))
        sed = 0d0
        allocate(tssp(1:ntime))
        allocate(time_hist(0:nhist))
        allocate(metal_hist(0:nhist))
        time_hist=0d0
        metal_hist=0d0
        allocate(mssp(1:nmetal))
        allocate(ssp(1:ntime,1:nmetal,1:nlambda))
        read(1)lssp
        read(1)tssp
        read(1)mssp
        read(1)ssp
        close(1)
    end subroutine read_ssp_file

    !---------------------------------------------------------------
    ! Subroutine: CHECK INPUT FILES
    !
    ! This routine checks that all required input files are in the
    ! expected location.
    !---------------------------------------------------------------
    subroutine check_input_files(check_hydro,check_info,check_amr,check_part)
        implicit none
        logical,intent(in)::check_hydro,check_info,check_amr,check_part

        integer :: ipos
        logical :: ok
        character(len=5)::nchar,ncharcpu,nchar1,nchar2
        character(len=128)::nomfich
        ! Check initialisations
        ipos=INDEX(repository,'output_')
        nchar1=repository1(ipos+7:ipos+13)
        nchar2=repository2(ipos+7:ipos+13)
        nchar=repository(ipos+7:ipos+13)
        read(nchar,*) outputID

        ! Check info file
        if (check_info) then
            nomfich=TRIM(repository)//'/info_'//TRIM(nchar)//'.txt'
            inquire(file=nomfich, exist=ok) ! verify input file
            if ( .not. ok ) then
                write(*,'(": Info file not found: ", A)') trim(nomfich)
                stop
            end if
        end if

        ! Check hydro file
        if (check_hydro) then
            nomfich=TRIM(repository)//'/hydro_'//TRIM(nchar)//'.out00001'
            inquire(file=nomfich, exist=ok) ! verify input file
            if ( .not. ok ) then
                write(*,'(": Hydro file not found: ", A)') trim(nomfich)
                stop
            end if
        end if

        ! Check amr file
        if (check_amr) then
            nomfich=TRIM(repository)//'/amr_'//TRIM(nchar)//'.out00001'
            inquire(file=nomfich, exist=ok) ! verify input file
            if ( .not. ok ) then
                write(*,'(": AMR file not found: ", A)') trim(nomfich)
                stop
            end if
        end if

        ! Check particle file
        if (check_part) then
            nomfich=TRIM(repository)//'/part_'//TRIM(nchar)//'.out00001'
            inquire(file=nomfich, exist=ok) ! verify input file
            if ( .not. ok ) then
                write(*,'(": Particle file not found: ", A)') trim(nomfich)
                stop
            end if
        end if
    end subroutine check_input_files

    !---------------------------------------------------------------
    ! Subroutine: COMPUTE EFFECTIVE LMAX
    !
    ! This subroutine computes the lmax to be used in the 
    ! calculations.
    !---------------------------------------------------------------
    subroutine compute_effective_lmax
        implicit none
        logical,parameter::re_up=.true.
        ! Maxram levels variables
        character(len=5)::nchar,ncharcpu
        character(len=128)::nomfich
        integer,dimension(:,:),allocatable::lev_ngridf,lev_ngridl
        integer::icpu,i
        integer::klevel,klevmax
        integer::allgridnum
        integer::lev_ncpu,lev_nbound,lev_nlevelmax
        logical::amrheader_read
        ! Fix lmax to default
        if(amr%lmax==0) amr%lmax=amr%nlevelmax
        if(amr%lmax.gt.amr%nlevelmax) amr%lmax=amr%nlevelmax
        if (maxram_levels) then       ! RAMSES entire grid definer
            klevmax=1
            amrheader_read=.false.
            CheckCPU: do icpu=1,amr%ncpu
                ! Open AMR file and skip header ###
                call title(icpu,ncharcpu)
                nomfich=trim(repository)//'/amr_'//trim(nchar)//'.out'//trim(ncharcpu)
                open(unit=10,file=nomfich,status='old',form='unformatted')
                ! First read initialises all the required variables
                if(amrheader_read) then
                    do i=1,21
                        read(10)
                    end do
                else
                    read(10)lev_ncpu
                    read(10)!ndim
                    read(10)!nx,ny,nz
                    read(10)lev_nlevelmax
                    read(10)!ngridmax
                    read(10)lev_nbound
                    do i=1,15
                        read(10)
                    end do
                    allocate(lev_ngridf(1:lev_ncpu+lev_nbound,1:lev_nlevelmax))
                    allocate(lev_ngridl(1:lev_ncpu,1:lev_nlevelmax))
                    amrheader_read=.true.
                end if
                ! Read grid numbers
                read(10)lev_ngridl
                lev_ngridf(1:lev_ncpu,1:lev_nlevelmax)=lev_ngridl
                close(10)
                ! Closed file for future use ######

                ! Review grid numbers in each level to find highest refinement in sim
                do klevel=1,amr%lmax
                    allgridnum = sum(lev_ngridf(:,klevel))
                    if ((allgridnum.gt.0).and.(klevel.gt.klevmax)) klevmax=klevel
                end do
            end do CheckCPU
            ! Reduce maximum level to highest AMR level with existing grids
            if (klevmax.lt.amr%lmax) then
                write(*,'(": Lmax ",I2.1," reduced to ",I2.1," by maxram")') amr%lmax, klevmax
                amr%lmax=klevmax
            else
                write(*,'(": Lmax (maxram) =    ",I7.1,". Lmaxsim    = ",I9.1)') amr%lmax, klevmax
            end if
            deallocate(lev_ngridf,lev_ngridl)
        else if (cosmo_levels) then   ! RAMSES cosmolevels definer
            maxcosmol=amr%nlevelmax-int(log(0.79/sim%aexp)/log(2.)+1)
            if (maxcosmol .lt. amr%lmax) then
                write(*,'(": Lmax ",I2.1," reduced to ",I2.1," by expansion")') amr%lmax, maxcosmol
                amr%lmax = maxcosmol
            else
                write(*,'(": Lmax (csmlev) =    ",I7.1)') amr%lmax
            end if
            if (re_up) then
                amr%lmax = amr%lmax + 1
                write(*,'(": Increasing lmax back 1 level up:",I3)') amr%lmax
            end if
        else
            write(*,'(": Maxlevel finder inactive, lmax=",I2.1)') amr%lmax
        end if
    end subroutine compute_effective_lmax
    !---------------------------------------------------------------
    ! Subroutine: PREPARE CODE LUMINOSITY
    !
    ! Convert log lum into flux divided by code mass unit.
    !---------------------------------------------------------------
    subroutine prepare_code_luminosity
        implicit none
        integer :: ilambda,itime,imetal
        real(kind=8)::qfac
        real(kind=8)::Area10pc
        qfac=sim%unit_m/1d6 ! The 1d6 is required due for STARDUST99 format: Mass_SSP=10^6 Msun
        Area10pc =(4.0d0*pi*(10.0*pc2cm)**2)
        do itime=1,ntime
            do imetal=1,nmetal
                do ilambda=1,nlambda
                    ! 10 parsec luminosity (abs. magnitude)
                    ssp(itime,imetal,ilambda)=qfac*10d0**ssp(itime,imetal,ilambda)/Area10pc 
                end do
            end do
        end do
    end subroutine prepare_code_luminosity
    !---------------------------------------------------------------
    ! Subroutine: COMPUTE DUST SED
    !
    ! Compute dust opacity at the SED wavelength.
    !---------------------------------------------------------------
    subroutine compute_dust2SED
        implicit none
        integer :: ilambda,j
        real(dbl) :: lll

        allocate(kappa_dust(1:nlambda))

        ilambda_min=1
        do ilambda=1,nlambda
            j=ilambda_min
            do while(ldust(j)<lssp(ilambda))
                j=j+1
            end do
            lll=lssp(ilambda)
            if(j==1)then
                kappa_dust(ilambda)=0d0
            else
                kappa_dust(ilambda)=kdust(j-1)*(lll-ldust(j))/(ldust(j-1)-ldust(j))+ &
                    & kdust(j)*(lll-ldust(j-1))/(ldust(j)-ldust(j-1))
            end if
            ilambda_min=j
        end do
    end subroutine compute_dust2SED
    !---------------------------------------------------------------
    ! Subroutine: COMPUTE SSP FILTER
    !
    ! Compute ssp emission in the chosen filter.
    !---------------------------------------------------------------
    subroutine compute_SSP2filter
        implicit none
        integer :: iii,ilambda,itime,imetal
        real(dbl) :: lll,fff
        real(dbl) :: normfilter,nu

        allocate(intssp(1:ntime,1:nmetal))
        intssp=0d0
        ilambda_min=1
        do while(lssp(ilambda_min)<lambdafilter(ifilter,1))
            ilambda_min=ilambda_min+1
        end do
        ilambda_min=ilambda_min-1
        ilambda_max=nlambda
        do while(lssp(ilambda_max)>lambdafilter(ifilter,nlambdafilter(ifilter)))
            ilambda_max=ilambda_max-1
        end do
        ilambda_max=ilambda_max+1
        
        ! Normalize filter
        do ilambda=ilambda_min,ilambda_max
            lll=lssp(ilambda)
            iii=int(dble(nlambdafilter(ifilter)-1)*(lll-lambdafilter(ifilter,1))/ &
                    & (lambdafilter(ifilter,nlambdafilter(ifilter))-lambdafilter(ifilter,1)))
            iii=iii+1
            iii=max(iii,1)
            iii=min(iii,nlambdafilter(ifilter)-1)
            fff=transfilter(ifilter,iii)*(lll-lambdafilter(ifilter,iii+1)) &
                    & /(lambdafilter(ifilter,iii)-lambdafilter(ifilter,iii+1))+ &
                    & transfilter(ifilter,iii+1)*(lll-lambdafilter(ifilter,iii)) &
                    & /(lambdafilter(ifilter,iii+1)-lambdafilter(ifilter,iii))
            fff=max(fff,0d0)
            normfilter=normfilter+fff*(lssp(ilambda+1)-lssp(ilambda))*2.0/(lssp(ilambda+1)+lssp(ilambda))
        end do

        do itime=1,ntime
            do imetal=1,nmetal
                do ilambda=ilambda_min,ilambda_max
                    lll=lssp(ilambda)
                    iii=int(dble(nlambdafilter(ifilter)-1)*(lll-lambdafilter(ifilter,1))/ &
                        & (lambdafilter(ifilter,nlambdafilter(ifilter))-lambdafilter(ifilter,1)))
                    iii=iii+1
                    iii=max(iii,1)
                    iii=min(iii,nlambdafilter(ifilter)-1)
                    fff=transfilter(ifilter,iii)*(lll-lambdafilter(ifilter,iii+1)) &
                        & /(lambdafilter(ifilter,iii)-lambdafilter(ifilter,iii+1))+ &
                        & transfilter(ifilter,iii+1)*(lll-lambdafilter(ifilter,iii)) &
                        & /(lambdafilter(ifilter,iii+1)-lambdafilter(ifilter,iii))
                    fff=max(fff,0d0)
                    nu=clight/Armstrong2cm*2.0/(lssp(ilambda+1)+lssp(ilambda))
                    intssp(itime,imetal)=intssp(itime,imetal)+fff*ssp(itime,imetal,ilambda) &
                        & *(lssp(ilambda+1)-lssp(ilambda))/normfilter/nu
                    filter_contrib(ilambda)=fff*(lssp(ilambda+1)-lssp(ilambda))/normfilter/nu
                end do
            end do
        end do
    end subroutine compute_SSP2filter
    !---------------------------------------------------------------
    ! Subroutine: READ DUST CUBE
    !
    ! Read the provided dust mass density array.
    ! If not present, make it!
    !---------------------------------------------------------------
    subroutine read_dust_cube
        use export_amr
        use obs_instruments
        implicit none
        logical :: ok
        integer :: i,j,k,ii,jj,kk
        type(chunk_handler) :: chunk
        type(region) :: bbox,temp_bbox
        real(dbl) :: max_column
        real(dbl) :: xx,yy,zz
        type(vector) :: pos

        inquire(file=rhofile, exist=ok) ! verify input file
        if (.not. ok) then
            write(*,'(": Dust file missing, so computing and saving in ",A)') TRIM(rhofile)
            ! Setting up chunk
            chunk%nx=nx;chunk%ny=ny;chunk%nz=nz
            chunk%nvars = 1
            chunk%filt = filt_hydro
            call allocate_chunk_handler(chunk)
            chunk%varnames(1) = 'dust_density'
            ! Extracting unigrid cube
            call get_bounding_box(cam,bbox)
            call get_unigrid(repository,bbox,amr%lmax,.false.,chunk)
            allocate(rho(1:n1,1:n2,1:n3))
            rho = chunk%data(1,:,:,:)
            ! Save cube to file
            open(unit=10,file=rhofile,form='unformatted')
            write(10)chunk%nx,chunk%ny,chunk%nz
            write(10)bbox%xmin,bbox%xmax,bbox%ymin,bbox%ymax,bbox%zmin,bbox%zmax
            write(10)chunk%data(1,:,:,:)
            close(10)
        else
            ! Read data cube
            write(*,'(": Reading dust density in file ",A)') TRIM(rhofile)
            open(unit=20,file=rhofile,form='unformatted')
            read(20)n1,n2,n3
            read(20)temp_bbox%xmin,temp_bbox%xmax,temp_bbox%ymin,temp_bbox%ymax,temp_bbox%zmin,temp_bbox%zmax
            
            ! Check whether this box bounds the rotated region
            call get_bounding_box(cam,bbox)
            if(  (temp_bbox%xmin.ne.bbox%xmin).or.&
                &(temp_bbox%xmax.ne.bbox%xmax).or.&
                &(temp_bbox%ymin.ne.bbox%ymin).or.&
                &(temp_bbox%ymax.ne.bbox%ymax).or.&
                &(temp_bbox%zmin.ne.bbox%zmin).or.&
                &(temp_bbox%zmax.ne.bbox%zmax)) then
                write(*,'(": Dust file limits do not contain the required region ",A)') TRIM(rhofile)
                write(*,'(": Recomputing dust file… ",A)') TRIM(rhofile)
                close(20)
                ! Setting up chunk
                chunk%nx=nx;chunk%ny=ny;chunk%nz=nz
                chunk%nvars = 1
                chunk%filt = filt_hydro
                call allocate_chunk_handler(chunk)
                chunk%varnames(1) = 'dust_density'
                ! Extracting unigrid cube
                call get_bounding_box(cam,bbox)
                call get_unigrid(repository,bbox,amr%lmax,.false.,chunk)
                allocate(rho(1:n1,1:n2,1:n3))
                rho = chunk%data(1,:,:,:)
                ! Save cube to file
                open(unit=10,file=rhofile,form='unformatted')
                write(10)chunk%nx,chunk%ny,chunk%nz
                write(10)bbox%xmin,bbox%xmax,bbox%ymin,bbox%ymax,bbox%zmin,bbox%zmax
                write(10)chunk%data(1,:,:,:)
                close(10)
            else
                allocate(rho(1:n1,1:n2,1:n3))
                read(20)rho
                close(20)
            end if
        end if

        ! Now transform into the coordinates of the camera 
        ! and compute actual column density
        write(*,'(": Interpolating into rotated cube using interpolation method ",A)')TRIM(interp_mode)
        allocate(column(1:nx,1:ny,1:nz))
        call los_transformation(cam,trans_matrix)
        do i=1,nx
            xx = dble(i) * (cam%region_size(1)/nx) 
            do j=1,ny
                yy = dble(j) * (cam%region_size(2)/ny)
                do k=1,nz
                    zz = dble(k) * (cam%far_cut_depth + cam%distance)/nz
                    pos = (/xx,yy,zz/)
                    call rotate_vector(pos,transpose(trans_matrix))
                    pos = pos + cam%centre
                    ii = int(dble(nx)*(pos%x-bbox%xmin)/(bbox%xmax-bbox%xmin)) + 1
                    jj = int(dble(ny)*(pos%y-bbox%ymin)/(bbox%ymax-bbox%ymin)) + 1
                    kk = int(dble(nz)*(pos%z-bbox%zmin)/(bbox%zmax-bbox%zmin)) + 1
                    if (TRIM(interp_mode).eq.'nearest') then
                        column(i,j,k) = rho(ii,jj,kk)
                    end if
                end do
            end do
        end do

        write(*,'(": Computing column density")')
        max_column = 0d0
        column(:,:,1) = column(:,:,1) * dzpix
        do i=1,nx
            do j=1,ny
                do k=2,nz
                    column(i,j,k) = column(i,j,k-1)+column(i,j,k) * dzpix
                    max_column = max(max_column,column(i,j,k))
                end do
            end do
        end do
        column = column * (sim%unit_d * sim%unit_l)
        write(*,'(": Maximum dust column density ",E13.5)') max_column * (sim%unit_d * sim%unit_l)/1.66d-24
        
    end subroutine read_dust_cube
    !---------------------------------------------------------------
    ! Subroutine: COMPUTE FILTER IMAGE
    !
    ! Compute image in the chosen filter while reading
    ! particles in chunks
    !---------------------------------------------------------------
    subroutine compute_filter_image
        implicit none
        integer::npart_proc
        integer::icpu
        real(dbl)::paint_val
        real(dbl)::ageYr,tbornYr
        
        mtot = 0.0d0
        npart_proc = 0

        cpuloop: do k=1,amr%ncpu_read
            icpu = amr%cpu_list(k)
            call title(icpu,ncharcpu)
            call check_stars
            if (.not.okstars) cycle
            call get_partdata

            if (TRIM(kernel_type).eq."proximity") then
                call measure_particles_proximity
            end if
            dust_opacity = 1.0d0
            StarsLoop: do i=1,npart2
                call check_okpart(parts(i),okpart)
                weight = 1.0d0

                ! For the required stellar particles
                ProjectStar: if (okpart) then
                    npart_proc = npart_proc + 1
                    call getpartvalue(bbox,parts(i),'age',ageYr,dx)
                    call getpartvalue(bbox,parts(i),'birth_date',tbornYr,dx)
                    ageYr = ageYr * 1d9; tbornYr = tbornYr * 1d9
                    
                end if ProjectStar
            end do StarsLoop
        end do cpuloop

    end subroutine compute_filter_image
end module sunset
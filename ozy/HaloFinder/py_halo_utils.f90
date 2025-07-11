module py_halo_utils

public

contains

  
  !*************************************************************************
  subroutine get_nb_halos(file,nhalo,nsubhalo)
    
    implicit none 

    character(1000),intent(in)  :: file
    integer(kind=4),intent(out) :: nhalo
    integer(kind=4),intent(out) :: nsubhalo

    open(unit=11,file=file,status='old',form='unformatted',action='read')
    read(11) !nbodies_full_sim
    read(11) !massp
    read(11) !aexp
    read(11) !omega_t
    read(11) !age_univ
    read(11) nhalo,nsubhalo
    close(11)

    return

  end subroutine get_nb_halos

  !*************************************************************************
  subroutine read_all_halos_with_contam(file,halos,n)

    implicit none 

    character(1000),intent(in) :: file
    integer(kind=4),intent(in) :: n
    real(kind=4),intent(out)   :: halos(n,33)
    
    integer(kind=4) :: i,nh,ns
    integer(kind=4) :: np,num,ts,lev,hosth,hosts,nbsub,nextsub,contam
    real(kind=4)    :: mass,x,y,z,vx,vy,vz,lx,ly,lz,r,a,b,c,ek,ep,et,spin,rvir,mvir,tvir,cvel,rho_0,r_c

    open(unit=11,file=file,status='old',form='unformatted')
    read(11) !nbodies_full_sim
    read(11) !massp
    read(11) !aexp
    read(11) !omega_t
    read(11) !age_univ
    read(11) nh,ns
    if (n > nh+ns) stop
    do i = 1,n
       read(11) np
       read(11) ! skip part id's
       read(11) num
       read(11) ts
       read(11) lev,hosth,hosts,nbsub,nextsub
       read(11) mass
       read(11) x,y,z
       read(11) vx,vy,vz
       read(11) lx,ly,lz
       read(11) r,a,b,c
       read(11) ek,ep,et
       read(11) spin
       read(11) rvir,mvir,tvir,cvel
       read(11) rho_0,r_c
       read(11) contam
       halos(i,:) = (/real(np,4),real(num,4),real(ts,4),real(lev,4),real(hosth,4),real(hosts,4),real(nbsub,4),&
            & real(nextsub,4),mass,x,y,z,vx,vy,vz,lx,ly,lz,r,a,b,c,&
            & ek,ep,et,spin,rvir,mvir,tvir,cvel,rho_0,r_c,real(contam,4)/)
    end do
    close(11)

    return
    
  end subroutine read_all_halos_with_contam
  
  !*************************************************************************
  subroutine read_all_halos(file,halos,n)

    implicit none 

    character(1000),intent(in) :: file
    integer(kind=4),intent(in) :: n
    real(kind=4),intent(out)   :: halos(n,32)
    
    integer(kind=4) :: i,nh,ns
    integer(kind=4) :: np,num,ts,lev,hosth,hosts,nbsub,nextsub,contam
    real(kind=4)    :: mass,x,y,z,vx,vy,vz,lx,ly,lz,r,a,b,c,ek,ep,et,spin,rvir,mvir,tvir,cvel,rho_0,r_c

    open(unit=11,file=file,status='old',form='unformatted')
    read(11) !nbodies_full_sim
    read(11) !massp
    read(11) !aexp
    read(11) !omega_t
    read(11) !age_univ
    read(11) nh,ns
    if (n > nh+ns) stop
    do i = 1,n
       read(11) np
       read(11) ! skip part id's
       read(11) num
       read(11) ts
       read(11) lev,hosth,hosts,nbsub,nextsub
       read(11) mass
       read(11) x,y,z
       read(11) vx,vy,vz
       read(11) lx,ly,lz
       read(11) r,a,b,c
       read(11) ek,ep,et
       read(11) spin
       read(11) rvir,mvir,tvir,cvel
       read(11) rho_0,r_c
       halos(i,:) = (/real(np,4),real(num,4),real(ts,4),real(lev,4),real(hosth,4),real(hosts,4),real(nbsub,4),&
            & real(nextsub,4),mass,x,y,z,vx,vy,vz,lx,ly,lz,r,a,b,c,&
            & ek,ep,et,spin,rvir,mvir,tvir,cvel,rho_0,r_c/)
    end do
    close(11)

    return

  end subroutine read_all_halos

  !*************************************************************************
  subroutine read_all_galaxies(file,halos,n)

    ! Galaxy maker version of read_all_halos

    implicit none 

    character(1000),intent(in) :: file
    integer(kind=4),intent(in) :: n
    real(kind=4),intent(out)   :: halos(n,32)
    
    integer(kind=4) :: i,nh,ns
    integer(kind=4) :: np,num,ts,lev,hosth,hosts,nbsub,nextsub,contam
    real(kind=4)    :: mass,r,a,b,c,ek,ep,et,spin,rvir,mvir,tvir,cvel,rho_0,r_c
    real(kind=8)    :: x,y,z,vx,vy,vz,lx,ly,lz

    open(unit=11,file=file,status='old',form='unformatted')
    read(11) !nbodies_full_sim
    read(11) !massp
    read(11) !aexp
    read(11) !omega_t
    read(11) !age_univ
    read(11) nh,ns
    if (n > nh+ns) stop
    do i = 1,n
       read(11) np
       read(11) ! skip part id's
       read(11) ! skip stellar particle mass
       read(11) ! skip stellar particle age
       read(11) ! skip stellar particle metallicity
       read(11) ! skip part id's
       read(11) num
       read(11) ts
       read(11) lev,hosth,hosts,nbsub,nextsub
       read(11) mass
       read(11) x,y,z
       read(11) vx,vy,vz
       read(11) lx,ly,lz
       read(11) r,a,b,c
       read(11) ek,ep,et
       read(11) spin
       read(11) rvir,mvir,tvir,cvel
       read(11) rho_0,r_c
       halos(i,:) = (/real(np,4),real(num,4),real(ts,4),real(lev,4),real(hosth,4),real(hosts,4),real(nbsub,4),&
            & real(nextsub,4),mass,real(x,4),real(y,4),real(z,4),real(vx,4),real(vy,4),real(vz,4), &
            & real(lx,4),real(ly,4),real(lz,4),r,a,b,c,&
            & ek,ep,et,spin,rvir,mvir,tvir,cvel,rho_0,r_c/)
    end do
    close(11)

    return

  end subroutine read_all_galaxies

  !***********************************************************************
  subroutine get_sfrs(starfile, lookback_myr, nGals, sfrs)

  ! Extract star formation rates for each galaxy and return in Msun/yr
  ! Parameters:
  ! starfile     => File containing all stellar particles per galaxy
  ! nGals        => Number of galaxies
  ! lookback_myr => lookback time for measuring SFR, in Myr
  ! sfrs         <= SFRs [Msun/yr] 
  !-----------------------------------------------------------------------
    implicit none 
    character(1000),intent(in) :: starfile
    integer,intent(in) :: nGals
    real,intent(in) :: lookback_myr
    real(kind=4),intent(out)   :: sfrs(nGals)
    integer::iGal, galNr, nGals2, iStar, nStars
    real(kind=8)::mTot
    real::sAge, sMass, sZ, sPosx, sPosy, sPosz, sMassInit
  !-----------------------------------------------------------------------
    ! Open file generated by halofinder containing stars in each galaxy
    open(unit=11,file=starfile,status='old',form='formatted')
    read(11, '(I12)') nGals2                        !   Number of galaxies
    if(nGals2 .ne. nGals) then
       print*,'nGals!=nGals2!!!'
       print*,'nGals = ',nGals
       print*,'nGals2 = ',nGals2
       return
    endif

    sfrs(:) = 0.
    do iGal=1,nGals
       write (*, "(A, I5, A, I5, A)", advance='no') &
            ' HaloFinder/f90/py_halo_utils.f90: Processing galaxy '     &
            ,iGal,' / ',nGals, char(13)
       read(11, '(I12, I12)') galNr, nStars
       mTot=0d0
       do iStar=1,nStars
          read(11, '(7(e14.6,1x))') &
               sAge, sMass, sZ, sPosx, sPosy, sPosz, sMassInit
          if(sAge .le. 0d0) then
             print*,'Negative age in HalFinder/f90/py_halo_utils.f90->get_sfrs() = ',sAge
             sAge = 1d-8
          endif
          sAge = sAge/1d6                                 !     age in Myr
          if(sAge .lt. lookback_myr) then
             sMassInit = sMassInit*1d11                   !   mass in Msun
             sfrs(iGal) = sfrs(iGal) + sMassInit
          endif
       end do ! Loop over stars
    end do ! Loop over galaxies
    write(*,*) char(13)
    close(11)
    sfrs(:) = sfrs(:) / lookback_myr / 1d6 ! SFR in Msun/yr

  end subroutine get_sfrs

  !***********************************************************************
  subroutine get_ages_metallicities(starfile, nGals, ages_Myr, Zs)

  ! Extract mass-weighted average stellar ages and metalicities for
  ! galaxies  and return
  ! Parameters:
  ! starfile     => File containing all stellar particles per galaxy
  ! nGals        => Number of galaxies
  ! ages_Myr    <=  Mass-weighted average stellar ages
  ! Zs          <=  Mass-weighted average metallicities
  !-----------------------------------------------------------------------
    implicit none 
    character(1000),intent(in) :: starfile
    integer,intent(in) :: nGals
    real(kind=4),intent(out)   :: ages_Myr(nGals), Zs(nGals)
    integer::iGal, galNr, nGals2, iStar, nStars
    real(kind=8)::mTot
    real::sAge, sMass, sZ, sPosx, sPosy, sPosz, sMassInit
  !-----------------------------------------------------------------------
    ! Open file generated by halofinder containing stars in each galaxy
    open(unit=11,file=starfile,status='old',form='formatted')
    read(11, '(I12)') nGals2                        !   Number of galaxies
    if(nGals2 .ne. nGals) then
       print*,'nGals!=nGals2!!!'
       print*,'nGals = ',nGals
       print*,'nGals2 = ',nGals2
       return
    endif

    ages_Myr(:) = 0.
    Zs(:) = 0.
    do iGal=1,nGals
       write (*, "(A, I5, A, I5, A)", advance='no') &
            ' HaloFinder/f90/py_halo_utils.f90: Processing galaxy '     &
            ,iGal,' / ',nGals, char(13)
       read(11, '(I12, I12)') galNr, nStars
       mTot=0d0
       do iStar=1,nStars
          read(11, '(7(e14.6,1x))') &
               sAge, sMass, sZ, sPosx, sPosy, sPosz, sMassInit
          if(sAge .le. 0d0) then
             print*,'Negative age in HalFinder/f90/py_halo_utils.f90->get_ages_metallicities() = ',sAge
             sAge = 1d-8
          endif
          mTot = mTot + sMass
          ages_Myr(iGal) = ages_Myr(iGal) + sAge * sMass
          Zs(iGal)   = Zs(iGal)   + sZ   * sMass 
       end do ! Loop over stars
       ages_Myr(iGal) = ages_Myr(iGal) / mtot / 1d6  ! age in Myr
       Zs(iGal)   = Zs(iGal)   / mTot
    end do ! Loop over galaxies
    write(*,*) char(13)
    close(11)

  end subroutine get_ages_metallicities

end module py_halo_utils

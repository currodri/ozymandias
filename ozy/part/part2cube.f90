subroutine integratecube(repository,ageweight,periodic,star,xmin,xmax,ymin,ymax,zmin,zmax,&
      & mtot,nparttot)
  !--------------------------------------------------------------------------
  ! Ce programme calcule la carte de densite surfacique projetee
  ! des particules de matiere noire d'une simulation RAMSES. 
  ! Version F90 par R. Teyssier le 01/04/01.
  !--------------------------------------------------------------------------
  implicit none
  
  real(KIND=8) :: xmin,xmax,ymin,ymax,zmin,zmax
  character(LEN=128) :: repository
  logical::periodic,star,ageweight
  ! Output variables
  integer :: nparttot
  real(KIND=8)::mtot

!f2py intent(in) repository,ageweight,periodic,star,xmin,xmax,ymin,ymax,zmin,zmax
!f2py intent(out) mtot
!f2py intent(out) nparttot

  integer::ncpu,ndim,npart,ngrid,n,i,j,k,icpu,ipos,n_frw,nstar,npartcount
  integer::ncpu2,npart2,ndim2,levelmin,levelmax,ilevel,iii
  integer::nx=0,ny=0,nz=0,ix,iy,iz,ixp1,iyp1,izp1,idim,jdim,kdim,ncpu_read
  real(KIND=8)::ddx,ddy,ddz,dex,dey,dez,t,time,time_tot,time_simu,weight
  real(KIND=8)::aexp,omega_m,omega_l,omega_b,omega_k,h0,unit_l,unit_t,unit_d
  integer::imin,imax,jmin,jmax,kmin,kmax,lmin
  real(KIND=8)::xxmin,xxmax,yymin,yymax,zzmin,zzmax,dx,dy,dz,deltax
  real(KIND=8),dimension(:),allocatable::aexp_frw,hexp_frw,tau_frw,t_frw
  real(KIND=8),dimension(:,:),allocatable::x
  real(KIND=8),dimension(:)  ,allocatable::m,age
  character(LEN=1)::proj='z'
  character(LEN=5)::nchar,ncharcpu
  character(LEN=80)::ordering
  character(LEN=80)::GMGM
  character(LEN=128)::nomfich
  logical::ok,ok_part
  integer::impi,ndom,bit_length,maxdom
  integer,dimension(1:8)::idom,jdom,kdom,cpu_min,cpu_max
  real(KIND=8),dimension(1:8)::bounding_min,bounding_max
  real(KIND=8)::dkey,order_min,dmax
  real(kind=8),dimension(:),allocatable::bound_key
  logical,dimension(:),allocatable::cpu_read
  integer,dimension(:),allocatable::cpu_list

  !-----------------------------------------------
  ! Lecture du fichier particules au format RAMSES
  !-----------------------------------------------
  ipos=INDEX(repository,'output_')
  nchar=repository(ipos+7:ipos+13)
  nomfich=TRIM(repository)//'/part_'//TRIM(nchar)//'.out00001'
  inquire(file=nomfich, exist=ok) ! verify input file 
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found.'
     stop
  endif

  nomfich=TRIM(repository)//'/info_'//TRIM(nchar)//'.txt'
  inquire(file=nomfich, exist=ok) ! verify input file 
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found.'
     stop
  endif
  open(unit=10,file=nomfich,form='formatted',status='old')
  ! read(10,'("ncpu        =",I11)')ncpu
  read(10,'(A13,I11)')GMGM,ncpu
  ! read(10,'("ndim        =",I11)')ndim
  read(10,'(A13,I11)')GMGM,ndim
  ! read(10,'("levelmin    =",I11)')levelmin
  read(10,'(A13,I11)')GMGM,levelmin
  ! read(10,'("levelmax    =",I11)')levelmax
  read(10,'(A13,I11)')GMGM,levelmax
  read(10,*)
  read(10,*)
  read(10,*)

  read(10,*)
  ! read(10,'("time        =",E23.15)')t
  read(10,'(A13,E23.15)')GMGM,t

!  read(10,'("aexp        =",E23.15)')aexp
  read(10,'(A13,E23.15)')GMGM,aexp

!  read(10,'("H0          =",E23.15)')h0
  read(10,'(A13,E23.15)')GMGM,h0

  !read(10,'("omega_m     =",E23.15)')omega_m
  read(10,'(A13,E23.15)')GMGM,omega_m

  !read(10,'("omega_l     =",E23.15)')omega_l
  read(10,'(A13,E23.15)')GMGM,omega_l

  !read(10,'("omega_k     =",E23.15)')omega_k
  read(10,'(A13,E23.15)')GMGM,omega_k

  !read(10,'("omega_b     =",E23.15)')omega_b
  read(10,'(A13,E23.15)')GMGM,omega_b

  !read(10,'("unit_l      =",E23.15)')unit_l
  read(10,'(A13,E23.15)')GMGM,unit_l

  !read(10,'("unit_d      =",E23.15)')unit_d
  read(10,'(A13,E23.15)')GMGM,unit_d

  !read(10,'("unit_t      =",E23.15)')unit_t
  read(10,'(A13,E23.15)')GMGM,unit_t

  read(10,*)

  ! read(10,'("ordering type=",A80)'),ordering
  read(10,'(A14,A80)')GMGM,ordering
  write(*,'(" ordering type=",A20)'),TRIM(ordering)
  read(10,*)
  allocate(cpu_list(1:ncpu))
  if(TRIM(ordering).eq.'hilbert')then
     allocate(bound_key(0:ncpu))
     allocate(cpu_read(1:ncpu))
     cpu_read=.false.
     do impi=1,ncpu
        read(10,'(I8,1X,E23.15,1X,E23.15)')i,bound_key(impi-1),bound_key(impi)
     end do
  endif
  close(10)


  !-----------------------
  ! Cosmological model
  !-----------------------
  ! Allocate look-up tables
  n_frw=1000
  allocate(aexp_frw(0:n_frw),hexp_frw(0:n_frw))
  allocate(tau_frw(0:n_frw),t_frw(0:n_frw))
  
  ! Compute Friedman model look up table
  write(*,*)'Computing Friedman model'
  call friedman(dble(omega_m),dble(omega_l),dble(omega_k), &
       & 1.d-6,1.d-3,aexp_frw,hexp_frw,tau_frw,t_frw,n_frw,time_tot)
  
  ! Find neighboring expansion factors
  i=1
  do while(aexp_frw(i)>aexp.and.i<n_frw)
     i=i+1
  end do
  ! Interploate time
  time_simu=t_frw(i)*(aexp-aexp_frw(i-1))/(aexp_frw(i)-aexp_frw(i-1))+ &
       & t_frw(i-1)*(aexp-aexp_frw(i))/(aexp_frw(i-1)-aexp_frw(i))
  write(*,*)'Age simu=',(time_tot+time_simu)/(h0*1d5/3.08d24)/(365.*24.*3600.*1d9)

  !-----------------------
  ! Map parameters
  !-----------------------
  if(nx==0)then
     nx=2**levelmin
  endif
  if(ny==0)then
     ny=nx
  end if
  if(nz==0)then
     nz=nx
  end if
  write(*,*)'time=',t
  idim=1
  jdim=2
  kdim=3
  xxmin=xmin ; xxmax=xmax
  yymin=ymin ; yymax=ymax
  zzmin=zmin ; zzmax=zmax
  dx=(xxmax-xxmin)/dble(nx)
  dy=(yymax-yymin)/dble(ny)
  dz=(zzmax-zzmin)/dble(nz)

  if(TRIM(ordering).eq.'hilbert')then

     dmax=max(xmax-xmin,ymax-ymin,zmax-zmin)
     do ilevel=1,levelmax
        deltax=0.5d0**ilevel
        if(deltax.lt.dmax)exit
     end do
     lmin=ilevel
     bit_length=lmin-1
     maxdom=2**bit_length
     imin=0; imax=0; jmin=0; jmax=0; kmin=0; kmax=0
     if(bit_length>0)then
        imin=int(xmin*dble(maxdom))
        imax=imin+1
        jmin=int(ymin*dble(maxdom))
        jmax=jmin+1
        kmin=int(zmin*dble(maxdom))
        kmax=kmin+1
     endif
     
     dkey=(dble(2**(levelmax+1)/dble(maxdom)))**ndim
     ndom=1
     if(bit_length>0)ndom=8
     idom(1)=imin; idom(2)=imax
     idom(3)=imin; idom(4)=imax
     idom(5)=imin; idom(6)=imax
     idom(7)=imin; idom(8)=imax
     jdom(1)=jmin; jdom(2)=jmin
     jdom(3)=jmax; jdom(4)=jmax
     jdom(5)=jmin; jdom(6)=jmin
     jdom(7)=jmax; jdom(8)=jmax
     kdom(1)=kmin; kdom(2)=kmin
     kdom(3)=kmin; kdom(4)=kmin
     kdom(5)=kmax; kdom(6)=kmax
     kdom(7)=kmax; kdom(8)=kmax
     
     do i=1,ndom
        if(bit_length>0)then
           call hilbert3d(idom(i),jdom(i),kdom(i),order_min,bit_length,1)
        else
           order_min=0.0d0
        endif
        bounding_min(i)=(order_min)*dkey
        bounding_max(i)=(order_min+1.0D0)*dkey
     end do
     cpu_min=0; cpu_max=0
     do impi=1,ncpu
        do i=1,ndom
           if (   bound_key(impi-1).le.bounding_min(i).and.&
                & bound_key(impi  ).gt.bounding_min(i))then
              cpu_min(i)=impi
           endif
           if (   bound_key(impi-1).lt.bounding_max(i).and.&
                & bound_key(impi  ).ge.bounding_max(i))then
              cpu_max(i)=impi
           endif
        end do
     end do
     
     ncpu_read=0
     do i=1,ndom
        do j=cpu_min(i),cpu_max(i)
           if(.not. cpu_read(j))then
              ncpu_read=ncpu_read+1
              cpu_list(ncpu_read)=j
              cpu_read(j)=.true.
           endif
        enddo
     enddo
  else
     ncpu_read=ncpu
     do j=1,ncpu
        cpu_list(j)=j
     end do
  end  if

  npart=0
  do k=1,ncpu_read
     icpu=cpu_list(k)
     call title(icpu,ncharcpu)
     nomfich=TRIM(repository)//'/part_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
     open(unit=1,file=nomfich,status='old',form='unformatted')
     read(1)ncpu2
     read(1)ndim2
     read(1)npart2
     read(1)
     read(1)nstar
     close(1)
     npart=npart+npart2
  end do
  write(*,*)'Found ',npart,' particles'
  if(nstar>0)then
     if(star)then
        write(*,*)'Keeping star particles.'
     else
        write(*,*)'Discard star particles.'
     endif
  endif

  !-----------------------------------------------
  ! Compute projected mass using CIC smoothing
  !----------------------------------------------
  mtot=0.0d0
  nparttot=0
  do k=1,ncpu_read
     icpu=cpu_list(k)
     call title(icpu,ncharcpu)
     nomfich=TRIM(repository)//'/part_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
     open(unit=1,file=nomfich,status='old',form='unformatted')
   !   write(*,*)'Processing file '//TRIM(nomfich)
     read(1)ncpu2
     read(1)ndim2
     read(1)npart2
     read(1)
     read(1)
     read(1)
     read(1)
     read(1)
     allocate(m(1:npart2))
     if(nstar>0)allocate(age(1:npart2))
     allocate(x(1:npart2,1:ndim2))
     ! Read position
     do i=1,ndim
        read(1)m
        x(1:npart2,i)=m
     end do
     ! Skip velocity
     do i=1,ndim
        read(1)m
     end do
     ! Read mass
     read(1)m
     if(nstar>0)then
        read(1) ! Skip identity
        read(1) ! Skip level
        read(1)age
     endif
     close(1)
     if(periodic)then
     do i=1,npart2
        ok_part=(x(i,1)>=xmin.and.x(i,1)<=xmax.and. &
             &   x(i,2)>=ymin.and.x(i,2)<=ymax.and. &
             &   x(i,3)>=zmin.and.x(i,3)<=zmax)

        if(nstar>0)then
           if(star)then
              ok_part=ok_part.and.(age(i).ne.0.0d0)
              if(ageweight)then
                 iii=1
                 do while(tau_frw(iii)>age(i).and.iii<n_frw)
                    iii=iii+1
                 end do
                 ! Interploate time
                 time=t_frw(iii)*(age(i)-tau_frw(iii-1))/(tau_frw(iii)-tau_frw(iii-1))+ &
                      & t_frw(iii-1)*(age(i)-tau_frw(iii))/(tau_frw(iii-1)-tau_frw(iii))
                 time=(time_simu-time)/(h0*1d5/3.08d24)/(365.*24.*3600.*1d9)
                 if(time>0.01)then
                    weight=(time/0.01)**(-0.7)
                 endif
              endif
           else
              ok_part=ok_part.and.(age(i).eq.0.0d0)
           endif
        endif

        if(ok_part)then
           ddx=(x(i,idim)-xxmin)/dx
           ddy=(x(i,jdim)-yymin)/dy
           ddz=(x(i,kdim)-zzmin)/dz
           ix=ddx
           iy=ddy
           iz=ddz
           ddx=ddx-ix
           ddy=ddy-iy
           ddz=ddz-iz
           dex=1.0-ddx
           dey=1.0-ddy
           dez=1.0-ddz
           if(ix<0)ix=ix+nx
           if(ix>=nx)ix=ix-nx
           if(iy<0)iy=iy+ny
           if(iy>=ny)iy=iy-ny
           if(iz<0)iz=iz+nz
           if(iz>=nz)iz=iz-nz
           ixp1=ix+1
           iyp1=iy+1
           izp1=iz+1
           if(ixp1<0)ixp1=ixp1+nx
           if(ixp1>=nx)ixp1=ixp1-nx
           if(iyp1<0)iyp1=iyp1+ny
           if(iyp1>=ny)iyp1=iyp1-ny
           if(izp1<0)izp1=izp1+nz
           if(izp1>=nz)izp1=izp1-nz
           mtot=mtot+m(i)
           nparttot=nparttot+1
        end if
     end do
     else
     do i=1,npart2
        weight=1.0
        ok_part=(x(i,1)>=xmin.and.x(i,1)<=xmax.and. &
             &   x(i,2)>=ymin.and.x(i,2)<=ymax.and. &
             &   x(i,3)>=zmin.and.x(i,3)<=zmax)

        if(nstar>0)then
           if(star)then
              ok_part=ok_part.and.(age(i).ne.0.0d0)
              if(ageweight)then
                 iii=1
                 do while(tau_frw(iii)>age(i).and.iii<n_frw)
                    iii=iii+1
                 end do
                 ! Interpolate time
                 time=t_frw(iii)*(age(i)-tau_frw(iii-1))/(tau_frw(iii)-tau_frw(iii-1))+ &
                      & t_frw(iii-1)*(age(i)-tau_frw(iii))/(tau_frw(iii-1)-tau_frw(iii))
                 time=(time_simu-time)/(h0*1d5/3.08d24)/(365.*24.*3600.*1d9)
                 if(time>0.01)then
                    weight=(time/0.01)**(-0.7)
                 endif
              endif
           else
              ok_part=ok_part.and.(age(i).eq.0.0d0)
           endif
        endif

        if(ok_part)then
           ddx=(x(i,idim)-xxmin)/dx
           ddy=(x(i,jdim)-yymin)/dy
           ddz=(x(i,kdim)-zzmin)/dz
           ix=ddx
           iy=ddy
           iz=ddz
           ddx=ddx-ix
           ddy=ddy-iy
           ddz=ddz-iz
           dex=1.0-ddx
           dey=1.0-ddy
           dez=1.0-ddz
           ixp1=ix+1
           iyp1=iy+1
           izp1=iz+1
           if(ix>=0.and.ix<nx.and.iy>=0.and.iy<ny.and.iz>=0.and.iz<nz)then
              mtot=mtot+m(i)
              nparttot=nparttot+1
           end if
        end if
     end do
     endif
     deallocate(x,m)
     if(nstar>0)deallocate(age)
  end do
!   write(*,*)'Total mass=',mtot
!   write(*,*)'n Particles selected=',nparttot
  return
end subroutine integratecube
!================================================================
!================================================================
!================================================================
!================================================================
subroutine integratesphere(repository,ageweight,periodic,star,xcenter,ycenter,zcenter,radius,&
      & mtot,nparttot)
  !--------------------------------------------------------------------------
  ! Ce programme calcule la carte de densite surfacique projetee
  ! des particules de matiere noire d'une simulation RAMSES. 
  ! Version F90 par R. Teyssier le 01/04/01.
  !--------------------------------------------------------------------------
  implicit none
  
  real(KIND=8) :: xcenter,ycenter,zcenter,radius
  character(LEN=128) :: repository
  logical::periodic,star,ageweight
  ! Output variables
  integer :: nparttot
  real(KIND=8)::mtot

!f2py intent(in) repository,ageweight,periodic,star,xcenter,ycenter,zcenter,radius
!f2py intent(out) mtot
!f2py intent(out) nparttot

  integer::ncpu,ndim,npart,ngrid,n,i,j,k,icpu,ipos,n_frw,nstar,npartcount
  integer::ncpu2,npart2,ndim2,levelmin,levelmax,ilevel,iii
  integer::nx=0,ny=0,nz=0,ix,iy,iz,ixp1,iyp1,izp1,idim,jdim,kdim,ncpu_read
  real(KIND=8)::ddx,ddy,ddz,dex,dey,dez,t,time,time_tot,time_simu,weight
  real(KIND=8)::aexp,omega_m,omega_l,omega_b,omega_k,h0,unit_l,unit_t,unit_d
  integer::imin,imax,jmin,jmax,kmin,kmax,lmin
  real(KIND=8)::xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1
  real(KIND=8)::xxmin,xxmax,yymin,yymax,zzmin,zzmax,dx,dy,dz,deltax
  real(KIND=8),dimension(:),allocatable::aexp_frw,hexp_frw,tau_frw,t_frw
  real(KIND=8),dimension(:,:),allocatable::x
  real(KIND=8),dimension(:),allocatable::m,age,r
  real(KIND=8),dimension(1:3) :: center
  character(LEN=1)::proj='z'
  character(LEN=5)::nchar,ncharcpu
  character(LEN=80)::ordering
  character(LEN=80)::GMGM
  character(LEN=128)::nomfich
  logical::ok,ok_part
  integer::impi,ndom,bit_length,maxdom
  integer,dimension(1:8)::idom,jdom,kdom,cpu_min,cpu_max
  real(KIND=8),dimension(1:8)::bounding_min,bounding_max
  real(KIND=8)::dkey,order_min,dmax
  real(kind=8),dimension(:),allocatable::bound_key
  logical,dimension(:),allocatable::cpu_read
  integer,dimension(:),allocatable::cpu_list

  !-----------------------------------------------
  ! Lecture du fichier particules au format RAMSES
  !-----------------------------------------------
  ipos=INDEX(repository,'output_')
  nchar=repository(ipos+7:ipos+13)
  nomfich=TRIM(repository)//'/part_'//TRIM(nchar)//'.out00001'
  inquire(file=nomfich, exist=ok) ! verify input file 
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found.'
     stop
  endif

  nomfich=TRIM(repository)//'/info_'//TRIM(nchar)//'.txt'
  inquire(file=nomfich, exist=ok) ! verify input file 
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found.'
     stop
  endif
  open(unit=10,file=nomfich,form='formatted',status='old')
  ! read(10,'("ncpu        =",I11)')ncpu
  read(10,'(A13,I11)')GMGM,ncpu
  ! read(10,'("ndim        =",I11)')ndim
  read(10,'(A13,I11)')GMGM,ndim
  ! read(10,'("levelmin    =",I11)')levelmin
  read(10,'(A13,I11)')GMGM,levelmin
  ! read(10,'("levelmax    =",I11)')levelmax
  read(10,'(A13,I11)')GMGM,levelmax
  read(10,*)
  read(10,*)
  read(10,*)

  read(10,*)
  ! read(10,'("time        =",E23.15)')t
  read(10,'(A13,E23.15)')GMGM,t

!  read(10,'("aexp        =",E23.15)')aexp
  read(10,'(A13,E23.15)')GMGM,aexp

!  read(10,'("H0          =",E23.15)')h0
  read(10,'(A13,E23.15)')GMGM,h0

  !read(10,'("omega_m     =",E23.15)')omega_m
  read(10,'(A13,E23.15)')GMGM,omega_m

  !read(10,'("omega_l     =",E23.15)')omega_l
  read(10,'(A13,E23.15)')GMGM,omega_l

  !read(10,'("omega_k     =",E23.15)')omega_k
  read(10,'(A13,E23.15)')GMGM,omega_k

  !read(10,'("omega_b     =",E23.15)')omega_b
  read(10,'(A13,E23.15)')GMGM,omega_b

  !read(10,'("unit_l      =",E23.15)')unit_l
  read(10,'(A13,E23.15)')GMGM,unit_l

  !read(10,'("unit_d      =",E23.15)')unit_d
  read(10,'(A13,E23.15)')GMGM,unit_d

  !read(10,'("unit_t      =",E23.15)')unit_t
  read(10,'(A13,E23.15)')GMGM,unit_t

  read(10,*)

  ! read(10,'("ordering type=",A80)'),ordering
  read(10,'(A14,A80)')GMGM,ordering
  write(*,'(" ordering type=",A20)'),TRIM(ordering)
  read(10,*)
  allocate(cpu_list(1:ncpu))
  if(TRIM(ordering).eq.'hilbert')then
     allocate(bound_key(0:ncpu))
     allocate(cpu_read(1:ncpu))
     cpu_read=.false.
     do impi=1,ncpu
        read(10,'(I8,1X,E23.15,1X,E23.15)')i,bound_key(impi-1),bound_key(impi)
     end do
  endif
  close(10)


  !-----------------------
  ! Cosmological model
  !-----------------------
  ! Allocate look-up tables
  n_frw=1000
  allocate(aexp_frw(0:n_frw),hexp_frw(0:n_frw))
  allocate(tau_frw(0:n_frw),t_frw(0:n_frw))
  
  ! Compute Friedman model look up table
  write(*,*)'Computing Friedman model'
  call friedman(dble(omega_m),dble(omega_l),dble(omega_k), &
       & 1.d-6,1.d-3,aexp_frw,hexp_frw,tau_frw,t_frw,n_frw,time_tot)
  
  ! Find neighboring expansion factors
  i=1
  do while(aexp_frw(i)>aexp.and.i<n_frw)
     i=i+1
  end do
  ! Interploate time
  time_simu=t_frw(i)*(aexp-aexp_frw(i-1))/(aexp_frw(i)-aexp_frw(i-1))+ &
       & t_frw(i-1)*(aexp-aexp_frw(i))/(aexp_frw(i-1)-aexp_frw(i))
  write(*,*)'Age simu=',(time_tot+time_simu)/(h0*1d5/3.08d24)/(365.*24.*3600.*1d9)

  !-----------------------
  ! Map parameters
  !-----------------------
  if(nx==0)then
     nx=2**levelmin
  endif
  if(ny==0)then
     ny=nx
  end if
  if(nz==0)then
     nz=nx
  end if
  write(*,*)'time=',t
  idim=1
  jdim=2
  kdim=3
  center=(/xcenter,ycenter,zcenter/)
  xmin=xcenter-radius ; xmax=xcenter+radius
  ymin=ycenter-radius ; ymax=ycenter+radius
  zmin=zcenter-radius ; zmax=zcenter+radius
  xxmin=xmin ; xxmax=xmax
  yymin=ymin ; yymax=ymax
  zzmin=zmin ; zzmax=zmax
  dx=(xxmax-xxmin)/dble(nx)
  dy=(yymax-yymin)/dble(ny)
  dz=(zzmax-zzmin)/dble(nz)

  if(TRIM(ordering).eq.'hilbert')then

     dmax=max(xmax-xmin,ymax-ymin,zmax-zmin)
     do ilevel=1,levelmax
        deltax=0.5d0**ilevel
        if(deltax.lt.dmax)exit
     end do
     lmin=ilevel
     bit_length=lmin-1
     maxdom=2**bit_length
     imin=0; imax=0; jmin=0; jmax=0; kmin=0; kmax=0
     if(bit_length>0)then
        imin=int(xmin*dble(maxdom))
        imax=imin+1
        jmin=int(ymin*dble(maxdom))
        jmax=jmin+1
        kmin=int(zmin*dble(maxdom))
        kmax=kmin+1
     endif
     
     dkey=(dble(2**(levelmax+1)/dble(maxdom)))**ndim
     ndom=1
     if(bit_length>0)ndom=8
     idom(1)=imin; idom(2)=imax
     idom(3)=imin; idom(4)=imax
     idom(5)=imin; idom(6)=imax
     idom(7)=imin; idom(8)=imax
     jdom(1)=jmin; jdom(2)=jmin
     jdom(3)=jmax; jdom(4)=jmax
     jdom(5)=jmin; jdom(6)=jmin
     jdom(7)=jmax; jdom(8)=jmax
     kdom(1)=kmin; kdom(2)=kmin
     kdom(3)=kmin; kdom(4)=kmin
     kdom(5)=kmax; kdom(6)=kmax
     kdom(7)=kmax; kdom(8)=kmax
     
     do i=1,ndom
        if(bit_length>0)then
           call hilbert3d(idom(i),jdom(i),kdom(i),order_min,bit_length,1)
        else
           order_min=0.0d0
        endif
        bounding_min(i)=(order_min)*dkey
        bounding_max(i)=(order_min+1.0D0)*dkey
     end do
     cpu_min=0; cpu_max=0
     do impi=1,ncpu
        do i=1,ndom
           if (   bound_key(impi-1).le.bounding_min(i).and.&
                & bound_key(impi  ).gt.bounding_min(i))then
              cpu_min(i)=impi
           endif
           if (   bound_key(impi-1).lt.bounding_max(i).and.&
                & bound_key(impi  ).ge.bounding_max(i))then
              cpu_max(i)=impi
           endif
        end do
     end do
     
     ncpu_read=0
     do i=1,ndom
        do j=cpu_min(i),cpu_max(i)
           if(.not. cpu_read(j))then
              ncpu_read=ncpu_read+1
              cpu_list(ncpu_read)=j
              cpu_read(j)=.true.
           endif
        enddo
     enddo
  else
     ncpu_read=ncpu
     do j=1,ncpu
        cpu_list(j)=j
     end do
  end  if

  npart=0
  do k=1,ncpu_read
     icpu=cpu_list(k)
     call title(icpu,ncharcpu)
     nomfich=TRIM(repository)//'/part_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
     open(unit=1,file=nomfich,status='old',form='unformatted')
     read(1)ncpu2
     read(1)ndim2
     read(1)npart2
     read(1)
     read(1)nstar
     close(1)
     npart=npart+npart2
  end do
  write(*,*)'Found ',npart,' particles'
  if(nstar>0)then
     if(star)then
        write(*,*)'Keeping star particles.'
     else
        write(*,*)'Discard star particles.'
     endif
  endif

  !-----------------------------------------------
  ! Compute projected mass using CIC smoothing
  !----------------------------------------------
  mtot=0.0d0
  nparttot=0
  do k=1,ncpu_read
     icpu=cpu_list(k)
     call title(icpu,ncharcpu)
     nomfich=TRIM(repository)//'/part_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
     open(unit=1,file=nomfich,status='old',form='unformatted')
   !   write(*,*)'Processing file '//TRIM(nomfich)
     read(1)ncpu2
     read(1)ndim2
     read(1)npart2
     read(1)
     read(1)
     read(1)
     read(1)
     read(1)
     allocate(m(1:npart2))
     allocate(r(1:npart2))
     if(nstar>0)allocate(age(1:npart2))
     allocate(x(1:npart2,1:ndim2))
     ! Read position
     do i=1,ndim
        read(1)m
        x(1:npart2,i)=m
     end do
     ! Compute radial distance
     call radial_distance(ndim,center,x,npart2,r)
     deallocate(x)
     ! Skip velocity
     do i=1,ndim
        read(1)m
     end do
     ! Read mass
     read(1)m
     if(nstar>0)then
        read(1) ! Skip identity
        read(1) ! Skip level
        read(1)age
     endif
     close(1)
     if(periodic)then
     do i=1,npart2
        ok_part=(r(i) <= radius)
        if(nstar>0)then
           if(star)then
              ok_part=ok_part.and.(age(i).ne.0.0d0)
              if(ageweight)then
                 iii=1
                 do while(tau_frw(iii)>age(i).and.iii<n_frw)
                    iii=iii+1
                 end do
                 ! Interploate time
                 time=t_frw(iii)*(age(i)-tau_frw(iii-1))/(tau_frw(iii)-tau_frw(iii-1))+ &
                      & t_frw(iii-1)*(age(i)-tau_frw(iii))/(tau_frw(iii-1)-tau_frw(iii))
                 time=(time_simu-time)/(h0*1d5/3.08d24)/(365.*24.*3600.*1d9)
                 if(time>0.01)then
                    weight=(time/0.01)**(-0.7)
                 endif
              endif
           else
              ok_part=ok_part.and.(age(i).eq.0.0d0)
           endif
        endif

        if(ok_part)then
           mtot=mtot+m(i)
           nparttot=nparttot+1
        end if
     end do
     else
     do i=1,npart2
        weight=1.0
        ok_part=(r(i) <= radius)

        if(nstar>0)then
           if(star)then
              ok_part=ok_part.and.(age(i).ne.0.0d0)
              if(ageweight)then
                 iii=1
                 do while(tau_frw(iii)>age(i).and.iii<n_frw)
                    iii=iii+1
                 end do
                 ! Interpolate time
                 time=t_frw(iii)*(age(i)-tau_frw(iii-1))/(tau_frw(iii)-tau_frw(iii-1))+ &
                      & t_frw(iii-1)*(age(i)-tau_frw(iii))/(tau_frw(iii-1)-tau_frw(iii))
                 time=(time_simu-time)/(h0*1d5/3.08d24)/(365.*24.*3600.*1d9)
                 if(time>0.01)then
                    weight=(time/0.01)**(-0.7)
                 endif
              endif
           else
              ok_part=ok_part.and.(age(i).eq.0.0d0)
           endif
        endif

        if(ok_part)then
            mtot=mtot+m(i)
            nparttot=nparttot+1
        end if
     end do
     endif
     deallocate(r,m)
     if(nstar>0)deallocate(age)
  end do
!   write(*,*)'Total mass=',mtot
!   write(*,*)'n Particles selected=',nparttot
  return
end subroutine integratesphere
!=======================================================================
subroutine title(n,nchar)
!=======================================================================
  implicit none
  integer::n
  character*5::nchar

  character*1::nchar1
  character*2::nchar2
  character*3::nchar3
  character*4::nchar4
  character*5::nchar5

  if(n.ge.10000)then
     write(nchar5,'(i5)') n
     nchar = nchar5
  elseif(n.ge.1000)then
     write(nchar4,'(i4)') n
     nchar = '0'//nchar4
  elseif(n.ge.100)then
     write(nchar3,'(i3)') n
     nchar = '00'//nchar3
  elseif(n.ge.10)then
     write(nchar2,'(i2)') n
     nchar = '000'//nchar2
  else
     write(nchar1,'(i1)') n
     nchar = '0000'//nchar1
  endif

end subroutine title
!================================================================
!================================================================
!================================================================
!================================================================
subroutine radial_distance(ndim,center,x,npart, r)
   implicit none
   integer :: ndim,j,npart
   real(KIND=8),dimension(1:ndim) :: center
   real(KIND=8),dimension(1:npart,1:ndim) :: x
   real(KIND=8),dimension(1:npart) :: r
   do j=1,npart
      r(j)=sqrt(sum((x(j,:)-center)**2))
   end do
end subroutine radial_distance
!================================================================
!================================================================
!================================================================
!================================================================
subroutine hilbert3d(x,y,z,order,bit_length,npoint)
  implicit none
  integer     ,INTENT(IN)                     ::bit_length,npoint
  integer     ,INTENT(IN) ,dimension(1:npoint)::x,y,z
  real(kind=8),INTENT(OUT),dimension(1:npoint)::order

  logical,dimension(0:3*bit_length-1)::i_bit_mask
  logical,dimension(0:1*bit_length-1)::x_bit_mask,y_bit_mask,z_bit_mask
  integer,dimension(0:7,0:1,0:11)::state_diagram
  integer::i,ip,cstate,nstate,b0,b1,b2,sdigit,hdigit

  if(bit_length>bit_size(bit_length))then
     write(*,*)'Maximum bit length=',bit_size(bit_length)
     write(*,*)'stop in hilbert3d'
     stop
  endif

  state_diagram = RESHAPE( (/   1, 2, 3, 2, 4, 5, 3, 5,&
                            &   0, 1, 3, 2, 7, 6, 4, 5,&
                            &   2, 6, 0, 7, 8, 8, 0, 7,&
                            &   0, 7, 1, 6, 3, 4, 2, 5,&
                            &   0, 9,10, 9, 1, 1,11,11,&
                            &   0, 3, 7, 4, 1, 2, 6, 5,&
                            &   6, 0, 6,11, 9, 0, 9, 8,&
                            &   2, 3, 1, 0, 5, 4, 6, 7,&
                            &  11,11, 0, 7, 5, 9, 0, 7,&
                            &   4, 3, 5, 2, 7, 0, 6, 1,&
                            &   4, 4, 8, 8, 0, 6,10, 6,&
                            &   6, 5, 1, 2, 7, 4, 0, 3,&
                            &   5, 7, 5, 3, 1, 1,11,11,&
                            &   4, 7, 3, 0, 5, 6, 2, 1,&
                            &   6, 1, 6,10, 9, 4, 9,10,&
                            &   6, 7, 5, 4, 1, 0, 2, 3,&
                            &  10, 3, 1, 1,10, 3, 5, 9,&
                            &   2, 5, 3, 4, 1, 6, 0, 7,&
                            &   4, 4, 8, 8, 2, 7, 2, 3,&
                            &   2, 1, 5, 6, 3, 0, 4, 7,&
                            &   7, 2,11, 2, 7, 5, 8, 5,&
                            &   4, 5, 7, 6, 3, 2, 0, 1,&
                            &  10, 3, 2, 6,10, 3, 4, 4,&
                            &   6, 1, 7, 0, 5, 2, 4, 3 /), &
                            & (/8 ,2, 12 /) )

  do ip=1,npoint

     ! convert to binary
     do i=0,bit_length-1
        x_bit_mask(i)=btest(x(ip),i)
        y_bit_mask(i)=btest(y(ip),i)
        z_bit_mask(i)=btest(z(ip),i)
     enddo

     ! interleave bits
     do i=0,bit_length-1
        i_bit_mask(3*i+2)=x_bit_mask(i)
        i_bit_mask(3*i+1)=y_bit_mask(i)
        i_bit_mask(3*i  )=z_bit_mask(i)
     end do

     ! build Hilbert ordering using state diagram
     cstate=0
     do i=bit_length-1,0,-1
        b2=0 ; if(i_bit_mask(3*i+2))b2=1
        b1=0 ; if(i_bit_mask(3*i+1))b1=1
        b0=0 ; if(i_bit_mask(3*i  ))b0=1
        sdigit=b2*4+b1*2+b0
        nstate=state_diagram(sdigit,0,cstate)
        hdigit=state_diagram(sdigit,1,cstate)
        i_bit_mask(3*i+2)=btest(hdigit,2)
        i_bit_mask(3*i+1)=btest(hdigit,1)
        i_bit_mask(3*i  )=btest(hdigit,0)
        cstate=nstate
     enddo

     ! save Hilbert key as double precision real
     order(ip)=0.
     do i=0,3*bit_length-1
        b0=0 ; if(i_bit_mask(i))b0=1
        order(ip)=order(ip)+dble(b0)*dble(2)**i
     end do

  end do

end subroutine hilbert3d
!================================================================
!================================================================
!================================================================
!================================================================
!================================================================
!================================================================
!================================================================
!================================================================
subroutine friedman(O_mat_0,O_vac_0,O_k_0,alpha,axp_min, &
     & axp_out,hexp_out,tau_out,t_out,ntable,age_tot)

  implicit none
  integer::ntable
  real(kind=8)::O_mat_0, O_vac_0, O_k_0
  real(kind=8)::alpha,axp_min,age_tot
  real(kind=8),dimension(0:ntable)::axp_out,hexp_out,tau_out,t_out
  ! ######################################################!
  ! This subroutine assumes that axp = 1 at z = 0 (today) !
  ! and that t and tau = 0 at z = 0 (today).              !
  ! axp is the expansion factor, hexp the Hubble constant !
  ! defined as hexp=1/axp*daxp/dtau, tau the conformal    !
  ! time, and t the look-back time, both in unit of 1/H0. !
  ! alpha is the required accuracy and axp_min is the     !
  ! starting expansion factor of the look-up table.       !
  ! ntable is the required size of the look-up table.     !
  ! ######################################################!
  real(kind=8)::axp_tau, axp_t
  real(kind=8)::axp_tau_pre, axp_t_pre
  real(kind=8)::dadtau, dadt
  real(kind=8)::dtau,dt
  real(kind=8)::tau,t
  integer::nstep,nout,nskip

!!$  if( (O_mat_0+O_vac_0+O_k_0) .ne. 1.0D0 )then
!!$     write(*,*)'Error: non-physical cosmological constants'
!!$     write(*,*)'O_mat_0,O_vac_0,O_k_0=',O_mat_0,O_vac_0,O_k_0
!!$     write(*,*)'The sum must be equal to 1.0, but '
!!$     write(*,*)'O_mat_0+O_vac_0+O_k_0=',O_mat_0+O_vac_0+O_k_0
!!$     stop
!!$  end if

  axp_tau = 1.0D0
  axp_t = 1.0D0
  tau = 0.0D0
  t = 0.0D0
  nstep = 0
  
  do while ( (axp_tau .ge. axp_min) .or. (axp_t .ge. axp_min) ) 
     
     nstep = nstep + 1
     dtau = alpha * axp_tau / dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)
     axp_tau_pre = axp_tau - dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)*dtau/2.d0
     axp_tau = axp_tau - dadtau(axp_tau_pre,O_mat_0,O_vac_0,O_k_0)*dtau
     tau = tau - dtau
     
     dt = alpha * axp_t / dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
     axp_t_pre = axp_t - dadt(axp_t,O_mat_0,O_vac_0,O_k_0)*dt/2.d0
     axp_t = axp_t - dadt(axp_t_pre,O_mat_0,O_vac_0,O_k_0)*dt
     t = t - dt
     
  end do

  age_tot=-t
  write(*,666)-t
  666 format(' Age of the Universe (in unit of 1/H0)=',1pe10.3)

  nskip=nstep/ntable
  
  axp_t = 1.d0
  t = 0.d0
  axp_tau = 1.d0
  tau = 0.d0
  nstep = 0
  nout=0
  t_out(nout)=t
  tau_out(nout)=tau
  axp_out(nout)=axp_tau
  hexp_out(nout)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau

  do while ( (axp_tau .ge. axp_min) .or. (axp_t .ge. axp_min) ) 
     
     nstep = nstep + 1
     dtau = alpha * axp_tau / dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)
     axp_tau_pre = axp_tau - dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)*dtau/2.d0
     axp_tau = axp_tau - dadtau(axp_tau_pre,O_mat_0,O_vac_0,O_k_0)*dtau
     tau = tau - dtau

     dt = alpha * axp_t / dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
     axp_t_pre = axp_t - dadt(axp_t,O_mat_0,O_vac_0,O_k_0)*dt/2.d0
     axp_t = axp_t - dadt(axp_t_pre,O_mat_0,O_vac_0,O_k_0)*dt
     t = t - dt
     
     if(mod(nstep,nskip)==0)then
        nout=nout+1
        t_out(nout)=t
        tau_out(nout)=tau
        axp_out(nout)=axp_tau
        hexp_out(nout)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau
     end if

  end do
  t_out(ntable)=t
  tau_out(ntable)=tau
  axp_out(ntable)=axp_tau
  hexp_out(ntable)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau

end subroutine friedman

function dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0) 
  real(kind=8)::dadtau,axp_tau,O_mat_0,O_vac_0,O_k_0
  dadtau = axp_tau*axp_tau*axp_tau *  &
       &   ( O_mat_0 + &
       &     O_vac_0 * axp_tau*axp_tau*axp_tau + &
       &     O_k_0   * axp_tau )
  dadtau = sqrt(dadtau)
  return
end function dadtau

function dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
  real(kind=8)::dadt,axp_t,O_mat_0,O_vac_0,O_k_0
  dadt   = (1.0D0/axp_t)* &
       &   ( O_mat_0 + &
       &     O_vac_0 * axp_t*axp_t*axp_t + &
       &     O_k_0   * axp_t )
  dadt = sqrt(dadt)
  return
end function dadt
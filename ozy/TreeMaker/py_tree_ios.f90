module py_tree_ios

  public

contains
  
  
  subroutine get_Nsteps(TreeDir,TreeNum,nsteps)
    implicit none
    character(1000),intent(in)  :: TreeDir
    integer(kind=4),intent(in)  :: TreeNum
    integer(kind=4),intent(out) :: nsteps
    character(1000) :: filename
    write(filename,'(a,a,i3.3,a)') trim(TreeDir),'/tstep_file_',TreeNum,'.001'  ! assume there is a single tree-file 
    open(unit=12,file=filename,status='old',form='unformatted')
    read(12) nsteps
    close(12)
    return
  end subroutine get_Nsteps


  subroutine get_nIDs(TreeDir,TreeNum,nIDs)
    implicit none
    character(1000),intent(in)  :: TreeDir
    integer(kind=4),intent(in)  :: TreeNum
    integer(kind=4),intent(out) :: nIDs
    character(1000) :: filename
    integer(kind=4) :: i,j
    write(filename,'(a,a,i3.3,a)') trim(TreeDir),'/tree_file_',TreeNum,'.001'  ! assume there is a single tree-file 
    open(unit=12,file=filename,status='old',form='unformatted')
    read(12) i,nIDs,j
    close(12)
    return
  end subroutine get_nIDs

  
  subroutine get_nProps(TreeDir,TreeNum,nProps)
    implicit none
    character(1000),intent(in)  :: TreeDir
    integer(kind=4),intent(in)  :: TreeNum
    integer(kind=4),intent(out) :: nProps
    character(1000) :: filename
    integer(kind=4) :: i
    write(filename,'(a,a,i3.3,a)') trim(TreeDir),'/props_',TreeNum,'.001'  ! assume there is a single tree-file 
    open(unit=12,file=filename,status='old',form='unformatted')
    read(12) i,nProps
    close(12)
    return
  end subroutine get_nProps

  
  subroutine read_timestep_props(TreeDir,TreeNum,nsteps,nhalos,aexp,age_univ)

    implicit none

    character(1000),intent(in)  :: TreeDir
    integer(kind=4),intent(in)  :: nsteps,TreeNum
    integer(kind=4),intent(out) :: nhalos(nsteps)
    real(kind=4),intent(out)    :: aexp(nsteps),age_univ(nsteps)
    character(1000) :: filename
    integer(kind=4) :: i

    write(filename,'(a,a,i3.3,a)') trim(TreeDir),'/tstep_file_',TreeNum,'.001'  ! assume there is a single tree-file 
    open(unit=12,file=filename,status='old',form='unformatted')
    read(12) ! skip nsteps
    read(12) nhalos
    read(12) aexp
    read(12) age_univ
    close(12)
    
    return
    
  end subroutine read_timestep_props
  
    
  subroutine read_tree_struct(TreeDir,TreeNum,nsteps,nhalos,nh,nIDs,IDs)
    
    implicit none

    ! nh is the sum of nhalos ... 
    
    character(1000),intent(in) :: TreeDir
    integer(kind=4),intent(in) :: nsteps,nIDs,TreeNum,nh
    integer(kind=4),intent(in) :: nhalos(nsteps)
    integer(kind=8),intent(out) :: IDs(nh,nIDs)
    integer(kind=4) :: i,j,ts,hlast
    character(1000) :: filename

    write(filename,'(a,a,i3.3,a)') trim(TreeDir),'/tree_file_',TreeNum,'.001'  ! assume there is a single tree-file 
    open(unit=12,file=filename,status='old',form='unformatted')
    read(12) ! skip nsteps,nIDs,nIndex
    read(12) ! skip nhalos
    hlast = 1
    do ts = 1,nsteps
       if (nhalos(ts) > 0) then 
          read(12) ((IDs(j,i),i=1,nIDs),j=hlast,hlast+nhalos(ts)-1)
          read(12) ! skip indexes
          hlast = hlast+nhalos(ts)
       end if
    end do
    close(12)
    return

  end subroutine read_tree_struct
  
  
  subroutine read_props(TreeDir,TreeNum,nsteps,nhalos,nh,nProps,props)
    
    implicit none

    ! nh is the sum of nhalos ... 
    
    character(1000),intent(in) :: TreeDir
    integer(kind=4),intent(in) :: nsteps,nh,TreeNum,nprops
    integer(kind=4),intent(in) :: nhalos(nsteps)
    real(kind=4),intent(out)   :: props(nh,nProps)
    integer(kind=4) :: i,j,ts,hlast
    character(1000) :: filename

    write(filename,'(a,a,i3.3,a)') trim(TreeDir),'/props_',TreeNum,'.001'  ! assume there is a single tree-file 
    open(unit=12,file=filename,status='old',form='unformatted')
    read(12) ! skip nsteps,nprops
    read(12) ! skip nhalos
    hlast = 1
    do ts = 1,nsteps
       if (nhalos(ts) > 0) then 
          read(12) ((props(j,i),i=1,nProps),j=hlast,hlast+nhalos(ts)-1)
          hlast = hlast+nhalos(ts)
       end if
    end do
    close(12)
    return

  end subroutine read_props
  

  


end module py_tree_ios


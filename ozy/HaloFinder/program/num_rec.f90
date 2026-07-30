! These are subroutines from numerical recipes therefore we refer to it
! for details on what they do 
!***********************************************************************
subroutine qromo(func,omm,oml,a,b,ss,choose)

  use halo_defs, only: dp
  
  implicit none

  integer(kind=4),parameter :: jmax=14,k=5,km=k-1,jmaxp=jmax+1
  integer(kind=4)           :: j
  real(dp),parameter    :: eps=1e-6
  real(dp)              :: s(jmaxp),h(jmaxp),a,b,ss,func,dss,omm,oml
  real(dp) :: zero = 0.d0, one = 1.d0
  external choose      

  h(1) = one

  do j = 1,jmax
     call choose(func,omm,oml,a,b,s(j),j)
     if (j .ge. k) then
        call polint(h(j-km),s(j-km),k,zero,ss,dss)
        if (abs(dss) .lt. eps*abs(ss)) return
     endif
     s(j+1) = s(j)
     h(j+1) = h(j)/9.
  enddo

  stop '> too many steps in qromo'

end subroutine qromo

!***********************************************************************
subroutine polint(xa,ya,n,x,y,dy)
  use halo_defs, only: dp
  implicit none

  integer(kind=4),parameter :: nmax=10
  integer(kind=4)           :: i,m,n,ns  
  real(dp)              :: xa(n),ya(n),c(nmax),d(nmax),x,y,dif,dift
  real(dp)              :: ho,hpa,den,w,dy

  ns  = 1
  dif = abs(x-xa(1))
  do i = 1,n 
     dift = abs(x-xa(i))
     if (dift .lt. dif) then
        ns  = i
        dif = dift
     endif
     c(i) = ya(i)
     d(i) = ya(i)
  enddo

  y  = ya(ns)
  ns = ns-1
  do m = 1,n-1
     do i = 1,n-m
        ho  = xa(i)-x
        hpa = xa(i+m)-x
        w   = c(i+1)-d(i)
        den = ho-hpa
        if (den .eq. 0.) stop
        den  = w/den
        d(i) = hpa*den
        c(i) = ho*den
     enddo
     if (2*ns .lt. n-m)then
        dy = c(ns+1)
     else
        dy = d(ns)
        ns = ns-1
     endif
     y = y+dy
  enddo

  return

end subroutine polint

!***********************************************************************
subroutine midpnt(func,omm,oml,a,b,s,n)
  use halo_defs, only: dp
  implicit none

  real(dp)    ::  func,a,b,s,x,ddel,del,sum,omm,oml
  integer(kind=4) ::  n,tnm,it,j
  save it
  external func

  if (n .eq. 1) then
     s    = (b-a)*func(0.5*(a+b),omm,oml)
     it   = 1
  else
     tnm  = it
     del  = (b-a)/(3.*tnm)
     ddel = del+del
     x    = a+0.5*del
     sum  = 0.
     do j = 1,it
        sum = sum+func(x,omm,oml)
        x   = x+ddel
        sum = sum+func(x,omm,oml)
        x   = x+del
     enddo
     s  = (s+(b-a)*sum/tnm)/3.
     it = 3*it
  endif

  return

end subroutine midpnt

!***********************************************************************
subroutine locate(xx,n,x,j)
  use halo_defs, only: dp
  implicit none

  integer(kind=4) ::  n
  real(dp)    ::  xx(n),x
  integer(kind=4) ::  j,jl,ju,jm

  jl = 0
  ju = n+1
  do while (ju-jl .gt. 1) 
     jm = (ju+jl)/2
     if ((xx(n) .gt. xx(1)) .eqv. (x .gt. xx(jm))) then
        jl = jm
     else
        ju = jm
     endif
  enddo
  j = jl

  return

end subroutine locate

!***********************************************************************
subroutine jacobi(a,np,d,v,nrot)
  use halo_defs, only: dp
  implicit none

  integer(kind=4)           :: np,nrot
  integer(kind=4),parameter :: nmax = 500
  real(dp)              :: a(np,np),d(np),v(np,np)
  integer(kind=4)           :: i,ip,iq,j
  real(dp)              :: c,g,h,s,sm,t,tau,theta,tresh,b(nmax),z(nmax)

  v = 0.
  do ip = 1,np
     v(ip,ip) = 1.
  end do
  z = 0.
  do ip = 1,np
     b(ip) = a(ip,ip)
     d(ip) = b(ip)
  end do
  nrot = 0
  do i = 1,50
     sm = 0.
     do ip = 1,np-1
        do iq = ip+1,np
           sm = sm+abs(a(ip,iq))
        end do
     end do
     if (sm .eq. 0.) return
     if (i .lt. 4) then
        tresh = 0.2*sm/np**2
     else
        tresh = 0.
     endif
     do ip = 1,np-1
        do iq = ip+1,np
           g = 100.*abs(a(ip,iq))
           if ((i .gt. 4) .and. (abs(d(ip))+g .eq. abs(d(ip))) .and. (abs(d(iq))+g .eq. abs(d(iq)))) then
              a(ip,iq) = 0.
           else if (abs(a(ip,iq)) .gt. tresh) then
              h = d(iq)-d(ip)
              if (abs(h)+g .eq. abs(h)) then
                 t     = a(ip,iq)/h
              else
                 theta = 0.5*h/a(ip,iq)
                 t     = 1./(abs(theta)+sqrt(1.+theta**2))
                 if (theta .lt. 0.0) t = -t
              endif
              c        = 1./sqrt(1+t**2)
              s        = t*c
              tau      = s/(1.+c)
              h        = t*a(ip,iq)
              z(ip)    = z(ip)-h
              z(iq)    = z(iq)+h
              d(ip)    = d(ip)-h
              d(iq)    = d(iq)+h
              a(ip,iq) = 0.0
              do j = 1,ip-1
                 g       = a(j,ip)
                 h       = a(j,iq)
                 a(j,ip) = g-s*(h+g*tau)
                 a(j,iq) = h+s*(g-h*tau)
              end do
              do j = ip+1,iq-1
                 g       = a(ip,j)
                 h       = a(j,iq)
                 a(ip,j) = g-s*(h+g*tau)
                 a(j,iq) = h+s*(g-h*tau)
              end do
              do j = iq+1,np
                 g       = a(ip,j)
                 h       = a(iq,j)
                 a(ip,j) = g-s*(h+g*tau)
                 a(iq,j) = h+s*(g-h*tau)
              end do
              do j = 1,np
                 g       = v(j,ip)
                 h       = v(j,iq)
                 v(j,ip) = g-s*(h+g*tau)
                 v(j,iq) = h+s*(g-h*tau)
              end do
              nrot = nrot+1
           endif
        enddo
     enddo
     do ip = 1,np
        b(ip) = b(ip)+z(ip)
        d(ip) = b(ip)
        z(ip) = 0.0
     end do
  end do

  stop 'too many iterations in jacobi'

  return

end subroutine jacobi

!***********************************************************************
subroutine indexx(n,arr,indx)
  use halo_defs, only: dp
  implicit none

  integer(kind=4)           :: n,indx(n)
  real(dp)              :: arr(n)
  integer(kind=4),parameter :: m=7,nstack=50
  integer(kind=4)           :: i,indxt,ir,itemp,j,jstack,k,l,istack(nstack)
  real(dp)              :: a

  do j = 1,n
     indx(j) = j
  enddo

  jstack = 0
  l      = 1
  ir     = n
1 if (ir-l .lt. m) then
     do j = l+1,ir
        indxt = indx(j)
        a     = arr(indxt)
        do i = j-1,1,-1
           if (arr(indx(i)) .le. a) goto 2
           indx(i+1) = indx(i)
        enddo
        i         = 0
2       indx(i+1) = indxt
     enddo
     if (jstack .eq. 0) return
     ir     = istack(jstack)
     l      = istack(jstack-1)
     jstack = jstack-2
  else
     k         = (l+ir)/2
     itemp     = indx(k)
     indx(k)   = indx(l+1)
     indx(l+1) = itemp
     if (arr(indx(l+1)) .gt. arr(indx(ir))) then
        itemp     = indx(l+1)
        indx(l+1) = indx(ir)
        indx(ir)  = itemp
     endif
     if (arr(indx(l)) .gt. arr(indx(ir))) then
        itemp    = indx(l)
        indx(l)  = indx(ir)
        indx(ir) = itemp
     endif
     if (arr(indx(l+1)) .gt. arr(indx(l))) then
        itemp     = indx(l+1)
        indx(l+1) = indx(l)
        indx(l)   = itemp
     endif
     i     = l+1
     j     = ir
     indxt = indx(l)
     a     = arr(indxt)
3    continue
     i     = i+1
     if (arr(indx(i)) .lt. a) goto 3
4    continue
     j     = j-1
     if (arr(indx(j)) .gt. a) goto 4
     if (j .lt. i) goto 5
     itemp   = indx(i)
     indx(i) = indx(j)
     indx(j) = itemp
     goto 3
5    continue
     indx(l) = indx(j)
     indx(j) = indxt
     jstack  = jstack+2
     if (jstack .gt. nstack) stop 'nstack too small in indexx'
     if (ir-i+1 .ge. j-l) then
        istack(jstack)   = ir
        istack(jstack-1) = i
        ir               = j-1
     else
        istack(jstack)   = j-1
        istack(jstack-1) = l
        l                = i
     endif
  endif
  goto 1

end subroutine indexx

!***********************************************************************
function rf(x,y,z)
  use halo_defs, only: dp
  implicit none 

  real(dp)           :: rf,x,y,z
  real(dp),parameter :: errtol=.08,tiny=1.5e-38,big=3.e37,third=1./3.,c1=1./24.,c2=.1,c3=3./44.,c4=1./14.
  real(dp)           :: alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt

  if (min(x,y,z) < 0.0 .or. min(x+y,x+z,y+z) < tiny .or. max(x,y,z) > big) stop 'invalid arguments in rf'

  xt    = x
  yt    = y
  zt    = z
1 continue
  sqrtx = sqrt(xt)
  sqrty = sqrt(yt)
  sqrtz = sqrt(zt)
  alamb = sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
  xt    = .25*(xt+alamb)
  yt    = .25*(yt+alamb)
  zt    = .25*(zt+alamb)
  ave   = third*(xt+yt+zt)
  delx  = (ave-xt)/ave
  dely  = (ave-yt)/ave
  delz  = (ave-zt)/ave
  if (max(abs(delx),abs(dely),abs(delz)) .gt. errtol) goto 1
  e2    = delx*dely-delz**2
  e3    = delx*dely*delz
  rf    = (1.+(c1*e2-c2-c3*e3)*e2+c4*e3)/sqrt(ave)

  return

end function rf

!***********************************************************************
function cubic(a1,a2,a3,a4)
  use halo_defs, only: dp
!     Function for evaluating a cubic equation.  In the case where
!     there is one real root, this routine returns it.  If there are three
!     roots, it returns the smallest positive one.
!
!     The cubic equation is a1*x**3 + a2*x**2 + a3*x + a4 = 0.
!
!     Ref: Tallarida, R. J., "Pocket Book of Integrals and Mathematical
!     Formulas," 2nd ed., Boca Raton: CRC Press 1992, pp. 8-9.

  implicit none

  real(dp) :: a1,a2,a3,a4,cubic
  real(dp) :: a,b,d,r,theta,pi,ar,ai,y,y1,y2,y3,p,q

  if (a1 .eq. real(0.0d0,kind=dp)) then
     stop '> ERROR: Quadratic/linear equation passed to cubic function'
  end if

  p = a3/a1 - (a2/a1)**2/3.0d0
  q = a4/a1 - a2*a3/a1**2/3.0d0 + 2.0d0*(a2/a1)**3/27.0d0
  d = p**3/27.0d0 + q**2/4.0d0

  if (d .gt. 0.0d0) then
     a = -0.5*q + sqrt(d)
     b = -0.5*q - sqrt(d)
     if (a .gt. 0.0) then
        a = a**(1.0d0/3.0d0)
     else
        a = - (-a)**(1.0d0/3.0d0)
     end if
     if (b .gt. 0.0) then
        b = b**(1.d0/3.d0)
     else
        b = - (-b)**(1.d0/3.d0)
     end if
     y = a + b
  else
     ar    = -0.5d0*q
     ai    = sqrt(-d)
     r     = (ar**2+ai**2)**(1.0d0/6.0d0)
!    atan2 is arctan
     theta = atan2(ai,ar)
     y1    = 2.0d0 * r * cos(theta/3.0d0) - a2/a1/3.0d0
     y     = y1
     pi    = 4.0d0*atan(1.0d0)
     y2    = 2.0d0 * r * cos(theta/3.0d0+2.0d0*pi/3.0d0) - a2/a1/3.0d0
     if (y .lt. 0.0d0 .or. (y2 .gt. 0.0d0 .and. y2 .lt. y)) y = y2
     y3    = 2.0d0 * r * cos(theta/3.0d0-2.d0*pi/3.0d0) - a2/a1/3.0d0
     if (y .lt. 0.0d0 .or. (y3 .gt. 0.0d0 .and. y3 .lt. y)) y = y3
  end if

  cubic = y

  return

end function cubic

!********************************************************************
subroutine int_indexx(n,arr,indx)

  integer :: n,indx(n),M,NSTACK
  integer :: arr(n)
  parameter (M=7,NSTACK=50)
  integer :: i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
  integer :: a
  do j=1,n
     indx(j)=j
  enddo
  jstack=0
  l=1
  ir=n
1 if(ir-l.lt.M)then
     do j=l+1,ir
        indxt=indx(j)
        a=arr(indxt)
        do i=j-1,1,-1
           if(arr(indx(i)).le.a)goto 2
           indx(i+1)=indx(i)
        enddo
        i=0
2       indx(i+1)=indxt
     enddo
     if(jstack==0)return
     ir=istack(jstack)
     l=istack(jstack-1)
     jstack=jstack-2
  else
     k=(l+ir)/2
     itemp=indx(k)
     indx(k)=indx(l+1)
     indx(l+1)=itemp
     if(arr(indx(l+1)).gt.arr(indx(ir)))then
        itemp=indx(l+1)
        indx(l+1)=indx(ir)
        indx(ir)=itemp
     endif
     if(arr(indx(l)).gt.arr(indx(ir)))then
        itemp=indx(l)
        indx(l)=indx(ir)
        indx(ir)=itemp
     endif
     if(arr(indx(l+1)).gt.arr(indx(l)))then
        itemp=indx(l+1)
        indx(l+1)=indx(l)
        indx(l)=itemp
     endif
     i=l+1
     j=ir
     indxt=indx(l)
     a=arr(indxt)
3    continue
     i=i+1
     if(arr(indx(i)).lt.a)goto 3
4    continue
     j=j-1
     if(arr(indx(j)).gt.a)goto 4
     if(j.lt.i)goto 5
     itemp=indx(i)
     indx(i)=indx(j)
     indx(j)=itemp
     goto 3
5    indx(l)=indx(j)
     indx(j)=indxt
     jstack=jstack+2
     if(jstack.gt.nstack) stop 'nstack too small in indexx'
     if(ir-i+1.ge.j-l)then
        istack(jstack)=ir
        istack(jstack-1)=i
        ir=j-1
     else
        istack(jstack)=j-1
        istack(jstack-1)=l
        l=i
     endif
  endif
  goto 1
end subroutine int_indexx



!********************************************************************
FUNCTION ran2(idum)
  use halo_defs, only: dp
  integer :: idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
  real(dp) :: ran2,AM,EPS,RNMX
  parameter (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
 &   IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, &
 &   NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
  integer :: idum2,j,k,iv(NTAB),iy
  save iv,iy,idum2
  data idum2/123456789/, iv/NTAB*0/, iy/0/

  if (idum.le.0) then
     idum=max(-idum,1)
     idum2=idum
     do j=NTAB+8,1,-1
        k=idum/IQ1
        idum=IA1*(idum-k*IQ1)-k*IR1
        if (idum.lt.0) idum=idum+IM1
        if (j.le.NTAB) iv(j)=idum
     enddo
     iy=iv(1)
  endif
  k=idum/IQ1
  idum=IA1*(idum-k*IQ1)-k*IR1
  if (idum.lt.0) idum=idum+IM1
  k=idum2/IQ2
  idum2=IA2*(idum2-k*IQ2)-k*IR2
  if (idum2.lt.0) idum2=idum2+IM2
  j=1+iy/NDIV
  iy=iv(j)-idum2
  iv(j)=idum
  if(iy.lt.1)iy=iy+IMM1
  ran2=min(AM*iy,RNMX)
end FUNCTION ran2

!********************************************************************
SUBROUTINE mrqmin(x,y,sig,ndata,a,ia,ma,covar,alpha,nca,chisq,funcs,alamda)
  use halo_defs, only: dp
  
  implicit none 
  
  INTEGER      :: ma,nca,ndata,ia(ma),MMAX
  REAL(dp) :: alamda,chisq,a(ma),alpha(nca,nca),covar(nca,nca),sig(ndata),x(ndata),y(ndata)
  PARAMETER (MMAX=20)
  !    USES covsrt,gaussj,mrqcof
  INTEGER      :: j,k,l,m,mfit
  REAL(dp) :: ochisq,atry(MMAX),beta(MMAX),da(MMAX)
  SAVE ochisq,atry,beta,da,mfit
  external funcs  
  
  if(alamda.lt.0.)then
     mfit=0
     do j=1,ma
        if (ia(j).ne.0) mfit=mfit+1
     enddo
     alamda=0.001
     call mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,nca,chisq,funcs)
     ochisq=chisq
     do j=1,ma
        atry(j)=a(j)
     enddo
  endif
  j=0
  do l=1,ma
     if(ia(l).ne.0) then
        j=j+1
        k=0
        do m=1,ma
           if(ia(m).ne.0) then
              k=k+1
              covar(j,k)=alpha(j,k)
           endif
        enddo
        covar(j,j)=alpha(j,j)*(1.+alamda)
        da(j)=beta(j)
     endif
  enddo
  call gaussj(covar,mfit,nca,da,1,1)
  if(alamda.eq.0.)then
     call covsrt(covar,nca,ma,ia,mfit)
     return
  endif
  j=0
  do l=1,ma
     if(ia(l).ne.0) then
        j=j+1
        atry(l)=a(l)+da(j)
     endif
  enddo

  ! in our case, the concentration parameter must be > 0.0 or the scale radius r_s becomes infinite
  ! in practice we do not allow r_s to become greater than r_vir, i.e. c must be greater than 1.0
  j=0
  do l=1,ma
     if(ia(l).ne.0) then
        j=j+1
        atry(l)= max(atry(l),1.0)
        da(j)  = atry(l) - a(l) 
     endif
  enddo
  
!  write(132,*) 'before mrqcof',a(1),atry(1),ia(1),da(1)
  call mrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,nca,chisq,funcs)
!  write(132,'(a15,1x,2(f6.3,1x),i1,2(e15.6,1x))') 'c,ochisq,chisq:',a(1),atry(1),ia(1),ochisq,chisq
  if(chisq.lt.ochisq)then
     alamda=0.1*alamda
     ochisq=chisq
     j=0
     do l=1,ma
        if(ia(l).ne.0) then
           j=j+1
           k=0
           do m=1,ma
              if(ia(m).ne.0) then
                 k=k+1
                 alpha(j,k)=covar(j,k)
              endif
           enddo
           beta(j)=da(j)
           a(l)=atry(l)
        endif
     enddo
  else
     alamda=10.*alamda
     chisq=ochisq
  endif
  return
  
end SUBROUTINE mrqmin

!*********************************************************************
SUBROUTINE mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,nalp,chisq,funcs)
   use halo_defs, only: dp
      
  implicit none

  INTEGER        :: ma,nalp,ndata,ia(ma),MMAX
  REAL(dp)   :: chisq,a(ma),alpha(nalp,nalp),beta(ma),sig(ndata),x(ndata),y(ndata)
  PARAMETER (MMAX=20)
  INTEGER        :: mfit,i,j,k,l,m
  REAL(dp)   :: dy,sig2i,wt,ymod,dyda(MMAX)
  external funcs 
  
  mfit=0
  do  j=1,ma
     if (ia(j).ne.0) mfit=mfit+1
  enddo
  
  do j=1,mfit
     do k=1,j
        alpha(j,k)=0.
     enddo
     beta(j)=0.
  enddo
  
  chisq=0.
  do i=1,ndata
     call funcs(x(i),a,ymod,dyda,ma)
!     write(132,'(a6,1x,i3,1x,5(e15.8,1x))') 'funcs:',i,x(i),y(i),sig(i),ymod,dyda(1)
     sig2i=1./(sig(i)*sig(i))
     dy=y(i)-ymod
     j=0
     do l=1,ma
        if(ia(l).ne.0) then
           j=j+1
           wt=dyda(l)*sig2i
           k=0
           do m=1,l
              if(ia(m).ne.0) then
                 k=k+1
                 alpha(j,k)=alpha(j,k)+wt*dyda(m)
              endif
           enddo
           beta(j)=beta(j)+dy*wt
        endif
     enddo
     chisq=chisq+dy*dy*sig2i
  enddo
  
  do  j=2,mfit
     do  k=1,j-1
        alpha(k,j)=alpha(j,k)
     enddo
  enddo
  
  return
  
end SUBROUTINE mrqcof

!*********************************************************************
SUBROUTINE gaussj(a,n,np,b,m,mp)
   use halo_defs, only: dp

  implicit none 

  INTEGER      :: m,mp,n,np,NMAX
  REAL(dp) :: a(np,np),b(np,mp)
  PARAMETER (NMAX=50)
  INTEGER      :: i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
  REAL(dp) :: big,dum,pivinv

  do j=1,n
     ipiv(j)=0
  enddo

  do i=1,n

     big=0.
     do j=1,n
        if(ipiv(j).ne.1)then
           do k=1,n
              if (ipiv(k).eq.0) then
                 if (abs(a(j,k)).ge.big)then
                    big=abs(a(j,k))
                    irow=j
                    icol=k
                 endif
              else if (ipiv(k).gt.1) then
                 pause 'singular matrix in gaussj'
              endif
           enddo
        endif
     enddo

     ipiv(icol)=ipiv(icol)+1
     if (irow.ne.icol) then
        do l=1,n
           dum=a(irow,l)
           a(irow,l)=a(icol,l)
           a(icol,l)=dum
        enddo
        do l=1,m
           dum=b(irow,l)
           b(irow,l)=b(icol,l)
           b(icol,l)=dum
        enddo
     endif
     indxr(i)=irow
     indxc(i)=icol
     if (a(icol,icol).eq.0.) pause 'singular matrix in gaussj'
     pivinv=1./a(icol,icol)
     a(icol,icol)=1.

     do l=1,n
        a(icol,l)=a(icol,l)*pivinv
     enddo

     do l=1,m
        b(icol,l)=b(icol,l)*pivinv
     enddo

     do ll=1,n
        if(ll.ne.icol)then
           dum=a(ll,icol)
           a(ll,icol)=0.
           do l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
           enddo
           do l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
           enddo
        endif
     enddo
     
  enddo

  do l=n,1,-1
     if(indxr(l).ne.indxc(l))then
        do k=1,n
           dum=a(k,indxr(l))
           a(k,indxr(l))=a(k,indxc(l))
           a(k,indxc(l))=dum
        enddo
     endif
  enddo
  
  return

end SUBROUTINE gaussj

!**********************************************************************************
 
SUBROUTINE covsrt(covar,npc,ma,ia,mfit)
   use halo_defs, only: dp

  implicit none 

  INTEGER      :: ma,mfit,npc,ia(ma)
  REAL(dp) :: covar(npc,npc)
  INTEGER      :: i,j,k
  REAL(dp) :: swap

  do i=mfit+1,ma
     do j=1,i
        covar(i,j)=0.
        covar(j,i)=0.
     enddo
  enddo

  k=mfit
  do j=ma,1,-1
     if(ia(j).ne.0)then
        do i=1,ma
           swap=covar(i,k)
           covar(i,k)=covar(i,j)
           covar(i,j)=swap
        enddo
        do i=1,ma
           swap=covar(k,i)
           covar(k,i)=covar(j,i)
           covar(j,i)=swap
        enddo
        k=k-1
     endif
  enddo  
  
  return

end SUBROUTINE covsrt


 


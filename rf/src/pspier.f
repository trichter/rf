C file: pspier.f
C
C  The program calculates piercing points of P-to-S converted phases
C  under the station at a given depth.
C  P wave slownesses (in sec/deg) and conversion depths (km) are read
C  from a file (DSLOW.STX)
C  The velocity model: /home/st/yuan/model/iasp91.dat
C  Output file: PSPIER.STX
C
C  Copyright:
C  Xiaohui Yuan, April 1997
C  Tom Richter, January 2011: changes for use with f2py
C                             use model ../data/iasp91.dat
C
C  MIT license
C
      subroutine pspier(rpier,plat,plon,n,depth,slat,slon,slow,azim)

      
      integer x
      integer n
      real*8 depth
      real*8 slat(n)
      real*8 slon(n)
      real*8 slow(n)
      real*8 azim(n)
      real*8 rpier(n)
      real*8 plat(n)
      real*8 plon(n)
      dimension depco(500),vpco(500),vsco(500),layco(500)
      dimension dep(2000),vp(2000),vs(2000),hk(2000),ro(2000)
c	dimension slat(10000),slon(10000),slow(10000),azim(10000)
c	dimension rpier(10000),plat(10000),plon(10000)

Cf2py intent(in) n
Cf2py intent(in) depth
Cf2py intent(in) slat
Cf2py intent(in) slon
Cf2py intent(in) slow
Cf2py intent(in) azim
Cf2py intent(out) rpier
Cf2py intent(out) plat
Cf2py intent(out) plon

      pi=3.1415927
      arc=pi/180.0
      depthmax = depth
      nray = n
c	reading conversion depth and P wave slowness from file
C   open(13,file='DSLOW.STX',status='old')
C   read(13,*) nray,depthmax
C   do 5 i=1,nray
C5  read(13,*) slat(i),slon(i),slow(i),azim(i)
C   close(13)

c      write(6,91)nray,depthmax
c91    format(i5,' rays,','  conversion at depth',f7.1,'km.')
c      write(6,*) 'output to file PSPIER.STX'

c	reading velocity model from file
      open(11,file='../data/iasp91.dat',status='old')
      read(11,*)
      read(11,*)
      nlayco=0
      do i=1,500
        read(11,*) depco(i),vpco(i),vsco(i),layco(i)
        if (abs(depco(i)).lt.0.001.and.abs(vpco(i)).lt.0.001) then
          if (abs(vsco(i)).lt.0.001.and.abs(layco(i)).lt.0.001) then
               goto 11
          end if
        end if
        nlayco=nlayco+1
      end do
11    close(11)
c
c	transformation of gradient model to homogeneous layered model
c
      lay=0
      do i=1,nlayco
        if(layco(i).eq.0) goto 20
        vp0=vpco(i-1)
        vs0=vsco(i-1)
        dep0=depco(i-1)
        devp=(vpco(i)-vpco(i-1))/layco(i)
        devs=(vsco(i)-vsco(i-1))/layco(i)
        dedep=(depco(i)-depco(i-1))/layco(i)
        do j=1,layco(i)
          lay=lay+1
          vp(lay)=vp0+devp*(j-0.5)
          vs(lay)=vs0+devs*(j-0.5)
          hk(lay)=dedep
          dep(lay)=dep0+dedep*(j-1)
          depm=dep(lay)+dedep
          if (depm.ge.depthmax) then
            hk(lay)=hk(lay)-(depm-depthmax)
            goto 21
          end if
        end do
20    x = 0
      end do
c  write out to file
c	open(12,file='model.lay',status='unknown')
c	write(12,*) 'layered model'
c	write(12,*) '     depth      vp        vs        h'
c	do 30 i=1,lay
c30	write(12,'(4f10.2)') dep(i),vp(i),vs(i),hk(i)
c	close(12)

c
c  Earth Flattening Approximation
c
21    call efa(lay,dep,hk,vp,vs,ro,1)
c
c  Calculation of piercing points
c
c   loop for rays
      do kray=1,nray
        rpier(kray)=0

c   loop for model layers
        do i=1,lay
          pv=slow(kray)*vs(i)/111.2
          rpier(kray)=rpier(kray)+hk(i)*pv/sqrt(1-pv*pv)
        end do
        alpha=azim(kray)*arc
        rdeg=rpier(kray)/111.2
        plat(kray)=slat(kray)+rdeg*cos(alpha)
        plon(kray)=slon(kray)+rdeg*sin(alpha)/cos(slat(kray)*arc)
      enddo
c
c  Output
c
c      write(6,97)
c97    format(3x,'trace',3x,'slat',6x,'slon',7x,'r',8x,'plat',6x,'plon')
c      do i=1,nray
c        write(6,99) i,slat(i),slon(i),rpier(i),plat(i),plon(i)
c      end do


c    open(14,file='PSPIER.STX',status='unknown')
c    do 200 i=1,nray
c200	write(14,99) i,slat(i),slon(i),rpier(i),plat(i),plon(i)
c    close(14)
c99    format(1x,i5,2f10.3,f10.1,2f12.5)
C    stop
      return
      end

C ********************************************

      subroutine efa(n,dep,h,vp,vs,ro,i)
c
c--- earth flattenning approximation for gradient model
c--- i=1(from spherical to flat), i.ne.1 (inverse)
c--- modefied from kosarev's "efal" in "bip.f"
c

      dimension dep(n),h(n),vp(n),vs(n),ro(n)

Cf2py intent(in) n
Cf2py intent(in,out) dep
Cf2py intent(in,out) h
Cf2py intent(in,out) vp
Cf2py intent(in,out) vs
Cf2py intent(in,out) ro

      if(i.eq.1) then
        za=0.
        zz=0.
        do k=1,n
c      write(*,*)'k=',k,' y1',y1,' dep,vp,vs,ro,h'
c     &,dep(k),vp(k),vs(k),ro(k),h(k)
              z1k=zz+h(k)*0.5
              zz=zz+h(k)
              y1=6370./(6370.-z1k)
              yz=6370./(6370.-zz)
              zm=6370.*alog(yz)
              h(k)=zm-za
              za=zm
              vp(k)=vp(k)*y1
              vs(k)=vs(k)*y1
              ro(k)=ro(k)/y1
              if (k.gt.1) dep(k)=dep(k-1)+h(k-1)
c      write(*,*)'k=',k,' y1',y1,' dep,vp,vs,ro,h'
c     &,dep(k),vp(k),vs(k),ro(k),h(k)
        end do
      else
        za=0.
        zz=0.
        do k=1,n
              z1k=zz+h(k)*0.5
              if(k.eq.n) z1k=zz
              zz=zz+h(k)
              z1k=6370.*(1.-exp(-z1k/6370.))
              zm=6370.*(1.-exp(-zz/6370.))
              y1=6370./(6370.-z1k)
              h(k)=zm-za
              za=zm
              vp(k)=vp(k)/y1
              vs(k)=vs(k)/y1
              ro(k)=ro(k)*y1
              if (k.gt.1) dep(k)=dep(k-1)+h(k-1)
        end do  
      end if
      return
      end
C end file pspier.f
C  sppier.f
C
C  The program calculates piercing points of S-to-P converted phases
C  under the station.
C  S wave slownesses (in sec/deg) and conversion depths (km) are read
C  from a file (DSLOW.STX)
C  The velocity model is in a seperate file, data are in the order of 
C  "depth vp vs no_of_layers", starting at surface (iasp91.dat).
C  Output file: SPPIER.STX
C
C  Xiaohui Yuan, April 1997, Aug. 2003
C
C  Tom Richter, 2013: changes for use with f2py
C                     use model ../data/iasp91.dat
C
C  MIT license
C
      subroutine sppier(rpier,plat,plon,n,depth,slat,slon,slow,azim)

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
C dimension slat(500),slon(500),slow(500),azim(500)
C dimension rpier(500),plat(500),plon(500)

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
c
c reading conversion depth and P wave slowness from file
C open(13,file='DSLOW.STX',status='old')
C read(13,*) nray,depthmax
C do 5 i=1,nray
C5  read(13,*) slat(i),slon(i),slow(i),azim(i)
C close(13)

C write(6,91)nray,depthmax
C91 format(i5,' rays,','  conversion at depth',f7.1,'km.')
C write(6,*) 'output to file SPPIER.STX'

c reading velocity model from file
      open(11,file='../data/iasp91.dat',status='old')
      read(11,*)
      read(11,*)
      nlayco=0
      do 10 i=1,500
      read(11,*) depco(i),vpco(i),vsco(i),layco(i)
      if (abs(depco(i)).lt.0.001.and.abs(vpco(i)).lt.0.001) then
         if (abs(vsco(i)).lt.0.001.and.abs(layco(i)).lt.0.001) then
           goto 11
         end if
      end if
      nlayco=nlayco+1
10    continue
11    close(11)

c
c transformation of gradient model to homogeneous layered model
c
      lay=0
      do 20 i=1,nlayco
      if(layco(i).eq.0) goto 20
      vp0=vpco(i-1)
      vs0=vsco(i-1)
      dep0=depco(i-1)
      devp=(vpco(i)-vpco(i-1))/layco(i)
      devs=(vsco(i)-vsco(i-1))/layco(i)
      dedep=(depco(i)-depco(i-1))/layco(i)
      do 19 j=1,layco(i)
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
19    continue
20    continue
21    continue
c  write out to file
C open(12,file='/home/st/felix/model/iasp91.dat',status='unknown')
C write(12,*) 'layered model'
C write(12,*) '     depth      vp        vs        h'
C do 30 i=1,lay
C30 write(12,'(4f10.2)') dep(i),vp(i),vs(i),hk(i)
C close(12)

c
c  Earth Flattening Approximation
c
      call efa(lay,dep,hk,vp,vs,ro,1)

c
c  Calculation of piercing points
c
c loop for rays
      do 110 kray=1,nray
      rpier(kray)=0

c loop for model layers
      do 100 i=1,lay
      pv=slow(kray)*vp(i)/111.2
      rpier(kray)=rpier(kray)+hk(i)*pv/sqrt(1-pv*pv)
100   continue

      alpha=azim(kray)*arc
      rdeg=rpier(kray)/111.2
      plat(kray)=slat(kray)+rdeg*cos(alpha)
      plon(kray)=slon(kray)+rdeg*sin(alpha)/cos(slat(kray)*arc)
110   continue
c
c  Output
c
C write(6,97)
C97 format(3x,'trace',3x,'slat',6x,'slon',7x,'r',8x,'prclat'
C     &  ,6x,'prclon')         
C do 201 i=1,nray
C201  write(6,99) i,slat(i),slon(i),rpier(i),plat(i),plon(i)

C open(14,file='SPPIER.STX',status='unknown')
C do 200 i=1,nray
C200  write(14,99) i,slat(i),slon(i),rpier(i),plat(i),plon(i)
C close(14)
C99 format(1x,i5,2f10.3,f10.1,2f12.5)
C stop
      return
      end

C********************************************
C
C      subroutine efa(n,dep,h,vp,vs,ro,i)
Cc
Cc--- earth flattenning approximation for gradient model
Cc--- i=1(from spherical to flat), i.ne.1 (inverse)
Cc--- modefied from kosarev's "efal" in "bip.f"
Cc
C      dimension dep(n),h(n),vp(n),vs(n),ro(n)
C      if(i.eq.1) then
C         za=0.
C         zz=0.
C         do k=1,n
Cc      write(*,*)'k=',k,' y1',y1,' dep,vp,vs,ro,h'
Cc     &,dep(k),vp(k),vs(k),ro(k),h(k)
C              z1k=zz+h(k)*0.5
C              zz=zz+h(k)
C              y1=6370./(6370.-z1k)
C              yz=6370./(6370.-zz)
C              zm=6370.*alog(yz)
C              h(k)=zm-za
C              za=zm
C              vp(k)=vp(k)*y1
C              vs(k)=vs(k)*y1
C              ro(k)=ro(k)/y1
C              if (k.gt.1) dep(k)=dep(k-1)+h(k-1)
Cc      write(*,*)'k=',k,' y1',y1,' dep,vp,vs,ro,h'
Cc     &,dep(k),vp(k),vs(k),ro(k),h(k)
C         end do
C      else
C         za=0.
C         zz=0.
C         do k=1,n
C              z1k=zz+h(k)*0.5
C              if(k.eq.n) z1k=zz
C              zz=zz+h(k)
C              z1k=6370.*(1.-exp(-z1k/6370.))
C              zm=6370.*(1.-exp(-zz/6370.))
C              y1=6370./(6370.-z1k)
C              h(k)=zm-za
C              za=zm
C              vp(k)=vp(k)/y1
C              vs(k)=vs(k)/y1
C              ro(k)=ro(k)*y1
C              if (k.gt.1) dep(k)=dep(k-1)+h(k-1)
C         end do  
C       end if
C       return
C       end

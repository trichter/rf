C  psmout.f
C
C  The program correct moveout for Ps/Ppps/Ppss/ Psss phases.
C  Reference slowness is 6.4 s/deg.
C  It is assumed that Ps and the other multiples have the same slowness
C  with the corresponding P phase which may have big deviations for 
C  multiples at deep depth.
C  It uses Earth-Flattening-Approximation.
C
C  Input data are ASCII written from SeismicHandler with the header: 
C  SLOWNESS" (TMP_PSMOUT.ASC)
C
c  Reference earth model: /home/st/yuan/model/iasp91.dat
c                         or a local file.
C
C  Output: binary files in SeismicHandler format.
C  (PSMOUT.Q??, PPPSMOUT.Q??, PPSSMOUT.Q??, PSSSMOUT.Q??)
C
C  Copyright:
C  Xiaohui Yuan, Sep. 1997.
C  XY, July 2003
C  XY, 02.11.2004
C  Tom Richter, January 2011: changes for use with f2py
C                             use model ../data/iasp91.dat
C
C  MIT license
C
      subroutine psmout(ntrace,nsample,x,slow,tbeg,tend,delt,itype,y)

      integer ntrace
      integer nsample
      integer itype

      real*8 tbeg
      real*8 tend
      real*8 delt
      real*8 delta
      real*8 slow(ntrace)
C      real*8 x(nsample)

      dimension ttmo(401,2001),tt(2001),maxt(401),maxy(10000)
      dimension x(ntrace, nsample), y(ntrace, nsample)
C      dimension x(20000), y(20000),ysum(20000),tref(20000)
      dimension ysum(20000),tref(20000)
      character model*50
      common /mod/ model
      common /slow/ slow1,slow2,dslow
      common /time/ time1,time2,dtime

Cf2py intent(in) ntrace
Cf2py intent(in) nsample
Cf2py intent(in) x
Cf2py intent(in) tbeg
Cf2py intent(in) tend
Cf2py intent(in) delt
Cf2py intent(in) slow
Cf2py intent(in) itype
Cf2py intent(out) y

      length = nsample
      delta=delt
c
c  open input and output data file and read some values
c

C    open(31,file='TMP_PSMOUT.ASC',status='old')
C    read(31,'(a)') name
C    read(31,*) ntrace,tbeg,tend,delt,itype
      numt=int((tend-tbeg)/delt+0.5)+1

      model='../data/iasp91.dat'
    
      do 4 j=1,numt
4       tref(j)=tbeg+(j-1)*delt

c
c   open output ascii files for different phases:
c     itype=1: Ps
c     itype=2: P2p1s
c     itype=3: P1p2s
c     itype=4: P3s
c
C23456789012345678901234567890123456789012345678901234567890123456789012
C    if (itype.eq.1) then
C      open(32,file='PSMOUT.QBN',status='unknown'
C    &             ,form='unformatted',access='direct',recl=4)
C    else if (itype.eq.2) then
C      open(32,file='PPPSMOUT.QBN',status='unknown'
C    &             ,form='unformatted',access='direct',recl=4)
C    else if (itype.eq.3) then
C      open(32,file='PPSSMOUT.QBN',status='unknown'
C    &             ,form='unformatted',access='direct',recl=4)
C    else if (itype.eq.4) then
C      open(32,file='PSSSMOUT.QBN',status='unknown'
C    &             ,form='unformatted',access='direct',recl=4)
C    else
C      write(*,*) 'error: wrong phase type'
C      stop
C    endif

c
c  moveout table
c
      call table(itype,nslow,maxt,ttmo)
c
      do 6 j=1,numt
6       ysum(j)=0
      do 8 i=1,ntrace
8       maxy(i)=numt
c
c  Loop for traces
c
      maxtr=0
      do 100 itr=1,ntrace
        write(*,96)itr
96      format('  trace: ',i5)

c
c  read in data of one trace
c
C5       read(31,*) ch,delta
C        if(ch.ne.'DELTA:') goto 5
C        read(31,*) ch,length
C        read(31,*) ch,slowness
C        read(31,*) (x(j),j=1,length)
        slowness=slow(itr)


        if (slowness.gt.slow2) then
          write(*,89) slow2
89          format(' slowness exceeds',f10.3,', neglected')
          goto 100
        endif
        maxtr=maxtr+1
c
c  moveout table for this trace
c
        is1=int((slowness-slow1)/dslow)+1
        if(is1.eq.nslow) is1=is1-1
        is2=is1+1
        si1=slow1+(is1-1)*dslow
        si2=si1+dslow
        maxtable=maxt(is2)
        do 10 j=1,maxtable
C        tt(j)=ttmo(is1,j) + (ttmo(is2,j)-ttmo(is1,j))*(slowness-si1)/(si2-si1)
        tt(j)=(ttmo(is2,j)-ttmo(is1,j))*(slowness-si1)/(si2-si1)
        tt(j)=ttmo(is1,j) + tt(j) 
10      continue
c
c  moveout correction
c
        do 20 j=1,numt
          if (tref(j).le.0) then
            treal=tref(j)
          else
            i1=int((tref(j)-time1)/dtime)+1
            i2=i1+1
            if (i2.gt.maxtable) then
                write(*,911) treal
911             format('  warning: maximum Pds delay time: ',f10.3,'s')
                maxy(itr)=j-1
                goto 21
            endif
            ti1=time1+(i1-1)*dtime
            ti2=ti1+dtime
            ttx=tt(i1)+(tref(j)-ti1)*(tt(i2)-tt(i1))/(ti2-ti1)
            treal=tref(j)+ttx
          end if

          i1=int((treal-tbeg)/delta+0.5)+1

          if (i1.gt.length) then
            write(*,*) '  warning: data not enough length'
            maxy(itr)=j-1
            goto 21
          else if (i1.eq.length) then
            i1=i1-1
          else
            i2=i1+1
            ti1=tbeg+(i1-1)*delta
            ti2=ti1+delta
C            y(j)=x(i1)+(treal-ti1)*(x(i2)-x(i1))/(ti2-ti1)
            y(itr,j)=(treal-ti1)*(x(itr,i2)-x(itr,i1))
            y(itr,j)=x(itr,i1)+y(itr,j)/(ti2-ti1)
C            ysum(j)=ysum(j)+y(j)
          endif
20      continue
21      continue
c
c  output to file
c
        do 30 j=1,numt
C          if (j.gt.maxy(itr)) y(j)=0
          if (j.gt.maxy(itr)) y(j,itr)=0
30      continue
C        do 40 j=1,numt
C          irecnum=j+(itr-1)*numt
C40      write(32,rec=irecnum) y(j)

100   continue

C    stop
      return
      end


********************************************************************

      subroutine table(itype,nslow,maxt,ttpsmo)
      common /slow/ slow1,slow2,dslow
      common /time/ time1,time2,dtime
      common /mod/ model
      dimension depco(500),vpco(500),vsco(500),layco(500)
      dimension dep(5000),vp(5000),vs(5000),hk(5000),ro(5000)
      dimension slowness(401),depth(2001),ttps(401,2001),ttps0(2001)
      dimension t0(2001),ttpsmo(401,2001),maxdep(401),maxt(401)
      character model*50
c
c  parameters for moveout table
c
      slow1=0
      slow2=15
      dslow=0.1
      nslow=int((slow2-slow1)/dslow)+1

      depth1=0
      depth2=2000
      ddepth=2
      ndepth=int((depth2-depth1)/ddepth)+1

      time1=0
      time2=250
      dtime=0.2
      ntime=int((time2-time1)/dtime)+1

      do 3 i=1,nslow
3     maxt(i)=ntime

c
c  read in velocity model
c
      open(11,file=model,status='old')
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

      maxdep0=0
      do 301 i=1,nslow
301   maxdep(i)=0

c
c  Loop for depth
c
      do 200 ldep=1,ndepth
      depthmax=depth1+ddepth*(ldep-1)
      depth(ldep)=depthmax
c
c  transformation of gradient model to homogeneous layered model
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
c
c  Earth Flattening Approximation
c
      call efa(lay,dep,hk,vp,vs,ro,1)
c
c  for reference slowness p=6.4
c
      slow=6.4/111.195
      ttp=0
      tts=0
      xp=0
      xs=0

      do 80 i=1,lay
      pvp=slow*vp(i)
      pvs=slow*vs(i)

c   check for post critical angle
      if (pvp.ge.1) goto 201

      h=hk(i)
      termp=1/sqrt(1-pvp*pvp)
      terms=1/sqrt(1-pvs*pvs)
      ttp=ttp+h*termp/vp(i)
      tts=tts+h*terms/vs(i)
      xp=xp+h*pvp*termp
      xs=xs+h*pvs*terms
80    continue
c
c   computes delay time for different phase type:
c     itype=1: Ps
c     itype=2: P2p1s
c     itype=3: P1p2s
c     itype=4: P3s
c
      if (itype.eq.1) then
        ttps0(ldep)=tts-ttp-(xs-xp)*slow
      else if (itype.eq.2) then
        ttps0(ldep)=tts+ttp-(xs+xp)*slow
      else if (itype.eq.3) then
        ttps0(ldep)=2*tts-2*xs*slow
      else if (itype.eq.4) then
        ttps0(ldep)=3*tts-ttp-(3*xs-xp)*slow
      else
        write(*,*) 'error in subroutine: wrong phase type'
        stop
      endif

      maxdep0=maxdep0+1

c
c  Main calculation, loop for slowness
c
      do 101 lslo=1,nslow
      slowness(lslo)=slow1+dslow*(lslo-1)
      slow=slowness(lslo)/111.195
      ttp=0
      tts=0
      xp=0
      xs=0

      do 100 i=1,lay
      pvp=slow*vp(i)
      pvs=slow*vs(i)

c   check for post critical angle
      if (pvp.ge.1) goto 101

      h=hk(i)
      termp=1/sqrt(1-pvp*pvp)
      terms=1/sqrt(1-pvs*pvs)
      ttp=ttp+h*termp/vp(i)
      tts=tts+h*terms/vs(i)
      xp=xp+h*pvp*termp
      xs=xs+h*pvs*terms
100   continue
c
c   computes delay time for different phase type:
c     itype=1: Ps
c     itype=2: P2p1s
c     itype=3: P1p2s
c     itype=4: P3s
c
      if (itype.eq.1) then
        ttps(lslo,ldep)=tts-ttp-(xs-xp)*slow
      else if (itype.eq.2) then
        ttps(lslo,ldep)=tts+ttp-(xs+xp)*slow
      else if (itype.eq.3) then
        ttps(lslo,ldep)=2*tts-2*xs*slow
      else if (itype.eq.4) then
        ttps(lslo,ldep)=3*tts-ttp-(3*xs-xp)*slow
      else
        write(*,*) 'error in subroutine: wrong phase type'
        stop
      endif
      maxdep(lslo)=maxdep(lslo)+1
101   continue
200   continue
201   continue
c
c  differences of delay times
c
      do 120 i=1,nslow
      do 120 j=1,maxdep(i)
      ttps(i,j)=ttps(i,j)-ttps0(j)
120   continue
c
c  Interpolation for an equal-sampling reference delay-time array.
c
      do 130 i=1,nslow

      ind=1

      do 129 j=1,ntime
      t0(j)=time1+(j-1)*dtime

141   continue
      if (ttps0(ind).ge.t0(j)) then
        if (ind.eq.1) then
          ttpsmo(i,j)=ttps(i,1)
        else if (ind.gt.maxdep(i)) then
          maxt(i)=j-1
          goto 130
        else
          ind1=ind-1
          ttpsmo(i,j)=(t0(j)-ttps0(ind1))*(ttps(i,ind)-ttps(i,ind1))
          ttpsmo(i,j)=ttps(i,ind1)+ttpsmo(i,j)/(ttps0(ind)-ttps0(ind1))
        end if
          else if (ind.lt.maxdep(i)) then
            ind=ind+1
            goto 141
          else
            maxt(i)=j-1
            goto 130
      end if

129   continue
130   continue

      return
      end

C************************************************************
C
C      subroutine efa(n,dep,h,vp,vs,ro,i)
Cc
Cc--- earth flattenning approximation for gradient model
Cc--- i=1(from spherical to flat), i.ne.1 (inverse)
Cc--- modefied from kosarev's "efal" in "bip.f"
Cc
C      dimension dep(n),h(n),vp(n),vs(n),ro(n)
C      if(i.eq.1) then
C        za=0.
C        zz=0.
C        do k=1,n
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
C        end do
C      else
C        za=0.
C        zz=0.
C        do k=1,n
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
C        end do  
C      end if
C      return
C      end

!---------------------------------------------------------------------!
! version V1T1
!---------------------------------------------------------------------!
	program swf_new
cc                                                                     cc
cc Contains the subroutines of the piecewise perturbation methods      cc
cc to solve the radial Schroedinger eq. for complex wavenumber.        cc
cc Written by                                                          cc
cc                   L. Gr. Ixaru, M. Rizea,                           cc
cc  Institute of Physics and Nuclear Engineering "Horia Hulubei"       cc
cc                            Bucharest                                cc  
cc                               and                                   cc
cc                            T. Vertse                                cc
cc Institute of Nuclear Research of the Hungarian Academy of Sciences  cc
cc                             Debrecen                                cc
cc                                                                     cc
	parameter(mxrd=1000)
	implicit double complex (a-h,o-q,s-z)
	implicit real*8 (r)
	common/lj/lsp,jsp
	dimension scawf(mxrd)
	open(unit=7,file='Out/cwf.out',status='unknown')
cc----------------------------------------
cc defines the single particle hamiltonian
	call potdat
cc--------------------------------
cc reads the single particle state
	call spdat(isp,esp)
cc----------------
cc calculates w.f.
	call swf(isp,esp,scawf) !
	end
cc=======================================================================
	subroutine swf(ii,espi,scawf)
	parameter(mxrd=1000,mni=1500,mnia=300,nmax=mxrd)
	implicit double complex (a-h,o-q,s-z)
	implicit real*8 (r)
	common/lj/lsp,jsp
	common/iz/iz
	common/nstep/rstep,nstep
	common/rmax/rmax
	common/poti1/v1s,v10,v11,v12
	common/poti2/pas(mni),v(mni),v1(mni),v2(mni),intf
	common/poti3/apas(mnia),av(mnia),av1(mnia),av2(mnia),intaf
	common/rmatch/rmatch
	common/nbo/nbo
	common/par/par(3)
	common/abscis/rxeq(0:mxrd),rx(0:mxrd),numeq,numneq
	common/rc1/rc1
	dimension scawf(mxrd)
	dimension par12(2),par23(2)
	dimension y(0:nmax),yp(0:nmax)
	dimension dy1(0:nmax),dy2(0:nmax)
	external spot,pot,apot,aspot
c------------------------------------------------------------------------
cc First,the ends of the intervals i1,i2 and the value of the assymptotic
cc point must be settled. They are stored in par(i),i=1,2,3,resp.
	nbo=5
	rrxmin=nbo*rstep
	rixmin=0.d0
	par(1)=dcmplx(rrxmin,rixmin)
	rrxmax=rmax
	rixmax=0.d0
	par(2)=dcmplx(rrxmax,rixmax)
	rrxass=rmax+5.d0
	rixass=5.d0
	par(3)=dcmplx(rrxass,rixass)
cc------------------------------------------------------
cc generating potential terms and their second derivates
c	call potgen
cc--------------------------------------------------------------------
cc The program now generates the the weighted potential values v1s,v10,
cc v11 and v12 which are further necessary for the integration on i1.
	call spi1(par(1),v1s,v10,v11,v12)
cc----------------------------------------------------------------------
cc Data which are further necessary for the integration on i2,viz. 
cc Partition of i2 (vector pas), the three weighted potential values
cc (vectors v,v1 and v2) and the total number of the intervals intf 
cc are now generated. The user should provide the maximal value allowed
cc for intf. This governs the dimension of the vectors pas, v, v1 and v2.
cc In output, intf contains the number of intervals which are actually
cc used. 
	rtol=1.d-12
	if(iz.eq.2) rtol=1.d-9
	par12(1)=par(1)
	par12(2)=par(2)
	rl=lsp
        intf=mni
        call sli(rtol,pas,v,v1,v2,par12,rl,spot,pot,intf)
cc The same data are now constructed to cover the ray between par(2)
cc and par(3), i.e. in the assymptotic region. they are stored in vectors
cc apas,av,av1 and av2 and in integer intaf.
	par23(1)=par(2)
	par23(2)=par(3)
        intaf=mnia
        call sli(rtol,apas,av,av1,av2,par23,rl,aspot,apot,intaf)
cc------------------------------------------
cc numeq = number of equidistant mesh points
cc rxeq(i) = equidistant mesh points
    	numeq=nstep
	rxeq(0)=0.d0
	do int=1,numeq
	  rxeq(int)=int*rstep
	enddo
cc----------------------------------------------------
cc numneq = number of non-equidistant mesh-points
cc rx(i),i=0,1,..numneq = non-equidistant mesh-points
	numneq=intf+nbo
	rx(0)=0.d0
	do int=1,nbo
	  rx(int)=int*rstep
	enddo
        x=par(1)
        do 16 int=1,intf
	  x=x+pas(int)
	  rx(int+nbo)=x
 16	continue
cc----------------------------------------------
cc The value of the matching interval is settled
        x=par(1)
        do 8 int=1,intf
	  x=x+pas(int)
	  if(abs(x).gt.rmatch) goto 9
 8	continue
 9	match=int !index for the maching radius
cc---------------------------------------------------
cc Calculation of the scattering states for igam=1 or
cc discrete states for igam=-2,-1,0
	if(ii.eq.1)then
	  call scatt(match,ii,espi,scawf)
	else
cc The procedure to locate the eigenvalues follows. The user should first
cc give a first guess for the complex k or the energy esp
	  ck2=rc1*espi
	  ck=cwn(ii,ck2)
c	  print*,'sfindev'
	  call sfindev(ck,match,swght,gap,dwght,rtol,ierr)
	  if(ierr.ne.0) stop 'ierr.ne.0, I stop'
	  espi=ck*ck/rc1
cc We now have the eigenvalue in the output ck, and also the values of
cc three parameters (swght, gap, gwght) which are further used in the 
cc computation of the normalized eigenfunction.
cc
cc The calculation of the normalized eigenfunction follows.
cc First the eigenfunction is calculated at the original, nonequidistant
cc mesh. Later it is interpolated and written out at equidistant
cc mesh points.
cc
cc nbo is the number of equidistant mesh points in the interval I1.
cc The values of the eigenfunctions and its derivative wrt. r
cc at points rx(i), i=1,2,..,nbo,...,numneq are stored
cc in vectors y  and yp respectively.
cc
	  call snorm(nbo,ck,gap,match,swght,dwght,y,yp)
cc Write out at the original mesh
cc Interpolate and write out at equidistant mesh
	  call eqmesh(y(nbo),dy1(1),yp(nbo),yp(numneq),dy2)
cc Normalised solution follows
	  do ir=1,numeq
	    if(ir.le.nbo) then
	      scawf(ir)=y(ir)
	    else
	      scawf(ir)=dy2(ir)
	    endif
	  enddo
	endif
cc-------------------------------------------------------------
	write(7,901) lsp,jsp,espi
	write(*,901) lsp,jsp,espi
 901	format(1x,'l,2*j,e -> ',2i5,2e18.7/)
	write(7,*)' ----- Wave Functions -----'
	write(7,900)
 900	format(1x,'The wave function is in wf.dat file'/)
	open(unit=3,file='Out/wf.dat',status='unknown')
	do ir=1,nstep
	   write(3,222) ir*rstep,scawf(ir),scawf(ir)**2
	enddo
 222	format(f10.3,2(2e15.5,2x))
	close(unit=3)
	return
	end
cc=======================================================================
        subroutine snorm(nbo,ck,gap,match,swght,dwght,y,yp)
cc                                                                      cc
cc This computes the normalized eigenfunction associated with the       cc
cc eigenvalue ck                                                        cc
cc                                                                      cc
cc Input parameters:                                                    cc
cc nbo                Number of equidistant points on i1 at which the   cc
cc                    solution is needed                                cc
cc ck                 The eigenvalue, as resulted from sbr.sfindev      cc
cc gap,swght, dwght   Three double complex numbers as resulted from     cc
cc                    sbr. sfindev                                      cc
cc                                                                      cc
cc Output parameters                                                    cc
cc y                  double complex vector                             cc
cc                                                                      cc
        implicit double complex(a-h,o-q,s-z)
        implicit double precision(r)        
        parameter(mni=1500,mnia=300)
        common/poti1/v1s,v10,v11,v12
        common/poti2/pas(mni),v(mni),v1(mni),v2(mni),intf
        common/poti3/apas(mnia),av(mnia),av1(mnia),av2(mnia),intaf
	common/par/par(3)
	common/lj/lsp,jsp
        dimension y(0:*),yp(0:*)
        dimension yk(10),ykp(10)
	rl=lsp
        y(0)=0.d0
        qsd=swght-dwght*gap*gap
        sqsd=sqrt(qsd)
        swght=(1.d0,0.d0)/sqsd  
        dwght=gap/sqsd
	rh=dreal(par(1))/nbo
        call soli1(rh,rl,y(1),yp(1),yk,ykp,ck,
     *v1s,v10,v11,v12,nbo)
        do 1 n=1,nbo
	  y(n)=y(n)*swght
	  yp(n)=yp(n)*swght
	  yk(n)=yk(n)*swght
	  ykp(n)=ykp(n)*swght
 1	continue
        isens=1
        do 25 int=1,match
	  int1=int+nbo
	  y(int1)=y(int1-1)
	  yp(int1)=yp(int1-1)
	  call slix(isens,int,int,ck,pas,v,v1,v2,y(int1),yp(int1),
     *yk(nbo),ykp(nbo))
 25	continue
        call soli3(ck,dy,dyp,dyk,dykp)
        isens=-1
        int1=intf+nbo
        y(int1)=dy*dwght
        yp(int1)=dyp*dwght
        dyk=dyk*dwght
        dykp=dykp*dwght
        do 28 int=intf,match+1,isens
	  int1=int+nbo
	  y(int1-1)=y(int1)
	  yp(int1-1)=yp(int1)
	  call slix(isens,int,int,ck,pas,v,v1,v2,y(int1-1),
     *yp(int1-1),dyk,dykp)
 28	continue
        return
        end
cc=======================================================================
        subroutine soli3(ck,dy,dyp,dyk,dykp)
cc                                                                      cc     
cc This computes the outgoing solution (dy),its derivative with respect cc
cc to x (dy), the derivative of y against k (dyk) and the derivative of cc
cc dy against k (dykp) at point par(3), the l.e.end of interval i2.     cc
cc                                                                      cc
cc Input parameter                                                      cc
cc ck   double complex value for momentum k                             cc
cc                                                                      cc
cc other data are communicated through 'common' blocks                  cc
cc                                                                      cc
        implicit double complex(a-h,o-q,s-z)
        implicit double precision(r)
        parameter(mnia=300)
        common/poti3/apas(mnia),av(mnia),av1(mnia),av2(mnia),intaf
	common/par/par(3)
	common/lj/lsp,jsp
        isens=-1
        coul=aspot(par(3))
        l=lsp
        call assol(par(3),ck,l,coul,dy,dyp,dyk,dykp)
        call slix(isens,intaf,1,ck,apas,av,av1,av2,dy,dyp,dyk,dykp)
        return
        end
cc=========================================================================
cc                                                                       cc
cc	subroutine wepc                                                  cc
cc Purpose:                                                              cc
cc compute the limit of a sequence of complex numbers                    cc
cc by the epsilon algorithm of p.wynn                                    cc
cc                                                                       cc
cc Method:                                                               cc
cc the epsilon algorithm of p.wynn                                       cc
cc                                                                       cc
cc Parameters:                                                           cc
cc Input:                                                                cc
cc - s =    successive terms of the sequence to be accelerated           cc
cc - inp =  number of calls made to the routine.                         cc
cc          must be set to zero before first entry                       cc
cc          and unchanged between calls.                                 cc
cc - rtol = relative tolerance                                           cc
cc - v    = working vector of dimension maxt                             cc
cc - maxt = the maximum number of terms of the sequence                  cc
cc                                                                       cc
cc Output:                                                               cc
cc - cr=    the result (the rightmost term of the epsilon array)         cc
cc - wr=    the precedent estimation of the limit (changed after         cc
cc          each successful calculation of the rightmost element         cc
cc          from the epsilon array belonging to an odd column)           cc
cc - rerr = the estimated relative error (a zero value means             cc
cc          that this estimation cannot be obtained)                     cc
cc - ier=   error flag                                                   cc
cc          ier=0 - the tolerance was reached                            cc
cc          ier=1 - the tolerance was not reached in the frame of        cc
cc                  supplyed terms                                       cc
cc          ier=2 - two successive terms are too close each other. the   cc
cc                  calculation cannot continue. the last estimation of  cc
cc                  the limit is returned as result.                     cc
cc                                                                       cc
C============SUBROUTINE WEPC====================
        subroutine wepc(s,inp,rtol,v,maxt,cr,wr,rerr,ier)
        implicit double complex(a-h,o-q,s-z)
        implicit double precision(r)
        dimension v(maxt)
        data oner/(1.d0,0.d0)/,zero/(0.d0,0.d0)/,rsmall/1.d-7/
*********************************************************************
* reps is approximatively the smallest positive number of the machine
        data reps/1.0d-60/
*********************************************************************
        if(inp.ne.0) goto 1
        cr=s
        v(1)=s
        rerr=0.d0
        ier=1
        goto 5
1       k=0
        a1=s
        dif=a1-v(1)
        if(abs(dif).lt.reps) goto 8
        a0=oner/dif
        a2=a1
        a1=a0
        k=k+1
3       if(k.eq.inp) goto 2
        dif=a1-v(k+1)
        if(abs(dif).ge.reps) goto 7
8       ier=2
        rerr=0.d0
        goto 5
7       a0=v(k)+oner/dif
        v(k)=a2
        a2=a1
        a1=a0
        k=k+1
        goto 3
2       v(k)=a2
        v(k+1)=a1
        cr=v(k+1)
        ier=1
        if((inp/2)*2.eq.inp) goto 4
        if(inp.eq.1) goto 5
        rrvk=dreal(v(k))
        rivk=imag(v(k))
        rrvk2=dreal(v(k-2))
        rivk2=imag(v(k-2))
        rrtol=rtol
        if(abs(rrvk).gt.rsmall) rrtol=rtol*abs(rrvk)
        ritol=rtol
        if(abs(rivk).gt.rsmall) ritol=rtol*abs(rivk)
        if(abs(rrvk-rrvk2).le.rrtol.and.abs(rivk-rivk2).le.ritol) then
        ier=0
        cr=v(k)
        endif
        goto 5
4       rrcr=dreal(cr)
        ricr=imag(cr)
        rrvk1=dreal(v(k-1))
        rivk1=imag(v(k-1))
        rrtol=rtol
        if(abs(rrcr).gt.rsmall) rrtol=rtol*abs(rrcr)
        ritol=rtol
        if(abs(ricr).gt.rsmall) ritol=rtol*abs(ricr)
        if(abs(rrcr-rrvk1).le.rrtol.and.abs(ricr-rivk1).le.ritol) ier=0
5       inp=inp+1
        if(ier.eq.2) then
        cr=wr
        return
        endif
        if((inp/2)*2.eq.inp) return
        wr=cr
        return
        end
cc=======================================================================
	subroutine assol(x,pe,l,coul,u,up,uk,ukp)
cc                                                            cc
cc This sbr computes the outgoing solution o_l = g_l + i f_l  cc
cc of the coulomb equation using an asymptotic series.        cc
cc the summation of the asymptotic series is done by two      cc
cc summation methods (wynn and levin) and the most accurate   cc
cc solution is retained.                                      cc
cc                                                            cc
cc Called routines: clevin,wepc                               cc
cc                                                            cc
cc Input parameters:                                          cc
cc x = radius                                                 cc
cc pe = k-value                                               cc
cc l = angular momentum                                       cc
cc coul = 2*k*eta, where eta is the sommerfeld parameter.     cc
cc                                                            cc
cc Output parameters:                                         cc
cc u = solution                                               cc
cc up = d(u)/dr                                               cc
cc uk = d(u)/dk                                               cc
cc ukp= d(d(u)/dk)/dr                                         cc
cc                                                            cc
	implicit double precision(r)
	implicit double complex(a-h,o-q,s-z)
	parameter (max=30,max1=max-1)
        dimension p(4,2,max),s(4,2),q(4,2),sq(4,2),qq(4,2),sqq(4,2)
     *  ,qw(4,2),sqw(4,2)
	dimension vvv(max)
      	dimension arup(0:max),arlo(0:max)
	data rtol/1.d-15/,reps/1.d-12/,rsmall/1.d-7/
	alpha=0.5d0*coul
 	eta=alpha/pe
        or=dcmplx(1.d0,0.d0)
        oi=dcmplx(0.d0,1.d0)
	oz=dcmplx(0.d0,0.d0)
	nmax=max
	nm1=nmax-1
	beta=or
        z=x*pe
        f1=alpha/(pe*z)
        f2=or-f1
	f3=oi*alpha/(pe*pe)
	f4=f3*pe
        f5=-.5d0*alpha/(pe*pe)
	zlog2=log(2.d0*z)
	s(1,1)=exp(oi*(z-eta*zlog2))
	s(2,1)=oi*f2*s(1,1)
	s(3,1)=(oi*f1/z-f2*f2)*s(1,1)
	s(4,1)=(-f1*(2.d0*oi/z+3.d0*f2)/z-oi*f2*f2*f2)*s(1,1)
	fs=f3*zlog2
	s(1,2)=fs*s(1,1)
	s(2,2)=(f3/z+oi*f2*fs)*s(1,1)
	s(3,2)=(-f1*(oi/z+2.d0*f2)/pe+(oi*f1/z-f2*f2)*fs)*s(1,1)
	w1=-.5d0*oi*(f4+l+or)*(f4-l)
	w2=f5*(2.d0*f4+or)
	p(1,1,1)=w1
	p(2,1,1)=-p(1,1,1)
	p(3,1,1)=2.d0*p(1,1,1)
	p(4,1,1)=-6.d0*p(1,1,1)
	p(1,2,1)=w2
	p(2,2,1)=-p(1,2,1)
        p(3,2,1)=2.d0*p(1,2,1)
	do 1 n=2,nmax
	  w1=-.5d0*oi*(f4+l+n)*(f4-l+n-or)/float(n)
	  w2=f5*(2.d0*(f4+n)-or)/float(n)
	  p(1,1,n)=w1*p(1,1,n-1)
	  p(2,1,n)=-n*p(1,1,n)
	  p(3,1,n)=-(n+or)*p(2,1,n)
	  p(4,1,n)=-(n+2.d0*or)*p(3,1,n)
	  p(1,2,n)=w1*p(1,2,n-1)+w2*p(1,1,n-1)
	  p(2,2,n)=-n*p(1,2,n)
	  p(3,2,n)=-(n+or)*p(2,2,n)
 1	continue
	zi=or/z
	zc=zi
	do 11 n=1,nmax
	  do 10 i=1,4
c	    do 10 j=1,2
	    do 101 j=1,2
c	      if(i+j.gt.5) goto 10
	      if(i+j.gt.5) goto 101
	      p(i,j,n)=p(i,j,n)*zc
 101	  continue
 10	  continue
	  zc=zc*zi
 11	continue
	do 20 i=1,4
c	  do 20 j=1,2
	  do 201 j=1,2
c	    if(i+j.gt.5) goto 20
	    if(i+j.gt.5) goto 201
	    inp=0
	    ss=oz
	    do 31 n=1,nm1
	      term=p(i,j,n)
	      ss=ss+term
 31	    continue
	    sum=oz
	    do 30 n=1,nm1
	      term=p(i,j,n)
	      sum=sum+term
	      sofn=sum
	      zofn=p(i,j,n+1)
	      if(abs(zofn).le.rtol) goto 21
c  convergence acceleration by levin	
	      call clevin(sofn,zofn,beta,inp,arup,arlo,max,estlim)
	      if(inp.eq.0) then
		estlio=estlim
		goto 23
	      endif
	      rres=dreal(estlim)
	      ries=imag(estlim)
	      rreso=dreal(estlio)
	      rieso=imag(estlio)
	      rreps=reps
	      if(abs(rres).gt.rsmall) rreps=reps*abs(rres)
	      rieps=reps
	      if(abs(ries).gt.rsmall) rieps=reps*abs(ries)
        if(abs(rres-rreso).le.rreps.and.abs(ries-rieso).le.rieps) 
     *  goto 21     
	      estlio=estlim
c  epsilon algorithm by wynn
 23	      call wepc(sum,inp,reps,vvv,max,cr,wr,rerr,ier)
	      if(ier.ne.1) goto 21
 30	    continue
c	write(*,*) 'tol is not achieved'
 21	    continue 
	    q(i,j)=wr
	    qq(i,j)=ss
	    qw(i,j)=estlim
 201	continue
 20	continue
        zc=zi
	do 40 i=2,4
	  do 41 j=1,2
	    if(i+j.gt.5) goto 41
	    q(i,j)=q(i,j)*zc
	    qq(i,j)=qq(i,j)*zc
	    qw(i,j)=qw(i,j)*zc
 41	  continue
	  zc=zc*zi
 40	continue
	q(1,1)=q(1,1)+or
	qq(1,1)=qq(1,1)+or
	qw(1,1)=qw(1,1)+or
	do 4 i=1,4
	  do 5 j=1,2
	    if(i+j.gt.5) goto 5
	    sq(i,j)=dcmplx(0.d0,0.d0)
	    sqq(i,j)=dcmplx(0.d0,0.d0)
	    sqw(i,j)=dcmplx(0.d0,0.d0)
	    do 6 k=1,i
c	      do 6 m=1,j
	      do 61 m=1,j
c		if(k+m.gt.5) goto 6
		if(k+m.gt.5) goto 61
		do 7 l1=1,i
c		  do 7 n=1,j
		  do 71 n=1,j
c		    if(l1+n.gt.5) goto 7
		    if(l1+n.gt.5) goto 71
c		    if(k+l1-i.ne.1) goto 7
		    if(k+l1-i.ne.1) goto 71
c		    if(m+n-j.ne.1) goto 7
		    if(m+n-j.ne.1) goto 71
		    factor=or
		    if(i-k.eq.1.or.i-l1.eq.1) factor=(i-1)*or
		    sq(i,j)=sq(i,j)+factor*s(k,m)*q(l1,n)
		    sqq(i,j)=sqq(i,j)+factor*s(k,m)*qq(l1,n)
		    sqw(i,j)=sqw(i,j)+factor*s(k,m)*qw(l1,n)
 71		  continue
 7		  continue
 61	    continue
 6	    continue
 5	  continue
 4	continue
	u=sq(1,1)
	up=pe*sq(2,1)
        upp=pe*pe*sq(3,1)
	v=sqq(1,1)
	vp=pe*sqq(2,1)
        vpp=pe*pe*sqq(3,1)
	w=sqw(1,1)
	wp=pe*sqw(2,1)
        wpp=pe*pe*sqw(3,1)
	uk=sq(1,2)+x*sq(2,1)
	ukp=sq(2,1)+pe*sq(2,2)+z*sq(3,1)
c	ukpp=pe*(2.d0*sq(3,1)+pe*sq(3,2)+z*sq(4,1))
	vk=sqq(1,2)+x*sqq(2,1)
	vkp=sqq(2,1)+pe*sqq(2,2)+z*sqq(3,1)
c	vkpp=pe*(2.d0*sqq(3,1)+pe*sqq(3,2)+z*sqq(4,1))
	wk=sqw(1,2)+x*sqw(2,1)
	wkp=sqw(2,1)+pe*sqw(2,2)+z*sqw(3,1)
c	wkpp=pe*(2.d0*sqw(3,1)+pe*sqw(3,2)+z*sqw(4,1))
	pot=(l*(l+1.d0)/x+2.d0*eta*pe)/x-pe*pe
	dpotk=-2.d0*pe
	du=upp-pot*u
	dv=vpp-pot*v
	dw=wpp-pot*w
	rdv=abs(dv)
	rdu=abs(du)
	rdw=abs(dw)
	rdmin=min(rdv,rdu,rdw)
	if(rdmin.eq.rdu) goto 25
	if(rdmin.eq.rdv) then
	  u=v
	  up=vp
	  uk=vk
	  ukp=vkp
	else
	  u=w
	  up=wp
	  uk=wk
	  ukp=wkp
	endif
 25	quady2=.5d0*(ukp*u-up*uk)/pe
c
c	quad2 is the integral from infinity to x of u*u; u is the
c	outgoing coulomb function
c
	return
	end
cc========================================================================
cc   subroutine clevin(sofn,rofn,beta,n,arup,arlo,larray,estlim)        cc
cc                                                                      cc
cc                                                                      cc
cc description :                                                        cc
cc -------------                                                        cc
cc                                                                      cc
cc subroutine clevin computes levin's transformation.                   cc
cc (d. levin (j. comput. math., vol. b3, 1973, 371-388 ) using the      cc
cc partial sums s(n) and terms a(n) of a given series as input data.    cc
cc                                                                      cc
cc here, the definition of levin's sequence transformation according to cc
cc eq. (7.1-7) of e. j. weniger, comput. phys. rep., vol. 10 (1989),    cc
cc 189 - 371, is used.                                                  cc
cc                                                                      cc
cc in levin's sequence transformation it is assumed that the partial    cc
cc sums s(n) of the series to be transformed can be written as          cc
cc                                                                      cc
cc s(n) = s + r(n),                                                     cc
cc                                                                      cc
cc with s being the limit or antilimit of this series, and that the     cc
cc exact remainder r(n) can be approximated by a remainder estimate     cc
cc omega(n) in such a way that the ratio                                cc
cc                                                                      cc
cc [ s(n) - s ] / omega(n)                                              cc
cc                                                                      cc
cc can be expressed as a by a poincare-type asymptotic expansion in     cc
cc inverse powers of (n + beta).                                        cc
cc                                                                      cc
cc levin suggested some simple remainder estimates omega(n) which can   cc
cc be computed be determined from at most two terms a(n) of the series  cc
cc to be transformed. with the help of these remainder estimates the    cc
cc following variants of levin's sequence transformations result        cc
cc (compare section 7.3 of e. j. weniger, comput. phys. rep.,           cc
cc  vol. 10 (1989), 189):                                               cc
cc                                                                      cc
cc (n+beta) * a(n)                     : u transformation, eq. (7.3-5)  cc
cc a(n)                                : t transformation, eq. (7.3-7)  cc
cc a(n+1)                              : d transformation, eq. (7.3-7)  cc
cc a(n + 1) * a(n) / [a(n) - a(n + 1)] : v transformation, eq. (7.3-11) cc
cc                                                                      cc
cc it is always assumed that the index n of the terms and the partial   cc
cc sums s(n) satisfies n >= 0.                                          cc
cc                                                                      cc
cc both numerator and denominator of this transformation are calculated cc
cc with the help of a simple extension of a 2-dimensional nonlinear     cc
cc 3-term recurrence formula which was originally derived by            cc
cc t.fessler, w.f.ford and w.f.smith                                    cc
cc ( acm trans.math.software ,vol. 9 ( 1983 ), pp.346-354 ).            cc
cc                                                                      cc
cc in this subroutine a variant of the so-called moving lozenge tech-   cc
cc nique ( p.wynn, r.f.t.i.-chiffres, vol.9, pp.327-362 ( 1966 )) is    cc
cc used, i.e., only one counterdiagonal of the numerator and denomina-  cc
cc tor tables have to be stored at a time. consequently, only two       cc
cc 1-dimensional arrays of appropriate length are needed in this        cc
cc subroutine.                                                          cc
cc in each call this subroutine calculates new counterdiagonals of the  cc
cc numerator and denominator tables and overwrites the previous entries.cc
cc in addition, a new estimate of the limit is calculated which is the  cc
cc ratio of the array elements arup(0) and arlo(0). hence, for every    cc
cc new sequence element s-of-n this subroutine has to be called again.  cc
cc                                                                      cc
cc input parameters :                                                   cc
cc ----------------                                                     cc
cc                                                                      cc
cc sofn   : element s(n) of the sequence to be accelerated              cc
cc                                                                      cc
cc rofn   : estimate of the leading term of the poincare-type asympto-  cc
cc          tic expansion of the remainder r(n) in inverse powers of    cc
cc          (n+beta), i.e. omega(n).                                    cc
cc          it is tacitly assumed that rofn is different from zero.     cc
cc                                                                      cc
cc beta   : parameter contained in the inverse powers                   cc
cc          ( n + beta ) ** ( -j )                                      cc
cc          which describe the higher order contributions in the        cc
cc          poincare-type expansion of the remainder r(n).              cc
cc                                                                      cc
cc n      : number of the last sequence element (n has to satisfy the   cc
cc          inequality n .ge. 0 )                                       cc
cc                                                                      cc
cc arup   : 1-dimensional array to store the actual counterdiagonal of  cc
cc          the numerator table.                                        cc
cc                                                                      cc
cc arlo   : 1-dimensional array to store the actual counterdiagonal of  cc
cc          the denominator table.                                      cc
cc                                                                      cc
cc larray : length of the 1-dimensional arrays arup and arlo.           cc
cc                                                                      cc
cc output parameter :                                                   cc
cc ----------------                                                     cc
cc                                                                      cc
cc estlim : approximation to the limit of the sequence s(n) which is    cc
cc          to be accelerated.                                          cc
cc          we have here estlim = arup(0)/arlo(0).                      cc
cc                                                                      cc
cc machine-dependent parameters:                                        cc
cc ----------------------------                                         cc
cc                                                                      cc
cc huge,tiny : huge should be set close to but not identical with the   cc
cc             largest floating point number representable on the com-  cc
cc             puter, and tiny should be set close to but not identical cc
cc             with the smallest floating point number representable    cc
cc             on the computer.                                         cc
cc                                                                      cc
cc                                                                      cc
	subroutine clevin(sofn,rofn,beta,n,arup,arlo,larray,estlim)
	implicit double complex(a-h,o-z)
	double precision huge, tiny, one
	dimension arup(0:larray),arlo(0:larray)
	parameter ( huge = 1.d+60 , tiny = 1.d-60 , one = 1.d0 )
	zone=dcmplx(one,0.d0)
c
c  check whether n is less than zero
c
	if (n.lt.0) then
	  write(*,1000)
 1000	  format('*** error in subroutine clevin ***'/
     1         '=================================')
	  write(6,1010) n
 1010	  format('illegal input parameter n which must be .ge. zero'/
     1         'here we have n = ',i4)
	  stop
c
c  check whether the arrays arlo and arup are large enough to store
c  all the elements required
c
	else if (n.gt.larray) then
	  write(*,1000)
	  write(*,1020) n,larray
 1020	  format('illegal input parameters n and larray'/
     1  'we must have n .le. larray'/
     2  'here we have n = ',i4,' and larray = ',i4)
	  return
	end if
c
c  this is reached if the general levin transformation has to be
c  calculated.
c
c  first it has to be checked whether n .ge. 1 holds, i.e., whether we
c  proceed with the calculation of the next counter diagonals of the
c  two arrays or whether we have to initialize a new calculation.
c  this is done if n = 0 holds.
c
      arup(n)   = sofn / rofn
      arlo(n)   = zone  / rofn
c
c  check whether n .gt. 0 holds. in that case the recursive computation
c  of the numerators and denominators of the general levin transformatio
c  is to be continued.
c
      if ( n .gt. 0 ) then
        arup(n-1) = arup(n) - arup(n-1)
        arlo(n-1) = arlo(n) - arlo(n-1)
c
c  check whether n .gt. 1 holds. in that case the actual recursion
c  of the numerators and denominators starts.
c
        if (n.gt.1) then
          bn1     = beta + float(n-1)
          bn2     = beta + float(n)
          coef    = bn1 / bn2
          do 10 j = 2,n
            fact  = (beta+float(n-j)) * coef**(j-2) / bn2
c
c  perform the calculation of the new elements of the two tables.
c
            arup(n-j)  = arup(n-j+1) - fact*arup(n-j)
            arlo(n-j)  = arlo(n-j+1) - fact*arlo(n-j)
10        continue
        end if
      end if
c
c  compute the estimate for the limit of the sequence.
c
        if (abs(arlo(0)).lt.tiny) then
          estlim = dcmplx(huge,huge)
        else
          estlim = arup(0)/arlo(0)
        end if
      return
      end
cc=======================================================================
        subroutine sfindev(ck,match,swght,gap,dwght,rtol,ierr)
cc                                                                      cc
cc This computes the eigenvalue of the boundary value problem.          cc
cc Newton method is used in the iteration                               cc
cc                                                                      cc
cc Input parameter                                                      cc
cc match   The interval at the r.h.end of which the solutions           cc
cc         should be matched                                            cc
cc                                                                      cc
cc Input/output parameter:                                              cc
cc ck      In input it furnishes a guess value. In output, it contains  cc
cc         the resultant eigenvalue.                                    cc
cc                                                                      cc
cc Output parameters:                                                   cc
cc swght,gap,dgwht   Double complex numbers which are further on used   cc
cc                   for the normalization of the eigenfunction         cc
cc                                                                      cc
        implicit double complex(a-h,o-q,s-z)
        implicit double precision(r)        
        parameter(mni=1500,mnia=300)
	parameter(maxnewit=20)
        common/poti1/v1s,v10,v11,v12
        common/poti2/pas(mni),v(mni),v1(mni),v2(mni),intf
        common/poti3/apas(mnia),av(mnia),av1(mnia),av2(mnia),intaf
	common/par/par(3)
	common/lj/lsp,jsp
	dimension y(1),yp(1),yk(1),ykp(1)
        newit=0
	rl=dfloat(lsp)
c	print*,'soli1'
 1	call soli1(dreal(par(1)),rl,y,yp,yk,ykp,ck,v1s,v10,v11,
     *v12,1)    
	sy=y(1)
	syp=yp(1)
c	print*,' sy',sy
c	print*,' syp',syp
	syk=yk(1)
	sykp=ykp(1)
        isens=1
        call slix(isens,1,match,ck,pas,v,v1,v2,sy,syp,syk,sykp)
        call soli3(ck,dy,dyp,dyk,dykp)
        isens=-1
        if(intf.ge.match+1)call slix(isens,
     *intf,match+1,ck,pas,v,v1,v2,dy,dyp,dyk,dykp)
        phi=sy*dyp-dy*syp
        dphi=syk*dyp+sy*dykp-dyk*syp-dy*sykp
        dphi=(dphi*sy*dy-phi*(syk*dy+sy*dyk))/(sy*dy)
        dk=-phi/dphi
        ck=ck+dk
        newit=newit+1
	if(abs(dk/ck).le.rtol) then
	  ierr=0
	  goto 2
	else
	  if(newit.le.maxnewit) goto 1
	  ierr=1
	  goto 2
	endif
 2	swght=.5d0*(syp*syk-sykp*sy)/ck
        dwght=.5d0*(dyp*dyk-dykp*dy)/ck
        gap=sy/dy
        return
        end
cc=======================================================================
	subroutine scatt(match,ii,espi,scawf)
	parameter(mxrd=1000,nmax=mxrd)
	implicit double precision(r)
	implicit double complex (a-h,o-q,s-z)
	common/rc1/rc1
	common/nbo/nbo
	common/abscis/rxeq(0:mxrd),rx(0:mxrd),numeq,numneq
	dimension sy(0:nmax),syp(0:nmax)
	dimension dy1(0:nmax),dy2(0:nmax),dyp2(0:nmax),dyp1(0:nmax)
	dimension scawf(mxrd)
        ck2=rc1*espi
        ck=cwn(ii,ck2)
cc Complex k is in ck and calculation of the normalized scattering
cc solution by calling subroutine scatwf as follows. 
cc You will have the wf and its derivative in vectors sy and syp. 
cc The vectors dy1, dyp1, dy2 and dyp2 play here the role of
cc working vectors.
        call scatwf(nbo,ck,match,sy,syp,dy1,dyp1,dy2,dyp2)
cc
cc Now the wavefunction at equidistant partition is generated.
cc
        call eqmesh(sy(nbo),dy1(1),syp(nbo),syp(numneq),dy2)
        do ir=1,numeq
	  if(ir.le.nbo) then
	    scawf(ir)=sy(ir)
	  else
	    scawf(ir)=dy2(ir)
	 endif
	enddo
cc
cc Set correct phase of the radial scattering wave function
cc
        cphase=(1.d0,0.d0)
        if(dreal(scawf(1)).lt.0.d0)cphase=(-1.d0,0.d0)
        do k=1,numeq
	  scawf(k)=scawf(k)*cphase
        enddo
	return
	end
cc=======================================================================
        subroutine eqmesh(y,b,ypmin,ypmax,sol)
	parameter(mxrd=1000)
        implicit double complex (a-h,o-q,s-z)
        implicit double precision (r)
	dimension y(0:*),b(1),sol(0:*)
	common/abscis/rxeq(0:mxrd),rx(0:mxrd),numeq,numneq
	common/nbo/nbo
        call splinc(rx(nbo),y(0),numneq-nbo,ypmin,ypmax,b)
	do 1 int=nbo+1,numeq
	  call spintc(rx(nbo),y(0),b,numneq-nbo,rxeq(int),sol(int))
 1	continue
	return
	end
cc=======================================================================
	subroutine spintc(rxa,ya,y2a,n,rx,y)
c
c  spline interpolation using the spline coefficients in y2a      
c
	implicit double complex (a-h,o-q,s-z)
	implicit double precision (r)
	dimension rxa(n),ya(n),y2a(n)
	klo=1
	khi=n
 1	if (khi-klo.gt.1) then
	  k=(khi+klo)/2
	  if(rxa(k).gt.rx)then
	    khi=k
	  else
	    klo=k
	  endif
	  goto 1
	endif
	rh=rxa(khi)-rxa(klo)
	if (rh.eq.0.d0) stop 'bad xa input.'
	ra=(rxa(khi)-rx)/rh
	rb=(rx-rxa(klo))/rh
	y=ra*ya(klo)+rb*ya(khi)+
     *((ra**3-ra)*y2a(klo)+(rb**3-rb)*y2a(khi))*(rh**2)/6.d0
	return
	end
cc=======================================================================
        subroutine scatwf(nbo,ck,match,sy,syp,dy1,dyp1,dy2,dyp2)
	parameter(mni=1500)
        implicit double complex(a-h,o-q,s-z)
        implicit double precision(r)        
        dimension cfc(0:3),cgc(0:3),cgcp(0:3),cfcp(0:3),sig(0:3)
	dimension sy(0:*),syp(0:*),dy1(0:*),dyp1(0:*),dy2(0:*),
     *  dyp2(0:*)
	dimension syk(mni),sykp(mni),x1(2)
        common/poti1/v1s,v10,v11,v12
        common/poti2/pas(mni),v(mni),v1(mni),v2(mni),intf
	common/lj/lsp,jsp
	common/par/par(3)
	common/rc1/rc1
	common/retak/retak
	rpi=4.d0*atan(1.d0)
	ci=(0.d0,1.d0)
	rh=dreal(par(1))/nbo
	rl=lsp
	call soli1(rh,rl,sy(1),syp(1),syk,sykp,ck,v1s,
     * v10,v11,v12,nbo)
        isens=1
cc Integrates outwards and store in sy
	do 1 int=1,match
	  int1=int+nbo
	  sy(int1)=sy(int1-1)
	  syp(int1)=syp(int1-1)
	  call slix(isens,int,int,ck,pas,v,v1,v2,sy(int1),syp(int1),
     *  syk(int1),sykp(int1))
 1	continue
	crho=ck*par(2)
	eta=retak/ck
	mode1=1
	kfn=0
	zlp=rl
	nlp=1
	ifail=1
	iprint=1
c
c	calculating coulomb functions at complex arguments
c
	call wclbes(crho,eta,zlp,nlp,cfc,cgc,cfcp,cgcp,sig,
     *  kfn,mode1,ifail,iprint)
c O=G+i*F and O'
        out=cgc(0)+ci*cfc(0)
        outp=cgcp(0)+ci*cfcp(0)
c dF/dr & dO/dr
	cfcp(0)=ck*cfcp(0)
        outp=outp*ck
	intf1=nbo+intf
c integrate inward F(r) and store in dy1
	dy1(intf1)=cfc(0)
	dyp1(intf1)=cfcp(0)
	dyk=(0.d0,0.d0)
	dykp=(0.d0,0.d0)
	isens=-1
	do 2 int=intf,match+1,-1
	  intf1=nbo+int
	  dy1(intf1-1)=dy1(intf1)
	  dyp1(intf1-1)=dyp1(intf1)
	  if(intf.ge.match+1)call slix(isens,
     *int,int,ck,pas,v,v1,v2,dy1(intf1-1),dyp1(intf1-1),dyk,dykp)
 2	continue
c dG/dr
c integrate inward G(r) and store in dy2
	cgcp(0)=ck*cgcp(0)
	intf1=nbo+intf
	dy2(intf1)=cgc(0)
	dyp2(intf1)=cgcp(0)
	dyk=(0.d0,0.d0)
	dykp=(0.d0,0.d0)
	do 3 int=intf,match+1,-1
	  intf1=nbo+int
	  dy2(intf1-1)=dy2(intf1)
	  dyp2(intf1-1)=dyp2(intf1)
	  if(intf.ge.match+1)call slix(isens,
     *int,int,ck,pas,v,v1,v2,dy2(intf1-1),dyp2(intf1-1),dyk,dykp)
 3	continue
	intf1=intf1-1
	x1(1)=(sy(intf1)*dyp2(intf1)/dy2(intf1)-syp(intf1))/
     +        (dy1(intf1)*dyp2(intf1)/dy2(intf1)-dyp1(intf1))
	x1(2)=(sy(intf1)*dyp1(intf1)/dy1(intf1)-syp(intf1))/
     +        (dy2(intf1)*dyp1(intf1)/dy1(intf1)-dyp2(intf1))
	a_u=0.5d0*(x1(2)+ci*x1(1))
	b_u=0.5d0*(x1(2)-ci*x1(1))
c below is the S-matrix of H_0
	smatrix=-b_u/a_u
c below is the corresponding phase shift of H_0
	delta=(0.d0,-0.5d0)*log(smatrix)
	om=(0.d0,-1.d0)*delta
	tmp1=exp(om)
	renfac=rc1
	ecm=ck*ck/renfac
	xn=sqrt(ck/(rpi*ecm))*(0.d0,0.5d0)
	coef=xn/a_u*tmp1
	nbmat=nbo+match
	nbintf=nbo+intf
	do 4 int=nbmat,nbintf
	  sy(int)=x1(1)*dy1(int)+x1(2)*dy2(int)
	  syp(int)=x1(1)*dyp1(int)+x1(2)*dyp2(int)
 4	continue
	sy(0)=(0.d0,0.d0)
	do 5 int=1,nbintf
	  sy(int)=sy(int)*coef
	  syp(int)=syp(int)*coef
 5	continue
	return
	end
cc==========================================================================
	subroutine soli1(rh,rl,u,up,ue,uep,pe,vc,v0,v1,v2,nbo)
cc                                                                        cc
cc This computes the regular solution at a series of points on i1.        cc
cc                                                                        cc
cc Input parameters                                                       cc
cc rl          Double precision number for the angular momentum l         cc
cc vc,v0,v1,v2 Double complex numbers as resulted from sbr.spi1           cc
cc pe          Double complex number representing the momentum k          cc
cc rh          Double precision number for the r.h.end of i1              cc
cc nbo         Number of equidistant points on i1 at which the solution   cc
cc             is required.                                               cc
cc                                                                        cc
cc Output parameters                                                      cc
cc u           Double complex vector. in u(i) the value of the solution   cc
cc             at point x_i=i*rh/nbo is found. i ranges from 1 to nbo     cc
cc up          The same for the first derivative of u                     cc
cc ue          The same for the first derivative of u against k           cc
cc uep         The same for the first derivative of up against k          cc
cc                                                                        cc
	implicit double precision(r)
	implicit double complex(a-h,o-q,s-z)
	dimension u(1),up(1),ue(1),uep(1)
	dimension f(4),ca(65),cap(65),a(65),ap(65),cfu(65),cfue(65)
	r=nbo*rh
	e=pe*pe
	l=rl
	rlmem=rl
	f(1)=r*vc
	f(2)=(v0-pe*pe)*r*r
	fp=-2.d0*pe*r*r
	f(3)=v1*r*r*r
	f(4)=v2*r*r*r*r
	cz=dcmplx(0.d0,0.d0)
	cu=dcmplx(1.d0,0.d0)
	iqd=0
	iqu=0
	iql=0
	do 20 i=1,65
	 cfu(i)=cz
	 cfue(i)=cz
	 a(i)=cz
c 20	ap(i)=cz
		ap(i)=cz
  20	continue
	a(2)=cu
	cfu(2)=cu
	cpot=l*(l+1)/(r*r)+vc/r+(v0-pe*pe)+v1*r+v2*r*r
	do 10 kpert=1,15
	 do 1 iq=iqd-1,iqu+2
	  iqf=iq+3
	  ca(iqf)=cz
	  cap(iqf)=cz
	  ipu=2
	  ipl=iq-iqu
	  if(ipl.le.-1) ipl=-1
	  ipl=-1
	  ipu=2
	  if(ipu.ge.iq) ipu=iq
	  do 2 ip=ipl,ipu
	   ca(iqf)=ca(iqf)+f(ip+2)*a(iq-ip+2)
	   cap(iqf)=cap(iqf)+f(ip+2)*ap(iq-ip+2)
 2	  continue
	  cap(iqf)=cap(iqf)+fp*a(iq+2)
 1	 continue
	 do 3 iq=1,iqu+4
	  ap(iq)=cz
c 3	 a(iq)=cz
	 a(iq)=cz
  3	 continue
	 iql=iql+1
	 iqu=iqu+4
	 do 4 iq=iql,iqu
	  den=(iq+l+1)*(iq+l)-l*(l+1)
	  a(iq+2)=ca(iq+1)/den
	  ap(iq+2)=cap(iq+1)/den
	  cfu(iq+2)=cfu(iq+2)+a(iq+2)
	  cfue(iq+2)=cfue(iq+2)+ap(iq+2)               
 4	 continue
	 y=cu
	 yp=(l+1)*cu
	 yp=yp/r
	 ye=cz
	 yep=cz
	 ypp=(l+1)*l*cu
	 ypp=ypp/(r*r)
	 yepp=cz
	 do 11 iq=1,iqu
	  iqc=iq+2
	  su=(iq+l+1)*cfu(iqc)/r
	  sue=(iq+l+1)*cfue(iqc)/r
	  y=y+cfu(iqc)
	  yp=yp+su
	  ypp=ypp+(iq+l)*su/r
	  ye=ye+cfue(iqc)
	  yep=yep+sue
	  yepp=yepp+(iq+l)*sue/r
 11	 continue
	 rtol=1.d-15
	 crty=cpot*y
	 crtye=-2.d0*pe*y+cpot*ye
	 dy=ypp-crty
	 dye=yepp-crtye
	 if(abs(dy).le.rtol*abs(ypp).and.abs(dye).le.rtol*abs
     1(yepp)) goto 15
 10	continue
 15	continue            
	do 16 n=1,nbo
	 us=cz
	 ups=cz
	 upps=cz
	 ues=cz
	 ueps=cz
	 uepps=cz
	 rl=n*rh/r
c	 print*,'rl',rl
	 do 19 in=1,iqu
	  iq=iqu-in+1
	  iqc=iq+2
	  su=(iq+l+1)*cfu(iqc)
	  sue=(iq+l+1)*cfue(iqc)
	  us=rl*(cfu(iqc)+us)
c	  print*,'us',us
	  ups=rl*(su+ups)
	  upps=rl*((iq+l)*su+upps)
	  ues=rl*(cfue(iqc)+ues)
	  ueps=rl*(sue+ueps)
	  uepps=rl*((iq+l)*sue+uepps)
 19	 continue
	 u(n)=cu+us
c	 print*,'n,u(n)=cu+us',n,u(n)
	 up(n)=(l+1)*cu+ups
	 upp=l*(l+1)*cu+upps
	 ue(n)=ues
	 uep(n)=ueps
	 uepp=uepps
	 ril=1.d0
	 do 21 i=1,l+1
	  ril=ril*rl
 21	 continue
	 rl=n*rh
	 u(n)=u(n)*ril
c	 print*,'n,ril,u(n)',n,ril,u(n)
c	 print*,'n,u(n)',n,u(n)
	 up(n)=up(n)*ril/rl
c	 print*,'n,up(n)',n,up(n)
c	 print*
	 upp=upp*ril/(rl*rl)
	 ue(n)=ue(n)*ril
	 uep(n)=uep(n)*ril/rl
	 uepp=uepp*ril/(rl*rl)
	 cpot=l*(l+1)/(rl*rl)+vc/rl+v0-e+v1*rl+v2*rl*rl
	 crty=u(n)*cpot
	 dy=upp-crty
	 crtye=-2.d0*pe*u(n)+cpot*ue(n)
	 dye=uepp-crtye
 16	continue                  
	rl=rlmem
	return
	end
cc=============================================================================
        subroutine slix(isens,inti,intf,ck,pas,v,v1,v2,y,yp,yk,ykp)
cc                                                                           cc
cc This sbr. propagates the solution along a prescribed interval             cc
cc                                                                           cc
cc Input parameters                                                          cc
cc isens        Sense of the propagation; isens=1 for forward while          cc
cc              isens=-1 for backward propagation.                           cc
cc inti         The label of the interval from which the propagation begins. cc
cc intf         The label of the last interval to be covered                 cc
cc              note: inti<intf requires isens=1,and so the solution is      cc
cc              propagated forwards from the l.h. end of inti up to the      cc
cc              r.h.end of intf;                                             cc
cc              inti>intf requires isens=-1; the solution is then propagated cc
cc              backwards from the r.h.end of inti up to the l.h.end of intf.cc
cc              One can also consider the propagation on a single interval.  cc
cc              Then inti=intf and, if isens is settled 1, the propagation   cc
cc ck           The double complex value of the momentum k                   cc
cc pas,v,v1,v2  Double precision vect. with the values resulted from sbr sli cc
cc                                                                           cc
cc input/output parameters                                                   cc
cc y,yp,yk,ykp  Four double complex parameters. In input, they furnish the   cc
cc              initial values for the solution (y), its first derivative    cc
cc              against x (yp), its first derivative against k (yk)          cc
cc              and the first derivative of yp against k (ykp).              cc
cc              In output these parameters contain the value of the same     cc
cc              quantities at the other end of the domain.                   cc
cc                                                                           cc
        implicit double precision(r)
        implicit double complex(a-h,o-q,s-z)
        dimension pas(1),v(1),v1(1),v2(1)
        e=ck*ck
        pe=ck
        y0=y
        yp0=yp
        ye0=yk
        yep0=ykp
	do 2 n=inti,intf,isens
	  h=pas(n)
	  eh=e*h*h
	  call tmcpm2(eh,v(n),v1(n),v2(n),csi,eta0,eta1,eta2,eta3,
     1eta4,eta5,eta6,t11,t12,t21,t22)
	  if(isens.lt.0) goto 4
	  y1=t11*y0+t12*yp0*h
	  yp1=t21*y0/h+t22*yp0
	  ye1=t11*ye0+t12*yep0*h-pe*h*h*(y0*eta0+yp0*h*eta1)
	  yep1=t21*ye0/h+t22*yep0-pe*h*(y0*(csi+eta0)+yp0*h*eta0)
	  goto 5
 4	  y1=t22*y0-t12*yp0*h
	  yp1=-t21*y0/h+t11*yp0   
	  ye1=t22*ye0-t12*yep0*h-pe*h*h*(y0*eta0-yp0*h*eta1)
	  yep1=-t21*ye0/h+t11*yep0+pe*h*(y0*(csi+eta0)-yp0*h*eta0)
 5	  y0=y1
	  yp0=yp1
	  ye0=ye1
	  yep0=yep1
 2	continue
        y=y0
        yp=yp0
        yk=ye0
        ykp=yep0
        return
        end
cc=======================================================================
	subroutine tmcpm2(eh,v,v1,v2,csi,eta0,eta1,eta2,eta3,
     1eta4,eta5,eta6,t11,t12,t21,t22)
c
c       calculation of the elements t11,..,t22 of the transfer
c         matrix t which is used to
c       advance the solution and its derivative from a point
c       r to r+h i.e. outward
c
c     (   y(r+h)  )      ( t11     t12*h  )   ( y(r) )
c     (           )   =  (                ) * (      )
c     (   y'(r+h) )      ( t21/h   t22    )   ( y'(r))
c
c       or from r to r-h i.e. inward
c
c     (   y(r-h)  )      ( t22     -t12*h  )   ( y(r) )
c     (           )   =  (                 ) * (      )
c     (   y'(r-h) )      ( t21/h   t11     )   ( y'(r))
c
	implicit double precision(r)
	implicit double complex(a-h,o-q,s-z)
        z=v-eh
        call cgebas(z,csi,eta0,eta1,eta2,eta3,eta4,eta5,eta6)
c
c       corrections of the 1-st order perturbation theory
c
        corr1=v1*eta1
        corr2=v2*eta2
        v11=v1*v1/6.d0
        v12=2.d0*v1*v2
        v22=.3d0*v2*v2
        auxv=v11+v22/3.d0
        auxc=auxv*eta2
c
c       corrections of the 2-nd order perturbation theory
c
        corr11=-auxc+(v12-v22)*eta3
        corr21=-auxv*eta1-(7.d0*v11+2.d0*v22)*eta2+19.d0*v22*eta3
        corr12=-auxv*eta3+3.d0*v22*eta4
        corr22=-auxc-(v12+v22)*eta3
c
c       calculation of the transfer matrix by adding corrections 
c       of the order 0,1 and 2 together
c
        t11=csi-corr1+corr11
        t21=z*(eta0+corr2)+corr21
        t12=eta0-corr2+corr12
        t22=csi+corr1+corr22
        return
        end
cc=======================================================================
      subroutine cgebas(z,csi,eta0,eta1,eta2,eta3,eta4,eta5,eta6)
c
c      subroutine to calculate the quantities csi,..,eta6 defined in p.198
c      of ref.1.  by eq. (3.4.15)-(3.4.18b)
c
	implicit double precision(r)
	implicit double complex(a-h,o-q,s-z)
	rtz=1.d0
	ru=abs(z)
	if(ru.lt.rtz) goto 72
	squ=sqrt(z)
	squ=squ*dcmplx(0.d0,1.d0)
	eta0=sin(squ)/squ
	csi=cos(squ)
 71	eta1=(csi-eta0)/z
	eta2=(eta0-3.d0*eta1)/z
	eta3=(eta1-5.d0*eta2)/z
	eta4=(eta2-7.d0*eta3)/z
	eta5=(eta3-9.d0*eta4)/z
	eta6=(eta4-11.d0*eta5)/z
	return
 72	eta6=(269583552.d4+z*(89861184.d3+z*(1321488.d3+z*(11592.d3+
     1z*(69.d3+z*(3.d2+z))))))/3643017329952.d5
	eta5=(14018344704.d2+z*(539167104.d2+z*(89861184.d1+z*(
     1880992.d1+z*(5796.d1+z*(276.d0+z))))))/14572069319808.d3
	eta4=z*eta6+11.d0*eta5
	eta3=z*eta5+9.d0*eta4
	eta2=z*eta4+7.d0*eta3
	eta1=z*eta3+5.d0*eta2
	eta0=z*eta2+3.d0*eta1
	csi=z*eta1+eta0
	return
	end
cc=======================================================================
	complex*16 function cwn(i,e)
	implicit real*8 (a-h,o-z)
cc This function calculates the complex wavenumber cwn from the
cc square of the complex energy  times 2m/h**2 (k**2=e*2m/h**2)
	complex*16 e
	cwn=sqrt(e)
	if(dimag(cwn).eq.0)return
	if(dimag(cwn).gt.0)then
	  if(i.eq.-1)return
	  cwn=-cwn
	endif
	if(i.eq.-1)cwn=-cwn
	return
	end	
cc=======================================================================
        subroutine sli(rtol,pas,v,v1,v2,par,rlorb,spot,pot,intf)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc                                                                      cc
cc This sbr. computes the partition and the associated vectors          cc
cc (v,v1,v2) to be further used by sbr. slix which solves the           cc
cc initial value problem for the radial schroedinger equation on i2     cc
cc and on some subinterval of i3.                                       cc
cc The equation to be solved is                                         cc
cc                      y" = (v(x)- k**2)*y                             cc
cc where the complex function v(x) has the form                         cc
cc                                                                      cc
cc                v(x) = l*(l +1)/x**2 +spot(x)/x +pot(x)               cc
cc                                                                      cc
cc Input data(to be declared in the calling program):                   cc
cc rtol   Accuracy requested in the results; rtol is a real double      cc
cc        precision number                                              cc
cc par    Double complex vector with two components. par(1) and par(2)  cc
cc        are the ends of the x domain of the system. The domain is     cc
cc        supposed to be the segment in the complex plane which joints  cc
cc        par(1) and par(2).                                            cc
cc rlorb  double precision number for the orbital momentum l.           cc
cc                                                                      cc
cc Externals furnished by the user:                                     cc
cc spot   Runction of the form                                          cc
cc                  double complex function spot(x)                     cc
cc        which furnishes spot(x), for any given double complex x.      cc
cc pot    Runction of the form                                          cc
cc                  double complex function pot(x)                      cc
cc        which allows calculating  pot(x) for any given x              cc
cc                                                                      cc
cc input/utput parameter:                                               cc
cc intf   In input this furnishes the maximum number of intervals       cc
cc        assigned to the partition. In output it furnishes the         cc
cc        total number of intervals of the resultant partition          cc
cc        consistent with the required tolerance rtol. If               cc
cc        it is found that this number exceeds the input value          cc
cc        the program displayes the message                             cc
cc                   "insufficient number of steps."                    cc
cc pas    Double complex vector; pas(int) is the stepsize of the        cc
cc        int-th interval of the partition.                             cc
cc v      Double complex vector. Value v(int), represents the gross     cc
cc        part (conveniently normalized) of the potential on the        cc
cc        interval int.                                                 cc
cc v1     Double complex vector. Value v1(int) represents the           cc
cc        (conveniently normalized) linear correction in the            cc
cc        potential on the interval int.                                cc
cc v2     The same as v1 for the quadratic correction.                  cc
cc                                                                      cc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        implicit double complex (a-h,o-q,s-z)
        implicit double precision (r)
        dimension par(2),pas(1),v(1),v1(1),v2(1)
	external spot,pot
        rlotol=rtol**(7.d0/6.d0)
        rl=(1.d-4)**(7.d0/6.d0)
        fac=(1.d2*rl/rlotol)**(1.d0/7.d0)
        nm=fac
        rnm=nm
        rld=abs(par(2)-par(1))
        versor=(par(2)-par(1))/rld
        int=0
        x1=par(1)
        pas(1)=versor
        ie=0
 10	int=int+1
        icont=0
 12	continue
	call sprel(x1,pas(int),v(int),v1(int),v2(int),spot,pot,rlorb)
        icont=icont+1
        if(ie.eq.1) goto 15
        rv1=abs(v1(int)*v1(int))
        rhi=abs(pas(int))
        rv1=rv1/(rhi*rhi*rhi*rhi*rhi*rhi)
        rvbar=1.d-1
        rw=rv1*rvbar
        rdh=.1d0*rhi
        call hsize(rw,rdh,rh,rl)
        rdist=abs(x1+versor*rh-par(1))
        if(rdist.gt.rld-0.01*rh) ie=1
        if(int.eq.intf) ie=1
        if(int.eq.intf) write(*,*)' insufficient number of steps' 
        if(int.eq.intf) write(7,*)' insufficient number of steps' 
        if(ie.eq.1) pas(int)=par(2)-x1
        if(ie.eq.1) goto 12
	if(icont.ge.3) goto 16
        if(abs(rh-rhi)/rhi.lt..11d0) goto 16
        pas(int)=versor*rh
        goto 12
 16	x1=x1+pas(int)
        pas(int+1)=pas(int)
        goto 10
 15	intf1=int
        ia=0
        if(abs(pas(intf1)).lt.1.d-12)ia=1
        if(ia.eq.1)intf1=intf1-1
        if(ia.eq.1)pas(intf1)=pas(intf1)+pas(intf1+1)
        intf2=intf1*nm
        if(intf2.gt.intf)write(*,*)'insufficient number of steps'
        if(intf2.gt.intf)write(7,*)'insufficient number of steps'
        if(intf2.gt.intf)stop
        intf=intf1
        do 1 int=1,intf
	  inm=(int-1)*nm+1
	  fac=pas(int)/rnm
	  do 2 i=inm,inm+nm-1
	    v(i)=fac
 2	  continue
 1	continue
        intf=intf*nm
        do 3 int=1,intf
	  pas(int)=v(int)
 3	continue
        x1=par(1)
        do 20 k=1,intf
	  call sprel(x1,pas(k),v(k),v1(k),v2(k),spot,pot,rlorb)
	  x2=x1+pas(k)
	  x1=x2
 20	continue
        return
        end
cc=======================================================================
        subroutine hsize(rw,rdh,rh,rlotol)
        implicit double precision (r)
        i=1
 1	i=i+1
        if(i.ge.22) return
        m=0
        rh=rdh*(i-1)
        rh2=rh*rh
        rh6=rh2*rh2*rh2
        rw8=rw*rh6*rh2
        if(rw8/126.d1.le.rlotol) m=m+1
        if(rw8*rh/1134.d1.le.rlotol) m=m+1
        if(rw8/(9.d1*rh).le.rlotol) m=m+1
        if(m.eq.3) goto 1
        return
        end
cc=======================================================================
	subroutine sprel(x1,h,v,v1,v2,spot,pot,rlorb)
	implicit double complex (a-h,o-q,s-z) 
	implicit double precision (r)
	dimension ppot(6),pspot(6)
	external spot,pot
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  this sbr. computes  v,v1,v2  on a single interval
c  whose ends are x1 and x1+h.          
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        rl=rlorb*(rlorb+1)
        call ssixval(pot,x1,h,ppot)
        call ssixval(spot,x1,h,pspot)
        do 3 m=1,6
	  x=x1+.2d0*h*(m-1)
	  ppot(m)=ppot(m)+pspot(m)/x+rl/(x*x)
 3	continue
        call sfitpot(ppot,h,v,v1,v2)
	return
	end
cc=======================================================================
        subroutine spi1(x1,vc,v0,v1,v2)
        implicit double complex (a-h,o-q,s-z)
        implicit double precision (r)
        dimension cpfs(6),cpf(6)
        external spot,pot
        data czero/(0.d0,0.d0)/
        call ssixval(spot,czero,x1,cpfs)
        call sfitpot(cpfs,x1,vs0,vs1,vs2)
        call ssixval(pot,czero,x1,cpf)
        call sfitpot(cpf,x1,vn0,vn1,vn2)
        call copot(vs0,vs1,vs2,vn0,vn1,vn2,x1,vts,vt0,vt1,vt2)
        vc=vts/x1
        v0=vt0/(x1*x1)
        v1=vt1/(x1*x1*x1)
        v2=vt2/(x1*x1*x1*x1)
        return
        end
cc=======================================================================
        subroutine copot(vs0,vs1,vs2,v0,v1,v2,h,vts,vt0,vt1,vt2)
        implicit double complex (a-h,o-q,s-z)
        data two/(2.d0,0.d0)/,thr/(3.d0,0.d0)/,fou/(4.d0,0.d0)/,
     *  twe/(12.d0,0.d0)/
        vts=(vs0-two*(vs1-vs2))/h
        vt0=v0-two*(v1-v2)+fou*(vs1-thr*vs2)/h
        vt1=fou*(v1-thr*v2)+twe*vs2/h
        vt2=twe*v2
        return
        end
cc=======================================================================
        subroutine sfitpot(v,h,v0,v1,v2)
        implicit double complex (a-h,o-q,s-z)
        dimension v(6)
        hh=h*h
	v0=hh*(19.d0*(v(1)+v(6))+75.d0*(v(2)+v(5))+5.d1*(v(3)+v(4)))/
     1288.d0
	v1=hh*(111.d0*(v(6)-v(1))+425.d0*(v(5)-v(2))
     1-15.d1*(v(4)-v(3)))/1344.d0
	v2=hh*25.d0*(5.d0*(v(1)+v(6))+6.d0*(v(2)+v(5))
     1-11.d0*(v(3)+v(4)))/1008.d0
	return
        end        
cc=======================================================================
        subroutine ssixval(func,x,h,pf)
        implicit double precision (r)
        implicit double complex (a-h,o-q,s-z)
        dimension pf(6)
        external func
        do 1 k=1,6
	  pf(k)=func(x+.2d0*h*(k-1))
 1	continue
        return
        end
cc=======================================================================
cc evaluated the singular and non-singular part of the potential and
cc their second derivative
	subroutine potgen
	parameter(mxrd=1000)
	implicit double complex (a-h,o-q,s-z)
	implicit real*8 (r)
	common/nstep/rstep,nstep
	dimension rxx(mxrd)
	dimension ysps(mxrd),yspn(mxrd),bsps(mxrd),bspn(mxrd)
	nsps=nstep+1
	nspn=nstep+1
cc singular and non-singular potential for r<rmax
	write(7,*)'The files spot.dat and nspot.dat store the singular'
	write(7,*)'and non-singular potential(chequear el denominador)'
	open(unit=8,file='Out/spot.dat',status='unknown')
	open(unit=9,file='Out/nspot.dat',status='unknown')
	do i=1,nsps
	 rxx(i)=0.05d0+i*rstep
	 ysps(i)=spot(dcmplx(rxx(i))) !singular potential
	 yspn(i)=pot(dcmplx(rxx(i))) !non-singular potential
	 write(8,900) rxx(i),ysps(i)/rxx(i)
	 write(9,900) rxx(i),yspn(i)
	enddo
 900	format(1x,f7.3,2e12.3)
	close(unit=8)
	close(unit=9)
cc calculates the second derivative (bsps) of ysps to be used in the
cc interpolation subroutine spintc
	ypmin=(ysps(2)-ysps(1))/rstep
	ypmax=(ysps(nsps)-ysps(nsps-1))/rstep
	call splinc(rxx,ysps,nsps,ypmin,ypmax,bsps)
cc calculates the second derivative (bspn) of yspn to be used in the
cc interpolation subroutine spintc
	ypmin=(yspn(2)-yspn(1))/rstep
	ypmax=(yspn(nspn)-yspn(nspn-1))/rstep
	call splinc(rxx,yspn,nspn,ypmin,ypmax,bspn)
	return
	end
cc=======================================================================
cc calculation of the spline coefficients for knots points x(i),y(i),i=1,..,n
cc and store them in y2(1),..,y2(n)
	subroutine splinc(rx,y,n,yp1,ypn,y2)
        parameter(mni=1500,mnia=300)
	implicit double complex (a-h,o-q,s-z)
	implicit double precision (r)
	dimension rx(n),y(n),y2(n),u(mni)
	if (dreal(yp1).gt..99e30) then
	   y2(1)=0.d0
	   u(1)=0.d0
	else
	   y2(1)=-0.5d0
	   u(1)=(3.d0/(rx(2)-rx(1)))*((y(2)-y(1))/(rx(2)-rx(1))-yp1)
	endif
	do 11 i=2,n-1
	   rsig=(rx(i)-rx(i-1))/(rx(i+1)-rx(i-1))
	   p=rsig*y2(i-1)+2.d0
	   y2(i)=(rsig-1.d0)/p
	   u(i)=(6.d0*((y(i+1)-y(i))/(rx(i+1)-rx(i))-(y(i)-y(i-1))
     *      /(rx(i)-rx(i-1)))/(rx(i+1)-rx(i-1))-rsig*u(i-1))/p
 11	continue
	if (dreal(ypn).gt..99e30) then
	   qn=0.d0
	   un=0.d0
	else
	   qn=0.5d0
	   un=(3.d0/(rx(n)-rx(n-1)))*(ypn-(y(n)-y(n-1))/(rx(n)-rx(n-1)))
	endif
	y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
	do 12 k=n-1,1,-1
	   y2(k)=y2(k)*y2(k+1)+u(k)
 12	continue
	return
	end
cc=========================================================
cc non-singular potential terms
	double complex function pot(zr)
	implicit double complex (a-h,o-q,s-z)
	implicit double precision (r)
	common/rbigr/rbigr,rdiff,rv0,rvso
	common/iz/iz
	common/rzt/rzt
	common/rc1/rc1
cc W-S term
	if(rdiff.lt.1.e-2)then
	  vvv=rv0
	  if(abs(zr).gt.rbigr) vvv=0.d0
	else
	  e1=exp((rbigr-zr)/rdiff)
	  vvv=rv0*e1/(1.d0+e1)
	endif
cc Coulomb term
	vcou=0.d0
	if(iz.gt.1) then
	  cost=1.43986d0*rzt
	  if(abs(zr).le.rbigr) vcou=cost/2.d0*(3.d0/rbigr-zr*zr/rbigr**3)
	endif
	pot=rc1*(-vvv+vcou)
	return
	end
cc===============================================================
cc calculation of the singular potential which behaves as 1/r
	double complex function spot(zr)
	implicit double complex (a-h,o-q,s-z)
	implicit double precision (r)
	common/lj/lsp,jsp
	common/iz/iz
	common/rbigr/rbigr,rdiff,rv0,rvso
	common/rzt/rzt
	common/rc1/rc1
cc spin-orbit term
	if(jsp.eq.0)then
	  als=0.d0
	  vls=0.d0
	else
	  e1=exp((rbigr-zr)/rdiff)
	  als=0.25d0*(jsp*(jsp+2)/4.d0-lsp*(lsp+1)-0.75d0)
	  vls=als*4.d0/rdiff*e1/(1.d0+e1)**2*rvso
	endif
	vcou=0.d0
cc point Coulomb term outside
	if(iz.gt.1)then
	  cost=1.43986d0*rzt
	  if(abs(zr).gt.rbigr) vcou=cost
	endif
	spot=rc1*(vcou-vls)
	return
	end
cc=======================================================================
	double complex function aspot(zr)
cc singular potential in the asymptotic region
	implicit double complex (a-h,o-q,s-z)
	implicit double precision (r)
	common/iz/iz
	common/rzt/rzt
	common/rbigr/rbigr,rdiff,rv0,rvso
	common/rc1/rc1
cc Coulomb term
	vcou=0.d0
        if(iz.gt.1)then
	  cost=1.43986d0*rzt
	  if(dreal(zr).gt.rbigr) vcou=cost
        endif
	aspot=rc1*vcou
	return
	end
cc=======================================================================
	double complex function apot(zr)
cc non-singular term in the asymptotic region: none
	implicit double complex (a-h,o-q,s-z)
	implicit double precision (r)
	apot=(0.d0,0.d0)
	return
	end
cc=======================================================================
cc reads the single particle state
	subroutine spdat(isp,esp)
	implicit double complex (a-h,o-q,s-z)
	implicit real*8 (r)
	common/lj/lsp,jsp
	open(unit=1,file='Read/sp.dat',status='old')
	read(1,*)isp,lsp,jsp,rener,renei
	close(unit=1)
	if(isp.lt.-2.or.isp.gt.1)stop 'isp should be -2,-1,0 or 1'
	esp=dcmplx(rener,renei)
	if(isp.eq.-2) write(7,*) ' Anti-bound state readed'
	if(isp.eq.-1) write(7,*) ' Bound state readed'
	if(isp.eq.0)  write(7,*) ' Resonant state readed'
	if(isp.eq.1)  write(7,*) ' Scattering state readed'
	write(7,901) lsp,jsp,esp
 901	format(1x,' l,2*j,e -> ',2i5,2f10.3/)
	return
	end
cc=======================================================================
	subroutine potdat
	parameter(mxrd=300)
	implicit double complex (a-h,o-q,s-z)
	implicit real*8 (r)
	common/nstep/rstep,nstep
	common/rmax/rmax
	common/iz/iz
	common/rbigr/rbigr,rdiff,rv0,rvso
	common/rzt/rzt
	common/rc1/rc1
	common/rmatch/rmatch
	common/retak/retak 
	open(unit=5,file='Read/pot.dat',status='old')
	read(5,*) rmax,rstep
	read(5,*) rmatch
	read(5,*) iz
	read(5,*) rat,rzt,rv0,r0,rdiff,rvso
	close(unit=5)
	nstep=rmax/rstep+0.1
	if(nstep.gt.mxrd) stop ' nstep.gt.mxrd. I stop'
	rbigr=r0*(rat**(1./3.))
	if(iz.ne.1.and.iz.ne.2)  stop 'iz should be 1(n) or 2(p). I stop'
	if(iz.eq.1) then
	  rap=1.0089d0 !neutron mass
	else
	  rap=1.0075d0 !proton mass
	endif
	red2m=rat*rap/(rat+rap)
	rc1=0.04783258d0*red2m !k^2=c1*e,0.478325=2/hbar^2 in (Mev*u*fm^2)^{-1}
	rzet=rzt
	if(iz.eq.1)rzet=0.d0
	retak=0.71993d0*rzet*rc1 !coulomb parameter*wave number
	write(7,1799)
 1799	format(1x,'SINGLE-PARTICLE HAMILTONIAN read from pot.dat')
	write(7,1800) rmax
 1800	format(1x,'Maximun radius at which u(r) is evaluated -> ',f7.3)
	write(7,1801) rstep
 1801	format(1x,'Step at which u(r) is evaluated -> ',f7.3)
	write(7,1802) nstep
 1802	format(1x,'Mesh point number -> ',i5)
	write(7,1804) rmatch
 1804	format(1x,'Matching radius -> ',f7.3)
	write(7,1805) iz
 1805	format(1x,'1(2) for neutrons(protons) -> ',i4)
	write(7,1809)
 1809	format(1x,'W-S parameters: at,zt,v0,r0,a,vso ')
	write(7,1806) rat,rzt,rv0,r0,rdiff,rvso
 1806	format(1x,6f7.3)
	write(7,1808) rbigr
 1808	format(1x,'Nuclear radius -> ',f7.3)
	write(7,1810) rc1
 1810	format(1x,'Convertion Factor c1: k^2=c1*e -> ',e15.5)
	write(7,1803) retak
 1803	format(1x,'Coulomb parameter*wave number -> ',e15.5/)
	return
	end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

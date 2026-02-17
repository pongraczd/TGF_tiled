#include "fintrf.h"  
C "fintrf.h": MATLAB/FORTRAN interface header file
C using the MEX file 'mexFunction' to call the Fortan subroutine 'One_polyhedron' 
C-----------------------------------------------------------------------
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
C
C MEX file: mexFunction
C nlhs: number of input
C plhs: input variant
C nrhs: number of output
C prhs: output variant
C
C-----------------------------------------------------------------------
C     Declarations 
      implicit none
C     mexFunction arguments:
      integer nlhs, nrhs, num(100)
      mwPointer plhs(*), prhs(*)

C     Function declarations:
      mwPointer mxGetPr, mxGetM, mxGetN
      mwPointer mxCreateDoubleMatrix
C mxGetPr: Pointer to the first element of real data
C mxGetM: return the number of rows
C mxGetN: return the number of columns
C mxCreateDoubleMatrix: create an m-by-n mxArray.
C
C     Array information, pointers to input/output mxArrays: 
      mwPointer x_pr, y_pr,poly_prx,poly_pry,poly_prz,poly_pre
      mwPointer m, n, m_out, n_out
      mwSize size,size_out    
      real*8 x(100),y(1000000,10)
	  real*8 poly_x(10000000),poly_y(10000000),poly_z(10000000)
        m_out=1000000
		n_out=10
		size_out=m_out*n_out
C-----------------------------------------------------------------------
C obtaining the size of input array
      m = mxGetM(prhs(1))
      n = mxGetN(prhs(1))
      size = m * n
      if(size .gt. 100) then
         call mexErrMsgIdAndTxt ('MATLAB:InputTooBig',
     +              'Input #1: number of elements exceeds buffer') 
      endif
C
C obtain input data from Matlab
      x_pr = mxGetPr(prhs(1))
C transform input data to double pricision data, and give values to Fortran variant x
      call mxCopyPtrToReal8(x_pr,x,size)
      num=int(x)
	  
c obtain size of prhs
      m = mxGetM(prhs(2))
      n = mxGetN(prhs(2))
      size = m * n
	       poly_prx = mxGetPr(prhs(2))
	      call mxCopyPtrToReal8(poly_prx,poly_x,size)
		   poly_pry = mxGetPr(prhs(3))
	      call mxCopyPtrToReal8(poly_pry,poly_y,size)
		   poly_prz = mxGetPr(prhs(4))
	      call mxCopyPtrToReal8(poly_prz,poly_z,size) 
C call Fortran function "One_polyhedron"
      call One_polyhedron(num,poly_x,poly_y,poly_z,y)

C send calculated results from Fortran to Matlab
C creat a double pricision array
      plhs(1) = mxCreateDoubleMatrix(m_out,n_out,0) 
      y_pr = mxGetPr(plhs(1))
      call mxCopyReal8ToPtr(y,y_pr,size_out)      

      return
      end

      subroutine One_polyhedron(num,poly_x,poly_y,poly_z,T_g
     $ )

	implicit real*8 (a-h,o-z)
c
c	Edit number of planes (nop) and number of vertices (nov)
c
	parameter(nop=7,nov=8,noed=100)

	dimension e(3,3)
	dimension noe(nop)
	dimension x(nov),y(nov),z(nov),h(nop),s(nop),cs(nop,3)
	dimension xo(nop),yo(nop),zo(nop)
	dimension xx(nop),yy(nop),zz(nop)
	dimension xxx(nop,noed),yyy(nop,noed),zzz(nop,noed)
	dimension k(nop,noed)
	dimension g(nop,noed,3)
	dimension plnorm(nop,3),segnorm(nop,noed,3)
	dimension dcnp(nop,3),dcnpq(nop,noed,3)
	dimension hh(nop,noed),ss(nop,noed)
	dimension cspq(nop,noed,3)
	dimension vij(3,3),x_temp(3*nov)
	dimension poly_x(10000000),poly_y(10000000),poly_z(10000000)
	dimension T_g(1000000,10)
	integer num,p_j,p_i

	snorm(xipq,yipq,zipq)=dsqrt(xipq*xipq+yipq*yipq+zipq*zipq)
	fln(aa,bb,cc,dd)=dlog((aa+bb)/(cc+dd))
	fan(ee,f,gg,t,yf,w)=datan((ee*t)/(f*w))-datan((ee*gg)/(f*yf))


	noe=(/3, 3, 4, 4, 4, 4, 4/)

  
	k(1,1:noe(1))= (/1, 2, 3/)
	k(2,1:noe(2))= (/3, 4, 1/)
	k(3,1:noe(3))= (/6, 5, 8, 7/)
	k(4,1:noe(4))= (/8, 5 ,1, 4/)
	k(5,1:noe(5))= (/6, 7, 3, 2/)
	k(6,1:noe(6))= (/5, 6, 2, 1/)
	k(7,1:noe(7))= (/7, 8, 4, 3/)
   


c
c	Vector basis (ei) situated at P
c

	e(1,1)=1
	e(1,2)=0
	e(1,3)=0
	e(2,1)=0
	e(2,2)=1
	e(2,3)=0
	e(3,1)=0
	e(3,2)=0
	e(3,3)=1

c      do 2001 i=1,10
c         T_g(i)=0.0
c2001  continue
c
c	Constants
c
	gc=6.67259e-08
	d=1.0
c    d=1.0
	pi=4.*atan(1.)

c
c
c	Building vectors Gij // Equation (17)
c
c
      do 2100 p_j=1,num
c	     read(20,*,iostat=ios) x_temp(1:3*nov)
         do 2000 p_i=1,nov		 
c	     x(p_i)=x_temp((p_i-1)*3+1)
c		 y(p_i)=x_temp((p_i-1)*3+2)
c		 z(p_i)=x_temp((p_i-1)*3+3)
         x(p_i)=poly_x((p_i-1)*num+p_j)
		 y(p_i)=poly_y((p_i-1)*num+p_j)
		 z(p_i)=poly_z((p_i-1)*num+p_j)
2000	 continue

 
  	do 70 i=1,nop
	do 70 j=1,noe(i)

	if(j.eq.noe(i)) then 

	g(i,j,1)=x(k(i,1))-x(k(i,j))
	g(i,j,2)=y(k(i,1))-y(k(i,j))
	g(i,j,3)=z(k(i,1))-z(k(i,j))


	else


	g(i,j,1)=x(k(i,j+1))-x(k(i,j))
	g(i,j,2)=y(k(i,j+1))-y(k(i,j))
	g(i,j,3)=z(k(i,j+1))-z(k(i,j))

	end if


70	continue

c
c
c	Computation of plane Normals (Np) // Equation (18)
c
c

	do 80 i=1,nop
	call cross (g(i,1,1),g(i,1,2),g(i,1,3),g(i,2,1),g(i,2,2),
     $g(i,2,3),plnorm(i,1),plnorm(i,2),plnorm(i,3))
	qnorm=snorm(plnorm(i,1),plnorm(i,2),plnorm(i,3))
	plnorm(i,1)=plnorm(i,1)/qnorm
	plnorm(i,2)=plnorm(i,2)/qnorm
	plnorm(i,3)=plnorm(i,3)/qnorm


c
c	Direction cosines computation
c
c
c	Notice: We computed a normed vector, i.e. Np's 
c	        x, y and z components are at the same time
c	        its direction cosines with the three axes.
c

	dcnp(i,1)=plnorm(i,1)
	dcnp(i,2)=plnorm(i,2)
	dcnp(i,3)=plnorm(i,3)


80	continue


c
c
c	cos(Np,ei), cos(npq,ei)
c
c
c	cs(i,j) = dcnp(i,j)
c
c	and
c
c	cspq(i,j,m) = dcnpq(i,j,m), m = 1,2,3
c
c

	do 23 i=1,nop
	do 23 j=1,3
	qnorm1=snorm(plnorm(i,1),plnorm(i,2),plnorm(i,3))
	qnorm2=snorm(e(j,1),e(j,2),e(j,3))
	call dot (plnorm(i,1),plnorm(i,2),plnorm(i,3),e(j,1),e(j,2),
     $e(j,3),dotnpi)
	cs(i,j)=dotnpi/(qnorm1*qnorm2)


23	continue

c
c
c	Computation of line segment normals (npq) // Equation (19)
c
c

	do 90 i=1,nop
	do 90 j=1,noe(i)
	call cross (g(i,j,1),g(i,j,2),g(i,j,3),plnorm(i,1),
     $plnorm(i,2),plnorm(i,3),segnorm(i,j,1),
     $segnorm(i,j,2),segnorm(i,j,3))
	qnorm=snorm(segnorm(i,j,1),segnorm(i,j,2),segnorm(i,j,3))
	segnorm(i,j,1)=segnorm(i,j,1)/qnorm
	segnorm(i,j,2)=segnorm(i,j,2)/qnorm
	segnorm(i,j,3)=segnorm(i,j,3)/qnorm


c
c	Direction cosines computation
c
c
c	Notice: We are computing a normed vector, i.e. npq's 
c	        x, y and z components are at the same time
c	        its direction cosines with the three axes.
c
	dcnpq(i,j,1)=segnorm(i,j,1)
	dcnpq(i,j,2)=segnorm(i,j,2)
	dcnpq(i,j,3)=segnorm(i,j,3)


90	continue


	
	do 24 i=1,nop
	do 24 j=1,noe(i)
 	do 24 m=1,3
	qnorm1=snorm(segnorm(i,j,1),segnorm(i,j,2),segnorm(i,j,3))
	qnorm2=snorm(e(m,1),e(m,2),e(m,3))
	call dot (segnorm(i,j,1),segnorm(i,j,2),segnorm(i,j,3),e(m,1),
     $e(m,2),e(m,3),dotnpq)
	cspq(i,j,m)=dotnpq/(qnorm1*qnorm2)


24	continue

c
c
c	Distances of P from planes Sp ( h(i) ) 
c	by subroutine 'planedist'
c
c	xo(i), yo(i) and zo(i) are the intersections of Sp(i)
c	with axes x, y and z respectively
c
c	and 
c
c	plane normal orientation ( s(i) ) // Equation (20)
c
c

	do 60 i=1,nop
	call planedist (x,y,z,k(i,1),k(i,2),k(i,3),h(i),xo(i),yo(i),
     $zo(i))


60	continue



	do 41 i=1,nop
	call dot (plnorm(i,1),plnorm(i,2),plnorm(i,3),-x(k(i,1)),
     $-y(k(i,1)),-z(k(i,1)),plnrmor)
	if (plnrmor.gt.0) then
	s(i)=-1
	else if (plnrmor.lt.0) then
	s(i)=1
	else if (plnrmor.eq.0) then
	s(i)=0
	end if


41	continue





c
c
c	Computation of the projections of P on the 'p' planes
c	of the polyhedron. [P'(x',y',z')]
c
c	xx(i), yy(i), zz(i) // Equation (21)
c
c

	do 31 i=1,nop

	if (xo(i).ge.0.and.yo(i).ge.0.and.zo(i).ge.0) then
	xx(i)=abs(dcnp(i,1)*h(i))
	yy(i)=abs(dcnp(i,2)*h(i))
	zz(i)=abs(dcnp(i,3)*h(i))

	else if (xo(i).lt.0.and.yo(i).lt.0.and.zo(i).lt.0) then

	 if (dcnp(i,1).gt.0) then
	 xx(i)=-dcnp(i,1)*h(i)
	 else
	 xx(i)=dcnp(i,1)*h(i)
	 end if

	 if (dcnp(i,2).gt.0) then
	 yy(i)=-dcnp(i,2)*h(i)
	 else
	 yy(i)=dcnp(i,2)*h(i)
	 end if

	 if (dcnp(i,3).gt.0) then
	 zz(i)=-dcnp(i,3)*h(i)
	 else
	 zz(i)=dcnp(i,3)*h(i)
	 end if

	else if (xo(i).lt.0.and.yo(i).lt.0) then

	 if (dcnp(i,1).gt.0) then
	 xx(i)=-dcnp(i,1)*h(i)
	 else
	 xx(i)=dcnp(i,1)*h(i)
	 end if

	 if (dcnp(i,2).gt.0) then
	 yy(i)=-dcnp(i,2)*h(i)
	 else
	 yy(i)=dcnp(i,2)*h(i)
	 end if

	zz(i)=abs(dcnp(i,3)*h(i))

	else if (yo(i).lt.0.and.zo(i).lt.0) then
	xx(i)=abs(dcnp(i,1)*h(i))

	 if (dcnp(i,2).gt.0) then
	 yy(i)=-dcnp(i,2)*h(i)
	 else
	 yy(i)=dcnp(i,2)*h(i)
	 end if

	 if (dcnp(i,3).gt.0) then
	 zz(i)=-dcnp(i,3)*h(i)
	 else
	 zz(i)=dcnp(i,3)*h(i)
	 end if

	else if (xo(i).lt.0.and.zo(i).lt.0) then

	 if (dcnp(i,1).gt.0) then
	 xx(i)=-dcnp(i,1)*h(i)
	 else
	 xx(i)=dcnp(i,1)*h(i)
	 end if

	yy(i)=abs(dcnp(i,2)*h(i))

	 if (dcnp(i,3).gt.0) then
	 zz(i)=-dcnp(i,3)*h(i)
	 else
	 zz(i)=dcnp(i,3)*h(i)
	 end if

	else if (xo(i).lt.0) then

	 if (dcnp(i,1).gt.0) then
	 xx(i)=-dcnp(i,1)*h(i)
	 else
	 xx(i)=dcnp(i,1)*h(i)
	 end if

	yy(i)=abs(dcnp(i,2)*h(i))
	zz(i)=abs(dcnp(i,3)*h(i))

	else if (yo(i).lt.0) then
	xx(i)=abs(dcnp(i,1)*h(i))

	 if (dcnp(i,2).gt.0) then
	 yy(i)=-dcnp(i,2)*h(i)
	 else
	 yy(i)=dcnp(i,2)*h(i)
	 end if

	zz(i)=abs(dcnp(i,3)*h(i))

	else if (zo(i).lt.0) then
	xx(i)=abs(dcnp(i,1)*h(i))
	yy(i)=abs(dcnp(i,2)*h(i))

	 if (dcnp(i,3).gt.0) then
	 zz(i)=-dcnp(i,3)*h(i)
	 else
	 zz(i)=dcnp(i,3)*h(i)
	 end if

	end if


31	continue




c
c	segment normal orientation	
c
c	senrmor = npq (dot) Gpq(ij-->P') // Equation (22)
c

	do 42 i=1,nop
	do 42 j=1,noe(i)
	call dot (segnorm(i,j,1),segnorm(i,j,2),segnorm(i,j,3),
     $xx(i)-x(k(i,j)),yy(i)-y(k(i,j)),zz(i)-z(k(i,j)),senrmor)
	if (abs(senrmor).lt.1e-10) then
	ss(i,j)=0
	go to 42
	else if (senrmor.gt.0) then
	ss(i,j)=-1
	else if (senrmor.lt.0) then
	ss(i,j)=1
	else if (senrmor.eq.0) then
	ss(i,j)=0
	end if


42	continue



c
c
c	Distances of P' from line segments Gpq (hpq) 
c
c	hh(i,j)
c
c	Computation of the projections of P' on each
c	line segment Gpq, i.e. P''(x'',y''',z''')
c	by calling subroutine 'segmproj'
c
c	xxx(i,j), yyy(i,j), zzz(i,j) // Equations (23) - (25)
c
c

	do 32 i=1,nop
	do 32 j=1,noe(i)

	if (j.eq.noe(i)) then

	if (ss(i,j).eq.0) then
	hh(i,j)=0
	xxx(i,j)=xx(i)
	yyy(i,j)=yy(i)
	zzz(i,j)=zz(i)
	go to 131
	end if

	call segmproj (x(k(i,j)),y(k(i,j)),z(k(i,j)),
     $x(k(i,1)),y(k(i,1)),z(k(i,1)),
     $xx(i),yy(i),zz(i),xseg,yseg,zseg)
	hh(i,j)=snorm(xseg-xx(i),yseg-yy(i),zseg-zz(i))


	xxx(i,j)=xseg
	yyy(i,j)=yseg
	zzz(i,j)=zseg
	go to 131

	end if

	if (ss(i,j).eq.0) then
	hh(i,j)=0
	xxx(i,j)=xx(i)
	yyy(i,j)=yy(i)
	zzz(i,j)=zz(i)
	go to 131
	end if

	call segmproj (x(k(i,j)),y(k(i,j)),z(k(i,j)),
     $x(k(i,j+1)),y(k(i,j+1)),z(k(i,j+1)),
     $xx(i),yy(i),zz(i),xseg,yseg,zseg)
	hh(i,j)=snorm(xseg-xx(i),yseg-yy(i),zseg-zz(i))


	xxx(i,j)=xseg
	yyy(i,j)=yseg
	zzz(i,j)=zseg

131	continue




32	continue



c
c
c	Computation of V and Vi
c
c

	v=0
	vx=0
	vy=0
	vz=0

	do 100 i=1,nop

	sum1=0
	sum2=0
	singarea=0
	singsegm=0
	singvert=0
	
	do 200 j=1,noe(i)




	thisvert=0

	gpqnorm=snorm(g(i,j,1),g(i,j,2),g(i,j,3))

c	
c	check if P' lays inside the polygon Np(i)
c

	if (ss(i,j).eq.1) then
	singarea=singarea+1
	end if

c
c	check if P' lays inside the segment Gpq 
c
c	or
c
c	on the vertex (i,j)
c
c	Condition controlled is  
c
c	ss(i,j) = 0
c

	if (ss(i,j).eq.0) then



	 if(j.eq.noe(i)) then 

	 gx1=xx(i)-x(k(i,1))
	 gx2=xx(i)-x(k(i,j))
	 gy1=yy(i)-y(k(i,1))
	 gy2=yy(i)-y(k(i,j))
	 gz1=zz(i)-z(k(i,1))
	 gz2=zz(i)-z(k(i,j))

	 else  

	 gx2=xx(i)-x(k(i,j))
	 gx1=xx(i)-x(k(i,j+1))
	 gy2=yy(i)-y(k(i,j))
	 gy1=yy(i)-y(k(i,j+1))
	 gz2=zz(i)-z(k(i,j))
	 gz1=zz(i)-z(k(i,j+1))

	 end if

	 enorm1=snorm(gx1,gy1,gz1)
	 enorm2=snorm(gx2,gy2,gz2)




c
c	check for the Vertex (i,j) 
c

	  if (enorm1.eq.0) then

	  singvert=singvert+1

	  thisvert=thisvert+1

	   	if (j.eq.noe(i)) then

	   g1norm=snorm(g(i,j,1),g(i,j,2),g(i,j,3))
	   g2norm=snorm(g(i,1,1),g(i,1,2),g(i,1,3))

	call dot(-g(i,j,1),-g(i,j,2),-g(i,j,3),g(i,1,1),g(i,1,2),
     $g(i,1,3),gdot)

	     		if (gdot.eq.0) then

	     theta=pi/2

	     		else

	     theta=dacos(gdot/(g1norm*g2norm))

	     		end if

	   	else

	   g1norm=snorm(g(i,j,1),g(i,j,2),g(i,j,3))
	   g2norm=snorm(g(i,j+1,1),g(i,j+1,2),g(i,j+1,3))

	call dot(-g(i,j,1),-g(i,j,2),-g(i,j,3),g(i,j+1,1),g(i,j+1,2),
     $g(i,j+1,3),gdot)

	     		if (gdot.eq.0) then

	     theta=pi/2

	     		else

	     theta=dacos(gdot/(g1norm*g2norm))

	     		end if

	   	end if

	  end if


	  if (enorm2.eq.0) then

	  singvert=singvert+1

	  thisvert=thisvert+1

	   	if (j.eq.1) then

	   g1norm=snorm(g(i,noe(i),1),g(i,noe(i),2),g(i,noe(i),3))
	   g2norm=snorm(g(i,1,1),g(i,1,2),g(i,1,3))

	call dot(-g(i,noe(i),1),-g(i,noe(i),2),-g(i,noe(i),3),
     $g(i,1,1),g(i,1,2),g(i,1,3),gdot)
     
	     		if (gdot.eq.0) then

	     theta=pi/2

	     		else

	     theta=dacos(gdot/(g1norm*g2norm))

	     		end if

	   	else

	   g1norm=snorm(g(i,j-1,1),g(i,j-1,2),g(i,j-1,3))
	   g2norm=snorm(g(i,j,1),g(i,j,2),g(i,j,3))

	call dot(-g(i,j-1,1),-g(i,j-1,2),-g(i,j-1,3),g(i,j,1),g(i,j,2),
     $g(i,j,3),gdot)

	     		if (gdot.eq.0) then

	     theta=pi/2

	     		else

	     theta=dacos(gdot/(g1norm*g2norm))

	     		end if

	   	end if

	  end if




	      if (enorm1.lt.gpqnorm.and.enorm2.lt.gpqnorm) then

	      singsegm=singsegm+1

	      end if




	end if









c
c
c	d1, d2 (i.e. s1, s2) are the distances of P'' (projection
c	of P' on Gpq) from the segments' Gpq vertices.
c	P'' is taken
c	as the origin of a 1-dimensional local coordinate system
c	on Gpq.
c	f1, f2 (i.e. r1, r2) are the 3d-distances of the vertices 
c	of Gpq from field point P.
c
c						



	if (j.eq.noe(i)) then

	d1=snorm(xxx(i,j)-x(k(i,j)),yyy(i,j)-y(k(i,j)),
     $zzz(i,j)-z(k(i,j)))
	d2=snorm(xxx(i,j)-x(k(i,1)),yyy(i,j)-y(k(i,1)),
     $zzz(i,j)-z(k(i,1)))
	f1=snorm(x(k(i,j)),y(k(i,j)),z(k(i,j)))
	f2=snorm(x(k(i,1)),y(k(i,1)),z(k(i,1)))

	if (abs(d1-f1).lt.1e-10.and.abs(d2-f2).lt.1e-10) then
	go to 1968
	end if


	else 

	d1=snorm(xxx(i,j)-x(k(i,j)),yyy(i,j)-y(k(i,j)),
     $zzz(i,j)-z(k(i,j)))
	d2=snorm(xxx(i,j)-x(k(i,j+1)),yyy(i,j)-y(k(i,j+1)),
     $zzz(i,j)-z(k(i,j+1)))
	f1=snorm(x(k(i,j)),y(k(i,j)),z(k(i,j)))
	f2=snorm(x(k(i,j+1)),y(k(i,j+1)),z(k(i,j+1)))


	if (abs(d1-f1).lt.1e-10.and.abs(d2-f2).lt.1e-10) then
	go to 1968
	end if



	end if



c	
c	check if P'' is located inside the segment Gpq
c

	if (d1.lt.gpqnorm.and.d2.lt.gpqnorm) then
	s1=-d1
	s2=d2
	r1=f1
	r2=f2
	go to 1972
	end if

	if (d1.lt.d2) then

	s1=d1
	s2=d2
	r1=f1
	r2=f2

	else if (d2.lt.d1) then

	s1=-d1
	s2=-d2
	r1=f1
	r2=f2

	end if


	go to 1972






1968	if (d1.lt.d2) then

	s1=d1
	s2=d2
	r1=f1
	r2=f2

	else if (d2.lt.d1) then

	s1=-d1
	s2=-d2
	r1=-f1
	r2=-f2

	else if (d1.eq.d2) then

	s1=-d1
	s2=d2
	r1=-f1
	r2=f2

	end if


1972	continue


	if (h(i).eq.0) then

	   if (thisvert.gt.0) then
	   c1=0
	   c2=0
	   go to 1973
	   end if


	if ((s1+s2).lt.1e-10.and.(r1+r2).lt.1e-10) then
	c1=0
	c2=0
	go to 1973
	end if

	c1=fln(s2,r2,s1,r1)
	c2=0

	
	go to 1973

	end if


	if (thisvert.gt.0) then
	c1=0
	else
	c1=fln(s2,r2,s1,r1)
	end if

	if (hh(i,j).eq.0) then
	c2=0
	else
	c2=fan(h(i),hh(i,j),s1,s2,r1,r2)
	end if

1973	continue

	sum1=sum1+ss(i,j)*hh(i,j)*c1
	sum2=sum2+ss(i,j)*c2



200	continue


c
c	If P' lays inside the polygon Ap then
c	ALL of ss(i,j) = 1. In that case sing = noe(i)
c	and the following correction has to be applied
c
c	Ap = lim (Ap) = Ap - 2*pi*h(i)
c	       R->0
c
c	If P' lays on an edge Gpq the above correction term is
c	
c	                   - pi*h(i)
c
c	and finally when P' is located on a vertex of Ap we get
c
c	                   - theta*h(i) 
c	where
c
c	theta = arccos{ (g1 (dot) g2) / (|g1||g2|) }, g1 and g2 are
c	the two vectors Gpq meeting at the vertex.
c
c	theta has already been computed (see above)
c


	if (singarea.eq.noe(i)) then
	area=2*pi*h(i)
	else
	area=0
	end if

	if (singsegm.gt.0) then
	edge=pi*h(i)
	else 
	edge=0
	end if

	if (singvert.gt.0) then
	vertex=theta*h(i)
	else
	vertex=0
	end if


c     
c     Final computation of Equations (10) and (11)
c


	v=v+s(i)*h(i)*(sum1+h(i)*sum2-area-edge-vertex)
	vx=vx+cs(i,1)*(sum1+h(i)*sum2-area-edge-vertex)
	vy=vy+cs(i,2)*(sum1+h(i)*sum2-area-edge-vertex)
	vz=vz+cs(i,3)*(sum1+h(i)*sum2-area-edge-vertex)
    



100	continue

      T_g(p_j,1)=v
      T_g(p_j,2)=vx
      T_g(p_j,3)=vy
      T_g(p_j,4)=vz
c
c
c	Computation of Vij
c
c
	do 901 i=1,3
	do 901 j=1,3
	vij(i,j)=0
901	continue

	do 400 i=1,3
	do 500 j=1,3


		do 600 m=1,nop

		sum1=0
		sum2=0

		singarea=0
		singsegm=0
		singvert=0



		do 700 n=1,noe(m)

		thisvert=0

		gpqnorm=snorm(g(m,n,1),g(m,n,2),g(m,n,3))

c	
c		check if P' lays inside the polygon Np(m)
c

		if (ss(m,n).eq.1) then
		singarea=singarea+1
		end if



c
c		check if P' lays inside the segment Gpq
c
c		or
c
c		on the vertex (m,n)
c
c		Condition controlled is  
c
c		ss(m,n) = 0
c

		if (ss(m,n).eq.0) then

		if (n.eq.noe(m)) then 

		gx1=xx(m)-x(k(m,1))
		gx2=xx(m)-x(k(m,n))
		gy1=yy(m)-y(k(m,1))
		gy2=yy(m)-y(k(m,n))
		gz1=zz(m)-z(k(m,1))
		gz2=zz(m)-z(k(m,n))

		else  

		gx2=xx(m)-x(k(m,n))
		gx1=xx(m)-x(k(m,n+1))
		gy2=yy(m)-y(k(m,n))
		gy1=yy(m)-y(k(m,n+1))
		gz2=zz(m)-z(k(m,n))
		gz1=zz(m)-z(k(m,n+1))

		end if

		enorm1=snorm(gx1,gy1,gz1)
		enorm2=snorm(gx2,gy2,gz2)



c
c		check if P' lays on the Vertex (m,n) 
c

		if (enorm1.eq.0) then

		singvert=singvert+1

		thisvert=thisvert+1

		if (n.eq.noe(m)) then

		g1norm=snorm(g(m,n,1),g(m,n,2),g(m,n,3))
		g2norm=snorm(g(m,1,1),g(m,1,2),g(m,1,3))

	


	call dot(-g(m,n,1),-g(m,n,2),-g(m,n,3),g(m,1,1),g(m,1,2),
     $g(m,1,3),gdot)

			if (gdot.eq.0) then

			theta=pi/2

			else

			theta=dacos(gdot/(g1norm*g2norm))

			end if

		else

		g1norm=snorm(g(m,n,1),g(m,n,2),g(m,n,3))
		g2norm=snorm(g(m,n+1,1),g(m,n+1,2),g(m,n+1,3))




	call dot(-g(m,n,1),-g(m,n,2),-g(m,n,3),g(m,n+1,1),g(m,n+1,2),
     $g(m,n+1,3),gdot)

			if (gdot.eq.0) then

			theta=pi/2

			else

			theta=dacos(gdot/(g1norm*g2norm))

			end if

		end if

		end if



		if (enorm2.eq.0) then

		singvert=singvert+1

		thisvert=thisvert+1

		if (n.eq.1) then

	g1norm=snorm(g(m,noe(m),1),g(m,noe(m),2),g(m,noe(m),3))
	g2norm=snorm(g(m,1,1),g(m,1,2),g(m,1,3))


	call dot(-g(m,noe(m),1),-g(m,noe(m),2),-g(m,noe(m),3),
     $g(m,1,1),g(m,1,2),g(m,1,3),gdot)

			if (gdot.eq.0) then

			theta=pi/2

			else

			theta=dacos(gdot/(g1norm*g2norm))

			end if

		else

	g1norm=snorm(g(m,n-1,1),g(m,n-1,2),g(m,n-1,3))
	g2norm=snorm(g(m,n,1),g(m,n,2),g(m,n,3))


	call dot(-g(m,n-1,1),-g(m,n-1,2),-g(m,n-1,3),g(m,n,1),g(m,n,2),
     $g(m,n,3),gdot)

			if (gdot.eq.0) then

			theta=pi/2

			else

			theta=dacos(gdot/(g1norm*g2norm))

			end if

		end if

		end if

		if (enorm1.lt.gpqnorm.and.enorm2.lt.gpqnorm) then

		singsegm=singsegm+1

		end if

		end if





c
c
c	d1, d2 (i.e. s1, s2) are the distances of P'' (projection
c	of P' on Gpq) from the segments' Gpq vertices.
c	P'' is taken
c	as the origin of a 1-dimensional local coordinate system
c	on Gpq.
c	f1, f2 (i.e. r1, r2) are the 3d-distances of the vertices 
c	of Gpq from filed point P.
c
c



		if (n.eq.noe(m)) then

		d1=snorm(xxx(m,n)-x(k(m,n)),yyy(m,n)-
     $y(k(m,n)),
     $zzz(m,n)-z(k(m,n)))
		d2=snorm(xxx(m,n)-x(k(m,1)),yyy(m,n)-
     $y(k(m,1)),
     $zzz(m,n)-z(k(m,1)))
		f1=snorm(x(k(m,n)),y(k(m,n)),z(k(m,n)))
		f2=snorm(x(k(m,1)),y(k(m,1)),z(k(m,1)))


	if (abs(d1-f1).lt.1e-10.and.abs(d2-f2).lt.1e-10) then
	go to 1958
	end if

		else

		d1=snorm(xxx(m,n)-x(k(m,n)),yyy(m,n)-
     $y(k(m,n)),
     $zzz(m,n)-z(k(m,n)))
		d2=snorm(xxx(m,n)-x(k(m,n+1)),
     $yyy(m,n)-y(k(m,n+1)),zzz(m,n)-z(k(m,n+1)))
		f1=snorm(x(k(m,n)),y(k(m,n)),z(k(m,n)))
		f2=snorm(x(k(m,n+1)),y(k(m,n+1)),z(k(m,n+1)))


	if (abs(d1-f1).lt.1e-10.and.abs(d2-f2).lt.1e-10) then
	go to 1958
	end if

		end if



c	
c	check if P'' is located inside the segment Gpq
c

	if (d1.lt.gpqnorm.and.d2.lt.gpqnorm) then
	s1=-d1
	s2=d2
	r1=f1
	r2=f2
	go to 1982
	end if

c
c	1-D Coordinate system orientation
c

	if (d1.lt.d2) then

	s1=d1
	s2=d2
	r1=f1
	r2=f2

	else if (d2.lt.d1) then

	s1=-d1
	s2=-d2
	r1=f1
	r2=f2

	end if

	go to 1982

1958	if (d1.lt.d2) then

	s1=d1
	s2=d2
	r1=f1
	r2=f2

	else if (d2.lt.d1) then

	s1=-d1
	s2=-d2
	r1=-f1
	r2=-f2

	end if


1982	continue

	if (h(m).eq.0) then

	   if (thisvert.gt.0) then
	   c1=0
	   c2=0
	   go to 1983
	   end if

	c1=fln(s2,r2,s1,r1)
	c2=0

	go to 1983

	end if


	if (thisvert.gt.0) then
	c1=0
	else
	c1=fln(s2,r2,s1,r1)
	end if

	if (hh(m,n).eq.0) then
	c2=0
	else
	c2=fan(h(m),hh(m,n),s1,s2,r1,r2)
	end if

1983	continue




		sum1=sum1+cspq(m,n,j)*c1
		sum2=sum2+ss(m,n)*c2



700	continue



c
c	If P' lays inside the polygon Ap then
c	ALL of ss(m,n) = 1. In that case sing = noe(m)
c	and the following correction has to be applied
c
c	Ap = lim (Ap) = Ap - 2*pi*h(m)
c	       R->0
c
c	If P' lays on an edge Gpq the above correction term is
c	
c	                   - pi*h(m)
c
c	and finally when P' is located on a vertex of Ap we get
c
c	                   - theta*h(m) 
c	where
c
c	theta = arccos{(g1 (dot) g2) / (|g1||g2|)}, g1 and g2 are
c	the two vectors Gpq meeting at the vertex.
c


	if (singarea.eq.noe(m)) then
	area=2*pi
	else
	area=0
	end if

	if (singsegm.gt.0) then
	edge=pi
	else 
	edge=0
	end if

	if (singvert.gt.0) then
	vertex=theta
	else
	vertex=0
	end if


c     
c     Final computation of Equations (12)
c

	
	vij(i,j)=vij(i,j)+cs(m,i)*(sum1+s(m)*cs(m,j)*(sum2-area
     $-edge-vertex))



600	continue

500	continue
400	continue
  
        
        T_g(p_j,5)=vij(1,1)
		T_g(p_j,6)=vij(1,2)
		T_g(p_j,7)=vij(1,3)
		T_g(p_j,8)=vij(2,2)
		T_g(p_j,9)=vij(2,3)
		T_g(p_j,10)=vij(3,3)
2100   continue 
        	
      close(11)
      close(12)
      close(13)
	  close(20)

1001	format(3i4)
1002	format(7(4i2)/)
      return
	end


c
c
c	subroutines
c
c
	subroutine planedist (x,y,z,kp1,kp2,kp3,hp,xax,yax,zax)
c
c	computes the distance of P(xp,yp,zp) from 
c	the plane a*x+b*y+c*z+d=0
c	for our case (xp=yp=zp=0) it is
c	h = d / sqrt(a^2+b^2+c^2)
c
	implicit real*8 (a-h,o-z)
	dimension x(8),y(8),z(8)
	a=(y(kp2)-y(kp1))*(z(kp3)-z(kp1))-
     $(y(kp3)-y(kp1))*(z(kp2)-z(kp1))
	b=(z(kp2)-z(kp1))*(x(kp3)-x(kp1))-
     $(z(kp3)-z(kp1))*(x(kp2)-x(kp1))
	c=(x(kp2)-x(kp1))*(y(kp3)-y(kp1))-
     $(x(kp3)-x(kp1))*(y(kp2)-y(kp1))
	d=-a*x(kp1)-b*y(kp1)-c*z(kp1)
	hp=abs(d/dsqrt(a*a+b*b+c*c))

	if (z(kp1).eq.z(kp2).and.z(kp2).eq.z(kp3)) then
	xax=0
	yax=0
	zax=-d/c
	else if (y(kp1).eq.y(kp2).and.y(kp2).eq.y(kp3)) then
	xax=0
	yax=-d/b
	zax=0
	else if (x(kp1).eq.x(kp2).and.x(kp2).eq.x(kp3)) then
	xax=-d/a
	yax=0
	zax=0
	end if

	if (a.eq.0) then
	xax=0
	else
	xax=-d/a
	end if

	if (b.eq.0) then
	yax=0
	else
	yax=-d/b
	end if

	if (c.eq.0) then
	zax=0
	else
	zax=-d/c
	end if



	return
	end


	subroutine segmproj(x1,y1,z1,x2,y2,z2,xo,yo,zo,xs,ys,zs)
c
c	Computes the projections of P'(xx(i),yy(i),zz(i)) from
c	each line segment lpq, i.e. the coordinate values
c	P''(xxx(i,j),yyy(i,j),zzz(i,j)). P'' is considered to lie
c	on the intersection of 3 planes (see notes), namely
c	plane XoX1X2, plane XX1X2 and the plane perpendicular to
c	vector X1X2. Building the respective plane equations, we get
c	a 3x3 system of equations whose solution gives the 
c	vector xxx. The equations (23) - (25) build this system.
c
	implicit real*8 (a-h,o-z)
	r1=x2-x1
	r2=y2-y1
	r3=z2-z1
	e1=x1-xo
	e2=y1-yo
	e3=z1-zo
	a1=x2-x1
	b1=y2-y1
	c1=z2-z1
	d1=a1*xo+b1*yo+c1*zo
 	call cross (e1,e2,e3,r1,r2,r3,a2,b2,c2)	
	call cross (a2,b2,c2,r1,r2,r3,a3,b3,c3)
	call dot (a2,b2,c2,xo,yo,zo,d2)
	call dot (a3,b3,c3,x1,y1,z1,d3)
	det=a1*b2*c3+b1*c2*a3+c1*a2*b3-a3*b2*c1-b3*c2*a1-c3*a2*b1
	dx=d1*b2*c3+b1*c2*d3+c1*d2*b3-d3*b2*c1-b3*c2*d1-c3*d2*b1
	dy=a1*d2*c3+d1*c2*a3+c1*a2*d3-a3*d2*c1-d3*c2*a1-c3*a2*d1
	dz=a1*b2*d3+b1*d2*a3+d1*a2*b3-a3*b2*d1-b3*d2*a1-d3*a2*b1
	xs=dx/det
	ys=dy/det
	zs=dz/det
	return
	end




	subroutine cross (a1,a2,a3,b1,b2,b3,cr1,cr2,cr3)
c
c	Cross product A x B of vectors A=(A1,A2,A3)u and
c	B=(B1,B2,B3)u
c
	implicit real*8 (a-h,o-z)
	cr1=a2*b3-a3*b2
	cr2=a3*b1-a1*b3
	cr3=a1*b2-a2*b1
	return
	end


	subroutine dot (a1,a2,a3,b1,b2,b3,cdot)
c
c	Dot product A . B of vectors A=(A1,A2,A3)u and
c	B=(B1,B2,B3)u
c
	implicit real*8 (a-h,o-z)
	cdot=a1*b1+a2*b2+a3*b3
	return
	end

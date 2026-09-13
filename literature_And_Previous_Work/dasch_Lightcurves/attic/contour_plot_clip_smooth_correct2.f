C Copyright the President and Fellows of Harvard College.
C Licensed under the MIT License

C Make 2D contour/greyscale plot (mean[Z]) from irregularly spaced X,Y,Z data points
C 1. Bin the data into a grid (bin_value=average of member values)
C 2. vizualize the array
C This version contour_plot_clip.f uses an iterative sigma-clipping routine to 
C derive mean, median & sigma. Default settings are 3 iterations, removing 3sigma outliers
C Parameter "clip" specifies number of clipping iterations. If non integer (ie. M.N), 
C the part following the decimal point is treated as a stopping condition, and 
C the clipping-code will keep iterating until fractional improvement in sigma is less than 0.N
C In contour_plot_clip_smooth.f, a weighting function is used to sample data from an
C arbitrarily sized "boxcar" onto a finely spaced grid. -thus smoothing out fluctuations
C and interpolating-over sparse areas. The smoothing scale is set by parameter "scale"
C which is specified in the same units as X & Y.
C 16Jul2007: Correct the input Z values using the smooth clipped-median error model.
C
C Jan 29, 2008: Edward J. Los  Increase the size of ID to avoid truncation.
C                              Correct calculation of sigma1
C                              Initialize the stop variable.
C                               
 
	character*100 infile,xtext,ytext,plotdev,outfile
	integer lim,lim2,lim3,nx,ny,ncon
	integer nvalues,goodvals,goodval2
	parameter (lim=1e7,lim2=1000,lim3=100)
        character*16 ID(lim)

	integer CHECK(lim),NPOUT(lim)
	real X(lim),Y(lim),Z(lim),TR(6),CON(lim3)	
	real MAP(lim2,lim2),TMP(lim),ZOUT(lim),ERROUT(lim)
	real AMED(lim2,lim2),AMEAN(lim2,lim2),APOINT(lim2,lim2),ASIGMA(lim2,lim2)
	real x1,x2,y1,y2,dx,dy
	real xmin,xmax,ymin,ymax,temp
	real medmin,medmax,sigmax,sigmin,npmin,npmax,dz
	real ave,med,rms, clip,scale,dist	
	real mapmin,mapmax

	read(*,'(A100)') infile
        read(*,*) nx,ny
	read(*,'(A100)') xtext
	read(*,'(A100)') ytext
	read(*,*) ncon
	read(*,'(A100)') plotdev
	read(*,*) clip
	read(*,*) scale
 	read(*,'(A100)') outfile

C Process the X and Y limits
	xmin= 1.0e6
 	xmax= -1.0e6
	ymin= 1.0e6
	ymax= -1.0e6

	open(1,file=infile,status='old')

	do i=1,lim
	  read(1,*,end=10,err=10) ID(i),X(i),Y(i),Z(i)
	  if(X(i).gt.xmax) xmax=X(i)
	  if(Y(i).gt.ymax) ymax=Y(i)
	  if(X(i).lt.xmin) xmin=X(i)
	  if(Y(i).lt.ymin) ymin=Y(i)
	enddo
10	np=i-1
        write(*,*) np,' points read'
	dx=(xmax-xmin)/real(nx)
	dy=(ymax-ymin)/real(ny)

C Slightly extend the X,Y limits to allow for co-ordinate round-off errors
	xmin=xmin-dx/10.0
	ymin=ymin-dy/10.0
	xmax=xmax+dx/10.0
        ymax=ymax+dy/10.0
	dx=(xmax-xmin)/real(nx)
        dy=(ymax-ymin)/real(ny)



C Report the input values for checking
        write(*,*) 'infile=',infile
        write(*,*) 'xmin=',xmin,' xmax=',xmax
        write(*,*) 'ymin=',ymin,' ymax=',ymax
        write(*,*) 'nx=',nx,' ny=',ny
 	write(*,*) 'grid-spacing: X=',dx,' Y=',dy
	write(*,*) 'smoothing scale=',scale
        write(*,*) 'xtext=',xtext,' ytext=',ytext,' ncon=',ncon
        write(*,*) 'plotdev=',plotdev
	write(*,*) 'clip=',clip	


C Tranformation Vector for plotting.
	TR(1)=xmin-dx/2
        TR(2)=dx
        TR(3)=0
        TR(4)=ymin-dy/2
        TR(5)=0
        TR(6)=dy


C Initialize the 2D Arrays to zero
	do i=1,lim2
		do j=1,lim2
		  AMED(i,j)=0.0
		  AMEAN(i,j)=0.0
		  ASIGMA(i,j)=0.0
		  APOINT(i,j)=0.0
		  MAP(i,j)=0.0
		enddo
	enddo
	do i=1,lim
		TMP(i)=0.0
		CHECK(i)=0
	enddo
		
C Sample (bin) the points, using a boxcar of radius "scale".
	write(*,*) 'Creating map of input data...'
	do i=1,nx
		x1=xmin+(real(i-1)*dx)
		do j=1,ny
			y1=ymin+(real(j-1)*dy)
			nvalues=0
			do k=1,np
				dist= sqrt((x1-X(k))**2 + (y1-Y(k))**2)
				if(dist.lt.scale) then

c				   print*,'i=',i,' j=',j,' k=',k,' dist=',dist,' X Y Z =',X(k),Y(k),Z(k)

					nvalues=nvalues+1
                                        TMP(nvalues)=Z(k)
				endif		
			enddo
c                       print*,nvalues,3.0,clip,ave,med,rms
			ave=0.0
			med=0.0
			rms=0.0
		        call sigmaclip (TMP,nvalues,goodvals,3.0,clip,ave,med,rms)  
			APOINT(i,j)=real(goodvals)
			AMED(i,j)=med
			ASIGMA(i,j)=rms
			AMEAN(i,j)=ave
c			print*,i,j,nvalues,goodvals,med,rms,ave
		enddo
	enddo

C Correct the input catalog values Z(i) --> ZOUT(i)
C use the smooth-clipped-median gridpoint closest to each actual datapoint.
C a more sophisticated treatment would interpolate from the grid. 
C However if the grid is fine enough there will be no significant difference.
C on this second run through the grid, we can replace bad cells, 
C If any grid-cell had no good points, replace its MEDIAN, MEAN & SIGMA values 
C with the average of the neighboring cells' values. Do not modify the APOINT array
C which provides a record of the quality of each grid point.

        print*,'Correcting input data...'
	do i=1,nx
		x1=xmin+(real(i-1)*dx)
		x2=x1+dx
		do j=1,ny
			y1=ymin+(real(j-1)*dy)
                       	y2=y1+dy
			if(APOINT(i,j).eq.0.0) then
				AMED(i,j)=(AMED(i-1,j)+AMED(i+1,j)+AMED(i,j-1)+AMED(i,j+1))/4.0
				AMEAN(i,j)=(AMEAN(i-1,j)+AMEAN(i+1,j)+AMEAN(i,j-1)+AMEAN(i,j+1))/4.0
				ASIGMA(i,j)=(ASIGMA(i-1,j)+ASIGMA(i+1,j)+ASIGMA(i,j-1)+ASIGMA(i,j+1))/4.0
			endif
			do k=1,np
				if(X(k).ge.x1.and.X(k).lt.x2.and.Y(k).ge.y1.and.Y(k).lt.y2) then
				  ZOUT(k)=AMED(i,j)
				  ERROUT(k)=ASIGMA(i,j)
				  NPOUT(k)=APOINT(i,j)
				  CHECK(k)=CHECK(k) + 1
				endif
			enddo
		enddo
	enddo
C Read back the CHECK array to make sure all points have been corrected exactly once.
	do i=1,np
c	        print*,i,ZOUT(i),ERROUT(i),NPOUT(i),CHECK(i)

		if(CHECK(i).ne.1) write(*,*) 'ERROR: checksum=',CHECK(i),' for datapoint ',i
	enddo


C Create a map of the corrected data, for comparison with the input data. Plot them with identical scaling.
	write(*,*) 'Creating map of corrected data...'
	 do i=1,lim
		TMP(i)=0.0	
	 enddo
	 do i=1,nx
                x1=xmin+(real(i-1)*dx)
                do j=1,ny
                        y1=ymin+(real(j-1)*dy)
                        nvalues=0
                        do k=1,np
                                dist= sqrt((x1-X(k))**2 + (y1-Y(k))**2)
                                if(dist.lt.scale) then
                                        nvalues=nvalues+1
                                        TMP(nvalues)=ZOUT(k)
                                endif
                        enddo
                        med=0.0
           call sigmaclip(TMP,nvalues,goodval2,3.0,clip,ave,med,rms)
                        MAP(i,j)=med
c			print*,i,j,med
                enddo
        enddo

C Write out the corrected data to a file
	open (unit=2,file=outfile,status='unknown') 
	write(*,*) 'Writing corrected data to file...'
	do i=1,np
		write(2,100) ID(i),X(i),Y(i),Z(i),ZOUT(i),ERROUT(i),NPOUT(i)
	enddo
100	format (1X,A20,4X,F12.8,4X,F12.8,4X,F10.6,4X,F9.6,4X,F9.6,4X,I9)

C Find high and low values for plotting.
	medmin=1.0e6
	medmax=-1.0e6
	sigmin=1.0e6
	sigmax=-1.0e6
	npmin=1.0e6
	npmax=-1.0e6
	mapmin=1.0e6
	mapmax=-1.0e6
	do i=1,nx
                do j=1,ny
	                if(AMED(i,j).gt.medmax) medmax=AMED(i,j)
			if(AMED(i,j).lt.medmin) medmin=AMED(i,j) 
	                if(ASIGMA(i,j).gt.sigmax) sigmax=ASIGMA(i,j)
			if(ASIGMA(i,j).lt.sigmin) sigmin=ASIGMA(i,j) 
	                if(APOINT(i,j).gt.npmax) npmax=APOINT(i,j)
			if(APOINT(i,j).lt.npmin) npmin=APOINT(i,j) 
			if(MAP(i,j).gt.mapmax) mapmax=MAP(i,j)
                        if(MAP(i,j).lt.mapmin) mapmin=MAP(i,j)			 
               enddo
        enddo 
	print*,'mapmin',mapmin,'mapmax',mapmax
	print*,'medmin',medmin,'medmax',medmax
c	mapmin=medmin
c	mapmax=medmax

C set contour levels
        if(ncon.gt.1) then
	  dz=(medmax-medmin)/real(ncon)
	  write(*,*) "Contour Levels for input data-map"
	  do i=1,ncon
		CON(i)=medmin+(real(i)*dz)	
		write(*,*)CON(i)
	  enddo
	  write(*,*) "------------------------"
	endif


C PLOT THE RESULTS!
        xmin= xmin-dx/2.
        ymin= ymin-dy/2.
        xmax= xmax+dx/2.
        ymax= ymax+dy/2.
        write (*,*) xmin,ymin,xmax,ymax
	write (*,*) medmax, medmin
C open graphics device:
        CALL PGOPEN(plotdev)
        call pgslw (4)
        call pgsch (1.5)

C subdivide the page into 4 panels
	call PGSUBP(2,2)

C Panel(1) : Clipped Median 
	call PGPANL(2,2)
        call PGENV(xmin, xmax, ymin, ymax,  0,  0)
	call PGMTXT ('B',2.5,0.5,0.5,xtext)
	call PGMTXT ('L',2.0,0.5,0.5,ytext)
	call PGMTXT ('T',3.0,0.0,0.0,plotdev)
	call PGMTXT ('T',1.0,0.5,0.5,'Clipped Median Map of Dataset')
	call PGMTXT ('T',1.0,1.0,0.0,'Median')
	call PGGRAY (AMED,lim2,lim2,1,nx,1,ny,medmax,medmin,TR)
        call PGWEDG('RG', 1.0, 4.0, medmax, medmin, '')
        call PGSLS(1)
        if(ncon.gt.1) then
           ncon = ncon * -1
           write(*,*) 'ncon=',ncon 
           call PGCONT (AMED,lim2,lim2,1,nx,1,ny,CON,ncon,TR)
        endif

C Panel(2) : Clipped SIGMA 
	call PGPANL(1,2)
        call PGENV(xmin, xmax, ymin, ymax,  0,  0)
        call PGMTXT ('B',2.5,0.5,0.5,xtext)
        call PGMTXT ('L',2.0,0.5,0.5,ytext)
	call PGMTXT ('T',1.0,0.5,0.5,'Scatter of Calibration Points')
        call PGMTXT ('T',1.0,1.0,0.0,'Sigma')
        call PGGRAY (ASIGMA,lim2,lim2,1,nx,1,ny,sigmax,sigmin,TR)
        call PGWEDG('RG', 1.0, 4.0, sigmax, sigmin, '')

C Panel(3) : Number of good points contributing to each grid-cell
        call PGPANL(2,1)
        call PGENV(xmin, xmax, ymin, ymax,  0,  0)
        call PGMTXT ('B',2.5,0.5,0.5,xtext)
        call PGMTXT ('L',2.0,0.5,0.5,ytext)
	call PGMTXT ('T',1.0,0.5,0.5,'Number of Points Contributing to each Grid-Cell')         
	call PGMTXT ('T',1.0,1.0,-1.5,'N')
        call PGGRAY (APOINT,lim2,lim2,1,nx,1,ny,npmax,npmin,TR)
        call PGWEDG('RG', 1.0, 4.0, npmax, npmin, '')

C Panel(4) : Map of corrected datapoints, Median
        call PGPANL(1,1)
        call PGENV(xmin, xmax, ymin, ymax,  0,  0)
        call PGMTXT ('B',2.5,0.5,0.5,xtext)
        call PGMTXT ('L',2.0,0.5,0.5,ytext)
        call PGMTXT ('T',1.0,0.5,0.5,'Corrected Data Map')
        call PGMTXT ('T',1.0,1.0,0.0,'Median')
        call PGGRAY (MAP,lim2,lim2,1,nx,1,ny,mapmax,mapmin,TR)
        call PGWEDG('RG', 1.0, 4.0, mapmax, mapmin, '')


	call PGIDEN

	call pgend
   

	END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C  sigmaclip.f
C  compute the median and RMS of a set of numbers
C  using an iterative sigma clipping, arguments are
C  the Nsigma and Niterations for clipping
C  If Niter is non integer, then the decimal 
C  part is treated as as a stopping condition
C  and the code will iterate until the fractional 
C  improvement in sigma is less that 0.N
C  The routine returns the clipped statistical values
C  and the array X modified to contain only the "good" points.
C23456
        subroutine sigmaclip (x,npoints,ngood,nsigma,
     +		  iter,mean,median,sigma)
	integer npoints,lim,ngood,niter
	real x(npoints),nsigma,sigma1,sigma,mean1,mean,median1,median,sum,sumsq
	real min, max, stop,iter,temp(npoints)

C initialize values of certain variables:
	ngood=npoints
        stop = 0.0
	if (iter-real(int(iter)).gt.0.0) then
		niter = 1000
		stop = iter-real(int(iter))
	else
	  	niter = int(iter)
	endif 

c First pass to get raw mean
	sum = 0.0
	do i = 1,npoints
          sum = sum + x(i)
	enddo
10      npoints = i-1
   	mean1 = sum/real(npoints)
	
c Get the raw Median
	call PIKSRT(npoints,x)  
        median1 = x(npoints/2)

c Second pass to get raw RMS or standard eviation
        sumsq = 0.0
   	do i = 1,npoints
	  sumsq = sumsq + (x(i)-mean1)**2    
	enddo
        sigma1 = sqrt(sumsq/real(npoints-1))

c	print*,'iteration 0',' ngood',ngood,' mean',mean1,' sigma',sigma1
	
C Now iterate the requested number of times, to obtain a cleaned up "clipped" mean and sigma
	sigma = sigma1	
        mean = mean1
	median = median1

	do i = 1,niter

		mean1 = mean
		sigma1 = sigma
		min = mean - real(nsigma)*sigma
		max = mean + real(nsigma)*sigma
		ngood = 0 
		sum = 0.0
		sumsq = 0.0
		do j = 1,npoints
			if(x(j).gt.min.and.x(j).lt.max) then
				sum = sum + x(j)
				ngood = ngood + 1
 				temp(ngood)=x(j)		
			endif					
		enddo
		mean = sum/real(ngood)
                call PIKSRT(ngood,temp)
                median = temp(ngood/2)
		do j = 1,npoints
                        if(x(j).gt.min.and.x(j).lt.max) then
				sumsq = sumsq + (x(j)-mean)**2
                        endif
                enddo
		sigma = sqrt(sumsq/real(ngood-1))
		
c		if ( (mean1-mean).lt.stop ) go to 100 
		if ( (sigma1-sigma).lt.stop ) go to 100 
 	enddo
100	mean=mean
        sigma=sigma
        median=median

c Make sure that weird numbers don't get spat out.
        if(ngood.eq.0)then
		mean=0.0
		sigma=0.0
		median=0.0
	endif

c finally replace the input array X with the array of goodpoints TEMP
	do i=1,npoints
		x(i)=0.0
	enddo
	do i=1,ngood
		x(i)=temp(i)
	enddo

c print helpful de-bugging information
c	print*,'npoints=',npoints,' ngoodpoints=',ngood,' niter=',niter
c	print*,'mean1=',mean1,' mean=',mean
c	print*,'median1=',median1,' median=',median
c	print*,'sigma1=',sigma1,' sigma=',sigma
c	print*,'good points:'
c       print*,x
c	print*,''

        return
	end


!*****************************************************
!* Sorts an array ARR of length N in ascending order *
!* by straight insertion.                            *
!* ------------------------------------------------- *
!* INPUTS:                                           *
!*	    N	  size of table ARR                  *
!*          ARR	  table to be sorted                 *
!* OUTPUT:                                           *
!*	    ARR   table sorted in ascending order    *
!*                                                   *
!* NOTE: Straight insertion is a N² routine and      *
!*       should only be used for relatively small    *
!*       arrays (N<100).                             *
!*****************************************************         
	SUBROUTINE PIKSRT(N,ARR)
  	real ARR(N)
  	do j=2, N
    		a=ARR(j)
    		do i=j-1,1,-1
      			if (ARR(i)<=a) goto 10
      			ARR(i+1)=ARR(i)
    		end do
		i=0
10  		ARR(i+1)=a
  	end do
  	return
	END

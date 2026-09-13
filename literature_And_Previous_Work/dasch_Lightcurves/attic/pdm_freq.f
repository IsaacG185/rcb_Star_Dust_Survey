C Copyright the President and Fellows of Harvard College.
C Licensed under the MIT License

* pdm analysis for the XTE PCA data
* assume data is in TIME,FLUX,ERROR form.
*2345678
*       Sep 15 2007 - Make failures silent to avoid swamping output in batch
*                     jobs.
*                     Remove "** ERROR:  No points in bin in PERIOD_PDM."
*                     Remove ** ERROR:  Zero variance in PERIOD_PDM."
*                     Remove "subroutine PERIOD_PDM fallen over"
*
*                     Remove pause statements.
*                     Correct input format to include the sequence number
*                     Add check for insufficient points
c Mar 24, 2008  Edward J. Los
c    Magnitude has moved from column 3 to column 5
*
	character*50 lightcurve,outname,plate
	integer nbins,mdata,ntrials,lim,iflag,ierror
	parameter(lim=1e6)
	double precision width,TIME(lim),FLUX(lim),PHASE(lim),
     +	XPDM,trial,p1,p2,res,fmean,fsdev,fadev,fvar,error,
     +  RESULTS(lim),PERIODS(lim),tspan,BEST,bestP

	print*,'PDM analysis -input data should be 3 columnns T F Err'
	PRINT*,'SGTL 18/10/2001'
	print*,'this version is linear in frequency space.'
	print*
	PRINT*,'lightcurve filename?'
	read(*,'(A50)') lightcurve
10	open(unit=1,file=lightcurve,status='old',iostat=ierror)
	if(ierror.eq.1)then
		print*,'input lightcurve not found - give alternative filename'
                RETURN
c		read(*,'(A50)') lightcurve
c		go to 10
	end if
	print*,'file open' 
	do i=1,lim
                read(1,*,end=100) error,TIME(i),error,error,FLUX(i),error,error,plate
        enddo
100	mdata=i-1
	print*,mdata,' points found - spanning',TIME(1),TIME(mdata)
        if(mdata.lt.2) then
		print*,'insufficient input points'
                RETURN
        end if
	tspan=TIME(mdata)-TIME(1)
	print*,tspan,'time units'
* set analysis parameters -period limits, resolution etc.....

	print*,'number of PDM bins per trial period?'
	read*,nbins
*	nbins=10
	width=1.0/DBLE(nbins)
	print*,'period range and frequency resolution?'
	read*,p1,p2,res
	if(p2.gt.tspan/2.0)then
		p2=tspan/2.0
		print*,'Maximum trial period is 1/2*data length: reset to',p2
c		pause
	end if
	f1=1./p2
	f2=1./p1
	ntrials=DINT((f2-f1)/res)
	print*,'nbins',nbins,'width',width,'p1',p1,'p2',p2,'res',res,'ntrials',ntrials
	call PERIOD_MOMENT(FLUX,mdata,fmean,fadev,fsdev,fvar)
	print*,'debug: period_moment completed - back in main prog'
*loop over trial periods - fold and find pdm:
	BEST=100.
	do i=1,ntrials
		trial=1./(f1+((i-1)*res))
c	print*,'trial period:',trial
		PERIODS(i)=1./trial
		call PHASE_FOLD(TIME,PHASE,mdata,trial)	
*	call GRAPH(PHASE,FLUX,mdata,0)
		call PERIOD_PDM(PHASE,FLUX,mdata,nbins,width,fvar,XPDM,lim,iflag)
		if(iflag.eq.1)then
c			print*,'subroutine PERIOD_PDM fallen over'
                        XPDM = 100;
c			pause
		end if
 		RESULTS(i)=XPDM
c		print*,'PDM Statistic:',XPDM
		if(XPDM.lt.BEST)then
			BEST=XPDM
			bestP=trial
		end if
	enddo			
** now plot the results!
	print*,'Best Period:',bestP,' Best PDM: ',BEST,' points: ',mdata
C	call GRAPH(PERIODS,RESULTS,ntrials,1)

*finally write a file
	print*,'name for results file?'
	read(*,'(A50)')outname
	open(unit=2,file=outname,status='new')
	do i=1,ntrials
		write(2,*) PERIODS(i),RESULTS(i)
	enddo
	print*,'*Finished*'
	END

      SUBROUTINE PERIOD_PDM(XDATA, YDATA, NDATA, NBIN, WBIN, VARI, PDM, 
     :                      MAXPTS, IFAIL)
 
C===========================================================================
C This routine calculates the PDM statistic of XDATA(NDATA), YDATA(NDATA),
C following the method described by Stellingwerf (1978, APJ, 224, 953).
C Basically, the data is folded on a trial frequency and split up into
C NBIN samples, each of width WBIN. The PDM statistic is then calculated
C by estimating the mean and variance of each sample, calculating the
C overall variance for all of the samples and then dividing this by the
C variance of the whole dataset. If there is an error IFAIL = 1, otherwise
C IFAIL = 0.
C
C Written by Vikram Singh Dhillon @LPO 3-March-1993.
C===========================================================================
 
      IMPLICIT NONE
 
C-----------------------------------------------------------------------------
C PERIOD_PDM declarations.
C-----------------------------------------------------------------------------
 
      INTEGER NDATA, NBIN, IFAIL, MAXPTS, MAXNPTS, I, J, K
      PARAMETER (MAXNPTS=1e6)
      INTEGER NPTS(MAXNPTS)
      DOUBLE PRECISION XDATA(NDATA), YDATA(NDATA)
      DOUBLE PRECISION WBIN, VARI, PDM
      DOUBLE PRECISION BINWID, CENBIN, MINBIN, MAXBIN
      DOUBLE PRECISION NOM, DENOM
      DOUBLE PRECISION AVE, ADEV, SDEV
      DOUBLE PRECISION SAMPLE(MAXNPTS), SVAR(MAXNPTS)
      double precision MEANS(MAXNPTS)
      BYTE BELL
      DATA BELL/7/
 
C-----------------------------------------------------------------------------
C Test maximum array dimensions.
C-----------------------------------------------------------------------------
 
      IF ( MAXNPTS.NE.MAXPTS ) THEN
         WRITE (*, *) BELL
         WRITE (*, *) '** ERROR: MAXNPTS not equal to MAXPTS.'
         WRITE (*, *) '** ERROR: Change dimension of MAXNPTS'
         WRITE (*, *) '** ERROR: in subroutine PERIOD_PDM.'
         RETURN
      END IF
 
C-----------------------------------------------------------------------------
C Loop through the NBIN samples.
C-----------------------------------------------------------------------------
 
      DO 100 I = 1, NBIN
 
C-----------------------------------------------------------------------------
C Determine the limits of the bin.
C-----------------------------------------------------------------------------
 
         BINWID = 1./DBLE(NBIN)
         CENBIN = (DBLE(I)*BINWID) - (0.5*BINWID)
         MINBIN = CENBIN - (0.5*WBIN)
         MAXBIN = CENBIN + (0.5*WBIN)
 
C-----------------------------------------------------------------------------
C Loop through the phases and store those that fall in the bin.
C-----------------------------------------------------------------------------
 
         NPTS(I) = 0
         DO 50 K = 1, 3
            DO 20 J = 1, NDATA
               IF ( K.EQ.1 ) THEN
                  IF ( XDATA(J)-1..LE.MAXBIN ) THEN
                     IF ( XDATA(J)-1..GT.MINBIN ) THEN
                        NPTS(I) = NPTS(I) + 1
                        SAMPLE(NPTS(I)) = YDATA(J)
                     END IF
                  END IF
               ELSE IF ( K.EQ.2 ) THEN
                  IF ( XDATA(J).LE.MAXBIN ) THEN
                     IF ( XDATA(J).GT.MINBIN ) THEN
                        NPTS(I) = NPTS(I) + 1
                        SAMPLE(NPTS(I)) = YDATA(J)
                     END IF
                  END IF
               ELSE IF ( K.EQ.3 ) THEN
                  IF ( XDATA(J)+1..LE.MAXBIN ) THEN
                     IF ( XDATA(J)+1..GT.MINBIN ) THEN
                        NPTS(I) = NPTS(I) + 1
                        SAMPLE(NPTS(I)) = YDATA(J)
                     END IF
                  END IF
               END IF
 20         CONTINUE
 50      CONTINUE
         
C-----------------------------------------------------------------------------
C Calculate the mean and variance of the sample and store the variance.
C-----------------------------------------------------------------------------
         
         IF ( NPTS(I).GT.1 ) THEN
          CALL PERIOD_MOMENT(SAMPLE, NPTS(I), AVE, ADEV, SDEV, SVAR(I))
          MEANS(I)=AVE 
	ELSE
            SVAR(I) = 0.
         END IF
*	print*,'variance',SVAR(I)
*	print*,'mean',AVE
 100  CONTINUE
* 	call GRAPH(XDATA,MEANS,NBIN,0)
C-----------------------------------------------------------------------------
C Calculate the overall variance for all of the samples and divide it by the
C variance of the whole dataset. The result is the PDM statistic.
C-----------------------------------------------------------------------------
 
      NOM = 0.
      DENOM = 0.
      DO 200 I = 1, NBIN
         NOM = NOM + (DBLE(NPTS(I)-1.)*SVAR(I))
         DENOM = DENOM + DBLE(NPTS(I))
 200  CONTINUE
      IF ( DENOM.EQ.0. ) THEN
c        WRITE (*, *) BELL
c        WRITE (*, *) '** ERROR:  No points in bin in PERIOD_PDM.'
         IFAIL = 1
         GO TO 300
      ELSE IF ( VARI.EQ.0. ) THEN
c        WRITE (*, *) BELL
c        WRITE (*, *) '** ERROR:  Zero variance in PERIOD_PDM.'
         IFAIL = 1
         GO TO 300
      END IF
      DENOM = DENOM - NBIN
      PDM = (NOM/DENOM)/VARI
 
C-----------------------------------------------------------------------------
C And return.
C-----------------------------------------------------------------------------
 
      IFAIL = 0
 300  CONTINUE
      RETURN
      END

C-----------------------------------------------------------------------------

      SUBROUTINE PERIOD_MOMENT(DATA, N, AVE, ADEV, SDEV, VAR)
 
C===============================================================================
C Given an array of DATA of length N, this routine returns its mean AVE,
C average deviation ADEV, standard deviation SDEV and variance VAR.
C
C Adapted from Numerical Recipes by Vikram Singh Dhillon @Sussex 29-April-1992.
C===============================================================================
 
C-------------------------------------------------------------------------------
C PERIOD_MOMENT declarations.
C-------------------------------------------------------------------------------

      INTEGER     N,J
      DOUBLE PRECISION        AVE,ADEV,SDEV,VAR
      DOUBLE PRECISION        S,DATA(N)
      BYTE BELL
      DATA BELL/7/
 
      IF ( N.LE.1 ) THEN
         WRITE (*, *) BELL
         WRITE (*, *) '** ERROR: N must be at least 2 in PERIOD_MOMENT.'
         GO TO 300
      END IF
 
C-------------------------------------------------------------------------------
C First pass to get the mean.
C-------------------------------------------------------------------------------
 
      S = 0.0
      DO 100 J = 1, N
         S = S + DATA(J)
 100  CONTINUE
      AVE = S/N
 
C-------------------------------------------------------------------------------
C Second pass to get the first (absolute) and second moments of the deviation
C from the mean.
C-------
-----------------------------------------------------------------------
 
      ADEV = 0.0
      VAR = 0.0
      DO 200 J = 1, N
         S = DATA(J) - AVE
         ADEV = ADEV +DABS(S)
         VAR = VAR + S*S
 200  CONTINUE
 
C-------------------------------------------------------------------------------
C Put the pieces together according to the conventional definitions.
C-------------------------------------------------------------------------------
 
      ADEV = ADEV/N
      VAR = VAR/DBLE(N-1)
      SDEV = DSQRT(VAR)
 
 300  CONTINUE
      RETURN
      END


**-----------------------------------------------------------------------------


	SUBROUTINE PHASE_FOLD(T,PH,ndata,period)
* simple folding subroutine -transform TIME(NDATA) to PHASE(NDATA)

        double precision T(ndata),PH(ndata),period,zeropt,cycles
        
        zeropt=T(1)
*	print*,'zeropoint:',zeropt,'period:',period
        do i=1,ndata
                cycles=(T(i)-zeropt)/period
                PH(i)=cycles-DBLE(DINT(cycles))
        enddo 
        return
        end   

**----------------------------------------------------------------------------

C	SUBROUTINE GRAPH(DX,DY,n,flag)
C	
C	character*10 display,c1
C	double precision DX(n),DY(n)
C	real X(n),Y(n),ymin,ymax
C	integer flag
C
C	ymin=100.
C	ymax=0.
C	do i=1,n
C		X(i)=SNGL(DX(i))
C		Y(i)=SNGL(DY(i))
CC		if(Y(i).lt.ymin)ymin=Y(i)
C		if(Y(i).gt.ymax)ymax=Y(i)
C	enddo
C
C	display='/xserve'
C10	call PGOPEN(display)
C	call PGENV( X(1),X(n),ymin,ymax,0,0)
C	call PGLAB('Trial Frequency /Hz','PDM Statistic',' ')
CC	call PGPT(n,X,Y,6)
C	if(flag.eq.1)call PGLINE(n,X,Y)
C	call PGCLOS
C	print*,'(H)ardcopy or e(X)it?'
C	read*,c1
C	if(c1.eq.'h')then
C		display='?'
C		go to 10	
CC	end if
C	return
C	END

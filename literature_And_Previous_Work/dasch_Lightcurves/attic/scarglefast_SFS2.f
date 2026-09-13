C Copyright the President and Fellows of Harvard College.
C Licensed under the MIT License

* FAST Lomb-Scargle for large datasets ( >10e6 points)  
* uses Vik Dhillon's PERIOD code -converted into double precision 
* plus pgplot and some of my subroutines

* Silas Laycock 
* Southampton
* July 2001
*2345678

C comments added Nov 6th 2003:
c Input is the lightcurve gx_LR1_bcor_faxb_pcu.txt
c parameters are fixed: period range is P1=0.5sec, P2=1000sec
C resolution is 0.00001 Hz
c these could be read in from the script rather than set in the program
C they are set to the values used in work presented in the big paper.
c This program is identical to scarglefast_3a.f which is for manual use.

* July 21 2004
* PARAMETERS FED IN BY REDUCTION SCRIPT. THEY ARE NO LONGER FIXED.
* VARIANCE AND M INDEPENDANT FREQUENCIES ARE PRINTED OUT

c
c Sep 17, 2007
c    Account for the sequence number in the input file
c    No need to skip the first line of the input file
c    Handle insufficient points gracefully
c Mar 24, 2008  Edward J. Los
c    Magnitude has moved from column 3 to column 5
c

	character*100 lightcurve
        character*100 plate
	integer	ndata,ierror,iINFO,lim1,lim2,nfreq
	parameter(lim1=1e7,lim2=1e7)
	double precision TIME(lim1),FLUX(lim1),FREQ(lim2),POW(lim2)
        double precision minF,maxF,intF,crap,P1,P2,TIMEDEL,mfreq
     	double precision variance, sig,s99,s95,s90
	
c	print*,'          ***** SCARGLEFAST *****'
c	print*,'fast evaluation of the Lomb-Scargle statistic '
c       print*,'using the method described by Press & Rybicki ' 
c       print*,'1989 Astrophysical Journal 338, 277-280.'
c	print*,'Now in Double Precision. Upto 1 million datapoints'
c	print*,'with peak detection. Interactive PGPLOT' 
c	print*,'display with significance estimation'
c	print*,'Silas Laycock          Southampton  July 25th 2001'
c	print*

	

c	print*,'lightcurve?'
	read(*,*) lightcurve
	open(unit=1,file='scarglefast_logfile',status='new')
        open(unit=2,file=lightcurve,status='old',iostat=ierror)
	if(ierror.ne.0)then
                write(1,*)'***error code',ierror,'reading file***'
        end if
	open(unit=4,file='VARMNT',status='unknown')

* skip the first line
c	read(2,*)
	do i=1,lim1
		read(2,*,end=100) crap,TIME(i),crap,crap,FLUX(i),crap,crap,plate
	enddo
100	ndata=i-1
	write(1,*) ndata,'datapoints'
	TIMEDEL=TIME(ndata)-TIME(1)
	write(1,*),'spanning',TIMEDEL,'seconds'

	read(*,*) P1,P2,intF
c 	P1=0.5
c	P2=1000
c	intF=0.00001
	write(1,*),'period range',P1,'-',P2
        write(1,*),'frequency resolution',intF
	iINFO=0
		minF=1.0D+00/P2
		maxF=1.0D+00/P1
	call PERIOD_SCARGLE(TIME,FLUX,ndata,FREQ,POW,lim2,nfreq,ierror,
     +                      minF,maxF,intF,iINFO,variance)     
	
        if (ierror.eq.1) THEN
           WRITE (*, *) '** ERROR: PERIOD_SCARGLE failed'
           RETURN        

        END IF
	open(unit=3,file='scarglefast_pds.txt',status='new')
	do i=1,nfreq
		write(3,*) FREQ(i),POW(i)
	enddo	
	write(1,*),nfreq,'frequencies'
	write(1,*),'PDS written to scarglefast_pds.txt'

c       Estimate number of independant frequencies (M)
	mfreq=2.0D+00*nfreq*intF*TIMEDEL

c       Compute 99% confidence level
        sig=2.0D+00*dble(nfreq)*intF*dble(TIMEDEL)
        s99=DLOG(sig/0.01D+00)
        s95=DLOG(sig/0.05D+00)
        s90=DLOG(sig/0.1D+00)        


	write(4,*) 'variance: ',real(variance)
	write(4,*) 'datapoints: ',ndata	
	write(4,*) 'independent_freq: ',real(mfreq)   
	write(4,*) 'Timespan: ',TIMEDEL
        write(4,*) '99%_level:',s99
        write(4,*) '95%_level:',s95
        write(4,*) '90%_level:',s90
        


	write(1,*),'FINISHED OK'
	END


 
      SUBROUTINE PERIOD_SCARGLE(X, Y, N, WK1, WK2, NWK, NOUT, IFAIL, 
     :                          FMIN, FMAX, FINT, INFO, VAR)
 
C==============================================================================
C Subroutine for fast evaluation of the Lomb-Scargle statistic using the method
C described by Press, W. H. & Rybicki, G. B., 1989, Astrophysical Journal, 338,
C 277-280. Given N data points with abscissas X (which need not be equally
C spaced) and ordinates Y, this routine fills array WK1 with a sequence of
C NOUT increasing frequencies (not angular frequencies) from FMIN up to FMAX
C (with frequency interval FMIN) and fills array WK2 with the values of
C the Lomb-Scargle normalized periodogram at those frequencies. The arrays X and
C Y are not altered. NWK, the dimension of WK1 and WK2, must be large enough
C for intermediate work space, or an error results. If any errors occur, IFAIL
C is set to 1, otherwise IFAIL = 0. If INFO=1 then loop information is output
C to the screen.
C
C Adapted for PERIOD by Vikram Singh Dhillon @Sussex 29-April-1992.
C
C==============================================================================
 
C------------------------------------------------------------------------------
C Number of interpolation points per 1/4 cycle of highest frequency.
C------------------------------------------------------------------------------

      INTEGER MACC 
      PARAMETER (MACC=4)
 
C------------------------------------------------------------------------------
C PERIOD_SCARGLE declarations.
C------------------------------------------------------------------------------

      INTEGER     NOUT,N,NWK,IFAIL,INFO,K,IEL,LOOP
      INTEGER     ISTEP,NDIM,J,NFREQ,NFREQT
      DOUBLE PRECISION        FMIN,FMAX,FINT
      DOUBLE PRECISION        AVE,ADEV,SDEV,VAR 
      DOUBLE PRECISION        X(N),Y(N),WK1(NWK),WK2(NWK)
      DOUBLE PRECISION        XMIN,XMAX,XDIF
      DOUBLE PRECISION        FNDIM,OFAC,HIFAC,FAC,CK,CKK,DEN,DF
      DOUBLE PRECISION        CTERM,STERM,HYPO,HC2WT,HS2WT,CWT,SWT 
      BYTE BELL
      DATA BELL/7/
      DATA ISTEP/50/
 
* 	print*,'in period_scargle'
C------------------------------------------------------------------------------
C Compute the mean, variance, and range of the data.
C------------------------------------------------------------------------------
      
      CALL PERIOD_MOMENT(Y, N, AVE, ADEV, SDEV, VAR)
      IF ( VAR.EQ.0. ) THEN
         WRITE (*, *) BELL
         WRITE (*, *) '** ERROR: Zero variance in PERIOD_SCARGLE.'
         IFAIL = 1
         GO TO 600
      END IF
      IFAIL = 0
      XMIN = X(1)
      XMAX = XMIN
      DO 100 J = 2, N
         IF ( X(J).LT.XMIN ) XMIN = X(J)
         IF ( X(J).GT.XMAX ) XMAX = X(J)
 100  CONTINUE
      XDIF = XMAX - XMIN
 
C------------------------------------------------------------------------------
C Compute OFAC and HIFAC.
C------------------------------------------------------------------------------
 
      OFAC = 1.0D+00/(FINT*XDIF)
      HIFAC = (FMAX/FINT)/(0.5*OFAC*N)
      NOUT = 0.5*OFAC*HIFAC*N
 
 	print*,'OFAC=',OFAC
 	print*,'HIFAC=',HIFAC
 	print*,'no. of Frequencies=',NOUT
C------------------------------------------------------------------------------
C Size the FFT as next power of 2 above NFREQT.
C------------------------------------------------------------------------------
 
      NFREQT = OFAC*HIFAC*N*MACC
 
      NFREQ = 64
 200  CONTINUE
      IF ( NFREQ.LT.NFREQT ) THEN
         NFREQ = NFREQ*2
         GO TO 200
      END IF
      NDIM = 2*NFREQ
      IF ( NDIM.GT.NWK ) THEN
         WRITE (*, *) BELL
         WRITE (*, *) 
     :               '** ERROR: Workspaces too small in PERIOD_SCARGLE.'
         IFAIL = 1
         GO TO 600
      END IF
      IFAIL = 0
 
C------------------------------------------------------------------------------
C Zero the workspaces.
C------------------------------------------------------------------------------
      
      DO 300 J = 1, NDIM
         WK1(J) = 0.
         WK2(J) = 0.
 300  CONTINUE
      FAC = NDIM/(XDIF*OFAC)
      FNDIM = NDIM
 
C------------------------------------------------------------------------------
C Extirpolate the data into the workspaces.
C------------------------------------------------------------------------------
      
      DO 400 J = 1, N
         CK = 1.0d+00 + DMOD((X(J)-XMIN)*FAC, FNDIM)
         CKK = 1.0d+00 + DMOD(2.*(CK-1.), FNDIM)
         CALL PERIOD_SPREAD(Y(J)-AVE, WK1, NDIM, CK, MACC, IFAIL)
         IF ( IFAIL.EQ.1 ) GO TO 600
         CALL PERIOD_SPREAD(1.0d+00, WK2, NDIM, CKK, MACC, IFAIL)
         IF ( IFAIL.EQ.1 ) GO TO 600
 400  CONTINUE
 
C------------------------------------------------------------------------------
C Take the Fast Fourier Transforms.
C------------------------------------------------------------------------------
      
      CALL PERIOD_REALFT(WK1, NFREQ, 1)
      CALL PERIOD_REALFT(WK2, NFREQ, 1)
      DF = 1.0d+00/(XDIF*OFAC)
      K = 3
 
C------------------------------------------------------------------------------
C Compute the Lomb-Scargle value for each frequency between FMIN, FMAX
C with step size FINT.
C------------------------------------------------------------------------------
      
      IEL = 0
      LOOP = 1
      DO 500 J = 1, NOUT
         IF ( (J*DF).LT.FMIN ) THEN
            K = K + 2
         ELSE
            HYPO = DSQRT(WK2(K)**2+WK2(K+1)**2)
            HC2WT = 0.5*WK2(K)/HYPO
            HS2WT = 0.5*WK2(K+1)/HYPO
            CWT = DSQRT(0.5+HC2WT)
            SWT = DSIGN(SQRT(0.5-HC2WT), HS2WT)
            DEN = 0.5*N + HC2WT*WK2(K) + HS2WT*WK2(K+1)
            CTERM = (CWT*WK1(K)+SWT*WK1(K+1))**2/DEN
            STERM = (CWT*WK1(K+1)-SWT*WK1(K))**2/(N-DEN)
            IEL = IEL + 1
            WK1(IEL) = J*DF
            WK2(IEL) = (CTERM+STERM)/(2.*VAR)
            IF ( INFO.EQ.1 ) THEN
               IF ( IEL.EQ.(LOOP*ISTEP) ) THEN
                  WRITE (*, 99001) WK1(IEL), WK2(IEL)
99001             FORMAT ('+Frequency =', D12.6, 
     :                    ',  Lomb-Scargle Statistic =', D12.6)
                  LOOP = LOOP + 1
               END IF
            END IF
            K = K + 2
         END IF
 500  CONTINUE
      NOUT = IEL
 
 600  CONTINUE
      RETURN
      END
 
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
 
* 	print*,'in period_moment'
 
      IF ( N.LE.1 ) THEN
         WRITE (*, *) BELL
         WRITE (*, *) '** ERROR: N must be at least 2 in PERIOD_MOMENT.'
         VAR = 0.0
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
C-------------------------------------------------------------------------------
 
      ADEV = 0.0
      VAR = 0.0
      DO 200 J = 1, N
         S = DATA(J) - AVE
         ADEV = ADEV + ABS(S)
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
 
C==============================================================================
 
      SUBROUTINE PERIOD_SPREAD(Y, YY, N, X, M, IFAIL)
 
C==============================================================================
C Given an array YY of length N, extirpolate (spread) a value Y into M actual
C array elements that best approximate the "fictional" (ie. possibly
C non-integer) array element number X. The weights used are coefficients of
C the Lagrange interpolating polynomial.
C
C Adapted for PERIOD by Vikram Singh Dhillon @Sussex 29-April-1992.
C==============================================================================
 
C------------------------------------------------------------------------------
C PERIOD_SPREAD declarations.
C------------------------------------------------------------------------------

      INTEGER M,IFAIL,J,N,IX,IHI,ILO 
      DOUBLE PRECISION YY(N), NFAC(10)
      BYTE BELL
      DOUBLE PRECISION FAC,NDEN,X,Y
      DATA BELL/7/
      DATA NFAC/1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880/	
 	
      IF ( M.GT.10 ) THEN
         WRITE (*, *) BELL
         WRITE (*, *) 
     :           '** ERROR: Factorial table too small in PERIOD_SPREAD.'
         IFAIL = 1
         GO TO 200
      END IF
      IFAIL = 0
      IX = X
      IF ( X.EQ.DBLE(IX) ) THEN
         YY(IX) = YY(IX) + Y
      ELSE
         ILO = MIN(MAX(IDINT(X-0.5d+00*DBLE(M)+1.0d+00),1), N-M+1)
         IHI = ILO + M - 1
         NDEN = NFAC(M)
         FAC = X - ILO
         DO 50 J = ILO + 1, IHI
            FAC = FAC*(X-J)
 50      CONTINUE
         YY(IHI) = YY(IHI) + Y*FAC/(NDEN*(X-IHI))
         DO 100 J = IHI - 1, ILO, -1
            NDEN = (NDEN/(J+1-ILO))*(J-IHI)
            YY(J) = YY(J) + Y*FAC/(NDEN*(X-J))
 100     CONTINUE
      END IF
 
 200  CONTINUE
      RETURN
      END
 
      SUBROUTINE PERIOD_REALFT(DATA, N, ISIGN)
 
C=============================================================================
C Calculate the Fourier Transform of a set of 2N real-valued data points.
C Replaces this data (which is stored in array DATA) by the positive frequency
C half of its complex Fourier Transform. The real-valued first and last
C components of the complex transform are returned as elements DATA(1) and
C DATA(2) respectively. N must be a power of 2. This routine also calculates
C the inverse transform of a complex data array if it is the transform of
C real data. (Result in this case must be multipled by 1/N.)
C
C Adapted from Numerical Recipes by Vikram Singh Dhillon @Sussex 1-July-1992.
C=============================================================================
 
C-----------------------------------------------------------------------------
C PERIOD_REALFT declarations. Double precision for ALL trigonometric variables.
C-----------------------------------------------------------------------------

      INTEGER          N,N2P3,I1,I2,I3,I4,ISIGN,I
      DOUBLE PRECISION WR, WI, WPR, WPI, WTEMP, THETA,C1,C2
      DOUBLE PRECISION             DATA(2*N)
      DOUBLE PRECISION WRS,WIS,H1R,H2R,H1I,H2I
     
C-----------------------------------------------------------------------------
C Initialize the recurrence.
C-----------------------------------------------------------------------------
 
      THETA = 3.141592653589793D0/DBLE(N)
      C1 = 0.5
      IF ( ISIGN.EQ.1 ) THEN
         C2 = -0.5
 
C-----------------------------------------------------------------------------
C The forward transform is here.
C-----------------------------------------------------------------------------
 
         CALL PERIOD_FOUR1(DATA, N, +1)
      ELSE
 
C-----------------------------------------------------------------------------
C Otherwise set up for an inverse transform.
C-----------------------------------------------------------------------------
 
         C2 = 0.5
         THETA = -THETA
      END IF
      WPR = -2.0D0*DSIN(0.5D0*THETA)**2
      WPI = DSIN(THETA)
      WR = 1.0D0 + WPR
      WI = WPI
      N2P3 = 2*N + 3
 
C-----------------------------------------------------------------------------
C Case I = 1 done separately below.
C-----------------------------------------------------------------------------
 
      DO 100 I = 2, N/2 + 1
         I1 = 2*I - 1
         I2 = I1 + 1
         I3 = N2P3 - I2
         I4 = I3 + 1
         WRS = DBLE(WR)
         WIS = DBLE(WI)
 
C-----------------------------------------------------------------------------
C The two separate transforms are separated out of Z.
C-----------------------------------------------------------------------------
 
         H1R = C1*(DATA(I1)+DATA(I3))
         H1I = C1*(DATA(I2)-DATA(I4))
         H2R = -C2*(DATA(I2)+DATA(I4))
         H2I = C2*(DATA(I1)-DATA(I3))
 
C-----------------------------------------------------------------------------
C Here they are recombined to form the true transform of the original
C real data.
C-----------------------------------------------------------------------------
 
         DATA(I1) = H1R + WRS*H2R - WIS*H2I
         DATA(I2) = H1I + WRS*H2I + WIS*H2R
         DATA(I3) = H1R - WRS*H2R + WIS*H2I
         DATA(I4) = -H1I + WRS*H2I + WIS*H2R
 
C-----------------------------------------------------------------------------
C The recurrence.
C-----------------------------------------------------------------------------
 
         WTEMP = WR
         WR = WR*WPR - WI*WPI + WR
         WI = WI*WPR + WTEMP*WPI + WI

 100  CONTINUE
      IF ( ISIGN.EQ.1 ) THEN
         H1R = DATA(1)
         DATA(1) = H1R + DATA(2)
 
C-----------------------------------------------------------------------------
C Squeeze the first and last data together to get them all within the
C original array.
C-----------------------------------------------------------------------------
 
         DATA(2) = H1R - DATA(2)

      ELSE

         H1R = DATA(1)
         DATA(1) = C1*(H1R+DATA(2))
         DATA(2) = C1*(H1R-DATA(2))
 
C-----------------------------------------------------------------------------
C This is the inverse transform for the case ISIGN = -1.
C-----------------------------------------------------------------------------
 
         CALL PERIOD_FOUR1(DATA, N, -1)

      END IF
 
      RETURN
      END
C____________________________________________________________________________

      SUBROUTINE PERIOD_FOUR1(DATA, NN, ISIGN)
 
C==============================================================================
C Replaces data by its discrete Fourier transform, if ISIGN is input as 1; or
C replaces DATA by NN times its inverse discrete Fourier transform, if ISIGN
C is input as -1. DATA is a complex array of length NN or, equivalently, a real
C array of length 2*NN. NN MUST be an integer power of 2 (this is not checked
C for!).
C
C Adapted from Numerical Recipes by Vikram Singh Dhillon @Sussex 11-Feb-1992.
C==============================================================================
 
C------------------------------------------------------------------------------
C Double precision for the trigonometric recurrences.
C------------------------------------------------------------------------------

      INTEGER          I,J,N,ISIGN,NN,M,MMAX,ISTEP 
      DOUBLE PRECISION WR, WI, WPR, WPI, WTEMP, THETA
      DOUBLE PRECISION  DATA(2*NN),TEMPR,TEMPI

      N = 2*NN
      J = 1
 
C------------------------------------------------------------------------------
C This is the bit-reversal section of the routine.
C------------------------------------------------------------------------------
 
      DO 100 I = 1, N, 2
         IF ( J.GT.I ) THEN
 
C------------------------------------------------------------------------------
C Exchange the two complex numbers.
C------------------------------------------------------------------------------
 
            TEMPR = DATA(J)
            TEMPI = DATA(J+1)
            DATA(J) = DATA(I)
            DATA(J+1) = DATA(I+1)
            DATA(I) = TEMPR
            DATA(I+1) = TEMPI
         END IF
         M = N/2
 50      CONTINUE
         IF ( (M.GE.2) .AND. (J.GT.M) ) THEN
            J = J - M
            M = M/2
            GO TO 50
         END IF
         J = J + M
 100  CONTINUE
 
C------------------------------------------------------------------------------
C Here begins the Danielson-Lanczos section of the routine.
C------------------------------------------------------------------------------
 
      MMAX = 2
 
C------------------------------------------------------------------------------
C Outer loop executed log_2(NN) times.
C------------------------------------------------------------------------------
 
 200  CONTINUE
      IF ( N.GT.MMAX ) THEN
         ISTEP = 2*MMAX
 
C------------------------------------------------------------------------------
C Initialise for the trigonometric recurrence.
C------------------------------------------------------------------------------
 
         THETA = 6.28318530717959D0/(ISIGN*MMAX)
         WPR = -2.D0*DSIN(0.5D0*THETA)**2
         WPI = DSIN(THETA)
         WR = 1.D0
         WI = 0.D0
 
C------------------------------------------------------------------------------
C Here are the two nested inner loops.
C------------------------------------------------------------------------------
 
         DO 250 M = 1, MMAX, 2
            DO 220 I = M, N, ISTEP
               
C------------------------------------------------------------------------------
C This is the Danielson-Lanczos formula.
C------------------------------------------------------------------------------
 
               J = I + MMAX
               TEMPR = (WR)*DATA(J) - (WI)*DATA(J+1)
               TEMPI = (WR)*DATA(J+1) + (WI)*DATA(J)
               DATA(J) = DATA(I) - TEMPR
               DATA(J+1) = DATA(I+1) - TEMPI
               DATA(I) = DATA(I) + TEMPR
               DATA(I+1) = DATA(I+1) + TEMPI
 220        CONTINUE
 
C------------------------------------------------------------------------------
C Trigonometric recurrence.
C------------------------------------------------------------------------------
 
            WTEMP = WR
            WR = WR*WPR - WI*WPI + WR
            WI = WI*WPR + WTEMP*WPI + WI
 250     CONTINUE
         MMAX = ISTEP
 
C------------------------------------------------------------------------------
C Not yet done.
C------------------------------------------------------------------------------
 
         GO TO 200
 
C------------------------------------------------------------------------------
C All done.
C------------------------------------------------------------------------------
 
      END IF
 
      RETURN
      END

*******************************************************************************

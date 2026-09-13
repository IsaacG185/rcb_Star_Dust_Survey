C Copyright the President and Fellows of Harvard College.
C Licensed under the MIT License

c read a power spectrum and find peaks
c select those that are above a given power 
c and write them to the output
c along with period uncertainty, amplitude and date information
c for every detected period, compute the significance and false-alarm probability
c for (1) a blind search, (2) a specified period range
	
	
	implicit none
	character*50 infile,outfile
	integer i,j,lim,nf,ibig,nfr,N
	parameter (lim=1e7)
	double precision  A,B,C,sig,M,Z,norm,thresh,p1,p2,Pmax,Fmax
	double precision F(lim),P(lim),df,period,exptim,Zthresh
	double precision U,fap,over,err1,err2,err3,err4,err5,err5a,err5b
	double precision Mr,sigr,fapr,freq,f1,f2,dw,amp,Z5

c first read the inputs from parameter file

	read(*,'(A50)') infile
	read(*,'(A50)') outfile
	read*,M,norm,N,exptim,thresh,p1,p2

c exptime is the full span of the data, N is number of points in lightcurve.
c M is number of independent frequencies in periodogram
c thresh is the threshold for considering a peak to be real, in units of significance e.g. 95%. 
	
c calculate Zthresh
c the power at which the required significance threshold is reached
	sig=1.0d+00 - (thresh/100.0d+00)
	Zthresh=-1.0*DLOG(sig/M)
	write(*,*) 'significance theshold',thresh,' power=',Zthresh

c read the data	
	open(unit=1,file=infile,status='old')
	open(unit=2,file=outfile,status='unknown')
	write(2,*) 'Period  Amplitude  Perr2(Kovacs)  Perr3(res)  Perr5(delchi)  Signif  FAP'

	Pmax=-1.0d-06
	do i=1,lim
		read(1,*,end=10,err=10) F(i),P(i)
		if(P(i).gt.Pmax)then
			Pmax=P(i)
			Fmax=F(i)
		endif
	enddo
10	nf=i-1
	df=F(11)-F(10)
	f1=1.0d+00/p2
	f2=1.0d+00/p1 
	over=nf/M
	print*,'frequency bins in PDS=',nf,' independent frequencies=',M,' over-sampling factor=',over,' maximum power',Pmax
	print*,'maximum power',Pmax,' at period=',1.0d+00/Fmax

c Loop over the power spectrum. (1) Find the peaks
c				(2) Determine significance
c				(4) Estimate amplitude
c				(3) Estimate period uncertainty
	do i=2,nf-1
C(1) FIND THE PEAKS
                A=P(i-1)
                B=P(i)
                C=P(i+1)
                if(B.gt.A.and.B.gt.C.and.B.gt.Zthresh)then
			Z=P(i)
			freq=F(i)
			period=1.0D+00/F(i)
			ibig=i

C TRANSFORM POWER TO AMPLITUDE (Lomb-Scargle variance normalized method)
C with amplitude defined as the peak-to-trough difference.
			amp=4.0D+00*DSQRT( (Z*norm)/N )

			write(*,*) 'PEAK FOUND',' period=',period,' power=',Z,' Amplitude=',amp


C(2) EVALUATE BLIND SEARCH SIGNIFICANCE   
        		fap=M*DEXP(-1.0d+00*Z)
c			print*,fap
			if(fap.gt.0.2D+00)then
c			    print*,'(1) using full eqn'
			    fap= 1.0D+00 - DEXP(-1.0D+00*Z)
			    fap= 1.0D+00 - fap**M
			endif
			sig= (1.0D+00 - fap)*100.0D+00


c Estimate the Error on the Period
c Version 1: periodogram resolution
c			err1= period - (1.0d+00/(freq+df))
c	 		print*,'err1',err1	

c Version 2: Kovacs 1981 formula: 
c dw = 3*PI*var/2*T*A*sqrt(N)		
c where var is supposed to be variance of noise alone.
			dw=4.71238898d+00 * (dsqrt(norm/N)) * (1.0d+00/(exptim*amp))
			err2=dw/6.283185307d+00
			err2=period - (1.0d+00/(freq+err2))
c       		print*,'err2',err2


c Version3: periodogram resolution modified by the oversampling factor
			err3= period - (1.0d+00/(freq+(df*over)))
			err3= err3/2.0d+00
c	 		print*,'err3',err3

c version 5: Consider the power to be equivalent to a delta-chi2 statistic, 
c hence uncertainty is the range overwhich power changes by 0.5
c this requires rexamining the power spectrum. 
c the peak power is P=Z, held in array element ibig, at frequency=freq 
c i.e. f,Z = F(ibig),P(ibig)
c call the Z-0.5 power Z5	
	
			Z5=Z-0.5d+00
c			print*,'Z',Z
c			print*,'Z5',Z5
			do j=ibig,ibig+1000
			    if(P(j).lt.Z5)then
				err5a=F(j)
c				print*,'err5a',err5a
c                       	print*,'power',P(j)
c                       	print*,'element ibig+',j
				goto 50			
			    endif
			enddo		      
50			do j=ibig,ibig-1000,-1
                	   if(P(j).lt.Z5)then  
                        	err5b=F(j)
c				print*,'err5b',err5b
c				print*,'power',P(j)
c				print*,'element ibig-',j
                        	goto 60
			   endif
        		enddo
60			err5=((1.0d+00/err5b)-(1.0d+00/err5a))/2.0d+00


C WRITE THE OUTPUT
		write(2,200) sngl(period),sngl(amp),sngl(err2),sngl(err3),sngl(err5),sngl(sig),sngl(fap)
200		format(4x,F10.6,4x,F10.6,4x,F10.6,4x,F10.6,4x,F10.6,4x,F10.6,4x,F10.6,4x,F10.6,4x,F10.6)

		endif
	enddo
	end








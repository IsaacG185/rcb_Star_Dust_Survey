C Copyright the President and Fellows of Harvard College.
C Licensed under the MIT License

* modified 6/AUG/2002 SGTL

* Folding program to fold without binning.

	character*50 datafile, fileout
	INTEGER LIMIT,points,icycle 
	PARAMETER( limit=1e6)
	double precision Z,P,cycle, F(1:limit),T(1:limit),E(1:limit),PHASE(1:limit)


c	print*,'PRODUCE EPOCH-FOLDED LIGHTCURVES'
c	print*,'Enter Lightcurve, Folding period & Zeropoint'
c	print*,' ( 0 for first datapoint)'
	
	read (*,'(A50)') datafile
	open (unit=1,file=datafile,status='old')
	read*,P,Z

	do i=1,limit
		read(1,*,end=2,err=2) T(i),F(i),E(i)
	enddo
2	points=i-1

c  calculate mean flux level just out of interest.
	ftotal=0.
	do 4 i=1,points	
		ftotal=ftotal+F(i)
4   	continue		
	fmean=ftotal/points
c	print*,'with Mean Flux level :',fmean			
	 
	
* now fold the data into two arrays according to the period.
	if(Z.eq.0.)Z=T(1)
	do 20 i=1,points
		cycle=(T(i)-Z)/P
		icycle=IDINT(cycle)
		PHASE(i)=cycle-DBLE(icycle)
		if ( PHASE(i) .lt. 0.0D+00 ) PHASE(i) = 1.0d+00 + PHASE(i)
20	continue

C repeat the arrays to a second cycle for plotting
	do 22 i=1,points*2
         	if ( i.gt.points ) then
			PHASE(i)=PHASE(i-points) + 1.0d+00
			F(i)=F(i-points)
			T(i)=T(i-points)
			E(i)=0.0d+00
		endif
22	continue

	do i=1,points*2
		write(*,99) PHASE(i),T(i),F(i),E(i)		
	enddo

99	format (4X,F10.6,4X,F18.6,4X,F10.6,4X,G12.6)

	end

C Copyright the President and Fellows of Harvard College.
C Licensed under the MIT License

C Assign spatial bin based on plate_dist, ra, sec and lookup table
C Both files have a 2-line header

C the lookup table defines the spatial bins both in terms of min/max values for 
C plate_distance and ra/dec. Thus circular and rectangular bins can be accomodated
C in a single script. (pre-supposing that the plates/bins are always oriented N/S etc.)

	character*50 list,lookuptable
	integer lbins,lpts,nbins,npts,check
	parameter (lbins=100, lpts=1e6)
	real XMIN(lbins),XMAX(lbins)
        real RAMIN(lbins),RAMAX(lbins),DECMIN(lbins),DECMAX(lbins)
        real BIN(lbins),VALUE(lbins),x,ra,dec


	read(*,'(A50)')list
	read(*,'(A50)')lookuptable

	open (1,file=list,status='old')
 	open (2,file=lookuptable,status='old')

C skip over the 2-line headers
	read(1,*) 
	read(1,*)
	read(2,*)
        read(2,*)
C read the lookup table into arrays
  	do i=1,lbins
		read(2,*,end=10,err=10) XMIN(i),XMAX(i),RAMIN(i),RAMAX(i),DECMIN(i),DECMAX(i),BIN(i),VALUE(i)
	enddo 
10	nbins=i-1

C read the points and assign a bin and value to each one.
C note that the outer bin is always an edge-trim, 
C which can cut into the radial extent of some of the inner bins. 
C Thus stars in the outer bin do not have meet the minimum plate_distance criterion.
        write(*,*) 'SPATIALBINVAL MAGLIMITVAL'
	write(*,*) '------------- -----------'
	do i=1,lpts
	   read(1,*,end=20,err=20) x,ra,dec
	   check = 0
	   do j=1,nbins
c              write(*,*) XMIN(j),XMAX(j),RAMIN(j),RAMAX(j),DECMIN(j),DECMAX(j) 
              if (x.ge.XMIN(j).and.x.lt.XMAX(j)) then
                 if(ra.ge.RAMIN(j).and.ra.lt.RAMAX(j)) then
			if (dec.ge.DECMIN(j).and.dec.lt.DECMAX(j)) then
         	           write(*,*) BIN(j), VALUE(j)
c                          write(*,*) x,ra,dec 
		           check=check+1
		        endif
		 endif
	      endif
	   enddo
	   if (check.eq.0) then 
	     if (ra.ge.RAMIN(nbins).and.ra.lt.RAMAX(nbins).and.dec.ge.DECMIN(nbins).and.dec.lt.DECMAX(nbins)) then
                write(*,*) BIN(nbins), VALUE(nbins)
		check=check+1
	     endif	
	   endif
	   if (check.eq.0) write(*,*) 99, 99
c           write(*,*) x,ra,dec
	enddo
20	npts=i-1

	end

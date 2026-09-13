% Copyright the President and Fellows of Harvard College.
% Licensed under the MIT License

% Revised defect filters for DASCH
% Sumin Tang, June 12, 2012; June 24, 2012; June 28, 2012
% 
% Compared with the 2008 version, the following improvements are made:
% 1. 2-sigma threshold revised to be 3-sigma (2.5sigma for theta and ellipticity)
% 2. theta vs ellipticity: corrected for discontinuity in theta (0 deg = 180 deg)
% 3. revised defination of real objects (position drad within 2-sigma)
% 4. revised binning in both spatial bins and parameter bins (more dynamic now)
% 5. add cut on FLUX_MAX vs ISO0 (both x and y-axis)
% 
% Input columns from SExtractor:
% NUMBER X_IMAGE Y_IMAGE FLUX_MAX MAG_ISO FWHM_WORLD THETA_J2000 ELLIPTICITY ISO0 plate_dist drad
% Output file: *_defectflag.db
%              2 columns: NUMBER flag
% Overall performance over 16 example plates (stars near the edges are not counted): 
% >90% 2-sigma drad stars are correctly flagged as good (flag=1)
% >90% 4-sigma drad stars (mostly dubious) are correctly flagged as bad (flag=0)

% Jun  4, 2008 Edward J. Los - modify to work under Octave, but only for debugging a C-language conversion
% Jul 14, 2008 Edward J. Los - x and y ranges must be computed from x3 and y3, not x and y
%                            - flag and num must be outside conditionals.

tic
clear
startTime = time();
debug_on_warning(1);

if (nargin() < 1) 
   printf("Usage: octave annular9.m rootname\n");
   %exit;
   rootname="mc39048_01_01r270ww"
else
   rootname = char(argv()(1));
endif
matchdir = getenv("DASCH_MATCH");

for kkk=1:1
    %name="ac05094_00_01ww"; % takes 439 seconds
    %name="mc34766_00_01r90ww";  % takes 574 seconds
    name = rootname;
    numbercheck = -1;

    plate=name;
    datfile = [matchdir '/' rootname '_defect_input.db'];
    numpoints = 0;
    fid = fopen(datfile,'rt');
    if (fid >= 0)     
        while ((txt = fgetl(fid)) != -1)
          numpoints++;
        endwhile
        fclose(fid);
    else
       printf("ERROR: opening %s\n",datfile);
       continue;
    endif
    if (numpoints < 3) 
       printf("ERROR: numpoints is %d for %s\n",numpoints,datfile);
       continue;
    endif
    numpoints = numpoints - 2;
    NUMBER(numpoints) = 0;
    X_IMAGE(numpoints) = 0;
    Y_IMAGE(numpoints) = 0;
    FLUX_MAX(numpoints) = 0;
    MAG_ISO(numpoints) = 0;
    FWHM_WORLD(numpoints) = 0;
    THETA_J2000(numpoints) = 0;
    ELLIPTICITY(numpoints) = 0;
    ISO0(numpoints) = 0;
    plate_dist(numpoints) = 0;
    drad(numpoints) = 0;
    count = 1;
    iii = 0;
    fid = fopen(datfile,'rt');
    txt = fgetl(fid);
    txt = fgetl(fid);
    while ((txt = fgetl(fid)) != -1)
       iii++;
       

       [NUMBER(iii),X_IMAGE(iii),Y_IMAGE(iii),FLUX_MAX(iii),MAG_ISO(iii),FWHM_WORLD(iii),THETA_J2000(iii),ELLIPTICITY(iii),ISO0(iii),plate_dist(iii),drad(iii),count] = sscanf(txt,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f","C");
       if (count != 11)
         printf("Line %d count %d\n",iii,count);
          printf("NUMBER[%d] %s in %s\n",iii,NUMBER(iii),txt);
         iii--;
       endif
    endwhile
    fclose(fid);
    if (iii != numpoints)
        printf("ERROR: numpoints %d is not iii %d\n",numpoints,iii);
        continue
    endif
    endTime = time() - startTime;
    printf("numpoints %d read in from %s at %d seconds\n",numpoints,datfile,endTime);   

    x=X_IMAGE;
    y=Y_IMAGE;
    
    datfile = [matchdir  '/' name '_defect_all.db'];
    numpoints2 = 0;
    fid = fopen(datfile,'rt');
    if (fid >= 0)     
        while ((txt = fgetl(fid)) != -1)
          numpoints2++;
        endwhile
        fclose(fid);
    else
       printf("ERROR: opening %s\n",datfile);
       continue;
    endif
    if (numpoints2 < 3) 
       printf("ERROR: numpoints2 is %d for %s\n",numpoints2,datfile);
       continue;
    endif
    numpoints2 = numpoints2 - 2;
    NUMBER3(numpoints2) = 0;
    x3(numpoints2) = 0;
    y3(numpoints2) = 0;
    FLUX_MAX3(numpoints2) = 0;
    MAG_ISO3(numpoints2) = 0;
    FWHM_WORLD3(numpoints2) = 0;
    THETA_J20003(numpoints2) = 0;
    ELLIPTICITY3(numpoints2) = 0;
    ISO03(numpoints2) = 0;
    plate_dist3(numpoints2) = 0;
    count = 1;
    iii = 0;
    fid = fopen(datfile,'rt');
    txt = fgetl(fid);
    txt = fgetl(fid);
    while ((txt = fgetl(fid)) != -1)
       iii++;
       

       [NUMBER3(iii),x3(iii),y3(iii),FLUX_MAX3(iii),MAG_ISO3(iii),FWHM_WORLD3(iii),THETA_J20003(iii),ELLIPTICITY3(iii),ISO03(iii),plate_dist3(iii),count] = sscanf(txt,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f","C");
       if (count != 10)
         printf("Line %d count %d\n",iii,count);
          printf("NUMBER3[%d] %s in %s\n",iii,NUMBER3(iii),txt);
         iii--;
       endif
    endwhile
    fclose(fid);
    if (iii != numpoints2)
        printf("ERROR: numpoints2 %d is not iii %d\n",numpoints2,iii);
        continue
    endif
    endTime = time() - startTime;
    printf("numpoints2 %d read in from %s at %d seconds\n",numpoints2,datfile,endTime);   

    x=X_IMAGE;
    y=Y_IMAGE;
    
     
    % changed: make the peak of the distribution near the center
    ixa=find(abs(THETA_J2000)<10);
    ixb=find(abs(THETA_J2000)>80);
    THETA=THETA_J2000;   % [-90 90] if centered around 0  
    
    THETA3=THETA_J20003;   % [-90 90] if centered around 0  
    %printf("ixaCount %d ixbCount %d\n",length(ixa),length(ixb));
    if length(ixa)<length(ixb)
        ix=find(THETA_J2000<0);
        THETA(ix)=THETA_J2000(ix)+180;   % [0 180] if centered around +-90
        
        ix3=find(THETA_J20003<0);
        THETA3(ix3)=THETA_J20003(ix3)+180; 
    end
   

   
    % to avoid crash for some problematic objects
    ix=find(FLUX_MAX<=0);
    FLUX_MAX(ix)=1;
    ix=find(ISO0<=0);
    ISO0(ix)=1;
    ix=find(FWHM_WORLD<=0);
    FWHM_WORLD(ix)=1/3600;
    
    iy=find(drad<90);
    meddrad=median(drad(iy));
    stddrad=std(drad(iy));
    
    %printf("drad median is %f from %d points\n",median(drad(find(drad<100))),length(find(drad<100)))

    ixr=find(drad<=meddrad+2*stddrad); % ---- changed: matches within 2-sigma
    ixf=find(drad>=meddrad+4*stddrad); % ---- changed: 4-sigma outliers in position
    %save('los3.tmp','ixr','ixf');

    % devided into 5x5 bins (ensure at least 200 points per bin), then do the local filtering     
    kxr=size(ixr(:),1);
    kall=size(X_IMAGE(:),1);
    nums=0;
    flags=0;
    
    printf("kxr is %d, kall is %d maximum plate distance is %f\n",kxr,kall,max(plate_dist3));

    nbin=min(10,round((kxr/1000).^0.5));  % nbin x nbin bins
    xrange=max(x3)-min(x3);
    yrange=max(y3)-min(y3);
    printf("overall  X: Min %f Max %f range %f  Y: Min %f Max %f range %f\n",min(x3),max(x3),xrange,min(y3),max(y3),yrange);

    x_bin = linspace(min(x3)+0.05*xrange,max(x3)-0.05*xrange,nbin+1);  % larger bins for the edge bins; equal effective area for each bin
    y_bin = linspace(min(y3)+0.05*yrange,max(y3)-0.05*yrange,nbin+1);
    
    x_bin(1)=min(x3);
    x_bin(nbin+1)=max(x3)+1;
    y_bin(1)=min(y3);
    y_bin(nbin+1)=max(y3)+1;

    printf("nbin is %d\n",nbin);
    %save('los4.tmp','x_bin','y_bin');
    for iii=1:nbin % 1:nbin
        for jjj=1:nbin %1:nbin
            endTime = time() - startTime;
            %printf("processing bin %d %d  X: Min %f Max %f range %f  Y: Min %f Max %f range %f  at %d seconds\n",iii,jjj,x_bin(iii),x_bin(iii+1),xrange,y_bin(jjj),y_bin(jjj+1),yrange,endTime);
            ixbin=find(x3>=x_bin(iii) & x3<x_bin(iii+1) & y3>=y_bin(jjj) & y3<y_bin(jjj+1)); % all objects



            % all objects 
            num=NUMBER3(ixbin);
            mag=MAG_ISO3(ixbin);  % instrumental ISO magnitude
            fwhm=log10(3600*FWHM_WORLD3(ixbin)); % log10 FWHM in arcseconds
            theta=THETA3(ixbin);
            e=ELLIPTICITY3(ixbin);
            fluxmax=log10(FLUX_MAX3(ixbin));
            iso0=log10(ISO03(ixbin));
            flag=zeros(size(mag));
            nums=[nums;num'];

            ixcheck = find(numbercheck == NUMBER3(ixbin));
            if length(ixcheck) > 0
                printf("NUMBER %d found in iii %d jjj %d index %d\n",NUMBER3(ixbin(ixcheck)),iii,jjj,ixcheck(1));
            endif
            ixbin2=     find(x>=x_bin(iii)&x<x_bin(iii+1)&y>=y_bin(jjj)&y<y_bin(jjj+1)&drad<90                     & x-min(x3)>0.05*xrange & max(x3)-x>0.05*xrange & y-min(y3)>0.05*yrange& max(y3)-y>0.05*yrange & plate_dist<0.85*max(plate_dist3));
            if (length(ixbin2) > 0)
                % cut off 10% edge for real stars
                ixbinr= find(x>=x_bin(iii)&x<x_bin(iii+1)&y>=y_bin(jjj)&y<y_bin(jjj+1)&drad<=2*median(drad(ixbin2))& x-min(x3)>0.05*xrange & max(x3)-x>0.05*xrange & y-min(y3)>0.05*yrange & max(y3)-y>0.05*yrange & plate_dist<0.85*max(plate_dist3));




                % cut off 20% edge for dubious stars
                ixbinf1=find(x>=x_bin(iii)&x<x_bin(iii+1)&y>=y_bin(jjj)&y<y_bin(jjj+1)&drad>4*median(drad(ixbin2)) & x-min(x3)>0.1*xrange & max(x3)-x>0.1*xrange & y-min(y3)>0.1*yrange & max(y3)-y>0.1*yrange & plate_dist<0.7*max(plate_dist3)); % 4-sigma objects
    
 
                
                if length(ixbinr)>50   % changed: put the requirement inside the loop for bins
                    % real objects
                    magr=MAG_ISO(ixbinr);  % instrumental ISO magnitude
                    fwhmr=log10(3600*FWHM_WORLD(ixbinr)); % log10 FWHM in arcseconds
                    thetar=THETA(ixbinr);
                    er=ELLIPTICITY(ixbinr);
                    fluxmaxr=log10(FLUX_MAX(ixbinr));
                    iso0r=log10(ISO0(ixbinr));
                    
                    % fake objects as 4-sigma object
                    magf1=MAG_ISO(ixbinf1);  % instrumental ISO magnitude
                    fwhmf1=log10(3600*FWHM_WORLD(ixbinf1)); % log10 FWHM in arcseconds
                    thetaf1=THETA(ixbinf1);
                    ef1=ELLIPTICITY(ixbinf1);
                    fluxmaxf1=log10(FLUX_MAX(ixbinf1));
                    iso0f1=log10(ISO0(ixbinf1));                
                    
                    
                    % filter 1-2: mag vs fwhm, ellipticity (theta does not depend on mag), 15 bins in mag
                    %             two iterations of 4-sigma clipping, 3-sigma rejection
                    % filter 3: 3-sigma clipping on theta, twice; 4-sigma on ellipticity
                    % filter 4: fluxmax vs iso0, 10 bins in mag
                    
                    
                    % change the binning to prctile; each bin contains similar number of stars
                    nmag=min(30,round(length(magr)/30));
                    maggrid=sort(prctile(magr', sort([0:100*15/length(magr):80/nmag, 100/nmag:100/nmag:100])));
                    maggrid(2:length(maggrid)-1) += 0.0001;                   

                    %percentvector = sort([0:100*15/length(magr):80/nmag, 100/nmag:100/nmag:100]);
                    %for mmm=1:length(maggrid)
                    %    printf("index %d, percent %f, vector %f\n",mmm-1,percentvector(mmm),maggrid(mmm));
                    %end

                    
                    ixfa=0;
                    ixfb=0;
                    clear fwhmc rmsfwhmc ec rmsec
                    for i=1:length(maggrid)-1
                        ix=find(magr>=maggrid(i)&magr<maggrid(i+1));

                        %X_IMAGEX=X_IMAGE(ixbinr);
                        %Y_IMAGEX=Y_IMAGE(ixbinr);
                        %NUMBERX=NUMBER(ixbinr);
                        %X_IMAGE4=X_IMAGEX(ix);
                        %Y_IMAGE4=Y_IMAGEX(ix);
                        %NUMBER4=NUMBERX(ix);
                        %magr4=magr(ix);
                        %debugfile = sprintf('%s/losixf2%d%d_%dm.tmp',matchdir,iii,jjj,i);
                        %fid = fopen(debugfile,'w');
                        %fprintf(fid,'NUMBER\tX_IMAGE\tY_IMAGE\tMAG_ISO\tcount\n');
                        %fprintf(fid,'------\t-------\t-------\t-------\t-----\n');
                        %for lll=1:length(ix)
                        %   fprintf(fid,'%d\t%f\t%f\t%f\t%d\n',NUMBER4(lll),X_IMAGE4(lll),Y_IMAGE4(lll),magr4(lll),(10000*iii)+(100*jjj)+i);
                        %end
                        %fclose(fid);     


                        maggridOK(i) = length(ix);
                        %if ((iii == 3) && (jjj == 3)) 
                        %   printf("Bin %d %d, ixbincount %5d, ixbin2count %5d, ixbinrcount %5d ixbinf1count %5d\n",
                        %        iii,jjj,length(ixbin),length(ixbin2),length(ixbinr),length(ixbinf1));
                        %endif
    
                        %printf("FWHM index length %d for iii %d jjj %d i %d\n",length(ix),iii,jjj,i);
                        if (maggridOK(i) > 0)
        
                            fwhmc0=mean(fwhmr(ix));
                            rmsfwhmc0=min(std(fwhmr(ix)),0.2);
                            ix0=find(abs(fwhmr(ix)-fwhmc0)<=3*rmsfwhmc0); % 3-sigma clipping
                            maggridOK(i) = length(ix0);
                            if (maggridOK(i) > 0)
                              rmsfwhmc0=std(fwhmr(ix(ix0)));
                              ix0=find(abs(fwhmr(ix)-fwhmc0)<=3*rmsfwhmc0); % 3-sigma clipping
                              maggridOK(i) = length(ix0);
                            endif
                            if (maggridOK(i) > 0)
                              rmsfwhmc0=std(fwhmr(ix(ix0)));
                              ix0=find(abs(fwhmr(ix)-fwhmc0)<=3*rmsfwhmc0); % 3-sigma clipping
                              maggridOK(i) = length(ix0);
                            endif
                            if (maggridOK(i) > 0)
                              fwhmc(i)=mean(fwhmr(ix(ix0)));
                              rmsfwhmc(i)=std(fwhmr(ix(ix0)));

                              %printf("Bin %d %d %d fwhmc %f rmsfwhmc %f ix0 %d\n",iii,jjj,i,fwhmc(i),rmsfwhmc(i),length(ix0));


                              ixa=find(mag>=maggrid(i)&mag<maggrid(i+1)&abs(fwhm-fwhmc(i))<3*rmsfwhmc(i));
                              %if ((iii == 3) && (jjj == 3)) 
                              %   printf("Bin %d %d, i %d fwhm %f rmsfwhmc %f, length(fwhm) %d length ixa %d\n",iii,jjj,i,fwhmc,rmsfwhmc,length(fwhm),length(ixa));
                              %endif
                              ixfa=[ixfa;ixa'];
                            endif
                            ec0=mean(er(ix));                    
                            rmsec0=std(er(ix));                    
                            ix0=find(abs(er(ix)-ec0)<=3*rmsec0); % 3-sigma clipping
                            rmsec0=std(er(ix(ix0)));
                            ix0=find(abs(er(ix)-ec0)<=3*rmsec0); % 3-sigma clipping
                            rmsec0=std(er(ix(ix0)));
                            ix0=find(abs(er(ix)-ec0)<=3*rmsec0); % 3-sigma clipping
                            ec(i)=mean(er(ix(ix0)));
                            rmsec(i)=std(er(ix(ix0)));

                            %printf("Bin %d %d %d ec %f rmsec %f ix0 %d\n",iii,jjj,i,ec(i),rmsec(i),length(ix0));

                            ixb=find(mag>=maggrid(i)&mag<maggrid(i+1)&abs(e-ec(i))<3*rmsec(i));
                            %if ((iii == 3) && (jjj == 3)) 
                            %  printf("Bin %d %d, i %d e %f rmsec %f, length(e) %d length ixb %d\n",iii,jjj,i,ec,rmsec,length(e),length(ixb));
                            %endif
                            ixfb=[ixfb;ixb'];      
                            %printf("Bin %d %d, i -1 fwhm %f rmsfwhmc %f, length(fwhm) %d length ixa %d\n",iii,jjj,fwhmc,rmsfwhmc,length(fwhm),length(ixa));
                            %if ((iii == 3) && (jjj == 3)) 
                            %   printf("Bin %d %d, i -1 e %f rmsec %f, length(e) %d length ixb %d\n",iii,jjj,ec,rmsec,length(e),length(ixb));
                            %endif
                        else
                            printf("FWHM index is zero for iii %d jjj %d i %d\n",iii,jjj,i);
                        endif
                    end
                    ixf1=ixfa(find(ixfa));
                    ixf2=ixfb(find(ixfb));
                    figure(1,'visible','off');
                    subplot(2,2,1)
                    plot(magr,fwhmr,'.','markersize',4)
                    hold on
                    if (length(magf1) > 0) 
                        plot(magf1,fwhmf1,'ro', 'markersize',2)
                    endif            
                    for i=1:length(maggrid)-1
                        if (maggridOK(i) > 0)
                          myx=[maggrid(i) maggrid(i+1)];
                          myy=[1 1]*3*rmsfwhmc(i)+fwhmc(i);
                          plot(myx,myy,'k-')
                          myy=-[1 1]*3*rmsfwhmc(i)+fwhmc(i);
                          plot(myx,myy,'k-')
                       endif
                    end
                    if length(ixcheck) > 0
                        printf("NUMBER %d MAG_ISO3 %f log10(FWHM_WORLD3) %f\n",NUMBER3(ixbin(ixcheck)),MAG_ISO3(ixbin(ixcheck)),log10(3600*FWHM_WORLD3(ixbin(ixcheck))));
                        plot(MAG_ISO3(ixbin(ixcheck)),log10(3600*FWHM_WORLD3(ixbin(ixcheck))),'go','markersize',8);
                    endif
                    xlabel('MAG\_ISO','Fontsize',16);
                    ylabel('log(FWHM) (arcsec)','Fontsize',16);
                    title([' ' plate ', [' num2str(iii) ',' num2str(jjj) ']'],'Fontsize',16);
                    
                    subplot(2,2,2)
                    plot(magr,er,'.','markersize',4)
                    hold on
                    plot(magf1,ef1,'ro', 'markersize',2)                
                    for i=1:length(maggrid)-1
                        myx=[maggrid(i) maggrid(i+1)];
                        myy=[1 1]*3*rmsec(i)+ec(i);
                        plot(myx,myy,'k-')
                        myy=-[1 1]*3*rmsec(i)+ec(i);
                        plot(myx,myy,'k-')
                    end
                    if length(ixcheck) > 0
                        printf("NUMBER %d MAG_ISO3 %f ELLIPTICITY3 %f\n",NUMBER3(ixbin(ixcheck)),MAG_ISO3(ixbin(ixcheck)),ELLIPTICITY3(ixbin(ixcheck)));
                        plot(MAG_ISO3(ixbin(ixcheck)),ELLIPTICITY3(ixbin(ixcheck)),'go','markersize',8);
                    endif
                    xlabel('MAG\_ISO','Fontsize',16);
                    ylabel('Ellipticity','Fontsize',16);
    
                    
                    % filter 3: 3-sigma clipping on theta, 2.5-sigma on ellipticity
    
  


                    if (length(thetar) == 0) 
                       thetac0 = 0;
                    else
                       thetac0=median(thetar);
                    endif
    
                    % changed: define distance to the median(thetar)
                    thetadist=abs(theta-thetac0);
                    ix=find(thetadist>90);
                    thetadist(ix)=180-thetadist(ix);
                    
                    thetadistr=abs(thetar-thetac0);
                    ix=find(thetadistr>90);
                    thetadistr(ix)=180-thetadistr(ix);
                    rmsthetac0=(sum(thetadistr.^2)/length(thetadistr))^0.5; % revised std, to reflect the continuous distribution of theta, i.e. 0==180
                    
                    ix0=find(thetadistr<=2*rmsthetac0); % 2-sigma clipping
                    rmsthetac0=(sum(thetadistr(ix0).^2)/length(ix0))^0.5;
                    ix0=find(thetadistr<=2*rmsthetac0); % 2-sigma clipping
                    rmsthetac0=(sum(thetadistr(ix0).^2)/length(ix0))^0.5;
                    ix0=find(thetadistr<=2*rmsthetac0); % 2-sigma clipping
                    rmsthetac=(sum(thetadistr(ix0).^2)/length(ix0))^0.5;
                    %printf("Bin %d %d thetac0 %f rmsthetac0 %f, length(ix0) %d\n",iii,jjj,thetac0,rmsthetac,length(ix0));
                    %if ((iii == 3) && (jjj == 3))
                    %    zzz = thetar(ix0);
                    %    save('los5.tmp','thetar','ix0','zzz','thetac0','rmsthetac0');
                    %    printf("Bin %d %d thetac0 %f rmsthetac0 %f, length(ix0) %d\n",iii,jjj,thetac0,rmsthetac0,length(ix0));
                    %endif
                    %if ((iii == 3) && (jjj == 3))
                    %   printf("Bin %d %d thetac0 %f rmsthetac0 %f, length(ix0) %d\n",iii,jjj,thetac0,rmsthetac0,length(ix0));
                    %endif
                    
                    if (length(er) == 0) 
                        ec0 = 0;
                    else 
                        ec0=median(er);
                    endif
                    rmsec0=std(er);
                    ix0=find(abs(er-ec0)<=2*rmsec0); % 2-sigma clipping
                    %if ((iii == 3) && (jjj == 3))
                    %    zzz = er(ix0);
                    %    save('los6.tmp','er','ix0','zzz','ec0','rmsec0');
                    %   printf("Bin %d %d ec0 %f rmsec0 %f, length(ix0) %d\n",iii,jjj,ec0,rmsec0,length(ix0));
                    %endif
                    rmsec0=std(er(ix0));
                    ix0=find(abs(er-ec0)<=2*rmsec0); % 2-sigma clipping
                    %if ((iii == 3) && (jjj == 3))
                    %   printf("Bin %d %d ec0 %f rmsec0 %f, length(ix0) %d\n",iii,jjj,ec0,rmsec0,length(ix0));
                    %endif
                    if (length(ix0) == 0)
                      ec = 0;
                      rmsec = 0;
                    else
                      ec=mean(er(ix0));
                      rmsec=std(er(ix0));
                    endif
                    %printf("Bin %d %d ec %f rmsec %f, length(ix0) %d\n",iii,jjj,ec,rmsec,length(ix0));
                    ixf3=find(thetadist<2.5*rmsthetac&abs(e-ec)<2.5*rmsec);   % 2.5 sigma threshold 
                    %if ((iii == 3) && (jjj == 3)) 
                    %     printf("Bin %d %d, thetac0 %f rmsthetac %f, ec %f, rmsec %f, theta length %d e length %d ixf3 length %d\n",iii,jjj,thetac0,rmsthetac,ec,rmsec,length(theta),length(e),length(ixf3));
                    %endif
    
                    
                    subplot(2,2,3)
                    plot(thetar,er,'.','markersize',4)
                    hold on
                    plot(thetaf1,ef1,'ro', 'markersize',2)               
                    plot([min(thetar) max(thetar)],[1 1]*rmsec*2.5+ec,'k--')
                    plot([min(thetar) max(thetar)],-[1 1]*rmsec*2.5+ec,'k--')
                    plot([1 1]*thetac0+2.5*rmsthetac,[0 1],'k--')
                    plot([1 1]*thetac0-2.5*rmsthetac,[0 1],'k--')
                    if length(ixcheck) > 0
                        printf("NUMBER %d THETA3 %f ELLIPTICITY3 %f\n",NUMBER3(ixbin(ixcheck)),THETA3(ixbin(ixcheck)),ELLIPTICITY3(ixbin(ixcheck)));
                        plot(THETA3(ixbin(ixcheck)),ELLIPTICITY3(ixbin(ixcheck)),'go','markersize',8);
                    endif
                    xlabel('\theta','Fontsize',16);
                    ylabel('Ellipticity','Fontsize',16);
                    ylim([0 1])
                    xlim([min(thetar) max(thetar)])
                    
                    
                    % filter 4: fluxmax vs iso0, 15 bins in iso0
                    % binning changed to be that each bin contains similar number of stars
                    niso0=min(30,round(length(iso0r)/30));
                    iso0grid=sort(prctile(iso0r', sort(100-[0:100*15/length(magr):80/niso0, 100/niso0:100/niso0:100])));
                    iso0grid(2:length(iso0grid)-1) += 0.0001;                   
                    
                    %percentvector = sort(100-[0:100*15/length(magr):80/niso0, 100/niso0:100/niso0:100]);
                    %for mmm=1:length(iso0grid)
                    %    printf("index %d, percent %f, vector %f\n",mmm-1,percentvector(mmm),iso0grid(mmm));
                    %end
                    


                    ixfd=0;
                    clear fluxmaxc rmsfluxmaxc
                    for i=1:length(iso0grid)-1

                        %if (iso0grid(i+1) < iso0grid(i)) 
                        %   printf("Inversion at i %d iso0grid %f %f\n",i,iso0grid(i),iso0grid(i+1));
                        %   isog0grid(i+1) = iso0grid(i);
                        %endif
                        ix=find(iso0r>=iso0grid(i)&iso0r<iso0grid(i+1));


                        X_IMAGEX=X_IMAGE(ixbinr);
                        Y_IMAGEX=Y_IMAGE(ixbinr);
                        NUMBERX=NUMBER(ixbinr);
                        X_IMAGE4=X_IMAGEX(ix);
                        Y_IMAGE4=Y_IMAGEX(ix);
                        NUMBER4=NUMBERX(ix);
                        iso0r4=iso0r(ix);
                        debugfile = sprintf('%s/losixf4_%d%d_%dm.tmp',matchdir,iii,jjj,i);
                        fid = fopen(debugfile,'w');
                        fprintf(fid,'NUMBER\tX_IMAGE\tY_IMAGE\tISO0\tcount\n');
                        fprintf(fid,'------\t-------\t-------\t----\t-----\n');
                        for lll=1:length(ix)
                           fprintf(fid,'%d\t%f\t%f\t%f\t%d\n',NUMBER4(lll),X_IMAGE4(lll),Y_IMAGE4(lll),iso0r4(lll),(10000*iii)+(100*jjj)+i);
                        end
                        fclose(fid);     



                        iso0gridOK(i) = length(ix);
                        %printf("FLUX_MAX index length %d for iii %d jjj %d i %d\n",length(ix),iii,jjj,i);
                        if (iso0gridOK(i) > 0) 
                            fluxmaxc0=mean(fluxmaxr(ix));
                            rmsfluxmaxc0=max(std(fluxmaxr(ix)),0.3);
                            ix0=find(abs(fluxmaxr(ix)-fluxmaxc0)<=2*rmsfluxmaxc0); % 2-sigma clipping
                            iso0gridOK(i) = length(ix0);
                            if length(ix0) > 0
                               rmsfluxmaxc0=std(fluxmaxr(ix(ix0)));
                               ix0=find(abs(fluxmaxr(ix)-fluxmaxc0)<=2*rmsfluxmaxc0); % 2-sigma clipping
                               iso0gridOK(i) = length(ix0);
                            endif
                            if length(ix0) > 0
                               rmsfluxmaxc0=std(fluxmaxr(ix(ix0)));
                               ix0=find(abs(fluxmaxr(ix)-fluxmaxc0)<=3*rmsfluxmaxc0); % 3-sigma clipping
                               iso0gridOK(i) = length(ix0);
                            endif
                            if (length(ix0) > 0) 
                                fluxmaxc(i)=mean(fluxmaxr(ix(ix0)));
                                rmsfluxmaxc(i)=std(fluxmaxr(ix(ix0)));
                                if i>1
                                    rmsfluxmaxc(i)=min(std(fluxmaxr(ix(ix0))),abs(fluxmaxc(i)-fluxmaxc(i-1))+rmsfluxmaxc(i-1));
                               end

                               %printf("Bin %d %d %d fluxmaxc %f rmsfluxmaxc %f ix0 %d\n",iii,jjj,i,fluxmaxc(i),rmsfluxmaxc(i),length(ix0));

                               ixd=find(iso0>=iso0grid(i)&iso0<iso0grid(i+1)& abs(fluxmax-fluxmaxc(i))<3*rmsfluxmaxc(i));
                               %if ((iii == 3) && (jjj == 3)) 
                               %     printf("Bin %d %d, i %d iso0c %f rmsiso0c %f, length(iso0) %d length ixd %d\n",iii,jjj,i,iso0c,rmsiso0c,length(iso0),length(ixd));
                               %endif
                               ixfd=[ixfd;ixd'];
                               %if (length(find(ixd==16))) > 0
                               %   i
                               %   ixd
                               %endif
                          endif
                        else 
                            %printf("FLUX_MAX index is zero for iii %d jjj %d i %d\n",iii,jjj,i);
                        endif
                    end
                    ixf4=ixfd(find(ixfd));
                     
                    
                    % bin in FLUX_MAX
                    niso0=min(30,round(length(iso0r)/30));
                    fluxmaxgrid=sort(prctile(fluxmaxr', sort(100-[0:100*15/length(magr):80/niso0, 100/niso0:100/niso0:100])));                                
                    fluxmaxgrid(2:length(fluxmaxgrid)-1) += 0.0001;                   


                    ixfe=0;
                    clear iso0c rmsiso0c
                    for i=1:length(fluxmaxgrid)-1
                        ix=find(fluxmaxr>=fluxmaxgrid(i)&fluxmaxr<fluxmaxgrid(i+1));





                        iso0cOK(i) = length(ix);
                        if (iso0cOK(i) > 0) 
                           iso0c0=mean(iso0r(ix));
                           rmsiso0c0=max(std(iso0r(ix)),0.3);
                           ix0=find(abs(iso0r(ix)-iso0c0)<=2*rmsiso0c0); % 2-sigma clipping
                           iso0cOK(i) = length(ix0);
                        endif
                        if (iso0cOK(i) > 0) 
                            rmsiso0c0=std(iso0r(ix(ix0)));
                            ix0=find(abs(iso0r(ix)-iso0c0)<=2*rmsiso0c0); % 2-sigma clipping
                            iso0cOK(i) = length(ix0);
                        endif
                        if (iso0cOK(i) > 0) 
                          rmsiso0c0=std(iso0r(ix(ix0)));
                          ix0=find(abs(iso0r(ix)-iso0c0)<=3*rmsiso0c0); % 3-sigma clipping
                          iso0cOK(i) = length(ix0);
                        endif
                        if (iso0cOK(i) > 0) 
                           iso0c(i)=mean(iso0r(ix(ix0)));
                           rmsiso0c(i)=min(std(iso0r(ix(ix0))),1);

                           %printf("Bin %d %d %d iso0c %f rmsiso0c %f ix0 %d\n",iii,jjj,i,iso0c(i),rmsiso0c(i),length(ix0));

                           ixe=find(fluxmax>=fluxmaxgrid(i)&fluxmax<fluxmaxgrid(i+1)& abs(iso0-iso0c(i))<3*rmsiso0c(i));
                           ixfe=[ixfe;ixe'];
                        endif
                    end
                    ixf5=ixfe(find(ixfe));
                    
                    
                    
                    subplot(2,2,4)
                    plot(iso0r,fluxmaxr,'.','markersize',4)
                    hold on
                    if length(iso0f1) > 0 
                        plot(iso0f1,fluxmaxf1,'ro', 'markersize',2)   
                    endif  
                    for i=1:length(iso0grid)-1
                        if iso0gridOK(i) > 0 
    
                           myx=[iso0grid(i) iso0grid(i+1)];
                           myy=[1 1]*3*rmsfluxmaxc(i)+fluxmaxc(i);
                           plot(myx,myy,'k-')
                           myy=-[1 1]*3*rmsfluxmaxc(i)+fluxmaxc(i);
                           plot(myx,myy,'k-')
                        endif
                    end
                    for i=1:length(fluxmaxgrid)-1
                        if (iso0cOK(i) > 0)
                            myy=[fluxmaxgrid(i) fluxmaxgrid(i+1)];
                            myx=[1 1]*3*rmsiso0c(i)+iso0c(i);
                            plot(myx,myy,'k-')
                            myx=-[1 1]*3*rmsiso0c(i)+iso0c(i);
                            plot(myx,myy,'k-')
                        endif
                    end
                        plot(iso0f1,fluxmaxf1,'ro', 'markersize',2)   
                    if length(ixcheck) > 0
                        printf("NUMBER %d log10(ISO03) %f log10(FLUX_MAX3) %f\n",NUMBER3(ixbin(ixcheck)),log10(ISO03(ixbin(ixcheck))),log10(FLUX_MAX3(ixbin(ixcheck))));
                        plot(log10(ISO03(ixbin(ixcheck))),log10(FLUX_MAX3(ixbin(ixcheck))),'go','markersize',8);
                    endif
                    ylabel('FLUX\_MAX','Fontsize',16);
                    xlabel('log(ISO0) (pixel^2)','Fontsize',16);
     
                    print([plate '_' num2str(iii) '-' num2str(jjj) '_plot.eps'],'-color');
                    %%%saveas(gcf,[plate '_' num2str(iii) '-' num2str(jjj) '_plot.eps'],'epsc2')
                    close
                    % only keep the objects passed all the tests:
                    printf("Bin %d %d, ixf1 %5d ixf2 %5d ixf3 %5d ixf4 %5d ixf5 %5d\n",iii,jjj,length(ixf1),length(ixf2),length(ixf3),length(ixf4),length(ixf5));
                    
                    %save(['los' num2str(iii) '_' num2str(jjj) '.tmp'],'ixf1','ixf2','ixf3','ixf4');
                    ixfall=[ixf1;ixf2;ixf3';ixf4;ixf5];
                    [m,n]=hist(ixfall,1:size(mag(:),1));
                    ix=find(m==5);
                    
                    flag(ix)=1;
                    
                end



            end
            flags=[flags;flag'];
        end
    end        
    kall=length(nums)-1;
    output=[nums(2:kall+1) flags(2:kall+1)];
    fid = fopen([matchdir  '/' name '_defectflag.db'],'w');
    fprintf(fid,'%7.0f %1.0f\n',output');
    fclose(fid);
    endTime = time() - startTime;
    printf("Creating summary at %d seconds\n",endTime);        
    % get 2 ratios: real pass ratio, fake filtered out ratio, only for MAG_ISO<max(MAG_ISO)-1:
    % match output and NUMBER drad
    numsall=nums(2:kall+1);
    flagsall=flags(2:kall+1);
    % [C,ia,ib] = intersect(A,B) also returns index vectors ia and ib, such that C = A(ia) and C = B(ib)
    [C,Ia,Ib] = intersect(numsall,NUMBER); % ia: ranked by numsall; ib: ranked by NUMBER
    Ya=numsall(Ia);
    Yb=NUMBER(Ib);
    drads=drad(Ib);
    % ERROR check for plate_dist3
    platedist=plate_dist(Ib);
    magisos=MAG_ISO(Ib);
    x4=x(Ib);
    y4=y(Ib);
    flagss=flagsall(Ia);
    output2=[Yb drads flagss'];
%   save('los2.tmp','output2');


    
    % define fake as 4-sigma
    ixsr=find(drads<=2*median(drads(find(drads<90))) & magisos<max(magisos)-1& x4-min(x4)>0.05*xrange & max(x4)-x4>0.05*xrange & y4-min(y4)>0.05*yrange & max(y4)-y4>0.05*yrange & platedist<0.85*max(plate_dist3));
    ixsf=find(drads>4*median(drads(find(drads<90)))& magisos<max(magisos)-1& x4-min(x4)>0.1*xrange & max(x4)-x4>0.1*xrange & y4-min(y4)>0.1*yrange & max(y4)-y4>0.1*yrange & platedist<0.7*max(plate_dist3));
    flagsr=flagss(ixsr);
    flagsf=flagss(ixsf);
%    save('los1.tmp','flagsr','flagsf','ixsr','ixsf');
    ratio1=size(find(flagsr),1)/size(ixsr(:),1);     % 2sigma stars kept
    if (size(ixsf(:),1) > 0)
        ratio2=size(find(flagsf==0),1)/size(ixsf(:),1);  % 4sigma stars flagged as defects
    else
        ratio2 = 9999;
    endif

    ratio3=size(ixsf(:),1)/size(ixsr(:),1);
    numreal=size(ixsr(:),1);
    numfake=size(ixsf(:),1);
    
    
    % plot 4-sigma object -----
    ix4sig=find( drads>4*median(drads(find(drads<90))) &  magisos<max(magisos)-1 ); 
    plot(x4(ix4sig),y4(ix4sig),'.','markersize',4)
    hold on
    plot([1 1]*min(x3)+0.1*range(x),[min(y3) max(y3)],'g')
    plot([1 1]*max(x3)-0.1*range(x),[min(y3) max(y3)],'g')
    plot([min(x3) max(x3)],[1 1]*min(y3)+0.1*range(y),'g')
    plot([min(x3) max(x3)],[1 1]*max(y3)-0.1*range(y),'g')
    s=0:0.01:2*pi;
    s1=((max(x3)^2+max(y3)^2)^0.5/2)*0.7*cos(s)+(max(x3)+min(x3))/2;
    s2=((max(x3)^2+max(y3)^2)^0.5/2)*0.7*sin(s)+(max(y3)+min(y3))/2;
    plot(s1,s2,'g')
    xlabel('X\_IMAGE','fontsize',18)
    ylabel('Y\_IMAGE','fontsize',18)
    xlim([0 max(x3)])
    ylim([0 max(y3)])
    title([plate ', ' num2str(length(ix4sig)) '; ' num2str(numfake) ' 4\sigma drad objects (' num2str(length(x)) ';' num2str(numreal) ' matched)'],'fontsize',12)   
    set(gca,'linewidth',1,'fontsize',16)
    daspect([1 1 1]) % equal scaling along each axis
    print([plate '-4sigma.eps'],'-color');
    %%%saveas(gcf,[plate '-4sigma.eps'],'epsc2')
    daspect('auto') % sets the data aspect ratio mode to auto.
    close
    


    endTime = time() - startTime;
    printf("Done at %d seconds\n",endTime);        

    ratios=[ratio1' ratio2' ratio3' numreal' numfake'];
    fid = fopen([matchdir  '/' plate '_defectflag_ratios.db'],'w');
    fprintf(fid,'2sigma-in-ratio  4sigma-out-ratio n_4sigma/n_2sigma n_2sigma n_4sigma \n');
    fprintf(fid,'%2.4f %2.4f %2.4f %7.0f %7.0f\n',ratios');
    fclose(fid);     
    
end

 
toc


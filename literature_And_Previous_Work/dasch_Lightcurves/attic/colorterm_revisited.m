% Copyright the President and Fellows of Harvard College.
% Licensed under the MIT License

% Do colorterm fitting before rlowess
%    Final verson
% Sumin Tang, July 26, 2008 
% ----------------------------------------
% 1. do a rough fitting, do a correlation between dmag and color, to decide the starting point
% 2. use metropolis algorithm to sampling color-term
% 3. stop when n>100 or ?
% 4. fit the curve, determine the color-term
% This is just after dividing the plate into 9 bins, and before the annular9.m
% The input file should be /dasch/raid008/Pipeline/catalog9bin/<plate>_a1-9.db
% But I don't know how to handle the format of input files <plate>_a1-9.db
% So instead, here I'm using the output files /dasch/Pipeline/bin9/<plate>_a1-9.out
% Only the following 5 columns are indeed used:
% iso stdmag color dra ddec
%
% Aug 25, 2008 Edward J Los - Port to Octave
% Jan  7, 2008 Edward J Los - use linear fit if only two points are available for 
%                             curve-fitting


tic
clear
debugSuffix2 = "";
scripts = getenv("DASCH_SCRIPTS");
printf("scripts is %s\n",scripts);
if (length(scripts) < 1) 
  exit;
endif
catalogbin = getenv("DASCH_CATALOGBIN");
printf("catalogbin is %s\n",catalogbin);
if (length(catalogbin) < 1) 
  exit;
endif
    oldformat = 0;



qualifier = "";

chdir(scripts);
pwd();    % this prevents a synchronization problem.


  listfile = argv()(1);
  if (nargin() == 2) 
      qualifier =  ['_' char(argv()(2))];
  endif
  printf("List file is %s qualifier is %s\n",listfile,qualifier);
  
  
  
  fid = fopen(listfile,'rt');
  iii = 1;
  B = cell(1);
  while ((txt = fgetl(fid)) != -1)
     B(iii) = txt;
     iii++;
  endwhile
  ksize = iii - 1;
  fclose(fid);
  


for kkk=1:ksize
    name=char(B(kkk));
  for bin=1:9 % 1:9
    chdir(catalogbin)
    pwd();
    printf("oldformat is %d for %s bin %d\n",oldformat,name,bin);
    if (oldformat == 1) 
       fid = fopen(['/dasch/junk/colorterm/input/' name qualifier '_a' num2str(bin) debugSuffix2 '.out'],'r');
    else % oldformat
       fid = fopen([name qualifier '_a' num2str(bin) debugSuffix2 '.db'],'r');
    endif % oldformat
    if (fid == -1)
       printf("Failed to open %s\n",[name qualifier '_a' num2str(bin) debugSuffix2 '.db']);
       continue;
    else
       if (oldformat == 1) 
          numpoints = 0;
          while ((txt = fgetl(fid)) != -1)
             numpoints++;
          endwhile
          printf("numpoints is %d for %s\n",numpoints,['/dasch/junk/colorterm/input/' name qualifier '_a' num2str(bin) debugSuffix2 '.out']);
       endif % oldformat
       fclose(fid);
    endif
    if (oldformat == 1) 
      b1 = cell(numpoints,1);
      a1(numpoints) = 0;
      a2(numpoints) = 0;
      a3(numpoints) = 0;
      a4(numpoints) = 0;
      a5(numpoints) = 0;
      a6(numpoints) = 0;
      a7(numpoints) = 0;
      a8(numpoints) = 0;
      a9(numpoints) = 0;
      a10(numpoints) = 0;
      a11(numpoints) = 0;
      a12(numpoints) = 0;
      count = 1;
      iii = 0;
       fid = fopen(['/dasch/junk/colorterm/input/' name qualifier '_a' num2str(bin) debugSuffix2 '.out'],'r');
 
      while ((txt = fgetl(fid)) != -1)
        iii++;
        % keep a1 around to avoid renumbering all of the indices 
        a1(iii) = 0;
        [b1(iii,1),a2(iii),a3(iii),a4(iii),a5(iii),a6(iii),a7(iii),a8(iii),a9(iii),a10(iii),a11(iii),a12(iii),count] = sscanf(txt,"%s %f %f %f %f %f %f %f %f %f %f %f","C");
        if (count != 12) 
           printf("Line %d count %d\n",iii,count);
           printf("b1[%d] %s\n",iii,b1(iii,1));
           iii--;
        endif
     endwhile
     fclose(fid);
     printf("Finished read of %s of size %8d %8d\n",name,numpoints,iii);
     stdmag = a4';
     color = a5';
     iso = a6';
     ddec = a11';
     dra = a12';

     clear b1 a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 a12

    else % oldformat
      load([name qualifier '_a' num2str(bin) debugSuffix2 '.db'],['output' num2str(bin) 'REF'],['output' num2str(bin)]);
      command = ['input = output' num2str(bin) ';'];
      eval(command);
      clear(['output' num2str(bin)]);
      command = ['inputREF = output' num2str(bin) 'REF;'];
      eval(command);
      clear(['outputREF' num2str(bin)]);
   
      stdmag = input(:,4);
      color = input(:,5);
      iso = input(:,6);
      ddec = input(:,11);
      dra = input(:,12);
    endif %oldformat
    
    % For MS stars, B-R is about [-0.5,2];
    % only use stars with stdmag between 4 and 19 (to avoid crazy magnitudes)
    % and use drad filter - use absolute value!
    ix=find(stdmag>4 & stdmag<19 & color>-0.5 & color<2 & abs(ddec)<3*median(abs(ddec)) & abs(dra)<3*median(abs(dra)));  
    
    stdmag0=stdmag(ix);
    iso0=iso(ix);
    color0=color(ix);
    
    [stdmag,ir]=sort(stdmag0); % sort by stdmag
    iso=iso0(ir);
    color=color0(ir);
    
    % 1st smoothing
    isos=smooth(stdmag,iso,0.1,'rlowess'); % 1.75 sec in my laptop
    
    % remove 3 sigma outliers, twice; 4-sigma, once
    rms2=(iso-isos).^2; % rms^2
    rms2s=std(rms2);    % standard deviation of rms^2
    ix=find(rms2<3^2*rms2s);  % remove 3-sigma ouliters
    iso3=iso(ix);
    iso3s=isos(ix);
    stdmag3=stdmag(ix);
    color3=color(ix);        
    
    rms3=(iso3-iso3s).^2;
    rms3s=std(iso3-iso3s).^2;
    ix=find(rms3<3^2*rms3s);     % remove 3-sigma ouliters   
    iso3=iso3(ix);
    iso3s=iso3s(ix);
    stdmag3=stdmag3(ix);
    color3=color3(ix);             
    
    rms3=(iso3-iso3s).^2;
    rms3s=std(iso3-iso3s).^2;
    ix=find(rms3<4^2*rms3s);     % remove 4-sigma ouliters   
    iso4=iso3(ix);
    iso4s=iso3s(ix);
    stdmag4=stdmag3(ix);
    color4=color3(ix);   
    
    % derive a roughly 3-sigma limiting mag
    rms4=(iso4-iso4s).^2;
    rms4s=smooth(stdmag4,rms4,0.2);
    rms0=(iso4-max(iso4s)).^2;
    rms0s=smooth(stdmag4,rms0,0.2);
    ix=find(rms4s>rms0s/16); % 16: 4 sigma limiting magnitude
    stdpointa=stdmag4(min(ix));
    isopointa=iso4s(min(ix));
    
    % only use stars brighter than limiting mag
    ix=find(stdmag4<stdpointa & iso4<isopointa+2*rms3s.^0.5);
    stdmag=stdmag4(ix);
    iso=iso4(ix);
    iso5s=iso4s(ix);
    color=color4(ix); 
    
    
    if size(stdmag(:),1)>100
        % CAUTION: systematic differences between rms(iso) and rms(stdmag)
        % for two dnr plates, when using rms(iso), the color-term is about -0.6 -0.7 
        
        
        % get rid of identical stdmag; if it is indentical, assign a random deviation to it
        % twice, to aviod fake randn produce same value
        dfi=diff(stdmag);
        k2=0;
        ix=find(dfi==0);
        stdmag(ix)=stdmag(ix)+randn(size(ix))*0.01;        
        dfi=diff(stdmag);
        k2=0;
        ix=find(dfi==0);
        stdmag(ix)=stdmag(ix)+randn(size(ix))*0.01;        
        
        isos=smooth(stdmag,iso,0.1); % to speed up, only use regular smooth; since already sigma-clipped, it is ok.
        
        
        % interp to get the daschmag for each iso value
        ix=find(iso<max(isos)&iso>min(isos));
        iso2=iso(ix);
        stdmag2=stdmag(ix);
        color2=color(ix);
        stdmags2=interp1(isos,stdmag,iso2,'cubic');
        dmag=stdmags2-stdmag2;
        myrms(1)=std(dmag); % Get an overall rms for cterm=0.
        cterm(1)=0;        
        
        coef=corr2(color2,dmag);
        if abs(coef)<0.5;
            cterm(2)=0;
            myrms(2)=myrms(1);
        else  % do a correlation to determine the second cterm value
            % linear fit; use bisector. Ref: Isobe et al. 1990, ApJ, 364, 104
            xdata=color2;
            ydata=dmag;
            xmean=mean(xdata);
            ymean=mean(ydata);
            
            Sxx=sum((xdata-xmean).^2);
            Syy=sum((ydata-ymean).^2);
            Sxy=sum((xdata-xmean).*(ydata-ymean));
            beta1=Sxy/Sxx; % OLS(X|Y)
            beta2=Syy/Sxy; % OLS(Y|X)
            
            beta3=(beta1*beta2 -1 + ((1+beta1^2)*(1+beta2^2))^0.5 )/(beta1+beta2);
            cterm(2)=beta3;
            
            refmag=stdmag+color*cterm(2);            
            isos=smooth(refmag,iso,0.1); % to speed up, only use regular smooth; since already sigma-clipped, it is ok.
            ix=find(iso<max(isos)&iso>min(isos));
            iso2=iso(ix);
            refmag2=refmag(ix);
            refmags2=interp1(isos,refmag,iso2,'cubic');
            dmag=refmags2-refmag2;
            myrms(2)=std(dmag); % Get an overall rms for cterm=beta3.
            
        end
        
        % Sampling: metropolis algorithm
        sigma=0.15;
        kdE=0.004;
        cterm0=cterm(2);
        myrms0=myrms(2);
        for km=3:50
            cterm(km)=cterm0+sigma*randn(1);
            refmag=stdmag+color*cterm(km);
            
            isos=smooth(refmag,iso,0.1); % to speed up, only use regular smooth; since already sigma-clipped, it is ok.
            
            ix=find(iso<max(isos)&iso>min(isos));
            iso2=iso(ix);
            refmag2=refmag(ix);
            refmags2=interp1(isos,refmag,iso2,'cubic');
            dmag=refmags2-refmag2;
            myrms(km)=std(dmag); % Get an overall rms for the given cterm
            
            dmyrms=myrms(km)-myrms0;
            if dmyrms<0
                cterm0=cterm(km);
                myrms0=myrms(km);
            else
                p=exp(-dmyrms/kdE);
                if rand(1)>p
                    cterm0=cterm(km-1);
                    myrms0=myrms(km-1);
                else
                    cterm0=cterm(km);
                    myrms0=myrms(km);
                end
            end
            
        end
        
        
        % fit the cterm vs myrms to derive the color-term
        minrms=min(myrms);
        mincterm=cterm(find(myrms==minrms));
        mincterm=mincterm(1); % to avoid crash by 2 minimum
        ix=find(abs(cterm-mincterm)<0.5 & myrms<1);   
        if (length(ix) > 2) 
          keff=size(ix(:),1);
          xdata=cterm(ix);
          ydata=myrms(ix);
          f=inline('c(1)*(x-c(2)).^2+c(3)','c','x');
          x0=[0.4,mincterm,minrms];
         [c,norm,resi,exitflag,output,lambda,jacobian]=lsqcurvefit(f,x0,xdata,ydata);
          
          x=xdata;
          covcolor=inv(jacobian'*jacobian);
          sigcolor=(norm/(keff-6)*covcolor).^0.5;
          errcolor=sigcolor(2,2);
          
          % good fit: points at both sides; c(1)>0; abs(c(2)-mincterm)<0.2; error<0.1;
          colorterm=c(2);
          kleft=size(find(cterm<colorterm),2);
          kright=size(find(cterm>colorterm),2);
          if min(kleft,kright)>3 & c(1)>0 & abs(colorterm-mincterm)<0.2 & errcolor<0.1
              flag=1;
          else
              flag=0;
              colorterm=cterm(2);
          end
          
          xfit=min(xdata):0.01:max(xdata)+0.01;
          yfit=c(1)*(xfit-c(2)).^2+c(3);
          
          eval('cd /dasch/junk/colorterm/octave2')
          pwd
          % h=axes('Fontsize',14);
          figure(1,'visible','off');
          newplot();
          plot(xdata,ydata,'b+',xfit,yfit,'r')
          xlabel('color-term','fontsize',18)
          ylabel('Photometric fitting RMS','fontsize',18)
          title([name ', colorterm = ' num2str(colorterm,3) ', flag=' num2str(flag)],'fontsize',16)
  %        saveas(gcf,[name '.colorterm.eps'],'epsc2')
          print([name '_a'  num2str(bin) '.out.colorterm.eps'])  % print -dpsc2 test     
          close

        else
           flag=0;
           errcolor = myrms(2);
           colorterm=cterm(2);
        endif

        y=[colorterm errcolor flag];
        eval('cd /dasch/junk/colorterm/octave2')
        pwd
        fid = fopen([name '_a'  num2str(bin)  '.colorterm.txt'],'w');    
        fprintf(fid,'%1.2f  %1.2f  %1.0f\n',y);
        fclose(fid);
        
        
        clear           covcolor    errcolor    iso2        k2          mincterm    refmag      rms4        stdmags2    
        clear REF         cterm       exitflag    iso3        kdE         minrms      refmag2     rms4s       stdpointa   
        clear c           cterm0      extinction  iso3s       keff        myrms       refmags2    rmsout      x           
        clear coef        ddec        f           iso4                    myrms0      resi        sigcolor    x0          
        clear color       dec         flag        iso4s       kleft              rms0        sigma       x_image     
        clear color0      dfi         flags       iso5s       km          norm        rms0s       stdmag      xdata       
        clear color2      dmag        h           isopointa   kright      number      rms2        stdmag0     xfit        
        clear color3      dmagout     ir          isos                output      rms2s       stdmag2     y_image     
        clear color4      dmyrms      iso         ix          lambda      p           rms3        stdmag3     ydata       
        clear colorterm   dra         iso0        jacobian    magout      ra          rms3s       stdmag4     yfit        
        clear numpoints iii fid count txt y
        
        %who -variables
    end
  end
end

toc

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
% Jan 12, 2009 Edward J Los - Integrate this script into the pipeline
% Jan 20, 2009 Edward J Los - Correct errors for insufficient bin points
% Feb  3, 2009 Edward J Los - Add an extra decimal point to the colorterm output
% Mar 18, 2009 Edward J Los - Perform colorterm calculations for bins 1-3, 4-6, and 7-9
%                           - Reject any Skymap and Tycho 2 objects
%                           - Correct tests in the length of ix

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
do_plots = getenv("DASCH_PLOT");
printf("do_plots is %s\n",do_plots);
if (length(do_plots) < 1) 
  exit;
endif



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
  
    chdir(catalogbin)
    pwd();


for kkk=1:ksize
    name=char(B(kkk));
    for bin=1:9 
       outcmd = ['rm -f ' name qualifier '_a'  num2str(bin) '.xcolorterm.eps'];
       system(outcmd);
       outcmd = ['rm -f ' name qualifier '_a'  num2str(bin)  '.xcolorterm.txt' ];
       system(outcmd);

    end

  for gbin=1:3:9 % 1:9
    clear input inputREF
    input = [];
    for sbin=gbin:gbin+2
        printf("colorterm %s%s sbin %d\n",name,qualifier,sbin);
        fid = fopen([name qualifier '_a' num2str(sbin) debugSuffix2 '.db'],'r');
        if (fid == -1)
          printf("Failed to open %s\n",[name qualifier '_a' num2str(sbin) debugSuffix2 '.db']);
          continue;
        else
          fclose(fid);
        endif
        load([name qualifier '_a' num2str(sbin) debugSuffix2 '.db'],['output' num2str(sbin) 'REF'],['output' num2str(sbin)]);
        command = ['sinput = output' num2str(sbin) ';'];
        eval(command);
        clear(['output' num2str(sbin)]);
        command = ['sinputREF = output' num2str(sbin) 'REF;'];
        eval(command);
        clear(['outputREF' num2str(sbin)]);
        %printf("length is %d for bin %d\n",length(sinput),sbin);

        if (length(input) == 0)
          input = sinput;
          inputREF = sinputREF;
        else
          lengthold = length(input);
          lengthnew = length(sinput);
          if (lengthnew > 0) 
            input(lengthold+1:lengthold+lengthnew,:) = sinput(1:lengthnew,:);
            inputREF(lengthold+1:lengthold+lengthnew) = sinputREF(1:lengthnew);
          endif
        endif
        clear lengthold lengthnew sinput sinputREF
      end
    
      
      if (length(input) == 0) 
        printf("consolidated length is %d for %s%s bin %d\n",length(input),name,qualifier,gbin);
        continue
      endif

      stdmag = input(:,4);
      color = input(:,5);
      iso = input(:,6);
      ddec = input(:,11);
      dra = input(:,12);
      bflags = input(:,8);
    
    % For MS stars, B-R is about [-0.5,2];
    % only use stars with stdmag between 4 and 19 (to avoid crazy magnitudes)
    % and use drad filter - use absolute value!
    % The bitand operation excludes Tycho2 and Skymap stars from the test set.  See pipelineutils.h
%    ix=find((bitand(bflags,61440) == 0) & (stdmag>4) & (stdmag<19) & (color>-0.5) & (color<2) & (abs(ddec)<3*median(abs(ddec))) & (abs(dra)<3*median(abs(dra))));  
    printf("ERROR: reversed sense of bitand(bflags,61440)\n");
    ix=find((bitand(bflags,61440) != 0) & (stdmag>4) & (stdmag<19) & (color>-0.5) & (color<2) & (abs(ddec)<3*median(abs(ddec))) & (abs(dra)<3*median(abs(dra))));  

    if (length(ix) <= 0) 
      printf("consolidated length ix (1) is %d for %s%s bin %d\n",length(ix),name,qualifier,gbin);
    endif

    if (length(ix) > 0)
    
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
    endif
    if (length(ix) > 0)
      iso3=iso(ix);
      iso3s=isos(ix);
      stdmag3=stdmag(ix);
      color3=color(ix);        
    
      rms3=(iso3-iso3s).^2;
      rms3s=std(iso3-iso3s).^2;
      ix=find(rms3<3^2*rms3s);     % remove 3-sigma ouliters   
    endif
    if (length(ix) > 0)
      iso3=iso3(ix);
      iso3s=iso3s(ix);
      stdmag3=stdmag3(ix);
      color3=color3(ix);             
    
      rms3=(iso3-iso3s).^2;
      rms3s=std(iso3-iso3s).^2;
      ix=find(rms3<4^2*rms3s);     % remove 4-sigma ouliters   
    endif
    if (length(ix) > 0)
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
    endif
    if (length(ix) > 0)
      stdpointa=stdmag4(min(ix));
      isopointa=iso4s(min(ix));
    
      % only use stars brighter than limiting mag
      ix=find(stdmag4<stdpointa & iso4<isopointa+2*rms3s.^0.5);
    endif
    if (length(ix) > 0)
      stdmag=stdmag4(ix);
      iso=iso4(ix);
      iso5s=iso4s(ix);
      color=color4(ix); 
    else
      stdmag=stdmag(ix);
    endif
    
    if (size(stdmag(:),1) <= 100) 
      printf("consolidated length stdmag is %d length ix %d for %s%s bin %d\n",length(stdmag),length(ix),name,qualifier,gbin);
    endif    

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
        
        if (length(ix) <= 2)
          printf("consolidated length second ix is %d for %s%s bin %d\n",length(ix),name,qualifier,gbin);
        endif

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
              flag=2;
          else
              flag=1;
              colorterm=cterm(2);
          end
          
          if strcmp(do_plots,"YES") 
            xfit=min(xdata):0.01:max(xdata)+0.01;
            yfit=c(1)*(xfit-c(2)).^2+c(3);
            
            % h=axes('Fontsize',14);
            figure(1,'visible','off');
            newplot();
            plot(xdata,ydata,'b+',xfit,yfit,'r')
            xlabel('color-term','fontsize',18)
            ylabel('Photometric fitting RMS','fontsize',18)
            title([name qualifier ', colorterm = ' num2str(colorterm,3) ', flag=' num2str(flag)],'fontsize',16)
    %        saveas(gcf,[name qualifier '.colorterm.eps'],'epsc2')
            print([name qualifier '_a'  num2str(gbin) '.xcolorterm.eps'])  % print -dpsc2 test     
            close
          %printf("Writing %s\n",[name qualifier '_a'  num2str(gbin) '.xcolorterm.eps']); 
          endif
        else
           flag=1;
           errcolor = myrms(2);
           colorterm=cterm(2);
        endif

        y=[colorterm errcolor flag];
        fid = fopen([name qualifier '_a'  num2str(gbin)  '.xcolorterm.txt'],'w');    
        fprintf(fid,'%1.3f  %1.3f  %1.0f\n',y);
        fclose(fid);
        
    end
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
    clear numpoints iii fid txt y
    clear bflags
        
    %who -variables
  end
end

toc

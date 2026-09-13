% Copyright the President and Fellows of Harvard College.
% Licensed under the MIT License

% Usage:
%
%   annular9.m <solutionNumber> <list> [qualifier]
%
% Inputs:
%
%   <listfile>  --  list of {name}, one per line; must be absolute path
%   $CATALOGBIN/{name}{solstr}{qual}.colorterm_marker.txt
%   $CATALOGBIN/{name}{solstr}{qual}_a{bin}.colorterm.txt
%   $CATALOGBIN/{name}{solstr}{qual}_a{bin}.db
%
% Outputs:
%
%   $BINOUTPUT/{name}{solstr}{qual}_a{bin}.grid
%   $BINOUTPUT/{name}{solstr}{qual}_a{bin}.para
%   $BINOUTPUT/{name}{solstr}{qual}_a{bin}.out
%   $BINOUTPUT/{name}{solstr}{qual}_a{bin}_a.eps  (maybe)
%   $BINOUTPUT/{name}{solstr}{qual}_a{bin}_b.eps  (maybe)
%   $BINOUTPUT/{name}{solstr}{qual}.annular9_marker.txt
%
% Environment:
%
%   $DASCH_BINOUTPUT
%   $DASCH_CATALOGBIN
%   $DASCH_NUMBINS  --  1 or 9
%   $DASCH_PLOT  --  either "YES" or not
%   $DASCH_SCRIPTS  -- for getting smooth.m

tic
clear

startTime = time();

pkg load image
pkg load struct
pkg load io
pkg load statistics
pkg load optim

debugSuffix = "";
debugSuffix2 = "";

smooth_colorterm = 0;
if (smooth_colorterm == 1)
  printf("ERROR: smooth_colorterm is enabled\n");
end

scripts = getenv("DASCH_SCRIPTS");
printf("scripts is %s\n",scripts);
if (length(scripts) < 1)
  exit;
end

catalogbin = getenv("DASCH_CATALOGBIN");
printf("catalogbin is %s\n",catalogbin);
if (length(catalogbin) < 1)
  exit;
end

binoutput = getenv("DASCH_BINOUTPUT");
printf("binoutput is %s\n",binoutput);
if (length(binoutput) < 1)
  exit;
end

do_plots = getenv("DASCH_PLOT");
printf("do_plots is %s\n",do_plots);
if (length(do_plots) < 1)
  exit;
end

num_bins = getenv("DASCH_NUMBINS");
printf("num_bins is %s\n",num_bins);
if (length(num_bins) < 1)
  exit;
end

qualifier = "";

addpath(scripts);

chdir(scripts);
pwd();    % this prevents a synchronization problem.

if (nargin() < 2)
  printf("Usage: octave annular9.m <solutionNumber> <list> [qualifier]\n");
  exit;
end

solutionNumber = argv()(1);
listfile = argv()(2);
if (nargin() == 3)
  qualifier =  ['_' char(argv()(3))];
end

if (char(solutionNumber) == "0")
  solutionString = "";
else
  solutionString = ['_s' char(solutionNumber)];
end

printf("List file is %s qualifier is %s%s\n",char(listfile),solutionString,qualifier);

fid = fopen(char(listfile),'rt');
iii = 1;
B = cell(1);

while ((txt = fgetl(fid)) != -1)
  B(iii) = txt;
  iii++;
end

fclose(fid);

ksize=length(B);

for kkk=1:ksize
  mosaicStartTime = time();
  name=char(B(kkk));

  % Initialize the RNG explicitly so that we can provide repeatable results.
  % Based on experimentation, the different rand*() functions seem to each
  % maintain their own random state, so you can't (safely) just seed one. For
  % copy-paste friendliness, we initialize all of them.
  rand("state", unicode2native(name, "UTF-8"));
  randn("state", unicode2native(name, "UTF-8"));
  rande("state", unicode2native(name, "UTF-8"));
  randg("state", unicode2native(name, "UTF-8"));
  randp("state", unicode2native(name, "UTF-8"));

  chdir(binoutput)
  pwd();

  lastcolorterm = 0;
  lasterrorcolor = 99;
  failcolorflag=0;

  [_err, _msg] = unlink([name solutionString qualifier '.annular9_marker.txt']);

  for bin=1:9
    [_err, _msg] = unlink([name solutionString qualifier '_a' num2str(bin) debugSuffix '.para']);
    [_err, _msg] = unlink([name solutionString qualifier '_a' num2str(bin) debugSuffix '.out']);
    [_err, _msg] = unlink([name solutionString qualifier '_a' num2str(bin) debugSuffix '.grid']);
    [_err, _msg] = unlink([name solutionString qualifier '_a' num2str(bin) debugSuffix '_a.eps']);
    [_err, _msg] = unlink([name solutionString qualifier '_a' num2str(bin) debugSuffix '_b.eps']);
  end

  chdir(catalogbin)
  pwd();

  fid = fopen([name solutionString qualifier '.colorterm_marker.txt'],'r');
  if (fid == -1)
    printf("ERRORD: Failed to open colorterm marker %s\n",[name solutionString qualifier '.colorterm_marker.txt']);
  else
    %printf("Found colorterm marker %s\n",[name solutionString qualifier '.colorterm_marker.txt']);
    fclose(fid);
    [_err, _msg] = unlink([name solutionString qualifier '.colorterm_marker.txt']);
  end

  for bin=1:9 % 1:9
    binStartTime = time();
    chdir(catalogbin)
    pwd();

    % Get the colorterm correction
    if (smooth_colorterm == 1)
      fid = fopen([name solutionString qualifier '_a'  num2str(bin)  '.smoothcolorterm.txt'],'r');
    else
      fid = fopen([name solutionString qualifier '_a'  num2str(bin)  '.colorterm.txt'],'r');
    end

    colorterm = lastcolorterm;
    errorcolor = lasterrorcolor;
    colorflag=failcolorflag;

    if (fid == -1)
      if (smooth_colorterm == 1)
        printf("WARNING: Failed to open %s\n",[name solutionString qualifier '_a'  num2str(bin)  '.smoothcolorterm.txt']);
      else
        printf("WARNING: Failed to open %s\n",[name solutionString qualifier '_a'  num2str(bin)  '.colorterm.txt']);
      end
    else
      txt = fgetl(fid);

      if (txt != -1)
        if (smooth_colorterm == 1)
          [oldcolorterm,errorcolor,colorflag,colorterm,count] = sscanf(txt,"%f %f %d %f","C");
          if (count != 4)
            printf("Decode count only %d from %s\n",count,[name solutionString qualifier '_a'  num2str(bin)  '.colorterm.txt']);
            colorterm = lastcolorterm;
            errorcolor = lasterrorcolor;
            colorflag=failcolorflag;
          else
            failcolorflag=3; % successful read, use this value for any subsequent failed bins
          end
        else
          [colorterm,errorcolor,colorflag,count] = sscanf(txt,"%f %f %d","C");
          if (count != 3)
            printf("Decode count only %d from %s\n",count,[name solutionString qualifier '_a'  num2str(bin)  '.colorterm.txt']);
            colorterm = lastcolorterm;
            errorcolor = lasterrorcolor;
            colorflag=failcolorflag;
          else
            failcolorflag=3; % successful read, use this value for any subsequent failed bins
          end
        end
      else
        printf("Failed to read a line from %s\n",[name solutionString qualifier '_a'  num2str(bin)  '.colorterm.txt']);
      end

      fclose(fid);
    end

    lastcolorterm = colorterm;
    lasterrorcolor = errorcolor;

    fid = fopen([name solutionString qualifier '_a' num2str(bin) debugSuffix2 '.db'],'r');
    if (fid == -1)
      printf("Failed to open %s\n",[name solutionString qualifier '_a' num2str(bin) debugSuffix2 '.db']);
      continue;
    else
      fclose(fid);
    end

    load([name solutionString qualifier '_a' num2str(bin) debugSuffix2 '.db'],['output' num2str(bin) 'REF'],['output' num2str(bin)]);
    command = ['input = output' num2str(bin) ';'];
    eval(command);
    clear(['output' num2str(bin)]);
    command = ['inputREF = output' num2str(bin) 'REF;'];
    eval(command);
    clear(['outputREF' num2str(bin)]);

    chdir(binoutput)
    pwd();

    % Adjust stdmag for color
    if (colorflag > 0)
      % Adjust color entries by COLOR_RATIO when VMag is substituted for FpgMag
      % ix = find(bitand(input(:,8),32768) != 0);
      % if (length(ix) > 0)
      %   input(ix,5) = input(ix,5) * COLOR_RATIO;
      % end

      ix = find(input(:,5) < 90.0); % ignore erroneous colors
      if (length(ix) > 0)
        input(ix,4) = input(ix,4) + (input(ix,5)*colorterm);
      end
    end

    % ------ input1: remove stdmag>19 (99.9), and the bad points at the bright end
    ksize1 =  size(input,1);
    stdpointa = 0;
    refmag4crit = 0;

    if (ksize1 > 10)
      [sref,ix]=sort(input(:,4));
      input=input(ix,:);
      inputREF=inputREF(ix,:);
      mref=input(2,4);
      mref2=input(4,4);
      ksize1 = length(input(:,4));
    end

    if (ksize1 > 10)
      ix2=find(input(1:10,6)>(min(input(:,6))+max(input(:,6)))/3);
      input(ix2,6)=99;
      ix=find(input(:,4)<19&input(:,4)>mref2-2 &input(:,4)>mref-1 & input(:,6)<0);
      kxa=size(ix,1);
      input1(1:kxa,:)=input(ix(1:kxa),:);
      input1REF(1:kxa,:)=inputREF(ix(1:kxa),:);

      % ------ input2: sort by stdmag
      [s,ir]=sort(input1(:,4));  % sort by stdmag
      input2=input1(ir,:);
      input2REF=input1REF(ir,:);
      refmag2=input2(:,4);       % refmag is stdmag in this whole script
      iso2=input2(:,6);
      clear input input1
      clear inputREF input1REF

      % -------- use mag as x, iso as y
      isosa = smooth(refmag2, iso2, 0.25, 'rlowess');
      % plot(refmag2,iso2,'b.','markersize',5)
      % hold on

      size2 = size(iso2, 1);
      size1 = round(size2 * 0.1);

      if size1 <= 1
        printf("setting size1 %d to 2\n", size1);
        size1 = 2;
      end

      isosb = smooth(refmag2(1:size1), iso2(1:size1), 0.1, 'rlowess');

      size3 = round(size2 * 0.4);
      isosc = smooth(refmag2(1:size3), iso2(1:size3), 0.2, 'rlowess');

      iso2s=[
        isosb(1:ceil(size1/8)); ...
        isosb(ceil(size1/8) + 1:ceil(size1/6)) * 0.5 + isosc(ceil(size1/8) + 1:ceil(size1/6)) * 0.5; ...
        isosc(ceil(size1/6) + 1:ceil(size1/2)); ...
        isosc(ceil(size1/2) + 1:ceil(size3/3)) * 0.5 + isosa(ceil(size1/2) + 1:ceil(size3/3)) * 0.5; ...
        isosa(ceil(size3/3) + 1:size2)
      ];

      iso3s = smooth(refmag2, iso2s, 0.03);

      ix = find(abs(iso3s - iso2s) > median(abs(iso3s - iso2s)) * 4);  % revised - throw out delta(iso)>median*4
      iso2s(ix) = iso3s(ix);
      iso2s(ceil(size1 / 7):ceil(size3 / 2)) = iso3s(ceil(size1 / 7):ceil(size3 / 2));
      iso4s = smooth(refmag2, iso2s, 0.04);

      iso2s(ceil(size1/4):ceil(size3/2)) = iso4s(ceil(size1/4):ceil(size3/2));
      iso5s = smooth(refmag2, iso2s, 0.04);

      iso2s(ceil(size1/4):ceil(size3/2)) = iso5s(ceil(size1/4):ceil(size3/2));
      iso6s = smooth(refmag2, iso2s, 0.03);

      iso2s(ceil(size1/8):ceil(size3/2)) = iso6s(ceil(size1/8):ceil(size3/2)); % Bugfix of Apr 26, 2011
      iso2s = smooth(refmag2, iso2s, 0.02);

      iso2s = smooth(refmag2, iso2s, 0.02);  % revised - double smooth to restrain non-monotonicity ripples

      % ------ determine the critical point a
      %  remove 4 sigma outliers
      rms2=(iso2-iso2s).^2;
      rms2s=std(iso2-iso2s).^2;  % standard deviation
      ix=find(rms2<4^2*rms2s);

      if (size(ix,1) > 2)
        input3=input2(ix,:);
        refmag3=refmag2(ix);
        iso3=iso2(ix);
        iso3s=iso2s(ix);

        %  remove 4 sigma outliers again
        rms3=(iso3-iso3s).^2;
        rms3s=std(iso3-iso3s).^2;
        ix=find(rms3<4^2*rms3s);
        input4=input3(ix,:);
        refmag4=refmag3(ix);
        iso4=iso3(ix);
        iso4s=iso3s(ix);

        % determine the flux limit
        rms4=(iso4-iso4s).^2;
        rms4s=smooth(refmag4,rms4,0.2);
        % Remove: rms0=(iso4-max(iso4s)).^2;
        miniso4s = min(iso4s);
        ix =  find(iso4s < (miniso4s+0.01));
        maxix = max(ix);
        maxiso4s = max(iso4s(maxix:size(iso4s,1)));
        rms0=(iso4-maxiso4s).^2;

        rms0s=smooth(refmag4,rms0,0.2);
        ix=find(rms4s>rms0s/4); % 4: 2 sigma limiting magnitude

        if (size(ix,1) > 0)
          stdpoint1=refmag4(min(ix));

          iy1=find(rms4s>=min(max(rms4s),5*median(rms4s)) & refmag4>(stdpoint1+min(refmag4))/2);
          if (size(iy1,1) > 0)
            stdpointb=refmag4(min(iy1));
          else
            stdpointb= stdpoint1;
          end

          iy2=find(rms4s>3*median(rms4s) & refmag4>(stdpoint1+min(refmag4))/2);
          ky2=ceil(size(iy2(:),1)/2);

          if ky2>1
            stdpointc=median(refmag4(iy2));

            lowlimit = round(median(iy2))-ky2;
            if (lowlimit < 1)
              lowlimit = 1;
            end

            highlimit = round(median(iy2))+ky2;
            if (highlimit > size(refmag4,1))
              highlimit = size(refmag4,1);
            end

            refmagtest=refmag4(lowlimit:highlimit);
            rmstest=rms4s(lowlimit:highlimit);
            iy2=find(rms4s>=max(rmstest) & refmag4>(stdpoint1+min(refmag4))/2);

            if (length(iy2) > 0)
              stdpointc=refmag4(min(iy2));
            else
              stdpointc=stdpointb;
            end
          else
            stdpointc=stdpointb;
          end

          % sort iso4  (Begin revision of Oct 20, 2008)
          [iso4sort,ir]=sort(iso4);
          refmag4sort=refmag4(ir);
          % get the median refmag value for the faintest objects
          % to avoid wrong limiting mag when there are very few points at the faint end
          size4=size(iso4(:),1);

          if (size4 > 10)
            ix=find(iso4>iso4sort(size4-10)-0.6);
            stdpointd=median(refmag4(ix));   % median refmag value for the faintest objects in iso4 with iso>max(iso)-0.6
            stdpointe=median(refmag4sort(size4-floor(size4/100)-10:size4));  % median refmag for the faintest 1% objects

            stdpoint2=min(min(stdpointb,stdpointc), min(stdpointd,stdpointe));
            stdpointa=min(stdpoint1,stdpoint2);
          else
            stdpointa=stdpoint1;
          end

          isopointa=min(iso4s(find(abs(refmag4-stdpointa)==min(abs(refmag4-stdpointa)))));
          rms4smax = max(rms4s);

          clear stdpoint1 stdpoint2 stdpointb stdpointc stdpointd stdpointe size4 iy1 iy2 refmagtest rmstest lowlimit highlimit maxix maxiso4s miniso4s
          % (end revision of Oct 20, 2008)

          if (stdpointa > 0)
            ix=find(refmag4<stdpointa);
            ksize1=size(ix(:),1);
          else
            ksize1=-1;
          end
        else
          ksize1=-8;
        end
      else
        ksize1=-7;
      end
    else
      ksize1=-2;
      refmag2 = 0; % avoid an exception when printing statistics
    end

    if (ksize1 > 100)
      % ---- to avoid non-monotonicity
      iso3s=smooth(refmag4,iso4s,0.03);
      ix=find(abs(iso3s-iso4s)>median(abs(iso3s-iso4s))*4);  % throw out delta(iso)>median*4
      iso4s(ix)=iso3s(ix);

      % Abandon all data for reference magnitudes brighter than the minimum iso4s

      iso4sminimum = min(iso4s);
      ix = find(iso4s == iso4sminimum);
      ixmax = ix(length(ix));
      if (ixmax <= size(iso4s(:),1))
        iso4s=iso4s(ixmax:size(iso4s(:),1));
        refmag4=refmag4(ixmax:size(refmag4(:),1));
      end

      clear ixmax iso4sminimum;

      ix=find([diff(iso4s);1]<=0 & refmag4<stdpointa);  % only consider points brighter than limiting magnitude
      if (length(ix) > 0)
        diffix=diff(ix);

        % if there is negative slope in the bright end, abandon the data
        if ix(1)==1
          ix3=find(diffix>1);

          if (length(ix3) > 0)
            iso4s=iso4s(ix3(1)+1:size(iso4s(:),1));
            refmag4=refmag4(ix3+1:size(refmag4(:),1));
          else
            iso4s=iso4s(length(ix)+1:size(iso4s(:),1));
            refmag4=refmag4(length(ix)+1:size(refmag4(:),1));
          end
        end
      end

      ix=find([diff(iso4s);1]<=0 & refmag4<stdpointa);  % only consider points brighter than limiting magnitude
      kx=size(ix(:),1);

      if kx>2
        iso3s=smooth(refmag4,iso4s,0.04);
        ix=find(abs(iso3s-iso4s)>median(abs(iso3s-iso4s))*3.5);  % throw out delta(iso)>median*3.5
        iso4s(ix)=iso3s(ix);
      end

      ix=find([diff(iso4s);1]<=0 & refmag4<stdpointa);  % only consider points brighter than limiting magnitude
      kx=size(ix(:),1);
      diffix=diff(ix);
      ix2=find(diffix==1);
      kx2=size(ix2(:),1);

      % interpolate the ripple points with the value of smoothed curve isoss
      loopcounter = 0;

      while ((kx2>2) && (loopcounter < 100))
        loopcounter = loopcounter+1;
        nr=0;   % number of ripples
        nr2=0;
        krnum=1;
        clear kstart kend

        for i=1:kx-2
          if (diffix(i)==1) && (diffix(i+1)>1)
            nr=nr+1;   %  number of ripples
            kend(nr)=ix(i+1);    % end ix of the ripple
          end

          if diffix(1)>1
            if (diffix(i)>1) && (diffix(i+1)==1)
              nr2=nr2+1;
              kstart(nr2)=ix(i+1);   %  start ix of the ripple
            end
          elseif diffix(1)==1
            kstart(1)=ix(1);
            if (diffix(i)>1) && (diffix(i+1)==1)
              nr2=nr2+1;
              kstart(nr2+1)=ix(i+1);  %   start ix of the ripple
            end
          end
        end

        nr3=min(nr,nr2);

        if nr3>0
          kstart=kstart(1:nr3);
          kend=kend(1:nr3);
          klength_half=ceil((kend-kstart+1)/2);   % the half length of the ripple

          for i=1:nr3
            ix=find(abs(iso4s-iso4s(kstart(i)))==min(iso4s-iso4s(kstart(i))));

            if size(ix(:),1)>0
              if (ix(1)>kend) && (refmag4(ix(1))<stdpointa)
                kstart0(i)=ix(1)+1;
              else
                kstart0(i)=max(kstart(i)-klength_half(i),1);
              end
            else
              kstart0(i)=max(kstart(i)-klength_half(i),1);
            end

            ix=find(abs(iso4s-iso4s(kend(i)))==min(iso4s-iso4s(kend(i))));

            if size(ix(:),1)>0
              if (ix(1)<kstart) && (refmag4(ix(1))<stdpointa)
                kend0(i)=ix(1);
              else
                kend0(i)=min(kend(i)+klength_half(i),size(refmag4(:),1));
              end
            else
              kend0(i)=min(kend(i)+klength_half(i),size(refmag4(:),1));
            end

            iso4s(kstart0(i):kend0(i))=iso4s(kstart0(i))+(iso4s(kend0(i))-iso4s(kstart0(i)))*(0:1:(kend0(i)-kstart0(i)))/(kend0(i)-kstart0(i));
          end
        end

        iso3s=smooth(refmag4,iso4s,0.03);
        ix=find(abs(iso3s-iso4s)>median(abs(iso3s-iso4s))*3);  % revised - throw out delta(iso)>median*3
        iso4s(ix)=iso3s(ix);

        ix=find([diff(iso4s);1]<0 & refmag4<stdpointa);  % only consider points brighter than limiting magnitude
        kx=size(ix(:),1);
        diffix=diff(ix);
        ix2=find(diffix==1);
        kx2=size(ix2(:),1);
      end

      if (loopcounter >= 100)
        ksize1 = -6; % flag for infinite loop
      end

      clear iso3s ix kend krnum kx nr2 diffix ix2 kend0 kstart kx2 nr3 ix3 klength_half kstart0 nr
      clear loopcounter;

      % still define these two variables in order to be consistent with the later part of the script
      refmag4crit=min(refmag4);
      iso4scrit=min(iso4s);
    end

    % refmag4, iso4, iso4s
    if (ksize1>100) && ((stdpointa-refmag4crit)>1.5)
      % ---- output grid, bright end extrapolate to 6, faint end extrapolate to max(iso)

      % 1. slope of extrapolate the bright end
      %             magmin=min(refmag4);
      %             magmax=min(refmag4)+1.5;
      magmin=refmag4crit;    % use refmag4crit instead of min(refmag4) to guarantee monotonicity
      magmax=refmag4crit+1.5;
      ixmin=find(abs(refmag4-magmin)==min(abs(refmag4-magmin)));
      isomin1=iso4s(ixmin(1));
      magmin1=refmag4(ixmin(1));
      ixmax=find(abs(refmag4-magmax)==min(abs(refmag4-magmax)));
      isomax1=iso4s(ixmax(1));
      magmax1=refmag4(ixmax(1));

      % Make sure that we do not have the same point for maximum and minimum

      index1 = ixmax(1);
      while ((magmax1 == magmin) && (ixmax < length(refmag4)))
        index1++;
        isomax1 = iso4s(index1);
        magmax1 = refmag4(index1);
      end

      if (magmax1 == magmin1)
        ksize1 = -4;  % Can not go on without a divide by zero error
      end
    end

    if (ksize1>100) && ((stdpointa-refmag4crit)>1.5)
      slope1=(isomax1-isomin1)/(magmax1-magmin1);

      % 2. slope of extrapolate the faint end
      magmin=stdpointa-0.3;
      magmax=stdpointa+0.2;
      ixmin=find(abs(refmag4-magmin)==min(abs(refmag4-magmin)));
      isomin1=iso4s(ixmin(1));
      magmin1=refmag4(ixmin(1));
      ixmax=find(abs(refmag4-magmax)==min(abs(refmag4-magmax)));
      isomax1=iso4s(max(ixmax));
      magmax1=refmag4(max(ixmax));

      if (magmax1 == magmin1)
        magmin=stdpointa-0.75;
        magmax=stdpointa;
        ixmin=find(abs(refmag4-magmin)==min(abs(refmag4-magmin)));
        isomin1=iso4s(ixmin(1));
        magmin1=refmag4(ixmin(1));
        ixmax=find(abs(refmag4-magmax)==min(abs(refmag4-magmax)));
        isomax1=iso4s(ixmax(1));
        magmax1=refmag4(ixmax(1));
      end

      if (magmax1 == magmin1)
        ksize1 = -5;  % Can not go on without a divide by zero error
      end
    end

    if (ksize1>100) && ((stdpointa-refmag4crit)>1.5)
      slope3=(isomax1-isomin1)/(magmax1-magmin1);

      % ---- grid
      gridmax=stdpointa+slope3*(max(iso2)-isopointa);
      gridmin=min(6,min(refmag4)-0.1);
      maggrid=gridmin:0.01:gridmax;

      % ----------- bright part
      %            ix1=find(maggrid<min(refmag4));
      ix1=find(maggrid<refmag4crit);
      kx1=length(ix1);
      %            isogrid(1:kx1)=min(iso4s)+slope1*(maggrid(ix1)-min(refmag4));
      isogrid(1:kx1)=iso4scrit+slope1*(maggrid(ix1)-refmag4crit);
      flaggrid(1:kx1)=0;  % 0 means bright end extrapolate
      % plot(maggrid(1:kx1),isogrid,'g')

      % ---------- interp1 the central part
      %            ix2=find(maggrid>=min(refmag4)&maggrid<=stdpointa);
      ix2=find(maggrid>=refmag4crit&maggrid<=stdpointa);
      kx2=length(ix2);
      % isogrid(kx1+1:kx1+kx2)=interp1(refmag4,iso4s,refmag4(kx1+1:kx1+kx2),'cubic');
      % ---- remove the non-distinct part
      dfi=diff(refmag4);
      k2=0;

      for i=1:size(dfi(:),1)  % remove identical iso, interp1: The data abscissae should be distinct
        if dfi(i)>0
          k2=k2+1;
          xi3(k2)=refmag4(i);
          yi3(k2)=iso4s(i);
        end
      end

      if k2>10
        isogrid(kx1+1:kx1+kx2)=interp1(xi3,yi3,maggrid(kx1+1:kx1+kx2),'cubic','extrap');
        isogrid(kx1+1:kx1+kx2)=smooth(maggrid(kx1+1:kx1+kx2),isogrid(kx1+1:kx1+kx2));
        flaggrid(kx1+1:kx1+kx2)=1;
      else
        isogrid(kx1+1:kx1+kx2)=maggrid(kx1+1:kx1+kx2);
        flaggrid(kx1+1:kx1+kx2)=-2;     % -2 means no points in the central part
      end

      clear xi3 yi3
      % plot(maggrid(kx1+1:kx1+kx2),isogrid(kx1+1:kx1+kx2),'k')

      isogrid(1:kx1)=isogrid(kx1+1)+slope1*(maggrid(1:kx1)-maggrid(kx1+1));

      % ----------- faint part
      ix3=find(maggrid>stdpointa);
      kx3=length(ix3);
      isogrid(kx1+kx2+1:kx1+kx2+kx3)=isopointa+slope3*(maggrid(ix3)-stdpointa);
      flaggrid(kx1+kx2+1:kx1+kx2+kx3)=-1;  % -1 means faint end extrapolate
      % plot(maggrid(kx1+kx2+1:kx1+kx2+kx3),isogrid(kx1+kx2+1:kx1+kx2+kx3),'c')

      isogridsmooth=smooth(maggrid,isogrid);
      ix=find(abs(isogridsmooth-isogrid')<max(abs(isogridsmooth-isogrid'))/1.5);  % ------throw out delta(iso)>max/1.5
      isogridnew=interp1(maggrid(ix),isogrid(ix),maggrid,'cubic','extrap');
      isogrid=isogridnew;
      clear isogridnew isogridsmooth

      % add additional judgement for turn over point (begin revision of Oct 20, 2008)
      ix=find(maggrid>stdpointa-2 & maggrid<stdpointa); % consider the faint end 2 mag brighter than the previous limiting mag
      maggrid2=maggrid(ix);
      isogrid2=isogrid(ix);
      isodiff=diff(diff(isogrid2));
      ix=find(isodiff<0); % turn over point

      if size(ix(:),1)>0
        turnpoint=maggrid2(ix(1));
        ix=find(abs(maggrid2-turnpoint)<0.5);
        maggrid3=maggrid2(ix);
        isogrid3=isogrid2(ix);
        slope4=(max(isogrid3)-min(isogrid3))/(max(maggrid3)-min(maggrid3));
        isogrid4=min(isogrid3)+slope4*(maggrid2-min(maggrid3));
        ix2=find(isogrid4-isogrid2>0.1); % turn over at least 0.1 in iso

        if size(ix2(:),1)>1
          stdpointa2=maggrid2(ix2(1));
          stdpointa=min(stdpointa, stdpointa2);
          % re-interpolate the faint part after getting the new stdpointa
          ix3=find(maggrid>stdpointa);
          kx3=length(ix3);
          kx2 = length(isogrid) - kx1 -kx3;
          isogrid(kx1+kx2+1:kx1+kx2+kx3)=isopointa+slope3*(maggrid(ix3)-stdpointa);
          flaggrid(kx1+kx2+1:kx1+kx2+kx3)=-1;  % -1 means faint end extrapolate

          isogridsmooth=smooth(maggrid,isogrid);
          ix=find(abs(isogridsmooth-isogrid')<max(abs(isogridsmooth-isogrid'))/1.5);  % ------throw out delta(iso)>max/1.5
          isogridnew=interp1(maggrid(ix),isogrid(ix),maggrid,'cubic','extrap');
          isogrid=isogridnew;
          clear isogridnew isogridsmooth
        end
      end

      clear isodiff isogrid2 isogrid4 ix2 maggrid3 turnpoint isogrid3 ix maggrid2 slope4 stdpointa2
      % (end revision of Oct 20, 2008)

      % ----- grid done except error: maggrid, isogrid, flaggrid

      % - Do the output file: iso2, fitted mag
      isoout=iso2;
      magout=interp1(isogrid,maggrid,isoout,'cubic','extrap');
      % plot(magout,isoout,'r.')
      dmagout=magout-refmag2;
      myout1=[isoout,magout,dmagout];
      % ---- output file done

      % remove flux limit and 5 sigma outliers
      ix=find(magout<=stdpointa & dmagout.^2<5^2*std(dmagout).^2);
      myout2=myout1(ix,:);
      % remove 5 sigma outliers again
      ix=find(myout2(:,3).^2<5^2*std(myout2(:,3))^2);
      myout3=myout2(ix,:);

      if (length(myout3) > 0)
        % --- calculate the error
        rms5=myout3(:,3).^2;
        rms5s=smooth(myout3(:,2),rms5,0.2);

        % ---- remove 4 sigma outliers
        ix=find((rms5<4^2*rms5s | rms5<4^2*median(rms5s)));
        myout4=myout3(ix,:);

        RMS2b=std(myout4(:,3));

        rms6=myout4(:,3).^2;
        errsa=smooth(myout4(:,2),rms6,0.1);
        ix=find(errsa<(RMS2b/4)^2);
        errsa(ix)=(RMS2b/4)^2;
        errs=smooth(myout4(:,2),errsa,0.02);  % error^2
        % errs=smooth(myout4(:,2),errs,0.05);  % error^2
        ix=find(errs<(RMS2b/4)^2);
        errs(ix)=(RMS2b/4)^2;

        griderr(1:kx1+kx2+kx3)=99;  % 99 means N/A, extrapolate
        % interpolate grid error
        fi2=[myout4(:,2),errs.^0.5];
        [s,ir]=sort(fi2(:,1));
        fi2=fi2(ir,:);
        xi2=fi2(:,1);  % ref mag
        yi2=fi2(:,2);  % err in mag
        dfi=diff(xi2);
        k2=0;

        for i=1:size(dfi(:),1)  % remove identical iso, interp1: The data abscissae should be distinct
          if dfi(i)>0
            k2=k2+1;
            xi3(k2)=xi2(i);
            yi3(k2)=yi2(i);
          end
        end

        if (k2>10) && (length(errs)>10)
          griderr(kx1+1:kx1+kx2)=interp1(xi3,yi3, maggrid(kx1+1:kx1+kx2),'cubic','extrap');
          griderr2(kx1+1:kx1+kx2)=smooth(maggrid(kx1+1:kx1+kx2),griderr(kx1+1:kx1+kx2),0.15);
        end

        clear fi2 dfi gridiso3 gridmagiso3 xi3 yi3

        % Now find rmsout from griderr
        ixrms = find((magout > refmag4crit) & (magout  < stdpointa));
        rmsout(1:length(magout)) = 99;
        rmsout(ixrms) = interp1(maggrid,griderr,magout(ixrms),'cubic','extrap');

        % ---------------- output files:

        myout=input2;  % [REF; ra; dec; stdmag; color; iso;  NUMBER; BFLAGS; X_IMAGE; Y_IMAGE; ddec; dra; extinction; spatial_bin; local_bin]
        myout(:,16)=magout;  % fitted mag for the given star  myout1=[iso2,magout,dmagout];
        myout(:,17)=dmagout;
        myout(:,18) = rmsout';
        outsize = length(myout(:,1));
        fid = fopen([name solutionString qualifier '_a' num2str(bin) debugSuffix '.out'],'wt');

        for outindex = 1:outsize
          fprintf(fid,"%s %f %f %f %f %f %d %d %f %f %f %f %f %f %f %f %f %f\n",input2REF(outindex,:),myout(outindex,2:18));
        end

        fclose(fid);

        % The upper limit is the maggrid value where maggrid+griderr == stdpointa
        upper = stdpointa - griderr;
        ixupper = find((flaggrid == 1) & ((stdpointa-maggrid) < griderr));

        if (length(ixupper) > 0)
          upper_limit = upper(min(ixupper));
          upper_limit = min(upper_limit,stdpointa - RMS2b);
        else
          upper_limit = stdpointa - RMS2b;
        end

        clear upper ixupper

        % limiting_mag, rms, n1, n2, max_bright_mag, max_bright_iso, limiting_iso, upper_limit, colorterm ,errorcolor, colorflag
        mypara=[stdpointa,RMS2b,length(refmag2),ksize1,refmag4crit,iso4scrit,isopointa,upper_limit,colorterm,errorcolor,colorflag];
        fid = fopen([name solutionString qualifier '_a' num2str(bin) debugSuffix '.para'],'wt');
        fprintf(fid,"%f %f %.0f %.0f %f %f %f %f %f %f %d\n",mypara(:));
        fclose(fid);

        % --------------- output grids
        mygrid=[maggrid', isogrid', griderr', flaggrid'];
        outsize = length(mygrid(:,1));
        fid = fopen([name solutionString qualifier '_a' num2str(bin) debugSuffix '.grid'],'wt');

        for outindex = 1:outsize
          fprintf(fid,"%f %f %f %.0f\n",mygrid(outindex,:));
        end

        fclose(fid);

        % ---------- figure2
        if strcmp(do_plots,"YES")
          if strcmp(version(),"3.8.2")
            graphics_toolkit("gnuplot");
          end

          ix=find(griderr<90);

          if (ix > 2)
            if strcmp(version(),"2.9.8")
              oneplot();
              __gnuplot_set__ term postscript eps color

              command = ['__gnuplot_set__ output  """' name solutionString qualifier '_a'  num2str(bin) '_b' debugSuffix '.eps' '"""'];
              eval(command);

              __gnuplot_set__ pointsize 5
              __gnuplot_set__ key right
              __gnuplot_set__ key top
              plotlabel  = ['GSC2.3.2 B mag  (ISO Smoothing  RMS = ' num2str(RMS2b,3) ')' ];
              xlabel(plotlabel);
              ylabel("Error (mag) or ISO - B mag");
              plottitle = [name solutionString qualifier debugSuffix ' Photometry Calibration, bin ' num2str(bin) '/9'  ];
              title(plottitle);
              axislimit = [min(maggrid),max(maggrid),0,max(isogrid-maggrid-min(isogrid-maggrid))+0.3];
              axis(axislimit);

              plot(maggrid,isogrid-maggrid-min(isogrid-maggrid),'b;ISO-B;',...
                    maggrid(ix),10*griderr(ix),'r;10*error;');

              closeplot;
            else % 3.0.3 plot
              figure(1,'visible','off');
              newplot();
              h2=plot(maggrid,isogrid-maggrid-min(isogrid-maggrid),'b',maggrid(ix),10*griderr(ix),'r','linewidth',1.5);
              ylabel('Error (mag) or ISO - B mag','Fontsize',18)
              xlabel('GSC2.2 B mag','Fontsize',18)
              ylim([0,max(isogrid-maggrid-min(isogrid-maggrid))+0.3])
              title([name solutionString qualifier ' Photometry Calibration, bin ' num2str(bin) '/9'],'Fontsize',18)
              legend('ISO-B','10*error','northeast outside')
              text(min(maggrid)+1,(max(isogrid-maggrid-min(isogrid-maggrid))+0.3)*0.9,['ISO Smoothing  RMS = ' num2str(RMS2b,3) ],'Fontsize',15, 'Color', 'k')
              print([name solutionString qualifier '_a'  num2str(bin) '_b.eps'],'-color')  % print -dpsc2 test
              close

            end % 3.03 plot
          end % length(ix) > 2
        end % do_plot

        % -------------------------- plot figure 1
        ix=find(griderr<90);

        if (ix > 2)
          if strcmp(do_plots,"YES")
            if strcmp(version(),"2.9.8")
              __gnuplot_set__ term postscript eps color

              command = ['__gnuplot_set__ output  """' name solutionString qualifier '_a'  num2str(bin) '_a' debugSuffix '.eps' '"""'];
              eval(command);

              __gnuplot_set__ pointsize 5
              multiplot(1,2);
              plottitle = [name solutionString qualifier debugSuffix ' Photometry Calibration, bin ' num2str(bin) '/9' ];
              title(plottitle);
              ylabel("Instrumental ISO");
              axislimit = [6,17,floor(min(iso2)),ceil(max(iso2))];
              axis(axislimit,"labely");
                    mplot(refmag2,iso2,"b.",
                    maggrid(1:kx1),isogrid(1:kx1),'g',
                    maggrid(kx1+1:kx1+kx2),isogrid(kx1+1:kx1+kx2),'r',...
                    maggrid(kx1+kx2+1:kx1+kx2+kx3),isogrid(kx1+kx2+1:kx1+kx2+kx3),'c',...
                    [1,1]*stdpointa,[min(iso2),max(iso2)],'m-');

              ix=find(griderr<90);
              title("");
              ylabel("Residuals (mag)");
              plotlabel  = [ 'GSC2.3.2 B mag  (ISO Smoothing  RMS = ' num2str(RMS2b,3) ' Limiting mag:  ' num2str(stdpointa,4) ')' ];
              xlabel(plotlabel);
              axislimit = [6,17,-2,2];
              axis(axislimit,"label");
              mplot(refmag2,-refmag2+magout,'b.',...
                    refmag2,randn(size(iso2))*0,'m-',...
                    maggrid(ix),griderr(ix),'r',...
                    maggrid(ix),-griderr(ix),'r');

              closeplot;
            else % 3.0.3 plot
              figure(1,'visible','off');
              newplot();

              subplot(2,1,1)
              %%%subplot('Position',[0.1 0.5 0.8 0.4])
              h1=plot(refmag2,iso2,'b.','markersize',4);
              hold on
              ylabel('Instrumental ISO','Fontsize',18)
              title([name solutionString qualifier ' Photometry Calibration, bin ' num2str(bin) '/9'],'Fontsize',18)
              set(gca,'XTick',7:16)
              set(gca,'ytick',ceil(min(iso2)):2:ceil(max(iso2)))
              xlim([6 17])
              ylim([floor(min(iso2)) ceil(max(iso2))])
              hold on
              plot(maggrid(1:kx1),isogrid(1:kx1),'g',maggrid(kx1+1:kx1+kx2),isogrid(kx1+1:kx1+kx2),'r',...
                  maggrid(kx1+kx2+1:kx1+kx2+kx3),isogrid(kx1+kx2+1:kx1+kx2+kx3),'c','linewidth',1.5)
              plot([1,1]*stdpointa,[min(iso2),max(iso2)],'k--','linewidth',1.5)
              text(6.5,ceil(max(iso2))-1.5,['Limiting mag:  ' num2str(stdpointa,4)],'Fontsize',15, 'Color', 'k')

              subplot(2,1,2)
              %%%subplot('Position',[0.1 0.1 0.8 0.4])
              plot(refmag2,-refmag2+magout,'b.','markersize',5); % ,refmag6,-refmag6+mag7s,'b.',-refmag6+starmag2(ix6)'
              set(gca,'ytick',0.1*(-20:5:+20))
              xlim([6 17])
              ylim([-2,2])
              ylabel('Residuals (mag)','Fontsize',18)
              xlabel('GSC2.2 B mag','Fontsize',18)
              hold on
              plot(refmag2,randn(size(iso2))*0,'k','linewidth',1)
              ix=find(griderr<90);
              plot(maggrid(ix),griderr(ix),'r',maggrid(ix),-griderr(ix),'r','linewidth',1)
              text(6.5,1.5,['ISO Smoothing  RMS = ' num2str(RMS2b,3) ],'Fontsize',15, 'Color', 'k')

              print([name solutionString qualifier '_a'  num2str(bin) '_a.eps'],'-color')  % print -dpsc2 test
              close
            end % 3.03 plot
          end
        else
          printf("length(ix) %d is too small for bin %d, name %s%s%s\n",length(ix),bin,name,solutionString,qualifier);
        end % length(ix) > 2
      else
        printf("ERROR: Failed to find rms for bin %d, name %s%s%s\n",bin,name,solutionString,qualifier);
      end
    else
      printf("ksize1 %d is too small or stdpoint %f is less than refmag4crit %f+1.5 for bin %d, name %s%s%s\n",ksize1,stdpointa,refmag4crit,bin,name,solutionString,qualifier);
    end

    endTime = time() - binStartTime;
    printf("annular9 completed bin %2d, %5d seconds, %5d refmag2 size, %5d ksize1 name %s%s%s\n",bin,endTime,length(refmag2),ksize1,name,solutionString,qualifier);

    command = ['clear output' num2str(bin) 'REF'];
    eval(command);

    clear RMS2b      flaggrid   input3     iso3s      isomax1    isosb      ixmax      kx3        magmin1    myout3     refmag4    rms3s      s          stdpointa
    clear ans        griderr    input4     iso4       isomin1    isosc      ixmin      kxa        magout     myout4     rms0       rms4       size1      stdpointb
    clear            gridmax    ir         iso4s      isoout     ix         k2         maggrid    mygrid     mypara     rms0s      rms4s      size2      xi2
    clear dmagout    gridmin    iso2       iso5s      isopointa  ix1        ksize1     magmax     myout                 rms2       rms5       size3      yi2
    clear errs       i          iso2s      iso6s      isopointb  ix2        kx1        magmax1    myout1     refmag2    rms2s      rms5s      slope1
    clear errsa      input2     iso3       isogrid    isosa      ix3        kx2        magmin     myout2     refmag3    rms3       rms6       slope3
    clear errsb      diso
    clear errsb      diso iso4scrit refmag4crit kx
    clear mref plotlabel  fid mref2   plottitle axislimit   iii  sref
    clear bin command                    outindex         txt   outsize
    clear input2REF
    clear  isoout1 magout1 dmagout1 tmpcrit tmpisoscrit
    clear ixrms rmsout
    clear iso4sort outcmd ky2 refmag4sort
    clear colorterm errorcolor colorflag

    %who -variables
    if (num_bins == "1")
      printf("annular9.m Single bin case\n");
      break;
    end
  end

  % Write a marker file to tell ingest_matlab2.csh that annular9.m completed without crashing
  fid = fopen([binoutput '/' name solutionString qualifier '.annular9_marker.txt'],'w');
  fprintf(fid,'%s%s%s\n',name,solutionString,qualifier);
  fclose(fid);

  endTime = time() - mosaicStartTime;
  printf("annular9.m time is %5d seconds for %s%s%s\n",endTime,name,solutionString,qualifier);
end

endTime = time()-startTime;
printf("annular9.m total time %d seconds\n",endTime);

toc

%  Jun 10, 2007 Sumin Tang    - original version
%  Oct 05, 2007 Edward J. Los - create inputREF to handle multiple non-numeric characters in the REF field.
%  Oct 09, 2007 Edward J. Los - Test size of stdpointa before calculating ksize1
%                               to avoid the following error in find:
%                               error: mx_el_lt: nonconformant arguments (..., op2 is 0x0)
%  Oct 22, 2007 Edward J. Los - force size1 to be at least 2 to guarantee that smooth.m will see a vector
%                               and not come up with a negative span.
%  Nov  4, 2007 Edward J. Los - Test for input number of stars being too small.
%  Nov  6, 2007 Edward J. Los - Remove directory references and use environment variables instead
%  Nov 29, 2007 Sumin Tang    - force the fitting curve to be monotonical in both *.out and *.grid
%                               to fix the non-monotonical problem
%  Dec 11, 2007 Edward J. Los - add NUMBER and FLAGS field. Filter out anything with FLAGS > 128
%  Dec 20, 2007 Edward J. Los - Filter out anything with FLAGS >= 65536
%  Dec 22, 2007 Edward J. Los - Add completion message to help find crashes
%  Dec 27, 2007 Edward J. Los - Avoid divide by zero error (flag it with ksize1 = -4 and -5)
%                               Handle non-existent input files.
%  Dec 29, 2007 Edward J. Los - Support recalibrated GSC magnitudes by adding a qualifier to the file names
%  Jan 18, 2008 Edward J. Los - Make plotting dependent on an environment variable
%                               Make binning dependent on the mosaic pixel location rather than ra and dec.
%                               (This script just passes through X_IMAGE and Y_IMAGE)
%                               Remove filtering with the FLAGS field
%                               Add refmag4crit, iso4scrit, and isopointa to the parameters block
%                               Add rmsout which is griderr interpolated for each output star
%                               Correct refmag4crit array size check
%  Feb  8, 2008 Sumin Tang    - Non-monotonicity revision
%  Feb 12, 2008 Edward J. Los - Always plot figure 1
%                               Add checks for null arrays and infinite loops
%                               Throw out all data less than the minimum iso4s
%                               Correct the bright slope calculation
%  Feb 18, 2008 Edward J. Los - Handle empty input bins
%  Feb 22, 2008 Edward J. Los - Handle single bin case flagged by the DASCH_NUMBINS environment value
%                               Delete any old files
%  Feb 24, 2008 Sumin Tang    - Non-monotonicity revisions
%  Mar  5, 2008 Edward J. Los - Add a column for extinction
%  Mar 28, 2008 Edward J. Los - Remove MAG_APER, ISOAREA, and plate_dist for performance
%  Mar 30, 2008 Edward J. Los - drad should be dra
%  Jul  1, 2008 Edward J. Los - Change GSC2.2 to GSC2.3.2
%  Aug  7, 2008 Sumin Tang    - Revise limiting magnitude calculation.
%  Aug  8, 2008 Edward J. Los - Correct octave crashes
%  Oct 20, 2008 Sumin Tang    - Revise limiting magnitude calculation (Memorandum of Oct 8, 2008)
%  Oct 21, 2008 Edward J. Los - Add upper limit code
%  Oct 22, 2008 Edward J. Los - Correct bug in limiting magnitude calculation.
%  Oct 29, 2008 Edward J. Los - Correct negative index bug in limiting magnitude calculation
%  Nov  3, 2008 Edward J. Los - Correct zero size myout3 when calculating rms
%                             - Add Sumin Tang s upper_limit revision
%  Dec 28, 2008 Edward J. Los - Convert to Octave 3.0.3 plotting
%                               Correct stdpointb calculation if iy1 is zero
%  Jan 12, 2008 Edward J. Los - Correct stdpointc calculation
%  Mar 10, 2009 Edward J. Los - Add experimental smoothed colorterm support
%  Mar 23, 2009 Edward J. Los - If a colorterm is not available for a bin, use one from an innermost bin.  Use flag 3 for this one
%  Apr  6, 2009 Edward J. Los - Add COLOR_RATIO
%  Aug 20, 2009 Edward J. Los - Correct plot output for pass2 cases
%  Sep 15, 2009 Edward J. Los - Test the colorterm marker file and report a crash
%  Sep 29, 2009 Edward J. Los - Add multiple exposure support
%  Mar  5, 2010 Edward J. Los - Add bin-specific performance information
%  Mar 10, 2010 Edward J. Los - Correct an exception when ksize1=-1
%  Mar 15, 2010 Edward J. Los - Write a marker file to detect crashes
%  Aug 27, 2010 Edward J. Los - Update to Octave 3.2.3
%  Oct 19, 2010 Edward J. Los - Conditionally turn off plotting
%  Mar 14, 2010 Edward J. Los - Add spatial_bin and local_bin columns
%  Apr 26, 2011 Edward J. Los - Correct smoothing bug with previously unused iso6s variable
%  Sep  7, 2011 Edward J. Los - Back out deprecated COLOR_RATIO
%  Jan 23, 2012 Edward J. Los - fix short-circuit operators
%                             - avoid using "best" legend position
%  Apr 15, 2014 Edward J. Los - Mark miscellaneous errors deferred with ERRORD
%  Dec  8, 2015 Edward J. Los - Support Octave 3.8.2

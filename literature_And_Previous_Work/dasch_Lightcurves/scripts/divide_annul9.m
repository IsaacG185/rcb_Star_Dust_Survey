% Copyright the President and Fellows of Harvard College.
% Licensed under the MIT License

% divide the stars in the catalog to be 9 groups; also output a whole plate catalog
% matched outputs in /ice5/dasch/photom_wcs/*_tnx/match_*_tnx_u.db:
%
%  Jul 11, 2007  Sumin Tang    - original version
%  Oct 05, 2007  Edward J. Los - create inputREF and output*REF arrays to handle multiple non-numeric characters in the REF field.
%                                The original REF field becomes a dummy zero to avoid changing subscripts.
%  Nov  6, 2007  Edward J. Los - Remove directory references and use environment variables instead
%  Dec 11, 2007  Edward J. Los - Add NUMBER and FLAGS columns
%  Dec 29, 2007  Edward J. Los - Support recalibrated GSC magnitudes by adding a qualifier to the file names
%  Jan 18, 2008  Edward J. Los - Make plotting dependent on an environment variable
%                                Make binning dependent on the mosaic pixel location rather than ra and dec.
%                                Remove filtering by deviation.
%  Feb 15, 2008  Edward J. Los - Restore filtering by deviation
%  Feb 18, 2008  Edward J. Los - Avoid plotting empty bins
%  Feb 22, 2008  Edward J. Los - Handle single bin case flagged by the DASCH_NUMBINS environment value
%                                Delete any old files
%  Feb 27, 2008  Edward J. Los - Handle empty file case
%  Mar  1, 2008  Edward J. Los - Replace sgolay curve fitting with rlowess curve fitting and reduce span by a factor of 100
%  Mar  5, 2008  Edward J. Los - Add a column for extinction
%  Mar  7, 2008  Edward J. Los - Backout the change to rlowess
%                              - Correct exceptions if deviation point rejection returns no points
%  Mar 28, 2008  Edward J. Los - Remove MAG_APER, ISOAREA, and plate_dist for performance
%  Mar 30, 2008  Edward J. Los - drad should really be dra
%                              - read file in native mode
%  Dec 28, 2008  Edward J. Los - Convert plotting to Octave 3.0.3
%  Jan 20, 2009  Edward J. Los - Avoid crash with empty bins which are below the horizon.
%  Feb 27, 2008  Edward J. Los - Delete stale files before proceeding
%  Sep 29, 2009  Edward J. Los - Add multiple exposure support
%  Aug 27, 2010  Edward J. Los - Update to Octave 3.2.3
%  Sep 13, 2010  Edward J. Los - Add a marker file to detect crashes
%  Sep 20, 2010  Edward J. Los - Correct position of marker file write
%  Oct  4, 2010  Edward J. Los - Do not change directories to avoid script load issues
%  Feb 22, 2011  Edward J. Los - Use binarray.m for the calculation of bins.
%  Mar 14, 2011  Edward J. Los - Back out binarray.m, spatial bins and local bins are computed now in prepare_octave.c
%  Apr 13, 2011  Edward J. Los - Write marker files before error exits.
%  Apr 26, 2011  Edward J. Los - Remove initial drad filtering, now covered in other parts of the pipeline
%  Dec  8, 2015  Edward J. Los - Support Octave 3.8.2
%  Oct 19, 2018  Edward J. Los - Support Octave 4.4.1 (use char() fopen())
%
% use 16 columns + the size of the mosaic in pixels
% REF ra_2    dec_2   Stdmag  color   MAG_ISO   NUMBER BFLAGS  X_IMAGE Y_IMAGE ddec dra  extinction spatial_bin local_bin width height scale
%   1   2       3        4      5        6        7      8        9       10    11   12      13         14          15     16     17    18

% should/can also work out a way to devide plates with arbitrary b0/a0 ratio
% eval('cd E:\Plates\GSC2\annular\test');
% name='mc00355_10_01ww';  % mc00355_10_01ww; mc17343_01_01r90ww; mc21438_12_01r90ww

tic
clear
startTime = time();
debugSuffix = "";
qualifier = "";

scripts = getenv("DASCH_SCRIPTS");
printf("scripts is %s\n",scripts);
if (length(scripts) < 1)
  exit;
endif

scratch = getenv("DASCH_SCRATCH");
printf("scratch is %s\n", scratch);
if (length(scratch) < 1)
  exit;
endif

catalogall = getenv("DASCH_CATALOGALL");
printf("catalogall is %s\n",catalogall);
if (length(catalogall) < 1)
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

num_bins = getenv("DASCH_NUMBINS");
printf("num_bins is %s\n",num_bins);
if (length(num_bins) < 1)
  exit;
endif

use_binarray = 1;

addpath(scripts);

chdir(scripts);
pwd();  % this prevents a synchronization problem

if (nargin() < 2)
   printf("Usage: octave divide_annul9.m solutionNumber list [qualifier]\n");
   exit;
endif
solutionNumber = argv()(1);
listfile = argv()(2);
if (nargin() == 3)
    qualifier =  ['_' char(argv()(3))];
else
    qualifier = "";
endif
if (char(solutionNumber) == "0")
    solutionString = "";
else
    solutionString = ['_s' char(solutionNumber)];
endif

printf("List file is %s qualifier is %s use_binarray %d\n",char(listfile),qualifier,use_binarray);
fid = fopen(char(listfile),'rt');
iii = 1;
B = cell(1);
while ((txt = fgetl(fid)) != -1)
   B(iii) = txt;
   iii++;
endwhile

fclose(fid);

ksize=length(B);

for kkk=1:ksize   % 1:ksize  or 249:249
    mosaicStartTime = time();
    name=char(B(kkk));
    printf("Processing %s%s\n",name,solutionString);
    % name2=[name solutionString '_tnx'];
    % input=[REF; ra; dec; stdmag; color; iso;aper;area; platedist; ddec; dra]';

    numpoints = 0;
    printf("Opening %s\n",[catalogall '/' name solutionString qualifier '.db']);
    fid = fopen([catalogall '/' name solutionString qualifier '.db'],'rt');
    if (fid >= 0)

       while ((txt = fgetl(fid)) != -1)
          numpoints++;
       endwhile
       fclose(fid);
    else
       printf("ERROR: divide_annul9 opening %s \n",[catalogall '/' name solutionString qualifier '.db']);
       % write a marker file so the shell does not keep retrying
       filename = [scratch '/' name solutionString qualifier '.divide_annul9.txt'];
       fid = fopen(filename,'w');
       fprintf(fid,'%s%s%s\n',name,solutionString,qualifier);
       fclose(fid);
       continue
    endif
    if (numpoints < 1)
       printf("ERROR: divide_annul9 numpoints is %d for %s%s%s \n",numpoints,name,solutionString,qualifier);
       % write a marker file so the shell does not keep retrying
       filename = [scratch '/' name solutionString qualifier '.divide_annul9.txt'];
       fid = fopen(filename,'w');
       fprintf(fid,'%s%s%s\n',name,solutionString,qualifier);
       fclose(fid);
       continue
    endif

     load([catalogall '/' name solutionString qualifier '.db'],'inputREF','inputx');
     iii = size(inputx,1);


    if iii>0
        width = inputx(1,16);
        height = inputx(1,17);
        scale = inputx(1,18);
        input = inputx(:,1:15);

        clear inputx

        printf("Finished creation of input array\n");

        a0=min(width/2,height/2);
        b0=max(width/2,height/2);

        x_center = 0.5 + (width * 0.5);
        y_center = 0.5 + (height * 0.5);

        printf("Finished read of %s%s%s of size %8d no reduction\n",name,solutionString,qualifier,numpoints);

        save([catalogall '/' name solutionString qualifier '_a0' debugSuffix  '.db'],'inputREF','input')

        x_image=input(:,9);
        y_image=input(:,10);

        if (use_binarray == 1)
           spatial_bin = input(:,14);
        else
           platedist= (((x_image-x_center).^2)+((y_image-y_center).^2)).^0.5;

           min_x_image = 0.5;
           max_x_image = width+0.5;
           min_y_image = 0.5;
           max_y_image = height+0.5;

           % %---- Divide into 9 annulars with similar area
           % a8=max(platedist)*0.8;
           a5=0.95*2*a0*(b0/a0/2/pi).^0.5;
           dr=1.05*2*a0*(4*b0/a0-1)/30/(2*pi*b0/a0)^0.5;
           a6=a5+1.1*dr;
           a7=a5+2.35*dr;
           a8=a5+5*dr;
           a4=a5*(4/5)^(0.5);
           a3=a5*(3/5)^(0.5);
           a2=a5*(2/5)^(0.5);
           a1=a5*(1/5)^(0.5);
        endif

        outcmd = ['rm -f ' catalogbin '/' name solutionString qualifier '_a1' debugSuffix '.db'];
        system(outcmd);
        outcmd = ['rm -f ' catalogbin '/' name solutionString qualifier '_a2' debugSuffix '.db'];
        system(outcmd);
        outcmd = ['rm -f ' catalogbin '/' name solutionString qualifier '_a3' debugSuffix '.db'];
        system(outcmd);
        outcmd = ['rm -f ' catalogbin '/' name solutionString qualifier '_a4' debugSuffix '.db'];
        system(outcmd);
        outcmd = ['rm -f ' catalogbin '/' name solutionString qualifier '_a5' debugSuffix '.db'];
        system(outcmd);
        outcmd = ['rm -f ' catalogbin '/' name solutionString qualifier '_a6' debugSuffix '.db'];
        system(outcmd);
        outcmd = ['rm -f ' catalogbin '/' name solutionString qualifier '_a7' debugSuffix '.db'];
        system(outcmd);
        outcmd = ['rm -f ' catalogbin '/' name solutionString qualifier '_a8' debugSuffix '.db'];
        system(outcmd);
        outcmd = ['rm -f ' catalogbin '/' name solutionString qualifier '_a9' debugSuffix '.db'];
        system(outcmd);
        outcmd = ['rm -f ' catalogbin '/' name solutionString qualifier '_9bins_a' debugSuffix '.eps'];
        system(outcmd);

        if (num_bins == "1")
          printf("divide_annul9.m Single bin case\n");
          output1=input;
          output1REF = inputREF;
          save([catalogbin '/' name solutionString qualifier '_a1' debugSuffix '.db'],'output1REF','output1')  % output catalog for the innermost annular

        else
          if (use_binarray == 1)
             ix1=find(spatial_bin == 1);
             kx1=size(ix1,1);
             if (kx1 > 0)
               output1(1:kx1,:)=input(ix1(1:kx1),:);
               output1REF(1:kx1,:) = inputREF(ix1(1:kx1),:);
               save([catalogbin '/' name solutionString qualifier '_a1' debugSuffix '.db'],'output1REF','output1')  % output catalog for the innermost annular
             else
               output1=input(1,:); % dummy point to avoid plotting error
             endif

             ix2=find(spatial_bin == 2);
             kx2=size(ix2,1);
             if (kx2 > 0)
               output2(1:kx2,:)=input(ix2(1:kx2),:);
               output2REF(1:kx2,:) = inputREF(ix2(1:kx2),:);
               save([catalogbin '/' name solutionString qualifier '_a2' debugSuffix '.db'],'output2REF','output2')
             else
               output2=input(1,:); % dummy point to avoid plotting error
             endif

             ix3=find(spatial_bin == 3);
             kx3=size(ix3,1);
             if (kx3 > 0)
               output3(1:kx3,:)=input(ix3(1:kx3),:);
               output3REF(1:kx3,:) = inputREF(ix3(1:kx3),:);
               save([catalogbin '/' name solutionString qualifier '_a3' debugSuffix '.db'],'output3REF','output3')
             else
               output3=input(1,:); % dummy point to avoid plotting error
             endif

             ix4=find(spatial_bin == 4);
             kx4=size(ix4,1);
             if (kx4 > 0)
               output4(1:kx4,:)=input(ix4(1:kx4),:);
               output4REF(1:kx4,:) = inputREF(ix4(1:kx4),:);
               save([catalogbin '/' name solutionString qualifier '_a4' debugSuffix '.db'],'output4REF','output4')
             else
               output4=input(1,:); % dummy point to avoid plotting error
             endif

             % trim off 10% of the edge
             dedge1=0.1;
             dedge2=0.055;

             ix5=find(spatial_bin == 5);
             kx5=size(ix5,1);
             if (kx5 > 0)
               output5(1:kx5,:)=input(ix5(1:kx5),:);
               output5REF(1:kx5,:) = inputREF(ix5(1:kx5),:);
               save([catalogbin '/' name solutionString qualifier '_a5' debugSuffix '.db'],'output5REF','output5')
             else
               output5=input(1,:); % dummy point to avoid plotting error
             endif

             ix6=find(spatial_bin == 6);
             kx6=size(ix6,1);
             if (kx6 > 0)
               output6(1:kx6,:)=input(ix6(1:kx6),:);
               output6REF(1:kx6,:) = inputREF(ix6(1:kx6),:);
               save([catalogbin '/' name solutionString qualifier '_a6' debugSuffix '.db'],'output6REF','output6')
             else
               output6=input(1,:); % dummy point to avoid plotting error
             endif

             ix7=find(spatial_bin == 7);
             kx7=size(ix7,1);
             if (kx7 > 0)
               output7(1:kx7,:)=input(ix7(1:kx7),:);
               output7REF(1:kx7,:) = inputREF(ix7(1:kx7),:);
               save([catalogbin '/' name solutionString qualifier '_a7' debugSuffix '.db'],'output7REF','output7')
             else
               output7=input(1,:); % dummy point to avoid plotting error
             endif

             ix8=find(spatial_bin == 8);
             kx8=size(ix8,1);
             if (kx8 > 0)
               output8(1:kx8,:)=input(ix8(1:kx8),:);
               output8REF(1:kx8,:) = inputREF(ix8(1:kx8),:);
               save([catalogbin '/' name solutionString qualifier '_a8' debugSuffix '.db'],'output8REF','output8')
             else
               output8=input(1,:); % dummy point to avoid plotting error
             endif

             ix9=find(spatial_bin == 9);
             kx9=size(ix9,1);
             if (kx9 > 0)
               output9(1:kx9,:)=input(ix9(1:kx9),:);
               output9REF(1:kx9,:) = inputREF(ix9(1:kx9),:);
               save([catalogbin '/' name solutionString qualifier '_a9' debugSuffix '.db'],'output9REF','output9')
             else
               output9=input(1,:); % dummy point to avoid plotting error
             endif

          else

             ix1=find(platedist<a1);
             kx1=size(ix1,1);
             if (kx1 > 0)
               output1(1:kx1,:)=input(ix1(1:kx1),:);
               output1(1:kx1,14) = 1;
               output1REF(1:kx1,:) = inputREF(ix1(1:kx1),:);
               save([catalogbin '/' name solutionString qualifier '_a1' debugSuffix '.db'],'output1REF','output1')  % output catalog for the innermost annular
             else
               output1=input(1,:); % dummy point to avoid plotting error
             endif

             ix2=find(platedist<a2&platedist>=a1);
             kx2=size(ix2,1);
             if (kx2 > 0)
               output2(1:kx2,:)=input(ix2(1:kx2),:);
               output2(1:kx2,14) = 2;
               output2REF(1:kx2,:) = inputREF(ix2(1:kx2),:);
               save([catalogbin '/' name solutionString qualifier '_a2' debugSuffix '.db'],'output2REF','output2')
             else
               output2=input(1,:); % dummy point to avoid plotting error
             endif

             ix3=find(platedist<a3&platedist>=a2);
             kx3=size(ix3,1);
             if (kx3 > 0)
               output3(1:kx3,:)=input(ix3(1:kx3),:);
               output3(1:kx3,14) = 3;
               output3REF(1:kx3,:) = inputREF(ix3(1:kx3),:);
               save([catalogbin '/' name solutionString qualifier '_a3' debugSuffix '.db'],'output3REF','output3')
             else
               output3=input(1,:); % dummy point to avoid plotting error
             endif

             ix4=find(platedist<a4&platedist>=a3);
             kx4=size(ix4,1);
             if (kx4 > 0)
               output4(1:kx4,:)=input(ix4(1:kx4),:);
               output4(1:kx4,14) = 4;
               output4REF(1:kx4,:) = inputREF(ix4(1:kx4),:);
               save([catalogbin '/' name solutionString qualifier '_a4' debugSuffix '.db'],'output4REF','output4')
             else
               output4=input(1,:); % dummy point to avoid plotting error
             endif

             % trim off 10% of the edge
             dedge1=0.1;
             dedge2=0.055;

             ix5=find(platedist<a5&platedist>=a4&(x_image<max_x_image-dedge1*a0)&(x_image>min_x_image+dedge1*a0)&(y_image<max_y_image-dedge1*a0)&(y_image>min_y_image+dedge1*a0));
             kx5=size(ix5,1);
             if (kx5 > 0)
               output5(1:kx5,:)=input(ix5(1:kx5),:);
               output5(1:kx5,14) = 5;
               output5REF(1:kx5,:) = inputREF(ix5(1:kx5),:);
               save([catalogbin '/' name solutionString qualifier '_a5' debugSuffix '.db'],'output5REF','output5')
             else
               output5=input(1,:); % dummy point to avoid plotting error
             endif

             ix6=find((platedist<a6&platedist>=a5)&(x_image<max_x_image-dedge1*a0)&(x_image>min_x_image+dedge1*a0)&(y_image<max_y_image-dedge1*a0)&(y_image>min_y_image+dedge1*a0));
             kx6=size(ix6,1);
             if (kx6 > 0)
               output6(1:kx6,:)=input(ix6(1:kx6),:);
               output6(1:kx6,14) = 6;
               output6REF(1:kx6,:) = inputREF(ix6(1:kx6),:);
               save([catalogbin '/' name solutionString qualifier '_a6' debugSuffix '.db'],'output6REF','output6')
             else
               output6=input(1,:); % dummy point to avoid plotting error
             endif

             ix7=find(platedist<a7&platedist>=a6&(x_image<max_x_image-dedge1*a0)&(x_image>min_x_image+dedge1*a0)&(y_image<max_y_image-dedge1*a0)&(y_image>min_y_image+dedge1*a0));
             kx7=size(ix7,1);
             if (kx7 > 0)
               output7(1:kx7,:)=input(ix7(1:kx7),:);
               output7(1:kx7,14) = 7;
               output7REF(1:kx7,:) = inputREF(ix7(1:kx7),:);
               save([catalogbin '/' name solutionString qualifier '_a7' debugSuffix '.db'],'output7REF','output7')
             else
               output7=input(1,:); % dummy point to avoid plotting error
             endif

             ix8=find(platedist<a8 & platedist>=a7 &(x_image<max_x_image-dedge1*a0)&(x_image>min_x_image+dedge1*a0)&(y_image<max_y_image-dedge1*a0)&(y_image>min_y_image+dedge1*a0));
             kx8=size(ix8,1);
             if (kx8 > 0)
               output8(1:kx8,:)=input(ix8(1:kx8),:);
               output8(1:kx8,14) = 8;
               output8REF(1:kx8,:) = inputREF(ix8(1:kx8),:);
               save([catalogbin '/' name solutionString qualifier '_a8' debugSuffix '.db'],'output8REF','output8')
             else
               output8=input(1,:); % dummy point to avoid plotting error
             endif

             ix9=find((platedist>=a8|(x_image>=max_x_image-dedge1*a0)|(x_image<=min_x_image+dedge1*a0)|(y_image>=max_y_image-dedge1*a0)|(y_image<=min_y_image+dedge1*a0)));
             kx9=size(ix9,1);
             if (kx9 > 0)
               output9(1:kx9,:)=input(ix9(1:kx9),:);
               output9(1:kx9,14) = 9;
               output9REF(1:kx9,:) = inputREF(ix9(1:kx9),:);
               save([catalogbin '/' name solutionString qualifier '_a9' debugSuffix '.db'],'output9REF','output9')
             else
               output9=input(1,:); % dummy point to avoid plotting error
             endif
          endif
          if strcmp(do_plots,"YES")
            if strcmp(version(),"2.9.8")
              __gnuplot_set__ term postscript eps color

              command = ['__gnuplot_set__ output  """' catalogbin '/' name solutionString qualifier '_9bins_a' debugSuffix '.eps' '"""'];
              eval(command);
              __gnuplot_set__ pointsize 5
              __gnuplot_set__ key right
              __gnuplot_set__ key top
              __gnuplot_set__ key outside
              __gnuplot_set__ xlabel """RA"""
              __gnuplot_set__ ylabel """DEC"""


              bin1 = ['g.;bin 1 ' num2str(kx1) ';'];
              bin2 = ['m.;bin 2 ' num2str(kx2) ';'];
              bin3 = ['6.;bin 3 ' num2str(kx3) ';'];
              bin4 = ['r.;bin 4 ' num2str(kx4) ';'];
              bin5 = ['g.;bin 5 ' num2str(kx5) ';'];
              bin6 = ['m.;bin 6 ' num2str(kx6) ';'];
              bin7 = ['6.;bin 7 ' num2str(kx7) ';'];
              bin8 = ['c.;bin 8 ' num2str(kx8) ';'];
              bin9 = ['b.;bin 9 ' num2str(kx9) ';'];
              plottitle = [catalogbin '/' name solutionString qualifier ' 9 bins'];
              title(plottitle)
              axis("tight","equal")
              commaflag = 0;
              outcmd = ['plot('];
              if (kx9 > 0)
                 if (commaflag > 0)
                    outcmd = [outcmd ','];
                 endif
                 outcmd = [outcmd 'output9(:,2),output9(:,3),bin9'];
                 commaflag = 1;
              endif
              if (kx8 > 0)
                 if (commaflag > 0)
                    outcmd = [outcmd ','];
                 endif
                 outcmd = [outcmd 'output8(:,2),output8(:,3),bin8'];
                 commaflag = 1;
              endif
              if (kx7 > 0)
                 if (commaflag > 0)
                    outcmd = [outcmd ','];
                 endif
                 outcmd = [outcmd 'output7(:,2),output7(:,3),bin7'];
                 commaflag = 1;
              endif
              if (kx6 > 0)
                 if (commaflag > 0)
                    outcmd = [outcmd ','];
                 endif
                 outcmd = [outcmd 'output6(:,2),output6(:,3),bin6'];
                 commaflag = 1;
              endif
              if (kx5 > 0)
                 if (commaflag > 0)
                    outcmd = [outcmd ','];
                 endif
                 outcmd = [outcmd 'output5(:,2),output5(:,3),bin5'];
                 commaflag = 1;
              endif
              if (kx4 > 0)
                 if (commaflag > 0)
                    outcmd = [outcmd ','];
                 endif
                 outcmd = [outcmd 'output4(:,2),output4(:,3),bin4'];
                 commaflag = 1;
              endif
              if (kx3 > 0)
                 if (commaflag > 0)
                    outcmd = [outcmd ','];
                 endif
                 outcmd = [outcmd 'output3(:,2),output3(:,3),bin3'];
                 commaflag = 1;
              endif
              if (kx2 > 0)
                 if (commaflag > 0)
                    outcmd = [outcmd ','];
                 endif
                 outcmd = [outcmd 'output2(:,2),output2(:,3),bin2'];
                 commaflag = 1;
              endif
              if (kx1 > 0)
                 if (commaflag > 0)
                    outcmd = [outcmd ','];
                 endif
                 outcmd = [outcmd 'output1(:,2),output1(:,3),bin1'];
                 commaflag = 1;
              endif
              outcmd = [outcmd ');'];

              if (commaflag > 0)
                 eval(outcmd);
              end
              closeplot;
          else % 3.0.3 plot
              if strcmp(version(),"3.8.2")
                  graphics_toolkit("gnuplot");
              endif
              figure(1,'visible','off');
              h=plot(output1(:,2),output1(:,3),'g+','markersize',1,output2(:,2),output2(:,3),'m+','markersize',1,output3(:,2),output3(:,3),'y+','markersize',1,output4(:,2),output4(:,3),'r+','markersize',1,...
                  output5(:,2),output5(:,3),'g+','markersize',1,output6(:,2),output6(:,3),'m+','markersize',1,output7(:,2),output7(:,3),'y+','markersize',1,output8(:,2),output8(:,3),'c+','markersize',1,...
                  output9(:,2),output9(:,3),'b+','markersize',1);
              legend(['bin 1, ' num2str(kx1)],['bin 2, ' num2str(kx2)],['bin 3, ' num2str(kx3)],['bin 4, ' num2str(kx4)],...
                  ['bin 5, ' num2str(kx5)],['bin 6, ' num2str(kx6)],['bin 7, ' num2str(kx7)],['bin 8, ' num2str(kx8)],...
                  ['bin 9, ' num2str(kx9)],'Location','NorthEastOutside')
              xlabel('RA','Fontsize',18)
              ylabel('DEC','Fontsize',18)
              title([name solutionString qualifier ' 9 bins'],'Fontsize',18)
              axis equal tight

              print([catalogbin '/' name solutionString qualifier '_9bins_a.eps'],'-color');
              close
            endif
          endif
        end % num_bins

        endTime = time() - mosaicStartTime;
        printf("divide_annul9.m time is %5d seconds for %s%s%s\n",endTime,name,solutionString,qualifier);


        clear outcmd commaflag
        clear a0         a5         b0         dedge2     input      ix3        ix8        kx3        kx8        min_x_image     output2    output7    sdec
        clear a1         a6         ddec       dr         ix         ix4        ix9        kx4        kx9        min_y_image      output3    output8    sra
        clear a2         a7         ddeviate   dra        ix0        ix5        kx0        kx5        max_x_image      name2      output4    output9
        clear a3         a8         x_image        drandom    ix1        ix6        kx1        kx6        max_y_image      output0    output5    platedist
        clear a4         ans        dedge1     h          ix2        ix7        kx2        kx7        meandev    output1    output6    y_image
        clear width height x_center y_center ra dec
        clear bin1 bin2 bin3 bin4 bin5 bin6 bin7 bin8 bin9 count fid
        clear iii numpoints plottitle txt txt1 mosaicStartTime
        clear command  b1
        clear inputREF output0REF output1REF output2REF output3REF output4REF output5REF output6REF output7REF output8REF output9REF
        %who -variables
    end

    % Write a marker file to tell annular9.m that colorterm.m completed without crashing
    filename = [scratch '/' name solutionString qualifier '.divide_annul9.txt'];
    fid = fopen(filename,'w');
    fprintf(fid,'%s%s%s\n',name,solutionString,qualifier);
    fclose(fid);

    clear name
end

endTime = time()-startTime;
printf("divide_annul9.m total time %d seconds\n",endTime);

toc

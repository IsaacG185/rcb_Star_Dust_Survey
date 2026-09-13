% Copyright the President and Fellows of Harvard College.
% Licensed under the MIT License

% Variable search in the Kepler fiedl
% output a list of variable candidates
% study the id table parameters of ASAS variables
% Sumin Tang, March 12,2012
%
% Input files: only ngood>=10 stars are included
%     324 ASAS_id_table_matched_ngood10.db
%     123 ASAS_id_table_matched_ngood10_Vamp04.db  % with V_amp>0.4
%   27979 id_table_keplerfield_2011_12_15b_ngood10.db
%
%  Original name: /dasch3/stang/kepler/kepler_variable_ASAS.m
%  Sep  6, 2012 Convert to Octave
%  Oct 26, 2012 Correct APASS designation
%  Jan 10, 2012 Add nearbyObjects

clear
startTime = time();
%debug_on_warning(1);
%eval('cd /home/scanner/junk/variable_search')

doPlots = 1;
keplerchips = 0;
fname='matched_apass';
fqualifier='_matched_apass'
xmaglabel='APASS V  mag'
%inputname0='id_table_apass2_ngood10.db';
%inputname0='id_table_apass_ngood10_2012_08_29.db';
inputname0='los11.db';
keplerchips = 0;

%fname='NONE';
%fqualifier='_NONE_apass'
%xmaglabel='APASS V mag'
%inputname0='id_table_apass2_ngood10.db';
%inputname0='id_table_apass_ngood10_2012_08_29.db';


%keplerchips = 1
%fname='Kepler';
%fqualifier='_Kepler_ASAS'
%xmaglabel='KIC g mag'
%inputname0='id_table_keplerfield_2011_12_15b_ngood10.db';

keplerchip = 'kepler_ra_dec_4_seasons.db';
%inputname2='ASAS_id_table_matched_ngood10.db';
%inputname3='ASAS_id_table_matched_ngood10_Vamp04.db';
%inputdir='/home/scanner/Pipeline/ingest';
%outputdir='/home/scanner/Pipeline/ingest/plots';
%outputdir2='/home/scanner/Pipeline/ingest/output';

inputdir='/dasch/Pipeline/ingest';
outputdir='/dasch/Pipeline/ingest/plots';
outputdir2='/dasch/Pipeline/ingest/output';

%[REF0 drad0 dradB0 max_drad0 max_dradB0 yrbegin0 yrend0 npoints0 min_local0 max_local0 range_local0 min_local20 max_local20 range_local20 median_local0 rms_local0 ngood0 ngoodB0 clip_median_local0 clip_rms_local0 clip_ngood0 median_iso0 range_iso0 error_bar_factor0 Malmquist_factor0 Malmquist_factorB0 Damon_factor0 Sextractor_Blend0 nblend0 nNonDamonBlue0 medNonDamonBlue0 nonDamonBluerms0 nDamonBlue0 medDamonBlue0 damonBlueRms0 magvslimitingcorr0 magvsracorr0 magvsdeccorr0 rmsdradrms20 rarms0 decrms0 nburst0 nburst20 nburst30 nburst40 ndip0 ndip20 ndip30 ndip40 ndev30 ndev20 adjacentburstdip0 adjacentburstdip20 adjacentburstdip30 lightcurverms10 lightcurverms20 lightcurverms30 lightcurverms40 lightcurverms50 slope_all0 slope_all_err0 nslope_all0 slope_600 slope_60_err0 nslope_600 rms_factor0 dmagcatalog0 Stdmag0 color0 ra0 declination0 MAGFlag0 gscclass0 VFlag0 RaPM0 DecPM0 keplerField0 versionId0]...
%    =textread(inputname0,'%q %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %q %f','headerlines',2);

    numpoints = 0;
    printf("Opening %s\n",[inputdir '/' inputname0]); 
    fid = fopen([inputdir '/' inputname0],'rt');
    if (fid >= 0) 
        
       while ((txt = fgetl(fid)) != -1)
          numpoints++;
       endwhile
       fclose(fid);
    else
       printf("ERROR:  opening %s \n",[inputdir '/' inputname0 '.db']); 
       exit
    endif
    numpoints = numpoints - 2;
    if (numpoints < 1) 
       printf("ERROR:  numpoints is %d \n");
       % write a marker file so the shell does not keep retrying
       exit
    endif
 

    
    REF9 = cell(numpoints,1);
    drad9(numpoints) = 0;
    dradB9(numpoints) = 0;
    max_drad9(numpoints) = 0;
    max_dradB9(numpoints) = 0;
    yrbegin9(numpoints) = 0;
    yrend9(numpoints) = 0;
    npoints9(numpoints) = 0;
    min_local9(numpoints) = 0;
    max_local9(numpoints) = 0;
    range_local9(numpoints) = 0;
    min_local29(numpoints) = 0;
    max_local29(numpoints) = 0;
    range_local29(numpoints) = 0;
    median_local9(numpoints) = 0;
    rms_local9(numpoints) = 0;
    ngood9(numpoints) = 0;
    ngoodB9(numpoints) = 0;
    clip_median_local9(numpoints) = 0;
    clip_rms_local9(numpoints) = 0;
    clip_ngood9(numpoints) = 0;
    median_iso9(numpoints) = 0;
    range_iso9(numpoints) = 0;
    error_bar_factor9(numpoints) = 0;
    Malmquist_factor9(numpoints) = 0;
    Malmquist_factorB9(numpoints) = 0;
    Damon_factor9(numpoints) = 0;
    Sextractor_Blend9(numpoints) = 0;
    nblend9(numpoints) = 0;
    nNonDamonBlue9(numpoints) = 0;
    medNonDamonBlue9(numpoints) = 0;
    nonDamonBluerms9(numpoints) = 0;
    nDamonBlue9(numpoints) = 0;
    medDamonBlue9(numpoints) = 0;
    damonBlueRms9(numpoints) = 0;
    magvslimitingcorr9(numpoints) = 0;
    magvsracorr9(numpoints) = 0;
    magvsdeccorr9(numpoints) = 0;
    rmsdradrms29(numpoints) = 0;
    rarms9(numpoints) = 0;
    decrms9(numpoints) = 0;
    nburst9(numpoints) = 0;
    nburst29(numpoints) = 0;
    nburst39(numpoints) = 0;
    nburst49(numpoints) = 0;
    ndip9(numpoints) = 0;
    ndip29(numpoints) = 0;
    ndip39(numpoints) = 0;
    ndip49(numpoints) = 0;
    ndev39(numpoints) = 0;
    ndev29(numpoints) = 0;
    adjacentburstdip9(numpoints) = 0;
    adjacentburstdip29(numpoints) = 0;
    adjacentburstdip39(numpoints) = 0;
    lightcurverms19(numpoints) = 0;
    lightcurverms29(numpoints) = 0;
    lightcurverms39(numpoints) = 0;
    lightcurverms49(numpoints) = 0;
    lightcurverms59(numpoints) = 0;
    slope_all9(numpoints) = 0;
    slope_all_err9(numpoints) = 0;
    nslope_all9(numpoints) = 0;
    slope_609(numpoints) = 0;
    slope_60_err9(numpoints) = 0;
    nslope_609(numpoints) = 0;
    rms_factor9(numpoints) = 0;
    dmagcatalog9(numpoints) = 0;
    Stdmag9(numpoints) = 0;
    color9(numpoints) = 0;
    ra9(numpoints) = 0;
    declination9(numpoints) = 0;
    MAGFlag9(numpoints) = 0;
    gscclass9(numpoints) = 0;
    VFlag9(numpoints) = 0;
    RaPM9(numpoints) = 0;
    DecPM9(numpoints) = 0;
    keplerField9 = cell(numpoints,1);
    versionId9(numpoints) = 0;
    nearbyObjects9 = cell(numpoints,1);
    lon9(numpoints) = 0;
    lat9(numpoints) = 0;
    
      count = 1;
    iii = 0;
    fid = fopen([inputdir '/' inputname0],'rt');
    txt = fgetl(fid);
    txt = fgetl(fid);
    while ((txt = fgetl(fid)) != -1)
       iii++;
       

[REF9(iii),drad9(iii),dradB9(iii),max_drad9(iii),max_dradB9(iii),yrbegin9(iii),yrend9(iii),npoints9(iii),min_local9(iii),max_local9(iii),range_local9(iii),min_local29(iii),max_local29(iii),range_local29(iii),median_local9(iii),rms_local9(iii),ngood9(iii),ngoodB9(iii),clip_median_local9(iii),clip_rms_local9(iii),clip_ngood9(iii),median_iso9(iii),range_iso9(iii),error_bar_factor9(iii),Malmquist_factor9(iii),Malmquist_factorB9(iii),Damon_factor9(iii),Sextractor_Blend9(iii),nblend9(iii),nNonDamonBlue9(iii),medNonDamonBlue9(iii),nonDamonBluerms9(iii),nDamonBlue9(iii),medDamonBlue9(iii),damonBlueRms9(iii),magvslimitingcorr9(iii),magvsracorr9(iii),magvsdeccorr9(iii),rmsdradrms29(iii),rarms9(iii),decrms9(iii),nburst9(iii),nburst29(iii),nburst39(iii),nburst49(iii),ndip9(iii),ndip29(iii),ndip39(iii),ndip49(iii),ndev39(iii),ndev29(iii),adjacentburstdip9(iii),adjacentburstdip29(iii),adjacentburstdip39(iii),lightcurverms19(iii),lightcurverms29(iii),lightcurverms39(iii),lightcurverms49(iii),lightcurverms59(iii),slope_all9(iii),slope_all_err9(iii),nslope_all9(iii),slope_609(iii),slope_60_err9(iii),nslope_609(iii),rms_factor9(iii),dmagcatalog9(iii),Stdmag9(iii),color9(iii),ra9(iii),declination9(iii),MAGFlag9(iii),gscclass9(iii),VFlag9(iii),RaPM9(iii),DecPM9(iii),keplerField9(iii),versionId9(iii),nearbyObjects9(iii),lon9(iii),lat9(iii),count] = sscanf(txt,"%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%s\t%f\t%s\t%f\t%f","C");
       if (count != 81)
         printf("Line %d count %d\n",iii,count);
          printf("REF9[%d] %s in %s\n",iii,REF9(iii),txt);
         iii--;
       endif
    endwhile
    fclose(fid);
    if (iii != numpoints)
        printf("ERROR: numpoints %d is not iii %d\n",numpoints,iii);
        continue
    endif
    endTime = time() - startTime;
    printf("numpoints %d read in from %s at %d seconds\n",numpoints,[inputdir '/' inputname0],endTime);   
  

%[REF2 drad2 dradB2 max_drad2 max_dradB2 yrbegin2 yrend2 npoints2 min_local2 max_local2 range_local2 min_local22 max_local22 range_local22 median_local2 rms_local2 ngood2 ngoodB2 clip_median_local2 clip_rms_local2 clip_ngood2 median_iso2 range_iso2 error_bar_factor2 Malmquist_factor2 Malmquist_factorB2 Damon_factor2 Sextractor_Blend2 nblend2 nNonDamonBlue2 medNonDamonBlue2 nonDamonBluerms2 nDamonBlue2 medDamonBlue2 damonBlueRms2 magvslimitingcorr2 magvsracorr2 magvsdeccorr2 rmsdradrms22 rarms2 decrms2 nburst2 nburst22 nburst32 nburst42 ndip2 ndip22 ndip32 ndip42 ndev32 ndev22 adjacentburstdip2 adjacentburstdip22 adjacentburstdip32 lightcurverms12 lightcurverms22 lightcurverms32 lightcurverms42 lightcurverms52 slope_all2 slope_all_err2 nslope_all2 slope_602 slope_60_err2 nslope_602 rms_factor2 dmagcatalog2 Stdmag2 color2 ra2 declination2 MAGFlag2 gscclass2 VFlag2 RaPM2 DecPM2 keplerField2 versionId2 No2 ID2 ASASIDV2 ASASIDI2 tMASSID2 RAd2 DEC2 V2 I2 V_I2 J2 J_H2 H_K2 Type2 Period2 V_amp2 I_amp2]...
%    =textread(inputname2,'%q %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %q %f %f %q %q %q %q %f %f %f %f %f %f %f %f %q %f %f %f','headerlines',2);

%[REF3 drad3 dradB3 max_drad3 max_dradB3 yrbegin3 yrend3 npoints3 min_local3 max_local3 range_local3 min_local23 max_local23 range_local23 median_local3 rms_local3 ngood3 ngoodB3 clip_median_local3 clip_rms_local3 clip_ngood3 median_iso3 range_iso3 error_bar_factor3 Malmquist_factor3 Malmquist_factorB3 Damon_factor3 Sextractor_Blend3 nblend3 nNonDamonBlue3 medNonDamonBlue3 nonDamonBluerms3 nDamonBlue3 medDamonBlue3 damonBlueRms3 magvslimitingcorr3 magvsracorr3 magvsdeccorr3 rmsdradrms23 rarms3 decrms3 nburst3 nburst23 nburst33 nburst43 ndip3 ndip23 ndip33 ndip43 ndev33 ndev23 adjacentburstdip3 adjacentburstdip23 adjacentburstdip33 lightcurverms13 lightcurverms23 lightcurverms33 lightcurverms43 lightcurverms53 slope_all3 slope_all_err3 nslope_all3 slope_603 slope_60_err3 nslope_603 rms_factor3 dmagcatalog3 Stdmag3 color3 ra3 declination3 MAGFlag3 gscclass3 VFlag3 RaPM3 DecPM3 keplerField3 versionId3 No3 ID3 ASASIDV3 ASASIDI3 tMASSID3 RAd3 DEC3 V3 I3 V_I3 J3 J_H3 H_K3 Type3 Period3 V_amp3 I_amp3]...
%    =textread(inputname3,'%q %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %q %f %f %q %q %q %q %f %f %f %f %f %f %f %f %q %f %f %f','headerlines',2);

%arrays for galactic bin limits
latMax = [90.000 85.000 85.000 85.000 85.000 85.000 85.000 75.000 75.000 75.000 75.000 75.000 75.000 75.000 75.000 75.000 75.000 75.000 75.000 65.000 65.000 65.000 65.000 65.000 65.000 65.000 65.000 65.000 65.000 65.000 65.000 65.000 65.000 65.000 65.000 65.000 65.000 55.000 55.000 55.000 55.000 55.000 55.000 55.000 55.000 55.000 55.000 55.000 55.000 55.000 55.000 55.000 55.000 55.000 55.000 55.000 55.000 55.000 55.000 55.000 ];
latMin = [85.000 75.000 75.000 75.000 75.000 75.000 75.000 65.000 65.000 65.000 65.000 65.000 65.000 65.000 65.000 65.000 65.000 65.000 65.000 55.000 55.000 55.000 55.000 55.000 55.000 55.000 55.000 55.000 55.000 55.000 55.000 55.000 55.000 55.000 55.000 55.000 55.000 45.000 45.000 45.000 45.000 45.000 45.000 45.000 45.000 45.000 45.000 45.000 45.000 45.000 45.000 45.000 45.000 45.000 45.000 45.000 45.000 45.000 45.000 45.000 ];
lonMax = [360.000 60.000 120.000 180.000 240.000 300.000 360.000 30.000 60.000 90.000 120.000 150.000 180.000 210.000 240.000 270.000 300.000 330.000 360.000 20.000 40.000 60.000 80.000 100.000 120.000 140.000 160.000 180.000 200.000 220.000 240.000 260.000 280.000 300.000 320.000 340.000 360.000 15.652 31.304 46.957 62.609 78.261 93.913 109.565 125.217 140.870 156.522 172.174 187.826 203.478 219.130 234.783 250.435 266.087 281.739 297.391 313.043 328.696 344.348 360.000 ];
lonMin = [0.000 0.000 60.000 120.000 180.000 240.000 300.000 0.000 30.000 60.000 90.000 120.000 150.000 180.000 210.000 240.000 270.000 300.000 330.000 0.000 20.000 40.000 60.000 80.000 100.000 120.000 140.000 160.000 180.000 200.000 220.000 240.000 260.000 280.000 300.000 320.000 340.000 0.000 15.652 31.304 46.957 62.609 78.261 93.913 109.565 125.217 140.870 156.522 172.174 187.826 203.478 219.130 234.783 250.435 266.087 281.739 297.391 313.043 328.696 344.348 ];



%for iteration=10:85
for iteration = 1:69
    %who 

    if (iteration == 1)
       %This iteration  uses everything
       ixfilter = find(ra9 > -500.);
       prefix = "all_";
    endif
    if ((iteration >= 2) && (iteration <= 9))
       angle = 105.0-(10*iteration);
       ixfilter = find(lat9 >= angle);
       prefix = ['lat_' num2str(angle) "_"];
    endif
    if (iteration > 9)
       topangle = latMax(iteration-9);
       bottomangle = latMin(iteration-9);
       rightangle = lonMax(iteration-9);
       leftangle = lonMin(iteration-9);
       ixfilter = find((lat9 >= bottomangle) & (lat9 < topangle) & (lon9 >= leftangle) & (lon9 <= rightangle));
       printf("Galactic bin %d top %f bottom %f left %f right %f\n",iteration-9,topangle,bottomangle,leftangle,rightangle);
       prefix = ['bin_' num2str(iteration-9) "_"];
 
    endif


    endTime = time() - startTime;
    printf("iteration %d prefix %s filter length %d at %d seconds\n",iteration,prefix,length(ixfilter),endTime);
    if (length(ixfilter) <= 1) 
      continue;
    endif


    REF0 = REF9(ixfilter);
    nearbyObjects0 = nearbyObjects9(ixfilter);
    drad0 = drad9(ixfilter);
    dradB0 = dradB9(ixfilter);
    max_drad0 = max_drad9(ixfilter);
    max_dradB0 = max_dradB9(ixfilter);
    yrbegin0 = yrbegin9(ixfilter);
    yrend0 = yrend9(ixfilter);
    npoints0 = npoints9(ixfilter);
    min_local0 = min_local9(ixfilter);
    max_local0 = max_local9(ixfilter);
    range_local0 = range_local9(ixfilter);
    min_local20 = min_local29(ixfilter);
    max_local20 = max_local29(ixfilter);
    range_local20 = range_local29(ixfilter);
    median_local0 = median_local9(ixfilter);
    rms_local0 = rms_local9(ixfilter);
    ngood0 = ngood9(ixfilter);
    ngoodB0 = ngoodB9(ixfilter);
    clip_median_local0 = clip_median_local9(ixfilter);
    clip_rms_local0 = clip_rms_local9(ixfilter);
    clip_ngood0 = clip_ngood9(ixfilter);
    median_iso0 = median_iso9(ixfilter);
    range_iso0 = range_iso9(ixfilter);
    error_bar_factor0 = error_bar_factor9(ixfilter);
    Malmquist_factor0 = Malmquist_factor9(ixfilter);
    Malmquist_factorB0 = Malmquist_factorB9(ixfilter);
    Damon_factor0 = Damon_factor9(ixfilter);
    Sextractor_Blend0 = Sextractor_Blend9(ixfilter);
    nblend0 = nblend9(ixfilter);
    nNonDamonBlue0 = nNonDamonBlue9(ixfilter);
    medNonDamonBlue0 = medNonDamonBlue9(ixfilter);
    nonDamonBluerms0 = nonDamonBluerms9(ixfilter);
    nDamonBlue0 = nDamonBlue9(ixfilter);
    medDamonBlue0 = medDamonBlue9(ixfilter);
    damonBlueRms0 = damonBlueRms9(ixfilter);
    magvslimitingcorr0 = magvslimitingcorr9(ixfilter);
    magvsracorr0 = magvsracorr9(ixfilter);
    magvsdeccorr0 = magvsdeccorr9(ixfilter);
    rmsdradrms20 = rmsdradrms29(ixfilter);
    rarms0 = rarms9(ixfilter);
    decrms0 = decrms9(ixfilter);
    nburst0 = nburst9(ixfilter);
    nburst20 = nburst29(ixfilter);
    nburst30 = nburst39(ixfilter);
    nburst40 = nburst49(ixfilter);
    ndip0 = ndip9(ixfilter);
    ndip20 = ndip29(ixfilter);
    ndip30 = ndip39(ixfilter);
    ndip40 = ndip49(ixfilter);
    ndev30 = ndev39(ixfilter);
    ndev20 = ndev29(ixfilter);
    adjacentburstdip0 = adjacentburstdip9(ixfilter);
    adjacentburstdip20 = adjacentburstdip29(ixfilter);
    adjacentburstdip30 = adjacentburstdip39(ixfilter);
    lightcurverms10 = lightcurverms19(ixfilter);
    lightcurverms20 = lightcurverms29(ixfilter);
    lightcurverms30 = lightcurverms39(ixfilter);
    lightcurverms40 = lightcurverms49(ixfilter);
    lightcurverms50 = lightcurverms59(ixfilter);
    slope_all0 = slope_all9(ixfilter);
    slope_all_err0 = slope_all_err9(ixfilter);
    nslope_all0 = nslope_all9(ixfilter);
    slope_600 = slope_609(ixfilter);
    slope_60_err0 = slope_60_err9(ixfilter);
    nslope_600 = nslope_609(ixfilter);
    rms_factor0 = rms_factor9(ixfilter);
    dmagcatalog0 = dmagcatalog9(ixfilter);
    Stdmag0 = Stdmag9(ixfilter);
    color0 = color9(ixfilter);
    ra0 = ra9(ixfilter);
    declination0 = declination9(ixfilter);
    MAGFlag0 = MAGFlag9(ixfilter);
    gscclass0 = gscclass9(ixfilter);
    VFlag0 = VFlag9(ixfilter);
    RaPM0 = RaPM9(ixfilter);
    DecPM0 = DecPM9(ixfilter);
    keplerField0 = keplerField9(ixfilter);
    versionId0 = versionId9(ixfilter);
    nearbyObjects0 = nearbyObjects9(ixfilter);
    lon0 = lon9(ixfilter);
    lat0 = lat9(ixfilter);

if (keplerchips == 0) 
    RACCD(1) = min(ra0);
    RACCD(2) = min(ra0)-0.1;
    RACCD(3) = max(ra0)+0.1;
    RACCD(4) = max(ra0);
    DECCCD(1) = min(declination0);
    DECCCD(2) = max(declination0);
    DECCCD(3) = min(declination0)+0.1;
    DECCCD(4) = max(declination0)+0.1;
else
    %[mod out row column RAa DECa RAb DECb RAc DECc RAd DECd RAe DECe]...
    %    =textread('keplerchip','%f %f %f %f %f %f %f %f %f %f %f %f %f %f','headerlines',5);
    %i=1:5:416;
    %RACCD=RAa(i);
    %DECCCD=DECa(i);
    numkeplerchip = 0;
    printf("Opening %s\n",[inputdir '/' keplerchip]); 
    fid = fopen([inputdir '/' keplerchip],'rt');
    if (fid >= 0) 
        
       while ((txt = fgetl(fid)) != -1)
          numkeplerchip++;
       endwhile
       fclose(fid);
    else
       printf("ERROR:  opening %s \n",[inputdir '/' keplerchip '.db']); 
       exit
    endif
    numkeplerchip = numkeplerchip - 5;
    if (numkeplerchip < 1) 
       printf("ERROR:  numkeplerchip is %d \n");
       % write a marker file so the shell does not keep retrying
       exit
    endif
 

    
    mod(numkeplerchip) = 0;
    out(numkeplerchip) = 0;
    row(numkeplerchip) = 0;
    column(numkeplerchip) = 0;
    RAa(numkeplerchip) = 0;
    DECa(numkeplerchip) = 0;
    RAb(numkeplerchip) = 0;
    DECb(numkeplerchip) = 0;
    RAc(numkeplerchip) = 0;
    DECc(numkeplerchip) = 0;
    RAd(numkeplerchip) = 0;
    DECd(numkeplerchip) = 0;
    RAe(numkeplerchip) = 0;
    DECe(numkeplerchip) = 0;
    
      count = 1;
    iii = 0;
    fid = fopen([inputdir '/' keplerchip],'rt');
    txt = fgetl(fid);
    txt = fgetl(fid);
    txt = fgetl(fid);
    txt = fgetl(fid);
    txt = fgetl(fid);
    while ((txt = fgetl(fid)) != -1)
       iii++;
       

[mod(iii),out(iii),row(iii),column(iii),RAa(iii),DECa(iii),RAb(iii),DECb(iii),RAc(iii),DECc(iii),RAd(iii),DECd(iii),RAe(iii),DECe(iii),count] = sscanf(txt,"%f %f %f %f %f %f %f %f %f %f %f %f %f %f","C");
       if (count != 14)
         printf("Line %d count %d\n",iii,count);
          printf("REF9[%d] %s in %s\n",iii,REF9(iii),txt);
         iii--;
       endif

    endwhile
    fclose(fid);
    if (iii != numkeplerchip)
        printf("ERROR: numkeplerchip %d is not iii %d\n",numkeplerchip,iii);
        continue
    endif
    endTime = time() - startTime;
    printf("numkeplerchip %d read in from %s at %d seconds\n",numkeplerchip,[inputdir '/' keplerchip],endTime);   
    i=1:5:416;
    RACCD=RAa(i);
    DECCCD=DECa(i);


endif
eval(['cd ' outputdir])


    figure(1,'visible','off');

% ------Plot the distribution of bins
h=axes('Fontsize',16);
plot(ra0, declination0,'.','markersize',2)
hold on
%plot(ra3,declination3,'ro','markersize',5)
xlabel('RA','fontsize',18)
ylabel('Dec','fontsize',18)
daspect([1 1 1]) % equal scaling along each axis
%legend('all DASCH stars','ASAS variables')

% --- plot the CCD chips
kccd=length(RACCD)/4;
RAchipa=zeros(kccd,1);
RAchipb=zeros(kccd,1);
RAchipc=zeros(kccd,1);
RAchipd=zeros(kccd,1);
DECchipa=zeros(kccd,1);
DECchipb=zeros(kccd,1);
DECchipc=zeros(kccd,1);
DECchipd=zeros(kccd,1);
for i=1:kccd
    RAchip=RACCD(i*4-3:i*4);
    DECchip=DECCCD(i*4-3:i*4);
    ixa=find(RAchip==min(RAchip));
    ixb=find(DECchip==max(DECchip));
    ixc=find(RAchip==max(RAchip));
    ixd=find(DECchip==min(DECchip));
    
    if (keplerchips == 0)
      RAchipa(i)=RAchip(ixa(1))-0.1;
      DECchipa(i)=DECchip(ixa(1));
      RAchipb(i)=RAchip(ixb(1));
      DECchipb(i)=DECchip(ixb(1))+0.1;
      RAchipc(i)=RAchip(ixc(1))+0.1;
      DECchipc(i)=DECchip(ixc(1));
      RAchipd(i)=RAchip(ixd(1));
      DECchipd(i)=DECchip(ixd(1))-0.1;
    else
      RAchipa(i)=RAchip(ixa)-0.1;
      DECchipa(i)=DECchip(ixa);
      RAchipb(i)=RAchip(ixb);
      DECchipb(i)=DECchip(ixb)+0.1;
      RAchipc(i)=RAchip(ixc)+0.1;
      DECchipc(i)=DECchip(ixc);
      RAchipd(i)=RAchip(ixd);
      DECchipd(i)=DECchip(ixd)-0.1;

    endif
    
    plot([RAchipa(i) RAchipb(i)],[DECchipa(i) DECchipb(i)],'g');
    plot([RAchipb(i) RAchipc(i)],[DECchipb(i) DECchipc(i)],'g');
    plot([RAchipc(i) RAchipd(i)],[DECchipc(i) DECchipd(i)],'g');
    plot([RAchipa(i) RAchipd(i)],[DECchipa(i) DECchipd(i)],'g');
    hold on
    text(RAchipc(i)+0.3,DECchipc(i)+0.1,num2str(i),'color','k');
end

set(gca,'XDir','reverse')
title([num2str(length(ra0)) ' stars in the ' fname ' field with ngood>=10'],'fontsize',18)
saveas(gcf,[prefix 'ra_dec_grid_' fname '.eps'],'epsc2')
close

daspect('auto') % sets the data aspect ratio mode to auto.
close

% DASCH stars on the chip which are not likely fake variables
% criteria, not 3-sigma outliers of Malmquist_factor0,
% Malmquist_factorB0, magvslimitingcorr0, rarms0,
% decrms0, magvsracorr0, magvsdeccorr0
xp=Malmquist_factor0;
ixf=find(xp<20);  % remove 99
if (length(ixf) > 0) 
  medxp=median(xp(ixf));
  sigxp=std(xp(ixf));
  ixf=find(abs(xp-medxp)<4*sigxp);  % 4-sigma clipping
  if (length(ixf) > 0)
    medmf=median(xp(ixf));
    sigmf=std(xp(ixf));
  else 
     printf("WARNING: No Malmquist_factor measurements available\n");
     medmf=0;
     sigmf=0;
  endif
else
  printf("WARNING: No Malmquist_factor measurements available\n");
  medmf=0;
  sigmf=0;
endif
xp=Malmquist_factorB0;
if (length(ixf) > 0) 
  ixf=find(xp<20);
  medxp=median(xp(ixf));
  sigxp=std(xp(ixf));
  ixf=find(abs(xp-medxp)<4*sigxp);
  medmfB=median(xp(ixf));
  sigmfB=std(xp(ixf));
else
  printf("WARNING: No Malmquist_factorB measurements available\n");
  medmfB=0;
  sigmfB=0;
endif        

xp=magvslimitingcorr0;
ixf=find(xp<20);  % remove 99
if (length(ixf) > 0) 
  medxp=median(xp(ixf));
  sigxp=std(xp(ixf));
  ixf=find(abs(xp-medxp)<4*sigxp);  % 4-sigma clipping
  medmaglim=median(xp(ixf));
  sigmaglim=std(xp(ixf));
else
  printf("WARNING: No magvslimitingcorr0 measurements available\n");
  medmaglim=0;
  sigmaglim=0;  
endif        

xp=rarms0;
ixf=find(xp<20);  % remove 99
if (length(ixf) > 0) 
  medxp=median(xp(ixf));
  sigxp=std(xp(ixf));
  ixf=find(abs(xp-medxp)<4*sigxp);  % 4-sigma clipping
  medrarms=median(xp(ixf));
  sigrarms=std(xp(ixf));
else
  printf("WARNING: No rarms measurements available\n");
  medrarms=0;
  sigrarms=0;
endif        


xp=decrms0;
if (length(ixf) > 0) 
  ixf=find(xp<20);  % remove 99
  medxp=median(xp(ixf));
  sigxp=std(xp(ixf));
  ixf=find(abs(xp-medxp)<4*sigxp);  % 4-sigma clipping
  meddecrms=median(xp(ixf));
  sigdecrms=std(xp(ixf));
else
  printf("WARNING: No decrarms measurements available\n");
  meddecrms=0;
  sigdecrms=0;
endif        



xp=magvsracorr0;
ixf=find(xp<20);  % remove 99
if (length(ixf) > 0) 
  medxp=median(xp(ixf));
  sigxp=std(xp(ixf));
  ixf=find(abs(xp-medxp)<4*sigxp);  % 4-sigma clipping
  medmagra=median(xp(ixf));
  sigmagra=std(xp(ixf));
else
  printf("WARNING: No magvsracorr0 measurements available\n");
  medmagra=0;
  sigmagra=0;
endif        

xp=magvsdeccorr0;
ixf=find(xp<20);  % remove 99
if (length(ixf) > 0) 
  medxp=median(xp(ixf));
  sigxp=std(xp(ixf));
  ixf=find(abs(xp-medxp)<4*sigxp);  % 4-sigma clipping
  medmagdec=median(xp(ixf));
  sigmagdec=std(xp(ixf));
else
  printf("WARNING: No magvsdeccorr0 measurements available\n");
  medmagdec=0;
  sigmagdec=0;
endif        


for ii=1:kccd  % compare ID table parameters on each CCD chip
    if (keplerchips == 0)
       CCDstring = [prefix ];
    else
       CCDstring = [prefix 'CCD_' num2str(ii)];
    endif
    x1=RAchipa(ii);
    y1=DECchipa(ii);
    x2=RAchipb(ii);
    y2=DECchipb(ii);
    k1=(y2-y1)/(x2-x1);
    b1=(x2*y1-x1*y2)/(x2-x1);
    
    x1=RAchipc(ii);
    y1=DECchipc(ii);
    x2=RAchipb(ii);
    y2=DECchipb(ii);
    k2=(y2-y1)/(x2-x1);
    b2=(x2*y1-x1*y2)/(x2-x1);
    
    x1=RAchipc(ii);
    y1=DECchipc(ii);
    x2=RAchipd(ii);
    y2=DECchipd(ii);
    k3=(y2-y1)/(x2-x1);
    b3=(x2*y1-x1*y2)/(x2-x1);
    
    x1=RAchipa(ii);
    y1=DECchipa(ii);
    x2=RAchipd(ii);
    y2=DECchipd(ii);
    k4=(y2-y1)/(x2-x1);
    b4=(x2*y1-x1*y2)/(x2-x1);
    
    % DASCH stars on the chip which are not likely fake variables
    % criteria, not 3-sigma outliers of Malmquist_factor0,
    % Malmquist_factorB0, magvslimitingcorr0, rarms0,
    % decrms0, magvsracorr0, magvsdeccorr0
    % some real variables may also have this effect (only detected on deeper plates, with a positive mf factor and magvslimcorr)
    ix=find(declination0<k1*ra0+b1 & declination0<k2*ra0+b2 & declination0>k3*ra0+b3 & declination0>k4*ra0+b4...
        & (rarms0-medrarms<4*sigrarms | rarms0==99) & ( decrms0-meddecrms<4*sigdecrms | decrms0==99)...
        & (abs(magvsracorr0-medmagra)<4*sigmagra | magvsracorr0==99) & (abs(magvsdeccorr0-medmagdec)<4*sigmagdec | magvsdeccorr0==99)... 
        & (Malmquist_factor0-medmf>-4*sigmf | Malmquist_factor0==99 )...
       & ( Malmquist_factorB0-medmfB>-4*sigmfB | Malmquist_factorB0==99)...
       & ( magvslimitingcorr0-medmaglim>-4*sigmaglim |magvslimitingcorr0==99)); ...
    printf("Filter outliers before %d after %d\n",length(declination0),length(ix));
    REF=REF0(ix);
    nearby=nearbyObjects0(ix);
    ra=ra0(ix);
    dec=declination0(ix);
    stdmag=Stdmag0(ix);
    ngood=ngood0(ix);
    ngoodB=ngoodB0(ix);
    % ID table parameters
    lcrms1=lightcurverms10(ix);
    lcrms2a=lightcurverms20(ix);
    lcrms3a=lightcurverms30(ix);
    lcrms4a=lightcurverms40(ix);
    lcrms5a=lightcurverms50(ix);
    % if lcrms2-5<0.05, set them to be 0.05
    lcrms2=max(0.05,lcrms2a);
    lcrms3=max(0.05,lcrms3a);
    lcrms4=max(0.05,lcrms4a);
    lcrms5=max(0.05,lcrms5a);
    
    rmslocal=rms_local0(ix);
    
    range1=range_local0(ix);
    range2=range_local20(ix);
    nadj1=adjacentburstdip0(ix);
    nadj2=adjacentburstdip20(ix);
    nadj3=adjacentburstdip30(ix);
    nburst1=nburst0(ix);
    nburst2=nburst20(ix);
    nburst3=nburst30(ix);
    nburst4=nburst40(ix);
    ndip1=ndip0(ix);
    ndip2=ndip20(ix);
    ndip3=ndip30(ix);
    ndip4=ndip40(ix);
    ndev2=ndev20(ix);
    ndev3=ndev30(ix);
    
    % ASAS variables with Vamp>0.4 mag on the chip
%    ix=find(declination3<k1*ra3+b1 & declination3<k2*ra3+b2 & declination3>k3*ra3+b3 & declination3>k4*ra3+b4);
%    REFv=REF3(ix);
%    rav=ra3(ix);
%    decv=declination3(ix);
%    stdmagv=Stdmag3(ix);
%    ngoodv=ngood3(ix);
%    ngoodBv=ngoodB3(ix);
%    typev=Type3(ix);
%    periodv=Period3(ix);
%    vampv=V_amp3(ix);
%    iampv=I_amp3(ix);
%    % ID table parameters
%    lcrms1v=lightcurverms13(ix);
%    lcrms2va=lightcurverms23(ix);
%    lcrms3va=lightcurverms33(ix);
%    lcrms4va=lightcurverms43(ix);
%    lcrms5va=lightcurverms53(ix);
%    % if lcrms2-5<0.05, set them to be 0.05
%    lcrms2v=max(0.05,lcrms2va);
%    lcrms3v=max(0.05,lcrms3va);
%    lcrms4v=max(0.05,lcrms4va);
%    lcrms5v=max(0.05,lcrms5va);
%    
%    rmslocalv=rms_local3(ix);
%    
%    range1v=range_local3(ix);
%    range2v=range_local23(ix);
%    nadj1v=adjacentburstdip3(ix);
%    nadj2v=adjacentburstdip23(ix);
%    nadj3v=adjacentburstdip33(ix);
%    nburst1v=nburst3(ix);
%    nburst2v=nburst23(ix);
%    nburst3v=nburst33(ix);
%    nburst4v=nburst43(ix);
%    ndip1v=ndip3(ix);
%    ndip2v=ndip23(ix);
%    ndip3v=ndip33(ix);
%    ndip4v=ndip43(ix);
%    ndev2v=ndev23(ix);
%    ndev3v=ndev33(ix);
    
    kcandidate=0;
    clear REFcandidate;
    clear REFnearby;
    clear REFreason;
    
    eval(['cd ' outputdir])
    % ------------------------------------------
    % --- Part 1. lightcurve rms and amplitude
    % ------------------------------------------
    % --- 1.0 rmslocal
    dN=8;  % number of bins in ngood and gmag % -------- need to be revised
    myvar=rmslocal;                             % --------
%    myvarv=rmslocalv;                           % --------
    myvarname='rmslocal';                       % --------
    nsig=5; % nsig-sigma outliers             % --------
    rmsthr=0.25; % require >=0.25 mag
    
    xvalue=ngood;
%   xvaluev=ngoodv;
    subplot(2,1,1)
    semilogx(xvalue,myvar,'.')
    hold on
%   semilogx(xvaluev,myvarv,'ro','markersize',8)
    ixvalue=prctile(xvalue',0:100/dN:100)';  % Y = prctile(X,p) calculates a value that is greater than p percent of the values in X. The values of p must lie in the interval [0 100].
    clear medvalue stdvalue stdvalue1 stdvalue2 kxoutlier
    for i=1:length(ixvalue)-1
      if(length(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1))) > 1) 
        medvalue(i)=median(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1)));
        stdvalue(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1)));
        stdvalue1(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue(i))); % 6-sigma clipping
        stdvalue2(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue1(i))); % 6-sigma clipping
        % stdvalue1(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue2(i))); % 6-sigma clipping
        % stdvalue2(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue1(i))); % 6-sigma clipping
        
        ixoutlier=find( xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & myvar>medvalue(i)+nsig*stdvalue2(i) & myvar>=rmsthr);  % nsig-sigma outlier
        semilogx([ixvalue(i) ixvalue(i+1)],[1 1]*medvalue(i)+nsig*stdvalue2(i),'g--')
        kxoutlier(i)=length(ixoutlier);
        if kxoutlier(i)>0
            REFcandidate(kcandidate+1:kcandidate+kxoutlier(i))=REF(ixoutlier);
            REFnearby(kcandidate+1:kcandidate+kxoutlier(i))=nearby(ixoutlier);
            tmpreason = cell(length(ixoutlier),1);
            tmpreason(:) = myvarname;
            REFreason(kcandidate+1:kcandidate+kxoutlier(i))=tmpreason;
            kcandidate=kcandidate+kxoutlier(i);
            semilogx(xvalue(ixoutlier),myvar(ixoutlier),'gd','markersize',12);
        end
      else
         medvalue(i) = 0;
      end
    end

    set(gca,'linewidth',1,'fontsize',12)
    xlabel('ngood','fontsize',16)
    ylabel(myvarname,'fontsize',16)
    title([CCDstring ': ' num2str(nsig) '\sigma outliers of ' myvarname],'fontsize',16)
    
    
    iz=find(stdmag<20);     xvalue=stdmag(iz);    myvar=myvar(iz);  REFmag=REF(iz); nearbymag = nearby(iz);
%   xvaluev=stdmagv;
    subplot(2,1,2)
    plot(xvalue,myvar,'.')
    hold on
%   plot(xvaluev,myvarv,'ro','markersize',8)
    xlim([6 15])
    ixvalue=prctile(xvalue',0:100/dN:100)';
    clear medvalue stdvalue stdvalue1 stdvalue2 kxoutlier
    for i=1:length(ixvalue)-1
      if(length(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1))) > 1) 
        medvalue(i)=median(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1)));
        stdvalue(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1)));
        stdvalue1(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue(i))); % 6-sigma clipping
        stdvalue2(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue1(i))); % 6-sigma clipping
        % stdvalue1(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue2(i))); % 6-sigma clipping
        % stdvalue2(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue1(i))); % 6-sigma clipping
        
        ixoutlier=find( xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & myvar>medvalue(i)+nsig*stdvalue2(i) & myvar>=rmsthr );  % nsig-sigma outlier
        plot([ixvalue(i) ixvalue(i+1)],[1 1]*medvalue(i)+nsig*stdvalue2(i),'g--')
        kxoutlier(i)=length(ixoutlier);
        if kxoutlier(i)>0
            REFcandidate(kcandidate+1:kcandidate+kxoutlier(i))=REFmag(ixoutlier);
            REFnearby(kcandidate+1:kcandidate+kxoutlier(i))=nearbymag(ixoutlier);
            tmpreason = cell(length(ixoutlier),1);
            tmpreason(:) = myvarname;
            REFreason(kcandidate+1:kcandidate+kxoutlier(i))=tmpreason;
            kcandidate=kcandidate+kxoutlier(i);
            plot(xvalue(ixoutlier),myvar(ixoutlier),'gd','markersize',12);
        end
      else
         medvalue(i) = 0;
      end
    end
    
    set(gca,'linewidth',1,'fontsize',12)
    xlabel(xmaglabel,'fontsize',16)
    ylabel(myvarname,'fontsize',16)
    if (doPlots == 1)
      saveas(gcf,[CCDstring '_' myvarname fqualifier '.eps'],'epsc2')
    endif    
    close
    %if sum(kxoutlier)>0
    %    myvarname
    %    REFcandidate
    %end
    
    % --- 1.1 lcrms1
    dN=8;  % number of bins in ngood and gmag % -------- need to be revised
    myvar=lcrms1;                             % --------
%   myvarv=lcrms1v;                           % --------
    myvarname='lcrms1';                       % --------
    nsig=5; % nsig-sigma outliers             % --------
    rmsthr=0.25; % require >=0.25 mag
    
    xvalue=ngood;
%   xvaluev=ngoodv;
    subplot(2,1,1)
    semilogx(xvalue,myvar,'.')
    hold on
%   semilogx(xvaluev,myvarv,'ro','markersize',8)
    
    ixvalue=prctile(xvalue',0:100/dN:100)';  % Y = prctile(X,p) calculates a value that is greater than p percent of the values in X. The values of p must lie in the interval [0 100].
    clear medvalue stdvalue stdvalue1 stdvalue2 kxoutlier
    for i=1:length(ixvalue)-1
      if(length(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1))) > 1) 
        medvalue(i)=median(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1)));
        medvalue(i)=median(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1)));
        stdvalue(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1)));
        if (length(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue(i))) > 1)
        stdvalue1(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue(i))); % 6-sigma clipping
        stdvalue2(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue1(i))); % 6-sigma clipping
        % stdvalue1(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue2(i))); % 6-sigma clipping
        % stdvalue2(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue1(i))); % 6-sigma clipping
        
        ixoutlier=find( xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & myvar>medvalue(i)+nsig*stdvalue2(i) & myvar>=rmsthr);  % nsig-sigma outlier
        semilogx([ixvalue(i) ixvalue(i+1)],[1 1]*medvalue(i)+nsig*stdvalue2(i),'g--')
        kxoutlier(i)=length(ixoutlier);
        if kxoutlier(i)>0
            REFcandidate(kcandidate+1:kcandidate+kxoutlier(i))=REF(ixoutlier);
            REFnearby(kcandidate+1:kcandidate+kxoutlier(i))=nearby(ixoutlier);
            tmpreason = cell(length(ixoutlier),1);
            tmpreason(:) = myvarname;
            REFreason(kcandidate+1:kcandidate+kxoutlier(i))=tmpreason;
            kcandidate=kcandidate+kxoutlier(i);
            semilogx(xvalue(ixoutlier),myvar(ixoutlier),'gd','markersize',12);
        end
        end
      else
         medvalue(i) = 0;
      end
    end
    
    set(gca,'linewidth',1,'fontsize',12)
    xlabel('ngood','fontsize',16)
    ylabel(myvarname,'fontsize',16)
    title([CCDstring ': ' num2str(nsig) '\sigma outliers of ' myvarname],'fontsize',16)
    
    
    iz=find(stdmag<20);     xvalue=stdmag(iz);    myvar=myvar(iz);  REFmag=REF(iz); nearbymag = nearby(iz);
%   xvaluev=stdmagv;
    subplot(2,1,2)
    plot(xvalue,myvar,'.')
    hold on
%   plot(xvaluev,myvarv,'ro','markersize',8)
    xlim([6 15])
    ixvalue=prctile(xvalue',0:100/dN:100)';
    clear medvalue stdvalue stdvalue1 stdvalue2 kxoutlier
    for i=1:length(ixvalue)-1
      if(length(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1))) > 1) 
        medvalue(i)=median(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1)));
        stdvalue(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1)));
        if (length(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue(i))) > 1)
        stdvalue1(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue(i))); % 6-sigma clipping
        stdvalue2(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue1(i))); % 6-sigma clipping
        % stdvalue1(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue2(i))); % 6-sigma clipping
        % stdvalue2(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue1(i))); % 6-sigma clipping
        
        ixoutlier=find( xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & myvar>medvalue(i)+nsig*stdvalue2(i)  & myvar>=rmsthr );  % nsig-sigma outlier
        plot([ixvalue(i) ixvalue(i+1)],[1 1]*medvalue(i)+nsig*stdvalue2(i),'g--')
        kxoutlier(i)=length(ixoutlier);
        if kxoutlier(i)>0
            REFcandidate(kcandidate+1:kcandidate+kxoutlier(i))=REFmag(ixoutlier);
            REFnearby(kcandidate+1:kcandidate+kxoutlier(i))=nearbymag(ixoutlier);
            tmpreason = cell(length(ixoutlier),1);
            tmpreason(:) = myvarname;
            REFreason(kcandidate+1:kcandidate+kxoutlier(i))=tmpreason;
            kcandidate=kcandidate+kxoutlier(i);
            plot(xvalue(ixoutlier),myvar(ixoutlier),'gd','markersize',12);
        end
        end
      else
         medvalue(i) = 0;
      end
    end
    
    set(gca,'linewidth',1,'fontsize',12)
    xlabel(xmaglabel,'fontsize',16)
    ylabel(myvarname,'fontsize',16)
    if (doPlots == 1)
      saveas(gcf,[CCDstring '_' myvarname fqualifier '.eps'],'epsc2')
    endif    
    close

    %if sum(kxoutlier)>0
    %    myvarname
    %    REFcandidate
    %end
    
    % --- 1.2 range1
    dN=8;  % number of bins in ngood and gmag % -------- need to be revised
    myvar=range1;                             % --------
%   myvarv=range1v;                           % --------
    myvarname='range1';                       % --------
    nsig=5; % nsig-sigma outliers             % --------
    
    
    xvalue=ngood;
%   xvaluev=ngoodv;
    subplot(2,1,1)
    semilogx(xvalue,myvar,'.')
    hold on
%   semilogx(xvaluev,myvarv,'ro','markersize',8)
    
    ixvalue=prctile(xvalue',0:100/dN:100)';  % Y = prctile(X,p) calculates a value that is greater than p percent of the values in X. The values of p must lie in the interval [0 100].
    clear medvalue stdvalue stdvalue1 stdvalue2 kxoutlier
    for i=1:length(ixvalue)-1
      if(length(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1))) > 1) 
        medvalue(i)=median(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1)));
        stdvalue(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1)));
        stdvalue1(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue(i))); % 6-sigma clipping
        stdvalue2(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue1(i))); % 6-sigma clipping
        % stdvalue1(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue2(i))); % 6-sigma clipping
        % stdvalue2(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue1(i))); % 6-sigma clipping
        
        ixoutlier=find( xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & myvar>medvalue(i)+nsig*stdvalue2(i) );  % nsig-sigma outlier
        semilogx([ixvalue(i) ixvalue(i+1)],[1 1]*medvalue(i)+nsig*stdvalue2(i),'g--')
        kxoutlier(i)=length(ixoutlier);
        if kxoutlier(i)>0
            REFcandidate(kcandidate+1:kcandidate+kxoutlier(i))=REF(ixoutlier);
            REFnearby(kcandidate+1:kcandidate+kxoutlier(i))=nearby(ixoutlier);
            tmpreason = cell(length(ixoutlier),1);
            tmpreason(:) = myvarname;
            REFreason(kcandidate+1:kcandidate+kxoutlier(i))=tmpreason;
            kcandidate=kcandidate+kxoutlier(i);
            semilogx(xvalue(ixoutlier),myvar(ixoutlier),'gd','markersize',12);
        end
      else
         medvalue(i) = 0;
      end
    end
    
    set(gca,'linewidth',1,'fontsize',12)
    xlabel('ngood','fontsize',16)
    ylabel(myvarname,'fontsize',16)
    title([CCDstring ': ' num2str(nsig) '\sigma outliers of ' myvarname],'fontsize',16)
    
    iz=find(stdmag<20);     xvalue=stdmag(iz);    myvar=myvar(iz);  REFmag=REF(iz); nearbymag = nearby(iz);
%   xvaluev=stdmagv;
    subplot(2,1,2)
    plot(xvalue,myvar,'.')
    hold on
%   plot(xvaluev,myvarv,'ro','markersize',8)
    xlim([6 15])
    ixvalue=prctile(xvalue',0:100/dN:100)';
    clear medvalue stdvalue stdvalue1 stdvalue2 kxoutlier
    for i=1:length(ixvalue)-1
      if(length(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1))) > 1) 
        medvalue(i)=median(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1)));
        stdvalue(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1)));
        stdvalue1(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue(i))); % 6-sigma clipping
        stdvalue2(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue1(i))); % 6-sigma clipping
        % stdvalue1(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue2(i))); % 6-sigma clipping
        % stdvalue2(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue1(i))); % 6-sigma clipping
        
        ixoutlier=find( xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & myvar>medvalue(i)+nsig*stdvalue2(i) );  % nsig-sigma outlier
        plot([ixvalue(i) ixvalue(i+1)],[1 1]*medvalue(i)+nsig*stdvalue2(i),'g--')
        kxoutlier(i)=length(ixoutlier);
        if kxoutlier(i)>0
            REFcandidate(kcandidate+1:kcandidate+kxoutlier(i))=REFmag(ixoutlier);
            REFnearby(kcandidate+1:kcandidate+kxoutlier(i))=nearbymag(ixoutlier);
            tmpreason = cell(length(ixoutlier),1);
            tmpreason(:) = myvarname;
            REFreason(kcandidate+1:kcandidate+kxoutlier(i))=tmpreason;
            kcandidate=kcandidate+kxoutlier(i);
            plot(xvalue(ixoutlier),myvar(ixoutlier),'gd','markersize',12);
        end
      else
         medvalue(i) = 0;
      end
    end
    
    set(gca,'linewidth',1,'fontsize',12)
    xlabel(xmaglabel,'fontsize',16)
    ylabel(myvarname,'fontsize',16)
    if (doPlots == 1)
      saveas(gcf,[CCDstring '_' myvarname fqualifier '.eps'],'epsc2')
    endif    
    close

    %if sum(kxoutlier)>0
    %    myvarname
    %    REFcandidate
    %end
    
    % --- 1.3 range2
    dN=8;  % number of bins in ngood and gmag % -------- need to be revised
    myvar=range2;                             % --------
%   myvarv=range2v;                           % --------
    myvarname='range2';                       % --------
    nsig=5; % nsig-sigma outliers             % --------
    
    
    xvalue=ngood;
%   xvaluev=ngoodv;
    subplot(2,1,1)
    semilogx(xvalue,myvar,'.')
    hold on
%   semilogx(xvaluev,myvarv,'ro','markersize',8)
    
    ixvalue=prctile(xvalue',0:100/dN:100)';  % Y = prctile(X,p) calculates a value that is greater than p percent of the values in X. The values of p must lie in the interval [0 100].
    clear medvalue stdvalue stdvalue1 stdvalue2 kxoutlier
    for i=1:length(ixvalue)-1
      if(length(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1))) > 1) 
        medvalue(i)=median(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1)));
        stdvalue(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1)));
        stdvalue1(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue(i))); % 6-sigma clipping
        stdvalue2(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue1(i))); % 6-sigma clipping
        % stdvalue1(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue2(i))); % 6-sigma clipping
        % stdvalue2(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue1(i))); % 6-sigma clipping
        
        ixoutlier=find( xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & myvar>medvalue(i)+nsig*stdvalue2(i) );  % nsig-sigma outlier
        semilogx([ixvalue(i) ixvalue(i+1)],[1 1]*medvalue(i)+nsig*stdvalue2(i),'g--')
        kxoutlier(i)=length(ixoutlier);
        if kxoutlier(i)>0
            REFcandidate(kcandidate+1:kcandidate+kxoutlier(i))=REF(ixoutlier);
            REFnearby(kcandidate+1:kcandidate+kxoutlier(i))=nearby(ixoutlier);
            tmpreason = cell(length(ixoutlier),1);
            tmpreason(:) = myvarname;
            REFreason(kcandidate+1:kcandidate+kxoutlier(i))=tmpreason;
            kcandidate=kcandidate+kxoutlier(i);
            semilogx(xvalue(ixoutlier),myvar(ixoutlier),'gd','markersize',12);
        end
      else
         medvalue(i) = 0;
      end
    end
    
    set(gca,'linewidth',1,'fontsize',12)
    xlabel('ngood','fontsize',16)
    ylabel(myvarname,'fontsize',16)
    title([CCDstring ': ' num2str(nsig) '\sigma outliers of ' myvarname],'fontsize',16)
    
    iz=find(stdmag<20);     xvalue=stdmag(iz);    myvar=myvar(iz);  REFmag=REF(iz); nearbymag = nearby(iz);
%   xvaluev=stdmagv;
    subplot(2,1,2)
    plot(xvalue,myvar,'.')
    hold on
%   plot(xvaluev,myvarv,'ro','markersize',8)
    xlim([6 15])
    ixvalue=prctile(xvalue',0:100/dN:100)';
    clear medvalue stdvalue stdvalue1 stdvalue2 kxoutlier
    for i=1:length(ixvalue)-1
      if(length(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1))) > 1) 
        medvalue(i)=median(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1)));
        stdvalue(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1)));
        if (length(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue(i))) > 1)
        stdvalue1(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue(i))); % 6-sigma clipping
        stdvalue2(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue1(i))); % 6-sigma clipping
        stdvalue1(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue2(i))); % 6-sigma clipping
        stdvalue2(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue1(i))); % 6-sigma clipping
        
        ixoutlier=find( xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & myvar>medvalue(i)+nsig*stdvalue2(i) );  % nsig-sigma outlier
        plot([ixvalue(i) ixvalue(i+1)],[1 1]*medvalue(i)+nsig*stdvalue2(i),'g--')
        kxoutlier(i)=length(ixoutlier);
        if kxoutlier(i)>0
            REFcandidate(kcandidate+1:kcandidate+kxoutlier(i))=REFmag(ixoutlier);
            REFnearby(kcandidate+1:kcandidate+kxoutlier(i))=nearbymag(ixoutlier);
            tmpreason = cell(length(ixoutlier),1);
            tmpreason(:) = myvarname;
            REFreason(kcandidate+1:kcandidate+kxoutlier(i))=tmpreason;
            kcandidate=kcandidate+kxoutlier(i);
            plot(xvalue(ixoutlier),myvar(ixoutlier),'gd','markersize',12);
        end
        end
      else
         medvalue(i) = 0;
      end
    end
    
    set(gca,'linewidth',1,'fontsize',12)
    xlabel(xmaglabel,'fontsize',16)
    ylabel(myvarname,'fontsize',16)
    if (doPlots == 1)
      saveas(gcf,[CCDstring '_' myvarname fqualifier '.eps'],'epsc2')
    endif    
    close

    %if sum(kxoutlier)>0
    %    myvarname
    %    REFcandidate
    %end
    
    % --- 1.4 lcrms1/lcrms2
    dN=8;  % number of bins in ngood and gmag % -------- need to be revised
    myvar=lcrms1./lcrms2;                             % --------
%   myvarv=lcrms1v./lcrms2v;                           % --------
    myvarname='lcrms1-to-lcrms2';                       % --------
    nsig=8; % nsig-sigma outliers             % --------
    
    
    xvalue=ngood;
%   xvaluev=ngoodv;
    subplot(2,1,1)
    semilogx(xvalue,myvar,'.')
    hold on
%   semilogx(xvaluev,myvarv,'ro','markersize',8)
    
    ixvalue=prctile(xvalue',0:100/dN:100)';  % Y = prctile(X,p) calculates a value that is greater than p percent of the values in X. The values of p must lie in the interval [0 100].
    clear medvalue stdvalue stdvalue1 stdvalue2 kxoutlier
    for i=1:length(ixvalue)-1
      if(length(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1))) > 1) 
        medvalue(i)=median(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1)));
        stdvalue(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1)));
        if (length(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue(i))) > 1)
        stdvalue1(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue(i))); % 6-sigma clipping
        stdvalue2(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue1(i))); % 6-sigma clipping
        % stdvalue1(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue2(i))); % 6-sigma clipping
        % stdvalue2(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue1(i))); % 6-sigma clipping
        
        ixoutlier=find( xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & myvar>medvalue(i)+nsig*stdvalue2(i) );  % nsig-sigma outlier
        semilogx([ixvalue(i) ixvalue(i+1)],[1 1]*medvalue(i)+nsig*stdvalue2(i),'g--')
        kxoutlier(i)=length(ixoutlier);
        if kxoutlier(i)>0
            REFcandidate(kcandidate+1:kcandidate+kxoutlier(i))=REF(ixoutlier);
            REFnearby(kcandidate+1:kcandidate+kxoutlier(i))=nearby(ixoutlier);
            tmpreason = cell(length(ixoutlier),1);
            tmpreason(:) = myvarname;
            REFreason(kcandidate+1:kcandidate+kxoutlier(i))=tmpreason;
            kcandidate=kcandidate+kxoutlier(i);
            semilogx(xvalue(ixoutlier),myvar(ixoutlier),'gd','markersize',12);
        end
        end
      else
         medvalue(i) = 0;
      end
    end
    
    set(gca,'linewidth',1,'fontsize',12)
    xlabel('ngood','fontsize',16)
    ylabel(myvarname,'fontsize',16)
    title([CCDstring ': ' num2str(nsig) '\sigma outliers of ' myvarname],'fontsize',16)
    
    iz=find(stdmag<20);     xvalue=stdmag(iz);    myvar=myvar(iz);  REFmag=REF(iz); nearbymag = nearby(iz);
%   xvaluev=stdmagv;
    subplot(2,1,2)
    plot(xvalue,myvar,'.')
    hold on
%   plot(xvaluev,myvarv,'ro','markersize',8)
    xlim([6 15])
    
    
    set(gca,'linewidth',1,'fontsize',12)
    xlabel(xmaglabel,'fontsize',16)
    ylabel(myvarname,'fontsize',16)
    if (doPlots == 1)
      saveas(gcf,[CCDstring '_' myvarname fqualifier '.eps'],'epsc2')
    endif
    close
    
    %if sum(kxoutlier)>0
    %    myvarname
    %    REFcandidate
    %end
    
    % --- 1.5 lcrms1/lcrms3
    dN=8;  % number of bins in ngood and gmag % -------- need to be revised
    myvar=lcrms1./lcrms3;                             % --------
%   myvarv=lcrms1v./lcrms3v;                           % --------
    myvarname='lcrms1-to-lcrms3';                       % --------
    nsig=8; % nsig-sigma outliers             % --------
    
    
    xvalue=ngood;
%   xvaluev=ngoodv;
    subplot(2,1,1)
    semilogx(xvalue,myvar,'.')
    hold on
%   semilogx(xvaluev,myvarv,'ro','markersize',8)
    
    ixvalue=prctile(xvalue',0:100/dN:100)';  % Y = prctile(X,p) calculates a value that is greater than p percent of the values in X. The values of p must lie in the interval [0 100].
    clear medvalue stdvalue stdvalue1 stdvalue2 kxoutlier
    for i=1:length(ixvalue)-1
      if(length(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1))) > 1) 
        medvalue(i)=median(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1)));
        stdvalue(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1)));
        if (length(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue(i))) > 1)
        stdvalue1(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue(i))); % 6-sigma clipping
        stdvalue2(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue1(i))); % 6-sigma clipping
        % stdvalue1(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue2(i))); % 6-sigma clipping
        % stdvalue2(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue1(i))); % 6-sigma clipping
        
        ixoutlier=find( xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & myvar>medvalue(i)+nsig*stdvalue2(i) );  % nsig-sigma outlier
        semilogx([ixvalue(i) ixvalue(i+1)],[1 1]*medvalue(i)+nsig*stdvalue2(i),'g--')
        kxoutlier(i)=length(ixoutlier);
        if kxoutlier(i)>0
            REFcandidate(kcandidate+1:kcandidate+kxoutlier(i))=REF(ixoutlier);
            REFnearby(kcandidate+1:kcandidate+kxoutlier(i))=nearby(ixoutlier);
            tmpreason = cell(length(ixoutlier),1);
            tmpreason(:) = myvarname;
            REFreason(kcandidate+1:kcandidate+kxoutlier(i))=tmpreason;
            kcandidate=kcandidate+kxoutlier(i);
            semilogx(xvalue(ixoutlier),myvar(ixoutlier),'gd','markersize',12);
        end
        end
      else
         medvalue(i) = 0;
      end
    end
    
    set(gca,'linewidth',1,'fontsize',12)
    xlabel('ngood','fontsize',16)
    ylabel(myvarname,'fontsize',16)
    title([CCDstring ': ' num2str(nsig) '\sigma outliers of ' myvarname],'fontsize',16)
    
    iz=find(stdmag<20);     xvalue=stdmag(iz);    myvar=myvar(iz);  REFmag=REF(iz); nearbymag = nearby(iz);
%   xvaluev=stdmagv;
    subplot(2,1,2)
    plot(xvalue,myvar,'.')
    hold on
%   plot(xvaluev,myvarv,'ro','markersize',8)
    xlim([6 15])
    
    
    set(gca,'linewidth',1,'fontsize',12)
    xlabel(xmaglabel,'fontsize',16)
    ylabel(myvarname,'fontsize',16)
    if (doPlots == 1)
      saveas(gcf,[CCDstring '_' myvarname fqualifier '.eps'],'epsc2')
    endif    
    close

    %if sum(kxoutlier)>0
    %    myvarname
    %    REFcandidate
    %end
    
    % --- 1.6 lcrms1/lcrms4
    dN=8;  % number of bins in ngood and gmag % -------- need to be revised
    myvar=lcrms1./lcrms4;                             % --------
%   myvarv=lcrms1v./lcrms4v;                           % --------
    myvarname='lcrms1-to-lcrms4';                       % --------
    nsig=6; % nsig-sigma outliers             % --------
    
    
    xvalue=ngood;
%   xvaluev=ngoodv;
    subplot(2,1,1)
    semilogx(xvalue,myvar,'.')
    hold on
%   semilogx(xvaluev,myvarv,'ro','markersize',8)
    
    ixvalue=prctile(xvalue',0:100/dN:100)';  % Y = prctile(X,p) calculates a value that is greater than p percent of the values in X. The values of p must lie in the interval [0 100].
    clear medvalue stdvalue stdvalue1 stdvalue2 kxoutlier
    for i=1:length(ixvalue)-1
      if(length(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1))) > 1) 
        medvalue(i)=median(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1)));
        stdvalue(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1)));
        if (length(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue(i))) > 1)
        stdvalue1(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue(i))); % 6-sigma clipping
        stdvalue2(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue1(i))); % 6-sigma clipping
        % stdvalue1(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue2(i))); % 6-sigma clipping
        % stdvalue2(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue1(i))); % 6-sigma clipping
        
        ixoutlier=find( xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & myvar>medvalue(i)+nsig*stdvalue2(i) );  % nsig-sigma outlier
        semilogx([ixvalue(i) ixvalue(i+1)],[1 1]*medvalue(i)+nsig*stdvalue2(i),'g--')
        kxoutlier(i)=length(ixoutlier);
        if kxoutlier(i)>0
            REFcandidate(kcandidate+1:kcandidate+kxoutlier(i))=REF(ixoutlier);
            REFnearby(kcandidate+1:kcandidate+kxoutlier(i))=nearby(ixoutlier);
            tmpreason = cell(length(ixoutlier),1);
            tmpreason(:) = myvarname;
            REFreason(kcandidate+1:kcandidate+kxoutlier(i))=tmpreason;
            kcandidate=kcandidate+kxoutlier(i);
            semilogx(xvalue(ixoutlier),myvar(ixoutlier),'gd','markersize',12);
        end
        end
      else
         medvalue(i) = 0;
      end
    end
    
    set(gca,'linewidth',1,'fontsize',12)
    xlabel('ngood','fontsize',16)
    ylabel(myvarname,'fontsize',16)
    title([CCDstring ': ' num2str(nsig) '\sigma outliers of ' myvarname],'fontsize',16)
    
 
    
    iz=find(stdmag<20);     xvalue=stdmag(iz);    myvar=myvar(iz);  REFmag=REF(iz); nearbymag = nearby(iz);
%   xvaluev=stdmagv;
    subplot(2,1,2)
    plot(xvalue,myvar,'.')
    hold on
%   plot(xvaluev,myvarv,'ro','markersize',8)
    xlim([6 15])
  
    set(gca,'linewidth',1,'fontsize',12)
    xlabel(xmaglabel,'fontsize',16)
    ylabel(myvarname,'fontsize',16)
    if (doPlots == 1)
      saveas(gcf,[CCDstring '_' myvarname fqualifier '.eps'],'epsc2')
    endif    
    close

    %if sum(kxoutlier)>0
    %    myvarname
    %    REFcandidate
    %end
    
    % --- 1.7 lcrms1/lcrms5
    dN=8;  % number of bins in ngood and gmag % -------- need to be revised
    myvar=lcrms1./lcrms5;                             % --------
%   myvarv=lcrms1v./lcrms5v;                           % --------
    myvarname='lcrms1-to-lcrms5';                       % --------
    nsig=6; % nsig-sigma outliers             % --------
    
    
    xvalue=ngood;
%   xvaluev=ngoodv;
    subplot(2,1,1)
    semilogx(xvalue,myvar,'.')
    hold on
%   semilogx(xvaluev,myvarv,'ro','markersize',8)
    
    ixvalue=prctile(xvalue',0:100/dN:100)';  % Y = prctile(X,p) calculates a value that is greater than p percent of the values in X. The values of p must lie in the interval [0 100].
    clear medvalue stdvalue stdvalue1 stdvalue2 kxoutlier
    for i=1:length(ixvalue)-1
      if(length(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1))) > 1) 
        medvalue(i)=median(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1)));
        stdvalue(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1)));
        if (length(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue(i))) > 1)
        stdvalue1(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue(i))); % 6-sigma clipping
        stdvalue2(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue1(i))); % 6-sigma clipping
        % stdvalue1(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue2(i))); % 6-sigma clipping
        % stdvalue2(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue1(i))); % 6-sigma clipping
        
        ixoutlier=find( xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & myvar>medvalue(i)+nsig*stdvalue2(i) );  % nsig-sigma outlier
        semilogx([ixvalue(i) ixvalue(i+1)],[1 1]*medvalue(i)+nsig*stdvalue2(i),'g--')
        kxoutlier(i)=length(ixoutlier);
        if kxoutlier(i)>0
            REFcandidate(kcandidate+1:kcandidate+kxoutlier(i))=REF(ixoutlier);
            REFnearby(kcandidate+1:kcandidate+kxoutlier(i))=nearby(ixoutlier);
            tmpreason = cell(length(ixoutlier),1);
            tmpreason(:) = myvarname;
            REFreason(kcandidate+1:kcandidate+kxoutlier(i))=tmpreason;
            kcandidate=kcandidate+kxoutlier(i);
            semilogx(xvalue(ixoutlier),myvar(ixoutlier),'gd','markersize',12);
        end
        end
      else
         medvalue(i) = 0;
      end
    end
    
    set(gca,'linewidth',1,'fontsize',12)
    xlabel('ngood','fontsize',16)
    ylabel(myvarname,'fontsize',16)
    title([CCDstring ': ' num2str(nsig) '\sigma outliers of ' myvarname],'fontsize',16)
    
    iz=find(stdmag<20);     xvalue=stdmag(iz);    myvar=myvar(iz);  REFmag=REF(iz); nearbymag = nearby(iz);
%   xvaluev=stdmagv;
    subplot(2,1,2)
    plot(xvalue,myvar,'.')
    hold on
%   plot(xvaluev,myvarv,'ro','markersize',8)
    xlim([6 15])
    
    
    set(gca,'linewidth',1,'fontsize',12)
    xlabel(xmaglabel,'fontsize',16)
    ylabel(myvarname,'fontsize',16)
    if (doPlots == 1)
      saveas(gcf,[CCDstring '_' myvarname fqualifier '.eps'],'epsc2')
    endif
    close
    
    %if sum(kxoutlier)>0
    %    myvarname
    %    REFcandidate
    %end
    
    % ------------------------------------------
    % --- Part 2. ndev
    % ------------------------------------------
    % --- 2.1 ndev2
    dN=8;  % number of bins in ngood and gmag % -------- need to be revised
    myvar=ndev2;                             % --------
%   myvarv=ndev2v;                           % --------
    myvarname='ndev2';                       % --------
    nsig=8; % nsig-sigma outliers             % --------
    nthr=20; % threshold of n             % --------
    
    
    xvalue=ngood;
%   xvaluev=ngoodv;
    subplot(2,1,1)
    semilogx(xvalue,myvar,'.')
    hold on
%   semilogx(xvaluev,myvarv,'ro','markersize',8)
    
    ixvalue=prctile(xvalue',0:100/dN:100)';  % Y = prctile(X,p) calculates a value that is greater than p percent of the values in X. The values of p must lie in the interval [0 100].
    clear medvalue stdvalue stdvalue1 stdvalue2 kxoutlier
    for i=1:length(ixvalue)-1
      if(length(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1))) > 1) 
        medvalue(i)=mean(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1)));
        stdvalue(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1)));
        if (length(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue(i))) > 1)
            stdvalue1(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue(i))); % 6-sigma clipping
            stdvalue2(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue1(i))); % 6-sigma clipping
            % stdvalue1(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue2(i))); % 6-sigma clipping
            % stdvalue2(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue1(i))); % 6-sigma clipping
            
            ixoutlier=find( xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & (myvar>medvalue(i)+nsig*stdvalue2(i) | myvar>=nthr ));  % nsig-sigma outlier or >=nthr
            semilogx([ixvalue(i) ixvalue(i+1)],[1 1]*medvalue(i)+nsig*stdvalue2(i),'g--')
            kxoutlier(i)=length(ixoutlier);
            if kxoutlier(i)>0
                REFcandidate(kcandidate+1:kcandidate+kxoutlier(i))=REF(ixoutlier);
                REFnearby(kcandidate+1:kcandidate+kxoutlier(i))=nearby(ixoutlier);
                tmpreason = cell(length(ixoutlier),1);
                tmpreason(:) = myvarname;
                REFreason(kcandidate+1:kcandidate+kxoutlier(i))=tmpreason;
                kcandidate=kcandidate+kxoutlier(i);
                semilogx(xvalue(ixoutlier),myvar(ixoutlier),'gd','markersize',12);
            end
        else
            medvalue(i) = 0;
        end
      else
         medvalue(i) = 0;
      end
    end
    
    set(gca,'linewidth',1,'fontsize',12)
    xlabel('ngood','fontsize',16)
    ylabel(myvarname,'fontsize',16)
    title([CCDstring ': ' num2str(nsig) '\sigma outliers of ' myvarname],'fontsize',16)
    
    iz=find(stdmag<20);     xvalue=stdmag(iz);    myvar=myvar(iz);  REFmag=REF(iz); nearbymag = nearby(iz);
%   xvaluev=stdmagv;
    subplot(2,1,2)
    plot(xvalue,myvar,'.')
    hold on
%   plot(xvaluev,myvarv,'ro','markersize',8)
    xlim([6 15])
    ixvalue=prctile(xvalue',0:100/dN:100)';
    clear medvalue stdvalue stdvalue1 stdvalue2 kxoutlier
    for i=1:length(ixvalue)-1
      if(length(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1))) > 1) 
        medvalue(i)=mean(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1)));
        stdvalue(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1)));
        if (length(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue(i))) > 1)
        stdvalue1(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue(i))); % 6-sigma clipping
        stdvalue2(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue1(i))); % 6-sigma clipping
        % stdvalue1(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue2(i))); % 6-sigma clipping
        % stdvalue2(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue1(i))); % 6-sigma clipping
        
        ixoutlier=find( xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & (myvar>medvalue(i)+nsig*stdvalue2(i) | myvar>=nthr));  % nsig-sigma outlier or >=nthr
        plot([ixvalue(i) ixvalue(i+1)],[1 1]*medvalue(i)+nsig*stdvalue2(i),'g--')
        kxoutlier(i)=length(ixoutlier);
        if kxoutlier(i)>0
            REFcandidate(kcandidate+1:kcandidate+kxoutlier(i))=REFmag(ixoutlier);
            REFnearby(kcandidate+1:kcandidate+kxoutlier(i))=nearbymag(ixoutlier);
            tmpreason = cell(length(ixoutlier),1);
            tmpreason(:) = myvarname;
            REFreason(kcandidate+1:kcandidate+kxoutlier(i))=tmpreason;
            kcandidate=kcandidate+kxoutlier(i);
            plot(xvalue(ixoutlier),myvar(ixoutlier),'gd','markersize',12);
        end
        end
      else
         medvalue(i) = 0;
      end
    end
    
    set(gca,'linewidth',1,'fontsize',12)
    xlabel(xmaglabel,'fontsize',16)
    ylabel(myvarname,'fontsize',16)
    if (doPlots == 1)
      saveas(gcf,[CCDstring '_' myvarname fqualifier '.eps'],'epsc2')
    endif    
    close

    %if sum(kxoutlier)>0
    %    myvarname
    %    REFcandidate
    %end
    
    % --- 2.2 ndev3
    dN=8;  % number of bins in ngood and gmag % -------- need to be revised
    myvar=ndev3;                             % --------
%   myvarv=ndev3v;                           % --------
    myvarname='ndev3';                       % --------
    nsig=8; % nsig-sigma outliers             % --------
    nthr=5; % threshold of n; all object with ndev3>=5 got picked             % --------
    nthrlow=3;  % all objects with ndec3<3 got rejected
    
    xvalue=ngood;
%   xvaluev=ngoodv;
    subplot(2,1,1)
    semilogx(xvalue,myvar,'.')
    hold on
%   semilogx(xvaluev,myvarv,'ro','markersize',8)
    
    ixvalue=prctile(xvalue',0:100/dN:100)';  % Y = prctile(X,p) calculates a value that is greater than p percent of the values in X. The values of p must lie in the interval [0 100].
    clear medvalue stdvalue stdvalue1 stdvalue2 kxoutlier
    for i=1:length(ixvalue)-1
      if(length(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1))) > 1) 
        medvalue(i)=mean(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1)));
        stdvalue(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1)));
        if (length(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue(i))) > 1)
        stdvalue1(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue(i))); % 6-sigma clipping
        if (length(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue1(i))) > 1)
          stdvalue2(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue1(i))); % 6-sigma clipping
          % stdvalue1(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue2(i))); % 6-sigma clipping
          % stdvalue2(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue1(i))); % 6-sigma clipping
          
          ixoutlier=find( xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & (myvar>medvalue(i)+nsig*stdvalue2(i) | myvar>=nthr ) & myvar>=nthrlow);  % nsig-sigma outlier or >=nthr
          semilogx([ixvalue(i) ixvalue(i+1)],[1 1]*medvalue(i)+nsig*stdvalue2(i),'g--')
          kxoutlier(i)=length(ixoutlier);
          if kxoutlier(i)>0
              REFcandidate(kcandidate+1:kcandidate+kxoutlier(i))=REF(ixoutlier);
              REFnearby(kcandidate+1:kcandidate+kxoutlier(i))=nearby(ixoutlier);
              tmpreason = cell(length(ixoutlier),1);
              tmpreason(:) = myvarname;
              REFreason(kcandidate+1:kcandidate+kxoutlier(i))=tmpreason;
              kcandidate=kcandidate+kxoutlier(i);
              semilogx(xvalue(ixoutlier),myvar(ixoutlier),'gd','markersize',12);
          end
          end
        else
           medvalue(i) = 0;
        end

      else
         medvalue(i) = 0;
      end
    end

    set(gca,'linewidth',1,'fontsize',12)
    xlabel('ngood','fontsize',16)
    ylabel(myvarname,'fontsize',16)
    title([CCDstring ': ' num2str(nsig) '\sigma outliers of ' myvarname],'fontsize',16)
    
    iz=find(stdmag<20);     xvalue=stdmag(iz);    myvar=myvar(iz);  REFmag=REF(iz); nearbymag = nearby(iz);
%   xvaluev=stdmagv;
    subplot(2,1,2)
    plot(xvalue,myvar,'.')
    hold on
%   plot(xvaluev,myvarv,'ro','markersize',8)
    xlim([6 15])
    ixvalue=prctile(xvalue',0:100/dN:100)';
    clear medvalue stdvalue stdvalue1 stdvalue2 kxoutlier
    for i=1:length(ixvalue)-1
      if(length(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1))) > 1) 
        medvalue(i)=mean(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1)));
        stdvalue(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1)));
        if (length(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue(i))) > 1)
        stdvalue1(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue(i))); % 6-sigma clipping
        if (length(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue1(i))) > 1)
          stdvalue2(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue1(i))); % 6-sigma clipping
          % stdvalue1(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue2(i))); % 6-sigma clipping
          % stdvalue2(i)=std(myvar(xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & abs(myvar-medvalue(i))<6*stdvalue1(i))); % 6-sigma clipping
          
          ixoutlier=find( xvalue>=ixvalue(i) & xvalue<ixvalue(i+1) & (myvar>medvalue(i)+nsig*stdvalue2(i) | myvar>=nthr) & myvar>=nthrlow);  % nsig-sigma outlier or >=nthr
          plot([ixvalue(i) ixvalue(i+1)],[1 1]*medvalue(i)+nsig*stdvalue2(i),'g--')
          kxoutlier(i)=length(ixoutlier);
          if kxoutlier(i)>0
              REFcandidate(kcandidate+1:kcandidate+kxoutlier(i))=REFmag(ixoutlier);
              REFnearby(kcandidate+1:kcandidate+kxoutlier(i))=nearbymag(ixoutlier);
              tmpreason = cell(length(ixoutlier),1);
              tmpreason(:) = myvarname;
              REFreason(kcandidate+1:kcandidate+kxoutlier(i))=tmpreason;
              kcandidate=kcandidate+kxoutlier(i);
              plot(xvalue(ixoutlier),myvar(ixoutlier),'gd','markersize',12);
          end
          end
        else
           medvalue(i) = 0;
        end
      else
         medvalue(i) = 0;
      end
    end
    
    set(gca,'linewidth',1,'fontsize',12)
    xlabel(xmaglabel,'fontsize',16)
    ylabel(myvarname,'fontsize',16)
    if (doPlots == 1)
      saveas(gcf,[CCDstring '_' myvarname fqualifier '.eps'],'epsc2')
    endif    
    close

    %if sum(kxoutlier)>0
    %    myvarname
    %    REFcandidate
    %end
    
    % ------------------------------------------
    % --- Part 3. nadj --- note I'm starting using nthr
    % ------------------------------------------
    % --- 3.1 nadj1
    myvar=nadj1;                             % --------
%   myvarv=nadj1v;                           % --------
    myvarname='nadj1';                       % --------
    ix=find(nadj1<mean(nadj1)+20*std(nadj1));
    if (length(ix) > 0)
    nthr=max(2,ceil(mean(nadj1(ix))+10*std(nadj1(ix))));  % threshold of n             % --------

    if (length(nthr) > 0)
        xvalue=ngood;
    %   xvaluev=ngoodv;
        subplot(2,1,1)
        semilogx(xvalue,myvar,'.')
        hold on
    %   semilogx(xvaluev,myvarv,'ro','markersize',8)
        
        ixoutlier=find(myvar>=nthr);
        kxoutlier=length(ixoutlier);
        if kxoutlier>0
            REFcandidate(kcandidate+1:kcandidate+kxoutlier)=REF(ixoutlier);
            REFnearby(kcandidate+1:kcandidate+kxoutlier)=nearby(ixoutlier);
            tmpreason = cell(length(ixoutlier),1);
            tmpreason(:) = myvarname;
            REFreason(kcandidate+1:kcandidate+kxoutlier)=tmpreason;
            kcandidate=kcandidate+kxoutlier;
            semilogx(xvalue(ixoutlier),myvar(ixoutlier),'gd','markersize',12);
        end
        
        set(gca,'linewidth',1,'fontsize',12)
        xlabel('ngood','fontsize',16)
        ylabel(myvarname,'fontsize',16)
        title([CCDstring ': ' myvarname '>= ' num2str(nthr)],'fontsize',16)
        
        xvalue=stdmag;
    %   xvaluev=stdmagv;
        subplot(2,1,2)
        plot(xvalue,myvar,'.')
        hold on
    %   plot(xvaluev,myvarv,'ro','markersize',8)
        xlim([6 15])
        if kxoutlier>0
            plot(xvalue(ixoutlier),myvar(ixoutlier),'gd','markersize',12);
        end
        set(gca,'linewidth',1,'fontsize',12)
        xlabel(xmaglabel,'fontsize',16)
        ylabel(myvarname,'fontsize',16)
        if (doPlots == 1)
          saveas(gcf,[CCDstring '_' myvarname fqualifier '.eps'],'epsc2')
        endif
        close

    end
    end
    %if sum(kxoutlier)>0
    %    myvarname
    %    REFcandidate
    %end
    
    % --- 3.2 nadj2
    myvar=nadj2;                             % --------
%   myvarv=nadj2v;                           % --------
    myvarname='nadj2';                       % --------
    nthr=1; % threshold of n             % --------
    
    xvalue=ngood;
%   xvaluev=ngoodv;
    subplot(2,1,1)
    semilogx(xvalue,myvar,'.')
    hold on
%   semilogx(xvaluev,myvarv,'ro','markersize',8)
    
    ixoutlier=find(myvar>=nthr);
    kxoutlier=length(ixoutlier);
    if kxoutlier>0
        REFcandidate(kcandidate+1:kcandidate+kxoutlier)=REF(ixoutlier);
        REFnearby(kcandidate+1:kcandidate+kxoutlier)=nearby(ixoutlier);
        tmpreason = cell(length(ixoutlier),1);
        tmpreason(:) = myvarname;
        REFreason(kcandidate+1:kcandidate+kxoutlier)=tmpreason;
        kcandidate=kcandidate+kxoutlier;
        semilogx(xvalue(ixoutlier),myvar(ixoutlier),'gd','markersize',12);
    end

    set(gca,'linewidth',1,'fontsize',12)
    xlabel('ngood','fontsize',16)
    ylabel(myvarname,'fontsize',16)
    title([CCDstring ': ' myvarname '>= ' num2str(nthr)],'fontsize',16)
    
    xvalue=stdmag;
%   xvaluev=stdmagv;
    subplot(2,1,2)
    plot(xvalue,myvar,'.')
    hold on
%   plot(xvaluev,myvarv,'ro','markersize',8)
    xlim([6 15])
    if kxoutlier>0
        plot(xvalue(ixoutlier),myvar(ixoutlier),'gd','markersize',12);
    end
    set(gca,'linewidth',1,'fontsize',12)
    xlabel(xmaglabel,'fontsize',16)
    ylabel(myvarname,'fontsize',16)
    if (doPlots == 1)
      saveas(gcf,[CCDstring '_' myvarname fqualifier '.eps'],'epsc2')
    endif    
    close

    %if sum(kxoutlier)>0
    %    myvarname
    %    REFcandidate
    %end
    
    % --- 3.3 nadj3
    myvar=nadj3;                             % --------
%   myvarv=nadj3v;                           % --------
    myvarname='nadj3';                       % --------
    nthr=1; % threshold of n             % --------
    
    xvalue=ngood;
%   xvaluev=ngoodv;
    subplot(2,1,1)
    semilogx(xvalue,myvar,'.')
    hold on
%   semilogx(xvaluev,myvarv,'ro','markersize',8)
    
    ixoutlier=find(myvar>=nthr);
    kxoutlier=length(ixoutlier);
    if kxoutlier>0
        REFcandidate(kcandidate+1:kcandidate+kxoutlier)=REF(ixoutlier);
        REFnearby(kcandidate+1:kcandidate+kxoutlier)=nearby(ixoutlier);
        tmpreason = cell(length(ixoutlier),1);
        tmpreason(:) = myvarname;
        REFreason(kcandidate+1:kcandidate+kxoutlier)=tmpreason;
        kcandidate=kcandidate+kxoutlier;
        semilogx(xvalue(ixoutlier),myvar(ixoutlier),'gd','markersize',12);
    end

    set(gca,'linewidth',1,'fontsize',12)
    xlabel('ngood','fontsize',16)
    ylabel(myvarname,'fontsize',16)
    title([CCDstring ': ' myvarname '>= ' num2str(nthr)],'fontsize',16)
    
    xvalue=stdmag;
%   xvaluev=stdmagv;
    subplot(2,1,2)
    plot(xvalue,myvar,'.')
    hold on
%   plot(xvaluev,myvarv,'ro','markersize',8)
    xlim([6 15])
    if kxoutlier>0
        plot(xvalue(ixoutlier),myvar(ixoutlier),'gd','markersize',12);
    end
    set(gca,'linewidth',1,'fontsize',12)
    xlabel(xmaglabel,'fontsize',16)
    ylabel(myvarname,'fontsize',16)
    if (doPlots == 1)
      saveas(gcf,[CCDstring '_' myvarname fqualifier '.eps'],'epsc2')
    endif
    close
    
    %if sum(kxoutlier)>0
    %    myvarname
    %    REFcandidate
    %end
    
    % ------------------------------------------
    % --- Part 4. nburst
    % ------------------------------------------
    % --- 4.1 nburst1
    myvar=nburst1;                             % --------
%   myvarv=nburst1v;                           % --------
    myvarname='nburst1';                       % --------
    nthr=2; % threshold of n             % --------
    
    xvalue=ngood;
%   xvaluev=ngoodv;
    subplot(2,1,1)
    semilogx(xvalue,myvar,'.')
    hold on
%   semilogx(xvaluev,myvarv,'ro','markersize',8)
    
    ixoutlier=find(myvar>=nthr);
    kxoutlier=length(ixoutlier);
    if kxoutlier>0
        REFcandidate(kcandidate+1:kcandidate+kxoutlier)=REF(ixoutlier);
        REFnearby(kcandidate+1:kcandidate+kxoutlier)=nearby(ixoutlier);
        tmpreason = cell(length(ixoutlier),1);
        tmpreason(:) = myvarname;
        REFreason(kcandidate+1:kcandidate+kxoutlier)=tmpreason;
        kcandidate=kcandidate+kxoutlier;
        semilogx(xvalue(ixoutlier),myvar(ixoutlier),'gd','markersize',12);
    end

    set(gca,'linewidth',1,'fontsize',12)
    xlabel('ngood','fontsize',16)
    ylabel(myvarname,'fontsize',16)
    title([CCDstring ': ' myvarname '>= ' num2str(nthr)],'fontsize',16)
    
    xvalue=stdmag;
%   xvaluev=stdmagv;
    subplot(2,1,2)
    plot(xvalue,myvar,'.')
    hold on
%   plot(xvaluev,myvarv,'ro','markersize',8)
    xlim([6 15])
    if kxoutlier>0
        plot(xvalue(ixoutlier),myvar(ixoutlier),'gd','markersize',12);
    end
    set(gca,'linewidth',1,'fontsize',12)
    xlabel(xmaglabel,'fontsize',16)
    ylabel(myvarname,'fontsize',16)
    if (doPlots == 1)
      saveas(gcf,[CCDstring '_' myvarname fqualifier '.eps'],'epsc2')
    endif
    close
    
    % --- 4.2 nburst2
    myvar=nburst2;                             % --------
%   myvarv=nburst2v;                           % --------
    myvarname='nburst2';                       % --------
    ix=find(nburst2<mean(nburst2)+20*std(nburst2));
    if (length(ix) > 0)
    nthr=max(2,ceil(mean(nburst2(ix))+10*std(nburst2(ix))));  % threshold of n             % --------
    
    if (length(nthr) > 0)
        xvalue=ngood;
    %   xvaluev=ngoodv;
        subplot(2,1,1)
        semilogx(xvalue,myvar,'.')
        hold on
    %   semilogx(xvaluev,myvarv,'ro','markersize',8)
        
        ixoutlier=find(myvar>=nthr);
        kxoutlier=length(ixoutlier);
        if kxoutlier>0
            REFcandidate(kcandidate+1:kcandidate+kxoutlier)=REF(ixoutlier);
            REFnearby(kcandidate+1:kcandidate+kxoutlier)=nearby(ixoutlier);
            tmpreason = cell(length(ixoutlier),1);
            tmpreason(:) = myvarname;
            REFreason(kcandidate+1:kcandidate+kxoutlier)=tmpreason;
            kcandidate=kcandidate+kxoutlier;
            semilogx(xvalue(ixoutlier),myvar(ixoutlier),'gd','markersize',12);
        end
        set(gca,'linewidth',1,'fontsize',12)
        xlabel('ngood','fontsize',16)
        ylabel(myvarname,'fontsize',16)
        title([CCDstring ': ' myvarname '>= ' num2str(nthr)],'fontsize',16)
        
        xvalue=stdmag;
    %   xvaluev=stdmagv;
        subplot(2,1,2)
        plot(xvalue,myvar,'.')
        hold on
    %   plot(xvaluev,myvarv,'ro','markersize',8)
        xlim([6 15])
        if kxoutlier>0
            plot(xvalue(ixoutlier),myvar(ixoutlier),'gd','markersize',12);
        end
        set(gca,'linewidth',1,'fontsize',12)
        xlabel(xmaglabel,'fontsize',16)
        ylabel(myvarname,'fontsize',16)
        if (doPlots == 1)
          saveas(gcf,[CCDstring '_' myvarname fqualifier '.eps'],'epsc2')
        endif
        close
    end    
    end
    % --- 4.3 nburst3
    myvar=nburst3;                             % --------
%   myvarv=nburst3v;                           % --------
    myvarname='nburst3';                       % --------
    ix=find(nburst3<mean(nburst3)+20*std(nburst3));
    if (length(ix) > 0)
    nthr=max(2,ceil(mean(nburst3(ix))+10*std(nburst3(ix))));  % threshold of n             % --------
    
    if (length(nthr) > 0)
        xvalue=ngood;
    %   xvaluev=ngoodv;
        subplot(2,1,1)
        semilogx(xvalue,myvar,'.')
        hold on
    %   semilogx(xvaluev,myvarv,'ro','markersize',8)
        
        ixoutlier=find(myvar>=nthr);
        kxoutlier=length(ixoutlier);
        if kxoutlier>0
            REFcandidate(kcandidate+1:kcandidate+kxoutlier)=REF(ixoutlier);
            REFnearby(kcandidate+1:kcandidate+kxoutlier)=nearby(ixoutlier);
            tmpreason = cell(length(ixoutlier),1);
            tmpreason(:) = myvarname;
            REFreason(kcandidate+1:kcandidate+kxoutlier)=tmpreason;
            kcandidate=kcandidate+kxoutlier;
            semilogx(xvalue(ixoutlier),myvar(ixoutlier),'gd','markersize',12);
        end

        set(gca,'linewidth',1,'fontsize',12)
        xlabel('ngood','fontsize',16)
        ylabel(myvarname,'fontsize',16)
        title([CCDstring ': ' myvarname '>= ' num2str(nthr)],'fontsize',16)
        
        xvalue=stdmag;
    %   xvaluev=stdmagv;
        subplot(2,1,2)
        plot(xvalue,myvar,'.')
        hold on
    %   plot(xvaluev,myvarv,'ro','markersize',8)
        xlim([6 15])
        if kxoutlier>0
            plot(xvalue(ixoutlier),myvar(ixoutlier),'gd','markersize',12);
        end
        set(gca,'linewidth',1,'fontsize',12)
        xlabel(xmaglabel,'fontsize',16)
        ylabel(myvarname,'fontsize',16)
        if (doPlots == 1)
          saveas(gcf,[CCDstring '_' myvarname fqualifier '.eps'],'epsc2')
        endif
        close
    end
    end
    % --- 4.4 nburst4
    myvar=nburst4;                             % --------
%   myvarv=nburst4v;                           % --------
    myvarname='nburst4';                       % --------
    ix=find(nburst4<mean(nburst4)+20*std(nburst4));
    if (length(ix) > 0)
    nthr=max(3,ceil(mean(nburst4(ix))+10*std(nburst4(ix))));  % threshold of n             % --------
    
    if (length(nthr) > 0)
        xvalue=ngood;
    %   xvaluev=ngoodv;
        subplot(2,1,1)
        semilogx(xvalue,myvar,'.')
        hold on
    %   semilogx(xvaluev,myvarv,'ro','markersize',8)
        
        ixoutlier=find(myvar>=nthr);
        kxoutlier=length(ixoutlier);
        if kxoutlier>0
            REFcandidate(kcandidate+1:kcandidate+kxoutlier)=REF(ixoutlier);
            REFnearby(kcandidate+1:kcandidate+kxoutlier)=nearby(ixoutlier);
            tmpreason = cell(length(ixoutlier),1);
            tmpreason(:) = myvarname;
            REFreason(kcandidate+1:kcandidate+kxoutlier)=tmpreason;
            kcandidate=kcandidate+kxoutlier;
            semilogx(xvalue(ixoutlier),myvar(ixoutlier),'gd','markersize',12);
        end
        set(gca,'linewidth',1,'fontsize',12)
        xlabel('ngood','fontsize',16)
        ylabel(myvarname,'fontsize',16)
        title([CCDstring ': ' myvarname '>= ' num2str(nthr)],'fontsize',16)
        
        xvalue=stdmag;
    %   xvaluev=stdmagv;
        subplot(2,1,2)
        plot(xvalue,myvar,'.')
        hold on
    %   plot(xvaluev,myvarv,'ro','markersize',8)
        xlim([6 15])
        if kxoutlier>0
            plot(xvalue(ixoutlier),myvar(ixoutlier),'gd','markersize',12);
        end
        set(gca,'linewidth',1,'fontsize',12)
        xlabel(xmaglabel,'fontsize',16)
        ylabel(myvarname,'fontsize',16)
        if (doPlots == 1)
          saveas(gcf,[CCDstring '_' myvarname fqualifier '.eps'],'epsc2')
        endif   
        close 
    end
    end
    %if sum(kxoutlier)>0
    %    myvarname
    %    REFcandidate
    %end
    
    % ------------------------------------------
    % --- Part 5. ndip
    % ------------------------------------------
    % --- 5.1 ndip1
    myvar=ndip1;                             % --------
%   myvarv=ndip1v;                           % --------
    myvarname='ndip1';                       % --------
    nthr=2; % threshold of n             % --------
    
    xvalue=ngood;
%   xvaluev=ngoodv;
    subplot(2,1,1)
    semilogx(xvalue,myvar,'.')
    hold on
%   semilogx(xvaluev,myvarv,'ro','markersize',8)
    
    ixoutlier=find(myvar>=nthr);
    kxoutlier=length(ixoutlier);
    if kxoutlier>0
        REFcandidate(kcandidate+1:kcandidate+kxoutlier)=REF(ixoutlier);
        REFnearby(kcandidate+1:kcandidate+kxoutlier)=nearby(ixoutlier);
        tmpreason = cell(length(ixoutlier),1);
        tmpreason(:) = myvarname;
        REFreason(kcandidate+1:kcandidate+kxoutlier)=tmpreason;
        kcandidate=kcandidate+kxoutlier;
        semilogx(xvalue(ixoutlier),myvar(ixoutlier),'gd','markersize',12);
    end
    set(gca,'linewidth',1,'fontsize',12)
    xlabel('ngood','fontsize',16)
    ylabel(myvarname,'fontsize',16)
    title([CCDstring ': ' myvarname '>= ' num2str(nthr)],'fontsize',16)
    
    xvalue=stdmag;
%   xvaluev=stdmagv;
    subplot(2,1,2)
    plot(xvalue,myvar,'.')
    hold on
%   plot(xvaluev,myvarv,'ro','markersize',8)
    xlim([6 15])
    if kxoutlier>0
        plot(xvalue(ixoutlier),myvar(ixoutlier),'gd','markersize',12);
    end
    set(gca,'linewidth',1,'fontsize',12)
    xlabel(xmaglabel,'fontsize',16)
    ylabel(myvarname,'fontsize',16)
    if (doPlots == 1)
      saveas(gcf,[CCDstring '_' myvarname fqualifier '.eps'],'epsc2')
    endif
    close
    
    % --- 5.2 ndip2
    myvar=ndip2;                             % --------
%   myvarv=ndip2v;                           % --------
    myvarname='ndip2';                       % --------
    ix=find(ndip2<mean(ndip2)+20*std(ndip2));
    if (length(ix) > 0) 
        nthr=max(2,ceil(mean(ndip2(ix))+10*std(ndip2(ix))));  % threshold of n             % --------
        
        if (length(nthr) > 0)
            xvalue=ngood;
        %   xvaluev=ngoodv;
            subplot(2,1,1)
            semilogx(xvalue,myvar,'.')
            hold on
        %   semilogx(xvaluev,myvarv,'ro','markersize',8)
            
            ixoutlier=find(myvar>=nthr);
            kxoutlier=length(ixoutlier);
            if kxoutlier>0
                REFcandidate(kcandidate+1:kcandidate+kxoutlier)=REF(ixoutlier);
                REFnearby(kcandidate+1:kcandidate+kxoutlier)=nearby(ixoutlier);
                tmpreason = cell(length(ixoutlier),1);
                tmpreason(:) = myvarname;
                REFreason(kcandidate+1:kcandidate+kxoutlier)=tmpreason;
                kcandidate=kcandidate+kxoutlier;
                semilogx(xvalue(ixoutlier),myvar(ixoutlier),'gd','markersize',12);
            end
            set(gca,'linewidth',1,'fontsize',12)
            xlabel('ngood','fontsize',16)
            ylabel(myvarname,'fontsize',16)
            title([CCDstring ': ' myvarname '>= ' num2str(nthr)],'fontsize',16)
            
            xvalue=stdmag;
        %   xvaluev=stdmagv;
            subplot(2,1,2)
            plot(xvalue,myvar,'.')
            hold on
        %   plot(xvaluev,myvarv,'ro','markersize',8)
            xlim([6 15])
            if kxoutlier>0
                plot(xvalue(ixoutlier),myvar(ixoutlier),'gd','markersize',12);
            end
            set(gca,'linewidth',1,'fontsize',12)
            xlabel(xmaglabel,'fontsize',16)
            ylabel(myvarname,'fontsize',16)
            if (doPlots == 1)
              saveas(gcf,[CCDstring '_' myvarname fqualifier '.eps'],'epsc2')
            endif
            close
        end
    end
    
    % --- 5.3 ndip3
    myvar=ndip3;                             % --------
%   myvarv=ndip3v;                           % --------
    myvarname='ndip3';                       % --------
    ix=find(ndip3<mean(ndip3)+20*std(ndip3));
    if (length(ix) > 0)
    nthr=max(2,ceil(mean(ndip3(ix))+10*std(ndip3(ix))));  % threshold of n             % --------
    
    if (length(nthr) > 0)
        xvalue=ngood;
    %   xvaluev=ngoodv;
        subplot(2,1,1)
        semilogx(xvalue,myvar,'.')
        hold on
    %   semilogx(xvaluev,myvarv,'ro','markersize',8)
        
        ixoutlier=find(myvar>=nthr);
        kxoutlier=length(ixoutlier);
        if kxoutlier>0
            REFcandidate(kcandidate+1:kcandidate+kxoutlier)=REF(ixoutlier);
            REFnearby(kcandidate+1:kcandidate+kxoutlier)=nearby(ixoutlier);
            tmpreason = cell(length(ixoutlier),1);
            tmpreason(:) = myvarname;
            REFreason(kcandidate+1:kcandidate+kxoutlier)=tmpreason;
            kcandidate=kcandidate+kxoutlier;
            semilogx(xvalue(ixoutlier),myvar(ixoutlier),'gd','markersize',12);
        end
        set(gca,'linewidth',1,'fontsize',12)
        xlabel('ngood','fontsize',16)
        ylabel(myvarname,'fontsize',16)
        title([CCDstring ': ' myvarname '>= ' num2str(nthr)],'fontsize',16)
        
        xvalue=stdmag;
    %   xvaluev=stdmagv;
        subplot(2,1,2)
        plot(xvalue,myvar,'.')
        hold on
    %   plot(xvaluev,myvarv,'ro','markersize',8)
        xlim([6 15])
        if kxoutlier>0
            plot(xvalue(ixoutlier),myvar(ixoutlier),'gd','markersize',12);
        end
        set(gca,'linewidth',1,'fontsize',12)
        xlabel(xmaglabel,'fontsize',16)
        ylabel(myvarname,'fontsize',16)
        if (doPlots == 1)
          saveas(gcf,[CCDstring '_' myvarname fqualifier '.eps'],'epsc2')
        endif
        close
    end
    end
    % --- 5.4 ndip4
    myvar=ndip4;                             % --------
%   myvarv=ndip4v;                           % --------
    myvarname='ndip4';                       % --------
    ix=find(ndip4<mean(ndip4)+20*std(ndip4));
    if (length(ix) > 0)
    nthr=max(3,ceil(mean(ndip4(ix))+10*std(ndip4(ix))));  % threshold of n             % --------
    
    if (length(nthr) > 0)
        xvalue=ngood;
    %   xvaluev=ngoodv;
        subplot(2,1,1)
        semilogx(xvalue,myvar,'.')
        hold on
    %   semilogx(xvaluev,myvarv,'ro','markersize',8)
        
        ixoutlier=find(myvar>=nthr);
        kxoutlier=length(ixoutlier);
        if kxoutlier>0
            REFcandidate(kcandidate+1:kcandidate+kxoutlier)=REF(ixoutlier);
            REFnearby(kcandidate+1:kcandidate+kxoutlier)=nearby(ixoutlier);
            tmpreason = cell(length(ixoutlier),1);
            tmpreason(:) = myvarname;
            REFreason(kcandidate+1:kcandidate+kxoutlier)=tmpreason;
            kcandidate=kcandidate+kxoutlier;
            semilogx(xvalue(ixoutlier),myvar(ixoutlier),'gd','markersize',12);
        end
        set(gca,'linewidth',1,'fontsize',12)
        xlabel('ngood','fontsize',16)
        ylabel(myvarname,'fontsize',16)
        title([CCDstring ': ' myvarname '>= ' num2str(nthr)],'fontsize',16)
        
        xvalue=stdmag;
    %   xvaluev=stdmagv;
        subplot(2,1,2)
        plot(xvalue,myvar,'.')
        hold on
    %   plot(xvaluev,myvarv,'ro','markersize',8)
        xlim([6 15])
        if kxoutlier>0
            plot(xvalue(ixoutlier),myvar(ixoutlier),'gd','markersize',12);
        end
        set(gca,'linewidth',1,'fontsize',12)
        xlabel(xmaglabel,'fontsize',16)
        ylabel(myvarname,'fontsize',16)
        if (doPlots == 1)
          saveas(gcf,[CCDstring '_' myvarname fqualifier '.eps'],'epsc2')
        endif
        close
    end
    end
    %if sum(kxoutlier)>0
    %    myvarname
    %    REFcandidate
    %end

    if (kcandidate > 0) 
      newREFreason = cell(length(REFcandidate),1);
      for (i = 1:length(REFcandidate))
        ixval=strcmp(REFcandidate,char(REFcandidate(i)));
        ix = find(ixval == 1);
        tmpreason = cell(length(ix),1);
        tmpreason(:) = char(REFreason(i));
        for (iy = 1:length(ix))
          if (index(char(newREFreason(ix(iy))),char(REFreason(i))) == 0) 
            if (length(char(newREFreason(ix(iy))) > 0))
               newREFreason(ix(iy)) = strcat(char(newREFreason(ix(iy))),";");
            end
            newREFreason(ix(iy)) = strcat(char(newREFreason(ix(iy))),char(REFreason(i)));
          end
        end
      
      end
  
  
      for (i = 1:length(REFcandidate))
          if (strcmp("X",REFnearby(i))) 
             newREFreason(i) = [char(REFcandidate(i))  " () " char(newREFreason(i))];
          else
             newREFreason(i) = [char(REFcandidate(i))  " (" char(REFnearby(i)) ") " char(newREFreason(i))];
          endif
      end
  
            % ttt = yyy(-1);
            % source /dasch/Pipeline/galactic_variable_search.m
      
      % --- remove duplicated entries and get the final output REF list
      eval(['cd ' outputdir2])
        myREF = unique(newREFreason);
        for i=1:length(myREF)
            fid = fopen([CCDstring 'candidate_REF' fqualifier '.txt'],'a');
            tmpREF = char(myREF{i});
            if (index(tmpREF,"APASS") == 1)
              signval = substr(tmpREF,13,1);
              if (signval == "1") 
                  sign = "+";
              else
                  sign = "-";
              endif
              tmpREF = ["APASS_J" substr(char(tmpREF),6,6) "." substr(char(tmpREF),12,1) sign substr(char(tmpREF),14,length(tmpREF)-13)];

            endif

            fprintf(fid,'%s \n',tmpREF);
            fclose(fid);
        end
    end    

end
clear adjacentburstdip0 adjacentburstdip20 adjacentburstdip30 b1 b2 b3 b4 CCDstring clip_median_local0 clip_ngood0 
clear clip_rms_local0 color0 damonBlueRms0 Damon_factor0 dec DECCCD DECchip DECchipa DECchipb DECchipc DECchipd declination0 
clear DecPM0 decrms0 dmagcatalog0 dN drad0 dradB0 error_bar_factor0 gscclass0 h i ii ix ixa ixb ixc ixd ixf ixoutlier ixval 
clear ixvalue iy iz k1 k2 k3 k4 kcandidate kccd keplerField0 kxoutlier lcrms1 lcrms2 lcrms2a lcrms3 lcrms3a lcrms4 lcrms4a 
clear lcrms5 lcrms5a lightcurverms10 lightcurverms20 lightcurverms30 lightcurverms40 lightcurverms50 MAGFlag0 magvsdeccorr0 
clear magvslimitingcorr0 magvsracorr0 Malmquist_factor0 Malmquist_factorB0 max_drad0 max_dradB0 max_local0 max_local20 
clear medDamonBlue0 meddecrms median_iso0 median_local0 medmagdec medmaglim medmagra medmf medmfB medNonDamonBlue0 medrarms 
clear medvalue medxp min_local0 min_local20 myREF myvar myvarname nadj1 nadj2 nadj3 nblend0 nburst0 nburst1 nburst2 nburst20 
clear nburst3 nburst30 nburst4 nburst40 nDamonBlue0 ndev2 ndev20 ndev3 ndev30 ndip0 ndip1 ndip2 ndip20 ndip3 ndip30 ndip4 
clear ndip40 newREFreason ngood ngood0 ngoodB ngoodB0 nNonDamonBlue0 nonDamonBluerms0 npoints0 nsig nslope_600 nslope_all0 
clear nthr nthrlow prefix ra ra0 RACCD RAchip RAchipa RAchipb RAchipc RAchipd range1 range2 range_iso0 range_local0 
clear range_local20 RaPM0 rarms0 REF REF0 REFcandidate REFmag REFreason rmsdradrms20 rms_factor0 rmslocal rms_local0 rmsthr 
clear Sextractor_Blend0 sigdecrms sigmagdec sigmaglim sigmagra sigmf sigmfB sign signval sigrarms sigxp slope_600 
clear slope_60_err0 slope_all0 slope_all_err0 stdmag Stdmag0 stdvalue stdvalue1 stdvalue2 tmpreason tmpREF versionId0 VFlag0 
clear x1 x2 xp xvalue y1 y2 yrbegin0 yrend0 lon0 lat0
clear nearbyObjects0 nearby REFnearby nearbymag
%who 
end % iteration loop
endTime = time() - startTime;
printf("Done at %d seconds\n",endTime);        
exit

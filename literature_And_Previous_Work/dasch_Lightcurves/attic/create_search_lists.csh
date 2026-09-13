# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

# create search lists in the proper format 
# Originally written by Silas Laycock
# Mar 10, 2008 Edward J. Los - used to assign GSC2.3.2 ID's to the set of stars. 
# Mar 24, 2008 Edward J. Los - assign the GSC RA and DEC instead of the original RA and DEC 
index -mb -n m44gsc.db ra%0.01 dec
index -mb -n CV.db ra%0.01 dec
index -mb -n ROSAT.db ra%0.01 dec
index -mb -n Ritter.db ra%0.01 dec
index -mb -n XMM.db ra%0.01 dec
index -mb -n XMMsimbad.db ra%0.01 dec
index -mb -n VSX.db ra%0.01 dec
index -mb -n M44.db ra%0.01 dec
index -mb -n ConstantStars.db ra%0.01 dec

search -j m44gsc.db -S2dd ra dec "00:00:20" < CV.db | column -a dra ddec drad | compute 'ddec=3600*(dec_2-dec_1); dra=3600*((ra_2-ra_1)*(cos(((dec_2+dec_1)/2)/57.29577951))); drad=sqrt(dra^2+ddec^2)' | sorttable -n drad | sorttable -u src_name| column src_name ra_2 dec_2 REF | sed 's/ra_2/ra/g' | sed 's/dec_2/dec/g' >! CV_REF.db

search -j m44gsc.db -S2dd ra dec "00:00:20" < Ritter.db | column -a dra ddec drad |compute 'ddec=3600*(dec_2-dec_1); dra=3600*((ra_2-ra_1)*(cos(((dec_2+dec_1)/2)/57.29577951))); drad=sqrt(dra^2+ddec^2)' | sorttable -n drad | sorttable -u src_name| column src_name ra_2 dec_2 REF | sed 's/ra_2/ra/g' | sed 's/dec_2/dec/g'> ! Ritter_REF.db

search -j m44gsc.db -S2dd ra dec "00:00:20" < ROSAT.db | column -a dra ddec drad | compute 'ddec=3600*(dec_2-dec_1); dra=3600*((ra_2-ra_1)*(cos(((dec_2+dec_1)/2)/57.29577951))); drad=sqrt(dra^2+ddec^2)' | sorttable -n drad |sorttable -u src_name| column src_name ra_2 dec_2 REF | sed 's/ra_2/ra/g' | sed 's/dec_2/dec/g'> ! ROSAT_REF.db

search -j m44gsc.db -S2dd ra dec "00:00:20" < XMM.db | column -a dra ddec drad |compute 'ddec=3600*(dec_2-dec_1); dra=3600*((ra_2-ra_1)*(cos(((dec_2+dec_1)/2)/57.29577951))); drad=sqrt(dra^2+ddec^2)' | sorttable -n drad | sorttable -u src_name| column src_name ra_2 dec_2 REF | sed 's/ra_2/ra/g' | sed 's/dec_2/dec/g'> ! XMM_REF.db

search -j m44gsc.db -S2dd ra dec "00:00:20" < XMMsimbad.db | column -a dra ddec drad |compute 'ddec=3600*(dec_2-dec_1); dra=3600*((ra_2-ra_1)*(cos(((dec_2+dec_1)/2)/57.29577951))); drad=sqrt(dra^2+ddec^2)' | sorttable -n drad | sorttable -u src_name| column src_name ra_2 dec_2 REF | sed 's/ra_2/ra/g' | sed 's/dec_2/dec/g' > ! XMMsimbad_REF.db

search -j m44gsc.db -S2dd ra dec "00:00:30" < VSX.db | column -a dra ddec drad | compute 'ddec=3600*(dec_2-dec_1); dra=3600*((ra_2-ra_1)*(cos(((dec_2+dec_1)/2)/57.29577951))); drad=sqrt(dra^2+ddec^2)' | sorttable -n drad | sorttable -u src_name | column src_name ra_2 dec_2 REF | sed 's/ra_2/ra/g' | sed 's/dec_2/dec/g' > ! VSX_REF.db

search -j m44gsc.db -S2dd ra dec "00:00:30" < M44.db | column -a dra ddec drad | compute 'ddec=3600*(dec_2-dec_1); dra=3600*((ra_2-ra_1)*(cos(((dec_2+dec_1)/2)/57.29577951))); drad=sqrt(dra^2+ddec^2)' | sorttable -n drad | sorttable -u src_name | column src_name ra_2 dec_2 REF | sed 's/ra_2/ra/g' | sed 's/dec_2/dec/g' > ! M44_REF.db

search -j m44gsc.db -S2dd ra dec "00:00:20" < ConstantStars.db | column -a dra ddec drad | compute 'ddec=3600*(dec_2-dec_1); dra=3600*((ra_2-ra_1)*(cos(((dec_2+dec_1)/2)/57.29577951))); drad=sqrt(dra^2+ddec^2)' | sorttable -n drad | sorttable -u src_name | column src_name ra_2 dec_2 REF | sed 's/ra_2/ra/g' | sed 's/dec_2/dec/g' > ! ConstantStars_REF.db

#column -i /data/bernie1/Plates/Pipeline_1.2/catalogs/gsc2_60by60/gsc2_within_10degrees_of_m44_Blt15.list -a REF | compute 'REF=src_name' | column src_name REF ra dec > General_REF.db

exit


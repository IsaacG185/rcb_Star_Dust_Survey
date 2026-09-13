#!/bin/csh
# Copyright the President and Fellows of Harvard College.
# Licensed under the MIT License

#  Sep  6, 2010 Edward J. Los - Convert build.csh to build on han.cfa.harvard.edu for the external DASCH website
 set lib64 = ""
gcc -g -c -O2 -I/usr/include/mysql pipelineutils.c  -o pipelineutils.o
gcc -g -c -O2  -I/usr/include/mysql  -I/home/elos/wcstools photometryutils.c    -D_FILE_OFFSET_BITS=64 -o photometryutils.o
gcc -g -c -O2 -I/usr/include/mysql RA2fpix.c  -o RA2fpix.o
gcc -g -c -O2  urlencode.c  -o urlencode.o
gcc -pedantic -Wall -g  -fPIC -O3 -g -DUSE_LIST_NODE_ALLOCATOR   -c -o kdtree.o kdtree.c

gcc -g -c -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64  -I/home/elos/install/include -I/usr/include/mysql -I/home/elos/wcstools searchgsc.c -o searchgsc.o

gcc -g -O2 -DNOVERBOSE  -DRUN_TIMER -DLIMIT_SHIFT -DCROP_BORDER -DSCANNER -DLOS_DEBUG -DFIX_TNX   -I/usr/include/mysql  -I/home/elos/wcstools/libwcs -I/home/elos/wcstools/libwcs/libwcs  -c -o scanread.o scanread.c

ar csr pipelineutils.a pipelineutils.o photometryutils.o  RA2fpix.o urlencode.o scanread.o searchgsc.o kdtree.o

echo "showflags"
gcc -g -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64   -I/usr/include/mysql  showflags.c   pipelineutils.a  /home/elos/install/lib/libwcs.a -lm -L/usr/lib64/mysql -lmysqlclient -o  showflags
cp showflags /proj/dasch/bin 

echo "votable"
gcc -g -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64  -I/usr/include/mysql -I/home/elos/wcstools/libwcs -I/home/elos/install/include  votable.c     pipelineutils.a  /home/elos/install/lib/libtable.a /home/elos/install/lib/libutil.a  /home/elos/install/lib/libwcs.a    -lm -L/usr/lib64/mysql -lmysqlclient -o  votable 
cp votable /proj/dasch/bin 

echo "web_query"
gcc -g -O2 -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -I/usr/include/mysql -I/home/elos/wcstools  web_query.c  pipelineutils.a  /home/elos/install/lib/libwcs.a  -lm  -lrt -lnsl -L/usr/lib64/mysql -lmysqlclient -o  web_query  
cp web_query /proj/dasch/bin


echo "web_point"
gcc -g -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64  -I/usr/include/mysql -I/home/elos/wcstools   -I/home/elos/install/include   web_point.c    pipelineutils.a  /home/elos/install/lib/libtable.a /home/elos/install/lib/libutil.a /home/elos/install/lib/libwcs.a   -lrt -lnsl  -lm  -L/usr/lib64/mysql -lmysqlclient -o  web_point
cp web_point /proj/dasch/bin

echo "dasch_scat"
gcc -g -O0   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64  -I/usr/include/mysql -I/home/elos/wcstools   -I/home/elos/install/include   dasch_scat.c     pipelineutils.a   /home/elos/install/lib/libwcs.a   -lm  -L/usr/lib64/mysql -lmysqlclient -lnsl -o  dasch_scat
cp dasch_scat /proj/dasch/bin


echo "web_plot"
gcc -g -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64  -I/usr/include/mysql -I/home/elos/wcstools    -I/home/elos/install/include web_plot.c     pipelineutils.a /home/elos/install/lib/libwcs.a   /home/elos/install/lib/libgd.a /home/elos/install/lib/libtable.a /home/elos/install/lib/libutil.a -lm   -lrt -lnsl  -L/usr/lib64/mysql -lmysqlclient -o  web_plot
cp web_plot /proj/dasch/bin

echo "web_authorize"
gcc -O2 -g -o web_authorize  -I/usr/include/mysql  -I/home/elos/wcstools    -I/home/elos/install/include web_authorize.c pipelineutils.a  /home/elos/install/lib/libwcs.a  /home/elos/install/lib/libutil.a  -L/usr/lib64/mysql -lmysqlclient -ldl -lm 
cp web_authorize /proj/dasch/bin


echo "dumpbin"
gcc -g -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1  -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64  -I/usr/include/mysql -I/home/elos/wcstools  -I/home/elos/install/include  dumpbin.c    /home/elos/install/lib/libtable.a /home/elos/install/lib/libutil.a pipelineutils.a  /home/elos/install/lib/libwcs.a  -lrt -lnsl -lm   -L/usr/lib64/mysql -lmysqlclient -o  dumpbin
cp dumpbin /proj/dasch/bin

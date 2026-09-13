        Build SELinux tools    Mon Aug 23 14:28:29 EDT 2021
This file is /home/scanner/Pipeline/selinux/readme.txt

[scanner@localhost Pipeline]$ rpm -q -a | grep -i "kernel"
kernel-headers-3.10.0-1127.13.1.el7.x86_64
kernel-tools-libs-3.10.0-1127.13.1.el7.x86_64
kernel-devel-3.10.0-229.el7.x86_64
kernel-tools-3.10.0-1127.13.1.el7.x86_64
abrt-addon-kerneloops-2.1.11-57.el7.centos.x86_64
kernel-3.10.0-1127.13.1.el7.x86_64
texlive-l3kernel-svn29409.SVN_4469-45.el7.noarch
kernel-3.10.0-229.el7.x86_64
kernel-devel-3.10.0-1127.13.1.el7.x86_64


/usr/src/kernels/3.10.0-1127.13.1.el7.x86_64
/usr/src/kernels/3.10.0-229.el7.x86_64

/usr/src/kernels/3.10.0-229.el7.x86_64/scripts/selinux/genheaders/genheaders.c
/usr/src/kernels/3.10.0-229.el7.x86_64/scripts/selinux/mdp/mdp.c


gcc

find /usr/src/kernels -type f  -exec grep -i "classmap" {} \; -print

/usr/src/kernels/3.10.0-1127.13.1.el7.x86_64/Kbuild/scripts/selinux/genheaders/genheaders.c
#include "classmap.h"
/usr/src/kernels/3.10.0-1127.13.1.el7.x86_64/Kbuild/scripts/selinux/mdp/mdp.c
#include "classmap.h"
/usr/src/kernels/3.10.0-1127.13.1.el7.x86_64/Kbuild/scripts/modinst/selinux/genheaders/genheaders.c
#include "classmap.h"
/usr/src/kernels/3.10.0-1127.13.1.el7.x86_64/Kbuild/scripts/modinst/selinux/mdp/mdp.c
$(obj)/flask.h: $(src)/include/classmap.h FORCE
/usr/src/kernels/3.10.0-1127.13.1.el7.x86_64/Kbuild/security/selinux/Makefile


[root@localhost kernels]# find /usr/src/kernels/ -type f -exec grep "classmap" {} \; -print
#include "classmap.h"
/usr/src/kernels/3.10.0-229.el7.x86_64/scripts/selinux/genheaders/genheaders.c
#include "classmap.h"
/usr/src/kernels/3.10.0-229.el7.x86_64/scripts/selinux/mdp/mdp.c
$(obj)/flask.h: $(src)/include/classmap.h FORCE
/usr/src/kernels/3.10.0-229.el7.x86_64/security/selinux/Makefile
#include "classmap.h"
/usr/src/kernels/3.10.0-1127.13.1.el7.x86_64/scripts/selinux/genheaders/genheaders.c
#include "classmap.h"
/usr/src/kernels/3.10.0-1127.13.1.el7.x86_64/scripts/selinux/mdp/mdp.c
#include "classmap.h"
/usr/src/kernels/3.10.0-1127.13.1.el7.x86_64/scripts/modinst/selinux/genheaders/genheaders.c
#include "classmap.h"
/usr/src/kernels/3.10.0-1127.13.1.el7.x86_64/scripts/modinst/selinux/mdp/mdp.c
$(obj)/flask.h: $(src)/include/classmap.h FORCE
/usr/src/kernels/3.10.0-1127.13.1.el7.x86_64/security/selinux/Makefile
#include "classmap.h"
/usr/src/kernels/3.10.0-1127.13.1.el7.x86_64/Kbuild/scripts/selinux/genheaders/genheaders.c
#include "classmap.h"
/usr/src/kernels/3.10.0-1127.13.1.el7.x86_64/Kbuild/scripts/selinux/mdp/mdp.c
#include "classmap.h"
/usr/src/kernels/3.10.0-1127.13.1.el7.x86_64/Kbuild/scripts/modinst/selinux/genheaders/genheaders.c
#include "classmap.h"
/usr/src/kernels/3.10.0-1127.13.1.el7.x86_64/Kbuild/scripts/modinst/selinux/mdp/mdp.c
$(obj)/flask.h: $(src)/include/classmap.h FORCE
/usr/src/kernels/3.10.0-1127.13.1.el7.x86_64/Kbuild/security/selinux/Makefile
[root@localhost kernels]# find . -type f -exec grep "flask" {} \; -print
Binary file /usr/src/kernels/3.10.0-229.el7.x86_64/scripts/selinux/genheaders/genheaders matches
/usr/src/kernels/3.10.0-229.el7.x86_64/scripts/selinux/genheaders/genheaders
	printf("usage: %s flask.h av_permissions.h\n", progname);
/usr/src/kernels/3.10.0-229.el7.x86_64/scripts/selinux/genheaders/genheaders.c
$(addprefix $(obj)/,$(selinux-y)): $(obj)/flask.h
quiet_cmd_flask = GEN     $(obj)/flask.h $(obj)/av_permissions.h
      cmd_flask = scripts/selinux/genheaders/genheaders $(obj)/flask.h $(obj)/av_permissions.h
targets += flask.h av_permissions.h
$(obj)/flask.h: $(src)/include/classmap.h FORCE
	$(call if_changed,flask)
/usr/src/kernels/3.10.0-229.el7.x86_64/security/selinux/Makefile
Binary file /usr/src/kernels/3.10.0-1127.13.1.el7.x86_64/scripts/selinux/genheaders/genheaders matches
/usr/src/kernels/3.10.0-1127.13.1.el7.x86_64/scripts/selinux/genheaders/genheaders
	printf("usage: %s flask.h av_permissions.h\n", progname);
/usr/src/kernels/3.10.0-1127.13.1.el7.x86_64/scripts/selinux/genheaders/genheaders.c
Binary file /usr/src/kernels/3.10.0-1127.13.1.el7.x86_64/scripts/modinst/selinux/genheaders/genheaders matches
/usr/src/kernels/3.10.0-1127.13.1.el7.x86_64/scripts/modinst/selinux/genheaders/genheaders
	printf("usage: %s flask.h av_permissions.h\n", progname);
/usr/src/kernels/3.10.0-1127.13.1.el7.x86_64/scripts/modinst/selinux/genheaders/genheaders.c
$(addprefix $(obj)/,$(selinux-y)): $(obj)/flask.h
quiet_cmd_flask = GEN     $(obj)/flask.h $(obj)/av_permissions.h
      cmd_flask = scripts/selinux/genheaders/genheaders $(obj)/flask.h $(obj)/av_permissions.h
targets += flask.h av_permissions.h
$(obj)/flask.h: $(src)/include/classmap.h FORCE
	$(call if_changed,flask)
/usr/src/kernels/3.10.0-1127.13.1.el7.x86_64/security/selinux/Makefile
Binary file /usr/src/kernels/3.10.0-1127.13.1.el7.x86_64/Kbuild/scripts/selinux/genheaders/genheaders matches
/usr/src/kernels/3.10.0-1127.13.1.el7.x86_64/Kbuild/scripts/selinux/genheaders/genheaders
	printf("usage: %s flask.h av_permissions.h\n", progname);
/usr/src/kernels/3.10.0-1127.13.1.el7.x86_64/Kbuild/scripts/selinux/genheaders/genheaders.c
Binary file /usr/src/kernels/3.10.0-1127.13.1.el7.x86_64/Kbuild/scripts/modinst/selinux/genheaders/genheaders matches
/usr/src/kernels/3.10.0-1127.13.1.el7.x86_64/Kbuild/scripts/modinst/selinux/genheaders/genheaders
	printf("usage: %s flask.h av_permissions.h\n", progname);
/usr/src/kernels/3.10.0-1127.13.1.el7.x86_64/Kbuild/scripts/modinst/selinux/genheaders/genheaders.c
$(addprefix $(obj)/,$(selinux-y)): $(obj)/flask.h
quiet_cmd_flask = GEN     $(obj)/flask.h $(obj)/av_permissions.h
      cmd_flask = scripts/selinux/genheaders/genheaders $(obj)/flask.h $(obj)/av_permissions.h
targets += flask.h av_permissions.h
$(obj)/flask.h: $(src)/include/classmap.h FORCE
	$(call if_changed,flask)
/usr/src/kernels/3.10.0-1127.13.1.el7.x86_64/Kbuild/security/selinux/Makefile
[root@localhost kernels]# 

[scanner@localhost selinux]$ ls -l  /usr/src/kernels/3.10.0-229.el7.x86_64/scripts/selinux/genheaders/genheaders 
-rwxr-xr-x. 1 root root 24280 Mar  6  2015 /usr/src/kernels/3.10.0-229.el7.x86_64/scripts/selinux/genheaders/genheaders
[scanner@localhost selinux]$ ls -l  /usr/src/kernels/3.10.0-1127.13.1.el7.x86_64/scripts/selinux/genheaders/genheaders 
-rwxr-xr-x. 1 root root 25368 Jun 23  2020 /usr/src/kernels/3.10.0-1127.13.1.el7.x86_64/scripts/selinux/genheaders/genheaders
[scanner@localhost selinux]$ ls -l  /usr/src/kernels/3.10.0-1127.13.1.el7.x86_64/scripts/modinst/selinux/genheaders/genheaders 
-rwxr-xr-x. 1 root root 25368 Aug  3 21:41 /usr/src/kernels/3.10.0-1127.13.1.el7.x86_64/scripts/modinst/selinux/genheaders/genheaders
[scanner@localhost selinux]$ ls -l  /usr/src/kernels/3.10.0-1127.13.1.el7.x86_64/Kbuild/scripts/selinux/genheaders/genheaders 
-rwxr-xr-x. 1 root root 25368 Aug  3 21:44 /usr/src/kernels/3.10.0-1127.13.1.el7.x86_64/Kbuild/scripts/selinux/genheaders/genheaders
[scanner@localhost selinux]$ ls -l  /usr/src/kernels/3.10.0-1127.13.1.el7.x86_64/Kbuild/scripts/modinst/selinux/genheaders/genheaders 
-rwxr-xr-x. 1 root root 25368 Aug  3 21:44 /usr/src/kernels/3.10.0-1127.13.1.el7.x86_64/Kbuild/scripts/modinst/selinux/genheaders/genheaders
[scanner@localhost selinux]$ 
[scanner@localhost selinux]$ 
[scanner@localhost selinux]$ md5sum  /usr/src/kernels/3.10.0-229.el7.x86_64/scripts/selinux/genheaders/genheaders 
8bff49fed897f7849294746102dd8631  /usr/src/kernels/3.10.0-229.el7.x86_64/scripts/selinux/genheaders/genheaders
[scanner@localhost selinux]$ md5sum  /usr/src/kernels/3.10.0-1127.13.1.el7.x86_64/scripts/selinux/genheaders/genheaders 
bf73f9ab9991c0cf79c881c796a7caaf  /usr/src/kernels/3.10.0-1127.13.1.el7.x86_64/scripts/selinux/genheaders/genheaders
[scanner@localhost selinux]$ md5sum  /usr/src/kernels/3.10.0-1127.13.1.el7.x86_64/scripts/modinst/selinux/genheaders/genheaders 
bf73f9ab9991c0cf79c881c796a7caaf  /usr/src/kernels/3.10.0-1127.13.1.el7.x86_64/scripts/modinst/selinux/genheaders/genheaders
[scanner@localhost selinux]$ md5sum  /usr/src/kernels/3.10.0-1127.13.1.el7.x86_64/Kbuild/scripts/selinux/genheaders/genheaders 
bf73f9ab9991c0cf79c881c796a7caaf  /usr/src/kernels/3.10.0-1127.13.1.el7.x86_64/Kbuild/scripts/selinux/genheaders/genheaders
[scanner@localhost selinux]$ md5sum  /usr/src/kernels/3.10.0-1127.13.1.el7.x86_64/Kbuild/scripts/modinst/selinux/genheaders/genheaders 
bf73f9ab9991c0cf79c881c796a7caaf  /usr/src/kernels/3.10.0-1127.13.1.el7.x86_64/Kbuild/scripts/modinst/selinux/genheaders/genheaders
[scanner@localhost selinux]$ 


cd /home/scanner/Pipeline/selinux/
[root@localhost selinux]# cp /usr/src/kernels/3.10.0-1127.13.1.el7.x86_64/Kbuild/security/selinux/Makefile ./Makefile2
[root@localhost selinux]# chgrp scanner Makefile2
[root@localhost selinux]# chown scanner Makefile2
[root@localhost selinux]# src="/home/scanner/Pipeline/selinux"
[root@localhost selinux]# echo "$src"
/home/scanner/Pipeline/selinux
cp /usr/src/kernels/3.10.0-1127.13.1.el7.x86_64/Kbuild/scripts/modinst/selinux/genheaders/genheaders /home/scanner/Pipeline/selinux/
cd /home/scanner/Pipeline/selinux/
chgrp scanner genheaders
chown scanner genheaders
./genheaders
usage: ./genheaders flask.h av_permissions.h
./genheaders flask.h av_permissions.h
exit
[root@localhost selinux]# exit
exit

ERROR: nothing in ./include

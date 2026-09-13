// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/*
 *  testtimer.c - test process kill ability
 *  gcc testtimer.c -g -o testtimer
 *
 */

#include <stdarg.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <signal.h>
#include <errno.h>
#include <unistd.h>
int main(int argc, char *argv[])
{ 

  struct timespec timeout;
  sigset_t waitset;
  siginfo_t waitinfo;
  time_t startTime;
  time_t curTime;
  int result;
  pid_t pgid;
  pid_t pid;
  pid_t ppid;

  time(&startTime);
  if (argc < 2) {
    printf("ERROR: there must be at least one argument\n");
    return(-1);
  }
  sigemptyset(&waitset);
  pgid = getpgid(0);
  pid = getpid();
  ppid = getppid();
#if 0
  sigaddset(&waitset,SIGKILL);
  sigaddset(&waitset,SIGSTOP);
  sigaddset(&waitset,SIGCONT);
  sigaddset(&waitset,SIGUSR1);
#endif
  timeout.tv_sec =  10;
  timeout.tv_nsec = 0;
  printf("testtimer beginning wait for %s\n",argv[1]);
  while(1) {

    result = sigtimedwait(&waitset,&waitinfo,&timeout);
    time(&curTime);
    curTime = curTime - startTime;
    fprintf(stderr,"testimer at %5d seconds pid %d, pgid %d, ppid %d for %s \n",curTime,pid,pgid,ppid,argv[1]);
    fflush(stderr);
#if 0
    printf("Timeout result %d, errno %d si_signo %d\n",result,errno,waitinfo.si_signo);
#endif
    if (result >= 0) {
      break;
    }

  }
  
  return(0);


}

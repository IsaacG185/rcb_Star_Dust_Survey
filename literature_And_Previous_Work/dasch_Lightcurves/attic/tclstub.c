// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* los.c
 *  Used when there are too many 3ware messages to read from mail
 *  gcc -g -c tclstub.c -o tclstub.o  
 */

#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
int TclFormatInt()
{
  printf("Illegal call to TclFormatInt\n");
  exit(-1);

}

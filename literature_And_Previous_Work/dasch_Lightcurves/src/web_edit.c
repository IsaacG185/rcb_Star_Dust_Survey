// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* web_edit.c
 *
 * This program is executed from the editstacks.php script on daschprivate to allow editing of the stacks
 *
 *  NOTE: the stackLocation format is unknown
 *  NOTE: boxnumbr is not supported
 *

 cc -ggdb -O2   -DLinux -DSTDC_HEADERS=1 -DCONFIG_BROKETS -DHAVE_DIRENT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_MEMCPY=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_FCNTL_H=1 -DHAVE_ALLOCA_H=1 -DHAVE_FLOAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_MEMCHR=1 -DHAVE_STRCHR=1 -DHAVE_STRERROR=1 -DHAVE_STRUCT_TM_TM_ZONE=1 -DHAVE_TM_ZONE=1 -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64  -I/usr/include/mysql -I/dasch/install/include web_edit.c pipelineutils.a -L /dasch/install/lib -lm  -L/usr/lib${lib64}/mysql  -lmysqlclient -ldl -pthread    -ltable -lutil  -lwcs plottransientstub.o -o web_edit
cp web_edit /dasch/install/bin

 *
 *
 * Jul  4, 2019 Edward J. Los   Initial version
 * Jul 30, 2019 Edward J. Los   Add '-t' editConditionVersionId and '-u' editEventsVersionId to handle deletion and undeletion
 * Aug 18, 2019 Edward J. Los   "Jacket JPEG","Plate JPEG","Flood", and "Scanned" can not be changed.
 *                              Validate remoteaddr
 * Aug 29, 2019 Edward J. Los   Handle plate number with leading zeros
 *
 *
 *
 *
 */


#include <math.h>
#include <time.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <ctype.h>

#include <table.h>

#include <mysql.h>

#include <libwcs/fitsfile.h>
#include <libwcs/wcs.h>

#include "scandb.h"
#include "pipelineutils.h"
#include "galaxyutils.h"
#include "photometryutils.h"
#include "daschunistd.h"

#define MAX_REMOTE_ADDR 25
#define MAX_MD5_LENGTH 26

extern char *catalogText[MAX_CATALOG_NUMBER];
char *urlencode(char *string, char *outstring, int outlength);


char *md5sumTable[] = {
  "g37FdU9QPPqu4JKf1Il05w",
  "cHj/sKxhtEWBHPMJYyMDAA",
  "VdOkd2/Wa0xwsQaVil0FSg",
  "shWLdovOCrk+VL5yW+ONfA",
  "Od8NGzeAOMyitU9pjshKbg",
  "+vITU31lDrc5WRVp6twJ9A",
  "73ibpQGG0bUSXJmFoDOEiw",
  "Sj+H22QqZXYC4pbfA/yE/w",
  "NSXU9dbbrvb6oAvll4MLgQ",
  "4o5Iq6zGF/akUbAPk32eZg",
};
int md5sumTableSize = sizeof(md5sumTable)/sizeof(char *);


/*  externally sourced md5sum algorithm.  Also used in CheckJPEG and CheckBacula  */
typedef struct MD5state
{
  unsigned int len;
  unsigned int state[4];
}MD5state;

typedef struct Table
{
  unsigned int	sin;	/* integer part of 4294967296 times abs(sin(i)) */
  unsigned char	x;	/* index into data block */
  unsigned char	rot;	/* amount to rotate left by */
} Table;

/*
 *	Rotate amounts used in the algorithm
 */
enum
  {
    S11=	7,
    S12=	12,
    S13=	17,
    S14=	22,

    S21=	5,
    S22=	9,
    S23=	14,
    S24=	20,

    S31=	4,
    S32=	11,
    S33=	16,
    S34=	23,

    S41=	6,
    S42=	10,
    S43=	15,
    S44=	21
  };


Table tab[] =
  {
    /* round 1 */
    { 0xd76aa478, 0, S11},
    { 0xe8c7b756, 1, S12},
    { 0x242070db, 2, S13},
    { 0xc1bdceee, 3, S14},
    { 0xf57c0faf, 4, S11},
    { 0x4787c62a, 5, S12},
    { 0xa8304613, 6, S13},
    { 0xfd469501, 7, S14},
    { 0x698098d8, 8, S11},
    { 0x8b44f7af, 9, S12},
    { 0xffff5bb1, 10, S13},
    { 0x895cd7be, 11, S14},
    { 0x6b901122, 12, S11},
    { 0xfd987193, 13, S12},
    { 0xa679438e, 14, S13},
    { 0x49b40821, 15, S14},

    /* round 2 */
    { 0xf61e2562, 1, S21},
    { 0xc040b340, 6, S22},
    { 0x265e5a51, 11, S23},
    { 0xe9b6c7aa, 0, S24},
    { 0xd62f105d, 5, S21},
    {  0x2441453, 10, S22},
    { 0xd8a1e681, 15, S23},
    { 0xe7d3fbc8, 4, S24},
    { 0x21e1cde6, 9, S21},
    { 0xc33707d6, 14, S22},
    { 0xf4d50d87, 3, S23},
    { 0x455a14ed, 8, S24},
    { 0xa9e3e905, 13, S21},
    { 0xfcefa3f8, 2, S22},
    { 0x676f02d9, 7, S23},
    { 0x8d2a4c8a, 12, S24},

    /* round 3 */
    { 0xfffa3942, 5, S31},
    { 0x8771f681, 8, S32},
    { 0x6d9d6122, 11, S33},
    { 0xfde5380c, 14, S34},
    { 0xa4beea44, 1, S31},
    { 0x4bdecfa9, 4, S32},
    { 0xf6bb4b60, 7, S33},
    { 0xbebfbc70, 10, S34},
    { 0x289b7ec6, 13, S31},
    { 0xeaa127fa, 0, S32},
    { 0xd4ef3085, 3, S33},
    {  0x4881d05, 6, S34},
    { 0xd9d4d039, 9, S31},
    { 0xe6db99e5, 12, S32},
    { 0x1fa27cf8, 15, S33},
    { 0xc4ac5665, 2, S34},

    /* round 4 */
    { 0xf4292244, 0, S41},
    { 0x432aff97, 7, S42},
    { 0xab9423a7, 14, S43},
    { 0xfc93a039, 5, S44},
    { 0x655b59c3, 12, S41},
    { 0x8f0ccc92, 3, S42},
    { 0xffeff47d, 10, S43},
    { 0x85845dd1, 1, S44},
    { 0x6fa87e4f, 8, S41},
    { 0xfe2ce6e0, 15, S42},
    { 0xa3014314, 6, S43},
    { 0x4e0811a1, 13, S44},
    { 0xf7537e82, 4, S41},
    { 0xbd3af235, 11, S42},
    { 0x2ad7d2bb, 2, S43},
    { 0xeb86d391, 9, S44},
  };
static unsigned char t64d[256];
static char t64e[64];

void init64(void)
{
  int c, i;

  memset(t64d, 255, 256);
  memset(t64e, '=', 64);
  i = 0;
  for(c = 'A'; c <= 'Z'; c++){
    t64e[i] = c;
    t64d[c] = i++;
  }
  for(c = 'a'; c <= 'z'; c++){
    t64e[i] = c;
    t64d[c] = i++;
  }
  for(c = '0'; c <= '9'; c++){
    t64e[i] = c;
    t64d[c] = i++;
  }
  t64e[i] = '+';
  t64d['+'] = i++;
  t64e[i] = '/';
  t64d['/'] = i;
}

/*
 *	input (unsigned int) into output (byte). Assumes len is
 *	a multiple of 4.
 */
void encode(unsigned char *output, unsigned int *input, unsigned int len)
{
  unsigned int x;
  unsigned char *e;

  for(e = output + len; output < e;) {
    x = *input++;
    *output++ = x;
    *output++ = x >> 8;
    *output++ = x >> 16;
    *output++ = x >> 24;
  }
}

/*
 *	decodes input (unsigned char) into output (unsigned int). Assumes len is
 *	a multiple of 4.
 */
void decode(unsigned int *output, unsigned char *input, unsigned int len)
{
  unsigned char *e;

  for(e = input+len; input < e; input += 4)
    *output++ = input[0] | (input[1] << 8) |
      (input[2] << 16) | (input[3] << 24);
}

/*
 *
 *  I require len to be a multiple of 64 for all but
 *  the last call
 */
MD5state* md5(unsigned char *p, unsigned int len, unsigned char *digest, MD5state *s)
{
  unsigned int a, b, c, d, tmp;
  unsigned int i, done;
  Table *t;
  unsigned char *end;
  unsigned int x[16];

  if(s == NULL){
    s = (MD5state *)calloc(sizeof(*s),1);
    if(s == NULL)
      return NULL;

    /* seed the state, these constants would look nicer big-endian */
    s->state[0] = 0x67452301;
    s->state[1] = 0xefcdab89;
    s->state[2] = 0x98badcfe;
    s->state[3] = 0x10325476;
  }
  s->len += len;

  i = len & 0x3f;
  if(i || len == 0){
    done = 1;

    /* pad the input, assume there's room */
    if(i < 56)
      i = 56 - i;
    else
      i = 120 - i;
    if(i > 0){
      memset(p + len, 0, i);
      p[len] = 0x80;
    }
    len += i;

    /* append the count */
    x[0] = s->len<<3;
    x[1] = s->len>>29;
    encode(p+len, x, 8);
  } else
    done = 0;

  for(end = p+len; p < end; p += 64){
    a = s->state[0];
    b = s->state[1];
    c = s->state[2];
    d = s->state[3];

    decode(x, p, 64);

    for(i = 0; i < 64; i++){
      t = tab + i;
      switch(i>>4){
      case 0:
        a += (b & c) | (~b & d);
        break;
      case 1:
        a += (b & d) | (c & ~d);
        break;
      case 2:
        a += b ^ c ^ d;
        break;
      case 3:
        a += c ^ (b | ~d);
        break;
      }
      a += x[t->x] + t->sin;
      a = (a << t->rot) | (a >> (32 - t->rot));
      a += b;

      /* rotate variables */
      tmp = d;
      d = c;
      c = b;
      b = a;
      a = tmp;
    }

    s->state[0] += a;
    s->state[1] += b;
    s->state[2] += c;
    s->state[3] += d;
  }

  /* return result */
  if(done){
    encode(digest, s->state, 16);
    free(s);
    return NULL;
  }
  return s;
}


int enc64(char *out, unsigned char *in, int n)
{
  int i;
  unsigned long b24;
  char *start = out;

  if(t64e[0] == 0)
    init64();
  for(i = n/3; i > 0; i--){
    b24 = (*in++)<<16;
    b24 |= (*in++)<<8;
    b24 |= *in++;
    *out++ = t64e[(b24>>18)];
    *out++ = t64e[(b24>>12)&0x3f];
    *out++ = t64e[(b24>>6)&0x3f];
    *out++ = t64e[(b24)&0x3f];
  }

  switch(n%3){
  case 2:
    b24 = (*in++)<<16;
    b24 |= (*in)<<8;
    *out++ = t64e[(b24>>18)];
    *out++ = t64e[(b24>>12)&0x3f];
    *out++ = t64e[(b24>>6)&0x3f];
    break;
  case 1:
    b24 = (*in)<<16;
    *out++ = t64e[(b24>>18)];
    *out++ = t64e[(b24>>12)&0x3f];
    *out++ = '=';
    break;
  }
  *out++ = '=';
  *out = 0;
  return out - start;
}


int main(int argc,char *argv[])
{
  char *argstr;
  char cmdchar;
  int errorFlag = 0;
  int plateNumber = 0;
  int newPlateNumber = 0;
  int insertPlateNumber = 0;
  int nvals;
  char series[MAX_SERIES_STRING];
  char newseries[MAX_SERIES_STRING];
  char insertseries[MAX_SERIES_STRING];
  char editdate[MAX_DATE_STRING];
  char newplate[MAX_PLATE_NAME];
  char tempplate[MAX_PLATE_NAME];
  char editevent[MAX_STACKS_BUFFER];

  int newvalidcondition = -1; /* -1 = unknown; 0 = "no" 1 = "yes" */
  int newpickedcomment = 0;
  char newpickednotes[MAX_STACKS_BUFFER];

  int newvalidevent = -1; /* -1 = unknown; 0 = "no" 1 = "yes" */
  int newpickedstatus = 0;
  char newstacklocation[MAX_STACKS_BUFFER];
  char neweventnotes[MAX_STACKS_BUFFER];
  char editorname[MAX_STACKS_BUFFER];
  char editcomment[MAX_STACKS_BUFFER];
  char logfile[MAX_STACKS_BUFFER];
  FILE *logHandle = NULL;
  char timestr[100];
  time_t startTime;
  struct tm *ptr;
  char *mysqlhost;
  MYSQL my_connection;
  MYSQL *pConnection = &my_connection;
  MYSQL my_phot_connection;
  MYSQL *pPhotConnection = &my_phot_connection;
  char *stacksusername;
  char *stackspassword;
  MYSQL my_stacks_connection;
  MYSQL *pStacksConnection = &my_stacks_connection;

  int condition_nrecs = 0;
  int conditionIndex;
  int editConditionVersionId = -1;
  int editConditionIndex = -1;
  int editEventsVersionId = -1;
  int editEventsIndex = -1;
  int curConditionVersionId = -1;
  int maxConditionVersionId = -1;
  int conditionChangedFlag = 0;
  int eventsChangedFlag = 0;
  int pickedStatusIndex;
  PPLATECONDITION pPlateCondition;
  PPLATECONDITION plateConditionTable = NULL;
  int events_nrecs = 0;
  int eventsIndex = 0;
  int curEventsVersionId = -1;
  int curEventsIndex = -1;
  int maxEventVersionId = -1;
  int nextVersionId = -1;
  PPLATEEVENTS pPlateEvents;
  PPLATEEVENTS plateEventsTable = NULL;

  int newcondition_nrecs = 0;
  int newconditionIndex;
  int maxNewConditionVersionId = -1;
  int curConditionIndex = -1;
  int pickedCommentIndex;
  PPLATECONDITION pPlateNewCondition;
  PPLATECONDITION plateNewConditionTable = NULL;
  int newevents_nrecs = 0;
  int neweventsIndex = 0;
  int maxNewEventVersionId = -1;
  PPLATEEVENTS pPlateNewEvents;
  PPLATEEVENTS plateNewEventsTable = NULL;

  char remoteaddr[MAX_REMOTE_ADDR];
  EXPOSURE exposureTable;
  PEXPOSURE pExposure = &exposureTable;
  int newPlateFlag = 0;
  int deletionFlag = 0;
  int restoreFlag = 0;    /* undelete */
  int multipleValidFlag = 0; /* Multiple valid entries for a given date */
  int multipleConditionValidFlag = 0;
  int multipleEventsValidFlag = 0;
  int commaFlag = 0;
  int commaFlag2 = 0;
  char queryString[MAX_QUERY_STRING];
  char tempString[MAX_QUERY_STRING];
  int res;
  char validString[MAX_QUERY_STRING];
  unsigned char digest[16];
  MD5state *s = NULL;
  int md5index;
  unsigned char testString[256*64];
  char pr64[MAX_MD5_LENGTH];
  int testLength;

  series[0] = 0;
  newseries[0] = 0;
  editdate[0] = 0;
  newplate[0] = 0;
  newpickednotes[0] = 0;
  newstacklocation[0] = 0;
  neweventnotes[0] = 0;
  editorname[0] = 0;
  editcomment[0] = 0;
  remoteaddr[0] = 0;
  logfile[0] = 0;
  editevent[0] = 0;
#if 0
  memset(testString,0,sizeof(testString));
  strcpy(testString,"xxxxx");
  testLength = strlen(testString);
  if (testLength > (MAX_MD5_LENGTH-1)) {
    printf("ERROR: MAX_MD5_LENGTH exceeded");
  } else {
    s = md5((unsigned char*)&testString[0],testLength,digest,s);
    enc64(pr64,digest,sizeof(digest));
    pr64[22] = '\0';  /* chop trailing == */
    printf("'%s'",pr64);
  }
  printf("\n");

  exit(-1);
#endif


  time(&startTime);
  ptr = localtime(&startTime);
  strftime(timestr,25,"%Y-%m-%dT%H:%M:%S", ptr);


  /* Loop through the arguments */
#if 0
  for (argctr = 0; argctr < argc; argctr++) {
    printf("argctr %d arg '%s'\n",argctr,argv[argctr]);
  }
#endif
  for (argv++; --argc > 0; argv++) {
    argstr = *argv;
    /* Decode arguments */
    if (argstr[0] != '-') {
      errorFlag = 1;
      printf("ERROR: unqualified argument %s argc: %d\n",argstr,argc);
    } else {
      while ((cmdchar = *++argstr) != 0) {
        switch(cmdchar) {
        case 's': /* series */
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for input file -%c\n",cmdchar);
            errorFlag = 1;
          } else {

            strncpy(series,*++argv,MAX_SERIES_STRING-1);
            if (strlen(series) >= MAX_SERIES_STRING-1) {
              printf("ERROR: MAX_SERIES_STRING too small for %s\n",*argv);
              return(-1);
            }
          }
          break;

        case 'n': /* plateNumber */
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for input file -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&plateNumber);
            if (nvals != 1) {
              printf("ERROR: Unable to decode plateNumber %s\n",*argv);
              errorFlag = 1;
            }
          }
          break;

        case 'd': /* date */
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for input file -%c\n",cmdchar);
            errorFlag = 1;
          } else {

            strncpy(editdate,*++argv,MAX_DATE_STRING-1);
            if (strlen(editdate) >= MAX_DATE_STRING-1) {
              printf("ERROR: MAX_DATE_STRING too small for %s\n",*argv);
              return(-1);
            }
          }
          break;

        case 't': /* editConditionVersionId */
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for input file -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&editConditionVersionId);
            if (nvals != 1) {
              printf("ERROR: Unable to decode editConditionVersionId %s\n",*argv);
              errorFlag = 1;
            }
            if (editConditionVersionId < 0) {
              editConditionVersionId = -1;
            }
          }
          break;

        case 'u': /* editEventsVersionId */
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for input file -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            nvals = sscanf(*++argv,"%d",&editEventsVersionId);
            if (nvals != 1) {
              printf("ERROR: Unable to decode editEventsVersionId %s\n",*argv);
              errorFlag = 1;
            }
            if (editEventsVersionId < 0) {
              editEventsVersionId = -1;
            }
          }
          break;


        case 'p': /* newplate */
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for input file -%c\n",cmdchar);
            errorFlag = 1;
          } else {

            strncpy(newplate,*++argv,MAX_PLATE_NAME-1);
            if (strlen(newplate) >= MAX_PLATE_NAME-1) {
              printf("ERROR: MAX_PLATE_NAME too small for %s\n",*argv);
              return(-1);
            }
          }
          break;

        case 'a': /* newvalidcondition */
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for input file -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            ++argv;
            if (strlen(*argv) > 0) {
              if (strcmp(*argv,"yes") == 0) {
                newvalidcondition = 1;
              } else if (strcmp(*argv,"no") == 0) {
                newvalidcondition = 0;
              } else {
                printf("ERROR: newvalidcondition %s is not 'yes' nor 'no'\n",*argv);
              }
            }
          }
          break;

        case 'b': /* newpickedcomment */
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for input file -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            ++argv;
            for (pickedCommentIndex = 0; pickedCommentIndex < PICKED_COMMENT_MAX; pickedCommentIndex++) {
              if (strstr(*argv,pickedCommentString[pickedCommentIndex]) != NULL) {
                newpickedcomment |= 1 << pickedCommentIndex;
              }
            }
          }
          break;

        case 'c': /* newpickednotes */
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for input file -%c\n",cmdchar);
            errorFlag = 1;
          } else {

            strncpy(newpickednotes,*++argv,MAX_STACKS_BUFFER-1);
            if (strlen(newpickednotes) >= MAX_STACKS_BUFFER-1) {
              printf("ERROR: MAX_STACKS_BUFFER too small for %s\n",*argv);
              return(-1);
            }
          }
          break;

        case 'j': /* newvalidevent */
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for input file -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            ++argv;
            if (strlen(*argv) > 0) {
              if (strcmp(*argv,"yes") == 0) {
                newvalidevent = 1;
              } else if (strcmp(*argv,"no") == 0) {
                newvalidevent = 0;
              } else {
                printf("ERROR: newvalidevent %s is not 'yes' nor 'no'\n",*argv);
              }
            }
          }
          break;


        case 'e': /* newpickedstatus */
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for input file -%c\n",cmdchar);
            errorFlag = 1;
          } else {
            ++argv;
            for (pickedStatusIndex = 0; pickedStatusIndex < PICKED_STATUS_MAX; pickedStatusIndex++) {
              if (strcmp(*argv,pickedStatusString[pickedStatusIndex]) == 0) {
                newpickedstatus = pickedStatusIndex;
                break;
              }
            }
          }
          break;

        case 'f': /* newstacklocation */
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for input file -%c\n",cmdchar);
            errorFlag = 1;
          } else {

            strncpy(newstacklocation,*++argv,MAX_STACKS_BUFFER-1);
            if (strlen(newstacklocation) >= MAX_STACKS_BUFFER-1) {
              printf("ERROR: MAX_STACKS_BUFFER too small for %s\n",*argv);
              return(-1);
            }
          }
          break;

        case 'g': /* neweventnotes */
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for input file -%c\n",cmdchar);
            errorFlag = 1;
          } else {

            strncpy(neweventnotes,*++argv,MAX_STACKS_BUFFER-1);
            if (strlen(neweventnotes) >= MAX_STACKS_BUFFER-1) {
              printf("ERROR: MAX_STACKS_BUFFER too small for %s\n",*argv);
              return(-1);
            }
          }
          break;


        case 'h': /* editorname */
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for input file -%c\n",cmdchar);
            errorFlag = 1;
          } else {

            strncpy(editorname,*++argv,MAX_STACKS_BUFFER-1);
            if (strlen(editorname) >= MAX_STACKS_BUFFER-1) {
              printf("ERROR: MAX_STACKS_BUFFER too small for %s\n",*argv);
              return(-1);
            }
          }
          break;

        case 'i': /* editcomment */
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for input file -%c\n",cmdchar);
            errorFlag = 1;
          } else {

            strncpy(editcomment,*++argv,MAX_STACKS_BUFFER-1);
            if (strlen(editcomment) >= MAX_STACKS_BUFFER-1) {
              printf("ERROR: MAX_STACKS_BUFFER too small for %s\n",*argv);
              return(-1);
            }
          }
          break;


        case 'k': /* remoteaddr */
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for input file -%c\n",cmdchar);
            errorFlag = 1;
          } else {

            strncpy(remoteaddr,*++argv,MAX_REMOTE_ADDR-1);
            if (strlen(remoteaddr) >= MAX_REMOTE_ADDR-1) {
              printf("ERROR: MAX_REMOTE_ADDR too small for %s\n",*argv);
              return(-1);
            }
          }
          break;

        case 'o': /* logfile */
          argc--;
          if (argc < 1) {
            printf("ERROR: Insufficient arguments for input file -%c\n",cmdchar);
            errorFlag = 1;
          } else {

            strncpy(logfile,*++argv,MAX_STACKS_BUFFER-1);
            if (strlen(logfile) >= MAX_STACKS_BUFFER-1) {
              printf("ERROR: MAX_STACKS_BUFFER too small for %s\n",*argv);
              return(-1);
            }
          }
          break;



        default:
          printf("ERROR: * illegal command -%c-\n",cmdchar);
          errorFlag = 1;

        }

      }

    }
  }

  if (logfile[0] == 0) {
    printf("ERROR: no log file specified\n");
    return(-1);
  }

  logHandle = fopen(logfile,"a+t");
  if (logHandle == NULL) {
    printf("ERROR opening the log file.  Other user may be active.  Retry later.\n");
    return(-1);
  }
  fprintf(logHandle,"%s '%s' '%s' ",timestr,remoteaddr,editorname);


  if (remoteaddr[0] == 0) {
    printf("ERROR: invalid remote host\n");
    fprintf(logHandle,"ERROR: invalid remote host\n");
    return(-1);
  }
  memset(testString,0,sizeof(testString));
  testLength = strlen(remoteaddr);
  if (testLength  > (MAX_MD5_LENGTH-1)) {
    printf("ERROR: invalid remote host length\n");
    fprintf(logHandle,"ERROR: invalid remote host length\n");
    return(-1);
  }
  strcpy((char *) testString,remoteaddr);
  testLength = strlen((char *) testString);
  s = md5((unsigned char*)&testString[0],testLength,digest,s);
  enc64(pr64,digest,sizeof(digest));
  pr64[22] = '\0';  /* chop trailing == */
  for (md5index = 0; md5index < md5sumTableSize; md5index++) {
    if (strcmp(pr64,md5sumTable[md5index]) == 0) {
      break;
    }
  }
  if (md5index >= md5sumTableSize) {
    printf("ERROR: unauthorized remote host '%s'\n",editcomment);
    fprintf(logHandle,"ERROR: unauthorized remote host '%s'\n",editcomment);
    return(-1);
  }


  if ((plateNumber <= 0) || (plateNumber >= MAX_PLATE_NUMBER)) {
    printf("ERROR: invalid plate number %d\n",plateNumber);
    fprintf(logHandle,"ERROR: invalid plate number %d\n",plateNumber);
    return(-1);
  }
  if ((editConditionVersionId < 0) && (editEventsVersionId < 0)) {
    printf("ERROR: No valid versionId found\n");
    fprintf(logHandle,"ERROR: No valid versionId found\n");
    return(-1);
  }


  if (strlen(editorname) <= 0) {
    printf("ERROR: No editor name specified\n");
    fprintf(logHandle,"ERROR: No editor name specified\n");
    return(-1);
  }

  if (ParseFilename2(newplate,newseries,&newPlateNumber) != 1)  {
    printf("ERROR: unable to parse new plate '%s'\n",newplate);
    fprintf(logHandle,"ERROR: unable to parse new plate '%s'\n",newplate);
    return(-1);
  }
  sprintf(tempplate,"%s%d",newseries,newPlateNumber);
  if (strcmp(tempplate,newplate) != 0) {
    sprintf(tempplate,"%s%05d",newseries,newPlateNumber);
    if (strcmp(tempplate,newplate) != 0) {
      printf("ERROR: extra characters in new plate '%s' (interpreted as '%s')\n",newplate,tempplate);
      fprintf(logHandle,"ERROR: extra characters in new plate '%s' (interpreted as '%s')\n",newplate,tempplate);
      return(-1);
    }
  }

  if (errorFlag != 0) {
    fprintf(logHandle,"ERROR: parse error\n");
    return(-1);
  }

  mysqlhost = getenv("DASCH_MYSQLHOST");
  if (mysqlhost == NULL) {
    printf("ERROR: DASCH_MYSQLHOST is not defined\n");
    fprintf(logHandle,"ERROR: DASCH_MYSQLHOST is not defined\n");
    return(-1);
  }

  // Connect to databases. These will abort the process if any unsolvable
  // problems occur.
  dasch_init_scandb(pConnection);
  dasch_init_photdb(pPhotConnection);

  stacksusername = getenv("DASCH_STACKS_USERNAME");
  if (stacksusername == NULL) {
    printf("DASCH_STACKS_USERNAME is not defined\n");
    fprintf(logHandle,"DASCH_STACKS_USERNAME is not defined\n");
    return(-1);
  }
  stackspassword = getenv("DASCH_STACKS_PASSWORD");
  if (stackspassword == NULL) {
    printf("DASCH_STACKS_PASSWORD is not defined\n");
    fprintf(logHandle,"DASCH_STACKS_PASSWORD is not defined\n");
    return(-1);
  }
  mysql_init(pStacksConnection);

  if (!mysql_real_connect(pStacksConnection,mysqlhost,stacksusername,stackspassword,"stacks",0,NULL,CLIENT_FOUND_ROWS)) {
    if (mysql_errno(pStacksConnection)) {
      printf("ERROR: MySQL error %d: %s\n",mysql_errno(pStacksConnection),mysql_error(pStacksConnection));
    }
    fprintf(logHandle,"ERROR: MySQL error %d: %s\n",mysql_errno(pStacksConnection),mysql_error(pStacksConnection));
    return(-1);
  }

  if (newvalidcondition >= 0) {
    condition_nrecs = GetPlateCondition(pStacksConnection,series,plateNumber,editdate,&plateConditionTable,logHandle);
  }
#if 0
  printf("%d condition records found\n",condition_nrecs);
#endif
  for (conditionIndex = 0; conditionIndex < condition_nrecs; conditionIndex++) {
    pPlateCondition = &plateConditionTable[conditionIndex];
    if ((editConditionVersionId >= 0) && (editConditionVersionId == pPlateCondition->versionId)) {
      editConditionIndex = conditionIndex;
      if ((pPlateCondition->valid > 0) && (newvalidcondition == 0)) {
        deletionFlag = 1;
      } else if ((pPlateCondition->valid == 0) && (newvalidcondition > 0)) {
        restoreFlag = 1;
      }
    }
    if (pPlateCondition->versionId > maxConditionVersionId) {
      maxConditionVersionId = pPlateCondition->versionId;
    }
    if (pPlateCondition->valid > 0) {
      if (curConditionVersionId < 0) {
        curConditionVersionId = pPlateCondition->versionId;
        curConditionIndex = conditionIndex;
      } else {
        multipleValidFlag = 1;
        multipleConditionValidFlag = 1;
      }
    }
  }


  if (newvalidevent >= 0) {
    events_nrecs = GetPlateEvents(pStacksConnection,series,plateNumber,editdate,&plateEventsTable,logHandle);
  }

#if 0
  printf("%d events records found\n",events_nrecs);
#endif
  for (eventsIndex = 0; eventsIndex < events_nrecs; eventsIndex++) {
    pPlateEvents = &plateEventsTable[eventsIndex];
    if ((editEventsVersionId >= 0) && (editEventsVersionId == pPlateEvents->versionId)) {
      editEventsIndex = conditionIndex;
      if ((pPlateEvents->valid > 0) && (newvalidevent == 0)) {
        deletionFlag = 1;
      } else if ((pPlateEvents->valid == 0) && (newvalidevent > 0)) {
        restoreFlag = 1;
      }
    }
    if (pPlateEvents->versionId > maxEventVersionId) {
      maxEventVersionId = pPlateEvents->versionId;
    }
    if (pPlateEvents->valid > 0) {
      if (curEventsVersionId < 0) {
        curEventsVersionId = pPlateEvents->versionId;
        curEventsIndex = eventsIndex;
      } else {
        multipleValidFlag = 1;
        multipleEventsValidFlag = 1;
      }
    }
  }


  InitSeriesTable(pConnection,pPhotConnection);
  strcpy(insertseries,series);
  insertPlateNumber = plateNumber;

  if ((strcmp(newseries,series) != 0) ||
      (newPlateNumber != plateNumber)) {
    newPlateFlag = 1;
    strcpy(insertseries,newseries);
    insertPlateNumber = newPlateNumber;
    if (deletionFlag != 0) {
      printf("ERROR: deletion allowed only for plate %s%05d, not new plate %s\n",series,plateNumber,newplate);
      fprintf(logHandle,"ERROR:  deletion allowed only for plate %s%05d, not new plate %s\n",series,plateNumber,newplate);
      return(-1);
    }
    if (restoreFlag != 0) {
      printf("ERROR: restore allowed only for plate %s%05d, not new plate %s\n",series,plateNumber,newplate);
      fprintf(logHandle,"ERROR: restore allowed only for plate %s%05d, not new plate %s\n",series,plateNumber,newplate);
      return(-1);
    }
    if (multipleValidFlag != 0) {
      printf("ERROR: rename allowed only for one valid entry for a given date for plate %s%05d, not new plate %s\n",series,plateNumber,newplate);
      fprintf(logHandle,"ERROR: rename allowed only for one valid entry for a given date for plate %s%05d, not new plate %s\n",series,plateNumber,newplate);
      return(-1);
    }


    if (GetSeriesId(newseries,0) < 1) {
      printf("ERROR: unknown new series %s\n", newseries);
      fprintf(logHandle,"ERROR: unknown new series %s\n", newseries);
      return(-1);
    }
    if (GetExposureInfo(pConnection,newseries,newPlateNumber,0,pExposure,0) != 1) {
      printf("ERROR: the new plate '%s' does not have a transcription\n",newplate);
      fprintf(logHandle,"ERROR: the new plate '%s' does not have a transcription\n",newplate);
      return(-1);
    }

    newcondition_nrecs = GetPlateCondition(pStacksConnection,newseries,newPlateNumber,editdate,&plateNewConditionTable,logHandle);
#if 0
    printf("%d newcondition records found\n",newcondition_nrecs);
#endif
    for (newconditionIndex = 0; newconditionIndex < newcondition_nrecs; newconditionIndex++) {
      pPlateNewCondition = &plateNewConditionTable[newconditionIndex];
      if (pPlateNewCondition->versionId > maxNewConditionVersionId) {
        maxNewConditionVersionId = pPlateNewCondition->versionId;
      }
      if (pPlateNewCondition->valid > 0) {
        printf("ERROR: a valid platecondition record with date %s already exists for new plate '%s'\n",editdate,newplate);
        fprintf(logHandle,"ERROR: a valid platecondition record with date %s already exists for new plate '%s'\n",editdate,newplate);
        return(-1);
      }
    }


    newevents_nrecs = GetPlateEvents(pStacksConnection,newseries,newPlateNumber,editdate,&plateNewEventsTable,logHandle);
#if 0
    printf("%d newevents records found\n",newevents_nrecs);
#endif
    for (neweventsIndex = 0; neweventsIndex < newevents_nrecs; neweventsIndex++) {
      pPlateNewEvents = &plateNewEventsTable[neweventsIndex];
      if (pPlateNewEvents->versionId > maxNewEventVersionId) {
        maxNewEventVersionId = pPlateNewEvents->versionId;
      }
      if (pPlateNewEvents->valid > 0) {
        printf("ERROR: a valid plateevents record with date %s already exists for new plate '%s'\n",editdate,newplate);
        fprintf(logHandle,"ERROR: a plateevents record with date %s already exists for new plate '%s'\n",editdate,newplate);
        return(-1);
      }
    }


  }

  if ((deletionFlag == 1) && (restoreFlag == 1)) {
    printf("ERROR: simultaneous restore and delete condition for %s%05d %s\n",series,plateNumber,editdate);
    fprintf(logHandle,"ERROR: simultaneous restore and delete condition for %s%05d %s\n",series,plateNumber,editdate);
    return(-1);
  }


  if (deletionFlag != 0) {
    strcpy(editevent,"Delete");
#if 0 /* Does not work for mc29946  because of a null date in one bad record but not the other */
    /* If we can delete a single record on both tables, then do it */
    if ((editConditionIndex < 0) && (multipleConditionValidFlag == 0) && (curConditionIndex >= 0)) {
      editConditionIndex = curConditionIndex;
      editConditionVersionId = curConditionVersionId;
    }
    if ((editEventsIndex < 0) && (multipleEventsValidFlag == 0) && (curEventsIndex >= 0)) {
      editEventsIndex = curEventsIndex;
      editEventsVersionId = curEventsVersionId;
    }
#endif
  } else if (restoreFlag != 0) {
    strcpy(editevent,"Restore");
    if (((editConditionIndex >= 0) && (curConditionIndex >= 0)) || (multipleConditionValidFlag == 1)) {
      printf("ERROR: invalid restore because valid condition records %d %d %d exist for %s%05d %s\n",editConditionVersionId,curConditionVersionId,multipleConditionValidFlag,series,plateNumber,editdate);
      fprintf(logHandle,"ERROR: invalid restore because valid condition records %d %d %d exist for %s%05d %s\n",editConditionVersionId,curConditionVersionId,multipleConditionValidFlag,series,plateNumber,editdate);
      return(-1);
    }
    if (((editEventsIndex >= 0) && (curEventsIndex >= 0)) || (multipleEventsValidFlag == 1)) {
      printf("ERROR: invalid restore because valid events records %d %d %d exist for %s%05d %s\n",editEventsVersionId,curEventsVersionId,multipleEventsValidFlag,series,plateNumber,editdate);
      fprintf(logHandle,"ERROR: invalid restore because valid events records %d %d %d exist for %s%05d %s\n",editEventsVersionId,curEventsVersionId,multipleEventsValidFlag,series,plateNumber,editdate);
      return(-1);
    }
  } else if (newPlateFlag != 0) {
    strcpy(editevent,"Rename");
  } else {
    strcpy(editevent,"Update");
    if (multipleValidFlag != 0) {
      printf("ERROR: update allowed only for one valid entry for a given date for plate %s%05d, not new plate %s\n",series,plateNumber,newplate);
      fprintf(logHandle,"ERROR: update allowed only for one valid entry for a given date for plate %s%05d, not new plate %s\n",series,plateNumber,newplate);
      return(-1);
    }
  }

  pPlateCondition = NULL;
  pPlateEvents = NULL;
  if ((deletionFlag == 0) && (restoreFlag == 0)) {
    if (curConditionVersionId >= 0) {
      pPlateCondition = &plateConditionTable[curConditionIndex];
      if (pPlateCondition->pickedComment != newpickedcomment) {
        conditionChangedFlag = 1;
      }
      if (strcmp(pPlateCondition->pickednotes,newpickednotes) != 0) {
        conditionChangedFlag = 1;
      }
    } else {
      if (newpickedcomment != 0) {
        conditionChangedFlag = 1;
      }
      if (strlen(newpickednotes) > 0) {
        conditionChangedFlag = 1;
      }
    }
    if (curEventsVersionId >= 0) {
      pPlateEvents = &plateEventsTable[curEventsIndex];
      if ((newpickedstatus == 0) &&
          /* change of Aug 18 2019: prevent changing externally curated pickedStatus values */
          ((strcmp(pickedStatusString[pPlateEvents->pickedStatus],"Jacket JPEG") == 0) ||
           (strcmp(pickedStatusString[pPlateEvents->pickedStatus],"Plate JPEG") == 0) ||
           (strcmp(pickedStatusString[pPlateEvents->pickedStatus],"Flood") == 0) ||
           (strcmp(pickedStatusString[pPlateEvents->pickedStatus],"Scanned") == 0))) {
        newpickedstatus = pPlateEvents->pickedStatus;
      }
      if (pPlateEvents->pickedStatus != newpickedstatus) {
        eventsChangedFlag = 1;
      }
      if (strcmp(pPlateEvents->stackLocation,newstacklocation) != 0) {
        eventsChangedFlag = 1;
      }
      if (strcmp(pPlateEvents->eventnotes,neweventnotes) != 0) {
        eventsChangedFlag = 1;
      }
    } else {
      if (newpickedstatus != 0) {
        eventsChangedFlag = 1;
      }
      if (strlen(newstacklocation) > 0) {
        eventsChangedFlag = 1;
      }
      if (strlen(neweventnotes) > 0) {
        eventsChangedFlag = 1;
      }
    }
    if ((newPlateFlag == 0)  && (conditionChangedFlag == 0) && (eventsChangedFlag == 0)) {
      printf("ERROR: no edit changes found for plate %s%05d \n",series,plateNumber);
      fprintf(logHandle,"ERROR: no change made for plate %s%05d \n",series,plateNumber);
      return(-1);
    }


    nextVersionId =   maxConditionVersionId;
    if (maxEventVersionId > nextVersionId) {
      nextVersionId = maxEventVersionId;
    }

    if (maxNewConditionVersionId > nextVersionId) { /* for rename */
      nextVersionId = maxNewConditionVersionId;
    }

    if (maxNewEventVersionId > nextVersionId) {
      nextVersionId = maxNewEventVersionId;
    }

    nextVersionId++;
    if (nextVersionId <= 0) {
      nextVersionId = 1;
    }
  } else {
    curConditionIndex = editConditionIndex;
    curConditionVersionId = editConditionVersionId;
    curEventsIndex = editEventsIndex;
    curEventsVersionId = editEventsVersionId;
  }



#if 0
  printf("ERROR: web_edit aborting transaction\n");
  exit(-1);
#endif
  fprintf(logHandle,"web_edit committing transaction\n");
  pPlateCondition = NULL;
  pPlateEvents = NULL;
  sprintf(queryString,"INSERT INTO edithistory (series,plateNumber,date,conditionVersionId,eventsVersionId,newVersionId,editdate,editorname,editevent) values ('%s',%d,'%s',%d,%d,%d,'%s','%s','%s');\n",
          series,
          plateNumber,
          timestr,
          curConditionVersionId,
          curEventsVersionId,
          nextVersionId,
          editdate,
          editorname,
          editevent);

  res = ExecuteQuery(pStacksConnection,queryString);
  if (res) {
    fprintf(logHandle,"%s '%s' '%s' ERROR executing %s",timestr,remoteaddr,editorname,queryString);
    return(-1);
  }
  fprintf(logHandle,"%s",queryString);

  commaFlag = 0;
  sprintf(queryString,"UPDATE edithistory SET ");
  if (newPlateFlag != 0) {
    sprintf(tempString,"newseries = '%s',newPlateNumber = %d",newseries,newPlateNumber);
    strcat(queryString,tempString);
    commaFlag = 1;
  }
  if (editcomment[0] != 0) {
    if (commaFlag == 1) {
      strcat(queryString,",");
    }
    sprintf(tempString,"editcomment = '%s'",editcomment);
    strcat(queryString,tempString);
    commaFlag = 1;
  }
  if (commaFlag == 1) {
    sprintf(tempString," WHERE series = '%s' and plateNumber = '%d' and date = '%s';\n",series,plateNumber,timestr);
    strcat(queryString,tempString);
    res = ExecuteQuery(pStacksConnection,queryString);
    if (res) {
      fprintf(logHandle,"%s '%s' '%s' ERROR executing %s",timestr,remoteaddr,editorname,queryString);
      return(-1);
    }
    fprintf(logHandle,"%s",queryString);
  }
  /* Write the new records */

  if ((deletionFlag == 0) && (restoreFlag == 0)) {
    if ((curConditionVersionId >= 0) || (conditionChangedFlag != 0)) {
      if ((newpickedcomment != 0) || (newpickednotes[0] != 0)) {


        sprintf(queryString,"INSERT INTO platecondition (series,plateNumber,date,versionId) values ('%s',%d,'%s',%d);\n",insertseries,insertPlateNumber,editdate,nextVersionId);
        res = ExecuteQuery(pStacksConnection,queryString);
        if (res) {
          fprintf(logHandle,"%s '%s' '%s' ERROR executing %s",timestr,remoteaddr,editorname,queryString);
          return(-1);
        }
        fprintf(logHandle,"%s",queryString);

        sprintf(queryString,"UPDATE platecondition SET ");
        commaFlag = 0;
        if (newpickedcomment != 0) {
          commaFlag = 1;
          commaFlag2 = 0;
          sprintf(tempString,"pickedComment='");
          for (pickedCommentIndex = 0; pickedCommentIndex < PICKED_COMMENT_MAX; pickedCommentIndex++) {
            if ((newpickedcomment & (1 << pickedCommentIndex)) != 0) {
              if (commaFlag2 != 0) {
                strcat(tempString,",");
              }
              commaFlag2 = 1;
              strcat(tempString,pickedCommentString[pickedCommentIndex]);
            }
          }
          strcat(tempString,"'");
          strcat(queryString,tempString);
        }
        if (newpickednotes[0] != 0) {
          if (commaFlag != 0) {
            strcat(queryString,",");
          }
          sprintf(tempString,"pickednotes='%s'",newpickednotes);
          strcat(queryString,tempString);
          commaFlag = 1;
        }
        if (commaFlag == 1) {
          sprintf(tempString," WHERE series = '%s' and plateNumber = '%d' and date = '%s' and versionId = %d;\n",insertseries,insertPlateNumber,editdate,nextVersionId);
          strcat(queryString,tempString);
          res = ExecuteQuery(pStacksConnection,queryString);
          if (res) {
            fprintf(logHandle,"%s '%s' '%s' ERROR executing %s",timestr,remoteaddr,editorname,queryString);
            return(-1);
          }
          fprintf(logHandle,"%s",queryString);
        }
      }

    }
    if ((curEventsVersionId >= 0) || (eventsChangedFlag != 0)) {
      sprintf(queryString,"INSERT INTO plateevents (series,plateNumber,date,versionId) values ('%s',%d,'%s',%d);\n",insertseries,insertPlateNumber,editdate,nextVersionId);
      res = ExecuteQuery(pStacksConnection,queryString);
      if (res) {
        fprintf(logHandle,"%s '%s' '%s' ERROR executing %s",timestr,remoteaddr,editorname,queryString);
        return(-1);
      }
      fprintf(logHandle,"%s",queryString);

      sprintf(queryString,"UPDATE plateevents SET ");
      commaFlag = 1;
      sprintf(tempString,"valid='yes',pickedStatus='%s'",pickedStatusString[newpickedstatus]);
      strcat(queryString,tempString);

      if (newstacklocation[0] != 0) {
        if (commaFlag != 0) {
          strcat(queryString,",");
        }
        sprintf(tempString,"stackLocation='%s'",newstacklocation);
        strcat(queryString,tempString);
        commaFlag = 1;
      }
      if (neweventnotes[0] != 0) {
        if (commaFlag != 0) {
          strcat(queryString,",");
        }
        sprintf(tempString,"eventnotes='%s'",neweventnotes);
        strcat(queryString,tempString);
        commaFlag = 1;
      }
      if (commaFlag == 1) {
        sprintf(tempString," WHERE series = '%s' and plateNumber = '%d' and date = '%s' and versionId = %d;\n",insertseries,insertPlateNumber,editdate,nextVersionId);
        strcat(queryString,tempString);
        res = ExecuteQuery(pStacksConnection,queryString);
        if (res) {
          fprintf(logHandle,"%s '%s' '%s' ERROR executing %s",timestr,remoteaddr,editorname,queryString);
          return(-1);
        }
        fprintf(logHandle,"%s",queryString);
      }


    }

  }




  /* Delete or restore the old records */
  if (restoreFlag) {
    strcpy(validString,"yes");
  } else {
    strcpy(validString,"no");
  }
  if (curConditionVersionId >= 0) {
    sprintf(queryString,"UPDATE platecondition SET valid = '%s' WHERE series = '%s' and plateNumber = %d and date = '%s' and versionId = '%d';\n",validString,series,plateNumber,editdate,curConditionVersionId);
    res = ExecuteQuery(pStacksConnection,queryString);
    if (res) {
      fprintf(logHandle,"%s '%s' '%s' ERROR executing %s",timestr,remoteaddr,editorname,queryString);
      return(-1);
    }
    fprintf(logHandle,"%s",queryString);

  }
  if (curEventsVersionId >= 0) {
    sprintf(queryString,"UPDATE plateevents SET valid = '%s' WHERE series = '%s' and plateNumber = %d and date = '%s' and versionId = '%d';\n",validString,series,plateNumber,editdate,curEventsVersionId);
    res = ExecuteQuery(pStacksConnection,queryString);
    if (res) {
      fprintf(logHandle,"%s '%s' '%s' ERROR executing %s",timestr,remoteaddr,editorname,queryString);
      return(-1);
    }
    fprintf(logHandle,"%s",queryString);
  }
  commaFlag = 0;







  if (plateConditionTable != NULL) {
    free(plateConditionTable);
    plateConditionTable = NULL;
  }
  if (plateEventsTable != NULL) {
    free(plateEventsTable);
    plateEventsTable = NULL;
  }
  if (plateNewConditionTable != NULL) {
    free(plateNewConditionTable);
    plateNewConditionTable = NULL;
  }
  if (plateNewEventsTable != NULL) {
    free(plateNewEventsTable);
    plateNewEventsTable = NULL;
  }
  mysql_close(pConnection);
  mysql_close(pPhotConnection);
  mysql_close(pStacksConnection);
  if (logHandle != NULL) {
    fclose(logHandle);
  }
  return(0);
}


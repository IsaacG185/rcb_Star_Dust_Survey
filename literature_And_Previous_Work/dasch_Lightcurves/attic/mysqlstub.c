// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/*
 *  mysqlstub.c - to allow compile against pipelineutils.a without MySQL
 *
 *  gcc -fPIC  -ggdb -O0   -D_FILE_OFFSET_BITS=64 -I/dasch/install/include  -I/usr/include/mysql  mysqlstub.c -c -o mysqlstub.o

 */
#include "mysql.h"
void STDCALL mysql_close(MYSQL *sock) {
  return;
}

unsigned int STDCALL mysql_errno(MYSQL *mysql) {
  return(0);
}

const char * STDCALL mysql_error(MYSQL *mysql) {
  return(0);
}

MYSQL_ROW	STDCALL mysql_fetch_row(MYSQL_RES *result) {
  return(0);
}

void		STDCALL mysql_free_result(MYSQL_RES *result) {
  return;
}

MYSQL *		STDCALL mysql_init(MYSQL *mysql) {
  return(0);
}

my_ulonglong STDCALL mysql_num_rows(MYSQL_RES *res) {
  return(0);
}

int		STDCALL mysql_query(MYSQL *mysql, const char *q) {
  return(0);
}

MYSQL_RES *     STDCALL mysql_store_result(MYSQL *mysql) {
  return(0);
}

MYSQL *		STDCALL mysql_real_connect(MYSQL *mysql, const char *host,
                                     const char *user,
                                     const char *passwd,
                                     const char *db,
                                     unsigned int port,
                                     const char *unix_socket,
                                     unsigned long clientflag) {
  return(0);
}


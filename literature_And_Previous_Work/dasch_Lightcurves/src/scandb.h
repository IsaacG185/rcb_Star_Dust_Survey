// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

#ifndef _DASCH_SCANDB_H
#define _DASCH_SCANDB_H 1

#include <mysql.h>

// This function will connect to the scanner SQL database. If a "too many
// connections" error is encountered, it will sleep for a random amount of time
// and retry, up to 5 times. If any other unhandled problem happens, the process
// will abort. The database connection parameters come from $DASCH_MYSQLHOST,
// $DASCH_USERNAME, and $DASCH_PASSWORD.
extern void dasch_init_scandb(MYSQL *conn);

#endif

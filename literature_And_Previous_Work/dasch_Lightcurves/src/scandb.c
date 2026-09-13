// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

#include <stdlib.h>  // getenv, srand, rand
#include <stdio.h> // fprintf, etc.
#include <errno.h> // errno
#include <string.h> // strerror
#include <time.h> // time
#include <unistd.h> // sleep

#include <mysql.h>
#include <errmsg.h> // CR_*
#include <mysqld_error.h> // ER_*


static int seeded_random = 0;
const static int N_ATTEMPTS = 5;


void
dasch_init_scandb(MYSQL *conn)
{
    char *host, *user, *pass;
    int attempt;
    unsigned int delay = 10;

    host = getenv("DASCH_MYSQLHOST");
    if (host == NULL) {
        fprintf(stderr, "fatal error: $DASCH_MYSQLHOST not defined\n");
        exit(1);
    }

    user = getenv("DASCH_USERNAME");
    if (user == NULL) {
        fprintf(stderr, "fatal error: $DASCH_USERNAME not defined\n");
        exit(1);
    }

    pass = getenv("DASCH_PASSWORD");
    if (pass == NULL) {
        fprintf(stderr, "fatal error: $DASCH_PASSWORD not defined\n");
        exit(1);
    }

    mysql_init(conn);

    for (attempt = 0; attempt < N_ATTEMPTS; attempt++) {
        unsigned int e;

        if (
            mysql_real_connect(
                conn,
                host,
                user,
                pass,
                "scanner",
                0,
                NULL,
                CLIENT_FOUND_ROWS
            )
        ) {
            return; // success!
        }

        e = mysql_errno(conn);

        if (e == 0) {
            fprintf(
                stderr,
                "fatal error: couldn\'t connect to scanner database: "
                "unspecified MySQL error\n"
            );
            exit(1);
        }

        if ((e == ER_CON_COUNT_ERROR || e == CR_SERVER_LOST) && attempt < N_ATTEMPTS - 1) {
            // error 1040 -- too many connections; error 2013 -- lost connection
            // to server. These are both observed to occur in cases where the
            // server is too busy, and it might help to sleep and try again.
            int d = delay;

            if (!seeded_random) {
                srand(time(NULL));
                seeded_random = 1;
            }

            d += rand() % 16;
            fprintf(
                stderr,
                "warning: trouble connecting to scanner database (error code %u); "
                "sleeping %u s and retrying\n",
                e,
                d
            );
            sleep(d);
            delay *= 2;
            continue;
        }

        // Otherwise, not an error we can handle. Abort.
        fprintf(
            stderr,
            "fatal error: couldn\'t connect to scanner database: "
            "%s (MySQL error %d)\n",
            mysql_error(conn),
            e
        );
        exit(1);
    }
}

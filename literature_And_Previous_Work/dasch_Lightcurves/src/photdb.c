// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

#include <stdlib.h>  // getenv
#include <stdio.h> // fprintf, etc.
#include <errno.h> // errno
#include <string.h> // strerror

#include <mysql.h>


static int magfiles_use_new_layout = 0;


void
dasch_init_photdb(MYSQL *conn)
{
    char *root, *host, *user, *pass, *dbstart, *dbname;
    char buf[1024];
    int n, overflow = 0;
    FILE *f;

    root = getenv("DASCH_PHOT_ROOT");
    if (!root || root == NULL) {
        fprintf(stderr, "fatal error: $DASCH_PHOT_ROOT not defined\n");
        exit(1);
    }

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

    n = snprintf(buf, sizeof(buf), "%s/dbname.txt", root);
    if (n < 1 || n >= sizeof(buf)) {
        fprintf(stderr, "fatal error: buffer overflow constructing dbname.txt path\n");
        exit(1);
    }

    f = fopen(buf, "r");
    if (f == NULL) {
        fprintf(stderr, "fatal error: failed to open `%s` for reading\n", buf);
        exit(1);
    }

    if (fgets(buf, sizeof(buf), f) == NULL) {
        char *msg;
        int err_num = errno;

        if (err_num == 0) {
            msg = "empty file";
        } else {
            msg = strerror(err_num);
        }

        fprintf(
            stderr,
            "fatal error: didn't get a line from $DASCH_PHOT_ROOT/dbname.txt: %s (error %d)\n",
            msg,
            err_num
        );
        fclose(f);
        exit(1);
    }

    fclose(f);

    dbstart = &buf[0];

    if (*dbstart == '*') {
        magfiles_use_new_layout = 1;
        dbstart++;
    }

    dbname = strtok(dbstart, " \t\n");
    if (dbname == NULL) {
        fprintf(stderr, "fatal error: no database name found in $DASCH_PHOT_ROOT/dbname.txt\n");
        exit(1);
    }

    // Finally we can get this started
    //
    // If needed, we could copy the retry code from scandb.c, but so far we
    // don't have lots of people trying to connect to this database at once.

    mysql_init(conn);

    if (
        !mysql_real_connect(
            conn,
            host,
            user,
            pass,
            dbname,
            0,
            NULL,
            CLIENT_FOUND_ROWS
        )
    ) {
        if (mysql_errno(conn)) {
            fprintf(
                stderr,
                "fatal error: couldn\'t connect to photometry database \"%s\": "
                "%s (MySQL error %d)\n",
                dbname,
                mysql_error(conn),
                mysql_errno(conn)
            );
        } else {
            fprintf(
                stderr,
                "fatal error: couldn\'t connect to photometry database \"%s\": "
                "unspecified MySQL error\n",
                dbname
            );
        }

        exit(1);
    }

    // For now (?), set up the legacy environment variables used to access the
    // photdb file trees. Note that putenv takes ownership of the buffer you
    // pass to it!

    n = snprintf(buf, sizeof(buf), "DASCH_PHOT_MAGNITUDES=%s/gsc2.3.2", root);
    if (n < 1 || n >= sizeof(buf) || putenv(strdup(buf))) {
        overflow = 1;
    }

    n = snprintf(buf, sizeof(buf), "DASCH_PHOT_MAGNITUDES1=%s/kepler", root);
    if (n < 1 || n >= sizeof(buf) || putenv(strdup(buf))) {
        overflow = 1;
    }

    n = snprintf(buf, sizeof(buf), "DASCH_PHOT_MAGNITUDES2=%s/apass", root);
    if (n < 1 || n >= sizeof(buf) || putenv(strdup(buf))) {
        overflow = 1;
    }

    n = snprintf(buf, sizeof(buf), "DASCH_PHOT_MAGNITUDES4=%s/atlas", root);
    if (n < 1 || n >= sizeof(buf) || putenv(strdup(buf))) {
        overflow = 1;
    }

    if (overflow) {
        fprintf(stderr, "fatal error: problem setting $DASCH_PHOT_MAGNITUDES path(s)\n");
        exit(1);
    }
}


int
dasch_photdb_magfiles_use_new_layout(void)
{
    return magfiles_use_new_layout;
}
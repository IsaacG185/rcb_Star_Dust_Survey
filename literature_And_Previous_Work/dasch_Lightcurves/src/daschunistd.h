// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

/* daschunistd.h
 *
 *  Because of the conflict between table.h and <unistd.h>, this file contains definitions normally appearing in unistd.h.  The file should be included in at least one instance of a program using unistd.h so that definition changes can be flagged
 *
 * Sep 15, 2015 Edward J. Los - Initial version
 */

#ifndef _DASCHUNISTD_H
#define _DASCHUNISTD_H 1
double exp10(double x); /* Not GNU standard. Use exp(x * ln(10.0)) otherwise. */
int close(int fd);
int unlink(const char *pathname);
off_t lseek(int fd, off_t offset, int whence);
ssize_t read(int fd, void *buf, size_t count);
ssize_t write(int fildes, const void *buf, size_t nbyte);
unsigned int sleep(unsigned int seconds);
#endif /* _DASCHUNISTD_H */

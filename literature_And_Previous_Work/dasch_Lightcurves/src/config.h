// Copyright the President and Fellows of Harvard College.
// Licensed under the MIT License

#ifndef __DASCH_PIPELINE_CONFIG_H__
#define __DASCH_PIPELINE_CONFIG_H__

/* Workarounds for the Starbase headers, which have the classic issue of
 * dependencies on should-be-private config.h macros. */

#ifndef __STDC__
#define __STDC__
#endif

#ifndef STDC_HEADERS
#define STDC_HEADERS 1
#endif

#define HAVE_MALLOC_H
#define HAVE_MEMCPY
#define HAVE_STDLIB_H

#endif
/*
//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 2002-2011, Sandia Corporation, Albuquerque, NM, USA.
// Under the terms of Contract DE-AC04-94AL85000, there is a
// non-exclusive license for use of this work by or on behalf of the
// U.S. Government.  Export of this program may require a license from
// the United States Government.
// 
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2011  Sandia Corporation
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $RCSfile: alloc.c,v $
//
// Purpose        : 
//
// Special Notes  : 
//
// Creator        : 
//
// Creation Date  : 
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.2 $
//
// Revision Date  : $Date: 2011/08/17 21:20:54 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------
*/

/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
**********/

/*
 * Memory alloction functions
 */

#include "spice.h"
#include "stdio.h"
#include "misc.h"
#include "util.h"
#ifdef DEBUG_MALLOC
#include "sm.h"
#endif

#include <stdio.h>

#ifndef HAS_BCOPY
#define bzero(s,n) memset(s,0,n)
#endif

/* Malloc num bytes and initialize to zero. Fatal error if the space can't
 * be malloc'd.   Return NULL for a request for 0 bytes.
 */
#undef SHARED_MEM

void bye_bye(i)
{
    printf ("inv = %d\n",1/i);
}
/*
*/
char *
tmalloc(num)
    int num;
{
    char *s;

    if (!num)
	return NULL;

#ifdef DEBUG_MALLOC
    s = sm_malloc((unsigned) num);
#else
    s = malloc((unsigned) num);
#endif
    if (!s) {
        fprintf(stderr, 
		"malloc: Internal Error: can't allocate %d bytes.\n", num);
        exit(EXIT_BAD);
    }

    bzero(s, num);

    return(s);
}

char *
trealloc(str, num)
    char *str;
    int num;
{
    char *s;

    if (!num) {
	if (str)
#ifdef SHARED_MEM
		FREE(str);
#else
#ifdef DEBUG_MALLOC
		sm_free(str);
#else
		free(str);
#endif
#endif
	return NULL;
    }
#ifdef SHARED_MEM
    if (!str)
	s = (char *) MALLOC_SM(num);
    else
        s = (char *) REALLOC_SM(str, (unsigned) num);
#else
    if (!str)
	s = tmalloc(num);
    else
#ifdef DEBUG_MALLOC
        s = sm_realloc(str, (unsigned) num);
#else
        s = realloc(str, (unsigned) num);
#endif
#endif
    if (!s) {
        fprintf(stderr, 
		"realloc: Internal Error: can't allocate %d bytes.\n", num);
        perror ("realloc");
        s = malloc((unsigned) num);
        bye_bye(0);
        fprintf (stderr, "From malloc of %d bytes: %lx\n",num,s);
        perror ("malloc");
        exit(EXIT_BAD);
    }
    return(s);
}

void
txfree(ptr)
	char	*ptr;
{
	if (ptr)
#ifdef SHARED_MEM
		FREE(ptr);
#else
#ifdef DEBUG_MALLOC
		sm_free(ptr);
#else
		free(ptr);
#endif
#endif
}

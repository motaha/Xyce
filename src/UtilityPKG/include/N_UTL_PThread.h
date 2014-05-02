//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
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
//-----------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_UTL_PThread.h,v $
//
// Purpose        : Stream buffer that performs indentationx
//
// Special Notes  :
//
// Creator        : David G. Baur  Raytheon  Sandia National Laboratories 1355 
//
// Creation Date  : 2013/04/18 18:01:27
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.2.2.2 $
//
// Revision Date  : $Date: 2014/03/03 18:29:29 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef Xyce_N_UTL_PThread_h
#define Xyce_N_UTL_PThread_h

#ifdef USE_THREADS
#include <pthread.h>
#endif

namespace Xyce {
namespace Util {

#ifdef USE_THREADS
typedef ::pthread_t xyce_pthread_t;
typedef ::pthread_attr_t xyce_pthread_attr_t;
typedef ::pthread_mutex_t xyce_pthread_mutex_t;

inline xyce_pthread_t xyce_pthread_self() 
{
  return ::pthread_self();
}

inline int xyce_pthread_mutex_lock(xyce_pthread_mutex_t *mutex) 
{
  return ::pthread_mutex_lock(mutex);
}

inline int xyce_pthread_mutex_unlock(xyce_pthread_mutex_t *mutex) 
{
  return ::pthread_mutex_unlock(mutex);
}

inline int xyce_pthread_create(xyce_pthread_t *__newthread, xyce_pthread_attr_t * __attr, void *(*__start_routine) (void *), void * __arg) 
{
  return ::pthread_create(__newthread, __attr, __start_routine, __arg);
}

inline int xyce_pthread_join(xyce_pthread_t __th, void **__thread_return) 
{
  return ::pthread_join(__th, __thread_return);
}

#else

typedef int xyce_pthread_t;
typedef int xyce_pthread_attr_t;
typedef int xyce_pthread_mutex_t;

inline xyce_pthread_t xyce_pthread_self() 
{
  return 1;
}

inline int xyce_pthread_mutex_lock(xyce_pthread_mutex_t *mutex) 
{
  return 0;
}

inline int xyce_pthread_mutex_unlock(xyce_pthread_mutex_t *mutex) 
{
  return 0;
}

inline int xyce_pthread_create(xyce_pthread_t *__newthread, xyce_pthread_attr_t * __attr, void *(*__start_routine) (void *), void * __arg) 
{
  (*__start_routine)(__arg);

  return 0;
}

inline int xyce_pthread_join(xyce_pthread_t __th, void **__thread_return) 
{
  return 0;
}

#endif

} // namespace Util
} // namespace Xyce

#endif // Xyce_N_UTL_PThread_h

//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2014 Sandia Corporation
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
// Filename       : $RCSfile: null_streambuf.h,v $
//
// Purpose        :
//
// Creator        : Scott A. Hutchinson, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/05/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.10.2.1 $
//
// Revision Date  : $Date: 2014/02/27 00:52:16 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef null_streambuf_h
#define null_streambuf_h

#ifdef HAVE_IOSTREAM
#include <iostream>
#else
#include <iostream.h>
#endif

#include <N_UTL_Xyce.h>

#if ! defined(streamsize)
#if ! defined(__SUNPRO_CC) || __SUNPRO_CC < 0x500
#define streamsize int
#else
#define streamsize streamsize
#endif
#endif

// Specialize the ANSI Standard C++ streambuf class that throws away everything
// given to it without generating an error.

class null_streambuf : public streambuf {
public:

  // Constructor
  null_streambuf();

  // Destructor
  virtual ~null_streambuf();

protected:

  // Called when output buffer is filled
  virtual int overflow(int c = EOF);

  // Sync is a no-op
  virtual int sync();

  // Setbuf is a no-op
  virtual streambuf * setbuf(char * s , streamsize n);

private:

  null_streambuf(const null_streambuf & ); // Not allowed
  null_streambuf & operator = (const null_streambuf & ); // Not allowed

  char buf[64]; // Throw away buffer
};

/*--------------------------------------------------------------------*/

#endif

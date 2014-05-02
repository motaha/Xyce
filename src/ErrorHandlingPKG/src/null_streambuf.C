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
// Filename       : $RCSfile: null_streambuf.C,v $
//
// Purpose        :
//
// Creator        : Scott A. Hutchinson, SNL, Parallel Computational Sciences
//
// Creation Date  : 06/02/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.7.2.1 $
//
// Revision Date  : $Date: 2014/02/27 00:52:17 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


#include <null_streambuf.h>

/*--------------------------------------------------------------------*/

null_streambuf::null_streambuf() : streambuf()
{
  setp( buf , buf + sizeof(buf) );
}

null_streambuf::~null_streambuf() {}

/*--------------------------------------------------------------------*/
/* Overflow */

int null_streambuf::overflow( int c )
{
  setp( buf , buf + sizeof(buf) );

  return c ;
}

/*--------------------------------------------------------------------*/

int null_streambuf::sync()
{
  return 0 ;
}

streambuf * null_streambuf::setbuf( char * s , streamsize n )
{
  return this ;
}

/*--------------------------------------------------------------------*/

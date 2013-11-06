//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2013  Sandia Corporation
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
// Filename       : $RCSfile: N_UTL_Xyce.h,v $
//
// Purpose        : This file is similar to Petra_Petra.h.  It  contains
//                  a number of Xyce-wide definitions and includes.
//
// Special Notes  :
//
// Creator        : Eric Keiter
//
// Creation Date  : 5/22/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.8.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:52 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef  _N_UTL_Xyce_H
#define  _N_UTL_Xyce_H

//DO NOT REMOVE
#include <string>
// Why?
// This is a terrible hack meant to tide the project over until misuse of
// "using namespace std" is completely removed.
// The reason?  On Windows, attempting to do a "using namespace std" before
// any include file has actually *defined* anything in that namespace is an
// error.  Because of that, it is a source of insanely cryptic windows build 
// errors when a Xyce developer inadvertently includes this header file as the 
// first line of any Xyce source file.
//  The correct solution to this issue is to REMOVE the dependence on this
//  using namespace thing, and to use "std::" where needed.  Then we can just
//  do away with N_UTL_Xyce.h entirely, as it serves no purpose other than to
//  obscure the using namespace.
using namespace std;

#endif


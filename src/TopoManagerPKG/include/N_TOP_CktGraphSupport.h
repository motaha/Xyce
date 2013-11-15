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

//-----------------------------------------------------------------------------
// Filename       : $RCSfile: N_TOP_CktGraphSupport.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/20/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.7.6.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:50 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef N_TOP_CktGraphSupport_h
#define N_TOP_CktGraphSupport_h 1

// ---------- Standard Includes ----------

#include <string>

// ---------- Xyce Includes ----------

#include <N_UTL_Xyce.h>

// ---------- Forward Declarations ----------

class N_TOP_CktGraphCreator;

//-----------------------------------------------------------------------------
// Class         : N_TOP_CktGraphSupport
// Purpose       : Factory class for CktGraphCreator
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
class N_TOP_CktGraphSupport
{

private:

  // Constructor (private).
  N_TOP_CktGraphSupport();

  // Copy constructor (private).
  N_TOP_CktGraphSupport(const N_TOP_CktGraphSupport & right);

  // Destructor (private).
  ~N_TOP_CktGraphSupport();

public:

  // Circuit graph factory.
  static N_TOP_CktGraphCreator * factory(const string & typeID, const int maxTries);

};

#endif

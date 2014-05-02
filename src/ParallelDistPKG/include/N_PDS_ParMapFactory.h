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

//-----------------------------------------------------------------------------
// Filename       : $RCSfile: N_PDS_ParMapFactory.h,v $
//
// Purpose        : Specification file for the parallel map abstract factory.
//
// Special Notes  : GoF Abstract Factory
//
// Creator        : Scott A. Hutchinson, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/08/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.14 $
//
// Revision Date  : $Date: 2014/02/24 23:49:25 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_PDS_ParMapFactory_h
#define Xyce_N_PDS_ParMapFactory_h

#include <vector>

#include <N_UTL_Xyce.h>
#include <N_PDS_fwd.h>

class N_PDS_ParMap;

//-----------------------------------------------------------------------------
// Class         : N_PDS_ParMapFactory
// Purpose       : Parallel map abstract factory (GoF).
// Special Notes :
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 03/08/00
//-----------------------------------------------------------------------------
class N_PDS_ParMapFactory
{

public:

  // Return a new ParMap
  static N_PDS_ParMap * create( int & numGlobalEntities,
                                const int & numLocalEntities,
                                const std::vector<int> & lbMap,
                                const int index_base = 0,
                                N_PDS_Comm * aComm = 0);

private:

  // Default constructor (private).
  N_PDS_ParMapFactory() { }

};

#endif

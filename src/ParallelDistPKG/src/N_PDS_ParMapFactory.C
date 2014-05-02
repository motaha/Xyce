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
// Filename       : $RCSfile: N_PDS_ParMapFactory.C,v $
//
// Purpose        : Implementation file for the parallel map abstract factory.
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
// Revision Number: $Revision: 1.13 $
//
// Revision Date  : $Date: 2014/02/24 23:49:25 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <N_PDS_ParMapFactory.h>
#include <N_PDS_ParMap.h>

#include <N_PDS_CommFactory.h>
#include <N_PDS_Comm.h>

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParMapFactory::create
// Purpose       : Creates a concrete N_PDS_ParMap object.  Takes a
//                 N_PDS_PetraComm object already instantiated.
// Special Notes : Part of a GoF Abstract Factory for creating parallel maps
//                 and associated communicators.
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 04/13/00
//-----------------------------------------------------------------------------
N_PDS_ParMap * N_PDS_ParMapFactory::create(
  		int & numGlobalEntities,
                const int & numLocalEntities,
                const std::vector<int> & lbMap,
                const int index_base,
  		N_PDS_Comm * aComm)
{
  return new N_PDS_ParMap(numGlobalEntities, numLocalEntities, lbMap, index_base, aComm);
}


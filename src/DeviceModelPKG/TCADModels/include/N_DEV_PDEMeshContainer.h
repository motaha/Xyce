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
// Filename       : $RCSfile: N_DEV_PDEMeshContainer.h,v $
//
// Purpose        : This is the base class for mesh processing/ownership.
//
// Special Notes  : Classes derived off of this one will be developed
//                  as needed.
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 04/21/02
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.8.2.1 $
//
// Revision Date  : $Date: 2014/02/26 20:16:31 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef Xyce_N_DEV_PDE_Mesh_Container_h
#define Xyce_N_DEV_PDE_Mesh_Container_h

#include <string>

#include <N_UTL_Misc.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : PDEMeshContainer
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 04/21/02
//-----------------------------------------------------------------------------
class PDEMeshContainer
{
public:
  PDEMeshContainer ();
  ~PDEMeshContainer ();

  virtual bool initializeMesh (const std::string & meshFileName);
};

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::PDEMeshContainer N_DEV_PDEMeshContainer;

#endif


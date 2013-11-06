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
// Filename       : $RCSfile: N_DEV_CCVS.C,v $
//
// Purpose        :
//
// Special Notes  : This device file is somewhat pointless.  CCVS devices
//                  are converted to Bsrc's behind the scenes, in the
//                  netlist parser.  As a result, this file never had
//                  to be fleshed out.  It should probably be removed.
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 02/28/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.12.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:38 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>



// N_DEV_CCVS
#include <N_DEV_CCVS.h>

namespace Xyce {
namespace Device {

// Class CCVS


//-----------------------------------------------------------------------------
// Function      : CCVS::CCVS
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
CCVS::CCVS(SolverState & ss1,
                       DeviceOptions & do1)
 : Device(ss1,do1)
{
}

//-----------------------------------------------------------------------------
// Function      : CCVS::CCVS
// Purpose       : copy constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
CCVS::CCVS(const CCVS &right)
  : Device(right)
{
}


//-----------------------------------------------------------------------------
// Function      : CCVS::~CCVS
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
CCVS::~CCVS()
{
}

// Additional Declarations

// Class CCVSInstance

//-----------------------------------------------------------------------------
// Function      : CCVSInstance::CCVSInstance
// Purpose       : copy constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
CCVSInstance::CCVSInstance(const CCVSInstance &right)
  :  DeviceInstance(right)
{

}

//-----------------------------------------------------------------------------
// Function      : CCVSInstance::~CCVSInstance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
CCVSInstance::~CCVSInstance()
{

}

// Additional Declarations

// Class CCVSModel

//-----------------------------------------------------------------------------
// Function      : CCVSModel::CCVSModel
// Purpose       : copy constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
CCVSModel::CCVSModel(const CCVSModel &right)
 : DeviceModel(right)
{
}


//-----------------------------------------------------------------------------
// Function      : CCVSModel::~CCVSModel
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
CCVSModel::~CCVSModel()
{
}

} // namespace Device
} // namespace Xyce

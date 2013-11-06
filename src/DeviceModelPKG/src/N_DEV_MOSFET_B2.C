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
// Filename       : $RCSfile: N_DEV_MOSFET_B2.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 02/28/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.11.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:38 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// N_DEV_MOSFET_B2
#include <N_DEV_MOSFET_B2.h>

namespace Xyce {
namespace Device {

// Class MOSFET_B2

//-----------------------------------------------------------------------------
// Function      : MOSFET_B2::MOSFET_B2
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
MOSFET_B2::MOSFET_B2 (SolverState & ss1,
                                  DeviceOptions & do1)
 : Device(ss1,do1)
{

}

//-----------------------------------------------------------------------------
// Function      : MOSFET_B2::MOSFET_B2
// Purpose       : copy constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
// MOSFET_B2::MOSFET_B2 (const MOSFET_B2 & right)
//  : Device(right)
// {

// }

//-----------------------------------------------------------------------------
// Function      : MOSFET_B2::~MOSFET_B2
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
MOSFET_B2::~MOSFET_B2 ()
{

}

// Additional Declarations

// Class MOSFET_B2Instance

//-----------------------------------------------------------------------------
// Function      : MOSFET_B2Instance::~MOSFET_B2Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
MOSFET_B2Instance::~MOSFET_B2Instance ()
{

}

//-----------------------------------------------------------------------------
// Function      : MOSFET_B2Model::~MOSFET_B2Model
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
MOSFET_B2Model::~MOSFET_B2Model ()
{
}

} // namespace Device
} // namespace Xyce


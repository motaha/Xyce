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
// Filename       : $RCSfile: N_DEV_DeviceTable.C,v $
//
// Purpose        : Describe the purpose of the contents of the file. If the
//                  contains the header file of a class, provide a clear
//                  description of the nature of the class.
//
// Special Notes  : Specify any "hidden" or subtle details of the class here.
//                  Portability details, error handling information, etc.
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 02/28/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.7.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:38 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


#include <N_DEV_DeviceTable.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : DeviceTable::DeviceTable
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
DeviceTable::DeviceTable(SolverState & ss1,
                                     DeviceOptions & do1)
 : Device(ss1,do1)
{

}

//-----------------------------------------------------------------------------
// Function      : DeviceTable::DeviceTable
// Purpose       : copy constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
DeviceTable::DeviceTable(const DeviceTable &right)
 : Device(right)
{

}


//-----------------------------------------------------------------------------
// Function      : DeviceTable::~DeviceTable
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
DeviceTable::~DeviceTable()
{

}

} // namespace Device
} // namespace Xyce

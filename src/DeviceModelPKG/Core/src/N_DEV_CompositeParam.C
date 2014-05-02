//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright(C) 2002-2014 Sandia Corporation
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
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
// Filename       : $RCSfile: N_DEV_CompositeParam.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Dave Shirley, PSSI
//
// Creation Date  : 05/05/05
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.27.2.2 $
//
// Revision Date  : $Date: 2014/03/14 19:48:40 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_DEV_CompositeParam.h>
#include <N_ERH_Message.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : CompositeParam::given
// Purpose       : Return whether param was given
// Special Notes :
// Scope         : public
// Creator       : Dave Baur, Raytheon
// Creation Date : 1/29/2014
//-----------------------------------------------------------------------------
///
/// given returns true if the value was specified in the netlist (not defaulted).
///
/// @param parameter_name    const reference to the name of the parameter
///
/// @return true if the value was specified in the netlist (was not defaulted)
bool CompositeParam::given(const std::string &parameter_name ) const
{
  ParameterMap::const_iterator it = getParameterMap().find(parameter_name);

  if (it == getParameterMap().end())
    Report::DevelFatal0() << "CompositeParam::Given: unrecognized param: " << parameter_name;

  return Xyce::Device::wasValueGiven(*this, (*it).second->getSerialNumber());
}

} // namespace Device
} // namespace Xyce

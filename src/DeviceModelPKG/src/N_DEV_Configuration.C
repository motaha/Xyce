//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2011  Sandia Corporation
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
// Filename       : $RCSfile: N_DEV_Configuration.C,v $
//
// Purpose        : 
//
// Special Notes  : 
//
// Creator        : David Baur
//
// Creation Date  : 04/18/2013
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.2.2.3 $
//
// Revision Date  : $Date: 2013/10/03 17:23:38 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <ostream>
#include <iomanip>

#include <N_DEV_Configuration.h>
#include <N_UTL_IndentStreamBuf.h>

namespace Xyce {
namespace Device {

std::ostream &outputConfiguration(std::ostream &os, const Configuration &configuration)
{
  os << "Nodes: " << configuration.numNodes << std::endl
     << "Optional Nodes: " << configuration.numOptionalNodes << std::endl
//     << "Model Required: " << (configuration.modelRequired ? "yes" : "no") << std::endl
     << "Primary Parameter: " << (configuration.primaryParameter.empty() ? "<none>" : configuration.primaryParameter) << std::endl
     << "Model Types: ";

  if (configuration.modelTypes.empty())
    os << "<none>";
  else {
    for (std::vector<std::string>::const_iterator it = configuration.modelTypes.begin(); it != configuration.modelTypes.end(); ++it) {
      if (it != configuration.modelTypes.begin())
        os << ", ";
      os << (*it);
    }
  }

  os << std::endl;

  return os;
}

} // namespace Device
} // namespace Xyce

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
// Filename       : $RCSfile: N_DEV_Configuration.h,v $
//
// Purpose        : Holds device configuration information
//
// Special Notes  :
//
// Creator        : David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
//
// Creation Date  : 2013/04/18 18:01:27
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.2.2.3 $
//
// Revision Date  : $Date: 2013/10/03 17:23:36 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Configuration_h
#define Xyce_N_DEV_Configuration_h

#include <iosfwd>
#include <string>
#include <vector>

namespace Xyce {
namespace Device {

/**
 * Class Configuration contains device instance configuration data.
 *
 */
class Configuration
{
  public:
    Configuration ()
      : numNodes(0),
        numOptionalNodes(0),
        numFillNodes(0),
        modelRequired(0),
        primaryParameter(""),
        modelTypes()
    {}

    Configuration (const Configuration &configuration)
      : numNodes(configuration.numNodes),
        numOptionalNodes(configuration.numOptionalNodes),
        numFillNodes(configuration.numFillNodes),
        modelRequired(configuration.modelRequired),
        primaryParameter(configuration.primaryParameter),
        modelTypes(configuration.modelTypes)
    {}

  public:
    int numNodes;
    int numOptionalNodes;
    int numFillNodes;
    int modelRequired;

  public:
    std::string primaryParameter;
    std::vector<std::string> modelTypes;
};

std::ostream &outputConfiguration(std::ostream &os, const Configuration &configuration);

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::Configuration Configuration;

#endif // Xyce_N_DEV_Configuration_h

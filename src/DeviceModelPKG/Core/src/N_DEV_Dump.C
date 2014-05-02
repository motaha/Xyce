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
// Filename       : $RCSfile: N_DEV_Dump.C,v $
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
// Revision Number: $Revision: 1.1 $
//
// Revision Date  : $Date: 2014/02/12 19:17:28 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <algorithm>
#include <iomanip>
#include <ostream>
#include <stdexcept>

#include <N_DEV_Dump.h>
#include <N_DEV_Configuration.h>
#include <N_UTL_IndentStreamBuf.h>
#include <N_UTL_Demangle.h>

namespace Xyce {
namespace Device {

struct type_compare_less
{
  bool operator()(const std::type_info *t0, const std::type_info *t1)
  {
    return t0->before(*t1);
  }
};

std::map<const std::type_info *, std::string, type_compare_less> s_typeInfoNameMap;

std::ostream &printTypeName(std::ostream &os, const std::type_info &type)
{
  if (s_typeInfoNameMap.empty()) {
    s_typeInfoNameMap[&typeid(bool)] = "boolean";
    s_typeInfoNameMap[&typeid(double)] = "double";
    s_typeInfoNameMap[&typeid(int)] = "int";
    s_typeInfoNameMap[&typeid(std::string)] = "string";
    s_typeInfoNameMap[&typeid(std::vector<double>)] = "double vector";
    s_typeInfoNameMap[&typeid(std::vector<int>)] = "int vector";
    s_typeInfoNameMap[&typeid(std::vector<std::string>)] = "string vector";
  }

  if (s_typeInfoNameMap[&type].empty())
    os << "composite"; // demangle(type.name());
  else
    os << s_typeInfoNameMap[&type];

  return os;
}

std::ostream &outputParameterMap(std::ostream &os, const ParameterMap &parameter_map);

std::ostream &outputDescriptor(std::ostream &os, const Descriptor &descriptor)
{
  if (&descriptor.getEntry()) {
    printTypeName(os, descriptor.getEntry().type());
  }

  if (!descriptor.getCompositeParametricData<void>()) {
    os << ", default ";
    descriptor.getEntry().print(os);
    if (descriptor.hasOriginalValueStored())
      os << ", original value managed, scaling enabled";
  }
  else {
    const ParametricData<void> &composite_parametric_data = *descriptor.getCompositeParametricData<void>();
    os << Util::push << std::endl;
    outputParameterMap(os, composite_parametric_data.getMap());
    os << Util::pop;
  }

  os << std::endl;

  return os;
}

std::ostream &outputParameterMap(std::ostream &os, const ParameterMap &parameter_map)
{
  for (ParameterMap::const_iterator it = parameter_map.begin(); it != parameter_map.end(); ++it) {
    os << (*it).first << ", ";

    outputDescriptor(os, *(*it).second);
  }

  return os;
}

std::ostream &operator<<(std::ostream &os, const Configuration &configuration)
{
  ParametricData<void> &instance_parameters = configuration.getInstanceParameters();
  ParametricData<void> &model_parameters = configuration.getModelParameters();

  os << "Configuration" << Xyce::Util::push << std::endl
     << "Name: " << configuration.getName() << std::endl
     << "Device Type: " << configuration.getDeviceTypeName() << std::endl
     << "Nodes: " << configuration.getNumNodes() << std::endl
     << "Optional Nodes: " << configuration.getNumOptionalNodes() << std::endl
     << "Fill Nodes: " << configuration.getNumFillNodes() << std::endl
     << "Model Required: " << (configuration.getModelRequired() ? "yes" : "no") << std::endl
     << "Linear Device: " << (configuration.getLinearDevice() ? "yes" : "no") << std::endl
     << "PDE Device: " << (configuration.getPDEDevice() ? "yes" : "no") << std::endl
     << "Primary Parameter: " << (configuration.getPrimaryParameter().empty() ? "<none>" : configuration.getPrimaryParameter()) << std::endl
     << "Instance Default Parameter: " << (configuration.getInstanceDefaultParameterName().empty() ? "<none>" : configuration.getInstanceDefaultParameterName()) << std::endl;

  os << "Model Types: ";
  for (std::vector<std::string>::const_iterator it = configuration.getModelTypeNames().begin(); it != configuration.getModelTypeNames().end(); ++it) {
    if (it != configuration.getModelTypeNames().begin())
      os << ", ";
    os << *it;
  }
  os << Xyce::Util::pop << std::endl
     << "Model Parameters" << Xyce::Util::push << std::endl;
  outputParameterMap(os, model_parameters.getMap());
  os << Xyce::Util::pop << std::endl
     <<  "Instance Parameters" << Xyce::Util::push << std::endl;
  outputParameterMap(os, instance_parameters.getMap());

  os << Xyce::Util::pop << std::endl
     << Xyce::Util::pop << std::endl;

  return os;
}

} // namespace Device
} // namespace Xyce

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
// Filename       : $RCSfile: N_DEV_LaTexDoc.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : David Baur
//
// Creation Date  : 3/20/2013
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.1.2.5 $
//
// Revision Date  : $Date: 2014/03/17 21:28:42 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

// All Xyce source files must include this!
#include <Xyce_config.h>

#include <string>
#include <map>
#include <vector>

// need this on Losedows for "std::tolower"
#ifdef HAVE_CCTYPE
#include <cctype>
#endif

#include <N_DEV_LaTexDoc.h>

#include <N_DEV_fwd.h>
#include <N_DEV_Units.h>
#include <N_DEV_Const.h>
#include <N_DEV_DeviceEntity.h>
#include <N_DEV_DeviceModel.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_SolverState.h>
#include <N_UTL_Misc.h>
#include <N_UTL_IndentStreamBuf.h>
#include <N_ERH_Message.h>

namespace Xyce {
namespace Device {

namespace {

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::escape
// Purpose       : Add escape backslash before underscore for LaTeX
// Special Notes :
// Scope         : private
// Creator       : Dave Shirley, PSSI
// Creation Date : 05/16/06
//-----------------------------------------------------------------------------
struct laTexProtect {
  laTexProtect(const std::string &t)
    : text(t)
  {}

  const std::string &text;
};


std::ostream &operator<<(std::ostream &os, const laTexProtect &laTex_protect)
{
  for (std::string::const_iterator cit = laTex_protect.text.begin(); cit != laTex_protect.text.end(); ++cit) {
    if ((*cit) == '_')
      os << "\\_\\-";
    else
      os << (*cit);
  }

  return os;
}


struct laTexClean {
  laTexClean(const std::string &t)
    : text(t)
  {}

  const std::string &text;
};


std::ostream &operator<<(std::ostream &os, const laTexClean &laTex_clean)
{
  for (std::string::const_iterator cit = laTex_clean.text.begin(); cit != laTex_clean.text.end(); ++cit) {
    if ((*cit) != '_')
      os << (*cit);
  }

  return os;
}


const std::string &categoryName(int category) {
  static std::vector<std::string> s_categoryList;

  if (s_categoryList.empty()) {
    s_categoryList.resize(CAT_MAX);
    s_categoryList[CAT_NONE] = "";
    s_categoryList[CAT_AC] = "AC Parameters";
    s_categoryList[CAT_BASIC] = "Basic Parameters";
    s_categoryList[CAT_BIN] = "Bin Parameters";
    s_categoryList[CAT_CAP] = "Capacitance Parameters";
    s_categoryList[CAT_CONTROL] = "Control Parameters";
    s_categoryList[CAT_CURRENT] = "Current Parameters";
    s_categoryList[CAT_DC] = "DC Parameters";
    s_categoryList[CAT_DEPENDENCY] = "Dependency Parameters";
    s_categoryList[CAT_DOPING] = "Doping Parameters";
    s_categoryList[CAT_FLICKER] = "Flicker Parameters";
    s_categoryList[CAT_GEOMETRY] = "Geometry Parameters";
    s_categoryList[CAT_INITIAL] = "Initial Condition Parameters";
    s_categoryList[CAT_MATERIAL] = "Material Parameters";
    s_categoryList[CAT_NQS] = "NQS Parameters";
    s_categoryList[CAT_RADP] = "Radiation Pulse Parameters";
    s_categoryList[CAT_RES] = "Resistance Parameters";
    s_categoryList[CAT_PROCESS] = "Process Parameters";
    s_categoryList[CAT_RF] = "RF Parameters";
    s_categoryList[CAT_RAD] = "Radiation Parameters";
    s_categoryList[CAT_TEMP] = "Temperature Parameters";
    s_categoryList[CAT_TUNNEL] = "Tunnelling Parameters";
    s_categoryList[CAT_VBI] = "Built-in Potential Lowering Parameters";
    s_categoryList[CAT_VOLT] = "Voltage Parameters";
    s_categoryList[CAT_ASYMRDS] = "Asymmetric and Bias-Dependent $R_{ds}$ Parameters";
    s_categoryList[CAT_IMPACT] = "Impact Ionization Current Parameters";
    s_categoryList[CAT_GDLEAKAGE] = "Gate-induced Drain Leakage Model Parameters";
  }

  return s_categoryList.size() ? s_categoryList[category] : s_categoryList[CAT_UNKNOWN];
}

void unitDescription(
  const std::string &   name,
  const ParameterMap &  parameter_map,
  const Descriptor &    descriptor,
  ParameterUnit &       unit,
  ParameterCategory &   category,
  std::string &         description)
{
  unit = descriptor.getUnit();
  if (unit == STANDARD)
  {
    for (const StdDescription *it = Units::descriptionTable; it != Units::descriptionTable + Units::descriptionTableSize; ++it) {
      if (name == (*it).Name)
      {
        description = (*it).Description;
        category = (*it).Category;
        unit = (*it).Unit;
        break;
      }
    }

    if (description.empty())
    {
      Report::UserWarning0() << "No description given for " << name << ".  Setting category and unit as UNKNOWN.";

      description = "Unspecified";
      category = CAT_UNKNOWN;
      unit = U_UNKNOWN;
    }
  }
  else
  {
    category = descriptor.getCategory();
    if (category == CAT_DEPENDENCY)
    {
      int num;
      std::string base_name = name.substr(1);
      if (name[0] == 'L')
      {
        description = "Length";
        num = 1;
      }
      else if (name[0] == 'W')
      {
        description = "Width";
        num = 1;
      }
      else if (name[0] == 'P')
      {
        description = "Cross-term";
        num = 2;
      }
      description += " dependence of ";
      description += base_name;

      ParameterMap::const_iterator it = parameter_map.find(base_name);
      if (it == parameter_map.end())
      {
        Report::DevelFatal().in("DeviceEntity::getUnitDescription")
          << "Base parameter not found for " << name << " in " <<  base_name;
      }

      else {
        const Descriptor &base_descriptor = *(*it).second;

        std::string base_description;
        ParameterUnit base_unit = U_INVALID;
        ParameterCategory base_category = CAT_UNKNOWN;

        unitDescription(base_name, parameter_map, base_descriptor, base_unit, base_category, base_description);
        if (base_unit != U_INVALID && base_unit != U_UNKNOWN)
        {
          for (int i = 0; i < num; ++i)
          {
            for (const UnitInfo *it = Units::unitTable; it != Units::unitTable + Units::unitTableSize; ++it)
            {
              if ((*it).Unit == base_unit)
              {
                unit = (*it).UnitM;
                if (unit == U_INVALID)
                {
                  Report::UserWarning0() << "Need unit for " << (*it).description << " times m";
                }
                break;
              }
            }
          }
        }
      }
    }
    else
    {
      description = descriptor.getDescription();
    }
  }
}

const UnitInfo &findUnit(const int unit)
{
  for (const UnitInfo *it = Units::unitTable; it != Units::unitTable + Units::unitTableSize; ++it)
    if ((*it).Unit == unit)
      return *it;
  return *(Units::unitTable);
}

std::ostream &documentParameter(std::ostream &os, const std::string &name, const int parameter_unit, const std::string &description, const Descriptor &descriptor)
{
  os << laTexProtect(name) << " & " << laTexProtect(description) << " & ";

  const UnitInfo &unit_info = findUnit(parameter_unit);
  os << unit_info.doc;
  os << " & ";

  if ((descriptor.isType<double>() && name == "TNOM") || name == "TEMP") 
  {
    double default_value = getDefaultValue<double>(descriptor);
    // os << default_value - CONSTCtoK;
    // os << CONSTREFTEMP - CONSTCtoK;
    if (default_value == 0.0)
    {
      os << "Ambient Temperature";
    }
    else if (default_value - CONSTCtoK < 0.0)
    {
      os << default_value;
    }
    else
    {
      os << default_value - CONSTCtoK;
    }
  }
  else
  {
    descriptor.getEntry().print(os);
  }

  os << " \\\\ \\hline" << std::endl;

  return os;
}

} // namespace <unnamed>


std::ostream &laTexDevice(std::ostream &os, const std::string &device_name, const int device_level, const int type, const std::string &device_description, const ParametricData<void> &parameters, OutputMode::Mode format)
{
  // Place caption on table and add index entry
  std::string device_description_lc = device_description;
  std::transform(device_description_lc.begin(), device_description_lc.end(), device_description_lc.begin(), (int (*)(int)) std::tolower);

  os << "\\index{" << laTexClean(device_description_lc) << "!" << (type ? "device model parameters" : "device instance parameters") << "}" << std::endl
     << "\\begin{DeviceParamTableGenerated}{" << laTexProtect(device_description) << " " << (type ? "Device Model Parameters" : "Device Instance Parameters") << "}"
     << "{" << device_name << "_" << device_level << (type ? "_Device_Model_Params" : "_Device_Instance_Params") << "}" << std::endl;

  // If catagorical listing, group by category
  if (format == OutputMode::DOC_CAT) {
    const ParameterMap &parameter_map = parameters.getMap();
    for (int category = CAT_NONE; category != CAT_MAX; ++category) {
      const std::string &header = categoryName(category);
      bool header_printed = false;

      for (ParameterMap::const_iterator it = parameter_map.begin(); it != parameter_map.end(); ++it) {
        std::string parameter_name = (*it).first;
        const Descriptor &descriptor = *(*it).second;

        if ((descriptor.getExpressionAccess() & ParameterType::NO_DOC) == 0) {
          std::string parameter_description;
          ParameterUnit parameter_unit = U_INVALID;
          ParameterCategory parameter_category = CAT_UNKNOWN;
          unitDescription((*it).first, parameter_map, descriptor, parameter_unit, parameter_category, parameter_description);

          if (parameter_category == category) {
            if (!header.empty() && !header_printed) {
              header_printed = true;

              os << std::endl
                 << "\\category{" << header << "}" << "\\\\ \\hline" << std::endl;
            }
            documentParameter(os, parameter_name, parameter_unit, parameter_description, descriptor);
          }
        }
      }
    }
  }
  else {
    const ParameterMap &parameter_map = parameters.getMap();
    for (ParameterMap::const_iterator it = parameter_map.begin(); it != parameter_map.end(); ++it) {
      std::string parameter_name = (*it).first;
      const Descriptor &descriptor = *(*it).second;

      if ((descriptor.getExpressionAccess() & ParameterType::NO_DOC) == 0) {

        std::string parameter_description;
        ParameterUnit parameter_unit = U_INVALID;
        ParameterCategory parameter_category = CAT_UNKNOWN;
        unitDescription((*it).first, parameter_map, descriptor, parameter_unit, parameter_category, parameter_description);

        documentParameter(os, parameter_name, parameter_unit, parameter_description, descriptor);
      }
    }
  }

  os << "\\end{DeviceParamTableGenerated}" << std::endl;

  return os;
}

} // namespace Device
} // namespace Xyce

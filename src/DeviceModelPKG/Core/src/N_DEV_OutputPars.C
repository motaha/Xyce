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
// Filename       : $RCSfile: N_DEV_OutputPars.C,v $
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
// Revision Number: $Revision: 1.4.2.10 $
//
// Revision Date  : $Date: 2013/12/08 23:38:51 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <string>
#include <map>
#include <vector>

#include <N_DEV_fwd.h>
#include <N_DEV_Units.h>
#include <N_DEV_Const.h>
#include <N_DEV_Factory.h>
#include <N_DEV_OutputPars.h>
#include <N_DEV_DeviceEntity.h>
#include <N_DEV_DeviceLevelKey.h>
#include <N_DEV_DeviceModel.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_SolverState.h>
#include <N_UTL_Misc.h>
#include <N_UTL_IndentStreamBuf.h>
#include <N_UTL_Demangle.h>

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
      os << "\\";
    os << (*cit);
  }

  return os;
}


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
    s_typeInfoNameMap[&typeid(CompositeMap)] = "composite";
  }

  if (s_typeInfoNameMap[&type].empty())
    os << demangle(type.name());
  else
    os << s_typeInfoNameMap[&type];

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
  const std::string &                           name,
  const ParametricData<void>::ParameterMap &    parameter_map,
  const Descriptor &                            descriptor,
  ParameterUnit &                               unit,
  ParameterCategory &                           category,
  std::string &                                 description)
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
      string msg("DeviceEntity::outputParams standard unit/description ");
      msg += "invoked for unrecognized name: ";
      msg += name;
      msg += ".  Setting category and unit as UNKNOWN.";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_WARNING, msg);

      description = string ("Unspecified");
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
      string base_name = name.substr(1);
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

      ParametricData<void>::ParameterMap::const_iterator it = parameter_map.find(base_name);
      if (it == parameter_map.end())
      {
        string msg("DeviceEntity::getUnitDescription: Base parameter not found for: ");
        msg += name;
        msg += " in ";
        msg += name;
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
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
                  string msg("DeviceEntity::getUnitDescription: need unit for ");
                  msg += (*it).description;
                  msg += " times m";
                  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_WARNING, msg);
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
  os << std::endl
     << "{\\tt " << laTexProtect(name) << "} & " << laTexProtect(description) << " & ";

  const UnitInfo &unit_info = findUnit(parameter_unit);
  os << unit_info.doc;
  os << " & ";

  if (descriptor.isType<double>() && name == "TNOM" || name == "TEMP") {
    double default_value = getDefaultValue<double>(descriptor);
    // os << default_value - CONSTCtoK;
    // os << CONSTREFTEMP - CONSTCtoK;
    if (default_value == 0.0)
      os << "Ambient Temperature";
    else if (default_value - CONSTCtoK < 0.0)
      os << default_value;
    else
      os << default_value - CONSTCtoK;
  }
  else
    descriptor.getEntry().print(os);

  os << " \\\\ \\hline" << endl;

  return os;
}

} // namespace <unnamed>


std::ostream &laTexDevice(std::ostream &os, const std::string &device_name, const int device_level, const int type, const std::string &device_description, const ParametricData<void> &parameters, OutputMode::Mode format)
{
  // Create table
  os << "% Column sizes are specified as percentages of the table width." << std::endl
     << "% If you need to tweak column widths, modify the coefficients" << std::endl
     << "% of \\hsize while keeping the sum of coefficients equal to one." << std::endl
     << "\\small" << std::endl
     << "\\begin{longtable}[H,t,b,p]{|>{\\setlength{\\hsize}{0.15\\hsize}}Y|" << std::endl
     << ">{\\setlength{\\hsize}{.65\\hsize}}Y|" << std::endl
     << ">{\\setlength{\\hsize}{.1\\hsize}}Y|" << std::endl
     << ">{\\setlength{\\hsize}{.1\\hsize}}Y|} " << std::endl << std::endl;

  // Place caption on table and add index entry
  std::string device_description_lc = device_description;
  std::transform(device_description_lc.begin(), device_description_lc.end(), device_description_lc.begin(), (int (*)(int)) std::tolower);

  os << std::endl
     << "\\caption{" << device_description << " " << (type ? "Device Model Parameters" : "Device Instance Parameters") << "."
     << "\\label{" << device_name << "_" << device_level << (type ? "_Device Model_Params" : "_Device Instance_Params") << "}}" << std::endl
     << "\\index{" << device_description_lc << "!" << (type ? "device model parameters" : "device instance parameters") << "} \\\\ \\hline " << std::endl;

  os << std::endl
     << "\\rowcolor{XyceDarkBlue}\\color{white}\\bf Parameter &" << std::endl
     << "\\color{white}\\bf Description & \\color{white}\\bf Units &" << std::endl
     << "\\color{white}\\bf Default \\endhead" << std::endl;

  // If catagorical listing, group by category
  if (format == OutputMode::DOC_CAT) {
    const ParametricData<void>::ParameterMap &parameter_map = parameters.getMap();
    for (int category = CAT_NONE; category != CAT_MAX; ++category) {
      const std::string &header = categoryName(category);
      bool header_printed = false;

      for (ParametricData<void>::ParameterMap::const_iterator it = parameter_map.begin(); it != parameter_map.end(); ++it) {
        std::string parameter_name = (*it).first;

        const Descriptor &descriptor = *(*it).second;

        std::string parameter_description;
        ParameterUnit parameter_unit = U_INVALID;
        ParameterCategory parameter_category = CAT_UNKNOWN;
        unitDescription((*it).first, parameter_map, descriptor, parameter_unit, parameter_category, parameter_description);

        if (parameter_category == category) {
          if (!header.empty() && !header_printed) {
            header_printed = true;

            os << std::endl
               << "\\multicolumn{4}{|c|}{\\color{XyceDarkBlue}\\em\\bfseries " << header << "}" << std::endl
               << "\\\\ \\hline\\hline" << std::endl;
          }

          documentParameter(os, parameter_name, parameter_unit, parameter_description, descriptor);
        }
      }
    }
  }
  else {
    const ParametricData<void>::ParameterMap &parameter_map = parameters.getMap();
    for (ParametricData<void>::ParameterMap::const_iterator it = parameter_map.begin(); it != parameter_map.end(); ++it) {
      std::string parameter_name = (*it).first;
      const Descriptor &descriptor = *(*it).second;

      std::string parameter_description;
      ParameterUnit parameter_unit = U_INVALID;
      ParameterCategory parameter_category = CAT_UNKNOWN;
      unitDescription((*it).first, parameter_map, descriptor, parameter_unit, parameter_category, parameter_description);

      documentParameter(os, parameter_name, parameter_unit, parameter_description, descriptor);
    }
  }

  os << std::endl  << "\\end{longtable}" << std::endl;

  return os;
}


std::ostream &outputParams(std::ostream &os, const ParametricData<void> &parametric_data, OutputMode::Mode mode)
{
  os << "Parameters" << Util::push << std::endl;
  outputParameterMap(os, parametric_data.getMap());
  os << Util::pop << std::endl;

  return os;

}


std::ostream &outputParameterMap(std::ostream &os, const ParametricData<void>::ParameterMap &parameter_map)
{
  for (ParametricData<void>::ParameterMap::const_iterator it = parameter_map.begin(); it != parameter_map.end(); ++it) {
    os << (*it).first << ", ";

    outputDescriptor(os, *(*it).second);
  }

  return os;
}

std::ostream &outputDescriptor(std::ostream &os, const Descriptor &descriptor)
{
  if (&descriptor.getEntry()) {
    printTypeName(os, descriptor.getEntry().type());
  }

  if (!descriptor.getCompositeParametricData<void>()) {
    os << ", default ";
    descriptor.getEntry().print(os);
    if (descriptor.getOriginalValueIndex() >= 0)
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

} // namespace Device
} // namespace Xyce

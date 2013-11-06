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
// Filename       : $RCSfile: N_DEV_Pars.C,v $
//
// Purpose        : 
//
// Special Notes  : 
//
// Creator        : David Baur
//
// Creation Date  : 3/28/2013
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.4.2.3 $
//
// Revision Date  : $Date: 2013/10/03 17:23:38 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <ostream>
#include <string>

#include <N_DEV_Pars.h>
#include <N_UTL_IndentStreamBuf.h>
#include <N_UTL_Demangle.h>

namespace Xyce {
namespace Device {

typedef std::map<std::pair<ParameterBase *, int>, double>  OriginalValueMap;
typedef std::map<std::pair<ParameterBase *, int>, bool>  GivenValueMap;

OriginalValueMap s_originalValueMap;                  ///< Map from device entity and original value index to original value
GivenValueMap s_givenValueMap;                        ///< Map from device entity and serial number to value given flag


/** 
 * Retrieve a parameter's original value
 *
 * @param entity device entity holding parameter
 * @param original_value_index index of original value
 *
 * @return original value of the parameter
 *
 * @date   Tue Aug  6 13:51:16 2013
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
 */
double getOriginalValue(ParameterBase *entity, int original_value_index)
{
  return s_originalValueMap[OriginalValueMap::key_type(entity, original_value_index)];
}


/** 
 * Set a parameter's original value
 *
 * @param entity device entity holding parameter
 * @param original_value_index index of original value
 * @param value value to be stored
 *
 * @date   Tue Aug  6 13:53:29 2013
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
 */
void setOriginalValue(ParameterBase *entity, int original_value_index, double value)
{
  s_originalValueMap[OriginalValueMap::key_type(entity, original_value_index)] = value;
}


/** 
 * Return true if a value was provided for the device
 *
 * @param entity device entity holding parameter
 * @param serial_number serial number of parameter
 *
 * @return true if a value was provided for the device
 *
 * @date   Tue Aug  6 13:54:25 2013
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
 */
bool wasValueGiven(ParameterBase *entity, int serial_number)
{
  return s_givenValueMap[GivenValueMap::key_type(entity, serial_number)];
}


/** 
 * Set the given value state of a parameter
 *
 * @param entity device entity holding parameter
 * @param serial_number serial number of parameter
 * @param value true if the value was given
 *
 * @date   Tue Aug  6 13:55:34 2013
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
 */
void setValueGiven(ParameterBase *entity, int serial_number, bool value)
{
  s_givenValueMap[GivenValueMap::key_type(entity, serial_number)] = value;
}

/** 
 * Adds an entry to the parameter name to descriptor map.
 *
 *
 * @param name Name or parameter to declare
 * @param descriptor Type information of parameter
 * @param parameter_data_class Typeinfo to display name of class on error
 * 
 * @date   Tue Aug  6 13:22:21 2013
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
 */
void ParametricData<void>::addDescriptor(const std::string &name, Descriptor *descriptor, const std::type_info &parameter_data_class)
{
  std::pair<ParameterMap::iterator, bool> result = map_.insert(ParameterMap::value_type(name, descriptor));

  if (!result.second) {
    std::ostringstream oss;

    oss << "Parameter " << name << " already added to class " << demangle(parameter_data_class.name());
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, oss.str());
  }
}

/** 
 * Report casting error when attempting to cast from from_type to to_type
 *
 * 
 * @param from_type Typeinfo casting from
 * @param to_type Typeinfo casting to
 *
 * @date   Tue Aug  6 13:40:56 2013
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
 */
void typeMismatch(const std::type_info &from_type, const::type_info &to_type)
{
  std::ostringstream oss;

  oss << "Attempting to cast parameter of type " << demangle(from_type.name()) << " to type " << demangle(to_type.name());

  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, oss.str());
}

/** 
 * Report error if both MIN_CAP and MIN_RES have been specified.
 * 
 * @param name Parameter name
 * @param expr_access Parameter expr access to verify
 * @param parameter_data_class Typeinfo to display name of class on error
 *
 * @author Dave Shirley, PSSI
 * @date   03/03/06
 */
void checkExprAccess(const std::string &name, ParameterType::ExprAccess &expr_access, const std::type_info &parameter_data_class)
{
  if ((expr_access & ParameterType::MIN_CAP) && (expr_access & ParameterType::MIN_RES))
  {
    std::ostringstream oss;

    oss << "Attempt to set MIN_CAP and MIN_RES on ParameterType::ExprAccess for parameter " << name << " in class " << parameter_data_class.name();

    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, oss.str());
  }
}

} // namespace Device
} // namespace Xyce

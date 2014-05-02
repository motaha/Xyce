//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2014 Sandia Corporation
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
// Filename       : $RCSfile: N_DEV_DeviceEntity.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 04/03/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.223.2.3 $
//
// Revision Date  : $Date: 2014/03/12 16:50:27 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <string>
#include <iostream>
#include <map>

#include <N_DEV_fwd.h>
#include <N_DEV_CompositeParam.h>
#include <N_DEV_Const.h>
#include <N_DEV_DeviceEntity.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_Message.h>
#include <N_DEV_Param.h>
#include <N_DEV_SolverState.h>
#include <N_UTL_BreakPoint.h>
#include <N_UTL_Expression.h>


namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::DeviceEntity
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 12/13/04
//-----------------------------------------------------------------------------
DeviceEntity::DeviceEntity(
  const char * const            entity_type,
  const std::string &           device_name,
  ParametricData<void> &        parametric_data,
  const SolverState &           solver_state,
  const DeviceOptions &         device_options,
  const std::string &           netlist_path,
  int                           netlist_line)
  : entityType_(entity_type),
    name_(device_name),
    parametricData_(parametric_data),
    solState_(solver_state),
    devOptions_(device_options),
    defaultParamName_(),
    netlistLocation_(netlist_path, netlist_line)
{}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::~DeviceEntity
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/30/00
//-----------------------------------------------------------------------------
DeviceEntity::~DeviceEntity()
{
  std::vector<Depend>::iterator d = dependentParams.begin();
  std::vector<Depend>::iterator end = dependentParams.end();
  for ( ; d != end; ++d)
  {
    delete d->expr;
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::scaleParam
//
// Purpose : Scales the original value of the specified parameter by the specified value.
// The parameter is never specified by the user so errors are developer caused.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/21/04
//-----------------------------------------------------------------------------
bool DeviceEntity::scaleParam( const std::string & paramName, double val, double val0)
{
  ParameterMap::const_iterator p_i = getParameterMap().find(paramName);
  if (p_i == getParameterMap().end())
  {
    DevelFatal(*this).in("DeviceEntity::scaleParam") << "Unrecognized parameter " << paramName;
    return false;
  }

  const Descriptor &param = *(*p_i).second;
  if (!param.hasOriginalValueStored())
  {
    DevelFatal(*this).in("DeviceEntity::scaleParam") << "Original value not available for parameter " << paramName;
    return false;
  }

  if (!param.isType<double>())
  {
    DevelFatal(*this).in("DeviceEntity::scaleParam") << "Can scale only double parameters, parameter " << paramName << " is not double";
    return false;
  }

  // Scale the parameter
  setValue<double, DeviceEntity>(*this, param, Xyce::Device::getOriginalValue(*this, param.getSerialNumber())*val + val0*(1.0-val));

  if (param.hasGivenMember())
    param.setGiven(*this, true);

  Xyce::Device::setValueGiven(*this, param.getSerialNumber(), true);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::scaleParam
//
// Purpose       : Scales the specified parameter by a specified value.
// The parameter is never specified by the user so errors are developer caused.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/21/04
//-----------------------------------------------------------------------------
bool DeviceEntity::scaleParam( const std::string & paramName, double val)
{
  ParameterMap::const_iterator p_i = getParameterMap().find(paramName);
  if (p_i == getParameterMap().end())
  {
    DevelFatal(*this).in("DeviceEntity::scaleParam") << "Unrecognized parameter " << paramName;
    return false;
  }

  const Descriptor &param = *(*p_i).second;
  if (!param.hasOriginalValueStored())
  {
    DevelFatal(*this).in("DeviceEntity::scaleParam") << "Original value not available for parameter " << paramName;
    return false;
  }

  if (!param.isType<double>())
  {
    DevelFatal(*this).in("DeviceEntity::scaleParam") << "Can scale only double parameters, parameter " << paramName << " is not double";
    return false;
  }

  // Scale the parameter
  param.value<double>(*this) = Xyce::Device::getOriginalValue(*this, param.getSerialNumber())*val;

  if (param.hasGivenMember())
    param.setGiven(*this, true);

  Xyce::Device::setValueGiven(*this, param.getSerialNumber(), true);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::scaleDefaultParam
// Purpose       :
// The parameter is never specified by the user so errors are developer caused.
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 50/19/05
//-----------------------------------------------------------------------------
bool DeviceEntity::scaleDefaultParam(double val)
{
  if (defaultParamName_.empty())
  {
    DevelFatal(*this).in("DeviceEntity::scaleDefaultParam") << "Device " << getName() << " does not have a default parameter";
    return false;
  }

  return scaleParam(defaultParamName_, val);
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::setParam
//
// Purpose       : This function loops over the vector of parameters, and
//                 sets the specified one (if found) to a specified value.
//
// Special Notes : This is kind of tricky, b/c some parameters are actually
//                 deep inside other classes (like a source class, for
//                 example)
//
//                 This function always returns a true, b/c there are many
//                 instances, (for example running in parallel), where one
//                 could set a param that didn't exist locally on
//                 processor.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/25/03
//-----------------------------------------------------------------------------
bool DeviceEntity::setParam(const std::string & paramName, double val)
{
  ParameterMap::const_iterator p_i = getParameterMap().find(paramName);
  if (p_i == getParameterMap().end())
    return false;

  if (isTempParam(paramName))
    val += CONSTCtoK;

  const Descriptor &param = *(*p_i).second;

  if (param.isType<double>())
    param.value<double>(*this) = val;
  else if (param.isType<int>())
    param.value<int>(*this) = static_cast <int> (val);
  else if (param.isType<long>())
    param.value<long>(*this) =  static_cast <long> (val);
  else if (param.isType<bool>())
    param.value<bool>(*this) = (val != 0);
  else
    DevelFatal0(*this) << "Illegal type for parameter " << paramName;

  if (param.hasGivenMember())
    param.setGiven(*this, true);

  Xyce::Device::setValueGiven(*this, param.getSerialNumber(), true);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::getParam
//
// Purpose       : returns the value of the requested param.
//
// Special Notes : This  function currently assumes that the requested
//                 param is a double-precision number.
//
//                 Parameters are not case-dependent.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/25/03
//-----------------------------------------------------------------------------
bool DeviceEntity::getParam ( const std::string & paramName, double & result)
{
  double val = 0.0;
  bool found = false;

  ParameterMap::const_iterator p_i = getParameterMap().find(paramName);
  if (p_i != getParameterMap().end())
  {
    found = true;
    const Descriptor &param = *(*p_i).second;
    if (param.isType<double>())
      val = param.value<double>(*this);
    else if (param.isType<int>())
      val = static_cast <double> (param.value<int>(*this));
    else if (param.isType<long>())
      val = static_cast <double> (param.value<long>(*this));
    else if (param.isType<bool>())
    {
      if (param.value<bool>(*this))
        val = 1;
      else
        val = 0;
    }
    else
    {
      DevelFatal(*this).in("DeviceEntity::getParam") << "Illegal type for parameter " << paramName;
    }
    if (isTempParam(paramName))
      val -= CONSTCtoK;
  }
  else
  {
    // If not recognized, just do nothing
  }
  result = val;

  return found;
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::setDefaultParam
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/06/03
//-----------------------------------------------------------------------------
bool DeviceEntity::setDefaultParam (double val)
{
  if (defaultParamName_.empty())
  {
    DevelFatal(*this).in("DeviceEntity::setDefaultParam") << getName() << " does not have a default parameter";
  }

  return setParam(defaultParamName_, val);
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::getDefaultParam
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/06/03
//-----------------------------------------------------------------------------
double DeviceEntity::getDefaultParam ()
{
  if (defaultParamName_.empty())
  {
    return 0.0;
  }

  double result = 0.0;

  getParam(defaultParamName_, result);

  return result;
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::setDependentParameter
// Purpose       : Add expression, param pairs for future updates
// Special Notes :  This is an overloaded method, used instead of the old
//                  monolithic one.
// Scope         : protected
// Creator       : Tom Russo
// Creation Date : 6 Nov 07
//-----------------------------------------------------------------------------
double DeviceEntity::setDependentParameter (Util::Param & par,
                                            double *res,
                                            ParameterType::ExprAccess depend)

{
  Depend dependentParam;
  setDependentParameter(par, dependentParam, depend);

  dependentParam.resultU.result = res;
  dependentParam.vectorIndex = -1;
  dependentParams.push_back(dependentParam);

  double rval;
  dependentParam.expr->evaluateFunction (rval);
  dependentParam.expr->set_sim_time( getSolverState().currTime );

  return rval;
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::setDependentParameter
// Purpose       : Add expression, param pairs for future updates
// Special Notes :  This is an overloaded method, used instead of the old
//                  monolithic one, and is specifically to set an element
//                  of a double vector.
// Scope         : protected
// Creator       : Tom Russo
// Creation Date : 6 Nov 07
//-----------------------------------------------------------------------------
double DeviceEntity::setDependentParameter (Util::Param & par,
                                            std::vector<double> *res,
                                            int ind,
                                            ParameterType::ExprAccess depend)

{
  Depend dependentParam;
  setDependentParameter(par,dependentParam, depend);

  dependentParam.resultU.resVec = res;
  dependentParam.vectorIndex = ind;
  dependentParams.push_back(dependentParam);

  double rval;
  dependentParam.expr->evaluateFunction (rval);
  dependentParam.expr->set_sim_time( getSolverState().currTime );

  return rval;
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::setDependentParameter
// Purpose       : Add expression, param pairs for future updates
// Special Notes :  This is a utility version, used by overloaded methods.
// Scope         : protected
// Creator       : Dave Shirley, PSSI
// Creation Date : 11/18/04
//-----------------------------------------------------------------------------

void DeviceEntity::setDependentParameter (Util::Param & par,
                                          Depend & dependentParam,
                                          ParameterType::ExprAccess depend)
{
  std::vector<std::string> instances, leads, names, variables;

  dependentParam.name = par.tag();
  if (isTempParam(par.tag()))
  {
    dependentParam.expr = new Util::Expression ("(" + par.stringValue() + ")+CONSTCtoK");
    dependentParam.expr->make_constant (std::string("CONSTCTOK"), CONSTCtoK);
  }
  else
  {
    dependentParam.expr = new Util::Expression (par.getValue<Util::Expression>());
  }

  names.clear();
  leads.clear();
  instances.clear();
  variables.clear();

  dependentParam.expr->get_names(XEXP_NODE, names);
  dependentParam.expr->get_names(XEXP_LEAD, leads);
  dependentParam.expr->get_names(XEXP_INSTANCE, instances);
  dependentParam.expr->get_names(XEXP_VARIABLE, variables);

  //std::vector<std::string>::iterator s;
  std::vector<std::string>::iterator iterS;

  if (!(depend & ParameterType::SOLN_DEP))
  {
    if (names.size() > 0 || instances.size() > 0)
    {
      UserError0(*this) << "Parameter " << par.tag() << " is not allowed to depend on voltage/current values";
      return;
    }
    if (depend & ParameterType::NO_DEP)
    {
      if (dependentParam.expr->get_num(XEXP_SPECIAL) > 0)
      {
        UserError0(*this) << "Parameter " << par.tag() << " is not allowed to depend on time";
        return;
      }
    }
  }

  if (leads.size() > 0)
  {
    char type;
    int index;
    for (std::vector<std::string>::const_iterator n_i=leads.begin(); n_i != leads.end(); ++n_i)
    {
      index = n_i->find_last_of(":");
      if (index == std::string::npos )
        type = (*n_i)[0];
      else
        type = (*n_i)[index+1];

      if (type != 'B' && type != 'E' && type != 'H')
      {
        UserError(*this) << "Illegal use of lead current specification in expression '" << dependentParam.expr->get_expression()
                         << "' in parameter " << par.tag();
      }
    }
    names.insert( names.end(), leads.begin(), leads.end() );
  }

  names.insert( names.end(), instances.begin(), instances.end() );

  dependentParam.lo_var = expVarNames.size();
  dependentParam.n_vars = names.size();
  dependentParam.vals.resize(dependentParam.n_vars);
  int expVarLen = dependentParam.lo_var+dependentParam.n_vars;
  expVarGIDs.resize(expVarLen);
  expVarLIDs.resize(expVarLen);
  expVarVals.resize(expVarLen);

  if (!variables.empty())
    names.insert( names.end(), variables.begin(), variables.end() );

  if ( !names.empty() )
  {
    // Order the names in the expression so that it agrees with the order
    // in names.
    dependentParam.expr->order_names( names );
  }
  for (int i=0 ; i<dependentParam.n_vars ; ++i)
    expVarNames.push_back(names[i]);

  if (dependentParam.n_vars > 0)
  {
    std::vector<double> zeros;
    zeros.resize(dependentParam.n_vars);
    for (int i=0 ; i<dependentParam.n_vars ; ++i)
      zeros[i] = 0;
    dependentParam.expr->set_vars(zeros);
  }

  dependentParam.global_params.clear();
  if (!variables.empty())
  {
    for (iterS=variables.begin() ; iterS!=variables.end() ; ++iterS)
    {
      if (getSolverState().global_params.find(*iterS) == getSolverState().global_params.end())
      {
        UserError0(*this) << "Global parameter " << *iterS << " not found";
      }
      else {
        dependentParam.expr->set_var(*iterS, getSolverState().global_params[*iterS]);
        dependentParam.global_params.push_back(*iterS);
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::updateDependentParameters
// Purpose       : Update values of parameters defined as expressions
// Special Notes :
// Scope         : protected
// Creator       : Dave Shirley, PSSI
// Creation Date : 03/15/05
//-----------------------------------------------------------------------------
bool DeviceEntity::updateDependentParameters(N_LAS_Vector & vars)
{
  std::vector<Depend>::iterator dpIter = dependentParams.begin();
  std::vector<Depend>::iterator end = dependentParams.end();
  double rval(0.0);
  int i, hi;
  bool changed = false;

  for ( ; dpIter != end ; ++dpIter)
  {
    if (dpIter->expr->set_sim_time( getSolverState().currTime ))
      changed = true;
    eVarVals.resize(dpIter->n_vars);
    if (dpIter->n_vars > 0)
    {
      hi = dpIter->lo_var+dpIter->n_vars;
      for (i=dpIter->lo_var ; i<hi ; ++i)
      {
        expVarVals[i] =vars[expVarLIDs[i]];
        eVarVals[i-dpIter->lo_var] = expVarVals[i];
      }
      if (dpIter->expr->set_vars(eVarVals))
        changed = true;
    }
    dpIter->expr->evaluateFunction (rval);
    if (dpIter->vectorIndex==-1)
      *(dpIter->resultU.result) = rval;
    else
      (*(dpIter->resultU.resVec))[dpIter->vectorIndex] = rval;
  }

  return changed;
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::updateGlobalParameters
// Purpose       : Update values of global parameters in expressions
// Special Notes :
// Scope         : protected
// Creator       : Dave Shirley, PSSI
// Creation Date : 11/17/05
//-----------------------------------------------------------------------------
bool DeviceEntity::updateGlobalParameters(std::map<std::string,double> & global_map)
{
  std::vector<Depend>::iterator dpIter = dependentParams.begin();
  std::vector<Depend>::iterator end = dependentParams.end();
  double rval;
  int i, hi;
  bool changed = false;

  for ( ; dpIter != end ; ++dpIter)
  {
    if (!dpIter->global_params.empty())
    {
      std::vector<std::string>::iterator gp=dpIter->global_params.begin();
      std::vector<std::string>::iterator gend=dpIter->global_params.end();
      for ( ; gp != gend; ++gp)
      {
        if (global_map.find(*gp) == global_map.end())
        {
          DevelFatal(*this).in("DeviceEntity::updateGlobalParameters") << "Failed to find global parameter " << *gp;
        }
        if (dpIter->expr->set_var(*gp, global_map[*gp]))
          changed = true;
      }
    }
  }

  return changed;
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::updateDependentParameters
// Purpose       : Update values of parameters defined as expressions
// Special Notes :
// Scope         : protected
// Creator       : Dave Shirley, PSSI
// Creation Date : 11/18/04
//-----------------------------------------------------------------------------
bool DeviceEntity::updateDependentParameters()
{
  double rval;
  bool changed = false;

  std::vector<Depend>::iterator dpIter = dependentParams.begin();
  std::vector<Depend>::iterator end = dependentParams.end();
  for ( ; dpIter != end; ++dpIter)
  {
    if (dpIter->expr->set_sim_time( getSolverState().currTime ))
      changed = true;
    dpIter->expr->evaluateFunction (rval);
    if (dpIter->vectorIndex == -1)
      *(dpIter->resultU.result) = rval;
    else
      (*(dpIter->resultU.resVec))[dpIter->vectorIndex] = rval;
  }

  return changed;
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::updateDependentParameters
// Purpose       : Update values of parameters defined as expressions with a
//                 specified temperature
// Special Notes :
// Scope         : protected
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/11/06
//-----------------------------------------------------------------------------
bool DeviceEntity::updateDependentParameters(double tempIn)
{
  double rval;
  bool changed = false;

  std::vector<Depend>::iterator dpIter = dependentParams.begin();
  std::vector<Depend>::iterator end = dependentParams.end();
  for ( ; dpIter != end; ++dpIter)
  {
    if (dpIter->expr->set_sim_time( getSolverState().currTime ) || dpIter->expr->set_temp(tempIn))
      changed = true;
    dpIter->expr->evaluateFunction (rval);
    if (dpIter->vectorIndex == -1)
      *(dpIter->resultU.result) = rval;
    else
      (*(dpIter->resultU.resVec))[dpIter->vectorIndex] = rval;

  }

  return changed;
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::getParamBreakpoints
// Purpose       : Add breakpoints caused by discontinuities in computed params
// Special Notes :
// Scope         : protected
// Creator       : Dave Shirley, PSSI
// Creation Date : 11/18/04
//-----------------------------------------------------------------------------
bool DeviceEntity::getParamBreakpoints( std::vector<Util::BreakPoint> & breakPointTimes )
{
  double bTime;

  std::vector<Depend>::iterator dpIter = dependentParams.begin();
  std::vector<Depend>::iterator end = dependentParams.end();
  for ( ; dpIter != end; ++dpIter)
  {
    bTime = dpIter->expr->get_break_time();
    if (bTime > getSolverState().currTime)
      breakPointTimes.push_back(bTime);
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::given
// Purpose       : Return whether param was given
// Special Notes :
// Scope         : protected
// Creator       : Dave Shirley, PSSI
// Creation Date : 11/23/04
//-----------------------------------------------------------------------------
bool DeviceEntity::given( const std::string & parameter_name ) const
{
  ParameterMap::const_iterator it = getParameterMap().find(parameter_name);
  if (it == getParameterMap().end())
    DevelFatal0(*this).in("DeviceEntity::given") << "Unrecognized parameter " << parameter_name;

  return Xyce::Device::wasValueGiven(*this, (*it).second->getSerialNumber());
}

//-----------------------------------------------------------------------------
// Function      : setParameters
// Purpose       : Set parameters according to a vector of params.  Used to
//                 set instance or model parameter sets to netlist values
// Special Notes :
// Scope         : protected
// Creator       : Dave Shirley, PSSI
// Creation Date : 11/20/04
//-----------------------------------------------------------------------------
void setParameters(DeviceEntity &entity, std::vector<Param>::const_iterator begin, std::vector<Param>::const_iterator end, const DeviceOptions &device_options)
{
  std::vector<std::string> composite_name_list;
  std::map<std::string, std::vector<CompositeParam *>, LessNoCase> composite_parameter_map;

  if (DEBUG_DEVICE && device_options.debugLevel > 0)
  {
    Xyce::dout() << std::endl << "In DeviceEntity::setParams, for " << entity.getEntityType()
                 << ": " << entity.getName() << " parameters are:" << std::endl;
    for (std::vector<Param>::const_iterator it = begin; it != end; ++it)
    {
      const Param &param = *it;
      Xyce::dout() << "Param = " << param.tag() << ", Type = ";
      int tmpType = param.getType();
      switch (tmpType)
      {
        case Util::STR:
          Xyce::dout() << "STR";
          break;
        case Util::DBLE:
          Xyce::dout() << "DBLE";
          break;
        case Util::INT:
          Xyce::dout() << "INT";
          break;
        case Util::LNG:
          Xyce::dout() << "LNG";
          break;
        case Util::EXPR:
          Xyce::dout() << "EXPR";
          break;
        case Util::BOOL:
          Xyce::dout() << "BOOL";
          break;
        case Util::STR_VEC:
          Xyce::dout() << "STR_VEC";
          break;
        case Util::INT_VEC:
          Xyce::dout() << "INT_VEC";
          break;
        case Util::DBLE_VEC:
          Xyce::dout() << "DBLE_VEC";
          break;
        case Util::DBLE_VEC_IND:
          Xyce::dout() << "DBLE_VEC_IND";
          break;
        case Util::COMPOSITE:
          Xyce::dout() << "COMPOSITE";
          break;
        default:
          Xyce::dout() << "Unknown";
      }
      Xyce::dout() << ", Value = " << param.stringValue();

      if (param.given())
      {
        Xyce::dout() << "  given=TRUE";
      }
      else
      {
        Xyce::dout() << "  given=FALSE";
      }

      if (param.default_val())
      {
        Xyce::dout() << "  default=TRUE" << std::endl;
      }
      else
      {
        Xyce::dout() << "  default=FALSE" << std::endl;
      }
    }
    Xyce::dout() << std::endl;
  }

  for (std::vector<Param>::const_iterator param_it = begin; param_it != end; ++param_it)
  {
    Param &param = const_cast<Param &>(*param_it);

    const std::string &tag = param.tag();

    // Is this parameter in the Entity?
    ParameterMap::const_iterator entity_parameter_it = entity.getParameterMap().find(tag);
    if (entity_parameter_it != entity.getParameterMap().end())
    {
      const Descriptor &descriptor = *(*entity_parameter_it).second;
      if (descriptor.hasGivenMember())
      {
        if (param.given())
        {
          descriptor.setGiven(entity, true);
        }
        else if (descriptor.getGiven(entity))
        {
          continue;
        }
      }

      Xyce::Device::setValueGiven(entity, descriptor.getSerialNumber(), param.given());
      if (param.given() || param.default_val())
      {
        if ( param.getType() == Util::EXPR )
        {
          if (descriptor.isType<double>())
          {
            double val = entity.setDependentParameter(param, &(descriptor.value<double>(entity)), descriptor.getExpressionAccess());
            param.setVal(val);
          }
          else if (descriptor.isType<std::vector<double> >())
          {
            int ind = (descriptor.value<std::vector<double> >(entity)).size();
            double val = entity.setDependentParameter (param, &(descriptor.value<std::vector<double> >(entity)), ind, descriptor.getExpressionAccess());
            (descriptor.value<std::vector<double> >(entity)).push_back(val);
          }
          else
          {
            DevelFatal(entity).in("DeviceEntity::setParams") << "Non double param " <<  tag << " cannot be set to expression";
          }
        }
        else
        {
          if (descriptor.isType<double>())
          {
            if (!param.isNumeric())
            {
              UserFatal(entity) << "Cannot convert parameter " << tag <<  " to a numeric value from " << param.stringValue();
            }

            descriptor.value<double>(entity) = param.getImmutableValue<double>();
            if (isTempParam(tag))
            {
              descriptor.value<double>(entity) += CONSTCtoK;
            }
            if (descriptor.hasOriginalValueStored())
            {
              Xyce::Device::setOriginalValue(entity, descriptor.getSerialNumber(), descriptor.value<double>(entity));
            }
          }
          else if (descriptor.isType<std::string>())
          {
            descriptor.value<std::string>(entity) = param.stringValue();
          }
          else if (descriptor.isType<int>())
          {
            if (!param.isInteger())
            {
              UserFatal(entity) << "Cannot convert parameter " << tag << " to an integer value from " << param.stringValue();
            }
            descriptor.value<int>(entity) = param.getImmutableValue<int>();
            if (descriptor.hasOriginalValueStored())
            {
              Xyce::Device::setOriginalValue(entity, descriptor.getSerialNumber(), static_cast<double> (descriptor.value<int>(entity)));
            }
          }
          else if (descriptor.isType<long>())
          {
            if (!param.isInteger())
            {
              UserFatal(entity) << "Cannot convert parameter " << tag << " to an integer value from " << param.stringValue();
            }
            descriptor.value<long>(entity) = param.getImmutableValue<long>();
            if (descriptor.hasOriginalValueStored())
            {
              Xyce::Device::setOriginalValue(entity, descriptor.getSerialNumber(), static_cast<double> (descriptor.value<long>(entity)));
            }
          }
          else if (descriptor.isType<bool>())
          {
            if (!param.isBool())
            {
              UserFatal(entity) << "Cannot convert parameter " << tag << " to a logical value from " << param.stringValue();
            }
            descriptor.value<bool>(entity) = param.getImmutableValue<bool>();
            if (descriptor.hasOriginalValueStored())
            {
              if (descriptor.value<bool>(entity))
              {
                Xyce::Device::setOriginalValue(entity, descriptor.getSerialNumber(), 1.0);
              }
              else
              {
                Xyce::Device::setOriginalValue(entity, descriptor.getSerialNumber(), 0.0);
              }
            }
          }
          else if (descriptor.isType<std::vector<int> >())
          {
            if (param.getType() == Util::INT_VEC)
            {
              (descriptor.value<std::vector<int> >(entity)) = param.getValue<std::vector<int> >();
            }
            else if (param.getType() == Util::INT)
            {
              (descriptor.value<std::vector<int> >(entity)).push_back(param.getImmutableValue<int>());
            }
          }
          else if (descriptor.isType<std::vector<double> >())
          {
            if (param.getType() == Util::DBLE_VEC)
            {
              (descriptor.value<std::vector<double> >(entity)) = param.getValue<std::vector<double> >();
            }
            else if (param.getType() == Util::DBLE)
            {
              (descriptor.value<std::vector<double> >(entity)).push_back(param.getImmutableValue<double>());
            }
          }
          else if (descriptor.isType<std::vector<std::string> >())
          {
            if (param.getType() == Util::STR_VEC)
            {
              (descriptor.value<std::vector<std::string> >(entity)) = param.getValue<std::vector<std::string> >();
            }
            else if (param.getType() == Util::STR)
            {
              (descriptor.value<std::vector<std::string> >(entity)).push_back(param.stringValue());
            }
          }
          else if (descriptor.isComposite())
          {
            composite_parameter_map[tag].clear();
            composite_name_list.push_back(tag);

            // Note: ERK.  This push-back is done to process the base-param tag of a vector composite.
            // For example, if the composite parameters are things like REGION0.XWIDTH, where REGION
            // is the base parameter tag, 0 is the index, and XWIDTH is the subcomponent, the
            // tag that should be pushed back is tagES=REGION.
            //
            // Note: ERK:  This function seems to implicitly rely on the base parameter always preceeding
            // subcomponent parameters.  So REGION (alone) should preceed REGION0.XWIDTH in  the
            // STL vector params that is passed into this function.  If it doesn't then vc_stat will
            // stay false and a fatal error will get thrown.
            if (DEBUG_DEVICE && device_options.debugLevel > 0)
            {
              Xyce::dout() << "pushing back composite " << tag << std::endl;
            }
          }
          else
          {
            DevelFatal(entity).in("DeviceEntity::setParams") << "Unknown type";
          }
        }
      }
    }

    // Is it a vector (why do nothing?)
    else if (param.stringValue() == "VECTOR")
    {
    }

    // Must be a composite
    else
    {
      bool vc_stat = false;

      std::string::size_type dot = tag.find_first_of('.');
      if (dot != std::string::npos)
      {
        for (std::vector<std::string>::const_iterator it = composite_name_list.begin(); it != composite_name_list.end(); ++it)
        {
          const std::string &composite_name = *it;

          if (tag.find(composite_name) == 0) // Tag starts with the composite parameter name
          {
            std::string param_name(tag.begin() + dot + 1, tag.end());

            vc_stat = true;

            int n = 0;
            {
              std::istringstream is(std::string(tag.begin() + composite_name.size(), tag.begin() + dot));
              is >> n;
            }

            if (param_name == "NAME")
            {
              if (n != composite_parameter_map[composite_name].size())
              {
                DevelFatal(entity).in("DeviceEntity::setParams") << "Error filling 'NAME' vector param " <<  composite_name;
              }

              std::string name = param.stringValue();
              CompositeParam *composite = entity.constructComposite(composite_name, name);
              composite_parameter_map[composite_name].push_back(composite);
              setDefaultParameters(*composite, composite->getParameterMap().begin(), composite->getParameterMap().end(), device_options);
            }
            else
            {
              if (n >= composite_parameter_map[composite_name].size())
              {
                DevelFatal(entity).in("DeviceEntity::setParams") << "Internal error in vector-composite, 'NAME' must come first " << composite_name;
              }
            }
            Xyce::Device::setParameters(*composite_parameter_map[composite_name][n], param_name, param);
          }
        }
      }
      if (!vc_stat)
      {
        UserError(entity) << "Undefined parameter " << tag << ", this parameter is in metadata, but not recognized in constructor";
      }
    }
  }

  if (!composite_name_list.empty())
  {
    for (std::vector<std::string>::const_iterator name_it = composite_name_list.begin(); name_it != composite_name_list.end(); ++name_it)
    {
      const std::string &vcs = *name_it;

      for (std::vector<CompositeParam *>::iterator it =  composite_parameter_map[vcs].begin(); it != composite_parameter_map[vcs].end(); ++it)
      {
        (*it)->processParams();
      }
    }
  }
}


//-----------------------------------------------------------------------------
// Function      : setParameters
// Purpose       : Set parameter according to input Param.  Used to
//                 set instance or model parameter sets to netlist values
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 05/06/05
//-----------------------------------------------------------------------------

void setParameters(CompositeParam &composite_param, const std::string & pName, const Param & ndParam )
{
  ParameterMap::const_iterator p_i = composite_param.getParameterMap().find(pName);
  if (p_i != composite_param.getParameterMap().end())
  {
    const Descriptor &p = *(*p_i).second;
    if (p.hasGivenMember())
    {
      if (ndParam.given())
        p.setGiven(composite_param, true);
      else if (p.getGiven(composite_param))
        return;
    }
    Xyce::Device::setValueGiven(composite_param, p.getSerialNumber(), ndParam.given());
    if (ndParam.given() || ndParam.default_val())
    {
      if (p.isType<double>())
      {
        p.value<double>(composite_param) = ndParam.getImmutableValue<double>();
        if (isTempParam(pName))
          p.value<double>(composite_param) += CONSTCtoK;
        if (p.hasOriginalValueStored())
          Xyce::Device::setOriginalValue(composite_param, p.getSerialNumber(), p.value<double>(composite_param));
      }
      else if (p.isType<std::string>())
      {
        p.value<std::string>(composite_param) = ndParam.stringValue();
      }
      else if (p.isType<int>())
      {
        p.value<int>(composite_param) = ndParam.getImmutableValue<int>();
        if (p.hasOriginalValueStored())
          Xyce::Device::setOriginalValue(composite_param, p.getSerialNumber(), static_cast<double>(p.value<int>(composite_param)));
      }
      else if (p.isType<long>())
      {
        p.value<long>(composite_param) = ndParam.getImmutableValue<long>();
        if (p.hasOriginalValueStored())
          Xyce::Device::setOriginalValue(composite_param, p.getSerialNumber(), static_cast<double>(p.value<long>(composite_param)));
      }
      else if (p.isType<bool>())
      {
        p.value<bool>(composite_param) = (ndParam.getImmutableValue<double>() != 0.0);
        if (p.hasOriginalValueStored())
        {
          if (p.value<bool>(composite_param))
            Xyce::Device::setOriginalValue(composite_param, p.getSerialNumber(), 1.0);
          else
            Xyce::Device::setOriginalValue(composite_param, p.getSerialNumber(), 0.0);
        }
      }
      else if (p.isType<std::vector<double> >())
      {
       (p.value<std::vector<double> >(composite_param)).push_back(ndParam.getImmutableValue<double>());
      }
      else if (p.isType<std::vector<std::string> >())
      {
        p.value<std::vector<std::string> >(composite_param).push_back(ndParam.stringValue());
      }
      else
      {
        Report::DevelFatal().in("CompositeParam::setParams") << "Unknown parameter type for " << pName;
      }
    }
  }
  else
  {
    Report::DevelFatal().in("CompositeParam::setParams") << "Undefined parameter " <<  pName << std::endl
                                                       << "This parameter is in metadata, but not recognized in constructor";
  }
}


void populateParams(const ParameterMap &parameter_map, std::vector<Param> &param_list, CompositeParamMap &composite_param_map)
{
  for (ParameterMap::const_iterator it = parameter_map.begin(); it != parameter_map.end(); ++it)
  {
    const Descriptor &param = *(*it).second;

    if (param.isType<double>())
    {
      if (param.getVec() == 0)
      {
        double val;
        if (isTempParam((*it).first))
          val = getDefaultValue<double>(param) - CONSTCtoK;
        else
          val = getDefaultValue<double>(param);
        param_list.push_back(Param((*it).first, val));
      }
      else if (param.getVec() > 0)
      {
        if (param.getVec() == 1)
        {
          // This converts parameters link IC1, IC2 to just IC and type vector
          std::string vPar((*it).first.substr(0, (*it).first.size()-1));
          param_list.push_back(Param(vPar, "VECTOR"));
        }
        // We will also output IC1, IC2 as type double so they can
        // be specified as individual elements
        // This allows TC=a, b to also be specified as TC1=a TC2=b
        double val = getDefaultValue<double>(param);
        param_list.push_back(Param((*it).first, val));
      }
    }
    else if (param.isType<bool>())
    {
      if (param.getVec() == 0)
        param_list.push_back(Param((*it).first, getDefaultValue<bool>(param)));
      else if (param.getVec() > 0)
      {
        if (param.getVec() == 1)
        {
          std::string vPar((*it).first.substr(0, (*it).first.size()-1));
          param_list.push_back(Param(vPar, "VECTOR"));
        }
        // We will also output IC1, IC2 as type double so they can
        // be specified as individual elements
        // This allows TC=a, b to also be specified as TC1=a TC2=b
        bool val = getDefaultValue<bool>(param);
        param_list.push_back(Param((*it).first, val));
      }
    }
    else if (param.isType<int>())
    {
      if (param.getVec() == 0)
        param_list.push_back(Param((*it).first, getDefaultValue<int>(param)));
      else if (param.getVec() > 0)
      {
        if (param.getVec() == 1)
        {
          std::string vPar((*it).first.substr(0, (*it).first.size()-1));
          param_list.push_back(Param(vPar, "VECTOR"));
        }
        int val = getDefaultValue<int>(param);
        param_list.push_back(Param((*it).first, val));
      }
    }
    else if (param.isType<std::string>())
    {
      if (param.getVec() == 0)
        param_list.push_back(Param((*it).first, getDefaultValue<std::string>(param)));
      else if (param.getVec() > 0)
      {
        if (param.getVec() == 1)
        {
          std::string vPar((*it).first.substr(0, (*it).first.size()-1));
          param_list.push_back(Param(vPar, "VECTOR"));
        }
        std::string val = getDefaultValue<std::string>(param);
        param_list.push_back(Param((*it).first, val));
      }
    }
    else if (param.isType<std::vector<std::string> >())
    {
      Param vc((*it).first, std::vector<std::string>()); // param.value<std::vector<std::string> >(entity));
      param_list.push_back(vc);
    }
    else if (param.isType<std::vector<double> >())
    {
      Param vc((*it).first, std::vector<double>()); // param.value<std::vector<double> >(entity));
      param_list.push_back(vc);
    }
    else if (param.isComposite())
    {
      Param vc2((*it).first, "VECTOR-COMPOSITE");
      vc2.setDefault(true);
      param_list.push_back(vc2);

      std::vector<Param> compositeParams;
      const ParametricData<CompositeParam> *c = param.getCompositeParametricData<CompositeParam>();

      if (c == 0)
      {
        Report::DevelFatal().in("populateParams") << "Vector-composite map for device type entity empty.";
      }

      // TODO: [DGB] I think when the Descriptor is refactored this will be clearer.  But this basically adds the
      //   type to the composite list with 'NAME' first.
      const ParametricData<CompositeParam> &d = *c;
      const ParameterMap &e = d.getMap();

      for (ParameterMap::const_iterator it4 = e.find("NAME"); it4 != e.end();) {
        const Descriptor &p = *(*it4).second;
        if (p.isType<double>())
          compositeParams.push_back(Param((*it4).first, getDefaultValue<double>(p)));
        else if (p.isType<bool>())
          compositeParams.push_back(Param((*it4).first, getDefaultValue<bool>(p)));
        else if (p.isType<int>())
          compositeParams.push_back(Param((*it4).first, getDefaultValue<int>(p)));
        else if (p.isType<std::string>())
          compositeParams.push_back(Param((*it4).first, getDefaultValue<std::string>(p)));
        if ((*it4).first == "NAME")
          it4 = e.begin();
        else
          it4++;
        if (it4 != e.end() && (*it4).first == "NAME")
          it4++;
      }

      composite_param_map[(*it).first] = compositeParams;
    }
    else
    {
//        Just skip these, like list of coupled inductors because not needed for metadata
//        std::string msg("DeviceEntity::getParams: Type not supported");
//        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
      Xyce::dout() << "In final else clause of DeviceEntity::getParams().";
      if( param.isType<std::vector<std::string> >() )
        Xyce::dout() << " type is STR_VEC ";
      if( param.isType<std::vector<double> >() )
        Xyce::dout() << " type is DBLE_VEC ";
      Xyce::dout() << it->first << " this item is NOT being added to default parameter list." << std::endl;
    }
  }
}


//-----------------------------------------------------------------------------
// Function      : printParameter
//
// Purpose       : This function finds a parameter in the par table, and then
//                 formats its value in a string.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/4/12
//-----------------------------------------------------------------------------
std::ostream &printParameter(std::ostream &os, const DeviceEntity &entity, const std::string &name, const Descriptor &param)
{
  if (param.isType<double>())
  {
    if (isTempParam(name))
    {
      os << param.value<double>(entity) - CONSTCtoK;
    }
    else
    {
      os << param.value<double>(entity);
    }
  }
  else if (param.isType<bool>())
  {
    os << param.value<bool>(entity);
  }
  else if (param.isType<int>())
  {
    os << param.value<int>(entity);
  }
  else if (param.isType<long>())
  {
    os << param.value<long>(entity);
  }
  else if (param.isType<std::string>())
  {
    os << param.value<std::string>(entity);
  }
  else if (param.isType<std::vector<std::string> >())
  {
    os << " (string vector) : ";
    os << "length = " << (param.value<std::vector<std::string> >(entity)).size();
    if ((param.value<std::vector<std::string> >(entity)).size() > 0)
    {
      os << " :";
      std::vector<std::string>::const_iterator iterStringVec = (param.value<std::vector<std::string> >(entity)).begin();
      for ( ; iterStringVec != (param.value<std::vector<std::string> >(entity)).end() ; ++iterStringVec)
      {
        os << "  " << *iterStringVec;
      }
    }
  }
  else if (param.isType<std::vector<int> >())
  {
    os << " (int vector) : ";
    os << "length = " << (param.value<std::vector<int> >(entity)).size();
    if ((param.value<std::vector<int> >(entity)).size() > 0)
    {
      os << " :";
      std::vector<int>::const_iterator v_i = (param.value<std::vector<int> >(entity)).begin();
      for ( ; v_i != (param.value<std::vector<int> >(entity)).end() ; ++v_i)
      {
        os << "  " << *v_i;
      }
    }
  }
  else if (param.isType<std::vector<double> >())
  {
    os << " (double vector) : ";
    os << "length = " << (param.value<std::vector<double> >(entity)).size();
    if ((param.value<std::vector<double> >(entity)).size() > 0)
    {
      os << " :";
      std::vector<double>::const_iterator iterDoubleVec = (param.value<std::vector<double> >(entity)).begin();
      for ( ; iterDoubleVec != (param.value<std::vector<double> >(entity)).end() ; ++iterDoubleVec)
      {
        os << "  " << *iterDoubleVec;
      }
    }
  }
  else if (param.isComposite())
  {
    os << " (composite) : " << (param.value<CompositeMap>(entity)).size();
    if ((param.value<CompositeMap>(entity)).size() > 0)
    {
      std::map<std::string, CompositeParam *>::const_iterator c_i = (param.value<CompositeMap>(entity)).begin();
      // (*c_i).second->printParams(os);
      printCompositeParameters(os, *(*c_i).second);
    }
  }
    
  return os;
}

std::ostream &printCompositeParameters(std::ostream &os, const CompositeParam &composite)
{
  for (ParameterMap::const_iterator it_parameter = composite.getParameterMap().find("NAME"); it_parameter != composite.getParameterMap().end() ; )
  {
    os << std::endl << "   " <<(*it_parameter).first << " ";
    const Descriptor &p = *(*it_parameter).second;
    if (p.isType<double>())
      os << "(double) : " << p.value<double>(composite);
    else if (p.isType<bool>())
      os << "(bool) : " << p.value<bool>(composite);
    else if (p.isType<int>())
      os << "(int) : " << p.value<int>(composite);
    else if (p.isType<long>())
      os << "(long) : " << p.value<long>(composite);
    else if (p.isType<std::string>())
      os << "(string) : " << p.value<std::string>(composite);
    else if (p.isType<std::vector<std::string> >())
    {
      const std::vector<std::string> &string_vector = p.value<std::vector<std::string> >(composite);

      for (std::vector<std::string>::const_iterator it = string_vector.begin(); it != string_vector.end(); ++it)
        os << "  " << *it;
    }
    else if (p.isType<std::vector<double> >())
    {
      const std::vector<double> &string_vector = p.value<std::vector<double> >(composite);

      for (std::vector<double>::const_iterator it = string_vector.begin(); it != string_vector.end(); ++it)
        os << "  " << *it;
    }

    if ((*it_parameter).first == "NAME")
      it_parameter = composite.getParameterMap().begin();
    else
      it_parameter++;

    if (it_parameter != composite.getParameterMap().end() && (*it_parameter).first == "NAME")
      it_parameter++;
  }

  return os;
}

// //-----------------------------------------------------------------------------
// // Function      : DeviceEntity::perturbSensParam
// // Purpose       :
// // Special Notes :
// // Scope         : public
// // Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// // Creation Date : 05/02/03
// //-----------------------------------------------------------------------------
// bool DeviceEntity::perturbSensParam (Param & ndParam)
// {
//   if ( !ndParam.isTimeDependent() )
//   {
//     double paramOld = ndParam.dVal();
//     getDeviceOptions().deviceSens_dp = getDeviceOptions().testJac_SqrtEta * (1.0 + fabs(ndParam.dVal()));
//     double paramNew = paramOld + getDeviceOptions().deviceSens_dp;

//     ndParam.setVal (paramNew);
//   }
//   else
//   {
//     UserWarning0(*this) << "Time dependent parameter (" << ndParam.uTag()
//                         << ") not allowed as a sensitivity parameter, continuing analysis without perturubing this parameter.";
//   }

//   if (Xyce::DEBUG_DEVICE && getDeviceOptions().debugLevel > 0)
//   {
//     Xyce::dout() << "DeviceEntity::perturbSensParams" << std::endl
//                  << "paramNew = " << paramNew <<std::endl
//                  << "paramOld = " << paramOld <<std::endl
//                  << "dp    = " << getDeviceOptions().deviceSens_dp << std::endl;
//   }

//   return true;
// }

} // namespace Device
} // namespace Xyce

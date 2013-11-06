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
// Revision Number: $Revision: 1.193.2.4 $
//
// Revision Date  : $Date: 2013/10/03 17:23:38 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------
#include <string>
#include <iostream>
#include <map>

// ----------   Xyce Includes   ----------
#include <N_DEV_fwd.h>
#include <N_DEV_DeviceEntity.h>
#include <N_DEV_DeviceModel.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_Units.h>
#include <N_DEV_Const.h>
#include <N_DEV_Param.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_CompositeParam.h>

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
  SolverState &         solver_state,
  DeviceOptions &       device_options,
  const std::string &   device_name,
  const std::string &   net_list_file_path,
  int                   net_list_file_line_number)
  : solState(solver_state),
    devOptions(device_options),
    commandLine(device_options.commandLine),
    paramNew(0.0),
    paramOld(0.0),
    name_(device_name),
    defaultParamName(""),
    netlistFileName_(net_list_file_path),
    lineNumber_(net_list_file_line_number)
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
  vector<Depend>::iterator d = dependentParams.begin();
  vector<Depend>::iterator end = dependentParams.end();
  for ( ; d != end; ++d)
  {
    delete d->expr;
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::processInstanceParams
// Purpose       : Generate error for virtual method
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 04/25/06
//-----------------------------------------------------------------------------
bool DeviceEntity::processInstanceParams(string s)
{
  string msg("DeviceEntity::processInstanceParams: method not found for: ");
  msg += getName();
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
  return true;
}

#if 0
//-----------------------------------------------------------------------------
// Function      : DeviceEntity::perturbSensParam
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/02/03
//-----------------------------------------------------------------------------
bool DeviceEntity::perturbSensParam (Param & ndParam)
{
  if ( !ndParam.isTimeDependent() )
  {
    paramOld = ndParam.dVal();
    getDeviceOptions().deviceSens_dp = getDeviceOptions().testJac_SqrtEta * (1.0 + fabs(ndParam.dVal()));
    paramNew = paramOld + getDeviceOptions().deviceSens_dp;

    ndParam.setVal (paramNew);
  }
  else
  {
    string msg("DeviceEntity::perturbSensParam\n");
    msg += "\t time dependent parameter (" + ndParam.uTag();
    msg += ") not allowed as a sensitivity parameter.\n";
    msg += "Continuing analysis without perturubing this parameter.";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_WARNING_0,msg);
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "DeviceEntity::perturbSensParams" << endl;
    cout << "paramNew = " << paramNew <<endl;
    cout << "paramOld = " << paramOld <<endl;
    cout << "dp    = " << getDeviceOptions().deviceSens_dp << endl;
  }
#endif
  return true;
}
#endif

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::scaleParam
//
// Purpose : Scales the original value of the specified parameter by the specified value.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/21/04
//-----------------------------------------------------------------------------
bool DeviceEntity::scaleParam( const std::string & paramName, double val, double val0)
{
  ExtendedString tmpName(paramName);
  tmpName.toUpper ();

  ParameterMap::const_iterator p_i = (*getPMap()).find(tmpName);
  if (p_i != (*getPMap()).end())
  {
    const Pars &param = static_cast<const Pars &>(*(*p_i).second);
    if (param.getOriginalValueIndex() >= 0)
    {
      if (param.isType<double>())
        setValue<double, DeviceEntity>(*this, param, getOriginalValue(this, param.getOriginalValueIndex())*val + val0*(1.0-val));
      else
      {
        string msg("DeviceEntity::scaleParam: can scale only double parameters: ");
        msg += paramName;
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
      }
    }
    else
    {
      string msg("DeviceEntity::scaleParam: original value not available for: ");
      msg += paramName;
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
    }
    if (param.hasGivenMember())
      param.setGiven(*this, true);
    setValueGiven(this, param.getSerialNumber(), true);
  }
  else
  {
    string msg("DeviceEntity::scaleParam: unrecognized param: ");
    msg += paramName;
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::scaleParam
//
// Purpose       : Scales the specified parameter by a specified value.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/21/04
//-----------------------------------------------------------------------------
bool DeviceEntity::scaleParam ( const std::string & paramName, double val)
{
  ParameterMap::const_iterator p_i = (*getPMap()).find(paramName);
  if (p_i != (*getPMap()).end())
  {
    const Pars &param = static_cast<const Pars &>(*(*p_i).second);
    if (param.getOriginalValueIndex() >= 0)
    {
      if (param.isType<double>())
        param.value<double>(*this) = getOriginalValue(this, param.getOriginalValueIndex())*val;
      else
      {
        string msg("DeviceEntity::scaleParam: can only scale double parameter: ");
        msg += paramName;
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
      }
    }
    else
    {
      string msg("DeviceEntity::scaleParam: original value not available for: ");
      msg += paramName;
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
    }
    if (param.hasGivenMember())
      param.setGiven(*this, true);
    setValueGiven(this, param.getSerialNumber(), true);
  }
  else
  {
    string msg("DeviceEntity::scaleParam: unrecognized param: ");
    msg += paramName;
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::scaleDefaultParam
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 50/19/05
//-----------------------------------------------------------------------------
bool DeviceEntity::scaleDefaultParam (double val)
{
  if (defaultParamName == "")
  {
    string msg("DeviceEntity::scaleDefaultParam. ");
    msg += getName() + " does not have a default parameter";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL, msg);
  }

  return scaleParam(defaultParamName, val);
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
bool DeviceEntity::setParam ( const string & paramName, double val)
{
  ExtendedString tmpName(paramName);
  tmpName.toUpper ();

  ParameterMap::const_iterator p_i = (*getPMap()).find(tmpName);
  if (p_i != (*getPMap()).end())
  {
    if (tmpName == "TEMP" || tmpName == "TNOM")
      val += CONSTCtoK;
    const Pars &param = static_cast<const Pars &>(*(*p_i).second);
    if (param.isType<double>())
      param.value<double>(*this) = val;
    else if (param.isType<int>())
      param.value<int>(*this) = static_cast <int> (val);
    else if (param.isType<long>())
      param.value<long>(*this) =  static_cast <long> (val);
    else if (param.isType<bool>())
      param.value<bool>(*this) = (val != 0);
    else
    {
      string msg("DeviceEntity::setParam: illegal type for param: ");
      msg += paramName;
      std::ostringstream oss;
      oss << "Error in " << netlistLocation() << "\n" << msg;
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, oss.str());
    }
    if (param.hasGivenMember())
      param.setGiven(*this, true);
    setValueGiven(this, param.getSerialNumber(), true);
  }
  else
  {
    return false;
  }

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
bool DeviceEntity::getParam ( const string & paramName, double & result)
{
  double val = 0.0;
  bool found = false;

  ExtendedString tmpName(paramName);
  tmpName.toUpper ();

  ParameterMap::const_iterator p_i = (*getPMap()).find(paramName);
  if (p_i != (*getPMap()).end())
  {
    found = true;
    const Pars &param = static_cast<const Pars &>(*(*p_i).second);
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
      string msg("DeviceEntity::getParam: illegal type for param: ");
      msg += paramName;
      std::ostringstream oss;
      oss << "Error in " << netlistLocation() << "\n" << msg;
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, oss.str());
    }
    if (tmpName == "TEMP" || tmpName == "TNOM")
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
  if (defaultParamName == "")
  {
    string msg("DeviceEntity::setDefaultParam. ");
    msg += getName() + " does not have a default parameter";
    std::ostringstream oss;
    oss << "Error in " << netlistLocation() << "\n" << msg;
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL, oss.str());
  }

  return setParam(defaultParamName, val);
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
  if (defaultParamName == "")
  {
    return 0.0;
  }

  double result = 0.0;

  getParam(defaultParamName, result);

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
double DeviceEntity::setDependentParameter (N_UTL_Param & par,
                                            double *res,
                                            ParameterType::ExprAccess depend)

{
  Depend dependentParam;
  setDependentParameter(par,dependentParam, depend);

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
double DeviceEntity::setDependentParameter (N_UTL_Param & par,
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

void DeviceEntity::setDependentParameter (N_UTL_Param & par,
                                          Depend & dependentParam,
                                          ParameterType::ExprAccess depend)
{
  vector<string> instances, leads, names, variables;

  dependentParam.name = par.tag();
  if (par.tag() == "TEMP" || par.tag() == "TNOM")
  {
    dependentParam.expr = new N_UTL_Expression ("(" + par.sVal() + ")+CONSTCtoK");
    dependentParam.expr->make_constant (string("CONSTCTOK"), CONSTCtoK);
  }
  else
  {
    dependentParam.expr = new N_UTL_Expression (par.eVal());
  }

  names.clear();
  leads.clear();
  instances.clear();
  variables.clear();

  dependentParam.expr->get_names(XEXP_NODE, names);
  dependentParam.expr->get_names(XEXP_LEAD, leads);
  dependentParam.expr->get_names(XEXP_INSTANCE, instances);
  dependentParam.expr->get_names(XEXP_VARIABLE, variables);

  //vector<string>::iterator s;
  vector<string>::iterator iterS;

  if (!(depend & ParameterType::SOLN_DEP))
  {
    if (names.size() > 0 || instances.size() > 0)
    {
      string msg("In device: " + getName() + ", Parameter: ");
      msg += par.tag();
      msg += " is not allowed to depend on voltage/current values";
      std::ostringstream oss;
      oss << "Error in " << netlistLocation() << "\n" << msg;
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL, oss.str());
    }
    if (depend & ParameterType::NO_DEP)
    {
      if (dependentParam.expr->get_num(XEXP_SPECIAL) > 0)
      {
        string msg("In device: " + getName() + ", Parameter: ");
        msg += par.tag();
        msg += " is not allowed to depend on time";
        std::ostringstream oss;
        oss << "Error in " << netlistLocation() << "\n" << msg;
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL, oss.str());
      }
    }
  }

  if (leads.size() > 0)
  {
    char type;
    int index;
    vector<string>::iterator n_i=leads.begin();
    vector<string>::iterator n_end=leads.end();
    for ( ; n_i != n_end; ++n_i)
    {
      index = n_i->find_last_of(":");
      if (index == string::npos )
        type = (*n_i)[0];
      else
        type = (*n_i)[index+1];
      if (type != 'B' && type != 'E' && type != 'H')
      {
        string msg("In device: " + getName() + ", Parameter: ");
        msg += par.tag();
        msg += ", Illegal use of lead current specification in expression: ";
        msg += dependentParam.expr->get_expression();
        std::ostringstream oss;
        oss << "Error in " << netlistLocation() << "\n" << msg;
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, oss.str());
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
    vector<double> zeros;
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
        string msg("In device: " + getName() + ", Global parameter: " + *iterS +" not found");
        std::ostringstream oss;
        oss << "Error in " << netlistLocation() << "\n" << msg;
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL, oss.str());
      }
      dependentParam.expr->set_var(*iterS,getSolverState().global_params[*iterS]);
      dependentParam.global_params.push_back(*iterS);
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
  vector<Depend>::iterator dpIter = dependentParams.begin();
  vector<Depend>::iterator end = dependentParams.end();
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
bool DeviceEntity::updateGlobalParameters(map<string,double> & global_map)
{
  vector<Depend>::iterator dpIter = dependentParams.begin();
  vector<Depend>::iterator end = dependentParams.end();
  double rval;
  int i, hi;
  bool changed = false;

  for ( ; dpIter != end ; ++dpIter)
  {
    if (!dpIter->global_params.empty())
    {
      vector<string>::iterator gp=dpIter->global_params.begin();
      vector<string>::iterator gend=dpIter->global_params.end();
      for ( ; gp != gend; ++gp)
      {
        if (global_map.find(*gp) == global_map.end())
        {
          string msg("Failed to find global param in map: ");
          msg += *gp;
          N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
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

  vector<Depend>::iterator dpIter = dependentParams.begin();
  vector<Depend>::iterator end = dependentParams.end();
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

  vector<Depend>::iterator dpIter = dependentParams.begin();
  vector<Depend>::iterator end = dependentParams.end();
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
bool DeviceEntity::getParamBreakpoints( vector<N_UTL_BreakPoint> & breakPointTimes )
{
  double bTime;

  vector<Depend>::iterator dpIter = dependentParams.begin();
  vector<Depend>::iterator end = dependentParams.end();
  for ( ; dpIter != end; ++dpIter)
  {
    bTime = dpIter->expr->get_break_time();
    if (bTime > getSolverState().currTime)
      breakPointTimes.push_back(bTime);
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::setDefaultParams
// Purpose       : Set parameters according to default values
// Special Notes : Normally called before setParams
// Scope         : protected
// Creator       : Dave Shirley, PSSI
// Creation Date : 12/14/04
//-----------------------------------------------------------------------------
void DeviceEntity::setDefaultParams ( )
{
// First, allocate and zero out original and given vals
  ParameterMap::const_iterator par_i = (*getPMap()).begin();
  ParameterMap::const_iterator par_end = (*getPMap()).end();
  for ( ; par_i != par_end; ++par_i)
  {
    Pars &param = static_cast<Pars &>(*(*par_i).second);
    if (param.hasGivenMember())
      param.setGiven(*this, false);

    if (param.isType<double>()) {
      if (param.getExpressionAccess() & MIN_RES)
      {
        setDefaultValue<double>(param, devOptions.minRes);
      }
      else if (param.getExpressionAccess() & MIN_CAP)
      {
        setDefaultValue<double>(param, devOptions.minCap);
      }
      param.value<double>(*this) = getDefaultValue<double>(param);
    }
    else if (param.isType<bool>())
      param.value<bool>(*this) = getDefaultValue<bool>(param);
    else if (param.isType<int>())
      param.value<int>(*this) = getDefaultValue<int>(param);
    else if (param.isType<long>())
      param.value<long>(*this) = getDefaultValue<long>(param);
    else if (param.isType<std::string>())
      param.value<std::string>(*this) = getDefaultValue<std::string>(param);
    else if (param.isType<std::vector<int> >())
      (param.value<std::vector<int> >(*this)).clear();
    else if (param.isType<std::vector<double> >())
      (param.value<std::vector<double> >(*this)).clear();
    else if (param.isType<std::vector<std::string> >())
      (param.value<std::vector<std::string> >(*this)).clear();
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::setParams
// Purpose       : Set parameters according to a vector of params.  Used to
//                 set instance or model parameter sets to netlist values
// Special Notes :
// Scope         : protected
// Creator       : Dave Shirley, PSSI
// Creation Date : 11/20/04
//-----------------------------------------------------------------------------
void DeviceEntity::setParams(vector<Param> & params)
{
  vector<string> vecCompStrings;
  map <string,vector<CompositeParam *> > compList;
  vector<CompositeParam *>::iterator compVecIter;
  bool vc_stat;
  int ipar=0;
  int iparSize=params.size();

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << endl << "In DeviceEntity::setParams, for ";
    if (dynamic_cast<DeviceModel * const>(this))
    {
      cout << "device model";
    }
    else if (dynamic_cast<DeviceInstance * const>(this))
    {
      cout << "device instance";
    }
    else
    {
      cout << "unknown entity";
    }
    cout << ": " << getName();
    if (lineNumber_ > 0)
    {
      cout << " from file: " << netlistFileName_ << " at line: " << lineNumber_;
    }
    cout << " parameters are:" << endl;
    for (ipar=0;ipar<iparSize;++ipar)
    {
      const Param & param = params[ipar];
      cout << "Param = " << param.tag() << ", Type = ";
      int tmpType = param.getType();
      switch (tmpType)
      {
        case STR:
          cout << "STR";
          break;
        case DBLE:
          cout << "DBLE";
          break;
        case INT:
          cout << "INT";
          break;
        case LNG:
          cout << "LNG";
          break;
        case EXPR:
          cout << "EXPR";
          break;
        case BOOL:
          cout << "BOOL";
          break;
        case STR_VEC:
          cout << "STR_VEC";
          break;
        case INT_VEC:
          cout << "INT_VEC";
          break;
        case DBLE_VEC:
          cout << "DBLE_VEC";
          break;
        case DBLE_VEC_IND:
          cout << "DBLE_VEC_IND";
          break;
        case COMPOSITE:
          cout << "COMPOSITE";
          break;
        default:
          cout << "Unknown";
      }
      cout << ", Value = " << param.sVal();

      if (param.given())
      {
        cout << "  given=TRUE";
      }
      else
      {
        cout << "  given=FALSE";
      }

      if (param.default_val())
      {
        cout << "  default=TRUE" << endl;
      }
      else
      {
        cout << "  default=FALSE" << endl;
      }
    }
    cout << endl;
  }
#endif

  for (ipar=0;ipar<iparSize;++ipar)
  {
    Param &param = params[ipar];
    ExtendedString tagES(param.tag());
    Param ndParam(param); // erk: not sure what ndParam is for...

    ParameterMap::const_iterator nParametricDataIter = (*getPMap()).find(tagES);
    if (nParametricDataIter != (*getPMap()).end())
    {
      const Pars &npar = static_cast<const Pars &>(*(*nParametricDataIter).second);
      if (npar.getExpressionAccess() & ParameterType::NO_INPUT)
      {
        string msg("DeviceEntity::setParams: parameter: ");
        msg += tagES;
        msg += " cannot be initialized";
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
      }
      if (npar.hasGivenMember())
      {
        if (ndParam.given())
        {
          npar.setGiven(*this, true);
        }
        else if (npar.getGiven(*this))
        {
          continue;
        }
      }
      setValueGiven(this, npar.getSerialNumber(), ndParam.given());
      if (ndParam.given() || ndParam.default_val())
      {
        if ( ndParam.getType() == EXPR )
        {
          if (npar.isType<double>())
          {
            double val = setDependentParameter (ndParam, &(npar.value<double>(*this)), npar.getExpressionAccess());
            param.setVal(val);
          }
          else if (npar.isType<std::vector<double> >())
          {
            int ind = (npar.value<std::vector<double> >(*this)).size();
            double val = setDependentParameter (ndParam, &(npar.value<std::vector<double> >(*this)), ind, npar.getExpressionAccess());
            (npar.value<std::vector<double> >(*this)).push_back(val);
          }
          else
          {
            string msg("DeviceEntity::setParams: non double param: ");
            msg += tagES;
            msg += " cannot be set to expression";
            N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
          }
        }
        else
        {
          if (npar.isType<double>())
          {
            if (!ndParam.isNumeric())
            {
              string msg("Cannot convert parameter: ");
              msg += tagES;
              msg += " to a numeric value, with input value of: ";
              msg += param.sVal();
              std::ostringstream oss;
              oss << "Error in " << netlistLocation() << "\n" << msg;
              N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, oss.str());
            }
            npar.value<double>(*this) = ndParam.dVal();
            if (tagES == "TNOM" || tagES == "TEMP")
            {
              npar.value<double>(*this) += CONSTCtoK;
            }
            if (npar.getOriginalValueIndex() >= 0)
            {
              setOriginalValue(this, npar.getOriginalValueIndex(), npar.value<double>(*this));
            }
          }
          else if (npar.isType<std::string>())
          {
            npar.value<std::string>(*this) = ndParam.sVal();
          }
          else if (npar.isType<int>())
          {
            if (!ndParam.isInteger())
            {
              string msg("Cannot convert parameter: ");
              msg += tagES;
              msg += " to an integer value, with input value of: ";
              msg += param.sVal();
              std::ostringstream oss;
              oss << "Error in " << netlistLocation() << "\n" << msg;
              N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, oss.str());
            }
            npar.value<int>(*this) = ndParam.iVal();
            if (npar.getOriginalValueIndex() >= 0)
            {
              setOriginalValue(this, npar.getOriginalValueIndex(), static_cast<double> (npar.value<int>(*this)));
            }
          }
          else if (npar.isType<long>())
          {
            if (!ndParam.isInteger())
            {
              string msg("Cannot convert parameter: ");
              msg += tagES;
              msg += " to an integer value, with input value of: ";
              msg += param.sVal();
              std::ostringstream oss;
              oss << "Error in " << netlistLocation() << "\n" << msg;
              N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, oss.str());
            }
            npar.value<long>(*this) = ndParam.lVal();
            if (npar.getOriginalValueIndex() >= 0)
            {
              setOriginalValue(this, npar.getOriginalValueIndex(), static_cast<double> (npar.value<long>(*this)));
            }
          }
          else if (npar.isType<bool>())
          {
            if (!ndParam.isBool())
            {
              string msg("Cannot convert parameter: " + tagES);
              msg += " to a logical value, with input value of: ";
              msg += param.sVal();
              std::ostringstream oss;
              oss << "Error in " << netlistLocation() << "\n" << msg;
              N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, oss.str());
            }
            npar.value<bool>(*this) = ndParam.bVal();
            if (npar.getOriginalValueIndex() >= 0)
            {
              if (npar.value<bool>(*this))
              {
                setOriginalValue(this, npar.getOriginalValueIndex(), 1.0);
              }
              else
              {
                setOriginalValue(this, npar.getOriginalValueIndex(), 0.0);
              }
            }
          }
          else if (npar.isType<std::vector<int> >())
          {
            if (ndParam.getType() == INT_VEC)
            {
              (npar.value<std::vector<int> >(*this)) = ndParam.iVecVal();
            }
            else if (ndParam.getType() == INT)
            {
              (npar.value<std::vector<int> >(*this)).push_back(ndParam.iVal());
            }
          }
          else if (npar.isType<std::vector<double> >())
          {
            if (ndParam.getType() == DBLE_VEC)
            {
              (npar.value<std::vector<double> >(*this)) = ndParam.dVecVal();
            }
            else if (ndParam.getType() == DBLE)
            {
              (npar.value<std::vector<double> >(*this)).push_back(ndParam.dVal());
            }
          }
          else if (npar.isType<std::vector<std::string> >())
          {
            if (ndParam.getType() == STR_VEC)
            {
              (npar.value<std::vector<std::string> >(*this)) = ndParam.sVecVal();
            }
            else if (ndParam.getType() == STR)
            {
              (npar.value<std::vector<std::string> >(*this)).push_back(ndParam.sVal());
            }
          }
          else if (npar.isType<CompositeMap>())
          {
            compList[tagES].clear();
            vecCompStrings.push_back(tagES);

#ifdef Xyce_DEBUG_DEVICE
            // Note: ERK.  This push-back is done to process the base-param tag of a vector composite.
            // For example, if the composite parameters are things like REGION0.XWIDTH, where REGION
            // is the base parameter tag, 0 is the index, and XWIDTH is the subcomponent, the
            // tag that should be pushed back is tagES=REGION.
            //
            // Note: ERK:  This function seems to implicitly rely on the base parameter always preceeding
            // subcomponent parameters.  So REGION (alone) should preceed REGION0.XWIDTH in  the
            // STL vector params that is passed into this function.  If it doesn't then vc_stat will
            // stay false and a fatal error will get thrown.
            if (getDeviceOptions().debugLevel > 0)
            {
              cout << "pushing back composite: tagES = " << tagES << endl;
            }
#endif
          }
          else
          {
            string msg("DeviceEntity::setParams: unknown type");
            N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
          }
        }
      }
    }
    else
    {
      if (param.sVal() == "VECTOR")
      {
      }
      else
      {
        vc_stat = false;
        std::string::size_type dot = tagES.find_first_of('.');
#ifdef Xyce_DEBUG_DEVICE
        if (getDeviceOptions().debugLevel > 0)
        {
          cout << "Inside of VECTOR-COMPOSITE: tagES = " << tagES ;
          cout << "  size of vecCompStrings = " << vecCompStrings.size();

          cout << "  Result of dot find: ";
          if (dot != string::npos)
          {
            cout << "Found it";
          }
          else
          {
            cout << "Not found";
          }
          cout << endl;
        }
#endif
        if (dot != string::npos)
        {
          int ivc=0;
          for (ivc=0;ivc<vecCompStrings.size();++ivc)
          {
            string & vcs = vecCompStrings[ivc];

            if (tagES.find(vcs) == 0)
            {
              string numString(tagES,(vcs).size(),dot-(vcs).size());
              string paramName(tagES,dot+1,tagES.size()-dot);

              vc_stat = true;

              int i=0, n=0;
              for (i=0 ; i<numString.size() ; ++i)
              {
                n *= 10;
                char c = numString[i];
                n += (static_cast<int> (c)) - (static_cast<int> ('0'));
              }
              if (paramName == "NAME")
              {
                if (n != compList[vcs].size())
                {
                  string msg("DeviceEntity::setParams: internal error filling 'NAME' vector");
                  msg += " (param: ";
                  msg += vcs;
                  msg += ")";
                  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
                }

                string tmpStr(ndParam.sVal());
                compList[vcs].push_back(constructComposite(vcs, tmpStr));
                compList[vcs][n]->setDefaultParams();
                Pars &nparComp = static_cast<Pars &>(*const_cast<ParameterMap &>(*getPMap())[vcs]);
                nparComp.value<CompositeMap>(*this)[tmpStr] = compList[vcs][n];
              }
              else
              {
                if (n >= compList[vcs].size())
                {
                  string msg("DeviceEntity::setParams: internal error in ");
                  msg += "vector-composite, 'NAME' must come first (param: ";
                  msg += vcs;
                  msg += ")";
                  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
                }
              }
              compList[vcs][n]->setParams(paramName, ndParam);
            }
          }
        }
        if (!vc_stat)
        {
          string msg("DeviceEntity::setParams: undefined parameter: ");
          msg += tagES;
          msg += "\nThis parameter is in metadata, but not recognized in constructor";
          N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
        }
      }
    }
  }
  if (!vecCompStrings.empty())
  {
    int ivc=0;
    for (ivc=0;ivc<vecCompStrings.size();++ivc)
    {
      string & vcs = vecCompStrings[ivc];

      compVecIter = compList[vcs].begin();
      for ( ; compVecIter != compList[vcs].end() ; ++compVecIter)
      {
        (*compVecIter)->processParams();
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::printFormattedOutputParam
//
// Purpose       : This function finds a parameter in the par table, and then
//                 formats its value in a string.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/4/12
//-----------------------------------------------------------------------------
std::ostream &DeviceEntity::printFormattedOutputParam(std::ostream &os, const string & paramName) const
{
  vector<int>::const_iterator v_i;
  string description;

  ExtendedString tmpName(paramName);
  tmpName.toUpper ();

  ParameterMap::const_iterator paramIter = (*getPMap()).find(tmpName);

  ExtendedString myName(getName());
  ExtendedString type("Unknown");
  string label;

  if (paramIter != (*getPMap()).end())
  {
    string fstring((*paramIter).first);

    const Pars &param = static_cast<const Pars &>(*(*paramIter).second);

    //if ( !(param.getExpressionAccess() & Pars::NO_INPUT) )
    {
      if (param.isType<double>())
      {
        if (fstring == "TNOM" || fstring == "TEMP")
        {
          os << param.value<double>(*this)-CONSTCtoK;
        }
        else
        {
          os << param.value<double>(*this);
        }
      }
      else if (param.isType<bool>())
      {
        os << param.value<bool>(*this);
      }
      else if (param.isType<int>())
      {
        os << param.value<int>(*this);
      }
      else if (param.isType<long>())
      {
        os << param.value<long>(*this);
      }
      else if (param.isType<std::string>())
      {
        os << param.value<std::string>(*this);
      }
      else if (param.isType<std::vector<std::string> >())
      {
        os << " (string vector) : ";
        os << "length = " << (param.value<std::vector<std::string> >(*this)).size();
        if ((param.value<std::vector<std::string> >(*this)).size() > 0)
        {
          os << " :";
          vector<string>::const_iterator iterStringVec;
          iterStringVec = (param.value<std::vector<std::string> >(*this)).begin();
          for ( ; iterStringVec != (param.value<std::vector<std::string> >(*this)).end() ; ++iterStringVec)
          {
            os << "  " << *iterStringVec;
          }
        }
      }
      else if (param.isType<std::vector<int> >())
      {
        os << " (int vector) : ";
        os << "length = " << (param.value<std::vector<int> >(*this)).size();
        if ((param.value<std::vector<int> >(*this)).size() > 0)
        {
          os << " :";
          v_i = (param.value<std::vector<int> >(*this)).begin();
          for ( ; v_i != (param.value<std::vector<int> >(*this)).end() ; ++v_i)
          {
            os << "  " << *v_i;
          }
        }
      }
      else if (param.isType<std::vector<double> >())
      {
        os << " (double vector) : ";
        os << "length = " << (param.value<std::vector<double> >(*this)).size();
        if ((param.value<std::vector<double> >(*this)).size() > 0)
        {
          os << " :";
          vector<double>::const_iterator iterDoubleVec;
          iterDoubleVec = (param.value<std::vector<double> >(*this)).begin();
          for ( ; iterDoubleVec != (param.value<std::vector<double> >(*this)).end() ; ++iterDoubleVec)
          {
            os << "  " << *iterDoubleVec;
          }
        }
      }
      else if (param.isType<CompositeMap>())
      {
        os << " (composite) : " << (param.value<CompositeMap>(*this)).size();
        if ((param.value<CompositeMap>(*this)).size() > 0)
        {
          map <string,CompositeParam *>::const_iterator c_i;
          c_i = (param.value<CompositeMap>(*this)).begin();
          (*c_i).second->printParams(os);
        }
      }
    }
  }
  return os;
}

//-----------------------------------------------------------------------------
// Function      : DeviceEntity::given
// Purpose       : Return whether param was given
// Special Notes :
// Scope         : protected
// Creator       : Dave Shirley, PSSI
// Creation Date : 11/23/04
//-----------------------------------------------------------------------------
bool DeviceEntity::given( const string & parameter_name )
{
  ParameterMap::const_iterator it = (*getPMap()).find(parameter_name);
  if (it == (*getPMap()).end())
  {
    string msg("DeviceEntity::Given: unrecognized param: ");
    msg += parameter_name;
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
  }

  return wasValueGiven(this, (*it).second->getSerialNumber());
}

void populateParams(const DeviceEntity::ParameterMap &parameter_map, vector<Param> & param_list, DeviceParamMap &composite_param_list)
{
  for (DeviceEntity::ParameterMap::const_iterator it = parameter_map.begin(); it != parameter_map.end(); ++it)
  {
    const DeviceEntity::Pars &param = static_cast<const DeviceEntity::Pars &>(*(*it).second);

    if (!(param.getExpressionAccess() & ParameterType::NO_INPUT))
    {
      if (param.isType<double>())
      {
        if (param.getVec() == 0)
        {
          double val;
          if ((*it).first == "TNOM" || (*it).first == "TEMP")
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
            string vPar((*it).first.substr(0, (*it).first.size()-1));
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
        else if (param.getVec() == 1)
        {
          string vPar((*it).first.substr(0, (*it).first.size()-1));
          param_list.push_back(Param(vPar, "VECTOR"));
        }
      }
      else if (param.isType<int>())
      {
        if (param.getVec() == 0)
          param_list.push_back(Param((*it).first, getDefaultValue<int>(param)));
        else if (param.getVec() == 1)
        {
          string vPar((*it).first.substr(0, (*it).first.size()-1));
          param_list.push_back(Param(vPar, "VECTOR"));
        }
      }
      else if (param.isType<std::string>())
      {
        if (param.getVec() == 0)
          param_list.push_back(Param((*it).first, getDefaultValue<std::string>(param)));
        else if (param.getVec() == 1)
        {
          string vPar((*it).first.substr(0, (*it).first.size()-1));
          param_list.push_back(Param(vPar, "VECTOR"));
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
      else if (param.isType<CompositeMap>())
      {
        Param vc2((*it).first, "VECTOR-COMPOSITE");
        vc2.setDefault(true);
        param_list.push_back(vc2);

        vector<Param> compositeParams;
        const ParametricData<CompositeParam> *c = param.getCompositeParametricData<CompositeParam>();

        if (c == 0)
        {
          string msg("Error: vector-composite map for device type entity empty.  You need to modify the device factory.");
          N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg );
        }

        // TODO: [DGB] I think when the Pars is refactored this will be clearer.  But this basically adds the
        //   type to the composite list with 'NAME' first.
        const ParametricData<CompositeParam> &d = *c;
        const ParametricData<CompositeParam>::ParameterMap &e = d.getMap();

        for (ParametricData<CompositeParam>::ParameterMap::const_iterator it4 = e.find("NAME"); it4 != e.end();) {
          const Descriptor &p = static_cast<const Descriptor &>(*(*it4).second);
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

        composite_param_list[(*it).first] = compositeParams;
      }
      else
      {
//        Just skip these, like list of coupled inductors because not needed for metadata
//        string msg("DeviceEntity::getParams: Type not supported");
//        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
        std::cout << "In final else clause of DeviceEntity::getParams().";
        if( param.isType<std::vector<std::string> >() )
          std::cout << " type is STR_VEC ";
        if( param.isType<std::vector<double> >() )
          std::cout << " type is DBLE_VEC ";
        std::cout << it->first << " this item is NOT being added to default parameter list." << std::endl;

      }
    }
  }
}

netlistLocation_ DeviceEntity::netlistLocation() const {
  return netlistLocation_(netlistFileName_, lineNumber_);
}

} // namespace Device
} // namespace Xyce

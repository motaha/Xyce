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
//  (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_DeviceEntity.h,v $
//
// Purpose        : This file contains the device entity base class.
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
// Revision Number: $Revision: 1.117.2.1 $
//
// Revision Date  : $Date: 2014/02/26 20:16:30 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------


#ifndef Xyce_N_DEV_DeviceEntity_h
#define Xyce_N_DEV_DeviceEntity_h

#include <iosfwd>
#include <map>
#include <string>
#include <vector>

#include <N_DEV_fwd.h>
#include <N_IO_fwd.h>
#include <N_UTL_fwd.h>
#include <N_UTL_NetlistLocation.h>
#include <N_DEV_Pars.h>

class N_LAS_Vector;

namespace Xyce {
namespace Device {

typedef std::map<std::string, std::vector<Param>, LessNoCase> CompositeParamMap;

void populateParams(const ParameterMap &parameter_map, std::vector<Param> &param_list, CompositeParamMap &composite_param_map);
void setParameters(CompositeParam &composite_param, const std::string &pName, const Param &ndParam );
void setParameters(DeviceEntity &entity, std::vector<Param>::const_iterator begin, std::vector<Param>::const_iterator end, const DeviceOptions &device_options);

std::ostream &printParameter(std::ostream &os, const DeviceEntity &entity, const std::string &name, const Descriptor &param);
std::ostream &printCompositeParameters(std::ostream &os, const CompositeParam &composite);

struct Depend
{
  std::string                 name;
  Util::Expression *          expr;
  union resUnion
  {
    double *                result;
    std::vector<double> *   resVec;
  } resultU;
  int                         vectorIndex;
  std::vector<double>         vals;
  std::vector<std::string>    global_params;
  int                         n_vars, lo_var;
};

//-----------------------------------------------------------------------------
// Class         : DeviceEntity
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/11/02
//-----------------------------------------------------------------------------
class DeviceEntity : public ParameterBase
{
public:
  DeviceEntity(
     const char * const        entity_type,
     const std::string &       device_name,
     ParametricData<void> &    parametric_data,
     const SolverState &       solver_state,
     const DeviceOptions &     device_options,
     const std::string &       netlist_path,
     int                       netlist_line);

private:
  DeviceEntity(const DeviceEntity &);
  DeviceEntity &operator=(const DeviceEntity &);

public:
  virtual ~DeviceEntity();

  virtual bool processParams() = 0;

  virtual bool processInstanceParams() = 0;

  virtual CompositeParam *constructComposite(const std::string &composite_name, const std::string &param_name) 
  {
    return NULL;
  }

  bool setDefaultParam(double val);
  double getDefaultParam();

  bool scaleParam(const std::string & paramName, double val, double val0);
  bool scaleParam(const std::string & paramName, double val);
  bool scaleDefaultParam(double val);

  bool setParam(const std::string & paramName, double val);
  bool getParam(const std::string & paramName, double & result);
  bool getParamBreakpoints( std::vector<Util::BreakPoint> & );

  bool updateDependentParameters(N_LAS_Vector & vars);
  bool updateDependentParameters(double temp_tmp);
  bool updateGlobalParameters(std::map<std::string, double> &);
  bool updateDependentParameters();

  double setDependentParameter(Util::Param &, double *, ParameterType::ExprAccess);
  double setDependentParameter(Util::Param &, std::vector<double> *, int , ParameterType::ExprAccess);
  void setDependentParameter(Util::Param & par, Depend & dependentParam, ParameterType::ExprAccess depend);

  void setDefaultParams() 
  {
    setDefaultParameters(*this, getParameterMap().begin(), getParameterMap().end(), devOptions_);
  }

  void setParams(const std::vector<Param> & params) 
  {
    setParameters(*this, params.begin(), params.end(), devOptions_);
  }

public:
  bool given(const std::string & parameter_name) const;

  const char *getEntityType() const 
  {
    return entityType_;
  }

  const std::string &getName() const 
  {
    return name_;
  }

  void setDefaultParamName(const std::string &default_param_name) 
  {
    defaultParamName_ = default_param_name;
  }

  const std::vector<Depend> &getDependentParams() 
  {
    return dependentParams;
  }

  const DeviceOptions &getDeviceOptions() const 
  {
    return devOptions_;
  }

  const SolverState &getSolverState() const 
  {
    return solState_;
  }

  const NetlistLocation &netlistLocation() const 
  {
    return netlistLocation_;
  }

  const ParameterMap &getParameterMap() const 
  {
    return parametricData_.getMap();
  }

private:
  void escape(std::string &) const;
  void checkDepend(ParameterType::ExprAccess &);

private:
  const char * const          entityType_;
  std::string                 name_;
  std::string                 defaultParamName_;
  ParametricData<void> &      parametricData_;
  NetlistLocation             netlistLocation_;

  const SolverState &         solState_;
  const DeviceOptions &       devOptions_;

protected:
  std::vector<Depend>         dependentParams;
  std::vector<int>            expVarGIDs;
  std::vector<int>            expVarLIDs;
  std::vector<std::string>    expVarNames;
  std::vector<double>         expVarVals;
  std::vector<double>         eVarVals;
};

struct DeviceEntityCmp : public std::binary_function<DeviceEntity, DeviceEntity, bool>

{
  bool operator()(const DeviceEntity &entity_0, const DeviceEntity &entity_1) const 
  {
    return less_nocase(entity_0.getName(), entity_1.getName());
  }
  bool operator()(const DeviceEntity *entity_0, const DeviceEntity *entity_1) const 
  {
    return less_nocase(entity_0->getName(), entity_1->getName());
  }
};

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::Depend sDepend;

#endif // Xyce_N_DEV_DeviceEntity_h

//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright(C) 2002-2013  Sandia Corporation
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
// Revision Number: $Revision: 1.102.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:37 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------


#ifndef Xyce_N_DEV_DeviceEntity_h
#define Xyce_N_DEV_DeviceEntity_h

// ---------- Standard Includes ----------
#include <iosfwd>
#include <list>
#include <map>
#include <string>
#include <vector>

// ------------- Xyce Includes ------------

#include <N_DEV_fwd.h>
#include <N_DEV_CompositeParam.h>
#include <N_DEV_Units.h>
#include <N_DEV_Pars.h>
#include <N_ERH_ErrorMgr.h>
#include <N_LAS_Vector.h>
#include <N_UTL_BreakPoint.h>
#include <N_UTL_Xyce.h>
#include <N_UTL_NoCase.h>

// ---------- Forward Declarations ----------
class N_IO_CmdParse;
class N_UTL_Expression;

namespace Xyce {
namespace Device {

struct netlistLocation_;

void populateParams(const ParametricData<DeviceEntity>::ParameterMap &parameter_map, std::vector<Param> & param_list, DeviceParamMap &composite_param_list);

struct Depend
{
    std::string                   name;
    N_UTL_Expression *            expr;
    union resUnion
    {
        double *                    result;
        std::vector<double> *       resVec;
    } resultU;
    int                           vectorIndex;
    std::vector<double>           vals;
    std::vector<std::string>      global_params;
    int                           n_vars, lo_var;
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
    typedef Descriptor Pars;
    typedef ParametricData<void>::ParameterMap ParameterMap;

    DeviceEntity(
      SolverState &       solver_state,
      DeviceOptions &     device_options,
      const std::string & device_name,
      const std::string & net_list_file_path,
      int                 net_list_file_line_number);

  private:
    DeviceEntity(const DeviceEntity &);
    DeviceEntity &operator=(const DeviceEntity &);

  public:
    virtual ~DeviceEntity();

    virtual bool processParams(std::string param = "") = 0;
    virtual bool processInstanceParams(std::string param = "");
    virtual CompositeParam *constructComposite(std::string & compositeName, std::string & paramName) {return NULL;}

    bool setDefaultParam(double val);
    double getDefaultParam();

    bool scaleParam(const std::string & paramName, double val, double val0);
    bool scaleParam(const std::string & paramName, double val);
    bool scaleDefaultParam(double val);

    bool setParam(const std::string & paramName, double val);
    bool getParam(const std::string & paramName, double & result);
    bool getParamBreakpoints( std::vector<N_UTL_BreakPoint> & );

    bool updateDependentParameters(N_LAS_Vector & vars);
    bool updateDependentParameters(double temp_tmp);
    bool updateGlobalParameters(map<std::string,double> &);
    bool updateDependentParameters();

#if 0
    bool perturbSensParam(Param & ndParam);
#endif

    std::ostream &printFormattedOutputParam(std::ostream &os, const std::string & paramName) const;

  protected:
    double setDependentParameter(N_UTL_Param &, double *, ParameterType::ExprAccess);
    double setDependentParameter(N_UTL_Param &, std::vector<double> *, int , ParameterType::ExprAccess);
    void setDependentParameter(N_UTL_Param & par, Depend & dependentParam, ParameterType::ExprAccess depend);

    void setDefaultParams();
    void setParams(std::vector<Param> & params);

  public:
    bool given(const std::string & parameter_name);

  public:
    const std::string &getName() const {
      return name_;
    }

    void setName(const std::string &name) {
      name_ = name;
    }

    const std::vector<Depend> &getDependentParams() {
      return dependentParams;
    }

    const DeviceOptions &getDeviceOptions() const {
      return devOptions;
    }

    DeviceOptions &getDeviceOptions() {
      return devOptions;
    }

    const SolverState &getSolverState() const {
      return solState;
    }

    SolverState &getSolverState() {
      return solState;
    }

    netlistLocation_ netlistLocation() const;

    // [DGB]  Hopefully htis crap goes away soon!
    virtual const ParametricData<void> &getMyParametricData() const = 0;

    const ParameterMap *getPMap() const {
      return reinterpret_cast<const ParameterMap *>(&getMyParametricData().getMap());
    }

  private:
    void escape(std::string &) const;
    void checkDepend(ParameterType::ExprAccess &);

  private:
    std::string                   name_;

    std::string                   netlistFileName_;
    int                           lineNumber_;

  protected:
    std::string                   defaultParamName;
    double                        paramNew;
    double                        paramOld;
    SolverState &                 solState;
    DeviceOptions &               devOptions;
    N_IO_CmdParse &               commandLine;

    std::vector<Depend>           dependentParams;
    std::vector<int>              expVarGIDs;
    std::vector<int>              expVarLIDs;
    std::vector<std::string>      expVarNames;
    std::vector<double>           expVarVals;
    std::vector<double>           eVarVals;
};

struct DeviceEntityCmp : public std::binary_function<DeviceEntity, DeviceEntity, bool>
{
    bool operator()(const DeviceEntity &entity_0, const DeviceEntity &entity_1) const {
      return less_nocase(entity_0.getName(), entity_1.getName());
    }
    bool operator()(const DeviceEntity *entity_0, const DeviceEntity *entity_1) const {
      return less_nocase(entity_0->getName(), entity_1->getName());
    }
};


struct netlistLocation_
{
  netlistLocation_(const std::string &file_path, int line_number)
    : filePath(file_path),
      lineNumber(line_number)
  {}

  const std::string &filePath;
  const int lineNumber;
};

inline std::ostream &operator<<(std::ostream &os, const netlistLocation_ &x) {
  os << "file " << x.filePath << " at or near line " << x.lineNumber;
  return os;
}

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::DeviceEntity N_DEV_DeviceEntity;
typedef Xyce::Device::Depend sDepend;

std::ostream &dumpConfiguration(std::ostream &os, const Configuration &configuration);
std::ostream &dumpParams(std::ostream &os, const N_DEV_DeviceEntity &entity);

#endif

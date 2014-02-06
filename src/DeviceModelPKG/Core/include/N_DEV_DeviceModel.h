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

//-----------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_DeviceModel.h,v $
//
// Purpose        : This file contains the device model base class.
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
// Revision Number: $Revision: 1.38.2.4 $
//
// Revision Date  : $Date: 2013/10/03 17:23:37 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_DeviceModel_h
#define Xyce_N_DEV_DeviceModel_h

#include <map>
#include <string>
#include <vector>

#include <N_DEV_fwd.h>
#include <N_DEV_DeviceEntity.h>
#include <N_DEV_Pars.h>

namespace Xyce {
namespace Device {

/** 
 * @class DeviceModel N_DEV_DeviceModel.h
 *
 * @author Eric Keiter, SNL, Parallel Computational Sciences
 * @date   4/03/00
 */
class DeviceModel : public DeviceEntity
{
  enum mType {TEMP, DOSE};
  enum iType {LIN, QUAD, PWL};
  enum fitType {LINEAR_FIT, LOG_FIT};

public:
    /** 
     * Add the parameter "TEMPMODEL" to the parametric_data.
     *
     * @param parametric_data 
     */
    template<class T>
    static void initThermalModel(ParametricData<T> &parametric_data) {
      parametric_data.addPar("TEMPMODEL", "NONE", false, ParameterType::NO_DEP, &DeviceModel::temperatureModel, NULL, U_NONE, CAT_CONTROL,
                             "Specification to type of parameter interpolation over temperature (see User's Guide section 5.3)");
    }

    /** 
     * Add the parameter "DOSEMODEL" to the parametric_data.
     *
     * @param parametric_data 
     */
    template<class T>
    static void initDoseModel(ParametricData<T> &parametric_data) {
      parametric_data.addPar ("DOSEMODEL", "NONE", false, ParameterType::NO_DEP, &DeviceModel::doseModel, NULL);
    }
    
  DeviceModel(const ModelBlock &model_block, SolverState &solver_state, DeviceOptions &device_options);

  virtual ~DeviceModel ();

private:
  DeviceModel();
  DeviceModel(const DeviceModel &);
  DeviceModel &operator=(const DeviceModel &);

public:
  void setModParams(std::vector<Param> params);

  virtual std::ostream &printOutInstances(std::ostream &os) const = 0;

    /** 
     * processParams 
     *
     * @param param 
     *
     * @return true if parameter processing was successful
     */
    virtual bool processParams(std::string param = "") = 0;

    /** 
     * processInstanceParams 
     *
     * @param param 
     *
     * @return true if parameter processing was successful
     */
    virtual bool processInstanceParams(std::string param = "") = 0;

  virtual bool clearTemperatureData () {return true;}

#if 0
  void perturbSensParam (Param & ndParam);
#endif

  void saveParams ();
  bool interpolateTNOM (double);
  bool interpolateDOSE (double);
  void restoreParams ();

  virtual bool getBinPrefixFlag () {return false;}

private:
  bool interpolated ();
  bool interpolate (double);

public:
  int getLevel() const {
    return level_;
  }

  void setLevel(int level) {
    level_ = level;
  }

  const std::string &getType() const {
    return type_;
  }

private:
  std::string                           type_;
  int                                   level_;
  std::string                           temperatureModel;
  std::string                           doseModel;
  mType                                 iModel;
  iType                                 iMethod;
  double                                base_temp;
  std::map<std::string, int>            fitMap;
  std::vector<double DeviceEntity::*>   fitParams;
  std::vector<double>                   oldParams;
  std::vector<double>                   base;
  std::vector< std::vector<double> >    fit;
  std::vector<double>                   min_par;
  std::vector<double>                   max_par;
  std::vector<fitType>                  parType;
};

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::DeviceModel N_DEV_DeviceModel;

#endif


// ----------------------------------------------------------------------------
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
// ----------------------------------------------------------------------------

//----------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_Resistor.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 02/28/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.195.2.5 $
//
// Revision Date  : $Date: 2013/10/03 17:23:38 $
//
// Current Owner  : $Author: tvrusso $
//----------------------------------------------------------------------------
#include <Xyce_config.h>

#include <N_DEV_Const.h>
#include <N_DEV_Factory.h>
#include <N_DEV_Resistor.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_MatrixLoadData.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>

namespace Xyce {
namespace Device {

/**
 * Configuration and parameter definitions for resistor instance
 *
 * @date   Tue Aug  6 08:42:40 2013
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
 */
template<>
ParametricData<Resistor::Instance>::ParametricData()
{
  // Device Configuration information
  setNumNodes(2);                                     ///< Device has two nodes
  setNumOptionalNodes(0);                             ///<   No optional nodes
  setNumFillNodes(0);                                 ///<   No fill nodes
  setModelRequired(0);                                ///< A model line is not required in the netlist
  setPrimaryParameter("R");                           ///< The first missing parameter name will be assumed to be a "R"
  addModelType("R");                                  ///< The model type is "R"

  // Create parameter definitions for parameter member variables
  addPar("R", 1000.0, false, ParameterType::TIME_DEP, &Resistor::Instance::R, U_OHM, CAT_NONE, "Resistance");
  addPar("L", 0.0, false, ParameterType::NO_DEP, &Resistor::Instance::length, U_METER, CAT_NONE, "Length");
  addPar("W", 0.0, false, ParameterType::NO_DEP, &Resistor::Instance::width, U_METER, CAT_NONE, "Width");
  addPar("TEMP", 0.0, false, ParameterType::TIME_DEP, &Resistor::Instance::temp, U_DEGC, CAT_NONE, "Temperature");

  addPar("TC1",   0.0, false,   ParameterType::NO_DEP, &Resistor::Instance::tempCoeff1, &Resistor::Instance::tempCoeff1Given, U_DEGCM1, CAT_NONE, "Linear Temperature Coefficient");
  addPar("TC2",   0.0, false,   ParameterType::NO_DEP, &Resistor::Instance::tempCoeff2, &Resistor::Instance::tempCoeff2Given, U_DEGCM2, CAT_NONE, "Quadratic Temperature Coefficient");
  makeVector("TC", 2);                                ///< Allow TC to be entered as a vector (TC=1,2)

  addPar("DTEMP",   0.0, false,   ParameterType::NO_DEP, &Resistor::Instance::dtemp, &Resistor::Instance::dtempGiven, U_DEGC, CAT_NONE, "Device Temperature -- For compatibility only. Parameter is NOT used");
}

/**
 * ParametricData<Resistor::Model>::ParametricData
 *
 *
 * @date   Tue Aug  6 08:42:32 2013
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
 *
 *
 * @return
 */
template<>
ParametricData<Resistor::Model>::ParametricData()
{
  // Create parameter definitions for parameter member variables
  addPar("TC1",   0.0, false,   ParameterType::NO_DEP, &Resistor::Model::tempCoeff1, U_DEGCM1, CAT_NONE, "Linear Temperature Coefficient");
  addPar("TC2",   0.0, false,   ParameterType::NO_DEP, &Resistor::Model::tempCoeff2, U_DEGCM2, CAT_NONE, "Quadratic Temperature Coefficient");
  addPar("RSH",   0.0, false,   ParameterType::NO_DEP, &Resistor::Model::sheetRes, U_OHM,  CAT_NONE, "Sheet Resistance");
  addPar("DEFW",  1.e-5, false, ParameterType::NO_DEP, &Resistor::Model::defWidth, U_METER,  CAT_NONE, "Default Instance Width");
  addPar("NARROW",0.0, false,   ParameterType::NO_DEP, &Resistor::Model::narrow, U_METER,  CAT_NONE, "Narrowing due to side etching");
  addPar("TNOM",  0.0, false,   ParameterType::NO_DEP, &Resistor::Model::tnom, U_DEGC, CAT_NONE, "Parameter Measurement Temperature");
}

namespace Resistor {

vector<vector<int> >
Instance::jacStamp;

/**
 * Configuration and parameter definitions of the resistor instance
 *
 * @return reference to the configuration and parameter definitions of the resistor instance singleton
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
 * @date   Mon Aug 12 09:15:58 2013
 */
ParametricData<Instance> &Instance::getParametricData() {
  static ParametricData<Instance> parMap;

  return parMap;
}

/**
 * Configuration and parameter definitions of the resistor model
 *
 * @return reference to the configuration and parameter definitions of the resistor model singleton
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
 * @date   Mon Aug 12 09:15:58 2013
 */
ParametricData<Model> &Model::getParametricData() {
  static ParametricData<Model> parMap;

  return parMap;
}

/**
 * Construct a resistor instance.
 *
 *
 * @note The parameter member variable initial values from the initializer if immediately replaced with the values from
 * the parameter definitions by the setDefaultParams() function.
 *
 * @note The matrix_load_data, solver_state and extern_data parameter are not used directly by the resistor instance,
 * but are passed on to the DeviceInstance.
 *
 * @param instance_block Instance information from parser
 * @param model Resistor model to add this instance to
 * @param matrix_load_data Solution matrix load data
 * @param solver_state Solution sover state
 * @param extern_data Solution external data
 * @param device_options Device options defined in netlist
 */
Instance::Instance(
  InstanceBlock &       instance_block,
  Model &               model,
  MatrixLoadData &      matrix_load_data,
  SolverState &         solver_state,
  ExternData &          extern_data,
  DeviceOptions &       device_options)
  : DeviceInstance(instance_block, matrix_load_data, solver_state, extern_data, device_options),
    model_(model),
    R(0.0),
    length(0.0),
    width(0.0),
    temp(device_options.temp.dVal()),
    tempCoeff1(0.0),
    tempCoeff2(0.0),
    dtemp(0.0),
    tempCoeff1Given(false),
    tempCoeff2Given(false),
    dtempGiven(false),
    G(0.0),
    i0(0.0),
    li_Pos(-1),
    li_Neg(-1),
    li_store_dev_i(0),
    APosEquPosNodeOffset(-1),
    APosEquNegNodeOffset(-1),
    ANegEquPosNodeOffset(-1),
    ANegEquNegNodeOffset(-1)
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  ,
    f_PosEquPosNodePtr(0),
    f_PosEquNegNodePtr(0),
    f_NegEquPosNodePtr(0),
    f_NegEquNegNodePtr(0)
#endif
{
  // Initialize DeviceInstance values
  numIntVars   = 0;                                   ///< Initialize number if internal nodes in DeviceInstance
  numExtVars   = 2;                                   ///< Initialize number if external nodes in DeviceInstance
  numStateVars = 0;                                   ///< Initialize number if state variables in DeviceInstance
  setNumStoreVars(0);                                 ///< Initialize number if store variables in DeviceInstance
  numLeadCurrentStoreVars = 1;                        ///< Initialize number if lead current variables in DeviceInstance

  defaultParamName = "R";                             ///< Default parameter is "R", resistance

  setName(instance_block.getName());                  ///< This instance name
  setModelName(model_.getName());                     ///< This model name

  // Initialize global resistor Jacobina stamp is it hasn't been
  if (jacStamp.empty())
  {
    jacStamp.resize(2);
    jacStamp[0].resize(2);
    jacStamp[1].resize(2);
    jacStamp[0][0] = 0;
    jacStamp[0][1] = 1;
    jacStamp[1][0] = 0;
    jacStamp[1][1] = 1;
  }

  // Set params to constant default values from parameter definition
  setDefaultParams();

  // Set params according to instance line and constant defaults from metadata
  setParams(instance_block.params);

  // Set any non-constant parameter defaults:
  if (!given("TEMP"))
    temp = device_options.temp.dVal();
  if (!given("W"))
    width = model_.defWidth;

  // Get temperature values from model is not given in instance
  if (!tempCoeff1Given)
    tempCoeff1 = model_.tempCoeff1;
  if (!tempCoeff2Given)
    tempCoeff2 = model_.tempCoeff2;

  // Calculate any parameters specified as expressions
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors
  if (!given("R"))
  {
    if (model_.given("RSH") && given("L") && (model_.sheetRes != 0) &&
        (length != 0))
    {
      R = model_.sheetRes * (length - model_.narrow)
          / (width - model_.narrow);
    }
    else
    {
      R = 1000.0;
      string msg="***********\n";
      msg += "Resistor: WARNING!  Resistance=0, ";
      msg += "set to default of 1000 ohms " + getName() + "\n";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_WARNING_0,msg);
    }
  }

  // Process the parameters to complete initialization
  processParams();
}

/**
 * Process parameter.
 *
 * @param param Parameter to be processed.
 *
 * @return true on success
 *
 * @author Eric Keiter, SNL, Parallel Computational Sciences
 * @date   6/03/02
 */
bool Instance::processParams(string param)
{
  // now set the temperature related stuff.
  return updateTemperature(temp);
}

/**
 * Register local IDs
 *
 * Register the local internal and external node IDs.
 *
 * @param intLIDVecRef internal local IDs
 * @param extLIDVecRef external local IDs
 *
 * @author Robert Hoekstra, SNL, Parallel Computational Sciences
 * @date   6/12/02
 */
void Instance::registerLIDs(
  const vector<int> & intLIDVecRef,
  const vector<int> & extLIDVecRef)
{
  string msg;

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline =
    "-----------------------------------------------------------------------------";
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << endl << dashedline << endl;
    cout << "  ResistorInstance::registerLIDs" << endl;
    cout << "  name = " << getName() << endl;
  }
#endif

  // Check if the size of the ID lists corresponds to the
  // proper number of internal and external variables.
  int numInt = intLIDVecRef.size();
  int numExt = extLIDVecRef.size();

  if (numInt != numIntVars)
  {
    msg = "ResistorInstance::registerLIDs:";
    msg += "numInt != numIntVars";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL,msg);
  }

  if (numExt != numExtVars)
  {
    msg = "ResistorInstance::registerLIDs:";
    msg += "numExt != numExtVars";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL,msg);
  }

  // Copy the local ID lists.
  intLIDVec = intLIDVecRef;                           ///< Set the internal local IDs in DeviceInstance
  extLIDVec = extLIDVecRef;                           ///< Set the external local IDs in DeviceInstance

  li_Pos = extLIDVec[0];
  li_Neg = extLIDVec[1];

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "  li_Pos_ = " << li_Pos << endl;
    cout << "  li_Neg_ = " << li_Neg << endl;
  }
#endif

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << dashedline << endl;
  }
#endif
}

/**
 * Register the local state IDs
 *
 * @note The resistor does not have any state vars.
 *
 * @param staLIDVecRef State variable local IDs
 *
 * @author Robert Hoekstra, SNL, Parallel Computational Sciences
 * @date   06/12/02
 */
void Instance::registerStateLIDs(const vector<int> & staLIDVecRef)
{
  string msg;

  // Check if the size of the ID lists corresponds to the
  // proper number of internal and external variables.
  int numSta = staLIDVecRef.size();

  if (numSta != numStateVars)
  {
    msg = "ResistorInstance::registerStateLIDs:";
    msg += "numSta != numStateVars";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL,msg);
  }
}


/**
 * Register the local store IDs
 *
 * One store var for device current.
 *
 * @param stoLIDVecRef Store variable local IDs
 *
 * @author Richard Schiek, Electrical Systems Modeling
 * @date   12/18/2012
 */
void Instance::registerStoreLIDs(const vector<int> & stoLIDVecRef)
{
  string msg;

  // Check if the size of the ID lists corresponds to the
  // proper number of internal and external variables.
  int numSto = stoLIDVecRef.size();

  if (numSto != getNumStoreVars())
  {
    msg = "ResistorInstance::registerStoreLIDs:";
    msg += "numSto != numStoreVars";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL,msg);
  }

  if (loadLeadCurrent)
  {
    li_store_dev_i = stoLIDVecRef[0];
  }
}

/**
 * Populates and returns the store name map.
 *
 * If the DeviceInstance::storeNameMap is empty, populate it first.
 *
 * @return reference to the DeviceInstance::storeNameMap
 *
 * @author Richard Schiek, Electrical Systems Modeling
 * @date   12/18/2012
 */
map<int,string> & Instance::getStoreNameMap()
{
  // set up the internal name map, if it hasn't been already.
  if (loadLeadCurrent && storeNameMap.empty ())
  {
    // change subcircuitname:devicetype_deviceName to
    // devicetype:subcircuitName:deviceName
    string modName(getName());
    spiceInternalName(modName);
    string tmpstr;
    tmpstr = modName+":DEV_I";
    storeNameMap[ li_store_dev_i ] = tmpstr;
  }

  return storeNameMap;
}

/**
 * Register the Jacobian local IDs
 *
 * @param jacLIDVec Jacobian local Ids
 *
 * @author Robert Hoekstra, SNL, Parallel Computational Sciences
 * @date   08/27/01
 */
void Instance::registerJacLIDs(const vector< vector<int> > & jacLIDVec)
{
  // Let DeviceInstance do its work.
  DeviceInstance::registerJacLIDs(jacLIDVec);

  // Store offsets of the components of the Jacobian of this instance
  APosEquPosNodeOffset = jacLIDVec[0][0];
  APosEquNegNodeOffset = jacLIDVec[0][1];
  ANegEquPosNodeOffset = jacLIDVec[1][0];
  ANegEquNegNodeOffset = jacLIDVec[1][1];
}

/**
 * Setup direct access pointer to solution matrix and vectors.
 *
 * The pointers to the martix are saved if Xyce_NONPOINTER_MATRIX_LOAD is defined in the preprocessor.
 *
 * @author Eric Keiter, SNL
 * @date   11/30/08
 */
void Instance::setupPointers()
{
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  N_LAS_Matrix & dFdx = *(extData.dFdxMatrixPtr);
  f_PosEquPosNodePtr = &(dFdx[li_Pos][APosEquPosNodeOffset]);
  f_PosEquNegNodePtr = &(dFdx[li_Pos][APosEquNegNodeOffset]);
  f_NegEquPosNodePtr = &(dFdx[li_Neg][ANegEquPosNodeOffset]);
  f_NegEquNegNodePtr = &(dFdx[li_Neg][ANegEquNegNodeOffset]);
#endif
}

/**
 * Update the intermediate variables
 *
 * @return true on success
 *
 * @author Eric R. Keiter, Dept. 9233.
 * @date   3/05/04
 */
bool Instance::updateIntermediateVars()
{
  bool bsuccess = true;

  double v_pos, v_neg;
  double * solVec = extData.nextSolVectorRawPtr;

  v_pos = solVec[li_Pos];
  v_neg = solVec[li_Neg];

  // Load RHS vector element for the positive circuit node KCL equ.
  i0 = (v_pos-v_neg)*G;

  return bsuccess;
}

/**
 * Update the primary state variables.
 *
 * @return true on success
 *
 * @author Eric Keiter, SNL, Parallel Computational Sciences
 * @date   01/29/01
 */
bool Instance::updatePrimaryState()
{
  return updateIntermediateVars();
}

/**
 * Load the DAE force vector
 *
 * Loads the F-vector contributions for a single resistor instance.
 *
 * @return true on success
 *
 * @author Eric Keiter, SNL, Parallel Computational Sciences
 * @date   01/24/03
 */
bool Instance::loadDAEFVector()
{
  double * fVec = extData.daeFVectorRawPtr;
  fVec[li_Pos] += i0;
  fVec[li_Neg] += -i0;

  if( loadLeadCurrent )
  {
    double * stoVec = extData.nextStoVectorRawPtr;
    stoVec[li_store_dev_i] = i0;
  }

  return true;
}

/**
 * Load the DAE the force derivative vector, dFdx
 *
 * Loads the F-vector contributions for a single resistor instance.
 *
 * @return true on success
 *
 * @author Eric Keiter, SNL, Parallel Computational Sciences
 * @date   03/05/04
 */
bool Instance::loadDAEdFdx()
{
  N_LAS_Matrix & dFdx = *(extData.dFdxMatrixPtr);
  dFdx[li_Pos][APosEquPosNodeOffset] += G;
  dFdx[li_Pos][APosEquNegNodeOffset] -= G;
  dFdx[li_Neg][ANegEquPosNodeOffset] -= G;
  dFdx[li_Neg][ANegEquNegNodeOffset] += G;
  return true;
}

/**
 * Update the temperature of the device
 *
 * @param temp_tmp temperature
 *
 * @return true on success
 *
 * @author Tom Russo, Component Information and Models
 * @date   02/27/01
 */
bool Instance::updateTemperature(const double & temp_tmp)
{
  bool bsuccess = true;
  double difference, factor;

  if (temp_tmp != -999.0)
    temp = temp_tmp;
  difference = temp - model_.tnom;
  factor = 1.0 + tempCoeff1*difference + tempCoeff2*difference*difference;

  if (R*factor != 0.0)
    G = 1.0/(R * factor);
  else
    G = 0.0;

  return bsuccess;
}

/**
 * Process model parameter
 *
 * @param param parameter to process
 *
 * @return true on success
 *
 * @author Eric Keiter, SNL, Parallel Computational Sciences
 * @date   6/03/02
 */
bool Model::processParams(string param)
{
  return true;
}

/**
 * Process the instance parameters of instance owned by this model
 *
 *
 * @param param
 *
 * @return
 *
 * @author Dave Shirely, PSSI
 * @date   03/23/06
 */
bool Model::processInstanceParams(string param)
{
  for (InstanceVector::const_iterator it = getInstanceVector().begin(); it != getInstanceVector().end(); ++it)
  {
    (*it)->processParams();
  }

  return true;
}

/**
 * Construct a resistor model
 *
 * @param model_block Model information from parser
 * @param solver_state Solution sover state
 * @param device_options Device options defined in netlist
 *
 * @author Eric Keiter, SNL, Parallel Computational Sciences
 * @date   5/16/00
 */
Model::Model(
  const ModelBlock &    model_block,
  SolverState &         solver_state,
  DeviceOptions &       device_options)
  : DeviceModel(model_block, solver_state, device_options),
    tempCoeff1(0.0),
    tempCoeff2(0.0),
    sheetRes(0.0),
    defWidth(10e-6),
    narrow(0.0),
    tnom(device_options.tnom)
{
  // Set params to constant default values.
  setDefaultParams();

  // Set params according to .model line and constant defaults from metadata.
  setModParams(model_block.params);

  // Set any non-constant parameter defaults.
  if (!given("TNOM"))
    tnom = device_options.tnom;

  // Calculate any parameters specified as expressions.
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors.
  processParams();
}

/**
 * Destroy this model.
 *
 * @author Eric Keiter, SNL, Parallel Computational Sciences
 * @date   3/16/00
 */
Model::~Model()
{
  // Destory all owned instances
  for (InstanceVector::const_iterator it = getInstanceVector().begin(); it != getInstanceVector().end(); ++it)
  {
    delete (*it);
  }
}

/**
 * Print instances owned by this model.
 *
 * For debugging
 *
 * @param os output stream
 *
 * @return reference to output stream
 *
 * @author Eric Keiter, SNL, Parallel Computational Sciences
 * @date   4/03/00
 */
std::ostream &Model::printOutInstances(std::ostream &os) const
{
  os << endl;
  os << "Number of Resistor Instances: " << instanceContainer.size() << endl;
  os << "    name     getModelName()  Parameters" << endl;

  int i = 0;
  for (InstanceVector::const_iterator it = instanceContainer.begin(); it != instanceContainer.end(); ++it)
  {
    os << "  " << i << ": " << (*it)->getName() << "\t";
    os << (*it)->getModelName();
    os << "\t\tR(Tnom) = " << (*it)->R;
    os << "\tG(T) = " << (*it)->G;
    os << endl;
    ++i;
  }

  os << endl;

  return os;
}

/**
 * Update state for all resistor instances, regardless of model.
 *
 * @param solVec solution vector
 * @param staVec state vector
 * @param stoVec store vector
 *
 * @return true on success
 *
 * @author Eric Keiter, SNL
 * @date   11/26/08
 */
bool Master::updateState(double * solVec, double * staVec, double * stoVec)
{
  for (InstanceVector::const_iterator it = getInstanceVector().begin(); it != getInstanceVector().end(); ++it)
  {
    Instance & ri = *(*it);

    double v_pos = solVec[ri.li_Pos];
    double v_neg = solVec[ri.li_Neg];

    // Load RHS vector element for the positive circuit node KCL equ.
    ri.i0 = (v_pos-v_neg)*ri.G;
  }

  return true;
}

/**
 * Load DAE vectors of all resistor instances, regardless of model
 *
 * @param solVec solution vector
 * @param fVec f vector
 * @param qVec q vector
 * @param storeLeadF store lead current f vector
 * @param storeLeadQ store lead current q vector
 *
 * @return true on success
 *
 * @author Eric Keiter, SNL
 * @date   11/26/08
 */
bool Master::loadDAEVectors (double * solVec, double * fVec, double *qVec,  double * storeLeadF, double * storeLeadQ)
{
  for (InstanceVector::const_iterator it = getInstanceVector().begin(); it != getInstanceVector().end(); ++it)
  {
    Instance & ri = *(*it);
    fVec[ri.li_Pos] += ri.i0;
    fVec[ri.li_Neg] += -ri.i0;
    if( ri.loadLeadCurrent )
    {
      storeLeadF[ri.li_store_dev_i] = ri.i0;
    }
  }
  return true;
}

/**
 * Load DAE matrices for all resistor instances, regardless of model
 *
 * @param dFdx derivative f matrix
 * @param dQdx derivative q matrix
 *
 * @return true on success
 *
 * @author Eric Keiter, SNL
 * @date   11/26/08
 */
bool Master::loadDAEMatrices(N_LAS_Matrix & dFdx, N_LAS_Matrix & dQdx)
{
  for (InstanceVector::const_iterator it = getInstanceVector().begin(); it != getInstanceVector().end(); ++it)
  {
    Instance & ri = *(*it);

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
    *(ri.f_PosEquPosNodePtr) += ri.G;
    *(ri.f_PosEquNegNodePtr) -= ri.G;
    *(ri.f_NegEquPosNodePtr) -= ri.G;
    *(ri.f_NegEquNegNodePtr) += ri.G;
#else
    dFdx[ri.li_Pos][ri.APosEquPosNodeOffset] += ri.G;
    dFdx[ri.li_Pos][ri.APosEquNegNodeOffset] -= ri.G;
    dFdx[ri.li_Neg][ri.ANegEquPosNodeOffset] -= ri.G;
    dFdx[ri.li_Neg][ri.ANegEquNegNodeOffset] += ri.G;
#endif
  }

  return true;
}

} // namespace Resistor
} // namespace Device
} // namespace Xyce

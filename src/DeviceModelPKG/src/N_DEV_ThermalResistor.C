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

//----------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_ThermalResistor.C,v $
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
// Revision Number: $Revision: 1.28.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:39 $
//
// Current Owner  : $Author: tvrusso $
//----------------------------------------------------------------------------
#include <Xyce_config.h>

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include <N_DEV_Const.h>
#include <N_DEV_ThermalResistor.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_MatrixLoadData.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>

#include <N_UTL_Expression.h>

namespace Xyce {
namespace Device {

template<>
ParametricData<ThermalResistor::Instance>::ParametricData()
{
    setNumNodes(2);
    setNumOptionalNodes(0);
    setNumFillNodes(0);
    setModelRequired(0);
    setPrimaryParameter("R");
    addModelType("R");

    // Set up double precision variables:
    addPar ("R", 1000.0, false, ParameterType::TIME_DEP,
      &ThermalResistor::Instance::R,
      NULL, U_OHM, CAT_NONE, "Resistance");

    addPar ("L", 0.0, false, ParameterType::NO_DEP,
      &ThermalResistor::Instance::length,
      NULL, U_METER, CAT_NONE, "Length of conductor");

    addPar ("W", 0.0, false, ParameterType::NO_DEP,
      &ThermalResistor::Instance::width,
      NULL, U_METER, CAT_NONE, "Width of conductor");

    addPar ("A", 0.0, false, ParameterType::NO_DEP,
      &ThermalResistor::Instance::area,
      NULL, U_METER2, CAT_NONE, "Area of conductor");

    addPar ("THERMAL_L", 0.0, false, ParameterType::NO_DEP,
      &ThermalResistor::Instance::thermalLength,
      NULL, U_METER, CAT_NONE, "Length of material thermally coupled to conductor");

    addPar ("THERMAL_A", 0.0, false, ParameterType::NO_DEP,
      &ThermalResistor::Instance::thermalArea,
      NULL, U_METER2, CAT_NONE, "Area of material thermally coupled to conductor");


     // This stuff is copied from the model:
    addPar ("RESISTIVITY", 0.0, false, ParameterType::TIME_DEP,
      &ThermalResistor::Instance::resistivity,
      NULL, U_OHMM, CAT_NONE, "Resistor material resistivity");

    addPar ("DENSITY", 0.0, false, ParameterType::TIME_DEP,
      &ThermalResistor::Instance::density,
      NULL, U_KGMM3, CAT_NONE, "Resistor material density (unused)");

    addPar ("HEATCAPACITY", 0.0, false, ParameterType::TIME_DEP,
      &ThermalResistor::Instance::heatCapacity,
      NULL, U_JMM3KM1, CAT_NONE, "Resistor material volumetric heat capacity");

    addPar ("THERMAL_HEATCAPACITY", 0.0, false, ParameterType::TIME_DEP,
      &ThermalResistor::Instance::thermalHeatCapacity,
      NULL, U_JMM3KM1, CAT_NONE, "Volumetric heat capacity of material thermally coupled to conductor");


    addPar ("TEMP", 0.0, false, ParameterType::TIME_DEP,
      &ThermalResistor::Instance::temp,
      NULL, U_DEGC, CAT_NONE, "Temperature");

    // Set up non-double precision variables:
    addPar ("OUTPUTINTVARS", false, false, ParameterType::NO_DEP,
      &ThermalResistor::Instance::outputInternalVarsFlag, NULL,
            U_NONE, CAT_CONTROL, "Debug Output switch");
}

template<>
ParametricData<ThermalResistor::Model>::ParametricData()
{
    addPar ("TC1",      0.0, false,   ParameterType::NO_DEP,
      &ThermalResistor::Model::tempCoeff1, NULL,
       U_DEGCM1, CAT_NONE, "Linear Temperature Coefficient");

    addPar ("TC2",      0.0, false,   ParameterType::NO_DEP,
      &ThermalResistor::Model::tempCoeff2, NULL,
       U_DEGCM2, CAT_NONE, "Quadratic Temperature Coefficient");

    addPar ("RSH",      0.0, false,   ParameterType::NO_DEP,
      &ThermalResistor::Model::sheetRes, NULL,
       U_OHM,  CAT_NONE, "Sheet Resistance");


    addPar ("RESISTIVITY", 0.0, false, ParameterType::TIME_DEP,
      &ThermalResistor::Model::resistivity, NULL,
       U_OHMM,  CAT_NONE,  "Resistor material resistivity");

    addPar ("DENSITY", 0.0, false, ParameterType::TIME_DEP,
      &ThermalResistor::Model::density, NULL,
       U_KGMM3,  CAT_NONE,  "Resistor material density (unused)");

    addPar ("HEATCAPACITY", 0.0, false, ParameterType::TIME_DEP,
      &ThermalResistor::Model::heatCapacity, NULL,
       U_JMM3KM1,  CAT_NONE,  "Resistor material volumetric heat capacity");

    addPar ("THERMAL_HEATCAPACITY", 0.0, false, ParameterType::TIME_DEP,
      &ThermalResistor::Model::thermalHeatCapacity, NULL,
       U_JMM3KM1,  CAT_NONE,  "Volumetric heat capacity of material thermally coupled to conductor");

    addPar ("DEFW",     1.e-5, false, ParameterType::NO_DEP,
      &ThermalResistor::Model::defWidth, NULL,
       U_METER,  CAT_NONE, "Default Instance Width");

    addPar ("NARROW",   0.0, false,   ParameterType::NO_DEP,
      &ThermalResistor::Model::narrow, NULL,
       U_METER,  CAT_NONE, "Narrowing due to side etching");

    addPar ("TNOM",     0.0, false,   ParameterType::NO_DEP,
      &ThermalResistor::Model::tnom, NULL,
       U_DEGC, CAT_NONE, "Parameter Measurement Temperature");
}

namespace ThermalResistor {

vector< vector<int> > Instance::jacStamp;



ParametricData<Instance> &Instance::getParametricData() {
  static ParametricData<Instance> parMap;

  return parMap;
}

ParametricData<Model> &Model::getParametricData() {
  static ParametricData<Model> parMap;

  return parMap;
}

// Class Instance
//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/03/02
//-----------------------------------------------------------------------------
bool Instance::processParams (string param)
{

  // now set the temperature related stuff.
  updateTemperature(temp);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : instance block constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
Instance::Instance(
  InstanceBlock & IB,
    Model & Riter,
    MatrixLoadData & mlData1,
    SolverState &ss1,
    ExternData  &ed1,
    DeviceOptions & do1)

  : DeviceInstance(IB, mlData1, ss1, ed1, do1),
    model_(Riter),
    R(0.0), G(0.0),
    i0(0.0),
    length(0.0),
    width(0.0),
    temp(getDeviceOptions().temp.dVal()),
    li_Pos(-1),
    li_Neg(-1),
    APosEquPosNodeOffset(-1),
    APosEquNegNodeOffset(-1),
    ANegEquPosNodeOffset(-1),
    ANegEquNegNodeOffset(-1),
    f_PosEquPosNodePtr(0),
    f_PosEquNegNodePtr(0),
    f_NegEquPosNodePtr(0),
    f_NegEquNegNodePtr(0),
    tempModelEnabled(false),
    outputInternalVarsFlag(false),
    li_TempState(-1),
    li_store_dev_i(-1)
{
  numIntVars   = 0;
  numExtVars   = 2;
  numStateVars = 0;
  numLeadCurrentStoreVars = 1; // one potential lead current 

  defaultParamName = "R";

  setName(IB.getName());
  setModelName(model_.getName());

  if( jacStamp.empty() )
  {
    jacStamp.resize(2);
    jacStamp[0].resize(2);
    jacStamp[1].resize(2);
    jacStamp[0][0] = 0;
    jacStamp[0][1] = 1;
    jacStamp[1][0] = 0;
    jacStamp[1][1] = 1;
  }


  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  // Set any non-constant parameter defaults:
  if (!given("TEMP"))
    temp = getDeviceOptions().temp.dVal();
  if (!given("W"))
    width = model_.defWidth;;

  // Nonzero value for numStateVars indicates that the self-consistent thermal
  // resistor model is being used.
  if (given("A") && given("L") &&
      ( (model_.given("HEATCAPACITY") && model_.given("RESISTIVITY"))
                || (given("HEATCAPACITY") && given("RESISTIVITY")) )
       )
  {
    numStateVars++;
    tempModelEnabled = true;
  }

  // If the instance parameters are NOT given, but the model parameters
  // ARE given, then copy the model params into the instance.  All the
  // work is actually done in the instance anyway,
  // but this stuff can be specified on the model level.
  if ( !(given("HEATCAPACITY")) && !(given("RESISTIVITY")) &&
      model_.given("HEATCAPACITY") && model_.given("RESISTIVITY")  )
  {
    resistivity = model_.resistivity;
    heatCapacity = model_.heatCapacity;
    thermalHeatCapacity = model_.thermalHeatCapacity;

    // copy over the dependent parameters.  For now, it only appears necessary
    // to copy the dependentParams vector, and not anything else like
    // the expVarLIDs vector.

    if (!(model_.dependentParams.empty()))
    {
      vector<sDepend> & model_dp = model_.dependentParams;
      int dpSize = model_dp.size();

      for (int i=0;i<dpSize;++i)
      {
        sDepend dpTmp;
        dpTmp.name = model_dp[i].name;
        dpTmp.vals = model_dp[i].vals;
        dpTmp.global_params = model_dp[i].global_params;
        dpTmp.n_vars = model_dp[i].n_vars;
        dpTmp.lo_var = model_dp[i].lo_var;
        dpTmp.vectorIndex = -1;

        // dpTmp needs to point to a copy of the original expression.
        dpTmp.expr = new N_UTL_Expression( *(model_dp[i].expr) );

        double *Dval;
        if (dpTmp.name=="RESISTIVITY")
        {
          //Dval = &Instance::resistivity;
          Dval = &resistivity;
          dpTmp.resultU.result = Dval;

        }

        if (dpTmp.name=="HEATCAPACITY")
        {
          Dval = &heatCapacity;
          dpTmp.resultU.result = Dval;
        }

        dependentParams.push_back(dpTmp);
      }
    }
  }

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:

  if (!given("R") && numStateVars == 0)
  {
    if (model_.given("RSH") && given("L") && (model_.sheetRes!=0) &&
        (length != 0))
    {
      R = model_.sheetRes * (length - model_.narrow)
        / (width - model_.narrow);
    }
    else
    {
      R = 1000;
      string msg="***********\n";
      msg += ": WARNING!  Resistance=0, ";
      msg += "set to default of 1000 ohms " + getName() + "\n";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_WARNING_0,msg);
    }
  }

  processParams ();
}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
Instance::~Instance ()
{
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/12/02
//-----------------------------------------------------------------------------
void Instance::registerLIDs( const vector<int> & intLIDVecRef,
                                           const vector<int> & extLIDVecRef )
{
  string msg;

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline =
  "-------------------------------------------------------------------------"
  "----";
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
    msg = "Instance::registerLIDs:";
    msg += "numInt != numIntVars";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
  }

  if (numExt != numExtVars)
  {
    msg = "Instance::registerLIDs:";
    msg += "numExt != numExtVars";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
  }

  // copy over the global ID lists.
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  li_Pos = extLIDVec[0];
  li_Neg = extLIDVec[1];

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 )
  {
    cout << "  li_Pos = " << li_Pos << endl;
    cout << "  li_Neg = " << li_Neg << endl;
  }
#endif

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 )
  {
    cout << dashedline << endl;
  }
#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       : Note that the resistor does not have any state vars.
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/12/02
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs(const vector<int> & staLIDVecRef )
{
  string msg;

  // Check if the size of the ID lists corresponds to the
  // proper number of internal and external variables.
  int numSta = staLIDVecRef.size();

  if (numSta != numStateVars)
  {
    msg = "Instance::registerLIDs:";
    msg += "numSta != numStateVars";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
  }

  // Copy over the global ID lists:
  staLIDVec = staLIDVecRef;

  if (numStateVars > 0)
    li_TempState = staLIDVec[0];

}

// ----------------------------------------------------------------------------
// Function      : Instance::registerStoreLIDs
// Purpose       : One store var for device current.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 04/24/2013
// ----------------------------------------------------------------------------
void Instance::registerStoreLIDs(const vector<int> & stoLIDVecRef)
{
  string msg;

  // Check if the size of the ID lists corresponds to the
  // proper number of internal and external variables.
  int numSto = stoLIDVecRef.size();

  if (numSto != getNumStoreVars())
  {
    msg = "Instance::registerStoreLIDs:";
    msg += "numSto != numStoreVars";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL,msg);
  }
  if( loadLeadCurrent )
  {
    li_store_dev_i = stoLIDVecRef[0];
  }
}

// ----------------------------------------------------------------------------
// Function      : Instance::getStoreNameMap
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 04/24/2013
// ----------------------------------------------------------------------------
map<int,string> & Instance::getStoreNameMap()
{
  // set up the internal name map, if it hasn't been already.
  if( loadLeadCurrent && storeNameMap.empty ())
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

//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/20/01
//-----------------------------------------------------------------------------
const vector< vector<int> > & Instance::jacobianStamp() const
{
  return jacStamp;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerJacLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/27/01
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const vector< vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );

  APosEquPosNodeOffset = jacLIDVec[0][0];
  APosEquNegNodeOffset = jacLIDVec[0][1];
  ANegEquPosNodeOffset = jacLIDVec[1][0];
  ANegEquNegNodeOffset = jacLIDVec[1][1];
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupPointers
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/12/08
//-----------------------------------------------------------------------------
void Instance::setupPointers ()
{
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  N_LAS_Matrix & dFdx = *(extData.dFdxMatrixPtr);
  f_PosEquPosNodePtr = &(dFdx[li_Pos][APosEquPosNodeOffset]);
  f_PosEquNegNodePtr = &(dFdx[li_Pos][APosEquNegNodeOffset]);
  f_NegEquPosNodePtr = &(dFdx[li_Neg][ANegEquPosNodeOffset]);
  f_NegEquNegNodePtr = &(dFdx[li_Neg][ANegEquNegNodeOffset]);
#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       : update intermediate variables for one diode instance
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, Dept. 9233.
// Creation Date : 3/05/04
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars ()
{
  double * solVec = extData.nextSolVectorRawPtr;

  if (tempModelEnabled)
  {
    N_LAS_Vector * staVectorPtr = extData.currStaVectorPtr;

    if (!getSolverState().dcopFlag)
    {
      if (li_TempState >= 0)
      {
        temp = (*staVectorPtr)[li_TempState];
        updateTemperature(temp);
      }
    }
  }

  double v_pos = solVec[li_Pos];
  double v_neg = solVec[li_Neg];

  // Load RHS vector element for the positive circuit node KCL equ.
  i0 = (v_pos-v_neg)*G;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/29/01
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState ()
{
  bool bsuccess = updateIntermediateVars ();

  if (tempModelEnabled)
  {
    double * staVec = extData.nextStaVectorRawPtr;
    if (li_TempState >= 0)
    {
      double dissipation = i0*i0*R;
      temp += dissipation*getSolverState().currTimeStep/(area*length*heatCapacity +
                            thermalArea*thermalLength*thermalHeatCapacity);
      staVec[li_TempState] = temp;
    }
  }

  return  bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::outputPlotFiles
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool Instance::outputPlotFiles ()
{
  bool bsuccess = true;

  if (tempModelEnabled && outputInternalVarsFlag)
  {
    cout.width(28); cout.precision(20); cout.setf(ios::scientific);
    N_LAS_Vector * sta1VectorPtr = extData.nextStaVectorPtr;
    N_LAS_Vector * sta2VectorPtr = extData.currStaVectorPtr;
    cout << "TEMP("<<getName()<<"):  " << getSolverState().currTime << "    "
        << ((*sta1VectorPtr)[li_TempState]-CONSTCtoK) << "    "
        << ((*sta2VectorPtr)[li_TempState]-CONSTCtoK)
        << endl;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 resistor  instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/24/03
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
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

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 resistor  instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  N_LAS_Matrix & dFdx = *(extData.dFdxMatrixPtr);

  dFdx[li_Pos][APosEquPosNodeOffset] += G;
  dFdx[li_Pos][APosEquNegNodeOffset] -= G;
  dFdx[li_Neg][ANegEquPosNodeOffset] -= G;
  dFdx[li_Neg][ANegEquNegNodeOffset] += G;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateTemperature
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, Component Information and Models
// Creation Date : 02/27/01
//-----------------------------------------------------------------------------
bool Instance::updateTemperature ( const double & temp_tmp)
{
  double difference, factor;

  if (tempModelEnabled)
  {
    updateDependentParameters(temp_tmp);
    R = resistivity * length / area;
    factor = 1;
  }
  else
  {
    difference = temp_tmp - model_.tnom;
    factor = 1.0 + (model_.tempCoeff1)*difference +
      (model_.tempCoeff2)*difference*difference;
  }

  if (R*factor != 0.0)
    G = 1.0/(R * factor);
  else
    G = 0.0;


  return true;
}

//-----------------------------------------------------------------------------
// Function      : Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/03/02
//-----------------------------------------------------------------------------
bool Model::processParams (string param)
{
  return true;
}

//----------------------------------------------------------------------------
// Function      : Model::processInstanceParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirely, PSSI
// Creation Date : 03/23/06
//----------------------------------------------------------------------------
bool Model::processInstanceParams(string param)
{
  vector<Instance*>::iterator iter;
  vector<Instance*>::iterator first = instanceContainer.begin();
  vector<Instance*>::iterator last  = instanceContainer.end();

  for (iter=first; iter!=last; ++iter)
  {
    (*iter)->processParams();
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : model block constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
Model::Model (const ModelBlock & MB,
                                                SolverState & ss1,
                                                DeviceOptions & do1)
  : DeviceModel(MB, ss1,do1),
    tempCoeff1(0.0),
    tempCoeff2(0.0),
    sheetRes(0.0),
    defWidth(10e-6),
    narrow(0.0),
    tnom(getDeviceOptions().tnom)
{

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to .model line and constant defaults from metadata:
  setModParams (MB.params);

  // Set any non-constant parameter defaults:
  if (!given("TNOM"))
    tnom = getDeviceOptions().tnom;

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  if (!given("THERMAL_HEATCAPACITY"))
    thermalHeatCapacity = heatCapacity;

  // calculate dependent (ie computed) params and check for errors:
  processParams();
}

//-----------------------------------------------------------------------------
// Function      : Model::~Model
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
Model::~Model ()
{
  vector<Instance*>::iterator iter;
  vector<Instance*>::iterator first = instanceContainer.begin();
  vector<Instance*>::iterator last  = instanceContainer.end();

  for (iter=first; iter!=last; ++iter)
  {
    delete (*iter);
  }
}

//-----------------------------------------------------------------------------
// Function      : Model::printOutInstances
// Purpose       : debugging tool.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/03/00
//-----------------------------------------------------------------------------
std::ostream &Model::printOutInstances(std::ostream &os) const
{
  vector<Instance*>::const_iterator iter;
  vector<Instance*>::const_iterator first = instanceContainer.begin();
  vector<Instance*>::const_iterator last  = instanceContainer.end();

  int i,isize;
  isize = instanceContainer.size();
  os << endl;
  os << "Number of Resistor Instances: " << isize << endl;
  os << "    name     getModelName()  Parameters" << endl;
  for (i=0, iter=first; iter!=last; ++iter, ++i)
  {
    os << "  " << i << ": " << (*iter)->getName() << "\t";
    os << (*iter)->getModelName();
    os << "\t\tR(Tnom) = " << (*iter)->R;
    os << "\tG(T) = " << (*iter)->G;
    os << endl;
  }

  os << endl;

  return os;
}

// ThermalResistor Master functions:

//-----------------------------------------------------------------------------
// Function      : Master::updateState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
bool Master::updateState (double * solVec, double * staVec, double * stoVec)
{
  bool bsuccess = true;

  // first loop over the models:
  for (ModelMap::const_iterator model_it = getModelMap().begin(); model_it != getModelMap().end(); ++model_it)
  {
    // loop over the instances for this model.
    InstanceVector::const_iterator first = (*model_it).second->instanceContainer.begin();
    InstanceVector::const_iterator last = (*model_it).second->instanceContainer.end();

    for (InstanceVector::const_iterator it = first; it != last; ++it)
    {
      Instance & ri = *(*it);

      bool btmp = ri.updateIntermediateVars ();
      bsuccess = bsuccess && btmp;

      if (ri.tempModelEnabled)
      {
        if (ri.li_TempState >= 0)
        {
          double dissipation = ri.i0*ri.i0*ri.R;
          double dt = ri.getSolverState().currTimeStep;
          ri.temp += dissipation*dt/(ri.area*ri.length*ri.heatCapacity +
                                ri.thermalArea*ri.thermalLength*ri.thermalHeatCapacity);
          staVec[ri.li_TempState] = ri.temp;
        }
      }
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Master::loadDAEVectors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
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

//-----------------------------------------------------------------------------
// Function      : Master::loadDAEMatrices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
bool Master::loadDAEMatrices (N_LAS_Matrix & dFdx, N_LAS_Matrix & dQdx)
{
#ifdef _OMP
#pragma omp parallel for
#endif
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

} // namespace ThermalResistor
} // namespace Device
} // namespace Xyce

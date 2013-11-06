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
// Filename       : $RCSfile: N_DEV_Capacitor.C,v $
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
// Revision Number: $Revision: 1.231.2.5 $
//
// Revision Date  : $Date: 2013/10/03 17:23:38 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

#include <N_UTL_Misc.h>

#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif

// ----------   Xyce Includes   ----------
#include <N_DEV_Const.h>
#include <N_DEV_Factory.h>
#include <N_DEV_Capacitor.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_MatrixLoadData.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>

#include <N_UTL_Expression.h>

namespace Xyce {
namespace Device {

template<>
ParametricData<Capacitor::Instance>::ParametricData()
{
  // Set up configuration constants:
  setNumNodes(2);
  setNumOptionalNodes(0);
  setNumFillNodes(0);
  setModelRequired(0);
  setPrimaryParameter("C");
  addModelType("C");

  // Set up double precision variables:
  addPar ("C", 1.e-6,false, ParameterType::SOLN_DEP, &Capacitor::Instance::C, NULL, U_FARAD, CAT_NONE, "Capacitance");
  addPar ("IC", 0.0, false, ParameterType::NO_DEP, &Capacitor::Instance::IC, &Capacitor::Instance::ICGiven, STANDARD, CAT_NONE, "");
  addPar ("L", 1.0, false, ParameterType::NO_DEP, &Capacitor::Instance::length, NULL, U_METER, CAT_NONE, "Semiconductor capacitor width");
  addPar ("W", 1.e-6,false, ParameterType::NO_DEP, &Capacitor::Instance::width, NULL, U_METER, CAT_NONE, "Semiconductor capacitor length");
  addPar ("AGE",0.0, false, ParameterType::NO_DEP, &Capacitor::Instance::age, NULL,U_HOUR, CAT_NONE, "Age of capacitor");
  addPar ("D", 0.0233, false, ParameterType::NO_DEP, &Capacitor::Instance::ageCoef, NULL, U_NONE, CAT_NONE, "Age degradation coefficient");
  addPar ("TEMP", 0.0, false,ParameterType::TIME_DEP, &Capacitor::Instance::temp, NULL, STANDARD, CAT_NONE, "");
  /* Genie 121412. Support TC for capacitor instances*/
  addPar ("TC1",   0.0, false,   ParameterType::NO_DEP, &Capacitor::Instance::tempCoeff1, &Capacitor::Instance::tempCoeff1Given, U_DEGCM1, CAT_NONE, "Linear Temperature Coefficient");
  addPar ("TC2",   0.0, false,   ParameterType::NO_DEP, &Capacitor::Instance::tempCoeff2, &Capacitor::Instance::tempCoeff2Given, U_DEGCM2, CAT_NONE, "Quadratic Temperature Coefficient");
  makeVector ("TC", 2);
}

template<>
ParametricData<Capacitor::Model>::ParametricData()
{
  // Set up double precision variables:
  addPar ("CJ", 0.0, false, ParameterType::NO_DEP, &Capacitor::Model::cj, NULL, U_FARADMM2, CAT_NONE, "Junction bottom capacitance");
  addPar ("CJSW", 0.0, false, ParameterType::NO_DEP, &Capacitor::Model::cjsw, NULL, U_FARADMM1, CAT_NONE, "Junction sidewall capacitance");
  addPar ("DEFW", 1.e-6, false, ParameterType::NO_DEP, &Capacitor::Model::defWidth, NULL, U_METER, CAT_NONE, "Default device width");
  addPar ("NARROW", 0.0, false, ParameterType::NO_DEP, &Capacitor::Model::narrow, NULL, U_METER, CAT_NONE, "Narrowing due to side etching");
  addPar ("TC1", 0.0, false, ParameterType::NO_DEP, &Capacitor::Model::tempCoeff1, NULL, STANDARD, CAT_NONE, "");
  addPar ("TC2", 0.0, false, ParameterType::NO_DEP, &Capacitor::Model::tempCoeff2, NULL, STANDARD, CAT_NONE, "");
  addPar ("TNOM", 0.0, false, ParameterType::NO_DEP, &Capacitor::Model::tnom, NULL, STANDARD, CAT_NONE, "");
}

namespace Capacitor {



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

  baseCap = C;
  if (!given("C") && given("AGE"))
  {
    string msg =
      "Age defined, but no base capacitance given.  Can't use age-aware with ";
    msg += "semiconductor capacitor options.";
    std::ostringstream oss;
    oss << "Error in " << netlistLocation() << "\n" << msg;
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL, oss.str());
  }

  // the age aware capacitor simply modifies the base capacitance.
  if (given("AGE") && age >= 1)
  {
    baseCap = baseCap*(1-ageCoef*log10(age));
  }

  if (!given("C") && !given("L"))
  {
    string msg = "Could find neither C parameter or L in instance.";
    std::ostringstream oss;
    oss << "Error in " << netlistLocation() << "\n" << msg;
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL, oss.str());
  }

  // Now we know we have either cap or length specified.
  if (!given("C"))
  {
#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      cout << "Doing semiconductor capacitor! "<< endl;
      cout << " cj = " << model_.cj << endl;
      cout << " cjsw = " << model_.cjsw << endl;
      cout << " width = " << width << endl;
      cout << " length = " << length << endl;
      cout << " narrow = " << model_.narrow << endl;
    }
#endif

    baseCap = C =
              model_.cj*(length-model_.narrow)*(width-model_.narrow) +
              2*model_.cjsw*(length+width-2*model_.narrow);
  }

  // If there are any time dependent parameters, set their values for
  // the current time.

  updateTemperature(temp);

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
  bool bsuccess = true;
  double difference, factor;

  difference = temp - model_.tnom;
  //factor = 1.0 + (model_.tempCoeff1)*difference +
  //(model_.tempCoeff2)*difference*difference;
  factor = 1.0 + tempCoeff1*difference +  //support specifying TC at the instance line
           tempCoeff2*difference*difference;
  C = baseCap*factor;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << "Instance::updateTemperature.  C = " << C << endl;
    cout << "   temp = " << temp << endl;
    cout << "   temp_tmp = " << temp_tmp << endl;
    cout << "   tnom = " << model_.tnom << endl;
    cout << "   difference = " << difference << endl;
    cout << "   tempCoeff1 = " << tempCoeff1 << endl;
    cout << "   tempCoeff2 = " << tempCoeff2 << endl;
    cout << "   baseCap = " << baseCap << endl;
    cout << "   factor = " << factor  << endl;
  }
#endif

  return bsuccess;
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
  InstanceBlock &      IB,
  Model & Citer,
  MatrixLoadData &     mlData1,
  SolverState &        ss1,
  ExternData &         ed1,
  DeviceOptions &            do1)
  : DeviceInstance(IB, mlData1, ss1, ed1, do1),
    model_(Citer),
    expNumVars(0),
    expPtr(0),
    baseCap(0.0),
    temp(getDeviceOptions().temp.dVal()),
    tempGiven(0),
    tempCoeff1(0.0),
    tempCoeff2(0.0),
    tempCoeff1Given(false),
    tempCoeff2Given(false),
    li_Pos(-1),
    li_Neg(-1),
    li_QState(-1),
    li_vcapState(-1),
    li_capState(-1),
    li_store_dev_i(-1),
    APosEquPosNodeOffset(-1),
    ANegEquPosNodeOffset(-1),
    APosEquNegNodeOffset(-1),
    ANegEquNegNodeOffset(-1),
    APosEquBraNodeOffset(-1),
    ANegEquBraNodeOffset(-1),
    ABraEquBraNodeOffset(-1),
    ABraEquPosNodeOffset(-1),
    ABraEquNegNodeOffset(-1),

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
    qPosEquPosNodePtr(0),
    qNegEquPosNodePtr(0),
    qPosEquNegNodePtr(0),
    qNegEquNegNodePtr(0),
    fPosEquBraNodePtr(0),
    fNegEquBraNodePtr(0),
    fBraEquBraNodePtr(0),
    fBraEquPosNodePtr(0),
    fBraEquNegNodePtr(0),
#endif

    ICGiven(false),
    solVarDepC(false),
    IC(0)

{
  numIntVars   = 0;
  numExtVars   = 2;
  numStateVars = 1;
  setNumStoreVars(0);
  numLeadCurrentStoreVars = 1; // lead current DEV_I

  devConMap[0] = 1;
  devConMap[1] = 2;

  defaultParamName = "C";

  setName(IB.getName());
  setModelName(model_.getName());

  if( jacStamp.empty() )
  {
    jacStamp_IC.resize(3);
    jacStamp_IC[0].resize(3);
    jacStamp_IC[1].resize(3);
    jacStamp_IC[2].resize(3);
    jacStamp_IC[0][0] = 0;
    jacStamp_IC[0][1] = 1;
    jacStamp_IC[0][2] = 2;
    jacStamp_IC[1][0] = 0;
    jacStamp_IC[1][1] = 1;
    jacStamp_IC[1][2] = 2;
    jacStamp_IC[2][0] = 0;
    jacStamp_IC[2][1] = 1;
    jacStamp_IC[2][2] = 2;

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

  //  Set any non-constant parameter defaults:

  if (!given("W"))
    width = model_.defWidth;
  if (!given("TEMP"))
    temp = getDeviceOptions().temp.dVal();

  if (!tempCoeff1Given)
    tempCoeff1=model_.tempCoeff1;
  if (!tempCoeff2Given)
    tempCoeff2=model_.tempCoeff2;



  // Handle case where capacitance is solution-variable dependent:
  if (dependentParams.size()>0)
  {
    vector<sDepend>::iterator d;
    vector<sDepend>::iterator begin=dependentParams.begin();
    vector<sDepend>::iterator end=dependentParams.end();

    for (d=begin; d!=end; ++d)
    {
      if (d->name != "C")
      {
        string msg="Error: Solution-variable-dependent parameter other than C detected for device ";
        msg += getName();
        std::ostringstream oss;
        oss << "Error in " << netlistLocation() << "\n" << msg;
        N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, oss.str());
      }
      else
      {
        expNumVars = d->n_vars;
        expPtr = d->expr;
        solVarDepC = true;
        // To do the proper integration of the charge, we need to save the
        // voltage drop, the old capacitance and
        // the derivatives of Q and C from the last step.
        numStateVars += 2+2*expNumVars;

        if (expPtr->getNumDdt() != 0)
        {
          string msg="Error: Solution-variable-dependent expression for device ";
          msg += getName();
          msg += " contains time derivatives.  Feature not supported.\n";
          std::ostringstream oss;
          oss << "Error in " << netlistLocation() << "\n" << msg;
          N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, oss.str());
        }

#ifdef Xyce_DEBUG_DEVICE
        if (getDeviceOptions().debugLevel > 0 )
        {
          cout << " Instance::Instance:\n";
          cout << "   name=" << getName() << endl;
          cout << " Found solution-dependent parameter C depending on "
               << expNumVars << " variables." << endl;
        }
#endif

        // We now need to extend the pos and neg rows of the jacstamps
        // to account for the additional dependencies:

        jacStamp[0].resize(2+expNumVars);
        jacStamp[1].resize(2+expNumVars);
        jacStamp_IC[0].resize(3+expNumVars);
        jacStamp_IC[1].resize(3+expNumVars);
        for (int i=0; i<expNumVars; ++i)
        {
          jacStamp[0][2+i]=2+i;
          jacStamp[1][2+i]=2+i;

          jacStamp_IC[0][3+i]=3+i;
          jacStamp_IC[1][3+i]=3+i;
        }

        // finally, allocate space to hold the derivatives of C w.r.t.
        // the variables it depends on:
        expVarDerivs.resize(expNumVars);
        // and LIDs for state vector
        li_dQdXState.resize(expNumVars);
        li_dCdXState.resize(expNumVars);

      }
    }
  }

  // Calculate any parameters specified as expressions:

  updateDependentParameters();

  // calculate dependent (ie computed) params:

  processParams ();

  // we're gonna have to fake a voltage source at the operating point
  if (ICGiven ) numIntVars = 1;
}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
Instance::~Instance()
{
}

// Additional Declarations

//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/20/02
//-----------------------------------------------------------------------------
void Instance::registerLIDs( const vector<int> & intLIDVecRef,
                             const vector<int> & extLIDVecRef)

{
  string msg;

  // Check if the size of the ID lists corresponds to the proper number of
  // internal and external variables.
  int numInt = intLIDVecRef.size();
  int numExt = extLIDVecRef.size();

  if (numInt != numIntVars)
  {
    msg = "Instance::registerLIDs:";
    msg += "numInt != numIntVars";
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL,msg);
  }

  if (numExt != numExtVars)
  {
    msg = "Instance::registerLIDs:";
    msg += "numExt != numExtVars";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
  }

  // Copy over the local ID lists:
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  // Now use these lists to obtain the indices into the linear algebra
  // entities.  This assumes an order.  For the matrix indices, first do the
  // rows.

  li_Pos = extLIDVec[0];
  li_Neg = extLIDVec[1];

  // For fake voltage source at operating point
  if (ICGiven)
  {
    li_Bra = intLIDVec[0];
  }

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline =
    "-----------------------------------------------------------------------------";
  if (getDeviceOptions().debugLevel > 0 )
  {
    cout << dashedline << endl;

    cout << "::registerLIDs:\n";
    cout << "  name = " << getName() << endl;

    cout << "\nsolution indices:\n";
    cout << "  li_Pos = "<< li_Pos << endl;
    cout << "  li_Neg = "<< li_Neg<< endl;
    if (ICGiven)
      cout << "  li_Bra = "<< li_Bra<< endl;

    cout << dashedline << endl;
  }
#endif

}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/20/02
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs( const vector<int> & staLIDVecRef)
{
  string msg;

  // Check if the size of the ID lists corresponds to the proper number of
  // internal and external variables.
  int numSta = staLIDVecRef.size();
  int i=0;

  if (numSta != numStateVars)
  {
    msg = "Instance::registerStateLIDs:";
    msg += "numSta != numStateVars";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
  }

  // Copy over the global ID lists:
  staLIDVec = staLIDVecRef;

  li_QState = staLIDVec[i++];

  // If the capacitance is voltage dependent, we have additional state  vars
  if (solVarDepC)
  {
    li_vcapState = staLIDVec[i++];
    li_capState = staLIDVec[i++];

    for (int j = 0; j<expNumVars; ++j)
    {
      li_dQdXState[j] = staLIDVec[i++];
    }

    for (int j = 0; j<expNumVars; ++j)
    {
      li_dCdXState[j] = staLIDVec[i++];
    }
  }

}


//-----------------------------------------------------------------------------
// Function      : Instance::registerStoreLIDs
// Purpose       : One store var for device current.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 01/23/2012
//-----------------------------------------------------------------------------
void Instance::registerStoreLIDs(const vector<int> & stoLIDVecRef )
{
  string msg;

  // Check if the size of the ID lists corresponds to the
  // proper number of internal and external variables.
  int numSto = stoLIDVecRef.size();

  if (numSto != getNumStoreVars())
  {
    msg = "Instance::registerStoreLIDs:";
    msg += "numSto != numStoreVars";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
  }
  if( loadLeadCurrent )
  {
    li_store_dev_i = stoLIDVecRef[0];
  }
}


//-----------------------------------------------------------------------------
// Function      : Instance::getIntNameMap
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/13/05
//-----------------------------------------------------------------------------
map<int,string> & Instance::getIntNameMap ()
{
  // set up the internal name map, if it hasn't been already.
  if (ICGiven)
  {
    if (intNameMap.empty ())
    {
      string tmpstr;
      tmpstr = getName()+"_branch";
      spiceInternalName (tmpstr);
      intNameMap[li_Bra] = tmpstr;
    }
  }

  return intNameMap;
}

//-----------------------------------------------------------------------------
// Function      : Instance::getStoreNameMap
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Electrical Systems Modeling
// Creation Date : 01/23/2013
//-----------------------------------------------------------------------------
map<int,string> & Instance::getStoreNameMap ()
{
  if( loadLeadCurrent && storeNameMap.empty () )
  {
    string tmpstr(getName());
    spiceInternalName (tmpstr);
    tmpstr = getName()+":DEV_I";
    storeNameMap[li_store_dev_i] = tmpstr;
  }

  return storeNameMap;
}


//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/27/02
//-----------------------------------------------------------------------------
const vector< vector<int> > & Instance::jacobianStamp() const
{
  if (ICGiven)
    return jacStamp_IC;
  else
    return jacStamp;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerJacLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/27/02
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const vector< vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );

  APosEquPosNodeOffset = jacLIDVec[0][0];
  APosEquNegNodeOffset = jacLIDVec[0][1];
  ANegEquPosNodeOffset = jacLIDVec[1][0];
  ANegEquNegNodeOffset = jacLIDVec[1][1];

  if (ICGiven)
  {
    APosEquBraNodeOffset = jacLIDVec[0][2];
    ANegEquBraNodeOffset = jacLIDVec[1][2];
    ABraEquPosNodeOffset = jacLIDVec[2][0];
    ABraEquNegNodeOffset = jacLIDVec[2][1];
    ABraEquBraNodeOffset = jacLIDVec[2][2];
  }
  // set offsets if we have a solution-variable dependent C
  if (solVarDepC)
  {
    int depVarsBaseIndex = 2;
    if (ICGiven)
    {
      depVarsBaseIndex=3;
    }

    APosEquDepVarOffsets.resize(expNumVars);
    ANegEquDepVarOffsets.resize(expNumVars);

    for ( int i=0; i<expNumVars; ++i)
    {
      APosEquDepVarOffsets[i] = jacLIDVec[0][depVarsBaseIndex+i];
      ANegEquDepVarOffsets[i] = jacLIDVec[1][depVarsBaseIndex+i];
    }
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << "-----------------------------------------\n";
    cout << "Instance::registerJacLIDs\n";
    cout << " APosEquPosNodeOffset: " << APosEquPosNodeOffset << endl;
    cout << " APosEquNegNodeOffset: " << APosEquNegNodeOffset << endl;
    cout << " ANegEquPosNodeOffset: " << ANegEquPosNodeOffset << endl;
    cout << " ANegEquNegNodeOffset: " << ANegEquNegNodeOffset << endl;
    cout << " APosEquBraNodeOffset: " << APosEquBraNodeOffset << endl;
    cout << " ANegEquBraNodeOffset: " << ANegEquBraNodeOffset << endl;
    cout << " ABraEquPosNodeOffset: " << ABraEquPosNodeOffset << endl;
    cout << " ABraEquNegNodeOffset: " << ABraEquNegNodeOffset << endl;
    cout << " ABraEquBraNodeOffset: " << ABraEquBraNodeOffset << endl;
    if (solVarDepC)
    {
      for ( int i=0; i<expNumVars; ++i)
      {
        cout << " APosEquDepVarOffsets["<<i<<"]: " << APosEquDepVarOffsets[i]
             << endl;
        cout << " ANegEquDepVarOffsets["<<i<<"]: " << ANegEquDepVarOffsets[i]
             << endl;
      }
    }
    cout << "-----------------------------------------\n";
  }
#endif

}

//-----------------------------------------------------------------------------
// Function      : Instance::setupPointers
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/30/08
//-----------------------------------------------------------------------------
void Instance::setupPointers ()
{
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  N_LAS_Matrix & dFdx = *(extData.dFdxMatrixPtr);
  N_LAS_Matrix & dQdx = *(extData.dQdxMatrixPtr);

  qPosEquPosNodePtr = &(dQdx[li_Pos][APosEquPosNodeOffset]);
  qPosEquNegNodePtr = &(dQdx[li_Pos][APosEquNegNodeOffset]);
  qNegEquPosNodePtr = &(dQdx[li_Neg][ANegEquPosNodeOffset]);
  qNegEquNegNodePtr = &(dQdx[li_Neg][ANegEquNegNodeOffset]);

  if (solVarDepC)
  {
    qPosEquDepVarsPtrs.resize(expNumVars);
    qNegEquDepVarsPtrs.resize(expNumVars);

    for (int i=0; i<expNumVars; ++i)
    {
      qPosEquDepVarsPtrs[i]=&(dQdx[li_Pos][APosEquDepVarOffsets[i]]);
      qNegEquDepVarsPtrs[i]=&(dQdx[li_Neg][ANegEquDepVarOffsets[i]]);
    }
  }

  // there are no contributions to the dFdx matrix from dependent C's, so
  // we don't bother with those pointers.

  if (ICGiven)
  {
    fPosEquBraNodePtr = &(dFdx[li_Pos][APosEquBraNodeOffset]);
    fNegEquBraNodePtr = &(dFdx[li_Neg][ANegEquBraNodeOffset]);
    fBraEquPosNodePtr = &(dFdx[li_Bra][ABraEquPosNodeOffset]);
    fBraEquNegNodePtr = &(dFdx[li_Bra][ABraEquNegNodeOffset]);
    fBraEquBraNodePtr = &(dFdx[li_Bra][ABraEquBraNodeOffset]);
  }
#endif
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
  double * solVec = extData.nextSolVectorRawPtr;
  double * staVec = extData.nextStaVectorRawPtr;
  double v_pos = solVec[li_Pos];
  double v_neg = solVec[li_Neg];

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << " ----------------------------------" << endl;
    cout << "Instance::updatePrimaryState:" << endl;
  }
#endif

  vcap = v_pos-v_neg;

  if( getSolverState().dcopFlag && ICGiven ) vcap = IC;

  if (!solVarDepC)
  {
    // Obtain the "current"  value for the charge stored in the capacitor.
    q0 = C*vcap;
    staVec[li_QState] = q0;
  }
  else
  {
    // For solution-variable dependent cap, can't use the same formulation.
    // The capacitance is *strictly* dQ/dV, so must integrate CdV to get
    // charge.  We do this by incrementally adding C'*deltaV as V changes.
    // C' is the average capacitance between this and the previous step.
    // Using the average assures charge conservation, at least when C is
    // a function of vcap alone.
    // When C is not a function of vcap alone, we have additional derivative
    // terms that must also be integrated.  Ick.

    if (getSolverState().dcopFlag)
    {
      q0 = vcap*C;
    }
    else
    {
      double * oldstaVector = extData.currStaVectorRawPtr;
      double oldC;
      double oldVcap;
      q0=oldstaVector[li_QState];
      oldC=oldstaVector[li_capState];
      oldVcap=oldstaVector[li_vcapState];

      q0 += 0.5*(oldC+C)*(vcap-oldVcap);
    }
    staVec[li_QState] = q0;
    staVec[li_vcapState]=vcap;
    staVec[li_capState]=C;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/29/01
//-----------------------------------------------------------------------------
bool Instance::updateSecondaryState ()
{
  double * staDerivVec = extData.nextStaDerivVectorRawPtr;
  double * staVec = extData.nextStaVectorRawPtr;
  double * stoVec = extData.nextStoVectorRawPtr;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << " ----------------------------------" << endl;
    cout << "Instance::updateSecondaryState:" << endl;
  }
#endif

  if (solVarDepC)
  {
    // In other devices, at this point we'd handle all the DDT issues in
    // the expression before doing the evaluation, then redo the computation
    // of C with the updated time derivatives.  But for now, we're going to
    // disallow use of ddt in C, so we don't have to.  This restriction is
    // enforced in the constructor, where we throw a fatal error if the user
    // gives us a ddt-dependent C.
    // So here, we're just calling evaluate to get the derivatives of C,
    // and discarding the actual value of C (which has already been stored)
    // This code is in updateSecondaryState merely for consistency with
    // other devices that put similar code here (and because it'll need to be
    // here when/if we implement "ddt" handling.
    double junk;
    expPtr->evaluate(junk,expVarDerivs);

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      cout << "  Derivatives of C w.r.t variables: " << endl;
      for (int i=0;i<expNumVars;++i)
      {
        cout << " expVarDerivs[ "<< i << " ] = " << expVarDerivs[i] << endl;
      }
    }
#endif

    // Now we have some trickiness because we are integrating CdV to get Q.
    // In order to get dQ/dX we have two cases:
    //   X is one of our capacitor's voltage nodes:  dQ/dX = C or -C depending
    //       on whether X is the pos or negative node
    //  X is NOT one of our nodes: dQ/dX = integral( dC/dX *dV)

    // For now, because we don't have an easy way to tell which of our
    // expression nodes is which, we'll just calculate all the dC/dX and dQ/dX
    // the same way.  When it comes time to assemble the jacobian, we'll use
    // the node/equation offsets to know when to skip adding in this component.

    // The logic here is similar to the computation of the charge itself.
    // We'll use "expVarDerivs" to hold the final dQ/dX values.


    // Need to save the dC/dX value for next step.
    for (int i=0;i<expNumVars; ++i)
    {
      staVec[li_dCdXState[i]] = expVarDerivs[i];
    }

    if (getSolverState().dcopFlag)
    {
      // dQ/dX is just dC/dX*vcap
      for (int i=0;i<expNumVars; ++i)
      {
        expVarDerivs[i] *= vcap;
      }
    }
    else
    {
      double * oldstaVector = extData.currStaVectorRawPtr;
      // otherwise we have to integrate

      // dQ/dx = olddQdX + .5*(olddCdX+newdCdX)*(vcap-oldvcap)

      for (int i=0; i< expNumVars; ++i)
      {
        expVarDerivs[i] = oldstaVector[li_dQdXState[i]]
                          + 0.5*(oldstaVector[li_dCdXState[i]]+expVarDerivs[i])*
                          (vcap-oldstaVector[li_vcapState]);

      }
    }
    // Regardless of whether it's dcop or not, expVarDerivs now contains
    // dQ/dX for all the X's.  It's WRONG if X is one of our voltage nodes,
    // so we have to be careful not to use it in that case.
    // Save to state:
    for (int i=0;i<expNumVars; ++i)
    {
      staVec[li_dQdXState[i]] = expVarDerivs[i];
    }

  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEQVector
//
// Purpose       : Loads the Q-vector contributions for a single
//                 capacitor instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/24/03
//-----------------------------------------------------------------------------
bool Instance::loadDAEQVector ()
{
  (*extData.daeQVectorPtr)[li_Pos] += q0;
  (*extData.daeQVectorPtr)[li_Neg] += -q0;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 capacitor instance.
//
// Special Notes : See the special notes for loadDAEFVector.
//
//                 For the capacitor this doesn't do anything, except in
//                 the case of IC= being specified.  In that case, then
//                 some extra stuff is contributed that doesn't have time
//                 derivatives.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/24/03
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  bool bsuccess = true;
  double Vpos = 0.0;
  double Vneg = 0.0;
  double v_tmp = 0.0;
  double * fVec = extData.daeFVectorRawPtr;

  if (ICGiven && getSolverState().dcopFlag)
  {
    // If we're doing the operating point and we have an initial condition,
    // get the current from the branch equation
    Vpos = (*extData.nextSolVectorPtr)[li_Pos];
    Vneg = (*extData.nextSolVectorPtr)[li_Neg];

    // load current into the F vector
    fVec[li_Pos] += (*extData.nextSolVectorPtr)[li_Bra];
    fVec[li_Neg] += -(*extData.nextSolVectorPtr)[li_Bra];

    if( loadLeadCurrent )
    {
      double * stoVec = extData.nextStoVectorRawPtr;
      stoVec[li_store_dev_i] = (*extData.nextSolVectorPtr)[li_Bra];
    }
  }

  // Initial condition stuff.
  v_tmp=0;
  if (ICGiven && getSolverState().dcopFlag)
  {
    v_tmp= (Vpos-Vneg-IC);
  }

  // Do this whenever there's a Branch equation, but only if there is one.
  // We'll be using 0 if we're not the OP.
  if (ICGiven)
  {
    fVec[li_Bra] += v_tmp;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx
//
// Purpose       : Loads the dQdx-matrix contributions for a single
//                 capacitor instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  if (!(ICGiven&& getSolverState().dcopFlag))
  {
    N_LAS_Matrix & dQdx = *(extData.dQdxMatrixPtr);
    dQdx[li_Pos][APosEquPosNodeOffset] += C;
    dQdx[li_Pos][APosEquNegNodeOffset] -= C;
    dQdx[li_Neg][ANegEquPosNodeOffset] -= C;
    dQdx[li_Neg][ANegEquNegNodeOffset] += C;


    // Remember the comments in updateSecondaryState:
    // expVarDerivs contains dQ/dX, but only when X is not one of our
    // nodal voltages.  If X *IS* one of our nodal voltages, dQ/dX is either
    // C or -C and is already handled above.  We need only do the stuff
    // below for the dependencies on voltages that are NOT our nodal
    // voltages.
    if (solVarDepC)
    {
      for (int i=0; i< expNumVars; ++i)
      {
        if ( (APosEquDepVarOffsets[i] != APosEquPosNodeOffset)
             && (APosEquDepVarOffsets[i] != APosEquNegNodeOffset))
        {
          dQdx[li_Pos][APosEquDepVarOffsets[i]] += expVarDerivs[i];
        }
        if ( (ANegEquDepVarOffsets[i] != ANegEquPosNodeOffset)
             && (ANegEquDepVarOffsets[i] != ANegEquNegNodeOffset))
        {
          dQdx[li_Neg][ANegEquDepVarOffsets[i]] -= expVarDerivs[i];
        }
      }
    }

  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 capacitor instance.
//
// Special Notes : See the special notes for loadDAEFVector.
//
//                 For the capacitor this doesn't do anything, unless IC=
//                 has been specified for an initial condition.  Then,
//                 there are extra equations that do not contain time
//                 derivatives.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  N_LAS_Matrix & dFdx = *(extData.dFdxMatrixPtr);

  if (ICGiven && getSolverState().dcopFlag)
  {
    // Special Jacobian if we're doing operating point when IC given
    dFdx[li_Pos][APosEquBraNodeOffset] +=  1.0;
    dFdx[li_Neg][ANegEquBraNodeOffset] += -1.0;
    dFdx[li_Bra][ABraEquPosNodeOffset] +=  1.0;
    dFdx[li_Bra][ABraEquNegNodeOffset] += -1.0;
  }
  else
  {
    if (ICGiven)
      dFdx[li_Bra][ABraEquBraNodeOffset] += 1.0;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setIC
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/10/02
//-----------------------------------------------------------------------------
bool Instance::setIC ()
{
  double Vic;
  double Vpos, Vneg;
  double * nextStaVector = extData.nextStaVectorRawPtr;
  double * currStaVector = extData.currStaVectorRawPtr;

  double * nextSolVector = extData.nextSolVectorRawPtr;
  double * currSolVector = extData.currSolVectorRawPtr;

  if (ICGiven)
  {
    Vic = IC;
    q0 = C*Vic;

    currStaVector[li_QState] = q0;
    nextStaVector[li_QState] = q0;

    Vneg = currSolVector[li_Neg];
    Vpos = Vneg + Vic;

    currSolVector[li_Pos] = Vpos;
    nextSolVector[li_Pos] = Vpos;
    currSolVector[li_Neg] = -Vpos;
    nextSolVector[li_Neg] = -Vpos;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::varTypes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/17/04
//-----------------------------------------------------------------------------
void Instance::varTypes( vector<char> & varTypeVec )
{
  if (ICGiven)
  {
    varTypeVec.resize(1);
    varTypeVec[0] = 'I';
  }
}

// Class Model

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

  if (!tnomGiven)
    tnom = getDeviceOptions().tnom;
  else
    tnom += CONSTCtoK; // if user-specified, assume it was in deg. C.

  // If there are any time dependent parameters, set their values for
  // the current time.

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
// Purpose       : block constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/17/00
//-----------------------------------------------------------------------------

Model::Model(const ModelBlock & MB,
             SolverState & ss1,
             DeviceOptions & do1)
  : DeviceModel(MB,ss1,do1),
    cj(0.0),
    cjsw(0.0),
    defWidth(10e-6),
    narrow(0),
    tempCoeff1(0.0),
    tempCoeff2(0.0),
    tnom(getDeviceOptions().tnom),
    tnomGiven(0)

{

  // Set params to constant default values :
  setDefaultParams ();

  // Set params according to .model line and constant defaults from metadata:
  setModParams (MB.params);

  // Set any non-constant parameter defaults:
  if (!given("TNOM"))
    tnom = getDeviceOptions().tnom;

  // Calculate any parameters specified as expressions:

  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:

  processParams ();
}

//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------

Model::~Model()
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
  os << "Number of capacitor instances: " << isize << endl;
  os << "    name\t\tmodelName\tParameters" << endl;

  for (i = 0, iter = first; iter != last; ++iter, ++i)
  {
    os << "  " << i << ": " << (*iter)->getName() << "\t";
    os << (*iter)->getModelName();
    os << "\t\tC = " << (*iter)->C;
    os << "\tIC = " << (*iter)->IC;
    os << endl;
  }

  os << endl;

  return os;
}

// Capacitor Master functions:

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
  for (InstanceVector::const_iterator it = getInstanceVector().begin(); it != getInstanceVector().end(); ++it)
  {
    Instance & ci = *(*it);

    double v_pos = solVec[ci.li_Pos];
    double v_neg = solVec[ci.li_Neg];
    ci.vcap = v_pos-v_neg;

    if( getSolverState().dcopFlag && ci.ICGiven )
    {
      ci.vcap = ci.IC;
    }

    if (!ci.solVarDepC)
    {
      // Obtain the "current"  value for the charge stored in the capacitor.
      ci.q0 = ci.C * ci.vcap;
      staVec[ci.li_QState] = ci.q0;
    }
    else
    {
      // fall back on old pre-turbo scheme
      bool tmpBool=true;
      tmpBool = ci.updatePrimaryState ();
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Master::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
bool Master::updateSecondaryState ( double * staDerivVec, double * stoVec )
{
  for (InstanceVector::const_iterator it = getInstanceVector().begin(); it != getInstanceVector().end(); ++it)
  {
    Instance & ci = *(*it);

    if (ci.solVarDepC)
    {
      // fall back on old-style code if we've got one of these
      ci.updateSecondaryState();
    }
  }

  return true;
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

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << " ----------------------------------" << endl;
    cout << " Master::loadDAEVectors: " << endl;
  }
#endif

  for (InstanceVector::const_iterator it = getInstanceVector().begin(); it != getInstanceVector().end(); ++it)
  {
    Instance & ci = *(*it);
    if (ci.ICGiven)
    {
      double Vpos (0.0), Vneg (0.0), v_tmp (0.0);

      // Initial condition
      if (getSolverState().dcopFlag)
      {
#ifdef Xyce_DEBUG_DEVICE
        if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
        {
          cout << " loading dcop F vector for cap " << ci.getName() << ":" << endl;
        }
#endif
        // If doing the DCOP and have IC=, get current from branch equation
        // ci.i0   = solVec[ci.li_Bra]; moved to CapacitorMaster::updateState() where stovec is passed in.
        Vpos    = solVec[ci.li_Pos];
        Vneg    = solVec[ci.li_Neg];
        fVec [ci.li_Pos] += solVec[ci.li_Bra];
        fVec [ci.li_Neg] += -solVec[ci.li_Bra];

        if( ci.loadLeadCurrent )
        {
          storeLeadF[ci.li_store_dev_i] = solVec[ci.li_Bra];
        }

        if( ci.loadLeadCurrent )
        {
          storeLeadF[ci.li_store_dev_i] = solVec[ci.li_Bra];
        }

        if( ci.loadLeadCurrent )
        {
          storeLeadF[ci.li_store_dev_i] = solVec[ci.li_Bra];
        }

#ifdef Xyce_DEBUG_DEVICE
        if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
        {
          cout << " f[ " << ci.li_Pos << " ] += " << solVec[ci.li_Bra]<< endl;
          cout << " f[ " << ci.li_Neg << " ] += " << -solVec[ci.li_Bra] << endl;
        }
#endif

        v_tmp= (Vpos-Vneg-ci.IC);
      }

      // Do this only if there's a Branch equation
      fVec[ci.li_Bra] += v_tmp;

#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
      {
        cout << " f[ " << ci.li_Bra << " ] += " << v_tmp << endl;
      }
#endif
    }

    qVec[ci.li_Pos] += ci.q0;
    qVec[ci.li_Neg] += -ci.q0;

    if( ci.loadLeadCurrent )
    {
      storeLeadQ[ci.li_store_dev_i] = ci.q0;
      storeLeadF[ci.li_store_dev_i] = 0;
    }

    if( ci.loadLeadCurrent )
    {
      storeLeadQ[ci.li_store_dev_i] = ci.q0;
      storeLeadF[ci.li_store_dev_i] = 0;
    }

    if( ci.loadLeadCurrent )
    {
      storeLeadQ[ci.li_store_dev_i] = ci.q0;
      storeLeadF[ci.li_store_dev_i] = 0;
    }

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      cout << " loading Q vector for cap " << ci.getName() << ":" << endl;
      cout << " q[ " << ci.li_Pos << " ] += " << ci.q0 << endl;
      cout << " q[ " << ci.li_Neg <<  " ] += " << -ci.q0 << endl;

    }
#endif

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

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << " ----------------------------------" << endl;
    cout << " Master::loadDAEMatrices: " << endl;
  }
#endif
  for (InstanceVector::const_iterator it = getInstanceVector().begin(); it != getInstanceVector().end(); ++it)
  {
    Instance & ci = *(*it);

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      cout << " loads for capacitor " << ci.getName() << endl;
    }
#endif

    if (ci.ICGiven && getSolverState().dcopFlag)
    {
      // Special Jacobian if we're doing operating point when IC given
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
      *(ci.fPosEquBraNodePtr) +=  1.0;
      *(ci.fNegEquBraNodePtr) += -1.0;
      *(ci.fBraEquPosNodePtr) +=  1.0;
      *(ci.fBraEquNegNodePtr) += -1.0;
#else
      dFdx[ci.li_Pos][ci.APosEquBraNodeOffset] +=  1.0;
      dFdx[ci.li_Neg][ci.ANegEquBraNodeOffset] += -1.0;
      dFdx[ci.li_Bra][ci.ABraEquPosNodeOffset] +=  1.0;
      dFdx[ci.li_Bra][ci.ABraEquNegNodeOffset] += -1.0;
#endif
    }
    else
    {
      if (ci.ICGiven)
      {
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
        *(ci.fBraEquBraNodePtr) += 1.0;
#else
        dFdx[ci.li_Bra][ci.ABraEquBraNodeOffset] += 1.0;
#endif
      }
    }

    if (!(ci.ICGiven&& getSolverState().dcopFlag))
    {
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
      *(ci.qPosEquPosNodePtr) += ci.C;
      *(ci.qPosEquNegNodePtr) -= ci.C;
      *(ci.qNegEquPosNodePtr) -= ci.C;
      *(ci.qNegEquNegNodePtr) += ci.C;

      if (ci.solVarDepC)
      {
        for (int i=0; i< ci.expNumVars; ++i)
        {
          // Similar to logic in loadDAEdQdX:
          if ((ci.qPosEquDepVarsPtrs[i] != ci.qPosEquPosNodePtr)
              && (ci.qPosEquDepVarsPtrs[i] != ci.qPosEquNegNodePtr))
          {
            *(ci.qPosEquDepVarsPtrs[i]) += ci.expVarDerivs[i];
          }

          if ((ci.qNegEquDepVarsPtrs[i] != ci.qNegEquPosNodePtr)
              && (ci.qNegEquDepVarsPtrs[i] != ci.qNegEquNegNodePtr))
          {
            *(ci.qNegEquDepVarsPtrs[i]) -= ci.expVarDerivs[i];
          }
#ifdef Xyce_DEBUG_DEVICE
          if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
          {
            if ((ci.qPosEquDepVarsPtrs[i] != ci.qPosEquPosNodePtr)
                && (ci.qPosEquDepVarsPtrs[i] != ci.qPosEquNegNodePtr))
            {
              cout << " q[pos][ " << ci.APosEquDepVarOffsets[i] << " ] += "
                   << ci.expVarDerivs[i] << endl;
            }
            if ((ci.qNegEquDepVarsPtrs[i] != ci.qNegEquPosNodePtr)
                && (ci.qNegEquDepVarsPtrs[i] != ci.qNegEquNegNodePtr))
            {
              cout << " q[neg][ " << ci.ANegEquDepVarOffsets[i] << " ] += "
                   << ci.expVarDerivs[i] << endl;
            }
          }
#endif

        }
      }
#else
      dQdx[ci.li_Pos][ci.APosEquPosNodeOffset] += ci.C;
      dQdx[ci.li_Pos][ci.APosEquNegNodeOffset] -= ci.C;
      dQdx[ci.li_Neg][ci.ANegEquPosNodeOffset] -= ci.C;
      dQdx[ci.li_Neg][ci.ANegEquNegNodeOffset] += ci.C;

      if (ci.solVarDepC)
      {
        for (int i=0; i< ci.expNumVars; ++i)
        {
          if ( (ci.APosEquDepVarOffsets[i] != ci.APosEquPosNodeOffset)
               && (ci.APosEquDepVarOffsets[i] != ci.APosEquNegNodeOffset))
          {
            dQdx[ci.li_Pos][ci.APosEquDepVarOffsets[i]] +=
              ci.expVarDerivs[i];
          }
          if ( (ci.ANegEquDepVarOffsets[i] != ci.ANegEquPosNodeOffset)
               && (ci.ANegEquDepVarOffsets[i] != ci.ANegEquNegNodeOffset))
          {
            dQdx[ci.li_Neg][ci.ANegEquDepVarOffsets[i]] -=
              ci.expVarDerivs[i];
          }
#ifdef Xyce_DEBUG_DEVICE
          if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
          {
            if ( (ci.APosEquDepVarOffsets[i] != ci.APosEquPosNodeOffset)
                 && (ci.APosEquDepVarOffsets[i] != ci.APosEquNegNodeOffset))
            {
              cout << " q[pos][ " << ci.APosEquDepVarOffsets[i] << " ] += "
                   << expVarDerivs[i] << endl;
            }
            if ( (ci.ANegEquDepVarOffsets[i] != ci.ANegEquPosNodeOffset)
                 && (ci.ANegEquDepVarOffsets[i] != ci.ANegEquNegNodeOffset))
            {
              cout << " q[neg][ " << ci.ANegEquDepVarOffsets[i] << " ] += "
                   << ci.expVarDerivs[i] << endl;
            }
          }
#endif
        }
      }

#endif
    }
  }

  return true;
}

} // namespace Capacitor
} // namespace Device
} // namespace Xyce

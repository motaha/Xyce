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
// Filename       : $RCSfile: N_DEV_Digital.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Dave Shirley, PSSI
//
// Creation Date  : 01/05/06
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.69.2.2 $
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
#include <iostream>

#include <sstream>

// ----------   Xyce Includes   ----------
#include <N_DEV_Const.h>
#include <N_DEV_Digital.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_MatrixLoadData.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>

namespace Xyce {
namespace Device {

template<>
ParametricData<Digital::Instance>::ParametricData()
{
  setNumNodes(2);
  setNumOptionalNodes(20);
  setNumFillNodes(0);
  setModelRequired(1);
  setPrimaryParameter("");
  addModelType("DIG");

  // Set up double precision variables:

  // Set up non-double precision variables:
  addPar ("IC1", false, false, ParameterType::NO_DEP,
          &Digital::Instance::ic1, NULL,
          U_LOGIC, CAT_NONE, "Vector of initial values for output(s)");
  addPar ("IC2", false, false, ParameterType::NO_DEP,
          &Digital::Instance::ic2, NULL,
          U_NONE, CAT_NONE, "");
  addPar ("IC3", false, false, ParameterType::NO_DEP,
          &Digital::Instance::ic3, NULL,
          U_NONE, CAT_NONE, "");
  makeVector ("IC", 3);
}

template<>
ParametricData<Digital::Model>::ParametricData()
{
  // Set up double precision variables:
  addPar ("VLO", 0., false, ParameterType::NO_DEP,
          &Digital::Model::vlo,
          NULL, U_VOLT, CAT_NONE, "Internal low state supply voltage");

  addPar ("VHI", 0., false, ParameterType::NO_DEP,
          &Digital::Model::vhi,
          NULL, U_VOLT, CAT_NONE, "Internal high state supply voltage");

  addPar ("VREF", 0., false, ParameterType::NO_DEP,
          &Digital::Model::vref,
          NULL, U_VOLT, CAT_NONE, "Internal reference voltage for inputs");

  addPar ("CLO", 1.e-6, false, ParameterType::NO_DEP,
          &Digital::Model::clo,
          NULL, U_FARAD, CAT_NONE, "Capacitance between output node and low reference");

  addPar ("CHI", 1.e-6, false, ParameterType::NO_DEP,
          &Digital::Model::chi,
          NULL, U_FARAD, CAT_NONE, "Capacitance between output node and high reference");

  addPar ("CLOAD", 1.e-6, false, ParameterType::NO_DEP,
          &Digital::Model::cload,
          NULL, U_FARAD, CAT_NONE, "Capacitance between input node and input reference");

  addPar ("RLOAD", 1000., false, ParameterType::NO_DEP,
          &Digital::Model::rload,
          NULL, U_OHM, CAT_NONE, "Resistance between input node and input reference");

  addPar ("S0RLO", 100., false, ParameterType::NO_DEP,
          &Digital::Model::s0rlo,
          NULL, U_OHM, CAT_NONE, "Low state resistance between output node and low reference");

  addPar ("S0RHI", 100., false, ParameterType::NO_DEP,
          &Digital::Model::s0rhi,
          NULL, U_OHM, CAT_NONE, "Low state resitance between output node and high reference");

  addPar ("S0TSW", 1.e-8, false, ParameterType::NO_DEP,
          &Digital::Model::s0tsw,
          NULL, U_SECOND, CAT_NONE, "Switching time transition to low state");

  addPar ("S0VLO", -1.5, false, ParameterType::NO_DEP,
          &Digital::Model::s0vlo,
          NULL, U_VOLT, CAT_NONE, "Minimum voltage to switch to low state");

  addPar ("S0VHI", 1.7, false, ParameterType::NO_DEP,
          &Digital::Model::s0vhi,
          NULL, U_VOLT, CAT_NONE, "Maximum voltage to switch to low state");

  addPar ("S1RLO", 100., false, ParameterType::NO_DEP,
          &Digital::Model::s1rlo,
          NULL, U_OHM, CAT_NONE, "High state resistance between output node and low reference");

  addPar ("S1RHI", 100., false, ParameterType::NO_DEP,
          &Digital::Model::s1rhi,
          NULL, U_OHM, CAT_NONE, "High state resistance between output node and high reference");

  addPar ("S1TSW", 1.e-8, false, ParameterType::NO_DEP,
          &Digital::Model::s1tsw,
          NULL, U_SECOND, CAT_NONE, "Switching time transition to high state");

  addPar ("S1VLO", 0.9, false, ParameterType::NO_DEP,
          &Digital::Model::s1vlo,
          NULL, U_VOLT, CAT_NONE, "Minimum voltage to switch to high state");

  addPar ("S1VHI", 7.0, false, ParameterType::NO_DEP,
          &Digital::Model::s1vhi,
          NULL, U_VOLT, CAT_NONE, "Maximum voltage to switch to high state");

  addPar ("DELAY", 1.e-8, false, ParameterType::NO_DEP,
          &Digital::Model::delay,
          NULL, U_SECOND, CAT_NONE, "Delay time of device");
}

namespace Digital {



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
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/05/06
//-----------------------------------------------------------------------------
bool Instance::processParams (string param)
{

  // If there are any time dependent parameters, set their values for
  // the current time.

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : instance block constructor
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/05/06
//-----------------------------------------------------------------------------
Instance::Instance(
  InstanceBlock & IB,
  Model & Diter,
  MatrixLoadData & mlData1,
  SolverState &ss1,
  ExternData  &ed1,
  DeviceOptions & do1)
  : DeviceInstance(IB, mlData1, ss1, ed1, do1),
    model_(Diter),
    li_Lo(-1),
    li_Hi(-1),
    li_Ref(-1),
    row_Lo(-1),
    row_Hi(-1),
    row_Ref(-1),
    breakTime(0.)
{
  // can't check for MPDE at this point because MPD_Manager may not exist yet.
  // for example with Xyce -param case.  So for now we can't catch this
  // error here.

  // #ifdef Xyce_MPDE
  //  string msg("::addInstance digital devices not compatible with MPDE");
  //  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL, std::ostringstream() << "Error in " << netlistLocation() << "\n" << msg);
  // #endif // Xyce_MPDE

  setName(IB.getName());
  int i, p1, p2;
  p1 = getName().find_first_of('%');
  p2 = getName().find_last_of('%');
  if (p1 != string::npos && p2 != string::npos && p2 > p1)
  {
    numExtVars = 0;
    if (!model_.given("VLO"))
    {
      li_Lo = 0;
      ++numExtVars;
    }
    if (!model_.given("VHI"))
    {
      li_Hi = 0;
      ++numExtVars;
    }
    if (!model_.given("VREF"))
    {
      li_Ref = 0;
      ++numExtVars;
    }
    string dev_type(getName().substr(p1+1,p2-p1-1));
    if (dev_type == "NOT")
    {
      numInput = 1;
      numOutput = 1;
      gate = NOT;
    }
    else if (dev_type == "AND")
    {
      numInput = 2;
      numOutput = 1;
      gate = AND;
    }
    else if (dev_type == "NAND")
    {
      numInput = 2;
      numOutput = 1;
      gate = NAND;
    }
    else if (dev_type == "OR")
    {
      numInput = 2;
      numOutput = 1;
      gate = OR;
    }
    else if (dev_type == "NOR")
    {
      numInput = 2;
      numOutput = 1;
      gate = NOR;
    }
    else if (dev_type == "ADD")
    {
      numInput = 3;
      numOutput = 2;
      gate = ADD;
    }
    else if (dev_type == "XOR")
    {
      numInput = 2;
      numOutput = 1;
      gate = XOR;
    }
    else if (dev_type == "NXOR")
    {
      numInput = 2;
      numOutput = 1;
      gate = NXOR;
    }
    //Genie 110812
    else if (dev_type == "DFF")
    {
      numInput = 4;  //PREB, CLRB, clock, data
      numOutput = 2; //Q, Q_bar
      gate = DFF;
    }
    else
    {
      string msg("Unknown digital device type: ");
      msg += dev_type;
      std::ostringstream oss;
      oss << "Error in " << netlistLocation() << "\n" << msg;
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL, oss.str());
    }
  }
  else if (getName().substr(0,7) == "Digital")
  {
    numExtVars = 0;
    numInput = 0;
    numOutput = 0;
  }
  else
  {
    string msg("Internal error in digital device name: ");
    msg += getName();
    std::ostringstream oss;
    oss << "Error in " << netlistLocation() << "\n" << msg;
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, oss.str());
  }

  int iBase = numExtVars;
  int oBase = numExtVars + numInput;
  numExtVars += numInput + numOutput;

////cout << "iBase = " << iBase << ", oBase = " << oBase << ", numExtVars = " << numExtVars << endl;


  if (numExtVars != IB.numExtVars)
  {
    ostringstream msg;
    msg << "Incorrect number of nodes in digital device: ";
    msg << getName();
    msg << ".  Found: ";
    msg << IB.numExtVars;
    msg << ", Should be: ";
    msg << numExtVars;
    msg << ".";
    std::ostringstream oss;
    oss << "Error in " << netlistLocation() << "\n" << msg.str();
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL, oss.str());
  }
  numIntVars   = 0;
  numStateVars = 4*numInput + 6*numOutput;

  li_Inp.resize(numInput);
  li_currentStateInp.resize(numInput);
  li_transitionTimeInp.resize(numInput);
  li_QinpState.resize(numInput);
  li_IinpState.resize(numInput);
  li_Out.resize(numOutput);
  li_currentStateOut.resize(numOutput);
  li_transitionTimeOut.resize(numOutput);
  li_QloState.resize(numOutput);
  li_IloState.resize(numOutput);
  li_QhiState.resize(numOutput);
  li_IhiState.resize(numOutput);

  qlo.resize(numOutput);
  ilo.resize(numOutput);
  vcaplo.resize(numOutput);
  qhi.resize(numOutput);
  ihi.resize(numOutput);
  vcaphi.resize(numOutput);
  qref.resize(numInput);
  iref.resize(numInput);
  vcapref.resize(numInput);
  rilo.resize(numOutput);
  rihi.resize(numOutput);
  riref.resize(numInput);
  currentOut.resize(numOutput);
  currentIn.resize(numInput);
  glo.resize(numOutput);
  ghi.resize(numOutput);

  qInp.resize(numInput);
  iInp.resize(numInput);
  vcapInp.resize(numInput);
  currentInp.resize(numInput);

  inpL.resize(numInput);
  iTime.resize(numInput);
  outL.resize(numOutput);
  oTime.resize(numOutput);

  //Genie 110812.
  //changeState.resize(numInput);


  // These are to store the indicies into the jacobian for the four element
  // stamps for the capacitor/resistors connected to the input/outputs. The
  // format is to have 6 int vectors for each stamp with format:
  // (row 1, col 1, col 2, row 2, col 1, col 2).  These will be replaced in
  // registerJacLIDs with the indicies into the actual jacobian matrix.

  li_jac_Ref.resize(numInput);
  li_jac_Hi.resize(numOutput);
  li_jac_Lo.resize(numOutput);

  devConMap.resize(numExtVars);
  for (i=0 ; i<numExtVars ; ++i)
    devConMap[i] = 1;

  setModelName(model_.getName());

  jacStamp.resize(numExtVars);
  int row = 0;
  if (li_Lo == 0)
  {
    row_Lo = row;
    jacStamp[row].push_back(row);
    for (i=0 ; i<numOutput ; ++i)
    {
      li_jac_Lo[i].push_back(row);
      li_jac_Lo[i].push_back(0);
      li_jac_Lo[i].push_back(jacStamp[row].size());
      li_jac_Lo[i].push_back(oBase+i);
      li_jac_Lo[i].push_back(jacStamp[oBase+i].size());
      jacStamp[row].push_back(oBase+i);
      jacStamp[oBase+i].push_back(row);
    }
    ++row;
  }
  if (li_Hi == 0)
  {
    row_Hi = row;
    jacStamp[row].push_back(row);
    for (i=0 ; i<numOutput ; ++i)
    {
      li_jac_Hi[i].push_back(row);
      li_jac_Hi[i].push_back(0);
      li_jac_Hi[i].push_back(jacStamp[row].size());
      li_jac_Hi[i].push_back(oBase+i);
      li_jac_Hi[i].push_back(jacStamp[oBase+i].size());
      jacStamp[row].push_back(oBase+i);
      jacStamp[oBase+i].push_back(row);
    }
    ++row;
  }
  if (li_Ref == 0)
  {
    row_Ref = row;
    jacStamp[row].push_back(row);
    for (i=0 ; i<numInput ; ++i)
    {
      li_jac_Ref[i].push_back(row);
      li_jac_Ref[i].push_back(0);
      li_jac_Ref[i].push_back(jacStamp[row].size());
      li_jac_Ref[i].push_back(iBase+i);
      li_jac_Ref[i].push_back(jacStamp[iBase+i].size());
      jacStamp[row].push_back(iBase+i);
      jacStamp[iBase+i].push_back(row);
    }
    ++row;
  }
  for (i=0 ; i<numInput ; ++i)
  {
    if (li_Ref == 0)
    {
      li_jac_Ref[i].push_back(jacStamp[iBase+i].size());
    }
    else
    {
      li_jac_Ref[i].push_back(iBase+i);
    }
    jacStamp[iBase+i].push_back(iBase+i);
  }
  for (i=0 ; i<numOutput ; ++i)
  {
    if (li_Lo == 0)
    {
      li_jac_Lo[i].push_back(jacStamp[oBase+i].size());
    }
    else
    {
      li_jac_Lo[i].push_back(oBase+i);
    }
    if (li_Hi == 0)
    {
      li_jac_Hi[i].push_back(jacStamp[oBase+i].size());
    }
    else
    {
      li_jac_Hi[i].push_back(oBase+i);
    }
    jacStamp[oBase+i].push_back(oBase+i);
  }


  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  //  Set any non-constant parameter defaults:

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params:
  processParams ();
}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/05/06
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
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/05/06
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

  int i=0, j;
  if (li_Lo == 0)
    li_Lo = extLIDVec[i++];
  if (li_Hi == 0)
    li_Hi = extLIDVec[i++];
  if (li_Ref == 0)
    li_Ref = extLIDVec[i++];
  for (j=0 ; j<numInput ; ++j)
    li_Inp[j] = extLIDVec[i++];
  for (j=0 ; j<numOutput ; ++j)
    li_Out[j] = extLIDVec[i++];

}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/05/06
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs( const vector<int> & staLIDVecRef)
{
  string msg;

  // Check if the size of the ID lists corresponds to the proper number of
  // internal and external variables.
  int numSta = staLIDVecRef.size();

  if (numSta != numStateVars)
  {
    msg = "Instance::registerStateLIDs:";
    msg += "numSta != numStateVars";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
  }

  // Copy over the global ID lists:
  staLIDVec = staLIDVecRef;

  int i, j=0;

  for (i=0 ; i<numInput ; ++i)
  {
    li_currentStateInp[i] = staLIDVec[j++];
    li_transitionTimeInp[i] = staLIDVec[j++];
    li_QinpState[i] = staLIDVec[j++];
    li_IinpState[i] = staLIDVec[j++];
  }

  for (i=0 ; i<numOutput ; ++i)
  {
    li_currentStateOut[i] = staLIDVec[j++];
    li_transitionTimeOut[i] = staLIDVec[j++];
    li_QloState[i] = staLIDVec[j++];
    li_IloState[i] = staLIDVec[j++];
    li_QhiState[i] = staLIDVec[j++];
    li_IhiState[i] = staLIDVec[j++];
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::getIntNameMap
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/05/06
//-----------------------------------------------------------------------------
map<int,string> & Instance::getIntNameMap ()
{
  return intNameMap;
}


//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/05/06
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
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/05/06
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const vector< vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );
  int i;
  int iBase = numExtVars - numInput - numOutput;
  int oBase = iBase + numInput;
  int lo_present, hi_present;

  if (row_Lo < 0)
    lo_present = 0;
  else
    lo_present = 1;
  if (row_Hi < 0)
    hi_present = 0;
  else
    hi_present = 1;

  if (row_Ref == -1)
  {
    for (i=0 ; i<numInput ; ++i)
    {
      li_jac_Ref[i].push_back(jacLIDVec[li_jac_Ref[i][0]][0]);
    }
  }
  else
  {
    for (i=0 ; i<numInput ; ++i)
    {
      li_jac_Ref[i][1] = jacLIDVec[li_jac_Ref[i][0]][li_jac_Ref[i][1]];
      li_jac_Ref[i][2] = jacLIDVec[li_jac_Ref[i][0]][li_jac_Ref[i][2]];
      li_jac_Ref[i][4] = jacLIDVec[li_jac_Ref[i][3]][li_jac_Ref[i][4]];
      li_jac_Ref[i][5] = jacLIDVec[li_jac_Ref[i][3]][li_jac_Ref[i][5]];
    }
  }

  if (row_Lo == -1)
  {
    for (i=0 ; i<numOutput ; ++i)
    {
      li_jac_Lo[i].push_back(jacLIDVec[li_jac_Lo[i][0]][hi_present]);
    }
  }
  else
  {
    for (i=0 ; i<numOutput ; ++i)
    {
      li_jac_Lo[i][1] = jacLIDVec[li_jac_Lo[i][0]][li_jac_Lo[i][1]];
      li_jac_Lo[i][2] = jacLIDVec[li_jac_Lo[i][0]][li_jac_Lo[i][2]];
      li_jac_Lo[i][4] = jacLIDVec[li_jac_Lo[i][3]][li_jac_Lo[i][4]];
      li_jac_Lo[i][5] = jacLIDVec[li_jac_Lo[i][3]][li_jac_Lo[i][5]];
    }
  }

  if (row_Hi == -1)
  {
    for (i=0 ; i<numOutput ; ++i)
    {
      li_jac_Hi[i].push_back(jacLIDVec[li_jac_Hi[i][0]][lo_present]);
    }
  }
  else
  {
    for (i=0 ; i<numOutput ; ++i)
    {
      li_jac_Hi[i][1] = jacLIDVec[li_jac_Hi[i][0]][li_jac_Hi[i][1]];
      li_jac_Hi[i][2] = jacLIDVec[li_jac_Hi[i][0]][li_jac_Hi[i][2]];
      li_jac_Hi[i][4] = jacLIDVec[li_jac_Hi[i][3]][li_jac_Hi[i][4]];
      li_jac_Hi[i][5] = jacLIDVec[li_jac_Hi[i][3]][li_jac_Hi[i][5]];
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/05/06
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState ()
{
  bool bsuccess = true;
  double v_poslo, v_poshi, v_posref, v_neg;
  double elapsed, time, frac, lastT;
  int currentState;
  double transitionTime;
  bool changeState = false; //Genie 110812
  bool clocking = false; //Genie 111212
  //vector<bool> changeState;  //Genie 110812
  bool toPrint = false; //Genie 022713

  int i;

  N_LAS_Vector & solVector = *(extData.nextSolVectorPtr);
  N_LAS_Vector & staVector = *(extData.nextStaVectorPtr);
  N_LAS_Vector & oldStaVector = *(extData.currStaVectorPtr);
  N_LAS_Vector & oldSolVector = *(extData.currSolVectorPtr);

  //Genie 113012
  ////cout <<"Digital::updatePrimaryState() is called." <<endl;

  // The convention in this device is to consider the 'positive' side of the
  // capacitors as the Vsrc or node that supplies the voltage, or for the output,
  // vref.  The 'positive' voltages are thus the same for all input/outputs
  if (li_Lo >= 0)
    v_poslo = solVector[li_Lo];
  else
    v_poslo = model_.vlo;
  if (li_Hi >= 0)               //Genie 112712. If VHI is not speficied, li_Hi >=0 and
    v_poshi = solVector[li_Hi];	// v_poshi is derived from CLOAD/RLOAD by solVector.
  else				// If VHI is specified, li_Hi == -1 and v_poshi is specified
    v_poshi = model_.vhi;	// by the netlist model card.
  if (li_Ref >= 0)
    v_posref = solVector[li_Ref];
  else
    v_posref = model_.vref;

  lastT = 0;

  //Genie 110812
  //changeState.resize(numInput);
  //for (i=0 ; i<numInput ; ++i)
  //     changeState[i] = false;

  //Genie 022713
  if (getSolverState().currTime > 3)
     toPrint = true;

  for (i=0 ; i<numInput ; ++i)
  {
    //initialize
    currentState = static_cast <int> (oldStaVector[li_currentStateInp[i]]);
    transitionTime = oldStaVector[li_transitionTimeInp[i]];
    //changeState[i] = false; //Genie 110812
    changeState = false; //Genie 022013 Clear the memory of changeState.
			//Attempt to debug DFF state not changing issue due to wrong oTime


    // obtain voltage drop accross the capacitor:

    v_neg = solVector[li_Inp[i]];

    time = getSolverState().currTime;
    elapsed = time - transitionTime;

    vcapref[i] = v_posref-v_neg;
    riref[i] = model_.gload*(vcapref[i]);

    // Obtain the "current"  value for the charge stored in the capacitors.
    qref[i] = model_.cload*vcapref[i];

    if (getSolverState().dcopFlag)
    {
      transitionTime = 0;
      if (-vcapref[i] < model_.s0vhi)
        currentState = 0;
      else
        currentState = 1;

      oldStaVector[li_currentStateInp[i]] = currentState;
      oldStaVector[li_transitionTimeInp[i]] = transitionTime;
    }

    //Genie 120312 Debug
    /*if (gate == NOT)
  {
   cout << "UpdatePrimaryState(NOT): dcopFlag = " << getSolverState().dcopFlag << ", s0vhi = " << (*M_iter)->s0vhi << "currentState = " << currentState << ", s1vlo = " << (*M_iter)->s1vlo << endl;
   cout << ", v_posref = " << v_posref << "vcapref[0] = " << vcapref[i] << ", v_neg = " << v_neg << endl;
  }*/
    //Genie 121812 clock Debug
    /*if (gate == DFF && i == 2)
  {
   cout << "UpdatePrimaryState(DFF): dcopFlag = " << getSolverState().dcopFlag << ", s0vhi = " << (*M_iter)->s0vhi << "currentState(in2) = " << currentState << ", s1vlo = " << (*M_iter)->s1vlo << endl;
   cout << ", vcapref[2] = " << vcapref[2] << endl;
   cout << ", oldStaVector[li_currentStateInp[2]] = " << static_cast <int> (oldStaVector[li_currentStateInp[2]]) << endl;
  }*/
   //Genie 021913 Debug
   /*
   if (gate == DFF && i == 3 && toPrint) //data debug
  {
   cout << "At simulation time " << getSolverState().currTime << endl;
   cout << "UpdatePrimaryState(DFF): dcopFlag = " << getSolverState().dcopFlag << ", s0vhi = " << (*M_iter)->s0vhi << "currentState(in3) = " << currentState << ", s1vlo = " << (*M_iter)->s1vlo << endl;
   cout << ", vcapref[3] = " << vcapref[3] << endl;
   cout << ", oldStaVector[li_currentStateInp[3]] = " << static_cast <int> (oldStaVector[li_currentStateInp[3]]) << endl;
  }*/
   //Genie 022813 Vrf debug
   if (gate == DFF && i == 1 && toPrint) //data debug
  {
   cout << "At simulation time " << getSolverState().currTime << endl;
   cout << "UpdatePrimaryState(DFF): dcopFlag = " << getSolverState().dcopFlag << ", s0vhi = " << model_.s0vhi << "currentState(in1) = " << currentState << ", s1vlo = " << model_.s1vlo << endl;
   cout << ", vcapref[1] = " << vcapref[1] << endl;
   cout << ", oldStaVector[li_currentStateInp[1]] = " << static_cast <int> (oldStaVector[li_currentStateInp[1]]) << endl;
  }


    iTime[i] = transitionTime;

    staVector[li_transitionTimeInp[i]] = transitionTime;

    if (currentState == 0)
    {
      inpL[i] = false;
      if (-vcapref[i] > model_.s0vhi)
      {
        if (-vcapref[i] > model_.s1vlo)
        {
          currentState = 1;
          changeState = true;
	  ////cout << "Gate Type " << gate << " changes state to 1" << endl;
	  if (gate == DFF && i == 2)  //Genie 111212: i==2 -> clk
          {   clocking = true;  // clock of  DFF changes state
            ////cout << "clock of DFF changes state to 1" << endl;
	  }
        }
      }
    }
    else
    {
      inpL[i] = true;
      if (-vcapref[i] < model_.s1vlo)
      {
        if (-vcapref[i] < model_.s0vhi)
        {
          currentState = 0;
          changeState = true;
	  ////cout << "Gate type "<< gate << " changes state to 0" << endl;
	  if (gate == DFF && i == 2)  //Genie 111212
	  {
            clocking = true;  // clock of  DFF changes state
            ////cout << "clock of DFF changes state to 0" << endl;
	  }
        }
      }
    }

    if (changeState)
    {
      double vOld, del;

      inpL[i] = (currentState == 1);
      if (li_Ref >= 0)
        vOld = oldSolVector[li_Ref];
      else
        vOld = model_.vref;
      vOld = oldSolVector[li_Inp[i]] - vOld;
      if (fabs(vcapref[i]+vOld) < 1.e-12)
        del = 0;
      else
      {
        if (inpL[i])
        {
          del = getSolverState().currTimeStep * (-vcapref[i] - model_.s0vhi)/(-vcapref[i] - vOld);
        }
        else
        {
          del = getSolverState().currTimeStep * (vcapref[i] + model_.s1vlo)/(vcapref[i] + vOld);
        }
      }
      staVector[li_transitionTimeInp[i]] = time - del;
      iTime[i] = time;
    }

    staVector[li_currentStateInp[i]] = currentState;
    staVector[li_QinpState[i]] = qref[i];

    if (iTime[i] > lastT)
      lastT = iTime[i];
  } // end for loop numInput i

  if (gate == NOT)
  {
    outL[0] = !inpL[0];
    //Genie 113012
    ////cout << "NOT: inpL[0] = " << inpL[0] << ", outL[0] = " << outL[0] << endl;
    oTime[0] = lastT+model_.delay;
  }
  else if (gate == AND)
  {
    outL[0] = inpL[0] & inpL[1];
    oTime[0] = lastT+model_.delay;
  }
  else if (gate == NAND)
  {
    outL[0] = !(inpL[0] & inpL[1]);
    oTime[0] = lastT+model_.delay;
  }
  else if (gate == OR)
  {
    outL[0] = inpL[0] | inpL[1];
    oTime[0] = lastT+model_.delay;
  }
  else if (gate == NOR)
  {
    outL[0] = !(inpL[0] | inpL[1]);
    oTime[0] = lastT+model_.delay;
  }
  else if (gate == ADD)
  {
    outL[0] = inpL[0] ^ inpL[1] ^ inpL[2];
    //outL[0] = (inpL[0] & inpL[1]) | (inpL[1] & inpL[2]) | (inpL[0] & inpL[2]);
    outL[1] = (inpL[0] & inpL[1]) | (inpL[1] & inpL[2]) | (inpL[0] & inpL[2]); //Genie, 101612. carry-out sum

    oTime[0] = lastT+model_.delay;
    oTime[1] = lastT+model_.delay;
  }
  else if (gate == XOR)
  {
    outL[0] = inpL[0] ^ inpL[1];
    oTime[0] = lastT+model_.delay;
  }
  else if (gate == NXOR)
  {
    outL[0] = !(inpL[0] ^ inpL[1]);
    oTime[0] = lastT+model_.delay;
  }
  //Genie
  else if (gate == DFF)
  { // DFF: in0: PREB, in1: CLRB, in2: clock, in3: data
    // DFF: out0: Q, out1: Q_bar
    // CD4013B: set = !PREB, reset = !CLRB
    // CD4013B: in0: set, in1: reset, in2: clock, in3: data
    // CD4013B: out0: Q, out1: Q_bar
    if (clocking && inpL[2] ==1) //clock rising edge 0->1
    {
        if (toPrint){
          cout << "At simulation time: " << getName() << ", " << getSolverState().currTime << endl;
	   cout << "clock rising edge, D = " << inpL[3] << ", set = " << !inpL[0] << ", reset = " << !inpL[1] << endl;
	}
        if (inpL[0] == 1 && inpL[1] == 1) //PREB = CLRB = 1
 	//if (inpL[0] == 0 && inpL[1] == 0) //S = R = 0
        //if (inpL[0] == 1 ) //Genie 022813. Brute force debug for...PREB = CLRB = 1
        {
            outL[0] = inpL[3];  //Q = D
            outL[1] = !(inpL[3]); //Q_bar = !D
	    if (toPrint) {
	        cout << "Set=Reset=0, " << "D = " << inpL[3] << ", Q = " << outL[0] << ", Q_bar = " << outL[1] << endl;
	    }
	}
    }
    else if (clocking && inpL[2] ==0) //clock falling edge 1->0
    {
        if (toPrint)
	     cout << "clock falling edge, D= " << inpL[3] << ", set = " << !inpL[0] << ", reset = " << !inpL[1] << endl;
 	if (inpL[0] == 1 && inpL[1] == 1) //PREB = CLRB = 1
	//if (inpL[0] == 0 && inpL[1] == 0) //S = R = 0
        {
            outL[0] = oldStaVector[li_currentStateOut[0]]; //no change
            outL[1] = oldStaVector[li_currentStateOut[1]]; //no change
             if (toPrint) {
	        cout << "Set=Reset=0, " << "D = " << inpL[3] << ", no change Q = " << outL[0] << ", Q_bar = " << outL[1] << endl;
	    }

	}
    }
    else // no clock change
    {
      if (inpL[0] == 1 && inpL[1] == 0) //PREB=1, CLRB = 0
	//if (inpL[0] == 0 && inpL[1] == 1) //S=0, R = 1
        {
            outL[0] = 0;
            outL[1] = 1;
	     if (toPrint) {
	        cout << "Set=0, Reset=1, " << "D = " << inpL[3] << ", Q = " << outL[0] << ", Q_bar = " << outL[1] << endl;
	    }
	}
	else if (inpL[0] == 0 && inpL[1] == 1) //PREB = 0, CLRB = 1
	//else if (inpL[0] == 1 && inpL[1] == 0) //S = 1, R = 0
        {
            outL[0] = 1;
            outL[1] = 0;
	     if (toPrint) {
	        cout << "Set=1, Reset=0, " << "D = " << inpL[3] << ", Q = " << outL[0] << ", Q_bar = " << outL[1] << endl;
	    }

	}
	else if (inpL[0] == 0 && inpL[1] == 0) //PREB = CLRB = 0
	//else if (inpL[0] == 1 && inpL[1] == 1) //S = R = 1
        {
            outL[0] = 1;
            outL[1] = 1;
	     if (toPrint) {
	        cout << "Set=Reset=1, " << "D = " << inpL[3] << ", Q = " << outL[0] << ", Q_bar = " << outL[1] << endl;
	    }
	}
    }
    oTime[0] = lastT+model_.delay;
    oTime[1] = lastT+model_.delay;
  }

  bool curr;
  breakTime = 0;
  for (i=0 ; i<numOutput ; ++i)
  {
    time = getSolverState().currTime;
    if (getSolverState().dcopFlag)
    {
      if (i == 0)
      {
        if (given("IC1"))
          outL[i] = ic1;
      }
      else if (i == 1)
      {
        if (given("IC2"))
          outL[i] = ic2;
      }
      else if (i == 2)
      {
        if (given("IC3"))
          outL[i] = ic3;
      }
      else
      {
        string msg("Insufficient initial conditions supported in digital device");
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
      }

      oldStaVector[li_currentStateOut[i]] = outL[i]?1:0;
      oldStaVector[li_transitionTimeOut[i]] = time;

      //Genie 022013 debug
      /*if (gate == DFF){
          cout << "DFF Output: At simulation time " << getSolverState().currTime << endl;
          cout << "oldStaVector[li_currentStateOut[" << i << "] = " << oldStaVector[li_currentStateOut[i]] << endl;
      }*/
    }

    //current logic state of output nodes
    currentState = static_cast <int> (oldStaVector[li_currentStateOut[i]]);
    transitionTime = oldStaVector[li_transitionTimeOut[i]];

    /*if (gate == NOT) //Genie
      cout << "Digital(NOT):: currentState (out) = " << currentState << ", outL[i] = " << outL[i] << endl; */
    //Genie
    if (gate == DFF && toPrint)
      cout << "Digital(DFF):: " << getName() << ", currentState (out) = " << currentState << ", i = " << i << ", outL[0] = " << outL[0] << ", outL[1] = " << outL[1] << ", oTime[i] = " << oTime[i] << endl;

    if (currentState == 1)
      curr = true;
    else
      curr = false;

    if (curr != outL[i]) //Genie 110812. This is executed when scopFlag is false
    {
      if (oTime[i] <= time)
      {
        currentState = 1-currentState;
        transitionTime = oTime[i];
      }
      else
      {
        if (breakTime == 0 || (breakTime > 0 && breakTime > oTime[i]))
        {
          breakTime = oTime[i];
        }
      }
    }

    staVector[li_currentStateOut[i]] = currentState;
    staVector[li_transitionTimeOut[i]] = transitionTime;

    // obtain voltage drop accross the capacitors:

    v_neg = solVector[li_Out[i]];

    elapsed = time - transitionTime;

    if (currentState == 0)
      elapsed /= model_.s0tsw;
    else if (currentState == 1)
      elapsed /= model_.s1tsw;
    else
    {
      string msg("Instance::updateSecondaryState: unrecognized state");
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
    }

    frac = exp(-elapsed); //Genie 112812. This line can be omitted.
    if (transitionTime == 0)
      frac = 0;
    else
    {
      if (elapsed > 1)
        frac = 0;
      else
      {
        // This is a simple linear transition.  Since there is a breakpoint at the start of
        // the transition it is OK to have a discontinuity there.
        frac = 1-elapsed;
      }
    }

    if (currentState == 0)
    {
      glo[i] = 1/(frac*model_.s1rlo + (1-frac)*model_.s0rlo);
      ghi[i] = 1/(frac*model_.s1rhi + (1-frac)*model_.s0rhi);
    }
    else
    {
      glo[i] = 1/(frac*model_.s0rlo + (1-frac)*model_.s1rlo);
      ghi[i] = 1/(frac*model_.s0rhi + (1-frac)*model_.s1rhi);
    }

    //Genie 113012
    /*if (gate == NOT)
  {
    cout << "Gate NOT new currentState (out)= " << currentState << ", frac = " << frac << endl;
    cout << "s0rlo = " << model_.s0rlo << ", s1rlo = " << model_.s1rlo << endl;
  }*/
    //Genie 021913
   if (gate == DFF && toPrint)
  {
    cout << "At simulation time = " << getSolverState().currTime << endl;
    cout << "Gate DFF new currentState (out[" << i << "]= " << currentState << endl;
  }

    rilo[i] = glo[i]*(v_poslo-v_neg);
    rihi[i] = ghi[i]*(v_poshi-v_neg);

   /*if (gate == NOT)
  {
   cout << "rilo[" << i << "] = " << rilo[i] << ", glo[i] = " << glo[i];
   cout << "v_poslo = " << v_poslo << "v_poshi = " << v_poshi << ", v_neg = " << v_neg << endl;
  }*/
   /*if (gate == DFF)
  {
   cout << "v_poslo = " << v_poslo << "v_poshi = " << v_poshi << ", v_neg = " << v_neg << endl;
  }*/
    vcaplo[i] = v_poslo-v_neg;
    vcaphi[i] = v_poshi-v_neg;

    // Obtain the "current"  value for the charge stored in the capacitors.
    ////qlo[i] = model_.clo*vcaphi[i]; //Genie 022113 Found bug?
    qlo[i] = model_.clo*vcaplo[i];  //Genie 022113. Fix
    qhi[i] = model_.chi*vcaphi[i];

   /*if (gate == NOT)
   {
    cout << "vcaplo = " << vcaplo[i] << ", vcaphi = " << vcaphi[i] << endl;
    cout << "qhi[i] = " << qhi[i] << ", chi = " << model_.chi << endl;
    cout << "qlo[i] = " << qlo[i] << ", clo = " << model_.clo << endl;
    cout << "model_.chi*vcaphi[i] = " << model_.chi*vcaphi[i] << endl;
   }*/
   if (gate == DFF && toPrint)
   {
    cout << "vcaplo[" << i << "] = "<< vcaplo[i] << ", vcaphi[" << i << "] = " << vcaphi[i] << endl;
   }

    staVector[li_QloState[i]] = qlo[i];
    staVector[li_QhiState[i]] = qhi[i];
  }// end for loop numOutput i

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/05/06
//-----------------------------------------------------------------------------
bool Instance::updateSecondaryState ()
{
  bool bsuccess = true;
  int i;

  N_LAS_Vector * staVectorPtr = extData.nextStaVectorPtr;

  // Now that the state vector for time=0 is up-to-date, get the
  // derivative with respect to time of the charge, to obtain the
  // best estimate for the current in the capacitors.

  for (i=0 ; i<numOutput ; ++i)
  {
    ilo[i] = (*extData.nextStaDerivVectorPtr)[li_QloState[i]];
    ihi[i] = (*extData.nextStaDerivVectorPtr)[li_QhiState[i]];

    currentOut[i] = ilo[i] + ihi[i] + rilo[i] + rihi[i];

    (*staVectorPtr)[li_IloState[i]] = ilo[i];
    (*staVectorPtr)[li_IhiState[i]] = ihi[i];
  }

  for (i=0 ; i<numInput ; ++i)
  {
    iref[i] = (*extData.nextStaDerivVectorPtr)[li_QinpState[i]];

    currentIn[i] = iref[i] + riref[i];

    (*staVectorPtr)[li_IinpState[i]] = iref[i];
  }


  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::getInstanceBreakPoints
// Purpose       : Add break point for anticipated digital output transition
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 02/07/06
//-----------------------------------------------------------------------------
bool Instance::getInstanceBreakPoints
( vector<N_UTL_BreakPoint> & breakPointTimes )
{
  if (breakTime > getSolverState().currTime)
  {
    breakPointTimes.push_back(breakTime);
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEQVector
//
// Purpose       : Loads the Q-vector contributions for a single
//                 digital instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/05/06
//-----------------------------------------------------------------------------
bool Instance::loadDAEQVector ()
{
  bool bsuccess = true;
  int i;

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline2("---------------------");
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << dashedline2 << endl;
    cout << "  Instance::loadDAEQVector" << endl;
    cout << "  name = " << getName() <<endl;
  }
#endif

  for (i=0 ; i<numOutput ; ++i)
  {
#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      cout << "  qlo[" << i << "] = " << qlo[i] << endl;
      cout << "  qhi[" << i << "] = " << qhi[i] << endl;
    }
#endif
    if (li_Lo >= 0)
    {

      (*extData.daeQVectorPtr)[li_Lo] += qlo[i];
    }
    if (li_Hi >= 0)
    {

      (*extData.daeQVectorPtr)[li_Hi] += qhi[i];
    }


    (*extData.daeQVectorPtr)[li_Out[i]] -= qlo[i];

    (*extData.daeQVectorPtr)[li_Out[i]] -= qhi[i];
  }

  for (i=0 ; i<numInput ; ++i)
  {
#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      cout << "  qref[" << i << "] = " << qref[i] << endl;
    }
#endif
    if (li_Ref >= 0)
    {

      (*extData.daeQVectorPtr)[li_Ref] += qref[i];
    }


    (*extData.daeQVectorPtr)[li_Inp[i]] -= qref[i];
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 digital instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/05/06
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  bool bsuccess = true;
  double coef = 0.0;
  double Vpos = 0.0;
  double Vneg = 0.0;
  double v_tmp = 0.0;
  int i;

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline2("---------------------");
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << dashedline2 << endl;
    cout << "  Instance::loadDAEFVector" << endl;
    cout << "  name = " << getName() <<endl;
  }
#endif

  for (i=0 ; i<numOutput ; ++i)
  {
#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      cout << "  rilo[" << i << "] = " << rilo[i] << endl;
      cout << "  rihi[" << i << "] = " << rihi[i] << endl;
    }
#endif
    if (li_Lo >= 0)
    {

      (*extData.daeFVectorPtr)[li_Lo] += rilo[i];
    }
    if (li_Hi >= 0)
    {

      (*extData.daeFVectorPtr)[li_Hi] += rihi[i];
    }


    (*extData.daeFVectorPtr)[li_Out[i]] -= rilo[i];

    (*extData.daeFVectorPtr)[li_Out[i]] -= rihi[i];
  }

  for (i=0 ; i<numInput ; ++i)
  {
#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      cout << "  riref[" << i << "] = " << riref[i] << endl;
    }
#endif
    if (li_Ref >= 0)
    {

      (*extData.daeFVectorPtr)[li_Ref] += riref[i];
    }


    (*extData.daeFVectorPtr)[li_Inp[i]] -= riref[i];
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << dashedline2 << endl;
  }
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx
//
// Purpose       : Loads the dQdx-matrix contributions for a single
//                 digital instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/05/06
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  bool bsuccess = true;
  int i;

  N_LAS_Matrix * dQdxMatPtr = extData.dQdxMatrixPtr;

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline2("---------------------");
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << dashedline2 <<endl;
    cout << "  Instance::loadDAEdQdx" << endl;
    cout << "  name = " << getName() << endl;
  }
#endif

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << "\nLoading DIGITAL dQdx matrix\n";
    cout << "Capacitance lo: " << model_.clo << endl;
    cout << "Capacitance hi: " << model_.chi << endl;
    cout << "Capacitance load: " << model_.cload << endl;
    cout << "DONE DIGITAL dQdx matrix LOAD\n";
  }
#endif

  for (i=0 ; i<numInput ; ++i)
  {
    if (row_Ref >= 0)
    {

      (*dQdxMatPtr)[li_Ref][li_jac_Ref[i][1]] += model_.cload;

      (*dQdxMatPtr)[li_Ref][li_jac_Ref[i][2]] -= model_.cload;

      (*dQdxMatPtr)[li_Inp[i]][li_jac_Ref[i][4]] -= model_.cload;

      (*dQdxMatPtr)[li_Inp[i]][li_jac_Ref[i][5]] += model_.cload;
    }
    else
    {

      (*dQdxMatPtr)[li_Inp[i]][li_jac_Ref[i][1]] += model_.cload;
    }
  }

  for (i=0 ; i<numOutput ; ++i)
  {
    if (row_Lo >= 0)
    {

      (*dQdxMatPtr)[li_Lo][li_jac_Lo[i][1]] += model_.clo;

      (*dQdxMatPtr)[li_Lo][li_jac_Lo[i][2]] -= model_.clo;

      (*dQdxMatPtr)[li_Out[i]][li_jac_Lo[i][4]] -= model_.clo;

      (*dQdxMatPtr)[li_Out[i]][li_jac_Lo[i][5]] += model_.clo;
    }
    else
    {

      (*dQdxMatPtr)[li_Out[i]][li_jac_Lo[i][1]] += model_.clo;
    }
    if (row_Hi >= 0)
    {

      (*dQdxMatPtr)[li_Hi][li_jac_Hi[i][1]] += model_.chi;

      (*dQdxMatPtr)[li_Hi][li_jac_Hi[i][2]] -= model_.chi;

      (*dQdxMatPtr)[li_Out[i]][li_jac_Hi[i][4]] -= model_.chi;

      (*dQdxMatPtr)[li_Out[i]][li_jac_Hi[i][5]] += model_.chi;
    }
    else
    {

      (*dQdxMatPtr)[li_Out[i]][li_jac_Hi[i][1]] += model_.chi;
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 digital instance.
//
// Special Notes : See the special notes for loadDAEFVector.
//                 For digital devices this is the contribution of the resistors
//
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/05/06
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  bool bsuccess = true;
  N_LAS_Matrix * dFdxMatPtr = extData.dFdxMatrixPtr;
  int i;

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline2("---------------------");
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << dashedline2 <<endl;
    cout << "  Instance::loadDAEdFdx" << endl;
  }
#endif

  for (i=0 ; i<numInput ; ++i)
  {
    if (row_Ref >= 0)
    {

      (*dFdxMatPtr)[li_Ref][li_jac_Ref[i][1]] += model_.gload;

      (*dFdxMatPtr)[li_Ref][li_jac_Ref[i][2]] -= model_.gload;

      (*dFdxMatPtr)[li_Inp[i]][li_jac_Ref[i][4]] -= model_.gload;

      (*dFdxMatPtr)[li_Inp[i]][li_jac_Ref[i][5]] += model_.gload;
    }
    else
    {

      (*dFdxMatPtr)[li_Inp[i]][li_jac_Ref[i][1]] += model_.gload;
    }
  }

  for (i=0 ; i<numOutput ; ++i)
  {
#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      cout << "  glo[" << i << "] = " << glo[i] << endl;
      cout << "  ghi[" << i << "] = " << ghi[i] << endl;
    }
#endif
    if (row_Lo >= 0)
    {

      (*dFdxMatPtr)[li_Lo][li_jac_Lo[i][1]] += glo[i];

      (*dFdxMatPtr)[li_Lo][li_jac_Lo[i][2]] -= glo[i];

      (*dFdxMatPtr)[li_Out[i]][li_jac_Lo[i][4]] -= glo[i];

      (*dFdxMatPtr)[li_Out[i]][li_jac_Lo[i][5]] += glo[i];
    }
    else
    {

      (*dFdxMatPtr)[li_Out[i]][li_jac_Lo[i][1]] += glo[i];
    }
    if (row_Hi >= 0)
    {

      (*dFdxMatPtr)[li_Hi][li_jac_Hi[i][1]] += ghi[i];

      (*dFdxMatPtr)[li_Hi][li_jac_Hi[i][2]] -= ghi[i];

      (*dFdxMatPtr)[li_Out[i]][li_jac_Hi[i][4]] -= ghi[i];

      (*dFdxMatPtr)[li_Out[i]][li_jac_Hi[i][5]] += ghi[i];
    }
    else
    {

      (*dFdxMatPtr)[li_Out[i]][li_jac_Hi[i][1]] += ghi[i];
    }
  }

  return bsuccess;
}

// Class Model
//-----------------------------------------------------------------------------
// Function      : Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/05/06
//-----------------------------------------------------------------------------
bool Model::processParams (string param)
{

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
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/05/06
//-----------------------------------------------------------------------------

Model::Model(const ModelBlock & MB,
             SolverState & ss1,
             DeviceOptions & do1)
  : DeviceModel(MB,ss1,do1)

{

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to .model line and constant defaults from metadata:
  setModParams (MB.params);

  // Set any non-constant parameter defaults:

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:
  if (rload == 0)
  {
    string msg("Zero load resistance in inputs");
    std::ostringstream oss;
    oss << "Error in " << netlistLocation() << "\n" << msg;
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL, oss.str());
  }
  gload = 1/rload;

  processParams ();
}

//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/05/06
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
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/05/06
//-----------------------------------------------------------------------------

std::ostream &Model::printOutInstances(std::ostream &os) const
{
  vector<Instance*>::const_iterator iter;
  vector<Instance*>::const_iterator first = instanceContainer.begin();
  vector<Instance*>::const_iterator last  = instanceContainer.end();

  int i,isize;

  isize = instanceContainer.size();
  os << endl;
  os << "Number of digital instances: " << isize << endl;
  os << "    name\t\tmodelName\tParameters" << endl;

  for (i = 0, iter = first; iter != last; ++iter, ++i)
  {
    os << "  " << i << ": " << (*iter)->getName() << "\t";
    os << (*iter)->getModelName();
    os << endl;
  }

  os << endl;

  return os;
}

} // namespace Digital
} // namespace Device
} // namespace Xyce

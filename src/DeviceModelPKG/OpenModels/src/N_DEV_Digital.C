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
// Revision Number: $Revision: 1.97.2.5 $
//
// Revision Date  : $Date: 2014/03/14 18:03:21 $
//
// Current Owner  : $Author: peshola $
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
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_Digital.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_Message.h>
#include <N_ERH_ErrorMgr.h>
#include <N_UTL_BreakPoint.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>

namespace Xyce {
namespace Device {


namespace Digital {

// enables debug output for just the digital devices
//static const int DEBUG_DEVICE = 1;

void Traits::loadInstanceParameters(ParametricData<Digital::Instance> &p)
{
// Set up double precision variables:

  // Set up non-double precision variables:
  p.addPar ("IC1", false, &Digital::Instance::ic1)
    .setUnit(U_LOGIC)
    .setDescription("Vector of initial values for output(s)");
  p.addPar ("IC2", false, &Digital::Instance::ic2);

  p.makeVector ("IC", 2);
}

void Traits::loadModelParameters(ParametricData<Digital::Model> &p)
{
  // Set up double precision variables:
  p.addPar ("VLO", 0., &Digital::Model::vlo)
    .setUnit(U_VOLT)
    .setDescription("Internal low state supply voltage");

  p.addPar ("VHI", 0., &Digital::Model::vhi)
    .setUnit(U_VOLT)
    .setDescription("Internal high state supply voltage");

  p.addPar ("VREF", 0., &Digital::Model::vref)
    .setUnit(U_VOLT)
    .setDescription("Internal reference voltage for inputs");

  p.addPar ("CLO", 1.e-6, &Digital::Model::clo)
    .setUnit(U_FARAD)
    .setDescription("Capacitance between output node and low reference");

  p.addPar ("CHI", 1.e-6, &Digital::Model::chi)
    .setUnit(U_FARAD)
    .setDescription("Capacitance between output node and high reference");

  p.addPar ("CLOAD", 1.e-6, &Digital::Model::cload)
    .setUnit(U_FARAD)
    .setDescription("Capacitance between input node and input reference");

  p.addPar ("RLOAD", 1000., &Digital::Model::rload)
    .setUnit(U_OHM)
    .setDescription("Resistance between input node and input reference");

  p.addPar ("S0RLO", 100., &Digital::Model::s0rlo)
    .setUnit(U_OHM)
    .setDescription("Low state resistance between output node and low reference");

  p.addPar ("S0RHI", 100., &Digital::Model::s0rhi)
    .setUnit(U_OHM)
    .setDescription("Low state resitance between output node and high reference");

  p.addPar ("S0TSW", 1.e-8, &Digital::Model::s0tsw)
    .setUnit(U_SECOND)
    .setDescription("Switching time transition to low state");

  p.addPar ("S0VLO", -1.5, &Digital::Model::s0vlo)
    .setUnit(U_VOLT)
    .setDescription("Minimum voltage to switch to low state");

  p.addPar ("S0VHI", 1.7, &Digital::Model::s0vhi)
    .setUnit(U_VOLT)
    .setDescription("Maximum voltage to switch to low state");

  p.addPar ("S1RLO", 100., &Digital::Model::s1rlo)
    .setUnit(U_OHM)
    .setDescription("High state resistance between output node and low reference");

  p.addPar ("S1RHI", 100., &Digital::Model::s1rhi)
    .setUnit(U_OHM)
    .setDescription("High state resistance between output node and high reference");

  p.addPar ("S1TSW", 1.e-8, &Digital::Model::s1tsw)
    .setUnit(U_SECOND)
    .setDescription("Switching time transition to high state");

  p.addPar ("S1VLO", 0.9, &Digital::Model::s1vlo)
    .setUnit(U_VOLT)
    .setDescription("Minimum voltage to switch to high state");

  p.addPar ("S1VHI", 7.0, &Digital::Model::s1vhi)
    .setUnit(U_VOLT)
    .setDescription("Maximum voltage to switch to high state");

  p.addPar ("DELAY", 1.e-8, &Digital::Model::delay)
    .setUnit(U_SECOND)
    .setDescription("Delay time of device");
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
bool Instance::processParams ()
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
  const Configuration & configuration,
  const InstanceBlock & IB,
  Model & Diter,
  const FactoryBlock &  factory_block)
  : DeviceInstance(IB, configuration.getInstanceParameters(), factory_block),
    model_(Diter),
    li_Lo(-1),
    li_Hi(-1),
    li_Ref(-1),
    row_Lo(-1),
    row_Hi(-1),
    row_Ref(-1),
    breakTime(0.)
{
  int i, tokenCount = 0, dev_numInputs = 0;
  tokenCount = count(getName().begin(),getName().end(),'%'); 
  if (tokenCount == 2 || tokenCount == 3)
  {
    int p1, p2, p3;
    p1 = getName().find_first_of('%');
    p2 = p1 + 1 + getName().substr(p1+1).find_first_of('%');
    // parse U devices with a variable number of inputs
    if (tokenCount == 3)
    {   
      p3 = getName().find_last_of('%');
      if (p3 != getName().size()-1){
        dev_numInputs = std::atoi(getName().substr(p3+1).c_str());
      }
    }
 
    numExtVars = 0;
    // Code required to support both Y-style and U-style digital devices.
    // Y digital devices are now deprecated.
    std::string dev_letter = getName().substr(p1-1,1);
    if (dev_letter == "U"){
      // For U devices, DPWR and DGND are always specified on the instance line.
      // Warning message if VHI, VLO or VREF are in the model card.
      li_Lo = 0;
      li_Hi = 0;
      li_Ref = 0;
      numExtVars += 2;
      if (model_.given("VLO"))
      { 
        UserWarning(*this)<< "VLO model parameter ignored for U digital device";
      }
      if (model_.given("VHI"))
      { 
        UserWarning(*this)<< "VHI model parameter ignored for U digital device";
      }
      if (model_.given("VREF"))
      { 
        UserWarning(*this)<< "VREF model parameter ignored for U digital device";
      }
    }
    else if (dev_letter == "Y")
    {
      // legacy code required to support VLO, VHI and VREF variables
      // being on the instance line rather than in model card in Y devices
      UserWarning(*this)<< "Y digital device (" << getName() << ") is deprecated. Consider using U device instead."; 
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
    }  
    else
    {
      UserError0(*this) << "Digital device letter must be Y or U: " << getName();
    }

    // Configure number of inputs/outputs for each device
    // Y devices are limited to 2 inputs.
    std::string dev_type(getName().substr(p1+1,p2-p1-1));
    if (dev_type == "NOT" || dev_type == "INV")
    {
      if (dev_type == "NOT")
      {
        UserWarning(*this)<< "NOT gate type (" << getName() << ") is deprecated. Consider using INV instead."; 
      }
      numInput = 1;
      numOutput = 1;
      gate = INV;
    }
    else if (dev_type == "AND")
    {
      (dev_letter == "Y") ? (numInput = 2) : (numInput = dev_numInputs);
      numOutput = 1;
      gate = AND;
    }
    else if (dev_type == "NAND")
    {
      (dev_letter == "Y") ? (numInput = 2) : (numInput = dev_numInputs);
      numOutput = 1;
      gate = NAND;
    }
    else if (dev_type == "OR")
    {
      (dev_letter == "Y") ? (numInput = 2) : (numInput = dev_numInputs);
      numOutput = 1;
      gate = OR;
    }
    else if (dev_type == "NOR")
    {
      (dev_letter == "Y") ? (numInput = 2) : (numInput = dev_numInputs);
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
    else if (dev_type == "DFF")
    {
      numInput = 4;  //PREB, CLRB, clock, data
      numOutput = 2; //Q, Q_bar
      gate = DFF;
    }
    // DLTCH device still under development
    /*else if (dev_type == "DLTCH")
    {
      numInput = 4;  // PREB, CLRB, Enable, Data
      numOutput = 2; // Q, Q_bar
      gate = DLTCH;
    }*/
    else if (dev_type == "BUF")
    {
      numInput = 1;
      numOutput = 1;
      gate = BUF;
    }
    else
    {
      UserError0(*this) << "Unknown digital device type " << dev_type;
    }
  }
  else
  {
    std::string msg("Internal error in digital device name: ");
    msg += getName();
    std::ostringstream oss;
    oss << "Error in " << netlistLocation() << "\n" << msg;
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, oss.str());
  }

  int iBase = numExtVars;
  int oBase = numExtVars + numInput;
  numExtVars += numInput + numOutput;

  if (DEBUG_DEVICE && getDeviceOptions().debugLevel > 0)
  {
    Xyce::dout() << "Digital Device " << getName() << " has iBase = " <<
        iBase << ", oBase = " << oBase << ", numExtVars = " <<
        numExtVars << std::endl;
  }

  // catch cases of AND, NAND, OR or NOR gate with only one input specified
  // or the number of nodes on the instance line does not match the 
  // (N) value specified as part of the gate type (e.g, AND(4))
  if ((gate == AND) || (gate == NAND) || (gate == OR) || (gate == NOR))
  {
    if (numInput == 1)
    {
      UserError0(*this) << "this device must have more than one input.";
    }
    if ( (dev_numInputs != 0) && (IB.numExtVars - iBase - numOutput != dev_numInputs) )
    {
      std::cout << IB.numExtVars << " " << dev_numInputs << " " << numExtVars << std::endl; 
      UserError0(*this) << "too few I/O nodes on instance line.";
    }
  }
  // catch case where gates with a fixed number of inputs have
  // wrong number of inputs
  if (numExtVars != IB.numExtVars)
  {
    UserError0(*this) << "Incorrect number of nodes in digital device "
         << ".  Found " << IB.numExtVars << ", should be " << numExtVars;
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

  jacStamp.resize(numExtVars);
  int row = 0;
  // Code required to support both Y-style and U-style digital devices.
  // Y digital devices are now deprecated.
  std::string dev_letter = getDeviceLetter();
  if (dev_letter == "U")
  {
    // digital power and digital ground node are always on
    // the instance line for a U device.  The low reference
    // voltage for inputs is assumed to be the same as the
    // digital ground node.
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
  
    row_Ref = row;
    for (i=0 ; i<numInput ; ++i)
    {
      li_jac_Ref[i].push_back(row);
      li_jac_Ref[i].push_back(0);
      li_jac_Ref[i].push_back(1);
      li_jac_Ref[i].push_back(iBase+i);
      li_jac_Ref[i].push_back(jacStamp[iBase+i].size());
      jacStamp[row].push_back(iBase+i);
      jacStamp[iBase+i].push_back(row);
    }
  }
  else if (dev_letter == "Y") 
  {
    // output low and output high reference voltages and input 
    // low reference voltage are optionally on the instance line 
    // for Y devices.
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
  }
  else
  {
    UserError0(*this) << "Digital device letter must be Y or U: " << getName();
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
void Instance::registerLIDs( const std::vector<int> & intLIDVecRef,
                             const std::vector<int> & extLIDVecRef)

{
  std::string msg;

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
  // Code required to support both Y-style and U-style digital devices.
  // Y digital devices are now deprecated
  std::string dev_letter = getDeviceLetter();
  if (dev_letter == "U")
  {
    // ordering on U-device instance line is dig_power_node (DPWR) then 
    // dig_ground_node (DGND).  Assume that input low-reference
    // voltage (VREF in Y devices) is equal to DGND.
    li_Hi = extLIDVec[i++];
    li_Lo = extLIDVec[i++];
    li_Ref = li_Lo;
  }
  else if (dev_letter == "Y")
  {
    // for Y devices, ordering of output high/low reference nodes (VHI/VLO)
    // on instance line is reversed.  Input low-reference voltage (VREF), 
    // VHI and VLO can either be on the instance line or in the model card.
    if (li_Lo == 0)
      li_Lo = extLIDVec[i++];
    if (li_Hi == 0)
      li_Hi = extLIDVec[i++];
    if (li_Ref == 0)
      li_Ref = extLIDVec[i++];
  }
  else
  {
    UserError0(*this) << "Digital device letter must be Y or U: " << getName();
  } 

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
void Instance::registerStateLIDs( const std::vector<int> & staLIDVecRef)
{
  std::string msg;

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
std::map<int,std::string> & Instance::getIntNameMap ()
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
const std::vector< std::vector<int> > & Instance::jacobianStamp() const
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
void Instance::registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );
  int i;
  int iBase = numExtVars - numInput - numOutput;
  int oBase = iBase + numInput;
  int lo_present, hi_present;

  (row_Lo < 0) ? (lo_present = 0) : (lo_present = 1);
  (row_Hi < 0) ? (hi_present = 0) : (hi_present = 1);

  // Code required to support both Y-style and U-style digital devices.
  // Y digital devices are now deprecated.
  std::string dev_letter = getDeviceLetter();
  
  if (dev_letter == "U")
  {
    if ( (row_Lo == -1) || (row_Hi == -1) || (row_Ref == -1) )
    {
       UserError0(*this) << "Internal error in Instance::registerJacLIDs() for " << getName();
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

      for (i=0 ; i<numOutput ; ++i)
      {
        li_jac_Lo[i][1] = jacLIDVec[li_jac_Lo[i][0]][li_jac_Lo[i][1]];
        li_jac_Lo[i][2] = jacLIDVec[li_jac_Lo[i][0]][li_jac_Lo[i][2]];
        li_jac_Lo[i][4] = jacLIDVec[li_jac_Lo[i][3]][li_jac_Lo[i][4]];
        li_jac_Lo[i][5] = jacLIDVec[li_jac_Lo[i][3]][li_jac_Lo[i][5]];
      }
  
      for (i=0 ; i<numOutput ; ++i)
      {
        li_jac_Hi[i][1] = jacLIDVec[li_jac_Hi[i][0]][li_jac_Hi[i][1]];
        li_jac_Hi[i][2] = jacLIDVec[li_jac_Hi[i][0]][li_jac_Hi[i][2]];
        li_jac_Hi[i][4] = jacLIDVec[li_jac_Hi[i][3]][li_jac_Hi[i][4]];
        li_jac_Hi[i][5] = jacLIDVec[li_jac_Hi[i][3]][li_jac_Hi[i][5]];
      }
    }
  }
  else if (dev_letter == "Y") 
  {
    // row_Ref == -1 means that VREF is in the model card
    if (row_Ref == -1)
    {
      for (i=0 ; i<numInput ; ++i)
      {
        li_jac_Ref[i].push_back(jacLIDVec[li_jac_Ref[i][0]][0]);
      }
    }
    else
    {
      // this for loop is identical for both U and Y devices
      for (i=0 ; i<numInput ; ++i)
      {
        li_jac_Ref[i][1] = jacLIDVec[li_jac_Ref[i][0]][li_jac_Ref[i][1]];
        li_jac_Ref[i][2] = jacLIDVec[li_jac_Ref[i][0]][li_jac_Ref[i][2]];
        li_jac_Ref[i][4] = jacLIDVec[li_jac_Ref[i][3]][li_jac_Ref[i][4]];
        li_jac_Ref[i][5] = jacLIDVec[li_jac_Ref[i][3]][li_jac_Ref[i][5]];
      }
    }
    
    // row_Lo == -1 means that VLO is in the model card
    if (row_Lo == -1)
    {
      for (i=0 ; i<numOutput ; ++i)
      {
        li_jac_Lo[i].push_back(jacLIDVec[li_jac_Lo[i][0]][hi_present]);
      }
    }
    else
    {
      // this for loop is identical for both U and Y devices
      for (i=0 ; i<numOutput ; ++i)
      {
        li_jac_Lo[i][1] = jacLIDVec[li_jac_Lo[i][0]][li_jac_Lo[i][1]];
        li_jac_Lo[i][2] = jacLIDVec[li_jac_Lo[i][0]][li_jac_Lo[i][2]];
        li_jac_Lo[i][4] = jacLIDVec[li_jac_Lo[i][3]][li_jac_Lo[i][4]];
        li_jac_Lo[i][5] = jacLIDVec[li_jac_Lo[i][3]][li_jac_Lo[i][5]];
      }
    }

    // row_Hi == -1 means that VHI is in the model card
    if (row_Hi == -1)
    {
      for (i=0 ; i<numOutput ; ++i)
      {
        li_jac_Hi[i].push_back(jacLIDVec[li_jac_Hi[i][0]][lo_present]);
      }
    }
    else
    {
      // this for loop is identical for both U and Y devices
      for (i=0 ; i<numOutput ; ++i)
      {
        li_jac_Hi[i][1] = jacLIDVec[li_jac_Hi[i][0]][li_jac_Hi[i][1]];
        li_jac_Hi[i][2] = jacLIDVec[li_jac_Hi[i][0]][li_jac_Hi[i][2]];
        li_jac_Hi[i][4] = jacLIDVec[li_jac_Hi[i][3]][li_jac_Hi[i][4]];
        li_jac_Hi[i][5] = jacLIDVec[li_jac_Hi[i][3]][li_jac_Hi[i][5]];
      }
    }
  }  
  else
  {
    UserError0(*this) << "Digital device letter must be Y or U: " << getName();
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
  //std::vector<bool> changeState;  //Genie 110812
  bool toPrint = false; //Genie 022713

  int i;

  N_LAS_Vector & solVector = *(extData.nextSolVectorPtr);
  N_LAS_Vector & staVector = *(extData.nextStaVectorPtr);
  N_LAS_Vector & oldStaVector = *(extData.currStaVectorPtr);
  N_LAS_Vector & oldSolVector = *(extData.currSolVectorPtr);

  // The convention in this device is to consider the 'positive' side of the
  // capacitors as the Vsrc or node that supplies the voltage, or for the
  // output, vref. The 'positive' voltages are thus the same for all
  // input/outputs
  (li_Lo >= 0) ? (v_poslo = solVector[li_Lo]) : (v_poslo = model_.vlo);

  // If VHI is not specified, li_Hi >=0 and v_poshi is derived from
  // CLOAD/RLOAD by solVector.  If VHI is specified, li_Hi == -1 and
  // v_poshi is specified by the netlist model card.
  (li_Hi >= 0) ? (v_poshi = solVector[li_Hi]) : (v_poshi = model_.vhi);
  (li_Ref >= 0) ? (v_posref = solVector[li_Ref]) : (v_posref = model_.vref);

  lastT = 0;

  for (i=0 ; i<numInput ; ++i)
  {
    //initialize
    currentState = static_cast <int> (oldStaVector[li_currentStateInp[i]]);
    transitionTime = oldStaVector[li_transitionTimeInp[i]];
    changeState = false; //Genie 022013 Clear the memory of changeState.

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
      (-vcapref[i] < model_.s0vhi) ? (currentState = 0) : (currentState = 1);

      oldStaVector[li_currentStateInp[i]] = currentState;
      oldStaVector[li_transitionTimeInp[i]] = transitionTime;
    }
 
    iTime[i] = transitionTime;

    staVector[li_transitionTimeInp[i]] = transitionTime;

    if (currentState == 0)
    {
      inpL[i] = false;
      if ((-vcapref[i] > model_.s0vhi) && (-vcapref[i] > model_.s1vlo))
      {
          currentState = 1;
          changeState = true;
          if (DEBUG_DEVICE && getDeviceOptions().debugLevel > 0)
	  {
	    Xyce::dout() << "Device " << getName() << " changed state from 0 to 1 at time " << getSolverState().currTime << std::endl;
	  }

          if (gate == DFF && i == 2)  // i==2 -> clk
          {   
            clocking = true;  // clock of  DFF changes state
	  }
      }
    }
    else
    {
      inpL[i] = true;
      if ((-vcapref[i] < model_.s1vlo) && (-vcapref[i] < model_.s0vhi))
      {
          currentState = 0;
          changeState = true;
          if (DEBUG_DEVICE && getDeviceOptions().debugLevel > 0)
	  {
	    Xyce::dout() << "Device " << getName() << " changed state from 1 to 0 at time " << getSolverState().currTime << std::endl;
	  }
	  if (gate == DFF && i == 2)  //Genie 111212
	  {
            clocking = true;  // clock of  DFF changes state
	  }
      }
    }

    if (changeState)
    {
      double vOld, del;

      inpL[i] = (currentState == 1);
      (li_Ref >= 0) ? (vOld = oldSolVector[li_Ref]) : (vOld = model_.vref);
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

  if (gate == INV)
  {
    outL[0] = !inpL[0];
    oTime[0] = lastT+model_.delay;
  }
  else if (gate == AND)
  {
    outL[0] = !(count(inpL.begin(),inpL.end(),false) > 0);
    oTime[0] = lastT+model_.delay;
  }
  else if (gate == NAND)
  {
    outL[0] = (count(inpL.begin(),inpL.end(),false) > 0);
    oTime[0] = lastT+model_.delay;
  }
  else if (gate == OR)
  {
    outL[0] = (count(inpL.begin(),inpL.end(),true) > 0);
    oTime[0] = lastT+model_.delay;
  }
  else if (gate == NOR)
  {
    outL[0] = !(count(inpL.begin(),inpL.end(),true) > 0);
    oTime[0] = lastT+model_.delay;
  }
  else if (gate == ADD)
  {
    outL[0] = inpL[0] ^ inpL[1] ^ inpL[2];
    // carry-out sum
    outL[1] = (inpL[0] & inpL[1]) | (inpL[1] & inpL[2]) | (inpL[0] & inpL[2]);

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
  else if (gate == BUF)
  {
    outL[0] = inpL[0];
    oTime[0] = lastT+model_.delay;
  }
  // DLTCH device still under development.  This else-if
  // clause not reachable in Xyce 6.1
  else if (gate == DLTCH)
  { // DLTCH: in0: PREB, in1: CLRB, in2: enable, in3: data
    // DLTCH: out0: Q, out1: Q_bar
    if ((inpL[0] == 1) && (inpL[1] == 0))
    {
      outL[0] = 0;
      outL[1] = 1;
      oTime[0] = lastT+model_.delay;
      oTime[1] = lastT+model_.delay;
    }
    else if ((inpL[0] == 0) && (inpL[1] == 1))
    {
      outL[0] = 1;
      outL[1] = 0;
      oTime[0] = lastT+model_.delay;
      oTime[1] = lastT+model_.delay;
    }
     else if ((inpL[0] == 0) && (inpL[1] == 0))
    { // this state is unstable, and needs further testing
      outL[0] = 1;
      outL[1] = 1;
      oTime[0] = lastT+model_.delay;
      oTime[1] = lastT+model_.delay;
    }
    else if (inpL[2] == 1)
    { // enable line, PREB and CLRB are TRUE
      outL[0] = inpL[3];
      outL[1] = !inpL[3];
      oTime[0] = lastT+model_.delay;
      oTime[1] = lastT+model_.delay;
    }
    else
    {
      // no op.  Keep outputs latched in current state
      // this statement deals with the startup condition,
      // when no IC's were specified.  However, it may cause
      // trouble when transitioning from the PREB/CLRB = 0
      // state
      outL[1] = !outL[0];
    }
  }
  else if (gate == DFF)
  { // DFF: in0: PREB, in1: CLRB, in2: clock, in3: data
    // DFF: out0: Q, out1: Q_bar
    // CD4013B: set = !PREB, reset = !CLRB
    // CD4013B: in0: set, in1: reset, in2: clock, in3: data
    // CD4013B: out0: Q, out1: Q_bar
    if (clocking && inpL[2] ==1) //clock rising edge 0->1
    {
        if (inpL[0] == 1 && inpL[1] == 1) //PREB = CLRB = 1
        {
            outL[0] = inpL[3];  //Q = D
            outL[1] = !(inpL[3]); //Q_bar = !D
	}
    }
    else if (clocking && inpL[2] ==0) //clock falling edge 1->0
    {
 	if (inpL[0] == 1 && inpL[1] == 1) //PREB = CLRB = 1
        {
            outL[0] = oldStaVector[li_currentStateOut[0]]; //no change
            outL[1] = oldStaVector[li_currentStateOut[1]]; //no change
	}
    }
    else // no clock change
    {
      if (inpL[0] == 1 && inpL[1] == 0) //PREB=1, CLRB = 0
        {
            outL[0] = 0;
            outL[1] = 1;
	}
	else if (inpL[0] == 0 && inpL[1] == 1) //PREB = 0, CLRB = 1
        {
            outL[0] = 1;
            outL[1] = 0;
	}
	else if (inpL[0] == 0 && inpL[1] == 0) //PREB = CLRB = 0
        {
            outL[0] = 1;
            outL[1] = 1;
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
      else
      {
        std::string msg("Insufficient initial conditions supported in digital device");
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
      }

      oldStaVector[li_currentStateOut[i]] = outL[i]?1:0;
      oldStaVector[li_transitionTimeOut[i]] = time;
    }

    //current logic state of output nodes
    currentState = static_cast <int> (oldStaVector[li_currentStateOut[i]]);
    transitionTime = oldStaVector[li_transitionTimeOut[i]];

    if (currentState == 1)
      curr = true;
    else
      curr = false;

    if (curr != outL[i]) // This is executed when scopFlag is false
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
      std::string msg("Instance::updateSecondaryState: unrecognized state");
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
        // This is a simple linear transition.  Since there is a
        // breakpoint at the start of the transition it is OK to
        // have a discontinuity there.
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

    rilo[i] = glo[i]*(v_poslo-v_neg);
    rihi[i] = ghi[i]*(v_poshi-v_neg);

    vcaplo[i] = v_poslo-v_neg;
    vcaphi[i] = v_poshi-v_neg;

    // Obtain the "current"  value for the charge stored in the capacitors.
    qlo[i] = model_.clo*vcaplo[i]; 
    qhi[i] = model_.chi*vcaphi[i];

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
( std::vector<N_UTL_BreakPoint> & breakPointTimes )
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

  if (DEBUG_DEVICE && getDeviceOptions().debugLevel > 2 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << subsection_divider << std::endl;
    Xyce::dout() << "  Instance::loadDAEQVector" << std::endl;
    Xyce::dout() << "  name = " << getName() <<std::endl;
  }

  for (i=0 ; i<numOutput ; ++i)
  {
    if (DEBUG_DEVICE && getDeviceOptions().debugLevel > 2 && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "  qlo[" << i << "] = " << qlo[i] << std::endl;
      Xyce::dout() << "  qhi[" << i << "] = " << qhi[i] << std::endl;
    }

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
    if (DEBUG_DEVICE && getDeviceOptions().debugLevel > 2 && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "  qref[" << i << "] = " << qref[i] << std::endl;
    }

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

  if (DEBUG_DEVICE && getDeviceOptions().debugLevel > 2 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << subsection_divider << std::endl;
    Xyce::dout() << "  Instance::loadDAEFVector" << std::endl;
    Xyce::dout() << "  name = " << getName() <<std::endl;
  }

  for (i=0 ; i<numOutput ; ++i)
  {
    if (DEBUG_DEVICE && getDeviceOptions().debugLevel > 2 && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "  rilo[" << i << "] = " << rilo[i] << std::endl;
      Xyce::dout() << "  rihi[" << i << "] = " << rihi[i] << std::endl;
    }

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
    if (DEBUG_DEVICE && getDeviceOptions().debugLevel > 2 && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "  riref[" << i << "] = " << riref[i] << std::endl;
    }

    if (li_Ref >= 0)
    {

      (*extData.daeFVectorPtr)[li_Ref] += riref[i];
    }


    (*extData.daeFVectorPtr)[li_Inp[i]] -= riref[i];
  }

  if (DEBUG_DEVICE && getDeviceOptions().debugLevel > 2 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << subsection_divider << std::endl;
  }

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

  if (DEBUG_DEVICE && getDeviceOptions().debugLevel > 2 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << subsection_divider <<std::endl;
    Xyce::dout() << "  Instance::loadDAEdQdx" << std::endl;
    Xyce::dout() << "  name = " << getName() << std::endl;
  }

  if (DEBUG_DEVICE && getDeviceOptions().debugLevel > 2 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "\nLoading DIGITAL dQdx matrix\n";
    Xyce::dout() << "Capacitance lo: " << model_.clo << std::endl;
    Xyce::dout() << "Capacitance hi: " << model_.chi << std::endl;
    Xyce::dout() << "Capacitance load: " << model_.cload << std::endl;
    Xyce::dout() << "DONE DIGITAL dQdx matrix LOAD\n";
  }

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

  if (DEBUG_DEVICE && getDeviceOptions().debugLevel > 2 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << subsection_divider <<std::endl;
    Xyce::dout() << "  Instance::loadDAEdFdx" << std::endl;
  }

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
    if (DEBUG_DEVICE && getDeviceOptions().debugLevel > 2 && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "  glo[" << i << "] = " << glo[i] << std::endl;
      Xyce::dout() << "  ghi[" << i << "] = " << ghi[i] << std::endl;
    }

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
bool Model::processParams ()
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
bool Model::processInstanceParams()
{

  std::vector<Instance*>::iterator iter;
  std::vector<Instance*>::iterator first = instanceContainer.begin();
  std::vector<Instance*>::iterator last  = instanceContainer.end();

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

Model::Model(
  const Configuration & configuration,
  const ModelBlock &    MB,
  const FactoryBlock &  factory_block)
  : DeviceModel(MB, configuration.getModelParameters(), factory_block)

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
    UserError0(*this) << "Zero load resistance in inputs";
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
  std::vector<Instance*>::iterator iter;
  std::vector<Instance*>::iterator first = instanceContainer.begin();
  std::vector<Instance*>::iterator last  = instanceContainer.end();

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
  std::vector<Instance*>::const_iterator iter;
  std::vector<Instance*>::const_iterator first = instanceContainer.begin();
  std::vector<Instance*>::const_iterator last  = instanceContainer.end();

  int i,isize;

  isize = instanceContainer.size();
  os << std::endl;
  os << "Number of digital instances: " << isize << std::endl;
  os << "    name\t\tmodelName\tParameters" << std::endl;

  for (i = 0, iter = first; iter != last; ++iter, ++i)
  {
    os << "  " << i << ": " << (*iter)->getName() << "\t";
    os << getName();
    os << std::endl;
  }

  os << std::endl;

  return os;
}

//-----------------------------------------------------------------------------
// Function      : Model::forEachInstance
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : David Baur
// Creation Date : 2/4/2014
//-----------------------------------------------------------------------------
/// Apply a device instance "op" to all instances associated with this
/// model
/// 
/// @param[in] op Operator to apply to all instances.
/// 
/// 
void Model::forEachInstance(DeviceInstanceOp &op) const /* override */ 
{
  for (std::vector<Instance *>::const_iterator it = instanceContainer.begin(); it != instanceContainer.end(); ++it)
    op(*it);
}


Device *Traits::factory(const Configuration &configuration, const FactoryBlock &factory_block)
{

  return new DeviceMaster<Traits>(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);
}

//-----------------------------------------------------------------------------
// Function      : Instance::getDeviceLetter ()
//
// Purpose       : Returns first letter of device name string (basically U or Y).
//                 Returns blank ("") if the function fails.
//
// Special Notes : 
//
// Scope         : public
// Creator       : Pete Sholander
// Creation Date : 03/11/14
//-----------------------------------------------------------------------------
std::string Instance::getDeviceLetter ()
{
  int p1 = getName().find_first_of('%');
  std::string dev_letter = (p1 != std::string::npos) ? getName().substr(p1-1,1) : "";
  return dev_letter;
}

void registerDevice()
{
  // NOT device is deprecated now
  // removed .registerDevice("dltch", 1) since that device is still 
  // under development
  Config<Traits>::addConfiguration()
    .registerDevice("inv", 1)
    .registerDevice("not",1)
    .registerDevice("and", 1)
    .registerDevice("nand", 1)
    .registerDevice("or", 1)
    .registerDevice("nor", 1)
    .registerDevice("add", 1)
    .registerDevice("xor", 1)
    .registerDevice("nxor", 1)
    .registerDevice("dff", 1)
    .registerDevice("buf", 1)
    .registerModelType("dig", 1);
}

} // namespace Digital
} // namespace Device
} // namespace Xyce

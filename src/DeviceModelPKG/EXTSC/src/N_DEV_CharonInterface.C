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
// Filename       : $RCSfile: N_DEV_CharonInterface.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/13/06
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.25 $
//
// Revision Date  : $Date: 2014/02/24 23:49:16 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------
#include <iostream>

#include "Teuchos_ParameterList.hpp"

// ----------   Xyce Includes   ----------
#include <N_DEV_CharonInterface.h>
#include <N_CIR_Xyce.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_SolverState.h>
#include <N_ERH_ErrorMgr.h>
#include <N_UTL_BreakPoint.h>

#include <N_TIA_TimeIntInfo.h>
#include <N_TIA_TwoLevelError.h>

// RPP: We have added a copy of a charon interface file to xyce to
// eliminate a circular dependency in nevada's build TPL process.
// This file must be kept up-to-date with Charon's otherwise they will
// not be able to link with each other.
#ifdef Xyce_CHARON
#include "Charon_CircuitInterface.h"
#endif

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : CharonInterface::CharonInterface
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/13/06
//-----------------------------------------------------------------------------
CharonInterface::CharonInterface(
  const DeviceOptions & do1,
  const std::string &   netlist,
  const SolverState &   ss1)
  : devOptions_(do1),
    inputFileName_(netlist),
    solState_(ss1)
{

}

//-----------------------------------------------------------------------------
// Function      : CharonInterface::CharonInterface
// Purpose       : copy constructor
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/13/06
//-----------------------------------------------------------------------------
CharonInterface::CharonInterface (const CharonInterface &right)
  : devOptions_(right.devOptions_),
    inputFileName_(right.inputFileName_),
    solState_(right.solState_)
{

}

//-----------------------------------------------------------------------------
// Function      : CharonInterface::~CharonInterface
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/13/06
//-----------------------------------------------------------------------------
CharonInterface::~CharonInterface()

{

}

//-----------------------------------------------------------------------------
// Function      : CharonInterface::initialize
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/13/2006
//-----------------------------------------------------------------------------
bool CharonInterface::initialize(N_PDS_Comm * comm)
{

#ifdef Xyce_DEBUG_DEVICE
  if (devOptions_.debugLevel > 0)
  {
    Xyce::dout() << "In CharonInterface::initialize" << std::endl;
  }
#endif

  input_list_ = Teuchos::rcp(new Teuchos::ParameterList);

  output_list_ = Teuchos::rcp(new Teuchos::ParameterList);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : CharonInterface::simulateStep
// Purpose       :
// Special Notes : The 'input vector' is a vector of voltages corresponding
//                 to the connected nodes.
//
//                 The 'output vector', is a vector of currents corresponding
//                 to the same nodes.
//
//                 The jacobian is mainly conductances - dI/dV.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/13/2006
//-----------------------------------------------------------------------------
bool CharonInterface::simulateStep
      ( const SolverState & solState,
        const std::map<std::string,double> & inputMap,
        std::vector<double> & outputVector,
        std::vector< std::vector<double> > & jacobian,
        N_TIA_TwoLevelError & tlError
      )
{

#ifdef Xyce_DEBUG_DEVICE
  if (devOptions_.debugLevel > 0)
  {
    Xyce::dout() << "In CharonInterface::simulateStep" << std::endl;
  }
#endif

  // RPP: need to get this from xyce time integrator...
  double currentTime = solState.currTime;
  double stepSize = solState.currTimeStep;

  input_list_->set("Current Time", currentTime);
  input_list_->set("Time Step Size", stepSize);
  input_list_->set("Time Step Type", "Backward Euler");

  // Tell charon whether we are in a steady state or transient mode
  if (solState_.dcopFlag)
    input_list_->set("Solve Type", "Steady State");
  else
    input_list_->set("Solve Type", "Transient");

#ifdef Xyce_CHARON

  charon::sc::CircuitInterface::getInstance().takeStep(input_list_,
						       inputMap,
						       output_list_,
						       outputVector,
						       jacobian);

#else
  std::string msg = "CharonInterface::simulateStep: Charon support has not been enabled.  Rebuild xyce with the flag --enable-charon.";
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL,msg);
  return true;
#endif

  // Get the output parameters
  if (output_list_->isParameter("Error Sum"))
    tlError.xErrorSum = output_list_->get<double>("Error Sum");
  else
    tlError.xErrorSum = 0.0;

  if (output_list_->isParameter("Inner Size"))
    tlError.innerSize = output_list_->get<double>("Number of DOF");
  else
    tlError.innerSize = 1.0;

  // Get the return status
  int return_status = -3;
  if (output_list_->isParameter("Charon Status"))
    return_status = output_list_->get<int>("Charon Status");
  if (return_status < 0)
    return false;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : CharonInterface::finalize
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/13/2006
//-----------------------------------------------------------------------------
bool CharonInterface::finalize ()
{
#ifdef Xyce_DEBUG_DEVICE
  if (devOptions_.debugLevel > 0)
  {
    Xyce::dout() << "In CharonInterface::finalize" << std::endl;
  }
#endif


  return true;
}

//-----------------------------------------------------------------------------
// Function      : CharonInterface::run
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/13/06
//-----------------------------------------------------------------------------
bool CharonInterface::run ()
{
#ifdef Xyce_DEBUG_DEVICE
  if (devOptions_.debugLevel > 0)
  {
    Xyce::dout() << "In CharonInterface::run" << std::endl;
  }
#endif

  return true;
}

//-----------------------------------------------------------------------------
// Function      : CharonInterface::homotopyStepSuccess
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
void CharonInterface::homotopyStepSuccess
      (const std::vector<std::string> & paramNames,
       const std::vector<double> & paramVals)
{

  return;
}

//-----------------------------------------------------------------------------
// Function      : CharonInterface::homotopyStepFailure
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/30/06
//-----------------------------------------------------------------------------
void CharonInterface::homotopyStepFailure ()
{

  return;
}

//-----------------------------------------------------------------------------
// Function      : CharonInterface::stepSuccess
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/13/06
//-----------------------------------------------------------------------------
void CharonInterface::stepSuccess (int analysis)
{

#ifdef Xyce_CHARON

  //Xyce::dout() << "ROGER ** in stepSuccess!!!!" << std::endl;
  bool is_active = true;
  charon::sc::CircuitInterface::getInstance().acceptTimeStep(is_active);

#else
  std::string msg = "CharonInterface::stepSuccess: Charon support has not been enabled.  Rebuild xyce with the flag --enable-charon.";
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL,msg);
#endif

  return;
}

//-----------------------------------------------------------------------------
// Function      : CharonInterface::stepFailure
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/13/06
//-----------------------------------------------------------------------------
void CharonInterface::stepFailure (int analysis)
{

  return;
}


//-----------------------------------------------------------------------------
// Function      : CharonInterface::getInitialQnorm
// Purpose       :
// Special Notes : no-op.
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/15/06
//-----------------------------------------------------------------------------
bool CharonInterface::getInitialQnorm (N_TIA_TwoLevelError & tle)
{
  return true;
}

} // namespace Device
} // namespace Xyce

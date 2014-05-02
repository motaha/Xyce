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
// Filename       : $RCSfile: N_NLS_NOX_Interface.C,v $
//
// Purpose        : Interface to Xyce vectors for NOX.
//
// Special Notes  :
//
// Creator        : Tammy Kolda, NLS, 8950
//
// Creation Date  : 01/31/02
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.158 $
//
// Revision Date  : $Date: 2014/02/24 23:49:25 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#ifndef HAVE_CONFIG_H
#define HAVE_CONFIG_H
#endif

#include <sstream>

#include "N_NLS_NOX_Interface.h"
#include "N_NLS_NOX_SharedSystem.h"
#include "N_NLS_NOX_Group.h"
#include "N_NLS_LOCA_Group.h"
#include "LOCA_GlobalData.H"
#include "LOCA_StatusTest_Wrapper.H"
#include "NOX_Solver_Factory.H"

#include "N_NLS_NOX_PseudoTransientTest.h"
#include "N_NLS_NOX_PseudoTransientSolver.h"

#include "N_NLS_NOX_XyceTests.h"
#include "N_NLS_NOX_AugmentLinSys.h"
#include "N_LAS_Solver.h"
#include "N_LAS_Builder.h"
#include "N_LAS_System.h"
#include "N_LAS_QueryUtil.h"
#include "N_LOA_Loader.h"
#include "N_IO_CmdParse.h"
#include "N_IO_OutputMgr.h"
#include "N_ANP_AnalysisInterface.h"

#ifdef Xyce_PARALLEL_MPI
#include <N_PDS_ParComm.h>
#else
#include <N_PDS_SerialComm.h>
#endif


// ----------   NOX Includes   ----------
#include "LOCA.H"

// -----------  Forward declarations  -------
using namespace N_NLS_NOX;

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::Interface::Interface
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
Interface::Interface(N_IO_CmdParse & cp) :
  N_NLS_NonLinearSolver(cp),
  dcParams_(DC_OP),
  transientParams_(TRANSIENT),
  hbParams_(HB_MODE),
  sharedSystemPtr_(0),
  mode_(DC_OP),
  lastParametersMode_(DC_OP),
  parametersMode_(DC_OP),
  usemode_(true),
  copiedGroupFlag_(false),
  DCOPused_(false),
  ICspecified_(false),
  NODESETspecified_(false),
  isFirstContinuationParam_(true),
  firstSolveComplete_(false),
  iParam_(0)
{
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::Interface::~Interface
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
Interface::~Interface()
{
  delete sharedSystemPtr_;
  if (!Teuchos::is_null(globalDataPtr_))
  {
    LOCA::destroyGlobalData(globalDataPtr_);
  }
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::Interface::setOptions
// Purpose       : Passes option block corresponding to "NONLIN" onto
//                 nonlinear solver. These parameters set convergence
//                 tolerances, the type of method used, and so on.
// Special Notes :
// Return Type   : boolean
//
// See also      : setTranOptions, setAnalysisMode
//
// - Input Arguments -
//
//    OB         : Option block containing options corresponding to
//                 "NONLIN" in the netlist.
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool Interface::setOptions(const N_UTL_OptionBlock& OB)
{
  return dcParams_.setOptions(OB);
}


//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::Interface::setTranOptions
//
// Purpose       : Passes option block corresponding to "NONLIN-TRAN" onto
//                 nonlinear solver. This affects the settings used when
//                 the mode is Transient.
//
// Special Notes :
// Return Type   : boolean
//
// See also      : setOptions, setAnalysisMode
//
// - Input Arguments -
//    OB         : Option block containing options corresponding to
//                 "NONLIN-TRAN" in the netlist.
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool Interface::setTranOptions(const N_UTL_OptionBlock& OB)
{
  return transientParams_.setOptions(OB);
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::Interface::setHBOptions
// Purpose       :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool Interface::setHBOptions(const N_UTL_OptionBlock& OB)
{
  return hbParams_.setOptions(OB);
}


//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::Interface::setLocaOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool Interface::setLocaOptions(const N_UTL_OptionBlock& OB)
{
  return dcParams_.setLocaOptions(OB);
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::Interface::setDCOPRestartOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool Interface::setDCOPRestartOptions(const N_UTL_OptionBlock& OB)
{
  DCOPspecified_ = true;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::Interface::setICOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool Interface::setICOptions(const N_UTL_OptionBlock& OB)
{
  ICspecified_ = true;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::Interface::setNodeSetOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool Interface::setNodeSetOptions(const N_UTL_OptionBlock& OB)
{
  NODESETspecified_ = true;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::Interface::initializeAll
//
// Purpose       : Called after all register and set functions.
//                 Once the various registrations have taken place,
//                 this function sets the remaining pointers.
//
// Special Notes :
// Derived Notes : This function also calls the base object initializeAll.
//
// Special Notes:  This function obtains the solution, temporary solution and
//                 f vectors from the LAS system class.
//
// Return Type   : boolean
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool Interface::initializeAll()
{
  // Use the base class initialization
  bool initok = N_NLS_NonLinearSolver::initializeAll();
  if (!initok)
  {
    return false;
  }

  // Set the processor ID for printing in the nonlinear solver
  // For now, processor 0 will output information to the screen
  //int myPID = ((**nextSolVectorPtrPtr_).epetraObj()).Comm().MyPID();
  int myPID = pdsMgrPtr_->getPDSComm()->procID();
  dcParams_.setOutputOptions(myPID, 0);
  transientParams_.setOutputOptions(myPID, 0);

  // get var type list
  std::vector<char> varTypeVec;
#if 0
  topologyRcp_->returnVarTypeVec( varTypeVec );
#endif
  hbParams_.setOutputOptions(myPID, 0);

  // Set up the status tests
#ifdef Xyce_NLS_MASKED_WRMS_NORMS
  bool testsok =
    dcParams_.createStatusTests(currSolVectorPtrPtr_, getLoader(), varTypeVec,
      lasSysPtr_->getNonTrivialDeviceMaskFlag(), lasSysPtr_->getDeviceMaskVector()) &&
    transientParams_.createStatusTests(currSolVectorPtrPtr_, getLoader(), varTypeVec,
      lasSysPtr_->getNonTrivialDeviceMaskFlag(), lasSysPtr_->getDeviceMaskVector()) &&
    hbParams_.createStatusTests(currSolVectorPtrPtr_, getLoader(),  varTypeVec,
      lasSysPtr_->getNonTrivialDeviceMaskFlag(), lasSysPtr_->getDeviceMaskVector());
#else
  bool testsok =
    dcParams_.createStatusTests(currSolVectorPtrPtr_, getLoader(), varTypeVec) &&
    transientParams_.createStatusTests(currSolVectorPtrPtr_, getLoader(), varTypeVec);
    hbParams_.createStatusTests(currSolVectorPtrPtr_, getLoader(), varTypeVec);
#endif

  if (!testsok)
    return false;

  // Set up any linear solver options
  // We only set the tolerance if adaptive forcing is being used.
  setAZ_Tol_DC = false;
  setAZ_Tol_Transient = false;
  if (!(dcParams_.getNoxParams()->sublist("Direction").sublist("Newton")
	  .get("Forcing Term Method","Constant") == std::string("Constant")))
  {
    setAZ_Tol_DC = true;
    //Xyce::dout() << "NOX is overriding DC LS tol." << std::endl;
  }
  if (!(transientParams_.getNoxParams()->sublist("Direction").sublist("Newton")
	  .get("Forcing Term Method", "Constant") == std::string("Constant")))
  {
    setAZ_Tol_Transient = true;
    //Xyce::dout() << "NOX is overriding Transient LS tol." << std::endl;
  }
  if (!(hbParams_.getNoxParams()->sublist("Direction").sublist("Newton")
	  .get("Forcing Term Method", "Constant") == std::string("Constant")))
  {
    setAZ_Tol_Transient = true;
    //Xyce::dout() << "NOX is overriding Transient LS tol." << std::endl;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::Interface::solve
//
// Purpose       : Reset all the counters and parameters and solve the
//                 nonlinear problem for this time step. The solution is
//                 stored in nextSolVector (obtained from the N_LAS_System
//                 registered above by registerLinearSystem).
//
// Special Notes : Should not be called until *after* initializeAll() has
//                 been called.
//
// Return Type   : Integer - postive for sucess, negative for failure.
// Special Notes :
// Scope         : public
// Creator       : Roger Pawlowski, SNL 9233
// Creation Date :
//-----------------------------------------------------------------------------
int Interface::solve (N_NLS_NonLinearSolver * nlsTmpPtr)
{

 try
 {

  // For base object
  N_NLS_NonLinearSolver::resetCountersAndTimers_();

#ifdef Xyce_DEBUG_NONLINEAR
  N_NLS_NonLinearSolver::setDebugFlags ();
#endif

  // Setup the status tests
  if (Teuchos::is_null(locaStatusTestPtr_))
  {
    locaDCOpStatusTestPtr_ =
      Teuchos::rcp(new LOCA::StatusTest::Wrapper(dcParams_.getStatusTests()));
    locaTransientStatusTestPtr_ =
      Teuchos::rcp(new LOCA::StatusTest::Wrapper(transientParams_.getStatusTests()));
    locaHBStatusTestPtr_ =
      Teuchos::rcp(new LOCA::StatusTest::Wrapper(hbParams_.getStatusTests()));

  }

  // Pick the parameter set to use.
  ParameterSet* paramsPtr;
  if ((usemode_) && (mode_ == TRANSIENT))
  {
    paramsPtr = &transientParams_;
    locaStatusTestPtr_ = locaTransientStatusTestPtr_;
    lastParametersMode_ = parametersMode_;
    parametersMode_ = TRANSIENT;
  }
  else if ((usemode_) && (mode_ == HB_MODE))
  {
    paramsPtr = &hbParams_;
    locaStatusTestPtr_ = locaHBStatusTestPtr_;
    lastParametersMode_ = parametersMode_;
    parametersMode_ = HB_MODE;
  }
  else
  {
    paramsPtr = &dcParams_;
    locaStatusTestPtr_ = locaDCOpStatusTestPtr_;
    lastParametersMode_ = parametersMode_;
    parametersMode_ = DC_OP;
  }

  if (Teuchos::is_null(globalDataPtr_))
    globalDataPtr_ = LOCA::createGlobalData(paramsPtr->getAllParams());

  // set the xyce return codes:
  paramsPtr->setStatusTestReturnCodes(retCodes_);

  // Set up the shared system (we have to redo this every time because
  // the object pointed to by nextSolVectorPtrPtr may have changed.
  //delete sharedSystemPtr_;
  if (sharedSystemPtr_ == 0)
    sharedSystemPtr_ = new SharedSystem(**nextSolVectorPtrPtr_,
					*rhsVectorPtr_,
					*jacobianMatrixPtr_,
					*NewtonVectorPtr_,
					*gradVectorPtr_,
					*lasSysPtr_,
					*this);
  else
    sharedSystemPtr_->reset(**nextSolVectorPtrPtr_,
			    *rhsVectorPtr_,
			    *jacobianMatrixPtr_,
			    *NewtonVectorPtr_,
			    *gradVectorPtr_,
			    *lasSysPtr_,
			    *this);

  //////////////////////////////////////////////////////////////////////////
  // erkeite: Group handling required by 2-level Newton:
  // Reset up the corresponding group as well
  if (nlsTmpPtr==0)
  {
    if (Teuchos::is_null(groupPtr_))
    {
       groupPtr_ = Teuchos::rcp(new N_NLS_LOCA::Group(globalDataPtr_,
						      *sharedSystemPtr_,
						      getLoader(),
						      *outMgrPtr_,
						      *anaIntPtr_)
				);
    }
    else
    {
      N_NLS_NOX::Vector tmpVec(**nextSolVectorPtrPtr_, *lasSysPtr_);
      groupPtr_->setX(tmpVec);
    }
  }
  else
  {
    copiedGroupFlag_ = true;
    Interface * nlsOtherPtr = dynamic_cast<Interface*>(nlsTmpPtr);
    groupPtr_ = nlsOtherPtr->getSolutionGroup();
  }
  // End of block needed by 2-level Newton.
  //////////////////////////////////////////////////////////////////////////

  int solverType = paramsPtr->getNoxSolverType();

  // Setting the nonContinuation flag is required to prevent incorrect
  // N_DEV_DeviceMgr::setParam calls.
  if (solverType==0)
  {
    groupPtr_->setNonContinuationFlag (true);
  }
  else
  {
    groupPtr_->setNonContinuationFlag (false);
  }

#ifdef Xyce_DEBUG_NONLINEAR
  Xyce::dout() << "solverType is " << solverType << std::endl;
#endif

  // (0) Standard Newton Method Solve (includes line search and
  // trust region based methods
  if (solverType == 0)
  {
    bool usedIC=false;
    bool usedNODESET=false;
    bool usedOP=false;

    // RPP: needed for tensor method.  DO NOT UNCOMMENT!
    // This somehow breaks the one-shot test circuit.
    // ERK:  it breaks it (and probably many other circuits) b/c voltage
    //       limiting is fragile and cannot handle extra (unexpected)
    //       F-loads.
    //groupPtr_->computeF();

    // Set up nox nonlinear solver.  solverPtr is only needed for
    // non-LOCA (i.e. just Newton's method, no continuation) solves.
    if (Teuchos::is_null(solverPtr_))
    {
#ifdef Xyce_DEBUG_NONLINEAR
      Xyce::dout() << "Creating new NLS solver b/c it is 0." <<std::endl;
#endif
      solverPtr_ = NOX::Solver::buildSolver(groupPtr_,
                                            paramsPtr->getStatusTests(),
                                            paramsPtr->getNoxParams());

      // If .IC, .NODESET or .OP have been specified, then set up the
      // augmented linear systems
      if ((usemode_) && (mode_ != TRANSIENT))
      {
        if (ICspecified_)
        {
          usedIC=icCont (paramsPtr);
        }
        else if (NODESETspecified_)
        {
          usedNODESET=nodesetCont0 (paramsPtr);
        }
        else
        {
          usedOP = opStartCont0 (paramsPtr);
        }
      }
    }
    else if ((usemode_) && (lastParametersMode_ != parametersMode_))
    {
#ifdef Xyce_DEBUG_NONLINEAR
      Xyce::dout() << "Creating new NLS solver b/c starting next phase, post-DC." <<std::endl;
#endif
      // Create a new solver for the next phase of the simulation.
      solverPtr_ = NOX::Solver::buildSolver(groupPtr_,
                                            paramsPtr->getStatusTests(),
                                            paramsPtr->getNoxParams());
    }
    else  // solverPtr is not null, and the mode didn't change since last call.
    {
#ifdef Xyce_DEBUG_NONLINEAR
      Xyce::dout() << "NOT Creating new NLS solver, just resetting." <<std::endl;
#endif
      solverPtr_->reset(groupPtr_->getX());
    }

    NOX::StatusTest::StatusType status = solverPtr_->solve();
    firstSolveComplete_ = true;

    if (usedIC)
    {
      Teuchos::RefCountPtr<N_NLS_NOX::AugmentLinSys> alsNull; // 0
      groupPtr_->setAugmentLinearSystem(false, alsNull);
    }
    // Send back the correct return code
#ifdef Xyce_DEBUG_NONLINEAR
    Xyce::dout() << "return code for transient params: " << transientParams_.getStatusTestReturnCode() <<std::endl;;
    Xyce::dout() << "return code for hb params: " << hbParams_.getStatusTestReturnCode() << std::endl;
    Xyce::dout() << "return code for dc params: " << dcParams_.getStatusTestReturnCode () <<std::endl;
#endif
    return (paramsPtr->getStatusTestReturnCode());
  }
  // (1) Natural Parameter Continuation
  else if (solverType == 1)
  {
    std::vector<std::string> pars;
    std::string con;
    int j, numParam = 0;
    double value;
    bool usedOP=false;
    bool usedNODESET=false;
    bool usedIC=false;

    if ((usemode_) && (mode_ != TRANSIENT))
    {
      if (ICspecified_)
      {
        usedIC=icCont (paramsPtr);
      }
      else if (NODESETspecified_)
      {
        usedNODESET=nodesetCont1 (paramsPtr);
      }
      else
      {
        usedOP = opStartCont1 (paramsPtr);
      }
    }

    Teuchos::RefCountPtr<Teuchos::ParameterList> locaList =
      paramsPtr->getLocaParams();

    // Create Parameter Vector and get the stepper parameter list.
    LOCA::ParameterVector locaPVec;
    Teuchos::ParameterList& stepperList = locaList->sublist("Stepper");
    Teuchos::ParameterList& stepSizeList = locaList->sublist("Step Size");

    while (paramsPtr->getVectorParam("CONPARAM", numParam, con))
    {
      pars.push_back(con);
      numParam++;
    }

    // Fail out if no parameters have been specified.
    if ( numParam == 0 ) {
      std::string message = "Using \"continuation=1\" requires a parameter to be set with the conparam keyword in the loca option block!";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL, message);
    }

    Teuchos::RefCountPtr<N_NLS_NOX::AugmentLinSys> als;
    Teuchos::RefCountPtr<N_NLS_NOX::AugmentLinSys> alsNull; // 0

    std::vector<double> minValue(numParam, 0.);
    std::vector<double> maxValue(numParam, 0.);
    std::vector<double> initialValue(numParam, 0.);
    std::vector<double> initialStepSize(numParam, 0.);
    std::vector<double> minStepSize(numParam, 0.);
    std::vector<double> maxStepSize(numParam, 0.);
    std::vector<double> aggressiveness(numParam, 0.);

    // Check the size of params to make sure users provided all required data
    std::vector<std::string> paramNames(0);
    paramNames.push_back("MINVALUE");
    paramNames.push_back("MAXVALUE");
    paramNames.push_back("INITIALVALUE");
    paramNames.push_back("INITIALSTEPSIZE");
    paramNames.push_back("MINSTEPSIZE");
    paramNames.push_back("MAXSTEPSIZE");
    paramNames.push_back("AGGRESSIVENESS");

    for (std::size_t p = 0 ; p < paramNames.size() ; ++p) {
      if ( paramsPtr->getVectorParamSize(paramNames[p]) != numParam) {
        std::string msg = "The parameter \"" +
          paramNames[p] +
          "\" must have a list of values with size equal to the numParamber of parameters specified in \"conparam\".";
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
      }
    }

    for (iParam_=0 ; iParam_<numParam ; ++iParam_)
    {
      paramsPtr->getVectorParam("MINVALUE", iParam_, value);
      minValue[iParam_] = value;

      paramsPtr->getVectorParam("MAXVALUE", iParam_, value);
      maxValue[iParam_] = value;

      paramsPtr->getVectorParam("INITIALVALUE", iParam_, value);
      initialValue[iParam_] = value;

      paramsPtr->getVectorParam("INITIALSTEPSIZE", iParam_, value);
      initialStepSize[iParam_] = value;

      paramsPtr->getVectorParam("MINSTEPSIZE", iParam_, value);
      minStepSize[iParam_] = value;

      paramsPtr->getVectorParam("MAXSTEPSIZE", iParam_, value);
      maxStepSize[iParam_] = value;

      paramsPtr->getVectorParam("AGGRESSIVENESS", iParam_, value);
      aggressiveness[iParam_] = value;

      locaPVec.addParameter (pars[iParam_], initialValue[iParam_]);

    }

    if (usedOP || usedNODESET)
    {
      const N_NLS_LOCA::Group & conLocaGrp =
        dynamic_cast<const N_NLS_LOCA::Group&>(solverPtr_->getSolutionGroup());
      *groupPtr_ = const_cast<N_NLS_LOCA::Group&>(conLocaGrp);
      solverPtr_ = Teuchos::null;
    }

    groupPtr_->setParams(locaPVec);

    LOCA::Abstract::Iterator::IteratorStatus locaStatus;

    for (iParam_=0 ; iParam_<numParam ; ++iParam_)
    {
      // Copy out the solution and use it in the next run
      if (iParam_ > 0)
      {
        groupPtr_->copy(*(stepperPtr_->getSolutionGroup()));
      }

      stepperList.set("Continuation Parameter", pars[iParam_]);
      stepperList.set("Initial Value", initialValue[iParam_]);
      stepperList.set("Max Value", maxValue[iParam_]);
      stepperList.set("Min Value", minValue[iParam_]);
      stepSizeList.set("Initial Step Size", initialStepSize[iParam_]);
      stepSizeList.set("Min Step Size", minStepSize[iParam_]);
      stepSizeList.set("Max Step Size", maxStepSize[iParam_]);
      stepSizeList.set("Aggressiveness", aggressiveness[iParam_]);

      for (j=0 ; j<iParam_ ; ++j)
      {
        locaPVec.setValue(pars[j], maxValue[j]);
      }
      for (j=iParam_ ; j<numParam ; ++j)
      {
        locaPVec.setValue(pars[j], initialValue[j]);
      }
      groupPtr_->setParams(locaPVec);

      if (iParam_==0)
      {
        if (!usedOP && !usedNODESET) // (usedOP and usedNODESET have already loaded F)
        {
          groupPtr_->computeF();
        }
      }

      // RPP 03/08/2006 If this is a special parameter used in voltage
      // node resistance (Spice's GMIN stepping) then we need to
      // intercept and handle the parameter in the solver (in the LOCA
      // Group) - the device package does nothing for this parameter.
      if (pars[iParam_] != "GSTEPPING")
      {
        // Do the continuation run
        resetStepper(globalDataPtr_, groupPtr_, locaStatusTestPtr_, paramsPtr->getAllParams());
        locaStatus = stepperPtr_->run();
        if (usedIC || usedNODESET)
        {
          groupPtr_->setAugmentLinearSystem(false, alsNull);
        }
      }
      else
      {
        if (iParam_ == 0 && DCOPused_)
        {
          std::string message = "'.dcop input=' and gstepping are incompatible";
          N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL, message);
        }
        Teuchos::RefCountPtr<N_NLS_NOX::AugmentLinSys> als =
          paramsPtr->createAugmentLinearSystem(lasSysPtr_);
        groupPtr_->setAugmentLinearSystem(true, als);

        // Do the continuation run
        resetStepper(globalDataPtr_, groupPtr_, locaStatusTestPtr_, paramsPtr->getAllParams());
        locaStatus = stepperPtr_->run();

        groupPtr_->setAugmentLinearSystem(false, alsNull);
      }

      // Kick out if continuation failed
      if (locaStatus != LOCA::Abstract::Iterator::Finished)
        return (-1);

      // Increment Param Number Tracking
      isFirstContinuationParam_ = false;
      firstSolveComplete_ = true;
    }
    return (paramsPtr->getStatusTestReturnCode());
  }
  // (2) Mosfet specific continuation
  else if (solverType == 2)
  {
    bool usedOP=false;
    bool usedNODESET=false;
    bool usedIC=false;
    Teuchos::RefCountPtr<N_NLS_NOX::AugmentLinSys> alsNull; // 0

    if ((usemode_) && (mode_ != TRANSIENT))
    {
      if (ICspecified_)
      {
        usedIC=icCont (paramsPtr);
      }
      else if (NODESETspecified_)
      {
        usedNODESET=nodesetCont1 (paramsPtr);
      }
      else
      {
        usedOP = opStartCont1 (paramsPtr);
      }
    }

    // check if use specified vector params.  If so, flag a warning.
    std::vector<std::string> pars;
    std::string con;
    int numParam;
    while (paramsPtr->getVectorParam("CONPARAM", numParam, con))
    {
      pars.push_back(con);
      numParam++;
    }

    if ( numParam > 1 ) {
      std::string message = "WARNING: Using \"continuation=2\" currently does not support vectorized loca parameters.  Only the first values in each list will be used.";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_WARNING, message);
    }

    Teuchos::RefCountPtr<Teuchos::ParameterList> locaList = paramsPtr->getLocaParams();

    // Create the continuation parameter names
    std::string gain = "mosfet:gainscale";
    std::string nonlinear = "mosfet:nltermscale";

    // Create Parameter Vector and get the stepper parameter list.
    LOCA::ParameterVector locaPVec;
    Teuchos::ParameterList& stepperList = locaList->sublist("Stepper");
    Teuchos::ParameterList& predictorList = locaList->sublist("Predictor");
    Teuchos::ParameterList& stepSizeList = locaList->sublist("Step Size");

    // Continuation solve from alpha2 (gain) = 0.0 -> 1.0
    //                    with alpha1 (nlterm) = 0.0 (constant)
    locaPVec.addParameter(gain, 0.0);
    locaPVec.addParameter(nonlinear, 0.0);
    groupPtr_->setParams(locaPVec);
    stepperList.set("Continuation Parameter", gain);

    stepperList.set("Initial Value", 0.0);
    stepperList.set("Max Value", 1.0);
    stepperList.set("Min Value",-1.0);

    stepSizeList.set("Initial Step Size", 0.2);
    stepSizeList.set("Min Step Size", 1.0e-4);
    stepSizeList.set("Max Step Size", 1.0);
    stepSizeList.set("Aggressiveness", 1.0);

    // assert the user-specified defaults, if any.
    dcParams_.applySavedLocaOptions();

    // Initialize parameters in xyce
    if (!usedOP && !usedNODESET) // (usedOP and usedNODESET have already loaded F)
    {
      groupPtr_->computeF();
    }

    // Do the continuation run
    iParam_ = 0;
    resetStepper(globalDataPtr_, groupPtr_, locaStatusTestPtr_, paramsPtr->getAllParams());
    LOCA::Abstract::Iterator::IteratorStatus locaStatus = stepperPtr_->run();

    // Kick out if continuation failed
    if (locaStatus != LOCA::Abstract::Iterator::Finished)
      return (-1);

    // Increment Param Number Tracking
    isFirstContinuationParam_ = false;
    firstSolveComplete_ = true;

    // Copy out the solution and use it in the next run
    groupPtr_->copy(*(stepperPtr_->getSolutionGroup()));

    // Continuation solve from alpha1 (nlterm) = 0.0 -> 1.0
    //                    with alpha2 (gain) = 1.0 (constant)
    stepperList.set("Continuation Parameter", nonlinear);

    stepperList.set("Initial Value", 0.0);
    stepperList.set("Max Value", 1.0);
    stepperList.set("Min Value",-1.0);

    stepSizeList.set("Initial Step Size", 0.2);
    stepSizeList.set("Min Step Size", 1.0e-4);
    stepSizeList.set("Max Step Size", 1.0);
    stepSizeList.set("Aggressiveness", 1.0);

    // assert the user-specified defaults, if any.
    dcParams_.applySavedLocaOptions();

    locaPVec.setValue(gain, 1.0);
    locaPVec.setValue(nonlinear, 0.0);
    groupPtr_->setParams(locaPVec);

    // Do the continuation run
    iParam_ = 1;
    resetStepper(globalDataPtr_, groupPtr_, locaStatusTestPtr_, paramsPtr->getAllParams());
    locaStatus = stepperPtr_->run();

    if (usedIC)
    {
      groupPtr_->setAugmentLinearSystem(false, alsNull);
    }

    // Return the solution status
    if (locaStatus != LOCA::Abstract::Iterator::Finished)
      return (-1);
    else
      return (paramsPtr->getStatusTestReturnCode());

  }
  else if (solverType == 3)  // GMIN stepping, simple specification
  {


    Teuchos::RefCountPtr<Teuchos::ParameterList> locaList = paramsPtr->getLocaParams();

    // Create the continuation parameter names
    std::string gmin = "GSTEPPING";

    // Create Parameter Vector and get the stepper parameter list.
    LOCA::ParameterVector locaPVec;
    Teuchos::ParameterList& stepperList = locaList->sublist("Stepper");
    Teuchos::ParameterList& predictorList = locaList->sublist("Predictor");
    Teuchos::ParameterList& stepSizeList = locaList->sublist("Step Size");

    // Continuation solve using Gmin stepping.
    locaPVec.addParameter(gmin, 0.0);
    groupPtr_->setParams(locaPVec);
    stepperList.set("Continuation Parameter", gmin);
    stepperList.set("Continuation Method", "Natural");

    stepSizeList.set("Method", "Adaptive");
    predictorList.set("Method", "Constant");

    stepperList.set("Initial Value", 4.0);
    stepperList.set("Min Value", -4.0);
    paramsPtr->set_gstepping_min_value (-4.0);
    stepperList.set("Max Value", 4.0);

    stepSizeList.set("Initial Step Size", -2.0);
    stepSizeList.set("Min Step Size", 1.0e-6);
    stepSizeList.set("Max Step Size", 1.0e+12);
    stepSizeList.set("Aggressiveness", 0.01);

    stepperList.set("Max Steps", 400);
    stepperList.set("Max Nonlinear Iterations", 20);

    // the following set of if-statemtents call functions that
    // allocate augmented linear systems for various scenarios.
    // It is important that the augmented systems get allocated
    // after the paramter (above) have been set.
    bool usedOP=false;
    bool usedNODESET=false;
    bool usedIC=false;
    Teuchos::RefCountPtr<N_NLS_NOX::AugmentLinSys> alsNull; // 0

    if ((usemode_) && (mode_ != TRANSIENT))
    {
      if (ICspecified_)
      {
        usedIC=icCont3 (paramsPtr);
      }
      else if (NODESETspecified_)
      {
        usedNODESET=nodesetCont1 (paramsPtr);
      }
      else
      {
        usedOP = opStartCont1 (paramsPtr);
      }
    }

    // Initialize parameters in xyce
    if (!usedOP && !usedNODESET) // (usedOP and usedNODESET have already loaded F)
    {
      groupPtr_->computeF();
    }

    // Do the continuation run
    iParam_ = 0;
    if (iParam_ == 0 && DCOPused_)
    {
      std::string message = "'.dcop input=' and gstepping are incompatible";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL, message);
    }
    if (!usedIC)
    {
      Teuchos::RefCountPtr<N_NLS_NOX::AugmentLinSys> als =
        paramsPtr->createAugmentLinearSystem(lasSysPtr_);
      groupPtr_->setAugmentLinearSystem(true, als);
    }

    // Do the continuation run
    resetStepper(globalDataPtr_, groupPtr_, locaStatusTestPtr_, paramsPtr->getAllParams());
    LOCA::Abstract::Iterator::IteratorStatus locaStatus = stepperPtr_->run();

    groupPtr_->setAugmentLinearSystem(false, alsNull);

    // Kick out if continuation failed
    if (locaStatus != LOCA::Abstract::Iterator::Finished)
      return (-1);
    else
      return (paramsPtr->getStatusTestReturnCode());
  }
  // (4) Mosfet specific continuation (ERK, new 2/21/2004)
  else if (solverType == 4)
  {

    Teuchos::RefCountPtr<Teuchos::ParameterList> locaList = paramsPtr->getLocaParams();

    // Create the continuation parameter names
    std::string gain = "mosfet:gainscale";
    std::string nonlinear = "mosfet:nltermscale";
    std::string size = "mosfet:sizescale";

    // Create Parameter Vector and get the stepper parameter list.
    LOCA::ParameterVector locaPVec;
    Teuchos::ParameterList& stepperList = locaList->sublist("Stepper");

    // Continuation solve from alpha2 (gain) = 0.0 -> 1.0
    //                    with alpha1 (nlterm) = 0.0 (constant)
    //                    sizeScale = 0.0 (constant)
    locaPVec.addParameter(gain, 0.0);
    locaPVec.addParameter(nonlinear, 0.0);
    locaPVec.addParameter(size, 0.0);
    groupPtr_->setParams(locaPVec);
    stepperList.set("Continuation Parameter", gain);
    stepperList.set("Initial Value", 0.0);
    stepperList.set("Max Value", 1.0);

    // Initialize parameters in xyce
    groupPtr_->computeF();

    // Do the continuation run
    resetStepper(globalDataPtr_, groupPtr_, locaStatusTestPtr_, paramsPtr->getAllParams());
    LOCA::Abstract::Iterator::IteratorStatus locaStatus = stepperPtr_->run();

    // Kick out if continuation failed
    if (locaStatus != LOCA::Abstract::Iterator::Finished)
      return (-1);

    // Increment Param Number Tracking
    isFirstContinuationParam_ = false;
    firstSolveComplete_ = true;

    // Copy out the solution and use it in the next run
    groupPtr_->copy(*(stepperPtr_->getSolutionGroup()));

    // Continuation solve from alpha1 (nlterm) = 0.0 -> 1.0
    //                    with alpha2 (gain) = 1.0 (constant)
    //                    with size   (scale) = 0.0 (constant)
    stepperList.set("Continuation Parameter", nonlinear);
    stepperList.set("Initial Value", 0.0);
    stepperList.set("Max Value", 1.0);
    locaPVec.setValue(gain, 1.0);
    locaPVec.setValue(nonlinear, 0.0);
    locaPVec.setValue(size, 0.0);
    groupPtr_->setParams(locaPVec);

    // Do the continuation run
    resetStepper(globalDataPtr_, groupPtr_, locaStatusTestPtr_, paramsPtr->getAllParams());
    locaStatus = stepperPtr_->run();

    // Kick out if continuation failed
    if (locaStatus != LOCA::Abstract::Iterator::Finished)
      return (-1);

    // Copy out the solution and use it in the next run
    groupPtr_->copy(*(stepperPtr_->getSolutionGroup()));

    // Continuation solve from alpha1 (nlterm) = 1.0 -> 1.0
    //                    with alpha2 (gain) = 1.0 (constant)
    //                    with size   (scale) = 0.0 -> 1.0
    stepperList.set("Continuation Parameter", size);
    stepperList.set("Initial Value", 0.0);
    stepperList.set("Max Value", 1.0);
    locaPVec.setValue(gain, 1.0);
    locaPVec.setValue(nonlinear, 1.0);
    locaPVec.setValue(size, 0.0);
    groupPtr_->setParams(locaPVec);

    // Do the continuation run
    resetStepper(globalDataPtr_, groupPtr_, locaStatusTestPtr_, paramsPtr->getAllParams());
    locaStatus = stepperPtr_->run();

    // Return the solution status
    if (locaStatus != LOCA::Abstract::Iterator::Finished)
      return (-1);
    else
      return (paramsPtr->getStatusTestReturnCode());

  }
  // (5) Mosfet:BSIM3:Inverter specific continuation (RPP, new 2/25/2004)
  else if (solverType == 5)
  {

    Teuchos::RefCountPtr<Teuchos::ParameterList> locaList = paramsPtr->getLocaParams();

    // Create the continuation parameter names
    //std::string gain = "mosfet:gainscale";
    std::string nonlinear = "mosfet:nltermscale";
    std::string size = "mosfet:sizescale";

    // Create Parameter Vector and get the stepper parameter list.
    LOCA::ParameterVector locaPVec;
    Teuchos::ParameterList& stepperList = locaList->sublist("Stepper");

    // Continuation solve from alpha2 (gain) = 0.0 -> 1.0
    //                    with alpha1 (nlterm) = 0.0 (constant)
    //                    sizeScale = 0.0 (constant)
    //locaPVec.addParameter(gain, 0.0);
    locaPVec.addParameter(nonlinear, 0.0);
    locaPVec.addParameter(size, 0.0);
    groupPtr_->setParams(locaPVec);
    stepperList.set("Continuation Parameter", nonlinear);
    stepperList.set("Initial Value", 0.0);
    stepperList.set("Max Value", 1.0);

    // Initialize parameters in xyce
    groupPtr_->computeF();

    // Do the continuation run
    resetStepper(globalDataPtr_, groupPtr_, locaStatusTestPtr_, paramsPtr->getAllParams());
    LOCA::Abstract::Iterator::IteratorStatus locaStatus = stepperPtr_->run();

    // Kick out if continuation failed
    if (locaStatus != LOCA::Abstract::Iterator::Finished)
      return (-1);

    // Increment Param Number Tracking
    isFirstContinuationParam_ = false;
    firstSolveComplete_ = true;

    // Copy out the solution and use it in the next run
    groupPtr_->copy(*(stepperPtr_->getSolutionGroup()));

    // Continuation solve from alpha1 (nlterm) = 0.0 -> 1.0
    //                    with alpha2 (gain) = 1.0 (constant)
    //                    with size   (scale) = 0.0 (constant)
    stepperList.set("Continuation Parameter", size);
    stepperList.set("Initial Value", 0.0);
    stepperList.set("Max Value", 1.0);
    //locaPVec.setValue(gain, 1.0);
    locaPVec.setValue(nonlinear, 1.0);
    locaPVec.setValue(size, 0.0);
    groupPtr_->setParams(locaPVec);

    // Do the continuation run
    resetStepper(globalDataPtr_, groupPtr_, locaStatusTestPtr_, paramsPtr->getAllParams());
    locaStatus = stepperPtr_->run();

    // Return the solution status
    if (locaStatus != LOCA::Abstract::Iterator::Finished)
      return (-1);
    else
      return (paramsPtr->getStatusTestReturnCode());

  }
  // (6) Mosfet:BSIM3:Inverter specific continuation (RPP, new 2/25/2004)
  else if (solverType == 6)
  {

    Teuchos::RefCountPtr<Teuchos::ParameterList> locaList = paramsPtr->getLocaParams();

    // Create the continuation parameter names
    std::string gain = "mosfet:gainscale";
    std::string nonlinear = "mosfet:nltermscale";
    std::string size = "mosfet:sizescale";

    // Create Parameter Vector and get the stepper parameter list.
    LOCA::ParameterVector locaPVec;
    Teuchos::ParameterList& stepperList = locaList->sublist("Stepper");

    // Continuation solve
    locaPVec.addParameter(gain, 1.0);
    locaPVec.addParameter(nonlinear, 1.0);
    locaPVec.addParameter(size, 0.0);
    groupPtr_->setParams(locaPVec);
    stepperList.set("Continuation Parameter", size);
    stepperList.set("Initial Value", 0.0);
    stepperList.set("Max Value", 1.0);

    // Initialize parameters in xyce
    groupPtr_->computeF();

    // Do the continuation run
    resetStepper(globalDataPtr_, groupPtr_, locaStatusTestPtr_, paramsPtr->getAllParams());
    LOCA::Abstract::Iterator::IteratorStatus locaStatus = stepperPtr_->run();

    // Return the solution status
    if (locaStatus != LOCA::Abstract::Iterator::Finished)
      return (-1);
    else
    {
      firstSolveComplete_ = true;
      return (paramsPtr->getStatusTestReturnCode());
    }

  }  // Block gainscale
  else if (solverType == 7)
  {

    // Get some initial objects
    Teuchos::RefCountPtr<Teuchos::ParameterList> locaList = paramsPtr->getLocaParams();
    Teuchos::ParameterList& stepperList = locaList->sublist("Stepper");
    LOCA::ParameterVector locaPVec;

    // Create storage for continuation objects
    int numGainBlocks = loaderPtr_->getHomotopyBlockSize();
    int numHomotopyContinuationRuns = 1 + numGainBlocks;
    std::vector<std::string> names(numHomotopyContinuationRuns);
    std::vector<double> initialVal(numHomotopyContinuationRuns);
    std::vector<double> finalVal(numHomotopyContinuationRuns);
    std::vector<double> minVal(numHomotopyContinuationRuns);
    std::vector<double> maxVal(numHomotopyContinuationRuns);

    // Set up continuation steps
    // ***************************************************
    // Changes vals below for continuation
    // ***************************************************
    for (int i = 0; i < numGainBlocks; ++i) {
      std::stringstream s;
      s << i;
      names[i] = "mosfet:gainscale_block_" + s.str();
      initialVal[i] = 0.0;
      finalVal[i] = 1.0;
    }
    names[numHomotopyContinuationRuns - 1] = "mosfet:nltermscale";
    initialVal[numHomotopyContinuationRuns - 1] = 0.0;
    finalVal[numHomotopyContinuationRuns - 1] = 1.0;
    // ***************************************************
    // ***************************************************

    // Ste up max/min bounds
    for (int i = 0; i < names.size(); ++i) {
      if (finalVal[i] > initialVal[i]) {
        minVal[i] = initialVal[i];
        maxVal[i] = finalVal[i];
      }
      else {
        minVal[i] = finalVal[i];
        maxVal[i] = initialVal[i];
      }
    }

    // Initialize loca parameter vector
    for (int i = 0; i < names.size(); ++i)
      locaPVec.addParameter(names[i], initialVal[i]);

    LOCA::Abstract::Iterator::IteratorStatus locaStatus;

    // Loop over the number of homotopy steps
    for (int hs = 0; hs < names.size(); ++hs) {
      for (int i = 0; i < names.size(); ++i) {
	if (i >= hs)
	  locaPVec.setValue(names[i], initialVal[i]);
	else
	  locaPVec.setValue(names[i], finalVal[i]);
      }
      groupPtr_->setParams(locaPVec);
      stepperList.set("Continuation Parameter", names[hs]);
      stepperList.set("Initial Value", initialVal[hs]);
      stepperList.set("Min Value", minVal[hs]);
      stepperList.set("Max Value", maxVal[hs]);

      // Initialize parameters in xyce
      groupPtr_->computeF();

      // Do the continuation run
      resetStepper(globalDataPtr_, groupPtr_, locaStatusTestPtr_, paramsPtr->getAllParams());
      locaStatus = stepperPtr_->run();

      // Kick out if continuation failed
      if (locaStatus != LOCA::Abstract::Iterator::Finished)
        return (-1);

      // Increment Param Number Tracking
      isFirstContinuationParam_ = false;
      firstSolveComplete_ = true;

      // Copy out the solution and use it in the next run
      groupPtr_->copy(*(stepperPtr_->getSolutionGroup()));
    }

    // Return converged solver code
    return (paramsPtr->getStatusTestReturnCode());

  }  // Test suite
  else if (solverType == 8)
  {

    // Get some initial objects
    Teuchos::RefCountPtr<Teuchos::ParameterList> locaList = paramsPtr->getLocaParams();
    Teuchos::ParameterList& stepperList = locaList->sublist("Stepper");
    LOCA::ParameterVector locaPVec;

    // Create storage for continuation objects
    int numHomotopyContinuationRuns = 1;
    std::vector<std::string> names(numHomotopyContinuationRuns);
    std::vector<double> initialVal(numHomotopyContinuationRuns);
    std::vector<double> finalVal(numHomotopyContinuationRuns);
    std::vector<double> minVal(numHomotopyContinuationRuns);
    std::vector<double> maxVal(numHomotopyContinuationRuns);

    // Set up continuation steps
    // ***************************************************
    // Changes vals below for continuation
    // ***************************************************
    names[0] = "mosfet:gainscale";
    //names[1] = "mosfet:nltermscale";
    //names[2] = "mosfet:sizescale";
    //names[2] = "mosfet:l";
    initialVal[0] = 0.0;
    //initialVal[1] = 0.0;
    //initialVal[2] = 0.0;
    finalVal[0] = 1.0;
    //finalVal[1] = 1.0;
    //finalVal[2] = 1.0;
    // ***************************************************
    std::string n1 = "mosfet:nltermscale";
    locaPVec.addParameter(n1, 0.0);
    //std::string n2 = "mosfet:sizescale";
    //locaPVec.addParameter(n2, 0.0);
    // ***************************************************
    // ***************************************************

    // Ste up max/min bounds
    for (int i = 0; i < names.size(); ++i) {
      if (finalVal[i] > initialVal[i]) {
        minVal[i] = initialVal[i];
        maxVal[i] = finalVal[i];
      }
      else {
        minVal[i] = finalVal[i];
        maxVal[i] = initialVal[i];
      }
    }

    // Initialize loca parameter vector
    for (int i = 0; i < names.size(); ++i)
      locaPVec.addParameter(names[i], initialVal[i]);

    LOCA::Abstract::Iterator::IteratorStatus locaStatus;

    LOCA::StatusTest::Wrapper test(paramsPtr->getStatusTests());

    // Loop over the number of homotopy steps
    for (int hs = 0; hs < names.size(); ++hs) {
      for (int i = 0; i < names.size(); ++i) {
        if (i >= hs)
          locaPVec.setValue(names[i], initialVal[i]);
        else
          locaPVec.setValue(names[i], finalVal[i]);
      }
      groupPtr_->setParams(locaPVec);
      stepperList.set("Continuation Parameter", names[hs]);
      stepperList.set("Initial Value", initialVal[hs]);
      stepperList.set("Min Value", minVal[hs]);
      stepperList.set("Max Value", maxVal[hs]);

      // Initialize parameters in xyce
      groupPtr_->computeF();

      // Do the continuation run
      resetStepper(globalDataPtr_, groupPtr_, locaStatusTestPtr_, paramsPtr->getAllParams());
      locaStatus = stepperPtr_->run();

      // Kick out if continuation failed
      if (locaStatus != LOCA::Abstract::Iterator::Finished)
        return (-1);

      // Increment Param Number Tracking
      isFirstContinuationParam_ = false;
      firstSolveComplete_ = true;

      // Copy out the solution and use it in the next run
      groupPtr_->copy(*(stepperPtr_->getSolutionGroup()));
    }

    // Return converged solver code
    return (paramsPtr->getStatusTestReturnCode());

  } // Pseudo Transient
  else if (solverType == 9)
  {
    Teuchos::RefCountPtr<NOX::StatusTest::Combo> ctest =
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));

    Teuchos::RefCountPtr<Teuchos::ParameterList> locaList =
      paramsPtr->getLocaParams();
    Teuchos::ParameterList& stepperList = locaList->sublist("Stepper");
    Teuchos::ParameterList& stepSizeList = locaList->sublist("Step Size");
    double initialStepSize =
      stepSizeList.get("Initial Step Size", 1.0e-3);
    double minStepSize = stepSizeList.get("Min Step Size", 1.0e-12);
    double maxStepSize = stepSizeList.get("Max Step Size", 1.0e4);
    Teuchos::RefCountPtr<Teuchos::ParameterList> noxList =
      paramsPtr->getNoxParams();

    // Create Pseudo Transient status tests.
    //paramsPtr->getStatusTests()
    Teuchos::RefCountPtr<NOX::StatusTest::MaxIters> mi =
      Teuchos::rcp(new NOX::StatusTest::MaxIters(stepperList.get("Max Steps", 200)));
    Teuchos::RefCountPtr<NOX::StatusTest::FiniteValue> fv =
      Teuchos::rcp(new NOX::StatusTest::FiniteValue);
    Teuchos::RefCountPtr<N_NLS_NOX::PseudoTransientTest> pt =
      Teuchos::rcp(new N_NLS_NOX::PseudoTransientTest(maxStepSize, 1.0e-8));

    ctest->addStatusTest(mi);
    ctest->addStatusTest(fv);
    ctest->addStatusTest(pt);

    // First solve - pseudo transient

    Teuchos::RefCountPtr<N_NLS_NOX::AugmentLinSys> als =
      paramsPtr->createAugmentLinearSystem(lasSysPtr_);

    solverPtr_ = Teuchos::rcp(new N_NLS_NOX::PseudoTransientBased(als,
								  groupPtr_,
								  ctest,
								  noxList,
								  initialStepSize,
								  minStepSize,
								  maxStepSize));

    NOX::StatusTest::StatusType status = solverPtr_->solve();
    firstSolveComplete_ = true;

    // RPP 3/7/2006: We don't care if pseudo transient solve fails the
    // solve() call above.  This is just to get an inital guess for
    // the corrector step in the next solve.  So we don't check the
    // status at this point.

    // Turn off pseudo transient in groups. (this is also done in the
    // pseudo transient solver, but just being safe - groups could be
    // different).
    Teuchos::RefCountPtr<N_NLS_NOX::AugmentLinSys> dummy;
    groupPtr_->setAugmentLinearSystem(false, dummy);

    // Second solve is the correct steady state solve
    solverPtr_ = NOX::Solver::buildSolver(groupPtr_,
					  paramsPtr->getStatusTests(),
					  paramsPtr->getNoxParams());
    status = solverPtr_->solve();

    // Send back the correct return code
    return (paramsPtr->getStatusTestReturnCode());
  }  // continuation = 4 + power node
  else if (solverType == 10)
  {

    // Get some initial objects
    Teuchos::RefCountPtr<Teuchos::ParameterList> locaList = paramsPtr->getLocaParams();
    Teuchos::ParameterList& stepperList = locaList->sublist("Stepper");
    LOCA::ParameterVector locaPVec;

    // Create storage for continuation objects
    int numHomotopyContinuationRuns = 4;
    std::vector<std::string> names(numHomotopyContinuationRuns);
    std::vector<double> initialVal(numHomotopyContinuationRuns);
    std::vector<double> finalVal(numHomotopyContinuationRuns);
    std::vector<double> minVal(numHomotopyContinuationRuns);
    std::vector<double> maxVal(numHomotopyContinuationRuns);

    // Set up continuation steps
    // ***************************************************
    // Changes vals below for continuation
    // ***************************************************
    names[0] = stepperList.get("Power Node", "VA:V0");
    names[1] = "mosfet:gainscale";
    names[2] = "mosfet:nltermscale";
    names[3] = "mosfet:sizescale";
    initialVal[0] = 0.0;
    initialVal[1] = 0.0;
    initialVal[2] = 0.0;
    initialVal[3] = 0.0;
    finalVal[0] = 1.0;
    finalVal[1] = 1.0;
    finalVal[2] = 1.0;
    finalVal[3] = 1.0;
    // ***************************************************
    //std::string n1 = "mosfet:nltermscale";
    //locaPVec.addParameter(n1, 0.0);
    //std::string n2 = "mosfet:sizescale";
    //locaPVec.addParameter(n2, 0.0);
    // ***************************************************
    // ***************************************************

    // Ste up max/min bounds
    for (int i = 0; i < names.size(); ++i) {
      if (finalVal[i] > initialVal[i]) {
        minVal[i] = initialVal[i];
        maxVal[i] = finalVal[i];
      }
      else {
        minVal[i] = finalVal[i];
        maxVal[i] = initialVal[i];
      }
    }

    // Initialize loca parameter vector
    for (int i = 0; i < names.size(); ++i)
      locaPVec.addParameter(names[i], initialVal[i]);

    LOCA::Abstract::Iterator::IteratorStatus locaStatus;

    // Loop over the number of homotopy steps
    for (int hs = 0; hs < names.size(); ++hs) {
      for (int i = 0; i < names.size(); ++i) {
        if (i >= hs)
          locaPVec.setValue(names[i], initialVal[i]);
        else
          locaPVec.setValue(names[i], finalVal[i]);
      }
      groupPtr_->setParams(locaPVec);
      stepperList.set("Continuation Parameter", names[hs]);
      stepperList.set("Initial Value", initialVal[hs]);
      stepperList.set("Min Value", minVal[hs]);
      stepperList.set("Max Value", maxVal[hs]);

      // Initialize parameters in xyce
      groupPtr_->computeF();

      // Do the continuation run
      resetStepper(globalDataPtr_, groupPtr_, locaStatusTestPtr_, paramsPtr->getAllParams());
      locaStatus = stepperPtr_->run();

      // Kick out if continuation failed
      if (locaStatus != LOCA::Abstract::Iterator::Finished)
          return (-1);

      firstSolveComplete_ = true;

      // Copy out the solution and use it in the next run
      groupPtr_->copy(*(stepperPtr_->getSolutionGroup()));
    }

    // Return converged solver code
    return (paramsPtr->getStatusTestReturnCode());

  }
  else if (solverType == 33)  // artificial parameter
  {

#ifdef Xyce_NOX_LOCA_ARTIFICIAL_HOMOTOPY_SUPPORT
    Teuchos::RefCountPtr<Teuchos::ParameterList> locaList = paramsPtr->getLocaParams();
    Teuchos::ParameterList& locaUtilsList = locaList->sublist("Utilities");

    Teuchos::RefCountPtr<LOCA::Homotopy::Group> hGrp =
      Teuchos::rcp(new LOCA::Homotopy::Group(*locaList, globalDataPtr_, groupPtr_, 1.0, 0.0));

    hGrp->computeF();

    locaList->sublist("Predictor").set("Secant", 0.999);
    locaList->sublist("Stepper").set("Max Value", 0.999);

    resetStepper(globalDataPtr_, hGrp, locaStatusTestPtr_, paramsPtr->getAllParams());

    LOCA::Abstract::Iterator::IteratorStatus locaStatus = stepperPtr_->run();
    firstSolveComplete_ = true;

    Teuchos::RefCountPtr<LOCA::Homotopy::Group> hGrp2 =
      Teuchos::rcp(new LOCA::Homotopy::Group(*locaList, globalDataPtr_, groupPtr_, 0.1, 1.0));

    locaList->sublist("Predictor").set("Secant", 0.999);
    locaList->sublist("Stepper").set("Initial Value", 0.999);
    locaList->sublist("Stepper").set("Max Value", 1.0);
    locaList->sublist("Step Size").set("Method", "Constant");
    locaList->sublist("Step Size").set("Initial Step Size", 0.0001);
    locaList->sublist("Step Size").set("Min Step Size", 0.0001);

    resetStepper(globalDataPtr_, hGrp2, locaStatusTestPtr_, paramsPtr->getAllParams());

    locaStatus = stepperPtr_->run();

    // Return the solution status
    if (locaStatus != LOCA::Abstract::Iterator::Finished)
      return (-1);
    else
      return (paramsPtr->getStatusTestReturnCode());

#else
    Xyce::Report::UserFatal0() << "Nonlinear Solver (NOX::Interface) Artificial parameter continuation requires "
                               << "building xyce with the define: -DXyce_NOX_LOCA_ARTIFICIAL_HOMOTOPY_SUPPORT to "
                               << "allow LOCA to augment the diagonal of Jacobian! Either rebuild Xyce or do not "
                               << "run Xyce with \"continuation=33\"";
#endif
  } // End of if (solverType == )

 } // try
 catch (const char* error_msg) {
   std::string nox_error = "NOX Error";
   std::string err_msg = std::string(error_msg);
   if (err_msg == nox_error) {
     const std::string message =
       "Caught a NOX Exception in N_NLS_NOX::Interface::solve()!";
     N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, message);
   }
   else // Otherwise, rethrow...
     throw;
 }
#ifndef Xyce_CHARON
 catch (const std::exception& e) {
   Xyce::dout() << e.what() << std::endl;
   const std::string message =
     "Caught std::exception in N_NLS_NOX::Interface::solve()!";
   N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, message);
 }
 catch (...) {
   const std::string message =
     "Caught Unknown Exception in N_NLS_NOX::Interface::solve()!";
   N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, message);
 }
#endif

  // Should never get this far
  return -1;
}


//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::Interface::opStartCont0
// Purpose       :
// Special Notes : returns true if DCOP restart is being used.
//
//                 The "found" variable indicates if any of the nodes specified
//                 in the dcop start file were found in this circuit.  If not,
//                 then don't bother with this.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Electrical and Microsystem Modeling
// Creation Date : 09/15/07
//-----------------------------------------------------------------------------
bool Interface::opStartCont0 (ParameterSet* paramsPtr)
{
  bool usedOP(false);

#ifdef Xyce_DEBUG_OP_START
  Xyce::dout() << "NOX_Interface:  Inside continuation=0 OP_START code (case 1)" << std::endl;
#endif
  if (!DCOPused_)
  {
    int found = 0;
    std::string icType;
    Xyce::NodeNamePairMap & op = outMgrPtr_->getICData(found, icType);
    Xyce::NodeNamePairMap & allNodes = outMgrPtr_->getAllNodes( );
#ifdef Xyce_PARALLEL_MPI
    N_PDS_Comm * pdsCommPtr = pdsMgrPtr_->getPDSComm();
#endif

    if (found > 0 && icType == "DCOP_RESTART")
    {
      DCOPused_ = true;
      usedOP = true;
      Teuchos::RefCountPtr<N_NLS_NOX::AugmentLinSys> als =
#ifdef Xyce_PARALLEL_MPI
      paramsPtr->createAugmentLinearSystem(lasSysPtr_, op, allNodes, pdsCommPtr);
#else
      paramsPtr->createAugmentLinearSystem(lasSysPtr_, op, allNodes);
#endif
      groupPtr_->setAugmentLinearSystem(true, als);
      NOX::StatusTest::StatusType status = solverPtr_->solve();
      // Create a new solver after performing the initial op_start solve.
      solverPtr_ = NOX::Solver::buildSolver(groupPtr_,
                                          paramsPtr->getStatusTests(),
                                          paramsPtr->getNoxParams());
      firstSolveComplete_ = true;
      groupPtr_->setAugmentLinearSystem(false, als);
      anaIntPtr_->completeOPStartStep();

#ifdef Xyce_DEBUG_OP_START
      // DNS: Uncommenting this will set debug output of linear system for every step
#ifdef Xyce_PARALLEL_MPI
      groupPtr_->setOutputLinear (&op, &allNodes, pdsCommPtr);
#else
      groupPtr_->setOutputLinear (&op, &allNodes);
#endif
#endif // debug op start

    }
  }
  return usedOP;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::Interface::opStartCont1
// Purpose       :
// Special Notes : returns true if DCOP restart is being used.
//
//                 The "found" variable indicates if any of the nodes specified
//                 in the dcop start file were found in this circuit.  If not,
//                 then don't bother with this.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Electrical and Microsystem Modeling
// Creation Date : 09/15/07
//-----------------------------------------------------------------------------
bool Interface::opStartCont1 (ParameterSet* paramsPtr)
{
  bool usedOP(false);

#ifdef Xyce_DEBUG_OP_START
  Xyce::dout() << "NOX_Interface:  Inside continuation=1 OP_START code (case 2)" << std::endl;
#endif
  if (!DCOPused_)
  {
    int found = 0;
    std::string icType;
    Xyce::NodeNamePairMap & op = outMgrPtr_->getICData(found, icType);
    Xyce::NodeNamePairMap & allNodes = outMgrPtr_->getAllNodes( );
#ifdef Xyce_PARALLEL_MPI
    N_PDS_Comm * pdsCommPtr = pdsMgrPtr_->getPDSComm();
#endif

    if (found > 0 && icType == "DCOP_RESTART")
    {
      DCOPused_ = true;
      usedOP = true;
      // Set up nox nonlinear solver
      if (Teuchos::is_null(solverPtr_))
      {
          solverPtr_ = NOX::Solver::buildSolver(groupPtr_,
					      paramsPtr->getStatusTests(),
					      paramsPtr->getNoxParams());
      }

      Teuchos::RefCountPtr<N_NLS_NOX::AugmentLinSys> als =
#ifdef Xyce_PARALLEL_MPI
      paramsPtr->createAugmentLinearSystem(lasSysPtr_, op, allNodes, pdsCommPtr);
#else
      paramsPtr->createAugmentLinearSystem(lasSysPtr_, op, allNodes);
#endif
      groupPtr_->setAugmentLinearSystem(true, als);
      NOX::StatusTest::StatusType status = solverPtr_->solve();

      firstSolveComplete_ = true;
      groupPtr_->setAugmentLinearSystem(false, als);
      anaIntPtr_->completeOPStartStep();

#ifdef Xyce_DEBUG_OP_START
      // DNS: Uncommenting this will set debug output of linear system for every step
#ifdef Xyce_PARALLEL_MPI
      groupPtr_->setOutputLinear (&op, &allNodes, pdsCommPtr);
#else
      groupPtr_->setOutputLinear (&op, &allNodes);
#endif
#endif // debug op start

    }
  }
  return usedOP;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::Interface::icCont
// Purpose       :
// Special Notes : returns true if IC is being used.
//
//                 The "found" variable indicates if any of the nodes specified
//                 in the dcop start file were found in this circuit.  If not,
//                 then don't bother with this.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Electrical and Microsystem Modeling
// Creation Date : 09/15/07
//-----------------------------------------------------------------------------
bool Interface::icCont (ParameterSet* paramsPtr)
{
  bool usedIC(false);

#ifdef Xyce_DEBUG_IC
  Xyce::dout() << "NOX_Interface:  Inside continuation=0 .IC code." << std::endl;
#endif
  int found = 0;
  std::string icType;
  Xyce::NodeNamePairMap & op = outMgrPtr_->getICData(found, icType);

  usedIC = (icType=="IC" && found > 0);
  if (usedIC)
  {
    bool useGminStepping=false;
    Teuchos::RefCountPtr<N_NLS_NOX::AugmentLinSys> als =
      paramsPtr->createAugmentLinearSystem(lasSysPtr_, op, useGminStepping);
    groupPtr_->setAugmentLinearSystem(true, als);
  }
  return usedIC;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::Interface::icCont3
// Purpose       : IC with gmin stepping
// Special Notes : returns true if IC is being used.
//
//                 The "found" variable indicates if any of the nodes specified
//                 in the dcop start file were found in this circuit.  If not,
//                 then don't bother with this.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 04/29/2012
//-----------------------------------------------------------------------------
bool Interface::icCont3 (ParameterSet* paramsPtr)
{
  bool usedIC(false);

#ifdef Xyce_DEBUG_IC
  Xyce::dout() << "NOX_Interface:  Inside continuation=3 .IC code." << std::endl;
#endif
  int found = 0;
  std::string icType = "";
  Xyce::NodeNamePairMap & op = outMgrPtr_->getICData(found, icType);

  usedIC = (icType=="IC" && found > 0);
  if (usedIC)
  {
    bool useGminStepping=true;
    Teuchos::RefCountPtr<N_NLS_NOX::AugmentLinSys> als =
      paramsPtr->createAugmentLinearSystem(lasSysPtr_, op, useGminStepping);
    groupPtr_->setAugmentLinearSystem(true, als);
  }
  return usedIC;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::Interface::nodesetCont0
// Purpose       :
// Special Notes : returns true if is being used.
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Electrical and Microsystem Modeling
// Creation Date : 09/15/07
//-----------------------------------------------------------------------------
bool Interface::nodesetCont0 (ParameterSet* paramsPtr)
{
  bool usedNODESET(false);

#ifdef Xyce_DEBUG_IC
  Xyce::dout() << "NOX_Interface:  Inside continuation=0 .NODESET code (case 1)" << std::endl;
#endif
  int found = 0;
  std::string icType;
  Xyce::NodeNamePairMap & op = outMgrPtr_->getICData(found, icType);
#ifdef Xyce_PARALLEL_MPI
  N_PDS_Comm * pdsCommPtr = pdsMgrPtr_->getPDSComm();
#endif

  usedNODESET = (icType=="NODESET" && found > 0);
  if (usedNODESET)
  {
    Teuchos::RefCountPtr<N_NLS_NOX::AugmentLinSys> als =
      paramsPtr->createAugmentLinearSystem(lasSysPtr_, op);

    groupPtr_->setAugmentLinearSystem(true, als);
    NOX::StatusTest::StatusType status = solverPtr_->solve();

    // Create a new solver after performing the initial nodeset solve.
    solverPtr_ = NOX::Solver::buildSolver(groupPtr_,
                                          paramsPtr->getStatusTests(),
                                          paramsPtr->getNoxParams());
    firstSolveComplete_ = true;
    groupPtr_->setAugmentLinearSystem(false, als);
    anaIntPtr_->completeOPStartStep();
  }
  return usedNODESET;
}


//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::Interface::nodesetCont1
// Purpose       :
// Special Notes : returns true if .NODESET is being used.
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Electrical and Microsystem Modeling
// Creation Date : 09/15/07
//-----------------------------------------------------------------------------
bool Interface::nodesetCont1 (ParameterSet* paramsPtr)
{
  bool usedNODESET(false);

#ifdef Xyce_DEBUG_IC
  Xyce::dout() << "NOX_Interface:  Inside continuation=1 .NODESET code (case 2)" << std::endl;
#endif

  int found = 0;
  std::string icType;
  Xyce::NodeNamePairMap & op = outMgrPtr_->getICData(found, icType);
#ifdef Xyce_PARALLEL_MPI
  N_PDS_Comm * pdsCommPtr = pdsMgrPtr_->getPDSComm();
#endif

  usedNODESET = (icType=="NODESET" && found > 0);
  if (usedNODESET)
  {
    // Set up nox nonlinear solver
    if (Teuchos::is_null(solverPtr_))
    {
      solverPtr_ = NOX::Solver::buildSolver(groupPtr_,
                                            paramsPtr->getStatusTests(),
                                            paramsPtr->getNoxParams());
    }

    Teuchos::RefCountPtr<N_NLS_NOX::AugmentLinSys> als =
      paramsPtr->createAugmentLinearSystem(lasSysPtr_, op);

    groupPtr_->setAugmentLinearSystem(true, als);
    NOX::StatusTest::StatusType status = solverPtr_->solve();

    firstSolveComplete_ = true;
    groupPtr_->setAugmentLinearSystem(false, als);
    anaIntPtr_->completeOPStartStep();
  }
  return usedNODESET;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::Interface::takeFirstSolveStep
// Purpose       : same as Interface::solve, except that solverPtr_->iterate is
//                 called instead of solverPtr_->solve.
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
int Interface::takeFirstSolveStep (N_NLS_NonLinearSolver * nlsTmpPtr)
{
  // For base object
  N_NLS_NonLinearSolver::resetCountersAndTimers_();

  // Pick the parameter set to use.
  ParameterSet* paramsPtr;
  if ((usemode_) && (mode_ == TRANSIENT))
    paramsPtr = &transientParams_;
  else if ((usemode_) && (mode_ == HB_MODE))
    paramsPtr = &hbParams_;
  else
    paramsPtr = &dcParams_;

  // set the xyce return codes:
  paramsPtr->setStatusTestReturnCodes(retCodes_);

  if (Teuchos::is_null(globalDataPtr_))
    globalDataPtr_ = LOCA::createGlobalData(paramsPtr->getAllParams());

  // Set up the shared system (we have to redo this every time because
  // the object pointed to by nextSolVectorPtrPtr may have changed.
  delete sharedSystemPtr_;
  sharedSystemPtr_ = new SharedSystem(**nextSolVectorPtrPtr_,
				      *rhsVectorPtr_,
				      *jacobianMatrixPtr_,
				      *NewtonVectorPtr_,
				      *gradVectorPtr_,
				      *lasSysPtr_,
				      *this);

  // Reset up the corresponding group as well
  //delete groupPtr_;
  if (nlsTmpPtr==0)
  {
    if (Teuchos::is_null(groupPtr_))
    {
       Xyce::dout() << "takeFirstSolveStep: allocating a new group!" << std::endl;
       groupPtr_ = Teuchos::rcp(new N_NLS_LOCA::Group(globalDataPtr_,
						      *sharedSystemPtr_,
						      getLoader(),
						      *outMgrPtr_,
						      *anaIntPtr_));
    }
    else
    {
       Xyce::dout() << "takeFirstSolveStep: using the old group!" << std::endl;
      N_NLS_NOX::Vector tmpVec(**nextSolVectorPtrPtr_, *lasSysPtr_);
      groupPtr_->setX(tmpVec);
    }
  }
  else
  {
    Xyce::dout() << "takeFirstSolveStep: copying over the passed group!" << std::endl;
    copiedGroupFlag_ = true;
    Interface * nlsOtherPtr = dynamic_cast<Interface*>(nlsTmpPtr);
    groupPtr_ = nlsOtherPtr->getSolutionGroup();
  }

  // Set up solver
  if (Teuchos::is_null(solverPtr_))
    solverPtr_ = NOX::Solver::buildSolver(groupPtr_,
				          paramsPtr->getStatusTests(),
			                  paramsPtr->getNoxParams());
  else
    solverPtr_->reset(groupPtr_->getX());

  // Solve
  NOX::StatusTest::StatusType status = solverPtr_->step();

  // Return the solution status
  return (status == NOX::StatusTest::Converged) ? 1 : -1;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::Interface::takeOneSolveStep
// Purpose       : same as Interface::takeFirstSolveStep, except that none of the
//                 set up stuff (like allocating the solverPtr) is done here.
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
int Interface::takeOneSolveStep ()
{
  // Solve
  NOX::StatusTest::StatusType status = solverPtr_->step();

  // Return the solution status
  return (status == NOX::StatusTest::Converged) ? 1 : -1;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::Interface::getNumIterations
// Purpose       :
// Special Notes :
// Return Type   : Integer (current number of nonlinear iterations)
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
int Interface::getNumIterations() const
{
  // Pick the parameter set to use.
  const ParameterSet* paramsPtr;
  if ((usemode_) && (mode_ == TRANSIENT))
    paramsPtr = &transientParams_;
  else if ((usemode_) && (mode_ == HB_MODE))
    paramsPtr = &hbParams_;
  else
    paramsPtr = &dcParams_;

  int solverType = paramsPtr->getNoxSolverType();

  if ((!Teuchos::is_null(solverPtr_)) && (solverType == 0))
    return solverPtr_->getNumIterations();
  else if ((!Teuchos::is_null(solverPtr_)) && (solverType == 1))
    return solverPtr_->getNumIterations();
  else if ((!Teuchos::is_null(solverPtr_)) && (solverType == 9))
    return solverPtr_->getNumIterations();
  else if ((!Teuchos::is_null(stepperPtr_)) && (solverType != 0))
  {
    return stepperPtr_->getSolver()->getNumIterations();
  }

  // Sometimes this is called before solve() itself, in which calse
  // the solverPtr_ has not yet been initialized, so we just return 0.
  return 0;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::Interface::getMaxNormF() const
// Purpose       :
// Special Notes :
// Return Type   : double (norm of F)
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
double Interface::getMaxNormF() const
{
  // Pick the parameter set to use.
  const ParameterSet* paramsPtr;
  if ((usemode_) && (mode_ == TRANSIENT))
  {
    paramsPtr = &transientParams_;
  }
  else if ((usemode_) && (mode_ == HB_MODE))
  {
    paramsPtr = &hbParams_;
  }
  else
  {
    paramsPtr = &dcParams_;
  }

  double maxNormF = paramsPtr->getMaxNormF();
  return maxNormF;
}


//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::Interface::getMaxNormFindex() const
// Purpose       :
// Special Notes :
// Return Type   : int (vector index norm of F)
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
int Interface::getMaxNormFindex() const
{
  // Pick the parameter set to use.
  const ParameterSet* paramsPtr;
  if ((usemode_) && (mode_ == TRANSIENT))
  {
    paramsPtr = &transientParams_;
  }
  else if ((usemode_) && (mode_ == HB_MODE))
  {
    paramsPtr = &hbParams_;
  }
  else
  {
    paramsPtr = &dcParams_;
  }

  int maxNormFindex = paramsPtr->getMaxNormFindex();
  return maxNormFindex;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::Interface::getDebugLevel
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 9/17/2007
//-----------------------------------------------------------------------------
int Interface::getDebugLevel() const
{
  const ParameterSet* paramsPtr;
  if ((usemode_) && (mode_ == TRANSIENT))
  {
    paramsPtr = &transientParams_;
  }
  else if ((usemode_) && (mode_ == HB_MODE))
  {
    paramsPtr = &hbParams_;
  }
  else
  {
    paramsPtr = &dcParams_;
  }

  return (paramsPtr->getDebugLevel());
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::Interface::getScreenOutputFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 9/17/2007
//-----------------------------------------------------------------------------
bool Interface::getScreenOutputFlag () const
{
  const ParameterSet* paramsPtr;
  if ((usemode_) && (mode_ == TRANSIENT))
  {
    paramsPtr = &transientParams_;
  }
  else if ((usemode_) && (mode_ == HB_MODE))
  {
    paramsPtr = &hbParams_;
  }
  else
  {
    paramsPtr = &dcParams_;
  }

  return (paramsPtr->getScreenOutputFlag());
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::Interface::getDebugMinTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 9/17/2007
//-----------------------------------------------------------------------------
double Interface::getDebugMinTime() const
{
  const ParameterSet* paramsPtr;
  if ((usemode_) && (mode_ == TRANSIENT))
  {
    paramsPtr = &transientParams_;
  }
  else if ((usemode_) && (mode_ == HB_MODE))
  {
    paramsPtr = &hbParams_;
  }
  else
  {
    paramsPtr = &dcParams_;
  }

  return (paramsPtr->getDebugMinTime());
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::Interface::getDebugMaxTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 9/17/2007
//-----------------------------------------------------------------------------
double Interface::getDebugMaxTime() const
{
  const ParameterSet* paramsPtr;
  if ((usemode_) && (mode_ == TRANSIENT))
  {
    paramsPtr = &transientParams_;
  }
  else if ((usemode_) && (mode_ == HB_MODE))
  {
    paramsPtr = &hbParams_;
  }
  else
  {
    paramsPtr = &dcParams_;
  }

  return (paramsPtr->getDebugMaxTime());
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::Interface::getDebugMinTimeStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 9/17/2007
//-----------------------------------------------------------------------------
int Interface::getDebugMinTimeStep() const
{
  const ParameterSet* paramsPtr;
  if ((usemode_) && (mode_ == TRANSIENT))
  {
    paramsPtr = &transientParams_;
  }
  else if ((usemode_) && (mode_ == HB_MODE))
  {
    paramsPtr = &hbParams_;
  }
  else
  {
    paramsPtr = &dcParams_;
  }

  return (paramsPtr->getDebugMinTimeStep());
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::Interface::getDebugMaxTimeStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 9/17/2007
//-----------------------------------------------------------------------------
int Interface::getDebugMaxTimeStep() const
{
  const ParameterSet* paramsPtr;
  if ((usemode_) && (mode_ == TRANSIENT))
  {
    paramsPtr = &transientParams_;
  }
  else if ((usemode_) && (mode_ == HB_MODE))
  {
    paramsPtr = &hbParams_;
  }
  else
  {
    paramsPtr = &dcParams_;
  }

  return (paramsPtr->getDebugMaxTimeStep());
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::Interface::getMMFormat
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 3/2/2011
//-----------------------------------------------------------------------------
bool Interface::getMMFormat () const
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::Interface::isFirstContinuationParam
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool Interface::isFirstContinuationParam() const
{
  return isFirstContinuationParam_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::Interface::isFirstSolveComplete
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool Interface::isFirstSolveComplete() const
{
  return firstSolveComplete_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::Interface::getContinuationStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
int Interface::getContinuationStep() const
{
  if (!Teuchos::is_null(stepperPtr_))
  {
    return stepperPtr_->getStepNumber();
  }
  return 0;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::Interface::getContinuationStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
int Interface::getParameterNumber() const
{
  return iParam_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::Interface::setAnalysisMode
//
// Purpose       : Specify the analysis mode to be used by the nonlinear
//                 solver in the next call to solve(). This *may* affect
//                 the parameters used by the solver.
//
// See Also      : setOptions, setTranOptions
//
// - Input Arguments -
//
//    mode       : Mode to be used in the next nonlinear solve.
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void Interface::setAnalysisMode(AnalysisMode mode)
{
  mode_ = mode;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::Interface::resetAll
// Purpose       : This is used when Xyce is doing a STEP loop, and
//                 needs to act like it is at the beginning of a transient
//                 simulation again, for the next parameter in the STEP loop.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
void Interface::resetAll (AnalysisMode mode)
{
  setAnalysisMode(mode);
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::Interface::copySolnVectors()
//
// Purpose       : To be called at the beginning of each nonlinear
//                 iteration. It equates the tmp and next vectors as
//                 well as some special hidden vectors.
//
// See Also      : setX0_
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool Interface::copySolnVectors()
{
  N_NLS_NonLinearSolver::setX0_();
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::Interface::getMatrixFreeFlag()
//
// Purpose       :
// See Also      :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, Ting Mei
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
bool Interface::getMatrixFreeFlag()
{
  return N_NLS_NonLinearSolver::getMatrixFreeFlag();
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::Interface::computeF()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool Interface::computeF()
{
  return N_NLS_NonLinearSolver::rhs_();
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::Interface::computeNewton
// Purpose       : Set up the parameters for the linear solver and then
//                 call newton_()
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool Interface::computeNewton(Teuchos::ParameterList& params)
{
  if (mode_ == DC_OP && setAZ_Tol_DC)
  {
    lasSolverPtr_->setTolerance(params.get("Tolerance", 1.0e-12));
  }
  else if (mode_ == TRANSIENT && setAZ_Tol_Transient)
  {
    lasSolverPtr_->setTolerance(params.get("Tolerance", 1.0e-12));
  }

  return N_NLS_NonLinearSolver::newton_();
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::Interface::computeJacobian()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool Interface::computeJacobian()
{
  bool status = N_NLS_NonLinearSolver::jacobian_();
  return status;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::Interface::applyJacobian(const N_LAS_Vector& input, N_LAS_Vector& result)
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool Interface::applyJacobian(const N_LAS_Vector& input, N_LAS_Vector& result)
{
  return N_NLS_NonLinearSolver::applyJacobian(input,result);
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::Interface::computeGradient()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool Interface::computeGradient()
{
  bool status = N_NLS_NonLinearSolver::gradient_();
  return status;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::Interface::getSolutionGroup
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
Teuchos::RefCountPtr<N_NLS_LOCA::Group> Interface::getSolutionGroup ()
{
  return groupPtr_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::Interface::getLoader
//
// Purpose       : LOCA needs access to loader to set parameters
//
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
N_LOA_Loader& Interface::getLoader() const
{
  return *(N_NLS_NonLinearSolver::loaderPtr_);
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::Interface::resetStepper
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void Interface::resetStepper(const Teuchos::RefCountPtr<LOCA::GlobalData>& gd,
			     const Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup>& initialGuess,
			     const Teuchos::RefCountPtr<NOX::StatusTest::Generic>& test,
			     const Teuchos::RefCountPtr<Teuchos::ParameterList>& p)
{
  stepperPtr_ =
    Teuchos::rcp(new LOCA::Stepper(gd, initialGuess, test, p));
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_NOX::Interface::getLocaFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool Interface::getLocaFlag ()
{
  // Pick the parameter set to use.
  const ParameterSet* paramsPtr;
  bool retCode;

  if ((usemode_) && (mode_ == TRANSIENT))
  {
    firstSolveComplete_ = false;
    paramsPtr = &transientParams_;
    int solverType = paramsPtr->getNoxSolverType();
    retCode = false;
    if (solverType != 0) retCode = true;
  }
  else
  {
    if ((usemode_) && (mode_ ==HB_MODE))
    {
      paramsPtr = &hbParams_;
    }
    else
    {
      paramsPtr = &dcParams_;
    }

    if (DCOPused_)
    {
      retCode = firstSolveComplete_;
    }
    else
    {
      int solverType = paramsPtr->getNoxSolverType();
      retCode=false;
      if (solverType != 0) retCode = true;
    }
  }

  return retCode;
}


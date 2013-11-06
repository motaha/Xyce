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
// Filename      : $RCSfile: N_ANP_AnalysisInterface.C,v $
// Purpose       :
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.36.2.3 $
// Revision Date  : $Date: 2013/10/03 17:23:31 $
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

#include <iostream>

// ----------   Xyce Includes   ----------

#include <N_ANP_AnalysisInterface.h>
#include <N_ANP_AnalysisManager.h>

#include <N_TIA_Assembler.h>
#include <N_TIA_DataStore.h>
#include <N_TIA_TimeIntInfo.h>
#include <N_TIA_TwoLevelError.h>

#include <N_UTL_OptionBlock.h>
#include <N_ERH_ErrorMgr.h>
#include <N_LOA_Loader.h>

#include <N_PDS_Manager.h>

#include <N_IO_CmdParse.h>

#include <N_UTL_BreakPoint.h>
#include <N_UTL_NoCase.h>

#include <N_TIA_TimeIntInfo.h>
#include <N_TIA_TwoLevelError.h>

#include <N_LOA_NonlinearEquationLoader.h>
#include <N_DEV_DeviceInterface.h>
#include <N_TOP_Topology.h>
#include <N_LAS_Builder.h>


//-----------------------------------------------------------------------------
// Function      : anpAnalysisModeToNLS
// Purpose       : Converte between N_NLS_Manager.h AnalysisMode enum and
//               : expanded ANP_AnalysisInterface.h ANP_Analysis_Mode enum.
// Special Notes :
// Scope         :
// Creator       : Todd Coffey, 1414
// Creation Date : 07/23/08
//-----------------------------------------------------------------------------
AnalysisMode anpAnalysisModeToNLS(ANP_Analysis_Mode mode)
{
  AnalysisMode outMode;
  if (mode == ANP_MODE_TRANSIENT)
  {
    outMode = TRANSIENT;
  }
  else if (mode == ANP_MODE_DC_OP)
  {
    outMode = DC_OP;
  }
  else if (mode == ANP_MODE_DC_SWEEP)
  {
    outMode = DC_SWEEP;
  }
  else if (mode == ANP_MODE_HB)
  {
    outMode = HB;
  }
  else
  {
    outMode = NUM_MODES; // Should be this be TRANSIENT?
  }
  return outMode;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::N_ANP_AnalysisInterface
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, 9233, Computational Sciences
// Creation Date : 11/10/00
//-----------------------------------------------------------------------------
N_ANP_AnalysisInterface::N_ANP_AnalysisInterface(N_IO_CmdParse & cp)
    : commandLine_(cp)
{
  anaManagerPtr_ = rcp(new N_ANP_AnalysisManager(commandLine_,this));
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::~N_ANP_AnalysisInterface
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/10/00
//-----------------------------------------------------------------------------
N_ANP_AnalysisInterface::~N_ANP_AnalysisInterface ()
{
  pkgOptMgrPtr_ = Teuchos::null;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::resetAll
// Purpose       : just like a destructor without the death
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, Rich Schiek, Ting Mei
// Creation Date : 7/24/08
//-----------------------------------------------------------------------------
void N_ANP_AnalysisInterface::resetAll()
{
  anaManagerPtr_->resetAll();
}


//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::registerTIAParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/06/00
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::registerTIAParams(const N_TIA_TIAParams & tiaParams_tmp)

{
  anaManagerPtr_->registerTIAParams(tiaParams_tmp);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::registerNonlinearEquationLoader
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd S. Coffey, 1414
// Creation Date : 07/23/08
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::registerNonlinearEquationLoader( N_LOA_NonlinearEquationLoader * nonlinearEquationLoaderPtr )
{
  anaManagerPtr_->registerNonlinearEquationLoader(nonlinearEquationLoaderPtr);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::registerDeviceInterface
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd S. Coffey, 1414
// Creation Date : 07/23/08
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::registerDeviceInterface( N_DEV_DeviceInterface * devInterfacePtr )
{
  anaManagerPtr_->registerDeviceInterface(devInterfacePtr);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::registerTopology
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd S. Coffey, 1414
// Creation Date : 07/23/08
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::registerTopology( N_TOP_Topology * topoMgrPtr )
{
  anaManagerPtr_->registerTopology( topoMgrPtr );
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::registerApplicationBuilder
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd S. Coffey, 1414
// Creation Date : 07/23/08
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::registerApplicationBuilder( N_LAS_Builder * appBuilderPtr )
{
  anaManagerPtr_->registerApplicationBuilder( appBuilderPtr );
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::registerPkgOptionsMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, 1437
// Creation Date : 10/21/08
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::registerPkgOptionsMgr( RCP<N_IO_PkgOptionsMgr> pkgOptPtr )
{
  pkgOptMgrPtr_ = pkgOptPtr;

  // this work was in the constructor, but now we don't know the pkgOptMgrPtr_ until
  // it is registered.  So, do this work now.
  string netListFile = "";
  if (commandLine_.getArgumentValue("netlist") != "")
  {
    netListFile = commandLine_.getArgumentValue("netlist");
  }

  pkgOptMgrPtr_->submitRegistration(
      "TIMEINT", netListFile, new N_ANP_AnalysisInterface_TranOptionsReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "TRAN", netListFile, new N_ANP_AnalysisInterface_TransAnalysisReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "DC", netListFile, new N_ANP_AnalysisInterface_DCAnalysisReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "OP", netListFile, new N_ANP_AnalysisInterface_OPAnalysisReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "STEP", netListFile, new N_ANP_AnalysisInterface_STEPAnalysisReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "OP_IO", netListFile, new N_ANP_AnalysisInterface_DCOPOptionsReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "SAVE", netListFile, new N_ANP_AnalysisInterface_SaveOptionsReg( this ) );

  // MPDE specific netlist lines
  pkgOptMgrPtr_->submitRegistration(
      "MPDE", netListFile, new N_ANP_AnalysisInterface_MPDE_AnalysisReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "MPDEINT", netListFile, new N_ANP_AnalysisInterface_MPDE_OptionsReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "TIMEINT-MPDE", netListFile, new N_ANP_AnalysisInterface_MPDE_TranMPDEOptionsReg( this ) );

  // HB Specific netlist lines
  pkgOptMgrPtr_->submitRegistration(
      "HB", netListFile, new N_ANP_AnalysisInterface_HB_AnalysisReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "HBINT", netListFile, new N_ANP_AnalysisInterface_HB_OptionsReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "LINSOL-HB", netListFile, new N_ANP_AnalysisInterface_HB_LinSolReg( this ) );

  // AC Specific netlist lines
  pkgOptMgrPtr_->submitRegistration(
      "AC", netListFile, new N_ANP_AnalysisInterface_AC_AnalysisReg( this ) );

  // MOR Specific netlist lines
  pkgOptMgrPtr_->submitRegistration(
      "MOR", netListFile, new N_ANP_AnalysisInterface_MOR_AnalysisReg( this ) );

  // MOR Specific netlist lines
  pkgOptMgrPtr_->submitRegistration(
      "MOR_OPTS", netListFile, new N_ANP_AnalysisInterface_MOR_OptionsReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "SENS", netListFile, new N_ANP_AnalysisInterface_SensOptionsReg( this ) );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::setTranAnalysisParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 9/27/00
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::setTranAnalysisParams(const N_UTL_OptionBlock &tiaParamsBlock)
{
  return anaManagerPtr_->setTranAnalysisParams(tiaParamsBlock);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::setDCAnalysisParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/6/00
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::setDCAnalysisParams(const N_UTL_OptionBlock & tiaParamsBlock)
{
  return anaManagerPtr_->setDCAnalysisParams(tiaParamsBlock);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::setACAnalysisParams
// Purpose       : Method to handle AC statements
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/18/07
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::setACAnalysisParams(const N_UTL_OptionBlock & ACParamBlock)
{
  return anaManagerPtr_->setACAnalysisParams(ACParamBlock);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::setOPAnalysisParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/15/07
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::setOPAnalysisParams(const N_UTL_OptionBlock & tiaParamsBlock)
{
  return anaManagerPtr_->setOPAnalysisParams(tiaParamsBlock);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::setSTEPAnalysisParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/30/03
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::setSTEPAnalysisParams(const N_UTL_OptionBlock & tiaParamsBlock)
{
  return anaManagerPtr_->setSTEPAnalysisParams(tiaParamsBlock);
}


//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::setSaveOptions
// Purpose       : Method to handle SAVE statements
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/18/07
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::setSaveOptions(const N_UTL_OptionBlock & tiaParamBlock)
{
  return anaManagerPtr_->setSaveOptions(tiaParamBlock);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::setDCOPRestartParams
// Purpose       : Method to handle DCOP restart statements
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/18/07
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::setDCOPRestartParams(const N_UTL_OptionBlock & tiaParamBlock)
{
  return anaManagerPtr_->setDCOPRestartParams(tiaParamBlock);
}


//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::setPauseTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 04/29/04
//-----------------------------------------------------------------------------
void N_ANP_AnalysisInterface::setPauseTime(const double pauseTime)
{
  anaManagerPtr_->setPauseTime(pauseTime);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::getPauseTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 04/29/04
//-----------------------------------------------------------------------------
double N_ANP_AnalysisInterface::getPauseTime()
{
  return(anaManagerPtr_->getPauseTime());
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::isPaused
//
// Purpose       : return the true if the simulation is currently paused
//
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Electrical and Microsystems Simulation
// Creation Date : 02/18/2008
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::isPaused()
{
  return(anaManagerPtr_->isPaused());
}
//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::resumeSimulation
// Purpose       : signal to control algorithm that simulation is continuing
//                 from previously paused simulation
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 04/29/04
//-----------------------------------------------------------------------------
void N_ANP_AnalysisInterface::resumeSimulation()
{
  anaManagerPtr_->resumeSimulation();
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::setTranOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 9/29/00
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::setTranOptions(const N_UTL_OptionBlock & OB)
{
  return anaManagerPtr_->setTranOptions(OB);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::setMPDEAnalysisParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 7/23/08
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::setMPDEAnalysisParams(const N_UTL_OptionBlock & OB)
{
  return anaManagerPtr_->setMPDEAnalysisParams(OB);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::setMPDEOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 7/23/08
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::setMPDEOptions(const N_UTL_OptionBlock & OB)
{
  return anaManagerPtr_->setMPDEOptions(OB);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::setHBAnalysisParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, Ting Mei
// Creation Date : 7/30/08
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::setHBAnalysisParams(const N_UTL_OptionBlock & OB)
{
  return anaManagerPtr_->setHBAnalysisParams(OB);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::setHBOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, Ting Mei
// Creation Date : 7/30/08
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::setHBOptions(const N_UTL_OptionBlock & OB)
{
  return anaManagerPtr_->setHBOptions(OB);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::setLinSol
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 07/12/2013
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::setLinSol(const N_UTL_OptionBlock & OB)
{
  return anaManagerPtr_->setLinSol(OB);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::setHBLinSol
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist
// Creation Date : 11/11/08
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::setHBLinSol(const N_UTL_OptionBlock & OB)
{
  return anaManagerPtr_->setHBLinSol(OB);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::setMORAnalysisParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist and Ting Mei
// Creation Date : 5/31/12
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::setMORAnalysisParams(const N_UTL_OptionBlock & OB)
{
  return anaManagerPtr_->setMORAnalysisParams(OB);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::setMOROptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist and Ting Mei
// Creation Date : 5/31/12
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::setMOROptions(const N_UTL_OptionBlock & OB)
{
  return anaManagerPtr_->setMOROptions(OB);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::setTRANMPDEOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 7/23/08
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::setTRANMPDEOptions(const N_UTL_OptionBlock & OB)
{
  return anaManagerPtr_->setTRANMPDEOptions(OB);
}


//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::setSensOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, Sandia
// Creation Date : 6/5/13
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::setSensOptions(const N_UTL_OptionBlock & OB)
{
  return anaManagerPtr_->setSensOptions(OB);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::registerElapsedTimer
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 02/24/06
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::registerElapsedTimer (N_UTL_Timer * et)
{
  return anaManagerPtr_->registerElapsedTimer(et);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::getBlockAnalysisFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/29/05
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::getBlockAnalysisFlag () const
{
  return anaManagerPtr_->getBlockAnalysisFlag();
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::getTransientFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/29/05
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::getTransientFlag () const
{
  return anaManagerPtr_->getTransientFlag();
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::initializeAll
// Purpose       : This function performs the final initializations.  Mainly,
//                 it initializes the two data container classes, DataStore and
//                 StepErrorControl.  It also registers the neccessary vectors
//                 with the LAS system class.
// Special Notes : This function should *only* be called after all the
//                 registrations have all been performed.  In particular, the
//                 N_LAS_System class *must* be registered before this function
//                 is called.
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/10/00
//-----------------------------------------------------------------------------

bool N_ANP_AnalysisInterface::initializeAll()
{
  return anaManagerPtr_->initializeAll();
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::registerLinearSystem
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/05/01
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::registerLinearSystem(N_LAS_System * lasSysPtr)
{
  anaManagerPtr_->registerLinearSystem(lasSysPtr);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::registerLoader
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/05/01
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::registerLoader(N_LOA_Loader * loaderPtr)
{
  return anaManagerPtr_->registerLoader(loaderPtr);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::registerOutputMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/05/01
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::registerOutputMgr(N_IO_OutputMgr * outputPtr)
{
  return anaManagerPtr_->registerOutputMgr(outputPtr);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::registerRestartMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/05/01
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::registerRestartMgr(N_IO_RestartMgr * restartPtr)
{
  return anaManagerPtr_->registerRestartMgr(restartPtr);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::registerNLSolver
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/05/01
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::registerNLSManager(N_NLS_Manager * nlsMgrPtr)
{
  return anaManagerPtr_->registerNLSManager(nlsMgrPtr);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::run
// Purpose       : Execute the top level control loop.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, 9233, Computational Sciences
// Creation Date : 11/10/00
//-----------------------------------------------------------------------------

bool N_ANP_AnalysisInterface::run()
{
  return anaManagerPtr_->run();
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::updateDerivs
// Purpose       : calls the time  int. method to update the  corrector
//                 derivatives.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/10/00
//-----------------------------------------------------------------------------

bool N_ANP_AnalysisInterface::updateDerivs()
{
  return anaManagerPtr_->updateDerivs ();
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::updateDivDiffs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/10/00
//-----------------------------------------------------------------------------

bool N_ANP_AnalysisInterface::updateDivDiffs()
{
  return anaManagerPtr_->updateDivDiffs();
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::getTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/27/00
//-----------------------------------------------------------------------------
double N_ANP_AnalysisInterface::getTime() const
{
  return anaManagerPtr_->getTime();
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::getCurrentTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Mems Modeling
// Creation Date : 8/21/2009
//-----------------------------------------------------------------------------
double N_ANP_AnalysisInterface::getCurrentTime() const
{
  return anaManagerPtr_->getCurrentTime();
}


//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::getFinalTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 9/04/01
//-----------------------------------------------------------------------------
double N_ANP_AnalysisInterface::getFinalTime() const
{
  return anaManagerPtr_->getFinalTime();
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::getInitialTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/04/02
//-----------------------------------------------------------------------------
double N_ANP_AnalysisInterface::getInitialTime() const
{
  return anaManagerPtr_->getInitialTime();
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::equateTmpVectors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/24/02
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::equateTmpVectors()
{
  return anaManagerPtr_->equateTmpVectors();
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::updateDerivsBlock
// Purpose       : calls the time  int. method to update the  corrector
//                 derivatives, but only does it a subset of the solution
//                 and state vectors.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/10/01
//-----------------------------------------------------------------------------
//inline
bool N_ANP_AnalysisInterface::updateDerivsBlock(const list<index_pair> & solGIDList,
                                      const list<index_pair> & staGIDList)
{
  return anaManagerPtr_->updateDerivsBlock(solGIDList, staGIDList);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::restartDataSize
// Purpose       : Gets the size of the restart data (bytes?).
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
int N_ANP_AnalysisInterface::restartDataSize( bool pack )
{
  return anaManagerPtr_->restartDataSize( pack );
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::dumpRestartData
// Purpose       : Dumps the restart data to a file.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::dumpRestartData(char * buf, int bsize, int & pos,
                                    N_PDS_Comm * comm, bool pack )
{
  return anaManagerPtr_->dumpRestartData(buf, bsize, pos, comm, pack);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::restoreRestartData
// Purpose       : Restores the restart data from a file.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::restoreRestartData(char * buf, int bsize, int & pos,
                                       N_PDS_Comm * comm, bool pack )
{
  return anaManagerPtr_->restoreRestartData(buf, bsize, pos, comm, pack);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::getSolnVarData
// Purpose       : Gets the solution variable data.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::getSolnVarData(const int & gid, vector<double> & varData)
{
  return anaManagerPtr_->getSolnVarData(gid, varData);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::getStateVarData
// Purpose       : Gets the state variable data.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::getStateVarData(const int & gid, vector<double> & varData)
{
  return anaManagerPtr_->getStateVarData(gid, varData);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::getStoreVarData
// Purpose       : Gets the store variable data.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::getStoreVarData(const int & gid, vector<double> & varData)
{
  return anaManagerPtr_->getStoreVarData(gid, varData);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::setSolnVarData
// Purpose       : Sets the solution variable data.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::setSolnVarData(const int & gid,
                                   const vector<double> & varData)
{
  return anaManagerPtr_->setSolnVarData(gid, varData);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::setStateVarData
// Purpose       : Sets the state variable data.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::setStateVarData(const int & gid,
                                    const vector<double> & varData)
{
  return anaManagerPtr_->setStateVarData(gid, varData);
}


//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::setStoreVarData
// Purpose       : Sets the state variable data.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::setStoreVarData(const int & gid,
                                    const vector<double> & varData)
{
  return anaManagerPtr_->setStoreVarData(gid, varData);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::setBeginningIntegrationFlag
// Purpose       : set beginning integration flag to true
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytems Modeling
// Creation Date : 03/5/2010
//-----------------------------------------------------------------------------
void N_ANP_AnalysisInterface::setBeginningIntegrationFlag(bool bif)
{
  anaManagerPtr_->setBeginningIntegrationFlag (bif);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::registerParallelServices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/30/01
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::registerParallelServices(N_PDS_Manager * pds)
{
  return anaManagerPtr_->registerParallelServices(pds);

}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::outputSummary
// Purpose       : Outputs summary information for time-integration function.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 07/30/02
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::outputSummary()
{
  bool bsuccess = true;

  if (!Teuchos::is_null(anaManagerPtr_))
  {
    bsuccess = anaManagerPtr_->printLoopInfo(0,0);
  }

  return bsuccess;

}


//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::setNextSolVectorPtr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 04/22/03
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::setNextSolVectorPtr (N_LAS_Vector * solVecPtr)
{
  return anaManagerPtr_->setNextSolVectorPtr (solVecPtr);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::loadRHS
// Purpose       : manages the RHS load.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 03/03/04
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::loadRHS ()
{
  return anaManagerPtr_->assemblerPtr->loadRHS();
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::loadJacobian
// Purpose       : Manages the jacobian load.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 03/03/04
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::loadJacobian ()
{
  return anaManagerPtr_->assemblerPtr->loadJacobian();
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::applyJacobian
// Purpose       : Manages the jacobian load.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 03/03/04
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::applyJacobian (const N_LAS_Vector& input, N_LAS_Vector& result)
{
  return anaManagerPtr_->assemblerPtr->applyJacobian(input,result);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::getResidualTime
// Purpose       : This should give the time for the most recent Residual load.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 03/03/04
//-----------------------------------------------------------------------------
double N_ANP_AnalysisInterface::getResidualTime()
{
  return anaManagerPtr_->assemblerPtr->residualTime_;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::getJacobianTime
// Purpose       : This should give the time for the most recent Jacobian load.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 03/03/04
//-----------------------------------------------------------------------------
double N_ANP_AnalysisInterface::getJacobianTime()
{
  return anaManagerPtr_->assemblerPtr->jacobianTime_;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::getTIAParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 11/05/04
//-----------------------------------------------------------------------------
N_TIA_TIAParams& N_ANP_AnalysisInterface::getTIAParams()
{
  return anaManagerPtr_->tiaParams;
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::simulationComplete
// Purpose       : return boolean signifying whether simulation is complete
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 05/06/2004
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::simulationComplete()
{
  return anaManagerPtr_->simulationComplete();
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::completeOPStartStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 06/28/2006
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::completeOPStartStep()
{
  return anaManagerPtr_->completeOPStartStep();
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::completeHomotopyStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, 9233
// Creation Date : 03/20/2006
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::completeHomotopyStep
    ( const vector<string> & paramNames,
      const vector<double> & paramVals,
      N_LAS_Vector * solnVecPtr )
{
  return anaManagerPtr_->completeHomotopyStep(paramNames, paramVals,solnVecPtr);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::failHomotopyStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, 9233
// Creation Date : 03/30/2006
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::failHomotopyStep ( )
{
  return anaManagerPtr_->failHomotopyStep ();
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::startTimeStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/20/2006
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::startTimeStep (const N_TIA_TimeIntInfo & tiInfo)
{
  return anaManagerPtr_->startTimeStep (tiInfo);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::runStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/06/2006
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::runStep
    (const N_TIA_TimeIntInfo & tiInfo, N_TIA_TwoLevelError & tlError)
{
  return anaManagerPtr_->runStep (tiInfo, tlError);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::startupSolvers
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/06/2006
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::startupSolvers ()
{
  return anaManagerPtr_->startupSolvers ();
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::finishSolvers
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/06/2006
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::finishSolvers ()
{
  return anaManagerPtr_->finishSolvers ();
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::provisionalStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/04/09
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::provisionalStep (double maxTimeStep,
      double &currTimeStep)
{
  return anaManagerPtr_->provisionalStep (maxTimeStep, currTimeStep);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::acceptProvisionalStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 04/10/09
//-----------------------------------------------------------------------------
void N_ANP_AnalysisInterface::acceptProvisionalStep ()
{
  anaManagerPtr_->acceptProvisionalStep ();
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::rejectProvisionalStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 04/10/09
//-----------------------------------------------------------------------------
void N_ANP_AnalysisInterface::rejectProvisionalStep ()
{
  anaManagerPtr_->rejectProvisionalStep ();
}

#ifdef Xyce_Dakota

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::getDakotaRunFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::getDakotaRunFlag()
{
  return anaManagerPtr_->getDakotaRunFlag();
}


//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::setDakotaRunFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL
// Creation Date :
//-----------------------------------------------------------------------------
void N_ANP_AnalysisInterface::setDakotaRunFlag( bool flag )
{
  anaManagerPtr_->setDakotaRunFlag( flag );
}


//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::getDakotaIteration
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL
// Creation Date :
//-----------------------------------------------------------------------------
int N_ANP_AnalysisInterface::getDakotaIteration()
{
    return anaManagerPtr_->getDakotaIteration();
}


//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::setDakotaIteration
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL
// Creation Date :
//-----------------------------------------------------------------------------
void N_ANP_AnalysisInterface::setDakotaIteration( int iterNumber )
{
  anaManagerPtr_->setDakotaIteration( iterNumber );
}

#endif

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::condutanceTest
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
void N_ANP_AnalysisInterface::condutanceTest ()
{
  anaManagerPtr_->conductanceTest ();
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::homotopyStepSuccess
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
void N_ANP_AnalysisInterface::homotopyStepSuccess
    ( const vector<string> & paramNames,
      const vector<double> & paramVals)
{
  anaManagerPtr_->homotopyStepSuccess ( paramNames, paramVals);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::homotopyStepFailure
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/30/06
//-----------------------------------------------------------------------------
void N_ANP_AnalysisInterface::homotopyStepFailure ()
{
  anaManagerPtr_->homotopyStepFailure ();
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::stepSuccess
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
void N_ANP_AnalysisInterface::stepSuccess (int analysis)
{
  anaManagerPtr_->stepSuccess (analysis);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::stepFailure
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
void N_ANP_AnalysisInterface::stepFailure (int analysis)
{
  anaManagerPtr_->stepFailure (analysis);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::getInitialQnorm
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::getInitialQnorm (N_TIA_TwoLevelError & tle)
{
  return anaManagerPtr_->getInitialQnorm (tle);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::getBreakPoints
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::getBreakPoints (vector<N_UTL_BreakPoint> &breakPointTimes)
{
  return anaManagerPtr_->getBreakPoints (breakPointTimes);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::getTimeIntInfo
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 02/11/07
//-----------------------------------------------------------------------------
void N_ANP_AnalysisInterface::getTimeIntInfo(N_TIA_TimeIntInfo & tiInfo)
{
  anaManagerPtr_->getTimeIntInfo (tiInfo);
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::silenceProgress
// Purpose       : shut up "percent complete" noise
// Special Notes :
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 11 Feb 2009
//-----------------------------------------------------------------------------
void N_ANP_AnalysisInterface::silenceProgress()
{
  anaManagerPtr_->silenceProgress();
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::enableProgress
// Purpose       : enable "percent complete" noise
// Special Notes :
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 11 Feb 2009
//-----------------------------------------------------------------------------
void N_ANP_AnalysisInterface::enableProgress()
{
  anaManagerPtr_->enableProgress();
}

//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::initializeTransientModel
// Purpose       : ModelEvaluator Interface
// Special Notes :
// Scope         : public
// Creator       : Coffey, Schiek, Mei
// Creation Date : 05/29/09
//-----------------------------------------------------------------------------
void N_ANP_AnalysisInterface::initializeTransientModel()
{
  anaManagerPtr_->initializeTransientModel();
}
//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::evalModel
// Purpose       : ModelEvaluator Interface
// Special Notes :
// Scope         : public
// Creator       : Coffey, Schiek, Mei
// Creation Date : 05/27/09
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::evalTransientModel(
      double t,
      N_LAS_Vector * SolVectorPtr,
      N_LAS_Vector * CurrSolVectorPtr,
      N_LAS_Vector * LastSolVectorPtr,
      N_LAS_Vector * StaVectorPtr,
      N_LAS_Vector * CurrStaVectorPtr,
      N_LAS_Vector * LastStaVectorPtr,
      N_LAS_Vector * StaDerivVectorPtr,
      N_LAS_Vector * StoVectorPtr,
      N_LAS_Vector * CurrStoVectorPtr,
      N_LAS_Vector * LastStoVectorPtr,
      N_LAS_Vector * stoLeadCurrQCompVectorPtr,
      N_LAS_Vector * QVectorPtr,
      N_LAS_Vector * FVectorPtr,
      N_LAS_Vector * dFdxdVpVectorPtr,
      N_LAS_Vector * dQdxdVpVectorPtr,
      N_LAS_Matrix * dQdxMatrixPtr,
      N_LAS_Matrix * dFdxMatrixPtr
      )
{
  return (
        anaManagerPtr_->evalTransientModel(
        t,
        SolVectorPtr,
        CurrSolVectorPtr,
        LastSolVectorPtr,
        StaVectorPtr,
        CurrStaVectorPtr,
        LastStaVectorPtr,
        StaDerivVectorPtr,
        StoVectorPtr,
        CurrStoVectorPtr,
        LastStoVectorPtr,
        stoLeadCurrQCompVectorPtr,
        QVectorPtr,
        FVectorPtr,
        dFdxdVpVectorPtr,
        dQdxdVpVectorPtr,
        dQdxMatrixPtr,
        dFdxMatrixPtr
        )
      );
}
//-----------------------------------------------------------------------------
// Function      : N_ANP_AnalysisInterface::evalTransientModelState
// Purpose       : ModelEvaluator Interface
// Special Notes :
// Scope         : public
// Creator       : Coffey, Schiek, Mei
// Creation Date : 05/29/09
//-----------------------------------------------------------------------------
bool N_ANP_AnalysisInterface::evalTransientModelState(
      double t,
      N_LAS_Vector * SolVectorPtr,
      N_LAS_Vector * StaVectorPtr,
      N_LAS_Vector * StoVectorPtr
      )
{
  return (
        anaManagerPtr_->evalTransientModelState(
        t,
        SolVectorPtr,
        StaVectorPtr,
        StoVectorPtr
        )
      );
}



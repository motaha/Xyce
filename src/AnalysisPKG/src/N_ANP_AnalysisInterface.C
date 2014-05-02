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

//-----------------------------------------------------------------------------
// Filename      : $RCSfile: N_ANP_AnalysisInterface.C,v $
// Purpose       :
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.46 $
// Revision Date  : $Date: 2014/02/24 23:49:12 $
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>

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
#include <N_NLS_Manager.h>

namespace Xyce {
namespace Analysis {

//-----------------------------------------------------------------------------
// Function      : anpAnalysisModeToNLS
// Purpose       : Converte between N_NLS_Manager.h AnalysisMode enum and
//               : expanded ANP_AnalysisInterface.h ANP_Analysis_Mode enum.
// Special Notes :
// Scope         :
// Creator       : Todd Coffey, 1414
// Creation Date : 07/23/08
//-----------------------------------------------------------------------------
AnalysisMode anpAnalysisModeToNLS(Analysis_Mode mode)
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
    outMode = HB_MODE;
  }
  else
  {
    outMode = NUM_MODES; // Should be this be TRANSIENT?
  }
  return outMode;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::AnalysisInterface
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, 9233, Computational Sciences
// Creation Date : 11/10/00
//-----------------------------------------------------------------------------
AnalysisInterface::AnalysisInterface(N_IO_CmdParse & cp)
    : commandLine_(cp)
{
  anaManagerPtr_ = rcp(new AnalysisManager(commandLine_,this));
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::~AnalysisInterface
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/10/00
//-----------------------------------------------------------------------------
AnalysisInterface::~AnalysisInterface ()
{}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::resetAll
// Purpose       : just like a destructor without the death
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, Rich Schiek, Ting Mei
// Creation Date : 7/24/08
//-----------------------------------------------------------------------------
void AnalysisInterface::resetAll()
{
  anaManagerPtr_->resetAll();
}


//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::registerTIAParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/06/00
//-----------------------------------------------------------------------------
bool AnalysisInterface::registerTIAParams(const N_TIA_TIAParams & tiaParams_tmp)

{
  anaManagerPtr_->registerTIAParams(tiaParams_tmp);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::registerNonlinearEquationLoader
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd S. Coffey, 1414
// Creation Date : 07/23/08
//-----------------------------------------------------------------------------
bool AnalysisInterface::registerNonlinearEquationLoader( N_LOA_NonlinearEquationLoader * nonlinearEquationLoaderPtr )
{
  anaManagerPtr_->registerNonlinearEquationLoader(nonlinearEquationLoaderPtr);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::registerDeviceInterface
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd S. Coffey, 1414
// Creation Date : 07/23/08
//-----------------------------------------------------------------------------
bool AnalysisInterface::registerDeviceInterface( N_DEV_DeviceInterface * devInterfacePtr )
{
  anaManagerPtr_->registerDeviceInterface(devInterfacePtr);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::registerTopology
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd S. Coffey, 1414
// Creation Date : 07/23/08
//-----------------------------------------------------------------------------
bool AnalysisInterface::registerTopology( N_TOP_Topology * topoMgrPtr )
{
  anaManagerPtr_->registerTopology( topoMgrPtr );
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::registerApplicationBuilder
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd S. Coffey, 1414
// Creation Date : 07/23/08
//-----------------------------------------------------------------------------
bool AnalysisInterface::registerApplicationBuilder( N_LAS_Builder * appBuilderPtr )
{
  anaManagerPtr_->registerApplicationBuilder( appBuilderPtr );
  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::registerPkgOptionsMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, 1437
// Creation Date : 10/21/08
//-----------------------------------------------------------------------------
bool AnalysisInterface::registerPkgOptionsMgr( N_IO_PkgOptionsMgr *pkgOptPtr )
{
  pkgOptMgrPtr_ = pkgOptPtr;

  // this work was in the constructor, but now we don't know the pkgOptMgrPtr_ until
  // it is registered.  So, do this work now.
  std::string netListFile = "";
  if (commandLine_.getArgumentValue("netlist") != "")
  {
    netListFile = commandLine_.getArgumentValue("netlist");
  }

  pkgOptMgrPtr_->submitRegistration(
      "TIMEINT", netListFile, new AnalysisInterface_TranOptionsReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "TRAN", netListFile, new AnalysisInterface_TransAnalysisReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "DC", netListFile, new AnalysisInterface_DCAnalysisReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "OP", netListFile, new AnalysisInterface_OPAnalysisReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "STEP", netListFile, new AnalysisInterface_STEPAnalysisReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "OP_IO", netListFile, new AnalysisInterface_DCOPOptionsReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "SAVE", netListFile, new AnalysisInterface_SaveOptionsReg( this ) );

  // MPDE specific netlist lines
  pkgOptMgrPtr_->submitRegistration(
      "MPDE", netListFile, new AnalysisInterface_MPDE_AnalysisReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "MPDEINT", netListFile, new AnalysisInterface_MPDE_OptionsReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "TIMEINT-MPDE", netListFile, new AnalysisInterface_MPDE_TranMPDEOptionsReg( this ) );

  // HB Specific netlist lines
  pkgOptMgrPtr_->submitRegistration(
      "HB", netListFile, new AnalysisInterface_HB_AnalysisReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "HBINT", netListFile, new AnalysisInterface_HB_OptionsReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "LINSOL-HB", netListFile, new AnalysisInterface_HB_LinSolReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "LINSOL", netListFile, new AnalysisInterface_LinSolReg( this ) );

  // AC Specific netlist lines
  pkgOptMgrPtr_->submitRegistration(
      "AC", netListFile, new AnalysisInterface_AC_AnalysisReg( this ) );

  // MOR Specific netlist lines
  pkgOptMgrPtr_->submitRegistration(
      "MOR", netListFile, new AnalysisInterface_MOR_AnalysisReg( this ) );

  // MOR Specific netlist lines
  pkgOptMgrPtr_->submitRegistration(
      "MOR_OPTS", netListFile, new AnalysisInterface_MOR_OptionsReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "SENS", netListFile, new AnalysisInterface_SensOptionsReg( this ) );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::setTranAnalysisParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 9/27/00
//-----------------------------------------------------------------------------
bool AnalysisInterface::setTranAnalysisParams(const N_UTL_OptionBlock &tiaParamsBlock)
{
  return anaManagerPtr_->setTranAnalysisParams(tiaParamsBlock);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::setDCAnalysisParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/6/00
//-----------------------------------------------------------------------------
bool AnalysisInterface::setDCAnalysisParams(const N_UTL_OptionBlock & tiaParamsBlock)
{
  return anaManagerPtr_->setDCAnalysisParams(tiaParamsBlock);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::setACAnalysisParams
// Purpose       : Method to handle AC statements
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/18/07
//-----------------------------------------------------------------------------
bool AnalysisInterface::setACAnalysisParams(const N_UTL_OptionBlock & ACParamBlock)
{
  return anaManagerPtr_->setACAnalysisParams(ACParamBlock);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::setOPAnalysisParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/15/07
//-----------------------------------------------------------------------------
bool AnalysisInterface::setOPAnalysisParams(const N_UTL_OptionBlock & tiaParamsBlock)
{
  return anaManagerPtr_->setOPAnalysisParams(tiaParamsBlock);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::setSTEPAnalysisParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/30/03
//-----------------------------------------------------------------------------
bool AnalysisInterface::setSTEPAnalysisParams(const N_UTL_OptionBlock & tiaParamsBlock)
{
  return anaManagerPtr_->setSTEPAnalysisParams(tiaParamsBlock);
}


//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::setSaveOptions
// Purpose       : Method to handle SAVE statements
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/18/07
//-----------------------------------------------------------------------------
bool AnalysisInterface::setSaveOptions(const N_UTL_OptionBlock & tiaParamBlock)
{
  return anaManagerPtr_->setSaveOptions(tiaParamBlock);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::setDCOPRestartParams
// Purpose       : Method to handle DCOP restart statements
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/18/07
//-----------------------------------------------------------------------------
bool AnalysisInterface::setDCOPRestartParams(const N_UTL_OptionBlock & tiaParamBlock)
{
  return anaManagerPtr_->setDCOPRestartParams(tiaParamBlock);
}


//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::setPauseTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 04/29/04
//-----------------------------------------------------------------------------
void AnalysisInterface::setPauseTime(const double pauseTime)
{
  anaManagerPtr_->setPauseTime(pauseTime);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::getPauseTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 04/29/04
//-----------------------------------------------------------------------------
double AnalysisInterface::getPauseTime()
{
  return(anaManagerPtr_->getPauseTime());
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::isPaused
//
// Purpose       : return the true if the simulation is currently paused
//
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Electrical and Microsystems Simulation
// Creation Date : 02/18/2008
//-----------------------------------------------------------------------------
bool AnalysisInterface::isPaused()
{
  return(anaManagerPtr_->isPaused());
}
//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::resumeSimulation
// Purpose       : signal to control algorithm that simulation is continuing
//                 from previously paused simulation
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 04/29/04
//-----------------------------------------------------------------------------
void AnalysisInterface::resumeSimulation()
{
  anaManagerPtr_->resumeSimulation();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::setTranOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 9/29/00
//-----------------------------------------------------------------------------
bool AnalysisInterface::setTranOptions(const N_UTL_OptionBlock & OB)
{
  return anaManagerPtr_->setTranOptions(OB);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::setMPDEAnalysisParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 7/23/08
//-----------------------------------------------------------------------------
bool AnalysisInterface::setMPDEAnalysisParams(const N_UTL_OptionBlock & OB)
{
  return anaManagerPtr_->setMPDEAnalysisParams(OB);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::setMPDEOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 7/23/08
//-----------------------------------------------------------------------------
bool AnalysisInterface::setMPDEOptions(const N_UTL_OptionBlock & OB)
{
  return anaManagerPtr_->setMPDEOptions(OB);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::setHBAnalysisParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, Ting Mei
// Creation Date : 7/30/08
//-----------------------------------------------------------------------------
bool AnalysisInterface::setHBAnalysisParams(const N_UTL_OptionBlock & OB)
{
  return anaManagerPtr_->setHBAnalysisParams(OB);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::setHBOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, Ting Mei
// Creation Date : 7/30/08
//-----------------------------------------------------------------------------
bool AnalysisInterface::setHBOptions(const N_UTL_OptionBlock & OB)
{
  return anaManagerPtr_->setHBOptions(OB);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::setLinSol
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 07/12/2013
//-----------------------------------------------------------------------------
bool AnalysisInterface::setLinSol(const N_UTL_OptionBlock & OB)
{
  return anaManagerPtr_->setLinSol(OB);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::setHBLinSol
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist
// Creation Date : 11/11/08
//-----------------------------------------------------------------------------
bool AnalysisInterface::setHBLinSol(const N_UTL_OptionBlock & OB)
{
  return anaManagerPtr_->setHBLinSol(OB);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::setMORAnalysisParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist and Ting Mei
// Creation Date : 5/31/12
//-----------------------------------------------------------------------------
bool AnalysisInterface::setMORAnalysisParams(const N_UTL_OptionBlock & OB)
{
  return anaManagerPtr_->setMORAnalysisParams(OB);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::setMOROptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist and Ting Mei
// Creation Date : 5/31/12
//-----------------------------------------------------------------------------
bool AnalysisInterface::setMOROptions(const N_UTL_OptionBlock & OB)
{
  return anaManagerPtr_->setMOROptions(OB);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::setTRANMPDEOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 7/23/08
//-----------------------------------------------------------------------------
bool AnalysisInterface::setTRANMPDEOptions(const N_UTL_OptionBlock & OB)
{
  return anaManagerPtr_->setTRANMPDEOptions(OB);
}


//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::setSensOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, Sandia
// Creation Date : 6/5/13
//-----------------------------------------------------------------------------
bool AnalysisInterface::setSensOptions(const N_UTL_OptionBlock & OB)
{
  return anaManagerPtr_->setSensOptions(OB);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::registerElapsedTimer
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 02/24/06
//-----------------------------------------------------------------------------
bool AnalysisInterface::registerElapsedTimer (N_UTL_Timer * et)
{
  return anaManagerPtr_->registerElapsedTimer(et);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::getBlockAnalysisFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/29/05
//-----------------------------------------------------------------------------
bool AnalysisInterface::getBlockAnalysisFlag () const
{
  return anaManagerPtr_->getBlockAnalysisFlag();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::getTransientFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/29/05
//-----------------------------------------------------------------------------
bool AnalysisInterface::getTransientFlag () const
{
  return anaManagerPtr_->getTransientFlag();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::initializeAll
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

bool AnalysisInterface::initializeAll()
{
  return anaManagerPtr_->initializeAll();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::registerLinearSystem
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/05/01
//-----------------------------------------------------------------------------
bool AnalysisInterface::registerLinearSystem(N_LAS_System * lasSysPtr)
{
  anaManagerPtr_->registerLinearSystem(lasSysPtr);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::registerLoader
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/05/01
//-----------------------------------------------------------------------------
bool AnalysisInterface::registerLoader(N_LOA_Loader * loaderPtr)
{
  return anaManagerPtr_->registerLoader(loaderPtr);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::registerOutputMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/05/01
//-----------------------------------------------------------------------------
bool AnalysisInterface::registerOutputMgr(N_IO_OutputMgr * outputPtr)
{
  return anaManagerPtr_->registerOutputMgr(outputPtr);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::registerRestartMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/05/01
//-----------------------------------------------------------------------------
bool AnalysisInterface::registerRestartMgr(N_IO_RestartMgr * restartPtr)
{
  return anaManagerPtr_->registerRestartMgr(restartPtr);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::registerNLSolver
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/05/01
//-----------------------------------------------------------------------------
bool AnalysisInterface::registerNLSManager(N_NLS_Manager * nlsMgrPtr)
{
  return anaManagerPtr_->registerNLSManager(nlsMgrPtr);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::run
// Purpose       : Execute the top level control loop.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, 9233, Computational Sciences
// Creation Date : 11/10/00
//-----------------------------------------------------------------------------

bool AnalysisInterface::run()
{
  return anaManagerPtr_->run();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::updateDerivs
// Purpose       : calls the time  int. method to update the  corrector
//                 derivatives.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/10/00
//-----------------------------------------------------------------------------

bool AnalysisInterface::updateDerivs()
{
  return anaManagerPtr_->updateDerivs ();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::updateDivDiffs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/10/00
//-----------------------------------------------------------------------------

bool AnalysisInterface::updateDivDiffs()
{
  return anaManagerPtr_->updateDivDiffs();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::getTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/27/00
//-----------------------------------------------------------------------------
double AnalysisInterface::getTime() const
{
  return anaManagerPtr_->getTime();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::getCurrentTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Mems Modeling
// Creation Date : 8/21/2009
//-----------------------------------------------------------------------------
double AnalysisInterface::getCurrentTime() const
{
  return anaManagerPtr_->getCurrentTime();
}


//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::getFinalTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 9/04/01
//-----------------------------------------------------------------------------
double AnalysisInterface::getFinalTime() const
{
  return anaManagerPtr_->getFinalTime();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::getInitialTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/04/02
//-----------------------------------------------------------------------------
double AnalysisInterface::getInitialTime() const
{
  return anaManagerPtr_->getInitialTime();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::equateTmpVectors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/24/02
//-----------------------------------------------------------------------------
bool AnalysisInterface::equateTmpVectors()
{
  return anaManagerPtr_->equateTmpVectors();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::updateDerivsBlock
// Purpose       : calls the time  int. method to update the  corrector
//                 derivatives, but only does it a subset of the solution
//                 and state vectors.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/10/01
//-----------------------------------------------------------------------------
//inline
bool AnalysisInterface::updateDerivsBlock(const std::list<index_pair> & solGIDList,
                                      const std::list<index_pair> & staGIDList)
{
  return anaManagerPtr_->updateDerivsBlock(solGIDList, staGIDList);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::restartDataSize
// Purpose       : Gets the size of the restart data (bytes?).
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
int AnalysisInterface::restartDataSize( bool pack )
{
  return anaManagerPtr_->restartDataSize( pack );
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::dumpRestartData
// Purpose       : Dumps the restart data to a file.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool AnalysisInterface::dumpRestartData(char * buf, int bsize, int & pos,
                                    N_PDS_Comm * comm, bool pack )
{
  return anaManagerPtr_->dumpRestartData(buf, bsize, pos, comm, pack);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::restoreRestartData
// Purpose       : Restores the restart data from a file.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool AnalysisInterface::restoreRestartData(char * buf, int bsize, int & pos,
                                       N_PDS_Comm * comm, bool pack )
{
  return anaManagerPtr_->restoreRestartData(buf, bsize, pos, comm, pack);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::getSolnVarData
// Purpose       : Gets the solution variable data.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool AnalysisInterface::getSolnVarData(const int & gid, std::vector<double> & varData)
{
  return anaManagerPtr_->getSolnVarData(gid, varData);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::getStateVarData
// Purpose       : Gets the state variable data.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool AnalysisInterface::getStateVarData(const int & gid, std::vector<double> & varData)
{
  return anaManagerPtr_->getStateVarData(gid, varData);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::getStoreVarData
// Purpose       : Gets the store variable data.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
bool AnalysisInterface::getStoreVarData(const int & gid, std::vector<double> & varData)
{
  return anaManagerPtr_->getStoreVarData(gid, varData);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::setSolnVarData
// Purpose       : Sets the solution variable data.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool AnalysisInterface::setSolnVarData(const int & gid,
                                   const std::vector<double> & varData)
{
  return anaManagerPtr_->setSolnVarData(gid, varData);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::setStateVarData
// Purpose       : Sets the state variable data.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 07/31/01
//-----------------------------------------------------------------------------
bool AnalysisInterface::setStateVarData(const int & gid,
                                    const std::vector<double> & varData)
{
  return anaManagerPtr_->setStateVarData(gid, varData);
}


//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::setStoreVarData
// Purpose       : Sets the state variable data.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool AnalysisInterface::setStoreVarData(const int & gid,
                                    const std::vector<double> & varData)
{
  return anaManagerPtr_->setStoreVarData(gid, varData);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::setBeginningIntegrationFlag
// Purpose       : set beginning integration flag to true
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytems Modeling
// Creation Date : 03/5/2010
//-----------------------------------------------------------------------------
void AnalysisInterface::setBeginningIntegrationFlag(bool bif)
{
  anaManagerPtr_->setBeginningIntegrationFlag (bif);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::registerParallelServices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/30/01
//-----------------------------------------------------------------------------
bool AnalysisInterface::registerParallelServices(N_PDS_Manager * pds)
{
  return anaManagerPtr_->registerParallelServices(pds);

}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::outputSummary
// Purpose       : Outputs summary information for time-integration function.
// Special Notes :
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 07/30/02
//-----------------------------------------------------------------------------
bool AnalysisInterface::outputSummary()
{
  bool bsuccess = true;

  if (!Teuchos::is_null(anaManagerPtr_))
  {
    bsuccess = anaManagerPtr_->printLoopInfo(0,0);
  }

  return bsuccess;

}


//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::setNextSolVectorPtr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 04/22/03
//-----------------------------------------------------------------------------
bool AnalysisInterface::setNextSolVectorPtr (N_LAS_Vector * solVecPtr)
{
  return anaManagerPtr_->setNextSolVectorPtr (solVecPtr);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::loadRHS
// Purpose       : manages the RHS load.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 03/03/04
//-----------------------------------------------------------------------------
bool AnalysisInterface::loadRHS ()
{
  return anaManagerPtr_->assemblerPtr->loadRHS();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::loadJacobian
// Purpose       : Manages the jacobian load.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 03/03/04
//-----------------------------------------------------------------------------
bool AnalysisInterface::loadJacobian ()
{
  return anaManagerPtr_->assemblerPtr->loadJacobian();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::applyJacobian
// Purpose       : Manages the jacobian load.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 03/03/04
//-----------------------------------------------------------------------------
bool AnalysisInterface::applyJacobian (const N_LAS_Vector& input, N_LAS_Vector& result)
{
  return anaManagerPtr_->assemblerPtr->applyJacobian(input,result);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::getResidualTime
// Purpose       : This should give the time for the most recent Residual load.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 03/03/04
//-----------------------------------------------------------------------------
double AnalysisInterface::getResidualTime()
{
  return anaManagerPtr_->assemblerPtr->residualTime_;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::getJacobianTime
// Purpose       : This should give the time for the most recent Jacobian load.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 03/03/04
//-----------------------------------------------------------------------------
double AnalysisInterface::getJacobianTime()
{
  return anaManagerPtr_->assemblerPtr->jacobianTime_;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::getTIAParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 11/05/04
//-----------------------------------------------------------------------------
N_TIA_TIAParams& AnalysisInterface::getTIAParams()
{
  return anaManagerPtr_->tiaParams;
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::simulationComplete
// Purpose       : return boolean signifying whether simulation is complete
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 05/06/2004
//-----------------------------------------------------------------------------
bool AnalysisInterface::simulationComplete()
{
  return anaManagerPtr_->simulationComplete();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::completeOPStartStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 06/28/2006
//-----------------------------------------------------------------------------
bool AnalysisInterface::completeOPStartStep()
{
  return anaManagerPtr_->completeOPStartStep();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::completeHomotopyStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, 9233
// Creation Date : 03/20/2006
//-----------------------------------------------------------------------------
bool AnalysisInterface::completeHomotopyStep
    ( const std::vector<std::string> & paramNames,
      const std::vector<double> & paramVals,
      N_LAS_Vector * solnVecPtr )
{
  return anaManagerPtr_->completeHomotopyStep(paramNames, paramVals,solnVecPtr);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::failHomotopyStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, 9233
// Creation Date : 03/30/2006
//-----------------------------------------------------------------------------
bool AnalysisInterface::failHomotopyStep ( )
{
  return anaManagerPtr_->failHomotopyStep ();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::startTimeStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/20/2006
//-----------------------------------------------------------------------------
bool AnalysisInterface::startTimeStep (const N_TIA_TimeIntInfo & tiInfo)
{
  return anaManagerPtr_->startTimeStep (tiInfo);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::runStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/06/2006
//-----------------------------------------------------------------------------
bool AnalysisInterface::runStep
    (const N_TIA_TimeIntInfo & tiInfo, N_TIA_TwoLevelError & tlError)
{
  return anaManagerPtr_->runStep (tiInfo, tlError);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::startupSolvers
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/06/2006
//-----------------------------------------------------------------------------
bool AnalysisInterface::startupSolvers ()
{
  return anaManagerPtr_->startupSolvers ();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::finishSolvers
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/06/2006
//-----------------------------------------------------------------------------
bool AnalysisInterface::finishSolvers ()
{
  return anaManagerPtr_->finishSolvers ();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::provisionalStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/04/09
//-----------------------------------------------------------------------------
bool AnalysisInterface::provisionalStep (double maxTimeStep,
      double &currTimeStep)
{
  return anaManagerPtr_->provisionalStep (maxTimeStep, currTimeStep);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::acceptProvisionalStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 04/10/09
//-----------------------------------------------------------------------------
void AnalysisInterface::acceptProvisionalStep ()
{
  anaManagerPtr_->acceptProvisionalStep ();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::rejectProvisionalStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 04/10/09
//-----------------------------------------------------------------------------
void AnalysisInterface::rejectProvisionalStep ()
{
  anaManagerPtr_->rejectProvisionalStep ();
}

#ifdef Xyce_Dakota

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::getDakotaRunFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool AnalysisInterface::getDakotaRunFlag()
{
  return anaManagerPtr_->getDakotaRunFlag();
}


//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::setDakotaRunFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL
// Creation Date :
//-----------------------------------------------------------------------------
void AnalysisInterface::setDakotaRunFlag( bool flag )
{
  anaManagerPtr_->setDakotaRunFlag( flag );
}


//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::getDakotaIteration
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL
// Creation Date :
//-----------------------------------------------------------------------------
int AnalysisInterface::getDakotaIteration()
{
    return anaManagerPtr_->getDakotaIteration();
}


//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::setDakotaIteration
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL
// Creation Date :
//-----------------------------------------------------------------------------
void AnalysisInterface::setDakotaIteration( int iterNumber )
{
  anaManagerPtr_->setDakotaIteration( iterNumber );
}

#endif

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::condutanceTest
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
void AnalysisInterface::condutanceTest ()
{
  anaManagerPtr_->conductanceTest ();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::homotopyStepSuccess
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
void AnalysisInterface::homotopyStepSuccess
    ( const std::vector<std::string> & paramNames,
      const std::vector<double> & paramVals)
{
  anaManagerPtr_->homotopyStepSuccess ( paramNames, paramVals);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::homotopyStepFailure
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/30/06
//-----------------------------------------------------------------------------
void AnalysisInterface::homotopyStepFailure ()
{
  anaManagerPtr_->homotopyStepFailure ();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::stepSuccess
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
void AnalysisInterface::stepSuccess (int analysis)
{
  anaManagerPtr_->stepSuccess (analysis);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::stepFailure
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
void AnalysisInterface::stepFailure (int analysis)
{
  anaManagerPtr_->stepFailure (analysis);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::getInitialQnorm
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
bool AnalysisInterface::getInitialQnorm (N_TIA_TwoLevelError & tle)
{
  return anaManagerPtr_->getInitialQnorm (tle);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::getBreakPoints
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
bool AnalysisInterface::getBreakPoints (std::vector<N_UTL_BreakPoint> &breakPointTimes)
{
  return anaManagerPtr_->getBreakPoints (breakPointTimes);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::getTimeIntInfo
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 02/11/07
//-----------------------------------------------------------------------------
void AnalysisInterface::getTimeIntInfo(N_TIA_TimeIntInfo & tiInfo)
{
  anaManagerPtr_->getTimeIntInfo (tiInfo);
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::silenceProgress
// Purpose       : shut up "percent complete" noise
// Special Notes :
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 11 Feb 2009
//-----------------------------------------------------------------------------
void AnalysisInterface::silenceProgress()
{
  anaManagerPtr_->silenceProgress();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::enableProgress
// Purpose       : enable "percent complete" noise
// Special Notes :
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 11 Feb 2009
//-----------------------------------------------------------------------------
void AnalysisInterface::enableProgress()
{
  anaManagerPtr_->enableProgress();
}

//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::initializeTransientModel
// Purpose       : ModelEvaluator Interface
// Special Notes :
// Scope         : public
// Creator       : Coffey, Schiek, Mei
// Creation Date : 05/29/09
//-----------------------------------------------------------------------------
void AnalysisInterface::initializeTransientModel()
{
  anaManagerPtr_->initializeTransientModel();
}
//-----------------------------------------------------------------------------
// Function      : AnalysisInterface::evalModel
// Purpose       : ModelEvaluator Interface
// Special Notes :
// Scope         : public
// Creator       : Coffey, Schiek, Mei
// Creation Date : 05/27/09
//-----------------------------------------------------------------------------
bool AnalysisInterface::evalTransientModel(
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
// Function      : AnalysisInterface::evalTransientModelState
// Purpose       : ModelEvaluator Interface
// Special Notes :
// Scope         : public
// Creator       : Coffey, Schiek, Mei
// Creation Date : 05/29/09
//-----------------------------------------------------------------------------
bool AnalysisInterface::evalTransientModelState(
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

} // namespace Analysis
} // namespace Xyce

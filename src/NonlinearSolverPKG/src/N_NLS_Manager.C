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
// Filename       : $RCSfile: N_NLS_Manager.C,v $
//
// Purpose        : Body for the nonlinear solver manager class implementation
//                  which will allow for the creation and selection of NLS
//                  algorithms.
//
// Special Notes  : GOF Strategy Pattern
//
// Creator        : Scott A. Hutchinson, SNL, Parallel Computational Sciences
//
// Creation Date  : 04/28/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.169.2.4 $
//
// Revision Date  : $Date $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ----------  Xyce Includes   ----------
#include <N_UTL_Misc.h>

#include <N_NLS_Manager.h>
#include <N_NLS_NOX_Interface.h>
#include <N_NLS_DampedNewton.h>
#include <N_NLS_TwoLevelNewton.h>
#include <N_NLS_Sensitivity.h>
#include <N_NLS_ConductanceExtractor.h>
#include <N_NLS_NonLinInfo.h>

#include <N_ERH_ErrorMgr.h>
#include <N_IO_CmdParse.h>
#include <N_IO_OutputMgr.h>
#include <N_UTL_Expression.h>

// ----------  Trilinos Includes   ----------

#include <Teuchos_RefCountPtr.hpp>

//-----------------------------------------------------------------------------
// Function      : N_NLS_Manager::N_NLS_Manager
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/01/00
//-----------------------------------------------------------------------------
N_NLS_Manager::N_NLS_Manager(N_IO_CmdParse & cp)
  : commandLine_(cp),
    twoLevelNewtonFlag_(false),
    noxFlag_(false),
    noxFlagInner_(false),
    noxFlagTransient_(false),
    anaIntPtr_(0x0),
    loaderPtr_(0x0),
    lasSysPtr_(0x0),
    rhsVecPtr_(0x0),
    setupSensFlag_(false),
    nlsSensitivityPtr_(0x0),
    outputPtr_(0x0),
    pdsMgrPtr_(0x0),
    nlsPtr_(0x0),
    conductanceExtractorPtr_(0x0),
    initializeAllFlag_(false),
    exprPtr(0x0)
{
  matrixFreeFlag_ = false;
}


//-----------------------------------------------------------------------------
// Function      : N_NLS_Manager::~N_NLS_Manager
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/01/00
//-----------------------------------------------------------------------------

N_NLS_Manager::~N_NLS_Manager()
{
  delete nlsPtr_;
  nlsPtr_ = 0x0;

  if (nlsSensitivityPtr_)
  {
    delete nlsSensitivityPtr_;
    nlsSensitivityPtr_ = 0x0;
  }

  if (conductanceExtractorPtr_)
  {
    delete conductanceExtractorPtr_;
    conductanceExtractorPtr_  = 0x0;
  }
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Manager::registerPkgOptionsMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, 1437
// Creation Date : 10/21/08
//-----------------------------------------------------------------------------
bool N_NLS_Manager::registerPkgOptionsMgr( RCP<N_IO_PkgOptionsMgr> pkgOptPtr )
{
  pkgOptMgrPtr_ = pkgOptPtr;
  string netListFile = "";
  if (commandLine_.getArgumentValue("netlist") != "")
  {
    netListFile = commandLine_.getArgumentValue("netlist");
  }
  pkgOptMgrPtr_->submitRegistration(
      "NONLIN", netListFile, new N_NLS_Manager_OptionsReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "NONLIN-TRAN", netListFile, new N_NLS_Manager_TranOptionsReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "NONLIN-HB", netListFile, new N_NLS_Manager_HBOptionsReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "LINSOL", netListFile, new N_NLS_Manager_LSOptionsReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "LOCA", netListFile, new N_NLS_Manager_LocaOptionsReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "SENS", netListFile, new N_NLS_Manager_SensOptionsReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "SENSITIVITY", netListFile, new N_NLS_Manager_SensitivityOptionsReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "NONLIN-TWOLEVEL", netListFile, new N_NLS_Manager_TwoLvlOptionsReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "NONLIN-TWOLEVEL-TRAN", netListFile, new N_NLS_Manager_TwoLvlTranOptionsReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "TIMEINT", netListFile, new N_NLS_Manager_TimeOptionsReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "DCOP", netListFile, new N_NLS_Manager_DCOPRestartOptionsReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "IC", netListFile, new N_NLS_Manager_ICOptionsReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "NODESET", netListFile, new N_NLS_Manager_NodeSetOptionsReg( this ) );
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Manager::registerParallelMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, Sandia
// Creation Date : 
//-----------------------------------------------------------------------------
bool N_NLS_Manager::registerParallelMgr (N_PDS_Manager * pdsMgrPtr )
{
  pdsMgrPtr_ = pdsMgrPtr;
  return (pdsMgrPtr_ != 0x0);
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Manager::setOptions
// Purpose       : Sets the nonlinear solver options.
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 9/29/00
//-----------------------------------------------------------------------------
bool N_NLS_Manager::setOptions(const N_UTL_OptionBlock & OB)
{
  optionBlockMap_["dcop"] = OB;
  outputPtr_->setEnableHomotopyFlag(true);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Manager::setTranOptions
// Purpose       : Sets the nonlinear solver options.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 9/05/01
//-----------------------------------------------------------------------------
bool N_NLS_Manager::setTranOptions(const N_UTL_OptionBlock & OB)
{
  optionBlockMap_["transient"] = OB;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Manager::setHBOptions
// Purpose       : Sets the nonlinear solver options.
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 02/03/2009
//-----------------------------------------------------------------------------
bool N_NLS_Manager::setHBOptions(const N_UTL_OptionBlock & OB)
{
  optionBlockMap_["hb"] = OB;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Manager::getHBOptions
// Purpose       : Sets the nonlinear solver options.
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 02/03/2009
//-----------------------------------------------------------------------------
bool N_NLS_Manager::getHBOptions(N_UTL_OptionBlock & HBOB)
{
  HBOB = optionBlockMap_["hb"];
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Manager::setLinSolOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 11/9/00
//-----------------------------------------------------------------------------
bool N_NLS_Manager::setLinSolOptions(const N_UTL_OptionBlock & OB)
{
  optionBlockMap_["petra"] = OB;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Manager::setLocaOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/02/03
//-----------------------------------------------------------------------------
bool N_NLS_Manager::setLocaOptions(const N_UTL_OptionBlock & OB)
{
  optionBlockMap_["loca"] = OB;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Manager::setTwoLevelLocaOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/02/03
//-----------------------------------------------------------------------------
bool N_NLS_Manager::setTwoLevelLocaOptions(const N_UTL_OptionBlock & OB)
{
  optionBlockMap_["twolevelloca"] = OB;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Manager::setTwoLevelOptions
// Purpose       : This option setter will create a new class structure.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/20/02
//-----------------------------------------------------------------------------
bool N_NLS_Manager::setTwoLevelOptions(const N_UTL_OptionBlock & OB)
{
  optionBlockMap_["twolevel"] = OB;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Manager::setTwoLevelTranOptions
// Purpose       : This option setter will create a new class structure.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/20/02
//-----------------------------------------------------------------------------
bool N_NLS_Manager::setTwoLevelTranOptions(const N_UTL_OptionBlock & OB)
{
  optionBlockMap_["twoleveltran"] = OB;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Manager::registerRHSVector
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/01/00
//-----------------------------------------------------------------------------
bool N_NLS_Manager::registerRHSVector(N_LAS_Vector* tmp_RHSVecPtr)
{
  rhsVecPtr_ = tmp_RHSVecPtr;
  return (rhsVecPtr_ != 0x0);
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Manager::registerLoader
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/01/00
//-----------------------------------------------------------------------------
bool N_NLS_Manager::registerLoader(N_LOA_Loader* tmp_LoaderPtr)
{
  loaderPtr_ = tmp_LoaderPtr;
  return (loaderPtr_ != 0x0);
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Manager::registerOutputMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 9/23/03
//-----------------------------------------------------------------------------
bool N_NLS_Manager::registerOutputMgr(N_IO_OutputMgr * outputPtr)
{
  outputPtr_ = outputPtr;
  return (outputPtr_ != 0x0);
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Manager::registerLinearSystem
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/09/00
//-----------------------------------------------------------------------------
bool N_NLS_Manager::registerLinearSystem(N_LAS_System* tmp_LasSysPtr)
{
  lasSysPtr_ = tmp_LasSysPtr;
  return (lasSysPtr_ != 0x0);
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Manager::registerPrecondFactory
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, SNL, Electrical & Microsystem Modeling
// Creation Date : 11/11/08
//-----------------------------------------------------------------------------
bool N_NLS_Manager::registerPrecondFactory(const RCP<N_LAS_PrecondFactory>& tmp_LasPrecPtr)
{
  lasPrecPtr_ = tmp_LasPrecPtr;
  return (!Teuchos::is_null(lasPrecPtr_));
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Manager::registerAnalysisInterface
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/03/00
//-----------------------------------------------------------------------------
bool N_NLS_Manager::registerAnalysisInterface(N_ANP_AnalysisInterface* tmp_anaIntPtr)
{
  anaIntPtr_ = tmp_anaIntPtr;
  return (anaIntPtr_ != 0x0);
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Manager::using Nox_
// Purpose       : This function determines if we are using NOX or not.
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/30/02
//-----------------------------------------------------------------------------
void N_NLS_Manager::usingNox_ ()
{
  noxFlag_ = true;
  noxFlagInner_ = true;
  noxFlagTransient_ = true;

  N_UTL_OptionBlock OBdcop = optionBlockMap_["dcop"];
  N_UTL_OptionBlock OBtran = optionBlockMap_["transient"];

  // scan for changes to the dcop values of noxFlag_
  for (list<N_UTL_Param>::const_iterator it_tpL = OBdcop.getParams().begin();
       it_tpL != OBdcop.getParams().end(); ++ it_tpL)
  {
    if (it_tpL->uTag() == "NOX")
    {
      noxFlag_ = it_tpL->iVal();
      noxFlagInner_ = noxFlag_;
    }
  }

  // now check for a nox flag on the .options twolevel line, if it exists.
  map<string,N_UTL_OptionBlock>::iterator obmIter = optionBlockMap_.find("twolevel") ;
  if ( obmIter != optionBlockMap_.end() )
  {
    //N_UTL_OptionBlock OBtwoLevel = optionBlockMap_["twolevel"];
    N_UTL_OptionBlock OBtwoLevel = obmIter->second;
    for (list<N_UTL_Param>::const_iterator it_tpL = OBtwoLevel.getParams().begin();
	 it_tpL != OBtwoLevel.getParams().end(); ++ it_tpL)
    {
      if (it_tpL->uTag() == "NOX")
      {
        noxFlagInner_ = it_tpL->iVal();
      }
    }
  }


  // scan for changes to the transient value for noxFlagTransient_
  for (list<N_UTL_Param>::const_iterator it_tpL = OBtran.getParams().begin();
       it_tpL != OBtran.getParams().end(); ++ it_tpL)
  {
    if (it_tpL->uTag() == "NOX")
    {
      noxFlagTransient_ = it_tpL->iVal();
    }
  }

  obmIter = optionBlockMap_.find("hb") ;
  if ( obmIter != optionBlockMap_.end() )
  {
    N_UTL_OptionBlock OBhb= obmIter->second;
    for (list<N_UTL_Param>::const_iterator it_tpL = OBhb.getParams().begin();
	 it_tpL != OBhb.getParams().end(); ++ it_tpL)
    {
      if (it_tpL->uTag() == "NOX")
      {
        noxFlag_ = it_tpL->iVal();
        noxFlagInner_ = noxFlag_;
      }
    }
  }

  // now check if the command line has specified nox.  The command line
  // overrides the netlist.
  if( commandLine_.getArgumentValue( "-nox" ) == "off" )
  {
    noxFlag_ = false;
    noxFlagInner_ = false;
    noxFlagTransient_ = false;
  }

#ifdef Xyce_DEBUG_NONLINEAR
  string msg = "noxFlag is: ";
  if (noxFlag_) msg += "TRUE";
  else          msg += "FALSE";
  msg += "\n";
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, msg);
  string msg2 = "noxFlagTransient is: ";
  if (noxFlagTransient_) msg2 += "TRUE";
  else          msg2 += "FALSE";
  msg2 += "\n";
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, msg2);
#endif

}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Manager::allocateSolver_
// Purpose       : This function determines which solver to allocate, and
//                 allocates it.
//
//                 Right now the possibilities are:
//
//                    N_NLS_DampedNewton
//                    N_NLS_NOXInterface
//                    N_NLS_TwoLevelNewton
//
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/28/02
//-----------------------------------------------------------------------------
bool N_NLS_Manager::allocateSolver_ ()
{
  bool bsuccess = true;
  bool bs1 = true;

  // determine if we are using NOX or not.
  usingNox_ ();

  // if ".options nonlin two-level" appeared in the netlist, then
  // allocate the two-level solver.  Otherwise allocate one of the single
  // level solvers.

  if ( (optionBlockMap_.find(    "twolevel") != optionBlockMap_.end()) ||
       (optionBlockMap_.find("twoleveltran") != optionBlockMap_.end())
     )
  {

    twoLevelNewtonFlag_ = true;
    nlsPtr_ = new N_NLS_TwoLevelNewton (noxFlag_, noxFlagInner_, commandLine_);
  }
  else
  {
    twoLevelNewtonFlag_ = false;
    if (noxFlag_)
    {
      if( nlsPtr_ != 0x0)
        delete nlsPtr_;
      nlsPtr_ = new N_NLS_NOX::Interface(commandLine_);
    }
    else
    {
      if( nlsPtr_ != 0x0)
        delete nlsPtr_;
      nlsPtr_ = new N_NLS_DampedNewton(commandLine_);
    }
  }

  // now register everything, now that the solver class is set up.
  bs1 = nlsPtr_->registerLinearSystem(lasSysPtr_);      bsuccess = bsuccess && bs1;
  bs1 = nlsPtr_->registerAnalysisInterface(anaIntPtr_); bsuccess = bsuccess && bs1;
  bs1 = nlsPtr_->registerLoader(loaderPtr_);            bsuccess = bsuccess && bs1;

  if (!Teuchos::is_null(lasPrecPtr_))
    bs1 = nlsPtr_->registerPrecondFactory(lasPrecPtr_); bsuccess = bsuccess && bs1;

  bs1 = nlsPtr_->registerOutputMgr (outputPtr_); bsuccess = bsuccess && bs1;
  bs1 = nlsPtr_->registerParallelMgr (pdsMgrPtr_); bsuccess = bsuccess && bs1;

  map<string,N_UTL_OptionBlock>::iterator obmIter = optionBlockMap_.find("dcop") ;

  if ( obmIter != optionBlockMap_.end() )
  {
    const N_UTL_OptionBlock OB = obmIter->second;
    bs1 = nlsPtr_->setOptions (OB);
    bsuccess = bsuccess && bs1;
  }

  obmIter = optionBlockMap_.find("transient") ;
  if ( obmIter != optionBlockMap_.end() )
  {
    const N_UTL_OptionBlock OB = obmIter->second;
    bs1 = nlsPtr_->setTranOptions (OB);
    bsuccess = bsuccess && bs1;
  }

  obmIter = optionBlockMap_.find("hb") ;
  if ( obmIter != optionBlockMap_.end() )
  {
    const N_UTL_OptionBlock OB = obmIter->second;
    bs1 = nlsPtr_->setHBOptions (OB);
    bsuccess = bsuccess && bs1;
  }

  obmIter = optionBlockMap_.find("loca") ;
  if ( obmIter != optionBlockMap_.end() )
  {
    const N_UTL_OptionBlock OB = obmIter->second;
    bs1 = nlsPtr_->setLocaOptions (OB);
    bsuccess = bsuccess && bs1;
  }

  obmIter = optionBlockMap_.find("twolevel") ;
  if ( obmIter != optionBlockMap_.end() )
  {
    const N_UTL_OptionBlock OB = obmIter->second;
    bs1 = nlsPtr_->setTwoLevelOptions (OB);
    bsuccess = bsuccess && bs1;
  }

  obmIter = optionBlockMap_.find("twoleveltran") ;
  if ( obmIter != optionBlockMap_.end() )
  {
    const N_UTL_OptionBlock OB = obmIter->second;
    bs1 = nlsPtr_->setTwoLevelTranOptions (OB);
    bsuccess = bsuccess && bs1;
  }

  obmIter = optionBlockMap_.find("twolevelloca") ;
  if ( obmIter != optionBlockMap_.end() )
  {
    const N_UTL_OptionBlock OB = obmIter->second;
    bs1 = nlsPtr_->setTwoLevelLocaOptions (OB);
    bsuccess = bsuccess && bs1;
  }

  obmIter = optionBlockMap_.find("petra") ;
  if ( obmIter != optionBlockMap_.end() )
  {
    const N_UTL_OptionBlock OB = obmIter->second;
    bs1 = nlsPtr_->setPetraOptions(OB);
    bsuccess = bsuccess && bs1;
  }

  obmIter = optionBlockMap_.find("dcop_restart") ;
  if ( obmIter != optionBlockMap_.end() )
  {
    const N_UTL_OptionBlock OB = obmIter->second;
    bs1 = nlsPtr_->setDCOPRestartOptions(OB);
    bsuccess = bsuccess && bs1;
  }

  obmIter = optionBlockMap_.find("ic") ;
  if ( obmIter != optionBlockMap_.end() )
  {
    const N_UTL_OptionBlock OB = obmIter->second;
    bs1 = nlsPtr_->setICOptions(OB);
    bsuccess = bsuccess && bs1;
  }

  obmIter = optionBlockMap_.find("nodeset") ;
  if ( obmIter != optionBlockMap_.end() )
  {
    const N_UTL_OptionBlock OB = obmIter->second;
    bs1 = nlsPtr_->setNodeSetOptions(OB);
    bsuccess = bsuccess && bs1;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Manager::initializeAll
// Purpose       : This can only be called after the linear system
//                 (N_LAS_System) has been registered, and after all the
//                 options have been set.
//
// Special Notes : This function obtains the solution, temporary solution and
//                 rhs vectors from the LAS system class.
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/12/00
//-----------------------------------------------------------------------------
bool N_NLS_Manager::initializeAll()
{
  bool bsuccess = true;
  bool bs1 = true;
  bool bs2 = true;

  bs1 = allocateSolver_ ();        bsuccess = bsuccess && bs1;
  nlsPtr_->setMatrixFreeFlag(matrixFreeFlag_);
  bs2 = nlsPtr_->initializeAll();  bsuccess = bsuccess && bs2;

  nlsPtr_->setReturnCodes (retCodes_);

  initializeAllFlag_ = true;

  if(conductanceExtractorPtr_ == 0x0)
  {
    conductanceExtractorPtr_ = new N_NLS_ConductanceExtractor (
      *nlsPtr_, *topPtr_, commandLine_);
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Manager::allocateTranSolver
// Purpose       : Allocate a different solver for transient problems if
//                 the user has requested a different one.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 12/22/08
//-----------------------------------------------------------------------------
void N_NLS_Manager::allocateTranSolver()
{
  bool bsuccess = true;
  bool bs1 = true;

  // only do a reallocation if the the solver type is changing from dcop to transient
  if( ((noxFlag_ == true) && (noxFlagTransient_ == false)) ||
      ((noxFlag_ == false) && (noxFlagTransient_ == true)) )
  {
    if (noxFlagTransient_)
    {
      if( nlsPtr_ != 0x0)
        delete nlsPtr_;
      nlsPtr_ = new N_NLS_NOX::Interface(commandLine_);
    }
    else
    {
      if( nlsPtr_ != 0x0)
        delete nlsPtr_;
      nlsPtr_ = new N_NLS_DampedNewton(commandLine_);
    }

    // now register everything, now that the solver class is set up.
    nlsPtr_->registerLinearSystem(lasSysPtr_);
    nlsPtr_->registerAnalysisInterface(anaIntPtr_);
    nlsPtr_->registerLoader(loaderPtr_);
    nlsPtr_->registerOutputMgr (outputPtr_);
    nlsPtr_->registerParallelMgr (pdsMgrPtr_); 

    nlsPtr_->setMatrixFreeFlag(matrixFreeFlag_);
    nlsPtr_->initializeAll();
    nlsPtr_->setReturnCodes (retCodes_);

    map<string,N_UTL_OptionBlock>::iterator obmIter = optionBlockMap_.find("transient") ;
    if ( obmIter != optionBlockMap_.end() )
    {
      //  const N_UTL_OptionBlock OB = optionBlockMap_["transient"];
      const N_UTL_OptionBlock OB = obmIter->second;
      bs1 = nlsPtr_->setTranOptions (OB);
      bsuccess = bsuccess && bs1;
    }

    initializeAllFlag_ = true;
    delete conductanceExtractorPtr_;
    conductanceExtractorPtr_ = new N_NLS_ConductanceExtractor (*nlsPtr_, *topPtr_, commandLine_);
  }
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Manager::isFirstContinuationParam
// Purpose       : This function returns true if is the first LOCA
//                 continuation param
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 7/06/05
//-----------------------------------------------------------------------------
bool N_NLS_Manager::isFirstContinuationParam()
{
  return nlsPtr_->isFirstContinuationParam ();
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Manager::isFirstSolveComplete
// Purpose       : This function incdicates if the initial solve has completed
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 6/22/06
//-----------------------------------------------------------------------------
bool N_NLS_Manager::isFirstSolveComplete()
{
  return nlsPtr_->isFirstSolveComplete ();
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Manager::getContinuationStep
// Purpose       : This function returns the value of the current LOCA
//                 continuation step.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 2/28/05
//-----------------------------------------------------------------------------
int N_NLS_Manager::getContinuationStep ()
{
  return nlsPtr_->getContinuationStep ();
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Manager::getLocaFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/09/05
//-----------------------------------------------------------------------------
bool N_NLS_Manager::getLocaFlag ()
{
  return  nlsPtr_->getLocaFlag ();
}
//-----------------------------------------------------------------------------
// Function      : N_NLS_Manager::getNumIterations
// Purpose       : This function returns the value of the current Newton
//                 iteration.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 2/23/01
//-----------------------------------------------------------------------------
int N_NLS_Manager::getNumIterations()
{
  return nlsPtr_->getNumIterations();
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Manager::getCouplingMode
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, Parallel Computational Sciences
// Creation Date : 12/05/02
//-----------------------------------------------------------------------------
int N_NLS_Manager::getCouplingMode ()
{
  return nlsPtr_->getCouplingMode ();
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Manager::getNonLinInfo
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, Parallel Computational Sciences
// Creation Date : 12/05/02
//-----------------------------------------------------------------------------
void N_NLS_Manager::getNonLinInfo (N_NLS_NonLinInfo & nlInfo)
{
  nlInfo.newtonIter    = nlsPtr_->getNumIterations();
  nlInfo.twoLevelNewtonCouplingMode  = nlsPtr_->getCouplingMode ();
  nlInfo.locaFlag      = nlsPtr_->getLocaFlag ();

  if (nlInfo.locaFlag)
  {
    nlInfo.continuationStep       = nlsPtr_->getContinuationStep ();
    nlInfo.firstContinuationParam = nlsPtr_->isFirstContinuationParam ();
    nlInfo.firstSolveComplete     = nlsPtr_->isFirstSolveComplete ();
  }
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Manager::getMaxNormF
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and Mems Modeling
// Creation Date : 9/28/2009
//-----------------------------------------------------------------------------
double N_NLS_Manager::getMaxNormF() const
{
  return nlsPtr_->getMaxNormF();
}


//-----------------------------------------------------------------------------
// Function      : N_NLS_Manager::getMaxNormFindex
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
int N_NLS_Manager::getMaxNormFindex () const
{
  return nlsPtr_->getMaxNormFindex ();
}

//-----------------------------------------------------------------------------
// Function      : solve
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Tamara Kolda, SNL, CSMR
// Creation Date : 1/30/02
//-----------------------------------------------------------------------------
int N_NLS_Manager::solve()
{
  int status = 0;

  status = nlsPtr_->solve();
  if (status >= 0)
  {
    if (Teuchos::is_null(exprPtr))
      exprPtr = Teuchos::rcp( new N_UTL_Expression (string("0")) );
    exprPtr->set_accepted_time();
  }

  return status;
}

//-----------------------------------------------------------------------------
// Function      : setAnalysisMode
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Tamara Kolda, SNL, CSMR
// Creation Date : 1/30/02
//-----------------------------------------------------------------------------
void N_NLS_Manager::setAnalysisMode(AnalysisMode mode)
{
  nlsPtr_->setAnalysisMode(mode);
}

//-----------------------------------------------------------------------------
// Function      : resetAll
// Purpose       : like setAnalysisMode, but will also result in NOX
//                 resetting a few extra things.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, 9233
// Creation Date : 1/30/02
//-----------------------------------------------------------------------------
void N_NLS_Manager::resetAll(AnalysisMode mode)
{
  nlsPtr_->resetAll (mode);
}

//-----------------------------------------------------------------------------
// Function      : getNumResidualLoads
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Tamara Kolda, SNL, CSMR
// Creation Date : 1/30/02
//-----------------------------------------------------------------------------
int N_NLS_Manager::getNumResidualLoads()
{
  return nlsPtr_->getNumResidualLoads();
}

//-----------------------------------------------------------------------------
// Function      : getNumJacobianLoads
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Tamara Kolda, SNL, CSMR
// Creation Date : 1/30/02
//-----------------------------------------------------------------------------
int N_NLS_Manager::getNumJacobianLoads()
{
  return nlsPtr_->getNumJacobianLoads();
}

//-----------------------------------------------------------------------------
// Function      : getNumLinearSolves
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Tamara Kolda, SNL, CSMR
// Creation Date : 1/30/02
//-----------------------------------------------------------------------------
int N_NLS_Manager::getNumLinearSolves()
{
  return nlsPtr_->getNumLinearSolves();
}

//-----------------------------------------------------------------------------
// Function      : getNumFailedLinearSolves()
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Tamara Kolda, SNL, CSMR
// Creation Date : 1/30/02
//-----------------------------------------------------------------------------
int N_NLS_Manager::getNumFailedLinearSolves()
{
  return nlsPtr_->getNumFailedLinearSolves();
}

//-----------------------------------------------------------------------------
// Function      : getTotalNumLinearIters
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences Department
// Creation Date : 3/20/02
//-----------------------------------------------------------------------------
unsigned int N_NLS_Manager::getTotalNumLinearIters()
{
  return nlsPtr_->getTotalNumLinearIters();
}

//-----------------------------------------------------------------------------
// Function      : getNumJacobianFactorizations
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Tamara Kolda, SNL, CSMR
// Creation Date : 1/30/02
//-----------------------------------------------------------------------------
int N_NLS_Manager::getNumJacobianFactorizations()
{
  return nlsPtr_->getNumJacobianFactorizations();
}

//-----------------------------------------------------------------------------
// Function      : getTotalLinearSolveTime
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Tamara Kolda, SNL, CSMR
// Creation Date : 1/30/02
//-----------------------------------------------------------------------------
double N_NLS_Manager::getTotalLinearSolveTime()
{
  return nlsPtr_->getTotalLinearSolveTime();
}

//-----------------------------------------------------------------------------
// Function      : getTotalResidualLoadTime
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Tamara Kolda, SNL, CSMR
// Creation Date : 1/30/02
//-----------------------------------------------------------------------------
double N_NLS_Manager::getTotalResidualLoadTime()
{
  return nlsPtr_->getTotalResidualLoadTime();
}

//-----------------------------------------------------------------------------
// Function      : getTotalJacobianLoadTime
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Tamara Kolda, SNL, CSMR
// Creation Date : 1/30/02
//-----------------------------------------------------------------------------
double N_NLS_Manager::getTotalJacobianLoadTime()
{
  return nlsPtr_->getTotalJacobianLoadTime();
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Manager::setSensOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/18/02
//-----------------------------------------------------------------------------
bool N_NLS_Manager::setSensOptions (const N_UTL_OptionBlock & OB)
{
  bool bsuccess = true;

  optionBlockMap_["sens"] = OB;
  return true;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Manager::setSensitivityOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/18/02
//-----------------------------------------------------------------------------
bool N_NLS_Manager::setSensitivityOptions (const N_UTL_OptionBlock & OB)
{
  bool bsuccess = true;

  optionBlockMap_["sensitivity"] = OB;
  return true;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Manager::enableSensitivity
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/30/03
//-----------------------------------------------------------------------------
bool N_NLS_Manager::enableSensitivity ()
{
  bool bsuccess = true;

  if (!setupSensFlag_)
  {
    bool b1 = setupSensitivity_ ();
    bsuccess = bsuccess && b1;
  }
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Manager::calcSensitivity
// Purpose       : This is the controller function for performing a direct
//                 sensitivity calculation.  It is generally called from
//                 the time integration package, as only that package really
//                 knows when to do it... (I may change this later.)
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/31/02
//-----------------------------------------------------------------------------
bool N_NLS_Manager::calcSensitivity ()
{
  bool bsuccess = true;

  if (!setupSensFlag_)
  {
    string msg = "Function N_NLS_Manager::calcSensitivity called,\n";
    msg += "but N_NLS_Manager::enableSensitivity must be called first.\n";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  }

  bsuccess = nlsSensitivityPtr_->solve ();

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Manager::setupSensitivity_
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/31/02
//-----------------------------------------------------------------------------
bool N_NLS_Manager::setupSensitivity_ ()
{
  bool bsuccess = true;
  bool bs1 = true;

  nlsSensitivityPtr_ = new N_NLS_Sensitivity (*nlsPtr_, *topPtr_, commandLine_);
  bs1 = nlsSensitivityPtr_->registerParallelMgr (pdsMgrPtr_); bsuccess = bsuccess && bs1;

  map<string,N_UTL_OptionBlock>::iterator obmIter = optionBlockMap_.find("sens") ;
  if ( obmIter != optionBlockMap_.end() )
  {
    const N_UTL_OptionBlock OB = obmIter->second;
    bs1 = nlsSensitivityPtr_->setOptions (OB);
    bsuccess = bsuccess && bs1;
  }

  obmIter = optionBlockMap_.find("sensitivity") ;
  if ( obmIter != optionBlockMap_.end() )
  {
    const N_UTL_OptionBlock OB = obmIter->second;
    bs1 = nlsSensitivityPtr_->setSensitivityOptions (OB);
    bsuccess = bsuccess && bs1;
  }

  setupSensFlag_ = true;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Manager::registerTopology
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/09/00
//-----------------------------------------------------------------------------
bool N_NLS_Manager::registerTopology (N_TOP_Topology * ptr)
{
  topPtr_ = ptr;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Manager::setReturnCodes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/02/03
//-----------------------------------------------------------------------------
void N_NLS_Manager::setReturnCodes (const N_NLS_ReturnCodes & retCodeTmp)
{
  retCodes_ = retCodeTmp;
  if (initializeAllFlag_)
  {
    nlsPtr_->setReturnCodes (retCodes_);
  }
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Manager::getReturnCodes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical and MEMS Modeling
// Creation Date : 9/28/2009
//-----------------------------------------------------------------------------
N_NLS_ReturnCodes N_NLS_Manager::getReturnCodes() const
{
  // the N_NLS_ReturnCodes structure is very simple, so just
  // return a copy of it.
  return retCodes_;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Manager::setTimeOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/29/05
//-----------------------------------------------------------------------------
bool N_NLS_Manager::setTimeOptions(const N_UTL_OptionBlock & OB)
{
  bool bsuccess = true;
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Manager::setDCOPRestartOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 09/17/07
//-----------------------------------------------------------------------------
bool N_NLS_Manager::setDCOPRestartOptions (const N_UTL_OptionBlock& OB )
{
  optionBlockMap_["dcop_restart"] = OB;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Manager::setICOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 09/17/07
//-----------------------------------------------------------------------------
bool N_NLS_Manager::setICOptions (const N_UTL_OptionBlock& OB )
{
  optionBlockMap_["ic"] = OB;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Manager::setNodeSetOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 09/25/07
//-----------------------------------------------------------------------------
bool N_NLS_Manager::setNodeSetOptions (const N_UTL_OptionBlock& OB )
{
  optionBlockMap_["nodeset"] = OB;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Manager::obtainConductances
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/04/06
//-----------------------------------------------------------------------------
bool N_NLS_Manager::obtainConductances (
        const map<string,double> & inputMap,
        vector<double> & outputVector,
        vector< vector<double> > & jacobian
    )
{
  bool bsuccess = true;

  bsuccess =
    conductanceExtractorPtr_->extract( inputMap, outputVector, jacobian );

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Manager::setMatrixFreeFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, Ting Mei
// Creation Date : 07/29/08
//-----------------------------------------------------------------------------
void N_NLS_Manager::setMatrixFreeFlag ( bool matrixFreeFlag )
{
  matrixFreeFlag_ = matrixFreeFlag;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_Manager::obtainConductances
//
// Purpose       : Same function as above, only this one is for ISO-devices
//                 instead of Vsrc's.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/04/06
//-----------------------------------------------------------------------------
bool N_NLS_Manager::obtainConductances (
        const string & isoName,
        vector< vector<double> > & jacobian )
{
  bool bsuccess = true;

  bsuccess =
    conductanceExtractorPtr_->extract( isoName, jacobian );

  return bsuccess;
}

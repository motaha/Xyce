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
// Filename      : $RCSfile: N_MPDE_Manager.C,v $
// Purpose       :
// Special Notes :
// Creator       : Eric R. Keiter, 9233, Computational Sciences
// Creation Date : 03/12/04
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.197 $
// Revision Date  : $Date: 2014/02/24 23:49:24 $
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


#define MPDE_OLD_TWO_PERIOD 1
// ---------- Standard Includes ----------
#include <iostream>

// ----------   Xyce Includes   ----------

#include <N_MPDE_Manager.h>
#include <N_MPDE_Builder.h>
#include <N_MPDE_Loader.h>

#include <N_MPDE_SawtoothLoader.h>
#include <N_MPDE_Discretization.h>
#include <N_MPDE_DeviceInterface.h>

#include <N_LAS_BlockVector.h>
#include <N_LAS_BlockMatrix.h>
#include <N_LAS_Builder.h>
#include <N_LAS_System.h>

#include <N_LAS_Vector.h>
#include <N_UTL_FFTInterface.hpp>
#include <N_UTL_Misc.h>

#include <N_ANP_AnalysisInterface.h>
#include <N_ANP_AnalysisManager.h>
#include <N_TIA_DataStore.h>
#include <N_TIA_MPDEInterface.h>
#include <N_TIA_TIAParams.h>

#include <N_NLS_Manager.h>

#include <N_TOP_Topology.h>

#include <N_PDS_Manager.h>

#include <N_IO_RestartMgr.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_CmdParse.h>

#include <N_UTL_OptionBlock.h>
#include <N_UTL_Timer.h>

#include <N_PDS_ParMap.h>

#include <N_LOA_Loader.h>
#include <N_LOA_NonlinearEquationLoader.h>

#include <N_ERH_ErrorMgr.h>

#include <N_LAS_PrecondFactory.h>

#include <Epetra_MultiVector.h>

using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::RCP;

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::N_MPDE_Manager
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
N_MPDE_Manager::N_MPDE_Manager(N_IO_CmdParse & cp)
 : debugLevel(0),
   commandLine_ (cp),
   tiaMPDEParams_(cp),
   //tiaParams_(cp),
   blockAnalysisFlag_(false),
   test_(false),
   size_(21),
   tranRunForSize_(false),
   maxCalcSize_(20),
   maxCalcSizeGiven_(false),
   fastSrc_("VIN"),
   fastSrcGiven_(false),
   oscOut_(""),
   oscOutGiven_(false),
   nonLteSteps_(0),
   nonLteStepsGiven_(false),
   period_(1.0),
   periodGiven_(false),
   startUpPeriods_(0),
   startUpPeriodsGiven_(false),
   startUpPeriodsFinished_(false),
   saveIcData_(false),
   fastTimeDisc_(0),
   fastTimeDiscOrder_(1),
   warpMPDE_(false),
   warpMPDEOSCOUT_(-1),
   warpPhase_(0),
   warpPhaseGiven_(false),
   warpPhaseCoeff_(0.0),
   warpPhaseCoeffGiven_(false),
   fftFlag_(false),
   icPer_(10),
   initialCondition_(0),
   mpdeOnFlag_(false),
   mpdeICFlag_(false),
   dcopExitFlag_(false),
   icExitFlag_(false),
   exitSawtoothStep_(-1)
{

  std::string netListFile = "";
  if (commandLine_.getArgumentValue("netlist") != "")
  {
    netListFile = commandLine_.getArgumentValue("netlist");
  }
  else
  {
    std::string msg="N_MPDE_Manager::N_MPDE_Manager ";
    msg += "No netlist specified. \n";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL, msg);
  }

}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::run
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
bool N_MPDE_Manager::run()
{
  bool returnValue=true;

  returnValue = initializeAll_();

  if( test_ )
  {
    runTests_();
  }
  else
  {
    //SHOULD CALL AN INITCOND CLASS HERE
    //mpdeOnFlag_ = false;
    setMPDEOnFlag( false );
    bool ret1 = runInitialCondition_();  returnValue = returnValue && ret1;

    //mpdeOnFlag_ = true;
    setMPDEOnFlag( true );
    bool ret2 = setupMPDEProblem_(); returnValue = returnValue && ret2;

    //SHOULD CALL AN PROBLEM CLASS HERE
    bool ret3 = runMPDEProblem_(); returnValue = returnValue && ret3;
  }

  return returnValue;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::setMPDEAnalysisParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 07/23/08
//-----------------------------------------------------------------------------
bool N_MPDE_Manager::setMPDEAnalysisParams(
  const N_UTL_OptionBlock & tiaParamsBlock)
{
  setBlockAnalysisFlag_( true );

  N_TIA_TIAParams tiaParams = tiaMPDEIfacePtr_->getTIAParams();

  std::list<N_UTL_Param>::const_iterator it_tp;
  std::list<N_UTL_Param>::const_iterator first = tiaParamsBlock.getParams().begin();
  std::list<N_UTL_Param>::const_iterator last = tiaParamsBlock.getParams().end();
  for (it_tp = first; it_tp != last; ++it_tp)
  {
    if (it_tp->uTag()      == "TSTART")
    {
      tiaParams.tStart = it_tp->getImmutableValue<double>();
      tiaParams.tStartGiven = true;
    }
    else if (it_tp->uTag() == "TSTOP")
    {
      tiaParams.finalTime   = it_tp->getImmutableValue<double>();
    }
    else if (it_tp->uTag() == "TSTEP")
    {
      tiaParams.userSpecified_startingTimeStep = it_tp->getImmutableValue<double>();
    }
    else if (it_tp->uTag() == "NOOP" ||
             it_tp->uTag() == "UIC")
    {
      tiaParams.NOOP = 1;
    }
    else if (it_tp->uTag() == "DTMAX")
    {
      tiaParams.maxTimeStep = it_tp->getImmutableValue<double>();
      tiaParams.maxTimeStepGiven = true;
#ifdef Xyce_DEBUG_ANALYSIS
      if (tiaParams.debugLevel > 0)
      {
        Xyce::dout() << "setting maxTimeStep = " << tiaParams.maxTimeStep << std::endl;
      }
#endif
    }

  }

  if (tiaParams.finalTime <= tiaParams.tStart || tiaParams.finalTime <= 0 || tiaParams.tStart < 0)
  {
    std::ostringstream ost;
    ost << " In N_MPDE_Manager::setMPDEAnalysisParams: " << std::endl;
    ost << "Final time of " << tiaParams.finalTime
        << " is earlier or same as start time of "
        << tiaParams.tStart << std::endl;
    ost << " Check netlist for invalid .MPDE specification " << std::endl;
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, ost.str());
  }

  // if starting time steps is baloney, set then to a default value.
  // ERK.  This trap is redudant with an identical one in step error
  // control.
  if ( tiaParams.userSpecified_startingTimeStep <= 0.0 )
  {
    tiaParams.userSpecified_startingTimeStep = 1.0e-10;
  }

#ifdef Xyce_DEBUG_ANALYSIS
  if (tiaParams.debugLevel > 0)
  {
    Xyce::dout() << std::endl
                 << Xyce::section_divider << std::endl
                 << "  Transient simulation parameters" << std::endl
                 << "  initial time = " << tiaParams.initialTime << std::endl
                 << "  final time   = " << tiaParams.finalTime << std::endl
                 << "  starting time step = " << tiaParams.userSpecified_startingTimeStep << std::endl
                 << "  tStart (time of first output << std::endl = " << tiaParams.tStart << std::endl
                 << (!tiaParams.NOOP ? "  NOOP/UIC is NOT set" : "  NOOP/UIC is set") << std::endl
                 << Xyce::section_divider << std::endl;
  }
#endif

  anaIntPtr_->registerTIAParams(tiaParams);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::setMPDEOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
bool N_MPDE_Manager::setMPDEOptions(const N_UTL_OptionBlock & OB)
{
  std::string msg;

  //MOST OF THIS SHOULD BE PUSHED TO INITCOND AND PROBLEM CLASSES
  std::list<N_UTL_Param>::const_iterator iterPL = OB.getParams().begin();
  std::list<N_UTL_Param>::const_iterator  endPL = OB.getParams().end();

  for( ; iterPL != endPL; ++iterPL )
  {
    ExtendedString tag = iterPL->tag();
    tag.toUpper();

    if ( tag == "N2" )
    {
      size_    = iterPL->getImmutableValue<int>();
    }
    else if ( tag == "AUTON2" )
    {
      tranRunForSize_ = iterPL->getImmutableValue<bool>();
    }
    else if ( tag == "AUTON2MAX" )
    {
      maxCalcSize_ = iterPL->getImmutableValue<int>();
      maxCalcSizeGiven_ = true;
    }
    else if (std::string(tag,0,6) == "OSCSRC" )
    {
      fastSrc_  = iterPL->stringValue();
      srcVec_.push_back(fastSrc_);
      fastSrcGiven_ = true;
    }
    else if ( tag == "NONLTESTEPS" )
    {
      nonLteSteps_ = iterPL->getImmutableValue<int>();
      nonLteStepsGiven_ = true;
    }
    else if ( tag == "STARTUPPERIODS" )
    {
      startUpPeriods_ = iterPL->getImmutableValue<int>();
      startUpPeriodsGiven_ = true;
    }
    else if( tag == "SAVEICDATA" )
    {
      saveIcData_ = true;
    }
    else if( tag == "OSCOUT" )
    {
      oscOut_   = iterPL->stringValue();
      oscOutGiven_ = true;
    }
    else if( tag == "PHASE" )
    {
      warpPhase_ = iterPL->getImmutableValue<int>();
      warpPhaseGiven_ = true;
    }
    else if( tag == "PHASECOEFF" )
    {
      warpPhaseCoeff_ = iterPL->getImmutableValue<double>();
      warpPhaseCoeffGiven_ = true;
    }
    else if( tag == "TEST" )
    {
      test_     = static_cast<bool> (iterPL->getImmutableValue<int>());
    }
    else if( tag == "T2" )
    {
      period_   = iterPL->getImmutableValue<double>();
      periodGiven_ = true;
    }
    else if( tag == "WAMPDE" )
    {
      warpMPDE_ = static_cast<bool>( iterPL->getImmutableValue<int>() );
    }
    else if( tag == "FREQDOMAIN" )
    {
      fftFlag_  = static_cast<bool>( iterPL->getImmutableValue<int>() );
    }
    else if( tag == "ICPER" )
    {
      icPer_    = iterPL->getImmutableValue<int>();
    }
    else if( tag == "IC" )
    {
      initialCondition_ = iterPL->getImmutableValue<int>();
    }
    else if( tag == "DIFF" )
    {
      fastTimeDisc_ = iterPL->getImmutableValue<int>();
    }
    else if( tag == "DIFFORDER" )
    {
      fastTimeDiscOrder_ = iterPL->getImmutableValue<int>();
    }
    else if (tag == "DCOPEXIT" )
    {
      dcopExitFlag_ = static_cast<bool>( iterPL->getImmutableValue<int>() );
    }
    else if (tag == "ICEXIT" )
    {
      icExitFlag_ = static_cast<bool>( iterPL->getImmutableValue<int>() );
    }
    else if (tag == "EXITSAWTOOTHSTEP" )
    {
      exitSawtoothStep_ = iterPL->getImmutableValue<int>();
    }
    else if (tag == "DEBUGLEVEL" )
    {
      debugLevel = iterPL->getImmutableValue<int>();
    }
    else
    {
      msg = "N_MPDE_Manager::setMPDEOptions";
      msg += " Unrecognized MPDEINT option: "+tag+"\n";
      N_ERH_ErrorMgr::report (N_ERH_ErrorMgr::USR_WARNING, msg);
    }
  }

  if( (fastTimeDisc_ == 1) && (fastTimeDiscOrder_ < 2 ) )
  {
    // if user asked for centered differences, then default to
    // second order if none was specified otherwise the Discretization
    // class will return backwards differences which may not be what
    // the user expects
    fastTimeDiscOrder_ = 2;
  }

  printParams_ ();

  return true;
}


//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::registerTranMPDEOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, 9233, Computational Sciences
// Creation Date : 08/20/07
//-----------------------------------------------------------------------------
bool N_MPDE_Manager::registerTranMPDEOptions(const N_UTL_OptionBlock & OB)
{
  saved_tranMPDEOB_ = OB;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::printParams_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, 9233, Computational Sciences
// Creation Date : 03/30/04
//-----------------------------------------------------------------------------
void N_MPDE_Manager::printParams_ ()
{
  Xyce::dout() << "\n" << std::endl;
  Xyce::dout() << Xyce::section_divider << std::endl;
  Xyce::dout() <<
                         "\n***** MPDE options:\n" << std::endl;

  Xyce::dout() <<
       "\tN2:\t\t\t" <<  size_ << std::endl;

  std::string msg;
  if (fastSrcGiven_)
  {
    for (unsigned int i1=0;i1<srcVec_.size();++i1)
    {
      msg = "\toscsrc:\t\t\t"+srcVec_[i1];
      Xyce::dout() <<msg << std::endl;
    }
  }
  else
  {
    msg = "\toscsrc:\t\t\tNot Specified";
    Xyce::dout() <<msg << std::endl;
  }

  if (oscOutGiven_)
  {
    msg = "\toscout:\t\t\t"+oscOut_;
  }
  else
  {
    msg = "\toscout:\t\t\tNot Specified";
  }
  Xyce::dout() <<msg << std::endl;

  Xyce::dout() <<
       "\tT2:\t\t\t" <<  period_ << std::endl;

  if( fastTimeDisc_ == 0 )
    msg = "\tdiscretization:\tBackward";
  else if( fastTimeDisc_ == 1 )
    msg = "\tdiscretization:\tCentered";
  else
    msg = "\tdiscretization:\tUNKNOWN";
  Xyce::dout() <<msg << std::endl;

  Xyce::dout() <<
       "\tDisc Order:\t\t" <<  fastTimeDiscOrder_ << std::endl;

  if (test_)
  {
    Xyce::dout() <<
      "\tTest Mode:" << std::endl;
  }
  else
  {
    Xyce::dout() <<
      "\tFull MPDE Mode(not test mode)" << std::endl;
  }

  if (warpMPDE_)
  {
    Xyce::dout() <<
      "\tWarpedMPDE:\t\tON" << std::endl;
  }
  else
  {
    Xyce::dout() <<
      "\tWarpedMPDE:\t\tOFF" << std::endl;
  }

  if (fftFlag_)
  {
    Xyce::dout() <<
      "\tFrequency domain:\tON" << std::endl;
  }
  else
  {
    Xyce::dout() <<
      "\tFrequency domain:\tOFF" << std::endl;
  }

  Xyce::dout() <<
       "\tICPER:\t\t\t" <<  icPer_ << std::endl;

  if (initialCondition_ == MPDE_IC_DCOP)
  {
    Xyce::dout() <<
       "\tInitial Condition:\tDCOP" << std::endl;
  }
  else if (initialCondition_ == MPDE_IC_SAWTOOTH)
  {
    Xyce::dout() <<
       "\tInitial Condition:\tSAWTOOTH" << std::endl;
  }

  Xyce::dout() << "\n" << std::endl;
  Xyce::dout() << Xyce::section_divider << std::endl;
  Xyce::dout() << "\n" << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::initializeAll_
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, 9233, Computational Sciences
// Creation Date : 03/12/04
//-----------------------------------------------------------------------------
bool N_MPDE_Manager::initializeAll_()
{
  bool returnValue = true;

#ifdef Xyce_DEBUG_MPDE
  Xyce::dout() << std::endl;
  Xyce::dout() << Xyce::section_divider << std::endl;
  Xyce::dout() <<
                          "  N_MPDE_Manager::initializeAll_" << std::endl;
#endif // Xyce_DEBUG_MPDE
  //Need this for IC to get "MPDE size" vectors
  //-----------------------------------------
#ifdef Xyce_DEBUG_MPDE
  Xyce::dout() << "N_MPDE_Manager::initializeAll_[Construct Builder]\n";
#endif // Xyce_DEBUG_MPDE

  // Set up the gid index for the osc out variable:
  initializeOscOut_();

  // may need to dynamically find size_ here

  //-----------------------------------------

  //Store copy of transient TIAParams for MPDE run
  tiaMPDEParams_ = anaIntPtr_->getTIAParams() ;

  //Register MPDE State with Dev Pkg
  //devInterfacePtr_->registerMPDEState( &state_ );

  //Set MPDE src
  //devInt_->setMPDESrc( fastSrc_ );
    
  std::vector<double> srcPeriods;
  //srcPeriods = devInterfacePtr_->registerFastSources( srcVec_ );
  srcPeriods = mpdeDeviceInterfacePtr_->getFastSourcePeriod( srcVec_ );

  // If running straight MPDE, the only option is to get the period
  // from the "oscsrc", or fast source.  If running warped, it is still
  // possible to get it from a fast source, but the T2 (period_) parameter
  // takes precedence.
  std::string msg("");
  if (periodGiven_ && fastSrcGiven_)
  {
    msg = "Note:  oscsrc is being ignored, as the user has set T2.\n";
    Xyce::dout() << msg << std::endl;
  }

  if (fastSrcGiven_  && !periodGiven_)
  {
    // don't assume the periods are all the same:
    period_ = srcPeriods[0];
    int numSrc = srcPeriods.size();
    for( int i=1; i<numSrc; i++ )
    {
      if( srcPeriods[i] > period_ )
      {
        period_ = srcPeriods[i];
      }
    }
  }

  if (warpMPDE_ && !(periodGiven_) && !(fastSrcGiven_))
  {
    msg = "The Warped MPDE needs for the user to set oscsrc or T2 (preferably T2).\n";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0,msg);
  }

  if (!warpMPDE_ && !(fastSrcGiven_))
  {
    msg = "The MPDE algorithm needs the user to set oscsrc.\n";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0,msg);
  }

  //
  // If it was requested, advance the solution a fixed number of startup periods
  //
  if( startUpPeriodsGiven_ )
  {
    bool startupPeriodsSuccess = runStartupPeriods_();
    if (!startupPeriodsSuccess)
    {
      std::string msg = "N_MPDE_Manager::initializeAll_().  initializeAll failed to calculate the startup periods.\n";
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL_0,msg);
    }
    returnValue = returnValue && startupPeriodsSuccess;
  }

  if( !tranRunForSize_ )
  {
    //Simple discrete time points for fast time scale
    //-----------------------------------------
    fastTimes_.resize(size_+1);
    double fastStep = period_/static_cast<double>(size_);
    int i=0;
    for( i = 0; i < size_; ++i )
    {
      fastTimes_[i] = static_cast<double>(i) * fastStep;
    }
    fastTimes_[size_] = period_;
    //-----------------------------------------
#ifdef Xyce_DEBUG_MPDE
    if( debugLevel > 0 )
    {
      Xyce::dout() << "MPDE Fast Step: " << fastStep << std::endl;
      for( i = 0; i < size_; ++i )
        Xyce::dout() << "fastTimes_["<<i<<"] = " << fastTimes_[i] << std::endl;
      Xyce::dout() << Xyce::section_divider << std::endl;
    }
#endif // Xyce_DEBUG_MPDE
  }
  else
  {
    returnValue = runTransientIC_();
    if (!returnValue)
    {
      std::string msg = "N_MPDE_Manager::initializeAll_().  initializeAll failed to compute the transient initial condition.\n";
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL_0,msg);
    }
    filterFastTimePoints_();
  }

  // now that we have size_, continue with the initialization

  mpdeDiscPtr_ = rcp(new N_MPDE_Discretization(
     static_cast<N_MPDE_Discretization::Type>(fastTimeDisc_), fastTimeDiscOrder_ ));

  mpdeBuilderPtr_ = rcp(new N_MPDE_Builder( rcp(this,false), size_, mpdeDiscPtr_, warpMPDE_ ));

#ifdef Xyce_DEBUG_MPDE
  Xyce::dout() << "N_MPDE_Manager::initializeAll_[Generate Maps]\n";
#endif // Xyce_DEBUG_MPDE
  mpdeBuilderPtr_->generateMaps( rcp( pdsMgrPtr_->getParallelMap( "SOLUTION" ), false) );
  mpdeBuilderPtr_->generateStateMaps( rcp( pdsMgrPtr_->getParallelMap( "STATE" ), false) );
  mpdeBuilderPtr_->generateStoreMaps( rcp( pdsMgrPtr_->getParallelMap( "STORE" ), false) );

  // Setup WaMPDE Phase Condition Object
  if (warpMPDE_)
  {
    // These three values are initialized during generateMaps in the MPDE Builder.
    int offset = mpdeBuilderPtr_->getMPDEOffset();
    int omegaGID = mpdeBuilderPtr_->getMPDEomegaGID();
    int phiGID = mpdeBuilderPtr_->getMPDEphiGID();
    int augProcID = mpdeBuilderPtr_->getMPDEaugProcID();

#ifdef Xyce_DEBUG_MPDE
    std::cout << "N_MPDE_Manager::initializeAll_[ offset = " << offset << ", omegaGID = " << omegaGID 
               << ", phiGID = " << phiGID << ", augProcID = " << augProcID << " ] " << std::endl;
#endif

    warpMPDEPhasePtr_ = rcp(
        new N_MPDE_WarpedPhaseCondition(warpPhase_, warpPhaseCoeff_, warpMPDEOSCOUT_, omegaGID, offset, size_, phiGID, augProcID)
        );
    mpdeBuilderPtr_->setWarpedPhaseCondition(warpMPDEPhasePtr_);
  }

  // Now that the builder is set up, allocate petra vectors, etc.
  mpdeICVectorPtr_ = rcp_dynamic_cast<N_LAS_BlockVector>(rcp(mpdeBuilderPtr_->createVector()));
  mpdeICStateVectorPtr_ = rcp_dynamic_cast<N_LAS_BlockVector>(rcp(mpdeBuilderPtr_->createStateVector()));
  mpdeICQVectorPtr_ = rcp_dynamic_cast<N_LAS_BlockVector>(rcp(mpdeBuilderPtr_->createVector()));
  mpdeICStoreVectorPtr_ = rcp_dynamic_cast<N_LAS_BlockVector>(rcp(mpdeBuilderPtr_->createStoreVector()));

#ifdef Xyce_DEBUG_MPDE
  Xyce::dout() << "N_MPDE_Manager::initializeAll_[Finished]\n";
#endif // Xyce_DEBUG_MPDE

  // set the fast source flag in the devices
  srcPeriods = mpdeDeviceInterfacePtr_->registerFastSources( srcVec_ );

  return returnValue;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::initializeOscOut_
// Purpose       : Sets up the warpMPDEOSCOUT_  global id (GID) index.
//
// Special Notes : This is called from initializeAll_.
//
//   This has to be done here  (rather than setMPDEOptions) because
//   it can only be done after the non-MPDE topology has been set up.
//
// Scope         : private
// Creator       : Eric Keiter, 1437, Computational Sciences
// Creation Date : 12/16/07
//-----------------------------------------------------------------------------
bool N_MPDE_Manager::initializeOscOut_()
{
  if (oscOutGiven_)
  {
    bool oscOutError(false);
    ExtendedString varName(oscOut_);
    varName.toUpper();
    //std::string varName(oscOut_);
    std::list<int> oscoutList, dummyList;
    char type1;
    int pos1 = varName.find("I(");
    int pos2 = varName.find("V(");
    int pos3 = varName.find(")");
    if (pos1!=(int)std::string::npos && pos3!=(int)std::string::npos) // this is a current variable
    {
      pos1+=2;
      int len = pos3-pos1;
      std::string tmpVar (varName,pos1,len);
#ifdef Xyce_DEBUG_MPDE
      Xyce::dout() << "tmpVar (for I-oscout) = " << tmpVar << std::endl;
#endif
      topoMgrPtr_->getNodeSVarGIDs(NodeID(tmpVar, Xyce::_DNODE), oscoutList, dummyList, type1);
    }
    else if (pos2 !=(int)std::string::npos && pos3!=(int)std::string::npos) // this is a voltage variable
    {
      pos2+=2;
      int len = pos3-pos2;
      std::string tmpVar (varName,pos2,len);
#ifdef Xyce_DEBUG_MPDE
      Xyce::dout() << "tmpVar (for V-oscout) = " << tmpVar << std::endl;
#endif
      topoMgrPtr_->getNodeSVarGIDs(NodeID(tmpVar, Xyce::_VNODE), oscoutList, dummyList, type1);
    }
    else
    {
      // do nothing, assume that the varName from the netlist requires no
      // modification.
      topoMgrPtr_->getNodeSVarGIDs(NodeID(varName,-1), oscoutList, dummyList, type1);
    }

    if (oscoutList.size()==1)
    {
      warpMPDEOSCOUT_=oscoutList.front();
#ifdef Xyce_DEBUG_MPDE
      if (warpMPDEOSCOUT_>=0)
        std::cout << "warpMPDEOSCOUT = " << warpMPDEOSCOUT_ << std::endl;
#endif // Xyce_DEBUG_MPDE
    }
#ifndef Xyce_PARALLEL_MPI
    else  // This is not an error in parallel.
    {
      oscOutError = true;
    }
#else
    int tmpOscOut = warpMPDEOSCOUT_;
    pdsMgrPtr_->getPDSComm()->maxAll( &tmpOscOut, &warpMPDEOSCOUT_, 1 );
    if ( warpMPDEOSCOUT_ < 0 ) 
    {
      oscOutError = true;
    }
#endif

    if (oscOutError)
    {
      std::string msg = "N_MPDE_Manager::setMPDEOptions";
      msg += " Unrecognized value for MPDE option:  oscOut_=";
      msg += oscOut_ + "\n";
      N_ERH_ErrorMgr::report (N_ERH_ErrorMgr::USR_WARNING, msg);
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::runInitialCondition_
// Purpose       : Execution of initial condition phase
// Special Notes : Produces an "MPDE" size initial value
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/19/04
//-----------------------------------------------------------------------------
bool N_MPDE_Manager::runInitialCondition_()
{
  bool returnValue=true;

  Xyce::lout() << " ***** Running MPDE initial conditions....\n" << std::endl;
#ifdef Xyce_DEBUG_MPDE
  Xyce::dout() << std::endl
               << Xyce::section_divider << std::endl
               << "  N_MPDE_Manager::runInitialCondition" << std::endl;
#endif // Xyce_DEBUG_MPDE
  int n2 = size_;
  std::string msg;

  RCP<N_TIA_DataStore> dsPtr = anaIntPtr_->getAnalysisMgr()->getTIADataStore();

  if (!warpMPDE_)
  {
    if( is_null(dcOpSolVecPtr_) )
    {
      // haven't done a dc op yet, so do it now
      runDCOP_();
    }

    // place dc op values in mpde block vectors
    mpdeICVectorPtr_->block(0) = *dcOpSolVecPtr_;
    mpdeICStateVectorPtr_->block(0) = *dcOpStateVecPtr_;
    mpdeICQVectorPtr_->block(0) = *dcOpQVecPtr_;
    mpdeICStoreVectorPtr_->block(0) = *dcOpStoreVecPtr_;
  }
  else
  {
    /* This initialization is wrong for WaMPDE.  
       These vectors aren't available yet.
    mpdeICVectorPtr_->block(0) = *(dsPtr->fastTimeSolutionVec[0]);
    mpdeICStateVectorPtr_->block(0) = *(dsPtr->fastTimeStateVec[0]);
    mpdeICQVectorPtr_->block(0) = *(dsPtr->fastTimeQVec[0]);
    mpdeICStoreVectorPtr_->block(0) = *(dsPtr->fastTimeStoreVec[0]);
       Initialize block 0 from the solution 
    */

    //std::cout << "Setting initial condition vectors for warped MPDE" << std::endl;
    
    mpdeICVectorPtr_->block(0) = *(dsPtr->currSolutionPtr);
    mpdeICStateVectorPtr_->block(0) = *(dsPtr->currStatePtr);
    mpdeICQVectorPtr_->block(0) = *(dsPtr->daeQVectorPtr);
    mpdeICStoreVectorPtr_->block(0) = *(dsPtr->currStorePtr);

    //std::cout << "Set initial condition vectors for warped MPDE" << std::endl;
  }
    
#ifdef Xyce_DEBUG_MPDE
  if (dcopExitFlag_) exit(0);
#endif // Xyce_DEBUG_MPDE

  switch (initialCondition_)
  {
  case MPDE_IC_DCOP:
    // 03/23/04 tscoffe:  Here's the meta algorithm for the
    // constant DCOP IC strategy:
    for (int i=1 ; i<n2 ; ++i)
    {
      //std::cout << "Computing the initial condition for IC block " << i << std::endl;
      mpdeState_.fastTime = fastTimes_[i];
      tiaMPDEIfacePtr_->runDCOP();
      mpdeICVectorPtr_->block(i) = *(tiaMPDEIfacePtr_->getFinalSolution());
      mpdeICStateVectorPtr_->block(i) = *(tiaMPDEIfacePtr_->getStateFinalSolution());
      mpdeICQVectorPtr_->block(i) = *(tiaMPDEIfacePtr_->getQVectorFinalSolution());
      mpdeICStoreVectorPtr_->block(i) = *(tiaMPDEIfacePtr_->getStoreFinalSolution());
    }
    
    //std::cout << "Done computing the initial condition" << std::endl;
    
    break;

  case MPDE_IC_SAWTOOTH:
    {// <----This bracket is neccessary for OSX compile. Don't remove.ERK

      // 03/23/04 tscoffe:  Here's the meta algorithm for the
      // "sawtooth" IC strategy:

      //Setting up specialized Loader
      N_MPDE_SawtoothLoader stLoader( mpdeState_ );
      stLoader.registerTIA( anaIntPtr_ );
      stLoader.registerAppLoader( appLoaderPtr_ );
      stLoader.setTimeShift( fastTimes_[0] );
      anaIntPtr_->registerLoader(&stLoader);

      // get the time integrator params, and save a copy.
      //tiaParams_ = tiaMPDEIfacePtr_->getTIAParams();
      N_TIA_TIAParams & tiaParams_ = tiaMPDEIfacePtr_->getTIAParams();
      N_TIA_TIAParams tiaParamsSave = tiaParams_;

      for (int i=1 ; i<n2 ; ++i)
      {
        // set DAE initial time = 0.0
        tiaParams_.initialTime = 0.0;
        // set DAE final time = h2
        double h2 = fastTimes_[i]-fastTimes_[i-1];
        tiaParams_.finalTime = h2;

        //tiaParams_.finalTime = fastTimes_[i];
        tiaParams_.pauseTime = tiaParams_.finalTime;
        // We need to do an initialization every time...
        tiaMPDEIfacePtr_->reInitialize();

#ifdef Xyce_DEBUG_MPDE
        Xyce::dout() << "Beginning SAWTOOTH IC: i = " << i
          << " initialTime = " << tiaParams_.initialTime
          << " finalTime = " << tiaParams_.finalTime << std::endl;
#endif // Xyce_DEBUG_MPDE

        // register new tiaParams with time integrator
        tiaMPDEIfacePtr_->registerTIAParams(tiaParams_); // redundant?

        // Set t = t + (i-1.0)*h2 in MPDE source in bhat
        stLoader.setTimeShift( fastTimes_[i-1] );

        // set DAE initial condition to mpdeICVectorPtr_[i-1]
        tiaMPDEIfacePtr_->setInitialCondition(&( mpdeICVectorPtr_->block(i-1)));
        tiaMPDEIfacePtr_->setStateInitialCondition(&( mpdeICStateVectorPtr_->block(i-1)));
        tiaMPDEIfacePtr_->setQVectorInitialCondition(&( mpdeICQVectorPtr_->block(i-1)));
        tiaMPDEIfacePtr_->setStoreInitialCondition(&( mpdeICStoreVectorPtr_->block(i-1)));
        // solve the DAE to time h2
        tiaMPDEIfacePtr_->runTransient();
        // put the dsPtr->currentSolutionPtr into mpdeICVectorPtr_[i]
        mpdeICVectorPtr_->block(i) = *(tiaMPDEIfacePtr_->getFinalSolution());
        mpdeICStateVectorPtr_->block(i) = *(tiaMPDEIfacePtr_->getStateFinalSolution());
        mpdeICStoreVectorPtr_->block(i) = *(tiaMPDEIfacePtr_->getStoreFinalSolution());
#ifdef Xyce_DEBUG_MPDE
        Xyce::dout() << "End SAWTOOTH IC: i = " << i <<std::endl;
        if (i==exitSawtoothStep_) exit(0);
#endif // Xyce_DEBUG_MPDE
      }
      anaIntPtr_->registerLoader(&*appLoaderPtr_);

      // restore the saved copy of time integration parameters.
      tiaMPDEIfacePtr_->registerTIAParams(tiaParamsSave);
    }
    break;

  case MPDE_IC_TRANSIENT:
    {
#ifdef Xyce_DEBUG_MPDE
      if( debugLevel > 0 )
        {
          for( unsigned int i=0; i<indicesUsed_.size(); i++ )
            Xyce::dout() << "indicesUsed_[ " << i << " ] = " << indicesUsed_[i] << std::endl;
        }
#endif

      // we did an initial transient run so pick up ic data from that
//      for (int i=1 ; i<n2 ; ++i)
      for (int i=0 ; i<n2 ; ++i)
      {
        RCP<N_TIA_DataStore> dsPtr = anaIntPtr_->getAnalysisMgr()->getTIADataStore();
#ifdef Xyce_DEBUG_MPDE
        if( debugLevel > 0 )
        {
          Xyce::dout() << "Loading initial condition data from time: fastTimes_["
            << i << "] = " << fastTimes_[i] << std::endl;
        }
        if ( debugLevel > 1 )
        {
          Xyce::dout() << "mpdeICVectorPtr_->block(" << i << ") = dsPtr->fastTimeSolutionVec[" << indicesUsed_[i] << "] = " << std::endl;
          (dsPtr->fastTimeSolutionVec[indicesUsed_[i]])->printPetraObject(Xyce::dout());
          Xyce::dout() << "mpdeICStateVectorPtr_->block(" << i << ") = dsPtr->fastTimeStateVec[" << indicesUsed_[i] << "] = " << std::endl;
          (dsPtr->fastTimeStateVec[indicesUsed_[i]])->printPetraObject(Xyce::dout());
          Xyce::dout() << "mpdeICQVectorPtr_->block(" << i << ") = dsPtr->fastTimeQVec[" << indicesUsed_[i] << "] = " << std::endl;
          (dsPtr->fastTimeQVec[indicesUsed_[i]])->printPetraObject(Xyce::dout());
          Xyce::dout() << "mpdeICStoreVectorPtr_->block(" << i << ") = dsPtr->fastTimeStoreVec[" << indicesUsed_[i] << "] = " << std::endl;
          (dsPtr->fastTimeStoreVec[indicesUsed_[i]])->printPetraObject(Xyce::dout());
        }
#endif // Xyce_DEBUG_MPDE

        mpdeICVectorPtr_->block(i) = *(dsPtr->fastTimeSolutionVec[indicesUsed_[i]]);
        mpdeICStateVectorPtr_->block(i) = *(dsPtr->fastTimeStateVec[indicesUsed_[i]]);
        mpdeICQVectorPtr_->block(i) = *(dsPtr->fastTimeQVec[indicesUsed_[i]]);
        mpdeICStoreVectorPtr_->block(i) = *(dsPtr->fastTimeStoreVec[indicesUsed_[i]]);
      }

    }
    break;

  case MPDE_IC_TRAN_MAP:
    {
      // at this point fastTimes_ is either setup
      // with fixed spaced points or a variable mesh.

      if( !tranRunForSize_ )
      {
        returnValue = runTransientIC_();
        // normally after runTransientIC_() we would filter the resulting
        // points down.  We only get this this part of the code
        // if one called this initial condition type, but wanted fixed
        // spaced points on the fast time scale.  So, we don't call
        // filterFastTimePoints_() here.
      }

      // use the last point for the endICvector
      RCP<N_TIA_DataStore> dsPtr = anaIntPtr_->getAnalysisMgr()->getTIADataStore();
      endIcSolVecPtr_ = rcp((dsPtr->fastTimeSolutionVec[indicesUsed_[n2-1]]));
      endIcStateVecPtr_ = rcp((dsPtr->fastTimeStateVec[indicesUsed_[n2-1]]));
      endIcQVecPtr_ = rcp((dsPtr->fastTimeQVec[indicesUsed_[n2-1]]));
      endIcStoreVecPtr_ = rcp((dsPtr->fastTimeStoreVec[indicesUsed_[n2-1]]));

      // now we linearlly interpolate between the data in
      // dcOpSolVecPtr_    and  endIcSolVecPtr_
      // dcOpStateVecPtr_  and  endIcStateVecPtr_
      // dcOpQVecPtr_      and  endIcQVecPtr_
      //
      // and, if enabled:
      // dcOpStoreVecPtr_  and  endIcStoreVecPtr_
      //
      // for each point in fastTimes_ and use that value
      // as an initial condition to go from t=fastTimes_[i] to period_
      //

      mpdeICVectorPtr_->block(0) = *endIcSolVecPtr_,
      mpdeICStateVectorPtr_->block(0) = *endIcStateVecPtr_;
      mpdeICQVectorPtr_->block(0) = *endIcQVecPtr_;
      mpdeICStoreVectorPtr_->block(0) = *endIcStoreVecPtr_;

      //Setting up specialized Loader
      N_MPDE_SawtoothLoader stLoader( mpdeState_ );
      stLoader.registerTIA( anaIntPtr_ );
      stLoader.registerAppLoader( appLoaderPtr_ );
      stLoader.setTimeShift( fastTimes_[0] );
      anaIntPtr_->registerLoader(&stLoader);

      // get the time integrator params, and save a copy.
      //tiaParams_ = tiaMPDEIfacePtr_->getTIAParams();
      N_TIA_TIAParams & tiaParams_ = tiaMPDEIfacePtr_->getTIAParams();
      N_TIA_TIAParams tiaParamsSave = tiaParams_;

      RCP<N_LAS_Vector> interpIcSolVecPtr = rcp(new N_LAS_Vector( *endIcSolVecPtr_ ));
      RCP<N_LAS_Vector> interpIcStateVecPtr = rcp(new N_LAS_Vector( *endIcStateVecPtr_ ));
      RCP<N_LAS_Vector> interpIcQVecPtr = rcp(new N_LAS_Vector( *endIcQVecPtr_ ));
      RCP<N_LAS_Vector> interpIcStoreVecPtr = rcp(new N_LAS_Vector( *endIcStoreVecPtr_ ));

      for (int i=1 ; i<n2 ; ++i)
      {
        // set DAE initial time = 0.0
        tiaParams_.initialTime = tiaMPDEParams_.initialTime + fastTimes_[n2-i];
        tiaParams_.finalTime = tiaMPDEParams_.initialTime + period_;
        tiaParams_.restartTimeStepScale = 1.0;  // try to take big steps.

        //tiaParams_.finalTime = fastTimes_[i];
        tiaParams_.pauseTime = tiaParams_.finalTime;
        // We need to do an initialization every time...
        tiaMPDEIfacePtr_->reInitialize();

        // here's were we interpolate to find the initial condition
        double fraction = fastTimes_[n2-i] / period_;

#ifdef Xyce_DEBUG_MPDE
        Xyce::dout() << "Beginning MPDE_IC_TRAN_MAP IC: i = " << i
          << " fraction = " << fraction
          << " initialTime = " << tiaParams_.initialTime
          << " finalTime = " << tiaParams_.finalTime << std::endl;
#endif // Xyce_DEBUG_MPDE

        // register new tiaParams with time integrator
        tiaMPDEIfacePtr_->registerTIAParams(tiaParams_); // redundant?

        // Set t = t + (i-1.0)*h2 in MPDE source in bhat
        stLoader.setTimeShift( -fastTimes_[n2-i] );

        interpIcSolVecPtr->putScalar( 0.0 );
        interpIcStateVecPtr->putScalar( 0.0 );
        interpIcQVecPtr->putScalar( 0.0 );
        interpIcStoreVecPtr->putScalar( 0.0 );
        interpIcSolVecPtr->linearCombo( (1.0-fraction), *dcOpSolVecPtr_, fraction, *endIcSolVecPtr_ );
        interpIcStateVecPtr->linearCombo( (1.0-fraction), *dcOpStateVecPtr_, fraction, *endIcStateVecPtr_ );
        interpIcQVecPtr->linearCombo( (1.0-fraction), *dcOpQVecPtr_, fraction, *endIcQVecPtr_ );
        interpIcStoreVecPtr->linearCombo( (1.0-fraction), *dcOpStoreVecPtr_, fraction, *endIcStoreVecPtr_ );

        // set DAE initial condition to mpdeICVectorPtr_[i-1]
        tiaMPDEIfacePtr_->setInitialCondition(interpIcSolVecPtr.get());
        tiaMPDEIfacePtr_->setStateInitialCondition(interpIcStateVecPtr.get());
        tiaMPDEIfacePtr_->setQVectorInitialCondition(interpIcQVecPtr.get());
        tiaMPDEIfacePtr_->setStoreInitialCondition(interpIcStoreVecPtr.get());

        // solve the DAE to time period_
        tiaMPDEIfacePtr_->runTransient();

        // put the dsPtr->currentSolutionPtr into mpdeICVectorPtr_[i]
        mpdeICVectorPtr_->block(i) = *(tiaMPDEIfacePtr_->getFinalSolution());
        mpdeICStateVectorPtr_->block(i) = *(tiaMPDEIfacePtr_->getStateFinalSolution());
        mpdeICQVectorPtr_->block(i) = *(tiaMPDEIfacePtr_->getQVectorFinalSolution());
        mpdeICStoreVectorPtr_->block(i) = *(tiaMPDEIfacePtr_->getStoreFinalSolution());
#ifdef Xyce_DEBUG_MPDE
        Xyce::dout() << "End MPDE_IC_TRAN_MAP IC: i = " << i <<std::endl;
#endif // Xyce_DEBUG_MPDE
      }

      // release these vectors from reference counting
      endIcSolVecPtr_.release();
      endIcStateVecPtr_.release();
      endIcQVecPtr_.release();
      endIcStoreVecPtr_.release();

      anaIntPtr_->registerLoader(&*appLoaderPtr_);

      // restore the saved copy of time integration parameters.
      tiaMPDEIfacePtr_->registerTIAParams(tiaParamsSave);
    }
    break;

  case MPDE_IC_TWO_PERIOD:
    {
 
      // The runTransientIC_ function would have done two periods of
      // the fast time source.  Now we need to look at the points in
      // the second period and match them up with the ones in the first
      // period and then do the inerpolation.

      // fastTimes_[] will already have been offset to run between 0
      // and one period_.  We can use indicesUsed_[] and the timeStep[]
      // array still in the dataStore to find solutions in the first
      // cycle that lie close to those in the fastTimes_[] cycle.

      RCP<N_TIA_DataStore> dsPtr = anaIntPtr_->getAnalysisMgr()->getTIADataStore();
      int numPoints = dsPtr->timeSteps.size();
      double startUpOffset = tiaMPDEParams_.initialTime;
      int lastJused = 0;
      for( int i=0; i < size_; i++ )
      {
        for( int j=lastJused; j < numPoints; j++ )
        {
          if( (fastTimes_[i] >= (dsPtr->timeSteps[j] - startUpOffset)) &&
              (fastTimes_[i] <  (dsPtr->timeSteps[j+1] - startUpOffset)) )
          {
            // found time points in the first period that bracket this one
#ifdef Xyce_DEBUG_MPDE
            if( debugLevel > 0 )
            {
              Xyce::dout() << "For fastTime_[" << i << "] = " << fastTimes_[i] <<
                " Found bracketing times : " << (dsPtr->timeSteps[j] - startUpOffset) << ", "
                << (dsPtr->timeSteps[j+1] - startUpOffset)
                << " at j = " << j;
            }
#endif

            // assume we're closest to the first point, but check if we're
            // really closest to the second one
            int interpolationPoint = j;
            if( fabs( fastTimes_[i] - (dsPtr->timeSteps[j] - startUpOffset)) >
                fabs( fastTimes_[i] - (dsPtr->timeSteps[j+1] - startUpOffset)) )
            {
              interpolationPoint = j+1;
            }
#ifdef Xyce_DEBUG_MPDE
            if( debugLevel > 0 )
            {
              Xyce::dout() << " closest to point " << interpolationPoint;
            }
#endif

            // now do the interpolation
            N_LAS_Vector * firstPeriodSolVecPtr = dsPtr->fastTimeSolutionVec[ interpolationPoint ];
            N_LAS_Vector * firstPeriodStateVecPtr = dsPtr->fastTimeStateVec[ interpolationPoint ];
            N_LAS_Vector * firstPeriodQVecPtr = dsPtr->fastTimeQVec[ interpolationPoint ];
            N_LAS_Vector * firstPeriodStoreVecPtr = dsPtr->fastTimeStoreVec[ interpolationPoint ];

            N_LAS_Vector * secondPeriodSolVecPtr = dsPtr->fastTimeSolutionVec[indicesUsed_[i]];
            N_LAS_Vector * secondPeriodStateVecPtr = dsPtr->fastTimeStateVec[indicesUsed_[i]];
            N_LAS_Vector * secondPeriodQVecPtr = dsPtr->fastTimeQVec[indicesUsed_[i]];
            N_LAS_Vector * secondPeriodStoreVecPtr = dsPtr->fastTimeStoreVec[indicesUsed_[i]];

            RCP<N_LAS_Vector> interpIcSolVecPtr = rcp(new N_LAS_Vector( *secondPeriodSolVecPtr ));
            RCP<N_LAS_Vector> interpIcStateVecPtr = rcp(new N_LAS_Vector( *secondPeriodStateVecPtr ));
            RCP<N_LAS_Vector> interpIcQVecPtr = rcp(new N_LAS_Vector( *secondPeriodQVecPtr ));
            RCP<N_LAS_Vector> interpIcStoreVecPtr = rcp(new N_LAS_Vector( *secondPeriodStoreVecPtr ));

            double fraction = fastTimes_[i] / period_;
#ifdef Xyce_DEBUG_MPDE
            if( debugLevel > 0 )
            {
              Xyce::dout() << " fraction = " << fraction << std::endl;
            }
#endif

            interpIcSolVecPtr->putScalar( 0.0 );
            interpIcStateVecPtr->putScalar( 0.0 );
            interpIcQVecPtr->putScalar( 0.0 );
            interpIcStoreVecPtr->putScalar( 0.0 );
            interpIcSolVecPtr->linearCombo( fraction, *firstPeriodSolVecPtr, (1.0-fraction), *secondPeriodSolVecPtr );
            interpIcStateVecPtr->linearCombo( fraction, *firstPeriodStateVecPtr, (1.0-fraction), *secondPeriodStateVecPtr );
            interpIcQVecPtr->linearCombo( fraction, *firstPeriodQVecPtr, (1.0-fraction), *secondPeriodQVecPtr );
            interpIcStoreVecPtr->linearCombo( fraction, *firstPeriodStoreVecPtr, (1.0-fraction), *secondPeriodStoreVecPtr );

            mpdeICVectorPtr_->block(i) = *interpIcSolVecPtr;
            mpdeICStateVectorPtr_->block(i) = *interpIcStateVecPtr;
            mpdeICQVectorPtr_->block(i) = *interpIcQVecPtr;
            mpdeICStoreVectorPtr_->block(i) = *interpIcStoreVecPtr;

            lastJused = interpolationPoint;
            break;  // break out of for(j..) loop.  We don't need to cycle more
          }
        }
      }

      // We've basically skipped a period with this initial condition.  So
      // have MPDE think it's starting a bit later.
      tiaMPDEParams_.initialTime += period_;

      // stil under development RLS
      //double error = checkPeriodicity_();
      //Xyce::dout() << "N_MPDE_Manager::runInitialConditon.  Based on two period anaysis, non-periodic error in this problem is " << error << std::endl;
    }
    break;

  default:
    msg = "N_MPDE_Manager::runInitialCondition_:";
    msg += "  Invalid IC option specified.\n";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL_0,msg);
    break;
  }

  // Assemble global block vector
  mpdeICVectorPtr_->assembleGlobalVector();
  mpdeICStateVectorPtr_->assembleGlobalVector();
  mpdeICQVectorPtr_->assembleGlobalVector();
  mpdeICStoreVectorPtr_->assembleGlobalVector();

  if (warpMPDE_)
  {
    //std::cout << "Setting omega = 1 into MPDE solution vector." << std::endl;
    // tscoffe/tmei 08/10/05:  put omega = 1 into MPDE solution vector
    // This is true for MPDE and WaMPDE because in WaMPDE we scale omega by T2.
    // tscoffe 12/14/06:  put phi(0) = omega(0) into MPDE solution vector.
    int omegaGID = warpMPDEPhasePtr_->getOmegaGID();
    int phiGID = warpMPDEPhasePtr_->getPhiGID();
    int omegaLID = mpdeICVectorPtr_->pmap()->globalToLocalIndex(omegaGID);
    int phiLID = mpdeICVectorPtr_->pmap()->globalToLocalIndex(phiGID);
    if (omegaLID >= 0)
      (*mpdeICVectorPtr_)[omegaLID] = 1.0;
    if (phiLID >= 0)
      (*mpdeICVectorPtr_)[phiLID] = 0.0;
  }

#ifdef Xyce_DEBUG_MPDE
  if ( debugLevel > 1 )
  {
    Xyce::dout() << "MPDE Initial Condition Solution!\n";
    mpdeICVectorPtr_->printPetraObject(std::cout);
    Xyce::dout() << "MPDE Initial Condition State Vector!\n";
    mpdeICStateVectorPtr_->printPetraObject(std::cout);
    Xyce::dout() << "MPDE Initial Condition Store Vector!\n";
    mpdeICStoreVectorPtr_->printPetraObject(std::cout);
  }
#endif // Xyce_DEBUG_MPDE

  // Now that the initial condition is over, turn off voltage limiting in
  // the device package (in case it is on).
  // tscoffe 1/17/07 We should experiment with turning this on after the state vector is blockized.
  //
  // ERK. 1/25/07.  I had to revamp voltlim in the device package to make it
  // work with MPDE, but I think it does now.
  //devInterfacePtr_->unsetVoltageLimiterFlag ();

#ifdef Xyce_DEBUG_MPDE
  Xyce::dout() << Xyce::section_divider << std::endl;
  if (icExitFlag_) exit(0);
#endif // Xyce_DEBUG_MPDE

  return returnValue;
}


//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::runDCop_()
// Purpose       : Runs the dc op problem
// Special Notes :
// Scope         : private
// Creator       : Richard Schiek, 1437, Electrical and Micro Modeling
// Creation Date : 09/17/07
//-----------------------------------------------------------------------------
bool N_MPDE_Manager::runDCOP_()
{
  bool success = true;

  // Perform DCOP
  success = tiaMPDEIfacePtr_->runDCOP();

  // store the dc op results in case we need them later
  dcOpSolVecPtr_ = rcp(new N_LAS_Vector( *(tiaMPDEIfacePtr_->getFinalSolution()) ));
  dcOpStateVecPtr_ = rcp(new N_LAS_Vector( *(tiaMPDEIfacePtr_->getStateFinalSolution()) ));
  dcOpStoreVecPtr_ = rcp(new N_LAS_Vector( *(tiaMPDEIfacePtr_->getStoreFinalSolution()) ));
  dcOpQVecPtr_ = rcp(new N_LAS_Vector( *(tiaMPDEIfacePtr_->getQVectorFinalSolution()) ));

  return success;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::runStartupPeriods_()
// Purpose       : Runs normal transient problem through the requested
//                 number of startup periods
// Special Notes :
// Scope         : private
// Creator       : Richard Schiek, 1437, Electrical and Micro Modeling
// Creation Date : 09/17/07
//-----------------------------------------------------------------------------
bool N_MPDE_Manager::runStartupPeriods_()
{
  bool returnValue = true;

  // try and advance time by startUpPeriods_ * period_

  // turn off the device manager's fast source flag
  //bool fastSourceSetFlagValue=devInterfacePtr_->getFastSourceSetFlag();
  //devInterfacePtr_->setFastSourceSetFlag(false);

  // need to advance time by startUpPeriods_ * period_
  N_TIA_TIAParams & tiaParams = tiaMPDEIfacePtr_->getTIAParams();
  N_TIA_TIAParams tiaParamsSave = tiaParams;

  // set DAE initial time = 0.0
  tiaParams.initialTime = 0.0;
  tiaParams.finalTime = startUpPeriods_ * period_;
  tiaParams.pauseTime = tiaParams.finalTime;
//  tiaParams.maxOrder = 5;
  tiaMPDEIfacePtr_->registerTIAParams(tiaParams); // redundant?

#ifdef Xyce_DEBUG_MPDE
  Xyce::dout() << "Advancing time through " << startUpPeriods_ << " startup periods"
     << " initialTime = " << tiaParams.initialTime
     << " finalTime = " << tiaParams.finalTime << std::endl;
#endif // Xyce_DEBUG_MPDE

  // tell the output manager to save this data as netlist.cir.startup.prn
  outMgrPtr_->setOutputFilenameSuffix(".startup");

  //runDCOP_();
  //tiaMPDEIfacePtr_->runTransient();

  if (warpMPDE_)
  {
    returnValue = tiaMPDEIfacePtr_->runTransient();
  }
  else
  {
    returnValue = tiaMPDEIfacePtr_->runTransientWithDCOP();
  }

  outMgrPtr_->finishOutput();

  // reset the output filename suffix
  outMgrPtr_->setOutputFilenameSuffix("");

  // put the dsPtr->currentSolutionPtr into dcOpSol and State Vec so that it
  // is used as our initial condition for the pending fast time scale runs
  dcOpSolVecPtr_ = rcp( new N_LAS_Vector( *(tiaMPDEIfacePtr_->getFinalSolution()) ));
  dcOpStateVecPtr_ = rcp( new N_LAS_Vector( *(tiaMPDEIfacePtr_->getStateFinalSolution()) ));
  dcOpQVecPtr_ = rcp( new N_LAS_Vector( *(tiaMPDEIfacePtr_->getQVectorFinalSolution()) ));
  dcOpStoreVecPtr_ = rcp( new N_LAS_Vector( *(tiaMPDEIfacePtr_->getStoreFinalSolution()) ));

  // tell mpde to start after this startup period
  tiaParamsSave.initialTime = startUpPeriods_ * period_;
  tiaMPDEParams_.initialTime = startUpPeriods_ * period_;
  tiaMPDEIfacePtr_->registerTIAParams(tiaParamsSave);
  startUpPeriodsFinished_ = true;

  return returnValue;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::runTransientIC_
// Purpose       : Conducts a regular transient run for MPDE initial conditions
// Special Notes :
// Scope         : private
// Creator       : Richard Schiek, 1437, Electrical and Micro Modeling
// Creation Date : 09/17/07
//-----------------------------------------------------------------------------
bool N_MPDE_Manager::runTransientIC_()
{
  bool returnValue = false;

#ifdef Xyce_DEBUG_MPDE
  Xyce::dout() << "N_MPDE_Manager::runTransientIC_() Doing initial transient run." << std::endl;
#endif
  // flag for N_TIA_ControlgAlgorithm that we're calculating an IC
  // this prevents extra DC op data from being printed.
  mpdeICFlag_ = true;

  if (!warpMPDE_)
  {
    if( is_null(dcOpSolVecPtr_) )
    {
      // haven't done a dc op yet, so do it now
      runDCOP_();
    }
  }
  
  // use an initial transient run to create a set of time points
  // for the fast time scale
  N_TIA_TIAParams & tiaParams = tiaMPDEIfacePtr_->getTIAParams();
  N_TIA_TIAParams tiaParamsSave = tiaParams;

  // now try the real thing
  tiaParams.initialTime = tiaMPDEParams_.initialTime;
#ifdef MPDE_OLD_TWO_PERIOD
  if( initialCondition_ == MPDE_IC_TWO_PERIOD )
  {
    // in this case we need to integrate forward 2 periods
    tiaParams.finalTime = tiaMPDEParams_.initialTime + 2.0 * period_;
  }
  else
  {
    tiaParams.finalTime = tiaMPDEParams_.initialTime + period_;
  }
#else
  tiaParams.finalTime = tiaMPDEParams_.initialTime + period_;
#endif

  tiaParams.pauseTime = tiaParams.finalTime;
  tiaParams.saveTimeStepsFlag = true;
//  tiaParams.maxOrder = 2;
//  tiaParams.integrationMethod = 7;
//  tiaParams.relErrorTol = 1e-3;

  // register new tiaParams with time integrator
  tiaMPDEIfacePtr_->registerTIAParams(tiaParams); // redundant?
//  anaIntPtr_->initializeAll();

//  if (!warpMPDE_)
  if (!is_null(dcOpSolVecPtr_))
  {
    tiaMPDEIfacePtr_->setInitialCondition(dcOpSolVecPtr_.get());
    tiaMPDEIfacePtr_->setStateInitialCondition(dcOpStateVecPtr_.get());
    tiaMPDEIfacePtr_->setQVectorInitialCondition(dcOpQVecPtr_.get());
    tiaMPDEIfacePtr_->setStoreInitialCondition(dcOpStoreVecPtr_.get());
  }

  outMgrPtr_->prepareOutput( Xyce::Analysis::ANP_MODE_TRANSIENT, std::vector<N_ANP_SweepParam>(), std::vector<N_ANP_SweepParam>() ); //, *dsPtr_->currSolutionPtr, *dsPtr_->currStatePtr, *dsPtr_->currStorePtr );
  
  {
    Xyce::IO::OutputMgr::ActiveOutput x(*outMgrPtr_);
    
    if (saveIcData_)
    {
      // tell the output manager to append future simulation data.  This
      // will let us keep the initial condition data
//      outMgrPtr_->setOutputFilenameSuffix( ".mpde_ic" );
      x.add(Xyce::IO::PrintType::MPDE_IC);
    }
    
    // success for this function is the success/failure of the transient run
    returnValue = tiaMPDEIfacePtr_->runTransient();
//  returnValue = tiaMPDEIfacePtr_->runTransientWithDCOP();

  }
    
#ifndef MPDE_OLD_TWO_PERIOD
  // If this is a two period initial condition calculation, make sure
  // the second period uses the same time steps as the first
  if( initialCondition_ == MPDE_IC_TWO_PERIOD )
  {
    if(saveIcData_)
    {
      // tell the output manager to append future simulation data.  This
      // will let us keep the initial condition data
      outMgrPtr_->setOutputFilenameSuffix( ".mpde_ic2" );
    }
    RCP<N_TIA_DataStore> dsPtr = anaIntPtr_->getAnalysisMgr()->getTIADataStore();
    int numPoints = dsPtr->timeSteps.size();

    // try using the 2-level machinery to take prescribed time steps.
    N_TIA_TimeIntInfo tiInfo;
    N_TIA_TwoLevelError tlError;

    for( int ts = 1; ts < numPoints; ts++)
    {
      tiaMPDEIfacePtr_->getTimeIntInfo( tiInfo );
      tiInfo.nextTime = dsPtr->timeSteps[ts] + period_;
      tiInfo.nextTimeStep = dsPtr->timeSteps[ts] - dsPtr->timeSteps[ts-1];
      tiInfo.finalTime = dsPtr->timeSteps[ts] + period_;
      tiInfo.beginIntegrationFlag = false;

      tiaMPDEIfacePtr_->runStep( tiInfo, tlError);
      tiaMPDEIfacePtr_->stepSuccess( 1 );
    }

    /*
      Xyce::dout() << "From one period there were " << numPoints << " saved time steps." << std::endl;
      numPoints = dsPtr->timeSteps.size();
      Xyce::dout() << "Now there are " << numPoints << " time steps saved." << std::endl;
      for( int ts=0; ts<numPoints; ts++ )
      {
      Xyce::dout() << "timeSteps[ " << ts << " ] = " << dsPtr->timeSteps[ts] << std::endl;
      }
    */
    // while runTransient() calls finishOutput(), runStep() does not.  So we
    // need to finish the output of the second part of the IC calculation
    outMgrPtr_->finishOutput();
  }
#endif

  // restore the saved copy of time integration parameters.
  tiaParamsSave.saveTimeStepsFlag = false;
  anaIntPtr_->registerLoader(&*appLoaderPtr_);
  tiaParamsSave.initialTime += period_;  // start MPDE problem after this transient init.
  tiaMPDEIfacePtr_->registerTIAParams(tiaParamsSave);
  mpdeICFlag_ = false;  // done calculating ic stuff
  return returnValue;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::filterFastTimePoints_()
// Purpose       : Tries to filter the fast time points from a transient run
//                 so that points are not too close together
// Special Notes :
// Scope         : private
// Creator       : Richard Schiek, 1437, Electrical and Micro Modeling
// Creation Date : 09/17/07
//-----------------------------------------------------------------------------
bool N_MPDE_Manager::filterFastTimePoints_()
{
  RCP<N_TIA_DataStore> dsPtr = anaIntPtr_->getAnalysisMgr()->getTIADataStore();
  int numPoints = dsPtr->timeSteps.size();
  int startIndex = 0;
  double extraPeriodShift = 0.0;

  // Note, if we integrated more then one period in the runTransientIC_()
  // function, the we'll try to pull in points only from the second period.
  if( initialCondition_ == MPDE_IC_TWO_PERIOD )
  {
    extraPeriodShift = period_;
    for(int i=0; i<numPoints; i++ )
    {
      if(dsPtr->timeSteps[i] > (tiaMPDEParams_.initialTime + period_))
      {
        startIndex = i;
        break;
      }
    }
  }

#ifdef Xyce_DEBUG_MPDE
  Xyce::dout() << "Initial transient run produced " << numPoints << " points." << std::endl;
#endif

  std::list<int> goodIndicies;
  int numGoodPoints = 0;
  int breakpoints = 0;          // need to keep track of how many breakpoints there are

  if( !maxCalcSizeGiven_  || (numPoints < maxCalcSize_ ) )
  {
    // keep all of the points from the transient run

    for( int i=startIndex; i < numPoints; i++ )
    {
      goodIndicies.push_back( i );
    }
    numGoodPoints = goodIndicies.size();

#ifdef Xyce_DEBUG_MPDE
    Xyce::dout() << " keeping all points for MPDE calculations." << std::endl;
#endif
  }
  else
  {
    // locate the set of points and indicies that are "good", i.e. those that
    // are not too close together
    const double relTimeTol = 0.005;
    const double absTimeDiffTol = 5.0e-9;

    // always keep first point
    goodIndicies.push_back(startIndex);
    int lastGoodIndex = startIndex;

    for( int i=startIndex + 1; i < numPoints; i++ )
    {
      // count up breakpoints
      if( dsPtr->timeStepsBreakpointFlag[i] == true )
      {
        breakpoints++;
      }

#ifdef Xyce_DEBUG_MPDE
      if( debugLevel > 0 )
      {
        Xyce::dout() << "\t\t timeStep[ " << i << " ] = " << dsPtr->timeSteps[i];
        if( dsPtr->timeStepsBreakpointFlag[i] == true )
        {
          Xyce::dout() << "  Breakpoint";
        }
        Xyce::dout() << std::endl;
      }
#endif

      double delta_t = (dsPtr->timeSteps[i] - dsPtr->timeSteps[lastGoodIndex]);
      if( ((delta_t > relTimeTol * period_) && (delta_t > absTimeDiffTol)) ||  dsPtr->timeStepsBreakpointFlag[i] == true )
      {
        // found a good point so save the index
        goodIndicies.push_back( i );
        lastGoodIndex = i;
      }
    }
    numGoodPoints = goodIndicies.size();

#ifdef Xyce_DEBUG_MPDE
    Xyce::dout() << " of which " << numGoodPoints << " fit tolerance requirements of abstol = "
      << absTimeDiffTol << " and reltol*period = " << relTimeTol * period_ <<std::endl;
#endif
  }

  if( (maxCalcSizeGiven_)  && (maxCalcSize_ < numGoodPoints ) )
  {
    // can't just copy all of the fast times,
    // filter down to maxCalcSizeGiven
    //vector<double> tmpTimePoints = dsPtr->timeSteps;
    fastTimes_.resize( maxCalcSize_ + 1);
    indicesUsed_.resize( maxCalcSize_ );
    // always grab the first and last point
    fastTimes_[0] = dsPtr->timeSteps[startIndex] - tiaMPDEParams_.initialTime - extraPeriodShift;
    indicesUsed_[0] = startIndex;
    fastTimes_[maxCalcSize_] = period_;

    // this is a little confusing but here is what we're trying to do
    // we have a list of indicies that are considered good points.
    // we need to sample this list a regular intervals to get it
    // into the max size requested by the users.  We also have
    // points that are breakpoints.  We must keep all the breakpoints.

    double sampleRate = (1.0*(numGoodPoints-1-breakpoints)) / (1.0*(maxCalcSize_-1));

    std::list<int>::iterator currentIndex = goodIndicies.begin();
    std::list<int>::iterator endIndex = goodIndicies.end();

    int indexCounter=0, k=0;
    while( currentIndex != endIndex )
    {
      int targetIndex = static_cast<int>( k * sampleRate + 0.5 );
      if( ((indexCounter == targetIndex) && (k < maxCalcSize_)) ||
        (dsPtr->timeStepsBreakpointFlag[*currentIndex] == true) )
      {
        indicesUsed_[k] = *currentIndex;
#ifdef Xyce_DEBUG_MPDE
        if( debugLevel > 0 )
        {
          Xyce::dout() << " indicesUsed_[" << k << "] = " << indicesUsed_[k] << std::endl;
        }
#endif
        fastTimes_[k] = dsPtr->timeSteps[*currentIndex] - tiaMPDEParams_.initialTime - extraPeriodShift;
        k++;
      }
      currentIndex++;
      indexCounter++;
    }
  }
  else
  {
    // number of good points is less than the requested number
    // of points so just copy them all over
    //
    fastTimes_.resize( numGoodPoints + 1 );
//    indiciesUsed_.resize( numGoodPoints );
//    fastTimes_[0] = dsPtr->timeSteps[startIndex] - tiaMPDEParams_.initialTime - extraPeriodShift;

    indicesUsed_.resize( numGoodPoints );
//    fastTimes_[0] = dsPtr->timeSteps[startIndex] - tiaMPDEParams_.initialTime - extraPeriodShift;
    fastTimes_[0] = 0.0;
    indicesUsed_[0] = startIndex;

    fastTimes_[numGoodPoints] = period_;

    std::list<int>::iterator currentIndex = goodIndicies.begin();
    std::list<int>::iterator endIndex = goodIndicies.end();
//    int k=0;
    int k=1;
    currentIndex++;
    while( currentIndex != endIndex )
    {
      indicesUsed_[k] = *currentIndex;
#ifdef Xyce_DEBUG_MPDE
      if( debugLevel > 0 )
      {
        Xyce::dout() << " indicesUsed_[" << k << "] = " << indicesUsed_[k] << std::endl;
      }
#endif
      fastTimes_[k] = dsPtr->timeSteps[*currentIndex] - tiaMPDEParams_.initialTime - extraPeriodShift;
      k++;
      currentIndex++;
    }
  }
  size_ = fastTimes_.size() - 1;

//#ifdef Xyce_DEBUG_MPDE
  Xyce::dout() << "MPDE: " << size_ << " fast time points added to the problem." << std::endl
    << "fast time point range: " << fastTimes_[0] << ", " << fastTimes_[size_] << std::endl;
  if( debugLevel > 0 )
  {
    Xyce::dout() << "MPDE: new fast times are:" << std::endl;
    int i=0;
    for( i = 0; i < size_; ++i )
      Xyce::dout() << "fastTimes_["<<i<<"] = " << fastTimes_[i]
        << " forward difference is " << (fastTimes_[i+1] - fastTimes_[i]) << std::endl;
    Xyce::dout() << "period = " << period_
      << " fastTimes_[numGoodPoints] = " << fastTimes_[numGoodPoints]
      << " extraPeriodShift = " << extraPeriodShift << std::endl;
  }
//#endif // Xyce_DEBUG_MPDE
  return true;
}


//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::setupMPDEProblem_
// Purpose       : Configures solvers, etc. for MPDE run
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/19/04
//-----------------------------------------------------------------------------
bool N_MPDE_Manager::setupMPDEProblem_()
{
#ifdef Xyce_DEBUG_MPDE
  Xyce::dout() << std::endl
               << Xyce::section_divider << std::endl
               << "  N_MPDE_Manager::setupMPDEProblem" << std::endl;
#endif // Xyce_DEBUG_MPDE

  Xyce::lout() << " ***** Setting up full MPDE problem....\n" << std::endl;

  //Destroy Solvers, etc. from IC phase
  //-----------------------------------------
  anaIntPtr_->resetAll();
  nlsMgrPtr_ = Teuchos::null;
  //-----------------------------------------

  //Finish setup of MPDE Builder
  //-----------------------------------------
  mpdeBuilderPtr_->generateGraphs
    ( *(pdsMgrPtr_->getMatrixGraph( "JACOBIAN" )),
      *(pdsMgrPtr_->getMatrixGraph( "JACOBIAN" )),
      *(pdsMgrPtr_->getMatrixGraph( "JACOBIAN" )) );
  //-----------------------------------------

  //Setup MPDE Loader
  //-----------------------------------------
  RCP<N_MPDE_Manager> mpdeMgrRCPtr = rcp(this,false);
  mpdeLoaderPtr_ = rcp(new N_MPDE_Loader( mpdeState_, mpdeDiscPtr_, mpdeMgrRCPtr, warpMPDEPhasePtr_ ));

  mpdeLoaderPtr_->registerAppLoader( appLoaderPtr_ );
  mpdeLoaderPtr_->registerMPDEDeviceInterface( mpdeDeviceInterfacePtr_ );
  mpdeLoaderPtr_->setFastTimes( fastTimes_ );
  mpdeLoaderPtr_->setPeriodFlags( nonPeriodicFlags );

  //Include these for the loader to use with devices
  mpdeLoaderPtr_->registerAppNextVec ( rcp(appBuilderPtr_->createVector()) );
  mpdeLoaderPtr_->registerAppCurrVec ( rcp(appBuilderPtr_->createVector()) );
  mpdeLoaderPtr_->registerAppLastVec ( rcp(appBuilderPtr_->createVector()) );

  mpdeLoaderPtr_->registerAppNextStaVec ( rcp(appBuilderPtr_->createStateVector()) );
  mpdeLoaderPtr_->registerAppCurrStaVec ( rcp(appBuilderPtr_->createStateVector()) );
  mpdeLoaderPtr_->registerAppLastStaVec ( rcp(appBuilderPtr_->createStateVector()) );
  mpdeLoaderPtr_->registerOmegadQdt2( rcp_dynamic_cast<N_LAS_BlockVector>(rcp(mpdeBuilderPtr_->createVector())) );

  mpdeLoaderPtr_->registerAppNextStoVec ( rcp(appBuilderPtr_->createStoreVector()) );
  mpdeLoaderPtr_->registerAppCurrStoVec ( rcp(appBuilderPtr_->createStoreVector()) );
  mpdeLoaderPtr_->registerAppLastStoVec ( rcp(appBuilderPtr_->createStoreVector()) );
  mpdeLoaderPtr_->registerAppStoLeadCurrQCompVec ( rcp(appBuilderPtr_->createStoreVector()) );

  mpdeLoaderPtr_->registerAppdQdx( rcp(appBuilderPtr_->createMatrix()) );
  mpdeLoaderPtr_->registerAppdFdx( rcp(appBuilderPtr_->createMatrix()) );

  mpdeLoaderPtr_->registerMPDEdQdx
    ( rcp_dynamic_cast<N_LAS_BlockMatrix>(rcp(mpdeBuilderPtr_->createMatrix())) );

  mpdeLoaderPtr_->registerMPDEdFdx
    ( rcp_dynamic_cast<N_LAS_BlockMatrix>(rcp(mpdeBuilderPtr_->createMatrix())) );

  //-----------------------------------------

  //Construct Solvers, etc. for MPDE Phase
  //-----------------------------------------
  nlsMgrPtr_ = rcp(new N_NLS_Manager(commandLine_));
  lasMPDESysPtr_ = rcp(new N_LAS_System());
  //-----------------------------------------

  // modify time integrator parameters for MPDE phase
  tiaMPDEParams_.userSpecified_startingTimeStep = period_;
  // if the user requested that a few steps be taken without LTE control,
  // then assert that option here.
  if( nonLteStepsGiven_ )
  {
    tiaMPDEParams_.errorAnalysisOption = 1;
    tiaMPDEParams_.errorAnalysisOptionResetIter = nonLteSteps_;
  }

  //Registration of Solvers, etc. for MPDE Phase
  //-----------------------------------------
  // Device package registrations.
  mpdeDeviceInterfacePtr_->registerNonlinearSolver( nlsMgrPtr_ );
  mpdeDeviceInterfacePtr_->registerAnalysisInterface( anaIntPtr_ );
  //let devices continue to use old lasSys object
  mpdeDeviceInterfacePtr_->registerLinearSystem( lasMPDESysPtr_ );

  //-----------------------------------------

  anaIntPtr_->registerTIAParams( tiaMPDEParams_ );
  anaIntPtr_->registerLinearSystem( &*lasMPDESysPtr_ );
  anaIntPtr_->registerNLSManager( &*nlsMgrPtr_ );
  anaIntPtr_->registerLoader( dynamic_cast<N_LOA_Loader*>( &*mpdeLoaderPtr_ ) );

  // Give NLS Manager the same old nonlinearEquationLoader as it just calls the TIA loader in newDAE
  nlsMgrPtr_->registerLoader( dynamic_cast<N_LOA_Loader*>( &*nonlinearEquationLoaderPtr_ ) );
  nlsMgrPtr_->registerLinearSystem( &*lasMPDESysPtr_ );
  nlsMgrPtr_->registerAnalysisInterface( &*anaIntPtr_ );
  nlsMgrPtr_->registerOutputMgr( &*outMgrPtr_ );
  nlsMgrPtr_->registerParallelMgr( &*pdsMgrPtr_ );

#ifdef Xyce_PARALLEL_MPI
  // Set reindexing for KLU
  N_UTL_OptionBlock LSOB = N_UTL_OptionBlock();
  LSOB.getParams().push_back( N_UTL_Param( "KLU_REINDEX", 1 ) );
  nlsMgrPtr_->setLinSolOptions( LSOB );
#endif

  //hack needed by TIA initialization currently
  mpdeBuilderPtr_->registerPDSManager( &*pdsMgrPtr_ );

  lasMPDESysPtr_->registerANPInterface( &*anaIntPtr_ );
  lasMPDESysPtr_->registerPDSManager( &*pdsMgrPtr_ );
  lasMPDESysPtr_->registerBuilder
    ( dynamic_cast<N_LAS_Builder*>(&*mpdeBuilderPtr_) );

  //need to cut out unnecessary stuff from this call for new dae
  lasMPDESysPtr_->initializeSystem();

  //-----------------------------------------

  //Initialization of Solvers, etc. for MPDE Phase
  //-----------------------------------------
  // register the saved MPDE transient options.
  anaIntPtr_->setTranOptions (saved_tranMPDEOB_);
  //Dummy call to setup time integrator for transient
  anaIntPtr_->setTranAnalysisParams( N_UTL_OptionBlock() );
  //need to cut out unnecessary stuff from this call for new dae
  anaIntPtr_->initializeAll();
  //need to cut out unnecessary stuff from this call for new dae
  nlsMgrPtr_->initializeAll();


#ifdef Xyce_DEBUG_MPDE
  Xyce::dout() << Xyce::section_divider << std::endl;
#endif // Xyce_DEBUG_MPDE
  return true;
}


//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::runMPDEProblem_
// Purpose       : Actual execution of MPDE problem
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/19/04
//-----------------------------------------------------------------------------
bool N_MPDE_Manager::runMPDEProblem_()
{
  bool returnValue = true;

#ifdef Xyce_DEBUG_MPDE
  Xyce::dout() << std::endl;
  Xyce::dout() << Xyce::section_divider << std::endl;
  Xyce::dout() <<
                          "  N_MPDE_Manager::runMPDEProblem_" << std::endl;
#endif // Xyce_DEBUG_MPDE
  tiaMPDEIfacePtr_->setInitialCondition(&*mpdeICVectorPtr_);
  tiaMPDEIfacePtr_->setStateInitialCondition(&*mpdeICStateVectorPtr_);
  tiaMPDEIfacePtr_->setQVectorInitialCondition(&*mpdeICQVectorPtr_);
  tiaMPDEIfacePtr_->setStoreInitialCondition(&*mpdeICStoreVectorPtr_);

  Xyce::lout() << " ***** Beginning full MPDE simulation....\n" << std::endl;

  // try to run the transient problem
  returnValue = tiaMPDEIfacePtr_->runTransient();
#ifdef Xyce_DEBUG_MPDE
  Xyce::dout() << Xyce::section_divider << std::endl;
#endif // Xyce_DEBUG_MPDE

  return returnValue;
}


//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::runTests_
// Purpose       : Testing MPDE objects
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, 9233, Computational Sciences
// Creation Date : 03/19/04
//-----------------------------------------------------------------------------
bool N_MPDE_Manager::runTests_()
{
  Xyce::dout() << "N_MPDE_Manager::runTests_\n";

  //Finish setup of MPDE Builder
  //-----------------------------------------
  Xyce::dout() << "N_MPDE_Manager::runTests_[Generate Graphs]\n";
  mpdeBuilderPtr_->generateGraphs
    ( *(pdsMgrPtr_->getMatrixGraph( "JACOBIAN" )),
      *(pdsMgrPtr_->getMatrixGraph( "JACOBIAN" )),
      *(pdsMgrPtr_->getMatrixGraph( "JACOBIAN" )) );

  Xyce::dout() << "N_MPDE_Manager::runTests_[Finished Graphs]\n";
  //-----------------------------------------

  //Setup MPDE Loader
  //-----------------------------------------
  Xyce::dout() << "N_MPDE_Manager::runTests_[Construct Loader]\n";
  RCP<N_MPDE_Manager> mpdeMgrRCPtr = rcp(this,false);
  mpdeLoaderPtr_ = rcp(new N_MPDE_Loader( mpdeState_, mpdeDiscPtr_, mpdeMgrRCPtr, warpMPDEPhasePtr_ ));
  mpdeLoaderPtr_->registerAppLoader( appLoaderPtr_ );
  mpdeLoaderPtr_->registerMPDEDeviceInterface( mpdeDeviceInterfacePtr_ );
  mpdeLoaderPtr_->setFastTimes( fastTimes_ );

  //Include these for the loader to use with devices
  mpdeLoaderPtr_->registerAppNextVec( rcp(appBuilderPtr_->createVector()) );
  mpdeLoaderPtr_->registerAppCurrVec( rcp(appBuilderPtr_->createVector()) );
  mpdeLoaderPtr_->registerAppLastVec( rcp(appBuilderPtr_->createVector()) );

  mpdeLoaderPtr_->registerAppNextStaVec( rcp(appBuilderPtr_->createStateVector()) );
  mpdeLoaderPtr_->registerAppCurrStaVec( rcp(appBuilderPtr_->createStateVector()) );
  mpdeLoaderPtr_->registerAppLastStaVec( rcp(appBuilderPtr_->createStateVector()) );
  mpdeLoaderPtr_->registerAppNextStoVec( rcp(appBuilderPtr_->createStoreVector()) );
  mpdeLoaderPtr_->registerAppCurrStoVec( rcp(appBuilderPtr_->createStoreVector()) );
  mpdeLoaderPtr_->registerAppLastStoVec( rcp(appBuilderPtr_->createStoreVector()) );
  mpdeLoaderPtr_->registerAppStoLeadCurrQCompVec( rcp(appBuilderPtr_->createStoreVector()) );
  mpdeLoaderPtr_->registerAppdQdx( rcp(appBuilderPtr_->createMatrix()) );
  mpdeLoaderPtr_->registerAppdFdx( rcp(appBuilderPtr_->createMatrix()) );
  Xyce::dout() << "N_MPDE_Manager::runTests_[Finished Loader]\n";
  //-----------------------------------------

  //Test Construction of "MPDE Size" Objects
  //-----------------------------------------
  Xyce::dout() << "N_MPDE_Manager::runTests_[Create Vectors]\n";
  RCP<N_LAS_Vector> Q = rcp(mpdeBuilderPtr_->createVector());
  RCP<N_LAS_Vector> F = rcp(mpdeBuilderPtr_->createVector());
  RCP<N_LAS_Vector> Res = rcp(mpdeBuilderPtr_->createVector());
  RCP<N_LAS_Vector> nextX = rcp(mpdeBuilderPtr_->createVector());
  RCP<N_LAS_Vector> currX = rcp(mpdeBuilderPtr_->createVector());
  RCP<N_LAS_Vector> lastX = rcp(mpdeBuilderPtr_->createVector());


  RCP<N_LAS_Vector> nextS = rcp(mpdeBuilderPtr_->createStateVector());
  RCP<N_LAS_Vector> currS = rcp(mpdeBuilderPtr_->createStateVector());
  RCP<N_LAS_Vector> lastS = rcp(mpdeBuilderPtr_->createStateVector());
  RCP<N_LAS_Vector> dSdt = rcp(mpdeBuilderPtr_->createStateVector());
  RCP<N_LAS_Vector> nextStore = rcp(mpdeBuilderPtr_->createStoreVector());
  RCP<N_LAS_Vector> currStore = rcp(mpdeBuilderPtr_->createStoreVector());
  RCP<N_LAS_Vector> lastStore = rcp(mpdeBuilderPtr_->createStoreVector());
  RCP<N_LAS_Vector> storeLeadCurrQComp = rcp(mpdeBuilderPtr_->createStoreVector());

  RCP<N_LAS_BlockVector> bQ = rcp_dynamic_cast<N_LAS_BlockVector>(Q);
  RCP<N_LAS_BlockVector> bF = rcp_dynamic_cast<N_LAS_BlockVector>(F);
  RCP<N_LAS_BlockVector> bRes = rcp_dynamic_cast<N_LAS_BlockVector>(Res);
  RCP<N_LAS_BlockVector> bNextX = rcp_dynamic_cast<N_LAS_BlockVector>(nextX);
  RCP<N_LAS_BlockVector> bCurrX = rcp_dynamic_cast<N_LAS_BlockVector>(currX);
  RCP<N_LAS_BlockVector> bLastX = rcp_dynamic_cast<N_LAS_BlockVector>(lastX);

  RCP<N_LAS_BlockVector> bNextS = rcp_dynamic_cast<N_LAS_BlockVector>(nextS);
  RCP<N_LAS_BlockVector> bCurrS = rcp_dynamic_cast<N_LAS_BlockVector>(currS);
  RCP<N_LAS_BlockVector> bLastS = rcp_dynamic_cast<N_LAS_BlockVector>(lastS);
  RCP<N_LAS_BlockVector> bdSdt = rcp_dynamic_cast<N_LAS_BlockVector>(dSdt);
  RCP<N_LAS_BlockVector> bNextStore = rcp_dynamic_cast<N_LAS_BlockVector>(nextStore);
  RCP<N_LAS_BlockVector> bCurrStore = rcp_dynamic_cast<N_LAS_BlockVector>(currStore);
  RCP<N_LAS_BlockVector> bLastStore = rcp_dynamic_cast<N_LAS_BlockVector>(lastStore);

  bNextX->printPetraObject(Xyce::dout());
  Xyce::dout() << "N_MPDE_Manager::runTests_[Finished Vectors]\n";

  Xyce::dout() << "N_MPDE_Manager::runTests_[Construct Matrices]\n";
  RCP<N_LAS_Matrix> dQdx = rcp(mpdeBuilderPtr_->createMatrix());
  RCP<N_LAS_Matrix> dFdx = rcp(mpdeBuilderPtr_->createMatrix());
  RCP<N_LAS_Matrix> Jac  = rcp(mpdeBuilderPtr_->createMatrix());

  RCP<N_LAS_BlockMatrix> bdQdx = rcp_dynamic_cast<N_LAS_BlockMatrix>(dQdx);
  RCP<N_LAS_BlockMatrix> bdFdx = rcp_dynamic_cast<N_LAS_BlockMatrix>(dFdx);
  RCP<N_LAS_BlockMatrix> bJac = rcp_dynamic_cast<N_LAS_BlockMatrix>(Jac);

  bJac->printPetraObject(Xyce::dout());
  Xyce::dout() << "N_MPDE_Manager::runTests_[Finished Matrices]\n";

  //-----------------------------------------

  //Test Load of MPDE Objects
  //-----------------------------------------
  Xyce::dout() << "N_MPDE_Manager::runTests_[Load Vectors]\n";
  mpdeLoaderPtr_->loadDAEVectors(
      &*nextX, &*currX, &*lastX,
      &*nextS, &*currS, &*lastS, &*dSdt,
      &*nextStore, &*currStore, &*lastStore, &*storeLeadCurrQComp, &*Q, &*F, 0, 0 );
  Xyce::dout() << "N_MPDE_Manager::runTests_[Load Matrices]\n";

  mpdeLoaderPtr_->loadDAEMatrices( &*nextX, &*nextS, &*dSdt, &*nextStore, &*dQdx, &*dFdx );
  Xyce::dout() << "N_MPDE_Manager::runTests_[Finished Loads]\n";

  bQ->printPetraObject(Xyce::dout());
  bF->printPetraObject(Xyce::dout());

  bdQdx->printPetraObject(Xyce::dout());
  bdFdx->printPetraObject(Xyce::dout());
  //-----------------------------------------

  //Pretend to solve DCOP problems
  //-----------------------------------------
  Res->linearCombo( 1.0, *Res, +1.0, *F );
  bRes->printPetraObject(Xyce::dout());

  Jac->add( *dFdx );
  bJac->printPetraObject(Xyce::dout());

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_MPDE_Manager::checkPeriodicity
// Purpose       : Return's an L2 norm of the difference of the initial
//                 condition signal's over two periods (if 2 were
//                 calculated) or just the endpoints if only one
//                 period was calculated.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, 1437, Electrical and Microsystem Simulation
// Creation Date : 06/13/2008
//-----------------------------------------------------------------------------
double  N_MPDE_Manager::checkPeriodicity_()
{
  double returnValue = 0.0;

  // get access to time history of simulation
  RCP<N_TIA_DataStore> dsPtr = anaIntPtr_->getAnalysisMgr()->getTIADataStore();
  int numPoints = dsPtr->timeSteps.size();

  // Note, if we integrated more then one period in the runTransientIC_()
  // function, the we'll try to pull in points only from the second period.
  if( initialCondition_ == MPDE_IC_TWO_PERIOD )
  {
    double startUpOffset = tiaMPDEParams_.initialTime - period_;
    int lastJused = 0;
    for( int i=0; i < size_; i++ )
    {
      for( int j=lastJused; j < numPoints; j++ )
      {
        if( (fastTimes_[i] >= (dsPtr->timeSteps[j] - startUpOffset)) &&
            (fastTimes_[i] <  (dsPtr->timeSteps[j+1] - startUpOffset)) )
        {
          // found time points in the first period that bracket this one
#ifdef Xyce_DEBUG_MPDE
          if( debugLevel > 0 )
          {
            Xyce::dout() << "For fastTime_[" << i << "] = " << fastTimes_[i] <<
              " Found bracketing times : " << (dsPtr->timeSteps[j] - startUpOffset) << ", "
              << (dsPtr->timeSteps[j+1] - startUpOffset)
              << " at j = " << j;
          }
#endif

          // assume we're closest to the first point, but check if we're
          // really closest to the second one
          int interpolationPoint = j;
          if( fabs( fastTimes_[i] - (dsPtr->timeSteps[j] - startUpOffset)) >
              fabs( fastTimes_[i] - (dsPtr->timeSteps[j+1] - startUpOffset)) )
          {
            interpolationPoint = j+1;
          }
#ifdef Xyce_DEBUG_MPDE
          if( debugLevel > 0 )
          {
            Xyce::dout() << " closest to point " << interpolationPoint << " index used = " << indicesUsed_[i] << std::endl;
          }
#endif
          N_LAS_Vector *thisPeriod = dsPtr->fastTimeSolutionVec[indicesUsed_[i]];
          N_LAS_Vector *lastPeriod = dsPtr->fastTimeSolutionVec[interpolationPoint];
          N_LAS_Vector scratchVec( *thisPeriod );
          scratchVec.linearCombo( 1.0, scratchVec, -1.0, *lastPeriod );
          scratchVec.infNorm(&returnValue);
          Xyce::dout() << i << " returnValue = " << returnValue << std::endl;
        }
      }
    }
  }
  return returnValue;
}


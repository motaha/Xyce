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
// Filename       : $RCSfile: N_DEV_DeviceMgr.C,v $
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
// Revision Number: $Revision: 1.557.2.6 $
//
// Revision Date  : $Date: 2013/10/03 17:23:38 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

#ifdef Xyce_DEBUG_DEVICE
#ifdef HAVE_CSTDIO
#include <cstdio>
#else
#include <stdio.h>
#endif
#endif

#include <algorithm>
#include <sstream>

// ----------   Xyce Includes   ----------
#include <N_DEV_Const.h>
#include <N_DEV_DeviceMgr.h>

#include <N_DEV_DeviceBld.h>

#include <N_DEV_PlaceHolder.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_Source.h>
#include <N_DEV_SourceData.h>
#include <N_DEV_Vsrc.h>
#include <N_DEV_DAC.h>
#include <N_DEV_ADC.h>
#include <N_DEV_Xygra.h>

#ifdef Xyce_EXTDEV
#include <N_DEV_ExternDevice.h>
#endif

#include <N_DEV_DeviceModel.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_DeviceOptions.h>

#include <N_LAS_System.h>
#include <N_LAS_Builder.h>
#include <N_LAS_MultiVector.h>
#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>


#include <N_ANP_AnalysisInterface.h>
#include <N_NLS_Manager.h>
#include <N_TIA_TimeIntInfo.h>
#include <N_NLS_NonLinInfo.h>
#include <N_TIA_TwoLevelError.h>

#include <N_UTL_Misc.h>
#include <N_UTL_OptionBlock.h>
#include <N_UTL_BreakPoint.h>
#include <N_UTL_Expression.h>

#include <N_IO_PkgOptionsMgr.h>

#include <N_ERH_ErrorMgr.h>

#include <N_MPDE_Manager.h>

#include <N_PDS_Manager.h>
#include <N_PDS_Comm.h>

#include <N_IO_CmdParse.h>
#include <N_IO_OutputMgr.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::factory
// Purpose       : factory function
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
DeviceMgr * DeviceMgr::factory (N_IO_CmdParse & cp)
{
   DeviceMgr * DM_ptr = new DeviceMgr(cp);
   DM_ptr->DevMgrPtr_ = DM_ptr;

   return DM_ptr;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::DeviceMgr
// Purpose       : constructor
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
DeviceMgr::DeviceMgr(N_IO_CmdParse & cp)
  : commandLine_(cp),
    devOptions_(cp),
    icLoads_(NULL),
    devSensPtr_ (NULL),
    jacobianLoadCalledBefore_(false),
    entityMapDone_ (false),
    numPDEDevices_(0),
    calledBeforeCSPI (false),
    numSensParams_(0),
    allDevsConverged_(false),
    sensFlag_(false),
    linearSystemFlag_(true),
    firstDependent(true),
    externalStateFlag_(false),
    parameterChanged_(false),
    breakPointInstancesInitialized(false),
    timeParamsProcessed_(0.0),
    devBuilder_(solState_, devOptions_),
    numThreads_(0),
    multiThreading_(false),
    nonTrivialDeviceMaskFlag(false),
    dotOpOutputFlag(false),
    numJacStaVectorPtr_(0),
    numJacSolVectorPtr_(0),
    numJacStoVectorPtr_(0),
    diagonalVectorPtr_(0)
{

  int i;
  dummy_ptr = PlaceHolder::factory(solState_,devOptions_);

  // instantiate the deviceAllocFlag array.
  for (i = 0; i < ModelType::NUMDEV; ++i)
  {
    deviceAllocFlag_[i] = 0;
    deviceUseFlag_[i]   = 0;
    deviceArray_[i]     = dummy_ptr;
    PDEDeviceFlag_[i]   = 0;
  }


  setUpDeviceIndexMap_();
  setUpMOSFETLevelMap_();
  setUpJFETLevelMap_();
  setUpDiodeLevelMap_();
  setUpBJTLevelMap_();
  setUpPDELevelMap_();
  setUpResistorLevelMap_();
  setUpIndLevelMap_();
  setUpNeuronLevelMap_();
  setUpNeuronPopLevelMap_();
  setUpSynapseLevelMap_();
  setUpRadLevelMap_();
  setUpDeviceModelTypeMap_();
  setUpPDEDeviceFlagArray_ ();

#ifdef Xyce_EXTDEV
  setUpPassThroughParamsMap_();
#endif

  // set the solution-device map pointer:   does not exist anymore...
  extData_.solDevInstMap = & solDevInstMap_;

  extData_.devMgrPtr = this;

#ifdef Xyce_DEBUG_DEVICE
  solState_.debugTimeFlag = true;
#endif
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::~DeviceMgr
// Purpose       : destructor
// Special Notes : De-allocates all the devices pointed  to by deviceArray
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
DeviceMgr::~DeviceMgr()
{
  delete numJacStaVectorPtr_;
  delete numJacSolVectorPtr_;
  delete numJacStoVectorPtr_;
  delete diagonalVectorPtr_;

  delete extData_.numJacRHSVectorPtr;
  delete extData_.numJacFVectorPtr;
  delete extData_.numJacQVectorPtr;
  delete extData_.perturbVectorPtr;
  delete extData_.numJacLoadFlagPtr;

  delete extData_.tmpdIdXPtr;
  delete extData_.tmpdQdXPtr;

  for (int i = 0; i < ModelType::NUMDEV; ++i)
  {
    if (deviceAllocFlag_[i] == 1)
    {
      delete deviceArray_[i];
      deviceAllocFlag_[i] = 0;
      deviceUseFlag_[i]   = 0;
    }
    PDEDeviceFlag_[i]   = 0;
  }

#ifdef Xyce_EXTDEV
  for (int i = 0; i < extDevIBPtrVec_.size(); ++i)
  {
    if (extDevIBPtrVec_[i] != NULL)
    {
      delete extDevIBPtrVec_[i];
      extDevIBPtrVec_[i] = 0;
    }
  }
#endif

  delete icLoads_;
  delete devSensPtr_;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::createDeviceByNetlistDeviceType
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Baur, SNL
// Creation Date : 04/18/2013
//-----------------------------------------------------------------------------
Device * DeviceMgr::createDeviceByNetlistDeviceType(const std::string &name, const int level)
{
  return devBuilder_.createDeviceByNetlistDeviceType(name, level);
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::registerPkgOptionsMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Mems Modeling
// Creation Date : 10/20/2008
//-----------------------------------------------------------------------------
bool DeviceMgr::registerPkgOptionsMgr( RCP<N_IO_PkgOptionsMgr> pkgOptPtr )
{
  pkgOptMgrPtr_ = pkgOptPtr;
  string netListFile("");
  if (commandLine_.getArgumentValue(string("netlist")) != "")
  {
    netListFile = commandLine_.getArgumentValue(string("netlist"));
  }
  pkgOptMgrPtr_->submitRegistration(
      "DEVICE", netListFile, new DeviceMgr_OptionsReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "SENS", netListFile, new DeviceMgr_SensOptionsReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "TIMEINT", netListFile, new DeviceMgr_TimeOptionsReg( this ) );

  // different analysis types.
  pkgOptMgrPtr_->submitRegistration(
      "TRAN", netListFile, new DeviceMgr_TransAnalysisReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "DC", netListFile, new DeviceMgr_DCAnalysisReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "OP", netListFile, new DeviceMgr_OPAnalysisReg( this ) );

  pkgOptMgrPtr_->submitRegistration(
      "STEP", netListFile, new DeviceMgr_STEPAnalysisReg( this ) );

  // MPDE specific netlist lines
  pkgOptMgrPtr_->submitRegistration(
      "MPDE", netListFile, new DeviceMgr_MPDE_AnalysisReg( this ) );

  // HB Specific netlist lines
  pkgOptMgrPtr_->submitRegistration(
      "HB", netListFile, new DeviceMgr_HB_AnalysisReg( this ) );

  // AC Specific netlist lines
  pkgOptMgrPtr_->submitRegistration(
      "AC", netListFile, new DeviceMgr_AC_AnalysisReg( this ) );

  // MOR Specific netlist lines
  pkgOptMgrPtr_->submitRegistration(
      "MOR", netListFile, new DeviceMgr_MOR_AnalysisReg( this ) );

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::registerSensParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/03/02
//-----------------------------------------------------------------------------
bool DeviceMgr::registerSensParams (const N_UTL_OptionBlock & OB)
{
  sensFlag_ = true;

#ifdef Xyce_DEBUG_DEVICE
  if (devOptions_.debugLevel > 0)
  {
    std::cout << "DeviceMgr::registerSensParams called!" <<std::endl;
  }
#endif

  // the devSensPtr will be deleted in the destructor.
  if (devSensPtr_==0)
  {
    devSensPtr_ = new DeviceSensitivities
      (devOptions_,extData_, solState_, *lasSysPtr_);
  }

  return devSensPtr_->registerSensParams (OB);
}


//-----------------------------------------------------------------------------
// Function      : N_DEV_DeviceMgr::registerLeadCurrentRequests
// Purpose       : this function is called from the output manager (through the
//                 device interface) to inform the device package of the devices
//                 for which lead currents have been requested.  The device
//                 manager will take care of doing isolated F and Q loads for
//                 these devices so the lead currents can be calculated
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical Systems Modeling
// Creation Date : 03/20/13
//-----------------------------------------------------------------------------
bool N_DEV_DeviceMgr::setLeadCurrentRequests(const std::set<std::string> & deviceNames )
{
  // this is called prior to fully constructing the devices.  So for now
  // save the list
  devicesNeedingLeadCurrentLoads_ = deviceNames;

#ifdef Xyce_DEBUG_DEVICE
  if (devOptions_.debugLevel > 0)
  {
    if( ! devicesNeedingLeadCurrentLoads_.empty() )
    {
      set<string>::iterator currentDeviceNameItr = devicesNeedingLeadCurrentLoads_.begin();
      set<string>::iterator endDeviceNameItr = devicesNeedingLeadCurrentLoads_.end();
      std::cout << "N_DEV_DeviceMgr::registerLeadCurrentRequests Devices for which lead currents were requested: ";
      while ( currentDeviceNameItr != endDeviceNameItr )
      {
        std::cout << *currentDeviceNameItr << "  ";
        currentDeviceNameItr++;
      }
      std::cout << std::endl;
    }
  }
#endif
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setTranAnalysisParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/14/11
//-----------------------------------------------------------------------------
bool DeviceMgr::setTranAnalysisParams (const N_UTL_OptionBlock & OB)
{
  solState_.TRANspecified = true;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setDCAnalysisParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/14/11
//-----------------------------------------------------------------------------
bool DeviceMgr::setDCAnalysisParams (const N_UTL_OptionBlock & OB)
{
  solState_.DCspecified = true;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setOPAnalysisParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/14/11
//-----------------------------------------------------------------------------
bool DeviceMgr::setOPAnalysisParams (const N_UTL_OptionBlock & OB)
{
  solState_.OPspecified = true;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setSTEPAnalysisParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/14/11
//-----------------------------------------------------------------------------
bool DeviceMgr::setSTEPAnalysisParams (const N_UTL_OptionBlock & OB)
{
  solState_.STEPspecified = true;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setMPDEAnalysisParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/14/11
//-----------------------------------------------------------------------------
bool DeviceMgr::setMPDEAnalysisParams (const N_UTL_OptionBlock & OB)
{
  solState_.MPDEspecified = true;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setHBAnalysisParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/14/11
//-----------------------------------------------------------------------------
bool DeviceMgr::setHBAnalysisParams (const N_UTL_OptionBlock & OB)
{
  solState_.HBspecified = true;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setACAnalysisParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/14/11
//-----------------------------------------------------------------------------
bool DeviceMgr::setACAnalysisParams (const N_UTL_OptionBlock & OB)
{
  solState_.ACspecified = true;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setMORAnalysisParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist and Ting Mei, SNL
// Creation Date : 5/30/12
//-----------------------------------------------------------------------------
bool DeviceMgr::setMORAnalysisParams (const N_UTL_OptionBlock & OB)
{
  solState_.MORspecified = true;
  return true;
}


//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getFastSourcePeriod
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/27/04
//-----------------------------------------------------------------------------
vector<double> DeviceMgr::getFastSourcePeriod (vector<string>& sourceNames)
{
  int numFastSrcs = sourceNames.size();

  //Setup return of source periods
  vector<double> srcPeriods(numFastSrcs);

  // Now loop over them, and mark them.
  for( int i = 0; i < numFastSrcs; ++i )
  {
    ExtendedString tmpName(sourceNames[i]);
    tmpName.toUpper();
    map <string,SourceInstance*>::iterator iterFS = indepSourcePtrMap_.find(tmpName);
    if ( iterFS != indepSourcePtrMap_.end() )
    {
      SourceInstance * SIPtr = iterFS->second;
      srcPeriods[i] = SIPtr->period();
    }
    else
    {
      string msg("DeviceMgr::getFastSourcePeriod ");
      msg += "Unable to find source: " + fastSourceNames_[i] + "\n";
      msg += "Potential names are: ";
      map <string,SourceInstance*>::iterator currentFS = indepSourcePtrMap_.begin();
      map <string,SourceInstance*>::iterator endFS = indepSourcePtrMap_.end();
      while( currentFS != endFS )
      {
        msg += (*currentFS).first;
        msg += " ";
        currentFS++;
      }

#ifdef Xyce_PARALLEL_MPI
#else
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL, msg);
#endif
    }
  }

  return srcPeriods;
}


//-----------------------------------------------------------------------------
// Function      : DeviceMgr::registerFastSources
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/27/04
//-----------------------------------------------------------------------------
vector<double> DeviceMgr::registerFastSources (vector<string>& sourceNames)
{
  int numFastSrcs = sourceNames.size();
  //Setup return of source periods
  vector<double> srcPeriods;

  if (numFastSrcs > 0)
  {
    srcPeriods.resize(numFastSrcs);

      // Default case, the sources are explicitely listed on the .option line
    fastSourceNames_.resize(numFastSrcs);
    copy( sourceNames.begin(), sourceNames.end(), fastSourceNames_.begin());

    // Now loop over them, and mark them.
    for( int i = 0; i < numFastSrcs; ++i )
    {
      ExtendedString tmpName(fastSourceNames_[i]);
      tmpName.toUpper();
      map <string,SourceInstance*>::iterator iterFS = indepSourcePtrMap_.find(tmpName);
      if ( iterFS != indepSourcePtrMap_.end() )
      {
        SourceInstance * SIPtr = iterFS->second;
        SIPtr->setFastSourceFlag (true);
        srcPeriods[i] = SIPtr->period();
      }
      else
      {
        string msg("DeviceMgr::registerFastSources ");
        msg += "Unable to find source: " + fastSourceNames_[i] + "\n";
        msg += "Potential names are: ";
        map <string,SourceInstance*>::iterator currentFS = indepSourcePtrMap_.begin();
        map <string,SourceInstance*>::iterator endFS = indepSourcePtrMap_.end();
        while( currentFS != endFS )
        {
          msg += (*currentFS).first;
          msg += " ";
          currentFS++;
        }

#ifdef Xyce_PARALLEL_MPI
#else
        N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL, msg);
#endif
      }
    }

  }
  else
  {
    // tscoffe/tmei 09/16/08
    // Special case:  Use all sources
    numFastSrcs = indepSourceInstancePtrVec_.size();
    srcPeriods.resize(numFastSrcs);
    for (int i=0 ; i<numFastSrcs ; ++i) {
      indepSourceInstancePtrVec_[i]->setFastSourceFlag(true);
      srcPeriods[i] = indepSourceInstancePtrVec_[i]->period();
    }
  }
  return srcPeriods;
}


//-----------------------------------------------------------------------------
// Function      : DeviceMgr::deRegisterFastSources
// Purpose       : reverses the effect of registerFastSources
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 06/12/2013
//-----------------------------------------------------------------------------
void DeviceMgr::deRegisterFastSources (vector<string>& sourceNames)
{
  int numFastSrcs = sourceNames.size();

  if (numFastSrcs > 0)
  {
      // Default case, the sources are explicitely listed on the .option line
    fastSourceNames_.resize(numFastSrcs);
    copy( sourceNames.begin(), sourceNames.end(), fastSourceNames_.begin());

    // Now loop over them, and mark them.
    for( int i = 0; i < numFastSrcs; ++i )
    {
      ExtendedString tmpName(fastSourceNames_[i]);
      tmpName.toUpper();
      map <string,SourceInstance*>::iterator iterFS = indepSourcePtrMap_.find(tmpName);
      if ( iterFS != indepSourcePtrMap_.end() )
      {
        SourceInstance * SIPtr = iterFS->second;
        SIPtr->setFastSourceFlag (false);
      }
      else
      {
        string msg("DeviceMgr::registerFastSources ");
        msg += "Unable to find source: " + fastSourceNames_[i] + "\n";
        msg += "Potential names are: ";
        map <string,SourceInstance*>::iterator currentFS = indepSourcePtrMap_.begin();
        map <string,SourceInstance*>::iterator endFS = indepSourcePtrMap_.end();
        while( currentFS != endFS )
        {
          msg += (*currentFS).first;
          msg += " ";
          currentFS++;
        }

#ifdef Xyce_PARALLEL_MPI
#else
        N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL, msg);
#endif
      }
    }

  }
  else
  {
    // Special case:  Use all sources
    numFastSrcs = indepSourceInstancePtrVec_.size();
    for (int i=0 ; i<numFastSrcs ; ++i) {
      indepSourceInstancePtrVec_[i]->setFastSourceFlag(false);
    }
  }
  return;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setFastSourceSetFlag
// Purpose       : traverse fast source list and remove any slow sources from
//                 the deviceArray
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 3/22/07
//-----------------------------------------------------------------------------
void DeviceMgr::deactivateSlowSources()
{
  // first back-up a copy of the deviceArray so we can edit out
  // the slow sources
  indepSourceInstanceBackupPtrVec_.resize( indepSourceInstancePtrVec_.size() );
  copy( indepSourceInstancePtrVec_.begin(), indepSourceInstancePtrVec_.end(),
    indepSourceInstanceBackupPtrVec_.begin());

  // erase the existing list of sources
  indepSourceInstancePtrVec_.clear();

  // now copy back only those that are fast sources
  vector<SourceInstance*>::iterator iter;
  vector<SourceInstance*>::iterator begin =indepSourceInstanceBackupPtrVec_.begin();
  vector<SourceInstance*>::iterator end =indepSourceInstanceBackupPtrVec_.end();
  for (iter=begin; iter!=end;++iter)
  {
    if( (*iter)->getFastSourceFlag() )
    {
      indepSourceInstancePtrVec_.push_back( *iter );
    }
  }

}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setFastSourceSetFlag
// Purpose       : restore any slow sources to the device array.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 3/22/07
//-----------------------------------------------------------------------------
void DeviceMgr::activateSlowSources()
{
  // restore the independent source list from backup
  indepSourceInstancePtrVec_.clear();
  indepSourceInstancePtrVec_.resize( indepSourceInstanceBackupPtrVec_.size() );
  copy( indepSourceInstanceBackupPtrVec_.begin(), indepSourceInstanceBackupPtrVec_.end(),
    indepSourceInstancePtrVec_.begin());
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setMPDEFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Parallel Computational Sciences
// Creation Date : 07/21/08
//-----------------------------------------------------------------------------
void DeviceMgr::setMPDEFlag( bool flagVal )
{
  solState_.mpdeOnFlag  = flagVal;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setBlockAnalysisFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Parallel Computational Sciences
// Creation Date : 07/21/08
//-----------------------------------------------------------------------------
void DeviceMgr::setBlockAnalysisFlag( bool flagVal )
{
  solState_.blockAnalysisFlag = flagVal;
  devOptions_.setBlockAnalysisFlag(flagVal);
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setFastTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Parallel Computational Sciences
// Creation Date : 07/21/08
//-----------------------------------------------------------------------------
void DeviceMgr::setFastTime( double timeVal )
{
  solState_.currFastTime = timeVal;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::initializeAll
// Purpose       : This function, via the LAS system class, sets up
//                 the pointers to the various linear algebra entities.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/12/00
//-----------------------------------------------------------------------------
bool DeviceMgr::initializeAll()
{
  bool bsuccess = true;

  extData_.lasSysPtr = lasSysPtr_;

  // nullify ptrs that are passed in at each step  (see the loadDAEVectors function args)
  extData_.nextSolVectorPtr = 0;
  extData_.currSolVectorPtr = 0;
  extData_.lastSolVectorPtr = 0;
  extData_.daeQVectorPtr    = 0;
  extData_.daeFVectorPtr    = 0;
  extData_.dFdxdVpVectorPtr = 0;
  extData_.dQdxdVpVectorPtr = 0;
  extData_.nextStaVectorPtr = 0;
  extData_.currStaVectorPtr = 0;
  extData_.lastStaVectorPtr = 0;
  extData_.storeLeadCurrQCompPtr = 0;
  extData_.nextStaDerivVectorPtr = 0;
  extData_.nextStoVectorPtr =
  extData_.currStoVectorPtr = 0;
  extData_.lastStoVectorPtr = 0;

#ifdef Xyce_DEBUG_DEVICE
  // get f vector pointer:
  extData_.fVectorPtr  = lasSysPtr_->getFVector();
  bsuccess = bsuccess && (extData_.fVectorPtr != 0);

  // create Jdxp vector pointer:
  extData_.JdxpVectorPtr = lasSysPtr_->getJDXPVector();
  bsuccess = bsuccess && (extData_.JdxpVectorPtr != 0);
#endif

#ifdef Xyce_DEBUG_VOLTLIM
  // get test matrix: (old DAE)
  extData_.JTestMatrixPtr = lasSysPtr_->getJacTestMatrix();
  bsuccess = bsuccess && (extData_.JTestMatrixPtr != 0);

  // get test matrix:
  extData_.FTestMatrixPtr = lasSysPtr_->getdFdxTestMatrix();
  bsuccess = bsuccess && (extData_.FTestMatrixPtr != 0);
  extData_.QTestMatrixPtr = lasSysPtr_->getdQdxTestMatrix();
  bsuccess = bsuccess && (extData_.QTestMatrixPtr != 0);

  // get dxVoltlim vector pointer:
  extData_.dxVoltlimVectorPtr = lasSysPtr_->getDxVoltlimVector() ;
  bsuccess = bsuccess && (extData_.dxVoltlimVectorPtr != 0);

  // create Jdx2 vector pointer: (old DAE)
  extData_.Jdx2VectorPtr = lasSysPtr_->getJDX2Vector();
  bsuccess = bsuccess && (extData_.Jdx2VectorPtr != 0);

  // create Jdx2 vector pointer:
  extData_.Fdx2VectorPtr = lasSysPtr_->getFDX2Vector();
  bsuccess = bsuccess && (extData_.Fdx2VectorPtr != 0);
  extData_.Qdx2VectorPtr = lasSysPtr_->getQDX2Vector();
  bsuccess = bsuccess && (extData_.Qdx2VectorPtr != 0);
#endif

  // get flag solution pointer-pointer:
  extData_.flagSolVectorPtr = lasSysPtr_->getFlagSolVector();
  bsuccess = bsuccess && (extData_.flagSolVectorPtr != 0);

  // get device mask pointer.
  extData_.deviceMaskVectorPtr = lasSysPtr_->getDeviceMaskVector ();
  bsuccess = bsuccess && (extData_.deviceMaskVectorPtr != 0);

  // create the temporary numerical jacobian vectors
  if (devOptions_.numericalJacobianFlag || devOptions_.testJacobianFlag || sensFlag_)
  {
    numJacStaVectorPtr_ = lasSysPtr_->builder().createStateVector();
    numJacSolVectorPtr_ = lasSysPtr_->builder().createVector();
    numJacStoVectorPtr_ = lasSysPtr_->builder().createStoreVector();

    extData_.numJacRHSVectorPtr = lasSysPtr_->builder().createVector();
    extData_.numJacFVectorPtr   = lasSysPtr_->builder().createVector();
    extData_.numJacQVectorPtr   = lasSysPtr_->builder().createVector();
    extData_.perturbVectorPtr   = lasSysPtr_->builder().createVector();
    extData_.numJacLoadFlagPtr  = lasSysPtr_->builder().createVector();
  }

  extData_.tmpdIdXPtr = lasSysPtr_->builder().createVector();
  extData_.tmpdQdXPtr = lasSysPtr_->builder().createVector();

  // create a diagonal vector to be used for 2-level
  diagonalVectorPtr_  = lasSysPtr_->builder().createVector();

  extData_.initializeAllFlag = true;

  // For Homotopy on block gainscale
  solState_.InitializeHomotopyBlockSize(devOptions_.numGainScaleBlocks);

#ifdef Xyce_SIZEOF
  int size = sizeof(*this);
  std::cout << "Size of device package after initializeAll  = " << size << std::endl;
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::resetForStepAnalysis
// Purpose       : 
// Special Notes : Some "resetForStep" functions (only HB so far) will 
//                 call dev->initializeAll.  So, this must be called first.
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/12/2013
//-----------------------------------------------------------------------------
void DeviceMgr::resetForStepAnalysis()
{
  delete numJacStaVectorPtr_;
  delete numJacSolVectorPtr_;
  delete numJacStoVectorPtr_;
  delete extData_.numJacRHSVectorPtr;
  delete extData_.numJacFVectorPtr;
  delete extData_.numJacQVectorPtr;
  delete extData_.perturbVectorPtr;
  delete extData_.numJacLoadFlagPtr;
  delete extData_.tmpdIdXPtr;
  delete extData_.tmpdQdXPtr;
  delete diagonalVectorPtr_;

  solState_.ltraTimeIndex = 0;
  solState_.ltraTimeHistorySize = 10;
  solState_.ltraTimePoints.resize(solState_.ltraTimeHistorySize);
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::createDevice
// Purpose       : This function creates a single device based on the passed
//                 index.
// Special Notes :
// Scope         : protected
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
bool DeviceMgr::createDeviceByModelType(const int model_type)
{
  bool itmp = true;

  if (model_type != 0)  // do not re-allocate the dummy.
  {
    if (deviceArray_[model_type] == NULL || deviceArray_[model_type] == dummy_ptr || deviceAllocFlag_[model_type] == 0)
    {
      deviceArray_[model_type] = devBuilder_.createDeviceByModelType(model_type);
      deviceAllocFlag_[model_type] = 1;

      devicePtrVec_.push_back(deviceArray_[model_type]);

      if (PDEDeviceFlag_[model_type] == 1)
      {
        pdeDevicePtrVec_.push_back(deviceArray_[model_type]);
      }
      else
      {
        nonPdeDevicePtrVec_.push_back(deviceArray_[model_type]);
      }

      itmp = true;
    }
  }

  return itmp;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::flushDevices
// Purpose       : This function deletes all the allocated devices, and resets
//                 the array of device pointers to all point to the dummy
//                 (placeholder) device class pointer.
// Special Notes :
// Scope         : protected
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
bool DeviceMgr::flushDevices ()
{
  int i;
  bool itmp = false;

  for (i = 0; i < ModelType::NUMDEV; ++i)
  {
    if (deviceAllocFlag_[i] == 1)
    {
      delete deviceArray_[i];
      deviceArray_    [i] = dummy_ptr;
      deviceAllocFlag_[i] = 0;
    }
  }

  itmp = true;
  return itmp;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getDeviceIndex
// Purpose       : This function returns the device type index for a given
//                 named string.  This assumes that the device names used
//                 in the simulation will obey the spice3f5 netlist language
//                 convention.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
int DeviceMgr::getDeviceIndex (const string & Name)
{
  int type = 0;

  DeviceIndexMap::const_iterator iterDI;
  iterDI = deviceIndexMap_.find(Name);
  if (iterDI != deviceIndexMap_.end())
  {
    return iterDI->second;
  }

  // first get the first letter of the name string.
  string letter(Name.begin(),Name.begin()+1);

  // Use the letter as an index into the deviceIndex map.
  if ( letter != "Y" )
  {
    type = deviceIndexMap_[letter];
  }
  else
  {
    // "Y" device requires special handling. The device type was embedded in
    // the device name to fit the way things worked for ordinary devices.
    // The format of the device name for a "Y" device is
    // "Y%<DeviceType>%<NetlistDeviceName>". Extract the device type and
    // use as the index into the deviceIndexMap map.
    string deviceType(Name.substr(Name.find_first_of("%")+1, Name.find_last_of("%")-2));
    type = deviceIndexMap_[deviceType];
  }

  return type;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getModelTypeIndex
// Purpose       : This function returns the device type index for a given
//                 named string.  The named string corrsponds to a model
//                 "type" name.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/06/00
//-----------------------------------------------------------------------------
int DeviceMgr::getModelTypeIndex(const string & ModelType)
{
  int type = 0;

  // Obtain the device index using the deviceModelTypeMap.
  type = deviceModelTypeMap_[ModelType];

  return type;
}

// //-----------------------------------------------------------------------------
// // Function      : DeviceMgr::getDeviceTypeOffset
// // Purpose       :
// // Special Notes :
// // Scope         : public
// // Creator       : Eric Keiter, SNL
// // Creation Date : 10/26/11
// //-----------------------------------------------------------------------------
// int DeviceMgr::getDeviceTypeOffset(const string & name, const int level, const string & modelType)
// {
//   ModelBlock MB;
//   MB.name = name;
//   MB.type = modelType;
//   MB.level = level;
//   return getDeviceTypeOffset(MB);
// }

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getDeviceTypeOffset
// Purpose       : Some devices (mainly MOSFETs) will have a different
//                 device array index, depending on what level is specified
//                 by their model statement.
//
//                 This version of the function determines the offset based
//                 on the model type and level.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
int DeviceMgr::getDeviceTypeOffset (const ModelBlock & MB)
{
  int offset = 0;

  ostringstream ost;

  // is this a mosfet?  if so the offset is not zero.
  ExtendedString tmpType(MB.type);
  tmpType.toUpper ();

  map<int,int>::iterator iterLevel;

  // offset for MOSFET levels:
  if (tmpType == "PMOS" || tmpType == "NMOS" || tmpType == "M")
  {
    iterLevel = MOSFETLevelMap_.find(MB.level);
    if ( iterLevel == MOSFETLevelMap_.end() )
    {
      ost <<"This Xyce build doesn't include";
      ost << " the MOSFET level=" << MB.level << " used by " << MB.name ;
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, ost.str() );
    }
    offset = iterLevel->second;
  }
  else if (tmpType == "PJF" || tmpType == "NJF" || tmpType == "J")
  {
    iterLevel = JFETLevelMap_.find(MB.level);
    if ( iterLevel == JFETLevelMap_.end() )
    {
      ost <<"This Xyce build doesn't include";
      ost << " the JFET level=" << MB.level << " used by " << MB.name ;
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, ost.str() );
    }
    offset = iterLevel->second;
  }
  else if (tmpType == "D")
  {
    iterLevel = DiodeLevelMap_.find(MB.level);
    if ( iterLevel == DiodeLevelMap_.end() )
    {
      ost <<"This Xyce build doesn't include";
      ost << " the DIODE level=" << MB.level << " used by " << MB.name ;
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, ost.str() );
    }
    offset = iterLevel->second;
  }
  else if (tmpType == "PNP" || tmpType == "NPN" || tmpType == "Q" || tmpType == "VBIC")
  {
    iterLevel = BJTLevelMap_.find(MB.level);
    if ( iterLevel == BJTLevelMap_.end() )
    {
      ost <<"This Xyce build doesn't include";
      ost << " the BJT level=" << MB.level << " used by " << MB.name ;
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, ost.str() );
    }
    offset = iterLevel->second;
  }
  // offset for PDE device levels.
  else if (tmpType == "ZOD" || tmpType == "Y%PDE%" || tmpType == "PDE")
  {
    iterLevel = PDELevelMap_.find(MB.level);
    if ( iterLevel == PDELevelMap_.end() )
    {
      ost <<"This Xyce build doesn't include";
      ost << " the PDE level=" << MB.level << " used by " << MB.name ;
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, ost.str() );
    }
    offset = iterLevel->second;
  }
  else if (tmpType == "R")
  {
    iterLevel = ResistorLevelMap_.find(MB.level);
    if ( iterLevel == ResistorLevelMap_.end() )
    {
      ost <<"This Xyce build doesn't include";
      ost << " the Resistor level=" << MB.level << " used by " << MB.name ;
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, ost.str() );
    }
    offset = iterLevel->second;
  }
  else if (tmpType == "CORE")
  {
    iterLevel = INDLevelMap_.find(MB.level);
    if ( iterLevel == INDLevelMap_.end() )
    {
      ost <<"This Xyce build doesn't include";
      ost << " the Mutual Inductor level=" << MB.level << " used by " << MB.name ;
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, ost.str() );
    }
    offset = iterLevel->second;
  }
  else if (tmpType == "NEURON")
  {
    iterLevel = NeuronLevelMap_.find(MB.level);
    if ( iterLevel == NeuronLevelMap_.end() )
    {
      ost <<"This Xyce build doesn't include";
      ost << " the Neuron level=" << MB.level << " used by " << MB.name ;
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, ost.str() );
    }
    offset = iterLevel->second;
  }
  else if (tmpType == "SYNAPSE")
  {
    iterLevel = SynapseLevelMap_.find(MB.level);
    if ( iterLevel == SynapseLevelMap_.end() )
    {
      ost <<"This Xyce build doesn't include";
      ost << " the Synapse level=" << MB.level << " used by " << MB.name ;
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, ost.str() );
    }
    offset = iterLevel->second;
  }
  else if (tmpType == "NEURONPOP")
  {
    iterLevel = NeuronPopLevelMap_.find(MB.level);
    if ( iterLevel == NeuronPopLevelMap_.end() )
    {
      ost <<"This Xyce build doesn't include";
      ost << " the NeuronPop level=" << MB.level << " used by " << MB.name ;
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, ost.str() );
    }
    offset = iterLevel->second;
  }
  else if (tmpType=="RAD" || tmpType == "NEUTRON")
  {
    iterLevel = radLevelMap_.find(MB.level);
    if ( iterLevel == radLevelMap_.end() )
    {
      ost <<"This Xyce build doesn't include";
      ost << " the rad level=" << MB.level << " used by " << MB.name ;
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, ost.str() );
    }
    offset = iterLevel->second;
  }

  return offset;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::addDeviceModel
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
bool DeviceMgr::addDeviceModel (const ModelBlock & MB)
{
  int baseType = getModelTypeIndex(MB.type);
  int offset = 0;

  if (MB.name != "") offset = getDeviceTypeOffset(MB);

  int type = baseType + offset;

  // Just in case, call the device creator function - if it has already
  // been allocated it won't create a redundant one.
  createDeviceByModelType(type);

  DeviceModel * dmPtr = deviceArray_[type]->addModel(MB);

  deviceUseFlag_[type] = 1;

  deviceModelNameMap_[MB.name] = type;
  deviceModelNameBaseTypeMap_[MB.name] = baseType;

  // add the various model vectors:
  if (dmPtr != 0)
  {
    modelPtrVec_.push_back(dmPtr);

    if ( baseType == ModelType::MOSFET1) mosfetModelPtrVec_.push_back(dmPtr);
    if ( baseType == ModelType::BJT    ) bjtModelPtrVec_.push_back(dmPtr);

    if ( type == ModelType::DIODE ) diodeModelPtrVec_.push_back(dmPtr);
    if ( type == ModelType::MOSFET_B3 ) bsim3ModelPtrVec_.push_back(dmPtr);
    if ( type == ModelType::MOSFET_B4 ) bsim4ModelPtrVec_.push_back(dmPtr);
    if ( type == ModelType::MOSFET_B3SOI ) bsimsoiModelPtrVec_.push_back(dmPtr);
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::verifyDeviceInstance
// Purpose       : This function verifies a device instance prior to
//                 instantiating.
//
//                 Theoretically, we could do this in addDeviceInstance() and
//                 not make a device if it fails some verification criteria (a
//                 resistor with zero resistance is the primary case here).
//                 However, later unlinking of redundant devices that are
//                 connected to just one node is difficult after
//                 addDeviceInstance() is called because the device instance
//                 pointer can be placed in many other containers
//                 It could be done, but this is a simpler first step to having the
//                 device manager be in charge of device verification -- rather
//                 than have it in toplogy or IO.
//
// Special Notes : return true if this device is ok to instantiate, false otherwise
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/18/2010
//-----------------------------------------------------------------------------
bool DeviceMgr::verifyDeviceInstance(InstanceBlock & IB)
{
  int  type=0;

  if (IB.getModelName() == "")
  {
    if( IB.getName().find_last_of(":") != string::npos )
    {
      type = getDeviceIndex(
        IB.getName().substr(IB.getName().find_last_of(":")+1,
              string::npos) );
    }
    else
    {
      type = getDeviceIndex( IB.getName() );
    }
  }
  else
  {
    type = deviceModelNameMap_[IB.getModelName()];
  }

  if (type == 0)
  {
     string msg("Device Manager::verifyDeviceInstance() could not correctly figure out what");
     msg += " type of device " + IB.getName() + " is.\n";
     if (IB.getModelName() != "")
     {
       msg += "Its model name is " + IB.getModelName() +", which was not found.\n";
     }
     else
     {
       msg += "Could not find this device name in the device index.\n";
     }
     N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
  }
  else if (IB.bsourceFlag)
  {
    // This is an E, F, G, or H source that is to be treated as a B source,
    // set its type now to BSRC.
    type = ModelType::BSRC;
  }

  // Check if this is a simple resistor, but with a resistance of zero. if so,
  // then return false so we don't make this device.
  if( (devOptions_.checkForZeroResistance) && (type == ModelType::RESISTOR ) )
  {
    const double zeroResistanceValue = devOptions_.zeroResistanceTol;
    // loop over the parameters
    vector<Param>::iterator currentParam = IB.params.begin();
    vector<Param>::iterator endParam = IB.params.end();
    while( currentParam != endParam )
    {
#ifdef Xyce_DEBUG_DEVICE
      //std::cout << "In addInstance:  Resistor parameter = " << currentParam->uTag()  << std::endl;
#endif
      if( (currentParam->uTag() == "R") )
      {
        if (currentParam->given())
        {
          vector<string> variables, specials;

          // check if this is a time-dependent, or variable-dependent expression.
          // If it is, then skip.
          Param * devPar = &(*(currentParam));
          N_UTL_Param * tmpPar = (dynamic_cast<N_UTL_Param*> (devPar));
          // only check if this is an expression-type parameter.
          if (tmpPar->getType() == EXPR)
          {
            N_UTL_Expression tmpExp = tmpPar->eVal();
            tmpExp.get_names(XEXP_VARIABLE, variables);
            tmpExp.get_names(XEXP_SPECIAL, specials);
          }

          if (specials.empty() && variables.empty())
          {
            if (fabs(currentParam->dVal()) < devOptions_.zeroResistanceTol)
            {
              // change device type to be the level 3 resistor
              return false;
            }
          }
        }
        break;
      }
      currentParam++;
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::addDeviceInstance
// Purpose       : addDeviceInstance will create a new instance of the
//                 designated device type.  This version of the function
//                 accepts a parameter list as one of the arguments,
//                 so it is assumed that a parameter instance will
//                 also have to be allocated for it.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
DeviceInstance * DeviceMgr::addDeviceInstance(InstanceBlock & IB)
{
  DeviceInstance * DI_ptr;

  int  type;
  int  baseType;

  // If the ModelName string is not a null string, use it to get the
  // device type index.
  // Otherwise, use the name of the device itself to determine the device.

  //if (IB.getModelName() == "") type = getDeviceIndex(IB.name);
  //else                    type = deviceModelNameMap_[IB.getModelName()];

  if (IB.getModelName() == "")
  {
    if( IB.getName().find_last_of(":") != string::npos )
    {
      type = getDeviceIndex(IB.getName().substr(IB.getName().find_last_of(":")+1, string::npos) );
    }
    else
    {
      type = getDeviceIndex( IB.getName() );
    }
    baseType = type;
  }
  else
  {
    type = deviceModelNameMap_[IB.getModelName()];
    baseType = deviceModelNameBaseTypeMap_[IB.getModelName()];
  }

  if (type == 0)
  {
     string msg("Device Manager could not correctly figure out what");
     msg += " type of device " + IB.getName() + " is.\n";
     if (IB.getModelName() != "")
     {
       msg += "Its model name is " + IB.getModelName() +", which was not found.\n";
     }
     else
     {
       msg += "Could not find this device name in the device index.\n";
     }
     N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
  }
  else if (IB.bsourceFlag)
  {
    // This is an E, F, G, or H source that is to be treated as a B source,
    // set its type now to BSRC.
    type = ModelType::BSRC;
  }


  // Check if this is a simple resistor, but with a resistance of zero.  If so,
  // then change the type to RESISTOR3.  This variant of the resistor acts like
  // a voltage souce with zero voltage difference.
  if( (devOptions_.checkForZeroResistance) && (type == ModelType::RESISTOR ) )
  {
    const double zeroResistanceValue = devOptions_.zeroResistanceTol;
    // loop over the parameters
    vector<Param>::iterator currentParam = IB.params.begin();
    vector<Param>::iterator endParam = IB.params.end();
    while( currentParam != endParam )
    {
      if( (currentParam->uTag() == "R") )
      {
        if (currentParam->given())
        {
          vector<string> variables, specials;

          // check if this is a time-dependent, or variable-dependent expression.
          // If it is, then skip.
          Param * devPar = &(*(currentParam));
          N_UTL_Param * tmpPar = (dynamic_cast<N_UTL_Param*> (devPar));
          // only check if this is an expression-type parameter.
          if (tmpPar->getType() == EXPR)
          {
            N_UTL_Expression tmpExp = tmpPar->eVal();
            tmpExp.get_names(XEXP_VARIABLE, variables);
            tmpExp.get_names(XEXP_SPECIAL, specials);
          }

          if (specials.empty() && variables.empty())
          {
            if (fabs(currentParam->dVal()) < devOptions_.zeroResistanceTol)
            {
              // change device type to be the level 3 resistor
              type = ModelType::RESISTOR3;
            }
          }
        }
        break;
      }
      currentParam++;
    }
  }

  // Just in case, call the device creator function - if it has already
  // been allocated it won't create a redundant one.

  createDeviceByModelType(type);

  string deviceTypeName;
  // Add an instance of this type.
  DI_ptr = deviceArray_[type]->addInstance(IB, mlData, solState_, extData_,devOptions_, deviceTypeName);

  // check if lead current were requested for this device
  //if( (devicesNeedingLeadCurrentLoads_.find( DI_ptr->getName() ) != devicesNeedingLeadCurrentLoads_.end() ) ||

  //string outputName = DI_ptr->getName();
  string outputName;
  setupIOName( DI_ptr->getName(), outputName);
  if( (devicesNeedingLeadCurrentLoads_.find( outputName ) != devicesNeedingLeadCurrentLoads_.end() ) ||
      (devOptions_.calculateAllLeadCurrents) )
  {
    DI_ptr->enableLeadCurrentCalc();
#ifdef Xyce_DEBUG_DEVICE
    if (devOptions_.debugLevel > 0)
    {
      std::cout << "N_DEV_DeviceMgr::addDeviceInstance Enabling lead current load for device \""
        << DI_ptr->getName()
        << "\" ->  \""
        << outputName
        << "\"" << std::endl;
    }
#endif
  }
#ifdef Xyce_DEBUG_DEVICE
  else
  {
    if (devOptions_.debugLevel > 0)
    {
      std::cout << "N_DEV_DeviceMgr::addDeviceInstance Cannot enable lead current load for device \""
        << DI_ptr->getName()
        << "\" ->  \""
        << outputName
        << "\""
        << std::endl;
    }
  }
#endif

  localDeviceCountMap_[deviceTypeName]++;
  deviceUseFlag_[type] = 1;

  linearSystemFlag_ = linearSystemFlag_ && deviceArray_[type]->isLinearDevice();

  solState_.PDESystemFlag = solState_.PDESystemFlag || PDEDeviceFlag_[type]==1;

#if 0
  std::cout << "In addDeviceInstance. PDESystemFlag = ";
  if ( solState_.PDESystemFlag ) std::cout << "TRUE" <<std::endl;
  else std::cout << "FALSE" <<std::endl;
#endif

  // Set up the instance vectors.  These are the main containers used in
  // the load procedures.  (rather than deviceArray_).
  instancePtrVec_.push_back(DI_ptr);

  // set up the list of pde device instances
  // and the list of non-pde devices instances.
  if (PDEDeviceFlag_[type]==1)
  {
    pdeInstancePtrVec_.push_back(DI_ptr);
  }
  else
  {
    nonPdeInstancePtrVec_.push_back(DI_ptr);
  }

  // set up the list of mosfet instances.
  if (baseType == ModelType::MOSFET1)
  {
    mosfetInstancePtrVec_.push_back(DI_ptr);
  }

  if (baseType == ModelType::BJT )
  {
    bjtInstancePtrVec_.push_back (DI_ptr);
  }


#ifdef Xyce_EXTDEV
  if (type == ModelType::EXTERN_DEVICE)
  {
    extDevInstancePtrVec_.push_back(dynamic_cast<ExternDevice::Instance*>(DI_ptr));
    extDevIBPtrVec_.push_back( new InstanceBlock( IB ) );
  }
#endif

  // set up the independent source map.
  if (type == ModelType::VSRC || type == ModelType::ISRC )
  {
    ExtendedString tmpName(IB.getName());
    tmpName.toUpper ();
    indepSourcePtrMap_[tmpName] = dynamic_cast<SourceInstance*>(DI_ptr);
    indepSourceInstancePtrVec_.push_back( dynamic_cast<SourceInstance*>(DI_ptr) );
  }

  if (type == ModelType::VSRC)
  {
    vsrcInstancePtrVec_.push_back(DI_ptr);
    ExtendedString tmpName(IB.getName());
    tmpName.toUpper ();
    vsrcInstancePtrMap_[tmpName] = dynamic_cast<Vsrc::Instance*>(DI_ptr);
  }

#if 0
  // ERK.  Commenting this out, as the PDE sources cannot be dynamically
  // casted to be SourceInstance classes.  They can't, because they
  // are not derrived from them.  The photocurrent capability in the PDE
  // devices is just an experiment, not a supported one, so I'll figure
  // out a fix for it later.
  if( type == TWO_D_PDE )
  {
    // a photocurrent addition to the two d pde can act like a souce
    // so add it to our list of independent sources
    indepSourceInstancePtrVec_.push_back( dynamic_cast<Source::Instance*>(DI_ptr) );
  }
#endif

  if ( DI_ptr->plotfileFlag () )
  {
    plotFileInstancePtrVec_.push_back(DI_ptr);
  }

  ExtendedString tmpDevName =  devOptions_.testJacDeviceName;
  tmpDevName.toUpper();
  // Set up the vector of devices subject to the jacobian test.
  if( DI_ptr->getName() == tmpDevName )
  {
    testJacDevicePtrVec_.push_back(DI_ptr);
  }

  return DI_ptr;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::deleteDeviceInstance
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/08/01
//-----------------------------------------------------------------------------
bool DeviceMgr::deleteDeviceInstance (const string & name)
{
  bool bsuccess = true;
  bool tmpBool = true;
  int type    = getDeviceIndex (name);

  if (deviceUseFlag_[type] == 1)
  {
    tmpBool = deviceArray_[type]->deleteInstance(name);
    bsuccess = bsuccess && tmpBool;
  }

  // note as this is written it ignores lots of other containers
  // this may be used for clean up at the end of a run, but
  // it is not sufficient to use during simulation setup.
  //
  // need to remove pointer to the instance "name" from other arrays
  // candidate lists are:
  //     vector <DeviceInstance*> instancePtrVec_;
  //     vector <DeviceInstance*> bpInstancePtrVec_; // instances with breakpoints functions
  //     vector <DeviceInstance*> pdeInstancePtrVec_;
  //     vector <DeviceInstance*> nonPdeInstancePtrVec_;
  //     vector <DeviceInstance*> mosfetInstancePtrVec_;
  //     vector <DeviceInstance*> vsrcInstancePtrVec_;
  //     vector <DeviceInstance*> bjtInstancePtrVec_;
  //     map <string, VsrcInstance*> vsrcInstancePtrMap_;
  //
  //     vector <DeviceInstance*> plotFileInstancePtrVec_;
  //
  //     map <string,SourceInstance*> indepSourcePtrMap_;
  //     vector <SourceInstance*> indepSourceInstancePtrVec_;

  string msg("DeviceMgr::deleteDeviceInstance:");
  msg += "  Not ready with the new boilerplate-free device package";
  N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL, msg);

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::printOutLists
// Purpose       : This function will send output to stdout of all the
//                 allocated model and instance lists.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/03/00
//-----------------------------------------------------------------------------
void DeviceMgr::printOutLists()
{
#ifdef Xyce_DEBUG_DEVICE
  if (devOptions_.debugLevel > 0)
  {
    int i;
    for (i = 1; i < ModelType::NUMDEV; ++i)
    {
      if (deviceUseFlag_[i] == 1) deviceArray_[i]->printOutModels(std::cout);
    }
  }
#endif
}

#ifdef Xyce_DEBUG_DEVICE
//-----------------------------------------------------------------------------
// Function      : DeviceMgr::debugOutput1
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 1/29/07
//-----------------------------------------------------------------------------
void DeviceMgr::debugOutput1()
{
   // dump the fvector and the jdxp vector to files.
  if (devOptions_.debugLevel > 3 && solState_.debugTimeFlag)
  {
    // To match the numbering scheme of the NLS debug output files,
    // it is neccessary to add +1 to the newton iterartion number.
    int outIter = solState_.newtonIter + 1;

    // f-vector
    char fn_fv[256]; for (int ich = 0; ich < 256; ++ich) fn_fv[ich] = 0;
    sprintf(fn_fv, "fvector.%03d.txt", outIter);

    // Note: this needs to change sign to match Spice.
    (extData_.fVectorPtr)->scale(-1.0);
    (extData_.fVectorPtr)->writeToFile(fn_fv);
    (extData_.fVectorPtr)->scale(-1.0);

    // jdxp-vector
    char fn_jdxp[256]; for (int ich = 0; ich < 256; ++ich) fn_jdxp[ich] = 0;
    sprintf(fn_jdxp, "Jdxp.%03d.txt", outIter);
    (extData_.JdxpVectorPtr)->writeToFile(fn_jdxp);

#ifdef Xyce_DEBUG_VOLTLIM
    // the voltlim dx vector.
    char fn_dxvl[256]; for (int ich = 0; ich < 256; ++ich) fn_dxvl[ich] = 0;
    sprintf(fn_dxvl, "dxVL.%03d.txt", outIter);
    (extData_.dxVoltlimVectorPtr)->writeToFile(fn_dxvl);

    // jdx2-vector
    char fn_jdx2[256]; for (int ich = 0; ich < 256; ++ich) fn_jdx2[ich] = 0;
    sprintf(fn_jdx2, "Jdx2.%03d.txt", outIter);
    (extData_.Jdx2VectorPtr)->writeToFile(fn_jdx2);
#endif

  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::debugOutput2
// Purpose       : new-dae version
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/19/08
//-----------------------------------------------------------------------------
void DeviceMgr::debugOutput2()
{
  // dump the fvector and the jdxp vector to files.
  if (devOptions_.debugLevel > -1 && solState_.debugTimeFlag)
  {
    // To match the numbering scheme of the NLS debug output files,
    // it is neccessary to add +1 to the newton iterartion number.
    int newtonIter = solState_.newtonIter + 1;
    int outputStepNumber = 0;

    if (solState_.tranopFlag)
    {
      outputStepNumber = 0;
    }
    else if (solState_.initTranFlag)
    {
      outputStepNumber = solState_.timeStepNumber+1;
    }
    else
    {
      outputStepNumber = solState_.timeStepNumber+1;
    }

#ifdef Xyce_DEBUG_VOLTLIM
    // the voltlim dx vector.
    char fn_dxvl[256]; for (int ich = 0; ich < 256; ++ich) fn_dxvl[ich] = 0;
    sprintf(fn_dxvl, "dxVL.%03d.%03d.txt", outputStepNumber, newtonIter);
    (extData_.dxVoltlimVectorPtr)->writeToFile(fn_dxvl);

    // Fdx2-vector
    char fn_fdx2[256]; for (int ich = 0; ich < 256; ++ich) fn_fdx2[ich] = 0;
    sprintf(fn_fdx2, "Fdx2.%03d.%03d.txt", outputStepNumber, newtonIter);
    (extData_.Fdx2VectorPtr)->writeToFile(fn_fdx2);

    // Qdx2-vector
    char fn_qdx2[256]; for (int ich = 0; ich < 256; ++ich) fn_qdx2[ich] = 0;
    sprintf(fn_qdx2, "Qdx2.%03d.%03d.txt", outputStepNumber, newtonIter);
    (extData_.Qdx2VectorPtr)->writeToFile(fn_qdx2);
#endif

  }
}
#endif

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setInitialGuess
// Purpose       : This is a function call that sets the initial guess for
// devices that have initial guesses.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/02/03
//-----------------------------------------------------------------------------
bool DeviceMgr::setInitialGuess (N_LAS_Vector * solVectorPtr)
{
  bool bsuccess = true;

  if (solVectorPtr != 0)
  {
    extData_.nextSolVectorPtr = solVectorPtr;

    // if two-level, and just the inner problem, only load the PDE devices.
    vector<DeviceInstance*>::iterator iter;
    vector<DeviceInstance*>::iterator begin;
    vector<DeviceInstance*>::iterator end;

    begin = pdeInstancePtrVec_.begin ();
    end = pdeInstancePtrVec_.end ();
    for (iter=begin; iter!=end;++iter)
    {
      bool tmpBool = (*iter)->setInitialGuess ();
      bsuccess = bsuccess && tmpBool;
    }
  }

  return bsuccess;
}

#ifdef Xyce_EXTDEV
//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setUpPassThroughParamsMap_()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date :
//-----------------------------------------------------------------------------
void DeviceMgr::setUpPassThroughParamsMap_()
{
  passThroughParamsMap_[ "MOSFET:GAINSCALE"       ] = 1;
  passThroughParamsMap_[ "MOSFET:GAIN"            ] = 1;
  passThroughParamsMap_[ "MOSFET:NLTERMSCALE"     ] = 1;
  passThroughParamsMap_[ "MOSFET:NLTERM"          ] = 1;
  passThroughParamsMap_[ "MOSFET_ALL:GAINSCALE"   ] = 1;
  passThroughParamsMap_[ "MOSFET_ALL:NLTERMSCALE" ] = 1;
  passThroughParamsMap_[ "MOSFET1:GAINSCALE"      ] = 1;
  passThroughParamsMap_[ "MOSFET1:NLTERMSCALE"    ] = 1;
  passThroughParamsMap_[ "MOSFET:W"               ] = 1;
  passThroughParamsMap_[ "MOSFET:L"               ] = 1;
  passThroughParamsMap_[ "MOSFET:SIZESCALE"       ] = 1;
  passThroughParamsMap_[ "MOSFET:TOX"             ] = 1;
  passThroughParamsMap_[ "TEMP"                   ] = 1;
  passThroughParamsMap_[ "BJT:NF"                 ] = 1;
  passThroughParamsMap_[ "BJT:NR"                 ] = 1;
  passThroughParamsMap_[ "BJT:EXPORD"             ] = 1;
  //passThroughParamsMap_[ "GSTEPPING"              ] = 1;
}
#endif

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setParam
//
// Purpose       : This function sets named parameters (name) to a
//                 specified value (val).
//
// Special Notes : Used for continuation calculations, as well as possibly
//                 intrusive sensitivity/optimization calculations.  It is
//                 assumed that this is called after everything (devices,
//                 solvers, etc.) is set up.
//
//                 The specified parameter can be either a natural or
//                 artificial parameter.
//
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/02/03
//-----------------------------------------------------------------------------
bool DeviceMgr::setParam (string & name, double val)
{
  bool bsuccess = true, success = true;
  vector<DeviceInstance*>::const_iterator iter;
  vector<DeviceInstance*>::const_iterator begin;
  vector<DeviceInstance*>::const_iterator end;

  vector<DeviceModel*>::const_iterator iterM;
  vector<DeviceModel*>::const_iterator beginM;
  vector<DeviceModel*>::const_iterator endM;

  ExtendedString tmpName(name);
  tmpName.toUpper ();

#ifdef Xyce_DEBUG_DEVICE
  if (devOptions_.debugLevel > 0 )
  {
    string netListFile("");
    if (commandLine_.getArgumentValue(string("netlist")) != "")
    {
      netListFile = commandLine_.getArgumentValue(string("netlist"));
    }

    std::cout << netListFile << "\t\t";
    std::cout << "DeviceMgr::setParam.  ";
    std::cout << name;
    std::cout << "  " << val;
    std::cout << std::endl;
  }
#endif

  // There are numerous special cases.  Check these first:
  // check if this is one of the designated artificial parameters.
  if (tmpName == "MOSFET:GAINSCALE" ||
      tmpName == "MOSFET:GAIN")
  {
    solState_.gainScale[0] = val;
    solState_.artParameterFlag = true;
  }
  else if (tmpName == "MOSFET:NLTERMSCALE" ||
           tmpName == "MOSFET:NLTERM")
  {
    solState_.nltermScale = val;
    solState_.artParameterFlag = true;
  }
#if 0
  else if (tmpName == "MOSFET_ALL:GAINSCALE")
  {
    solState_.gainScale[0] = val;
    solState_.mos1GainScale = val;
    solState_.mos3GainScale = val;
    solState_.artParameterFlag = true;
    solState_.mos1ArtParameterFlag = true;
    solState_.mos3ArtParameterFlag = true;
  }
  else if (tmpName == "MOSFET_ALL:NLTERMSCALE")
  {
    solState_.nltermScale = val;
    solState_.mos1NltermScale = val;
    solState_.mos3NltermScale = val;
    solState_.artParameterFlag = true;
    solState_.mos1ArtParameterFlag = true;
    solState_.mos3ArtParameterFlag = true;
  }
  else if (tmpName == "MOSFET_13:GAINSCALE")
  {
    solState_.mos1GainScale = val;
    solState_.mos3GainScale = val;
    solState_.mos1ArtParameterFlag = true;
    solState_.mos3ArtParameterFlag = true;
  }
  else if (tmpName == "MOSFET_13:NLTERMSCALE")
  {
    solState_.mos1NltermScale = val;
    solState_.mos3NltermScale = val;
    solState_.mos1ArtParameterFlag = true;
    solState_.mos3ArtParameterFlag = true;
  }
  else if (tmpName == "MOSFET1:GAINSCALE")
  {
    solState_.mos1GainScale = val;
    solState_.mos1ArtParameterFlag = true;
  }
  else if (tmpName == "MOSFET1:NLTERMSCALE")
  {
    solState_.mos1NltermScale = val;
    solState_.mos1ArtParameterFlag = true;
  }
  else if (tmpName == "MOSFET3:GAINSCALE")
  {
    solState_.mos3GainScale = val;
    solState_.mos3ArtParameterFlag = true;
  }
  else if (tmpName == "MOSFET3:NLTERMSCALE")
  {
    solState_.mos3NltermScale = val;
    solState_.mos3ArtParameterFlag = true;
  }
#endif
  else if (tmpName == "MOSFET:GAINSCALE_BLOCK_0")
  {
    solState_.gainScale[0] = val;
    solState_.artParameterFlag = true;
  }
  else if (tmpName == "MOSFET:GAINSCALE_BLOCK_1")
  {
    solState_.gainScale[1] = val;
    solState_.artParameterFlag = true;
  }
  else if (tmpName == "MOSFET:GAINSCALE_BLOCK_2")
  {
    solState_.gainScale[2] = val;
    solState_.artParameterFlag = true;
  }
  else if (tmpName == "MOSFET:GAINSCALE_BLOCK_3")
  {
    solState_.gainScale[3] = val;
    solState_.artParameterFlag = true;
  }
  else if (tmpName == "MOSFET:GAINSCALE_BLOCK_4")
  {
    solState_.gainScale[4] = val;
    solState_.artParameterFlag = true;
  }
  else if (tmpName == "MOSFET:GAINSCALE_BLOCK_5")
  {
    solState_.gainScale[5] = val;
    solState_.artParameterFlag = true;
  }
  else if (tmpName == "MOSFET:GAINSCALE_BLOCK_6")
  {
    solState_.gainScale[6] = val;
    solState_.artParameterFlag = true;
  }
  else if (tmpName == "MOSFET:GAINSCALE_BLOCK_7")
  {
    solState_.gainScale[7] = val;
    solState_.artParameterFlag = true;
  }
  else if (tmpName == "MOSFET:GAINSCALE_BLOCK_8")
  {
    solState_.gainScale[8] = val;
    solState_.artParameterFlag = true;
  }
  else if (tmpName == "MOSFET:GAINSCALE_BLOCK_9")
  {
    solState_.gainScale[9] = val;
    solState_.artParameterFlag = true;
  }
  else if (tmpName == "MOSFET:L")
  {
    solState_.sizeParameterFlag = true;

    double newL = val;

    // now loop over all the mosfet instances, and change the l.
    begin = mosfetInstancePtrVec_.begin ();
    end   = mosfetInstancePtrVec_.end ();
    for (iter=begin; iter!=end;++iter)
    {
      string nameL("l");
      success = (*iter)->setParam (nameL, newL);
      success = success && (*iter)->processParams ();
    }

  }
  else if (tmpName == "MOSFET:W")
  {
    solState_.sizeParameterFlag = true;
    double newW = val;

    // now loop over all the mosfet instances, and change the w.
    begin = mosfetInstancePtrVec_.begin ();
    end   = mosfetInstancePtrVec_.end ();
    for (iter=begin; iter!=end;++iter)
    {
      string nameW("w");
      success = (*iter)->setParam (nameW, newW);
      success = success && (*iter)->processParams ();
    }

  }
  else if (tmpName == "MOSFET:SIZESCALE")
  {
    solState_.sizeParameterFlag = true;
    solState_.sizeScale = val;

    // What we want at val = 0.0.  For all the lengths and widths to have a
    // ratio of:  L=5u W=175u  or L=5u W=270u.  Compromise.  L/W = 5/200,
    // with L actually being 5u.

    double length0 = devOptions_.length0;
    double width0  = devOptions_.width0;
    // now loop over all the mosfet instances, and change the l, w.
    begin = mosfetInstancePtrVec_.begin ();
    end   = mosfetInstancePtrVec_.end ();
    for (iter=begin; iter!=end;++iter)
    {
      string nameL("l");
      string nameW("w");

      success = (*iter)->scaleParam(nameL, solState_.sizeScale, length0);
      success = success || (*iter)->scaleParam(nameW, solState_.sizeScale, width0);
      success = success && (*iter)->processParams ();
    }

  }
  else if (tmpName == "MOSFET:TOX")
  {
    solState_.sizeParameterFlag = true;
    solState_.sizeScale = val;

    // What we want at val = 0.0.  For all the lengths and widths to have a
    // ratio of:  L=5u W=175u  or L=5u W=270u.  Compromise.  L/W = 5/200,
    // with L actually being 5u.

    double tox0    = devOptions_.tox0;

    // loop over all the models and change tox.
    beginM = mosfetModelPtrVec_.begin ();
    endM   = mosfetModelPtrVec_.end ();
    for (iterM=beginM; iterM!=endM;++iterM)
    {
      string nameTOX("tox");

      success = (*iterM)->scaleParam(nameTOX, solState_.sizeScale, tox0);

      success = success && (*iterM)->processParams ();
      success = success && (*iterM)->processInstanceParams ();
    }

  }
  else if (tmpName == "VSRCSCALE") // This scalar, val, is assumed to go from 0 to 1
  {
    int vsrcSize = vsrcInstancePtrVec_.size();
    for (int i=0;i<vsrcSize;++i)
    {
      bool s1 = vsrcInstancePtrVec_[i]->scaleDefaultParam (val);
      bool s2 = vsrcInstancePtrVec_[i]->processParams ();
      success = success && s1 && s2;
    }
  }
  else if (tmpName == "TEMP")
  {
    updateTemperature (val);
  }
  // if this is called, need to be running 2-level newton, and
  // need to have called "enablePDEContinuation" first.
  else if (tmpName == "PDEALPHA")
  {
    solState_.pdeAlpha = val; // not important - part of planned refactor.

    if (!solState_.PDEcontinuationFlag)
    {
      string msg("DeviceMgr::setParam:");
      msg += " tried to set pdeAlpha without first calling";
      msg += " enablePDEContinaution.";
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL, msg);
    }

    // This should be part of the inner level solve for a 2-level
    // algorithm.  Set this for the PDE devices.
    vector<DeviceInstance*>::iterator iter;
    vector<DeviceInstance*>::iterator begin;
    vector<DeviceInstance*>::iterator end;
    begin = instancePtrVec_.begin ();
    end = instancePtrVec_.end ();

    for (iter=begin; iter!=end;++iter)
    {
      (*iter)->setPDEContinuationAlpha (val);
    }
  }
  else if (tmpName == "PDEBETA")
  {
    solState_.PDEcontinuationFlag = true;

    vector<DeviceInstance*>::iterator iter;
    vector<DeviceInstance*>::iterator begin;
    vector<DeviceInstance*>::iterator end;
    begin = instancePtrVec_.begin ();
    end = instancePtrVec_.end ();

    for (iter=begin; iter!=end;++iter)
    {
      (*iter)->setPDEContinuationBeta (val);
    }
  }
  else if (tmpName == "PDECHARGEALPHA")
  {
    solState_.chargeAlpha = val;
    solState_.chargeHomotopy = true;
  }
  else if (tmpName == "BJT:BF")
  {
    double scale = val;
    string newBF("bf");

    // loop over all the models and change bf
    beginM = bjtModelPtrVec_.begin ();
    endM   = bjtModelPtrVec_.end ();
    for (iterM=beginM; iterM!=endM;++iterM)
    {
      success = (*iterM)->scaleParam(newBF, scale, 0.0);
      success = success && (*iterM)->processParams();
      success = success && (*iterM)->processInstanceParams();
    }
  }
  else if (tmpName == "BJT:NF")
  {
    double scale = val;
    string newNF("nf");

    // loop over all the models and change nf
    beginM = bjtModelPtrVec_.begin ();
    endM   = bjtModelPtrVec_.end ();
    for (iterM=beginM; iterM!=endM;++iterM)
    {
      success = (*iterM)->scaleParam(newNF, scale, 10.0);
      success = success && (*iterM)->processParams();
      success = success && (*iterM)->processInstanceParams();
    }
  }
  else if (tmpName == "BJT:NR")
  {
    double scale = val;
    string newNR("nr");

    // loop over all the models and change nr
    beginM = bjtModelPtrVec_.begin ();
    endM   = bjtModelPtrVec_.end ();
    for (iterM=beginM; iterM!=endM;++iterM)
    {
      success = (*iterM)->scaleParam(newNR, scale, 10.0);
      success = success && (*iterM)->processParams();
      success = success && (*iterM)->processInstanceParams();
    }
  }
  else if (tmpName == "BJT:EXPORD")
  {
    devOptions_.exp_order = val;
    solState_.bjtArtParameterFlag = true;
  }
  else if (tmpName == "GSTEPPING")
  {
    // Do nothing!!!! This is not a device variable but is used by the
    // AugmentLinSys.
  }
  else if (tmpName == "DIODE:N")
  {
    double scale = val;
    string newN("n");

    // loop over all the models and change n
    beginM = diodeModelPtrVec_.begin ();
    endM   = diodeModelPtrVec_.end ();
    for (iterM=beginM; iterM!=endM;++iterM)
    {
      success = (*iterM)->scaleParam(newN, scale, 10.0);
      success = success && (*iterM)->processParams();
      success = success && (*iterM)->processInstanceParams();
    }
  }
  else if (tmpName == "GMIN")
  {
    devOptions_.gmin = devOptions_.gmin_orig *val + devOptions_.gmin_init * (1.0-val);
  }
  else if (solState_.global_params.find(tmpName) != solState_.global_params.end())
  {
    if (solState_.global_params[tmpName] != val)
    {
      solState_.global_params[tmpName] = val;
      vector<DeviceEntity*>::iterator iter;
      vector<DeviceEntity*>::iterator begin = dependentPtrVec_.begin();
      vector<DeviceEntity*>::iterator end = dependentPtrVec_.end();
      for (iter=begin; iter!=end;++iter)
      {
        if ((*iter)->updateGlobalParameters (solState_.global_params));
        {
          (*iter)->processParams();
          (*iter)->processInstanceParams();
        }
      }
    }
  }
  else
  {
    // If not artificial, then search for the appropriate natural param(s).

    // the devSensPtr will be deleted in the destructor.
    if (devSensPtr_==0)
    {
      devSensPtr_ = new DeviceSensitivities
        (devOptions_,extData_, solState_, *lasSysPtr_);
    }

    string paramName;
    string entityName;
    devSensPtr_->stripParamName (name, paramName, entityName);

    DeviceEntity * dePtr = devSensPtr_->getDeviceEntity(name);

#ifdef Xyce_PARALLEL_MPI
    double foundParam = 0.0;
    double finalParam = 0.0;
#endif
    bool entityFound = (dePtr!=0)?true:false;

    if (entityFound)
    {
      bool found;
      if (paramName == "")
      {
        found = dePtr->setDefaultParam (val);
      }
      else
      {
        found = dePtr->setParam (paramName, val);
      }
      if (found)
      {
        dePtr->processParams (); // if this "entity" is a model, then need to
                                 // also do a  "processParams" on the related
                                 // instances.
        dePtr->processInstanceParams();
      }

#ifdef Xyce_PARALLEL_MPI
      foundParam = found?1.0:0.0;
#endif
      entityFound = found;
    }

#ifdef Xyce_PARALLEL_MPI
    N_PDS_Comm * pdsCommPtr = pdsMgrPtr_->getPDSComm();
    pdsCommPtr->barrier();
    pdsCommPtr->sumAll(&foundParam, &finalParam, 1);
    entityFound = (finalParam != 0.0)?true:false;
#endif

    if(!entityFound)
    {
      string msg("Unable to find parameter " + name + "\n");
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg);
    }
  }

   // Certain parameters should be passed through to the inner solve,
   // if there is an inner circuit problem.  The names of parameters
   // that should be passed through are stored in the map.
#ifdef Xyce_EXTDEV
  if ( passThroughParamsMap_.find(tmpName) != passThroughParamsMap_.end())
  {

    vector<ExternDevice::Instance*>::iterator iter;
    vector<ExternDevice::Instance*>::iterator begin;
    vector<ExternDevice::Instance*>::iterator end;
    begin = extDevInstancePtrVec_.begin ();
    end = extDevInstancePtrVec_.end ();

    for (iter=begin; iter!=end;++iter)
    {
      bool bs1 = (*iter)->setInternalParam (name, val);
      //bsuccess = bsuccess && bs1;
    }

  }
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getParam
// Purpose       : Returns the current value of a named parameter.
//
// Special Notes : This works in parallel.
//
//                 If the parameter is not found, this function sets val to
//                 0.0, and returns a "false" boolean.  It does not invoke the
//                 error handler.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 9/13/04
//-----------------------------------------------------------------------------
bool DeviceMgr::getParam (const string & name, double & val)
{
  val = 0.0;
  ExtendedString tmpName(name);
  tmpName.toUpper ();
  bool entityFound = false;


  // check if this is one of the designated artificial parameters.
  if (tmpName == "MOSFET:GAINSCALE")
  {
    val = solState_.gainScale[0];
    entityFound = true;
  }
  else if (tmpName == "MOSFET:NLTERMSCALE")
  {
    val = solState_.nltermScale;
    entityFound = true;
  }
  else if (tmpName == "MOSFET:L")
  {
    val = devOptions_.defl;
    entityFound = true;
  }
  else if (tmpName == "MOSFET:W")
  {
    val = devOptions_.defw;
    entityFound = true;
  }
  else if (tmpName == "TEMP")
  {
    val = devOptions_.temp.dVal();;
    val -= CONSTCtoK;
    entityFound = true;
  }
  else if (solState_.global_params.find(tmpName) !=
      solState_.global_params.end())
  {
    val = solState_.global_params[tmpName];
    entityFound = true;
  }
  else
  {
    // the devSensPtr will be deleted in the destructor.
    if (devSensPtr_==0)
    {
      devSensPtr_ = new DeviceSensitivities
        (devOptions_,extData_, solState_, *lasSysPtr_);
    }

    string paramName;
    string entityName;
    devSensPtr_->stripParamName (name, paramName, entityName);


    DeviceEntity * dePtr = devSensPtr_->getDeviceEntity(name);
    entityFound = (dePtr!=0)?true:false;

    if (entityFound)
    {
      entityFound = dePtr->getParam (paramName, val);
      if (paramName == "")
        entityFound = true;
    }

    if (!entityFound)
      outputMgrPtr_->getMeasureValue(tmpName, val, entityFound);

#ifdef Xyce_PARALLEL_MPI
    double foundParam = entityFound?1:0;
    double finalParam = 0.0;
    double globalVal;
    N_PDS_Comm * pdsCommPtr = pdsMgrPtr_->getPDSComm();
    pdsCommPtr->barrier();
    pdsCommPtr->sumAll(&foundParam, &finalParam, 1);
    if (finalParam != 0.0)
    {
      entityFound = true;
      pdsCommPtr->sumAll(&val, &globalVal, 1);
      val = globalVal/finalParam;
    }
#endif
  }

  return entityFound;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getParam
// Purpose       : Returns the current value of a named parameter.
//
// Special Notes : This works in parallel.
//
//                 This is different than the other "getParam" function,
//                 in that it generates a fatal error if the parameter is
//                 not found.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/26/03
//-----------------------------------------------------------------------------
double DeviceMgr::getParam (const string & name)
{
  bool bsuccess = false;
  double val = 0.0;
  bsuccess = getParam(name,val);

  if(!bsuccess)
  {
    string msg("Unable to find parameter " + name + "\n");
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg);
  }

  return val;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getVsrcLIDs
//
// Purpose       : Returns the LID of the voltage drop row from the named
//                 voltage source.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/26/03
//-----------------------------------------------------------------------------
bool DeviceMgr::getVsrcLIDs (
    string & srcName, int & li_Pos, int & li_Neg, int & li_Bra)
{
  ExtendedString tmpName(srcName);
  tmpName.toUpper();
  map <string, N_DEV_VsrcInstance*>::iterator iterVsrc = vsrcInstancePtrMap_.find(tmpName);
  if ( iterVsrc != vsrcInstancePtrMap_.end())
  {
    N_DEV_VsrcInstance * vsrcPtr = iterVsrc->second;
    vsrcPtr->getLIDs(li_Pos,li_Neg,li_Bra);
  }
  else
  {
    string msg("DeviceMgr::getVoltageDropRow ");
    msg += "Unable to find source: " + srcName + "\n";
#ifdef Xyce_PARALLEL_MPI
#else
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL, msg);
#endif
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::updateState
// Purpose       : This should be called prior to loadDAEVectors.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/17/05
//-----------------------------------------------------------------------------
bool DeviceMgr::updateState (
     N_LAS_Vector * nextSolVectorPtr,
     N_LAS_Vector * currSolVectorPtr,
     N_LAS_Vector * lastSolVectorPtr,
     N_LAS_Vector * nextStaVectorPtr,
     N_LAS_Vector * currStaVectorPtr,
     N_LAS_Vector * lastStaVectorPtr,
     N_LAS_Vector * nextStoVectorPtr,
     N_LAS_Vector * currStoVectorPtr,
     N_LAS_Vector * lastStoVectorPtr
     )
{
  bool bsuccess = true;
  bool tmpBool = true;

#ifdef Xyce_PARALLEL_MPI
  N_ERH_ErrorMgr::startSafeBarrier();
#endif

  tmpBool = setupSolverInfo_();
  bsuccess = bsuccess && tmpBool;

  // copy over the passed pointers:
  extData_.nextSolVectorPtr = nextSolVectorPtr;
  extData_.currSolVectorPtr = currSolVectorPtr;
  extData_.lastSolVectorPtr = lastSolVectorPtr;
  extData_.nextStaVectorPtr = nextStaVectorPtr;
  extData_.currStaVectorPtr = currStaVectorPtr;
  extData_.lastStaVectorPtr = lastStaVectorPtr;
  extData_.nextStoVectorPtr = nextStoVectorPtr;
  extData_.currStoVectorPtr = currStoVectorPtr;
  extData_.lastStoVectorPtr = lastStoVectorPtr;

#ifdef Xyce_PARALLEL_MPI
  extData_.nextSolVectorPtr->importOverlap();
#endif

  // Now reset the relevant RAW pointers:
  extData_.nextSolVectorRawPtr = &((*extData_.nextSolVectorPtr)[0]);
  extData_.currSolVectorRawPtr = &((*extData_.currSolVectorPtr)[0]);
  extData_.lastSolVectorRawPtr = &((*extData_.lastSolVectorPtr)[0]);
  extData_.nextStaVectorRawPtr = &((*extData_.nextStaVectorPtr)[0]);
  extData_.currStaVectorRawPtr = &((*extData_.currStaVectorPtr)[0]);
  extData_.lastStaVectorRawPtr = &((*extData_.lastStaVectorPtr)[0]);
  extData_.nextStoVectorRawPtr = &((*extData_.nextStoVectorPtr)[0]);
  extData_.currStoVectorRawPtr = &((*extData_.currStoVectorPtr)[0]);
  extData_.lastStoVectorRawPtr = &((*extData_.lastStoVectorPtr)[0]);

#ifdef Xyce_OLD_LOOPING
  // call all the intermediate vars loads:
  vector<DeviceInstance*>::iterator iter;
  vector<DeviceInstance*>::iterator begin;
  vector<DeviceInstance*>::iterator end;

  // if inner problem, only do the PDE loads.
  if (solState_.twoLevelNewtonCouplingMode==INNER_PROBLEM)
  {
    begin = pdeInstancePtrVec_.begin ();
    end   = pdeInstancePtrVec_.end ();
  }
  else
  {
    begin = instancePtrVec_.begin ();
    end   = instancePtrVec_.end ();
  }

  updateDependentParameters_();

  // updateIntermediateVars_();

  for (iter=begin; iter!=end;++iter)
  {
    tmpBool = (*iter)->updatePrimaryState ();
    bsuccess = bsuccess && tmpBool;
  }
#else

  updateDependentParameters_();
  // updateIntermediateVars_();

  // if inner problem, only do the PDE loads.
  if (solState_.twoLevelNewtonCouplingMode==INNER_PROBLEM)
  {
    // call all the intermediate vars loads:
    vector<Device*>::iterator iter;
    vector<Device*>::iterator begin = pdeDevicePtrVec_.begin ();
    vector<Device*>::iterator end = pdeDevicePtrVec_.end ();
    for (iter=begin; iter!=end;++iter)
    {
      tmpBool = (*iter)->updateState ( extData_.nextSolVectorRawPtr,
                                      extData_.nextStaVectorRawPtr, extData_.nextStoVectorRawPtr);
      bsuccess = bsuccess && tmpBool;
    }
  }
  else
  {
    int numDevices = devicePtrVec_.size();
    for(int i=0; i< numDevices; ++i)
    {
      bsuccess=devicePtrVec_.at(i)->updateState ( extData_.nextSolVectorRawPtr,
                                                 extData_.nextStaVectorRawPtr, extData_.nextStoVectorRawPtr);
    }
  }

#endif

#ifdef Xyce_EXTDEV
  updateExternalDevices_();
#endif

#ifdef Xyce_PARALLEL_MPI
  extData_.nextStaVectorPtr->importOverlap();
  extData_.nextStoVectorPtr->importOverlap();
#endif

#ifdef Xyce_PARALLEL_MPI
  N_ERH_ErrorMgr::safeBarrier(0);
#endif

  return true;
}


//-----------------------------------------------------------------------------
// Function      : DeviceMgr::loadDAEMatrices
// Purpose       : This function loads the various DAE related matrices to
//                 set up the following expression:
//
//                 residual:  f(x) = dQ/dt + F(x) - B(t) = 0
//
//                 jacobian:  J(x) = d(dQ/dt)dx + dFdx
//                                 = d(dQdx)dt + dFdx
//
// Special Notes : ERK.  This function is only called if using the new
//                 (new-DAE) integrator.  As such, it is also used for
//                 MPDE.  It is not used for the old (ODE) integrator.
//                 The corresponding function for the old integrator
//                 is loadJacobianMatrix
//
//                 Note that this function, unlike the loadJacobianMatrix
//                 function *requires* that vector/matrix pointers
//                 be passed in to the function.  When
//                 running in new-DAE or MPDE, we can't rely (much) on
//                 registered vectors. (set up in initializeAll).  In
//                 particular, MPDE cycles through different vector and
//                 matrix blocks, which each correspond to different
//                 "fast" time points, and for each block the state
//                 information is different.
//
//                 NOTE:  This function assumes that all the load matrices
//                 are zeroed out.  That is, dFdx, dQdx are all zeros.
//
//                 If that isn't true, then this function will not produce
//                 the correct answer.  The reason for doing this
//                 is MPDE.  For the warped case, the entire linear system
//                 needs to be zeroed out, but the full system includes
//                 phase equations that never get passed in here.
//                 That zeroing now happens upstream from here, in either
//                 the MPDE loader or the time integrator.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/03/04
//-----------------------------------------------------------------------------
bool DeviceMgr::loadDAEMatrices
  (N_LAS_Vector * tmpSolVectorPtr,
   N_LAS_Vector * tmpStateVectorPtr,
   N_LAS_Vector * tmpStateDerivVectorPtr,
   N_LAS_Vector * tmpStoreVectorPtr,
   N_LAS_Matrix * tmpdQdxMatrixPtr,
   N_LAS_Matrix * tmpdFdxMatrixPtr)
{
  bool bsuccess = true;
  bool tmpBool = true;

#ifdef Xyce_PARALLEL_MPI
  N_ERH_ErrorMgr::startSafeBarrier();
#endif

  // copy over the passed pointers:
  (extData_.nextSolVectorPtr) = tmpSolVectorPtr;
  bool resetRawMatrixPointers = true;
  if (
      (extData_.dQdxMatrixPtr == tmpdQdxMatrixPtr) &&
      (extData_.dFdxMatrixPtr == tmpdFdxMatrixPtr)
     )
  {
    resetRawMatrixPointers = false;
  }
  if (resetRawMatrixPointers)
  {
    extData_.dQdxMatrixPtr = tmpdQdxMatrixPtr;
    extData_.dFdxMatrixPtr = tmpdFdxMatrixPtr;
  }

  extData_.nextStaVectorPtr = tmpStateVectorPtr;
  extData_.nextStaDerivVectorPtr = tmpStateDerivVectorPtr;
  extData_.nextStoVectorPtr = tmpStoreVectorPtr;

  // setup the relevant RAW vector pointers:
  extData_.nextSolVectorRawPtr = &((*extData_.nextSolVectorPtr)[0]);
  extData_.nextStaVectorRawPtr = &((*extData_.nextStaVectorPtr)[0]);
  extData_.nextStaDerivVectorRawPtr = &((*extData_.nextStaDerivVectorPtr)[0]);
  extData_.nextStoVectorRawPtr = &((*extData_.nextStoVectorPtr)[0]);

//#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  // setup the relevant RAW matrix pointers (down in the devices that need them):
  if (resetRawMatrixPointers || solState_.blockAnalysisFlag)
  {
    this->setupRawMatrixPointers_();
  }
//#endif

  vector<DeviceInstance*>::iterator iter;
  vector<DeviceInstance*>::iterator begin;
  vector<DeviceInstance*>::iterator end;

  // if this is an "inner problem" phase of a Two-Level Newton
  // simulation, then only load the PDE devices.  Everything else just
  // gets "1" on the diagonal.
  //
  // Note, it is possible to just load 1's using a petra call, so I may
  // get rid of the "trivial" matrix stamp stuff soon.

  if (solState_.twoLevelNewtonCouplingMode==INNER_PROBLEM)
  {
#if 1
    begin = nonPdeInstancePtrVec_.begin();
    end   = nonPdeInstancePtrVec_.end();
    for (iter=begin; iter!=end;++iter)
    {
      (*iter)->loadTrivialDAE_FMatrixStamp ();
    }
#else
    diagonalVectorPtr_->putScalar(1.0);
    extData_.dFdxMatrixPtr->replaceDiagonal ( (*diagonalVectorPtr_) );

    begin = pdeInstancePtrVec_.begin();
    end   = pdeInstancePtrVec_.end();
    for (iter=begin; iter!=end;++iter)
    {
      (*iter)->zeroMatrixDiagonal(extData_.dFdxMatrixPtr);
    }
#endif

    begin = pdeInstancePtrVec_.begin();
    end   = pdeInstancePtrVec_.end();
    for (iter=begin; iter!=end;++iter)
    {
      tmpBool = (*iter)->loadDAEdQdx();bsuccess = bsuccess && tmpBool;
      tmpBool = (*iter)->loadDAEdFdx();bsuccess = bsuccess && tmpBool;
    }
  }
  // Else, do a normal analytical matrix load.
  else
  {
#ifdef Xyce_OLD_LOOPING
    begin = instancePtrVec_.begin();
    end   = instancePtrVec_.end();
    for (iter=begin; iter!=end;++iter)
    {
      tmpBool = (*iter)->loadDAEdQdx();bsuccess = bsuccess && tmpBool;
      tmpBool = (*iter)->loadDAEdFdx();bsuccess = bsuccess && tmpBool;
    }
#else
    int numDevices = devicePtrVec_.size();
    for(int i=0; i< numDevices; ++i)
    {
      bsuccess=devicePtrVec_.at(i)->loadDAEMatrices ( *(extData_.dFdxMatrixPtr) , *(extData_.dQdxMatrixPtr) );
    }
#endif
  }

  // Run jacobian diagnostic.
  if (devOptions_.testJacobianFlag &&
      (solState_.timeStepNumber >= devOptions_.testJacStartStep &&
       solState_.timeStepNumber <= devOptions_.testJacStopStep)
      )
  {
#ifdef Xyce_DEBUG_DEVICE
    devOptions_.debugLevel -= 2;
#endif
    // Test just the specified device(s), if the user specified any.
    if( (devOptions_.testJacDeviceNameGiven) )
    {
      begin = testJacDevicePtrVec_.begin();
      end   = testJacDevicePtrVec_.end();
    }
    else // Test all the devices:
    {
      begin = instancePtrVec_.begin();
      end   = instancePtrVec_.end();
    }

    for (iter=begin; iter!=end;++iter)
    {
      (*iter)->testDAEMatrices (nameVec_);
    }
#ifdef Xyce_DEBUG_DEVICE
    devOptions_.debugLevel += 2;
#endif
  }

#if 0
  // figure this stuff out later...
  // Add in terms for .IC initial conditions if DC OP
  if (solState_.twoLevelNewtonCouplingMode==INNER_PROBLEM)
  {
    if (icLoads_ != NULL && solState_.dcopFlag)
    {
      double diag = 10000.0;
      int size = icLoads_->size();
      for (int i = 0; i < size; ++i)
        extData_.JMatrixPtr->sumIntoRow( (*icLoads_)[i].first, 1, &diag,
                                         &((*icLoads_)[i].first) );
    }
  }
#endif

  //Tell Jacobian, fill is complete allowing accumulation if necessary
  extData_.dQdxMatrixPtr->fillComplete();
  extData_.dFdxMatrixPtr->fillComplete();

#ifdef Xyce_PARALLEL_MPI
  N_ERH_ErrorMgr::safeBarrier(0);
#endif

#ifdef Xyce_DEBUG_DEVICE
  if (devOptions_.debugLevel > 1 && solState_.debugTimeFlag)
  {
    int newtonIter = solState_.newtonIter;
    std::cout << "-----------------------------------------------" << std::endl;
    std::cout <<  "Q-matrix: nonlinear iteration = " << newtonIter << "\n";
    //tmpdQdxMatrixPtr->printPetraObject();
    extData_.dQdxMatrixPtr->printPetraObject();
    std::cout << std::endl;
    std::cout << "-----------------------------------------------" << std::endl;
    std::cout <<  "F-matrix: nonlinear iteration = " << newtonIter << "\n";
    //tmpdFdxMatrixPtr->printPetraObject();
    extData_.dFdxMatrixPtr->printPetraObject();
    std::cout << std::endl;
    std::cout << "-----------------------------------------------" << std::endl;
  }
#endif


  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::loadDAEVectors
// Purpose       : This function loads the various DAE related vectors to
//                 set up the following expression:
//
//                 f(x) = dQ/dt + F(x) - B(t) = 0
//
// Special Notes : ERK.  This function is only called if using the new
//                 (new-DAE) integrator.  As such, it is also used for
//                 MPDE.  It is not used for the old (ODE) integrator.
//                 The corresponding function for the old integrator
//                 is loadRHSVector.
//
//                 Note that this function, unlike the loadRHSVector
//                 function *requires* that vectors be passed in.  When
//                 running in new-DAE or MPDE, we can't rely (much) on
//                 registered vectors. (set up in initializeAll).  In
//                 particular, MPDE cycles through different vector and
//                 matrix blocks, which each correspond to different
//                 "fast" time points, and for each block the state
//                 information is different.
//
//                 NOTE:  This function assumes that all the load vectors
//                 are zeroed out.  That is, F, Q, B, dFdxdVp and dQdxVp
//                 are all zeros.
//
//                 If that isn't true, then this function will not produce
//                 the correct answer.  The reason for doing this
//                 is MPDE.  For the warped case, the entire linear system
//                 needs to be zeroed out, but the full system includes
//                 phase equations that never get passed in here.
//
//                 That zeroing now happens upstream from here, in either
//                 the MPDE loader or the time integrator.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/26/03
//-----------------------------------------------------------------------------
bool DeviceMgr::loadDAEVectors
  (N_LAS_Vector * tmpSolVectorPtr,
   N_LAS_Vector * tmpCurrSolVectorPtr,
   N_LAS_Vector * tmpLastSolVectorPtr,
   N_LAS_Vector * tmpStaVectorPtr,
   N_LAS_Vector * tmpCurrStaVectorPtr,
   N_LAS_Vector * tmpLastStaVectorPtr,
   N_LAS_Vector * tmpStaDerivVectorPtr,
   N_LAS_Vector * tmpStoVectorPtr,
   N_LAS_Vector * tmpCurrStoVectorPtr,
   N_LAS_Vector * tmpLastStoVectorPtr,
   N_LAS_Vector * tmpStoLeadCurrQCompVectorPtr,
   N_LAS_Vector * tmpQVectorPtr,
   N_LAS_Vector * tmpFVectorPtr,
   N_LAS_Vector * tmpdFdxdVpVectorPtr,
   N_LAS_Vector * tmpdQdxdVpVectorPtr)
{
  bool bsuccess = true;
  bool tmpBool = true;

#ifdef Xyce_PARALLEL_MPI
  N_ERH_ErrorMgr::startSafeBarrier();
#endif

  // copy over the passed pointers:
  extData_.nextSolVectorPtr = tmpSolVectorPtr;
  extData_.currSolVectorPtr = tmpCurrSolVectorPtr;
  extData_.lastSolVectorPtr = tmpLastSolVectorPtr;
  extData_.daeQVectorPtr    = tmpQVectorPtr;
  extData_.daeFVectorPtr    = tmpFVectorPtr;
  extData_.dFdxdVpVectorPtr = tmpdFdxdVpVectorPtr;
  extData_.dQdxdVpVectorPtr = tmpdQdxdVpVectorPtr;
  extData_.nextStaVectorPtr = tmpStaVectorPtr;
  extData_.currStaVectorPtr = tmpCurrStaVectorPtr;
  extData_.lastStaVectorPtr = tmpLastStaVectorPtr;
  extData_.storeLeadCurrQCompPtr = tmpStoLeadCurrQCompVectorPtr;
  extData_.nextStaDerivVectorPtr = tmpStaDerivVectorPtr;
  extData_.nextStoVectorPtr = tmpStoVectorPtr;
  extData_.currStoVectorPtr = tmpCurrStoVectorPtr;
  extData_.lastStoVectorPtr = tmpLastStoVectorPtr;

  // Make sure all boundary data is valid in the solution vector
#ifdef Xyce_PARALLEL_MPI
  extData_.nextSolVectorPtr->importOverlap();
  extData_.nextStaDerivVectorPtr->importOverlap();
#endif

  // Set up the relevant RAW Pointers:
  setupRawVectorPointers_ ();

#ifdef Xyce_OLD_LOOPING
  // call all the intermediate vars loads:
  vector<DeviceInstance*>::iterator iter;
  vector<DeviceInstance*>::iterator begin;
  vector<DeviceInstance*>::iterator end;

  // if inner problem, only do the PDE loads.
  if (solState_.twoLevelNewtonCouplingMode==INNER_PROBLEM)
  {
    begin = pdeInstancePtrVec_.begin ();
    end   = pdeInstancePtrVec_.end ();
  }
  else
  {
    begin = instancePtrVec_.begin ();
    end   = instancePtrVec_.end ();
  }

  //updateDependentParameters_();

  for (iter=begin; iter!=end;++iter)
  {
    tmpBool = (*iter)->updateSecondaryState();
    bsuccess = bsuccess && tmpBool;
  }
#else


  // call all the intermediate vars loads:
  vector<Device*>::iterator iter;
  vector<Device*>::iterator begin;
  vector<Device*>::iterator end;

  // if inner problem, only do the PDE loads.
  if (solState_.twoLevelNewtonCouplingMode==INNER_PROBLEM)
  {
    begin = pdeDevicePtrVec_.begin ();
    end   = pdeDevicePtrVec_.end ();
  }
  else
  {
    begin = devicePtrVec_.begin ();
    end   = devicePtrVec_.end ();
  }

  //updateDependentParameters_();

#ifndef Xyce_EXCLUDE_SECONDARY_STATE
  for (iter=begin; iter!=end;++iter)
  {
    tmpBool = (*iter)->updateSecondaryState (extData_.nextStaDerivVectorRawPtr, extData_.nextStoVectorRawPtr);
    bsuccess = bsuccess && tmpBool;
  }
#endif // Xyce_EXCLUDE_SECONDARY_STATE

#endif

#ifdef Xyce_PARALLEL_MPI
  extData_.nextStaVectorPtr->importOverlap();
  extData_.nextStoVectorPtr->importOverlap();
#endif


#ifdef Xyce_OLD_LOOPING
  for (iter=begin; iter!=end;++iter)
  {
    tmpBool=(*iter)->loadDAEQVector();bsuccess = bsuccess && tmpBool;
    tmpBool=(*iter)->loadDAEFVector();bsuccess = bsuccess && tmpBool;
  }
#else

  if (solState_.twoLevelNewtonCouplingMode==INNER_PROBLEM)
  {
    int numDevices = pdeDevicePtrVec_.size();
    for(int i=0; i< numDevices; ++i)
    {
      bsuccess=pdeDevicePtrVec_.at(i)->loadDAEVectors( extData_.nextSolVectorRawPtr,
                                                       extData_.daeFVectorRawPtr,
                                                       extData_.daeQVectorRawPtr,
                                                       extData_.nextStoVectorRawPtr,
                                                       extData_.storeLeadCurrQCompRawPtr);
    }
  }
  else
  {
    int numDevices = devicePtrVec_.size();
    for(int i=0; i< numDevices; ++i)
    {
      bsuccess=devicePtrVec_.at(i)->loadDAEVectors( extData_.nextSolVectorRawPtr,
                                                    extData_.daeFVectorRawPtr,
                                                    extData_.daeQVectorRawPtr,
                                                    extData_.nextStoVectorRawPtr,
                                                    extData_.storeLeadCurrQCompRawPtr);
      // std::cout << "Device i = " << i << "  " << devicePtrVec_.at(i)->name << std::endl;
      // std::cout.precision(15);
      // std::cout << "F = " << std::endl;
      // extData_.daeFVectorPtr->printPetraObject();
      // std::cout << "Q = " << std::endl;
      // extData_.daeQVectorPtr->printPetraObject();
    }
  }

#endif

  // dump to the screen:
#ifdef Xyce_DEBUG_DEVICE
  // note: this should eventually go away!
  if (devOptions_.debugLevel > 1 && solState_.debugTimeFlag)
  {
    int newtonIter = solState_.newtonIter;
    std::cout <<  "Q-vector: nonlinear iteration = " << newtonIter << "\n";
    extData_.daeQVectorPtr->printPetraObject();
    std::cout << std::endl;
    std::cout <<  "F-vector: nonlinear iteration = " << newtonIter << "\n";
    extData_.daeFVectorPtr->printPetraObject();
    std::cout << std::endl;

    if (devOptions_.voltageLimiterFlag)
    {
      std::cout << "\n\n  dFdxdVp vector: nonlinear iteration = " << newtonIter << "\n";
      extData_.dFdxdVpVectorPtr->printPetraObject();
      std::cout << std::endl;
      std::cout << "\n\n  dQdxdVp vector: nonlinear iteration = " << newtonIter << "\n";
      extData_.dQdxdVpVectorPtr->printPetraObject();
      std::cout << std::endl;
    }
  }

  debugOutput2();
#endif

  // Update parallel if necessary
  extData_.daeQVectorPtr->fillComplete();
  extData_.daeFVectorPtr->fillComplete();
  extData_.dFdxdVpVectorPtr->fillComplete();
  extData_.dQdxdVpVectorPtr->fillComplete();

#ifdef Xyce_PARALLEL_MPI
  N_ERH_ErrorMgr::safeBarrier(0);
#endif

#ifdef Xyce_SIZEOF
  int size = sizeof(*this);
  std::cout << "Size of device package after vector load = " << size << std::endl;
#endif

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::loadDeviceMask ()
// Purpose       : let devices set elements of a mask telling the time
//                 integrator what equations should be ignored in taking
//                 weighted norms for error control purposes.
// Special Notes : Devices should *only* zero internal variables, and then
//                 only those that absolutely should never be used to
//                 control step size (e.g. excess phase variables in BJTs)
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 01/18/07
//-----------------------------------------------------------------------------
bool DeviceMgr::loadDeviceMask()
{

  // call all the intermediate vars loads:
  vector<DeviceInstance*>::iterator iter;
  vector<DeviceInstance*>::iterator begin;
  vector<DeviceInstance*>::iterator end;

  begin = instancePtrVec_.begin ();
  end   = instancePtrVec_.end ();

  nonTrivialDeviceMaskFlag = false;

  for (iter=begin; iter!=end;++iter)
  {
    nonTrivialDeviceMaskFlag |= (*iter)->loadDeviceMask();
  }

  extData_.deviceMaskVectorPtr->fillComplete();

  // make sure the system's flag reflects ours:
  extData_.lasSysPtr->setNonTrivialDeviceMaskFlag(nonTrivialDeviceMaskFlag);

  return nonTrivialDeviceMaskFlag;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::addGlobalPar ()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 11/18/05
//-----------------------------------------------------------------------------
void DeviceMgr::addGlobalPar (N_UTL_Param & par)
{
  int pos;
  if (par.getType() == EXPR)
  {
    vector<string> variables, specials;
    vector<string>::iterator vs_i;
    double val;

    solState_.global_expressions.push_back(par.eVal());
    solState_.global_exp_names.push_back(par.uTag());
    pos = solState_.global_expressions.size()-1;
    N_UTL_Expression *e = &(solState_.global_expressions[pos]);
    e->get_names(XEXP_VARIABLE, variables);
    e->get_names(XEXP_SPECIAL, specials);
    if (!specials.empty())
    {
      e->set_sim_time(solState_.currTime);
    }
    if (!variables.empty())
    {
      for (vs_i=variables.begin() ; vs_i!=variables.end() ; ++vs_i)
        e->set_var(*vs_i,solState_.global_params[*vs_i]);
    }
    e->evaluateFunction (val);
    solState_.global_params[par.uTag()] = val;
  }
  else
  {
    solState_.global_params[par.uTag()] = par.dVal();
  }
}


//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getGlobalPar ()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schie, SNL, Electrical Systems Modeling
// Creation Date : 01/25/13
//-----------------------------------------------------------------------------
double DeviceMgr::getGlobalPar (const string & parName) const
{
  double retVal = 0;
  map<string,double>::const_iterator parLocItr = solState_.global_params.find(parName);
  if( parLocItr != solState_.global_params.end() )
  {
    // extract the value for return
    retVal = parLocItr->second;
  }
  else
  {
    string msg = "Could not find global parameter \"" + parName + "\"";
    N_ERH_ErrorMgr::report (N_ERH_ErrorMgr::USR_FATAL_0, msg);
  }
  return retVal;
}


//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setUpPDEDeviceFlagArray ()
// Purpose       : This function is intended to set which device types
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/18/02
//-----------------------------------------------------------------------------

bool DeviceMgr::setUpPDEDeviceFlagArray_ ()
{
  PDEDeviceFlag_[ModelType::DIODE_PDE]     = 1;
  PDEDeviceFlag_[ModelType::TWO_D_PDE]     = 1;

#ifdef Xyce_RAD_MODELS
  ExtendedModelLevels::setUpPDEDeviceFlagArray_ (PDEDeviceFlag_);
#endif

  return true;
}

#if 0
//-----------------------------------------------------------------------------
// Function      : DeviceMgr::loadPerturbationVector_ ()
// Purpose       : This function sets up the perturbation vector which
//                 is  needed in for the numerical Jacobian calculation.
//
// Special Notes : The formula is copied over from MPSalsa's num_jac_delta
//                 function in the file rf_fill_num_jac.c.
//
//                 Note that the bracket operator uses *local* indexing.
//
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/21/00
//-----------------------------------------------------------------------------

bool DeviceMgr::loadPerturbationVector_()
{
  int i;
  int isize = lasSysPtr_->getSolutionSize();
  for (i = 0; i < isize; ++i)
    (*(extData_.perturbVectorPtr))[i]  =
      solState_.numJacSqrtEta *
        (1.0 + fabs( (*(extData_.nextSolVectorPtr))[i]));

  return true;
}
#endif

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::updateIntermediateVars_
// Purpose       : This function calls updateIntermediateVars
//                 for the current devices.
// Special Notes :
// Scope         : private
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 12/14/2006
//-----------------------------------------------------------------------------
bool DeviceMgr::updateIntermediateVars_()
{
  bool bsuccess = true;

  vector<DeviceInstance*>::iterator iter;
  vector<DeviceInstance*>::iterator begin =instancePtrVec_.begin();
  vector<DeviceInstance*>::iterator end =instancePtrVec_.end();
  for (iter=begin; iter!=end;++iter)
  {
    (*iter)->updateIntermediateVars ();
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::updatePrimaryState_
// Purpose       : This function updates primary states
//                 for the present time step.
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/21/00
//-----------------------------------------------------------------------------
bool DeviceMgr::updatePrimaryState_()
{
  bool bsuccess = true;

  vector<DeviceInstance*>::iterator iter;
  vector<DeviceInstance*>::iterator begin =instancePtrVec_.begin();
  vector<DeviceInstance*>::iterator end =instancePtrVec_.end();
  for (iter=begin; iter!=end;++iter)
  {
    (*iter)->updatePrimaryState ();
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::updateSecondaryState_
// Purpose       : This function function updates secondary states
//                 for the present time step.
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/21/00
//-----------------------------------------------------------------------------
bool DeviceMgr::updateSecondaryState_()
{
  bool bsuccess = true;

  vector<DeviceInstance*>::iterator iter;
  vector<DeviceInstance*>::iterator begin =instancePtrVec_.begin();
  vector<DeviceInstance*>::iterator end =instancePtrVec_.end();
  for (iter=begin; iter!=end;++iter)
  {
    (*iter)->updateSecondaryState ();
  }

  return bsuccess;
}

//----------------------------------------------------------------------------
// Function       : DeviceMgr::updateDependentParameters_
// Purpose        : This function updates all dependent parameters for
//                  the current time step.
// Special Notes  : This was evolved from updateTimeDependentParameters_
// Scope          : private
// Creator        : Dave Shirley
// Creation Date  : 08/17/06
//----------------------------------------------------------------------------
bool DeviceMgr::updateDependentParameters_()
{
  bool bsuccess = true;
  bool tmpBool = true;
  N_LAS_Vector * solVectorPtr = extData_.nextSolVectorPtr;

  map<string,double> & gp = solState_.global_params;
  vector<N_UTL_Expression> & ge = solState_.global_expressions;

  if (timeParamsProcessed_ != solState_.currTime)
    parameterChanged_ = true;
  if (!ge.empty())
  {
    // Update global params for new time and other global params
    vector<string> variables;
    vector<string>::iterator vs_i;
    vector<string>::iterator vs_end;
    double val;
    bool changed;

    int pos = 0;
    vector<N_UTL_Expression>::iterator g_i = ge.begin();
    vector<N_UTL_Expression>::iterator g_end = ge.end();
    for ( ; g_i != g_end; ++g_i)
    {
      changed = false;
      g_i->get_names(XEXP_VARIABLE, variables);
      if (g_i->set_sim_time(solState_.currTime))
        changed = true;
      if (!variables.empty())
      {
        vs_i=variables.begin();
        vs_end=variables.end();
        for ( ; vs_i!=vs_end; ++vs_i)
        {
          if (g_i->set_var(*vs_i,gp[*vs_i]))
            changed = true;
        }
      }
      if (changed)
      {
        parameterChanged_ = true;
        g_i->evaluateFunction (val);
        gp[solState_.global_exp_names[pos]] = val;
      }
      ++pos;
    }
  }

  // do the models:
  if (firstDependent)
  {
    dependentPtrVec_.clear();
    firstDependent = false;

    vector<DeviceModel*>::iterator iterM;
    vector<DeviceModel*>::iterator beginM =modelPtrVec_.begin();
    vector<DeviceModel*>::iterator endM =modelPtrVec_.end();
    for (iterM=beginM; iterM!=endM;++iterM)
    {
      if (!(*iterM)->getDependentParams().empty())
      {
        dependentPtrVec_.push_back(static_cast<DeviceEntity *>(*iterM));
        tmpBool = (*iterM)->updateGlobalParameters(gp);
        bsuccess = bsuccess && tmpBool;
        tmpBool = (*iterM)->updateDependentParameters (*solVectorPtr);
        bsuccess = bsuccess && tmpBool;
        (*iterM)->processParams();
        (*iterM)->processInstanceParams();
      }
    }

    // do the instances
    vector<DeviceInstance*>::iterator iter;
    vector<DeviceInstance*>::iterator begin =instancePtrVec_.begin();
    vector<DeviceInstance*>::iterator end =instancePtrVec_.end();
    for (iter=begin; iter!=end;++iter)
    {
      if (!(*iter)->getDependentParams().empty())
      {
        dependentPtrVec_.push_back(static_cast<DeviceEntity *>(*iter));
        tmpBool = (*iter)->updateGlobalParameters (gp);
        bsuccess = bsuccess && tmpBool;
        tmpBool = (*iter)->updateDependentParameters (*solVectorPtr);
        bsuccess = bsuccess && tmpBool;
        (*iter)->processParams();
      }
    }
  }
  else
  {
    bool changed;
    vector<DeviceEntity*>::iterator iter;
    vector<DeviceEntity*>::iterator begin = dependentPtrVec_.begin();
    vector<DeviceEntity*>::iterator end = dependentPtrVec_.end();
    for (iter=begin; iter!=end;++iter)
    {
      changed = false;
      if (parameterChanged_)
      {
        tmpBool = (*iter)->updateGlobalParameters (gp);
        changed = changed || tmpBool;
        bsuccess = bsuccess && tmpBool;
      }
      tmpBool = (*iter)->updateDependentParameters (*solVectorPtr);
      changed = changed || tmpBool;
      bsuccess = bsuccess && tmpBool;
      if (changed)
      {
        (*iter)->processParams();
        (*iter)->processInstanceParams();
      }
    }
  }
  timeParamsProcessed_ = solState_.currTime;
  parameterChanged_ = false;

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : DeviceMgr::loadBVectorsforAC
// Purpose       : This function loads the B-vector contributions for sources.
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool DeviceMgr::loadBVectorsforAC(N_LAS_Vector * bVecRealPtr,
                                        N_LAS_Vector * bVecImagPtr)
{
  bool bsuccess = true;

//  setupSolverInfo_ ();

// copy over the passed pointers:
  extData_.bVecRealPtr = bVecRealPtr;
  extData_.bVecImagPtr = bVecImagPtr;

#ifdef Xyce_PARALLEL_MPI
  extData_.nextSolVectorPtr->importOverlap();
#endif
  // Now reset the relevant RAW pointers:
  extData_.bVecRealRawPtr = &((*extData_.bVecRealPtr)[0]);
  extData_.bVecImagRawPtr = &((*extData_.bVecImagPtr)[0]);

  vector<SourceInstance*>::iterator vIter;
  vector<SourceInstance*>::iterator vBegin =indepSourceInstancePtrVec_.begin();
  vector<SourceInstance*>::iterator vEnd =indepSourceInstancePtrVec_.end();
  for (vIter=vBegin; vIter!=vEnd;++vIter)
  {
    (*vIter)->loadBVectorsforAC(extData_.bVecRealRawPtr, extData_.bVecImagRawPtr);
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getBMatrixEntriesforMOR()
// Purpose       : This function obtains the indices for the B-vector contributions for sources.
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist and Ting Mei, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool DeviceMgr::getBMatrixEntriesforMOR(std::vector<int>& bMatEntriesVec,
                                              std::vector<int>& bMatPosEntriesVec)
{
  bool bsuccess = true;

  int lpos, lneg, lbra;
  vector<N_DEV_SourceInstance*>::iterator vIter;
  vector<N_DEV_SourceInstance*>::iterator vBegin =indepSourceInstancePtrVec_.begin();
  vector<N_DEV_SourceInstance*>::iterator vEnd =indepSourceInstancePtrVec_.end();

  for (vIter=vBegin; vIter!=vEnd;++vIter)
  {
     N_DEV_VsrcInstance* vsrc = dynamic_cast<N_DEV_VsrcInstance *>(*vIter);
     if (vsrc != 0)
     {
       vsrc->getLIDs(lpos, lneg, lbra);
       bMatEntriesVec.push_back( lbra );
       bMatPosEntriesVec.push_back( lpos );
     }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::updateSources
// Purpose       : This function function updates sources for the present
//                 time step.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/26/00
//-----------------------------------------------------------------------------
bool DeviceMgr::updateSources()
{
  bool bsuccess = true;

  setupSolverInfo_ ();

  vector<SourceInstance*>::iterator vIter;
  vector<SourceInstance*>::iterator vBegin =indepSourceInstancePtrVec_.begin();
  vector<SourceInstance*>::iterator vEnd =indepSourceInstancePtrVec_.end();
  for (vIter=vBegin; vIter!=vEnd;++vIter)
  {
    (*vIter)->updateSource();
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setICs
// Purpose       : This function function sets initial conditions for devices
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 9/13/00
//-----------------------------------------------------------------------------
bool DeviceMgr::setICs(
   N_LAS_Vector * tmpSolVectorPtr,
   N_LAS_Vector * tmpCurrSolVectorPtr,
   N_LAS_Vector * tmpLastSolVectorPtr,
   N_LAS_Vector * tmpStaVectorPtr,
   N_LAS_Vector * tmpCurrStaVectorPtr,
   N_LAS_Vector * tmpLastStaVectorPtr,
   N_LAS_Vector * tmpStaDerivVectorPtr,
   N_LAS_Vector * tmpStoVectorPtr,
   N_LAS_Vector * tmpCurrStoVectorPtr,
   N_LAS_Vector * tmpLastStoVectorPtr,
   N_LAS_Vector * tmpQVectorPtr,
   N_LAS_Vector * tmpFVectorPtr,
   N_LAS_Vector * tmpdFdxdVpVectorPtr,
   N_LAS_Vector * tmpdQdxdVpVectorPtr)
{
  bool bsuccess = true;

#ifdef Xyce_PARALLEL_MPI
  N_ERH_ErrorMgr::startSafeBarrier();
#endif

  // copy over the passed pointers:
  extData_.nextSolVectorPtr = tmpSolVectorPtr;
  extData_.currSolVectorPtr = tmpCurrSolVectorPtr;
  extData_.lastSolVectorPtr = tmpLastSolVectorPtr;
  extData_.daeQVectorPtr    = tmpQVectorPtr;
  extData_.daeFVectorPtr    = tmpFVectorPtr;
  extData_.dFdxdVpVectorPtr = tmpdFdxdVpVectorPtr;
  extData_.dQdxdVpVectorPtr = tmpdQdxdVpVectorPtr;
  extData_.nextStaVectorPtr = tmpStaVectorPtr;
  extData_.currStaVectorPtr = tmpCurrStaVectorPtr;
  extData_.lastStaVectorPtr = tmpLastStaVectorPtr;
  extData_.nextStaDerivVectorPtr = tmpStaDerivVectorPtr;
  extData_.nextStoVectorPtr = tmpStoVectorPtr;
  extData_.currStoVectorPtr = tmpCurrStoVectorPtr;
  extData_.lastStoVectorPtr = tmpLastStoVectorPtr;

  // Make sure all boundary data is valid in the solution vector
#ifdef Xyce_PARALLEL_MPI
  extData_.nextSolVectorPtr->importOverlap();
  extData_.nextStaDerivVectorPtr->importOverlap();
#endif

  // if IC's on devices are set, we need to ensure that the
  // raw pointers are up to date first.
  setupRawVectorPointers_ ();

  vector<DeviceInstance*>::iterator iter;
  vector<DeviceInstance*>::iterator begin =instancePtrVec_.begin();
  vector<DeviceInstance*>::iterator end =instancePtrVec_.end();
  for (iter=begin; iter!=end;++iter)
  {
    (*iter)->setIC ();
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getLinearSystemFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/17/00
//-----------------------------------------------------------------------------
bool DeviceMgr::getLinearSystemFlag()
{
  return linearSystemFlag_;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getVoltageLimiterFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/11/04
//-----------------------------------------------------------------------------
bool DeviceMgr::getVoltageLimiterFlag ()
{
  return devOptions_.voltageLimiterFlag;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getPDESystemFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/25/01
//-----------------------------------------------------------------------------
bool DeviceMgr::getPDESystemFlag()
{
  return solState_.PDESystemFlag;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::runParameterTests
// Purpose       : The  purpose of this function is to test out the processing
//                 of various device parameters.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/14/01
//-----------------------------------------------------------------------------
bool DeviceMgr::runParameterTests(string & deviceName)
{
  bool bsuccess = true;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::output
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/13/01
//-----------------------------------------------------------------------------
bool DeviceMgr::output ()
{
  bool bsuccess = true;
  bool tmpBool = true;

  vector<DeviceInstance*>::iterator iter;
  vector<DeviceInstance*>::iterator begin = plotFileInstancePtrVec_.begin ();
  vector<DeviceInstance*>::iterator end = plotFileInstancePtrVec_.end ();
  for (iter=begin; iter!=end;++iter)
  {
    tmpBool = (*iter)->outputPlotFiles ();
    bsuccess = bsuccess && tmpBool;
  }

  // only output .OP information once.
  if ( !dotOpOutputFlag && solState_.OPspecified )
  {
    dotOpOutput();
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::finishOutput
// Purpose       : Same as output, only this one forces the output.
//
// Special Notes : This function was neccessary with the implementation of
//                 outputInterval.  The final output, which needs to
//                 happen after the transient is over, won't neccessarily
//                 happen with outputInterval enabled.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/19/04
//-----------------------------------------------------------------------------
bool DeviceMgr::finishOutput ()
{
  bool bsuccess = true;
  bool tmpBool = true;

  solState_.forceFinalOutput = true;
  tmpBool = output ();
  bsuccess = bsuccess && tmpBool;
  solState_.forceFinalOutput = false;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::dotOpOutput
// Purpose       :
// Special Notes : This is a quick-and-dirty implementation, to get something
//                 working quickly.
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 5/3/12
//-----------------------------------------------------------------------------
void  DeviceMgr::dotOpOutput ()
{
  dotOpOutputFlag = true;

  const string dashedline("--------------------------------------------------"
    "---------------------------");
  //std::cout << dashedline << std::endl;
  string msg = dashedline + "\nOperating point information:";
  N_ERH_ErrorMgr::report (N_ERH_ErrorMgr::USR_INFO_0, msg);
  int i;
  for (i = 1; i < ModelType::NUMDEV; ++i)
  {
    if (deviceUseFlag_[i] == 1) deviceArray_[i]->printDotOpOutput(std::cout);
  }
  msg = dashedline;
  N_ERH_ErrorMgr::report (N_ERH_ErrorMgr::USR_INFO_0, msg);
  //outputParams (int mode, map <int, DeviceEntity *> & base )
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setGlobalFlags
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 11/2/12
//-----------------------------------------------------------------------------
void DeviceMgr::setGlobalFlags()
{
#ifdef Xyce_PARALLEL_MPI
  // just in case this has not been synchronized in parallel
  double pdeSys_glob=0.0;
  double pdeSys_par = solState_.PDESystemFlag?1.0:0.0;

  N_PDS_Comm *pdsCommPtr = pdsMgrPtr_->getPDSComm();
  pdsCommPtr->barrier();
  pdsCommPtr->sumAll(&pdeSys_par, &pdeSys_glob, 1);
  solState_.PDESystemFlag = (pdeSys_glob != 0.0)?true:false;

#if 0
  if (solState_.PDESystemFlag)
    std::cout << "PDESystemFlag = TRUE"<<std::endl;
  else
    std::cout << "PDESystemFlag = FALSE"<<std::endl;
#endif
#endif

  return;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getBreakPoints
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/08/01
//-----------------------------------------------------------------------------
bool DeviceMgr::getBreakPoints ( vector<N_UTL_BreakPoint> & breakPointTimes )
{
  vector<DeviceInstance*>::iterator iterI;
  vector<DeviceModel*>::iterator iterM;
  bool bsuccess = true;
  bool tmpBool = true;

  tmpBool = setupSolverInfo_();
  bsuccess = bsuccess && tmpBool;
  setupRawVectorPointers_ ();

  // For some devices we only need to set breakpoints caused by discontinuities
  // in their parameters:

  for (iterM = modelPtrVec_.begin() ; iterM != modelPtrVec_.end() ; ++iterM)
  {
    if (!(*iterM)->getDependentParams().empty())
    {
      tmpBool = (*iterM)->getParamBreakpoints(breakPointTimes);
      bsuccess = bsuccess && tmpBool;
    }
  }

  for (iterI = instancePtrVec_.begin() ; iterI != instancePtrVec_.end() ; ++iterI)
  {
    if (!(*iterI)->getDependentParams().empty())
    {
      tmpBool = (*iterI)->getParamBreakpoints(breakPointTimes);
      bsuccess = bsuccess && tmpBool;
    }
  }

  // Breakpoints for global params:
  vector<N_UTL_Expression>::iterator globalExp_i =
    solState_.global_expressions.begin();
  vector<N_UTL_Expression>::iterator globalExp_end =
    solState_.global_expressions.end();
  for ( ; globalExp_i != globalExp_end; ++globalExp_i)
  {
    double bTime = globalExp_i->get_break_time();
    if (bTime > solState_.currTime)
      breakPointTimes.push_back(bTime);
  }

  if (!breakPointInstancesInitialized)
  {
    vector<DeviceInstance*>::iterator beginI = instancePtrVec_.begin ();
    vector<DeviceInstance*>::iterator endI = instancePtrVec_.end ();

    for (iterI=beginI;iterI!=endI;++iterI)
    {
      // this function returns false if it is the base class, and true otherwise.
      bool functionSetup = (*iterI)->getInstanceBreakPoints ( breakPointTimes );
      if (functionSetup)
      {
        bpInstancePtrVec_.push_back(*iterI);
      }
    }
    breakPointInstancesInitialized = true;
  }
  else
  {
    vector<DeviceInstance*>::iterator beginI = bpInstancePtrVec_.begin ();
    vector<DeviceInstance*>::iterator endI = bpInstancePtrVec_.end ();
    for (iterI=beginI;iterI!=endI;++iterI)
    {
      bool functionSetup = (*iterI)->getInstanceBreakPoints ( breakPointTimes );
    }
  }

#ifdef Xyce_EXTDEV
  vector<N_DEV_ExternDeviceInstance*>::iterator iter;
  vector<N_DEV_ExternDeviceInstance*>::iterator begin;
  vector<N_DEV_ExternDeviceInstance*>::iterator end;
  begin = extDevInstancePtrVec_.begin ();
  end = extDevInstancePtrVec_.end ();
  for (iter=begin; iter!=end;++iter)
  {
    (*iter)->getBreakPoints ( breakPointTimes );
  }
#endif

  return bsuccess;

}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setupSolverInfo_
//
// Purpose       :
//
// Special Notes : This function gets called a lot, and this can be kind of
//                 confusing.  Probably, this function is being used to handle
//                 too many different types of data.
//
//                 For example, it gets called at the beginning
//                 of the "updateSources" function.  Why?  Because the sources
//                 need to know which step we are at, and/or what the current
//                 time is, to do their update properly.
//
//                 However, at the moment, updateSources also provides
//                 information that setupSolverInfo needs.  Only after the
//                 sources have been updated do we know if a sweep source has
//                 been reset.  And, the sweepSourceResetFlag  is used by
//                 setupSolverInfo, to set up the initJctFlag boolean.
//                 So it will need to be called at least one more
//                 time, at the beginning of the RHS load.  (which it is).
//
//                 Anyway, any of the functions that are called from the
//                 outside, such as:  updateSources, loadRHSVector,
//                 loadJacobianMatrix, etc.... have no way of knowing when,
//                 w.r.t. the solvers they are being called.  The only way
//                 to do this properly is to have each of them request
//                 the current solver state, before they go to do their
//                 work.  Hence, this function is called a lot.
//
//                 Unfortunately, this has led to a somewhat sloppy and
//                 confusing interface between the solvers and the
//                 device package.  I wanted to avoid having a lot of
//                 function arguments being passed around for each of
//                 these functions, in part because the calling code
//                 (NLS) doesn't know everything. - NLS knows about the newton
//                 step, but it doesn't know the time step, for example.
//
//                 At some point I hope to refactor this.
//
// Scope         : private
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/30/01
//-----------------------------------------------------------------------------
bool DeviceMgr::setupSolverInfo_ ()
{
  bool bsuccess = true;

  N_TIA_TimeIntInfo & tiInfo = solState_.tiInfo;
  N_NLS_NonLinInfo  & nlInfo = solState_.nlInfo;

  // Time integration info:
  anaIntPtr_->getTimeIntInfo(tiInfo);
  solState_.pdt                  = tiInfo.pdt;
  solState_.currTimeStep         = tiInfo.nextTimeStep;
  solState_.lastTimeStep         = tiInfo.currTimeStep;
  solState_.currTime             = tiInfo.nextTime;
  solState_.finalTime            = tiInfo.finalTime;
  solState_.startingTimeStep     = tiInfo.startingTimeStep;
  solState_.bpTol                = tiInfo.bpTol;
  solState_.currentOrder         = tiInfo.currentOrder; // BNB
  solState_.usedOrder            = tiInfo.usedOrder;  // BNB

  if (solState_.mpdeOnFlag == 1)
  {
    solState_.dcopFlag = 0;
    solState_.initTranFlag = 0;
    solState_.beginIntegrationFlag = 0;
  }
  else
  {
    solState_.dcopFlag             = tiInfo.dcopFlag;
    solState_.initTranFlag         = tiInfo.initTranFlag;
    solState_.beginIntegrationFlag = tiInfo.beginIntegrationFlag;
  }

  solState_.inputOPFlag          = tiInfo.inputOPFlag;
  solState_.acopFlag             = tiInfo.acopFlag;
  solState_.tranopFlag           = tiInfo.tranopFlag;
  solState_.transientFlag        = tiInfo.transientFlag;
  solState_.dcsweepFlag          = tiInfo.dcsweepFlag;
  solState_.sweepSourceResetFlag = tiInfo.sweepSourceResetFlag;

  solState_.timeStepNumber       = tiInfo.timeStepNumber;
//  solState_.initTranFlag         = tiInfo.initTranFlag;
//  solState_.beginIntegrationFlag = tiInfo.beginIntegrationFlag;

  solState_.doubleDCOPStep       = tiInfo.doubleDCOPStep;
  solState_.doubleDCOPEnabled    = tiInfo.doubleDCOPEnabled;
  solState_.stepLoopIter         = tiInfo.stepLoopIter;

  // Nonlinear solver info:
  nlsMgrPtr_->getNonLinInfo(nlInfo);
  solState_.newtonIter           = nlInfo.newtonIter;
  solState_.twoLevelNewtonCouplingMode         = nlInfo.twoLevelNewtonCouplingMode;

  // Get LOCA-specific information.  Note - in general, LOCA is only used for
  // steady state calculations - there isn't much point in using it
  // for transient.  A common case is one where LOCA is used for the
  // tranOP, but not for the subsequent transient phase.  The
  // locaEnabledFlag should switch from true to false under
  // that scenario, once the transient phase starts.
  solState_.locaEnabledFlag      = nlInfo.locaFlag;
  if (solState_.locaEnabledFlag) // if no LOCA, these are 0, true, respectively.
  {
    solState_.continuationStepNumber = nlInfo.continuationStep;
    solState_.firstContinuationParam = nlInfo.firstContinuationParam;
    solState_.firstSolveComplete     = nlInfo.firstSolveComplete;
  }
  // Done with LOCA information.

  // Setup the initialize junctions flag.
  // The initJct flag should only be true if we are at the first Newton step of
  // an initial point in the calculation.  Examples include:
  //     - 1st Newton step of the DCOP initialization for transient (tranOp)
  //     - 1st Newton step of first DC sweep step.
  //     - 1st Newton step of a sweep that has been reset.  That typically
  //         happens if the sweep is multi-dimensional, and the inner loop has
  //         cycled back to the beginning again.
  //
  bool resetFlag =  (solState_.timeStepNumber==0) || (solState_.sweepSourceResetFlag);

  // Do this if using LOCA for a DC or tranop calculation.
  if (solState_.dcopFlag && solState_.locaEnabledFlag)
  {
    resetFlag = resetFlag && (solState_.continuationStepNumber==0);
  }

  solState_.initJctFlag = ( (solState_.dcopFlag) &&
                            (solState_.newtonIter==0) &&
                             solState_.firstContinuationParam &&
                             !solState_.firstSolveComplete && resetFlag );

  // initFixFlag: try to mimic "MODEINITFIX" of SPICE.  This is set if:
  //   DCOP or TranOP
  //   Not first iteration
  //   Any device not converged

  solState_.initFixFlag = ( (solState_.dcopFlag) &&
                            !(allDevsConverged()) &&
                            (solState_.newtonIter!=0) &&
                             solState_.firstContinuationParam &&
                             !solState_.firstSolveComplete && resetFlag );

  if ( solState_.dcopFlag )
  {
    solState_.ltraTimeIndex = 0;
    solState_.ltraTimeHistorySize = 10;
    solState_.ltraTimePoints.resize(solState_.ltraTimeHistorySize);
  }

  // One final check.  See if the "external state" has been set.  If so,
  // check to see if it has the initJctFlag set.  If not, then we probably
  // shouldn't either.  The external state comes from a higher up level
  // in a multi-level newton solve.
  //
  // This should be made more detailed later.
  if (externalStateFlag_)
  {
    if (solState_.newtonIter==0 && solState_.dcopFlag)
    {
      solState_.initJctFlag = solStateExternal_.initJctFlag;
    }
  }

  // The first DCOP step of a "double DCOP" simulation is a special case,
  // in which the nonlinear poisson is solved in place of drift-diffusion
  // equations for the PDE devices.  For this initialization problem, the
  // circuit subproblem should not be included.
  if ((solState_.doubleDCOPEnabled) &&
      (solState_.dcopFlag) &&
      (solState_.doubleDCOPStep == 0))
  {
    solState_.twoLevelNewtonCouplingMode       = INNER_PROBLEM;
  }

  // if necessary, set up the names vector
  if (devOptions_.testJacobianFlag)
  {
    NodeNamePairMap & nodeNames = outputMgrPtr_->getAllNodes();
    int nodeNameSize = nodeNames.size();
    nameVec_.resize(nodeNameSize+1,"gnd");
    NodeNamePairMap::iterator mapI, mapEnd;
    mapEnd = nodeNames.end();
    mapI =  nodeNames.begin();
    for ( ; mapI != mapEnd ; ++mapI)
    {
      nameVec_[(*mapI).second.first] = (*mapI).first;
    }
  }

#ifdef Xyce_DEBUG_DEVICE
  if (devOptions_.debugLevel > 0 && solState_.debugTimeFlag)
  {
    std::cout << solState_;
  }

  solState_.debugTimeFlag =
  (solState_.currTime       >= devOptions_.debugMinTime &&
   solState_.currTime       <= devOptions_.debugMaxTime) &&
  (solState_.timeStepNumber >= devOptions_.debugMinTimestep &&
   solState_.timeStepNumber <= devOptions_.debugMaxTimestep);

#endif  // Xyce_DEBUG_DEVICE

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setupRawVectorPointers_
// Purpose       : set up raw pointers
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/25/08
//-----------------------------------------------------------------------------
bool DeviceMgr::setupRawVectorPointers_ ()
{
  if (extData_.daeQVectorPtr != 0)
  {
    extData_.daeQVectorRawPtr    = &((*extData_.daeQVectorPtr)[0]);
  }

  if (extData_.daeFVectorPtr != 0)
  {
    extData_.daeFVectorRawPtr    = &((*extData_.daeFVectorPtr)[0]);
  }

  if (extData_.dFdxdVpVectorPtr != 0)
  {
    extData_.dFdxdVpVectorRawPtr = &((*extData_.dFdxdVpVectorPtr)[0]);
  }

  if (extData_.dQdxdVpVectorPtr != 0)
  {
    extData_.dQdxdVpVectorRawPtr = &((*extData_.dQdxdVpVectorPtr)[0]);
  }

  if ( extData_.nextSolVectorPtr != 0 )
  {
    extData_.nextSolVectorRawPtr = &((*extData_.nextSolVectorPtr)[0]);
  }

  if ( extData_.currSolVectorPtr != 0 )
  {
    extData_.currSolVectorRawPtr = &((*extData_.currSolVectorPtr)[0]);
  }

  if ( extData_.lastSolVectorPtr != 0 )
  {
    extData_.lastSolVectorRawPtr = &((*extData_.lastSolVectorPtr)[0]);
  }

  if ( extData_.nextStaVectorPtr != 0 )
  {
    extData_.nextStaVectorRawPtr = &((*extData_.nextStaVectorPtr)[0]);
  }

  if ( extData_.currStaVectorPtr != 0 )
  {
    extData_.currStaVectorRawPtr = &((*extData_.currStaVectorPtr)[0]);
  }

  if ( extData_.lastStaVectorPtr != 0 )
  {
    extData_.lastStaVectorRawPtr = &((*extData_.lastStaVectorPtr)[0]);
  }

  if ( extData_.nextStaDerivVectorPtr != 0 )
  {
    extData_.nextStaDerivVectorRawPtr = &((*extData_.nextStaDerivVectorPtr)[0]);
  }

  if ( extData_.nextStoVectorPtr != 0 )
  {
    extData_.nextStoVectorRawPtr = &((*extData_.nextStoVectorPtr)[0]);
  }

  if ( extData_.currStoVectorPtr != 0 )
  {
    extData_.currStoVectorRawPtr = &((*extData_.currStoVectorPtr)[0]);
  }

  if ( extData_.lastStoVectorPtr != 0 )
  {
    extData_.lastStoVectorRawPtr = &((*extData_.lastStoVectorPtr)[0]);
  }

  if ( extData_.storeLeadCurrQCompPtr != 0 )
  {
    extData_.storeLeadCurrQCompRawPtr = &((*extData_.storeLeadCurrQCompPtr)[0]);
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setupRawMatrixPointers_
// Purpose       : set up raw pointers for matrices
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, 1414
// Creation Date : 06/03/09
//-----------------------------------------------------------------------------
bool DeviceMgr::setupRawMatrixPointers_ ()
{
    vector<DeviceInstance*>::iterator iter;
    vector<DeviceInstance*>::iterator begin;
    vector<DeviceInstance*>::iterator end;
    begin = instancePtrVec_.begin();
    end   = instancePtrVec_.end();
    for (iter=begin; iter!=end;++iter)
    {
      (*iter)->setupPointers();
    }
    return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getMaxTimeStepSize
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/31/01
//-----------------------------------------------------------------------------
double DeviceMgr::getMaxTimeStepSize ()
{
  double maxStep = devOptions_.defaultMaxTimeStep;

  vector<DeviceInstance*>::iterator iter;
  vector<DeviceInstance*>::iterator begin =instancePtrVec_.begin();
  vector<DeviceInstance*>::iterator end =instancePtrVec_.end();
  for (iter=begin; iter!=end;++iter)
  {
    double step = (*iter)->getMaxTimeStepSize ();
    SourceInstance * srcInst = dynamic_cast<SourceInstance*>(*iter);
    if( !srcInst || !srcInst->getFastSourceFlag() )
      maxStep = Xycemin( step, maxStep );
  }

  return maxStep;
}


//-----------------------------------------------------------------------------
// Function      : DeviceMgr::declareCurrentStepAsBreakpoint
// Purpose       : If during a device load, a device must act in a discontinuous
//                 fashion, let the analysis manager know that this step should
//                 be treated as such.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 3/5/2010
//-----------------------------------------------------------------------------
void DeviceMgr::declareCurrentStepAsBreakpoint()
{
  anaIntPtr_->setBeginningIntegrationFlag(true);
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::enablePDEContinuation
// Purpose       : This function turns on the Continuation flag, which lets
//                 the devices which are Continuation enabled know that they
//                 need to set up their variable parameters.
//
// Special Notes : Currently, only PDE devices can take advantage of this
//                 capability, and it is only used in the context of two-level
//                 Newton.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 10/22/02
//-----------------------------------------------------------------------------
int DeviceMgr::enablePDEContinuation ()
{
#ifdef Xyce_DEBUG_DEVICE
  std::cout << "DeviceMgr::enablePDEContinuation" << std::endl;
#endif

  bool bsuccess = true;
  solState_.PDEcontinuationFlag  = true;
  solState_.maxPDEContinuationSteps = 1;
  solState_.currPDEContinuationStep = 0;

  vector<DeviceInstance*>::iterator iter;
  vector<DeviceInstance*>::iterator begin =instancePtrVec_.begin();
  vector<DeviceInstance*>::iterator end =instancePtrVec_.end();

  for (iter=begin; iter!=end;++iter)
  {
    bool tmpSuccess = (*iter)->enablePDEContinuation();
    bsuccess = bsuccess && tmpSuccess;
  }

  // if any of the devices feels that it needs more than the specified
  // number of continuation steps, re-do, with the new step number.
  if (solState_.maxPDEContinuationSteps != 1)
  {
    for (iter=begin; iter!=end;++iter)
    {
      bool tmpSuccess = (*iter)->enablePDEContinuation();
      bsuccess = bsuccess && tmpSuccess;
    }
  }

  int returnedSteps;

  if (!bsuccess)  returnedSteps = -1;
  else            returnedSteps = solState_.maxPDEContinuationSteps;

  return returnedSteps;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::disablePDEContinuation
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 10/22/02
//-----------------------------------------------------------------------------
bool DeviceMgr::disablePDEContinuation ()
{
  bool bsuccess = true;
  solState_.PDEcontinuationFlag = false;

  vector<DeviceInstance*>::iterator iter;
  vector<DeviceInstance*>::iterator begin =instancePtrVec_.begin();
  vector<DeviceInstance*>::iterator end =instancePtrVec_.end();
  for (iter=begin; iter!=end;++iter)
  {
    bool tmpSuccess = (*iter)->disablePDEContinuation();
    bsuccess = bsuccess && tmpSuccess;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::calcPDESubProblemInfo
//
// Purpose       : Determines the number of PDE sub-problems.,
//
//                 This is mainly used/needed for 2-level problems.
//
//                 Also determines the number of interface nodes per
//                 sub-problem (ie number of electrodes on each device).
//
// Special Notes : Need to modify to work correctly in parallel, probably.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 12/03/02
//-----------------------------------------------------------------------------
bool DeviceMgr::calcPDESubProblemInfo ()
{
  bool bsuccess = true;

  numPDEDevices_ = pdeInstancePtrVec_.size ();

  // now set up numInterfaceNodes_;
  numInterfaceNodes_.resize(numPDEDevices_);

  for (int i=0;i<numPDEDevices_;++i)
  {
    numInterfaceNodes_[i] = pdeInstancePtrVec_[i]->getNumExtVars ();
  }

  calledBeforeCSPI = true;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getNumInterfaceNodes
// Purpose       : returns the vector calculaed in calcPDESubProblemInfo.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 12/03/02
//-----------------------------------------------------------------------------
void DeviceMgr::getNumInterfaceNodes (vector<int> & numINodes)
{
  if (!calledBeforeCSPI)
  {
    calcPDESubProblemInfo ();
  }

  int size = numINodes.size ();
  int size2 = numInterfaceNodes_.size();

  if (size < size2 ) numINodes.resize(size2);

  for (int i=0;i<size2;++i)
  {
    numINodes[i] = numInterfaceNodes_[i];
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::loadCouplingRHS
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 12/03/02
//-----------------------------------------------------------------------------
bool DeviceMgr::loadCouplingRHS (int iPDEDevice, int iElectrode, N_LAS_Vector * dfdvPtr)
{
  return pdeInstancePtrVec_[iPDEDevice]->loadDFDV(iElectrode,dfdvPtr);
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::calcCouplingTerms
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 12/03/02
//-----------------------------------------------------------------------------
bool DeviceMgr::calcCouplingTerms (int iPDEDevice, int iElectrode, const N_LAS_Vector * dxdvPtr)
{
  return pdeInstancePtrVec_[iPDEDevice]->calcConductance(iElectrode, dxdvPtr);
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::raiseDebugLevel
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 04/23/03
//-----------------------------------------------------------------------------
bool DeviceMgr::raiseDebugLevel(int increment)
{
#ifdef Xyce_DEBUG_DEVICE
  devOptions_.debugLevel += increment;
#endif
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getDACDeviceNames
// Purpose       : Get a list of names of all DACs in the circuit
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and MOdels
// Creation Date : 05/05/04
//-----------------------------------------------------------------------------
bool DeviceMgr::getDACDeviceNames(vector < string > & dacNames)
{

  bool bSuccess = true;

  dacNames.clear();

  // find the device index for DAC devices:

  int dacIndex = getDeviceIndex("Y%DAC%DUMMYDAC");

  // get the singleton for the DAC device
  Device* devicePtr = returnDevicePtr(dacIndex);

  // If it wasn't found, there were no DAC devices
  if (devicePtr == NULL || devicePtr == dummy_ptr)
  {
    string msg("DeviceMgr::getDACDeviceNames:");
    msg += "  No singleton found, circuit does not apparently have any DAC devices in it.\n";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_WARNING_0, msg);
    // this is not a failure, though
  }
  else
  {
    N_DEV_DACMaster* DACPtr = dynamic_cast<N_DEV_DACMaster*>(devicePtr);
    bSuccess = DACPtr->getDACDeviceNames( dacNames );
  }

  return bSuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getADCMap
// Purpose       : Get instance params for all ADCs in the circuit
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and MOdels
// Creation Date : 05/05/04
//-----------------------------------------------------------------------------
bool DeviceMgr::getADCMap(map<string,map<string,double> >&ADCMap)
{

  bool bSuccess = true;

  ADCMap.clear();

  // find the device index for ADC devices:

  int adcIndex = getDeviceIndex("Y%ADC%DUMMYADC");

  // get the singleton for the ADC device
  Device* devicePtr = returnDevicePtr(adcIndex);

  // If it wasn't found, there were no ADC devices
  if (devicePtr == NULL || devicePtr == dummy_ptr)
  {
    string msg("DeviceMgr::getADCMap:");
    msg += "  No singleton found, circuit does not apparently have any ADC devices in it.\n";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_WARNING_0, msg);
    // this is not a failure, though
  }
  else
  {
    N_DEV_ADCMaster* ADCPtr = dynamic_cast<N_DEV_ADCMaster*>(devicePtr);
    bSuccess = ADCPtr->getADCMap( ADCMap );
  }

  return bSuccess;
}

//----------------------------------------------------------------------------
// Function       : updateTimeVoltagePairs
// Purpose        : Update the DAC devices in a circuit by adding the set
//                  of time and voltage pairs built up on the "digital side"
//                  since the last update and by removing the time-voltage
//                  pairs for times that pre-date the given simulation time.
// Special Notes  : The current method for locating DAC instances
//                  works for the serial case only. The parallel
//                  case will be added in later modifications.
//                Function moved from N_CIR_Xyce class by TVR, 05/07/2004
// Scope          :
// Creator        : Tom Russo
// Creation Date  : 05/07/2004
//----------------------------------------------------------------------------
bool DeviceMgr::updateTimeVoltagePairs(
   map< string, vector< pair<double,double> >* > const & timeVoltageUpdateMap)
{
  bool bSuccess = true;

  // Get the device index for DAC devices. Note that Xyce's convention
  // for DAC device naming breaks the Netlist usual conventions.
  // Also note that the Device Instance name supplied here does
  // not have to be for an instance that really exists, just one
  // that meets the naming convention for the given type of device.
  int dacIndex = getDeviceIndex("Y%DAC%DUMMYDAC");

  // Get the singleton for the "DAC" type of device.
  Device* devicePtr = returnDevicePtr(dacIndex);

  if (devicePtr == NULL || devicePtr == dummy_ptr)
  {
    // not a failure, there are no DACs
  }
  else
  {
    // Get all of the instances associated with devicePtr.
    vector<DeviceInstance*> deviceInstances;
    devicePtr->getInstances(deviceInstances);

    // Update each of the DAC device instances in turn.
    vector<DeviceInstance*>::iterator iterI;
    vector<DeviceInstance*>::iterator beginI = deviceInstances.begin();
    vector<DeviceInstance*>::iterator endI = deviceInstances.end();
    for (iterI = beginI; iterI != endI; ++iterI)
    {
      N_DEV_DACInstance* dacInstancePtr = dynamic_cast<N_DEV_DACInstance*>(*iterI);
      // Get the name of the given DAC device instance.
      string const & dacName(dacInstancePtr->getName());
      // Try to find an entry for that DAC in the map.
      map< string, vector< pair<double,double> >* >::const_iterator mapIter;
      mapIter = timeVoltageUpdateMap.find(dacName);
      if (mapIter == timeVoltageUpdateMap.end())
      {
        // See if there is an entry for a stripped down version
        // of the DAC name.
        //string strippedDacName(dacName.substr(
        //   dacName.find_last_of("%") + 1, dacName.length() - 1 ));

        mapIter = timeVoltageUpdateMap.find(dacName);
        if (mapIter == timeVoltageUpdateMap.end())
        {
          string msg("In N_CIR_Xyce::updateTimeVoltagePairs: ");
          msg += "Failed to find a map entry for the DAC named ";
          msg += dacName;
          msg += "\n";
          N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_WARNING_0, msg);
          continue;
        }
      }

#ifdef Xyce_DEBUG_DEVICE
      std::cout << "DeviceMgr::updateTimeVoltagePairs: found DAC with name: "
           << dacName << "\n";
#endif

      // Update the time-voltage pairs for the given DAC instance.
      if (!dacInstancePtr->updateTVVEC(*(*mapIter).second))
      {
        string msg("In N_CIR_Xyce::updateTimeVoltagePairs: ");
        msg += "Failed to update the time-voltage pairs for the DAC named ";
        msg += dacName;
        msg += "\n";
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_WARNING_0, msg);
        continue;
      }
    }
  }

  return bSuccess;
}

//----------------------------------------------------------------------------
// Function      : DeviceMgr::setADCWidths
// Purpose       : Set the widths of every ADC instance given a map of names
//                 and widths
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 07/29/2002
//----------------------------------------------------------------------------
bool DeviceMgr::setADCWidths(map<string,int> const & ADCWidthMap)
{
  bool bSuccess = true;
#ifdef Xyce_DEBUG_DEVICE
  std::cout << "entering DeviceMgr::setADCWidths " << std::endl;
#endif

  // Get the device index for ADC devices. Note that Xyce's convention
  // for ADC device naming breaks the Netlist usual conventions.
  // Also note that the Device Instance name supplied here does
  // not have to be for an instance that really exists, just one
  // that meets the naming convention for the given type of device.
  int adcIndex = getDeviceIndex("Y%ADC%DUMMYADC");

  // Get the singleton for the "ADC" type of device.
  Device* devicePtr = returnDevicePtr(adcIndex);

  if (devicePtr == NULL || devicePtr == dummy_ptr)
  {
    // This is not a failure, just means there aren't any ADCs
  }
  else
  {
    // Get all of the instances associated with devicePtr.
    vector<DeviceInstance*> deviceInstances;
    devicePtr->getInstances(deviceInstances);

    // Update each of the ADC device instances in turn.
    vector<DeviceInstance*>::iterator iterI;
    vector<DeviceInstance*>::iterator beginI = deviceInstances.begin();
    vector<DeviceInstance*>::iterator endI = deviceInstances.end();
    for (iterI = beginI; iterI != endI; ++iterI)
    {
      N_DEV_ADCInstance* adcInstancePtr = dynamic_cast<N_DEV_ADCInstance*>(*iterI);
      // Get the name of the given ADC device instance.
      string const & adcName(adcInstancePtr->getName());
      // Try to find an entry for that ADC in the map.
      map< string, int >::const_iterator mapIter;
      mapIter = ADCWidthMap.find(adcName);
      if (mapIter == ADCWidthMap.end())
      {
        // See if there is an entry for a stripped down version
        // of the ADC name.
        //string strippedAdcName(adcName.substr(
        //   adcName.find_last_of("%") + 1, adcName.length() - 1 ));

        mapIter = ADCWidthMap.find(adcName);
        if (mapIter == ADCWidthMap.end())
        {
          string msg("In N_CIR_Xyce::setADCWidths: ");
          msg += "Failed to find a map entry for the ADC named ";
          msg += adcName;
          msg += "\n";
          N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_WARNING_0, msg);
          continue;
        }
      }

#ifdef Xyce_DEBUG_DEVICE
      std::cout << "DeviceMgr::setADCWidths found ADC with name: "
           << adcName << " and will set bit vector width to: "
           << (*mapIter).second << "\n";
#endif

      // Update the time-voltage pairs for the given DAC instance.
      if (!adcInstancePtr->setBitVectorWidth((*mapIter).second))
      {
        string msg("In N_CIR_Xyce::setADCWidths: ");
        msg += "Failed to update the bit vector width for the ADC named ";
        msg += adcName;
        msg += "\n";
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_WARNING_0, msg);
        continue;
      }
    }
  }

  return bSuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getTimeVoltagePairs
// Purpose       : Get time/voltage pairs from all ADCs in the circuit
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and MOdels
// Creation Date : 05/05/04
//-----------------------------------------------------------------------------
bool DeviceMgr::getTimeVoltagePairs(
   map<string, vector< pair<double,double> > >&TimeVoltageMap)
{
  bool bSuccess = true;

  TimeVoltageMap.clear();

  // find the device index for ADC devices:

  int adcIndex = getDeviceIndex("Y%ADC%DUMMYADC");

  // get the singleton for the ADC device
  Device* devicePtr = returnDevicePtr(adcIndex);

  // If it wasn't found, there were no ADC devices
  if (devicePtr == NULL || devicePtr == dummy_ptr)
  {
    // this is not a failure, though, and no point complaining
    //    string msg("DeviceMgr::getADCDevices:");
    //    msg += "  No singleton found, circuit does not apparently have any ADC devices in it.\n";
    //    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_WARNING_0, msg);

  }
  else
  {
    N_DEV_ADCMaster* ADCPtr = dynamic_cast<N_DEV_ADCMaster*>(devicePtr);
    bSuccess = ADCPtr->getTimeVoltagePairs( TimeVoltageMap );
  }

  return bSuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getDeviceNames
// Purpose       : Get a list of names of all devices of given type in the
//                 circuit
// Special Notes : "deviceType" is a string of the form that a device of the
//                  given type would have in the netlist, e.g. "R1" for a
//                  resistor or "Y%XYGRA%DUMMY" for a Xygra device.
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/25/08
//-----------------------------------------------------------------------------
bool DeviceMgr::getDeviceNames(const string & deviceType,
                                     vector < string > & deviceNames)
{

  bool bSuccess = false;

  deviceNames.clear();

  // find the device index for given type of devices:

  int deviceIndex = getDeviceIndex(deviceType);

  // get the singleton for the device
  Device* devicePtr = returnDevicePtr(deviceIndex);

  // If it wasn't found, there were no devices of this type
  if (devicePtr == NULL || devicePtr == dummy_ptr)
  {
    string msg("DeviceMgr::getDeviceNames:");
    msg += "  No singleton found, circuit does not apparently have any devices matching " + deviceType + "in it.\n";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_WARNING_0, msg);
    bSuccess = false;
    // this is not a failure, though
  }
  else
  {
    bSuccess = devicePtr -> getDeviceNames ( deviceNames );
  }

  return bSuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getXygraInstancePtr_
// Purpose       : Returns the pointer to a named Xygra device instance
// Special Notes :
// Scope         : PRIVATE
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/27/08
//-----------------------------------------------------------------------------
Xygra::Instance * DeviceMgr::getXygraInstancePtr_(const string & deviceName)
{
  map<string,Xygra::Instance *>::iterator mapIter;

  // See if we've looked this up before.
  mapIter = xygraPtrMap_.find(deviceName);

  if (mapIter == xygraPtrMap_.end())
  {

    // we haven't looked it up before, do so now.
    // find the device index for given type of devices:
    int deviceIndex = getDeviceIndex("Y%XYGRA%DUMMY");

    // get the singleton for the device
    Device* devicePtr = returnDevicePtr(deviceIndex);
    // If it wasn't found, there were no devices of this type
    if (devicePtr == NULL || devicePtr == dummy_ptr)
    {
      string msg("DeviceMgr::getXygraInstancePtr_:");
      msg += "  No singleton found, circuit does not apparently have any Xygra devices in it.\n";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL, msg);
      return 0; // never reached, just shut up compiler warnings
    }
    else
    {
      // we have some, so get 'em.
      vector<DeviceInstance*> deviceInstances;
      devicePtr->getInstances(deviceInstances);

      int theNamedInstanceIndex = -1;

      int numInstances = deviceInstances.size();
      for (int i = 0 ; i < numInstances; ++i)
      {
        if (deviceInstances[i]->getName() == deviceName)
        {
          theNamedInstanceIndex = i;
          break;
        }
      }
      if (theNamedInstanceIndex == -1)
      {
        string msg("DeviceMgr::getXygraInstancePtr_:");
        msg += "  No Xygra device named " + deviceName + " found.\n";
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL, msg);
      }
      Xygra::Instance* xygraInstancePtr = dynamic_cast<Xygra::Instance*>(deviceInstances[theNamedInstanceIndex]);
      xygraPtrMap_[deviceName] = xygraInstancePtr;
    }
  }

  // when we get here we've either already had the pointer in our map, or
  // we just put it there.  All other cases were fatal errors and would have
  // kept us from getting here.
  return xygraPtrMap_[deviceName];
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::xygraGetNumNodes
// Purpose       : Returns the number of nodes that a given Xygra instance
//                 has, given the name of that instance in the netlist.
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/27/08
//-----------------------------------------------------------------------------
int DeviceMgr::xygraGetNumNodes(const string & deviceName)
{
  Xygra::Instance * xygraInstancePtr = getXygraInstancePtr_(deviceName);
  return xygraInstancePtr -> getNumNodes();
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::xygraGetNumWindings
// Purpose       : Returns the number of windings that a given Xygra instance
//                 has, given the name of that instance in the netlist.
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/27/08
//-----------------------------------------------------------------------------
int DeviceMgr::xygraGetNumWindings(const string & deviceName)
{
  Xygra::Instance * xygraInstancePtr = getXygraInstancePtr_(deviceName);
  return xygraInstancePtr -> getNumWindings();
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::xygraGetCoilWindings
// Purpose       : Returns the number of windings that a given Xygra instance
//                 has in each coil, given the name of that instance in the
//                 netlist.
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 09/15/08
//-----------------------------------------------------------------------------
void DeviceMgr::xygraGetCoilWindings(const string & deviceName,
                                           vector<int> & cW)
{
  Xygra::Instance * xygraInstancePtr = getXygraInstancePtr_(deviceName);
  xygraInstancePtr -> getCoilWindings(cW);
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::xygraGetCoilNames
// Purpose       : Returns the names of each coil in a given Xygra instance
//                 given the name of that instance in the netlist.
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 09/29/08
//-----------------------------------------------------------------------------
void DeviceMgr::xygraGetCoilNames(const string & deviceName,
                                           vector<string> & cN)
{
  Xygra::Instance * xygraInstancePtr = getXygraInstancePtr_(deviceName);
  xygraInstancePtr -> getCoilNames(cN);
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::xygraSetConductances
// Purpose       : Sets th conductance matrix on the specified Xygra device
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/27/08
//-----------------------------------------------------------------------------
bool DeviceMgr::xygraSetConductances(const string & deviceName,
                                           const vector<vector<double> > & cM)
{
  Xygra::Instance * xygraInstancePtr = getXygraInstancePtr_(deviceName);
  return xygraInstancePtr -> setConductances(cM);
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::xygraSetK
// Purpose       : Sets the K matrix on the specified Xygra device
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/27/08
//-----------------------------------------------------------------------------
bool DeviceMgr::xygraSetK(const string & deviceName,
                                const vector<vector<double> > & kM,
                                const double t)
{
  Xygra::Instance * xygraInstancePtr = getXygraInstancePtr_(deviceName);
  return xygraInstancePtr -> setK(kM,t);
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::xygraSetSources
// Purpose       : Set the S vector on named device
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/27/08
//-----------------------------------------------------------------------------
bool DeviceMgr::xygraSetSources(const string & deviceName,
                                      const vector<double> & sV,
                                      const double t)
{
  Xygra::Instance * xygraInstancePtr = getXygraInstancePtr_(deviceName);
  return xygraInstancePtr -> setSources(sV,t);
}


//-----------------------------------------------------------------------------
// Function      : DeviceMgr::xygraGetVoltages
// Purpose       : Retrieve the voltages on nodes for named Xygra device
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/27/08
//-----------------------------------------------------------------------------
bool DeviceMgr::xygraGetVoltages(const string & deviceName,
                                       vector<double> & vN)
{
  Xygra::Instance * xygraInstancePtr = getXygraInstancePtr_(deviceName);
  return xygraInstancePtr -> getVoltages(vN);
}


//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getHomotopyBlockSize
// Purpose       : Returns the number of mosfet gainscale blocks (a value
//                 in the device options).
// Special Notes : Needed for block homotopy on gainscale.
// Scope         : public
// Creator       : Roger P. Pawlowski, SNL, Computational Sciences
// Creation Date : 01/26/2005
//-----------------------------------------------------------------------------
int DeviceMgr::getHomotopyBlockSize() const
{
  return devOptions_.numGainScaleBlocks;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::updateTemperature
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 10/26/04
//-----------------------------------------------------------------------------
bool DeviceMgr::updateTemperature (double val)
{
  bool bsuccess = true, success = true;

  // convert to kelvin:
  double Ctemp = val;
  double Ktemp = val + CONSTCtoK;

#ifdef Xyce_DEBUG_DEVICE
  if (devOptions_.debugLevel > 0 && solState_.debugTimeFlag)
  {
    std::cout << "In DeviceMgr::updateTemperature."
      "  new C temp = " << Ctemp << " K temp = " << Ktemp;
    std::cout << std::endl;
  }
#endif

  // First set the global temp.  This is used in each device if the "tempGiven"
  // variable is false.  This should be in Kelvin.
  devOptions_.temp.setVal(Ktemp);

  // loop over the bsim3 models and delete the size dep params.
  vector<DeviceModel *>::const_iterator iterM;
  vector<DeviceModel *>::const_iterator beginM = bsim3ModelPtrVec_.begin ();
  vector<DeviceModel *>::const_iterator endM = bsim3ModelPtrVec_.end ();
  for (iterM=beginM;iterM!=endM;++iterM)
  {
    (*iterM)->clearTemperatureData ();
  }

  // loop over the bsim4 models and delete the size dep params.
  beginM = bsim4ModelPtrVec_.begin ();
  endM = bsim4ModelPtrVec_.end ();
  for (iterM=beginM;iterM!=endM;++iterM)
  {
    (*iterM)->clearTemperatureData ();
  }

  // loop over the b3soi models and delete the size dep params.
  beginM = bsimsoiModelPtrVec_.begin ();
  endM = bsimsoiModelPtrVec_.end ();
  for (iterM=beginM;iterM!=endM;++iterM)
  {
    (*iterM)->clearTemperatureData ();
  }


  // Loop over all models, call processParams, with CTemp.
  // "XYCEADMS*TEMP is there to force Verilog devices, which might have
  // temperature dependence through "$temperature" instead of a "TEMP"
  // parameter, to work properly.  If so, they need the temperature set in
  // Kelvin.
  string tname("TEMP");
  string tname2("XYCEADMSMODTEMP");
  string tname3("XYCEADMSINSTTEMP");
  beginM = modelPtrVec_.begin();
  endM   = modelPtrVec_.end();
  for (iterM=beginM; iterM!=endM;++iterM)
  {
    success = (*iterM)->setParam (tname, Ctemp);
    success = (*iterM)->setParam (tname2, Ktemp) || success;
    success = success && (*iterM)->processParams ();
  }

  // Loop over device instances, and set the temperature.  This should be
  // in C, if going through processParams, and K if going through
  // the updateTemperature function.
  vector<DeviceInstance*>::const_iterator iter;
  vector<DeviceInstance*>::const_iterator begin;
  vector<DeviceInstance*>::const_iterator end;

  begin = instancePtrVec_.begin();
  end   = instancePtrVec_.end();
  for (iter=begin; iter!=end;++iter)
  {
    success = (*iter)->setParam (tname, Ctemp);
    success = (*iter)->setParam (tname3, Ktemp) || success;
    success = success && (*iter)->processParams ();
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::allDevsConverged
// Purpose       : Check whether any device has taken an action that renders
//                  normal convergence checks invalid (i.e. that the current
//                  step must be assumed unconverged).
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 03/22/05
//-----------------------------------------------------------------------------
bool DeviceMgr::allDevsConverged()
{
  bool allDevsConv = true;
#ifdef Xyce_PARALLEL_MPI
  double dc_par;
  double dc_glob;
#endif

  // if two-level, and just the inner problem, only check the PDE devices,
  // coz those are all we loaded.
  vector<DeviceInstance*>::iterator iter;
  vector<DeviceInstance*>::iterator begin;
  vector<DeviceInstance*>::iterator end;
  if (solState_.twoLevelNewtonCouplingMode==INNER_PROBLEM)
  {
    begin = pdeInstancePtrVec_.begin ();
    end = pdeInstancePtrVec_.end ();
  }
  else
  {
    begin = instancePtrVec_.begin ();
    end = instancePtrVec_.end ();
  }
  for (iter=begin; iter!=end;++iter)
  {
    bool tmpBool =  (*iter)->isConverged();
    allDevsConv = allDevsConv && tmpBool;
  }

#ifdef Xyce_PARALLEL_MPI
  // Take care of moving this info across all processors and getting
  // a global answer
  dc_glob=0.0;

  if (allDevsConv)
    dc_par=0.0;
  else
    dc_par=1.0;

  N_PDS_Comm *pdsCommPtr = pdsMgrPtr_->getPDSComm();
  pdsCommPtr->barrier();
  pdsCommPtr->sumAll(&dc_par, &dc_glob, 1);

  // If any processor has allDevConverged_ false, dc_glob will be nonzero.
  if (dc_glob == 0.0)
    allDevsConv = true;
  else
    allDevsConv = false;
#endif

#ifdef Xyce_DEBUG_DEVICE
  if (devOptions_.debugLevel > 0)
  {
    if (allDevsConv)
    {
      std::cout << "All devices converged!" << std::endl;
    }
    else
    {
      std::cout << "At least one device NOT converged!" << std::endl;
    }
  }
#endif

  return allDevsConv;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::innerDevsConverged
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/22/06
//-----------------------------------------------------------------------------
bool DeviceMgr::innerDevsConverged()
{
  bool innerDevsConv = true;

#ifdef Xyce_EXTDEV

#ifdef Xyce_PARALLEL_MPI
  double dc_par;
  double dc_glob;
#endif

  vector<N_DEV_ExternDeviceInstance*>::iterator iter;
  vector<N_DEV_ExternDeviceInstance*>::iterator begin;
  vector<N_DEV_ExternDeviceInstance*>::iterator end;
  begin = extDevInstancePtrVec_.begin ();
  end = extDevInstancePtrVec_.end ();

  for (iter=begin; iter!=end;++iter)
  {
    bool tmpFlag = (*iter)->isInnerSolveConverged();
    innerDevsConv = innerDevsConv && tmpFlag;
  }

#ifdef Xyce_PARALLEL_MPI
  // Take care of moving this info across all processors and getting
  // a global answer
  dc_glob=0.0;

  if (innerDevsConv)
    dc_par=0.0;
  else
    dc_par=1.0;

  N_PDS_Comm *pdsCommPtr = pdsMgrPtr_->getPDSComm();
  pdsCommPtr->barrier();
  pdsCommPtr->sumAll(&dc_par, &dc_glob, 1);

  // If any processor has innerDevConverged_ false, dc_glob will be nonzero.
  if (dc_glob == 0.0)
    innerDevsConv = true;
  else
    innerDevsConv = false;
#endif

#endif // ext-dev

  return innerDevsConv;
}

#ifdef Xyce_EXTDEV
//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setupExternalDevices
// Purpose       : In parallel, we need to setup all external devices
//                 and appropriately setup the list of instances
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Elec. & MicroSystems Modeling
// Creation Date : 03/10/06
//-----------------------------------------------------------------------------
bool DeviceMgr::setupExternalDevices()
{
#ifdef Xyce_PARALLEL_MPI
    N_PDS_Comm * pdsCommPtr = pdsMgrPtr_->getPDSComm();
    int procID = pdsCommPtr->procID();
    int numProc = pdsCommPtr->numProc();

    int numExt = extDevInstancePtrVec_.size();
    int numExtTotal = 0;
    pdsCommPtr->sumAll( &numExt, &numExtTotal, 1 );

    vector<ExternDevice::Instance*> tmpVec = extDevInstancePtrVec_;

    //resetup the instance vector to have a size totally the global
    //number of ext devices
    if( numExtTotal > 0 )
    {
      extDevInstancePtrVec_.resize( numExtTotal );

      int loc = 0;
      //      const int bufSize = 1000;
      //      char buf[bufSize];
      for( int proc = 0; proc < numProc; ++proc )
      {
        int cnt = 0;
        if( proc == procID ) cnt = numExt;
        pdsCommPtr->bcast( &cnt, 1, proc );

        for( int i = 0; i < cnt; ++i )
        {
          if( proc == procID )
          {
            int size = extDevIBPtrVec_[i]->packedByteCount();
            int bufSize = size+100;
            char *buf=new char[bufSize];
            pdsCommPtr->bcast( &size, 1, proc );
            int pos = 0;
            extDevIBPtrVec_[i]->pack( buf, bufSize, pos, pdsCommPtr );
            pdsCommPtr->bcast( buf, size, proc );
            extDevInstancePtrVec_[loc] = tmpVec[i];
            delete [] buf;
          }
          else
          {
            int size = 0;
            pdsCommPtr->bcast( &size, 1, proc );
            int bufSize = size+100;
            char *buf=new char[bufSize];
            pdsCommPtr->bcast( buf, size, proc );
            int pos = 0;
            InstanceBlock IB;
            IB.unpack( buf, bufSize, pos, pdsCommPtr );
            extDevInstancePtrVec_[loc] = addExtDeviceInstance_( IB );
            delete [] buf;
          }
          extDevInstancePtrVec_[loc]->setOwningProc(proc);
          extDevInstancePtrVec_[loc]->setComm(pdsCommPtr);
          ++loc;
        }
      }

      assert( loc == numExtTotal );
    }

#endif
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::updateExternalDevices_
// Purpose       : Do the actual solve of the external devices
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Elec. & MicroSystems Modeling
// Creation Date : 03/10/06
//-----------------------------------------------------------------------------
void DeviceMgr::updateExternalDevices_()
{
  vector<N_DEV_ExternDeviceInstance*>::iterator iter;
  vector<N_DEV_ExternDeviceInstance*>::iterator begin;
  vector<N_DEV_ExternDeviceInstance*>::iterator end;
  begin = extDevInstancePtrVec_.begin ();
  end = extDevInstancePtrVec_.end ();
  for (iter=begin; iter!=end;++iter)
  {
    (*iter)->runExternalDevice();
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::addExtDeviceInstance_
// Purpose       : adds an external device instance on the processors
//               : that don't actually own it.
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/16/06
//-----------------------------------------------------------------------------
N_DEV_ExternDeviceInstance *
DeviceMgr::addExtDeviceInstance_(InstanceBlock & IB)
{
  DeviceInstance * DI_ptr;
  N_DEV_ExternDeviceInstance * EDI_ptr;

  int  type;

  if (IB.getModelName() == "")
  {
    if( IB.getName().find_last_of(":") != string::npos )
    {
      type = getDeviceIndex(IB.getName().substr(IB.getName().find_last_of(":")+1, string::npos) );
    }
    else
    {
      type = getDeviceIndex( IB.getName() );
    }
  }
  else
  {
    type = deviceModelNameMap_[IB.getModelName()];
  }

  if (type == 0)
  {
     string msg("Device Manager could not correctly figure out what");
     msg += " type of device " + IB.getName() + " is.\n";
     if (IB.getModelName() != "")
     {
       msg += "Its model name is " + IB.getModelName() +", which was not found.\n";
     }
     else
     {
       msg += "Could not find this device name in the device index.\n";
     }
     N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
  }

  // Just in case, call the device creator function - if it has already
  // been allocated it won't create a redundant one.

  createDeviceByModelType(type);

  string deviceTypeName;
  // Add an instance of this type.
  DI_ptr = deviceArray_[type]->addInstance(IB, mlData, solState_, extData_,devOptions_, deviceTypeName);
  EDI_ptr = dynamic_cast<N_DEV_ExternDeviceInstance*>(DI_ptr);

  return EDI_ptr;
}
#endif // Xyce_EXTDEV

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::homotopyStepSuccess
// Purpose       :
// Special Notes : Needed for 2-level
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
void DeviceMgr::homotopyStepSuccess
      (const vector<string> & paramNames,
       const vector<double> & paramVals)
{
#ifdef Xyce_EXTDEV
  int numExt = extDevInstancePtrVec_.size();
  vector<N_DEV_ExternDeviceInstance*>::iterator iter;
  vector<N_DEV_ExternDeviceInstance*>::iterator begin;
  vector<N_DEV_ExternDeviceInstance*>::iterator end;
  begin = extDevInstancePtrVec_.begin ();
  end = extDevInstancePtrVec_.end ();
  for (iter=begin; iter!=end;++iter)
  {
    (*iter)->homotopyStepSuccess (paramNames, paramVals);
  }
#endif // Xyce_EXTDEV

  return;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::homotopyStepFailure
// Purpose       :
// Special Notes : Needed for 2-level
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/30/06
//-----------------------------------------------------------------------------
void DeviceMgr::homotopyStepFailure ()
{
#ifdef Xyce_EXTDEV
  int numExt = extDevInstancePtrVec_.size();
  vector<N_DEV_ExternDeviceInstance*>::iterator iter;
  vector<N_DEV_ExternDeviceInstance*>::iterator begin;
  vector<N_DEV_ExternDeviceInstance*>::iterator end;
  begin = extDevInstancePtrVec_.begin ();
  end = extDevInstancePtrVec_.end ();
  for (iter=begin; iter!=end;++iter)
  {
    (*iter)->homotopyStepFailure ();
  }
#endif // Xyce_EXTDEV

  return;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::stepSuccess
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
void DeviceMgr::stepSuccess (int analysis)
{
#ifdef Xyce_EXTDEV
  int numExt = extDevInstancePtrVec_.size();
  vector<N_DEV_ExternDeviceInstance*>::iterator iter;
  vector<N_DEV_ExternDeviceInstance*>::iterator begin;
  vector<N_DEV_ExternDeviceInstance*>::iterator end;
  begin = extDevInstancePtrVec_.begin ();
  end = extDevInstancePtrVec_.end ();
  for (iter=begin; iter!=end;++iter)
  {
    (*iter)->stepSuccess (analysis);
  }
#endif // Xyce_EXTDEV

  return;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::stepFailure
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
void DeviceMgr::stepFailure (int analysis)
{
#ifdef Xyce_EXTDEV
  int numExt = extDevInstancePtrVec_.size();
  vector<ExternDevice::Instance*>::iterator iter;
  vector<ExternDevice::Instance*>::iterator begin;
  vector<ExternDevice::Instance*>::iterator end;
  begin = extDevInstancePtrVec_.begin ();
  end = extDevInstancePtrVec_.end ();
  for (iter=begin; iter!=end;++iter)
  {
    (*iter)->stepFailure (analysis);
  }
#endif // Xyce_EXTDEV

  return;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::acceptStep
// Purpose       : Communicate to devices that the current step is accepted
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL
// Creation Date : 01/23/07
//-----------------------------------------------------------------------------
void DeviceMgr::acceptStep()
{
  vector<DeviceInstance*>::iterator iter;
  vector<DeviceInstance*>::iterator begin;
  vector<DeviceInstance*>::iterator end;

  // The time history for the LTRA device(s) has to be tracked
  // separately because it can be truncated if the user specifies that
  // option. This has to be called before
  // LTRAInstance::acceptStep(). Note that the DCOP is stored at
  // index zero and the first time step is stored at index 1 for the
  // time history vectors associated with the LTRA
  if (solState_.dcopFlag)
  {
    solState_.ltraTimeIndex = 0;
    solState_.ltraTimeHistorySize = 10;
    solState_.ltraTimePoints.resize(solState_.ltraTimeHistorySize);
  }
  else
  {
    solState_.ltraTimeIndex++;
    if (solState_.ltraTimeIndex >= solState_.ltraTimeHistorySize)
    {
      solState_.ltraTimeHistorySize += 10;
      solState_.ltraTimePoints.resize(solState_.ltraTimeHistorySize);
    }
    solState_.ltraTimePoints[solState_.ltraTimeIndex] = solState_.currTime;
  }

  // make sure we have the right solver state.  Putting this here
  // makes it essential that acceptStep be called BEFORE any changes have
  // been made to the solution or state vectors (e.g. rotating them) and
  // while "nextTime" in the solver still references the time for which the
  // solution vector is valid.

  bool tmpBool = setupSolverInfo_();
  solState_.acceptedTime = solState_.currTime;

  begin = instancePtrVec_.begin ();
  end   = instancePtrVec_.end ();

  for (iter=begin; iter!=end;++iter)
  {
    (*iter)->acceptStep();
  }

  // If the TRYCOMPACT option is set then the LTRA model will try to
  // compact the amount of data stored and speed up the convolutions. If
  // any of the LTRA instances request this, done in their acceptStep()
  // member function, then perform the time-step compaction here.
  if (solState_.ltraDoCompact)
  {
    solState_.ltraTimePoints[solState_.ltraTimeIndex-1] =
      solState_.ltraTimePoints[solState_.ltraTimeIndex];

    solState_.ltraTimeIndex--;

    // reset the flag for the next time step
    solState_.ltraDoCompact = false;
  }
}


//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getInitialQnorm
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/12/07
//-----------------------------------------------------------------------------
bool DeviceMgr::getInitialQnorm (vector <N_TIA_TwoLevelError> & tleVec)
{
  bool bsuccess = true;

#ifdef Xyce_EXTDEV
  int numExt = extDevInstancePtrVec_.size();

  tleVec.resize(numExt);

  vector<ExternDevice::Instance*>::iterator iter;
  vector<ExternDevice::Instance*>::iterator begin;
  vector<ExternDevice::Instance*>::iterator end;
  begin = extDevInstancePtrVec_.begin ();
  end = extDevInstancePtrVec_.end ();

  int iext=0;
  for (iter=begin; iter!=end;++iter,++iext)
  {
    bool bs1 = (*iter)->getInitialQnorm (tleVec[iext]);
    bsuccess = bsuccess && bs1;
  }
#endif // Xyce_EXTDEV

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getInnerLoopErrorSums
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
bool DeviceMgr::getInnerLoopErrorSums (
  vector <N_TIA_TwoLevelError> & tleVec)
{
  bool bsuccess = true;

#ifdef Xyce_EXTDEV
  int numExt = extDevInstancePtrVec_.size();

  tleVec.resize(numExt);

  vector<ExternDevice::Instance*>::iterator iter;
  vector<ExternDevice::Instance*>::iterator begin;
  vector<ExternDevice::Instance*>::iterator end;
  begin = extDevInstancePtrVec_.begin ();
  end = extDevInstancePtrVec_.end ();

  int iext=0;
  for (iter=begin; iter!=end;++iter,++iext)
  {
    bool bs1 = (*iter)->getInnerLoopErrorSum (tleVec[iext]);
    bsuccess = bsuccess && bs1;
  }
#endif // Xyce_EXTDEV

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::updateStateArrays()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
bool DeviceMgr::updateStateArrays()
{
  bool bsuccess = true;

#ifdef Xyce_EXTDEV

  vector<ExternDevice::Instance*>::iterator iter;
  vector<ExternDevice::Instance*>::iterator begin;
  vector<ExternDevice::Instance*>::iterator end;
  begin = extDevInstancePtrVec_.begin ();
  end = extDevInstancePtrVec_.end ();

  for (iter=begin; iter!=end;++iter)
  {
    bool bs1 = (*iter)->updateStateArrays();
    bsuccess = bsuccess && bs1;
  }
#endif // Xyce_EXTDEV

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : DeviceMgr::startTimeStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
bool DeviceMgr::startTimeStep ()
{
  bool bsuccess = true;
  bool tmpBool = setupSolverInfo_();

#ifdef Xyce_EXTDEV

  vector<ExternDevice::Instance*>::iterator iter;
  vector<ExternDevice::Instance*>::iterator begin;
  vector<ExternDevice::Instance*>::iterator end;
  begin = extDevInstancePtrVec_.begin ();
  end = extDevInstancePtrVec_.end ();

  for (iter=begin; iter!=end;++iter)
  {
    bool bs1 = (*iter)->startTimeStep ();
    bsuccess = bsuccess && bs1;
  }
#endif // Xyce_EXTDEV

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setExternalSolverState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
void DeviceMgr::setExternalSolverState (const SolverState & ss)
{
  externalStateFlag_ = true;
  solStateExternal_ = ss;
  return;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setUpDeviceModelNameMap
// Purpose       : This function initializes the deviceModelNameMap structure.
// Special Notes :
//                 This map sets up a relationship between the "type"
//                 parameter, which is one of the manadatory fields in a model
//                 statement, and indices into the device array.  The type
//                 paramter is *not* the *name* of a model, and it says nothing
//                 about the level of the model.  So it is possible that this
//                 map will not fully resolve the deviceArray index for a given
//                 model.
//
//                 For devices (such as MOSFETs) which require a
//                 level-dependent offset to the device array index, it is
//                 neccessary to call two functions to fully resolve the device
//                 array index:
//
//                 int index  = getModelTypeIndex (string & ModelType)
//                 int offset = getDeviceTypeOffset (string & ModelName)
//                 index += offset
//
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/06/00
//-----------------------------------------------------------------------------

void DeviceMgr::setUpDeviceModelTypeMap_()
{
  deviceModelTypeMap_["r"   ] = ModelType::RESISTOR;
  deviceModelTypeMap_["c"   ] = ModelType::CAPACITOR;
  deviceModelTypeMap_["l"   ] = ModelType::INDUCTOR;

#ifdef Xyce_OLD_SWITCH
  deviceModelTypeMap_["sw"  ] = ModelType::SW;
#else
  deviceModelTypeMap_["switch"]  = ModelType::SW;
  deviceModelTypeMap_["iswitch"] = ModelType::SW;
  deviceModelTypeMap_["vswitch"] = ModelType::SW;
#endif

  deviceModelTypeMap_["d"   ] = ModelType::DIODE;
  deviceModelTypeMap_["npn" ] = ModelType::BJT;
  deviceModelTypeMap_["pnp" ] = ModelType::BJT;
  deviceModelTypeMap_["njf" ] = ModelType::JFET;
  deviceModelTypeMap_["pjf" ] = ModelType::JFET;

  deviceModelTypeMap_["nmos"] = ModelType::MOSFET1;
  deviceModelTypeMap_["pmos"] = ModelType::MOSFET1;

  deviceModelTypeMap_["nmf" ] = ModelType::MESFET;
  deviceModelTypeMap_["pmf" ] = ModelType::MESFET;

  deviceModelTypeMap_["zod" ] = ModelType::DIODE_PDE;

#ifdef Xyce_EXTDEV
  deviceModelTypeMap_["ext"] = ModelType::EXTERN_DEVICE;
#endif

  deviceModelTypeMap_["ADC" ] = ModelType::ADC;
  deviceModelTypeMap_["DAC" ] = ModelType::DAC;

  deviceModelTypeMap_["mil"] = ModelType::MUTUAL_INDUCTOR_LINEAR;
  deviceModelTypeMap_["min"] = ModelType::MUTUAL_INDUCTOR_NONLINEAR;
  deviceModelTypeMap_["core"] = ModelType::MUTUAL_INDUCTOR_NONLINEAR;

  deviceModelTypeMap_["opamp" ] = ModelType::OPAMP;
  deviceModelTypeMap_["dig" ] = ModelType::DIGITAL;
  deviceModelTypeMap_["neuron"] = ModelType::NEURON;
  deviceModelTypeMap_["synapse"] = ModelType::SYNAPSE;
  deviceModelTypeMap_["neuronpop"] = ModelType::NEURONPOP;
  deviceModelTypeMap_["newd" ] = ModelType::NEW_DEVICE;
  deviceModelTypeMap_["xygra" ] = ModelType::XYGRA;
  deviceModelTypeMap_["vbic"  ] = ModelType::BJT;

  deviceModelTypeMap_["rom"  ] = ModelType::ROM;
  deviceModelTypeMap_["rxn"  ] = ModelType::RXNSET;

  deviceModelTypeMap_["ltra"] = ModelType::LTRA;

#ifdef Xyce_RAD_MODELS
  ExtendedModelLevels::setUpExtendedDeviceModelTypeMap(deviceModelTypeMap_);
#endif

}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setUpDeviceIndexMap
// Purpose       : This function initializes the deviceIndexMap structure.
//                 This map structure is used to determine, partially, the
//                 index into deviceArray appropriate for a given device
//                 instance name.  For devices that have multiple "levels", or
//                 whatever, the deviceModelNameMap and the MOSFETLevelMap are
//                 needed to further refine the index.
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/31/00
//-----------------------------------------------------------------------------
void DeviceMgr::setUpDeviceIndexMap_()
{
  // This list is of the "canonical" Spice3f5 devices.
  //      lower case map                     upper case map
  deviceIndexMap_["r"] = ModelType::RESISTOR;
  deviceIndexMap_["c"] = ModelType::CAPACITOR;
  deviceIndexMap_["l"] = ModelType::INDUCTOR;
  deviceIndexMap_["d"] = ModelType::DIODE;
  deviceIndexMap_["q"] = ModelType::BJT;
  deviceIndexMap_["j"] = ModelType::JFET;
  deviceIndexMap_["z"] = ModelType::MESFET;
  deviceIndexMap_["m"] = ModelType::MOSFET1;
  deviceIndexMap_["i"] = ModelType::ISRC;

  deviceIndexMap_["e"] = ModelType::VCVS;
  deviceIndexMap_["f"] = ModelType::BSRC;
  deviceIndexMap_["g"] = ModelType::VCCS;
  deviceIndexMap_["h"] = ModelType::BSRC;

  deviceIndexMap_["v"] = ModelType::VSRC;
  deviceIndexMap_["b"] = ModelType::BSRC;
  deviceIndexMap_["o"] = ModelType::LTRA;
  deviceIndexMap_["t"] = ModelType::TRA;
  deviceIndexMap_["s"] = ModelType::SW;

  // "Y" device types, breaks usual netlist syntax rules for devices.
  // Specifed as (for example):  YPDE deviceName node1 node2 ...
  //
  // In the above example, the device instance name coming in from
  // the parser will be Y%PDE%deviceName.

  // PDE based devices.
  deviceIndexMap_["pde"] = ModelType::DIODE_PDE;

#ifdef Xyce_EXTDEV
  deviceIndexMap_["ext"] = ModelType::EXTERN_DEVICE;
#endif

  //
  deviceIndexMap_["adc"] = ModelType::ADC;
  deviceIndexMap_["dac"] = ModelType::DAC;
  deviceIndexMap_["mil"] = ModelType::MUTUAL_INDUCTOR_LINEAR;
  deviceIndexMap_["min"] = ModelType::MUTUAL_INDUCTOR_NONLINEAR;
  deviceIndexMap_["opamp"] = ModelType::OPAMP;
  deviceIndexMap_["NOT"] = ModelType::DIGITAL;
  deviceIndexMap_["AND"] = ModelType::DIGITAL;
  deviceIndexMap_["NAND"] = ModelType::DIGITAL;
  deviceIndexMap_["OR"] = ModelType::DIGITAL;
  deviceIndexMap_["NOR"] = ModelType::DIGITAL;
  deviceIndexMap_["ADD"] = ModelType::DIGITAL;
  deviceIndexMap_["XOR"] = ModelType::DIGITAL;
  deviceIndexMap_["NXOR"] = ModelType::DIGITAL;
  deviceIndexMap_["DFF"] = ModelType::DIGITAL;
  deviceIndexMap_["ACC"] = ModelType::ACC;
  deviceIndexMap_["neuron"] = ModelType::NEURON;
  deviceIndexMap_["synapse"] = ModelType::SYNAPSE;
  deviceIndexMap_["neuronpop"] = ModelType::NEURONPOP;
  deviceIndexMap_["newd"] = ModelType::NEW_DEVICE;
  deviceIndexMap_["xygra"] = ModelType::XYGRA;
  deviceIndexMap_["rom"] = ModelType::ROM;
  deviceIndexMap_["rxn"] = ModelType::RXNSET;

#ifdef Xyce_RAD_MODELS
  ExtendedModelLevels::setUpExtendedDeviceIndexMap (deviceIndexMap_);
#endif
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setUpMESFETLevelMap
// Purpose       : This function initializes the MESFETLevelMap structure.
// Special Notes : These are integer offsets for the deviceArray, which
//                 are determined by what MESFET subtype is being used.
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/31/00
//-----------------------------------------------------------------------------

void DeviceMgr::setUpMESFETLevelMap_()

{
  modelLevelMap_[ModelType::MESFET] = &MESFETLevelMap_;

  // Note:  Do not add to the map until your model is implemented - one of
  // the parser tests in the code is to see if a level exists or not.

  MESFETLevelMap_[1]  = ModelType::MESFET      - ModelType::MESFET;
}


//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setUpMOSFETLevelMap
// Purpose       : This function initializes the MOSFETLevelMap structure.
// Special Notes : These are integer offsets for the deviceArray, which
//                 are determined by what MOSFET subtype is being used.
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/31/00
//-----------------------------------------------------------------------------

void DeviceMgr::setUpMOSFETLevelMap_()

{
  modelLevelMap_[ModelType::MOSFET1] = &MOSFETLevelMap_;

  // The first MOSFET class corresponds to the level=1 MOSFET, so all of the
  // offsets are relative to its index.

  // Note:  Do not add to the map until your model is implemented - one of
  // the parser tests in the code is to see if a level exists or not.

  MOSFETLevelMap_[1]  = ModelType::MOSFET1      - ModelType::MOSFET1;
  MOSFETLevelMap_[2]  = ModelType::MOSFET2      - ModelType::MOSFET1;
  MOSFETLevelMap_[3]  = ModelType::MOSFET3      - ModelType::MOSFET1;
  MOSFETLevelMap_[6]  = ModelType::MOSFET6      - ModelType::MOSFET1;
  MOSFETLevelMap_[9]  = ModelType::MOSFET_B3    - ModelType::MOSFET1;
  MOSFETLevelMap_[49]  = ModelType::MOSFET_B3    - ModelType::MOSFET1;

  MOSFETLevelMap_[10] = ModelType::MOSFET_B3SOI - ModelType::MOSFET1;
  MOSFETLevelMap_[57] = ModelType::MOSFET_B3SOI - ModelType::MOSFET1; // Hspice level number for bsoi

  MOSFETLevelMap_[14] = ModelType::MOSFET_B4    - ModelType::MOSFET1;
  MOSFETLevelMap_[54] = ModelType::MOSFET_B4    - ModelType::MOSFET1; // Hspice level number for b4

  MOSFETLevelMap_[18] = ModelType::VDMOS        - ModelType::MOSFET1;
  MOSFETLevelMap_[103] = ModelType::ADMS_PSP103 - ModelType::MOSFET1;
#ifdef Xyce_NONFREE_MODELS
  MOSFETLevelMap_[301] = ModelType::ADMS_EKV - ModelType::MOSFET1;
#endif

#ifdef Xyce_RAD_MODELS
  ExtendedModelLevels::setUpMOSFETLevelMap(MOSFETLevelMap_);
#endif
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setUpJFETLevelMap
// Purpose       : This function initializes the JFETLevelMap structure.
// Special Notes : These are integer offsets for the deviceArray, which
//                 are determined by what JFET subtype is being used.
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/31/00
//-----------------------------------------------------------------------------

void DeviceMgr::setUpJFETLevelMap_()
{
  modelLevelMap_[ModelType::JFET] = &JFETLevelMap_;

  // The first JFET class corresponds to the level=1 JFET, so all of the
  // offsets are relative to its index.

  // Note:  Do not add to the map until your model is implemented - one of
  // the parser tests in the code is to see if a level exists or not.

  JFETLevelMap_[1]  = ModelType::JFET      - ModelType::JFET;
  JFETLevelMap_[2]  = ModelType::JFET      - ModelType::JFET;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setUpDiodeLevelMap
// Purpose       : This function initializes the DiodeLevelMap structure.
// Special Notes : These are integer offsets for the deviceArray, which
//                 are determined by what diode subtype is being used.
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/31/00
//-----------------------------------------------------------------------------

void DeviceMgr::setUpDiodeLevelMap_()

{
  modelLevelMap_[ModelType::DIODE] = &DiodeLevelMap_;

  // The first diode class corresponds to the level=1 Diode, so all of the
  // offsets are relative to its index.

  // Note:  Do not add to the map until your model is implemented - one of
  // the parser tests in the code is to see if a level exists or not.

  DiodeLevelMap_[1] = ModelType::DIODE - ModelType::DIODE;
  DiodeLevelMap_[2] = ModelType::DIODE - ModelType::DIODE;
#ifdef Xyce_RAD_MODELS
  ExtendedModelLevels::setUpDiodeLevelMap(DiodeLevelMap_);
#endif
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setUpBJTLevelMap
// Purpose       : This function initializes the BJTLevelMap structure.
// Special Notes : These are integer offsets for the deviceArray, which
//                 are determined by what BJT subtype is being used.
// Scope         : private
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 9/25/02
//-----------------------------------------------------------------------------

void DeviceMgr::setUpBJTLevelMap_()

{
  modelLevelMap_[ModelType::BJT] = &BJTLevelMap_;

  // The first BJT class corresponds to the level=1 BJT, so all of the
  // offsets are relative to its index.

  // Note:  Do not add to the map until your model is implemented - one of
  // the parser tests in the code is to see if a level exists or not.

  BJTLevelMap_[1] = ModelType::BJT - ModelType::BJT;
  BJTLevelMap_[10] = ModelType::ADMS_VBIC - ModelType::BJT;
  BJTLevelMap_[23] = ModelType::ADMS_HBT_X - ModelType::BJT;

#ifdef Xyce_RAD_MODELS
  ExtendedModelLevels::setUpBJTLevelMap(BJTLevelMap_);
#endif
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setUpPDELevelMap
// Purpose       : This function initializes the PDELevelMap structure.
// Special Notes : These are integer offsets for the deviceArray, which
//                 are determined by what PDE subtype is being used.
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/31/00
//-----------------------------------------------------------------------------
void DeviceMgr::setUpPDELevelMap_()
{
  modelLevelMap_[ModelType::DIODE_PDE] = &PDELevelMap_;

  // The first PDE class corresponds to the level=1 PDE, so all of the
  // offsets are relative to its index.

  // Note:  Do not add to the map until your model is implemented - one of
  // the parser tests in the code is to see if a level exists or not.

  PDELevelMap_[1] = ModelType::DIODE_PDE   - ModelType::DIODE_PDE;
  PDELevelMap_[2] = ModelType::TWO_D_PDE   - ModelType::DIODE_PDE;
#ifdef Xyce_RAD_MODELS
  ExtendedModelLevels::setUpPDELevelMap(PDELevelMap_);
#endif
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setUpRadLevelMap
// Purpose       : This function initializes the radLevelMap structure.
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 2/11/08
//-----------------------------------------------------------------------------
void DeviceMgr::setUpRadLevelMap_()
{
#ifdef Xyce_RAD_MODELS
  ExtendedModelLevels::setUpRadLevelMap(modelLevelMap_, radLevelMap_);
#endif
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setUpResistorLevelMap
// Purpose       : This function initializes the ResistorLevelMap structure.
// Special Notes : These are integer offsets for the deviceArray, which
//                 are determined by what resistor subtype is being used.
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/31/00
//-----------------------------------------------------------------------------
void DeviceMgr::setUpResistorLevelMap_()
{
  modelLevelMap_[ModelType::RESISTOR] = &ResistorLevelMap_;

  // The first resistor class corresponds to the level=1 Resistor, so all of the
  // offsets are relative to its index.

  // Note:  Do not add to the map until your model is implemented - one of
  // the parser tests in the code is to see if a level exists or not.

  ResistorLevelMap_[1] = ModelType::RESISTOR - ModelType::RESISTOR;
  ResistorLevelMap_[2] = ModelType::THERMAL_RESISTOR - ModelType::RESISTOR;
  ResistorLevelMap_[3] = ModelType::RESISTOR3 - ModelType::RESISTOR;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setUpInductorLevelMap
// Purpose       : This function initializes the Indcutor Level Map structure.
// Special Notes : These are integer offsets for the deviceArray, which
//                 are determined by what nonlinear inductor subtype is being used.
// Scope         : private
// Creator       : Richard Schiek, SNL, Electrical and Microsystems Sim.
// Creation Date : 2/27/08
//-----------------------------------------------------------------------------
void DeviceMgr::setUpIndLevelMap_()
{
  modelLevelMap_[ModelType::MUTUAL_INDUCTOR_NONLINEAR] = &INDLevelMap_;

  // The first non-linear class corresponds to the level=1 non-linear mutual
  // inductor, so all of the offsets are relative to its index.

  // Note:  Do not add to the map until your model is implemented - one of
  // the parser tests in the code is to see if a level exists or not.

  INDLevelMap_[1] = ModelType::MUTUAL_INDUCTOR_NONLINEAR - ModelType::MUTUAL_INDUCTOR_NONLINEAR;
  INDLevelMap_[2] = ModelType::MUTUAL_INDUCTOR_NONLINEAR_2 - ModelType::MUTUAL_INDUCTOR_NONLINEAR;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setUpNeuronLevelMap
// Purpose       : This function initializes the Neuron Level Map structure.
// Special Notes : These are integer offsets for the deviceArray, which
//                 are determined by what nonlinear inductor subtype is being used.
// Scope         : private
// Creator       : Richard Schiek, SNL, Electrical and Microsystems Sim.
// Creation Date : 2/27/08
//-----------------------------------------------------------------------------
void DeviceMgr::setUpNeuronLevelMap_()
{
  modelLevelMap_[ModelType::NEURON] = &NeuronLevelMap_;

  // The first  class corresponds to the level=1 neuron
  // so all of the offsets are relative to its index.

  // Note:  Do not add to the map until your model is implemented - one of
  // the parser tests in the code is to see if a level exists or not.

  NeuronLevelMap_[1] = ModelType::NEURON   - ModelType::NEURON;
  NeuronLevelMap_[2] = ModelType::NEURON_2 - ModelType::NEURON;
  NeuronLevelMap_[3] = ModelType::NEURON_3 - ModelType::NEURON;
  NeuronLevelMap_[4] = ModelType::NEURON_4 - ModelType::NEURON;
  NeuronLevelMap_[5] = ModelType::NEURON_5 - ModelType::NEURON;
  NeuronLevelMap_[6] = ModelType::NEURON_6 - ModelType::NEURON;
  NeuronLevelMap_[7] = ModelType::NEURON_7 - ModelType::NEURON;
  NeuronLevelMap_[8] = ModelType::NEURON_8 - ModelType::NEURON;
  NeuronLevelMap_[9] = ModelType::NEURON_9 - ModelType::NEURON;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setUpNeuronPopLevelMap_
// Purpose       : This function initializes the Neuron Level Map structure.
// Special Notes : These are integer offsets for the deviceArray, which
//                 are determined by what nonlinear inductor subtype is being used.
// Scope         : private
// Creator       : Richard Schiek, SNL, Electrical and Microsystems Sim.
// Creation Date : 2/27/08
//-----------------------------------------------------------------------------
void DeviceMgr::setUpNeuronPopLevelMap_()
{
  modelLevelMap_[ModelType::NEURONPOP] = &NeuronPopLevelMap_;

  // The first  class corresponds to the level=1 neuron pop
  // so all of the offsets are relative to its index.

  // Note:  Do not add to the map until your model is implemented - one of
  // the parser tests in the code is to see if a level exists or not.

  NeuronPopLevelMap_[1] = ModelType::NEURONPOP   - ModelType::NEURONPOP;

}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setUpSynapseLevelMap
// Purpose       : This function initializes the Synapse Level Map structure.
// Special Notes : These are integer offsets for the deviceArray, which
//                 are determined by what nonlinear inductor subtype is being used.
// Scope         : private
// Creator       : Christy Warrender, SNL, Cognitive Modeling
// Creation Date : 11/18/10
//-----------------------------------------------------------------------------
void DeviceMgr::setUpSynapseLevelMap_()
{
  modelLevelMap_[ModelType::SYNAPSE] = &SynapseLevelMap_;

  // The first  class corresponds to the level=1 synapse
  // so all of the offsets are relative to its index.

  // Note:  Do not add to the map until your model is implemented - one of
  // the parser tests in the code is to see if a level exists or not.

  SynapseLevelMap_[1] = ModelType::SYNAPSE   - ModelType::SYNAPSE;
  SynapseLevelMap_[2] = ModelType::SYNAPSE_2 - ModelType::SYNAPSE;
  SynapseLevelMap_[3] = ModelType::SYNAPSE_3 - ModelType::SYNAPSE;
  SynapseLevelMap_[4] = ModelType::SYNAPSE_4 - ModelType::SYNAPSE;
}


//-----------------------------------------------------------------------------
// Function      : DeviceMgr::restartDataSize
// Purpose       :
// Special Notes : 
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 06/26/2013
//-----------------------------------------------------------------------------
int DeviceMgr::restartDataSize( bool pack )
{
  int numdoubles = solState_.ltraTimePoints.size();
  int numSize_t = 3;
  int count = sizeof(double) * (numdoubles);
  count += sizeof(size_t) * numSize_t;

  // bump up size for unpacked data.  This is an empirical multiplier.
  if( !pack ) 
  {
    count *= 3;
  }

  return count;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::dumpRestartData
// Purpose       : Output restart data.
// Special Notes : This function is called by the restart manager to output
//                 persistent data for the device package.  It should NOT
//                 include any data from individual devices, as that restart
//                 data is collected elsewhere.
// Scope         : 
// Creator       : Eric R. Keiter, SNL
// Creation Date : 06/26/2013
//-----------------------------------------------------------------------------
bool DeviceMgr::dumpRestartData
(char * buf, int bsize, int & pos, N_PDS_Comm * comm, bool pack )
{
  bool retval=true;

  if ( pack )
  {
    size_t size=solState_.ltraTimePoints.size();
    comm->pack( &(solState_.ltraTimeIndex), 1, buf, bsize, pos );
    comm->pack( &(solState_.ltraTimeHistorySize), 1, buf, bsize, pos );
    comm->pack( &(size), 1, buf, bsize, pos );
    comm->pack( &(solState_.ltraTimePoints[0]), size, buf, bsize, pos );
  }
  else
  {
    int count = restartDataSize( false );
    int startIndex = pos;
    for( int i = startIndex; i < (startIndex+count); ++i) buf[i] = ' ';

    size_t size=solState_.ltraTimePoints.size();
    ostringstream ost;
    ost.width(24);ost.precision(16);ost.setf(ios::scientific);
    ost << solState_.ltraTimeIndex << " ";
    ost << solState_.ltraTimeHistorySize << " ";
#ifdef Xyce_DEBUG_RESTART
    std::cout <<
      "DeviceMgr::getRestartData:  ltraTimeIndex = " << solState_.ltraTimeIndex <<std::endl;
    std::cout <<
      "DeviceMgr::getRestartData:  ltraTimeHistorySize = " << solState_.ltraTimeHistorySize <<std::endl;
#endif
    ost << size << " ";
    for (int i=0;i<size;i++)
    {
      ost << solState_.ltraTimePoints[i] << " ";
#ifdef Xyce_DEBUG_RESTART
    std::cout <<
      "DeviceMgr::dumpRestartData:  ltraTimePoints["<<i<<"] =" 
      << solState_.ltraTimePoints[i]<<std::endl;
#endif
    }

    string data( ost.str() );
    for( unsigned int i = 0; i < data.length(); ++i ) buf[startIndex+i] = data[i];

    // The line above copies the characters of the data string into buf,
    // but doesn't null-terminate buf.
    // it is essential to terminate the buffer with a null, or attempts
    // to construct a string object from it will get memory access problems.
    buf[startIndex+data.length()] = '\0';
    pos += data.length();
  }

  return retval;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::restoreRestartData
// Purpose       : Load restart data.
// Special Notes : 
// Scope         : 
// Creator       : Eric R. Keiter, SNL
// Creation Date : 06/26/2013
//-----------------------------------------------------------------------------
bool DeviceMgr::restoreRestartData
(char * buf, int bsize, int & pos, N_PDS_Comm * comm, bool pack )
{
  bool retval=true;

  if( pack )
  {
    comm->unpack(buf, bsize, pos, &(solState_.ltraTimeIndex), 1);
    comm->unpack(buf, bsize, pos, &(solState_.ltraTimeHistorySize), 1);
    size_t size=0;
    comm->unpack(buf, bsize, pos, &(size), 1);
    solState_.ltraTimePoints.resize(size);
    comm->unpack(buf, bsize, pos, &(solState_.ltraTimePoints[0]), size);
  }
  else
  {
    string str1(buf);
    int length = str1.size() - pos;
    string str2(str1,pos,length);

    istringstream ist( str2 );

    ist >> solState_.ltraTimeIndex;
    ist >> solState_.ltraTimeHistorySize;
#ifdef Xyce_DEBUG_RESTART
    std::cout <<
      "DeviceMgr::restoreRestartData:  ltraTimeIndex = " << solState_.ltraTimeIndex <<std::endl;
    std::cout <<
      "DeviceMgr::restoreRestartData:  ltraTimeHistorySize = " << solState_.ltraTimeHistorySize <<std::endl;
#endif
    size_t size=0;
    ist >> size;
    solState_.ltraTimePoints.resize(size);
    for (int i=0;i<size;++i)
    {
      ist >> solState_.ltraTimePoints[i];
#ifdef Xyce_DEBUG_RESTART
    std::cout <<
      "DeviceMgr::restoreRestartData:  ltraTimePoints["<<i<<"] =" 
      << solState_.ltraTimePoints[i]<<std::endl;
#endif
    }

    pos += ist.tellg();
  }

  return retval;
}


} // namespace Device
} // namespace Xyce

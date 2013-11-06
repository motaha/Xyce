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
// Filename       : $RCSfile: N_DEV_DeviceInterface.C,v $
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
// Revision Number: $Revision: 1.134.2.3 $
//
// Revision Date  : $Date: 2013/10/03 17:23:38 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ----------   Std Includes  ------------

#include <N_UTL_Misc.h>

#ifdef HAVE_ALGORITHM
#include <algorithm>
#else
#ifdef HAVE_ALGO_H
#include <algo.h>
#else
#error Must have either <algorithm> or <algo.h>!
#endif
#endif

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceInterface.h>
#include <N_DEV_DeviceMgr.h>
#include <N_ERH_ErrorMgr.h>
#include <N_UTL_BreakPoint.h>
#include <N_TIA_TwoLevelError.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::factory
// Purpose       : factory function
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/02
//-----------------------------------------------------------------------------
DeviceInterface * DeviceInterface::factory (N_IO_CmdParse & cp)
{
   DeviceInterface * DI_ptr = new DeviceInterface(cp);
   return DI_ptr;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::DeviceInterface
// Purpose       : constructor
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/02
//-----------------------------------------------------------------------------
DeviceInterface::DeviceInterface(N_IO_CmdParse & cp)
{
  devMgrPtr_ = DeviceMgr::factory(cp);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::~DeviceInterface
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/02
//-----------------------------------------------------------------------------
DeviceInterface::~DeviceInterface()
{
  delete devMgrPtr_;
}


//-----------------------------------------------------------------------------
// Function      : DeviceInterface::returnDevicePtr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/02
//-----------------------------------------------------------------------------
Device * DeviceInterface::returnDevicePtr (int Index)
{
  return devMgrPtr_->returnDevicePtr(Index);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::registerLinearSystem
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/02
//-----------------------------------------------------------------------------
bool DeviceInterface::registerLinearSystem (N_LAS_System * tmp_system_ptr)
{
  return devMgrPtr_->registerLinearSystem(tmp_system_ptr);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::registerAnalysisInterface
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/02
//-----------------------------------------------------------------------------
bool DeviceInterface::registerAnalysisInterface (N_ANP_AnalysisInterface * tmp_anaIntPtr)
{
  return devMgrPtr_->registerAnalysisInterface (tmp_anaIntPtr);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::registerOutputManager
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 08/15/06
//-----------------------------------------------------------------------------
bool DeviceInterface::registerOutputMgr (N_IO_OutputMgr * tmp_outputMgrPtr)
{
  return devMgrPtr_->registerOutputMgr(tmp_outputMgrPtr);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::registerParallelMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/29/01
//-----------------------------------------------------------------------------
bool DeviceInterface::registerParallelMgr (N_PDS_Manager * tmp_pdsMgrPtr )
{
  return devMgrPtr_->registerParallelMgr(tmp_pdsMgrPtr);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::registerNonlinearSolver
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/02
//-----------------------------------------------------------------------------
bool DeviceInterface::registerNonlinearSolver (N_NLS_Manager * tmp_nlsMgrPtr)
{
  return devMgrPtr_->registerNonlinearSolver(tmp_nlsMgrPtr);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::registerPkgOptionsMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Mems Modeling
// Creation Date : 10/20/2008
//-----------------------------------------------------------------------------
bool DeviceInterface::registerPkgOptionsMgr( RCP<N_IO_PkgOptionsMgr> pkgOptPtr )
{
  return devMgrPtr_->registerPkgOptionsMgr(pkgOptPtr);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::getFastSourcePeriod
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/27/04
//-----------------------------------------------------------------------------
vector<double> DeviceInterface::getFastSourcePeriod (vector<string>& sourceNames)
{
  return devMgrPtr_->getFastSourcePeriod (sourceNames);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::registerFastSources
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/27/04
//-----------------------------------------------------------------------------
vector<double> DeviceInterface::registerFastSources(vector<string>& sourceNames)
{
  return devMgrPtr_->registerFastSources (sourceNames);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::deRegisterFastSources
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 06/12/2013
//-----------------------------------------------------------------------------
void DeviceInterface::deRegisterFastSources(vector<string>& sourceNames)
{
  return devMgrPtr_->deRegisterFastSources (sourceNames);
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_DeviceInterface::registerLeadCurrentRequests
// Purpose       : this function is called from the output manager to inform the
//                 device package of the devices for which lead currents have
//                 been requested.  The device manager will take care of doing
//                 isolated F and Q loads for these devices so the lead currents
//                 can be calculated
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Electrical Systems Modeling
// Creation Date : 03/20/13
//-----------------------------------------------------------------------------
bool DeviceInterface::setLeadCurrentRequests( const set<string> & deviceNames )
{
  return devMgrPtr_->setLeadCurrentRequests( deviceNames );
}


//-----------------------------------------------------------------------------
// Function      : N_DEV_DeviceInterface::deactivateSlowSources
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/22/07
//-----------------------------------------------------------------------------
void DeviceInterface::deactivateSlowSources()
{
  devMgrPtr_->deactivateSlowSources();
}


//-----------------------------------------------------------------------------
// Function      : DeviceInterface::activateSlowSources
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/22/07
//-----------------------------------------------------------------------------
void DeviceInterface::activateSlowSources()
{
  devMgrPtr_->activateSlowSources();
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::setMPDEFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Parallel Computational Sciences
// Creation Date : 07/21/08
//-----------------------------------------------------------------------------
void DeviceInterface::setMPDEFlag( bool flagVal )
{
  devMgrPtr_->setMPDEFlag( flagVal );
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::setBlockAnalysisFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Parallel Computational Sciences
// Creation Date : 07/21/08
//-----------------------------------------------------------------------------
void DeviceInterface::setBlockAnalysisFlag( bool flagVal )
{
  devMgrPtr_->setBlockAnalysisFlag( flagVal );
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::setFastTime
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Parallel Computational Sciences
// Creation Date : 07/21/08
//-----------------------------------------------------------------------------
void DeviceInterface::setFastTime( double timeVal )
{
  devMgrPtr_->setFastTime( timeVal );
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::initializeAll
// Purpose       : This function, via the LAS system class, sets up
//                 the pointers to the various linear algebra entities.
// Special Notes :
// Scope         : protected
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/12/00
//-----------------------------------------------------------------------------
bool DeviceInterface::initializeAll()
{
  return devMgrPtr_->initializeAll();
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::resetForStepAnalysis
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/12/2013
//-----------------------------------------------------------------------------
void DeviceInterface::resetForStepAnalysis ()
{
  return devMgrPtr_->resetForStepAnalysis();
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::createDeviceByNetlistDeviceType
// Purpose       : 
// Special Notes :
// Scope         : protected
// Creator       : Dave Baur, SNL
// Creation Date : 04/18/2013
//-----------------------------------------------------------------------------
Device * DeviceInterface::createDeviceByNetlistDeviceType(const std::string &name, const int level)
{
  return devMgrPtr_->createDeviceByNetlistDeviceType(name, level);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::getNumSupportedDevices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/16/05
//-----------------------------------------------------------------------------
int DeviceInterface::getNumSupportedDevices ()
{
  return devMgrPtr_->getNumSupportedDevices ();
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::flushDevices
// Purpose       : This function deletes all the allocated devices, and resets
//                 the array of device pointers to all point to the dummy
//                 (placeholder) device class pointer.
// Special Notes :
// Scope         : protected
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/02
//-----------------------------------------------------------------------------
bool DeviceInterface::flushDevices ()
{
  return devMgrPtr_->flushDevices ();
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::addDeviceModel
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/02
//-----------------------------------------------------------------------------
bool DeviceInterface::addDeviceModel (const ModelBlock & MB)
{
  return devMgrPtr_->addDeviceModel(MB);
}


//-----------------------------------------------------------------------------
// Function      : DeviceInterface::verifyDeviceInstance
// Purpose       : This funciton verifies a device instance prior to instantiating.
//                 Theoretically, we could do this in addDeviceInstance() and
//                 not make a device if it fails some verification criteria (a
//                 resistor with zero resistance is the primary case here).  However,
//                 later unlinking of redundant devices that are connected to just
//                 one node is difficut after addDeviceInstance() is called because
//                 the device instance pointer can be placed in many other containers
//                 It could be done, but this is a simpler first step to having the
//                 device manager be in charge of device verification -- rather
//                 than have it in toplogy or io.
//
// Special Notes : return true if this device is ok to instantiate, false otherwise
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/18/2010
//-----------------------------------------------------------------------------
bool DeviceInterface::verifyDeviceInstance(InstanceBlock & IB)
{
  return devMgrPtr_->verifyDeviceInstance( IB );
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::addDeviceInstance
// Purpose       : addDeviceInstance will create a new instance of the
//                 designated device type.  This version of the function
//                 accepts a parameter list as one of the arguments,
//                 so it is assumed that a parameter instance will
//                 also have to be allocated for it.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/02
//-----------------------------------------------------------------------------
DeviceInstance * DeviceInterface::addDeviceInstance
                                        (InstanceBlock & IB)
{
  return devMgrPtr_->addDeviceInstance(IB);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::deleteDeviceInstance
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/02
//-----------------------------------------------------------------------------
bool DeviceInterface::deleteDeviceInstance (const string & name)
{
  return devMgrPtr_->deleteDeviceInstance(name);
}

// //----------------------------------------------------------------------------
// // Function       : DeviceInterface::getDeviceIndex
// // Purpose        :
// // Special Notes  :
// // Scope          : public
// // Creator        : Lon Waters
// // Creation Date  : 06/09/2003
// //----------------------------------------------------------------------------
// int DeviceInterface::getDeviceIndex(string const& type)
// {
//   return devMgrPtr_->getDeviceIndex(type);
// }

// //----------------------------------------------------------------------------
// // Function       : DeviceInterface::getDeviceTypeOffset
// // Purpose        :
// // Special Notes  :
// // Scope          : public
// // Creator        : Dave Shirley, PSSI
// // Creation Date  : 10/04/2005
// //----------------------------------------------------------------------------
// int DeviceInterface::getDeviceTypeOffset
//     (const string & name, const int level, const string & type)
// {
//   return devMgrPtr_->getDeviceTypeOffset(name, level, type);
// }

//----------------------------------------------------------------------------
// Function       : DeviceInterface::getDeviceCountMap
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Eric R. Keiter, SNL.
// Creation Date  : 05/07/2010
//----------------------------------------------------------------------------
const map<string,int> & DeviceInterface::getDeviceCountMap()
{
  return devMgrPtr_->getDeviceCountMap();
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::printOutLists
// Purpose       : This function will send output to stdout of all the
//                 allocated model and instance lists.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/02
//-----------------------------------------------------------------------------
void DeviceInterface::printOutLists()
{
  devMgrPtr_->printOutLists();
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::loadDeviceMask
// Purpose       : This function loads the mask elements for all
//                 of the circuit devices.
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 01/19/07
//-----------------------------------------------------------------------------
bool DeviceInterface::loadDeviceMask()
{
  return devMgrPtr_->loadDeviceMask();
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::setInitialGuess
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/02/03
//-----------------------------------------------------------------------------
bool DeviceInterface::setInitialGuess
  (N_LAS_Vector * solVectorPtr)
{
  return devMgrPtr_->setInitialGuess(solVectorPtr);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::setParam
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/02/03
//-----------------------------------------------------------------------------
bool DeviceInterface::setParam(string & name, double val)
{
  return devMgrPtr_->setParam (name,val);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::getParam
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/02/03
//-----------------------------------------------------------------------------
double DeviceInterface::getParam(const string & name)
{
  return devMgrPtr_->getParam (name);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::getParam
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/02/03
//-----------------------------------------------------------------------------
bool DeviceInterface::getParam(const string & name, double & val)
{
  return devMgrPtr_->getParam (name, val);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::getVsrcLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/06
//-----------------------------------------------------------------------------
bool DeviceInterface::getVsrcLIDs
    (string & srcName, int & li_Pos, int & li_Neg, int & li_Bra)
{
  return devMgrPtr_-> getVsrcLIDs (srcName,li_Pos,li_Neg,li_Bra);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::updateSources
// Purpose       : This function function updates sources for the present
//                 time step.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/02
//-----------------------------------------------------------------------------
bool DeviceInterface::updateSources()
{
  return devMgrPtr_->updateSources ();
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::setICs
// Purpose       : This function function sets initial conditions for devices
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/02
//-----------------------------------------------------------------------------
bool DeviceInterface::setICs(
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
  return devMgrPtr_->setICs( tmpSolVectorPtr,
                             tmpCurrSolVectorPtr,
                             tmpLastSolVectorPtr,
                             tmpStaVectorPtr,
                             tmpCurrStaVectorPtr,
                             tmpLastStaVectorPtr,
                             tmpStaDerivVectorPtr,
                             tmpStoVectorPtr,
                             tmpCurrStoVectorPtr,
                             tmpLastStoVectorPtr,
                             tmpQVectorPtr,
                             tmpFVectorPtr,
                             tmpdFdxdVpVectorPtr,
                             tmpdQdxdVpVectorPtr );
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::getLinearSystemFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/02
//-----------------------------------------------------------------------------
bool DeviceInterface::getLinearSystemFlag()
{
  return devMgrPtr_->getLinearSystemFlag ();
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::getVoltageLimterFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/11/04
//-----------------------------------------------------------------------------
bool DeviceInterface::getVoltageLimiterFlag ()
{
  return devMgrPtr_->getVoltageLimiterFlag ();
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::getPDESystemFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/02
//-----------------------------------------------------------------------------
bool DeviceInterface::getPDESystemFlag()
{
  return devMgrPtr_->getPDESystemFlag ();
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::runParameterTests
// Purpose       : The  purpose of this function is to test out the processing
//                 of various device parameters.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/02
//-----------------------------------------------------------------------------
bool DeviceInterface::runParameterTests(string & deviceName)
{
  return devMgrPtr_->runParameterTests(deviceName);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::output
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/02
//-----------------------------------------------------------------------------
bool DeviceInterface::output ()
{
  return devMgrPtr_->output ();
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::finishOutput
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 06/19/04
//-----------------------------------------------------------------------------
bool DeviceInterface::finishOutput ()
{
  return devMgrPtr_->finishOutput ();
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface:: dotOpOutput
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 5/3/12
//-----------------------------------------------------------------------------
void DeviceInterface::dotOpOutput ()
{
  return devMgrPtr_->dotOpOutput ();
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::setGlobalFlags
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 11/2/12
//-----------------------------------------------------------------------------
void DeviceInterface::setGlobalFlags()
{
  return devMgrPtr_->setGlobalFlags ();
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::getBreakPoints
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/02
//-----------------------------------------------------------------------------
bool DeviceInterface::getBreakPoints ( vector<N_UTL_BreakPoint> & breakPointTimes )
{

  return devMgrPtr_->getBreakPoints (breakPointTimes);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::getMaxTimeStepSize
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/02
//-----------------------------------------------------------------------------
double DeviceInterface::getMaxTimeStepSize ()
{
  return devMgrPtr_->getMaxTimeStepSize ();
}


//-----------------------------------------------------------------------------
// Function      : DeviceInterface::registerOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/02
//-----------------------------------------------------------------------------
bool DeviceInterface::registerOptions(const N_UTL_OptionBlock & OB)

{
  return devMgrPtr_->registerOptions(OB);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::registerSensParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 06/03/02
//-----------------------------------------------------------------------------
bool DeviceInterface::registerSensParams (const N_UTL_OptionBlock & OB)

{
  return devMgrPtr_->registerSensParams (OB);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::registerICLoads
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 03/05/02
//-----------------------------------------------------------------------------
bool DeviceInterface::registerICLoads( vector< pair<int,double> > * icLoads )
{
  return devMgrPtr_->registerICLoads(icLoads);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::enablePDEContinuation
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 10/22/02
//-----------------------------------------------------------------------------
int DeviceInterface::enablePDEContinuation ()
{
  return devMgrPtr_->enablePDEContinuation ();
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::disablePDEContinuation
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 10/22/02
//-----------------------------------------------------------------------------
bool DeviceInterface::disablePDEContinuation ()
{
  return devMgrPtr_->disablePDEContinuation ();
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::getNumInterfaceNodes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 12/03/02
//-----------------------------------------------------------------------------
void DeviceInterface::getNumInterfaceNodes (vector<int> & numINodes)
{
  devMgrPtr_->getNumInterfaceNodes (numINodes);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::loadCouplingRHS
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 12/03/02
//-----------------------------------------------------------------------------
bool DeviceInterface::loadCouplingRHS (int iPDEDevice, int iElectrode, N_LAS_Vector * dfdvPtr)
{
  return devMgrPtr_->loadCouplingRHS (iPDEDevice, iElectrode,dfdvPtr);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::calcCouplingTerms
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 12/03/02
//-----------------------------------------------------------------------------
bool DeviceInterface::calcCouplingTerms (int iSubProblem, int iElectrode, const N_LAS_Vector * dxdvPtr)
{
  return devMgrPtr_->calcCouplingTerms (iSubProblem, iElectrode, dxdvPtr);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::raiseDebugLevel
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 04/23/03
//-----------------------------------------------------------------------------
bool DeviceInterface::raiseDebugLevel(int increment)
{
  return devMgrPtr_->raiseDebugLevel(increment);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::loadDAEMatrices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 03/03/04
//-----------------------------------------------------------------------------
bool DeviceInterface::loadDAEMatrices (
     N_LAS_Vector * tmpSolVectorPtr,
     N_LAS_Vector * tmpStaVectorPtr,
     N_LAS_Vector * tmpStaDerivVectorPtr,
     N_LAS_Vector * tmpStoVectorPtr,
     N_LAS_Matrix * tmpdQdxMatrixPtr,
     N_LAS_Matrix * tmpdFdxMatrixPtr)
{
  return devMgrPtr_->loadDAEMatrices (
                                  tmpSolVectorPtr,
                                  tmpStaVectorPtr,
                                  tmpStaDerivVectorPtr,
                                  tmpStoVectorPtr,
                                  tmpdQdxMatrixPtr,
                                  tmpdFdxMatrixPtr);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::loadDAEVectors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 02/20/03
//-----------------------------------------------------------------------------
bool DeviceInterface::loadDAEVectors (
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
   N_LAS_Vector * tmpStoLeadCurrQCompVectorPtr,
   N_LAS_Vector * tmpQVectorPtr,
   N_LAS_Vector * tmpFVectorPtr,
   N_LAS_Vector * tmpdFdxdVpVectorPtr,
   N_LAS_Vector * tmpdQdxdVpVectorPtr)

{
  return devMgrPtr_->loadDAEVectors ( tmpSolVectorPtr,
                                      tmpCurrSolVectorPtr,
                                      tmpLastSolVectorPtr,
                                      tmpStaVectorPtr,
                                      tmpCurrStaVectorPtr,
                                      tmpLastStaVectorPtr,
                                      tmpStaDerivVectorPtr,
                                      tmpStoVectorPtr,
                                      tmpCurrStoVectorPtr,
                                      tmpLastStoVectorPtr,
                                      tmpStoLeadCurrQCompVectorPtr,
                                      tmpQVectorPtr,
                                      tmpFVectorPtr,
                                      tmpdFdxdVpVectorPtr,
                                      tmpdQdxdVpVectorPtr );
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::updateState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 11/17/05
//-----------------------------------------------------------------------------
bool DeviceInterface::updateState      (
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
  return devMgrPtr_->updateState (
    nextSolVectorPtr, currSolVectorPtr, lastSolVectorPtr,
    nextStaVectorPtr, currStaVectorPtr, lastStaVectorPtr,
    nextStoVectorPtr, currStoVectorPtr, lastStoVectorPtr
    );
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::setVoltageLimiterFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 04/20/04
//-----------------------------------------------------------------------------
void DeviceInterface::setVoltageLimiterFlag ()
{
  devMgrPtr_->setVoltageLimiterFlag ();
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::unsetVoltageLimiterFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Computational Sciences
// Creation Date : 04/20/04
//-----------------------------------------------------------------------------
void DeviceInterface::unsetVoltageLimiterFlag ()
{
  devMgrPtr_->unsetVoltageLimiterFlag ();
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::loadBVectorsforAC
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool DeviceInterface::loadBVectorsforAC (N_LAS_Vector * bVecRealPtr,
                        N_LAS_Vector * bVecImagPtr)
{
  return devMgrPtr_->loadBVectorsforAC (bVecRealPtr, bVecImagPtr);
}

bool DeviceInterface::getBMatrixEntriesforMOR(std::vector<int>& bMatEntriesVec,
                                                    std::vector<int>& bMatPosEntriesVec)
{
  return devMgrPtr_->getBMatrixEntriesforMOR( bMatEntriesVec, bMatPosEntriesVec );
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::getDACDeviceNames
// Purpose       : Get a list of names of all DACs in the circuit
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and MOdels
// Creation Date : 05/05/04
//-----------------------------------------------------------------------------
bool DeviceInterface::getDACDeviceNames(vector < string > & dacNames)
{
  return devMgrPtr_->getDACDeviceNames(dacNames);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::getADCMap
// Purpose       : Get instance params for all ADCs in the circuit
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and MOdels
// Creation Date : 05/05/04
//-----------------------------------------------------------------------------
bool DeviceInterface::getADCMap(map<string,map<string,double> >&ADCMap)
{
  return devMgrPtr_->getADCMap(ADCMap);
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
bool DeviceInterface::updateTimeVoltagePairs(
   map< string, vector< pair<double,double> >* > const & timeVoltageUpdateMap)
{
  return devMgrPtr_->updateTimeVoltagePairs(timeVoltageUpdateMap);
}

//----------------------------------------------------------------------------
// Function      : DeviceInterface::setADCWidths
// Purpose       : Set the widths of every ADC instance given a map of names
//                 and widths
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 07/29/2002
//----------------------------------------------------------------------------
bool DeviceInterface::setADCWidths(map<string,int> const & ADCWidthMap)
{
  return devMgrPtr_->setADCWidths(ADCWidthMap);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::getTimeVoltagePairs
// Purpose       : Get time/voltage pairs from all ADCs in the circuit
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and MOdels
// Creation Date : 05/05/04
//-----------------------------------------------------------------------------
bool DeviceInterface::getTimeVoltagePairs(
   map<string, vector< pair<double,double> > >&TimeVoltageMap)
{
  return devMgrPtr_->getTimeVoltagePairs(TimeVoltageMap);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::getDeviceNames
// Purpose       : Get a list of names of all devices of given type in the
//                 circuit
// Special Notes : "deviceType" is string of the form that an actual
//                 device in the netlist would have, e.g. "R1" for a resistor
//                 or "Y%XYGRA%DUMMY" for a Xygra device.
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/25/08
//-----------------------------------------------------------------------------
bool DeviceInterface::getDeviceNames(const string & deviceType,
                                           vector < string > & deviceNames)
{
  return devMgrPtr_->getDeviceNames(deviceType, deviceNames);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::xygraGetNumNodes
// Purpose       : Returns the number of nodes that a given Xygra instance
//                 has, given the name of that instance in the netlist.
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/27/08
//-----------------------------------------------------------------------------
int DeviceInterface::xygraGetNumNodes(const string & deviceName)
{
  return devMgrPtr_->xygraGetNumNodes(deviceName);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::xygraGetNumWindings
// Purpose       : Returns the number of windings that a given Xygra instance
//                 has, given the name of that instance in the netlist.
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/27/08
//-----------------------------------------------------------------------------
int DeviceInterface::xygraGetNumWindings(const string & deviceName)
{
  return devMgrPtr_->xygraGetNumWindings(deviceName);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::xygraGetCoilWindings
// Purpose       : Returns the number of windings that a given Xygra instance
//                 has in each coil, given the name of that instance in the
//                 netlist.
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/27/08
//-----------------------------------------------------------------------------
void DeviceInterface::xygraGetCoilWindings(const string & deviceName,
                                                vector<int> & cW)
{
  devMgrPtr_->xygraGetCoilWindings(deviceName,cW);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::xygraGetCoilNames
// Purpose       : Returns the names of each coil in a given Xygra instance
//                 given the name of that instance in the netlist.
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 09/29/08
//-----------------------------------------------------------------------------
void DeviceInterface::xygraGetCoilNames(const string & deviceName,
                                                vector<string> & cN)
{
  devMgrPtr_->xygraGetCoilNames(deviceName,cN);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::xygraSetConductances
// Purpose       : sets the conductance matrix on the given Xygra device
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/27/08
//-----------------------------------------------------------------------------
bool DeviceInterface::xygraSetConductances(const string & deviceName,
                                            const vector<vector<double> > &cM)
{
  return devMgrPtr_->xygraSetConductances(deviceName, cM);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::xygraSetK
// Purpose       : sets the K matrix on the given Xygra device
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/27/08
//-----------------------------------------------------------------------------
bool DeviceInterface::xygraSetK(const string & deviceName,
                                      const vector<vector<double> > &kM,
                                      const double t)
{
  return devMgrPtr_->xygraSetK(deviceName, kM,t);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::xygraSetSources
// Purpose       : sets the S vector on the given Xygra device
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/27/08
//-----------------------------------------------------------------------------
bool DeviceInterface::xygraSetSources(const string & deviceName,
                                            const vector<double> &sV,
                                            const double t)
{
  return devMgrPtr_->xygraSetSources(deviceName, sV, t);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::xygraGetVoltages
// Purpose       : retrieve nodal voltages from the given Xygra device
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 08/27/08
//-----------------------------------------------------------------------------
bool DeviceInterface::xygraGetVoltages(const string & deviceName,
                                             vector<double> &vN)
{
  return devMgrPtr_->xygraGetVoltages(deviceName, vN);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::getHomotopyBlockSize
// Purpose       : Returns the block size for block gainscale continuation.
// Special Notes : This call is virtual and should be overridden by
//                 the DeviceMgr
// Scope         : public
// Creator       : Roger Pawlowski, SNL, Computational Sciences
// Creation Date : 01/26/2005
//----------------------------------------------------------------------------
int DeviceInterface::getHomotopyBlockSize() const
{
  return devMgrPtr_->getHomotopyBlockSize();
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::addGlobalPar
// Purpose       : Add global parameter
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 11/18/2005
//----------------------------------------------------------------------------
void DeviceInterface::addGlobalPar(N_UTL_Param & par)
{
  devMgrPtr_->addGlobalPar(par);
}


//-----------------------------------------------------------------------------
// Function      : DeviceInterface::getGlobalPar
// Purpose       : get the value of a global parameter
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Electrical Systems Modeling
// Creation Date : 01/25/2013
//----------------------------------------------------------------------------
double DeviceInterface::getGlobalPar( const string & parName  ) const
{
  return devMgrPtr_->getGlobalPar( parName );
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
bool DeviceInterface::allDevsConverged()
{
  return devMgrPtr_->allDevsConverged();
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::innerDevsConverged
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 03/22/06
//-----------------------------------------------------------------------------
bool DeviceInterface::innerDevsConverged()
{
  return devMgrPtr_->innerDevsConverged();
}

#ifdef Xyce_EXTDEV
//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setupExternalDevices
// Purpose       : Setup external devices, especially critical for parallel
// Special Notes : needed for 2-level.
// Scope         : public
// Creator       : Robert Hoekstra, Elec. & MicroSystems Engring
// Creation Date : 03/11/06
//-----------------------------------------------------------------------------
void DeviceInterface::setupExternalDevices()
{
  devMgrPtr_->setupExternalDevices();
}
#endif

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::homotopyStepSuccess
// Purpose       :
// Special Notes : needed for 2-level.
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
void DeviceInterface::homotopyStepSuccess
      (const vector<string> & paramNames,
       const vector<double> & paramVals)
{
  return devMgrPtr_->homotopyStepSuccess (paramNames, paramVals);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::homotopyStepFailure
// Purpose       :
// Special Notes : needed for 2-level.
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/30/06
//-----------------------------------------------------------------------------
void DeviceInterface::homotopyStepFailure ()
{
  return devMgrPtr_->homotopyStepFailure ();
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::stepSuccess
// Purpose       :
// Special Notes : needed for 2-level.
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
void DeviceInterface::stepSuccess (int analysis)
{
  return devMgrPtr_->stepSuccess (analysis);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::stepFailure
// Purpose       :
// Special Notes : needed for 2-level.
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/12/06
//-----------------------------------------------------------------------------
void DeviceInterface::stepFailure (int analysis)
{
  return devMgrPtr_->stepFailure (analysis);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::acceptStep
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL
// Creation Date : 01/23/07
//-----------------------------------------------------------------------------
void DeviceInterface::acceptStep ()
{
  devMgrPtr_->acceptStep ();
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::getInitialQnorm
// Purpose       :
// Special Notes : needed for 2-level.
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/18/07
//-----------------------------------------------------------------------------
bool DeviceInterface::getInitialQnorm (vector<N_TIA_TwoLevelError> & tleVec )
{
  return devMgrPtr_->getInitialQnorm (tleVec);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::getInnerLoopErrorSums
// Purpose       :
// Special Notes : needed for 2-level.
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/15/06
//-----------------------------------------------------------------------------
bool DeviceInterface::getInnerLoopErrorSums (vector <N_TIA_TwoLevelError> & tleVec)
{
  return devMgrPtr_->getInnerLoopErrorSums (tleVec);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::updateStateArrays()
// Purpose       :
// Special Notes : needed for 2-level.
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/15/06
//-----------------------------------------------------------------------------
bool DeviceInterface::updateStateArrays()
{
  return devMgrPtr_->updateStateArrays();
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::startTimeStep
// Purpose       :
// Special Notes : needed for 2-level.
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
bool DeviceInterface::startTimeStep ()
{
  return devMgrPtr_->startTimeStep ();
}


//-----------------------------------------------------------------------------
// Function      : DeviceInterface::setExternalSolverState
// Purpose       :
// Special Notes : needed for 2-level.
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 03/20/06
//-----------------------------------------------------------------------------
void DeviceInterface::setExternalSolverState(const SolverState & ss)
{
  devMgrPtr_->setExternalSolverState (ss);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::updateSolverState
// Purpose       :
// Special Notes : needed for mixed signal
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 05/22/09
//-----------------------------------------------------------------------------
void DeviceInterface::updateSolverState ()
{
  bool tmp = devMgrPtr_->setupSolverInfo();
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::restartDataSize
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 06/26/2013
//-----------------------------------------------------------------------------
int DeviceInterface::restartDataSize(bool pack)
{
  return devMgrPtr_->restartDataSize(pack );
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::dumpRestartData
// Purpose       : Output restart data.
// Special Notes : 
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 06/26/2013
//-----------------------------------------------------------------------------
bool DeviceInterface::dumpRestartData
(char * buf, int bsize, int & pos, N_PDS_Comm * comm, bool pack )
{
  return devMgrPtr_->dumpRestartData(buf, bsize, pos, comm, pack );
}

//-----------------------------------------------------------------------------
// Function      : DeviceInterface::restoreRestartData
// Purpose       : Load restart data.
// Special Notes : 
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 06/26/2013
//-----------------------------------------------------------------------------
bool DeviceInterface::restoreRestartData
(char * buf, int bsize, int & pos, N_PDS_Comm * comm, bool pack )
{
  return devMgrPtr_->restoreRestartData(buf, bsize, pos, comm, pack);
}

} // namespace Device
} // namespace Xyce

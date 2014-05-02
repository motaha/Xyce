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
// Filename       : $RCSfile: N_DEV_DeviceInterface.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/05/02
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.145.2.1 $
//
// Revision Date  : $Date: 2014/02/26 20:16:30 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_DeviceInterface_h
#define Xyce_N_DEV_DeviceInterface_h

// ---------- Standard Includes ----------

#include <vector>
#include <map>
#include <string>

#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RefCountPtr;
using Teuchos::rcp;

#include <N_DEV_fwd.h>
#include <N_IO_fwd.h>
#include <N_UTL_fwd.h>
#include <N_ANP_fwd.h>

#include <N_UTL_Misc.h>
#include <N_UTL_BreakPoint.h>
#include <N_UTL_Param.h>

#include <N_PDS_Manager.h>

class N_TIA_TwoLevelError;

class N_NLS_Manager;

class N_LAS_System;
class N_LAS_Vector;
class N_LAS_Matrix;

class N_MPDE_Manager;

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : N_DEV_DeviceInterface
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/05/02
//-----------------------------------------------------------------------------

class DeviceInterface
{
  // functions:
public:
  static DeviceInterface * factory(N_IO_CmdParse & cp);

  ~DeviceInterface();

  // registration functions:
  bool registerLinearSystem(N_LAS_System * tmp_system_ptr);

  bool registerAnalysisInterface (N_ANP_AnalysisInterface * tmp_anaIntPtr);

  bool registerOutputMgr(N_IO_OutputMgr * tmp_outputMgrPtr);

  bool registerParallelMgr(N_PDS_Manager * tmp_pdsMgrPtr);

  bool registerNonlinearSolver (N_NLS_Manager * tmp_nlsMgrPtr);

  bool registerOptions(const Util::OptionBlock & OB);

  bool registerSensParams (const Util::OptionBlock & OB);

  bool registerICLoads( std::vector<std::pair<int,double> > * icLoads );

  bool registerPkgOptionsMgr( N_IO_PkgOptionsMgr *pkgOptPtr );

  const EntityTypeIdDeviceMap &getDeviceMap() const;
    
  // this function is called from the output manager to inform the
  // device package of the devices for which lead currents have been requested.
  // The device manager will take care of doing isolated F and Q loads
  // for these devices so the lead currents can be calculated
  bool setLeadCurrentRequests( const std::set<std::string> & deviceNames );

  // MPDE registrations
  std::vector<double> getFastSourcePeriod (std::vector<std::string>& sourceNames);
  std::vector<double> registerFastSources (std::vector<std::string> & sourceNames);
  void deRegisterFastSources (std::vector<std::string> & sourceNames);
  void deactivateSlowSources();
  void activateSlowSources();

  void setMPDEFlag( bool flagVal );
  void setBlockAnalysisFlag( bool flagVal );
  void setFastTime( double timeVal );

  // Initialization function, to be called after all registrations are
  // finished, and the linear system class is completely set up.
  bool initializeAll();

  void resetForStepAnalysis ();

  // Device accessor functions:
  bool addDeviceModel(const ModelBlock & MB);

  bool verifyDeviceInstance(InstanceBlock & IB);

  DeviceInstance * addDeviceInstance(InstanceBlock & IB);

  bool deleteDeviceInstance (const std::string & name);

  const std::map<std::string,int> & getDeviceCountMap ();
  void addDeviceToCount(const std::string & device_name);

  //    void  printOutLists();

  bool output ();
  bool finishOutput ();

  void  dotOpOutput ();

  // needed for parallel only:
  void setGlobalFlags ();

  // Load functions:
  bool loadDeviceMask();
  bool setInitialGuess   (N_LAS_Vector * solVectorPtr);

  bool setParam(std::string & name, double val);
  double getParamAndReduce(const std::string & name);
  bool getParamAndReduce(const std::string & name, double & val);
  double getParamNoReduce(const std::string & name) const;
  bool findParam(const std::string & name) const;

  bool getVsrcLIDs(std::string & srcName, int & li_Pos, int & li_Neg, int & li_Bra);

  bool updateSources();

  const DeviceSensitivities &getDeviceSensitivities() const;
  DeviceEntity *getDeviceEntity(const std::string &param) const;

  bool getLinearSystemFlag ();
  bool getVoltageLimiterFlag ();
  bool getPDESystemFlag ();

  // setup initial conditions on devices
  bool setICs (N_LAS_Vector * tmpSolVectorPtr,
               N_LAS_Vector * tmpCurrSolVectorPtr,
               N_LAS_Vector * tmpLastSolVectorPtr,
               N_LAS_Vector * tmpStaVectorPtr,
               N_LAS_Vector * tmpCurrStaVectorPtr,
               N_LAS_Vector * tmpLasStaVectorPtr,
               N_LAS_Vector * tmpStaDerivVectorPtr,
               N_LAS_Vector * tmpStoVectorPtr,
               N_LAS_Vector * tmpCurrStoVectorPtr,
               N_LAS_Vector * tmpLasStoVectorPtr,
               N_LAS_Vector * tmpQVectorPtr,
               N_LAS_Vector * tmpFVectorPtr,
               N_LAS_Vector * tmpdFdxdVpVectorPtr,
               N_LAS_Vector * tmpdQdxdVpVectorPtr);

  // Testing functions
  bool runParameterTests (std::string & deviceName);

  // time integration stuff:
  bool   getBreakPoints     ( std::vector<N_UTL_BreakPoint> & breakPointTimes );
  double getMaxTimeStepSize ();

  // Used by mixed signal:
  bool getDACDeviceNames(std::vector<std::string> & dacNames);
  bool getADCMap(std::map<std::string,std::map<std::string,double> >& ADCMap);
  bool updateTimeVoltagePairs(
     std::map<std::string,std::vector<std::pair<double,double> >* > const & timeVoltageUpdateMap);
  bool getTimeVoltagePairs(
     std::map<std::string,std::vector<std::pair<double,double> > > & timeVoltageUpdateMap);

  bool setADCWidths(std::map<std::string,int> const & ADCWidthMap);

  // Generic API calls (currently used only by Xygra)
  bool getDeviceNames(const std::string & deviceType,
                      std::vector<std::string> & deviceNames);

  int xygraGetNumNodes (const std::string & deviceName);
  int xygraGetNumWindings (const std::string & deviceName);
  void xygraGetCoilWindings (const std::string & deviceName,
                             std::vector<int> & cW);
  void xygraGetCoilNames (const std::string & deviceName,
                          std::vector<std::string> & cN);
  bool xygraSetConductances (const std::string & deviceName,
                             const std::vector<std::vector<double> > & cM);
  bool xygraSetK (const std::string & deviceName,
                  const std::vector<std::vector<double> > & kM,
                  const double t=0);
  bool xygraSetSources (const std::string & deviceName,
                        const std::vector<double> & sV,
                        const double t=0);
  bool xygraGetVoltages (const std::string & deviceName,
                         std::vector< double > & vN);

  // Two-level Newton and PDE-Continuation
  int  enablePDEContinuation ();
  bool disablePDEContinuation ();
  void getNumInterfaceNodes (std::vector<int> & numINodes);
  bool loadCouplingRHS (int iPDEDevice, int iElectrode, N_LAS_Vector * dfdvPtr);
  bool calcCouplingTerms (int iSubProblem, int iElectrode, const N_LAS_Vector * dxdvPtr);
  bool raiseDebugLevel (int increment);

  // load functions:
  bool loadDAEMatrices  (N_LAS_Vector * tmpSolVectorPtr,
                         N_LAS_Vector * tmpStaVectorPtr,
                         N_LAS_Vector * tmpStaDerivVectorPtr,
                         N_LAS_Vector * tmpStoVectorPtr,
                         N_LAS_Matrix * tmpdQdxMatrixPtr,
                         N_LAS_Matrix * tmpdFdxMatrixPtr);

  bool loadDAEVectors   (N_LAS_Vector * tmpSolVectorPtr,
                         N_LAS_Vector * tmpCurrSolVectorPtr,
                         N_LAS_Vector * tmpLastSolVectorPtr,
                         N_LAS_Vector * tmpStaVectorPtr,
                         N_LAS_Vector * tmpCurrStaVectorPtr,
                         N_LAS_Vector * tmpLasStaVectorPtr,
                         N_LAS_Vector * tmpStaDerivVectorPtr,
                         N_LAS_Vector * tmpStoVectorPtr,
                         N_LAS_Vector * tmpCurrStoVectorPtr,
                         N_LAS_Vector * tmpLasStoVectorPtr,
                         N_LAS_Vector * tmpStoLeadCurrQCompVectorPtr,
                         N_LAS_Vector * tmpQVectorPtr,
                         N_LAS_Vector * tmpFVectorPtr,
                         N_LAS_Vector * tmpdFdxdVpVectorPtr,
                         N_LAS_Vector * tmpdQdxdVpVectorPtr);

  bool updateState
  (N_LAS_Vector * nextSolVectorPtr,
   N_LAS_Vector * currSolVectorPtr,
   N_LAS_Vector * lastSolVectorPtr,
   N_LAS_Vector * nextStaVectorPtr,
   N_LAS_Vector * currStaVectorPtr,
   N_LAS_Vector * lastStaVectorPtr,
   N_LAS_Vector * nextStoVectorPtr,
   N_LAS_Vector * currStoVectorPtr,
   N_LAS_Vector * lastStoVectorPtr
   );

  bool loadBVectorsforAC (N_LAS_Vector * bVecRealPtr,
                          N_LAS_Vector * bVecImagPtr);

  bool getBMatrixEntriesforMOR(std::vector<int>& bMatEntriesVec, std::vector<int>& bMatPosEntriesVec);

  // voltlim doesn't work with MPDE, but does work for DCOP and the
  // various initial conditions used in MPDE.  Hence the need for these
  // functions.
  //    void setVoltageLimiterFlag ();
  void unsetVoltageLimiterFlag ();

  void setVoltageLimiterFlag (bool flagVal);

  int getHomotopyBlockSize() const;

  void addGlobalPar (Util::Param &);
  const double *findGlobalPar( const std::string & parName) const;
  double getGlobalPar( const std::string & parName ) const;

  // For convergence testing
  bool allDevsConverged();

  // For 2-level, need inner "solve" convergence.
  bool innerDevsConverged();

  // Functions needed for power node (2-level) algorithm):
#ifdef Xyce_EXTDEV
  void setupExternalDevices();
#endif

  void homotopyStepSuccess
  (const std::vector<std::string> & paramNames,
   const std::vector<double> & paramVals);

  void homotopyStepFailure ();

  void stepSuccess (int analysis);
  void stepFailure (int analysis);

  void acceptStep();

  bool getInitialQnorm (std::vector<N_TIA_TwoLevelError> & tleVec );
  bool getInnerLoopErrorSums (std::vector <N_TIA_TwoLevelError> & tleVec);

  bool updateStateArrays();
  bool startTimeStep ();
  void setExternalSolverState (const SolverState & ss);

  void updateSolverState ();

  int restartDataSize(bool pack);

  // Output restart data.
  bool dumpRestartData(char * buf, int bsize, int & pos, N_PDS_Comm * comm, bool pack );

  // Load restart data.
  bool restoreRestartData(char * buf, int bsize, int & pos, N_PDS_Comm * comm, bool pack );

protected:

private:
  DeviceInterface(N_IO_CmdParse & cp);
  DeviceInterface(const DeviceInterface &right);


public:

protected:

private:
  DeviceMgr     * devMgrPtr_;
};

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::DeviceInterface N_DEV_DeviceInterface;

#endif


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
// Revision Number: $Revision: 1.128.2.3 $
//
// Revision Date  : $Date: 2013/10/03 17:23:37 $
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

// ----------   Xyce Includes   ----------

#include <N_DEV_fwd.h>
#include <N_IO_fwd.h>
#include <N_UTL_fwd.h>
#ifdef Xyce_RAD_MODELS
#include <N_DEV_ExtendedModelTypes.h>
#else
#include <N_DEV_ModelTypes.h>
#endif
#include <N_UTL_Misc.h>
#include <N_UTL_BreakPoint.h>
#include <N_UTL_Param.h>
// to eliminate RCP compiler warnings, N_PDS_Manager header is being put here.
#include <N_PDS_Manager.h>

// ---------- Forward Declarations ----------

class N_ANP_AnalysisInterface;
class N_TIA_TwoLevelError;
class N_IO_PkgOptionsMgr;

class N_NLS_Manager;

class N_LAS_System;
class N_LAS_Vector;
class N_LAS_Matrix;

class N_MPDE_Manager;

class N_IO_CmdParse;

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

    // This function returns a base class device pointer.  The pointer should
    // point to a "type" level device class, which is (usual) two levels of
    // inheritance below the base device class.
    Device * returnDevicePtr(int Index);

    // These functions are used by the parser.
    Device * createDeviceByNetlistDeviceType(const std::string &name, const int level);

    int getNumSupportedDevices ();

    // registration functions:
    bool registerLinearSystem(N_LAS_System * tmp_system_ptr);

    bool registerAnalysisInterface (N_ANP_AnalysisInterface * tmp_anaIntPtr);

    bool registerOutputMgr(N_IO_OutputMgr * tmp_outputMgrPtr);

    bool registerParallelMgr(N_PDS_Manager * tmp_pdsMgrPtr);

    bool registerNonlinearSolver (N_NLS_Manager * tmp_nlsMgrPtr);

    bool registerOptions(const N_UTL_OptionBlock & OB);

    bool registerSensParams (const N_UTL_OptionBlock & OB);

    bool registerICLoads( vector< pair<int,double> > * icLoads );

    bool registerPkgOptionsMgr( RCP<N_IO_PkgOptionsMgr> pkgOptPtr );

    // this function is called from the output manager to inform the
    // device package of the devices for which lead currents have been requested.
    // The device manager will take care of doing isolated F and Q loads
    // for these devices so the lead currents can be calculated
    bool setLeadCurrentRequests( const set<string> & deviceNames );

    // MPDE registrations
    vector<double> getFastSourcePeriod (std::vector<std::string>& sourceNames);
    vector<double> registerFastSources (std::vector<std::string> & sourceNames);
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
    bool flushDevices();

    bool addDeviceModel(const ModelBlock & MB);

    bool verifyDeviceInstance(InstanceBlock & IB);

    DeviceInstance * addDeviceInstance(InstanceBlock & IB);

    bool deleteDeviceInstance (const std::string & name);

    const map<string,int> & getDeviceCountMap ();

    void  printOutLists();

    bool output ();
    bool finishOutput ();

    void  dotOpOutput ();

    // needed for parallel only:
    void setGlobalFlags ();

    // Load functions:
    bool loadDeviceMask();
    bool setInitialGuess   (N_LAS_Vector * solVectorPtr);

    bool setParam          (string & name, double val);
    double getParam        (const string & name);
    bool   getParam        (const string & name, double & val);
    bool   getVsrcLIDs     (string & srcName, int & li_Pos, int & li_Neg, int & li_Bra);

    bool updateSources();

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
    bool runParameterTests (string & deviceName);

    // time integration stuff:
    bool   getBreakPoints     ( vector<N_UTL_BreakPoint> & breakPointTimes );
    double getMaxTimeStepSize ();

   // Used by mixed signal:
  bool getDACDeviceNames(std::vector<std::string> & dacNames);
    bool getADCMap(map<string,map<string,double> >& ADCMap);
    bool updateTimeVoltagePairs(
       map<string,vector<pair<double,double> >* > const & timeVoltageUpdateMap);
    bool getTimeVoltagePairs(
       map<string,vector<pair<double,double> > > & timeVoltageUpdateMap);

    bool setADCWidths(map<string,int> const & ADCWidthMap);

    // Generic API calls (currently used only by Xygra)
  bool getDeviceNames(const std::string & deviceType,
                      std::vector<std::string> & deviceNames);

  int xygraGetNumNodes (const std::string & deviceName);
  int xygraGetNumWindings (const std::string & deviceName);
  void xygraGetCoilWindings (const std::string & deviceName,
                              vector<int> & cW);
  void xygraGetCoilNames (const std::string & deviceName,
                          std::vector<std::string> & cN);
  bool xygraSetConductances (const std::string & deviceName,
                               const vector< vector<double> > & cM);
  bool xygraSetK (const std::string & deviceName,
                    const vector< vector<double> > & kM,
                    const double t=0);
  bool xygraSetSources (const std::string & deviceName,
                          const vector<double> & sV,
                          const double t=0);
  bool xygraGetVoltages (const std::string & deviceName,
                           vector< double > & vN);

    // Two-level Newton and PDE-Continuation
    int  enablePDEContinuation ();
    bool disablePDEContinuation ();
    void getNumInterfaceNodes (vector<int> & numINodes);
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
    void setVoltageLimiterFlag ();
    void unsetVoltageLimiterFlag ();

    int getHomotopyBlockSize() const;

    void addGlobalPar (N_UTL_Param &);
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
       const vector<double> & paramVals);

    void homotopyStepFailure ();

    void stepSuccess (int analysis);
    void stepFailure (int analysis);

    void acceptStep();

    bool getInitialQnorm (vector<N_TIA_TwoLevelError> & tleVec );
    bool getInnerLoopErrorSums (vector <N_TIA_TwoLevelError> & tleVec);

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


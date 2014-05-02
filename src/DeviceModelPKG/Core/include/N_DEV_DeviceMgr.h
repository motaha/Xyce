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
// Filename       : $RCSfile: N_DEV_DeviceMgr.h,v $
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
// Revision Number: $Revision: 1.350.2.1 $
//
// Revision Date  : $Date: 2014/02/26 20:16:30 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_DeviceMgr_h
#define Xyce_N_DEV_DeviceMgr_h

#include <map>
#include <set>
#include <string>
#include <vector>

#include <N_DEV_fwd.h>
#include <N_IO_fwd.h>
#include <N_UTL_fwd.h>
#include <N_ANP_fwd.h>

#include <N_DEV_DeviceOptions.h>
#include <N_DEV_DeviceSupport.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_Pars.h>
#include <N_DEV_SolverState.h>
#include <N_IO_PkgOptionsMgr.h>

class N_LAS_Matrix;
class N_LAS_System;
class N_LAS_Vector;
class N_NLS_Manager;
class N_PDS_Manager;
class N_TIA_TwoLevelError;

namespace Xyce {
namespace Device {

void setupIOName(const std::string & name, std::string & outputName);

//-----------------------------------------------------------------------------
// Class         : DeviceMgr
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
class DeviceMgr
{
public:
  typedef std::vector<Device *> DeviceVector;
  typedef std::vector<DeviceInstance *> InstanceVector;
  typedef std::vector<DeviceModel *> ModelVector;
  typedef std::map<std::string, DeviceEntity *, LessNoCase> DeviceEntityMap;
  typedef std::map<std::string, ModelTypeId, LessNoCase> ModelTypeNameModelTypeIdMap;

  static DeviceMgr * factory(IO::CmdParse & cp);

private:
  DeviceMgr(IO::CmdParse & cp);                       ///< Only the factory can create a device manager

  DeviceMgr(const DeviceMgr &);                       ///< No copying
  DeviceMgr &operator=(const DeviceMgr &);            ///< No assignment

public:
  ~DeviceMgr();

  // registration functions:
  bool registerLinearSystem(N_LAS_System * tmp_system_ptr);

  bool registerAnalysisInterface(N_ANP_AnalysisInterface * tmp_anaIntPtr);

  bool registerOutputMgr (IO::OutputMgr * tmp_outputMgrPtr);

  bool registerParallelMgr(N_PDS_Manager * tmp_pdsMgrPtr);

  bool registerNonlinearSolver (N_NLS_Manager * tmp_nlsMgrPtr);

  bool registerICLoads( std::vector<std::pair<int,double> > * icLoads );

  bool registerPkgOptionsMgr( IO::PkgOptionsMgr *pkgOptPtr );

  // this function is called from the output manager (through the
  // device interface) to inform the device package of the devices for
  // which lead currents have been requested. The device manager will
  // take care of doing isolated F and Q loads for these devices so
  // the lead currents can be calculated
  bool setLeadCurrentRequests(const std::set<std::string> & deviceNames );

  // MPDE related registrations:
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

  int getHomotopyBlockSize() const;

  //    void  printOutLists();

  bool output ();
  bool finishOutput ();

  void  dotOpOutput ();

  // Load functions:
  bool setInitialGuess (N_LAS_Vector * solVectorPtr);
  bool loadDeviceMask();

  void debugOutput1();
  void debugOutput2();

  bool setParam(std::string & name, double val);
  double getParamAndReduce(const std::string & name);
  bool getParamAndReduce(const std::string & name, double & val);
  double getParamNoReduce(const std::string & name) const;
  bool findParam(const std::string & name) const;
  bool   getVsrcLIDs       (std::string & srcName, int & li_Pos, int & li_Neg, int & li_Bra);

  bool   updateTemperature (double val);

  bool updateSources();

  const EntityTypeIdDeviceMap &getDeviceMap() const 
  {
    return deviceMap_;
  }

  bool resetRHSLoadFlags (int index);

  const DeviceSensitivities &getDeviceSensitivities() const 
  {
    return *devSensPtr_;
  }

  const N_NLS_Manager &getNlsMgrPtr() const 
  {
    return *nlsMgrPtr_;
  }

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
               N_LAS_Vector * tmpLastStoVectorPtr,
               N_LAS_Vector * tmpQVectorPtr,
               N_LAS_Vector * tmpFVectorPtr,
               N_LAS_Vector * tmpdFdxdVpVectorPtr,
               N_LAS_Vector * tmpdQdxdVpVectorPtr);

  // Testing functions
  bool runParameterTests (std::string & deviceName);

  // time integration stuff:
  bool   getBreakPoints     ( std::vector<Util::BreakPoint> & breakPointTimes );
  double getMaxTimeStepSize ();
  void  declareCurrentStepAsBreakpoint();

  // Used by mixed signal:
  bool getDACDeviceNames(std::vector<std::string> & dacNames);
  bool getADCMap(std::map<std::string,std::map<std::string,double> >& ADCMap);
  bool updateTimeVoltagePairs(
     std::map<std::string,std::vector<std::pair<double,double> >* > const & timeVoltageUpdateMap);
  bool getTimeVoltagePairs(
     std::map<std::string,std::vector<std::pair<double,double> > > & timeVoltageUpdateMap);

  bool setADCWidths(std::map<std::string,int> const & ADCWidthMap);

  // Generic API calls (so far used only by Xygra interface)
  bool getDeviceNames(const std::string & model_type_name, std::vector<std::string> & deviceNames);

  int xygraGetNumNodes (const std::string & deviceName);
  int xygraGetNumWindings (const std::string & deviceName);
  void xygraGetCoilWindings (const std::string & deviceName,
                             std::vector<int> & cW);
  void xygraGetCoilNames (const std::string & deviceName,
                          std::vector<std::string> & cN);
  bool xygraSetConductances (const std::string & deviceName,
                             const std::vector<std::vector<double> > &cM );
  bool xygraSetSources (const std::string & deviceName,
                        const std::vector<double> &sV,
                        const double t=0);
  bool xygraSetK (const std::string & deviceName,
                  const std::vector<std::vector<double> > &kM,
                  const double t=0);
  bool xygraGetVoltages (const std::string & deviceName,
                         std::vector<double> & vN );

  // two-level newton and pde-continuation
  int  enablePDEContinuation ();
  bool disablePDEContinuation ();
  void getNumInterfaceNodes (std::vector<int> & numInterfaceNodes);
  bool loadCouplingRHS (int iPDEDevice, int iElectrode, N_LAS_Vector * dfdvPtr);
  bool calcCouplingTerms (int iSubProblem, int iElectrode, const N_LAS_Vector * dxdvPtr);
  bool raiseDebugLevel (int increment);

  bool calcPDESubProblemInfo ();

  // load functions:
  bool loadDAEMatrices  (N_LAS_Vector * tmpSolVectorPtr,
                         N_LAS_Vector * tmpStaVectorPtr,
                         N_LAS_Vector * tmpStaDerivVectorPtr,
                         N_LAS_Vector * tmpStoVectorPtr,
                         N_LAS_Matrix * tmpdQdxMatrixPtr,
                         N_LAS_Matrix * tmpdFdxMatrixPtr);

  bool loadDAEVectors   (N_LAS_Vector * tmpNextSolVectorPtr,
                         N_LAS_Vector * tmpCurrSolVectorPtr,
                         N_LAS_Vector * tmpLastSolVectorPtr,
                         N_LAS_Vector * tmpNextStaVectorPtr,
                         N_LAS_Vector * tmpCurrStaVectorPtr,
                         N_LAS_Vector * tmpLastStaVectorPtr,
                         N_LAS_Vector * tmpStaDerivVectorPtr,
                         N_LAS_Vector * tmpNextStoVectorPtr,
                         N_LAS_Vector * tmpCurrStoVectorPtr,
                         N_LAS_Vector * tmpLastStoVectorPtr,
                         N_LAS_Vector * tmpStoLeadCurrQCompVectorPtr,
                         N_LAS_Vector * tmpQVectorPtr,
                         N_LAS_Vector * tmpFVectorPtr,
                         N_LAS_Vector * tmpdFdxdVpVectorPtr,
                         N_LAS_Vector * tmpdQdxdVpVectorPtr);

  bool updateState      (N_LAS_Vector * nextSolVectorPtr,
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

  void setVoltageLimiterFlag ( bool flagVal );

  void addGlobalPar(Util::Param &);
  const double *findGlobalPar( const std::string & parName) const;
  double getGlobalPar( const std::string & parName ) const;


  // structs related to options registration
  struct DeviceMgr_OptionsReg : public IO::PkgOptionsReg
  {
    DeviceMgr_OptionsReg(DeviceMgr * mgr) : devMgr(mgr) {}

    bool operator()(const Util::OptionBlock & options)
    {return devMgr->registerOptions(options);}

    DeviceMgr * devMgr;
  };

  struct DeviceMgr_SensOptionsReg : public IO::PkgOptionsReg
  {
    DeviceMgr_SensOptionsReg(DeviceMgr * mgr) : devMgr(mgr) {}

    bool operator()(const Util::OptionBlock & options)
    {return devMgr->registerSensParams(options);}

    DeviceMgr * devMgr;
  };

  struct DeviceMgr_TimeOptionsReg : public IO::PkgOptionsReg
  {
    DeviceMgr_TimeOptionsReg(DeviceMgr * mgr) : devMgr(mgr) {}

    bool operator()(const Util::OptionBlock & options)
    {return devMgr->registerTimeOptions(options);}

    DeviceMgr * devMgr;
  };

  // .TRAN
  struct DeviceMgr_TransAnalysisReg : public IO::PkgOptionsReg
  {
    DeviceMgr_TransAnalysisReg(DeviceMgr * mgr) : Mgr(mgr) {}

    bool operator()(const Util::OptionBlock & options)
    {return Mgr->setTranAnalysisParams(options);}

    DeviceMgr * Mgr;
  };

  // .DC
  struct DeviceMgr_DCAnalysisReg : public IO::PkgOptionsReg
  {
    DeviceMgr_DCAnalysisReg(DeviceMgr * mgr) : Mgr(mgr) {}

    bool operator()(const Util::OptionBlock & options)
    {return Mgr->setDCAnalysisParams(options);}

    DeviceMgr * Mgr;
  };

  // .OP
  struct DeviceMgr_OPAnalysisReg : public IO::PkgOptionsReg
  {
    DeviceMgr_OPAnalysisReg(DeviceMgr * mgr) : Mgr(mgr) {}

    bool operator()(const Util::OptionBlock & options)
    {return Mgr->setOPAnalysisParams(options);}

    DeviceMgr * Mgr;
  };

  // .STEP
  struct DeviceMgr_STEPAnalysisReg : public IO::PkgOptionsReg
  {
    DeviceMgr_STEPAnalysisReg(DeviceMgr * mgr) : Mgr(mgr) {}

    bool operator()(const Util::OptionBlock & options)
    {return Mgr->setSTEPAnalysisParams(options);}

    DeviceMgr * Mgr;
  };

  // .MPDE
  struct DeviceMgr_MPDE_AnalysisReg : public IO::PkgOptionsReg
  {
    DeviceMgr_MPDE_AnalysisReg(DeviceMgr * devMgr) : devMgr_(devMgr) {}

    bool operator()(const Util::OptionBlock & options)
    {return devMgr_->setMPDEAnalysisParams(options);}

    DeviceMgr * devMgr_;
  };

  // .HB
  struct DeviceMgr_HB_AnalysisReg : public IO::PkgOptionsReg
  {
    DeviceMgr_HB_AnalysisReg(DeviceMgr * devMgr) : devMgr_(devMgr) {}

    bool operator()(const Util::OptionBlock & options)
    {return devMgr_->setHBAnalysisParams(options);}

    DeviceMgr * devMgr_;
  };

  // .AC
  struct DeviceMgr_AC_AnalysisReg : public IO::PkgOptionsReg
  {
    DeviceMgr_AC_AnalysisReg(DeviceMgr * devMgr) : devMgr_(devMgr) {}

    bool operator()(const Util::OptionBlock & options)
    {return devMgr_->setACAnalysisParams(options);}

    DeviceMgr * devMgr_;
  };

  // .MOR
  struct DeviceMgr_MOR_AnalysisReg : public IO::PkgOptionsReg
  {
    DeviceMgr_MOR_AnalysisReg(DeviceMgr * devMgr) : devMgr_(devMgr) {}

    bool operator()(const Util::OptionBlock & options)
    {return devMgr_->setMORAnalysisParams(options);}

    DeviceMgr * devMgr_;
  };

  // functions related to options registration
  bool registerOptions(const Util::OptionBlock & OB) {
    bool result = devOptions_.registerOptions(OB);
    devOptions_.applyCmdLineOptions(commandLine_);
    return result;
  }
  bool registerSensParams(const Util::OptionBlock & OB);
  bool registerTimeOptions(const Util::OptionBlock & OB){return true;}
  bool setTranAnalysisParams(const Util::OptionBlock & OB);
  bool setDCAnalysisParams(const Util::OptionBlock & OB);
  bool setOPAnalysisParams(const Util::OptionBlock & OB);
  bool setSTEPAnalysisParams(const Util::OptionBlock & OB);
  bool setMPDEAnalysisParams(const Util::OptionBlock & OB);
  bool setHBAnalysisParams(const Util::OptionBlock & OB);
  bool setACAnalysisParams(const Util::OptionBlock & OB);
  bool setMORAnalysisParams(const Util::OptionBlock & OB);

  // convergence:  allow devices to signal back to the solvers that
  // they've played some game that invalidates normal convergence tests,
  // and so the solution should be considered unconverged no matter how
  // small the various norms are.
  bool allDevsConverged();

  // Similar to allDevsConverged, but specific to "inner" devices of
  // 2-level solves.  They are handled slightly differently.
  bool innerDevsConverged();

  // Functions needed for power node (2-level) algorithm):
#ifdef Xyce_EXTDEV
  // for the parallel case, we need to give all the processors a copy
  // of the device so all the parallel synchronized calls such are
  // called by all processors together.
  bool setupExternalDevices();
#endif

  const std::map<std::string,int> &getDeviceCountMap() {
    return localDeviceCountMap_;
  }

  void addDeviceToCount(const std::string & device_name)
  {
    localDeviceCountMap_[device_name]++;
  }

  void addDeviceEntity(const std::string &param, DeviceEntity *entity) 
  {
    nameDevEntityMap_[param] = entity;
  }
    
  DeviceEntity *getDeviceEntity(const std::string &param) const;

  void homotopyStepSuccess
  (const std::vector<std::string> & paramNames,
   const std::vector<double> & paramVals);

  void homotopyStepFailure ();

  void stepSuccess (int analysis);
  void stepFailure (int analysis);

  void acceptStep();

  bool getInitialQnorm (std::vector<N_TIA_TwoLevelError> & tleVec );
  bool getInnerLoopErrorSums (std::vector<N_TIA_TwoLevelError> & tleVec);

  bool updateStateArrays();
  bool startTimeStep ();
  void setExternalSolverState (const SolverState & ss);

  int restartDataSize(bool pack);

  // Output restart data.
  bool dumpRestartData(char * buf, int bsize, int & pos, N_PDS_Comm * comm, bool pack );

  // Load restart data.
  bool restoreRestartData(char * buf, int bsize, int & pos, N_PDS_Comm * comm, bool pack );

private:
  Device *getDevice(EntityTypeId model_type_id) 
  {
    EntityTypeIdDeviceMap::iterator it = deviceMap_.find(model_type_id);
    return it == deviceMap_.end() ? 0 : (*it).second;
  }

  Device &getDeviceByModelType(const EntityTypeId model_type);
  EntityTypeId getModelGroup(const std::string &device_type_name);

  bool setupSolverInfo_ ();
  bool setupRawVectorPointers_ ();
  bool setupRawMatrixPointers_ ();


#ifdef Xyce_EXTDEV
  void setUpPassThroughParamsMap_();
#endif

  bool updateIntermediateVars_();
  bool updatePrimaryState_();
  bool updateSecondaryState_();

  bool updateDependentParameters_();

#ifdef Xyce_EXTDEV
  // Do the actual solve/calculation for the external devices
  void updateExternalDevices_();
  // add external devices for processors that don't own it
  ExternDevice::Instance * addExtDeviceInstance_(InstanceBlock & IB);
#endif

  Xygra::Instance * getXygraInstancePtr_(const std::string & deviceName);

  // attributes:

private:
  IO::CmdParse & commandLine_;                        ///< Command line

  // user-defined options:
  DeviceOptions devOptions_;

  EntityTypeIdDeviceMap   deviceMap_;

  DeviceSupport devSupport;  // doesn't appear to be used in the code at all.
  // an object owned by DeviceInstance looks like it's
  // used in all cases.

  bool allDevsConverged_;
  bool sensFlag_;
  bool linearSystemFlag_;

  bool firstDependent;

  bool jacobianLoadCalledBefore_;

  bool entityMapDone_;
  bool parameterChanged_;

  bool breakPointInstancesInitialized;

  double timeParamsProcessed_;

  ExternData externData_;

  N_LAS_Vector       *  numJacSolVectorPtr_;
  N_LAS_Vector       *  numJacStaVectorPtr_;
  N_LAS_Vector       *  numJacStoVectorPtr_;

  N_LAS_Vector       *  diagonalVectorPtr_;

  N_LAS_System        * lasSysPtr_;

  N_ANP_AnalysisInterface       * anaIntPtr_;

  IO::OutputMgr      * outputMgrPtr_;

  N_PDS_Manager       * pdsMgrPtr_;

  N_NLS_Manager       * nlsMgrPtr_;

  DeviceMgr     * DevMgrPtr_;

  IO::PkgOptionsMgr *pkgOptMgrPtr_;

  std::vector< std::pair<int,double> > *                icLoads_;

  ModelTypeNameModelTypeIdMap         modelTypeMap_;          ///< Model type name to model
  ModelTypeNameModelTypeIdMap         modelGroupMap_;         ///< Model type name to model group

  std::map<std::string,int>             localDeviceCountMap_;

  std::map<std::string, Xygra::Instance *> xygraPtrMap_;

  std::multimap<int,DeviceInstance *> solDevInstMap_;

  DeviceVector devicePtrVec_;
  DeviceVector pdeDevicePtrVec_;
  DeviceVector nonPdeDevicePtrVec_;

  InstanceVector instancePtrVec_;
  InstanceVector bpInstancePtrVec_; // instances with breakpoints functions
  InstanceVector pdeInstancePtrVec_;
  InstanceVector nonPdeInstancePtrVec_;
  InstanceVector mosfetInstancePtrVec_;
  InstanceVector vsrcInstancePtrVec_;
  InstanceVector bjtInstancePtrVec_;
  std::map<std::string, Vsrc::Instance *> vsrcInstancePtrMap_;

  InstanceVector plotFileInstancePtrVec_;

  std::map<std::string,SourceInstance *> indepSourcePtrMap_;
  std::vector<SourceInstance *> indepSourceInstancePtrVec_;

  // this is used to store the contents of the indepSourceInstancePtrVec_
  // during an mpde initialization where we'll remove slow sources from
  // the that vector so that they don't get updated
  std::vector<SourceInstance *> indepSourceInstanceBackupPtrVec_;

  std::set<std::string> devicesNeedingLeadCurrentLoads_;

#ifdef Xyce_EXTDEV
  std::vector<ExternDevice::Instance *> extDevInstancePtrVec_;
  std::vector<InstanceBlock *> extDevIBPtrVec_;

  std::map<std::string, int> passThroughParamsMap_;
#endif

  mutable DeviceEntityMap     nameDevEntityMap_;

  // vector of pointers to devices under test.
  InstanceVector testJacDevicePtrVec_;

  ModelVector modelPtrVec_;
  ModelVector mosfetModelPtrVec_;
  ModelVector bsim3ModelPtrVec_;
  ModelVector bsim4ModelPtrVec_;
  ModelVector bsimsoiModelPtrVec_;
  ModelVector bjtModelPtrVec_;
  ModelVector diodeModelPtrVec_;

  std::vector<DeviceEntity *> dependentPtrVec_;

  std::vector<int> numInterfaceNodes_;
  int numPDEDevices_;
  bool calledBeforeCSPI;

  int numThreads_;
  bool multiThreading_;

  // real time solver data:
  SolverState solState_;
  SolverState solStateExternal_;
  bool externalStateFlag_;

  // sensitivities:
  DeviceSensitivities * devSensPtr_;

  // MPDE fast source list
  std::vector<std::string> fastSourceNames_;

  // Device mask flag:
  bool nonTrivialDeviceMaskFlag;

  // .OP output flags.
  bool dotOpOutputFlag;

  // temporary jacobian load structures:
  MatrixLoadData mlData;

  // solution variable names vector:
  std::vector<std::string> nameVec_;

public:
  // needed for parallel only:
  void setGlobalFlags ();

  bool setupSolverInfo() 
  {
    return setupSolverInfo_();
  }
};

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::registerLinearSystem
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/13/00
//-----------------------------------------------------------------------------
inline bool DeviceMgr::registerLinearSystem (N_LAS_System * tmp_system_ptr)
{
  lasSysPtr_ = tmp_system_ptr;

  if (lasSysPtr_ != 0) return true;
  else                    return false;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::registerAnalysisInterface
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
inline bool DeviceMgr::registerAnalysisInterface (N_ANP_AnalysisInterface * tmp_anaIntPtr)
{
  anaIntPtr_              = tmp_anaIntPtr;

  if (anaIntPtr_ != 0) return true;
  else                    return false;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::registerOutputMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 8/15/06
//-----------------------------------------------------------------------------
inline bool DeviceMgr::registerOutputMgr (IO::OutputMgr * tmp_outputMgrPtr)
{
  outputMgrPtr_              = tmp_outputMgrPtr;

  if (outputMgrPtr_ != 0) return true;
  else                    return false;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::registerParallelMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/29/01
//-----------------------------------------------------------------------------
inline bool DeviceMgr::registerParallelMgr (N_PDS_Manager * tmp_pdsMgrPtr )
{
  pdsMgrPtr_              = tmp_pdsMgrPtr;

  if (pdsMgrPtr_ != 0 ) return true;
  else                     return false;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::registerNonlinearSolver
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 2/23/01
//-----------------------------------------------------------------------------
inline bool DeviceMgr::registerNonlinearSolver (N_NLS_Manager * tmp_nlsMgrPtr)
{
  nlsMgrPtr_              = tmp_nlsMgrPtr;

  if (nlsMgrPtr_ != 0) return true;
  else                    return false;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setVoltageLimiterFlag ()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/09/04
//-----------------------------------------------------------------------------
//inline void DeviceMgr::setVoltageLimiterFlag ()
//{
//  devOptions_.voltageLimiterFlag = true;
//}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::unsetVoltageLimiterFlag ()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/09/04
//-----------------------------------------------------------------------------
inline void DeviceMgr::unsetVoltageLimiterFlag ()
{
  devOptions_.voltageLimiterFlag = false;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::setVoltageLimiterFlag ()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 02/18/14
//-----------------------------------------------------------------------------
inline void DeviceMgr::setVoltageLimiterFlag ( bool flagVal )
{
  devOptions_.voltageLimiterFlag = flagVal;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::registerICLoads
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 04/03/01
//-----------------------------------------------------------------------------
inline bool DeviceMgr::registerICLoads( std::vector< std::pair<int,double> > * icLoads )
{
  return (icLoads_ = icLoads) != 0;
}

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::DeviceMgr N_DEV_DeviceMgr;

#endif // Xyce_N_DEV_DeviceMgr_h

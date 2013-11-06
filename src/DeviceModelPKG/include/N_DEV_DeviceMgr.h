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
// Revision Number: $Revision: 1.323.2.4 $
//
// Revision Date  : $Date: 2013/10/03 17:23:37 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_DeviceMgr_h
#define Xyce_N_DEV_DeviceMgr_h

// ---------- Standard Includes ----------

#include <vector>
#include <string>
#include <map>
#include <set>

#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RefCountPtr;
using Teuchos::rcp;

// ----------   Xyce Includes   ----------
#include <N_UTL_Misc.h>
#include <N_UTL_NoCase.h>
#include <N_DEV_fwd.h>
#include <N_IO_fwd.h>
#include <N_UTL_fwd.h>

#include <N_UTL_OptionBlock.h>

#include <N_IO_PkgOptionsMgr.h>

#include <N_DEV_Device.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_DeviceBld.h>
#include <N_DEV_DeviceEntity.h>
#include <N_DEV_DeviceSupport.h>

#include <N_DEV_DeviceSensitivities.h>
#include <N_UTL_BreakPoint.h>
#include <N_MPDE_Manager.h>

// ---------- Forward Declarations ----------
class N_ANP_AnalysisInterface;
class N_NLS_Manager;
class N_TIA_TwoLevelError;

class N_LAS_System;
class N_LAS_Matrix;
class N_LAS_MultiVector;
class N_LAS_Vector;

class N_PDS_Manager;

class N_ERH_ErrorMgr;

class colData;
class valData;

class N_IO_CmdParse;

#ifdef Xyce_RAD_MODELS
#include <N_DEV_ExtendedModelTypes.h>
#include <N_DEV_ExtendedModelLevels.h>
#else
#include <N_DEV_ModelTypes.h>
#endif

namespace Xyce {
namespace Device {

// ----------------------------------------------------------------------------
// Function      : setupIOName
//
// Purpose       : This function takes the device instance name and creates
//                 an appropriate "outputName" to be used for file outputs.
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/26/2013
// ----------------------------------------------------------------------------
inline void setupIOName (const string & name, string & outputName)
{
  // get rid of the y-device prefix.
  string firstString("Y%");
  string secondString("%");
  std::string::size_type pos1 = name.find_first_of(firstString);
  std::string::size_type pos2 = name.find_last_of(secondString);

  if (pos1 != std::string::npos && pos2 != std::string::npos)
  {
    string tmp1 = "";
    if (pos1 > 0) tmp1 = name.substr(0,pos1);
    string tmp2 = name.substr(pos2+1, name.length()-1);
    outputName = tmp1 + tmp2;
  }
  else
  {
    outputName = name;
  }

  return;
}

//-----------------------------------------------------------------------------
// Class         : DeviceMgr
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
using Xyce::Device::ModelType;

class DeviceMgr
{
  public:
    typedef std::map<std::string, int, LessNoCase> DeviceIndexMap;
    typedef std::map<std::string, int, LessNoCase> DeviceModelTypeMap;
    typedef std::vector<Device *> DeviceVector;

    static DeviceMgr * factory(N_IO_CmdParse & cp);

    ~DeviceMgr();

    const Device &getDeviceArray(int i) const {
      return *deviceArray_[i];
    }
    Device &getDeviceArray(int i) {
      return *deviceArray_[i];
    }

    const int &getDeviceAllocFlag(int i) const {
      return deviceAllocFlag_[i];
    }

    const N_NLS_Manager &getNlsMgrPtr() const {
      return *nlsMgrPtr_;
    }

    // This function returns a base class device pointer.  The pointer should
    // point to a "type" level device class, which is (usual) two levels of
    // inheritance below the base device class.
    Device * returnDevicePtr(int Index);

    // These functions are used by the parser.
    Device * createDeviceByNetlistDeviceType(const std::string &name, const int level);

    template <class T>
    ParametricData<T> *getParametricData(const std::string &name, const int level) {
      return 0;
    }

    int getNumSupportedDevices ();

    // registration functions:
    bool registerLinearSystem(N_LAS_System * tmp_system_ptr);

    bool registerAnalysisInterface(N_ANP_AnalysisInterface * tmp_anaIntPtr);

    bool registerOutputMgr (N_IO_OutputMgr * tmp_outputMgrPtr);

    bool registerParallelMgr(N_PDS_Manager * tmp_pdsMgrPtr);

    bool registerNonlinearSolver (N_NLS_Manager * tmp_nlsMgrPtr);

    const map<string,int> & getDeviceCountMap ();

    bool registerICLoads( vector< pair<int,double> > * icLoads );

    bool registerPkgOptionsMgr( RCP<N_IO_PkgOptionsMgr> pkgOptPtr );

    // this function is called from the output manager (through the
    // device interface) to inform the device package of the devices for
    // which lead currents have been requested. The device manager will
    // take care of doing isolated F and Q loads for these devices so
    // the lead currents can be calculated
    bool setLeadCurrentRequests(const std::set<std::string> & deviceNames );

    // MPDE related registrations:
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

    int getHomotopyBlockSize() const;

    void  printOutLists();

    bool output ();
    bool finishOutput ();

    void  dotOpOutput ();

    // Load functions:
    bool setInitialGuess (N_LAS_Vector * solVectorPtr);
    bool loadDeviceMask();

#ifdef Xyce_DEBUG_DEVICE
    void debugOutput1();
    void debugOutput2();
#endif

    bool   setParam          (string & name, double val);
    double getParam          (const string & name);
    bool   getParam          (const string & name, double & val);
    bool   getVsrcLIDs       (string & srcName, int & li_Pos, int & li_Neg, int & li_Bra);

    bool   updateTemperature (double val);

    bool updateSources();
    bool resetRHSLoadFlags (int index);

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
    bool runParameterTests (string & deviceName);

    // time integration stuff:
    bool   getBreakPoints     ( vector<N_UTL_BreakPoint> & breakPointTimes );
    double getMaxTimeStepSize ();
    void  declareCurrentStepAsBreakpoint();

    // Used by mixed signal:
    bool getDACDeviceNames(vector<string> & dacNames);
    bool getADCMap(map<string,map<string,double> >& ADCMap);
    bool updateTimeVoltagePairs(
      map<string,vector<pair<double,double> >* > const & timeVoltageUpdateMap);
    bool getTimeVoltagePairs(
      map<string,vector<pair<double,double> > > & timeVoltageUpdateMap);

    bool setADCWidths(map<string,int> const & ADCWidthMap);

    // Generic API calls (so far used only by Xygra interface)
    bool getDeviceNames(const string & deviceType,
                        vector<string> & deviceNames);

    int xygraGetNumNodes (const string & deviceName);
    int xygraGetNumWindings (const string & deviceName);
    void xygraGetCoilWindings (const string & deviceName,
                               vector<int> & cW);
    void xygraGetCoilNames (const string & deviceName,
                            vector<string> & cN);
    bool xygraSetConductances (const string & deviceName,
                               const vector<vector<double> > &cM );
    bool xygraSetSources (const string & deviceName,
                          const vector<double> &sV,
                          const double t=0);
    bool xygraSetK (const string & deviceName,
                    const vector<vector<double> > &kM,
                    const double t=0);
    bool xygraGetVoltages (const string & deviceName,
                           vector<double> & vN );

    // two-level newton and pde-continuation
    int  enablePDEContinuation ();
    bool disablePDEContinuation ();
    void getNumInterfaceNodes (vector<int> & numInterfaceNodes);
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
    void setVoltageLimiterFlag ();
    void unsetVoltageLimiterFlag ();
    void addGlobalPar(N_UTL_Param &);
    double getGlobalPar( const string & parName ) const;


    // structs related to options registration
    struct DeviceMgr_OptionsReg : public N_IO_PkgOptionsReg
    {
        DeviceMgr_OptionsReg(DeviceMgr * mgr) : devMgr(mgr) {}

        bool operator()(const N_UTL_OptionBlock & options)
        {return devMgr->registerOptions(options);}

        DeviceMgr * devMgr;
    };

    struct DeviceMgr_SensOptionsReg : public N_IO_PkgOptionsReg
    {
        DeviceMgr_SensOptionsReg(DeviceMgr * mgr) : devMgr(mgr) {}

        bool operator()(const N_UTL_OptionBlock & options)
        {return devMgr->registerSensParams(options);}

        DeviceMgr * devMgr;
    };

    struct DeviceMgr_TimeOptionsReg : public N_IO_PkgOptionsReg
    {
        DeviceMgr_TimeOptionsReg(DeviceMgr * mgr) : devMgr(mgr) {}

        bool operator()(const N_UTL_OptionBlock & options)
        {return devMgr->registerTimeOptions(options);}

        DeviceMgr * devMgr;
    };

    // .TRAN
    struct DeviceMgr_TransAnalysisReg : public N_IO_PkgOptionsReg
    {
        DeviceMgr_TransAnalysisReg(DeviceMgr * mgr) : Mgr(mgr) {}

        bool operator()(const N_UTL_OptionBlock & options)
        {return Mgr->setTranAnalysisParams(options);}

        DeviceMgr * Mgr;
    };

    // .DC
    struct DeviceMgr_DCAnalysisReg : public N_IO_PkgOptionsReg
    {
        DeviceMgr_DCAnalysisReg(DeviceMgr * mgr) : Mgr(mgr) {}

        bool operator()(const N_UTL_OptionBlock & options)
        {return Mgr->setDCAnalysisParams(options);}

        DeviceMgr * Mgr;
    };

    // .OP
    struct DeviceMgr_OPAnalysisReg : public N_IO_PkgOptionsReg
    {
        DeviceMgr_OPAnalysisReg(DeviceMgr * mgr) : Mgr(mgr) {}

        bool operator()(const N_UTL_OptionBlock & options)
        {return Mgr->setOPAnalysisParams(options);}

        DeviceMgr * Mgr;
    };

    // .STEP
    struct DeviceMgr_STEPAnalysisReg : public N_IO_PkgOptionsReg
    {
        DeviceMgr_STEPAnalysisReg(DeviceMgr * mgr) : Mgr(mgr) {}

        bool operator()(const N_UTL_OptionBlock & options)
        {return Mgr->setSTEPAnalysisParams(options);}

        DeviceMgr * Mgr;
    };

    // .MPDE
    struct DeviceMgr_MPDE_AnalysisReg : public N_IO_PkgOptionsReg
    {
        DeviceMgr_MPDE_AnalysisReg(DeviceMgr * devMgr) : devMgr_(devMgr) {}

        bool operator()(const N_UTL_OptionBlock & options)
        {return devMgr_->setMPDEAnalysisParams(options);}

        DeviceMgr * devMgr_;
    };

    // .HB
    struct DeviceMgr_HB_AnalysisReg : public N_IO_PkgOptionsReg
    {
        DeviceMgr_HB_AnalysisReg(DeviceMgr * devMgr) : devMgr_(devMgr) {}

        bool operator()(const N_UTL_OptionBlock & options)
        {return devMgr_->setHBAnalysisParams(options);}

        DeviceMgr * devMgr_;
    };

    // .AC
    struct DeviceMgr_AC_AnalysisReg : public N_IO_PkgOptionsReg
    {
        DeviceMgr_AC_AnalysisReg(DeviceMgr * devMgr) : devMgr_(devMgr) {}

        bool operator()(const N_UTL_OptionBlock & options)
        {return devMgr_->setACAnalysisParams(options);}

        DeviceMgr * devMgr_;
    };

    // .MOR
    struct DeviceMgr_MOR_AnalysisReg : public N_IO_PkgOptionsReg
    {
        DeviceMgr_MOR_AnalysisReg(DeviceMgr * devMgr) : devMgr_(devMgr) {}

        bool operator()(const N_UTL_OptionBlock & options)
        {return devMgr_->setMORAnalysisParams(options);}

        DeviceMgr * devMgr_;
    };

    // functions related to options registration
    bool registerOptions(const N_UTL_OptionBlock & OB) { return devOptions_.registerOptions(OB); }
    bool registerSensParams(const N_UTL_OptionBlock & OB);
    bool registerTimeOptions(const N_UTL_OptionBlock & OB){return true;}
    bool setTranAnalysisParams(const N_UTL_OptionBlock & OB);
    bool setDCAnalysisParams(const N_UTL_OptionBlock & OB);
    bool setOPAnalysisParams(const N_UTL_OptionBlock & OB);
    bool setSTEPAnalysisParams(const N_UTL_OptionBlock & OB);
    bool setMPDEAnalysisParams(const N_UTL_OptionBlock & OB);
    bool setHBAnalysisParams(const N_UTL_OptionBlock & OB);
    bool setACAnalysisParams(const N_UTL_OptionBlock & OB);
    bool setMORAnalysisParams(const N_UTL_OptionBlock & OB);

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

    int restartDataSize(bool pack);

    // Output restart data.
    bool dumpRestartData(char * buf, int bsize, int & pos, N_PDS_Comm * comm, bool pack );

    // Load restart data.
    bool restoreRestartData(char * buf, int bsize, int & pos, N_PDS_Comm * comm, bool pack );

  protected:

  private:
    bool createDeviceByModelType(const int model_type);
    bool createDeviceByModelLevel(const std::string &model_name, const int level);
    int getDeviceTypeOffset(const ModelBlock & MB);
    int getModelTypeIndex(const std::string & ModelType);
    int getDeviceIndex(const std::string & Name);

    //bool loadPerturbationVector_ ();
    bool setupSolverInfo_ ();
    bool setupRawVectorPointers_ ();
    bool setupRawMatrixPointers_ ();

    // copy constructor should be private: This class should be
    // a singleton.
    DeviceMgr(const DeviceMgr &right);
    DeviceMgr(N_IO_CmdParse & cp);

    void setUpDeviceIndexMap_();
    void setUpMESFETLevelMap_();
    void setUpMOSFETLevelMap_();
    void setUpJFETLevelMap_();
    void setUpDiodeLevelMap_();
    void setUpBJTLevelMap_();
    void setUpPDELevelMap_();
    void setUpResistorLevelMap_();
    void setUpIndLevelMap_();
    void setUpNeuronLevelMap_();
    void setUpNeuronPopLevelMap_();
    void setUpSynapseLevelMap_();
    void setUpRadLevelMap_();
    void setUpDeviceModelTypeMap_();

    bool setUpPDEDeviceFlagArray_ ();

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
    N_DEV_ExternDeviceInstance * addExtDeviceInstance_(InstanceBlock & IB);
#endif

    Xygra::Instance * getXygraInstancePtr_(const std::string & deviceName);

    // attributes:
  public:
    Device * dummy_ptr;

  protected:

  private:
    Device * deviceArray_[ModelType::NUMDEV];
    int deviceAllocFlag_[ModelType::NUMDEV];
    int deviceUseFlag_[ModelType::NUMDEV];
    int PDEDeviceFlag_[ModelType::NUMDEV];

    int numSensParams_;

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

    DeviceBuilder devBuilder_;

    ExternData extData_;

    N_LAS_Vector       *  numJacSolVectorPtr_;
    N_LAS_Vector       *  numJacStaVectorPtr_;
    N_LAS_Vector       *  numJacStoVectorPtr_;

    N_LAS_Vector       *  diagonalVectorPtr_;

    N_LAS_System        * lasSysPtr_;

    N_ANP_AnalysisInterface       * anaIntPtr_;

    N_IO_OutputMgr      * outputMgrPtr_;

    N_PDS_Manager       * pdsMgrPtr_;

    N_NLS_Manager       * nlsMgrPtr_;

    DeviceMgr     * DevMgrPtr_;

    RCP<N_IO_PkgOptionsMgr> pkgOptMgrPtr_;

    vector< pair<int,double> > *        icLoads_;
    DeviceIndexMap                      deviceIndexMap_;
    map<int, map<int,int> *>            modelLevelMap_;
    DeviceModelTypeMap                  deviceModelTypeMap_;
    map<string,int>                     deviceModelNameMap_;
    map<string,int>                     deviceModelNameBaseTypeMap_;
    map<string,int>                     localDeviceCountMap_;
    map<int,int>                        MESFETLevelMap_;
    map<int,int>                        MOSFETLevelMap_;
    map<int,int>                        JFETLevelMap_;
    map<int,int>                        DiodeLevelMap_;
    map<int,int>                        BJTLevelMap_;
    map<int,int>                        PDELevelMap_;
    map<int,int>                        radLevelMap_;
    map<int,int>                        ResistorLevelMap_;
    map<int,int>                        INDLevelMap_;
    map<int,int>                        NeuronLevelMap_;
    map<int,int>                        NeuronPopLevelMap_;
    map<int,int>                        SynapseLevelMap_;

    map<string, Xygra::Instance *> xygraPtrMap_;

    multimap <int,DeviceInstance*> solDevInstMap_;

    DeviceVector devicePtrVec_;
    DeviceVector pdeDevicePtrVec_;
    DeviceVector nonPdeDevicePtrVec_;

    vector <DeviceInstance*> instancePtrVec_;
    vector <DeviceInstance*> bpInstancePtrVec_; // instances with breakpoints functions
    vector <DeviceInstance*> pdeInstancePtrVec_;
    vector <DeviceInstance*> nonPdeInstancePtrVec_;
    vector <DeviceInstance*> mosfetInstancePtrVec_;
    vector <DeviceInstance*> vsrcInstancePtrVec_;
    vector <DeviceInstance*> bjtInstancePtrVec_;
    map <string, Vsrc::Instance*> vsrcInstancePtrMap_;

    vector <DeviceInstance*> plotFileInstancePtrVec_;

    map <string,SourceInstance*> indepSourcePtrMap_;
    vector <SourceInstance*> indepSourceInstancePtrVec_;

    // this is used to store the contents of the indepSourceInstancePtrVec_
    // during an mpde initialization where we'll remove slow sources from
    // the that vector so that they don't get updated
    vector <SourceInstance*> indepSourceInstanceBackupPtrVec_;

    set<string> devicesNeedingLeadCurrentLoads_;

#ifdef Xyce_EXTDEV
    vector <N_DEV_ExternDeviceInstance*> extDevInstancePtrVec_;
    vector<InstanceBlock*> extDevIBPtrVec_;

    map <string, int> passThroughParamsMap_;
#endif

    // vector of pointers to devices under test.
    vector <DeviceInstance*> testJacDevicePtrVec_;

    vector <DeviceModel*> modelPtrVec_;
    vector <DeviceModel*> mosfetModelPtrVec_;
    vector <DeviceModel*> bsim3ModelPtrVec_;
    vector <DeviceModel*> bsim4ModelPtrVec_;
    vector <DeviceModel*> bsimsoiModelPtrVec_;
    vector <DeviceModel*> bjtModelPtrVec_;
    vector <DeviceModel*> diodeModelPtrVec_;

    vector<DeviceEntity*> dependentPtrVec_;

    vector<int> numInterfaceNodes_;
    int numPDEDevices_;
    bool calledBeforeCSPI;

    int numThreads_;
    bool multiThreading_;

    // command line reference:
    N_IO_CmdParse & commandLine_;

    // real time solver data:
    SolverState solState_;
    SolverState solStateExternal_;
    bool externalStateFlag_;

    // user-defined options:
    DeviceOptions devOptions_;

    // sensitivities:
    N_DEV_DeviceSensitivities * devSensPtr_;

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

    bool setupSolverInfo() {
      return setupSolverInfo_();
    }
};

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::returnDevicePtr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
inline Device * DeviceMgr::returnDevicePtr (int Index)
{
  if (Index >= ModelType::NUMDEV) Index = ModelType::NUMDEV - 1;
  Device * devptr = deviceArray_[Index];
  return devptr;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::getNumSupportedDevices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/16/05
//-----------------------------------------------------------------------------
inline int DeviceMgr::getNumSupportedDevices ()
{
  return ModelType::NUMDEV;
}


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
inline bool DeviceMgr::registerOutputMgr (N_IO_OutputMgr * tmp_outputMgrPtr)
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
inline void DeviceMgr::setVoltageLimiterFlag ()
{
  devOptions_.voltageLimiterFlag = true;
}

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
//// Function      : DeviceMgr::getDeviceCountMap()
//// Purpose       :
//// Special Notes :
//// Scope         : public
//// Creator       : Eric R. Keiter, SNL
//// Creation Date : 05/07/2010
////-----------------------------------------------------------------------------
inline const map<string,int> & DeviceMgr::getDeviceCountMap()
{
  return localDeviceCountMap_;
}

//-----------------------------------------------------------------------------
// Function      : DeviceMgr::registerICLoads
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Computational Sciences
// Creation Date : 04/03/01
//-----------------------------------------------------------------------------
inline bool DeviceMgr::registerICLoads( vector< pair<int,double> > * icLoads )
{
  return (icLoads_ = icLoads) != 0;
}

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::DeviceMgr N_DEV_DeviceMgr;

#endif


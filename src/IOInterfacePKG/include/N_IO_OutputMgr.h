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
// Filename       : $RCSfile: N_IO_OutputMgr.h,v $
//
// Purpose        : Generate global id structures and proc maps
//                  and distribute nodes to processors
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 10/10/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.188.2.2 $
//
// Revision Date  : $Date: 2014/03/13 21:17:52 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_OutputMgr_h
#define Xyce_N_IO_OutputMgr_h

#include <iterator>
#include <list>
#include <set>
#include <string>
#include <time.h>
#include <vector>

#ifdef Xyce_USE_HDF5
#include <hdf5.h>
#endif

#include <N_ANP_AnalysisInterface.h>
#include <N_ANP_SweepParam.h>
#include <N_ANP_fwd.h>
#include <N_DEV_fwd.h>
#include <N_IO_FourierMgr.h>
#include <N_IO_MeasureManager.h>
#include <N_IO_Objective.h>
#include <N_IO_Outputter.h>
#include <N_IO_PkgOptionsMgr.h>
#include <N_IO_fwd.h>
#include <N_PDS_Comm.h>
#include <N_PDS_MPI.h>
#include <N_PDS_Serial.h>
#include <N_UTL_Misc.h>
#include <N_UTL_NoCase.h>
#include <N_UTL_OptionBlock.h>
#include <N_UTL_Xyce.h>
#include <N_UTL_fwd.h>

class N_LAS_BlockVector;
class N_LAS_Matrix;
class N_LAS_Vector;
class N_MPDE_Manager;

namespace Xyce {
namespace IO {

typedef std::map<std::string, Objective> ObjectiveMap;
typedef std::map<PrintType::PrintType, Outputter::Interface *> OutputterMap;
typedef std::vector<std::vector<Outputter::Interface *> > ActiveOutputterStack;
typedef std::map<OutputType::OutputType, PrintParameters> OutputParameterMap;

enum Domain {DOMAIN_TIME, DOMAIN_FREQUENCY};

//-----------------------------------------------------------------------------
// Function      : getWidthFromStaticIndex
//
// Purpose       : This function figures out the correct column width from the
//                 staticIndex_ and the delimiter_
//
// Special Notes :
//
// Creator       : Todd Coffey, 1414
// Creation Date : 9/19/08
//-----------------------------------------------------------------------------
int getWidthFromStaticIndex(int index, const std::string &delim);

void printLineDiagnostics();

//-----------------------------------------------------------------------------
// Class         : OutputMgr
// Purpose       :
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/10/00
//-----------------------------------------------------------------------------
class OutputMgr
{
public:
  class ActiveOutput
  {
  public:
    ActiveOutput(OutputMgr &output_manager)
      : outputManager_(output_manager)
    {
      outputManager_.pushActiveOutputters();
    }

    ~ActiveOutput()
    {
      outputManager_.finishOutput();
      outputManager_.popActiveOutputters();
    }

    void add(PrintType::PrintType print_type)
    {
      outputManager_.addActiveOutputter(print_type);
    }

  private:
    OutputMgr &     outputManager_;
  };

  typedef std::map<std::string, std::pair<int, std::ostream *> > OpenPathStreamMap;

  // Factory to generate instance of class
  static OutputMgr * factory(CmdParse & cp);

  // Destructor
  ~OutputMgr();

  // registration functions:

  bool registerTopology(Topo::Topology * topology_interface)
  {
    topology_ = topology_interface;
    return true;
  }

  bool registerDeviceInterface(Device::DeviceInterface * device_interface)
  {
    deviceInterface_ = device_interface;
    return true;
  }

  bool registerAnalysisInterface(N_ANP_AnalysisInterface * analysis_interface)
  {
    analysisInterface_ = analysis_interface;
    return true;
  }

  bool registerMPDEManager(N_MPDE_Manager * mpde_manager)
  {
    mpdeMgrPtr_ = mpde_manager;
    return true;
  }

  bool getOutputIntervals(double & initialInterval, std::vector<std::pair< double, double > > & intervalPairs);

  bool registerPkgOptionsMgr( PkgOptionsMgr &pkgOpt);
  bool registerParallelServices(N_PDS_Comm * tmp_pds_ptr);
  bool registerDCOPOptions(const Util::OptionBlock & option_block);
  bool registerDCOptions(const Util::OptionBlock & option_block);
  bool registerTranOptions(const Util::OptionBlock & option_block);
  bool registerMPDETranOptions(const Util::OptionBlock & option_block);
  bool registerHBOptions(const Util::OptionBlock & option_block);
  bool registerSTEPOptions(const Util::OptionBlock & option_block);
  bool registerOutputOptions(const Util::OptionBlock & option_block);
  bool registerDeviceOptions(const Util::OptionBlock & option_block);
  bool registerDistribNodes(const Util::OptionBlock & option_block);

  bool registerPRINTSet(const Util::OptionBlock & option_block);
  bool registerIC(const Util::OptionBlock & option_block);
  bool registerNodeSet(const Util::OptionBlock & option_block);
  bool registerSave(const Util::OptionBlock & option_block);
  bool registerLoad(const Util::OptionBlock & option_block);
  bool registerOP (const Util::OptionBlock & option_block);
  bool registerSens(const Util::OptionBlock & option_block);

  bool setRESULTParams(const Util::OptionBlock & option_block);
  bool setOBJECTIVEParams(const Util::OptionBlock & option_block);

  // Turns on output
  void enableOverrideRawOutput(const PrintParameters &print_parameters);
  void enableTransientOutput();
  void enableMPDEOutput();
  void enableACOutput();
  void enableDCSweepOutput();
  void enableHBOutput();
  void enableHomotopyOutput(Analysis::Analysis_Mode analysis_mode);
  void enableSensitivityOutput (Analysis::Analysis_Mode analysis_mode);
  void enableMPDEOutput(const CmdParse & commandLine, PrintParameters &print_parameters);

  void prepareOutput(
     Analysis::Analysis_Mode                         analysis_mode,
     const std::vector<N_ANP_SweepParam> &     step_sweep_parameters,
     const std::vector<N_ANP_SweepParam> &     dc_sweep_parameterds);

  void setSweepParameters(const std::vector<N_ANP_SweepParam> &step_sweep_parameters, const std::vector<N_ANP_SweepParam> &dc_sweep_parameters);

  void fixupNodeNames();
  void fixupPrintParameters(PrintParameters &print_parameters);

  // Runs specified output commands
  void output(
     const double & time,
     const int stepNumber,
     const int maxStep,
     const std::vector<N_ANP_SweepParam> & stepParamVec1,
     const int dcNumber,
     const int maxDC,
     const std::vector<N_ANP_SweepParam> & dcParamVec1,
     N_LAS_Vector * solnVecPtr,
     N_LAS_Vector * stateVecPtr,
     N_LAS_Vector * storeVecPtr,
     const std::vector<double> & objectiveVec,
     const std::vector<double> & dOdpVec, 
     const std::vector<double> & dOdpAdjVec,
     const std::vector<double> & scaled_dOdpVec, 
     const std::vector<double> & scaled_dOdpAdjVec,
     bool skipPrintLineOutput=false);

  void registerNodeDevNames(const std::set<std::string> * nodeNamesIn,
                            const std::map<std::string, RefCountPtr<Device::InstanceBlock> > * devNamesIn);
  void printLineDiagnostics();
  void delayedPrintLineDiagnostics ();

  bool setupInitialConditions ( N_LAS_Vector & solnVec, N_LAS_Vector & flagVec);
  void outputDCOP(N_LAS_Vector & solnVec);

  const Device::DeviceInterface &getDeviceInterface() const
  {
    return *deviceInterface_;
  }

  Device::DeviceInterface &getDeviceInterface()
  {
    return *deviceInterface_;
  }

  NodeNamePairMap &getICData( int &, std::string &);

  NodeNamePairMap &getAllNodes()
  {
    return allNodes_;
  }

  const NodeNamePairMap &getAllNodes() const
  {
    return allNodes_;
  }

  const NodeNamePairMap &getStateNodes() const
  {
    return stateNodes_;
  }

  const NodeNamePairMap &getStoreNodes() const
  {
    return storeNodes_;
  }

  // Runs specified RESULT output commands
  void outputRESULT (N_LAS_Vector * solnVecPtr, N_LAS_Vector * stateVecPtr, N_LAS_Vector * storeVecPtr);

  // // Finishes STEP/RESULT output
  void finishOutputSTEP ();

  void outputAC(double freq, const N_LAS_Vector * freqDomainSolnVecReal, const N_LAS_Vector * freqDomainSolnVecImaginary);
  void outputMPDE (double time, const N_LAS_Vector * solnVecPtr );
  void outputHomotopy(const std::vector<std::string> & parameter_names, const std::vector<double> & param_values, const N_LAS_Vector * solution_vector);
  void outputHB(
     const int                             stepNumber,
     const int                             maxStep,
     const std::vector<N_ANP_SweepParam> & stepParamVec1,
     const std::vector<double>& timePoints,
     const std::vector<double>& freqPoints,
     const N_LAS_BlockVector & timeDomainSolnVec,
     const N_LAS_BlockVector & freqDomainSolnVecReal,
     const N_LAS_BlockVector & freqDomainSolnVecImaginary);
  void outputMORTF(bool origSys, const double & freq, const Teuchos::SerialDenseMatrix<int, std::complex<double> >& H);

  void outputROM(
     const Teuchos::SerialDenseMatrix<int, double>& Ghat,
     const Teuchos::SerialDenseMatrix<int, double>& Chat,
     const Teuchos::SerialDenseMatrix<int, double>& Bhat,
     const Teuchos::SerialDenseMatrix<int, double>& Lhat);

  void outputROM(
     const N_LAS_Matrix& Ghat,
     const N_LAS_Matrix& Chat,
     const Teuchos::SerialDenseMatrix<int, double>& Bhat,
     const Teuchos::SerialDenseMatrix<int, double>& Lhat);

  void recordDakotaFileName (std::string &);
  void recordDakotaParNames (std::vector<std::string> &);
  void outputDakotaResults();

  void resetOutput();
  void finishOutput();

  // if any macro-simulation level functions were defined, finish them and output results
  void outputMacroResults();

  // routines to get/set variables that external programs such as Dakota want to query the value of.
  void finalizeResponseVars();
  bool registerResponseVars (const std::string & objString, RCP< std::vector< double > > varVectorPtr );
  void saveResponseVarValues( N_LAS_Vector * solnVecPtr );

  void getMeasureValue (const std::string &name, double &result, bool &found) const;

  // Additional Public Declarations

  std::vector<std::string> getVariableNames();

  void setExternalNetlistParams( std::vector<std::pair< std::string, std::string> > & externalNetlistParams );

  void remeasure();

  // set the suffix that will be used to generate the outputfile
  // baseFileName + suffix + ".prn"
  // This lets MPDE save a start-up file and an initial condition file
  // as well as the expected baseFileName.prn file
  void setOutputFilenameSuffix( std::string newSuffix )
  {
    filenameSuffix_ = newSuffix;
  };

  void setResponseFilename( std::string aFileName )
  {
    responseFileName_ = aFileName;
    responseFileNameGiven_ = true;
  };

  N_PDS_Comm * getCommPtr() const
  {
    return pdsCommPtr_;
  }

  int getProcID() const
  {
    return pdsCommPtr_->procID();
  }

  int getNumProcs() const
  {
    return pdsCommPtr_->numProc();
  }

  const std::string &getTitle() const
  {
    return title_;
  }

  void setTitle(const std::string &title)
  {
    title_ = title;
  }

  const std::string &getNetListFilename() const
  {
    return netListFilename_;
  }

  const CmdParse &getCommandLine() const
  {
    return commandLine_;
  }

  const std::vector<N_ANP_SweepParam> &getStepParamVec() const
  {
    return stepParamVec_;
  }

  int getStepLoopNumber() const
  {
    return stepLoopNumber_;
  }

  int getMaxParamSteps() const
  {
    return maxParamSteps_;
  }

  const std::vector<N_ANP_SweepParam> &getDCParamVec() const
  {
    return dcParamVec_;
  }

  int getDCLoopNumber() const
  {
    return dcLoopNumber_;
  }

  int getMaxDCSteps() const
  {
    return maxDCSteps_;
  }

  PrintType::PrintType getPrintType() const
  {
    return PRINTType_;
  }

  bool getTempSweepFlag() const
  {
    return tempSweepFlag_;
  }

  const ParameterList *getVariableList() const
  {
    return printParameters_ ? &printParameters_->variableList_ : 0;
  }

  double getPRINTDCvalue() const
  {
    return PRINTdcvalue_;
  }

  double getPRINTDCstart() const
  {
    return PRINTdcstart_;
  }

  double getPRINTDCstop() const
  {
    return PRINTdcstop_;
  }

  const N_ANP_AnalysisInterface *getAnaIntPtr() const
  {
    return analysisInterface_;
  }

  const std::string &getPRINTDCname() const
  {
    return PRINTdcname_;
  }

  bool getPrintEndOfSimulationLine() const
  {
    return printEndOfSimulationLine_;
  }

  bool getOutputVersionInRawFile() const
  {
    return outputVersionInRawFile_; 
  }


  void setEnableHomotopyFlag(bool value)
  {
    enableHomotopyFlag_ = value;
  }

  const Topo::Topology *getTopPtr() const
  {
    return topology_;
  }

  double getCircuitTime() const
  {
    return circuitTime_;
  }

  void setCircuitTime(double time)
  {
    circuitTime_ = time;
  }

  double getCircuitFrequency() const
  {
    return circuitFrequency_;
  }

  void setCircuitFrequency(double frequency)
  {
    circuitFrequency_ = frequency;
  }

  double getCircuitTemp() const
  {
    return circuitTemp_;
  }

  bool getSTEPEnabledFlag() const
  {
    return STEPEnabledFlag_;
  }

  const std::string &getFilenameSuffix() const
  {
    return filenameSuffix_;
  }

  Format::Format getFormat() const
  {
    return format_;
  }

  bool getNoIndex() const
  {
    return noIndex_;
  }

  bool getPRINTcsv() const
  {
    return format_ == Format::CSV;
  }

  bool getPRINTprobe() const
  {
    return format_ == Format::PROBE;
  }

  bool getPRINTtecplot() const
  {
    return format_ == Format::TECPLOT;
  }

  bool getDetailedDeviceFlag() const
  {
    return detailedDeviceCountFlag_;
  }

  const Measure::Manager &getMeasureManager() const
  {
    return measureManager_;
  }

  Measure::Manager &getMeasureManager()
  {
    return measureManager_;
  }

  const ObjectiveMap &getObjectiveMap() const
  {
    return objective_;
  }

  ObjectiveMap &getObjectiveMap()
  {
    return objective_;
  }

  const N_MPDE_Manager *getMpdeMgrPtr() const
  {
    return mpdeMgrPtr_;
  }

  int getCurrentOutputterIndex() const
  {
    return currentOutputter_->getIndex();
  }

  double getTemperature() const
  {
    return circuitTemp_;
  }

  double getTime() const
  {
    return circuitTime_;
  }

  double getFrequency() const
  {
    return circuitFrequency_;
  }

  double getStepSweep(size_t index) const 
  {
    return stepParamVec_[index].currentVal;
  }

  double getDCSweep(size_t index) const 
  {
    return dcParamVec_[index].currentVal;
  }

  void pushActiveOutputters() 
  {
    activeOutputterStack_.push_back(std::vector<Outputter::Interface *>());
  }

  void popActiveOutputters() 
  {
    activeOutputterStack_.pop_back();
  }

  void addActiveOutputter(PrintType::PrintType print_type) 
  {
    OutputterMap::iterator it = outputterMap_.find(print_type);

    if (!activeOutputterStack_.empty() && it != outputterMap_.end())
      activeOutputterStack_.back().push_back((*it).second);
  }

  void setCurrentOutputter(Outputter::Interface *outputter) 
  {
    currentOutputter_ = outputter;
  }

  Objective &getObjective(const std::string &varType) const 
  {
    ObjectiveMap::const_iterator it = objective_.find(varType);
    return const_cast<Objective &>((*it).second);
  }

  std::ostream *openFile(const std::string &path, std::ios_base::openmode mode);
  std::ostream *openFile(const std::string &path);
  std::ostream *openBinaryFile(const std::string &path);
  int closeFile(std::ostream *os);

private:
  // Copy constructor should be private: This class should be a singleton.
  OutputMgr( const OutputMgr & rhs );
  OutputMgr & operator=( const OutputMgr & rhs );

  // Default constructor, accessible only through singleton factory (private)
  OutputMgr( CmdParse & cp);

  void outputDCOP_restartFile ( N_LAS_Vector & solnVec);
  bool inputDCOP_( N_LAS_Vector & solnVec, N_LAS_Vector & flagVec);

  bool setupIC_or_NODESET( N_LAS_Vector & solnVec, N_LAS_Vector & flagVec,
                           bool & useFlag,
                           std::string & icType,
                           std::vector<Util::OptionBlock> & initBlockVec);

  void outputIC_or_NODESET (N_LAS_Vector & solnVec);

  std::set<std::string> deviceParamCheck;

  bool enableHDF5Output();
  bool updateHDF5Output( N_LAS_Vector * solnVecPtr);
  bool closeHDF5Output();

private:
  OutputterMap                outputterMap_;
  Outputter::Interface *      currentOutputter_;
  OutputParameterMap          outputParameterMap_;
  ActiveOutputterStack        activeOutputterStack_;

  // Output file name
  // a suffix that will be added to the output filename
  // to distinguish different types of output such as multiple iteration numbers
  // circuit.cir.1.prn, circuit.cir.2.prn or domains like circuit.cir.FD.prn
  // each derived type from OutputFileBase handles a particular
  // format of output like standard, probe, tecplot etc.
  // Store a vector of them so that more than one type can be output at the
  // same time
  std::string                 netListFilename_;
  std::string                 title_;
  std::string                 filenameSuffix_;

  std::vector< OutputFileBase *> outputHandlersVec_;

  std::string responseFileName_;
  bool responseFileNameGiven_;

  // Pointer to the ostream for .result output
  std::ostream *              resultStreamPtr_;

  // std::map<std::string,int>   localDeviceCountMap_;

  bool output_op_;
  bool input_op_;
  std::string output_op_file_;
  std::string input_op_file_;

  std::string saveFileLevel_;
  std::string saveFileType_;
  std::string saveOutputFile_;

  std::vector<std::string> stepParams_;
  std::vector<std::string> dcParams_;

  // For early check of .print argument validity
  const std::set<std::string> *nodeNames_;
  const std::map<std::string, RefCountPtr<Device::InstanceBlock> > *devNames_;

  // these variables used to be static:
  unsigned long staticIndex_;
  unsigned long staticMPDEIndex_;

  bool outputPrnAdded_;

  N_PDS_Comm *                pdsCommPtr_;
  Topo::Topology *            topology_;
  Device::DeviceInterface *     deviceInterface_;
  N_ANP_AnalysisInterface *   analysisInterface_;
  N_MPDE_Manager *            mpdeMgrPtr_;

  // print statement vars
  bool rawFlag_;
  bool tranFlag_;
  bool probeFlag_;
  bool DCSweepFlag_;
  bool ACFlag_;
  bool enableHomotopyFlag_;
  bool enableSensitivityFlag_;
  bool homotopyFlag_;
  bool sensitivityFlag_;
  bool HBFlag_;
  bool MPDEFlag_;
  bool PRINTflag_;

  PrintType::PrintType PRINTType_;

  Format::Format      format_;
  bool                noIndex_;

  PrintParameters *   printParameters_;
  PrintParameters     defaultPrintParameters_;

  Format::Format formatResult_;

  double PRINTdcstart_;
  double PRINTdcstop_;
  double PRINTdcvalue_;
  std::string PRINTdcname_;

  // initial condition vars
  bool ICflag_;
  bool saveFlag_;
  bool loadFlag_;
  bool outputOnceAlreadyFlag_;
  std::vector<Util::OptionBlock> ICblockVec_;
  bool nodesetflag_;
  std::vector<Util::OptionBlock> nodesetblockVec_;

  double initialOutputInterval_;
  std::vector<std::pair< double, double > > outputIntervalPairs_;

  std::set< std::string> distribVsrcSet_;
  std::set< std::string> distribVnodeSet_;

  bool STEPEnabledFlag_;
  int  STEPcounter_;
  bool RESULTinitialized_;

  bool printEndOfSimulationLine_;  // flag to indicate if user wants the "End of Xyce(TM)" line in the output.
  bool outputVersionInRawFile_;    // flag to indicate that Version should be output in the header of a RAW file.

  // response functions requested in the external params passed to Xyce.
  std::vector< std::pair< std::string, std::string> > variablesUsedInSimulation_ ;
  std::vector< std::pair< std::string, std::string> > responseFunctionsRequested_ ;
  std::vector< std::pair< std::string, std::string> > derivativeVariablesRequested_ ;
  std::vector< std::pair< std::string, std::string> > analysisComponentsRequested_ ;

  // response related variable lists
  Util::OpList responseVarList_;
  RCP< std::vector< double > > responseVarPtr_;
  std::vector<std::string> responseNames_;
  int numResponseVars_;

  bool tempSweepFlag_;
  bool outputCalledBefore_;

  // vector of expressions:
  std::vector<Util::ExpressionData*> expVector_;

  // vector of expressions from .RESULT
  std::vector<Util::ExpressionData*> resultVector_;

  ObjectiveMap objective_;

  // step loop information
  int stepLoopNumber_;
  int maxParamSteps_;
  std::vector<N_ANP_SweepParam> stepParamVec_;

  // dc loop information
  int dcLoopNumber_;
  int maxDCSteps_;
  std::vector<N_ANP_SweepParam> dcParamVec_;

  double circuitTime_;                        ///< transient circuit time:
  double circuitFrequency_;                   ///< ac current frequency
  double circuitTemp_;                        ///< circuit temperature

  // command line parser
  CmdParse & commandLine_;

  bool isStarredPrintLineProcessed;

  NodeNamePairMap allNodes_;
  NodeNamePairMap stateNodes_;
  NodeNamePairMap storeNodes_;
  NodeNamePairMap externNodes_;

  bool detailedDeviceCountFlag_;

  NodeNamePairMap opData_;
  int op_found_, total_soln_;

  int icType_;

  // HDF5 file support vars
  bool hdf5FileNameGiven_;
  bool hdf5HeaderWritten_;
  std::string hdf5FileName_;
  int hdf5IndexValue_;

#ifdef Xyce_USE_HDF5
  hid_t hdf5FileId_;
  hid_t hdf5PlistId_;
#endif  // Xyce_USE_HDF5

  // .sens output variables
  bool sensObjFuncGiven_;
  std::string sensObjFunction_;
  std::vector<std::string> sensParamNameVec_;
  std::vector<std::string> sensParFullNames_;



  Measure::Manager            measureManager_;
  FourierMgr                  fourierManager_;

  OpenPathStreamMap           openPathStreamMap_;
};

std::ostream &printGlobalDeviceCounts(std::ostream &os, Parallel::Machine comm, 
                                      const std::map<std::string,int> &local_device_count_map);
std::ostream &printLocalDeviceCount(std::ostream &os, Parallel::Machine comm, 
                                    const std::map<std::string, int> & local_device_count_map);

} // namespace IO
} // namespace Xyce

typedef Xyce::IO::OutputMgr N_IO_OutputMgr;

#endif

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
// Revision Number: $Revision: 1.165.2.7 $
//
// Revision Date  : $Date: 2013/10/03 17:23:42 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_OutputMgr_h
#define Xyce_N_IO_OutputMgr_h

// ---------- Standard Includes ----------
#include <list>
#include <string>
#include <vector>
#include <set>
#include <time.h>

#ifdef Xyce_USE_HDF5
#include <hdf5.h>
#endif

// ----------   Xyce Includes   ----------
#include <N_DEV_fwd.h>
#include <N_UTL_Xyce.h>
#include <N_UTL_Misc.h>
#include <N_UTL_OptionBlock.h>
#include <N_UTL_NoCase.h>

#include <N_IO_Outputter.h>
#include <N_IO_PkgOptionsMgr.h>
#include <N_ANP_SweepParam.h>
#include <N_ANP_AnalysisInterface.h>
#include <N_IO_Objective.h>
#include <N_IO_MeasureManager.h>
#include <N_IO_FourierMgr.h>

#include <N_PDS_Comm.h>

// ---------- Trilinos Includes ----------
#include <Teuchos_SerialDenseMatrix.hpp>

#ifdef Xyce_Dakota
#include <Epetra_SerialDenseVector.h>
#endif

// ---------- Forward Declarations ----------

class N_ANP_AnalysisInterface;
class N_IO_CmdParse;
class N_IO_MeasureBase;
class N_IO_OutputFileBase;
class N_LAS_BlockVector;
class N_LAS_Matrix;
class N_LAS_Vector;
class N_MPDE_Manager;
class N_TOP_Topology;
class N_UTL_Expression;
class N_UTL_ExpressionData;

namespace Xyce {
namespace IO {

typedef std::list<N_UTL_Param> ParameterList;
typedef std::map<std::string, N_IO_Objective> ObjectiveMap;
typedef std::map<PrintType::PrintType, Outputter::Interface *> OutputterMap;
typedef std::vector<std::vector<Outputter::Interface *> > ActiveOutputterStack;
typedef std::map<OutputType::OutputType, PrintParameters> OutputParameterMap;

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
// Class         : N_IO_OutputMgr
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

    enum Domain {DOMAIN_TIME, DOMAIN_FREQUENCY};

    // Factory to generate instance of class
    static OutputMgr * factory(N_IO_CmdParse & cp);

    // Destructor
    ~OutputMgr();

    // registration functions:

    bool registerTopology(N_TOP_Topology * topology_interface) {
      topPtr_ = topology_interface;
      return true;
    }

    bool registerDeviceInterface(N_DEV_DeviceInterface * device_interface) {
      devPtr_ = device_interface;
      return true;
    }

    bool registerAnalysisInterface(N_ANP_AnalysisInterface * analysis_interface) {
      anaIntPtr_ = analysis_interface;
      return true;
    }

    bool registerMPDEManager(N_MPDE_Manager * mpde_manager) {
      mpdeMgrPtr_ = mpde_manager;
      return true;
    }

    bool getOutputIntervals(double & initialInterval, vector < pair < double, double > > & intervalPairs);

    bool registerPkgOptionsMgr( N_IO_PkgOptionsMgr &pkgOpt);
    bool registerParallelServices(N_PDS_Comm * tmp_pds_ptr);
    bool registerDCOPOptions(const N_UTL_OptionBlock & option_block);
    bool registerDCOptions(const N_UTL_OptionBlock & option_block);
    bool registerTranOptions(const N_UTL_OptionBlock & option_block);
    bool registerMPDETranOptions(const N_UTL_OptionBlock & option_block);
    bool registerHBOptions(const N_UTL_OptionBlock & option_block);
    bool registerSTEPOptions(const N_UTL_OptionBlock & option_block);
    bool registerOutputOptions(const N_UTL_OptionBlock & option_block);
    bool registerDeviceOptions(const N_UTL_OptionBlock & option_block);
    bool registerDistribNodes(const N_UTL_OptionBlock & option_block);

    bool registerPRINTSet(const N_UTL_OptionBlock & option_block);
    bool registerIC(const N_UTL_OptionBlock & option_block);
    bool registerNodeSet(const N_UTL_OptionBlock & option_block);
    bool registerSave(const N_UTL_OptionBlock & option_block);
    bool registerLoad(const N_UTL_OptionBlock & option_block);
    bool registerOP (const N_UTL_OptionBlock & option_block);

    bool setRESULTParams(const N_UTL_OptionBlock & option_block);
    bool setOBJECTIVEParams(const N_UTL_OptionBlock & option_block);

    // Turns on output
    void enableOverrideRawOutput(const PrintParameters &print_parameters);
    void enableTransientOutput();
    void enableMPDEOutput();
    void enableACOutput();
    void enableDCSweepOutput();
    void enableHBOutput();
    void enableHomotopyOutput(ANP_Analysis_Mode analysis_mode);
    void enableMPDEOutput(const N_IO_CmdParse & commandLine, PrintParameters &print_parameters);

    void prepareOutput(
      ANP_Analysis_Mode                         analysis_mode,
      const std::vector<N_ANP_SweepParam> &     step_sweep_parameters,
      const std::vector<N_ANP_SweepParam> &     dc_sweep_parameterds);

    void setSweepParameters(const std::vector<N_ANP_SweepParam> &step_sweep_parameters, const std::vector<N_ANP_SweepParam> &dc_sweep_parameters);

    void fixupPrintParameters(Domain domain, PrintParameters &print_parameters, bool expandComplexTypes=true);

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
      bool skipPrintLineOutput=false);

    void printDeviceCounts ();
    void addDeviceToCount (std::string & devNameIn);
    void registerNodeDevNames(const std::set<std::string> * nodeNamesIn,
                              const std::map<std::string, RefCountPtr<N_DEV_InstanceBlock> > * devNamesIn);
    void printLineDiagnostics();
    void delayedPrintLineDiagnostics ();

    bool setupInitialConditions ( N_LAS_Vector & solnVec, N_LAS_Vector & flagVec);
    void outputDCOP(N_LAS_Vector & solnVec);

    NodeNamePairMap & getICData ( int &, std::string &);

    NodeNamePairMap & getAllNodes ( );

    // Runs specified RESULT output commands
    void outputRESULT (N_LAS_Vector * solnVecPtr, N_LAS_Vector * stateVecPtr, N_LAS_Vector * storeVecPtr);

    // // Finishes STEP/RESULT output
    void finishOutputSTEP ();

    void outputAC(double freq, const N_LAS_Vector * freqDomainSolnVecReal, const N_LAS_Vector * freqDomainSolnVecImaginary);
    void outputMPDE (double time, const N_LAS_Vector * solnVecPtr );
    void outputHomotopy(const std::vector<string> & parameter_names, const std::vector<double> & param_values, const N_LAS_Vector * solution_vector);
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

    void getMeasureValue (const std::string &name, double &result, bool &found);

    // Additional Public Declarations

    std::vector<std::string> getVariableNames();

    void setExternalNetlistParams( std::vector< pair< std::string, std::string> > & externalNetlistParams );

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

    // Utility functions for all the output functions.
    // set the context data for the N_UTL_Param object for later look up by getPrintValue()
    bool setParamContextType_(ParameterList::iterator iterParam);

    // // not all calling points for getPrintValue_ will have a state/store vector (i.e.
    // // in block analysis.  Thus, let these pointers default to null when not given.
    double getPrintValue(ParameterList::const_iterator it_tpL, const N_LAS_Vector * solnVecPtr,
                         const N_LAS_Vector * stateVecPtr=NULL, const N_LAS_Vector * storeVecPtr=NULL,
                         const N_LAS_Vector * solnVecImagPtr=NULL );

    N_PDS_Comm * getCommPtr() {
      return pdsCommPtr_;
    }

    int getProcID() const {
      return pdsCommPtr_->procID();
    }

    int getNumProcs() const {
      return pdsCommPtr_->numProc();
    }

    const std::string &getTitle() const {
      return title_;
    }

    void setTitle(const std::string &title) {
      title_ = title;
    }

    const std::string &getNetListFilename() const {
      return netListFilename_;
    }

    const N_IO_CmdParse &getCommandLine() const {
      return commandLine_;
    }

    const std::vector<N_ANP_SweepParam> &getStepParamVec() const {
      return stepParamVec_;
    }

    int getStepLoopNumber() const {
      return stepLoopNumber_;
    }

    int getMaxParamSteps() const {
      return maxParamSteps_;
    }

    const std::vector<N_ANP_SweepParam> &getDCParamVec() const {
      return dcParamVec_;
    }

    int getDCLoopNumber() const {
      return dcLoopNumber_;
    }

    int getMaxDCSteps() const {
      return maxDCSteps_;
    }

    PrintType::PrintType getPrintType() const {
      return PRINTType_;
    }

    bool getTempSweepFlag() const {
      return tempSweepFlag_;
    }

    const ParameterList *getVariableList() const {
      return printParameters_ ? &printParameters_->variableList_ : 0;
    }

    double getPRINTDCvalue() const {
      return PRINTdcvalue_;
    }

    double getPRINTDCstart() const {
      return PRINTdcstart_;
    }

    double getPRINTDCstop() const {
      return PRINTdcstop_;
    }

    const N_ANP_AnalysisInterface *getAnaIntPtr() const {
      return anaIntPtr_;
    }

    const std::string &getPRINTDCname() const {
      return PRINTdcname_;
    }

    bool getPrintEndOfSimulationLine() const {
      return printEndOfSimulationLine_;
    }

    void setEnableHomotopyFlag(bool value) {
      enableHomotopyFlag_ = value;
    }

    const N_TOP_Topology *getTopPtr() const {
      return topPtr_;
    }

    double getCircuitTime() const {
      return circuitTime_;
    }

    void setCircuitTime(double time) {
      circuitTime_ = time;
    }

    double getCircuitFrequency() const {
      return circuitFrequency_;
    }

    void setCircuitFrequency(double frequency) {
      circuitFrequency_ = frequency;
    }

    double getCircuitTemp() const {
      return circuitTemp_;
    }

    double getFilter() const {
      return filter_;
    }

    bool getSTEPEnabledFlag() const {
      return STEPEnabledFlag_;
    }

    const std::string &getFilenameSuffix() const {
      return filenameSuffix_;
    }

    Format::Format getFormat() const {
      return format_;
    }

    bool getNoIndex() const {
      return noIndex_;
    }

    bool getPRINTcsv() const {
      return format_ == Format::CSV;
    }

    bool getPRINTprobe() const {
      return format_ == Format::PROBE;
    }

    bool getPRINTtecplot() const {
      return format_ == Format::TECPLOT;
    }

    const ObjectiveMap &getObjectiveMap() const {
      return objective;
    }

    ObjectiveMap &getObjectiveMap() {
      return objective;
    }

    const N_MPDE_Manager *getMpdeMgrPtr() const {
      return mpdeMgrPtr_;
    }

    void pushActiveOutputters() {
      activeOutputterStack_.push_back(std::vector<Outputter::Interface *>());
    }

    void popActiveOutputters() {
      activeOutputterStack_.pop_back();
    }

    void addActiveOutputter(PrintType::PrintType print_type) {
      OutputterMap::iterator it = outputterMap_.find(print_type);

      if (it != outputterMap_.end())
        activeOutputterStack_.back().push_back((*it).second);
    }

    void setCurrentOutputter(Outputter::Interface *outputter) {
      currentOutputter_ = outputter;
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
    OutputMgr( N_IO_CmdParse & cp);

    void outputDCOP_restartFile ( N_LAS_Vector & solnVec);
    bool inputDCOP_( N_LAS_Vector & solnVec, N_LAS_Vector & flagVec);

    bool setupIC_or_NODESET( N_LAS_Vector & solnVec, N_LAS_Vector & flagVec,
                             bool & useFlag,
                             std::string & icType,
                             std::vector<N_UTL_OptionBlock> & initBlockVec);

    void outputIC_or_NODESET (N_LAS_Vector & solnVec);

    std::set<std::string> deviceParamCheck;

    void gatherGlobalDeviceCount_ (
      std::map<std::string,int> & globalDeviceCount , const std::map<std::string,int> & localDeviceCount );

    void formatPrintDeviceCount_ ( const std::map<std::string,int> & deviceCountMap_ , std::string & msg);

    void outputLocalDeviceCount_ ( const std::map<std::string,int> & localDeviceCount );

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
    // each derived type from N_IO_OutputFileBase handles a particular
    // format of output like standard, probe, tecplot etc.
    // Store a vector of them so that more than one type can be output at the
    // same time
    std::string                 netListFilename_;
    std::string                 title_;
    std::string                 filenameSuffix_;

    std::vector< N_IO_OutputFileBase *> outputHandlersVec_;

    std::string responseFileName_;
    bool responseFileNameGiven_;

    // Pointer to the ostream for .result output
    std::ostream *              resultStreamPtr_;

    // option to output values less than "filter_" as zero.
    double                      filter_;
    bool                        filterGiven_;

    std::map<std::string,int>   localDeviceCountMap_;

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
    const std::map<std::string, RefCountPtr<N_DEV_InstanceBlock> > *devNames_;

    // these variables used to be static:
    unsigned long staticIndex_;
    unsigned long staticMPDEIndex_;

    bool outputPrnAdded_;

    N_PDS_Comm *                pdsCommPtr_;
    N_TOP_Topology *            topPtr_;
    N_DEV_DeviceInterface *     devPtr_;
    N_ANP_AnalysisInterface *   anaIntPtr_;
    N_MPDE_Manager *            mpdeMgrPtr_;

    // print statement vars
    bool rawFlag_;
    bool tranFlag_;
    bool probeFlag_;
    bool DCSweepFlag_;
    bool ACFlag_;
    bool enableHomotopyFlag_;
    bool homotopyFlag_;
    bool HBFlag_;
    bool MPDEFlag_;
    bool PRINTflag_;

    PrintType::PrintType PRINTType_;

    Format::Format      format_;
    bool                noIndex_;

    PrintParameters *   printParameters_;

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
    std::vector<N_UTL_OptionBlock> ICblockVec_;
    bool nodesetflag_;
    std::vector<N_UTL_OptionBlock> nodesetblockVec_;

    double initialOutputInterval_;
    vector < pair < double, double > > outputIntervalPairs_;

    set < std::string > distribVsrcSet_;
    set < std::string > distribVnodeSet_;

    bool STEPEnabledFlag_;
    int  STEPcounter_;
    bool RESULTinitialized_;

    bool printEndOfSimulationLine_;  // flag to indicate if user wants the "End of Xyce(TM)" line in the output.

    // response functions requested in the external params passed to Xyce.
    std::vector< pair< std::string, std::string> > variablesUsedInSimulation_ ;
    std::vector< pair< std::string, std::string> > responseFunctionsRequested_ ;
    std::vector< pair< std::string, std::string> > derivativeVariablesRequested_ ;
    std::vector< pair< std::string, std::string> > analysisComponentsRequested_ ;

    // response related variable lists
    std::vector<ParameterList > responseVarList_;
    RCP< std::vector< double > > responseVarPtr_;
    std::vector<std::string> responseNames_;
    int numResponseVars_;

    bool tempSweepFlag_;
    bool outputCalledBefore_;

    // vector of expressions:
    std::vector<N_UTL_ExpressionData*> expVector_;

    // vector of expressions from .RESULT
    std::vector<N_UTL_ExpressionData*> resultVector_;

    ObjectiveMap objective;

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
    N_IO_CmdParse & commandLine_;

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

    N_IO_MeasureManager         measureManager_;
    N_IO_FourierMgr             fourierManager_;

    OpenPathStreamMap           openPathStreamMap_;
};

} // namespace IO
} // namespace Xyce

typedef Xyce::IO::OutputMgr N_IO_OutputMgr;

#endif

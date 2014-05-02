//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright(C) 2002-2014 Sandia Corporation
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
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
// Filename       : $RCSfile: N_IO_OutputMgr.C,v $
//
// Purpose        : Output Manager
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
// Revision Number: $Revision: 1.419.2.2 $
//
// Revision Date  : $Date: 2014/03/18 21:43:40 $
//
// Current Owner  : $Author: erkeite $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_UTL_Misc.h>

#include <iostream>
#include <fstream>
#include <sstream>

#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif

#include <N_ANP_AnalysisInterface.h>
#include <N_ANP_AnalysisManager.h>
#include <N_DEV_DeviceInterface.h>
#include <N_DEV_DeviceSensitivities.h>
#include <N_ERH_ErrorMgr.h>
#include <N_IO_CmdParse.h>
#include <N_IO_FourierMgr.h>
#include <N_IO_Objective.h>
#include <N_IO_Op.h>
#include <N_IO_OutputFileBase.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_OutputPrn.h>
#include <N_IO_mmio.h>
#include <N_LAS_BlockVector.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_System.h>
#include <N_LAS_Vector.h>
#include <N_MPDE_Manager.h>
#include <N_PDS_Comm.h>
#include <N_PDS_ParMap.h>
#include <N_PDS_Serial.h>
#include <N_PDS_MPI.h>
#include <N_TOP_Topology.h>
#include <N_UTL_Algorithm.h>
#include <N_UTL_Expression.h>
#include <N_UTL_ExpressionData.h>
#include <N_UTL_LogStream.h>
#include <N_UTL_Version.h>

#include <Teuchos_as.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_oblackholestream.hpp>

#undef HAVE_LIBPARMETIS
#include <EpetraExt_RowMatrixOut.h>
#include <Epetra_CrsMatrix.h>

#ifdef Xyce_USE_HDF5
#ifdef Xyce_PARALLEL_MPI
#include <N_PDS_ParComm.h>
#include <N_PDS_MPIComm.h>
#endif // Xyce_PARALLEL_MPI
#include <Epetra_Map.h>
#include <Epetra_Comm.h>
#include <Epetra_MultiVector.h>
#include <N_IO_OutputHDF5.h>
#endif // Xyce_USE_HDF5

namespace Xyce {
namespace IO {

void gatherGlobalDeviceCount(Parallel::Machine comm, std::map<std::string,int> &globalDeviceMap, const std::map<std::string,int> &local_device_count_map );
std::ostream &printDeviceCount(std::ostream &os, const std::map<std::string, int> &device_count_map);

//-----------------------------------------------------------------------------
// Function      : OutputMgr::factory
// Purpose       : factory function
//
// Special Notes : ERK.  10/16/2005.  This used to be a singleton(ie a
//                 static pointer was returned) but had to be changed
//                 so that the library version of Xyce would work
//                 correctly.
//
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/10/00
//-----------------------------------------------------------------------------
OutputMgr * OutputMgr::factory(CmdParse & cp)
{
  OutputMgr * OM_ptr = new OutputMgr(cp);
  return OM_ptr;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::OutputMgr
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/10/00
//-----------------------------------------------------------------------------
OutputMgr::OutputMgr(CmdParse & cp)
  : commandLine_(cp),
    ACFlag_(false),
    rawFlag_(false),
    tranFlag_(false),
    MPDEFlag_(false),
    DCSweepFlag_(false),
    HBFlag_(false),
    homotopyFlag_(false),
    sensitivityFlag_(false),
    enableHomotopyFlag_(false),
    enableSensitivityFlag_(false),
    format_(Format::STD),
    noIndex_(false),
    formatResult_(Format::STD),
    printParameters_(&defaultPrintParameters_),
    PRINTType_(PrintType::NONE),
    PRINTdcstart_(0.0),
    PRINTdcstop_(0.0),
    PRINTdcvalue_(0.0),
    PRINTdcname_(""),
    ICflag_(false),
    saveFlag_(false),
    loadFlag_(false),
    outputOnceAlreadyFlag_(false),
    nodesetflag_(false),
    stepLoopNumber_(0),
    maxParamSteps_(0),
    dcLoopNumber_(0),
    maxDCSteps_(0),
    circuitTime_(0.0),
    circuitTemp_(27.0),
    circuitFrequency_(0.0),
    output_op_(false),
    input_op_(false),
    output_op_file_(""),
    input_op_file_(""),
    saveFileLevel_("ALL"),
    saveFileType_(".NODESET"),
    saveOutputFile_(""),
    topology_(0),
    deviceInterface_(0),
    analysisInterface_(0),
    initialOutputInterval_(0.0),
    staticIndex_(0),
    netListFilename_(cp.getArgumentValue("netlist")),
    title_(cp.getArgumentValue("netlist")),
    measureManager_(*this),
    fourierManager_(*this),
    filenameSuffix_(""),            // set using an inlined accessor function
    responseFileName_(""),          // set using inlined accessor function
    responseFileNameGiven_(false),
    tempSweepFlag_(false),
    outputCalledBefore_(false),
    STEPEnabledFlag_(false),
    STEPcounter_(0),
    printEndOfSimulationLine_(true),
    outputVersionInRawFile_(false),
    RESULTinitialized_(false),
    resultStreamPtr_(0),
    op_found_(0),
    total_soln_(0),
    numResponseVars_(0),
    detailedDeviceCountFlag_(false),
    pdsCommPtr_(0),
    icType_(-1),
    hdf5FileNameGiven_(false),
    hdf5HeaderWritten_(false),
    hdf5IndexValue_(0),
    sensObjFunction_(""),
    sensObjFuncGiven_(false)
{
  if (commandLine_.getArgumentValue("-delim") == "TAB")
    defaultPrintParameters_.delimiter_ = "\t";
  else if (commandLine_.getArgumentValue("-delim") == "COMMA")
    defaultPrintParameters_.delimiter_ = ",";
  else
    defaultPrintParameters_.delimiter_ = commandLine_.getArgumentValue("-delim");

  defaultPrintParameters_.rawOverride_ = commandLine_.argExists("-r");

  // If the output file is specified on the command line it takes
  // precedence, reset the value of netListFilename_
  if (commandLine_.argExists("-o"))
  {
    defaultPrintParameters_.filename_ = commandLine_.getArgumentValue("-o");
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::~OutputMgr
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/10/00
//-----------------------------------------------------------------------------
OutputMgr::~OutputMgr()
{
  for (std::vector< OutputFileBase *>::iterator it = outputHandlersVec_.begin(); it != outputHandlersVec_.end(); ++it)
    delete (*it);

  for (OutputterMap::iterator it = outputterMap_.begin(); it != outputterMap_.end(); ++it)
    delete (*it).second;

  for (OpenPathStreamMap::iterator it = openPathStreamMap_.begin(); it != openPathStreamMap_.end(); ++it)
    delete (*it).second.second;

  for (Util::OpList::iterator it = responseVarList_.begin(); it != responseVarList_.end(); ++it)
    delete *it;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::openFile
// Purpose       : open named file in given mode, create stream
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 07/08/2013
//-----------------------------------------------------------------------------
std::ostream *OutputMgr::openFile(const std::string &path, std::ios_base::openmode mode)
{
  OpenPathStreamMap::iterator it = openPathStreamMap_.find(path);

  if (path == "CONSOLE")
    return &Xyce::dout();
  else if (it != openPathStreamMap_.end()) {
    ++(*it).second.first;
    return (*it).second.second;
  }
  else {
    std::ostream *os = new std::ofstream(path.c_str(), mode);
    openPathStreamMap_[path] = std::pair<int, std::ostream *>(1, os);

    if (!os->good())
    {
      Report::UserFatal0() << "Failure opening " << path;
    }

    return os;
  }

}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::openFile
// Purpose       : open named file for output only, create stream
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 07/08/2013
//-----------------------------------------------------------------------------
std::ostream *OutputMgr::openFile(const std::string &path)
{
  return openFile(path, std::ios_base::out);
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::openBinaryFile
// Purpose       : open named file in binary mode for output only, create 
//                 stream
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 07/08/2013
//-----------------------------------------------------------------------------
std::ostream *OutputMgr::openBinaryFile(const std::string &path)
{
  return openFile(path, std::ios_base::out | std::ios_base::binary);
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::closeFile
// Purpose       : Close given stream
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 07/08/2013
//-----------------------------------------------------------------------------
int OutputMgr::closeFile(std::ostream *os)
{
  if (os == &Xyce::dout())
    return 1;

  int open_count = 0;

  for (OpenPathStreamMap::iterator it = openPathStreamMap_.begin(); it != openPathStreamMap_.end(); ++it) {
    if ((*it).second.second == os) {
      open_count = --(*it).second.first;
      if (open_count == 0) {
        delete os;
        openPathStreamMap_.erase(it);
        break;
      }
    }
  }

  return open_count;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::prepareOutput
// Purpose       : checks to make sure requested values can be output
// Special Notes : Primary purpose is for error checking BEFORE calculation.
//                 Also allows any name translation for output, if needed.
//
//                 This function also now serves as an initialization for
//                 allNodes_.  allNodes_ needs to be set up after topology
//                 is done setting up.
//
// erkeite  8/10/2007:
// Note that this is currently(august 2007) one of 3 different functions
// that check the .print line for errors.  That is a bit confusing,
// unfortunately.  For commentary and some explanation of this function,
// as well as the other two, see the comments for the delayedPrintLineDiagnostics
// function.  This whole functionality really should be refactored, or
// at least made more transparent to other developers.
//
// Function renamed from check_output to prepareOutput and heavily reorganized
// by David Baur on 6/28/2013.
//
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 06/01/05
//-----------------------------------------------------------------------------
void OutputMgr::prepareOutput(
  Analysis::Analysis_Mode               analysis_mode,
  const std::vector<N_ANP_SweepParam> & step_sweep_parameters,
  const std::vector<N_ANP_SweepParam> & dc_sweep_parameters)
{
  fixupNodeNames();

  setSweepParameters(step_sweep_parameters, dc_sweep_parameters);

  // Setup rawfile if requested
  if (commandLine_.argExists("-r"))
  {
    if (activeOutputterStack_.empty())
      activeOutputterStack_.push_back(std::vector<Outputter::Interface *>());

    enableOverrideRawOutput(*printParameters_);
  }
  else {
    switch (analysis_mode) {
      case Analysis::ANP_MODE_DC_OP:
        break;

      case Analysis::ANP_MODE_DC_SWEEP:
        enableDCSweepOutput();
        break;

      case Analysis::ANP_MODE_TRANSIENT:
        enableTransientOutput();
        break;

      case Analysis::ANP_MODE_MPDE:
//      enableTransientOutput(step_sweep_parameters, dc_sweep_parameters);
        enableMPDEOutput();
        break;

      case Analysis::ANP_MODE_HB:
        enableHBOutput();
        break;

      case Analysis::ANP_MODE_AC:
        enableACOutput();
        break;

      case Analysis::ANP_MODE_MOR:
        break;

      default:     // silence warnings from clang about uncaught cases in enum
        break;
    }

    if (enableHomotopyFlag_)
    {
      enableHomotopyOutput(analysis_mode);
    }

    if (enableSensitivityFlag_)
    {
      enableSensitivityOutput(analysis_mode);
    }
  }

  if (outputterMap_.empty()) {
    Outputter::TimePrn *outputter = new Outputter::TimePrn(*this, *printParameters_);
    outputter->parse();
    outputterMap_[PrintType::TRAN] = outputter;
  }

  measureManager_.fixupMeasureParameters();
  fourierManager_.fixupFourierParameters();
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::registerParallelServices
// Purpose       : Registers N_PDS_Comm object for parallel communication
// Special Notes : inline
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/10/00
//-----------------------------------------------------------------------------
bool OutputMgr::registerParallelServices(N_PDS_Comm * tmp_pds_ptr)
{
  pdsCommPtr_ = tmp_pds_ptr;

  return pdsCommPtr_;
}

//-----------------------------------------------------------------------------
// Namespace     : Unnamed
// Purpose       : file-local scoped methods and data
// Special Notes :
// Creator       : dgbaur
// Creation Date : 06/28/2013
//-----------------------------------------------------------------------------
namespace {

//-----------------------------------------------------------------------------
// Class         : OutputMgr_STEPOptionsReg
// Purpose       : functor for registering STEP options
// Special Notes : Used by package manager submitRegistration method
// Creator       : dgbaur
// Creation Date : 06/28/2013
//-----------------------------------------------------------------------------
struct OutputMgr_STEPOptionsReg : public PkgOptionsReg
{
  OutputMgr_STEPOptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const Util::OptionBlock & options )
  { return outputManager_.registerSTEPOptions( options ); }

  OutputMgr &outputManager_;
};



//-----------------------------------------------------------------------------
// Class         : OutputMgr_DCOptionsReg
// Purpose       : functor for registering DC options
// Special Notes : Used by package manager submitRegistration method
// Creator       : dgbaur
// Creation Date : 06/28/2013
//-----------------------------------------------------------------------------
struct OutputMgr_DCOptionsReg : public PkgOptionsReg
{
  OutputMgr_DCOptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const Util::OptionBlock & options )
  { return outputManager_.registerDCOptions( options ); }

  OutputMgr &outputManager_;
};



//-----------------------------------------------------------------------------
// Class         : OutputMgr_DCOPOptionsReg
// Purpose       : functor for registering DCOP options
// Special Notes : Used by package manager submitRegistration method
// Creator       : dgbaur
// Creation Date : 06/28/2013
//-----------------------------------------------------------------------------
struct OutputMgr_DCOPOptionsReg : public PkgOptionsReg
{
  OutputMgr_DCOPOptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const Util::OptionBlock & options )
  { return outputManager_.registerDCOPOptions( options ); }

  OutputMgr &outputManager_;
};



//-----------------------------------------------------------------------------
// Class         : OutputMgr_TranOptionsReg
// Purpose       : functor for registering transient options
// Special Notes : Used by package manager submitRegistration method
// Creator       : dgbaur
// Creation Date : 06/28/2013
//-----------------------------------------------------------------------------
struct OutputMgr_TranOptionsReg : public PkgOptionsReg
{
  OutputMgr_TranOptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const Util::OptionBlock & options )
  { return outputManager_.registerTranOptions( options ); }

  OutputMgr &outputManager_;
};



//-----------------------------------------------------------------------------
// Class         : OutputMgr_TranOptionsReg
// Purpose       : functor for registering transient options
// Special Notes : Used by package manager submitRegistration method
// Creator       : dgbaur
// Creation Date : 06/28/2013
//-----------------------------------------------------------------------------
struct OutputMgr_MPDETranOptionsReg : public PkgOptionsReg
{
  OutputMgr_MPDETranOptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const Util::OptionBlock & options )
  { return outputManager_.registerMPDETranOptions( options ); }

  OutputMgr &outputManager_;
};



//-----------------------------------------------------------------------------
// Class         : OutputMgr_HBOptionsReg
// Purpose       : functor for registering HB options
// Special Notes : Used by package manager submitRegistration method
// Creator       : dgbaur
// Creation Date : 06/28/2013
//-----------------------------------------------------------------------------
struct OutputMgr_HBOptionsReg : public PkgOptionsReg
{
  OutputMgr_HBOptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const Util::OptionBlock & options )
  { return outputManager_.registerHBOptions( options ); }

  OutputMgr &outputManager_;
};



//-----------------------------------------------------------------------------
// Class         : OutputMgr_OptionsReg
// Purpose       : functor for registering Output options
// Special Notes : Used by package manager submitRegistration method
// Creator       : dgbaur
// Creation Date : 06/28/2013
//-----------------------------------------------------------------------------
struct OutputMgr_OptionsReg : public PkgOptionsReg
{
  OutputMgr_OptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const Util::OptionBlock & options )
  { return outputManager_.registerOutputOptions( options ); }

  OutputMgr &outputManager_;
};



//-----------------------------------------------------------------------------
// Class         : OutputMgr_DeviceOptionsReg
// Purpose       : functor for registering Device options
// Special Notes : Used by package manager submitRegistration method
// Creator       : dgbaur
// Creation Date : 06/28/2013
//-----------------------------------------------------------------------------
struct OutputMgr_DeviceOptionsReg : public PkgOptionsReg
{
  OutputMgr_DeviceOptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const Util::OptionBlock & options )
  { return outputManager_.registerDeviceOptions( options ); }

  OutputMgr &outputManager_;
};



//-----------------------------------------------------------------------------
// Class         : OutputMgr_PrintOptionsReg
// Purpose       : functor for registering Print options
// Special Notes : Used by package manager submitRegistration method
// Creator       : dgbaur
// Creation Date : 06/28/2013
//-----------------------------------------------------------------------------
struct OutputMgr_PrintOptionsReg : public PkgOptionsReg
{
  OutputMgr_PrintOptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const Util::OptionBlock & options )
  { return outputManager_.registerPRINTSet
      ( options ); }

  OutputMgr &outputManager_;
};



//-----------------------------------------------------------------------------
// Class         : OutputMgr_ICOptionsReg
// Purpose       : functor for registering IC options
// Special Notes : Used by package manager submitRegistration method
// Creator       : dgbaur
// Creation Date : 06/28/2013
//-----------------------------------------------------------------------------
struct OutputMgr_ICOptionsReg : public PkgOptionsReg
{
  OutputMgr_ICOptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const Util::OptionBlock & options )
  { return outputManager_.registerIC ( options ); }

  OutputMgr &outputManager_;
};



//-----------------------------------------------------------------------------
// Class         : OutputMgr_NodeSetOptionsReg
// Purpose       : functor for registering NodeSet options
// Special Notes : Used by package manager submitRegistration method
// Creator       : dgbaur
// Creation Date : 06/28/2013
//-----------------------------------------------------------------------------
struct OutputMgr_NodeSetOptionsReg : public PkgOptionsReg
{
  OutputMgr_NodeSetOptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const Util::OptionBlock & options )
  { return outputManager_.registerNodeSet( options ); }

  OutputMgr &outputManager_;
};



//-----------------------------------------------------------------------------
// Class         : OutputMgr_ResultOptionsReg
// Purpose       : functor for registering Result options
// Special Notes : Used by package manager submitRegistration method
// Creator       : dgbaur
// Creation Date : 06/28/2013
//-----------------------------------------------------------------------------
struct OutputMgr_ResultOptionsReg : public PkgOptionsReg
{
  OutputMgr_ResultOptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const Util::OptionBlock & options )
  { return outputManager_.setRESULTParams( options ); }

  OutputMgr &outputManager_;
};



//-----------------------------------------------------------------------------
// Class         : OutputMgr_ObjectiveOptionsReg
// Purpose       : functor for registering Objective options
// Special Notes : Used by package manager submitRegistration method
// Creator       : dgbaur
// Creation Date : 06/28/2013
//-----------------------------------------------------------------------------
struct OutputMgr_ObjectiveOptionsReg : public PkgOptionsReg
{
  OutputMgr_ObjectiveOptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const Util::OptionBlock & options )
  { return outputManager_.setOBJECTIVEParams( options ); }

  OutputMgr &outputManager_;
};



//-----------------------------------------------------------------------------
// Class         : OutputMgr_SaveOptionsReg
// Purpose       : functor for registering Save options
// Special Notes : Used by package manager submitRegistration method
// Creator       : dgbaur
// Creation Date : 06/28/2013
//-----------------------------------------------------------------------------
struct OutputMgr_SaveOptionsReg : public PkgOptionsReg
{
  OutputMgr_SaveOptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const Util::OptionBlock & options )
  { return outputManager_.registerSave( options ); }

  OutputMgr &outputManager_;
};



//-----------------------------------------------------------------------------
// Class         : OutputMgr_LoadOptionsReg
// Purpose       : functor for registering Load options
// Special Notes : Used by package manager submitRegistration method
// Creator       : dgbaur
// Creation Date : 06/28/2013
//-----------------------------------------------------------------------------
struct OutputMgr_LoadOptionsReg : public PkgOptionsReg
{
  OutputMgr_LoadOptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const Util::OptionBlock & options )
  { return outputManager_.registerLoad( options ); }

  OutputMgr &outputManager_;
};



//-----------------------------------------------------------------------------
// Class         : OutputMgr_OPOptionsReg
// Purpose       : functor for registering OP options
// Special Notes : Used by package manager submitRegistration method
// Creator       : dgbaur
// Creation Date : 06/28/2013
//-----------------------------------------------------------------------------
struct OutputMgr_OPOptionsReg : public PkgOptionsReg
{
  OutputMgr_OPOptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const Util::OptionBlock & options )
  { return outputManager_.registerOP ( options ); }

  OutputMgr &outputManager_;
};



//-----------------------------------------------------------------------------
// Class         : OutputMgr_MeasureOptionsReg
// Purpose       : functor for registering Measure options
// Special Notes : Used by package manager submitRegistration method
// Creator       : dgbaur
// Creation Date : 06/28/2013
//-----------------------------------------------------------------------------
struct OutputMgr_MeasureOptionsReg : public PkgOptionsReg
{
  OutputMgr_MeasureOptionsReg( Measure::Manager &mgr )
    : measureManager_(mgr)
  {}

  bool operator()( const Util::OptionBlock & options )
  { return measureManager_.addMeasure( options ); }

  Measure::Manager &measureManager_;
};

//-----------------------------------------------------------------------------
// Class         : OutputMgr_FourierOptionsReg
// Purpose       : functor for registering Fourier options
// Special Notes : Used by package manager submitRegistration method
// Creator       : dgbaur
// Creation Date : 06/28/2013
//-----------------------------------------------------------------------------
struct OutputMgr_FourierOptionsReg : public PkgOptionsReg
{
  OutputMgr_FourierOptionsReg( FourierMgr &mgr )
    : fourierManager_(mgr)
  {}

  bool operator()( const Util::OptionBlock & options )
  { return fourierManager_.addFourierAnalysis( options ); }

  FourierMgr &fourierManager_;
};

//-----------------------------------------------------------------------------
// Class         : OutputMgr_SensOptionsReg
// Purpose       : functor for registering sensitivity options
// Special Notes : Used by package manager submitRegistration method
// Creator       : Eric Keiter
// Creation Date : 02/10/2014
//-----------------------------------------------------------------------------
struct OutputMgr_SensOptionsReg : public PkgOptionsReg
{
  OutputMgr_SensOptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const Util::OptionBlock & options )
  { return outputManager_.registerSens( options ); }

  OutputMgr &outputManager_;
};

} // namespace <unnamed>

//-----------------------------------------------------------------------------
// Function      : OutputMgr::registerPkgOptionsMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, 1437
// Creation Date : 10/21/08
//-----------------------------------------------------------------------------
bool OutputMgr::registerPkgOptionsMgr(PkgOptionsMgr &pkgOpt)
{
  pkgOpt.submitRegistration(
    "DC", netListFilename_, new OutputMgr_DCOptionsReg(*this));

  pkgOpt.submitRegistration(
    "TRAN", netListFilename_, new OutputMgr_TranOptionsReg(*this));

  pkgOpt.submitRegistration(
    "MPDE", netListFilename_, new OutputMgr_MPDETranOptionsReg(*this));

  pkgOpt.submitRegistration(
    "HB", netListFilename_, new OutputMgr_HBOptionsReg(*this));

  pkgOpt.submitRegistration(
    "STEP", netListFilename_, new OutputMgr_STEPOptionsReg(*this));

  pkgOpt.submitRegistration(
    "OP_IO", netListFilename_, new OutputMgr_DCOPOptionsReg(*this));

  pkgOpt.submitRegistration(
    "OUTPUT", netListFilename_, new OutputMgr_OptionsReg(*this));

  pkgOpt.submitRegistration(
    "PRINT", netListFilename_, new OutputMgr_PrintOptionsReg(*this));

  pkgOpt.submitRegistration(
    "IC", netListFilename_, new OutputMgr_ICOptionsReg(*this));

  pkgOpt.submitRegistration(
    "NODESET", netListFilename_, new OutputMgr_NodeSetOptionsReg(*this));

  pkgOpt.submitRegistration(
    "RESULT", netListFilename_, new OutputMgr_ResultOptionsReg(*this));

  pkgOpt.submitRegistration(
    "OBJECTIVE", netListFilename_, new OutputMgr_ObjectiveOptionsReg(*this));

  pkgOpt.submitRegistration(
    "SAVE", netListFilename_, new OutputMgr_SaveOptionsReg(*this));

  pkgOpt.submitRegistration(
    "LOAD", netListFilename_, new OutputMgr_LoadOptionsReg(*this));

  pkgOpt.submitRegistration(
    "OP", netListFilename_, new OutputMgr_OPOptionsReg(*this));

  pkgOpt.submitRegistration(
    "DEVICE", netListFilename_, new OutputMgr_DeviceOptionsReg(*this));

  pkgOpt.submitRegistration(
    "MEASURE", netListFilename_, new OutputMgr_MeasureOptionsReg(measureManager_));

  pkgOpt.submitRegistration(
    "FOUR", netListFilename_, new OutputMgr_FourierOptionsReg(fourierManager_));

  pkgOpt.submitRegistration(
      "SENS", netListFilename_, new OutputMgr_SensOptionsReg(*this)); 

  return true;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::registerDCOPOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 05/08/06
//-----------------------------------------------------------------------------
bool OutputMgr::registerDCOPOptions(const Util::OptionBlock & option_block)
{
  for (ParameterList::const_iterator iterPL = option_block.getParams().begin(); iterPL != option_block.getParams().end(); ++iterPL)
  {
    if (iterPL->tag() == "INPUT")
    {
      input_op_ = true;
      input_op_file_ = iterPL->stringValue();
    }
    else if (iterPL->tag() == "OUTPUT")
    {
      output_op_ = true;
      output_op_file_ = iterPL->stringValue();
    }
    else if (iterPL->tag() == "TIME")
    {
      // do nothing, this will be handled in the time integrator.
    }
    else
    {
      Report::UserWarning0() << "Parameter " << iterPL->tag() << " not recognized in .DCOP command";
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::registerDCOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 08/21/06
//-----------------------------------------------------------------------------
bool OutputMgr::registerDCOptions(const Util::OptionBlock & option_block)
{
  for (ParameterList::const_iterator 
       iterPL = option_block.getParams().begin(); 
       iterPL != option_block.getParams().end(); ++iterPL)
  {
    if (iterPL->tag() == "PARAM")
    {
      dcParams_.push_back(iterPL->stringValue());
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::registerTranOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 09/25/07
//-----------------------------------------------------------------------------
bool OutputMgr::registerTranOptions(const Util::OptionBlock & OB)
{
//  PRINTType_ = PrintType::TRAN;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::registerMPDETranOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, Rich Schiek, Ting Mei
// Creation Date : 7/24/08
//-----------------------------------------------------------------------------
bool OutputMgr::registerMPDETranOptions(const Util::OptionBlock & OB)
{
//  PRINTType_ = PrintType::TRAN;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::registerHBOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey, Rich Schiek, Ting Mei
// Creation Date : 7/24/08
//-----------------------------------------------------------------------------
bool OutputMgr::registerHBOptions(const Util::OptionBlock & OB)
{
//  PRINTType_ = PrintType::HB;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::registerSTEPOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 08/21/06
//-----------------------------------------------------------------------------
bool OutputMgr::registerSTEPOptions(const Util::OptionBlock & option_block)
{
  for (ParameterList::const_iterator 
      iterPL = option_block.getParams().begin(); 
      iterPL != option_block.getParams().end(); 
      ++iterPL)
  {
    if (iterPL->tag() == "PARAM")
    {
      stepParams_.push_back(iterPL->stringValue());
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::registerDeviceOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 5/10/10
//-----------------------------------------------------------------------------
bool OutputMgr::registerDeviceOptions(const Util::OptionBlock & option_block)
{
  for (ParameterList::const_iterator 
      iter = option_block.getParams().begin(); 
      iter != option_block.getParams().end(); 
      ++iter)
  {
    if (iter->tag() == "DETAILED_DEVICE_COUNTS")
    {
      detailedDeviceCountFlag_ = iter->getImmutableValue<bool>();
    }
  }

  return true;
}


//-----------------------------------------------------------------------------
// Function      : OutputMgr::registerOutputOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/30/01
//-----------------------------------------------------------------------------
bool OutputMgr::registerOutputOptions(const Util::OptionBlock & OB)
{

  // Need to capture a variable length of output intervals.  The use case is:
  //
  // .OPTIONS OUTPUT INITIAL_INTERVAL=<interval> [<t0> <i0> [<t1> <i1>...]]
  //
  // additional options can appear before or after the INITIAL_INTERVAL as in
  //
  //.OPTIONS OUTPUT HDF5FILE=xxxx INITIAL_INTERVAL=<interval> [<t0> <i0> [<t1> <i1>...]]
  // or
  // .OPTIONS OUTPUT INITIAL_INTERVAL=<interval> [<t0> <i0> [<t1> <i1>...]] HDF5FILE=xxxx
  //

  ParameterList::const_iterator iterPL = OB.getParams().begin();
  while (iterPL != OB.getParams().end())
  {
    if (iterPL->tag() == "INITIAL_INTERVAL")
    {
      // start handling a list of intervals
      // set interval value
      initialOutputInterval_ = iterPL->getImmutableValue<double>();
      // look for optional time pairs.
      outputIntervalPairs_.clear();
      bool doneWithTimePairs = false;
      
      // need to point at next parameter
      ++iterPL;

      while ((iterPL != OB.getParams().end()) && !doneWithTimePairs)
      {
        if (iterPL->tag() == "TIME")
        {
          double t = iterPL->getImmutableValue<double>();
          ++iterPL;
          double iv = iterPL->getImmutableValue<double>();
          ++iterPL;

          outputIntervalPairs_.push_back(std::pair<double, double>(t, iv));
        }
        else
        {
          // didn't find a time pair so bail out of this loop.
          doneWithTimePairs=true;
        }
      }
    }
    else if (iterPL->tag()=="HDF5FILENAME")
    {
    // look for other option tags
      hdf5FileNameGiven_=true;
      hdf5FileName_=iterPL->stringValue();
      ++iterPL;    
    } 
    else if (iterPL->tag()=="PRINTENDOFSIMLINE")
    {
      // look for flag to turn off "End of Xyce(TM) Simulation" line
      printEndOfSimulationLine_=iterPL->getImmutableValue<bool>();
      ++iterPL;    
    }
    else if (iterPL->tag()=="OUTPUTVERSIONINRAWFILE")
    {
      // look for flag to toggle output of version in header of RAW file 
     outputVersionInRawFile_=iterPL->getImmutableValue<bool>();
     ++iterPL;    
    }
    else
    {
      // silently ignore?
      ++iterPL;    
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::registerDistribNodes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 11/16/01
//-----------------------------------------------------------------------------
bool OutputMgr::registerDistribNodes(const Util::OptionBlock & option_block)
{
  for (ParameterList::const_iterator 
      iterPL = option_block.getParams().begin(); 
      iterPL != option_block.getParams().end(); 
      ++iterPL)
  {
    distribVnodeSet_.insert(iterPL->tag());
    distribVsrcSet_.insert(iterPL->stringValue());
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::registerNodeDevNames
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 08/18/06
//-----------------------------------------------------------------------------
void OutputMgr::registerNodeDevNames(
    const std::set<std::string> * nodeNamesIn,
    const std::map<std::string, RefCountPtr<N_DEV_InstanceBlock> > * devNamesIn)
{
  nodeNames_ = nodeNamesIn;
  devNames_ = devNamesIn;
}


//-----------------------------------------------------------------------------
// Function      : OutputMgr::setExternalNetlistParams
// Purpose       : Called from owning class to set any external parameters
//                 set by the user.  Note, these parameter lists may also
//                 request response functions that the output manger must report.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical Systems Modeling
// Creation Date : 08/07/2012
//-----------------------------------------------------------------------------
void OutputMgr::setExternalNetlistParams(
  std::vector< std::pair< std::string, std::string> > & externalNetlistParams)
{

  //  externalParams can contain section names from dakota like
  // variables 2, responses 4, derivatives 4.  If we don't find any sections tags
  // then we can assume that all the parameters are just variable names to set.
  // if we find tags, then use the "responses" section to record what response functions
  // need to be reported.
  std::string sectionTag="variables";
  std::vector< std::pair< std::string, std::string > >::iterator externalParamsIter = externalNetlistParams.begin();
  std::vector< std::pair< std::string, std::string > >::iterator externalParamsEnd = externalNetlistParams.end();
  while (externalParamsIter != externalParamsEnd)
  {
    if (externalParamsIter->first == "variables")
    {
      // found a section marker, use it to set the section
      sectionTag="variables";
    }
    else if (externalParamsIter->first == "functions")
    {
      // found a section marker, use it to set the section
      sectionTag="functions";
    }
    else if (externalParamsIter->first == "derivative_variables")
    {
      // found a section marker, use it to set the section
      sectionTag="derivative_variables";
    }
    else if (externalParamsIter->first == "analysis_components")
    {
      // found a section marker, use it to set the section
      sectionTag="analysis_components";
    }
    else
    {
      // found data so put it in the right container
      if (sectionTag=="variables")
      {
        variablesUsedInSimulation_.push_back(*externalParamsIter);
      }
      else if (sectionTag=="functions")
      {
        responseFunctionsRequested_.push_back(*externalParamsIter);
      }
      else if (sectionTag=="derivative_variables")
      {
        derivativeVariablesRequested_.push_back(*externalParamsIter);
      }
      else if (sectionTag=="analysis_components")
      {
        analysisComponentsRequested_.push_back(*externalParamsIter);
      }
    }
    externalParamsIter++;
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::printLineDiagnostics
// Purpose       :
// Special Notes : erkeite:  8/10/2007
//
// For commentary and some explanation of this function, see the comments
// for the delayedPrintLineDiagnostics function.
//
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 08/18/06
//-----------------------------------------------------------------------------
void OutputMgr::printLineDiagnostics()
{
  if (!printParameters_ || printParameters_->variableList_.empty())
  {
    return;
  }

  std::map<std::string, bool> stringStat;
  std::map<std::string, bool> nodeStat;
  std::map<std::string, bool> instanceStat;

  for (ParameterList::const_iterator iterParam = printParameters_->variableList_.begin() ; iterParam != printParameters_->variableList_.end(); ++iterParam)
  {
    if (true) // (*iterParam).getSimContext() == Util::UNDEFINED)
    {
      std::string varType(iterParam->tag());

      bool done = false;
      std::vector<std::string> nodes;
      std::vector<std::string> instances;
      std::vector<std::string> leads;
      std::vector<std::string> strings;
      std::vector<std::string> special;
      if (Util::hasExpressionTag(*iterParam) )
      {
        // parameter starts with "{" but may not have been been parsed into an expression.
        // check if there is an underlying expression object with the parameter
        if (iterParam->getType() == Util::EXPR)
        {
          iterParam->getValue<Util::Expression>().get_names(XEXP_NODE, nodes);
          iterParam->getValue<Util::Expression>().get_names(XEXP_INSTANCE, instances);
          iterParam->getValue<Util::Expression>().get_names(XEXP_LEAD, leads);
          iterParam->getValue<Util::Expression>().get_names(XEXP_STRING, strings);
          iterParam->getValue<Util::Expression>().get_names(XEXP_SPECIAL, special);  // special returns vars like TIME
        }
        else if ( ((iterParam->getType()) == Util::DBLE) || ((iterParam->getType()) == Util::INT) )
        {
        }
        instances.insert(instances.end(), leads.begin(), leads.end());
        strings.insert(strings.end(), special.begin(), special.end());  // add specials to strings
      }
      else
      {
        int numIndices = iterParam->getImmutableValue<int>();
        if (numIndices > 0 &&(varType == "I" ||(varType.size() == 2 && varType[0] == 'I')))
        {
          // any devices found in this I(xxx) structure need to be communicated to the device manager
          // so that the lead currents can be calculated
          if (numIndices != 1)
          {
            Report::DevelFatal0() << "Only one device argument allowed in I() in .PRINT command";
          }
          ++iterParam;
          if (varType.size() == 2)
          {
            instances.push_back(iterParam->tag() + "{" + varType[1] + "}");
          }
          else
          {
            instances.push_back(iterParam->tag());
          }
        }
        else if (numIndices > 0 &&(varType == "V" ||((varType.size() == 2 || varType.size() == 3) && varType[0] == 'V')))
        {
          if (numIndices < 1 || numIndices > 2)
          {
            Report::DevelFatal0() << "Only one or two node arguments allowed in V() in .PRINT command";
          }
          for (; numIndices > 0 ; --numIndices)
          {
            ++iterParam;
            nodes.push_back(iterParam->tag());
          }
        }
        else if (numIndices > 0 && varType == "N")
        {
          // Don't attempt to check this type of variable.
          //
          // N-variables can only be fully resolved once the devices are allocated,
          // and variables which are internal to the devices have been assigned names.
          //(N-variables are the same as V-variables, except that they can include
          // internal vars).
          if (numIndices != 1)
          {
            Report::DevelFatal0() << "Only one device argument allowed in N() in .PRINT command";
          }
          ++iterParam;
        }
        else
        {
          strings.push_back(iterParam->tag());
        }
      }

      for (std::vector<std::string>::iterator 
          iter_s = strings.begin(); 
          iter_s != strings.end(); ++iter_s)
      {
        done = false;

        if (*iter_s == "TEMP" || *iter_s == "TIME" || 
            *iter_s == "FREQUENCY" || *iter_s == "INDEX")
          done = true;
        std::vector<std::string>::iterator step_i;
        std::vector<std::string>::iterator step_end;
        if (!done)
        {
          step_i = dcParams_.begin();
          step_end = dcParams_.end();
          // check if this is a step sweep value.
          for (; step_i != step_end ; ++step_i)
          {
            if (*step_i == *iter_s)
            {
              done = true;
              break;
            }
          }
        }

        double mValue;
        bool mFound;
        measureManager_.getMeasureValue(*iter_s, mValue, mFound);
        if (mFound)
          done = true;

        if (!done)
        {
          step_i = stepParams_.begin();
          step_end = stepParams_.end();
          // check if this is a step sweep value.
          for (; step_i != step_end ; ++step_i)
          {
            if (*step_i == *iter_s)
            {
              done = true;
              break;
            }
          }
        }

        if (!done)
        {
          double result;
          done = deviceInterface_->getParamAndReduce(*iter_s, result);

          // Note:  when this is called, the devices haven't been allocated
          // yet.  So, this particular check isn't reliable.  If the device and/or
          // parameter is not found by the device package via getParam, then
          // simply save it for a later diagnostic.
          if (!done)
          {
            deviceParamCheck.insert(*iter_s);
          }
          done = true;
        }
        stringStat[*iter_s] = done;
      }


      for (std::vector<std::string>::iterator 
          iter_s = nodes.begin(); 
          iter_s != nodes.end(); ++iter_s)
      {
        bool tmpBool = false;
        // ERK: in rare cases, nodeNames_ can be NULL.  The function 
        // registerNodeDevNames is only called
        // on processors with greater than zero devices.  Occasionally, Xyce 
        // will give a processor nothing to do.
        if (nodeNames_ !=0)
        {
          // allow * to pass
          tmpBool = (*iter_s == "*") || (nodeNames_->find(*iter_s) != nodeNames_->end());
        }
        nodeStat[*iter_s] = tmpBool;
      }

      for (std::vector<std::string>::iterator iter_s = instances.begin(); iter_s != instances.end(); ++iter_s)
      {
        if ((*iter_s).substr((*iter_s).size()-1, 1) == "}")
        {
          bool tmpBool = false;
          // ERK: in rare cases, devNames_ can be NULL.  The function 
          // registerNodeDevNames is only called
          // on processors with greater than zero devices.  Occasionally, Xyce 
          // will give a processor nothing to do.
          if (devNames_ !=0)
          {
            // allow * to pass
            tmpBool = (*iter_s == "*") || (devNames_->find((*iter_s).substr(0, (*iter_s).size()-3)) != devNames_->end());
          }
          instanceStat[*iter_s] = tmpBool;
        }
        else
        {
          bool tmpBool = false;
          // ERK: in rare cases, devNames_ can be NULL.  The function 
          // registerNodeDevNames is only called
          // on processors with greater than zero devices.  Occasionally, Xyce 
          // will give a processor nothing to do.
          if (devNames_ !=0)
          {
            // allow * to pass
            tmpBool = (*iter_s == "*") || (devNames_->find(*iter_s) != devNames_->end());
          }
          instanceStat[*iter_s] = tmpBool;
        }
      }
    }
  }

  // Sum results if parallel
  if (!pdsCommPtr_->isSerial())
  {
    std::map<std::string, bool>::iterator stat_i, stat_end;

    int totSize = nodeStat.size() + instanceStat.size();
    std::vector<int> stat(totSize, 0);
    std::vector<int> stat_g(totSize, 0);
    int i=0;
    stat_i = nodeStat.begin();
    stat_end = nodeStat.end();
    for (; stat_i != stat_end ; ++stat_i)
    {
      if ((*stat_i).second)
      {
        stat[i++] = 1;
      }
      else
      {
        stat[i++] = 0;
      }
    }
    stat_i = instanceStat.begin();
    stat_end = instanceStat.end();
    for (; stat_i != stat_end ; ++stat_i)
    {
      if ((*stat_i).second)
      {
        stat[i++] = 1;
      }
      else
      {
        stat[i++] = 0;
      }
    }
    pdsCommPtr_->sumAll(&stat[0], &stat_g[0], totSize);
    i = 0;
    stat_i = nodeStat.begin();
    stat_end = nodeStat.end();
    for (; stat_i != stat_end ; ++stat_i)
    {
      if (stat_g[i++] > 0)
      {
        (*stat_i).second = true;
      }
    }
    stat_i = instanceStat.begin();
    stat_end = instanceStat.end();
    for (; stat_i != stat_end ; ++stat_i)
    {
      if (stat_g[i++] > 0)
      {
        (*stat_i).second = true;
      }
    }
  }

  // Generate message
  std::ostringstream oss;
  int count = 0;

  for (std::map<std::string, bool>::iterator 
      stat_i = nodeStat.begin(); 
      stat_i != nodeStat.end(); ++stat_i)
  {
    if (!(*stat_i).second)
    {
      if (count != 0)
      {
        oss << ", ";
      }
      oss << "node " << (*stat_i).first;
      ++count;
    }
  }

  for (std::map<std::string, bool>::iterator 
      stat_i = instanceStat.begin(); 
      stat_i != instanceStat.end(); ++stat_i)
  {
    if (!(*stat_i).second)
    {
      if (count != 0)
      {
        oss << ", ";
      }
      if (((*stat_i).first)[((*stat_i).first).size()-1] == '}')
      {
        oss << "I" << ((*stat_i).first).substr(((*stat_i).first).size()-2, 1)
            << "("
            << ((*stat_i).first).substr(0, ((*stat_i).first).size()-3)
            << ")";
        ++count;
      }
      else
      {
        oss << "I(" << (*stat_i).first << ")";
        ++count;
      }
    }
  }

  for (std::map<std::string, bool>::iterator 
      stat_i = stringStat.begin(); 
      stat_i != stringStat.end(); ++stat_i)
  {
    if (!(*stat_i).second)
    {
      if (count != 0)
      {
        oss << ", ";
      }
      oss << (*stat_i).first;
      ++count;
    }
  }

  if (count != 0)
  {
    Report::UserError0().at(printParameters_->netlistLocation_) << "There " 
      << (count == 1 ? "was " : "were ") << count << " undefined symbol" 
      << (count == 1 ? "" : "s") << " in .PRINT command: " << oss.str();
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::delayedPrintLineDiagnostics
// Purpose       : Yet another function that checks the .print line.
//
// Special Notes : erkeite.  8/10/2007
//
//                 This function has to be called after the devices in
//                 the device package have all been allocated.  The function
//                 printLineDiagnostics is called earlier than that, so
//                 it doesn't handle device-specific parameter outputs.
//
//                 We have 3 completely different functions that do different
//                 aspects of checking the print line:
//
//                (1) printLineDiagnostics:  happens really early, as soon as
//                 topology is constructed, but before devivces are allocated.
//
//                (2) delayedPrintLineDiagnostics:  happens just a little later, to
//                 check .print line fields that are specific to device entities(ie
//                 aren't just I and V solution vars).
//
//                (3) check_output:  Happens at the beginning of the simulation,
//                 after everything is completely set up, but before any solves
//                 happen.  The advantage of this one is that it actually uses
//                 all the same function calls as the regular output call that
//                 will happen later.  So, it is probably the most thorough test
//                 of all.  However, for really big simulations, it happens relatively
//                 late in the setup.  For really large simulations, setup can take
//                 a long time, and it is often useful to get diagnostics earlier
//                 than that.
//
//                 Note from TVR 10/22/2013:  "check_output", referred to
//                 above, has been broken apart and renamed "prepareOutput".
//                 The notes from ERK above are therefore somewhat out of date.
//                 The refactor, though, is not complete, and the replacement
//                 code is not well documented.
//
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/17/07
//-----------------------------------------------------------------------------
void OutputMgr::delayedPrintLineDiagnostics()
{
  if (deviceParamCheck.empty())
    return;

  std::ostringstream oss;
  int count = 0;
  for (std::set<std::string>::iterator 
      iter = deviceParamCheck.begin(); 
      iter != deviceParamCheck.end(); ++iter)
  {
    std::string s = *iter;

    double result = 0.0;
    bool done = deviceInterface_->getParamAndReduce(s, result);
    if (!done)
    {
      if (count != 0)
        oss << ", ";
      oss << *iter;
      ++count;
    }
  }

  if (count != 0)
  {
    Report::UserFatal0() << "There " << (count == 1 ? "was " : "were ") 
      << count << " undefined symbol" << (count == 1 ? "" : "s") 
      << " in .PRINT command: " << oss.str();
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::getOutputIntervals
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/30/01
//-----------------------------------------------------------------------------
bool OutputMgr::getOutputIntervals(
    double & initialInterval, 
    std::vector< std::pair<double, double> > & intervalPairs)
{
  initialInterval = initialOutputInterval_;
  intervalPairs = outputIntervalPairs_;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_IO OutputMgr::setRESULTParams
// Purpose       : Sets the RESULT calculation parameters.
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/29/04
//-----------------------------------------------------------------------------
bool OutputMgr::setRESULTParams(const Util::OptionBlock & OB)
{
#ifdef Xyce_DEBUG_IO
  Xyce::dout() << "In N_IO OutputMgr::setRESULTParams" << std::endl;
#endif

  ParameterList::const_iterator it_tp;
  ParameterList::const_iterator it_param;
  ParameterList::const_iterator it_type;
  ParameterList::const_iterator first = OB.getParams().begin();
  ParameterList::const_iterator last = OB.getParams().end();

  std::string msg;

  // first check to see that there is only 1 PARAM set.  They need to be in
  // separate lines for this to work.
  int countPar = 0;
  for (it_tp = first; it_tp != last; ++it_tp)
  {
    if (it_tp->uTag() == "EXPRESSION")
    {
      ++countPar;
    }
  }
  if (countPar > 1)
  {
    Report::DevelFatal0() << "Only one expression per .RESULT command.  Each parameter needs its own .RESULT line.";
  }

  Util::ExpressionData * expDataPtr;

  for (it_tp = first; it_tp != last; ++it_tp)
  {
    Util::Param expParam = *it_tp;

    if (!expParam.hasExpressionValue())
    {
      Report::DevelFatal0() << "Parameter must be an expression in .RESULT command";
    }
    else
    {
      // expression should have already been resolved.  Check for this
      // case before we try and create a new expression.
      if (expParam.getType() == Util::EXPR)
      {
        expDataPtr = new Util::ExpressionData(expParam.getValue<Util::Expression>(), *this);
      }
      else
      {
        expDataPtr = new Util::ExpressionData(expParam.stringValue(), *this);
      }
      resultVector_.push_back(expDataPtr);
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_IO OutputMgr::setOBJECTIVEParams
// Purpose       : Sets the OBJECTIVE calculation parameters.
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 11/07/05
//-----------------------------------------------------------------------------
bool OutputMgr::setOBJECTIVEParams(const Util::OptionBlock & OB)
{
#ifdef Xyce_DEBUG_IO
  Xyce::dout() << "In N_IO OutputMgr::setOBJECTIVEParams" << std::endl;
#endif

  ParameterList::const_iterator it_tp;
  ParameterList::const_iterator it_param;
  ParameterList::const_iterator it_type;
  ParameterList::const_iterator first = OB.getParams().begin();
  ParameterList::const_iterator last = OB.getParams().end();

  std::string name("");

  for (it_tp = first; it_tp != last; ++it_tp)
  {
    if (it_tp->uTag() == "NAME")
    {
      name = it_tp->usVal();
      if (objective_.find(name) != objective_.end())
      {
        Report::DevelFatal0() << "Duplicate objective name " << name;
      }
      break;
    }
  }
  objective_[name].initialize(OB, this);

  return true;
}


//-----------------------------------------------------------------------------
// Function      : OutputMgr::registerPRINTSet
// Purpose       : registers set of variables to print for .PRINT
//
// Special Notes : ERK.  This function used to also output the header.
//                 Now, that functionality is down in outputPRINTHeader_.
//
//                 Also, this function allocates expression handlers for
//                 any expressions specified.
//
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/10/00
//-----------------------------------------------------------------------------
bool OutputMgr::registerPRINTSet(const Util::OptionBlock &print_block)
{
//  TraceIO trace__("bool OutputMgr::registerPRINTSet(const Util::OptionBlock &print_block)");

  PrintParameters print_parameters;

  print_parameters.netlistLocation_ = print_block.getNetlistLocation();
  
  ParameterList::const_iterator iterParam = print_block.getParams().begin();
  for (; iterParam != print_block.getParams().end(); ++iterParam)
  {

#ifdef Xyce_DEBUG_IO
    Xyce::dout() << "iterParam->tag = " << iterParam->tag() << std::endl;
#endif

    if (iterParam->tag() == "WIDTH") {
      print_parameters.streamWidth_ = iterParam->getImmutableValue<int>();
    }

    else if (iterParam->tag() == "TYPE") {
      std::string s = iterParam->stringValue();
      if (s == "TRAN")
        PRINTType_ = PrintType::TRAN;
      else if (s == "AC")
        PRINTType_ = PrintType::AC;
      else if (s == "HB")
        PRINTType_ = PrintType::HB;
      else if (s == "DC")
        PRINTType_ = PrintType::DC;
      else
      {
        Report::DevelFatal0() << "Unrecognized analysis type " << s;
      }
      print_parameters.printType_ = PRINTType_;
    }

    else if (iterParam->tag() == "PRECISION") {
      print_parameters.streamPrecision_ = iterParam->getImmutableValue<int>();
    }

    else if (iterParam->tag() == "FILTER")
    {
      print_parameters.filter_ = iterParam->getImmutableValue<double>();
    }
    else if (iterParam->tag() == "FORMAT")
    {
      std::string s = iterParam->stringValue();

      if (s == "STD")
        format_ = Format::STD;
      else if (s == "TECPLOT")
        format_ = Format::TECPLOT;
      else if (s == "PROBE")
        format_ = Format::PROBE;
      else if (s == "NOINDEX") {
        format_ = Format::STD;
        noIndex_ = true;
        print_parameters.index_ = false;
      }
      else if (s == "CSV") {
        format_ = Format::CSV;
        noIndex_ = true;
        print_parameters.index_ = false;
        print_parameters.delimiter_ = ",";
      }
      else if (s == "RAW")
        format_ = Format::RAW;
      else
      {
        Report::DevelFatal0() << "Unrecognized print format " << s;
      }
      print_parameters.format_ = format_;
    }

    else if (iterParam->tag() == "TIMEWIDTH")
    {
      print_parameters.timeWidth_ = iterParam->getImmutableValue<int>();
    }

    else if (iterParam->tag() == "TIMESCALEFACTOR")
    {
      print_parameters.outputTimeScaleFactor_ = iterParam->getImmutableValue<double>();
    }

    else if (iterParam->tag() == "FILE")
    {
      // netListFilename_ should be the default unless FILE was set to
      // something other than "" which is its default value.
      if (iterParam->stringValue() != "")
      {
        print_parameters.filename_ = iterParam->stringValue();
      }
    }

    else if (iterParam->tag() == "DELIMITER")
    {
      if (iterParam->stringValue() == "TAB") {
        print_parameters.delimiter_ = "\t";
      }

      else if (iterParam->stringValue() == "COMMA") {
        print_parameters.delimiter_ = ",";
      }

      else if (iterParam->stringValue() != "")
      {
        Report::UserWarning0() << "Invalid value of DELIMITER in .PRINT statment, ignoring";
      }
    }

    else
    {
      // This must be the first print variable.
      break;
    }
  }

  // Remaining parameters are variables
  print_parameters.variableList_.assign(iterParam, print_block.getParams().end());

  // Increase width to make sure that the columns do not run together
  if (print_parameters.streamWidth_ - print_parameters.streamPrecision_  < 9)
    print_parameters.streamWidth_ = print_parameters.streamPrecision_ + 9;

  // If the delimiter is specified on the command line it takes
  // precedence, reset the value of netListFilename_
  if (commandLine_.argExists("-delim"))
  {
    if (commandLine_.getArgumentValue("-delim") == "TAB")
      print_parameters.delimiter_ = "\t";
    else if (commandLine_.getArgumentValue("-delim") == "COMMA")
      print_parameters.delimiter_ = ",";
    else
      print_parameters.delimiter_ = commandLine_.getArgumentValue("-delim");
  }

  print_parameters.rawOverride_ = commandLine_.argExists("-r");

  // If the output file is specified on the command line it takes
  // precedence, reset the value of netListFilename_
  if (commandLine_.argExists("-o"))
  {
    print_parameters.filename_ = commandLine_.getArgumentValue("-o");
  }

  print_parameters.index_ = !noIndex_ 
                            && print_parameters.format_ != Format::PROBE 
                            && print_parameters.format_ != Format::TECPLOT
                            && print_parameters.format_ != Format::RAW 
                            && print_parameters.format_ != Format::RAW_ASCII;

  // Assemble the apropriate flavors of output variable lists based on the PRINT type and format
  if (PRINTType_ == PrintType::AC) {
    outputParameterMap_[OutputType::AC] = print_parameters;
    outputParameterMap_[OutputType::AC_IC] = print_parameters;
    outputParameterMap_[OutputType::HOMOTOPY] = print_parameters;
    outputParameterMap_[OutputType::SENS] = print_parameters;
    printParameters_ = &outputParameterMap_[OutputType::AC];
  }
  else if (PRINTType_ == PrintType::HB) {
    outputParameterMap_[OutputType::HB_FD] = print_parameters;
    outputParameterMap_[OutputType::HB_TD] = print_parameters;
    outputParameterMap_[OutputType::HB_IC] = print_parameters;
    outputParameterMap_[OutputType::HB_STARTUP] = print_parameters;
    printParameters_ = &outputParameterMap_[OutputType::HB_FD];
  }
  else if (PRINTType_ == PrintType::TRAN) {
    outputParameterMap_[OutputType::TRAN] = print_parameters;
    outputParameterMap_[OutputType::HOMOTOPY] = print_parameters;
    outputParameterMap_[OutputType::SENS] = print_parameters;
    printParameters_ = &outputParameterMap_[OutputType::TRAN];
  }
  else if (PRINTType_ == PrintType::DC) {
    outputParameterMap_[OutputType::DC] = print_parameters;
    outputParameterMap_[OutputType::HOMOTOPY] = print_parameters;
    outputParameterMap_[OutputType::SENS] = print_parameters;
    printParameters_ = &outputParameterMap_[OutputType::DC];
  }
  else {
    Report::UserError0() << "Unrecognized .PRINT type";
  }

  activeOutputterStack_.push_back(std::vector<Outputter::Interface *>());

#ifdef Xyce_DEBUG_IO
  Xyce::dout() << "netListFilename_ = " << netListFilename_ << std::endl;
  Xyce::dout() << " PRINTSET Registered:: finished" << std::endl;
#endif

  return true;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::registerIC
// Purpose       : registers set of variables to set for .IC or .DCVOLT
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 09/06/07
//-----------------------------------------------------------------------------
bool OutputMgr::registerIC(const Util::OptionBlock & OB)
{
#ifdef Xyce_DEBUG_IO
  Xyce::dout() << " IC Registered::" << std::endl;
#endif

  ICflag_ = true;

  // save the ICblock_.  This will be needed later.
  // We allow multiple .IC blocks in the netlist, so there needs to be an
  // STL vector of these.
  ICblockVec_.push_back(OB);
  int end = ICblockVec_.size()-1;

#ifdef Xyce_DEBUG_IO
  ParameterList::iterator iterParam = ICblockVec_[end].getParams().begin();
  ParameterList::iterator endParam  = ICblockVec_[end].getParams().end();

  while (iterParam != endParam)
  {
   Xyce::dout() << "iterParam->tag = " << iterParam->tag() << std::endl;
     ++iterParam;
  }
#endif

  return true;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::registerNodeSet
// Purpose       : registers set of variables to set for .NODESET
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 09/06/07
//-----------------------------------------------------------------------------
bool OutputMgr::registerNodeSet(const Util::OptionBlock & OB)
{
#ifdef Xyce_DEBUG_IO
 Xyce::dout() << " NODESET Registered::" << std::endl;
#endif

  nodesetflag_ = true;

  // save the nodesetblock_.  This will be needed later.
  // We allow multiple .nodeset blocks in the netlist, so there needs to be an
  // STL vector of these.
  nodesetblockVec_.push_back(OB);
  int end = nodesetblockVec_.size()-1;

#ifdef Xyce_DEBUG_IO
  ParameterList::iterator iterParam = nodesetblockVec_[end].getParams().begin();
  ParameterList::iterator endParam  = nodesetblockVec_[end].getParams().end();

  while (iterParam != endParam)
  {
   Xyce::dout() << "iterParam->tag = " << iterParam->tag() << std::endl;
     ++iterParam;
  }
#endif

  return true;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::registerSave
// Purpose       : registers set of variables to set for .SAVE.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/11/07
//-----------------------------------------------------------------------------
bool OutputMgr::registerSave(const Util::OptionBlock & OB)
{
#ifdef Xyce_DEBUG_IO
 Xyce::dout() << " SAVE Registered::" << std::endl;
#endif

  saveFlag_ = true;

  ExtendedString sval("");

  ParameterList::const_iterator iterPL = OB.getParams().begin();
  ParameterList::const_iterator iterPL_end = OB.getParams().end();

  while (iterPL != iterPL_end)
  {
#ifdef Xyce_DEBUG_IO
   Xyce::dout() << "iterPL->tag = " << iterPL->tag() << std::endl;
#endif
    if (iterPL->tag() == "TYPE")
    {
      sval = iterPL->stringValue();
      sval.toUpper();

      if (sval == "NODESET" || sval == ".NODESET")
      {
        saveFileType_ = ".NODESET";
      }
      else if (sval == "IC" || sval == ".IC")
      {
        saveFileType_ = ".IC";
      }
      else
      {
        Report::UserWarning0() << "Unrecognized type specified on .SAVE command.  Defaulting to .NODESET";
      }
    }
    else if (iterPL->tag() == "FILE")
    {
      saveOutputFile_ = iterPL->stringValue();
    }
    else if (iterPL->tag() == "TIME")
    {
      // do nothing, this will be handled in the time integrator.
    }
    else if (iterPL->tag() == "LEVEL")
    {
      sval = iterPL->stringValue();
      sval.toUpper();

      if (sval == "ALL")
      {
        // do nothing
      }
      else if (sval == "NONE")
      {
        // none means don't output anything.  Pretend the .SAVE line isn't in the netlist.
        saveFlag_ = false;
        saveFileLevel_ = "NONE";
      }
      else if (sval == "TOP")
      {
        Report::UserWarning0() << "LEVEL=TOP in .SAVE line not supported.  Ignoring. ";
      }
      else
      {
        Report::UserWarning0() << "Unrecognized LEVEL " << sval << " specified in .SAVE command";
      }
    }
    else
    {
      Report::UserWarning0() << "Parameter " << iterPL->tag() << " not recognized in .SAVE command";
    }


    ++iterPL;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::registerLoad
// Purpose       : registers set of variables to set for .SAVE.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/11/07
//-----------------------------------------------------------------------------
bool OutputMgr::registerLoad(const Util::OptionBlock & OB)
{
  loadFlag_ = true;

  Report::UserWarning0() << ".LOAD not supported yet.  Use .INCLUDE instead";

  return true;
}


//-----------------------------------------------------------------------------
// Function      : OutputMgr::registerOPOptions
// Purpose       : registers set of variables to set for .OP.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/15/07
//-----------------------------------------------------------------------------
bool OutputMgr::registerOP(const Util::OptionBlock & OB)
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::registerSens
// Purpose       : registers set of variables to set for .SENS.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 2/10/2014
//-----------------------------------------------------------------------------
bool OutputMgr::registerSens(const Util::OptionBlock & OB)
{
  bool bsuccess = true;
  std::list<N_UTL_Param>::const_iterator iter = OB.getParams().begin();
  std::list<N_UTL_Param>::const_iterator end   = OB.getParams().end();

  for ( ; iter != end; ++ iter)
  {
    if (iter->uTag() == "OBJFUNC")
    {
      ExtendedString func = iter->stringValue();
      func.toUpper();
      sensObjFunction_ = func;
      sensObjFuncGiven_ = true;
    }
    else if ( std::string( iter->uTag() ,0,5) == "PARAM") // this is a vector
    {
      ExtendedString tag = iter->stringValue();
      tag.toUpper();
      sensParamNameVec_.push_back(tag);
    }
    else
    {
      Xyce::Report::UserWarning() << iter->uTag() << " is not a recognized sensitivity solver option.\n" << std::endl;
    }
  }
  enableSensitivityFlag_ = true;


  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::enableHomotopyOutput
// Purpose       : turns on Homotopy output
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/10/00
//-----------------------------------------------------------------------------
void OutputMgr::enableHomotopyOutput(Analysis::Analysis_Mode analysis_mode)
{
  // prepare output manager to write
  if (homotopyFlag_)
  {
    Report::UserWarning0() << "Homotopyfile already initialized.  Contents may be overwritten.";
  }
  else
  {
    homotopyFlag_ = true;

    PrintParameters homotopy_print_parameters = outputParameterMap_[OutputType::HOMOTOPY];

    if (analysis_mode == Analysis::ANP_MODE_TRANSIENT)
      homotopy_print_parameters.variableList_.push_front(Util::Param("TIME", 0.0));
    if (homotopy_print_parameters.index_)
      homotopy_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));

    fixupPrintParameters(homotopy_print_parameters);

    Outputter::Interface *outputter;
    if (homotopy_print_parameters.format_ == Format::STD) {
      homotopy_print_parameters.defaultExtension_ = ".HOMOTOPY.prn";
      homotopy_print_parameters.extraExtension_ = ".HOMOTOPY.prn";
      outputter = new Outputter::HomotopyPrn(*this, homotopy_print_parameters);
    }
    else if (homotopy_print_parameters.format_ == Format::TECPLOT) {
      homotopy_print_parameters.defaultExtension_ = ".HOMOTOPY.dat";
      homotopy_print_parameters.extraExtension_ = ".HOMOTOPY.dat";
      outputter = new Outputter::HomotopyTecPlot(*this, homotopy_print_parameters);
    }
    // else if (homotopy_print_parameters.format_ == Format::CSV) {
    //   homotopy_print_parameters.defaultExtension_ = ".HOMOTOPY.csv";
    //   outputter                            = new Outputter::HomotopyTecplot(*this, homotopy_print_parameters);
    // }
    else if (homotopy_print_parameters.format_ == Format::PROBE) {
      homotopy_print_parameters.defaultExtension_ = ".HOMOTOPY.csd";
      homotopy_print_parameters.extraExtension_ = ".HOMOTOPY.csd";
      outputter = new Outputter::HomotopyProbe(*this, homotopy_print_parameters);
    }
    else
    {
      Report::UserWarning0() 
        << "Homotopy output cannot be written in " 
        << homotopy_print_parameters.format_ << " format, using standard format";

      outputter = new Outputter::HomotopyPrn(*this, homotopy_print_parameters);
    }

    outputter->parse();
    outputterMap_[PrintType::HOMOTOPY] = outputter;
    addActiveOutputter(PrintType::HOMOTOPY);
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::enableSensitivityOutput
// Purpose       : turns on sensitivity output
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 2/10/2014
//-----------------------------------------------------------------------------
void OutputMgr::enableSensitivityOutput (Analysis::Analysis_Mode analysis_mode)
{
  // prepare output manager to write
  if (sensitivityFlag_)
  {
    Report::UserWarning0() << "Sensfile already initialized.  Contents may be overwritten.";
  }
  else
  {
    sensitivityFlag_ = true;

    PrintParameters sensitivity_print_parameters = outputParameterMap_[OutputType::SENS];

    if (analysis_mode == Analysis::ANP_MODE_TRANSIENT)
      sensitivity_print_parameters.variableList_.push_front(Util::Param("TIME", 0.0));
    if (sensitivity_print_parameters.index_)
      sensitivity_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));

    fixupPrintParameters(sensitivity_print_parameters);

    Outputter::Interface *outputter;
    if (sensitivity_print_parameters.format_ == Format::STD) {
      sensitivity_print_parameters.defaultExtension_ = ".SENS.prn";
      sensitivity_print_parameters.extraExtension_ = ".SENS.prn";
      outputter = new Outputter::SensitivityPrn(*this, sensitivity_print_parameters);
    }
    else if (sensitivity_print_parameters.format_ == Format::TECPLOT) {
      sensitivity_print_parameters.defaultExtension_ = ".SENS.dat";
      sensitivity_print_parameters.extraExtension_ = ".SENS.dat";
      outputter = new Outputter::SensitivityTecPlot(*this, sensitivity_print_parameters);
    }
#if 0
    else if (sensitivity_print_parameters.format_ == Format::PROBE) {
      sensitivity_print_parameters.defaultExtension_ = ".SENS.csd";
      sensitivity_print_parameters.extraExtension_ = ".SENS.csd";
      outputter = new Outputter::SensitivityProbe(*this, sensitivity_print_parameters);
    }
#endif
    else
    {
      Report::UserWarning0() 
        << "Sensitivity output cannot be written in requested format, using standard format";

      sensitivity_print_parameters.defaultExtension_ = ".SENS.prn";
      sensitivity_print_parameters.extraExtension_ = ".SENS.prn";
      outputter = new Outputter::SensitivityPrn(*this, sensitivity_print_parameters);
    }

    outputter->parse();
    outputterMap_[PrintType::SENS] = outputter;
    addActiveOutputter(PrintType::SENS);
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::enableACoutput
// Purpose       : turns on AC output
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 5/22/2013
//-----------------------------------------------------------------------------
void OutputMgr::enableACOutput()
{
  // prepare output manager to write
  if (ACFlag_)
  {
    Report::UserWarning0() << "AC file already initialized.  Contents may be overwritten.";
  }

  else
  {
    ACFlag_ = true;

    PrintParameters freq_print_parameters = outputParameterMap_[OutputType::AC];
    PrintParameters ac_ic_print_parameters = outputParameterMap_[OutputType::AC_IC];

    if (freq_print_parameters.format_ != Format::PROBE)
    {
      freq_print_parameters.variableList_.push_front(Util::Param("FREQUENCY", 0.0));
    }
    if (freq_print_parameters.index_)
    {
      freq_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));
    }
    freq_print_parameters.expandComplexTypes_ = freq_print_parameters.format_ != Format::PROBE
                                                && freq_print_parameters.format_ != Format::RAW
                                                && freq_print_parameters.format_ != Format::RAW_ASCII;

    if (ac_ic_print_parameters.format_ != Format::PROBE)
    {
      ac_ic_print_parameters.variableList_.push_front(Util::Param("TIME", 0.0));
    }
    if (ac_ic_print_parameters.index_)
    {
      ac_ic_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));
    }

    fixupPrintParameters(freq_print_parameters);
    fixupPrintParameters(ac_ic_print_parameters);

    Outputter::Interface *outputter_fd;
    Outputter::Interface *outputter_prn;
    if (freq_print_parameters.format_ == Format::STD) {
      ac_ic_print_parameters.extraExtension_ = ".TD.prn";
      outputter_fd = new Outputter::FrequencyPrn(*this, freq_print_parameters);
      outputter_prn = new Outputter::TimePrn(*this, ac_ic_print_parameters);
    }
    else if (freq_print_parameters.format_ == Format::CSV) {
      ac_ic_print_parameters.extraExtension_ = ".TD.csv";
      outputter_fd = new Outputter::FrequencyCSV(*this, freq_print_parameters);
      outputter_prn = new Outputter::TimeCSV(*this, ac_ic_print_parameters);
    }
    else if (freq_print_parameters.format_ == Format::PROBE) {
      ac_ic_print_parameters.extraExtension_ = ".TD.csd";
      outputter_fd = new Outputter::FrequencyProbe(*this, freq_print_parameters);
      outputter_prn = new Outputter::TimeProbe(*this, ac_ic_print_parameters);
    }
    else if (freq_print_parameters.format_ == Format::TECPLOT) {
      ac_ic_print_parameters.extraExtension_ = ".TD.dat";
      outputter_fd = new Outputter::FrequencyTecPlot(*this, freq_print_parameters);
      outputter_prn = new Outputter::TimeTecPlot(*this, ac_ic_print_parameters);
    }
    else if (freq_print_parameters.format_ == Format::RAW) {
      ac_ic_print_parameters.extraExtension_ = ".TD.raw";
      if (commandLine_.argExists("-a")) {
        outputter_fd = new Outputter::FrequencyRawAscii(*this, freq_print_parameters);
        outputter_prn = new Outputter::TimeRawAscii(*this, ac_ic_print_parameters);
      }
      else {
        outputter_fd = new Outputter::FrequencyRaw(*this, freq_print_parameters);
        outputter_prn = new Outputter::TimeRaw(*this, ac_ic_print_parameters);
      }
    }
    else
    {
      Report::UserWarning0() << "AC output cannot be written in " << freq_print_parameters.format_ << " format, using standard format";

      outputter_fd = new Outputter::FrequencyPrn(*this, freq_print_parameters);
      outputter_prn = new Outputter::TimePrn(*this, ac_ic_print_parameters);
    }

    outputter_fd->parse();
    outputter_prn->parse();
    outputterMap_[PrintType::AC] = outputter_fd;
    outputterMap_[PrintType::AC_IC] = outputter_prn;
    addActiveOutputter(PrintType::AC);
    addActiveOutputter(PrintType::AC_IC);
  }
}


//-----------------------------------------------------------------------------
// Function      : OutputMgr::enableTransientOutput
// Purpose       : turns on Transient output
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 6/26/2013
//-----------------------------------------------------------------------------
void OutputMgr::enableTransientOutput()
{
  // prepare output manager to write
  if (tranFlag_)
  {
    Report::UserWarning0() << "Transient file already initialized.  Contents may be overwritten.";
  }

  else
  {
    tranFlag_ = true;

    PrintParameters transient_print_parameters = outputParameterMap_[OutputType::TRAN];

    if (transient_print_parameters.format_ != Format::PROBE)
      transient_print_parameters.variableList_.push_front(Util::Param("TIME", 0.0));
    if (transient_print_parameters.index_)
      transient_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));

    fixupPrintParameters(transient_print_parameters);

    Outputter::Interface *outputter;
    if (transient_print_parameters.format_ == Format::STD) {
      outputter = new Outputter::TimePrn(*this, transient_print_parameters);
    }
    else if (transient_print_parameters.format_ == Format::CSV) {
      outputter = new Outputter::TimeCSV(*this, transient_print_parameters);
    }
    else if (transient_print_parameters.format_ == Format::TECPLOT) {
      outputter = new Outputter::TimeTecPlot(*this, transient_print_parameters);
    }
    else if (transient_print_parameters.format_ == Format::PROBE) {
      outputter = new Outputter::TimeProbe(*this, transient_print_parameters);
    }
    else if (transient_print_parameters.format_ == Format::RAW) {
      if (commandLine_.argExists("-a"))
        outputter = new Outputter::TimeRawAscii(*this, transient_print_parameters);
      else
        outputter = new Outputter::TimeRaw(*this, transient_print_parameters);
    }
    else
    {
      Report::UserWarning0() << "Transient output cannot be written in " << transient_print_parameters.format_ << " format, using standard format";

      outputter = new Outputter::TimePrn(*this, transient_print_parameters);
    }

    outputter->parse();
    outputterMap_[PrintType::TRAN] = outputter;
    addActiveOutputter(PrintType::TRAN);
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::enableMPDEOutput
// Purpose       : turns on MPDE output
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 5/21/2013
//-----------------------------------------------------------------------------
void OutputMgr::enableMPDEOutput()
{
  // prepare output manager to write
  if (MPDEFlag_)
  {
    Report::UserWarning0() << "MPDE file already initialized.  Contents may be overwritten.";
  }

  else
  {
    MPDEFlag_ = true;

    PrintParameters mpde_print_parameters = outputParameterMap_[OutputType::TRAN];
    PrintParameters mpde_ic_print_parameters = outputParameterMap_[OutputType::TRAN];

    mpde_ic_print_parameters.defaultExtension_ = ".mpde_ic.prn";

    if (mpde_ic_print_parameters.format_ != Format::PROBE)
      mpde_ic_print_parameters.variableList_.push_front(Util::Param("TIME", 0.0));
    if (mpde_ic_print_parameters.index_)
      mpde_ic_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));

    fixupPrintParameters(mpde_print_parameters);
    fixupPrintParameters(mpde_ic_print_parameters);

    Outputter::Interface *outputter_mpde;
    Outputter::Interface *outputter_mpde_ic;
    if (mpde_print_parameters.format_ == Format::STD) {
      outputter_mpde = new Outputter::MPDEPrn(*this, mpde_print_parameters);
      outputter_mpde_ic = new Outputter::TimePrn(*this, mpde_ic_print_parameters);
    }
    // else if (format_ == Format::CSV) {
    //   outputter_mpde = new Outputter::CSV(*this, mpde_print_parameters);
    // }
    else if (mpde_print_parameters.format_ == Format::TECPLOT) {
      outputter_mpde = new Outputter::MPDETecPlot(*this, mpde_print_parameters);
      outputter_mpde_ic = new Outputter::TimeTecPlot(*this, mpde_ic_print_parameters);
    }
    // else if (format_ == Format::PROBE) {
    //   outputter_mpde = new Outputter::TimeProbe(*this, mpde_print_parameters);
    // }
    else
    {
      Report::UserWarning0() << "MPDE output cannot be written in " << mpde_print_parameters.format_ << " format, using standard format";

      outputter_mpde = new Outputter::MPDEPrn(*this, mpde_print_parameters);
      outputter_mpde_ic = new Outputter::TimePrn(*this, mpde_ic_print_parameters);
    }

    outputter_mpde->parse();
    outputter_mpde_ic->parse();
    outputterMap_[PrintType::MPDE] = outputter_mpde;
    outputterMap_[PrintType::MPDE_IC] = outputter_mpde_ic;
    addActiveOutputter(PrintType::MPDE);
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::enableDCSweepOutput
// Purpose       : Enable output for DC sweeps
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 6/26/2013
//-----------------------------------------------------------------------------
void OutputMgr::enableDCSweepOutput()
{
  // prepare output manager to write
  if (DCSweepFlag_)
  {
    Report::UserWarning0() << "DC file already initialized.  Contents may be overwritten.";
  }

  else
  {
    DCSweepFlag_ = true;

    PrintParameters dc_print_parameters = outputParameterMap_[OutputType::DC];

    if (dc_print_parameters.index_)
      dc_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));

    fixupPrintParameters(dc_print_parameters);

    Outputter::Interface *outputter_prn;
    if (dc_print_parameters.format_ == Format::STD) {
      dc_print_parameters.defaultExtension_ = ".prn";
      outputter_prn = new Outputter::TimePrn(*this, dc_print_parameters);
    }
    else if (dc_print_parameters.format_ == Format::CSV) {
      outputter_prn = new Outputter::TimeCSV(*this, dc_print_parameters);
    }
    else if (dc_print_parameters.format_ == Format::RAW) {
      outputter_prn = new Outputter::TimeRaw(*this, dc_print_parameters);
    }
    else if (dc_print_parameters.format_ == Format::TECPLOT) {
      outputter_prn = new Outputter::TimeTecPlot(*this, dc_print_parameters);
    }
    else if (dc_print_parameters.format_ == Format::PROBE) {
      outputter_prn = new Outputter::TimeProbe(*this, dc_print_parameters);
    }
    else
    {
      Report::UserWarning0() << "DC output cannot be written in " << dc_print_parameters.format_ << " format, using standard format";

      outputter_prn = new Outputter::TimePrn(*this, dc_print_parameters);
    }

    outputter_prn->parse();
    outputterMap_[PrintType::TRAN] = outputter_prn;
    addActiveOutputter(PrintType::TRAN);
  }
}


//-----------------------------------------------------------------------------
// Function      : OutputMgr::enableHBoutput
// Purpose       : turns on HB output
// Special Notes :
// Scope         : public
// Creator       : David Baur
// Creation Date : 5/21/2013
//-----------------------------------------------------------------------------
void OutputMgr::enableHBOutput()
{
  /*
  if (getNumProcs() > 1) {
    Report::UserFatal0() << "Parallel processing is not supported in HB analysis";
  }
  */

  // prepare output manager to write
  if (HBFlag_)
  {
    Report::UserWarning0() << "HBfile already initialized.  Contents may be overwritten.";
  }

  else
  {
    HBFlag_ = true;

    PrintParameters freq_print_parameters = outputParameterMap_[OutputType::HB_FD];
    if (freq_print_parameters.format_ != Format::PROBE)
      freq_print_parameters.variableList_.push_front(Util::Param("FREQUENCY", 0.0));
    if (freq_print_parameters.index_)
      freq_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));
    freq_print_parameters.expandComplexTypes_ = freq_print_parameters.format_ != Format::PROBE
                                                && freq_print_parameters.format_ != Format::RAW
                                                && freq_print_parameters.format_ != Format::RAW_ASCII;

    PrintParameters time_print_parameters = outputParameterMap_[OutputType::HB_TD];
    if (time_print_parameters.format_ != Format::PROBE)
      time_print_parameters.variableList_.push_front(Util::Param("TIME", 0.0));
    if (time_print_parameters.index_)
      time_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));

    PrintParameters hb_ic_print_parameters = outputParameterMap_[OutputType::HB_IC];
    if (hb_ic_print_parameters.format_ != Format::PROBE)
      hb_ic_print_parameters.variableList_.push_front(Util::Param("TIME", 0.0));
    if (hb_ic_print_parameters.index_)
      hb_ic_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));

    PrintParameters hb_startup_print_parameters = outputParameterMap_[OutputType::HB_STARTUP];
    if (hb_startup_print_parameters.format_ != Format::PROBE)
      hb_startup_print_parameters.variableList_.push_front(Util::Param("TIME", 0.0));
    if (hb_startup_print_parameters.index_)
      hb_startup_print_parameters.variableList_.push_front(Util::Param("INDEX", 0.0));

    fixupPrintParameters(freq_print_parameters);
    fixupPrintParameters(time_print_parameters);
    fixupPrintParameters(hb_ic_print_parameters);
    fixupPrintParameters(hb_startup_print_parameters);

    Outputter::Interface *outputter_hb;
    Outputter::Interface *outputter_init;
    Outputter::Interface *outputter_startup;
    if (freq_print_parameters.format_ == Format::STD) {
      hb_ic_print_parameters.defaultExtension_ = ".hb_ic.prn";
      hb_startup_print_parameters.defaultExtension_ = ".startup.prn";
      outputter_hb = new Outputter::HBPrn(*this, freq_print_parameters, time_print_parameters);
      outputter_init = new Outputter::TimePrn(*this, hb_ic_print_parameters);
      outputter_startup = new Outputter::TimePrn(*this, hb_startup_print_parameters);
    }
    else if (freq_print_parameters.format_ == Format::CSV) {
      hb_ic_print_parameters.defaultExtension_ = ".hb_ic.csv";
      hb_startup_print_parameters.defaultExtension_ = ".startup.csv";
      outputter_hb = new Outputter::HBCSV(*this, freq_print_parameters, time_print_parameters);
      outputter_init = new Outputter::TimeCSV(*this, hb_ic_print_parameters);
      outputter_startup = new Outputter::TimeCSV(*this, hb_startup_print_parameters);
    }
    else if (freq_print_parameters.format_ == Format::TECPLOT) {
      hb_ic_print_parameters.defaultExtension_ = ".hb_ic.dat";
      hb_startup_print_parameters.defaultExtension_ = ".startup.dat";
      outputter_hb = new Outputter::HBTecPlot(*this, freq_print_parameters, time_print_parameters);
      outputter_init = new Outputter::TimeTecPlot(*this, hb_ic_print_parameters);
      outputter_startup = new Outputter::TimeTecPlot(*this, hb_startup_print_parameters);
    }
    else
    {
      Report::UserWarning0() << "HB output cannot be written in " << freq_print_parameters.format_ << " format, using standard format";

      outputter_hb = new Outputter::HBPrn(*this, freq_print_parameters, time_print_parameters);
      outputter_init = new Outputter::TimePrn(*this, hb_ic_print_parameters);
    }

    outputter_hb->parse();
    outputter_init->parse();
    outputter_startup->parse();
    outputterMap_[PrintType::HB] = outputter_hb;
    outputterMap_[PrintType::HB_IC] = outputter_init;
    outputterMap_[PrintType::HB_STARTUP] = outputter_startup;
    addActiveOutputter(PrintType::HB);
    addActiveOutputter(PrintType::HB_IC);
  }
}


//-----------------------------------------------------------------------------
// Function      : OutputMgr::enableOverrideRawOutput
// Purpose       : turns on RAW output from the "-r" command line override
// Special Notes : "Override" raw is different from "format=raw" on the
//                 .print line.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 7/2/2013
//-----------------------------------------------------------------------------
void OutputMgr::enableOverrideRawOutput(const PrintParameters &print_parameters)
{
  // prepare output manager to write
  if (rawFlag_)
  {
    Report::UserWarning0() << "Rawfile already initialized.  Contents may be overwritten.";
  }

  else
  {
    rawFlag_ = true;

    PrintParameters raw_print_parameters = print_parameters;
    raw_print_parameters.filename_ = commandLine_.getArgumentValue("-r");

    fixupPrintParameters(raw_print_parameters);

    Outputter::Interface *outputter;
    if (commandLine_.argExists("-a"))
      outputter = new Outputter::OverrideRawAscii(*this, raw_print_parameters);
    else
      outputter = new Outputter::OverrideRaw(*this, raw_print_parameters);

    outputter->parse();
    outputterMap_[PrintType::RAW_OVERRIDE] = outputter;
    addActiveOutputter(PrintType::RAW_OVERRIDE);
  }
}


//-----------------------------------------------------------------------------
// Function      : OutputMgr::enableHDF5Output
// Purpose       : creates HDF5 output object and prepares for writing
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 03/28/2012
//-----------------------------------------------------------------------------
bool OutputMgr::enableHDF5Output()
{
  bool result = true;

#ifdef Xyce_USE_HDF5
  hdf5PlistId_ = H5Pcreate(H5P_FILE_ACCESS);

  // based on if we're running in serial or parallel, add parallel output
  // option to hdf5PlistId_
  if (!pdsCommPtr_->isSerial())
  {
#ifdef Xyce_PARALLEL_MPI
    N_PDS_ParComm *aParComm = dynamic_cast<N_PDS_ParComm *>(pdsCommPtr_);
    H5Pset_fapl_mpio(hdf5PlistId_, aParComm->mpiComm()->comm(), MPI_INFO_NULL);
    H5Pset_fclose_degree(hdf5PlistId_, H5F_CLOSE_STRONG);
#endif
  }
  std::string hdf5ExtendedFileName = hdf5FileName_ + filenameSuffix_ + ".h5d";
  hdf5FileId_ = H5Fcreate(hdf5ExtendedFileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, hdf5PlistId_);
#endif  // Xyce_USE_HDF5

  return result;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::updateHDF5Output
// Purpose       : writes curent solution data to the HDF5 output object
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 03/28/2012
//-----------------------------------------------------------------------------
bool OutputMgr::updateHDF5Output(
      N_LAS_Vector * solnVecPtr)
{
  bool result = true;

#ifdef Xyce_USE_HDF5
  int status=0;
  if (!hdf5HeaderWritten_)
  {
    hdf5HeaderWritten_=true;
    // write out solution var map

    // now get all the node names and GID's.  We need to create the dataset on all
    // processors and then write it distributed over all processors.
    //std::map<int, std::string> nodeRef;
    //std::vector<char> typeRef;
    //topology_->getRawData(nodeRef, typeRef);

    // in parallel we can't write an array of variale length node names.
    // Thus we need to make an array of fixed length strings the longest
    // of which is the maximum length of a node name.  This wastes a bit of
    // output space, but we only do this once in the solution var map.

    // first on each processor, find the max length of a node name.
    int maxNodeNameLength=0;
    NodeNamePairMap::iterator nameIter = allNodes_.begin();
    NodeNamePairMap::iterator nameEnd = allNodes_.end();

    for (; nameIter != nameEnd ; ++nameIter)
    {
      int nodeNameLength = nameIter->first.size();
      maxNodeNameLength = Xycemax(maxNodeNameLength, nodeNameLength);
    }

    // now have all the processors send the max length to proc 0.
    int globalMaxNodeNameLength;
    pdsCommPtr_->barrier();
    pdsCommPtr_->maxAll(&maxNodeNameLength, &globalMaxNodeNameLength, 1);

    globalMaxNodeNameLength++;  // add one for string terminator
    // now make a char[][] array of size [numLocalNodes][maxNodeNameLength]
    int numLocalNodes = solnVecPtr->localLength();
    char* nodeNameSpace = new char[ numLocalNodes * globalMaxNodeNameLength ];
    // could simplify this with a placement new operator(i.e. new nodeNameSpace nodeNameArray).
    char** nodeNameArray = new char * [ numLocalNodes ];
    for (int i=0; i< numLocalNodes; i++)
    {
      nodeNameArray[i] = nodeNameSpace + i * sizeof(char) * globalMaxNodeNameLength;//new char(globalMaxNodeNameLength);
    }
    // also need an array for global ID's
    int* nodeGIDArray = new int[ numLocalNodes ];

    // to fill the arrays in parallel, each processor needs to know how many
    // unknowns are on prior processors.  We could infer this from the GID values,
    // however that would not work in block analysis modes like MPDE, HB and AC
    // Thus, we will have the individual processes communicate this info .
    // This duplicates what could be done with an epetra multivector, if we could
    // construct one that held strings or char*.  We can't so for now this is what
    // we'll do.
    int* unknownsPerProc = new int[ getNumProcs() ];

    pdsCommPtr_->petraComm()->GatherAll(&numLocalNodes, unknownsPerProc, 1);
    std::vector<int> totalUnknownsPriorToProc;
    totalUnknownsPriorToProc.resize(getNumProcs());
    int sum = 0;
    for (int i=0; i<getNumProcs(); i++)
    {
      totalUnknownsPriorToProc[i]=sum;
      sum += unknownsPerProc[i];
    }

    // fill up the arrays
    nameIter = allNodes_.begin();
    for (int index=0; nameIter != nameEnd ; ++nameIter, index++)
    {
      strncpy(nodeNameArray[index], nameIter->first.c_str(), globalMaxNodeNameLength);
      nodeGIDArray[index] = solnVecPtr->pmap()->petraMap()->GID(nameIter->second.first);
    }

    // make the group for writing
    // in parallel, all processors must create groups and datasets.
    hid_t solVarMapGroup = H5Gcreate(hdf5FileId_, "SolVarMap", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (solVarMapGroup < 0)
      Xyce::dout() << "Error in making solVarMapGroup on procID_ " << getProcID() << std::endl;

    // make the file data space
    hsize_t nodeNamesDim[1] = {solnVecPtr->globalLength()};
    hid_t fileDataspace = H5Screate_simple(1,  nodeNamesDim, NULL);

    // make the memorySpace
    hsize_t nodeNamesLocalDim[1] = {solnVecPtr->localLength()};
    hid_t memorySpace = H5Screate_simple(1, nodeNamesLocalDim, NULL);

    // make the dataset
    hid_t fixedLenString = H5Tcreate(H5T_STRING, globalMaxNodeNameLength);
    hid_t nodeNameDataSet = H5Dcreate(solVarMapGroup, "VarNames", fixedLenString, fileDataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // select portion for writing
    hid_t fileDataspaceSelection = H5Dget_space(nodeNameDataSet);
    hsize_t offset[1] = {totalUnknownsPriorToProc[getProcID()]};
    hsize_t stride[1] = {1};
    hsize_t count[1] = {unknownsPerProc[getProcID()]};
    hsize_t block[1] = {1};
    H5Sselect_hyperslab(fileDataspaceSelection, H5S_SELECT_SET, offset, stride, count, block);

    // writing property list
    hid_t writePropertyList = H5Pcreate(H5P_DATASET_XFER);
    #ifdef Xyce_PARALLEL_MPI
    if (getNumProcs() > 1)
    {
      H5Pset_dxpl_mpio(writePropertyList, H5FD_MPIO_COLLECTIVE);
    }
    #endif

    status = H5Dwrite(nodeNameDataSet, fixedLenString, memorySpace, fileDataspaceSelection, writePropertyList, nodeNameSpace);

    // write the GID array
    hid_t gidDataSet = H5Dcreate(solVarMapGroup, "VarGID", H5T_NATIVE_INT, fileDataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sselect_hyperslab(fileDataspaceSelection, H5S_SELECT_SET, offset, stride, count, block);
    status = H5Dwrite(gidDataSet, H5T_NATIVE_INT, memorySpace, fileDataspaceSelection, writePropertyList, nodeGIDArray);
    status = H5Dclose(gidDataSet);

    delete [] nodeNameArray;
    delete [] nodeGIDArray;
    delete [] unknownsPerProc;

    status = H5Tclose(fixedLenString);
    if (status < 0)
      Xyce::dout() << "Error closing fixedLenString. on " << getProcID() << std::endl;
    status = H5Dclose(nodeNameDataSet);
    if (status < 0)
      Xyce::dout() << "Error closing nodeNameDataSet. on " << getProcID() << std::endl;
    status = H5Sclose(fileDataspace);
    if (status < 0)
      Xyce::dout() << "Error closing solVarMapDataspace. on " << getProcID() << std::endl;
    status = H5Gclose(solVarMapGroup);
    if (status < 0)
      Xyce::dout() << "Error closing solVarMapGroup. on " << getProcID() << std::endl;

    // now output the independent variables
    hid_t independentVarGroup = H5Gcreate(hdf5FileId_, "IndepVars", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (independentVarGroup < 0)
      Xyce::dout() << "Error in making independentVarGroup on getProcID() " << getProcID() << std::endl;

    // filespace & memspace
    hsize_t independentDataDim[1] = {1};
    hsize_t independentMaxDataDim[1] = {H5S_UNLIMITED};
    hid_t independentFileSpace = H5Screate_simple(1, independentDataDim, independentMaxDataDim);
    hid_t independentMemorySpace = H5Screate_simple(1, independentDataDim, independentMaxDataDim);

    // array will be unlimited in length because we don't know how many steps Xyce will take
    // So we've set the max dimension to H5S_UNLIMITED above. we need to set the size of the
    // chunk by which it will grow as well
    // now the dataset
    hid_t independentDataProperty_ = H5Pcreate(H5P_DATASET_CREATE);
    hsize_t chunkSize[1] = {1};
    H5Pset_chunk(independentDataProperty_, 1, chunkSize);

    // make the dataset
    hid_t independentVarDataSet = H5Dcreate(independentVarGroup, "Time", H5T_NATIVE_DOUBLE, independentFileSpace, H5P_DEFAULT, independentDataProperty_, H5P_DEFAULT);

    // write
    writePropertyList=H5Pcreate(H5P_DATASET_XFER);
    #ifdef Xyce_PARALLEL_MPI
    if (getNumProcs() > 1)
    {
      H5Pset_dxpl_mpio(writePropertyList, H5FD_MPIO_INDEPENDENT);
    }
    #endif
    if (getProcID() == 0)
    {
      H5Dwrite(independentVarDataSet, H5T_NATIVE_DOUBLE, independentMemorySpace, independentFileSpace, writePropertyList, &circuitTime_);
    }
    H5Sclose(independentFileSpace);
    H5Sclose(independentMemorySpace);
    H5Dclose(independentVarDataSet);
    H5Gclose(independentVarGroup);
    H5Pclose(independentDataProperty_);

    // now output the first solution vector
    hid_t solutionVecGroup = H5Gcreate(hdf5FileId_, "SolutionVec", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (solutionVecGroup < 0)
      Xyce::dout() << "Error in making solutionVecGroup on getProcID() " << getProcID() << std::endl;

    // filespace & memspace
    hsize_t dependentDataLocalDim[2] = {1, solnVecPtr->localLength()};
    hsize_t dependentDataGlobalDim[2] = {1, solnVecPtr->globalLength()};
    hsize_t dependentMaxDataDim[2] = {H5S_UNLIMITED, solnVecPtr->globalLength()};
    // remember, filespace is the size of the global data on the disk.
    hid_t dependentFileSpace = H5Screate_simple(2, dependentDataGlobalDim, dependentMaxDataDim);
    // memspace defines the space used locally on a given processor.
    hid_t dependentMemorySpace = H5Screate_simple(2, dependentDataLocalDim, dependentMaxDataDim);

    // array will be unlimited in length because we don't know how many steps Xyce will take
    // So we've set the max dimension to H5S_UNLIMITED above. we need to set the size of the
    // chunk by which it will grow as well
    // now the dataset
    hid_t dependentDataProperty_ = H5Pcreate(H5P_DATASET_CREATE);
    hsize_t solChunkSize[2] = {1, solnVecPtr->globalLength()};
    H5Pset_chunk(dependentDataProperty_, 2, solChunkSize);

    // make the dataset
    hid_t dependentVarDataSet = H5Dcreate(solutionVecGroup, "Solution", H5T_NATIVE_DOUBLE, dependentFileSpace, H5P_DEFAULT, dependentDataProperty_, H5P_DEFAULT);

    // set up hyperslab to define the relationship between memspace and filespace.
    hid_t dependentVarFilespace = H5Dget_space(dependentVarDataSet);
    hsize_t solVecStart[2] = {0, solnVecPtr->pmap()->petraMap()->GID(0) - solnVecPtr->pmap()->petraMap()->IndexBase()};
    //hsize_t solVecStart[2] = {0, solnVecPtr->pmap()->petraMap()->IndexBase()};
    hsize_t solVecStride[2] = {1, 1};
    hsize_t solVecCount[2] = {1, solnVecPtr->localLength()};
    hsize_t solVecBlock[2] = {1, 1};
    H5Sselect_hyperslab(dependentFileSpace, H5S_SELECT_SET, solVecStart, solVecStride, solVecCount, solVecBlock);

    writePropertyList=H5Pcreate(H5P_DATASET_XFER);
    #ifdef Xyce_PARALLEL_MPI
    if (getNumProcs() > 1)
    {
      H5Pset_dxpl_mpio(writePropertyList, H5FD_MPIO_COLLECTIVE);
    }
    #endif

    // write the data
    H5Dwrite(dependentVarDataSet, H5T_NATIVE_DOUBLE, dependentMemorySpace, dependentFileSpace, writePropertyList, &((*solnVecPtr)[0]));

    H5Pclose(writePropertyList);
    H5Sclose(dependentFileSpace);
    H5Sclose(dependentMemorySpace);
    H5Dclose(dependentVarDataSet);
    H5Gclose(solutionVecGroup);
  }
  else
  {
    // update  independent variables
    hid_t independentVarGroup = H5Gopen(hdf5FileId_, "IndepVars", H5P_DEFAULT);
    hid_t independentVarDataSet = H5Dopen(independentVarGroup, "Time", H5P_DEFAULT);

    // data sets have been written to once.  Now extend them with new data
    hsize_t newDim[1] = {hdf5IndexValue_+1};
    hsize_t coords[1] = {hdf5IndexValue_};

    hid_t writePropertyList=H5Pcreate(H5P_DATASET_XFER);
    #ifdef Xyce_PARALLEL_MPI
    if (getNumProcs() > 1)
    {
      H5Pset_dxpl_mpio(writePropertyList, H5FD_MPIO_INDEPENDENT);
    }
    #endif

    H5Dset_extent(independentVarDataSet, newDim);

    // make the memspace
    hsize_t independentDataDim[1] = {1};
    hsize_t independentMaxDataDim[1] = {H5S_UNLIMITED};
    hid_t independentVarMemSpace = H5Screate_simple(1, independentDataDim, independentMaxDataDim);
    //hid_t independentVarMemSpace = H5Dget_space(independentVarDataSet);
    // get filespace
    hid_t independentVarFilespace = H5Dget_space(independentVarDataSet);

    // indicate that we have only updated the last element in this space
    //H5Sselect_elements(independentVarFilespace, H5S_SELECT_SET, 1, coords);
    hsize_t start[1] = {hdf5IndexValue_};
    hsize_t stride[1] = {1};
    hsize_t count[1] = {1};
    hsize_t block[1] = {1};
    H5Sselect_hyperslab(independentVarFilespace, H5S_SELECT_SET, start, stride, count, block);
    // write the update to the data space

    if (getProcID() == 0)
    {
      //H5Dwrite(independentVarDataSet, H5T_NATIVE_DOUBLE, independentVarMemSpace, independentVarFilespace, writePropertyList, &circuitTime_-hdf5IndexValue_);
      H5Dwrite(independentVarDataSet, H5T_NATIVE_DOUBLE, independentVarMemSpace, independentVarFilespace, writePropertyList, &circuitTime_);
    }

    H5Dclose(independentVarDataSet);
    H5Gclose(independentVarGroup);
    H5Sclose(independentVarMemSpace);
    H5Sclose(independentVarFilespace);

    // update solution
    hid_t solutionVecGroup = H5Gopen(hdf5FileId_, "SolutionVec", H5P_DEFAULT);
    hid_t dependentVarDataSet = H5Dopen(solutionVecGroup, "Solution", H5P_DEFAULT);

    // data sets have been written to once.  Now extend them with new data
    hsize_t solutionNewDim[2] = {hdf5IndexValue_+1, solnVecPtr->globalLength() };
    hsize_t solutionCoords[2] = {hdf5IndexValue_, 0};
    H5Dset_extent(dependentVarDataSet, solutionNewDim);

    // filespace & memspace
    hsize_t dependentDataLocalDim[2] = {1, solnVecPtr->localLength()};
    hsize_t dependentDataGlobalDim[2] = {1, solnVecPtr->globalLength()};
    hsize_t dependentMaxDataDim[2] = {H5S_UNLIMITED, solnVecPtr->globalLength()};
    // remember, filespace is the size of the global data on the disk.
    // hid_t dependentFileSpace = H5Screate_simple(2, dependentDataGlobalDim, dependentMaxDataDim);
    // memspace defines the space used locally on a given processor.
    hid_t dependentMemorySpace = H5Screate_simple(2, dependentDataLocalDim, dependentMaxDataDim);
    // get filespace
    hid_t dependentFileSpace = H5Dget_space(dependentVarDataSet);

    // set up hyperslab to define the relationship between memspace and filespace.
    hsize_t solVecStart[2] = {hdf5IndexValue_, solnVecPtr->pmap()->petraMap()->GID(0) - solnVecPtr->pmap()->petraMap()->IndexBase()};
    //hsize_t solVecStart[2] = {0, solnVecPtr->pmap()->petraMap()->IndexBase()};
    hsize_t solVecStride[2] = {1, 1};
    hsize_t solVecCount[2] = {1, solnVecPtr->localLength()};
    hsize_t solVecBlock[2] = {1, 1};
    H5Sselect_hyperslab(dependentFileSpace, H5S_SELECT_SET, solVecStart, solVecStride, solVecCount, solVecBlock);

    writePropertyList=H5Pcreate(H5P_DATASET_XFER);
    #ifdef Xyce_PARALLEL_MPI
    if (getNumProcs() > 1)
    {
      H5Pset_dxpl_mpio(writePropertyList, H5FD_MPIO_COLLECTIVE);
    }
    #endif

    // write the data
    H5Dwrite(dependentVarDataSet, H5T_NATIVE_DOUBLE, dependentMemorySpace, dependentFileSpace, writePropertyList, &((*solnVecPtr)[0]));

    H5Pclose(writePropertyList);
    H5Sclose(dependentFileSpace);
    H5Sclose(dependentMemorySpace);
    H5Dclose(dependentVarDataSet);
    H5Gclose(solutionVecGroup);
  }

  // update index before return
  hdf5IndexValue_++;
#endif  // Xyce_USE_HDF5

  return result;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::closeHDF5Output
// Purpose       : closes HDF5 ouput object
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 03/28/2012
//-----------------------------------------------------------------------------
bool OutputMgr::closeHDF5Output()
{
  bool result = true;

#ifdef Xyce_USE_HDF5
  int status = 0;
  H5Pclose(hdf5PlistId_);
  status = H5Fclose(hdf5FileId_);
  if (status < 0)
      Xyce::dout() << "Error closing hdf5FileId_. on " << getProcID() << std::endl;
#endif  // Xyce_USE_HDF5

  return result;
}


//-----------------------------------------------------------------------------
// Function      : Xyce::IO::removeStarVariables
// Purpose       : Process V(*) and I(*) on .print line
// Special Notes : Replaces v(*) and i(*) that appears on the .print line
//                 with a list of all v() or i() variables as appropriate.
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date : 6/28/2013
//-----------------------------------------------------------------------------
void
removeStarVariables(ParameterList &variable_list, N_PDS_Comm &communicator, const NodeNamePairMap &all_nodes, const NodeNamePairMap &external_nodes)
{
  bool vStarFound = false;
  bool iStarFound = false;

  ParameterList::iterator vStarPosition = variable_list.end();
  ParameterList::iterator iStarPosition = variable_list.end();

  // remove the v(*) and i(*) from the .print line
  for (ParameterList::iterator it = variable_list.begin(); it != variable_list.end(); )
  {
    // process the * entries
    if ((*it).tag() == "*")
    {
      // move to the type
      --it;

      // remember type of replacement and location of *
      if ((*it).tag() == "V")
      {
        vStarFound = true;
        vStarPosition = it;
        --vStarPosition;
      }
      else if ((*it).tag() == "I")
      {
        iStarFound = true;
        iStarPosition = it;
        --iStarPosition;
      }

      // remove the v or i
      it = variable_list.erase(it);

      // remove the *
      it = variable_list.erase(it);
    }

    // move to next item on .print line
    else
    {
      ++it;
    }
  }

  // create list of * entries
  if (vStarFound || iStarFound)
  {
    int sB = 0;

    ParameterList vStarList;

    for (NodeNamePairMap::const_iterator it = external_nodes.begin(); it != external_nodes.end() ; ++it)
    {
      ExtendedString tmpStr((*it).first);
      tmpStr.toUpper();
// DEBUG_ELR:  do all curr end in "branch"?
      std::string::size_type pos = tmpStr.rfind("BRANCH");

      if (pos == std::string::npos && vStarFound)
      {
        vStarList.push_back(Util::Param("V", 1.0));
        vStarList.push_back(Util::Param(tmpStr, 0.0));

        if (!communicator.isSerial())
        {
          sB += tmpStr.size();
          sB += sizeof(int);
        }
      }
    }

    ParameterList iStarList;

    for (NodeNamePairMap::const_iterator iter_a = all_nodes.begin(); iter_a != all_nodes.end() ; ++iter_a)
    {
      ExtendedString tmpStr((*iter_a).first);
      tmpStr.toUpper();
      size_t pos = tmpStr.rfind("BRANCH");

// DEBUG_ELR: drop Y*branch; broken
//        else if (pos != std::string::npos && iStarFound)
      if (pos != std::string::npos && iStarFound && tmpStr[0] != 'Y')
      {
        tmpStr = tmpStr.substr(0, pos - 1);

        iStarList.push_back(Util::Param("I", 1.0));
        iStarList.push_back(Util::Param(tmpStr, 0.0));

        if (!communicator.isSerial())
        {
          sB += tmpStr.size();
          sB += sizeof(int);
        }
      }
    }

    if (!communicator.isSerial())
    {
      // setup buffers for exchange
      int localSize = sizeof(int) + sizeof(int) + sB;
      int bSize, p, size, length;

      communicator.maxAll(&localSize, &bSize,  1);

      // reserve memory for exchange
      char * tmpBuffer = new char[ bSize ];

      // sync the temporary list to the common print block
      if (communicator.procID() == 0)
      {
        sB = 0;

        // DEBUG_ELR: switch to gatherAll()
        // rx local lists from pN
        for (int i = 1; i < communicator.numProc() ; ++i)
        {
          p = 0;

          communicator.recv(tmpBuffer, bSize, i);

          communicator.unpack(tmpBuffer, bSize, p, &size, 1);

          for (int j = 0; j < size; ++j)
          {
            communicator.unpack(tmpBuffer, bSize, p, &length, 1);
            std::string tmpStr(tmpBuffer + p, length);
            p += length;

            Util::Param param;
            param.setTag("V");
            param.setVal(1.0);
            vStarList.push_back(param);

            param.setTag(tmpStr);
            param.setVal(0.0);
            vStarList.push_back(param);

            sB += length + sizeof(int);
          }

          communicator.unpack(tmpBuffer, bSize, p, &size, 1);

          for (int j = 0; j < size; ++j)
          {
            communicator.unpack(tmpBuffer, bSize, p, &length, 1);
            std::string tmpStr(tmpBuffer + p, length);
            p += length;

            Util::Param param;
            param.setTag("I");
            param.setVal(1.0);
            iStarList.push_back(param);

            param.setTag(tmpStr);
            param.setVal(0.0);
            iStarList.push_back(param);

            sB += length + sizeof(int);
          }
        }

        // free memory
        delete [] tmpBuffer;

        // allocate mem
        bSize += sB;
        tmpBuffer = new char[ bSize ];

        // pack new list
        p = 0;
        size = vStarList.size() / 2;
        ParameterList::iterator iter_p = vStarList.begin();
        ParameterList::iterator iter_p_end = vStarList.end();

        communicator.pack(&size, 1, tmpBuffer, bSize, p);

        // pack name sizes and chars
        for (; iter_p != iter_p_end; ++iter_p)
        {
          // move to id
          ++iter_p;
          size =(iter_p->tag()).length();
          communicator.pack(&size, 1, tmpBuffer, bSize, p);
          communicator.pack((iter_p->tag()).c_str(), size, tmpBuffer, bSize, p);
        }

        size = iStarList.size() / 2;
        iter_p = iStarList.begin();
        iter_p_end = iStarList.end();

        communicator.pack(&size, 1, tmpBuffer, bSize, p);

        // pack name sizes and chars
        for (; iter_p != iter_p_end; ++iter_p)
        {
          // move to id
          ++iter_p;
          size =(iter_p->tag()).length();
          communicator.pack(&size, 1, tmpBuffer, bSize, p);
          communicator.pack((iter_p->tag()).c_str(), size, tmpBuffer, bSize, p);
        }

        // tx merged list to pN
        communicator.bcast(&bSize, 1, 0);
        communicator.bcast(tmpBuffer, bSize, 0);

        // free memory
        delete [] tmpBuffer;
      }
      else
      {
        p = 0;
        size = vStarList.size() / 2;
        ParameterList::iterator iter_p = vStarList.begin();
        ParameterList::iterator iter_p_end = vStarList.end();

        communicator.pack(&size, 1, tmpBuffer, bSize, p);

        // pack name sizes and chars
        for (; iter_p != iter_p_end; ++iter_p)
        {
          // move to id
          ++iter_p;

          size =(iter_p->tag()).length();
          communicator.pack(&size, 1, tmpBuffer, bSize, p);
          communicator.pack((iter_p->tag()).c_str(), size, tmpBuffer, bSize, p);
        }

        size = iStarList.size() / 2;
        iter_p = iStarList.begin();
        iter_p_end = iStarList.end();

        communicator.pack(&size, 1, tmpBuffer, bSize, p);

        // pack name sizes and chars
        for (; iter_p != iter_p_end; ++iter_p)
        {
          // move to id
          ++iter_p;

          size =(iter_p->tag()).length();
          communicator.pack(&size, 1, tmpBuffer, bSize, p);
          communicator.pack((iter_p->tag()).c_str(), size, tmpBuffer, bSize, p);
        }

        // tx local list to p0
        communicator.send(tmpBuffer, p, 0);

        // free memory
        delete [] tmpBuffer;
        vStarList.clear();
        iStarList.clear();

        // rx merged list from p0
        communicator.bcast(&bSize, 1, 0);
        tmpBuffer = new char[ bSize ];

        communicator.bcast(tmpBuffer, bSize, 0);

        // unpack and store
        p = 0;
        communicator.unpack(tmpBuffer, bSize, p, &size, 1);

        for (int j = 0; j < size; ++j)
        {
          communicator.unpack(tmpBuffer, bSize, p, &length, 1);
          std::string tmpStr(tmpBuffer + p, length);
          p += length;

          Util::Param param;
          param.setTag("V");
          param.setVal(1.0);
          vStarList.push_back(param);

          param.setTag(tmpStr);
          param.setVal(0.0);
          vStarList.push_back(param);
        }

        communicator.unpack(tmpBuffer, bSize, p, &size, 1);

        for (int j = 0; j < size; ++j)
        {
          communicator.unpack(tmpBuffer, bSize, p, &length, 1);
          std::string tmpStr(tmpBuffer + p, length);
          p += length;

          Util::Param param;
          param.setTag("I");
          param.setVal(1.0);
          iStarList.push_back(param);

          param.setTag(tmpStr);
          param.setVal(0.0);
          iStarList.push_back(param);
        }

        // free memory
        delete [] tmpBuffer;
      }
    }

    // advance iterators
    ++vStarPosition;
    ++iStarPosition;

    // append temporary lists to print block, erasing temporary lists
    variable_list.splice(vStarPosition, vStarList);
    variable_list.splice(iStarPosition, iStarList);
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::setSweepParameters
// Purpose       : Copy .DC or .STEP sweep parameters, and set up flags if
//                 sweeping the TEMP variable.
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 6/28/2013
//-----------------------------------------------------------------------------

void
OutputMgr::setSweepParameters(const std::vector<N_ANP_SweepParam> &step_sweep_parameters, const std::vector<N_ANP_SweepParam> &dc_sweep_parameters)
{
  if (!dc_sweep_parameters.empty())
  {
    dcParamVec_ = dc_sweep_parameters;
  }

  if (!step_sweep_parameters.empty())
  {
    stepParamVec_ =  step_sweep_parameters;
  }

  // check if one of the sweep variables is temperature.
  //(this only needs to be checked one time.)
  if (!step_sweep_parameters.empty())
  {
    for (std::vector<N_ANP_SweepParam>::const_iterator iterP = stepParamVec_.begin();iterP != stepParamVec_.end(); ++iterP)
    {
      if (equal_nocase((*iterP).name, "TEMP"))
      {
        tempSweepFlag_ = true;
      }
    }
  }

  // check if one of the DC sweep variables is temperature.
  //(this only needs to be checked one time.)
  if (!dc_sweep_parameters.empty() && !outputCalledBefore_)
  {
    for (std::vector<N_ANP_SweepParam>::const_iterator iterP = dcParamVec_.begin();iterP != dcParamVec_.end(); ++iterP)
    {
      if (equal_nocase((*iterP).name, "TEMP"))
      {
        tempSweepFlag_ = true;
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::fixupNodeNames
// Purpose       : obtain names of nodes from topology
// Special Notes : Why is this named "fixup"?  It does no fixing?
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 11/15/2013
//-----------------------------------------------------------------------------
void OutputMgr::fixupNodeNames()
{
  topology_->getNodeNames(allNodes_);
  topology_->getStateNodeNames(stateNodes_);
  topology_->getStoreNodeNames(storeNodes_);
  topology_->getExternNodeNames(externNodes_);
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::fixupPrintParameters
// Purpose       : Perform some .print line checks and munging, primarily
//                 dealing with use of V(*) and I(*)
// Special Notes : Mostly a re-wrapping of work that used to be done in the
//                 check_output function, which was too sprawling to remain
//                 as a single function.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 6/26/2013
//-----------------------------------------------------------------------------

void OutputMgr::fixupPrintParameters(PrintParameters &print_parameters)
{

  // Handle v(*) and i(*) print line options, by replacing with complete list
  removeStarVariables(print_parameters.variableList_, *pdsCommPtr_, allNodes_, externNodes_);

  // setup hdf5 output if requested
  if (hdf5FileNameGiven_)
  {
#ifdef Xyce_USE_HDF5
    enableHDF5Output();
#endif // Xyce_USE_HDF5
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::output
// Purpose       : Runs specified output commands
//
// Special Notes : ERK.  I've refactored this so that it receives STEP
//                 and DC sweep parameter information(mainly name and value).
//                 This function is called from the time integrator, which has
//                 all this info.
//
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/10/00
//-----------------------------------------------------------------------------
void OutputMgr::output(
  const double &                        time,
  const int                             stepNumber,
  const int                             maxStep,
  const std::vector<N_ANP_SweepParam> & step_sweep_parameters,
  const int                             dcNumber,
  const int                             maxDC,
  const std::vector<N_ANP_SweepParam> & dc_sweep_parameters,
  N_LAS_Vector *                        solnVecPtr,
  N_LAS_Vector *                        stateVecPtr,
  N_LAS_Vector *                        storeVecPtr,
  const std::vector<double>           & objectiveVec,
  const std::vector<double>           & dOdpVec, 
  const std::vector<double>           & dOdpAdjVec,
  const std::vector<double>           & scaled_dOdpVec, 
  const std::vector<double>           & scaled_dOdpAdjVec,
  bool                                  skipPrintLineOutput)
{
  // copy over time:
  circuitTime_ = time;

  // copy over the step sweep information:
  stepLoopNumber_ = stepNumber;
  maxParamSteps_ = maxStep;
  if (maxParamSteps_ > 0)
  {
    STEPEnabledFlag_ = true;
  }

  // copy the new values into the locally owned vector:
  if (!step_sweep_parameters.empty())
  {
    stepParamVec_ = step_sweep_parameters;
  }

  // copy over the dc sweep information:
  dcLoopNumber_ = dcNumber;
  maxDCSteps_ = maxDC;

  if (!analysisInterface_->getBlockAnalysisFlag()) {
    // This error test should not be used in the MPDE case, as at least
    // one of the initial conditions that can be
    // used by MPDE is a DC sweep.
    if (maxDCSteps_ > 0 && PRINTType_ == PrintType::TRAN)
    {
      Report::DevelFatal0() << "Print type is inconsistent, maxDCSteps = " << maxDCSteps_;
    }
  }

  if (!dc_sweep_parameters.empty())
  {
    dcParamVec_ = dc_sweep_parameters;

    // For now just have PRINTdcvalue, etc just be the first parameter.
    const N_ANP_SweepParam &firstParam = dc_sweep_parameters.front();
    PRINTdcname_   = firstParam.name;
    PRINTdcvalue_  = firstParam.currentVal;
    if (firstParam.type == "LIST")
    {
      PRINTdcstart_  = firstParam.valList[0];
      int size1 = firstParam.valList.size();
      PRINTdcstop_   = firstParam.valList[size1-1];
    }
    else
    {
      PRINTdcstart_  = firstParam.startVal;
      PRINTdcstop_   = firstParam.stopVal;
    }
  }

  // Check for temperature:
  circuitTemp_ = deviceInterface_->getParamAndReduce("TEMP");

  // call on the measure manager to update any active measures
  // update .measure functions before outputting the .print line
  if (PRINTType_ == PrintType::TRAN)
  {
    measureManager_.updateTranMeasures(circuitTime_, solnVecPtr, stateVecPtr, storeVecPtr);
    fourierManager_.updateFourierData(circuitTime_, solnVecPtr);
  }
  else
  {
    measureManager_.updateDcMeasures(dcParamVec_, solnVecPtr, stateVecPtr, storeVecPtr);
  }

  // Needs to pass skipPrintLineOutput
  if (!skipPrintLineOutput)
  {
    if (!activeOutputterStack_.empty())
    {
      for (std::vector<Outputter::Interface *>::const_iterator 
          it = activeOutputterStack_.back().begin(); 
          it != activeOutputterStack_.back().end(); ++it)
      {
        (*it)->output(solnVecPtr, stateVecPtr, storeVecPtr);
      }
    }
  }

  if (sensitivityFlag_)
  {
    if (sensParFullNames_.empty())
    {
      // include the objective function in the output
      if (!objectiveVec.empty())
      {
        sensParFullNames_.push_back(sensObjFunction_);
      }

      // set up the parameter names
      std::vector<std::string> parNamesDirect = sensParamNameVec_;
      std::vector<std::string> parNamesDirectScaled = sensParamNameVec_;
      std::vector<std::string> parNamesAdjoint = sensParamNameVec_;
      std::vector<std::string> parNamesAdjointScaled = sensParamNameVec_;

      if (!(dOdpVec.empty()))
      {
        for (std::vector<std::string>::iterator itN = parNamesDirect.begin(); 
            itN != parNamesDirect.end(); ++itN)
        {
          (*itN) = "  d" + sensObjFunction_ + "/d(" + (*itN) + ")";
          (*itN) += std::string("_Dir");
           sensParFullNames_.push_back((*itN));
        }
      }

      if (!(scaled_dOdpVec.empty()))
      {
        for (std::vector<std::string>::iterator itN = parNamesDirectScaled.begin(); 
            itN != parNamesDirectScaled.end(); ++itN)
        {
          (*itN) = "  d" + sensObjFunction_ + "/d(" + (*itN) + ")";
          (*itN) += std::string("_Dir_scaled");
           sensParFullNames_.push_back((*itN));
        }
      }
   
      if(!(dOdpAdjVec.empty()))
      {
        for (std::vector<std::string>::iterator itN = parNamesAdjoint.begin(); 
            itN != parNamesAdjoint.end(); ++itN)
        {
          (*itN) = "  d" + sensObjFunction_ + "/d(" + (*itN) + ")";
          (*itN) += std::string("_Adj");
           sensParFullNames_.push_back((*itN));
        }
      }

 
      if(!(scaled_dOdpAdjVec.empty()))
      {
        for (std::vector<std::string>::iterator itN = parNamesAdjointScaled.begin(); 
            itN != parNamesAdjointScaled.end(); ++itN)
        {
          (*itN) = "  d" + sensObjFunction_ + "/d(" + (*itN) + ")";
          (*itN) += std::string("_Adj_scaled");
           sensParFullNames_.push_back((*itN));
        }
      }

    }

    if (!activeOutputterStack_.empty())
    {
      for (std::vector<Outputter::Interface *>::const_iterator 
          it = activeOutputterStack_.back().begin(); 
          it != activeOutputterStack_.back().end(); ++it)
      {
        (*it)->outputSensitivity(sensParFullNames_, 
            objectiveVec, 
            dOdpVec, dOdpAdjVec, 
            scaled_dOdpVec, scaled_dOdpAdjVec, 
            solnVecPtr, stateVecPtr, storeVecPtr);
      }
    }
  }

#ifdef Xyce_USE_HDF5
  if (hdf5FileNameGiven_)
  {
    updateHDF5Output(solnVecPtr);
  }
#endif // Xyce_USE_HDF5

  // if any variables are needed for a response function
  // save them now.
  if (numResponseVars_ != 0)
  {
    saveResponseVarValues(solnVecPtr);
  }

  // transient or dc values for the objective function call.
  double arg1 = 0.0;
  double arg2 = 0.0;
  if (PRINTType_ == PrintType::TRAN || analysisInterface_->getAnalysisMgr()->getTransientFlag())
  {
    arg1 = circuitTime_;
  }
  else
  {
    if (dcParamVec_.size() > 0)
    {
      arg1 = dcParamVec_[0].currentVal;
    }
    if (dcParamVec_.size() > 1)
    {
      arg2 = dcParamVec_[1].currentVal;
    }
  }

  // if there are any simulation level functions(objective objects)
  // that need data from this time step, send it to them as well
  for (ObjectiveMap::iterator ob = objective_.begin(); ob != objective_.end(); ++ob)
  {
    (*ob).second.save(arg1, arg2, solnVecPtr, stateVecPtr, storeVecPtr);
  }

  outputCalledBefore_ = true;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::outputDCOP
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
void OutputMgr::outputDCOP(N_LAS_Vector & solnVec)
{
  if (outputOnceAlreadyFlag_)
  {
    return;
  }

  if (output_op_)
  {
    outputDCOP_restartFile(solnVec);
  }
  else if (saveFlag_)
  {
    // check for .IC
    outputIC_or_NODESET(solnVec);
  }

  outputOnceAlreadyFlag_ = true;

  return;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::outputDCOP
// Purpose       : Output DC operating point solution vars for use in subsequent
//                 runs to speed up DCOP
// Special Notes :
// Scope         : private
// Creator       : Dave Shirley, PSSI
// Creation Date : 04/10/06
//-----------------------------------------------------------------------------
void OutputMgr::outputDCOP_restartFile(N_LAS_Vector & solnVec)
{
  NodeNamePairMap::iterator name_i, name_end;
  int ind;
  std::ofstream opOut;

  if (!output_op_)
    return;

  name_i = allNodes_.begin();
  name_end = allNodes_.end();
  for (; name_i != name_end ; ++name_i)
  {
    ind =(*name_i).second.first;
   (*name_i).second.second = solnVec[ind];
  }

  if (getProcID() == 0)
  {
    std::string outputFileName;
    if (output_op_file_ == "")
      outputFileName = netListFilename_ + ".op";
    else
      outputFileName = output_op_file_;
    opOut.open(outputFileName.c_str());

    name_i = allNodes_.begin();

    for (; name_i != name_end ; ++name_i)
    {
      opOut <<(*name_i).first << "   " <<(*name_i).second.second << std::endl;
    }
  }

#ifdef Xyce_PARALLEL_MPI
  name_i = allNodes_.begin();
  int bufSize, mySize = sizeof(int);
  for (; name_i != name_end ; ++name_i)
  {
    mySize += sizeof(int) +(*name_i).first.size() + sizeof(double);
  }
  pdsCommPtr_->maxAll(&mySize, &bufSize, 1);
  std::vector<char> buffer(bufSize, '\0');
  double value;
  int pos, len, size, numNodes;
  if (getProcID() == 0)
  {
    int i, j;
    int one = 1;
    int two = 2;
    pdsCommPtr_->send(&one, 1, 1);
    for (i=1 ; i<getNumProcs() ; ++i)
    {
      pdsCommPtr_->recv(&size, 1, i);
      pdsCommPtr_->iRecv(&buffer[0], size, i);
      pdsCommPtr_->send(&two, 1, i);
      if (i < getNumProcs() - 1)
        pdsCommPtr_->send(&one, 1, i+1);
      pdsCommPtr_->waitAll();
      pos = 0;
      pdsCommPtr_->unpack(&buffer[0], bufSize, pos, &numNodes, 1);
      for (j=0 ; j<numNodes ; ++j)
      {
        pdsCommPtr_->unpack(&buffer[0], bufSize, pos, &len, 1);
        std::string name(&(buffer[pos]), len);
        pos += len;
        pdsCommPtr_->unpack(&buffer[0], bufSize, pos, &value, 1);
        opOut << name << "   " << value << std::endl;
      }
    }
  }
  else
  {
    size = allNodes_.size();
    pos = 0;
    pdsCommPtr_->pack(&size, 1, &buffer[0], bufSize, pos);
    name_i = allNodes_.begin();
    for (; name_i != name_end ; ++name_i)
    {
      size =(*name_i).first.size();
      pdsCommPtr_->pack(&size, 1, &buffer[0], bufSize, pos);
      pdsCommPtr_->pack((*name_i).first.c_str(), size, &buffer[0], bufSize, pos);
      value =(*name_i).second.second;
      pdsCommPtr_->pack(&value, 1, &buffer[0], bufSize, pos);
    }
    if (pos != mySize)
    {
      Report:: UserFatal0() << "Internal error 1 in OutputMgr::outputDCOP";
    }
    int flag;
    pdsCommPtr_->recv(&flag, 1, 0);
    if (flag != 1)
    {
      Report:: UserFatal0() << "Internal error 3 in OutputMgr::outputDCOP";
    }
    pdsCommPtr_->send(&mySize, 1, 0);
    pdsCommPtr_->recv(&flag, 1, 0);
    if (flag != 2)
    {
      Report:: UserFatal0() << "Internal error 3 in OutputMgr::outputDCOP";
    }
    pdsCommPtr_->send(&buffer[0], mySize, 0);
  }
  pdsCommPtr_->barrier();
  if (getProcID() == 0)
#endif
    opOut.close();
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::outputIC_or_NODESET
// Purpose       : Outputs either an *.ic file, either in .ic or .nodeset format.
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL.
// Creation Date :
//-----------------------------------------------------------------------------
void OutputMgr::outputIC_or_NODESET(N_LAS_Vector & solnVec)
{
  NodeNamePairMap::iterator name_i, name_end;
  int ind;
  std::ofstream saveOutputStream;

  name_i = allNodes_.begin();
  name_end = allNodes_.end();

  for (; name_i != name_end ; ++name_i)
  {
    ind =(*name_i).second.first;
   (*name_i).second.second = solnVec[ind];
  }

  if (getProcID() == 0)
  {
    std::string outputFileName;
    if (saveOutputFile_  == "")
      outputFileName = netListFilename_ + ".ic";
    else
      outputFileName = saveOutputFile_;

    saveOutputStream.open(outputFileName.c_str());

    name_i = allNodes_.begin();

    for (; name_i != name_end ; ++name_i)
    {
      // make sure this is not a branch current:
      ExtendedString tmpName((*name_i).first);
      tmpName.toUpper();
      std::string::size_type pos = tmpName.rfind("BRANCH");
      if (pos == std::string::npos)
      {
        saveOutputStream << saveFileType_ << " V(";
        saveOutputStream <<(*name_i).first << ") = " <<(*name_i).second.second << std::endl;
      }
    }
  }

#ifdef Xyce_PARALLEL_MPI
  // compute size of buffer for this processor and number of voltage variables being sent
  int numNodes=0;  // this value will be used if getProcID() != 0
  name_i = allNodes_.begin();
  int bufSize, mySize = sizeof(int);
  for (; name_i != name_end ; ++name_i)
  {
    // make sure this is not a branch current, which will not be sent:
    ExtendedString tmpName((*name_i).first);
    tmpName.toUpper();
    std::string::size_type pos = tmpName.rfind("BRANCH");
    if (pos == std::string::npos)
    {
      mySize += sizeof(int) +(*name_i).first.size() + sizeof(double);
      ++numNodes;
    }
  }
  pdsCommPtr_->maxAll(&mySize, &bufSize, 1);
  std::vector<char> buffer(bufSize, '\0');
  double value;
  int pos, len, size;
  if (getProcID() == 0)
  {
    int i, j;
    int one = 1;
    int two = 2;
    pdsCommPtr_->send(&one, 1, 1);
    for (i=1 ; i<getNumProcs() ; ++i)
    {
      pdsCommPtr_->recv(&size, 1, i);
      pdsCommPtr_->iRecv(&buffer[0], size, i);
      pdsCommPtr_->send(&two, 1, i);
      if (i < getNumProcs() - 1)
        pdsCommPtr_->send(&one, 1, i+1);
      pdsCommPtr_->waitAll();
      pos = 0;
      pdsCommPtr_->unpack(&buffer[0], bufSize, pos, &numNodes, 1);
      for (j=0 ; j<numNodes ; ++j)
      {
        pdsCommPtr_->unpack(&buffer[0], bufSize, pos, &len, 1);
        std::string name(&(buffer[pos]), len);
        pos += len;
        pdsCommPtr_->unpack(&buffer[0], bufSize, pos, &value, 1);
        saveOutputStream << saveFileType_ << " V(";
        saveOutputStream << name << ") = " << value << std::endl;
      }
    }
  }
  else
  {
    pos = 0;
    pdsCommPtr_->pack(&numNodes, 1, &buffer[0], bufSize, pos);
    name_i = allNodes_.begin();
    for (; name_i != name_end ; ++name_i)
    {
      // make sure this is not a branch current, which will not be sent:
      ExtendedString tmpName((*name_i).first);
      tmpName.toUpper();
      std::string::size_type tmpPos = tmpName.rfind("BRANCH");
      if (tmpPos == std::string::npos)
      {
        size =(*name_i).first.size();
        pdsCommPtr_->pack(&size, 1, &buffer[0], bufSize, pos);
        pdsCommPtr_->pack((*name_i).first.c_str(), size, &buffer[0], bufSize, pos);
        value =(*name_i).second.second;
        pdsCommPtr_->pack(&value, 1, &buffer[0], bufSize, pos);
      }
    }
    if (pos != mySize)
    {
      Report:: UserFatal0() << "Internal error 1 in OutputMgr::outputIC_or_NODESET";
    }
    int flag;
    pdsCommPtr_->recv(&flag, 1, 0);
    if (flag != 1)
    {
      Report:: UserFatal0() << "Internal error 2 in OutputMgr::outputIC_or_NODESET";
    }
    pdsCommPtr_->send(&mySize, 1, 0);
    pdsCommPtr_->recv(&flag, 1, 0);
    if (flag != 2)
    {
      Report:: UserFatal0() << "Internal error 3 in OutputMgr::outputIC_or_NODESET";
    }
    pdsCommPtr_->send(&buffer[0], mySize, 0);
  }
  pdsCommPtr_->barrier();
  if (getProcID() == 0)
#endif
    saveOutputStream.close();
}


//-----------------------------------------------------------------------------
// Function      : OutputMgr::inputDCOP_
// Purpose       : Attempt to input better starting estimate of DC operating
//                 point from a previous run
// Special Notes :
// Scope         : private
// Creator       : Dave Shirley, PSSI
// Creation Date : 04/11/06
//-----------------------------------------------------------------------------
bool OutputMgr::inputDCOP_(N_LAS_Vector & nextSolnVec,
                                 N_LAS_Vector & flagSolnVec)
{
  int ind;
  double value;
  std::ifstream opIn;
  bool success = false;
  bool fileFound = false;
  int nodesMatched = 0;
  int totalNodes = 0;
  std::set<std::string> notFound;
  std::set<std::string> notMatched;
  std::set<std::string> matched;
  std::set<std::string>::iterator nm_i;
  std::set<std::string>::iterator nm_end;
  std::string inputFileName;

  totalNodes = allNodes_.size();
  if (!input_op_)
    return success;

  flagSolnVec.putScalar(-1.0);
  nextSolnVec.putScalar(0.0);

  NodeNamePairMap::iterator nodes_i = allNodes_.begin();
  NodeNamePairMap::iterator nodes_end = allNodes_.end();
  for (; nodes_i != nodes_end ; ++nodes_i)
  {
    flagSolnVec[(*nodes_i).second.first] = 0;
  }

#ifdef Xyce_PARALLEL_MPI
  int bufSize = 100000;
  std::vector<char> buffer(bufSize, '\0');
  std::vector<std::string> candidates;
  std::vector<int> candidateFound;
  std::vector<int> allFound;
  int numCandidates = 0;;
  int pos, len, i;

  if (getProcID() == 0)
  {
    pos = 0;
#endif
    if (input_op_file_ == "")
    {
      inputFileName = netListFilename_ + ".op";
    }
    else
    {
      inputFileName = input_op_file_;
    }
    opIn.open(inputFileName.c_str(), std::ios::in);

    if (opIn.good())
    {
     Xyce::dout() << "Reading in operating point initial estimate from: " << inputFileName << std::endl;
      fileFound = true;
      while (1)
      {
        std::string node("");
        opIn >> node;
        if (node == "")
          break;
        opIn >> value;

        if (allNodes_.find(node) != allNodes_.end())
        {
          ind = allNodes_[node].first;
          nextSolnVec[ind] = value;
          flagSolnVec[ind] = 1;
          if (matched.find(node) == matched.end())
          {
            matched.insert(node);
            nodesMatched++;
          }
        }
        else
        {
#ifdef Xyce_PARALLEL_MPI
          len = node.size();
          if (pos+sizeof(int)+len+sizeof(double) >= bufSize)
          {
            numCandidates = candidates.size();
            pdsCommPtr_->bcast(&pos, 1, 0);
            pdsCommPtr_->bcast(&numCandidates, 1, 0);
            pdsCommPtr_->bcast(&buffer[0], pos, 0);
            candidateFound.resize(numCandidates);
            allFound.resize(numCandidates);
            for (i=0 ; i<numCandidates ; i++)
            {
              candidateFound[i] = 0;
            }
            pdsCommPtr_->sumAll(&candidateFound[0], &allFound[0], numCandidates);
            for (i=0 ; i<numCandidates ; i++)
            {
              if (allFound[i] == 0)
                notFound.insert(candidates[i]);
            }
            pos = 0;
            candidates.clear();
          }
          candidates.push_back(node);
          pdsCommPtr_->pack(&len, 1, &buffer[0], bufSize, pos);
          pdsCommPtr_->pack(node.c_str(), len, &buffer[0], bufSize, pos);
          pdsCommPtr_->pack(&value, 1, &buffer[0], bufSize, pos);
#else
          notFound.insert(node);
#endif
        }
      }
      opIn.close();
    }
#ifdef Xyce_PARALLEL_MPI
    if (!fileFound)
      pos = -1;
    pdsCommPtr_->bcast(&pos, 1, 0);
    if (pos > 0)
    {
      numCandidates = candidates.size();
      pdsCommPtr_->bcast(&numCandidates, 1, 0);
      pdsCommPtr_->bcast(&buffer[0], pos, 0);
      candidateFound.resize(numCandidates);
      allFound.resize(numCandidates);
      for (i=0 ; i<numCandidates ; i++)
      {
        candidateFound[i] = 0;
      }
      pdsCommPtr_->sumAll(&candidateFound[0], &allFound[0], numCandidates);
      for (i=0 ; i<numCandidates ; i++)
      {
        if (allFound[i] == 0)
          notFound.insert(candidates[i]);
      }
      pos = 0;
      pdsCommPtr_->bcast(&pos, 1, 0);
    }
  }
  else
  {
    int bufLen;
    pdsCommPtr_->bcast(&bufLen, 1, 0);
    while (bufLen > 0)
    {
      pdsCommPtr_->bcast(&numCandidates, 1, 0);
      pdsCommPtr_->bcast(&buffer[0], bufLen, 0);
      if (bufLen > 0)
      {
        candidateFound.resize(numCandidates);
        allFound.resize(numCandidates);
        for (i=0 ; i<numCandidates ; i++)
        {
          candidateFound[i] = 0;
        }
        pos = 0;
        i = 0;
        while (pos < bufLen)
        {
          pdsCommPtr_->unpack(&buffer[0], bufSize, pos, &len, 1);
          std::string name(&(buffer[pos]), len);
          pos += len;
          pdsCommPtr_->unpack(&buffer[0], bufSize, pos, &value, 1);
          if (allNodes_.find(name) != allNodes_.end())
          {
            ind = allNodes_[name].first;
            nextSolnVec[ind] = value;
            flagSolnVec[ind] = 1;
            if (matched.find(name) == matched.end())
            {
              matched.insert(name);
              nodesMatched++;
            }
            candidateFound[i] = 1;
          }
          i++;
        }
        pdsCommPtr_->sumAll(&candidateFound[0], &allFound[0], numCandidates);
        pdsCommPtr_->bcast(&bufLen, 1, 0);
      }
    }
  }
  nextSolnVec.importOverlap();
  flagSolnVec.importOverlap();
  int myNodesMatched = nodesMatched;
  pdsCommPtr_->sumAll(&myNodesMatched, &nodesMatched, 1);
  int myTotalNodes = totalNodes;
  pdsCommPtr_->sumAll(&myTotalNodes, &totalNodes, 1);
#endif

  op_found_ = nodesMatched;
  total_soln_ = totalNodes;
  std::set<std::string>::iterator matched_i = matched.begin();
  std::set<std::string>::iterator matched_end = matched.end();
  for (; matched_i != matched_end ; ++matched_i)
  {
    opData_[*matched_i].first = allNodes_[*matched_i].first;
    opData_[*matched_i].second = nextSolnVec[opData_[*matched_i].first];
  }


  nodes_i = allNodes_.begin();
  for (; nodes_i != nodes_end ; ++nodes_i)
  {
    if (matched.find((*nodes_i).first) == matched.end())
      notMatched.insert((*nodes_i).first);
  }

#ifdef Xyce_PARALLEL_MPI
  int numBad = notMatched.size();
// Xyce::dout() << "PE: " << getProcID() << " has " << numBad << " mismatched nodes" << std::endl;
  int myFlag, minProc;
  int myBufSize = 0;
  if (numBad > 0 && getProcID() > 0)
  {
    myFlag = getProcID();
    myBufSize = numBad*sizeof(int);
    nm_i = notMatched.begin();
    nm_end = notMatched.end();
    for (; nm_i != nm_end ; ++nm_i)
      myBufSize +=(*nm_i).size();
  }
  else
  {
    myFlag = getNumProcs();
  }
  pdsCommPtr_->maxAll(&myBufSize, &bufSize, 1);
  if (bufSize > 0)
  {
    buffer.resize(bufSize);
    if (myBufSize > 0)
    {
      pos = 0;
      nm_i = notMatched.begin();
      nm_end = notMatched.end();
      for (; nm_i != nm_end ; ++nm_i)
      {
        len =(*nm_i).size();
        pdsCommPtr_->pack(&len, 1, &buffer[0], bufSize, pos);
        pdsCommPtr_->pack((*nm_i).c_str(), len, &buffer[0], bufSize, pos);
      }
    }
    pdsCommPtr_->minAll(&myFlag, &minProc, 1);
    while (minProc < getNumProcs())
    {
      if (getProcID() == minProc)
      {
        pdsCommPtr_->send(&myBufSize, 1, 0);
        pdsCommPtr_->send(&buffer[0], myBufSize, 0);
        myFlag = getNumProcs();
      }
      else if (getProcID() == 0)
      {
        pdsCommPtr_->recv(&myBufSize, 1, minProc);
        pdsCommPtr_->recv(&buffer[0], myBufSize, minProc);
        pos = 0;
        while (pos < myBufSize)
        {
          pdsCommPtr_->unpack(&buffer[0], bufSize, pos, &len, 1);
          std::string name(&(buffer[pos]), len);
          pos += len;
          notMatched.insert(name);
        }
      }
      pdsCommPtr_->minAll(&myFlag, &minProc, 1);
    }
  }

  if (getProcID() == 0)
  {
#endif
    if (fileFound)
    {
      if (totalNodes > nodesMatched)
       Xyce::dout() << "DCOP restart:  Initialized " << nodesMatched << " of a possible " << totalNodes << " nodes" << std::endl;
      else
       Xyce::dout() << "DCOP restart:  All " << totalNodes << " nodes initialized" << std::endl;

      if (notFound.size() > 0)
      {
       Xyce::dout() << "DCOP restart:  Nodes specified in " << inputFileName << " but not present in circuit:" << std::endl;
        nm_i = notFound.begin();
        nm_end = notFound.end();
        for (; nm_i != nm_end ; ++nm_i)
         Xyce::dout() << *nm_i << std::endl;
      }
      if (notMatched.size() > 0)
      {
       Xyce::dout() << "DCOP restart:  Nodes present in circuit, but not specified in " << inputFileName << ":" << std::endl;
        nm_i = notMatched.begin();
        nm_end = notMatched.end();
        for (; nm_i != nm_end ; ++nm_i)
         Xyce::dout() << *nm_i << std::endl;
      }
    }
#ifdef Xyce_PARALLEL_MPI
  }
#endif

  if (nodesMatched > 0)
    success = true;

  return success;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::getICData
// Purpose       : Provide OP data for LOCA
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 05/04/06
//-----------------------------------------------------------------------------
NodeNamePairMap & OutputMgr::getICData(int & op_found, std::string & icType)
{
  op_found = op_found_;
  if (icType_==0)
  {
    icType = "DCOP_RESTART";
  }
  else if (icType_==1)
  {
    icType = "IC";
  }
  else if (icType_==2)
  {
    icType = "NODESET";
  }
  else
  {
    icType = "";
  }

  return opData_;
}


//-----------------------------------------------------------------------------
// Function      : OutputMgr::setupIC_or_NODESET
// Purpose       :
// Special Notes : Assumes that allNodes_ is set up.
// Scope         : private
// Creator       : Eric Keiter
// Creation Date : 09/13/07
//-----------------------------------------------------------------------------
bool OutputMgr::setupIC_or_NODESET(N_LAS_Vector & nextSolnVec,
                                N_LAS_Vector & flagSolnVec,
                                bool & useFlag,
                                std::string & icType,
                                std::vector<Util::OptionBlock> & initBlockVec)
{
  int lid(0);
  bool success(false);

  std::set<std::string> notFoundInCkt;
  std::set<std::string> notSpecified;
  std::set<std::string> matched;
  std::set<std::string>::iterator nm_i;
  std::set<std::string>::iterator nm_end;

  int totalNodes(allNodes_.size());

  if (!useFlag) return success;

  flagSolnVec.putScalar(-1.0);
  nextSolnVec.putScalar(0.0);

  NodeNamePairMap::iterator nodes_i = allNodes_.begin();
  NodeNamePairMap::iterator nodes_end = allNodes_.end();
  for (; nodes_i != nodes_end ; ++nodes_i)
  {
    flagSolnVec[(*nodes_i).second.first] = 0;
  }

  // icblock loop
  int icBlockIndex(0);
  int icBlockEnd = initBlockVec.size();

  for (;icBlockIndex < icBlockEnd; ++icBlockIndex)
  {
    // param loop
    ParameterList::const_iterator itPar  = initBlockVec[icBlockIndex].getParams().begin();
    ParameterList::const_iterator endPar = initBlockVec[icBlockIndex].getParams().end();
    for (;itPar != endPar; ++itPar)
    {
      std::string node("");
      double value(0.0);

      // the first tag is always "V".  At this point I variables are not supported,
      // and there is an error trap for them in the OptionBlock.C file.
      ++itPar;
      node = itPar->tag();
      ++itPar;

      ExtendedString strValue(itPar->tag());
      if (strValue.isValue())
      {
        value = strValue.Value();
      }
      else
      {
        Report::UserFatal0() << "Problems processing " << icType << " values";
      }

      bool globalFound(false);
      bool localFound(false);
      NodeNamePairMap::iterator iterCI = allNodes_.find(node);
      if (iterCI != allNodes_.end())
      {
        lid = iterCI->second.first;
        nextSolnVec[lid] = value;
        flagSolnVec[lid] = 1;
        localFound = true;
      }

#ifdef Xyce_DEBUG_IC
      if (localFound)
      {
        Xyce::dout().width(10);Xyce::dout().precision(3);Xyce::dout().setf(std::ios::scientific);
        Xyce::dout()
#ifdef Xyce_PARALLEL_MPI
          << "procID="<<getProcID() << "  "
#endif
          << icType + " found, and set: V(" << node << "):  solution["<<lid<<"] = " << value << std::endl;
      }
#endif

#ifdef Xyce_PARALLEL_MPI
      // Take care of moving this info across all processors and getting
      // a global assessment of whether or not the specified node actually
      // exists in the ckt.
      double dGlobal=0.0;
      double dLocal=0.0;
      dLocal=localFound?1.0:0.0;
      pdsCommPtr_->barrier();
      pdsCommPtr_->sumAll(&dLocal, &dGlobal, 1);
      globalFound =(dGlobal != 0.0)?true:false;

#ifdef Xyce_DEBUG_IC
      if (globalFound)
        Xyce::dout() << "procID="<<getProcID() << "  globalFound = true" << std::endl;
      else
        Xyce::dout() << "procID="<<getProcID() << "  globalFound = false" << std::endl;
#endif

#else
      globalFound = localFound;
#endif

      if (globalFound)
      {
        success = true;

        if (matched.find(node) == matched.end())
        {
          matched.insert(node);
        }
      }
      else
      {
        notFoundInCkt.insert(node);
      }
    }
  }

#ifdef Xyce_PARALLEL_MPI
  nextSolnVec.importOverlap();
  flagSolnVec.importOverlap();
#endif

  op_found_ = matched.size();
  total_soln_ = totalNodes;
  std::set<std::string>::iterator matched_i = matched.begin();
  std::set<std::string>::iterator matched_end = matched.end();
  for (; matched_i != matched_end ; ++matched_i)
  {
    // only add this to the opData map if it is local to the processor.
    // The allNodes object only contains local stuff.  Otherwise there
    // isn't much point.
    NodeNamePairMap::iterator iterCI = allNodes_.find(*matched_i);
    //if (allNodes_.find(*matched_i) != allNodes_.end())
    if (iterCI != allNodes_.end())
    {
      opData_[*matched_i].first = iterCI->second.first;
      opData_[*matched_i].second = nextSolnVec[opData_[*matched_i].first];
    }
  }


#ifdef Xyce_DEBUG_IO
  // Do final accounting of which nodes have been set and which have not.

  // Identify nodes in the circuit which were not set.
  nodes_i = allNodes_.begin();
  for (; nodes_i != nodes_end ; ++nodes_i)
  {
    if (matched.find((*nodes_i).first) == matched.end())
      notSpecified.insert((*nodes_i).first);
  }

#ifdef Xyce_PARALLEL_MPI
  // ERK.  If running in parallel, make the notSpecified list a global list:
  int bufSize = 100000;
  std::vector<char> buffer(bufSize, '\0');
  int pos, len, size, numNodes;

  int numBad = notSpecified.size();
  int myFlag, minProc;
  int myBufSize = 0;
  if (numBad > 0 && getProcID() > 0)
  {
    myFlag = getProcID();
    myBufSize = numBad*sizeof(int);
    nm_i = notSpecified.begin();
    nm_end = notSpecified.end();
    for (; nm_i != nm_end ; ++nm_i)
    {
      myBufSize +=(*nm_i).size();
    }
  }
  else
  {
    myFlag = getNumProcs();
  }

  pdsCommPtr_->maxAll(&myBufSize, &bufSize, 1);
  if (bufSize > 0)
  {
    buffer.resize(bufSize);
    if (myBufSize > 0)
    {
      pos = 0;
      nm_i = notSpecified.begin();
      nm_end = notSpecified.end();
      for (; nm_i != nm_end ; ++nm_i)
      {
        len =(*nm_i).size();
        pdsCommPtr_->pack(&len, 1, &buffer[0], bufSize, pos);
        pdsCommPtr_->pack((*nm_i).c_str(), len, &buffer[0], bufSize, pos);
      }
    }
    pdsCommPtr_->minAll(&myFlag, &minProc, 1);
    while (minProc < getNumProcs())
    {
      if (getProcID() == minProc)
      {
        pdsCommPtr_->send(&myBufSize, 1, 0);
        pdsCommPtr_->send(&buffer[0], myBufSize, 0);
        myFlag = getNumProcs();
      }
      else if (getProcID() == 0)
      {
        pdsCommPtr_->recv(&myBufSize, 1, minProc);
        pdsCommPtr_->recv(&buffer[0], myBufSize, minProc);
        pos = 0;
        while (pos < myBufSize)
        {
          pdsCommPtr_->unpack(&buffer[0], bufSize, pos, &len, 1);
          std::string name(&(buffer[pos]), len);
          pos += len;
          notSpecified.insert(name);
        }
      }
      pdsCommPtr_->minAll(&myFlag, &minProc, 1);
    }
  }
#endif // mpi
#endif // debug io

  if (getProcID() == 0)
  {

#ifdef Xyce_DEBUG_IO
    if (totalNodes > matched.size())
    {
      Xyce::dout() << icType << ":  Initialized " << matched.size() << " of a possible " << totalNodes << " nodes" << std::endl;
    }
    else
    {
      Xyce::dout() << icType << ":  All " << totalNodes << " nodes initialized" << std::endl;
    }
#endif

    if (notFoundInCkt.size() > 0)
    {
      lout() << icType << ":  Nodes specified on ." << icType << " line, but not present in circuit. (ignoring):" << std::endl;

      nm_i = notFoundInCkt.begin();
      nm_end = notFoundInCkt.end();
      for (; nm_i != nm_end ; ++nm_i)
      {
        lout() << *nm_i << std::endl;
      }
    }

#ifdef Xyce_DEBUG_IO
    // Note: for a typical .IC line, this list of unspecified nodes will be
    // a very long list.  So, don't output except in debug mode.
    if (notSpecified.size() > 0)
    {
      dout() << icType << ":  Nodes present in circuit, but not specified on ." << icType << " line(ignoring):" << std::endl;

      nm_i = notSpecified.begin();
      nm_end = notSpecified.end();
      for (; nm_i != nm_end ; ++nm_i)
      {
        dout() << *nm_i << std::endl;
      }
    }
#endif // Xyce_DEBUG_IO

  }

  return success;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::setupInitialConditions
// Purpose       : This function is an umbrella function, under which the
//                 various types of initial conditions are potentially set up.
//                 These include, but aren't limited to .IC, .NODESET and
//                 .DCOP restart.
//
//                 In addition to setting these things up, it does some
//                 nominal checking to make sure that more than on IC hasn't
//                 been specified(ie DCOP restart and .IC can't both be set).
// Special Notes :
// Scope         : private
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 09/13/07
//-----------------------------------------------------------------------------
bool OutputMgr::setupInitialConditions
#if 0
  (N_LAS_Vector & solnVec, N_LAS_Vector & flagVec, int & icType)
#else
  (N_LAS_Vector & solnVec, N_LAS_Vector & flagVec)
#endif
{
    bool dcopRestart = false;
    bool dotIC = false;
    bool NODESET = false;
    icType_ = -1;

    if (input_op_ && ICflag_)
    {
      Report::UserFatal0() << "Cannot set both DCOP restart and .IC simultaneously.";
    }

    if (input_op_ && nodesetflag_)
    {
      Report::UserFatal0() << "Cannot set both DCOP restart and .NODESET simultaneously.";
    }

    if (ICflag_ && nodesetflag_)
    {
      Report::UserFatal0() << "Cannot set both .IC and .NODESET simultaneously.";
    }

    if (input_op_)
    {
      // check for dcop restart(the original)
      dcopRestart = inputDCOP_(solnVec, flagVec);
      if (dcopRestart) icType_ = 0;
    }
    else if (ICflag_)
    {
      // check for .IC
      std::string ictype("IC");
      dotIC = setupIC_or_NODESET(solnVec, flagVec, ICflag_, ictype, ICblockVec_);
      if (dotIC) icType_ = 1;
    }
    else if (nodesetflag_)
    {
    // check for .NODESET
      std::string ictype("NODESET");
      NODESET = setupIC_or_NODESET(solnVec, flagVec, nodesetflag_, ictype, nodesetblockVec_);
      if (NODESET) icType_ = 2;
    }

#if 0
    icType = icType_;
#endif

    return(dcopRestart || dotIC || NODESET);
}



#ifdef Xyce_Dakota
//-----------------------------------------------------------------------------
// Function      : OutputMgr::recordDakotaParNames
// Purpose       : Records parameter names for a Dakota optimization
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 03/23/06
//-----------------------------------------------------------------------------
void OutputMgr::recordDakotaParNames(std::vector<std::string> & names)
{
  responseNames_ = names;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::recordDakotaFileName
// Purpose       : Records input file name for a Dakota optimization
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 03/23/06
//-----------------------------------------------------------------------------
void OutputMgr::recordDakotaFileName(std::string & fileName)
{
  responseFileName_ = fileName;
  responseFileNameGiven_ = true;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::outputDakotaResults
// Purpose       : Outputs parameter results after Dakota optimization
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 03/23/06
//-----------------------------------------------------------------------------
void OutputMgr::outputDakotaResults()
{
  double val;
  int i, j, num = responseNames_.size();
  std::vector<bool> done(num, false);

  std::ofstream dakOut;
  dakOut.open(std::string(responseFileName_ + ".result").c_str());

  dakOut << "Dakota optimization results:" << std::endl;
  for (i=0 ; i<num ; ++i)
  {
    if (!done[i])
    {
      std::string::size_type cpos = responseNames_[i].find_first_of(':');
      if (cpos == std::string::npos)
      {
        deviceInterface_->getParamAndReduce(responseNames_[i], val);
        dakOut << i << " : " << responseNames_[i] << " = " << val << std::endl;
      }
      else
      {
        std::string prefix(responseNames_[i].substr(0, cpos));
        dakOut << "Results for " << prefix << ":" << std::endl;
        prefix += ":";
        int parNum = 0;
        for (j=i ; j<num ; ++j)
        {
          if (responseNames_[j].substr(0, cpos+1) == prefix)
          {
            done[j] = true;
            deviceInterface_->getParamAndReduce(responseNames_[j], val);
            if (parNum > 0)
              dakOut << "  ";
            if (parNum%4 == 0)
            {
              if (parNum > 0)
                dakOut << std::endl;
              dakOut << "+  ";
            }
            dakOut << responseNames_[j].substr(cpos+1) << " = " << val;
            parNum++;
          }
        }
        dakOut << std::endl;
      }
    }
  }
  dakOut.close();
}
#endif  // Xyce_Dakota

//-----------------------------------------------------------------------------
// Function      : OutputMgr::getVariableNames
// Purpose       :
// Special Notes :
// Scope         : private
// Creator       : Todd Coffey
// Creation Date : 07/28/09
//-----------------------------------------------------------------------------
std::vector<std::string> OutputMgr::getVariableNames()
{
  std::vector<std::string> names;
  for (ParameterList::const_iterator iterParam = printParameters_->variableList_.begin(); iterParam != printParameters_->variableList_.end(); ++iterParam)
  {
    names.push_back(iterParam->tag());
  }
  return names;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::outputAC
// Purpose       : .PRINT output for ac runs
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey
// Creation Date : 8/5/08
//-----------------------------------------------------------------------------
void OutputMgr::outputAC(
    double frequency, 
    const N_LAS_Vector * real_solution_vector, 
    const N_LAS_Vector *imaginary_solution_vector)
{
  circuitFrequency_ = frequency;

  measureManager_.updateAcMeasures(frequency, real_solution_vector, imaginary_solution_vector);
  if (!activeOutputterStack_.empty())
  {
    for (std::vector<Outputter::Interface *>::const_iterator 
        it = activeOutputterStack_.back().begin(); it != activeOutputterStack_.back().end(); ++it)
    {
      (*it)->outputAC(frequency, real_solution_vector, imaginary_solution_vector);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::outputHomotopy
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey
// Creation Date : 8/5/08
//-----------------------------------------------------------------------------
void OutputMgr::outputHomotopy(
    const std::vector<std::string> & parameter_names, 
    const std::vector<double> & param_values, 
    const N_LAS_Vector * solution_vector)
{
  if (!activeOutputterStack_.empty())
  {
    for (std::vector<Outputter::Interface *>::const_iterator 
        it = activeOutputterStack_.back().begin(); 
        it != activeOutputterStack_.back().end(); ++it)
    {
      (*it)->outputHomotopy(parameter_names, param_values, solution_vector);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::outputMORTF
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey
// Creation Date : 8/5/08
//-----------------------------------------------------------------------------
void OutputMgr::outputMORTF(
    bool origSystem, 
    const double & freq, 
    const Teuchos::SerialDenseMatrix<int, std::complex<double> >& H)
{
  if (!activeOutputterStack_.empty())
  {
    for (std::vector<Outputter::Interface *>::const_iterator 
        it = activeOutputterStack_.back().begin(); 
        it != activeOutputterStack_.back().end(); ++it)
    {
      (*it)->outputMORTF(origSystem, freq, H);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::outputROM
// Purpose       : Output reduced order model to file in Matrix Market format.
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist
// Creation Date : 6/7/12
//-----------------------------------------------------------------------------
void OutputMgr::outputROM(
      const Teuchos::SerialDenseMatrix<int, double>& Ghat,
      const Teuchos::SerialDenseMatrix<int, double>& Chat,
      const Teuchos::SerialDenseMatrix<int, double>& Bhat,
      const Teuchos::SerialDenseMatrix<int, double>& Lhat
     )
{
  if (getProcID() == 0)
  {
    // Eliminate the ".cir" suffix from the filename.
    // Leaving ".cir" in the netListFilename_ causes issues when reading the files back in.
    int pos = netListFilename_.find(".cir");
    std::string baseoutputname;
    if (pos > 0)
    {
      baseoutputname = netListFilename_.substr(0, pos);
    }
    else
    {
      baseoutputname = netListFilename_;
    }

    // Open files for Ghat, Chat, Bhat, and Lhat
    FILE *c_file, *g_file, *b_file, *l_file;
    MMIO::MM_typecode matcode;
    std::string cfile = baseoutputname + ".Chat";
    std::string gfile = baseoutputname + ".Ghat";
    std::string bfile = baseoutputname + ".Bhat";
    std::string lfile = baseoutputname + ".Lhat";
    c_file = fopen(cfile.c_str(), "w");
    g_file = fopen(gfile.c_str(), "w");
    b_file = fopen(bfile.c_str(), "w");
    l_file = fopen(lfile.c_str(), "w");
    if (c_file == NULL || g_file == NULL || b_file == NULL || l_file == NULL)
    {
      Report::DevelFatal0() << "Cannot open one of the ROM files for output: " 
        << cfile << ", " << gfile << ", " << bfile << ", " << lfile;
    }
    mm_initialize_typecode(&matcode);
    mm_set_matrix(&matcode);
    mm_set_array(&matcode);
    mm_set_general(&matcode);
    mm_set_real(&matcode);

    int ret = 0;

    // Write the headers
    ret = MMIO::mm_write_banner(g_file, matcode);
    ret = MMIO::mm_write_banner(c_file, matcode);
    ret = MMIO::mm_write_banner(b_file, matcode);
    ret = MMIO::mm_write_banner(l_file, matcode);

    // Write the matrix array sizes
    ret = MMIO::mm_write_mtx_array_size(g_file, Ghat.numRows(), Ghat.numCols());
    ret = MMIO::mm_write_mtx_array_size(c_file, Chat.numRows(), Chat.numCols());
    ret = MMIO::mm_write_mtx_array_size(b_file, Bhat.numRows(), Bhat.numCols());
    ret = MMIO::mm_write_mtx_array_size(l_file, Lhat.numRows(), Lhat.numCols());

    // Write Ghat
    for (int j=0; j<Ghat.numCols(); j++) {
      for (int i=0; i<Ghat.numRows(); i++) {
        fprintf(g_file, "%22.16e\n", Ghat(i, j));
      }
    }

    // Write Chat
    for (int j=0; j<Chat.numCols(); j++) {
      for (int i=0; i<Chat.numRows(); i++) {
        fprintf(c_file, "%22.16e\n", Chat(i, j));
      }
    }

    // Write Bhat
    for (int j=0; j<Bhat.numCols(); j++) {
      for (int i=0; i<Bhat.numRows(); i++) {
        fprintf(b_file, "%22.16e\n", Bhat(i, j));
      }
    }

    // Write Lhat
    for (int j=0; j<Lhat.numCols(); j++) {
      for (int i=0; i<Lhat.numRows(); i++) {
        fprintf(l_file, "%22.16e\n", Lhat(i, j));
      }
    }

    // Close the files
    fclose(g_file);
    fclose(c_file);
    fclose(b_file);
    fclose(l_file);
  } // end proc0 check
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::outputROM
// Purpose       : Output reduced order model to file in Matrix Market format.
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist
// Creation Date : 6/7/12
//-----------------------------------------------------------------------------
void OutputMgr::outputROM(
  const N_LAS_Matrix& Ghat,
  const N_LAS_Matrix& Chat,
  const Teuchos::SerialDenseMatrix<int, double>& Bhat,
  const Teuchos::SerialDenseMatrix<int, double>& Lhat)
{
  // Eliminate the ".cir" suffix from the filename.
  // Leaving ".cir" in the netListFilename_ causes issues when reading the files back in.
  int pos = netListFilename_.find(".cir");
  std::string baseoutputname;
  if (pos > 0)
  {
    baseoutputname = netListFilename_.substr(0, pos);
  }
  else
  {
    baseoutputname = netListFilename_;
  }

  // Create file string for Chat and Ghat
  std::string gfile = baseoutputname + ".Ghat";
  std::string cfile = baseoutputname + ".Chat";

  // Get Epetra_CrsMatrix objects from Ghat and Chat
  Epetra_CrsMatrix& epetraGhat =(const_cast<N_LAS_Matrix*>(&Ghat))->epetraObj();
  Epetra_CrsMatrix& epetraChat =(const_cast<N_LAS_Matrix*>(&Chat))->epetraObj();

  // Write out objects using EpetraExt
  EpetraExt::RowMatrixToMatrixMarketFile(gfile.c_str(), epetraGhat);
  EpetraExt::RowMatrixToMatrixMarketFile(cfile.c_str(), epetraChat);

  // Write out Bhat and Lhat.
  // NOTE:  Only do this on one processor if running in parallel.

  if (getProcID() == 0)
  {

    // Open files for Bhat, and Lhat
    FILE *b_file, *l_file;
    MMIO::MM_typecode matcode;

    std::string bfile = baseoutputname + ".Bhat";
    std::string lfile = baseoutputname + ".Lhat";

    b_file = fopen(bfile.c_str(), "w");
    l_file = fopen(lfile.c_str(), "w");
    if (b_file == NULL || l_file == NULL)
    {
      Report::DevelFatal0() << "Cannot open one of the ROM files for output: " 
        << bfile << ", " << lfile;
    }
    mm_initialize_typecode(&matcode);
    mm_set_matrix(&matcode);
    mm_set_array(&matcode);
    mm_set_general(&matcode);
    mm_set_real(&matcode);

    int ret = 0;

    // Write the headers
    ret = MMIO::mm_write_banner(b_file, matcode);
    ret = MMIO::mm_write_banner(l_file, matcode);

    // Write the matrix array sizes
    ret = MMIO::mm_write_mtx_array_size(b_file, Bhat.numRows(), Bhat.numCols());
    ret = MMIO::mm_write_mtx_array_size(l_file, Lhat.numRows(), Lhat.numCols());

    // Write Bhat
    for (int j=0; j<Bhat.numCols(); j++) {
      for (int i=0; i<Bhat.numRows(); i++) {
        fprintf(b_file, "%22.16e\n", Bhat(i, j));
      }
    }

    // Write Lhat
    for (int j=0; j<Lhat.numCols(); j++) {
      for (int i=0; i<Lhat.numRows(); i++) {
        fprintf(l_file, "%22.16e\n", Lhat(i, j));
      }
    }

    // Close the files
    fclose(b_file);
    fclose(l_file);
  } // end proc0 check
}


//-----------------------------------------------------------------------------
// Function      : OutputMgr::finishOutput
// Purpose       : Runs specified finish output commands
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 12/13/00
//-----------------------------------------------------------------------------
void OutputMgr::finishOutput()
{
    if (!activeOutputterStack_.empty())
    {
      for (std::vector<Outputter::Interface *>::const_iterator it = activeOutputterStack_.back().begin(); it != activeOutputterStack_.back().end(); ++it)
      {
        (*it)->finishOutput();
      }
    }

#ifdef Xyce_USE_HDF5
  if (hdf5FileNameGiven_)
  {
    closeHDF5Output();
  }
#endif // Xyce_USE_HDF5
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::resetOutput
// Purpose       : Call outputter reset functions
// Special Notes : In current implementation, few outputter classes actually
//                 do anything in their reset methods.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 6/28/2013
//-----------------------------------------------------------------------------
void OutputMgr::resetOutput()
{
  if (!activeOutputterStack_.empty())
  {
    for (std::vector<Outputter::Interface *>::const_iterator it = activeOutputterStack_.back().begin(); it != activeOutputterStack_.back().end(); ++it)
    {
      (*it)->resetOutput();
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::outputHB
// Purpose       : .print output for Harmonic Balance runs
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 6/5/2013
//-----------------------------------------------------------------------------
void OutputMgr::outputHB(
  const int                             stepNumber,
  const int                             maxStep,
  const std::vector<N_ANP_SweepParam> & step_sweep_parameters,
  const std::vector<double> &           timePoints,
  const std::vector<double> &           freqPoints,
  const N_LAS_BlockVector &             timeDomainSolnVec,
  const N_LAS_BlockVector &             freqDomainSolnVecReal,
  const N_LAS_BlockVector &             freqDomainSolnVecImaginary)
{
  // copy over the step sweep information:
  stepLoopNumber_ = stepNumber;
  maxParamSteps_ = maxStep;
  if (maxParamSteps_ > 0)
  {
    STEPEnabledFlag_ = true;
  }

  // copy the new values into the locally owned vector:
  std::vector<N_ANP_SweepParam>::const_iterator firstParam = step_sweep_parameters.begin();
  std::vector<N_ANP_SweepParam>::const_iterator lastParam = step_sweep_parameters.end();

  if (!(step_sweep_parameters.empty()))
  {
    stepParamVec_ = step_sweep_parameters;
  }

  if (!activeOutputterStack_.empty())
  {
    for (std::vector<Outputter::Interface *>::const_iterator it = activeOutputterStack_.back().begin(); it != activeOutputterStack_.back().end(); ++it)
    {
      (*it)->outputHB(timePoints, freqPoints, timeDomainSolnVec, freqDomainSolnVecReal, freqDomainSolnVecImaginary);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::outputMPDE
// Purpose       : .PRINT output for mpde runs
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 9/22/03
//-----------------------------------------------------------------------------
void OutputMgr::outputMPDE(double time, const N_LAS_Vector * solnVecPtr)
{
  if (!activeOutputterStack_.empty())
    for (std::vector<Outputter::Interface *>::const_iterator it = activeOutputterStack_.back().begin(); it != activeOutputterStack_.back().end(); ++it)
      (*it)->outputMPDE(time, solnVecPtr);
}



//-----------------------------------------------------------------------------
// Function      : OutputMgr::outputRESULT
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/30/04
//-----------------------------------------------------------------------------
void OutputMgr::outputRESULT(N_LAS_Vector * solnVecPtr, N_LAS_Vector * stateVecPtr, N_LAS_Vector * storeVecPtr)
{

#ifdef Xyce_DEBUG_IO
  Xyce::dout() << std::endl << "OutputMgr::outputRESULT" << std::endl;
#endif

  std::vector<N_ANP_SweepParam>::iterator iterParam;
  std::vector<N_ANP_SweepParam>::iterator firstParam=stepParamVec_.begin();
  std::vector<N_ANP_SweepParam>::iterator lastParam=stepParamVec_.end();

  Util::ExpressionData * expDataPtr = NULL;

  int ires=0;
  std::string delim("");
  int width = 17;
  int precision = 8;

  if (getProcID() == 0)
  {
    if (!RESULTinitialized_)
    {
      RESULTinitialized_ = true;

      std::string resultfilename("");
      if (netListFilename_  != "")
      {
        resultfilename = netListFilename_ + ".res";
      }
      else
      {
        resultfilename = "output.res";
      }
      resultStreamPtr_ = new std::ofstream(resultfilename.c_str());

      resultStreamPtr_->setf(std::ios::scientific);
      resultStreamPtr_->precision(precision);

      if (!noIndex_)
      {
       (*resultStreamPtr_) << "STEP" ;
      }
      if (delim == "")(*resultStreamPtr_) << "   ";

      for (iterParam=firstParam; iterParam != lastParam;++iterParam)
      {
        if (delim == "")
         (*resultStreamPtr_) << "           ";
        else
         (*resultStreamPtr_) << printParameters_->delimiter_;

       (*resultStreamPtr_) << iterParam->name;
      }

      for (ires=0;ires< resultVector_.size(); ++ires)
      {
        if (delim == "")
         (*resultStreamPtr_) << "           ";
        else
         (*resultStreamPtr_) << printParameters_->delimiter_;

        expDataPtr = resultVector_[ires];
        (*resultStreamPtr_) << delim << expDataPtr->getExpression();
      }

     (*resultStreamPtr_) << "\n";

    } // RESULTinitialized_

    if (delim == "") resultStreamPtr_->width(8);
    else             resultStreamPtr_->width(0);
    resultStreamPtr_->setf(std::ios::left, std::ios::adjustfield);
    if (!noIndex_)
    {
     (*resultStreamPtr_) << stepLoopNumber_;
    }

    if (delim == "") resultStreamPtr_->width(width);

    for (iterParam=firstParam; iterParam != lastParam;++iterParam)
    {
      if (delim=="") { resultStreamPtr_->width(width); }
      else           {(*resultStreamPtr_)<<delim; }

     (*resultStreamPtr_) << iterParam->currentVal;
    }
  } // getProcID()

  for (ires=0;ires< resultVector_.size(); ++ires)
  {
    if (getProcID() == 0)
    {
      if (delim=="") { resultStreamPtr_->width(width); }
      else           {(*resultStreamPtr_)<<delim; }
    }

    expDataPtr = resultVector_[ires];
    double result = expDataPtr->evaluate(solnVecPtr, stateVecPtr, storeVecPtr);
    if (getProcID() == 0)
    {
     (*resultStreamPtr_) << result;
    }
  }

  if (getProcID() == 0)
  {
   (*resultStreamPtr_) << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::finishOutputSTEP
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/30/04
//-----------------------------------------------------------------------------
void OutputMgr::finishOutputSTEP()
{
  if (getProcID() == 0)
  {
    if (!activeOutputterStack_.empty())
    {
      for (std::vector<Outputter::Interface *>::const_iterator it = activeOutputterStack_.back().begin(); it != activeOutputterStack_.back().end(); ++it)
      {
        (*it)->finishOutputStep();
      }
    }

    // Deal with the result file:
    if (resultStreamPtr_)
    {
      (*resultStreamPtr_) << "End of Xyce(TM) Parameter Sweep" << std::endl;
    }

    if (resultStreamPtr_ != &Xyce::dout() && resultStreamPtr_)
    {
      delete resultStreamPtr_;
      resultStreamPtr_ = 0;
      // must re-set this flag, or Dakota runs that use Result files
      // can get confused.
      RESULTinitialized_ = false;
    }
  } // getProcID()
}


//-----------------------------------------------------------------------------
// Function      : OutputMgr::outputMacroResults
// Purpose       : if any post-process analysis was specified(like objective or measure)
//                 output results.  Called after simulation is over.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL Electrical and Microsystem Modeling
// Creation Date : 11/05/08
//-----------------------------------------------------------------------------
void OutputMgr::outputMacroResults()
{
  if (!objective_.empty())
  {
    Xyce::dout() << std::endl
      << " ***** Analysis Functions ***** " << std::endl
      << std::endl;

    std::map<std::string, Objective>::iterator ob = objective_.begin();
    std::map<std::string, Objective>::iterator ob_end = objective_.end();
    for (; ob != ob_end ; ++ob)
    {
      std::string name((*ob).first);
      double val = objective_[name].evaluate();
      Xyce::dout() << name << " = " << val << std::endl;
    }
  }

  // This is a null output stream that helps ensure a function that needs to be called
  // on every processor in parallel only outputs on one processor.
  Teuchos::oblackholestream outputBHS;
  std::ofstream outputFileStream;

  // Output the Measure results only if .measure is being performed on any variables.
  if (measureManager_.isMeasureActive())
  {
    // Output the Measure results to Xyce::dout().
    // Make sure the function gets called on all processors, but only one outputs it.
    if (getProcID() == 0)
    {
      measureManager_.outputResults( Xyce::dout() );
    }
    else
    {
      measureManager_.outputResults( outputBHS );
    }

    // Output the Measure results to file.
    if (getProcID() == 0)
    {
      // Adding "0" to the end of this string for the step number, which was always 0
      // for .measure previously anyways.
      std::string filename = netListFilename_ + ".mt0";
      outputFileStream.open( filename.c_str() );
      measureManager_.outputResults( outputFileStream, false );
      outputFileStream.close();
    }
    else {
      measureManager_.outputResults( outputBHS, false );
    }
  }

  // Output the Fourier results to file if Fourier analysis is being performed.
  // Make sure the function gets called on all processors, but only one outputs it.
  if (fourierManager_.isFourierActive())
  {
    if (getProcID() == 0)
    {
      std::string filename = netListFilename_ + ".four";
      outputFileStream.open( filename.c_str() );
      fourierManager_.outputResults( outputFileStream );
      outputFileStream.close();
    }
    else {
      fourierManager_.outputResults( outputBHS );
    }
  }

  // if the response list is not empty, try to dump those results to a file
  // a big limitation here is that all responses must be measure functions.
  // need to make this more flexible to include all solution vars and
  // objectives and results.
  if (!responseFunctionsRequested_.empty())
  {
    std::ofstream responseOFS;
    std::string outputResponseFilename;
    if (responseFileNameGiven_)
    {
      // do need the suffix in this case as a name was specified
      outputResponseFilename = responseFileName_;
    }
    else
    {
      outputResponseFilename = "response.out";
    }

    responseOFS.open(outputResponseFilename.c_str());

    std::vector< std::pair< std::string, std::string> >::iterator currRespItr = responseFunctionsRequested_.begin();
    std::vector< std::pair< std::string, std::string> >::iterator endRespItr = responseFunctionsRequested_.end();
    while (currRespItr != endRespItr)
    {
      double respvalue = -1.0;
      bool found = false;
      // need to parse name from XXX_X:name So find last ":" and extract
      // remaining string as value that dakota is looking for.
      std::string tempName=currRespItr->first;
      std::string::size_type beginingOfVarName = tempName.find_last_of(":");
      //beginingOfVarName++;

      if (beginingOfVarName != std::string::npos)
      {
        int numChars =(currRespItr->first).length() - beginingOfVarName;
        tempName.assign(currRespItr->first, beginingOfVarName+1, numChars);
      }
      //Xyce::dout() << " looking for " << currRespItr->first << std::endl;
      ExtendedString es(tempName);
      measureManager_.getMeasureValue(es.toUpper(), respvalue, found);
      if (! found)
      {
        // look for value in .objective functions
        std::map<std::string, Objective>::iterator ob = objective_.begin();
        std::map<std::string, Objective>::iterator ob_end = objective_.end();
        for (; ob != ob_end ; ++ob)
        {
          ExtendedString name((*ob).first);
          //Xyce::dout() << " looking for " << es.toUpper() << " in " << name.toUpper() << std::endl;
          if( es.toUpper() == name.toUpper() )
          {
            respvalue = objective_[name].evaluate();
            found = true;
          }
        }
      }
      responseOFS << respvalue   << "   " << tempName << std::endl;
      currRespItr++;
    }
    responseOFS.close();
  }
}


//-----------------------------------------------------------------------------
// Function      : OutputMgr::finalizeResponseVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 04/04/06
//-----------------------------------------------------------------------------
void OutputMgr::finalizeResponseVars()
{
  std::vector<std::string>::iterator currentVar;
  std::vector<std::string>::iterator endVar;
  int expValue = 0;

  // save the resuts of any end point variables if there are any
  if (numResponseVars_ != 0)
  {
    // save the independant variable
    int varNumber = 0;
    if (dcParamVec_.empty())  // DNS: is this a reliable way to determine if transient?
    {
      responseVarPtr_->at(varNumber) = circuitTime_;
    }
    else
    {
      responseVarPtr_->at(varNumber) = dcParamVec_[ dcLoopNumber_ ].currentVal;
    }
    varNumber++;
  }
}



// routines to tell the output manager which variables external programs will
// need as output.  By default we'll only remember the last timepoint
// or dc step unless asked to track all history.

//-----------------------------------------------------------------------------
// Function      : OutputMgr::registerResponseVars
// Purpose       : Create an objective from string submitted by external program
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 04/04/06
//-----------------------------------------------------------------------------

bool OutputMgr::registerResponseVars(const std::string & objString, RCP<std::vector< double > > varVectorPtr)
{
  bool result = true;
  responseVarPtr_ = varVectorPtr;

  ExtendedString sVal(objString);
  ParameterList::iterator pl_i;
  Util::Param parameter;

  sVal.toUpper();
  if (sVal.size() < 3 || sVal[1] != '(' || sVal[sVal.size()-1] != ')') {
    Report::DevelFatal0() << "OutputMgr::registerResponseVars: response var not of format V() or I(): '" << objString << "'";
  }
  numResponseVars_++;

  ParameterList pList;

  parameter.setTag(sVal.substr(0, 1));
  parameter.setVal(1.0);
  pList.push_back(parameter);

  parameter.setTag(sVal.substr(2, sVal.size()-3));
  parameter.setVal(0.0);
  pList.push_back(parameter);

  makeOps(*this, pList.begin(), pList.end(), std::back_inserter(responseVarList_));

  return result;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::saveResponseVarValues
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL, Parallel Computational Sciences
// Creation Date : 08/11/04
//-----------------------------------------------------------------------------
void OutputMgr::saveResponseVarValues(N_LAS_Vector * solnVecPtr)
{
  int expValue = 0;

  // save the resuts of any end point variables if there are any
  ParameterList vlist;
  vlist.push_back(Util::Param("", 0));
  ParameterList::iterator it_vlist = vlist.begin();

  // save the independant variable
  int varNumber = 0;
  if (dcParamVec_.empty())  // DNS: is this a reliable way to determine if transient?
  {
    responseVarPtr_->at(varNumber) = circuitTime_;
  }
  else
  {
    responseVarPtr_->at(varNumber) = dcParamVec_[ dcLoopNumber_ ].currentVal;
  }
  varNumber++;

  // loop over the response variable list
  for (Util::OpList::const_iterator it = responseVarList_.begin(); it != responseVarList_.end(); ++it)
  {
    double result = getValue(getCommPtr()->comm(), *(*it), solnVecPtr, 0, 0, 0).real();

    responseVarPtr_->at(varNumber) = result;
    varNumber++;
  }
}

//-----------------------------------------------------------------------------
// Function      : getMeasureValue
// Purpose       : get a .measure value
// Special Notes :
// Creator       : Dave Shirley, PSSI
// Creation Date : 4/29/10
//-----------------------------------------------------------------------------
void OutputMgr::getMeasureValue(const std::string &name, double &result, bool &found) const
{
  found = false;
  measureManager_.getMeasureValue(name, result, found);
}


//-----------------------------------------------------------------------------
// Function      : OutputMgr::remeasure
// Purpose       : Recompute measure functions based on existing Xyce output
// Special Notes :
// Creator       : Rich Schiek
// Creation Date : 4/8/13
//-----------------------------------------------------------------------------
void OutputMgr::remeasure()
{
  Xyce::dout() << "In OutputMgr::remeasure " << std::endl;

  // get file to open for remeasure
  std::string existingDataFile = commandLine_.getArgumentValue("-remeasure");
  Xyce::dout() << "file to reprocess through measure functions. " << existingDataFile << std::endl;

  // open file for reading
  // just support PRN format for now
  OutputFileBase * fileToReMeasure = new N_IO_OutputPrn();
  if (!fileToReMeasure->openFileForRead(existingDataFile))
  {
    // open failed.  report error and exit remeasure routine

    return;
  }

  // load data-names & allocate space for a line
  std::vector<std::string> fileVarNames;
  if (!(fileToReMeasure->getOutputVarNames(fileVarNames)))
  {
    // reading var names failed.  report error and exit remeasure routine.

  }

  // this code will need to move to the output handler classes as it will be
  // specific to the file type(i.e. some formats store the entire solution
  // regardless of what is on the print line.
  fileToReMeasure->convertOutputNamesToSolVarNames(fileVarNames);

//   Xyce::dout() << "Original var names: " << std::endl;
//   for (int i=0; i<fileVarNames.size(); i++)
//   {
//     Xyce::dout() << "\"" << fileVarNames[i] << "\", ";
//   }
//   Xyce::dout() << std::endl;
//
//   fileToReMeasure->convertOutputNamesToSolVarNames(fileVarNames);
//
//   Xyce::dout() << "extracted var names: " << std::endl;
//   for (int i=0; i<fileVarNames.size(); i++)
//   {
//     Xyce::dout() << "\"" << fileVarNames[i] << "\", ";
//   }
//   Xyce::dout() << std::endl;

  // assume analysis type is DC(Xyce hasn't processed the analysis type yet
  // when this function is called.  In the next loop if we find a column of
  // data for TIME then we can treat this as a transient analysis.

  PRINTType_ = PrintType::DC;
  int timeIndex=0;
  // set up allNodes map to map var names to indices in the solution vector
  int numVars = fileVarNames.size();
  for (int i=0; i<numVars; i++)
  {
    allNodes_[fileVarNames[i]]=std::make_pair(i, 0);
    ExtendedString tmpStr(fileVarNames[i]);
    tmpStr.toUpper();
    allNodes_[tmpStr]=std::make_pair(i, 0);
    // while scanning fileVarNames look for "TIME" as a name in the 0th or 1st
    // column.  We will use this as a key to figure out if the output file is
    // transient or DC data(no support for .measure in other analysis types yet)
    if ((i<2) &&(fileVarNames[i]=="TIME"))
    {
      PRINTType_ = PrintType::TRAN;
      timeIndex = i;
    }
  }

  // create an N_LAS_Vector to hold the data from the file.  This is
  // needed to support getPrintValue_ for evaluating the measure functions.
  std::vector<int> lbMap;
  lbMap.resize(numVars);
  for (int i=0;i<numVars; i++)
  {
    lbMap[i]=i;
  }
  N_PDS_ParMap * aParMapPtr = new N_PDS_ParMap(numVars, numVars, lbMap, 0, pdsCommPtr_);
  N_LAS_Vector * varValuesVecPtr = new N_LAS_Vector(*aParMapPtr);
  varValuesVecPtr->putScalar(0);

  // set up context for items in .measure line

  // run though lines in the file calling update measure as we go.
  while (fileToReMeasure->getOutputNextVarValues(varValuesVecPtr))
  {
    if (PRINTType_ == PrintType::TRAN)
    {
      circuitTime_=(*varValuesVecPtr)[timeIndex];
      measureManager_.updateTranMeasures(circuitTime_, varValuesVecPtr, 0, 0 );
    }
    else
    {
      measureManager_.updateDcMeasures(dcParamVec_, varValuesVecPtr, 0, 0 );
    }
    varValuesVecPtr->putScalar(0);

  }

  delete varValuesVecPtr;
  delete aParMapPtr;

  // This is a null output stream that helps ensure a function that needs to be called
  // on every processor in parallel only outputs on one processor.
  Teuchos::oblackholestream outputBHS;
  std::ofstream outputFileStream;

  // Output the Measure results to Xyce::dout().
  // Make sure the function gets called on all processors, but only one outputs it.
  if (getProcID() == 0)
  {
    measureManager_.outputResults( Xyce::dout() );
  }
  else
  {
    measureManager_.outputResults( outputBHS );
  }

  // Output the Measure results to file.
  if (getProcID() == 0)
  {
    // Adding "0" to the end of this string for the step number, which was always 0
    // for .measure previously anyways.
    std::string filename = netListFilename_ + ".mt0";
    outputFileStream.open( filename.c_str() );
    measureManager_.outputResults( outputFileStream );
    outputFileStream.close();
  }
  else {
    measureManager_.outputResults( outputBHS );
  }

  fileToReMeasure->closeFileForRead();
}

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
int getWidthFromStaticIndex(int index, const std::string &delim)
{
  int W = 0;
  if (delim == "")
  {
    if (index < 10000000)
    {
      W = 8;
    }
    else
    {
      W = 1;
      unsigned long SI = index;
      while (SI)
      {
        SI /= 10;
        ++W;
      }
    }
  }
  return W;
}

//-----------------------------------------------------------------------------
// Function      : gatherGlobalDeviceCount
//
// Purpose       : In parallel, gathers the local device counts and sums them
//                 into the global device count map.
//
// Special Notes : In serial, just copies the local map into the global one.
//
// Scope         : private
// Creator       : Eric R. Keiter, SNL
// Creation Date : 05/10/10
//-----------------------------------------------------------------------------
void gatherGlobalDeviceCount(
  Parallel::Machine                     comm,
  std::map<std::string, int> &          globalDeviceMap,
  const std::map<std::string, int> &    localDeviceMap)
{
  if (Parallel::is_parallel_run(comm))
  {
    std::map<std::string, int>::const_iterator dc;
    std::map<std::string, int>::const_iterator dc_end = localDeviceMap.end();

    int i, len = 0;
    const int size = Parallel::size(comm);
    const int rank = Parallel::rank(comm);

    std::set<std::string> known;

    int lowestKnown = size;
    if (!localDeviceMap.empty())
    {
      lowestKnown = rank;
    }
    Parallel::AllReduce(comm, MPI_MIN, &lowestKnown, 1);
    dc = localDeviceMap.begin();
    while (lowestKnown < size)
    {
      std::string curr;
      if (lowestKnown == rank)
      {
        curr =(*dc).first;
        len = curr.size();
      }

      // comm.bcast(&len, 1, lowestKnown);
      Parallel::Broadcast(comm, &len, 1, lowestKnown);
      curr.resize(len);
      // comm.bcast(&curr[0], len, lowestKnown);
      Parallel::Broadcast(comm, &curr[0], len, lowestKnown);
      dc = localDeviceMap.find(curr);

      int numDev = 0;
      if (dc != dc_end)
      {
        numDev =(*dc).second;
      }
      known.insert(curr);

      int numDevTot = 0;
      // comm.sumAll(&numDev, &numDevTot, 1);
      Parallel::AllReduce(comm, MPI_SUM, &numDev, &numDevTot, 1);
      globalDeviceMap[curr] = numDevTot;
      lowestKnown = size;
      dc = localDeviceMap.begin();
      for (; dc != dc_end ; ++dc)
      {
        if (known.find((*dc).first) == known.end())
        {
          lowestKnown = rank;
          curr =(*dc).first;
          break;
        }
      }
      // i = lowestKnown;
      // comm.minAll(&i, &lowestKnown, 1);
      Parallel::AllReduce(comm, MPI_MIN, &lowestKnown, 1);
    }
  }
  else
  {
    globalDeviceMap = localDeviceMap;
  }
}

//-----------------------------------------------------------------------------
// Function      : printDeviceCount
//
// Purpose       : takes the(passed) device count map, and formats a string
//                 that can be output, either via the error handler or Xyce::dout().
//
// Special Notes :
//
// Scope         : private
// Creator       : Eric R. Keiter, SNL
// Creation Date : 05/10/10
//-----------------------------------------------------------------------------
std::ostream &printDeviceCount(std::ostream &os, const std::map<std::string, int> & device_count_map)
{
  int maxLen = 15;

  int totDev = 0;
  for (std::map<std::string, int>::const_iterator dc = device_count_map.begin(); dc != device_count_map.end(); ++dc)
  {
    int len = (*dc).first.size();
    if (len > maxLen)
      maxLen = len;
    totDev += (*dc).second;
  }

  int width = 0;
  for (int i = totDev; i != 0; i /= 10)
    width++;

  for (std::map<std::string, int>::const_iterator dc = device_count_map.begin(); dc != device_count_map.end(); ++dc)
  {
    int len = (*dc).first.size();
    os << "       " << (*dc).first;
    for (int i = 0; i < maxLen - len + 1 ; ++i)
      os << " ";
    os.width(width);
    os.setf(std::ios::right);
    os << (*dc).second << "\n";
  }
  os << "       ";
  for (int i = 0; i < maxLen + width + 1; ++i)
    os << "-";

  os << "\n       Total Devices";
  for (int i = 0; i < maxLen - 12; ++i)
  {
    os << " ";
  }
  os.width(width);
  os.setf(std::ios::right);
  os << totDev;

  return os;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::printDeviceCounts
// Purpose       : Print device summary
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 08/15/06
//-----------------------------------------------------------------------------
std::ostream &printGlobalDeviceCounts(std::ostream &os, Parallel::Machine comm, const std::map<std::string,int> &local_device_count_map)
{
  std::map<std::string,int> global_device_count_map;

  gatherGlobalDeviceCount(comm, global_device_count_map, local_device_count_map);
  printDeviceCount(os, global_device_count_map);

  return os;
}

//-----------------------------------------------------------------------------
// Function      : printLocalDeviceCount
//
// Purpose       : If running in parallel, this function can be used to output
//                 the local count on each processor.
//
// Special Notes : debug only.
//
// Scope         : private
// Creator       : Eric R. Keiter
// Creation Date : 05/10/10
//-----------------------------------------------------------------------------
std::ostream &printLocalDeviceCount(std::ostream &os, Parallel::Machine comm, const std::map<std::string, int> & local_device_count_map)
{
  const int size = Parallel::size(comm);
  const int rank = Parallel::rank(comm);

  for (int p = 0; p < size; ++p)
  {
    Parallel::Barrier(comm);
    if (p == rank)
    {
      os << "\n\tDevice Count for Processor " << p << std::endl;

      printDeviceCount(os, local_device_count_map);

      os << std::endl;
    }
  }

  return os;
}

} // namespace IO
} // namespace Xyce

//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright(C) 2002-2013  Sandia Corporation
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
// Revision Number: $Revision: 1.374.2.21 $
//
// Revision Date  : $Date: 2013/10/03 17:23:43 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>
// ---------- Standard Includes ----------

#include <N_UTL_Misc.h>

#include <iostream>
#include <fstream>

#include <sstream>

#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif

//------- Xyce Includes ----------
#include <N_IO_CmdParse.h>
#include <N_IO_OutputMgr.h>

#include <N_TOP_Topology.h>

#include <N_DEV_DeviceInterface.h>

#include <N_LAS_System.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>
#include <N_LAS_BlockVector.h>

#include <N_ANP_AnalysisInterface.h>
#include <N_ANP_AnalysisManager.h>
#include <N_MPDE_Manager.h>

#include <N_PDS_Comm.h>
#include <N_PDS_ParMap.h>
#include <N_ERH_ErrorMgr.h>

#include <N_UTL_Version.h>
#include <N_UTL_Expression.h>
#include <N_UTL_LogStream.h>
#include <N_IO_Objective.h>
#include <N_UTL_ExpressionData.h>
#include <N_IO_mmio.h>

//------ Misc Includes --------
#include <Teuchos_as.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_oblackholestream.hpp>

#include <EpetraExt_RowMatrixOut.h>
#include <Epetra_CrsMatrix.h>

#include <N_IO_OutputFileBase.h>
#include <N_IO_OutputPrn.h>

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


//------ Extern Declarations --------

namespace Xyce {
namespace IO {

namespace {

void fixupExpressions(ParameterList::iterator begin, ParameterList::iterator end, OutputMgr &output_manager)
{
  // set up expressions.  This needs to be done on each processor.
  // expressions should have been fully resolved by this point.  We
  // can check that by look at the tag() and value() of the parameter.
  // if the N_UTL_Param has been fully resolved, then both the
  // hasExpressionTag() and hasExpressionValue() will be true.
  for (ParameterList::iterator it = begin ; it != end; ++it)
  {
    if ((*it).hasExpressionTag()) // check for expression.
    {
      if ((*it).getType() == EXPR)
      {
        (*it).setSimContextAndData(EXPRESSION, rcp(new N_UTL_ExpressionData((*it).ePtr(), &output_manager)));
      }
      else if( ((*it).getType() == DBLE) || ((*it).getType() == INT) )
      {
        (*it).setSimContextAndData(CONSTANT, (*it).dVal() );
      }
      else
      {
        (*it).setSimContextAndData(EXPRESSION, rcp(new N_UTL_ExpressionData((*it).tag(), &output_manager)));
      }
    }
  }
}

} // namespace <unnnamed>

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
OutputMgr * OutputMgr::factory(N_IO_CmdParse & cp)
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
OutputMgr::OutputMgr(N_IO_CmdParse & cp)
  : commandLine_(cp),
    ACFlag_(false),
    rawFlag_(false),
    tranFlag_(false),
    MPDEFlag_(false),
    DCSweepFlag_(false),
    HBFlag_(false),
    homotopyFlag_(false),
    enableHomotopyFlag_(false),
    format_(Format::STD),
    noIndex_(false),
    formatResult_(Format::STD),
    printParameters_(0),
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
  topPtr_(0),
  devPtr_(0),
  anaIntPtr_(0),
  initialOutputInterval_(0.0),
  filter_(0.0),
  filterGiven_(false),
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
  hdf5IndexValue_(0)
{
  activeOutputterStack_.push_back(std::vector<Outputter::Interface *>());
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
  for (std::vector< N_IO_OutputFileBase *>::iterator it = outputHandlersVec_.begin(); it != outputHandlersVec_.end(); ++it)
    delete (*it);

  for (OutputterMap::iterator it = outputterMap_.begin(); it != outputterMap_.end(); ++it)
    delete (*it).second;

  for (OpenPathStreamMap::iterator it = openPathStreamMap_.begin(); it != openPathStreamMap_.end(); ++it)
    delete (*it).second.second;
}

std::ostream *OutputMgr::openFile(const std::string &path, std::ios_base::openmode mode)
{
  OpenPathStreamMap::iterator it = openPathStreamMap_.find(path);

  if (path == "CONSOLE")
    return &std::cout;
  else if (it != openPathStreamMap_.end()) {
    ++(*it).second.first;
    return (*it).second.second;
  }
  else {
    std::ostream *os = new std::ofstream(path.c_str(), mode);
    openPathStreamMap_[path] = std::pair<int, std::ostream *>(1, os);

    if (!os->good())
    {
      std::string msg("Failure opening " + path + "\n");
      N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_FATAL_0, msg);
    }

    return os;
  }

}

std::ostream *OutputMgr::openFile(const std::string &path)
{
  return openFile(path, ios_base::out);
}

std::ostream *OutputMgr::openBinaryFile(const std::string &path)
{
  return openFile(path, ios_base::out | ios_base::binary);
}

int OutputMgr::closeFile(std::ostream *os)
{
  if (os == &std::cout)
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
// Function      : OutputMgr::check_output
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
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 06/01/05
//-----------------------------------------------------------------------------
void OutputMgr::prepareOutput(
  ANP_Analysis_Mode                     analysis_mode,
  const std::vector<N_ANP_SweepParam> & step_sweep_parameters,
  const std::vector<N_ANP_SweepParam> & dc_sweep_parameters)
{
  setSweepParameters(step_sweep_parameters, dc_sweep_parameters);

  // Setup rawfile if requested
  if (commandLine_.argExists("-r"))
  {
    enableOverrideRawOutput(*printParameters_);
  }
  else {
    switch (analysis_mode) {
      case ANP_MODE_DC_OP:
        break;

      case ANP_MODE_DC_SWEEP:
        enableDCSweepOutput();
        break;

      case ANP_MODE_TRANSIENT:
        enableTransientOutput();
        break;

      case ANP_MODE_MPDE:
//      enableTransientOutput(step_sweep_parameters, dc_sweep_parameters);
        enableMPDEOutput();
        break;

      case ANP_MODE_HB:
        enableHBOutput();
        break;

      case ANP_MODE_AC:
        enableACOutput();
        break;

      case ANP_MODE_MOR:
        break;
    }

    if (enableHomotopyFlag_)
      enableHomotopyOutput(analysis_mode);
  }

  if (outputterMap_.empty()) {
    Outputter::TimePrn *outputter = new Outputter::TimePrn(*this, *printParameters_);
    outputter->parse();
    outputterMap_[PrintType::TRAN] = outputter;
  }
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

namespace {
struct OutputMgr_STEPOptionsReg : public N_IO_PkgOptionsReg
{
  OutputMgr_STEPOptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const N_UTL_OptionBlock & options )
  { return outputManager_.registerSTEPOptions( options ); }

  OutputMgr &outputManager_;
};

struct OutputMgr_DCOptionsReg : public N_IO_PkgOptionsReg
{
  OutputMgr_DCOptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const N_UTL_OptionBlock & options )
  { return outputManager_.registerDCOptions( options ); }

  OutputMgr &outputManager_;
};

struct OutputMgr_DCOPOptionsReg : public N_IO_PkgOptionsReg
{
  OutputMgr_DCOPOptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const N_UTL_OptionBlock & options )
  { return outputManager_.registerDCOPOptions( options ); }

  OutputMgr &outputManager_;
};

struct OutputMgr_TranOptionsReg : public N_IO_PkgOptionsReg
{
  OutputMgr_TranOptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const N_UTL_OptionBlock & options )
  { return outputManager_.registerTranOptions( options ); }

  OutputMgr &outputManager_;
};

struct OutputMgr_MPDETranOptionsReg : public N_IO_PkgOptionsReg
{
  OutputMgr_MPDETranOptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const N_UTL_OptionBlock & options )
  { return outputManager_.registerMPDETranOptions( options ); }

  OutputMgr &outputManager_;
};

struct OutputMgr_HBOptionsReg : public N_IO_PkgOptionsReg
{
  OutputMgr_HBOptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const N_UTL_OptionBlock & options )
  { return outputManager_.registerHBOptions( options ); }

  OutputMgr &outputManager_;
};

struct OutputMgr_OptionsReg : public N_IO_PkgOptionsReg
{
  OutputMgr_OptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const N_UTL_OptionBlock & options )
  { return outputManager_.registerOutputOptions( options ); }

  OutputMgr &outputManager_;
};

struct OutputMgr_DeviceOptionsReg : public N_IO_PkgOptionsReg
{
  OutputMgr_DeviceOptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const N_UTL_OptionBlock & options )
  { return outputManager_.registerDeviceOptions( options ); }

  OutputMgr &outputManager_;
};

struct OutputMgr_PrintOptionsReg : public N_IO_PkgOptionsReg
{
  OutputMgr_PrintOptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const N_UTL_OptionBlock & options )
  { return outputManager_.registerPRINTSet
      ( options ); }

  OutputMgr &outputManager_;
};

struct OutputMgr_ICOptionsReg : public N_IO_PkgOptionsReg
{
  OutputMgr_ICOptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const N_UTL_OptionBlock & options )
  { return outputManager_.registerIC ( options ); }

  OutputMgr &outputManager_;
};

struct OutputMgr_NodeSetOptionsReg : public N_IO_PkgOptionsReg
{
  OutputMgr_NodeSetOptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const N_UTL_OptionBlock & options )
  { return outputManager_.registerNodeSet( options ); }

  OutputMgr &outputManager_;
};

struct OutputMgr_ResultOptionsReg : public N_IO_PkgOptionsReg
{
  OutputMgr_ResultOptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const N_UTL_OptionBlock & options )
  { return outputManager_.setRESULTParams( options ); }

  OutputMgr &outputManager_;
};

struct OutputMgr_ObjectiveOptionsReg : public N_IO_PkgOptionsReg
{
  OutputMgr_ObjectiveOptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const N_UTL_OptionBlock & options )
  { return outputManager_.setOBJECTIVEParams( options ); }

  OutputMgr &outputManager_;
};

struct OutputMgr_SaveOptionsReg : public N_IO_PkgOptionsReg
{
  OutputMgr_SaveOptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const N_UTL_OptionBlock & options )
  { return outputManager_.registerSave( options ); }

  OutputMgr &outputManager_;
};

struct OutputMgr_LoadOptionsReg : public N_IO_PkgOptionsReg
{
  OutputMgr_LoadOptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const N_UTL_OptionBlock & options )
  { return outputManager_.registerLoad( options ); }

  OutputMgr &outputManager_;
};

struct OutputMgr_OPOptionsReg : public N_IO_PkgOptionsReg
{
  OutputMgr_OPOptionsReg( OutputMgr &mgr )
    : outputManager_(mgr)
  {}

  bool operator()( const N_UTL_OptionBlock & options )
  { return outputManager_.registerOP ( options ); }

  OutputMgr &outputManager_;
};

struct OutputMgr_MeasureOptionsReg : public N_IO_PkgOptionsReg
{
  OutputMgr_MeasureOptionsReg( N_IO_MeasureManager &mgr )
    : measureManager_(mgr)
  {}

  bool operator()( const N_UTL_OptionBlock & options )
  { return measureManager_.addMeasure( options ); }

  N_IO_MeasureManager &measureManager_;
};

struct OutputMgr_FourierOptionsReg : public N_IO_PkgOptionsReg
{
  OutputMgr_FourierOptionsReg( N_IO_FourierMgr &mgr )
    : fourierManager_(mgr)
  {}

  bool operator()( const N_UTL_OptionBlock & options )
  { return fourierManager_.addFourierAnalysis( options ); }

  N_IO_FourierMgr &fourierManager_;
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
bool OutputMgr::registerPkgOptionsMgr(N_IO_PkgOptionsMgr &pkgOpt)
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
bool OutputMgr::registerDCOPOptions(const N_UTL_OptionBlock & option_block)
{
  for (ParameterList::const_iterator iterPL = option_block.getParams().begin(); iterPL != option_block.getParams().end(); ++iterPL)
  {
    if (iterPL->tag() == "INPUT")
    {
      input_op_ = true;
      input_op_file_ = iterPL->sVal();
    }
    else if (iterPL->tag() == "OUTPUT")
    {
      output_op_ = true;
      output_op_file_ = iterPL->sVal();
    }
    else if (iterPL->tag() == "TIME")
    {
      // do nothing, this will be handled in the time integrator.
    }
    else
    {
      std::string msg("Parameter not recognized in .DCOP line: ");
      msg += iterPL->tag();
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_WARNING_0, msg);
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
bool OutputMgr::registerDCOptions(const N_UTL_OptionBlock & option_block)
{
  for (ParameterList::const_iterator iterPL = option_block.getParams().begin(); iterPL != option_block.getParams().end(); ++iterPL)
    if (iterPL->tag() == "PARAM")
      dcParams_.push_back(iterPL->sVal());

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
bool OutputMgr::registerTranOptions(const N_UTL_OptionBlock & OB)
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
bool OutputMgr::registerMPDETranOptions(const N_UTL_OptionBlock & OB)
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
bool OutputMgr::registerHBOptions(const N_UTL_OptionBlock & OB)
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
bool OutputMgr::registerSTEPOptions(const N_UTL_OptionBlock & option_block)
{
  for (ParameterList::const_iterator iterPL = option_block.getParams().begin(); iterPL != option_block.getParams().end(); ++iterPL)
    if (iterPL->tag() == "PARAM")
      stepParams_.push_back(iterPL->sVal());

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
bool OutputMgr::registerDeviceOptions(const N_UTL_OptionBlock & option_block)
{
  for (ParameterList::const_iterator iter = option_block.getParams().begin(); iter != option_block.getParams().end(); ++iter)
    if (iter->tag() == "DETAILED_DEVICE_COUNTS")
      detailedDeviceCountFlag_ = iter->bVal();

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
bool OutputMgr::registerOutputOptions(const N_UTL_OptionBlock & OB)
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
      initialOutputInterval_ = iterPL->dVal();
      // look for optional time pairs.
      outputIntervalPairs_.clear();
      bool doneWithTimePairs = false;
      while ((iterPL != OB.getParams().end()) && !doneWithTimePairs)
      {
        if (iterPL->tag() == "TIME")
        {
          double t = iterPL->dVal();
          ++iterPL;
          double iv = iterPL->dVal();
          ++iterPL;
          outputIntervalPairs_.push_back(pair<double, double>(t, iv));
        }
        else
        {
          // didn't find a time pair so bail out of this loop.
          doneWithTimePairs=true;
        }
      }
    }

    // look for other option tags
    if (iterPL->tag()=="HDF5FILENAME")
    {
      hdf5FileNameGiven_=true;
      hdf5FileName_=iterPL->sVal();
    }

    // look for flag to turn off "End of Xyce(TM) Simulation" line
    if (iterPL->tag()=="PRINTENDOFSIMLINE")
    {
      printEndOfSimulationLine_=iterPL->bVal();
    }
    ++iterPL;
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
bool OutputMgr::registerDistribNodes(const N_UTL_OptionBlock & option_block)
{
  for (ParameterList::const_iterator iterPL = option_block.getParams().begin(); iterPL != option_block.getParams().end(); ++iterPL)
  {
    distribVnodeSet_.insert(iterPL->tag());
    distribVsrcSet_.insert(iterPL->sVal());
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
void OutputMgr::registerNodeDevNames(const std::set<std::string> * nodeNamesIn,
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
void OutputMgr::setExternalNetlistParams(std::vector< pair< std::string, std::string> > & externalNetlistParams)
{

  //  externalParams can contain section names from dakota like
  // variables 2, responses 4, derivatives 4.  If we don't find any sections tags
  // then we can assume that all the parameters are just variable names to set.
  // if we find tags, then use the "responses" section to record what response functions
  // need to be reported.
  std::string sectionTag="variables";
  std::vector< pair< std::string, std::string > >::iterator externalParamsIter = externalNetlistParams.begin();
  std::vector< pair< std::string, std::string > >::iterator externalParamsEnd = externalNetlistParams.end();
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
    return;

  std::map<std::string, bool> stringStat;
  std::map<std::string, bool> nodeStat;
  std::map<std::string, bool> instanceStat;

  for (ParameterList::const_iterator iterParam = printParameters_->variableList_.begin() ; iterParam != printParameters_->variableList_.end(); ++iterParam)
  {
    if ((*iterParam).getSimContext() == UNDEFINED)
    {
      std::string varType(iterParam->tag());

      bool done = false;
      std::vector<std::string> nodes;
      std::vector<std::string> instances;
      std::vector<std::string> leads;
      std::vector<std::string> strings;
      std::vector<std::string> special;
      if (iterParam->hasExpressionTag() )
      {
        // parameter starts with "{" but may not have been been parsed into an expression.
        // check if there is an underlying expression object with the parameter
        if (iterParam->getType() == EXPR)
        {
          iterParam->ePtr()->get_names(XEXP_NODE, nodes);
          iterParam->ePtr()->get_names(XEXP_INSTANCE, instances);
          iterParam->ePtr()->get_names(XEXP_LEAD, leads);
          iterParam->ePtr()->get_names(XEXP_STRING, strings);
          iterParam->ePtr()->get_names(XEXP_SPECIAL, special);  // special returns vars like TIME
        }
        else if ( ((iterParam->getType()) == DBLE) || ((iterParam->getType()) == INT) )
        {
          // an expression that has been resolved to a constant, usually on another processor
          // (*iterParam).setSimContextAndData( CONSTANT, iterParam->dVal());
        }
        instances.insert(instances.end(), leads.begin(), leads.end());
        strings.insert(strings.end(), special.begin(), special.end());  // add specials to strings
      }
      else
      {
        int numIndices = iterParam->iVal();
        if (numIndices > 0 &&(varType == "I" ||(varType.size() == 2 && varType[0] == 'I')))
        {
          // any devices found in this I(xxx) structure need to be communicated to the device manager
          // so that the lead currents can be calculated
          if (numIndices != 1)
          {
            std::string msg("Only one device argument allowed in I() in .print");
            N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
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
            std::string msg("Only one or two node arguments allowed in V() in .print");
            N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
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
            std::string msg("Only one device argument allowed in N() in .print");
            N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
          }
          ++iterParam;
        }
        else
        {
          strings.push_back(iterParam->tag());
        }
      }

      for (std::vector<std::string>::iterator iter_s = strings.begin(); iter_s != strings.end(); ++iter_s)
      {
        done = false;

        if (*iter_s == "TEMP" || *iter_s == "TIME" || *iter_s == "FREQUENCY" || *iter_s == "INDEX")
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
          done = devPtr_->getParam(*iter_s, result);

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


      for (std::vector<std::string>::iterator iter_s = nodes.begin(); iter_s != nodes.end(); ++iter_s)
      {
        bool tmpBool = false;
        // ERK: in rare cases, nodeNames_ can be NULL.  The function registerNodeDevNames is only called
        // on processors with greater than zero devices.  Occasionally, Xyce will give a processor
        // nothing to do.
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
          // ERK: in rare cases, devNames_ can be NULL.  The function registerNodeDevNames is only called
          // on processors with greater than zero devices.  Occasionally, Xyce will give a processor
          // nothing to do.
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
          // ERK: in rare cases, devNames_ can be NULL.  The function registerNodeDevNames is only called
          // on processors with greater than zero devices.  Occasionally, Xyce will give a processor
          // nothing to do.
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
        stat[i++] = 1;
      else
        stat[i++] = 0;
    }
    stat_i = instanceStat.begin();
    stat_end = instanceStat.end();
    for (; stat_i != stat_end ; ++stat_i)
    {
      if ((*stat_i).second)
        stat[i++] = 1;
      else
        stat[i++] = 0;
    }
    pdsCommPtr_->sumAll(&stat[0], &stat_g[0], totSize);
    i = 0;
    stat_i = nodeStat.begin();
    stat_end = nodeStat.end();
    for (; stat_i != stat_end ; ++stat_i)
    {
      if (stat_g[i++] > 0)
        (*stat_i).second = true;
    }
    stat_i = instanceStat.begin();
    stat_end = instanceStat.end();
    for (; stat_i != stat_end ; ++stat_i)
    {
      if (stat_g[i++] > 0)
        (*stat_i).second = true;
    }
  }

  // Generate message
  std::map<std::string, bool>::iterator stat_i, stat_end;
  std::string msg("Undefined item(s) in .print: ");
  int mlen = msg.size();
  stat_i = nodeStat.begin();
  stat_end = nodeStat.end();
  for (; stat_i != stat_end ; ++stat_i)
  {
    if (!(*stat_i).second)
    {
      msg += "V(" +(*stat_i).first + ") ";
    }
  }
  stat_i = instanceStat.begin();
  stat_end = instanceStat.end();
  for (; stat_i != stat_end ; ++stat_i)
  {
    if (!(*stat_i).second)
    {
      if (((*stat_i).first)[((*stat_i).first).size()-1] == '}')
      {
        msg += "I" +((*stat_i).first).substr(((*stat_i).first).size()-2, 1);
        msg += "(";
        msg +=((*stat_i).first).substr(0, ((*stat_i).first).size()-3);
        msg += ") ";
      }
      else
      {
        msg += "I(" +(*stat_i).first + ") ";
      }
    }
  }
  stat_i = stringStat.begin();
  stat_end = stringStat.end();
  for (; stat_i != stat_end ; ++stat_i)
  {
    if (!(*stat_i).second)
    {
      msg +=(*stat_i).first + " ";
    }
  }
  if (mlen < msg.size())
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL_0, msg);
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
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/17/07
//-----------------------------------------------------------------------------
void OutputMgr::delayedPrintLineDiagnostics()
{
  if (deviceParamCheck.empty())
    return;

  std::string badParams;
  for (std::set<std::string>::iterator iter = deviceParamCheck.begin(); iter != deviceParamCheck.end(); ++iter)
  {
    std::string s(*iter);
    double result = 0.0;
    bool done = devPtr_->getParam(s , result);
    if (!done)
    {
      badParams += *iter + "  ";
    }
  }

  if (badParams != "")
  {
    std::string msg("Undefined parameter(s) in .print: ");
    msg += badParams;
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL_0, msg);
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
bool OutputMgr::getOutputIntervals(double & initialInterval, std::vector< pair<double, double> > & intervalPairs)
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
bool OutputMgr::setRESULTParams(const N_UTL_OptionBlock & OB)
{
#ifdef Xyce_DEBUG_IO
  cout << "In N_IO OutputMgr::setRESULTParams" << std::endl;
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
    msg =  "N_IO OutputMgr::setRESULTParams\n";
    msg += "You have more than one result expression on a single line.\n";
    msg += "Each parameter needs its own line.\n";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  }

  N_UTL_ExpressionData * expDataPtr;

  for (it_tp = first; it_tp != last; ++it_tp)
  {
    N_UTL_Param expParam = *it_tp;

    if (!expParam.hasExpressionValue())
    {
      msg = "Non-Expression specified in .RESULT\n";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
    }
    else
    {
      // expression should have already been resolved.  Check for this
      // case before we try and create a new expression.
      if (expParam.getType() == EXPR)
      {
        expDataPtr = new N_UTL_ExpressionData(expParam.ePtr(), this);
      }
      else
      {
        expDataPtr = new N_UTL_ExpressionData(expParam.sVal(), this);
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
bool OutputMgr::setOBJECTIVEParams(const N_UTL_OptionBlock & OB)
{
#ifdef Xyce_DEBUG_IO
  cout << "In N_IO OutputMgr::setOBJECTIVEParams" << std::endl;
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
      if (objective.find(name) != objective.end())
      {
        std::string msg("Duplicate objective name: ");
        msg += name;
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
      }
      break;
    }
  }
  objective[name].initialize(OB, this);

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
bool OutputMgr::registerPRINTSet(const N_UTL_OptionBlock &print_block)
{
#ifdef Xyce_DEBUG_IO
  cout << " PRINTSET Registered::" << std::endl;
#endif

  PrintParameters print_parameters;

  ParameterList::const_iterator iterParam = print_block.getParams().begin();
  for (; iterParam != print_block.getParams().end(); ++iterParam)
  {

#ifdef Xyce_DEBUG_IO
    cout << "iterParam->tag = " << iterParam->tag() << std::endl;
#endif

    if (iterParam->tag() == "WIDTH") {
      print_parameters.streamWidth_ = iterParam->iVal();
    }

    else if (iterParam->tag() == "TYPE") {
      std::string s = iterParam->sVal();
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
        std::string msg("Unrecognized analysis type, " + s);
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
      }
      print_parameters.printType_ = PRINTType_;
    }

    else if (iterParam->tag() == "PRECISION") {
      print_parameters.streamPrecision_ = iterParam->iVal();
    }

    else if (iterParam->tag() == "FILTER")
    {
      filterGiven_ = true;
      filter_ = iterParam->dVal();
    }
    else if (iterParam->tag() == "FORMAT")
    {
      std::string s = iterParam->sVal();

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
        std::string msg("Unrecognized print format, " + s);
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
      }
      print_parameters.format_ = format_;
    }

    else if (iterParam->tag() == "TIMEWIDTH")
    {
      print_parameters.timeWidth_ = iterParam->iVal();
    }

    else if (iterParam->tag() == "TIMESCALEFACTOR")
    {
      print_parameters.outputTimeScaleFactor_ = iterParam->dVal();
    }

    else if (iterParam->tag() == "FILE")
    {
      // netListFilename_ should be the default unless FILE was set to
      // something other than "" which is its default value.
      if (iterParam->sVal() != "")
      {
        print_parameters.filename_ = iterParam->sVal();
      }
    }

    else if (iterParam->tag() == "DELIMITER")
    {
      if (iterParam->sVal() == "TAB") {
        print_parameters.delimiter_ = "\t";
      }

      else if (iterParam->sVal() == "COMMA") {
        print_parameters.delimiter_ = ",";
      }

      else if (iterParam->sVal() != "")
      {
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_WARNING_0,
                                "Invalid value of DELIMITER in .PRINT statment, ignoring");
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

  print_parameters.index_ = !noIndex_ && print_parameters.format_ != Format::PROBE && print_parameters.format_ != Format::TECPLOT
                            && print_parameters.format_ != Format::RAW && print_parameters.format_ != Format::RAW_ASCII;

  // Assemble the apropriate flavors of output variable lists based on the PRINT type and format
  if (PRINTType_ == PrintType::AC) {
    outputParameterMap_[OutputType::AC] = print_parameters;
    outputParameterMap_[OutputType::AC_IC] = print_parameters;
    outputParameterMap_[OutputType::HOMOTOPY] = print_parameters;
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
    printParameters_ = &outputParameterMap_[OutputType::TRAN];
  }
  else if (PRINTType_ == PrintType::DC) {
    outputParameterMap_[OutputType::DC] = print_parameters;
    outputParameterMap_[OutputType::HOMOTOPY] = print_parameters;
    printParameters_ = &outputParameterMap_[OutputType::DC];
  }
  else {
    std::string msg("Unrecognized PRINTType_");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL_0, msg);
  }

#ifdef Xyce_DEBUG_IO
  cout << "netListFilename_ = " << netListFilename_ << std::endl;
  cout << " PRINTSET Registered:: finished" << std::endl;
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
bool OutputMgr::registerIC(const N_UTL_OptionBlock & OB)
{
#ifdef Xyce_DEBUG_IO
  cout << " IC Registered::" << std::endl;
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
    cout << "iterParam->tag = " << iterParam->tag() << std::endl;
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
bool OutputMgr::registerNodeSet(const N_UTL_OptionBlock & OB)
{
#ifdef Xyce_DEBUG_IO
  cout << " NODESET Registered::" << std::endl;
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
    cout << "iterParam->tag = " << iterParam->tag() << std::endl;
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
bool OutputMgr::registerSave(const N_UTL_OptionBlock & OB)
{
#ifdef Xyce_DEBUG_IO
  cout << " SAVE Registered::" << std::endl;
#endif

  saveFlag_ = true;

  ExtendedString sval("");

  ParameterList::const_iterator iterPL = OB.getParams().begin();
  ParameterList::const_iterator iterPL_end = OB.getParams().end();

  while (iterPL != iterPL_end)
  {
#ifdef Xyce_DEBUG_IO
    cout << "iterPL->tag = " << iterPL->tag() << std::endl;
#endif
    if (iterPL->tag() == "TYPE")
    {
      sval = iterPL->sVal();
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
        std::string msg("Unrecognized type specified on .SAVE line.  Defaulting to .NODESET");
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_WARNING_0, msg);
      }
    }
    else if (iterPL->tag() == "FILE")
    {
      saveOutputFile_ = iterPL->sVal();
    }
    else if (iterPL->tag() == "TIME")
    {
      // do nothing, this will be handled in the time integrator.
    }
    else if (iterPL->tag() == "LEVEL")
    {
      sval = iterPL->sVal();
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
        std::string msg("LEVEL=TOP in .SAVE line not supported.  Ignoring. ");
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_WARNING_0, msg);
      }
      else
      {
        std::string msg("Unrecognized LEVEL specified in .SAVE line: ");
        msg += sval;
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_WARNING_0, msg);
      }
    }
    else
    {
      std::string msg("Parameter not recognized in .SAVE line: ");
      msg += iterPL->tag();
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_WARNING_0, msg);
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
bool OutputMgr::registerLoad(const N_UTL_OptionBlock & OB)
{
#ifdef Xyce_DEBUG_IO
  cout << " Load Registered::" << std::endl;
#endif

  loadFlag_ = true;
  std::string msg(".LOAD not supported yet.  Use .INCLUDE instead");
  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_WARNING_0, msg);

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
bool OutputMgr::registerOP(const N_UTL_OptionBlock & OB)
{
#ifdef Xyce_DEBUG_IO
  cout << " OP Registered::" << std::endl;
#endif

  //  std::string msg(".OP is still under construction and not fully supported yet. ");
  //  N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_WARNING_0, msg);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::enableHomotopyOutput
// Purpose       : turns on Homotopy output
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/10/00
//-----------------------------------------------------------------------------
void OutputMgr::enableHomotopyOutput(ANP_Analysis_Mode analysis_mode)
{
  // prepare output manager to write
  if (homotopyFlag_)
  {
    std::string msg("homotopyfile already initialized.  Contents may be overwritten.\n");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_WARNING_0, msg);
  }

  else
  {
    homotopyFlag_ = true;

    PrintParameters homotopy_print_parameters = outputParameterMap_[OutputType::HOMOTOPY];

    if (analysis_mode == ANP_MODE_TRANSIENT)
      homotopy_print_parameters.variableList_.push_front(N_UTL_Param("TIME", 0.0));
    if (homotopy_print_parameters.index_)
      homotopy_print_parameters.variableList_.push_front(N_UTL_Param("INDEX", 0.0));

    fixupPrintParameters(DOMAIN_TIME, homotopy_print_parameters);

    Outputter::Interface *outputter;
    if (homotopy_print_parameters.format_ == Format::STD) {
      outputter = new Outputter::HomotopyPrn(*this, homotopy_print_parameters);
    }
    // else if (format_ == Format::TECPLOT) {
    //   outputter = new Outputter::HomotopyTecPlot(*this, print_parameters);
    // }
    else
    {
      string msg("Homotopy output can only use tecplot or gnuplot format.");
      msg += "  Resetting to gnuplot format.\n";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_WARNING, msg);

      outputter = new Outputter::HomotopyPrn(*this, homotopy_print_parameters);
    }

    outputter->parse();
    outputterMap_[PrintType::HOMOTOPY] = outputter;
    addActiveOutputter(PrintType::HOMOTOPY);
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::enableACoutput
// Purpose       : turns on AC output
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/10/00
//-----------------------------------------------------------------------------
void OutputMgr::enableACOutput()
{
  if (getNumProcs() > 1) {
    std::string msg("Parallel processing is not supported in A/C analysis");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL_0, msg);
  }

  // prepare output manager to write
  if (ACFlag_)
  {
    std::string msg("ACfile already initialized.  Contents may be overwritten.\n");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_WARNING_0, msg);
  }

  else
  {
    ACFlag_ = true;

    PrintParameters freq_print_parameters = outputParameterMap_[OutputType::AC];
    PrintParameters ac_ic_print_parameters = outputParameterMap_[OutputType::AC_IC];

    if (freq_print_parameters.format_ != Format::PROBE)
    {
      freq_print_parameters.variableList_.push_front(N_UTL_Param("FREQUENCY", 0.0));
    }
    if (freq_print_parameters.index_)
    {
      freq_print_parameters.variableList_.push_front(N_UTL_Param("INDEX", 0.0));
    }
    
    if (ac_ic_print_parameters.format_ != Format::PROBE)
    {
      ac_ic_print_parameters.variableList_.push_front(N_UTL_Param("TIME", 0.0));
    }
    if (ac_ic_print_parameters.index_)
    {
      ac_ic_print_parameters.variableList_.push_front(N_UTL_Param("INDEX", 0.0));
    }

    if( (freq_print_parameters.format_ == Format::PROBE) ||
        (freq_print_parameters.format_ == Format::RAW)   ||
        (freq_print_parameters.format_ == Format::RAW_ASCII) )
    {
      // for probe and raw output formats the output has a complex 
      // type so we don't need to expand V(a)
      // into Re(V(a)) Im(V(a)).  Thus we override the default 
      // parameter expandComplexTypes to false.
      fixupPrintParameters(DOMAIN_FREQUENCY, freq_print_parameters, false);
    }
    else
    {
      fixupPrintParameters(DOMAIN_FREQUENCY, freq_print_parameters);
    }
    fixupPrintParameters(DOMAIN_TIME, ac_ic_print_parameters);

    Outputter::Interface *outputter_fd;
    Outputter::Interface *outputter_prn;
    if (freq_print_parameters.format_ == Format::STD) {
      outputter_fd = new Outputter::FrequencyPrn(*this, freq_print_parameters);
      outputter_prn = new Outputter::TimePrn(*this, ac_ic_print_parameters);
    }
    else if (freq_print_parameters.format_ == Format::CSV) {
      outputter_fd = new Outputter::FrequencyCSV(*this, freq_print_parameters);
      outputter_prn = new Outputter::TimeCSV(*this, ac_ic_print_parameters);
    }
    else if (freq_print_parameters.format_ == Format::PROBE) {
      outputter_fd = new Outputter::FrequencyProbe(*this, freq_print_parameters);
      outputter_prn = new Outputter::TimeProbe(*this, ac_ic_print_parameters);
    }
    else if (freq_print_parameters.format_ == Format::TECPLOT) {
      outputter_fd = new Outputter::FrequencyTecPlot(*this, freq_print_parameters);
      outputter_prn = new Outputter::TimeTecPlot(*this, ac_ic_print_parameters);
    }
    else if (freq_print_parameters.format_ == Format::RAW) {
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
      string msg("AC output can only use tecplot or gnuplot format.");
      msg += "  Resetting to gnuplot format.\n";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_WARNING, msg);

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
// Function      : OutputMgr::enableACoutput
// Purpose       : turns on AC output
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/10/00
//-----------------------------------------------------------------------------
void OutputMgr::enableTransientOutput()
{
  // prepare output manager to write
  if (tranFlag_)
  {
    std::string msg("ACfile already initialized.  Contents may be overwritten.\n");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_WARNING_0, msg);
  }

  else
  {
    tranFlag_ = true;

    PrintParameters transient_print_parameters = outputParameterMap_[OutputType::TRAN];

    if (transient_print_parameters.format_ != Format::PROBE)
      transient_print_parameters.variableList_.push_front(N_UTL_Param("TIME", 0.0));
    if (transient_print_parameters.index_)
      transient_print_parameters.variableList_.push_front(N_UTL_Param("INDEX", 0.0));

    fixupPrintParameters(DOMAIN_TIME, transient_print_parameters);

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
      string msg("AC output can only use tecplot or gnuplot format.");
      msg += "  Resetting to gnuplot format.\n";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_WARNING, msg);

      outputter = new Outputter::TimePrn(*this, transient_print_parameters);
    }

    outputter->parse();
    outputterMap_[PrintType::TRAN] = outputter;
    addActiveOutputter(PrintType::TRAN);
  }
}

void OutputMgr::enableMPDEOutput()
{
  // prepare output manager to write
  if (MPDEFlag_)
  {
    std::string msg("ACfile already initialized.  Contents may be overwritten.\n");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_WARNING_0, msg);
  }

  else
  {
    MPDEFlag_ = true;

    PrintParameters mpde_print_parameters = outputParameterMap_[OutputType::TRAN];
    PrintParameters mpde_ic_print_parameters = outputParameterMap_[OutputType::TRAN];

    mpde_ic_print_parameters.extension_ = ".mpde_ic.prn";

    if (mpde_ic_print_parameters.format_ != Format::PROBE)
      mpde_ic_print_parameters.variableList_.push_front(N_UTL_Param("TIME", 0.0));
    if (mpde_ic_print_parameters.index_)
      mpde_ic_print_parameters.variableList_.push_front(N_UTL_Param("INDEX", 0.0));

    fixupPrintParameters(DOMAIN_TIME, mpde_print_parameters);
    fixupPrintParameters(DOMAIN_TIME, mpde_ic_print_parameters);

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
      string msg("MPDE output can only use tecplot or gnuplot format.");
      msg += "  Resetting to gnuplot format.\n";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_WARNING, msg);

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


void OutputMgr::enableDCSweepOutput()
{
  // prepare output manager to write
  if (DCSweepFlag_)
  {
    std::string msg("ACfile already initialized.  Contents may be overwritten.\n");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_WARNING_0, msg);
  }

  else
  {
    DCSweepFlag_ = true;

    PrintParameters dc_print_parameters = outputParameterMap_[OutputType::DC];

    if (dc_print_parameters.index_)
      dc_print_parameters.variableList_.push_front(N_UTL_Param("INDEX", 0.0));

    fixupPrintParameters(DOMAIN_TIME, dc_print_parameters);

    Outputter::Interface *outputter_prn;
    if (dc_print_parameters.format_ == Format::STD) {
      dc_print_parameters.extension_ = ".prn";
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
      string msg("AC output can only use tecplot or gnuplot format.");
      msg += "  Resetting to gnuplot format.\n";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_WARNING, msg);

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
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/10/00
//-----------------------------------------------------------------------------
void OutputMgr::enableHBOutput()
{
  // prepare output manager to write
  if (HBFlag_)
  {
    std::string msg("HBfile already initialized.  Contents may be overwritten.\n");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_WARNING_0, msg);
  }

  else
  {
    HBFlag_ = true;

    PrintParameters freq_print_parameters = outputParameterMap_[OutputType::HB_FD];
    if (freq_print_parameters.format_ != Format::PROBE)
      freq_print_parameters.variableList_.push_front(N_UTL_Param("FREQUENCY", 0.0));
    if (freq_print_parameters.index_)
      freq_print_parameters.variableList_.push_front(N_UTL_Param("INDEX", 0.0));

    PrintParameters time_print_parameters = outputParameterMap_[OutputType::HB_TD];
    if (time_print_parameters.format_ != Format::PROBE)
      time_print_parameters.variableList_.push_front(N_UTL_Param("TIME", 0.0));
    if (time_print_parameters.index_)
      time_print_parameters.variableList_.push_front(N_UTL_Param("INDEX", 0.0));

    PrintParameters hb_ic_print_parameters = outputParameterMap_[OutputType::HB_IC];
    if (hb_ic_print_parameters.format_ != Format::PROBE)
      hb_ic_print_parameters.variableList_.push_front(N_UTL_Param("TIME", 0.0));
    if (hb_ic_print_parameters.index_)
      hb_ic_print_parameters.variableList_.push_front(N_UTL_Param("INDEX", 0.0));

    PrintParameters hb_startup_print_parameters = outputParameterMap_[OutputType::HB_STARTUP];
    if (hb_startup_print_parameters.format_ != Format::PROBE)
      hb_startup_print_parameters.variableList_.push_front(N_UTL_Param("TIME", 0.0));
    if (hb_startup_print_parameters.index_)
      hb_startup_print_parameters.variableList_.push_front(N_UTL_Param("INDEX", 0.0));

    fixupPrintParameters(DOMAIN_FREQUENCY, freq_print_parameters);
    fixupPrintParameters(DOMAIN_TIME, time_print_parameters);
    fixupPrintParameters(DOMAIN_TIME, hb_ic_print_parameters);
    fixupPrintParameters(DOMAIN_TIME, hb_startup_print_parameters);

    Outputter::Interface *outputter_hb;
    Outputter::Interface *outputter_init;
    Outputter::Interface *outputter_startup;
    if (freq_print_parameters.format_ == Format::STD) {
      hb_ic_print_parameters.extension_ = ".hb_ic.prn";
      hb_startup_print_parameters.extension_ = ".startup.prn";
      outputter_hb = new Outputter::HBPrn(*this, freq_print_parameters, time_print_parameters);
      outputter_init = new Outputter::TimePrn(*this, hb_ic_print_parameters);
      outputter_startup = new Outputter::TimePrn(*this, hb_startup_print_parameters);
    }
    else if (freq_print_parameters.format_ == Format::CSV) {
      hb_ic_print_parameters.extension_ = ".hb_ic.csv";
      hb_startup_print_parameters.extension_ = ".startup.csv";
      outputter_hb = new Outputter::HBCSV(*this, freq_print_parameters, time_print_parameters);
      outputter_init = new Outputter::TimeCSV(*this, hb_ic_print_parameters);
      outputter_startup = new Outputter::TimeCSV(*this, hb_startup_print_parameters);
    }
    else if (freq_print_parameters.format_ == Format::TECPLOT) {
      hb_ic_print_parameters.extension_ = ".hb_ic.dat";
      hb_startup_print_parameters.extension_ = ".startup.dat";
      outputter_hb = new Outputter::HBTecPlot(*this, freq_print_parameters, time_print_parameters);
      outputter_init = new Outputter::TimeTecPlot(*this, hb_ic_print_parameters);
      outputter_startup = new Outputter::TimeTecPlot(*this, hb_startup_print_parameters);
    }
    else
    {
      string msg("HB output can only use tecplot or gnuplot format.");
      msg += "  Resetting to gnuplot format.\n";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_WARNING, msg);

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
// Function      : OutputMgr::enableRawOutput
// Purpose       : turns on RAW output
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/10/00
//-----------------------------------------------------------------------------
void OutputMgr::enableOverrideRawOutput(const PrintParameters &print_parameters)
{
  // prepare output manager to write
  if (rawFlag_)
  {
    std::string msg("Rawfile already initialized.  Contents may be overwritten.\n");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_WARNING_0, msg);
  }

  else
  {
    rawFlag_ = true;

    PrintParameters raw_print_parameters = print_parameters;
    raw_print_parameters.filename_ = commandLine_.getArgumentValue("-r");

    fixupPrintParameters(DOMAIN_FREQUENCY, raw_print_parameters);

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
    //topPtr_->getRawData(nodeRef, typeRef);

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

    globalMaxNodeNameLength+1;  // add one for string terminator
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
      std::cout << "Error in making solVarMapGroup on procID_ " << getProcID() << std::endl;

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
      std::cout << "Error closing fixedLenString. on " << getProcID() << std::endl;
    status = H5Dclose(nodeNameDataSet);
    if (status < 0)
      std::cout << "Error closing nodeNameDataSet. on " << getProcID() << std::endl;
    status = H5Sclose(fileDataspace);
    if (status < 0)
      std::cout << "Error closing solVarMapDataspace. on " << getProcID() << std::endl;
    status = H5Gclose(solVarMapGroup);
    if (status < 0)
      std::cout << "Error closing solVarMapGroup. on " << getProcID() << std::endl;

    // now output the independent variables
    hid_t independentVarGroup = H5Gcreate(hdf5FileId_, "IndepVars", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (independentVarGroup < 0)
      std::cout << "Error in making independentVarGroup on getProcID() " << getProcID() << std::endl;

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
      std::cout << "Error in making solutionVecGroup on getProcID() " << getProcID() << std::endl;

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
      std::cout << "Error closing hdf5FileId_. on " << getProcID() << std::endl;
#endif  // Xyce_USE_HDF5

  return result;
}


void
removeStarVariables(ParameterList &variable_list, N_PDS_Comm &communicator, const NodeNamePairMap &allNodes_, const NodeNamePairMap &externNodes_)
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

    for (NodeNamePairMap::const_iterator it = externNodes_.begin(); it != externNodes_.end() ; ++it)
    {
      ExtendedString tmpStr((*it).first);
      tmpStr.toUpper();
// DEBUG_ELR:  do all curr end in "branch"?
      std::string::size_type pos = tmpStr.rfind("BRANCH");

      if (pos == std::string::npos && vStarFound)
      {
        vStarList.push_back(N_UTL_Param("V", 1.0));
        vStarList.push_back(N_UTL_Param(tmpStr, 0.0));

        if (!communicator.isSerial())
        {
          sB += tmpStr.size();
          sB += sizeof(int);
        }
      }
    }

    ParameterList iStarList;

    for (NodeNamePairMap::const_iterator iter_a = allNodes_.begin(); iter_a != allNodes_.end() ; ++iter_a)
    {
      ExtendedString tmpStr((*iter_a).first);
      tmpStr.toUpper();
      size_t pos = tmpStr.rfind("BRANCH");

// DEBUG_ELR: drop Y*branch; broken
//        else if (pos != std::string::npos && iStarFound)
      if (pos != std::string::npos && iStarFound && tmpStr[0] != 'Y')
      {
        tmpStr = tmpStr.substr(0, pos - 1);

        iStarList.push_back(N_UTL_Param("I", 1.0));
        iStarList.push_back(N_UTL_Param(tmpStr, 0.0));

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

            N_UTL_Param param;
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

            N_UTL_Param param;
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

          N_UTL_Param param;
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

          N_UTL_Param param;
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


void OutputMgr::fixupPrintParameters(Domain domain, PrintParameters &print_parameters, bool expandComplexTypes)
{
  fixupExpressions(print_parameters.variableList_.begin(), print_parameters.variableList_.end(), *this);

  topPtr_->getNodeNames(allNodes_);
  topPtr_->getStateNodeNames(stateNodes_);
  topPtr_->getStoreNodeNames(storeNodes_);
  topPtr_->getExternNodeNames(externNodes_);

  removeStarVariables(print_parameters.variableList_, *pdsCommPtr_, allNodes_, externNodes_);

  // setup hdf5 output if requeted
  if (hdf5FileNameGiven_)
  {
#ifdef Xyce_USE_HDF5
    enableHDF5Output();
#endif // Xyce_USE_HDF5
  }

  for (ParameterList::iterator it = print_parameters.variableList_.begin(); it != print_parameters.variableList_.end(); )
  {
    if ((*it).getSimContext() == UNDEFINED)
    {
      bool result = setParamContextType_(it);
      if (!result)
      {
        std::string msg("Can't find context for print variable ");
        msg += (*it).tag() + " " + (*it).tag() + " " + (*it).tag();
        // msg is not quite complete at this stage.  it could be that it has two tags
        // need better error report in this case.
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL_0, msg);

      }
    }

    // If frequency domain and time domain 'V' specifies copy and make this one real and next one imaginary
    if (expandComplexTypes && ((*it).getSimContext() == SOLUTION_VAR && domain == DOMAIN_FREQUENCY)) {
      ParameterList::iterator next_it = it;
      ++next_it;

      print_parameters.variableList_.insert(next_it, (*it));
      (*it).setTag(std::string("Re(") + (*it).tag() + ")");
      (*it).setSimContextAndData(SOLUTION_VAR_REAL);
      ++it;
      (*it).setTag(std::string("Im(") + (*it).tag() + ")");
      (*it).setSimContextAndData(SOLUTION_VAR_IMAG);
    }

    // If frequency domain and time domain 'V' specifies copy and make this one real and next one imaginary
    if (expandComplexTypes && ((*it).getSimContext() == VOLTAGE_DIFFERENCE && domain == DOMAIN_FREQUENCY)) {
      ParameterList::iterator next_it = it;
      ++next_it;

      print_parameters.variableList_.insert(next_it, (*it));
      (*it).setTag(std::string("Re(") + (*it).tag() + ")");
      (*it).setSimContextAndData(VOLTAGE_DIFFERENCE_REAL);
      ++it;
      (*it).setTag(std::string("Im(") + (*it).tag() + ")");
      (*it).setSimContextAndData(VOLTAGE_DIFFERENCE_IMAG);
    }

    // in seting the context of the PRINTblock_ params, some of them will no
    // longer be needed. for example V(a) is stored in two parameters one for "V"
    // and the second for "a".  After setParamContextType_ is called, the first
    // param will have all the data needed to evalueate V(a) and the second one will
    // be labeled as NODE_OR_DEVICE_NAME.  If it is, then we can safely erase it
    // from the PRINTblock_.params
    if ((*it).getSimContext() == NODE_OR_DEVICE_NAME)
    {
      it = print_parameters.variableList_.erase(it);
    }
    else
    {
      ++it;
    }
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
  bool                                  skipPrintLineOutput)
{
#ifdef Xyce_PARALLEL_MPI
  N_ERH_ErrorMgr::startSafeBarrier();
#endif

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

  if (!anaIntPtr_->getBlockAnalysisFlag()) {
    // This error test should not be used in the MPDE case, as at least
    // one of the initial conditions that can be
    // used by MPDE is a DC sweep.
    if (maxDCSteps_ > 0 && PRINTType_ == PrintType::TRAN)
    {
      std::string msg("Print type is inconsistent.");
      msg += " maxDCSteps = ";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg, maxDCSteps_);
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
  circuitTemp_ = devPtr_->getParam("TEMP");

  // Needs to pass skipPrintLineOutput
  if (!skipPrintLineOutput)
  {
    if (!activeOutputterStack_.empty())
      for (std::vector<Outputter::Interface *>::const_iterator it = activeOutputterStack_.back().begin(); it != activeOutputterStack_.back().end(); ++it)
        (*it)->output(solnVecPtr, stateVecPtr, storeVecPtr);
  }

#ifdef Xyce_USE_HDF5
  if (hdf5FileNameGiven_)
  {
    updateHDF5Output(solnVecPtr);
  }
#endif // Xyce_USE_HDF5

  // call on the measure manager to update any active measures
  RCP< N_LAS_Vector > solnVecRCP(solnVecPtr, false);
  if (PRINTType_ == PrintType::TRAN)
  {
    measureManager_.updateTranMeasures(circuitTime_, solnVecRCP);
    fourierManager_.updateFourierData(circuitTime_, solnVecRCP);
  }
  else
  {
    measureManager_.updateDcMeasures(dcParamVec_, solnVecRCP);
  }

  // if any variables are needed for a response function
  // save them now.
  if (numResponseVars_ != 0)
  {
    saveResponseVarValues(solnVecPtr);
  }

  // transient or dc values for the objective function call.
  double arg1 = 0.0;
  double arg2 = 0.0;
  if (PRINTType_ == PrintType::TRAN || anaIntPtr_->getAnalysisMgr()->getTransientFlag())
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
  for (ObjectiveMap::iterator ob = objective.begin(); ob != objective.end(); ++ob)
    (*ob).second.save(arg1, arg2, solnVecPtr, stateVecPtr, storeVecPtr);

  outputCalledBefore_ = true;

#ifdef Xyce_PARALLEL_MPI
  N_ERH_ErrorMgr::safeBarrier(0);
#endif
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::printDeviceCounts
// Purpose       : Print device summary
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 08/15/06
//-----------------------------------------------------------------------------
void OutputMgr::printDeviceCounts()
{
  std::string msg;
  std::map<std::string,int> globalDeviceCountMap;

  if (localDeviceCountMap_.empty())
  {
    // get device count from device package.
    const std::map<std::string, int> & localDeviceCount = devPtr_->getDeviceCountMap();
    gatherGlobalDeviceCount_(globalDeviceCountMap , localDeviceCount);
    formatPrintDeviceCount_ (globalDeviceCountMap, msg);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_OUT_0, msg);

    if (detailedDeviceCountFlag_)
    {
      outputLocalDeviceCount_(localDeviceCount);
    }
  }
  else
  {
    gatherGlobalDeviceCount_(globalDeviceCountMap , localDeviceCountMap_);
    formatPrintDeviceCount_ (globalDeviceCountMap, msg);
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_OUT_0, msg);

    if (detailedDeviceCountFlag_)
    {
      outputLocalDeviceCount_(localDeviceCountMap_);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::gatherGlobalDeviceCount_
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
void OutputMgr::gatherGlobalDeviceCount_(
      std::map<std::string, int> & globalDeviceCount ,
      const std::map<std::string, int> & localDeviceCount)
{

  if (!pdsCommPtr_->isSerial())
  {
    std::map<std::string, int>::const_iterator dc;
    std::map<std::string, int>::const_iterator dc_end = localDeviceCount.end();

    int i, len = 0;
    int procID = pdsCommPtr_->procID();
    int numProc = pdsCommPtr_->numProc();
    int lowestKnown;
    int numDev, numDevTot;
    std::set<std::string> known;
    std::string curr;

    lowestKnown = numProc;
    if (!localDeviceCount.empty())
    {
      lowestKnown = procID;
    }
    i = lowestKnown;
    pdsCommPtr_->minAll(&i, &lowestKnown, 1);
    dc = localDeviceCount.begin();
    while (lowestKnown < numProc)
    {
      if (lowestKnown == procID)
      {
        curr =(*dc).first;
        len = curr.size();
      }
      pdsCommPtr_->bcast(&len, 1, lowestKnown);
      curr.resize(len);
      pdsCommPtr_->bcast(&curr[0], len, lowestKnown);
      dc = localDeviceCount.find(curr);
      numDev = 0;
      if (dc != dc_end)
      {
        numDev =(*dc).second;
      }
      known.insert(curr);
      pdsCommPtr_->sumAll(&numDev, &numDevTot, 1);
      globalDeviceCount[curr] = numDevTot;
      lowestKnown = numProc;
      dc = localDeviceCount.begin();
      for (; dc != dc_end ; ++dc)
      {
        if (known.find((*dc).first) == known.end())
        {
          lowestKnown = procID;
          curr =(*dc).first;
          break;
        }
      }
      i = lowestKnown;
      pdsCommPtr_->minAll(&i, &lowestKnown, 1);
    }
  }
  else
  {
    globalDeviceCount = localDeviceCount;
  }
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::formatPrintDeviceCount_
//
// Purpose       : takes the(passed) device count map, and formats a string
//                 that can be output, either via the error handler or std::cout.
//
// Special Notes :
//
// Scope         : private
// Creator       : Eric R. Keiter, SNL
// Creation Date : 05/10/10
//-----------------------------------------------------------------------------
void OutputMgr::formatPrintDeviceCount_(const std::map<std::string, int> & deviceCountMap_, std::string & msg)
{
  int i, len = 0, maxLen = 15, totDev = 0, width = 0;
  std::map<std::string, int>::const_iterator dc = deviceCountMap_.begin();
  std::map<std::string, int>::const_iterator dc_end = deviceCountMap_.end();
  msg = "";

  for (; dc != dc_end ; ++dc)
  {
    len =(*dc).first.size();
    if (len > maxLen)
      maxLen = len;
    totDev +=(*dc).second;
  }
  i = totDev;
  while (i > 0)
  {
    width++;
    i /= 10;
  }

  dc = deviceCountMap_.begin();
  for (; dc != dc_end ; ++dc)
  {
    len =(*dc).first.size();
    msg += "       " +(*dc).first;
    for (i=0 ; i<maxLen-len+1 ; ++i)
      msg += " ";
  std::ostringstream ost;
    ost.width(width);
    ost.setf(ios::right);
    ost <<(*dc).second;
    msg += ost.str();
    msg += "\n";
  }
  msg += "       ";
  for (i=0 ; i<maxLen+width+1 ; ++i)
    msg += "-";
  msg += "\n       Total Devices";
  for (i=0 ; i<maxLen-12 ; ++i)
  {
    msg += " ";
  }
  std::ostringstream ost;
  ost.width(width);
  ost.setf(ios::right);
  ost << totDev;
  msg += ost.str();
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::outputLocalDeviceCount_
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
void OutputMgr::outputLocalDeviceCount_(const std::map<std::string, int> & localDeviceCount)
{
  int masterRank = 0;
  int numProcs = pdsCommPtr_->numProc();
  int thisProc = pdsCommPtr_->procID();

  pdsCommPtr_->barrier();

  for (int p = 0; p < numProcs; ++p)
  {
    pdsCommPtr_->barrier();
    if (p==thisProc)
    {
      std::string msg1="\n\tDevice Count for proc=";
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_OUT, msg1, p);

      std::string msg="";
      formatPrintDeviceCount_(localDeviceCount , msg);
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_OUT, msg);
    }
  }
  pdsCommPtr_->barrier();

}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::addDeviceToCount
// Purpose       : Add a device for device count report
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 08/15/06
//-----------------------------------------------------------------------------
void OutputMgr::addDeviceToCount(std::string & devNameIn)
{
  localDeviceCountMap_[devNameIn]++;

  return;
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
  ofstream opOut;

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
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL,
        std::string("Internal error 1 in OutputMgr::outputDCOP"));
    }
    int flag;
    pdsCommPtr_->recv(&flag, 1, 0);
    if (flag != 1)
    {
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL,
        std::string("Internal error 2 in OutputMgr::outputDCOP"));
    }
    pdsCommPtr_->send(&mySize, 1, 0);
    pdsCommPtr_->recv(&flag, 1, 0);
    if (flag != 2)
    {
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL,
        std::string("Internal error 3 in OutputMgr::outputDCOP"));
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
  ofstream saveOutputStream;

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
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL,
        std::string("Internal error 1 in OutputMgr::outputIC_or_NODESET"));
    }
    int flag;
    pdsCommPtr_->recv(&flag, 1, 0);
    if (flag != 1)
    {
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL,
        std::string("Internal error 2 in OutputMgr::outputIC_or_NODESET"));
    }
    pdsCommPtr_->send(&mySize, 1, 0);
    pdsCommPtr_->recv(&flag, 1, 0);
    if (flag != 2)
    {
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL,
        std::string("Internal error 3 in OutputMgr::outputIC_or_NODESET"));
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
  ifstream opIn;
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
    opIn.open(inputFileName.c_str(), ios::in);

    if (opIn.good())
    {
      cout << "Reading in operating point initial estimate from: " << inputFileName << std::endl;
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
//  cout << "PE: " << getProcID() << " has " << numBad << " mismatched nodes" << std::endl;
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
        cout << "DCOP restart:  Initialized " << nodesMatched << " of a possible " << totalNodes << " nodes" << std::endl;
      else
        cout << "DCOP restart:  All " << totalNodes << " nodes initialized" << std::endl;

      if (notFound.size() > 0)
      {
        cout << "DCOP restart:  Nodes specified in " << inputFileName << " but not present in circuit:" << std::endl;
        nm_i = notFound.begin();
        nm_end = notFound.end();
        for (; nm_i != nm_end ; ++nm_i)
          cout << *nm_i << std::endl;
      }
      if (notMatched.size() > 0)
      {
        cout << "DCOP restart:  Nodes present in circuit, but not specified in " << inputFileName << ":" << std::endl;
        nm_i = notMatched.begin();
        nm_end = notMatched.end();
        for (; nm_i != nm_end ; ++nm_i)
          cout << *nm_i << std::endl;
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
// Function      : OutputMgr::getAllNodes
// Purpose       : provides a map of node names.  Needed by a variety of
//                 things, including LOCA.
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 05/25/06
//-----------------------------------------------------------------------------
NodeNamePairMap & OutputMgr::getAllNodes()
{
  return allNodes_;
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
                                std::vector<N_UTL_OptionBlock> & initBlockVec)
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
      // and there is an error trap for them in the N_IO_OptionBlock.C file.
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
        std::string msg("Problems processing" + icType +  " values");
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL_0, msg);
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
        cout.width(10);cout.precision(3);cout.setf(ios::scientific);
        cout
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
        cout << "procID="<<getProcID() << "  globalFound = true" << std::endl;
      else
        cout << "procID="<<getProcID() << "  globalFound = false" << std::endl;
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
      cout << icType << ":  Initialized " << matched.size() << " of a possible " << totalNodes << " nodes" << std::endl;
    }
    else
    {
      cout << icType << ":  All " << totalNodes << " nodes initialized" << std::endl;
    }
#endif

    if (notFoundInCkt.size() > 0)
    {
      std::string msg1(icType + ":  Nodes specified on ." + icType +" line, but not present in circuit. (ignoring):");
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_OUT_0, msg1);

      nm_i = notFoundInCkt.begin();
      nm_end = notFoundInCkt.end();
      for (; nm_i != nm_end ; ++nm_i)
      {
        std::string msg2(*nm_i);
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_OUT_0, msg2);
      }
    }

#ifdef Xyce_DEBUG_IO
    // Note: for a typical .IC line, this list of unspecified nodes will be
    // a very long list.  So, don't output except in debug mode.
    if (notSpecified.size() > 0)
    {
      std::string msg1(icType + ":  Nodes present in circuit, but not specified on ." + icType + " line(ignoring):");
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_OUT_0, msg1);

      nm_i = notSpecified.begin();
      nm_end = notSpecified.end();
      for (; nm_i != nm_end ; ++nm_i)
      {
        std::string msg2(*nm_i);
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_OUT_0, msg2);
      }
    }
#endif // debug io

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
      std::string msg("Cannot set both DCOP restart and .IC simultaneously.");
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL_0, msg);
    }

    if (input_op_ && nodesetflag_)
    {
      std::string msg("Cannot set both DCOP restart and .NODESET simultaneously.");
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL_0, msg);
    }

    if (ICflag_ && nodesetflag_)
    {
      std::string msg("Cannot set both .IC and .NODESET simultaneously.");
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL_0, msg);
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

  ofstream dakOut;
  dakOut.open(std::string(responseFileName_ + ".result").c_str());

  dakOut << "Dakota optimization results:" << std::endl;
  for (i=0 ; i<num ; ++i)
  {
    if (!done[i])
    {
      std::string::size_type cpos = responseNames_[i].find_first_of(':');
      if (cpos == std::string::npos)
      {
        devPtr_->getParam(responseNames_[i], val);
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
            devPtr_->getParam(responseNames_[j], val);
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
#endif

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
// Function      : OutputMgr::setParamContextType_
// Purpose       : This routine looks at the N_UTL_Param object passed in
//                 and sets its context type.  This makes looking up the
//                 param's value faster as once it is done, no further
//                 context searches are needed.  Much of this code originally
//                 resided in getPrintValue_() but is now here to streamline
//                 the lookup process.
// Special Notes :
// Scope         : private
// Creator       : Rich Schiek, Electrical Systems Modeling, SNL
// Creation Date : 11/6/12
//-----------------------------------------------------------------------------
bool OutputMgr::setParamContextType_(ParameterList::iterator originalParamItr)
{
  bool foundAndSetContext=false;
  int index1 = -1;
  int index2 = -1;
  // rather than change the status of the passed in originalParamItr,
  // make a copy for advancement called nextParamItr.
  // in some cases we need to look ahead to the next one or two
  // iterators in the ParameterList objects.  Once all the nodes have
  // their context set, we can remove the extra nodes that were used for
  // holding part of the context.
  ParameterList::iterator nextParamItr = originalParamItr;

  std::ostringstream ost;
  ost << getProcID();
  string pN(ost.str());

#ifdef Xyce_DEBUG_IO
    cout << " Proc " <<  pN  << ": In setParamContextType_ with:";
    cout << originalParamItr->tag() << " : " << originalParamItr->sVal() << std::endl;
#endif

  // varType should be I, V, N, B, H as in I(name), V(name), or V(name, name)
  string varType(originalParamItr->tag());
  char type1, type2;

  double result=0;
#ifdef Xyce_PARALLEL_MPI
  double final=0;
#endif

  bool resultFound = false;

  // If this is "TEMP", we know where to find it.
  if (varType == "TEMP")
  {
    originalParamItr->setSimContextAndData(TEMPERATURE);
    foundAndSetContext=true;
  }

  // if this is TIME.  We know where to find it.
  if (varType == "TIME")
  {
    originalParamItr->setSimContextAndData(TIME_VAR);
    foundAndSetContext=true;
  }

  // if this is FREQ.  We know where to find it.
  if (varType == "FREQUENCY" || varType == "FREQUENCY")
  {
    originalParamItr->setSimContextAndData(FREQUENCY);
    foundAndSetContext=true;
  }

  // if this is INDEX.  We know where to find it.
  if (varType == "INDEX")
  {
    originalParamItr->setSimContextAndData(INDEX);
    foundAndSetContext=true;
  }

  // check if this is a step sweep value.
  if (!foundAndSetContext)
  {
    int numParams = stepParamVec_.size();
    for (int i=0; i<numParams; i++)
    {
      if (varType == stepParamVec_.at(i).name)
      {
        originalParamItr->setSimContextAndData(STEP_SWEEP_VAR, i);
        foundAndSetContext=true;
        break;
      }
    }
  }

  // check if this is a dc sweep value.
  if (!foundAndSetContext)
  {
    int numParams = dcParamVec_.size();
    for (int i=0; i<numParams; i++)
    {
      if (varType == dcParamVec_.at(i).name)
      {
        originalParamItr->setSimContextAndData(DC_SWEEP_VAR, i);
        foundAndSetContext=true;
        break;
      }
    }
  }

  // Last easy check. Global params
  if (!foundAndSetContext)
  {
    if (varType == "GLOBAL_PARAMETER")
    {
      originalParamItr->setSimContextAndData(GLOBAL_PARAMETER);
      foundAndSetContext=true;
    }
  }

#ifdef Xyce_DEBUG_IO
  cout << "\nsetParamContextType_: varType = " << varType;
  if (!foundAndSetContext)
    cout << "\tIS NOT a Sweep variable\n";
  else
    cout << "\tIS a Sweep variable\n";
  cout << std::endl;
#endif

  // erkeite note:  12/16/2007.  This needs to be handled a lot more cleanly.
  bool thisIsAnExpression = false;
  bool thisIsA_V_or_I_Var = false;
  bool thisIsProbablyALeadCurrent = false;
  bool thisIsANodeVar = false;
  bool iteraterAdvanced = false;

  if (!foundAndSetContext)
  {
    // now check for expressions and/or solution variables:
    //
    // If we get to this point, then this is not a sweep variable.
    // If not a sweep value, other valid print variables
    // include '{', 'I', 'V' and 'N'.  The current variables(for
    // BJT lead currents, etc) can be 'IB', 'IC' etc.
    // and for AC, we also have VR, VI, VM and VDB
    thisIsAnExpression = originalParamItr->hasExpressionTag();

    if (!thisIsAnExpression && originalParamItr->iVal() > 0 &&
       ((varType[0] == 'I') ||(varType[0] == 'V') ) )
    {
      thisIsA_V_or_I_Var = true;
    }

    if (!thisIsAnExpression && originalParamItr->iVal() > 0 && varType=="N")
    {
      thisIsANodeVar = true;
    }
  }

#ifdef Xyce_DEBUG_IO
  cout << "setParamContextType_: varType = " << varType;
  if (thisIsAnExpression)
    cout << "\tIS an expression.";
  else
    cout << "\tIS NOT an expression.";
  cout << std::endl;
  cout << "setParamContextType_: varType = " << varType;
  if (thisIsA_V_or_I_Var)
    cout << "\tIS a voltage or current variable(either solution or lead current)";
  else
    cout << "\tIS NOT a voltage or current variable(either solution or lead current)";
  cout << std::endl;
  if (thisIsANodeVar)
    cout << "\tIS a node variable";
  else
    cout << "\tIS NOT a node variable";
  cout << std::endl;
#endif

  // If the print value is an expression, process here.
  if (!foundAndSetContext && thisIsAnExpression)
  {
#ifdef Xyce_DEBUG_IO
    cout << "\nsetParamContextType_: processing expression" << std::endl;
#endif
    // This context is set in registerPRINTSet() because at that time expressions are flesed out
    std::cout << "OutputMgr::setParamContextType_ In expression block with unresolved context!" << std::endl;
    foundAndSetContext=true;
  }

  // One wrinkle in this is that the list of N_UTL_Param objects stores
  // I(name), V(name) and N(name) as two sequential N_UTL_Param ojects one for the
  // I, V or N and the following one for the "name".
  // Likewise V(name1, name2) uses three N_UTL_Param objects, one for V
  // and two for name1 and name2.
  // These following N_UTL_Param objects which hold names are needed once to
  // find indicies but not after that.  Thus there is an enum NODE_OR_DEVICE_NAME
  // to categorize those N_UTL_Param objects that can be dropped after context
  // resolution(or at least ignored when later processing things in getPrintValue()
  string nodeName("");    // node name of originalParamItr
  string nodeName2("");                        // node name of 2nd paramItr(needed for V(a, b))

  // Advance the parameter iterator and set up the node name(s).
  int numIndices = 0;
  if (!foundAndSetContext &&(thisIsA_V_or_I_Var || thisIsANodeVar))
  {
    // calling this only makes sense for non expressions so
    // wait to call it until we're sure we don't have an expression.
    string printVarName;
    numIndices = originalParamItr->iVal();
    ++nextParamItr;
    nodeName = nextParamItr->tag();
    // we know the context of the nextParamItr at this point so set it.
    nextParamItr->setSimContextAndData(NODE_OR_DEVICE_NAME);
    if (numIndices==2)
    {
      ++nextParamItr;
      nodeName2 = nextParamItr->tag();
      // we know the context of the nextParamItr at this point so set it.
      nextParamItr->setSimContextAndData(NODE_OR_DEVICE_NAME);
      printVarName =originalParamItr->tag() + "(" + nodeName + "," + nodeName2 + ")";
    }
    else
    {
      printVarName =originalParamItr->tag() + "(" + nodeName + ")";
    }
    originalParamItr->setTag(printVarName);

    iteraterAdvanced = true;
  }

  // If this is a node value('N'), attempt determine the value.
  // Note, a N(node) value is the same as a V(node) value, it is just less
  // error-checked, and thus can refer to internal device variables, not just
  // user-specified voltage nodes.  N() can also access values in the state
  // and store vectors.
  // if (!foundAndSetContext &&(thisIsANodeVar || thisIsA_V_or_I_Var)) -- not a good check here
  if (!foundAndSetContext && thisIsANodeVar)
  {
    NodeNamePairMap::iterator iterCI = allNodes_.find(nodeName);
    if (iterCI != allNodes_.end())
    {
      // object is part of a solution var
      int ind = iterCI->second.first;
      originalParamItr->setSimContextAndData(SOLUTION_VAR, ind);
      foundAndSetContext = true;
      //std::cout << " found context as sol-var at index " << ind << std::endl;
    }
    else
    {
      // check state vars.
      iterCI = stateNodes_.find(nodeName);
      if (iterCI != stateNodes_.end())
      {
        int ind = iterCI->second.first;
        originalParamItr->setSimContextAndData(STATE_VAR, ind);
        foundAndSetContext = true;
      }
      else
      {
        // check store vars
        iterCI = storeNodes_.find(nodeName);
        if (iterCI != storeNodes_.end())
        {
          int ind = iterCI->second.first;
          originalParamItr->setSimContextAndData(STORE_VAR, ind);
          foundAndSetContext = true;
        }
      }
    }
  }


  // Determine value if this a 'I' or a 'V' variable.  If it is a V-variable, it is
  // definately a solution variable and not a lead current.
  if (!foundAndSetContext && thisIsA_V_or_I_Var)
  {
    list<int> svGIDList1, svGIDList2, dummyList;
    svGIDList1.clear();
    svGIDList2.clear();
    dummyList.clear();

    if (varType[0] == 'I')
    {
      std::string::size_type pos = nodeName.find_last_of(":");
      if (pos == std::string::npos)
        pos = 0;
      else
        ++pos;
      // The level 3 resistor is really just a voltage source with a zero volt difference
      // Thus, it has a real branch current that is in the solution vector.  By putting
      // nodeName[pos] =='R' in the following conditional, a later block of code will
      // be executed(just a few lines down that starts:)
      //   if ((!(svGIDList1.empty()) || nodeName == "0") && !thisIsProbablyALeadCurrent)
      // If this is a level 3 resistor, then its branch current will be found in that
      // block of code.  If it's not a level 3 resistor, then its current will not be found
      // and a later block of code for finding lead currents will be executed because
      // resultFound will still be false.  Thus, it's safe to add 'R' to this if statement.
      //
      if (nodeName[pos] == 'V' || nodeName[pos] == 'L' ||
          nodeName[pos] == 'B' || nodeName[pos] == 'E' || nodeName[pos] == 'R')
      {
        thisIsProbablyALeadCurrent = false;
      }
      else
      {
        thisIsProbablyALeadCurrent = true;
      }
    }
    else
    {
      thisIsProbablyALeadCurrent = false;
    }

#ifdef Xyce_DEBUG_IO
    cout << "\nsetParamContextType_: doing the old stuff" << std::endl;
#endif

    if ((varType == "I" && distribVsrcSet_.count(nodeName))  ||
        (varType[0] == 'V' && distribVnodeSet_.count(nodeName)))
    {
      nodeName += "__" + pN;
    }

    // Need to get the solution variables from topology with the known node type.
    if (varType == "I")
      topPtr_->getNodeSVarGIDs(NodeID(nodeName, _DNODE), svGIDList1, dummyList, type1);
    if (varType[0] == 'V')
      topPtr_->getNodeSVarGIDs(NodeID(nodeName, _VNODE), svGIDList1, dummyList, type1);

#ifdef Xyce_DEBUG_IO
    cout << "nodeName = "<< nodeName << " nodeName2 = " << nodeName2 <<std::endl;
    list<int>::iterator isvg = svGIDList1.begin();
    for (;isvg!=svGIDList1.end();++isvg)
    {
      cout << "svGID: " << *isvg << std::endl;
    }
#endif

    if (numIndices == 2)
    {
      if (distribVnodeSet_.count(nodeName))
      {
        nodeName2 += "__" + pN;
      }
      // Only "V" allows for numIndices to be larger than 1.
      topPtr_->getNodeSVarGIDs(NodeID(nodeName2, _VNODE), svGIDList2, dummyList, type2);
    }

    // this is an awkward check.  In parallel if a node is not on the local processor
    // then svGIDList1 will be empty.  If the node tag was "V" and the next nodeName was "0"
    // then the user wanted the ground node, which isn't calculated but is zero.
    // Thus, what this conditional is asking is:
    //   NOT(var found in local solution vec OR V(GROUND)) AND NOT a Lead Current
    //if ((!(svGIDList1.empty()) || nodeName == "0") && !thisIsProbablyALeadCurrent)
    // would this a better check?  looking for "0" will default to an index of -1 as will svGIDList1.empty==true.
    //    if (!thisIsProbablyALeadCurrent)
    // no because if svGIDList is empty then type is not correctly set and we can't make good decisions.

    if ((!(svGIDList1.empty()) || nodeName == "0") && !thisIsProbablyALeadCurrent)
    {
      if (varType[0] == 'V')
      {
        if (!svGIDList1.empty())
          index1 = *(svGIDList1.begin());

        if (numIndices == 1)
        {
          if (type1 != 'V')
          {
            string msg("Voltage Output requested for device node(" + nodeName + ") by PRINT statement\n");
            N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
          }
          // before we set the context need to check for additional
          // letters after V as in VR, VI, VM and VDB as they change
          // how we process the voltage.
          if (varType == "V" )
          {
            originalParamItr->setSimContextAndData(SOLUTION_VAR, index1);
            foundAndSetContext = true;
          }
          else if (varType == "VR" )
          {
            originalParamItr->setSimContextAndData(SOLUTION_VAR_REAL, index1);
            foundAndSetContext = true;
          }
          else if (varType == "VI" )
          {
            originalParamItr->setSimContextAndData(SOLUTION_VAR_IMAG, index1);
            foundAndSetContext = true;
          }
          else if (varType == "VM" )
          {
            originalParamItr->setSimContextAndData(SOLUTION_VAR_MAG, index1);
            foundAndSetContext = true;
          }
          else if (varType == "VP" )
          {
            originalParamItr->setSimContextAndData(SOLUTION_VAR_PHASE, index1);
            foundAndSetContext = true;
          }
          else if (varType == "VDB" )
          {
            originalParamItr->setSimContextAndData(SOLUTION_VAR_DB, index1);
            foundAndSetContext = true;
          }
        }
        else if (numIndices == 2)
        {
          // check to ensure that User asked for a voltage as that's the
          // only thing that takes two args
          if (type1 != 'V')
          {
            string msg("Voltage Output requested for device node(" + nodeName2 + ") by PRINT statement\n");
            N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
          }

          // user requested a voltage difference V(node1, node2)
          // report V(node1) and -V(node2) depending on what's on
          // this processor
          if (svGIDList2.size() > 0)
            index2 = *(svGIDList2.begin());

          originalParamItr->setSimContextAndData(VOLTAGE_DIFFERENCE, index1, index2);
          foundAndSetContext = true;
        }
        else
        {
          // In the serial case, we know for sure there was a request for
          // output at a non-existent voltage node.
          string msg("Output requested for non-existent node \"" + nodeName + nodeName2);
          msg += "\" by PRINT statement\n";
          N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
        }
      }
      else if (varType == "I")
      {
#ifdef Xyce_DEBUG_IO
        cout <<"Processing varType == I" << std::endl;
#endif
        char ch_type = nodeName[0];
        if (nodeName.find_last_of(":") != string::npos)
          ch_type = nodeName[nodeName.find_last_of(":")+1];
        if (type1 != 'D')
        {
          string msg("Current Output requested for voltage node(" + nodeName + ") by PRINT statement\n");
          N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
        }

        if (!svGIDList1.empty())
          index1 = *(svGIDList1.begin());

        // need this conditional as this could be I(Rmyresistor).  Only the level 3 resistor has
        // a current in the solution vector.  I'll need this check untill lead currents are properly
        // handled by the store vector
        if (index1 != -1)
        {
          originalParamItr->setSimContextAndData(SOLUTION_VAR, index1);
          foundAndSetContext = true;
        }
#ifdef Xyce_DEBUG_IO
        cout <<"index1 = " << index1 << "  foundAndSetContext = " << foundAndSetContext << std::endl;
#endif
      }
    }
    else if (svGIDList1.empty() && varType[0] == 'V' && numIndices == 2)
    {
      // This case is *extremely* unfortunate, but is a result of the
      // design of this function being very fragile with respect to parallel
      // implementation.
      //
      // When we have V(a,b) (the only case where 2 indices are legal),
      // and we're paralell, and on a processor where A is a node stored
      // off-processor, the previous catch-all case misses it, and we would
      // otherwise fall off and set this node to return "constant".  This
      // is a serious flaw in the logic and makes V(A,B) be completely wrong
      // in parallel.  It also makes V(A,B) invisible to all processors that
      // don't own "A" --- meaning that all the clumsy logic in getPrintValue
      // won't catch this case.  Unless we do this little game right here.
      //
      // Here, if any processor has V(A,B) and A is off processor, we
      // *STILL* set its context as VOLTAGE_DIFFERENCE, and fudge the
      // indices.
      //
      // Logic inside getPrintValue check these special cases and do The Right
      // Thing.
      //
      // The *entire* handling of the print line needs to be refactored with
      // parallelism in mind from the outset, and must not remain a hodge-podge
      // of serial code with parallelism band-aided in.

      index1 = -1;
      index2 = -1;
      if (!svGIDList2.empty())
        index2 = *(svGIDList2.begin());

      originalParamItr->setSimContextAndData(VOLTAGE_DIFFERENCE, index1, index2);
      foundAndSetContext = true;
    }

#ifdef Xyce_DEBUG_IO
    cout << "procID = " << getProcID() << ":";
    cout << " nodeName = " << nodeName << ":";
    cout << " nodeName2 = " << nodeName2 << ":";
    cout << " index1 = " << index1;
    cout << " result = " << result << std::endl;
#endif
  } // !resultFound &&(expression || solution var)

  bool globalFound;

#ifdef Xyce_PARALLEL_MPI
  // need to check in parallel that any node not found locally
  // was found on another processor.

  // double dblFound = resultFound?1.0:0.0;
  double dblFound = 0.0;
  if (foundAndSetContext)
    dblFound = 1.0;
  double globalDblFound = 0.0;
  pdsCommPtr_->maxAll(&dblFound, &globalDblFound, 1);
  int globalContexFound=UNDEFINED;
  int localContext = originalParamItr->getSimContext();
  pdsCommPtr_->maxAll(&localContext, &globalContexFound, 1);
  // resultFound =(globalDblFound!=0.0)?true:false;
  if (globalDblFound != 0.0)
  {
    if (foundAndSetContext==false)
    {
      // a valid context was found on another processor,
      // so set this context to CONSTANT and value of zero
      if ( static_cast<SimulatorVariableContext>(globalContexFound) == SOLUTION_VAR )
      {
        // this is a hack to handle complex types from getPrintValue to in parallel 
        // In parallel all procs need to call getPrintValue() twice to get the real
        // and imaginary parts of a solution var.  But, if the context on processors 
        // that don't own the var is set to CONSTANT then they won't know to make 
        // two calls to getPrintValue() and Xyce will hang waiting on parallel 
        // communication. The real solution to this is to have getPrintValue return
        // complex types in one call.  But that is for a later release
        originalParamItr->setSimContextAndData(static_cast<SimulatorVariableContext>(globalContexFound));
      }
      else
      {
        // a safe context assignment for all other cases.
        originalParamItr->setSimContextAndData(CONSTANT, 0.0);
      }
      foundAndSetContext=true;
    }
  }
  else
  {
    foundAndSetContext=false;
  }
#endif

  // Handle lead currents, including PDE lead currents.
  // Note; the devPtr_->getParam will work in parallel, so this can happen
  // after the maxall, above.
  if (!foundAndSetContext && thisIsA_V_or_I_Var)
  {
    string parTmp("");
    string parTmp2("");
    if (!iteraterAdvanced)
    {
      ++nextParamItr;
      nodeName = nextParamItr->tag();
    }

    if (varType == "I" ||(varType.size() == 2 && varType[0] == 'I'))
    {
      // nodeName could be circuit_context:DeviceTypeDeviceName while
      // internally it should be DeviceType:circuit_context:DeviceName.
      // generate the modified nodename here
      string modifiedName;
      std::string::size_type lastColonInName = nodeName.find_last_of(":");
      //std::cout << "==>lastColonInName = " << lastColonInName << std::endl;
      if ((lastColonInName != string::npos) &&(lastColonInName+1 < nodeName.length()))
      {
        string::iterator deviceName = nodeName.begin()+lastColonInName+1;
        string::iterator namePrefixEnd = nodeName.begin()+lastColonInName;
        modifiedName.append(deviceName, deviceName+1);
        modifiedName.append(":");
        modifiedName.append(nodeName.begin(), namePrefixEnd+1);
        modifiedName.append(deviceName+1, nodeName.end());
        //std::cout << "modifiedName = " << modifiedName <<std::endl;
      }
      else
      {
        modifiedName = nodeName;
      }
      // could be a device lead current "DEV_I" or a branch current.
      // so we don't have to duplicate solution vars(branch currents) in the
      // store vector, look for each type.
      parTmp = modifiedName + ":DEV_" + varType;  // if it is in the state/store vec.
      parTmp2 = modifiedName + "_BRANCH";         // if it is in the solution vec.
    }
    else
    {
      // check for v(name) in allNodes_
      // end up here when topology isn't completed(so svGIDList1.empty() == true)
      NodeNamePairMap::iterator iterCI = allNodes_.find(nodeName);
      if (iterCI != allNodes_.end())
      {
        // object is part of a solution var
        int ind = iterCI->second.first;
        originalParamItr->setSimContextAndData(SOLUTION_VAR, ind);
        foundAndSetContext = true;
        //std::cout << " found context as sol-var at index " << ind << std::endl;
      }

    }

    // this if block allows for spaces in YPDE names as in I1(YPDE NAME)
    // we try to find devices based on parTmp in the following blocks of code,
    // so do this modification now.
    std::string::size_type space = parTmp.find_first_of(" ");
    if (space != std::string::npos)
    {
      if (space == 4 && parTmp.substr(0, 4) == "YPDE")
      {
        parTmp.replace(4, 1, "%");
        parTmp.insert(1, "%");
      }
    }

    // search the store vector for current.
    NodeNamePairMap::iterator iterCI;
#ifdef Xyce_DEBUG_IO
    std::cout << "allNodes_: " << std::endl;
    iterCI = allNodes_.begin();
    while (iterCI != allNodes_.end())
    {
      std::cout << "allNode-map \"" << iterCI->first << "\" =(" << iterCI->second.first
        << ", " << iterCI->second.second << ")" << std::endl;
      iterCI++;
    }

    std::cout << "stateNodes_: " << std::endl;
    iterCI = stateNodes_.begin();
    while (iterCI != stateNodes_.end())
    {
      std::cout << "stateNode-map \"" << iterCI->first << "\" =(" << iterCI->second.first
      << ", " << iterCI->second.second << ")" << std::endl;
      iterCI++;
    }

    std::cout << "storeNodes_: " << std::endl;
    iterCI = storeNodes_.begin();
    while (iterCI != storeNodes_.end())
    {
      std::cout << "storeNode-map \"" << iterCI->first << "\" =(" << iterCI->second.first
      << ", " << iterCI->second.second << ")" << std::endl;
      iterCI++;
    }

    std::cout << " parTmp = \"" << parTmp << "\"" << std::endl;
    std::cout << " parTmp2 = \"" << parTmp2 << "\"" << std::endl;
#endif // Xyce_DEBUG_IO

    iterCI = storeNodes_.find(parTmp);
    if (iterCI != storeNodes_.end())
    {
      int ind = iterCI->second.first;
      originalParamItr->setSimContextAndData(STORE_VAR, ind);
      nextParamItr->setSimContextAndData(NODE_OR_DEVICE_NAME);
      foundAndSetContext = true;
#ifdef Xyce_DEBUG_IO
      std::cout << " found store var context for " << parTmp << " index " << ind << std::endl;
#endif
    }

    if (!foundAndSetContext)
    {
      iterCI = allNodes_.find(parTmp2);
      if (iterCI != allNodes_.end())
      {

        int ind = iterCI->second.first;
        originalParamItr->setSimContextAndData(SOLUTION_VAR, ind);
        nextParamItr->setSimContextAndData(NODE_OR_DEVICE_NAME);
#ifdef Xyce_DEBUG_IO
        std::cout << " found solution var context for " << parTmp2 << " index " << ind << std::endl;
#endif
        foundAndSetContext = true;
      }
    }

    // in parallel we may have found the variable in the state or store vectors
    // sync up the foundAndSetContext flag before checking with the device
    // package for a device parameter.
#ifdef Xyce_PARALLEL_MPI
    double foundParam = foundAndSetContext?1:0;
    double finalParam = 0.0;
    double globalVal;
    pdsCommPtr_->barrier();
    pdsCommPtr_->sumAll(&foundParam, &finalParam, 1);
    if (!foundAndSetContext &&(finalParam != 0.0))
    {
      // context was found on another processor.
      // it's safe to set this context to the store vec.
      // and the index to -1.
      foundAndSetContext=true;
      originalParamItr->setSimContextAndData(STORE_VAR, -1);
      nextParamItr->setSimContextAndData(NODE_OR_DEVICE_NAME);
#ifdef Xyce_DEBUG_IO
      std::cout << " found solution/state var context off proc for " << parTmp << " index -1" << std::endl;
#endif
    }
#endif

    if (!foundAndSetContext)
    {
      // this is confusing.  While the solution vector and state/store vector's
      // use maps with modified device names(specifically where the device
      // type is always first as in "D:subcircuitname:devciename" as apposed to
      // "subcircuitname:Ddevicename", the device manager does not use the modified
      // device name to find a device.  So, when we set up the device name below
      // use the "nodeName" rather than the "modifiedDeviceName".
      parTmp = nodeName + ":DEV_" + varType;
      // have to repeat this check for spaces as in I(YPDE NAME)
      std::string::size_type space = parTmp.find_first_of(" ");
      if (space != std::string::npos)
      {
        if (space == 4 && parTmp.substr(0, 4) == "YPDE")
        {
          parTmp.replace(4, 1, "%");
          parTmp.insert(1, "%");
        }
      }
      originalParamItr->setSimContextAndData(DEVICE_PARAMETER, parTmp);
      nextParamItr->setSimContextAndData(NODE_OR_DEVICE_NAME);

      // This starts to querry the device manager to test if
      // this parameter is real, or if it need further qualification

#ifdef Xyce_DEBUG_IO
      std::cout << " Checking device manager for " << parTmp << std::endl;
#endif
      foundAndSetContext = devPtr_->getParam(parTmp, result);
      if (!foundAndSetContext)
      {
        string ppde("Y%PDE%" + parTmp);
        originalParamItr->setSimContextAndData(DEVICE_PARAMETER, ppde);
        foundAndSetContext = devPtr_->getParam(ppde, result);
      }
#ifdef Xyce_DEBUG_IO
      if (foundAndSetContext)
      {
        std::cout << " found device parameter context for " << parTmp << std::endl;
      }
#endif
    }
  }
  globalFound = foundAndSetContext;

  // If the parameter is still "not found", then try to get it from
  // the device package.  I added this as an afterthought. ERK.
  if (!globalFound)
  {
    if (thisIsA_V_or_I_Var)
    {
      if (varType == "V")
      {
        string msg("VOLTAGE Output requested for non-existent node \"" + nodeName);
        msg += "\" by PRINT or Expression statement\n";
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
      }
      else if (varType == "I" ||(varType.size() == 2 && varType[0] == 'I'))
      {
        string msg("CURRENT Output requested for unsupported or non-existent device \"" + nodeName);
        msg += "\" by PRINT or Expression statement\n";
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
      }
    }
    else
    {
      // possibly it is a global param or other known param
      globalFound = devPtr_->getParam(varType, result);
      if (globalFound)
      {
        originalParamItr->setSimContextAndData(DEVICE_PARAMETER, varType);
        foundAndSetContext=true;
      }
    }
  }

  // Handle objective vars and measure vars.
  if (!globalFound && !thisIsA_V_or_I_Var)
  {
    if (objective.find(varType) != objective.end())
    {
      originalParamItr->setSimContextAndData(OBJECTIVE_FUNCTION, varType);
      foundAndSetContext=true;
    }
    else  // not an objective so potentially a measure var.
    {
      bool mFound;
      measureManager_.getMeasureValue(varType, result, mFound);
      if (!mFound) {
        string msg("Can't find print variable "+varType+" "+nodeName);
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL_0, msg);
      }
      else
      {
        originalParamItr->setSimContextAndData(MEASURE_FUNCTION, varType);
        foundAndSetContext=true;
      }
    }
  }

  return foundAndSetContext;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::getPrintValue_
// Purpose       :
// Special Notes : The design of this method is defective.
//                 We are doing an MPI reduce operation on Every. Single.
//                 Element. on the print line, one at a time.
//                 This is RIPE for refactor, as what we should be doing
//                 is accumulating every processor's print element contributions
//                 into a vector, and then doing a reduce on the whole
//                 vector at once.  Further, some print options, e.g. V(A,B)
//                 may require that we do a SUM across processors to get the
//                 correct result.  This can be done with a custom reduction
//                 operator rather than a generic maxAll vs. sumAll.  But
//                 again, this requires a refactor of the outputMgr class
//
// Scope         : private
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/11/04
//-----------------------------------------------------------------------------
double OutputMgr::getPrintValue(ParameterList::const_iterator iterParam, const N_LAS_Vector * solnVecPtr, const N_LAS_Vector * stateVecPtr, const N_LAS_Vector * storeVecPtr, const N_LAS_Vector * solnVecImagPtr)
{
  double result=0;
  bool resultFound = false;
  bool requireSum = false;

  SimulatorVariableContext theParamsContext = iterParam->getSimContext();

  switch (theParamsContext)
  {
    case INDEX:
      result = currentOutputter_->getIndex();
      resultFound = true;
      break;

    case CONSTANT:
      result = iterParam->getFixedValue();
      resultFound = true;
      break;

    case TEMPERATURE:
      result = circuitTemp_;
      resultFound = true;
      break;

    case TIME_VAR:
      result = circuitTime_; // *print_parameters.outputTimeScaleFactor_;
      resultFound = true;
      break;

    case FREQUENCY:
      result = circuitFrequency_;
      resultFound = true;
      break;

    case STEP_SWEEP_VAR:
    {
      int sweepVarIndex = iterParam->getExtraIndex1();
      result = stepParamVec_.at(sweepVarIndex).currentVal;
      resultFound = true;
    }
    break;

    case DC_SWEEP_VAR:
    {
      int dcSweepVarIndex = iterParam->getExtraIndex1();
      result = dcParamVec_.at(dcSweepVarIndex).currentVal;
      resultFound = true;
    }
    break;

    case GLOBAL_PARAMETER:
      result = devPtr_->getGlobalPar(iterParam->sVal());
      resultFound = true;
      break;

    case EXPRESSION:
      result = iterParam->getExpressionDataPointer()->evaluate(solnVecPtr, stateVecPtr, storeVecPtr);
      if( !(iterParam->getExpressionDataPointer()->numUnresolvedStringsChecked) )
      {
        // first time this expression has been resolved.  Check for unresolved symbols 
        // emit a fatal error if no expression in parallel has zero unresolved symbols 
        int numUnresolvedSymbols = iterParam->getExpressionDataPointer()->getNumUnresolvedStrings();
        
#ifdef Xyce_PARALLEL_MPI
        pdsCommPtr_->barrier();
        int minNumUnresolvedSymbols=0;
        pdsCommPtr_->minAll(&numUnresolvedSymbols, &minNumUnresolvedSymbols, 1);
        numUnresolvedSymbols=minNumUnresolvedSymbols;
#endif

        if (numUnresolvedSymbols > 0 )
        {
          string msg("Can't resolve all symbols in expression variable ");
          msg += iterParam->getExpressionDataPointer()->expPtr->get_expression();
          N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_FATAL_0,msg);
        }
        // set numUnresolvedStringsChecked flag to true so we only do this extra work once.
        iterParam->getExpressionDataPointer()->numUnresolvedStringsChecked = true;
      }
      resultFound = true;
      break;

    case SOLUTION_VAR:
    {
      int solVarIndex = iterParam->getExtraIndex1();
      if (solVarIndex == -1)
      {
        // off processor result, so return zero
        result = 0.0;
      }
      else
      {
        //result =(*solnVecPtr)[solVarIndex];
        result = solnVecPtr->getElementByGlobalIndex(solVarIndex);
      }
      resultFound = true;
    }
    break;

    case VOLTAGE_DIFFERENCE:
    {
      int solVarIndex1 = iterParam->getExtraIndex1();
      int solVarIndex2 = iterParam->getExtraIndex2();

      // the following may not be necessary as the epetra-map for the solution vector may just
      // return zero for an index of "-1".
      if (solVarIndex1 == -1 && solVarIndex2 == -1)
      {
        // off processor result, so return zero
        result = 0.0;
      }
      else if (solVarIndex1 == -1)
      {
        // solVarIndex1 is off proc.
        result = -solnVecPtr->getElementByGlobalIndex(solVarIndex2);
      }
      else if (solVarIndex2 == -1)
      {
        // solVarIndex2 is off proc.
        result = solnVecPtr->getElementByGlobalIndex(solVarIndex1);
      }
      else
      {
        result = solnVecPtr->getElementByGlobalIndex(solVarIndex1) - solnVecPtr->getElementByGlobalIndex(solVarIndex2);
      }
      resultFound=true;
      requireSum = true;
    }
    break;

    case STATE_VAR:
    {
      int stateVarIndex = iterParam->getExtraIndex1();
      // stateVecPtr can be NULL so check.
      if (stateVecPtr != NULL)
      {
        if (stateVarIndex == -1)
        {
          // off processor result, so return zero
          result = 0.0;
        }
        else
        {
          result =(*stateVecPtr)[stateVarIndex]; // stateVecPtr->getElementByGlobalIndex(stateVarIndex);
        }
      }
      resultFound = true;
    }
    break;

    case STORE_VAR:
    {
      int storeVarIndex = iterParam->getExtraIndex1();
      // stoVecPtr can be NULL so check.
      if (storeVecPtr != NULL)
      {
        if (storeVarIndex == -1)
        {
          // off processor result, so return zero
          result = 0.0;
        }
        else
        {
          result = (*storeVecPtr)[storeVarIndex]; // storeVecPtr->getElementByGlobalIndex(storeVarIndex);
        }
      }
      resultFound = true;
    }
    break;

    case DEVICE_PARAMETER:
    {
      string varName(iterParam->getQualifiedParameterOrFunctionName());
      result = devPtr_->getParam(varName);
      resultFound=true;
    }
    break;

    case OBJECTIVE_FUNCTION:
    {
      // this looks out of place here.  It should be in it's own update
      // function like measureManager_.update...
      string varType = iterParam->getQualifiedParameterOrFunctionName();
      if (objective[varType].var1.empty() && objective[varType].var2.empty())
      {
        // saving single values(potentially all values) as there isn't
        // any external data to tell us what needs to be saved
        result = objective[varType].save(solnVecPtr, stateVecPtr, storeVecPtr);
        resultFound=true;
      }
      else
      {
        // use user supplied external data to save only the simulation
        // results we really need.
        double v1=0.0;
        double v2=0.0;
        ParameterList vlist;
        vlist.push_back(N_UTL_Param("", 0));
        ParameterList::iterator it_vlist = vlist.begin();
        if (!objective[varType].var1.empty())
        {
          it_vlist->setTag(objective[varType].var1);
          v1 = getPrintValue(it_vlist, solnVecPtr, stateVecPtr, storeVecPtr);
        }
        if (!objective[varType].var2.empty())
        {
          it_vlist->setTag(objective[varType].var2);
          v2 = getPrintValue(it_vlist, solnVecPtr);
        }
        result = objective[varType].save(v1, v2, solnVecPtr, stateVecPtr, storeVecPtr);
        resultFound=true;
      }
    }
    break;

    case MEASURE_FUNCTION:
    {
      bool mFound=false;
      measureManager_.getMeasureValue(iterParam->getQualifiedParameterOrFunctionName(), result, mFound);
      resultFound=mFound;
    }
    break;

    // these are for solution variables under AC analysis where
    // solution is complex and could be printed as real, imaginary,
    // magnitude or phase
    case SOLUTION_VAR_REAL:
    {
      int solVarIndex = iterParam->getExtraIndex1();
      if (solVarIndex == -1)
      {
        // off processor result, so return zero
        result = 0.0;
      }
      else
      {
        result = solnVecPtr->getElementByGlobalIndex(solVarIndex);
      }
      resultFound = true;
    }
    break;

    case SOLUTION_VAR_IMAG:
    {
      int solVarIndex = iterParam->getExtraIndex1();
      if (solVarIndex == -1 ||(solnVecImagPtr == 0))
      {
        // off processor result, so return zero
        result = 0.0;
      }
      else
      {
        result = solnVecImagPtr->getElementByGlobalIndex(solVarIndex);
      }
      resultFound = true;
    }
    break;

    case SOLUTION_VAR_MAG:
    {
      int solVarIndex = iterParam->getExtraIndex1();
      if (solVarIndex == -1 )
      {
        // off processor result, so return zero
        result = 0.0;
      }
      else
      {
        double realComp = solnVecPtr->getElementByGlobalIndex(solVarIndex);
        double imagComp = 0.0;
        if (solnVecImagPtr != 0 )
          imagComp = solnVecImagPtr->getElementByGlobalIndex(solVarIndex);
        result =  sqrt(pow(realComp, 2) + pow(imagComp, 2) );
      }
      resultFound = true;
    }
    break;

    case SOLUTION_VAR_PHASE:
    {
      int solVarIndex = iterParam->getExtraIndex1();
      if (solVarIndex == -1 )
      {
        // off processor result, so return zero
        result = 0.0;
      }
      else
      {
        double realComp = solnVecPtr->getElementByGlobalIndex(solVarIndex);
        double imagComp = 0.0;
        if (solnVecImagPtr != 0 )
          imagComp = solnVecImagPtr->getElementByGlobalIndex(solVarIndex);
        // phase in radians is atan2(imag(z), real(z)).
        result =  atan2(imagComp, realComp );
      }
      resultFound = true;
    }
    break;

    case SOLUTION_VAR_DB:
    {
      int solVarIndex = iterParam->getExtraIndex1();
      if (solVarIndex == -1 )
      {
        // off processor result, so return zero
        result = 0.0;
      }
      else
      {
        double realComp = solnVecPtr->getElementByGlobalIndex(solVarIndex);
        double imagComp = 0.0;
        if (solnVecImagPtr != 0 )
          imagComp = solnVecImagPtr->getElementByGlobalIndex(solVarIndex);
        // DB = 20*log10(magnitude)
        result =  20*log10(sqrt(pow(realComp, 2) + pow(imagComp, 2)) );
      }
      resultFound = true;
    }
    break;

    case VOLTAGE_DIFFERENCE_REAL:
    {
      int solVarIndex1 = iterParam->getExtraIndex1();
      int solVarIndex2 = iterParam->getExtraIndex2();
      // the following may not be necessary as the epetra-map for the solution vector may just
      // return zero for an index of "-1".
      if (solVarIndex1 == -1 && solVarIndex2 == -1)
      {
        // off processor result, so return zero
        result = 0.0;
      }
      else if (solVarIndex1 == -1)
      {
        // solVarIndex1 is off proc.
        result = -solnVecPtr->getElementByGlobalIndex(solVarIndex2);
      }
      else if (solVarIndex2 == -1)
      {
        // solVarIndex2 is off proc.
        result = solnVecPtr->getElementByGlobalIndex(solVarIndex1);
      }
      else
      {
        result = solnVecPtr->getElementByGlobalIndex(solVarIndex1) - solnVecPtr->getElementByGlobalIndex(solVarIndex2);
      }
      resultFound=true;
      requireSum = true;
    }
    break;

    case VOLTAGE_DIFFERENCE_IMAG:
    {
      int solVarIndex1 = iterParam->getExtraIndex1();
      int solVarIndex2 = iterParam->getExtraIndex2();
      // the following may not be necessary as the epetra-map for the solution vector may just
      // return zero for an index of "-1".
      if (solVarIndex1 == -1 && solVarIndex2 == -1)
      {
        // off processor result, so return zero
        result = 0.0;
      }
      else if (solVarIndex1 == -1)
      {
        // solVarIndex1 is off proc.
        result = -solnVecImagPtr->getElementByGlobalIndex(solVarIndex2);
      }
      else if (solVarIndex2 == -1)
      {
        // solVarIndex2 is off proc.
        result = solnVecImagPtr->getElementByGlobalIndex(solVarIndex1);
      }
      else
      {
        result = solnVecImagPtr->getElementByGlobalIndex(solVarIndex1) - solnVecImagPtr->getElementByGlobalIndex(solVarIndex2);
      }
      resultFound=true;
      requireSum = true;
    }
    break;

    case NODE_OR_DEVICE_NAME:
      // we can ignore this becasue it's the name associated with a prior
      // I(), V() or N() param.
      break;

    case UNDEFINED:
    default:
      break;
  }

  // check filter option
  if (filterGiven_ &&(fabs(result) < filter_))
  {
    result = 0.0;
  }

  // If running in parallel, make value of "result" be the same on
  // each processor.
#ifdef Xyce_PARALLEL_MPI
  pdsCommPtr_->barrier();

  double final=0;

  // HACK ALERT:  Voltage differences are special, and could need to be a
  // combination of data from off-processor rather than a global max/min.
  // Those special cases require us to use sumAll instead.
  //  This is RIPE for refactor and must be undertaken after release 6.0!
  if (!requireSum)
  {
    pdsCommPtr_->maxAll(&result, &final, 1);
    if (final == 0.0) pdsCommPtr_->minAll(&result, &final, 1);
  }
  else
  {
    pdsCommPtr_->sumAll(&result, &final, 1);
  }

  result = final;

#endif

  return result;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::outputAC
// Purpose       : .PRINT output for ac runs
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey
// Creation Date : 8/5/08
//-----------------------------------------------------------------------------
void OutputMgr::outputAC(double frequency, const N_LAS_Vector * real_solution_vector, const N_LAS_Vector *imaginary_solution_vector)
{
  circuitFrequency_ = frequency;

  if (!activeOutputterStack_.empty())
    for (std::vector<Outputter::Interface *>::const_iterator it = activeOutputterStack_.back().begin(); it != activeOutputterStack_.back().end(); ++it)
      (*it)->outputAC(frequency, real_solution_vector, imaginary_solution_vector);
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::outputAC
// Purpose       : .PRINT output for ac runs
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey
// Creation Date : 8/5/08
//-----------------------------------------------------------------------------
void OutputMgr::outputHomotopy(const std::vector<string> & parameter_names, const std::vector<double> & param_values, const N_LAS_Vector * solution_vector)
{
  if (!activeOutputterStack_.empty())
    for (std::vector<Outputter::Interface *>::const_iterator it = activeOutputterStack_.back().begin(); it != activeOutputterStack_.back().end(); ++it)
      (*it)->outputHomotopy(parameter_names, param_values, solution_vector);
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::outputAC
// Purpose       : .PRINT output for ac runs
// Special Notes :
// Scope         : public
// Creator       : Todd Coffey
// Creation Date : 8/5/08
//-----------------------------------------------------------------------------
void OutputMgr::outputMORTF(bool origSystem, const double & freq, const Teuchos::SerialDenseMatrix<int, std::complex<double> >& H)
{
  if (!activeOutputterStack_.empty())
    for (std::vector<Outputter::Interface *>::const_iterator it = activeOutputterStack_.back().begin(); it != activeOutputterStack_.back().end(); ++it)
      (*it)->outputMORTF(origSystem, freq, H);
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
    string baseoutputname;
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
    N_IO_MMIO::MM_typecode matcode;
    string cfile = baseoutputname + ".Chat";
    string gfile = baseoutputname + ".Ghat";
    string bfile = baseoutputname + ".Bhat";
    string lfile = baseoutputname + ".Lhat";
    c_file = fopen(cfile.c_str(), "w");
    g_file = fopen(gfile.c_str(), "w");
    b_file = fopen(bfile.c_str(), "w");
    l_file = fopen(lfile.c_str(), "w");
    if (c_file == NULL || g_file == NULL || b_file == NULL || l_file == NULL)
    {
      string msg="Error: Cannot open one of the ROM files for output: " + cfile + ", " + gfile + ", " + bfile + ", " + lfile;
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
    }
    mm_initialize_typecode(&matcode);
    mm_set_matrix(&matcode);
    mm_set_array(&matcode);
    mm_set_general(&matcode);
    mm_set_real(&matcode);

    int ret = 0;

    // Write the headers
    ret = N_IO_MMIO::mm_write_banner(g_file, matcode);
    ret = N_IO_MMIO::mm_write_banner(c_file, matcode);
    ret = N_IO_MMIO::mm_write_banner(b_file, matcode);
    ret = N_IO_MMIO::mm_write_banner(l_file, matcode);

    // Write the matrix array sizes
    ret = N_IO_MMIO::mm_write_mtx_array_size(g_file, Ghat.numRows(), Ghat.numCols());
    ret = N_IO_MMIO::mm_write_mtx_array_size(c_file, Chat.numRows(), Chat.numCols());
    ret = N_IO_MMIO::mm_write_mtx_array_size(b_file, Bhat.numRows(), Bhat.numCols());
    ret = N_IO_MMIO::mm_write_mtx_array_size(l_file, Lhat.numRows(), Lhat.numCols());

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
  string baseoutputname;
  if (pos > 0)
  {
    baseoutputname = netListFilename_.substr(0, pos);
  }
  else
  {
    baseoutputname = netListFilename_;
  }

  // Create file string for Chat and Ghat
  string gfile = baseoutputname + ".Ghat";
  string cfile = baseoutputname + ".Chat";

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
    N_IO_MMIO::MM_typecode matcode;

    string bfile = baseoutputname + ".Bhat";
    string lfile = baseoutputname + ".Lhat";

    b_file = fopen(bfile.c_str(), "w");
    l_file = fopen(lfile.c_str(), "w");
    if (b_file == NULL || l_file == NULL)
    {
      string msg="Error: Cannot open one of the ROM files for output: " + bfile + ", " + lfile;
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
    }
    mm_initialize_typecode(&matcode);
    mm_set_matrix(&matcode);
    mm_set_array(&matcode);
    mm_set_general(&matcode);
    mm_set_real(&matcode);

    int ret = 0;

    // Write the headers
    ret = N_IO_MMIO::mm_write_banner(b_file, matcode);
    ret = N_IO_MMIO::mm_write_banner(l_file, matcode);

    // Write the matrix array sizes
    ret = N_IO_MMIO::mm_write_mtx_array_size(b_file, Bhat.numRows(), Bhat.numCols());
    ret = N_IO_MMIO::mm_write_mtx_array_size(l_file, Lhat.numRows(), Lhat.numCols());

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
      for (std::vector<Outputter::Interface *>::const_iterator it = activeOutputterStack_.back().begin(); it != activeOutputterStack_.back().end(); ++it)
        (*it)->finishOutput();

#ifdef Xyce_USE_HDF5
  if (hdf5FileNameGiven_)
  {
    closeHDF5Output();
  }
#endif // Xyce_USE_HDF5
}

void OutputMgr::resetOutput()
{
  if (!activeOutputterStack_.empty())
    for (std::vector<Outputter::Interface *>::const_iterator it = activeOutputterStack_.back().begin(); it != activeOutputterStack_.back().end(); ++it)
      (*it)->resetOutput();
}

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
  vector <N_ANP_SweepParam>::const_iterator firstParam = step_sweep_parameters.begin();
  vector <N_ANP_SweepParam>::const_iterator lastParam = step_sweep_parameters.end();

  if (!(step_sweep_parameters.empty()))
  {
    stepParamVec_ = step_sweep_parameters;
  }

  if (!activeOutputterStack_.empty())
    for (std::vector<Outputter::Interface *>::const_iterator it = activeOutputterStack_.back().begin(); it != activeOutputterStack_.back().end(); ++it)
      (*it)->outputHB(timePoints, freqPoints, timeDomainSolnVec, freqDomainSolnVecReal, freqDomainSolnVecImaginary);
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
  cout << std::endl << "OutputMgr::outputRESULT" << std::endl;
#endif

  std::vector<N_ANP_SweepParam>::iterator iterParam;
  std::vector<N_ANP_SweepParam>::iterator firstParam=stepParamVec_.begin();
  std::vector<N_ANP_SweepParam>::iterator lastParam=stepParamVec_.end();

  N_UTL_ExpressionData * expDataPtr = NULL;

  int ires=0;
  string delim("");
  int width = 17;
  int precision = 8;

  if (getProcID() == 0)
  {
    if (!RESULTinitialized_)
    {
      RESULTinitialized_ = true;

      string resultfilename("");
      if (netListFilename_  != "")
      {
        resultfilename = netListFilename_ + ".res";
      }
      else
      {
        resultfilename = "output.res";
      }
      resultStreamPtr_ = new ofstream(resultfilename.c_str());

      resultStreamPtr_->setf(ios::scientific);
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
       (*resultStreamPtr_) << delim << expDataPtr->expression;
      }

     (*resultStreamPtr_) << "\n";

    } // RESULTinitialized_

    if (delim == "") resultStreamPtr_->width(8);
    else             resultStreamPtr_->width(0);
    resultStreamPtr_->setf(ios::left, ios::adjustfield);
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
      for (std::vector<Outputter::Interface *>::const_iterator it = activeOutputterStack_.back().begin(); it != activeOutputterStack_.back().end(); ++it)
        (*it)->finishOutputStep();

    // Deal with the result file:
    if (resultStreamPtr_)
    {
      (*resultStreamPtr_) << "End of Xyce(TM) Parameter Sweep" << std::endl;
    }

    if (resultStreamPtr_ != &cout && resultStreamPtr_)
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
// Purpose       : if any macro analysis was specified(like objective or measure)
//                 output results
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL Electrical and Microsystem Modeling
// Creation Date : 11/05/08
//-----------------------------------------------------------------------------
void OutputMgr::outputMacroResults()
{
#ifdef Xyce_DEBUG_IO
  if (!objective.empty())
  {
    std::cout << std::endl
      << " ***** Analysis Functions ***** " << std::endl
      << std::endl;

    std::map<string, N_IO_Objective>::iterator ob = objective.begin();
    std::map<string, N_IO_Objective>::iterator ob_end = objective.end();
    for (; ob != ob_end ; ++ob)
    {
      string name((*ob).first);
      double val = objective[name].evaluate();
      cout << name << " = " << val << std::endl;
    }
  }
#endif

  // This is a null output stream that helps ensure a function that needs to be called
  // on every processor in parallel only outputs on one processor.
  Teuchos::oblackholestream outputBHS;
  std::ofstream outputFileStream;

  // Output the Measure results only if .measure is being performed on any variables.
  if (measureManager_.isMeasureActive())
  {
    // Output the Measure results to std::cout.
    // Make sure the function gets called on all processors, but only one outputs it.
    if (getProcID() == 0)
    {
      measureManager_.outputResults( std::cout );
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
    ofstream responseOFS;
    string outputResponseFilename;
    if (responseFileNameGiven_)
    {
      // do need the suffix in this case as a name was specified
      outputResponseFilename = responseFileName_;
    }
    else
    {
      outputResponseFilename = "response.out"; // + filenameSuffix_;
    }

    responseOFS.open(outputResponseFilename.c_str());

    std::vector< pair< string, string> >::iterator currRespItr = responseFunctionsRequested_.begin();
    std::vector< pair< string, string> >::iterator endRespItr = responseFunctionsRequested_.end();
    while (currRespItr != endRespItr)
    {
      double respvalue = 0.0;
      bool found = false;
      // need to parse name from XXX_X:name So find last ":" and extract
      // remaining string as value that dakota is looking for.
      string tempName=currRespItr->first;
      std::string::size_type beginingOfVarName = tempName.find_last_of(":");
      //beginingOfVarName++;

      if (beginingOfVarName != std::string::npos)
      {
        int numChars =(currRespItr->first).length() - beginingOfVarName;
        tempName.assign(currRespItr->first, beginingOfVarName+1, numChars);
      }
      ExtendedString es(tempName);
      std::cout << "Calling getMeasureValue with " << es.toUpper() << std::endl;
      measureManager_.getMeasureValue(es.toUpper(), respvalue, found);
      if (! found)
        responseOFS << "not found" << "   " << currRespItr->first << std::endl;
      else
        responseOFS << respvalue   << "   " << currRespItr->first << std::endl;
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
  N_UTL_Param parameter;

  sVal.toUpper();
  if (sVal.size() < 3 || sVal[1] != '(' || sVal[sVal.size()-1] != ')') {
    std::string msg("OutputMgr::registerResponseVars: response var not of format V() or I(): '");
    msg += objString;
    msg += "'";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, msg);
  }
  numResponseVars_++;
  responseVarList_.resize(numResponseVars_);
  ParameterList & pList = responseVarList_[numResponseVars_-1];
  parameter.setTag(sVal.substr(0, 1));
  parameter.setVal(1.0);
  pList.push_back(parameter);
  parameter.setTag(sVal.substr(2, sVal.size()-3));
  parameter.setVal(0.0);
  pList.push_back(parameter);
  pl_i = pList.begin();
  setParamContextType_(pl_i);

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
  std::vector<ParameterList >::iterator currentVar;
  std::vector<ParameterList >::iterator endVar;
  ParameterList::iterator pl_i;
  int expValue = 0;

  // save the resuts of any end point variables if there are any
  ParameterList vlist;
  vlist.push_back(N_UTL_Param("", 0));
  ParameterList::iterator it_vlist = vlist.begin();

  // save the independant variable
  int varNumber = 0;
  if (dcParamVec_.empty())  // DNS: is this a reliable way to determine if transient?
  {
    responseVarPtr_->at(varNumber) = circuitTime_;
  }
  else
  {
    responseVarPtr_->at(varNumber) =
      dcParamVec_[ dcLoopNumber_ ].currentVal;
  }
  varNumber++;

  // loop over the response variable list
  currentVar = responseVarList_.begin();
  endVar = responseVarList_.end();
  for (; currentVar != endVar ; ++currentVar)
  {
    pl_i =(*currentVar).begin();
    double result = getPrintValue(pl_i, solnVecPtr);

    //std::cout << result << std::endl;

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
void OutputMgr::getMeasureValue(const std::string &name, double &result, bool &found)
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
  std::cout << "In OutputMgr::remeasure " << std::endl;

  // get file to open for remeasure
  std::string existingDataFile = commandLine_.getArgumentValue("-remeasure");
  std::cout << "file to reprocess through measure functions. " << existingDataFile << std::endl;

  // open file for reading
  // just support PRN format for now
  N_IO_OutputFileBase * fileToReMeasure = new N_IO_OutputPrn();
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

//   std::cout << "Original var names: " << std::endl;
//   for (int i=0; i<fileVarNames.size(); i++)
//   {
//     std::cout << "\"" << fileVarNames[i] << "\", ";
//   }
//   std::cout << std::endl;
//
//   fileToReMeasure->convertOutputNamesToSolVarNames(fileVarNames);
//
//   std::cout << "extracted var names: " << std::endl;
//   for (int i=0; i<fileVarNames.size(); i++)
//   {
//     std::cout << "\"" << fileVarNames[i] << "\", ";
//   }
//   std::cout << std::endl;

  // assume analysis type is DC(Xyce hasn't processed the analysis type yet
  // when this function is called.  In the next loop if we find a column of
  // data for TIME then we can treat this as a transient analysis.

  PRINTType_ = PrintType::DC;
  int timeIndex=0;
  // set up allNodes map to map var names to indices in the solution vector
  int numVars = fileVarNames.size();
  for (int i=0; i<numVars; i++)
  {
    allNodes_[fileVarNames[i]]=make_pair(i, 0);
    ExtendedString tmpStr(fileVarNames[i]);
    tmpStr.toUpper();
    allNodes_[tmpStr]=make_pair(i, 0);
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
    RCP< N_LAS_Vector > varValuesVecRCP(varValuesVecPtr, false);
    if (PRINTType_ == PrintType::TRAN)
    {
      circuitTime_=(*varValuesVecPtr)[timeIndex];
      measureManager_.updateTranMeasures(circuitTime_, varValuesVecRCP);
    }
    else
    {
      measureManager_.updateDcMeasures(dcParamVec_, varValuesVecRCP);
    }
    varValuesVecPtr->putScalar(0);

  }

  delete varValuesVecPtr;
  delete aParMapPtr;

  // This is a null output stream that helps ensure a function that needs to be called
  // on every processor in parallel only outputs on one processor.
  Teuchos::oblackholestream outputBHS;
  std::ofstream outputFileStream;

  // Output the Measure results to std::cout.
  // Make sure the function gets called on all processors, but only one outputs it.
  if (getProcID() == 0)
  {
    measureManager_.outputResults( std::cout );
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

} // namespace IO
} // namespace Xyce

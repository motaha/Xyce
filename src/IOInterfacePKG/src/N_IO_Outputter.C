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
// Filename       : $RCSfile: N_IO_Outputter.C,v $
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
// Revision Number: $Revision: 1.16.2.21 $
//
// Revision Date  : $Date: 2013/10/03 17:23:43 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>
// ---------- Standard Includes ----------

#include <N_UTL_Misc.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include <N_PDS_Comm.h>
#include <N_IO_Outputter.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_CmdParse.h>

#include <N_ANP_AnalysisInterface.h>
#include <N_ANP_AnalysisManager.h>
#include <N_MPDE_Manager.h>

#include <N_UTL_Version.h>

#include <N_LAS_Vector.h>
#include <N_LAS_BlockVector.h>

namespace Xyce {
namespace IO {

namespace { // <unnamed>

//-----------------------------------------------------------------------------
// Function      : OutputMgr::getTimeDateStamp
// Purpose       : get current date and time and format for .PRINT output
// Special Notes : inline
// Scope         : private
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
std::string getTimeDateStamp()
{
  const time_t now = time( NULL);
  char timeDate[ 40 ];

  // format for output
  strftime( timeDate, 40, "TIME='%I:%M:%S %p' DATE='%b %d, %Y' ", localtime( &now));

  return std::string( timeDate);
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::getTecplotTimeDateStamp
// Purpose       : Get current date and time and format for .PRINT output
// Special Notes : tecplot version of getTimeDateStamp.
// Scope         : private
// Creator       : Eric Keiter, SNL
// Creation Date : 9/6/04
//-----------------------------------------------------------------------------
std::string getTecplotTimeDateStamp()
{
  const time_t now = time( NULL);
  char timeDate[ 40 ];

  // format for output
  strftime( timeDate, 40, "TIME= \" %I:%M:%S %p %b %d, %Y \" ", localtime( &now));

  return std::string( timeDate);
}

void tecplotFreqHeader(OutputMgr &output_manager, const PrintParameters &print_parameters, std::ostream &stream, int counter);
void fixupColumns(const OutputMgr &output_manager, PrintParameters &print_parameters);

std::ostream &printHeader(std::ostream &os, const PrintParameters &print_parameters);
std::ostream &printHeader(std::ostream &os, const Table::ColumnList &column_list, const std::string &delimiter);
std::ostream &printHeader(std::ostream &os, const Table::Column &column);
std::ostream &printValue(std::ostream &os, const Table::Column &column, const std::string &delimiter, const int column_index, double value);

std::string outputFilename(const PrintParameters &print_parameters, const std::string &netListFilename)
{
  if (!print_parameters.filename_.empty())
  {
    return print_parameters.filename_ + print_parameters.suffix_;
  }
  else
  {
    return netListFilename + print_parameters.suffix_ + print_parameters.extension_;
  }
}

} // namespace <unnamed>

namespace Outputter {

TimePrn::TimePrn(OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    firstTimePrint_(true),
    index_(0),
    outFilename_(),
    suffix_(),
    outStreamPtr_(0),
    headerPrintCalls_(0)
{
  if (printParameters_.extension_.empty())
    printParameters_.extension_ = ".prn";

  fixupColumns(outputManager_, printParameters_);
}

TimePrn::~TimePrn()
{
  outputManager_.closeFile(outStreamPtr_);
}

void TimePrn::doParse()
{
  outFilename_ = outputFilename(printParameters_, outputManager_.getNetListFilename());
}

void TimePrn::doOpen()
{
  if (outputManager_.getProcID() == 0 && outStreamPtr_ == 0)
    outStreamPtr_ = outputManager_.openFile(outFilename_);
}

void TimePrn::timeHeader()
{
  if (outputManager_.getProcID() == 0 && headerPrintCalls_ == 0) {
    printHeader(*outStreamPtr_, printParameters_);
  }
}

void TimePrn::doOutputTime(const N_LAS_Vector * solnVecPtr, const N_LAS_Vector * stateVecPtr, const N_LAS_Vector * storeVecPtr)
{
  outputManager_.setCurrentOutputter(this);

  if (outputManager_.getProcID() == 0)
  {
    if (firstTimePrint_) // Setup Output Stream and Print Out Header
    {
      doOpen();

      timeHeader();

      ++headerPrintCalls_;

      index_ = 0;

      firstTimePrint_ = false;
    }
  } // procID

  std::ostream &os = *outStreamPtr_;
  outputManager_.getCommPtr()->barrier();

  int column_index = 0;
  for (ParameterList::const_iterator it = printParameters_.variableList_.begin() ; it != printParameters_.variableList_.end(); ++it, ++column_index)
  {
    double result = outputManager_.getPrintValue(it, solnVecPtr, stateVecPtr, storeVecPtr);
    if ((*it).getSimContext() == TIME_VAR)
      result *= printParameters_.outputTimeScaleFactor_;
    
    if (outputManager_.getProcID() == 0)
      printValue(os, printParameters_.table_.columnList_[column_index], printParameters_.delimiter_, column_index, result);

    outputManager_.getCommPtr()->barrier();
  }

  ++index_;

  if (outputManager_.getProcID() == 0)
    os << std::endl;
}

void TimePrn::doFinishOutput()
{
  if (outputManager_.getProcID() == 0)
  {

    if (outStreamPtr_)
    {
      if (!outputManager_.getSTEPEnabledFlag() && outputManager_.getPrintEndOfSimulationLine())
      {
        (*outStreamPtr_) << "End of Xyce(TM) Simulation" << std::endl;
        // set to zero so that if there is more output, as in MPDE,
        // a new header will get printed as well.  Can't do this
        // at the top of this function as it would break .STEP's use
        // of this function.
        headerPrintCalls_ = 0;
      }
    }

    if (!outputManager_.getSTEPEnabledFlag())
    {
      outputManager_.closeFile(outStreamPtr_);
      outStreamPtr_ = 0;
    }
  } // procID

  firstTimePrint_ = true;
  index_ = 0;
}

void
TimePrn::doFinishOutputStep()
{
  // Deal with the *prn file:
  if (outStreamPtr_)
  {
    if (outputManager_.getFormat() != Format::PROBE && outputManager_.getPrintEndOfSimulationLine())
      (*outStreamPtr_) << "End of Xyce(TM) Parameter Sweep" << std::endl;

    outputManager_.closeFile(outStreamPtr_);
    outStreamPtr_ = 0;
  }
}

FrequencyPrn::FrequencyPrn(OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    firstTimePrint_(true),
    outFilename_(),
    suffix_(),
    outStreamPtr_(0),
    headerPrintCalls_(0),
    stepCount_(0)
{
  if (printParameters_.extension_.empty())
    printParameters_.extension_ = ".FD.prn";

  fixupColumns(outputManager_, printParameters_);
}

FrequencyPrn::~FrequencyPrn()
{
  outputManager_.closeFile(outStreamPtr_);
}

void FrequencyPrn::doParse()
{
  outFilename_ = outputFilename(printParameters_, outputManager_.getNetListFilename());
}

void FrequencyPrn::doOpen()
{
  if (outputManager_.getProcID() == 0 && outStreamPtr_ == 0)
    outStreamPtr_ = outputManager_.openFile(outFilename_);
}

void
FrequencyPrn::doOutputFrequency(double frequency, const N_LAS_Vector *real_solution_vector, const N_LAS_Vector *imaginary_solution_vector)
{
  outputManager_.setCurrentOutputter(this);

  if (outputManager_.getProcID() == 0)
  {
    if (firstTimePrint_) //Setup Output Stream and Print Out Header
    {
      doOpen();

      index_ = 0;

      printHeader(*outStreamPtr_, printParameters_);

      ++stepCount_;
      firstTimePrint_ = false;
    }
  }

  std::ostream &os = *outStreamPtr_;

  int column_index = 0;
  for (ParameterList::const_iterator it = printParameters_.variableList_.begin(); it != printParameters_.variableList_.end(); ++it, ++column_index)
  {
    double varValue = outputManager_.getPrintValue(it, real_solution_vector, NULL, NULL, imaginary_solution_vector);

    if (outputManager_.getProcID() == 0)
      printValue(os, printParameters_.table_.columnList_[column_index], printParameters_.delimiter_, column_index, varValue);

    outputManager_.getCommPtr()->barrier();
  }

  ++index_;

  if (outputManager_.getProcID() == 0)
    os << std::endl;
}

void FrequencyPrn::doFinishOutput()
{
  if (outputManager_.getProcID() == 0)
  {
    if (outStreamPtr_)
    {
      if (!outputManager_.getSTEPEnabledFlag() && outputManager_.getPrintEndOfSimulationLine())
      {
        (*outStreamPtr_) << "End of Xyce(TM) Simulation" << std::endl;
        headerPrintCalls_ = 0;
      }

      if (!outputManager_.getSTEPEnabledFlag())
      {
        outputManager_.closeFile(outStreamPtr_);
        outStreamPtr_ = 0;
      }
    }
  } // procID

  firstTimePrint_ = true;
}

void
FrequencyPrn::doFinishOutputStep()
{
  // Deal with the *prn file:
  if (outStreamPtr_)
  {
    if (outputManager_.getFormat() != Format::PROBE && outputManager_.getPrintEndOfSimulationLine())
      (*outStreamPtr_) << "End of Xyce(TM) Parameter Sweep" << std::endl;

    outputManager_.closeFile(outStreamPtr_);
    outStreamPtr_ = 0;
  }
}

TimeCSV::TimeCSV(OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    firstTimePrint_(true),
    outFilename_(),
    suffix_(),
    outStreamPtr_(0),
    headerPrintCalls_(0)
{
  if (printParameters_.extension_.empty())
    printParameters_.extension_ = ".csv";

  fixupColumns(outputManager_, printParameters_);
}

TimeCSV::~TimeCSV()
{
  outputManager_.closeFile(outStreamPtr_);
}


void TimeCSV::doParse()
{
  outFilename_ = outputFilename(printParameters_, outputManager_.getNetListFilename());
}

void TimeCSV::doOpen()
{
  if (outputManager_.getProcID() == 0 && outStreamPtr_ == 0)
    outStreamPtr_ = outputManager_.openFile(outFilename_);
}

void TimeCSV::timeHeader()
{
  if (outputManager_.getProcID() == 0)
  {
    if (headerPrintCalls_ == 0)
    {
      printHeader(*outStreamPtr_, printParameters_);
    }
  } // procID

  ++headerPrintCalls_;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::outputPRINT_
// Purpose       : .PRINT output
// Special Notes :
// Scope         : private
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/10/00
//-----------------------------------------------------------------------------
void TimeCSV::doOutputTime(const N_LAS_Vector * solnVecPtr, const N_LAS_Vector * stateVecPtr, const N_LAS_Vector * storeVecPtr)
{
  outputManager_.setCurrentOutputter(this);
  double time = outputManager_.getCircuitTime();

  if (outputManager_.getProcID() == 0)
  {

    if (firstTimePrint_) // Setup Output Stream and Print Out Header
    {
      doOpen();

      index_ = 0;

      timeHeader();
      firstTimePrint_ = false;
    }
  } // procID

  std::ostream &os = *outStreamPtr_;
  outputManager_.getCommPtr()->barrier();

  int column_index = 0;
  for (ParameterList::const_iterator it = printParameters_.variableList_.begin() ; it != printParameters_.variableList_.end(); ++it, ++column_index)
  {
    double result = outputManager_.getPrintValue(it, solnVecPtr, stateVecPtr, storeVecPtr);
    if ((*it).getSimContext() == TIME_VAR)
      result *= printParameters_.outputTimeScaleFactor_;

    if (outputManager_.getProcID() == 0)
      printValue(os, printParameters_.table_.columnList_[column_index], printParameters_.delimiter_, column_index, result);

    outputManager_.getCommPtr()->barrier();
  }

  ++index_;

  if (outputManager_.getProcID() == 0)
    os << std::endl;
}

void TimeCSV::doFinishOutput()
{
  if (outputManager_.getProcID() == 0)
  {
    if (outStreamPtr_)
    {
      if (!outputManager_.getSTEPEnabledFlag())
      {
        outputManager_.closeFile(outStreamPtr_);
        outStreamPtr_ = 0;
      }
    }
  } // procID

  firstTimePrint_ = true;
}

void
TimeCSV::doFinishOutputStep()
{
  if (outStreamPtr_) {
    outputManager_.closeFile(outStreamPtr_);
    outStreamPtr_ = 0;
  }
}

FrequencyCSV::FrequencyCSV(OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    firstTimePrint_(true),
    outFilename_(),
    suffix_(),
    outStreamPtr_(0),
    headerPrintCalls_(0),
    stepCount_(0)
{
  if (printParameters_.extension_.empty())
    printParameters_.extension_ = ".FD.csv";

  fixupColumns(outputManager_, printParameters_);
}

FrequencyCSV::~FrequencyCSV()
{
  outputManager_.closeFile(outStreamPtr_);
}


void FrequencyCSV::doParse()
{
  outFilename_ = outputFilename(printParameters_, outputManager_.getNetListFilename());
}

void FrequencyCSV::doOpen()
{
  if (outputManager_.getProcID() == 0 && outStreamPtr_ == 0)
    outStreamPtr_ = outputManager_.openFile(outFilename_);
}

void
FrequencyCSV::doOutputFrequency(double frequency, const N_LAS_Vector *real_solution_vector, const N_LAS_Vector *imaginary_solution_vector)
{
  outputManager_.setCurrentOutputter(this);

  if (outputManager_.getProcID() == 0)
  {
    if (firstTimePrint_) //Setup Output Stream and Print Out Header
    {
      doOpen();

      printHeader(*outStreamPtr_, printParameters_);

      index_ = 0;
      ++stepCount_;
      firstTimePrint_ = false;
    }
  }

  std::ostream &os = *outStreamPtr_;

  int column_index = 0;
  for (ParameterList::const_iterator it = printParameters_.variableList_.begin(); it != printParameters_.variableList_.end(); ++it, ++column_index)
  {
    double varValue = outputManager_.getPrintValue(it, real_solution_vector, NULL, NULL, imaginary_solution_vector);

    if (outputManager_.getProcID() == 0)
      printValue(os, printParameters_.table_.columnList_[column_index], printParameters_.delimiter_, column_index, varValue);

    outputManager_.getCommPtr()->barrier();
  }

  ++index_;

  if (outputManager_.getProcID() == 0)
    os << std::endl;
}

void FrequencyCSV::doFinishOutput()
{
  if (outputManager_.getProcID() == 0)
  {
    if (outStreamPtr_) {

      if (!outputManager_.getSTEPEnabledFlag())
      {
        outputManager_.closeFile(outStreamPtr_);
        outStreamPtr_ = 0;
      }
    }
  }

  firstTimePrint_ = true;
}

void
FrequencyCSV::doFinishOutputStep()
{
  // Deal with the *prn file:
  if (outStreamPtr_)
  {
    if (outputManager_.getFormat() != Format::PROBE && outputManager_.getPrintEndOfSimulationLine())
      (*outStreamPtr_) << "End of Xyce(TM) Parameter Sweep" << std::endl;
  }

  outputManager_.closeFile(outStreamPtr_);
  outStreamPtr_ = 0;
}

TimeTecPlot::TimeTecPlot(OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    firstTimePrint_(true),
    outFilename_(),
    suffix_(),
    outStreamPtr_(0),
    headerPrintCalls_(0)
{
  if (printParameters_.extension_.empty())
    printParameters_.extension_ = ".dat";

  fixupColumns(outputManager_, printParameters_);
}

TimeTecPlot::~TimeTecPlot()
{
  outputManager_.closeFile(outStreamPtr_);
}


void TimeTecPlot::doParse()
{
  outFilename_ = outputFilename(printParameters_, outputManager_.getNetListFilename());
}

void TimeTecPlot::doOpen()
{
  if (outputManager_.getProcID() == 0 && outStreamPtr_ == 0)
    outStreamPtr_ = outputManager_.openFile(outFilename_);
}

void TimeTecPlot::timeHeader()
{
  std::ostream &os = *outStreamPtr_;

  index_ = 0;

  if (outputManager_.getProcID() == 0)
  {
    int tecplotHeaderPrecision_ = 2;

    if (headerPrintCalls_ == 0)
    {
      os << "TITLE = \"" << outputManager_.getNetListFilename() << " - " << outputManager_.getTitle() << "\", " << std::endl;
      os << "\tVARIABLES = ";

      // output the user-specified solution vars:
      for (ParameterList::const_iterator iterParam = printParameters_.variableList_.begin() ; iterParam != printParameters_.variableList_.end(); ++iterParam)
        os << "\" " << (*iterParam).tag() << "\" " << std::endl;

      // output some AUXDATA
      os << "DATASETAUXDATA " << getTecplotTimeDateStamp() << std::endl;

      if (!outputManager_.getTempSweepFlag())
      {
        os.setf(ios::scientific);
        os.precision( tecplotHeaderPrecision_);
        os << "DATASETAUXDATA TEMP = \"" << outputManager_.getCircuitTemp() << " \"" << std::endl;
      }
    } // print header calls=0

    os << "ZONE F=POINT ";

    if (outputManager_.getStepParamVec().empty())
    {
      os << "T=\"Xyce data\" ";
    }
    else
    {
      os.setf(ios::scientific);
      os.precision( tecplotHeaderPrecision_);
      os << "T= \" ";
      for (std::vector<N_ANP_SweepParam>::const_iterator iterParam = outputManager_.getStepParamVec().begin(); iterParam != outputManager_.getStepParamVec().end(); ++iterParam)
      {
        os << " " << iterParam->name << " = " << iterParam->currentVal;
      }
      os << "\" ";
    }
    os << std::endl;

    os.setf(ios::scientific);
    os.precision(printParameters_.streamPrecision_);

    // put in the various sweep parameters as auxdata:
    if (!outputManager_.getStepParamVec().empty()) {
      for (std::vector<N_ANP_SweepParam>::const_iterator iterParam = outputManager_.getStepParamVec().begin(); iterParam != outputManager_.getStepParamVec().end(); ++iterParam)
      {
        // convert any ":" or "%" in the name to a "_", so as not to confuse tecplot.
        std::string tmpName(iterParam->name);
        replace(tmpName.begin(), tmpName.end(), '%', '_');
        replace(tmpName.begin(), tmpName.end(), ':', '_');
        os << "AUXDATA " << tmpName << " = " << "\" " << iterParam->currentVal << "\" ";
      }
      os << std::endl;
    }

    os.setf(ios::left, ios::adjustfield);

  } // procID

  ++headerPrintCalls_;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::outputPRINT_
// Purpose       : .PRINT output
// Special Notes :
// Scope         : private
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/10/00
//-----------------------------------------------------------------------------
void TimeTecPlot::doOutputTime(const N_LAS_Vector * solnVecPtr, const N_LAS_Vector * stateVecPtr, const N_LAS_Vector * storeVecPtr)
{
  outputManager_.setCurrentOutputter(this);
  double time = outputManager_.getCircuitTime();

  bool firstColPrinted = false;

  if (outputManager_.getProcID() == 0)
  {
    if (firstTimePrint_) // Setup Output Stream and Print Out Header
    {
      doOpen();

      timeHeader();

      firstTimePrint_ = false;
    }
  } // procID

  outputManager_.getCommPtr()->barrier();

  int i = 1;
  for (ParameterList::const_iterator iterParam = printParameters_.variableList_.begin() ; iterParam != printParameters_.variableList_.end(); ++iterParam, ++i)
  {
    double result = outputManager_.getPrintValue(iterParam, solnVecPtr, stateVecPtr, storeVecPtr);
    if ((*iterParam).getSimContext() == TIME_VAR)
      result *= printParameters_.outputTimeScaleFactor_;

    if (outputManager_.getProcID() == 0)
    {
      (*outStreamPtr_) << result << " ";
    } // procID

    outputManager_.getCommPtr()->barrier();
  }

  if (outputManager_.getProcID() == 0)
  {
    (*outStreamPtr_) << std::endl;
  } // procID
}

void TimeTecPlot::doFinishOutput()
{
  if (outputManager_.getProcID() == 0)
  {
    if (outStreamPtr_)
    {
      headerPrintCalls_ = 0;

      if (!outputManager_.getSTEPEnabledFlag())
      {
        outputManager_.closeFile(outStreamPtr_);
        outStreamPtr_ = 0;
      }
    }
  } // procID

  firstTimePrint_ = true;
}

void
TimeTecPlot::doFinishOutputStep()
{
  // Deal with the *tecplot file:
  if (outStreamPtr_)
  {
    outputManager_.closeFile(outStreamPtr_);
    outStreamPtr_ = 0;
  }
}

FrequencyTecPlot::FrequencyTecPlot(OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_(),
    outStreamPtr_(0),
    stepCount_(0),
    firstTime_(true),
    index_(0)
{
  if (printParameters_.extension_.empty())
    printParameters_.extension_ = ".dat";

  fixupColumns(outputManager_, printParameters_);
}

FrequencyTecPlot::~FrequencyTecPlot()
{
  outputManager_.closeFile(outStreamPtr_);
}


void FrequencyTecPlot::doParse()
{
  outFilename_ = outputFilename(printParameters_, outputManager_.getNetListFilename());
}

void FrequencyTecPlot::doOpen()
{
  if (outputManager_.getProcID() == 0 && outStreamPtr_ == 0)
    outStreamPtr_ = outputManager_.openFile(outFilename_);
}

void FrequencyTecPlot::frequencyHeader()
{
  index_ = 0;

  // STD header
  // Freq Domain Headers
  tecplotFreqHeader(outputManager_, printParameters_, *outStreamPtr_, stepCount_);

  outStreamPtr_->setf(ios::scientific);
  outStreamPtr_->precision(printParameters_.streamPrecision_);
  outStreamPtr_->setf(ios::left, ios::adjustfield);

  ++stepCount_;
}

void
FrequencyTecPlot::doOutputFrequency(double frequency, const N_LAS_Vector *real_solution_vector, const N_LAS_Vector *imaginary_solution_vector)
{
  outputManager_.setCurrentOutputter(this);
  bool firstColPrinted = false;
  int index = 0;

  std::ostringstream ost;
  ost << outputManager_.getProcID();
  string pN(ost.str());

  if (outputManager_.getProcID() == 0)
  {
    if (firstTime_) //Setup Output Stream and Print Out Header
    {
      doOpen();

      frequencyHeader();
      firstTime_ = false;
    }
    else // If running tecplot, the header needs to be revised each time.
         // Don't know how yet.
    {

    }

    if (printParameters_.delimiter_ == "")
    {
      outStreamPtr_->width(printParameters_.streamWidth_);
    }

    // Output a TAB.
    if (printParameters_.delimiter_ != "")
    {
      (*outStreamPtr_) << printParameters_.delimiter_;
    }
  } // procID


  if (outputManager_.getProcID() == 0)
  {
    if (printParameters_.delimiter_ == "")
    {
      outStreamPtr_->width(printParameters_.streamWidth_);
    }
    else
    {
      outStreamPtr_->width(0);
      if (firstColPrinted) {
        (*outStreamPtr_) << printParameters_.delimiter_;
      }
    }

    (*outStreamPtr_) << frequency;
    firstColPrinted = true;

  }

  // periodic time-domain steady-state output
  {
    for (ParameterList::const_iterator iterParam = printParameters_.variableList_.begin(); iterParam != printParameters_.variableList_.end(); ++iterParam)
    {
      double varValue = outputManager_.getPrintValue(iterParam, real_solution_vector, NULL, NULL, imaginary_solution_vector);

      if (outputManager_.getProcID() == 0)
      {
        (*outStreamPtr_) << varValue << " ";
      }
      outputManager_.getCommPtr()->barrier();

    } // end of output variable loop.
  }

  if (outputManager_.getProcID() == 0)
  {
    (*outStreamPtr_) << std::endl;
  }
}

void FrequencyTecPlot::doFinishOutput()
{
  firstTime_ = true;
}

void
FrequencyTecPlot::doFinishOutputStep()
{
}

TimeProbe::TimeProbe(OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    firstTimePrint_(true),
    printCount_(0),
    outFilename_(),
    suffix_(),
    outStreamPtr_(0),
    headerPrintCalls_(0)
{
  if (printParameters_.extension_.empty())
    printParameters_.extension_ = ".csd";

  fixupColumns(outputManager_, printParameters_);
}

TimeProbe::~TimeProbe()
{}


void TimeProbe::doParse()
{
  outFilename_ = outputFilename(printParameters_, outputManager_.getNetListFilename());
}

void TimeProbe::doOpen()
{
  if (outputManager_.getProcID() == 0 && outStreamPtr_ == 0)
  {
    outStreamPtr_ = outputManager_.openFile(outFilename_);
  }
}

void TimeProbe::timeHeader()
{
  std::ostream &os = *outStreamPtr_;

  if (outputManager_.getProcID() == 0)
  {
    // count the number of output variables.
    printCount_ = 0;
    for (ParameterList::const_iterator iterParam2 = printParameters_.variableList_.begin() ; iterParam2 != printParameters_.variableList_.end(); ++iterParam2)
    {
      if (iterParam2->getSimContext() != UNDEFINED)
        ++printCount_;
    }

    os << "#H" << std::endl;
    os << "SOURCE='Xyce' VERSION='"
       << N_UTL_Version::getShortVersionString() << "'" << std::endl;
    os << "TITLE='* " << outputManager_.getNetListFilename() << "'" << std::endl;

    os.setf(ios::scientific);
    os.precision(0); // streamPrecision_);
    if (outputManager_.getStepParamVec().empty())
    {
      os << "SUBTITLE='Xyce data";
    }
    else
    {
      std::vector<N_ANP_SweepParam>::const_iterator iterParam;
      std::vector<N_ANP_SweepParam>::const_iterator firstParam=outputManager_.getStepParamVec().begin();
      std::vector<N_ANP_SweepParam>::const_iterator lastParam=outputManager_.getStepParamVec().end();

      os << "SUBTITLE='Step param";
      for (iterParam=firstParam;iterParam!=lastParam;++iterParam)
      {
        os << " " << iterParam->name << " = ";
        os << iterParam->currentVal;
      }
    }

    os << " ' " << std::endl;;

    // set the time/date stamp
    os << getTimeDateStamp();
    os.setf(ios::scientific);
    os.precision(printParameters_.streamPrecision_);
    os << "TEMPERATURE='" << outputManager_.getCircuitTemp();
    os << "'" << std::endl;

    if (printParameters_.printType_ == PrintType::TRAN)
      os <<
        "ANALYSIS='Transient Analysis' SERIALNO='12345'" <<  std::endl;
    else
      os << "ANALYSIS='DC Sweep' " <<
        "SERIALNO='12345'" <<  std::endl;

    os << "ALLVALUES='NO' COMPLEXVALUES='NO' " <<
      "NODES='" << printCount_ << "'" << std::endl;

    if (printParameters_.printType_ == PrintType::TRAN)
    {
      os << "SWEEPVAR='Time' SWEEPMODE='VAR_STEP'" <<
        std::endl;
    }
    else
    {
      os << "SWEEPVAR='";
      os << outputManager_.getPRINTDCname();
      os << "' SWEEPMODE=";
      if ((outputManager_.getDCParamVec().size() > 0) && (outputManager_.getDCParamVec()[0].type == "LIST"))
      {
        os << "'LIST'" << std::endl;
      }
      else
      {
        os << "'VAR_STEP'" << std::endl;
      }
    }

    if (printParameters_.printType_ == PrintType::TRAN)
    {
      os << "XBEGIN='" << outputManager_.getAnaIntPtr()->getInitialTime()
         << "' XEND='" << outputManager_.getAnaIntPtr()->getFinalTime() << "'" << std::endl;
    }
    else
    {
      os << "XBEGIN='" << outputManager_.getPRINTDCstart()
         << "' XEND='" << outputManager_.getPRINTDCstop() << "'" << std::endl;
    }

    os << "FORMAT='0 VOLTSorAMPS;EFLOAT : "
       << "NODEorBRANCH;NODE  '  " << std::endl;
    os << "DGTLDATA='NO'";

    int dcSize = outputManager_.getDCParamVec().size();
    if (dcSize > 1)
    {
      os << "  ";
      for (int idc=1;idc<dcSize;++idc)
      {
        os << "SWEEP" << idc+1 << "PARM='";
        os << outputManager_.getDCParamVec()[idc].name;
        os << "' ";
        os << "SWEEP" << idc+1 << "VALUE='";
        os.setf(ios::scientific);
        os.precision(printParameters_.streamPrecision_);
        os << outputManager_.getDCParamVec()[idc].currentVal;
        os << "' ";
        os << std::endl;
      }
    }
    else
    {
      os << std::endl;
    }
    os << "#N" << std::endl;

    int i = 0;
    for (ParameterList::const_iterator iterParam = printParameters_.variableList_.begin() ; iterParam != printParameters_.variableList_.end(); ++iterParam, ++i)
    {
      if (i > 18)
      {
        i = 0;
        os << std::endl;
      }
      os << "'" << (*iterParam).tag() << "' ";
    }
    if (i != 0)
      os << std::endl;

    os.setf(ios::scientific);
    os.precision(printParameters_.streamPrecision_);
    os.setf(ios::left, ios::adjustfield);

  } // procID

  ++headerPrintCalls_;
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::outputPRINT_
// Purpose       : .PRINT output
// Special Notes :
// Scope         : private
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/10/00
//-----------------------------------------------------------------------------
void TimeProbe::doOutputTime(const N_LAS_Vector * solnVecPtr, const N_LAS_Vector * stateVecPtr, const N_LAS_Vector * storeVecPtr)
{
  outputManager_.setCurrentOutputter(this);
  double time = outputManager_.getCircuitTime();

  if (outputManager_.getProcID() == 0)
  {

    if (firstTimePrint_) // Setup Output Stream and Print Out Header
    {
      doOpen();

      timeHeader();

      index_ = 0;

      firstTimePrint_ = false;
    }

    std::ostream &os = *outStreamPtr_;

    os.width( 0);
    if (printParameters_.printType_ == PrintType::TRAN)
      os << "#C " << time << " " << printCount_ << std::endl;
    else
      os << "#C " << outputManager_.getPRINTDCvalue() << " " << printCount_ << std::endl;
  } // procID

  outputManager_.getCommPtr()->barrier();

  std::ostream &os = *outStreamPtr_;

  int i = 1;
  for (ParameterList::const_iterator iterParam = printParameters_.variableList_.begin() ; iterParam != printParameters_.variableList_.end(); ++iterParam, ++i)
  {
    double result = outputManager_.getPrintValue(iterParam, solnVecPtr, stateVecPtr, storeVecPtr);
    if ((*iterParam).getSimContext() == TIME_VAR)
      result *= printParameters_.outputTimeScaleFactor_;

    if (outputManager_.getProcID() == 0)
    {
      outStreamPtr_->width(0);
      os << result << ":" << i << "   ";
      if ((i/5)*5 == i)os << std::endl;
    } // procID

    outputManager_.getCommPtr()->barrier();
  }

  if (outputManager_.getProcID() == 0)
  {
    if (((i - 1)/5)*5 != (i - 1))
      os << std::endl;
  }
}

void TimeProbe::doFinishOutput()
{
  if (outputManager_.getProcID() == 0)
  {

    if (outStreamPtr_)
    {
      // print the last line of the *probe file.
      (*outStreamPtr_) << "#;" << std::endl;
    }

    if (!outputManager_.getSTEPEnabledFlag())
    {
      outputManager_.closeFile(outStreamPtr_);
      outStreamPtr_ = 0;
    }
  } // procID

  firstTimePrint_ = true;
}

void
TimeProbe::doFinishOutputStep()
{
  outputManager_.closeFile(outStreamPtr_);
  outStreamPtr_ = 0;
}

FrequencyProbe::FrequencyProbe(OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    firstTimePrint_(true),
    printCount_(0),
    outFilename_(),
    suffix_(),
    outStreamPtr_(0),
    headerPrintCalls_(0)
{
  if (printParameters_.extension_.empty())
    printParameters_.extension_ = ".csd";

  fixupColumns(outputManager_, printParameters_);
}

FrequencyProbe::~FrequencyProbe()
{}


void FrequencyProbe::doParse()
{
  outFilename_ = outputFilename(printParameters_, outputManager_.getNetListFilename());
}

void FrequencyProbe::doOpen()
{
  if (outputManager_.getProcID() == 0 && outStreamPtr_ == 0)
  {
    outStreamPtr_ = outputManager_.openFile(outFilename_);
  }
}

void FrequencyProbe::frequencyHeader()
{
  std::ostream &os = *outStreamPtr_;

  if (outputManager_.getProcID() == 0)
  {
    // count the number of output variables.
    printCount_ = 0;
    for (ParameterList::const_iterator iterParam2 = printParameters_.variableList_.begin() ; iterParam2 != printParameters_.variableList_.end(); ++iterParam2)
    {
      if (iterParam2->getSimContext() != UNDEFINED)
        ++printCount_;
    }

    os << "#H" << std::endl;
    os << "SOURCE='Xyce' VERSION='"
       << N_UTL_Version::getShortVersionString() << "'" << std::endl;
    os << "TITLE='* " << outputManager_.getNetListFilename() << "'" << std::endl;

    os.setf(ios::scientific);
    os.precision(0); // streamPrecision_);
    if (outputManager_.getStepParamVec().empty())
    {
      os << "SUBTITLE='Xyce data";
    }
    else
    {
      std::vector<N_ANP_SweepParam>::const_iterator iterParam;
      std::vector<N_ANP_SweepParam>::const_iterator firstParam=outputManager_.getStepParamVec().begin();
      std::vector<N_ANP_SweepParam>::const_iterator lastParam=outputManager_.getStepParamVec().end();

      os << "SUBTITLE='Step param";
      for (iterParam=firstParam;iterParam!=lastParam;++iterParam)
      {
        os << " " << iterParam->name << " = ";
        os << iterParam->currentVal;
      }
    }

    os << " ' " << std::endl;;

    // set the time/date stamp
    os << getTimeDateStamp();
    os.setf(ios::scientific);
    os.precision(printParameters_.streamPrecision_);
    os << "TEMPERATURE='" << outputManager_.getCircuitTemp();
    os << "'" << std::endl;

    os << "ANALYSIS='AC Sweep' " <<
      "SERIALNO='12345'" <<  std::endl;

    os << "ALLVALUES='NO' COMPLEXVALUES='YES' " <<
      "NODES='" << printCount_ << "'" << std::endl;

    os << "SWEEPVAR='";
    string varName = outputManager_.getPRINTDCname();
    if( varName == "" )
    { 
      varName="FREQ";
    }
    os << varName;
    os << "' SWEEPMODE=";
    if ( (!outputManager_.getDCParamVec().empty()) && (outputManager_.getDCParamVec()[0].type == "LIST") )
    {
      os << "'LIST'" << std::endl;
    }
    else
    {
      os << "'VAR_STEP'" << std::endl;
    }


    os << "XBEGIN='" << outputManager_.getPRINTDCstart()
       << "' XEND='" << outputManager_.getPRINTDCstop() << "'" << std::endl;

    os << "FORMAT='0 VOLTSorAMPS;EFLOAT : "
       << "NODEorBRANCH;NODE  '  " << std::endl;
    os << "DGTLDATA='NO'";

    int dcSize = outputManager_.getDCParamVec().size();
    if (dcSize > 1)
    {
      os << "  ";
      for (int idc=1;idc<dcSize;++idc)
      {
        os << "SWEEP" << idc+1 << "PARM='";
        os << outputManager_.getDCParamVec()[idc].name;
        os << "' ";
        os << "SWEEP" << idc+1 << "VALUE='";
        os.setf(ios::scientific);
        os.precision(printParameters_.streamPrecision_);
        os << outputManager_.getDCParamVec()[idc].currentVal;
        os << "' ";
        os << std::endl;
      }
    }
    else
    {
      os << std::endl;
    }
    os << "#N" << std::endl;

    int i = 0;
    for (ParameterList::const_iterator iterParam = printParameters_.variableList_.begin() ; iterParam != printParameters_.variableList_.end(); ++iterParam, ++i)
    {
      os << "'" << (*iterParam).tag() << "' ";
      if (i > 3)
      {
        i = 0;
        os << std::endl;
      }
    }

    if (i != 0)
      os << std::endl;

    os << flush;

    os.setf(ios::scientific);
    os.precision(printParameters_.streamPrecision_);
    os.setf(ios::left, ios::adjustfield);

  } // procID

  ++headerPrintCalls_;
}


void
FrequencyProbe::doOutputFrequency(double frequency, const N_LAS_Vector *real_solution_vector, const N_LAS_Vector *imaginary_solution_vector)
{
  if (outputManager_.getProcID() == 0)
  {
    if (firstTimePrint_) // Setup Output Stream and Print Out Header
    {
      doOpen();

      frequencyHeader();

      index_ = 0;
      firstTimePrint_ = false;
    }

    std::ostream &os = *outStreamPtr_;

    os.width( 0);
    os << "#C " << frequency << " " << printCount_ << std::endl;
  } // procID

  outputManager_.getCommPtr()->barrier();

  std::ostream &os = *outStreamPtr_;

  int i = 1;
  for (ParameterList::const_iterator iterParam = printParameters_.variableList_.begin() ; iterParam != printParameters_.variableList_.end(); ++iterParam, ++i)
  {
    // if context type of the param is SOLUTION_VAR then it is complex and 
    // we need to get both the real and imaginary part for output.
    // Any other context type is just real and the imaginary part will be zero
    double varValue = outputManager_.getPrintValue(iterParam, real_solution_vector, NULL, NULL, imaginary_solution_vector);
    double varValueIm = 0.0;
    if( iterParam->getSimContext() == SOLUTION_VAR )
    {
      // note we fake out getPrintValue to return the imaginary component by passing it in as the real vector.
      // not a great idea, but getPrintValue() will be rewritten to return complex types later.
      varValueIm = outputManager_.getPrintValue(iterParam, imaginary_solution_vector, NULL, NULL, imaginary_solution_vector);
    }

    if (outputManager_.getProcID() == 0)
    {
      os << varValue << "/" << varValueIm << ":" << i << "   ";
      if ((i/5)*5 == i)os << std::endl;
    }

    outputManager_.getCommPtr()->barrier();
  }

  if (outputManager_.getProcID() == 0)
  {
    if (((i - 1)/5)*5 != (i - 1))
      os << std::endl;
  }
}

void FrequencyProbe::doFinishOutput()
{
  if (outputManager_.getProcID() == 0)
  {

    if (outStreamPtr_)
    {
      // print the last line of the *probe file.
      (*outStreamPtr_) << "#;" << std::endl;
    }

    if (!outputManager_.getSTEPEnabledFlag())
    {
      outputManager_.closeFile(outStreamPtr_);
      outStreamPtr_ = 0;
    }
  } // procID

  firstTimePrint_ = true;
}

void
FrequencyProbe::doFinishOutputStep()
{
  outputManager_.closeFile(outStreamPtr_);
  outStreamPtr_ = 0;
}

HBPrn::HBPrn(OutputMgr &output_manager, const PrintParameters &freq_print_parameters, const PrintParameters &time_print_parameters)
  : outputManager_(output_manager),
    freqPrintParameters_(freq_print_parameters),
    timePrintParameters_(time_print_parameters),
    stepCount_(0),
    index_(0),
    firstTimeHB_(true),
    timeFilename_(),
    freqFilename_(),
    timeStreamPtr_(0),
    freqStreamPtr_(0)
{
  if (timePrintParameters_.extension_.empty())
    timePrintParameters_.extension_ = ".HB.TD.prn";

  if (freqPrintParameters_.extension_.empty())
    freqPrintParameters_.extension_ = ".HB.FD.prn";

  fixupColumns(outputManager_, timePrintParameters_);
  fixupColumns(outputManager_, freqPrintParameters_);
}

HBPrn::~HBPrn()
{
  outputManager_.closeFile(timeStreamPtr_);
  outputManager_.closeFile(freqStreamPtr_);
}


void HBPrn::doParse()
{
  timeFilename_ = outputFilename(timePrintParameters_, outputManager_.getNetListFilename());
  freqFilename_ = outputFilename(freqPrintParameters_, outputManager_.getNetListFilename());
}

void HBPrn::doOpen()
{
  if (outputManager_.getProcID() == 0 && timeStreamPtr_ == 0)
    timeStreamPtr_ = outputManager_.openFile(timeFilename_);
  if (outputManager_.getProcID() == 0 && freqStreamPtr_ == 0)
    freqStreamPtr_ = outputManager_.openFile(freqFilename_);
}

void HBPrn::doOutputHB(
  const std::vector<double>&    timePoints,
  const std::vector<double>&    freqPoints,
  const N_LAS_BlockVector &     timeDomainSolnVec,
  const N_LAS_BlockVector &     freqDomainSolnVecReal,
  const N_LAS_BlockVector &     freqDomainSolnVecImaginary)
{
  outputManager_.setCurrentOutputter(this);
  int index = 0;

  int blockCount = timeDomainSolnVec.blockCount();

  if (outputManager_.getProcID() == 0)
  {

    if (firstTimeHB_) //Setup Output Stream and Print Out Header
    {
      doOpen();

      printHeader(*timeStreamPtr_, timePrintParameters_);
      printHeader(*freqStreamPtr_, freqPrintParameters_);
      index_ = 0;

      ++stepCount_;

      firstTimeHB_ = false;
    }
  } // procID

  std::ostream &freq_os = *freqStreamPtr_;
  std::ostream &time_os = *timeStreamPtr_;

  // Get the indices:

  // Loop over the time points of the N_LAS_BlockVecor:
  for (int iblock = 0; iblock < blockCount; ++iblock)
  {
    outputManager_.setCircuitTime(timePoints[iblock]);
    outputManager_.setCircuitFrequency(freqPoints[iblock]);

    N_LAS_Vector * solnVecPtr = &(timeDomainSolnVec.block(iblock));
    N_LAS_Vector * realVecPtr = &(freqDomainSolnVecReal.block(iblock));
    N_LAS_Vector * imagVecPtr = &(freqDomainSolnVecImaginary.block(iblock));

    { // periodic time-domain steady-state output
      int column_index = 0;
      for (ParameterList::const_iterator it = timePrintParameters_.variableList_.begin(); it != timePrintParameters_.variableList_.end(); ++it, ++column_index)
      {
        double result = outputManager_.getPrintValue(it, solnVecPtr);

        if (outputManager_.getProcID() == 0)
          printValue(time_os, timePrintParameters_.table_.columnList_[column_index], timePrintParameters_.delimiter_, column_index, result);

        outputManager_.getCommPtr()->barrier();
      } // end of output variable loop.
    } // periodic time-domain steady-state output

    { // Fourier coefficient output
      int column_index = 0;
      for (ParameterList::const_iterator it = freqPrintParameters_.variableList_.begin(); it != freqPrintParameters_.variableList_.end(); ++it, ++column_index)
      {
        // state and store vec are not available in this context, but we must
        // pass in both the real and imaginary vectors
        double varValue = outputManager_.getPrintValue(it, realVecPtr, NULL, NULL, imagVecPtr);

        if (outputManager_.getProcID() == 0)
          printValue(freq_os, freqPrintParameters_.table_.columnList_[column_index], freqPrintParameters_.delimiter_, column_index, varValue);

        outputManager_.getCommPtr()->barrier();
      }
    } // Fourier coefficient output

    if (outputManager_.getProcID() == 0) {
      freq_os << std::endl;
      time_os << std::endl;
    }

    ++index_;
  }
}

void HBPrn::doFinishOutput()
{
  if (timeStreamPtr_)
  {
    if (!outputManager_.getSTEPEnabledFlag() && outputManager_.getPrintEndOfSimulationLine())
    {
      (*timeStreamPtr_) << "End of Xyce(TM) Simulation" << std::endl;
    }
  }

  if (freqStreamPtr_)
  {
    if (!outputManager_.getSTEPEnabledFlag() && outputManager_.getPrintEndOfSimulationLine())
    {
      (*freqStreamPtr_) << "End of Xyce(TM) Simulation" << std::endl;
    }
  }

  if (!outputManager_.getSTEPEnabledFlag())
  {
    outputManager_.closeFile(timeStreamPtr_);
    timeStreamPtr_ = 0;
    outputManager_.closeFile(freqStreamPtr_);
    freqStreamPtr_ = 0;

    firstTimeHB_ = true;
  }
}

void
HBPrn::doFinishOutputStep()
{
  // Deal with the *prn file:
  if (timeStreamPtr_)
  {
    if (outputManager_.getFormat() != Format::PROBE && outputManager_.getPrintEndOfSimulationLine())
      (*timeStreamPtr_) << "End of Xyce(TM) Parameter Sweep" << std::endl;
  }
  if (freqStreamPtr_)
  {
    if (outputManager_.getFormat() != Format::PROBE && outputManager_.getPrintEndOfSimulationLine())
      (*freqStreamPtr_) << "End of Xyce(TM) Parameter Sweep" << std::endl;
  }

  outputManager_.closeFile(timeStreamPtr_);
  timeStreamPtr_ = 0;
  outputManager_.closeFile(freqStreamPtr_);
  freqStreamPtr_ = 0;
}

HBCSV::HBCSV(OutputMgr &output_manager, const PrintParameters &freq_print_parameters, const PrintParameters &time_print_parameters)
  : outputManager_(output_manager),
    freqPrintParameters_(freq_print_parameters),
    timePrintParameters_(time_print_parameters),
    stepCount_(0),
    index_(0),
    firstTimeHB_(true),
    timeFilename_(),
    freqFilename_(),
    timeStreamPtr_(0),
    freqStreamPtr_(0)
{
  if (timePrintParameters_.extension_.empty())
    timePrintParameters_.extension_ = ".HB.TD.csv";

  if (freqPrintParameters_.extension_.empty())
    freqPrintParameters_.extension_ = ".HB.FD.csv";

  fixupColumns(outputManager_, timePrintParameters_);
  fixupColumns(outputManager_, freqPrintParameters_);
}

HBCSV::~HBCSV()
{
  outputManager_.closeFile(timeStreamPtr_);
  outputManager_.closeFile(freqStreamPtr_);
}


void HBCSV::doParse()
{
  timeFilename_ = outputFilename(timePrintParameters_, outputManager_.getNetListFilename());
  freqFilename_ = outputFilename(freqPrintParameters_, outputManager_.getNetListFilename());
}

void HBCSV::doOpen()
{
  if (outputManager_.getProcID() == 0 && timeStreamPtr_ == 0)
    timeStreamPtr_ = outputManager_.openFile(timeFilename_);
  if (outputManager_.getProcID() == 0 && freqStreamPtr_ == 0)
    freqStreamPtr_ = outputManager_.openFile(freqFilename_);
}

void HBCSV::doOutputHB(
  const std::vector<double>&    timePoints,
  const std::vector<double>&    freqPoints,
  const N_LAS_BlockVector &     timeDomainSolnVec,
  const N_LAS_BlockVector &     freqDomainSolnVecReal,
  const N_LAS_BlockVector &     freqDomainSolnVecImaginary)
{
  outputManager_.setCurrentOutputter(this);
  int index = 0;

  int blockCount = timeDomainSolnVec.blockCount();

  if (outputManager_.getProcID() == 0)
  {

    if (firstTimeHB_) //Setup Output Stream and Print Out Header
    {
      doOpen();

      printHeader(*timeStreamPtr_, timePrintParameters_);
      printHeader(*freqStreamPtr_, freqPrintParameters_);
      index_ = 0;

      ++stepCount_;

      firstTimeHB_ = false;
    }
  } // procID

  std::ostream &freq_os = *freqStreamPtr_;
  std::ostream &time_os = *timeStreamPtr_;

  // Get the indices:

  // Loop over the time points of the N_LAS_BlockVecor:
  for (int iblock = 0; iblock < blockCount; ++iblock)
  {
    outputManager_.setCircuitTime(timePoints[iblock]);
    outputManager_.setCircuitFrequency(freqPoints[iblock]);

    N_LAS_Vector * solnVecPtr = &(timeDomainSolnVec.block(iblock));
    N_LAS_Vector * realVecPtr = &(freqDomainSolnVecReal.block(iblock));
    N_LAS_Vector * imagVecPtr = &(freqDomainSolnVecImaginary.block(iblock));

    { // periodic time-domain steady-state output
      int column_index = 0;
      for (ParameterList::const_iterator it = timePrintParameters_.variableList_.begin(); it != timePrintParameters_.variableList_.end(); ++it, ++column_index)
      {
        double result = outputManager_.getPrintValue(it, solnVecPtr);

        if (outputManager_.getProcID() == 0)
          printValue(time_os, timePrintParameters_.table_.columnList_[column_index], timePrintParameters_.delimiter_, column_index, result);

        outputManager_.getCommPtr()->barrier();
      } // end of output variable loop.
    } // periodic time-domain steady-state output

    { // Fourier coefficient output
      int column_index = 0;
      for (ParameterList::const_iterator it = freqPrintParameters_.variableList_.begin(); it != freqPrintParameters_.variableList_.end(); ++it, ++column_index)
      {
        // state and store vec are not available in this context, but we must
        // pass in both the real and imaginary vectors
        double varValue = outputManager_.getPrintValue(it, realVecPtr, NULL, NULL, imagVecPtr);

        if (outputManager_.getProcID() == 0)
          printValue(freq_os, freqPrintParameters_.table_.columnList_[column_index], freqPrintParameters_.delimiter_, column_index, varValue);

        outputManager_.getCommPtr()->barrier();
      }
    } // Fourier coefficient output

    if (outputManager_.getProcID() == 0) {
      freq_os << std::endl;
      time_os << std::endl;
    }

    ++index_;
  }
}

void HBCSV::doFinishOutput()
{
  if (timeStreamPtr_)
  {
    if (!outputManager_.getSTEPEnabledFlag() && outputManager_.getPrintEndOfSimulationLine())
    {
      (*timeStreamPtr_) << "End of Xyce(TM) Simulation" << std::endl;
    }
  }

  if (freqStreamPtr_)
  {
    if (!outputManager_.getSTEPEnabledFlag() && outputManager_.getPrintEndOfSimulationLine())
    {
      (*freqStreamPtr_) << "End of Xyce(TM) Simulation" << std::endl;
    }
  }

  if (!outputManager_.getSTEPEnabledFlag())
  {
    outputManager_.closeFile(timeStreamPtr_);
    timeStreamPtr_ = 0;
    outputManager_.closeFile(freqStreamPtr_);
    freqStreamPtr_ = 0;

    firstTimeHB_ = true;
  }
}

void
HBCSV::doFinishOutputStep()
{
  // Deal with the *prn file:
  if (timeStreamPtr_)
  {
    if (outputManager_.getFormat() != Format::PROBE && outputManager_.getPrintEndOfSimulationLine())
      (*timeStreamPtr_) << "End of Xyce(TM) Parameter Sweep" << std::endl;
  }
  if (freqStreamPtr_)
  {
    if (outputManager_.getFormat() != Format::PROBE && outputManager_.getPrintEndOfSimulationLine())
      (*freqStreamPtr_) << "End of Xyce(TM) Parameter Sweep" << std::endl;
  }

  outputManager_.closeFile(timeStreamPtr_);
  timeStreamPtr_ = 0;

  outputManager_.closeFile(freqStreamPtr_);
  freqStreamPtr_ = 0;
}

HBTecPlot::HBTecPlot(OutputMgr &output_manager, const PrintParameters &freq_print_parameters, const PrintParameters &time_print_parameters)
  : outputManager_(output_manager),
    freqPrintParameters_(freq_print_parameters),
    timePrintParameters_(time_print_parameters),
    stepCount_(0),
    index_(0),
    firstTimeHB_(true),
    timeFilename_(),
    freqFilename_(),
    timeStreamPtr_(0),
    freqStreamPtr_(0)
{
  if (timePrintParameters_.extension_.empty())
    timePrintParameters_.extension_ = ".HB.TD.dat";

  if (freqPrintParameters_.extension_.empty())
    freqPrintParameters_.extension_ = ".HB.FD.dat";

  fixupColumns(outputManager_, timePrintParameters_);
  fixupColumns(outputManager_, freqPrintParameters_);
}

HBTecPlot::~HBTecPlot()
{
  outputManager_.closeFile(timeStreamPtr_);
  outputManager_.closeFile(freqStreamPtr_);
}


void HBTecPlot::doParse()
{
  timeFilename_ = outputFilename(timePrintParameters_, outputManager_.getNetListFilename());
  freqFilename_ = outputFilename(freqPrintParameters_, outputManager_.getNetListFilename());
}

void HBTecPlot::doOpen()
{
  if (outputManager_.getProcID() == 0 && timeStreamPtr_ == 0)
    timeStreamPtr_ = outputManager_.openFile(timeFilename_);
  if (outputManager_.getProcID() == 0 && freqStreamPtr_ == 0)
    freqStreamPtr_ = outputManager_.openFile(freqFilename_);
}

void HBTecPlot::doOutputHB(
  const std::vector<double>&    timePoints,
  const std::vector<double>&    freqPoints,
  const N_LAS_BlockVector &     timeDomainSolnVec,
  const N_LAS_BlockVector &     freqDomainSolnVecReal,
  const N_LAS_BlockVector &     freqDomainSolnVecImaginary)
{
  int index = 0;

  int blockCount = timeDomainSolnVec.blockCount();

  std::ostringstream ost;
  ost << outputManager_.getProcID();
  string pN(ost.str());

  if (outputManager_.getProcID() == 0)
  {

    if (firstTimeHB_) //Setup Output Stream and Print Out Header
    {
      doOpen();

      tecplotTimeHBHeader(*timeStreamPtr_);
      tecplotFreqHeader(outputManager_, freqPrintParameters_, *freqStreamPtr_, stepCount_);

      ++stepCount_;

      firstTimeHB_ = false;
    }
  } // procID

  std::ostream &freq_os = *freqStreamPtr_;
  std::ostream &time_os = *timeStreamPtr_;

  // Loop over the time points of the N_LAS_BlockVecor:
  for (int iblock = 0; iblock < blockCount; ++iblock)
  {
    outputManager_.setCircuitTime(timePoints[iblock]);
    outputManager_.setCircuitFrequency(freqPoints[iblock]);

    N_LAS_Vector *solnVecPtr = &(timeDomainSolnVec.block(iblock));
    N_LAS_Vector *realVecPtr = &(freqDomainSolnVecReal.block(iblock));
    N_LAS_Vector *imagVecPtr = &(freqDomainSolnVecImaginary.block(iblock));

    { // periodic time-domain steady-state output
      int column_index = 0;
      for (ParameterList::const_iterator it = timePrintParameters_.variableList_.begin(); it != timePrintParameters_.variableList_.end(); ++it, ++column_index)
      {
        double result = outputManager_.getPrintValue(it, solnVecPtr);

        if (outputManager_.getProcID() == 0)
          printValue(time_os, timePrintParameters_.table_.columnList_[column_index], timePrintParameters_.delimiter_, column_index, result);

        outputManager_.getCommPtr()->barrier();
      } // end of output variable loop.
    } // periodic time-domain steady-state output

    { // Fourier coefficient output
      int column_index = 0;
      for (ParameterList::const_iterator it = freqPrintParameters_.variableList_.begin(); it != freqPrintParameters_.variableList_.end(); ++it, ++column_index)
      {
        double varValue = outputManager_.getPrintValue(it, realVecPtr, NULL, NULL, imagVecPtr);
        if (outputManager_.getProcID() == 0)
          printValue(freq_os, freqPrintParameters_.table_.columnList_[column_index], freqPrintParameters_.delimiter_, column_index, varValue);

        outputManager_.getCommPtr()->barrier();
      }
    } // Fourier coefficient output

    if (outputManager_.getProcID() == 0) {
      freq_os << std::endl;
      time_os << std::endl;
    }

    ++index_;
  } // time scale loop.
}

void HBTecPlot::doFinishOutput()
{
  if (timeStreamPtr_)
  {
    if (!outputManager_.getSTEPEnabledFlag() && outputManager_.getPrintEndOfSimulationLine())
    {
      (*timeStreamPtr_) << "End of Xyce(TM) Simulation" << std::endl;
    }
  }

  if (freqStreamPtr_)
  {
    if (!outputManager_.getSTEPEnabledFlag() && outputManager_.getPrintEndOfSimulationLine())
    {
      (*freqStreamPtr_) << "End of Xyce(TM) Simulation" << std::endl;
    }
  }

  if (!outputManager_.getSTEPEnabledFlag())
  {

    outputManager_.closeFile(timeStreamPtr_);
    timeStreamPtr_ = 0;
    outputManager_.closeFile(freqStreamPtr_);
    freqStreamPtr_ = 0;

    firstTimeHB_ = true;
  }
}

void
HBTecPlot::doFinishOutputStep()
{
  // Deal with the *dat file:
  if (timeStreamPtr_)
  {
    if (outputManager_.getFormat() != Format::PROBE && outputManager_.getPrintEndOfSimulationLine())
      (*timeStreamPtr_) << "End of Xyce(TM) Parameter Sweep" << std::endl;
  }
  if (freqStreamPtr_)
  {
    if (outputManager_.getFormat() != Format::PROBE && outputManager_.getPrintEndOfSimulationLine())
      (*freqStreamPtr_) << "End of Xyce(TM) Parameter Sweep" << std::endl;
  }

  outputManager_.closeFile(timeStreamPtr_);
  timeStreamPtr_ = 0;
  outputManager_.closeFile(freqStreamPtr_);
  freqStreamPtr_ = 0;

}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::tecplotTimeHBHeader_
// Purpose       : header for tecplot. Time Domain HB(default)
// Special Notes :
// Scope         : private
// Creator       : Todd Coffey, Ting Mei
// Creation Date : 7/31/08
//-----------------------------------------------------------------------------
void HBTecPlot::tecplotTimeHBHeader( ostream & stream)
{
  if (stepCount_ == 0)
  {
    stream << " TITLE                            = \" Xyce Time Domain HB data, " << outputManager_.getNetListFilename() << "\", " << std::endl;
    stream << "\tVARIABLES                       = \"TIME \" " << std::endl;

    // output the user-specified solution vars:
    for (ParameterList::const_iterator iterParam = timePrintParameters_.variableList_.begin() ; iterParam != timePrintParameters_.variableList_.end(); ++iterParam)
    {
      stream << "\" " << (*iterParam).tag() << "\" " << std::endl;
    }
  }

  // output some AUXDATA
  stream << "DATASETAUXDATA ";
  stream << getTecplotTimeDateStamp();
  stream << std::endl;
  stream << "ZONE F=POINT\n";

  if (outputManager_.getStepParamVec().empty())
  {
    stream << " T=\"Xyce data\" ";
  }
  else
  {
    std::vector<N_ANP_SweepParam>::const_iterator iterParam;
    std::vector<N_ANP_SweepParam>::const_iterator firstParam=outputManager_.getStepParamVec().begin();
    std::vector<N_ANP_SweepParam>::const_iterator lastParam=outputManager_.getStepParamVec().end();

    stream << " T= \" ";
    for (iterParam=firstParam;iterParam!=lastParam;++iterParam)
    {
      int tecplotHeaderPrecision_ = 2;
      stream.setf(ios::scientific);
      stream.precision( tecplotHeaderPrecision_);
      stream << " " << iterParam->name << " = ";
      stream << iterParam->currentVal;
    }
    stream << "\" ";
  }

  stream << std::endl;
  stream << flush;
}

MPDEPrn::MPDEPrn(OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_(),
    outStreamPtr_(0),
    stepCount_(0),
    firstTimeMPDE_(true),
    n1_(0),
    n2_(0)
{
  if (printParameters_.extension_.empty())
    printParameters_.extension_ = ".MPDE.prn";

  printParameters_.table_.addColumn("TIME1", printParameters_.streamWidth_, printParameters_.streamPrecision_, Table::JUSTIFICATION_RIGHT);
  printParameters_.table_.addColumn("TIME2", printParameters_.streamWidth_, printParameters_.streamPrecision_, Table::JUSTIFICATION_RIGHT);

  fixupColumns(outputManager_, printParameters_);
}

MPDEPrn::~MPDEPrn()
{
  outputManager_.closeFile(outStreamPtr_);
}


void MPDEPrn::doParse()
{
  // if (outputManager_.getSTEPEnabledFlag())
  // {
  //   std::ostringstream num;
  //   num << stepCount_;
  //   outFilename_ = outputManager_.getNetListFilename() + ".STEP" + num.str() + ".prn";
  //   ++stepCount_;
  // }
  // else
  // {
    outFilename_ = outputFilename(printParameters_, outputManager_.getNetListFilename());
  // }
}

void MPDEPrn::doOpen()
{
  if (outputManager_.getProcID() == 0 && outStreamPtr_ == 0)
    outStreamPtr_ = outputManager_.openFile(outFilename_);
}

void MPDEPrn::mpdeHeader()
{}

void MPDEPrn::doOutputTime(const N_LAS_Vector *solution_vector, const N_LAS_Vector *state_vector, const N_LAS_Vector *store_vector)
{}

void MPDEPrn::doOutputMPDE(double time, const N_LAS_Vector *solution_vector)
{
  outputManager_.setCurrentOutputter(this);

  const N_LAS_BlockVector * blockSolVecPtr = dynamic_cast<const N_LAS_BlockVector*>(solution_vector);
  if (blockSolVecPtr == 0)
    return;

  int blockCount = blockSolVecPtr->blockCount();
  n2_ = blockCount; // fast time points.
  ++n1_;            // slow time points.  increments by one each call.

  // get the fast time points:
  const std::vector<double> & fastTimes = outputManager_.getMpdeMgrPtr()->getFastTimePoints();

  if (outputManager_.getProcID() == 0)
  {
    if (firstTimeMPDE_) //Setup Output Stream and Print Out Header
    {
      doOpen();

      mpdeHeader();

      firstTimeMPDE_ = false;
    }
  } // procID

  std::ostream &os = *outStreamPtr_;
  
  // Loop over the fast time points of the N_LAS_BlockVecor:
  for (int iblock=0;iblock<n2_+1;++iblock)
  {
    N_LAS_Vector * solnVecPtr;

    if (iblock == n2_)
    {
      solnVecPtr = &(blockSolVecPtr->block(0));
    }
    else
    {
      solnVecPtr = &(blockSolVecPtr->block(iblock));
    }

    if (outputManager_.getProcID() == 0)
    {
      //-------------------------------------
      // Get the 2 time values first.
      //-------------------------------------
      double first  = time;
      double second = fastTimes[iblock];

      // time 1:
      printValue(os, printParameters_.table_.columnList_[0], printParameters_.delimiter_, 0, first);

      // time 2:
      printValue(os, printParameters_.table_.columnList_[1], printParameters_.delimiter_, 1, second);
    }

    int column_index = 2;
    for (ParameterList::const_iterator iterParam = printParameters_.variableList_.begin(); iterParam != printParameters_.variableList_.end(); ++iterParam, ++column_index)
    {
      double result = outputManager_.getPrintValue(iterParam, solnVecPtr);

      if (outputManager_.getProcID() == 0)
      {
        printValue(os, printParameters_.table_.columnList_[column_index], printParameters_.delimiter_, column_index, result);
      }

      outputManager_.getCommPtr()->barrier();

    }

    if (outputManager_.getProcID() == 0)
    {
      (*outStreamPtr_) << std::endl;
    }
  } // fast time scale loop.
}

void MPDEPrn::doOutputHomotopy(const std::vector<std::string> & parameter_names, const std::vector<double> & parameter_values, const N_LAS_Vector * solution_vector)
{}

void MPDEPrn::doOutputMORTF(bool origSystem, const double & freq, const Teuchos::SerialDenseMatrix<int, std::complex<double> >& H)
{}

void MPDEPrn::doOutputHB(
  const std::vector<double>&    timePoints,
  const std::vector<double>&    freqPoints,
  const N_LAS_BlockVector &     timeDomainSolnVec,
  const N_LAS_BlockVector &     freqDomainSolnVecReal,
  const N_LAS_BlockVector &     freqDomainSolnVecImaginary)
{
}

void
MPDEPrn::doOutputFrequency(double frequency, const N_LAS_Vector *real_solution_vector, const N_LAS_Vector *imaginary_solution_vector)
{}


void MPDEPrn::doResetOutput()
{}

void MPDEPrn::doFinishOutput()
{
  outputManager_.closeFile(outStreamPtr_);
  outStreamPtr_ = 0;
}

void
MPDEPrn::doFinishOutputStep()
{
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::stdTimeMPDEHeader_
// Purpose       : header for std. Time Domain MPDE(default)
// Special Notes :
// Scope         : private
// Creator       : Todd Coffey, Ting Mei
// Creation Date : 7/31/08
//-----------------------------------------------------------------------------
void MPDEPrn::stdTimeMPDEHeader( std::ostream & stream)
{
}

MPDETecPlot::MPDETecPlot(OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_(),
    outStreamPtr_(0),
    stepCount_(0),
    firstTimeMPDE_(true),
    n1_(0),
    n2_(0)
{
}

MPDETecPlot::~MPDETecPlot()
{
  outputManager_.closeFile(outStreamPtr_);
}


void MPDETecPlot::doParse()
{
  ///////////////////////////////////////////////////////////////////////
  if (outputManager_.getSTEPEnabledFlag())
  {
    std::ostringstream num;
    num << stepCount_;
    outFilename_ = outputManager_.getNetListFilename() + ".STEP" + num.str() + ".MPDE.prn";
    ++stepCount_;
  }
  else
  {
    outFilename_ = outputManager_.getNetListFilename() + ".MPDE.prn";
  }
  ///////////////////////////////////////////////////////////////////////
}

void MPDETecPlot::doOpen()
{
  if (outputManager_.getProcID() == 0 && outStreamPtr_ == 0)
    outStreamPtr_ = outputManager_.openFile(outFilename_);
}

void MPDETecPlot::doOutputHeader()
{
  ParameterList::const_iterator iterParam = printParameters_.variableList_.begin();
  ParameterList::const_iterator iterParam_end = printParameters_.variableList_.end();

  while ( iterParam != iterParam_end &&
          iterParam->tag() != "I" &&
          iterParam->tag() != "V" &&
          iterParam->tag() != "VSTART" &&
          !(iterParam->hasExpressionTag()))
  {
    ++iterParam;
  }

  (*outStreamPtr_) << " TITLE = \" Xyce MPDE data, " << outputManager_.getNetListFilename() << "\", " << std::endl;
  (*outStreamPtr_) << "\tVARIABLES = \"T1(sec) \", \"T2(sec)\", \n";

  // output the user-specified solution vars:
  while (iterParam != printParameters_.variableList_.end())
  {
    (*outStreamPtr_) << "\" ";
    (*outStreamPtr_) << (*iterParam).tag();

    ++iterParam;
    (*outStreamPtr_) << "\" " << std::endl;
  }
  // output some AUXDATA
  (*outStreamPtr_) << "DATASETAUXDATA ";
  (*outStreamPtr_) << getTecplotTimeDateStamp();
  (*outStreamPtr_) << std::endl;

  (*outStreamPtr_) << "ZONE I=" << n2_ + 1<<", ";
  (*outStreamPtr_) << " J=" << n1_ << ", ";
  (*outStreamPtr_) << " F=POINT\n";

  (*outStreamPtr_) << std::endl;

  (*outStreamPtr_) << flush;

  outStreamPtr_->setf(ios::scientific);
  outStreamPtr_->precision(80); // streamPrecision_);
  outStreamPtr_->setf(ios::left, ios::adjustfield);
}

void MPDETecPlot::doOutputTime(const N_LAS_Vector *solution_vector, const N_LAS_Vector *state_vector, const N_LAS_Vector *store_vector)
{}

void MPDETecPlot::doOutputMPDE(double time, const N_LAS_Vector *solution_vector)
{
  outputManager_.setCurrentOutputter(this);
  const N_LAS_BlockVector * blockSolVecPtr = dynamic_cast<const N_LAS_BlockVector*>(solution_vector);

  int blockCount = blockSolVecPtr->blockCount();
  n2_ = blockCount; // fast time points.
  ++n1_;            // slow time points.  increments by one each call.

  // get the fast time points:
  const std::vector<double> & fastTimes = outputManager_.getMpdeMgrPtr()->getFastTimePoints();

  std::ostringstream ost;
  ost << outputManager_.getProcID();
  string pN(ost.str());

  if (outputManager_.getProcID() == 0)
  {
    if (firstTimeMPDE_) //Setup Output Stream and Print Out Header
    {
      doOpen();

      doOutputHeader();

      firstTimeMPDE_ = false;
    }
    else // If running tecplot, the header needs to be revised each time.
      // Don't know how yet.
    {

    }

    if (printParameters_.delimiter_ == "")
      outStreamPtr_->width(printParameters_.streamWidth_);

    // Output a TAB.
    if (printParameters_.delimiter_ != "")
      (*outStreamPtr_) << printParameters_.delimiter_;
  } // procID

  ParameterList::const_iterator iterParam = printParameters_.variableList_.begin();
  ParameterList::const_iterator begin_tpL;
  ParameterList::const_iterator last = printParameters_.variableList_.end();

  while ( iterParam != last &&
          iterParam->tag() != "I" &&
          iterParam->tag() != "V" &&
          iterParam->tag() != "VSTART" &&
          !(iterParam->hasExpressionTag()))
  {
    ++iterParam;
  }
  begin_tpL = iterParam;

  // Get the indices:

  // Loop over the fast time points of the N_LAS_BlockVecor:
  for (int iblock=0;iblock<n2_+1;++iblock)
  {
    N_LAS_Vector * solnVecPtr;

    if (iblock == n2_)
    {
      solnVecPtr = &(blockSolVecPtr->block(0));
    }
    else
    {
      solnVecPtr = &(blockSolVecPtr->block(iblock));
    }

    if (outputManager_.getProcID() == 0)
    {
      //-------------------------------------
      // Get the 2 time values first.
      //-------------------------------------
      double first  = 0.0;
      double second = 0.0;

      second = fastTimes[iblock];
      first  = time;

      // time 1:
      if (printParameters_.delimiter_ == "")
        outStreamPtr_->width(printParameters_.streamWidth_);
      else
      {
        outStreamPtr_->width(0);
        if (printParameters_.delimiter_ != "")
          (*outStreamPtr_) << printParameters_.delimiter_;
      }

      (*outStreamPtr_) << first;

      // time 2:
      if (printParameters_.delimiter_ == "")
        outStreamPtr_->width(printParameters_.streamWidth_);
      else
      {
        outStreamPtr_->width(0);
        if (printParameters_.delimiter_ != "")
          (*outStreamPtr_) << printParameters_.delimiter_;
      }

      (*outStreamPtr_) << second;
    }
//
    int i;
    for (i = 1, iterParam=begin_tpL ; iterParam != last; ++iterParam, ++i)
    {
      double result;
      result = 0.0;
      result = outputManager_.getPrintValue(iterParam, solnVecPtr);

      if (outputManager_.getProcID() == 0)
      {
        if (printParameters_.delimiter_ == "")
          outStreamPtr_->width(printParameters_.streamWidth_);
        else
        {
          outStreamPtr_->width(0);
          if (printParameters_.delimiter_ != "")
            (*outStreamPtr_) << printParameters_.delimiter_;
        }
        (*outStreamPtr_) << result;
      }

      outputManager_.getCommPtr()->barrier();

    } // end of output variable loop.

    if (outputManager_.getProcID() == 0)
    {
      (*outStreamPtr_) << "\n";
    }

  } // fast time scale loop.

  if (outputManager_.getProcID() == 0)
  {
    (*outStreamPtr_) << std::endl << flush;
  }

}

void MPDETecPlot::doOutputHomotopy(const std::vector<std::string> & parameter_names, const std::vector<double> & parameter_values, const N_LAS_Vector * solution_vector)
{}

void MPDETecPlot::doOutputMORTF(bool origSystem, const double & freq, const Teuchos::SerialDenseMatrix<int, std::complex<double> >& H)
{}

void MPDETecPlot::doOutputHB(
  const std::vector<double>&    timePoints,
  const std::vector<double>&    freqPoints,
  const N_LAS_BlockVector &     timeDomainSolnVec,
  const N_LAS_BlockVector &     freqDomainSolnVecReal,
  const N_LAS_BlockVector &     freqDomainSolnVecImaginary)
{}

void
MPDETecPlot::doOutputFrequency(double frequency, const N_LAS_Vector *real_solution_vector, const N_LAS_Vector *imaginary_solution_vector)
{}

void MPDETecPlot::doResetOutput()
{}

void MPDETecPlot::doFinishOutput()
{
  outputManager_.closeFile(outStreamPtr_);
  outStreamPtr_= 0;
}

void
MPDETecPlot::doFinishOutputStep()
{
}

//-----------------------------------------------------------------------------
// Function      : OutputMgr::stdTimeMPDEHeader_
// Purpose       : header for std. Time Domain MPDE(default)
// Special Notes :
// Scope         : private
// Creator       : Todd Coffey, Ting Mei
// Creation Date : 7/31/08
//-----------------------------------------------------------------------------
void MPDETecPlot::stdTimeMPDEHeader( std::ostream & stream)
{
}


HomotopyPrn::HomotopyPrn(OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_(),
    outStreamPtr_(0),
    stepCount_(0),
    index_(0),
    printCount_(0),
    firstTimeHomotopy_(true)
{
  fixupColumns(outputManager_, printParameters_);
}

HomotopyPrn::~HomotopyPrn()
{
  outputManager_.closeFile(outStreamPtr_);
}


void HomotopyPrn::doParse()
{
  outFilename_ = outputManager_.getNetListFilename() + ".HOMOTOPY.prn";
}

void HomotopyPrn::doOpen()
{
  if (outputManager_.getProcID() == 0 && outStreamPtr_ == 0)
  {
    outStreamPtr_ = outputManager_.openFile(outFilename_);
  }
}

void HomotopyPrn::homotopyHeader(const std::vector<std::string> & parameter_names, const std::vector<double> & param_values, const N_LAS_Vector * solution_vector)
{
  std::ostream &os = *outStreamPtr_;

  if (columnList_.empty()) {
    Table::Justification justification = printParameters_.delimiter_.empty() ? Table::JUSTIFICATION_CENTER :  Table::JUSTIFICATION_NONE;

    for (std::vector<std::string>::const_iterator it = parameter_names.begin(); it != parameter_names.end(); ++it)
      columnList_.push_back(Table::Column((*it), std::ios_base::scientific, printParameters_.streamWidth_, printParameters_.streamPrecision_, justification));
  }

  index_ = 0;

  if (stepCount_ == 0)
  {
    int column_index = 0;
    for (Table::ColumnList::const_iterator it = printParameters_.table_.columnList_.begin(); it != printParameters_.table_.columnList_.end(); ++it, ++column_index) {
      if (it != printParameters_.table_.columnList_.begin())
        os << (printParameters_.delimiter_.empty() ? " " : printParameters_.delimiter_);

      if (column_index == 1) {
        for (Table::ColumnList::const_iterator it2 = columnList_.begin(); it2 != columnList_.end(); ++it2) {
          if (it2 != columnList_.begin())
            os << printParameters_.delimiter_;
          printHeader(os, (*it2));
        }
      }

      printHeader(os, (*it));
    }

    os << std::endl;
  }

  ++stepCount_;
}

void HomotopyPrn::doOutputHomotopy(const std::vector<std::string> & parameter_names, const std::vector<double> & parameter_values, const N_LAS_Vector * solution_vector)
{
  outputManager_.setCurrentOutputter(this);

  double tmpTime = outputManager_.getAnaIntPtr()->getTime();

  if (outputManager_.getProcID() == 0)
  {

    if (firstTimeHomotopy_) //Setup Output Stream and Print Out Header
    {
      doOpen();

      homotopyHeader(parameter_names, parameter_values, solution_vector);

      firstTimeHomotopy_ = false;
    }
  }

  std::ostream &os = *outStreamPtr_;

  int column_index = 0;
  for (ParameterList::const_iterator it = printParameters_.variableList_.begin(); it != printParameters_.variableList_.end(); ++it, ++column_index)
  {
    double result = outputManager_.getPrintValue(it, solution_vector);

    if (outputManager_.getProcID() == 0) {
      if (column_index == 1)
        for (int i = 0; i < parameter_values.size(); ++i)
          printValue(os, columnList_[i], printParameters_.delimiter_, 1, parameter_values[i]);

      printValue(os, printParameters_.table_.columnList_[column_index], printParameters_.delimiter_, column_index, result);
    }

    outputManager_.getCommPtr()->barrier();
  }

  if (outputManager_.getProcID() == 0)
    os << std::endl;

  ++index_;
}

void HomotopyPrn::doResetOutput()
{}

void HomotopyPrn::doFinishOutput()
{
  firstTimeHomotopy_ = true;
}

void
HomotopyPrn::doFinishOutputStep()
{
  // close the homotopy file.
  if (outStreamPtr_)
  {
    (*outStreamPtr_) << "End of Xyce(TM) Homotopy Simulation" << std::endl;
  }

  outputManager_.closeFile(outStreamPtr_);
  outStreamPtr_ = 0;
}

HomotopyTecPlot::HomotopyTecPlot(OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_(),
    outStreamPtr_(0),
    stepCount_(0),
    index_(0),
    printCount_(0),
    firstTimeHomotopy_(true)
{
}

HomotopyTecPlot::~HomotopyTecPlot()
{
  outputManager_.closeFile(outStreamPtr_);
}


void HomotopyTecPlot::doParse()
{
  string finalSuffix = ".prn";
  if (outputManager_.getFormat() == Format::TECPLOT)
  {
    finalSuffix = ".dat";
  }
  else if (outputManager_.getFormat() == Format::CSV)
  {
    finalSuffix = ".csv";
  }
  else if (outputManager_.getFormat() == Format::PROBE)
  {
    finalSuffix = ".csd";
  }

  outFilename_ = outputManager_.getNetListFilename() + ".HOMOTOPY" + finalSuffix;
}

void HomotopyTecPlot::doOpen()
{
  if (outputManager_.getProcID() == 0 && outStreamPtr_ == 0)
  {
    outStreamPtr_ = outputManager_.openFile(outFilename_);
  }
}

void HomotopyTecPlot::doOutputHeader(const std::vector<std::string> & parameter_names, const std::vector<double> & param_values, const N_LAS_Vector * solution_vector)
{
  outputManager_.setCurrentOutputter(this);
  index_ = 0;

  ParameterList::const_iterator iterParam = printParameters_.variableList_.begin();
  ParameterList::const_iterator last = printParameters_.variableList_.end();


  while ( iterParam != last && iterParam->tag() != "I" &&
          iterParam->tag() != "V" &&
          iterParam->tag() != "VSTART" &&
          !(iterParam->hasExpressionTag()))
  {
    ++iterParam;
  }

  if (stepCount_ == 0)
  {
    (*outStreamPtr_) << " TITLE = \" Xyce homotopy data, " << outputManager_.getNetListFilename() << "\", " << std::endl;
    (*outStreamPtr_) << "\tVARIABLES = ";
    if (!outputManager_.getNoIndex())
    {
      (*outStreamPtr_) << "\"Index \" " << std::endl;
    }
    if (printParameters_.printType_ == PrintType::TRAN)
    {
      (*outStreamPtr_) << "\"TIME \" " << std::endl;
    }

    // output the continuation parameters:
    std::vector<std::string>::const_iterator iter_name;
    for (iter_name=parameter_names.begin(); iter_name!= parameter_names.end();
         ++iter_name)
    {
      (*outStreamPtr_) << "\" ";
      (*outStreamPtr_) << *iter_name;
      (*outStreamPtr_) << "\" " << std::endl;
    }

    // output the user-specified solution vars:
    while (iterParam != printParameters_.variableList_.end())
    {
      (*outStreamPtr_) << "\" ";
      (*outStreamPtr_) << (*iterParam).tag();
      (*outStreamPtr_) << "\" " << std::endl;
      ++iterParam;
    }
  }

  // output some AUXDATA
  (*outStreamPtr_) << "DATASETAUXDATA ";
  (*outStreamPtr_) << getTecplotTimeDateStamp();
  (*outStreamPtr_) << std::endl;

  (*outStreamPtr_) << "ZONE F=POINT";


  if (outputManager_.getStepParamVec().empty())
  {
    (*outStreamPtr_) << " T=\"Xyce data\" ";
  }
  else
  {
    std::vector<N_ANP_SweepParam>::const_iterator iterParam;
    std::vector<N_ANP_SweepParam>::const_iterator firstParam=outputManager_.getStepParamVec().begin();
    std::vector<N_ANP_SweepParam>::const_iterator lastParam=outputManager_.getStepParamVec().end();

    (*outStreamPtr_) << " T= \" ";
    for (iterParam=firstParam;iterParam!=lastParam;++iterParam)
    {
      int tecplotHeaderPrecision_ = 2;
      outStreamPtr_->setf(ios::scientific);
      outStreamPtr_->precision( tecplotHeaderPrecision_);
      (*outStreamPtr_) << " " << iterParam->name << " = ";
      (*outStreamPtr_) << iterParam->currentVal;
    }
    (*outStreamPtr_) << "\" ";
  }

  (*outStreamPtr_) << std::endl;
  (*outStreamPtr_) << flush;


  (*outStreamPtr_) << flush;

  outStreamPtr_->setf(ios::scientific);
  outStreamPtr_->precision(printParameters_.streamPrecision_);
  outStreamPtr_->setf(ios::left, ios::adjustfield);


  ++stepCount_;
}

void HomotopyTecPlot::doOutputTime(const N_LAS_Vector *solution_vector, const N_LAS_Vector *state_vector, const N_LAS_Vector *store_vector)
{}

void HomotopyTecPlot::doOutputMPDE(double time, const N_LAS_Vector *solution_vector)
{}

void HomotopyTecPlot::doOutputHomotopy(const std::vector<std::string> & parameter_names, const std::vector<double> & parameter_values, const N_LAS_Vector * solution_vector)
{
  outputManager_.setCurrentOutputter(this);
  std::ostringstream ost;
  ost << outputManager_.getProcID();
  string pN(ost.str());

  double tmpTime = outputManager_.getAnaIntPtr()->getTime();

  if (outputManager_.getProcID() == 0)
  {

    if (firstTimeHomotopy_) //Setup Output Stream and Print Out Header
    {
      doOpen();
      doOutputHeader(parameter_names, parameter_values, solution_vector);

      firstTimeHomotopy_ = false;
    }

    if (printParameters_.delimiter_ == "")
      outStreamPtr_->width(8);
    else
      outStreamPtr_->width(0);
    (*outStreamPtr_) << index_++;

    if (printParameters_.delimiter_ == "")
      outStreamPtr_->width(printParameters_.streamWidth_);

    if (printParameters_.printType_ == PrintType::TRAN)
    {
      // Output a TAB.
      if (printParameters_.delimiter_ != "")
        (*outStreamPtr_) << printParameters_.delimiter_;
      (*outStreamPtr_) << tmpTime;
    }

    //-------------------------------------
    //HOMOTOPY PARAM VALUE OUTPUT GOES HERE
    //-------------------------------------

    for (int iparam=0;iparam < parameter_values.size(); ++iparam)
    {

      if (printParameters_.delimiter_ == "")
        outStreamPtr_->width(printParameters_.streamWidth_);
      else
      {
        outStreamPtr_->width(0);
        if (printParameters_.delimiter_ != "")
          (*outStreamPtr_) << printParameters_.delimiter_;
      }

      (*outStreamPtr_) << parameter_values[iparam];
    }
  } // procID

  ParameterList::const_iterator iterParam = printParameters_.variableList_.begin();
  ParameterList::const_iterator last = printParameters_.variableList_.end();

  bool foundFirstParam=false;
  while ( iterParam != last &&  !foundFirstParam)
  {
    if (iterParam->getSimContext() != UNDEFINED)
    {
      // found the beginning so exit this loop
      foundFirstParam=true;
    }
    else
    {
      ++iterParam;
    }
  }

  int i;
  for (i = 1; iterParam != last; ++iterParam, ++i)
  {
    double result;
    result = 0.0;
    result = outputManager_.getPrintValue(iterParam, solution_vector);

    if (outputManager_.getProcID() == 0)
    {
      if (printParameters_.delimiter_ == "")
        outStreamPtr_->width(printParameters_.streamWidth_);
      else
      {
        outStreamPtr_->width(0);
        if (printParameters_.delimiter_ != "")
          (*outStreamPtr_) << printParameters_.delimiter_;
      }
      (*outStreamPtr_) << result;
    }

    outputManager_.getCommPtr()->barrier();
  }

  if (outputManager_.getProcID() == 0)
    (*outStreamPtr_) << std::endl;
}

void HomotopyTecPlot::doOutputMORTF(bool origSystem, const double & freq, const Teuchos::SerialDenseMatrix<int, std::complex<double> >& H)
{}

void HomotopyTecPlot::doOutputHB(
  const std::vector<double>&    timePoints,
  const std::vector<double>&    freqPoints,
  const N_LAS_BlockVector &     timeDomainSolnVec,
  const N_LAS_BlockVector &     freqDomainSolnVecReal,
  const N_LAS_BlockVector &     freqDomainSolnVecImaginary)
{
}

void
HomotopyTecPlot::doOutputFrequency(double frequency, const N_LAS_Vector *real_solution_vector, const N_LAS_Vector *imaginary_solution_vector)
{}


void HomotopyTecPlot::doResetOutput()
{}

void HomotopyTecPlot::doFinishOutput()
{
  firstTimeHomotopy_ = true;
}

void
HomotopyTecPlot::doFinishOutputStep()
{
  // close the homotopy file.
  if (outStreamPtr_)
  {
    (*outStreamPtr_) << "End of Xyce(TM) Homotopy Simulation" << std::endl;
  }

  outputManager_.closeFile(outStreamPtr_);
  outStreamPtr_ = 0;
}

HomotopyProbe::HomotopyProbe(OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_(),
    outStreamPtr_(0),
    stepCount_(0),
    index_(0),
    printCount_(0),
    firstTimeHomotopy_(true)
{
}

HomotopyProbe::~HomotopyProbe()
{
  outputManager_.closeFile(outStreamPtr_);
}


void HomotopyProbe::doParse()
{
  string finalSuffix = ".prn";
  if (outputManager_.getFormat() == Format::TECPLOT)
  {
    finalSuffix = ".dat";
  }
  else if (outputManager_.getFormat() == Format::CSV)
  {
    finalSuffix = ".csv";
  }
  else if (outputManager_.getFormat() == Format::PROBE)
  {
    finalSuffix = ".csd";
  }

  outFilename_ = outputManager_.getNetListFilename() + ".HOMOTOPY" + finalSuffix;
}

void HomotopyProbe::doOpen()
{
  if (outputManager_.getProcID() == 0 && outStreamPtr_ == 0)
  {
    outStreamPtr_ = outputManager_.openFile(outFilename_);
  }
}

void HomotopyProbe::doOutputHeader(const std::vector<std::string> & parameter_names, const std::vector<double> & param_values, const N_LAS_Vector * solution_vector)
{
  index_ = 0;

  ParameterList::const_iterator iterParam = printParameters_.variableList_.begin();
  ParameterList::const_iterator last = printParameters_.variableList_.end();


  while ( iterParam != last && iterParam->tag() != "I" &&
          iterParam->tag() != "V" &&
          iterParam->tag() != "VSTART" &&
          !(iterParam->hasExpressionTag()))
  {
    ++iterParam;
  }

  printCount_ = 0;
  ParameterList::const_iterator iterParam2;
  for (iterParam2 = iterParam; iterParam2 != last; ++iterParam2)
  {
    if (iterParam2->tag() == "I" || iterParam2->tag() == "V" ||
        iterParam2->hasExpressionTag() )
      ++printCount_;
  }

  (* outStreamPtr_) << "#H" << std::endl;

  (*outStreamPtr_) << "SOURCE='Xyce' VERSION='"
                   << N_UTL_Version::getShortVersionString() << "'" << std::endl;

  (*outStreamPtr_) << "TITLE='* " << outputManager_.getNetListFilename() << "'" << std::endl;
  (*outStreamPtr_) << "SUBTITLE='spice probe data'" << std::endl;

  // set the time/date stamp
  (*outStreamPtr_) << getTimeDateStamp();
  outStreamPtr_->setf(ios::scientific);
  outStreamPtr_->precision(printParameters_.streamPrecision_);
  (*outStreamPtr_) << "TEMPERATURE='" << outputManager_.getCircuitTemp();
  (*outStreamPtr_) << "'" << std::endl;

  if (printParameters_.printType_ == PrintType::TRAN)
    (*outStreamPtr_) << "ANALYSIS='Transient Analysis' SERIALNO='12345'" <<  std::endl;
  else
    (*outStreamPtr_) << "ANALYSIS='DC transfer characteristic' " <<
      "SERIALNO='12345'" <<  std::endl;

  (*outStreamPtr_) << "ALLVALUES='NO' COMPLEXVALUES='NO' " <<
    "NODES='" << printCount_ << "'" << std::endl;

  if (printParameters_.printType_ == PrintType::TRAN)
  {
    (*outStreamPtr_) << "SWEEPVAR='Time' SWEEPMODE='VAR_STEP'" <<
      std::endl;
  }
  else
  {
    (*outStreamPtr_) << "SWEEPVAR='Voltage' SWEEPMODE='VAR_STEP'" <<
      std::endl;
  }

  // This line assumes that we're doing a homotopy that goes from 0 to 1.
  // This will never be a transient output.
  (*outStreamPtr_) << "XBEGIN='0.0'  XEND='1.0'"<<std::endl;

  (*outStreamPtr_) << "FORMAT='0 VOLTSorAMPS;EFLOAT : "
                   << "NODEorBRANCH;NODE  '  " << std::endl;

  (*outStreamPtr_) << "DGTLDATA='NO'" << std::endl;

  (*outStreamPtr_) << "#N" << std::endl;

  // print the continuation parameter names:
  int i;
  std::vector<std::string>::const_iterator iter_name;
  for (i=0, iter_name=parameter_names.begin(); iter_name!= parameter_names.end();
       ++iter_name, ++i)
  {
    (*outStreamPtr_) << "'";
    (*outStreamPtr_) << *iter_name;
    (*outStreamPtr_) << "' ";

    if (i > 3)
    {
      i = 0;
      (*outStreamPtr_) << std::endl;
    }
  }

  // print output variable names:
  for (; iterParam != last; ++i)
  {
    (*outStreamPtr_) << "'";
    (*outStreamPtr_) << (*iterParam).tag();
    (*outStreamPtr_) << "' ";

    ++iterParam;

    if (i > 3)
    {
      i = 0;
      (*outStreamPtr_) << std::endl;
    }
  }
  if (i != 0)(*outStreamPtr_) << std::endl;


  (*outStreamPtr_) << flush;

  outStreamPtr_->setf(ios::scientific);
  outStreamPtr_->precision(printParameters_.streamPrecision_);
  outStreamPtr_->setf(ios::left, ios::adjustfield);


  ++stepCount_;
}

void HomotopyProbe::doOutputTime(const N_LAS_Vector *solution_vector, const N_LAS_Vector *state_vector, const N_LAS_Vector *store_vector)
{}

void HomotopyProbe::doOutputMPDE(double time, const N_LAS_Vector *solution_vector)
{}

void HomotopyProbe::doOutputHomotopy(const std::vector<std::string> & parameter_names, const std::vector<double> & parameter_values, const N_LAS_Vector * solution_vector)
{
  outputManager_.setCurrentOutputter(this);
  std::ostringstream ost;
  ost << outputManager_.getProcID();
  string pN(ost.str());

  double tmpTime = outputManager_.getAnaIntPtr()->getTime();

  if (outputManager_.getProcID() == 0)
  {

    if (firstTimeHomotopy_) //Setup Output Stream and Print Out Header
    {
      doOpen();
      doOutputHeader(parameter_names, parameter_values, solution_vector);

      firstTimeHomotopy_ = false;
    }

    outStreamPtr_->width( 0);
    if (printParameters_.printType_ == PrintType::TRAN)
    {
      (*outStreamPtr_) << "#C " << tmpTime << " ";
      (*outStreamPtr_) << printCount_ << std::endl;
    }
    else
    {
      (*outStreamPtr_) << "#C " << outputManager_.getPRINTDCvalue() << " ";
      (*outStreamPtr_) << printCount_ << std::endl;
    }

    //-------------------------------------
    //HOMOTOPY PARAM VALUE OUTPUT GOES HERE
    //-------------------------------------

    for (int iparam=0;iparam < parameter_values.size(); ++iparam)
    {

      if (printParameters_.delimiter_ == "")
        outStreamPtr_->width(printParameters_.streamWidth_);
      else
      {
        outStreamPtr_->width(0);
        if (printParameters_.delimiter_ != "")
          (*outStreamPtr_) << printParameters_.delimiter_;
      }

      (*outStreamPtr_) << parameter_values[iparam];
    }
  } // procID

  ParameterList::const_iterator iterParam = printParameters_.variableList_.begin();
  ParameterList::const_iterator last = printParameters_.variableList_.end();

  bool foundFirstParam=false;
  while ( iterParam != last &&  !foundFirstParam)
  {
    if (iterParam->getSimContext() != UNDEFINED)
    {
      // found the beginning so exit this loop
      foundFirstParam=true;
    }
    else
    {
      ++iterParam;
    }
  }

  int i;
  for (i = 1; iterParam != last; ++iterParam, ++i)
  {
    double result;
    result = 0.0;
    result = outputManager_.getPrintValue(iterParam, solution_vector);

    if (outputManager_.getProcID() == 0)
    {
      if (printParameters_.delimiter_ == "")
        outStreamPtr_->width(printParameters_.streamWidth_);
      else
      {
        outStreamPtr_->width(0);
        if (printParameters_.delimiter_ != "")
          (*outStreamPtr_) << printParameters_.delimiter_;
      }
      (*outStreamPtr_) << result;
    }

    outputManager_.getCommPtr()->barrier();
  }

  if (outputManager_.getProcID() == 0)
    (*outStreamPtr_) << std::endl;
}

void HomotopyProbe::doOutputMORTF(bool origSystem, const double & freq, const Teuchos::SerialDenseMatrix<int, std::complex<double> >& H)
{}

void HomotopyProbe::doOutputHB(
  const std::vector<double>&    timePoints,
  const std::vector<double>&    freqPoints,
  const N_LAS_BlockVector &     timeDomainSolnVec,
  const N_LAS_BlockVector &     freqDomainSolnVecReal,
  const N_LAS_BlockVector &     freqDomainSolnVecImaginary)
{
}

void
HomotopyProbe::doOutputFrequency(double frequency, const N_LAS_Vector *real_solution_vector, const N_LAS_Vector *imaginary_solution_vector)
{}


void HomotopyProbe::doResetOutput()
{}

void HomotopyProbe::doFinishOutput()
{
  firstTimeHomotopy_ = true;
}

void
HomotopyProbe::doFinishOutputStep()
{
  // close the homotopy file.
  if (outStreamPtr_)
  {
    (*outStreamPtr_) << "#;" << std::endl;
  }

  outputManager_.closeFile(outStreamPtr_);
  outStreamPtr_ = 0;
}

TimeRaw::TimeRaw(OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_(""),
    outStreamPtr_( NULL),
    numPoints_( 0),
    numPointsPos_( 0),
    outputRAWTitleAndDate_(false)
{
  if (printParameters_.extension_.empty())
    printParameters_.extension_ = ".raw";

  fixupColumns(outputManager_, printParameters_);
}

TimeRaw::~TimeRaw() {
  outputManager_.closeFile(outStreamPtr_);
}

void TimeRaw::doParse() {
  // prepare output manager to write
  outFilename_ = outputFilename(printParameters_, outputManager_.getNetListFilename());
}


void TimeRaw::doOpen()
{
  if (outputManager_.getProcID() == 0 && outStreamPtr_ == 0)
    outStreamPtr_ = outputManager_.openBinaryFile(outFilename_);
}

void TimeRaw::timeHeader()
{
  std::vector< pair< std::string, char > > nT;
  NodeNamePairMap::iterator name_i, name_end;
  name_end = outputManager_.getAllNodes().end();

  std::ostream &os = *outStreamPtr_;

  if (outputManager_.getProcID() == 0)
  {
    if (!outputRAWTitleAndDate_)
    {
      // in multi-plot RAW files, the Title and Date are only output once
      // afer that each plot gets the same header.  Thus, only do this section
      // once.

      outputRAWTitleAndDate_=true;

      os << "Title: " << outputManager_.getTitle() << std::endl;

      // create formatted timestamp
      const time_t now = time( NULL);
      char timeDate[ 40 ];
      strftime( timeDate, 40, "%a %b %d %I:%M:%S %Y", localtime( &now));
      os << "Date: " << timeDate << std::endl;
    }

    // format plot name
    std::ostringstream plotName;

    if (!outputManager_.getStepParamVec().empty())
    {
      plotName << "Step Analysis: Step " << outputManager_.getStepLoopNumber() + 1
               << " of " << outputManager_.getMaxParamSteps()
               << " params: ";
      std::vector<N_ANP_SweepParam>::const_iterator currItr = outputManager_.getStepParamVec().begin();
      std::vector<N_ANP_SweepParam>::const_iterator endItr = outputManager_.getStepParamVec().end();
      while ( currItr != endItr)
      {
        plotName << " name = " << currItr->name
                 << " value = " << currItr->currentVal << "  ";
        currItr++;
      }
    }

    // while there is a doOutputHeaderAC(), AC can call this function
    // if it outputting a DC operating point.  Thus it's included
    // in this if statement.
    if (printParameters_.printType_ == PrintType::TRAN)
    {
      plotName << "Transient Analysis";
    }
    else if ( printParameters_.printType_ == PrintType::AC)
    {
      plotName << "DC operating point";
    }
    else
    {
      plotName << "DC transfer characteristic";
    }

    os << "Plotname: " << plotName.str() << std::endl;

    // format the flags
    std::string flags("real");
    os << "Flags: " << flags << std::endl;

  } // end proc0 check

  // prepare header for partial dump
  if (outputManager_.getProcID() == 0)
  {
    int numVars = 0;
    switch (printParameters_.printType_) {
      case PrintType::TRAN:
      case PrintType::AC:
        break;

      case PrintType::DC:
        ++numVars; // sweep
        break;
    }

    numVars += printParameters_.variableList_.size();

    // format number of internal and external variables included here + time
    os << "No. Variables: " << numVars << std::endl;
  }

  outputManager_.getCommPtr()->barrier();

  if (outputManager_.getProcID() == 0)
  {
    // format total number of data points & remember the streampos
    os << "No. Points: ";
    numPointsPos_ = outStreamPtr_->tellp(); // <= 344 due to 80 char line limit
    os << "                  " << std::endl; // this will be overwritten

    // for full compatability with spice3 raw file, output version number here
    os << "Version: " << N_UTL_Version::getFullVersionString() << std::endl;

    // NOTE: dimensions, command, options, and scale not considered

    // The variable list appears next.  The format is:
    //         [tab](index) [tab](name) [tab](type)
    // Note that the(index) corresponds to the rawfile, not the soln vec.

    // write variable names for select data points
    std::string tmpNodeName, tmpType;

    int i = 0;
    // write the variable information
    os << "Variables:" << std::endl;

    switch (printParameters_.printType_) {
      case PrintType::TRAN:
        break;

      case PrintType::AC:
        break;

      case PrintType::DC:
        os << "\t" << 0 << "\t" << "sweep\tvoltage\n";
        ++i;
        break;
    }

    for (ParameterList::const_iterator iterParam = printParameters_.variableList_.begin() ; iterParam != printParameters_.variableList_.end(); ++iterParam, ++i)
    {
      // set the type
      if (iterParam->hasExpressionTag()) { tmpType = "expression"; }
      else if (iterParam->tag() == "INDEX")  { }
      else if (iterParam->tag() == "TIME")  { tmpType = "time"; }
      else if (iterParam->tag() == "FREQUENCY")  { tmpType = "frequency"; }
      else if (iterParam->tag()[0] == 'I')  { tmpType = "current";    }
      else if (iterParam->tag()[0] == 'V')  { tmpType = "voltage";    }
      else                              { tmpType = "unknown";    }

      // write the header line
      os << "\t" << i
         << "\t" << (*iterParam).tag()
         << "\t" << tmpType
         << "\n";

      // NOTE: other types, and params & dims, not considered
    }

    // write data marker
    // this string is actually ignored, but the pair of EOL markers is used
    os << "Binary:" << std::endl;
  }
}

void TimeRaw::doOutputTime(const N_LAS_Vector * solnVecPtr, const N_LAS_Vector * stateVecPtr, const N_LAS_Vector * storeVecPtr)
{
  outputManager_.setCurrentOutputter(this);

  // write rawfile header on first call
  if (numPoints_ == 0)
  {
    doOpen();

    timeHeader();
  }

  std::ostream &os = *outStreamPtr_;

  // file IO on proc 0 only
  if (outputManager_.getProcID() == 0)
  {
    // output the time/step values
    switch (printParameters_.printType_) {
      case PrintType::TRAN:
        break;

      case PrintType::AC:
        break;

      case PrintType::DC:
      {
        double dc_sweep = outputManager_.getPRINTDCvalue();
        outStreamPtr_->write((char *) &dc_sweep, sizeof(double));
      }
      break;
    }
  }

  outputManager_.getCommPtr()->barrier();

  // select values to write from .PRINT line if FORMAT=RAW
  for (ParameterList::const_iterator iterParam = printParameters_.variableList_.begin() ; iterParam != printParameters_.variableList_.end(); ++iterParam)
  {
    // retrieve values from all procs
    double tmp = outputManager_.getPrintValue( iterParam, solnVecPtr);
    if ((*iterParam).getSimContext() == TIME_VAR)
      tmp *= printParameters_.outputTimeScaleFactor_;

    // file IO only on proc 0
    if (outputManager_.getProcID() == 0)
    {
      // write binary data to rawfile
      outStreamPtr_->write((char *) &tmp , sizeof(tmp));
    } // end proc0
  } // end for

  outputManager_.getCommPtr()->barrier();

  // keep track of number of datapoints
  ++numPoints_;
}

void TimeRaw::doFinishOutput()
{
  if (outputManager_.getProcID() == 0)
  {
    if (outStreamPtr_ && numPoints_ != 0) {

      // need to move file pointer back to header and
      // write out the number of points.
      long currentFilePos = outStreamPtr_->tellp();

      // locate the position for number of points
      outStreamPtr_->seekp( numPointsPos_);

      // overwrite blank space with value
      (*outStreamPtr_) << numPoints_;

      // move file pointer to the end again.
      outStreamPtr_->seekp( currentFilePos);
    }
  }

  // reset numPoints_ as it is used as a flag to print the header.
  numPoints_ = 0;
}

void
TimeRaw::doFinishOutputStep()
{}

FrequencyRaw::FrequencyRaw(OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_(""),
    outStreamPtr_( NULL),
    numPoints_( 0),
    numPointsPos_( 0),
    outputRAWTitleAndDate_(false)
{
  if (printParameters_.extension_.empty())
    printParameters_.extension_ = ".raw";

  fixupColumns(outputManager_, printParameters_);
}

FrequencyRaw::~FrequencyRaw() {
  outputManager_.closeFile(outStreamPtr_);
}

void FrequencyRaw::doParse() {
  // prepare output manager to write
  numPoints_ = 0;

  outFilename_ = outputFilename(printParameters_, outputManager_.getNetListFilename());
}


void FrequencyRaw::doOpen()
{
  if (outputManager_.getProcID() == 0 && outStreamPtr_ == 0)
    outStreamPtr_ = outputManager_.openBinaryFile(outFilename_);
}


void FrequencyRaw::frequencyHeader()
{
  // This routine was split off of doOutputHeader because an
  // AC analysis can output both real data "doOutput" and complex
  // data via doOutputAC.  Each of these doOutput routines must
  // call a doHeader routine which can specify the right underlying
  // data type.  So rather than add more conditional logic we chose
  // to separate out the functionality for the different use cases.
  int numVars;
  std::vector< pair< std::string, char > > nT;
  NodeNamePairMap::iterator name_i, name_end;
  name_end = outputManager_.getAllNodes().end();

  if (outputManager_.getProcID() == 0)
  {
    if (!outputRAWTitleAndDate_)
    {
      // in multi-plot RAW files, the Title and Date are only output once
      // afer that each plot gets the same header.  Thus, only do this section
      // once.

      outputRAWTitleAndDate_=true;

      (*outStreamPtr_) << "Title: " << outputManager_.getTitle() << std::endl;

      // create formatted timestamp
      const time_t now = time( NULL);
      char timeDate[ 40 ];
      strftime( timeDate, 40, "%a %b %d %I:%M:%S %Y", localtime( &now));
      (*outStreamPtr_) << "Date: " << timeDate << std::endl;
    }

    // format plot name
    std::ostringstream plotName;

    if (!outputManager_.getStepParamVec().empty())
    {
      plotName << "Step Analysis: Step " << outputManager_.getStepLoopNumber() + 1
               << " of " << outputManager_.getMaxParamSteps()
               << " params: ";
      std::vector<N_ANP_SweepParam>::const_iterator currItr = outputManager_.getStepParamVec().begin();
      std::vector<N_ANP_SweepParam>::const_iterator endItr = outputManager_.getStepParamVec().end();
      while ( currItr != endItr)
      {
        plotName << " name = " << currItr->name
                 << " value = " << currItr->currentVal << "  ";
        currItr++;
      }
    }
    if (!outputManager_.getDCParamVec().empty())
    {
      plotName << "DC Sweep: Step " << outputManager_.getDCLoopNumber() + 1
               << " of " << outputManager_.getMaxDCSteps()
               << " params: ";
      std::vector<N_ANP_SweepParam>::const_iterator currItr = outputManager_.getDCParamVec().begin();
      std::vector<N_ANP_SweepParam>::const_iterator endItr = outputManager_.getDCParamVec().end();
      while ( currItr != endItr)
      {
        plotName << " name = " << currItr->name
                 << " value = " << currItr->currentVal << "  ";
        currItr++;
      }
    }
    if (printParameters_.printType_ == PrintType::TRAN)
    {
      plotName << "Transient Analysis";
    }
    else if ( printParameters_.printType_ == PrintType::AC)
    {
      plotName << "AC Analysis";
    }
    else
    {
      plotName << "DC transfer characteristic";
    }

    (*outStreamPtr_) << "Plotname: " << plotName.str() << std::endl;

    // format the flags
    std::string flags("complex");
    (*outStreamPtr_) << "Flags: " << flags << std::endl;

  } // end proc0 check

  // prepare header for partial dump
  if (outputManager_.getProcID() == 0)
  {
    numVars = 0;
    for (ParameterList::const_iterator iterParam = printParameters_.variableList_.begin() ; iterParam != printParameters_.variableList_.end(); ++iterParam)
    {
      // count params
      ++numVars;

      // skip over node names
      if (iterParam->tag() == "I" || iterParam->tag() == "V")
      {
        ++iterParam;
      }
    }
  } // end proc0 check

    // format number of internal and external variables included here + time
  (*outStreamPtr_) << "No. Variables: " << numVars << std::endl;

  outputManager_.getCommPtr()->barrier();

  if (outputManager_.getProcID() == 0)
  {
    // format total number of data points & remember the streampos
    (*outStreamPtr_) << "No. Points: ";
    numPointsPos_ = outStreamPtr_->tellp(); // <= 344 due to 80 char line limit
    (*outStreamPtr_) << "                  " << std::endl; // this will be overwritten

    // for full compatability with spice3 raw file, output version number here
    (*outStreamPtr_) << "Version: " << N_UTL_Version::getFullVersionString() << std::endl;

    // NOTE: dimensions, command, options, and scale not considered

    // write the variable information
    (*outStreamPtr_) << "Variables:" << std::endl;

    // The variable list appears next.  The format is:
    //         [tab](index) [tab](name) [tab](type)
    // Note that the(index) corresponds to the rawfile, not the soln vec.

    // write variable names for select data points
    std::string tmpNodeName, tmpType;
    int i = 0;
    // add timestep header info(Spice3f5 style)
    if (printParameters_.printType_ == PrintType::TRAN)
    {
    }
    else if (printParameters_.printType_ == PrintType::AC)
    {
    }
    else
    {
      (*outStreamPtr_) << "\t" << 0 << "\t";
      (*outStreamPtr_) << "sweep\tvoltage\n";
      ++i;
    }

    for (ParameterList::const_iterator iterParam = printParameters_.variableList_.begin() ; iterParam != printParameters_.variableList_.end(); ++iterParam, ++i)
    {
      // set the type
      if (iterParam->hasExpressionTag()) { tmpType = "expression"; }
      else if (iterParam->tag() == "INDEX")  { }
      else if (iterParam->tag() == "TIME")  { tmpType = "time"; }
      else if (iterParam->tag() == "FREQUENCY")  { tmpType = "frequency"; }
      else if (iterParam->tag()[0] == 'I')  { tmpType = "current";    }
      else if (iterParam->tag()[0] == 'V')  { tmpType = "voltage";    }
      else                              { tmpType = "unknown";    }

      (*outStreamPtr_) << "\t" << i
                       << "\t" << (*iterParam).tag()
                       << "\t" << tmpType
                       << "\n";

      // NOTE: other types, and params & dims, not considered
    }

    // write data marker
    // this string is actually ignored, but the pair of EOL markers is used
    (*outStreamPtr_) << "Binary:" << std::endl;
  } // end proc0 check
}



void FrequencyRaw::doFinishOutput()
{
  if (outputManager_.getProcID() == 0)
  {
    if (outStreamPtr_ && numPoints_ != 0) {
      // need to move file pointer back to header and
      // write out the number of points.
      long currentFelePost = outStreamPtr_->tellp();

      // locate the position for number of points
      outStreamPtr_->seekp( numPointsPos_);

      // overwrite blank space with value
      (*outStreamPtr_) << numPoints_;

      // move file pointer to the end again.
      outStreamPtr_->seekp( currentFelePost);
    }
  }

  // reset numPoints_ as it is used as a flag to print the header.
  numPoints_ = 0;
}

void
FrequencyRaw::doOutputFrequency(double frequency, const N_LAS_Vector *real_solution_vector, const N_LAS_Vector *imaginary_solution_vector)
{
  outputManager_.setCurrentOutputter(this);

  if (!outStreamPtr_)
  {
    doOpen();

    numPoints_ = 0;

    frequencyHeader();
  }

  outputManager_.getCommPtr()->barrier();

  // select values to write from .PRINT line if FORMAT=RAW
  for (ParameterList::const_iterator iterParam = printParameters_.variableList_.begin() ; iterParam != printParameters_.variableList_.end(); ++iterParam)
  {
    // if context type of the param is SOLUTION_VAR then it is complex and 
    // we need to get both the real and imaginary part for output.
    // Any other context type is just real and the imaginary part will be zero
    double varValue = outputManager_.getPrintValue(iterParam, real_solution_vector, NULL, NULL, imaginary_solution_vector);
    double varValueIm = 0.0;
    if( iterParam->getSimContext() == SOLUTION_VAR )
    {
      // note we fake out getPrintValue to return the imaginary component by passing it in as the real vector.
      // not a great idea, but getPrintValue() will be rewritten to return complex types later.
      varValueIm = outputManager_.getPrintValue(iterParam, imaginary_solution_vector, NULL, NULL, imaginary_solution_vector);
    }
    if (outputManager_.getProcID() == 0)
    {
      outStreamPtr_->write((char *)&varValue , sizeof( double));
      outStreamPtr_->write((char *)&varValueIm , sizeof( double));
    }
  }

  outputManager_.getCommPtr()->barrier();

  // keep track of number of datapoints
  ++numPoints_;
}

void
FrequencyRaw::doFinishOutputStep()
{}

TimeRawAscii::TimeRawAscii(OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_( ""),
    outStreamPtr_( NULL),
    numPoints_(0),
    numPointsPos_( 0),
    outputRAWTitleAndDate_(false)
{
  if (printParameters_.extension_.empty())
    printParameters_.extension_ = ".raw";

  fixupColumns(outputManager_, printParameters_);
}

TimeRawAscii::~TimeRawAscii() {
  outputManager_.closeFile(outStreamPtr_);
}

void TimeRawAscii::doParse() {
  outFilename_ = outputFilename(printParameters_, outputManager_.getNetListFilename());
}

void TimeRawAscii::doOpen()
{
  if (outputManager_.getProcID() == 0 && outStreamPtr_ == 0) {
    outStreamPtr_ = outputManager_.openFile(outFilename_);
    
    // set output value characteristics
    outStreamPtr_->setf( ios::scientific);
    outStreamPtr_->precision(8);
    outStreamPtr_->setf( ios::left, ios::adjustfield);
  }
}

void TimeRawAscii::timeHeader()
{
  int numVars;
  std::vector< pair< std::string, char > > nT;
  NodeNamePairMap::iterator name_i, name_end;
  name_end = outputManager_.getAllNodes().end();

  if (outputManager_.getProcID() == 0)
  {
    if (!outputRAWTitleAndDate_)
    {
      // in multi-plot RAW files, the Title and Date are only output once
      // afer that each plot gets the same header.  Thus, only do this section
      // once.

      outputRAWTitleAndDate_ = true;

      (*outStreamPtr_) << "Title: " << outputManager_.getTitle() << std::endl;

      // create formatted timestamp
      const time_t now = time( NULL);
      char timeDate[ 40 ];
      strftime( timeDate, 40, "%a %b %d %I:%M:%S %Y", localtime( &now));
      (*outStreamPtr_) << "Date: " << timeDate << std::endl;
    }

    // format plot name
    std::ostringstream plotName;

    if (!outputManager_.getStepParamVec().empty())
    {
      plotName << "Step Analysis: Step " << outputManager_.getStepLoopNumber() + 1
               << " of " << outputManager_.getMaxParamSteps()
               << " params: ";
      std::vector<N_ANP_SweepParam>::const_iterator currItr = outputManager_.getStepParamVec().begin();
      std::vector<N_ANP_SweepParam>::const_iterator endItr = outputManager_.getStepParamVec().end();
      while ( currItr != endItr)
      {
        plotName << " name = " << currItr->name
                 << " value = " << currItr->currentVal << "  ";
        currItr++;
      }
    }

    // while there is a doOutputHeaderAC(), AC can call this function
    // if it outputting a DC operating point.  Thus it's included
    // in this if statement.
    if (printParameters_.printType_ == PrintType::TRAN)
    {
      plotName << "Transient Analysis";
    }
    else if ( printParameters_.printType_ == PrintType::AC)
    {
      plotName << "DC operating point";
    }
    else
    {
      plotName << "DC transfer characteristic";
    }

    (*outStreamPtr_) << "Plotname: " << plotName.str() << std::endl;

    // format the flags
    std::string flags("real");
    (*outStreamPtr_) << "Flags: " << flags << std::endl;

  } // end proc0 check

  // prepare header for partial dump
  if (outputManager_.getProcID() == 0)
  {
    numVars = 0;
    if (printParameters_.printType_ == PrintType::TRAN)
    {
    }
    else if (printParameters_.printType_ == PrintType::AC)
    {
    }
    else
    {
      ++numVars;
    }
    for (ParameterList::const_iterator iterParam = printParameters_.variableList_.begin() ; iterParam != printParameters_.variableList_.end(); ++iterParam)
    {
      // count params
      ++numVars;

      // skip over node names
      if (iterParam->tag() == "I" || iterParam->tag() == "V")
      {
        ++iterParam;
      }
    }
  } // end proc0 check

  outputManager_.getCommPtr()->barrier();

  if (outputManager_.getProcID() == 0)
  {
    // format number of internal and external variables included here + time
    (*outStreamPtr_) << "No. Variables: " << numVars << std::endl;

    // format total number of data points & remember the streampos
    (*outStreamPtr_) << "No. Points: ";
    numPointsPos_ = outStreamPtr_->tellp(); // <= 344 due to 80 char line limit
    (*outStreamPtr_) << "                  " << std::endl; // this will be overwritten

    // for full compatability with spice3 raw file, output version number here
    (*outStreamPtr_) << "Version: " << N_UTL_Version::getFullVersionString() << std::endl;

    // NOTE: dimensions, command, options, and scale not considered

    // write the variable information
    (*outStreamPtr_) << "Variables:" << std::endl;

    // The variable list appears next.  The format is:
    //         [tab](index) [tab](name) [tab](type)
    // Note that the(index) corresponds to the rawfile, not the soln vec.

    // write variable names for select data points
    std::string tmpNodeName, tmpType;
    int i = 0;
    if (printParameters_.printType_ == PrintType::TRAN)
    {
    }
    else if (printParameters_.printType_ == PrintType::AC)
    {
    }
    else
    {
      // add timestep header info(Spice3f5 style)
      (*outStreamPtr_) << "\t" << 0 << "\t";

      (*outStreamPtr_) << "sweep\tvoltage\n";
      ++i;
    }

    for (ParameterList::const_iterator iterParam = printParameters_.variableList_.begin() ; iterParam != printParameters_.variableList_.end(); ++iterParam, ++i)
    {
      // set the type
      if (iterParam->hasExpressionTag()) { tmpType = "expression"; }
      else if (iterParam->tag() == "INDEX")  { }
      else if (iterParam->tag() == "TIME")  { tmpType = "time"; }
      else if (iterParam->tag() == "FREQUENCY")  { tmpType = "frequency"; }
      else if (iterParam->tag()[0] == 'I')  { tmpType = "current";    }
      else if (iterParam->tag()[0] == 'V')  { tmpType = "voltage";    }
      else                              { tmpType = "unknown";    }

      // write the header line
      (*outStreamPtr_) << "\t" << i
                       << "\t" << (*iterParam).tag()
                       << "\t" << tmpType
                       << "\n";

      // NOTE: other types, and params & dims, not considered
    }

    // write data marker
    (*outStreamPtr_) << "Values:" << std::endl;
  } // end proc0 check
}

void TimeRawAscii::doOutputTime(const N_LAS_Vector * solnVecPtr, const N_LAS_Vector * stateVecPtr, const N_LAS_Vector * storeVecPtr)
{
  outputManager_.setCurrentOutputter(this);

  // write rawfile header on first call
  if (numPoints_ == 0)
  {
    doOpen();

    timeHeader();
  }

  // file IO on proc 0 only
  if (outputManager_.getProcID() == 0)
  {
    // write the index to ascii rawfile
    (*outStreamPtr_) << numPoints_;

  }

  outputManager_.getCommPtr()->barrier();

  // select values to write from .PRINT line if FORMAT=RAW
  for (ParameterList::const_iterator iterParam = printParameters_.variableList_.begin() ; iterParam != printParameters_.variableList_.end(); ++iterParam)
  {
    // retrieve values from all procs
    double tmp = outputManager_.getPrintValue( iterParam, solnVecPtr);
    if ((*iterParam).getSimContext() == TIME_VAR)
      tmp *= printParameters_.outputTimeScaleFactor_;

    // file IO only on proc 0
    if (outputManager_.getProcID() == 0)
    {
      // write formatted values to rawfile
      (*outStreamPtr_) << "\t" << tmp << "\n";

    } // end proc0

  } // end for


  if (outputManager_.getProcID() == 0)
  {
    // write newline to rawfile
    (*outStreamPtr_) << std::endl;
  }

  outputManager_.getCommPtr()->barrier();

  // keep track of number of datapoints
  ++numPoints_;
}
void TimeRawAscii::doFinishOutput()
{
  if (outputManager_.getProcID() == 0)
  {
    if (outStreamPtr_ && numPoints_ != 0) {
      // need to move file pointer back to header and
      // write out the number of points.
      long currentFelePost = outStreamPtr_->tellp();

      // locate the position for number of points
      outStreamPtr_->seekp( numPointsPos_);

      // overwrite blank space with value
      (*outStreamPtr_) << numPoints_;

      // move file pointer to the end again.
      outStreamPtr_->seekp( currentFelePost);
    }
  }

  // reset numPoints_ as it is used as a flag
  // by outputRAW_() to print the header.
  numPoints_=0;
}


void
TimeRawAscii::doFinishOutputStep()
{}


FrequencyRawAscii::FrequencyRawAscii(OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_( ""),
    outStreamPtr_( NULL),
    numPoints_(0),
    numPointsPos_( 0),
    outputRAWTitleAndDate_(false)
{
  if (printParameters_.extension_.empty())
    printParameters_.extension_ = ".raw";

  fixupColumns(outputManager_, printParameters_);
}

FrequencyRawAscii::~FrequencyRawAscii() {
  outputManager_.closeFile(outStreamPtr_);
}

void FrequencyRawAscii::doParse() {
  outFilename_ = outputFilename(printParameters_, outputManager_.getNetListFilename());
}

void FrequencyRawAscii::doOpen()
{
  if (outputManager_.getProcID() == 0 && outStreamPtr_ == 0) {
    outStreamPtr_ = outputManager_.openFile(outFilename_);
    
    // set output value characteristics
    outStreamPtr_->setf( ios::scientific);
    outStreamPtr_->precision(8);
    outStreamPtr_->setf( ios::left, ios::adjustfield);
  }
}

void FrequencyRawAscii::frequencyHeader()
{
  // This routine was split off of doOutputHeader because an
  // AC analysis can output both real data "doOutput" and complex
  // data via doOutputAC.  Each of these doOutput routines must
  // call a doHeader routine which can specify the right underlying
  // data type.  So rather than add more conditional logic we chose
  // to separate out the functionality for the different use cases.
  int numVars;
  std::vector< pair< std::string, char > > nT;
  NodeNamePairMap::iterator name_i, name_end;
  name_end = outputManager_.getAllNodes().end();

  if (outputManager_.getProcID() == 0)
  {
    if (!outputRAWTitleAndDate_)
    {
      // in multi-plot RAW files, the Title and Date are only output once
      // afer that each plot gets the same header.  Thus, only do this section
      // once.

      outputRAWTitleAndDate_ = true;

      (*outStreamPtr_) << "Title: " << outputManager_.getTitle() << std::endl;

      // create formatted timestamp
      const time_t now = time( NULL);
      char timeDate[ 40 ];
      strftime( timeDate, 40, "%a %b %d %I:%M:%S %Y", localtime( &now));
      (*outStreamPtr_) << "Date: " << timeDate << std::endl;
    }

    // format plot name
    std::ostringstream plotName;

    if (!outputManager_.getStepParamVec().empty())
    {
      plotName << "Step Analysis: Step " << outputManager_.getStepLoopNumber() + 1
               << " of " << outputManager_.getMaxParamSteps()
               << " params: ";
      std::vector<N_ANP_SweepParam>::const_iterator currItr = outputManager_.getStepParamVec().begin();
      std::vector<N_ANP_SweepParam>::const_iterator endItr = outputManager_.getStepParamVec().end();
      while ( currItr != endItr)
      {
        plotName << " name = " << currItr->name
                 << " value = " << currItr->currentVal << "  ";
        currItr++;
      }
    }
    if (!outputManager_.getDCParamVec().empty())
    {
      plotName << "DC Sweep: Step " << outputManager_.getDCLoopNumber() + 1
               << " of " << outputManager_.getMaxDCSteps()
               << " params: ";
      std::vector<N_ANP_SweepParam>::const_iterator currItr = outputManager_.getDCParamVec().begin();
      std::vector<N_ANP_SweepParam>::const_iterator endItr = outputManager_.getDCParamVec().end();
      while ( currItr != endItr)
      {
        plotName << " name = " << currItr->name
                 << " value = " << currItr->currentVal << "  ";
        currItr++;
      }
    }
    if (printParameters_.printType_ == PrintType::TRAN)
    {
      plotName << "Transient Analysis";
    }
    else if ( printParameters_.printType_ == PrintType::AC)
    {
      plotName << "AC Analysis";
    }
    else
    {
      plotName << "DC transfer characteristic";
    }

    (*outStreamPtr_) << "Plotname: " << plotName.str() << std::endl;

    // format the flags
    std::string flags("complex");
    (*outStreamPtr_) << "Flags: " << flags << std::endl;

  } // end proc0 check

  // prepare header for partial dump
  if (outputManager_.getProcID() == 0)
  {
    numVars = 0;
    if (printParameters_.printType_ == PrintType::TRAN)
    {
    }
    else if (printParameters_.printType_ == PrintType::AC)
    {
    }
    else
    {
      ++numVars;
    }
    for (ParameterList::const_iterator iterParam = printParameters_.variableList_.begin() ; iterParam != printParameters_.variableList_.end(); ++iterParam)
    {
      // count params
      ++numVars;

      // skip over node names
      if (iterParam->tag() == "I" || iterParam->tag() == "V")
      {
        ++iterParam;
      }
    }

    // format number of internal and external variables included here + time
    (*outStreamPtr_) << "No. Variables: " << numVars << std::endl;

  } // end proc0 check

  outputManager_.getCommPtr()->barrier();

  if (outputManager_.getProcID() == 0)
  {
    // format total number of data points & remember the streampos
    (*outStreamPtr_) << "No. Points: ";
    numPointsPos_ = outStreamPtr_->tellp(); // <= 344 due to 80 char line limit
    (*outStreamPtr_) << "                  " << std::endl; // this will be overwritten

    // for full compatability with spice3 raw file, output version number here
    (*outStreamPtr_) << "Version: " << N_UTL_Version::getFullVersionString() << std::endl;

    // NOTE: dimensions, command, options, and scale not considered

    // write the variable information
    (*outStreamPtr_) << "Variables:" << std::endl;

    // The variable list appears next.  The format is:
    //         [tab](index) [tab](name) [tab](type)
    // Note that the(index) corresponds to the rawfile, not the soln vec.

    // write variable names for select data points
    std::string tmpNodeName, tmpType;
    int i = 0;

    // add timestep header info(Spice3f5 style)
    if (printParameters_.printType_ == PrintType::TRAN)
    {
    }
    else if (printParameters_.printType_ == PrintType::AC)
    {
    }
    else
    {
      (*outStreamPtr_) << "\t" << 0 << "\t";
      (*outStreamPtr_) << "sweep\tvoltage\n";
      ++i;
    }

    for (ParameterList::const_iterator iterParam = printParameters_.variableList_.begin() ; iterParam != printParameters_.variableList_.end(); ++iterParam, ++i)
    {
      // set the type
      if (iterParam->hasExpressionTag()) { tmpType = "expression"; }
      else if (iterParam->tag() == "INDEX")  { }
      else if (iterParam->tag() == "TIME")  { tmpType = "time"; }
      else if (iterParam->tag() == "FREQUENCY")  { tmpType = "frequency"; }
      else if (iterParam->tag()[0] == 'I')  { tmpType = "current";    }
      else if (iterParam->tag()[0] == 'V')  { tmpType = "voltage";    }
      else                              { tmpType = "unknown";    }

      // write the header line
      (*outStreamPtr_) << "\t" << i
                       << "\t" << (*iterParam).tag()
                       << "\t" << tmpType
                       << "\n";

      // NOTE: other types, and params & dims, not considered
    }

    // write data marker
    (*outStreamPtr_) << "Values:" << std::endl;
  } // end proc0 check
}

void
FrequencyRawAscii::doOutputFrequency(double frequency, const N_LAS_Vector *real_solution_vector, const N_LAS_Vector *imaginary_solution_vector)
{
  outputManager_.setCurrentOutputter(this);

  if (!outStreamPtr_)
  {
    doOpen();

    numPoints_ = 0;

    frequencyHeader();
  }

  // file IO on proc 0 only
  if (outputManager_.getProcID() == 0)
  {

    // write the index to ascii rawfile
    (*outStreamPtr_) << numPoints_;

  }

  outputManager_.getCommPtr()->barrier();

  // select values to write from .PRINT line if FORMAT=RAW
  for (ParameterList::const_iterator iterParam = printParameters_.variableList_.begin() ; iterParam != printParameters_.variableList_.end(); ++iterParam)
  {
    // if context type of the param is SOLUTION_VAR then it is complex and 
    // we need to get both the real and imaginary part for output.
    // Any other context type is just real and the imaginary part will be zero
    double varValue = outputManager_.getPrintValue(iterParam, real_solution_vector, NULL, NULL, imaginary_solution_vector);
    double varValueIm = 0.0;
    if( iterParam->getSimContext() == SOLUTION_VAR )
    {
      // note we fake out getPrintValue to return the imaginary component by passing it in as the real vector.
      // not a great idea, but getPrintValue() will be rewritten to return complex types later.
      varValueIm = outputManager_.getPrintValue(iterParam, imaginary_solution_vector, NULL, NULL, imaginary_solution_vector);
    }
    // file IO only on proc 0
    if (outputManager_.getProcID() == 0)
    {
      (*outStreamPtr_) << "\t"  << varValue << ", " << varValueIm << "\n";
    }
  }

  if (outputManager_.getProcID() == 0)
  {
    // write newline to rawfile
    (*outStreamPtr_) << std::endl;
  }

  // sync after writing data for step
  outputManager_.getCommPtr()->barrier();

  // keep track of number of datapoints
  ++numPoints_;
}

void FrequencyRawAscii::doFinishOutput()
{
  if (outputManager_.getProcID() == 0)
  {
    if (outStreamPtr_) {

      // need to move file pointer back to header and
      // write out the number of points.
      long currentFelePost = outStreamPtr_->tellp();

      // locate the position for number of points
      outStreamPtr_->seekp( numPointsPos_);

      // overwrite blank space with value
      (*outStreamPtr_) << numPoints_;

      // move file pointer to the end again.
      outStreamPtr_->seekp( currentFelePost);
    }
  }

  // reset numPoints_ as it is used as a flag
  // by outputRAW_() to print the header.
  numPoints_=0;
}


void
FrequencyRawAscii::doFinishOutputStep()
{}


OverrideRaw::OverrideRaw(OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_( ""),
    outStreamPtr_( NULL),
    numPoints_( 0),
    numPointsPos_( 0),
    outputRAWTitleAndDate_(false)
{
  if (printParameters_.extension_.empty())
    printParameters_.extension_ = ".raw";
}

OverrideRaw::~OverrideRaw() {
  outputManager_.closeFile(outStreamPtr_);
}

void OverrideRaw::doParse() {
  outFilename_ = outputFilename(printParameters_, outputManager_.getNetListFilename());
}

void OverrideRaw::doOpen()
{
  if (outputManager_.getProcID() == 0 && outStreamPtr_ == 0){
    outStreamPtr_ = outputManager_.openBinaryFile(outFilename_);

    // set output value characteristics
    outStreamPtr_->setf( ios::scientific);
    outStreamPtr_->precision(8);
    outStreamPtr_->setf( ios::left, ios::adjustfield);
  }
}


void OverrideRaw::timeHeader()
{
  int numVars;
  std::vector< pair< std::string, char > > nT;
  NodeNamePairMap::iterator name_i, name_end;
  name_end = outputManager_.getAllNodes().end();

  if (outputManager_.getProcID() == 0)
  {
    if (!outputRAWTitleAndDate_)
    {
      // in multi-plot RAW files, the Title and Date are only output once
      // afer that each plot gets the same header.  Thus, only do this section
      // once.

      outputRAWTitleAndDate_=true;

      (*outStreamPtr_) << "Title: " << outputManager_.getTitle() << std::endl;

      // create formatted timestamp
      const time_t now = time( NULL);
      char timeDate[ 40 ];
      strftime( timeDate, 40, "%a %b %d %I:%M:%S %Y", localtime( &now));
      (*outStreamPtr_) << "Date: " << timeDate << std::endl;
    }

    // format plot name
    std::ostringstream plotName;

    if (!outputManager_.getStepParamVec().empty())
    {
      plotName << "Step Analysis: Step " << outputManager_.getStepLoopNumber() + 1
               << " of " << outputManager_.getMaxParamSteps()
               << " params: ";
      std::vector<N_ANP_SweepParam>::const_iterator currItr = outputManager_.getStepParamVec().begin();
      std::vector<N_ANP_SweepParam>::const_iterator endItr = outputManager_.getStepParamVec().end();
      while ( currItr != endItr)
      {
        plotName << " name = " << currItr->name
                 << " value = " << currItr->currentVal << "  ";
        currItr++;
      }
    }

    // while there is a doOutputHeaderAC(), AC can call this function
    // if it outputting a DC operating point.  Thus it's included
    // in this if statement.
    if (printParameters_.printType_ == PrintType::TRAN)
    {
      plotName << "Transient Analysis";
    }
    else if ( printParameters_.printType_ == PrintType::AC)
    {
      plotName << "DC operating point";
    }
    else
    {
      plotName << "DC transfer characteristic";
    }

    (*outStreamPtr_) << "Plotname: " << plotName.str() << std::endl;

    // format the flags
    std::string flags("real");
    (*outStreamPtr_) << "Flags: " << flags << std::endl;

  } // end proc0 check

  // prepare header for partial dump
  // get node types
  std::vector< char > typeRefs;
  outputManager_.getTopPtr()->returnVarTypeVec( typeRefs);

  // count data points
  numVars = outputManager_.getAllNodes().size();

#ifdef Xyce_PARALLEL_MPI

  int pos, s, e, len, final;
  char t;

  // store the total num of vars
  outputManager_.getCommPtr()->barrier();
  outputManager_.getCommPtr()->sumAll( &numVars, &final, 1);
  numVars = final;

  // limiting names to not longer than 64 chars +1(size) + 1(type)
  int bSize = sizeof( int) +( final *( 65 + sizeof( int)));
  char * b = new char [ bSize ];

  if (outputManager_.getProcID() == 0)
  {

#endif

    // store local names
    for ( name_i = outputManager_.getAllNodes().begin(); name_i != name_end ; ++name_i)
    {
      nT.push_back( pair< std::string, char>((*name_i).first,
                                             typeRefs[(*name_i).second.first]));
    }

#ifdef Xyce_PARALLEL_MPI

    // receive data points from procN in order
    for ( int p = 1; p < outputManager_.getNumProcs(); ++p)
    {
      pos = 0;

      outputManager_.getCommPtr()->recv( &s, 1, p);
      outputManager_.getCommPtr()->recv( b, s, p);
      outputManager_.getCommPtr()->unpack( b, bSize, pos, &e, 1);

      for ( int j = 0; j < e; ++j)
      {
        outputManager_.getCommPtr()->unpack( b, bSize, pos, &len, 1);

        std::string name( &(b[pos]), len);
        pos += len * sizeof( char);

        outputManager_.getCommPtr()->unpack( b, bSize, pos, &t, 1);

        // store remote names
        nT.push_back( pair< std::string, char >( name, t));
      }
    }
  } // end proc0 check

  else
  {
    pos = 0;

    e = outputManager_.getAllNodes().size();
    s = sizeof( int);
    outputManager_.getCommPtr()->pack( &e, 1, b, bSize, pos);

    for ( name_i = outputManager_.getAllNodes().begin(); name_i != name_end ; ++name_i)
    {
      len =((*name_i).first).size();
      s += sizeof( int);
      outputManager_.getCommPtr()->pack( &len, 1, b, bSize, pos);

      s += len * sizeof( char);
      outputManager_.getCommPtr()->pack(((*name_i).first).c_str(), len, b, bSize, pos);

      t = typeRefs[(*name_i).second.first];
      s += sizeof( char);
      outputManager_.getCommPtr()->pack( &t, 1, b, bSize, pos);
    }

    outputManager_.getCommPtr()->send( &s, 1, 0);
    outputManager_.getCommPtr()->send( b, s, 0);
  } // end procN check

    // clean up
  if (b != NULL)
  {
    free( b);
  }

#endif

  // format number of internal and external variables included here + time
  if (outputManager_.getProcID() == 0)
    (*outStreamPtr_) << "No. Variables: " << numVars + 1 << std::endl;

  outputManager_.getCommPtr()->barrier();

  if (outputManager_.getProcID() == 0)
  {
    // format total number of data points & remember the streampos
    (*outStreamPtr_) << "No. Points: ";
    numPointsPos_ = outStreamPtr_->tellp(); // <= 344 due to 80 char line limit
    (*outStreamPtr_) << "                  " << std::endl; // this will be overwritten

    // for full compatability with spice3 raw file, output version number here
    (*outStreamPtr_) << "Version: " << N_UTL_Version::getFullVersionString() << std::endl;

    // NOTE: dimensions, command, options, and scale not considered

    // The variable list appears next.  The format is:
    //         [tab](index) [tab](name) [tab](type)
    // Note that the(index) corresponds to the rawfile, not the soln vec.

    // write variable names for select data points
    std::string tmpNodeName, tmpType;
    // write the variable information
    (*outStreamPtr_) << "Variables:" << std::endl;

    // add timestep header info(Spice3f5 style)
    (*outStreamPtr_) << "\t" << 0 << "\t";

    if (printParameters_.printType_ == PrintType::TRAN)
    {
      (*outStreamPtr_) << "time\ttime\n";
    }
    else if (printParameters_.printType_ == PrintType::AC)
    {
      (*outStreamPtr_) << "frequency\tfrequency\n";
    }
    else
    {
      (*outStreamPtr_) << "sweep\tvoltage\n";
    }

    std::string::size_type uPos;
    for ( int i = 0; i < numVars; ++i)
    {
      tmpNodeName =(nT[i]).first;

      // format is:  [tab](index) [tab](name) [tab](type)
      //   the index corresponds to the rawfile, not the soln vec
      if (strspn( tmpNodeName.c_str(), "0123456789") == tmpNodeName.size())
      {
        // sconvert, spice3f5, & chilespice wrap numeric voltage node names in V()
        tmpNodeName = "V(" +(nT[i]).first + ")";
      }

      uPos = tmpNodeName.rfind( "_", tmpNodeName.size()); if (uPos != std::string::npos)
                                                          {
                                                            tmpNodeName.replace( uPos, 1, "#");
                                                          }

      (*outStreamPtr_) << "\t" << i + 1 << "\t" << tmpNodeName << "\t";

      if ((nT[i]).second == 'I' ||(nT[i]).second == 'i')
      {
        (*outStreamPtr_) << "current\n" ;
      }
      else
      {
        (*outStreamPtr_) << "voltage\n";
      }

      // NOTE: other types, and params & dims, not considered
    }

    // write data marker
    // this string is actually ignored, but the pair of EOL markers is used
    (*outStreamPtr_) << "Binary:" << std::endl;
  } // end proc0 check
}

void OverrideRaw::frequencyHeader()
{
  // This routine was split off of doOutputHeader because an
  // AC analysis can output both real data "doOutput" and complex
  // data via doOutputAC.  Each of these doOutput routines must
  // call a doHeader routine which can specify the right underlying
  // data type.  So rather than add more conditional logic we chose
  // to separate out the functionality for the different use cases.
  int numVars;
  std::vector< pair< std::string, char > > nT;
  NodeNamePairMap::iterator name_i, name_end;
  name_end = outputManager_.getAllNodes().end();

  if (outputManager_.getProcID() == 0)
  {
    if (!outputRAWTitleAndDate_)
    {
      // in multi-plot RAW files, the Title and Date are only output once
      // afer that each plot gets the same header.  Thus, only do this section
      // once.

      outputRAWTitleAndDate_=true;

      (*outStreamPtr_) << "Title: " << outputManager_.getTitle() << std::endl;

      // create formatted timestamp
      const time_t now = time( NULL);
      char timeDate[ 40 ];
      strftime( timeDate, 40, "%a %b %d %I:%M:%S %Y", localtime( &now));
      (*outStreamPtr_) << "Date: " << timeDate << std::endl;
    }

    // format plot name
    std::ostringstream plotName;

    if (!outputManager_.getStepParamVec().empty())
    {
      plotName << "Step Analysis: Step " << outputManager_.getStepLoopNumber() + 1
               << " of " << outputManager_.getMaxParamSteps()
               << " params: ";
      std::vector<N_ANP_SweepParam>::const_iterator currItr = outputManager_.getStepParamVec().begin();
      std::vector<N_ANP_SweepParam>::const_iterator endItr = outputManager_.getStepParamVec().end();
      while ( currItr != endItr)
      {
        plotName << " name = " << currItr->name
                 << " value = " << currItr->currentVal << "  ";
        currItr++;
      }
    }
    if (!outputManager_.getDCParamVec().empty())
    {
      plotName << "DC Sweep: Step " << outputManager_.getDCLoopNumber() + 1
               << " of " << outputManager_.getMaxDCSteps()
               << " params: ";
      std::vector<N_ANP_SweepParam>::const_iterator currItr = outputManager_.getDCParamVec().begin();
      std::vector<N_ANP_SweepParam>::const_iterator endItr = outputManager_.getDCParamVec().end();
      while ( currItr != endItr)
      {
        plotName << " name = " << currItr->name
                 << " value = " << currItr->currentVal << "  ";
        currItr++;
      }
    }
    if (printParameters_.printType_ == PrintType::TRAN)
    {
      plotName << "Transient Analysis";
    }
    else if ( printParameters_.printType_ == PrintType::AC)
    {
      plotName << "AC Analysis";
    }
    else
    {
      plotName << "DC transfer characteristic";
    }

    (*outStreamPtr_) << "Plotname: " << plotName.str() << std::endl;

    // format the flags
    std::string flags("complex");
    (*outStreamPtr_) << "Flags: " << flags << std::endl;

  } // end proc0 check

  // prepare header for partial dump
  // get node types
  std::vector< char > typeRefs;
  outputManager_.getTopPtr()->returnVarTypeVec( typeRefs);

  // count data points
  numVars = outputManager_.getAllNodes().size();

#ifdef Xyce_PARALLEL_MPI

  int pos, s, e, len, final;
  char t;

  // store the total num of vars
  outputManager_.getCommPtr()->barrier();
  outputManager_.getCommPtr()->sumAll( &numVars, &final, 1);
  numVars = final;

  // limiting names to not longer than 64 chars +1(size) + 1(type)
  int bSize = sizeof( int) +( final *( 65 + sizeof( int)));
  char * b = new char [ bSize ];

  if (outputManager_.getProcID() == 0)
  {

#endif

    // store local names
    for ( name_i = outputManager_.getAllNodes().begin(); name_i != name_end ; ++name_i)
    {
      nT.push_back( pair< std::string, char>((*name_i).first,
                                             typeRefs[(*name_i).second.first]));
    }

#ifdef Xyce_PARALLEL_MPI

    // receive data points from procN in order
    for ( int p = 1; p < outputManager_.getNumProcs(); ++p)
    {
      pos = 0;

      outputManager_.getCommPtr()->recv( &s, 1, p);
      outputManager_.getCommPtr()->recv( b, s, p);
      outputManager_.getCommPtr()->unpack( b, bSize, pos, &e, 1);

      for ( int j = 0; j < e; ++j)
      {
        outputManager_.getCommPtr()->unpack( b, bSize, pos, &len, 1);

        std::string name( &(b[pos]), len);
        pos += len * sizeof( char);

        outputManager_.getCommPtr()->unpack( b, bSize, pos, &t, 1);

        // store remote names
        nT.push_back( pair< std::string, char >( name, t));
      }
    }
  } // end proc0 check

  else
  {
    pos = 0;

    e = outputManager_.getAllNodes().size();
    s = sizeof( int);
    outputManager_.getCommPtr()->pack( &e, 1, b, bSize, pos);

    for ( name_i = outputManager_.getAllNodes().begin(); name_i != name_end ; ++name_i)
    {
      len =((*name_i).first).size();
      s += sizeof( int);
      outputManager_.getCommPtr()->pack( &len, 1, b, bSize, pos);

      s += len * sizeof( char);
      outputManager_.getCommPtr()->pack(((*name_i).first).c_str(), len, b, bSize, pos);

      t = typeRefs[(*name_i).second.first];
      s += sizeof( char);
      outputManager_.getCommPtr()->pack( &t, 1, b, bSize, pos);
    }

    outputManager_.getCommPtr()->send( &s, 1, 0);
    outputManager_.getCommPtr()->send( b, s, 0);
  } // end procN check

    // clean up
  if (b != NULL)
  {
    free( b);
  }

#endif

  // format number of internal and external variables included here + time
  if (outputManager_.getProcID() == 0)
    (*outStreamPtr_) << "No. Variables: " << numVars + 1 << std::endl;

  outputManager_.getCommPtr()->barrier();

  if (outputManager_.getProcID() == 0)
  {
    // format total number of data points & remember the streampos
    (*outStreamPtr_) << "No. Points: ";
    numPointsPos_ = outStreamPtr_->tellp();
    (*outStreamPtr_) << "                  " << std::endl; // this will be overwritten

    // for full compatability with spice3 raw file, output version number here
    (*outStreamPtr_) << "Version: " << N_UTL_Version::getFullVersionString() << std::endl;

    // NOTE: dimensions, command, options, and scale not considered

    // write the variable information
    (*outStreamPtr_) << "Variables:" << std::endl;

    // The variable list appears next.  The format is:
    //         [tab](index) [tab](name) [tab](type)
    // Note that the(index) corresponds to the rawfile, not the soln vec.

    // write variable names for select data points
    std::string tmpNodeName, tmpType;
    // add timestep header info(Spice3f5 style)
    (*outStreamPtr_) << "\t" << 0 << "\t";

    if (printParameters_.printType_ == PrintType::TRAN)
    {
      (*outStreamPtr_) << "time\ttime\n";
    }
    else if (printParameters_.printType_ == PrintType::AC)
    {
      (*outStreamPtr_) << "frequency\tfrequency\n";
    }
    else
    {
      (*outStreamPtr_) << "sweep\tvoltage\n";
    }

    std::string::size_type uPos;
    for ( int i = 0; i < numVars; ++i)
    {
      tmpNodeName =(nT[i]).first;

      // format is:  [tab](index) [tab](name) [tab](type)
      //   the index corresponds to the rawfile, not the soln vec
      if (strspn( tmpNodeName.c_str(), "0123456789") == tmpNodeName.size())
      {
        // sconvert, spice3f5, & chilespice wrap numeric voltage node names in V()
        tmpNodeName = "V(" +(nT[i]).first + ")";
      }

      uPos = tmpNodeName.rfind( "_", tmpNodeName.size()); if (uPos != std::string::npos)
                                                          {
                                                            tmpNodeName.replace( uPos, 1, "#");
                                                          }

      (*outStreamPtr_) << "\t" << i + 1 << "\t" << tmpNodeName << "\t";

      if ((nT[i]).second == 'I' ||(nT[i]).second == 'i')
      {
        (*outStreamPtr_) << "current\n" ;
      }
      else
      {
        (*outStreamPtr_) << "voltage\n";
      }

      // NOTE: other types, and params & dims, not considered
    }

    // write data marker
    // this string is actually ignored, but the pair of EOL markers is used
    (*outStreamPtr_) << "Binary:" << std::endl;
  } // end proc0 check
}


void OverrideRaw::doOutputTime(const N_LAS_Vector * solnVecPtr, const N_LAS_Vector * stateVecPtr, const N_LAS_Vector * storeVecPtr)
{
  outputManager_.setCurrentOutputter(this);

  // write rawfile header on first call
  if (numPoints_ == 0)
  {
    doOpen();

    timeHeader();
  }

  // file IO on proc 0 only
  if (outputManager_.getProcID() == 0)
  {
    double tmp = 0.0;

// output the time/step values
    if (printParameters_.printType_ == PrintType::TRAN)
    {
      tmp = outputManager_.getCircuitTime();
      tmp *= printParameters_.outputTimeScaleFactor_;
    }
    else
    {
      tmp = outputManager_.getPRINTDCvalue();
    }
    outStreamPtr_->write((char *)&tmp , sizeof( double));
  }

  outputManager_.getCommPtr()->barrier();

  NodeNamePairMap::iterator name_i, name_end;
  name_end = outputManager_.getAllNodes().end();

#ifdef Xyce_PARALLEL_MPI

  int bSize, pos, s;
  int e = outputManager_.getAllNodes().size();

  outputManager_.getCommPtr()->barrier();
  outputManager_.getCommPtr()->maxAll( &e, &bSize, 1);

  bSize = bSize * sizeof( double);
  char * b = new char [ bSize ];

  if (outputManager_.getProcID() == 0)
  {

#endif

    // dump proc0 data points
    for ( name_i = outputManager_.getAllNodes().begin(); name_i != name_end ; ++name_i)
    {
      double tmp =(*solnVecPtr)[(*name_i).second.first];
      if (fabs( tmp) < outputManager_.getFilter())
      {
        tmp = 0.0;
      }

      // write binary data to rawfile
      outStreamPtr_->write((char *)&tmp , sizeof( double));
    }

#ifdef Xyce_PARALLEL_MPI

    // receive data points from procN in order
    for ( int p = 1; p < outputManager_.getNumProcs(); ++p)
    {
      double tmp;
      outputManager_.getCommPtr()->recv( &e, 1, p);
      s = e * sizeof( double);
      pos = 0;
      outputManager_.getCommPtr()->recv( b, s, p);
      for ( int j = 0; j < e; j++)
      {
        outputManager_.getCommPtr()->unpack( b, bSize, pos, &tmp, 1);

        // write binary data to rawfile
        outStreamPtr_->write((char *)&tmp , sizeof( double));
      }
    }

  } // end proc0 work

  else
  {
    // pack and send local data points to proc0
    s = e * sizeof( double);
    pos = 0;

    for ( name_i = outputManager_.getAllNodes().begin(); name_i != name_end ; ++name_i)
    {
      double tmp =(*solnVecPtr)[(*name_i).second.first];
      if (fabs( tmp) < outputManager_.getFilter())
      {
        tmp = 0.0;
      }

      outputManager_.getCommPtr()->pack( &tmp, 1, b, bSize, pos);
    }

    outputManager_.getCommPtr()->send( &e, 1, 0);
    outputManager_.getCommPtr()->send( b, s, 0);
  }

  // clean up
  delete b;
#endif

  outputManager_.getCommPtr()->barrier();

  // keep track of number of datapoints
  ++numPoints_;
}

void OverrideRaw::doResetOutput()
{}

void OverrideRaw::doFinishOutput()
{
  if (outputManager_.getProcID() == 0)
  {
    // need to move file pointer back to header and
    // write out the number of points.
    long currentFelePost = outStreamPtr_->tellp();

    // locate the position for number of points
    outStreamPtr_->seekp( numPointsPos_);

    // overwrite blank space with value
    (*outStreamPtr_) << numPoints_;

    // move file pointer to the end again.
    outStreamPtr_->seekp( currentFelePost);
  }

  // reset numPoints_ as it is used as a flag
  // by outputRAW_() to print the header.
  numPoints_=0;
}

void
OverrideRaw::doOutputFrequency(double frequency, const N_LAS_Vector *real_solution_vector, const N_LAS_Vector *imaginary_solution_vector)
{
  outputManager_.setCurrentOutputter(this);
  double tmp;

  if (numPoints_ == 0)
  {
    doOpen();

    frequencyHeader();
  }

  outputManager_.getCommPtr()->barrier();

  tmp=frequency;

  // file IO only on proc 0
  if (outputManager_.getProcID() == 0)
  {
    outStreamPtr_->write((char *)&tmp , sizeof( double));
    tmp = 0.0;
    outStreamPtr_->write((char *)&tmp , sizeof( double));
  }

  NodeNamePairMap::iterator name_i, name_end;
  name_end = outputManager_.getAllNodes().end();

#ifdef Xyce_PARALLEL_MPI

  int bSize, pos, s;
  int e = outputManager_.getAllNodes().size();

  outputManager_.getCommPtr()->barrier();
  outputManager_.getCommPtr()->maxAll( &e, &bSize, 1);

  bSize = bSize * sizeof( double);
  char * b = new char [ bSize ];

  if (outputManager_.getProcID() == 0)
  {

#endif

    // dump proc0 data points
    for ( name_i = outputManager_.getAllNodes().begin(); name_i != name_end ; ++name_i)
    {
      tmp =(*real_solution_vector)[(*name_i).second.first];
      double tmpImg =(*imaginary_solution_vector)[(*name_i).second.first];
      if (fabs( tmp) < outputManager_.getFilter())
      {
        tmp = 0.0;
      }
      if (fabs( tmpImg) < outputManager_.getFilter())
      {
        tmpImg = 0.0;
      }

      // write formatted values to rawfile
      if (outputManager_.getProcID() == 0)
      {
        outStreamPtr_->write((char *)&tmp , sizeof( double));
        outStreamPtr_->write((char *)&tmpImg , sizeof( double));
      }
    }

#ifdef Xyce_PARALLEL_MPI

    // receive data points from procN in order
    for ( int p = 1; p < outputManager_.getNumProcs(); ++p)
    {
      outputManager_.getCommPtr()->recv( &e, 1, p);
      s = e * sizeof( double);
      pos = 0;
      outputManager_.getCommPtr()->recv( b, s, p);
      for ( int j = 0; j < e; j++)
      {
        outputManager_.getCommPtr()->unpack( b, bSize, pos, &tmp, 1);

        // write formatted values to rawfile
        if (outputManager_.getProcID() == 0)
        {
          outStreamPtr_->write((char *)&tmp , sizeof( double));
          tmp = 0;
          outStreamPtr_->write((char *)&tmp , sizeof( double));
        }
      }
    }

  } // end proc0 work

  else
  {
    // pack and send local data points to proc0
    s = e * sizeof( double);
    pos = 0;

    for ( name_i = outputManager_.getAllNodes().begin(); name_i != name_end ; ++name_i)
    {
      tmp =(*real_solution_vector)[(*name_i).second.first];
      if (fabs( tmp) < outputManager_.getFilter())
      {
        tmp = 0.0;
      }

      outputManager_.getCommPtr()->pack( &tmp, 1, b, bSize, pos);
    }

    outputManager_.getCommPtr()->send( &e, 1, 0);
    outputManager_.getCommPtr()->send( b, s, 0);
  }

  // clean up
  delete b;
#endif

  outputManager_.getCommPtr()->barrier();

  // keep track of number of datapoints
  ++numPoints_;
}

void OverrideRaw::doOutputHB(
  const std::vector<double>&    timePoints,
  const std::vector<double>&    freqPoints,
  const N_LAS_BlockVector &     timeDomainSolnVec,
  const N_LAS_BlockVector &     freqDomainSolnVecReal,
  const N_LAS_BlockVector &     freqDomainSolnVecImaginary)
{}

void OverrideRaw::doOutputMPDE(double time, const N_LAS_Vector *solution_vector)
{}

void OverrideRaw::doOutputHomotopy(const std::vector<std::string> & parameter_names, const std::vector<double> & parameter_values, const N_LAS_Vector * solution_vector)
{}

void OverrideRaw::doOutputMORTF(bool origSystem, const double & freq, const Teuchos::SerialDenseMatrix<int, std::complex<double> >& H)
{}

void
OverrideRaw::doFinishOutputStep()
{}


OverrideRawAscii::OverrideRawAscii(OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_( ""),
    outStreamPtr_( NULL),
    numPoints_( 0),
    numPointsPos_( 0),
    outputRAWTitleAndDate_(false)
{
  if (printParameters_.extension_.empty())
    printParameters_.extension_ = ".raw";
}

OverrideRawAscii::~OverrideRawAscii() {
  outputManager_.closeFile(outStreamPtr_);
}

void OverrideRawAscii::doParse() {
  outFilename_ = outputFilename(printParameters_, outputManager_.getNetListFilename());
}

void OverrideRawAscii::doOpen()
{
  if (outputManager_.getProcID() == 0 && outStreamPtr_ == 0) {
    outStreamPtr_ = outputManager_.openFile(outFilename_);

    // set output value characteristics
    outStreamPtr_->setf( ios::scientific);
    outStreamPtr_->precision(8);
    outStreamPtr_->setf( ios::left, ios::adjustfield);
  }
}

void OverrideRawAscii::timeHeader()
{
  int numVars;
  std::vector< pair< std::string, char > > nT;
  NodeNamePairMap::iterator name_i, name_end;
  name_end = outputManager_.getAllNodes().end();

  if (outputManager_.getProcID() == 0)
  {
    if (!outputRAWTitleAndDate_)
    {
      // in multi-plot RAW files, the Title and Date are only output once
      // afer that each plot gets the same header.  Thus, only do this section
      // once.

      outputRAWTitleAndDate_ = true;

      (*outStreamPtr_) << "Title: " << outputManager_.getTitle() << std::endl;

      // create formatted timestamp
      const time_t now = time( NULL);
      char timeDate[ 40 ];
      strftime( timeDate, 40, "%a %b %d %I:%M:%S %Y", localtime( &now));
      (*outStreamPtr_) << "Date: " << timeDate << std::endl;
    }

    // format plot name
    std::ostringstream plotName;

    if (!outputManager_.getStepParamVec().empty())
    {
      plotName << "Step Analysis: Step " << outputManager_.getStepLoopNumber() + 1
               << " of " << outputManager_.getMaxParamSteps()
               << " params: ";
      std::vector<N_ANP_SweepParam>::const_iterator currItr = outputManager_.getStepParamVec().begin();
      std::vector<N_ANP_SweepParam>::const_iterator endItr = outputManager_.getStepParamVec().end();
      while ( currItr != endItr)
      {
        plotName << " name = " << currItr->name
                 << " value = " << currItr->currentVal << "  ";
        currItr++;
      }
    }
    if (!outputManager_.getDCParamVec().empty())
    {
      plotName << "DC Sweep: Step " << outputManager_.getDCLoopNumber() + 1
               << " of " << outputManager_.getMaxDCSteps()
               << " params: ";
      std::vector<N_ANP_SweepParam>::const_iterator currItr = outputManager_.getDCParamVec().begin();
      std::vector<N_ANP_SweepParam>::const_iterator endItr = outputManager_.getDCParamVec().end();
      while ( currItr != endItr)
      {
        plotName << " name = " << currItr->name
                 << " value = " << currItr->currentVal << "  ";
        currItr++;
      }
    }
    // while there is a doOutputHeaderAC(), AC can call this function
    // if it outputting a DC operating point.  Thus it's included
    // in this if statement.
    if (printParameters_.printType_ == PrintType::TRAN)
    {
      plotName << "Transient Analysis";
    }
    else if ( printParameters_.printType_ == PrintType::AC)
    {
      plotName << "DC operating point";
    }
    else
    {
      plotName << "DC transfer characteristic";
    }

    (*outStreamPtr_) << "Plotname: " << plotName.str() << std::endl;

    // format the flags
    std::string flags("real");
    (*outStreamPtr_) << "Flags: " << flags << std::endl;

  } // end proc0 check

  // prepare header for partial dump
  // get node types
  std::vector< char > typeRefs;
  outputManager_.getTopPtr()->returnVarTypeVec( typeRefs);

  // count data points
  numVars = outputManager_.getAllNodes().size();

#ifdef Xyce_PARALLEL_MPI

  int pos, s, e, len, final;
  char t;

  // store the total num of vars
  outputManager_.getCommPtr()->barrier();
  outputManager_.getCommPtr()->sumAll( &numVars, &final, 1);
  numVars = final;

  // limiting names to not longer than 64 chars +1(size) + 1(type)
  int bSize = sizeof( int) +( final *( 65 + sizeof( int)));
  char * b = new char [ bSize ];

  if (outputManager_.getProcID() == 0)
  {

#endif
    ++numVars; // SWEEP, TIME or FREQUENCY

    // store local names
    for ( name_i = outputManager_.getAllNodes().begin(); name_i != name_end ; ++name_i)
    {
      nT.push_back( pair< std::string, char>((*name_i).first,
                                             typeRefs[(*name_i).second.first]));
    }

#ifdef Xyce_PARALLEL_MPI

    // receive data points from procN in order
    for ( int p = 1; p < outputManager_.getNumProcs(); ++p)
    {
      pos = 0;

      outputManager_.getCommPtr()->recv( &s, 1, p);
      outputManager_.getCommPtr()->recv( b, s, p);
      outputManager_.getCommPtr()->unpack( b, bSize, pos, &e, 1);

      for ( int j = 0; j < e; ++j)
      {
        outputManager_.getCommPtr()->unpack( b, bSize, pos, &len, 1);

        std::string name( &(b[pos]), len);
        pos += len * sizeof( char);

        outputManager_.getCommPtr()->unpack( b, bSize, pos, &t, 1);

        // store remote names
        nT.push_back( pair< std::string, char >( name, t));
      }
    }
  } // end proc0 check

  else
  {
    pos = 0;

    e = outputManager_.getAllNodes().size();
    s = sizeof( int);
    outputManager_.getCommPtr()->pack( &e, 1, b, bSize, pos);

    for ( name_i = outputManager_.getAllNodes().begin(); name_i != name_end ; ++name_i)
    {
      len =((*name_i).first).size();
      s += sizeof( int);
      outputManager_.getCommPtr()->pack( &len, 1, b, bSize, pos);

      s += len * sizeof( char);
      outputManager_.getCommPtr()->pack(((*name_i).first).c_str(), len, b, bSize, pos);

      t = typeRefs[(*name_i).second.first];
      s += sizeof( char);
      outputManager_.getCommPtr()->pack( &t, 1, b, bSize, pos);
    }

    outputManager_.getCommPtr()->send( &s, 1, 0);
    outputManager_.getCommPtr()->send( b, s, 0);
  } // end procN check

    // clean up
  if (b != NULL)
  {
    free( b);
  }

#endif

  outputManager_.getCommPtr()->barrier();

  if (outputManager_.getProcID() == 0)
  {
    // format number of internal and external variables included here + time
    (*outStreamPtr_) << "No. Variables: " << numVars << std::endl;

    // format total number of data points & remember the streampos
    (*outStreamPtr_) << "No. Points: ";
    numPointsPos_ = outStreamPtr_->tellp(); // <= 344 due to 80 char line limit
    (*outStreamPtr_) << "                  " << std::endl; // this will be overwritten

    // for full compatability with spice3 raw file, output version number here
    (*outStreamPtr_) << "Version: " << N_UTL_Version::getFullVersionString() << std::endl;

    // NOTE: dimensions, command, options, and scale not considered

    // write the variable information
    (*outStreamPtr_) << "Variables:" << std::endl;

    // The variable list appears next.  The format is:
    //         [tab](index) [tab](name) [tab](type)
    // Note that the(index) corresponds to the rawfile, not the soln vec.

    // write variable names for select data points
    std::string tmpNodeName, tmpType;
    // add timestep header info(Spice3f5 style)
    (*outStreamPtr_) << "\t" << 0 << "\t";

    if (printParameters_.printType_ == PrintType::TRAN)
    {
      (*outStreamPtr_) << "time\ttime\n";
    }
    else if (printParameters_.printType_ == PrintType::AC)
    {
      (*outStreamPtr_) << "frequency\tfrequency\n";
    }
    else
    {
      (*outStreamPtr_) << "sweep\tvoltage\n";
    }


    std::string::size_type uPos;
    for ( int i = 0; i < numVars - 1; ++i)
    {
      tmpNodeName =(nT[i]).first;

      // format is:  [tab](index) [tab](name) [tab](type)
      //   the index corresponds to the rawfile, not the soln vec
      if (strspn( tmpNodeName.c_str(), "0123456789") == tmpNodeName.size())
      {
        // sconvert, spice3f5, & chilespice wrap numeric voltage node names in V()
        tmpNodeName = "V(" +(nT[i]).first + ")";
      }

      uPos = tmpNodeName.rfind( "_", tmpNodeName.size()); if (uPos != std::string::npos)
                                                          {
                                                            tmpNodeName.replace( uPos, 1, "#");
                                                          }

      (*outStreamPtr_) << "\t" << i + 1 << "\t" << tmpNodeName << "\t";

      if ((nT[i]).second == 'I' ||(nT[i]).second == 'i')
      {
        (*outStreamPtr_) << "current\n" ;
      }
      else
      {
        (*outStreamPtr_) << "voltage\n";
      }

      // NOTE: other types, and params & dims, not considered
    }

    // write data marker
    if (outputManager_.getProcID() == 0)
      (*outStreamPtr_) << "Values:" << std::endl;
  } // end proc0 check
}


void OverrideRawAscii::frequencyHeader()
{
  // This routine was split off of doOutputHeader because an
  // AC analysis can output both real data "doOutput" and complex
  // data via doOutputAC.  Each of these doOutput routines must
  // call a doHeader routine which can specify the right underlying
  // data type.  So rather than add more conditional logic we chose
  // to separate out the functionality for the different use cases.
  int numVars;
  std::vector< pair< std::string, char > > nT;
  NodeNamePairMap::iterator name_i, name_end;
  name_end = outputManager_.getAllNodes().end();

  if (outputManager_.getProcID() == 0)
  {
    if (!outputRAWTitleAndDate_)
    {
      // in multi-plot RAW files, the Title and Date are only output once
      // afer that each plot gets the same header.  Thus, only do this section
      // once.

      outputRAWTitleAndDate_ = true;

      (*outStreamPtr_) << "Title: " << outputManager_.getTitle() << std::endl;

      // create formatted timestamp
      const time_t now = time( NULL);
      char timeDate[ 40 ];
      strftime( timeDate, 40, "%a %b %d %I:%M:%S %Y", localtime( &now));
      (*outStreamPtr_) << "Date: " << timeDate << std::endl;
    }

    // format plot name
    std::ostringstream plotName;

    if (!outputManager_.getStepParamVec().empty())
    {
      plotName << "Step Analysis: Step " << outputManager_.getStepLoopNumber() + 1
               << " of " << outputManager_.getMaxParamSteps()
               << " params: ";
      std::vector<N_ANP_SweepParam>::const_iterator currItr = outputManager_.getStepParamVec().begin();
      std::vector<N_ANP_SweepParam>::const_iterator endItr = outputManager_.getStepParamVec().end();
      while ( currItr != endItr)
      {
        plotName << " name = " << currItr->name
                 << " value = " << currItr->currentVal << "  ";
        currItr++;
      }
    }
    if (!outputManager_.getDCParamVec().empty())
    {
      plotName << "DC Sweep: Step " << outputManager_.getDCLoopNumber() + 1
               << " of " << outputManager_.getMaxDCSteps()
               << " params: ";
      std::vector<N_ANP_SweepParam>::const_iterator currItr = outputManager_.getDCParamVec().begin();
      std::vector<N_ANP_SweepParam>::const_iterator endItr = outputManager_.getDCParamVec().end();
      while ( currItr != endItr)
      {
        plotName << " name = " << currItr->name
                 << " value = " << currItr->currentVal << "  ";
        currItr++;
      }
    }
    if (printParameters_.printType_ == PrintType::TRAN)
    {
      plotName << "Transient Analysis";
    }
    else if ( printParameters_.printType_ == PrintType::AC)
    {
      plotName << "AC Analysis";
    }
    else
    {
      plotName << "DC transfer characteristic";
    }

    (*outStreamPtr_) << "Plotname: " << plotName.str() << std::endl;

    // format the flags
    std::string flags("complex");
    (*outStreamPtr_) << "Flags: " << flags << std::endl;

  } // end proc0 check

  // prepare header for partial dump
  // get node types
  std::vector< char > typeRefs;
  outputManager_.getTopPtr()->returnVarTypeVec( typeRefs);

  // count data points
  numVars = outputManager_.getAllNodes().size();

#ifdef Xyce_PARALLEL_MPI

  int pos, s, e, len, final;
  char t;

  // store the total num of vars
  outputManager_.getCommPtr()->barrier();
  outputManager_.getCommPtr()->sumAll( &numVars, &final, 1);
  numVars = final;

  // limiting names to not longer than 64 chars +1(size) + 1(type)
  int bSize = sizeof( int) +( final *( 65 + sizeof( int)));
  char * b = new char [ bSize ];

  if (outputManager_.getProcID() == 0)
  {

#endif
    // store local names
    for ( name_i = outputManager_.getAllNodes().begin(); name_i != name_end ; ++name_i)
    {
      nT.push_back( pair< std::string, char>((*name_i).first,
                                             typeRefs[(*name_i).second.first]));
    }

#ifdef Xyce_PARALLEL_MPI

    // receive data points from procN in order
    for ( int p = 1; p < outputManager_.getNumProcs(); ++p)
    {
      pos = 0;

      outputManager_.getCommPtr()->recv( &s, 1, p);
      outputManager_.getCommPtr()->recv( b, s, p);
      outputManager_.getCommPtr()->unpack( b, bSize, pos, &e, 1);

      for ( int j = 0; j < e; ++j)
      {
        outputManager_.getCommPtr()->unpack( b, bSize, pos, &len, 1);

        std::string name( &(b[pos]), len);
        pos += len * sizeof( char);

        outputManager_.getCommPtr()->unpack( b, bSize, pos, &t, 1);

        // store remote names
        nT.push_back( pair< std::string, char >( name, t));
      }
    }
  } // end proc0 check

  else
  {
    pos = 0;

    e = outputManager_.getAllNodes().size();
    s = sizeof( int);
    outputManager_.getCommPtr()->pack( &e, 1, b, bSize, pos);

    for ( name_i = outputManager_.getAllNodes().begin(); name_i != name_end ; ++name_i)
    {
      len =((*name_i).first).size();
      s += sizeof( int);
      outputManager_.getCommPtr()->pack( &len, 1, b, bSize, pos);

      s += len * sizeof( char);
      outputManager_.getCommPtr()->pack(((*name_i).first).c_str(), len, b, bSize, pos);

      t = typeRefs[(*name_i).second.first];
      s += sizeof( char);
      outputManager_.getCommPtr()->pack( &t, 1, b, bSize, pos);
    }

    outputManager_.getCommPtr()->send( &s, 1, 0);
    outputManager_.getCommPtr()->send( b, s, 0);
  } // end procN check

    // clean up
  if (b != NULL)
  {
    free( b);
  }

#endif

  // format number of internal and external variables included here + time
  if (outputManager_.getProcID() == 0)
    (*outStreamPtr_) << "No. Variables: " << numVars + 1 << std::endl;

  outputManager_.getCommPtr()->barrier();

  if (outputManager_.getProcID() == 0)
  {
    // format total number of data points & remember the streampos
    (*outStreamPtr_) << "No. Points: ";
    numPointsPos_ = outStreamPtr_->tellp(); // <= 344 due to 80 char line limit
    (*outStreamPtr_) << "                  " << std::endl; // this will be overwritten

    // for full compatability with spice3 raw file, output version number here
    (*outStreamPtr_) << "Version: " << N_UTL_Version::getFullVersionString() << std::endl;

    // NOTE: dimensions, command, options, and scale not considered

    // write the variable information
    (*outStreamPtr_) << "Variables:" << std::endl;

    // The variable list appears next.  The format is:
    //         [tab](index) [tab](name) [tab](type)
    // Note that the(index) corresponds to the rawfile, not the soln vec.

    // write variable names for select data points
    std::string tmpNodeName, tmpType;
    if (printParameters_.printType_ == PrintType::TRAN)
    {
      (*outStreamPtr_) << "\t0\ttime\ttime\n";
    }
    else if (printParameters_.printType_ == PrintType::AC)
    {
      (*outStreamPtr_) << "\t0\tfrequency\tfrequency\n";
    }
    else
    {
      (*outStreamPtr_) << "\t0\tsweep\tvoltage\n";
    }

    std::string::size_type uPos;
    for ( int i = 0; i < numVars; ++i)
    {
      tmpNodeName =(nT[i]).first;

      // format is:  [tab](index) [tab](name) [tab](type)
      //   the index corresponds to the rawfile, not the soln vec
      if (strspn( tmpNodeName.c_str(), "0123456789") == tmpNodeName.size())
      {
        // sconvert, spice3f5, & chilespice wrap numeric voltage node names in V()
        tmpNodeName = "V(" +(nT[i]).first + ")";
      }

      uPos = tmpNodeName.rfind( "_", tmpNodeName.size()); if (uPos != std::string::npos)
                                                          {
                                                            tmpNodeName.replace( uPos, 1, "#");
                                                          }

      (*outStreamPtr_) << "\t" << i + 1 << "\t" << tmpNodeName << "\t";

      if ((nT[i]).second == 'I' ||(nT[i]).second == 'i')
      {
        (*outStreamPtr_) << "current\n" ;
      }
      else
      {
        (*outStreamPtr_) << "voltage\n";
      }

      // NOTE: other types, and params & dims, not considered
    }

    // write data marker
    (*outStreamPtr_) << "Values:" << std::endl;
  } // end proc0 check
}

void OverrideRawAscii::doOutputTime(const N_LAS_Vector * solnVecPtr, const N_LAS_Vector * stateVecPtr, const N_LAS_Vector * storeVecPtr)
{
  outputManager_.setCurrentOutputter(this);

  // write rawfile header on first call

  if (numPoints_ == 0)
  {
    doOpen();

    timeHeader();
  }

  std::ostream &os = *outStreamPtr_;

  // file IO on proc 0 only
  if (outputManager_.getProcID() == 0)
  {
    double tmp;

    // write the index to ascii rawfile
    os << numPoints_;

    // output the time/step values
    if ( printParameters_.printType_ == PrintType::TRAN) {
      tmp = outputManager_.getCircuitTime();
      tmp *= printParameters_.outputTimeScaleFactor_;
    }
    else
      tmp = outputManager_.getPRINTDCvalue();

    os << "\t" << tmp << "\n";
  }

  outputManager_.getCommPtr()->barrier();

  NodeNamePairMap::iterator name_i, name_end;
  name_end = outputManager_.getAllNodes().end();

#ifdef Xyce_PARALLEL_MPI

  int bSize, pos, s;
  int e = outputManager_.getAllNodes().size();

  outputManager_.getCommPtr()->barrier();
  outputManager_.getCommPtr()->maxAll( &e, &bSize, 1);

  bSize = bSize * sizeof( double);
  char * b = new char [ bSize ];

  if (outputManager_.getProcID() == 0)
  {

#endif

    // dump proc0 data points
    for ( name_i = outputManager_.getAllNodes().begin(); name_i != name_end ; ++name_i)
    {
      double tmp =(*solnVecPtr)[(*name_i).second.first];
      if (fabs( tmp) < outputManager_.getFilter())
      {
        tmp = 0.0;
      }

      // write formatted values to rawfile
      os << "\t" << tmp << "\n";
    }

#ifdef Xyce_PARALLEL_MPI

    // receive data points from procN in order
    for ( int p = 1; p < outputManager_.getNumProcs(); ++p)
    {
      outputManager_.getCommPtr()->recv( &e, 1, p);
      s = e * sizeof( double);
      pos = 0;
      outputManager_.getCommPtr()->recv( b, s, p);
      for ( int j = 0; j < e; j++)
      {
        double tmp = 0.0;
        outputManager_.getCommPtr()->unpack( b, bSize, pos, &tmp, 1);

        // write formatted values to rawfile
        os << "\t" << tmp << "\n";
      }
    }

  } // end proc0 work

  else
  {
    // pack and send local data points to proc0
    s = e * sizeof( double);
    pos = 0;

    for ( name_i = outputManager_.getAllNodes().begin(); name_i != name_end ; ++name_i)
    {
      double tmp = (*solnVecPtr)[(*name_i).second.first];
      if (fabs( tmp) < outputManager_.getFilter())
      {
        tmp = 0.0;
      }

      outputManager_.getCommPtr()->pack( &tmp, 1, b, bSize, pos);
    }

    outputManager_.getCommPtr()->send( &e, 1, 0);
    outputManager_.getCommPtr()->send( b, s, 0);
  }

  // clean up
  delete b;
#endif

  if (outputManager_.getProcID() == 0)
  {
    // write newline to rawfile
    os << std::endl;
  }

  outputManager_.getCommPtr()->barrier();

  // keep track of number of datapoints
  ++numPoints_;
}

void
OverrideRawAscii::doOutputFrequency(double frequency, const N_LAS_Vector *real_solution_vector, const N_LAS_Vector *imaginary_solution_vector)
{
  outputManager_.setCurrentOutputter(this);
  double tmp;

  if (numPoints_ == 0)
  {
    doOpen();

    frequencyHeader();
  }

  // file IO on proc 0 only
  if (outputManager_.getProcID() == 0)
  {
    // write the index to ascii rawfile
    (*outStreamPtr_) << numPoints_;
  }

  outputManager_.getCommPtr()->barrier();

  tmp=frequency;

  // file IO only on proc 0
  if (outputManager_.getProcID() == 0)
  {
    (*outStreamPtr_) << "\t" << frequency << ", " << 0.0 << "\n";
  }

  NodeNamePairMap::iterator name_i, name_end;
  name_end = outputManager_.getAllNodes().end();

#ifdef Xyce_PARALLEL_MPI

  int bSize, pos, s;
  int e = outputManager_.getAllNodes().size();

  outputManager_.getCommPtr()->barrier();
  outputManager_.getCommPtr()->maxAll( &e, &bSize, 1);

  bSize = bSize * sizeof( double);
  char * b = new char [ bSize ];

  if (outputManager_.getProcID() == 0)
  {

#endif

    // dump proc0 data points
    for ( name_i = outputManager_.getAllNodes().begin(); name_i != name_end ; ++name_i)
    {
      tmp =(*real_solution_vector)[(*name_i).second.first];
      double tmpImg =(*imaginary_solution_vector)[(*name_i).second.first];
      if (fabs( tmp) < outputManager_.getFilter())
      {
        tmp = 0.0;
      }
      if (fabs( tmpImg) < outputManager_.getFilter())
      {
        tmpImg = 0.0;
      }

      // write formatted values to rawfile
      (*outStreamPtr_) << "\t" << tmp << ", " << tmpImg << "\n";
    }

#ifdef Xyce_PARALLEL_MPI

    // receive data points from procN in order
    for ( int p = 1; p < outputManager_.getNumProcs(); ++p)
    {
      outputManager_.getCommPtr()->recv( &e, 1, p);
      s = e * sizeof( double);
      pos = 0;
      outputManager_.getCommPtr()->recv( b, s, p);
      for ( int j = 0; j < e; j++)
      {
        outputManager_.getCommPtr()->unpack( b, bSize, pos, &tmp, 1);

        // write formatted values to rawfile
        (*outStreamPtr_) << "\t" << tmp << "\n";
      }
    }

  } // end proc0 work

  else
  {
    // pack and send local data points to proc0
    s = e * sizeof( double);
    pos = 0;

    for ( name_i = outputManager_.getAllNodes().begin(); name_i != name_end ; ++name_i)
    {
      tmp =(*real_solution_vector)[(*name_i).second.first];
      if (fabs( tmp) < outputManager_.getFilter())
      {
        tmp = 0.0;
      }

      outputManager_.getCommPtr()->pack( &tmp, 1, b, bSize, pos);
    }

    outputManager_.getCommPtr()->send( &e, 1, 0);
    outputManager_.getCommPtr()->send( b, s, 0);
  }

  // clean up
  delete b;
#endif

  if (outputManager_.getProcID() == 0)
  {
    // write newline to rawfile
    (*outStreamPtr_) << std::endl;
  }

  // sync after writing data for step
  outputManager_.getCommPtr()->barrier();

  // keep track of number of datapoints
  ++numPoints_;
}

void OverrideRawAscii::doOutputHB(
  const std::vector<double>&    timePoints,
  const std::vector<double>&    freqPoints,
  const N_LAS_BlockVector &     timeDomainSolnVec,
  const N_LAS_BlockVector &     freqDomainSolnVecReal,
  const N_LAS_BlockVector &     freqDomainSolnVecImaginary)
{}

void OverrideRawAscii::doOutputMPDE(double time, const N_LAS_Vector *solution_vector)
{}

void OverrideRawAscii::doOutputHomotopy(const std::vector<std::string> & parameter_names, const std::vector<double> & parameter_values, const N_LAS_Vector * solution_vector)
{}

void OverrideRawAscii::doOutputMORTF(bool origSystem, const double & freq, const Teuchos::SerialDenseMatrix<int, std::complex<double> >& H)
{}

void OverrideRawAscii::doResetOutput()
{}

void OverrideRawAscii::doFinishOutput()
{
  if (outputManager_.getProcID() == 0)
  {
    if (outStreamPtr_) {
      // need to move file pointer back to header and
      // write out the number of points.
      long currentFelePost = outStreamPtr_->tellp();

      // locate the position for number of points
      outStreamPtr_->seekp( numPointsPos_);

      // overwrite blank space with value
      (*outStreamPtr_) << numPoints_;

      // move file pointer to the end again.
      outStreamPtr_->seekp( currentFelePost);
    }
  }

  // reset numPoints_ as it is used as a flag
  // by outputRAW_() to print the header.
  numPoints_=0;
}


void
OverrideRawAscii::doFinishOutputStep()
{}

MOR::MOR(OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    firstTimePrint_(true),
    index_(0),
    outFilename_(),
    suffix_(),
    outStreamPtr_(0),
    headerPrintCalls_(0)
{
  if (printParameters_.extension_.empty())
    printParameters_.extension_ = ".prn";
}

MOR::~MOR()
{
  if (outStreamPtr_ != &std::cout)
    delete outStreamPtr_;
}

void MOR::doParse()
{
  if (openOrig_)
  {
    outFilename_ = outputManager_.getNetListFilename() + ".Orig.FD.prn";
  }
  else
  {
    outFilename_ = outputManager_.getNetListFilename() + ".Red.FD.prn";
  }
}

void MOR::doOpen()
{
  if (stepCount_ == 0)
  {
    outStreamPtr_ = new std::ofstream(outFilename_.c_str());
  }
}

//-----------------------------------------------------------------------------
// Function      : MOR::outputMORHeaders
// Purpose       : .PRINT file header output for MOR runs.
//
// Special Notes : This function assumes the "regsiterPRINTSet" was already
//                 called, so the PRINTblock_ and format_ are already
//                 set up.
//
// Scope         : public
// Creator       : Heidi Thornquist, Ting Mei
// Creation Date : 5/25/12
//-----------------------------------------------------------------------------
void MOR::outputMORHeaders(int numPorts)
{
  index_ = 0;

#ifdef Xyce_DEBUG_IO
  cout << "MOR Header output:" << std::endl;
#endif

  // hardwire the AC print format to tecplot or gnuplot(not PROBE!)
  if (format_ != Format::TECPLOT && format_ != Format::STD)
  {
    format_ = Format::STD;
    string msg("MOR output can only use tecplot or gnuplot format.");
    msg += "  Resetting to gnuplot format.\n";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_WARNING, msg);
  }

  if (format_ == Format::TECPLOT) // tecplot header
  {
    // STD header
    // Freq Domain Headers
    tecplotFreqMORHeader_(*outStreamPtr_, stepCount_, numPorts);
  }
  else
  {
    // STD header
    // Freq Domain Headers
    stdFreqMORHeader_(*outStreamPtr_, numPorts);
  }

  outStreamPtr_->setf(ios::scientific);
  outStreamPtr_->precision(printParameters_.streamPrecision_);
  outStreamPtr_->setf(ios::left, ios::adjustfield);

  ++stepCount_;
}



//-----------------------------------------------------------------------------
// Function      : MOR::outputMORTF
// Purpose       : .PRINT output for mor runs for original and reduced system
// Special Notes :
// Scope         : public
// Creator       : Heidi Thornquist, Ting Mei
// Creation Date : 5/25/12
//-----------------------------------------------------------------------------
void MOR::doOutputMORTF(bool origSystem, const double & freq, const Teuchos::SerialDenseMatrix<int, std::complex<double> >& H)
{
  bool firstColPrinted = false;

#ifdef Xyce_DEBUG_IO
  cout << "Begin outputMORTF" << std::endl;
#endif

  std::ostringstream ost;
  ost << outputManager_.getProcID();

  if (outputManager_.getProcID() == 0)
  {
    if (firstTimePrint_) //Setup Output Stream and Print Out Header
    {
      doOpen();

      outputMORHeaders(H.numRows());

      firstTimePrint_ = false;
    }

    if (printParameters_.delimiter_ == "")
    {
      outStreamPtr_->width(printParameters_.streamWidth_);
    }

    // Output a TAB.
    if (printParameters_.delimiter_ != "")
    {
      (*outStreamPtr_) << printParameters_.delimiter_;
    }
  } // procID

  // Get the indices:
  // First print out the Index:
  if (index_ && format_ != Format::TECPLOT)
  {
    int W = getWidthFromStaticIndex(index_, printParameters_.delimiter_);
    outStreamPtr_->width(W);
    (*outStreamPtr_) << index_++;
    firstColPrinted = true;
  }

  if (outputManager_.getProcID() == 0)
  {
    //-------------------------------------
    // Get the time value first.
    //-------------------------------------
    if (printParameters_.delimiter_ == "")
    {
      outStreamPtr_->width(printParameters_.streamWidth_);
    }
    else
    {
      outStreamPtr_->width(0);
      if (firstColPrinted) {
        (*outStreamPtr_) << printParameters_.delimiter_;
      }
    }

    (*outStreamPtr_) << freq;
    firstColPrinted = true;
  }

  {
    int N = H.numRows();
    for (int i=0; i<N; ++i)
    {
      for (int j=0; j<N; ++j)
      {
        if (outputManager_.getProcID() == 0)
        {
          if (printParameters_.delimiter_ == "")
          {
            outStreamPtr_->width(printParameters_.streamWidth_);
          }
          else
          {
            outStreamPtr_->width(0);
            if (printParameters_.delimiter_ != "") {
              (*outStreamPtr_) << printParameters_.delimiter_;
            }
          }
          (*outStreamPtr_) << H(i, j).real();

          if (printParameters_.delimiter_ == "")
          {
            outStreamPtr_->width(printParameters_.streamWidth_);
          }
          else
          {
            outStreamPtr_->width(0);
            if (printParameters_.delimiter_ != "") {
              (*outStreamPtr_) << printParameters_.delimiter_;
            }
          }
          (*outStreamPtr_) << H(i, j).imag();

        }
        outputManager_.getCommPtr()->barrier();
      }
    }
  }

  if (outputManager_.getProcID() == 0)
  {
    (*outStreamPtr_) << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : MOR::stdFreqMORHeader_
// Purpose       : header for std. Freq Domain MOR(default)
// Special Notes :
// Scope         : private
// Creator       : Heidi Thornquist, Ting Mei
// Creation Date : 5/31/12
//-----------------------------------------------------------------------------
void MOR::stdFreqMORHeader_(ostream & stream, int numPorts)
{
  Table table;
  Table::Justification freqJust = Table::JUSTIFICATION_LEFT;
  if (index_)
  {
    table.addColumn("Index", printParameters_.timeWidth_- 1, printParameters_.streamPrecision_, Table::JUSTIFICATION_LEFT);
    freqJust = Table::JUSTIFICATION_CENTER;
  }
  table.addColumn("FREQ", printParameters_.streamWidth_, printParameters_.streamPrecision_, freqJust);

  // Assume that the input/output ports are the same, making a square H(s)
  for (int i=0; i<numPorts; ++i)
  {
    for (int j=0; j<numPorts; ++j)
    {
      std::ostringstream realName;
      std::ostringstream imaginaryName;
      realName << "Re(H(" << i << ", " << j << "))";
      table.addColumn(realName.str(), printParameters_.streamWidth_, printParameters_.streamPrecision_, Table::JUSTIFICATION_CENTER);
      imaginaryName << "Im(H(" << i << ", " << j << "))";
      table.addColumn(imaginaryName.str(), printParameters_.streamWidth_, printParameters_.streamPrecision_, Table::JUSTIFICATION_CENTER);
    }
  }
  // std::vector<std::string> paddedHeaderNames = stdHeaderMangler(table, printParameters_.delimiter_);
  // int N = paddedHeaderNames.size();
  // for (int i=0 ; i < N ; ++i)
  // {
  //   stream << paddedHeaderNames[i];
  // }
  stream << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : MOR::tecplotFreqMORHeader_
// Purpose       : header for tecplot. Freq Domain for MOR(default)
// Special Notes :
// Scope         : private
// Creator       : Heidi Thornquist and Ting Mei
// Creation Date : 5/31/12
//-----------------------------------------------------------------------------
void MOR::tecplotFreqMORHeader_(ostream & stream, int counter, int numPorts)
{
  if (counter == 0)
  {
    stream << " TITLE = \" Xyce Frequency Domain data, " << outputManager_.getNetListFilename() << "\", " << std::endl;
    stream << "\tVARIABLES = \"FREQ \" " << std::endl;

    // Assume that the input/output ports are the same, making a square H(s)
    for (int i=0; i<numPorts; ++i)
    {
      for (int j=0; j<numPorts; ++j)
      {
        std::ostringstream realName;
        std::ostringstream imaginaryName;
        realName << "Re(H(" << i << ", " << j << "))";
        imaginaryName << "Im(H(" << i << ", " << j << "))";
        stream << "\" ";
        stream << realName.str();
        stream << "\" " << std::endl;
        stream << "\" ";
        stream << imaginaryName.str();
        stream << "\" " << std::endl;
      }
    }

    // output some AUXDATA
    stream << "DATASETAUXDATA ";
    stream << getTecplotTimeDateStamp();
    stream << std::endl;
    stream << "ZONE F=POINT\n";
  }
  stream << std::endl;
  stream << flush;
}

//-----------------------------------------------------------------------------
// Function      : MOR::outputPRINT_
// Purpose       : .PRINT FORMAT=STD time domain output
// Special Notes :
// Scope         : private
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 10/10/00
//-----------------------------------------------------------------------------
void MOR::doOutputTime(const N_LAS_Vector * solnVecPtr, const N_LAS_Vector * stateVecPtr, const N_LAS_Vector * storeVecPtr)
{
}

void
MOR::doOutputFrequency(double frequency, const N_LAS_Vector *real_solution_vector, const N_LAS_Vector *imaginary_solution_vector)
{
}

void MOR::doOutputHB(
  const std::vector<double>&    timePoints,
  const std::vector<double>&    freqPoints,
  const N_LAS_BlockVector &     timeDomainSolnVec,
  const N_LAS_BlockVector &     freqDomainSolnVecReal,
  const N_LAS_BlockVector &     freqDomainSolnVecImaginary)
{}

void MOR::doOutputMPDE(double time, const N_LAS_Vector *solution_vector)
{}

void MOR::doOutputHomotopy(const std::vector<std::string> & parameter_names, const std::vector<double> & parameter_values, const N_LAS_Vector * solution_vector)
{}

void MOR::doFinishOutput()
{}

void MOR::doResetOutput() {
  firstTimePrint_ = true;
  stepCount_ = 0;
}

void MOR::doFinishOutputStep()
{
}

} // namespace outputter

namespace { // <unnamed>

//-----------------------------------------------------------------------------
// Function      : MOR::tecplotFreqHeader_
// Purpose       : header for tecplot. Freq Domain(default)
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter
// Creation Date : 4/11/12
//-----------------------------------------------------------------------------
void tecplotFreqHeader(OutputMgr &output_manager, const PrintParameters &print_parameters, std::ostream &os, int counter)
{
  if (counter == 0)
  {
    os << " TITLE = \" Xyce Frequency Domain data, " << output_manager.getNetListFilename() << "\", " << std::endl;
    os << "\tVARIABLES = \"FREQ \" " << std::endl;

    // output the user-specified solution vars:
    for (ParameterList::const_iterator iterParam = print_parameters.variableList_.begin() ; iterParam != print_parameters.variableList_.end(); ++iterParam)
    {
      os << "\" ";
      os << "Re(" << (*iterParam).tag() << ")";
      os << "\" " << std::endl;
      os << "\" ";
      os << "Im(" << (*iterParam).tag() << ")";
      os << "\" " << std::endl;
    }
  }

  // output some AUXDATA
  os << "DATASETAUXDATA ";
  os << getTecplotTimeDateStamp();
  os << std::endl;
  os << "ZONE F=POINT\n";

  if (output_manager.getStepParamVec().empty())
  {
    os << " T=\"Xyce data\" ";
  }
  else
  {
    std::vector<N_ANP_SweepParam>::const_iterator iterParam;
    std::vector<N_ANP_SweepParam>::const_iterator firstParam=output_manager.getStepParamVec().begin();
    std::vector<N_ANP_SweepParam>::const_iterator lastParam=output_manager.getStepParamVec().end();

    os << " T= \" ";
    for (iterParam=firstParam;iterParam!=lastParam;++iterParam)
    {
      int tecplotHeaderPrecision_ = 2;
      os.setf(ios::scientific);
      os.precision( tecplotHeaderPrecision_);
      os << " " << iterParam->name << " = ";
      os << iterParam->currentVal;
    }
    os << "\" ";
  }

  os << std::endl;
}

std::ostream &printHeader(std::ostream &os, const PrintParameters &print_parameters)
{
  return printHeader(os, print_parameters.table_.columnList_, print_parameters.delimiter_);
}


std::ostream &printHeader(std::ostream &os, const Table::ColumnList &column_list, const std::string &delimiter)
{
  for (Table::ColumnList::const_iterator it = column_list.begin(); it != column_list.end(); ++it) {
    if (it != column_list.begin())
      os << (delimiter.empty() ? " " : delimiter);

    printHeader(os, (*it));
  }

  os << std::endl;

  return os;
}


std::ostream &printHeader(std::ostream &os, const Table::Column &column)
{
  std::string name = column.name_;
  if (name == "INDEX")
    name = "Index";

  size_t left_padding = 0;
  size_t right_padding = 0;

  // if (column.width_ < name.size())
  //   column.width_ = name.size();

  if (column.width_ > name.size()) {
    switch (column.justification_) {
      case Table::JUSTIFICATION_LEFT:
        right_padding = column.width_ - left_padding - name.size();
        break;
      case Table::JUSTIFICATION_CENTER:
        left_padding = (column.width_ - column.name_.size())/2;
        right_padding = column.width_ - left_padding - name.size();
        break;
      case Table::JUSTIFICATION_RIGHT:
        left_padding = column.width_ - column.name_.size();
        break;
    }
  }

  os << std::setw(left_padding) << "" << std::setw(0) << name << std::setw(right_padding) << "";

  return os;
}

void fixupColumns(const OutputMgr &output_manager, PrintParameters &print_parameters)
{
  Table::Justification justification = print_parameters.delimiter_.empty() ? Table::JUSTIFICATION_CENTER :  Table::JUSTIFICATION_NONE;

  for (ParameterList::const_iterator it = print_parameters.variableList_.begin() ; it != print_parameters.variableList_.end(); ++it)
  {
    switch ((*it).getSimContext())
    {
      case INDEX:
        print_parameters.table_.addColumn("INDEX", ios_base::fixed, 5, 0, Table::JUSTIFICATION_LEFT);
        break;

      case TIME_VAR:
        print_parameters.table_.addColumn("TIME", print_parameters.streamWidth_, print_parameters.streamPrecision_, justification);
        break;

      case FREQUENCY:
        print_parameters.table_.addColumn("FREQ", print_parameters.streamWidth_, print_parameters.streamPrecision_, justification);
        break;

      default:
        print_parameters.table_.addColumn((*it).tag(), print_parameters.streamWidth_, print_parameters.streamPrecision_, justification);
        break;
    }
  }
}

std::ostream &printValue(std::ostream &os, const Table::Column &column, const std::string &delimiter, const int column_index, double value)
{
  if (delimiter.empty()) {
    if (column_index != 0)
      os <<  " ";
    os << std::resetiosflags(ios_base::floatfield) << std::setiosflags(column.format_)
       << std::resetiosflags(ios_base::adjustfield) << std::setiosflags(column.justification_ == Table::JUSTIFICATION_LEFT ? ios_base::left : ios_base::right)
       << std::setprecision(column.precision_) << std::setw(column.width_)
       << value;
  }
  else {
    if (column_index != 0)
      os << delimiter;
    os << std::resetiosflags(ios_base::floatfield) << std::setiosflags(column.format_) << std::setw(0) << std::setprecision(column.precision_) << value;
  }

  return os;
}

} // namespace <unnamed>
} // namespace IO
} // namespace Xyce

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
// Filename       : $RCSfile: N_IO_Outputter.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Dave Baur
//
// Creation Date  :
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.53.2.3 $
//
// Revision Date  : $Date: 2014/03/18 21:43:40 $
//
// Current Owner  : $Author: erkeite $
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
#include <N_IO_Op.h>

#include <N_ANP_AnalysisInterface.h>
#include <N_ANP_AnalysisManager.h>
#include <N_MPDE_Manager.h>

#include <N_UTL_Version.h>

#include <N_LAS_Vector.h>
#include <N_LAS_BlockVector.h>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Namespace     : Unnamed
// Purpose       : file-local scoped methods and data
// Special Notes :
// Creator       : dgbaur
// Creation Date : 06/28/2013
//-----------------------------------------------------------------------------
namespace {

//-----------------------------------------------------------------------------
// Function      : deleteList
// Purpose       : templated function to delete items from container given
//                 begin and end iterators
// Special Notes :
// Scope         : file-local
// Creator       : David Baur, Raytheon
// Creation Date : 11/21/2013
//-----------------------------------------------------------------------------
template <class II>
void deleteList(II begin, II end)
{
  for (; begin != end; ++begin)
    delete *begin;
}

//-----------------------------------------------------------------------------
// Function      : getTimeDateStamp
// Purpose       : get current date and time and format for .PRINT output
// Special Notes :
// Scope         : file-local
// Creator       : David Baur, Raytheon
// Creation Date : 06/14/2013
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
// Function      : getTecplotTimeDateStamp
// Purpose       : Get current date and time and format for .PRINT output
// Special Notes : tecplot version of getTimeDateStamp.
// Scope         : file-local
// Creator       : David Baur, Raytheon
// Creation Date : 6/14/2013
//-----------------------------------------------------------------------------
std::string getTecplotTimeDateStamp()
{
  const time_t now = time( NULL);
  char timeDate[ 40 ];

  // format for output
  strftime( timeDate, 40, "TIME= \" %I:%M:%S %p %b %d, %Y \" ", localtime( &now));

  return std::string( timeDate);
}

//-----------------------------------------------------------------------------
//The following are declarations of functions in the unnamed namespace
//that are defined in a second namespace{} block at the end of this
//file
//-----------------------------------------------------------------------------
void tecplotFreqHeader(OutputMgr &output_manager,
                       const PrintParameters &print_parameters,
                       const Util::OpList &op_list, std::ostream &os,
                       int counter);
void fixupColumns(const OutputMgr &output_manager,
                  PrintParameters &print_parameters, Util::OpList &op_list);

std::ostream &printHeader(std::ostream &os,
                          const PrintParameters &print_parameters);
std::ostream &printHeader(std::ostream &os,
                          const Table::ColumnList &column_list,
                          const std::string &delimiter);
std::ostream &printHeader(std::ostream &os, const Table::Column &column);
std::ostream &printValue(std::ostream &os, const Table::Column &column,
                         const std::string &delimiter, const int column_index,
                         double value);
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// Function      : outputFilename
// Purpose       : return a string for the name of the results output file
// Special Notes : uses the name on the .print line as the base filename if
//                 given, otherwise uses the netlist name as the base.
//                 Adds suffixes or extensions as provided by analysis type.
// Scope         : file-local
// Creator       : David Baur, Raytheon
// Creation Date : 6/14/2013
//-----------------------------------------------------------------------------
std::string outputFilename(const PrintParameters &print_parameters, const std::string &netListFilename)
{
  if (!print_parameters.filename_.empty())
  {
    return print_parameters.filename_ + print_parameters.suffix_ + print_parameters.extraExtension_;
  }
  else
  {
    return netListFilename + print_parameters.suffix_ + print_parameters.defaultExtension_;
  }
}

//-----------------------------------------------------------------------------
// Function      : filter
// Purpose       : Applies a filter to a double value, returning 0 if the
//                 absolute value of the data is less than the filter.
// Special Notes :
// Scope         : file-local
// Creator       : David Baur, Raytheon
// Creation Date : 11/25/2013
//-----------------------------------------------------------------------------
inline double filter(double value, double filter)
{
  return std::abs(value) < filter ? 0.0 : value;
}

} // namespace <unnamed>

namespace Outputter {

//-----------------------------------------------------------------------------
// Class         : TimePrn
// Purpose       : Outputter class for transient runs and "PRN" (Standard)
//                 output format
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Function      : TimePrn::TimePrn
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 06/07/2013
//-----------------------------------------------------------------------------
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
  if (printParameters_.defaultExtension_.empty())
    printParameters_.defaultExtension_ = ".prn";

  fixupColumns(outputManager_, printParameters_, opList_);
}

//-----------------------------------------------------------------------------
// Function      : TimePrn::~TimePrn
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 06/07/2013
//-----------------------------------------------------------------------------
TimePrn::~TimePrn()
{
  outputManager_.closeFile(outStreamPtr_);

  deleteList(opList_.begin(), opList_.end());
}

//-----------------------------------------------------------------------------
// Function      : TimePrn::doParse
// Purpose       : create the output file name based on print parameters and
//                 netlist file name.
// Special Notes : Does no parsing whatsoever.  Vestigal name from when this
//                 method still read the command line args itself.
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
void TimePrn::doParse()
{
  outFilename_ = outputFilename(printParameters_, outputManager_.getNetListFilename());
}

//-----------------------------------------------------------------------------
// Function      : TimePrn::doOpen
// Purpose       : Open the output file
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
void TimePrn::doOpen()
{
  if (outputManager_.getProcID() == 0 && outStreamPtr_ == 0)
  {
    outStreamPtr_ = outputManager_.openFile(outFilename_);
  }
}

//-----------------------------------------------------------------------------
// Function      : TimePrn::timeHeader
// Purpose       : output the header
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013 ?
//-----------------------------------------------------------------------------
void TimePrn::timeHeader()
{
  if (outputManager_.getProcID() == 0 && headerPrintCalls_ == 0)
  {
    printHeader(*outStreamPtr_, printParameters_);
  }
}

//-----------------------------------------------------------------------------
// Function      : TimePrn::doOutputTime
// Purpose       : Output the current data at a time point
// Special Notes :
// Scope         : public
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
void TimePrn::doOutputTime(
    const N_LAS_Vector * solnVecPtr, 
    const N_LAS_Vector * stateVecPtr, 
    const N_LAS_Vector * storeVecPtr)
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
  for (Util::OpList::const_iterator it = opList_.begin() ; it != opList_.end(); ++it, ++column_index)
  {
    double result = getValue(outputManager_.getCommPtr()->comm(), *(*it), solnVecPtr, 0, stateVecPtr, storeVecPtr).real();
    result = filter(result, printParameters_.filter_);
    if ((*it)->opType() == Util::TIME_VAR)
      result *= printParameters_.outputTimeScaleFactor_;

    if (outputManager_.getProcID() == 0)
      printValue(os, printParameters_.table_.columnList_[column_index], printParameters_.delimiter_, column_index, result);

    outputManager_.getCommPtr()->barrier();
  }

  ++index_;

  if (outputManager_.getProcID() == 0)
    os << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : TimePrn::doFinishOutput
// Purpose       : Output the footer, close stream
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
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

//-----------------------------------------------------------------------------
// Function      : TimePrn::doFinishOutputStep
// Purpose       : output footer and close stream for parameter sweep file
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
void TimePrn::doFinishOutputStep()
{
  // Deal with the *prn file:
  if (outStreamPtr_)
  {
    if ( outputManager_.getPrintEndOfSimulationLine() )
    {
      (*outStreamPtr_) << "End of Xyce(TM) Parameter Sweep" << std::endl;
    }

    outputManager_.closeFile(outStreamPtr_);
    outStreamPtr_ = 0;
  }
}


//-----------------------------------------------------------------------------
// Class         : FrequencyPrn
// Purpose       : Outputter class for frequency domain data in "PRN" (Standard)
//                 output format
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : FrequencyPrn::FrequencyPrn
// Purpose       : constructor
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
FrequencyPrn::FrequencyPrn(OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    firstTimePrint_(true),
    outFilename_(),
    suffix_(),
    outStreamPtr_(0),
    stepCount_(0)
{
  if (printParameters_.defaultExtension_.empty())
    printParameters_.defaultExtension_ = ".FD.prn";

  fixupColumns(outputManager_, printParameters_, opList_);
}

//-----------------------------------------------------------------------------
// Function      : FrequencyPrn::~FrequencyPrn
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
FrequencyPrn::~FrequencyPrn()
{
  outputManager_.closeFile(outStreamPtr_);

  deleteList(opList_.begin(), opList_.end());
}

//-----------------------------------------------------------------------------
// Function      : FrequencyPrn::doParse
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void FrequencyPrn::doParse()
{
  outFilename_ = outputFilename(printParameters_, outputManager_.getNetListFilename());
}

//-----------------------------------------------------------------------------
// Function      : FrequencyPrn::doOpen
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void FrequencyPrn::doOpen()
{
  if (outputManager_.getProcID() == 0 && outStreamPtr_ == 0)
  {
    outStreamPtr_ = outputManager_.openFile(outFilename_);
  }
}

//-----------------------------------------------------------------------------
// Function      : FrequencyPrn::doOutputFrequency
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void FrequencyPrn::doOutputFrequency(double frequency, const N_LAS_Vector *real_solution_vector, const N_LAS_Vector *imaginary_solution_vector)
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
  for (Util::OpList::const_iterator it = opList_.begin(); it != opList_.end(); ++it, ++column_index)
  {
    double result = getValue(outputManager_.getCommPtr()->comm(), *(*it), real_solution_vector, imaginary_solution_vector, 0, 0).real();
    if (outputManager_.getProcID() == 0)
      printValue(os, printParameters_.table_.columnList_[column_index], printParameters_.delimiter_, column_index, result);

    outputManager_.getCommPtr()->barrier();
  }

  ++index_;

  if (outputManager_.getProcID() == 0)
    os << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : FrequencyPrn::doFinishOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void FrequencyPrn::doFinishOutput()
{
  if (outputManager_.getProcID() == 0)
  {
    if (outStreamPtr_)
    {
      if (!outputManager_.getSTEPEnabledFlag() && outputManager_.getPrintEndOfSimulationLine())
      {
        firstTimePrint_ = true;
        (*outStreamPtr_) << "End of Xyce(TM) Simulation" << std::endl;
        firstTimePrint_ = true;
      }

      if (!outputManager_.getSTEPEnabledFlag())
      {
        outputManager_.closeFile(outStreamPtr_);
        outStreamPtr_ = 0;
      }
    }
  } // procID

  index_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : FrequencyPrn::doFinishOutputStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void FrequencyPrn::doFinishOutputStep()
{
  // Deal with the *prn file:
  if (outStreamPtr_)
  {
    if ( outputManager_.getPrintEndOfSimulationLine() )
    {
      (*outStreamPtr_) << "End of Xyce(TM) Parameter Sweep" << std::endl;
    }

    outputManager_.closeFile(outStreamPtr_);
    outStreamPtr_ = 0;
  }
}

//-----------------------------------------------------------------------------
// Class         : TimeCSV
// Purpose       : Outputter class for transient runs and "CSV"
//                 (comma separarated variable) output format
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : TimeCSV::TimeCSV
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
TimeCSV::TimeCSV(OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    firstTimePrint_(true),
    outFilename_(),
    suffix_(),
    outStreamPtr_(0),
    headerPrintCalls_(0)
{
  if (printParameters_.defaultExtension_.empty())
    printParameters_.defaultExtension_ = ".csv";

  fixupColumns(outputManager_, printParameters_, opList_);
}

//-----------------------------------------------------------------------------
// Function      : TimeCSV::~TimeCSV
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
TimeCSV::~TimeCSV()
{
  outputManager_.closeFile(outStreamPtr_);

  deleteList(opList_.begin(), opList_.end());
}


//-----------------------------------------------------------------------------
// Function      : TimeCSV::doParse
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void TimeCSV::doParse()
{
  outFilename_ = outputFilename(printParameters_, outputManager_.getNetListFilename());
}

//-----------------------------------------------------------------------------
// Function      : TimeCSV::doOpen
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void TimeCSV::doOpen()
{
  if (outputManager_.getProcID() == 0 && outStreamPtr_ == 0)
  {
    outStreamPtr_ = outputManager_.openFile(outFilename_);
  }
}

//-----------------------------------------------------------------------------
// Function      : TimeCSV::timeHeader
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
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
void TimeCSV::doOutputTime(
    const N_LAS_Vector * solnVecPtr, 
    const N_LAS_Vector * stateVecPtr, 
    const N_LAS_Vector * storeVecPtr)
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
  for (Util::OpList::const_iterator it = opList_.begin() ; it != opList_.end(); ++it, ++column_index)
  {
    double result = getValue(outputManager_.getCommPtr()->comm(), *(*it), solnVecPtr, 0, stateVecPtr, storeVecPtr).real();
    result = filter(result, printParameters_.filter_);
    if ((*it)->opType() == Util::TIME_VAR)
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

//-----------------------------------------------------------------------------
// Function      : TimeCSV::doFinishOutputStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void TimeCSV::doFinishOutputStep()
{
  if (outStreamPtr_)
  {
    outputManager_.closeFile(outStreamPtr_);
    outStreamPtr_ = 0;
  }
}

//-----------------------------------------------------------------------------
// Class         : FrequencyCSV
// Purpose       : Outputter class for frequency-domain runs and "CSV"
//                 (comma separarated variable) output format
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : FrequencyCSV::FrequencyCSV
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
FrequencyCSV::FrequencyCSV(OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    firstTimePrint_(true),
    outFilename_(),
    suffix_(),
    outStreamPtr_(0),
    stepCount_(0)
{
  if (printParameters_.defaultExtension_.empty())
    printParameters_.defaultExtension_ = ".FD.csv";

  fixupColumns(outputManager_, printParameters_, opList_);
}

//-----------------------------------------------------------------------------
// Function      : FrequencyCSV::~FrequencyCSV
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
FrequencyCSV::~FrequencyCSV()
{
  outputManager_.closeFile(outStreamPtr_);

  deleteList(opList_.begin(), opList_.end());
}


//-----------------------------------------------------------------------------
// Function      : FrequencyCSV::doParse
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void FrequencyCSV::doParse()
{
  outFilename_ = outputFilename(printParameters_, outputManager_.getNetListFilename());
}

//-----------------------------------------------------------------------------
// Function      : FrequencyCSV::doOpen
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void FrequencyCSV::doOpen()
{
  if (outputManager_.getProcID() == 0 && outStreamPtr_ == 0)
  {
    outStreamPtr_ = outputManager_.openFile(outFilename_);
  }
}

//-----------------------------------------------------------------------------
// Function      : FrequencyCSV::doOutputFrequency
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void FrequencyCSV::doOutputFrequency(double frequency, const N_LAS_Vector *real_solution_vector, const N_LAS_Vector *imaginary_solution_vector)
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
  for (Util::OpList::const_iterator it = opList_.begin(); it != opList_.end(); ++it, ++column_index)
  {
    double result = getValue(outputManager_.getCommPtr()->comm(), *(*it), real_solution_vector, imaginary_solution_vector, 0, 0).real();
    if (outputManager_.getProcID() == 0)
      printValue(os, printParameters_.table_.columnList_[column_index], printParameters_.delimiter_, column_index, result);

    outputManager_.getCommPtr()->barrier();
  }

  ++index_;

  if (outputManager_.getProcID() == 0)
    os << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : FrequencyCSV::doFinishOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void FrequencyCSV::doFinishOutput()
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
  }

  firstTimePrint_ = true;
}

//-----------------------------------------------------------------------------
// Function      : FrequencyCSV::doFinishOutputStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void FrequencyCSV::doFinishOutputStep()
{
  // Deal with the *prn file:
  if (outStreamPtr_)
  {
    if ( outputManager_.getPrintEndOfSimulationLine() )
    {
      (*outStreamPtr_) << "End of Xyce(TM) Parameter Sweep" << std::endl;
    }
  }

  outputManager_.closeFile(outStreamPtr_);
  outStreamPtr_ = 0;
}

//-----------------------------------------------------------------------------
// Class         : TimeTecPlot
// Purpose       : Outputter class for transient runs, tecplot output format
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : TimeTecPlot::TimeTecPlot
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
TimeTecPlot::TimeTecPlot(OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    firstTimePrint_(true),
    outFilename_(),
    suffix_(),
    outStreamPtr_(0),
    headerPrintCalls_(0)
{
  if (printParameters_.defaultExtension_.empty())
    printParameters_.defaultExtension_ = ".dat";

  fixupColumns(outputManager_, printParameters_, opList_);
}

//-----------------------------------------------------------------------------
// Function      : TimeTecPlot::~TimeTecPlot
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
TimeTecPlot::~TimeTecPlot()
{
  outputManager_.closeFile(outStreamPtr_);

  deleteList(opList_.begin(), opList_.end());
}

//-----------------------------------------------------------------------------
// Function      : TimeTecPlot::doParse
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void TimeTecPlot::doParse()
{
  outFilename_ = outputFilename(printParameters_, outputManager_.getNetListFilename());
}

//-----------------------------------------------------------------------------
// Function      : TimeTecPlot::doOpen
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void TimeTecPlot::doOpen()
{
  if (outputManager_.getProcID() == 0 && outStreamPtr_ == 0)
  {
    outStreamPtr_ = outputManager_.openFile(outFilename_);
  }
}

//-----------------------------------------------------------------------------
// Function      : TimeTecPlot::timeHeader
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
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
      for (Util::OpList::const_iterator it = opList_.begin() ; it != opList_.end(); ++it)
        os << "\" " << (*it)->getName() << "\" " << std::endl;

      // output some AUXDATA
      os << "DATASETAUXDATA " << getTecplotTimeDateStamp() << std::endl;

      if (!outputManager_.getTempSweepFlag())
      {
        os.setf(std::ios::scientific);
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
      os.setf(std::ios::scientific);
      os.precision( tecplotHeaderPrecision_);
      os << "T= \" ";
      for (std::vector<N_ANP_SweepParam>::const_iterator it = outputManager_.getStepParamVec().begin(); it != outputManager_.getStepParamVec().end(); ++it)
      {
        os << " " << it->name << " = " << it->currentVal;
      }
      os << "\" ";
    }
    os << std::endl;

    os.setf(std::ios::scientific);
    os.precision(printParameters_.streamPrecision_);

    // put in the various sweep parameters as auxdata:
    if (!outputManager_.getStepParamVec().empty())
    {
      for (std::vector<N_ANP_SweepParam>::const_iterator it = outputManager_.getStepParamVec().begin(); it != outputManager_.getStepParamVec().end(); ++it)
      {
        // convert any ":" or "%" in the name to a "_", so as not to confuse tecplot.
        std::string name(it->name);
        std::replace(name.begin(), name.end(), '%', '_');
        std::replace(name.begin(), name.end(), ':', '_');
        os << "AUXDATA " << name << " = " << "\" " << it->currentVal << "\" ";
      }
      os << std::endl;
    }

    os.setf(std::ios::left, std::ios::adjustfield);

  } // procID

  ++headerPrintCalls_;
}

//-----------------------------------------------------------------------------
// Function      : TimeTecPlot::doOutputTime
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void TimeTecPlot::doOutputTime(
    const N_LAS_Vector * solnVecPtr, 
    const N_LAS_Vector * stateVecPtr, 
    const N_LAS_Vector * storeVecPtr)
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
  for (Util::OpList::const_iterator it = opList_.begin() ; it != opList_.end(); ++it, ++i)
  {
    double result = getValue(outputManager_.getCommPtr()->comm(), *(*it), solnVecPtr, 0, stateVecPtr, storeVecPtr).real();
    result = filter(result, printParameters_.filter_);
    if ((*it)->opType() == Util::TIME_VAR)
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

//-----------------------------------------------------------------------------
// Function      : TimeTecPlot::doFinishOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void TimeTecPlot::doFinishOutput()
{
  if (outputManager_.getProcID() == 0)
  {
    if (outStreamPtr_)
    {
      if (!outputManager_.getSTEPEnabledFlag() && outputManager_.getPrintEndOfSimulationLine() )
      {
        (*outStreamPtr_) << "End of Xyce(TM) Simulation" << std::endl;
        outputManager_.closeFile(outStreamPtr_);
        outStreamPtr_ = 0;
        headerPrintCalls_ = 0;
      }
    }
  } // procID

  firstTimePrint_ = true;
}

//-----------------------------------------------------------------------------
// Function      : TimeTecPlot:: doFinishOutputStep()
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void TimeTecPlot::doFinishOutputStep()
{
  // Deal with the *tecplot file:
  if (outStreamPtr_)
  {
    if (outputManager_.getPrintEndOfSimulationLine())
    {
      (*outStreamPtr_) << "End of Xyce(TM) Parameter Sweep" << std::endl;
    }

    outputManager_.closeFile(outStreamPtr_);
    outStreamPtr_ = 0;
  }
}

//-----------------------------------------------------------------------------
// Class         : FrequencyTecPlot
// Purpose       : Outputter class for frequency domain runs, tecplot output
//                 format
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : FrequencyTecPlot::FrequencyTecPlot
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
FrequencyTecPlot::FrequencyTecPlot(OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_(),
    outStreamPtr_(0),
    stepCount_(0),
    firstTime_(true),
    index_(0)
{
  if (printParameters_.defaultExtension_.empty())
    printParameters_.defaultExtension_ = ".FD.dat";

  fixupColumns(outputManager_, printParameters_, opList_);
}

//-----------------------------------------------------------------------------
// Function      : FrequencyTecPlot::~FrequencyTecPlot()
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
FrequencyTecPlot::~FrequencyTecPlot()
{
  outputManager_.closeFile(outStreamPtr_);

  deleteList(opList_.begin(), opList_.end());
}


//-----------------------------------------------------------------------------
// Function      : FrequencyTecPlot::doParse()
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void FrequencyTecPlot::doParse()
{
  outFilename_ = outputFilename(printParameters_, outputManager_.getNetListFilename());
}

//-----------------------------------------------------------------------------
// Function      : FrequencyTecPlot::doOpen()
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void FrequencyTecPlot::doOpen()
{
  if (outputManager_.getProcID() == 0 && outStreamPtr_ == 0)
  {
    outStreamPtr_ = outputManager_.openFile(outFilename_);
  }
}

//-----------------------------------------------------------------------------
// Function      : FrequencyTecPlot::frequencyHeader
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void FrequencyTecPlot::frequencyHeader()
{
  index_ = 0;

  // STD header
  // Freq Domain Headers
  tecplotFreqHeader(outputManager_, printParameters_, opList_, *outStreamPtr_, stepCount_);

  outStreamPtr_->setf(std::ios::scientific);
  outStreamPtr_->precision(printParameters_.streamPrecision_);
  outStreamPtr_->setf(std::ios::left, std::ios::adjustfield);

  ++stepCount_;
}

//-----------------------------------------------------------------------------
// Function      : FrequencyTecPlot::doOutputFrequency
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void FrequencyTecPlot::doOutputFrequency(double frequency, const N_LAS_Vector *real_solution_vector, const N_LAS_Vector *imaginary_solution_vector)
{
  outputManager_.setCurrentOutputter(this);
  bool firstColPrinted = false;

  if (outputManager_.getProcID() == 0)
  {
    if (firstTime_) //Setup Output Stream and Print Out Header
    {
      doOpen();

      frequencyHeader();
      firstTime_ = false;
    }

  } // procID

  std::ostream &os = *outStreamPtr_;

  if (outputManager_.getProcID() == 0)
  {
    if (printParameters_.delimiter_ == "")
    {
      os.width(printParameters_.streamWidth_);
    }

    // Output a TAB.
    if (printParameters_.delimiter_ != "")
    {
      os << printParameters_.delimiter_;
    }

    if (printParameters_.delimiter_ == "")
    {
      os.width(printParameters_.streamWidth_);
    }
    else
    {
      os.width(0);
      if (firstColPrinted)
      {
        os << printParameters_.delimiter_;
      }
    }

    firstColPrinted = true;
  }

  // periodic time-domain steady-state output
  for (Util::OpList::const_iterator it = opList_.begin(); it != opList_.end(); ++it)
  {
    double result = getValue(outputManager_.getCommPtr()->comm(), *(*it), real_solution_vector, imaginary_solution_vector, 0, 0).real();
    if (outputManager_.getProcID() == 0)
    {
      os << result << " ";
    }
    outputManager_.getCommPtr()->barrier();

  } // end of output variable loop.

  if (outputManager_.getProcID() == 0)
  {
    os << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : FrequencyTecPlot::doFinishOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void FrequencyTecPlot::doFinishOutput()
{
  if (outputManager_.getProcID() == 0)
  {
    if (outStreamPtr_)
    {
      if (!outputManager_.getSTEPEnabledFlag() && outputManager_.getPrintEndOfSimulationLine())
      {
        (*outStreamPtr_) << "End of Xyce(TM) Simulation" << std::endl;
      }

      if (!outputManager_.getSTEPEnabledFlag())
      {
        outputManager_.closeFile(outStreamPtr_);
        outStreamPtr_ = 0;
      }
    }
  } // procID

  firstTime_ = true;
}

//-----------------------------------------------------------------------------
// Function      : FrequencyTecPlot::doFinishOutputStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void FrequencyTecPlot::doFinishOutputStep()
{
  if (outStreamPtr_)
  {
    if ( outputManager_.getPrintEndOfSimulationLine() )
    {
      (*outStreamPtr_) << "End of Xyce(TM) Parameter Sweep" << std::endl;
    }

    outputManager_.closeFile(outStreamPtr_);
    outStreamPtr_ = 0;
  }
}

//-----------------------------------------------------------------------------
// Class         : TimeProbe
// Purpose       : Outputter class for transient runs, Probe output
//                 format (PSpice-compatibility output)
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : TimeProbe::TimeProbe
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
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
  if (printParameters_.defaultExtension_.empty())
    printParameters_.defaultExtension_ = ".csd";

  fixupColumns(outputManager_, printParameters_, opList_);
}

//-----------------------------------------------------------------------------
// Function      : TimeProbe::~TimeProbe
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
TimeProbe::~TimeProbe()
{
  outputManager_.closeFile(outStreamPtr_);

  deleteList(opList_.begin(), opList_.end());
}


//-----------------------------------------------------------------------------
// Function      : TimeProbe::doParse
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void TimeProbe::doParse()
{
  outFilename_ = outputFilename(printParameters_, outputManager_.getNetListFilename());
}

//-----------------------------------------------------------------------------
// Function      : TimeProbe::doOpen
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void TimeProbe::doOpen()
{
  if (outputManager_.getProcID() == 0 && outStreamPtr_ == 0)
  {
    outStreamPtr_ = outputManager_.openFile(outFilename_);
  }
}

//-----------------------------------------------------------------------------
// Function      : TimeProbe::timeHeader
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void TimeProbe::timeHeader()
{
  std::ostream &os = *outStreamPtr_;

  if (outputManager_.getProcID() == 0)
  {
    // count the number of output variables.
    printCount_ = 0;
    for (Util::OpList::const_iterator it = opList_.begin() ; it != opList_.end(); ++it)
    {
      if ((*it)->opType() != Util::UNDEFINED)
        ++printCount_;
    }

    os << "#H" << std::endl;
    os << "SOURCE='Xyce' VERSION='"
       << N_UTL_Version::getShortVersionString() << "'" << std::endl;
    os << "TITLE='* " << outputManager_.getNetListFilename() << "'" << std::endl;

    os.setf(std::ios::scientific);
    os.precision(0); // streamPrecision_);
    if (outputManager_.getStepParamVec().empty())
    {
      os << "SUBTITLE='Xyce data";
    }
    else
    {
      os << "SUBTITLE='Step param";
      for (std::vector<N_ANP_SweepParam>::const_iterator it = outputManager_.getStepParamVec().begin(); it != outputManager_.getStepParamVec().end(); ++it)
      {
        os << " " << it->name << " = " << it->currentVal;
      }
    }

    os << " ' " << std::endl;;

    // set the time/date stamp
    os << getTimeDateStamp();
    os.setf(std::ios::scientific);
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
        os.setf(std::ios::scientific);
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
    for (Util::OpList::const_iterator it = opList_.begin() ; it != opList_.end(); ++it, ++i)
    {
      if (i > 18)
      {
        i = 0;
        os << std::endl;
      }
      os << "'" << (*it)->getName() << "' ";
    }
    if (i != 0)
      os << std::endl;

    os.setf(std::ios::scientific);
    os.precision(printParameters_.streamPrecision_);
    os.setf(std::ios::left, std::ios::adjustfield);

  } // procID

  ++headerPrintCalls_;
}

//-----------------------------------------------------------------------------
// Function      : TimeProbe::doOutputTime
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void TimeProbe::doOutputTime(
    const N_LAS_Vector * solnVecPtr, 
    const N_LAS_Vector * stateVecPtr, 
    const N_LAS_Vector * storeVecPtr)
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
  for (Util::OpList::const_iterator it = opList_.begin() ; it != opList_.end(); ++it, ++i)
  {
    double result = getValue(outputManager_.getCommPtr()->comm(), *(*it), solnVecPtr, 0, stateVecPtr, storeVecPtr).real();
    if ((*it)->opType() == Util::TIME_VAR)
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

//-----------------------------------------------------------------------------
// Function      : TimeProbe::doFinishOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
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

//-----------------------------------------------------------------------------
// Function      : TimeProbe::doFinishOutputStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void TimeProbe::doFinishOutputStep()
{
  outputManager_.closeFile(outStreamPtr_);
  outStreamPtr_ = 0;
}

//-----------------------------------------------------------------------------
// Class         : FrequencyProbe
// Purpose       : Outputter class for frequency-domain runs, Probe output
//                 format (PSpice-compatibility output)
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : FrequencyProbe::FrequencyProbe
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
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
  if (printParameters_.defaultExtension_.empty())
    printParameters_.defaultExtension_ = ".csd";

  fixupColumns(outputManager_, printParameters_, opList_);
}

//-----------------------------------------------------------------------------
// Function      : FrequencyProbe::~FrequencyProbe
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
FrequencyProbe::~FrequencyProbe()
{
  outputManager_.closeFile(outStreamPtr_);

  deleteList(opList_.begin(), opList_.end());
}


//-----------------------------------------------------------------------------
// Function      : FrequencyProbe::doParse
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void FrequencyProbe::doParse()
{
  outFilename_ = outputFilename(printParameters_, outputManager_.getNetListFilename());
}

//-----------------------------------------------------------------------------
// Function      : FrequencyProbe::doOpen
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void FrequencyProbe::doOpen()
{
  if (outputManager_.getProcID() == 0 && outStreamPtr_ == 0)
  {
    outStreamPtr_ = outputManager_.openFile(outFilename_);
  }
}

//-----------------------------------------------------------------------------
// Function      : FrequencyProbe::frequencyHeader
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void FrequencyProbe::frequencyHeader()
{
  std::ostream &os = *outStreamPtr_;

  if (outputManager_.getProcID() == 0)
  {
    // count the number of output variables.
    printCount_ = 0;
    for (Util::OpList::const_iterator it = opList_.begin() ; it != opList_.end(); ++it)
    {
      if ((*it)->opType() != Util::UNDEFINED)
        ++printCount_;
    }

    os << "#H" << std::endl;
    os << "SOURCE='Xyce' VERSION='"
       << N_UTL_Version::getShortVersionString() << "'" << std::endl;
    os << "TITLE='* " << outputManager_.getNetListFilename() << "'" << std::endl;

    os.setf(std::ios::scientific);
    os.precision(0); // streamPrecision_);
    if (outputManager_.getStepParamVec().empty())
    {
      os << "SUBTITLE='Xyce data";
    }
    else
    {
      os << "SUBTITLE='Step param";
      for (std::vector<N_ANP_SweepParam>::const_iterator it = outputManager_.getStepParamVec().begin(); it != outputManager_.getStepParamVec().end(); ++it)
      {
        os << " " << it->name << " = " << it->currentVal;
      }
    }

    os << " ' " << std::endl;;

    // set the time/date stamp
    os << getTimeDateStamp();
    os.setf(std::ios::scientific);
    os.precision(printParameters_.streamPrecision_);
    os << "TEMPERATURE='" << outputManager_.getCircuitTemp();
    os << "'" << std::endl;

    os << "ANALYSIS='AC Sweep' " <<
      "SERIALNO='12345'" <<  std::endl;

    os << "ALLVALUES='NO' COMPLEXVALUES='YES' " <<
      "NODES='" << printCount_ << "'" << std::endl;

    os << "SWEEPVAR='";
    std::string varName = outputManager_.getPRINTDCname();
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
        os.setf(std::ios::scientific);
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
    for (Util::OpList::const_iterator it = opList_.begin() ; it != opList_.end(); ++it, ++i)
    {
      os << "'" << (*it)->getName() << "' ";
      if (i > 3)
      {
        i = 0;
        os << std::endl;
      }
    }

    if (i != 0)
      os << std::endl;

    os.flush();

    os.setf(std::ios::scientific);
    os.precision(printParameters_.streamPrecision_);
    os.setf(std::ios::left, std::ios::adjustfield);

  } // procID

  ++headerPrintCalls_;
}

//-----------------------------------------------------------------------------
// Function      : FrequencyProbe::doOutputFrequency
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void FrequencyProbe::doOutputFrequency(double frequency, const N_LAS_Vector *real_solution_vector, const N_LAS_Vector *imaginary_solution_vector)
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
  for (Util::OpList::const_iterator it = opList_.begin() ; it != opList_.end(); ++it, ++i)
  {
    complex result = getValue(outputManager_.getCommPtr()->comm(), *(*it), real_solution_vector, imaginary_solution_vector, 0, 0);
    if (outputManager_.getProcID() == 0)
    {
      os << result.real() << "/" << result.imag() << ":" << i << "   ";
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

//-----------------------------------------------------------------------------
// Function      : FrequencyProbe::doFinishOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
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

//-----------------------------------------------------------------------------
// Function      : FrequencyProbe::doFinishOutputStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void FrequencyProbe::doFinishOutputStep()
{
  outputManager_.closeFile(outStreamPtr_);
  outStreamPtr_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : HBPrn::HBPrn
// Purpose       : Outputter for HB runs, standard (PRN) output
//                 format
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
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
  if (timePrintParameters_.defaultExtension_.empty())
    timePrintParameters_.defaultExtension_ = ".HB.TD.prn";

  if (freqPrintParameters_.defaultExtension_.empty())
    freqPrintParameters_.defaultExtension_ = ".HB.FD.prn";

  fixupColumns(outputManager_, timePrintParameters_, timeOpList_);
  fixupColumns(outputManager_, freqPrintParameters_, freqOpList_);
}

//-----------------------------------------------------------------------------
// Function      : HBPrn::~HBPrn
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
HBPrn::~HBPrn()
{
  outputManager_.closeFile(timeStreamPtr_);
  outputManager_.closeFile(freqStreamPtr_);

  deleteList(timeOpList_.begin(), timeOpList_.end());
  deleteList(freqOpList_.begin(), freqOpList_.end());
}


//-----------------------------------------------------------------------------
// Function      : HBPrn::doParse
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void HBPrn::doParse()
{
  timeFilename_ = outputFilename(timePrintParameters_, outputManager_.getNetListFilename());
  freqFilename_ = outputFilename(freqPrintParameters_, outputManager_.getNetListFilename());
}

//-----------------------------------------------------------------------------
// Function      : HBPrn::doOpen
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void HBPrn::doOpen()
{
  if (outputManager_.getProcID() == 0 && timeStreamPtr_ == 0)
  {
    timeStreamPtr_ = outputManager_.openFile(timeFilename_);
  }
  if (outputManager_.getProcID() == 0 && freqStreamPtr_ == 0)
  {
    freqStreamPtr_ = outputManager_.openFile(freqFilename_);
  }
}

//-----------------------------------------------------------------------------
// Function      : HBPrn::doOutputHB
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
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
      for (Util::OpList::const_iterator it = timeOpList_.begin(); it != timeOpList_.end(); ++it, ++column_index)
      {
        double result = getValue(outputManager_.getCommPtr()->comm(), *(*it), solnVecPtr, 0, 0, 0).real();
        if (outputManager_.getProcID() == 0)
          printValue(time_os, timePrintParameters_.table_.columnList_[column_index], timePrintParameters_.delimiter_, column_index, result);

        outputManager_.getCommPtr()->barrier();
      } // end of output variable loop.
    } // periodic time-domain steady-state output

    { // Fourier coefficient output
      int column_index = 0;
      for (Util::OpList::const_iterator it = freqOpList_.begin(); it != freqOpList_.end(); ++it, ++column_index)
      {
        // state and store vec are not available in this context, but we must
        // pass in both the real and imaginary vectors
        double result = getValue(outputManager_.getCommPtr()->comm(), *(*it), realVecPtr, imagVecPtr, 0, 0).real();
        if (outputManager_.getProcID() == 0)
          printValue(freq_os, freqPrintParameters_.table_.columnList_[column_index], freqPrintParameters_.delimiter_, column_index, result);

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

//-----------------------------------------------------------------------------
// Function      : HBPrn::doFinishOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
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

//-----------------------------------------------------------------------------
// Function      : HBPrn::doFinishOutputStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void HBPrn::doFinishOutputStep()
{
  // Deal with the *prn file:
  if (timeStreamPtr_)
  {
    if ( outputManager_.getPrintEndOfSimulationLine() )
    {
      (*timeStreamPtr_) << "End of Xyce(TM) Parameter Sweep" << std::endl;
    }
  }
  if (freqStreamPtr_)
  {
    if ( outputManager_.getPrintEndOfSimulationLine() )
    {
      (*freqStreamPtr_) << "End of Xyce(TM) Parameter Sweep" << std::endl;
    }
  }

  outputManager_.closeFile(timeStreamPtr_);
  timeStreamPtr_ = 0;
  outputManager_.closeFile(freqStreamPtr_);
  freqStreamPtr_ = 0;
}

//-----------------------------------------------------------------------------
// Class         : HBCSV
// Purpose       : Outputter class for HB runs, CSV (comma separated) output
//                 format
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : HBCSV::HBCSV
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
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
  if (timePrintParameters_.defaultExtension_.empty())
    timePrintParameters_.defaultExtension_ = ".HB.TD.csv";

  if (freqPrintParameters_.defaultExtension_.empty())
    freqPrintParameters_.defaultExtension_ = ".HB.FD.csv";

  fixupColumns(outputManager_, timePrintParameters_, timeOpList_);
  fixupColumns(outputManager_, freqPrintParameters_, freqOpList_);
}

//-----------------------------------------------------------------------------
// Function      : HBCSV::~HBCSV
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
HBCSV::~HBCSV()
{
  outputManager_.closeFile(timeStreamPtr_);
  outputManager_.closeFile(freqStreamPtr_);

  deleteList(timeOpList_.begin(), timeOpList_.end());
  deleteList(freqOpList_.begin(), freqOpList_.end());
}


//-----------------------------------------------------------------------------
// Function      : HBCSV::doParse
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void HBCSV::doParse()
{
  timeFilename_ = outputFilename(timePrintParameters_, outputManager_.getNetListFilename());
  freqFilename_ = outputFilename(freqPrintParameters_, outputManager_.getNetListFilename());
}

//-----------------------------------------------------------------------------
// Function      : HBCSV::doOpen
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void HBCSV::doOpen()
{
  if (outputManager_.getProcID() == 0 && timeStreamPtr_ == 0)
  {
    timeStreamPtr_ = outputManager_.openFile(timeFilename_);
  }
  if (outputManager_.getProcID() == 0 && freqStreamPtr_ == 0)
  {
    freqStreamPtr_ = outputManager_.openFile(freqFilename_);
  }
}

//-----------------------------------------------------------------------------
// Function      : HBCSV::doOutputHB
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
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
      for (Util::OpList::const_iterator it = timeOpList_.begin(); it != timeOpList_.end(); ++it, ++column_index)
      {
        double result = getValue(outputManager_.getCommPtr()->comm(), *(*it), solnVecPtr, 0, 0, 0).real();
        if (outputManager_.getProcID() == 0)
          printValue(time_os, timePrintParameters_.table_.columnList_[column_index], timePrintParameters_.delimiter_, column_index, result);

        outputManager_.getCommPtr()->barrier();
      } // end of output variable loop.
    } // periodic time-domain steady-state output

    { // Fourier coefficient output
      int column_index = 0;
      for (Util::OpList::const_iterator it = freqOpList_.begin(); it != freqOpList_.end(); ++it, ++column_index)
      {
        // state and store vec are not available in this context, but we must
        // pass in both the real and imaginary vectors
        double result = getValue(outputManager_.getCommPtr()->comm(), *(*it), realVecPtr, imagVecPtr, 0, 0).real();
        if (outputManager_.getProcID() == 0)
          printValue(freq_os, freqPrintParameters_.table_.columnList_[column_index], freqPrintParameters_.delimiter_, column_index, result);

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

//-----------------------------------------------------------------------------
// Function      : HBCSV::doFinishOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
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

//-----------------------------------------------------------------------------
// Function      : HBCSV::doFinishOutputStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void HBCSV::doFinishOutputStep()
{
  // Deal with the *prn file:
  if (timeStreamPtr_)
  {
    if ( outputManager_.getPrintEndOfSimulationLine() )
    {
      (*timeStreamPtr_) << "End of Xyce(TM) Parameter Sweep" << std::endl;
    }
  }
  if (freqStreamPtr_)
  {
    if ( outputManager_.getPrintEndOfSimulationLine() )
    {
      (*freqStreamPtr_) << "End of Xyce(TM) Parameter Sweep" << std::endl;
    }
  }

  outputManager_.closeFile(timeStreamPtr_);
  timeStreamPtr_ = 0;

  outputManager_.closeFile(freqStreamPtr_);
  freqStreamPtr_ = 0;
}

//-----------------------------------------------------------------------------
// Class         : HBTecPlot
// Purpose       : Outputter class for HB runs, TecPlot output
//                 format
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : HBTecPlot::HBTecPlot
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
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
  if (timePrintParameters_.defaultExtension_.empty())
    timePrintParameters_.defaultExtension_ = ".HB.TD.dat";

  if (freqPrintParameters_.defaultExtension_.empty())
    freqPrintParameters_.defaultExtension_ = ".HB.FD.dat";

  fixupColumns(outputManager_, timePrintParameters_, timeOpList_);
  fixupColumns(outputManager_, freqPrintParameters_, freqOpList_);
}

//-----------------------------------------------------------------------------
// Function      : HBTecPlot::~HBTecPlot
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
HBTecPlot::~HBTecPlot()
{
  outputManager_.closeFile(timeStreamPtr_);
  outputManager_.closeFile(freqStreamPtr_);

  deleteList(timeOpList_.begin(), timeOpList_.end());
  deleteList(freqOpList_.begin(), freqOpList_.end());
}

//-----------------------------------------------------------------------------
// Function      : HBTecPlot::doParse
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void HBTecPlot::doParse()
{
  timeFilename_ = outputFilename(timePrintParameters_, outputManager_.getNetListFilename());
  freqFilename_ = outputFilename(freqPrintParameters_, outputManager_.getNetListFilename());
}

//-----------------------------------------------------------------------------
// Function      : HBTecPlot::doOpen
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void HBTecPlot::doOpen()
{
  if (outputManager_.getProcID() == 0 && timeStreamPtr_ == 0)
  {
    timeStreamPtr_ = outputManager_.openFile(timeFilename_);
  }
  if (outputManager_.getProcID() == 0 && freqStreamPtr_ == 0)
  {
    freqStreamPtr_ = outputManager_.openFile(freqFilename_);
  }
}

//-----------------------------------------------------------------------------
// Function      : HBTecPlot::doOutputHB
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
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
  std::string pN(ost.str());

  if (outputManager_.getProcID() == 0)
  {

    if (firstTimeHB_) //Setup Output Stream and Print Out Header
    {
      doOpen();

      tecplotTimeHBHeader(*timeStreamPtr_);
      tecplotFreqHeader(outputManager_, freqPrintParameters_, freqOpList_, *freqStreamPtr_, stepCount_);

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
      for (Util::OpList::const_iterator it = timeOpList_.begin(); it != timeOpList_.end(); ++it, ++column_index)
      {
        double result = getValue(outputManager_.getCommPtr()->comm(), *(*it), solnVecPtr, 0, 0, 0).real();
        if (outputManager_.getProcID() == 0)
          printValue(time_os, timePrintParameters_.table_.columnList_[column_index], timePrintParameters_.delimiter_, column_index, result);

        outputManager_.getCommPtr()->barrier();
      } // end of output variable loop.
    } // periodic time-domain steady-state output

    { // Fourier coefficient output
      int column_index = 0;
      for (Util::OpList::const_iterator it = freqOpList_.begin(); it != freqOpList_.end(); ++it, ++column_index)
      {
        double result = getValue(outputManager_.getCommPtr()->comm(), *(*it), realVecPtr, imagVecPtr, 0, 0).real();
        if (outputManager_.getProcID() == 0)
          printValue(freq_os, freqPrintParameters_.table_.columnList_[column_index], freqPrintParameters_.delimiter_, column_index, result);

        outputManager_.getCommPtr()->barrier();
      }
    } // Fourier coefficient output

    if (outputManager_.getProcID() == 0)
    {
      freq_os << std::endl;
      time_os << std::endl;
    }

    ++index_;
  } // time scale loop.
}

//-----------------------------------------------------------------------------
// Function      : HBTecPlot::doFinishOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
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
  else
  {
    firstTimeHB_ = true;
  }
}

//-----------------------------------------------------------------------------
// Function      : HBTecPlot::doFinishOutputStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void HBTecPlot::doFinishOutputStep()
{
  // Deal with the *dat file:
  if (timeStreamPtr_)
  {
    if ( outputManager_.getPrintEndOfSimulationLine() )
    {
      (*timeStreamPtr_) << "End of Xyce(TM) Parameter Sweep" << std::endl;
    }
  }
  if (freqStreamPtr_)
  {
    if ( outputManager_.getPrintEndOfSimulationLine() )
    {
      (*freqStreamPtr_) << "End of Xyce(TM) Parameter Sweep" << std::endl;
    }
  }

  outputManager_.closeFile(timeStreamPtr_);
  timeStreamPtr_ = 0;
  outputManager_.closeFile(freqStreamPtr_);
  freqStreamPtr_ = 0;

}

//-----------------------------------------------------------------------------
// Function      : HBTecPlot::tecplotTimeHBHeader
// Purpose       : header for tecplot. Time Domain HB(default)
// Special Notes :
// Scope         :
// Creator       : Todd Coffey, Ting Mei
// Creation Date : 7/31/08
//-----------------------------------------------------------------------------
void HBTecPlot::tecplotTimeHBHeader( std::ostream & stream)
{
  std::ostream &os = *timeStreamPtr_;

  static const int tecplotHeaderPrecision_ = 2;
  os.setf(std::ios::scientific);
  os.precision( tecplotHeaderPrecision_);


  if (stepCount_ == 0)
  {
    os << " TITLE      = \" Xyce Time Domain HB data, " << outputManager_.getNetListFilename() << "\", " << std::endl;

    // output the user-specified solution vars:
    os << "\tVARIABLES = ";
    for (Util::OpList::const_iterator it = timeOpList_.begin() ; it != timeOpList_.end(); ++it)
    {
      os << "\" " << (*it)->getName() << "\" " << std::endl;
    }

    os << "DATASETAUXDATA " << getTecplotTimeDateStamp() << std::endl;
    if (!outputManager_.getTempSweepFlag())
    {
      stream << "DATASETAUXDATA TEMP = \"" << outputManager_.getCircuitTemp() <<
 " \"" << std::endl;
    }
  }

  os << "ZONE F=POINT  ";

  if (outputManager_.getStepParamVec().empty())
  {
    os << " T=\"Xyce data\" ";
  }
  else
  {
    os << " T= \" ";
    for (std::vector<N_ANP_SweepParam>::const_iterator it = outputManager_.getStepParamVec().begin(); it != outputManager_.getStepParamVec().end(); ++it)
    {
      os << " " << it->name << " = " << it->currentVal;
    }
    os << "\" ";
  }

  os << std::endl;

  // put in the various sweep parameters as auxdata:
  if (!outputManager_.getStepParamVec().empty())
  {
    for (std::vector<N_ANP_SweepParam>::const_iterator iterParam = outputManager_.getStepParamVec().begin();
        iterParam != outputManager_.getStepParamVec().end();
        ++iterParam)
    {
      // convert any ":" or "%" in the name to a "_", so as not to confuse tecplot.
      std::string tmpName(iterParam->name);
      replace(tmpName.begin(), tmpName.end(), '%', '_');
      replace(tmpName.begin(), tmpName.end(), ':', '_');
      os << "AUXDATA " << tmpName << " = " << "\" " << iterParam->currentVal
 << "\" ";
    }
    os << std::endl;
  }

  os << std::flush;
}

//-----------------------------------------------------------------------------
// Class         : MPDEPrn
// Purpose       : Outputter class for MPDE runs, standard (PRN) output
//                 format
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : MPDEPrn::MPDEPrn
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
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
  if (printParameters_.defaultExtension_.empty())
    printParameters_.defaultExtension_ = ".MPDE.prn";

  printParameters_.table_.addColumn("TIME1", printParameters_.streamWidth_, printParameters_.streamPrecision_, Table::JUSTIFICATION_RIGHT);
  printParameters_.table_.addColumn("TIME2", printParameters_.streamWidth_, printParameters_.streamPrecision_, Table::JUSTIFICATION_RIGHT);

  fixupColumns(outputManager_, printParameters_, opList_);
}

//-----------------------------------------------------------------------------
// Function      : MPDEPrn::~MPDEPrn
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
MPDEPrn::~MPDEPrn()
{
  outputManager_.closeFile(outStreamPtr_);

  deleteList(opList_.begin(), opList_.end());
}


//-----------------------------------------------------------------------------
// Function      : MPDEPrn::doParse
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void MPDEPrn::doParse()
{
  outFilename_ = outputFilename(printParameters_, outputManager_.getNetListFilename());
}

//-----------------------------------------------------------------------------
// Function      : MPDEPrn::doOpen
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void MPDEPrn::doOpen()
{
  if (outputManager_.getProcID() == 0 && outStreamPtr_ == 0)
  {
    outStreamPtr_ = outputManager_.openFile(outFilename_);
  }
}

//-----------------------------------------------------------------------------
// Function      : MPDEPrn::mpdeHeader
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void MPDEPrn::mpdeHeader()
{}

//-----------------------------------------------------------------------------
// Function      : MPDEPrn::doOutputMPDE
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
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
    for (Util::OpList::const_iterator it = opList_.begin(); it != opList_.end(); ++it, ++column_index)
    {
      double result = getValue(outputManager_.getCommPtr()->comm(), *(*it), solnVecPtr, 0, 0, 0).real();
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

//-----------------------------------------------------------------------------
// Function      : MPDEPrn::doFinishOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void MPDEPrn::doFinishOutput()
{
  outputManager_.closeFile(outStreamPtr_);
  outStreamPtr_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : MPDEPrn::doFinishOutputStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void MPDEPrn::doFinishOutputStep()
{
}

//-----------------------------------------------------------------------------
// Function      : MPDEPrn::stdTimeMPDEHeader
// Purpose       : header for std. Time Domain MPDE(default)
// Special Notes :
// Scope         :
// Creator       : Todd Coffey, Ting Mei
// Creation Date : 7/31/08
//-----------------------------------------------------------------------------
void MPDEPrn::stdTimeMPDEHeader( std::ostream & stream)
{
}

//-----------------------------------------------------------------------------
// Class         : MPDETecPlot
// Purpose       : Outputter class for MPDE runs, TecPlot output
//                 format
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : MPDETecPlot::MPDETecPlot
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
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

//-----------------------------------------------------------------------------
// Function      : MPDETecPlot::~MPDETecPlot
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
MPDETecPlot::~MPDETecPlot()
{
  outputManager_.closeFile(outStreamPtr_);

  deleteList(opList_.begin(), opList_.end());
}


//-----------------------------------------------------------------------------
// Function      : MPDETecPlot::doParse
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void MPDETecPlot::doParse()
{
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
}

//-----------------------------------------------------------------------------
// Function      : MPDETecPlot::doOpen
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void MPDETecPlot::doOpen()
{
  if (outputManager_.getProcID() == 0 && outStreamPtr_ == 0)
  {
    outStreamPtr_ = outputManager_.openFile(outFilename_);
  }
}

//-----------------------------------------------------------------------------
// Function      : MPDETecPlot::doOutputHeader
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void MPDETecPlot::doOutputHeader()
{
  (*outStreamPtr_) << " TITLE = \" Xyce MPDE data, " << outputManager_.getNetListFilename() << "\", " << std::endl
                   << "\tVARIABLES = \"T1(sec) \", \"T2(sec)\", " << std::endl;

  // output the user-specified solution vars:
  for (Util::OpList::const_iterator it = opList_.begin(); it != opList_.end(); ++it)
  {
    (*outStreamPtr_) << "\" "<< (*it)->getName() << "\" " << std::endl;
  }

  // output some AUXDATA
  (*outStreamPtr_) << "DATASETAUXDATA " << getTecplotTimeDateStamp() << std::endl
                   << "ZONE I=" << n2_ + 1 << ", " << " J=" << n1_ << ", " << " F=POINT\n" << std::endl;

  outStreamPtr_->setf(std::ios::scientific);
  outStreamPtr_->precision(printParameters_.streamPrecision_);
  outStreamPtr_->setf(std::ios::left, std::ios::adjustfield);
}

//-----------------------------------------------------------------------------
// Function      : MPDETecPlot::doOutputMPDE
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
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
  std::string pN(ost.str());

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
    {
      outStreamPtr_->width(printParameters_.streamWidth_);
    }

    // Output a TAB.
    if (printParameters_.delimiter_ != "")
    {
      (*outStreamPtr_) << printParameters_.delimiter_;
    }
  } // procID

  Util::OpList::const_iterator iterParam = opList_.begin();
  Util::OpList::const_iterator begin_tpL;
  Util::OpList::const_iterator last = opList_.end();

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
      {
        outStreamPtr_->width(printParameters_.streamWidth_);
      }
      else
      {
        outStreamPtr_->width(0);
        if (printParameters_.delimiter_ != "")
          (*outStreamPtr_) << printParameters_.delimiter_;
      }

      (*outStreamPtr_) << first;

      // time 2:
      if (printParameters_.delimiter_ == "")
      {
        outStreamPtr_->width(printParameters_.streamWidth_);
      }
      else
      {
        outStreamPtr_->width(0);
        if (printParameters_.delimiter_ != "")
        {
          (*outStreamPtr_) << printParameters_.delimiter_;
        }
      }

      (*outStreamPtr_) << second;
    }
//
    int i;
    for (i = 1, iterParam=begin_tpL ; iterParam != last; ++iterParam, ++i)
    {
      double result = getValue(outputManager_.getCommPtr()->comm(), *(*iterParam), solnVecPtr, 0, 0, 0).real();
      if (outputManager_.getProcID() == 0)
      {
        if (printParameters_.delimiter_ == "")
        {
          outStreamPtr_->width(printParameters_.streamWidth_);
        }
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
    (*outStreamPtr_) << std::endl;
    (*outStreamPtr_).flush();
  }
}

//-----------------------------------------------------------------------------
// Function      : MPDETecPlot::doFinishOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void MPDETecPlot::doFinishOutput()
{
  if (outStreamPtr_)
  {
    if (!outputManager_.getSTEPEnabledFlag() && outputManager_.getPrintEndOfSimulationLine())
    {
      (*outStreamPtr_) << "End of Xyce(TM) Simulation" << std::endl;
    }
  }

  if (!outputManager_.getSTEPEnabledFlag())
  {
    outputManager_.closeFile(outStreamPtr_);
    outStreamPtr_= 0;
  }
}

//-----------------------------------------------------------------------------
// Function      : MPDETecPlot::doFinishOutputStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void MPDETecPlot::doFinishOutputStep()
{
  // Deal with the *dat file:
  if (outStreamPtr_)
  {
    if ( outputManager_.getPrintEndOfSimulationLine() )
    {
      (*outStreamPtr_) << "End of Xyce(TM) Parameter Sweep" << std::endl;
    }
  }

  outputManager_.closeFile(outStreamPtr_);
  outStreamPtr_= 0;
}

//-----------------------------------------------------------------------------
// Function      : MPDETecPlot::stdTimeMPDEHeader
// Purpose       : header for std. Time Domain MPDE(default)
// Special Notes :
// Scope         : private
// Creator       : Todd Coffey, Ting Mei
// Creation Date : 7/31/08
//-----------------------------------------------------------------------------
void MPDETecPlot::stdTimeMPDEHeader( std::ostream & stream)
{
}

//-----------------------------------------------------------------------------
// Class         : HomotopyPrn
// Purpose       : Outputter class for homotopy output, standard (PRN) output
//                 format
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : HomotopyPrn::HomotopyPrn
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
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
  fixupColumns(outputManager_, printParameters_, opList_);
}

//-----------------------------------------------------------------------------
// Function      : HomotopyPrn::~HomotopyPrn
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
HomotopyPrn::~HomotopyPrn()
{
  outputManager_.closeFile(outStreamPtr_);

  deleteList(opList_.begin(), opList_.end());
}


//-----------------------------------------------------------------------------
// Function      : HomotopyPrn::doParse
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void HomotopyPrn::doParse()
{
  outFilename_ = outputFilename(printParameters_, outputManager_.getNetListFilename());
}

//-----------------------------------------------------------------------------
// Function      : HomotopyPrn::doOpen
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void HomotopyPrn::doOpen()
{
  if (outputManager_.getProcID() == 0 && outStreamPtr_ == 0)
  {
    outStreamPtr_ = outputManager_.openFile(outFilename_);
  }
}

//-----------------------------------------------------------------------------
// Function      : HomotopyPrn::homotopyHeader
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void HomotopyPrn::homotopyHeader(const std::vector<std::string> & parameter_names, 
    const std::vector<double> & param_values, const N_LAS_Vector * solution_vector)
{
  std::ostream &os = *outStreamPtr_;

  if (columnList_.empty())
  {
    Table::Justification justification = printParameters_.delimiter_.empty() ? 
      Table::JUSTIFICATION_CENTER :  Table::JUSTIFICATION_NONE;

    for (std::vector<std::string>::const_iterator it = parameter_names.begin(); 
        it != parameter_names.end(); ++it)
    {
      columnList_.push_back(Table::Column((*it), std::ios_base::scientific, 
            printParameters_.streamWidth_, printParameters_.streamPrecision_, justification));
    }
  }

  index_ = 0;

  int homotopyParamStartIndex=1;
  if (printParameters_.index_==false) // if noindex, then use 0 for start of homotopy params, otherwise 1
  {
    homotopyParamStartIndex=0;
  }

  if (stepCount_ == 0)
  {
    int column_index = 0;
    for (Table::ColumnList::const_iterator it = printParameters_.table_.columnList_.begin(); 
        it != printParameters_.table_.columnList_.end(); ++it, ++column_index)
    {
      if (it != printParameters_.table_.columnList_.begin())
      {
        os << (printParameters_.delimiter_.empty() ? " " : printParameters_.delimiter_);
      }

      if (column_index == homotopyParamStartIndex)
      {
        for (Table::ColumnList::const_iterator it2 = columnList_.begin(); it2 != columnList_.end(); ++it2)
        {
          if (it2 != columnList_.begin())
          {
            os << printParameters_.delimiter_;
          }
          printHeader(os, (*it2));
        }
      }

      printHeader(os, (*it));
    }

    os << std::endl;
  }

  ++stepCount_;
}

//-----------------------------------------------------------------------------
// Function      : HomotopyPrn::doOutputHomotopy
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
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

  int homotopyParamStartIndex=1;
  if (printParameters_.index_==false) // if noindex, then use 0 for start of homotopy params, otherwise 1
  {
    homotopyParamStartIndex=0;
  }

  int column_index = 0;
  for (Util::OpList::const_iterator it = opList_.begin(); it != opList_.end(); ++it, ++column_index)
  {
    double result = getValue(outputManager_.getCommPtr()->comm(), *(*it), solution_vector, 0, 0, 0).real();
    if (outputManager_.getProcID() == 0)
    {
      if (column_index == homotopyParamStartIndex)
      {
        for (int i = 0; i < parameter_values.size(); ++i)
        {
          printValue(os, columnList_[i], printParameters_.delimiter_, 1, parameter_values[i]);
        }
      }

      printValue(os, printParameters_.table_.columnList_[column_index], 
          printParameters_.delimiter_, column_index, result);
    }

    outputManager_.getCommPtr()->barrier();
  }

  if (outputManager_.getProcID() == 0)
  {
    os << std::endl;
  }

  ++index_;
}

//-----------------------------------------------------------------------------
// Function      : HomotopyPrn::doFinishOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void HomotopyPrn::doFinishOutput()
{
  if (outputManager_.getProcID() == 0)
  {
    if (outStreamPtr_)
    {
      if (!outputManager_.getSTEPEnabledFlag() && outputManager_.getPrintEndOfSimulationLine())
      {
        (*outStreamPtr_) << "End of Xyce(TM) Homotopy Simulation" << std::endl;
      }
    }

    if (!outputManager_.getSTEPEnabledFlag())
    {
      outputManager_.closeFile(outStreamPtr_);
      outStreamPtr_ = 0;
    }
  } // procID
  firstTimeHomotopy_ = true;
}

//-----------------------------------------------------------------------------
// Function      : HomotopyPrn::doFinishOutputStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void HomotopyPrn::doFinishOutputStep()
{
  // close the homotopy file.
  if (outStreamPtr_)
  {
    if ( outputManager_.getPrintEndOfSimulationLine() )
    {
      (*outStreamPtr_) << "End of Xyce(TM) Homotopy Simulation" << std::endl;
    }
  }

  outputManager_.closeFile(outStreamPtr_);
  outStreamPtr_ = 0;
}

//-----------------------------------------------------------------------------
// Class         : HomotopyTecPlot
// Purpose       : Outputter class for homotopy output, TecPlot output
//                 format
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : HomotopyTecPlot::HomotopyTecPlot
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
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
  fixupColumns(outputManager_, printParameters_, opList_);
}

//-----------------------------------------------------------------------------
// Function      : HomotopyTecPlot::~HomotopyTecPlot
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
HomotopyTecPlot::~HomotopyTecPlot()
{
  outputManager_.closeFile(outStreamPtr_);
  deleteList(opList_.begin(), opList_.end());
}

//-----------------------------------------------------------------------------
// Function      : HomotopyTecPlot::doParse
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void HomotopyTecPlot::doParse()
{
  outFilename_ = outputFilename(printParameters_, outputManager_.getNetListFilename());
}

//-----------------------------------------------------------------------------
// Function      : HomotopyTecPlot::doOpen
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void HomotopyTecPlot::doOpen()
{
  if (outputManager_.getProcID() == 0 && outStreamPtr_ == 0)
  {
    outStreamPtr_ = outputManager_.openFile(outFilename_);
  }
}

//-----------------------------------------------------------------------------
// Function      : HomotopyTecPlot::doOutputHeader
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void HomotopyTecPlot::doOutputHeader(const std::vector<std::string> & parameter_names, 
    const std::vector<double> & param_values, 
    const N_LAS_Vector * solution_vector)
{
  if (columnList_.empty())
  {
    Table::Justification justification = printParameters_.delimiter_.empty() ? 
      Table::JUSTIFICATION_CENTER :  Table::JUSTIFICATION_NONE;

    for (std::vector<std::string>::const_iterator it = parameter_names.begin(); 
        it != parameter_names.end(); ++it)
    {
      columnList_.push_back(Table::Column((*it), std::ios_base::scientific, 
        printParameters_.streamWidth_, printParameters_.streamPrecision_, justification));
    }
  }

  std::ostream &os = *outStreamPtr_;

  index_ = 0;
  ParameterList::const_iterator iterParam = printParameters_.variableList_.begin();
  ParameterList::const_iterator last = printParameters_.variableList_.end();

  if (stepCount_ == 0)
  {
    os << " TITLE = \" Xyce homotopy data, " << outputManager_.getNetListFilename() << "\", " << std::endl;
    os << "\tVARIABLES = ";

    // output the continuation parameters:
    std::vector<std::string>::const_iterator iter_name;
    for (iter_name = parameter_names.begin(); iter_name!= parameter_names.end(); ++iter_name)
    {
      os << "\" ";
      os << *iter_name;
      os << "\" " << std::endl;
    }

    // output the user-specified solution vars:
    for (Util::OpList::const_iterator it = opList_.begin(); it != opList_.end(); ++it)
    {
      os << "\" " << (*it)->getName() << "\" " << std::endl;
    }
  }

  // output some AUXDATA
  os << "DATASETAUXDATA ";
  os << getTecplotTimeDateStamp();
  os << std::endl;

  os << "ZONE F=POINT";


  if (outputManager_.getStepParamVec().empty())
  {
    os << " T=\"Xyce data\" ";
  }
  else
  {
    os << " T= \" ";
    for (std::vector<N_ANP_SweepParam>::const_iterator it = outputManager_.getStepParamVec().begin(); 
         it != outputManager_.getStepParamVec().end(); ++it)
    {
      static const int tecplotHeaderPrecision = 2;
      os.setf(std::ios::scientific);
      os.precision(tecplotHeaderPrecision);
      os << " " << it->name << " = " << it->currentVal;
    }
    os << "\" ";
  }

  os << std::endl;

  os.setf(std::ios::scientific);
  os.precision(printParameters_.streamPrecision_);
  os.setf(std::ios::left, std::ios::adjustfield);

  ++stepCount_;
}

//-----------------------------------------------------------------------------
// Function      : HomotopyTecPlot::doOutputHomotopy
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void HomotopyTecPlot::doOutputHomotopy(const std::vector<std::string> & parameter_names, 
    const std::vector<double> & parameter_values, const N_LAS_Vector * solution_vector)
{
  outputManager_.setCurrentOutputter(this);


  double tmpTime = outputManager_.getAnaIntPtr()->getTime();

  if (outputManager_.getProcID() == 0) 
  {
    if (firstTimeHomotopy_) //Setup Output Stream and Print Out Header
    {
      doOpen();
      doOutputHeader(parameter_names, parameter_values, solution_vector);
      firstTimeHomotopy_ = false;
    }
  }

  std::ostream &os = *outStreamPtr_;

  int column_index = 0;
  for (Util::OpList::const_iterator it = opList_.begin(); 
      it != opList_.end(); ++it, ++column_index)
  {
    double result = 
      getValue(outputManager_.getCommPtr()->comm(), *(*it), solution_vector, 0, 0, 0).real();
    if (outputManager_.getProcID() == 0)
    {
      if (column_index == 0)
      {
        for (int i = 0; i < parameter_values.size(); ++i)
        {
          printValue(os, columnList_[i], printParameters_.delimiter_, 1, parameter_values[i]);
        }
      }

      printValue(os, printParameters_.table_.columnList_[column_index], 
          printParameters_.delimiter_, column_index, result);
    }

    outputManager_.getCommPtr()->barrier();
  }

  if (outputManager_.getProcID() == 0)
  {
    os << std::endl;
  }

  ++index_;
}

//-----------------------------------------------------------------------------
// Function      : HomotopyTecPlot::doFinishOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void HomotopyTecPlot::doFinishOutput()
{
  if (outputManager_.getProcID() == 0)
  {
    if (outStreamPtr_)
    {
      if (!outputManager_.getSTEPEnabledFlag() && outputManager_.getPrintEndOfSimulationLine())
      {
        (*outStreamPtr_) << "End of Xyce(TM) Homotopy Simulation" << std::endl;
      }
    }

    if (!outputManager_.getSTEPEnabledFlag())
    {
      outputManager_.closeFile(outStreamPtr_);
      outStreamPtr_ = 0;
    }
  } // procID
  firstTimeHomotopy_ = true;
}

//-----------------------------------------------------------------------------
// Function      : HomotopyTecPlot::doFinishOutputStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void HomotopyTecPlot::doFinishOutputStep()
{
  // close the homotopy file.
  if (outStreamPtr_)
  {
    if ( outputManager_.getPrintEndOfSimulationLine() )
    {
      (*outStreamPtr_) << "End of Xyce(TM) Homotopy Simulation" << std::endl;
    }
  }

  outputManager_.closeFile(outStreamPtr_);
  outStreamPtr_ = 0;
}

//-----------------------------------------------------------------------------
// Class         : HomotopyProbe
// Purpose       : Outputter class for homotopy output, Probe output
//                 format
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : HomotopyProbe::HomotopyProbe
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
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

//-----------------------------------------------------------------------------
// Function      : HomotopyProbe::~HomotopyProbe
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
HomotopyProbe::~HomotopyProbe()
{
  outputManager_.closeFile(outStreamPtr_);

  deleteList(opList_.begin(), opList_.end());
}

//-----------------------------------------------------------------------------
// Function      :
// Purpose       : HomotopyProbe::doParse
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
HomotopyProbe::doParse()
{
  outFilename_ = outputFilename(printParameters_, outputManager_.getNetListFilename());
}

//-----------------------------------------------------------------------------
// Function      : HomotopyProbe::doOpen
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void  HomotopyProbe::doOpen()
{
  if (outputManager_.getProcID() == 0 && outStreamPtr_ == 0)
  {
    outStreamPtr_ = outputManager_.openFile(outFilename_);
  }
}

//-----------------------------------------------------------------------------
// Function      : HomotopyProbe::doOutputHeader
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void HomotopyProbe::doOutputHeader(const std::vector<std::string> & parameter_names, const std::vector<double> & param_values, const N_LAS_Vector * solution_vector)
{
  std::ostream &os = *outStreamPtr_;

  index_ = 0;

  printCount_ = opList_.size();

  os << "#H" << std::endl;

  os << "SOURCE='Xyce' VERSION='"
                   << N_UTL_Version::getShortVersionString() << "'" << std::endl;

  os << "TITLE='* " << outputManager_.getNetListFilename() << "'" << std::endl;
  os << "SUBTITLE='spice probe data'" << std::endl;

  // set the time/date stamp
  os << getTimeDateStamp();
  os.setf(std::ios::scientific);
  os.precision(printParameters_.streamPrecision_);
  os << "TEMPERATURE='" << outputManager_.getCircuitTemp();
  os << "'" << std::endl;

  if (printParameters_.printType_ == PrintType::TRAN)
    os << "ANALYSIS='Transient Analysis' SERIALNO='12345'" <<  std::endl;
  else
    os << "ANALYSIS='DC transfer characteristic' " << "SERIALNO='12345'" <<  std::endl;

  os << "ALLVALUES='NO' COMPLEXVALUES='NO' " << "NODES='" << printCount_ << "'" << std::endl;

  if (printParameters_.printType_ == PrintType::TRAN)
  {
    os << "SWEEPVAR='Time' SWEEPMODE='VAR_STEP'" << std::endl;
  }
  else
  {
    os << "SWEEPVAR='Voltage' SWEEPMODE='VAR_STEP'" <<std::endl;
  }

  // This line assumes that we're doing a homotopy that goes from 0 to 1.
  // This will never be a transient output.
  os << "XBEGIN='0.0'  XEND='1.0'"<<std::endl;

  os << "FORMAT='0 VOLTSorAMPS;EFLOAT : " << "NODEorBRANCH;NODE  '  " << std::endl;

  os << "DGTLDATA='NO'" << std::endl;

  os << "#N" << std::endl;

  // print the continuation parameter names:
  int i = 0;
  for (std::vector<std::string>::const_iterator iter_name = parameter_names.begin(); iter_name != parameter_names.end(); ++iter_name, ++i)
  {
    os << "'" << *iter_name << "' ";

    if (i > 3)
    {
      i = 0;
      os << std::endl;
    }
  }

  // print output variable names:
  for (Util::OpList::const_iterator it = opList_.begin(); it != opList_.end(); ++it, ++i)
  {
    os << "'" << (*it)->getName() << "' ";

    if (i > 3)
    {
      i = 0;
      os << std::endl;
    }
  }
  if (i != 0)
    os << std::endl;


  os.flush();

  os.setf(std::ios::scientific);
  os.precision(printParameters_.streamPrecision_);
  os.setf(std::ios::left, std::ios::adjustfield);


  ++stepCount_;
}

//-----------------------------------------------------------------------------
// Function      : HomotopyProbe::doOutputHomotopy
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void  HomotopyProbe::doOutputHomotopy(const std::vector<std::string> & parameter_names, const std::vector<double> & parameter_values, const N_LAS_Vector * solution_vector)
{
  std::ostream &os = *outStreamPtr_;

  outputManager_.setCurrentOutputter(this);

  double tmpTime = outputManager_.getAnaIntPtr()->getTime();

  if (outputManager_.getProcID() == 0)
  {

    if (firstTimeHomotopy_) //Setup Output Stream and Print Out Header
    {
      doOpen();
      doOutputHeader(parameter_names, parameter_values, solution_vector);

      firstTimeHomotopy_ = false;
    }

    os.width( 0);
    if (printParameters_.printType_ == PrintType::TRAN)
    {
      os << "#C " << tmpTime << " ";
      os << printCount_ << std::endl;
    }
    else
    {
      os << "#C " << outputManager_.getPRINTDCvalue() << " ";
      os << printCount_ << std::endl;
    }

    //-------------------------------------
    //HOMOTOPY PARAM VALUE OUTPUT GOES HERE
    //-------------------------------------

    for (int iparam=0;iparam < parameter_values.size(); ++iparam)
    {

      if (printParameters_.delimiter_ == "")
      {
        os.width(printParameters_.streamWidth_);
      }
      else
      {
        os.width(0);
        if (printParameters_.delimiter_ != "")
          os << printParameters_.delimiter_;
      }

      os << parameter_values[iparam];
    }
  } // procID

  Util::OpList::const_iterator iterParam = opList_.begin();
  Util::OpList::const_iterator last = opList_.end();

  int i;
  for (i = 1; iterParam != last; ++iterParam, ++i)
  {
    double result = getValue(outputManager_.getCommPtr()->comm(), *(*iterParam), solution_vector, 0, 0, 0).real();
    if (outputManager_.getProcID() == 0)
    {
      if (printParameters_.delimiter_ == "")
      {
        os.width(printParameters_.streamWidth_);
      }
      else
      {
        os.width(0);
        if (printParameters_.delimiter_ != "")
          os << printParameters_.delimiter_;
      }
      os << result;
    }

    outputManager_.getCommPtr()->barrier();
  }

  if (outputManager_.getProcID() == 0)
    os << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : HomotopyProbe::doFinishOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void  HomotopyProbe::doFinishOutput()
{
  firstTimeHomotopy_ = true;
}

//-----------------------------------------------------------------------------
// Function      : HomotopyProbe::doFinishOutputStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
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



//-----------------------------------------------------------------------------
// Class         : SensitivityPrn
// Purpose       : Outputter class for sensitivity output, standard (PRN) output
//                 format
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : SensitivityPrn::SensitivityPrn
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
SensitivityPrn::SensitivityPrn(OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_(),
    outStreamPtr_(0),
    stepCount_(0),
    index_(0),
    printCount_(0),
    firstTimeSensitivity_(true),
    headerPrintCalls_(0)
{
  if (printParameters_.defaultExtension_.empty())
    printParameters_.defaultExtension_ = "SENS.prn";

  fixupColumns(outputManager_, printParameters_, opList_);
}

//-----------------------------------------------------------------------------
// Function      : SensitivityPrn::~SensitivityPrn
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
SensitivityPrn::~SensitivityPrn()
{
  outputManager_.closeFile(outStreamPtr_);

  deleteList(opList_.begin(), opList_.end());
}


//-----------------------------------------------------------------------------
// Function      : SensitivityPrn::doParse
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
void SensitivityPrn::doParse()
{
  outFilename_ = outputFilename(printParameters_, outputManager_.getNetListFilename());
}

//-----------------------------------------------------------------------------
// Function      : SensitivityPrn::doOpen
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
void SensitivityPrn::doOpen()
{
  if (outputManager_.getProcID() == 0 && outStreamPtr_ == 0)
  {
    outStreamPtr_ = outputManager_.openFile(outFilename_);
  }
}

//-----------------------------------------------------------------------------
// Function      : SensitivityPrn::sensitivityHeader
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
void SensitivityPrn::sensitivityHeader(
    const std::vector<std::string> & parameter_names)
{
  std::ostream &os = *outStreamPtr_;

  if (columnList_.empty())
  {
    Table::Justification justification = 
      printParameters_.delimiter_.empty() ? Table::JUSTIFICATION_CENTER :  Table::JUSTIFICATION_NONE;

    for (std::vector<std::string>::const_iterator 
        it = parameter_names.begin(); it != parameter_names.end(); ++it)
    {
      columnList_.push_back(Table::Column((*it), std::ios_base::scientific, 
            printParameters_.streamWidth_, 
            printParameters_.streamPrecision_, 
            justification));
    }
  }

  index_ = 0;

  if (stepCount_ == 0)
  {
    int column_index = 0;
    for (Table::ColumnList::const_iterator 
        it = printParameters_.table_.columnList_.begin(); 
        it != printParameters_.table_.columnList_.end(); 
        ++it, ++column_index)
    {
      if (it != printParameters_.table_.columnList_.begin())
      {
        os << (printParameters_.delimiter_.empty() ? " " : printParameters_.delimiter_);
      }
      printHeader(os, (*it));
    }

    for (Table::ColumnList::const_iterator it2 = columnList_.begin(); it2 != columnList_.end(); ++it2)
    {
      if (it2 != columnList_.begin())
      {
        os << printParameters_.delimiter_;
      }
      printHeader(os, (*it2));
    }
    os << std::endl;
  }

  ++stepCount_;
}

//-----------------------------------------------------------------------------
// Function      : SensitivityPrn::doOutputSensitivity
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
void SensitivityPrn::doOutputSensitivity(
    const std::vector<std::string> & parameter_names, 
    const std::vector<double> & objective_values, 
    const std::vector<double> & direct_values, 
    const std::vector<double> & adjoint_values,
    const std::vector<double> & scaled_direct_values, 
    const std::vector<double> & scaled_adjoint_values,
    const N_LAS_Vector *solution_vector,
    const N_LAS_Vector *state_vector, 
    const N_LAS_Vector *store_vector)
{
  outputManager_.setCurrentOutputter(this);

  double tmpTime = outputManager_.getAnaIntPtr()->getTime();

  if (outputManager_.getProcID() == 0)
  {

    if (firstTimeSensitivity_) //Setup Output Stream and Print Out Header
    {
      doOpen();

      sensitivityHeader(parameter_names);

      firstTimeSensitivity_ = false;

      ++headerPrintCalls_;

      index_ = 0;
    }
  }

  std::ostream &os = *outStreamPtr_;
  outputManager_.getCommPtr()->barrier();

  int column_index = 0;
  for (Util::OpList::const_iterator it = opList_.begin(); it != opList_.end(); ++it, ++column_index)
  {
     double result = getValue(outputManager_.getCommPtr()->comm(), *(*it), 
         solution_vector, 0, state_vector, store_vector).real();

    result = filter(result, printParameters_.filter_);

    if ((*it)->opType() == Util::TIME_VAR)
      result *= printParameters_.outputTimeScaleFactor_;

    if (outputManager_.getProcID() == 0)
      printValue(os, printParameters_.table_.columnList_[column_index], 
          printParameters_.delimiter_, column_index, result);

    outputManager_.getCommPtr()->barrier();
  }

  for (int i = 0; i < objective_values.size(); ++i)
  {
    if (outputManager_.getProcID() == 0)
    {
      printValue(os, columnList_[i], printParameters_.delimiter_, 1, objective_values[i]);
    }
    outputManager_.getCommPtr()->barrier();
  }

  for (int i = 0; i < direct_values.size(); ++i)
  {
    if (outputManager_.getProcID() == 0)
    {
      printValue(os, columnList_[i], printParameters_.delimiter_, 1, direct_values[i]);
    }
    outputManager_.getCommPtr()->barrier();
  }

  for (int i = 0; i < scaled_direct_values.size(); ++i)
  {
    if (outputManager_.getProcID() == 0)
    {
      printValue(os, columnList_[i], printParameters_.delimiter_, 1, scaled_direct_values[i]);
    }
    outputManager_.getCommPtr()->barrier();
  }

  for (int i = 0; i < adjoint_values.size(); ++i)
  {
    if (outputManager_.getProcID() == 0)
    {
      printValue(os, columnList_[i], printParameters_.delimiter_, 1, adjoint_values[i]);
    }
    outputManager_.getCommPtr()->barrier();
  }

  for (int i = 0; i < scaled_adjoint_values.size(); ++i)
  {
    if (outputManager_.getProcID() == 0)
    {
      printValue(os, columnList_[i], printParameters_.delimiter_, 1, scaled_adjoint_values[i]);
    }
    outputManager_.getCommPtr()->barrier();
  }

  if (outputManager_.getProcID() == 0)
  {
    os << std::endl;
  }

  ++index_;
}

//-----------------------------------------------------------------------------
// Function      : SensitivityPrn::doFinishOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
void SensitivityPrn::doFinishOutput()
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

  firstTimeSensitivity_ = true;
  index_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : SensitivityPrn::doFinishOutputStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
void SensitivityPrn::doFinishOutputStep()
{
  // close the sensitivity file.
  if (outStreamPtr_)
  {
    if ( outputManager_.getPrintEndOfSimulationLine() )
    {
      (*outStreamPtr_) << "End of Xyce(TM) Sensitivity Simulation" << std::endl;
    }
  }

  outputManager_.closeFile(outStreamPtr_);
  outStreamPtr_ = 0;
}


//-----------------------------------------------------------------------------
// Class         : SensitivityTecPlot
// Purpose       : Outputter class for sensitivity output, TecPlot output
//                 format
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : SensitivityTecPlot::SensitivityTecPlot
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
SensitivityTecPlot::SensitivityTecPlot(OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_(),
    outStreamPtr_(0),
    stepCount_(0),
    index_(0),
    printCount_(0),
    firstTimeSensitivity_(true),
    headerPrintCalls_(0)
{
  if (printParameters_.defaultExtension_.empty())
    printParameters_.defaultExtension_ = "SENS.dat";

  fixupColumns(outputManager_, printParameters_, opList_);
}

//-----------------------------------------------------------------------------
// Function      : SensitivityTecPlot::~SensitivityTecPlot
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
SensitivityTecPlot::~SensitivityTecPlot()
{
  outputManager_.closeFile(outStreamPtr_);

  deleteList(opList_.begin(), opList_.end());
}

//-----------------------------------------------------------------------------
// Function      : SensitivityTecPlot::doParse
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
void SensitivityTecPlot::doParse()
{
  outFilename_ = outputFilename(printParameters_, outputManager_.getNetListFilename());
}

//-----------------------------------------------------------------------------
// Function      : SensitivityTecPlot::doOpen
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
void SensitivityTecPlot::doOpen()
{
  if (outputManager_.getProcID() == 0 && outStreamPtr_ == 0)
  {
    outStreamPtr_ = outputManager_.openFile(outFilename_);
  }
}

//-----------------------------------------------------------------------------
// Function      : SensitivityTecPlot::sensitivityHeader
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
void SensitivityTecPlot::sensitivityHeader(
    const std::vector<std::string> & parameter_names)
{
  std::ostream &os = *outStreamPtr_;

  if (columnList_.empty())
  {
    Table::Justification justification = 
      printParameters_.delimiter_.empty() ? Table::JUSTIFICATION_CENTER :  Table::JUSTIFICATION_NONE;

    for (std::vector<std::string>::const_iterator 
        it = parameter_names.begin(); it != parameter_names.end(); ++it)
    {
      columnList_.push_back(Table::Column((*it), std::ios_base::scientific, 
            printParameters_.streamWidth_, 
            printParameters_.streamPrecision_, 
            justification));
    }
  }


  index_ = 0;
  if (outputManager_.getProcID() == 0)
  {
    int tecplotHeaderPrecision_ = 2;

    if (stepCount_ == 0)
    {
      os << "TITLE = \"" << outputManager_.getNetListFilename() << " - " << outputManager_.getTitle() << "\", " << std::endl;
      os << "\tVARIABLES = ";

      // output the user-specified solution vars:
      for (Util::OpList::const_iterator 
          it = opList_.begin() ; it != opList_.end(); ++it)
      {
        os << "\" " << (*it)->getName() << "\" " << std::endl;
      }

      // output the .sens params
      std::vector<std::string>::const_iterator iter_name;
      for (iter_name = parameter_names.begin(); iter_name!= parameter_names.end(); ++iter_name)
      {
        os << "\" ";
        os << *iter_name;
        os << "\" " << std::endl;
      }  

      // output some AUXDATA
      os << "DATASETAUXDATA " << getTecplotTimeDateStamp() << std::endl;

      if (!outputManager_.getTempSweepFlag())
      {
        os.setf(std::ios::scientific);
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
      os.setf(std::ios::scientific);
      os.precision( tecplotHeaderPrecision_);
      os << "T= \" ";
      for (std::vector<N_ANP_SweepParam>::const_iterator 
          it = outputManager_.getStepParamVec().begin(); 
          it != outputManager_.getStepParamVec().end(); ++it)
      {
        os << " " << it->name << " = " << it->currentVal;
      }
      os << "\" ";
    }
    os << std::endl;

    os.setf(std::ios::scientific);
    os.precision(printParameters_.streamPrecision_);

    // put in the various sweep parameters as auxdata:
    if (!outputManager_.getStepParamVec().empty())
    {
      for (std::vector<N_ANP_SweepParam>::const_iterator 
          it = outputManager_.getStepParamVec().begin(); 
          it != outputManager_.getStepParamVec().end(); ++it)
      {
        // convert any ":" or "%" in the name to a "_", so as not to confuse tecplot.
        std::string name(it->name);
        std::replace(name.begin(), name.end(), '%', '_');
        std::replace(name.begin(), name.end(), ':', '_');
        os << "AUXDATA " << name << " = " << "\" " << it->currentVal << "\" ";
      }
      os << std::endl;
    }

    os.setf(std::ios::left, std::ios::adjustfield);

  } // procID

  ++stepCount_;
}

//-----------------------------------------------------------------------------
// Function      : SensitivityTecPlot::doOutputSensitivity
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
void SensitivityTecPlot::doOutputSensitivity(
    const std::vector<std::string> & parameter_names, 
    const std::vector<double> & objective_values, 
    const std::vector<double> & direct_values, 
    const std::vector<double> & adjoint_values,
    const std::vector<double> & scaled_direct_values, 
    const std::vector<double> & scaled_adjoint_values,
    const N_LAS_Vector *solution_vector,
    const N_LAS_Vector *state_vector, 
    const N_LAS_Vector *store_vector)
{
  outputManager_.setCurrentOutputter(this);
  double tmpTime = outputManager_.getAnaIntPtr()->getTime();

  if (outputManager_.getProcID() == 0)
  {

    if (firstTimeSensitivity_) //Setup Output Stream and Print Out Header
    {
      doOpen();

      sensitivityHeader(parameter_names);

      firstTimeSensitivity_ = false;

      ++headerPrintCalls_;

      index_ = 0;
    }
  }

  std::ostream &os = *outStreamPtr_;
  outputManager_.getCommPtr()->barrier();

  int column_index = 0;
  for (Util::OpList::const_iterator it = opList_.begin(); it != opList_.end(); ++it, ++column_index)
  {
     double result = getValue(outputManager_.getCommPtr()->comm(), *(*it), 
         solution_vector, 0, state_vector, store_vector).real();

    result = filter(result, printParameters_.filter_);

    if ((*it)->opType() == Util::TIME_VAR)
      result *= printParameters_.outputTimeScaleFactor_;

    if (outputManager_.getProcID() == 0)
      printValue(os, printParameters_.table_.columnList_[column_index], 
          printParameters_.delimiter_, column_index, result);

    outputManager_.getCommPtr()->barrier();
  }

  for (int i = 0; i < objective_values.size(); ++i)
  {
    if (outputManager_.getProcID() == 0)
    {
      printValue(os, columnList_[i], printParameters_.delimiter_, 1, objective_values[i]);
    }
    outputManager_.getCommPtr()->barrier();
  }

  for (int i = 0; i < direct_values.size(); ++i)
  {
    if (outputManager_.getProcID() == 0)
    {
      printValue(os, columnList_[i], printParameters_.delimiter_, 1, direct_values[i]);
    }
    outputManager_.getCommPtr()->barrier();
  }

  for (int i = 0; i < scaled_direct_values.size(); ++i)
  {
    if (outputManager_.getProcID() == 0)
    {
      printValue(os, columnList_[i], printParameters_.delimiter_, 1, scaled_direct_values[i]);
    }
    outputManager_.getCommPtr()->barrier();
  }

  for (int i = 0; i < adjoint_values.size(); ++i)
  {
    if (outputManager_.getProcID() == 0)
    {
      printValue(os, columnList_[i], printParameters_.delimiter_, 1, adjoint_values[i]);
    }
    outputManager_.getCommPtr()->barrier();
  }

  for (int i = 0; i < scaled_adjoint_values.size(); ++i)
  {
    if (outputManager_.getProcID() == 0)
    {
      printValue(os, columnList_[i], printParameters_.delimiter_, 1, scaled_adjoint_values[i]);
    }
    outputManager_.getCommPtr()->barrier();
  }

  if (outputManager_.getProcID() == 0)
  {
    os << std::endl;
  }

  ++index_;
}

//-----------------------------------------------------------------------------
// Function      : SensitivityTecPlot::doFinishOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
void SensitivityTecPlot::doFinishOutput()
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

  firstTimeSensitivity_ = true;
  index_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : SensitivityTecPlot::doFinishOutputStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
void SensitivityTecPlot::doFinishOutputStep()
{
  // close the sensitivity file.
  if (outStreamPtr_)
  {
    if ( outputManager_.getPrintEndOfSimulationLine() )
    {
      (*outStreamPtr_) << "End of Xyce(TM) Sensitivity Simulation" << std::endl;
    }
  }

  outputManager_.closeFile(outStreamPtr_);
  outStreamPtr_ = 0;
}

//-----------------------------------------------------------------------------
// Class         : TimeRaw
// Purpose       : Outputter class for transient output, rawfile output
//                 format
// Special Notes : Invoked by "FORMAT=raw" on .print line, not -r on command
//                 line.  -r is handled by the "OverrideRaw" classes.
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : TimeRaw::TimeRaw
// Purpose       : Constructor
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
TimeRaw::TimeRaw(OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_(""),
    outStreamPtr_( NULL),
    numPoints_( 0),
    numPointsPos_( 0),
    outputRAWTitleAndDate_(false)
{
  if (printParameters_.defaultExtension_.empty())
    printParameters_.defaultExtension_ = ".raw";

  fixupColumns(outputManager_, printParameters_, opList_);
}

//-----------------------------------------------------------------------------
// Function      : TimeRaw::~TimeRaw
// Purpose       : destructor
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
TimeRaw::~TimeRaw()
{
  outputManager_.closeFile(outStreamPtr_);

  deleteList(opList_.begin(), opList_.end());
}

//-----------------------------------------------------------------------------
// Function      : TimeRaw::doParse
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void  TimeRaw::doParse() {
  outFilename_ = outputFilename(printParameters_, outputManager_.getNetListFilename());
}


//-----------------------------------------------------------------------------
// Function      : TimeRaw::doOpen
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void  TimeRaw::doOpen()
{
  if (outputManager_.getProcID() == 0 && outStreamPtr_ == 0)
    outStreamPtr_ = outputManager_.openBinaryFile(outFilename_);
}

//-----------------------------------------------------------------------------
// Function      : TimeRaw::timeHeader
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void  TimeRaw::timeHeader()
{
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

      outputRAWTitleAndDate_ = true;

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
      for (std::vector<N_ANP_SweepParam>::const_iterator it = outputManager_.getStepParamVec().begin();
          it != outputManager_.getStepParamVec().end(); ++it)
      {
        plotName << " name = " << it->name << " value = " << it->currentVal << "  ";
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
  int numVars = 0;
  if (printParameters_.printType_ == PrintType::DC)
    ++numVars;

  numVars += opList_.size();

  // format number of internal and external variables included here + time
  if (outputManager_.getProcID() == 0)
  {
    os << "No. Variables: " << numVars << std::endl;

    // format total number of data points & remember the streampos
    os << "No. Points: ";
    numPointsPos_ = outStreamPtr_->tellp(); // <= 344 due to 80 char line limit
    os << "                  " << std::endl; // this will be overwritten

    if (outputManager_.getOutputVersionInRawFile() )
    {
      // spice3 does not output the version number of the simulator in the 
      // in the raw file header.  Optionally let one output the version 
      // if whatever program is going to read the file expects it
      os << "Version: " << N_UTL_Version::getFullVersionString() << std::endl;
    }

    // NOTE: dimensions, command, options, and scale not considered

    // write the variable information
    os << "Variables:" << std::endl;

    // The variable list appears next.  The format is:
    //         [tab](index) [tab](name) [tab](type)
    // Note that the(index) corresponds to the rawfile, not the soln vec.

    // write variable names for select data points
    std::string tmpNodeName, tmpType;
    int i = 0;
    if (printParameters_.printType_ == PrintType::DC)
    {
      os << "\t" << 0 << "\t" << "sweep\tvoltage\n";
      ++i;
    }

    for (Util::OpList::const_iterator it = opList_.begin() ; it != opList_.end(); ++it, ++i)
    {
      // set the type
      if (Util::hasExpressionTag((*it)->getName())) { tmpType = "expression"; }
      else if ((*it)->getName() == "INDEX")  { }
      else if ((*it)->getName() == "TIME")  { tmpType = "time"; }
      else if ((*it)->getName() == "FREQUENCY")  { tmpType = "frequency"; }
      else if ((*it)->getName()[0] == 'I')  { tmpType = "current";    }
      else if ((*it)->getName()[0] == 'V')  { tmpType = "voltage";    }
      else                              { tmpType = "unknown";    }

      // write the header line
      os << "\t" << i
         << "\t" << (*it)->getName()
         << "\t" << tmpType
         << "\n";

      // NOTE: other types, and params & dims, not considered
    }

    // write data marker
    // this string is actually ignored, but the pair of EOL markers is used
    os << "Binary:" << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : TimeRaw::doOutputTime
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void  TimeRaw::doOutputTime(
    const N_LAS_Vector * solnVecPtr, 
    const N_LAS_Vector * stateVecPtr, 
    const N_LAS_Vector * storeVecPtr)
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
    switch (printParameters_.printType_)
    {
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

  // select values to write from .PRINT line if FORMAT=RAW
  for (Util::OpList::const_iterator it = opList_.begin() ; it != opList_.end(); ++it)
  {
    // retrieve values from all procs
    double result = getValue(outputManager_.getCommPtr()->comm(),  *(*it), solnVecPtr, 0, 0, 0).real();
    if ((*it)->opType() == Xyce::Util::TIME_VAR)
      result *= printParameters_.outputTimeScaleFactor_;

    // file IO only on proc 0
    if (outputManager_.getProcID() == 0)
    {
      outStreamPtr_->write((char *) &result , sizeof(result));
    } // end proc0
  } // end for

  outputManager_.getCommPtr()->barrier();

  // keep track of number of datapoints
  ++numPoints_;
}

//-----------------------------------------------------------------------------
// Function      : TimeRaw::doFinishOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void  TimeRaw::doFinishOutput()
{
  if (outputManager_.getProcID() == 0)
  {
    if (outStreamPtr_ && numPoints_ != 0)
    {

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

//-----------------------------------------------------------------------------
// Class         : FrequencyRaw
// Purpose       : Outputter class for frequency-domain output, rawfile output
//                 format
// Special Notes : Invoked by "FORMAT=raw" on .print line, not -r on command
//                 line.  -r is handled by the "OverrideRaw" classes.
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : FrequencyRaw::FrequencyRaw
// Purpose       : Constructor
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
FrequencyRaw::FrequencyRaw(OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_(""),
    outStreamPtr_( NULL),
    numPoints_( 0),
    numPointsPos_( 0),
    outputRAWTitleAndDate_(false)
{
  if (printParameters_.defaultExtension_.empty())
    printParameters_.defaultExtension_ = ".raw";

  fixupColumns(outputManager_, printParameters_, opList_);
}

//-----------------------------------------------------------------------------
// Function      : FrequencyRaw::~FrequencyRaw
// Purpose       : destructor
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
FrequencyRaw::~FrequencyRaw()
{
  outputManager_.closeFile(outStreamPtr_);

  deleteList(opList_.begin(), opList_.end());
}

//-----------------------------------------------------------------------------
// Function      : FrequencyRaw::doParse
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void  FrequencyRaw::doParse()
{
  // prepare output manager to write
  numPoints_ = 0;

  outFilename_ = outputFilename(printParameters_, outputManager_.getNetListFilename());
}


//-----------------------------------------------------------------------------
// Function      : FrequencyRaw::doOpen
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void  FrequencyRaw::doOpen()
{
  if (outputManager_.getProcID() == 0 && outStreamPtr_ == 0)
    outStreamPtr_ = outputManager_.openBinaryFile(outFilename_);
}


//-----------------------------------------------------------------------------
// Function      : FrequencyRaw::frequencyHeader
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void  FrequencyRaw::frequencyHeader()
{
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
      for (std::vector<N_ANP_SweepParam>::const_iterator it = outputManager_.getStepParamVec().begin();
          it != outputManager_.getStepParamVec().end(); ++it)
      {
        plotName << " name = " << it->name << " value = " << it->currentVal << "  ";
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

    os << "Plotname: " << plotName.str() << std::endl;

    // format the flags
    std::string flags("complex");
    os << "Flags: " << flags << std::endl;

  } // end proc0 check

  int numVars = 0;
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
  numVars += opList_.size();

  if (outputManager_.getProcID() == 0)
  {
    // format number of internal and external variables included here + time
    os << "No. Variables: " << numVars << std::endl;

    // format total number of data points & remember the streampos
    os << "No. Points: ";
    numPointsPos_ = outStreamPtr_->tellp(); // <= 344 due to 80 char line limit
    os << "                  " << std::endl; // this will be overwritten

    if (outputManager_.getOutputVersionInRawFile() )
    {
      // spice3 does not output the version number of the simulator in the 
      // in the raw file header.  Optionally let one output the version 
      // if whatever program is going to read the file expects it
      os << "Version: " << N_UTL_Version::getFullVersionString() << std::endl;
    }

    // NOTE: dimensions, command, options, and scale not considered

    // write the variable information
    os << "Variables:" << std::endl;

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
      os << "\t" << 0 << "\t" << "sweep\tvoltage\n";
      ++i;
    }

    for (Util::OpList::const_iterator it = opList_.begin() ; it != opList_.end(); ++it, ++i)
    {
      std::string tmpNodeName, tmpType;
      // set the type
      if (Util::hasExpressionTag((*it)->getName())) { tmpType = "expression"; }
      else if ((*it)->getName() == "INDEX")  { }
      else if ((*it)->getName() == "TIME")  { tmpType = "time"; }
      else if ((*it)->getName() == "FREQUENCY")  { tmpType = "frequency"; }
      else if ((*it)->getName()[0] == 'I')  { tmpType = "current";    }
      else if ((*it)->getName()[0] == 'V')  { tmpType = "voltage";    }
      else                              { tmpType = "unknown";    }

      // write the header line
      os << "\t" << i
         << "\t" << (*it)->getName()
         << "\t" << tmpType
         << "\n";

      // NOTE: other types, and params & dims, not considered
    }

    // write data marker
    // this string is actually ignored, but the pair of EOL markers is used
    os << "Binary:" << std::endl;
  } // end proc0 check
}



//-----------------------------------------------------------------------------
// Function      : FrequencyRaw::doFinishOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void  FrequencyRaw::doFinishOutput()
{
  if (outputManager_.getProcID() == 0)
  {
    if (outStreamPtr_ && numPoints_ != 0)
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
  }

  // reset numPoints_ as it is used as a flag to print the header.
  numPoints_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : FrequencyRaw::doOutputFrequency
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
FrequencyRaw::doOutputFrequency(double frequency, const N_LAS_Vector *real_solution_vector, const N_LAS_Vector *imaginary_solution_vector)
{
  outputManager_.setCurrentOutputter(this);

  if (outputManager_.getProcID() == 0)
  {
    if (!outStreamPtr_)
    {
      doOpen();

      numPoints_ = 0;

      frequencyHeader();
    }
  }

  outputManager_.getCommPtr()->barrier();

  // select values to write from .PRINT line if FORMAT=RAW
  for (Util::OpList::const_iterator it = opList_.begin() ; it != opList_.end(); ++it)
  {
    complex result = getValue(outputManager_.getCommPtr()->comm(), *(*it), real_solution_vector, imaginary_solution_vector, 0, 0);
    if (outputManager_.getProcID() == 0)
    {
      double realPart=result.real();
      double imagPart=result.imag();
      outStreamPtr_->write((char *) &realPart, sizeof( double));
      outStreamPtr_->write((char *) &imagPart, sizeof( double));
    }
  }

  outputManager_.getCommPtr()->barrier();

  // keep track of number of datapoints
  ++numPoints_;
}

//-----------------------------------------------------------------------------
// Class         : TimeRawAscii
// Purpose       : Outputter class for transient output, rawfile output
//                 format, ascii version
// Special Notes : Invoked by "FORMAT=raw" on .print line, not -r on command
//                 line.  -r is handled by the "OverrideRaw" classes.
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : TimeRawAscii::TimeRawAscii
// Purpose       : Constructor
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
TimeRawAscii::TimeRawAscii(OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_( ""),
    outStreamPtr_( NULL),
    numPoints_(0),
    numPointsPos_( 0),
    outputRAWTitleAndDate_(false)
{
  if (printParameters_.defaultExtension_.empty())
    printParameters_.defaultExtension_ = ".raw";

  fixupColumns(outputManager_, printParameters_, opList_);
}

//-----------------------------------------------------------------------------
// Function      : TimeRawAscii::~TimeRawAscii
// Purpose       : Destructor
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
TimeRawAscii::~TimeRawAscii()
{
  outputManager_.closeFile(outStreamPtr_);

  deleteList(opList_.begin(), opList_.end());
}

//-----------------------------------------------------------------------------
// Function      : TimeRawAscii::doParse
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void  TimeRawAscii::doParse()
{
  outFilename_ = outputFilename(printParameters_, outputManager_.getNetListFilename());
}

//-----------------------------------------------------------------------------
// Function      : TimeRawAscii::doOpen
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void  TimeRawAscii::doOpen()
{
  if (outputManager_.getProcID() == 0 && outStreamPtr_ == 0)
  {
    outStreamPtr_ = outputManager_.openFile(outFilename_);

    // set output value characteristics
    outStreamPtr_->setf( std::ios::scientific);
    outStreamPtr_->precision(8);
    outStreamPtr_->setf( std::ios::left, std::ios::adjustfield);
  }
}

//-----------------------------------------------------------------------------
// Function      : TimeRawAscii::timeHeader
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void  TimeRawAscii::timeHeader()
{
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

      outputRAWTitleAndDate_ = true;

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
      for (std::vector<N_ANP_SweepParam>::const_iterator it = outputManager_.getStepParamVec().begin(); it != outputManager_.getStepParamVec().end(); ++it)
      {
        plotName << " name = " << it->name << " value = " << it->currentVal << "  ";
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
  int numVars = 0;
  if (printParameters_.printType_ == PrintType::DC)
    ++numVars;

  numVars += opList_.size();

  if (outputManager_.getProcID() == 0)
  {
    // format number of internal and external variables included here + time
    os << "No. Variables: " << numVars << std::endl;

    // format total number of data points & remember the streampos
    os << "No. Points: ";
    numPointsPos_ = outStreamPtr_->tellp(); // <= 344 due to 80 char line limit
    os << "                  " << std::endl; // this will be overwritten

    if (outputManager_.getOutputVersionInRawFile() )
    {
      // spice3 does not output the version number of the simulator in the 
      // in the raw file header.  Optionally let one output the version 
      // if whatever program is going to read the file expects it
      os << "Version: " << N_UTL_Version::getFullVersionString() << std::endl;
    }

    // NOTE: dimensions, command, options, and scale not considered

    // write the variable information
    os << "Variables:" << std::endl;

    // The variable list appears next.  The format is:
    //         [tab](index) [tab](name) [tab](type)
    // Note that the(index) corresponds to the rawfile, not the soln vec.

    // write variable names for select data points
    std::string tmpNodeName, tmpType;
    int i = 0;
    if (printParameters_.printType_ == PrintType::DC)
    {
      // add timestep header info(Spice3f5 style)
      os << "\t" << 0 << "\t" << "sweep\tvoltage\n";
      ++i;
    }

    for (Util::OpList::const_iterator it = opList_.begin() ; it != opList_.end(); ++it, ++i)
    {
      // set the type
      if (Util::hasExpressionTag((*it)->getName())) { tmpType = "expression"; }
      else if ((*it)->getName() == "INDEX")  { }
      else if ((*it)->getName() == "TIME")  { tmpType = "time"; }
      else if ((*it)->getName() == "FREQUENCY")  { tmpType = "frequency"; }
      else if ((*it)->getName()[0] == 'I')  { tmpType = "current";    }
      else if ((*it)->getName()[0] == 'V')  { tmpType = "voltage";    }
      else                              { tmpType = "unknown";    }

      // write the header line
      os << "\t" << i
         << "\t" << (*it)->getName()
         << "\t" << tmpType
         << "\n";

      // NOTE: other types, and params & dims, not considered
    }

    // write data marker
    os << "Values:" << std::endl;
  } // end proc0 check
}

//-----------------------------------------------------------------------------
// Function      : TimeRawAscii::doOutputTime
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void  TimeRawAscii::doOutputTime(
    const N_LAS_Vector * solnVecPtr, 
    const N_LAS_Vector * stateVecPtr, 
    const N_LAS_Vector * storeVecPtr)
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
    // write the index to ascii rawfile
    os << numPoints_;
  }

  // select values to write from .PRINT line if FORMAT=RAW
  for (Util::OpList::const_iterator it = opList_.begin() ; it != opList_.end(); ++it)
  {
    // retrieve values from all procs
    double result = getValue(outputManager_.getCommPtr()->comm(), *(*it), solnVecPtr, 0, 0, 0).real();
    if ((*it)->opType() == Xyce::Util::TIME_VAR)
      result *= printParameters_.outputTimeScaleFactor_;

    // file IO only on proc 0
    if (outputManager_.getProcID() == 0)
    {
      os << "\t" << result << "\n";
    }
  }

  if (outputManager_.getProcID() == 0)
  {
    os << std::endl;
  }

  outputManager_.getCommPtr()->barrier();

  // keep track of number of datapoints
  ++numPoints_;
}
//-----------------------------------------------------------------------------
// Function      : TimeRawAscii::doFinishOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void  TimeRawAscii::doFinishOutput()
{
  if (outputManager_.getProcID() == 0)
  {
    if (outStreamPtr_ && numPoints_ != 0)
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
  }

  // reset numPoints_ as it is used as a flag
  // by outputRAW_() to print the header.
  numPoints_=0;
}


//-----------------------------------------------------------------------------
// Function      : TimeRawAscii::doFinishOutputStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
TimeRawAscii::doFinishOutputStep()
{
}

//-----------------------------------------------------------------------------
// Class         : FrequencyRawAscii
// Purpose       : Outputter class for frequency-domain output, rawfile output
//                 format, ascii version
// Special Notes : Invoked by "FORMAT=raw" on .print line, not -r on command
//                 line.  -r is handled by the "OverrideRaw" classes.
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : FrequencyRawAscii::FrequencyRawAscii
// Purpose       : Constructor
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
FrequencyRawAscii::FrequencyRawAscii(OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_( ""),
    outStreamPtr_( NULL),
    numPoints_(0),
    numPointsPos_( 0),
    outputRAWTitleAndDate_(false)
{
  if (printParameters_.defaultExtension_.empty())
    printParameters_.defaultExtension_ = ".raw";

  fixupColumns(outputManager_, printParameters_, opList_);
}

//-----------------------------------------------------------------------------
// Function      : FrequencyRawAscii::~FrequencyRawAscii
// Purpose       : destructor
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
FrequencyRawAscii::~FrequencyRawAscii()
{
  outputManager_.closeFile(outStreamPtr_);

  deleteList(opList_.begin(), opList_.end());
}

//-----------------------------------------------------------------------------
// Function      : FrequencyRawAscii::doParse
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void  FrequencyRawAscii::doParse()
{
  outFilename_ = outputFilename(printParameters_, outputManager_.getNetListFilename());
}

//-----------------------------------------------------------------------------
// Function      : FrequencyRawAscii::doOpen
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void  FrequencyRawAscii::doOpen()
{
  if (outputManager_.getProcID() == 0 && outStreamPtr_ == 0)
  {
    outStreamPtr_ = outputManager_.openFile(outFilename_);

    // set output value characteristics
    outStreamPtr_->setf( std::ios::scientific);
    outStreamPtr_->precision(8);
    outStreamPtr_->setf( std::ios::left, std::ios::adjustfield);
  }
}

//-----------------------------------------------------------------------------
// Function      : FrequencyRawAscii::frequencyHeader
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void  FrequencyRawAscii::frequencyHeader()
{
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

      outputRAWTitleAndDate_ = true;

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
      for (std::vector<N_ANP_SweepParam>::const_iterator it = outputManager_.getStepParamVec().begin(); it != outputManager_.getStepParamVec().end(); ++it)
      {
        plotName << " name = " << it->name << " value = " << it->currentVal << "  ";
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

    os << "Plotname: " << plotName.str() << std::endl;

    // format the flags
    std::string flags("complex");
    os << "Flags: " << flags << std::endl;
  }

  int numVars = 0;
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
  numVars += opList_.size();

  if (outputManager_.getProcID() == 0)
  {
    // format number of internal and external variables included here + time
    os << "No. Variables: " << numVars << std::endl;

    // format total number of data points & remember the streampos
    os << "No. Points: ";
    numPointsPos_ = outStreamPtr_->tellp(); // <= 344 due to 80 char line limit
    os << "                  " << std::endl; // this will be overwritten

    if (outputManager_.getOutputVersionInRawFile() )
    {
      // spice3 does not output the version number of the simulator in the 
      // in the raw file header.  Optionally let one output the version 
      // if whatever program is going to read the file expects it
      os << "Version: " << N_UTL_Version::getFullVersionString() << std::endl;
    }


    // NOTE: dimensions, command, options, and scale not considered

    // write the variable information
    os << "Variables:" << std::endl;

    // The variable list appears next.  The format is:
    //         [tab](index) [tab](name) [tab](type)
    // Note that the(index) corresponds to the rawfile, not the soln vec.

    // write variable names for select data points
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
      os << "\t" << 0 << "\t" << "sweep\tvoltage\n";
      ++i;
    }

    for (Util::OpList::const_iterator it = opList_.begin() ; it != opList_.end(); ++it, ++i)
    {
      std::string tmpNodeName, tmpType;
      // set the type
      if (Util::hasExpressionTag((*it)->getName())) { tmpType = "expression"; }
      else if ((*it)->getName() == "INDEX")  { }
      else if ((*it)->getName() == "TIME")  { tmpType = "time"; }
      else if ((*it)->getName() == "FREQUENCY")  { tmpType = "frequency"; }
      else if ((*it)->getName()[0] == 'I')  { tmpType = "current";    }
      else if ((*it)->getName()[0] == 'V')  { tmpType = "voltage";    }
      else                              { tmpType = "unknown";    }

      // write the header line
      os << "\t" << i
         << "\t" << (*it)->getName()
         << "\t" << tmpType
         << "\n";

      // NOTE: other types, and params & dims, not considered
    }

    // write data marker
    os << "Values:" << std::endl;
  } // end proc0 check
}

//-----------------------------------------------------------------------------
// Function      : FrequencyRawAscii::doOutputFrequency
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
FrequencyRawAscii::doOutputFrequency(double frequency, const N_LAS_Vector *real_solution_vector, const N_LAS_Vector *imaginary_solution_vector)
{
  outputManager_.setCurrentOutputter(this);
  if (outputManager_.getProcID() == 0)
  {
    if (!outStreamPtr_)
    {
      doOpen();

      numPoints_ = 0;

      frequencyHeader();
    }
  }

  // file IO on proc 0 only
  if (outputManager_.getProcID() == 0)
  {

    // write the index to ascii rawfile
    (*outStreamPtr_) << numPoints_;

  }

  outputManager_.getCommPtr()->barrier();

  // select values to write from .PRINT line if FORMAT=RAW
  for (Util::OpList::const_iterator it = opList_.begin() ; it != opList_.end(); ++it)
  {
    complex result = getValue(outputManager_.getCommPtr()->comm(), *(*it), real_solution_vector, imaginary_solution_vector, 0, 0);
    if (outputManager_.getProcID() == 0)
    {
      (*outStreamPtr_) << "\t"  << result.real() << ", " << result.imag() << "\n";
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

//-----------------------------------------------------------------------------
// Function      : FrequencyRawAscii::doFinishOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void  FrequencyRawAscii::doFinishOutput()
{
  if (outputManager_.getProcID() == 0)
  {
    if (outStreamPtr_)
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
  }

  // reset numPoints_ as it is used as a flag
  // by outputRAW_() to print the header.
  numPoints_=0;
}


//-----------------------------------------------------------------------------
// Function      : FrequencyRawAscii::doFinishOutputStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
FrequencyRawAscii::doFinishOutputStep()
{
}

//-----------------------------------------------------------------------------
// Class         : OverrideRaw
// Purpose       : Outputter class for rawfile output format
// Special Notes : Invoked by -r on command line.  FORMAT=RAW on print line
//                 is handled by the TimeRaw* and FrequencyRaw* classes
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : OverrideRaw::OverrideRaw
// Purpose       : Constructor
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
OverrideRaw::OverrideRaw(OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_( ""),
    outStreamPtr_( NULL),
    numPoints_( 0),
    numPointsPos_( 0),
    outputRAWTitleAndDate_(false)
{
  if (printParameters_.defaultExtension_.empty())
    printParameters_.defaultExtension_ = ".raw";
}

//-----------------------------------------------------------------------------
// Function      : OverrideRaw::~OverrideRaw
// Purpose       : destructor
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
OverrideRaw::~OverrideRaw()
{
  outputManager_.closeFile(outStreamPtr_);

  deleteList(opList_.begin(), opList_.end());
}

//-----------------------------------------------------------------------------
// Function      : OverrideRaw::doParse
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void  OverrideRaw::doParse()
{
  outFilename_ = outputFilename(printParameters_, outputManager_.getNetListFilename());
}

//-----------------------------------------------------------------------------
// Function      : OverrideRaw::doOpen
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void  OverrideRaw::doOpen()
{
  if (outputManager_.getProcID() == 0 && outStreamPtr_ == 0)
  {
    outStreamPtr_ = outputManager_.openBinaryFile(outFilename_);

    // set output value characteristics
    outStreamPtr_->setf( std::ios::scientific);
    outStreamPtr_->precision(8);
    outStreamPtr_->setf( std::ios::left, std::ios::adjustfield);
  }
}


//-----------------------------------------------------------------------------
// Function      : OverrideRaw::timeHeader
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void  OverrideRaw::timeHeader()
{
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

      outputRAWTitleAndDate_ = true;

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
      for (std::vector<N_ANP_SweepParam>::const_iterator it = outputManager_.getStepParamVec().begin(); it != outputManager_.getStepParamVec().end(); ++it)
      {
        plotName << " name = " << it->name << " value = " << it->currentVal << "  ";
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
  // get node types
  std::vector< char > typeRefs;
  outputManager_.getTopPtr()->returnVarTypeVec( typeRefs);

  // count data points
  int numVars = outputManager_.getAllNodes().size();

  std::vector< std::pair< std::string, char > > nT;

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
      nT.push_back( std::pair< std::string, char>((*name_i).first,
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
        nT.push_back( std::pair< std::string, char >( name, t));
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
    os << "No. Variables: " << numVars + 1 << std::endl;

  outputManager_.getCommPtr()->barrier();

  if (outputManager_.getProcID() == 0)
  {
    // format total number of data points & remember the streampos
    os << "No. Points: ";
    numPointsPos_ = outStreamPtr_->tellp(); // <= 344 due to 80 char line limit
    os << "                  " << std::endl; // this will be overwritten

    if (outputManager_.getOutputVersionInRawFile() )
    {
      // spice3 does not output the version number of the simulator in the 
      // in the raw file header.  Optionally let one output the version 
      // if whatever program is going to read the file expects it
      os << "Version: " << N_UTL_Version::getFullVersionString() << std::endl;
    }

    // NOTE: dimensions, command, options, and scale not considered

    // The variable list appears next.  The format is:
    //         [tab](index) [tab](name) [tab](type)
    // Note that the(index) corresponds to the rawfile, not the soln vec.

    // write variable names for select data points
    std::string tmpType;
    // write the variable information
    os << "Variables:" << std::endl;

    // add timestep header info(Spice3f5 style)
    os << "\t" << 0 << "\t";

    if (printParameters_.printType_ == PrintType::TRAN)
    {
      os << "time\ttime\n";
    }
    else if (printParameters_.printType_ == PrintType::AC)
    {
      os << "frequency\tfrequency\n";
    }
    else
    {
      os << "sweep\tvoltage\n";
    }

    std::string::size_type uPos;
    for ( int i = 0; i < numVars; ++i)
    {
      std::string tmpNodeName =(nT[i]).first;

      // format is:  [tab](index) [tab](name) [tab](type)
      //   the index corresponds to the rawfile, not the soln vec
      if (strspn( tmpNodeName.c_str(), "0123456789") == tmpNodeName.size())
      {
        // sconvert, spice3f5, & chilespice wrap numeric voltage node names in V()
        tmpNodeName = "V(" +(nT[i]).first + ")";
      }

      uPos = tmpNodeName.rfind( "_", tmpNodeName.size());
      if (uPos != std::string::npos)
      {
        tmpNodeName.replace( uPos, 1, "#");
      }

      os << "\t" << i + 1 << "\t" << tmpNodeName << "\t";

      if ((nT[i]).second == 'I' ||(nT[i]).second == 'i')
      {
        os << "current\n" ;
      }
      else
      {
        os << "voltage\n";
      }

      // NOTE: other types, and params & dims, not considered
    }

    // write data marker
    // this string is actually ignored, but the pair of EOL markers is used
    os << "Binary:" << std::endl;
  } // end proc0 check
}

//-----------------------------------------------------------------------------
// Function      : OverrideRaw::frequencyHeader
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void  OverrideRaw::frequencyHeader()
{
  // This routine was split off of doOutputHeader because an
  // AC analysis can output both real data "doOutput" and complex
  // data via doOutputAC.  Each of these doOutput routines must
  // call a doHeader routine which can specify the right underlying
  // data type.  So rather than add more conditional logic we chose
  // to separate out the functionality for the different use cases.
  int numVars;
  std::vector< std::pair< std::string, char > > nT;
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
      for (std::vector<N_ANP_SweepParam>::const_iterator it = outputManager_.getStepParamVec().begin(); it != outputManager_.getStepParamVec().end(); ++it)
      {
        plotName << " name = " << it->name << " value = " << it->currentVal << "  ";
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
      nT.push_back( std::pair< std::string, char>((*name_i).first,
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
        nT.push_back( std::pair< std::string, char >( name, t));
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

    if (outputManager_.getOutputVersionInRawFile() )
    {
      // spice3 does not output the version number of the simulator in the 
      // in the raw file header.  Optionally let one output the version 
      // if whatever program is going to read the file expects it
      (*outStreamPtr_) << "Version: " << N_UTL_Version::getFullVersionString() << std::endl;
    }

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

      uPos = tmpNodeName.rfind( "_", tmpNodeName.size());
      if (uPos != std::string::npos)
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


//-----------------------------------------------------------------------------
// Function      : OverrideRaw::doOutputTime
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void  OverrideRaw::doOutputTime(
    const N_LAS_Vector * solnVecPtr, 
    const N_LAS_Vector * stateVecPtr, 
    const N_LAS_Vector * storeVecPtr)
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
      double result = (*solnVecPtr)[(*name_i).second.first];
      result = filter(result, printParameters_.filter_);

      // write binary data to rawfile
      outStreamPtr_->write((char *)&result , sizeof( double));
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
      double result =(*solnVecPtr)[(*name_i).second.first];
      result = filter(result, printParameters_.filter_);

      outputManager_.getCommPtr()->pack( &result, 1, b, bSize, pos);
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

//-----------------------------------------------------------------------------
// Function      : OverrideRaw::doResetOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void  OverrideRaw::doResetOutput()
{
}

//-----------------------------------------------------------------------------
// Function      : OverrideRaw::doFinishOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void  OverrideRaw::doFinishOutput()
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

//-----------------------------------------------------------------------------
// Function      : OverrideRaw::doOutputFrequency
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
OverrideRaw::doOutputFrequency(double frequency, const N_LAS_Vector *real_solution_vector, const N_LAS_Vector *imaginary_solution_vector)
{
  outputManager_.setCurrentOutputter(this);

  if (numPoints_ == 0)
  {
    doOpen();

    frequencyHeader();
  }

  outputManager_.getCommPtr()->barrier();

  double result = frequency;
  // file IO only on proc 0
  if (outputManager_.getProcID() == 0)
  {
    outStreamPtr_->write((char *) &result, sizeof(double));
    result = 0.0;
    outStreamPtr_->write((char *) &result, sizeof(double));
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
      complex result = complex((*real_solution_vector)[(*name_i).second.first], (*imaginary_solution_vector)[(*name_i).second.first]);

      // write formatted values to rawfile
      if (outputManager_.getProcID() == 0)
      {
        double realPart=result.real();
        double imagPart=result.imag();
        outStreamPtr_->write((char *) &realPart, sizeof(double));
        outStreamPtr_->write((char *) &imagPart, sizeof(double));
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
        double result = 0.0;
        outputManager_.getCommPtr()->unpack( b, bSize, pos, &result, 1);

        // write formatted values to rawfile
        if (outputManager_.getProcID() == 0)
        {
          outStreamPtr_->write((char *) &result, sizeof(double));
          result = 0.0;
          outStreamPtr_->write((char *) &result, sizeof(double));
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
      double result =(*real_solution_vector)[(*name_i).second.first];
      result = filter(result, printParameters_.filter_);
      outputManager_.getCommPtr()->pack( &result, 1, b, bSize, pos);
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

//-----------------------------------------------------------------------------
// Function      : OverrideRaw::doOutputHB
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void  OverrideRaw::doOutputHB(
  const std::vector<double>&    timePoints,
  const std::vector<double>&    freqPoints,
  const N_LAS_BlockVector &     timeDomainSolnVec,
  const N_LAS_BlockVector &     freqDomainSolnVecReal,
  const N_LAS_BlockVector &     freqDomainSolnVecImaginary)
{
}

//-----------------------------------------------------------------------------
// Function      : OverrideRaw::doOutputMPDE
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void  OverrideRaw::doOutputMPDE(double time, const N_LAS_Vector *solution_vector)
{
}

//-----------------------------------------------------------------------------
// Function      : OverrideRaw::doOutputHomotopy
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void  OverrideRaw::doOutputHomotopy(const std::vector<std::string> & parameter_names, const std::vector<double> & parameter_values, const N_LAS_Vector * solution_vector)
{
}

//-----------------------------------------------------------------------------
// Function      : OverrideRaw::doOutputMORTF
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void  OverrideRaw::doOutputMORTF(bool origSystem, const double & freq, const Teuchos::SerialDenseMatrix<int, std::complex<double> >& H)
{
}

//-----------------------------------------------------------------------------
// Function      : OverrideRaw::doFinishOutputStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
OverrideRaw::doFinishOutputStep()
{
}


//-----------------------------------------------------------------------------
// Class         : OverrideRawAscii
// Purpose       : Outputter class for rawfile output format, ascii version
// Special Notes : Invoked by -r on command line.  FORMAT=RAW on print line
//                 is handled by the TimeRaw* and FrequencyRaw* classes
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : OverrideRawAscii::OverrideRawAscii
// Purpose       : Constructor
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
OverrideRawAscii::OverrideRawAscii(OutputMgr &output_manager, const PrintParameters &print_parameters)
  : outputManager_(output_manager),
    printParameters_(print_parameters),
    outFilename_( ""),
    outStreamPtr_( NULL),
    numPoints_( 0),
    numPointsPos_( 0),
    outputRAWTitleAndDate_(false)
{
  if (printParameters_.defaultExtension_.empty())
    printParameters_.defaultExtension_ = ".raw";
}

//-----------------------------------------------------------------------------
// Function      : OverrideRawAscii::~OverrideRawAscii
// Purpose       : destructor
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
OverrideRawAscii::~OverrideRawAscii()
{
  outputManager_.closeFile(outStreamPtr_);

  deleteList(opList_.begin(), opList_.end());
}

//-----------------------------------------------------------------------------
// Function      : OverrideRawAscii::doParse
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void  OverrideRawAscii::doParse()
{
  outFilename_ = outputFilename(printParameters_, outputManager_.getNetListFilename());
}

//-----------------------------------------------------------------------------
// Function      : OverrideRawAscii::doOpen
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void  OverrideRawAscii::doOpen()
{
  if (outputManager_.getProcID() == 0 && outStreamPtr_ == 0)
  {
    outStreamPtr_ = outputManager_.openFile(outFilename_);

    // set output value characteristics
    outStreamPtr_->setf( std::ios::scientific);
    outStreamPtr_->precision(8);
    outStreamPtr_->setf( std::ios::left, std::ios::adjustfield);
  }
}

//-----------------------------------------------------------------------------
// Function      : OverrideRawAscii::timeHeader
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void  OverrideRawAscii::timeHeader()
{
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
      for (std::vector<N_ANP_SweepParam>::const_iterator it = outputManager_.getStepParamVec().begin(); it != outputManager_.getStepParamVec().end(); ++it)
      {
        plotName << " name = " << it->name << " value = " << it->currentVal << "  ";
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
  int numVars = outputManager_.getAllNodes().size();
  std::vector< std::pair< std::string, char > > nT;

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
      nT.push_back( std::pair< std::string, char>((*name_i).first,
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
        nT.push_back( std::pair< std::string, char >( name, t));
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

    if (outputManager_.getOutputVersionInRawFile() )
    {
      // spice3 does not output the version number of the simulator in the 
      // in the raw file header.  Optionally let one output the version 
      // if whatever program is going to read the file expects it
      (*outStreamPtr_) << "Version: " << N_UTL_Version::getFullVersionString() << std::endl;
    }

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

      uPos = tmpNodeName.rfind( "_", tmpNodeName.size());
      if (uPos != std::string::npos)
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


//-----------------------------------------------------------------------------
// Function      : OverrideRawAscii::frequencyHeader
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void  OverrideRawAscii::frequencyHeader()
{
  // This routine was split off of doOutputHeader because an
  // AC analysis can output both real data "doOutput" and complex
  // data via doOutputAC.  Each of these doOutput routines must
  // call a doHeader routine which can specify the right underlying
  // data type.  So rather than add more conditional logic we chose
  // to separate out the functionality for the different use cases.
  int numVars;
  std::vector< std::pair< std::string, char > > nT;
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
      for (std::vector<N_ANP_SweepParam>::const_iterator it = outputManager_.getStepParamVec().begin();
          it != outputManager_.getStepParamVec().end(); ++it)
      {
        plotName << " name = " << it->name << " value = " << it->currentVal << "  ";
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
      nT.push_back( std::pair< std::string, char>((*name_i).first,
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
        nT.push_back( std::pair< std::string, char >( name, t));
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

    if (outputManager_.getOutputVersionInRawFile() )
    {
      // spice3 does not output the version number of the simulator in the 
      // in the raw file header.  Optionally let one output the version 
      // if whatever program is going to read the file expects it
      (*outStreamPtr_) << "Version: " << N_UTL_Version::getFullVersionString() << std::endl;
    }

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

      uPos = tmpNodeName.rfind( "_", tmpNodeName.size());
      if (uPos != std::string::npos)
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

//-----------------------------------------------------------------------------
// Function      : OverrideRawAscii::doOutputTime
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void  OverrideRawAscii::doOutputTime(
    const N_LAS_Vector * solnVecPtr, 
    const N_LAS_Vector * stateVecPtr, 
    const N_LAS_Vector * storeVecPtr)
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
    if ( printParameters_.printType_ == PrintType::TRAN)
    {
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
      double result = (*solnVecPtr)[(*name_i).second.first];
      result = filter(result, printParameters_.filter_);

      os << "\t" << result << "\n";
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
      double result = (*solnVecPtr)[(*name_i).second.first];
      result = filter(result, printParameters_.filter_);

      outputManager_.getCommPtr()->pack( &result, 1, b, bSize, pos);
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

//-----------------------------------------------------------------------------
// Function      : OverrideRawAscii::doOutputFrequency
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
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
      complex result  = complex((*real_solution_vector)[(*name_i).second.first], (*imaginary_solution_vector)[(*name_i).second.first]);

      (*outStreamPtr_) << "\t" << result.real() << ", " << result.imag() << "\n";
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
      double result = (*real_solution_vector)[(*name_i).second.first];
      result = filter(result, printParameters_.filter_);

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

//-----------------------------------------------------------------------------
// Function      : OverrideRawAscii::doOutputHB
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void  OverrideRawAscii::doOutputHB(
  const std::vector<double>&    timePoints,
  const std::vector<double>&    freqPoints,
  const N_LAS_BlockVector &     timeDomainSolnVec,
  const N_LAS_BlockVector &     freqDomainSolnVecReal,
  const N_LAS_BlockVector &     freqDomainSolnVecImaginary)
{
}

//-----------------------------------------------------------------------------
// Function      : OverrideRawAscii::doOutputMPDE
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void  OverrideRawAscii::doOutputMPDE(double time, const N_LAS_Vector *solution_vector)
{
}

//-----------------------------------------------------------------------------
// Function      : OverrideRawAscii::doOutputHomotopy
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void  OverrideRawAscii::doOutputHomotopy(const std::vector<std::string> & parameter_names, const std::vector<double> & parameter_values, const N_LAS_Vector * solution_vector)
{
}

//-----------------------------------------------------------------------------
// Function      : OverrideRawAscii::doOutputMORTF
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void  OverrideRawAscii::doOutputMORTF(bool origSystem, const double & freq, const Teuchos::SerialDenseMatrix<int, std::complex<double> >& H)
{
}

//-----------------------------------------------------------------------------
// Function      : OverrideRawAscii::doResetOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void  OverrideRawAscii::doResetOutput()
{
}

//-----------------------------------------------------------------------------
// Function      : OverrideRawAscii::doFinishOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void  OverrideRawAscii::doFinishOutput()
{
  if (outputManager_.getProcID() == 0)
  {
    if (outStreamPtr_)
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
  }

  // reset numPoints_ as it is used as a flag
  // by outputRAW_() to print the header.
  numPoints_=0;
}


//-----------------------------------------------------------------------------
// Function      : OverrideRawAscii::doFinishOutputStep
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
OverrideRawAscii::doFinishOutputStep()
{
}

//-----------------------------------------------------------------------------
// Class         : MOR
// Purpose       : Outputter class for MOR runs
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date : 6/7/2013
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Function      : MOR::MOR
// Purpose       : Constructor
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
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
  if (printParameters_.defaultExtension_.empty())
    printParameters_.defaultExtension_ = ".prn";
}

//-----------------------------------------------------------------------------
// Function      : MOR::~MOR
// Purpose       : destructor
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
MOR::~MOR()
{
  if (outStreamPtr_ != &Xyce::dout())
  {
    delete outStreamPtr_;
  }

  deleteList(opList_.begin(), opList_.end());
}

//-----------------------------------------------------------------------------
// Function      : MOR::doParse
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void  MOR::doParse()
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

//-----------------------------------------------------------------------------
// Function      : MOR::doOpen
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void  MOR::doOpen()
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
void  MOR::outputMORHeaders(int numPorts)
{
  index_ = 0;

#ifdef Xyce_DEBUG_IO
  Xyce::dout() << "MOR Header output:" << std::endl;
#endif

  // hardwire the AC print format to tecplot or gnuplot(not PROBE!)
  if (format_ != Format::TECPLOT && format_ != Format::STD)
  {
    Report::UserWarning() << "MOR output can only use tecplot format, resetting to gnuplot format.";

    format_ = Format::STD;
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

  outStreamPtr_->setf(std::ios::scientific);
  outStreamPtr_->precision(printParameters_.streamPrecision_);
  outStreamPtr_->setf(std::ios::left, std::ios::adjustfield);

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
void  MOR::doOutputMORTF(bool origSystem, const double & freq, const Teuchos::SerialDenseMatrix<int, std::complex<double> >& H)
{
  bool firstColPrinted = false;

#ifdef Xyce_DEBUG_IO
  Xyce::dout() << "Begin outputMORTF" << std::endl;
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
      if (firstColPrinted)
      {
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
            if (printParameters_.delimiter_ != "")
            {
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
            if (printParameters_.delimiter_ != "")
            {
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
void  MOR::stdFreqMORHeader_(std::ostream & stream, int numPorts)
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
// Scope         :
// Creator       : Heidi Thornquist and Ting Mei
// Creation Date : 5/31/12
//-----------------------------------------------------------------------------
void  MOR::tecplotFreqMORHeader_(std::ostream & stream, int counter, int numPorts)
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
  stream.flush();
}

//-----------------------------------------------------------------------------
// Function      : MOR::doResetOutput
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void  MOR::doResetOutput()
{
  firstTimePrint_ = true;
  stepCount_ = 0;
}

} // namespace outputter

//-----------------------------------------------------------------------------
// Namespace     : Unnamed
// Purpose       : file-local scoped methods and data
// Special Notes : This second unnamed namespace block defines some functions
//                 that had only been declared in the previous block
// Creator       : dgbaur
// Creation Date : 06/28/2013
//-----------------------------------------------------------------------------
namespace { // <unnamed>

//-----------------------------------------------------------------------------
// Function      : tecplotFreqHeader
// Purpose       : header for tecplot. Freq Domain(default)
// Special Notes :
// Scope         : file-local
// Creator       : Eric Keiter
// Creation Date : 4/11/12
//-----------------------------------------------------------------------------
void tecplotFreqHeader(OutputMgr &output_manager, const PrintParameters &print_parameters, const Util::OpList &op_list, std::ostream &os, int counter)
{
  static const int tecplotHeaderPrecision_ = 2;
  os.setf(std::ios::scientific);
  os.precision( tecplotHeaderPrecision_);

  if (counter == 0)
  {
    os << " TITLE = \" Xyce Frequency Domain data, " << output_manager.getNetListFilename() << "\", " << std::endl;
    os << "\tVARIABLES = ";

    // output the user-specified solution vars:
    for (Util::OpList::const_iterator it = op_list.begin() ; it != op_list.end(); ++it)
    {
      os << "\" ";
      if ( (*it)->getName() == "FREQUENCY" )
      {
        os << "FREQ";
      }
      else
      {
        os << (*it)->getName() ;
      }
      os << "\" " << std::endl;
    }
    os << "DATASETAUXDATA ";
    os << getTecplotTimeDateStamp();
    os << std::endl;

    if (!output_manager.getTempSweepFlag())
    {
      os << "DATASETAUXDATA TEMP = \"" << output_manager.getCircuitTemp() << " \"" << std::endl;
    }
  }

  // output some AUXDATA
  os << "ZONE F=POINT  ";

  if (output_manager.getStepParamVec().empty())
  {
    os << " T=\"Xyce data\" ";
  }
  else
  {
    os << " T= \" ";
    for (std::vector<N_ANP_SweepParam>::const_iterator it = output_manager.getStepParamVec().begin(); it != output_manager.getStepParamVec().end(); ++it)
    {
      os << " " << it->name << " = " << it->currentVal;
    }
    os << "\" ";
  }

  os << std::endl;

  // put in the various sweep parameters as auxdata:
  if (!output_manager.getStepParamVec().empty())
  {
    for (std::vector<N_ANP_SweepParam>::const_iterator iterParam = output_manager.getStepParamVec().begin();
    iterParam != output_manager.getStepParamVec().end();
    ++iterParam)
    {
      // convert any ":" or "%" in the name to a "_", so as not to confuse tecplot.
      std::string tmpName(iterParam->name);
      replace(tmpName.begin(), tmpName.end(), '%', '_');
      replace(tmpName.begin(), tmpName.end(), ':', '_');
      os << "AUXDATA " << tmpName << " = " << "\" " << iterParam->currentVal << "\" ";
    }
    os << std::endl;
  }

  os << std::flush; 
}

//-----------------------------------------------------------------------------
// Function      : printHeader
// Purpose       : Given print parameters and a stream, print the header
// Special Notes : top level function
// Scope         : file-local
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
std::ostream &printHeader(std::ostream &os, const PrintParameters &print_parameters)
{
  return printHeader(os, print_parameters.table_.columnList_, print_parameters.delimiter_);
}


//-----------------------------------------------------------------------------
// Function      : printHeader
// Purpose       : Given stream, column list, and delimiter, print header
// Special Notes :
// Scope         : file-local
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
std::ostream &printHeader(std::ostream &os, const Table::ColumnList &column_list, const std::string &delimiter)
{
  for (Table::ColumnList::const_iterator it = column_list.begin(); it != column_list.end(); ++it)
  {
    if (it != column_list.begin())
      os << (delimiter.empty() ? " " : delimiter);

    printHeader(os, (*it));
  }

  os << std::endl;

  return os;
}


//-----------------------------------------------------------------------------
// Function      : printHeader
// Purpose       : print a single column of the header on the given stream.
// Special Notes :
// Scope         : file-local
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
std::ostream &printHeader(std::ostream &os, const Table::Column &column)
{
  std::string name = column.name_;
  if (name == "INDEX")
    name = "Index";

  size_t left_padding = 0;
  size_t right_padding = 0;

  // if (column.width_ < name.size())
  //   column.width_ = name.size();

  if (column.width_ > name.size())
  {
    switch (column.justification_)
    {
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

//-----------------------------------------------------------------------------
// Function      : createOps
// Purpose       : given an output manager and a (print) parameter list, and
//                 a back_inserter iterator for an OpList, generate all the
//                 "Ops" needed to obtain the values in the print parameter
//                 list, and add to the end of the OpList associated
//                 with the iterator.
// Special Notes : If we're in frequency domain, solution var access
//                 such as V(A) gets expanded into two ops, one for real part,
// Scope         : file-local
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void
createOps(const OutputMgr &output_manager, bool expandComplexTypes, const NetlistLocation &netlist_location, ParameterList::iterator begin, ParameterList::iterator end, std::back_insert_iterator<Util::OpList> inserter)
{
  Util::OpList tempOpList;

  makeOps(output_manager, netlist_location, begin, end, std::back_inserter(tempOpList));

  for (Util::OpList::const_iterator it = tempOpList.begin(); it != tempOpList.end(); ++it)
  {
    if (expandComplexTypes && (*it)->opType() == Util::SOLUTION_VAR)
    {
      std::string solutionName = (*it)->getName();
      int index = -1;
      const SolutionOp *op = dynamic_cast<const SolutionOp *>(*it);
      if (op)
        index = op->index_;

      delete *it;
      *inserter++ = new SolutionRealOp("Re(" + solutionName + ")", index);
      *inserter++ = new SolutionImaginaryOp("Im(" + solutionName + ")", index);
    }
    else if (expandComplexTypes && (*it)->opType() == Util::VOLTAGE_DIFFERENCE)
    {
      std::string solutionName = (*it)->getName();
      int index1 = -1;
      int index2 = -1;
      const VoltageDifferenceOp *op =
        dynamic_cast<const VoltageDifferenceOp *>(*it);
      if (op)
      {
        index1 = op->index1_;
        index2 = op->index2_;
      }

      delete *it;
      *inserter++ = new VoltageDifferenceRealOp("Re(" + solutionName + ")",
                                                index1, index2);
      *inserter++ = new VoltageDifferenceImaginaryOp("Im(" + solutionName + ")",
                                                     index1, index2);
    }
    else
      *inserter++ = *it;
  }
}

//-----------------------------------------------------------------------------
// Function      : fixupColumns
// Purpose       : given an output manager and a (print) parameter list, and
//                 an OpList, fill the list with the ops needed via createOps,
//                 and add additional output columns for things like
//                 index, time, frequency.  Set formatting options for these
//                 additional columns, set delimiter.
// Special Notes :
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
void fixupColumns(const OutputMgr &output_manager, PrintParameters &print_parameters, Util::OpList &op_list)
{
  createOps(output_manager, print_parameters.expandComplexTypes_, print_parameters.netlistLocation_, print_parameters.variableList_.begin(), print_parameters.variableList_.end(), std::back_inserter(op_list));

  Table::Justification justification = print_parameters.delimiter_.empty() ? Table::JUSTIFICATION_CENTER :  Table::JUSTIFICATION_NONE;

  for (Util::OpList::const_iterator it = op_list.begin() ; it != op_list.end(); ++it)
  {
    switch ((*it)->opType())
    {
      case Util::INDEX:
        print_parameters.table_.addColumn("INDEX", std::ios_base::fixed, 5, 0, Table::JUSTIFICATION_LEFT);
        break;

      case Util::TIME_VAR:
        print_parameters.table_.addColumn("TIME", print_parameters.streamWidth_, print_parameters.streamPrecision_, justification);
        break;

      case Util::FREQUENCY:
        print_parameters.table_.addColumn("FREQ", print_parameters.streamWidth_, print_parameters.streamPrecision_, justification);
        break;

      default:
        print_parameters.table_.addColumn((*it)->getName(), print_parameters.streamWidth_, print_parameters.streamPrecision_, justification);
        break;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : printValue
// Purpose       :
// Special Notes :
// Scope         : file-local
// Creator       : David Baur, Raytheon
// Creation Date :
//-----------------------------------------------------------------------------
std::ostream &printValue(std::ostream &os, const Table::Column &column, const std::string &delimiter, const int column_index, double value)
{
  if (delimiter.empty())
  {
    if (column_index != 0)
      os <<  " ";
    os << std::resetiosflags(std::ios_base::floatfield) << std::setiosflags(column.format_)
       << std::resetiosflags(std::ios_base::adjustfield) 
       << std::setiosflags(column.justification_ == Table::JUSTIFICATION_LEFT ? std::ios_base::left : std::ios_base::right)
       << std::setprecision(column.precision_) << std::setw(column.width_)
       << value;
  }
  else
  {
    if (column_index != 0)
      os << delimiter;
    os << std::resetiosflags(std::ios_base::floatfield) << std::setiosflags(column.format_) 
      << std::setw(0) << std::setprecision(column.precision_) << value;
  }

  return os;
}

} // namespace <unnamed>
} // namespace IO
} // namespace Xyce

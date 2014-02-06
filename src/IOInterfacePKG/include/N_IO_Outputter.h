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
// Filename       : $RCSfile: N_IO_Outputter.h,v $
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
// Revision Number: $Revision: 1.9.2.10 $
//
// Revision Date  : $Date: 2014/01/10 00:08:15 $
//
// Current Owner  : $Author: erkeite $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_Outputter_h
#define Xyce_N_IO_Outputter_h


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
#include <N_UTL_Demangle.h>

#include <N_IO_Outputter.h>
#include <N_IO_PkgOptionsMgr.h>
#include <N_ANP_SweepParam.h>
#include <N_IO_Objective.h>


// ---------- Trilinos Includes ----------
#include <Teuchos_SerialDenseMatrix.hpp>
#ifdef Xyce_Dakota
#include <Epetra_SerialDenseVector.h>
#endif

class N_ANP_AnalysisInterface;
class N_IO_CmdParse;
class N_IO_MeasureBase;
class N_IO_MeasureManager;
class N_IO_OutputFileBase;
class N_LAS_BlockVector;
class N_LAS_Matrix;
class N_LAS_Vector;
class N_MPDE_Manager;
class N_PDS_Comm;
class N_TOP_Topology;
class N_UTL_Expression;
class N_UTL_ExpressionData;

namespace Xyce {
namespace IO {

typedef std::list<N_UTL_Param> ParameterList;


//-----------------------------------------------------------------------------
// Class         : N_IO_Table
//
// Purpose       : This struct manages the header data
//
// Special Notes :
//
// Creator       : Todd Coffey, 1414
// Creation Date : 9/19/08
//-----------------------------------------------------------------------------
struct Table
{
    // Enum for column justification (left/center/right) in std header output
    enum Justification
    {
      JUSTIFICATION_LEFT,
      JUSTIFICATION_CENTER,
      JUSTIFICATION_RIGHT,
      JUSTIFICATION_NONE
    };

    struct Column
    {
        Column()
          : name_(),
            format_(ios_base::scientific),
            width_(17),
            precision_(9),
            justification_(JUSTIFICATION_LEFT)
        {}

        Column(const Column &column)
          : name_(column.name_),
            format_(column.format_),
            width_(column.width_),
            precision_(column.precision_),
            justification_(column.justification_)
        {}

        Column(std::string name, ios_base::fmtflags format, int width, int precision, Justification justification)
          : name_(name),
            format_(format),
            width_(width),
            precision_(precision),
            justification_(justification)
        {}

        std::string             name_;
        ios_base::fmtflags      format_;
        int                     width_;
        int                     precision_;
        Justification           justification_;
    };

    typedef std::vector<Column> ColumnList;

    Table()
    {}

    Table(const Table &table)
      : columnList_(table.columnList_.begin(), table.columnList_.end())
    {}

    Table &operator=(const Table &table) {
      columnList_.assign(table.columnList_.begin(), table.columnList_.end());

      return *this;
    }

    virtual ~Table()
    {}

    void addColumn(std::string name, ios_base::fmtflags format, int width, int precision, Justification justification)
    {
      columnList_.push_back(Column(name, format, width, precision, justification));
    }

    void addColumn(std::string name, int width, int precision, Justification justification)
    {
      columnList_.push_back(Column(name, ios_base::scientific, width, precision, justification));
    }

    ColumnList          columnList_;
};


namespace Format {

enum Format {STD, TECPLOT, PROBE, CSV, RAW, RAW_ASCII};

}

namespace PrintType {

enum PrintType {NONE, DC, TRAN, AC, AC_IC, HB, HB_IC, HB_STARTUP, HOMOTOPY, MPDE, MPDE_IC, RAW_OVERRIDE};

}

namespace OutputType {

enum OutputType {DC, TRAN, AC, AC_IC, HB_FD, HB_TD, HB_IC, HB_STARTUP, DCOP, HOMOTOPY, MPDE};

}

struct PrintParameters
{
    PrintParameters()
      : filename_(),
        suffix_(),
        extension_(),
        rawOverride_(false),
        printType_(PrintType::NONE),
        format_(Format::STD),
        index_(true),
        variableList_(),
        table_(),
        streamWidth_(17),
        streamPrecision_(9),
        timeWidth_(8),
        delimiter_(),
        outputTimeScaleFactor_(1.0)
    {}

    PrintParameters(const PrintParameters &print_parameters)
      : filename_(print_parameters.filename_),
        suffix_(print_parameters.suffix_),
        extension_(print_parameters.extension_),
        printType_(print_parameters.printType_),
        format_(print_parameters.format_),
        rawOverride_(print_parameters.rawOverride_),
        index_(print_parameters.index_),
        variableList_(print_parameters.variableList_.begin(), print_parameters.variableList_.end()),
        table_(print_parameters.table_),
        streamWidth_(print_parameters.streamWidth_),
        streamPrecision_(print_parameters.streamPrecision_),
        timeWidth_(print_parameters.timeWidth_),
        delimiter_(print_parameters.delimiter_),
        outputTimeScaleFactor_(print_parameters.outputTimeScaleFactor_)
    {}

    PrintParameters &operator=(const PrintParameters &print_parameters)
    {
      filename_ = print_parameters.filename_;
      suffix_ = print_parameters.suffix_;
      extension_ = print_parameters.extension_;
      printType_ = print_parameters.printType_;
      format_ = print_parameters.format_;
      index_ = print_parameters.index_;
      rawOverride_ = print_parameters.rawOverride_;
      variableList_.assign(print_parameters.variableList_.begin(), print_parameters.variableList_.end());
      table_ = print_parameters.table_;
      streamWidth_ = print_parameters.streamWidth_;
      streamPrecision_ = print_parameters.streamPrecision_;
      timeWidth_ = print_parameters.timeWidth_;
      delimiter_ = print_parameters.delimiter_;
      outputTimeScaleFactor_ = print_parameters.outputTimeScaleFactor_;

      return *this;
    }

    virtual ~PrintParameters()
    {}

    std::string                 filename_;
    std::string                 suffix_;
    std::string                 extension_;

    bool                        rawOverride_;
    bool                        index_;
    PrintType::PrintType        printType_;
    Format::Format              format_;
    ParameterList               variableList_;
    Table                       table_;
    int                         streamWidth_;
    int                         streamPrecision_;
    int                         timeWidth_;
    std::string                 delimiter_;

    double                      outputTimeScaleFactor_;         ///< output in something other than seconds (such as milli-seconds)
};

namespace Outputter {

class Interface
{
    const static int debug = false;

  public:
    void output(const N_LAS_Vector *solution_vector, const N_LAS_Vector *state_vector, const N_LAS_Vector *store_vector) {
      if (debug) std::cout << demangle(typeid(*this).name()) << " doOutput" << std::endl;

      doOutputTime(solution_vector, state_vector, store_vector);
    }

    void finishOutput() {
      if (debug) std::cout << demangle(typeid(*this).name()) << " doFinishOutput" << std::endl;

      doFinishOutput();
    }

    void resetOutput() {
      if (debug) std::cout << demangle(typeid(*this).name()) << " doResetOutput" << std::endl;

      doResetOutput();
    }

    void outputAC(double frequency, const N_LAS_Vector *real_solution_vector, const N_LAS_Vector *imaginary_solution_vector) {
      if (debug) std::cout << demangle(typeid(*this).name()) << " doOutputAC" << std::endl;

      doOutputFrequency(frequency, real_solution_vector, imaginary_solution_vector);
    }

    void outputHB (
      const std::vector<double>& timePoints,
      const std::vector<double>& freqPoints,
      const N_LAS_BlockVector & timeDomainSolnVec,
      const N_LAS_BlockVector & freqDomainSolnVecReal,
      const N_LAS_BlockVector & freqDomainSolnVecImaginary)
    {
      if (debug) std::cout << demangle(typeid(*this).name()) << " doOutputHB" << std::endl;

      doOutputHB(timePoints, freqPoints, timeDomainSolnVec, freqDomainSolnVecReal, freqDomainSolnVecImaginary);
    }

    void outputMPDE(double time, const N_LAS_Vector *solution_vector) {
      if (debug) std::cout << demangle(typeid(*this).name()) << " doOutputMPDE" << std::endl;

      doOutputMPDE(time, solution_vector);
    }

    void outputHomotopy(const std::vector<std::string> & parameter_names, const std::vector<double> & parameter_values, const N_LAS_Vector * solution_vector) {
      if (debug) std::cout << demangle(typeid(*this).name()) << " doOutputHomotopy" << std::endl;

      doOutputHomotopy(parameter_names, parameter_values, solution_vector);
    }

    void parse() {
      if (debug) std::cout << demangle(typeid(*this).name()) << " doParse" << std::endl;

      doParse();
    }

    void finishOutputStep() {
      if (debug) std::cout << demangle(typeid(*this).name()) << " doFinishOutputStep" << std::endl;

      doFinishOutputStep();
    }

    void outputMORTF(bool origSystem, const double & freq, const Teuchos::SerialDenseMatrix<int, std::complex<double> >& H) {
      if (debug) std::cout << demangle(typeid(*this).name()) << " doFinishOutputStep" << std::endl;

      doOutputMORTF(origSystem, freq, H);
    }

    virtual double getIndex() const = 0;

    virtual ~Interface()
    {}

  protected:
    virtual void doParse() = 0;
    virtual void doOpen() = 0;
    virtual void doOutputTime(const N_LAS_Vector *solution_vector, const N_LAS_Vector *state_vector, const N_LAS_Vector *store_vector) = 0;
    virtual void doResetOutput() = 0;
    virtual void doFinishOutput() = 0;
    virtual void doFinishOutputStep() = 0;

    virtual void doOutputFrequency(double frequency, const N_LAS_Vector *real_solution_vector, const N_LAS_Vector *imaginary_solution_vector) = 0;
    virtual void doOutputHB(const std::vector<double>& timePoints, const std::vector<double>& freqPoints,
                            const N_LAS_BlockVector & timeDomainSolnVec, const N_LAS_BlockVector & freqDomainSolnVecReal,
                            const N_LAS_BlockVector & freqDomainSolnVecImaginary) = 0;
    virtual void doOutputMPDE(double time, const N_LAS_Vector *solution_vector) = 0;
    virtual void doOutputHomotopy(const std::vector<std::string> & parameter_names, const std::vector<double> & param_values, const N_LAS_Vector * solution_vector) = 0;

    virtual void doOutputMORTF(bool origSystem, const double & freq, const Teuchos::SerialDenseMatrix<int, std::complex<double> >& H)= 0;
};

class TimeInterface : public Interface
{
  public:
    virtual void doResetOutput() {}

    virtual void doOutputFrequency(double frequency, const N_LAS_Vector *real_solution_vector, const N_LAS_Vector *imaginary_solution_vector) {}

    virtual void doOutputHB(const std::vector<double>& timePoints, const std::vector<double>& freqPoints,
                            const N_LAS_BlockVector & timeDomainSolnVec, const N_LAS_BlockVector & freqDomainSolnVecReal,
                            const N_LAS_BlockVector & freqDomainSolnVecImaginary) {}

    virtual void doOutputMPDE(double time, const N_LAS_Vector *solution_vector) {}

    virtual void doOutputHomotopy(const std::vector<std::string> & parameter_names, const std::vector<double> & param_values, const N_LAS_Vector * solution_vector) {}

    virtual void doOutputMORTF(bool origSystem, const double & freq, const Teuchos::SerialDenseMatrix<int, std::complex<double> >& H) {}
};

class FrequencyInterface : public Interface
{
  public:

    virtual void doResetOutput() {}

    virtual void doOutputTime(const N_LAS_Vector *solution_vector, const N_LAS_Vector *state_vector, const N_LAS_Vector *store_vector) {}

    virtual void doOutputHB(const std::vector<double>& timePoints, const std::vector<double>& freqPoints,
                            const N_LAS_BlockVector & timeDomainSolnVec, const N_LAS_BlockVector & freqDomainSolnVecReal,
                            const N_LAS_BlockVector & freqDomainSolnVecImaginary) {}

    virtual void doOutputMPDE(double time, const N_LAS_Vector *solution_vector) {}

    virtual void doOutputHomotopy(const std::vector<std::string> & parameter_names, const std::vector<double> & param_values, const N_LAS_Vector * solution_vector) {}

    virtual void doOutputMORTF(bool origSystem, const double & freq, const Teuchos::SerialDenseMatrix<int, std::complex<double> >& H) {}
};

class HBInterface : public Interface
{
  public:

    virtual void doResetOutput() {}

    virtual void doOutputFrequency(double frequency, const N_LAS_Vector *real_solution_vector, const N_LAS_Vector *imaginary_solution_vector) {}

    virtual void doOutputTime(const N_LAS_Vector *solution_vector, const N_LAS_Vector *state_vector, const N_LAS_Vector *store_vector) {}

    virtual void doOutputMPDE(double time, const N_LAS_Vector *solution_vector) {}

    virtual void doOutputHomotopy(const std::vector<std::string> & parameter_names, const std::vector<double> & param_values, const N_LAS_Vector * solution_vector) {}

    virtual void doOutputMORTF(bool origSystem, const double & freq, const Teuchos::SerialDenseMatrix<int, std::complex<double> >& H) {}
};

class TimePrn : public TimeInterface
{
  public:
    TimePrn(OutputMgr &output_manager, const PrintParameters &print_parameters);

    virtual ~TimePrn();

  private:
    TimePrn(const TimePrn &);
    TimePrn &operator=(const TimePrn &);

  public:
    void setOutputFilenameSuffix(const std::string &suffix) {
      suffix_ = suffix;
    }

    virtual double getIndex() const {
      return index_;
    }

    virtual void doParse();
    virtual void doOpen();
    virtual void doOutputTime(const N_LAS_Vector *solution_vector, const N_LAS_Vector *state_vector, const N_LAS_Vector *store_vector);
    virtual void doFinishOutput();
    virtual void doFinishOutputStep();

    void timeHeader();

  private:
    OutputMgr &         outputManager_;
    PrintParameters     printParameters_;
    bool                firstTimePrint_;
    int                 index_;
    std::string         outFilename_;
    std::string         suffix_;
    std::ostream *      outStreamPtr_;
    int                 headerPrintCalls_;
    int                 stepCount_;
};

class FrequencyPrn : public FrequencyInterface
{
  public:
    FrequencyPrn(OutputMgr &output_manager, const PrintParameters &print_parameters);

    virtual ~FrequencyPrn();

  private:
    FrequencyPrn(const FrequencyPrn &);
    FrequencyPrn &operator=(const FrequencyPrn &);

  public:
    void setOutputFilenameSuffix(const std::string &suffix) {
      suffix_ = suffix;
    }

    virtual double getIndex() const {
      return index_;
    }

    virtual void doParse();
    virtual void doOpen();
    virtual void doFinishOutput();
    virtual void doFinishOutputStep();

    virtual void doOutputFrequency(double frequency, const N_LAS_Vector *real_solution_vector, const N_LAS_Vector *imaginary_solution_vector);

  private:
    OutputMgr &         outputManager_;
    PrintParameters     printParameters_;
    bool                firstTimePrint_;
    int                 index_;
    std::string         outFilename_;
    std::string         suffix_;
    std::ostream *      outStreamPtr_;
    int                 stepCount_;
};

class TimeCSV : public TimeInterface
{
  public:
    TimeCSV(OutputMgr &output_manager, const PrintParameters &print_parameters);

    virtual ~TimeCSV();

  private:
    TimeCSV(const TimeCSV &);
    TimeCSV &operator=(const TimeCSV &);

  public:
    void setOutputFilenameSuffix(const std::string &suffix) {
      suffix_ = suffix;
    }

    virtual double getIndex() const {
      return index_;
    }

    virtual void doParse();
    virtual void doOpen();
    virtual void doOutputTime(const N_LAS_Vector *solution_vector, const N_LAS_Vector *state_vector, const N_LAS_Vector *store_vector);
    virtual void doFinishOutput();
    virtual void doFinishOutputStep();

    void timeHeader();

  private:
    OutputMgr &         outputManager_;
    PrintParameters     printParameters_;
    bool                firstTimePrint_;
    std::string         outFilename_;
    std::string         suffix_;
    std::ostream *      outStreamPtr_;
    int                 headerPrintCalls_;
    int                 index_;
    int                 stepCount_;
};

class FrequencyCSV : public FrequencyInterface
{
  public:
    FrequencyCSV(OutputMgr &output_manager, const PrintParameters &print_parameters);

    virtual ~FrequencyCSV();

  private:
    FrequencyCSV(const FrequencyCSV &);
    FrequencyCSV &operator=(const FrequencyCSV &);

  public:
    void setOutputFilenameSuffix(const std::string &suffix) {
      suffix_ = suffix;
    }

    virtual double getIndex() const {
      return index_;
    }

    virtual void doParse();
    virtual void doOpen();
    virtual void doFinishOutput();
    virtual void doFinishOutputStep();

    virtual void doOutputFrequency(double frequency, const N_LAS_Vector *real_solution_vector, const N_LAS_Vector *imaginary_solution_vector);

  private:
    OutputMgr &         outputManager_;
    PrintParameters     printParameters_;
    bool                firstTimePrint_;
    std::string         outFilename_;
    std::string         suffix_;
    std::ostream *      outStreamPtr_;
    int                 headerPrintCalls_;
    int                 index_;
    int                 stepCount_;
};

class TimeTecPlot : public TimeInterface
{
  public:
    TimeTecPlot(OutputMgr &output_manager, const PrintParameters &print_parameters);

    virtual ~TimeTecPlot();

  private:
    TimeTecPlot(const TimeTecPlot &);
    TimeTecPlot &operator=(const TimeTecPlot &);

  public:
    void probeHeader( std::ostream &ostreamPtr );
    void tecplotHeader( std::ostream &ostreamPtr );
    void stdHeader( std::ostream &ostreamPtr );
    void setOutputFilenameSuffix(const std::string &suffix) {
      suffix_ = suffix;
    }

    virtual double getIndex() const {
      return index_;
    }

    virtual void doParse();
    virtual void doOpen();
    virtual void doOutputTime(const N_LAS_Vector *solution_vector, const N_LAS_Vector *state_vector, const N_LAS_Vector *store_vector);
    virtual void doFinishOutput();
    virtual void doFinishOutputStep();

    void timeHeader();

  private:
    OutputMgr &         outputManager_;
    PrintParameters     printParameters_;
    bool                firstTimePrint_;
    std::string         outFilename_;
    std::string         suffix_;
    std::ostream *      outStreamPtr_;
    int                 headerPrintCalls_;
    int                 index_;
};

struct FrequencyTecPlot : public FrequencyInterface
{
    FrequencyTecPlot(OutputMgr &output_manager, const PrintParameters &print_parameters);

    virtual ~FrequencyTecPlot();

  private:
    FrequencyTecPlot(const FrequencyTecPlot &);
    FrequencyTecPlot &operator=(const FrequencyTecPlot &);

  public:
    virtual double getIndex() const {
      return index_;
    }

    virtual void doParse();
    virtual void doOpen();
    virtual void doFinishOutput();
    virtual void doFinishOutputStep();

    virtual void doOutputFrequency(double frequency, const N_LAS_Vector *real_solution_vector, const N_LAS_Vector *imaginary_solution_vector);

    void frequencyHeader();

  private:
    OutputMgr &         outputManager_;
    PrintParameters     printParameters_;
    std::string         outFilename_;
    std::ostream *      outStreamPtr_;
    int                 stepCount_;
    bool                firstTime_;
    unsigned long       index_;
};

struct OverrideRaw : public Interface
{
    OverrideRaw(OutputMgr &output_manager, const PrintParameters &print_parameters);

    virtual ~OverrideRaw();

  private:
    OverrideRaw(const OverrideRaw &);
    OverrideRaw &operator=(const OverrideRaw &);

  public:
    virtual double getIndex() const {
      return 0.0;
    }

    virtual void doParse();
    virtual void doOpen();
    virtual void doOutputTime(const N_LAS_Vector *solution_vector, const N_LAS_Vector *state_vector, const N_LAS_Vector *store_vector);
    virtual void doResetOutput();
    virtual void doFinishOutput();
    virtual void doFinishOutputStep();

    virtual void doOutputFrequency(double frequency, const N_LAS_Vector *real_solution_vector, const N_LAS_Vector *imaginary_solution_vector);
    virtual void doOutputHB(const std::vector<double>& timePoints, const std::vector<double>& freqPoints,
                            const N_LAS_BlockVector & timeDomainSolnVec, const N_LAS_BlockVector & freqDomainSolnVecReal,
                            const N_LAS_BlockVector & freqDomainSolnVecImaginary);
    virtual void doOutputMPDE(double time, const N_LAS_Vector *solution_vector);
    virtual void doOutputHomotopy(const std::vector<std::string> & parameter_names, const std::vector<double> & param_values, const N_LAS_Vector * solution_vector);
    virtual void doOutputMORTF(bool origSystem, const double & freq, const Teuchos::SerialDenseMatrix<int, std::complex<double> >& H);

    // raw headers output the data type: real or complex.
    // thus we need separate routines to handle if the the doOutputHeader call
    // originates form doOutput or doOutputAC
    void timeHeader();
    void frequencyHeader();

  private:
    OutputMgr &         outputManager_;
    PrintParameters     printParameters_;
    std::string         outFilename_;
    int                 numPoints_;
    long                numPointsPos_;
    std::ostream *      outStreamPtr_;
    bool                outputRAWTitleAndDate_;
};

struct TimeRaw : public TimeInterface
{
    TimeRaw(OutputMgr &output_manager, const PrintParameters &print_parameters);

    virtual ~TimeRaw();

  private:
    TimeRaw(const TimeRaw &);
    TimeRaw &operator=(const TimeRaw &);

  public:
    virtual double getIndex() const {
      return 0.0;
    }

    virtual void doParse();
    virtual void doOpen();
    virtual void doOutputTime(const N_LAS_Vector *solution_vector, const N_LAS_Vector *state_vector, const N_LAS_Vector *store_vector);
    virtual void doFinishOutput();
    virtual void doFinishOutputStep();

    // raw headers output the data type: real or complex.
    // thus we need separate routines to handle if the the doOutputHeader call
    // originates form doOutput or doOutputAC
    void timeHeader();

  private:
    OutputMgr &         outputManager_;
    PrintParameters     printParameters_;
    std::string         outFilename_;
    int                 numPoints_;
    long                numPointsPos_;
    std::ostream *      outStreamPtr_;
    bool                outputRAWTitleAndDate_;
};

struct FrequencyRaw : public FrequencyInterface
{
    FrequencyRaw(OutputMgr &output_manager, const PrintParameters &print_parameters);

    virtual ~FrequencyRaw();

  private:
    FrequencyRaw(const FrequencyRaw &);
    FrequencyRaw &operator=(const FrequencyRaw &);

  public:
    virtual double getIndex() const {
      return 0.0;
    }

    virtual void doParse();
    virtual void doOpen();
    virtual void doFinishOutput();
    virtual void doFinishOutputStep();

    virtual void doOutputFrequency(double frequency, const N_LAS_Vector *real_solution_vector, const N_LAS_Vector *imaginary_solution_vector);

    // raw headers output the data type: real or complex.
    // thus we need separate routines to handle if the the doOutputHeader call
    // originates form doOutput or doOutputAC
    void frequencyHeader();

  private:
    OutputMgr &         outputManager_;
    PrintParameters     printParameters_;
    std::string         outFilename_;
    int                 numPoints_;
    long                numPointsPos_;
    std::ostream *      outStreamPtr_;
    bool                outputRAWTitleAndDate_;
};

struct TimeRawAscii : public TimeInterface
{
    TimeRawAscii(OutputMgr &output_manager, const PrintParameters &print_parameters);

    virtual ~TimeRawAscii();

  private:
    TimeRawAscii(const TimeRawAscii &);
    TimeRawAscii &operator=(const TimeRawAscii &);

  public:
    virtual double getIndex() const {
      return 0.0;
    }

    virtual void doParse();
    virtual void doOpen();
    virtual void doOutputTime(const N_LAS_Vector *solution_vector, const N_LAS_Vector *state_vector, const N_LAS_Vector *store_vector);
    virtual void doFinishOutput();
    virtual void doFinishOutputStep();

    // raw headers output the data type: real or complex.
    // thus we need separate routines to handle if the the doOutputHeader call
    // originates form doOutput or doOutputAC
    void timeHeader();

  private:
    OutputMgr &         outputManager_;
    PrintParameters     printParameters_;
    std::string         outFilename_;
    int                 numPoints_;
    long                numPointsPos_;
    std::ostream *      outStreamPtr_;
    bool                printAll_;
    bool                outputRAWTitleAndDate_;
};

struct FrequencyRawAscii : public FrequencyInterface
{
    FrequencyRawAscii(OutputMgr &output_manager, const PrintParameters &print_parameters);

    virtual ~FrequencyRawAscii();

  private:
    FrequencyRawAscii(const FrequencyRawAscii &);
    FrequencyRawAscii &operator=(const FrequencyRawAscii &);

  public:
    virtual double getIndex() const {
      return 0.0;
    }

    virtual void doParse();
    virtual void doOpen();
    virtual void doFinishOutput();
    virtual void doFinishOutputStep();

    virtual void doOutputFrequency(double frequency, const N_LAS_Vector *real_solution_vector, const N_LAS_Vector *imaginary_solution_vector);

    // raw headers output the data type: real or complex.
    // thus we need separate routines to handle if the the doOutputHeader call
    // originates form doOutput or doOutputAC
    void frequencyHeader();


  private:
    OutputMgr &         outputManager_;
    PrintParameters     printParameters_;
    std::string         outFilename_;
    int                 numPoints_;
    long                numPointsPos_;
    std::ostream *      outStreamPtr_;
    bool                printAll_;
    bool                outputRAWTitleAndDate_;
};

struct OverrideRawAscii : public Interface
{
    OverrideRawAscii(OutputMgr &output_manager, const PrintParameters &print_parameters);

    virtual ~OverrideRawAscii();

  private:
    OverrideRawAscii(const OverrideRawAscii &);
    OverrideRawAscii &operator=(const OverrideRawAscii &);

  public:
    virtual double getIndex() const {
      return 0.0;
    }

    virtual void doParse();
    virtual void doOpen();
    virtual void doOutputTime(const N_LAS_Vector *solution_vector, const N_LAS_Vector *state_vector, const N_LAS_Vector *store_vector);
    virtual void doResetOutput();
    virtual void doFinishOutput();
    virtual void doFinishOutputStep();

    virtual void doOutputFrequency(double frequency, const N_LAS_Vector *real_solution_vector, const N_LAS_Vector *imaginary_solution_vector);
    virtual void doOutputHB(const std::vector<double>& timePoints, const std::vector<double>& freqPoints,
                            const N_LAS_BlockVector & timeDomainSolnVec, const N_LAS_BlockVector & freqDomainSolnVecReal,
                            const N_LAS_BlockVector & freqDomainSolnVecImaginary);
    virtual void doOutputMPDE(double time, const N_LAS_Vector *solution_vector);
    virtual void doOutputHomotopy(const std::vector<std::string> & parameter_names, const std::vector<double> & param_values, const N_LAS_Vector * solution_vector);
    virtual void doOutputMORTF(bool origSystem, const double & freq, const Teuchos::SerialDenseMatrix<int, std::complex<double> >& H);

    // raw headers output the data type: real or complex.
    // thus we need separate routines to handle if the the doOutputHeader call
    // originates form doOutput or doOutputAC
    void timeHeader();
    void frequencyHeader();


  private:
    OutputMgr &         outputManager_;
    PrintParameters     printParameters_;
    std::string         outFilename_;
    int                 numPoints_;
    long                numPointsPos_;
    std::ostream *      outStreamPtr_;
    bool                printAll_;
    bool                outputRAWTitleAndDate_;
};


struct TimeProbe : public TimeInterface
{
    TimeProbe(OutputMgr &output_manager, const PrintParameters &print_parameters);

    virtual ~TimeProbe();

  private:
    TimeProbe(const TimeProbe &);
    TimeProbe &operator=(const TimeProbe &);

  public:
    virtual double getIndex() const {
      return index_;
    }

    virtual void doParse();
    virtual void doOpen();
    virtual void doOutputTime(const N_LAS_Vector *solution_vector, const N_LAS_Vector *state_vector, const N_LAS_Vector *store_vector);
    virtual void doFinishOutput();
    virtual void doFinishOutputStep();


    void timeHeader();

  private:
    OutputMgr &         outputManager_;
    PrintParameters     printParameters_;
    bool                firstTimePrint_;
    int                 printCount_;
    std::string         outFilename_;
    std::string         suffix_;
    std::ostream *      outStreamPtr_;
    int                 headerPrintCalls_;
    int                 index_;
};

struct FrequencyProbe : public FrequencyInterface
{
    FrequencyProbe(OutputMgr &output_manager, const PrintParameters &print_parameters);

    virtual ~FrequencyProbe();

  private:
    FrequencyProbe(const FrequencyProbe &);
    FrequencyProbe &operator=(const FrequencyProbe &);

  public:
    virtual double getIndex() const {
      return index_;
    }

    virtual void doParse();
    virtual void doOpen();
    virtual void doFinishOutput();
    virtual void doFinishOutputStep();

    virtual void doOutputFrequency(double frequency, const N_LAS_Vector *real_solution_vector, const N_LAS_Vector *imaginary_solution_vector);

    void frequencyHeader();

  private:
    OutputMgr &         outputManager_;
    PrintParameters     printParameters_;
    bool                firstTimePrint_;
    int                 printCount_;
    std::string         outFilename_;
    std::string         suffix_;
    std::ostream *      outStreamPtr_;
    int                 headerPrintCalls_;
    int                 index_;
};

struct HBPrn : public HBInterface
{
    HBPrn(OutputMgr &output_manager, const PrintParameters &freq_print_parameters, const PrintParameters &time_print_parameters);

    virtual ~HBPrn();

  private:
    HBPrn(const HBPrn &);
    HBPrn &operator=(const HBPrn &);

  public:

    virtual double getIndex() const {
      return index_;
    }

    virtual void doParse();
    virtual void doOpen();
    virtual void doFinishOutput();
    virtual void doFinishOutputStep();
    virtual void doOutputHB(const std::vector<double>& timePoints, const std::vector<double>& freqPoints,
                            const N_LAS_BlockVector & timeDomainSolnVec, const N_LAS_BlockVector & freqDomainSolnVecReal,
                            const N_LAS_BlockVector & freqDomainSolnVecImaginary);

    void doOutputHeader();

  private:
    OutputMgr &         outputManager_;
    PrintParameters     freqPrintParameters_;
    PrintParameters     timePrintParameters_;
    int                 stepCount_;
    int                 index_;
    bool                firstTimeHB_;
    std::string         timeFilename_;
    std::string         freqFilename_;
    std::ostream *      timeStreamPtr_;
    std::ostream *      freqStreamPtr_;
};

struct HBCSV : public HBInterface
{
    HBCSV(OutputMgr &output_manager, const PrintParameters &freq_print_parameters, const PrintParameters &time_print_parameters);

    virtual ~HBCSV();

  private:
    HBCSV(const HBCSV &);
    HBCSV &operator=(const HBCSV &);

  public:

    virtual double getIndex() const {
      return index_;
    }

    virtual void doParse();
    virtual void doOpen();
    virtual void doFinishOutput();
    virtual void doFinishOutputStep();
    virtual void doOutputHB(const std::vector<double>& timePoints, const std::vector<double>& freqPoints,
                            const N_LAS_BlockVector & timeDomainSolnVec, const N_LAS_BlockVector & freqDomainSolnVecReal,
                            const N_LAS_BlockVector & freqDomainSolnVecImaginary);

  private:
    OutputMgr &         outputManager_;
    PrintParameters     freqPrintParameters_;
    PrintParameters     timePrintParameters_;
    int                 stepCount_;
    int                 index_;
    bool                firstTimeHB_;
    std::string         timeFilename_;
    std::string         freqFilename_;
    std::ostream *      timeStreamPtr_;
    std::ostream *      freqStreamPtr_;
};

struct HBTecPlot : public HBInterface
{
    HBTecPlot(OutputMgr &output_manager, const PrintParameters &freq_print_parameters, const PrintParameters &time_print_parameters);

    virtual ~HBTecPlot();

  private:
    HBTecPlot(const HBTecPlot &);
    HBTecPlot &operator=(const HBTecPlot &);

  private:
    void tecplotTimeHBHeader( ostream & stream);

  public:
    virtual double getIndex() const {
      return 0.0;
    }

    virtual void doParse();
    virtual void doOpen();
    virtual void doFinishOutput();
    virtual void doFinishOutputStep();

    virtual void doOutputHB(const std::vector<double>& timePoints, const std::vector<double>& freqPoints,
                            const N_LAS_BlockVector & timeDomainSolnVec, const N_LAS_BlockVector & freqDomainSolnVecReal,
                            const N_LAS_BlockVector & freqDomainSolnVecImaginary);

  private:
    OutputMgr &         outputManager_;
    PrintParameters     freqPrintParameters_;
    PrintParameters     timePrintParameters_;
    int                 stepCount_;
    int                 index_;
    bool                firstTimeHB_;
    std::string         timeFilename_;
    std::string         freqFilename_;
    std::ostream *      timeStreamPtr_;
    std::ostream *      freqStreamPtr_;
};

struct MPDEPrn : public Interface
{
    MPDEPrn(OutputMgr &output_manager, const PrintParameters &print_parameters);

    virtual ~MPDEPrn();

  private:
    MPDEPrn(const MPDEPrn &);
    MPDEPrn &operator=(const MPDEPrn &);

  public:

    void stdTimeMPDEHeader(std::ostream & stream );

    virtual double getIndex() const {
      return 0.0;
    }

    virtual void doParse();
    virtual void doOpen();
    virtual void doOutputTime(const N_LAS_Vector *solution_vector, const N_LAS_Vector *state_vector, const N_LAS_Vector *store_vector);
    virtual void doResetOutput();
    virtual void doFinishOutput();
    virtual void doFinishOutputStep();

    virtual void doOutputFrequency(double frequency, const N_LAS_Vector *real_solution_vector, const N_LAS_Vector *imaginary_solution_vector);
    virtual void doOutputHB(const std::vector<double>& timePoints, const std::vector<double>& freqPoints,
                            const N_LAS_BlockVector & timeDomainSolnVec, const N_LAS_BlockVector & freqDomainSolnVecReal,
                            const N_LAS_BlockVector & freqDomainSolnVecImaginary);
    virtual void doOutputMPDE(double time, const N_LAS_Vector *solution_vector);
    virtual void doOutputHomotopy(const std::vector<std::string> & parameter_names, const std::vector<double> & param_values, const N_LAS_Vector * solution_vector);
    virtual void doOutputMORTF(bool origSystem, const double & freq, const Teuchos::SerialDenseMatrix<int, std::complex<double> >& H);

    void mpdeHeader();

  private:
    OutputMgr &         outputManager_;
    PrintParameters     printParameters_;
    std::string         outFilename_;
    std::ostream *      outStreamPtr_;
    int                 stepCount_;
    bool                firstTimeMPDE_;
    int                 n1_;
    int                 n2_;
};


struct MPDETecPlot : public Interface
{
    MPDETecPlot(OutputMgr &output_manager, const PrintParameters &print_parameters);

    virtual ~MPDETecPlot();

  private:
    MPDETecPlot(const MPDETecPlot &);
    MPDETecPlot &operator=(const MPDETecPlot &);

  public:

    void stdTimeMPDEHeader(std::ostream & stream );

    virtual double getIndex() const {
      return 0.0;
    }

    virtual void doParse();
    virtual void doOpen();
    virtual void doOutputTime(const N_LAS_Vector *solution_vector, const N_LAS_Vector *state_vector, const N_LAS_Vector *store_vector);
    virtual void doResetOutput();
    virtual void doFinishOutput();
    virtual void doFinishOutputStep();

    virtual void doOutputFrequency(double frequency, const N_LAS_Vector *real_solution_vector, const N_LAS_Vector *imaginary_solution_vector);
    virtual void doOutputHB(const std::vector<double>& timePoints, const std::vector<double>& freqPoints,
                            const N_LAS_BlockVector & timeDomainSolnVec, const N_LAS_BlockVector & freqDomainSolnVecReal,
                            const N_LAS_BlockVector & freqDomainSolnVecImaginary);
    virtual void doOutputMPDE(double time, const N_LAS_Vector *solution_vector);
    virtual void doOutputHomotopy(const std::vector<std::string> & parameter_names, const std::vector<double> & param_values, const N_LAS_Vector * solution_vector);
    virtual void doOutputMORTF(bool origSystem, const double & freq, const Teuchos::SerialDenseMatrix<int, std::complex<double> >& H);

    void doOutputHeader();

  private:
    OutputMgr &         outputManager_;
    PrintParameters     printParameters_;
    std::string         outFilename_;
    std::ostream *      outStreamPtr_;
    int                 stepCount_;
    bool                firstTimeMPDE_;
    int                 n1_;
    int                 n2_;
};


struct HomotopyPrn : public Interface
{
    HomotopyPrn(OutputMgr &output_manager, const PrintParameters &print_parameters);

    virtual ~HomotopyPrn();

  private:
    HomotopyPrn(const HomotopyPrn &);
    HomotopyPrn &operator=(const HomotopyPrn &);

  public:

    virtual double getIndex() const {
      return index_;
    }

    virtual void doParse();
    virtual void doOpen();
    virtual void doResetOutput();
    virtual void doFinishOutput();
    virtual void doFinishOutputStep();

    virtual void doOutputTime(const N_LAS_Vector *solution_vector, const N_LAS_Vector *state_vector, const N_LAS_Vector *store_vector) {}
    virtual void doOutputFrequency(double frequency, const N_LAS_Vector *real_solution_vector, const N_LAS_Vector *imaginary_solution_vector) {}
    virtual void doOutputHB(const std::vector<double>& timePoints, const std::vector<double>& freqPoints,
                            const N_LAS_BlockVector & timeDomainSolnVec, const N_LAS_BlockVector & freqDomainSolnVecReal,
                            const N_LAS_BlockVector & freqDomainSolnVecImaginary) {}
    virtual void doOutputMPDE(double time, const N_LAS_Vector *solution_vector) {}
    virtual void doOutputMORTF(bool origSystem, const double & freq, const Teuchos::SerialDenseMatrix<int, std::complex<double> >& H) {}

    virtual void doOutputHomotopy(const std::vector<std::string> & parameter_names, const std::vector<double> & param_values, const N_LAS_Vector * solution_vector);

    void homotopyHeader(const std::vector<std::string> & parameter_names, const std::vector<double> & param_values, const N_LAS_Vector * solution_vector);

  private:
    OutputMgr &         outputManager_;
    PrintParameters     printParameters_;
    std::string         outFilename_;
    std::ostream *      outStreamPtr_;
    int                 stepCount_;
    unsigned long       index_;
    int                 printCount_;
    bool                firstTimeHomotopy_;
    Table::ColumnList   columnList_;
};

struct HomotopyTecPlot : public Interface
{
    HomotopyTecPlot(OutputMgr &output_manager, const PrintParameters &print_parameters);

    virtual ~HomotopyTecPlot();

  private:
    HomotopyTecPlot(const HomotopyTecPlot &);
    HomotopyTecPlot &operator=(const HomotopyTecPlot &);

  public:

    virtual double getIndex() const {
      return 0.0;
    }

    virtual void doParse();
    virtual void doOpen();
    virtual void doOutputTime(const N_LAS_Vector *solution_vector, const N_LAS_Vector *state_vector, const N_LAS_Vector *store_vector);
    virtual void doResetOutput();
    virtual void doFinishOutput();
    virtual void doFinishOutputStep();

    virtual void doOutputFrequency(double frequency, const N_LAS_Vector *real_solution_vector, const N_LAS_Vector *imaginary_solution_vector);
    virtual void doOutputHB(const std::vector<double>& timePoints, const std::vector<double>& freqPoints,
                            const N_LAS_BlockVector & timeDomainSolnVec, const N_LAS_BlockVector & freqDomainSolnVecReal,
                            const N_LAS_BlockVector & freqDomainSolnVecImaginary);
    virtual void doOutputMPDE(double time, const N_LAS_Vector *solution_vector);
    virtual void doOutputHomotopy(const std::vector<std::string> & parameter_names, const std::vector<double> & param_values, const N_LAS_Vector * solution_vector);
    virtual void doOutputMORTF(bool origSystem, const double & freq, const Teuchos::SerialDenseMatrix<int, std::complex<double> >& H);

    void doOutputHeader(const std::vector<std::string> & parameter_names, const std::vector<double> & param_values, const N_LAS_Vector * solution_vector);

  private:
    OutputMgr &         outputManager_;
    PrintParameters     printParameters_;
    std::string         outFilename_;
    std::ostream *      outStreamPtr_;
    int                 stepCount_;
    unsigned long       index_;
    int                 printCount_;
    bool                firstTimeHomotopy_;
};

struct HomotopyProbe : public Interface
{
    HomotopyProbe(OutputMgr &output_manager, const PrintParameters &print_parameters);

    virtual ~HomotopyProbe();

  private:
    HomotopyProbe(const HomotopyProbe &);
    HomotopyProbe &operator=(const HomotopyProbe &);

  public:

    virtual double getIndex() const {
      return 0.0;
    }

    virtual void doParse();
    virtual void doOpen();
    virtual void doOutputTime(const N_LAS_Vector *solution_vector, const N_LAS_Vector *state_vector, const N_LAS_Vector *store_vector);
    virtual void doResetOutput();
    virtual void doFinishOutput();
    virtual void doFinishOutputStep();

    virtual void doOutputFrequency(double frequency, const N_LAS_Vector *real_solution_vector, const N_LAS_Vector *imaginary_solution_vector);
    virtual void doOutputHB(const std::vector<double>& timePoints, const std::vector<double>& freqPoints,
                            const N_LAS_BlockVector & timeDomainSolnVec, const N_LAS_BlockVector & freqDomainSolnVecReal,
                            const N_LAS_BlockVector & freqDomainSolnVecImaginary);
    virtual void doOutputMPDE(double time, const N_LAS_Vector *solution_vector);
    virtual void doOutputHomotopy(const std::vector<std::string> & parameter_names, const std::vector<double> & param_values, const N_LAS_Vector * solution_vector);
    virtual void doOutputMORTF(bool origSystem, const double & freq, const Teuchos::SerialDenseMatrix<int, std::complex<double> >& H);

    void doOutputHeader(const std::vector<std::string> & parameter_names, const std::vector<double> & param_values, const N_LAS_Vector * solution_vector);

  private:
    OutputMgr &         outputManager_;
    PrintParameters     printParameters_;
    std::string         outFilename_;
    std::ostream *      outStreamPtr_;
    int                 stepCount_;
    unsigned long       index_;
    int                 printCount_;
    bool                firstTimeHomotopy_;
};

class MOR : public Interface
{
  public:
    MOR(OutputMgr &output_manager, const PrintParameters &print_parameters);

    virtual ~MOR();

  private:
    MOR(const MOR &);
    MOR &operator=(const MOR &);

  public:
    void setOutputFilenameSuffix(const std::string &suffix) {
      suffix_ = suffix;
    }

    virtual double getIndex() const {
      return index_;
    }

    virtual void doParse();
    virtual void doOpen();
    virtual void doOutputTime(const N_LAS_Vector *solution_vector, const N_LAS_Vector *state_vector, const N_LAS_Vector *store_vector);
    virtual void doResetOutput();
    virtual void doFinishOutput();
    virtual void doFinishOutputStep();

    virtual void doOutputFrequency(double frequency, const N_LAS_Vector *real_solution_vector, const N_LAS_Vector *imaginary_solution_vector);
    virtual void doOutputHB(const std::vector<double>& timePoints, const std::vector<double>& freqPoints,
                            const N_LAS_BlockVector & timeDomainSolnVec, const N_LAS_BlockVector & freqDomainSolnVecReal,
                            const N_LAS_BlockVector & freqDomainSolnVecImaginary);
    virtual void doOutputMPDE(double time, const N_LAS_Vector *solution_vector);
    virtual void doOutputHomotopy(const std::vector<std::string> & parameter_names, const std::vector<double> & param_values, const N_LAS_Vector * solution_vector);
    virtual void doOutputMORTF(bool origSystem, const double & freq, const Teuchos::SerialDenseMatrix<int, std::complex<double> >& H);

    void openMORFiles(bool openOrig);
    void outputMORHeaders(int numPorts);
    void tecplotFreqMORHeader_(ostream & stream, int counter, int numPorts);
    void stdFreqMORHeader_(ostream & stream, int numPorts);

  private:
    OutputMgr &         outputManager_;
    PrintParameters     printParameters_;
    bool                firstTimePrint_;
    int                 index_;
    Format::Format      format_;
    std::string         outFilename_;
    std::string         suffix_;
    std::ostream *      outStreamPtr_;
    int                 headerPrintCalls_;
    int                 stepCount_;
    bool                openOrig_;
};

} // namespace outputter

} // namespace IO
} // namespace Xyce

#endif // Xyce_N_IO_Outputter_h

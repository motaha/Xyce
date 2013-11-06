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
// Filename       : $RCSfile: N_IO_OutputFileBase.h,v $
//
// Purpose        : Base class for handling file output of simulation results.
//
// Special Notes  :
//
// Creator        : Richard Schiek, Electrical Systems Modeling, Sandia National Laboratories
//
// Creation Date  : 12/05/12
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.2.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:42 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_OutputFileBase_h
#define Xyce_N_IO_OutputFileBase_h

// ----------   Standard Includes   ----------
#include <string>
#include <list>
#include <vector>
#include <ostream>
#include <fstream>

// ----------   Xyce Includes   ----------
#include <N_UTL_Xyce.h>
#include <N_UTL_Misc.h>
#include <N_UTL_Param.h>
#include <N_UTL_OptionBlock.h>
#include <N_LAS_Vector.h>
#include <N_LAS_BlockVector.h>
#include <N_ANP_SweepParam.h>
#include <N_IO_OutputMgr.h>

// ---------- Forward Declarations ----------

//-------------------------------------------------------------------------
// Class         : N_IO_OutputFileBase
// Purpose       : Base class for output of simulation results to a file.
// Special Notes : 
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//-------------------------------------------------------------------------
class N_IO_OutputFileBase
{
 
public:
  N_IO_OutputFileBase();
  ~N_IO_OutputFileBase();
  
  // these functions have basic implementations that 
  // should work for most applications.  Thus they don't always 
  // need an override.
  virtual void openFile( string basename, string simulationSuffix="");
  void setFileSuffix( string newSuffix ) {fileSuffix_ = newSuffix;};
  virtual void closeOutput();
  
  // These functions will depend on the output format.  Thus, 
  // they emit errors if the base clase version is called.
  virtual void outputHeader();
  
  void outputHOMOTOPYHeader( const vector<string> & paramNames);
  
  virtual void outputDC( 
      const int dcNumber, 
      const int maxDC, 
      const vector<N_ANP_SweepParam> & dcParamVec1,
      N_LAS_Vector * solnVecPtr,
      N_LAS_Vector * stateVecPtr,
      N_LAS_Vector * storeVecPtr );
      
  virtual void outputTran(
      const double & time, 
      N_LAS_Vector * solnVecPtr,
      N_LAS_Vector * stateVecPtr,
      N_LAS_Vector * storeVecPtr );
    
  virtual void outputStep(
      const int stepNumber, 
      const int maxStep, 
      const vector<N_ANP_SweepParam> & stepParamVec1,
      N_LAS_Vector * solnVecPtr,
      N_LAS_Vector * stateVecPtr,
      N_LAS_Vector * storeVecPtr );
      
  virtual void outputAC( 
      const double & freq,
      N_LAS_Vector * freqDomainSolnVecReal,
      N_LAS_Vector * freqDomainSolnVecImaginary);
      
  void outputHOMOTOPY( const vector<string> & paramNames,
      const vector<double> & paramVals,
      N_LAS_Vector * solnVecPtr );
  
  virtual void outputMPDE(const double & time, N_LAS_Vector * solnVecPtr );
  
  virtual void outputHB(
      const N_LAS_BlockVector & timeDomainSolnVec, 
      const N_LAS_BlockVector & freqDomainSolnVecReal,
      const N_LAS_BlockVector & freqDomainSolnVecImaginary);
  
  virtual void outputMOR();
  
  virtual void finishOutput();
  
  // these functions are intended to let Xyce re-read a simulation output 
  // file and then recalculate output metrics in .measure() statements 
  // without re-running the original simulation. 
  virtual bool openFileForRead( string filename );
  virtual bool getOutputVarNames( vector< string > & varNames );
  virtual bool getOutputNextVarValues( N_LAS_Vector * varValues );
  virtual bool closeFileForRead();
  
  void convertOutputNamesToSolVarNames( vector< string > & varNames);
  
protected:  
  ostream * ostreamPtr_;           // the output stream to use.
  string outputFileBaseName_;      // the string name of the output file 
  string fileSuffix_;              // filename suffix if needed.
  string simulationSuffix_;        // any string that should be added as a suffix on the filename before the fileSuffix_ 
                                      
  string fullFileName_;            // should be  outputFileBaseName_ + simulationSuffixNumber_ + fileSuffix_
  
  string fileFormatName_;       // This is the name of the derived class file format. 
                                // Any virtual functions not overridden by the derived 
                                // class can make an error citing what wasn't supported.
                                // as in "Probe formated output does not currently support MPDE output."
  
  bool appendOutputFlag_;       // default to false, but if true, then append new data to output file.
  
  // for re-reading existing files
  ifstream * istreamPtr_;
  string inputFileBaseName_;
};

#endif // Xyce_N_IO_OutputFileBase_h

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
// Filename       : $RCSfile: N_IO_OutputPrn.h,v $
//
// Purpose        : Base class for handling file output of simulation results.
//
// Special Notes  :
//
// Creator        : Richard Schiek, Electrical Systems Modeling, Sandia National Laboratories
//
// Creation Date  : 12/06/12
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.7 $
//
// Revision Date  : $Date: 2014/02/24 23:49:20 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_OutputPrn_h
#define Xyce_N_IO_OutputPrn_h

// ----------   Standard Includes   ----------

// ----------   Xyce Includes   ----------
#include<N_IO_OutputFileBase.h>

// ---------- Forward Declarations ----------

class N_IO_OutputPrn : public N_IO_OutputFileBase
{
  public:

    N_IO_OutputPrn();
    ~N_IO_OutputPrn();

    // these functions are intended to let Xyce re-read a simulation output
    // file and then recalculate output metrics in .measure() statements
    // without re-running the original simulation.
    bool getOutputVarNames( std::vector< std::string > & varNames );
    bool getOutputNextVarValues( N_LAS_Vector * varValues );


// These functions will depend on the output format.  Thus,
// they emit errors if the base clase version is called.
    virtual void outputHeader()
    {}

    virtual void outputDC(
      const int dcNumber,
      const int maxDC,
      const std::vector<N_ANP_SweepParam> & dcParamVec1,
      N_LAS_Vector * solnVecPtr,
      N_LAS_Vector * stateVecPtr,
      N_LAS_Vector * storeVecPtr )
    {}

    virtual void outputTran(
      const double & time,
      N_LAS_Vector * solnVecPtr,
      N_LAS_Vector * stateVecPtr,
      N_LAS_Vector * storeVecPtr )
    {}

    virtual void outputStep(
      const int stepNumber,
      const int maxStep,
      const std::vector<N_ANP_SweepParam> & stepParamVec1,
      N_LAS_Vector * solnVecPtr,
      N_LAS_Vector * stateVecPtr,
      N_LAS_Vector * storeVecPtr )
    {}

    virtual void outputAC(
      const double & freq,
      N_LAS_Vector * freqDomainSolnVecReal,
      N_LAS_Vector * freqDomainSolnVecImaginary)
    {}

    virtual void outputMPDE(const double & time, N_LAS_Vector * solnVecPtr )
    {}

    virtual void outputHB(
      const N_LAS_BlockVector & timeDomainSolnVec,
      const N_LAS_BlockVector & freqDomainSolnVecReal,
      const N_LAS_BlockVector & freqDomainSolnVecImaginary)
    {}

    virtual void outputMOR()
    {}

    virtual void finishOutput()
    {}
};

#endif // Xyce_N_IO_OutputPrn_h

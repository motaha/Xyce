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
// Filename       : $RCSfile: N_ANP_OutputMgrAdapter.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date  : 01/24/08
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.26.2.3 $
//
// Revision Date  : $Date: 2013/10/03 17:23:30 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_ANP_OutputMgrAdapter_h
#define Xyce_N_ANP_OutputMgrAdapter_h

#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RefCountPtr;
using Teuchos::rcp;
#include <Teuchos_SerialDenseMatrix.hpp>

// ----------   Xyce Includes   ----------
#include <N_UTL_Xyce.h>
#include <N_UTL_NoCase.h>
#include <N_IO_OutputMgr.h>

// ---------- Forward Declarations ----------
class N_ANP_SweepParam;
class N_LAS_Vector;
class N_PDS_Comm;

//-------------------------------------------------------------------------
// Class         : N_ANP_OutputMgrAdapter
// Purpose       : Inteface class for the output manager
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//-------------------------------------------------------------------------
class N_ANP_OutputMgrAdapter
{

  public:
    N_ANP_OutputMgrAdapter( );

    virtual ~N_ANP_OutputMgrAdapter() {}

    void registerOutputMgr( N_IO_OutputMgr * outputMgrPtr )
    {
      outputMgrRCPtr_ = rcp( outputMgrPtr, false );
    }

    void setStepParamVec( const RefCountPtr< vector<N_ANP_SweepParam> > & paramVec )
    {
      stepParamVecRCPtr_ = paramVec;
    }

    void setDCParamVec( const RefCountPtr< vector<N_ANP_SweepParam> > & paramVec )
    {
      dcParamVecRCPtr_ = paramVec;
    }

    // accessor methods
    int getStepAnalysisStepNumber()
    { return stepAnalysisStepNumber_; }

    int getStepAnalysisMaxSteps()
    { return stepAnalysisMaxSteps_; }

    int getDCAnalysisStepNumber()
    { return dcAnalysisStepNumber_; }

    int getDCAnalysisMaxSteps()
    { return dcAnalysisMaxSteps_; }

#ifdef Xyce_PARALLEL_MPI
    N_PDS_Comm * getCommPtr () {return outputMgrRCPtr_->getCommPtr();}
#endif

    // this is only used to construct ExpressionData objects so they have a way
    // to connect to the OutputMgr.  Need to refactor this so they can take
    // an N_ANP_OutputMgrAdapter reference or RCP.
    N_IO_OutputMgr * getOutputMgrPtr() { return outputMgrRCPtr_.getRawPtr(); }

    void setStepAnalysisStepNumber( int num )
    { stepAnalysisStepNumber_ = num; }

    void setStepAnalysisMaxSteps( int num )
    { stepAnalysisMaxSteps_ = num; }

    void setDCAnalysisStepNumber( int num )
    { dcAnalysisStepNumber_ = num; }

    void setDCAnalysisMaxSteps( int num )
    { dcAnalysisMaxSteps_ = num; }

    void check_output(ANP_Analysis_Mode analysis_mode)
    {
      outputMgrRCPtr_->prepareOutput(analysis_mode, *stepParamVecRCPtr_, *dcParamVecRCPtr_);
    }

    void tranOutput(double time, N_LAS_Vector & currSolutionPtr,
      N_LAS_Vector & stateVecPtr, N_LAS_Vector & storeVecPtr, bool skipPrintLineOutput=false )
    {
      outputMgrRCPtr_->output(time,
          stepAnalysisStepNumber_, stepAnalysisMaxSteps_, *stepParamVecRCPtr_,
          dcAnalysisStepNumber_, dcAnalysisMaxSteps_, *dcParamVecRCPtr_,
          & currSolutionPtr, & stateVecPtr, & storeVecPtr, skipPrintLineOutput);
    }

    void dcOutput( int dcStepNumber, N_LAS_Vector & currSolutionPtr, N_LAS_Vector & stateVecPtr, N_LAS_Vector & storeVecPtr )
    {
      outputMgrRCPtr_->output(0.0,
          stepAnalysisStepNumber_, stepAnalysisMaxSteps_, *stepParamVecRCPtr_,
          dcStepNumber, dcAnalysisMaxSteps_, *dcParamVecRCPtr_,
          & currSolutionPtr, & stateVecPtr, & storeVecPtr);
    }


    void outputRESULT( N_LAS_Vector & currSolutionPtr, N_LAS_Vector & currStatePtr, N_LAS_Vector & currStorePtr )
    {
      outputMgrRCPtr_->outputRESULT( & currSolutionPtr, & currStatePtr, & currStorePtr );
    }

    void finishOutputSTEP()
    {
      outputMgrRCPtr_->finishOutputSTEP ();
    }

    void finishOutput()
    {
      outputMgrRCPtr_->finishOutput();
    }

    bool setupInitialConditions ( N_LAS_Vector & solnVec, N_LAS_Vector & flagVec)
    {
      return outputMgrRCPtr_->setupInitialConditions( solnVec, flagVec);
    }

    void outputDCOP( N_LAS_Vector & currSolutionPtr )
    {
      outputMgrRCPtr_->outputDCOP( currSolutionPtr );
    }


    void outputMPDE ( double time, const N_LAS_Vector & solnVecPtr )
    {
      outputMgrRCPtr_->outputMPDE(time, & solnVecPtr);
    }

    void outputHB (
        const vector< double > & timePoints, const vector< double > & freqPoints,
        const N_LAS_BlockVector & timeDomainSolnVec, const N_LAS_BlockVector & freqDomainSolnVecReal,
        const N_LAS_BlockVector & freqDomainSolnVecImaginary)
    {
      outputMgrRCPtr_->outputHB(
          stepAnalysisStepNumber_, stepAnalysisMaxSteps_, *stepParamVecRCPtr_,
          timePoints, freqPoints, timeDomainSolnVec, freqDomainSolnVecReal, freqDomainSolnVecImaginary);
    }

    void outputAC (double freq, const N_LAS_Vector & solnVecRealPtr, const N_LAS_Vector & solnVecImaginaryPtr)
    {
      outputMgrRCPtr_->outputAC(freq, & solnVecRealPtr, & solnVecImaginaryPtr);
    }

    void outputMORTF ( bool origSys, const double & freq, const Teuchos::SerialDenseMatrix<int, std::complex<double> >& H )
    {
      outputMgrRCPtr_->outputMORTF( origSys, freq, H );
    }

    void resetOutputMORTF()
    {
      outputMgrRCPtr_->resetOutput();
    }

    void outputROM( const Teuchos::SerialDenseMatrix<int, double>& Ghat, const Teuchos::SerialDenseMatrix<int, double>& Chat,
                    const Teuchos::SerialDenseMatrix<int, double>& Bhat, const Teuchos::SerialDenseMatrix<int, double>& Lhat )
    {
      outputMgrRCPtr_->outputROM( Ghat, Chat, Bhat, Lhat );
    }

    void outputROM( const N_LAS_Matrix& Ghat, const N_LAS_Matrix& Chat,
                    const Teuchos::SerialDenseMatrix<int, double>& Bhat,
                    const Teuchos::SerialDenseMatrix<int, double>& Lhat )
    {
      outputMgrRCPtr_->outputROM( Ghat, Chat, Bhat, Lhat );
    }

    bool getOutputIntervals(double & initialInterval, vector < pair < double, double > > * intervalPairs)
    {
      return outputMgrRCPtr_->getOutputIntervals( initialInterval, *intervalPairs );
    }


    void outputHomotopy( const vector<string> & paramNames, const vector<double> & paramVals, N_LAS_Vector & solnVecPtr )
    {
      outputMgrRCPtr_->outputHomotopy ( paramNames, paramVals, & solnVecPtr );
    }

    Xyce::NodeNamePairMap & getAllNodes ( )
    {
      return outputMgrRCPtr_->getAllNodes ();
    }

  private:
    RefCountPtr< N_IO_OutputMgr > outputMgrRCPtr_;
    RefCountPtr< vector<N_ANP_SweepParam> > stepParamVecRCPtr_;
    RefCountPtr< vector<N_ANP_SweepParam> > dcParamVecRCPtr_;

    int stepAnalysisStepNumber_;
    int stepAnalysisMaxSteps_;
    int dcAnalysisStepNumber_;
    int dcAnalysisMaxSteps_;

};

#endif


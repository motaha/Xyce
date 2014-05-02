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
// Revision Number: $Revision: 1.38 $
//
// Revision Date  : $Date: 2014/02/24 23:49:12 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_ANP_OutputMgrAdapter_h
#define Xyce_N_ANP_OutputMgrAdapter_h

#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RefCountPtr;
using Teuchos::rcp;
#include <Teuchos_SerialDenseMatrix.hpp>

#include <N_UTL_Xyce.h>
#include <N_ANP_fwd.h>
#include <N_PDS_fwd.h>

#include <N_IO_OutputMgr.h>

class N_LAS_Vector;

namespace Xyce {
namespace Analysis {

//-------------------------------------------------------------------------
// Class         : OutputMgrAdapter
// Purpose       : Inteface class for the output manager
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//-------------------------------------------------------------------------
class OutputMgrAdapter
{

  public:
    OutputMgrAdapter( );

    virtual ~OutputMgrAdapter() {}

    void registerOutputMgr( N_IO_OutputMgr * outputMgrPtr )
    {
      outputManager_ = outputMgrPtr;
    }

    void setStepParamVec( const RefCountPtr< std::vector<SweepParam> > & paramVec )
    {
      stepParamVecRCPtr_ = paramVec;
    }

    void setDCParamVec( const RefCountPtr< std::vector<SweepParam> > & paramVec )
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

    N_PDS_Comm * getCommPtr () {
      return outputManager_->getCommPtr();
    }

    // this is only used to construct ExpressionData objects so they have a way
    // to connect to the OutputMgr.  Need to refactor this so they can take
    // an OutputMgrAdapter reference or RCP.
    N_IO_OutputMgr * getOutputMgrPtr() {
      return outputManager_;
    }

    void setStepAnalysisStepNumber( int num )
    { stepAnalysisStepNumber_ = num; }

    void setStepAnalysisMaxSteps( int num )
    { stepAnalysisMaxSteps_ = num; }

    void setDCAnalysisStepNumber( int num )
    { dcAnalysisStepNumber_ = num; }

    void setDCAnalysisMaxSteps( int num )
    { dcAnalysisMaxSteps_ = num; }

    void check_output(Xyce::Analysis::Analysis_Mode analysis_mode)
    {
      outputManager_->prepareOutput(analysis_mode, *stepParamVecRCPtr_, *dcParamVecRCPtr_);
    }

    void tranOutput(double time, N_LAS_Vector & currSolutionPtr,
      N_LAS_Vector & stateVecPtr, N_LAS_Vector & storeVecPtr, 
      std::vector<double> & objectiveVec_, 
      std::vector<double> & dOdpVec_, 
      std::vector<double> & dOdpAdjVec_,
      std::vector<double> & scaled_dOdpVec_, 
      std::vector<double> & scaled_dOdpAdjVec_,
      bool skipPrintLineOutput=false )
    {
      outputManager_->output(time,
          stepAnalysisStepNumber_, stepAnalysisMaxSteps_, *stepParamVecRCPtr_,
          dcAnalysisStepNumber_, dcAnalysisMaxSteps_, *dcParamVecRCPtr_,
          & currSolutionPtr, & stateVecPtr, & storeVecPtr, objectiveVec_,
          dOdpVec_, dOdpAdjVec_, scaled_dOdpVec_, scaled_dOdpAdjVec_,
          skipPrintLineOutput);
    }

    void dcOutput( 
        int dcStepNumber, 
        N_LAS_Vector & currSolutionPtr, N_LAS_Vector & stateVecPtr, N_LAS_Vector & storeVecPtr,
        std::vector<double> & objectiveVec_, 
        std::vector<double> & dOdpVec_, 
        std::vector<double> & dOdpAdjVec_, 
        std::vector<double> & scaled_dOdpVec_, 
        std::vector<double> & scaled_dOdpAdjVec_)
    {
      outputManager_->output(0.0,
          stepAnalysisStepNumber_, stepAnalysisMaxSteps_, *stepParamVecRCPtr_,
          dcStepNumber, dcAnalysisMaxSteps_, *dcParamVecRCPtr_,
          & currSolutionPtr, & stateVecPtr, & storeVecPtr, objectiveVec_,
            dOdpVec_, dOdpAdjVec_, scaled_dOdpVec_, scaled_dOdpAdjVec_);
    }


    void outputRESULT( N_LAS_Vector & currSolutionPtr, N_LAS_Vector & currStatePtr, N_LAS_Vector & currStorePtr )
    {
      outputManager_->outputRESULT( & currSolutionPtr, & currStatePtr, & currStorePtr );
    }

    void finishOutputSTEP()
    {
      outputManager_->finishOutputSTEP ();
    }

    void finishOutput()
    {
      outputManager_->finishOutput();
    }

    bool setupInitialConditions ( N_LAS_Vector & solnVec, N_LAS_Vector & flagVec)
    {
      return outputManager_->setupInitialConditions( solnVec, flagVec);
    }

    void outputDCOP( N_LAS_Vector & currSolutionPtr )
    {
      outputManager_->outputDCOP( currSolutionPtr );
    }


    void outputMPDE ( double time, const N_LAS_Vector & solnVecPtr )
    {
      outputManager_->outputMPDE(time, & solnVecPtr);
    }

    void outputHB (
        const std::vector< double > & timePoints, const std::vector< double > & freqPoints,
        const N_LAS_BlockVector & timeDomainSolnVec, const N_LAS_BlockVector & freqDomainSolnVecReal,
        const N_LAS_BlockVector & freqDomainSolnVecImaginary, const N_LAS_BlockVector & timeDomainStoreVec,
        const N_LAS_BlockVector & freqDomainStoreVecReal, const N_LAS_BlockVector & freqDomainStoreVecImaginary)
    {
      outputManager_->outputHB(
          stepAnalysisStepNumber_, stepAnalysisMaxSteps_, *stepParamVecRCPtr_,
          timePoints, freqPoints, timeDomainSolnVec, freqDomainSolnVecReal, freqDomainSolnVecImaginary);
    }

    void outputAC (double freq, const N_LAS_Vector & solnVecRealPtr, const N_LAS_Vector & solnVecImaginaryPtr)
    {
      outputManager_->outputAC(freq, & solnVecRealPtr, & solnVecImaginaryPtr);
    }

    void outputMORTF ( bool origSys, const double & freq, const Teuchos::SerialDenseMatrix<int, std::complex<double> >& H )
    {
      outputManager_->outputMORTF( origSys, freq, H );
    }

    void resetOutputMORTF()
    {
      outputManager_->resetOutput();
    }

    void outputROM( const Teuchos::SerialDenseMatrix<int, double>& Ghat, const Teuchos::SerialDenseMatrix<int, double>& Chat,
                    const Teuchos::SerialDenseMatrix<int, double>& Bhat, const Teuchos::SerialDenseMatrix<int, double>& Lhat )
    {
      outputManager_->outputROM( Ghat, Chat, Bhat, Lhat );
    }

    void outputROM( const N_LAS_Matrix& Ghat, const N_LAS_Matrix& Chat,
                    const Teuchos::SerialDenseMatrix<int, double>& Bhat,
                    const Teuchos::SerialDenseMatrix<int, double>& Lhat )
    {
      outputManager_->outputROM( Ghat, Chat, Bhat, Lhat );
    }

    bool getOutputIntervals(double & initialInterval, std::vector<std::pair< double, double > > * intervalPairs)
    {
      return outputManager_->getOutputIntervals( initialInterval, *intervalPairs );
    }


    void outputHomotopy( const std::vector<std::string> & paramNames, const std::vector<double> & paramVals, N_LAS_Vector & solnVecPtr )
    {
      outputManager_->outputHomotopy ( paramNames, paramVals, & solnVecPtr );
    }

    Xyce::NodeNamePairMap & getAllNodes ( )
    {
      return outputManager_->getAllNodes ();
    }

  private:
    N_IO_OutputMgr *            outputManager_;
    
    RefCountPtr< std::vector<SweepParam> > stepParamVecRCPtr_;
    RefCountPtr< std::vector<SweepParam> > dcParamVecRCPtr_;

    int stepAnalysisStepNumber_;
    int stepAnalysisMaxSteps_;
    int dcAnalysisStepNumber_;
    int dcAnalysisMaxSteps_;

};

} // namespace Analysis
} // namespace Xyce

typedef Xyce::Analysis::OutputMgrAdapter N_ANP_OutputMgrAdapter;

#endif // Xyce_N_ANP_OutputMgrAdapter_h


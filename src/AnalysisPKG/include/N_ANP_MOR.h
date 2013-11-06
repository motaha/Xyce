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
// Filename       : $RCSfile: N_ANP_MOR.h,v $
//
// Purpose        : MOR analysis class
//
// Special Notes  : Specify any "hidden" or subtle details of the class here.
//                  Portability details, error handling information, etc.
//
// Creator        : Ting Mei   
//
// Creation Date  : 01/11
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.8.2.2 $
// Revision Date  : $Date: 2013/10/03 17:23:30 $
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_ANP_MOR_h
#define Xyce_N_ANP_MOR_h

#include <Teuchos_RCP.hpp>
using Teuchos::RCP;
using Teuchos::rcp;
#include <Teuchos_SerialDenseMatrix.hpp>


// ----------   Xyce Includes   ----------
#include <N_ANP_AnalysisBase.h>
#include <N_UTL_FixedQueue.h>

// ---------- Forward Declarations ----------
class N_ANP_AnalysisManager;
class N_UTL_ExpressionData;
class N_LAS_Matrix;
class N_LAS_Vector;
class N_LAS_BlockMatrix;
class N_LAS_BlockVector;
class Amesos_BaseSolver;
class Epetra_LinearProblem;

//-------------------------------------------------------------------------
// Class         : N_ANP_MOR 
// Purpose       : MOR analysis class
// Special Notes : 
// Creator       : Ting Mei
// Creation Date : 05/24/11
//-------------------------------------------------------------------------
class N_ANP_MOR: public N_ANP_AnalysisBase 
{
  public:
    N_ANP_MOR( N_ANP_AnalysisManager * anaManagerPtr );

    ~N_ANP_MOR();
    
    bool setAnalysisParams(const N_UTL_OptionBlock & paramsBlock);

    bool run();
    bool init();
    bool reduceSystem();
    bool evalOrigTransferFunction();
    bool evalRedTransferFunction();

    bool processSuccessfulStep(bool origSys);
    bool processFailedStep();
    bool finish();
    bool handlePredictor();

    void setDCOPFlag(bool flag) { dcopFlag_ = flag; }
    bool getDCOPFlag() { return dcopFlag_; }

  private:
    
    bool isPaused; 
    bool dcopFlag_;               // true if this is a DCOP calculation.
    int morEvalSize_;
    int numPorts_;
 
    list < int > morEvalFailures_; 
    std::vector<std::string> portList_;
 
    double stepMult_;
    double fStep_;
    double currentFreq_; 
    double s0_;

    int setupSweepParam_();

    bool updateCurrentFreq_(int stepNumber);
 
    bool createOrigLinearSystem_();
    bool createRedLinearSystem_();
    
    bool updateOrigLinearSystemFreq_();
    bool updateRedLinearSystemFreq_();

    bool solveOrigLinearSystem_();
    bool solveRedLinearSystem_();

    bool sparsifyRedSystem_();

    // Original system
    RCP<N_LAS_Matrix> CPtr_;
    RCP<N_LAS_Matrix> GPtr_;
    RCP<N_LAS_Matrix> sCpG_MatrixPtr_;
    RCP<N_LAS_Matrix> redCPtr_, redGPtr_;
    RCP<N_LAS_MultiVector> RPtr_, BPtr_, VPtr_;
    std::vector<int> bMatEntriesVec_, bMatPosEntriesVec_;

    // Original system, real-equivalent form
    RCP<N_LAS_BlockMatrix> sCpG_REFMatrixPtr_;
    RCP<N_LAS_BlockVector> REFBPtr_;
    RCP<N_LAS_BlockVector> REFXPtr_; // Store solution from Amesos here.
    //RCP<N_LAS_BlockVector> LPtr_;  

    // Reduced system
    Teuchos::SerialDenseMatrix<int, double> redC_;
    Teuchos::SerialDenseMatrix<int, double> redG_;
    Teuchos::SerialDenseMatrix<int, double> redB_;
    Teuchos::SerialDenseMatrix<int, double> redL_;  // redL_ != redB_

    // Reduced system, real-equivalent form
    Teuchos::SerialDenseMatrix<int, double> sCpG_redMatrix_, sCpG_tmpMatrix_;
    Teuchos::SerialDenseMatrix<int, double> ref_redB_;

    // Transfer functions
    Teuchos::SerialDenseMatrix<int, std::complex<double> > origH_;
    Teuchos::SerialDenseMatrix<int, std::complex<double> > redH_;

    // Original system solver objects
    RCP<Amesos_BaseSolver> blockSolver_, origSolver_;
    RCP<Epetra_LinearProblem> blockProblem_, origProblem_;
};

#endif

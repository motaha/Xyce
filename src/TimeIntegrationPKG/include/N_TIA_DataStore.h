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
// Filename      : $RCSfile: N_TIA_DataStore.h,v $
//
// Purpose       : This file handles the class that defines the data arrays
//                 needed for the time integration algorithms.
//
// Special Notes :
//
// Creator       : Buddy Watts
//
// Creation Date : 6/1/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.84 $
//
// Revision Date  : $Date: 2014/02/24 23:49:26 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_TIA_DATA_STORE_H
#define Xyce_N_TIA_DATA_STORE_H

// ---------- Standard Includes ----------

#include <list>

// ----------   Xyce Includes   ----------
#include <N_TIA_TIAParams.h>
#include <N_TIA_TwoLevelError.h>

class N_LAS_MultiVector;
class N_LAS_Vector;
class N_LAS_Matrix;
class N_LAS_System;

//-----------------------------------------------------------------------------
// Class         : N_TIA_DataStore
// Purpose       : This is the class for defining data arrays needed in the
//                 time integration algorithms.
// Special Notes :
// Creator       : Buddy Watts, SNL
// Creation Date : 6/01/00
//-----------------------------------------------------------------------------
class N_TIA_DataStore
{
  public:
    N_TIA_DataStore(N_TIA_TIAParams * tiaPtr, N_LAS_System * lsPtr);
    N_TIA_DataStore(const N_TIA_DataStore& rhs);
    virtual ~N_TIA_DataStore();

    // TIA Arrays (pointers) for  Integration Solution Process:
    unsigned int solutionSize;
    unsigned int stateSize;

    // temporary vectors:
    N_LAS_Vector * tmpSolVectorPtr;
    N_LAS_Vector * tmpStaVectorPtr;
    N_LAS_Vector * tmpStaDerivPtr;
    N_LAS_Vector * tmpStaDivDiffPtr;

    N_LAS_Vector * tmpStoVectorPtr;

    // Predictors
    N_LAS_Vector * xn0Ptr;

    // Solutions:
    N_LAS_Vector * currSolutionPtr;
    N_LAS_Vector * lastSolutionPtr;
    N_LAS_Vector * oldeSolutionPtr;
    N_LAS_Vector * nextSolutionPtr;
    N_LAS_Vector * flagSolutionPtr;

    N_LAS_Vector * savedNextSolutionPtr;

    // States:
    N_LAS_Vector * currStatePtr;
    N_LAS_Vector * lastStatePtr;
    N_LAS_Vector * oldeStatePtr;
    N_LAS_Vector * nextStatePtr;

    // Storage:
    N_LAS_Vector * currStorePtr;
    N_LAS_Vector * lastStorePtr;
    N_LAS_Vector * oldeStorePtr;
    N_LAS_Vector * nextStorePtr;
    // for lead current calculations.  F component is 
    // held in the store vector, Q component is here
    N_LAS_Vector * currStoreLeadCurrQCompPtr;
    N_LAS_Vector * lastStoreLeadCurrQCompPtr;
    N_LAS_Vector * oldStoreLeadCurrQCompPtr;
    N_LAS_Vector * nextStoreLeadCurrQCompPtr;

    // Derivatives of Solutions:
    N_LAS_Vector * currSolutionDerivPtr;
    N_LAS_Vector * lastSolutionDerivPtr;
    N_LAS_Vector * oldeSolutionDerivPtr;
    N_LAS_Vector * nextSolutionDerivPtr;

    // Derivatives of States:
    N_LAS_Vector * currStateDerivPtr;
    N_LAS_Vector * lastStateDerivPtr;
    N_LAS_Vector * oldeStateDerivPtr;
    N_LAS_Vector * nextStateDerivPtr;
    
    // Derivatives of Store for lead curent calculations 
    N_LAS_Vector * currStoreLeadCurrQCompDerivPtr;
    N_LAS_Vector * lastStoreLeadCurrQCompDerivPtr;
    N_LAS_Vector * oldStoreLeadCurrQCompDerivPtr;
    N_LAS_Vector * nextStoreLeadCurrQCompDerivPtr;

    // Scaled Divided Differences
    N_LAS_Vector * currSolutionDivDiffPtr;
    N_LAS_Vector * lastSolutionDivDiffPtr;
    N_LAS_Vector * oldeSolutionDivDiffPtr;
    N_LAS_Vector * nextSolutionDivDiffPtr;

    // Scaled Divided Differences
    N_LAS_Vector * currStateDivDiffPtr;
    N_LAS_Vector * lastStateDivDiffPtr;
    N_LAS_Vector * oldeStateDivDiffPtr;
    N_LAS_Vector * nextStateDivDiffPtr;

    // Error Vectors
    N_LAS_Vector * errWtVecPtr;
    N_LAS_Vector * absErrTolPtr;
    N_LAS_Vector * relErrTolPtr;

    // Jacobian and RHS
    N_LAS_Matrix * JMatrixPtr;
    N_LAS_Vector * RHSVectorPtr;
#ifdef Xyce_DEBUG_DEVICE
    N_LAS_Vector * JdxpVectorPtr; // only used in the device manager for debug purposes these days.
#endif

    // NonLinear Solution Vectors
    N_LAS_Vector * newtonCorrectionPtr;

    // Mask for error norms (to allow some equations not to take part in
    // weighted norms)
    N_LAS_Vector * deviceMaskPtr;
    // TVR: I toyed with the idea of having this flag here, but went with 
    // keeping it in the LAS_System instead --- we call an accessor method
    // to set and get the flag when we create or use the mask.  
    // flag showing whether mask is trivial or not
    //   bool nonTrivialDeviceMask;
    
    // this is a simple vector indicating var types.  Now
    // we just handle V and I differently, but we could do more with this
    std::vector<char> varTypeVec;
    
    // To remove conditionals from setErrorWtVector() we'll create
    // lists of indexes of unknows that are handled in different ways
    std::vector<int> indexVVars;
    std::vector<int> indexIVars;
    std::vector<int> indexMaskedVars;
    int numVVars, numIVars, numMaskedVars;
    bool indexVecsInitialized;
    
    // limiter flag:
    bool limiterFlag;

    // 2-level information:
    std::vector<N_TIA_TwoLevelError> innerErrorInfoVec;

    // new-DAE data (originally from the new-DAE derrived class)
    // Error Vectors
    N_LAS_Vector * qErrWtVecPtr;

    // DAE formulation vectors
    N_LAS_Vector * daeQVectorPtr;
    N_LAS_Vector * daeFVectorPtr;

    // voltage limiting vectors
    N_LAS_Vector * dFdxdVpVectorPtr;
    N_LAS_Vector * dQdxdVpVectorPtr;

    // DAE formulation matrices
    N_LAS_Matrix * dQdxMatrixPtr;
    N_LAS_Matrix * dFdxMatrixPtr;

    // HB temporary Matvec storage vectors
    N_LAS_Vector * dQdxVecVectorPtr;
    N_LAS_Vector * dFdxVecVectorPtr;

    // History arrays
    std::vector<N_LAS_Vector*> xHistory;
    std::vector<N_LAS_Vector*> qHistory;
    std::vector<N_LAS_Vector*> sHistory;    // state history
    std::vector<N_LAS_Vector*> stoHistory;  // store history
    std::vector<N_LAS_Vector*> stoLeadCurrQCompHistory;  // store history for lead current Q component.

    // Predictors
    N_LAS_Vector * qn0Ptr;
    N_LAS_Vector * qpn0Ptr;

    N_LAS_Vector * sn0Ptr;
    N_LAS_Vector * spn0Ptr;
  
    N_LAS_Vector * ston0Ptr;
    N_LAS_Vector * stopn0Ptr;
    
    N_LAS_Vector * stoQCn0Ptr;
    N_LAS_Vector * stoQCpn0Ptr;

    // Nonlinear solution vector:
    N_LAS_Vector * qNewtonCorrectionPtr;
    N_LAS_Vector * sNewtonCorrectionPtr;
    N_LAS_Vector * stoNewtonCorrectionPtr;
    N_LAS_Vector * stoLeadCurrQCompNewtonCorrectionPtr;

    // Step-size selection temporary vectors
    N_LAS_Vector * delta_x;
    N_LAS_Vector * delta_q;

    // Temporary vectors for WaMPDE interpolation
    N_LAS_Vector * tmpXn0APtr;
    N_LAS_Vector * tmpXn0BPtr;

    // These are for MPDE fast time scale points
    std::vector<double> timeSteps;
    std::vector<bool> timeStepsBreakpointFlag;
    std::vector<N_LAS_Vector*> fastTimeSolutionVec;
    std::vector<N_LAS_Vector*> fastTimeStateVec;
    std::vector<N_LAS_Vector*> fastTimeQVec;
    std::vector<N_LAS_Vector*> fastTimeStoreVec;
    //std::vector<N_LAS_Vector*> pastQVecHistory;

    // these are placeholders until the transient sensitivites are set 
    // up properly
    std::vector<double> objectiveVec_; 
    std::vector<double> dOdpVec_; 
    std::vector<double> dOdpAdjVec_;
    std::vector<double> scaled_dOdpVec_; 
    std::vector<double> scaled_dOdpAdjVec_;

  protected:

  private:


    // DataStore Functions
  public:
    void initializeDataArrays();
    void enableOrderOneStart();

    void updateSolDataArrays();
    bool updateStateDataArrays();

    void outputSolDataArrays(std::ostream &os);
    virtual void setConstantHistory();
    virtual void setZeroHistory();
    virtual void setErrorWtVector();
    virtual double WRMS_errorNorm();

    virtual double partialErrorNormSum();
    virtual double partialQErrorNormSum();

    virtual double partialSum_m1(int currentOrder);
    virtual double partialSum_m2(int currentOrder);

    virtual double globalLength ();
    void computeDividedDifferences();

    void computeDivDiffsBlock(const std::list<index_pair> & solGIDList,
                              const std::list<index_pair> & staGIDList);

    void printOutPointers ();
    bool equateTmpVectors ();
    bool usePreviousSolAsPredictor ();

    void outputPredictedSolution(std::ostream &os);
    void outputPredictedDerivative(std::ostream &os);

    virtual void stepLinearCombo ();

    virtual double partialSum_p1(int currentOrder, int maxOrder);
    virtual double partialSum_q1();

    virtual double delta_x_errorNorm_m1();
    virtual double delta_x_errorNorm_m2();
    virtual double delta_x_errorNorm_p1();
    virtual double delta_x_errorNorm_q1();

    virtual bool getSolnVarData( const int & gid, std::vector<double> & varData );
    virtual bool getStateVarData( const int & gid, std::vector<double> & varData );
    virtual bool setSolnVarData( const int & gid, const std::vector<double> & varData );
    virtual bool setStateVarData( const int & gid, const std::vector<double> & varData );
    virtual bool getStoreVarData( const int & gid, std::vector<double> & varData );
    virtual bool setStoreVarData( const int & gid, const std::vector<double> & varData );

    N_LAS_System * lasSysPtr;

    bool setNextSolVectorPtr (N_LAS_Vector * solVecPtr);
    bool unsetNextSolVectorPtr ();
  
    bool resetAll ();
    bool resetFastTimeData ();

  protected:

    N_TIA_TIAParams * tiaParamsPtr_;
    double          * dataBlockPtr;
    N_TIA_DataStore & operator=(const N_TIA_DataStore& rhs);

    bool nextSolPtrSwitched;
    

  private:
    N_TIA_DataStore();

};

#endif // Xyce_N_TIA_DATA_STORE_H


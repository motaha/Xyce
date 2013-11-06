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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_NLS_NOX_ParameterSet.h,v $
//
// Purpose        : Specification file which declares an interface common to
//                  all supported nonlinear solver algorithms.  The Manager
//                  class uses this interface to call a concrete algorithm.
//
// Special Notes  : This is the "Strategy" class in the Strategy design
//                  pattern.
//
// Creator        : Scott A. Hutchinson, SNL, Parallel Computational Sciences
//
// Creation Date  : 04/28/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.40.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:47 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef Xyce_N_NLS_NOX_ParameterSet_h
#define Xyce_N_NLS_NOX_ParameterSet_h

// ---------- Standard Includes ----------

#include <vector>

// ----------   Xyce Includes   ----------
#include <N_UTL_Misc.h>
#include <N_UTL_NoCase.h>

// ----------   NOX Includes   ----------

#include "NOX.H"
#include "Teuchos_RefCountPtr.hpp"

// ---------- Forward Declarations ----------

class N_LAS_Vector;
class N_LAS_System;
class N_NLS_ReturnCodes;
class N_LOA_Loader;
class N_PDS_Comm;
namespace N_NLS_NOX {
  class AugmentLinSys;
}
namespace NOX {
  namespace Status {
    class Test;
  }
}


//-----------------------------------------------------------------------------
// Class         : N_NLS_NonLinearSolver
// Purpose       : Nonlinear Solver Abstract Class
// Creator       : Tammy Kolda, SNL, 8950
// Creation Date : 2/5/02
//-----------------------------------------------------------------------------

namespace N_NLS_NOX {
class ParameterSet {

public:

  ParameterSet(AnalysisMode mode);
  ~ParameterSet();
  bool setOptions(const N_UTL_OptionBlock& OB);
  bool setLocaOptions(const N_UTL_OptionBlock& OB, bool saveCopy=true);
  bool applySavedLocaOptions()
  {
    bool ret=true;
    if (savedLocaOptions_)
    {
      ret=setLocaOptions(savedLocaOB_, false);
    }
    return ret;
  }

  bool setOutputOptions(int myPID, int outputProcess);
#ifdef Xyce_NLS_MASKED_WRMS_NORMS
  bool createStatusTests(N_LAS_Vector** currSolnVectorPtrPtr,
			 N_LOA_Loader& l, vector<char> & varTypeVec,
			 bool nonTrivialDeviceMaskFlag, N_LAS_Vector * maskVectorPtr);
#else
  bool createStatusTests(N_LAS_Vector** currSolnVectorPtrPtr,
			 N_LOA_Loader& l, vector<char> & varTypeVec);
#endif
  Teuchos::RefCountPtr<NOX::StatusTest::Generic> getStatusTests();
  bool getVectorParam (const string &, int, double &);
  bool getVectorParam (const string &, int, string &);
  int getVectorParamSize(const string& vectorName);
  int getStatusTestReturnCode() const;
  void setStatusTestReturnCodes (const N_NLS_ReturnCodes & retCodesTmp);

  Teuchos::RefCountPtr<Teuchos::ParameterList> getAllParams();
  Teuchos::RefCountPtr<Teuchos::ParameterList> getNoxParams();
  Teuchos::RefCountPtr<Teuchos::ParameterList> getLocaParams();
  Teuchos::RefCountPtr<Teuchos::ParameterList> getDebugParams();
  int getNoxSolverType() const;

#ifdef Xyce_DEBUG_NONLINEAR
  inline int  getDebugLevel() const
  {
    return debugLevel_;
  };

  inline int getDebugMinTimeStep() const
  {
    return debugMinTimeStep_;
  };

  inline int getDebugMaxTimeStep() const
  {
    return debugMaxTimeStep_;
  };

  inline double getDebugMinTime() const
  {
    return debugMinTime_;
  };

  inline double getDebugMaxTime() const
  {
    return debugMaxTime_;
  };

  inline bool getScreenOutputFlag () const
  {
    return screenOutputFlag_;
  };
#endif

  double getMaxNormF() const;
  int getMaxNormFindex () const;

  inline bool isParamsSet() const
  {
    return isParamsSet_;
  };

  inline void set_gstepping_min_value (double val)
  {
    gstepping_min_value_=val;
  }

  inline void set_gstepping_minimum_conductance (double val)
  {
    gstepping_minimum_conductance_ = val;
  }

  Teuchos::RefCountPtr<N_NLS_NOX::AugmentLinSys>
    createAugmentLinearSystem(N_LAS_System* ls) const;

  // Create augmented linear system, DCOP restart version
#ifdef Xyce_PARALLEL_MPI
  Teuchos::RefCountPtr<N_NLS_NOX::AugmentLinSys>
    createAugmentLinearSystem(N_LAS_System* ls,
       Xyce::NodeNamePairMap & op,
       Xyce::NodeNamePairMap & allNodes, N_PDS_Comm * pdsCommPtr) const;
#else
  Teuchos::RefCountPtr<N_NLS_NOX::AugmentLinSys>
    createAugmentLinearSystem(N_LAS_System* ls,
       Xyce::NodeNamePairMap & op,
       Xyce::NodeNamePairMap & allNodes) const;
#endif

  // Create augmented linear system, IC version
    Teuchos::RefCountPtr<N_NLS_NOX::AugmentLinSys>
    createAugmentLinearSystem(N_LAS_System* ls,
                              Xyce::NodeNamePairMap & op,
                              bool gminStepping=false) const;

private:

  void unsupportedOption_(const string& tag);
  bool parseOptionBlock_(const N_UTL_OptionBlock& OB);

private:

  Teuchos::RefCountPtr<Teuchos::ParameterList> allParams_;
  Teuchos::ParameterList& noxParams_;
  Teuchos::ParameterList& locaParams_;
  Teuchos::ParameterList& debugParams_;
  Teuchos::ParameterList statusTestParams_;

  map<string, vector<N_UTL_Param> > vectorParams;

  // Combo of all tests
  Teuchos::RefCountPtr<NOX::StatusTest::Combo> comboPtr_;

  // Vector containing all the tests we create
  vector< Teuchos::RefCountPtr<NOX::StatusTest::Generic> > tests_;

  bool isParamsSet_;
  bool isStatusTestsSet_;

  AnalysisMode mode_;

  // For special Xyce specific directions (solvers 7 and 8).
  //NOX::Parameter::DirectionConstructorT<N_NLS_NOX::FastNewton> dirFN;
  //NOX::Parameter::DirectionConstructorT<N_NLS_NOX::NewtonSDCombo> dirNSDCombo;

  int noxSolver;

  // This object acts as a factory for the augment system strategy.
  enum VoltageListType {
    VLT_DOFS,
    VLT_Node,
    VLT_None
  };

  VoltageListType voltageListType_;

  double voltageScaleFactor_;

  // Minimum parameter value used in gstepping
  double gstepping_min_value_;

  // Evil, evil gmin stepping hack, a residual conductance that should almost
  // always be zero.
  double gstepping_minimum_conductance_;

  bool savedLocaOptions_;
  N_UTL_OptionBlock savedLocaOB_;

#ifdef Xyce_DEBUG_NONLINEAR
  // debug info:
  int  debugLevel_;
  int debugMinTimeStep_;
  int debugMaxTimeStep_;
  double debugMinTime_;
  double debugMaxTime_;
  bool screenOutputFlag_;
#endif
};
} // namespace N_NLS_NOX

#endif


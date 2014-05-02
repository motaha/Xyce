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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_NLS_NOX_Interface.h,v $
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
// Revision Number: $Revision: 1.53 $
//
// Revision Date  : $Date: 2014/02/24 23:49:24 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef Xyce_N_NLS_NOX_Interface_h
#define Xyce_N_NLS_NOX_Interface_h

#include <N_UTL_Misc.h>
#include <N_PDS_fwd.h>
#include <N_IO_fwd.h>

#include <N_NLS_Manager.h>	// defines AnalysisMode
#include <N_NLS_NonLinearSolver.h>
#include <N_NLS_NOX_ParameterSet.h>

class N_PDS_Manager;

namespace NOX {
  namespace Parameter {
    class List;
  }
  namespace StatusTest {
    class Generic;
  }
}

namespace LOCA {
  class Stepper;
  namespace MultiContinuation {
    class AbstractGroup;
  }
  class GlobalData;
  namespace StatusTest {
    class Wrapper;
  }
}

namespace N_NLS_NOX {
  class SharedSystem;
  class Group;
  class Vector;
}

namespace N_NLS_LOCA {
  class Group;
}

class N_LOA_Loader;

//-----------------------------------------------------------------------------
// Class         : N_NLS_NonLinearSolver
// Purpose       : Nonlinear Solver Abstract Class
// Creator       : Tammy Kolda, SNL, 8950
// Creation Date : 2/5/02
//-----------------------------------------------------------------------------

namespace N_NLS_NOX {
class Interface : public N_NLS_NonLinearSolver
{

public:

  Interface(N_IO_CmdParse & cp);
  ~Interface();

  bool setOptions(const N_UTL_OptionBlock& OB);
  bool setTranOptions(const N_UTL_OptionBlock& OB);
  bool setHBOptions(const N_UTL_OptionBlock& OB);

  
  bool setLocaOptions(const N_UTL_OptionBlock& OB);
  bool setDCOPRestartOptions(const N_UTL_OptionBlock& OB);
  bool setICOptions(const N_UTL_OptionBlock& OB);
  bool setNodeSetOptions(const N_UTL_OptionBlock& OB);
  bool initializeAll();
  int solve (N_NLS_NonLinearSolver * nlsTmpPtr = NULL);

  int takeFirstSolveStep (N_NLS_NonLinearSolver * nlsTmpPtr = NULL);
  int takeOneSolveStep ();

  Teuchos::RefCountPtr<N_NLS_LOCA::Group> getSolutionGroup ();

  int getNumIterations() const;
  double getMaxNormF() const;
  int getMaxNormFindex() const;

  int getDebugLevel() const;
  bool getScreenOutputFlag () const;
  double getDebugMinTime() const;
  double getDebugMaxTime() const;
  int getDebugMinTimeStep() const;
  int getDebugMaxTimeStep() const;
  bool getMMFormat () const;

  // Returns the continuation step number if available.
  int getContinuationStep() const;

  // Returns the parameter number:
  int getParameterNumber() const;

  // Returns true if this is the first continuation param
  bool isFirstContinuationParam() const;

  // Returns true if this is the first solve has been completed
  bool isFirstSolveComplete() const;
  bool getLocaFlag ();
  void setAnalysisMode(AnalysisMode mode);
  void resetAll (AnalysisMode mode);
  bool copySolnVectors();

  // Returns flag for Matrix free loads
  bool getMatrixFreeFlag();

  bool computeF();

  bool computeJacobian();

  bool applyJacobian(const N_LAS_Vector& input, N_LAS_Vector& result);

  bool computeNewton(Teuchos::ParameterList& p);

  bool computeGradient();

  N_LOA_Loader& getLoader() const;

protected:
  // Resets the stepper by destroying and reallocating.
  void resetStepper(const Teuchos::RefCountPtr<LOCA::GlobalData>& gd,
		    const Teuchos::RefCountPtr<LOCA::MultiContinuation::AbstractGroup>& initialGuess,
		    const Teuchos::RefCountPtr<NOX::StatusTest::Generic>& test,
		    const Teuchos::RefCountPtr<Teuchos::ParameterList>& p);

  // Functions for DC_OP restart and other initial condition options.
  bool opStartCont0 (ParameterSet* paramsPtr);
  bool opStartCont1 (ParameterSet* paramsPtr);

  bool icCont (ParameterSet* paramsPtr);
  bool icCont3 (ParameterSet* paramsPtr);

  bool nodesetCont0 (ParameterSet* paramsPtr);
  bool nodesetCont1 (ParameterSet* paramsPtr);
private:

  // Parameters for DC_OP
  N_NLS_NOX::ParameterSet dcParams_;

  bool DCOPused_;
  bool DCOPspecified_;
  bool ICspecified_;
  bool NODESETspecified_;

  // Parameters for Transient
  N_NLS_NOX::ParameterSet transientParams_;

  N_NLS_NOX::ParameterSet hbParams_;

  // Shared system
  N_NLS_NOX::SharedSystem* sharedSystemPtr_;

  // Global data for loca groups
  Teuchos::RefCountPtr<LOCA::GlobalData> globalDataPtr_;

  // LOCA Wrapper Status Tests
  Teuchos::RefCountPtr<LOCA::StatusTest::Wrapper> locaTransientStatusTestPtr_;
  Teuchos::RefCountPtr<LOCA::StatusTest::Wrapper> locaDCOpStatusTestPtr_;
  Teuchos::RefCountPtr<LOCA::StatusTest::Wrapper> locaStatusTestPtr_;

  Teuchos::RefCountPtr<LOCA::StatusTest::Wrapper> locaHBStatusTestPtr_;
  
  // Nox group
  Teuchos::RefCountPtr<N_NLS_LOCA::Group> groupPtr_;

  // NOX Solver
  Teuchos::RefCountPtr<NOX::Solver::Generic> solverPtr_;

  // LOCA Stepper
  Teuchos::RefCountPtr<LOCA::Stepper> stepperPtr_;

  // Current analysis mode
  AnalysisMode mode_;

  // Whether or not we should use the current analysis mode
  bool usemode_;

  // save the parameters mode.
  AnalysisMode lastParametersMode_;
  AnalysisMode parametersMode_;

  bool copiedGroupFlag_;

  // Keep track of whether to set the linear solver tolerance
  // (Only do this if adaptive forcing is on)
  bool setAZ_Tol_DC;
  bool setAZ_Tol_Transient;

  //are we on the first LOCA continuation parameter?
  bool isFirstContinuationParam_;

  //is first solve completed? 
  bool firstSolveComplete_;

  //parameter index
  int iParam_;
};
} // namespace N_NLS_NOX

#endif


//-------------------------------------------------------------------------
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
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_NLS_NOX_PseudoTransientSolver.h,v $
//
// Purpose        : 
//
// Special Notes  : 
//
// Creator        : 
//
// Creation Date  : 
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.15 $
//
// Revision Date  : $Date: 2014/02/24 23:49:24 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------


#ifndef NOX_SOLVER_PSEUDOTRANSIENTSOLVER_H
#define NOX_SOLVER_PSEUDOTRANSIENTSOLVER_H

#include <N_UTL_Misc.h>

#include "NOX_Solver_Generic.H"	         // base class

#include "NOX_LineSearch_Generic.H"      // class data element
#include "NOX_Direction_Generic.H"       // class data element

#include "NOX_Solver_PrePostOperator.H"  // class data element
#include "NOX_Utils.H"		         // class data element
#include "NOX_StatusTest_FiniteValue.H"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"     

// Forward declarations
namespace N_NLS_LOCA {
  class Group;
}
namespace NOX {
  namespace StatusTest {
    class Generic;
  }
  class GlobalData;
}
namespace N_NLS_NOX {
  class AugmentLinSys;
}


namespace N_NLS_NOX {

  class PseudoTransientBased : public NOX::Solver::Generic {

public:

  //! Constructor
  PseudoTransientBased(const Teuchos::RefCountPtr<N_NLS_NOX::AugmentLinSys>& als,
		       const Teuchos::RefCountPtr<NOX::Abstract::Group>& grp, 
		       const Teuchos::RefCountPtr<NOX::StatusTest::Generic>& tests, 
		       const Teuchos::RefCountPtr<Teuchos::ParameterList>& params,
		       double initialStepSize,
		       double minStepSize,
		       double maxStepSize);

  //! Destructor
  virtual ~PseudoTransientBased();

  virtual void reset(const NOX::Abstract::Vector& initial_guess);
  virtual void reset(const NOX::Abstract::Vector& initial_guess,
                     const Teuchos::RCP<NOX::StatusTest::Generic>& test);
  virtual NOX::StatusTest::StatusType getStatus();
  virtual NOX::StatusTest::StatusType step();
  virtual NOX::StatusTest::StatusType solve();
  virtual const NOX::Abstract::Group& getSolutionGroup() const;
  virtual const NOX::Abstract::Group& getPreviousSolutionGroup() const;
  virtual int getNumIterations() const;
  virtual const Teuchos::ParameterList& getList() const;

  //! Return the line search step size from the current iteration
  virtual double getStepSize() const;

  //! Return the pseudo transient step size.
  virtual double getPseudoTransientStepSize() const;

protected:
  
  //! Print out initialization information and calcuation the RHS.
  virtual void init();

  //! Prints the current iteration information.
  virtual void printUpdate();

protected:
  
  //! Global Data.
  Teuchos::RefCountPtr<NOX::GlobalData> globalData;

  //! RCP to the strategy for augmenting the linear system.
  Teuchos::RefCountPtr<N_NLS_NOX::AugmentLinSys> augmentLSStrategy;

  //! Current solution.
  Teuchos::RefCountPtr<NOX::Abstract::Group> solnPtr;		

  //! Previous solution pointer. 
  /*! We have both a pointer and a reference because we need to create
    a DERIVED object and then want to have a reference to it. */
  Teuchos::RefCountPtr<NOX::Abstract::Group> oldSolnPtr;	
  //! Previous solution reference.
  NOX::Abstract::Group& oldSoln;	

  //! Current search direction.pointer.
  /*! We have both a pointer and a reference because we need to create
    a DERIVED object and then want to have a reference to it. */
  Teuchos::RefCountPtr<NOX::Abstract::Vector> dirPtr;
  //! Current search direction.reference.
  NOX::Abstract::Vector& dir;

  //! Stopping test.
  Teuchos::RefCountPtr<NOX::StatusTest::Generic> testPtr;		

  //! Input parameters.
  Teuchos::RefCountPtr<Teuchos::ParameterList> paramsPtr;	

  //! Utils
  NOX::Utils& utils;

  //! Linesearch. 
  Teuchos::RCP<NOX::LineSearch::Generic> lineSearch; 

  //! %Search %Direction. 
  Teuchos::RCP<NOX::Direction::Generic> direction; 

  //! Current step.
  double step_;			

  //! Number of nonlinear iterations.
  int nIter;                    

  //! %Status of nonlinear solver.
  NOX::StatusTest::StatusType status;

  //! Pointer to a user defined NOX::Abstract::PrePostOperator object.
  NOX::Solver::PrePostOperator prePostOperator;

  double initialStepSize_;
  double minStepSize_;
  double maxStepSize_;
  double stepSize_;
  double previousStepSize_;
  double scaleFactor_;

  N_NLS_LOCA::Group* group_;
  N_NLS_LOCA::Group* previousGroup_;
  
  NOX::StatusTest::FiniteValue fvTest_;

  //! Type of check to use for status tests.
  NOX::StatusTest::CheckType checkType;

};
} // namespace NOX

#endif


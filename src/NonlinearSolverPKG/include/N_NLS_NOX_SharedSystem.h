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
// Filename       : $RCSfile: N_NLS_NOX_SharedSystem.h,v $
//
// Purpose        : Interface to let multiple N_NLS::NOX::Group's
//                  share a single system of RHS Vector, Jacobian
//                  matrix, Newton vector, and gradient vector.
//                  Closely related to the NOX::Epetra::SharedSystem
//                  class.
//
// Special Notes  :
//
// Creator        : Tammy Kolda, NLS, 8950
//
// Creation Date  : 01/31/02
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.20.6.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:47 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef Xyce_N_NLS_NOX_SharedSystem_h
#define Xyce_N_NLS_NOX_SharedSystem_h

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

// ----------   NOX Includes   ----------

#include "N_NLS_NOX_Vector.h"
#include "N_NLS_NOX_Group.h"

// ---------- Forward Declarations ----------

class N_LAS_Vector;
class N_LAS_Matrix;
class Ifpack_IlukGraph;
class Ifpack_CrsRiluk;

namespace N_NLS_NOX {
  class Interface;
}

namespace NOX {
  namespace Parameter {
    class List;
  }
}

//-----------------------------------------------------------------------------
// Class         : N_NLS::NOX::SharedSystem
//
// Purpose       :
//
//      Interface to let multiple N_NLS::NOX::Group's share the
//      vectors and matrices in the Xyce nonlinear solver.
//
//      Closely related conceptually to the
//      NOX::Epetra::SharedJacobian class.
//
// Creator       : Tammy Kolda, SNL, 8950
//
// Creation Date : 1/31/02
//-----------------------------------------------------------------------------

namespace N_NLS_NOX {
class SharedSystem {

public:

  //---------------------------------------------------------------------------
  // Function      : SharedSystem
  // Purpose       : Constructor. Creates a shared system containing
  //                 the soln vector, the previous solution vector,
  //                 the RHS vector, the Newton vector, the Jacobian
  //                 matrix, and a reference to the interface used to
  //                 call the evaluation functions.
  //---------------------------------------------------------------------------
  SharedSystem(N_LAS_Vector& x,
	       N_LAS_Vector& f,
	       N_LAS_Matrix& jacobian,
	       N_LAS_Vector& newton,
	       N_LAS_Vector& gradient,
	       N_LAS_System& lasSys,
	       N_NLS_NOX::Interface& interface);

  //---------------------------------------------------------------------------
  // Function      : ~SharedSystem
  // Purpose       : Destructor.
  //---------------------------------------------------------------------------
  ~SharedSystem();


  //---------------------------------------------------------------------------
  // Function      : reset
  // Purpose       : reset the Xyce fill objects - pointers may have changed!
  //---------------------------------------------------------------------------
  void reset(N_LAS_Vector& x,
	     N_LAS_Vector& f,
	     N_LAS_Matrix& jacobian,
	     N_LAS_Vector& newton,
	     N_LAS_Vector& gradient,
	     N_LAS_System& lasSys,
	     N_NLS_NOX::Interface& interface);

  //---------------------------------------------------------------------------
  // Function      : isJacobianOwner
  // Purpose       : Verify that the group pointed to by grp is owner of the
  //                 Jacobian matrix.
  //---------------------------------------------------------------------------
  inline bool isJacobianOwner(const Group* grp) const
  {
    return (grp == ownerOfJacobian_);
  };

  //---------------------------------------------------------------------------
  // Function      : areStateVectors
  // Purpose       : To compute a Jacobian, the state vectors must be 
  //                 updated with respect to the solution in the group.  
  //                 However, the state vectors are updated ONLY during 
  //                 calls to compute the residual.  This method checks 
  //                 to see if the state vectors still correspond to this 
  //                 group.  Returns true if state vectors are correct.
  //---------------------------------------------------------------------------
  inline bool areStateVectors(const Group* grp) const
  {
    return (grp == ownerOfStateVectors_);
  };

  //---------------------------------------------------------------------------
  // Function      : computeF
  // Purpose       : Compute the F corresponding to the current
  //                 primary solution vector. Makes the primary
  //                 solution vector owner in to the owner of the F.
  //---------------------------------------------------------------------------
  bool computeF(const Vector& solution, Vector& F, const Group* grp);

  //---------------------------------------------------------------------------
  // Function      : computeJacobian
  // Purpose       : Compute the Jacobian corresponding to the current
  //                 primary solution vector. 
  //---------------------------------------------------------------------------
  bool computeJacobian(Group* grp);

  //---------------------------------------------------------------------------
  // Function      : computeNewton
  // Purpose       : Compute the Newton corresponding to the current
  //                 primary solution vector. 
  //---------------------------------------------------------------------------
  bool computeNewton(const Vector& F, Vector& Newton,
		     Teuchos::ParameterList& params);

  //---------------------------------------------------------------------------
  // Function      : computeGradient
  // Purpose       : Compute the Gradient corresponding to the current
  //                 primary solution vector. 
  //---------------------------------------------------------------------------
  bool computeGradient(const Vector& F, Vector& Gradient);

  bool applyJacobian(const Vector& input, Vector& result) const;

  bool applyJacobianTranspose(const Vector& input, Vector& result) const;

  N_NLS_NOX::Vector& getSolutionVector();

  // Take ownership of const Jacobian.
  const N_LAS_Matrix& getJacobian() const;

  // Take ownership of Jacobian and get a reference to it.
  N_LAS_Matrix& getJacobian(const Group* grp);

  // Take ownership of the state vectors.
  void getStateVectors(const Group* grp);

  // Get a pointer to the N_LAS_System object
  N_LAS_System* getLasSystem();

#ifdef Xyce_DEBUG_NONLINEAR
  // Use for debugging (corresponding to the ones in N_NLS_NonLinearSolver).
  void debugOutput1 (N_LAS_Matrix & jacobian, N_LAS_Vector & rhs);
  void debugOutput3 (N_LAS_Vector & dxVector, N_LAS_Vector & xVector);
#endif

  // Preconditioning objects for the Group::applyRightPreconditioning method
  bool computePreconditioner();
  bool deletePreconditioner();
  bool applyRightPreconditioning(bool useTranspose, 
				 Teuchos::ParameterList& params,
				 const Vector& input, 
				 Vector& result);

  // This is used to construct vectors in the group.  We used to
  // clone a time integrator vector but when the DC Op point fails,
  // somewhere (I have no idea where) it decides to delete vectors
  // before the Nonlinear solver can delete theirs.  This causes
  // a seg fault.
  // 
  N_NLS_NOX::Vector* cloneSolutionVector() const;

  // Take ownership of const newton vector
  const N_NLS_NOX::Vector & getNewtonVector() const;

  void printSoln() {xyceSolnPtr_->print();}
  void printRes() {xyceFPtr_->print();}

private:

  // Views of xyce objects used in the fill
  Vector* xyceSolnPtr_;                     // Solution vector
  Vector* xyceFPtr_;                        // Residual Vector
  N_LAS_Matrix* xyceJacobianPtr_;           // Jacobian matrix
  Vector* xyceNewtonPtr_;                   // Newton Vector
  Vector* xyceGradientPtr_;                 // gradient Vector
  N_LAS_System* xyceLasSysPtr_;             // LAS System
  N_NLS_NOX::Interface* xyceInterfacePtr_;  // Nonlinear Solver Interface

  // Flag for Matrix Free Loads tscoffe/tmei 07/29/08
  bool matrixFreeFlag_;

  const Group* ownerOfJacobian_;
  const Group* ownerOfStateVectors_;

  // Ifpack preconditioning objects for applyRightPreconditioning method
  mutable Ifpack_IlukGraph* ifpackGraphPtr_;
  mutable Ifpack_CrsRiluk* ifpackPreconditionerPtr_;

}; // class SharedSystem
} // namespace N_NLS_NOX

#endif // Xyce_N_NLS_NOX_SharedSystem_h


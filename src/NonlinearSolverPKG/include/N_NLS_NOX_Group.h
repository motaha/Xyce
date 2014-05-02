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
// Filename       : $RCSfile: N_NLS_NOX_Group.h,v $
//
// Purpose        : Interface to Xyce for NOX groups.
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
// Revision Number: $Revision: 1.19 $
//
// Revision Date  : $Date: 2014/02/24 23:49:24 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef Xyce_N_NLS_NOX_Group_h
#define Xyce_N_NLS_NOX_Group_h

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

// ----------   NOX Includes   ----------

#include "NOX_Abstract_Group.H"
#include "Teuchos_RefCountPtr.hpp"

// ---------- Forward Declarations ----------
namespace N_NLS_NOX {
  class Vector;
  class SharedSystem;
}
namespace NOX {
  namespace Abstract {
    class Vector;
    class Group;
  }
  namespace Parameter {
    class List;
  }
}
class Ifpack_IlukGraph;
class Ifpack_CrsRiluk;


// N_NLS_NOX namespace is for all NOX-related classes in the Xyce
// Nonlinear Solver Package
namespace N_NLS_NOX {


//-----------------------------------------------------------------------------
// Class         : N_NLS::NOX::Group
//
// Purpose       :
//
//      NOX Group Interface for Xyce
//
// Creator       : Tammy Kolda, SNL, 8950
//
// Creation Date : 2/1/02
//-----------------------------------------------------------------------------

class Group : public virtual NOX::Abstract::Group {

public:

  Group(SharedSystem& s);
  Group(const Group& source, NOX::CopyType type = NOX::DeepCopy);

  ~Group();

  NOX::Abstract::Group& operator=(const Group& source);
  NOX::Abstract::Group& operator=(const NOX::Abstract::Group& source);

  void setX(const Vector& input);
  void setX(const NOX::Abstract::Vector& input);

  void computeX(const Group& grp, const Vector& d, double step);
  void computeX(const NOX::Abstract::Group& grp, const NOX::Abstract::Vector& d, double step);

  NOX::Abstract::Group::ReturnType computeF();
  NOX::Abstract::Group::ReturnType computeJacobian();
  NOX::Abstract::Group::ReturnType computeGradient();
  NOX::Abstract::Group::ReturnType computeNewton(Teuchos::ParameterList& params);
  NOX::Abstract::Group::ReturnType applyJacobian(const Vector& input, Vector& result) const;
  NOX::Abstract::Group::ReturnType applyJacobian(const NOX::Abstract::Vector& input, NOX::Abstract::Vector& result) const;

  NOX::Abstract::Group::ReturnType applyJacobianTranspose(const Vector& input, Vector& result) const;
  NOX::Abstract::Group::ReturnType applyJacobianTranspose(const NOX::Abstract::Vector& input, NOX::Abstract::Vector& result) const;
  NOX::Abstract::Group::ReturnType 
    applyJacobianInverse(Teuchos::ParameterList& params, 
				const Vector& input, Vector& result) const;
  NOX::Abstract::Group::ReturnType 
    applyJacobianInverse(Teuchos::ParameterList& params, 
				const NOX::Abstract::Vector& input,
				NOX::Abstract::Vector& result) const;
  
  NOX::Abstract::Group::ReturnType 
    applyRightPreconditioning(bool useTranspose,
				     Teuchos::ParameterList& params,
				     const Vector& input, 
				     Vector& result) const;
  NOX::Abstract::Group::ReturnType 
    applyRightPreconditioning(bool useTranspose,
				     Teuchos::ParameterList& params,
				     const NOX::Abstract::Vector& input, 
				     NOX::Abstract::Vector& result) const;

  bool isF() const;
  bool isJacobian() const;
  bool isGradient() const;
  bool isNewton() const;

  const NOX::Abstract::Vector& getX() const;
  const NOX::Abstract::Vector& getF() const;
  double getNormF() const;
  const NOX::Abstract::Vector& getGradient() const;
  const NOX::Abstract::Vector& getNewton() const;

  Teuchos::RefCountPtr<NOX::Abstract::Group> 
    clone(NOX::CopyType type = NOX::DeepCopy) const;

protected:

  // resets the isValid flags to false
  void resetIsValid_();

  // Throws an error
  void throwError(std::string method, std::string message) const;


protected:

  // Reference to the shared Newton system
  SharedSystem* sharedSystemPtr_;

  // NOX Vectors for storing values
  Teuchos::RefCountPtr<N_NLS_NOX::Vector> xVecPtr_;
  N_NLS_NOX::Vector& xVec_;
  Teuchos::RefCountPtr<N_NLS_NOX::Vector> fVecPtr_;
  N_NLS_NOX::Vector& fVec_;
  Teuchos::RefCountPtr<N_NLS_NOX::Vector> newtonVecPtr_;
  Teuchos::RefCountPtr<N_NLS_NOX::Vector> gradVecPtr_;

  // Booleans for tracking whether or not these values have been
  // computed for the currect x.
  bool isValidF_;
  bool isValidJacobian_;
  bool isValidGradient_;
  bool isValidNewton_;
  mutable bool isValidPreconditioner_;

  // Value of the 2-Norm of F
  double normF_;

  // Flag to determine if the solver should refactor the 
  // preconditioner (iterative) or Jacobian (direct).
  mutable bool haveSolverFactors_;

}; // class SharedSystem
} // namespace N_NLS_NOX

#endif // Xyce_N_NLS_NOX_SharedSystem_h


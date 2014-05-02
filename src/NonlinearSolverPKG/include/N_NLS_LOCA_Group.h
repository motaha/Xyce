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
// Filename       : $RCSfile: N_NLS_LOCA_Group.h,v $
//
// Purpose        : Interface to Xyce for LOCA groups.
//
// Special Notes  :
//
// Creator        : Roger Pawlowski, NLS, 9233
//
// Creation Date  : 02/17/03
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.26 $
//
// Revision Date  : $Date: 2014/02/24 23:49:24 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef Xyce_N_NLS_LOCA_Group_h
#define Xyce_N_NLS_LOCA_Group_h

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#ifdef Xyce_PARALLEL_MPI
#include "N_PDS_ParMap.h"
#include "N_PDS_ParComm.h"
#endif
#include <N_IO_fwd.h>
#include <N_UTL_fwd.h>
#include <N_ANP_fwd.h>

// ----------   NOX Includes   ----------

#include "N_NLS_NOX_Group.h"       // base class
#include "LOCA_Abstract_Group.H"   // base class
#include "LOCA_Parameter_Vector.H" // data member
#include "LOCA_DerivUtils.H"       // data member
#include "N_LAS_Vector.h"          // data member
#include "Teuchos_RefCountPtr.hpp" // data member

// ---------- Forward Declarations ----------

namespace N_NLS_NOX {
  class Vector;
  class SharedSystem;
  class AugmentLinSys;
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

class N_LOA_Loader;
class N_LAS_Matrix;

// N_NLS_LOCA namespace is for all LOCA-related classes in the Xyce
// Nonlinear Solver Package
namespace N_NLS_LOCA {


//-----------------------------------------------------------------------------
// Class         : N_NLS::LOCA::Group
//
// Purpose       :
//
//      NOX Group Interface for Xyce
//
// Creator       : Roger Pawlowski, SNL, 9233
//
// Creation Date : 2/17/03
//-----------------------------------------------------------------------------

class Group : public N_NLS_NOX::Group, public LOCA::Abstract::Group {

public:

  //! Basic Constructor
    Group(Teuchos::RefCountPtr<LOCA::GlobalData> globalData,
	  N_NLS_NOX::SharedSystem& s, N_LOA_Loader& l, N_IO_OutputMgr& o,
	  N_ANP_AnalysisInterface & t);

  //! Copy Constructor
  Group(const Group& source, NOX::CopyType type = NOX::DeepCopy);

  //! Destructor
  ~Group();

  //! Assignment Operator
  NOX::Abstract::Group& operator=(const NOX::Abstract::Group& source);

  //! Assignment Operator
  N_NLS_NOX::Group& operator=(const N_NLS_NOX::Group& source);

  //! Assignment Operator
  LOCA::Abstract::Group& operator=(const LOCA::Abstract::Group& source);

  //! Assignment Operator
  N_NLS_LOCA::Group& operator=(const N_NLS_LOCA::Group& source);

  //! Special LOCA assignment operator
  void 	copy (const NOX::Abstract::Group &source);

  //! Cloning function
  Teuchos::RefCountPtr<NOX::Abstract::Group>
    clone(NOX::CopyType type = NOX::DeepCopy) const;

  //! Overloaded function evluation routine
  NOX::Abstract::Group::ReturnType computeF();

  //! Overloaded Jacobian evaluation routine
  NOX::Abstract::Group::ReturnType computeJacobian();

  void setParams(const LOCA::ParameterVector& p);

  const LOCA::ParameterVector& getParams() const;

  void setParam(int paramID, double value);

  double getParam(int paramID) const;

  void setParam(std::string paramID, double value);

  double getParam(std::string paramID) const;

  void setScaleVec(const NOX::Abstract::Vector& s);

  const NOX::Abstract::Vector& getScaleVec() const;

  NOX::Abstract::Group::ReturnType
    augmentJacobianForHomotopy(double conParamValue);

  void printSolution (const double conParam) const;

  void printSolution (const NOX::Abstract::Vector &x,
		      const double conParam) const;

  void stepFailed    ();
  void stepSucceeded ();

  // Pseudo Transient methods
  void setAugmentLinearSystem(bool enable,
		  const Teuchos::RefCountPtr<N_NLS_NOX::AugmentLinSys>& ls);

  // Continuation flag accessors
  void setNonContinuationFlag (bool value);
  bool getNonContinuationFlag ();

  // dcop restart functions
    void setOutputLinear (Xyce::NodeNamePairMap * op,
                          Xyce::NodeNamePairMap * allNodes
#ifdef Xyce_PARALLEL_MPI
       , N_PDS_Comm * pdsCommPtr);
#else
       );
#endif

private:

  // dcop restart data.
  bool outputLinear_;
  int serialNumber_;
  std::map<int, double> oldSol_;
  Xyce::NodeNamePairMap *op_;
  Xyce::NodeNamePairMap *allNodes_;
#ifdef Xyce_PARALLEL_MPI
  N_PDS_Comm * pdsCommPtr_;
#endif

  // dcop restart function
  void outputLinearSystem_ (N_LAS_Matrix* jacobian,
                           N_LAS_Vector* solution,
                           N_LAS_Vector* residual_vector);

  //! Keep a reference to the loader to set parameters.
  N_LOA_Loader& loader;

  //! For output to a file we need xyce's output manager.
  N_IO_OutputMgr& outputMgr;

  //! need xyce's time integration manager.
  N_ANP_AnalysisInterface & anaInt;

  //! Parameter vector container
  LOCA::ParameterVector params;

  //! Utilities for computing derivatives
  LOCA::DerivUtils derivUtils;

  //! Temporary vector used for homotopy calculation
  N_LAS_Vector* tmpVectorPtr;

  //! LOCA Scaling Vector
  const NOX::Abstract::Vector* scalingVecPtr;

  // Objects for Pseudo transient continuation
  bool useAugmentLinSys_;
  Teuchos::RefCountPtr<N_NLS_NOX::AugmentLinSys> augmentLSStrategy_;

  // Flag to indicate if this is a traditional newton solve or not.
  bool nonContinuationSolve_;

}; // class N_NLS_LOCA::Group
} // namespace N_NLS_LOCA

#endif // Xyce_N_NLS_LOCA_Group_h


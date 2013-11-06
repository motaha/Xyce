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
// Filename       : $RCSfile: N_NLS_LOCA_Group.C,v $
//
// Purpose        : Interface to Xyce for LOCA continuation routines.
//
// Special Notes  :
//
// Creator        : Roger Pawlowski, SNL 9233
//
// Creation Date  : 02/17/03
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.44.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:48 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include <N_UTL_Misc.h>

#include "N_NLS_LOCA_Group.h"
//#include "N_NLS_NOX.h"
//#include "N_NLS_NOX_Vector.h"
#include "N_NLS_NOX_SharedSystem.h"
#include "N_NLS_NOX_AugmentLinSys.h"
#include "N_LAS_Vector.h"
#include "N_LAS_Matrix.h"
#include "N_LAS_System.h"
#include "N_LAS_Builder.h"
#include "N_ERH_ErrorMgr.h"
#include "N_LOA_Loader.h"
#include "N_IO_OutputMgr.h"
#include "N_ANP_AnalysisInterface.h"

// ----------   NOX Includes   ----------
#include "Teuchos_ParameterList.hpp"
#include "NOX_Abstract_Vector.H"
#include "LOCA_Parameter_Vector.H"

//using namespace N_NLS_LOCA;
using namespace N_NLS_NOX;

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::Group
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
N_NLS_LOCA::Group::Group(Teuchos::RefCountPtr<LOCA::GlobalData> gd,
			 N_NLS_NOX::SharedSystem& s, N_LOA_Loader& l,
			 N_IO_OutputMgr& o, N_ANP_AnalysisInterface & t) :
  N_NLS_NOX::Group(s),
  LOCA::Abstract::Group(gd),
  loader(l),
  outputMgr(o),
  anaInt(t),
  derivUtils(gd),
  tmpVectorPtr(0),
  scalingVecPtr(0),
  useAugmentLinSys_(false),
  outputLinear_(false),
  serialNumber_(0),
  op_(0),
  allNodes_(0),
#ifdef Xyce_PARALLEL_MPI
  pdsCommPtr_(0),
#endif
  nonContinuationSolve_(true)
{
  tmpVectorPtr = sharedSystemPtr_->getLasSystem()->builder().createVector();
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::Group
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
N_NLS_LOCA::Group::Group(const Group& source, NOX::CopyType type) :
  N_NLS_NOX::Group(source, type),
  LOCA::Abstract::Group(source, type),
  loader(source.loader),
  outputMgr(source.outputMgr),
  anaInt(source.anaInt),
  params(source.params),
  derivUtils(source.derivUtils),
  tmpVectorPtr(0),
  scalingVecPtr(source.scalingVecPtr),
  useAugmentLinSys_(source.useAugmentLinSys_),
  outputLinear_(source.outputLinear_),
  serialNumber_(source.serialNumber_),
  oldSol_(source.oldSol_),
  op_(source.op_),
  allNodes_(source.allNodes_),
#ifdef Xyce_PARALLEL_MPI
  pdsCommPtr_(source.pdsCommPtr_),
#endif
  augmentLSStrategy_(source.augmentLSStrategy_),
  nonContinuationSolve_(source.nonContinuationSolve_)
{
  tmpVectorPtr = sharedSystemPtr_->getLasSystem()->builder().createVector();
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::~Group
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
N_NLS_LOCA::Group::~Group()
{
  delete tmpVectorPtr;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::operator=
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
NOX::Abstract::Group&
N_NLS_LOCA::Group::operator=(const NOX::Abstract::Group& source)
{
  return operator=(dynamic_cast<const N_NLS_LOCA::Group&>(source));
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::operator=
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
N_NLS_NOX::Group&
N_NLS_LOCA::Group::operator=(const N_NLS_NOX::Group& source)
{
  return operator=(dynamic_cast<const N_NLS_LOCA::Group&>(source));
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::operator=
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
LOCA::Abstract::Group&
N_NLS_LOCA::Group::operator=(const LOCA::Abstract::Group& source)
{
  return operator=(dynamic_cast<const N_NLS_LOCA::Group&>(source));
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::operator=
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
N_NLS_LOCA::Group&
N_NLS_LOCA::Group::operator=(const N_NLS_LOCA::Group& source)
{
  N_NLS_NOX::Group::operator=(source);
  params = source.params;
  derivUtils = source.derivUtils;
  if (source.scalingVecPtr != 0)
    scalingVecPtr = source.scalingVecPtr;
  useAugmentLinSys_ = source.useAugmentLinSys_;
  augmentLSStrategy_ = source.augmentLSStrategy_;
  return *this;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::copy
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void N_NLS_LOCA::Group::copy(const NOX::Abstract::Group &source)
{
  *this = source;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::clone
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
Teuchos::RefCountPtr<NOX::Abstract::Group> N_NLS_LOCA::Group::
clone(NOX::CopyType type) const
{
  Teuchos::RefCountPtr<N_NLS_LOCA::Group> ptr =
    Teuchos::rcp(new Group(*this, type));
  return ptr;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::computeF()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
NOX::Abstract::Group::ReturnType N_NLS_LOCA::Group::computeF()
{
  // Set continuation parameters, if neccessary
  // Note:  If the tranOP solve was a continuation, but the
  // transient is traditional Newton, make SURE that the
  // setParam calls are not happening for transient!
  if (!nonContinuationSolve_)
  {
    for (int i = 0; i < params.length(); ++i) {
      string label = params.getLabel(i);
      loader.setParam(label, params.getValue(i));

      if (label == "GSTEPPING" && useAugmentLinSys_)
        augmentLSStrategy_->setProgressVariable(params.getValue(i));

    }
  }

  NOX::Abstract::Group::ReturnType status = N_NLS_NOX::Group::computeF();

  if (useAugmentLinSys_)
    augmentLSStrategy_->augmentResidual(xVec_.getNativeVectorPtr(),
					fVec_.getNativeVectorPtr());

  return status;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::computeJacobian()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
NOX::Abstract::Group::ReturnType N_NLS_LOCA::Group::computeJacobian()
{
  // Set continuation parameters, if neccessary
  // Note:  If the tranOP solve was a continuation, but the
  // transient is traditional Newton, make SURE that the
  // setParam calls are not happening for transient!
  if (!nonContinuationSolve_)
  {
    for (int i = 0; i < params.length(); ++i) {
      string label = params.getLabel(i);
      loader.setParam(label, params.getValue(i));

      if (label == "GSTEPPING" && useAugmentLinSys_)
        augmentLSStrategy_->setProgressVariable(params.getValue(i));

    }
  }

  NOX::Abstract::Group::ReturnType status =
    N_NLS_NOX::Group::computeJacobian();

  // Augment jacobian for pseudo transient if enabled
  if (useAugmentLinSys_) {
    N_LAS_Matrix& jacobian =
      const_cast<N_LAS_Matrix&>(sharedSystemPtr_->getJacobian());
    augmentLSStrategy_->augmentJacobian(&jacobian);
  }

  if (outputLinear_)
  {
    N_LAS_Matrix& jacobian =
      const_cast<N_LAS_Matrix&>(sharedSystemPtr_->getJacobian());
    cout << "After computeJacobian, linear system is:" << endl;
    outputLinearSystem_ (&jacobian, xVec_.getNativeVectorPtr(),
                                   fVec_.getNativeVectorPtr());
  }

#ifdef Xyce_DEBUG_NONLINEAR
  N_LAS_Matrix& jacobian =
    const_cast<N_LAS_Matrix&>(sharedSystemPtr_->getJacobian());
  sharedSystemPtr_->debugOutput1( jacobian,
        (*(fVec_.getNativeVectorPtr())));
#endif

  return status;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::setOutputLinear
// Purpose       :
// Special Notes : for DCOP restart (debugging)
// Scope         : public
// Creator       : Dave Shirley
// Creation Date : 2006
//-----------------------------------------------------------------------------
void N_NLS_LOCA::Group::setOutputLinear (Xyce::NodeNamePairMap * op,
                                         Xyce::NodeNamePairMap * allNodes
#ifdef Xyce_PARALLEL_MPI
       , N_PDS_Comm * pdsCommPtr)
#else
       )
#endif
{
  op_ = op;
  allNodes_ = allNodes;
#ifdef Xyce_PARALLEL_MPI
  pdsCommPtr_ = pdsCommPtr;
#endif
  outputLinear_ = true;
  serialNumber_ = 0;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::outputLinearSystem_
// Purpose       :
// Special Notes : for DCOP restart (debugging)
// Scope         : private
// Creator       : Dave Shirly
// Creation Date : 2006
//-----------------------------------------------------------------------------
void N_NLS_LOCA::Group::outputLinearSystem_ (N_LAS_Matrix* jacobian,
                         N_LAS_Vector* solPtr_,
                         N_LAS_Vector* resPtr_)
{
  Xyce::NodeNamePairMap::iterator op_i;
  Xyce::NodeNamePairMap::iterator op_end;
  int i, row, global_row;
  int num;
  vector<int> col;
  vector<double> val;
  map<int,string> rowOut;
  int rowLen, GID;

  if (!outputLinear_)
    return;

#ifdef Xyce_PARALLEL_MPI
  N_PDS_ParMap * pmap_;
  pmap_ = resPtr_->pmap();
  int procID = pdsCommPtr_->procID();
#endif
  op_i = allNodes_->begin();
  op_end = allNodes_->end();
  for ( ; op_i != op_end ; ++op_i)
  {
   ostringstream s;
    row = (*op_i).second.first;
#ifdef Xyce_PARALLEL_MPI
    global_row = pmap_->localToGlobalIndex(row);
#else
    global_row = row;
#endif
    s << "Global: " << global_row << " : " << (*op_i).first << "  Row: " << global_row;
    s << " Value: " << (*solPtr_)[row];
    if (serialNumber_ > 0)
    {
      s << endl;
      s << "  Delta Value: " << (*solPtr_)[row] - oldSol_[row];
      s << endl;
    }
    oldSol_[row] = (*solPtr_)[row];
    s << " Residual: " << (*resPtr_)[row];
#ifdef Xyce_PARALLEL_MPI
    s << "  proc: " << procID;
#endif
    s << endl;
    rowLen = jacobian->getLocalRowLength(row);
    col.resize(rowLen);
    val.resize(rowLen);
    jacobian->getRowCopy(global_row, rowLen, rowLen, &val[0], &col[0]);
    for (i=0 ; i<rowLen ; i++)
    {
      if (i>1 && i%10 == 0)
      s << endl;
      GID = col[i];
      s << "  " << GID << "(" << val[i] << ")";
    }
    rowOut[global_row] = s.str();
  }
  serialNumber_++;
  map<int,string>::iterator row_i;
  map<int,string>::iterator row_end = rowOut.end();
  row_i = rowOut.begin();
  string str;
#ifdef Xyce_PARALLEL_MPI
  int numG;
  int pos, posG;
  int len;
  string buf;
#endif
  int big=2000000000;
  for ( ; ; ++row_i )
  {
    if (row_i == row_end)
    {
      str = "";
      num = big;
    }
    else
    {
      str = (*row_i).second;
      num = (*row_i).first;
    }
#ifdef Xyce_PARALLEL_MPI
    numG = -1;
    while (numG != num)
    {
      pdsCommPtr_->minAll (&num, &numG, 1);
      if (numG == big)
        break;
      if (num == numG)
        pos = procID;
      else
        pos = 0;
      pdsCommPtr_->sumAll (&pos, &posG, 1);
      if (procID == 0)
      {
        if (posG != 0)
        {
          pdsCommPtr_->recv(&len, 1, posG);
          buf.resize(len);
          pdsCommPtr_->recv(&buf[0], len, posG);
          cout << buf << endl;
        }
        else
          cout << str << endl;
      }
      else if (procID == posG)
      {
        len = str.size();
        pdsCommPtr_->send(&len, 1, 0);
        pdsCommPtr_->send(str.c_str(), len, 0);
      }
    }
    if (numG == big)
      break;
#else
    cout << str << endl;
    if (num == big)
      break;
#endif
  }
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::setParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void N_NLS_LOCA::Group::setParams(const LOCA::ParameterVector& p)
{
  N_NLS_NOX::Group::resetIsValid_();
  params = p;
  return;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::getParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
const LOCA::ParameterVector& N_NLS_LOCA::Group::getParams() const
{
  return params;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::setParam
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void N_NLS_LOCA::Group::setParam(int paramID, double value)
{
  N_NLS_NOX::Group::resetIsValid_();
  params.setValue(paramID, value);
  return;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::getParam
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
double N_NLS_LOCA::Group::getParam(int paramID) const
{
  return params.getValue(paramID);
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::setParam
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void N_NLS_LOCA::Group::setParam(string paramID, double value)
{
  N_NLS_NOX::Group::resetIsValid_();
  params.setValue(paramID, value);
  return;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::getParam
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
double N_NLS_LOCA::Group::getParam(string paramID) const
{
  return params.getValue(paramID);
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::setScaleVec
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void N_NLS_LOCA::Group::setScaleVec(const NOX::Abstract::Vector& s)
{
  scalingVecPtr = &s;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::getScaleVec()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
const NOX::Abstract::Vector& N_NLS_LOCA::Group::getScaleVec() const
{
  if (scalingVecPtr == 0) {
    const string message = "ERROR: N_NLS_LOCA::Group::getScaleVec() - scaling vector not set!";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, message);
  }

  return (*scalingVecPtr);
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::augmentJacobianForHomotopy
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
NOX::Abstract::Group::ReturnType
N_NLS_LOCA::Group::augmentJacobianForHomotopy(double conParamValue)
{

  N_LAS_Matrix& jacobian =
    const_cast<N_LAS_Matrix&>(sharedSystemPtr_->getJacobian());

  //jacobian.printPetraObject();
  jacobian.scale(conParamValue);
  jacobian.getDiagonal(*tmpVectorPtr);
  (*tmpVectorPtr).addScalar(1.0 - conParamValue);
  jacobian.replaceDiagonal(*tmpVectorPtr);
  //jacobian.printPetraObject();

  return NOX::Abstract::Group::Ok;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::printSolution
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void N_NLS_LOCA::Group::printSolution (const double conParam) const
{
  // ERK:  This is a KLUDGE!  I put it here because the "printSolution"
  // functions were the only place in this class that were called
  // only after a continuation step was successful.
  //
  // Note: homotopy output is now called from the time integrator,
  // so it is this function still results in the the solution being
  // printed.
  //
  // Note:  This will be replaced by "stepSuccess" (see below), once
  // stepSuccess completely works.
  anaInt.completeHomotopyStep
        (params.getNamesVector(),
			   params.getValuesVector(),
			   xVec_.getNativeVectorPtr());

  return;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::printSolution
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void N_NLS_LOCA::Group::printSolution (const NOX::Abstract::Vector &x,
				       const double conParam) const
{
  const N_NLS_NOX::Vector& noxVector =
    dynamic_cast<const N_NLS_NOX::Vector&>(x);

  // ERK: This is a KLUDGE!  I put it here because the "printSolution"
  // functions were the only place in this class that were called
  // only after a continuation step was successful.
  //
  // Note: homotopy output is now called from the time integrator,
  // so it is this function still results in the the solution being
  // printed.
  //
  // Note:  This will be replaced by "stepSuccess" (see below), once
  // stepSuccess completely works.
  anaInt.completeHomotopyStep
        (params.getNamesVector(),
			   params.getValuesVector(),
			   xVec_.getNativeVectorPtr());

  return;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::stepFailed
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void N_NLS_LOCA::Group::stepFailed    ()
{
#ifdef Xyce_DEBUG_NONLINEAR
  cout << "In N_NLS_LOCA::Group::stepFailed" << endl;
#endif
  anaInt.failHomotopyStep ();
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::stepSucceeded
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void N_NLS_LOCA::Group::stepSucceeded ()
{
#ifdef Xyce_DEBUG_NONLINEAR
  cout << "In N_NLS_LOCA::Group::stepSucceeded" << endl;
#endif
#ifdef Xyce_UPDATED_LOCA
  // Do nothing.  Unfortunately, this function doesn't get called
  // on the initial solve.  (for continuation param = initial value)
  // So, for now, still have to use the printSolution function, which
  // is called every time.
    const string message = "ERROR: N_NLS_LOCA::Group::stepSucceeded is not fully supported yet!  Recompile without Xyce_UPDATED_LOCA defined\n";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, message);
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::setAugmentLinearSystem
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void N_NLS_LOCA::Group::setAugmentLinearSystem(bool enable_value,
		const Teuchos::RefCountPtr<N_NLS_NOX::AugmentLinSys>& ls)
{
  useAugmentLinSys_ = enable_value;
  augmentLSStrategy_ = ls;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::setNonContinuationFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void N_NLS_LOCA::Group::setNonContinuationFlag (bool value)
{
  nonContinuationSolve_ = value;
}

//-----------------------------------------------------------------------------
// Function      : N_NLS_LOCA::Group::getNonContinuationFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool N_NLS_LOCA::Group::getNonContinuationFlag ()
{
  return nonContinuationSolve_;
}

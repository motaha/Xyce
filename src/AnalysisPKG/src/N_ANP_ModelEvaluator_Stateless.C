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
// Filename      : $RCSfile: N_ANP_ModelEvaluator_Stateless.C,v $
// Purpose       : This file supports the integration of Xyce with
//                 Rythmos.
// Special Notes :
// Creator       : Todd Coffey, 1414
// Creation Date : mid 2009
//
// Revision Information:
// ---------------------

#include <Xyce_config.h>

// ----------   Xyce Includes   ----------
#include <N_DEV_fwd.h>
#undef HAVE_LIBPARMETIS
#include <N_ANP_ModelEvaluator_Stateless.h>
#include <N_LAS_Vector.h>
#include <N_LAS_BlockVector.h>
#include <N_LAS_BlockMatrix.h>
#include <N_PDS_ParMap.h>
// ----------   Trilinos Includes   ----------
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_Operator.h>
#include <Epetra_CrsMatrix.h>

using Teuchos::rcp;

namespace Xyce {
namespace Analysis {

RCP<ModelEvaluator_Stateless> N_ANP_modelEvaluator_Stateless()
{
    return rcp(new ModelEvaluator_Stateless());
}


RCP<ModelEvaluator_Stateless> N_ANP_modelEvaluator_Stateless(
    const RCP<ModelEvaluator>& xyceME)
{
    RCP<ModelEvaluator_Stateless>
    xyceMES = N_ANP_modelEvaluator_Stateless();
    xyceMES->set_XyceModelEvaluator(xyceME);
    return xyceMES;
}


ModelEvaluator_Stateless::ModelEvaluator_Stateless()
  :isInitialized_(false)
{
}


ModelEvaluator_Stateless::~ModelEvaluator_Stateless()
{
}

void ModelEvaluator_Stateless::set_XyceModelEvaluator(
  const RCP<Analysis::ModelEvaluator>& xyceME)
{
    TEUCHOS_ASSERT(xyceME->isInitialized());
    xyceME_ = xyceME;
    tempStateVector_ = rcp(new Epetra_Vector(*xyceME_->get_g_map(0)));
    tempStateDotVector_ = rcp(new Epetra_Vector(*xyceME_->get_g_map(0)));
    tempVoltLimQVector_ = rcp(new Epetra_Vector(*xyceME_->get_x_map()));
    tempVoltLimFVector_ = rcp(new Epetra_Vector(*xyceME_->get_x_map()));
    tempFVector_ = rcp(new Epetra_Vector(*xyceME_->get_f_map()));
    tempWOperator_ = xyceME_->create_W();
    this->setupInOutArgs_();
}

void ModelEvaluator_Stateless::setupInOutArgs_()
{
  if (!isInitialized_) {
    Np_ = 0;
    Ng_ = 0;
    EpetraExt::ModelEvaluator::InArgsSetup inArgs;
    inArgs.setSupports(IN_ARG_t,true);
    inArgs.setSupports(IN_ARG_x,true);
    inArgs.setSupports(IN_ARG_x_dot,true);
    inArgs.setSupports(IN_ARG_alpha,true);
    inArgs.setSupports(IN_ARG_beta,true);
    inArgs.set_Np(Np_);
    inArgs_ = inArgs;

    EpetraExt::ModelEvaluator::OutArgsSetup outArgs;
    outArgs.setSupports(OUT_ARG_f,true);
    outArgs.setSupports(OUT_ARG_W,true);
    outArgs.set_Np_Ng(Np_,Ng_);
    outArgs_ = outArgs;

    isInitialized_ = true;
  }
}


EpetraExt::ModelEvaluator::InArgs ModelEvaluator_Stateless::createInArgs() const
{
  TEUCHOS_ASSERT(isInitialized_);
  return inArgs_;
}


EpetraExt::ModelEvaluator::OutArgs ModelEvaluator_Stateless::createOutArgs() const
{
  TEUCHOS_ASSERT(isInitialized_);
  return outArgs_;
}


RCP<const Epetra_Map> ModelEvaluator_Stateless::get_x_map() const
{
  TEUCHOS_ASSERT(isInitialized_);
  return xyceME_->get_x_map();
}


RCP<const Epetra_Map> ModelEvaluator_Stateless::get_f_map() const
{
  TEUCHOS_ASSERT(isInitialized_);
  return xyceME_->get_f_map();
}


RCP<const Epetra_Map> ModelEvaluator_Stateless::get_p_map(int p) const
{
  TEUCHOS_ASSERT(isInitialized_);
  return Teuchos::null;
}


RCP<const Epetra_Map> ModelEvaluator_Stateless::get_g_map(int p) const
{
  TEUCHOS_ASSERT(isInitialized_);
  return Teuchos::null;
}


void ModelEvaluator_Stateless::evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const
{

  TEUCHOS_ASSERT(isInitialized_);
  double t = inArgs.get_t();
  RCP<const Epetra_Vector> x = inArgs.get_x().assert_not_null();
  RCP<const Epetra_Vector> x_dot = inArgs.get_x_dot().assert_not_null();
  RCP<Epetra_Vector> f = outArgs.get_f();
  RCP<Epetra_Operator> W = outArgs.get_W();

  // First we evaluate the state vector:
  {
    EpetraExt::ModelEvaluator::InArgs xyceInArgs = xyceME_->createInArgs();
    EpetraExt::ModelEvaluator::OutArgs xyceOutArgs = xyceME_->createOutArgs();

    xyceInArgs.set_t(t);
    xyceInArgs.set_x(x);
    xyceInArgs.set_x_dot(x_dot);
    xyceOutArgs.set_g(0,tempStateVector_);
    xyceME_->evalModel(xyceInArgs,xyceOutArgs);
  }

  {
    // Fabricate a StateDot Vector
    // tempStateDotVector_ = tempStateVector_
    tempStateDotVector_->Update(1.0,*tempStateVector_,0.0);
  }

  if (Teuchos::is_null(f)) {
    // If the client did not ask for f, then we have to load it into a
    // temporary vector because Xyce does not support loading f
    // separately from W.
    f = tempFVector_;
  }
  if (Teuchos::is_null(W)) {
    // If the client did not ask for W, then we have to load it into a
    // temporary vector because Xyce does not support loading f
    // separately from W.
    W = tempWOperator_;
  }

  // Second, we evaluate f and W:
  {
    EpetraExt::ModelEvaluator::InArgs xyceInArgs = xyceME_->createInArgs();
    EpetraExt::ModelEvaluator::OutArgs xyceOutArgs = xyceME_->createOutArgs();
    xyceInArgs.set_t(t);
    xyceInArgs.set_x(x);
    xyceInArgs.set_x_dot(x_dot);
    xyceInArgs.set_alpha(inArgs.get_alpha());
    xyceInArgs.set_beta(inArgs.get_beta());
    xyceInArgs.set_p(0,tempStateVector_);
    xyceInArgs.set_p(1,tempStateDotVector_);
    xyceOutArgs.set_g(1,tempVoltLimQVector_);
    xyceOutArgs.set_g(2,tempVoltLimFVector_);
    xyceOutArgs.set_f(f);
    xyceOutArgs.set_W(W);
    xyceME_->evalModel(xyceInArgs,xyceOutArgs);
  }
}


RCP<Epetra_Operator> ModelEvaluator_Stateless::create_W() const
{
  TEUCHOS_ASSERT(isInitialized_);
  return xyceME_->create_W();
}

} // namespace Analysis
} // namespace Xyce

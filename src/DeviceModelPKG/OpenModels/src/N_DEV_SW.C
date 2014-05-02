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

//----------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_SW.C,v $
//
// Purpose        :
//
// Special Notes  :  "Creator and Creation Date" are actually the dates
//                   this file was created.  When it was created it was
//                   just a placeholder.
//                   Actual implementation of the Xyce voltage controlled
//                   switch occurred on 5/22/2001, and was done by Tom Russo,
//                   SNL, Component Information and Models.
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 02/28/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.157.2.2 $
//
// Revision Date  : $Date: 2014/03/06 23:33:43 $
//
// Current Owner  : $Author: tvrusso $
//----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_SW.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_Message.h>
#include <N_ERH_ErrorMgr.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>

#include <N_UTL_Expression.h>

namespace Xyce {
namespace Device {


namespace SW {


void Traits::loadInstanceParameters(ParametricData<SW::Instance> &p)
{
// Set up double precision variables:
    p.addPar ("CONTROL", 0.0, false, ParameterType::SOLN_DEP,
        &SW::Instance::CONTROL,
        NULL, U_NONE,CAT_NONE,"");

    // Set up non-double precision variables:
    p.addPar ("ON", false,  &SW::Instance::ON);
    p.addPar ("OFF", false, &SW::Instance::OFF);
}

void Traits::loadModelParameters(ParametricData<SW::Model> &p)
{
    p.addPar ("RON",      1.0,   false, ParameterType::NO_DEP,
        &SW::Model::RON,
        NULL,U_OHM,CAT_NONE,"On resistance");

    p.addPar ("ROFF",     1.0e6, false, ParameterType::NO_DEP,
        &SW::Model::ROFF,
        NULL,U_OHM,CAT_NONE,"Off resistance");

    p.addPar ("VON",      1.0,   false, ParameterType::NO_DEP,
        &SW::Model::VON,
        NULL,U_VOLT,CAT_NONE,"On voltage");

    p.addPar ("VOFF",     0.0,   false, ParameterType::NO_DEP,
        &SW::Model::VOFF,
        NULL,U_VOLT,CAT_NONE,"Off voltage");

    p.addPar ("ION",      0.001, false, ParameterType::NO_DEP,
        &SW::Model::ION,
        NULL,U_AMP,CAT_NONE,"On current");

    p.addPar ("IOFF",     0.0,   false, ParameterType::NO_DEP,
        &SW::Model::IOFF,
        NULL,U_AMP,CAT_NONE,"Off current");

    p.addPar ("ON",       1.0,   false, ParameterType::NO_DEP,
        &SW::Model::ON,
        NULL,U_NONE,CAT_NONE,"On control value");

    p.addPar ("OFF",      0.0,   false, ParameterType::NO_DEP,
        &SW::Model::OFF,
        NULL,U_NONE,CAT_NONE,"Off control value");
}


// Class Instance
//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 03/16/05
//-----------------------------------------------------------------------------

bool Instance::processParams ()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock & IB,
  Model & SWiter,
  const FactoryBlock &  factory_block)
  : DeviceInstance(IB, configuration.getInstanceParameters(), factory_block),
    model_(SWiter),
    R(0.0),
    ON(false),
    G(0.0),
    SW_STATE(0.0),
    switch_state(0.0),
    li_switch_state(-1),
    li_Pos(-1),
    li_Neg(-1),
    li_store_dev_i(-1),
    APosEquPosNodeOffset(-1),
    APosEquNegNodeOffset(-1),
    ANegEquPosNodeOffset(-1),
    ANegEquNegNodeOffset(-1),
    fPosEquPosNodePtr(0),
    fPosEquNegNodePtr(0),
    fNegEquPosNodePtr(0),
    fNegEquNegNodePtr(0),
    Exp_ptr(0)
{
  numIntVars = 0;
  numExtVars = 2;
  numStateVars = 1;
  numLeadCurrentStoreVars = 1; // for device lead current if needed.

  jacStamp.resize(2);
  jacStamp[0].resize(2);
  jacStamp[0][0]=0;
  jacStamp[0][1]=1;
  jacStamp[1].resize(2);
  jacStamp[1][0]=0;
  jacStamp[1][1]=1;

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  // Set any non-constant parameter defaults:
  if (given("OFF"))
  {
    if (given("ON"))
    {
      UserError0(*this) << "Cannot specify both 'on' and off' for switch";
    }
    ON = !OFF;
  }
  if (!given("CONTROL"))
  {
    UserError0(*this) << "Must specify 'control' for switch";
  }

  std::vector<sDepend>::iterator d;
  std::vector<sDepend>::iterator begin = dependentParams.begin();
  std::vector<sDepend>::iterator end = dependentParams.end();

  for  (d = begin ; d != end ; ++d)
  {
    if (d->name == "CONTROL")
    {
      expNumVars = d->n_vars;
      expBaseVar = d->lo_var;
      Exp_ptr = d->expr;

      expNumDdt = Exp_ptr->getNumDdt();
      ddtVals.resize(expNumDdt);
      li_ddt.resize(expNumDdt);
      numStateVars += expNumDdt;

      jacStamp[0].resize(2+expNumVars);
      jacStamp[1].resize(2+expNumVars);
      for( int i = 0; i < expNumVars; ++i )
      {
        jacStamp[0][2+i] = 2+i;
        jacStamp[1][2+i] = 2+i;
      }
      expVarDerivs.resize(expNumVars);
      myVarVals.resize(expNumVars);
    }
  }

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:
  processParams ();
}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
Instance::~Instance ()
{
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/20/02
//-----------------------------------------------------------------------------
void Instance::registerLIDs( const std::vector<int> & intLIDVecRef,
                                     const std::vector<int> & extLIDVecRef)
{
  AssertLIDs(intLIDVecRef.size() == numIntVars);
  AssertLIDs(extLIDVecRef.size() == numExtVars);

  // copy over the global ID lists.
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  // now use these lists to obtain the indices into the
  // linear algebra entities.  This assumes an order.
  // For the matrix  indices, first do the rows.

  li_Pos = extLIDVec[0];
  li_Neg = extLIDVec[1];

}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       : Note that the SW does not have any state vars.
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/20/02
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs( const std::vector<int> & staLIDVecRef )
{
  AssertLIDs(staLIDVecRef.size() == numStateVars);

  // copy over the global ID lists:
  staLIDVec = staLIDVecRef;
  li_switch_state = staLIDVec[0];
  for (int i=0 ; i<expNumDdt ; ++i)
  {
    li_ddt[i] = staLIDVecRef[i+1];
  }
}


//-----------------------------------------------------------------------------
// Function      : N_DEV_SWInstance::registerStoreLIDs
// Purpose       : One store var for device current.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 04/05/2013
//-----------------------------------------------------------------------------
void N_DEV_SWInstance::registerStoreLIDs(const std::vector<int> & stoLIDVecRef )
{
  AssertLIDs(stoLIDVecRef.size() == getNumStoreVars());

  if( loadLeadCurrent )
  {
    li_store_dev_i = stoLIDVecRef[0];
  }
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_SWInstance::getStoreNameMap
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 04/05/2013
//-----------------------------------------------------------------------------
std::map<int,std::string> & N_DEV_SWInstance::getStoreNameMap ()
{
  // set up the internal name map, if it hasn't been already.
  if( loadLeadCurrent && storeNameMap.empty ())
  {
    // change subcircuitname:devicetype_deviceName to
    // devicetype:subcircuitName:deviceName
    std::string modName(getName());
    spiceInternalName(modName);
    std::string tmpstr;
    tmpstr = modName+":DEV_I";
    storeNameMap[ li_store_dev_i ] = tmpstr;
  }

  return storeNameMap;
}


//-----------------------------------------------------------------------------
// Function      : Instance::getDepSolnVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/06/01
//-----------------------------------------------------------------------------
const std::vector<std::string> & Instance::getDepSolnVars()
{
  return DeviceInstance::getDepSolnVars();
}


//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 09/2/02
//-----------------------------------------------------------------------------
const std::vector< std::vector<int> > & Instance::jacobianStamp() const
{
  return jacStamp;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerJacLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 09/2/02
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );

  APosEquPosNodeOffset = jacLIDVec[0][0];
  APosEquNegNodeOffset = jacLIDVec[0][1];
  ANegEquPosNodeOffset = jacLIDVec[1][0];
  ANegEquNegNodeOffset = jacLIDVec[1][1];

  APosEquControlNodeOffset.resize( expNumVars );
  ANegEquControlNodeOffset.resize( expNumVars );
  for( int i = 0; i < expNumVars; ++i )
  {
    APosEquControlNodeOffset[i] = jacLIDVec[0][2+i];
    ANegEquControlNodeOffset[i] = jacLIDVec[1][2+i];
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupPointers
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/30/08
//-----------------------------------------------------------------------------
void Instance::setupPointers ()
{
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  N_LAS_Matrix & dFdx = *(extData.dFdxMatrixPtr);
  N_LAS_Matrix & dQdx = *(extData.dQdxMatrixPtr);

  fPosEquPosNodePtr = &(dFdx[li_Pos][APosEquPosNodeOffset]);
  fPosEquNegNodePtr = &(dFdx[li_Pos][APosEquNegNodeOffset]);
  fNegEquPosNodePtr = &(dFdx[li_Neg][ANegEquPosNodeOffset]);
  fNegEquNegNodePtr = &(dFdx[li_Neg][ANegEquNegNodeOffset]);

  fPosEquControlNodePtr.resize( expNumVars );
  fNegEquControlNodePtr.resize( expNumVars );
  for( int i = 0; i < expNumVars; ++i )
  {
    fPosEquControlNodePtr[i] = &(dFdx[li_Pos][APosEquControlNodeOffset[i]]);
    fNegEquControlNodePtr[i] = &(dFdx[li_Neg][ANegEquControlNodeOffset[i]]);
  }
#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/29/01
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState ()
{
  double * staVec = extData.nextStaVectorRawPtr;
  bool bsuccess = updateIntermediateVars ();

  //  obtain the current value of the switch state
  switch_state = SW_STATE;

  staVec[li_switch_state] = switch_state;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/29/01
//-----------------------------------------------------------------------------
bool Instance::updateSecondaryState ()
{
  double current_state;
  double control;
  int i;

  double * solVec = extData.nextSolVectorRawPtr;

  // Get time derivatives from time integrator, and evaluate expression to get
  // derivatives with respect to independent quantities

  if (expNumDdt > 0)
  {
    double * staDerivVec = extData.nextStaDerivVectorRawPtr;

    for (i=0 ; i<expNumDdt ; ++i)
    {
      ddtVals[i] = staDerivVec[li_ddt[i]];
    }
    Exp_ptr->setDdtDerivs(ddtVals);
  }
  // Evaluate Expression with corrected time derivative values

  Exp_ptr->evaluate( expVal, expVarDerivs );
  control = expVal;

  // This is not really correct, an interim hack.  This is supposed
  // to be where we deal with the specification of ON or OFF from the
  // netlist.
  if (getSolverState().initJctFlag)
  {
    if (ON)
      current_state = 1;
    else
      current_state = 0;
  }
  else
  {
    current_state = (control-model_.OFF)*model_.dInv;
  }

  v_pos = solVec[li_Pos];
  v_neg = solVec[li_Neg];

  if (current_state >= 1.0)
  {
    R = model_.RON;
    G = 1.0/R;
    for (i=0 ; i<expNumVars ; ++i)
      expVarDerivs[i] = 0;
  }
  else if ( current_state <= 0.0)
  {
    R = model_.ROFF;
    G = 1.0/R;
    for (i=0 ; i<expNumVars ; ++i)
      expVarDerivs[i] = 0;
  }
  else
  {
    current_state = 2*current_state - 1;
    G = exp(-model_.Lm - 0.75*model_.Lr*current_state +
		0.25*model_.Lr*current_state*current_state*current_state);
    R = 1.0/G;
    for (i=0 ; i<expNumVars ; ++i)
    {
      expVarDerivs[i] = G * (1.5 * (current_state*current_state-1) * model_.Lr *
                                model_.dInv * expVarDerivs[i]);
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       : update intermediate variables for one switch instance
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, Component Information and Models
// Creation Date : 1/10/01
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars ()
{
  // Get values of the arguments for ddt() calls in expression so that the derivatives
  // can be determined by the time integration class
  if (expNumDdt > 0)
  {
    double * staVec = extData.nextStaVectorRawPtr;

    Exp_ptr->getDdtVals (ddtVals);
    for (int i=0 ; i<expNumDdt ; ++i)
    {
      staVec[li_ddt[i]] = ddtVals[i];
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 resistor  instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 11/13/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  double * fVec = extData.daeFVectorRawPtr;

  // load RHS vector element for the positive circuit node KCL equ.
  double coef = (v_pos-v_neg)*G;

  fVec[li_Pos] += coef;
  fVec[li_Neg] += -coef;
  if( loadLeadCurrent )
  {
    double * stoVec = extData.nextStoVectorRawPtr;
    stoVec[li_store_dev_i] = coef;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 resistor  instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 11/13/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  N_LAS_Matrix & dFdx = *(extData.dFdxMatrixPtr);

  dFdx[li_Pos][APosEquPosNodeOffset] += G;
  dFdx[li_Pos][APosEquNegNodeOffset] -= G;
  dFdx[li_Neg][ANegEquPosNodeOffset] -= G;
  dFdx[li_Neg][ANegEquNegNodeOffset] += G;

  if( expNumVars )
  {
    for( int i = 0; i < expNumVars; ++i )
    {
      dFdx[li_Pos][APosEquControlNodeOffset[i]] += (v_pos-v_neg) * expVarDerivs[i];
      dFdx[li_Neg][ANegEquControlNodeOffset[i]] -= (v_pos-v_neg) * expVarDerivs[i];
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
Model::Model(
  const Configuration & configuration,
  const ModelBlock &    MB,
  const FactoryBlock &  factory_block)
  : DeviceModel(MB, configuration.getModelParameters(), factory_block),
  dtype(1),
  RON(0.0),
  ROFF(0.0),
  ON(0.0),
  OFF(0.0)
{
  if (getType() != "")
  {
    if (getType() == "SWITCH" ) {
      dtype = 1;
    }
    else if (getType() == "ISWITCH") {
      dtype = 2;
    }
    else if (getType() == "VSWITCH") {
      dtype = 3;
    }
    else
    {
      UserError0(*this) << "Unrecognized model type " << getType();
    }
  }


  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to .model line and constant defaults from metadata:
  setModParams (MB.params);

  // Set any non-constant parameter defaults:

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:

  if (dtype == 2)
  {
    if (!given("ON"))
      ON = ION;
    if (!given("OFF"))
      OFF = IOFF;
  }
  else if (dtype == 3)
  {
    if (!given("ON"))
      ON = VON;
    if (!given("OFF"))
      OFF = VOFF;
  }

  processParams();
}

//-----------------------------------------------------------------------------
// Function      : Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 03/16/05
//-----------------------------------------------------------------------------
bool Model::processParams ()
{
  double del;
  Lm = log (sqrt(RON*ROFF));
  Lr = log (RON/ROFF);

  del = ON-OFF;

  if (del < 0 && del > -1e-12)
    del = -1e-12;
  if (del >= 0 && del < 1e-12)
    del = 1e-12;
  dInv = 1.0/del;

  return true;
}

//----------------------------------------------------------------------------
// Function      : Model::processInstanceParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirely, PSSI
// Creation Date : 03/23/06
//----------------------------------------------------------------------------
bool Model::processInstanceParams()
{

  std::vector<Instance*>::iterator iter;
  std::vector<Instance*>::iterator first = instanceContainer.begin();
  std::vector<Instance*>::iterator last  = instanceContainer.end();

  for (iter=first; iter!=last; ++iter)
  {
    (*iter)->processParams();
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Model::~Model
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
Model::~Model ()
{
  std::vector<Instance*>::iterator iter;
  std::vector<Instance*>::iterator first = instanceContainer.begin();
  std::vector<Instance*>::iterator last  = instanceContainer.end();

  for (iter=first; iter!=last; ++iter)
    delete (*iter);
}

//-----------------------------------------------------------------------------
// Function      : Model::printOutInstances
// Purpose       : debugging tool.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/03/00
//-----------------------------------------------------------------------------
std::ostream &Model::printOutInstances(std::ostream &os) const
{
  std::vector<Instance*>::const_iterator iter;
  std::vector<Instance*>::const_iterator first = instanceContainer.begin();
  std::vector<Instance*>::const_iterator last  = instanceContainer.end();

  int i;
  os << std::endl;
  os << "    name     model name  Parameters" << std::endl;
  for (i=0, iter=first; iter!=last; ++iter, ++i)
  {
    os << "  " << i << ": " << (*iter)->getName() << "      ";
    os << getName();
    os << "    R = " << (*iter)->R;
    os << "  G = " << (*iter)->G;
    os << "  State = " << (*iter)->SW_STATE;
    os << std::endl;
  }

  os << std::endl;

  return os;
}

//-----------------------------------------------------------------------------
// Function      : Model::forEachInstance
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : David Baur
// Creation Date : 2/4/2014
//-----------------------------------------------------------------------------
/// Apply a device instance "op" to all instances associated with this
/// model
/// 
/// @param[in] op Operator to apply to all instances.
/// 
/// 
void Model::forEachInstance(DeviceInstanceOp &op) const /* override */ 
{
  for (std::vector<Instance *>::const_iterator it = instanceContainer.begin(); it != instanceContainer.end(); ++it)
    op(*it);
}


//-----------------------------------------------------------------------------
// SW Master functions:
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Function      : Master::updateState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
bool Master::updateState (double * solVec, double * staVec, double * stoVec)
{
  bool bsuccess = true;

  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & si = *(*it);

    bool btmp = si.updateIntermediateVars ();
    bsuccess = bsuccess && btmp;

    //  obtain the current value of the switch state
    si.switch_state = si.SW_STATE;
    staVec[si.li_switch_state] = si.switch_state;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Master::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
bool Master::updateSecondaryState ( double * staDerivVec, double * stoVec )
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & si = *(*it);

    double current_state;
    double control;
    double * solVec = si.extData.nextSolVectorRawPtr;

    // Get time derivatives from time integrator, and evaluate expression to get
    // derivatives with respect to independent quantities

    if (si.expNumDdt > 0)
    {
      for (int i=0 ; i<si.expNumDdt ; ++i)
      {
        si.ddtVals[i] = staDerivVec[si.li_ddt[i]];
      }
      si.Exp_ptr->setDdtDerivs(si.ddtVals);
    }
    // Evaluate Expression with corrected time derivative values

    si.Exp_ptr->evaluate( si.expVal, si.expVarDerivs );
    control = si.expVal;

    // This is not really correct, an interim hack.  This is supposed
    // to be where we deal with the specification of ON or OFF from the
    // netlist.
    if (getSolverState().initJctFlag)
    {
      if (si.ON)
        current_state = 1;
      else
        current_state = 0;
    }
    else
    {
      current_state = (control-si.getModel().OFF)*si.getModel().dInv;
    }

    si.v_pos = solVec[si.li_Pos];
    si.v_neg = solVec[si.li_Neg];

    if (current_state >= 1.0)
    {
      si.R = si.getModel().RON;
      si.G = 1.0/si.R;
      for (int i=0 ; i<si.expNumVars ; ++i)
        si.expVarDerivs[i] = 0;
    }
    else if ( current_state <= 0.0)
    {
      si.R = si.getModel().ROFF;
      si.G = 1.0/si.R;
      for (int i=0 ; i<si.expNumVars ; ++i)
        si.expVarDerivs[i] = 0;
    }
    else
    {
      current_state = 2*current_state - 1;
      si.G = exp(-si.getModel().Lm - 0.75*si.getModel().Lr*current_state +
      0.25*si.getModel().Lr*current_state*current_state*current_state);
      si.R = 1.0/si.G;
      for (int i=0 ; i<si.expNumVars ; ++i)
      {
        si.expVarDerivs[i] = si.G * (1.5 * (current_state*current_state-1) * si.getModel().Lr *
                                  si.getModel().dInv * si.expVarDerivs[i]);
      }
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Master::loadDAEVectors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
bool Master::loadDAEVectors (double * solVec, double * fVec, double *qVec,  double * storeLeadF, double * storeLeadQ)
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
     Instance & si = *(*it);
    // F-vector:
    double coef = (si.v_pos-si.v_neg)*si.G;
    fVec[si.li_Pos] += coef;
    fVec[si.li_Neg] += -coef;
    if( si.loadLeadCurrent )
    {
      storeLeadF[si.li_store_dev_i] = coef;
    }
    // Q-vector:
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Master::loadDAEMatrices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
bool Master::loadDAEMatrices (N_LAS_Matrix & dFdx, N_LAS_Matrix & dQdx)
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
     Instance & si = *(*it);
#ifndef Xyce_NONPOINTER_MATRIX_LOAD

    *si.fPosEquPosNodePtr += si.G;
    *si.fPosEquNegNodePtr -= si.G;
    *si.fNegEquPosNodePtr -= si.G;
    *si.fNegEquNegNodePtr += si.G;

    if( si.expNumVars )
    {
      for( int j = 0; j < si.expNumVars; ++j )
      {
        *si.fPosEquControlNodePtr[j] += (si.v_pos-si.v_neg) * si.expVarDerivs[j];
        *si.fNegEquControlNodePtr[j] -= (si.v_pos-si.v_neg) * si.expVarDerivs[j];
      }
    }
#else

    dFdx[si.li_Pos][si.APosEquPosNodeOffset] += si.G;
    dFdx[si.li_Pos][si.APosEquNegNodeOffset] -= si.G;
    dFdx[si.li_Neg][si.ANegEquPosNodeOffset] -= si.G;
    dFdx[si.li_Neg][si.ANegEquNegNodeOffset] += si.G;

    if( si.expNumVars )
    {
      for( int i = 0; i < si.expNumVars; ++i )
      {
        dFdx[si.li_Pos][si.APosEquControlNodeOffset[i]] += (si.v_pos-si.v_neg) * si.expVarDerivs[i];
        dFdx[si.li_Neg][si.ANegEquControlNodeOffset[i]] -= (si.v_pos-si.v_neg) * si.expVarDerivs[i];
      }
    }
#endif
  }
  return true;
}

Device *Traits::factory(const Configuration &configuration, const FactoryBlock &factory_block)
{

  return new Master(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);
}

void registerDevice()
{
  Config<Traits>::addConfiguration()
    .registerDevice("s", 1)

#ifdef Xyce_OLD_SWITCH
    .registerModelType("sw", 1)
#else
    .registerModelType("switch", 1)
    .registerModelType("iswitch", 1)
    .registerModelType("vswitch", 1)
#endif
    ;
}

} // namespace SW
} // namespace Device
} // namespace Xyce

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
// Filename       : $RCSfile: N_DEV_Inductor.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 02/28/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.269.2.4 $
//
// Revision Date  : $Date: 2014/03/06 23:33:43 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------
#include <Xyce_config.h>

// ---------- Standard Includes ----------
#include <N_UTL_Misc.h>

#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif

#include <algorithm>
#include <fstream>
#include <set>
#include <vector>

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_Inductor.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_Message.h>
#include <N_ERH_ErrorMgr.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>

namespace Xyce {
namespace Device {


namespace Inductor {


void Traits::loadInstanceParameters(ParametricData<Inductor::Instance> &p)
{
  p.addPar("L",    0.0, &Inductor::Instance::baseL)
    .setExpressionAccess(ParameterType::TIME_DEP)
    .setUnit(U_HENRY)
    .setDescription("Inductance");

  p.addPar("IC",   0.0, &Inductor::Instance::IC)
    .setGivenMember(&Inductor::Instance::ICGiven)
    .setUnit(U_AMP)
    .setDescription("Initial current through device");

  p.addPar("TEMP", 0.0, &Inductor::Instance::temp)
    .setExpressionAccess(ParameterType::TIME_DEP)
    .setGivenMember(&Inductor::Instance::tempGiven)
    .setUnit(U_DEGC)
    .setCategory(CAT_MATERIAL)
    .setDescription("Device temperature");

  p.addPar("TC1", 0.0, &Inductor::Instance::tempCoeff1)
    .setGivenMember(&Inductor::Instance::tempCoeff1Given)
    .setUnit(U_DEGCM1)
    .setDescription("Linear Temperature Coefficient");

  p.addPar("TC2", 0.0, &Inductor::Instance::tempCoeff2)
    .setGivenMember(&Inductor::Instance::tempCoeff2Given)
    .setUnit(U_DEGCM2)
    .setDescription("Quadratic Temperature Coefficient");

  // This call tells the parameter handling code that TC can be specified
  // as a vector with up to two elements as in TC=a,b.  It then translates
  // TC=a,b into TC1=a TC2=b.  Likewise, TC=a will translate into TC1=a
  p.makeVector ("TC", 2);
}

void Traits::loadModelParameters(ParametricData<Inductor::Model> &p)
{
  // Set up double precision variables:
  p.addPar("L", 1.0, &Inductor::Model::L)
    .setDescription("Inductance Multiplier");

  p.addPar("IC", 0.0, &Inductor::Model::IC)
    .setUnit(U_AMP)
    .setDescription("Initial current through device");

  p.addPar("TNOM", 27.0, &Inductor::Model::tnom)
    .setUnit(U_DEGC)
    .setCategory(CAT_MATERIAL)
    .setDescription("Reference temperature");

  p.addPar("TC1",0.0, &Inductor::Model::tempCoeff1)
    .setUnit(U_DEGCM1)
    .setCategory(CAT_MATERIAL)
    .setDescription("First order temperature coeff.");

  p.addPar("TC2", 0.0, &Inductor::Model::tempCoeff2)
    .setUnit(U_DEGCM2)
    .setCategory(CAT_MATERIAL)
    .setDescription("Second order temperature coeff.");
}

//
// static class member inits
//
std::vector< std::vector<int> > Instance::jacStamp_BASE;

// Class Instance
//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/03/02
//-----------------------------------------------------------------------------
bool Instance::processParams()
{
  // If there are any time dependent parameters, set their values at for
  // the current time.

  // now set the temperature related stuff.
  updateTemperature(temp);
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateTemperature
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, Component Information and Models
// Creation Date : 02/27/01
//-----------------------------------------------------------------------------
bool Instance::updateTemperature ( const double & temp)
{
  double difference = temp - model_.tnom;
  double factor = model_.L*(1.0 + tempCoeff1*difference +
                         tempCoeff2*difference*difference);
  L = baseL*factor;
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
  const InstanceBlock & instance_block,
  Model &               model,
  const FactoryBlock &  factory_block)
  : DeviceInstance(instance_block, configuration.getInstanceParameters(), factory_block),
    L(0),
    IC(0),
    ICGiven(false),
    model_(model),
    baseL(0.0),
    temp(getDeviceOptions().temp.getImmutableValue<double>()),
    tempGiven(0),
    tempCoeff1(0.0),
    tempCoeff2(0.0),
    tempCoeff1Given(false),
    tempCoeff2Given(false),
    li_fstate(-1),
    li_Pos(-1),
    li_Neg(-1),
    li_Bra(-1),
    ABraEquPosNodeOffset(-1),
    ABraEquNegNodeOffset(-1),
    ABraEquBraVarOffset(-1),
    APosEquBraVarOffset(-1),
    ANegEquBraVarOffset(-1)
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  ,
    fPosEquBraVarPtr(0),
    fNegEquBraVarPtr(0),
    fBraEquPosNodePtr(0),
    fBraEquNegNodePtr(0),
    fBraEquBraVarPtr(0),
    qBraEquBraVarPtr(0)
#endif
{
  numExtVars   = 2;
  numIntVars   = 1;
  numStateVars = 1;

  if( jacStamp_BASE.empty() )
  {
    jacStamp_BASE.resize(3);

    jacStamp_BASE[0].resize(1);
    jacStamp_BASE[0][0] = 2;

    jacStamp_BASE[1].resize(1);
    jacStamp_BASE[1][0] = 2;

    jacStamp_BASE[2].resize(3);
    jacStamp_BASE[2][0] = 0;
    jacStamp_BASE[2][1] = 1;
    jacStamp_BASE[2][2] = 2;
  }

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (instance_block.params);

  // Set any non-constant parameter defaults:
  if (!given("L"))
  {
    UserError0(*this) << "Could not find L parameter in instance.";
  }

  if (!given("TEMP"))
    temp = getDeviceOptions().temp.getImmutableValue<double>();

  if (!tempCoeff1Given)
    tempCoeff1=model_.tempCoeff1;
  if (!tempCoeff2Given)
    tempCoeff2=model_.tempCoeff2;

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:
  processParams ();

  // set up numIntVars:
  numIntVars = 1;

  // set up numStateVars:
  numStateVars = 2;
}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
Instance::~Instance()
{
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

  fPosEquBraVarPtr  = &(dFdx[li_Pos][APosEquBraVarOffset]);
  fNegEquBraVarPtr  = &(dFdx[li_Neg][ANegEquBraVarOffset]);
  fBraEquPosNodePtr = &(dFdx[li_Bra][ABraEquPosNodeOffset]);
  fBraEquNegNodePtr = &(dFdx[li_Bra][ABraEquNegNodeOffset]);
  fBraEquBraVarPtr  = &(dFdx[li_Bra][ABraEquBraVarOffset]);

  qBraEquBraVarPtr = &(dQdx[li_Bra][ABraEquBraVarOffset]);

#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/21/02
//-----------------------------------------------------------------------------
void Instance::registerLIDs(const std::vector<int> & intLIDVecRef,
                                          const std::vector<int> & extLIDVecRef)
{
  AssertLIDs(intLIDVecRef.size() == numIntVars);
  AssertLIDs(extLIDVecRef.size() == numExtVars);

  // copy over the global ID lists.
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  // Now use these lists to obtain the indices into the
  // linear algebra entities.  This assumes an order.
  // For the matrix  indices, first do the rows.

  li_Pos = extLIDVec[0];
  li_Neg = extLIDVec[1];
  li_Bra = intLIDVec[0];

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    Xyce::dout() << section_divider << std::endl;

    Xyce::dout() << "::registerLIDs:\n";
    Xyce::dout() << "  name = " << getName() << std::endl;

    Xyce::dout() << "\nlocal solution indices:\n";
    Xyce::dout() << "  li_Pos = "<< li_Pos << std::endl;
    Xyce::dout() << "  li_Neg = "<< li_Neg << std::endl;
    Xyce::dout() << "  li_Bra = "<< li_Bra << std::endl;

    Xyce::dout() << section_divider << std::endl;
  }
#endif

}

//-----------------------------------------------------------------------------
// Function      : Instance::getIntNameMap
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/13/05
//-----------------------------------------------------------------------------
std::map<int,std::string> & Instance::getIntNameMap ()
{
  // set up the internal name map, if it hasn't been already.
  if (intNameMap.empty ())
  {
    // set up internal name map
    std::string tmpstr(getName()+"_branch");
    spiceInternalName (tmpstr);
    intNameMap[ li_Bra ] = tmpstr;
  }

  return intNameMap;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/22/02
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs( const std::vector<int> & staLIDVecRef )
{
  AssertLIDs(staLIDVecRef.size() == numStateVars);

  // copy over the global ID lists.
  staLIDVec = staLIDVecRef;

  li_fstate = staLIDVec[0];
}

//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL
// Creation Date : 08/21/02
//-----------------------------------------------------------------------------
const std::vector< std::vector<int> > & Instance::jacobianStamp() const
{
  return jacStamp_BASE;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerJacLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL
// Creation Date : 08/27/02
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );

  APosEquBraVarOffset = jacLIDVec[0][0];
  ANegEquBraVarOffset = jacLIDVec[1][0];
  ABraEquPosNodeOffset = jacLIDVec[2][0];
  ABraEquNegNodeOffset = jacLIDVec[2][1];
  ABraEquBraVarOffset = jacLIDVec[2][2];
}

//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/01/01
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState ()
{
  double * solVec = extData.nextSolVectorRawPtr;
  double * staVec = extData.nextStaVectorRawPtr;
  double current = solVec[li_Bra];
  if( (getSolverState().dcopFlag) && ICGiven )
    current = IC;

  f0 = L*current;
  staVec[li_fstate] = f0;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/01/01
//-----------------------------------------------------------------------------
bool Instance::updateSecondaryState ()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEQVector
//
// Purpose       : Loads the Q-vector contributions for a single
//                 instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/26/03
//-----------------------------------------------------------------------------
bool Instance::loadDAEQVector ()
{
  double * qVec = extData.daeQVectorRawPtr;

  qVec[li_Bra] += f0;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/26/03
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  double current;
  double coef;

  double * solVec = extData.nextSolVectorRawPtr;
  double * fVec = extData.daeFVectorRawPtr;

  double v_pos = solVec[li_Pos];
  double v_neg = solVec[li_Neg];
  double vind = v_pos-v_neg;

  // In the case that an initial condition is specified, the inductor is set up
  // like a current source for the DC operating point. We don't deal with the
  // node voltages in that case, so set coef to 0.
  if (getSolverState().dcopFlag && ICGiven)
  {
    current = IC;
    coef = 0.0;
  }
  else
  {
    current = solVec[li_Bra];
    coef = -vind;
  }

  // load the current into the two KCL rhs vector elements
  fVec[li_Pos] += current;
  fVec[li_Neg] += -current;
  fVec[li_Bra] += coef;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx
//
// Purpose       : Loads the Q-vector contributions for a single
//                 instance.
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  N_LAS_Matrix & dQdx = *(extData.dQdxMatrixPtr);
  dQdx[li_Bra][ABraEquBraVarOffset] += L;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  N_LAS_Matrix & dFdx = *(extData.dFdxMatrixPtr);
  if ( getSolverState().dcopFlag && ICGiven )
  {
    // In the case that an initial condition is specified for an
    // inductor, the DC op should be set up like a current source just
    // for the operating point calculation.
    dFdx[li_Pos][APosEquBraVarOffset]  += 0.0;
    dFdx[li_Neg][ANegEquBraVarOffset]  += 0.0;
    dFdx[li_Bra][ABraEquPosNodeOffset] += 0.0;
    dFdx[li_Bra][ABraEquNegNodeOffset] += 0.0;
    dFdx[li_Bra][ABraEquBraVarOffset]  += 1.0;
  }
  else
  {
    dFdx[li_Pos][APosEquBraVarOffset]  += 1.0;
    dFdx[li_Neg][ANegEquBraVarOffset]  -= 1.0;
    dFdx[li_Bra][ABraEquPosNodeOffset] -= 1.0;
    dFdx[li_Bra][ABraEquNegNodeOffset] += 1.0;
    dFdx[li_Bra][ABraEquBraVarOffset]  += 0.0;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setIC
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 1/11/02
//-----------------------------------------------------------------------------
bool Instance::setIC ()
{
  int i_bra_sol;
  int i_f_state;
  double * nextStaVector = extData.nextStaVectorRawPtr;
  double * currStaVector = extData.currStaVectorRawPtr;

  double * nextSolVector = extData.nextSolVectorRawPtr;
  double * currSolVector = extData.currSolVectorRawPtr;

  if (ICGiven)
  {
    f0 = L*IC;
    currStaVector[li_fstate] = f0;
    nextStaVector[li_fstate] = f0;

    currSolVector[li_Bra] = IC;
    nextSolVector[li_Bra] = IC;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::varTypes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 2/17/04
//-----------------------------------------------------------------------------
void Instance::varTypes( std::vector<char> & varTypeVec )
{
  varTypeVec.resize(1);
  varTypeVec[0] = 'I';
}

//-----------------------------------------------------------------------------
// Function      : Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/03/02
//-----------------------------------------------------------------------------
bool Model::processParams ()
{
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
// Function      : Model::Model
// Purpose       : block constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
Model::Model(
  const Configuration & configuration,
  const ModelBlock &    MB,
  const FactoryBlock &  factory_block)
  : DeviceModel(MB, configuration.getModelParameters(), factory_block),
    L(0.0),
    IC(0.0),
    tempCoeff1(0.0),
    tempCoeff2(0.0),
    tnom(getDeviceOptions().tnom),
    tnomGiven(0)
{

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to .model line and constant defaults from metadata:
  setModParams (MB.params);

  // Set any non-constant parameter defaults:
  if (!given("TNOM"))
    tnom = getDeviceOptions().tnom;

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:

  processParams ();
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
  {
    delete (*iter);
  }

}

//-----------------------------------------------------------------------------
// Function      : Model::printOutInstances
// Purpose       : debugging tool.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/04/00
//-----------------------------------------------------------------------------
std::ostream &Model::printOutInstances(std::ostream &os) const
{
  std::vector<Instance*>::const_iterator iter;
  std::vector<Instance*>::const_iterator first = instanceContainer.begin();
  std::vector<Instance*>::const_iterator last  = instanceContainer.end();

  int i, isize;
  isize = instanceContainer.size();

  os << std::endl;
  os << "Number of Inductor instances: " << isize << std::endl;
  os << "    name=\t\tmodelName\tParameters" << std::endl;
  for (i=0, iter=first; iter!=last; ++iter, ++i)
  {
    os << "  " << i << ": " << (*iter)->getName() << "\t";
    os << getName();
    os << "\t\tL = " << (*iter)->L;
    os << "\tIC = " << (*iter)->IC;
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


// Inductor Master functions:

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
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & inst = *(*it);

    double current = solVec[inst.li_Bra];

    if( (getSolverState().dcopFlag) && inst.ICGiven )
      current = inst.IC;

    inst.f0 = inst.L*current;
    staVec[inst.li_fstate] = inst.f0;
  }

  return true;
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
     Instance & inst = *(*it);

    double current = 0.0;
    double coef = 0.0;

    double v_pos = solVec[inst.li_Pos];
    double v_neg = solVec[inst.li_Neg];
    double vind = v_pos-v_neg;

    // In the case that an initial condition is specified, the inductor is set up
    // like a current source for the DC operating point. We don't deal with the
    // node voltages in that case, so set coef to 0.
    if (getSolverState().dcopFlag && inst.ICGiven)
    {
      current = inst.IC;
      solVec[inst.li_Bra] = current;
      coef = 0.0;
    }
    else
    {
      current = solVec[inst.li_Bra];
      coef = -vind;
    }

    // load the current into the two KCL rhs vector elements
    fVec[inst.li_Pos] += current;
    fVec[inst.li_Neg] += -current;
    fVec[inst.li_Bra] += coef;
    qVec[inst.li_Bra] += inst.f0;
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
     Instance & inst = *(*it);

     if ( getSolverState().dcopFlag && inst.ICGiven )
    {
      // In the case that an initial condition is specified for an
      // inductor, the DC op should be set up like a current source just
      // for the operating point calculation.
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
      *inst.fBraEquBraVarPtr  += 1.0;
#else
      dFdx[inst.li_Bra][inst.ABraEquBraVarOffset]  += 1.0;
#endif
    }
    else
    {
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
      *inst.fPosEquBraVarPtr  += 1.0;
      *inst.fNegEquBraVarPtr  -= 1.0;
      *inst.fBraEquPosNodePtr -= 1.0;
      *inst.fBraEquNegNodePtr += 1.0;
#else
      dFdx[inst.li_Pos][inst.APosEquBraVarOffset]  += 1.0;
      dFdx[inst.li_Neg][inst.ANegEquBraVarOffset]  -= 1.0;
      dFdx[inst.li_Bra][inst.ABraEquPosNodeOffset] -= 1.0;
      dFdx[inst.li_Bra][inst.ABraEquNegNodeOffset] += 1.0;
#endif
    }

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
    *inst.qBraEquBraVarPtr += inst.L;
#else
    dQdx[inst.li_Bra][inst.ABraEquBraVarOffset] += inst.L;
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
    .registerDevice("l", 1)
    .registerModelType("l", 1);
}

} // namespace Inductor
} // namespace Device
} // namespace Xyce

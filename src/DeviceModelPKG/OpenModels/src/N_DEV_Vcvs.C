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
// Filename       : $RCSfile: N_DEV_Vcvs.C,v $
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
// Revision Number: $Revision: 1.131.2.2 $
//
// Revision Date  : $Date: 2014/03/06 23:33:43 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------
#include <Xyce_config.h>

// ---------- Standard Includes ----------
#ifdef HAVE_CSTDIO
#include <cstdio>
#else
#include <stdio.h>
#endif

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_Vcvs.h>
#include <N_DEV_Message.h>
#include <N_ERH_ErrorMgr.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>

namespace Xyce {
namespace Device {


namespace Vcvs {


void Traits::loadInstanceParameters(ParametricData<Vcvs::Instance> &p)
{
// Set up double precision variables:
    p.addPar ("G", 0.0, false, ParameterType::NO_DEP,
        &Vcvs::Instance::Gain,
        NULL,U_NONE,CAT_NONE,"Gain");
}

void Traits::loadModelParameters(ParametricData<Vcvs::Model> &p)
{
}



std::vector< std::vector<int> > Instance::jacStamp;

//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : "instance block" constructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 12/20/00
//-----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock &         IBref,
  Model &                       Viter,
  const FactoryBlock &          factory_block)
  : DeviceInstance(IBref, configuration.getInstanceParameters(), factory_block),
    IB(IBref),
    model_(Viter),
    Gain(1.0),
    li_Pos(-1),
    li_Neg(-1),
    li_Bra(-1),
    li_ContPos(-1),
    li_ContNeg(-1),
    ABraEquPosNodeOffset(-1),
    ABraEquNegNodeOffset(-1),
    ABraEquContPosNodeOffset(-1),
    ABraEquContNegNodeOffset(-1),
    APosEquBraVarOffset(-1),
    ANegEquBraVarOffset(-1),

    f_BraEquPosNodePtr(0),
    f_BraEquNegNodePtr(0),
    f_BraEquContPosNodePtr(0),
    f_BraEquContNegNodePtr(0),
    f_PosEquBraVarPtr(0),
    f_NegEquBraVarPtr(0)
{
  numIntVars   = 1;
  numExtVars   = 4;
  numStateVars = 0;

  std::string prefix(getName() + "-");
  std::string name2;

  name2 = prefix + "I_branch";


  if( jacStamp.empty() )
  {
    jacStamp.resize(5);
    jacStamp[0].resize(1);
    jacStamp[0][0]=4;
    jacStamp[1].resize(1);
    jacStamp[1][0]=4;
    jacStamp[4].resize(4);
    jacStamp[4][0]=0;
    jacStamp[4][1]=1;
    jacStamp[4][2]=2;
    jacStamp[4][3]=3;
  }


  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  // Set any non-constant parameter defaults:
  if (!given("G"))
  {
    UserError0(*this) << "Could not find Gain parameter in instance.";
  }
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
// Creator       : Robert Hoekstra, SNL, Computational Science
// Creation Date : 6/21/02
//-----------------------------------------------------------------------------
void Instance::registerLIDs ( const std::vector<int> & intLIDVecRef,
                                        const std::vector<int> & extLIDVecRef )
{
  std::string msg;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    Xyce::dout() << std::endl << section_divider << std::endl;
    Xyce::dout() << "  VcvsInstance::registerLIDs" << std::endl;
    Xyce::dout() << "  name = " << getName() << std::endl;
  }
#endif

  // Check if the size of the ID lists corresponds to the
  // proper number of internal and external variables.
  int numInt = intLIDVecRef.size();
  int numExt = extLIDVecRef.size();

  if (numInt != numIntVars)
  {
    msg = "Instance::registerLIDs:";
    msg += "numInt != numIntVars";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
  }

  if (numExt != numExtVars)
  {
    msg = "Instance::registerLIDs:";
    msg += "numExt != numExtVars";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
  }

  // copy over the global ID lists.
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  // Now use these lists to obtain the indices into the
  // linear algebra entities.  This assumes an order.
  // For the matrix  indices, first do the rows.

  li_Pos = extLIDVec[0];
  li_Neg = extLIDVec[1];
  li_ContPos = extLIDVec[2];
  li_ContNeg = extLIDVec[3];

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    Xyce::dout() << "  li_Pos = " << li_Pos << std::endl;
    Xyce::dout() << "  li_Neg = " << li_Neg << std::endl;
    Xyce::dout() << "  li_ContPos = " << li_ContPos << std::endl;
    Xyce::dout() << "  li_ContNeg = " << li_ContNeg << std::endl;
  }
#endif

  li_Bra = intLIDVec[0];

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
    Xyce::dout() << "  li_Bra = " << li_Bra << std::endl;
#endif

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
    Xyce::dout() << section_divider << std::endl;
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
    // set up the internal names map:
    std::string tmpstr(getName()+"_branch");
    spiceInternalName (tmpstr);
    intNameMap[ li_Bra ] = tmpstr;
  }

  return intNameMap;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateGIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/21/02
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs( const std::vector<int> & staLIDVecRef )
{
  AssertLIDs(staLIDVecRef.size() == numStateVars);
}

//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 9/2/02
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
// Creation Date : 9/2/02
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );

  APosEquBraVarOffset = jacLIDVec[0][0];
  ANegEquBraVarOffset = jacLIDVec[1][0];
  ABraEquPosNodeOffset = jacLIDVec[4][0];
  ABraEquNegNodeOffset = jacLIDVec[4][1];
  ABraEquContPosNodeOffset = jacLIDVec[4][2];
  ABraEquContNegNodeOffset = jacLIDVec[4][3];
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupPointers
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/12/08
//-----------------------------------------------------------------------------
void Instance::setupPointers ()
{
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  N_LAS_Matrix & dFdx = *(extData.dFdxMatrixPtr);

  f_PosEquBraVarPtr = &(dFdx[li_Pos][ APosEquBraVarOffset ]);
  f_NegEquBraVarPtr = &(dFdx[li_Neg][ ANegEquBraVarOffset ]);
  f_BraEquPosNodePtr = &(dFdx[li_Bra][ ABraEquPosNodeOffset ]);
  f_BraEquNegNodePtr = &(dFdx[li_Bra][ ABraEquNegNodeOffset ]);
  f_BraEquContPosNodePtr = &(dFdx[li_Bra][ ABraEquContPosNodeOffset ]);
  f_BraEquContNegNodePtr = &(dFdx[li_Bra][ ABraEquContNegNodeOffset ]);
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
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 Vcvs instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/19/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  double * fVec = extData.daeFVectorRawPtr;
  double * solVec = extData.nextSolVectorRawPtr;

  double v_pos = solVec[li_Pos];
  double v_neg = solVec[li_Neg];
  double v_cont_pos = solVec[li_ContPos];
  double v_cont_neg = solVec[li_ContNeg];

  double i_bra = solVec[li_Bra];

  fVec[li_Pos] += i_bra;
  fVec[li_Neg] += -i_bra;

  double src = Gain * ( v_cont_pos - v_cont_neg ) - v_pos + v_neg;
  fVec[li_Bra] += -src;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 resistor  instance.
//
// Special Notes : The F-vector is an algebraic constaint.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/19/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  N_LAS_Matrix & dFdx = *(extData.dFdxMatrixPtr);

  dFdx[li_Pos][APosEquBraVarOffset] += 1.0;
  dFdx[li_Neg][ANegEquBraVarOffset] -= 1.0;

  dFdx[li_Bra][ABraEquPosNodeOffset] += 1.0;
  dFdx[li_Bra][ABraEquNegNodeOffset] -= 1.0;
  dFdx[li_Bra][ABraEquContPosNodeOffset] -= Gain;
  dFdx[li_Bra][ABraEquContNegNodeOffset] += Gain;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::varTypes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/17/04
//-----------------------------------------------------------------------------
void Instance::varTypes( std::vector<char> & varTypeVec )
{
  varTypeVec.resize(1);
  varTypeVec[0] = 'I';
}

// Class Model

//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : "Model block" constructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 12/21/00
//-----------------------------------------------------------------------------
Model::Model(
  const Configuration & configuration,
  const ModelBlock &    MB,
  const FactoryBlock &  factory_block)
  : DeviceModel(MB, configuration.getModelParameters(), factory_block)
{
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


// Vcvs Master functions:

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
#ifdef _OMP
#pragma omp parallel for
#endif
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & vi = *(*it);

    double v_pos = solVec[vi.li_Pos];
    double v_neg = solVec[vi.li_Neg];
    double v_cont_pos = solVec[vi.li_ContPos];
    double v_cont_neg = solVec[vi.li_ContNeg];

    double i_bra = solVec[vi.li_Bra];

    fVec[vi.li_Pos] += i_bra;
    fVec[vi.li_Neg] += -i_bra;

    double src = vi.Gain * ( v_cont_pos - v_cont_neg ) - v_pos + v_neg;
    fVec[vi.li_Bra] += -src;

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
    Instance & vi = *(*it);

#ifndef Xyce_NONPOINTER_MATRIX_LOAD

    *vi.f_PosEquBraVarPtr += 1.0;
    *vi.f_NegEquBraVarPtr -= 1.0;

    *vi.f_BraEquPosNodePtr += 1.0;
    *vi.f_BraEquNegNodePtr -= 1.0;
    *vi.f_BraEquContPosNodePtr -= vi.Gain;
    *vi.f_BraEquContNegNodePtr += vi.Gain;
#else

    dFdx[vi.li_Pos][vi.APosEquBraVarOffset] += 1.0;
    dFdx[vi.li_Neg][vi.ANegEquBraVarOffset] -= 1.0;

    dFdx[vi.li_Bra][vi.ABraEquPosNodeOffset] += 1.0;
    dFdx[vi.li_Bra][vi.ABraEquNegNodeOffset] -= 1.0;
    dFdx[vi.li_Bra][vi.ABraEquContPosNodeOffset] -= vi.Gain;
    dFdx[vi.li_Bra][vi.ABraEquContNegNodeOffset] += vi.Gain;
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
    .registerDevice("e", 1)
    .registerModelType("e", 1);
}

} // namespace Vcvs
} // namespace Device
} // namespace Xyce

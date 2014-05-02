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
// Filename       : $RCSfile: N_DEV_Resistor3.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
//
// Creation Date : 2/1/10
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.40.2.3 $
//
// Revision Date  : $Date: 2014/03/06 23:33:43 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------
#include <Xyce_config.h>

#include <memory>

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_Resistor3.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_Message.h>
#include <N_ERH_ErrorMgr.h>

#include <N_DEV_Resistor.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>
#include <N_UTL_BreakPoint.h>

namespace Xyce {
namespace Device {
namespace Resistor3 {

void Traits::loadInstanceParameters(ParametricData<Resistor3::Instance> &p)
{
  p.addPar ("R", 1000.0, false,  ParameterType::TIME_DEP,
            &Resistor3::Instance::R,    NULL,
            U_OHM, CAT_NONE, "Resistance");

  p.addPar ("L", 0.0, false,     ParameterType::NO_DEP,
            &Resistor3::Instance::length, NULL,
            U_METER, CAT_NONE, "Length");

  p.addPar ("W", 0.0, false,     ParameterType::NO_DEP,
            &Resistor3::Instance::width,  NULL,
            U_METER, CAT_NONE, "Width");

  p.addPar ("TEMP", 0.0, false,  ParameterType::TIME_DEP,
            &Resistor3::Instance::temp, NULL,
            U_DEGC, CAT_NONE, "Device temperature");

  p.addPar ("TC1",   0.0, false,   ParameterType::NO_DEP,
            &Resistor3::Instance::tempCoeff1,
            &Resistor3::Instance::tempCoeff1Given,
            U_DEGCM1, CAT_NONE, "Linear Temperature Coefficient");

  p.addPar ("TC2",   0.0, false,   ParameterType::NO_DEP,
            &Resistor3::Instance::tempCoeff2,
            &Resistor3::Instance::tempCoeff2Given,
            U_DEGCM2, CAT_NONE, "Quadratic Temperature Coefficient");

  // This call tells the parameter handling code that TC can be specified
  // as a vector with up to two elements as in TC=a,b.  It then translates
  // TC=a,b into TC1=a TC2=b.  Likewise, TC=a will translate into TC1=a
  p.makeVector ("TC", 2);

  p.addPar ("DTEMP",   0.0, false,   ParameterType::NO_DEP,
            &Resistor3::Instance::dtemp,
            &Resistor3::Instance::dtempGiven,
            U_DEGC, CAT_NONE, "Device Temperature -- For compatibility only. Parameter is NOT used");
}

void Traits::loadModelParameters(ParametricData<Resistor3::Model> &p)
{}

std::vector< std::vector<int> > Instance::jacStamp;
std::vector< std::vector<int> > Instance::jacStampPDE;


// Class Instance
//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : "instance block" constructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/1/10
//-----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock &         IB,
  Model &                       Viter,
  const FactoryBlock &          factory_block)
  : DeviceInstance(IB, configuration.getInstanceParameters(), factory_block),
    model_(Viter),
    srcCurrent(0.0),
    srcVoltage(0.0),
    srcDrop(0.0),
    srcBC(0.0),
    
    li_Pos(-1),
    li_Neg(-1),
    li_Bra(-1),
    tempCoeff1(0.0),
  tempCoeff2(0.0),
  dtemp(0.0),
  tempCoeff1Given(false),
  tempCoeff2Given(false),
  dtempGiven(false),

  ABraEquPosNodeOffset(-1),
  ABraEquNegNodeOffset(-1),
  APosEquBraVarOffset(-1),
  ANegEquBraVarOffset(-1),
  ABraEquBraVarOffset(-1),
  APosEquPosNodeOffset(-1),
  ANegEquNegNodeOffset(-1),

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  fBraEquPosNodePtr(0),
  fBraEquNegNodePtr(0),
  fPosEquBraVarPtr(0),
  fNegEquBraVarPtr(0),
  fPosEquPosNodePtr(0),
  fNegEquNegNodePtr(0),
  fBraEquBraVarPtr(0),
#endif

  scale(1.0),
  nlstep(-1),
  source(0.0),
  v_pos(0.0),
  v_neg(0.0),
  i_bra(0.0)
{
  numIntVars   = 1;
  numExtVars   = 2;
  numStateVars = 0;

  if( jacStamp.empty() )
  {
    jacStamp.resize(3);
    jacStamp[0].resize(1);
    jacStamp[0][0] = 2;
    jacStamp[1].resize(1);
    jacStamp[1][0] = 2;
    jacStamp[2].resize(2);
    jacStamp[2][0] = 0;
    jacStamp[2][1] = 1;

    // PDE supporting stamp.  This includes diagonal elements, needed by the
    // 2-level Newton.
    jacStampPDE.resize(3);
    jacStampPDE[0].resize(2);
    jacStampPDE[0][0] = 0;
    jacStampPDE[0][1] = 2;
    jacStampPDE[1].resize(2);
    jacStampPDE[1][0] = 1;
    jacStampPDE[1][1] = 2;
    jacStampPDE[2].resize(3);
    jacStampPDE[2][0] = 0;
    jacStampPDE[2][1] = 1;
    jacStampPDE[2][2] = 2;
  }

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  processParams();

  // Calculate any parameters specified as expressions:
  updateDependentParameters();
  processParams();
}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/1/10
//-----------------------------------------------------------------------------
Instance::~Instance ()
{
}

//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/1/10
//-----------------------------------------------------------------------------
bool Instance::processParams ()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/1/10
//-----------------------------------------------------------------------------
void Instance::registerLIDs ( const std::vector<int> & intLIDVecRef,
	                                const std::vector<int> & extLIDVecRef)
{
  std::string msg;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 )
  {
    Xyce::dout() << std::endl << section_divider << std::endl;
    Xyce::dout() << "  Instance::registerLIDs" << std::endl;
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
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL_0,msg);
  }

  if (numExt != numExtVars)
  {
    msg = "Instance::registerLIDs:";
    msg += "numExt != numExtVars";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL_0,msg);
  }

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
  if (getDeviceOptions().debugLevel > 0 )
  {
    Xyce::dout() << "  li_Pos = " << li_Pos << std::endl;
    Xyce::dout() << "  li_Neg = " << li_Neg << std::endl;
    Xyce::dout() << "  li_Bra = " << li_Bra << std::endl;
    Xyce::dout() << section_divider << std::endl;
  }
#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::getIntNameMap
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/1/10
//-----------------------------------------------------------------------------
std::map<int,std::string> & Instance::getIntNameMap ()
{
  // set up the internal name map, if it hasn't been already.
  if (intNameMap.empty ())
  {
    // set up the internal name map
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
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/1/10
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
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/1/10
//-----------------------------------------------------------------------------
const std::vector< std::vector<int> > & Instance::jacobianStamp() const
{
  if (getSolverState().PDESystemFlag)
    return jacStampPDE;

  return jacStamp;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerJacLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/1/10
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );
  if (getSolverState().PDESystemFlag)
  {
    APosEquBraVarOffset  = jacLIDVec[0][1];
    ANegEquBraVarOffset  = jacLIDVec[1][1];
    ABraEquPosNodeOffset = jacLIDVec[2][0];
    ABraEquNegNodeOffset = jacLIDVec[2][1];
  }
  else
  {
    APosEquBraVarOffset  = jacLIDVec[0][0];
    ANegEquBraVarOffset  = jacLIDVec[1][0];
    ABraEquPosNodeOffset = jacLIDVec[2][0];
    ABraEquNegNodeOffset = jacLIDVec[2][1];
  }

}

//-----------------------------------------------------------------------------
// Function      : Instance::setupPointers
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/1/10
//-----------------------------------------------------------------------------
void Instance::setupPointers ()
{
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  N_LAS_Matrix & dFdx = *(extData.dFdxMatrixPtr);

  fPosEquBraVarPtr = &(dFdx[li_Pos][APosEquBraVarOffset]);
  fNegEquBraVarPtr = &(dFdx[li_Neg][ANegEquBraVarOffset]);
  fBraEquPosNodePtr = &(dFdx[li_Bra][ABraEquPosNodeOffset]);
  fBraEquNegNodePtr = &(dFdx[li_Bra][ABraEquNegNodeOffset]);
#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/1/10
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars ()
{
  double * solVec = extData.nextSolVectorRawPtr;

  // Get the value for the source.
  source = 0.0;

  // get the value for v_pos, v_neg, i_bra
  v_pos = solVec[li_Pos];
  v_neg = solVec[li_Neg];
  i_bra = solVec[li_Bra];

  srcCurrent = i_bra;
  srcDrop    = (v_pos-v_neg);
  srcBC      = source;
  srcVoltage = srcDrop-srcBC;

  if( getDeviceOptions().scale_src != 0.0 )
  {
    srcCurrent *= scale;

    // first newton step, generate new scaling
    //if( nlsMgrPtr->getNonLinearIter() != nlstep )
    if( getSolverState().newtonIter != nlstep )
    {
      //nlstep = nlsMgrPtr->getNonLinearIter();
      nlstep = getSolverState().newtonIter;
      double new_scale = fabs(i_bra) * scale;

      scale = Xycemax( new_scale, getDeviceOptions().scale_src );
    }

    srcVoltage *= scale;
    srcDrop    *= scale;
    srcBC      *= scale;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/1/10
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState ()
{
  bool bsuccess = updateIntermediateVars ();
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 vsrc instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/1/10
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  // Using values determined during the loadRHS function call.
  // (srcCurrent, srcVoltage).

  double * fVec = extData.daeFVectorRawPtr;

  fVec[li_Pos] += srcCurrent;
  fVec[li_Neg] += -srcCurrent;
  fVec[li_Bra] += srcDrop;
  fVec[li_Bra] -= srcBC;

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
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/1/10
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  N_LAS_Matrix & dFdx = *(extData.dFdxMatrixPtr);

  dFdx[li_Pos][APosEquBraVarOffset] += 1.0;
  dFdx[li_Neg][ANegEquBraVarOffset] -= 1.0;
  dFdx[li_Bra][ABraEquPosNodeOffset] += 1.0;
  dFdx[li_Bra][ABraEquNegNodeOffset] -= 1.0;

  return true;
}

// end of new-DAE functions

//-----------------------------------------------------------------------------
// Function      : Instance::getMaxTimeStepSize
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/1/10
//-----------------------------------------------------------------------------
double Instance::getMaxTimeStepSize  ()
{
  double maxStep = 1.0e99;
  return maxStep;
}

//-----------------------------------------------------------------------------
// Function      : Instance::varTypes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/1/10
//-----------------------------------------------------------------------------
void Instance::varTypes( std::vector<char> & varTypeVec )
{
  varTypeVec.resize(1);
  varTypeVec[0] = 'I';
}

// Class Model

//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : copy constructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/1/10
//-----------------------------------------------------------------------------
Model::Model(
  const Configuration & configuration,
  const ModelBlock &    MB,
  const FactoryBlock &  factory_block)
  : DeviceModel(MB, configuration.getModelParameters(), factory_block),
    DC_TRAN (0)
{
}

//-----------------------------------------------------------------------------
// Function      : Model::~Model
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/1/10
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

// Additional Declarations

//-----------------------------------------------------------------------------
// Function      : Model::printOutInstances
// Purpose       : debugging tool.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/1/10
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



//-----------------------------------------------------------------------------
// Vsrc Master functions:
//-----------------------------------------------------------------------------


//-----------------------------------------------------------------------------
// Function      : Master::updateState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/1/10
//-----------------------------------------------------------------------------
bool Master::updateState (double * solVec, double * staVec, double * stoVec)
{

  // // first loop over the models:
  // for (ModelMap::const_iterator model_it = getModelMap().begin(); model_it != getModelMap().end(); ++model_it)
  // {
  //   // loop over the instances for this model.
  //   InstanceVector::const_iterator first = (*model_it).second->instanceContainer.begin();
  //   InstanceVector::const_iterator last = (*model_it).second->instanceContainer.end();

  //   for (InstanceVector::const_iterator it = first; it != last; ++it)
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
    {
      Instance & vi = *(*it);

      // Get the value for the source.
      vi.source = 0.0;

      // get the value for v_pos, v_neg, i_bra
      vi.v_pos = solVec[vi.li_Pos];
      vi.v_neg = solVec[vi.li_Neg];
      vi.i_bra  = solVec[vi.li_Bra];

      vi.srcCurrent = vi.i_bra;
      vi.srcDrop    = (vi.v_pos-vi.v_neg);
      vi.srcBC      = vi.source;
      vi.srcVoltage = vi.srcDrop-vi.srcBC;

      if( getDeviceOptions().scale_src != 0.0 )
      {
        vi.srcCurrent *= vi.scale;

        // first newton step, generate new scaling
        if( getSolverState().newtonIter != vi.nlstep )
        {
          vi.nlstep = getSolverState().newtonIter;
          double new_scale = fabs(vi.i_bra) * vi.scale;

          vi.scale = Xycemax( new_scale, getDeviceOptions().scale_src );
        }

        vi.srcVoltage *= vi.scale;
        vi.srcDrop    *= vi.scale;
        vi.srcBC      *= vi.scale;
      }
    }
//  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Master::loadDAEVectors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/1/10
//-----------------------------------------------------------------------------
bool Master::loadDAEVectors (double * solVec, double * fVec, double *qVec,  double * storeLeadF, double * storeLeadQ)
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & vi = *(*it);


    fVec[vi.li_Pos] += vi.srcCurrent;

    fVec[vi.li_Neg] += -vi.srcCurrent;

    fVec[vi.li_Bra] += vi.srcDrop;

    fVec[vi.li_Bra] -= vi.srcBC;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Master::loadDAEMatrices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/1/10
//-----------------------------------------------------------------------------
bool Master::loadDAEMatrices (N_LAS_Matrix & dFdx, N_LAS_Matrix & dQdx)
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & vi = *(*it);

#ifndef Xyce_NONPOINTER_MATRIX_LOAD

    *(vi.fPosEquBraVarPtr) += 1.0;

    *(vi.fNegEquBraVarPtr) -= 1.0;

    *(vi.fBraEquPosNodePtr) += 1.0;

    *(vi.fBraEquNegNodePtr) -= 1.0;
#else

    dFdx[vi.li_Pos][vi.APosEquBraVarOffset] += 1.0;

    dFdx[vi.li_Neg][vi.ANegEquBraVarOffset] -= 1.0;

    dFdx[vi.li_Bra][vi.ABraEquPosNodeOffset] += 1.0;

    dFdx[vi.li_Bra][vi.ABraEquNegNodeOffset] -= 1.0;
#endif

  }

  return true;
}

Device *
Traits::factory(const Configuration &configuration, const FactoryBlock &factory_block)
{
  return new Master(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);
}

void registerDevice() 
{
  Config<Traits>::addConfiguration()
    .registerDevice("r", 3)
    .registerModelType("r", 3);
}

} // namespace Resistor3
} // namespace Device
} // namespace Xyce

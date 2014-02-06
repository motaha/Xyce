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
// Revision Number: $Revision: 1.22.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:38 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------
#include <Xyce_config.h>

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include <N_DEV_Resistor3.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_MatrixLoadData.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>
#include <N_UTL_BreakPoint.h>

namespace Xyce {
namespace Device {

template<>
ParametricData<Resistor3::Instance>::ParametricData()
{
    setNumNodes(2);
    setNumOptionalNodes(0);
    setNumFillNodes(0);
    setModelRequired(0);
    setPrimaryParameter("R");
    addModelType("R");

    // Set up double precision variables:
    addPar ("R", 1000.0, false,  ParameterType::TIME_DEP,
        &Resistor3::Instance::R,    NULL,
        U_OHM, CAT_NONE, "Resistance");

    addPar ("L", 0.0, false,     ParameterType::NO_DEP,
        &Resistor3::Instance::length, NULL,
        U_METER, CAT_NONE, "Length");

    addPar ("W", 0.0, false,     ParameterType::NO_DEP,
        &Resistor3::Instance::width,  NULL,
        U_METER, CAT_NONE, "Width");

    addPar ("TEMP", 0.0, false,  ParameterType::TIME_DEP,
        &Resistor3::Instance::temp, NULL,
        U_DEGC, CAT_NONE, "Temperature");

    addPar ("TC1",   0.0, false,   ParameterType::NO_DEP,
      &Resistor3::Instance::tempCoeff1,
      &Resistor3::Instance::tempCoeff1Given,
      U_DEGCM1, CAT_NONE, "Linear Temperature Coefficient");

    addPar ("TC2",   0.0, false,   ParameterType::NO_DEP,
      &Resistor3::Instance::tempCoeff2,
      &Resistor3::Instance::tempCoeff2Given,
      U_DEGCM2, CAT_NONE, "Quadratic Temperature Coefficient");

    // This call tells the parameter handling code that TC can be specified
    // as a vector with up to two elements as in TC=a,b.  It then translates
    // TC=a,b into TC1=a TC2=b.  Likewise, TC=a will translate into TC1=a
    makeVector ("TC", 2);

    addPar ("DTEMP",   0.0, false,   ParameterType::NO_DEP,
      &Resistor3::Instance::dtemp,
      &Resistor3::Instance::dtempGiven,
      U_DEGC, CAT_NONE, "Device Temperature -- For compatibility only. Parameter is NOT used");
}

template<>
ParametricData<Resistor3::Model>::ParametricData()
{
}

namespace Resistor3 {

vector< vector<int> > Instance::jacStamp;
vector< vector<int> > Instance::jacStampPDE;



ParametricData<Instance> &Instance::getParametricData() {
  static ParametricData<Instance> parMap;

  return parMap;
}

ParametricData<Model> &Model::getParametricData() {
  static ParametricData<Model> parMap;

  return parMap;
}


// Class Instance
//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : "instance block" constructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/1/10
//-----------------------------------------------------------------------------
Instance::Instance(InstanceBlock & IB,
			                     Model & Viter,
                                       MatrixLoadData & mlData1,
                                       SolverState &ss1,
                                       ExternData  &ed1,
                                       DeviceOptions & do1)

  : DeviceInstance(IB,mlData1,ss1,ed1,do1),
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

  // calculate dependent (ie computed) params and check for errors:

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
bool Instance::processParams (string param)
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
void Instance::registerLIDs ( const vector<int> & intLIDVecRef,
	                                const vector<int> & extLIDVecRef)
{
  string msg;

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline(
      "-----------------------------------------------------------------------------");
  if (getDeviceOptions().debugLevel > 0 )
  {
    cout << endl << dashedline << endl;
    cout << "  Instance::registerLIDs" << endl;
    cout << "  name = " << getName() << endl;
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
    cout << "  li_Pos = " << li_Pos << endl;
    cout << "  li_Neg = " << li_Neg << endl;
    cout << "  li_Bra = " << li_Bra << endl;
    cout << dashedline << endl;
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
map<int,string> & Instance::getIntNameMap ()
{
  // set up the internal name map, if it hasn't been already.
  if (intNameMap.empty ())
  {
    // set up the internal name map
    string tmpstr(getName()+"_branch");
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
void Instance::registerStateLIDs( const vector<int> & staLIDVecRef )
{
  string msg;

  // Check if the size of the ID lists corresponds to the
  // proper number of internal and external variables.
  int numSta = staLIDVecRef.size();

  if (numSta != numStateVars)
  {
    msg = "Instance::registerStateLIDs:";
    msg += "numSTa != numStateVars";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL_0,msg);
  }

}

//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystems Modeling
// Creation Date : 2/1/10
//-----------------------------------------------------------------------------
const vector< vector<int> > & Instance::jacobianStamp() const
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
void Instance::registerJacLIDs( const vector< vector<int> > & jacLIDVec )
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
void Instance::varTypes( vector<char> & varTypeVec )
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
Model::Model (const ModelBlock & MB,
                                        SolverState & ss1,
                                        DeviceOptions & do1)
  : DeviceModel(MB,ss1,do1),
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
  vector<Instance*>::iterator iter;
  vector<Instance*>::iterator first = instanceContainer.begin();
  vector<Instance*>::iterator last  = instanceContainer.end();

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
  vector<Instance*>::const_iterator iter;
  vector<Instance*>::const_iterator first = instanceContainer.begin();
  vector<Instance*>::const_iterator last  = instanceContainer.end();

  int i;
  os << endl;
  os << "    name     getModelName()  Parameters" << endl;
  for (i=0, iter=first; iter!=last; ++iter, ++i)
  {
    os << "  " << i << ": " << (*iter)->getName() << "      ";
    os << (*iter)->getModelName();
    os << endl;
  }

  os << endl;

  return os;
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

  // first loop over the models:
  for (ModelMap::const_iterator model_it = getModelMap().begin(); model_it != getModelMap().end(); ++model_it)
  {
    // loop over the instances for this model.
    InstanceVector::const_iterator first = (*model_it).second->instanceContainer.begin();
    InstanceVector::const_iterator last = (*model_it).second->instanceContainer.end();

    for (InstanceVector::const_iterator it = first; it != last; ++it)
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
  }

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
  for (InstanceVector::const_iterator it = getInstanceVector().begin(); it != getInstanceVector().end(); ++it)
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
  for (InstanceVector::const_iterator it = getInstanceVector().begin(); it != getInstanceVector().end(); ++it)
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

} // namespace Resistor3
} // namespace Device
} // namespace Xyce

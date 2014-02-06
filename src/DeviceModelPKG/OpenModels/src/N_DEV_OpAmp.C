//--------------------------------------------------------------------------
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
// Filename       : N_DEV_OpAmp.C
//
// Purpose        : An Ideal OpAmp Model
//
// Special Notes  : This model assumes infinite gain
//                   and unbounded output voltages
//
// Creator        : Brian Fett, SNL
//
// Creation Date  : 07/27/05
//
//-------------------------------------------------------------------------
#include <Xyce_config.h>

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include <N_DEV_OpAmp.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_MatrixLoadData.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>

namespace Xyce {
namespace Device {

template<>
ParametricData<OpAmp::Instance>::ParametricData()
{
    setNumNodes(3);
    setNumOptionalNodes(0);
    setNumFillNodes(0);
    setModelRequired(1);
    addModelType("OPAMP");

    // Set up double precision variables:
    addPar ("FAKEPARAM", 0.0, false, ParameterType::NO_DEP,
            &OpAmp::Instance::FAKEPARAM,
        NULL, U_NONE,CAT_NONE,"");
}

template<>
ParametricData<OpAmp::Model>::ParametricData()
{
#if 0
    // Set up double precision variables:
    addPar ("FAKEPARAM",0.0,false,Pars::NO_DEP,
        &OpAmp::Model::FAKEPARAM,
        NULL, U_NONE,CAT_NONE,"");
#endif
}

namespace OpAmp {

vector< vector<int> > Instance::jacStamp;



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
// Purpose       : instance block constructor
// Special Notes :
// Scope         : public
// Creator       : Brian Fett, SNL
// Creation Date : 08/05/07
//-----------------------------------------------------------------------------
Instance::Instance(InstanceBlock & IB,
                                    Model & iter,
                                    MatrixLoadData & mlData1,
                                    SolverState &ss1,
                                    ExternData  &ed1,
                                    DeviceOptions & do1)

  : DeviceInstance(IB,mlData1,ss1,ed1,do1),
    model_(iter),
  outCurrent(0.0),
  deltaVoltage(0.0),

  FAKEPARAM(0.0),

  li_Pos(-1),
  li_Neg(-1),
  li_Out(-1),
  li_Bra(-1),

  ABraEquPosNodeOffset(-1),
  ABraEquNegNodeOffset(-1),
  AOutEquBraVarOffset(-1),

  v_pos(0.0),
  v_neg(0.0),
  v_out(0.0),
  i_bra(0.0)
{
  numIntVars   = 1;
  numExtVars   = 3;
  numStateVars = 0;

  if( jacStamp.empty() )
  {
    jacStamp.resize(4);    //    V1 V2 V3 I
    jacStamp[2].resize(1); // N1
    jacStamp[2][0] = 3;    // N2
    jacStamp[3].resize(2); // N3          X
    jacStamp[3][0] = 0;    // Br X  X
    jacStamp[3][1] = 1;
  }


  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  // calculate dependent (ie computed) params and check for errors:

}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Brian Fett
// Creation Date : 07/28/05
//-----------------------------------------------------------------------------
Instance::~Instance ()
{
}
//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Brian Fett, SNL
// Creation Date : 08/05/05
//-----------------------------------------------------------------------------
bool Instance::processParams (string param)
{

  return true;
}


// Additional Declarations

//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Brian Fett, SNL
// Creation Date : 08/05/05
//-----------------------------------------------------------------------------
void Instance::registerLIDs ( const vector<int> & intLIDVecRef,
	                                const vector<int> & extLIDVecRef)
{

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline =
      "-----------------------------------------------------------------------------";
  if (getDeviceOptions().debugLevel > 0 )
  {
    cout << endl << dashedline << endl;
    cout << "  OpAmpInstance::registerLIDs" << endl;
    cout << "  name = " << getName() << endl;
  }
#endif

  // Check if the size of the ID lists corresponds to the
  // proper number of internal and external variables.
  int numInt = intLIDVecRef.size();
  int numExt = extLIDVecRef.size();

  if (numInt != numIntVars)
  {
    string msg = "Instance::registerLIDs:";
    msg += "numInt != numIntVars";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL_0,msg);
  }

  if (numExt != numExtVars)
  {
    string msg = "Instance::registerLIDs:";
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
  li_Out = extLIDVec[2];
  li_Bra = intLIDVec[0];

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 )
  {
    cout << "  li_Pos = " << li_Pos << endl;
    cout << "  li_Neg = " << li_Neg << endl;
    cout << "  li_Out = " << li_Out << endl;
    cout << "  li_Bra = " << li_Bra << endl;
    cout << dashedline << endl;
  }
#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/20/02
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs( const vector<int> & staLIDVecRef )
{
  // Check if the size of the ID lists corresponds to the
  // proper number of internal and external variables.
  int numSta = staLIDVecRef.size();

  if (numSta != numStateVars)
  {
    string msg = "Instance::registerStateLIDs:";
    msg += "numSTa != numStateVars";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL_0,msg);
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/21/02
//-----------------------------------------------------------------------------
const vector< vector<int> > & Instance::jacobianStamp() const
{
  return jacStamp;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerJacLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Brian Fett, SNL
// Creation Date : 08/05/05
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const vector< vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );

  AOutEquBraVarOffset  = jacLIDVec[2][0];
  ABraEquPosNodeOffset = jacLIDVec[3][0];
  ABraEquNegNodeOffset = jacLIDVec[3][1];
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Brian Fett, SNL
// Creation Date : 08/05/05
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars ()
{
  N_LAS_Vector * solVectorPtr = extData.nextSolVectorPtr;

  // get the value for v_pos, v_neg, v_out, i_bra

  v_pos = (*solVectorPtr)[li_Pos];
  v_neg = (*solVectorPtr)[li_Neg];
  v_out = (*solVectorPtr)[li_Out];
  i_bra = (*solVectorPtr)[li_Bra];

  outCurrent = i_bra;
  deltaVoltage = v_pos - v_neg;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    const string dashedline2 = "---------------------";
    cout << dashedline2 << endl;
    cout << "  Instance::updateIntermediateVars" << endl;
    cout << "  name = " << getName() <<endl;
    cout << "  v_pos  = " << v_pos << endl;
    cout << "  v_neg  = " << v_neg << endl;
    cout << "  v_out  = " << v_out << endl;
    cout << "  i_bra  = " << i_bra << endl;
    cout << endl;
    cout << dashedline2 << endl;
  }
#endif

  return true;
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
#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << "  Instance::updatePrimaryState" << endl;
  }
#endif
  bool bsuccess = true;
  bsuccess = updateIntermediateVars ();
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 OpAmp instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Brian Fett, SNL
// Creation Date : 08/05/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  // bool bsuccess = true;
  // int i_pos_rhs, i_neg_rhs, i_bra_rhs;

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline2 = "---------------------";
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << dashedline2 << endl;
    cout << "  Instance::loadDAEFVector" << endl;
    cout << "  name       = " << getName() <<endl;
    cout << "  Output Current = " << outCurrent << endl;
    cout << "  Delta  Voltage = " << deltaVoltage << endl;
  }
#endif

  // Using values determined during the loadRHS function call.
  // (outCurrent, deltaVoltage).

  //(*extData.daeFVectorPtr)[li_Pos] += 0;
  //(*extData.daeFVectorPtr)[li_Neg] += 0;

  (*extData.daeFVectorPtr)[li_Out] += outCurrent;

  (*extData.daeFVectorPtr)[li_Bra] += deltaVoltage;

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
// Creator       : Brian Fett, SNL
// Creation Date : 08/02/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  // bool bsuccess = true;

  N_LAS_Matrix * dFdxMatPtr = extData.dFdxMatrixPtr;

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline2 = "---------------------";

  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << dashedline2 << endl;
    cout << "name = " << getName() << endl;
    cout << "\nOPAMP dFdx LOADS\n";
    cout << "Pos,Bra: " << li_Out << ",";
    cout << AOutEquBraVarOffset << ": " << 1.0 << endl;
    cout << "Bra,Pos: " << li_Bra << ",";
    cout << ABraEquPosNodeOffset << ": " << 1.0 << endl;
    cout << "Bra,Neg: " << li_Bra << ",";
    cout << ABraEquNegNodeOffset << ": " << -1.0 << endl;
    cout << "DONE OPAMP dFdx LOAD\n";
  }
#endif

  (*dFdxMatPtr)[li_Out][AOutEquBraVarOffset] += 1.0;
  (*dFdxMatPtr)[li_Bra][ABraEquPosNodeOffset] += 1.0;
  (*dFdxMatPtr)[li_Bra][ABraEquNegNodeOffset] -= 1.0;

  return true;
}

// end of new-DAE functions

// Class Model

//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : model block constructor
// Special Notes :
// Scope         : public
// Creator       : Brian Fett, SNL
// Creation Date : 08/05/05
//-----------------------------------------------------------------------------
Model::Model (const ModelBlock & MB,
                                        SolverState & ss1,
                                        DeviceOptions & do1)
  : DeviceModel(MB,ss1,do1),
  FAKEPARAM(0.0)
{
}

//-----------------------------------------------------------------------------
// Function      : Model::~Model
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/06/00
//-----------------------------------------------------------------------------
Model::~Model ()
{
  vector<Instance*>::iterator iter = instanceContainer.begin();
  vector<Instance*>::iterator last  = instanceContainer.end();

  for ( ; iter!=last; ++iter)
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
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/03/00
//-----------------------------------------------------------------------------
std::ostream &Model::printOutInstances(std::ostream &os) const
{
  vector<Instance*>::const_iterator iter = instanceContainer.begin();
  vector<Instance*>::const_iterator last  = instanceContainer.end();

  os << endl;
  os << "    name     getModelName()  Parameters" << endl;
  for (int i=0; iter!=last; ++iter, ++i)
  {
    os << "  " << i << ": " << (*iter)->getName() << "      ";
    os << (*iter)->getModelName();
    os << endl;
  }

  os << endl;

  return os;
}

} // namespace OpAmp
} // namespace Device
} // namespace Xyce

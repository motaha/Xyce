//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2011  Sandia Corporation
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
// Filename       : $RCSfile: N_DEV_Synapse.C,v $
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
// Revision Number: $Revision: 1.21.2.1 $
//
// Revision Date  : $Date: 2013/10/03 17:23:33 $
//
// Current Owner  : $Author: tvrusso $
//----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
//
#include <N_DEV_Const.h>
#include <N_DEV_Synapse.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_MatrixLoadData.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>

namespace Xyce {
namespace Device {
template<>
ParametricData<Synapse::Instance>::ParametricData()
{
  // Set up configuration constants:
  setNumNodes(2);
  setNumOptionalNodes(0);
  setNumFillNodes(0);
  setModelRequired(0);
  addModelType("SYNAPSE");
}

template<>
ParametricData<Synapse::Model>::ParametricData()
{
  // Default values are taken from Destexhe 98, values for AMPA receptor
  addPar ("GMAX", 1.0E-9, false, ParameterType::NO_DEP,
          &Synapse::Model::gMax, NULL,
          U_OHMM1, CAT_NONE, "Maximal Synaptic Conductance");

  addPar ("EREV", 0.0, false, ParameterType::NO_DEP,
          &Synapse::Model::eRev, NULL,
          U_VOLT, CAT_NONE, "Reversal Potential");

  addPar ("ALPHA", 1.1E-6, false, ParameterType::NO_DEP,
          &Synapse::Model::alpha, NULL,
          U_SECM1, CAT_NONE, "Forward rate constant for receptor opening");

  addPar ("BETA", 190.0, false, ParameterType::NO_DEP,
          &Synapse::Model::beta, NULL,
          U_SECM1, CAT_NONE, "Backward rate constant for receptor opening");

  addPar ("VP", 0.002, false, ParameterType::NO_DEP,
          &Synapse::Model::vP, NULL,
          U_VOLT, CAT_NONE, "Presynaptic voltage at which neurotransmitter concentration is half-maximal");

  // KP should NOT be 0
  addPar ("KP", 0.005, false, ParameterType::NO_DEP,
          &Synapse::Model::kP, NULL,
          U_VOLT, CAT_NONE, "Steepness parameter for neurotransmitter concentration");

  addPar ("TMAX", 0.001, false, ParameterType::NO_DEP,
          &Synapse::Model::tMax, NULL,
          U_NONE, CAT_NONE, "Maximal neurotransmitter concentration");

  // NOTE:  not including concentration units - TMAX should be in moles/liter, and
  // alpha has concentration in the denominator; they cancel out
}

namespace Synapse {

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
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/03/02
//-----------------------------------------------------------------------------
bool Instance::processParams (string param)
{
  // now set the temperature related stuff.
  //updateTemperature(temp);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : instance block constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
Instance::Instance(
  InstanceBlock & IB,
  Model & Riter,
  MatrixLoadData & mlData1,
  SolverState &ss1,
  ExternData  &ed1,
  DeviceOptions & do1)

  : DeviceInstance(IB, mlData1, ss1, ed1, do1),
    model_(Riter),
    li_Prev(-1),
    li_Post(-1),
    APostEquPostNodeOffset(-1),
    APostEquRNodeOffset(-1),
    AREquPostNodeOffset(-1),
    AREquRNodeOffset(-1),
    f_PostEquPostNodePtr(0),
    f_PostEquRNodePtr(0),
    f_REquPostNodePtr(0),
    f_REquRNodePtr(0)
{
  numIntVars   = 1;
  numExtVars   = 2;
  numStateVars = 0;

  setName(IB.getName());
  setModelName(model_.getName());

  if( jacStamp.empty() )
  {
    // Vpre, Vpost, r
    jacStamp.resize(3);
    jacStamp[0].resize(0);
    jacStamp[1].resize(2);
    jacStamp[2].resize(2);
    jacStamp[1][0] = 1;  // Vpost  PostPostVar
    jacStamp[1][1] = 2;  // r      PostRVar
    jacStamp[2][0] = 0;  // Vpre   RVpre
    jacStamp[2][1] = 2;  // r      RR
  }


  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  // Set any non-constant parameter defaults:

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
// Creation Date : 6/12/02
//-----------------------------------------------------------------------------
void Instance::registerLIDs( const vector<int> & intLIDVecRef,
                             const vector<int> & extLIDVecRef )
{
  string msg;

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline =
    "-------------------------------------------------------------------------"
    "----";
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << endl << dashedline << endl;
    cout << "  SynapseInstance::registerLIDs" << endl;
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

  li_Prev = extLIDVec[0];
  li_Post = extLIDVec[1];

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 )
  {
    cout << "  li_Prev = " << li_Prev << endl;
    cout << "  li_Post = " << li_Post << endl;
  }
#endif

  li_rVar = intLIDVec[0];

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 )
  {
    cout << dashedline << endl;
  }
#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::getIntNameMap
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Christina Warrender, SNL, Cognitive Modeling
// Creation Date : 10/17/11
//-----------------------------------------------------------------------------
map<int,string> & Instance::getIntNameMap ()
{
  // set up the internal name map, if it hasn't been already.
  if (intNameMap.empty ())
  {
    string tmpstr;

    tmpstr = getName() + "_r";
    spiceInternalName (tmpstr);
    intNameMap[ li_rVar ] = tmpstr;
  }
  return intNameMap;
}



//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       : Note that the Synapse does not have any state vars.
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/12/02
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs(const vector<int> & staLIDVecRef )
{
  string msg;

  // Check if the size of the ID lists corresponds to the
  // proper number of internal and external variables.
  int numSta = staLIDVecRef.size();

  if (numSta != numStateVars)
  {
    msg = "Instance::registerLIDs:";
    msg += "numSta != numStateVars";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/20/01
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
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/27/01
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const vector< vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );

  APostEquPostNodeOffset = jacLIDVec[1][0];
  APostEquRNodeOffset = jacLIDVec[1][1];
  AREquPostNodeOffset = jacLIDVec[2][0];
  AREquRNodeOffset = jacLIDVec[2][1];
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
  f_PostEquPostNodePtr = &(dFdx[li_Post][APostEquPostNodeOffset]);
  f_PostEquRNodePtr = &(dFdx[li_Post][APostEquRNodeOffset]);
  f_REquPostNodePtr = &(dFdx[li_rVar][AREquPostNodeOffset]);
  f_REquRNodePtr = &(dFdx[li_rVar][AREquRNodeOffset]);
#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       : update intermediate variables for one diode instance
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, Dept. 9233.
// Creation Date : 3/05/04
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars()
{
  bool bsuccess = true;

  N_LAS_Vector *  solVecPtr = extData.nextSolVectorPtr;

  double vPre  = (*solVecPtr)[li_Prev];
  double vPost  = (*solVecPtr)[li_Post];
  double rVar  = (*solVecPtr)[li_rVar];

  // Load RHS vector element for the positive circuit node KCL equ.
  {
    Sacado::Fad::SFad<double,2> vPostVar( 2, 0, vPost );
    Sacado::Fad::SFad<double,2> rVarS( 2, 1, rVar);

    // parameters
    Sacado::Fad::SFad<double,2> gMaxVar( model_.gMax );
    Sacado::Fad::SFad<double,2> eRevVar( model_.eRev );

    Sacado::Fad::SFad<double,2> resultFad;
    resultFad = PostCurrentEqu( vPostVar, rVarS, gMaxVar, eRevVar );
    ipost = resultFad.val();
    didVpost = resultFad.dx(0);
    didr = resultFad.dx(1);
  }
  {
    Sacado::Fad::SFad<double,2> vPreVar( 2, 0, vPre );
    Sacado::Fad::SFad<double,2> rVarS( 2, 1, rVar);

    // parameters
    Sacado::Fad::SFad<double,2> alphaVar( model_.alpha );
    Sacado::Fad::SFad<double,2> betaVar( model_.beta );
    Sacado::Fad::SFad<double,2> tMaxVar( model_.tMax );
    Sacado::Fad::SFad<double,2> vPVar( model_.vP );
    Sacado::Fad::SFad<double,2> kPVar( model_.kP );

    Sacado::Fad::SFad<double,2> resultFad;
    resultFad = rEquF( vPreVar, rVarS, alphaVar, betaVar, tMaxVar, vPVar, kPVar);
    rFval = resultFad.val();
    drFdVpre = resultFad.dx(0);
    drFdr = resultFad.dx(1);
  }
  return bsuccess;
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
  bool bsuccess = updateIntermediateVars();
  return  bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/29/01
//-----------------------------------------------------------------------------
bool Instance::updateSecondaryState()
{
  return  true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEQVector
//
// Purpose       : Loads the Q-vector contributions for a single
//                 synapse
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 8/31/2010
//-----------------------------------------------------------------------------
bool Instance::loadDAEQVector ()
{
  N_LAS_Vector *  qVecPtr = extData.daeQVectorPtr;
  N_LAS_Vector *  solVecPtr = extData.nextSolVectorPtr;
  double rVar  = (*solVecPtr)[li_rVar];
  (*qVecPtr)[li_rVar] -= rVar;

  return true;
}
//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 Synapse  instance.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/24/03
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  N_LAS_Vector *  fVecPtr = extData.daeFVectorPtr;
  (*fVecPtr)[li_Prev] += 0.0;
  (*fVecPtr)[li_Post] += ipost;
  (*fVecPtr)[li_rVar] += rFval;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx
//
// Purpose       : Loads the dQdx-matrix contributions for a single
//                 synapse instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  N_LAS_Matrix & dQdx = *(extData.dQdxMatrixPtr);
  dQdx[li_rVar][AREquRNodeOffset] += -1.0;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
// Purpose       : Loads the F-vector contributions for a single
//                 Synapse  instance.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  N_LAS_Matrix & dFdx = *(extData.dFdxMatrixPtr);
  dFdx[li_Post][APostEquPostNodeOffset] += didVpost;
  dFdx[li_Post][APostEquRNodeOffset]    += didr;
  dFdx[li_rVar][AREquPostNodeOffset]    += drFdVpre;
  dFdx[li_rVar][AREquRNodeOffset]       += drFdr;

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
bool Instance::updateTemperature ( const double & temp_tmp)
{
  bool bsuccess = true;
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/03/02
//-----------------------------------------------------------------------------
bool Model::processParams (string param)
{
  if (kP == 0.0)
    kP = 0.001;

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
bool Model::processInstanceParams(string param)
{

  vector<Instance*>::iterator iter;
  vector<Instance*>::iterator first = instanceContainer.begin();
  vector<Instance*>::iterator last  = instanceContainer.end();

  for (iter=first; iter!=last; ++iter)
  {
    (*iter)->processParams();
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : model block constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
Model::Model (const ModelBlock & MB,
              SolverState & ss1,
              DeviceOptions & do1)
  : DeviceModel(MB, ss1,do1),
    gMax(0.0),
    eRev(0.0),
    alpha(0.0),
    beta(0.0),
    vP(0.0),
    kP(0.0),
    tMax(0.0)
{

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to .model line and constant defaults from metadata:
  setModParams (MB.params);

  // Set any non-constant parameter defaults:

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:
  processParams();
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
  vector<Instance*>::iterator iter;
  vector<Instance*>::iterator first = instanceContainer.begin();
  vector<Instance*>::iterator last  = instanceContainer.end();

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
  vector<Instance*>::const_iterator iter;
  vector<Instance*>::const_iterator first = instanceContainer.begin();
  vector<Instance*>::const_iterator last  = instanceContainer.end();

  int i,isize;
  isize = instanceContainer.size();
  os << endl;
  os << "Number of Synapse Instances: " << isize << endl;
  os << "    name     getModelName()  Parameters" << endl;
  for (i=0, iter=first; iter!=last; ++iter, ++i)
  {
    os << "  " << i << ": " << (*iter)->getName() << "\t";
    os << (*iter)->getModelName();
    os << endl;
  }

  os << endl;
}

// Synapse Master functions:

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
  for (InstanceVector::const_iterator it = getInstanceVector().begin(); it != getInstanceVector().end(); ++it)
  {
    (*it)->updateIntermediateVars();
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

} // namespace Synapse
} // namespace Device
} // namespace Xyce

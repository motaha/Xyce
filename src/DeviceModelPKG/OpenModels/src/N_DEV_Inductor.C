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
// Revision Number: $Revision: 1.243.2.3 $
//
// Revision Date  : $Date: 2013/10/03 17:23:38 $
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
#include <N_DEV_Inductor.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_MatrixLoadData.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>

namespace Xyce {
namespace Device {

template<>
ParametricData<Inductor::Instance>::ParametricData()
{
  setNumNodes(2);
  setNumOptionalNodes(0);
  setNumFillNodes(0);
  setModelRequired(0);
  setPrimaryParameter("L");
  addModelType("L");

  // Set up double precision variables:
  addPar ("L",    0.0, false, ParameterType::TIME_DEP,
          &Inductor::Instance::baseL,
          NULL,
          U_HENRY, CAT_NONE, "Inductance");

  addPar ("IC",   0.0, false,   ParameterType::NO_DEP,
          &Inductor::Instance::IC,
          &Inductor::Instance::ICGiven,
          U_AMP, CAT_NONE, "Initial current through device");

  addPar ("TEMP", 0.0, false, ParameterType::TIME_DEP,
          &Inductor::Instance::temp,
          &Inductor::Instance::tempGiven,
          U_DEGC, CAT_MATERIAL, "Temperature");

  //Genie 121412. Make Inductor support TC
  addPar ("TC1",   0.0, false,   ParameterType::NO_DEP,
          &Inductor::Instance::tempCoeff1,
          &Inductor::Instance::tempCoeff1Given,
          U_DEGCM1, CAT_NONE, "Linear Temperature Coefficient");

  addPar ("TC2",   0.0, false,   ParameterType::NO_DEP,
          &Inductor::Instance::tempCoeff2,
          &Inductor::Instance::tempCoeff2Given,
          U_DEGCM2, CAT_NONE, "Quadratic Temperature Coefficient");

  // This call tells the parameter handling code that TC can be specified
  // as a vector with up to two elements as in TC=a,b.  It then translates
  // TC=a,b into TC1=a TC2=b.  Likewise, TC=a will translate into TC1=a
  makeVector ("TC", 2);
}

template<>
ParametricData<Inductor::Model>::ParametricData()
{
  // Set up double precision variables:
  addPar ("L", 1.0, false, ParameterType::NO_DEP,
          &Inductor::Model::L,
          NULL,
          U_NONE, CAT_NONE, "Inductance Multiplier");

  addPar ("IC", 0.0, false, ParameterType::NO_DEP,
          &Inductor::Model::IC,
          NULL,
          U_AMP, CAT_NONE, "Initial current through device");

  addPar ("TNOM",  27.0,false, ParameterType::NO_DEP,
          &Inductor::Model::tnom,
          NULL,
          U_DEGC, CAT_MATERIAL, "Reference temperature");

  addPar ("TC1",0.0, false, ParameterType::NO_DEP,
          &Inductor::Model::tempCoeff1,
          NULL,
          U_DEGCM1, CAT_MATERIAL, "First order temperature coeff.");

  addPar ("TC2",0.0, false, ParameterType::NO_DEP,
          &Inductor::Model::tempCoeff2,
          NULL,
          U_DEGCM2, CAT_MATERIAL, "Second order temperature coeff.");
}

namespace Inductor {
//
// static class member inits
//
vector< vector<int> > Instance::jacStamp_BASE;



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
bool Instance::processParams(string param)
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
  //double factor = model_.L*(1.0 + (model_.tempCoeff1)*difference +
  //                       (model_.tempCoeff2)*difference*difference);
  //support specifying TC at the instance line
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
Instance::Instance(InstanceBlock & IB,
                                Model & Iiter,
                                MatrixLoadData & mlData1,
                                SolverState &ss1,
                                ExternData  &ed1,
                                DeviceOptions & do1)

  : DeviceInstance (IB, mlData1, ss1, ed1, do1),
    L(0),
    IC(0),
    ICGiven(false),
    model_(Iiter),
    baseL(0.0),
    temp(getDeviceOptions().temp.dVal()),
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
  setParams (IB.params);

  // Set any non-constant parameter defaults:
  if (!given("L"))
  {
    string msg("Could not find L parameter in instance.");
    std::ostringstream oss;
    oss << "Error in " << netlistLocation() << "\n" << msg;
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, oss.str());
  }

  if (!given("TEMP"))
    temp = getDeviceOptions().temp.dVal();

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
void Instance::registerLIDs(const vector<int> & intLIDVecRef,
                                          const vector<int> & extLIDVecRef)
{
  string msg;
  string tmpstr;

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
  li_Bra = intLIDVec[0];

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline(
"-----------------------------------------------------------------------------");
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << dashedline << endl;

    cout << "::registerLIDs:\n";
    cout << "  name = " << getName() << endl;

    cout << "\nlocal solution indices:\n";
    cout << "  li_Pos = "<< li_Pos << endl;
    cout << "  li_Neg = "<< li_Neg << endl;
    cout << "  li_Bra = "<< li_Bra << endl;

    cout << dashedline << endl;
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
map<int,string> & Instance::getIntNameMap ()
{
  // set up the internal name map, if it hasn't been already.
  if (intNameMap.empty ())
  {
    // set up internal name map
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
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/22/02
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
    msg += "numSta != numStateVars";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
  }

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
const vector< vector<int> > & Instance::jacobianStamp() const
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
void Instance::registerJacLIDs( const vector< vector<int> > & jacLIDVec )
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

  // obtain current accross the inductor
  double current = solVec[li_Bra];

  if( (getSolverState().dcopFlag) && ICGiven )
    current = IC;

  // obtain the "current" value for the flux stored in the inductor
  f0 = L*current;

  // place this value for the charge in the state vector.
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
    // obtain the "current" value for the flux stored in the inductor
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
void Instance::varTypes( vector<char> & varTypeVec )
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
bool Model::processParams (string param)
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
// Purpose       : block constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
Model::Model (const ModelBlock & MB,
                                                SolverState & ss1,
                                                DeviceOptions & do1)
  : DeviceModel(MB,ss1,do1),
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
// Creation Date : 4/04/00
//-----------------------------------------------------------------------------
std::ostream &Model::printOutInstances(std::ostream &os) const
{
  vector<Instance*>::const_iterator iter;
  vector<Instance*>::const_iterator first = instanceContainer.begin();
  vector<Instance*>::const_iterator last  = instanceContainer.end();

  int i, isize;
  isize = instanceContainer.size();

  os << endl;
  os << "Number of Inductor instances: " << isize << endl;
  os << "    name=\t\tmodelName\tParameters" << endl;
  for (i=0, iter=first; iter!=last; ++iter, ++i)
  {
    os << "  " << i << ": " << (*iter)->getName() << "\t";
    os << (*iter)->getModelName();
    os << "\t\tL = " << (*iter)->L;
    os << "\tIC = " << (*iter)->IC;
    os << endl;
  }

  os << endl;

  return os;
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
  for (InstanceVector::const_iterator it = getInstanceVector().begin(); it != getInstanceVector().end(); ++it)
  {
    Instance & inst = *(*it);

    // obtain current accross the inductor
    double current = solVec[inst.li_Bra];

    if( (getSolverState().dcopFlag) && inst.ICGiven )
      current = inst.IC;

    // obtain the "current" value for the flux stored in the inductor
    inst.f0 = inst.L*current;

    // place this value for the charge in the state vector.
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
  for (InstanceVector::const_iterator it = getInstanceVector().begin(); it != getInstanceVector().end(); ++it)
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
  for (InstanceVector::const_iterator it = getInstanceVector().begin(); it != getInstanceVector().end(); ++it)
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

} // namespace Inductor
} // namespace Device
} // namespace Xyce

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

//-----------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_Bsrc.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Robert J Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 06/05/01
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.168.2.3 $
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

#ifdef HAVE_CSTDIO
#include <cstdio>
#else
#include <stdio.h>
#endif

// ----------   Xyce Includes   ----------
#include <N_DEV_Bsrc.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_MatrixLoadData.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>

#include <N_UTL_Expression.h>
#include <N_UTL_BreakPoint.h>

namespace Xyce {
namespace Device {

template<>
ParametricData<Bsrc::Instance>::ParametricData()
{
  // Set up configuration constants:
  setNumNodes(2);
  setNumOptionalNodes(0);
  setNumFillNodes(0);
  setModelRequired(0);

  // Set up double precision variables:
  addPar ("I", 0.0, false, ParameterType::SOLN_DEP,
          &Bsrc::Instance::I,
          NULL, U_AMP, CAT_NONE, "Current for current source");

  addPar ("V", 0.0, false, ParameterType::SOLN_DEP,
          &Bsrc::Instance::V,
          NULL, U_VOLT, CAT_NONE, "Voltage for voltage source");
}

template<>
ParametricData<Bsrc::Model>::ParametricData()
{}

namespace Bsrc {



ParametricData<Instance> &Instance::getParametricData() {
  static ParametricData<Instance> parMap;

  return parMap;
}

ParametricData<Model> &Model::getParametricData() {
  static ParametricData<Model> parMap;

  return parMap;
}

#define Xyce_NONPOINTER_MATRIX_LOAD 1

// Class Instance
//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : "instance block" constructor
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/05/01
//-----------------------------------------------------------------------------
Instance::Instance(InstanceBlock & IB,
                   Model & BMiter,
                   MatrixLoadData & mlData1,
                   SolverState &ss1,
                   ExternData  &ed1,
                   DeviceOptions & do1)
  : DeviceInstance(IB, mlData1, ss1, ed1, do1),
    model_(BMiter),
    IB(IB),
    li_Pos(-1),
    li_Neg(-1),
    li_Bra(-1),
    li_store_branch(-1),
    ABraEquPosNodeOffset(-1),
    ABraEquNegNodeOffset(-1),
    APosEquBraVarOffset(-1),
    ANegEquBraVarOffset(-1),
    fBraEquPosNodePtr(0),
    fBraEquNegNodePtr(0),
    fPosEquBraVarPtr(0),
    fNegEquBraVarPtr(0),
    Exp_ptr(0),
    isVSRC(false),
    scale(1.0),
    nlstep(-1),
    expNumVars(0),
    expBaseVar(0),
    expNumDdt(0),
    expVal(0)
{
  numIntVars   = 1;
  numExtVars   = 2;
  numStateVars = 0;

  setName(IB.getName());
  setModelName(getName());


  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  // Set any non-constant parameter defaults:
  if (given("I") && !given("V"))
  {
    isVSRC = false;
    // current source doesn't have current as part of the solution vector
    // so store it in the store vector
    setNumStoreVars(1);
  }
  else if (!given("I") && given("V"))
  {
    isVSRC = true;
  }
  else
  {
    string msg("Instance::Instance: Must supply one of V= or I=");
    std::ostringstream oss;
    oss << "Error in " << netlistLocation() << "\n" << msg;
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, oss.str());
  }

  if (isVSRC)
  {
    numIntVars = 1;
  }
  else
  {
    numIntVars = 0;
  }

  vector<sDepend>::iterator d;
  vector<sDepend>::iterator begin = dependentParams.begin();
  vector<sDepend>::iterator end = dependentParams.end();

  for  (d = begin ; d != end ; ++d)
  {
    if (d->name == "I" || d->name == "V")
    {
      expNumVars = d->n_vars;
      expBaseVar = d->lo_var;
      Exp_ptr = d->expr;

      expNumDdt = Exp_ptr->getNumDdt();
      ddtVals.resize(expNumDdt);
      li_ddt.resize(expNumDdt);
      numStateVars += expNumDdt;

      expVarDerivs.resize(expNumVars);
      myVarVals.resize(expNumVars);
      break;
    }
  }

  if( jacStamp.empty() )
  {
    if( isVSRC )
    {
      jacStamp.resize(3);
      jacStamp[0].resize(1);
      jacStamp[0][0]=2;
      jacStamp[1].resize(1);
      jacStamp[1][0]=2;
      jacStamp[2].resize(2+expNumVars);
      jacStamp[2][0]=0;
      jacStamp[2][1]=1;
      for( int i = 0; i < expNumVars; ++i )
        jacStamp[2][i+2] = i+3;
    }
    else
    {
      jacStamp.resize( 2 );
      jacStamp[0].resize(expNumVars);
      jacStamp[1].resize(expNumVars);
      for( int i = 0; i < expNumVars; ++i )
      {
        jacStamp[0][i] = i+2;
        jacStamp[1][i] = i+2;
      }
    }
  }

  // Calculate any parameters specified as expressions:

  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:

  processParams();

}

//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/30/03
//-----------------------------------------------------------------------------
bool Instance::processParams (string param)
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/05/01
//-----------------------------------------------------------------------------
Instance::~Instance ()
{
}

// Additional Declarations
//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/21/02
//-----------------------------------------------------------------------------
void Instance::registerLIDs( const vector<int> & intLIDVecRef,
                             const vector<int> & extLIDVecRef )
{
  string msg;

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline =
    "-----------------------------------------------------------------------------";
  cout << endl << dashedline << endl;
  cout << "  BsrcInstance::registerLIDs" << endl;
  cout << "  name = " << getName() << endl;
#endif

  // Check if the size of the ID lists corresponds to the
  // proper number of internal and external variables.
  int numInt = intLIDVecRef.size();
  int numExt = extLIDVecRef.size();

  if ( numInt != numIntVars )
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

#ifdef Xyce_DEBUG_DEVICE
  cout << "  li_Pos = " << li_Pos << endl;
  cout << "  li_Neg = " << li_Neg << endl;
#endif

  if( isVSRC )
  {
    li_Bra = intLIDVec[0];

#ifdef Xyce_DEBUG_DEVICE
    cout << "  li_Bra = " << li_Bra << endl;
#endif
  }

#ifdef Xyce_DEBUG_DEVICE
  cout << dashedline << endl;
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
  if (isVSRC)
  {
    if (intNameMap.empty ())
    {
      // set up internal name map
      string tmpstr(getName()+"_branch");
      spiceInternalName (tmpstr);
      intNameMap[ li_Bra ] = tmpstr;
    }
  }
  return intNameMap;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStoreLIDs
// Purpose       : One store var for device current if this is a current source
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 01/17/2013
//-----------------------------------------------------------------------------
void Instance::registerStoreLIDs(const vector<int> & stoLIDVecRef )
{
  if (!isVSRC)
  {
    string msg;

    // Check if the size of the ID lists corresponds to the
    // proper number of internal and external variables.
    int numSto = stoLIDVecRef.size();

    if (numSto != getNumStoreVars())
    {
      msg = "ResistorInstance::registerStoreLIDs:";
      msg += "numSto != numStoreVars";
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
    }
    li_store_branch = stoLIDVecRef[0];
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::getStoreNameMap
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 12/18/2012
//-----------------------------------------------------------------------------
map<int,string> & Instance::getStoreNameMap ()
{
  if (!isVSRC)
  {
    // set up the internal name map, if it hasn't been already.
    if (storeNameMap.empty ())
    {
      // change subcircuitname:devicetype_deviceName to
      // devicetype:subcircuitName:deviceName
      string modName(getName());
      spiceInternalName(modName);
      string tmpstr;
      tmpstr = modName+":DEV_I";
      storeNameMap[ li_store_branch ] = tmpstr;
    }
  }
  return storeNameMap;
}


//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/21/02
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs( const vector<int> & staLIDVecRef )
{
  string msg;

  // Check if the size of the ID lists corresponds to the
  // proper number of internal and external variables.

  if (staLIDVecRef.size() != numStateVars ||
      li_ddt.size() != expNumDdt ||
      numStateVars != expNumDdt)
  {
    msg = "Instance::registerStateLIDs:";
    msg += "numSta != numStateVars";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
  }
  for (int i=0 ; i<expNumDdt ; ++i)
  {
    li_ddt[i] = staLIDVecRef[i];
  }
}


//-----------------------------------------------------------------------------
// Function      : Instance::getDepSolnVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/06/01
//-----------------------------------------------------------------------------
const vector<string> & Instance::getDepSolnVars()
{
  return DeviceInstance::getDepSolnVars();
}


//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 9/2/02
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
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 9/2/02
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs(
  const vector< vector<int> > & jacLIDVec)
{
  DeviceInstance::registerJacLIDs( jacLIDVec );
  if( isVSRC )
  {
    APosEquBraVarOffset = jacLIDVec[0][0];
    ANegEquBraVarOffset = jacLIDVec[1][0];
    ABraEquPosNodeOffset = jacLIDVec[2][0];
    ABraEquNegNodeOffset = jacLIDVec[2][1];
    ABraEquExpVarOffsets.resize( expNumVars );
    for( int i = 0; i < expNumVars; ++i )
    {
      ABraEquExpVarOffsets[i] = jacLIDVec[2][i+2];
    }
  }
  else
  {
    APosEquExpVarOffsets.resize( expNumVars );
    ANegEquExpVarOffsets.resize( expNumVars );
    for( int i = 0; i < expNumVars; ++i )
    {
      APosEquExpVarOffsets[i] = jacLIDVec[0][i];
      ANegEquExpVarOffsets[i] = jacLIDVec[1][i];
    }

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

  if( isVSRC )
  {
    fPosEquBraVarPtr = &(dFdx[li_Pos][APosEquBraVarOffset]);
    fNegEquBraVarPtr = &(dFdx[li_Neg][ANegEquBraVarOffset]);
    fBraEquPosNodePtr = &(dFdx[li_Bra][ABraEquPosNodeOffset]);
    fBraEquNegNodePtr = &(dFdx[li_Bra][ABraEquNegNodeOffset]);

    fBraEquExpVarPtrs.resize( expNumVars );
    for( int i = 0; i < expNumVars; ++i )
    {
      fBraEquExpVarPtrs[i] = &(dFdx[li_Bra][ ABraEquExpVarOffsets[i] ]);
    }
  }
  else
  {
    fPosEquExpVarPtrs.resize( expNumVars );
    fNegEquExpVarPtrs.resize( expNumVars );
    for( int i = 0; i < expNumVars; ++i )
    {
      fPosEquExpVarPtrs[i]  = &(dFdx[li_Pos][ APosEquExpVarOffsets[i] ]);
      fNegEquExpVarPtrs[i]  = &(dFdx[li_Neg][ ANegEquExpVarOffsets[i] ]);
    }
  }

#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
//
// Purpose       : Calls the expression handler to evaluate the expression
//                 and various derivatives.  These quantities are needed
//                 for the vector and matrix loads.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 06/05/01
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars ()
{
  if (expNumVars == 0)
  {
    if (isVSRC)
    {
      expVal = V;
    }
    else
    {
      expVal = I;
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/05/01
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState ()
{
  bool bsuccess=updateIntermediateVars ();

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

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/05/01
//-----------------------------------------------------------------------------
bool Instance::updateSecondaryState ()
{
  // Get time derivatives from time integrator, and evaluate expression to get
  // derivatives with respect to independent quantities

  if (expNumDdt > 0)
  {
    double * staDerivVec = extData.nextStaDerivVectorRawPtr;

    for (int i=0 ; i<expNumDdt ; ++i)
    {
      ddtVals[i] = staDerivVec[li_ddt[i]];
    }
    Exp_ptr->setDdtDerivs(ddtVals);
  }
  // Evaluate Expression with corrected time derivative values
  if (expNumVars != 0)
  {
    Exp_ptr->evaluate( expVal, expVarDerivs);
  }

  // Test derivatives, if too big, zero out
  for (int i = 0; i < expNumVars; ++i)
  {
    double maxMag = 1.0e+10;
    if (expVarDerivs[i] > maxMag || expVarDerivs[i] < -maxMag)
    {
      cout << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
      cout << "WARNING: Expression derivative exceeded magnitude of ";
      cout << maxMag << "\n";
      cout << "         Value being reduced\n";
      cout << "Bsource is : " << getName() << "\n";
      cout << " Expression derivative evaluates to " << expVarDerivs[i];
      cout << "\n";
      cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
      expVarDerivs[i] = (expVarDerivs[i] > 0) ? maxMag : -maxMag;
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 vsrc instance.
//
// Special Notes : See the special notes for loadDAEFVector.
//
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 04/27/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  double source(0.0), v_pos(0.0), v_neg(0.0), i_bra(0.0);
  double * solVec = extData.nextSolVectorRawPtr;
  double * fVec = extData.daeFVectorRawPtr;
  double * stoVec = extData.nextStoVectorRawPtr;

  source = expVal;

  //VSRC or ISRC
  if (isVSRC)
  {
    // get the value for v_pos, v_neg, i_bra
    v_pos = solVec[li_Pos];
    v_neg = solVec[li_Neg];
    i_bra = solVec[li_Bra];

    double c_tmp = i_bra;
    double v_tmp = (v_pos-v_neg-source);

    fVec[li_Pos] +=  c_tmp;
    fVec[li_Neg] += -c_tmp;
    fVec[li_Bra] +=  v_tmp;
  }
  else
  {
    fVec[li_Pos] += source;
    fVec[li_Neg] += -source;

    stoVec[li_store_branch] = source;
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
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 04/27/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  N_LAS_Matrix & dFdx = *(extData.dFdxMatrixPtr);

  double coef = 1.0;

  if( isVSRC )
  {
    if( getDeviceOptions().scale_src != 0.0 )
    {
      coef *= scale;
    }

    dFdx[li_Pos][APosEquBraVarOffset]  += coef;
    dFdx[li_Neg][ANegEquBraVarOffset]  -= coef;
    dFdx[li_Bra][ABraEquPosNodeOffset] += coef;
    dFdx[li_Bra][ABraEquNegNodeOffset] -= coef;

    for( int i = 0; i < expNumVars; ++i )
    {
      dFdx[li_Bra][ABraEquExpVarOffsets[i]] -= expVarDerivs[i];
    }
  }
  else
  {
    if( expNumVars )
    {
      for( int i = 0; i < expNumVars; ++i )
      {
        dFdx[li_Pos][APosEquExpVarOffsets[i]] += expVarDerivs[i];
        dFdx[li_Neg][ANegEquExpVarOffsets[i]] -= expVarDerivs[i];
      }
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::varTypes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/17/02
//-----------------------------------------------------------------------------
void Instance::varTypes( vector<char> & varTypeVec )
{
  if( !isVSRC )
  {
    varTypeVec.resize(1);
    varTypeVec[0] = 'I';
  }
}

// Class Model

//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/05/01
//-----------------------------------------------------------------------------
Model::Model (const ModelBlock & MB,
              SolverState & ss1,
              DeviceOptions & do1)
  : DeviceModel(MB,ss1,do1)
{
}

//-----------------------------------------------------------------------------
// Function      : Model::~Model
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/05/01
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
// Function      : DeviceModel::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 06/03/02
//-----------------------------------------------------------------------------
bool Model::processParams(std::string param)
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceModel::processInstanceParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 03/23/06
//-----------------------------------------------------------------------------
bool Model::processInstanceParams(std::string param)
{
  return true;
}

// Additional Declarations

//-----------------------------------------------------------------------------
// Function      : Model::printOutInstances
// Purpose       : debugging tool.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/05/01
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
// Bsrc Master functions:

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
    Instance & bi =  *(*it);
    if (bi.expNumVars == 0)
    {
      if (bi.isVSRC)
      {
        bi.expVal = bi.V;
      }
      else
      {
        bi.expVal = bi.I;
        stoVec[bi.li_store_branch]=bi.I;
      }
    }
    // Get values of the arguments for ddt() calls in expression so that the derivatives
    // can be determined by the time integration class
    if (bi.expNumDdt > 0)
    {
      bi.Exp_ptr->getDdtVals (bi.ddtVals);
      for (int j=0 ; j<bi.expNumDdt ; ++j)
      {
        staVec[bi.li_ddt[j]] = bi.ddtVals[j];
      }
    }
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
  for (InstanceVector::const_iterator it = getInstanceVector().begin(); it != getInstanceVector().end(); ++it)
  {
    Instance & bi =  *(*it);

    // Get time derivatives from time integrator, and evaluate expression to get
    // derivatives with respect to independent quantities

    if (bi.expNumDdt > 0)
    {
      for (int j=0 ; j<bi.expNumDdt ; ++j)
      {
        bi.ddtVals[j] = staDerivVec[bi.li_ddt[j]];
      }
      bi.Exp_ptr->setDdtDerivs(bi.ddtVals);
    }
    // Evaluate Expression with corrected time derivative values
    if (bi.expNumVars != 0)
    {
      bi.Exp_ptr->evaluate( bi.expVal, bi.expVarDerivs);
      if (!bi.isVSRC)
      {
        stoVec[bi.li_store_branch]=bi.expVal;
      }
    }

    // Test derivatives, if too big, zero out
    for (int k = 0; k < bi.expNumVars; ++k)
    {
      double maxMag = 1.0e+10;
      if (bi.expVarDerivs[k] > maxMag || bi.expVarDerivs[k] < -maxMag)
      {
        cout << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
        cout << "WARNING: Expression derivative exceeded magnitude of ";
        cout << maxMag << "\n";
        cout << "         Value being reduced\n";
        cout << "Bsource is : " << bi.getName() << "\n";
        cout << " Expression derivative evaluates to " << bi.expVarDerivs[k];
        cout << "\n";
        cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
        bi.expVarDerivs[k] = (bi.expVarDerivs[k] > 0) ? maxMag : -maxMag;
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
  for (InstanceVector::const_iterator it = getInstanceVector().begin(); it != getInstanceVector().end(); ++it)
  {
    Instance & bi = *(*it);
    double v_pos(0.0), v_neg(0.0), i_bra(0.0);
    double source = bi.expVal;

    //VSRC or ISRC
    if (bi.isVSRC)
    {
      // get the value for v_pos, v_neg, i_bra
      v_pos = solVec[bi.li_Pos];
      v_neg = solVec[bi.li_Neg];
      i_bra = solVec[bi.li_Bra];

      double c_tmp = i_bra;
      double v_tmp = (v_pos-v_neg-source);

      fVec[bi.li_Pos] +=  c_tmp;
      fVec[bi.li_Neg] += -c_tmp;
      fVec[bi.li_Bra] +=  v_tmp;
    }
    else
    {
      fVec[bi.li_Pos] += source;
      fVec[bi.li_Neg] += -source;
    }
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
    Instance & bi = *(*it);
    double coef = 1.0;

    if( bi.isVSRC )
    {
      if( getDeviceOptions().scale_src != 0.0 )
      {
        coef *= bi.scale;
      }

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
      *bi.fPosEquBraVarPtr  += coef;
      *bi.fNegEquBraVarPtr  -= coef;
      *bi.fBraEquPosNodePtr += coef;
      *bi.fBraEquNegNodePtr -= coef;

      for( int j = 0; j < bi.expNumVars; ++j )
      {
        *bi.fBraEquExpVarPtrs[j] -= bi.expVarDerivs[j];
      }
#else
      dFdx[bi.li_Pos][bi.APosEquBraVarOffset]  += coef;
      dFdx[bi.li_Neg][bi.ANegEquBraVarOffset]  -= coef;
      dFdx[bi.li_Bra][bi.ABraEquPosNodeOffset] += coef;
      dFdx[bi.li_Bra][bi.ABraEquNegNodeOffset] -= coef;

      for( int j = 0; j < bi.expNumVars; ++j )
      {
        dFdx[bi.li_Bra][bi.ABraEquExpVarOffsets[j]] -= bi.expVarDerivs[j];
      }
#endif
    }
    else
    {
      if( bi.expNumVars )
      {
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
        for( int k = 0; k < bi.expNumVars; ++k )
        {
          *bi.fPosEquExpVarPtrs[k] += bi.expVarDerivs[k];
          *bi.fNegEquExpVarPtrs[k] -= bi.expVarDerivs[k];
        }
#else
        for( int j = 0; j < bi.expNumVars; ++j )
        {
          dFdx[bi.li_Pos][bi.APosEquExpVarOffsets[j]] += bi.expVarDerivs[j];
          dFdx[bi.li_Neg][bi.ANegEquExpVarOffsets[j]] -= bi.expVarDerivs[j];
        }
#endif
      }
    }
  }
  return true;
}

} // namespace Bsrc
} // namespace Device
} // namespace Xyce

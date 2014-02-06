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
// Filename       : $RCSfile: N_DEV_MutIndNonLin2.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Rich Schiek, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/21/2005
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.28.2.3 $
//
// Revision Date  : $Date: 2013/10/03 17:23:38 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// #define Xyce_NO_MUTIND_MASK 1
// ---------- Standard Includes ----------

#include <N_UTL_Misc.h>

#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif

#include <fstream>
#include <algorithm>
#include <vector>
#include <set>

// ----------   Xyce Includes   ----------

#include <N_DEV_MutIndNonLin2.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_MatrixLoadData.h>

//This contains important constants like permitivity of free space
#include <N_DEV_Const.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>

#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RefCountPtr;
using Teuchos::rcp;

namespace Xyce {
namespace Device {

template<>
ParametricData<MutIndNonLin2::Instance>::ParametricData()
{
  getConfigTable().primaryParameter = "";
  getConfigTable().modelTypes.clear();
    // Set up configuration constants:
    setNumNodes(2);
    setNumOptionalNodes(0);
    setNumFillNodes(0);
    setModelRequired(1);
    addModelType("CORE");

    // Set up double precision variables:
    addPar ("COUP_VAL",   1.0, false, ParameterType::NO_DEP,
      &MutIndNonLin2::Instance::mutualCup,
      &MutIndNonLin2::Instance::mutualCupGiven,
      U_NONE, CAT_NONE, "Coupling coefficient");

    addPar ("NONLINEARCOUPLING", 0.0, false, ParameterType::NO_DEP,
      &MutIndNonLin2::Instance::nonlinFlag,
      &MutIndNonLin2::Instance::nonlinFlagGiven,
      U_NONE, CAT_NONE, "Nonlinear coupling flag");

    // Set up non-double precision variables:
    addPar ("COUPLEDMutIndNonLin", std::vector<string>(), false, ParameterType::NO_DEP,
            &MutIndNonLin2::Instance::inductorNames, NULL, U_NONE, CAT_NONE, "" );
    addPar ("COUPLEDINDUCTANCE", std::vector<double>(), false, ParameterType::NO_DEP,
            &MutIndNonLin2::Instance::inductorInductances, NULL, U_NONE, CAT_NONE, "");
    addPar ("NODE1", std::vector<string>(), false, ParameterType::NO_DEP,
            &MutIndNonLin2::Instance::inductorsNode1, NULL, U_NONE, CAT_NONE, "");
    addPar ("NODE2", std::vector<string>(), false, ParameterType::NO_DEP,
            &MutIndNonLin2::Instance::inductorsNode2, NULL, U_NONE, CAT_NONE, "");
    addPar ("COUPLING", std::vector<double>(), false, ParameterType::NO_DEP,
            &MutIndNonLin2::Instance::couplingCoefficient, NULL, U_NONE, CAT_NONE, "Coupling coefficient");
    addPar ("COUPLEDINDUCTOR", std::vector<string>(), false, ParameterType::NO_DEP,
            &MutIndNonLin2::Instance::couplingInductor, NULL, U_NONE, CAT_NONE, "");
  }

template<>
ParametricData<MutIndNonLin2::Model>::ParametricData()
{
    addPar ("A", 1000.0, false, ParameterType::NO_DEP,
      &MutIndNonLin2::Model::A,
    NULL, U_AMPMM1, CAT_MATERIAL, "Thermal energy parameter");

    addPar ("AREA", 0.1, false, ParameterType::NO_DEP,
      &MutIndNonLin2::Model::Area,
    NULL, U_CM2, CAT_GEOMETRY, "Mean magnetic cross-sectional area");

    addPar ("ALPHA", 5.0e-5, false, ParameterType::NO_DEP,
      &MutIndNonLin2::Model::Alpha,
    NULL, U_NONE, CAT_GEOMETRY, "Domain coupling parameter");

    addPar ("BETAH", 0.0001, false, ParameterType::NO_DEP,
      &MutIndNonLin2::Model::BetaH,
    NULL, U_NONE, CAT_NONE, "Modeling constant");

    addPar ("BETAM", 3.125e-5, false, ParameterType::NO_DEP,
      &MutIndNonLin2::Model::BetaM,
    NULL, U_NONE, CAT_NONE, "Modeling constant");

    addPar ("C", 0.2, false, ParameterType::NO_DEP,
      &MutIndNonLin2::Model::C,
    NULL, U_NONE, CAT_MATERIAL, "Domain flesing parameter");

    addPar ("DELV", 0.1, false, ParameterType::NO_DEP,
      &MutIndNonLin2::Model::DeltaV,
    NULL, U_VOLT, CAT_NONE, "Smoothing coefficient for voltage difference over first inductor");

    addPar ("GAP", 0.0, false, ParameterType::NO_DEP,
      &MutIndNonLin2::Model::Gap,
    NULL, U_CM, CAT_GEOMETRY, "Effective air gap");

    addPar ("K", 500.0, false, ParameterType::NO_DEP,
      &MutIndNonLin2::Model::Kirr,
    NULL, U_AMPMM1, CAT_MATERIAL, "Domain anisotropy parameter");

    addPar ("KIRR", 500.0, false, ParameterType::NO_DEP,
      &MutIndNonLin2::Model::Kirr,
    NULL, U_AMPMM1, CAT_MATERIAL, "Domain anisotropy parameter");

    addPar ("MS", 1.0e+6, false, ParameterType::NO_DEP,
      &MutIndNonLin2::Model::Ms,
    NULL, U_AMPMM1, CAT_MATERIAL, "Saturation magnetization");

    addPar ("LEVEL", 0.0, false, ParameterType::NO_DEP,
      &MutIndNonLin2::Model::LevelIgnored,
    NULL, U_NONE, CAT_NONE, "for pspice compatibility -- ignored");

    addPar ("PACK", 0.0, false, ParameterType::NO_DEP,
      &MutIndNonLin2::Model::PackIgnored,
    NULL, U_NONE, CAT_NONE, "for pspice compatibility -- ignored");

    addPar ("PATH", 1.0, false, ParameterType::NO_DEP,
      &MutIndNonLin2::Model::Path,
    NULL, U_CM, CAT_GEOMETRY, "Total mean magnetic path");

    addPar ("VINF", 1.0, false, ParameterType::NO_DEP,
      &MutIndNonLin2::Model::Vinf,
    NULL, U_VOLT, CAT_NONE, "Smoothing coefficient for voltage difference over first inductor");

    addPar ("TNOM", 27.0, false, ParameterType::NO_DEP,
      &MutIndNonLin2::Model::tnom,
    NULL, U_DEGC, CAT_MATERIAL, "Reference temperature");

    addPar ("TC1", 0.0, false, ParameterType::NO_DEP,
      &MutIndNonLin2::Model::tempCoeff1,
    NULL, U_NONE, CAT_MATERIAL, "First order temperature coeff.");

    addPar ("TC2", 0.0, false, ParameterType::NO_DEP,
      &MutIndNonLin2::Model::tempCoeff2,
    NULL, U_NONE, CAT_MATERIAL, "Second order temperature coeff.");

    addPar ("PZEROTOL", 0.1, false, ParameterType::NO_DEP,
      &MutIndNonLin2::Model::pZeroTol,
    NULL, U_NONE, CAT_NONE, "Tolerance for nonlinear zero crossing");

    addPar ("MVARSCALING", 1.0, false, ParameterType::NO_DEP,
      &MutIndNonLin2::Model::mVarScaling,
    NULL, U_NONE, CAT_NONE, "M-variable scaling.");

    addPar ("RVARSCALING", 1.0, false, ParameterType::NO_DEP,
      &MutIndNonLin2::Model::rVarScaling,
    NULL, U_NONE, CAT_NONE, "R-variable scaling");

    addPar ("MEQNSCALING", 1.0, false, ParameterType::NO_DEP,
      &MutIndNonLin2::Model::mEqScaling,
    NULL, U_NONE, CAT_NONE, "M-equation scaling");

    addPar ("REQNSCALING", 1.0, false, ParameterType::NO_DEP,
      &MutIndNonLin2::Model::rEqScaling,
    NULL, U_NONE, CAT_NONE, "R-equation scaling");

    // Set up non-double precision variables:
    addPar ("OUTPUTSTATEVARS", 0, false, ParameterType::NO_DEP,
      &MutIndNonLin2::Model::outputStateVars, NULL, U_NONE, CAT_NONE, "Flag to save state variables" );
    addPar ("INCLUDEDELTAM", 0, false, ParameterType::NO_DEP,
      &MutIndNonLin2::Model::includeDeltaM,
      &MutIndNonLin2::Model::includeDeltaMGiven, U_NONE, CAT_NONE, "Flag to make M calculation implicit" );
    addPar ("USERKINTEGRATION", 0, false, ParameterType::NO_DEP,
      &MutIndNonLin2::Model::useRKIntegration,
      &MutIndNonLin2::Model::useRKIntegrationGiven, U_NONE, CAT_NONE, "Flag to use 4th order Runge-Kutta integration for dM/dH" );
    addPar ("USESTATEDERIV", 0, false, ParameterType::NO_DEP,
      &MutIndNonLin2::Model::useStateDeriv, NULL, U_NONE, CAT_NONE, "Flag to use state vector for derivatives" );
    addPar ("VOLTAGELIMITERFLAG", 0, false, ParameterType::NO_DEP,
      &MutIndNonLin2::Model::voltageLimiterFlag,
      NULL, U_NONE, CAT_NONE, "Flag to use voltage limiting on Mag and R internal variables" );
    addPar ("MAGLIMITTHRES", 0.1, false, ParameterType::NO_DEP,
      &MutIndNonLin2::Model::magLimitThres,
      NULL, U_NONE, CAT_NONE, "Threshold over which newton interation changes in Mag are limited." );
    addPar ("RLIMITTHRES", 0.1, false, ParameterType::NO_DEP,
      &MutIndNonLin2::Model::rLimitThres,
      NULL, U_NONE, CAT_NONE, "Threshold over which newton interation changes in R are limited." );
    addPar ("FACTORMS", 0, false, ParameterType::NO_DEP,
      &MutIndNonLin2::Model::factorMS,
      NULL, U_NONE, CAT_NONE, "Flag to save state variables" );
}

namespace MutIndNonLin2 {




ParametricData<Instance> &Instance::getParametricData() {
  static ParametricData<Instance> parMap;

  return parMap;
}

ParametricData<Model> &Model::getParametricData() {
  static ParametricData<Model> parMap;

  return parMap;
}

//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
Instance::Instance(InstanceBlock & IB,
                                Model & Iiter,
                                MatrixLoadData & mlData1,
                                SolverState &ss1,
                                ExternData  &ed1,
                                DeviceOptions & do1)
  : DeviceInstance (IB, mlData1, ss1, ed1, do1),
    model_(Iiter),
    temp(getDeviceOptions().temp.dVal()),
    P(0.0),
    dP_dM(0.0),
    dP_dBranchCurrentSum(0.0),
    dP_dV1Pos(0.0),
    dP_dV1Neg(0.0),
    branchCurrentSum(0.0),
    mEquFval(0.0),
    MagVar(0.0),
    oldBranchCurrentSum(0.0),
    MagVarUpdate(0.0),
    lastMagUpdate(0.0),
    PPreviousStep(0.0),
    includeDeltaM(false),
    useRKIntegration(false),
    outputStateVarsFlag( false )
{
#ifdef Xyce_DEBUG_DEVICE
  cout << "In Instance constructor" << endl;
#endif

  // for a simple case of 2 leads, we have 3 internal vars (I_branch, H, M)
  // and one state var (I_branch)
  numExtVars   = 2;
  numIntVars   = 3;
  numStateVars = 0;
  tempGiven    = false;

  setName(IB.getName());
  const int ibev = IB.numExtVars;
  const int ibiv = IB.numIntVars;
  setModelName(model_.getName());


#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
    {
      cout << "Instance::Instance() " << endl;
    }
#endif

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  // these two vectors are used for 4th order RK
  // the vectors are one shorted than one would think you would need because
  // current step is held by local variables (thus, this just data from steps
  // n-1, n-2 and n-3
  branchCurrentSumHistory.resize(3);
  PFunctionHistory.resize(3);

  // Set any non-constant parameter defaults:
#ifdef Xyce_DEBUG_DEVICE
  //if (getDeviceOptions().debugLevel > 0)
  //{
    //cout << "\n\nInstance Params:\n";
//    outputParams(0);
  //}
#endif

  // if the model card askes for the delta M calculation to be implicit, then
  // change the includeDeltaM flag
  if( model_.includeDeltaMGiven && (model_.includeDeltaM > 0))
  {
    includeDeltaM = true;
  }

  if( model_.useRKIntegrationGiven && (model_.useRKIntegration > 0))
  {
    useRKIntegration = true;
  }

  // now load the instance data vector
  for( int i=0; i<inductorNames.size(); ++i )
  {
    InductorInstanceData * inductorData = new InductorInstanceData();
    inductorData->name = inductorNames[i];
    inductorData->L = inductorInductances[i];
    inductorData->baseL = inductorInductances[i];
    inductorData->ICGiven = false;
    inductorData->inductorCurrentOffsets.resize( inductorNames.size() );

    instanceData.push_back( inductorData );
  }
  numInductors = instanceData.size();

  inductorCurrents.resize( numInductors );
  inductanceVals.resize( numInductors );
  LOI.resize( numInductors );
  LO.resize( numInductors );
  for( int i=0; i<numInductors; ++i)
  {
    LO[i].resize( numInductors );
  }
  dHe_dI.resize(numInductors);
  dManp_dI.resize(numInductors);
  ddelM_dI.resize(numInductors);
  dMirrp_dI.resize(numInductors);
  dP_dI.resize( numInductors );

  updateInductanceMatrix();

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:
  processParams ();

  // if the user has requested output of the internal vars H and M
  // then open a file for that output.
  if( model_.outputStateVars > 0 )
  {
    outputStateVarsFlag = true;
    std::string filename( "Inductor_" );
    filename.append( getName() );
    filename.append( ".dat" );
    // convert any embeded ':' or '%' characters to '_'
    replace( filename.begin(), filename.end(), '%', '_' );
    replace( filename.begin(), filename.end(), ':', '_' );

    outputFileStreamPtr = rcp( new ofstream() );
    outputFileStreamPtr->open( filename.c_str() );
    if( !(*outputFileStreamPtr) )
    {
      string msg("Instance constructor.\n");
      msg += "\tCould not open file for output of state variables. name =" + getName();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
    }
    (*outputFileStreamPtr).setf(ios::scientific, ios::floatfield );
    (*outputFileStreamPtr).width(20);
    (*outputFileStreamPtr).precision(12);
  }

  // update internal/external/state variable counts
  numExtVars = 2*numInductors;       // 2 nodes Vin, Vout per inductor
  numIntVars = numInductors;     // branch current per inductor
  // if we're including the deltaM and deltaHapp variables, then there will be two more internal vars
  if(includeDeltaM)
  {
    numIntVars+=1;
  }

  // allocate space for InductorOffsets
  deltaMEquInductorOffsets.resize(numInductors);

  // set up the jacobian stamp
  // for an individual inductor with the two interal variables would be
  //
  //          V1   V2   Ib
  //  kcl1               1
  //  kcl2              -1
  //  branch  1    -1   L/dt
  //
  //  for a collection of these, the internal variable, branch equations,
  //  must be at the end of a given stamp row as well as the internal
  //  vars for M and R in this non-linear version.
  //
  //  So for N inductors the samp is:
  //
  //          V1  V2  V3  V4 ... V2N  I1  I2  ... IN   M
  //  kcl1                             1
  //  kcl2                            -1
  //  kcl3                                 1
  //  kcl4                                -1
  //  branch1 1   -1                 L/dt  c  ... c
  //  branch2          1  -1          c  L/dt ... c
  //  delta M                          x    x  ... x    x
  //
  //  where "c" is an induced current change and "x" are
  //  values which must be computed.

  jacStamp.resize( numExtVars + numIntVars);

  for( int i=0; i< numInductors; ++i )
  {
    //
    // allocate space
    //
    // kcl V+ node
    jacStamp[2*i].resize(1);
    // kcl V- node
    jacStamp[2*i+1].resize(1);
    if( i == 0 )
    {
      jacStamp[2*numInductors].resize(numInductors + 2);
    }
    else
    {
      jacStamp[2*numInductors + i].resize(numInductors + 2);
    }

    //
    // fill in dependency
    //
    // kcl V+ node
    jacStamp[2*i  ][0] = 2*numInductors + i;
    // kcl V- node
    jacStamp[2*i+1][0] = 2*numInductors + i;

    if( i==0 )
    {
      jacStamp[2*numInductors ][0] = 0;
      jacStamp[2*numInductors ][1] = 1;
      for( int j=0; j<numInductors; ++j )
      {
        jacStamp[2*numInductors][j+2] = 2*numInductors + j;
      }
    }
    else
    {
      jacStamp[2*numInductors + i][0] = 2*i;
      jacStamp[2*numInductors + i][1] = 2*i + 1;
      for( int j=0; j<numInductors; ++j )
      {
        jacStamp[2*numInductors + i][j+2] = 2*numInductors + j;
      }
    }
  }

  if(includeDeltaM)
  {
    // now the deltaM equation
    jacStamp[ 3*numInductors ].resize(numInductors + 1);
    // deltaM offsets to I_1 ... I_n and deltaM
    for(int i=0; i<numInductors; ++i)
    {
      jacStamp[ 3*numInductors ][i]=2*numInductors+i;
    }
    jacStamp[ 3*numInductors ][numInductors] = 3*numInductors;

  }
#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "Instance::Instance----------" << endl;
    cout << "numExtVars = " << numExtVars << ", " << ibev << endl
      << "numIntVars = " << numIntVars << ", " << ibiv << endl
      << "numStateVars = " << numStateVars << endl
      << "numInductors = " << numInductors << endl
      << "jacStamp = " << endl;
    for( int i = 0; i<jacStamp.size(); ++i )
    {
      cout << "jacStamp[ " << i << " ] = { ";
      for( int j=0; j<jacStamp[i].size(); ++j)
      {
        cout << jacStamp[i][j];
        if( j != ( jacStamp[i].size() -1 ) )
        {
          cout << ", ";
        }
      }
      cout << " }" << endl;
    }
  }
#endif
}


//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
Instance::~Instance()
{
  // Close output file if we opened one
  if( outputStateVarsFlag && outputFileStreamPtr->is_open() )
  {
    outputFileStreamPtr->close();
    if( !(*outputFileStreamPtr) )
    {
      string msg("Instance destructor.\n");
      msg += "\tCould not close file for output of state variables. name =" + getName();
         N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
    }
    // since outputFileStreamPtr is a ref counted pointer
    // we don't need to delete it.
  }

  vector< InductorInstanceData* >::iterator currentInductor = instanceData.begin();
  vector< InductorInstanceData* >::iterator endInductor = instanceData.end();
  for ( ; currentInductor != endInductor ; ++currentInductor)
  {
    if (*currentInductor != NULL)
    {
      delete *currentInductor;
      *currentInductor = NULL;
    }
  }
  instanceData.clear();
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
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
  // get the current values of the inductances and currentOffsets
  vector< InductorInstanceData* >::iterator currentInductor = instanceData.begin();
  vector< InductorInstanceData* >::iterator endInductor = instanceData.end();
  int i = 0;
  int j = 0;
  while( currentInductor != endInductor )
  {
    (*currentInductor)->li_Pos = extLIDVec[ i++ ];
    (*currentInductor)->li_Neg = extLIDVec[ i++ ];
    (*currentInductor)->li_Branch = intLIDVec[ j++ ];
    currentInductor++;
  }

  if(includeDeltaM)
  {
    // now get deltaHapp and deltaM
    //li_deltaHappVar = intLIDVec[ j++ ];
    li_deltaMagVar  = intLIDVec[ j++ ];
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "Instance::registerLIDs------------------------" << endl;
    currentInductor = instanceData.begin();
    i=0;
    while( currentInductor != endInductor )
    {
      cout << "Inductor [ " << i++ << " ] "
           << "   li_Pos = " << (*currentInductor)->li_Pos
           << "   li_Neg = " << (*currentInductor)->li_Neg
           << "   li_Branch = " << (*currentInductor)->li_Branch << endl;
      currentInductor++;
    }
    if(includeDeltaM)
    {
      cout << " li_deltaMagVar = " << li_deltaMagVar << endl;
           //<< " li_deltaHappVar = " << li_deltaHappVar << endl;
    }
  }
#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::getIntNameMap
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
map<int,string> & Instance::getIntNameMap()
{
  // set up the internal name map, if it hasn't been already.
  if (intNameMap.empty ())
  {
    string baseString(getName() + "_");
    string tempString;
    vector< InductorInstanceData* >::iterator currentInductor = instanceData.begin();
    vector< InductorInstanceData* >::iterator endInductor = instanceData.end();
    int i = 0;
    int j = 0;
    while( currentInductor != endInductor )
    {
      tempString = baseString + (*currentInductor)->name +"_branch";
      spiceInternalName (tempString);
      intNameMap[ (*currentInductor)->li_Branch ] = tempString;
      currentInductor++;
    }
  }

  return intNameMap;
}


//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
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

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "Instance::registerStateLIDs-------------------" << endl;
  }
#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
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
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const vector< vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );
#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "Instance::registerJacLIDs ----------------------------" << endl;

    cout << "jacLIDVec = " << endl;
    for( int i = 0; i<jacStamp.size(); ++i )
    {
      cout << "jacLIDVec[ " << i << " ] = { ";
      for( int j=0; j<jacLIDVec[i].size(); ++j)
      {
        cout << jacLIDVec[i][j];
        if( j != ( jacLIDVec[i].size() -1 ) )
        {
          cout << ", ";
        }
      }
      cout << " }" << endl;
    }
  }
#endif

  vector< InductorInstanceData* >::iterator currentInductor = instanceData.begin();
  vector< InductorInstanceData* >::iterator endInductor = instanceData.end();
  // int numInductors = instanceData.size();  // don't need this as it's defined at class level
  int i = 0;
  while( currentInductor != endInductor )
  {
    (*currentInductor)->APosEquBraVarOffset  = jacLIDVec[ 2*i     ][ 0 ];
    (*currentInductor)->ANegEquBraVarOffset  = jacLIDVec[ 2*i + 1 ][ 0 ];
    (*currentInductor)->vPosOffset = jacLIDVec[ 2*numInductors + i ][ 0 ];
    (*currentInductor)->vNegOffset = jacLIDVec[ 2*numInductors + i ][ 1 ];

    (*currentInductor)->ABraEquPosNodeOffset = jacLIDVec[ 2*numInductors + i ][ 0  ];
    (*currentInductor)->ABraEquNegNodeOffset = jacLIDVec[ 2*numInductors + i ][ 1  ];
    for( int j=0; j<numInductors; ++j )
    {
      if( i == j )
      {
        (*currentInductor)->ABraEquBraVarOffset  = jacLIDVec[ 2*numInductors + i ][ j + 2  ];
      }
      (*currentInductor)->inductorCurrentOffsets[ j ] = jacLIDVec[ 2*numInductors + i ][ j + 2  ];
    }

    //NoMag (*currentInductor)->magOffset = jacLIDVec[ 2*numInductors + i ][ numInductors + 2 + extraOffset ];
    currentInductor++;
    i++;
  }

  if(includeDeltaM)
  {
    // now get the deltaM
    for( int i=0; i<numInductors; i++)
    {
      deltaMEquInductorOffsets[i] = jacLIDVec[ 3*numInductors ][ i ];
    }
    deltaMEquDeltaMOffset = jacLIDVec[ 3*numInductors ][ numInductors ];
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    currentInductor = instanceData.begin();
    i=0;
    while( currentInductor != endInductor )
    {
      cout << "Inductor [ " << i << " ] " << (*currentInductor)->name << endl
           << "   APosEquBraVarOffset = " << (*currentInductor)->APosEquBraVarOffset << endl
           << "   ANegEquBraVarOffset = " << (*currentInductor)->ANegEquBraVarOffset << endl
           << "   vPosOffset = " << (*currentInductor)->vPosOffset << endl
           << "   vNegOffset = " << (*currentInductor)->vNegOffset << endl
           << "   ABraEquPosNodeOffset = " << (*currentInductor)->ABraEquPosNodeOffset << endl
           << "   ABraEquNegNodeOffset = " << (*currentInductor)->ABraEquNegNodeOffset << endl
           << "   ABraEquBraVarOffset = " << (*currentInductor)->ABraEquBraVarOffset << endl
           << "   magOffset = " << (*currentInductor)->magOffset << endl;
      cout << "\tInductor branch offsets = { ";
      for( int j=0; j<numInductors ; ++j )
      {
        cout << (*currentInductor)->inductorCurrentOffsets[ j ] << ", ";
      }
      cout << "} " << endl;
      i++;
      currentInductor++;
    }
    cout << endl;

    if(includeDeltaM)
    {
      cout << "deltaMEquInductorOffsets = ";
      for( int i=0; i<numInductors; i++ )
      {
        cout << deltaMEquInductorOffsets[i] << " ";
      }
      cout //<< "deltaMEquDeltaHappOffset = " << deltaMEquDeltaHappOffset
        << " deltaMEquDeltaMOffset = " << deltaMEquDeltaMOffset << endl;
    }
  }
#endif
}


//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
bool Instance::processParams(string param)
{
  // now set the temperature related stuff.
  if (tempGiven)
  {
    updateTemperature(temp);
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateTemperature
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
bool Instance::updateTemperature ( const double & temp)
{
  bool bsuccess = true;

  // current temp difference from reference temp.
  double difference = temp - model_.tnom;

  vector< InductorInstanceData* >::iterator currentData = instanceData.begin();
  while( currentData != instanceData.end() )
  {
    double factor = 1.0 + (model_.tempCoeff1)*difference +
                          (model_.tempCoeff2)*difference*difference;
    (*currentData)->L = ((*currentData)->baseL)*factor;
    currentData++;
  }

  // now that the inductances have changed we need to update the matrix.
  updateInductanceMatrix();

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       : updates a set of common variables used by RHS and jacobian
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars ()
{
  bool bsuccess = true;
#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0  && getSolverState().debugTimeFlag)
  {
    cout << "Instance::updateIntermediateVars" << endl;
  }
#endif
  N_LAS_Vector & solVector = *(extData.nextSolVectorPtr);
  N_LAS_Vector & lastSolVector = *(extData.lastSolVectorPtr);

  // some parameters in the model class that we will use often
  const double A      = model_.A;
  const double Alpha  = model_.Alpha;
  const double Area   = model_.Area;
  const double BetaH  = model_.BetaH;
  const double BetaM  = model_.BetaM;
  const double C      = model_.C;
  const double DeltaV = model_.DeltaV;
  const double Gap    = model_.Gap;
  const double Ms     = model_.Ms;
  const double Kirr   = model_.Kirr;
  const double Path   = model_.Path;
  const double Vinf   = model_.Vinf;

  double latestMag; //NoMag  = solVector[ li_MagVar ];

  //sum of currents through the inductors
  branchCurrentSum = 0.0;
  for(int i=0; i<numInductors; i++)
  {
    branchCurrentSum += solVector[instanceData[i]->li_Branch] * inductanceVals[ i ];
  }

  latestMag = MagVar + MagVarUpdate;

  // used in voltage drop over first inductor
  double V1Pos = solVector[(instanceData[0])->li_Pos];
  double V1Neg = solVector[(instanceData[0])->li_Neg];

  double qV = (DeltaV / Vinf) * (V1Pos - V1Neg);

  double tanh_qV = 0.0;
  if (fabs(qV) < CONSTTANH_THRESH)
  {
    tanh_qV = tanh(qV);
  }
  else if (qV < 0.0)
  {
    tanh_qV = -1.0;
  }
  else
  {
    tanh_qV = 1.0;
  }

  double Happ = branchCurrentSum / Path;

#ifdef MS_FACTORING2
  double H = Happ - (Gap / Path) * latestMag * Ms;
  double He = H + Alpha * latestMag * Ms;
#else
  double H = Happ - (Gap / Path) * latestMag;
  double He = H + Alpha * latestMag;
#endif

  double Heo = BetaH*A;

  // terms that come up frequently
  double gap_path = Gap / Path;
  double He2 = He*He;
  double Heo2 = Heo*Heo;
  double sq_Heo2He2 = sqrt(Heo2 + He2);

  double delM0 = model_.BetaM * Ms;
  double Man = Ms * He / ( A + sq_Heo2He2 );
#ifdef MS_FACTORING2
  double delM = Man - latestMag*Ms;
#else
  double delM = Man - latestMag;
#endif

  double delM2 = delM*delM;
  double delM02 = delM0*delM0;
  double sq_delM02delM2 = sqrt( delM02 + delM2 );

  double Pold = P;

  //Mirrp = (delM * tanh_qV + sq_delM02delM2 ) / (2*( Kirr- Alpha * Ms * sq_delM02delM2));
  //Manp =  (A + Heo2/sq_Heo2He2) / pow(A + sq_Heo2He2, 2.0);
  //P = ( C * Manp + (1 - C) * Mirrp) / (1 + (gap_path - Alpha) * C * Ms * Manp + gap_path * (1-C) * Ms * Mirrp);
 #ifdef MS_FACTORING2
  double Mirrp = (delM * tanh_qV + sq_delM02delM2 ) / (2*( Kirr- Alpha * sq_delM02delM2));
  double Manp =  Ms * (A + Heo2/sq_Heo2He2) / pow(A + sq_Heo2He2, 2.0);
  P = ( C * Manp + (1 - C) * Mirrp) / ((1 + (gap_path - Alpha) * C * Manp + gap_path * (1-C) * Mirrp)*Ms);
 #else
  double Mirrp = (delM * tanh_qV + sq_delM02delM2 ) / (2*( Kirr- Alpha * sq_delM02delM2));
  double Manp =  Ms*(A + Heo2/sq_Heo2He2) / pow(A + sq_Heo2He2, 2.0);
  P = ( C * Manp + (1 - C) * Mirrp) / (1 + (gap_path - Alpha) * C * Manp + gap_path * (1-C) * Mirrp);
#endif

  /*
  // this works to force tigher convergence tolerance on P and thus delta M as well.
  // doesn't help in addressing the difference in M between the level 1 and level 2
  // mutual inductors.  Must be in the integration method then.
  origFlag = true;
  //std::cout << "newton itr, time, P - Pold = " << getSolverState().newtonIter << ", " << getSolverState().currTime << ", " << (P-Pold) << ", " << (branchCurrentSum - oldBranchCurrentSum);
  if (fabs(P - Pold) > 0.005 )
  {
    // P changed too much (due to changes in M and I)
    // take this as a limiting step
    // note that we have started limiting
    origFlag = false;
    //std::cout << "limiting ";
    //P = Pold;
  }
  //std::cout << std::endl;
  */

  // at this point we have P so now we can update mag.
  /*
  The problem is that if deltaM is too big, then we need to shrink the time step.  One
  way to control this is to set a max time step.  But what we really need to do is calculate
  deltaM and then if it's over some fraction of Ms then turn on the limiting flag (or bail on the step
  but I think turning on limiting is safer and if we hit maxItter with it on then we'll get that step
  rejected.
  */

  if( useRKIntegration )
  {
    // use 4th order runga-kutta to estimate MagVarUpdate
    double stepLen = branchCurrentSumHistory[0] + branchCurrentSumHistory[1] + branchCurrentSumHistory[2] + (branchCurrentSum - oldBranchCurrentSum);
    MagVarUpdate = stepLen * (  PFunctionHistory[0] +
                    2*PFunctionHistory[1] +
                    2*PFunctionHistory[2] +
                      P) / 6;
  }
  else
  {
    // forward euler method
    MagVarUpdate = P * (branchCurrentSum - oldBranchCurrentSum) / model_.Path;

    // trap
    double MagVarUpdateWithTrap = 0.5 * (P + PPreviousStep) * (branchCurrentSum - oldBranchCurrentSum) / model_.Path;
    origFlag = true;
    if( fabs( MagVarUpdate ) > 0.25 * Ms )
    {
      // step was too big, so
      // turn on limiting
      origFlag = false;
    }
  }

  // try dH/dt to update mag
//   double MagVarUpdateByTime = P * (branchCurrentSum - oldBranchCurrentSum) / (model_.Path);
//   cout << "time, update, updatebyTime = " << getSolverState().currTime << ", " << MagVarUpdate << ", " << MagVarUpdateByTime << endl;
//   if( fabs(MagVarUpdateByTime) < fabs(MagVarUpdate) )
//   {
//     MagVarUpdate = MagVarUpdateByTime;
//   }

  latestMag = MagVar + MagVarUpdate;

  if(includeDeltaM)
  {
    // in this case we're scaling MagVarUpdate by Ms because it's being solved
    // with the full system and big changes in
    MagVarUpdate /=model_.Ms;
  }

  double dP_Denom = 1.0 + (gap_path - Alpha)*C*Manp + gap_path * (1.0-C) * Mirrp;

  // now get dP_dI for each inductor
  for( int i=0; i<numInductors; i++)
  {

    // old way

    //    dHe_dI = (1/lt)*N; % vector
    dHe_dI[ i ] = inductanceVals[ i ] / Path;

    //    dManp_dI = ( ( (-Ms*Heo^2*(Heo^2+He^2)^(-3/2)*He*dHe_dI) * (A + sqrt(Heo^2+He^2))^2 ) - ...
    //                ( ((A+sqrt(Heo^2+He^2))*(Heo^2+He^2)^(-1/2)*2*He*dHe_dI) * (Ms*(A + Heo^2/(sqrt(Heo^2+He^2)))) ) ...
    //              ) / ...
    //          (A + sqrt(Heo^2+He^2))^4; % vector

    /*
    dManp_dI[i] = dHe_dI[i] * ( ( (-Ms * Heo2 * pow((Heo2+He2),(-3.0/2.0)) * He) * pow((A + sqrt(Heo2 + He2)),2.0) ) -
                    ( ( (A + sqrt(Heo2 + He2)) * (1/sqrt(Heo2 + He2)) * 2 * He) *
                      (Ms * (A + Heo2/(sqrt(Heo2 + He2)))) ) ) /
                  pow((A + sqrt(Heo2 + He2)),4.0);
    */
    dManp_dI[i] = ( -Ms * He / (pow(A + sq_Heo2He2, 2.0)*sq_Heo2He2)) *
                   ( (Heo2 / (Heo2 + He2)) + (2.0*(A + Heo2 / sq_Heo2He2)/(A+sq_Heo2He2)) ) * dHe_dI[i];
    // ddelM_dI = ( (Ms*dHe_dI) * (A+sqrt(Heo^2+He^2)) - ...
    //            (Ms*He) * ((Heo^2+He^2)^(-1/2)*He*dHe_dI) ...
    //          ) / ...
    //      (A + sqrt(Heo^2+He^2))^2; % vector

    //ddelM_dI[i] = (Ms / (A + sq_Heo2He2)) * (1.0 - He2/sq_Heo2He2) * dHe_dI[i];
    //ddelM_dI[i] = (Ms / (A + sq_Heo2He2)) * (1.0 - He2/((A + sq_Heo2He2)*sq_Heo2He2)) * dHe_dI[i] - dM_dI;
    ddelM_dI[i] = (Ms / (A + sq_Heo2He2)) * (1.0 - He2/((A + sq_Heo2He2)*sq_Heo2He2)) * dHe_dI[i];
    // dMirrp_dI = ( (ddelM_dI*tanh(qV)+(delM0^2+delM^2)^(-1/2)*delM*ddelM_dI) * (2*(Kirr-alpha*sqrt(delM0^2+delM^2))) - ...
    //                (-alpha*(delM0^2+delM^2)^(-1/2)*2*delM*ddelM_dI) * (delM*tanh(qV)+sqrt(delM0^2+delM^2)) ...
    //      ) / ...
    //      (2*(Kirr-alpha*sqrt(delM0^2+delM^2)))^2; % vector
    dMirrp_dI[i] = (1.0/(2.0*(Kirr - Alpha*sq_delM02delM2))) *
                   (tanh_qV + delM/sq_delM02delM2 +
                     (2.0*Alpha*delM*(delM*tanh_qV +
                       sq_delM02delM2)/(2.0*(Kirr-Alpha*sq_delM02delM2)*sq_delM02delM2))) * ddelM_dI[i];

    // dP_dI = ( (c*dManp_dI+(1-c)*dMirrp_dI) * (1+(lg/lt-alpha)*c*Manp+(lg/lt)*(1-c)*Mirrp) - ...
    //            ((lg/lt-alpha)*c*dManp_dI + (lg/lt)*(1-c)*dMirrp_dI) * (c*Manp+(1-c)*Mirrp) ...
    //      ) / ...
    //      (1+(lg/lt-alpha)*c*Manp+(lg/lt)*(1-c)*Mirrp)^2; % vector
    dP_dI[i] = (1.0/dP_Denom) * (C * dManp_dI[i] + (1.0-C) * dMirrp_dI[i]) -
          ( (C*Manp + (1.0-C)*Mirrp)/pow(dP_Denom,2.0) ) *
            ( (gap_path - Alpha)*C*dManp_dI[i] + gap_path*(1.0-C)*dMirrp_dI[i] );

   /* not a reliable way to get dP_dI. Need to check why.
   const int numDeriv = 2;
   Sacado::Fad::SFad<double,4> MagVarFad( numDeriv, 0, latestMag );
   Sacado::Fad::SFad<double,4> branchCurrentSumVarFad( numDeriv, 1, solVector[instanceData[i]->li_Branch] * inductanceVals[ i ] );
   Sacado::Fad::SFad<double,4> V1PosVarFad( numDeriv, 2, V1Pos );
   Sacado::Fad::SFad<double,4> V1NegVarFad( numDeriv, 3, V1Neg );

   Sacado::Fad::SFad<double,4> resultFad;
   resultFad = Pcalc( MagVarFad, branchCurrentSumVarFad, V1PosVarFad, V1NegVarFad );
   dP_dI[i] = resultFad.dx(1)/inductanceVals[i];


   P = resultFad.val();
   dP_dM = resultFad.dx(0);
   dP_dBranchCurrentSum = resultFad.dx(1);
   dP_dV1Pos = resultFad.dx(2);
   dP_dV1Neg = resultFad.dx(2);
   */

  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateInductanceMatrix()
// Purpose       : A matrix of inductances is used often enough that it
//                 calculated and stored as a member variable here
//                 If and inductance ever changes say from updating
//                 the temperature or a parameter udpate, then this
//                 routine must be called again.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
void Instance::updateInductanceMatrix()
{
  vector< InductorInstanceData* >::iterator
    currentInductor = instanceData.begin();
  vector< InductorInstanceData* >::iterator
    endInductor = instanceData.end();

  // collec the inductances
  int i=0;
  while( currentInductor != endInductor )
  {
    inductanceVals[ i ] = (*currentInductor)->L;
    i++;
    currentInductor++;
  }

  double Area = model_.Area;
  double Path = model_.Path;

  // compute the inductance matrix
  for( i=0; i<numInductors; ++i)
  {
    for( int j=0; j<numInductors; ++j)
    {
      // 4.0e-7 * M_PI is a magnetic constant, the permeability of free space [Henries/m]
      LO[i][j] = mutualCup * 4.0e-7 * M_PI * (Area / Path) * inductanceVals[i] * inductanceVals[j];
    }
  }

}

//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState()
{
  bool bsuccess = true;
#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << "Instance::updatePrimaryState---------------" << endl
         << "\tname = " << getName() << endl;
  }
#endif
  // udate dependent parameters
  updateIntermediateVars ();

#if 0
  // don't need to do this as we're not using the state vector
  N_LAS_Vector & solVector = *(extData.nextSolVectorPtr);
  N_LAS_Vector & staVector = *(extData.nextStaVectorPtr);

  vector< InductorInstanceData* >::iterator currentInductor = instanceData.begin();
  vector< InductorInstanceData* >::iterator endInductor = instanceData.end();
  int i = 0;
  while( currentInductor != endInductor )
  {
    double current = solVector[ ( (*currentInductor)->li_Branch) ];
    if( (getSolverState().dcopFlag) && ((*currentInductor)->ICGiven) )
    {
      current = (*currentInductor)->IC;
    }
    // place this value for the charge in the state vector.
    staVector[((*currentInductor)->li_currentState)] = current;
    currentInductor++;
    i++;
  }
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
bool Instance::updateSecondaryState ()
{
  bool bsuccess = true;
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::acceptStep
// Purpose       : This function updates the value of MagVar
//
// Scope         : public
// Creator       : Rich Schiek, SNL, Electrical Systems Modeling
// Creation Date : 01/25/2011
//-----------------------------------------------------------------------------
void Instance::acceptStep()
{
  if (!getSolverState().dcopFlag)
  {
   if(includeDeltaM)
   {
     MagVar += MagVarUpdate*model_.Ms;
   }
   else
   {
     MagVar += MagVarUpdate;
   }
   lastMagUpdate = MagVarUpdate;
   PPreviousStep = P;
   if( fabs(MagVar) > 2*model_.Ms )
   {
     MagVar = 0.0;
   }

   if( useRKIntegration )
   {
     // fill in history for RK integration of dM/dH
     for(int i=0; i<2; i++)
     {
       branchCurrentSumHistory[i] = branchCurrentSumHistory[i+1];
       PFunctionHistory[i] = PFunctionHistory[i+1];
     }
     branchCurrentSumHistory[2] = branchCurrentSum-oldBranchCurrentSum;
     PFunctionHistory[2] = PPreviousStep;
   }
   oldBranchCurrentSum = branchCurrentSum;
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDeviceMask
//
// Purpose       : Loads the zero elements of the device mask
//
// Special Notes : elements of the error vector associated with zero
//                 elements of the mask will not be included in weighted
//                 norms by the time integrator.
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 01/19/07
//-----------------------------------------------------------------------------
bool Instance::loadDeviceMask ()
{
  bool returnVal=false;
#if 0
  // ifndef Xyce_NO_MUTIND_MASK
  N_LAS_Vector * maskVectorPtr = extData.deviceMaskVectorPtr;

  (*maskVectorPtr)[li_MagVar] = 0.0;
  returnVal = true;
#endif
  return (returnVal);
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
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
bool Instance::loadDAEQVector ()
{
  bool bsuccess = true;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << "Instance::loadDAEQVector------------------------" << endl
         << "\tname = " << getName() << endl;
  }
#endif

  N_LAS_Vector * daeQVecPtr = extData.daeQVectorPtr;
  N_LAS_Vector & solVector = *(extData.nextSolVectorPtr);

  // update LOI -- the following product
  // I = column vector of currents
  // L = row vector of inductances
  // LO = matrix = mutualCup * sqrt( L' * L )
  // LOI = column vector = mutualCup * sqrt( L' * L ) * I
  // LOI[1] = mutualCup * sqrt(L[1]*L[1])*I[1]) +
  //          mutualCup * sqrt(L[1]*L[2])*I[2]) + ...
  //          mutualCup * sqrt(L[1]*L[n])*I[n])

  vector< InductorInstanceData* >::iterator currentInductor = instanceData.begin();
  vector< InductorInstanceData* >::iterator endInductor = instanceData.end();
  int i = 0;
  while( currentInductor != endInductor )
  {
    if( (getSolverState().dcopFlag) && (*currentInductor)->ICGiven == true )
    {
      inductorCurrents[ i ] = (*currentInductor)->IC;
    }
    else
    {
      inductorCurrents[ i ] = solVector[ (*currentInductor)->li_Branch ];
    }
    i++;
    currentInductor++;
  }

  for( i = 0; i < numInductors; ++i )
  {
    LOI[ i ] = 0;
    for( int j = 0; j < numInductors; ++j )
    {
      LOI[i] += LO[i][j] * inductorCurrents[j];
    }
  }

  // loop over each inductor and load it's Q vector components
  // and each inductor's contribution to the R equ.
  currentInductor = instanceData.begin();
  endInductor = instanceData.end();
  i = 0;
  while( currentInductor != endInductor )
  {

    (*daeQVecPtr)[((*currentInductor)->li_Branch)] += LOI[ i ];
    i++;
    currentInductor++;
  }

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 instance.
//
// Special Notes : See the special notes for loadDAEFVector.
//
//                 Same as loadRHS, but without the capacitor
//                 currents.
//
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  bool bsuccess=true;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << "Instance::loadDAEFVector------------------------" << endl
         << "\tname = " << getName() << endl;
  }
#endif

  N_LAS_Vector * daeFVecPtr = extData.daeFVectorPtr;
  N_LAS_Vector & solVector = *(extData.nextSolVectorPtr);

  double Gap = model_.Gap;
  double Path = model_.Path;

  // for the M equation
/*

  NoMag (*daeFVecPtr)[li_MagVar] += mEquFval;
*/

  // used in scaling the branch equation;
  double mid = 1.0 + (1.0 - (Gap / Path))*P;

  // loop over each inductor and load it's F vector components
  vector< InductorInstanceData* >::iterator currentInductor = instanceData.begin();
  vector< InductorInstanceData* >::iterator endInductor = instanceData.end();
  int i=0;
  while( currentInductor != endInductor )
  {
    double current   = solVector[(*currentInductor)->li_Branch];
    double vNodePos  = solVector[(*currentInductor)->li_Pos];
    double vNodeNeg  = solVector[(*currentInductor)->li_Neg];


    (*daeFVecPtr)[((*currentInductor)->li_Pos)]    +=  current;

    (*daeFVecPtr)[((*currentInductor)->li_Neg)]    += -current;

    (*daeFVecPtr)[((*currentInductor)->li_Branch)] += -((vNodePos - vNodeNeg)/mid);

    currentInductor++;
    i++;
  }

  if(includeDeltaM)
  {
    // the deltaHapp equation
    //(*daeFVecPtr)[li_deltaHappVar] += solVector[li_deltaHappVar];
    //(*daeFVecPtr)[li_deltaHappVar] -= HappVarUpdate;

    // the deltaM equation
    (*daeFVecPtr)[li_deltaMagVar] += solVector[li_deltaMagVar];
    (*daeFVecPtr)[li_deltaMagVar] -= MagVarUpdate;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx
//
// Purpose       : Loads the Q-vector contributions for a single
//                 instance.
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  bool bsuccess = true;
  int i;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << "Instance::loadDAEdQdx-----------------------" << endl
         << "\tname = " << getName() << endl;
  }
#endif
  N_LAS_Matrix * dQdxMatPtr = extData.dQdxMatrixPtr;

  // loop over each inductor and load it's Q vector components
  vector< InductorInstanceData* >::iterator currentInductor = instanceData.begin();
  vector< InductorInstanceData* >::iterator endInductor = instanceData.end();
  i = 0;
  while( currentInductor != endInductor )
  {
    for( int j=0; j<numInductors; ++j )
    {
      (*dQdxMatPtr)[((*currentInductor)->li_Branch)]
                   [(*currentInductor)->inductorCurrentOffsets[j]] += LO[i][j];
    }
    i++;
    currentInductor++;
  }

  return bsuccess;
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
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  bool bsuccess = true;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << "Instance::loadDAEdFdx----------------------" << endl
         << "\tname = " << getName() << endl;
  }
#endif

  N_LAS_Vector & solVector = *(extData.nextSolVectorPtr);
  N_LAS_Vector & lastSolVector = *(extData.lastSolVectorPtr);
  N_LAS_Matrix * dFdxMatPtr = extData.dFdxMatrixPtr;

  // pull these parameters up from the model class to make it easier
  // to view the equations.
  const double Gap = model_.Gap;
  const double Path = model_.Path;

  // terms for the M equation
  /* NoMag
  (*dFdxMatPtr)[ li_MagVar ][ mEquMOffset ]    += 1.0;
  for( int i=0; i<numInductors; i++ )
  {
    (*dFdxMatPtr)[ li_MagVar ][ mEquInductorOffsets[i] ] = inductorInductances[i] * deltaBranchCurrentSum * dP_dBranchCurrentSum;
  }
  */

  // loop over each inductor and load it's dFdx components
#ifdef MS_FACTORING2
  double mid = 1.0 + (1.0 - (Gap/Path))*P*model_.Ms;
#else
  double mid = 1.0 + (1.0 - (Gap/Path))*P;
#endif

  vector< InductorInstanceData* >::iterator currentInductor = instanceData.begin();
  vector< InductorInstanceData* >::iterator endInductor = instanceData.end();
  while( currentInductor != endInductor )
  {
    // do the normal work for an inductor
    (*dFdxMatPtr)[((*currentInductor)->li_Pos)]   [((*currentInductor)->APosEquBraVarOffset)]  +=  1.0;
    (*dFdxMatPtr)[((*currentInductor)->li_Neg)]   [((*currentInductor)->ANegEquBraVarOffset)]  += -1.0;
    (*dFdxMatPtr)[((*currentInductor)->li_Branch)][((*currentInductor)->ABraEquPosNodeOffset)] += -1.0/mid;
    (*dFdxMatPtr)[((*currentInductor)->li_Branch)][((*currentInductor)->ABraEquNegNodeOffset)] +=  1.0/mid;

    double delV = solVector[(*currentInductor)->li_Pos] - solVector[(*currentInductor)->li_Neg];

    for( int j = 0; j<numInductors; ++j )
    {

      (*dFdxMatPtr)[((*currentInductor)->li_Branch)][(*currentInductor)->inductorCurrentOffsets[j]] +=
        delV * (1.0 - (Gap/Path)) * dP_dI[j]/(mid*mid);

    /* getting the dP_dI this way doesn't work very well.
     (*dFdxMatPtr)[((*currentInductor)->li_Branch)][(*currentInductor)->inductorCurrentOffsets[j]] +=
        delV * (1.0 - (Gap/Path)) * dP_dBranchCurrentSum * inductorInductances[j] / (mid*mid);
    */
    }

    /*NoMag
    (*dFdxMatPtr)[(*currentInductor)->li_Branch][(*currentInductor)->magOffset]  += delV * (1.0 - (Gap/Path)) * dP_dM/(mid*mid);

    (*dFdxMatPtr)[(*currentInductor)->li_Branch][(*currentInductor)->vPosOffset] += delV * (1.0 - (Gap/Path)) * dP_dV1Pos/(mid*mid);

    (*dFdxMatPtr)[(*currentInductor)->li_Branch][(*currentInductor)->vNegOffset] += delV * (1.0 - (Gap/Path)) * dP_dV1Neg/(mid*mid);
    */
    currentInductor++;
  }

  if(includeDeltaM)
  {
    // the deltaHapp equation
    //(*dFdxMatPtr)[li_deltaHappVar][deltaHappEquDeltaHappOffset] = 1.0;

    // the deltaM equation
    (*dFdxMatPtr)[li_deltaMagVar][deltaMEquDeltaMOffset] = 1.0;
    //(*dFdxMatPtr)[li_deltaMagVar][deltaMEquDeltaHappOffset] = -P/model_.Ms;

    for( int i=0; i<numInductors; i++ )
    {
      //(*dFdxMatPtr)[li_deltaHappVar][deltaHappEquInductorOffsets[i]] = -inductanceVals[i] / Path;
      (*dFdxMatPtr)[li_deltaMagVar][deltaMEquInductorOffsets[i]] =
        -((inductanceVals[i]*(solVector[ instanceData[i]->li_Branch ] - lastSolVector[ instanceData[i]->li_Branch ] ) * dP_dI[i])
        + P*inductanceVals[i])/(model_.Path*model_.Ms );
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::outputPlotFiles
// Purpose       : If requested by the use in the model statement,
//                 this routine outputs values of the internal
//                 state variables M, H and R to a file
//                 named "Inductor_name.dat".  File is opened
//                 and closed in the contructor and destructor.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
bool Instance::outputPlotFiles ()
{
  bool bsuccess = true;
  if( outputStateVarsFlag && (*outputFileStreamPtr) )
  {
    N_LAS_Vector & solVector = *(extData.nextSolVectorPtr);

#ifdef MS_FACTORING2
  double latestMag = MagVar*model_.Ms;
#else
  double latestMag = MagVar;
#endif

  if( includeDeltaM )
  {
  (*outputFileStreamPtr)
      << getSolverState().currTime << "  "
      << latestMag
      << endl;
  }
  else
  {
    (*outputFileStreamPtr)
      << getSolverState().currTime << "  "
      << latestMag
      << endl;
  }

  }
  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : Instance::setIC
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
bool Instance::setIC ()
{
  int i_bra_sol;
  int i_f_state;

  bool bsuccess = true;
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::varTypes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
void Instance::varTypes( vector<char> & varTypeVec )
{
  varTypeVec.resize(numInductors);
  for(int i=0; i<numInductors; i++)
  {
    varTypeVec[i] = 'I';
  }
}



//-----------------------------------------------------------------------------
// Function      : Instance::Pcalc
// Purpose       :
  // this is a templated function for a complicated term P(M,I_1... I_n) that relates
  // the magnetic saturation of the mutual indcutor to the individual currents
  // through the inductors.  We'll need dP_dM and this tempated function automates
  // that calculation via Sacado
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 10/13/2011
//-----------------------------------------------------------------------------
template <typename ScalarT>
ScalarT Instance::Pcalc(
    const ScalarT & Mag, const ScalarT & CurrentSum,
    const ScalarT & Vpos, const ScalarT & Vneg)
   // independent variable M
   // independent variable Sum(n_i * I_i) windings * current through each inductor
   // independent variable Vpos, Vneg -- voltage drop over first inductor
{
     // some parameters in the model class that we will use often
     const double A      = model_.A;
     const double Alpha  = model_.Alpha;
     const double Area   = model_.Area;
     const double BetaH  = model_.BetaH;
     const double BetaM  = model_.BetaM;
     const double C      = model_.C;
     const double DeltaV = model_.DeltaV;
     const double Gap    = model_.Gap;
     const double Ms     = model_.Ms;
     const double Kirr   = model_.Kirr;
     const double Path   = model_.Path;
     const double Vinf   = model_.Vinf;

     ScalarT qV = (DeltaV / Vinf) * (Vpos - Vneg);

     ScalarT tanh_qV = 0.0;
     if (fabs(qV) < CONSTTANH_THRESH)
     {
       ScalarT tanh_qV = tanh(qV);
     }
     else if (qV < 0.0)
     {
       tanh_qV = -1.0;
     }
     else
     {
       tanh_qV = 1.0;
     }

     ScalarT Happ = CurrentSum / Path;

#ifdef MS_FACTORING2
     ScalarT H = Happ - (Gap / Path) * Mag * Ms;
     ScalarT He = H + Alpha * Mag * Ms;
#else
     ScalarT H = Happ - (Gap / Path) * Mag;
     ScalarT He = H + Alpha * Mag;
#endif

     ScalarT Heo = BetaH*A;

     // terms that come up frequently
     ScalarT gap_path = Gap / Path;
     ScalarT He2 = He*He;
     ScalarT Heo2 = Heo*Heo;
     ScalarT sq_Heo2He2 = sqrt(Heo2 + He2);

     ScalarT delM0 = model_.BetaM * Ms;
     ScalarT Man = Ms * He / ( A + sq_Heo2He2 );
#ifdef MS_FACTORING2
     ScalarT delM = Man - Mag*Ms;
#else
     ScalarT delM = Man - Mag;
#endif

     ScalarT delM2 = delM*delM;
     ScalarT delM02 = delM0*delM0;
     ScalarT sq_delM02delM2 = sqrt( delM02 + delM2 );

     //Mirrp = (delM * tanh_qV + sq_delM02delM2 ) / (2*( Kirr- Alpha * Ms * sq_delM02delM2));
     //Manp =  (A + Heo2/sq_Heo2He2) / pow(A + sq_Heo2He2, 2.0);
     //P = ( C * Manp + (1 - C) * Mirrp) / (1 + (gap_path - Alpha) * C * Ms * Manp + gap_path * (1-C) * Ms * Mirrp);
#ifdef MS_FACTORING2
     ScalarT Mirrp = (delM * tanh_qV + sq_delM02delM2 ) / (2*( Kirr- Alpha * sq_delM02delM2));
     ScalarT Manp =  Ms * (A + Heo2/sq_Heo2He2) / pow(A + sq_Heo2He2, 2.0);
     ScalarT Pval = ( C * Manp + (1 - C) * Mirrp) / ((1 + (gap_path - Alpha) * C * Manp + gap_path * (1-C) * Mirrp)*Ms);
#else
     ScalarT Mirrp = (delM * tanh_qV + sq_delM02delM2 ) / (2*( Kirr- Alpha * sq_delM02delM2));
     ScalarT Manp =  Ms*(A + Heo2/sq_Heo2He2) / pow(A + sq_Heo2He2, 2.0);
     ScalarT Pval = ( C * Manp + (1 - C) * Mirrp) / (1 + (gap_path - Alpha) * C * Manp + gap_path * (1-C) * Mirrp);
#endif

// matlab code
/*
    V = Vpos-Vneg;

    qV = deltav*V/vinf;
%    qV = real(V); % scalar (analytic extension for complex step derivative check)

    delM0 = betaM*Ms; % scalar, modeling parameter
    Happ = (1/lt)*sum(I.*N'); % scalar
    H = Happ-(lg/lt)*M; % scalar
    He = H + alpha*M; % scalar
%    qHe = deltaHe*He/Heinf;
    Heo = betaH*A; % scalar, modeling parameter
%    Man = Ms*He/(A+He*tanh(qHe)); % scalar
    Man = Ms*He/(A+sqrt(Heo^2+He^2)); % scalar
    delM = Man - M; % scalar
%    qdelM = deltadelM*delM/delMinf;
%    Mirrp = (delM*tanh(qV)+delM*tanh(qdelM))/(2*(Kirr-alpha*delM*tanh(qdelM))); % scalar
    Mirrp = (delM*tanh(qV)+sqrt(delM0^2+delM^2))/(2*(Kirr-alpha*sqrt(delM0^2+delM^2))); % scalar
%    Manp = (Ms*(A - He^2*(1-tanh(qHe)^2)*deltaHe/Heinf))/((A+He*tanh(qHe))^2); % scalar
    Manp = (Ms*(A + (Heo^2)/sqrt(Heo^2+He^2)))/((A+sqrt(Heo^2+He^2))^2); % scalar
    P = (c*Manp+(1-c)*Mirrp)/(1+(lg/lt-alpha)*c*Manp+(lg/lt)*(1-c)*Mirrp); % scalar
*/

     return Pval;
}  // end of Pcalc() function


//-----------------------------------------------------------------------------
// Function      : Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
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
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
Model::Model (const ModelBlock & MB,
                                                SolverState & ss1,
                                                DeviceOptions & do1)
  : DeviceModel(MB,ss1,do1),
    A(0.0),
    Alpha(0.0),
    Area(0.0),
    BetaH(0.0),
    BetaM(0.0),
    C(0.0),
    DeltaV(0.0),
    Gap(0.0),
    Kirr(0.0),
    Ms(0.0),
    Path(0.0),
    Vinf(0.0),
    tempCoeff1(0.0),
    tempCoeff2(0.0),
    outputStateVars(0),
    useRKIntegration(0),
    includeDeltaM(0),
    tnom(do1.tnom)
{
  setLevel(2);


  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to .model line and constant defaults from metadata:
  setModParams (MB.params);

  // scale gap, path and area from cm and cm^2 to m and m^2
  Gap *= 1.0e-2;
  Path *= 1.0e-2;
  Area *= 1.0e-4;

  // Set any non-constant parameter defaults:

  if (!given("TNOM"))
  {
    tnom = getDeviceOptions().tnom;
  }

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
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
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

// additional Declarations

//-----------------------------------------------------------------------------
// Function      : Model::printOutInstances
// Purpose       : debugging tool.
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
std::ostream &Model::printOutInstances(std::ostream &os) const
{
  vector<Instance*>::const_iterator iter;
  vector<Instance*>::const_iterator first = instanceContainer.begin();
  vector<Instance*>::const_iterator last  = instanceContainer.end();

  int i, isize;
  isize = instanceContainer.size();

  os << endl;
  os << "Number of MutIndNonLin instances: " << isize << endl;
  os << "    name=\t\tmodelName\tParameters" << endl;
  for (i=0, iter=first; iter!=last; ++iter, ++i)
  {
    os << "  " << i << ": " << (*iter)->getName() << "\t";
    os << (*iter)->getModelName();
    os << endl;
  }

  os << endl;

  return os;
}

} // namespace MutIndNonLin2
} // namespace Device
} // namespace Xyce

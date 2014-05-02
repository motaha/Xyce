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
// Filename       : $RCSfile: N_DEV_MutIndNonLin.C,v $
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
// Revision Number: $Revision: 1.150.2.3 $
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

#include <fstream>
#include <algorithm>
#include <vector>
#include <set>

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_MutIndNonLin.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_Message.h>
#include <N_ERH_ErrorMgr.h>

//This contains important constants like permitivity of free space
#include <N_DEV_Const.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>

#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RefCountPtr;
using Teuchos::rcp;

namespace Xyce {
namespace Device {


namespace MutIndNonLin {


void Traits::loadInstanceParameters(ParametricData<MutIndNonLin::Instance> &p)
{
    // Set up configuration constants:
// Set up double precision variables:
    p.addPar ("COUP_VAL",   1.0, false, ParameterType::NO_DEP,
        &MutIndNonLin::Instance::mutualCup,
        &MutIndNonLin::Instance::mutualCupGiven,
        U_NONE, CAT_NONE, "Coupling coefficient");

    p.addPar ("NONLINEARCOUPLING", 0.0, false, ParameterType::NO_DEP,
        &MutIndNonLin::Instance::nonlinFlag,
        &MutIndNonLin::Instance::nonlinFlagGiven,
        U_NONE, CAT_NONE, "Nonlinear coupling flag");

    // Set up non-double precision variables:
    p.addPar ("COUPLEDMutIndNonLin", std::vector<std::string>(), false, ParameterType::NO_DEP,
            &MutIndNonLin::Instance::inductorNames, NULL, U_NONE, CAT_NONE, "" );
    p.addPar ("COUPLEDINDUCTANCE", std::vector<double>(), false, ParameterType::NO_DEP,
            &MutIndNonLin::Instance::inductorInductances, NULL, U_NONE, CAT_NONE, "");
    p.addPar ("NODE1", std::vector<std::string>(), false, ParameterType::NO_DEP,
            &MutIndNonLin::Instance::inductorsNode1, NULL, U_NONE, CAT_NONE, "");
    p.addPar ("NODE2", std::vector<std::string>(), false, ParameterType::NO_DEP,
            &MutIndNonLin::Instance::inductorsNode2, NULL, U_NONE, CAT_NONE, "");
    p.addPar ("COUPLING", std::vector<double>(), false, ParameterType::NO_DEP,
            &MutIndNonLin::Instance::couplingCoefficient, NULL, U_NONE, CAT_NONE, "Coupling coefficient");
    p.addPar ("COUPLEDINDUCTOR", std::vector<std::string>(), false, ParameterType::NO_DEP,
            &MutIndNonLin::Instance::couplingInductor, NULL, U_NONE, CAT_NONE, "");
}

void Traits::loadModelParameters(ParametricData<MutIndNonLin::Model> &p)
{
    p.addPar ("A", 1000.0, false, ParameterType::NO_DEP,
      &MutIndNonLin::Model::A,
      NULL, U_AMPMM1, CAT_MATERIAL, "Thermal energy parameter");

    p.addPar ("AREA", 0.1, false, ParameterType::NO_DEP,
      &MutIndNonLin::Model::Area,
      NULL, U_CM2, CAT_GEOMETRY, "Mean magnetic cross-sectional area");

    p.addPar ("ALPHA", 5.0e-5, false, ParameterType::NO_DEP,
      &MutIndNonLin::Model::Alpha,
      NULL, U_NONE, CAT_GEOMETRY, "Domain coupling parameter");

    p.addPar ("BETAH", 0.0001, false, ParameterType::NO_DEP,
      &MutIndNonLin::Model::BetaH,
      NULL, U_NONE, CAT_NONE, "Modeling constant");

    p.addPar ("BETAM", 3.125e-5, false, ParameterType::NO_DEP,
      &MutIndNonLin::Model::BetaM,
      NULL, U_NONE, CAT_NONE, "Modeling constant");

    p.addPar ("C", 0.2, false, ParameterType::NO_DEP,
      &MutIndNonLin::Model::C,
      NULL, U_NONE, CAT_MATERIAL, "Domain flexing parameter");

    //p.addPar ("DELVSCALING", 1.0e3, false, ParameterType::NO_DEP,
    p.addPar ("DELVSCALING", 1.0e-1, false, ParameterType::NO_DEP,
      &MutIndNonLin::Model::DeltaVScaling,
      NULL, U_VOLT, CAT_NONE, "Smoothing coefficient for voltage difference over first inductor");

    p.addPar ("CONSTDELVSCALING", true, false, ParameterType::NO_DEP,
      &MutIndNonLin::Model::UseConstantDeltaVScaling,
      NULL, U_VOLT, CAT_NONE, "Use constant scaling factor to smooth voltage difference over first inductor");

    p.addPar ("GAP", 0.0, false, ParameterType::NO_DEP,
      &MutIndNonLin::Model::Gap,
      NULL, U_CM, CAT_GEOMETRY, "Effective air gap");

    p.addPar ("K", 500.0, false, ParameterType::NO_DEP,
      &MutIndNonLin::Model::Kirr,
      NULL, U_AMPMM1, CAT_MATERIAL, "Domain anisotropy parameter");

    p.addPar ("KIRR", 500.0, false, ParameterType::NO_DEP,
      &MutIndNonLin::Model::Kirr,
      NULL, U_AMPMM1, CAT_MATERIAL, "Domain anisotropy parameter");

    p.addPar ("MS", 1.0e+6, false, ParameterType::NO_DEP,
      &MutIndNonLin::Model::Ms,
      NULL, U_AMPMM1, CAT_MATERIAL, "Saturation magnetization");

    p.addPar ("LEVEL", 0.0, false, ParameterType::NO_DEP,
      &MutIndNonLin::Model::LevelIgnored,
      NULL, U_NONE, CAT_NONE, "for pspice compatibility -- ignored");

    p.addPar ("PACK", 0.0, false, ParameterType::NO_DEP,
      &MutIndNonLin::Model::PackIgnored,
      NULL, U_NONE, CAT_NONE, "for pspice compatibility -- ignored");

    p.addPar ("PATH", 1.0, false, ParameterType::NO_DEP,
      &MutIndNonLin::Model::Path,
      NULL, U_CM, CAT_GEOMETRY, "Total mean magnetic path");

    p.addPar ("TNOM", 27.0, false, ParameterType::NO_DEP,
      &MutIndNonLin::Model::tnom,
      NULL, U_DEGC, CAT_MATERIAL, "Reference temperature");

    p.addPar ("TC1", 0.0, false, ParameterType::NO_DEP,
      &MutIndNonLin::Model::tempCoeff1,
      NULL, U_NONE, CAT_MATERIAL, "First order temperature coeff.");

    p.addPar ("TC2", 0.0, false, ParameterType::NO_DEP,
      &MutIndNonLin::Model::tempCoeff2,
      NULL, U_NONE, CAT_MATERIAL, "Second order temperature coeff.");

    p.addPar ("PZEROTOL", 0.1, false, ParameterType::NO_DEP,
      &MutIndNonLin::Model::pZeroTol,
      NULL, U_NONE, CAT_NONE, "Tolerance for nonlinear zero crossing");

    p.addPar ("MVARSCALING", 1.0, false, ParameterType::NO_DEP,
      &MutIndNonLin::Model::mVarScaling,
      NULL, U_NONE, CAT_NONE, "M-variable scaling.");

    p.addPar ("RVARSCALING", 1.0, false, ParameterType::NO_DEP,
      &MutIndNonLin::Model::rVarScaling,
      NULL, U_NONE, CAT_NONE, "R-variable scaling");

    p.addPar ("MEQNSCALING", 1.0, false, ParameterType::NO_DEP,
      &MutIndNonLin::Model::mEqScaling,
      NULL, U_NONE, CAT_NONE, "M-equation scaling");

    p.addPar ("REQNSCALING", 1.0, false, ParameterType::NO_DEP,
      &MutIndNonLin::Model::rEqScaling,
      NULL, U_NONE, CAT_NONE, "R-equation scaling");

    // Set up non-double precision variables:
    p.addPar ("OUTPUTSTATEVARS", 0, false, ParameterType::NO_DEP,
      &MutIndNonLin::Model::outputStateVars,
      NULL, U_NONE, CAT_NONE, "Flag to save state variables" );
    p.addPar ("FACTORMS", 0, false, ParameterType::NO_DEP,
      &MutIndNonLin::Model::factorMS,
      NULL, U_NONE, CAT_NONE, "Flag to save state variables" );
    p.addPar ("BHSIUNITS", 0, false, ParameterType::NO_DEP,
      &MutIndNonLin::Model::BHSiUnits,
      NULL, U_NONE, CAT_NONE, "Flag to report B and H in SI units" );
}


// Class Instance

//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock &         IB,
  Model &                       Iiter,
  const FactoryBlock &          factory_block)
  : DeviceInstance(IB, configuration.getInstanceParameters(), factory_block),
    model_(Iiter),
    temp(getDeviceOptions().temp.getImmutableValue<double>()),
    outputStateVarsFlag( false ),
    maxVoltageDrop(1.0e-10)
{
  scalingRHS = 1.0;
  numExtVars   = 2;
  numIntVars   = 3;

  numStateVars = 2;
  setNumStoreVars(3);

  tempGiven    = false;

  const int ibev = IB.numExtVars;
  const int ibiv = IB.numIntVars;

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

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

  // set up the device connectivity map
  // each simple inductor in this mutual inductor
  // is maked as a connection (given a common, non-zero
  // value in devConMap)
  devConMap.resize(2*numInductors);
  for(int i=0; i<numInductors; i++)
  {
    devConMap[i] = devConMap[i+1] = (i+1);
  }

  mEquInductorOffsets.resize( numInductors );
  rEquInductorOffsets.resize( numInductors );
  inductorCurrents.resize( numInductors );
  inductanceVals.resize( numInductors );
  LOI.resize( numInductors );
  LO.resize( numInductors );
  for( int i=0; i<numInductors; ++i)
  {
    LO[i].resize( numInductors );
  }

  // set up the device connectivity map
  // each simple inductor in this mutual inductor
  // is maked as a connection (given a common, non-zero
  // value in devConMap)
  devConMap.resize(2*numInductors);
  for(int i=0; i<numInductors; i++)
  {
    devConMap[i] = devConMap[i+1] = (i+1);
  }

  updateInductanceMatrix();

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:
  processParams ();

  // if the user has requested output of the state variables M, H and R
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

    outputFileStreamPtr = rcp( new std::ofstream() );
    outputFileStreamPtr->open( filename.c_str() );
    if( !(*outputFileStreamPtr) )
    {
      std::string msg("Instance constructor.\n");
      msg += "\tCould not open file for output of state variables. name =" + getName();
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
    }
    (*outputFileStreamPtr).setf(std::ios::scientific, std::ios::floatfield );
    (*outputFileStreamPtr).width(20);
    (*outputFileStreamPtr).precision(12);
  }

  // size some vectors needed in loadDAEdFdx
  dHe_dI.resize( numInductors );
  dManp_dI.resize( numInductors );
  ddelM_dI.resize( numInductors );
  dMirrp_dI.resize( numInductors );
  dP_dI.resize( numInductors );

  // update internal/external/state variable counts
  numExtVars = 2*numInductors;
  numIntVars = numInductors + 2;
  numStateVars = 5;  // extra state variables for M and dM/dt and R, H & B
  //numStateVars += 2*numInductors;  individual inductors no longer need state / store space

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
  //          V1  V2  V3  V4 ... V2N  I1  I2  ... IN  M  R
  //  kcl1                             1
  //  kcl2                            -1
  //  kcl3                                 1
  //  kcl4                                -1
  //  branch1 1   -1                 L/dt  c  ... c   x
  //  branch2 x    x   1  -1          c  L/dt ... c   x
  //  M equ   x    x                  x    x  ... x   x  x
  //  R equ                           x    x  ... x      x
  //
  //  where "c" is an induced current change and "x" are
  //  values which must be computed.

  jacStamp.resize( 3 * numInductors + 2);

  for( int i=0; i< numInductors; ++i )
  {
    //
    // allocate space
    //
    // kcl V+ node
    jacStamp[2*i].resize(1);
    // kcl V- node
    jacStamp[2*i+1].resize(1);
    // branch node -- every branch needs to access the first
    // inductor's V+ and V-, so they all contribute there
    if( i == 0 )
    {
      jacStamp[2*numInductors].resize(numInductors + 3);
    }
    else
    {
      jacStamp[2*numInductors + i].resize(numInductors + 5);
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
      jacStamp[2*numInductors][numInductors+2] = 3*numInductors;
    }
    else
    {
      jacStamp[2*numInductors + i][0] = 0;
      jacStamp[2*numInductors + i][1] = 1;
      jacStamp[2*numInductors + i][2] = 2*i;
      jacStamp[2*numInductors + i][3] = 2*i + 1;
      for( int j=0; j<numInductors; ++j )
      {
        jacStamp[2*numInductors + i][j+4] = 2*numInductors + j;
      }
      jacStamp[2*numInductors + i][numInductors+4] = 3*numInductors;
    }
  }

  // now the M equation
  jacStamp[ 3*numInductors    ].resize(numInductors + 4);
  // and the R equations
  jacStamp[ 3*numInductors + 1].resize(numInductors + 1);

  // M offsets for V+ and V- on first inductor
  jacStamp[ 3*numInductors    ][0] = 0;
  jacStamp[ 3*numInductors    ][1] = 1;

  // M and R offsets to each inductor's branch equ.
  for(int i=0; i<numInductors; ++i)
  {
    jacStamp[ 3*numInductors     ][i+2]=2*numInductors+i;
    jacStamp[ 3*numInductors + 1 ][i]=2*numInductors+i;
  }

  // M offsets to M and R
  jacStamp[ 3*numInductors    ][numInductors + 2] = 3*numInductors;
  jacStamp[ 3*numInductors    ][numInductors + 3] = 3*numInductors + 1;

  // R offsets to R
  jacStamp[ 3*numInductors + 1][numInductors ]    = 3*numInductors + 1;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    Xyce::dout() << "Instance::Instance----------" << std::endl;
    Xyce::dout() << "numExtVars = " << numExtVars << ", " << ibev << std::endl
      << "numIntVars = " << numIntVars << ", " << ibiv << std::endl
      << "numStateVars = " << numStateVars << std::endl
      << "numInductors = " << numInductors << std::endl
      << "jacStamp = " << std::endl;
    for( int i = 0; i<jacStamp.size(); ++i )
    {
      Xyce::dout() << "jacStamp[ " << i << " ] = { ";
      for( int j=0; j<jacStamp[i].size(); ++j)
      {
        Xyce::dout() << jacStamp[i][j];
        if( j != ( jacStamp[i].size() -1 ) )
        {
          Xyce::dout() << ", ";
        }
      }
      Xyce::dout() << " }" << std::endl;
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
      std::string msg("Instance destructor.\n");
      msg+= "\tCould not close file for output of state variables. name =" + getName();
         N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
    }
    // since outputFileStreamPtr is a ref counted pointer
    // we don't need to delete it.
  }

  std::vector< InductorInstanceData* >::iterator currentInductor = instanceData.begin();
  std::vector< InductorInstanceData* >::iterator endInductor = instanceData.end();
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
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
bool Instance::processParams()
{
  // now set the temperature related stuff.
  if (tempGiven)
  {
    updateTemperature(temp);
  }

  return true;
}


//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
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
  // get the current values of the inductances and currentOffsets
  std::vector< InductorInstanceData* >::iterator currentInductor = instanceData.begin();
  std::vector< InductorInstanceData* >::iterator endInductor = instanceData.end();
  int i = 0;
  int j = 0;
  while( currentInductor != endInductor )
  {
    (*currentInductor)->li_Pos = extLIDVec[ i++ ];
    (*currentInductor)->li_Neg = extLIDVec[ i++ ];
    (*currentInductor)->li_Branch = intLIDVec[ j++ ];
    currentInductor++;
  }

  // now get the M and R local id's
  li_MagVar = intLIDVec[ j++ ];
  li_RVar   = intLIDVec[ j++ ];

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    Xyce::dout() << "Instance::registerLIDs------------------------" << std::endl;
    currentInductor = instanceData.begin();
    i=0;
    while( currentInductor != endInductor )
    {
      Xyce::dout() << "Inductor [ " << i++ << " ] "
           << "   li_Pos = " << (*currentInductor)->li_Pos
           << "   li_Neg = " << (*currentInductor)->li_Neg
           << "   li_Branch = " << (*currentInductor)->li_Branch << std::endl;
      currentInductor++;
    }
    Xyce::dout() << " li_MagVar = " << li_MagVar << std::endl
         << " li_RVar = " << li_RVar << std::endl;
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
std::map<int,std::string> & Instance::getIntNameMap()
{
  // set up the internal name map, if it hasn't been already.
  if (intNameMap.empty ())
  {
    std::string baseString(getName() + "_");
    std::string tempString;
    std::vector< InductorInstanceData* >::iterator
      currentInductor = instanceData.begin();
    std::vector< InductorInstanceData* >::iterator
      endInductor = instanceData.end();
    int i = 0;
    int j = 0;
    while( currentInductor != endInductor )
    {
      tempString = baseString + (*currentInductor)->name +"_branch";
      spiceInternalName (tempString);
      intNameMap[ (*currentInductor)->li_Branch ] = tempString;
      currentInductor++;
    }

    tempString = baseString + "m";
    intNameMap[ li_MagVar ] = tempString;
    tempString = baseString + "r";
    intNameMap[ li_RVar ] = tempString;
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
void Instance::registerStateLIDs( const std::vector<int> & staLIDVecRef )
{
  AssertLIDs(staLIDVecRef.size() == numStateVars);

  // copy over the global ID lists.
  staLIDVec = staLIDVecRef;
  int i = 0;

  li_MagVarState = staLIDVec[i++];
  li_MagVarDerivState = staLIDVec[i++];

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    Xyce::dout() << "Instance::registerStateLIDs-------------------" << std::endl;

    Xyce::dout() << "li_MagVarState = " << li_MagVarState << std::endl
      << "li_MagVarDerivState = " << li_MagVarDerivState << std::endl
      ;
  }
#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStoreLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL
// Creation Date : 8/17/2012
//-----------------------------------------------------------------------------
void Instance::registerStoreLIDs(const std::vector<int> & stoLIDVecRef )
{
  AssertLIDs(stoLIDVecRef.size() == getNumStoreVars());

  // copy over the global ID lists.
  stoLIDVec = stoLIDVecRef;

  li_RVarStore = stoLIDVec[0];
  li_BVarStore = stoLIDVec[1];
  li_HVarStore = stoLIDVec[2];
}

//-----------------------------------------------------------------------------
// Function      : Instance::getStateNameMap
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Electrical Systems Modeling
// Creation Date : 08/01/2012
//-----------------------------------------------------------------------------
std::map<int,std::string> & Instance::getStateNameMap()
{
  // set up the internal name map, if it hasn't been already.
  if (stateNameMap.empty ())
  {
    std::string baseString(getName() + "_");
    std::string tempString;
    tempString = baseString + "m";
    stateNameMap[ li_MagVarState ] = tempString;
    tempString = baseString + "dmdt";
    stateNameMap[ li_MagVarDerivState ] = tempString;
  }

  return stateNameMap;
}


//-----------------------------------------------------------------------------
// Function      : Instance::getStoreNameMap
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Electrical Systems Modeling
// Creation Date : 08/01/2012
//-----------------------------------------------------------------------------
std::map<int,std::string> & Instance::getStoreNameMap()
{

  // set up the internal name map, if it hasn't been already.
  if (storeNameMap.empty ())
  {
    std::string baseString(getName() + "_");
    std::string tempString;
    tempString = baseString + "r";
    storeNameMap[ li_RVarStore ] = tempString;
    tempString = baseString + "b";
    storeNameMap[ li_BVarStore ] = tempString;
    tempString = baseString + "h";
    storeNameMap[ li_HVarStore ] = tempString;
  }
  return storeNameMap;
}

//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
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
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );
#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    Xyce::dout() << "Instance::registerJacLIDs ----------------------------" << std::endl;

    Xyce::dout() << "jacLIDVec = " << std::endl;
    for( int i = 0; i<jacStamp.size(); ++i )
    {
      Xyce::dout() << "jacLIDVec[ " << i << " ] = { ";
      for( int j=0; j<jacLIDVec[i].size(); ++j)
      {
        Xyce::dout() << jacLIDVec[i][j];
        if( j != ( jacLIDVec[i].size() -1 ) )
        {
          Xyce::dout() << ", ";
        }
      }
      Xyce::dout() << " }" << std::endl;
    }
  }
#endif

  std::vector< InductorInstanceData* >::iterator currentInductor = instanceData.begin();
  std::vector< InductorInstanceData* >::iterator endInductor = instanceData.end();
  // int numInductors = instanceData.size();  // don't need this as it's defined at class level
  int i = 0;
  while( currentInductor != endInductor )
  {
    (*currentInductor)->APosEquBraVarOffset  = jacLIDVec[ 2*i     ][ 0 ];
    (*currentInductor)->ANegEquBraVarOffset  = jacLIDVec[ 2*i + 1 ][ 0 ];
    (*currentInductor)->vPosOffset = jacLIDVec[ 2*numInductors + i ][ 0 ];
    (*currentInductor)->vNegOffset = jacLIDVec[ 2*numInductors + i ][ 1 ];
    int extraOffset = 2;
    if( i == 0)
    {
      extraOffset = 0;
    }
    (*currentInductor)->ABraEquPosNodeOffset = jacLIDVec[ 2*numInductors + i ][ 0 + extraOffset ];
    (*currentInductor)->ABraEquNegNodeOffset = jacLIDVec[ 2*numInductors + i ][ 1 + extraOffset ];
    for( int j=0; j<numInductors; ++j )
    {
      if( i == j )
      {
        (*currentInductor)->ABraEquBraVarOffset  = jacLIDVec[ 2*numInductors + i ][ j + 2 + extraOffset ];
      }
      (*currentInductor)->inductorCurrentOffsets[ j ] = jacLIDVec[ 2*numInductors + i ][ j + 2 + extraOffset ];
    }
    (*currentInductor)->magOffset = jacLIDVec[ 2*numInductors + i ][ numInductors + 2 + extraOffset ];
    currentInductor++;
    i++;
  }

  // now get the M equation offsets
  mEquVPosOffset = jacLIDVec[ 3*numInductors ][0];
  mEquVNegOffset = jacLIDVec[ 3*numInductors ][1];
  for( i=0; i<numInductors; ++i )
  {
    mEquInductorOffsets[i] = jacLIDVec[ 3*numInductors ][ i + 2];
  }
  mEquMOffset = jacLIDVec[ 3*numInductors ][ numInductors + 2 ];
  mEquROffset = jacLIDVec[ 3*numInductors ][ numInductors + 3 ];

  // now get the R equation offsets
  for( i=0; i<numInductors; ++i )
  {
    rEquInductorOffsets[i] = jacLIDVec[ 3*numInductors + 1 ][ i ];
  }
  rEquROffset = jacLIDVec[ 3*numInductors + 1 ][ numInductors ];

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    currentInductor = instanceData.begin();
    i=0;
    while( currentInductor != endInductor )
    {
      Xyce::dout() << "Inductor [ " << i << " ] " << (*currentInductor)->name << std::endl
           << "   APosEquBraVarOffset = " << (*currentInductor)->APosEquBraVarOffset << std::endl
           << "   ANegEquBraVarOffset = " << (*currentInductor)->ANegEquBraVarOffset << std::endl
           << "   vPosOffset = " << (*currentInductor)->vPosOffset << std::endl
           << "   vNegOffset = " << (*currentInductor)->vNegOffset << std::endl
           << "   ABraEquPosNodeOffset = " << (*currentInductor)->ABraEquPosNodeOffset << std::endl
           << "   ABraEquNegNodeOffset = " << (*currentInductor)->ABraEquNegNodeOffset << std::endl
           << "   ABraEquBraVarOffset = " << (*currentInductor)->ABraEquBraVarOffset << std::endl
           << "   magOffset = " << (*currentInductor)->magOffset << std::endl;
      Xyce::dout() << "\tInductor branch offsets = { ";
      for( int j=0; j<numInductors ; ++j )
      {
        Xyce::dout() << (*currentInductor)->inductorCurrentOffsets[ j ] << ", ";
      }
      Xyce::dout() << "} " << std::endl;
      i++;
      currentInductor++;
    }

    Xyce::dout() << "mEquVPosOffset = " << mEquVPosOffset << "\tmEquVNegOffset = " << mEquVNegOffset << std::endl;
    Xyce::dout() << "mEquInductorOffsets = ";
    for(i=0;i<numInductors; ++i)
    {
      Xyce::dout() << mEquInductorOffsets[i] << ", ";
    }
    Xyce::dout() << std::endl
      << "mEquMOffset = " << mEquMOffset << "\tmEquROffset = " << mEquROffset  << std::endl;

    Xyce::dout() << "rEquInductorOffsets = ";
    for(i=0;i<numInductors; ++i)
    {
      Xyce::dout() << rEquInductorOffsets[i] << ", ";
    }
    Xyce::dout() << std::endl
      << "rEquROffset = " << rEquROffset << std::endl;
  }
#endif
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

  std::vector< InductorInstanceData* >::iterator currentData = instanceData.begin();
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
    Xyce::dout() << "Instance::updateIntermediateVars " << std::endl;
  }
#endif

  N_LAS_Vector & solVector = *(extData.nextSolVectorPtr);

  // some parameters in the model class that we will use often
  const double A      = model_.A;
  const double Alpha  = model_.Alpha;
  const double Area   = model_.Area;
  const double BetaH  = model_.BetaH;
  const double BetaM  = model_.BetaM;
  const double C      = model_.C;
  const double DeltaVScaling = model_.DeltaVScaling;
  const double Gap    = model_.Gap;
  const double Ms     = model_.Ms;
  const double Kirr   = model_.Kirr;
  const double Path   = model_.Path;

  const double mVarScaling = model_.mVarScaling;
  const double rVarScaling = model_.rVarScaling;
  const double mEqScaling = model_.mEqScaling;
  const double rEqScaling = model_.rEqScaling;

  // calculate the voltage drop over the first inductor
  // as this is needed later
  double Vpos = solVector[(instanceData[0])->li_Pos];
  double Vneg = solVector[(instanceData[0])->li_Neg];

  // voltage drop over first inductor.
  double voltageDrop= Vpos - Vneg;

  // only update maxVoltageDrop when system has converged or we may
  // get wildly wrong values.
  N_LAS_Vector & lastSolVector = *(extData.currSolVectorPtr);
  double lastVoltageDrop = lastSolVector[(instanceData[0])->li_Pos] - lastSolVector[(instanceData[0])->li_Neg];
  if ( (getSolverState().newtonIter == 0) && (fabs(lastVoltageDrop) > maxVoltageDrop) )
  {
    maxVoltageDrop=fabs(lastVoltageDrop);
    Xyce::dout() << std::endl << " maxVoltageDrop = " << maxVoltageDrop << std::endl;
  }

  // approximate the sgn( voltageDrop ) with
  // tanh ( scalefactor * voltageDrop / maxVoltageDrop )
  if( model_.UseConstantDeltaVScaling )
  {
    qV = DeltaVScaling * voltageDrop;
  }
  else
  {
    qV = DeltaVScaling * voltageDrop / maxVoltageDrop;
  }

  double tanh_qV = 0.0;

  if ( (fabs(qV) < CONSTTANH_THRESH) )
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

  Happ = 0.0;
  std::vector< InductorInstanceData* >::iterator currentInductor = instanceData.begin();
  std::vector< InductorInstanceData* >::iterator endInductor = instanceData.end();
  int il=0;
  while( currentInductor != endInductor )
  {
    Happ += solVector[(*currentInductor)->li_Branch] * inductanceVals[ il ];
    il++;
    currentInductor++;
  }
  Happ /= Path;

  double latestMag = mVarScaling * solVector[ li_MagVar ];
  if( model_.factorMS )
  {
    latestMag *= Ms;
  }

  double H = Happ - (Gap / Path) * latestMag;

  He = H + Alpha * latestMag;

  Heo = BetaH*A;

  // terms that come up frequently
  const double gap_path = Gap / Path;
  const double He2 = He*He;
  const double Heo2 = Heo*Heo;
  const double sq_Heo2He2 = sqrt(Heo2 + He2);

  delM0 = model_.BetaM * Ms;
  double Man = Ms * He / ( A + sq_Heo2He2 );
  delM = Man - latestMag;

  // terms that come up frequently
  const double delM2 = delM*delM;
  const double delM02 = delM0*delM0;
  const double sq_delM02delM2 = sqrt( delM02 + delM2 );

  if( model_.factorMS )
  {
    Mirrp = (delM * tanh_qV + sq_delM02delM2 ) / (2*( Kirr- Alpha * sq_delM02delM2));
    Manp =  Ms * (A + Heo2/sq_Heo2He2) / pow(A + sq_Heo2He2, 2.0);
    P = ( C * Manp + (1 - C) * Mirrp) / ((1 + (gap_path - Alpha) * C * Manp + gap_path * (1-C) * Mirrp)*Ms);
  }
  else
  {
    Mirrp = (delM * tanh_qV + sq_delM02delM2 ) / (2*( Kirr- Alpha * sq_delM02delM2));
    Manp =  Ms * (A + Heo2/sq_Heo2He2) / pow(A + sq_Heo2He2, 2.0);
    P = ( C * Manp + (1 - C) * Mirrp) / (1 + (gap_path - Alpha) * C * Manp + gap_path * (1-C) * Mirrp);
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0  && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "\tA = " << A << std::endl
         << "\tArea = " << Area << std::endl
         << "\tPath = " << Path << std::endl
         << "\tGap = " << Gap << std::endl
         << "\tC = " << C << std::endl
         << "\tVpos = " << Vpos << std::endl
         << "\tVneg = " << Vneg << std::endl
         << "\tvoltageDrop = " << voltageDrop << std::endl
         << "\tqV = " << qV << std::endl
         << "\tdelM0 = " << delM0 << std::endl
         << "\tHapp = " << Happ
         << "\tlatestMag = " << latestMag
         << "\tlatestR = " <<  rVarScaling * solVector[ li_RVar ] << std::endl
         << "\tHe = " << He << std::endl
         << "\tH = " << H << std::endl
         << "\tHeo = " << Heo << std::endl
         << "\tMan = " << Man << std::endl
         << "\tdelM = " << delM << std::endl
         << "\tMirrp = " << Mirrp << std::endl
         << "\tManp = " << Manp << std::endl
         << "\tP  = " << P << std::endl
         << "\tgetSolverState().newtonIter = " << getSolverState().newtonIter << std::endl
         << std::endl;
  }
#endif

  // now calculate important derivative quantities

  double dHe_dM =  ((Alpha - gap_path) * mVarScaling);

  double dManp_dM = ( -Ms * He / (pow(A + sq_Heo2He2, 2.0)*sq_Heo2He2)) *
                    ( (Heo2 / (Heo2 + He2)) + (2.0*(A + Heo2 / sq_Heo2He2)/(A+sq_Heo2He2)) ) * dHe_dM;

  double ddelM_dM = ( dHe_dM*Ms/(A + sq_Heo2He2) ) * (1.0 - He2 / ((A + sq_Heo2He2)*sq_Heo2He2)) - mVarScaling;

  double dMirrp_dM = (1.0/(2.0*(Kirr - Alpha*sq_delM02delM2))) *
                     (tanh_qV + delM/sq_delM02delM2 +
                       (2.0*Alpha*delM*(delM*tanh_qV + sq_delM02delM2)
                     /(2.0*(Kirr-Alpha*sq_delM02delM2)*sq_delM02delM2))) * ddelM_dM;

  double dP_Denom=0.0;
  if( model_.factorMS )
  {
    dP_Denom = 1.0 + (gap_path - Alpha)*C*Manp + gap_path * (1.0-C) * Mirrp;

    dP_dM = (1.0/dP_Denom) * (C * dManp_dM + (1.0-C) * dMirrp_dM) -
              ( (C*Manp + (1.0-C)*Mirrp)/pow(dP_Denom,2.0) ) *
                ( (gap_path - Alpha)*C*dManp_dM + gap_path*(1.0-C)*dMirrp_dM );
    dP_dM /= Ms;
  }
  else
  {
    dP_Denom = 1.0 + (gap_path - Alpha)*C*Manp + gap_path * (1.0-C) * Mirrp;

    dP_dM = (1.0/dP_Denom) * (C * dManp_dM + (1.0-C) * dMirrp_dM) -
              ( (C*Manp + (1.0-C)*Mirrp)/pow(dP_Denom,2.0) ) *
                ( (gap_path - Alpha)*C*dManp_dM + gap_path*(1.0-C)*dMirrp_dM );
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "\tA = " << A << std::endl
      << "\tAlpha = " << Alpha << std::endl
      << "\tC = " << C << std::endl
      << "\tGap = " << Gap << std::endl
      << "\tMs = " << Ms << std::endl
      << "\tKirr = " << Kirr << std::endl
      << "\tPath = " << Path << std::endl
      << "\tHe2 = " << He2 << std::endl
      << "\tHeo2 = " << Heo2 << std::endl
      << "\tdelM2 = " << delM2 << std::endl
      << "\tdelM02 = " << delM02 << std::endl
      << "\tdHe_dM = " << dHe_dM << std::endl
      << "\tdManp_dM = " << dManp_dM << std::endl
      << "\tddelM_dM = " << ddelM_dM << std::endl
      << "\tdMirrp_dM = " << dMirrp_dM << std::endl
      << "\tdP_dM = " << dP_dM << std::endl
      << "\tdenom 1+(1-lg/lt)P = " << (1+(1-Gap/Path)*P) << std::endl;
  }
#endif

  //    % Now find (dP/dI_i): (this is nearly identical to dP/dM)
  currentInductor = instanceData.begin();

  for( int i=0; i<numInductors; ++i )
  {

    dHe_dI[ i ] = inductanceVals[ i ] / Path;
    dManp_dI[i] = ( -Ms * He / (pow(A + sq_Heo2He2, 2.0)*sq_Heo2He2)) *
                   ( (Heo2 / (Heo2 + He2)) + (2.0*(A + Heo2 / sq_Heo2He2)/(A+sq_Heo2He2)) ) * dHe_dI[i];
    ddelM_dI[i] = (Ms / (A + sq_Heo2He2)) * (1.0 - He2/((A + sq_Heo2He2)*sq_Heo2He2)) * dHe_dI[i];
    dMirrp_dI[i] = (1.0/(2.0*(Kirr - Alpha*sq_delM02delM2))) *
                   (tanh_qV + delM/sq_delM02delM2 +
                     (2.0*Alpha*delM*(delM*tanh_qV +
                       sq_delM02delM2)/(2.0*(Kirr-Alpha*sq_delM02delM2)*sq_delM02delM2))) * ddelM_dI[i];
    dP_dI[i] = (1.0/dP_Denom) * (C * dManp_dI[i] + (1.0-C) * dMirrp_dI[i]) -
          ( (C*Manp + (1.0-C)*Mirrp)/pow(dP_Denom,2.0) ) *
            ( (gap_path - Alpha)*C*dManp_dI[i] + gap_path*(1.0-C)*dMirrp_dI[i] );

  if( model_.factorMS )
  {
    dP_dI[i] /= Ms;
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
       Xyce::dout() << "\tdHe_dI[ " << i << " ] =" << dHe_dI[ i ] << std::endl
            << "\tdManp_dI[ " << i << " ] = " << dManp_dI[i] << std::endl
            << "\tddelM_dI[ " << i << " ] = " << ddelM_dI[i] << std::endl
            << "\tMirrp_dI[ " << i << " ] = " << dMirrp_dI[i] << std::endl
            << "\tdP_dI[ " << i << " ] = " << dP_dI[i] << std::endl;
    }
#endif
    currentInductor++;
  }

  // Now find (dP/dV_1):
  double dMirrp_dVp = (delM * DeltaVScaling * (1.0-pow(tanh_qV,2.0))) /
                      (2.0 * (Kirr - Alpha * sq_delM02delM2));
  double dMirrp_dVn = -dMirrp_dVp;

  dP_dVp = (1.0/dP_Denom) * ((1.0-C) * dMirrp_dVp) -
            ( (C*Manp + (1.0-C)*Mirrp)/pow(dP_Denom,2.0) ) * (  gap_path*(1.0-C)*dMirrp_dVp );

  if( model_.factorMS )
  {
    dP_dVp /= Ms;
  }

  dP_dVn = -dP_dVp;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "\tdMirrp_dVp = " << dMirrp_dVp << std::endl
      << "\tdMirrp_dVn = " << dMirrp_dVn << std::endl
      << "\tdP_dVp = " << dP_dVp << std::endl
      << "\tdP_dVn = " << dP_dVn << std::endl;
  }
#endif

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
  std::vector< InductorInstanceData* >::iterator
    currentInductor = instanceData.begin();
  std::vector< InductorInstanceData* >::iterator
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
bool Instance::updatePrimaryState ()
{
  bool bsuccess = true;
#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "Instance::updatePrimaryState---------------" << std::endl
         << "\tname = " << getName() << std::endl;
  }
#endif

  updateIntermediateVars ();

  N_LAS_Vector & solVector = *(extData.nextSolVectorPtr);
  N_LAS_Vector & staVector = *(extData.nextStaVectorPtr);
  N_LAS_Vector & stoVector = *(extData.nextStoVectorPtr);
  double mVarScaling = model_.mVarScaling;
  double rVarScaling = model_.rVarScaling;

  // place current values of mag, H and R in state vector
  // must unscale them as the rest of the class assumes
  // that these aren't scaled yet.
  staVector[ li_MagVarState ] = solVector[ li_MagVar ];
  stoVector[ li_RVarStore ] = solVector[ li_RVar ];

  double latestMag = mVarScaling * solVector[ li_MagVar ];
  if( model_.factorMS )
  {
    latestMag *= model_.Ms;
  }

  // B and H are quantities that we can calculate from M and R.  Store them in the state vector
  // for output if the user requests it.
  stoVector[ li_HVarStore ] = model_.HCgsFactor * (Happ  - (model_.Gap / model_.Path) * latestMag);
  stoVector[ li_BVarStore ] = model_.BCgsFactor * (4.0e-7 * M_PI * (stoVector[ li_HVarStore ] + latestMag));

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
#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "Instance::updateSecondaryState-------------" << std::endl
         << "\tname = " << getName() << std::endl;
  }
#endif

  N_LAS_Vector & staVector = *(extData.nextStaVectorPtr);
  N_LAS_Vector & staDerivVec = *(extData.nextStaDerivVectorPtr);

  // copy derivitive of Mag from result vector into state vector
  staVector[ li_MagVarDerivState ] = staDerivVec[ li_MagVarState ];

  return bsuccess;
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
  double mVarScaling = model_.mVarScaling;
  double rVarScaling = model_.rVarScaling;
  double mEqScaling = model_.mEqScaling;
  double rEqScaling = model_.rEqScaling;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "Instance::loadDAEQVector------------------------" << std::endl
         << "\tname = " << getName() << std::endl;
  }
#endif

  N_LAS_Vector & staVector = *(extData.nextStaVectorPtr);
  N_LAS_Vector & solVector = *(extData.nextSolVectorPtr);
  double * qVec = extData.daeQVectorRawPtr;

  // update LOI -- the following product
  // I = column vector of currents
  // L = row vector of inductances
  // LO = matrix = mutualCup * sqrt( L' * L )
  // LOI = column vector = mutualCup * sqrt( L' * L ) * I
  // LOI[1] = mutualCup * sqrt(L[1]*L[1])*I[1]) +
  //          mutualCup * sqrt(L[1]*L[2])*I[2]) + ...
  //          mutualCup * sqrt(L[1]*L[n])*I[n])

  std::vector< InductorInstanceData* >::iterator currentInductor = instanceData.begin();
  std::vector< InductorInstanceData* >::iterator endInductor = instanceData.end();
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

    qVec[((*currentInductor)->li_Branch)] += LOI[ i ];

    double current = inductorCurrents[ i ];
    double windings = (*currentInductor)->L;

    qVec[ li_RVar ] += rEqScaling * current * windings;
    i++;
    currentInductor++;
  }

  double latestMag = mVarScaling * staVector[ li_MagVarState ];

  // M equation
  if(!getSolverState().dcopFlag)
  {

    qVec[ li_MagVar ] += mEqScaling * latestMag;
  }

  return bsuccess;
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
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  bool bsuccess=true;
  double mVarScaling = model_.mVarScaling;
  double rVarScaling = model_.rVarScaling;
  double mEqScaling = model_.mEqScaling;
  double rEqScaling = model_.rEqScaling;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "Instance::loadDAEFVector------------------------" << std::endl
         << "\tname = " << getName() << std::endl;
  }
#endif

  N_LAS_Vector & staVector = *(extData.nextStaVectorPtr);
  N_LAS_Vector & solVector = *(extData.nextSolVectorPtr);
  N_LAS_Vector & staDerivVec = *(extData.nextStaDerivVectorPtr);
  N_LAS_Vector & stoVector = *(extData.nextStoVectorPtr);

  double * fVec = extData.daeFVectorRawPtr;

  double latestR = rVarScaling * stoVector[ li_RVarStore ];

  if(getSolverState().dcopFlag)
  {
    //enforce R = 0 in dc op
    latestR = 0.0;
  }

  // for the M equation

  fVec[li_MagVar] -= mEqScaling *  P * latestR / (model_.Path);

  // if |P| is near zero, then the M equation becomes dM/dt = 0, or M is
  // constant.  In this case we'll add a diagonal element for M so that
  // sole dM/dt element in dQ/dX doesn't cause a time step too small error
  // Since P is normally very large, we'll test for |P| <= 1.0.
  if( fabs( P ) <= model_.pZeroTol )
  {

    fVec[li_MagVar] -= mVarScaling * staVector[ li_MagVarState ];
  }

  // for the R equation
  fVec[li_RVar] -= rEqScaling * rVarScaling * stoVector[ li_RVarStore ];

  // used in scaling the branch equation;
  double mid=1.0;
  if( model_.factorMS )
  {
    mid = 1.0 + (1.0 - ((model_.Gap) / (model_.Path)))*P*(model_.Ms);
  }
  else
  {
    mid = 1.0 + (1.0 - ((model_.Gap) / (model_.Path)))*P;
  }

  // loop over each inductor and load it's F vector components
  std::vector< InductorInstanceData* >::iterator currentInductor = instanceData.begin();
  std::vector< InductorInstanceData* >::iterator endInductor = instanceData.end();
  int i=0;
  while( currentInductor != endInductor )
  {
    double current   = solVector[(*currentInductor)->li_Branch];
    double vNodePos  = solVector[(*currentInductor)->li_Pos];
    double vNodeNeg  = solVector[(*currentInductor)->li_Neg];


    fVec[((*currentInductor)->li_Pos)]    +=  scalingRHS * current;

    fVec[((*currentInductor)->li_Neg)]    += -scalingRHS * current;

    fVec[((*currentInductor)->li_Branch)] += -((vNodePos - vNodeNeg)/mid);
    double windings = (*currentInductor)->L;

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "  Inductor = " << (*currentInductor)->name
           << " li_Pos = " << (*currentInductor)->li_Pos
           << " li_Neg = " << (*currentInductor)->li_Neg
           << " li_Branch = " << (*currentInductor)->li_Branch
           << "\tPos/Neg current*windings = " << scalingRHS*current*windings
           << "\tBranch = " << ((vNodePos - vNodeNeg)/mid)
           << std::endl;
    }
#endif
    currentInductor++;
    i++;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx
//
// Purpose       : Loads the Q-vector contributions for a single instance.
//
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  bool bsuccess = true;

  double mVarScaling = model_.mVarScaling;
  double rVarScaling = model_.rVarScaling;
  double mEqScaling = model_.mEqScaling;
  double rEqScaling = model_.rEqScaling;

  int i;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "Instance::loadDAEdQdx-----------------------" << std::endl
         << "\tname = " << getName() << std::endl;
  }
#endif

  N_LAS_Matrix * dQdxMatPtr = extData.dQdxMatrixPtr;

  // update M equation
  if(!getSolverState().dcopFlag)
  {

    (*dQdxMatPtr)[li_MagVar][mEquMOffset] += mEqScaling * mVarScaling;
  }

  // update the R equation
  for( i = 0; i< numInductors; i++ )
  {

    (*dQdxMatPtr)[li_RVar][rEquInductorOffsets[i] ] += rEqScaling * inductanceVals[i];
  }


  // loop over each inductor and load it's Q vector components
  std::vector< InductorInstanceData* >::iterator currentInductor = instanceData.begin();
  std::vector< InductorInstanceData* >::iterator endInductor = instanceData.end();
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
// Purpose       : Loads the F-vector contributions for a single instance.
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
  double mVarScaling = model_.mVarScaling;
  double rVarScaling = model_.rVarScaling;
  double mEqScaling = model_.mEqScaling;
  double rEqScaling = model_.rEqScaling;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "Instance::loadDAEdFdx----------------------" << std::endl
         << "\tname = " << getName() << std::endl;
  }
#endif

  N_LAS_Vector & staVector = *(extData.nextStaVectorPtr);
  N_LAS_Vector & solVector = *(extData.nextSolVectorPtr);
  N_LAS_Vector & staDerivVec = *(extData.nextStaDerivVectorPtr);
  N_LAS_Vector & stoVector = *(extData.nextStoVectorPtr);
  N_LAS_Matrix * dFdxMatPtr = extData.dFdxMatrixPtr;

  // udate dependent parameters
  //updateIntermediateVars();

  // pull these parameters up from the model class to make it easier
  // to view the equations.
  const double Gap = model_.Gap;
  const double Path = model_.Path;

  // terms that come up frequently
  double latestR = rVarScaling * stoVector[ li_RVarStore ];

  // terms for the M equation
  if(!getSolverState().dcopFlag)
  {

    (*dFdxMatPtr)[ li_MagVar ][ mEquMOffset ]    -= mEqScaling * dP_dM * latestR / Path;   // d/dM

    (*dFdxMatPtr)[ li_MagVar ][ mEquROffset ]    -= mEqScaling * P * rVarScaling / Path;   // d/dR


    (*dFdxMatPtr)[ li_MagVar ][ mEquVPosOffset ] -= mEqScaling * dP_dVp * latestR / Path;  // d/dV_+

    (*dFdxMatPtr)[ li_MagVar ][ mEquVNegOffset ] -= mEqScaling * dP_dVn * latestR / Path;  // d/dV_-
    for( int i = 0; i<numInductors; ++i)
    {

      (*dFdxMatPtr)[ li_MagVar ][mEquInductorOffsets[i] ] -=
         mEqScaling * dP_dI[i] * latestR / Path;                     // d/dI_i;
    }
  }
  else
  {
    // the above load for the M equation is basically zero in the dc op.  We
    // need something on the diagonal for M to make the matrix non-singular

    (*dFdxMatPtr)[ li_MagVar ][ mEquMOffset ]    += getSolverState().pdt;
  }

  // if |P| is near zero, then the M equation becomes dM/dt = 0, or M is
  // constant.  In this case we'll add a diagonal element for M so that
  // sole dM/dt element in dQ/dX doesn't cause a time step too small error
  // Since P is normally very large, we'll test for |P| <= 1.0.
  if( fabs( P ) <= model_.pZeroTol )
  {

    (*dFdxMatPtr)[ li_MagVar ][ mEquMOffset ] += 1.0;
  }

  // update the R equation

  (*dFdxMatPtr)[ li_RVar ][rEquROffset] -= rEqScaling * rVarScaling;

  // loop over each inductor and load it's dFdx components
  double mid=0.0;
  if( model_.factorMS )
  {
    mid = 1.0 + (1.0 - ((model_.Gap) / (model_.Path)))*P*(model_.Ms);
  }
  else
  {
    mid = 1.0 + (1.0 - ((model_.Gap) / (model_.Path)))*P;
  }

  std::vector< InductorInstanceData* >::iterator currentInductor = instanceData.begin();
  std::vector< InductorInstanceData* >::iterator endInductor = instanceData.end();
  while( currentInductor != endInductor )
  {
    // do the normal work for an inductor

    (*dFdxMatPtr)[((*currentInductor)->li_Pos)]   [((*currentInductor)->APosEquBraVarOffset)]  +=  scalingRHS;

    (*dFdxMatPtr)[((*currentInductor)->li_Neg)]   [((*currentInductor)->ANegEquBraVarOffset)]  += -scalingRHS;

    (*dFdxMatPtr)[((*currentInductor)->li_Branch)][((*currentInductor)->ABraEquPosNodeOffset)] += -1.0/mid;

    (*dFdxMatPtr)[((*currentInductor)->li_Branch)][((*currentInductor)->ABraEquNegNodeOffset)] +=  1.0/mid;

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      Xyce::dout()
       << "(*currentInductor)->li_Pos = " << (*currentInductor)->li_Pos << std::endl
       << "(*currentInductor)->li_Neg = " << (*currentInductor)->li_Neg << std::endl
       << "(*currentInductor)->li_Branch = " << (*currentInductor)->li_Branch << std::endl
       << "(*currentInductor)->APosEquBraVarOffset = " << (*currentInductor)->APosEquBraVarOffset << std::endl
       << "(*currentInductor)->ANegEquBraVarOffset = " << (*currentInductor)->ANegEquBraVarOffset << std::endl
       << "(*currentInductor)->ABraEquPosNodeOffset = " << (*currentInductor)->ABraEquPosNodeOffset << std::endl
       << "(*currentInductor)->ABraEquNegNodeOffset = " << (*currentInductor)->ABraEquNegNodeOffset << std::endl
       << "(*dFdxMatPtr)["<<((*currentInductor)->li_Pos)<<"]   ["<<((*currentInductor)->APosEquBraVarOffset)<<"] =  " << scalingRHS << std::endl
       << "(*dFdxMatPtr)["<<((*currentInductor)->li_Neg)<<"]   ["<<((*currentInductor)->ANegEquBraVarOffset)<<"]  =  " << -scalingRHS << std::endl
       << "(*dFdxMatPtr)["<<((*currentInductor)->li_Branch)<<"]["<<((*currentInductor)->ABraEquPosNodeOffset)<<"] = " << -1/mid << std::endl
       << "(*dFdxMatPtr)["<<((*currentInductor)->li_Branch)<<"]["<<((*currentInductor)->ABraEquNegNodeOffset)<<"] = " << 1/mid << std::endl;
    }
#endif

    double delV = solVector[(*currentInductor)->li_Pos] - solVector[(*currentInductor)->li_Neg];

    for( int j = 0; j<numInductors; ++j )
    {

      (*dFdxMatPtr)[((*currentInductor)->li_Branch)][(*currentInductor)->inductorCurrentOffsets[j]] +=
        delV * (1.0 - (Gap/Path)) * dP_dI[j]/(mid*mid);

#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
        {
          Xyce::dout() << "(*dFdxMatPtr)[((*currentInductor)->li_Branch)][(*currentInductor)->inductorCurrentOffsets[j]] =  " << delV * (1 - (Gap/Path)) * dP_dI[j]/(mid*mid) << std::endl;
        }
#endif
    }

    (*dFdxMatPtr)[(*currentInductor)->li_Branch][(*currentInductor)->magOffset]  += delV * (1.0 - (Gap/Path)) * dP_dM/(mid*mid);

    (*dFdxMatPtr)[(*currentInductor)->li_Branch][(*currentInductor)->vPosOffset] += delV * (1.0 - (Gap/Path)) * dP_dVp/(mid*mid);

    (*dFdxMatPtr)[(*currentInductor)->li_Branch][(*currentInductor)->vNegOffset] += delV * (1.0 - (Gap/Path)) * dP_dVn/(mid*mid);

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "(*dFdxMatPtr)[(*currentInductor)->li_Branch][(*currentInductor)->magOffset] =  " << delV * (1 - (Gap/Path)) * dP_dM/(mid*mid) << std::endl
       << "(*dFdxMatPtr)[(*currentInductor)->li_Branch][(*currentInductor)->vPosOffset]  =  " << delV * (1 - (Gap/Path)) * dP_dVp/(mid*mid) << std::endl
       << "(*dFdxMatPtr)[(*currentInductor)->li_Branch][(*currentInductor)->vNegOffset] = " << delV * (1 - (Gap/Path)) * dP_dVn/(mid*mid) << std::endl;
    }
#endif
    currentInductor++;
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
    N_LAS_Vector & staVector = *(extData.nextStaVectorPtr);
    N_LAS_Vector & stoVector = *(extData.nextStoVectorPtr);
    double mVarScaling = model_.mVarScaling;
    double rVarScaling = model_.rVarScaling;

    double latestMag = mVarScaling * staVector[ li_MagVarState ];
    if( model_.factorMS )
    {
      latestMag *= model_.Ms;
    }
    double latestR   = rVarScaling * stoVector[ li_RVarStore ];
    (*outputFileStreamPtr)
      << getSolverState().currTime << "  "
      << latestMag << "\t  "
      << latestR << "\t "
      << staVector[ li_BVarStore ] << "\t "
      << staVector[ li_HVarStore ]
      << std::endl;

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
void Instance::varTypes( std::vector<char> & varTypeVec )
{
  varTypeVec.resize(numInductors+2);
  for(int i=0; i<numInductors; i++)
  {
    varTypeVec[i] = 'I';
  }
  // I don't know what should be used for non I,V vars.
  varTypeVec[numInductors] = 'M';
  varTypeVec[numInductors+1] = 'R';
}

//-----------------------------------------------------------------------------
// Function      : Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
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
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 03/21/2005
//-----------------------------------------------------------------------------
Model::Model(
  const Configuration & configuration,
  const ModelBlock &    MB,
  const FactoryBlock &  factory_block)
  : DeviceModel(MB, configuration.getModelParameters(), factory_block),
    A(0.0),
    Alpha(0.0),
    Area(0.0),
    BetaH(0.0),
    BetaM(0.0),
    C(0.0),
    DeltaVScaling(0.0),
    Gap(0.0),
    Kirr(0.0),
    Ms(0.0),
    Path(0.0),
    tempCoeff1(0.0),
    tempCoeff2(0.0),
    outputStateVars(0),
    factorMS(0),
    BCgsFactor( 10000.0 ),
    HCgsFactor( 0.012566370614359 ),  // 4 pi / 1000
    UseConstantDeltaVScaling( true ),
    tnom(getDeviceOptions().tnom)
{
  setLevel(1);


  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to .model line and constant defaults from metadata:
  setModParams (MB.params);

  // scale gap, path and area from cm and cm^2 to m and m^2
  Gap *= 1.0e-2;
  Path *= 1.0e-2;
  Area *= 1.0e-4;

  if( BHSiUnits != 0 )
  {
    // user requested SI units over the default of CGS units.  Change
    // conversion factor to unity.
    BCgsFactor=1.0;
    HCgsFactor=1.0;
  }

  // Set any non-constant parameter defaults:
  // when Ms factoring is off, scaling of M/R is still needed.
  if( factorMS == 0 )
  {
    mVarScaling=1.0e3;
    rVarScaling=1.0e3;
    mEqScaling=1.0e-3;
    rEqScaling=1.0e-3;
  }


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
  std::vector<Instance*>::iterator iter;
  std::vector<Instance*>::iterator first = instanceContainer.begin();
  std::vector<Instance*>::iterator last  = instanceContainer.end();

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
  std::vector<Instance*>::const_iterator iter;
  std::vector<Instance*>::const_iterator first = instanceContainer.begin();
  std::vector<Instance*>::const_iterator last  = instanceContainer.end();

  int i, isize;
  isize = instanceContainer.size();

  os << std::endl;
  os << "Number of MutIndNonLin instances: " << isize << std::endl;
  os << "    name=\t\tmodelName\tParameters" << std::endl;
  for (i=0, iter=first; iter!=last; ++iter, ++i)
  {
    os << "  " << i << ": " << (*iter)->getName() << "\t";
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
// Function      : Master::updateState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/25/08
//-----------------------------------------------------------------------------
bool Master::updateState (double * solVec, double * staVec, double * stoVec)
{
  bool bsuccess = true;
  bool tmpBool = true;

  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    tmpBool = (*it)->updatePrimaryState ();
    bsuccess = bsuccess && tmpBool;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Master::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/25/08
//-----------------------------------------------------------------------------
bool Master::updateSecondaryState (double * staDerivVec, double * stoVec)
{
  bool bsuccess = true;
  bool tmpBool = true;

  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    tmpBool = (*it)->updateSecondaryState ();
    bsuccess = bsuccess && tmpBool;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Master::loadDAEVectors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/25/08
//-----------------------------------------------------------------------------
bool Master::loadDAEVectors (double * solVec, double * fVec, double *qVec,  double * storeLeadF, double * storeLeadQ)
{
  bool bsuccess = true;
  bool tmpBool = true;

  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    tmpBool = (*it)->loadDAEFVector();
    bsuccess = bsuccess && tmpBool;
    tmpBool = (*it)->loadDAEQVector();
    bsuccess = bsuccess && tmpBool;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Master::loadDAEMatrices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/25/08
//-----------------------------------------------------------------------------
bool Master::loadDAEMatrices (N_LAS_Matrix & dFdx, N_LAS_Matrix & dQdx)
{
  bool bsuccess = true;
  bool tmpBool = true;

  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    tmpBool = (*it)->loadDAEdFdx ();
    bsuccess = bsuccess && tmpBool;
    tmpBool = (*it)->loadDAEdQdx ();
    bsuccess = bsuccess && tmpBool;
  }

  return bsuccess;
}

Device *Traits::factory(const Configuration &configuration, const FactoryBlock &factory_block)
{

  return new Master(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);
}

void registerDevice()
{
  Config<Traits>::addConfiguration()
    .registerDevice("min", 1)
    .registerModelType("min", 1)
    .registerModelType("core", 1);
}

} // namespace MutIndNonLin
} // namespace Device
} // namespace Xyce

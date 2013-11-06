// ----------------------------------------------------------------------------
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
// ----------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_RxnSet.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL
//
// Creation Date  : 02/09/08
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.14.2.3 $
// Revision Date  : $Date: 2013/10/03 17:23:38 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

#include <N_UTL_Misc.h>
#include <iostream>

#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif

#ifdef Xyce_DEBUG_DEVICE
#ifdef HAVE_CSTDIO
#include <cstdio>
#else
#include <stdio.h>
#endif
#endif


// ----------    Xyce Includes  ----------
#include <N_DEV_Const.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_Region.h>
#include <N_DEV_RegionData.h>
#include <N_DEV_SpecieSource.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_MultiVector.h>
#include <N_LAS_Vector.h>

#include <N_DEV_RxnSet.h>
#include <N_DEV_ReactionNetwork.h>

#include <N_UTL_Expression.h>

namespace Xyce {
namespace Device {

template<>
ParametricData<RxnSet::Instance>::ParametricData()
{
  setNumNodes(2);
  setNumOptionalNodes(0);
  setNumFillNodes(0);
  setModelRequired(1);
  addModelType("PN");
  addModelType("NP");
  addModelType("RXN");

  addPar ("TEMP", 0.0, false, ParameterType::TIME_DEP, &RxnSet::Instance::TEMP, NULL, STANDARD, CAT_NONE, "");

  // user-specified scaling vars:
  addPar ("X0", 1.0e-7, false, ParameterType::NO_DEP, &RxnSet::Instance::x0_user, NULL, U_NONE, CAT_NONE, "Length Scalar.");
  addPar ("C0", 1.0e+12, false, ParameterType::NO_DEP, &RxnSet::Instance::C0_user, NULL, U_NONE, CAT_NONE, "Concentration Scalar.");
  addPar ("t0", 1.0e-6, false, ParameterType::NO_DEP, &RxnSet::Instance::t0_user, NULL, U_NONE, CAT_NONE, "Time Scalar.");
  addPar ("outputXscalar", 1.0, false, ParameterType::NO_DEP, &RxnSet::Instance::outputXscalar, NULL, U_NONE, CAT_NONE, "Scalar for X axis in tecplot file (default unit is cm)");
  addPar ("OUTPUTINTERVAL", 0.0, false, ParameterType::NO_DEP, &RxnSet::Instance::outputInterval, NULL, U_NONE, CAT_NONE, "Output Interval(sec)");

  // Set up map for non-double precision variables:

  addPar ("COPYRXN", false, false, ParameterType::NO_DEP, &RxnSet::Instance::reactionFileCopyFlag, NULL, U_LOGIC, CAT_NONE, "Flag for processing chemistry file only once, then copied to other mesh points.");
  addPar ("SCALERXN", true, false, ParameterType::NO_DEP, &RxnSet::Instance::useScaledVariablesFlag, NULL, U_LOGIC, CAT_NONE, "Flag for applying scaling to the reaction equations.");
  addPar ("DIFFUSION", false, false, ParameterType::NO_DEP, &RxnSet::Instance::diffusionFlag, &RxnSet::Instance::diffusionFlagGiven, U_LOGIC, CAT_NONE, "Flag for enabling lattice defect diffusion.");
  addPar ("TRANSPORT", false, false, ParameterType::NO_DEP, &RxnSet::Instance::transportFlag, &RxnSet::Instance::transportFlagGiven, U_LOGIC, CAT_NONE, "Flag for enabling lattice defect diffusion.  Identical to DIFFUSION flag, above.  Do not set both!");
  addPar ("EXCLUDENOSOURCE", true, false, ParameterType::NO_DEP, &RxnSet::Instance::excludeNoSourceRegionsFlag, &RxnSet::Instance::excludeNoSourceRegionsFlagGiven, U_LOGIC, CAT_NONE, "Flag for excluding regions that are outside of source region from computing defect reaction equations.  This is a speed optimization.  Turning it on will NOT change the answer");
  addPar ("COLUMNREORDER", false, false, ParameterType::NO_DEP, &RxnSet::Instance::columnReorderingFlag, NULL, U_LOGIC, CAT_NONE, "Debug Flag for turning on/off column reordering.");
  addPar ("OUTPUTREGION", 1, false, ParameterType::NO_DEP, &RxnSet::Instance::outputRegion, NULL); addPar ("TECPLOTLEVEL", 0, false, ParameterType::NO_DEP, &RxnSet::Instance::tecplotLevel, NULL, U_NONE, CAT_NONE, "Integer number to determine type of tecplot output.  0=no output.  1=single time-dependent file, with each time step in a different zone.");
  addPar ("DIRICHLETBC", false, false, ParameterType::NO_DEP, &RxnSet::Instance::dirichletBCFlag, NULL, U_LOGIC, CAT_NONE, "Flag for using Dirichlet boundary conditions.");
}

template<>
ParametricData<RxnSet::Model>::ParametricData()
{
  addPar ("TNOM", 0.0, false, ParameterType::NO_DEP, &RxnSet::Model::TNOM, NULL, U_DEGC, CAT_UNKNOWN,
          "Parameter measurement temperature");

  // rxn stuff:

  addPar ("XLO", 1.0e-5, false, ParameterType::NO_DEP, &RxnSet::Model::xlo, NULL, U_CM, CAT_UNKNOWN, "Left edge of integration volume.");
  addPar ("XHI", 3.0e-4, false, ParameterType::NO_DEP, &RxnSet::Model::xhi, NULL, U_CM, CAT_UNKNOWN, "Right edge of integration volume");
  addPar ("XLO_SOURCE",1.0e-5, false, ParameterType::NO_DEP, &RxnSet::Model::xlo_source, &RxnSet::Model::xlo_sourceGiven, U_CM, CAT_UNKNOWN, "Left edge of source region");
  addPar ("XHI_SOURCE",3.0e-4, false, ParameterType::NO_DEP, &RxnSet::Model::xhi_source, &RxnSet::Model::xhi_sourceGiven, U_CM, CAT_UNKNOWN, "Right edge of source region");
  addPar ("MASTERSOURCE", 0.0, false, ParameterType::TIME_DEP, &RxnSet::Model::masterSource, NULL, STANDARD, CAT_NONE, "");
  addPar ("REACTION_FILE",string("NOFILE"),false,ParameterType::NO_DEP, &RxnSet::Model::rxnFileName, NULL, U_NONE, CAT_NONE, "");
  addPar ("NUMBER_REGIONS", 0, false, ParameterType::NO_DEP, &RxnSet::Model::userNumRegions, NULL, U_NONE, CAT_NONE, "Number of mesh points.");

  DeviceModel::initThermalModel(*this);

  addComposite("DOPINGPROFILES", DopeInfo::getParametricData(), &RxnSet::Model::regionMap);
  addComposite("REGION", DopeInfo::getParametricData(), &RxnSet::Model::regionMap);
  addComposite("SOURCELIST", SpecieSource::getParametricData(), &RxnSet::Model::defectSourceMap);
}

namespace RxnSet {



ParametricData<Instance> &Instance::getParametricData() {
  static ParametricData<Instance> parMap;

  return parMap;
}

ParametricData<Model> &Model::getParametricData() {
  static ParametricData<Model> parMap;

  return parMap;
}

//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/09/08
//-----------------------------------------------------------------------------
bool Instance::processParams (string param)
{
  updateTemperature(TEMP);
  return true;
}

//----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : instance block constructor
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter
// Creation Date : 02/09/08
//----------------------------------------------------------------------------

Instance::Instance(InstanceBlock & IB,
                   Model & it_MB,
                   MatrixLoadData & mlData1,
                   SolverState &ss1,
                   ExternData  &ed1,
                   DeviceOptions & do1)
  : DevicePDEInstance(IB,mlData1, ss1, ed1, do1),
    model_(it_MB),
    reactionFileCopyFlag(false),
    haveAnyReactions(false),
    useScaledVariablesFlag(true),
    useDopingArrayData(false),
    outputInterval(0.0),
    outputIndex(0),
    lastOutputTime(-10.0),
    outputRegion(0),
    tecplotLevel(0),
    callsOTEC(0),
    callsOTECcarrier(0),
    TEMP(300.0),
    newtonIterOld (0),

    li_Pos(-1),
    li_Neg(-1),

    outputXscalar(1.0),
    excludeNoSourceRegionsFlag(true),
    excludeNoSourceRegionsFlagGiven(false),
    transportFlagGiven(false),
    transportFlag(false),
    diffusionFlagGiven(false),
    diffusionFlag(false),
    dirichletBCFlag(false),
    columnReorderingFlag(false),
    xloIndex(-1),
    xhiIndex(-1)
{
  numIntVars   = 0;
  numExtVars   = 2;
  numStateVars = 10;

  devConMap.resize(2);
  devConMap[0] = 1;
  devConMap[1] = 1;


  // Set up mapping from param names to class variables:

  // if (parMap.empty())
  // {
  //   setNumNodes(2);
  //   setNumOptionalNodes(0);
  //   setNumFillNodes(0);
  //   setModelRequired(1);
  //   addModelType("PN");
  //   addModelType("NP");
  //   addModelType("RXN");

  //   addPar ("TEMP", 0.0, false, ParameterType::TIME_DEP,
  //           static_cast <double DeviceEntity:: *> (&Instance::TEMP),
  //           NULL, STANDARD, CAT_NONE, "");

  //   // user-specified scaling vars:
  //   addPar ("X0", 1.0e-7, false, ParameterType::NO_DEP,
  //           static_cast <double DeviceEntity:: *> (&Instance::x0_user),
  //           NULL, U_NONE, CAT_NONE, "Length Scalar.");

  //   addPar ("C0", 1.0e+12, false, ParameterType::NO_DEP,
  //           static_cast <double DeviceEntity:: *> (&Instance::C0_user),
  //           NULL, U_NONE, CAT_NONE, "Concentration Scalar.");

  //   addPar ("t0", 1.0e-6, false, ParameterType::NO_DEP,
  //           static_cast <double DeviceEntity:: *> (&Instance::t0_user),
  //           NULL, U_NONE, CAT_NONE, "Time Scalar.");

  //   addPar ("outputXscalar", 1.0, false, ParameterType::NO_DEP,
  //           static_cast <double DeviceEntity:: *> (&Instance::outputXscalar),
  //           NULL, U_NONE, CAT_NONE, "Scalar for X axis in tecplot file (default unit is cm)");

  //   addPar ("OUTPUTINTERVAL", 0.0, false, ParameterType::NO_DEP,
  //           static_cast <double DeviceEntity:: *> (&Instance::outputInterval),
  //           NULL, U_NONE, CAT_NONE, "Output Interval(sec)");

  //   // Set up map for non-double precision variables:

  //   addPar ("COPYRXN", false, false, ParameterType::NO_DEP,
  //           static_cast <bool DeviceEntity:: *> (&Instance::reactionFileCopyFlag), NULL,
  //           U_LOGIC, CAT_NONE, "Flag for processing chemistry file only once, then copied to other mesh points.");

  //   addPar ("SCALERXN", true, false, ParameterType::NO_DEP,
  //           static_cast <bool DeviceEntity:: *> (&Instance::useScaledVariablesFlag), NULL,
  //           U_LOGIC, CAT_NONE, "Flag for applying scaling to the reaction equations.");

  //   addPar ("DIFFUSION", false, false, ParameterType::NO_DEP,
  //           static_cast <bool DeviceEntity:: *> (&Instance::diffusionFlag),
  //           static_cast <bool DeviceEntity:: *> (&Instance::diffusionFlagGiven),
  //           U_LOGIC, CAT_NONE, "Flag for enabling lattice defect diffusion.");

  //   addPar ("TRANSPORT", false, false, ParameterType::NO_DEP,
  //           static_cast <bool DeviceEntity:: *> (&Instance::transportFlag),
  //           static_cast <bool DeviceEntity:: *> (&Instance::transportFlagGiven),
  //           U_LOGIC, CAT_NONE, "Flag for enabling lattice defect diffusion.  Identical to DIFFUSION flag, above.  Do not set both!");

  //   addPar ("EXCLUDENOSOURCE", true, false, ParameterType::NO_DEP,
  //           static_cast <bool DeviceEntity:: *> (&Instance::excludeNoSourceRegionsFlag),
  //           static_cast <bool DeviceEntity:: *> (&Instance::excludeNoSourceRegionsFlagGiven),
  //           U_LOGIC, CAT_NONE, "Flag for excluding regions that are outside of source region from computing defect reaction equations.  This is a speed optimization.  Turning it on will NOT change the answer");

  //   addPar ("COLUMNREORDER", false, false, ParameterType::NO_DEP,
  //           static_cast <bool DeviceEntity:: *> (&Instance::columnReorderingFlag), NULL,
  //           U_LOGIC, CAT_NONE, "Debug Flag for turning on/off column reordering.");

  //   addPar ("OUTPUTREGION", 1, false, ParameterType::NO_DEP,
  //           static_cast <int DeviceEntity:: *> (&Instance::outputRegion), NULL);

  //   addPar ("TECPLOTLEVEL", 0, false, ParameterType::NO_DEP,
  //           static_cast <int DeviceEntity:: *> (&Instance::tecplotLevel), NULL,
  //           U_NONE, CAT_NONE, "Integer number to determine type of tecplot output.  0=no output.  1=single time-dependent file, with each time step in a different zone.");

  //   addPar ("DIRICHLETBC", false, false, ParameterType::NO_DEP,
  //           static_cast <bool DeviceEntity:: *> (&Instance::dirichletBCFlag), NULL,
  //           U_LOGIC, CAT_NONE, "Flag for using Dirichlet boundary conditions.");
  // }

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  // Set any non-constant parameter defaults:
  if (!given("TEMP"))
    TEMP = getDeviceOptions().temp.dVal();


  if (diffusionFlagGiven && !transportFlagGiven)
  {
    transportFlag = diffusionFlag;
  }

  if (diffusionFlagGiven && transportFlagGiven) // both given
  {
    string msg = "Both transportFlag and diffusionFlag set in " + getName() + ".  Using transportFlag";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
  }

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  setupScalingVars ();

  // Handle the reactions.
  setupMeshUniform ();
  allocateRegions ();
  scaleMesh ();

  // Do the rest of the defect chemistry setup, which is independent of whether
  // or not the model or instance spefication has been used.
  initializeChemistry ();

  int numRegions=regVec.size();

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    vector<RegionData*> & rdVec1 = model_.regionDataVec;
    if (!(rdVec1.empty()))
    {
      cout << "Model Region Data vector:" << endl;
      for (int ireg=0;ireg<numRegions;++ireg)
      {
        RegionData & rd = *(rdVec1[ireg]);
        cout << ireg << ":  "<< rd;
      }
    }
  }
#endif

  // calculate dependent (ie computed) params:
  processParams ();
  setupFluxVec ();
  setupJacStamp ();
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupMeshUniform
// Purpose       : Sets up xVec and dxVec (unscaled).
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/27/08
//-----------------------------------------------------------------------------
void Instance::setupMeshUniform ()
{
  if ( model_.userNumRegions > 0)
  {
    int size = model_.userNumRegions;
    double xmin = model_.xlo;
    double xmax = model_.xhi;

    double dn = static_cast<double>(size-1);
    double dx = (xmax-xmin)/dn;
    double xtmp = xmin;

    xVec.resize(size,0.0);
    dxVec.resize(size,0.0);

    xloStencilVec.resize(size,0);
    xhiStencilVec.resize(size,0);

    int i=0;
    for (; i<size;++i)
    {
      xVec[i] = xtmp;
      xtmp += dx;
    }

    for (i=0;i<size-1;++i)
    {
      dxVec[i] = xVec[i+1]-xVec[i];
    }

    if (size-2 >= 0)
    {
      dxVec[size-1] = dxVec[size-2];
    }

    // set up the edge indices
    xloIndex=0;
    xhiIndex=size-1;
  }
  else
  {
    // no-op
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::allocateRegions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/03/08
//-----------------------------------------------------------------------------
void Instance::allocateRegions ()
{
  map <string, DopeInfo *> diMap = model_.dopeInfoMap;

  if ( model_.userNumRegions > 0 && model_.rxnFileName != "NOFILE" )
  {
    //setupDopingArrays();
    useDopingArrayData = true;

    // Now set up corresponding region Data classes, to match the region classes.
    // mesh stuff will be fixed later also:
    int numReg = model_.userNumRegions;

    vector<RegionData*> * rdVecPtr(0);
    rdVecPtr = &(model_.regionDataVec);
    // if (*rdVecPtr) is not empty, then the contents came from the input file and
    // should be deleted, as this specification overrides.
    (*rdVecPtr).clear();
    int iReg=0;
    for (iReg=0;iReg!=numReg;++iReg)
    {
      RegionData * regDataPtr = new RegionData ();
      double xtmp = xVec[iReg];

      vector<char> tmpname;
      tmpname.resize( getName().size()+10 );
      sprintf( &(tmpname[0]), "%s_%03d",getName().c_str(),iReg);
      string regName(&(tmpname[0]));

      vector<char> tmpname2;
      tmpname2.resize( outputName.size()+10 );
      sprintf( &(tmpname2[0]), "%s_%03d",outputName.c_str(),iReg);
      string outRegName(&(tmpname2[0]));

      regDataPtr->xloc = xtmp;
      regDataPtr->name = regName;
      regDataPtr->outName = outRegName;
      regDataPtr->reactionFile = model_.rxnFileName;

      (*rdVecPtr).push_back(regDataPtr);
    }

    if (reactionFileCopyFlag)
    {
      // Allocate a single reaction network class, and use it to parse the reactions
      // file.  This reaction network will then be copied into each reaction region
      // as it is constructed.
      N_DEV_ReactionNetwork tmpReactions1;
      N_DEV_ReactionNetwork tmpReactions2;

      tmpReactions1.setApplySources(true);
      tmpReactions1.setReactionNetworkFromFile(model_.rxnFileName);

      tmpReactions2.setApplySources(false);
      tmpReactions2.setReactionNetworkFromFile(model_.rxnFileName);

      // Now that the region data class have been set up, allocate the regions.
      for (iReg=0;iReg!=numReg;++iReg)
      {
        Region * regPtr;
        bool sourceOn(true);
        sourceOn = (xVec[iReg] >= model_.xlo_source) &&
                   (xVec[iReg] <= model_.xhi_source);

        if (sourceOn)
        {
          regPtr = new Region
                   ( *((*rdVecPtr)[iReg]),
                     devOptions,
                     solState,
                     tmpReactions1);
        }
        else
        {
          // if no source and transport is off, then no point in
          // solving reaction network at this mesh cell.
          if (transportFlag)
          {
            regPtr = new Region
                     ( *((*rdVecPtr)[iReg]),
                       devOptions,
                       solState,
                       tmpReactions2);
          }
          else
          {
            ((*rdVecPtr)[iReg])->doNothing=true;

            regPtr = new Region
                     ( *((*rdVecPtr)[iReg]),
                       devOptions,
                       solState,
                       tmpReactions2);
          }
        }
        regVec.push_back(regPtr);
      }
    }
    else
    {
      // Now that the region data class have been set up, allocate the regions.
      for (iReg=0;iReg!=numReg;++iReg)
      {
        bool sourceOn(true);
        sourceOn = (xVec[iReg] >= model_.xlo_source) &&
                   (xVec[iReg] <= model_.xhi_source);

        if (excludeNoSourceRegionsFlag)
        {
          if (!transportFlag && !sourceOn)
          {
            ((*rdVecPtr)[iReg])->doNothing=true;
          }
        }

        Region * regPtr = new Region
                          ( *((*rdVecPtr)[iReg]),
                            devOptions,
                            solState, sourceOn);

        regVec.push_back(regPtr);
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupScalingVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/04/08
//-----------------------------------------------------------------------------
void Instance::setupScalingVars ()
{
  scalingVars.t0 = 1.0e-6;          // time scaling (s)
  scalingVars.C0 = 1.0e+12;         // concentration scaling (cm^-3);
  scalingVars.x0  = 1.0e-7;         // distance scaling (cm)

  // time scaling (s)
  if (given("t0"))
  {
    scalingVars.t0 = t0_user;
  }

  // concentration scaling (cm^-3);
  if (given("C0"))
  {
    scalingVars.C0 = C0_user;
  }

  // distance scaling (cm)
  if (given("X0"))
  {
    scalingVars.x0  = x0_user;
  }

  double rx0 = 1.0/scalingVars.x0;

  // area scaling (cm^2)
  scalingVars.a0 = scalingVars.x0*scalingVars.x0;

  // diffusion coefficient scaling (cm^2/s)
  //scalingVars.D0  = scalingVars.t0*rx0*rx0;
  scalingVars.D0  = scalingVars.a0/scalingVars.t0;

  // recombination rate scaling (cm^-3/s)
  //scalingVars.R0  = scalingVars.D0*scalingVars.C0/(scalingVars.x0*scalingVars.x0);
  scalingVars.R0  = scalingVars.C0/scalingVars.t0;

  scalingVars.rR0 = 1.0/scalingVars.R0;

  // rate constant scaling.   k0 = 1/(C0*t0) = cm^3/sec
  scalingVars.rk0 = scalingVars.C0*scalingVars.t0;
  scalingVars.rt0 = 1.0/scalingVars.t0;
  scalingVars.k0 = 1.0/scalingVars.rk0;
}

//-----------------------------------------------------------------------------
// Function      : Instance::initializeChemistry
//
// Purpose       : This function takes care of some of the initial setup for
//                 defect chemistry.  It is independent of the model or
//                 instance specification for the chemistry, and thus gets
//                 its own function.  The other functions are
//                 setupUserSpecifiedRegions (for instance based).
//
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/02/08
//-----------------------------------------------------------------------------
void Instance::initializeChemistry ()
{
#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    // debug outputs.  Check what is in the defect map.
    map <string,CompositeParam *>::iterator itSource = model_.defectSourceMap.begin();
    map <string,CompositeParam *>::iterator itSourceEnd = model_.defectSourceMap.end ();

    cout << "----------------------------------------------" << endl;
    cout << "Instance::initializeChemistry ():" << endl;
    for (; itSource!=itSourceEnd; ++itSource)
    {
      string speciesName (itSource->first);
      cout << "speciesName = " << speciesName << endl;
    }
    cout << "----------------------------------------------" << endl;
  }
#endif

  vector<RegionData*> * rdVecPtr(0);
  rdVecPtr = &(model_.regionDataVec);

  int numRegions=regVec.size();
  for (int ireg=0;ireg<numRegions;++ireg)
  {
    regVec[ireg]->initializeReactionNetwork(scalingVars);

    if (useDopingArrayData)
    {
      map <string, DopeInfo *> diMap = model_.dopeInfoMap;
      if (!(diMap.empty()))
      {
        // Push in some initial values for various species.
        // This includes doping species like boron- (BM) as well as E and H.
        map<string, DopeInfo *>::iterator iter;
        map<string, DopeInfo *>::iterator start = diMap.begin();
        map<string, DopeInfo *>::iterator end   = diMap.end();

        for ( iter = start; iter != end; ++iter )
        {
          DopeInfo & di = *(iter->second);

          bool reactantExist = false;
          reactantExist = regVec[ireg]->reactantExist(di.speciesName);

          if (reactantExist)
          {
            double dopeValue = di.splintDopeVec[ireg];
            regVec[ireg]->setInitialCondition(di.speciesName, dopeValue);
          }
        }
      }
      else
      {
#if 0
        // This is the original doping/carrier intialization.
        // It is hardwired to boron and phosphorus, so
        // it should be phased out and/or deprecated.

        if (regVec[ireg]->reactantExist("BM"))
        {
          regVec[ireg]->setInitialCondition("BM", (*rdVecPtr)[ireg]->Boron_Concentration);
        }

        if (regVec[ireg]->reactantExist("PP"))
        {
          regVec[ireg]->setInitialCondition("PP", (*rdVecPtr)[ireg]->Phosphorus_Concentration);
        }
#endif
      }

      if (model_.given("MASTERSOURCE"))
      {
        bool sourceOn(true);
        double xlos = model_.xlo_source *(useScaledVariablesFlag?(1.0/scalingVars.x0):1.0);
        double xhis = model_.xhi_source *(useScaledVariablesFlag?(1.0/scalingVars.x0):1.0);

        sourceOn = (xVec[ireg] >= xlos) && (xVec[ireg] <= xhis);

        if (sourceOn)
        {
          map <string,CompositeParam *>::iterator itSource = model_.defectSourceMap.begin();
          map <string,CompositeParam *>::iterator itSourceEnd = model_.defectSourceMap.end ();

          for (; itSource!=itSourceEnd; ++itSource)
          {
            string speciesName (itSource->first);
            regVec[ireg]->addMasterSource(speciesName);
          }
        }
      }

      if (useScaledVariablesFlag)
      {
        regVec[ireg]->scaleVariables ();
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupJacStamp
// Purpose       :
// Special Notes :
//
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/09/08
//-----------------------------------------------------------------------------
void Instance::setupJacStamp ()
{
  // Now it's time to do the jacobian stamp.
  // Set up full jacstamp
  int localPosIndex = 0;
  int localNegIndex = 1;
  int indexCounter=localNegIndex+1;

  jacStamp.resize(2);

  jacStamp[localPosIndex].resize(2);  // Base'
  jacStamp[localPosIndex][localPosIndex] = localPosIndex;     // B'-B' (0)
  jacStamp[localPosIndex][localNegIndex] = localNegIndex;     // B'-E' (1)

  jacStamp[localNegIndex].resize(2);  // Emitter'
  jacStamp[localNegIndex][localPosIndex] = localPosIndex;     // E'-B' (0)
  jacStamp[localNegIndex][localNegIndex] = localNegIndex;     // E'-E' (1)

  // Rxn model blocks are next.
  regLastIndexVec.resize(regVec.size(), -1); // index at the end of the stamp.
  regFirstReactantIndexVec.resize(regVec.size(), -1); // index of first reactant
  regNumSpecieVec.resize(regVec.size(), 0);

  int ireg=0;

  // Call each region to have it augment its share of the jacstamp.
  int numRegions=regVec.size();
  for (ireg=0;ireg<numRegions;++ireg)
  {
    int concentrationSize = regVec[ireg]->getNumSpecies();
    if (concentrationSize != 0)
    {
      // communicate down to the region the voltage inputs indices.
      vector<int> voltageNodeColDep(2,-1);
      voltageNodeColDep[0]=localPosIndex; // pos (previously base) voltage
      voltageNodeColDep[1]=localNegIndex; // neg (perviously emitter) voltage

      int firstReactant=-1;
      int lastIndex=-1;
      regVec[ireg]->setupJacStamp(jacStamp,voltageNodeColDep, firstReactant, lastIndex);
      regFirstReactantIndexVec[ireg] = firstReactant;
      regLastIndexVec[ireg] = lastIndex;

      numIntVars += regVec[ireg]->getNumIntVars();
      numStateVars += concentrationSize;
    }
  }

  // Do diffusion. for all mobile species.
  int numSpecies = thVec.size();
  int isp=0;

  for (ireg=0;ireg<numRegions;++ireg)
  {
    if (!(regVec[ireg]->getDoNothingFlag()) )
    {
      regNumSpecieVec[ireg] = numSpecies;
    }
    else
    {
      regNumSpecieVec[ireg] = 0;
    }
#ifdef Xyce_DEBUG_DEVICE

    cout << regVec[ireg]->getName ()
         << "   numSpecies = " << regNumSpecieVec[ireg] << endl;
#endif
  }

  if (transportFlag && numRegions > 1)
  {
    for (isp=0;isp<numSpecies;++isp)
    {
      if (!(thVec[isp].transportFlag)) continue;

      // If this is the base or emitter, then add a single column to the stamp.
      // Otherwise, if this is the BE point (in the middle) then add 2 columns,
      // for 2 neighbors.
      for (ireg=0;ireg<numRegions;++ireg)
      {
        int row = regFirstReactantIndexVec[ireg] + isp;
        if (row < 0) continue;

        int rowSize = jacStamp[row].size();
        if (ireg==0)
        {
          int col = regFirstReactantIndexVec[1] + isp;
          jacStamp[row].resize(rowSize+1);
          jacStamp[row][rowSize] = col;
        }
        else if (ireg==numRegions-1)
        {
          int col = regFirstReactantIndexVec[numRegions-2] + isp;
          jacStamp[row].resize(rowSize+1);
          jacStamp[row][rowSize] = col;
        }
        else
        {
          int col1 = regFirstReactantIndexVec[ireg+1] + isp;
          int col2 = regFirstReactantIndexVec[ireg-1] + isp;
          jacStamp[row].resize(rowSize+2);
          jacStamp[row][rowSize  ] = col1;
          jacStamp[row][rowSize+1] = col2;
        }
      }
    }
  }

  // Just put in a dense-ish row for these terms for now.
  // These represent the dependence of generation-recombination
  // currents, calculated from regions, on the species
  // densities of those same regions.
  {
    //
    // do the base-prime row and the emit-prime row.
    int posRow = 0;
    int posSize = jacStamp[posRow].size();
    int negRow = 1;
    int negSize = jacStamp[negRow].size();
    for (ireg=0;ireg<numRegions;++ireg)
    {
      int posIndex = regFirstReactantIndexVec[ireg];
      int numSpecie = regNumSpecieVec[ireg];

      if (posIndex != -1)
      {
        int posRowSize = jacStamp[posRow].size();
        jacStamp[posRow].resize(posRowSize+numSpecie);

        int negRowSize = jacStamp[negRow].size();
        jacStamp[negRow].resize(negRowSize+numSpecie);

        for (isp=0;isp<numSpecie;++isp)
        {
          int icol = posRowSize+isp;
          jacStamp[posRow][icol] = posIndex+isp ;

          icol = negRowSize+isp;
          jacStamp[negRow][icol] = posIndex+isp ;
        }
      }
    }
  }

  // set up normal jacMap for when all resistances nonzero
  // If nothing is remapped, this amounts to a null operation when the
  // map is used later.  The maps become important when we start
  // remapping nodes because of zero lead resistances
  jacMap.clear();
  jacMap2.clear();
  jacMap.resize(jacStamp.size());
  jacMap2.resize(jacStamp.size());

  // This is only done once per model, and size is small --- not gonna
  // bother with temporary vars just to keep the .size() call out of the
  // loops

  int mapSize = jacMap.size();
  for (int i=0;i<mapSize;++i)
  {
    jacMap[i]=i;
    jacMap2[i].resize(jacStamp[i].size());
    for (int j=0;j<jacStamp[i].size();++j)
    {
      jacMap2[i][j] = j;
    }
  }

  // Now fix the ordering of the columns in the jacStamp.  If the columns in each row
  // are not in ascending order, then the jacStampMap calls below (for removing
  // variables) will not work correctly.
  //
  // ERK: 01/13/10:  this flag doesn't appear to be needed for this model,
  // as this model doesn't contain any optional prime nodes.  This block
  // of code was copied over from another model, which *does* need it.
  // As such, the columnReorderFlag is "false" by default, whereas it is
  // "true" in the orginal device from which this was copied.
  if (columnReorderingFlag)
  {
    vector< vector<int> > tempStamp_eric;
    vector< vector<int> > tempMap2_eric;
    jacStampMap_fixOrder(jacStamp, jacMap2, tempStamp_eric, tempMap2_eric);
    jacStamp = tempStamp_eric;
    jacMap2 = tempMap2_eric;
  }

  // Now we have to selectively map away bits and pieces of the
  // jacobian stamp based on absent resistances
#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > -50 && getSolverState().debugTimeFlag)
  {
    cout << "jacStamp Before removing terminal nodes:"<<endl;
    outputJacStamp(jacStamp);
    cout << "jacMap2  Before removing terminal nodes:"<<endl;
    outputJacStamp(jacMap2 );
  }
#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::scaleMesh
// Purpose       :
// Special Notes : This should be called after the doping-profile based
//                 regions are set up.
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/27/08
//-----------------------------------------------------------------------------
void Instance::scaleMesh ()
{
  vector<RegionData*> * rdVecPtr(0);
  rdVecPtr = &(model_.regionDataVec);
  int i=0;
  int size = (*rdVecPtr).size();
  if (useScaledVariablesFlag)
  {
    for (i=0;i<size;++i)
    {
      xVec[i] *= (1.0/scalingVars.x0);
    }
  }

  for (i=0;i<size-1;++i)
  {
    dxVec[i] = xVec[i+1]-xVec[i];
  }

  if (size-2 >= 0)
  {
    dxVec[size-1] = dxVec[size-2];
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 1)
  {
    for (i=0;i<size;++i)
    {
      cout << "Scaled Mesh:  xVec["<<i<<"] = " << xVec[i] << endl;
    }
    cout << endl;
    for (i=0;i<size-1;++i)
    {
      cout << "Scaled Mesh:  dxVec["<<i<<"] = " << dxVec[i] << endl;
    }
    cout << endl;
  }
#endif

}

//-----------------------------------------------------------------------------
// Function      : Instance::setupFluxVec
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 10/27/08
//-----------------------------------------------------------------------------
void Instance::setupFluxVec ()
{
  vector<RegionData*> * rdVecPtr(0);
  rdVecPtr = &(model_.regionDataVec);
  if ( !((*rdVecPtr).empty()) )
  {
    // do all the mobile species
    int numSpecies = thVec.size();
    int isp=0;
    int size = (*rdVecPtr).size();
    for (isp=0;isp<numSpecies;++isp)
    {
      thVec[isp].fluxVec.resize(size-1,0.0);
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::outputPlotFiles
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/09/08
//-----------------------------------------------------------------------------
bool Instance::outputPlotFiles ()
{
  bool bsuccess = true;
  bool bs1 = true;
  bool skipOutput = false;

  // If using output interval, check if enough time has passed to do
  // another output.  (only applies for transient - not DCOP).
  if ( !(getSolverState().dcopFlag) && !(getSolverState().forceFinalOutput) && given("OUTPUTINTERVAL") )
  {
    double outMult = static_cast<double> (outputIndex);
    double nextOutputTime = outMult * outputInterval;

    if (nextOutputTime > getSolverState().currTime)
    {
      skipOutput = true;
    }
  }

  // If this is a "forced" final output, make sure that it didn't already output.
  // This can happen if the output interval is an exact multiple of the
  // total simulation time.
  if (getSolverState().forceFinalOutput &&
      getSolverState().currTime==lastOutputTime) skipOutput=true;

  if (skipOutput) return bsuccess;
  ++outputIndex;
  lastOutputTime = getSolverState().currTime;


  // output individual tecplot files for each region.
  // These are time-dependent files, with no spatial depedence.
  int numRegions = regVec.size();
  if ( (tecplotLevel==2 && numRegions > 1) || (tecplotLevel==1 && numRegions==1) )
  {
    for (int ireg=0;ireg<numRegions;++ireg)
    {
      if (regVec[ireg]->haveAnyReactions())
      {
        bs1 = regVec[ireg]->outputTecplot ();
        bsuccess = bsuccess && bs1;
      }
    }
  }

  // output a single, spatially-dependent file(s).
  // Don't bother if there's only one (or less) regions.
  if (tecplotLevel==1 && numRegions > 1)
  {
    bs1 = outputTecplot ();
    bsuccess = bsuccess && bs1;
  }
  else if (tecplotLevel==3 && numRegions > 1)
  {
    bs1 = output2DTecplot ();
    bsuccess = bsuccess && bs1;
  }

  if (tecplotLevel > 0 && given("OUTPUTREGION"))
  {
    int numRegions = regVec.size();
    if (outputRegion < numRegions && outputRegion >= 0)
    {
      if (regVec[outputRegion]->haveAnyReactions())
      {
        bs1 = regVec[outputRegion]->outputTecplot ();
        bsuccess = bsuccess && bs1;
      }
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::outputTecplot
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/09/08
//-----------------------------------------------------------------------------
bool Instance::outputTecplot ()
{
  bool bsuccess = true;

  int i(0);
  int NX = regVec.size();

  char filename[256];
  for(i=0;i<256;++i) filename[i] = static_cast<char>(0);

  sprintf(filename,"%s.dat",outputName.c_str());
  double time = getSolverState().currTime;
  FILE *fp1(NULL);

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "  Instance::outputTecplot.  filename = " << string(filename) <<endl;
  }
#endif

  if (callsOTEC <= 0)
  {
    fp1 = fopen(filename,"w");
    fprintf(fp1,
            " TITLE = \"Spatially Dependent defect data for compact rxn device: %s  time = %20.12e seconds.\",\n",
            outputName.c_str(),time);
  }
  else
  {
    fp1 = fopen(filename,"a");
  }

  int rSize=(regVec[0])->getNumSpecies();
  int cSize=(regVec[0])->getNumConstants();

  if (callsOTEC <= 0)
  {
    fprintf(fp1,"%s","\tVARIABLES = \"X \",\n");

    Region & reg = (*(regVec[0]));
    for (int iconst=0;iconst<cSize;++iconst)
    {
      fprintf(fp1, "\t    \"%s \",\n" , ( reg.getConstantsName(iconst)).c_str());
    }

    for (int ispec=0;ispec<rSize;++ispec)
    {
      fprintf(fp1,"\t    \"%s \",\n", ( reg.getSpeciesName(ispec)).c_str());
    }
  }

  fprintf(fp1,"\tDATASETAUXDATA %s\n", tecplotTimeDateStamp().c_str() );
  fprintf(fp1,"\tZONE F=POINT,I=%d", NX);

  if (getSolverState().dcopFlag)
  {
    fprintf(fp1,"  T = \"DCOP step = %d\" \n", callsOTEC);
  }
  else
  {
    fprintf(fp1,"  T = \"time step = %d, time=%20.12e\" AUXDATA time = \"%20.12e seconds\" \n", callsOTEC , time, time);
  }

  vector<RegionData*> * rdVecPtr(0);
  rdVecPtr = &(model_.regionDataVec);

  if (NX <= 1)
  {
    fprintf(fp1,"  %20.12e", time);
  }

  int numprint=0;
  for (i=0;i<NX;++i)
  {
    double val(0.0);
    val = ((*rdVecPtr)[i])->xloc * outputXscalar;
    fprintf(fp1,"  %20.12e", val);

    Region & reg = (*(regVec[i]));

    for (int iconst=0;iconst<cSize;++iconst)
    {
      val = reg.getConstantsVal(iconst) ;
      fprintf(fp1,"  %20.12e", val);

      if (numprint >= 6)
      {
        fprintf(fp1,"%s","\n"); numprint = 0;
      }
      else
      {
        numprint++;
      }
    }

    for (int ispec=0;ispec<rSize;++ispec)
    {
      val = 0.0;
      if (reg.haveAnyReactions())
      {
        val = reg.getSpeciesVal(ispec);
      }
      fprintf(fp1,"  %20.12e", val);

      if (numprint >= 6)
      {
        fprintf(fp1,"%s","\n"); numprint = 0;
      }
      else
      {
        numprint++;
      }
    }
  }

  fprintf(fp1,"%s","\n");

  ++callsOTEC;
  fclose(fp1);

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::output2DTecplot
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 10/14/09
//-----------------------------------------------------------------------------
bool Instance::output2DTecplot ()
{
  bool bsuccess = true;

  int i(0);
  int NX = regVec.size();

  char filename[256];
  for(i=0;i<256;++i) filename[i] = static_cast<char>(0);

  sprintf(filename,"%s.dat",outputName.c_str());
  double time = getSolverState().currTime;
  FILE *fp1(NULL);

  if (callsOTEC <= 0)
  {
    fp1 = fopen(filename,"w");
    fprintf(fp1,
            " TITLE = \"Spatially Dependent defect data for compact rxn device: %s  time = %20.12e seconds.\",\n",
            outputName.c_str(),time);
  }
  else
  {
    fp1 = fopen(filename,"a");
  }

  int rSize=(regVec[0])->getNumSpecies();
  int cSize=(regVec[0])->getNumConstants();

  if (callsOTEC <= 0)
  {
    fprintf(fp1,"%s","\tVARIABLES = \"X \",\n");
    fprintf(fp1,"%s","\t    \"Time(sec) \",\n");

    Region & reg = (*(regVec[0]));
    for (int iconst=0;iconst<cSize;++iconst)
    {
      fprintf(fp1, "\t    \"%s \",\n" , ( reg.getConstantsName(iconst)).c_str());
    }

    for (int ispec=0;ispec<rSize;++ispec)
    {
      fprintf(fp1,"\t    \"%s \",\n", ( reg.getSpeciesName(ispec)).c_str());
    }

    fprintf(fp1,"\tZONE F=POINT,I=%d", NX);
  }

  fprintf(fp1,"%s"," \n");

  vector<RegionData*> * rdVecPtr(0);
  rdVecPtr = &(model_.regionDataVec);

  if (NX <= 1)
  {
    fprintf(fp1,"  %20.12e", time);
  }

  int numprint=0;
  for (i=0;i<NX;++i)
  {
    double val(0.0);
    val = ((*rdVecPtr)[i])->xloc;
    fprintf(fp1,"  %20.12e", val);

    fprintf(fp1,"  %20.12e", time);

    Region & reg = (*(regVec[i]));

    for (int iconst=0;iconst<cSize;++iconst)
    {
      val = reg.getConstantsVal(iconst) ;
      fprintf(fp1,"  %20.12e", val);

      if (numprint >= 6)
      {
        fprintf(fp1,"%s","\n"); numprint = 0;
      }
      else
      {
        numprint++;
      }
    }

    for (int ispec=0;ispec<rSize;++ispec)
    {
      val = reg.getSpeciesVal(ispec);
      fprintf(fp1,"  %20.12e", val);

      if (numprint >= 6)
      {
        fprintf(fp1,"%s","\n"); numprint = 0;
      }
      else
      {
        numprint++;
      }
    }
  }

  fprintf(fp1,"%s","\n");

  ++callsOTEC;
  fclose(fp1);

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::outputCarrierDensities
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/09/08
//-----------------------------------------------------------------------------
bool Instance::outputCarrierDensities ()
{
  bool bsuccess = true;

  int i(0);
  int NX = regVec.size();
  char filename[256];   for(i=0;i<256;++i) filename[i] = static_cast<char>(0);
  sprintf(filename,"%scarrier.dat",outputName.c_str());

  FILE *fp1(NULL);
  fp1 = fopen(filename,"w");
  int cSize=(regVec[0])->getNumConstants();

  vector<RegionData*> * rdVecPtr(0);
  rdVecPtr = &(model_.regionDataVec);

  for (i=0;i<NX;++i)
  {
    double val(0.0);
    val = ((*rdVecPtr)[i])->xloc;
    fprintf(fp1,"  %20.12e", val);

    Region & reg = (*(regVec[i]));

    for (int iconst=0;iconst<cSize;++iconst)
    {
      val = reg.getConstantsVal(iconst) ;
      fprintf(fp1,"  %20.12e", val);
    }
    fprintf(fp1,"%s","\n");
  }

  ++callsOTECcarrier;
  fclose(fp1);

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/09/08
//-----------------------------------------------------------------------------
Instance::~Instance()
{
  if (!(regVec.empty()))
  {
    int size = regVec.size();
    int i=0;
    for (i=0;i<size;++i)
    {
      if (regVec[i] != 0)
      {
        delete regVec[i];
        regVec[i] = 0;
      }
    }
  }

  // Loop over the dopeInfoMap (if it is not empty) and delete its contents.
  if (!(dopeInfoMap.empty()))
  {
    map <string,DopeInfo *>::iterator iter;
    map <string,DopeInfo *>::iterator begin = dopeInfoMap.begin();
    map <string,DopeInfo *>::iterator end   = dopeInfoMap.end  ();

    for(iter=begin;iter!=end;++iter)
    {
      if (iter->second != 0) delete iter->second;
    }
  }

  regVec.clear();
  dopeInfoMap.clear();
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 02/09/08
//-----------------------------------------------------------------------------
void Instance::registerLIDs( const vector<int> & intLIDVecRef,
                             const vector<int> & extLIDVecRef )
{
  string msg;

  // Check if the size of the ID lists corresponds to the proper number of
  // internal and external variables
  int numInt = intLIDVecRef.size();
  int numExt = extLIDVecRef.size();

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline =
    "----------------------------------------------------------------------------";

  if (getDeviceOptions().debugLevel > 0)
  {
    cout << dashedline << endl;
    cout << "  Instance::registerLIDs" <<endl;
    cout << "  name = " << getName() << endl;
    cout << "  number of internal variables: " << numInt << endl;
    cout << "  number of external variables: " << numExt << endl;
    cout << "  numIntVars = " << numIntVars << endl;
    cout << "  numExtVars = " << numExtVars << endl;
  }
#endif

  if (numInt != numIntVars)
  {
    msg = "Instance::registerLIDs:";
    msg += "numInt != numIntVars";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
  }

  if (numExt != numExtVars)
  {
    msg = "Instance::registerLIDs:";
    msg += "numExt != numExtVars";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
  }

  // Copy over the global ID lists.
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "  Internal LID List" << endl;
    for( int i = 0; i < intLIDVec.size(); ++i )
      cout << "  " << intLIDVec[i] << endl;
    cout << "  External LID List" << endl;
    for( int i = 0; i < extLIDVec.size(); ++i )
      cout << "  " << extLIDVec[i] << endl;
  }
#endif

  // Use these lists to obtain the indices into the linear algebra entities.
  // First do external variables:
  int extIndex = 0;
  int intIndex = 0;

  li_Pos = extLIDVec[extIndex++];
  li_Neg = extLIDVec[extIndex++];

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "  li_Pos = " << li_Pos << endl;
    cout << "  li_Neg = " << li_Neg << endl;
  }
#endif

  // Register the LIDs for each reaction region.
  int numRegions = regVec.size();
  for (int ireg=0;ireg<numRegions;++ireg)
  {
    regVec[ireg]->registerLIDs(intLIDVec, extLIDVec,intIndex);
  }

  // Do diffusion. for all mobile species.
  int i=0;
  int ispec=0;
  int numSpecies = thVec.size();
  int size = regVec.size();

  if (transportFlag)
  {
    for (ispec=0;ispec<numSpecies;++ispec)
    {
      string speciesName = regVec[0]->getSpeciesName(ispec);
      thVec[ispec].specie_id.resize(size,-1);
      for (i=0;i<size;++i)
      {
        thVec[ispec].specie_id[i] = (regVec[i])->getSpeciesLID (speciesName);
      }
    }
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
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
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/09/08
//-----------------------------------------------------------------------------
map<int,string> & Instance::getIntNameMap ()
{
  // set up the internal name map, if it hasn't been already.
  if (intNameMap.empty ())
  {
    // If we have internal nodes, then we need to do this for debugging
    // set up the internal name map:
    int numRegions = regVec.size();
    for (int ireg=0;ireg<numRegions;++ireg)
    {
      regVec[ireg]->augmentNameMap(intNameMap, *this);
    }
  }
  return intNameMap;
}


//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 02/09/08
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs( const vector<int> & staLIDVecRef )
{
  string msg;

  // Check if the size of the ID list corresponds to the proper number of state
  // variables.
  int numSta = staLIDVecRef.size();

  if (numSta != numStateVars)
  {
    msg = "Instance::registerStateLIDs:";
    msg += "numSta != numStateVars";
    N_ERH_ErrorMgr::report (N_ERH_ErrorMgr::DEV_FATAL, msg);
  }

  // Copy over the global ID lists:
  staLIDVec = staLIDVecRef;

  int i=0;

  int numRegions = regVec.size();
  for (int ireg=0;ireg<numRegions;++ireg)
  {
    regVec[ireg]->registerStateLIDs(staLIDVec,i);
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
// Creator       : Eric R. Keiter, SNL
// Creation Date : 02/09/08
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
// Creator       : Eric R. Keiter, SNL
// Creation Date : 02/09/08
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const vector< vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );

  vector<int> &map=jacMap;
  vector< vector<int> > &map2=jacMap2;
  int posI = 0;
  int negI = 1;

  APosEquPosNodeOffset  = jacLIDVec[map[posI]][map2[posI][posI]];
  APosEquNegNodeOffset  = jacLIDVec[map[posI]][map2[posI][negI]];

  ANegEquPosNodeOffset  = jacLIDVec[map[negI]][map2[negI][posI]];
  ANegEquNegNodeOffset  = jacLIDVec[map[negI]][map2[negI][negI]];

  // Do the reaction species block
  int numRegions = regVec.size();
  for (int ireg=0;ireg<numRegions;++ireg)
  {
    regVec[ireg]->registerJacLIDs(jacLIDVec, jacMap, jacMap2);
  }


#if 0
  // For now diffusion can only be set up in the jacobian matrix
  // via setRow function calls.

  // If this is the base or emitter, then add a single column to the stamp.
  // Otherwise, if this is the BE point (in the middle) then add 2 columns,
  // for 2 neighbors.
  for (int ireg=0;ireg<numRegions;++ireg)
  {
    int row = regPosIndexVec[ireg] + regV0subIndexVec[ireg];
    if (row < 0) continue;

    int rowSize = jacStamp[row].size();
    if (ireg==0)
    {
      int col = regPosIndexVec[1] + regV0subIndexVec[1];
      jacStamp[row].resize(rowSize+1);
      jacStamp[row][rowSize] = col;
    }
    else if (ireg==numRegions-1)
    {
      int col = regPosIndexVec[numRegions-2] + regV0subIndexVec[numRegions-2];
      jacStamp[row].resize(rowSize+1);
      jacStamp[row][rowSize] = col;
    }
    else
    {
      int col1 = regPosIndexVec[ireg+1] + regV0subIndexVec[ireg+1];
      int col2 = regPosIndexVec[ireg-1] + regV0subIndexVec[ireg-1];
      jacStamp[row].resize(rowSize+2);
      jacStamp[row][rowSize+1] = col1;
      jacStamp[row][rowSize+2] = col2;
    }
  }
#endif

  for (int iReg=0;iReg<numRegions;++iReg)
  {
    if (regVec[iReg]->haveAnyReactions())
    {
      haveAnyReactions=true;
    }
  }

}

//-----------------------------------------------------------------------------
// Function      : Instance::setupPointers()
//
// Purpose       : Sets up raw pointers for optimized matrix loads.
//
// Special Notes : ERK: This function is primarily concerned with determining
//                 matrix pointers associated with coupling the reaction model
//                 to the terminal nodes.  I've had a lot of difficulty
//                 getting the LID offsets correct for those terms.  To preserve
//                 load speed, I've decided to go directly to raw pointers
//                 for those terms.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 05/18/10
//-----------------------------------------------------------------------------
void Instance::setupPointers()
{
  int numRegions = regVec.size();

  N_LAS_Matrix & dFdxMat = *(extData.dFdxMatrixPtr);

  APosEqu_SpeciesPtr.resize(numRegions);
  ANegEqu_SpeciesPtr.resize(numRegions);

  APosEqu_ConstPtr.resize(numRegions);
  ANegEqu_ConstPtr.resize(numRegions);

  for (int ireg=0;ireg<numRegions;++ireg)
  {
    N_LAS_Matrix & dFdx = *(extData.dFdxMatrixPtr);
    N_LAS_Matrix & dQdx = *(extData.dQdxMatrixPtr);
    regVec[ireg]->setupPointers (dFdx,dQdx);
  }

}

//-----------------------------------------------------------------------------
// Function      : Instance::updateTemperature
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 02/09/08
//-----------------------------------------------------------------------------
bool Instance::updateTemperature( const double & temp )
{

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << "Start Instance::updateTemperature" << endl;
    cout << "temp = "<<temp << endl;
  }
#endif
  if( temp != -999.0 ) TEMP = temp;

  if(model_.interpolateTNOM(TEMP) )
  {
    // some params may have changed during interpolation
    model_.processParams();
  }

  //Generation of temperature based factors
  double TNOM  = model_.TNOM;

  // update the rate constants for the defect reactions
  int numRegions = regVec.size();
  for (int i=0;i<numRegions;++i)
  {
    regVec[i]->setRateConstants(TEMP);
  }

  // all Species diffusion:
  if (thVec.empty() && !(regVec.empty()) )
  {
    int numSpecies = regVec[0]->getNumSpecies();
    thVec.reserve(numSpecies);
    for (int ispec=0;ispec<numSpecies;++ispec)
    {
      string speciesName = regVec[0]->getSpeciesName(ispec);

      double Dtmp = (regVec[0])->getDiffusionCoefficient (ispec,TEMP);

      thVec.push_back(N_DEV_TransportHelper(Dtmp,speciesName));

#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 0)
      {
        if (useScaledVariablesFlag)
        {
          cout << "Vacancy Diffusion: scaling=TRUE ";
          cout << "D_"<<speciesName<<" = " << Dtmp << endl;
        }
        else
        {
          cout << "Vacancy Diffusion: scaling=FALSE ";
          cout << "D_"<<speciesName<<" = " << Dtmp << endl;
        }
      }
#endif
    }

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0)
    {
      for (int ispec=0;ispec<numSpecies;++ispec)
      {
        cout << "Vacancy Diffusion: D_" << thVec[ispec].name << "  transportFlag = ";
        if (thVec[ispec].transportFlag) cout << "TRUE";
        else cout << "FALSE";
        cout << endl;
      }
    }
#endif

  }

  return true;
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
// Creator       : Eric R. Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 02/09/08
//-----------------------------------------------------------------------------
bool Instance::loadDeviceMask ()
{
  bool returnVal=false;
  N_LAS_Vector & maskVector = *(extData.deviceMaskVectorPtr);

  int numRegions = regVec.size();
  if (numRegions > 0)
  {
    for (int ireg=0;ireg<numRegions;++ireg)
    {
      bool tmpval = regVec[ireg]->loadDeviceMask (maskVector);

      // if any of these return true, then set return val to true.
      if (tmpval) returnVal = true;
    }
  }

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
// Creator       : Eric R. Keiter, SNL
// Creation Date : 02/09/08
//-----------------------------------------------------------------------------
bool Instance::loadDAEQVector ()
{
  double * daeQVec = extData.daeQVectorRawPtr;
  double * dQdxdVp = extData.dQdxdVpVectorRawPtr;

  // Do the reaction terms.  There will be no voltage limiting
  // contributions here.
  int numRegions = regVec.size();
  for (int ireg=0;ireg<numRegions;++ireg)
  {
    regVec[ireg]->loadDAEQVector (daeQVec);
  }

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
// Creator       : Eric R. Keiter, SNL
// Creation Date : 02/09/08
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  double * daeFVec = extData.daeFVectorRawPtr;
  double * dFdxdVp = extData.dFdxdVpVectorRawPtr;
  double * solVector = extData.nextSolVectorRawPtr;

  double vbe_diff = 0.0;

  // Now do the reaction terms.  There will be no voltage limiting
  // contributions here.
  int numRegions = regVec.size();
  for (int ireg=0;ireg<numRegions;++ireg)
  {
    // Internal reaction terms are handled by the region class.
    regVec[ireg]->loadDAEFVector (daeFVec);
    regVec[ireg]->loadDAEdFdxdV(dFdxdVp, vbe_diff);
  }

  if (transportFlag)
  {
    if (!getSolverState().dcopFlag)
    {
      int numSpecies = thVec.size();
      int isp=0;
      for (isp=0;isp<numSpecies;++isp)
      {
        if (!(thVec[isp].transportFlag)) continue;

        int i=0;
        int size = regVec.size();

        vector<int> & specie_id = thVec[isp].specie_id;

        double dx1 = dxVec[0];
        double xloVal = (thVec[isp].fluxVec[0]-thVec[isp].flux_bc1)/dx1;
        daeFVec[specie_id[0]] += xloVal;

        for (i=1;i<size-1;++i)
        {
          double fluxDif = (thVec[isp].fluxVec[i]-thVec[isp].fluxVec[i-1]);
          double aveDx = (dxVec[i-1]+dxVec[i])*0.5;
          daeFVec[specie_id[i]] += (fluxDif)/aveDx;
        }

        double dx2 = dxVec[size-2];
        double xhiVal = (thVec[isp].flux_bc2-thVec[isp].fluxVec[size-2])/dx2;
        daeFVec[specie_id[size-1]] += xhiVal;
      } // end of species loop.
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx
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
// Creator       : Eric R. Keiter, SNL
// Creation Date : 02/09/08
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  bool bsuccess = true;
  bool tmpBool = true;
  N_LAS_Matrix & dQdxMat = *(extData.dQdxMatrixPtr);

#ifdef Xyce_DEBUG_DEVICE
  char tmpstr[256];
  const string dashedline2 = "---------------------";

  if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
  {
    cout << dashedline2 << endl;
    cout << "Rxn dQdx load:" << endl;
    cout << "  name = " << getName() << endl;
  }
#endif

  tmpBool = loadQMatrix (dQdxMat);
  bsuccess = bsuccess && tmpBool;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
    cout << dashedline2 << endl;
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadQMatrix
//
// Purpose       : Loads the Q-vector contributions for a single
//                 instance.
//
// Special Notes : See the special notes for loadDAEdQdxMatrix
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 10/19/08
//-----------------------------------------------------------------------------
bool Instance::loadQMatrix (N_LAS_Matrix & dQdxMat)
{
  bool bsuccess = true;

  // Finally, we need to add in the block of jacobian entries for the
  // reactions
  int numRegions = regVec.size();
  for (int ireg=0;ireg<numRegions;++ireg)
  {
    regVec[ireg]->loadDAEdQdx (dQdxMat);
  }
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 instance.
//
// Special Notes : This is an algebraic constaint.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 02/09/08
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  bool bsuccess = true;
  bool tmpBool = true;

  N_LAS_Matrix & dFdxMat = *(extData.dFdxMatrixPtr);

#ifdef Xyce_DEBUG_DEVICE
  char tmpstr[256];
  const string dashedline2 = "---------------------";

  if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
  {
    cout << dashedline2 << endl;
    cout << "Rxn dFdx load:" << endl;
    cout << "  name = " << getName() << endl;
  }
#endif

  tmpBool = loadFMatrix (dFdxMat);
  bsuccess = bsuccess && tmpBool;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
  {
    cout << dashedline2 << endl;
  }
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadFMatrix
//
// Purpose       : Loads the F-vector contributions for a single
//                 instance.
//
// Special Notes : See the special notes for loadDAEdFdxMatrix
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 10/19/08
//-----------------------------------------------------------------------------
bool Instance::loadFMatrix (N_LAS_Matrix & dFdxMat)
{
  bool bsuccess=true;

  // Finally, we need to add in the block of jacobian entries for the
  // reactions
  int numRegions = regVec.size();
  int cSize(0);
  int rSize(0);
  if (numRegions > 0)
  {
    cSize = (regVec[0])->getNumConstants();
    rSize = (regVec[0])->getNumSpecies();
    for (int ireg=0;ireg<numRegions;++ireg)
    {
      regVec[ireg]->loadDAEdFdx(dFdxMat);
    }
    if (cols.size() < rSize) cols.resize(rSize,0);
    if (vals.size() < cols.size()) vals.resize(cols.size(),0.0);
  }

  // This set of loads should be updated to use the bracket operators.
  if (transportFlag)
  {
    if (!getSolverState().dcopFlag)
    {
      int numSpecies = thVec.size();
      int isp=0;
      for (isp=0;isp<numSpecies;++isp)
      {
        if (!(thVec[isp].transportFlag)) continue;

        int i=0;
        int size = regVec.size();

        double DiffC = thVec[isp].D_specie;
        vector<int> & specie_id = thVec[isp].specie_id;

        for (i=0;i<size;++i)
        {
          int row = specie_id[i];
          if (i==0)
          {
            double aveDx = (dxVec[i]);
            double coef = DiffC/(aveDx*dxVec[i]);
            int count = 2;

            cols[0] = specie_id[i  ]; vals[0] = thVec[isp].bcScale1 * coef;
            cols[1] = specie_id[i+1]; vals[1] =-coef;

            bool bs1 = dFdxMat.sumIntoLocalRow (row, count, &vals[0], &cols[0]);
          }
          else if (i==(size-1))
          {
            double aveDx = (dxVec[i-1]);
            double coef = DiffC/(aveDx*dxVec[i-1]);
            int count = 2;
            cols[0] = specie_id[i-1]; vals[0] =-coef;
            cols[1] = specie_id[i  ]; vals[1] = thVec[isp].bcScale2 * coef;

            bool bs1 = dFdxMat.sumIntoLocalRow (row, count, &vals[0], &cols[0]);
          }
          else
          {
            double aveDx = 0.5*(dxVec[i-1] + dxVec[i]);
            double coef1 = DiffC/(aveDx*dxVec[i-1]);
            double coef2 = DiffC/(aveDx*dxVec[i  ]);
            double coefSum = coef1+coef2;
            int count = 3;

            cols[0] = specie_id[i-1]; vals[0] =-coef1;
            cols[1] = specie_id[i  ]; vals[1] = coefSum;
            cols[2] = specie_id[i+1]; vals[2] =-coef2;

            bool bs1 = dFdxMat.sumIntoLocalRow (row, count, &vals[0], &cols[0]);
          }
        } // for loop over regions
      } // for loop over species.
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 02/09/08
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState()
{
  bool bsuccess = true;

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline2 = "---------------------";

  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << dashedline2 << endl;
    cout << "  Start Instance::updatePrimaryState\n";
    cout << "  name = " << getName() << endl;
  }
#endif

  // Do the bulk of the work in updateIntermediateVars:
  updateIntermediateVars();

  // From this point onward, the function does nothing if running MPDE.
  double * staVector = extData.nextStaVectorRawPtr;
  double * currStaVector = extData.currStaVectorRawPtr;

  // We need time derivatives of concentrations from time integrator to
  // calculate RHS for old DAE, so save the concentrations in state.
  int numRegions = regVec.size();
  for (int ireg=0;ireg<numRegions;++ireg)
  {
    if (regVec[ireg]->haveAnyReactions())
    {
      int rSize = regVec[ireg]->getNumSpecies();
      for (int i=0;i<rSize;++i)
      {
        staVector[regVec[ireg]->getStateConcentrationLID(i)] = regVec[ireg]->getStateConcentration(i);
      }
    }
  }

  // if this is the first newton step of the first time step
  // of the transient simulation, we need to enforce that the
  // time derivatives w.r.t. charge are zero.  This is to maintain 3f5
  // compatibility.  ERK.

  if (!(getSolverState().dcopFlag) && (getSolverState().initTranFlag) && getSolverState().newtonIter==0)
  {
    // Also need to force these derivatives to be zero
    numRegions = regVec.size();
    for (int ireg=0;ireg<numRegions;++ireg)
    {
      if (regVec[ireg]->haveAnyReactions())
      {
        int rSize = regVec[ireg]->getNumSpecies();
        for (int i=0;i<rSize;++i)
        {
          currStaVector[regVec[ireg]->getStateConcentrationLID(i)] = regVec[ireg]->getStateConcentration(i);
        }
      }
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 02/09/08
//-----------------------------------------------------------------------------
bool Instance::updateSecondaryState()
{
  double * staDeriv = extData.nextStaDerivVectorRawPtr;

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline2 = "---------------------";

  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << dashedline2 << endl;
    cout << "  Start Instance::updateSecondaryState\n";
    cout << "  name = " << getName() << endl;
  }
#endif

  // We need time derivatives of concentrations from time integrator
  int numRegions = regVec.size();
  for (int ireg=0;ireg<numRegions;++ireg)
  {
    regVec[ireg]->updateSecondaryState (staDeriv);
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter
// Creation Date : 02/09/08
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars()
{
  double * solVec = extData.nextSolVectorRawPtr;
  double * oldSolVec = extData.currSolVectorRawPtr;

  double * staVec = extData.nextStaVectorRawPtr;
  double * currStaVec = extData.currStaVectorRawPtr;

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline2 = "---------------------";
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << dashedline2 << endl;
    cout << "  In Instance::updateIntermediateVars" << endl;
    cout << "  name = " << getName() << endl;
  }
#endif

  // Reset this to show we're starting from the originals
  origFlag = true;

  // update the "old" newton iteration number.
  if (getSolverState().newtonIter != 0 && getSolverState().newtonIter != newtonIterOld)
  {
    newtonIterOld = getSolverState().newtonIter;
  }

  double time=getSolverState().currTime;

  // For reaction kinetics, all we have to do here is copy the concentrations
  // out of the solution vector and into the temp storage
  int numRegions = regVec.size();
  for (int ireg=0;ireg<numRegions;++ireg)
  {
    regVec[ireg]->setMasterSourceValue(model_.masterSource);
    regVec[ireg]->updateIntermediateVars (solVec, oldSolVec, time);
  }


  int numSpecies = thVec.size();
  int isp=0;
  if (transportFlag)
  {
    for (isp=0;isp<numSpecies;++isp)
    {
      if (!(thVec[isp].transportFlag)) continue;

      int i=0;
      int size = thVec[isp].fluxVec.size(); // = numRegions-1
      double DiffC = thVec[isp].D_specie;
      for (i=0;i<size;++i)
      {
        double n2 = solVec[thVec[isp].specie_id[i  ]];
        double n1 = solVec[thVec[isp].specie_id[i+1]];

        thVec[isp].fluxVec[i] =  DiffC*(n2-n1)/dxVec[i];
      }

      // this imposes a direchlet BC (n=0) condition on all species at the ends.
      if (dirichletBCFlag)
      {
        double n1(0.0); // i+1
        double n2(0.0); // i

        n1 = solVec[thVec[isp].specie_id[0]];
        n2 = 0.0; // i=-1
        thVec[isp].flux_bc1 =  DiffC*(n2-n1)/dxVec[0];

        n1 = 0.0;  // i=max.  Remember "size" here is # fluxes, not # regions (=fluxes+1).
        n2 = solVec[thVec[isp].specie_id[size]];
        thVec[isp].flux_bc2 =  DiffC*(n2-n1)/dxVec[size];

        thVec[isp].bcScale1=2.0;
        thVec[isp].bcScale2=2.0;
      }
      else
      {
        thVec[isp].flux_bc1 = 0.0;
        thVec[isp].flux_bc2 = 0.0;
      }
    }
  }

  return true;
}


//-----------------------------------------------------------------------------
// Function      : Instance::getInstanceBreakPoints
// Purpose       : This function adds break points to a vector of breakpoints.
//
//                 It does not bother to check them in any way, or put them
//                 in order.  It only adds them in.
//
// Special Notes : Breakpoints other than those from the photocurrent
//                 in this device are all generated by the reaction network
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 02/09/08
//-----------------------------------------------------------------------------
bool Instance::getInstanceBreakPoints(
  vector<N_UTL_BreakPoint> &breakPointTimes)
{
  int numRegions = regVec.size();
  double junk;
  for (int ireg=0;ireg<numRegions;++ireg)
  {
    regVec[ireg]->setSimTime(getSolverState().currTime);
    junk=regVec[ireg]->getBreakTime();
    breakPointTimes.push_back(junk);
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Model::constructComposite
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 02/09/08
//-----------------------------------------------------------------------------
CompositeParam *Model::constructComposite
(string & cName, string & pName)
{
#ifdef Xyce_DEBUG_DEVICE
  //cout << "Model::constructComposite  name = " << name << "   cName = " << cName <<endl;
#endif

  if (cName == "DOPINGPROFILES" || cName == "REGION")
  {
    DopeInfo *doping = new DopeInfo();
    dopeInfoMap[pName] = doping;
    return (static_cast<CompositeParam *> (doping));
  }
  else if (cName == "SOURCELIST")
  {
    SpecieSource *sources = new SpecieSource();
    defectSourceMap[pName] = sources;
    return (static_cast<CompositeParam *> (sources));
  }
  else
  {
    string msg =
      "Model::constructComposite: unrecognized composite name: ";
    msg += cName;
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
  }

  return NULL;
}

//-----------------------------------------------------------------------------
// Function      : Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/09/08
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
// Creator       : Eric R. Keiter, SNL
// Creation Date : 02/09/08
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
// Purpose       : modelblock constructor
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date : 02/09/08
//-----------------------------------------------------------------------------
Model::Model(const ModelBlock & MB,
             SolverState & ss1,
             DeviceOptions & do1)
  : DevicePDEModel(MB,ss1,do1),
    TNOM(300.0),

    userNumRegions(0),
    rxnFileName("NOFILE"),
    xlo(1.0e-5),
    xhi(3.0e-4),

    xlo_source(1.0e-5),
    xhi_source(3.0e-4),
    xlo_sourceGiven(false),
    xhi_sourceGiven(false),

    masterSource(0.0)
{

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to .model line and constant defaults from metadata:
  setModParams (MB.params);

  // Set any non-constant parameter defaults:
  if (!given("TNOM"))
    TNOM = getDeviceOptions().tnom;

  // if the specification makes no sense, then turn it off:
  if (xlo_source >= xhi_source)
  {
    xlo_sourceGiven = false;
    xhi_sourceGiven = false;
    string msg = "XLO_SOURCE >= XHI_SOURCE.  Ignoring, and using a spatially uniform source";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_WARNING, msg);
  }

  // If the source region was not specified(or turned off), make sure it
  // extends through the entire integration volume:
  if (!xlo_sourceGiven)
  {
    xlo_source = xlo;
  }

  if (!xhi_sourceGiven)
  {
    xhi_source = xhi;
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
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/09/08
//-----------------------------------------------------------------------------
Model::~Model()
{
  vector<Instance*>::iterator iter;
  vector<Instance*>::iterator first = instanceContainer.begin();
  vector<Instance*>::iterator last  = instanceContainer.end();

  for (iter=first; iter!=last; ++iter)
  {
    delete (*iter);
  }

  if (!(regionDataVec.empty()))
  {
    int size = regionDataVec.size();
    int i=0;
    for (i=0;i<size;++i)
    {
      if (regionDataVec[i] != 0)
      {
        delete regionDataVec[i];
        regionDataVec[i] = 0;
      }
    }
  }

  // Loop over the dopeInfoMap (if it is not empty) and delete its contents.
  if (!(dopeInfoMap.empty()))
  {
    map <string,DopeInfo *>::iterator iter;
    map <string,DopeInfo *>::iterator begin = dopeInfoMap.begin();
    map <string,DopeInfo *>::iterator end   = dopeInfoMap.end  ();

    for(iter=begin;iter!=end;++iter)
    {
      if (iter->second != 0) delete iter->second;
    }
  }

  // Do the same for the defect source map.
  if (!(defectSourceMap.empty()))
  {
    map <string,CompositeParam *>::iterator iter;
    map <string,CompositeParam *>::iterator begin = defectSourceMap.begin();
    map <string,CompositeParam *>::iterator end   = defectSourceMap.end  ();

    for(iter=begin;iter!=end;++iter)
    {
      if (iter->second != 0)
      {
        delete iter->second;
        iter->second=0;
      }
    }
    defectSourceMap.clear();
  }

}

//-----------------------------------------------------------------------------
// Function      : Model::printOutInstances
// Purpose       : debugging tool.
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date : 02/09/08
//-----------------------------------------------------------------------------
std::ostream &Model::printOutInstances(std::ostream &os) const
{
  vector<Instance*>::const_iterator iter;
  vector<Instance*>::const_iterator first = instanceContainer.begin();
  vector<Instance*>::const_iterator last  = instanceContainer.end();

  int i;
  os << "     name     getModelName()  Parameters" << endl;
  for (i = 0, iter = first; iter != last; ++iter, ++i)
  {
    os << "  " << i << ": " << (*iter)->getName() << "        ";
    os << (*iter)->getModelName();

    os << endl;
    os << "  TEMP  = " << (*iter)->TEMP  << endl;

    os << endl;
  }

  os << endl;

  return os;
}

} // namespace RxnSet
} // namespace Device
} // namespace Xyce

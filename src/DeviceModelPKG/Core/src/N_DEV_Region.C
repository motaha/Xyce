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
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_Region.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL
//
// Creation Date  : 07/19/06
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.17 $
//
// Revision Date  : $Date: 2014/02/24 23:49:15 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

#include <N_UTL_Misc.h>
#include <N_UTL_Xyce.h>

#ifdef Xyce_DEBUG_DEVICE
#include <iostream>
#endif

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
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_RegionData.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_MultiVector.h>
#include <N_LAS_Vector.h>

#include <N_DEV_Region.h>
#include <N_DEV_Specie.h>
#include <N_DEV_ScalingVars.h>
#include <N_UTL_Expression.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : Region::Region
// Purpose       : constructor 1
// Special Notes : creates reaction region with default reaction network
// Scope         : public
// Creator       :
// Creation Date : 8/24/06
//-----------------------------------------------------------------------------
Region::Region( RegionData & rd,
                const DeviceOptions & devOp,
                const SolverState & solst,
                bool sourceOn)
  : name(rd.name),
    outputName(rd.outName),
    devOptions(devOp),
    solState(solst),
    regData(rd),
    useScaledVariablesFlag(false),
    variablesScaledFlag(false),
    rateConstantsScaledFlag(false),
    callsOTEC(0),
    baseReactionIndex(-1),
    x0(0.0),
    a0(0.0),
    C0(0.0),
    D0(0.0),
    u0(0.0),
    R0(0.0),
    rR0(0.0),
    t0(0.0),
    k0(0.0),
    rt0(0.0),
    rk0(0.0),
    outputBefore1(false),
    outputBefore2(false)
{
  theReactions.setApplySources(sourceOn);
  createDefaultReactionNetwork( rd.reactionFile );
}

//-----------------------------------------------------------------------------
// Function      : Region::Region
// Purpose       : constructor 2
// Special Notes : creates reaction region with copied reaction network
//                 for now does nothing but call initializers, but moved
//                 out of header file in case we want it to be more than a
//                 simple inlined deal.
// Scope         : public
// Creator       :
// Creation Date : 7/29/06
//-----------------------------------------------------------------------------
Region::Region ( RegionData & rd,
                 const DeviceOptions & devOp,
                 const SolverState & solst,
                 ReactionNetwork & reactionNet)
  
  : name(rd.name),
    outputName(rd.outName),
    devOptions(devOp),
    solState(solst),
    regData(rd),
    useScaledVariablesFlag(false),
    variablesScaledFlag(false),
    rateConstantsScaledFlag(false),
    callsOTEC(0),
    baseReactionIndex(-1),
    x0(0.0),
    a0(0.0),
    C0(0.0),
    D0(0.0),
    u0(0.0),
    R0(0.0),
    rR0(0.0),
    t0(0.0),
    k0(0.0),
    rt0(0.0),
    rk0(0.0),
    outputBefore1(false),
    outputBefore2(false),
    theReactions(reactionNet)
{
#ifdef Xyce_DEBUG_DEVICE
  if (devOptions.debugLevel > 0)
  {
    Xyce::dout() << theReactions << std::endl;
  }
#endif
}

//-----------------------------------------------------------------------------
// Function      : Region::~Region
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date : 7/29/06
//-----------------------------------------------------------------------------
Region::~Region ()
{
  theReactions.clear();
}

//-----------------------------------------------------------------------------
// Function      : Region::createDefaultRectionNetwork
// Purpose       : like it says
// Special Notes :
// Scope         : private
// Creator       :
// Creation Date : 7/29/06
//-----------------------------------------------------------------------------
void
Region::createDefaultReactionNetwork(const std::string &reactionSpecFile)
{
  theReactions.clear();

  theReactions.setReactionNetworkFromFile(reactionSpecFile);

#ifdef Xyce_DEBUG_DEVICE
  if (devOptions.debugLevel > 0)
  {
    Xyce::dout() << theReactions << std::endl;
  }
#endif

}

//-----------------------------------------------------------------------------
// Function      : Region::initializeReactionNetwork
// Purpose       : Initialize the 0-d reaction network equations
// Special Notes : The reaction network itself has already been created when
//                 the region was set up, here we just initialize some
//                 properties of it.
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 3/29/06
//-----------------------------------------------------------------------------
void Region::initializeReactionNetwork(ScalingVars & sv)
{

  int cSize=theReactions.getNumConstants();
  int sSize=theReactions.getNumSpecies();

  setupScalingVars(sv);

  theReactions.setScaleParams(C0,t0,x0);

  // We will need to get the initial conditions in some other way, too.
  // this presumes we have all these specific named species.  We might have
  // more, or fewer.
  if (cSize > 0)
  {
    theConstantConcentrations.resize(cSize,0.0);
  }

  // initial condition
  // These are just this way for testing purposes.
  if (sSize > 0)
  {
    initialConcentrations.resize(sSize,0.0);
  }

  // Estimate the electron and hole concentrations based on the
  // doping levels.  For now we are just considering the exact center
  // of the depletion region.

  if (sSize > 0 ) // don't initialize any concentrations at all if there
                  // are no variable species!
  {
    int icSize = theReactions.getNumInitialConditions();
    int i;
    for (i=0; i<icSize; ++i)
    {
      std::pair<std::string,double> theIC=theReactions.getInitialCondition(i);
      int specIndex=theReactions.getReactantNum(theIC.first);
      if (specIndex < 0)
        theConstantConcentrations[-(specIndex+1)] = theIC.second;
      else
        initialConcentrations[specIndex] = theIC.second;
    }
  }

}

//-----------------------------------------------------------------------------
// Function      : Region::setInitialCondition
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date : 3/29/06
//-----------------------------------------------------------------------------
void Region::setInitialCondition(const std::string & reactantname, const double val)
{
  int specIndex=theReactions.getReactantNum(reactantname);
  initialConcentrations[specIndex] = val;
}

//-----------------------------------------------------------------------------
// Function      : Region::addSource
// Purpose       : add a source term to the reaction network
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/4/06
//-----------------------------------------------------------------------------
void Region::addSource(std::string speciesName, Util::Expression *expr)
{
  theReactions.addSourceTerm(speciesName, expr);
}

//-----------------------------------------------------------------------------
// Function      : Region::addMasterSource
// Purpose       : add a master source term to the reaction network
// Special Notes : "master" source terms are those that depend on a single
//                  master source value that is set each time we evaluate
//                  the RHS.
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/4/06
//-----------------------------------------------------------------------------
void Region::addMasterSource(std::string speciesName)
{
  theReactions.addMasterSourceTerm(speciesName);
}

//-----------------------------------------------------------------------------
// Function      : Region::setRateConstants
// Purpose       : Set rate constants for reactions defined in
//                 initializeReactionNetwork.  This is meant to be called
//                 from updateTemperature.
// Special Notes : Cheesy hard-coded version for initial testing
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 3/29/06
//-----------------------------------------------------------------------------
void Region::setRateConstants(double T)
{
  theReactions.setRateConstantsFromCalc(T);

#ifdef Xyce_DEBUG_DEVICE

  if (!outputBefore1)
  {
    if (devOptions.debugLevel > 0)
      Xyce::dout() << theReactions << std::endl;
    outputBefore1 = true;
  }

#endif

  // if the variables have been scaled, then the rate constants should be
  // as well:
  if (variablesScaledFlag)
  {
    scaleRateConstants ();
  }

#ifdef Xyce_DEBUG_DEVICE
  if (variablesScaledFlag && !outputBefore2)
  {
    if (devOptions.debugLevel > 0)
    {
      Xyce::dout() << "Scaled Reactions: " << std::endl;
      Xyce::dout() << theReactions << std::endl;
    }
    outputBefore2=true;
  }
#endif
}

//----------------------------------------------------------------------------
// Function      : Region::setupScalingVars
// Purpose       :
//
// Special Notes : Most of these aren't used.  The important ones are
//                 T0 and C0.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 6/01/06
//----------------------------------------------------------------------------
void Region::setupScalingVars (ScalingVars & sv)
{
  t0 = sv.t0;
  C0 = sv.C0;
  x0 = sv.x0;
  double rx0 = 1.0/x0;

  a0 = sv.a0;

  D0  = sv.D0;
  R0  = sv.R0;

  rR0 = 1.0/R0;

  rk0 = C0*t0;
  rt0 = 1.0/t0;
  k0 = 1.0/rk0;

#ifdef Xyce_DEBUG_DEVICE
  if (devOptions.debugLevel > -2)
  {
    Xyce::dout() << "N_DEV_Region::setupScalingVars:" << std::endl;
    Xyce::dout().width(20);
    Xyce::dout().precision(12);
    Xyce::dout() << " " << std::endl;
    Xyce::dout() << "t0 = " << t0 << std::endl;
    Xyce::dout() << "x0 = " << x0 << std::endl;
    Xyce::dout() << "a0 = " << a0 << std::endl;
    Xyce::dout() << "C0 = " << C0 << std::endl;
    Xyce::dout() << "D0 = " << D0 << std::endl;
    Xyce::dout() << "R0 = " << R0 << std::endl;
    Xyce::dout() << "k0 = " << k0 << std::endl;
    Xyce::dout() << "rk0 = " << rk0 << std::endl;
  }
#endif

}

//----------------------------------------------------------------------------
// Function      : Region::scaleVariables
// Purpose       : Applies the scaling variables.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 6/01/06
//----------------------------------------------------------------------------
void Region::scaleVariables ()
{
  int i=0;
  int size=0;

  // vector of constant concentrations (species held fixed)
  size = theConstantConcentrations.size();
  for (i=0;i<size;++i)
  {
    theConstantConcentrations[i] /= C0;
  }

  // working storage for communicating between updateIntermediateVars
  // and updatePrimaryState
  size = tempConcentrations.size();
  for (i=0;i<size;++i)
  {
    tempConcentrations[i] /= C0;
  }

  // initial conditions
  size = initialConcentrations.size();
  for (i=0;i<size;++i)
  {
    initialConcentrations[i] /= C0;
  }

  variablesScaledFlag = true;
}

//----------------------------------------------------------------------------
// Function      : Region::unscaleVariables
// Purpose       : Restores things to normal units.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 6/01/06
//----------------------------------------------------------------------------
void Region::unscaleVariables ()
{
  int i=0;
  int size=0;

  // vector of constant concentrations (species held fixed)
  size = theConstantConcentrations.size();
  for (i=0;i<size;++i)
  {
    theConstantConcentrations[i] *= C0;
  }

  // working storage for communicating between updateIntermediateVars
  // and updatePrimaryState
  size = tempConcentrations.size();
  for (i=0;i<size;++i)
  {
    tempConcentrations[i] *= C0;
  }

  // initial conditions
  size = initialConcentrations.size();
  for (i=0;i<size;++i)
  {
    initialConcentrations[i] *= C0;
  }

  variablesScaledFlag = false;
}

//----------------------------------------------------------------------------
// Function      : Region::scaleRateConstants
// Purpose       :
// Special Notes : This function assumes that the rate constants are
//                 in units of cm^3/s, need to be multiplied by C0*t0.
//
//                 Some rate constants include e- or h+ densities in them.
//                 ie, are k*Ne or k*Nh, rather than just k.  For these
//                 rate constants, the scalar multiple should just be t0.
//
//                 To understand this consider the expression:
//
//                 k * Na * Nb, where k is a rate and Na and Nb are
//                 the concentrations of reactants.  If the variables have
//                 been scaled already, then Na and Nb are already scaled
//                 to dimensionless units.  (ie Na /= C0).
//
//                 k, being in cm^3/s, needs to be multiplied by t0/C0.  So,
//                 the total scaling multiplied by the expression is:
//
//                 (rate scaling) * (concentration scaling)^2
//                 (t0*C0) * 1/C0 * 1/C0 = t0/C0.
//
//                 This is (correctly) the same as the scaling applied
//                 to dNa/dt, which will also be multiplied by t0/C0.
//
//                 In the case of an implicit e- density, then k is actually
//                 k*Ne, in units of s^-1.  Then the expression is:
//
//                 kNe * Nb.  The scalar  will be:
//
//                 t0 * 1/C0 = t0/C0, which is correct.
//
//                 So the punchline is:
//
//                 (1) rates including Ne or Nh, multiply by t0.
//                 (2) rates not including Ne or Nh, multiply by rk0=C0*t0.
//
// For now, just like the other reaction-specific functions, this function
// is hardwired to a specific reaction network.  It will have to be
// generalized later.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 6/01/06
//----------------------------------------------------------------------------
void Region::scaleRateConstants ()
{
  theReactions.scaleRateConstantsFromCalc();

  rateConstantsScaledFlag = true;
}

//----------------------------------------------------------------------------
// Function      : Region::unscaleRateConstants
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date : 6/01/06
//----------------------------------------------------------------------------
void Region::unscaleRateConstants ()
{
  theReactions.unscaleRateConstantsFromCalc();

  rateConstantsScaledFlag = false;
}

//-----------------------------------------------------------------------------
// Function      : Region::outputTecplot
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 7/19/06
//-----------------------------------------------------------------------------
bool Region::outputTecplot ()
{
  bool bsuccess = true;

  int i;
  char filename[32];   for(i=0;i<32;++i) filename[i] = static_cast<char>(0);

  sprintf(filename,"%s.dat",outputName.c_str());

  double time = solState.currTime;
  double step = solState.currTimeStep;

#ifdef Xyce_DEBUG_DEVICE
  if (devOptions.debugLevel > 0 && solState.debugTimeFlag)
  {
    Xyce::dout() << std::endl;
    Xyce::dout() << section_divider << std::endl;
    Xyce::dout() << "In N_DEV_DiodePDEInstance::outputTecplot.  filename = ";
    Xyce::dout() << std::string(filename);
    Xyce::dout() << " time = " << time;
    Xyce::dout() << " step = " << step;
    Xyce::dout() << std::endl;
  }
#endif

  FILE *fp1;

  if (callsOTEC <= 0)
  {
    fp1 = fopen(filename,"w");
  }
  else
  {
    fp1 = fopen(filename,"a");
  }

  // If this is the first call, print the header.
  int rSize=theReactions.getNumSpecies();
  int cSize=theReactions.getNumConstants();
  if (callsOTEC <= 0)
  {
    fprintf(fp1,
   " TITLE = \"Data for reaction model: %s.\",\n", name.c_str());
    fprintf(fp1,"%s","\tVARIABLES = \"TIME \",\n");

    for (i=0;i<cSize;++i)
    {
      fprintf(fp1, "\t    \"%s \",\n" , (theReactions.getConstantsName(i)).c_str());
    }

    for (i=0;i<rSize;++i)
    {
      fprintf(fp1,"\t    \"%s \",\n", (theReactions.getSpeciesName(i)).c_str());
    }
    fprintf(fp1,"%s","\tZONE F=POINT\n");
  }

  // Now print the data:
  fprintf(fp1,"  %20.12e", solState.currTime);


  for (i=0;i<cSize;++i)
  {
    double printVal =
      theConstantConcentrations[i] *((variablesScaledFlag)?(C0):(1.0));
    fprintf(fp1,"  %20.12e", printVal);
  }

  for (i=0;i<rSize;++i)
  {
    double printVal =
      tempConcentrations[i]*((variablesScaledFlag)?(C0):(1.0));
    fprintf(fp1,"  %20.12e", printVal);
  }

  fprintf(fp1,"%s","\n");

  ++callsOTEC;
  fclose(fp1);

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Region::setupJacStamp
//
// Purpose       : adds reaction region's information to jac stamp.
//
// Special Notes : Modified to return the base index.
//
// Scope         : public
// Creator       : Eric R. Keiter
// Creation Date : 7/19/06
//-----------------------------------------------------------------------------
void Region::setupJacStamp ( std::vector< std::vector<int> > & jacStamp,
                                      std::vector<int> &colDep,
                                      int & firstReactant,
                                      int & lastIndex )
{
  int concentrationSize = theReactions.getNumSpecies();

  if (concentrationSize != 0 && !(regData.doNothing) )  // then we have jacStamp modifications to do
  {
    if (colDep.size() != 2)
    {
      Report::DevelFatal0() << "Region::setupJacStamp ERROR:  colDep != 2 ";
    }

    int stampSize=jacStamp.size();
    baseReactionIndex=stampSize;
    jacStamp.resize(stampSize+concentrationSize);

    // This simply sets a full square block, irresepective of how the
    // network actually couples variables.  First cut, possibly inefficient,
    // probably not worth doing more specifically unless there's a huge
    // reaction network with many species and very sparse matrix.
    for (int i=0; i<concentrationSize; ++i)
    {
#if 0
      jacStamp[baseReactionIndex+i].resize(concentrationSize+2);
#else
      jacStamp[baseReactionIndex+i].resize(concentrationSize);
#endif
      for (int j=0;j<concentrationSize; ++j)
      {
        jacStamp[baseReactionIndex+i][j]=baseReactionIndex+j;
      }

#if 0
      // This is why colDep must be of size 2. (one for e-, one for h+)
      // add the colDep colummns:
      jacStamp[baseReactionIndex+i][concentrationSize  ]=colDep[0];
      jacStamp[baseReactionIndex+i][concentrationSize+1]=colDep[1];
#endif
    }
  }
  firstReactant = baseReactionIndex;
  lastIndex = jacStamp.size()-1;
}

//-----------------------------------------------------------------------------
// Function      : Region::registerLIDs
// Purpose       :
// Special Notes : updates the intIndex
// Scope         : public
// Creator       : Eric R. Keiter
// Creation Date : 7/19/06
//-----------------------------------------------------------------------------
void Region::registerLIDs( const std::vector<int> & intLIDVec,
                                    const std::vector<int> & extLIDVec,
                                    int & intIndex)
{

  if ( !(regData.doNothing) )
  {
    if (baseReactionIndex != -1)
    {
      int rSize = theReactions.getNumSpecies();
      li_Concentrations.clear();
      li_Concentrations.resize(rSize);
#ifdef Xyce_DEBUG_DEVICE
      if (devOptions.debugLevel > 0)
      {
        Xyce::dout() << "  li_Concentrations: " <<  std::endl;
      }
#endif
      for (int i=0;i<rSize;++i)
      {
        li_Concentrations[i] = intLIDVec[intIndex++];

#ifdef Xyce_DEBUG_DEVICE
        if (devOptions.debugLevel > 0)
        {
          Xyce::dout() << "  li_Concentrations["<< i << "] = " << li_Concentrations[i] << std::endl;
        }
#endif
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Region::augmentNameMap
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter
// Creation Date : 7/20/06
//-----------------------------------------------------------------------------
void Region::augmentNameMap (
    std::map<int,std::string> & intNameMap,
    DeviceInstance & di)
{

  std::string tmpstr;

  if ( !(regData.doNothing) )
  {
    if (baseReactionIndex != -1)
    {
      int rSize=theReactions.getNumSpecies();
      for (int i=0;i<rSize;++i)
      {
        std::string speciesName=theReactions.getSpeciesName(i);
        tmpstr = name+"_Conc_"+speciesName;
        di.spiceInternalName(tmpstr);
        intNameMap[li_Concentrations[i]] = tmpstr;
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Region::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter
// Creation Date : 7/20/06
//-----------------------------------------------------------------------------
void Region::registerStateLIDs
   (const std::vector<int> & staLIDVec, int & i)
{
  if (baseReactionIndex != -1)
  {
    int rSize = theReactions.getNumSpecies();
    li_state_Concentrations.resize(rSize);
    for (int j=0;j<rSize;++j)
    {
      li_state_Concentrations[j]=staLIDVec[i++];
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Region::registerJacLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter
// Creation Date : 7/20/06
//-----------------------------------------------------------------------------
void Region::registerJacLIDs
  ( const std::vector< std::vector<int> > & jacLIDVec,
    const std::vector<int> &map,
    const std::vector< std::vector<int> > &map2
  )
{
  if ( !(regData.doNothing) )
  {
    if (baseReactionIndex != -1)
    {
      int rSize = theReactions.getNumSpecies();
      AConcentrationEquConcentrationNodeOffsets.clear();
      AConcentrationEquAuxNodeOffsets.clear();

      AConcentrationEquConcentrationNodeOffsets.resize(rSize);
      AConcentrationEquAuxNodeOffsets.resize(rSize);

      for (int i=0;i<rSize;++i)
      {
        AConcentrationEquConcentrationNodeOffsets[i].resize(rSize);
        for (int j=0;j<rSize;++j)
        {
          AConcentrationEquConcentrationNodeOffsets[i][j]=
            jacLIDVec[map[baseReactionIndex+i]][map2[baseReactionIndex+i][j]];
        }

#if 0
        int bRI = baseReactionIndex;
        AConcentrationEquAuxNodeOffsets[i].resize(2);
        for (int j=0; j<2; ++j)
        {
          AConcentrationEquAuxNodeOffsets[i][j]=
            jacLIDVec[map[bRI+i]][map2[bRI+i][rSize+j]];
        }
#endif
      }
    }
  }

}

//-----------------------------------------------------------------------------
// Function      : Region::setupPointers()
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 5/19/10
//-----------------------------------------------------------------------------
void Region::setupPointers (N_LAS_Matrix & dfdx, N_LAS_Matrix & dqdx)
{
  if ( !(regData.doNothing) )
  {
    int rSize = theReactions.getNumSpecies();

    dfdxConcEquConcVarPtrs.clear();
    dfdxConcEquConcVarPtrs.resize(rSize);

    dqdxConcEquConcVarPtrs.clear();
    dqdxConcEquConcVarPtrs.resize(rSize);

#if 0
    dfdxConcEquAuxVarPtrs.clear();
    dfdxConcEquAuxVarPtrs.resize(rSize);

    dqdxConcEquAuxVarPtrs.clear();
    dqdxConcEquAuxVarPtrs.resize(rSize);
#endif

    for (int i=0;i<rSize;++i)
    {
      dfdxConcEquConcVarPtrs[i].resize(rSize,0);
      dqdxConcEquConcVarPtrs[i].resize(rSize,0);

      for (int j=0;j<rSize; ++j)
      {
        int lidRow = li_Concentrations[i];
        int lidCol = li_Concentrations[j];

        double * dfdxPtr = dfdx.returnRawEntryPointer (lidRow, lidCol);
        dfdxConcEquConcVarPtrs[i][j] = dfdxPtr;

        double * dqdxPtr = dqdx.returnRawEntryPointer (lidRow, lidCol);
        dqdxConcEquConcVarPtrs[i][j] = dqdxPtr;
      }

#if 0
      dfdxConcEquAuxVarPtrs[i].resize(2,0);
      dqdxConcEquAuxVarPtrs[i].resize(2,0);
      for (int j=0;j<2; ++j)
      {
        int lidRow = li_Concentrations[i];
        int coloffset = AConcentrationEquAuxNodeOffsets[i][j];

        dfdxConcEquAuxVarPtrs[i][j] = &(dfdx[lidRow][coloffset]);
        dqdxConcEquAuxVarPtrs[i][j] = &(dqdx[lidRow][coloffset]);
      }
#endif
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Region::getDoNothingFlag
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 5/25/10
//-----------------------------------------------------------------------------
bool Region::getDoNothingFlag ()
{
  return regData.doNothing;
}

//-----------------------------------------------------------------------------
// Function      : Region::getBreakTime
// Purpose       : Return next breakpoint time from reaction network
// Special Notes :
//
//
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 8/4/06
//-----------------------------------------------------------------------------
double Region::getBreakTime()
{
  return (theReactions.getBreakpointTime());
}

//-----------------------------------------------------------------------------
// Function      : Region::updateIntermediateVars
// Purpose       :
// Special Notes :
//
//    For reaction kinetics, all we have to do here is copy the concentrations
//    out of the solution vector and into the temp storage
//
// Scope         : public
// Creator       : Eric R. Keiter
// Creation Date : 7/19/06
//-----------------------------------------------------------------------------
void Region::updateIntermediateVars
    ( double * solVector, double * oldSolVector, double time)
{
  if (baseReactionIndex != -1)
  {
    int rSize = theReactions.getNumSpecies();
    if (tempConcentrations.size() != rSize)
    {
      tempConcentrations.clear();
      tempConcentrations.resize(rSize,0.0);
    }

    for (int i=0;i<rSize;++i)
    {
      tempConcentrations[i] = solVector[li_Concentrations[i]];
    }

    // Now update the imposed time-dependent generation rates:
    double scalar= ((variablesScaledFlag)?(t0/C0):(1.0));
    theReactions.setSimTime(time);
    theReactions.setSourceScaleFac(scalar);

    // Now compute ddt and jac for later use in the various loads
    if (ddt.size() != rSize)
    {
      ddt.clear();
      ddt.resize(rSize,0.0);
    }
    else
    {
      for (int i=0;i<rSize;++i)
      {
        ddt[i]=0.0;
      }
    }
    theReactions.getDdt( tempConcentrations, theConstantConcentrations, ddt);

    // Jacobian
    if (tempJac.size() != rSize)
    {
      tempJac.clear();
      tempJac.resize(rSize);
      for (int i=0; i<rSize; ++i)
      {
        tempJac[i].resize(rSize,0.0);
      }
    }
    else
    {
      for (int i=0;i<rSize;++i)
      {
        for (int j=0;j<rSize;++j)
        {
          tempJac[i][j]=0.0;
        }
      }
    }

    theReactions.getJac( tempConcentrations, theConstantConcentrations, tempJac);
  }
}

//-----------------------------------------------------------------------------
// Function      : Region::loadDAEQVector
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter
// Creation Date : 7/21/06
//-----------------------------------------------------------------------------
bool Region::loadDAEQVector (double * qVec)
{

  // ERK:  The Q-vector must be loaded! Not loading during DCOP is a mistake!
  //
  // The Q-vector does not get summed into the residual during the DCOP.  However,
  // it has to be computed during the DCOP, or else dQ/dt will be enormous in the
  // first time step out of the DCOP.  (it will be (Q_1 - 0.0)/dt rather than (Q_1-Q_0)/dt.
  //
  // dQdt matrix does get summed into the jacobian, however, and this is because the
  // for some devices (particularly the stand-alone capacitor), excluding dQdx will result
  // in a singular jacobian.
  //
  // Note that this term will appear to cause the numerical jacobian test to fail for dQdx
  // during the DCOP.  This is not a meaningful failure.  The dQdx jacobian test should be turned
  // off for the DCOP phase.

  //if (!solState.dcopFlag)
  {
    if (baseReactionIndex != -1)
    {
      int rSize=theReactions.getNumSpecies();

      for (int i=0;i<rSize;++i)
      {
        double val = tempConcentrations[i];
        if (variablesScaledFlag)
        {
          val *= t0; // need to double-check this.
        }
        qVec[li_Concentrations[i]] += val;
      }
    }
  }

  return true;
}


//-----------------------------------------------------------------------------
// Function      : Region::loadDAEFVector
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter
// Creation Date : 7/21/06
//-----------------------------------------------------------------------------
bool Region::loadDAEFVector (double * fVec)
{
  if (baseReactionIndex != -1)
  {
    int rSize=theReactions.getNumSpecies();

    if (solState.dcopFlag || regData.doNothing)
    {
      // satisfy initial condition at dcop
      for (int i=0;i<rSize;++i)
      {
        fVec[li_Concentrations[i]] += (tempConcentrations[i] - initialConcentrations[i]);
      }
    }
    else
    {
      // ddt from reaction rate should be equal to ddt from time integrator
      for (int i=0;i<rSize;++i)
      {
        fVec[li_Concentrations[i]] += (-ddt[i]);
      }
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Region::loadDAEdFdxdV
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter
// Creation Date : 10/26/06
//-----------------------------------------------------------------------------
bool Region::loadDAEdFdxdV (double * dfdxdv,double vdiff)
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Region::loadDAEdQdx
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter
// Creation Date : 7/21/06
//-----------------------------------------------------------------------------
bool Region::loadDAEdQdx (N_LAS_Matrix & dqdx)
{

  if (baseReactionIndex != -1)
  {
    if ( !(solState.dcopFlag || regData.doNothing) )
    {
      double val = ((variablesScaledFlag)?(t0):(1.0));
      int rSize = theReactions.getNumSpecies();
      for (int i=0;i<rSize; ++i)
      {
        *(dqdxConcEquConcVarPtrs[i][i]) += val;
      }
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Region::loadDAEdFdx
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter
// Creation Date : 7/21/06
//-----------------------------------------------------------------------------
bool Region::loadDAEdFdx (N_LAS_Matrix & dfdx)
{
  if (baseReactionIndex != -1)
  {
    int rSize = theReactions.getNumSpecies();
    if ( solState.dcopFlag || regData.doNothing )
    {
      // jacobian for dcop is just ones on the diagonals, nothing else
      for (int i=0; i<rSize; ++i)
      {
        *(dfdxConcEquConcVarPtrs[i][i]) += 1.0;
      }
    }
    else
    {
      for (int i=0;i<rSize; ++i)
      {
        for (int j=0;j<rSize; ++j)
        {
          double val = (-tempJac[i][j]);
          *(dfdxConcEquConcVarPtrs[i][j]) += val;
        }
      }
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Region::updateSecondaryState
// Purpose       : Obtains time derivatives of concentrations from
//                 time integrator
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter
// Creation Date : 7/21/06
//-----------------------------------------------------------------------------
bool Region::updateSecondaryState (double * staDeriv)
{
  if (baseReactionIndex != -1)
  {
   int rSize = theReactions.getNumSpecies();
    if (tempConcentrationDerivs.size() != rSize)
    {
      tempConcentrationDerivs.clear();
      tempConcentrationDerivs.resize(rSize,0.0);
    }

    for (int i=0;i<rSize;++i)
    {
      tempConcentrationDerivs[i] = staDeriv[li_state_Concentrations[i]];

      // If the variablesScaledFlag is true, then the concentrations
      // are already scaled, but the time is not.
      if (variablesScaledFlag)
      {
        tempConcentrationDerivs[i] *= t0;
      }
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Region::haveAnyReactions()
// Purpose       :
// Special Notes :
//
// TVR: This one is a legacy of how I was keeping track of when
// the network was initialized and jacobian set up.  All the places
// where this baseReactionIndex != -1 occur were really only checking
// that the thing had been initialized.
//
// ERK:  updated to include a test for "do nothing" flag.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 05/22/10
//-----------------------------------------------------------------------------
bool Region::haveAnyReactions()
{
  bool retVal = (baseReactionIndex != -1);

  if (regData.doNothing)
  {
    retVal = false;
  }

  return retVal;
}

//-----------------------------------------------------------------------------
// Function      : Region::getNumIntVars
//
// Purpose       : Return the number of internal variables (solution vars)
//                 this region is adding to the system.
//
// Special Notes : The number of internal variables is comprised of the
//                 number of species in the reaction network.
//
// Scope         : public
//
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 10/24/06
//-----------------------------------------------------------------------------
int Region::getNumIntVars()
{
  int hInt=0;
  int numSpeciesVars=0;

  if (!(regData.doNothing))
  {
    numSpeciesVars=getNumSpecies();
  }

  return (numSpeciesVars+hInt);
}

//-----------------------------------------------------------------------------
// Function      : Region::loadDeviceMask
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Electrical and Microsystems Modeling
// Creation Date : 5/15/07
//-----------------------------------------------------------------------------
bool Region::loadDeviceMask (N_LAS_Vector & mask)
{
  return true;
}

} // namespace Device
} // namespace Xyce

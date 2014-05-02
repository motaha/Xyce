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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_Neuron3.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Richard Schiek, Electrical and Microsytem Modeling
//
// Creation Date  : 06/10/09
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.33.2.2 $
//
// Revision Date  : $Date: 2014/03/06 23:33:43 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------
#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif

#include <N_UTL_Misc.h>

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_Neuron3.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_Message.h>
#include <N_ERH_ErrorMgr.h>

#include <N_DEV_Neuron.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>

namespace Xyce {
namespace Device {


namespace Neuron3 {


void Traits::loadInstanceParameters(ParametricData<Neuron3::Instance> &p)
{
// Set up map for double precision variables:
  p.addPar ("R", 0.0, false, ParameterType::NO_DEP,
          &Neuron3::Instance::rInt,
          &Neuron3::Instance::rIntGiven, U_OHMMM1, CAT_NONE, "Intracellular resistivity");

  p.addPar ("A", 0.0, false, ParameterType::NO_DEP,
          &Neuron3::Instance::radius,
          &Neuron3::Instance::radiusGiven, U_METER, CAT_NONE, "Segment radius");

  p.addPar ("L", 0.0, false, ParameterType::NO_DEP,
          &Neuron3::Instance::length,
          &Neuron3::Instance::lengthGiven, U_METER, CAT_NONE, "Cable length");

  p.addPar ("RPS", 1.0e-6, false, ParameterType::NO_DEP,
          &Neuron3::Instance::rIntPrevious,
          &Neuron3::Instance::rIntPreviousGiven, U_OHMMM1, CAT_NONE, "Previous segment, intracellular resistivity");

  p.addPar ("APS", 0.0, false, ParameterType::NO_DEP,
          &Neuron3::Instance::radiusPrevious,
          &Neuron3::Instance::radiusPreviousGiven, U_METER, CAT_NONE, "Previous segment, segment radius");

  p.addPar ("LPS", 0.0, false, ParameterType::NO_DEP,
          &Neuron3::Instance::lengthPrevious,
          &Neuron3::Instance::lengthPreviousGiven, U_METER, CAT_NONE, "Previous segment length");

  p.addPar ("RNS", 1.0e-6, false, ParameterType::NO_DEP,
          &Neuron3::Instance::rIntNext,
          &Neuron3::Instance::rIntNextGiven, U_OHMMM1, CAT_NONE, "Next segment, intracellular resistivity");

  p.addPar ("ANS", 0.0, false, ParameterType::NO_DEP,
          &Neuron3::Instance::radiusNext,
          &Neuron3::Instance::radiusNextGiven, U_METER, CAT_NONE, "Next segment, segment radius");

  p.addPar ("LNS", 0.0, false, ParameterType::NO_DEP,
          &Neuron3::Instance::lengthNext,
          &Neuron3::Instance::lengthNextGiven, U_METER, CAT_NONE, "Next segment length");

  // add non-double types
  p.addPar("N", 0, false, ParameterType::NO_DEP,  &Neuron3::Instance::nSeg ,
          &Neuron3::Instance::nSegGiven , U_NONE, CAT_NONE, "Number of segments" );
}

void Traits::loadModelParameters(ParametricData<Neuron3::Model> &p)
{
  // Set up map for normal (double) param variables:
  p.addPar ("CMEM", 0.0, false, ParameterType::NO_DEP,
          &Neuron3::Model::cMem,
          &Neuron3::Model::cMemGiven,
          U_FARAD, CAT_NONE, "Membrane capacitance");

  p.addPar ("GMEM", 0.0, false, ParameterType::NO_DEP,
          &Neuron3::Model::gMem,
          &Neuron3::Model::gMemGiven,
          U_OHMM1, CAT_NONE, "Membrane conductance");

  p.addPar ("VREST", 0.0, false, ParameterType::NO_DEP,
          &Neuron3::Model::vRest,
          &Neuron3::Model::vRestGiven,
          U_VOLT, CAT_NONE, "Resting potential");

  p.addPar ("EK", 0.0, false, ParameterType::NO_DEP,
          &Neuron3::Model::eK,
          &Neuron3::Model::eKGiven,
          U_VOLT, CAT_NONE, "Potassium resting potential");

  p.addPar ("GK", 0.0, false, ParameterType::NO_DEP,
          &Neuron3::Model::gK,
          &Neuron3::Model::gKGiven,
          U_OHMM1, CAT_NONE, "Potassium base conductance");

  p.addPar ("ENA", 0.0, false, ParameterType::NO_DEP,
          &Neuron3::Model::eNa,
          &Neuron3::Model::eNaGiven,
          U_VOLT, CAT_NONE, "Sodium resting potential");

  p.addPar ("GNA", 0.0, false, ParameterType::NO_DEP,
          &Neuron3::Model::gNa,
          &Neuron3::Model::gNaGiven,
          U_OHMM1, CAT_NONE, "Sodium base conductance");

  p.addPar ("R", 0.0, false, ParameterType::NO_DEP,
          &Neuron3::Model::rInt,
          &Neuron3::Model::rIntGiven,
          U_OHMMM1, CAT_NONE, "Intracellular resistivity");

  p.addPar ("A", 0.0, false, ParameterType::NO_DEP,
          &Neuron3::Model::radius,
          &Neuron3::Model::radiusGiven,
          U_METER, CAT_NONE, "Segment radius");

  p.addPar ("L", 0.0, false, ParameterType::NO_DEP,
          &Neuron3::Model::length,
          &Neuron3::Model::lengthGiven,
          U_METER, CAT_NONE, "Cable length");

  p.addPar ("RPS", 1.0e-6, false, ParameterType::NO_DEP,
          &Neuron3::Model::rIntPrevious,
          &Neuron3::Model::rIntPreviousGiven,
          U_OHMMM1, CAT_NONE, "Previous segment, intracellular resistivity");

  p.addPar ("APS", 0.0, false, ParameterType::NO_DEP,
          &Neuron3::Model::radiusPrevious,
          &Neuron3::Model::radiusPreviousGiven,
          U_METER, CAT_NONE, "Previous segment, segment radius");

  p.addPar ("LPS", 0.0, false, ParameterType::NO_DEP,
          &Neuron3::Model::lengthPrevious,
          &Neuron3::Model::lengthPreviousGiven,
          U_METER, CAT_NONE, "Previous segment length");

  p.addPar ("RNS", 1.0e-6, false, ParameterType::NO_DEP,
          &Neuron3::Model::rIntNext,
          &Neuron3::Model::rIntNextGiven,
          U_OHMMM1, CAT_NONE, "Next segment, intracellular resistivity");

  p.addPar ("ANS", 0.0, false, ParameterType::NO_DEP,
          &Neuron3::Model::radiusNext,
          &Neuron3::Model::radiusNextGiven,
          U_METER, CAT_NONE, "Next segment, segment radius");

  p.addPar ("LNS", 0.0, false, ParameterType::NO_DEP,
          &Neuron3::Model::lengthNext,
          &Neuron3::Model::lengthNextGiven,
          U_METER, CAT_NONE, "Next segment length");

  // add non-double types
  p.addPar("N", 0, false, ParameterType::NO_DEP,  &Neuron3::Model::nSeg ,
          &Neuron3::Model::nSegGiven , U_NONE, CAT_NONE, "Number of segments" );
}



//
// static class member inits
//


// Class Instance
//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Instance::processParams()
{
  // If there are any time dependent parameters, set their values at for
  // the current time.

  // now set the temperature related stuff.
  //updateTemperature(temp);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateTemperature
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Instance::updateTemperature ( const double & temp)
{
  bool bsuccess = true;
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock &         IB,
  Model &                       Miter,
  const FactoryBlock &          factory_block)
  : DeviceInstance(IB, configuration.getInstanceParameters(), factory_block),
    model_(Miter),
    rInt(0.0),
    radius(0.0),
    length(0.0),
    segArea(0.0),
    nSeg(0),
    rIntGiven(false),
    radiusGiven(false),
    lengthGiven(false),
    nSegGiven(false),
    rIntPrevious(0.0),
    radiusPrevious(0.0),
    lengthPrevious(0.0),
    rIntNext(0.0),
    radiusNext(0.0),
    lengthNext(0.0),
    rIntPreviousGiven(false),
    radiusPreviousGiven(false),
    lengthPreviousGiven(false),
    rIntNextGiven(false),
    radiusNextGiven(false),
    lengthNextGiven(false),
    kcl1Fvalue(0.0),
    kcl2Fvalue(0.0),
    dkcl1F_dVin(0.0),
    dkcl1F_dVs0(0.0),
    dkcl2F_dVout(0.0),
    dkcl2F_dVsn(0.0),
    li_Pos(0),
    li_Neg(0),
    APosEquPosNodeOffset(0),
    APosEquNextNodeOffset(0),
    ANegEquNegNodeOffset(0),
    ANegEquLastNodeOffset(0)
{
  // we have pased the model and instance parameters so now we can calculate
  // the number of internal vars
  numExtVars   = 2;  // input and output voltage

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line
  setParams (IB.params);

  // Set any non-constant parameter defaults:

  //if (!given("TEMP"))
  //  temp = getDeviceOptions().temp.dVal();

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:
  processParams ();

  // pull unspecified params out of the model if they weren't specified here
  if( !nSegGiven && model_.nSegGiven )
  {
    nSeg = model_.nSeg;
    nSegGiven = true;
  }
  if( !rIntGiven && model_.rIntGiven )
  {
    rInt = model_.rInt;
    rIntGiven = true;
  }
  if( !radiusGiven && model_.radiusGiven )
  {
    radius = model_.radius;
    radiusGiven = true;
  }
  if( !lengthGiven && model_.lengthGiven )
  {
    length = model_.length;
    lengthGiven = true;
  }

  // if nSeg is still unknown then estimate it via lambda-d rule (TO DO)
  if( !nSegGiven )
  {
    nSeg = 10;
  }

  // now that nSeg, length and radius are set calculate segment area
  segArea = 2.0 * M_PI * radius * length / nSeg;

  // now we can calculate the number of internal vars
  numIntVars   = nSeg * 4;   // ion channel vars + one internal voltage node var for each segment
  numStateVars = nSeg * 2;
  int numVars = numExtVars + numIntVars;


  //
  // i(n) - I(n)/A - g(n,n+1) * (V(n+1) - V(n)) - g(n,n-1) * (V(n-1) - V(n)) + Cm dV(n)/d(t) = 0
  //
  // Vin   : g(0,1) * (V(1) - Vin ) = 0
  // Vout  : g(n,n-1) * (V(n-1) - Vout) = 0
  // Vnode : i(n) - I(n)/A - g(n,n+1) * (V(n+1) - V(n)) - g(n,n-1) * (V(n-1) - V(n)) + Cm dV(n)/d(t) = 0
  //        plus node supporting equations (a, b, m)
  //
  // jacobian format
  //             Vin      Vout      V1   n   m   h   V2   n   m   h    V(nSeg)   n   m   h
  // kcl Vin    -g(0,1)            g(0,1)
  // kcl Vout           -g(n,n-1)                                      g(n,n-1)
  // kcl V1     yes                yes yes yes yes  yes
  // a                             yes yes
  // b                             yes     yes
  // m                             yes         yes
  // kcl V2                        yes               yes yes yes yes    yes
  // a 2                                             yes yes
  // b 2                                             yes     yes
  // m 2                                             yes         yes
  // kcl VnSeg            yes                        yes                 yes     yes yes yes
  // a nSeg                                                              yes     yes
  // b nSeg                                                              yes         yes
  // m nSeg                                                              yes             yes
  //
  // jacobian element count by row:
  // vin:  2
  // vout: 2
  // v1:   6
  // a1:   2
  // b1:   2
  // m1:   2
  // v2:   6
  // a2:   2
  // b2:   2
  // m2:   2
  // vnSeg:6
  // a:    2
  // b:    2
  // m:    2
  //
  // set up jacStamp
  if( jacStamp.empty() )       // redundant as jacStamp is not static for this device
  {                            // it can't be as each cable may have a different number of nodes
    jacStamp.resize(numVars);
    jacStamp[0].resize(2);
    jacStamp[0][0] = 0;
    jacStamp[0][1] = 1;
    jacStamp[1].resize(2);
    jacStamp[1][0] = 1;
    jacStamp[1][1] = numVars - 4;
    for( int i=2; i<numVars; i+=4)
    {
      jacStamp[i].resize(6);
      if( i == 2 )
      {
        jacStamp[i][0] = 0;
      }
      else
      {
        jacStamp[i][0] = i-4;
      }
      jacStamp[i][1] = i;
      jacStamp[i][2] = i+1;
      jacStamp[i][3] = i+2;
      jacStamp[i][4] = i+3;
      if( i==(numVars-4) )
      {
        jacStamp[i][5] = 1;
      }
      else
      {
        jacStamp[i][5] = i+4;
      }
      jacStamp[i+1].resize(2);
      jacStamp[i+1][0] = i;
      jacStamp[i+1][1] = i+1;
      jacStamp[i+2].resize(2);
      jacStamp[i+2][0] = i;
      jacStamp[i+2][1] = i+2;
      jacStamp[i+3].resize(2);
      jacStamp[i+3][0] = i;
      jacStamp[i+3][1] = i+3;
    }

  }

  /*
  // print out jacStamp
  int numRows = jacStamp.size();
  for( int i=0; i< numRows; i++ )
  {
  int numCol = jacStamp[i].size();
  Xyce::dout() << "jacStamp[ " << i << " ] = { ";
  for(int j=0; j<numCol; j++)
  {
  Xyce::dout() << jacStamp[i][j] << "  ";
  }
  Xyce::dout() << " } " <<  std::endl;
  }
  Xyce::dout() << std::endl;
  */

  // calculate segment conductivities used in load calls:
  // Note: conductivity to the previous and next segment is not symmetric if the segment length and/or radius is not are not equal
  // the full formula is:
  //   g(n,n') = radius * (radius')^2 / ( rInt segLength * ( segLength * (radius')^2 + segLength' * (radius)^2 ) )
  //
  // equation 6.30 Theoretical neuroscience: computational and mathematical modeling of neural systems, Peter Dayan and L. F. Abbot 2001
  // If we allow variable segment lengths and radii then we'll need to update this calculation
  //
  // Note: the above equation has the wrong units unless rInt (which is Ohm Length) is rLong (which is Ohm/Length).
  //
  gForward.resize(nSeg);
  gBackward.resize(nSeg);
  double segLength = length / nSeg;
  double rLong = rInt / (M_PI * radius * radius);  // longitudinal resistivity (ohm/length)
  gBackward[0] = radius * (radiusPrevious * radiusPrevious) / (rLong * segLength * ( segLength * radiusPrevious * radiusPrevious + lengthPrevious * radius * radius ));
  gForward[0] = radius * (radius * radius) / (rLong * segLength * ( segLength * radius * radius + segLength * radius * radius ));
  gBackward[nSeg-1] = radius * (radius * radius) / (rLong * segLength * ( segLength * radius * radius + segLength * radius * radius ));
  gForward[nSeg-1] = radius * (radiusNext * radiusNext) / (rLong * segLength * ( segLength * radiusNext * radiusNext + lengthNext * radius * radius ));
  for(int i=1; i<(nSeg-1); i++)
  {
    gBackward[i] = radius * (radius * radius) / (rLong * segLength * ( segLength * radius * radius + segLength * radius * radius ));
    gForward[i] = gBackward[i];
  }
  /*
    gBackward[0] = radius * (radiusPrevious * radiusPrevious) / (rInt * segLength * ( segLength * radiusPrevious * radiusPrevious + lengthPrevious * radius * radius ));
    gForward[0] = radius * (radius * radius) / (rInt * segLength * ( segLength * radius * radius + segLength * radius * radius ));
    gBackward[nSeg-1] = radius * (radius * radius) / (rInt * segLength * ( segLength * radius * radius + segLength * radius * radius ));
    gForward[nSeg-1] = radius * (radiusNext * radiusNext) / (rInt * segLength * ( segLength * radiusNext * radiusNext + lengthNext * radius * radius ));
    for(int i=1; i<(nSeg-1); i++)
    {
    gBackward[i] = radius * (radius * radius) / (rInt * segLength * ( segLength * radius * radius + segLength * radius * radius ));
    gForward[i] = gBackward[i];
    }
  */
  // allocate space for load and jacobian terms per segment
  // variable indecies loads
  li_Vol.resize(nSeg);
  li_nPro.resize(nSeg);
  li_mPro.resize(nSeg);
  li_hPro.resize(nSeg);
  li_KCurrentState.resize(nSeg);
  li_NaCurrentState.resize(nSeg);
  segFvalue.resize(nSeg);
  segQvalue.resize(nSeg);
  segNEquFvalue.resize(nSeg);
  segNEquQvalue.resize(nSeg);
  segMEquFvalue.resize(nSeg);
  segMEquQvalue.resize(nSeg);
  segHEquFvalue.resize(nSeg);
  segHEquQvalue.resize(nSeg);

  // jacobian elements
  segF_dVp.resize(nSeg);
  segF_dV.resize(nSeg);
  segF_dVn.resize(nSeg);
  segF_dn.resize(nSeg);
  segF_dm.resize(nSeg);
  segF_dh.resize(nSeg);
  segQ_dV.resize(nSeg);
  dnF_dV.resize(nSeg);
  dnF_dn.resize(nSeg);
  dnQ_dn.resize(nSeg);
  dmF_dV.resize(nSeg);
  dmF_dm.resize(nSeg);
  dmQ_dm.resize(nSeg);
  dhF_dV.resize(nSeg);
  dhF_dh.resize(nSeg);
  dhQ_dh.resize(nSeg);

  // state vars
  potassiumCurrent.resize(nSeg);
  sodiumCurrent.resize(nSeg);

}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
Instance::~Instance()
{
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
void Instance::registerLIDs(const std::vector<int> & intLIDVecRef,
                            const std::vector<int> & extLIDVecRef)
{
  AssertLIDs(intLIDVecRef.size() == numIntVars);
  AssertLIDs(extLIDVecRef.size() == numExtVars);

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    Xyce::dout() << std::endl << section_divider << std::endl;
    Xyce::dout() << "  Instance::registerLIDs" << std::endl;
    Xyce::dout() << "  name = " << getName() << std::endl;
  }
#endif

  // copy over the global ID lists.
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  li_Pos  = extLIDVec[0];
  li_Neg  = extLIDVec[1];
  for( int i=0, j=0; i<nSeg; i++, j+=4)
  {
    li_Vol[i]  = intLIDVec[j];
    li_nPro[i] = intLIDVec[j+1];
    li_mPro[i] = intLIDVec[j+2];
    li_hPro[i] = intLIDVec[j+3];
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 )
  {
    Xyce::dout() << "  li_Pos = " << li_Pos << std::endl
         << "  li_Neg = " << li_Neg << std::endl;
    for( int i=0; i<nSeg; i++ )
    {
      Xyce::dout() << "  li_Vol[ " << i << " ] = " << li_Vol[i] << std::endl
           << "  li_nPro[ " << i << " ] = " << li_nPro[i] << std::endl
           << "  li_mPro[ " << i << " ] = " << li_mPro[i] << std::endl
           << "  li_hPro[ " << i << " ] = " << li_hPro[i] << std::endl;
    }
  }
#endif

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 )
  {
    Xyce::dout() << section_divider << std::endl;
  }
#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::getIntNameMap
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
std::map<int,std::string> & Instance::getIntNameMap ()
{
  // set up the internal name map, if it hasn't been already.
  if (intNameMap.empty ())
  {
    std::string tmpstr;
    for( int i=0; i<nSeg; i++)
    {
      std::ostringstream segNumber;
      segNumber << i;

      tmpstr = getName() + "_" + "V" + segNumber.str();
      spiceInternalName (tmpstr);
      intNameMap[ li_Vol[i] ] = tmpstr;

      tmpstr = getName() + "_" + "N" + segNumber.str();
      spiceInternalName (tmpstr);
      intNameMap[ li_nPro[i] ] = tmpstr;

      tmpstr = getName() + "_" + "M" + segNumber.str();
      spiceInternalName (tmpstr);
      intNameMap[ li_mPro[i] ] = tmpstr;

      tmpstr = getName() + "_" + "H" + segNumber.str();
      spiceInternalName (tmpstr);
      intNameMap[ li_hPro[i] ] = tmpstr;
    }
  }
  return intNameMap;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs( const std::vector<int> & staLIDVecRef )
{
  AssertLIDs(staLIDVecRef.size() == numStateVars);

  // copy over the global ID lists.
  staLIDVec = staLIDVecRef;

  for( int i=0, j=0; i<nSeg; i++, j+=2)
  {
    li_KCurrentState[i] = staLIDVec[j];
    li_NaCurrentState[i] = staLIDVec[j+1];
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
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Instance::loadDeviceMask ()
{
  bool returnVal=false;
  N_LAS_Vector * maskVectorPtr = extData.deviceMaskVectorPtr;

//   Xyce::dout() << "Masking n, m and h" << std::endl;
//   (*maskVectorPtr)[li_nPro] = 0.0;
//   (*maskVectorPtr)[li_mPro] = 0.0;
//   (*maskVectorPtr)[li_hPro] = 0.0;
//   returnVal = true;

  return (returnVal);
}

//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
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
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );

  // external terminals
  APosEquPosNodeOffset = jacLIDVec[0][0];
  APosEquNextNodeOffset = jacLIDVec[0][1];
  ANegEquNegNodeOffset = jacLIDVec[1][0];
  ANegEquLastNodeOffset = jacLIDVec[1][1];
  /*
    Xyce::dout() << "APosEquPosNodeOffset = " << APosEquPosNodeOffset << std::endl;
    Xyce::dout() << "APosEquNextNodeOffset = " << APosEquNextNodeOffset << std::endl;
    Xyce::dout() << "ANegEquNegNodeOffset = " << ANegEquNegNodeOffset << std::endl;
    Xyce::dout() << "ANegEquLastNodeOffset = " << ANegEquLastNodeOffset << std::endl;
  */

  // internal variables
  SegVEqnVpreOffset.resize(nSeg);
  SegVEqnVsegOffset.resize(nSeg);
  SegVEqnVnexOffset.resize(nSeg);
  SegVEqnNOffset.resize(nSeg);
  SegVEqnMOffset.resize(nSeg);
  SegVEqnHOffset.resize(nSeg);
  NEquVNodeOffset.resize(nSeg);
  NEquNNodeOffset.resize(nSeg);
  MEquVNodeOffset.resize(nSeg);
  MEquMNodeOffset.resize(nSeg);
  HEquVNodeOffset.resize(nSeg);
  HEquHNodeOffset.resize(nSeg);
  for(int i=0, j=2; i<nSeg; i++, j+=4 )
  {
    SegVEqnVpreOffset[i] = jacLIDVec[j][0];
    SegVEqnVsegOffset[i] = jacLIDVec[j][1];
    SegVEqnNOffset[i] = jacLIDVec[j][2];
    SegVEqnMOffset[i] = jacLIDVec[j][3];
    SegVEqnHOffset[i] = jacLIDVec[j][4];
    SegVEqnVnexOffset[i] = jacLIDVec[j][5];

    NEquVNodeOffset[i] = jacLIDVec[j+1][0];
    NEquNNodeOffset[i] = jacLIDVec[j+1][1];
    MEquVNodeOffset[i] = jacLIDVec[j+2][0];
    MEquMNodeOffset[i] = jacLIDVec[j+2][1];
    HEquVNodeOffset[i] = jacLIDVec[j+3][0];
    HEquHNodeOffset[i] = jacLIDVec[j+3][1];
    /*
      Xyce::dout() <<  SegVEqnVpreOffset[i] << ", "
      << SegVEqnVsegOffset[i] << ", "
      << SegVEqnNOffset[i] << ", "
      << SegVEqnMOffset[i] << ", "
      << SegVEqnHOffset[i] << ", "
      << SegVEqnVnexOffset[i] << ", "
      << NEquVNodeOffset[i] << ", "
      << NEquNNodeOffset[i] << ", "
      << MEquVNodeOffset[i] << ", "
      << MEquMNodeOffset[i] << ", "
      << HEquVNodeOffset[i] << ", "
      << HEquHNodeOffset[i] << std::endl;
    */
  }
}


//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars()
{
  //Xyce::dout() << "Instance::updateIntermediateVars()" << std::endl;

  bool bsuccess = true;

  // here we take the current solutions for V1, V2, n, m and h
  // and use those to calculate all the terms needed for the next
  // load cycle (F, Q, dFdX, dQdX)
  N_LAS_Vector * solVectorPtr = extData.nextSolVectorPtr;

  double vIn = (*solVectorPtr)[li_Pos];
  double vOut = (*solVectorPtr)[li_Neg];

  // take care of the input and output nodes as they are different
  kcl1Fvalue = gForward[0] * ((*solVectorPtr)[li_Vol[0]] - vIn );
  dkcl1F_dVin = -gForward[0];
  dkcl1F_dVs0 =  gForward[0];
  kcl2Fvalue = gBackward[nSeg-1] * ((*solVectorPtr)[li_Vol[nSeg-1]] - vOut );
  dkcl2F_dVout = -gBackward[nSeg-1];
  dkcl2F_dVsn = gBackward[nSeg-1];
  //Xyce::dout() << "end update: " << vIn << ", " << vOut << ", " << gForward[0] << ", " << gBackward[nSeg-1] << std::endl;

  // loop over segments getting all the load and jacobian terms for each segment
  for( int i=0; i<nSeg; i++ )
  {
    // for this segment get the values of the local vars
    double vSeg  = (*solVectorPtr)[li_Vol[i]];
    double vNext = 0.0;
    if (i == (nSeg - 1))
    {
      vNext = vOut;
    }
    else
    {
      vNext = (*solVectorPtr)[li_Vol[i+1]];
    }
    double vPrev = 0.0;
    if (i == 0 )
    {
      vPrev = vIn;
    }
    else
    {
      vPrev = (*solVectorPtr)[li_Vol[i-1]];
    }
    double nVarSeg = (*solVectorPtr)[li_nPro[i]];
    double mVarSeg = (*solVectorPtr)[li_mPro[i]];
    double hVarSeg = (*solVectorPtr)[li_hPro[i]];

    // do the voltage equation for this node
    // get function and derivative values as we go.
    {
      // F part
      // use scoping to avoid lots of similar variable names
      const int numDeriv = 6;
      Sacado::Fad::SFad<double,6> vVar( numDeriv, 0, vSeg );
      Sacado::Fad::SFad<double,6> vVpr( numDeriv, 1, vPrev );
      Sacado::Fad::SFad<double,6> vVne( numDeriv, 2, vNext );
      Sacado::Fad::SFad<double,6> nVar( numDeriv, 3, nVarSeg );
      Sacado::Fad::SFad<double,6> mVar( numDeriv, 4, mVarSeg );
      Sacado::Fad::SFad<double,6> hVar( numDeriv, 5, hVarSeg );
      // parameters
      Sacado::Fad::SFad<double,6> gPrev( gBackward[i] );
      Sacado::Fad::SFad<double,6> gNext( gForward[i] );
      Sacado::Fad::SFad<double,6> gMemVar( model_.gMem * segArea );
      Sacado::Fad::SFad<double,6> vRestVar( model_.vRest );
      Sacado::Fad::SFad<double,6> gKVar( model_.gK * segArea );
      Sacado::Fad::SFad<double,6> eKVar( model_.eK );
      Sacado::Fad::SFad<double,6> gNaVar( model_.gNa * segArea );
      Sacado::Fad::SFad<double,6> eNaVar( model_.eNa );
      //Xyce::dout() << "segment update: " << i << ", " << vSeg << ", " << vPrev << ", " << vNext << ", " << gBackward[i] << ", " << gForward[i] << std::endl;
      // compute the value and derivative terms for KCL 1 F
      Sacado::Fad::SFad<double,6> resultFad;
      resultFad = kcl1EquF( vVar, vVpr, vVne, nVar, mVar, hVar, gPrev, gNext, gMemVar, vRestVar, gKVar, eKVar, gNaVar, eNaVar );
      segFvalue[i] = resultFad.val();
      segF_dV[i]   = resultFad.dx(0);
      segF_dVp[i]  = resultFad.dx(1);
      segF_dVn[i]  = resultFad.dx(2);
      segF_dn[i]   = resultFad.dx(3);
      segF_dm[i]   = resultFad.dx(4);
      segF_dh[i]   = resultFad.dx(5);
    }
    {
      // Q part
      const int numDeriv = 2;
      Sacado::Fad::SFad<double,2> vVar( numDeriv, 0, vSeg );

      // parameters
      Sacado::Fad::SFad<double,2> cMemVar( model_.cMem * segArea );

      Sacado::Fad::SFad<double,2> resultFad;
      resultFad    = kcl1EquQ( vVar, cMemVar );
      segQvalue[i] = resultFad.val();
      segQ_dV[i]   = resultFad.dx(0);

    }

    // n - equation
    {
      const int numDeriv = 2;
      Sacado::Fad::SFad<double,2> vVar( numDeriv, 0, vSeg );
      Sacado::Fad::SFad<double,2> nVar( numDeriv, 1, nVarSeg );
      // parameter
      Sacado::Fad::SFad<double,2> vRestVar( model_.vRest );

      Sacado::Fad::SFad<double,2> resultFad = nEquF( vVar, nVar, vRestVar);
      segNEquFvalue[i] = resultFad.val();
      dnF_dV[i]        = resultFad.dx(0);
      dnF_dn[i]        = resultFad.dx(1);
    }
    {
      const int numDeriv = 1;
      Sacado::Fad::SFad<double,1> nVar( numDeriv, 0, nVarSeg );

      Sacado::Fad::SFad<double,1> resultFad = nEquQ( nVar );
      segNEquQvalue[i] = resultFad.val();
      dnQ_dn[i]        = resultFad.dx(0);
    }

    // m - equation
    {
      const int numDeriv = 2;
      Sacado::Fad::SFad<double,2> vVar( numDeriv, 0, vSeg );
      Sacado::Fad::SFad<double,2> mVar( numDeriv, 1, mVarSeg );
      // parameter
      Sacado::Fad::SFad<double,2> vRestVar( model_.vRest );

      Sacado::Fad::SFad<double,2> resultFad = mEquF( vVar, mVar, vRestVar );
      segMEquFvalue[i] = resultFad.val();
      dmF_dV[i]        = resultFad.dx(0);
      dmF_dm[i]        = resultFad.dx(1);
    }
    {
      const int numDeriv = 1;
      Sacado::Fad::SFad<double,1> mVar( numDeriv, 0, mVarSeg );

      Sacado::Fad::SFad<double,1> resultFad = mEquQ( mVar );
      segMEquQvalue[i] = resultFad.val();
      dmQ_dm[i]        = resultFad.dx(0);
    }

    // h - equation
    {
      const int numDeriv = 2;
      Sacado::Fad::SFad<double,2> vVar( numDeriv, 0, vSeg );
      Sacado::Fad::SFad<double,2> hVar( numDeriv, 1, hVarSeg );
      // parameter
      Sacado::Fad::SFad<double,2> vRestVar( model_.vRest );

      Sacado::Fad::SFad<double,2> resultFad = hEquF( vVar, hVar, vRestVar );
      segHEquFvalue[i] = resultFad.val();
      dhF_dV[i]        = resultFad.dx(0);
      dhF_dh[i]        = resultFad.dx(1);
    }
    {
      const int numDeriv = 1;
      Sacado::Fad::SFad<double,1> hVar( numDeriv, 0, hVarSeg );

      Sacado::Fad::SFad<double,1> resultFad = hEquQ( hVar );
      segHEquQvalue[i] = resultFad.val();
      dhQ_dh[i]        = resultFad.dx(0);
    }

    // calculate potassium current
    potassiumCurrent[i] = model_.gK * pow(nVarSeg, 4.0) * (vSeg - model_.eK);

    // calculate sodium current
    sodiumCurrent[i] = model_.gNa * pow(mVarSeg, 3.0) * hVarSeg * (vSeg - model_.eNa);

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << "Segment " << i << std::endl
                << "vPrev = " << vPrev << std::endl
                << "vSeg = " << vSeg << std::endl
                << "vNext = " << vNext << std::endl
                << "n, m, h = " << nVarSeg << ", " << mVarSeg << ", " << hVarSeg << std::endl
                << "segFvalue =  " << segFvalue[i] << std::endl
                << "segQvalue =  " << segQvalue[i] << std::endl
                << "segNEquFvalue = " << segNEquFvalue[i] << std::endl
                << "segNEquQvalue = " << segNEquQvalue[i] << std::endl
                << "segMEquFvalue = " << segMEquFvalue[i] << std::endl
                << "segMEquQvalue = " << segMEquQvalue[i] << std::endl
                << "segHEquFvalue = " << segHEquFvalue[i] << std::endl
                << "segHEquQvalue = " << segHEquQvalue[i] << std::endl
                << std::endl;

    }
#endif


  }



  return bsuccess;
}
//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState()
{
  bool bsuccess = true;

  updateIntermediateVars ();

  N_LAS_Vector * solVectorPtr = extData.nextSolVectorPtr;
  N_LAS_Vector * staVectorPtr = extData.nextStaVectorPtr;

  for( int i=0; i<nSeg; i++)
  {
    (*staVectorPtr)[li_KCurrentState[i]]  = potassiumCurrent[i];
    (*staVectorPtr)[li_NaCurrentState[i]] = sodiumCurrent[i];
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Instance::updateSecondaryState ()
{
  bool bsuccess = true;
  N_LAS_Vector * staVectorPtr = extData.nextStaVectorPtr;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEQVector
//
// Purpose       : Loads the Q-vector contributions
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Instance::loadDAEQVector ()
{
  bool bsuccess = true;

  N_LAS_Vector * daeQVecPtr = extData.daeQVectorPtr;

  for( int i=0; i<nSeg ; i++)
  {
    (*daeQVecPtr)[li_Vol[i]]  += segQvalue[i];
    (*daeQVecPtr)[li_nPro[i]] += segNEquQvalue[i];
    (*daeQVecPtr)[li_mPro[i]] += segMEquQvalue[i];
    (*daeQVecPtr)[li_hPro[i]] += segHEquQvalue[i];
  }

  return bsuccess;
}



//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 Neuron 3 instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  bool bsuccess=true;

  N_LAS_Vector * daeFVecPtr = extData.daeFVectorPtr;

  (*daeFVecPtr)[li_Pos]  += kcl1Fvalue;
  (*daeFVecPtr)[li_Neg]  += kcl2Fvalue;

  for( int i=0; i<nSeg ; i++)
  {
    (*daeFVecPtr)[li_Vol[i]]  += segFvalue[i];
    (*daeFVecPtr)[li_nPro[i]] += segNEquFvalue[i];
    (*daeFVecPtr)[li_mPro[i]] += segMEquFvalue[i];
    (*daeFVecPtr)[li_hPro[i]] += segHEquFvalue[i];
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx
//
// Purpose       : Loads the Q-vector contributions for a single
//                 Neuron 3 instance.
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  bool bsuccess = true;

  N_LAS_Matrix * dQdxMatPtr = extData.dQdxMatrixPtr;

  for( int i=0; i<nSeg ; i++)
  {
    (*dQdxMatPtr)[li_Vol[i]][SegVEqnVsegOffset[i]] += segQ_dV[i];
    (*dQdxMatPtr)[li_nPro[i]][NEquNNodeOffset[i]]  += dnQ_dn[i];
    (*dQdxMatPtr)[li_mPro[i]][MEquMNodeOffset[i]]  += dmQ_dm[i];
    (*dQdxMatPtr)[li_hPro[i]][HEquHNodeOffset[i]]  += dhQ_dh[i];
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 Neuron 3 instance.
//
// Special Notes : This is an algebraic constaint.
//
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  bool bsuccess = true;

  N_LAS_Matrix * dFdxMatPtr = extData.dFdxMatrixPtr;

  (*dFdxMatPtr)[li_Pos][APosEquPosNodeOffset] += dkcl1F_dVin;
  (*dFdxMatPtr)[li_Pos][APosEquNextNodeOffset] += dkcl1F_dVs0;

  (*dFdxMatPtr)[li_Neg][ANegEquNegNodeOffset] += dkcl2F_dVout;
  (*dFdxMatPtr)[li_Neg][ANegEquLastNodeOffset] += dkcl2F_dVsn;

  for( int i=0; i<nSeg ; i++)
  {
    (*dFdxMatPtr)[li_Vol[i]][SegVEqnVpreOffset[i]] += segF_dVp[i];
    (*dFdxMatPtr)[li_Vol[i]][SegVEqnVsegOffset[i]] += segF_dV[i];
    (*dFdxMatPtr)[li_Vol[i]][SegVEqnVnexOffset[i]] += segF_dVn[i];
    (*dFdxMatPtr)[li_Vol[i]][SegVEqnNOffset[i]]    += segF_dn[i];
    (*dFdxMatPtr)[li_Vol[i]][SegVEqnMOffset[i]]    += segF_dm[i];
    (*dFdxMatPtr)[li_Vol[i]][SegVEqnHOffset[i]]    += segF_dh[i];

    (*dFdxMatPtr)[li_nPro[i]][NEquVNodeOffset[i]]  += dnF_dV[i];
    (*dFdxMatPtr)[li_nPro[i]][NEquNNodeOffset[i]]  += dnF_dn[i];
    (*dFdxMatPtr)[li_mPro[i]][MEquVNodeOffset[i]]  += dmF_dV[i];
    (*dFdxMatPtr)[li_mPro[i]][MEquMNodeOffset[i]]  += dmF_dm[i];
    (*dFdxMatPtr)[li_hPro[i]][HEquVNodeOffset[i]]  += dhF_dV[i];
    (*dFdxMatPtr)[li_hPro[i]][HEquHNodeOffset[i]]  += dhF_dh[i];
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setIC
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Instance::setIC ()
{
  bool bsuccess = true;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::varTypes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
void Instance::varTypes( std::vector<char> & varTypeVec )
{
  //varTypeVec.resize(1);
  //varTypeVec[0] = 'I';
}


//-----------------------------------------------------------------------------
// Function      : Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
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
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
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
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
Model::Model(
  const Configuration & configuration,
  const ModelBlock &    MB,
  const FactoryBlock &  factory_block)
  : DeviceModel(MB, configuration.getModelParameters(), factory_block),
    rInt(0.0),
    radius(0.0),
    length(0.0),
    rIntGiven(false),
    radiusGiven(false),
    lengthGiven(false)
{

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to .model line and constant defaults from metadata:
  setModParams (MB.params);

  // Set any non-constant parameter defaults:
  //if (!given("TNOM"))
  //  tnom = getDeviceOptions().tnom;

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
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
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
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
std::ostream &Model::printOutInstances(std::ostream &os) const
{
  std::vector<Instance*>::const_iterator iter;
  std::vector<Instance*>::const_iterator first = instanceContainer.begin();
  std::vector<Instance*>::const_iterator last  = instanceContainer.end();

  int i, isize;
  isize = instanceContainer.size();

  os << std::endl;
  os << "Number of Neuron instances: " << isize << std::endl;
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


Device *Traits::factory(const Configuration &configuration, const FactoryBlock &factory_block)
{

  return new DeviceMaster<Traits>(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);
}

void registerDevice()
{
  Config<Traits>::addConfiguration()
    .registerDevice("neuron", 3)
    .registerModelType("neuron", 3);
}

} // namespace Neuron3
} // namespace Device
} // namespace Xyce

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
// Filename       : $RCSfile: N_DEV_Neuron6.C,v $
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
// Revision Number: $Revision: 1.37.2.1 $
//
// Revision Date  : $Date: 2013/10/03 17:23:33 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

#include <N_UTL_Misc.h>

// ----------   Xyce Includes   ----------
#include <N_DEV_Neuron6.h>
#include <N_DEV_Neuron_CommonEquations.h>
#include <N_DEV_MembranePassive.h>
#include <N_DEV_MembraneHH.h>
#include <N_DEV_MembraneCS.h>
#include <N_DEV_MembraneUserDefined.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_MatrixLoadData.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>

namespace Xyce {
namespace Device {

template<>
ParametricData<Neuron6::Instance>::ParametricData()
{
  setNumNodes(2);
  setNumOptionalNodes(0);
  setNumFillNodes(0);
  setModelRequired(1);
  setPrimaryParameter("");
  addModelType("NEURON");

  // Set up map for double precision variables:
  addPar ("R", 1.0, false, ParameterType::NO_DEP,
          &Neuron6::Instance::rInt,
          &Neuron6::Instance::rIntGiven,
          U_OHMM, CAT_NONE, "Intracellular resistivity");
  // typical values 1-3 kOhm-mm; use 1 Kohm-mm which is the same as 1 Ohm-m
  addPar ("A", 0.00025, false, ParameterType::NO_DEP,
          &Neuron6::Instance::radius,
          &Neuron6::Instance::radiusGiven,
          U_METER, CAT_NONE, "Segment radius");
  // 250 microns, based on NEURON default diameter of 500 microns
  addPar ("L", 0.0001, false, ParameterType::NO_DEP,
          &Neuron6::Instance::length,
          &Neuron6::Instance::lengthGiven,
          U_METER, CAT_NONE, "Cable length");
  // 100 microns, based on NEURON default length

  // add non-double types
  addPar("N", 1, false, ParameterType::NO_DEP,  &Neuron6::Instance::nSeg ,
          &Neuron6::Instance::nSegGiven ,
         U_NONE, CAT_NONE, "Number of segments" );
}

template<>
ParametricData<Neuron6::Model>::ParametricData()
{
  // Set up map for double precision variables:
  addPar ("CMEM", 0.0, false, ParameterType::NO_DEP,
          &Neuron6::Model::cMem,
          &Neuron6::Model::cMemGiven,
          U_FARADMM2, CAT_NONE, "Membrane capacitance");

  addPar ("GMEM", 0.0, false, ParameterType::NO_DEP,
          &Neuron6::Model::gMem,
          &Neuron6::Model::gMemGiven,
          U_OHMM1MM2, CAT_NONE, "Membrane conductance");

  addPar ("VREST", 0.0, false, ParameterType::NO_DEP,
          &Neuron6::Model::vRest,
          &Neuron6::Model::vRestGiven,
          U_VOLT, CAT_NONE, "Resting potential");

  addPar ("EK", 0.0, false, ParameterType::NO_DEP,
          &Neuron6::Model::eK,
          &Neuron6::Model::eKGiven,
          U_VOLT, CAT_NONE, "Potassium resting potential");

  addPar ("GK", 0.0, false, ParameterType::NO_DEP,
          &Neuron6::Model::gK,
          &Neuron6::Model::gKGiven,
          U_OHMM1MM2, CAT_NONE, "Potassium base conductance");

  addPar ("ENA", 0.0, false, ParameterType::NO_DEP,
          &Neuron6::Model::eNa,
          &Neuron6::Model::eNaGiven,
          U_VOLT, CAT_NONE, "Sodium resting potential");

  addPar ("GNA", 0.0, false, ParameterType::NO_DEP,
          &Neuron6::Model::gNa,
          &Neuron6::Model::gNaGiven,
          U_OHMM1MM2, CAT_NONE, "Sodium base conductance");

  addPar ("EA", 0.0, false, ParameterType::NO_DEP,
          &Neuron6::Model::eA,
          &Neuron6::Model::eAGiven,
          U_VOLT, CAT_NONE, "a-current rest potential");

  addPar ("GA",0.0, false, ParameterType::NO_DEP,
          &Neuron6::Model::gA,
          &Neuron6::Model::gAGiven,
          U_OHMM1MM2, CAT_NONE, "a-current base conductance");

  addPar ("ECA", 0.0, false, ParameterType::NO_DEP,
          &Neuron6::Model::eCa,
          &Neuron6::Model::eCaGiven,
          U_VOLT, CAT_NONE, "Calcium rest potential");

  addPar ("GCA", 0.0, false, ParameterType::NO_DEP,
          &Neuron6::Model::gCa,
          &Neuron6::Model::gCaGiven,
          U_OHMM1MM2, CAT_NONE, "Calcium base conductance");

  addPar ("EKCA", 0.0, false, ParameterType::NO_DEP,
          &Neuron6::Model::eKCa,
          &Neuron6::Model::eKCaGiven,
          U_VOLT, CAT_NONE, "Potassium-calcium rest potential");

  addPar ("GKCA", 0.0, false, ParameterType::NO_DEP,
          &Neuron6::Model::gKCa,
          &Neuron6::Model::gKCaGiven,
          U_OHMM1MM2, CAT_NONE, "Potassium-calcium base conductance");

  addPar ("CAINIT", 0.0, false, ParameterType::NO_DEP,
          &Neuron6::Model::CaInit,
          &Neuron6::Model::CaInitGiven,
          U_MOLAR, CAT_NONE, "initial intra-cellular calcium concentration");

  addPar ("CAGAMMA", 0.0, false, ParameterType::NO_DEP,
          &Neuron6::Model::CaGamma,
          &Neuron6::Model::CaGammaGiven,
          U_NONE, CAT_NONE, "calcium current to concentration multiplier");

  addPar ("CATAU", 0.0, false, ParameterType::NO_DEP,
          &Neuron6::Model::CaTau,
          &Neuron6::Model::CaTauGiven,
          U_SECOND, CAT_NONE, "calcium removal time constant");

  addPar ("R", 1.0, false, ParameterType::NO_DEP,
          &Neuron6::Model::rInt,
          &Neuron6::Model::rIntGiven,
          U_OHMM, CAT_NONE, "Intracellular resistivity");
  // typical values 1-3 kOhm-mm; use 1 Kohm-mm which is the same as 1 Ohm-m

  addPar ("A", 0.00025, false, ParameterType::NO_DEP,
          &Neuron6::Model::radius,
          &Neuron6::Model::radiusGiven,
          U_METER, CAT_NONE, "Segment radius");
  // 250 microns, based on NEURON default diameter of 500 microns

  addPar ("L", 0.0001, false, ParameterType::NO_DEP,
          &Neuron6::Model::length,
          &Neuron6::Model::lengthGiven,
          U_METER, CAT_NONE, "Cable length");
  // 100 microns, based on NEURON default length

  addPar ("I", 0.0, false, ParameterType::SOLN_DEP,
          &Neuron6::Model::I,
          NULL, U_AMP, CAT_NONE, "Current for user-defined current equation");


  // add non-double types
  addPar("IONCHANNELMODEL", "", false, ParameterType::NO_DEP,  &Neuron6::Model::ionChannelModel ,
          &Neuron6::Model::ionChannelModelGiven , U_NONE, CAT_NONE, "Neuron6::Model to use for ion channels" );
  addPar("N", 1, false, ParameterType::NO_DEP,  &Neuron6::Model::nSeg ,
          &Neuron6::Model::nSegGiven , U_NONE, CAT_NONE, "Number of segments" );

  addPar("MM_CURRENT", vector<string>(), false, ParameterType::NO_DEP,  &Neuron6::Model::membraneCurrentEqus,
          & Neuron6::Model::membraneCurrentEqusGiven, U_NONE, CAT_NONE, "Contribution to membrane current" );
  addPar("MM_INDVARS", vector<string>(), false, ParameterType::NO_DEP,  &Neuron6::Model::membraneIndpVars,
          & Neuron6::Model::membraneIndpVarsGiven, U_NONE, CAT_NONE, "Independant variables for ion channel equations" );
  addPar("MM_INDFEQUS", vector<string>(), false, ParameterType::NO_DEP,  &Neuron6::Model::membraneIndpFEqus,
          & Neuron6::Model::membraneIndpFEqusGiven, U_NONE, CAT_NONE, "Independant variables: F equations" );
  addPar("MM_INDQEQUS", vector<string>(), false, ParameterType::NO_DEP,  &Neuron6::Model::membraneIndpQEqus,
          & Neuron6::Model::membraneIndpQEqusGiven, U_NONE, CAT_NONE, "Independant variables: Q equations" );
  addPar("MM_FUNCTIONS", vector<string>(), false, ParameterType::NO_DEP,  &Neuron6::Model::membraneFunctions,
          & Neuron6::Model::membraneFunctionsGiven, U_NONE, CAT_NONE, "Functions for membrane Neuron6::Model" );
  addPar("MM_PARAMETERS", vector<string>(), false, ParameterType::NO_DEP,  &Neuron6::Model::membraneParameters,
          & Neuron6::Model::membraneParametersGiven, U_NONE, CAT_NONE, "Parameters for membrane Neuron6::Model" );
}

namespace Neuron6 {



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
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
Instance::Instance(InstanceBlock & IB,
                   Model & Miter,
                   MatrixLoadData & mlData1,
                   SolverState &ss1,
                   ExternData  &ed1,
                   DeviceOptions & do1)
  : DeviceInstance (IB, mlData1, ss1, ed1, do1),
    model_(Miter),
    rInt(0.0),
    radius(0.0),
    length(0.0),
    segArea(0.0),
    gSeg(0.0),
    nSeg(0),
    rIntGiven(false),
    radiusGiven(false),
    lengthGiven(false),
    nSegGiven(false),
    numIntVarsPerSegment(0),
    numStateVarsPerSegment(0),
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
  setName(IB.getName());
  setModelName(model_.getName());


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

  // if nSeg is still unknown then estimate it via lambda-d rule
  if (!nSegGiven)
  {
    // equations from The Neuron Book, pgs 122-123
    //     says d_lambda value of 0.1 is adequate for most purposes;
    //     use a smaller value if membrane time constant (tau_m) is shorter than 8 ms
    //     but generally uses directly specify nSeg if the default rule doesn't apply
    double d_lambda = 0.1;
    double f = 100.0;	// frequency
    // In equations from Neuron book, C was given in uF; we use F
    double cMem = model_.cMem * 1.0e6;
    // NEURON version of lambda_f equation took d in microns, other
    // distances in cm, and returned lambda_f in microns:
    //   lambda_f = 1.0e5 * sqrt(2*radius/(4*M_PI*f*rInt*cMem));
    //   nSeg = int((length/(d_lambda*lambda_f)+0.9)/2)*2 + 1;
    // I modified the coefficient in the lambda_f equation
    // to keep distance units in cm
    // (should also work in other units as long as the units are consistent)
    double lambda_f = 1.0e3 * sqrt(2*radius/(4*M_PI*f*rInt*cMem));
    nSeg = int((length/(d_lambda*lambda_f)+0.9)/2)*2 + 1;
  }

  // now that nSeg, length and radius are set calculate segment area
  segArea = 2.0 * M_PI * radius * length / nSeg;

  // now we can calculate the number of internal vars
  // by default we'll have a voltage at each node (nSeg)
  // and no state vars.  If the user has an ion channel on, then we'll add in vars for that
  // ask the membrane model for the number of vars it has.
  numIntVarsPerSegment = model_.membraneModel_->numVars();
  numStateVarsPerSegment = 0;

  /*
    if( model_.ConnorStevensOn_ ) {
    numIntVarsPerSegment += 9;
    numStateVarsPerSegment += 2;   // two currents per segment
    }
  */

  numIntVars = numIntVarsPerSegment*nSeg;
  numStateVars = numStateVarsPerSegment*nSeg;

  // total up number of vars.
  int numVars = numExtVars + numIntVars;

  //
  // i(n) - I(n)/A - g(n,n+1) * (V(n+1) - V(n)) - g(n,n-1) * (V(n-1) - V(n)) + Cm dV(n)/d(t) = 0 LHS needs leak term memG ( vSeg - vRest )
  //
  // Vin   : i(n) - I(n)/A - g(n,n+1) * (V(n+1) - V(n)) + Cm dV(n)/d(t) = 0
  // Vout  : i(n) - I(n)/A - g(n,n-1) * (V(n-1) - V(n)) + Cm dV(n)/d(t) = 0
  // Vnode : i(n) - I(n)/A - g(n,n+1) * (V(n+1) - V(n)) - g(n,n-1) * (V(n-1) - V(n)) + Cm dV(n)/d(t) = 0
  //        plus node supporting equations (a, b, m)
  //
  // jacobian format for just the passive cable
  //             Vin      Vout      V1        V(nSeg)
  // kcl Vin    yes               g(0,1)
  // kcl Vout             yes         g(n,n-1)
  // kcl V1     yes                yes        yes
  //

  // if the membrane model includes additional internal variables, the above changes to something
  // along these lines:
  //             Vin      Vout      V1     x1     y1    V(nSeg)     x(nSeg)     y(nSeg)
  // kcl Vin    yes               g(0,1)
  // kcl Vout             yes         g(n,n-1)
  // kcl V1     yes                yes     yes   yes     yes
  // note that V1's dependence on internal variables for segment 1 comes before dependence on
  // Vnext, except for the last segment, in which case Vnext is Vout

  // set up jacStamp.  This is dependant on the membrane model.  The only part this
  // constructor really knows about is the external variables Vin and vOut and internal node
  // voltages

  if( jacStamp.empty() )       // redundant as jacStamp is not static for this device
  {                            // it can't be as each cable may have a different number of nodes
    jacStamp.resize(numVars);
    jacStamp[0].resize(2);
    jacStamp[0][0] = 0;                               // Vin
    jacStamp[0][1] = 2;                               // Vseg[0]
    jacStamp[1].resize(2);
    jacStamp[1][0] = 1;                               // Vout
    jacStamp[1][1] = numVars - numIntVarsPerSegment;  // VnSeg[nSeg-1]

    // now loop over each segment and have the membrane model fill that instanceContainer
    for( int i=0; i<nSeg; i++ )
    {
      // the cable model should take care of the Vpre, V, Vnext dependence as that
      // is part of the cable equation.  Let the membraneModel_ handle what happens
      // at the membrane level
      int offset = numExtVars + i * numIntVarsPerSegment;	// row for vSeg
      jacStamp[offset].resize( numIntVarsPerSegment + 2 );   // + 2 for Vin and Vout

      // have to handle a few special cases which can alter the ordering of Vprev, Vseg, Vnext
      // for each segment number, save the offsets for the previous, current, and next segment so
      // we don't have to rethink these special cases every time.
      if( nSeg == 1 )
      {
        // here the ordering is Vin  Vout  Vseg
        prevMap[i] = 0;
        nextMap[i] = 1;
        segMap[i] = 2;
        jacStamp[offset][prevMap[i]] = 0;                                     // Vin
        jacStamp[offset][nextMap[i]] = 1;                                     // Vout
        jacStamp[offset][segMap[i]] = 2;                                      // Vseg

      }
      else if( i==0 )
      {
        // ordering is Vin Vseg (seg int vars)  Vnext
        prevMap[i] = 0;
        segMap[i] = 1;
        nextMap[i] = numIntVarsPerSegment+1;
        jacStamp[offset][prevMap[i]] = 0;                                     // Vin
        jacStamp[offset][segMap[i]] = offset;                                 // Vseg
        jacStamp[offset][nextMap[i]] = offset + numIntVarsPerSegment;         // Vnext
      }
      else if( i==(nSeg-1) )
      {
        // ordering is Vout Vprev Vseg
        nextMap[i] = 0;
        prevMap[i] = 1;
        segMap[i] = 2;
        jacStamp[offset][nextMap[i]] = 1;                                     // Vout
        jacStamp[offset][prevMap[i]] = offset - numIntVarsPerSegment;         // Vprev
        jacStamp[offset][segMap[i]] = offset;                                 // Vseg
      }
      else
      {
        // ordering is Vprev Vseg (seg int vars) Vnext
        prevMap[i] = 0;
        segMap[i] = 1;
        nextMap[i] = numIntVarsPerSegment+1;
        jacStamp[offset][prevMap[i]] = offset - numIntVarsPerSegment;         // Vprev
        jacStamp[offset][segMap[i]] = offset;                                 // Vseg
        jacStamp[offset][nextMap[i]] = offset + numIntVarsPerSegment;         // Vnext
      }

      // pass the membraneModel_ enough information for it to construct its part of the jacobian
      model_.membraneModel_->setJacStamp( numExtVars, i, segMap[i], jacStamp );
    }

  }

  /*
  // print out jacStamp
  std::cout << "jacStamp for Neuron6" << std::endl;
  int numRows = jacStamp.size();
  for( int i=0; i< numRows; i++ )
  {
  int numCol = jacStamp[i].size();
  std::cout << "jacStamp[ " << i << " ] = { ";
  for(int j=0; j<numCol; j++)
  {
  std::cout << jacStamp[i][j] << "  ";
  }
  std::cout << " } " <<  std::endl;
  }
  std::cout << std::endl;
  */



  // calculate segment conductivities used in load calls:
  // equation 6.30 Theoretical neuroscience: computational and mathematical modeling of neural systems, Peter Dayan and L. F. Abbot 2001
  //   g(n,n') = radius * (radius')^2 / ( rInt segLength * ( segLength * (radius')^2 + segLength' * (radius)^2 ) )
  // this equation is in terms of current/area, in which case conductivity to the previous and next
  // segment is not symmetric if the segment length and/or radius is not are not equal
  // But since we work in current rather than current density, this is not a problem
  double segLength = length / nSeg;
  double rLong = rInt * segLength / (M_PI * radius * radius);  // longitudinal resistance (ohm)

  // watch out for divide-by-0 cases; if resistivity is close to 0, just set conductance to some large value
  // TODO:  Really should be testing for anything close enough to zero to cause problems.  Is there a better 'large' value to use?
  if (rLong == 0.0)
  {
    gSeg = 1000.0;
  }
  else
  {
    gSeg = 1.0 / rLong;
  }

  // variable indecies loads
  if( model_.ConnorStevensOn_ )
  {
    li_Vol.resize(nSeg);
    li_nPro.resize(nSeg);
    li_mPro.resize(nSeg);
    li_hPro.resize(nSeg);
    li_aPro.resize(nSeg);
    li_bPro.resize(nSeg);
    li_MPro.resize(nSeg);
    li_HPro.resize(nSeg);
    li_cPro.resize(nSeg);
    li_CaPro.resize(nSeg);
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
    segAEquFvalue.resize(nSeg);
    segAEquQvalue.resize(nSeg);
    segBEquFvalue.resize(nSeg);
    segBEquQvalue.resize(nSeg);
    segM_EquFvalue.resize(nSeg);
    segM_EquQvalue.resize(nSeg);
    segH_EquFvalue.resize(nSeg);
    segH_EquQvalue.resize(nSeg);
    segCEquFvalue.resize(nSeg);
    segCEquQvalue.resize(nSeg);
    segCaEquFvalue.resize(nSeg);
    segCaEquQvalue.resize(nSeg);

    // jacobian elements
    segF_dVp.resize(nSeg);
    segF_dV.resize(nSeg);
    segF_dVn.resize(nSeg);
    segF_dn.resize(nSeg);
    segF_dm.resize(nSeg);
    segF_dh.resize(nSeg);
    segF_da.resize(nSeg);
    segF_db.resize(nSeg);
    segF_dM.resize(nSeg);
    segF_dH.resize(nSeg);
    segF_dc.resize(nSeg);
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
    daF_dV.resize(nSeg);
    daF_da.resize(nSeg);
    daQ_da.resize(nSeg);
    dbF_dV.resize(nSeg);
    dbF_db.resize(nSeg);
    dbQ_db.resize(nSeg);
    dMF_dV.resize(nSeg);
    dMF_dM.resize(nSeg);
    dMQ_dM.resize(nSeg);
    dHF_dV.resize(nSeg);
    dHF_dH.resize(nSeg);
    dHQ_dH.resize(nSeg);
    dcF_dV.resize(nSeg);
    dcF_dc.resize(nSeg);
    dcF_dCa.resize(nSeg);
    dcQ_dc.resize(nSeg);
    dCaF_dV.resize(nSeg);
    dCaF_dM.resize(nSeg);
    dCaF_dH.resize(nSeg);
    dCaF_dCa.resize(nSeg);
    dCaQ_dCa.resize(nSeg);

    // state vars
    potassiumCurrent.resize(nSeg);
    sodiumCurrent.resize(nSeg);
  }
  else
  {
    /*
      li_Vol.resize(nSeg);
      segFvalue.resize(nSeg);
      segQvalue.resize(nSeg);

      // jacobian elements
      segF_dVp.resize(nSeg);
      segF_dV.resize(nSeg);
      segF_dVn.resize(nSeg);
      segQ_dV.resize(nSeg);
    */
  }

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
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Instance::processParams(string param)
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
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
void Instance::registerLIDs(const vector<int> & intLIDVecRef,
                            const vector<int> & extLIDVecRef)
{
  string msg;

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline =
    "-------------------------------------------------------------------------"
    "----";
  if (getDeviceOptions().debugLevel > 0)
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

  li_Pos  = extLIDVec[0];
  li_Neg  = extLIDVec[1];

  // resize our storage location for the internal vars.
  li_internalVars.resize( numIntVars );

  // now copy in the local ids
  for( int i=0; i<numIntVars; i++ )
  {
    li_internalVars[i] = intLIDVec[i];
  }
  /*
    if( model_.ConnorStevensOn_ )
    {
    for( int i=0, j=0; i<nSeg; i++, j+=10)
    {
    li_Vol[i]  = intLIDVec[j];
    li_nPro[i] = intLIDVec[j+1];
    li_mPro[i] = intLIDVec[j+2];
    li_hPro[i] = intLIDVec[j+3];
    li_aPro[i] = intLIDVec[j+4];
    li_bPro[i] = intLIDVec[j+5];
    li_MPro[i] = intLIDVec[j+6];
    li_HPro[i] = intLIDVec[j+7];
    li_cPro[i] = intLIDVec[j+8];
    li_CaPro[i] = intLIDVec[j+9];
    }
    }
    else
    {
    for( int i=0, j=0; i<nSeg; i++, j+=1)
    {
    li_Vol[i]  = intLIDVec[j];
    }
    }

    #ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 )
    {
    cout << "  li_Pos = " << li_Pos << endl
    << "  li_Neg = " << li_Neg << endl;
    for( int i=0; i<nSeg; i++ )
    {
    cout << "  li_Vol[ " << i << " ] = " << li_Vol[i] << endl;

    if( model_.ConnorStevensOn_ )
    {
    std::cout << "  li_nPro[ " << i << " ] = " << li_nPro[i] << endl
    << "  li_mPro[ " << i << " ] = " << li_mPro[i] << endl
    << "  li_hPro[ " << i << " ] = " << li_hPro[i] << endl
    << "  li_aPro[ " << i << " ] = " << li_aPro[i] << endl
    << "  li_bPro[ " << i << " ] = " << li_bPro[i] << endl
    << "  li_MPro[ " << i << " ] = " << li_MPro[i] << endl
    << "  li_HPro[ " << i << " ] = " << li_HPro[i] << endl
    << "  li_cPro[ " << i << " ] = " << li_cPro[i] << endl
    << "  li_CaPro[ " << i << " ] = " << li_CaPro[i] << endl;
    }
    }
    }
    #endif

    #ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 )
    {
    cout << dashedline << endl;
    }
    #endif
  */
}

//-----------------------------------------------------------------------------
// Function      : Instance::getIntNameMap
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
map<int,string> & Instance::getIntNameMap ()
{
  // set up the internal name map, if it hasn't been already.
  if (intNameMap.empty ())
  {
    string tmpstr;
    for( int i=0; i<nSeg; i++)
    {
      ostringstream segNumber;
      segNumber << i;
      string segNumStr = segNumber.str();

      tmpstr = getName() + "_" + "V" + segNumStr;
      spiceInternalName (tmpstr);
      intNameMap[ li_internalVars[i*numIntVarsPerSegment] ] = tmpstr;

      if( model_.ConnorStevensOn_ )
      {
        tmpstr = getName() + "_" + "N" + segNumStr;
        spiceInternalName (tmpstr);
        intNameMap[ li_nPro[i] ] = tmpstr;

        tmpstr = getName() + "_" + "M" + segNumStr;
        spiceInternalName (tmpstr);
        intNameMap[ li_mPro[i] ] = tmpstr;

        tmpstr = getName() + "_" + "H" + segNumStr;
        spiceInternalName (tmpstr);
        intNameMap[ li_hPro[i] ] = tmpstr;

        tmpstr = getName() + "_" + "A" + segNumStr;
        spiceInternalName (tmpstr);
        intNameMap[ li_aPro[i] ] = tmpstr;

        tmpstr = getName() + "_" + "B" + segNumStr;
        spiceInternalName (tmpstr);
        intNameMap[ li_bPro[i] ] = tmpstr;

        tmpstr = getName() + "_" + "M_" + segNumStr;
        spiceInternalName (tmpstr);
        intNameMap[ li_MPro[i] ] = tmpstr;

        tmpstr = getName() + "_" + "H_" + segNumStr;
        spiceInternalName (tmpstr);
        intNameMap[ li_HPro[i] ] = tmpstr;

        tmpstr = getName() + "_" + "C" + segNumStr;
        spiceInternalName (tmpstr);
        intNameMap[ li_cPro[i] ] = tmpstr;

        tmpstr = getName() + "_" + "Ca" + segNumStr;
        spiceInternalName (tmpstr);
        intNameMap[ li_CaPro[i] ] = tmpstr;
      }

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

  if( model_.ConnorStevensOn_ )
  {
    for( int i=0, j=0; i<nSeg; i++, j+=2)
    {
      li_KCurrentState[i] = staLIDVec[j];
      li_NaCurrentState[i] = staLIDVec[j+1];
    }
  }
  else
  {
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

//   std::cout << "Masking n, m and h" << std::endl;
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
const vector< vector<int> > & Instance::jacobianStamp() const
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
void Instance::registerJacLIDs( const vector< vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );

  // resize our storage location and store the results

  int numRows = jacLIDVec.size();
  jacobianOffsets.resize( numRows );
  for( int i=0; i< numRows; i++ )
  {
    int numCol = jacLIDVec[i].size();
    jacobianOffsets[i].resize( numCol );
    for( int j=0; j< numCol; j++ )
    {
      jacobianOffsets[i][j] = jacLIDVec[i][j];
    }
  }

  // external terminals
  APosEquPosNodeOffset = jacLIDVec[0][0];
  APosEquNextNodeOffset = jacLIDVec[0][1];
  ANegEquNegNodeOffset = jacLIDVec[1][0];
  ANegEquLastNodeOffset = jacLIDVec[1][1];
  /*
    std::cout << "APosEquPosNodeOffset = " << APosEquPosNodeOffset << std::endl;
    std::cout << "APosEquNextNodeOffset = " << APosEquNextNodeOffset << std::endl;
    std::cout << "ANegEquNegNodeOffset = " << ANegEquNegNodeOffset << std::endl;
    std::cout << "ANegEquLastNodeOffset = " << ANegEquLastNodeOffset << std::endl;
  */

  if( model_.ConnorStevensOn_ )
  {
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
    AEquVNodeOffset.resize(nSeg);
    AEquANodeOffset.resize(nSeg);
    BEquVNodeOffset.resize(nSeg);
    BEquBNodeOffset.resize(nSeg);
    M_EquVNodeOffset.resize(nSeg);
    M_EquM_NodeOffset.resize(nSeg);
    H_EquVNodeOffset.resize(nSeg);
    H_EquH_NodeOffset.resize(nSeg);
    CEquVNodeOffset.resize(nSeg);
    CEquCNodeOffset.resize(nSeg);
    CEquCaNodeOffset.resize(nSeg);
    CaEquVNodeOffset.resize(nSeg);
    CaEquM_NodeOffset.resize(nSeg);
    CaEquH_NodeOffset.resize(nSeg);
    CaEquCaNodeOffset.resize(nSeg);

    for(int i=0, j=2; i<nSeg; i++, j+=10 )
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
      AEquVNodeOffset[i] = jacLIDVec[j+4][0];
      AEquANodeOffset[i] = jacLIDVec[j+4][1];
      BEquVNodeOffset[i] = jacLIDVec[j+5][0];
      BEquBNodeOffset[i] = jacLIDVec[j+5][1];
      M_EquVNodeOffset[i] = jacLIDVec[j+6][0];
      M_EquM_NodeOffset[i] = jacLIDVec[j+6][1];
      H_EquVNodeOffset[i] = jacLIDVec[j+7][0];
      H_EquH_NodeOffset[i] = jacLIDVec[j+7][1];
      CEquVNodeOffset[i] = jacLIDVec[j+8][0];
      CEquCNodeOffset[i] = jacLIDVec[j+8][1];
      CEquCaNodeOffset[i] = jacLIDVec[j+8][2];
      CaEquVNodeOffset[i] = jacLIDVec[j+9][0];
      CaEquM_NodeOffset[i] = jacLIDVec[j+9][1];
      CaEquH_NodeOffset[i] = jacLIDVec[j+9][2];
      CaEquCaNodeOffset[i] = jacLIDVec[j+9][3];

      /*
        std::cout <<  SegVEqnVpreOffset[i] << ", "
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
  else
  {
    // internal variables
    SegVEqnVpreOffset.resize(nSeg);
    SegVEqnVsegOffset.resize(nSeg);
    SegVEqnVnexOffset.resize(nSeg);

    for(int i=0, j=2; i<nSeg; i++, j+=numIntVarsPerSegment )
    {
      // std::cout << " i = " << i << " j = " << j << " jacLIDVec[ " << j << " ].size() = " << jacLIDVec[j].size() << std::endl;
      SegVEqnVpreOffset[i] = jacLIDVec[j][0];
      SegVEqnVsegOffset[i] = jacLIDVec[j][1];
      SegVEqnVnexOffset[i] = jacLIDVec[j][2];

      /*
        std::cout <<  SegVEqnVpreOffset[i] << ", "
        << SegVEqnVsegOffset[i] << ", "
        << SegVEqnVnexOffset[i] << std::endl;
      */
    }
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
bool Instance::updateIntermediateVars ()
{
  bool bsuccess = true;

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
bool Instance::updatePrimaryState ()
{
  bool bsuccess = true;

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

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEQVector
//
// Purpose       : Loads the Q-vector contributions for a single
//                 Neuron 6 instance.
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

  N_LAS_Vector * solVectorPtr = extData.nextSolVectorPtr;
  N_LAS_Vector * daeQVecPtr = extData.daeQVectorPtr;

  // no Q component for the cable component of this devcie

  // now let the membrane model load it's Q component (memCap dV/dt, etc.)
  for( int i=0; i<nSeg ; i++)
  {
    model_.membraneModel_->loadDAEQVector ( i, li_internalVars, solVectorPtr, daeQVecPtr, segArea );
    /*
      (*daeQVecPtr)[li_Vol[i]]  += segQvalue[i];
      if( model_.ConnorStevensOn_ )
      {
      (*daeQVecPtr)[li_nPro[i]] += segNEquQvalue[i];
      (*daeQVecPtr)[li_mPro[i]] += segMEquQvalue[i];
      (*daeQVecPtr)[li_hPro[i]] += segHEquQvalue[i];
      (*daeQVecPtr)[li_aPro[i]] += segAEquQvalue[i];
      (*daeQVecPtr)[li_bPro[i]] += segBEquQvalue[i];
      (*daeQVecPtr)[li_MPro[i]] += segM_EquQvalue[i];
      (*daeQVecPtr)[li_HPro[i]] += segH_EquQvalue[i];
      (*daeQVecPtr)[li_cPro[i]] += segCEquQvalue[i];
      (*daeQVecPtr)[li_CaPro[i]] += segCaEquQvalue[i];
      }
    */
  }

  return bsuccess;
}



//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 Neuron 6 instance.
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

  N_LAS_Vector * solVectorPtr = extData.nextSolVectorPtr;
  N_LAS_Vector * daeFVecPtr = extData.daeFVectorPtr;

  double vIn = (*solVectorPtr)[li_Pos];
  double vOut = (*solVectorPtr)[li_Neg];

  // take care of the input and output nodes as they are different
  (*daeFVecPtr)[li_Pos]  += -2.0 * gSeg * ((*solVectorPtr)[li_internalVars[0]] - vIn );
  (*daeFVecPtr)[li_Neg]  += -2.0 * gSeg * ((*solVectorPtr)[li_internalVars[(nSeg-1)*numIntVarsPerSegment]] - vOut );

  for( int i=0; i<nSeg ; i++)
  {
    // for this segment get the values of the local vars
    double vSeg  = (*solVectorPtr)[li_internalVars[i*numIntVarsPerSegment]];
    double vNext = 0.0;
    double gNext = 0.0;
    if (i == (nSeg - 1))
    {
      vNext = vOut;
      gNext = gSeg * 2.0;
    }
    else
    {
      vNext = (*solVectorPtr)[li_internalVars[(i+1)*numIntVarsPerSegment]];
      gNext = gSeg;
    }
    double vPrev = 0.0;
    double gPrev = 0.0;
    if (i == 0 )
    {
      vPrev = vIn;
      gPrev = gSeg * 2.0;
    }
    else
    {
      vPrev = (*solVectorPtr)[li_internalVars[(i-1)*numIntVarsPerSegment]];
      gPrev = gSeg;
    }

    // li_internalVars lists numIntVarsPerSegment for each segment.  V for each segment will
    // be first.  So, V's offset in li_internalVars is i * numIntVarsPerSegment

    (*daeFVecPtr)[li_internalVars[i*numIntVarsPerSegment]] += - gPrev * (vPrev - vSeg) - gNext * (vNext - vSeg);

    model_.membraneModel_->loadDAEFVector ( i, li_internalVars, solVectorPtr, daeFVecPtr, segArea );

  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx
//
// Purpose       : Loads the Q-vector contributions for a single
//                 Neuron 6 instance.
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  bool bsuccess = true;

  N_LAS_Vector * solVectorPtr = extData.nextSolVectorPtr;
  N_LAS_Matrix * dQdxMatPtr = extData.dQdxMatrixPtr;

  for( int i=0; i<nSeg ; i++)
  {
    // let the membrane model load it's part
    model_.membraneModel_->loadDAEdQdx ( i, segMap[i], li_internalVars, jacobianOffsets, solVectorPtr, dQdxMatPtr, segArea );

    /*
      if( model_.ConnorStevensOn_ )
      {
      (*dQdxMatPtr)[li_nPro[i]][NEquNNodeOffset[i]]  += dnQ_dn[i];
      (*dQdxMatPtr)[li_mPro[i]][MEquMNodeOffset[i]]  += dmQ_dm[i];
      (*dQdxMatPtr)[li_hPro[i]][HEquHNodeOffset[i]]  += dhQ_dh[i];
      (*dQdxMatPtr)[li_aPro[i]][AEquANodeOffset[i]] += daQ_da[i];
      (*dQdxMatPtr)[li_bPro[i]][BEquBNodeOffset[i]] += dbQ_db[i];
      (*dQdxMatPtr)[li_MPro[i]][M_EquM_NodeOffset[i]] += dMQ_dM[i];
      (*dQdxMatPtr)[li_HPro[i]][H_EquH_NodeOffset[i]] += dHQ_dH[i];
      (*dQdxMatPtr)[li_cPro[i]][CEquCNodeOffset[i]] += dcQ_dc[i];
      (*dQdxMatPtr)[li_CaPro[i]][CaEquCaNodeOffset[i]] += dCaQ_dCa[i];
      }
    */
  }

  return bsuccess;
}



//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 Neuron 6 instance.
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

  N_LAS_Vector * solVectorPtr = extData.nextSolVectorPtr;
  N_LAS_Matrix * dFdxMatPtr = extData.dFdxMatrixPtr;

  (*dFdxMatPtr)[li_Pos][APosEquPosNodeOffset]  +=  2.0 * gSeg;
  (*dFdxMatPtr)[li_Pos][APosEquNextNodeOffset] += -2.0 * gSeg;

  (*dFdxMatPtr)[li_Neg][ANegEquNegNodeOffset]  +=  2.0 * gSeg;
  (*dFdxMatPtr)[li_Neg][ANegEquLastNodeOffset] += -2.0 * gSeg;

  for( int i=0; i<nSeg ; i++)
  {
    int offset = i * numIntVarsPerSegment;

    int row = numExtVars + i * numIntVarsPerSegment;

    double gPrev = gSeg;
    double gNext = gSeg;
    if (i == 0)
    {
      gPrev = gSeg * 2.0;
    }
    if (i == (nSeg-1))
    {
      gNext = gSeg * 2.0;
    }

    (*dFdxMatPtr)[li_internalVars[offset]][jacobianOffsets[row][prevMap[i]]]           +=  -gPrev;         // Vpre
    (*dFdxMatPtr)[li_internalVars[offset]][jacobianOffsets[row][segMap[i]]]     +=  gPrev + gNext;  // Vseg
    (*dFdxMatPtr)[li_internalVars[offset]][jacobianOffsets[row][nextMap[i]]] +=  -gNext;         // Vnext

    // now let the membrane model load it's part
    model_.membraneModel_->loadDAEdFdx ( i, segMap[i], li_internalVars, jacobianOffsets, solVectorPtr, dFdxMatPtr, segArea );

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
void Instance::varTypes( vector<char> & varTypeVec )
{
  //varTypeVec.resize(1);
  //varTypeVec[0] = 'I';
}


//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : block constructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
Model::Model (const ModelBlock & MB,
              SolverState & ss1,
              DeviceOptions & do1)
  : DeviceModel(MB,ss1,do1),
    cMem(0.0),
    gMem(0.0),
    vRest(0.0),
    eNa(0.0),
    gNa(0.0),
    eK(0.0),
    gK(0.0),
    eA(0.0),
    gA(0.0),
    eCa(0.0),
    gCa(0.0),
    eKCa(0.0),
    gKCa(0.0),
    CaInit(0.0),
    CaGamma(0.0),
    CaTau(0.0),
    rInt(0.0),
    radius(0.0),
    length(0.0),
    nSeg(0.0),
    rIntGiven(false),
    radiusGiven(false),
    lengthGiven(false),
    nSegGiven(false),
    ionChannelModelGiven(false),
    cMemGiven(false),
    gMemGiven(false),
    vRestGiven(false),
    eNaGiven(false),
    gNaGiven(false),
    eKGiven(false),
    gKGiven(false),
    eAGiven(false),
    gAGiven(false),
    eCaGiven(false),
    gCaGiven(false),
    eKCaGiven(false),
    gKCaGiven(false),
    CaInitGiven(false),
    CaGammaGiven(false),
    CaTauGiven(false),
    hodgenHuxleyOn_(false),
    ConnorStevensOn_(false),
    sodiumOn_(false),
    potassiumOn_(false),
    aCurrentOn_(false),
    calciumOn_(false),
    membraneIndpVarsGiven(false),
    membraneIndpFEqusGiven(false),
    membraneIndpQEqusGiven(false)
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

  // check the specified model type and allocate the needed membrane model
  if( ionChannelModelGiven )
  {
    // should change the case on ionChannelModel and simplify this set of clauses
    if( ionChannelModel == "passive" || ionChannelModel == "PASSIVE" )
    {
      membraneModel_ = rcp( new MembranePassive( solState, cMem, gMem, vRest) );

    }
    else if( ionChannelModel == "hh" || ionChannelModel == "HH" )
    {
      membraneModel_ = rcp( new MembraneHH( solState, cMem, gMem, vRest, eK, gK, eNa, gNa) );
    }
    else if( ionChannelModel == "cs" || ionChannelModel == "CS")
    {
      membraneModel_ = rcp( new MembraneCS( solState ) );
    }
    else if( ionChannelModel == "ud" || ionChannelModel == "UD")
    {
      membraneModel_ = rcp( new MembraneUserDefined( solState, cMem, gMem, vRest,
                                                           membraneCurrentEqus, membraneIndpVars, membraneIndpFEqus, membraneIndpQEqus, membraneFunctions, membraneParameters) );
    }
    else
    {
      // issue error and stop
      string msg("Model: unknown ion channel model given \"");
      msg += ionChannelModel;
      msg += "\"";
      N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_FATAL, msg );
    }

  }
  else
  {
    // no model given so assume passive
    membraneModel_ = rcp( new MembranePassive( solState, cMem, gMem, vRest) );
  }
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
// Function      : Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
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
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
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
// Function      : Model::printOutInstances
// Purpose       : debugging tool.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
std::ostream &Model::printOutInstances(std::ostream &os) const
{
  vector<Instance*>::const_iterator iter;
  vector<Instance*>::const_iterator first = instanceContainer.begin();
  vector<Instance*>::const_iterator last  = instanceContainer.end();

  int i, isize;
  isize = instanceContainer.size();

  os << endl;
  os << "Number of Neuron instances: " << isize << endl;
  os << "    name=\t\tmodelName\tParameters" << endl;
  for (i=0, iter=first; iter!=last; ++iter, ++i)
  {
    os << "  " << i << ": " << (*iter)->getName() << "\t";
    os << (*iter)->getModelName();
    os << endl;
  }

  os << endl;
}

} // namespace Neuron6
} // namespace Device
} // namespace Xyce

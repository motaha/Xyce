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

//-----------------------------------------------------------------------------
// Filename       : $RCSfile$
//
// Purpose        : Neuron classes.
//
// Special Notes  :
//
// Creator        : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date  : 06/10/09
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision$
//
// Revision Date  : $Date$
//
// Current Owner  : $Author$
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Neuron4_h
#define Xyce_N_DEV_Neuron4_h

#include <Sacado.hpp>

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>

#ifdef HAVE_MATH_H
#include <math.h>
#endif

namespace Xyce {
namespace Device {
namespace Neuron4 {

// ---------- Forward Declarations ----------
class Model;

//-----------------------------------------------------------------------------
// Class         : Instance
// Purpose       : This is class refers to a single instance of the
//                 Neuron device.  It has two nodes associated with it, a
//                 positive and a negative node.   See the NeuronInstance
//                 class for a more detailed explanation.
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
class Instance : public DeviceInstance
{
  friend class ParametricData<Instance>;
  friend class Model;

public:
  static ParametricData<Instance> &getParametricData();

  virtual const ParametricData<void> &getMyParametricData() const {
    return getParametricData();
  }

  Instance(InstanceBlock & IB,
           Model & Miter,
           MatrixLoadData & mlData1,
           SolverState &ss1,
           ExternData  &ed1,
           DeviceOptions & do1);

  ~Instance();

private:
  Instance(const Instance &);
  Instance &operator=(const Instance &);

public:
  void registerLIDs( const vector<int> & intLIDVecRef,
                     const vector<int> & extLIDVecRef );
  void registerStateLIDs( const vector<int> & staLIDVecRef );

  map<int,string> & getIntNameMap ();
  bool loadDeviceMask();
  const vector< vector<int> > & jacobianStamp() const;
  void registerJacLIDs( const vector< vector<int> > & jacLIDVec );

  bool processParams (string param = "");
  bool updateTemperature(const double & temp_tmp);

  bool updateIntermediateVars ();
  bool updatePrimaryState ();
  bool updateSecondaryState ();
  bool setIC ();

  void varTypes( vector<char> & varTypeVec );

  // load functions, residual:
  bool loadDAEQVector ();
  bool loadDAEFVector ();

  void auxDAECalculations ();

  // load functions, Jacobian:
  bool loadDAEdQdx ();
  bool loadDAEdFdx ();

private:
  // These functions represents the equations that need to be solved
  // for this device.  Since Xyce loads an F and Q contribution, the
  // equations are broken up into their F and Q components.  Thus there
  // is a kcl1EquF() and kcl1EquQ().  Automatic differentiation will
  // be used to generate all the derivatives of these equations for the
  // dF/dX and dQ/dX loads

  // first we list some utility functions for calculating coefficients
  //
  // These functions expect V to be in milli-volts and then return values that
  // are in 1/ms.  Thus the extra factor's of 1000 here and there

  // potassium current, functions for activator equation
  template <typename ScalarT>
  static ScalarT alphaN( const ScalarT Vin)
  {
    ScalarT vScaled = 1000.0 * Vin;  // convert voltage to milli-volts
    ScalarT r = 1000.0 * (0.02 * (vScaled + 45.7)) / (1.0 - std::exp(-0.1*(vScaled+45.7)));
    // result.  the 1000 factor is to change from 1/ms to 1/s
    return r;
  }

  template <typename ScalarT>
  static ScalarT betaN( const ScalarT Vin)
  {
    ScalarT vScaled = 1000.0 * Vin;  // convert voltage to milli-volts
    ScalarT r = 1000.0 * 0.25 * std::exp( -0.0125 * (vScaled + 55.7));
    // result.  the 1000 factor is to change from 1/ms to 1/s
    return r;
  }

  // sodium current, functions for activator equation
  template <typename ScalarT>
  static ScalarT alphaM( const ScalarT Vin)
  {
    ScalarT vScaled = 1000.0 * Vin;  // convert voltage to milli-volts
    ScalarT r = 1000.0 * (0.38 * (vScaled + 29.7)) / (1.0 - std::exp(-0.1*(vScaled+29.7)));
    // result.  the 1000 factor is to change from 1/ms to 1/s
    return r;
  }

  template <typename ScalarT>
  static ScalarT betaM( const ScalarT Vin)
  {
    ScalarT vScaled = 1000.0 * Vin;  // convert voltage to milli-volts
    ScalarT r = 1000.0 * 15.2 * std::exp( -0.0556 * (vScaled + 54.7));
    // result.  the 1000 factor is to change from 1/ms to 1/s
    return r;
  }

  template <typename ScalarT>
  static ScalarT alphaH( const ScalarT Vin)
  {
    ScalarT vScaled = 1000.0 * Vin;  // convert voltage to milli-volts
    ScalarT r = 1000.0 * 0.266 * std::exp( -0.05 * (vScaled + 48.0));
    // result.  the 1000 factor is to change from 1/ms to 1/s
    return r;
  }

  template <typename ScalarT>
  static ScalarT betaH( const ScalarT Vin)
  {
    ScalarT vScaled = 1000.0 * Vin;  // convert voltage to milli-volts
    ScalarT r = 1000.0 * 3.8 / (1.0 + std::exp(-0.1*(vScaled+18.0)));
    // result.  the 1000 factor is to change from 1/ms to 1/s
    return r;
  }

  // a-current functions
  template <typename ScalarT>
  static ScalarT aInf( const ScalarT Vin)
  {
    ScalarT vScaled = 1000.0 * Vin;  // convert voltage to milli-volts
    ScalarT r = std::pow( ((0.0761 * std::exp(0.0314 * (vScaled+94.22))) / (1.0+std::exp(0.0346*(vScaled+1.17)))), 1.0/3.0);
    return r;
  }

  template <typename ScalarT>
  static ScalarT aTau( const ScalarT Vin)
  {
    ScalarT vScaled = 1000.0 * Vin;  // convert voltage to milli-volts
    ScalarT r = (0.3632 + 1.158 / (1.0 + std::exp(0.0497 * (vScaled + 55.96)))) / 1000.0;
    return r;
  }

  template <typename ScalarT>
  static ScalarT bInf( const ScalarT Vin)
  {
    ScalarT vScaled = 1000.0 * Vin;  // convert voltage to milli-volts
    ScalarT r = std::pow( (1.0 / (1.0 + std::exp(0.0688*(vScaled+53.3)))), 4.0);
    return r;
  }

  template <typename ScalarT>
  static ScalarT bTau( const ScalarT Vin)
  {
    ScalarT vScaled = 1000.0 * Vin;  // convert voltage to milli-volts
    ScalarT r = (1.24 + 2.678 / (1.0 + std::exp(0.0624 * (vScaled + 50.0)))) / 1000.0;
    return r;
  }

  // transient calcium functions
  template <typename ScalarT>
  static ScalarT M_Inf( const ScalarT Vin)
  {
    ScalarT vScaled = 1000.0 * Vin;  // convert voltage to milli-volts
    ScalarT r = 1.0/(1.0 + std::exp(-(vScaled+57)/6.2));
    return r;
  }

  template <typename ScalarT>
  static ScalarT M_Tau( const ScalarT Vin)
  {
    ScalarT vScaled = 1000.0 * Vin;  // convert voltage to milli-volts
    ScalarT r = (0.612 + 1.0/(std::exp(-(vScaled+132)/16.7) + std::exp((vScaled+16.8)/18.2)) ) / 1000.0;
    return r;
  }

  template <typename ScalarT>
  static ScalarT H_Inf( const ScalarT Vin)
  {
    ScalarT vScaled = 1000.0 * Vin;  // convert voltage to milli-volts
    ScalarT r = 1.0 / (1.0 + std::exp((vScaled+81)/4.0));
    return r;
  }

  template <typename ScalarT>
  static ScalarT H_Tau( const ScalarT Vin)
  {
    ScalarT vScaled = 1000.0 * Vin;  // convert voltage to milli-volts
    ScalarT r;
    if( vScaled < -80.0 )
    {
      r = std::exp( (vScaled + 467)/66.6 ) / 1000.0;
    }
    else
    {
      r = ( 28.0 + std::exp(-(vScaled+22.0)/10.5)) / 1000.0;
    }
    return r;
  }

  // Calcium dependent Potassium conductances
  template <typename ScalarT>
  static ScalarT C_Inf( const ScalarT Vin, const ScalarT CaConc)
  {
    ScalarT vScaled = 1000.0 * Vin;  // convert voltage to milli-volts
    ScalarT r = (CaConc / (CaConc + 3.0)) * (1.0 / (1.0 + std::exp(-(vScaled+28.3)/12.6 )));
    return r;
  }

  template <typename ScalarT>
  static ScalarT C_Tau( const ScalarT Vin)
  {
    ScalarT vScaled = 1000.0 * Vin;  // convert voltage to milli-volts
    ScalarT r = (90.3 - 75.1/(1.0 + std::exp(-(vScaled+46)/22.7))) / 1000.0;
    return r;
  }


  // now the device equations
  // KCL equation 1
  template <typename ScalarT>
  static ScalarT kcl1EquF( const ScalarT& VSeg, const ScalarT& VSegP, const ScalarT & VSegN, const ScalarT& n, const ScalarT& m, const ScalarT& h,
                           const ScalarT& a, const ScalarT& b, const ScalarT& MC, const ScalarT& HC, const ScalarT& CC,
                           const ScalarT& gPrev, const ScalarT& gNext,
                           const ScalarT& memG, const ScalarT& restV, const ScalarT& Kg, const ScalarT& Ke, const ScalarT& NaG, const ScalarT& NaE,
                           const ScalarT& Ag, const ScalarT& Ae, const ScalarT& CaTg, const ScalarT& CaE, const ScalarT& KCaG)
  {
    ScalarT powN = n * n * n * n;
    ScalarT powM = m * m * m;
    ScalarT powA = a * a * a;
    ScalarT powMC = MC * MC;
    ScalarT powCC = CC * CC * CC * CC;
    ScalarT r = memG * (VSeg - restV) + Kg * powN * (VSeg - Ke ) + NaG * powM * h * (VSeg - NaE )
                + Ag * powA * b * (VSeg - Ae) + CaTg * powMC * HC * (VSeg - CaE) + KCaG * powCC * (VSeg - Ke) - gNext * (VSegN - VSeg) - gPrev * (VSegP - VSeg);
    return r;
  }

  template <typename ScalarT>
  static ScalarT kcl1EquQ( const ScalarT& VSeg, const ScalarT& memC )
  {
    ScalarT r = memC * VSeg;
    return r;
  }

#if 0
  // KCL equation 2 -- -1 * equation 1 because of device symmetry
  template <typename ScalarT>
  static ScalarT kcl2EquF( const ScalarT& Vn1, const ScalarT& Vn2, const ScalarT& n, const ScalarT& m, const ScalarT& h,
                           const ScalarT& a, const ScalarT& b, const ScalarT& MC, const ScalarT& HC, const ScalarT& CC,
                           const ScalarT& memG, const ScalarT& restV, const ScalarT& Kg, const ScalarT& Ke, const ScalarT& NaG, const ScalarT& NaE,
                           const ScalarT& Ag, const ScalarT& Ae, const ScalarT& CaTg, const ScalarT& CaE, const ScalarT& KCaG)
  {
    ScalarT powN = n * n * n * n;
    ScalarT powM = m * m * m;
    ScalarT powA = a * a * a;
    ScalarT powMC = MC * MC;
    ScalarT powCC = CC * CC * CC * CC;
    ScalarT r = -1.0 * (memG * (Vn1 - Vn2 - restV) + Kg * powN * (Vn1 - Vn2 - Ke ) + NaG * powM * h * (Vn1 - Vn2 - NaE )
                        + Ag * powA * b * (Vn1 - Vn2 - Ae) + CaTg * powMC * HC * (Vn1 - Vn2 - CaE) + KCaG * powCC * (Vn1 - Vn2 - Ke) );
    return r;
  }

  template <typename ScalarT>
  static ScalarT kcl2EquQ( const ScalarT& Vn1, const ScalarT& Vn2, const ScalarT& memC )
  {
    ScalarT r = -1.0 * memC * (Vn1 - Vn2);
    return r;
  }
#endif

  // n conservation equation
  template <typename ScalarT>
  static ScalarT nEquF( const ScalarT& Vn1, const ScalarT& n, const ScalarT& Vrest )
  {
    ScalarT vDiff = Vn1; // - Vrest;
    ScalarT alpha = alphaN<ScalarT>(vDiff);
    ScalarT beta = betaN<ScalarT>(vDiff);
    ScalarT r = (alpha + beta) * n - alpha;
    return r;
  }

  template <typename ScalarT>
  static ScalarT nEquQ( const ScalarT& n )
  {
    ScalarT r = n;
    return r;
  }

  // m conservation equation
  template <typename ScalarT>
  static ScalarT mEquF( const ScalarT& Vn1, const ScalarT& m, const ScalarT& Vrest )
  {
    ScalarT vDiff = Vn1; // - Vrest;
    ScalarT alpha = alphaM<ScalarT>(vDiff);
    ScalarT beta = betaM<ScalarT>(vDiff);
    ScalarT r = (alpha + beta) * m - alpha;
    return r;
  }

  template <typename ScalarT>
  static ScalarT mEquQ( const ScalarT& m )
  {
    ScalarT r = m;
    return r;
  }

  // h conservation equation
  template <typename ScalarT>
  static ScalarT hEquF( const ScalarT& Vn1, const ScalarT& h, const ScalarT& Vrest )
  {
    ScalarT vDiff = Vn1; // - Vrest;
    ScalarT alpha = alphaH<ScalarT>(vDiff);
    ScalarT beta = betaH<ScalarT>(vDiff);
    ScalarT r = (alpha + beta) * h - alpha;
    return r;
  }

  template <typename ScalarT>
  static ScalarT hEquQ( const ScalarT& h )
  {
    ScalarT r = h;
    return r;
  }

  // a conservation equation
  template <typename ScalarT>
  static ScalarT aEquF( const ScalarT& Vn1, const ScalarT& a, const ScalarT& Vrest )
  {
    ScalarT vDiff = Vn1; // - Vrest;
    ScalarT Inf = aInf<ScalarT>(vDiff);
    ScalarT Tau = aTau<ScalarT>(vDiff);
    ScalarT r = (a - Inf)/Tau;
    return r;
  }

  template <typename ScalarT>
  static ScalarT aEquQ( const ScalarT& a )
  {
    ScalarT r = a;
    return r;
  }

  // b conservation equation
  template <typename ScalarT>
  static ScalarT bEquF( const ScalarT& Vn1, const ScalarT& b, const ScalarT& Vrest )
  {
    ScalarT vDiff = Vn1; // - Vrest;
    ScalarT Inf = bInf<ScalarT>(vDiff);
    ScalarT Tau = bTau<ScalarT>(vDiff);
    ScalarT r = (b - Inf)/Tau;
    return r;
  }

  template <typename ScalarT>
  static ScalarT bEquQ( const ScalarT& b )
  {
    ScalarT r = b;
    return r;
  }

  // M conservation equation
  template <typename ScalarT>
  static ScalarT M_EquF( const ScalarT& Vn1, const ScalarT& M, const ScalarT& Vrest )
  {
    ScalarT vDiff = Vn1; // - Vrest;
    ScalarT Inf = M_Inf<ScalarT>(vDiff);
    ScalarT Tau = M_Tau<ScalarT>(vDiff);
    ScalarT r = (M - Inf)/Tau;
    return r;
  }

  template <typename ScalarT>
  static ScalarT M_EquQ( const ScalarT& M )
  {
    ScalarT r = M;
    return r;
  }

  // H conservation equation
  template <typename ScalarT>
  static ScalarT H_EquF( const ScalarT& Vn1, const ScalarT& H, const ScalarT& Vrest )
  {
    ScalarT vDiff = Vn1; // - Vrest;
    ScalarT Inf = H_Inf<ScalarT>(vDiff);
    ScalarT Tau = H_Tau<ScalarT>(vDiff);
    ScalarT r = (H - Inf)/Tau;
    return r;
  }

  template <typename ScalarT>
  static ScalarT H_EquQ( const ScalarT& H )
  {
    ScalarT r = H;
    return r;
  }

  // C conservation equation
  template <typename ScalarT>
  static ScalarT C_EquF( const ScalarT& Vn1, const ScalarT& C, const ScalarT& CaConc, const ScalarT& Vrest )
  {
    ScalarT vDiff = Vn1; // - Vrest;
    ScalarT Inf = C_Inf<ScalarT>(vDiff, CaConc);
    ScalarT Tau = C_Tau<ScalarT>(vDiff);
    ScalarT r = (C - Inf)/Tau;
    return r;
  }

  template <typename ScalarT>
  static ScalarT C_EquQ( const ScalarT& C )
  {
    ScalarT r = C;
    return r;
  }

  // Calcium conservation equation
  template <typename ScalarT>
  static ScalarT Ca_EquF( const ScalarT& Vn1, const ScalarT& MC, const ScalarT& HC, const ScalarT& Ca,
                          const ScalarT& CaTg, const ScalarT& CaE, const ScalarT& CaGamma, const ScalarT& CaTau )
  {
    ScalarT r = CaGamma * CaTg * MC * MC * HC * (Vn1 - CaE) + Ca / CaTau;
    return r;
  }

  template <typename ScalarT>
  static ScalarT Ca_EquQ( const ScalarT& Ca)
  {
    ScalarT r = Ca;
    return r;
  }

public:
  // iterator reference to the Neuron model which owns this instance.
  // Getters and setters
  Model &getModel() {
    return model_;
  }

private:

  Model &       model_;         //< Owning model

  // model level parameters that can be overridden at the instance level
  double rInt;     // intracellular resistivity
  double radius;   // Segment radius
  double length;   // cable length (segment length = length/nSeg)
  double segArea;  // segment area (derrived from radius, length and nSeg)
  int    nSeg;     // number of segments
  bool rIntGiven;
  bool radiusGiven;
  bool lengthGiven;
  bool nSegGiven;

  // the cable equation is dependent on parameters from the previous and next segment of other
  // cables attached to the input/ouput nodes (i.e. when we branch a cable).  This lets
  // one set the parameters for the previous/next cable segment.
  double rIntPrevious;
  double radiusPrevious;
  double lengthPrevious;
  double rIntNext;
  double radiusNext;
  double lengthNext;
  bool   rIntPreviousGiven;
  bool   radiusPreviousGiven;
  bool   lengthPreviousGiven;
  bool   rIntNextGiven;
  bool   radiusNextGiven;
  bool   lengthNextGiven;

  // conductance between segments -- calculated from radius, rInt and length and number of segments
  // gBackward connects to previous node. gForward connects to the next node
  vector<double> gBackward;
  vector<double> gForward;

  // derrived quantities computed in updateIntermediateVars
  // and used in the load functions (no q terms on the external nodes)
  double kcl1Fvalue;
  double kcl2Fvalue;
  // internal segments
  vector<double> segFvalue;
  vector<double> segQvalue;
  vector<double> segNEquFvalue, segNEquQvalue;
  vector<double> segMEquFvalue, segMEquQvalue;
  vector<double> segHEquFvalue, segHEquQvalue;
  vector<double> segAEquFvalue, segAEquQvalue;
  vector<double> segBEquFvalue, segBEquQvalue;
  vector<double> segM_EquFvalue, segM_EquQvalue;
  vector<double> segH_EquFvalue, segH_EquQvalue;
  vector<double> segCEquFvalue, segCEquQvalue;
  vector<double> segCaEquFvalue, segCaEquQvalue;

  // jacobian terms
  double dkcl1F_dVin, dkcl1F_dVs0;
  double dkcl2F_dVout, dkcl2F_dVsn;
  // internal equations
  vector<double> segF_dVp, segF_dV, segF_dVn, segF_dn, segF_dm, segF_dh, segF_da, segF_db, segF_dM, segF_dH, segF_dc;
  vector<double> segQ_dV;
  vector<double> dnF_dV, dnF_dn, dnQ_dn;
  vector<double> dmF_dV, dmF_dm, dmQ_dm;
  vector<double> dhF_dV, dhF_dh, dhQ_dh;
  vector<double> daF_dV, daF_da, daQ_da;
  vector<double> dbF_dV, dbF_db, dbQ_db;
  vector<double> dMF_dV, dMF_dM, dMQ_dM;
  vector<double> dHF_dV, dHF_dH, dHQ_dH;
  vector<double> dcF_dV, dcF_dc, dcF_dCa, dcQ_dc;
  vector<double> dCaF_dV, dCaF_dM, dCaF_dH, dCaF_dCa, dCaQ_dCa;

  // state variables
  vector<double> potassiumCurrent;
  vector<double> sodiumCurrent;

  // local state indices (offsets)
  vector<int> li_KCurrentState;
  vector<int> li_NaCurrentState;

  // local solution indices (offsets)
  int li_Pos;      // local index to positive node on this device
  int li_Neg;      // local index to negative node on this device
  // local solution indices for internal vars (variable number of these)
  vector<int> li_Vol;      // local index to segment voltage
  vector<int> li_nPro;     // local index to n promoter value (Na current)
  vector<int> li_mPro;     // local index to m promoter value (K current)
  vector<int> li_hPro;     // local index to h promoter value (K current)
  vector<int> li_aPro;     // local index to a promoter value
  vector<int> li_bPro;     // local index to a promoter value
  vector<int> li_MPro;     // local index to a promoter value
  vector<int> li_HPro;     // local index to a promoter value
  vector<int> li_cPro;     // local index to a promoter value
  vector<int> li_CaPro;     // local index to a promoter value

  // Matrix equation index variables:

  // Offset variables corresponding to the above declared indices.
  int APosEquPosNodeOffset, APosEquNextNodeOffset;
  int ANegEquNegNodeOffset, ANegEquLastNodeOffset;
  vector<int> SegVEqnVpreOffset;
  vector<int> SegVEqnVsegOffset;
  vector<int> SegVEqnVnexOffset;
  vector<int> SegVEqnNOffset;
  vector<int> SegVEqnMOffset;
  vector<int> SegVEqnHOffset;
  vector<int> SegVEqnAOffset;
  vector<int> SegVEqnBOffset;
  vector<int> SegVEqnM_Offset;
  vector<int> SegVEqnH_Offset;
  vector<int> SegVEqnCOffset;
  vector<int> NEquVNodeOffset;
  vector<int> NEquNNodeOffset;
  vector<int> MEquVNodeOffset;
  vector<int> MEquMNodeOffset;
  vector<int> HEquVNodeOffset;
  vector<int> HEquHNodeOffset;
  vector<int> AEquVNodeOffset;
  vector<int> AEquANodeOffset;
  vector<int> BEquVNodeOffset;
  vector<int> BEquBNodeOffset;
  vector<int> M_EquVNodeOffset;
  vector<int> M_EquM_NodeOffset;
  vector<int> H_EquVNodeOffset;
  vector<int> H_EquH_NodeOffset;
  vector<int> CEquVNodeOffset;
  vector<int> CEquCNodeOffset;
  vector<int> CEquCaNodeOffset;
  vector<int> CaEquVNodeOffset;
  vector<int> CaEquM_NodeOffset;
  vector<int> CaEquH_NodeOffset;
  vector<int> CaEquCaNodeOffset;

  vector< vector<int> > jacStamp;
};

//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 06/10/09
//-----------------------------------------------------------------------------
class Model : public DeviceModel
{
  typedef std::vector<Instance *> InstanceVector;

  friend class ParametricData<Model>;
  friend class Instance;

public:
  static ParametricData<Model> &getParametricData();

  virtual const ParametricData<void> &getMyParametricData() const {
    return getParametricData();
  }

  Model(const ModelBlock & MB,
        SolverState & ss1,
        DeviceOptions & do1);
  ~Model();

private:
  Model();
  Model(const Model &);
  Model &operator=(const Model &);

public:
  virtual std::ostream &printOutInstances(std::ostream &os) const;

  bool processParams (string param = "");
  bool processInstanceParams (string param = "");

private:

  // parameter variables
  double cMem;     // membrane capacitance
  double gMem;     // membrane conductance
  double vRest;    // resting potential
  double eNa;      // sodium rest potential
  double gNa;      // sodium base conductance
  double eK;       // potassium rest potential
  double gK;       // potassium base conductance
  double eA;       // a-current rest potential
  double gA;       // a-current base conductance
  double eCa;      // Calcium rest potential
  double gCa;      // Calcium base conductance
  double eKCa;     // potassium-calcium rest potential
  double gKCa;     // potassium-calcium base conductance
  double CaInit;  // initial intra-cellular calcium concentration
  double CaGamma;  // calcium current to concentration multiplier
  double CaTau;    // calcium removal time constant
  double rInt;     // intracellular resistivity
  double radius;   // Segment radius
  double length;   // cable length (segment length = length/nSeg)
  int    nSeg;     // number of segments

  // the cable equation is dependent on parameters from the previous and next segment of other
  // cables attached to the input/ouput nodes (i.e. when we branch a cable).  This lets
  // one set the parameters for the previous/next cable segment.
  double rIntPrevious;
  double radiusPrevious;
  double lengthPrevious;
  double rIntNext;
  double radiusNext;
  double lengthNext;
  bool   rIntPreviousGiven;
  bool   radiusPreviousGiven;
  bool   lengthPreviousGiven;
  bool   rIntNextGiven;
  bool   radiusNextGiven;
  bool   lengthNextGiven;

  // flags that parameters were given
  bool rIntGiven;
  bool radiusGiven;
  bool lengthGiven;
  bool nSegGiven;

  // flags that parameters were given
  bool cMemGiven;
  bool gMemGiven;
  bool vRestGiven;
  bool eNaGiven;
  bool gNaGiven;
  bool eKGiven;
  bool gKGiven;
  bool eAGiven;
  bool gAGiven;
  bool eCaGiven;
  bool gCaGiven;
  bool eKCaGiven;
  bool gKCaGiven;
  bool CaInitGiven;
  bool CaGammaGiven;
  bool CaTauGiven;


public:
  InstanceVector &getInstanceVector() {
    return instanceContainer;
  }

  const InstanceVector &getInstanceVector() const {
    return instanceContainer;
  }

private:
  vector<Instance*> instanceContainer;
};

} // namespace Neuron2
} // namespace Device
} // namespace Xyce

typedef Xyce::Device::Neuron4::Instance N_DEV_NeuronInstance4;
typedef Xyce::Device::Neuron4::Model N_DEV_NeuronModel4;

#endif

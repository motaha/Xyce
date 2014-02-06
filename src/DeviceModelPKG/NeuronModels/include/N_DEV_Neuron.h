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
// Filename       : $RCSfile: N_DEV_Neuron.h,v $
//
// Purpose        : Neuron classes.
//
// Special Notes  :
//
// Creator        : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date  : 01/02/08
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.21.2.1 $
//
// Revision Date  : $Date: 2013/10/03 17:23:33 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Neuron_h
#define Xyce_N_DEV_Neuron_h

#include <Sacado.hpp>

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceTemplate.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>

#ifdef HAVE_MATH_H
#include <math.h>
#endif

namespace Xyce {
namespace Device {
namespace Neuron {

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
// Creation Date : 01/02/08
//-----------------------------------------------------------------------------
class Instance : public DeviceInstance
{
  friend class ParametricData<Instance>;
  friend class Model;
  friend class Master;

public:
  static vector< vector<int> > jacStamp;
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
  // These functions represent the equations that need to be solved
  // for this device.  Since Xyce loads an F and Q contribution, the
  // equations are broken up into their F and Q components.  Thus there
  // is a kcl1EquF() and kcl1EquQ().  Automatic differentiation will
  // be used to generate all the derivatives of these equations for the
  // dF/dX and dQ/dX loads

  // first we list some utility functions for calculating coefficients.
  // alpha and beta equations are taken from Koch.
  // They're generally functions of membrane voltage; here, the membrane voltage is
  // the difference between Vn1 and Vn2.
  // Also, in the Koch formulation, the voltage is relative to the resting potential;
  // equations adjusted here to accommodate a nonzero resting potential vRest.
  // These functions expect V to be in milli-volts and then return values that
  // are in 1/ms.  Thus the extra factor's of 1000 here and there
  //
  // potassium current, functions for activator equation
  template <typename ScalarT>
  static ScalarT alphaN( const ScalarT & Vn1, const ScalarT & Vn2, const ScalarT & Vrest)
  {
    ScalarT vDiff = 1000.0 * (Vn1 - Vn2 - Vrest);  // convert voltage to milli-volts
    // and shift to account for nonzero Vrest
    ScalarT r;
    if ((vDiff > 9.99) && (vDiff < 10.01) )
    {
      r = 1.0/(10.0 * ( std::exp( (10.0 - vDiff)/10.0 )));
    }
    else
    {
      r = (10.0 - vDiff) /
          (100.0 * ( std::exp( (10.0 - vDiff)/10.0 ) - 1.0 ));
    }
    r *= 1000.0; // change from 1/ms to 1/s
    return r;
  }

  template <typename ScalarT>
  static ScalarT betaN( const ScalarT & Vn1, const ScalarT & Vn2, const ScalarT & Vrest)
  {
    ScalarT vDiff = 1000.0 * (Vn1 - Vn2 - Vrest);
    ScalarT r = 0.125 * std::exp( -vDiff/80.0 );
    r *= 1000.0; // change from 1/ms to 1/s
    return r;
  }

  // sodium current, functions for activator equation
  template <typename ScalarT>
  static ScalarT alphaM( const ScalarT & Vn1, const ScalarT  & Vn2, const ScalarT & Vrest)
  {
    ScalarT vDiff = 1000.0 * (Vn1 - Vn2 - Vrest);
    ScalarT r;
    if ((vDiff > 24.99) && (vDiff < 25.01) )
    {
      r = (1.0) /
          (( std::exp( (25.0 - vDiff)/10.0 )));
    }
    else
    {
      r = (25.0 - vDiff) /
          (10.0 * ( std::exp( (25.0 - vDiff)/10.0 ) - 1.0 ));
    }
    r *= 1000.0; // change from 1/ms to 1/s
    return r;
  }

  template <typename ScalarT>
  static ScalarT betaM( const ScalarT & Vn1, const ScalarT & Vn2, const ScalarT & Vrest)
  {
    ScalarT vDiff = 1000.0 * (Vn1 - Vn2 - Vrest);
    ScalarT r = 4.0 * std::exp( -vDiff/18.0 );
    r *= 1000.0; // change from 1/ms to 1/s

    return r;
  }

  template <typename ScalarT>
  static ScalarT alphaH( const ScalarT & Vn1, const ScalarT & Vn2, const ScalarT & Vrest)
  {
    ScalarT vDiff = 1000.0 * (Vn1 - Vn2 - Vrest);
    ScalarT r = 0.07 * std::exp( -vDiff/20.0 );
    r *= 1000.0; // change from 1/ms to 1/s
    return r;
  }

  template <typename ScalarT>
  static ScalarT betaH( const ScalarT & Vn1, const ScalarT & Vn2, const ScalarT & Vrest)
  {
    ScalarT vDiff = 1000.0 * (Vn1 - Vn2 - Vrest);
    ScalarT r = 1.0 / ( std::exp( (30.0 - vDiff)/10.0 ) + 1.0 );
    r *= 1000.0; // change from 1/ms to 1/s
    return r;
  }

  // now the device equations
  // KCL equation 1
  template <typename ScalarT>
  static ScalarT kcl1EquF( const ScalarT& Vn1, const ScalarT& Vn2, const ScalarT& n, const ScalarT& m, const ScalarT& h,
                           const ScalarT& memG, const ScalarT& leakE, const ScalarT& Kg, const ScalarT& Ke, const ScalarT& NaG, const ScalarT& NaE )
  {
    ScalarT powN = n * n * n * n;
    ScalarT powM = m * m * m;
    ScalarT r = memG * (Vn1 - Vn2 - leakE) + Kg * powN * (Vn1 - Vn2 - Ke ) + NaG * powM * h * (Vn1 - Vn2 - NaE );
    return r;
  }

  template <typename ScalarT>
  static ScalarT kcl1EquQ( const ScalarT& Vn1, const ScalarT& Vn2, const ScalarT& memC )
  {
    ScalarT r = memC * (Vn1 - Vn2);
    return r;
  }

  // KCL equation 2 -- -1 * equation 1 because of device symmetry
  template <typename ScalarT>
  static ScalarT kcl2EquF( const ScalarT& Vn1, const ScalarT& Vn2, const ScalarT& n, const ScalarT& m, const ScalarT& h,
                           const ScalarT& memG, const ScalarT& leakE, const ScalarT& Kg, const ScalarT& Ke, const ScalarT& NaG, const ScalarT& NaE )
  {
    ScalarT powN = n * n * n * n;
    ScalarT powM = m * m * m;
    ScalarT r = -1.0*(memG * (Vn1 - Vn2 - leakE) + Kg * powN * (Vn1 - Vn2 - Ke ) + NaG * powM * h * (Vn1 - Vn2 - NaE ));
    return r;
  }

  template <typename ScalarT>
  static ScalarT kcl2EquQ( const ScalarT& Vn1, const ScalarT& Vn2, const ScalarT& memC )
  {
    ScalarT r = -1.0 * memC * (Vn1 - Vn2);
    return r;
  }

  // n conservation equation
  template <typename ScalarT>
  static ScalarT nEquF( const ScalarT& Vn1, const ScalarT& Vn2, const ScalarT& n, const ScalarT& Vrest )
  {
    ScalarT r = alphaN<ScalarT>( Vn1, Vn2, Vrest ) * (1.0 - n ) - betaN<ScalarT>( Vn1, Vn2, Vrest ) * n;
    return r;
  }

  template <typename ScalarT>
  static ScalarT nEquQ( const ScalarT& n )
  {
    ScalarT r = -n;
    return r;
  }

  // m conservation equation
  template <typename ScalarT>
  static ScalarT mEquF( const ScalarT& Vn1, const ScalarT& Vn2, const ScalarT& m, const ScalarT& Vrest )
  {
    ScalarT r = alphaM<ScalarT>( Vn1, Vn2, Vrest ) * (1.0 - m ) - betaM<ScalarT>( Vn1, Vn2, Vrest ) * m;
    return r;
  }

  template <typename ScalarT>
  static ScalarT mEquQ( const ScalarT& m )
  {
    ScalarT r = -m;
    return r;
  }

  // h conservation equation
  template <typename ScalarT>
  static ScalarT hEquF( const ScalarT& Vn1, const ScalarT& Vn2, const ScalarT& h, const ScalarT& Vrest )
  {
    ScalarT r = alphaH<ScalarT>( Vn1, Vn2, Vrest ) * (1.0 - h ) - betaH<ScalarT>( Vn1, Vn2, Vrest ) * h;
    return r;
  }

  template <typename ScalarT>
  static ScalarT hEquQ( const ScalarT& h )
  {
    ScalarT r = -h;
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

  // derrived quantities computed in updateIntermediateVars
  // and used in the load functions
  double kcl1Fvalue, kcl1Qvalue;
  double kcl2Fvalue, kcl2Qvalue;
  double nEquFvalue, nEquQvalue;
  double mEquFvalue, mEquQvalue;
  double hEquFvalue, hEquQvalue;
  double dkcl1F_dV1, dkcl1F_dV2, dkcl1F_dn, dkcl1F_dm, dkcl1F_dh, dkcl1Q_dV1, dkcl1Q_dV2;
  double dkcl2F_dV1, dkcl2F_dV2, dkcl2F_dn, dkcl2F_dm, dkcl2F_dh, dkcl2Q_dV1, dkcl2Q_dV2;
  double dnF_dV1, dnF_dV2, dnF_dn, dnQ_dn;
  double dmF_dV1, dmF_dV2, dmF_dm, dmQ_dm;
  double dhF_dV1, dhF_dV2, dhF_dh, dhQ_dh;

  // state variables
  double potassiumCurrent;
  double sodiumCurrent;

  // local state indices (offsets)
  int li_KCurrentState;
  int li_NaCurrentState;

  // local solution indices (offsets)
  int li_Pos;      // local index to positive node on this device
  int li_Neg;      // local index to negative node on this device
  int li_nPro;     // local index to n promoter value (Na current)
  int li_mPro;     // local index to m promoter value (K current)
  int li_hPro;     // local index to h promoter value (K current)

  // Matrix equation index variables:

  // Offset variables corresponding to the above declared indices.
  int APosEquPosNodeOffset;
  int APosEquNegNodeOffset;
  int APosEquNNodeOffset;
  int APosEquMNodeOffset;
  int APosEquHNodeOffset;

  int ANegEquPosNodeOffset;
  int ANegEquNegNodeOffset;
  int ANegEquNNodeOffset;
  int ANegEquMNodeOffset;
  int ANegEquHNodeOffset;

  int ANEquPosNodeOffset;
  int ANEquNegNodeOffset;
  int ANEquNNodeOffset;

  int AMEquPosNodeOffset;
  int AMEquNegNodeOffset;
  int AMEquMNodeOffset;

  int AHEquPosNodeOffset;
  int AHEquNegNodeOffset;
  int AHEquHNodeOffset;
};

//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/02/08
//-----------------------------------------------------------------------------
class Model : public DeviceModel
{
  typedef std::vector<Instance *> InstanceVector;

  friend class ParametricData<Model>;
  friend class Instance;
  friend class Master;

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
  double gMem;     // membrane conductance of leak current
  double eLeak;    // reversal potential of leak current
  double eNa;      // sodium reversal potential
  double gNa;      // sodium base conductance
  double eK;       // potassium reversal potential
  double gK;       // potassium base conductance
  double vRest;    // resting potential

  // flags that parameters were given
  bool cMemGiven;
  bool gMemGiven;
  bool eLeakGiven;
  bool eNaGiven;
  bool gNaGiven;
  bool eKGiven;
  bool gKGiven;
  bool vRestGiven;


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


//-----------------------------------------------------------------------------
// Class         : Master
// Purpose       :
// Special Notes :
// Creator       : Richard Schiek, SNL, Parallel Computational Sciences
// Creation Date : 06/01/12
//-----------------------------------------------------------------------------
class Master : public Xyce::Device::DeviceTemplate<Model, Instance>
{
  friend class Instance;
  friend class Model;

public:
  Master (
    const string dn,
    const string cn,
    const string dmName,
    LinearDevice linearDev,
    SolverState & ss1,
    DeviceOptions & do1)
    : Xyce::Device::DeviceTemplate<Model, Instance>(
      dn, cn, dmName, linearDev, ss1, do1)
  {

  }

  virtual bool updateState (double * solVec, double * staVec, double * stoVec);
};

} // namespace Neuron
} // namespace Device
} // namespace Xyce

typedef Xyce::Device::Neuron::Instance N_DEV_NeuronInstance;
typedef Xyce::Device::Neuron::Model N_DEV_NeuronModel;
typedef Xyce::Device::Neuron::Master N_DEV_NeuronMaster;

#endif

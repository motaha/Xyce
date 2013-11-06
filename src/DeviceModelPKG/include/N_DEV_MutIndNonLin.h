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
// Filename       : $RCSfile: N_DEV_MutIndNonLin.h,v $
//
// Purpose        : Non-Linear Mutual Inductor classes.
//
// Special Notes  :
//
// Creator        : Richard Schiek, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/21/05
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.68.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:37 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_MutIndNonLin_h
#define Xyce_N_DEV_MutIndNonLin_h

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>

// defines simple container InductorInstanceData
#include <N_DEV_MutIndLin.h>

#include <Teuchos_RefCountPtrDecl.hpp>
#include <fstream>

namespace Xyce {
namespace Device {
namespace MutIndNonLin {

// ---------- Forward Declarations ----------
class Model;

//-----------------------------------------------------------------------------
// Class         : Instance
// Purpose       : This is class refers to a single instance of the nonlinear
//                 mutual inductor device.
// Special Notes :
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 3/21/05
//-----------------------------------------------------------------------------

class Instance : public DeviceInstance
{
  friend class ParametricData<Instance>;
  friend class Model;
  friend class Master;

public:
  static ParametricData<Instance> &getParametricData();

  virtual const ParametricData<void> &getMyParametricData() const {
    return getParametricData();
  }

  Instance(InstanceBlock & IB,
                         Model & Iiter,
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
  map<int,string> & getIntNameMap();
  map<int,string> & getStateNameMap();
  map<int,string> & getStoreNameMap();

  void registerStateLIDs( const vector<int> & staLIDVecRef );
  void registerStoreLIDs( const vector<int> & staLIDVecRef );

  const vector< vector<int> > & jacobianStamp() const;
  void registerJacLIDs( const vector< vector<int> > & jacLIDVec );

  bool processParams (string param = "");
  bool updateTemperature(const double & temp_tmp);
  void updateInductanceMatrix();
  bool updateIntermediateVars ();
  bool updatePrimaryState ();
  bool updateSecondaryState ();
  bool setIC ();

  bool plotfileFlag () {return true;}

  void varTypes( vector<char> & varTypeVec );

  // load functions, residual:
  bool loadDAEQVector ();
  bool loadDAEFVector ();

  // DAE load functions, Jacobian:
  bool loadDAEdQdx ();
  bool loadDAEdFdx ();

  void auxDAECalculations ();

  bool outputPlotFiles ();

  // iterator reference to the inductor model which owns this instance.
  // Getters and setters
  Model &getModel() {
    return model_;
  }

private:

  Model &       model_;         //< Owning model

  // This container bundles up the physical data for each inductor
  // involved in this mutual inductor
  int numInductors;
  double L;
  vector< InductorInstanceData* > instanceData;

  // These vectors let the new param options load and set inductor data
  // the parser passes all of these to us
  vector< string > inductorNames;
  vector< double > inductorInductances;
  vector< string > inductorsNode1;
  vector< string > inductorsNode2;
  // and here's the list of ones we are coupling
  vector< string > couplingInductor;
  vector< double > couplingCoefficient;
  //vector< vector< double > > mutualCouplingCoef;

  // local indices for extra equations
  int li_MagVar;
  int li_RVar;

  // offsets in the jacobian
  int mEquVPosOffset, mEquVNegOffset;
  vector< int > mEquInductorOffsets;
  int mEquMOffset, mEquROffset;

  int rEquROffset;
  vector< int > rEquInductorOffsets;

  // state variable for mag, h and r
  int li_MagVarState;
  int li_MagVarDerivState;
  int li_RVarStore;
  int li_BVarStore;
  int li_HVarStore;

  double nonlinFlag;   // flag created by parser.  Don't need it but must read it in
  bool nonlinFlagGiven;
  double mutualCup;    // mutaul coupling value
  bool mutualCupGiven;

  vector< double > inductanceVals;      // the inductances of the inductors
  vector< vector< double > > LO;        // L' * L (matrix)
  vector< double > inductorCurrents;    // currents through inductors (col vec.)
  vector< double > LOI;                 // LO * I (col vector).

  double temp;         // temperature of this instance
  bool tempGiven;      // flag if temp was given

  // intermediate values needed for nonlinear model
  double qV;
  double delM0;        // modeling constant
  double Happ;
  double He;
  double Heo;
  double delM;
  double Mirrp;
  double Manp;
  double P;
  double dP_dM;
  double dP_dVp;
  double dP_dVn;

  // variables for limiting of non-linear, internal vars
  double Mag_orig;
  double Rvar_orig;

  // these vectors are used repeadily in loadDAEdFdx
  // so rather than create and destroy them on each call
  // we will allocate them in the constructor.
  vector< double > dHe_dI;
  vector< double > dManp_dI;
  vector< double > ddelM_dI;
  vector< double > dMirrp_dI;
  vector< double > dP_dI;

  // scaling crontrol for new equations
  double scalingRHS;   // scaling for loading DAE or RHS components

  vector< vector<int> > jacStamp;

  // output stream for output of internal state if requested by user
  Teuchos::RefCountPtr< ofstream > outputFileStreamPtr;
  bool outputStateVarsFlag;
};

//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
// Special Notes :
// Creator       : Rich Schiek, SNL, Parallel Computational Sciences
// Creation Date : 3/21/05
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

  // Data Members for Associations

public:
  InstanceVector &getInstanceVector() {
    return instanceContainer;
  }

  const InstanceVector &getInstanceVector() const {
    return instanceContainer;
  }

private:
  vector<Instance*> instanceContainer;

private:

  double A;                  // Thermal energy parameter (amp/m)
  double Alpha;              // domain coupling parameter (dimensionless)
  double Area;               // mean magnetic cross-sectional area (m^2)
  double BetaH;              // modeling constant (dimensionless)
  double BetaM;              // modeling constant (dimensionless)
  double C;                  // domain flesing parameter (dimensionless)
  double DeltaV;             // smoothing coefficient for V_1 in tanh
  double Gap;                // effective air gap (m)
  double Kirr;               // domain anisotropy parameter (amp/m)
  double Ms;                 // saturation magnetization (amp/m)
  double LevelIgnored;       // for pspice compatibility -- ignored
  double PackIgnored;        // for pspice compatibility -- ignored
  double Path;               // total mean magnetic path (m)
  double Vinf;               // smoothing coefficient for V+1 in tanh
  double tempCoeff1;         // first order temperature coeff.
  double tempCoeff2;         // second order temperature coeff.
  double tnom;               // reference temperature
  double pZeroTol;           // absolute value below which to consider P=0
  double HCgsFactor;         // conversion factor to change H from SI units to CGS units (H/M to Oersted)
  double BCgsFactor;         // conversion factor to change B form SI units (Tesla) to CGS units (Gauss)
  double mVarScaling;        // scaling for M variable
  double rVarScaling;        // scaling for R variable
  double mEqScaling;         // scaling for M equation
  double rEqScaling;         // scaling for r equation
  int outputStateVars;       // flag indicating if user wants M,H and R output
  int factorMS;              // flag to factor Ms out of M
  int BHSiUnits;             // flag to indicate that B and H should be output in SI units. Default is CGS
                             // units for output while SI units are used for calculations.

  // flags indicating if temperature parameters were given
  bool tc1Given;
  bool tc2Given;
  bool tnomGiven;
};

//-----------------------------------------------------------------------------
// Class         : Master
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
class Master : public Xyce::Device::DeviceTemplate<Model, Instance>
{
  public:
    Master (
      const std::string &dn,
      const std::string &cn,
      const std::string &dmName,
           LinearDevice linearDev,
           SolverState & ss1,
           DeviceOptions & do1)
      : Xyce::Device::DeviceTemplate<Model, Instance>(
           dn, cn, dmName, linearDev, ss1, do1)
    {

    }

    virtual bool updateState (double * solVec, double * staVec, double * stoVec);
    virtual bool updateSecondaryState (double * staDeriv, double * stoVec);

    // load functions, residual:
    virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec, double * storeLeadF, double * storeLeadQ);
    //virtual bool loadDAEQVector (double * solVec, double * qVec);
    //virtual bool loadDAEFVector (double * solVec, double * fVec);

    // load functions, Jacobian:
    virtual bool loadDAEMatrices (N_LAS_Matrix & dFdx, N_LAS_Matrix & dQdx);
    //virtual bool loadDAEdQdx (N_LAS_Matrix & dQdx);
    //virtual bool loadDAEdFdx (N_LAS_Matrix & dFdx);

    friend class Instance;
    friend class Model;
};

} // namespace MutIndNonLin
} // namespace Device
} // namespace Xyce

typedef Xyce::Device::MutIndNonLin::Instance N_DEV_MutIndNonLinInstance;
typedef Xyce::Device::MutIndNonLin::Model N_DEV_MutIndNonLinModel;
typedef Xyce::Device::MutIndNonLin::Master N_DEV_MutIndNonLinMaster;

#endif

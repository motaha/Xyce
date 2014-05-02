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

//-----------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_MutIndNonLin2.h,v $
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
// Revision Number: $Revision: 1.26.2.1 $
//
// Revision Date  : $Date: 2014/02/26 20:16:30 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_MutIndNonLin2_h
#define Xyce_N_DEV_MutIndNonLin2_h

// ----------   Xyce Includes   ----------
#include <N_DEV_Configuration.h>
#include <N_DEV_Const.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>
#include <N_DEV_MutIndLin.h>

#include <N_DEV_MutIndNonLin.h>

#include <Teuchos_RefCountPtrDecl.hpp>
#include <Sacado.hpp>
#include <fstream>

namespace Xyce {
namespace Device {
namespace MutIndNonLin2 {

class Model;
class Instance;

struct Traits : public DeviceTraits<Model, Instance, MutIndNonLin::Traits>
{
  static const char *name() {return "Nonlinear Mutual Inductor";}
  static const char *deviceTypeName() {return "K level 2";}
  static const int numNodes() {return 2;}
  static const bool modelRequired() {return true;}
  static const bool isLinearDevice() {return false;}

  static Device *factory(const Configuration &configuration, const FactoryBlock &factory_block);
  static void loadModelParameters(ParametricData<Model> &model_parameters);
  static void loadInstanceParameters(ParametricData<Instance> &instance_parameters);
};

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
  friend class Traits;
    
public:

  Instance(
     const Configuration &       configuration,
     const InstanceBlock &     IB,
     Model &                   Iiter,
     const FactoryBlock &      factory_block);

  ~Instance();

private:
  Instance(const Instance &);
  Instance &operator=(const Instance &);

public:
  void registerLIDs( const std::vector<int> & intLIDVecRef,
                     const std::vector<int> & extLIDVecRef );
  std::map<int,std::string> & getIntNameMap();
  void registerStateLIDs( const std::vector<int> & staLIDVecRef );

  const std::vector< std::vector<int> > & jacobianStamp() const;
  void registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec );

  bool processParams ();
  bool updateTemperature(const double & temp_tmp);
  void updateInductanceMatrix();
  bool updateIntermediateVars ();
  bool updatePrimaryState ();
  bool updateSecondaryState ();
  bool loadDeviceMask();
  bool setIC ();

  bool plotfileFlag () {return true;}

  void varTypes( std::vector<char> & varTypeVec );

  void acceptStep();

  // load functions, residual:
  bool loadDAEQVector ();
  bool loadDAEFVector ();

  // load functions, Jacobian:
  bool loadDAEdQdx ();
  bool loadDAEdFdx ();

  void auxDAECalculations ();

  bool outputPlotFiles ();

  // iterator reference to the inductor model which owns this instance.
  // Getters and setters
  Model &getModel() 
  {
    return model_;
  }

private:

  Model &       model_;         //< Owning model

  // This container bundles up the physical data for each inductor
  // involved in this mutual inductor
  int numInductors;
  double L;
  std::vector< InductorInstanceData* > instanceData;

  // These vectors let the new param options load and set inductor data
  // the parser passes all of these to us
  std::vector< std::string > inductorNames;
  std::vector< double > inductorInductances;
  std::vector< std::string > inductorsNode1;
  std::vector< std::string > inductorsNode2;
  // and here's the list of ones we are coupling
  std::vector< std::string > couplingInductor;
  std::vector< double > couplingCoefficient;
  //std::vector< std::vector< double > > mutualCouplingCoef;

  // local indices for extra equations
  //int li_deltaHappVar;
  int li_deltaMagVar;

  // offsets in the jacobian
  //int deltaHappEquDeltaHappOffset;
  //std::vector< int > deltaHappEquInductorOffsets;

  int deltaMEquDeltaMOffset;
  //int deltaMEquDeltaHappOffset;
  std::vector< int > deltaMEquInductorOffsets;

  double nonlinFlag;   // flag created by parser.  Don't need it but must read it in
  bool nonlinFlagGiven;
  double mutualCup;    // mutaul coupling value
  bool mutualCupGiven;

  std::vector< double > inductanceVals;      // the inductances of the inductors
  std::vector< std::vector< double > > LO;        // L' * L (matrix)
  std::vector< double > inductorCurrents;    // currents through inductors (col vec.)
  std::vector< double > LOI;                 // LO * I (col vector).

  double temp;         // temperature of this instance
  bool tempGiven;      // flag if temp was given

  double branchCurrentSum;
  double deltaBranchCurrentSum;
  double P;
  double PPreviousStep;
  double dP_dM;
  double dP_dBranchCurrentSum;
  double dP_dV1Pos;
  double dP_dV1Neg;
  double mEquFval;
  std::vector<double> dHe_dI;
  std::vector<double> dManp_dI;
  std::vector<double> ddelM_dI;
  std::vector<double> dMirrp_dI;
  std::vector<double> dP_dI;
  double MagVar;
  double oldBranchCurrentSum;  // last branch current sum.
  double MagVarUpdate;
  double lastMagUpdate;
  bool useRKIntegration;
  std::vector<double> branchCurrentSumHistory;  // these two vectors are used for 4th order RK
  std::vector<double> PFunctionHistory;         // integration of dM/dH.

  bool includeDeltaM;  // flag to include deltaM in solution

  // non static jacStamp as each mutual inductor will have a variable number of components
  std::vector< std::vector<int> > jacStamp;

  // output stream for output of internal state if requested by user
  Teuchos::RefCountPtr< std::ofstream > outputFileStreamPtr;
  bool outputStateVarsFlag;

  // this is a templated function for a complicated term P(M,I_1... I_n) that relates
  // the magnetic saturation of the mutual indcutor to the individual currents
  // through the inductors.  We'll need dP_dM and this tempated function automates
  // that calculation via Sacado
  template <typename ScalarT>
  ScalarT Pcalc( const ScalarT & Mag, const ScalarT & CurrentSum, const ScalarT & Vpos, const ScalarT & Vneg);
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
  friend class Traits;
    
public:
  Model(
     const Configuration &       configuration,
     const ModelBlock &        MB,
     const FactoryBlock &      factory_block);
  ~Model();

private:
  Model();
  Model(const Model &);
  Model &operator=(const Model &);

public:
  virtual void forEachInstance(DeviceInstanceOp &op) const /* override */;
    
  virtual std::ostream &printOutInstances(std::ostream &os) const;

  bool processParams ();
  bool processInstanceParams ();

  // Data Members for Associations

public:
  void addInstance(Instance *instance) 
  {
    instanceContainer.push_back(instance);
  }

  InstanceVector &getInstanceVector() 
  {
    return instanceContainer;
  }

  const InstanceVector &getInstanceVector() const 
  {
    return instanceContainer;
  }

private:
  std::vector<Instance*> instanceContainer;

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
  double mVarScaling;        // scaling for M variable
  double rVarScaling;        // scaling for R variable
  double mEqScaling;         // scaling for M equation
  double rEqScaling;         // scaling for r equation
  int outputStateVars;       // flag indicating if user wants M,H and R output

  int factorMS;              // flag to factor Ms out of M (not used in level 2)
  int includeDeltaM;         // flag to make delta M calculation implicit
  int useRKIntegration;      // flag to use 4th order runga-kutta for dM/dH integration
  int useStateDeriv;         // flag to use state vector for dH/dt calculation
  int voltageLimiterFlag;    // flag indicating that we should use limiting on internal vars Mag and R
  double magLimitThres;      // iteration threshold overwhich changes in Mag var are limited
  double rLimitThres;        // iteration threshold over which changes in R var are limited

  // flags indicating if temperature parameters were given
  bool tc1Given;
  bool tc2Given;
  bool tnomGiven;
  bool includeDeltaMGiven;
  bool useRKIntegrationGiven;
};

void registerDevice();

} // namespace MutIndNonLin2
} // namespace Device
} // namespace Xyce

typedef Xyce::Device::MutIndNonLin2::Instance N_DEV_MutIndNonLin2Instance;
typedef Xyce::Device::MutIndNonLin2::Model N_DEV_MutIndNonLin2Model;

#endif

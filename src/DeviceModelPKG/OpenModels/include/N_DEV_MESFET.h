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
// Filename       : $RCSfile: N_DEV_MESFET.h,v $
//
// Purpose        : MESFET transistor classes.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 02/28/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.51.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:37 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_MESFET_h
#define Xyce_N_DEV_MESFET_h

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceTemplate.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>

namespace Xyce {
namespace Device {
namespace MESFET {

// ---------- Forward Declarations -------
class Model;

//--------------- Class Instance ----------------------------------
//-----------------------------------------------------------------------------
// Class         : Instance
// Purpose       :
// Special Notes :
// Creator       : pmc
// Creation Date : 11/16/2003
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
                        Model & Miter,
                        MatrixLoadData & mlData1,
                        SolverState &ss1,
                        ExternData  &ed1,
                        DeviceOptions &do1);
    ~Instance();

private:
  Instance(const Instance &);
  Instance &operator=(const Instance &);

public:
    void registerLIDs( const vector<int> & intLIDVecRef,
                        const vector<int> & extLIDVecRef );
    void registerStateLIDs( const vector<int> & staLIDVecRef );
    void registerStoreLIDs(const vector<int> & stoLIDVecRef);

    map<int,string> & getIntNameMap ();
    map<int,string> & getStoreNameMap ();

    const vector< vector<int> > & jacobianStamp() const;
    void registerJacLIDs( const vector< vector<int> > & jacLIDVec );

    // load functions, residual:
    bool loadDAEQVector ();
    bool loadDAEFVector ();

    // load functions, Jacobian:
    bool loadDAEdQdx ();
    bool loadDAEdFdx ();

    void setupPointers ();

    bool updateIntermediateVars ();
    bool updatePrimaryState ();

    bool processParams (string param = "");
    bool updateTemperature(const double & temp_tmp);

    inline bool isConverged();

public:
  // Getters and setters
  Model &getModel() {
    return model_;
  }

private:
    static vector< vector<int> > jacStamp_DC_SC;
    static vector< vector<int> > jacStamp_DC;
    static vector< vector<int> > jacStamp_SC;
    static vector< vector<int> > jacStamp;

    static vector<int> jacMap_DC_SC;
    static vector<int> jacMap_DC;
    static vector<int> jacMap_SC;
    static vector<int> jacMap;

    static vector< vector<int> > jacMap2_DC_SC;
    static vector< vector<int> > jacMap2_DC;
    static vector< vector<int> > jacMap2_SC;
    static vector< vector<int> > jacMap2;


  Model &       model_;         //< Owning model

    bool limitedFlag;      // device convergence.
    bool   off;            // non-zero indicates device is off for dc analysis
    int    ic;             // initial VDS, VGS vector
    double area;           // area factor
    double ic_vds;         // initial D-S voltage
    double ic_vgs;         // initial G-S voltage
    double temp;           // instance temperature
    double drainCond;      // drain conductance
    double sourceCond;     // source conductance
    double tCGS;           // temp. corr. gate-source capacitance
    double tCGD;           // temp. corr. gate-drain capacitance
    double tIS;            // temp. corr. saturation current
    double tPB;            // temp. corr. gate potential
    double tMESb;          // temp. corr. dopinng tail parameter
    double tBeta;          // temp. corr. Beta
    double tvt0;           // temp. corr. vt0
    double tLambda;        // temp. corr. Lambda
    double tAlpha;         // saturation voltage parameter
    double tRD;            // temp. corr. drain resistance
    double tRS;            // temp. corr. source resistance
    double vt;             // set in updateTemperature to CONSTKoverQ*temp

    int dNode;             // number of drain node
    int gNode;             // number of gate node
    int sNode;             // number of source node
    int dpNode;            // internal drain node
    int spNode;            // internal source node
    double Vgs;            // voltage G-S
    double Vgd;            // voltage G-D
    double gm;             // transconductance
    double gds;            // dc conductance D-S
    double ggs;            // dc conductance G-S
    double ggd;            // dc conductance G-D
    double cdrain;         // channel current
    double cd;             // drain current
    double cg;             // total dc gate current
    double cgd;            // dc gate-drain current
    double corDepCap;      // joining point of fwd bias dep. cap.
    double vcrit;          // critical voltage
    double f1;             // coeff. of capacitance polynomial exponent
    double f2;
    double f3;
    double Bfac;           // doping profile parameter
    double p;              // power dissipated by the MESFET

    double GateSourceOverlapCap ;
    double GateDrainOverlapCap ;

    // Solution variables and intermediate quantities
    // drain,source,gate, drainprime and sourceprime voltages
    double Vd;
    double Vs;
    double Vg;
    double Vdp;
    double Vsp;

    // voltage drops between pairs of nodes
    double Vddp; // drain-drain'
    double Vssp; // source-source'
    double Vgsp; // gate-source'
    double Vgdp; // gate-drain'
    double Vdpsp; //drop across channel

    // vector local indices
    int li_Drain;
    int li_DrainPrime;
    int li_Source;
    int li_SourcePrime;
    int li_Gate;

    /////////////////////////////////////////////////////////////////
    //  Jacobian matrix indices:
    //  This is a 5x5 matrix block, of which 15 entries are nonzero:
    //
    // ----------------------------------------------------
    // | #NZ     |       |                                |
    // | entries |       |  V_d   V_g   V_s   V_d'  V_s'  |
    // ----------------------------------------------------
    // |    2    | KCL_d |   a                  b         |
    // |    3    | KCL_g |         c            e     f   |
    // |    2    | KCL_s |               g            h   |
    // |    4    | KCL_d'|   m     n            p     q   |
    // |    4    | KCL_s'|         r     s      u     v   |
    // ----------------------------------------------------
    //     15 total

    ////////////////////////////////////////////////////////////////////
    // Offset variables corresponding to the above declared indices.

    // Jacobian Matrix Offset:

 // V_d Row:
    int ADrainEquDrainNodeOffset;             // a
    int ADrainEquDrainPrimeNodeOffset;        // b

  // V_g Row:
    int AGateEquGateNodeOffset;               // c
    int AGateEquDrainPrimeNodeOffset;         // d
    int AGateEquSourcePrimeNodeOffset;        // e

  // V_s Row:
    int ASourceEquSourceNodeOffset;           // f
    int ASourceEquSourcePrimeNodeOffset;      // g

  // V_d' Row:
    int ADrainPrimeEquDrainNodeOffset;        // h
    int ADrainPrimeEquGateNodeOffset;         // i
    int ADrainPrimeEquDrainPrimeNodeOffset;   // j
    int ADrainPrimeEquSourcePrimeNodeOffset;  // k

  // V_s' Row:
    int ASourcePrimeEquGateNodeOffset;        // l
    int ASourcePrimeEquSourceNodeOffset;      // m
    int ASourcePrimeEquDrainPrimeNodeOffset;  // n
    int ASourcePrimeEquSourcePrimeNodeOffset; // o


#ifndef Xyce_NONPOINTER_MATRIX_LOAD
    ////////////////////////////////////////////////////////////////////
    // dFdx Matrix Ptr:
 // V_d Row:
    double * f_DrainEquDrainNodePtr;             // a
    double * f_DrainEquDrainPrimeNodePtr;        // b

  // V_g Row:
    double * f_GateEquGateNodePtr;               // c
    double * f_GateEquDrainPrimeNodePtr;         // d
    double * f_GateEquSourcePrimeNodePtr;        // e

  // V_s Row:
    double * f_SourceEquSourceNodePtr;           // f
    double * f_SourceEquSourcePrimeNodePtr;      // g

  // V_d' Row:
    double * f_DrainPrimeEquDrainNodePtr;        // h
    double * f_DrainPrimeEquGateNodePtr;         // i
    double * f_DrainPrimeEquDrainPrimeNodePtr;   // j
    double * f_DrainPrimeEquSourcePrimeNodePtr;  // k

  // V_s' Row:
    double * f_SourcePrimeEquGateNodePtr;        // l
    double * f_SourcePrimeEquSourceNodePtr;      // m
    double * f_SourcePrimeEquDrainPrimeNodePtr;  // n
    double * f_SourcePrimeEquSourcePrimeNodePtr; // o

    ////////////////////////////////////////////////////////////////////
    // dQdx Matrix Ptr:

 // V_d Row:
    double * q_DrainEquDrainNodePtr;             // a
    double * q_DrainEquDrainPrimeNodePtr;        // b

  // V_g Row:
    double * q_GateEquGateNodePtr;               // c
    double * q_GateEquDrainPrimeNodePtr;         // d
    double * q_GateEquSourcePrimeNodePtr;        // e

  // V_s Row:
    double * q_SourceEquSourceNodePtr;           // f
    double * q_SourceEquSourcePrimeNodePtr;      // g

  // V_d' Row:
    double * q_DrainPrimeEquDrainNodePtr;        // h
    double * q_DrainPrimeEquGateNodePtr;         // i
    double * q_DrainPrimeEquDrainPrimeNodePtr;   // j
    double * q_DrainPrimeEquSourcePrimeNodePtr;  // k

  // V_s' Row:
    double * q_SourcePrimeEquGateNodePtr;        // l
    double * q_SourcePrimeEquSourceNodePtr;      // m
    double * q_SourcePrimeEquDrainPrimeNodePtr;  // n
    double * q_SourcePrimeEquSourcePrimeNodePtr; // o
#endif

    ////////////////////////////////////////////////////////////////////
    // 3f5 State Variables & related quantities:
    // voltage drops
    double vgs;
    double vgd;
    double vds;

    // "original" versions of various voltage drop variables:
    double vgs_orig;
    double vgd_orig;
    double vds_orig;

    // "old" versions of various voltage drop variables:
    double vgs_old;
    double vgd_old;
    double vds_old;

    int newtonIterOld;

    int mode;       // device mode : 1 = normal, -1 = inverse

    //gate-source capacitor
    double capgs;   // value
    double qgs;     // charge
    double cqgs;    // trans. G-S current
    double gcgs;    // conductance, d(cqgs)/dVgs

    // gate-drain capacitor
    double capgd;   // value
    double qgd;     // charge
    double cqgd;    // trans. G-D current
    double gcgd;    // conductance, d(cqgd)/dVgd

    // current through source and drain resistors
    double Isource;
    double Idrain;

     // local indices
    int li_store_vgs;
    int li_store_vgd;
    int li_store_dev_id;
    int li_store_dev_ig;
    int li_store_dev_is;

    int li_state_qgs;
    int li_state_gcgs;
    int li_state_qgd;
    int li_state_gcgd;
};

//------------------- Class Model ---------------------------------
//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
// Special Notes :
// Creator       : pmc
// Creation Date : 11/16/2003
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

    //user-specified parameters
      double AF;      // flicker noise exponent, default=1
      double B;       // doping tail parameter, default=1
      double ALPHA;   // saturation voltage paramter, default=2
      double BETA;    // transconductance paramter (A/V^2), default=1.0E-4
      double CGS;     // zero-bias G-S junction capacitance (F), default=5pf
      double CGD;     // zero-bias G-D junction capacitance (F), default=1pf
      double FC;      // coeff. for forward-bias depletion cap., default=0.5
      double IS;      // gate junction saturation current (A), default=1.0E-14
      double KF;      // flicker noise coefficient, default=0
      double LAMBDA;  // channel-length modulation parameter (1/V), default=0
      double PB;      // gate junction potential (V), default=1
      double RD;      // drain  ohmic resistance (ohms), default=0
      double RS;      // source  ohmic resistance (ohms), default=0
      double TNOM;    // parameter measurement temperature (C), default=27
      double VTO;     // threshhold voltage (V), default=-2.0
      double fNcoef;  // default=0;
      double fNexp;   // default=1;


      int dtype;      //  device type:  1 = NJF, -1 = PJF
};

//-----------------------------------------------------------------------------
// Function      : Instance:isConverged ()
// Purpose       : Return whether a MOSFET device has done something that
//                  should be interpreted as invalidating other convergence
//                  tests
//                  In case of mesfet, just do it if the limiter function
//                  pnjlim changed anything.
//                  This actually agrees with how the Check flag
//                  is used in Spice3F5 mesfetload.c
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 03/22/05
//-----------------------------------------------------------------------------
inline bool Instance::isConverged()
{
  return (!limitedFlag);
}

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

    virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec, double * storeLeadF, double * storeLeadQ);
    virtual bool loadDAEMatrices (N_LAS_Matrix & dFdx, N_LAS_Matrix & dQdx);

    friend class Instance;
    friend class Model;
};

} // namespace MESFET
} // namespace Device
} // namespace Xyce

typedef Xyce::Device::MESFET::Instance N_DEV_MESFETInstance;
typedef Xyce::Device::MESFET::Model N_DEV_MESFETModel;
typedef Xyce::Device::MESFET::Master N_DEV_MESFETMaster;

#endif

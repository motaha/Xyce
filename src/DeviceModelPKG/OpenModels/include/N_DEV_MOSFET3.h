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
// Filename       : $RCSfile: N_DEV_MOSFET3.h,v $
//
// Purpose        : Level 3 Metal-oxide-semiconductor field effect transistor
//                  (MOSFET) classes.
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
// Revision Number: $Revision: 1.97.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:37 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_MOSFET3_h
#define Xyce_N_DEV_MOSFET3_h

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceTemplate.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>
#include <N_DEV_DeviceBlock.h>

namespace Xyce {
namespace Device {
namespace MOSFET3 {

// ---------- Forward Declarations -------
class Model;

//-----------------------------------------------------------------------------
// Class         : Instance
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
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
             DeviceOptions & do1);

    ~Instance();

private:
  Instance(const Instance &);
  Instance &operator=(const Instance &);

public:
    void registerLIDs( const vector<int> & intLIDVecRef,
                       const vector<int> & extLIDVecRef );
    void registerStateLIDs( const vector<int> & staLIDVecRef);
    void registerStoreLIDs(const vector<int> & stoLIDVecRef);

    map<int,string> & getIntNameMap ();
    map<int,string> & getStoreNameMap ();

    const vector< vector<int> > & jacobianStamp() const;
    void registerJacLIDs( const vector< vector<int> > & jacLIDVec );

    bool processParams (string param = "");

    bool updateTemperature(const double & temp_tmp);
    bool updateIntermediateVars ();
    bool updatePrimaryState ();

    // load functions, residual:
    bool loadDAEQVector ();
    bool loadDAEFVector ();

    // load functions, Jacobian:
    bool loadDAEdQdx ();
    bool loadDAEdFdx ();

    void setupPointers();

    // Additional Public Declarations
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

  protected:
  private:
    // instance variables ripped -- bleeding and without anesthetic -- from
    // 3f5, with obvious modifications to names (remove MOS3 prefix)
    int states;     // index into state table for this device
    int dNode;  // number of the gate node of the mosfet
    int gNode;  // number of the gate node of the mosfet
    int sNode;  // number of the source node of the mosfet
    int bNode;  // number of the bulk node of the mosfet
    int dNodePrime; // number of the internal drain node of the mosfet
    int sNodePrime; // number of the internal source node of the mosfet


    bool OFF;                   // device initialized OFF (vbs=vgs=vds=0)

    double l;   // the length of the channel region
    double w;   // the width of the channel region
    double drainArea;   // the area of the drain diffusion
    double sourceArea;  // the area of the source diffusion
    double drainSquares;    // the length of the drain in squares
    double sourceSquares;   // the length of the source in squares
    double drainPerimeter;
    double sourcePerimeter;
    double sourceConductance;   //conductance of source(or 0):set in setup
    double drainConductance;    //conductance of drain(or 0):set in setup
    double temp;    // operating temperature of this instance
    double numberParallel; // number simulated parallel mosfets

    double tTransconductance;   // temperature corrected transconductance
    double tSurfMob;            // temperature corrected surface mobility
    double tPhi;                // temperature corrected Phi
    double tVto;                // temperature corrected Vto
    double tSatCur;             // temperature corrected saturation Cur.
    double tSatCurDens; // temperature corrected saturation Cur. density
    double tCbd;                // temperature corrected B-D Capacitance
    double tCbs;                // temperature corrected B-S Capacitance
    double tCj;         // temperature corrected Bulk bottom Capacitance
    double tCjsw;       // temperature corrected Bulk side Capacitance
    double tBulkPot;    // temperature corrected Bulk potential
    double tDepCap;     // temperature adjusted transition point in
                            // the cureve matching Fc * Vj
    double tVbi;        // temperature adjusted Vbi

    double icVBS;   // initial condition B-S voltage
    double icVDS;   // initial condition D-S voltage
    double icVGS;   // initial condition G-S voltage
    double von;
    double vdsat;
    double sourceVcrit; // vcrit for pos. vds
    double drainVcrit;  // vcrit for neg. vds
    double cd;
    double cbs;
    double cbd;
    double gmbs;
    double gm;
    double gds;
    double gbd;
    double gbs;
    double capbd;
    double capbs;
    double Cbd;
    double Cbdsw;
    double Cbs;
    double Cbssw;
    double f2d;
    double f3d;
    double f4d;
    double f2s;
    double f3s;
    double f4s;
    int mode;       // device mode : 1 = normal, -1 = inverse
    double mode_low;
    double mode_high;

    bool limitedFlag;  // for convergence testing.

    bool IC_GIVEN;

    //end of 3f5 outtakes

    ////////////////////////////////////////////////////////////////////
    // these are intermediate variables added to the instance class instead
    // of leaving them to be calculated repeatedly in the load function

    // some caluclated quantities
    double EffectiveLength;
    double DrainSatCur;
    double SourceSatCur;
    double GateSourceOverlapCap;
    double GateDrainOverlapCap;
    double GateBulkOverlapCap;
    double OxideCap;

    // Solution variables and intermediate quantities
    // drain,source,gate, bulk, drainprime and sourceprime voltages
    double Vd;
    double Vs;
    double Vg;
    double Vb;
    double Vdp;
    double Vsp;
    // voltage drops between pairs of nodes
    double Vddp; // drain-drain'
    double Vssp; // source-source'
    double Vbsp; // bulk-source'
    double Vbdp; // bulk-drain'
    double Vgsp; // gate-source'
    double Vgdp; // gate-drain'
    double Vgb;  //gate-bulk
    double Vdpsp; //drop across channel

    // the gate-drain voltage drop isn't actually a state variable, but it
    // is calculated at the same time and in the same manner as the state
    // vars.  So here we go, sticking it in the instance class.
    double vgd;

    // Some stuff from mos3temp that were local vars but used elsewhere
    double vt;   // set in updateTemperature to CONSTKoverQ*temp


    // the variables capgs, capgd and capgb are the raw output of
    // qmeyer.  They get massaged into total capacitances in
    // updateIntermediateVars, and get used in updatePrimaryState to get
    // charges on the capacitors.

    double Capgs;   // total gate-source capacitance
    double Capgd;   // total gate-drain capacitance
    double Capgb;   // total gate-bulk capacitance

    // current through source and drain resistors
    double Isource;
    double Idrain;
    
    double cdrain;  // the channel current shouldn't be a local variable in */
		    // updateIntermediateVars!

    // these are calculated in loadRHS and used in the jacobian load
    double Gm,Gmbs;  // we do this so we don't really need the xnrm/xrev vars
    double revsum;   // described in comments at the end of
    double nrmsum;   // updateIntermediateVars (uIVB in remaining comments)
    double cdreq;

    // end of intermediate variables that aren't state variables
    ////////////////////////////
    //
    // vector indices
    int li_Drain;
    int li_DrainPrime;
    int li_Source;
    int li_SourcePrime;
    int li_Gate;
    int li_Bulk;

    ////////////////////////////////////////////////////////////////////
    // The following verbatim from Level=1, which has the same jacobian
    // structure
    ////////////////////////////////////////////////////////////////////
    //  Jacobian matrix indices:
    //  This is a 6x6 matrix block, of which 22 entries are nonzero:
    //
    // ---------------------------------------------------------
    // | #NZ     |       |                                     |
    // | entries |       |  V_d   V_g   V_s   V_b   V_d'  V_s' |
    // ---------------------------------------------------------
    // |    2    | KCL_d |   a                       b         |
    // |    4    | KCL_g |         c           d     e     f   |
    // |    2    | KCL_s |               g                 h   |
    // |    4    | KCL_b |         i           j     k     l   |
    // |    5    | KCL_d'|   m     n           o     p     q   |
    // |    5    | KCL_s'|         r     s     t     u     v   |
    // ---------------------------------------------------------
    //     22 total

    ////////////////////////////////////////////////////////////////////
    // Offset variables corresponding to the above declared indices.

    // Jacobian Matrix Offset:

  // V_d Row:
    int ADrainEquDrainNodeOffset;             // a
    int ADrainEquDrainPrimeNodeOffset;        // b

  // V_g Row:
    int AGateEquGateNodeOffset;               // c
    int AGateEquBulkNodeOffset;               // d
    int AGateEquDrainPrimeNodeOffset;         // e
    int AGateEquSourcePrimeNodeOffset;        // f

  // V_s Row:
    int ASourceEquSourceNodeOffset;           // g
    int ASourceEquSourcePrimeNodeOffset;      // h

  // V_b Row:
    int ABulkEquGateNodeOffset;               // i
    int ABulkEquBulkNodeOffset;               // j
    int ABulkEquDrainPrimeNodeOffset;         // k
    int ABulkEquSourcePrimeNodeOffset;        // l

  // V_d' Row:
    int ADrainPrimeEquDrainNodeOffset;        // m
    int ADrainPrimeEquGateNodeOffset;         // n
    int ADrainPrimeEquBulkNodeOffset;         // o
    int ADrainPrimeEquDrainPrimeNodeOffset;   // p
    int ADrainPrimeEquSourcePrimeNodeOffset;  // q

  // V_s' Row:
    int ASourcePrimeEquGateNodeOffset;        // r
    int ASourcePrimeEquSourceNodeOffset;      // s
    int ASourcePrimeEquBulkNodeOffset;        // t
    int ASourcePrimeEquDrainPrimeNodeOffset;  // u
    int ASourcePrimeEquSourcePrimeNodeOffset; // v

  // Jacobian Matrix Pointers:

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  // F-vector pointers:
  // V_d Row:
    double * f_DrainEquDrainNodePtr;             // a
    double * f_DrainEquDrainPrimeNodePtr;        // b

  // V_g Row:
    double * f_GateEquGateNodePtr;               // c
    double * f_GateEquBulkNodePtr;               // d
    double * f_GateEquDrainPrimeNodePtr;         // e
    double * f_GateEquSourcePrimeNodePtr;        // f

  // V_s Row:
    double * f_SourceEquSourceNodePtr;           // g
    double * f_SourceEquSourcePrimeNodePtr;      // h

  // V_b Row:
    double * f_BulkEquGateNodePtr;               // i
    double * f_BulkEquBulkNodePtr;               // j
    double * f_BulkEquDrainPrimeNodePtr;         // k
    double * f_BulkEquSourcePrimeNodePtr;        // l

  // V_d' Row:
    double * f_DrainPrimeEquDrainNodePtr;        // m
    double * f_DrainPrimeEquGateNodePtr;         // n
    double * f_DrainPrimeEquBulkNodePtr;         // o
    double * f_DrainPrimeEquDrainPrimeNodePtr;   // p
    double * f_DrainPrimeEquSourcePrimeNodePtr;  // q

  // V_s' Row:
    double * f_SourcePrimeEquGateNodePtr;        // r
    double * f_SourcePrimeEquSourceNodePtr;      // s
    double * f_SourcePrimeEquBulkNodePtr;        // t
    double * f_SourcePrimeEquDrainPrimeNodePtr;  // u
    double * f_SourcePrimeEquSourcePrimeNodePtr; // v

  // Q-vector pointers:
  // V_d Row:
    double * q_DrainEquDrainNodePtr;             // a
    double * q_DrainEquDrainPrimeNodePtr;        // b

  // V_g Row:
    double * q_GateEquGateNodePtr;               // c
    double * q_GateEquBulkNodePtr;               // d
    double * q_GateEquDrainPrimeNodePtr;         // e
    double * q_GateEquSourcePrimeNodePtr;        // f

  // V_s Row:
    double * q_SourceEquSourceNodePtr;           // g
    double * q_SourceEquSourcePrimeNodePtr;      // h

  // V_b Row:
    double * q_BulkEquGateNodePtr;               // i
    double * q_BulkEquBulkNodePtr;               // j
    double * q_BulkEquDrainPrimeNodePtr;         // k
    double * q_BulkEquSourcePrimeNodePtr;        // l

  // V_d' Row:
    double * q_DrainPrimeEquDrainNodePtr;        // m
    double * q_DrainPrimeEquGateNodePtr;         // n
    double * q_DrainPrimeEquBulkNodePtr;         // o
    double * q_DrainPrimeEquDrainPrimeNodePtr;   // p
    double * q_DrainPrimeEquSourcePrimeNodePtr;  // q

  // V_s' Row:
    double * q_SourcePrimeEquGateNodePtr;        // r
    double * q_SourcePrimeEquSourceNodePtr;      // s
    double * q_SourcePrimeEquBulkNodePtr;        // t
    double * q_SourcePrimeEquDrainPrimeNodePtr;  // u
    double * q_SourcePrimeEquSourcePrimeNodePtr; // v
#endif

    ////////////////////////////////////////////////////////////////////
    // 3f5 State Variables & related quantities:
    // voltage drops
    double vbd;
    double vbs;
    double vgs;
    double vds;

    // "original" versions of various voltage drop variables:
    double vgs_orig;
    double vds_orig;
    double vbs_orig;
    double vbd_orig;
    double vgd_orig;

    // "old" versions of various voltage drop variables:
    double vgs_old;
    double vds_old;
    double vbs_old;
    double vbd_old;
    double vgd_old;


    // meyer capacitances
    //gate-source capacitor
    double capgs; //value
    double qgs;   // charge
    // gate-drain capacitor
    double capgd; //value
    double qgd;   //charge
    //gate-bulk capacitor
    double capgb; //value
    double qgb;   //charge

    // diode capacitances
    double qbd; // bulk-drain capacitor charge
    double qbs;  // bulk-source capacitor charge

    // indices into the state vector.
    int li_store_vbd;
    int li_store_vbs;
    int li_store_vgs;
    int li_store_vds;
    int li_store_von;

    // place in store vec for lead currents.
    int li_store_dev_id;
    int li_store_dev_ig;
    int li_store_dev_is;
    int li_store_dev_ib;
    
    int li_state_capgs;
    int li_state_capgd;
    int li_state_capgb;

    int li_state_qgs;
    int li_state_qgd;
    int li_state_qgb;

    int li_state_qbd;
    int li_state_qbs;

    int blockHomotopyID; // For homotopy
    double randomPerturb; // For homotopy
};


//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
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

    int dtype;                  // device type : 1 = nmos,  -1 = pmos
    double model_l;             // the length of the channel region
    double model_w;             // the width of the channel region
    double tnom;                // temperature at which parameters measured
    double latDiff;
    double jctSatCurDensity;    // input - use tSatCurDens
    double jctSatCur;           // input - use tSatCur instead
    double drainResistance;
    double sourceResistance;
    double sheetResistance;
    double transconductance;    // input - use tTransconductance
    double gateSourceOverlapCapFactor;
    double gateDrainOverlapCapFactor;
    double gateBulkOverlapCapFactor;
    double oxideCapFactor;
    double vt0;                 // input - use tVto
    double capBD;               // input - use tCbs
    double capBS;               // input - use tCbd
    double bulkCapFactor;       // input - use tCj
    double sideWallCapFactor;   // input - use tCjsw
    double bulkJctPotential;    // input - use tBulkPot
    double bulkJctBotGradingCoeff;
    double bulkJctSideGradingCoeff;
    double fwdCapDepCoeff;
    double phi;                 // input - use tPhi
    double gamma;
    double substrateDoping;
    int gateType;
    double surfaceStateDensity;
    double oxideThickness;
    double surfaceMobility;     // input - use tSurfMob
    double surfaceMobility0;
    double eta;
    double junctionDepth;
    double coeffDepLayWidth;    // xd
    double narrowFactor;        // delta
    double delta;               // input delta
    double fastSurfaceStateDensity; // nfs
    double theta;               // theta
    double maxDriftVel;         // vmax
    double alpha;               // alpha
    double kappa;               // kappa
    double fNcoef;
    double fNexp;

    bool capBDGiven ;
    bool capBSGiven ;
    bool bulkCapFactorGiven ;
    bool sideWallCapFactorGiven   ;

    // These variables were used as temporaries in mos3temp, but since
    // the calculations in mos3temp are split between the model block
    // constructor and the function updateTemperature, we need them to be
    // model variables.

    double fact1;
    double vtnom;
    double egfet1;
    double pbfact1;
};

//-----------------------------------------------------------------------------
// Function      : Instance:isConverged ()
// Purpose       : Return whether a MOSFET device has done something that
//                  should be interpreted as invalidating other convergence
//                  tests
//                  In case of mos3, just do it if the limiter function
//                  pnjlim changed anything.
//                  This actually agrees with how the Check flag
//                  is used in Spice3F5 mos3load.c
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

    // new DAE stuff:
    // new DAE load functions, residual:
    virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec, double * storeLeadF, double * storeLeadQ);

    // new DAE load functions, Jacobian:
    virtual bool loadDAEMatrices (N_LAS_Matrix & dFdx, N_LAS_Matrix & dQdx);

    friend class Instance;
    friend class Model;
};

} // namespace MOSFET3
} // namespace Device
} // namespace Xyce

typedef Xyce::Device::MOSFET3::Instance N_DEV_MOSFET3Instance;
typedef Xyce::Device::MOSFET3::Model N_DEV_MOSFET3Model;
typedef Xyce::Device::MOSFET3::Master N_DEV_MOSFET3Master;

#endif
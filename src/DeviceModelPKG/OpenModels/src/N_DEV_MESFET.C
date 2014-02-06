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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_MESFET.C,v $
//
// Purpose        :
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
// Revision Number: $Revision: 1.92.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:38 $
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

#ifdef HAVE_CSTDIO
#include <cstdio>
#else
#include <stdio.h>
#endif

// ----------   Xyce Includes   ----------
#include <N_DEV_Const.h>
#include <N_DEV_MESFET.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_MatrixLoadData.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>

namespace Xyce {
namespace Device {

template<>
ParametricData<MESFET::Instance>::ParametricData()
{
    // Set up configuration constants:
    setNumNodes(3);
    setNumOptionalNodes(0);
    setNumFillNodes(0);
    setModelRequired(1);
    addModelType("NMF");
    addModelType("PMF");

    // Set up double precision variables:
    addPar ("TEMP", 0.0, false, ParameterType::TIME_DEP,
      &MESFET::Instance::temp,
      NULL, STANDARD, CAT_NONE, "");

    addPar ("AREA", 1.0, false,   ParameterType::NO_DEP,
      &MESFET::Instance::area,
      NULL, U_METER2, CAT_GEOMETRY, "device area");
}

template<>
ParametricData<MESFET::Model>::ParametricData()
{
    // Set up double precision variables:
    addPar ("AF", 1.0, false, ParameterType::NO_DEP,
      &MESFET::Model::AF,
      NULL, U_NONE, CAT_FLICKER, "Flicker noise exponent");

    addPar ("B", 0.3, false, ParameterType::NO_DEP,
      &MESFET::Model::B,
      NULL, U_VOLTM1, CAT_PROCESS, "Doping tail parameter");

    addPar ("BETA", 2.5e-3, false, ParameterType::NO_DEP,
      &MESFET::Model::BETA,
      NULL, U_AMPVM2, CAT_PROCESS, "Transconductance parameter");

    addPar ("ALPHA", 2.0, false, ParameterType::NO_DEP,
      &MESFET::Model::ALPHA,
      NULL, U_VOLTM1, CAT_PROCESS, "Saturation voltage parameter");

    addPar ("CGS", 0.0, false, ParameterType::MIN_CAP,
      &MESFET::Model::CGS,
      NULL, U_FARAD, CAT_CAP, "Zero-bias gate-source junction capacitance");

    addPar ("CGD", 0.0, false, ParameterType::MIN_CAP,
      &MESFET::Model::CGD,
      NULL, U_FARAD, CAT_CAP, "Zero-bias gate-drain junction capacitance");

    addPar ("FC", 0.5, false, ParameterType::NO_DEP,
      &MESFET::Model::FC,
      NULL, U_FARAD, CAT_CAP, "Coefficient for forward-bias depletion capacitance");

    addPar ("IS", 1e-14, false, ParameterType::NO_DEP,
      &MESFET::Model::IS,
      NULL, U_AMP, CAT_CURRENT, "Gate junction saturation current");

    addPar ("KF", 0.05, false, ParameterType::NO_DEP,
      &MESFET::Model::KF,
      NULL, U_NONE, CAT_FLICKER, "Flicker noise coefficient");

    addPar ("LAMBDA",0.0, false, ParameterType::NO_DEP,
      &MESFET::Model::LAMBDA,
      NULL, U_VOLTM1, CAT_VOLT, "Channel length modulation");

    addPar ("PB", 1.0, false, ParameterType::NO_DEP,
      &MESFET::Model::PB,
      NULL, U_VOLT, CAT_VOLT, "Gate junction potential");

    addPar ("RD", 0.0, false, ParameterType::MIN_RES,
      &MESFET::Model::RD,
      NULL, U_OHM, CAT_RES, "Drain ohmic resistance");

    addPar ("RS", 0.0, false, ParameterType::MIN_RES,
      &MESFET::Model::RS,
      NULL, U_OHM, CAT_RES, "Source ohmic resistance");

    addPar ("TNOM", 0.0, false, ParameterType::NO_DEP,
      &MESFET::Model::TNOM,
      NULL, STANDARD, CAT_NONE, "");

    addPar ("VTO", 0.0, false, ParameterType::NO_DEP,
      &MESFET::Model::VTO,
      NULL, U_VOLT, CAT_VOLT, "Threshold voltage");

    DeviceModel::initThermalModel(*this);
}

namespace MESFET {

vector< vector<int> > Instance::jacStamp_DC_SC;
vector< vector<int> > Instance::jacStamp_DC;
vector< vector<int> > Instance::jacStamp_SC;
vector< vector<int> > Instance::jacStamp;

vector<int> Instance::jacMap_DC_SC;
vector<int> Instance::jacMap_DC;
vector<int> Instance::jacMap_SC;
vector<int> Instance::jacMap;

vector< vector<int> > Instance::jacMap2_DC_SC;
vector< vector<int> > Instance::jacMap2_DC;
vector< vector<int> > Instance::jacMap2_SC;
vector< vector<int> > Instance::jacMap2;



ParametricData<Instance> &Instance::getParametricData() {
  static ParametricData<Instance> parMap;

  return parMap;
}

ParametricData<Model> &Model::getParametricData() {
  static ParametricData<Model> parMap;

  return parMap;
}

//------------------- Class Model ---------------------------------
//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : model block constructor
// Special Notes :
// Scope         : public
// Creator       : pmc
// Creation Date : 11/16/2003
//-----------------------------------------------------------------------------
Model::Model (const ModelBlock & MB,
                                        SolverState & ss1,
                                        DeviceOptions & do1)
  : DeviceModel(MB, ss1,do1),
        AF(1.0),
        B(0.3),
        ALPHA(2.0),
        BETA(2.5e-3),
        CGS(0.0),
        CGD(0.0),
        FC(0.5),
        IS(1.0e-14),
        KF(0.0),
        LAMBDA(0.0),
        PB(1.0),
        RD(0.0),
        RS(0.0),
        TNOM(CONSTREFTEMP),
        VTO(-2.0),
        fNcoef(0.0),
        fNexp(1.0),
        dtype(CONSTNMOS)
{
  if (getType() != "")
  {
    if (getType() == "NMF") {
      dtype = CONSTNMOS;
    }
    else if (getType() == "PMF") {
      dtype = CONSTPMOS;
    }
    else
    {
      string msg =
        ":: model block constructor:\n";
      msg += "Could not recognize the type for model ";
      msg += getName();
      std::ostringstream oss;
      oss << "Error in " << netlistLocation() << "\n" << msg;
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, oss.str());
    }
  }


  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to .model line and constant defaults from metadata:
  setModParams (MB.params);

  // Set any non-constant parameter defaults:
  if (!given("TNOM"))
    TNOM = getDeviceOptions().tnom;

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
// Creation Date : 3/16/00
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


//-----------------------------------------------------------------------------
// Function      : Model::printOutInstances
// Purpose       : debugging tool.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
std::ostream &Model::printOutInstances(std::ostream &os) const
{
  vector<Instance*>::const_iterator iter;
  vector<Instance*>::const_iterator first = instanceContainer.begin();
  vector<Instance*>::const_iterator last  = instanceContainer.end();

  int i,isize;
  isize = instanceContainer.size();
  os << endl;
  os << "Number of MESFET Instances: " << isize << endl;
  os << "    name     getModelName()  Parameters" << endl;

  for (i=0, iter=first; iter!=last; ++iter, ++i)
  {
    os << "  " << i << ": " << (*iter)->getName() << "\t";
    os << (*iter)->getModelName();
    os << endl;
  }

  os << endl;

  return os;
}

//-----------------------------------------------------------------------------
// Function      : Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : pmc
// Creation Date : 11/16/2003
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
// Creator       : Dave Shirely, PSSI
// Creation Date : 03/23/06
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

//------------------------ Class Instance -------------------------
//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : instance block constructor
// Special Notes :
// Scope         : public
// Creator       : pmc
// Creation Date : 11/16/2003
//-----------------------------------------------------------------------------
Instance::Instance(InstanceBlock & IB,
    Model & Miter,
    MatrixLoadData & mlData1,
    SolverState &ss1,
    ExternData  &ed1,
    DeviceOptions & do1)
  : DeviceInstance(IB, mlData1, ss1, ed1, do1),
    model_(Miter),
    limitedFlag(false),
    off(0),
    ic(0),
    area(1.0),
    ic_vds(0.0),
    ic_vgs(0.0),
    temp(getDeviceOptions().temp.dVal()),
    sourceCond(0.0),
    drainCond(0.0),
    tCGS(0.0),
    tCGD(0.0),
    tIS(0.0),
    tPB(0.0),
    tBeta(0.0),
    tvt0(0.0),
    tLambda(0.0),
    tAlpha(0.0),
    tRD(0.0),
    tRS(0.0),
    tMESb(0.0),
    Bfac(0.0),
    dNode(0),
    gNode(0),
    sNode(0),
    dpNode(0),
    spNode(0),
    Vgs(0.0),
    Vgd(0.0),
    gm(0.0),
    gds(0.0),
    ggs(0.0),
    ggd(0.0),
    p(0.0),
    // Solution variables and intermediate quantities
    // drain,source,gate, drainprime and sourceprime voltages
    Vd(0.0),
    Vs(0.0),
    Vg(0.0),
    Vdp(0.0),
    Vsp(0.0),
    // vector local indices
    li_Drain(-1),
    li_DrainPrime(-1),
    li_Source(-1),
    li_SourcePrime(-1),
    li_Gate(-1),
  // Jacobian Matrix
   // Jacobian Matrix Offset:
  // V_d Row:
    ADrainEquDrainNodeOffset(-1),
    ADrainEquDrainPrimeNodeOffset(-1),
  // V_g Row:
    AGateEquGateNodeOffset(-1),
    AGateEquDrainPrimeNodeOffset(-1),
    AGateEquSourcePrimeNodeOffset(-1),
  // V_s Row:
    ASourceEquSourceNodeOffset(-1),
    ASourceEquSourcePrimeNodeOffset(-1),
  // V_d' Row:
    ADrainPrimeEquDrainNodeOffset(-1),
    ADrainPrimeEquGateNodeOffset(-1),
    ADrainPrimeEquDrainPrimeNodeOffset(-1),
    ADrainPrimeEquSourcePrimeNodeOffset(-1),
 // V_s' Row:
    ASourcePrimeEquGateNodeOffset(-1),
    ASourcePrimeEquSourceNodeOffset(-1),
    ASourcePrimeEquDrainPrimeNodeOffset(-1),
    ASourcePrimeEquSourcePrimeNodeOffset(-1),

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
   // dFdx Matrix Ptr:
  // V_d Row:
    f_DrainEquDrainNodePtr(0),
    f_DrainEquDrainPrimeNodePtr(0),
  // V_g Row:
    f_GateEquGateNodePtr(0),
    f_GateEquDrainPrimeNodePtr(0),
    f_GateEquSourcePrimeNodePtr(0),
  // V_s Row:
    f_SourceEquSourceNodePtr(0),
    f_SourceEquSourcePrimeNodePtr(0),
  // V_d' Row:
    f_DrainPrimeEquDrainNodePtr(0),
    f_DrainPrimeEquGateNodePtr(0),
    f_DrainPrimeEquDrainPrimeNodePtr(0),
    f_DrainPrimeEquSourcePrimeNodePtr(0),
 // V_s' Row:
    f_SourcePrimeEquGateNodePtr(0),
    f_SourcePrimeEquSourceNodePtr(0),
    f_SourcePrimeEquDrainPrimeNodePtr(0),
    f_SourcePrimeEquSourcePrimeNodePtr(0),

   // dQdx Matrix Ptr:
  // V_d Row:
    q_DrainEquDrainNodePtr(0),
    q_DrainEquDrainPrimeNodePtr(0),
  // V_g Row:
    q_GateEquGateNodePtr(0),
    q_GateEquDrainPrimeNodePtr(0),
    q_GateEquSourcePrimeNodePtr(0),
  // V_s Row:
    q_SourceEquSourceNodePtr(0),
    q_SourceEquSourcePrimeNodePtr(0),
  // V_d' Row:
    q_DrainPrimeEquDrainNodePtr(0),
    q_DrainPrimeEquGateNodePtr(0),
    q_DrainPrimeEquDrainPrimeNodePtr(0),
    q_DrainPrimeEquSourcePrimeNodePtr(0),
 // V_s' Row:
    q_SourcePrimeEquGateNodePtr(0),
    q_SourcePrimeEquSourceNodePtr(0),
    q_SourcePrimeEquDrainPrimeNodePtr(0),
    q_SourcePrimeEquSourcePrimeNodePtr(0),
#endif
    vgs(0.0),
    vgd(0.0),
    vgs_old(0.0),
    vgd_old(0.0),
    vds_old(0.0),
    vgs_orig(0.0),
    vgd_orig(0.0),

    capgs(0.0),
    qgs(0.0),
    cqgs(0.0),
    capgd(0.0),
    qgd(0.0),
    cqgd(0.0),
    mode(1),
    // local indices
    li_store_vgs(-1),
    li_store_vgd(-1),
    li_store_dev_id(-1),
    li_store_dev_ig(-1),
    li_store_dev_is(-1),
    li_state_qgs(-1),
    li_state_gcgs(-1),
    li_state_qgd(-1),
    li_state_gcgd(-1)
{
  setName(IB.getName());
  setModelName(model_.getName());

  numIntVars   = 2;
  numExtVars   = 3;
  numStateVars = 4;
  setNumStoreVars(2);
  numLeadCurrentStoreVars = 3; // lead currents drain, gate and source

  devConMap.resize(3);
  devConMap[0] = 1;
  devConMap[1] = 2;
  devConMap[2] = 1;

  if( jacStamp.empty() )
  {
    // stamp for RS!=0, RD!=0
    jacStamp_DC_SC.resize(5);
    jacStamp_DC_SC[0].resize(2);  // Drain row
    jacStamp_DC_SC[0][0]=0;       // d-d
    jacStamp_DC_SC[0][1]=3;       // d-d'
    jacStamp_DC_SC[1].resize(3);  // Gate row
    jacStamp_DC_SC[1][0]=1;       // g-g
    jacStamp_DC_SC[1][1]=3;       // g-d'
    jacStamp_DC_SC[1][2]=4;       // g-s'
    jacStamp_DC_SC[2].resize(2);  // Source row
    jacStamp_DC_SC[2][0]=2;       // s-s
    jacStamp_DC_SC[2][1]=4;       // s-s'
    jacStamp_DC_SC[3].resize(4);  // Drain' row
    jacStamp_DC_SC[3][0]=0;       // d'-d
    jacStamp_DC_SC[3][1]=1;       // d'-g
    jacStamp_DC_SC[3][2]=3;       // d'-d'
    jacStamp_DC_SC[3][3]=4;       // d'-s'
    jacStamp_DC_SC[4].resize(4);  // Source' row
    jacStamp_DC_SC[4][0]=1;       // s'-g
    jacStamp_DC_SC[4][1]=2;       // s'-s
    jacStamp_DC_SC[4][2]=3;       // s'-d'
    jacStamp_DC_SC[4][3]=4;       // s'-s'

    jacMap_DC_SC.clear();
    jacStampMap(jacStamp_DC_SC, jacMap_DC_SC, jacMap2_DC_SC,
                jacStamp_DC,    jacMap_DC, jacMap2_DC, 4, 2, 5);

    jacStampMap(jacStamp_DC_SC, jacMap_DC_SC, jacMap2_DC_SC,
                jacStamp_SC,    jacMap_SC, jacMap2_SC, 3, 0, 5);

    jacStampMap(jacStamp_DC, jacMap_DC, jacMap2_DC,
                jacStamp,    jacMap, jacMap2, 3, 0, 5);

  }

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  // Set any non-constant parameter defaults:
  if (!given("TEMP"))
    temp = getDeviceOptions().temp.dVal();

  updateDependentParameters();

  // Calculate any parameters specified as expressions:
  processParams ();

  // process source/drain series resistance
  drainCond = 0;
  if (model_.RD != 0)
    drainCond = area/model_.RD;
  sourceCond = 0;
  if (model_.RS != 0)
    sourceCond = area/model_.RS;

  numIntVars = (((sourceCond == 0.0)?0:1)+((drainCond == 0.0)?0:1));

}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
Instance::~Instance ()
{
}

//----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : pmc
// Creation Date : 11/16/2003
//----------------------------------------------------------------------------
void Instance::registerLIDs( const vector<int> & intLIDVecRef,
                                       const vector<int> & extLIDVecRef )
{
  string msg;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
  const string dashedline =
  "-------------------------------------------------------------------------"
  "----";
    cout << endl << dashedline << endl;
    cout << "  Instance::registerLIDs" << endl;
    cout << "  name = " << getName() << endl;
  }
#endif

  // Check if the size of the ID lists corresponds to the
  // proper number of internal and external variables.
  int numInt = intLIDVecRef.size();
  int numExt = extLIDVecRef.size();

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << "  number of internal variables: " << numInt << endl;
    cout << "  number of external variables: " << numExt << endl;
  }
#endif

  numIntVars = (((sourceCond == 0.0)?0:1)+((drainCond == 0.0)?0:1));

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

  // now use these lists to obtain the indices into the
  // linear algebra entities.  This assumes an order.
  // For the matrix  indices, first do the rows.

  li_Drain = extLIDVec[0];
  li_Gate  = extLIDVec[1];
  li_Source = extLIDVec[2];

  int intLoc = 0;

  if( drainCond )
    li_DrainPrime = intLIDVec[intLoc++];
  else
    li_DrainPrime = li_Drain;

  if( sourceCond )
    li_SourcePrime = intLIDVec[intLoc];
  else
    li_SourcePrime = li_Source;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
  const string dashedline =
    "------------------------------------------------------------------------"
    "-----";
    cout << "\n variable local indices:\n";
    cout << "  li_Drain       = " << li_Drain << endl;
    cout << "  li_DrainPrime  = " << li_DrainPrime << endl;
    cout << "  li_Source      = " << li_Source << endl;
    cout << "  li_SourcePrime = " << li_SourcePrime << endl;
    cout << "  li_Gate        = " << li_Gate << endl;

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
// Creation Date : 05/13/05
//-----------------------------------------------------------------------------
map<int,string> & Instance::getIntNameMap ()
{
  // set up the internal name map, if it hasn't been already.
  if (intNameMap.empty ())
  {
    // set up the internal name map:
    string tmpstr;

    if (drainCond != 0.0)
    {
      tmpstr = getName()+"_drainprime";
      spiceInternalName (tmpstr);
      intNameMap[ li_DrainPrime ] = tmpstr;
    }

    if (sourceCond != 0.0)
    {
      tmpstr = getName()+"_sourceprime";
      spiceInternalName (tmpstr);
      intNameMap[ li_SourcePrime ] = tmpstr;
    }
  }

  return intNameMap;
}


//-----------------------------------------------------------------------------
// Function      : N_DEV_MESFETInstance::getStoreNameMap
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 4/4/2013
//-----------------------------------------------------------------------------
map<int,string> & N_DEV_MESFETInstance::getStoreNameMap ()
{
  // set up the internal name map, if it hasn't been already.
  if( loadLeadCurrent && storeNameMap.empty ())
  {
    // change subcircuitname:devicetype_deviceName to
    // devicetype:subcircuitName:deviceName
    string modName(getName());
    spiceInternalName(modName);
    string tmpstr;
    tmpstr = modName+":DEV_ID";
    storeNameMap[ li_store_dev_id ] = tmpstr;
    tmpstr = modName+":DEV_IG";
    storeNameMap[ li_store_dev_ig ] = tmpstr;
    tmpstr = modName+":DEV_IS";
    storeNameMap[ li_store_dev_is ] = tmpstr;
  }

  return storeNameMap;
}


//----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//----------------------------------------------------------------------------
void Instance::registerStateLIDs(const vector<int> & staLIDVecRef)
{
  string msg;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
  const string dashedline =
    "------------------------------------------------------------------------"
    "-----";
    cout << endl;
    cout << dashedline << endl;
    cout << "  In Instance::registerStateLIDs\n\n";
    cout << "  name             = " << getName() << endl;
  }
#endif

  // Check if the size of the ID lists corresponds to the proper number of
  // internal and external variables.
  int numSta = staLIDVecRef.size();

  if (numSta != numStateVars)
  {
    msg = "Instance::registerStateLIDs:";
    msg += "numSta != numStateVars";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
    cout << "  Number of State LIDs: " << numSta << endl;
#endif

  // Copy over the global ID lists:
  staLIDVec = staLIDVecRef;

  int lid=0;
  li_state_qgs  = staLIDVec[lid++];
  li_state_gcgs = staLIDVec[lid++];

  li_state_qgd  = staLIDVec[lid++];
  li_state_gcgd = staLIDVec[lid++];


#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
  const string dashedline =
    "------------------------------------------------------------------------"
    "-----";
    cout << "  State local indices:" << endl;
    cout << endl;

    cout << "  li_state_qgs       = " << li_state_qgs << endl;
    cout << "  li_state_gcgs      = " << li_state_gcgs;
    cout << "  li_state_qgd       = " << li_state_qgd;
    cout << "  li_state_gcgd      = " << li_state_gcgd << endl;;

    cout << dashedline << endl;
  }
#endif

}

//----------------------------------------------------------------------------
// Function      : Instance::registerStoreLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 12/9/11
//----------------------------------------------------------------------------
void Instance::registerStoreLIDs(const vector<int> & stoLIDVecRef)
{
#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
  const string dashedline =
    "------------------------------------------------------------------------"
    "-----";
    cout << endl;
    cout << dashedline << endl;
    cout << "  In Instance::registerStoreLIDs\n\n";
    cout << "  name             = " << getName() << endl;
  }
#endif

  // Check if the size of the ID lists corresponds to the proper number of
  // internal and external variables.
  int numSto = stoLIDVecRef.size();

  if (numSto != getNumStoreVars())
  {
    string msg = "Instance::registerStoreLIDs:";
    msg += "numSto != numStoreVars";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
    cout << "  Number of Store LIDs: " << numSto << endl;
#endif

  // Copy over the global ID lists:
  stoLIDVec = stoLIDVecRef;

  int lid=0;
  li_store_vgs = stoLIDVec[lid++];
  li_store_vgd = stoLIDVec[lid++];

  if( loadLeadCurrent )
  {
    li_store_dev_id = stoLIDVec[lid++];
    li_store_dev_ig = stoLIDVec[lid++];
    li_store_dev_is = stoLIDVec[lid++];
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
  const string dashedline =
    "------------------------------------------------------------------------"
    "-----";
    cout << "  Store local indices:" << endl;
    cout << endl;
    cout << "  li_store_vgs       = " << li_store_vgs;
    cout << "  li_store_vgd       = " << li_store_vgd;
    cout << dashedline << endl;
  }
#endif
}

//----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//----------------------------------------------------------------------------
const vector< vector<int> > & Instance::jacobianStamp() const
{
  if( drainCond != 0.0 && sourceCond != 0.0 )
    return jacStamp_DC_SC;
  else if( drainCond != 0.0 && sourceCond == 0.0 )
    return jacStamp_DC;
  else if( drainCond == 0.0 && sourceCond != 0.0 )
    return jacStamp_SC;
  else if( drainCond == 0.0 && sourceCond == 0.0 )
    return jacStamp;
  else
    return jacStamp;
}

//----------------------------------------------------------------------------
// Function      : Instance::registerJacLIDs
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//----------------------------------------------------------------------------
void Instance::registerJacLIDs( const vector< vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );
  vector<int> map;
  vector< vector<int> > map2;

  if (drainCond != 0.0)
  {
    if (sourceCond != 0.0)
    {
      map = jacMap_DC_SC;
      map2 = jacMap2_DC_SC;
    }
    else
    {
      map = jacMap_DC;
      map2 = jacMap2_DC;
    }
  }
  else
  {
    if (sourceCond != 0.0)
    {
      map = jacMap_SC;
      map2 = jacMap2_SC;
    }
    else
    {
      map = jacMap;
      map2 = jacMap2;
    }
  }

  ADrainEquDrainNodeOffset             = jacLIDVec[map[0]][map2[0][0]];
  ADrainEquDrainPrimeNodeOffset        = jacLIDVec[map[0]][map2[0][1]];

  AGateEquGateNodeOffset               = jacLIDVec[map[1]][map2[1][0]];
  AGateEquDrainPrimeNodeOffset         = jacLIDVec[map[1]][map2[1][1]];
  AGateEquSourcePrimeNodeOffset        = jacLIDVec[map[1]][map2[1][2]];

  ASourceEquSourceNodeOffset           = jacLIDVec[map[2]][map2[2][0]];
  ASourceEquSourcePrimeNodeOffset      = jacLIDVec[map[2]][map2[2][1]];

  ADrainPrimeEquDrainNodeOffset        = jacLIDVec[map[3]][map2[3][0]];
  ADrainPrimeEquGateNodeOffset         = jacLIDVec[map[3]][map2[3][1]];
  ADrainPrimeEquDrainPrimeNodeOffset   = jacLIDVec[map[3]][map2[3][2]];
  ADrainPrimeEquSourcePrimeNodeOffset  = jacLIDVec[map[3]][map2[3][3]];

  ASourcePrimeEquGateNodeOffset        = jacLIDVec[map[4]][map2[4][0]];
  ASourcePrimeEquSourceNodeOffset      = jacLIDVec[map[4]][map2[4][1]];
  ASourcePrimeEquDrainPrimeNodeOffset  = jacLIDVec[map[4]][map2[4][2]];
  ASourcePrimeEquSourcePrimeNodeOffset = jacLIDVec[map[4]][map2[4][3]];
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupPointers
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/30/08
//-----------------------------------------------------------------------------
void Instance::setupPointers ()
{
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  N_LAS_Matrix & dFdx = *(extData.dFdxMatrixPtr);
  N_LAS_Matrix & dQdx = *(extData.dQdxMatrixPtr);

  // F-matrix:
  f_DrainEquDrainNodePtr             = 	&(dFdx[li_Drain][ADrainEquDrainNodeOffset]);
  f_DrainEquDrainPrimeNodePtr        = 	&(dFdx[li_Drain][ADrainEquDrainPrimeNodeOffset]);

  f_GateEquGateNodePtr               = 	&(dFdx[li_Gate][AGateEquGateNodeOffset]);
  f_GateEquDrainPrimeNodePtr         = 	&(dFdx[li_Gate][AGateEquDrainPrimeNodeOffset]);
  f_GateEquSourcePrimeNodePtr        = 	&(dFdx[li_Gate][AGateEquSourcePrimeNodeOffset]);

  f_SourceEquSourceNodePtr           = 	&(dFdx[li_Source][ASourceEquSourceNodeOffset]);
  f_SourceEquSourcePrimeNodePtr      = 	&(dFdx[li_Source][ASourceEquSourcePrimeNodeOffset]);

  f_DrainPrimeEquDrainNodePtr        = 	&(dFdx[li_DrainPrime][ADrainPrimeEquDrainNodeOffset]);
  f_DrainPrimeEquGateNodePtr         = 	&(dFdx[li_DrainPrime][ADrainPrimeEquGateNodeOffset]);
  f_DrainPrimeEquDrainPrimeNodePtr   = 	&(dFdx[li_DrainPrime][ADrainPrimeEquDrainPrimeNodeOffset]);
  f_DrainPrimeEquSourcePrimeNodePtr  = 	&(dFdx[li_DrainPrime][ADrainPrimeEquSourcePrimeNodeOffset]);

  f_SourcePrimeEquGateNodePtr        = 	&(dFdx[li_SourcePrime][ASourcePrimeEquGateNodeOffset]);
  f_SourcePrimeEquSourceNodePtr      = 	&(dFdx[li_SourcePrime][ASourcePrimeEquSourceNodeOffset]);
  f_SourcePrimeEquDrainPrimeNodePtr  = 	&(dFdx[li_SourcePrime][ASourcePrimeEquDrainPrimeNodeOffset]);
  f_SourcePrimeEquSourcePrimeNodePtr = 	&(dFdx[li_SourcePrime][ASourcePrimeEquSourcePrimeNodeOffset]);

  // Q-matrix:
  q_DrainEquDrainNodePtr             = 	&(dQdx[li_Drain][ADrainEquDrainNodeOffset]);
  q_DrainEquDrainPrimeNodePtr        = 	&(dQdx[li_Drain][ADrainEquDrainPrimeNodeOffset]);

  q_GateEquGateNodePtr               = 	&(dQdx[li_Gate][AGateEquGateNodeOffset]);
  q_GateEquDrainPrimeNodePtr         = 	&(dQdx[li_Gate][AGateEquDrainPrimeNodeOffset]);
  q_GateEquSourcePrimeNodePtr        = 	&(dQdx[li_Gate][AGateEquSourcePrimeNodeOffset]);

  q_SourceEquSourceNodePtr           = 	&(dQdx[li_Source][ASourceEquSourceNodeOffset]);
  q_SourceEquSourcePrimeNodePtr      = 	&(dQdx[li_Source][ASourceEquSourcePrimeNodeOffset]);

  q_DrainPrimeEquDrainNodePtr        = 	&(dQdx[li_DrainPrime][ADrainPrimeEquDrainNodeOffset]);
  q_DrainPrimeEquGateNodePtr         = 	&(dQdx[li_DrainPrime][ADrainPrimeEquGateNodeOffset]);
  q_DrainPrimeEquDrainPrimeNodePtr   = 	&(dQdx[li_DrainPrime][ADrainPrimeEquDrainPrimeNodeOffset]);
  q_DrainPrimeEquSourcePrimeNodePtr  = 	&(dQdx[li_DrainPrime][ADrainPrimeEquSourcePrimeNodeOffset]);

  q_SourcePrimeEquGateNodePtr        = 	&(dQdx[li_SourcePrime][ASourcePrimeEquGateNodeOffset]);
  q_SourcePrimeEquSourceNodePtr      = 	&(dQdx[li_SourcePrime][ASourcePrimeEquSourceNodeOffset]);
  q_SourcePrimeEquDrainPrimeNodePtr  = 	&(dQdx[li_SourcePrime][ASourcePrimeEquDrainPrimeNodeOffset]);
  q_SourcePrimeEquSourcePrimeNodePtr = 	&(dQdx[li_SourcePrime][ASourcePrimeEquSourcePrimeNodeOffset]);

#endif
}

//----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : pmc
// Creation Date : 11/16/2003
//----------------------------------------------------------------------------
bool Instance::updatePrimaryState ()
{
  double * staVec = extData.nextStaVectorRawPtr;
  bool bsuccess = updateIntermediateVars ();

  double * stoVec = extData.nextStoVectorRawPtr;
  stoVec[li_store_vgs] = vgs;
  stoVec[li_store_vgd] = vgd;
  staVec[li_state_qgs] = qgs;
  staVec[li_state_qgd] = qgd;

  return  bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : pmc
// Creation Date : 11/16/2003
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars ()
{
  double * solVec = extData.nextSolVectorRawPtr;
  double * currStaVec = extData.currStaVectorRawPtr;

  int    dtype;
  double csat, betap;
  double vgst, vgdt;
  double evgs, evgd;
  double sarg, vtf;;

// from the spice jfet
  double czgd, czgs;
  double czgdf2, czgsf2;
  double fcpb2;
  double twop;
  int    icheck, ichk1;

// for the Shockley version
//  double A, B, C, B12, C12, D, Vdsat;
//  double delta;
  double prod, denom, invdenom, afact, lfact;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    const string dashedline2 = "---------------------";
    cout << dashedline2 << endl;
    cout <<"  Instance::updateIntermediateVars.\n"<<endl;
    cout <<"  name = " << getName() << endl;
    cout <<"  Model name = " << model_.getName() << endl;
    cout <<"  dtype is " << model_.dtype << endl;
    cout << endl;
    cout.width(25); cout.precision(17); cout.setf(ios::scientific);
  }
#endif

  icheck = 1;
  dtype  = model_.dtype;

  //  we need our solution variables for any of this stuff
  Vd  = 0.0;
  Vs  = 0.0;
  Vg  = 0.0;
  Vdp = 0.0;
  Vsp = 0.0;

  Vd  = solVec[li_Drain];
  Vg  = solVec[li_Gate];
  Vs  = solVec[li_Source];
  Vsp = solVec[li_SourcePrime];
  Vdp = solVec[li_DrainPrime];

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << " " << endl;
    cout << " Vg  = " << Vg << endl;
    cout << " Vd  = " << Vd << endl;
    cout << " Vs  = " << Vs << endl;
    cout << " Vdp = " << Vdp << endl;
    cout << " Vsp = " << Vsp << endl;
  }
#endif

  // now we need voltage drops
  Vddp  = Vd - Vdp;
  Vssp  = Vs - Vsp;
  Vgsp  = Vg - Vsp;
  Vgdp  = Vg - Vdp;
  Vdpsp = Vdp - Vsp;

  // Now the things that the 3f5 code really uses
  vgs = dtype * Vgsp;
  vgd = dtype * Vgdp;
  vds = vgs-vgd;

  origFlag = 1;
  limitedFlag = false;
  vgs_orig = vgs;
  vgd_orig = vgd;
  vds_orig = vds;

  if (getSolverState().newtonIter == 0)
  {
    newtonIterOld=0;
    if (getSolverState().initJctFlag && getDeviceOptions().voltageLimiterFlag)
    {
      if (getSolverState().inputOPFlag)
      {
        N_LAS_Vector * flagSolVectorPtr = extData.flagSolVectorPtr;
        if ((*flagSolVectorPtr)[li_Drain] == 0 || (*flagSolVectorPtr)[li_Gate] == 0 ||
            (*flagSolVectorPtr)[li_Source] == 0 || (*flagSolVectorPtr)[li_SourcePrime] ||
            (*flagSolVectorPtr)[li_DrainPrime] )
        {
          vgs = 0;
          vgd = 0;
          vds = vgs-vgd;
        }
      }
      else
      {
        vgs = 0;
        vgd = 0;
        vds = vgs-vgd;
      }
    }
    if (!(getSolverState().dcopFlag)||(getSolverState().locaEnabledFlag && getSolverState().dcopFlag))
    {
      double * currStoVec = extData.currStoVectorRawPtr;
      vgs_old = currStoVec[li_store_vgs];
      vgd_old = currStoVec[li_store_vgd];
    }
    else
    { // there is no history
      vgs_old = vgs;
      vgd_old = vgd;
    }
  }
  else
  {
    double *stoVec = extData.nextStoVectorRawPtr;
    vgs_old = stoVec[li_store_vgs];
    vgd_old = stoVec[li_store_vgd];
  }

  // SPICE-type Voltage Limiting
  ///////////////////////////////

  if (getDeviceOptions().voltageLimiterFlag)
  {
#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      cout << " before limiting: " << endl;
      cout << " vgs = " << vgs <<  "   vgs_old = " << vgs_old << endl;
      cout << " vgd = " << vgd <<  "   vgd_old = " << vgd_old << endl;
    }
#endif

    ichk1=1;
    vgs = devSupport.pnjlim(vgs, vgs_old, vt, vcrit, &icheck);
    vgd = devSupport.pnjlim(vgd, vgd_old, vt, vcrit, &ichk1);

    if (ichk1 == 1) {icheck=1;}
    if (icheck == 1) limitedFlag=true;

    vgs = devSupport.fetlim(vgs, vgs_old, tvt0);
    vgd = devSupport.fetlim(vgd, vgd_old, tvt0);
    vds = vgs-vgd;

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      cout << " After limiting: " << endl;
      cout << " vgs = " << vgs << endl;
      cout << " vgd = " << vgd << endl;
      cout << " " << endl;
    }
#endif
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << "vgs   = " << vgs << endl;
    cout << "vgd   = " << vgd << endl;
    cout << "vds   = " << vds << endl;
    cout << "Vddp  = " << Vddp << endl;
    cout << "Vssp  = " << Vssp << endl;
    cout << "Vgsp  = " << Vgsp << endl;
    cout << "Vgdp  = " << Vgdp << endl;
    cout << "Vdpsp = " << Vdpsp << endl;
    cout << " " << endl;
  }
#endif
  // Now set the origFlag
  if (vgs_orig != vgs || vds_orig != vds || vgd_orig != vgd) origFlag = 0;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    if (origFlag == 0)
    {
      cout << " Something modified the voltages. " << endl;
      cout << " Voltage       before                    after                   diff " << endl;
      cout << " vgs  " << vgs_orig << "  " << vgs << " " << vgs-vgs_orig << endl;
      cout << " vgd  " << vgd_orig << "  " << vgd << " " << vgd-vgd_orig << endl;
      cout << " vds  " << vds_orig << "  " << vds << " " << vds-vds_orig << endl;
      cout << " " << endl;
    }
  }
#endif

  // update the "old" variables:
  if (getSolverState().newtonIter != 0 && getSolverState().newtonIter != newtonIterOld)
  {
    newtonIterOld = getSolverState().newtonIter;
  }

  //
  //  the following block of code evaluates the dc current and its
  //  derivatives and the charges associated with the gate and
  //  channel
  //

  // vt set in updateTemperature
  vtf = 5.0*vt;
  csat = tIS;
  if (vgs <= -vtf)
  {
    ggs = -csat/vgs + getDeviceOptions().gmin;
    cg  = ggs*vgs;
  }
  else
  {
    evgs = exp(vgs/vt);
    ggs = csat*evgs/vt + getDeviceOptions().gmin;
    cg = csat*(evgs-1) + getDeviceOptions().gmin*vgs;
  }
  if (vgd <= -vtf)
  {
    ggd = -csat/vgd + getDeviceOptions().gmin;
    cgd = ggd*vgd;
  }
  else
  {
    evgd = exp(vgd/vt);
    ggd = csat*evgd/vt  + getDeviceOptions().gmin;
    cgd = csat*(evgd-1) + getDeviceOptions().gmin*vgd;
  }
  cg = cg + cgd;

  // 3f5 does this simple stuff
  if (vds >= 0)
    mode = 1;
  else
    mode = -1;

  if (vds >= 0)  // normal mode
  {
    vgst = vgs-tvt0;
    if (vgst <= 0)
    {
      //
      //   normal mode, cutoff region
      //
      cdrain = 0;
      gm = 0;
      gds = 0;
    }
    else
    {
      prod = 1 + tLambda*vds;
      betap = tBeta*prod;
      denom = 1 + tMESb*vgst;
      invdenom = 1/denom;
      if (vds >= ( 3/tAlpha ) )
      {
        //
        //   normal mode, saturation region
        //
        cdrain = betap*vgst*vgst*invdenom;
        gm = betap*vgst*(1 + denom)*invdenom*invdenom;
        gds = tLambda*tBeta*vgst*vgst*invdenom;
      }
      else
      {
        //
        //   normal mode, linear region
        //
        afact = 1 - tAlpha*vds/3;
        lfact = 1 - afact*afact*afact;
        cdrain = betap*vgst*vgst*invdenom*lfact;
        gm = betap*vgst*(1 + denom)*invdenom*invdenom*lfact;
        gds = tBeta*vgst*vgst*invdenom*(tAlpha*afact*afact*prod + lfact*tLambda);
      }
    }
  }
  else   // inverse mode
  {
    vgdt = vgd - tvt0;
    if (vgdt <= 0)
    {
      //
      //   inverse mode, cutoff region
      //
      cdrain = 0;
      gm = 0;
      gds = 0;
    }
    else
    {
      //
      //   inverse mode, saturation region
      //
      prod = 1 - tLambda*vds;
      betap = tBeta*prod;
      denom = 1 + tMESb*vgdt;
      invdenom = 1/denom;
      if ( -vds >= 3/tAlpha )
      {
        cdrain = -betap*vgdt*vgdt*invdenom;
        gm = -betap*vgdt*(1 + denom)*invdenom*invdenom;
        gds = tLambda*tBeta*vgdt*vgdt*invdenom - gm;
      }
      else
      {
        //
        //  inverse mode, linear region
        //
        afact = 1 + tAlpha*vds/3;
        lfact = 1 - afact*afact*afact;
        cdrain = -betap*vgdt*vgdt*invdenom*lfact;
        gm = -betap*vgdt*(1 + denom)*invdenom*invdenom*lfact;
        gds = tBeta*vgdt*vgdt*invdenom*(tAlpha*afact*afact*prod
                                                    + lfact*tLambda) - gm;
      }
    }
  }
  cd = cdrain-cgd;

  //
  //     charge storage elements
  //
  twop  = 2.0*tPB;
  fcpb2 = corDepCap*corDepCap;
  czgs  = tCGS;
  czgd  = tCGD;
  if(czgs != 0)
  {
    czgsf2=czgs/f2;
    if (vgs < corDepCap)
    {
      sarg=sqrt(1-vgs/tPB);
      qgs = twop*czgs*(1-sarg);
      capgs=czgs/sarg;
    }
    else
    {
      qgs = czgs*f1 + czgsf2*(f3 *(vgs - corDepCap)
           +(vgs*vgs - fcpb2)/(2*twop));
      capgs=czgsf2*(f3 + vgs/twop);
    }
  }
  else
  {
    qgs=0.0;
    capgs=0.0;
  }

  if(czgd != 0)
  {
    czgdf2=czgd/f2;
    if (vgd < corDepCap)
    {
      sarg=sqrt(1-vgd/tPB);
      qgd = twop*czgd*(1-sarg);
      capgd=czgd/sarg;
    }
    else
    {
      qgd = czgd*f1 + czgdf2*( f3*(vgd - corDepCap)
          +(vgd*vgd - fcpb2)/(2*twop) );
      capgd=czgdf2*(f3 + vgd/twop);
    }
  }
  else
  {
    qgd=0.0;
    capgd=0.0;
  }

  Idrain  = drainCond  * Vddp;
  Isource = sourceCond * Vssp;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << " Done with Instance::updateIntermediateVars." << endl;
    cout << "  mode    = " << mode << endl;
    cout << "  tBeta   = " << tBeta << endl;
    cout << "  Idrain  = " << Idrain << endl;
    cout << "  Isource = " << Isource << endl;
    cout << "  gds     = " << gds << endl;
    cout << "  gm      = " << gm << endl;
  }
#endif

  /// CURRENTS to load into RHS:

  // so at this point:

  // current out of drain is
  // Idrain

  // current out of gate:
  // dtype*( d/dt(qgs) + d/dt(qgd) )

  //  the current *out of* the source should be simply
  // Isource

  // current out of drain' is
  // -Idrain - dtype*( d/dt(qgd) -  cdrain )

  // the current out of the source' is
  //  -Isource - dtype*( d/dt(qgs) +  cdrain )

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEQVector
//
// Purpose       : Loads the Q-vector contributions for a single
//                 voltage source instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) + B(t) = 0
//
//                 The "Q" vector contains charges and fluxes, mostly.
//                 The voltage source will not make any contributions to Q,
//                 so this function does nothing.
//
//    from updateSecondaryState:
//
//    ggd = ggd + capgd*(getSolverState().pdt);
//    ggs = ggs + capgs*(getSolverState().pdt);
//
//    // Sum the capacitor currents into the DC currents.
//    cg = cg + cqgs + cqgd;
//    cd  = cd  - cqgd;
//    cgd = cgd + cqgd;
//
//    So:
//
//    replace ggd with capgd.
//    replace ggs with capgs
//
//    replace cg with qgs+qgd
//    replace cd with -qgd
//    replace cgd with qgd.
//
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 06/01/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEQVector ()
{
  double * qVec = extData.daeQVectorRawPtr;
  double * dQdxdVp = extData.dQdxdVpVectorRawPtr;

  // set up the final load variables:
  int Dtype = model_.dtype;
  double ceqgd = Dtype*(qgd);
  double ceqgs = Dtype*(((qgs+qgd)-qgd));
  double cdreq = Dtype*(((-qgd)+qgd));

  double ceqgd_Jdxp = -Dtype*(capgd*(vgd-vgd_orig));
  double ceqgs_Jdxp = -Dtype*(capgs*(vgs-vgs_orig));
  double cdreq_Jdxp = 0.0;

  qVec[li_Gate       ] += ( ceqgs+ceqgd);
  qVec[li_DrainPrime ] -= (-cdreq+ceqgd);
  qVec[li_SourcePrime] -= ( cdreq+ceqgs);

  if (!origFlag)
  {
    dQdxdVp[li_Gate       ] -= ( ceqgs_Jdxp+ceqgd_Jdxp);
    dQdxdVp[li_DrainPrime ] += (-cdreq_Jdxp+ceqgd_Jdxp);
    dQdxdVp[li_SourcePrime] += ( cdreq_Jdxp+ceqgs_Jdxp);
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 MESFET instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 06/01/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  double * fVec = extData.daeFVectorRawPtr;
  double * dFdxdVp = extData.dFdxdVpVectorRawPtr;

  // set up the final load variables:
  int Dtype = model_.dtype;
  double ceqgd = Dtype*(cgd);
  double ceqgs = Dtype*((cg-cgd));
  double cdreq = Dtype*((cd+cgd));

  double ceqgd_Jdxp = -Dtype*(ggd*(vgd-vgd_orig));
  double ceqgs_Jdxp = -Dtype*(ggs*(vgs-vgs_orig));
  double cdreq_Jdxp = -Dtype*(gds*(vds-vds_orig)+gm*(vgs-vgs_orig));

  // optional load resistors:
  if (drainCond  != 0.0)
  {
    fVec[li_Drain ] += Idrain;
  }
  if (sourceCond != 0.0)
  {
    fVec[li_Source] += Isource;
  }

  fVec[li_Gate       ] += (ceqgs+ceqgd);
  fVec[li_DrainPrime ] -= (Idrain +(-cdreq+ceqgd));
  fVec[li_SourcePrime] -= (Isource+(cdreq+ceqgs));

  if (!origFlag)
  {
    dFdxdVp[li_Gate       ] -= ( ceqgs_Jdxp+ceqgd_Jdxp);
    dFdxdVp[li_DrainPrime ] += (-cdreq_Jdxp+ceqgd_Jdxp);
    dFdxdVp[li_SourcePrime] += ( cdreq_Jdxp+ceqgs_Jdxp);
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx
//
// Purpose       : Loads the Q-vector contributions for a single
//                 MESFET instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 06/01/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  N_LAS_Matrix & dQdx = *(extData.dQdxMatrixPtr);

  dQdx[li_Gate       ][AGateEquGateNodeOffset        ] += capgd+capgs;
  dQdx[li_Gate       ][AGateEquDrainPrimeNodeOffset        ] -= capgd;
  dQdx[li_Gate       ][AGateEquSourcePrimeNodeOffset       ] -= capgs;
  dQdx[li_DrainPrime ][ADrainPrimeEquGateNodeOffset        ] -= capgd;
  dQdx[li_DrainPrime ][ADrainPrimeEquDrainPrimeNodeOffset  ] += capgd;
  dQdx[li_SourcePrime][ASourcePrimeEquGateNodeOffset       ] -= capgs;
  dQdx[li_SourcePrime][ASourcePrimeEquSourcePrimeNodeOffset] += capgs;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 resistor  instance.
//
// Special Notes : The F-vector is an algebraic constaint.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 06/01/05
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  N_LAS_Matrix & dFdx = *(extData.dFdxMatrixPtr);

  dFdx[li_Drain][ADrainEquDrainNodeOffset] += drainCond;
  dFdx[li_Drain][ADrainEquDrainPrimeNodeOffset] -= drainCond;

  dFdx[li_Gate][AGateEquGateNodeOffset] += ggd+ggs;
  dFdx[li_Gate][AGateEquDrainPrimeNodeOffset] -= ggd;
  dFdx[li_Gate][AGateEquSourcePrimeNodeOffset] -= ggs;

  dFdx[li_Source][ASourceEquSourceNodeOffset] += sourceCond;
  dFdx[li_Source][ASourceEquSourcePrimeNodeOffset] -= sourceCond;

  dFdx[li_DrainPrime][ADrainPrimeEquDrainNodeOffset] -= drainCond;
  dFdx[li_DrainPrime][ADrainPrimeEquGateNodeOffset] += gm-ggd;
  dFdx[li_DrainPrime][ADrainPrimeEquDrainPrimeNodeOffset] +=
     drainCond+gds+ggd;
  dFdx[li_DrainPrime][ADrainPrimeEquSourcePrimeNodeOffset] += -gds-gm;

  dFdx[li_SourcePrime][ASourcePrimeEquGateNodeOffset] -= gm+ggs;
  dFdx[li_SourcePrime][ASourcePrimeEquSourceNodeOffset] -= sourceCond;
  dFdx[li_SourcePrime][ASourcePrimeEquDrainPrimeNodeOffset] -= gds;
  dFdx[li_SourcePrime][ASourcePrimeEquSourcePrimeNodeOffset]
     += sourceCond+gds+gm+ggs;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateTemperature
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
bool Instance::updateTemperature ( const double & temp_tmp)
{
  bool bsuccess = true;
  double tnom, ratio;
  double arg, arg1;
  double ratio1;
  double fact1, fact2;
  double kt, kt1;
  double vtnom;
  double egfet, egfet1;
  double pbfact;
  double cjfact, cjfact1;
  double gmanew, gmaold;
  double pbo;
  double xfc;
  double Pb;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
//    const string dashedline2 = "---------------------";
//    cout << dashedline2 << endl;
    cout << "  Instance::Begin of updateTemperature. \n";
    cout << "  name = " << getName() << endl;
    cout << endl;
  }
#endif

  // first set the instance temperature to the new temperature:
  if (temp_tmp != -999.0) temp = temp_tmp;
  if (model_.interpolateTNOM(temp))
  {
    // make sure interpolation doesn't take any resistance negative
    if(model_.RD < 0) model_.RD = 0;
    if(model_.RS < 0) model_.RS = 0;

    // some params may have changed during interpolation
    // model_.processParams();
  }

  Pb   = model_.PB;
  tnom = model_.TNOM;
  ratio = temp/tnom;

  //  first do the model stuff

  vtnom  = tnom*CONSTKoverQ;
  fact1  = tnom/CONSTREFTEMP;
  kt1    = CONSTboltz*tnom;
  egfet1 = 1.16 - (7.02e-4*tnom*tnom)/(tnom + 1108);
  arg1   = -egfet1/(2.0*kt1) + 1.1150877/(CONSTboltz*2.0*CONSTREFTEMP);
  pbfact = -2.0*vtnom*(1.5*log(fact1) + CONSTQ*arg1);
  pbo    = (Pb - pbfact)/fact1;
  gmaold = (Pb - pbo)/pbo;
  cjfact = 1.0/(1.0 + 0.5*(4e-4*(tnom - CONSTREFTEMP) - gmaold));

  if(model_.FC >.95) {
      cout << "Depletion cap. coeff. FC too large, limited to .95" <<
      cout << endl;
      model_.FC = .95;
  }
  xfc = log(1.0 - model_.FC);
  f2  = exp(1.5*xfc);
  f3  = 1.0 - 1.5*model_.FC;
  // skip  bFac

  //  now do the instance stuff

  vt = temp*CONSTKoverQ;
  kt = temp*CONSTboltz;
  fact2 = temp/CONSTREFTEMP;
  ratio1 = ratio - 1.0;
  tIS = model_.IS*exp(ratio1*1.11/vt)*area;

  tCGS  = model_.CGS*cjfact*area;
  tCGD  = model_.CGD*cjfact*area;
  egfet = 1.16 - (7.02e-4*temp*temp)/(temp + 1108);
  arg   = -egfet/(2.0*kt) + 1.1150877/(CONSTboltz*2.0*CONSTREFTEMP);
  pbfact = -2.0*vt*(1.5*log(fact2) + CONSTQ*arg);
  tPB    = fact2*pbo + pbfact;
  gmanew = (tPB - pbo)/pbo;
  cjfact1 = 1.0 + 0.5*(4e-4*(temp - CONSTREFTEMP) - gmanew);
  tCGS *= cjfact1;
  tCGD *= cjfact1;

  corDepCap = model_.FC*tPB;
  f1    = tPB*(1.0 - exp((0.5)*xfc))/(0.5);
  vcrit = vt * log(vt/(CONSTroot2 * tIS));

  // the following parameters have no temperature dependence in Spice 3f5
  //
  tBeta   =  model_.BETA*area;   // transconductance parameter
  tvt0    =  model_.VTO;         // threshold voltage
  tLambda =  model_.LAMBDA;      // channel-length modulation
  tAlpha  =  model_.ALPHA;       // saturation voltage parameter
  tRD     =  model_.RD/area;     // drain ohmic resistance
  tRS     =  model_.RS/area;     // source ohmic resistance
  tMESb   =  model_.B;           // dopinng tail parameter

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << "temp   = "<< temp << endl;
    cout << "tnom   = " << tnom << endl;
    cout << "ratio  = " << ratio << endl;
    cout << "vt     = " << vt << endl;
    cout << "kt     = " << kt << endl;
    cout << "fact2  = " << fact2 << endl;
    cout << "egfet  = " << egfet << endl;
    cout << "arg    = " << arg << endl;
    cout << "pbfact = " << pbfact << endl;
    cout << "PB     = " << Pb << endl;
    cout << "pbo    = " << pbo << endl;
    cout << "f2     = " << f2 << endl;
    cout << "f3     = " << f3 << endl;
    cout << "corDepCap= " << corDepCap << endl;
    cout << "tBeta   = " << tBeta << endl;
    cout << "tvt0    = " << tvt0 << endl;
    cout << "tPB     = " << tPB << endl;
    cout << "tMESb   = " << tMESb << endl;
    cout << "tLambda = " << tLambda << endl;
    //cout << "tTheta  = " << tTheta << endl;
    cout << "tRD     = " << tRD << endl;
    cout << "tRS     = " << tRS << endl;
    cout << " " << endl;
  }
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
bool Instance::processParams (string param)
{
  updateTemperature(temp);

  return true;
}

// MESFET Master functions:

//-----------------------------------------------------------------------------
// Function      : Master::updateState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
bool Master::updateState (double * solVec, double * staVec, double * stoVec)
{
  bool bsuccess = true;

  for (InstanceVector::const_iterator it = getInstanceVector().begin(); it != getInstanceVector().end(); ++it)
  {
    Instance & ji = *(*it);

    bool btmp = ji.updateIntermediateVars ();
    bsuccess = bsuccess && btmp;

    double * stoVec = ji.extData.nextStoVectorRawPtr;
    stoVec[ji.li_store_vgs] = ji.vgs;
    stoVec[ji.li_store_vgd] = ji.vgd;
    staVec[ji.li_state_qgs] = ji.qgs;
    staVec[ji.li_state_qgd] = ji.qgd;
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Master::loadDAEVectors
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
bool Master::loadDAEVectors (double * solVec, double * fVec, double *qVec,  double * storeLeadF, double * storeLeadQ)
{
  for (InstanceVector::const_iterator it = getInstanceVector().begin(); it != getInstanceVector().end(); ++it)
  {
    Instance & ji = *(*it);

    // F-vector:
    double * dFdxdVp = ji.extData.dFdxdVpVectorRawPtr;

    // set up the final load variables:
    int Dtype = ji.getModel().dtype;
    double f_ceqgd = Dtype*(ji.cgd);
    double f_ceqgs = Dtype*((ji.cg-ji.cgd));
    double f_cdreq = Dtype*((ji.cd+ji.cgd));

    double f_ceqgd_Jdxp = -Dtype*(ji.ggd*(ji.vgd-ji.vgd_orig));
    double f_ceqgs_Jdxp = -Dtype*(ji.ggs*(ji.vgs-ji.vgs_orig));
    double f_cdreq_Jdxp = -Dtype*(ji.gds*(ji.vds-ji.vds_orig)+ji.gm*(ji.vgs-ji.vgs_orig));

    // optional load resistors:
    if (ji.drainCond  != 0.0)
    {
      fVec[ji.li_Drain ] += ji.Idrain;
    }
    if (ji.sourceCond != 0.0)
    {
      fVec[ji.li_Source] += ji.Isource;
    }
    fVec[ji.li_Gate       ] += (f_ceqgs+f_ceqgd);
    fVec[ji.li_DrainPrime ] -= (ji.Idrain +(-f_cdreq+f_ceqgd));
    fVec[ji.li_SourcePrime] -= (ji.Isource+(f_cdreq+f_ceqgs));

    if (!ji.origFlag)
    {
      dFdxdVp[ji.li_Gate       ] -= ( f_ceqgs_Jdxp+f_ceqgd_Jdxp);
      dFdxdVp[ji.li_DrainPrime ] += (-f_cdreq_Jdxp+f_ceqgd_Jdxp);
      dFdxdVp[ji.li_SourcePrime] += ( f_cdreq_Jdxp+f_ceqgs_Jdxp);
    }

    // Q-vector:
    double * dQdxdVp = ji.extData.dQdxdVpVectorRawPtr;

    // set up the final load variables:
    double q_ceqgd = Dtype*(ji.qgd);
    double q_ceqgs = Dtype*(((ji.qgs+ji.qgd)-ji.qgd));
    double q_cdreq = Dtype*(((-ji.qgd)+ji.qgd));

    double q_ceqgd_Jdxp = -Dtype*(ji.capgd*(ji.vgd-ji.vgd_orig));
    double q_ceqgs_Jdxp = -Dtype*(ji.capgs*(ji.vgs-ji.vgs_orig));
    double q_cdreq_Jdxp = 0.0;

    qVec[ji.li_Gate       ] += ( q_ceqgs+q_ceqgd);
    qVec[ji.li_DrainPrime ] -= (-q_cdreq+q_ceqgd);
    qVec[ji.li_SourcePrime] -= ( q_cdreq+q_ceqgs);

    if (!ji.origFlag)
    {
      dQdxdVp[ji.li_Gate       ] -= ( q_ceqgs_Jdxp+q_ceqgd_Jdxp);
      dQdxdVp[ji.li_DrainPrime ] += (-q_cdreq_Jdxp+q_ceqgd_Jdxp);
      dQdxdVp[ji.li_SourcePrime] += ( q_cdreq_Jdxp+q_ceqgs_Jdxp);
    }

    if( ji.loadLeadCurrent )
    {
      if (ji.drainCond != 0.0)
      {
        storeLeadF[ji.li_store_dev_id] = ji.Idrain;
      }
      else
      {
        storeLeadF[ji.li_store_dev_id] = -(ji.Idrain +(-f_cdreq+f_ceqgd));
        storeLeadQ[ji.li_store_dev_id] = -(-q_cdreq+q_ceqgd);
      }
      if (ji.sourceCond != 0.0)
      {
        storeLeadF[ji.li_store_dev_is] = ji.Isource;
      }
      else
      {
        storeLeadF[ji.li_store_dev_is] = -(ji.Isource+(f_cdreq+f_ceqgs));
        storeLeadQ[ji.li_store_dev_is] = -( q_cdreq+q_ceqgs);
      }
      storeLeadF[ji.li_store_dev_ig] = (f_ceqgs+f_ceqgd);
      storeLeadQ[ji.li_store_dev_ig] = (q_ceqgs+q_ceqgd);
    }

  }

  return true;
}

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
//-----------------------------------------------------------------------------
// Function      : Master::loadDAEMatrices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/12/08
//-----------------------------------------------------------------------------
bool Master::loadDAEMatrices (N_LAS_Matrix & dFdx, N_LAS_Matrix & dQdx)
{
  for (InstanceVector::const_iterator it = getInstanceVector().begin(); it != getInstanceVector().end(); ++it)
  {
    Instance & ji = *(*it);

    // F-matrix:

    *ji.f_DrainEquDrainNodePtr += ji.drainCond;

    *ji.f_DrainEquDrainPrimeNodePtr -= ji.drainCond;


    *ji.f_GateEquGateNodePtr += ji.ggd+ji.ggs;

    *ji.f_GateEquDrainPrimeNodePtr -= ji.ggd;

    *ji.f_GateEquSourcePrimeNodePtr -= ji.ggs;


    *ji.f_SourceEquSourceNodePtr += ji.sourceCond;

    *ji.f_SourceEquSourcePrimeNodePtr -= ji.sourceCond;


    *ji.f_DrainPrimeEquDrainNodePtr -= ji.drainCond;

    *ji.f_DrainPrimeEquGateNodePtr += ji.gm-ji.ggd;

    *ji.f_DrainPrimeEquDrainPrimeNodePtr += ji.drainCond+ji.gds+ji.ggd;

    *ji.f_DrainPrimeEquSourcePrimeNodePtr += -ji.gds-ji.gm;


    *ji.f_SourcePrimeEquGateNodePtr -= ji.gm+ji.ggs;

    *ji.f_SourcePrimeEquSourceNodePtr -= ji.sourceCond;

    *ji.f_SourcePrimeEquDrainPrimeNodePtr -= ji.gds;

    *ji.f_SourcePrimeEquSourcePrimeNodePtr += ji.sourceCond+ji.gds+ji.gm+ji.ggs;

    // Q-matrix:

    *ji.q_GateEquGateNodePtr         += ji.capgd+ji.capgs;

    *ji.q_GateEquDrainPrimeNodePtr         -= ji.capgd;

    *ji.q_GateEquSourcePrimeNodePtr        -= ji.capgs;

    *ji.q_DrainPrimeEquGateNodePtr         -= ji.capgd;

    *ji.q_DrainPrimeEquDrainPrimeNodePtr   += ji.capgd;

    *ji.q_SourcePrimeEquGateNodePtr        -= ji.capgs;

    *ji.q_SourcePrimeEquSourcePrimeNodePtr += ji.capgs;
  }

  return true;
}

#else
//-----------------------------------------------------------------------------
// Function      : Master::loadDAEMatrices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 12/12/08
//-----------------------------------------------------------------------------
bool Master::loadDAEMatrices (N_LAS_Matrix & dFdx, N_LAS_Matrix & dQdx)
{
  int sizeInstances = instanceContainer_.size();
  for (int i=0; i<sizeInstances; ++i)
  {
    Instance & ji = *(instanceContainer_.at(i));

    // dFdx matrix:

    dFdx[ji.li_Drain][ji.ADrainEquDrainNodeOffset] += ji.drainCond;

    dFdx[ji.li_Drain][ji.ADrainEquDrainPrimeNodeOffset] -= ji.drainCond;


    dFdx[ji.li_Gate][ji.AGateEquGateNodeOffset] += ji.ggd+ji.ggs;

    dFdx[ji.li_Gate][ji.AGateEquDrainPrimeNodeOffset] -= ji.ggd;

    dFdx[ji.li_Gate][ji.AGateEquSourcePrimeNodeOffset] -= ji.ggs;


    dFdx[ji.li_Source][ji.ASourceEquSourceNodeOffset] += ji.sourceCond;

    dFdx[ji.li_Source][ji.ASourceEquSourcePrimeNodeOffset] -= ji.sourceCond;


    dFdx[ji.li_DrainPrime][ji.ADrainPrimeEquDrainNodeOffset] -= ji.drainCond;

    dFdx[ji.li_DrainPrime][ji.ADrainPrimeEquGateNodeOffset] += ji.gm-ji.ggd;

    dFdx[ji.li_DrainPrime][ji.ADrainPrimeEquDrainPrimeNodeOffset] +=
      ji.drainCond+ji.gds+ji.ggd;

    dFdx[ji.li_DrainPrime][ji.ADrainPrimeEquSourcePrimeNodeOffset] += -ji.gds-ji.gm;


    dFdx[ji.li_SourcePrime][ji.ASourcePrimeEquGateNodeOffset] -= ji.gm+ji.ggs;

    dFdx[ji.li_SourcePrime][ji.ASourcePrimeEquSourceNodeOffset] -= ji.sourceCond;

    dFdx[ji.li_SourcePrime][ji.ASourcePrimeEquDrainPrimeNodeOffset] -= ji.gds;

    dFdx[ji.li_SourcePrime][ji.ASourcePrimeEquSourcePrimeNodeOffset]
      += ji.sourceCond+ji.gds+ji.gm+ji.ggs;

    // dQdx matrix:

    dQdx[ji.li_Gate       ][ji.AGateEquGateNodeOffset        ] += ji.capgd+ji.capgs;

    dQdx[ji.li_Gate       ][ji.AGateEquDrainPrimeNodeOffset        ] -= ji.capgd;

    dQdx[ji.li_Gate       ][ji.AGateEquSourcePrimeNodeOffset       ] -= ji.capgs;

    dQdx[ji.li_DrainPrime ][ji.ADrainPrimeEquGateNodeOffset        ] -= ji.capgd;

    dQdx[ji.li_DrainPrime ][ji.ADrainPrimeEquDrainPrimeNodeOffset  ] += ji.capgd;

    dQdx[ji.li_SourcePrime][ji.ASourcePrimeEquGateNodeOffset       ] -= ji.capgs;

    dQdx[ji.li_SourcePrime][ji.ASourcePrimeEquSourcePrimeNodeOffset] += ji.capgs;
  }

  return true;
}
#endif

} // namespace MESFET
} // namespace Device
} // namespace Xyce

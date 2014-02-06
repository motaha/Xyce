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
// Filename       : $RCSfile: N_DEV_ISRC.C,v $
//
// Purpose        : Independent current source
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
// Revision Number: $Revision: 1.142.2.3 $
//
// Revision Date  : $Date: 2013/10/03 17:23:38 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------
#ifdef HAVE_CSTDIO
#include <cstdio>
#else
#include <stdio.h>
#endif

// ----------   Xyce Includes   ----------
#include <N_DEV_ISRC.h>
#include <N_DEV_SourceData.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_MatrixLoadData.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>

namespace Xyce {
namespace Device {

template<>
ParametricData<ISRC::Instance>::ParametricData()
{
    // Set up configuration constants:
    setNumNodes(2);
    setNumOptionalNodes(0);
    setNumFillNodes(0);
    setModelRequired(0);
    setPrimaryParameter("DCV0");

    // Set up double precision variables:
    // DC parameters
    addPar ("DCV0",       0.0, true, NO_DEP,
      &ISRC::Instance::DCV0,
      NULL, U_VOLT, CAT_NONE, "DC Current");

    // Pulse parameters
    addPar ("V0",         0.0, false, NO_DEP,
      &ISRC::Instance::par0,
      NULL, U_AMP, CAT_NONE, "Offset Current");

    addPar ("V1",         0.0, false, NO_DEP,
      &ISRC::Instance::par0,
      NULL, U_AMP, CAT_NONE, "Initial Current");

    addPar ("V2",         0.0, false, NO_DEP,
      &ISRC::Instance::par1,
      NULL, U_AMP, CAT_NONE, "Pulsed Current");

    addPar ("TD",         0.0, false, NO_DEP,
      &ISRC::Instance::par2,
      NULL, U_SECOND, CAT_NONE, "Delay");

    addPar ("TR",         0.0, false, NO_DEP,
      &ISRC::Instance::par3,
      NULL, U_SECOND, CAT_NONE, "Rise Time");

    addPar ("TF",         0.0, false, NO_DEP,
      &ISRC::Instance::par4,
      NULL, U_SECOND, CAT_NONE, "Fall Time");

    addPar ("PW",         0.0, false, NO_DEP,
      &ISRC::Instance::par5,
      NULL, U_SECOND, CAT_NONE, "Pulse Width");

    addPar ("PER",        0.0, false, NO_DEP,
      &ISRC::Instance::par6,
      NULL, U_SECOND, CAT_NONE, "Period");

    addPar ("SF",        0.0, false, NO_DEP,
      &ISRC::Instance::par7,
      NULL, U_NONE, CAT_NONE, "Scale Factor -- smooth pulse only");

    // Sin parameters
    addPar ("VA",         0.0, false, NO_DEP,
      &ISRC::Instance::par1,
      NULL, U_AMP, CAT_NONE, "Amplitude");

    addPar ("FREQ",       0.0, false, NO_DEP,
      &ISRC::Instance::par3,
      NULL, U_SECM1, CAT_NONE, "Frequency");

    addPar ("THETA",      0.0, false, NO_DEP,
      &ISRC::Instance::par4,
      NULL, U_NONE, CAT_NONE, "Theta");

    addPar ("PHASE",      0.0, false, NO_DEP,
      &ISRC::Instance::par5,
      NULL, U_NONE, CAT_NONE, "Phase");


    // Exp parameters
    addPar ("TD1",        0.0, false, NO_DEP,
      &ISRC::Instance::par2,
      NULL, U_SECOND, CAT_NONE, "Rise Delay Time");

    addPar ("TAU1",       0.0, false, NO_DEP,
      &ISRC::Instance::par3,
      NULL, U_SECOND, CAT_NONE, "Rise Time Constant");

    addPar ("TD2",        0.0, false, NO_DEP,
      &ISRC::Instance::par4,
      NULL, U_SECOND, CAT_NONE, "Fall Delay Time");

    addPar ("TAU2",       0.0, false, NO_DEP,
      &ISRC::Instance::par5,
      NULL, U_SECOND, CAT_NONE, "Fall Time Constant");

    // AC parameters
    addPar ("ACMAG",         0.0, false, NO_DEP,
      &ISRC::Instance::ACMAG,
      NULL, U_VOLT, CAT_NONE, "Amplitude");

    addPar ("ACPHASE",      0.0, false, NO_DEP,
      &ISRC::Instance::ACPHASE,
      NULL, U_NONE, CAT_NONE, "Phase");

    // SFFM parameters
    addPar ("FC",         0.0, false, NO_DEP,
      &ISRC::Instance::par2,
      NULL, U_SECM1, CAT_NONE, "Carrier Frequency");

    addPar ("FS",         0.0, false, NO_DEP,
      &ISRC::Instance::par4,
      NULL, U_SECM1, CAT_NONE, "Signal Frequency");

    addPar ("MDI",        0.0, false, NO_DEP,
      &ISRC::Instance::par3,
      NULL, U_NONE, CAT_NONE, "Modulation Index");


    // PWL params
    addPar ("R",          0.0, false, NO_DEP,
      &ISRC::Instance::REPEATTIME,
      NULL, U_SECOND, CAT_NONE, "Repeat Time");

    addPar ("REPEATTIME", 0.0, false, NO_DEP,
      &ISRC::Instance::REPEATTIME,
      NULL, U_SECOND, CAT_NONE, "Repeat Time");

    addPar ("T",          0.0, false, NO_DEP,
      &ISRC::Instance::T,
      NULL, U_SECOND, CAT_NONE, "Time");  // time-voltage pairs

    addPar ("V",          0.0, false, NO_DEP,
      &ISRC::Instance::V,
      NULL, U_AMP, CAT_NONE, "Current"); // time-voltage pairs

    // Set up non-double precision variables:
    addPar ("TRANSIENTSOURCETYPE", (int)_DC_DATA, false, NO_DEP,
      &ISRC::Instance::TRANSIENTSOURCETYPE,
      &ISRC::Instance::TRANSIENTSOURCETYPEgiven,
      U_NONE, CAT_NONE, "" );

    addPar ("ACSOURCETYPE", (int) _AC_DATA, false, NO_DEP,
      &ISRC::Instance::ACSOURCETYPE,
      &ISRC::Instance::ACSOURCETYPEgiven,
      U_NONE, CAT_NONE, "" );

    addPar ("DCSOURCETYPE", (int) _DC_DATA, false, NO_DEP,
      &ISRC::Instance::DCSOURCETYPE,
      &ISRC::Instance::DCSOURCETYPEgiven,
      U_NONE, CAT_NONE, "" );


    addPar ("NUM", 0, false, NO_DEP, &ISRC::Instance::NUM,
        NULL, U_NONE, CAT_NONE, "" );
    addPar ("REPEAT", 0, false, NO_DEP, &ISRC::Instance::REPEAT,
        NULL, U_NONE, CAT_NONE, "" );
}

template<>
ParametricData<ISRC::Model>::ParametricData()
{
}

namespace ISRC {
vector< vector<int> > Instance::jacStamp;

ParametricData<Instance> &Instance::getParametricData() {
  static ParametricData<Instance> parMap;

  return parMap;
}

ParametricData<Model> &Model::getParametricData() {
  static ParametricData<Model> parMap;

  return parMap;
}

//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : instance block constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/06/00
//-----------------------------------------------------------------------------
Instance::Instance(InstanceBlock & IB,
                   Model & model,
                   MatrixLoadData & mlData1,
                   SolverState &ss1,
                   ExternData  &ed1,
                   DeviceOptions & do1)
  : SourceInstance(IB, mlData1, ss1, ed1, do1),
    model_(model),
    li_Pos(-1),
    li_Neg(-1),
    li_store_dev_i(-1),
    DCV0(0.0),
    par0(0.0),
    par1(0.0),
    par2(0.0),
    par3(0.0),
    par4(0.0),
    par5(0.0),
    par6(0.0),
    par7(0.0),
    REPEATTIME(),
    T(0.0),
  V(0.0),
  ACMAG(1.0),
  ACPHASE(0.0),
  NUM(0),
  REPEAT(0),
  TRANSIENTSOURCETYPE(_DC_DATA),
  TRANSIENTSOURCETYPEgiven(false),
  ACSOURCETYPE(_AC_DATA),
  ACSOURCETYPEgiven(false),
  DCSOURCETYPE(_AC_DATA),
  DCSOURCETYPEgiven(false),
  gotParams(false)
{
  numIntVars = 0;
  numExtVars = 2;
  numStateVars = 0;
  numLeadCurrentStoreVars = 1; // lead current DEV_I

  if( jacStamp.empty() )
    jacStamp.resize(2);

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  // Set any non-constant parameter defaults:

  if (solState.ACspecified && ACSOURCETYPEgiven)
  {
    acData_ptr = new ACData (IB.params,solState,devOptions);
  }

  if (DCSOURCETYPEgiven) // this will always be given, if the source spec was valid.
  {
    dcData_ptr = new ConstData (IB.params,solState,devOptions);
  }

  if (solState.HBspecified || TRANSIENTSOURCETYPEgiven)
  {
    switch (TRANSIENTSOURCETYPE)
    {
      case _SIN_DATA:
        Data_ptr = new SinData (IB.params,solState,devOptions);
        break;

      case _EXP_DATA:
        Data_ptr = new ExpData (IB.params,solState,devOptions);
        break;

      case _PULSE_DATA:
        Data_ptr = new PulseData (IB.params,solState,devOptions);
        break;

      case _PWL_DATA:
        Data_ptr = new PWLinData (IB.params,solState,devOptions);
        break;

      case _SFFM_DATA:
        Data_ptr = new SFFMData (IB.params,solState,devOptions);
        break;

      case _DC_DATA:
        Data_ptr = new ConstData (IB.params,solState,devOptions);
        break;

      case _SMOOTH_PULSE_DATA:
        Data_ptr = new SmoothPulseData (IB.params,solState,devOptions);
        break;

//	    case _AC_DATA:
//	      Data_ptr = new ACData (IB.params,solState,devOptions);
//	      break;

      default:
        string msg = "Instance::Instance(IB)\n";
        msg += "\tCannot identify source data type for " + getName() + ".";
        std::ostringstream oss;
        oss << "Error in " << netlistLocation() << "\n" << msg;
        N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, oss.str());
        break;
    }
  }

  //defaultParamName = Data_ptr->defaultParamName_;
  defaultParamName = "DCV0";
  processParams();

  // Calculate any parameters specified as expressions:

  updateDependentParameters();
  processParams();

  // calculate dependent (ie computed) params and check for errors:


}

//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/06/00
//-----------------------------------------------------------------------------
bool Instance::processParams (string param)
{
  if (gotParams)
  {
    if (dcData_ptr != 0)
    {
      dcData_ptr->setParams (&DCV0);
    }
    if (acData_ptr != 0)
    {
      acData_ptr->setParams (&ACMAG);
    }
    if (Data_ptr != 0)
    {
      Data_ptr->setParams(&par0);
    }
  }
  else
  {
    if (dcData_ptr != 0)
    {
      dcData_ptr->getParams (&DCV0);
    }
    if (acData_ptr != 0)
    {
      acData_ptr->getParams (&ACMAG);
    }
    if (Data_ptr != 0)
    {
      Data_ptr->getParams(&par0);
    }
    gotParams = true;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/06/00
//-----------------------------------------------------------------------------
Instance::~Instance ()
{
  if (Data_ptr != 0)
  {
    delete Data_ptr;
  }
  if (acData_ptr != 0)
  {
    delete acData_ptr;
  }
  if (dcData_ptr != 0)
  {
    delete dcData_ptr;
  }
}

// additional Declarations
//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/20/02
//-----------------------------------------------------------------------------
void Instance::registerLIDs( const vector<int> & intLIDVecRef,
	                               const vector<int> & extLIDVecRef)
{
  string msg;

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline(
      "-----------------------------------------------------------------------------");
  if (devOptions.debugLevel > 0 )
  {
    cout << endl << dashedline << endl;
    cout << "  ISRCInstance::registerLIDs" << endl;
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

  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  // Now use these lists to obtain the indices into the
  // linear algebra entities.  This assumes an order.    Note that
  // for a current  source, there  will be no Jacobian entries.

  li_Pos = extLIDVec[0];
  li_Neg = extLIDVec[1];

#ifdef Xyce_DEBUG_DEVICE
  if (devOptions.debugLevel > 0 )
  {
    cout << "  li_Pos = " << li_Pos << endl;
    cout << "  li_Neg = " << li_Neg << endl;
    cout << dashedline << endl;
  }
#endif

}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/20/02
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs(const vector<int> & staLIDVecRef )
{
  string msg;

  // Check if the size of the ID lists corresponds to the
  // proper number of internal and external variables.
  int numSta = staLIDVecRef.size();

  if (numSta != numStateVars)
  {
    msg = "Instance::registerStateLIDs:";
    msg += "numSTa != numStateVars";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStoreLIDs
// Purpose       : One store var for device current.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 03/27/2013
//-----------------------------------------------------------------------------
void Instance::registerStoreLIDs(const vector<int> & stoLIDVecRef )
{
  string msg;

  // Check if the size of the ID lists corresponds to the
  // proper number of internal and external variables.
  int numSto = stoLIDVecRef.size();

  if (numSto != getNumStoreVars())
  {
    msg = "Instance::registerStoreLIDs:";
    msg += "numSto != numStoreVars";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
  }
  if( loadLeadCurrent )
  {
    li_store_dev_i = stoLIDVecRef[0];
  }
}


//-----------------------------------------------------------------------------
// Function      : Instance::getStoreNameMap
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 03/27/2013
//-----------------------------------------------------------------------------
map<int,string> & Instance::getStoreNameMap ()
{
  // set up the internal name map, if it hasn't been already.
  if( loadLeadCurrent && storeNameMap.empty ())
  {
    // change subcircuitname:devicetype_deviceName to
    // devicetype:subcircuitName:deviceName
    string modName(getName());
    spiceInternalName(modName);
    string tmpstr;
    tmpstr = modName+":DEV_I";
    storeNameMap[ li_store_dev_i ] = tmpstr;
  }

  return storeNameMap;
}


//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 9/5/2
//-----------------------------------------------------------------------------
const vector< vector<int> > & Instance::jacobianStamp() const
{
  return jacStamp;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadBVectorsforAC
//
// Purpose       : Loads the B-vector contributions for a single
//                 isrc instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Ting Mei, SNL
// Creation Date : 5/11
//-----------------------------------------------------------------------------
bool Instance::loadBVectorsforAC(double * bVecReal, double * bVecImag )
{
  if (acData_ptr != 0)
  {
    bool flag = true;
    acData_ptr->setRealFlag(flag);

    acData_ptr->updateSource ();
    double source = acData_ptr->returnSource();

    bVecReal[li_Pos] -= source;
    bVecReal[li_Neg] += source;

    flag = false;
    acData_ptr->setRealFlag(flag);

    acData_ptr->updateSource ();
    source = acData_ptr->returnSource();

    bVecImag[li_Pos] -= source;
    bVecImag[li_Neg] += source;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/29/01
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState ()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
// Purpose       : Loads the F-vector contributions for a single
//                 ISRC instance.
// Special Notes : Does nothing.
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  bool bsuccess = true;

  double * fVec = extData.daeFVectorRawPtr;

  // get the source value:
  SourceData *dataPtr = dcData_ptr; // by default assume the DC value.
  if ( (solState.HBspecified || solState.tranopFlag || solState.transientFlag) && Data_ptr != 0 )
  {
    dataPtr = Data_ptr;
  }

  double source = 0.0;
  if (dataPtr != 0)
  {
    source = dataPtr->returnSource();
  }
  fVec[li_Pos] += source;
  fVec[li_Neg] -= source;

  if( loadLeadCurrent )
  {
    double * stoVec = extData.nextStoVectorRawPtr;
    stoVec[li_store_dev_i] = source;
  }

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : Model::Model
// Purpose       : copy constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/06/00
//-----------------------------------------------------------------------------
Model::Model (const ModelBlock & MB,
                                        SolverState & ss1,
                                        DeviceOptions & do1)
  : DeviceModel(MB,ss1,do1)
{

}

//-----------------------------------------------------------------------------
// Function      : Model::~Model
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/06/00
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

//----------------------------------------------------------------------------
// Function      : Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Lon Waters
// Creation Date : 07/29/2002
//----------------------------------------------------------------------------
bool Model::processParams(string param)
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

//-----------------------------------------------------------------------------
// Function      : Model::printOutInstances
// Purpose       : debugging tool.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/03/00
//-----------------------------------------------------------------------------
std::ostream &Model::printOutInstances(std::ostream &os) const
{
#ifdef Xyce_DEBUG_DEVICE

  vector<Instance*>::const_iterator iter;
  vector<Instance*>::const_iterator first = instanceContainer.begin();
  vector<Instance*>::const_iterator last  = instanceContainer.end();

  int i;
  os << endl;
  os << "    name     modelName  Parameters" << endl;
  for (i=0, iter=first; iter!=last; ++iter, ++i)
  {
    os << "  " << i << ": " << (*iter)->getName() << "      ";
    os << (*iter)->getModelName();
    os << endl;
    if ( (*iter)->Data_ptr != 0 )
    {
      (*iter)->Data_ptr->printOutParams ();
    }

    if ( (*iter)->dcData_ptr != 0 )
    {
      (*iter)->dcData_ptr->printOutParams ();
    }

    if ( (*iter)->acData_ptr != 0 )
    {
      (*iter)->acData_ptr->printOutParams ();
    }
  }

  os << endl;
#endif
  return os;
}



//-----------------------------------------------------------------------------
// ISRC Master functions:
//-----------------------------------------------------------------------------

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
  return true;
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
    Instance & inst = *(*it);

    SourceData *dataPtr = inst.dcData_ptr; // by default assume the DC value.
    if ( (getSolverState().HBspecified || getSolverState().tranopFlag || getSolverState().transientFlag) && inst.Data_ptr != 0 )
    {
      dataPtr = inst.Data_ptr;
    }

    double source = 0.0;
    if (dataPtr != 0)
    {
      source = dataPtr->returnSource();
    }
    fVec[inst.li_Pos] += source;
    fVec[inst.li_Neg] -= source;
    if( inst.loadLeadCurrent )
    {
      storeLeadF[inst.li_store_dev_i] = source;
    }

  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Master::loadDAEMatrices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
bool Master::loadDAEMatrices (N_LAS_Matrix & dFdx, N_LAS_Matrix & dQdx)
{
  return true;
}

} // namespace Resistor
} // namespace Device
} // namespace Xyce

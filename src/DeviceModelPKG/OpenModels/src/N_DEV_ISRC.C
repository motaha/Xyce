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
// Revision Number: $Revision: 1.169.2.2 $
//
// Revision Date  : $Date: 2014/03/06 23:33:43 $
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
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_ISRC.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_SourceData.h>
#include <N_DEV_Message.h>
#include <N_ERH_ErrorMgr.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>

namespace Xyce {
namespace Device {
namespace ISRC {

void Traits::loadInstanceParameters(ParametricData<ISRC::Instance> &p)
{
  // DC parameters
  p.addPar ("DCV0", 0.0, &ISRC::Instance::DCV0)
    .setOriginalValueStored(true)
    .setUnit(U_VOLT)
    .setDescription("DC Current");

  // Pulse parameters
  p.addPar ("V0", 0.0, &ISRC::Instance::par0)
    .setUnit(U_AMP)
    .setDescription("Offset Current");

  p.addPar ("V1", 0.0, &ISRC::Instance::par0)
    .setUnit(U_AMP)
    .setDescription("Initial Current");

  p.addPar ("V2", 0.0, &ISRC::Instance::par1)
    .setUnit(U_AMP)
    .setDescription("Pulsed Current");

  p.addPar ("TD", 0.0, &ISRC::Instance::par2)
    .setUnit(U_SECOND)
    .setDescription("Delay");

  p.addPar ("TR", 0.0, &ISRC::Instance::par3)
    .setUnit(U_SECOND)
    .setDescription("Rise Time");

  p.addPar ("TF", 0.0, &ISRC::Instance::par4)
    .setUnit(U_SECOND)
    .setDescription("Fall Time");

  p.addPar ("PW", 0.0, &ISRC::Instance::par5)
    .setUnit(U_SECOND)
    .setDescription("Pulse Width");

  p.addPar ("PER", 0.0, &ISRC::Instance::par6)
    .setUnit(U_SECOND)
    .setDescription("Period");

  p.addPar ("SF", 0.0, &ISRC::Instance::par7)
    .setDescription("Scale Factor -- smooth pulse only");

  // Sin parameters
  p.addPar ("VA", 0.0, &ISRC::Instance::par1)
    .setUnit(U_AMP)
    .setDescription("Amplitude");

  p.addPar ("FREQ", 0.0, &ISRC::Instance::par3)
    .setUnit(U_SECM1)
    .setDescription("Frequency");

  p.addPar ("THETA", 0.0, &ISRC::Instance::par4)
    .setDescription("Theta");

  p.addPar ("PHASE", 0.0, &ISRC::Instance::par5)
    .setDescription("Phase");

  // Exp parameters
  p.addPar ("TD1", 0.0, &ISRC::Instance::par2)
    .setUnit(U_SECOND)
    .setDescription("Rise Delay Time");

  p.addPar ("TAU1", 0.0, &ISRC::Instance::par3)
    .setUnit(U_SECOND)
    .setDescription("Rise Time Constant");

  p.addPar ("TD2", 0.0, &ISRC::Instance::par4)
    .setUnit(U_SECOND)
    .setDescription("Fall Delay Time");

  p.addPar ("TAU2", 0.0, &ISRC::Instance::par5)
    .setUnit(U_SECOND)
    .setDescription("Fall Time Constant");

  // AC parameters
  p.addPar ("ACMAG", 0.0, &ISRC::Instance::ACMAG)
    .setUnit(U_VOLT)
    .setDescription("Amplitude");

  p.addPar ("ACPHASE", 0.0, &ISRC::Instance::ACPHASE)
    .setDescription("Phase");

  // SFFM parameters
  p.addPar ("FC", 0.0, &ISRC::Instance::par2)
    .setUnit(U_SECM1)
    .setDescription("Carrier Frequency");

  p.addPar ("FS", 0.0, &ISRC::Instance::par4)
    .setUnit(U_SECM1)
    .setDescription("Signal Frequency");

  p.addPar ("MDI", 0.0, &ISRC::Instance::par3)
    .setDescription("Modulation Index");

  // PWL params
  p.addPar ("R", 0.0, &ISRC::Instance::REPEATTIME)
    .setUnit(U_SECOND)
    .setDescription("Repeat Time");

  p.addPar ("REPEATTIME", 0.0, &ISRC::Instance::REPEATTIME)
    .setUnit(U_SECOND)
    .setDescription("Repeat Time");

  p.addPar ("T", 0.0, &ISRC::Instance::T)
    .setUnit(U_SECOND)
    .setDescription("Time");  // time-voltage pairs

  p.addPar ("V", 0.0, &ISRC::Instance::V)
    .setUnit(U_AMP)
    .setDescription("Current"); // time-voltage pairs

  // Set up non-double precision variables:
  p.addPar ("TRANSIENTSOURCETYPE", (int)_DC_DATA, &ISRC::Instance::TRANSIENTSOURCETYPE)
    .setGivenMember(&ISRC::Instance::TRANSIENTSOURCETYPEgiven);

  p.addPar ("ACSOURCETYPE", (int) _AC_DATA, &ISRC::Instance::ACSOURCETYPE)
    .setGivenMember(&ISRC::Instance::ACSOURCETYPEgiven);

  p.addPar ("DCSOURCETYPE", (int) _DC_DATA, &ISRC::Instance::DCSOURCETYPE)
    .setGivenMember(&ISRC::Instance::DCSOURCETYPEgiven);

  p.addPar ("NUM", 0, &ISRC::Instance::NUM);

  p.addPar ("REPEAT", false, &ISRC::Instance::REPEAT);
}

void Traits::loadModelParameters(ParametricData<ISRC::Model> &p)
{
}


std::vector< std::vector<int> > Instance::jacStamp;

//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : instance block constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/06/00
//-----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock & instance_block,
  Model &               model,
  const FactoryBlock &  factory_block)
  : SourceInstance(instance_block, configuration.getInstanceParameters(), factory_block),
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
  REPEAT(false),
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
  setParams (instance_block.params);

  // Set any non-constant parameter defaults:

  if (getSolverState().ACspecified && ACSOURCETYPEgiven)
  {
    acData_ptr = new ACData (instance_block.params,getSolverState(),getDeviceOptions());
  }

  if (DCSOURCETYPEgiven) // this will always be given, if the source spec was valid.
  {
    dcData_ptr = new ConstData (instance_block.params,getSolverState(),getDeviceOptions());
  }

  if (getSolverState().HBspecified || TRANSIENTSOURCETYPEgiven)
  {
    switch (TRANSIENTSOURCETYPE)
    {
      case _SIN_DATA:
        Data_ptr = new SinData (instance_block.params,getSolverState(),getDeviceOptions());
        break;

      case _EXP_DATA:
        Data_ptr = new ExpData (instance_block.params,getSolverState(),getDeviceOptions());
        break;

      case _PULSE_DATA:
        Data_ptr = new PulseData (instance_block.params,getSolverState(),getDeviceOptions());
        break;

      case _PWL_DATA:
        Data_ptr = new PWLinData (instance_block.params,getSolverState(),getDeviceOptions());
        break;

      case _SFFM_DATA:
        Data_ptr = new SFFMData (instance_block.params,getSolverState(),getDeviceOptions());
        break;

      case _DC_DATA:
        Data_ptr = new ConstData (instance_block.params,getSolverState(),getDeviceOptions());
        break;

      case _SMOOTH_PULSE_DATA:
        Data_ptr = new SmoothPulseData (instance_block.params,getSolverState(),getDeviceOptions());
        break;

//	    case _AC_DATA:
//	      Data_ptr = new ACData (instance_block.params,getSolverState(),getDeviceOptions());
//	      break;

      default:
        UserError0(*this) << "Cannot identify source data type for " << getName();
        break;
    }
  }

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
bool Instance::processParams ()
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
void Instance::registerLIDs( const std::vector<int> & intLIDVecRef,
	                               const std::vector<int> & extLIDVecRef)
{
  AssertLIDs(intLIDVecRef.size() == numIntVars);
  AssertLIDs(extLIDVecRef.size() == numExtVars);

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 )
  {
    Xyce::dout() << std::endl << section_divider << std::endl;
    Xyce::dout() << "  ISRCInstance::registerLIDs" << std::endl;
    Xyce::dout() << "  name = " << getName() << std::endl;
  }
#endif

  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  // Now use these lists to obtain the indices into the
  // linear algebra entities.  This assumes an order.    Note that
  // for a current  source, there  will be no Jacobian entries.

  li_Pos = extLIDVec[0];
  li_Neg = extLIDVec[1];

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 )
  {
    Xyce::dout() << "  li_Pos = " << li_Pos << std::endl;
    Xyce::dout() << "  li_Neg = " << li_Neg << std::endl;
    Xyce::dout() << section_divider << std::endl;
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
void Instance::registerStateLIDs(const std::vector<int> & staLIDVecRef )
{
  AssertLIDs(staLIDVecRef.size() == numStateVars);
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStoreLIDs
// Purpose       : One store var for device current.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 03/27/2013
//-----------------------------------------------------------------------------
void Instance::registerStoreLIDs(const std::vector<int> & stoLIDVecRef )
{
  AssertLIDs(stoLIDVecRef.size() == getNumStoreVars());

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
std::map<int,std::string> & Instance::getStoreNameMap ()
{
  // set up the internal name map, if it hasn't been already.
  if( loadLeadCurrent && storeNameMap.empty ())
  {
    // change subcircuitname:devicetype_deviceName to
    // devicetype:subcircuitName:deviceName
    std::string modName(getName());
    spiceInternalName(modName);
    std::string tmpstr;
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
const std::vector< std::vector<int> > & Instance::jacobianStamp() const
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
  if ( (getSolverState().HBspecified || getSolverState().tranopFlag || getSolverState().transientFlag) && Data_ptr != 0 )
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
Model::Model(
  const Configuration & configuration,
  const ModelBlock &    MB,
  const FactoryBlock &  factory_block)
  : DeviceModel(MB, configuration.getModelParameters(), factory_block)
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
  std::vector<Instance*>::iterator iter;
  std::vector<Instance*>::iterator first = instanceContainer.begin();
  std::vector<Instance*>::iterator last  = instanceContainer.end();

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
bool Model::processParams()
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

  std::vector<Instance*>::const_iterator iter;
  std::vector<Instance*>::const_iterator first = instanceContainer.begin();
  std::vector<Instance*>::const_iterator last  = instanceContainer.end();

  int i;
  os << std::endl;
  os << "    name     modelName  Parameters" << std::endl;
  for (i=0, iter=first; iter!=last; ++iter, ++i)
  {
    os << "  " << i << ": " << (*iter)->getName() << "      ";
    os << getName();
    os << std::endl;
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

  os << std::endl;
#endif
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
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
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

Device *Traits::factory(const Configuration &configuration, const FactoryBlock &factory_block)
{

  return new Master(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);
}

void registerDevice()
{
  Config<Traits>::addConfiguration()
    .registerDevice("i", 1);
}

} // namespace Resistor
} // namespace Device
} // namespace Xyce

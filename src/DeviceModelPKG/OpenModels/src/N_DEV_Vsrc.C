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
// Filename       : $RCSfile: N_DEV_Vsrc.C,v $
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
// Revision Number: $Revision: 1.226.2.2 $
//
// Revision Date  : $Date: 2014/03/06 23:33:43 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------
#include <Xyce_config.h>

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_SourceData.h>
#include <N_DEV_Vsrc.h>
#include <N_DEV_Message.h>
#include <N_ERH_ErrorMgr.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>
#include <N_UTL_BreakPoint.h>

namespace Xyce {
namespace Device {


namespace Vsrc {


void Traits::loadInstanceParameters(ParametricData<Vsrc::Instance> &p)
{
    // Set up double precision variables:
    // DC parameters
    p.addPar ("DCV0",         0.0, true, ParameterType::NO_DEP,
      &Vsrc::Instance::DCV0,
      NULL, U_VOLT, CAT_NONE, "DC Voltage");

    // Pulse parameters
    p.addPar ("V0",         0.0, false, ParameterType::NO_DEP,
      &Vsrc::Instance::par0,
      NULL, U_VOLT, CAT_NONE, "Offset Voltage");

    p.addPar ("V1",         0.0, false, ParameterType::NO_DEP,
      &Vsrc::Instance::par0,
      NULL, U_VOLT, CAT_NONE, "Initial Voltage");

    p.addPar ("V2",         0.0, false, ParameterType::NO_DEP,
      &Vsrc::Instance::par1,
      NULL, U_VOLT, CAT_NONE, "Pulsed Voltage");

    p.addPar ("TD",         0.0, false, ParameterType::NO_DEP,
      &Vsrc::Instance::par2,
      NULL, U_SECOND, CAT_NONE, "Delay");

    p.addPar ("TR",         0.0, false, ParameterType::NO_DEP,
      &Vsrc::Instance::par3,
      NULL, U_SECOND, CAT_NONE, "Rise Time");

    p.addPar ("TF",         0.0, false, ParameterType::NO_DEP,
      &Vsrc::Instance::par4,
      NULL, U_SECOND, CAT_NONE, "Fall Time");

    p.addPar ("PW",         0.0, false, ParameterType::NO_DEP,
      &Vsrc::Instance::par5,
      NULL, U_SECOND, CAT_NONE, "Pulse Width");

    p.addPar ("PER",        0.0, false, ParameterType::NO_DEP,
      &Vsrc::Instance::par6,
      NULL, U_SECOND, CAT_NONE, "Period");

    p.addPar ("SF",        0.0, false, ParameterType::NO_DEP,
      &Vsrc::Instance::par7,
      NULL, U_NONE, CAT_NONE, "Scale Factor -- smooth pulse only");

    // Sin parameters
    p.addPar ("VA",         0.0, false, ParameterType::NO_DEP,
      &Vsrc::Instance::par1,
      NULL, U_VOLT, CAT_NONE, "Amplitude");

    p.addPar ("FREQ",       0.0, false, ParameterType::NO_DEP,
      &Vsrc::Instance::par3,
      NULL, U_SECM1, CAT_NONE, "Frequency");

    p.addPar ("THETA",      0.0, false, ParameterType::NO_DEP,
      &Vsrc::Instance::par4,
      NULL, U_NONE, CAT_NONE, "Theta");

    p.addPar ("PHASE",      0.0, false, ParameterType::NO_DEP,
      &Vsrc::Instance::par5,
      NULL, U_NONE, CAT_NONE, "Phase");

    // Exp parameters
    p.addPar ("TD1",        0.0, false, ParameterType::NO_DEP,
      &Vsrc::Instance::par2,
      NULL, U_SECOND, CAT_NONE, "Rise Delay Time");

    p.addPar ("TAU1",       0.0, false, ParameterType::NO_DEP,
      &Vsrc::Instance::par3,
      NULL, U_SECOND, CAT_NONE, "Rise Time Constant");

    p.addPar ("TD2",        0.0, false, ParameterType::NO_DEP,
      &Vsrc::Instance::par4,
      NULL, U_SECOND, CAT_NONE, "Fall Delay Time");

    p.addPar ("TAU2",       0.0, false, ParameterType::NO_DEP,
      &Vsrc::Instance::par5,
      NULL, U_SECOND, CAT_NONE, "Fall Time Constant");

// AC parameters
    p.addPar ("ACMAG",         0.0, false, ParameterType::NO_DEP,
      &Vsrc::Instance::ACMAG,
      NULL, U_VOLT, CAT_NONE, "Amplitude");

    p.addPar ("ACPHASE",      0.0, false, ParameterType::NO_DEP,
      &Vsrc::Instance::ACPHASE,
      NULL, U_NONE, CAT_NONE, "Phase");

    // SFFM parameters
    p.addPar ("FC",         0.0, false, ParameterType::NO_DEP,
      &Vsrc::Instance::par2,
      NULL, U_SECM1, CAT_NONE, "Carrier Frequency");

    p.addPar ("FS",         0.0, false, ParameterType::NO_DEP,
      &Vsrc::Instance::par4,
      NULL, U_SECM1, CAT_NONE, "Signal Frequency");

    p.addPar ("MDI",        0.0, false, ParameterType::NO_DEP,
      &Vsrc::Instance::par3,
      NULL, U_NONE, CAT_NONE, "Modulation Index");

    // PWL params
    p.addPar ("R",          0.0, false, ParameterType::NO_DEP,
      &Vsrc::Instance::REPEATTIME,
      NULL, U_SECOND, CAT_NONE, "Repeat Time");

    p.addPar ("REPEATTIME", 0.0, false, ParameterType::NO_DEP,
      &Vsrc::Instance::REPEATTIME,
      NULL, U_SECOND, CAT_NONE, "Repeat Time");

    p.addPar ("T",          0.0, false, ParameterType::NO_DEP,
      &Vsrc::Instance::T,
      NULL, U_SECOND, CAT_NONE, "Time");  // time-voltage pairs

    p.addPar ("V",          0.0, false, ParameterType::NO_DEP,
      &Vsrc::Instance::V,
      NULL, U_VOLT, CAT_NONE, "Voltage"); // time-voltage pairs

    // Set up exceptions (ie variables that are not doubles):
    p.addPar ("TRANSIENTSOURCETYPE", (int) _DC_DATA, false, ParameterType::NO_DEP,
      &Vsrc::Instance::TRANSIENTSOURCETYPE,
      &Vsrc::Instance::TRANSIENTSOURCETYPEgiven,
      U_NONE, CAT_NONE, "" );

    p.addPar ("ACSOURCETYPE", (int) _AC_DATA, false, ParameterType::NO_DEP,
      &Vsrc::Instance::ACSOURCETYPE,
      &Vsrc::Instance::ACSOURCETYPEgiven,
      U_NONE, CAT_NONE, "" );

    p.addPar ("DCSOURCETYPE", (int) _DC_DATA, false, ParameterType::NO_DEP,
      &Vsrc::Instance::DCSOURCETYPE,
      &Vsrc::Instance::DCSOURCETYPEgiven,
      U_NONE, CAT_NONE, "" );

    p.addPar ("NUM",        0, false, ParameterType::NO_DEP,
      &Vsrc::Instance::NUM,
      NULL, U_NONE, CAT_NONE, "" );

    p.addPar ("REPEAT", false, false, ParameterType::NO_DEP,
      &Vsrc::Instance::REPEAT,
      NULL, U_NONE, CAT_NONE, "" );
}

void Traits::loadModelParameters(ParametricData<Vsrc::Model> &p)
{
}


std::vector< std::vector<int> > Instance::jacStamp;
std::vector< std::vector<int> > Instance::jacStampPDE;

// Class Instance
//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : "instance block" constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 04/06/00
//-----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock &         IB,
  Model &                       Viter,
  const FactoryBlock &          factory_block)
  : SourceInstance(IB, configuration.getInstanceParameters(), factory_block),
    model_(Viter),
  srcCurrent(0.0),
  srcVoltage(0.0),
  srcDrop(0.0),
  srcBC(0.0),

  scale(1.0),
  nlstep(-1),
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

  li_Pos(-1),
  li_Neg(-1),
  li_Bra(-1),

  gotParams(false),

  ABraEquPosNodeOffset(-1),
  ABraEquNegNodeOffset(-1),
  APosEquBraVarOffset(-1),
  ANegEquBraVarOffset(-1),
  ABraEquBraVarOffset(-1),
  APosEquPosNodeOffset(-1),
  ANegEquNegNodeOffset(-1),

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  fBraEquPosNodePtr(0),
  fBraEquNegNodePtr(0),
  fPosEquBraVarPtr(0),
  fNegEquBraVarPtr(0),
  fPosEquPosNodePtr(0),
  fNegEquNegNodePtr(0),
  fBraEquBraVarPtr(0),
#endif

  source(0.0),
  v_pos(0.0),
  v_neg(0.0),
  i_bra(0.0)
{
  numIntVars   = 1;
  numExtVars   = 2;
  numStateVars = 0;

  if( jacStamp.empty() )
  {
    jacStamp.resize(3);
    jacStamp[0].resize(1);
    jacStamp[0][0] = 2;
    jacStamp[1].resize(1);
    jacStamp[1][0] = 2;
    jacStamp[2].resize(2);
    jacStamp[2][0] = 0;
    jacStamp[2][1] = 1;

    // PDE supporting stamp.  This includes diagonal elements, needed by the
    // 2-level Newton.
    jacStampPDE.resize(3);
    jacStampPDE[0].resize(2);
    jacStampPDE[0][0] = 0;
    jacStampPDE[0][1] = 2;
    jacStampPDE[1].resize(2);
    jacStampPDE[1][0] = 1;
    jacStampPDE[1][1] = 2;
    jacStampPDE[2].resize(3);
    jacStampPDE[2][0] = 0;
    jacStampPDE[2][1] = 1;
    jacStampPDE[2][2] = 2;
  }

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  // Set any non-constant parameter defaults:
  if (getSolverState().ACspecified && ACSOURCETYPEgiven)
  {
    acData_ptr = new ACData (IB.params,getSolverState(),getDeviceOptions());
  }

  if (DCSOURCETYPEgiven) // this will always be given, if the source spec was valid.
  {
    dcData_ptr = new ConstData (IB.params,getSolverState(),getDeviceOptions());
  }

  if (getSolverState().HBspecified || TRANSIENTSOURCETYPEgiven)
  {
    switch (TRANSIENTSOURCETYPE)
    {
      case _SIN_DATA:
        Data_ptr = new SinData (IB.params,getSolverState(),getDeviceOptions());
        break;

      case _EXP_DATA:
        Data_ptr = new ExpData (IB.params,getSolverState(),getDeviceOptions());
        break;

      case _PULSE_DATA:
        Data_ptr = new PulseData (IB.params,getSolverState(),getDeviceOptions());
        break;

      case _PWL_DATA:
        Data_ptr = new PWLinData (IB.params,getSolverState(),getDeviceOptions());
        break;

      case _SFFM_DATA:
        Data_ptr = new SFFMData (IB.params,getSolverState(),getDeviceOptions());
        break;

      case _DC_DATA:
        Data_ptr = new ConstData (IB.params,getSolverState(),getDeviceOptions());
        break;

      case _SMOOTH_PULSE_DATA:
        Data_ptr = new SmoothPulseData (IB.params,getSolverState(),getDeviceOptions());
        break;

//      case _AC_DATA:
//        Data_ptr = new ACData (IB.params,getSolverState(),getDeviceOptions());
 //        break;

      default:
        UserFatal0(*this) << "Cannot identify source data type for " << getName();
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
// Creation Date : 07/26/03
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

// Additional Declarations

//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/20/02
//-----------------------------------------------------------------------------
void Instance::registerLIDs ( const std::vector<int> & intLIDVecRef,
	                                const std::vector<int> & extLIDVecRef)
{
  std::string msg;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 )
  {
    Xyce::dout() << std::endl << section_divider << std::endl;
    Xyce::dout() << "  VsrcInstance::registerLIDs" << std::endl;
    Xyce::dout() << "  name = " << getName() << std::endl;
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
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL_0,msg);
  }

  if (numExt != numExtVars)
  {
    msg = "Instance::registerLIDs:";
    msg += "numExt != numExtVars";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL_0,msg);
  }

  // copy over the global ID lists.
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  // Now use these lists to obtain the indices into the
  // linear algebra entities.  This assumes an order.
  // For the matrix  indices, first do the rows.

  li_Pos = extLIDVec[0];
  li_Neg = extLIDVec[1];
  li_Bra = intLIDVec[0];

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 )
  {
    Xyce::dout() << "  li_Pos = " << li_Pos << std::endl;
    Xyce::dout() << "  li_Neg = " << li_Neg << std::endl;
    Xyce::dout() << "  li_Bra = " << li_Bra << std::endl;
    Xyce::dout() << section_divider << std::endl;
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
std::map<int,std::string> & Instance::getIntNameMap ()
{
  // set up the internal name map, if it hasn't been already.
  if (intNameMap.empty ())
  {
    // set up the internal name map
    std::string tmpstr(getName()+"_branch");
    spiceInternalName (tmpstr);
    intNameMap[ li_Bra ] = tmpstr;
  }

  return intNameMap;
}



//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/20/02
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs( const std::vector<int> & staLIDVecRef )
{
  AssertLIDs(staLIDVecRef.size() == numStateVars);
}

//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/21/02
//-----------------------------------------------------------------------------
const std::vector< std::vector<int> > & Instance::jacobianStamp() const
{
  if (getSolverState().PDESystemFlag)
    return jacStampPDE;

  return jacStamp;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerJacLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/27/02
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );
  if (getSolverState().PDESystemFlag)
  {
    APosEquBraVarOffset  = jacLIDVec[0][1];
    ANegEquBraVarOffset  = jacLIDVec[1][1];
    ABraEquPosNodeOffset = jacLIDVec[2][0];
    ABraEquNegNodeOffset = jacLIDVec[2][1];
  }
  else
  {
    APosEquBraVarOffset  = jacLIDVec[0][0];
    ANegEquBraVarOffset  = jacLIDVec[1][0];
    ABraEquPosNodeOffset = jacLIDVec[2][0];
    ABraEquNegNodeOffset = jacLIDVec[2][1];
  }

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

  fPosEquBraVarPtr = &(dFdx[li_Pos][APosEquBraVarOffset]);
  fNegEquBraVarPtr = &(dFdx[li_Neg][ANegEquBraVarOffset]);
  fBraEquPosNodePtr = &(dFdx[li_Bra][ABraEquPosNodeOffset]);
  fBraEquNegNodePtr = &(dFdx[li_Bra][ABraEquNegNodeOffset]);
#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, 9233, Computational Sciences
// Creation Date : 03/05/04
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars ()
{
  double * solVec = extData.nextSolVectorRawPtr;

  // Get the value for the source.
  SourceData *dataPtr = dcData_ptr; // by default assume the DC value.
  if ( (getSolverState().HBspecified || getSolverState().tranopFlag || getSolverState().transientFlag) && Data_ptr != 0 )
  {
    dataPtr = Data_ptr;
  }

  if (dataPtr != 0)
  {
    source = dataPtr->returnSource();
  }
  else
  {
    source = 0.0;
  }

  // get the value for v_pos, v_neg, i_bra
  v_pos = solVec[li_Pos];
  v_neg = solVec[li_Neg];
  i_bra = solVec[li_Bra];

  srcCurrent = i_bra;
  srcDrop    = (v_pos-v_neg);
  srcBC      = source;
  srcVoltage = srcDrop-srcBC;

  if( getDeviceOptions().scale_src != 0.0 )
  {
    srcCurrent *= scale;

    // first newton step, generate new scaling
    //if( nlsMgrPtr->getNonLinearIter() != nlstep )
    if( getSolverState().newtonIter != nlstep )
    {
      //nlstep = nlsMgrPtr->getNonLinearIter();
      nlstep = getSolverState().newtonIter;
      double new_scale = fabs(i_bra) * scale;

      scale = Xycemax( new_scale, getDeviceOptions().scale_src );
    }

    srcVoltage *= scale;
    srcDrop    *= scale;
    srcBC      *= scale;
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
  bool bsuccess = updateIntermediateVars ();
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 vsrc instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  double * fVec = extData.daeFVectorRawPtr;

  fVec[li_Pos] += srcCurrent;
  fVec[li_Neg] += -srcCurrent;
  fVec[li_Bra] += srcDrop;
  fVec[li_Bra] -= srcBC;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadBVectorsforAC
//
// Purpose       : Loads the B-vector contributions for a single
//                 vsrc instance.
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
    source = acData_ptr->returnSource();
    srcBC = source;

    if( getDeviceOptions().scale_src != 0.0 )
    {
      srcBC *= scale;
    }

    bVecReal[li_Bra] += srcBC;

    flag = false;
    acData_ptr->setRealFlag(flag);

    acData_ptr->updateSource ();
    source = acData_ptr->returnSource();
    srcBC = source;

    if( getDeviceOptions().scale_src != 0.0 )
    {
      srcBC *= scale;
    }

    bVecImag[li_Bra] += srcBC;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 resistor  instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/04
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  N_LAS_Matrix & dFdx = *(extData.dFdxMatrixPtr);

  dFdx[li_Pos][APosEquBraVarOffset] += 1.0;
  dFdx[li_Neg][ANegEquBraVarOffset] -= 1.0;
  dFdx[li_Bra][ABraEquPosNodeOffset] += 1.0;
  dFdx[li_Bra][ABraEquNegNodeOffset] -= 1.0;

  return true;
}

// end of new-DAE functions

//-----------------------------------------------------------------------------
// Function      : Instance::getMaxTimeStepSize
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/23/03
//-----------------------------------------------------------------------------
double Instance::getMaxTimeStepSize  ()
{
  double maxStep = 1.0e+100;
  if (Data_ptr != 0)
  {
    maxStep = Data_ptr->getMaxTimeStepSize ();
  }
  return maxStep;
}

//-----------------------------------------------------------------------------
// Function      : Instance::varTypes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/17/04
//-----------------------------------------------------------------------------
void Instance::varTypes( std::vector<char> & varTypeVec )
{
  varTypeVec.resize(1);
  varTypeVec[0] = 'I';
}

// Class Model

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
  : DeviceModel(MB, configuration.getModelParameters(), factory_block),
    DC_TRAN (0)
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

// Additional Declarations

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
  std::vector<Instance*>::const_iterator iter;
  std::vector<Instance*>::const_iterator first = instanceContainer.begin();
  std::vector<Instance*>::const_iterator last  = instanceContainer.end();

  int i;
  os << std::endl;
  os << "    name     model name  Parameters" << std::endl;
  for (i=0, iter=first; iter!=last; ++iter, ++i)
  {
    os << "  " << i << ": " << (*iter)->getName() << "      ";
    os << getName();
    os << std::endl;
  }

  os << std::endl;

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
// Vsrc Master functions:
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
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & vi        = *(*it);
    // Get the value for the source.
    SourceData *dataPtr  = vi.dcData_ptr; // by default assume the DC value.
    if ( (getSolverState().HBspecified || getSolverState().tranopFlag || getSolverState().transientFlag) && vi.Data_ptr != 0 )
    {
      dataPtr            = vi.Data_ptr;
    }

    if (dataPtr != 0)
    {
      vi.source          = dataPtr->returnSource();
    }
    else
    {
      vi.source          = 0.0;
    }

    // get the value for v_pos, v_neg, i_bra
    vi.v_pos             = solVec[vi.li_Pos];
    vi.v_neg             = solVec[vi.li_Neg];
    vi.i_bra             = solVec[vi.li_Bra];

    vi.srcCurrent        = vi.i_bra;
    vi.srcDrop           = (vi.v_pos-vi.v_neg);
    vi.srcBC             = vi.source;
    vi.srcVoltage        = vi.srcDrop-vi.srcBC;

    if( getDeviceOptions().scale_src != 0.0 )
    {
      vi.srcCurrent     *= vi.scale;

      // first newton step, generate new scaling
      if( getSolverState().newtonIter != vi.nlstep )
      {
        vi.nlstep        = getSolverState().newtonIter;
        double new_scale = fabs(vi.i_bra) * vi.scale;

        vi.scale         = Xycemax( new_scale, getDeviceOptions().scale_src );
      }

      vi.srcVoltage     *= vi.scale;
      vi.srcDrop        *= vi.scale;
      vi.srcBC          *= vi.scale;
    }
  }

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
    Instance & vi = *(*it);
    fVec[vi.li_Pos] += vi.srcCurrent;
    fVec[vi.li_Neg] += -vi.srcCurrent;
    fVec[vi.li_Bra] += vi.srcDrop;
    fVec[vi.li_Bra] -= vi.srcBC;
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
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & vi = *(*it);

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
    *(vi.fPosEquBraVarPtr) += 1.0;
    *(vi.fNegEquBraVarPtr) -= 1.0;
    *(vi.fBraEquPosNodePtr) += 1.0;
    *(vi.fBraEquNegNodePtr) -= 1.0;
#else
    dFdx[vi.li_Pos][vi.APosEquBraVarOffset] += 1.0;
    dFdx[vi.li_Neg][vi.ANegEquBraVarOffset] -= 1.0;
    dFdx[vi.li_Bra][vi.ABraEquPosNodeOffset] += 1.0;
    dFdx[vi.li_Bra][vi.ABraEquNegNodeOffset] -= 1.0;
#endif
  }
  return true;
}

Device *
Traits::factory(const Configuration &configuration, const FactoryBlock &factory_block)
{

  return new Master(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);
}

void registerDevice()
{
  Config<Traits>::addConfiguration()
    .registerDevice("v", 1);
}

} // namespace Vsrc
} // namespace Device
} // namespace Xyce

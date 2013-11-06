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

//----------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_CSW.C,v $
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
// Revision Number: $Revision: 1.26.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:38 $
//
// Current Owner  : $Author: tvrusso $
//----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------
#if 0

#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif

// ----------   Xyce Includes   ----------
#include <N_DEV_CSW.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_DeviceOptions.h>

#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>

namespace Xyce {
namespace Device {

vector< vector<int> > CSWInstance::jacStamp;

//-----------------------------------------------------------------------------
// Function      : CSW::CSW
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
CSW::CSW(SolverState & ss1,
                    DeviceOptions & do1)
 : Device(ss1,do1)
{
  name="switch";
  linearDeviceFlag = 1;
}

//-----------------------------------------------------------------------------
// Function      : CSW::CSW
// Purpose       : copy constructor
// Special Notes :
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
CSW::CSW(const CSW& right)
  : Device(right)
{

}

//-----------------------------------------------------------------------------
// Function      : CSW::~CSW
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
CSW::~CSW()
{
  vector<CSWModel*>::iterator iterM;
  vector<CSWModel*>::iterator firstM = modelContainer.begin ();
  vector<CSWModel*>::iterator lastM  = modelContainer.end ();

  // loop over models:
  for (iterM = firstM; iterM != lastM; ++iterM)
  {
    delete (*iterM);
  }
}

//-----------------------------------------------------------------------------
// Function      : CSW::addModel
// Purpose       : This function adds a model to the list of
//                 switch models.
//
// Special Notes : This is the  "public" version of the
//                 function, meaning that it returns a bool to indicate
//                 success or failure (rather than a pointer to the model).
//                 Also, it checks to see if the model being added is
//                 already part of the list.  If the model already exists,
//                 it generates an error.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/31/00
//-----------------------------------------------------------------------------
DeviceModel * CSW::addModel (const ModelBlock & MB)
{
  // Check to make sure that there doesn't already exist a model
  // of this name.
  vector<CSWModel*>::iterator iter;
  vector<CSWModel*>::iterator first = modelContainer.begin ();
  vector<CSWModel*>::iterator last  = modelContainer.end ();

  int icomp = 0;
  for (iter=first; iter!= last; ++iter)
    if ((*iter)->name == MB.name) icomp = 1;

  string msg;
  if (icomp == 1)
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "\n");
    msg = "\nCSW::addModel:\n";
    msg += "\tAttempted to add model that already exists.  ";
    msg += MB.name;
    msg += "\nIgnoring all but the first definition.\n";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_WARNING,msg);
    return 0;
  }
  else
  {
#ifdef Xyce_DEBUG_DEVICE
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_INFO_0, "\n");
    msg = "\nCSW::addModel:\n";
    msg += MB.name + " does not exist.  adding to list.";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_DEBUG_0,msg);
#endif
  }

  vector<CSWModel*>::iterator SWM_iter;
  return addModel_internal(MB,SWM_iter);

  //if (SWM_iter == modelContainer.end()) return 0;
  //else                            return 1;

}

//-----------------------------------------------------------------------------
// Function      : CSW::addModel_internal
// Purpose       : This function adds a model to the list of
//                 SW(switch) models.
//
// Special Notes : This is the "internal", private version of the function.
//                 It is called by both the public version (see above),
//                 and also by the addInstance function, when neccessary.
//
//                 This function does not check to see if the model already
//                 exists - it assumes that check has already been performed.
//                 Also, it returns a pointer to the model, which is needed
//                 by the addInstance function.
//
//
// Scope         : private
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/03/00
//-----------------------------------------------------------------------------
DeviceModel * CSW::addModel_internal
(const ModelBlock & MB, vector<CSWModel*>::iterator & SWM_iter)
{
  CSWModel * SWMptr = new CSWModel(MB,solState,devOptions);

  // add to the SW model list.
  modelContainer.push_back(SWMptr);
  SWM_iter = modelContainer.end()-1;

  return SWMptr;
}

//-----------------------------------------------------------------------------
// Function      : CSW::addInstance
// Purpose       : This function adds an instance to the list of
//                 SW instances.
//
// Special Notes : All device instances are contained in one of many lists.
//                 Each switch model owns exactly one list of instances.
//                 If an instance is added which does not reference a model,
//                 or that references a model that does not exist, then a
//                 new model has to be added.
//
//                 For now, this function will only add a model for instances
//                 that do not specifically reference one.  If an instance
//                 references a non-existant model, this function will
//                 return a fatal error.  (This function should only be called
//                 from the device manager function, addDeviceInstance, and
//                 that function will not call this one if the instance
//                 refers to a non-existant model.)
//
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/2/00
//-----------------------------------------------------------------------------
DeviceInstance * CSW::addInstance (
  InstanceBlock & IB,
  MatrixLoadData & mlData1,
  SolverState &ss1,
  ExternData  &ed1,
  DeviceOptions & do1)
{
  CSWInstance * SWIptr =
    new CSWInstance(IB,mlData1,ss1,ed1,do1);

  vector<CSWModel*>::iterator iter;
  vector<CSWModel*>::iterator first = modelContainer.begin ();
  vector<CSWModel*>::iterator last  = modelContainer.end ();

  ModelBlock MB;
  string modelName = IB.modelName;

  if (modelName == "")  // If NULL model name, use the default model.
    modelName = "CSWModel1";
  else
    modelName = IB.modelName;

  int ifound = 0;
  for (iter=first; iter != last; ++iter)
    if ((*iter)->name == modelName)
    {  ifound = 1; SWIptr->SWM_iter = iter; break; }

  if (ifound == 0)  // Instance references a non-existant model.
  {
    string msg = "CSW::addInstance could not find model ";
    msg += modelName + " which is referenced by instance " + IB.name;
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
  }

  SWIptr->name      = IB.name;
  SWIptr->modelName = modelName;

  (*(SWIptr->SWM_iter))->instanceContainer.push_back(SWIptr);
  return SWIptr;
}

//-----------------------------------------------------------------------------
// Function      : CSW::printOutModels
// Purpose       : debugging tool
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/03/00
//-----------------------------------------------------------------------------
void CSW::printOutModels ()
{

  cout << endl;

  vector<CSWModel*>::iterator iter;
  vector<CSWModel*>::iterator first = modelContainer.begin ();
  vector<CSWModel*>::iterator last  = modelContainer.end ();

  int isize = modelContainer.size();

  const string dashedline =
"-----------------------------------------------------------------------------";
  cout << dashedline << endl;
  cout << "Number of SW models: " << isize << endl;

  int i;
  for (iter=first,i=0; iter!=last; ++iter,++i)
  {
    cout << i << ": name = " << (*iter)->name << " type = " << (*iter)->type;
    cout << " level = " << (*iter)->level << endl;
    (*iter)->printOutInstances();
  }
  cout << dashedline << endl;

}

//-----------------------------------------------------------------------------
// Function      : CSW::getModelPointer
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/01/00
//-----------------------------------------------------------------------------
DeviceModel * CSW::getModelPointer (
  const string & modelNameRef, bool & bsuccess)
{
  DeviceModel * modelPtr = NULL;

  vector<CSWModel*>::iterator iter;
  vector<CSWModel*>::iterator begin = modelContainer.begin();
  vector<CSWModel*>::iterator end = modelContainer.end();

  for (iter=begin; iter!=end; ++iter)
  {
    ExtendedString tmpname((*iter)->name);
    if (tmpname == modelNameRef)
    { modelPtr = (*iter); bsuccess = true;}
  }

  return modelPtr;
}

//-----------------------------------------------------------------------------
// Function      : CSW::deleteInstance
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/08/01
//-----------------------------------------------------------------------------
bool CSW::deleteInstance (const string & tmpName)
{
  bool bsuccess = false;

  vector<CSWModel*>::iterator iterM;
  vector<CSWModel*>::iterator firstM = modelContainer.begin();
  vector<CSWModel*>::iterator lastM  = modelContainer.end();

  vector<CSWInstance*>::iterator iterI;
  vector<CSWInstance*>::iterator firstI;
  vector<CSWInstance*>::iterator lastI;

  // loop over models:
  for (iterM = firstM; iterM != lastM; ++iterM)
  {
    // loop over instances:
    firstI = (*iterM)->instanceContainer.begin();
    lastI  = (*iterM)->instanceContainer.end();

    for (iterI = firstI; iterI != lastI; ++iterI)
    {
      ExtendedString nameArg(tmpName);  nameArg.toUpper ();
      ExtendedString nameIter((*iterI)->name); nameIter.toUpper ();
      if (nameArg == nameIter)
      {
        bsuccess = true;
        delete (*iterI);
        (*iterM)->instanceContainer.erase(iterI);
        return bsuccess;
      }
    }
  }

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : CSWInstance::CSWInstance
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
CSWInstance::CSWInstance(
  InstanceBlock & IB,
  MatrixLoadData & mlData1,
  SolverState &ss1,
  ExternData  &ed1,
  DeviceOptions & do1)

  : DeviceInstance(IB,mlData1,ss1,ed1,do1),
    R(0.0),
    ON(false),
    G(0.0),
    SW_STATE(0.0),
    switch_state(0.0),
    li_switch_state(-1),
    switch_state_initialized(false),
    li_Pos(-1),
    li_Neg(-1),
    li_PosCntl(-1),
    li_NegCntl(-1),
    APosEquPosNodeOffset(-1),
    APosEquNegNodeOffset(-1),
    ANegEquPosNodeOffset(-1),
    ANegEquNegNodeOffset(-1)
{
  numIntVars = 0;
  numExtVars = 4;

  numStateVars = 1;

  vector<Param>::const_iterator iter;
  vector<Param>::const_iterator begin = IB.params.begin();
  vector<Param>::const_iterator end   = IB.params.end();

  for (iter=begin; iter !=  end; ++iter)
  {
    ExtendedString tagES(iter->tag());
    tagES.toUpper();

    if (tagES == "ON")   { if (iter->given() && iter->iVal()==1) ON=true; }
    else if (tagES == "OFF") { if (iter->given() && iter->iVal()==1) ON=false;}
#ifdef Xyce_DEBUG_DEVICE
    cout << " found tagES="<< tagES << endl;
    cout << " given is " ;
    if (iter->given())
      cout << "true" << endl;
    else
      cout << "false" << endl;

#endif
  }


  if( jacStamp.empty() )
  {
    jacStamp.resize(4);
    jacStamp[0].resize(2);
    jacStamp[0][0]=0;
    jacStamp[0][1]=1;
    jacStamp[1].resize(2);
    jacStamp[1][0]=0;
    jacStamp[1][1]=1;
  }
}

//-----------------------------------------------------------------------------
// Function      : CSWInstance::CSWInstance
// Purpose       : copy constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
CSWInstance::CSWInstance (const CSWInstance & right)
  : DeviceInstance(right),
    R(right.R),
    ON(right.ON),
    G(right.G),
    SW_STATE(right.SW_STATE),
    switch_state(right.switch_state),
    li_switch_state(right.li_switch_state),
    switch_state_initialized(right.switch_state_initialized),
    SWM_iter(right.SWM_iter),
    li_Pos(right.li_Pos),
    li_Neg(right.li_Neg),
    li_PosCntl(right.li_PosCntl),
    li_NegCntl(right.li_NegCntl),
    APosEquPosNodeOffset(right.APosEquPosNodeOffset),
    APosEquNegNodeOffset(right.APosEquNegNodeOffset),
    ANegEquPosNodeOffset(right.ANegEquPosNodeOffset),
    ANegEquNegNodeOffset(right.ANegEquNegNodeOffset)
{

}

//-----------------------------------------------------------------------------
// Function      : CSWInstance::CSWInstance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
CSWInstance::~CSWInstance ()
{

}

//-----------------------------------------------------------------------------
// Function      : CSWInstance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/20/02
//-----------------------------------------------------------------------------
void CSWInstance::registerLIDs( const vector<int> & intLIDVecRef,
                                     const vector<int> & extLIDVecRef)
{
  string msg;

  // Check if the size of the ID lists corresponds to the
  // proper number of internal and external variables.
  int numInt = intLIDVecRef.size();
  int numExt = extLIDVecRef.size();

  if (numInt != numIntVars)
  {
    msg = "CSWInstance::registerLIDs:";
    msg += "numInt != numIntVars";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
  }

  if (numExt != numExtVars)
  {
    msg = "CSWInstance::registerLIDs:";
    msg += "numExt != numExtVars";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
  }

  // copy over the global ID lists.
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  // now use these lists to obtain the indices into the
  // linear algebra entities.  This assumes an order.
  // For the matrix  indices, first do the rows.

  li_Pos = extLIDVec[0];
  li_Neg = extLIDVec[1];
  li_PosCntl = extLIDVec[2];
  li_NegCntl = extLIDVec[3];

}

//-----------------------------------------------------------------------------
// Function      : CSWInstance::registerStateLIDs
// Purpose       : Note that the SW does not have any state vars.
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/20/02
//-----------------------------------------------------------------------------
void CSWInstance::registerStateLIDs( const vector<int> & staLIDVecRef )
{
  string msg;

  // Check if the size of the ID lists corresponds to the
  // proper number of internal and external variables.
  int numSta = staLIDVecRef.size();

  if (numSta != numStateVars)
  {
    msg = "CSWInstance::registerLIDs:";
    msg += "numSta != numStateVars";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
  }

  // copy over the global ID lists:
  staLIDVec = staLIDVecRef;
  li_switch_state = staLIDVec[0];
}

//-----------------------------------------------------------------------------
// Function      : CSWInstance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 09/2/02
//-----------------------------------------------------------------------------
const vector< vector<int> > & CSWInstance::jacobianStamp() const
{
  return jacStamp;
}

//-----------------------------------------------------------------------------
// Function      : CSWInstance::registerJacLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 09/2/02
//-----------------------------------------------------------------------------
void CSWInstance::registerJacLIDs( const vector< vector<int> > & jacLIDVec )
{
  APosEquPosNodeOffset = jacLIDVec[0][0];
  APosEquNegNodeOffset = jacLIDVec[0][1];
  ANegEquPosNodeOffset = jacLIDVec[1][0];
  ANegEquNegNodeOffset = jacLIDVec[1][1];
}

//-----------------------------------------------------------------------------
// Function      : CSWInstance::updatePrimaryStateBlock
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/29/01
//-----------------------------------------------------------------------------
bool CSWInstance::updatePrimaryStateBlock ()
{
  double * staVec = extData.nextStaVectorRawPtr;

  bool bsuccess = updateIntermediateVarsBlock ();

  //  obtain the current value of the switch state
  switch_state = SW_STATE;
  // we keep this around so updateIntermediateVarsBlock knows whether to bother
  // getting the "last" state
  switch_state_initialized=true;

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline2 = "---------------------";
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << dashedline2 << endl;
    cout << "  CSW::udpatePrimaryStateBlock\n";
    cout << "  name   = " << name << endl;
    cout << "  switch_state = " << switch_state << endl;
  }
#endif

  staVec[li_switch_state] = switch_state;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : CSWInstance::updateSecondaryStateBlock
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 01/29/01
//-----------------------------------------------------------------------------
bool CSWInstance::updateSecondaryStateBlock ()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : CSWInstance::updateIntermediateVarsBlock
// Purpose       : update intermediate variables for one switch instance
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, Component Information and Models
// Creation Date : 1/10/01
//-----------------------------------------------------------------------------
bool CSWInstance::updateIntermediateVarsBlock()
{
  bool bsuccess = true;

  double v_pos_cntl, v_neg_cntl;
  double current_state;
  double v_ctrl;

  double * solVec = extData.nextSolVectorRawPtr;

#ifdef Xyce_DEBUG_DEVICE
  const string dashedline2 = "---------------------";
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << endl << dashedline2 << endl;
    cout << "  name = " << name << endl;
  }
#endif

  // get solution vector indices:
  v_pos_cntl = solVec[li_PosCntl];
  v_neg_cntl = solVec[li_NegCntl];

  v_ctrl = v_pos_cntl - v_neg_cntl;
  // This is not really correct, an interim hack.  This is supposed
  // to be where we deal with the specification of ON or OFF from the
  // netlist.
  if ((getSolverState().dcopFlag && getSolverState().newtonIter==0))
    {
      if (ON)
	current_state = 1;
      else
	current_state = 0;
#ifdef Xyce_DEBUG_DEVICE
      cout << " SETTING INITIAL CONDITION TO " << current_state << endl;
#endif
    }
  else
    current_state = (v_ctrl-(*SWM_iter)->VOFF)*(*SWM_iter)->VdInv;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << "  v_ctrl is " << v_ctrl << endl;
    cout << "  VdInv (von-voff) is " << (*SWM_iter)->VdInv << endl;
    cout << "  current_state (vc-voff)/(von-voff) is ";
    cout << current_state << endl;
  }
#endif

  if (current_state >= 1.0)
  {
    R = (*SWM_iter)->RON;
    G = 1.0/R;
  }
  else if ( current_state <= 0.0)
  {
    R = (*SWM_iter)->ROFF;
    G = 1.0/R;
  }
  else
  {
	current_state = 2*current_state - 1;
    G = exp(-(*SWM_iter)->Lm - 0.75*(*SWM_iter)->Lr*current_state +
		0.25*(*SWM_iter)->Lr*current_state*current_state*current_state);
    R = 1.0/G;

//    calculations for derivatives with respect to v_ctrl
//    K = -(*SWM_iter)->Lm - 0.75*(*SWM_iter)->Lr*current_state +
//        0.25*(*SWM_iter)->Lr*current_state*current_state*current_state;
//    G = exp(K);
//    dGdv_ctrl = K * exp(K) *
//                ( 0.75 * (*SWM_iter)->Lr * ( current_state*current_state - 1 ) ) *
//                (*SWM_iter)->VdInv;
//    dIdv_ctrl = dGdv_ctrl * ( v_pos - v_neg );

  }

  v_pos = solVec[li_Pos];
  v_neg = solVec[li_Neg];

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    cout << " v_pos=" << v_pos;
    cout << " v_neg=" << v_neg << endl;
    cout << "  G = " << G << endl;
    cout << "  R = " << R << endl;
  }
#endif

  return bsuccess;

}

//-----------------------------------------------------------------------------
// Function      : CSWModel::CSWModel
// Purpose       : copy constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
CSWModel::CSWModel (const CSWModel & right)
  : DeviceModel(right),
    RON(right.RON),
    ROFF(right.ROFF),
    VON(right.VON),
    VOFF(right.VOFF),
    Lm(right.Lm),
    Lr(right.Lr),
    VdInv(right.VdInv)
{

}

//-----------------------------------------------------------------------------
// Function      : CSWModel::CSWModel
// Purpose       : copy constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
CSWModel::CSWModel (const ModelBlock & MB,
                                    SolverState & ss1,
                                    DeviceOptions & do1)
  : DeviceModel(MB,ss1,do1),
  RON(0.0),
  ROFF(0.0),
  VON(0.0),
  VOFF(0.0)
{
  double del;

  vector<Param>::const_iterator iter;
  vector<Param>::const_iterator begin = MB.params.begin();
  vector<Param>::const_iterator end   = MB.params.end();

  for (iter=begin; iter !=  end;  ++iter)
  {
    ExtendedString tagES(iter->tag());
    tagES.toUpper();

    if (tagES == "RON")		RON	=	iter->dVal();
    if (tagES == "ROFF")	ROFF	=	iter->dVal();
    if (tagES == "VON")		VON	=	iter->dVal();
    if (tagES == "VOFF")	VOFF	=	iter->dVal();
  }

  del = VON-VOFF;
  if (del < 0 && del > -1e-12)
    del = -1e-12;
  if (del >= 0 && del < 1e-12)
    del = 1e-12;
  VdInv = 1.0/del;
  Lm = log (sqrt(RON*ROFF));
  Lr = log (RON/ROFF);
}

//-----------------------------------------------------------------------------
// Function      : CSWModel::~CSWModel
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
CSWModel::~CSWModel ()
{
  vector<CSWInstance*>::iterator iter;
  vector<CSWInstance*>::iterator first = instanceContainer.begin();
  vector<CSWInstance*>::iterator last  = instanceContainer.end();

  for (iter=first; iter!=last; ++iter)
    delete (*iter);
}

//-----------------------------------------------------------------------------
// Function      : CSWModel::printOutInstances
// Purpose       : debugging tool.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/03/00
//-----------------------------------------------------------------------------
void CSWModel::printOutInstances ()
{
  vector<CSWInstance*>::iterator iter;
  vector<CSWInstance*>::iterator first = instanceContainer.begin();
  vector<CSWInstance*>::iterator last  = instanceContainer.end();

  int i;
  cout << endl;
  cout << "    name     modelName  Parameters" << endl;
  for (i=0, iter=first; iter!=last; ++iter, ++i)
  {
    cout << "  " << i << ": " << (*iter)->name << "      ";
    cout << (*iter)->modelName;
    cout << "    R = " << (*iter)->R;
    cout << "  G = " << (*iter)->G;
    cout << "  State = " << (*iter)->SW_STATE;
    cout << endl;
  }

  cout << endl;
}

} // namespace Device
} // namespace Xyce

#endif // if 0


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
// Filename       : $RCSfile: N_DEV_DeviceInstance.C,v $
//
// Purpose        : Implementation of the base device instance class.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/30/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.169.2.4 $
//
// Revision Date  : $Date: 2013/10/03 17:23:38 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------
#include <Xyce_config.h>

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_ExternData.h>
#include <N_ERH_ErrorMgr.h>
#include <N_LAS_Matrix.h>
#include <N_LAS_Vector.h>
#include <N_DEV_Const.h>
#include <N_DEV_DeviceMgr.h>
#include <N_DEV_ExternDevice.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_NumericalJacobian.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::DeviceInstance
// Purpose       : instance block constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/30/00
//-----------------------------------------------------------------------------
DeviceInstance::DeviceInstance(
  InstanceBlock &       instance_block,
  MatrixLoadData &      matrix_load_data,
  SolverState &         solver_state,
  ExternData &          extern_data,
  DeviceOptions &       device_options)
  : DeviceEntity(solver_state, device_options, instance_block.getName(), instance_block.netlistFileName_, instance_block.lineNumber_),
    modelName_(instance_block.getModelName()),
    psLoaded(false),
    ssLoaded(false),
    rhsLoaded(false),
    origFlag(true),
    PDEDeviceFlag(false),
    mlData(matrix_load_data),
    cols(matrix_load_data.cols),
    vals(matrix_load_data.vals),
    numIntVars(0),
    numExtVars(2),
    numStateVars(0),
    numStoreVars(0),
    numLeadCurrentStoreVars(0),
    configuredForLeadCurrent(false),
    loadLeadCurrent(false),
    mergeRowColChecked(false),
    extData(extern_data),
    numJacPtr(NULL)
{
  devConMap.resize(2);
  devConMap[0] = 1;
  devConMap[1] = 1;
  numJacPtr = new NumericalJacobian (matrix_load_data, solver_state, extern_data, device_options);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::~DeviceInstance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/30/00
//-----------------------------------------------------------------------------
DeviceInstance::~DeviceInstance()
{
  delete numJacPtr;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::getDepSolnVars
//
// Purpose       : This function configures a device to an auxiliary F & Q vector
//                 so that lead currents can be calculated for this device.c.
//
// Special Notes : This must be called soon after the constructor call
//                 before the store vector is allocated.
//
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 03/22/13
//-----------------------------------------------------------------------------
void DeviceInstance::enableLeadCurrentCalc()
{
  if(!configuredForLeadCurrent)
  {
    // indicated that this device is now configured for lead current calculation
    // this avoids claiming too much space in the store vector if this function
    // is called more than once for a device.
    configuredForLeadCurrent = true;

    // set device instance flag to indicate the need to load lead current
    // data into store F & Q vectors
    loadLeadCurrent = true;

    // request additional space in store vector for lead current vars
    numStoreVars = numStoreVars +  numLeadCurrentStoreVars;
  }
  return;
}


//-----------------------------------------------------------------------------
// Function      : N_DEV_DeviceInstance::getDepSolnVars
//
// Purpose       : Topology uses this method to check for late dependencies
//                 due to such things as Expressions in the B-src.
//
// Special Notes : Returns empty list for devices that use this base method
//
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/05/01
//-----------------------------------------------------------------------------
const vector<string> & DeviceInstance::getDepSolnVars()
{
  return expVarNames;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::registerDepSolnLIDs
// Purpose       : Allows registration of LIDs of nodes and instances that
//                 appear in expressions that occur in the device.
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 03/15/05
//-----------------------------------------------------------------------------
void DeviceInstance::registerDepSolnLIDs
(const vector< vector<int> > & depSolnLIDVecRef)
{
  int size = expVarLIDs.size();
  if (size != depSolnLIDVecRef.size())
  {
    string msg = "DeviceInstance::registerDepSolnLIDs: Inconsistent ";
    msg += "number of LIDs returned from topology";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL, msg);
  }
  for (int i = 0; i < size; ++i)
  {
    if (depSolnLIDVecRef[i].size() != 1)
    {
      string msg = "Problem with value for: " + expVarNames[i];
      msg += ".  This may be an incorrect usage of a lead current in place of a current ";
      msg += "through a voltage source.";
      std::ostringstream oss;
      oss << "Error in " << netlistLocation() << "\n" << msg;
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, oss.str());
    }
    expVarLIDs[i] = depSolnLIDVecRef[i][0];
  }

  return;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::registerDepSolnGIDs
//
// Purpose       : This function allows global ID's to be registered
//                 with a device.  These global ID's refer to the global
//                 indices for the solution vector.  Given these ID's, the
//                 device can then also determine the (row,col) ID's needed
//                 to load the  jacobian matrix.
//
// Special Notes : The method is used for late resolution of variables
//                 of dependency such as for B-source.  Does nothing for
//                 devices using this base class method.
//
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/05/01
//-----------------------------------------------------------------------------
void DeviceInstance::registerDepSolnGIDs(
  const vector< vector<int> > & varList )
{
  int size = expVarGIDs.size();
  for (int i = 0; i < size; ++i)
  {
    expVarGIDs[i] = varList[i][0];
  }
  return;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::registerJacLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date :
//-----------------------------------------------------------------------------
void DeviceInstance::registerJacLIDs( const vector< vector<int> > & jacLIDVec )
{
  if (getDeviceOptions().numericalJacobianFlag || getDeviceOptions().testJacobianFlag)
  {
    devJacLIDs = jacLIDVec;
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::getDepSolnGIDVec
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 03/17/05
//-----------------------------------------------------------------------------
void DeviceInstance::getDepSolnGIDVec( vector<int> & depGIDVec )
{
  depGIDVec = expVarGIDs;
  return;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::getDepStateVars
//
// Purpose       : Topology uses this method to check for late dependencies
//                 due to such things as Expressions in the B-src.
//
// Special Notes : Returns empty list for devices that use this base method
//
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/05/01
//-----------------------------------------------------------------------------
const vector<string> & DeviceInstance::getDepStateVars()
{
  static vector<string> emptyList;
  emptyList.clear();
  return emptyList;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::registerDepStateGIDs
//
// Purpose       : This function allows global ID's to be registered
//                 with a device.  These global ID's refer to the global
//                 indices for the solution vector.  Given these ID's, the
//                 device can then also determine the (row,col) ID's needed
//                 to load the  jacobian matrix.
//
// Special Notes : The method is used for late resolution of variables
//                 of dependency such as for B-source.  Does nothing for
//                 devices using this base class method.
//
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/05/01
//-----------------------------------------------------------------------------
void DeviceInstance::registerDepStateGIDs(
  const vector< vector<int> > & varList )
{
  if( varList.size() != 0 )
  {
    string msg;
    msg = "DeviceInstance::registerDepStateGIDs\n";
    msg += "\tCall to registerDepStateGIDs for a device which doesn't use it.";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::getDepStoreVars
//
// Purpose       : Topology uses this method to check for late dependencies
//                 due to such things as Expressions in the B-src.
//
// Special Notes : Returns empty list for devices that use this base method
//
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
const vector<string> & DeviceInstance::getDepStoreVars()
{
  static vector<string> emptyList;
  emptyList.clear();
  return emptyList;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::registerDepStoreGIDs
//
// Purpose       : This function allows global ID's to be registered
//                 with a device.  These global ID's refer to the global
//                 indices for the solution vector.  Given these ID's, the
//                 device can then also determine the (row,col) ID's needed
//                 to load the  jacobian matrix.
//
// Special Notes : The method is used for late resolution of variables
//                 of dependency such as for B-source.  Does nothing for
//                 devices using this base class method.
//
// Scope         : public
// Creator       : Eric Keiter
// Creation Date :
//-----------------------------------------------------------------------------
void DeviceInstance::registerDepStoreGIDs(
  const vector< vector<int> > & varList )
{
  if( varList.size() != 0 )
  {
    string msg;
    msg = "DeviceInstance::registerDepStoreGIDs\n";
    msg += "\tCall to registerDepStoreGIDs for a device which doesn't use it.";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::getIndexPairList
//
// Purpose       : This function
//
// Special Notes :
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/26/00
//-----------------------------------------------------------------------------
bool DeviceInstance::getIndexPairList(list<index_pair> & iplRef)
{
  iplRef = indexPairList;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::testDAEMatrices
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 12/15/06
//-----------------------------------------------------------------------------
bool DeviceInstance::testDAEMatrices (vector<string> & nameVec)
{
  bool bsuccess = true;

  // if necessary, consolodate the LIDs vector.
  if (devLIDs.empty())
  {
    devLIDs = extLIDVec;
    devLIDs.insert(devLIDs.end(), intLIDVec.begin(), intLIDVec.end());
    devLIDs.insert(devLIDs.end(), expVarLIDs.begin(), expVarLIDs.end());
  }

  bsuccess = numJacPtr-> testDAEMatrices (*this, nameVec);

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::trivialStampLoader
// Purpose       : This function contains most of the original
//                 loadTrivialMatrixStamp function.  See comments above.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/25/05
//-----------------------------------------------------------------------------
bool DeviceInstance::trivialStampLoader (N_LAS_Matrix * matPtr)
{
  int count = 0;
  int i;
  int row;

  vector<int>::const_iterator firstVar;
  vector<int>::const_iterator lastVar;
  vector<int>::const_iterator iterVar;
  int localRows = matPtr->getLocalNumRows();

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << endl;
    cout << "Loading trivial stamp for " << getName() << endl;
  }
#endif

  if (cols.size() < 1) cols.resize(1);
  if (vals.size() < 1) vals.resize(1);

  for (i=0;i<2;++i)
  {
    // do external vars first.
    if (i==0)
    {
      firstVar = extLIDVec.begin ();
      lastVar  = extLIDVec.end ();
    }
    // then do internal vars.
    else
    {
      firstVar = intLIDVec.begin ();
      lastVar  = intLIDVec.end ();
    }

    for (iterVar=firstVar; iterVar!=lastVar; ++iterVar)
    {
      row = *iterVar;

#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 0)
        cout << "matrix row = " << row << endl;
#endif
      if (row < 0)
      {
#ifdef Xyce_DEBUG_DEVICE
        if (getDeviceOptions().debugLevel > 0) cout << "\tNOT loading this one - too small" << endl;
#endif
        continue;
      }
#if 0
      if (row >= localRows)
      {
#ifdef Xyce_DEBUG_DEVICE
        if (getDeviceOptions().debugLevel > 0) cout << "\tNOT loading this one - too big" << endl;
#endif
        continue;
      }
#endif

#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 0) cout << "\tloading this one" << endl;
#endif
      count = 1;
      vals[0] = 1.0;
      cols[0] = row;

      //matPtr->putLocalRow(row, count, &vals[0], &cols[0]);
      matPtr->replaceLocalRow(row, count, &vals[0], &cols[0]);
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::loadTrivialDAE_FMatrixStamp
// Purpose       : See loadTrivialMatrixStamp - this is the same thing, except
//                 for the new-DAE F-Matrix.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/25/05
//-----------------------------------------------------------------------------
bool DeviceInstance::loadTrivialDAE_FMatrixStamp ()
{
  return trivialStampLoader (extData.dFdxMatrixPtr);
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::zeroMatrixDiagonal
// Purpose       : puts zeros into on matrix diagonal, but just for this
//                 device.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 05/25/05
//-----------------------------------------------------------------------------
bool DeviceInstance::zeroMatrixDiagonal (N_LAS_Matrix * matPtr)
{
  int count = 0;
  int i;
  int row;

  vector<int>::const_iterator firstVar;
  vector<int>::const_iterator lastVar;
  vector<int>::const_iterator iterVar;
  int localRows = matPtr->getLocalNumRows();

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    cout << endl;
    cout << "Zeroing the matrix diagonal for " << getName() << endl;
  }
#endif

  if (cols.size() < 1) cols.resize(1);
  if (vals.size() < 1) vals.resize(1);

  for (i=0;i<2;++i)
  {
    // do external vars first.
    if (i==0)
    {
      firstVar = extLIDVec.begin ();
      lastVar  = extLIDVec.end ();
    }
    // then do internal vars.
    else
    {
      firstVar = intLIDVec.begin ();
      lastVar  = intLIDVec.end ();
    }

    for (iterVar=firstVar; iterVar!=lastVar; ++iterVar)
    {
      row = *iterVar;

#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 0)
        cout << "matrix row = " << row << endl;
#endif
      if (row < 0) continue;
      if (row >= localRows) continue;

#ifdef Xyce_DEBUG_DEVICE
      if (getDeviceOptions().debugLevel > 0)
        cout << "\tloading this one" << endl;
#endif
      count = 1;
      vals[0] = 0.0;
      cols[0] = row;

      matPtr->putLocalRow(row, count, &vals[0], &cols[0]);
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::spiceInternalName
// Purpose       : Converts a conventional "Xyce style" internal variable
//                 name to a spice style name.
// Special Notes :
//                 Example:
//                 chilespice:    d:subckt1:test1_internal
//                 xyce:          xsubckt1:dtest1_internal
//
//                 This function will not remove subcircuit "x" characters,
//                 but does move the first character of the local name
//                 to the beginning of the full string.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 09/12/01
//-----------------------------------------------------------------------------
void DeviceInstance::spiceInternalName (string & tmpname)
{
  size_t pos=tmpname.find_last_of(":");
  if (pos != string::npos)
  {
    tmpname=tmpname.substr(pos+1,1)+":"+tmpname.substr(0,pos+1)+tmpname.substr(pos+2);
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::enablePDEContinuation
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/23/03
//-----------------------------------------------------------------------------
bool DeviceInstance::enablePDEContinuation()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::disablePDEContinuation
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/23/03
//-----------------------------------------------------------------------------
bool DeviceInstance::disablePDEContinuation()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::setPDEContinuationAlpha
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/23/03
//-----------------------------------------------------------------------------
void DeviceInstance::setPDEContinuationAlpha (double alpha)
{
  return;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::setPDEContinuationBeta
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/28/04
//-----------------------------------------------------------------------------
void DeviceInstance::setPDEContinuationBeta (double beta)
{
  return;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::setInitialGuess
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/23/03
//-----------------------------------------------------------------------------
bool DeviceInstance::setInitialGuess ()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::getMaxTimeStepSize
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/23/03
//-----------------------------------------------------------------------------
double DeviceInstance::getMaxTimeStepSize  ()
{
  return getDeviceOptions().defaultMaxTimeStep;
}


//-----------------------------------------------------------------------------
// Function      : DeviceInstance::jacStampMap
// Purpose       : Compute Jacobian Stamp and Map for devices that can have merged nodes
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 2/13/04
//-----------------------------------------------------------------------------
void DeviceInstance::jacStampMap
(vector< vector<int> > & stamp_parent,
 vector<int> & map_parent,
 vector< vector<int> > & map2_parent,
 vector< vector<int> > & stamp,
 vector<int> & map,
 vector< vector<int> > & map2,
 int from, int to, int original_size)
{
  int i, j, t_index, f_index, t_index2, f_index2, p_size, p_row, f_mod;
  int map_2_from;
  bool new_col;
  vector< vector<int> > fill;
  vector<int> dup;
  vector< vector<int> > map2_tmp;

  if (from <= to)
  {
    string msg = "Internal Error 1 in DeviceInstance::jacStampMap. from <= to.";
    //msg += " from index = " << from;
    //msg += " to index = " << to;
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
  }

  if (map_parent.size() == 0)
  {
    map_parent.resize(original_size);
    map2_parent.resize(original_size);
    for (i=0 ; i<original_size ; ++i)
    {
      map_parent[i] = i;
      map2_parent[i].resize(stamp_parent[i].size());
      for (j=0 ; j<stamp_parent[i].size() ; ++j)
      {
        map2_parent[i][j] = j;
      }
    }
  }
  map2.resize(original_size);

  // This is to merge the column that is being eliminated into the column it is merged.
  // If extra elements are present then we must increment the map2 value for later
  // elements.  There are multiple cases depending what is populated.
  for (i=0 ; i<original_size ; ++i)
  {
    f_index = -1;
    t_index = -1;
    p_row = map_parent[i];
    for (j=0 ; j<stamp_parent[p_row].size() ; ++j)
    {
      if (stamp_parent[p_row][j] == from)
        f_index = j;
      if (stamp_parent[p_row][j] == to)
        t_index = j;
    }
    map2[i].resize(map2_parent[i].size());
    f_index2 = -1;
    t_index2 = -1;
    for (j=0 ; j<map2_parent[i].size() ; ++j)
    {
      map2[i][j] = map2_parent[i][j];
      if (stamp_parent[p_row][map2[i][j]] == from)
        f_index2 = j;
      if (stamp_parent[p_row][map2[i][j]] == to)
        t_index2 = j;
    }
    if (f_index >= 0)
    {
      if (t_index >= 0)
      {
        for (j=0 ; j<map2[i].size() ; ++j)
        {
          if (map2[i][j] > f_index)
            map2[i][j]--;
        }
        if (f_index2 >= 0 && t_index2 >= 0)
          map2[i][f_index2] = map2[i][t_index2];
      }
      else
      {
        for (j=0 ; j<map2[i].size() ; ++j)
        {
          if (stamp_parent[p_row][map2[i][j]] > to && stamp_parent[p_row][map2[i][j]] < from)
            ++map2[i][j];
        }
        t_index = 0;
        for (j=0 ; j<stamp_parent[p_row].size() ; ++j)
        {
          if (to > stamp_parent[p_row][j])
            ++t_index;
          else
            break;
        }
        if (f_index2 >= 0)
          map2[i][f_index2] = t_index;
      }
    }
  }
  map.resize(original_size);
  p_size = stamp_parent.size();

  // This is to merge the row that is being eliminated into the row it is merged into
  // if extra elements are present then we must increment the map2 value for later
  // elements in the row
  map2_tmp = map2;

  for (i=0 ; i<stamp_parent[from].size() && stamp_parent[from][i] < p_size-1 ; ++i)
  {
    new_col = false;
    for (j=0 ; j<stamp_parent[to].size() ; ++j)
    {
      if (j == 0)
        new_col = true;
      if (stamp_parent[from][i] == stamp_parent[to][j])
      {
        new_col = false;
        break;
      }
    }
    if (new_col)
    {
      for (j=0 ; j<map2[to].size() ; ++j)
      {
        if (stamp_parent[to][map2_tmp[to][j]] > stamp_parent[from][i])
        {
          ++map2[to][j];
        }
      }
    }
  }

  stamp.resize(p_size-1);

  dup.resize(original_size);
  f_mod = from;
  for (i=1 ; i<original_size ; ++i)
  {
    dup[i] = -1;
    for (j=0 ; j<i ; ++j)
    {
      if (map_parent[i] == map_parent[j])
      {
        dup[i] = j;
        if (i<=f_mod)
          ++f_mod;
        break;
      }
    }
  }

  for (i=0 ; i<f_mod ; ++i)
    map[i] = map_parent[i];
  map[f_mod] = map[to];
  for (i=f_mod+1 ; i<original_size ; ++i)
    map[i] = map_parent[i]-1;

  for (i=1 ; i<original_size ; ++i)
  {
    if (dup[i] >= 0)
      map[i] = map[dup[i]];
  }


  // Now that we know where the row originally came from, we can do any renumbering
  // needed in the map2 source row
  map2_tmp = map2;
  map_2_from = from;
  if (map[from] != map[to]) {
    for (i=from+1 ; i<original_size ; ++i) {
      if (map[i] == map[to]) {
        map_2_from = i;
        break;
      }
    }
    if (map_2_from == from) {
      string msg = "Internal Error 2 in DeviceInstance::jacStampMap";
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
    }
  }

  for (i=0 ; i<stamp_parent[to].size() && stamp_parent[to][i] < p_size-1 ; ++i)
  {
    new_col = false;
    for (j=0 ; j<stamp_parent[from].size() ; ++j)
    {
      if (j == 0)
        new_col = true;
      if (stamp_parent[to][i] == stamp_parent[from][j])
      {
        new_col = false;
        break;
      }
    }
    if (new_col)
    {
      for (j=0 ; j<map2[map_2_from].size() ; ++j)
      {
        if (stamp_parent[from][map2_tmp[map_2_from][j]] > stamp_parent[to][i])
        {
          ++map2[map_2_from][j];
        }
      }
    }
  }

  fill.resize(p_size);
  for (i=0 ; i<p_size ; ++i)
  {
    fill[i].resize(p_size);
    for (j=0 ; j<p_size ; ++j)
      fill[i][j] = 0;
    for (j=0 ; j<stamp_parent[i].size() ; ++j)
      fill[i][stamp_parent[i][j]] = 1;
  }
  for (i=0 ; i<p_size ; ++i)
  {
    fill[to][i] += fill[from][i];
  }
  for (i=0 ; i<p_size ; ++i)
  {
    fill[i][to] += fill[i][from];
  }
  for (i=from ; i<p_size-1 ; ++i)
  {
    for (j=0 ; j<p_size ; ++j)
      fill[i][j] = fill[i+1][j];
    for (j=0 ; j<p_size ; ++j)
      fill[j][i] = fill[j][i+1];
  }
  for (i=0 ; i<p_size-1 ; ++i)
  {
    stamp[i].clear();
    for (j=0 ; j<p_size-1 ; ++j)
    {
      if (fill[i][j] > 0)
        stamp[i].push_back(j);
    }
  }
  return;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::jacStampMap_fixOrder
//
// Purpose       : This function corrects the compressed row column array so
//                 that the column indices are in ascending order.
//
// Special Notes : The reason for this function is that the
//                 DeviceInstance::jacStampMap function
//                 implicitly requires an ordered jacStamp to work correctly.
//                 Some devices (particularly devices that have meshes),
//                 will start out with an non-ordered stamp, at least in the
//                 column entries.
//
//                 Note, this does require the "map" argument, because,
//                 unlike the function DeviceInstnace::jacStamMap,
//                 because this function doesn't change the row ordering,
//                 or remove or merge any rows.
//
//                 This function only changes the column ordering in the
//                 compressed row form of the jacStamp.  It thus requires
//                 modifications to map2, which is essentially a column map.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 2/01/08
//-----------------------------------------------------------------------------
void DeviceInstance::jacStampMap_fixOrder
(vector< vector<int> > & stamp_parent,
 vector< vector<int> > & map2_parent,
 vector< vector<int> > & stamp,
 vector< vector<int> > & map2)
{
  int i, j;
  int current_size = stamp_parent.size();

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
  {
    cout << "------------------------------------------------------" << endl;
    cout << "Begin DeviceInstance::jacStampMap_fixOrder."<<endl;
    cout << "------------------------------------------------------" << endl;
  }
#endif

  // if this is the first time this function is called (for a particular device), then
  // allocate the map and set up their trivial contents.
  if (map2_parent.size() == 0)
  {
    map2_parent.resize(current_size);
    for (i=0 ; i<current_size ; ++i)
    {
      map2_parent[i].resize(stamp_parent[i].size());
      for (j=0 ; j<stamp_parent[i].size() ; ++j)
      {
        map2_parent[i][j] = j;
      }
    }
  }

  stamp.clear();
  map2.clear();

  // To make this simple, start out with a full, dense stamp.
  vector < vector <int> > denseStamp(current_size);
  for (i=0;i<current_size;++i)
  {
    denseStamp[i].resize(current_size,-1);

    for (j=0;j<stamp_parent[i].size();++j)
    {
      int denseCol = stamp_parent[i][j];
      denseStamp[i][denseCol] = j;
    }
  }

  // At this point, the denseStamp has been set up.  Now use it to re-create the
  // compressed-row stamp.  By simply looping over the dense stamp, the column order
  // in the compressed row stamp will automatically be ascending.
  // map2 is set up here as well, by pulling the values out that we previously put into
  // dense stamp.
  stamp.resize(current_size);
  map2.resize(current_size);
  for (i=0;i<current_size;++i)
  {
    for (j=0;j<current_size;++j)
    {
      int colMapIndex = denseStamp[i][j];
      if (colMapIndex!=-1)
      {
        stamp[i].push_back(j);
      }
    }

    int stampRowLength=stamp[i].size();
    map2[i].resize(stampRowLength, 0);

    int k=0;
    for (j=0;j<current_size;++j)
    {
      int colMapIndex = denseStamp[i][j];
      if (colMapIndex!=-1 && colMapIndex < stampRowLength)
      {
        map2[i][colMapIndex] = k;
        ++k;
      }
    }
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 1 && getSolverState().debugTimeFlag)
  {
    cout << "From inside of DeviceInstance::jacStampMap_fixOrder:"<<endl;
    cout << "The original parent stamp is:" << endl;
    outputJacStamp(stamp_parent);
    cout << "The new reduced stamp is:" << endl;
    outputJacStamp(stamp);
    cout << "The dense stamp is:" << endl;
    outputJacStamp(denseStamp);

    cout << "The new map is:" << endl;
    outputJacStamp(map2);

    cout << "------------------------------------------------------" << endl;
    cout << "End DeviceInstance::jacStampMap_fixOrder."<<endl;
    cout << "------------------------------------------------------" << endl;
  }
#endif

  return;
}
//-----------------------------------------------------------------------------
// Function      : DeviceInstance::outputJacStamp
// Purpose       : Output jacStamp (for debugging)
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/19/06
//-----------------------------------------------------------------------------
void DeviceInstance::outputJacStamp(const vector<vector<int> > & jac)
{
  int i,j;
  for (i=0 ; i<jac.size() ; ++i)
  {
    cout << "Row: " << i << " ::";
    for (j=0 ; j<jac[i].size() ; ++j)
      cout << "  " << jac[i][j];
    cout << endl;
  }
  cout << endl;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::outputJacMaps
// Purpose       : Output jacMap and jacMap2 (for debugging)
// Special Notes :
// Scope         : public
// Creator       : Keith Santarelli, Electrical & Microsystems Modeling
// Creation Date : 02/20/08
//-----------------------------------------------------------------------------
void DeviceInstance::outputJacMaps(const vector<int>  & jacMap,
					 const vector<vector<int> > & jacMap2)
{
  int i,j;

  for (i=0 ; i<jacMap.size() ; ++i)
  {
    cout << "Row " << i << ": ";
    for (j=0; j < jacMap2[i].size(); j++)
    {
      cout << jacMap[i]<< "," << jacMap2[i][j] << " ";
    }
    cout << endl;
  }

  cout << endl;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::registerGIDData
// Purpose       : Insert GID data into 'indexPairList' object
//
// Special Notes : This information is neccessary for the numerical Jacobian.
//                 It is only called once.
//
//                 The numerical jacobian may get confused by duplicate
//                 matrix entries, in that it might load them 2x.  For that
//                 reason, this function checks for duplicates.
//
// Scope         : public
// Creator       : Rob Hoekstra
// Creation Date : 12/13/04
//-----------------------------------------------------------------------------
void DeviceInstance::registerGIDData(
  const vector<int> & counts,
  const vector<int> & GIDs,
  const vector< vector<int> > & jacGIDs )
{

  // the test for the PDE system flag only needs to be here until the
  // "diagonal" matrix load is fixed.
  //if (getDeviceOptions().numericalJacobianFlag || getSolverState().PDESystemFlag)
  if (getDeviceOptions().numericalJacobianFlag)
  {
    indexPairList.clear();

    int extSize = counts[0];
    int intSize = counts[1];
    int expSize = counts[2];
    int size = GIDs.size();
    int i;

    map<int,int> testIndexMap;

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      cout << "DeviceInstance::registerGIDData for " << getName() << endl;
      cout << "  extSize      = " << extSize  << endl;
      cout << "  intSize      = " << intSize  << endl;
      cout << "  expSize      = " << expSize  << endl;
      cout << "  GIDs.size    = " << size << endl;
      cout << "  jacGIDs.size = " << jacGIDs.size () << endl;
      cout << endl;
    }
#endif

    // Copy out external gids:
    extGIDList.clear ();
    for (i = 0; i < extSize; ++i )
    {
      if ( testIndexMap.find(GIDs[i]) == testIndexMap.end() )
      {
        extGIDList.push_back( index_pair( GIDs[i], 1 ) );
        testIndexMap[GIDs[i]] = 1;
      }
    }

    testIndexMap.clear ();

    // Copy out internal gids:
    intGIDList.clear ();
    for (; i<intSize+extSize;++i)
    {
      if ( testIndexMap.find(GIDs[i]) == testIndexMap.end() )
      {
        intGIDList.push_back( index_pair( GIDs[i], 1 ) );
        testIndexMap[GIDs[i]] = 1;
      }
    }

    testIndexMap.clear ();

    // Copy out the exp var gid's, if they exist.  These will
    // only exist in devices which depend on expressions, like the Bsrc.
    expVarGIDs.clear ();
    for (; i<intSize+extSize+expSize;++i)
    {
      if ( testIndexMap.find(GIDs[i]) == testIndexMap.end() )
      {
        expVarGIDs.push_back( GIDs[i] );
        testIndexMap[GIDs[i]] = 1;
      }
    }

    // Now copy the exp var GIDs into the extVarGID's.
    // add the contents of expVarGIDs to the extGIDListRef.
    // This is done because the numerical jacobian treats expression
    // GIDs the same as external (nodal) GIDs.
    int expS = expVarGIDs.size();
    for (int itmp = 0; itmp < expS; ++itmp)
    {
      extGIDList.push_back( index_pair( expVarGIDs[itmp], 1 ) );
    }

    testIndexMap.clear ();

    // do the index pairs for the jacobian matrix
    indexPairList.clear ();
    for(i = 0; i < jacGIDs.size () ; ++i )
    {

      if ( testIndexMap.find(GIDs[i]) == testIndexMap.end() )
      {
        testIndexMap[GIDs[i]] = 1;

        map<int,int> testJMap;
        testJMap.clear ();
        int length = jacGIDs[i].size();
        for( int j = 0; j < length; ++j )
        {
          if ( testJMap.find(jacGIDs[i][j]) == testJMap.end () )
          {
            indexPairList.push_back( index_pair( GIDs[i], jacGIDs[i][j] ) );
            testJMap[jacGIDs[i][j]] = 1;
          }
        }
      }
    }

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      vector<int>::const_iterator begin = GIDs.begin ();
      vector<int>::const_iterator end   = GIDs.end   ();
      vector<int>::const_iterator iterGID;

      cout << " Complete GIDs :" << endl;
      for (iterGID=begin;iterGID!=end ;++iterGID)
      {
        cout << "\tgid=" << *iterGID;
        cout << endl;
      }
      cout << endl;

      list<index_pair>::iterator first;
      list<index_pair>::iterator last;
      list<index_pair>::iterator iter;

      cout << " intGIDList :" << endl;
      first = intGIDList.begin ();
      last  = intGIDList.end   ();
      for (iter=first;iter!=last;++iter)
      {
        cout << "\tgid=" <<iter->row;
        cout << endl;
      }
      cout << endl;

      cout << " extGIDList :" << endl;
      first = extGIDList.begin ();
      last  = extGIDList.end   ();
      for (iter=first;iter!=last;++iter)
      {
        cout << "\tgid=" <<iter->row;
        cout << endl;
      }
      cout << endl;


      cout << " indexPairList :" << endl;
      first = indexPairList.begin();
      last  = indexPairList.end  ();
      for (iter=first;iter!=last;++iter)
      {
        cout << "  row=" <<iter->row;
        cout << "  col=" <<iter->col;
        cout << endl;
      }
      cout << endl << endl;
    }
#endif

  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::getNumNodes
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 9/16/05
//-----------------------------------------------------------------------------
int DeviceInstance::getNumNodes() const
{
  const Configuration &cTab = getMyParametricData().getConfigTable();

  return cTab.numNodes;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::getNumOptionalNodes
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 9/16/05
//-----------------------------------------------------------------------------
int DeviceInstance::getNumOptionalNodes() const
{
  const Configuration &cTab = getMyParametricData().getConfigTable();

  return cTab.numOptionalNodes;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::getNumFillNodes
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 9/16/05
//-----------------------------------------------------------------------------
int DeviceInstance::getNumFillNodes() const
{
  const Configuration &cTab = getMyParametricData().getConfigTable();

  return cTab.numFillNodes;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::getPrimaryParameter
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 9/16/05
//-----------------------------------------------------------------------------
void DeviceInstance::getPrimaryParameter(string & primaryParameter) const
{
  const Configuration &cTab = getMyParametricData().getConfigTable();

  primaryParameter = cTab.primaryParameter;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::getModelTypes
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 9/16/05
//-----------------------------------------------------------------------------
void DeviceInstance::getModelTypes(vector<string> & modTypesVector) const
{
  const Configuration &cTab = getMyParametricData().getConfigTable();

  modTypesVector = cTab.modelTypes;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::updateTemperature
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/30/00
//-----------------------------------------------------------------------------
bool DeviceInstance::updateTemperature(const double & temp_tmp)
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 06/03/02
//-----------------------------------------------------------------------------
bool DeviceInstance::processParams (string param)
{
  string msg;
  msg = "DeviceInstance::processParams does not exist";
  msg += "  for this device: " + getName() + "\n";
  N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);

  return true;
}

} // namespace Device
} // namespace Xyce

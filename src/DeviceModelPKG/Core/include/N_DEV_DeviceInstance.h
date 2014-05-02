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
// Filename       : $RCSfile: N_DEV_DeviceInstance.h,v $
//
// Purpose        : This file contains the device instance base class.
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
// Revision Number: $Revision: 1.179.2.1 $
//
// Revision Date  : $Date: 2014/02/26 20:16:30 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_DeviceInstance_h
#define Xyce_N_DEV_DeviceInstance_h

#include <list>
#include <map>
#include <string>
#include <vector>

#include <N_DEV_fwd.h>
#include <N_DEV_DeviceEntity.h>
#include <N_DEV_DeviceSupport.h>
#include <N_UTL_Misc.h>

class N_LAS_Matrix;
class N_LAS_MultiVector;
class N_LAS_Vector;
class N_LAS_Solver;

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : DeviceInstance
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
class DeviceInstance  : public DeviceEntity
{
private:
  DeviceInstance();

public:
  DeviceInstance(
     const InstanceBlock &     instance_block,
     ParametricData<void> &    parametric_data,
     const FactoryBlock &      factory_block);

  virtual ~DeviceInstance();

private:
  DeviceInstance(const DeviceInstance &);
  DeviceInstance &operator=(const DeviceInstance &);

public:
  // This function configures the device to request space in the store
  // vector for lead current calculations.  It must be called soon
  // after the constructor call before the store vector is allocated.
  virtual void enableLeadCurrentCalc();

  virtual void registerGIDs(
     const std::list<index_pair> & intGIDListRef,
     const std::list<index_pair> & extGIDListRef ) {}

  virtual void registerStateGIDs(const std::list<index_pair> & staGIDListRef) {}

  virtual void registerStoreGIDs(const std::list<index_pair> & stoGIDListRef) {}

  virtual void registerLIDs( const std::vector<int> & intLIDVecRef,
                             const std::vector<int> & extLIDVecRef ) {}
  virtual void registerStateLIDs( const std::vector<int> & staLIDVecRef ) {}

  virtual void registerStoreLIDs( const std::vector<int> & stoLIDVecRef ) {}

  virtual const std::vector<std::string> & getDepSolnVars();
  virtual void registerDepSolnGIDs( const std::vector< std::vector<int> > & varList );
  virtual const std::vector<std::string> & getDepStateVars();
  virtual void registerDepStateGIDs( const std::vector< std::vector<int> > & varList );
  virtual const std::vector<std::string> & getDepStoreVars();
  virtual void registerDepStoreGIDs( const std::vector< std::vector<int> > & varList );

  virtual void registerDepSolnLIDs( const std::vector< std::vector<int> > & depSolnLIDVecRef );
  virtual void registerDepStateLIDs( const std::vector< std::vector<int> > & depStaLIDVecRef ) {}
  virtual void registerDepStoreLIDs( const std::vector< std::vector<int> > & depStoLIDVecRef ) {}

  virtual const std::vector< std::vector<int> > & jacobianStamp() const 
  {
    static std::vector< std::vector<int> > dummy;
    return dummy;
  }

  virtual void registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec );

  virtual void registerGIDData(
     const std::vector<int> & counts,
     const std::vector<int> & GIDs,
     const std::vector< std::vector<int> > & jacGIDs );

  virtual void setupPointers()
  {}

  virtual void getDepSolnGIDVec( std::vector<int> & depGIDVec );

  virtual bool getIndexPairList(std::list<index_pair> & iplRef);

  virtual bool getInstanceBreakPoints (std::vector<N_UTL_BreakPoint> &breakPointTimes);

  virtual bool updateSource ();

  virtual bool processParams ();

  virtual bool processInstanceParams () 
  {
    return true;
  }

  virtual bool updateTemperature(const double & temp_tmp);

  virtual bool isConverged();

  virtual bool testDAEMatrices (std::vector<std::string> & nameVec);

  virtual bool loadTrivialDAE_FMatrixStamp ();
  bool trivialStampLoader (N_LAS_Matrix * matPtr);
  bool zeroMatrixDiagonal (N_LAS_Matrix * matPtr);

  virtual bool updateIntermediateVars () = 0;
  virtual bool updatePrimaryState ();
  virtual bool updateSecondaryState ();
  virtual bool setIC ();

  // This indicates if the device has functions that can output plot files for
  // internal variables.
  virtual bool plotfileFlag () {return false;}

  // load zeros into mask for equations that should not be used
  // to compute error estimates.  Return true if any zeros set.
  // Default implementation just does nothing (leaves everything 1.0)
  virtual bool loadDeviceMask() {return false;}

  // tell device instance that current solution has been accepted at
  // current time.  Most devices don't care, but the transmission line
  // does.
  virtual void acceptStep() {}

  // new DAE functions:
  virtual bool loadDAEQVector ()=0;
  virtual bool loadDAEFVector ()=0;

  virtual bool loadDAEdQdx ()=0;
  virtual bool loadDAEdFdx ()=0;

  int getNumIntVars() const 
  {
    return numIntVars;
  }

  int getNumExtVars() const 
  {
    return numExtVars;
  }

  int getNumStateVars() const 
  {
    return numStateVars;
  }

  int getNumStoreVars() const 
  {
    return numStoreVars;
  }

  void setNumStoreVars(int num_store_vars) 
  {
    numStoreVars = num_store_vars;
  }

  virtual void getDevConMap(std::vector<int> &);

  virtual DeviceState * getInternalState();
  virtual bool setInternalState( const DeviceState & state );

  virtual bool loadDFDV(int iElectrode, N_LAS_Vector * dfdvPtr);
  virtual bool calcConductance (int iElectrode, const N_LAS_Vector * dxdvPtr);

  // internal name map stuff:
  virtual std::map<int, std::string> & getIntNameMap ();
  virtual std::map<int, std::string> & getStateNameMap ();
  virtual std::map<int, std::string> & getStoreNameMap ();

  void spiceInternalName (std::string & tmpname);

  virtual bool outputPlotFiles () {return true;}

  // two level newton and PDE-continuation
  virtual bool enablePDEContinuation();
  virtual bool disablePDEContinuation();
  virtual void setPDEContinuationAlpha (double alpha);
  virtual void setPDEContinuationBeta  (double beta );

  virtual bool setInitialGuess ();
  virtual double getMaxTimeStepSize  ();

  virtual void varTypes( std::vector<char> & varTypeVec ) {}

protected:
  void jacStampMap(
     std::vector< std::vector<int> > & stamp_parent,
     std::vector<int> & map_parent,
     std::vector< std::vector<int> > & map2_parent,
     std::vector< std::vector<int> > & stamp,
     std::vector<int> & map,
     std::vector< std::vector<int> > & map2,
     int from, int to, int original_size);

  void jacStampMap_fixOrder(
     std::vector< std::vector<int> > & stamp_parent,
     std::vector< std::vector<int> > & map2_parent,
     std::vector< std::vector<int> > & stamp,
     std::vector< std::vector<int> > & map2);

  void outputJacStamp(const std::vector<std::vector<int> > & jac);
  void outputJacMaps(const std::vector<int>  & jacMap, const std::vector<std::vector<int> > & jacMap2);

public:
  bool getOrigFlag() const 
  {
    return origFlag;
  }

  void setOrigFlag(bool origFlag_local) 
  {
    origFlag = origFlag_local;
  }

  const std::vector<int> &getDevLIDs() const 
  {
    return devLIDs;
  }

  const std::vector<std::vector<int> > &getDevJacLIDs() const 
  {
    return devJacLIDs;
  }

  const std::vector<int> &getStaLIDVec() const 
  {
    return staLIDVec;
  }

  bool getMergeRowColChecked() const 
  {
    return mergeRowColChecked;
  }
  void setMergeRowColChecked(bool mergeRowColChecked_local) 
  {
    mergeRowColChecked = mergeRowColChecked_local;
  }

  const MatrixLoadData &getMatrixLoadData() const 
  {
    return mlData;
  }

  MatrixLoadData &getMatrixLoadData() 
  {
    return mlData;
  }

private:
  MatrixLoadData &      mlData;

protected:
  const ExternData & extData;
  std::list<index_pair> intGIDList;
  std::list<index_pair> extGIDList;
  std::list<index_pair> indexPairList;

  std::list<index_pair> staGIDList;

  std::vector<int> intLIDVec;
  std::vector<int> extLIDVec;

  std::vector<int> staLIDVec;
  std::vector<int> stoLIDVec;

  // devLIDs is a combined LID vector, containing int, ext, and expVar ID's.
  std::vector<int> devLIDs;
  std::vector< std::vector<int> > devJacLIDs;

  std::map<int,std::string> intNameMap;
  std::map<int,std::string> stateNameMap;
  std::map<int,std::string> storeNameMap;

  // device support class: (limiter functions, etc.)
  DeviceSupport devSupport;

private:
  bool configuredForLeadCurrent;  // flag which indicates that numStoreVars already includes numLeadCurrentStoreVars

public:
  std::vector<int> & cols;
  std::vector<double> & vals;

  NumericalJacobian * numJacPtr;

  bool psLoaded;
  bool ssLoaded;
  bool rhsLoaded;

  bool origFlag;

  int numIntVars;
  int numExtVars;
  int numStateVars;
  int numStoreVars;

  int numLeadCurrentStoreVars;
  bool loadLeadCurrent;           // flag indicating that we want to load lead current data during F & Q load

  std::vector<int> devConMap;

  bool mergeRowColChecked;
};

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 11/30/00
//-----------------------------------------------------------------------------
inline bool DeviceInstance::updateSecondaryState ()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::setIC
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 07/22/03
//-----------------------------------------------------------------------------
inline bool DeviceInstance::setIC ()
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::getInternalState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 09/02/01
//-----------------------------------------------------------------------------
inline DeviceState * DeviceInstance::getInternalState()
{
  return NULL;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::getIntNameMap
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 09/07/01
//-----------------------------------------------------------------------------
inline std::map<int, std::string> & DeviceInstance::getIntNameMap ()
{
  return intNameMap;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::getStateNameMap
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 07/31/2012
//-----------------------------------------------------------------------------
inline std::map<int, std::string> & DeviceInstance::getStateNameMap ()
{
  return stateNameMap;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::getStoreNameMap
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 07/31/2012
//-----------------------------------------------------------------------------
inline std::map<int, std::string> & DeviceInstance::getStoreNameMap ()
{
  return storeNameMap;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::getDevConMap
// Purpose       : Get connectivity map for leads.  Zero means a lead is
//                 connected to ground.  Other values indicate subsets of
//                 leads that have sonnection to each other.  Example would
//                 be a mosfet which would have 1 for gate and 2 for drain
//                 and source and zero for bulk, assuming bulk is grounded
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 06/20/05
//-----------------------------------------------------------------------------
inline void DeviceInstance::getDevConMap(std::vector<int> & conMap)
{
  conMap = devConMap;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::loadDFDV
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 12/03/02
//-----------------------------------------------------------------------------
inline bool DeviceInstance::loadDFDV(int iElectrode, N_LAS_Vector * dfdvPtr)
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::calcConductance
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 12/03/02
//-----------------------------------------------------------------------------
inline bool DeviceInstance::calcConductance (int iElectrode, const N_LAS_Vector * dxdvPtr)
{
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance:isConverged ()
// Purpose       : Return whether a device has done something that should
//                  be interpreted as invalidating other convergence tests
//                  (i.e. that means this step should not be considered
//                   converged even if norms are good)
//                 Since origFlag is set to true by the DeviceInstance
//                 constructor, this is a suitable base class method for
//                 almost all devices.  Devices with more complex convergence
//                 issues can override.
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Component Information and Models
// Creation Date : 03/22/05
//-----------------------------------------------------------------------------
inline bool DeviceInstance::isConverged()
{
  return origFlag;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::getInstanceBreakPoints
// Purpose       : virtual function for obtaining breakpoints from a device.
//
// Special Notes : No-op for the base class version.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 2/05/06
//-----------------------------------------------------------------------------
inline bool DeviceInstance::getInstanceBreakPoints(std::vector<N_UTL_BreakPoint> &breakPointTimes)
{
  return false;
}

//-----------------------------------------------------------------------------
// Function      : DeviceInstance::updateSource
// Purpose       : virtual function for obtaining breakpoints from a device.
//
// Special Notes : No-op for the base class version.
//
// Scope         : public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 2/05/06
//-----------------------------------------------------------------------------
inline bool DeviceInstance::updateSource ()
{
  return true;
}

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::DeviceInstance N_DEV_DeviceInstance;

#endif


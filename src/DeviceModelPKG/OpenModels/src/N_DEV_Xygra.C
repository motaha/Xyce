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
// Filename       : $RCSfile: N_DEV_Xygra.C,v $
//
// Purpose        :  Provide a linking device for coupling to Alegra
//                   simulations.
//
// Special Notes  :  Alegra uses API calls to set conductance and source
//                   information used to construct F (or RHS), and to retrieve
//                   nodal voltages.
//
// Creator        : Tom Russo, SNL, Electrical and Microsystems Modeling
//
// Creation Date  : 08/18/08
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.42.2.2 $
//
// Revision Date  : $Date: 2014/03/06 23:33:44 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------
#include <Xyce_config.h>

// ---------- Standard Includes ----------
#include <N_UTL_Misc.h>

#include <sstream>

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_DeviceMaster.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_Xygra.h>
#include <N_DEV_Message.h>
#include <N_ERH_ErrorMgr.h>

#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>

namespace Xyce {
namespace Device {

template<>
ParametricData<XygraCoilData>::ParametricData()
{
  addPar ("NAME", "COIL0", &XygraCoilData::name);
  addPar ("NUMWINDINGS", 1, &XygraCoilData::numWindings);
}

ParametricData<XygraCoilData> &XygraCoilData::getParametricData() {
  static ParametricData<XygraCoilData> parMap;

  return parMap;
}


//-----------------------------------------------------------------------------
// Function      : XygraCoilData::XygraCoilData
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 9/11/2008
//-----------------------------------------------------------------------------
XygraCoilData::XygraCoilData()
  : CompositeParam(getParametricData()),
    name("coil"),
    numWindings(1)
{}

//-----------------------------------------------------------------------------
// Function      : XygraCoilData::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 9/11/2008
//-----------------------------------------------------------------------------
void XygraCoilData::processParams ()
{
}

#ifdef Xyce_DEBUG_DEVICE
//-----------------------------------------------------------------------------
// Function      : operator<<
// Purpose       : "<<" operator for XygraCoilData
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 9/11/2008
//-----------------------------------------------------------------------------
std::ostream & operator<<(std::ostream & os, const XygraCoilData & xcd)
{
  os << " XygraCoilData for: name = " << xcd.getName() <<
    " numWindings=" << xcd.getNumWindings() <<
    std::endl;

  return os;
}
#endif




namespace Xygra {


void Traits::loadInstanceParameters(ParametricData<Xygra::Instance> &p)
{
  // Set up configuration constants:
p.addComposite("COIL", XygraCoilData::getParametricData(), &Xygra::Instance::coilDataMap);
}

void Traits::loadModelParameters(ParametricData<Xygra::Model> &p)
{}


// Class Instance
//-----------------------------------------------------------------------------
// Function      : Instance::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo
// Creation Date : 8/18/2008
//-----------------------------------------------------------------------------
bool Instance::processParams()
{
  // If there are any time dependent parameters, set their values at for
  // the current time.

  // now set the temperature related stuff.
  //updateTemperature(temp);

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateTemperature
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, Component Information and Models
// Creation Date : 8/18/2008
//-----------------------------------------------------------------------------
bool Instance::updateTemperature ( const double & temp)
{
  bool bsuccess = true;
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::Instance
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/18/2008
//-----------------------------------------------------------------------------
Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock &         IB,
  Model &                       Miter,
  const FactoryBlock &          factory_block)

  : DeviceInstance(IB, configuration.getInstanceParameters(), factory_block),
    model_(Miter),
    t0_(0),
    t1_(0)
{
  numExtVars   = IB.numExtVars;  // we have as many as were specified on the
                                 // instance line
  numIntVars   = 0;
  numStateVars = 0;

  // Start with both vectors empty.  We'll only make them the right size
  // when they get set from setConductances or setSources, so we can tell
  // whether to bother using them.  If not the right size, no contribution
  // from those terms.
  theSourceVector_.clear();
  theConductanceMatrix_.clear();
  theKMatrix_.clear();
  k0_.clear();
  k1_.clear();
  s0_.clear();
  s1_.clear();


  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to instance line and constant defaults from metadata:
  setParams (IB.params);

  // Set any non-constant parameter defaults:

  //if (!given("TEMP"))
  //  temp = getDeviceOptions().temp.dVal();

  // Calculate any parameters specified as expressions:
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors:
  processParams ();


  int totalNumIntVars=0;
  totalNumWindings = 0;
  if (!coilDataVec.empty())
  {
#ifdef Xyce_DEBUG_DEVICE
    Xyce::dout() << " We were given a coil spec." << std::endl;
    Xyce::dout() << "   There were " << coilDataVec.size() << " coils." << std::endl;
#endif

    // numExtVars is apparently 0 for the default device (which is never
    // used except for parameter output) and its associated coilDataVec always
    // has one.  So skip this test if numExtVars==0.

    if ( numExtVars != 0 && numExtVars != coilDataVec.size()*2)
    {
      UserError0(*this) << "Xygra Device " << getName() << "has "<< coilDataVec.size() << " coils  and "  << numExtVars << " external nodes.  Number of external nodes should be twice number of coils.";
    }

    nCoils=coilDataVec.size();
    nWindings.resize(nCoils);
    coilNames.resize(nCoils);
    for ( int i = 0; i < nCoils; ++i)
    {
#ifdef Xyce_DEBUG_DEVICE
      Xyce::dout() << " Coil["<< i << "] name is " << coilDataVec[i]->getName() << ", has " << coilDataVec[i]->getNumWindings() << " windings. " << std::endl;
#endif // Xyce_DEBUG_DEVICE
      totalNumIntVars += coilDataVec[i]->getNumWindings() - 1;
      nWindings[i] = coilDataVec[i]->getNumWindings();
      coilNames[i] = coilDataVec[i]->getName();
      totalNumWindings += nWindings[i];
    }

#ifdef Xyce_DEBUG_DEVICE
    Xyce::dout() << " We would have " << totalNumIntVars << " internal vars " << std::endl;
#endif
  }
  else
  {
    // if no coil spec given, each pair of external nodes corresponds to a
    // one-winding coil.
    nCoils=numExtVars/2;
    nWindings.resize(nCoils,1);
    totalNumWindings = nCoils;
    if (nCoils*2 != numExtVars)
    {
      std::ostringstream ost;

      ost << "Instance::Instance:";
      ost << "Number of nodes given to device " << getName() << "is not even."
          << std::endl;
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, ost.str() );
    }

  }
  // we'll use these in mapping between winding currents/dV's and nodal
  // contributions
  coilExtStart.resize(nCoils);
  coilIntStart.resize(nCoils);
  windingNodes.resize(totalNumWindings);

  // set up numIntVars:
  numIntVars = totalNumIntVars;

  // Now that we have computed our numIntVars, it's OK to
  // set up jacStamp.  For now, we'll have to assume a full (dense) Jacobian
  // maybe there's a way to get this done later, but I think not, because
  // topology needs to know it too soon, and we won't know the real structure
  // until it's too late to use it.
  //
  setupJacStamp_();

  // set up numStateVars:
  numStateVars = 0;
}

//-----------------------------------------------------------------------------
// Function      : Instance::~Instance
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/18/2008
//-----------------------------------------------------------------------------
Instance::~Instance()
{
}

//-----------------------------------------------------------------------------
// Function      : Instance::setupJacStamp_
// Purpose       : Utility function to set up the jacobian stamp
// Special Notes : Just makes a dense jacobian stamp for now
// Scope         : private
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/18/08
//-----------------------------------------------------------------------------
void Instance::setupJacStamp_()
{
  int numNodes = numExtVars+numIntVars;
  jacStamp_.resize(numNodes);
  for (int i=0; i<numNodes; ++i)
  {
    jacStamp_[i].resize(numNodes);
    for (int j=0; j<numNodes; ++j)
    {
      jacStamp_[i][j]=j;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::setConductances
// Purpose       : take matrix of conductances, copy into our local array
// Special Notes : The matrix is dense
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/18/08
//-----------------------------------------------------------------------------
bool Instance::setConductances(const std::vector< std::vector<double> > &cM )
{
  int numNodes = numExtVars+numIntVars;
  // Sanity check:
  if (cM.size() != numNodes)
  {
    std::ostringstream ost;

    ost << "Instance::setConductances:";
    ost << "  Input matrix passed to device " << getName()
        << " (" << cM.size()
        << ") does not have number of rows required by netlist specification ("
        << numNodes << ")." << std::endl;
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, ost.str());
  }
  theConductanceMatrix_.resize(numNodes);
  for (int i=0; i<numNodes; ++i)
  {
    if (cM[i].size() != numNodes)
    {
      std::ostringstream ost;

      ost << "Instance::setConductances:";
      ost << "  row " << i << "of matrix passed to device " << getName()
          << " has " << cM[i].size()
          << " columns instead of number required by netlist specification ("
          << numNodes << ")." << std::endl;
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, ost.str());
    }

#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << std::endl << subsection_divider << std::endl;
      Xyce::dout() << " Device " << getName()  << " setConductances called for time " << getSolverState().currTime << std::endl;
      Xyce::dout() << std::endl << subsection_divider << std::endl;
    }
#endif

    theConductanceMatrix_[i].resize(numNodes);
    theConductanceMatrix_[i]=cM[i];
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setK
// Purpose       : take matrix of K values, copy into our local array
// Special Notes : The matrix is dense
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/18/08
//-----------------------------------------------------------------------------
bool Instance::setK(const std::vector< std::vector<double> > &kM ,
                    const double t)
{
  std::vector<std::vector<double> > * kPtr;

  // Decide which of our two K's we're setting
  if (t==0  || t == getSolverState().currTime)
  {
    kPtr = &k0_;
    t0_=t;
  }
  else
  {
    kPtr = &k1_;
    k0_=k1_; // copy old "future" K to "current" K, because we're resetting
    // future
    t0_=t1_;  // save old "future" time
    t1_=t;   // set new one
  }

  // Sanity check:
  if (kM.size() != totalNumWindings)
  {
    std::ostringstream ost;

    ost << "Instance::setK:";
    ost << "  Input matrix passed to device " << getName()
        << " (" << kM.size()
        << ") does not have number of rows required by netlist specification ("
        << totalNumWindings << ")." << std::endl;
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, ost.str());
  }

  for (int i=0; i<totalNumWindings; ++i)
  {
    if (kM[i].size() != totalNumWindings)
    {
      std::ostringstream ost;

      ost << "Instance::setK:";
      ost << "  row " << i << "of matrix passed to device " << getName()
          << " has " << kM[i].size()
          << " columns instead of number required by netlist specification ("
          << totalNumWindings << ")." << std::endl;
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, ost.str());
    }
  }

  // Sanity check is done, the matrix we're given fits our needs.  Now just
  // copy it.
  (*kPtr) = kM;

  // If this is the first call, we have to set K1 properly.
  if (t==0)
  {
    k1_ = k0_; // initialize k1
    t1_=t0_;
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << std::endl << subsection_divider << std::endl;
    Xyce::dout() << " Device " << getName()  << " setK called for time " << getSolverState().currTime << std::endl;
    Xyce::dout() << std::endl << subsection_divider << std::endl;
  }
#endif
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::getVoltages
// Purpose       : Copy the values of this devices' nodal voltages
//                 into the provided array.
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/18/08
//-----------------------------------------------------------------------------
bool Instance::getVoltages( std::vector< double > &nV )
{
  bool bsuccess=true;
  N_LAS_Vector * solVectorPtr = extData.nextSolVectorPtr;
  int numNodes = numExtVars + numIntVars;
  nV.resize(numNodes);

  // We will output these in a different order than we have them stored,
  // so that all of a coil's nodes are in a predictable order:
  // Coil1 positive
  // [coil1 internal nodes]
  // coil1 negative
  // coil2 positive [...etc....]
  int offset = 0;
  for (int coil=0; coil < nCoils; ++coil)
  {
    nV[offset++] = (*solVectorPtr)[li_Nodes_[coilExtStart[coil]]];
    for (int winding=0; winding< nWindings[coil]-1; ++winding)
    {
      nV[offset++] = (*solVectorPtr)[li_Nodes_[coilIntStart[coil]+winding]];
    }
    nV[offset++] = (*solVectorPtr)[li_Nodes_[coilExtStart[coil]+1]];
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::getCoilWindings
// Purpose       : return individual winding counts for each coil in provided
//                 array
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 9/15/08
//-----------------------------------------------------------------------------
void Instance::getCoilWindings( std::vector< int > &cW )
{
  cW = nWindings;
  return;
}

//-----------------------------------------------------------------------------
// Function      : Instance::getCoilNames
// Purpose       : return individual names for each coil in provided
//                 array
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 9/29/08
//-----------------------------------------------------------------------------
void Instance::getCoilNames( std::vector< std::string > &cN )
{
  cN = coilNames;
  return;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setSources
// Purpose       : take vector of source terms, copy into our local array
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/18/08
//-----------------------------------------------------------------------------
bool Instance::setSources(const std::vector< double > &sV,
                          const double t)
{
  int numNodes = numExtVars + numIntVars;

  std::vector<double> * sPtr;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << std::endl << subsection_divider << std::endl;
    Xyce::dout() << " Device " << getName()  << " setSources called for time " << getSolverState().currTime << " with value t= " << t << std::endl;
    Xyce::dout() << std::endl << subsection_divider << std::endl;
  }
#endif

  // Decide which of our two S's we're setting
  if (t==0  || t == getSolverState().currTime)
  {
    sPtr = &s0_;
    t0_=t;
#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << " setting source s0 for t=" << t << std::endl;
    }
#endif
  }
  else
  {
    sPtr = &s1_;
    s0_=s1_; // copy old "future" S to "current" S, because we're resetting
    // future
    t0_=t1_;  // save old "future" time
    t1_=t;   // set new one
#ifdef Xyce_DEBUG_DEVICE
    if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
    {
      Xyce::dout() << " setting source s1 for t=" << t << std::endl;
    }
#endif
  }

  // Sanity check:
  if (sV.size() != totalNumWindings)
  {
    std::ostringstream ost;

    ost << "Instance::setSources:";
    ost << "  Input vector passed to device " << getName()
        << " (" << sV.size()
        << ") does not have number of rows required by netlist specification ("
        << totalNumWindings << ")." << std::endl;
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, ost.str());
  }

  // Sanity check passes, copy the vector
  (*sPtr) = sV;
  // If this is the first call, we have to set S1 properly.
  if (t==0)
  {
    s1_ = s0_; // initialize k1
    t1_=t0_;
  }


  return true;
}

//-----------------------------------------------------------------------------
// Function      : Instance::interpolateSandK_
// Purpose       : use s0,k0,t0, s1, k1, t1 to interpolate S and K to current
//                 time
// Special Notes :
// Scope         : private
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 9/15/2008
//-----------------------------------------------------------------------------
void Instance::interpolateSandK_()
{
  // if S has been set at all:
  if (!s0_.empty());
  {
    theSourceVector_=s0_;
    if (getSolverState().currTime > t0_ && !s1_.empty() && t0_ != t1_)
    {
      double fac=(getSolverState().currTime - t0_)/(t1_-t0_);
      for (int i=0;i<totalNumWindings; ++i)
      {
        theSourceVector_[i] += (s1_[i]-s0_[i])*fac;
#ifdef Xyce_DEBUG_DEVICE
        if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
        {
          Xyce::dout() << " interpolated source for t=" << getSolverState().currTime << std::endl;
          Xyce::dout() << " source vector is: " << std::endl;
          Xyce::dout() << " s0["<<i<<"] = " << s0_[i] << std::endl;
          Xyce::dout() << " s1["<<i<<"] = " << s1_[i] << std::endl;
          Xyce::dout() << " fac = " << fac << std::endl;
          Xyce::dout() << " s["<<i<<"] = " << theSourceVector_[i] << std::endl;
        }
#endif
      }
    }
  }

  // if K has been set at all:
  if (!k0_.empty());
  {
    theKMatrix_=k0_;
    if (getSolverState().currTime > t0_ && !k1_.empty() && t0_ != t1_)
    {
      double fac=(getSolverState().currTime - t0_)/(t1_-t0_);
      for (int i=0;i<totalNumWindings; ++i)
      {
        for (int j=0;i<totalNumWindings; ++i)
        {
          theKMatrix_[i][j] += (k1_[i][j]-k0_[i][j])*fac;
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/18/2008
//-----------------------------------------------------------------------------
void Instance::registerLIDs(const std::vector<int> & intLIDVecRef,
                            const std::vector<int> & extLIDVecRef)
{
  AssertLIDs(intLIDVecRef.size() == numIntVars);
  AssertLIDs(extLIDVecRef.size() == numExtVars);

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0)
  {
    Xyce::dout() << std::endl << section_divider << std::endl;
    Xyce::dout() << "  Instance::registerLIDs" << std::endl;
    Xyce::dout() << "  name = " << getName() << std::endl;
  }
#endif

  int numNodes = numExtVars + numIntVars;

  // copy over the global ID lists.
  intLIDVec = intLIDVecRef;
  extLIDVec = extLIDVecRef;

  li_Nodes_.resize(numNodes);

  int extVar=0;
  int intVar=0;
  int theVar=0;

  // First comes the external nodes
  for (int coil=0; coil<nCoils; ++coil)
  {
    coilExtStart[coil] = theVar;
    li_Nodes_[theVar++] = extLIDVec[extVar++];
    li_Nodes_[theVar++] = extLIDVec[extVar++];
  }
  // now the internals

  for (int coil=0; coil<nCoils; ++coil)
  {
    coilIntStart[coil] = theVar;
    for( int i=0; i< nWindings[coil]-1; ++i)
    {
      li_Nodes_[theVar++] = intLIDVec[intVar++];
    }
  }

  // This doesn't really belong in registerLIDs, but since it's only here that
  // we've calculated everything we need to know, and we need this stuff
  // later, this is a good place to do it.

  int globalWinding=0;
  for (int coil=0; coil<nCoils; ++coil)
  {
    for ( int coilWinding=0; coilWinding<nWindings[coil]; coilWinding++)
    {
      int posNode;
      int negNode;
      if (coilWinding==0)
      {
        posNode=coilExtStart[coil];
      }
      else
      {
        posNode=coilIntStart[coil]+(coilWinding-1);
      }
      if (coilWinding==nWindings[coil]-1)
      {
        negNode=coilExtStart[coil]+1;
      }
      else
      {
        negNode=coilIntStart[coil]+coilWinding;
      }
      windingNodes[globalWinding++] = std::pair<int,int>(posNode,negNode);
    }
  }

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 )
  {
    for (int i=0; i<numNodes; ++i)
    {
      Xyce::dout() << "  li_Nodes_[" <<i << "] = " << li_Nodes_[i] << std::endl;
    }

    for (int winding=0; winding<totalNumWindings; winding++)
    {
      Xyce::dout() << "Winding " << winding << " between node "
           << windingNodes[winding].first << " and "
           << windingNodes[winding].second << std::endl;
    }
  }
#endif

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 )
  {
    Xyce::dout() << section_divider << std::endl;
  }
#endif
}

//-----------------------------------------------------------------------------
// Function      : Instance::getIntNameMap
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/18/2008
//-----------------------------------------------------------------------------
std::map<int,std::string> & Instance::getIntNameMap ()
{
  // set up the internal name map, if it hasn't been already.

  if (intNameMap.empty() && numIntVars != 0)
  {
    std::string tmpstr;
    // set up the internal name map
    std::ostringstream ost;

    for (int coil = 0; coil< nCoils; ++coil)
    {
      if (nWindings[coil]>1)
      {
        // this coil has internal nodes
        for (int nodeOffset=1; nodeOffset<nWindings[coil]; ++nodeOffset)
        {
          int localIndex=li_Nodes_[coilIntStart[coil]+nodeOffset-1];
          ost.str("");
          ost << getName() << "_coil" << coil << "_Internal"<<nodeOffset;
          tmpstr = ost.str();
          spiceInternalName(tmpstr);
          intNameMap[localIndex] = tmpstr;
        }
      }
    }
  }


  return intNameMap;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerStateLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/18/2008
//-----------------------------------------------------------------------------
void Instance::registerStateLIDs( const std::vector<int> & staLIDVecRef )
{
  AssertLIDs(staLIDVecRef.size() == numStateVars);

  // copy over the global ID lists.
  staLIDVec = staLIDVecRef;
}

//-----------------------------------------------------------------------------
// Function      : Instance::jacobianStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/18/2008
//-----------------------------------------------------------------------------
const std::vector< std::vector<int> > & Instance::jacobianStamp() const
{
  return jacStamp_;
}

//-----------------------------------------------------------------------------
// Function      : Instance::registerJacLIDs
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/18/2008
//-----------------------------------------------------------------------------
void Instance::registerJacLIDs( const std::vector< std::vector<int> > & jacLIDVec )
{
  DeviceInstance::registerJacLIDs( jacLIDVec );
  int numNodes=numExtVars+numIntVars;
  A_Equ_NodeOffsets_.resize(numNodes);
  for (int equ=0; equ < numNodes; ++equ)
  {
    A_Equ_NodeOffsets_[equ].resize(numNodes);
    for (int node=0; node < numNodes; ++node)
    {
      A_Equ_NodeOffsets_[equ][node] = jacLIDVec[equ][node];
    }
  }

}

//-----------------------------------------------------------------------------
// Function      : Instance::updatePrimaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/18/2008
//-----------------------------------------------------------------------------
bool Instance::updatePrimaryState()
{
  bool bsuccess = true;

#ifdef Xyce_DEBUG_DEVICE
  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << "Instance::updatePrimaryState" <<std::endl;
  }
#endif
  bsuccess = updateIntermediateVars();
  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::updateSecondaryState
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/18/2008
//-----------------------------------------------------------------------------
bool Instance::updateSecondaryState ()
{
  bool bsuccess = true;

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : Instance::updateIntermediateVars
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 9/12/2008
//-----------------------------------------------------------------------------
bool Instance::updateIntermediateVars ()
{

  bool bsuccess=true;
  int numNodes=numIntVars+numExtVars;
  N_LAS_Vector * solVectorPtr = extData.nextSolVectorPtr;


#ifdef Xyce_DEBUG_DEVICE

  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << subsection_divider << std::endl;
    Xyce::dout() << "Instance::updateIntermediateVars "<<std::endl;
    Xyce::dout() << "  name = " << getName() << std::endl;
  }
#endif

  if (getSolverState().newtonIter == 0) // start of new timestep
  {
    interpolateSandK_();
  }

  if (solutionVars.size() != numNodes)
  {
    // set up all the Fad vectors with the right sizes
    // We only need to do this the first time around
    solutionVars.resize(numNodes);
    fContributions.resize(numNodes);
    dV.resize(totalNumWindings);
    windingCurrents.resize(totalNumWindings);

    // Now, we have to loop over all those vectors and do "resize" operations
    // on every element in them --- this sets the number of derivatives
    // on the "DFad" type
    for (int node=0 ; node<numNodes ; ++node)
    {
      solutionVars[node].resize(numNodes);
      fContributions[node].resize(numNodes);
    }
    for (int winding=0 ; winding<totalNumWindings ; ++winding)
    {
      dV[winding].resize(numNodes);
      windingCurrents[winding].resize(numNodes);
    }
  }
  // initialize all the contributions to zero (automatically zeros derivs)
  for (int node=0; node<numNodes; ++node)
  {
    fContributions[node] = 0.0;
  }

  // same for the windings:
  for (int winding=0; winding<totalNumWindings; ++winding)
  {
    windingCurrents[winding] = 0.0;
  }


  // Extract solution variables and set as DFad independent variables:
  // (diff sets the row of the derivatives matrix to the identity row
  for (int node=0; node<numNodes; ++node)
  {
    solutionVars[node] = (*solVectorPtr)[li_Nodes_[node]];
    solutionVars[node].diff(node,numNodes);
#ifdef Xyce_DEBUG_DEVICE
    Xyce::dout() << "solutionVar[" << node << "] = " << solutionVars[node] << std::endl;
#endif
  }

  // Now compute winding dV's:
  for (int winding=0; winding<totalNumWindings; ++winding)
  {
    int posNode=windingNodes[winding].first;
    int negNode=windingNodes[winding].second;
    dV[winding] = solutionVars[posNode]-solutionVars[negNode];
  }

  // Now compute winding currents:
  // I_i=S_i + sum(j=0,totalNumWindings,K_ij*dV_j)
  for (int windingI=0; windingI<totalNumWindings; ++windingI)
  {
    // we might nave no source vector
    if (!theSourceVector_.empty())
      windingCurrents[windingI] = theSourceVector_[windingI];

    // we might nave no K matrix
    if (!theKMatrix_.empty())
    {
      for (int windingJ=0; windingJ<totalNumWindings; ++windingJ)
      {
        windingCurrents[windingI] += theKMatrix_[windingI][windingJ]*dV[windingJ];
      }
    }
  }

  // Now assemble the contributions for the nodes:
  for (int winding=0; winding<totalNumWindings; ++winding)
  {
    int posNode=windingNodes[winding].first;
    int negNode=windingNodes[winding].second;
#ifdef Xyce_DEBUG_DEVICE
    Xyce::dout() <<  "Winding " << winding << " adding " << windingCurrents[winding]
         << "to  contribution to node " << posNode << " and subtracting same "
         << " from node " << negNode << std::endl;
#endif
    fContributions[posNode] += windingCurrents[winding];
    fContributions[negNode] -= windingCurrents[winding];
  }

#ifdef Xyce_DEBUG_DEVICE
  Xyce::dout() << " Contributions for device " << getName() << std::endl;
  for (int node=0; node<numNodes; ++node)
  {
    Xyce::dout() << "   F[" << node << "] = " << fContributions[node] << std::endl;
  }
  Xyce::dout() << "Winding potential drops and currents: " << std::endl;
  for (int winding=0; winding<totalNumWindings; ++winding)
  {
    Xyce::dout() << "  dV[" << winding << "] = " << dV[winding] << std::endl;
    Xyce::dout() << "   I[" << winding << "] = " << windingCurrents[winding] << std::endl;
  }
#endif

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEQVector
//
// Purpose       : Loads the Q-vector contributions for a single
//                 instance.
//
// Special Notes : The "Q" vector is part of a standard DAE formalism in
//                 which the system of equations is represented as:
//
//                 f(x) = dQ(x)/dt + F(x) - B(t) = 0
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/18/2008
//-----------------------------------------------------------------------------
bool Instance::loadDAEQVector ()
{
  bool bsuccess = true;

  N_LAS_Vector * solVectorPtr = extData.nextSolVectorPtr;
  N_LAS_Vector * staVectorPtr = extData.nextStaVectorPtr;

  return bsuccess;
}



//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEFVector
//
// Purpose       : Loads the F-vector contributions for a single
//                 instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/18/2008
//-----------------------------------------------------------------------------
bool Instance::loadDAEFVector ()
{
  bool bsuccess=true;

  N_LAS_Vector * daeFVecPtr = extData.daeFVectorPtr;
  int numNodes = numExtVars+numIntVars;

#ifdef Xyce_DEBUG_DEVICE

  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << subsection_divider << std::endl;
    Xyce::dout() << "Instance::loadDAEFVector "<<std::endl;
    Xyce::dout() << "  name = " << getName() << std::endl;
  }
#endif

  for (int node=0; node<numNodes; ++node)
  {
    (*daeFVecPtr)[li_Nodes_[node]] += fContributions[node].val();
#ifdef Xyce_DEBUG_DEVICE
    Xyce::dout() << "  fVec[" << node << "] += "
         << fContributions[node].val() << std::endl;
    Xyce::dout() << "   loaded into local ID " << li_Nodes_[node] << std::endl;
#endif
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdQdx
//
// Purpose       : Loads the Q-vector contributions for a single
//                 instance.
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/18/2008
//-----------------------------------------------------------------------------
bool Instance::loadDAEdQdx ()
{
  bool bsuccess = true;

  //  N_LAS_Matrix * dQdxMatPtr = extData.dQdxMatrixPtr;

  return bsuccess;
}



//-----------------------------------------------------------------------------
// Function      : Instance::loadDAEdFdx ()
//
// Purpose       : Loads the F-vector contributions for a single
//                 instance.
//
// Special Notes :
//
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/18/2008
//-----------------------------------------------------------------------------
bool Instance::loadDAEdFdx ()
{
  bool bsuccess = true;

  N_LAS_Matrix * dFdxMatPtr = extData.dFdxMatrixPtr;
  N_LAS_Vector * staVectorPtr = extData.nextStaVectorPtr;
  N_LAS_Vector * solVectorPtr = extData.nextSolVectorPtr;

  int numNodes = numExtVars+numIntVars;
#ifdef Xyce_DEBUG_DEVICE

  if (getDeviceOptions().debugLevel > 0 && getSolverState().debugTimeFlag)
  {
    Xyce::dout() << subsection_divider << std::endl;
    Xyce::dout() << "Instance::loadDAEdFdx "<<std::endl;
    Xyce::dout() << "  name = " << getName() << std::endl;
  }
#endif

#if 0
  // Problem:  In orginal Xygra testing, I tried to have multi-winding,
  // source-only coils.  This lead to the problem where internal nodes were
  // completely floating, and lead to singular jacobians.  My initial reaction
  // was to populate the diagonals of the  jacobian for source-only
  // coils.  This is The Wrong Thing To Do.
  // It breaks the completely acceptable case of source-only single-winding
  // "coils" such as those needed by Emphasis, leading to nonlinear convergence
  // failures on what should be a completely linear problem.
  //
  // There should be something sane to do for multi-winding source-only
  // coils, but these are probably not going to be an important issue for now.
  // So punt it.
  if (theKMatrix_.empty())
  {
    for (int equ=0; equ < numNodes; ++equ)
    {
      (*dFdxMatPtr)[li_Nodes_[equ]][A_Equ_NodeOffsets_[equ][equ]] += 1.0;
    }
  }
  else
#endif
  {
    // Using Sacado differentiation, set the dFdX terms
    for (int equ=0; equ < numNodes; ++equ)
    {
      for (int node=0; node < numNodes; ++node)
      {
        (*dFdxMatPtr)[li_Nodes_[equ]][A_Equ_NodeOffsets_[equ][node]]
          += fContributions[equ].dx(node);
#ifdef Xyce_DEBUG_DEVICE
        Xyce::dout() << " dFdX[" << equ << "]["<<node<<"] += "
             << fContributions[equ].dx(node) << std::endl;
#endif
      }
    }
  }

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::setIC
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/18/2008
//-----------------------------------------------------------------------------
bool Instance::setIC ()
{
  bool bsuccess = true;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Instance::varTypes
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/18/2008
//-----------------------------------------------------------------------------
void Instance::varTypes( std::vector<char> & varTypeVec )
{
  //varTypeVec.resize(1);
  //varTypeVec[0] = 'I';
}

//-----------------------------------------------------------------------------
// Function      : Instance::constructComposite
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 9/11/2008
//-----------------------------------------------------------------------------
CompositeParam *Instance::constructComposite(const std::string & cName, const std::string & pName)
{
  if (cName == "COIL")
  {
    XygraCoilData *xcd = new XygraCoilData ();
    coilDataVec.push_back(xcd);
    return (static_cast<CompositeParam *> (xcd));
  }
  else
  {
    std::string msg =
      "Instance::constructComposite: unrecognized composite name: ";
    msg += cName;
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL,msg);
  }
  // never reached
  return NULL;
}

//-----------------------------------------------------------------------------
// Function      : Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/18/2008
//-----------------------------------------------------------------------------
bool Model::processParams ()
{
  return true;
}

//----------------------------------------------------------------------------
// Function      : Model::processInstanceParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/18/2008
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
// Function      : Model::Model
// Purpose       : block constructor
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/18/2008
//-----------------------------------------------------------------------------
Model::Model(
  const Configuration & configuration,
  const ModelBlock &    MB,
  const FactoryBlock &  factory_block)
  : DeviceModel(MB, configuration.getModelParameters(), factory_block)
{

  // // Set up mapping from param names to class variables:
  // if (parMap.empty())
  // {
  //   // Set up double precision variables:

  // }

  // Set params to constant default values:
  setDefaultParams ();

  // Set params according to .model line and constant defaults from metadata:
  setModParams (MB.params);

  // Set any non-constant parameter defaults:
  //if (!given("TNOM"))
  //  tnom = getDeviceOptions().tnom;

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
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/18/2008
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

//-----------------------------------------------------------------------------
// Function      : Model::printOutInstances
// Purpose       : debugging tool.
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/18/2008
//-----------------------------------------------------------------------------
std::ostream &Model::printOutInstances(std::ostream &os) const
{
  std::vector<Instance*>::const_iterator iter;
  std::vector<Instance*>::const_iterator first = instanceContainer.begin();
  std::vector<Instance*>::const_iterator last  = instanceContainer.end();

  int i, isize;
  isize = instanceContainer.size();

  os << std::endl;
  os << "Number of Xygra instances: " << isize << std::endl;
  os << "    name=\t\tmodelName\tParameters" << std::endl;
  for (i=0, iter=first; iter!=last; ++iter, ++i)
  {
    os << "  " << i << ": " << (*iter)->getName() << "\t";
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


Device *Traits::factory(const Configuration &configuration, const FactoryBlock &factory_block)
{

  return new DeviceMaster<Traits>(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);
}

void registerDevice()
{
  Config<Traits>::addConfiguration()
    .registerDevice("xygra", 1)
    .registerModelType("xygra", 1);
}

} // namespace Xygra
} // namespace Device
} // namespace Xyce

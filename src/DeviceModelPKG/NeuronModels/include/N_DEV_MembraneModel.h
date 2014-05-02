//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2011  Sandia Corporation
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
// Filename       : $RCSfile: N_DEV_MembraneModel.h,v $
//
// Purpose        : Neuron classes.
//
// Special Notes  :
//
// Creator        : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date  : 08/11/10
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.10.2.1 $
//
// Revision Date  : $Date: 2014/02/26 20:16:30 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_MembraneModel_h
#define Xyce_N_DEV_MembraneModel_h

#include <vector>

// ----------   Xyce Includes   ----------
#include <N_DEV_fwd.h>


// ---------- Forward Declarations ----------
class N_LAS_Vector;
class N_LAS_Matrix;

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : N_DEV_MembraneModel
// Purpose       : This is class is a virtual base class to define the base
//                 interface for membrane modeling equations.  Classes
//                 that derive from this will define the equations that
//                 need to be solved for various ion currents
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 08/11/2010
//-----------------------------------------------------------------------------
class MembraneModel
{
public:
  MembraneModel(const SolverState & ss1)
    : solState(ss1),
      numIndependentVars_(0),
      numExternalVars_(2)
  {}
    
  ~MembraneModel() {}

  int numVars() { return numIndependentVars_; }

  virtual void setJacStamp( int numExtVars, int segmentNumber, int vOffset, std::vector< std::vector< int > > & segmentJacStamp ) {}
  virtual void loadDAEQVector( int segmentNumber, std::vector< int > & lidIndexVector, N_LAS_Vector * solnVecPtr, N_LAS_Vector * daeQVecPtr, double segArea) {}
  virtual void loadDAEFVector( int segmentNumber, std::vector< int > & lidIndexVector, N_LAS_Vector * solnVecPtr, N_LAS_Vector * daeFVecPtr, double segArea) {}
  virtual void loadDAEdQdx( int segmentNumber, int vOffset, std::vector< int > & lidIndexVector, std::vector< std::vector< int > > & jacobianOffsets, N_LAS_Vector * solnVecPtr, N_LAS_Matrix * dQdxMatPtr, double segArea) {}
  virtual void loadDAEdFdx( int segmentNumber, int vOffset, std::vector< int > & lidIndexVector, std::vector< std::vector< int > > & jacobianOffsets, N_LAS_Vector * solnVecPtr, N_LAS_Matrix * dFdxMatPtr, double segArea) {}

  int numIndependentVars_;

  // these are modeling parameters common to most model types
  // this is some data duplication between the devcie model the membrane model
  // need to think of a cleaner way to do this
  double cMem_;
  double gMem_;
  double vRest_;

  const int numExternalVars_;  // always assume that the owning cable equation has two external vars (a Vin and Vout)
  // if this assumption needs to be changed we'll have this constant to show where we made such an assumption.

  const SolverState & solState;  // this is here incase a model needs to be aware of simulator data like dcopflag or currtime.
};

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::MembraneModel N_DEV_MembraneModel;

#endif


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
// Filename       : $RCSfile: N_DEV_MembranePassive.h,v $
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

#ifndef Xyce_N_DEV_MembranePassive_h
#define Xyce_N_DEV_MembranePassive_h


// ----------   Xyce Includes   ----------
#include <N_DEV_MembraneModel.h>


// ---------- Forward Declarations ----------

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : MembranePassive
// Purpose       : This is class defines a passive ion channel membrane with
//                 a membrane capacitance and resistance
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 08/11/2010
//-----------------------------------------------------------------------------
class MembranePassive : public MembraneModel
{
public:
  MembranePassive(const SolverState & ss1, double cMem, double gMem, double vRest);
  ~MembranePassive() {}

  void setJacStamp( int numExtVars, int segmentNumber, int vOffset, std::vector< std::vector< int > > & segmentJacStamp );
  void loadDAEQVector( int segmentNumber, std::vector< int > & lidIndexVector, N_LAS_Vector * solnVecPtr, N_LAS_Vector * daeQVecPtr, double segArea);
  void loadDAEFVector( int segmentNumber, std::vector< int > & lidIndexVector, N_LAS_Vector * solnVecPtr, N_LAS_Vector * daeFVecPtr, double segArea);
  void loadDAEdQdx( int segmentNumber, int vOffset, std::vector< int > & lidIndexVector, std::vector< std::vector< int > > & jacobianOffsets, N_LAS_Vector * solnVecPtr, N_LAS_Matrix * dQdxMatPtr, double segArea);
  void loadDAEdFdx( int segmentNumber, int vOffset, std::vector< int > & lidIndexVector, std::vector< std::vector< int > > & jacobianOffsets, N_LAS_Vector * solnVecPtr, N_LAS_Matrix * dFdxMatPtr, double segArea);

  // constitutive parameters
  double cMem_;     // membrane capacitance
  double gMem_;     // membrane conductance
  double vRest_;    // membrane rest voltage
};

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::MembranePassive N_DEV_MembranePassive;

#endif

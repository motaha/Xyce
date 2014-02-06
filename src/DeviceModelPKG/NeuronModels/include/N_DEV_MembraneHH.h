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
// Filename       : $RCSfile: N_DEV_MembraneHH.h,v $
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
// Revision Number: $Revision: 1.7.2.1 $
//
// Revision Date  : $Date: 2013/10/03 17:23:33 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_MembraneHH_h
#define Xyce_N_DEV_MembraneHH_h

#include <N_DEV_MembraneModel.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : MembraneHH
// Purpose       : This is class defines a membrane with Hodgkin-Huxley equations
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 08/11/2010
//-----------------------------------------------------------------------------
class MembraneHH : public MembraneModel
{
  public:
    MembraneHH(SolverState & ss1, double cMem, double gMem, double vRest, double eK, double gK, double eNa, double gNa);
    ~MembraneHH() {}

    void setJacStamp( int numExtVars, int segmentNumber, int vOffset, vector< vector< int > > & segmentJacStamp );
    void loadDAEQVector( int segmentNumber, vector< int > & lidIndexVector, N_LAS_Vector * solnVecPtr, N_LAS_Vector * daeQVecPtr, double segArea);
    void loadDAEFVector( int segmentNumber, vector< int > & lidIndexVector, N_LAS_Vector * solnVecPtr, N_LAS_Vector * daeFVecPtr, double segArea);
    void loadDAEdQdx( int segmentNumber, int vOffset, vector< int > & lidIndexVector, vector< vector< int > > & jacobianOffsets, N_LAS_Vector * solnVecPtr, N_LAS_Matrix * dQdxMatPtr, double segArea);
    void loadDAEdFdx( int segmentNumber, int vOffset, vector< int > & lidIndexVector, vector< vector< int > > & jacobianOffsets, N_LAS_Vector * solnVecPtr, N_LAS_Matrix * dFdxMatPtr, double segArea);

    double cMem_;     // membrane capacitance
    double gMem_;     // membrane conductance
    double vRest_;    // membrane rest voltage
    double eK_;
    double gK_;
    double eNa_;
    double gNa_;
};

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::MembraneHH N_DEV_MembraneHH;

#endif

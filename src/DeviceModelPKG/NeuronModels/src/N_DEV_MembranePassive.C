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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_MembranePassive.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Richard Schiek, Electrical and Microsytem Modeling
//
// Creation Date  : 08/11/2010
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.12 $
//
// Revision Date  : $Date: 2014/01/23 16:19:04 $
//
// Current Owner  : $Author: dgbaur $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include<iostream>

// ---------- Standard Includes ----------

#include <N_UTL_Misc.h>

// ----------   Xyce Includes   ----------
#include <N_DEV_MembranePassive.h>
#include <N_LAS_Vector.h>
#include <N_LAS_Matrix.h>
#include <N_DEV_SolverState.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : MembranePassive::MembranePassive
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 08/11/2010
//-----------------------------------------------------------------------------
MembranePassive::MembranePassive(const SolverState & ss1, double cMem, double gMem, double vRest)
  : MembraneModel(ss1),
    gMem_(gMem),
    cMem_(cMem),
    vRest_(vRest)
{
  // passive membrane just has voltage as its unknown variable
  // so set up numIndependentVars_ for that
  numIndependentVars_ = 1;
}

//-----------------------------------------------------------------------------
// Function      : MembranePassive::setJacStamp
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 08/11/2010
//-----------------------------------------------------------------------------
void MembranePassive::setJacStamp( int numExtVars, int segmentNumber, int vOffset, std::vector< std::vector< int > > & segmentJacStamp )
{
  // In a passive cable the membrane is just two passive elements, a capacitor and a resistor
  // thus the membrane current, I = f(Vsegment).

  /*
  int offset = numExtVars + numIndependentVars_*segmentNumber;

  int jacobianRowSize = segmentJacStamp[offset].size();
  */
  /* need to handle these in a better way
     currently they can confuse parameter testing when dummy devices are created.
  if( jacobianRowSize == 1 )
  {
    if( segmentJacStamp[ offset ][0] != offset )
    {
      Xyce::dout() << "Potential error in MembranePassive::setJacStamp().  segmentJacStamp[ " << offset << " ][0] != " << offset << std::endl;
    }
  }
  else if( jacobianRowSize > 1 )
  {
    if( segmentJacStamp[ offset ][1] != offset )
    {
      Xyce::dout() << "Potential error in MembranePassive::setJacStamp().  segmentJacStamp[ " << offset << " ][1] != " << offset << std::endl;
    }
  }
  */

}

//-----------------------------------------------------------------------------
// Function      : MembranePassive::loadDAEQVector
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 08/11/2010
//-----------------------------------------------------------------------------
void MembranePassive::loadDAEQVector( int segmentNumber, std::vector< int > & lidIndexVector,  N_LAS_Vector * solnVecPtr, N_LAS_Vector * daeQVecPtr, double segArea)
{
  // Each segment will have numIndependentVars_ with segment voltage being the first
  // so, the cMem dV/dt term will be at segmentNumber * numIndependentVars_.
  // in the case of the passive cable numIndependentVars_=1.
  (*daeQVecPtr)[lidIndexVector[segmentNumber]] += cMem_ * segArea * (*solnVecPtr)[lidIndexVector[segmentNumber]];
}

//-----------------------------------------------------------------------------
// Function      : MembranePassive::loadDAEFVector
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 08/11/2010
//-----------------------------------------------------------------------------
void MembranePassive::loadDAEFVector( int segmentNumber, std::vector< int > & lidIndexVector, N_LAS_Vector * solnVecPtr, N_LAS_Vector * daeFVecPtr, double segArea)
{
  // Each segment will have numIndependentVars_ with segment voltage being the first
  // so, the cMem dV/dt term will be at segmentNumber * numIndependentVars_.
  // in the case of the passive cable numIndependentVars_=1.
  (*daeFVecPtr)[lidIndexVector[segmentNumber]] += gMem_ * segArea * ((*solnVecPtr)[lidIndexVector[segmentNumber]] - vRest_ );
}

//-----------------------------------------------------------------------------
// Function      : MembranePassive::loadDAEdQdx
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 08/11/2010
//-----------------------------------------------------------------------------
void MembranePassive::loadDAEdQdx( int segmentNumber, int vOffset, std::vector< int > & lidIndexVector, std::vector< std::vector< int > > & jacobianOffsets, N_LAS_Vector * solnVecPtr, N_LAS_Matrix * dQdxMatPtr, double segArea)
{
   // while lidIndexVector lists LID's for just the segment variables (just V in the case
   // of a passive cable). The jacobianOffsets includes the Vin and Vout as the first
   // two variables.  Thus, there is a constant offset of 2 for everything in jacobianOffsets

   // And, as in the Q and F load functions,  Each segment will have numIndependentVars_ with segment voltage being the first
   // so, the cMem dV/dt term will be at segmentNumber * numIndependentVars_.
   // in the case of the passive cable numIndependentVars_=1.

   int row = numExternalVars_ + segmentNumber;       // numExternalVars_ a contant of 2 assumed in MembraneModel base class

   (*dQdxMatPtr)[lidIndexVector[segmentNumber]][jacobianOffsets[row][vOffset]] += cMem_ * segArea;
}

//-----------------------------------------------------------------------------
// Function      : MembranePassive::loadDAEdFdx
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsytem Modeling
// Creation Date : 08/11/2010
//-----------------------------------------------------------------------------
void MembranePassive::loadDAEdFdx( int segmentNumber, int vOffset, std::vector< int > & lidIndexVector, std::vector< std::vector< int > > & jacobianOffsets, N_LAS_Vector * solnVecPtr, N_LAS_Matrix * dFdxMatPtr, double segArea)
{
   // while lidIndexVector lists LID's for just the segment variables (just V in the case
   // of a passive cable). The jacobianOffsets includes the Vin and Vout as the first
   // two variables.  Thus, there is a constant offset of 2 for everything in jacobianOffsets

   // And, as in the Q and F load functions,  Each segment will have numIndependentVars_ with segment voltage being the first
   // so, the cMem dV/dt term will be at segmentNumber * numIndependentVars_.
   // in the case of the passive cable numIndependentVars_=1.

   int row = numExternalVars_ + segmentNumber;       // numExternalVars_ a contant of 2 assumed in MembraneModel base class

   (*dFdxMatPtr)[lidIndexVector[segmentNumber]][jacobianOffsets[row][vOffset]] += gMem_ * segArea;

}

} // namespace Device
} // namespace Xyce

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

//-----------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_bcData.h,v $
//
// Purpose        : This file contains the classes neccessary for a PDE
//                  based diode simulation.
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
// Revision Number: $Revision: 1.9.6.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:36 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_bcData_h
#define Xyce_N_DEV_bcData_h

// ----------   Xyce Includes   ----------

// ---------- Forward Declarations -------
#include <iosfwd>

//-----------------------------------------------------------------------------
// Class         : bcData
//
// Purpose       : Boundary condition info for a single boundary.
//
// Special Notes : Kind of like the deviceInterfaceNode class from the 2D
//                 device.
//
// Creator       : Eric Keiter
// Creation Date : 04/01/03
//-----------------------------------------------------------------------------
class bcData
{
  private:
  protected:
  public:
    bcData ():
     eName(""),
     nName(""),
     type("ntype"),
     index(-1),
     Vequ(0.0),
     VequGiven(false),
     Vckt(0.0), 
     Vckt_orig(0.0),
     gid(-1), 
     lid(-1),          // local ID of the circuit node.
     lidOffset(-1),    // diagonal col in the ckt node matrix row. DMA only.
     lidlidPtr(0),     // diagonal matrix pointer.
     colArray(),
     li_colArray(),
     crossOffsets(),
     given(false),
     area(1.0),
     areaGiven(false),
     currentSum(0.0), elecCurrent(0.0), holeCurrent(0.0),
     stateC(-1), stateC_owned(true), 
     li_stateC(-1),
     meshIndex(0),
     neighborNode(1),
     dIdVckt(0.0),
     dFdVckt(),
     dIdX(),
     dIdXcols(),
     dIdXoffset(),
     numCrossTerms(0),
     Vbc(0.0),
     nnbc(0.0),
     npbc(0.0),
     dxdvPtr(0),
     dxdvAllocated(false),
     Vckt_old(0.0), Vckt_final(0.0), 
     Vckt_delta(0.0),
     Vckt_deltaC(0.0),
     Vckt_ramp(0.0), Vckt_ramp_old(0.0),
     displCurrent(0.0),
     material ("neutral"),
     materialGiven(false),
     oxideBndryFlag(false)
      {};

  private:
  protected:
  public:
    // electrode name
    string eName;     // device electrode name.  Should be in all upper case.
    string nName;     // circuit node name.
    string type;      // string to indicate ntype or ptype.

    int    index;     // index w.r.t. the order it was listed in the netlist.

    // Equilibrium voltage
    double Vequ;
    bool   VequGiven;

    // Circuit node voltage
    double Vckt;

    // Vckt, before voltage limiting.
    double Vckt_orig;

    int    gid;       // global ID of the circuit node.  This is the global
                      // solution vector index.  (used like a Vrowarray entry)

    int lid;            // local ID of the circuit node.
    int lidOffset;      // diagonal col in the ckt node matrix row. DMA only.
    double * lidlidPtr; // lid matrix diagonal pointer.

    vector<int> colArray;    // matrix col array for the KCL equation.
    vector<int> li_colArray; // matrix col array for the KCL equation, dma
    
    vector<int> crossOffsets; // columns for the other ckt nodes.

    bool   given;     // did the user specify this node or not.  If not, the
                      // device electrode with be (probably) connected to gnd.

    // geometrical stuff:
    double area;                // total area for the edge.
    bool   areaGiven;           // given boolean.

    // state variable stuff:
    double currentSum;    // sum of currents, to be stored as a state variable.
    double elecCurrent;   // electron current
    double holeCurrent;   // hole     current

    int    stateC;        // state variable index, current.
    bool   stateC_owned;

    //local id's (offsets)
    int li_stateC;

    // mesh index of this electrode:
    int meshIndex;
    int neighborNode;  // node neighboring this boundary.

    // information needed for the 2-level Newton.
    double dIdVckt; // derivative of currentSum w.r.t. Vckt.

    vector<double> dFdVckt;         // deriv. of residual w.r.t. Vckt.

    vector<double> dIdX;     // deriv. of currentSum w.r.t. PDE solution vars.
    vector<int>    dIdXcols; // nodes neighboring this boundary.
                             // (this may be the same as colArray, in which
			     // case one of them will be removed...)

    vector<int>    dIdXoffset; // If running with DMA, use this instead of
                               // dIdXcols.

    int numCrossTerms;

    // boundary conditions to be imposed on V,n and p.
    double Vbc;
    double nnbc;
    double npbc;

    // dxdv vector:
    N_LAS_Vector * dxdvPtr;
    bool dxdvAllocated;

    // These BC variables are only used for Continuation NL solves.
    // first 3 are w.r.t change in Vckt between ckt iterations.  Big change
    double Vckt_old;    // ckt value from the previous Newton solve.
    double Vckt_final;  // eventual ckt voltage value for the current solve.
    double Vckt_delta;  // total change in the Vckt. (Vckt_final-Vckt_old)

    // next 3, w.r.t change between continuation solves. small change
    double Vckt_deltaC;  // incremental, intermediate change in the Vckt.
    double Vckt_ramp;   // current intermediate value of Vckt used during
                        // continuation.
    double Vckt_ramp_old;   // old intermediate value of Vckt used during
                            // continuation.

    double displCurrent;

    // material information:
    string material;
    bool   materialGiven;
    bool   oxideBndryFlag;

    vector<int> volIndices;

};
// inline functions
//-----------------------------------------------------------------------------
// Function      : bcData::operator<<
// Purpose       : "<<" operator
// Special Notes :
// Scope         : public
// Creator       : Eric R. Keiter, 9233, SNL, Parallel Computational Sciences
// Creation Date : 11/04/10
//-----------------------------------------------------------------------------

inline ostream & operator<<(ostream & os, const bcData & bc)
{
  os << "------------------------------------------\n";
  os << "electrode name = " << bc.eName << "\n";
  os << "node      name = " << bc.nName << "\n";
  os << "         index = " << bc.index << "\n";

  os << " Vequ = " << bc.Vequ << "\n";
  if (bc.VequGiven)
  {
  os << " VequGiven = true\n";
  }
  else
  {
  os << " VequGiven = false\n";
  }

  os << " Vckt = " << bc.Vckt << "\n";
  os << " Vckt_orig = " << bc.Vckt_orig << "\n";
  os << " gid = " << bc.gid << "\n";
  os << " lid = " << bc.lid << "\n";
  os << " lidOffset = " << bc.lidOffset << "\n";
  os << " given = " << bc.given << "\n";
  os << " area = " << bc.area << "\n";

  if (bc.areaGiven)
  {
  os << " areaGiven = true\n";
  }
  else
  {
  os << " areaGiven = false\n";
  }

  os << " currentSum = " << bc.currentSum << "\n";
  os << " elecCurrent = " << bc.elecCurrent << "\n";
  os << " holeCurrent = " << bc.holeCurrent << "\n";

  os <<  " stateC = " << bc.stateC << "\n";
  os <<  " stateC_owned = " << bc.stateC_owned << "\n";
  os <<  " li_stateC = " << bc.li_stateC << "\n";
  os <<  " meshIndex = " << bc.meshIndex << "\n";
  os <<  " neighborNode = " << bc.neighborNode << "\n";

  os <<  " dIdVckt = " << bc.dIdVckt << "\n";

  os << " numCrossTerms = " << bc.numCrossTerms << "\n";

  os << " Vbc = " << bc.Vbc << "\n";
  os << " nnbc = " << bc.nnbc << "\n";
  os << " npbc = " << bc.npbc << "\n";


  os << "material = " << bc.material << "\n";

  if (bc.materialGiven)
  {
  os << "materialGiven = true\n";
  }
  else
  {
  os << "materialGiven = false\n";
  }
  if (bc.oxideBndryFlag)
  {
  os << "oxideBndrFlag = true\n";
  }
  else
  {
  os << "oxideBndrFlag = false\n";
  }

  os << " Vckt_old = " << bc.Vckt_old << "\n";
  os << " Vckt_final = " << bc.Vckt_final << "\n";
  os << " Vckt_delta = " << bc.Vckt_delta << "\n";
  os << " Vckt_deltaC = " << bc.Vckt_deltaC << "\n";
  os << " Vckt_ramp = " << bc.Vckt_ramp << "\n";
  os << " Vckt_ramp_old = " << bc.Vckt_ramp_old << "\n";
  os << " displCurrent = " << bc.displCurrent << "\n";

  //N_LAS_Vector * dxdvPtr;
  if (bc.dxdvAllocated)
  {
    os << " dxdvAllocated = true\n";
  }
  else
  {
    os << " dxdvAllocated = false\n";
  }

  int i=0;
  for (i=0;i<bc.dFdVckt.size();++i)
  {
    cout << " dFdVckt["<<i<<"] = " << bc.dFdVckt[i] << "\n";
  }
  for (i=0;i<bc.dIdX.size();++i)
  {
    cout << " dIdX["<<i<<"] = " << bc.dIdX[i] << "\n";
  }
  for (i=0;i<bc.dIdXcols.size();++i)
  {
    cout << " dIdXcols["<<i<<"] = " << bc.dIdXcols[i] << "\n";
  }
  for (i=0;i<bc.dIdXoffset.size();++i)
  {
    cout << " dIdXoffset["<<i<<"] = " << bc.dIdXoffset[i] << "\n";
  }
  for (i=0;i<bc.colArray.size();++i)
  {
    cout << " collArray["<<i<<"] = " << bc.colArray[i] << "\n";
  }
  for (i=0;i<bc.li_colArray.size();++i)
  {
    cout << " li_collArray["<<i<<"] = " << bc.li_colArray[i] << "\n";
  }
  for (i=0;i<bc.crossOffsets.size();++i)
  {
    cout << " crossOffsets["<<i<<"] = " << bc.crossOffsets[i] << "\n";
  }

  os << "------------------------------------------\n";
  os << endl;

  return os;
}

#endif
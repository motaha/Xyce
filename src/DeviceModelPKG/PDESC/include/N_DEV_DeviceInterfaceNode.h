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
// Filename       : $RCSfile: N_DEV_DeviceInterfaceNode.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 07/08/03
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.10.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:35 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_DEVICE_INTERFACE_NODE_h
#define Xyce_N_DEV_DEVICE_INTERFACE_NODE_h

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : DeviceInterfaceNode
// Purpose       : This class contains information about a single circuit node
//                 which is connected to a single device electrode.
// Special Notes :
// Creator       : Eric Keiter
// Creation Date : 04/23/02
//-----------------------------------------------------------------------------
class DeviceInterfaceNode
{
  private:
  protected:
  public:
    DeviceInterfaceNode () :
     eName(""), nName(""),index(-1), gid(-1), given(false),
     firstMeshNodeIndex(-1),
     stateC(-1), stateC_owned(true), currentSum(0.0), area(0.0),
     chargeSum(0.0),
     Vckt(0.0), Vckt_ramp(0.0), Vckt_ramp_old(0.0),
     Vckt_old(0.0), Vckt_final(0.0), Vckt_deltaC(0.0), Vckt_delta(0.0),
     Vckt_orig(0.0),
     dIdVckt(0.0),
     dQdVckt(0.0),
     neumannBCFlagV(false),
     neumannBCFlagN(false),
     neumannBCFlagP(false),
     numBoundaryPoints(0),
     dxdvPtr(NULL),
     dxdvAllocated(false),
     numCrossTerms(0),
     material ("neutral"),
     materialGiven(false),
     oxideBndryFlag(false),
     oxthick(0.0),
     oxcharge(0.0)
     {
       Vcol.reserve(20);
       Ncol.reserve(20);
       Pcol.reserve(20);
       areaVector.reserve(20);
     };

  private:
  protected:
  public:
    string eName;     // device electrode name.  Should be in all upper case.
    string nName;     // circuit node name.

    int    index;     // index w.r.t. the order it was listed in the netlist.

    int   labelIndex; // index into the label array.

    int    gid;       // global ID of the circuit node.  This is the global
                      // solution vector index.  (used like a Vrowarray entry)

    int lid;          // local ID of the circuit node.
    int lidOffset;    // diagonal col in the ckt node matrix row. DMA only.

    vector<int> crossOffsets; // columns for the other ckt nodes.

    bool   given;     // did the user specify this node or not.  If not, the
                      // device electrode with be (probably) connected to gnd.


    int firstMeshNodeIndex;
                      // the mesh index for a single point on the
                      // mesh that is part of this electrode.

    // The "gid" Vcol, Ncol, and Pcol are only really used during the setup
    // phase, when setting up (row, col) pairs.
    // As such, they might not be needed.
    vector<int> Vcol;  // column array
    vector<int> Ncol;  // column array
    vector<int> Pcol;  // column array

    // geometrical stuff:
    double area;                // total area for the edge.
    vector<double> areaVector;  // area for each edge node

    // state variable stuff:
    double currentSum;    // sum of currents, to be stored as a state variable.
    int    stateC;        // state variable index, current.
    bool   stateC_owned;

    // capacitance related stuff
    double chargeSum;     // total charge on the electrode.

    //local id's (offsets)
    int li_stateC;

    // number of mesh points for this boundary.
    int numBoundaryPoints;

    // from the circuit node:
    double Vckt;

    // Vckt, before voltage limiting.
    double Vckt_orig;

    // derivative information needed for the 2-level Newton.
    double dIdVckt; // derivative of currentSum w.r.t. Vckt.
    double dQdVckt; // derivative of chargeSum w.r.t. Vckt.

    vector<double> dFdVckt;         // deriv. of residual w.r.t. Vckt.
    vector<int>    neighborNodes;   // nodes neighboring this boundary.

    vector<double> dQdX;     // deriv. of chargeSum w.r.t. PDE solution vars.
    vector<double> dIdX;     // deriv. of currentSum w.r.t. PDE solution vars.
    vector<int>    dIdXcols; // nodes neighboring this boundary.

    vector<int>    dIdXoffset; // If running with DMA, use this instead of
                               // dIdXcols.

    int numCrossTerms;

    // equilibrium voltage at this electrode (if the applied voltage is
    // zero, there is still an internal voltage drop).
    vector<double> VequVec;

    // boundary conditions to be imposed on V,n and p.
    vector<double> VbcVec;
    vector<double> nnbcVec;
    vector<double> npbcVec;

    // dxdv vector:
    N_LAS_Vector * dxdvPtr;
    bool dxdvAllocated;

    // this map is between mesh node ID's and the indices of VbcVec
    map<int,int> meshGlobalToLocal;

    // These BC variables are only used for Continuation NL solves.
    // first 3 are w.r.t change in Vckt between ckt iterations.  Big change
    double Vckt_old;    // ckt value from the previous Newton solve.
    double Vckt_final;  // eventual ckt voltage value for the current solve.
    double Vckt_delta;  // total change in the Vckt. (Vckt_final-Vckt_old)

    // next 3, w.r.t change between continuation solves. small change
    double Vckt_deltaC; // incremental, intermediate change in the Vckt.
    double Vckt_ramp;   // current intermediate value of Vckt used during
                        // continuation.
    double Vckt_ramp_old;   // old intermediate value of Vckt used during
                            // continuation.

    // neumann BC flag;
    bool neumannBCFlagV;
    bool neumannBCFlagN;
    bool neumannBCFlagP;

    // material information:
    string material;
    bool   materialGiven;
    bool   oxideBndryFlag;
    double oxthick;
    double oxcharge;
};


} // namespace Device
} // namespace Xyce

typedef Xyce::Device::DeviceInterfaceNode N_DEV_DeviceInterfaceNode;

#endif // Xyce_N_DEV_DEVICE_INTERFACE_NODE_h



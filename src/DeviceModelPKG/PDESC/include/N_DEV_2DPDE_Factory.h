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
// Filename       : $RCSfile: N_DEV_2DPDE_Factory.h,v $
//
// Purpose        : This file contains the classes neccessary for a 2D PDE
//                  based simulation.  MOSFETs, BJTs, Diodes, etc.
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 11/14/01
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.4.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:35 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_2DPDE_Factory_h
#define Xyce_N_DEV_2DPDE_Factory_h

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceTemplate.h>
#include <N_DEV_2DPDE.h>

//-----------------------------------------------------------------------------
// Class         : N_DEV_2DPDE
// Purpose       : Handles 2-dimensional PDE devices.
//
// Special Notes :
//
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 2/03/06
//-----------------------------------------------------------------------------
struct N_DEV_2DPDE
{
  //---------------------------------------------------------------------------
  // Function      : factory
  // Purpose       : This is the factory function for the 2DPDE class.
  //
  // Special Notes : ERK.  10/16/2005.  This used to be a singleton (ie a
  //                 static pointer was returned) but had to be changed
  //                 so that the library version of Xyce would work
  //                 correctly.
  //
  // Scope         : public
  // Creator       : Eric Keiter, SNL, Parallel Computational Sciences
  // Creation Date : 3/16/00
  //---------------------------------------------------------------------------
  static N_DEV_Device * factory(N_DEV_SolverState & ss1,
                                N_DEV_DeviceOptions & do1)
  {
     const string name("2D PDE Device");
     const string className("N_DEV_2DPDE");
     const string defaultModelName ("PDE level 2");

     // no required parameters:
     vector< pair<string,double> > nameVec;


     // Set up a "default" instance block:
     N_DEV_Param p ("NX", 11, true);
     N_DEV_InstanceBlock IB(defaultModelName);
     IB.params.push_back(p);
     p.set("NODE", "DEFAULT");
     p.setGiven(true);
     IB.params.push_back(p);
     p.set("NODE0.NAME","COLLECTOR");
     IB.params.push_back(p);
     p.set("NODE1.NAME","BASE");
     IB.params.push_back(p);
     p.set("NODE2.NAME","EMITTER");
     IB.params.push_back(p);

     p.set("REGION", "DEFAULT");
     IB.params.push_back(p);
     p.set("REGION0.NAME", "DEFAULT_0");
     IB.params.push_back(p);
     p.set("REGION1.NAME", "DEFAULT_1");
     IB.params.push_back(p);
     p.set("REGION0.FUNCTION", "STEP");
     IB.params.push_back(p);
     p.set("REGION1.FUNCTION", "STEP");
     IB.params.push_back(p);
     p.set("REGION0.TYPE", "NTYPE");
     IB.params.push_back(p);
     p.set("REGION1.TYPE", "PTYPE");
     IB.params.push_back(p);
     p.set("REGION0.XLOC",0.0);
     IB.params.push_back(p);
     p.set("REGION1.XLOC",2.5e-4);
     IB.params.push_back(p);
     p.set("REGION0.FLATX",-1);
     IB.params.push_back(p);
     p.set("REGION1.FLATX",1);
     IB.params.push_back(p);


     N_DEV_Device * devptr =
       new Xyce::Device::DeviceTemplate<N_DEV_2DPDEModel, N_DEV_2DPDEInstance>
                   ( name,
                     className,
                     defaultModelName,
                     IB,
                     false,  // not linear
                     true,   // requires a model
                     ss1,do1);

     return devptr;
  }
};

#endif

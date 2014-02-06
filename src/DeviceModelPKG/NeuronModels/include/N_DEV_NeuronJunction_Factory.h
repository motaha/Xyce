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
// Filename       : $RCSfile: N_DEV_NeuronJunction_Factory.h,v $
//
// Purpose        : NeuronJunction classes
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
// Revision Number: $Revision: 1.3.2.1 $
//
// Revision Date  : $Date: 2013/10/03 17:23:33 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_NeuronJunctionFactory_h
#define Xyce_N_DEV_NeuronJunctionFactory_h

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceTemplate.h>
#include <N_DEV_NeuronJunction.h>

//-----------------------------------------------------------------------------
// Class         : N_DEV_NeuronJunction
// Purpose       : Handles NeuronJunctions
//
// Special Notes :
//
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 2/03/06
//-----------------------------------------------------------------------------
struct N_DEV_NeuronJunction
{
  //---------------------------------------------------------------------------
  // Function      : factory
  // Purpose       : This is the factory function for the NeuronJunction class.
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
     const string name("NeuronJunction");
     const string className("N_DEV_NeuronJunction");
     const string defaultModelName("YNEURONJUNCTION level 1");

     vector< pair<string,double> > nameVec;

     N_DEV_Device * devptr =
       new N_DEV_NeuronJunctionMaster
       //new Xyce::Device::DeviceTemplate<N_DEV_NeuronJunctionModel, N_DEV_NeuronJunctionInstance>
                   ( name,
                     className,
                     defaultModelName,
                     nameVec,
                     true,    // linear
                     false,   // model not required
                     ss1,do1);

     return devptr;
  }
};

#endif


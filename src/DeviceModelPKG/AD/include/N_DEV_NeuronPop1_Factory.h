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
// Filename       : $RCSfile: N_DEV_NeuronPop1_Factory.h,v $
//
// Purpose        : Neuron classes.
//
// Special Notes  :
//
// Creator        : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date  : 06/10/09
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.2.2.1 $
//
// Revision Date  : $Date: 2013/10/03 17:23:33 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_NeuronPop1_Factory_h
#define Xyce_N_DEV_NeuronPop1_Factory_h

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceTemplate.h>
#include <N_DEV_NeuronPop1.h>

//-----------------------------------------------------------------------------
// Class         : N_DEV_NeuronPop1
// Purpose       : Handles Neuron populations
//
// Special Notes :
//
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 02/23/11
//-----------------------------------------------------------------------------
struct N_DEV_NeuronPop1
{
  //---------------------------------------------------------------------------
  // Function      : factory
  // Purpose       : This is the factory function for the Neuron population class.
  //
  // Special Notes : ERK.  10/16/2005.  This used to be a singleton (ie a
  //                 static pointer was returned) but had to be changed
  //                 so that the library version of Xyce would work
  //                 correctly.
  //
  // Scope         : public
  // Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
  // Creation Date : 02/23/11
  //---------------------------------------------------------------------------
  static N_DEV_Device * factory(N_DEV_SolverState & ss1,
                                N_DEV_DeviceOptions & do1)
  {
     const string name("NeuronPopulation");
     const string className("N_DEV_NeuronPop1");
     const string defaultModelName("YNEURONPOP level 1");

     //pair <string,double> p(string("Neuron"), 1.0);
     vector< pair<string,double> > nameVec; //(1, p);

     N_DEV_Device * devptr =
       new Xyce::Device::DeviceTemplate<N_DEV_NeuronPopModel1, N_DEV_NeuronPopInstance1>
                   ( name,
                     className,
                     defaultModelName,
                     nameVec,
                     true,   // linear
                     true,   // model required
                     ss1,do1);

     return devptr;
  }
};

#endif

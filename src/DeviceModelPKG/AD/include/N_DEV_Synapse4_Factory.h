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
// Filename       : $RCSfile: N_DEV_Synapse4_Factory.h,v $
//
// Purpose        : Synapse4 classes
//
// Special Notes  :
//
// Creator        : Christina Warrender, SNL
//
// Creation Date  : 10/12/11
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

#ifndef Xyce_N_DEV_Synapse4_Factory_h
#define Xyce_N_DEV_Synapse4_Factory_h

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceTemplate.h>
#include <N_DEV_Synapse4.h>

//-----------------------------------------------------------------------------
// Class         : N_DEV_Synapse4
// Purpose       : Handles Synapse4s
//
// Special Notes :
//
// Creator       : Christina Warrender, SNL
// Creation Date : 10/12/11
//-----------------------------------------------------------------------------
struct N_DEV_Synapse4
{
  //---------------------------------------------------------------------------
  // Function      : factory
  // Purpose       : This is the factory function for the Synapse4 class.
  //
  // Special Notes : ERK.  10/16/2005.  This used to be a singleton (ie a
  //                 static pointer was returned) but had to be changed
  //                 so that the library version of Xyce would work
  //                 correctly.
  //
  // Scope         : public
  // Creator       : Christina Warrender, SNL
  // Creation Date : 10/12/11
  //---------------------------------------------------------------------------
  static N_DEV_Device * factory(N_DEV_SolverState & ss1,
                                N_DEV_DeviceOptions & do1)
  {
     const string name("Synapse");
     const string className("N_DEV_Synapse4");
     const string defaultModelName("YSYNAPSE level 4");

     vector< pair<string,double> > nameVec;

     N_DEV_Device * devptr =
       new N_DEV_SynapseMaster4
       //new Xyce::Device::DeviceTemplate<N_DEV_SynapseModel4, N_DEV_SynapseInstance4>
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


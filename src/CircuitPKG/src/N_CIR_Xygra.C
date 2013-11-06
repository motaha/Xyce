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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_CIR_Xygra.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Tom Russo, SNL, Electrical and Microsystems Modeling
//
// Creation Date  : 8/21/08
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.3.6.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:32 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

// ---------- Standard Includes ----------
// ----------   Xyce Includes   ----------
#include <N_CIR_Xygra.h>
#include <N_DEV_DeviceInterface.h>

//-----------------------------------------------------------------------------
// Function      : N_DEV_Xygra::N_DEV_Xygra
// Purpose       : ctor
// Special Notes : 
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/21/08
//-----------------------------------------------------------------------------

N_CIR_Xygra::N_CIR_Xygra()
  : N_CIR_Xyce()
{
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_Xygra::~N_DEV_Xygra
// Purpose       : dtor
// Special Notes : 
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/21/08
//-----------------------------------------------------------------------------

N_CIR_Xygra::~N_CIR_Xygra()
{
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_Xygra::getDeviceNames
// Purpose       : get all names of devices of specified type in the netlist
// Special Notes : "deviceType" takes a string of the form the device would
//                 have when instantiated in a netlist, e.g. "R1" for a resistor
//                 or "Y%XYGRA%DUMMY" for a Xygra device
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/25/08
//-----------------------------------------------------------------------------
bool N_CIR_Xygra::getDeviceNames(const string & deviceType, 
                      vector<string> & deviceNames)
{
  return devIntPtr_->getDeviceNames(deviceType,deviceNames);
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_Xygra::xygraGetNumNodes
// Purpose       : get number of nodes in specified Xygra device
// Special Notes : 
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/25/08
//-----------------------------------------------------------------------------
int N_CIR_Xygra::xygraGetNumNodes(const string & deviceName)
{
  return devIntPtr_->xygraGetNumNodes(deviceName);
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_Xygra::xygraGetNumWindings
// Purpose       : get number of Windings in specified Xygra device
// Special Notes : 
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/25/08
//-----------------------------------------------------------------------------
int N_CIR_Xygra::xygraGetNumWindings(const string & deviceName)
{
  return devIntPtr_->xygraGetNumWindings(deviceName);
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_Xygra::xygraGetCoilWindings
// Purpose       : get number of Windings in each coil of specified Xygra
//                 device
// Special Notes : 
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 9/15/08
//-----------------------------------------------------------------------------
void N_CIR_Xygra::xygraGetCoilWindings(const string & deviceName,
                                      vector<int> & cW)
{
  devIntPtr_->xygraGetCoilWindings(deviceName,cW);
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_Xygra::xygraGetCoilNames
// Purpose       : get names of each coil of specified Xygra
//                 device
// Special Notes : 
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 9/29/08
//-----------------------------------------------------------------------------
void N_CIR_Xygra::xygraGetCoilNames(const string & deviceName,
                                      vector<string> & cN)
{
  devIntPtr_->xygraGetCoilNames(deviceName,cN);
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_Xygra::xygraSetConductances
// Purpose       : Set the conductance matrix for a named Xygra device
// Special Notes : 
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/26/08
//-----------------------------------------------------------------------------
bool N_CIR_Xygra::xygraSetConductances(const string & deviceName, 
                                  const vector<vector<double> > & cM)
{
  return devIntPtr_->xygraSetConductances(deviceName, cM);
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_Xygra::xygraSetK
// Purpose       : Set the K matrix for a named Xygra device
// Special Notes : 
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/26/08
//-----------------------------------------------------------------------------
bool N_CIR_Xygra::xygraSetK(const string & deviceName, 
                            const vector<vector<double> > & kM,
                            const double t)
{
  return devIntPtr_->xygraSetK(deviceName, kM, t);
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_Xygra::xygraSetSources
// Purpose       : Set the S vector for a named Xygra device
// Special Notes : 
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 9/18/2008
//-----------------------------------------------------------------------------
bool N_CIR_Xygra::xygraSetSources(const string & deviceName, 
                                  const vector<double > & sV,
                                  const double t)
{
  return devIntPtr_->xygraSetSources(deviceName, sV, t);
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_Xygra::xygraGetVoltages
// Purpose       : retrieve nodal voltages for a named Xygra device
// Special Notes : Order of output voltages is 
//                 (coil1+,internal nodes for coil 1, coil1-,coil2+...)
// Scope         : public
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 9/18/2008
//-----------------------------------------------------------------------------
bool N_CIR_Xygra::xygraGetVoltages(const string & deviceName, 
                                   vector<double > & vN)
{
  return devIntPtr_->xygraGetVoltages(deviceName, vN);
}


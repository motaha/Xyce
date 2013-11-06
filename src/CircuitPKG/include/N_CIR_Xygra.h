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
// Filename       : $RCSfile: N_CIR_Xygra.h,v $
//
// Purpose        : Provide a class for Xyce/Alegra coupling
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
// Revision Date  : $Date: 2013/10/03 17:23:31 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef Xyce_N_CIR_XYGRA_H
#define Xyce_N_CIR_XYGRA_H

// ---------- Standard Includes ----------
// ----------   Xyce Includes   ----------
#include <N_CIR_Xyce.h>

//-----------------------------------------------------------------------------
// Class         : N_CIR_Xygra
// Purpose       :
// Special Notes :  
//                 
// Creator       : Tom Russo, SNL, Electrical and Microsystems Modeling
// Creation Date : 8/21/08
//-----------------------------------------------------------------------------
class N_CIR_Xygra : public N_CIR_Xyce
{

 public:
  N_CIR_Xygra();
  virtual ~N_CIR_Xygra();
  // and of course all of N_CIR_Xyce's base class public and protected 
  // methods...

  bool getDeviceNames(const string & deviceType, 
                      vector<string> & deviceNames);
  int xygraGetNumNodes(const string & deviceName);
  int xygraGetNumWindings(const string & deviceName);
  void xygraGetCoilWindings(const string & deviceName,
                            vector<int> & cW);
  void xygraGetCoilNames(const string & deviceName,
                            vector<string> & cN);
  bool xygraSetConductances(const string & deviceName, 
                            const vector< vector<double> > & cM);
  bool xygraSetK(const string & deviceName, 
                 const vector< vector<double> > & kM,
                 const double t=0);
  bool xygraSetSources(const string & deviceName, 
                       const vector< double > & sV,
                       const double t=0);
  bool xygraGetVoltages(const string & deviceName, 
                        vector< double > & vN);
} ;

#endif

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
// Filename       : $RCSfile: N_UTL_Version.h,v $
//
// Purpose        : set version string  
// Special Notes  : 
//
// Creator        : Eric Rankin
//
// Creation Date  : 
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.4.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:52 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------



#ifndef N_UTL_VERSION_H
#define N_UTL_VERSION_H



// ---------- Standard Includes -----------------------------------------------
#include <string>

using namespace std;

//-----------------------------------------------------------------------------
// Class         : N_UTL_Version
// Purpose       : Wrapper for retrieving Xyce version string info.
// Special Notes : 
// Creator       : 
// Creation Date : 
//-----------------------------------------------------------------------------
class N_UTL_Version 
{

public:

  // get full banner string for Xyce version
  static string getFullVersionString();
  // get the maj-min-rev number for Xyce version
  static string getShortVersionString();

  static string getBuildVariant();

  static string getCapabilities();

  static string getLicense();
  
private:
  
  // no need to create objects
  N_UTL_Version();
  
};
 
#endif 

//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2014 Sandia Corporation
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
// Filename       : $RCSfile: N_TOP_System.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 02/24/03
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.7 $
//
// Revision Date  : $Date: 2014/02/24 23:49:27 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef N_TOP_System_h
#define N_TOP_System_h 1

#include <iosfwd>
#include <map>
#include <string>

#include <N_UTL_Xyce.h>

class N_PDS_Manager;

namespace Xyce {
namespace Topo {

//-----------------------------------------------------------------------------
// Class         : Xyce::Topo::System
// Purpose       : Abstract System Interface
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 2/24/03
//-----------------------------------------------------------------------------
class System
{

 public:

  // Destructor
  virtual ~System() {}

  virtual System * clone() = 0;

  virtual Graph & graph() = 0;

  virtual Directory & directory( const string & type = "" ) = 0;

  virtual NodeFactory & nodeFactory() = 0;

  virtual MappingInterface & mapping() = 0;

 protected:

  //Default Constructor protected so only derived types can call it
  System() {}

 private:

  //Copy and Assignment private so only clone can be used
  System( const System & );
  System & operator=( const System & );

  friend ostream & operator<<( ostream & os, const System & );
};

} //namespace Topo
} //namespace Xyce

typedef Xyce::Topo::System N_TOP_System;

#endif

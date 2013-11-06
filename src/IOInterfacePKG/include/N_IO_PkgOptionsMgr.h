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
// Filename       : $RCSfile: N_IO_PkgOptionsMgr.h,v $
//
// Purpose        : Setup to register input options with pkg's that request them
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 01/28/03
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.8.6.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:42 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_IO_PkgOptionsMgr_h
#define Xyce_N_IO_PkgOptionsMgr_h

// ---------- Standard Includes ----------

#include <string>
#include <map>

#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RefCountPtr;
using Teuchos::rcp;

// ----------   Xyce Includes   ----------

#include <N_UTL_Xyce.h>
#include <N_UTL_Misc.h>

#include <N_UTL_OptionBlock.h>

#include <N_ERH_ErrorMgr.h>

// ---------- Forward Declarations ----------


//-----------------------------------------------------------------------------
// Class         : N_IO_PkgOptionsReg
// Purpose       : Abstract IF for pkg option registration
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 1/28/03
//-----------------------------------------------------------------------------
struct N_IO_PkgOptionsReg
{
  virtual ~N_IO_PkgOptionsReg() {}

  virtual bool operator()( const N_UTL_OptionBlock & options ) = 0;
};

//-----------------------------------------------------------------------------
// Class         : N_IO_PkgOptionsMgr
// Purpose       :
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 1/28/03
//-----------------------------------------------------------------------------
class N_IO_PkgOptionsMgr
{

 public:
  N_IO_PkgOptionsMgr();
  ~N_IO_PkgOptionsMgr();

  // registration functions:
  bool submitRegistration ( const string & regName, 
                            const string & circName, 
                            N_IO_PkgOptionsReg * reg );

  bool submitOptions( const N_UTL_OptionBlock & options, const string & circName );
  
 private:

  typedef pair<const string,N_IO_PkgOptionsReg*> RegistrationsValueType ;
  typedef pair<const string,N_UTL_OptionBlock> OptionsValueType;

  typedef multimap<string,N_IO_PkgOptionsReg*> RegistrationsMap;
  typedef multimap<string,N_UTL_OptionBlock> OptionsMap;

  typedef RegistrationsMap::iterator RegistrationsIterator;
  typedef OptionsMap::iterator OptionsIterator;

  typedef map <string, RegistrationsMap > circRegsMap;
  typedef map <string, OptionsMap> circOptsMap;
  typedef circRegsMap::iterator circRegsIterator;
  typedef circOptsMap::iterator circOptsIterator;

  circRegsMap circSpecificRegs_;
  circOptsMap circSpecificOpts_;

};

#endif


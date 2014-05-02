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
// Revision Number: $Revision: 1.14.2.1 $
//
// Revision Date  : $Date: 2014/02/26 20:42:38 $
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

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Class         : PkgOptionsReg
// Purpose       : Abstract IF for pkg option registration
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 1/28/03
//-----------------------------------------------------------------------------
struct PkgOptionsReg
{
  virtual ~PkgOptionsReg() {}

  virtual bool operator()( const Util::OptionBlock & options ) = 0;
};

//-----------------------------------------------------------------------------
// Class         : PkgOptionsMgr
// Purpose       :
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 1/28/03
//-----------------------------------------------------------------------------
class PkgOptionsMgr
{

public:
  PkgOptionsMgr();
  ~PkgOptionsMgr();

  // registration functions:
  bool submitRegistration ( const std::string & regName, 
                            const std::string & circName, 
                            PkgOptionsReg * reg );

  bool submitOptions( const Util::OptionBlock & options, const std::string & circName );
  
private:

  typedef std::pair<const std::string,PkgOptionsReg*> RegistrationsValueType ;
  typedef std::pair<const std::string,Util::OptionBlock> OptionsValueType;

  typedef std::multimap<std::string,PkgOptionsReg*> RegistrationsMap;
  typedef std::multimap<std::string,Util::OptionBlock> OptionsMap;

  typedef RegistrationsMap::iterator RegistrationsIterator;
  typedef OptionsMap::iterator OptionsIterator;

  typedef std::map<std::string, RegistrationsMap > circRegsMap;
  typedef std::map<std::string, OptionsMap> circOptsMap;
  typedef circRegsMap::iterator circRegsIterator;
  typedef circOptsMap::iterator circOptsIterator;

  circRegsMap circSpecificRegs_;
  circOptsMap circSpecificOpts_;

};

} // namespace IO
} // namespace Xyce

typedef Xyce::IO::PkgOptionsMgr N_IO_PkgOptionsMgr;
typedef Xyce::IO::PkgOptionsReg N_IO_PkgOptionsReg;

#endif


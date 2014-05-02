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
// Filename       : $RCSfile: N_IO_PkgOptionsMgr.C,v $
//
// Purpose        : Used to send input options to registered pkg's
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL
//
// Creation Date  : 1/28/2003
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.21 $
//
// Revision Date  : $Date: 2014/02/24 23:49:22 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <N_IO_PkgOptionsMgr.h>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Function      : PkgOptionsMgr::PkgOptionsMgr
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, SNL
// Creation Date : 10/21/2008
//-----------------------------------------------------------------------------
PkgOptionsMgr::PkgOptionsMgr()
{
}

//-----------------------------------------------------------------------------
// Function      : PkgOptionsMgr::~PkgOptionsMgr
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL
// Creation Date : 1/28/2003
//-----------------------------------------------------------------------------
PkgOptionsMgr::~PkgOptionsMgr()
{
  circRegsIterator iterCR = circSpecificRegs_.begin();
  circRegsIterator  endCR = circSpecificRegs_.end  ();

  for ( ; iterCR != endCR; ++iterCR)
  {
    RegistrationsIterator iterR = iterCR->second.begin();
    RegistrationsIterator endR = iterCR->second.end();
    for( ; iterR != endR; ++iterR )
    {
      delete iterR->second;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : PkgOptionsMgr::submitRegistration
// Purpose       : Register an object for recieving options
//
// Special Notes : The design of this function and the other 'submit' function
//                 is intended to be robust enough that it doesn't matter
//                 which call is done first.  If the options are submitted
//                 first, then the submitRegistration function will send
//                 them out to the various registration functors.  If the
//                 registration functors are submitted first, then when the
//                 options are submitted, they'll be sent out.
//
// Scope         : public
// Creator       : Robert Hoekstra, SNL
// Creation Date : 1/28/2003
//-----------------------------------------------------------------------------
bool PkgOptionsMgr::submitRegistration( const std::string & regName,
                                             const std::string & circName,
                                             N_IO_PkgOptionsReg * reg )
{

  // if a registration for this particular circuit hasn't happened yet, then
  // create one.
  if (circSpecificRegs_.find(circName) ==  circSpecificRegs_.end() )
  {
    RegistrationsMap registrations;
    circSpecificRegs_[circName] = registrations;
  }

  circSpecificRegs_[circName].insert( RegistrationsValueType(regName,reg) );

  if (circSpecificOpts_.find(circName) != circSpecificOpts_.end() )
  {
    if( circSpecificOpts_[circName].count( regName ) )
    {
      OptionsIterator iterO = circSpecificOpts_[circName].lower_bound( regName );
      OptionsIterator  endO = circSpecificOpts_[circName].upper_bound( regName );
      for( ; iterO != endO; ++iterO )
      {
        (*reg)( iterO->second );
      }
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : PkgOptionsMgr::submitOptions
// Purpose       : Register an option block to be given to pkg's
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL
// Creation Date : 1/28/2003
//-----------------------------------------------------------------------------
bool PkgOptionsMgr::submitOptions( const Util::OptionBlock & options, const std::string & circName )
{

  // if a registration for this particular circuit hasn't happened yet, then
  // create one.
  if (circSpecificOpts_.find(circName) ==  circSpecificOpts_.end() )
  {
    OptionsMap tmpOpt;
    circSpecificOpts_[circName] = tmpOpt;
  }

  Util::OptionBlock newOptions( options );
  std::string oname(newOptions.getName());

  Util::OptionBlock::ParameterList::iterator iterO = newOptions.begin();
  Util::OptionBlock::ParameterList::iterator  endO = newOptions.end();
  for( ; iterO != endO; ++iterO )
  {
    if( iterO->tag() == "SINGLETON" && iterO->getImmutableValue<int>() )
    {
      if( circSpecificOpts_[circName].count( oname ) )
      {
        std::string msg("Only 1 option line of type: " + oname + " allowed!\n");
        N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL_0, msg );
      }
      else
      {
        //strip out parameter befor passing on
        iterO = newOptions.getParams().erase( iterO );
      }
    }
  }

  circSpecificOpts_[circName].insert( OptionsValueType( oname, newOptions ) );

  if (circSpecificRegs_.find(circName) !=  circSpecificRegs_.end() )
  {
    if( circSpecificRegs_[circName].count( oname ) )
    {
      RegistrationsIterator iterR = circSpecificRegs_[circName].lower_bound( oname );
      RegistrationsIterator endR = circSpecificRegs_[circName].upper_bound( oname );
      for( ; iterR != endR; ++iterR )
      {
        (*iterR->second)( newOptions );
      }
    }
  }

  return true;
}

} // namespace IO
} // namespace Xyce

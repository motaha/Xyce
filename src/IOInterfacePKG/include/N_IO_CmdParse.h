//-------------------------------------------------------------------------
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
//-------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Filename      : N_IO_CmdParse.h
//
// Purpose       : This file defines the command line parser class.
//
// Special Notes :
//
//
// Creator       : Eric Keiter
//
// Creation Date : 06/17/99
//
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.36.2.2 $
//
// Revision Date  : $Date: 2014/03/10 16:15:21 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#ifndef _N_IO_CMDPARSE_
#define _N_IO_CMDPARSE_

#include <iosfwd>
#include <string>
#include <map>

#include <N_IO_fwd.h>
#include <N_DEV_fwd.h>
#include <N_UTL_Xyce.h>

class N_PDS_Manager;

namespace Xyce {
namespace IO {

void usage(std::ostream &os);

// KRS, 12/10/07: the following functions are related to some private variables that actually have nothing to do with
// parsing the command line.  In order to make Xyce produce netlist files that contain resistors between ground and
// either a node that has no DC path to ground or a node that is only connected to one device terminal, some information
// needs to be passed between the I/O package and the topology package.  Since both packages have access to the command
// line arguments, I'm kludging the information I need to pass here.

class HangingResistor {
public:
  HangingResistor()
    : netlistCopy_(false),
      oneTerm_(false),
      noDCPath_(false),
      oneTermRes_(""),
      noDCPathRes_("")
  {}

  HangingResistor(const HangingResistor &right)
    : netlistCopy_(right.netlistCopy_),
      oneTerm_(right.oneTerm_),
      noDCPath_(right.noDCPath_),
      oneTermRes_(right.oneTermRes_),
      noDCPathRes_(right.noDCPathRes_)
  {}

  bool getNetlistCopy() 
  {
    return netlistCopy_;
  }

  bool getOneTerm() 
  {
    return oneTerm_;
  }

  bool getNoDCPath() 
  {
    return noDCPath_;
  }

  std::string getOneTermRes() 
  {
    return oneTermRes_;
  }

  std::string getNoDCPathRes() 
  {
    return noDCPathRes_;
  }

  void setNetlistCopy(bool netlistCopyArg) 
  {
    netlistCopy_ = netlistCopyArg;
  }

  void setOneTerm(bool oneTermArg) 
  {
    oneTerm_ = oneTermArg;
  }

  void setNoDCPath(bool noDCPathArg) 
  {
    noDCPath_ = noDCPathArg;
  }

  void setOneTermRes(const std::string & oneTermResArg) 
  {
    oneTermRes_ = oneTermResArg;
  }

  void setNoDCPathRes(const std::string & noDCPathResArg) 
  {
    noDCPathRes_ = noDCPathResArg;
  }

private:
  //KRS, 12/10/07:  private variables introduced to make copies of netlists
  //which automatically add resistors between "dangling" nodes and ground
  bool        netlistCopy_;
  bool        oneTerm_;
  bool        noDCPath_;
  std::string oneTermRes_;
  std::string noDCPathRes_;
};


//-----------------------------------------------------------------------------
// Class         : CmdParse
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/14/01
//-----------------------------------------------------------------------------
class CmdParse
{
public:
  CmdParse();

  virtual ~CmdParse();

  CmdParse(const CmdParse &);

private:
  CmdParse &operator=(const CmdParse &);

public:
  bool registerParallelMgr( N_PDS_Manager * parMgr ) 
  {
    parMgr_ = parMgr;

    return true;
  }

  int argc() const 
  {
    return iargs;
  }

  char **argv() const 
  {
    return cargs;
  }

  void setCommandArgs ();

  int parseCommandLine(int iargs, char **cargs);

  bool argExists(const std::string & arg_tmp) const;

  std::string getArgumentValue(const std::string & argumentName) const;

  // Returns the number of cmd line args.
  int numArgs();

  void printArgMap();

  void getArg(int i, std::string & arg);

  void setNetlist(const std::string & newNetlist);

  HangingResistor &getHangingResistor()
  {
    return hangingResistor_;
  }

private:
  bool isSwitchArg(const std::string &arg);
  bool isStringValuedArg(const std::string &arg);

  void handleParameterOutputs_(const std::string & tmp, int &i, int &ia, char** ca);

private :
  int         iargs;
  char **     cargs;
  bool        allocatedCargs_;

  N_PDS_Manager *     parMgr_;

  std::map<std::string, int>          swArgs;
  std::map<std::string, std::string>  stArgs;

  std::map<std::string, int>          argIndex;

  HangingResistor             hangingResistor_;
};

} // namespace IO
} // namespace Xyce

typedef Xyce::IO::CmdParse N_IO_CmdParse;

#endif

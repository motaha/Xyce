//-------------------------------------------------------------------------
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
// Revision Number: $Revision: 1.27.2.3 $
//
// Revision Date  : $Date: 2013/10/03 17:23:42 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef _N_IO_CMDPARSE_
#define _N_IO_CMDPARSE_

#include <string>
#include <map>
#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RefCountPtr;
using Teuchos::rcp;
using Teuchos::rcpFromRef;

#include <N_DEV_fwd.h>
#include <N_UTL_Xyce.h>
#include <N_PDS_Manager.h>

//-----------------------------------------------------------------------------
// Class         : N_IO_CmdParse
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/14/01
//-----------------------------------------------------------------------------
class N_IO_CmdParse
{

public:

  // constructor
  N_IO_CmdParse();

  // special copy constructor
  N_IO_CmdParse(N_IO_CmdParse & right);

  // destructor
  virtual ~N_IO_CmdParse();

  void setCommandArgs ();
  bool registerParallelMgr( N_PDS_Manager * parMgr )
  {
    // Using rcpFromRef to create a NON-OWNING RCP using an undefined type.
    // Note:  N_PDS_Manager is forward declared above.
    RefCountPtr<N_PDS_Manager> parMgrRCPtr = rcpFromRef(*parMgr);
    registerParallelMgr(parMgrRCPtr);
    return true;
  }

  bool registerParallelMgr( RefCountPtr<N_PDS_Manager> parMgr )
  {
    parMgr_ = parMgr;
    return true;
  }

  void parseCommandLine(int iargs, char **cargs);

//  bool argExists(const char *arg_tmp) const;
  bool argExists(const string & arg_tmp) const;

  string getArgumentValue(const string & argumentName) const;

  // Returns the number of cmd line args.
  int numArgs();

  void printArgMap();

  void getArg(int i, string & arg);

  void setNetlist(string & newNetlist);

  //KRS, 12/10/07:  the following functions are related to some private
  //variables that actually have nothing to do with parsing the command line.
  //In order to make Xyce produce netlist files that contain resistors between
  //ground and either a node that has no DC path to ground or a node that is
  //only connected to one device terminal, some information needs to be passed
  //between the I/O package and the topology package.  Since both packages have
  //access to the command line arguments, I'm kludging the information I need
  //to pass here.

  bool getNetlistCopy();
  bool getOneTerm();
  bool getNoDCPath();
  string getOneTermRes();
  string getNoDCPathRes();

  void setNetlistCopy(bool netlistCopyArg);
  void setOneTerm(bool oneTermArg);
  void setNoDCPath(bool noDCPathArg);
  void setOneTermRes(const string & oneTermResArg);
  void setNoDCPathRes(const string & noDCPathResArg);

protected:

private :

  int iargs;
  char ** cargs;

  std::string usage_;

  RefCountPtr<N_PDS_Manager> parMgr_;

  map<string,int> swArgs;
  map<string,string> stArgs;

  map<string,int> argIndex;

  bool isSwitchArg(string arg);
  bool isStringValuedArg(string arg);

  void handleParameterOutputs_ (string & tmp, int &i, int &ia, char** ca);

  //KRS, 12/10/07:  private variables introduced to make copies of netlists
  //which automatically add resistors between "dangling" nodes and ground
  bool netlistCopy_;
  bool oneTerm_;
  bool noDCPath_;
  bool allocatedCargs_;
  string oneTermRes_;
  string noDCPathRes_;

  friend class Xyce::Device::XyceInterface;

};

#endif

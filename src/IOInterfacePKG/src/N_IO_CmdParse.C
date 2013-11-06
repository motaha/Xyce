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
// Filename      : N_IO_CmdParse.C
//
// Purpose       : This file contains the command line parser class.
//
// Special Notes :
//
// Creator       : Eric Keiter
//
// Creation Date : 06/17/99
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.101.2.8 $
//
// Revision Date  : $Date: 2013/10/03 17:23:43 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

#include <iostream>
#include <sstream>

#ifdef HAVE_STRINGS_H
#include <strings.h>
#else
#ifdef HAVE_STRING_H
// On windows, we don't have this and cannot use strcasecmp, and need to use
// _stricmp instead.  Groan.
#include <string.h>
#else
error FIXME: Neither strings.h (POSIX) nor string.h (Windows) found.  Cannot find a case-insensitive string compare to use.
#endif  // HAVE_STRING_H
#endif  // HAVE_STRINGS_H

// ----------   Xyce Includes   ----------

#include <N_ERH_ErrorMgr.h>

#include <N_IO_CmdParse.h>

#include <N_UTL_Misc.h>
#include <N_UTL_Expression.h>

#include <N_PDS_Comm.h>

#include <N_DEV_fwd.h>
#include <N_DEV_DeviceInterface.h>
#include <N_DEV_DeviceMgr.h>
#include <N_DEV_Device.h>
#include <N_DEV_DeviceModel.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_OutputPars.h>
#include <N_DEV_MatrixLoadData.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_Factory.h>
#include <N_UTL_LogStream.h>
#include <N_UTL_IndentStreamBuf.h>

#include <N_UTL_Version.h>
#include <N_UTL_IndentStreamBuf.h>

struct DeviceEntityCmp : public std::binary_function<N_DEV_DeviceEntity, N_DEV_DeviceEntity, bool>
{
  static bool less_nocase(const std::string &lhs, const std::string &rhs) {
#ifdef HAVE_STRINGS_H
    return strcasecmp(lhs.c_str(), rhs.c_str()) < 0;
#else
    return _stricmp(lhs.c_str(), rhs.c_str()) < 0;
#endif
  }

  bool operator()(const N_DEV_DeviceEntity &entity_0, const N_DEV_DeviceEntity &entity_1) const {
    return less_nocase(entity_0.getName(), entity_1.getName());
  }
  bool operator()(const N_DEV_DeviceEntity *entity_0, const N_DEV_DeviceEntity *entity_1) const {
    return less_nocase(entity_0->getName(), entity_1->getName());
  }
};

struct DeviceCmp : public std::binary_function<N_DEV_Device, N_DEV_Device, bool>
{
  static bool less_nocase(const std::string &lhs, const std::string &rhs) {
#ifdef HAVE_STRINGS_H
    return strcasecmp(lhs.c_str(), rhs.c_str()) < 0;
#else
    return _stricmp(lhs.c_str(), rhs.c_str()) < 0;
#endif
  }
};

namespace {
void tolower(std::string &s)
{
  for (string::iterator it = s.begin(); it != s.end(); ++it)
    (*it) = ((*it) < 'A'  || (*it) > 'Z') ? (*it) : (*it) - 'A' + 'a';
}

}

//-----------------------------------------------------------------------------
// Function      : N_IO_CmdParse::N_IO_CmdParse
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Eric Keiter
// Creation Date : 02/14/01
//-----------------------------------------------------------------------------
N_IO_CmdParse::N_IO_CmdParse()
: iargs(0),
  cargs(0),
  netlistCopy_(false),
  oneTerm_(false),
  noDCPath_(false),
  oneTermRes_(""),
  noDCPathRes_(""),
  allocatedCargs_(false)
{
  setCommandArgs();

  // Set up argument list for usage statement;
  std::ostringstream oss;

  oss << "Usage: Xyce [arguments] netlist\n\n"
      << "Arguments:\n"
      << "  -h                          print usage and exit\n"
      << "  -v                          print version info and exit\n"
      << "  -capabilities               print compiled-in options and exit\n"
      << "  -license                    print license and exit\n"
      << "  -param [Device [level [<inst|mod>]]] print a terse summary of model and/or device parameters\n"
      << "  -info [Device [level [<inst|mod>]]]  print a description of model and/or device parameters\n"
      << "  -syntax                     check netlist syntax and exit\n"
      << "  -norun                      netlist syntax and topology and exit\n"
      << "  -gui                        gui file output\n"
      << "  -jacobian_test              jacobian matrix diagnostic\n"
      //      << "  -test                       run unit tests\n"
      << "  -delim <TAB|COMMA|string>   set the output file field delimiter\n"
      << "  -o <file>                   place the results into <file>\n"
      << "  -l <file>                   place the log output into <file>\n"
      << "  -remeasure [existing Xyce output file] recompute .measure() results with existing data\n"

#ifndef Xyce_PARALLEL_MPI
      << "  -r <file>                   generate a rawfile named <file> in binary format\n"
      << "  -a                          use with  -r <file>  to output in ascii format\n"

#endif
      << "  -nox <on|off>               NOX nonlinear solver usage\n"

      << "  -linsolv <solver>           force usage of specific linear solver\n"

      << "  -maxord <1..5>              maximum time integration  order\n"
    //<< "  -method <1..4>              time integration method (old-dae only)\n"

#ifdef Xyce_Dakota
      << "  -dakota <dakota input file> dakota input file for this simulation\n"
#endif
      << "  -prf <param file name>      specify a file with simulation parameters\n"
      << "  -rsf <response file name>   specify a file to save simulation responses functions.\n"

#ifdef Xyce_DEBUG_NONLINEAR
      << "  -ndl <0..N>                 set the nonlinear solver debug level (overrides netlist)\n"
#endif
#ifdef Xyce_DEBUG_TIME
      << "  -tdl <0..N>                 set the time integration debug level\n"
#endif
#ifdef Xyce_DEBUG_DEVICE
      << "  -ddl <0..N>                 set the device debug level (overrides netlist)\n"
      << "  -sdl <0..N>                 set the device sensitivity debug level\n"
#endif

#ifdef Xyce_TOTALVIEW_BOGON
      << "  -mpichtv                    required for totalview \n"
#endif
      << std::endl;

  usage_ = oss.str();
}

//-----------------------------------------------------------------------------
// Function      : N_IO_CmdParse::N_IO_CmdParse
// Purpose       : copy constructor
// Special Notes :
// Scope         : Public
// Creator       : Eric Keiter
// Creation Date : 01/07/07
//-----------------------------------------------------------------------------
N_IO_CmdParse::N_IO_CmdParse(N_IO_CmdParse & right) :
  iargs(right.iargs),
  swArgs(right.swArgs),
  stArgs(right.stArgs),
  argIndex(right.argIndex),
  usage_(right.usage_),
  parMgr_(right.parMgr_),
  netlistCopy_(right.netlistCopy_),
  oneTerm_(right.oneTerm_),
  noDCPath_(right.noDCPath_),
  oneTermRes_(right.oneTermRes_),
  noDCPathRes_(right.noDCPathRes_)
{
  // if we use the copy constructor, then we will own
  // cargs and must delete it. There should be a better
  // way to do this with ref counted pointers, but I'm still
  // working on it. RLS
  allocatedCargs_ = true;
  cargs = new char*[iargs];
  int i=0;
  int j=0;
  for (i=0;i<iargs;++i)
  {
    if (right.cargs[i] == NULL)
    {
      cargs[i] = NULL;
      continue;
    }

    string tmpString(right.cargs[i]);
    int size = tmpString.size()+2;
    cargs[i] = new char[size];
    for (j=0;j<size;++j) cargs[i][j] = 0;

    sprintf(cargs[i], "%s", tmpString.c_str());
  }
}

//-----------------------------------------------------------------------------
// Function      : N_IO_CmdParse::N_IO_CmdParse
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Eric Keiter
// Creation Date : 02/14/01
//-----------------------------------------------------------------------------
N_IO_CmdParse::~N_IO_CmdParse()
{
  if(allocatedCargs_)
  {
    // by defalt this class will just hold on to the char** pointer
    // from the os for cargs.  However, if we use the copy constructor
    // then we actually allocate and populate cargs.  If this is the case
    // then we need to delete it here
    for( int i=0; i<iargs; i++)
    {
      if( cargs[i] )
        delete [] cargs[i];
    }
    if( cargs )
      delete [] cargs;
  }
}


//-----------------------------------------------------------------------------
// Function      : N_IO_CmdParse::setCommandArgs
// Purpose       : Initialize command line argumentws
// Special Notes :
// Scope         : Public
// Creator       : Dave Shirley, PSSi
// Creation Date : 01/24/11
//-----------------------------------------------------------------------------
void N_IO_CmdParse::setCommandArgs ()
{
  stArgs.clear();
  swArgs.clear();
  argIndex.clear();

  // Set the maps containing the expected arguments. There are two types
  // of arguments switch arguments and string valued arguments.
  swArgs[ "-h" ] = 0;           // Help option, print usage and exit.
  swArgs[ "-test" ] = 0;
  swArgs[ "-v" ] = 0;
  swArgs[ "-capabilities" ] = 0;
  swArgs[ "-license" ] = 0;
  swArgs[ "-syntax" ] = 0;
  swArgs[ "-norun" ] = 0;
  swArgs[ "-gui" ] = 0;
  swArgs[ "-jacobian_test" ] = 0;

  stArgs[ "netlist" ] = "";     // The input netlist will be a special case.
  stArgs[ "-delim" ] = "";
  stArgs[ "-o" ] = "";
  stArgs[ "-l" ] = "";          // Output log information to a file.

  stArgs[ "-r" ] = "";          // Output binary rawfile.
  swArgs[ "-a" ] = 0;           // Use ascii instead of binary in rawfile output

  stArgs[ "-nox" ] = "";

  stArgs[ "-linsolv" ] = "";

  stArgs[ "-maxord" ] = "";
  //stArgs[ "-method" ] = "";

#ifdef Xyce_Dakota
  stArgs[ "-dakota" ] = "";     // specify a dakota input file to couple with this simulation
#endif
  stArgs[ "-prf" ] = "";        // specify a parameter input file to set runtime params from a file
  stArgs[ "-rsf" ] = "";        // specify a response output file to save results to a file

#ifdef Xyce_DEBUG_NONLINEAR
  stArgs[ "-ndl" ] = "";
#endif
#ifdef Xyce_DEBUG_TIME
  stArgs[ "-tdl" ] = "";
#endif
#ifdef Xyce_DEBUG_DEVICE
  stArgs[ "-ddl" ] = "";
  stArgs[ "-sdl" ] = "";
#endif
#ifdef Xyce_TOTALVIEW_BOGON
  swArgs[ "-mpichtv" ] = 0;
#endif
  stArgs[ "-remeasure" ] = "";  // specify a existing data file on which Xyce will recompute .measure functions.
  return;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_CmdParse::setNetlist
// Purpose       : This function is used for 2-level newton solves.
//
// Special Notes : For 2-level newton, the inner solve is handled by a new
//                 allocation of N_CIR_Xyce, which needs to be passed in
//                 a set of command line args.  Those command line args
//                 need to be identical to the ones passed in by the user,
//                 except that the netlist needs to be different, as the
//                 inner problem is described in a different file.
//
//                 This function is called after a copy of the original
//                 N_IO_CmdParse class has been made.
//
//
// Scope         : Public
// Creator       : Eric Keiter
// Creation Date : 01/07/07
//-----------------------------------------------------------------------------
void N_IO_CmdParse::setNetlist(string & newNetlist)
{
  int netIndex = 0;

  if ( argIndex.find("netlist") != argIndex.end())
  {
    netIndex = argIndex["netlist"];
  }
  else
  {
    string errorMsg("N_IO_CmdParse::setNetlist. Unable to find netlist argument.");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, errorMsg);
  }


#if 0
  //++netIndex;

  if ( netIndex >= iargs )
  {
    // Unexectedly ran out of arguments on the command line.
    string errorMsg("N_IO_CmdParse::setNetlist.  Did not find previous netlist setting.");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, errorMsg);
  }
  else if (cargs[netIndex][0] == '-')
  {
    // Error if we ran into another option here.
    string errorMsg("N_IO_CmdParse::setNetlist.  Expected option value,");
    errorMsg += " but found option: " + string(cargs[netIndex]) + "\n";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, errorMsg);
  }
  else // found it!
#endif

  {
    delete [] cargs[netIndex];

    int newSize = newNetlist.size()+2;
    cargs[netIndex] = new char[newSize];
    for (int i=0;i<newSize;++i)  cargs[netIndex][i] = 0;

    sprintf(cargs[netIndex], "%s", newNetlist.c_str());

    stArgs["netlist"] = newNetlist;
  }

}

//-----------------------------------------------------------------------------
// Function      : N_IO_CmdParse::parseCommandLine
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Lon Waters
// Creation Date : 03/13/02
//-----------------------------------------------------------------------------
void N_IO_CmdParse::parseCommandLine(int ia, char** ca)
{
  iargs = ia;
  cargs = ca;

  N_PDS_Comm & comm = (*parMgr_->getPDSComm());

  bool isSerial = comm.isSerial();
  int procID = comm.procID();

  bool allExit = false;

#ifdef Xyce_PARALLEL_MPI
  N_ERH_ErrorMgr::startSafeBarrier();   //  All procs call (1)
#endif
  setCommandArgs();

  if( !procID )
  {
    for (int i = 1; i < ia; ++i)
    {
      // MPICH sometimes creates arguments that get replaced by NULLs by
      // MPI_INIT :-P.  They are always after the arguments that we care about,
      // so stop command line parsing at that point.
      if (ca[i] == NULL) break;

      string tmpArg(ca[i]);

      argIndex[tmpArg] = i;

      if (tmpArg[0] == '-')
      {
        if( tmpArg == "-h" )
        {
          cout << usage_ << endl;
          allExit = true;
          break;
        }
        else if( tmpArg == "-v" )
        {
          cout << N_UTL_Version::getFullVersionString() << endl;
          allExit = true;
        }
        else if( tmpArg == "-capabilities" )
        {
          cout << N_UTL_Version::getCapabilities() << endl;
          allExit = true;
        }
        else if( tmpArg == "-license" )
        {
          cout << N_UTL_Version::getLicense() << endl;
          allExit = true;
        }
        else if( tmpArg == "-param" || tmpArg == "-info" ||
                 tmpArg == "-doc"   || tmpArg == "-doc_cat" )
        {
          handleParameterOutputs_(tmpArg, i, ia, ca);
          allExit=true;
        }
        else if( tmpArg == "-gui")
        {
           string msg("Xyce setting up problem");
           N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::GUI_PROGRESS, msg);
           swArgs[tmpArg] = i;
        }
        else if (argExists(tmpArg))
        {
          string errorMsg("More than one \"" + tmpArg);
          errorMsg += "\" command line argument found.";
          errorMsg += "  Using last one.\n";
          cerr << endl << errorMsg;
          // Note: following command results in seg fault, don't know why.
          //N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_WARNING, errorMsg );
        }
        else if (isSwitchArg(tmpArg))
        {
          swArgs[tmpArg] = i;
        }
        else if (isStringValuedArg(tmpArg))
        {
          ++i;

          if ( i >= ia )
          {
            // Unexectedly ran out of arguments on the command line.
            cerr << usage_ << endl;
            string errorMsg("Did not find required value for option ");
            errorMsg += string(ca[i-1]) + "\n";
            N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL_0, errorMsg);
          }
          else if (ca[i][0] == '-')
          {
            // Error if we ran into another option here.
            string errorMsg("Expected option value,");
            errorMsg += " but found option: " + string(ca[i]) + "\n";
            N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL_0, errorMsg);
          }
          else
          {
            stArgs[tmpArg] = ca[i];
          }
        }
        else
        {
          // Invalid option, stop here.
          cerr << usage_ << endl;
          string errorMsg("Invalid option given: " + tmpArg + "\n");
          N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL_0, errorMsg);
        }
      }
      else
      {
        if (stArgs["netlist"] == "")
        {
          // Assume this field is the input netlist.
          stArgs["netlist"] = tmpArg;
          argIndex["netlist"] = i;
        }
        else
        {
          // Already found netlist, report error and terminate.
          // or this could be the case of a Dakota call to Xyce of
          // the format Xyce mycircuit.cir -prf params.in results.out
          // The more correct way would be to use the "-rsf <filename>"
          // to specify the results file.
          // This is a bit of a hack, but if the extra netlist starts with
          // res* (as in results* or response*) then assume it is the response
          // file and save it with the tag "-rsf" (provided that tag has not
          // been used).
          if ((stArgs["-rsf"]=="") &&
              ((tmpArg.find_first_of("res") != std::string::npos)
                || (tmpArg.find_first_of("RES") != std::string::npos)) )
          {
            // this is likely a response file name.  Store it as such
            stArgs["-rsf"]=tmpArg;
          }
          else
          {
            cerr << usage_ << endl;
            string errorMsg("Found second netlist on command line: " + tmpArg + "\n");
  #ifndef Xyce_TOTALVIEW_BOGON
            N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL_0, errorMsg);
  #else
            N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_WARNING_0, errorMsg);
  #endif
          }
        }
      }
    }

    if( !allExit && ( stArgs["netlist"] == "" ) )
    {
      // No netlist found on command line.
      cerr << usage_ << endl;
      string errorMsg("Netlist not found on command line\n");
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL_0, errorMsg);
    }
  }

#ifdef Xyce_PARALLEL_MPI
  N_ERH_ErrorMgr::safeBarrier(0);   //  All procs call (1)
#endif

  // Push cmds to other procs for parallel
  if( !isSerial )
  {
    int size, length, val, pos;
    int bsize = 500;

    char * buf = new char[bsize];

    if( !procID )
    {
      // broadcast termination flag to other procs and shutdown locally
      if( allExit )
      {
        size = -1;
        comm.bcast( &size, 1, 0 );
        Xyce_exit( 0 );
      }

      size = swArgs.size();
      comm.bcast( &size, 1, 0 );

      map<string,int>::iterator it_siM = swArgs.begin();
      map<string,int>::iterator end_siM = swArgs.end();
      for( ; it_siM != end_siM; ++it_siM )
      {
        length = it_siM->first.size();
        comm.bcast( &length, 1, 0 );
        pos = 0;
        comm.pack( it_siM->first.c_str(), length, buf, bsize, pos );
        comm.bcast( buf, length, 0 );

        val = it_siM->second;
        comm.bcast( &val, 1, 0 );
      }

      int size = stArgs.size();
      comm.bcast( &size, 1, 0 );

      map<string,string>::iterator it_ssM = stArgs.begin();
      map<string,string>::iterator end_ssM = stArgs.end();
      for( ; it_ssM != end_ssM; ++it_ssM )
      {
        length = it_ssM->first.size();
        comm.bcast( &length, 1, 0 );
        pos = 0;
        comm.pack( it_ssM->first.c_str(), length, buf, bsize, pos );
        comm.bcast( buf, length, 0 );

        length = it_ssM->second.size();
        comm.bcast( &length, 1, 0 );
        pos = 0;
        comm.pack( it_ssM->second.c_str(), length, buf, bsize, pos );
        comm.bcast( buf, length, 0 );
      }
    }

    // receive cmds from proc 0
    else
    {
      string str1, str2;

      comm.bcast( &size, 1, 0 );

      // shutdown nodes individually
      if( -1 == size )
      {
        Xyce_exit( 0 );
      }

      for( int i = 0; i < size; ++i )
      {
        comm.bcast( &length, 1, 0 );
        comm.bcast( buf, length, 0 );
        str1 = string( buf, length );

        comm.bcast( &val, 1, 0 );

        swArgs[str1] = val;
      }

      comm.bcast( &size, 1, 0 );
      for( int i = 0; i < size; ++i )
      {
        comm.bcast( &length, 1, 0 );
        comm.bcast( buf, length, 0 );
        str1 = string( buf, length );

        comm.bcast( &length, 1, 0 );
        comm.bcast( buf, length, 0 );
        str2 = string( buf, length );

        stArgs[str1] = str2;
      }
    }

    delete [] buf;

#ifdef Xyce_DEBUG_PARALLEL_CMDLINE
    cout << "*************************\n";
    cout << "PROC ID: " << procID << endl;
    cout << "Switched Args\n";
    cout << "-------------\n";
    map<string,int>::iterator it_siM = swArgs.begin();
    map<string,int>::iterator end_siM = swArgs.end();
    for( ; it_siM != end_siM; ++it_siM )
      cout << " " << it_siM->first << ": " << it_siM->second << endl;
    cout << "-------------\n";
    cout << "String Args\n";
    cout << "-------------\n";
    map<string,string>::iterator it_ssM = stArgs.begin();
    map<string,string>::iterator end_ssM = stArgs.end();
    for( ; it_ssM != end_ssM; ++it_ssM )
      cout << " " << it_ssM->first << ": " << it_ssM->second << endl;
    cout << "-------------\n";
    cout << "*************************\n";
#endif
  }

  // check for shutdown
  if( allExit && isSerial ) Xyce_exit( 0 );
}

//-----------------------------------------------------------------------------
// Function      : N_IO_CmdParse::handleParameterOutputs_
//
// Purpose       : Handles command line arguments like -doc, to output model
//                 parameters either to stdout or to files, depending on which
//                 option is chosen.
//
// Special Notes : erkeite: This function is called from parseCommandLine.  I
//                 moved this code out of that function because that function had
//                 gotten too long and confusing.
//
// Scope         : Public
// Creator       : Eric Keiter
// Creation Date : 10/31/07
//-----------------------------------------------------------------------------
void N_IO_CmdParse::handleParameterOutputs_ (string & tmpArg, int &i, int &ia, char** ca)
{
  Xyce::Device::OutputMode::Mode format = Xyce::Device::OutputMode::DEFAULT;
  if (tmpArg == "-param")
    format = Xyce::Device::OutputMode::PARAM;
  else if (tmpArg == "-info")
    format = Xyce::Device::OutputMode::INFO;
  else if (tmpArg == "-doc")
    format = Xyce::Device::OutputMode::DOC;
  else if (tmpArg == "-doc_cat")
    format = Xyce::Device::OutputMode::DOC_CAT;

  std::string option_device_name;
  int option_device_level = -1;
  bool print_model = true;
  bool print_instance=true;

  if (i < ia - 1 && ca[i + 1][0] != '-')
  {
    option_device_name = ca[i + 1];
    if (option_device_name.size() > 1 && tolower(option_device_name[0]) != 'y')
    {
      option_device_name = "";
    }
    else
    {
      if (tolower(option_device_name[0]) == 'y')
        option_device_name = option_device_name.substr(1);
      i++;
      if (i < ia - 1 && ca[i + 1][0] != '-')
      {
        option_device_level = atoi(ca[i + 1]);
        if (option_device_level != 0)
        {
          i++;
          if (i < ia - 1 && ca[i + 1][0] != '-')
          {
            if (tolower(ca[i + 1][0]) == 'm')
              print_instance = false;
            else if (tolower(ca[i + 1][0]) == 'i')
              print_model = false;
            else
              i--;
            i++;
          }
        }
      }
    }
  }

  Xyce::Device::Registry2::KeyVector keys = Xyce::Device::getXyceInstanceRegistry().getDerivedKey(typeid(Xyce::Device::Tom::Signature));

  N_DEV_DeviceMgr *dIntPtr = N_DEV_DeviceMgr::factory (*this);

  typedef map<pair<string, int>, N_DEV_Device *> DeviceMap;

  if (format == Xyce::Device::OutputMode::PARAM) {
    static Xyce::Util::indent_streambuf loutStreambuf(std::cout.rdbuf());
    static std::ostream lout(&loutStreambuf);

    for (Xyce::Device::Registry2::KeyVector::const_iterator it = keys.begin(); it != keys.end(); ++it)
    {
      Xyce::Device::ParametricData<void> &instance_parameters = Xyce::Device::Tom::execute(Xyce::Device::getXyceInstanceRegistry(), (*it))();
      Xyce::Device::ParametricData<void> &model_parameters = Xyce::Device::Tom::execute(Xyce::Device::getXyceModelRegistry(), (*it))();

      Xyce::lout() << (*it).first << " level " << (*it).second << Xyce::Util::push << std::endl
                   <<  "Model" << Xyce::Util::push << std::endl;

      Xyce::lout() << "Parameters" << Xyce::Util::push << std::endl;
      Xyce::Device::outputParameterMap(Xyce::lout(), model_parameters.getMap());
      Xyce::lout() << Xyce::Util::pop << std::endl;

      Xyce::lout() << Xyce::Util::pop << std::endl;

      Xyce::lout() <<  "Instance" << Xyce::Util::push << std::endl;

      Xyce::lout() << "Configuration" << Xyce::Util::push << std::endl;
      Xyce::Device::outputConfiguration(Xyce::lout(), instance_parameters.getConfigTable());
      Xyce::lout() << Xyce::Util::pop << std::endl;

      Xyce::lout() << "Parameters" << Xyce::Util::push << std::endl;
      Xyce::Device::outputParameterMap(Xyce::lout(), instance_parameters.getMap());
      Xyce::lout() << Xyce::Util::pop << std::endl;

      Xyce::lout() << Xyce::Util::pop << std::endl;

      Xyce::lout() << Xyce::Util::pop << std::endl;
    }
  }
  else {
    for (Xyce::Device::Registry2::KeyVector::const_iterator it = keys.begin(); it != keys.end(); ++it)
    {
      const Xyce::Device::DeviceLevelKey &device_key = (*it);

      std::string device_name = (*it).first;
      const int device_level = (*it).second;

      device_name[0] = toupper(device_name[0]);

      if ((option_device_name.empty() || Xyce::equal_nocase(option_device_name, device_name)) && (option_device_level == -1 || option_device_level == device_level)) {
        Xyce::Device::SolverState solState;
        Xyce::Device::DeviceOptions devOptions(*this);

        Xyce::Device::Device *device = Xyce::Device::Bob::create(Xyce::Device::getXyceRegistry2(), device_key)(solState, devOptions);
        std::string device_description = device->getName();

        Xyce::Device::ParametricData<void> &instance_parameters = Xyce::Device::Tom::execute(Xyce::Device::getXyceInstanceRegistry(), device_key)();
        Xyce::Device::ParametricData<void> &model_parameters = Xyce::Device::Tom::execute(Xyce::Device::getXyceModelRegistry(), device_key)();

        if (device->getName() == "Behavioral Digital")
          device_name = "Digital";
        
        if (print_instance && !instance_parameters.getMap().empty()) {
          std::ostringstream path;
          path << device_name << "_" << device_level << "_Device_Instance"
               << (format == Xyce::Device::OutputMode::DOC_CAT ? "_Category" : "")
               << "_Params.tex";

          std::ofstream os(path.str().c_str(), ios_base::out);

          laTexDevice(os, device_name, device_level, 0, device_description, instance_parameters, format);
        }

        if (print_model && !model_parameters.getMap().empty()) {
          std::ostringstream path;
          path << device_name << "_" << device_level << "_Device_Model"
               << (format == Xyce::Device::OutputMode::DOC_CAT ? "_Category" : "")
               << "_Params.tex";

          std::ofstream os(path.str().c_str(), ios_base::out);

          laTexDevice(os, device_name, device_level, 1, device_description, model_parameters, format);
        }

        delete device;
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : N_IO_CmdParse::numArgs
// Purpose       : Returns the number of cmd line args.
// Special Notes :
// Scope         : Public
// Creator       : Eric Keiter
// Creation Date : 02/19/01
//-----------------------------------------------------------------------------
int N_IO_CmdParse::numArgs()
{
  return iargs;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_CmdParse::printArgMap
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Eric Keiter
// Creation Date : 02/19/01
//-----------------------------------------------------------------------------
void N_IO_CmdParse::printArgMap()
{
  map<string,string>::iterator iter;
  map<string,string>::iterator begin = stArgs.begin();
  map<string,string>::iterator end   = stArgs.end();

  cout << endl << "Command Line Argument Map:" << endl;
  cout << endl;

  for (iter=begin;iter!=end;++iter)
  {
    cout << "   map[ ";
    cout << (iter->first);
    cout << " ] = ";
    cout << (iter->second) << endl;
  }
  cout << endl;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_CmdParse::argExists
// Purpose       : This function returns true if the specified
//                 argument exists on the command line. It returns
//                 false if either the specified argument does not exist
//                 on the command line or there is no such option.
// Special Notes :
// Scope         : Public
// Creator       : Eric Keiter
// Creation Date : 02/19/01
//-----------------------------------------------------------------------------
bool N_IO_CmdParse::argExists(const string & arg_tmp) const
{
  map<string,int>::const_iterator it = swArgs.find(arg_tmp);
  if (it != swArgs.end() && (*it).second != 0)
    return true;
  else {
    map<string,string>::const_iterator it = stArgs.find(arg_tmp);
    if (it == stArgs.end())
      return false;
    else
      return (*it).second != "";
  }
}

//-----------------------------------------------------------------------------
// Function      : N_IO_CmdParse::getArgumentValue
// Purpose       : This function returns the value of the specified
//                 string argument
// Special Notes :
// Scope         : Public
// Creator       : Lon Waters
// Creation Date : 03/14/2002
//-----------------------------------------------------------------------------
string N_IO_CmdParse::getArgumentValue(const string & argumentName) const
{
  map<string,string>::const_iterator it = stArgs.find(argumentName);

  return it == stArgs.end() ? "" : (*it).second;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_CmdParse::isSwitchArg
// Purpose       : This function returns true if the specified argument
//                 is a switch argument
// Special Notes :
// Scope         : Private
// Creator       : Lon Waters
// Creation Date : 03/13/2002
//-----------------------------------------------------------------------------
bool N_IO_CmdParse::isSwitchArg(string arg)
{
  bool bsuccess = false;

  if ( swArgs.count(arg) ) bsuccess = true;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_CmdParse::isStringValuedArg
// Purpose       : This function returns true if the specified argument
//                 is a string valued argument
// Special Notes :
// Scope         : Private
// Creator       : Lon Waters
// Creation Date : 03/13/2002
//-----------------------------------------------------------------------------
bool N_IO_CmdParse::isStringValuedArg(string arg)
{
  bool bsuccess = false;

  if ( stArgs.count(arg) ) bsuccess = true;

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : N_IO_CmdParse::getArg
// Purpose       :
// Special Notes :
// Scope         : Private
// Creator       : Eric Keiter
// Creation Date : 10/24/2006
//-----------------------------------------------------------------------------
void N_IO_CmdParse::getArg(int i, string & arg)
{
  arg.clear();
  if (cargs[i] != NULL)
  {
    arg = string(cargs[i]);
  }
}

//-----------------------------------------------------------------------------
// Function      : N_IO_CmdParse::getNetlistCopy()
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Keith Santarelli
// Creation Date : 12/10/2007
//-----------------------------------------------------------------------------
bool N_IO_CmdParse::getNetlistCopy()
{
  return netlistCopy_;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_CmdParse::getOneTerm()
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Keith Santarelli
// Creation Date : 12/10/2007
//-----------------------------------------------------------------------------
bool N_IO_CmdParse::getOneTerm()
{
  return oneTerm_;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_CmdParse::getNoDCPath()
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Keith Santarelli
// Creation Date : 12/10/2007
//-----------------------------------------------------------------------------
bool N_IO_CmdParse::getNoDCPath()
{
  return noDCPath_;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_CmdParse::getOneTermRes()
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Keith Santarelli
// Creation Date : 12/10/2007
//-----------------------------------------------------------------------------
string N_IO_CmdParse::getOneTermRes()
{
  return oneTermRes_;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_CmdParse::getNoDCPathRes()
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Keith Santarelli
// Creation Date : 12/10/2007
//-----------------------------------------------------------------------------
string N_IO_CmdParse::getNoDCPathRes()
{
  return noDCPathRes_;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_CmdParse::setNetlistCopy()
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Keith Santarelli
// Creation Date : 12/10/2007
//-----------------------------------------------------------------------------
void N_IO_CmdParse::setNetlistCopy(bool netlistCopyArg)
{
  netlistCopy_ = netlistCopyArg;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_CmdParse::setOneTerm()
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Keith Santarelli
// Creation Date : 12/10/2007
//-----------------------------------------------------------------------------
void N_IO_CmdParse::setOneTerm(bool oneTermArg)
{
  oneTerm_ = oneTermArg;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_CmdParse::setNoDCPath()
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Keith Santarelli
// Creation Date : 12/10/2007
//-----------------------------------------------------------------------------
void N_IO_CmdParse::setNoDCPath(bool noDCPathArg)
{
  noDCPath_ = noDCPathArg;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_CmdParse::setOneTermRes()
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Keith Santarelli
// Creation Date : 12/10/2007
//-----------------------------------------------------------------------------
void N_IO_CmdParse::setOneTermRes(const string & oneTermResArg)
{
  oneTermRes_ = oneTermResArg;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_CmdParse::setNoDCPathRes()
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Keith Santarelli
// Creation Date : 12/10/2007
//-----------------------------------------------------------------------------
void N_IO_CmdParse::setNoDCPathRes(const string & noDCPathResArg)
{
  noDCPathRes_ = noDCPathResArg;
}

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
// Revision Number: $Revision: 1.132.2.5 $
//
// Revision Date  : $Date: 2014/03/12 17:25:19 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>
#include <sstream>
#include <cctype>

#include <N_IO_CmdParse.h>

#include <N_DEV_Configuration.h>
#include <N_DEV_DeviceMgr.h>
#include <N_DEV_Dump.h>
#include <N_DEV_LaTexDoc.h>
#include <N_DEV_RegisterDevices.h>
#include <N_ERH_Message.h>
#include <N_ERH_Progress.h>
#include <N_PDS_Comm.h>
#include <N_PDS_Manager.h>
#include <N_UTL_IndentStreamBuf.h>
#include <N_UTL_LogStream.h>
#include <N_UTL_NoCase.h>
#include <N_UTL_Version.h>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Function      : usage
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : David Baur
// Creation Date : 3/7/2014
//-----------------------------------------------------------------------------
///
/// Print Xyce usage
///
/// @param os   Stream to write usage information to
///
void usage(std::ostream &os)
{
  os << "Usage: Xyce [arguments] netlist\n\n"
     << "Arguments:\n"
     << "  -b                          batch mode flag for spice compatibility (ignored)\n"
     << "  -h                          print usage and exit\n"
     << "  -v                          print version info and exit\n"
     << "  -capabilities               print compiled-in options and exit\n"
     << "  -license                    print license and exit\n"
     << "  -param [device [level [<inst|mod>]]] print a terse summary of model and/or device parameters\n"
     << "  -info [device [level [<inst|mod>]]] output latex tables of model and device parameters to files\n"
     << "  -doc [device [level [<inst|mod>]]] output latex tables of model and device parameters to files\n"
     << "  -doc_cat [device [level [<inst|mod>]]] output latex tables of model and device parameters to files\n"
     << "  -syntax                     check netlist syntax and exit\n"
     << "  -norun                      netlist syntax and topology and exit\n"
     << "  -namesfile                  output internal names file and exit\n"
     << "  -gui                        gui file output\n"
     << "  -jacobian_test              jacobian matrix diagnostic\n"
     << "  -delim <TAB|COMMA|string>   set the output file field delimiter\n"
     << "  -o <path>                   place the results into <path>\n"
     << "  -l <path>                   place the log output into <path>, \"cout\" to log to stdout\n"
     << "  -per-processor              create log file for each procesor, add .<n>.<r> to log path\n"
     << "  -remeasure [existing Xyce output file] recompute .measure() results with existing data\n"
     << "  -nox <on|off>               NOX nonlinear solver usage\n"
     << "  -linsolv <solver>           force usage of specific linear solver\n"
     << "  -maxord <1..5>              maximum time integration  order\n"
     << "  -prf <param file name>      specify a file with simulation parameters\n"
     << "  -rsf <response file name>   specify a file to save simulation responses functions.\n"

#ifndef Xyce_PARALLEL_MPI
     << "  -r <file>                   generate a rawfile named <file> in binary format\n"
     << "  -a                          use with  -r <file>  to output in ascii format\n"
#endif

#ifdef HAVE_DLFCN_H
     << "  -plugin                     load device plugin\n"
#endif

#ifdef Xyce_Dakota
     << "\nDakota build arguments:\n"
     << "  -dakota <dakota input file> dakota input file for this simulation\n"
#endif

#if defined Xyce_DEBUG_NONLINEAR || defined Xyce_DEBUG_TIME || defined Xyce_DEBUG_DEVICE || defined Xyce_TOTALVIEW_BOGON
     << "\nDebug arguments:\n"
#endif

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
}

//-----------------------------------------------------------------------------
// Function      : CmdParse::CmdParse
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Eric Keiter
// Creation Date : 02/14/01
//-----------------------------------------------------------------------------
CmdParse::CmdParse()
: iargs(0),
  cargs(0),
  allocatedCargs_(false)
{
  setCommandArgs();
}

//-----------------------------------------------------------------------------
// Function      : CmdParse::CmdParse
// Purpose       : copy constructor
// Special Notes :
// Scope         : Public
// Creator       : Eric Keiter
// Creation Date : 01/07/07
//-----------------------------------------------------------------------------
CmdParse::CmdParse(const CmdParse & right) :
  iargs(right.iargs),
  swArgs(right.swArgs),
  stArgs(right.stArgs),
  argIndex(right.argIndex),
  parMgr_(right.parMgr_)
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

    std::string tmpString(right.cargs[i]);
    int size = tmpString.size()+2;
    cargs[i] = new char[size];
    for (j=0;j<size;++j) cargs[i][j] = 0;

    sprintf(cargs[i], "%s", tmpString.c_str());
  }
}

//-----------------------------------------------------------------------------
// Function      : CmdParse::CmdParse
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Eric Keiter
// Creation Date : 02/14/01
//-----------------------------------------------------------------------------
CmdParse::~CmdParse()
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
// Function      : CmdParse::setCommandArgs
// Purpose       : Initialize command line arguments
// Special Notes :
// Scope         : Public
// Creator       : Dave Shirley, PSSi
// Creation Date : 01/24/11
//-----------------------------------------------------------------------------
void CmdParse::setCommandArgs()
{
  stArgs.clear();
  swArgs.clear();
  argIndex.clear();

  stArgs[ "netlist" ] = "";     // The input netlist will be a special case.

  // Set the maps containing the expected arguments. There are two types
  // of arguments switch arguments and string valued arguments.
  swArgs[ "-b" ] = 0;
  swArgs[ "-h" ] = 0;           // Help option, print usage and exit.
  swArgs[ "-test" ] = 0;
  swArgs[ "-v" ] = 0;
  swArgs[ "-capabilities" ] = 0;
  swArgs[ "-license" ] = 0;
  swArgs[ "-syntax" ] = 0;
  swArgs[ "-norun" ] = 0;
  swArgs[ "-namesfile" ] = 0;
  swArgs[ "-gui" ] = 0;
  swArgs[ "-jacobian_test" ] = 0;

  stArgs[ "-delim" ] = "";
  stArgs[ "-o" ] = "";
  stArgs[ "-l" ] = "";          // Output log information to a file.
  swArgs[ "-per-processor" ] = 0;
  swArgs[ "-messy-cout" ] = 0;
  stArgs[ "-remeasure" ] = "";  // specify a existing data file on which Xyce will recompute .measure functions.
  stArgs[ "-nox" ] = "";
  stArgs[ "-linsolv" ] = "";
  stArgs[ "-maxord" ] = "";
  stArgs[ "-prf" ] = "";        // specify a parameter input file to set runtime params from a file
  stArgs[ "-rsf" ] = "";        // specify a response output file to save results to a file
  stArgs[ "-r" ] = "";          // Output binary rawfile.
  swArgs[ "-a" ] = 0;           // Use ascii instead of binary in rawfile output

#ifdef HAVE_DLFCN_H
  stArgs[ "-plugin" ] = "";
#endif

#ifdef Xyce_Dakota
  stArgs[ "-dakota" ] = "";     // specify a dakota input file to couple with this simulation
#endif

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
  return;
}

//-----------------------------------------------------------------------------
// Function      : CmdParse::setNetlist
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
//                 CmdParse class has been made.
//
//
// Scope         : Public
// Creator       : Eric Keiter
// Creation Date : 01/07/07
//-----------------------------------------------------------------------------
void CmdParse::setNetlist(const std::string & newNetlist)
{
  int netIndex = 0;

  if ( argIndex.find("netlist") != argIndex.end())
  {
    netIndex = argIndex["netlist"];
  }
  else
  {
    Report::DevelFatal0().in("CmdParse::setNetlist") << "Unable to find netlist argument.";
  }


#if 0
  //++netIndex;

  if ( netIndex >= iargs )
  {
    // Unexectedly ran out of arguments on the command line.
    Report::DevelFatal0().in("CmdParse::setNetlist") << "Did not find previous netlist setting.";
  }
  else if (cargs[netIndex][0] == '-')
  {
    // Error if we ran into another option here.
    Report::DevelFatal0().in("CmdParse::setNetlist") << "Expected option value, but found option " << cargs[netIndex];
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
// Function      : CmdParse::parseCommandLine
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Lon Waters
// Creation Date : 03/13/02
//-----------------------------------------------------------------------------
int CmdParse::parseCommandLine(int ia, char** ca)
{
  iargs = ia;
  cargs = ca;

  N_PDS_Comm & comm = (*parMgr_->getPDSComm());

  bool isSerial = comm.isSerial();
  int procID = comm.procID();

  bool allExit = false;

  setCommandArgs();

  if (procID == 0)
  {
    for (int i = 1; i < ia; ++i)
    {
      // MPICH sometimes creates arguments that get replaced by NULLs by
      // MPI_INIT :-P.  They are always after the arguments that we care about,
      // so stop command line parsing at that point.
      if (ca[i] == NULL) break;

      std::string tmpArg(ca[i]);

      argIndex[tmpArg] = i;

      if (tmpArg[0] == '-')
      {
        if( tmpArg == "-h" )
        {
          usage(Xyce::lout());
          allExit = true;
          break;
        }
        else if( tmpArg == "-v" )
        {
          Xyce::lout() << N_UTL_Version::getFullVersionString() << std::endl;
          allExit = true;
        }
        else if( tmpArg == "-capabilities" )
        {
          Xyce::lout() << N_UTL_Version::getCapabilities() << std::endl;
          allExit = true;
        }
        else if( tmpArg == "-license" )
        {
          Xyce::lout() << N_UTL_Version::getLicense() << std::endl;
          allExit = true;
        }
        else if (tmpArg == "-plugin")
        {
          ++i;

          if ( i >= ia )
          {
            // Unexectedly ran out of arguments on the command line.
            Report::UserError0() << "Did not find required value for option " << ca[i-1];
            continue;
          }
          else if (ca[i][0] == '-')
          {
            // Error if we ran into another option here.
            Report::UserError0() << "Expected option value, but found option " << ca[i];
            continue;
          }
          else
          {
            stArgs[tmpArg] = ca[i];
          }

          const std::string plugin = ca[i];

          for (std::string::size_type i = 0, j = plugin.find_first_of(", "); i != std::string::npos; i = (j == std::string::npos ? j : j + 1), j = plugin.find_first_of(", ", i))
             Device::registerPlugin(plugin.substr(i, j).c_str());
        }
        else if( tmpArg == "-param" || tmpArg == "-info" ||
                 tmpArg == "-doc"   || tmpArg == "-doc_cat" )
        {
          handleParameterOutputs_(tmpArg, i, ia, ca);
          allExit = true;
        }
        else if( tmpArg == "-gui")
        {
           Report::signalProgress("Xyce setting up problem");
           swArgs[tmpArg] = i;
        }
        else if (argExists(tmpArg))
        {
          Xyce::lout() << "More than one \"" << tmpArg << "\" command line argument found, using last one" << std::endl;
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
            Report::UserError0() << "Did not find required value for option " << ca[i-1];
          }
          else if (ca[i][0] == '-')
          {
            // Error if we ran into another option here.
            Report::UserError0() << "Expected option value, but found option " << ca[i];
          }
          else
          {
            stArgs[tmpArg] = ca[i];
          }
        }
        else
        {
          // Invalid option, stop here.
          Xyce::lout() << "Invalid option given: " << tmpArg << std::endl;
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
            Xyce::lout() << "Found second netlist on command line: " << tmpArg << std::endl;
          }
        }
      }
    }

    if( !allExit && ( stArgs["netlist"] == "" ) )
    {
      IO::usage(Xyce::lout());
      Report::UserError0() << "Netlist not found on command line";
    }
  }

  N_ERH_ErrorMgr::safeBarrier(comm.comm());   //  All procs call (1)

  // Push cmds to other procs for parallel
  if( !isSerial )
  {
    int size, length, val, pos;
    int bsize = 500;

    char * buf = new char[bsize];

    if (procID == 0)
    {
      // broadcast termination flag to other procs and shutdown locally
      if (allExit)
      {
        size = -1;
        comm.bcast( &size, 1, 0 );
        return -1;
      }

      size = swArgs.size();
      comm.bcast( &size, 1, 0 );

      std::map<std::string,int>::iterator it_siM = swArgs.begin();
      std::map<std::string,int>::iterator end_siM = swArgs.end();
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

      std::map<std::string,std::string>::iterator it_ssM = stArgs.begin();
      std::map<std::string,std::string>::iterator end_ssM = stArgs.end();
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
      std::string str1, str2;

      comm.bcast( &size, 1, 0 );

      // shutdown nodes individually
      if( -1 == size )
      {
        return -1;
      }

      for( int i = 0; i < size; ++i )
      {
        comm.bcast( &length, 1, 0 );
        comm.bcast( buf, length, 0 );
        str1 = std::string( buf, length );

        comm.bcast( &val, 1, 0 );

        swArgs[str1] = val;
      }

      comm.bcast( &size, 1, 0 );
      for( int i = 0; i < size; ++i )
      {
        comm.bcast( &length, 1, 0 );
        comm.bcast( buf, length, 0 );
        str1 = std::string( buf, length );

        comm.bcast( &length, 1, 0 );
        comm.bcast( buf, length, 0 );
        str2 = std::string( buf, length );

        stArgs[str1] = str2;
      }
    }

    delete [] buf;

#ifdef Xyce_DEBUG_PARALLEL_CMDLINE
    Xyce::dout() << "*************************\n"
                 << "PROC ID: " << procID << std::endl
                 << "Switched Args\n"
                 << "-------------\n";
    std::map<std::string,int>::iterator it_siM = swArgs.begin();
    std::map<std::string,int>::iterator end_siM = swArgs.end();
    for( ; it_siM != end_siM; ++it_siM )
      Xyce::dout() << " " << it_siM->first << ": " << it_siM->second << std::endl;
    Xyce::dout() << std::endl
                 << "String Args\n"
                 << "-------------\n";
    std::map<std::string,std::string>::iterator it_ssM = stArgs.begin();
    std::map<std::string,std::string>::iterator end_ssM = stArgs.end();
    for( ; it_ssM != end_ssM; ++it_ssM )
      Xyce::dout() << " " << it_ssM->first << ": " << it_ssM->second << std::endl;
    Xyce::dout() << "*************************\n";
#endif
  }

  // return -1 for shutdown
  return allExit && isSerial ? -1 : 0;
}

//-----------------------------------------------------------------------------
// Function      : CmdParse::handleParameterOutputs_
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
void CmdParse::handleParameterOutputs_(const std::string &tmpArg, int &i, int &ia, char** ca)
{
  Device::OutputMode::Mode format = Device::OutputMode::DEFAULT;
  if (tmpArg == "-param")
    format = Device::OutputMode::PARAM;
  else if (tmpArg == "-info")
    format = Device::OutputMode::INFO;
  else if (tmpArg == "-doc")
    format = Device::OutputMode::DOC;
  else if (tmpArg == "-doc_cat")
    format = Device::OutputMode::DOC_CAT;

  std::string option_device_name;
  int option_device_level = -1;
  bool print_model = true;
  bool print_instance = true;

  if (i < ia - 1 && ca[i + 1][0] != '-')
  {
    option_device_name = ca[i + 1];
    if (option_device_name.size() > 1 && std::tolower(option_device_name[0]) != 'y')
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

  typedef std::map<std::pair<std::string, int>, Device::Device *> DeviceMap;

  if (format == Device::OutputMode::PARAM) {
    for (Device::Configuration::ConfigurationMap::const_iterator it = Device::Configuration::getConfigurationMap().begin(); it != Device::Configuration::getConfigurationMap().end(); ++it)
    {
      Xyce::dout() << (*it).first << Xyce::Util::push << std::endl
                   << *(*it).second << Xyce::Util::pop << std::endl;
    }
  }
  else {
    for (Device::Configuration::ConfigurationMap::const_iterator it = Device::Configuration::getConfigurationMap().begin(); it != Device::Configuration::getConfigurationMap().end(); ++it)
    {
      const NameLevelKey &device_key = (*it).first;

      std::string device_name = (*it).first.first;
      const int device_level = (*it).first.second;

      device_name[0] = toupper(device_name[0]);

      if ((option_device_name.empty() || Xyce::equal_nocase(option_device_name, device_name)) && (option_device_level == -1 || option_device_level == device_level)) {
        Device::Configuration &configuration = *(*it).second;

        std::string device_description = configuration.getName();
        Device::ParametricData<void> &instance_parameters = configuration.getInstanceParameters();
        Device::ParametricData<void> &model_parameters = configuration.getModelParameters();

        if (configuration.getName() == "Behavioral Digital")
          device_name = "Digital";

        if (print_instance && !instance_parameters.getMap().empty()) {
          std::ostringstream path;
          path << device_name << "_" << device_level << "_Device_Instance"
               << (format == Device::OutputMode::DOC_CAT ? "_Category" : "")
               << "_Params.tex";

          std::ofstream os(path.str().c_str(), std::ios_base::out);
          os << "% This table was generated by Xyce:" << std::endl
             << "%   Xyce " << (format == Device::OutputMode::DOC_CAT ? "-doc_cat " : "-doc ")
             << device_name << " " << device_level << std::endl
             << "%" << std::endl;

          laTexDevice(os, device_name, device_level, 0, device_description, instance_parameters, format);
        }

        if (print_model && !model_parameters.getMap().empty()) {
          std::ostringstream path;
          path << device_name << "_" << device_level << "_Device_Model"
               << (format == Device::OutputMode::DOC_CAT ? "_Category" : "")
               << "_Params.tex";

          std::ofstream os(path.str().c_str(), std::ios_base::out);

          os << "% This table was generated by Xyce:" << std::endl
             << "%   Xyce " << (format == Device::OutputMode::DOC_CAT ? "-doc_cat " : "-doc ")
             << device_name << " " << device_level << std::endl
             << "%" << std::endl;

            laTexDevice(os, device_name, device_level, 1, device_description, model_parameters, format);
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : CmdParse::numArgs
// Purpose       : Returns the number of cmd line args.
// Special Notes :
// Scope         : Public
// Creator       : Eric Keiter
// Creation Date : 02/19/01
//-----------------------------------------------------------------------------
int CmdParse::numArgs()
{
  return iargs;
}

//-----------------------------------------------------------------------------
// Function      : CmdParse::printArgMap
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Eric Keiter
// Creation Date : 02/19/01
//-----------------------------------------------------------------------------
void CmdParse::printArgMap()
{
  std::map<std::string,std::string>::iterator iter;
  std::map<std::string,std::string>::iterator begin = stArgs.begin();
  std::map<std::string,std::string>::iterator end   = stArgs.end();

  Xyce::dout() << std::endl << "Command Line Argument Map:" << std::endl;
  Xyce::dout() << std::endl;

  for (iter=begin;iter!=end;++iter)
  {
    Xyce::dout() << "   map[ ";
    Xyce::dout() << (iter->first);
    Xyce::dout() << " ] = ";
    Xyce::dout() << (iter->second) << std::endl;
  }
  Xyce::dout() << std::endl;
}

//-----------------------------------------------------------------------------
// Function      : CmdParse::argExists
// Purpose       : This function returns true if the specified
//                 argument exists on the command line. It returns
//                 false if either the specified argument does not exist
//                 on the command line or there is no such option.
// Special Notes :
// Scope         : Public
// Creator       : Eric Keiter
// Creation Date : 02/19/01
//-----------------------------------------------------------------------------
bool CmdParse::argExists(const std::string & arg_tmp) const
{
  std::map<std::string,int>::const_iterator it = swArgs.find(arg_tmp);
  if (it != swArgs.end() && (*it).second != 0)
    return true;
  else {
    std::map<std::string,std::string>::const_iterator it = stArgs.find(arg_tmp);
    if (it == stArgs.end())
      return false;
    else
      return (*it).second != "";
  }
}

//-----------------------------------------------------------------------------
// Function      : CmdParse::getArgumentValue
// Purpose       : This function returns the value of the specified
//                 string argument
// Special Notes :
// Scope         : Public
// Creator       : Lon Waters
// Creation Date : 03/14/2002
//-----------------------------------------------------------------------------
std::string CmdParse::getArgumentValue(const std::string & argumentName) const
{
  std::map<std::string,std::string>::const_iterator it = stArgs.find(argumentName);

  return it == stArgs.end() ? "" : (*it).second;
}

//-----------------------------------------------------------------------------
// Function      : CmdParse::isSwitchArg
// Purpose       : This function returns true if the specified argument
//                 is a switch argument
// Special Notes :
// Scope         : Private
// Creator       : Lon Waters
// Creation Date : 03/13/2002
//-----------------------------------------------------------------------------
bool CmdParse::isSwitchArg(const std::string &arg)
{
  bool bsuccess = false;

  if ( swArgs.count(arg) ) bsuccess = true;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : CmdParse::isStringValuedArg
// Purpose       : This function returns true if the specified argument
//                 is a string valued argument
// Special Notes :
// Scope         : Private
// Creator       : Lon Waters
// Creation Date : 03/13/2002
//-----------------------------------------------------------------------------
bool CmdParse::isStringValuedArg(const std::string &arg)
{
  bool bsuccess = false;

  if ( stArgs.count(arg) ) bsuccess = true;

  return bsuccess;
}


//-----------------------------------------------------------------------------
// Function      : CmdParse::getArg
// Purpose       :
// Special Notes :
// Scope         : Private
// Creator       : Eric Keiter
// Creation Date : 10/24/2006
//-----------------------------------------------------------------------------
void CmdParse::getArg(int i, std::string & arg)
{
  arg.clear();
  if (cargs[i] != NULL)
  {
    arg = std::string(cargs[i]);
  }
}

} // namespace IO
} // namespace Xyce

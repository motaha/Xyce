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
// Filename       : $RCSfile: N_IO_DistributionTool.h,v $
//
// Purpose        : Declares the DistributionTool class.  Distribution tool
//                  buffers and distributes circuit blocks (and related data
//                  such as option blocks, metadata, etc) for/during parsing.
//
// Special Notes  :
//
// Creator        : Eric Rankin, SNL
//
// Creation Date  : 03/12/2003
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.37.6.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:42 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------


#ifndef N_IO_DISTRIBUTIONTOOL_H
#define N_IO_DISTRIBUTIONTOOL_H


//---------- Standard Includes ------------------------------------------------
#include <vector>
#include <string>

#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RefCountPtr;
using Teuchos::rcp;

//---------- Xyce Includes ----------------------------------------------------
#include <N_IO_CircuitContext.h>
#include <N_IO_SpiceSeparatedFieldTool.h>


//---------- Forward Declarations ---------------------------------------------
class N_PDS_Comm;
class N_IO_CircuitBlock;
class N_IO_CmdParse;
class N_IO_PkgOptionsMgr;


//-----------------------------------------------------------------------------
// Class          : N_IO_DistributionTool
// Purpose        : Buffers and distributes circuit blocks (and related data
//                  such as option blocks, metadata, etc) for/during parsing.
//-----------------------------------------------------------------------------
class N_IO_DistributionTool
{

public:
  // ctor
  N_IO_DistributionTool( N_IO_CircuitBlock * cktBlk, 
                         N_IO_CmdParse & cp, 
                         N_PDS_Comm * pdsCommPtr );
  
  // dtor
  ~N_IO_DistributionTool();

  // prepare proc 0 to distribute, others receive
  bool start();
  
  // check for name collisions between nodes and devices
  void checkNodeDevConflicts();

  // set total device count and determine device/proc ratio
  void deviceCount( int devices );

  // send circuit context to all procs
  bool circuitContext( N_IO_CircuitContext * const circuitContexts );

   // stage options for tx
  bool circuitOptions( const list<N_UTL_OptionBlock> & options );

void updateCircuitOptions( const list<N_UTL_OptionBlock> & updateOptionList );

  // Send a circuit device line to current proc
  bool circuitDeviceLine(
   vector< N_IO_SpiceSeparatedFieldTool::StringToken > & deviceLine );
  void endDeviceLines();

  // change current subcircuit context after a new subcircuit is started
  bool circuitStart( string const & subcircuitName,
                     list<string> const & nodes,
                     string const & prefix,
                     vector<N_DEV_Param> const & params );

  void setFileName ( string const & fileNameIn );
  
  // change current subcircuit context to previous context
  bool circuitEnd();
  
  // clean up after reaching the end of the entire circuit
  bool done();
  
  // Method to register the package options manager
  bool registerPkgOptionsMgr( RCP<N_IO_PkgOptionsMgr> pkgOptPtr );

  
private:
  // register options  
  bool registerCircuitOptions();
  
 
#ifdef PMPDE  
  // setup distribution mode for MPDE
  void setupMPDE();
#endif


#ifdef Xyce_PARALLEL_MPI

  // pack subcircuit data into buffer; helper for circuitStart()
  bool packSubcircuitData( string const & subcircuitName,
                           list<string> const & nodes,
                           string const & prefix,
                           vector<N_DEV_Param> const & params );
  
  // send options, metatdata, and context to all procs when able
  bool broadcastGlobalData();

  // send buffer from proc 0
  void send(int size = -1);
  
  // receive all data from proc 0
  bool receive();
  
#ifdef PMPDE
  // unpack and parse netlist lines
  void processReceivedLines( vector<char *>& bufs, vector<int>& bufSize );
#endif
  
  // send circuit context blocks to all procs
  void bcastCircuitContext(); 
  
  // receive circuit context from proc 0
  bool receiveCircuitContext();

  // send option blocks to all procs
  void bcastCircuitOptions();  

  // receive option blocks from proc 0
  bool receiveCircuitOptions();

#endif
   

  // parser circuit block
  N_IO_CircuitBlock * cktBlk_;
  
  // package options manager
  RCP<N_IO_PkgOptionsMgr> pkgOptMgrPtr_;
    
  // proc that will receive the next tx
  int currProc_;

#ifdef PMPDE
  // flag MPDE state
  bool usingMPDE_;
#endif  
  
#ifdef Xyce_PARALLEL_MPI  

  // parallel services
  N_PDS_Comm * pdsCommPtr_;
  
  // total number of available procs for distribution 
  int numProcs_;

  // max number of devices to tx to each proc
  int procDeviceCount_;
  
  // total number of devices;
  int devices_;

  // keep track to txmitted devices for the current proc
  int deviceLinesSent_;

  // length of buffer
  int charBufferSize_;

  // length of data packed in buffer
  int charBufferPos_;
  
  // tx/rx buffer
  char * charBuffer_;
  
  // subcircuit data
  vector<string> subcircuitNames_;
  vector< list<string> > subcircuitNodes_;
  vector<string> subcircuitPrefixes_;
  vector< vector<N_DEV_Param> > subcircuitParams_;
  
  // flags for managing global data
  bool circuitContextReady_;

  // flag for adjusting buffered items limit
  bool suppress_bufmgr_;
#endif
  bool circuitOptionsReady_;

  // global data
  N_IO_CircuitContext * circuitContexts_;
  list< N_UTL_OptionBlock > options_;
//  vector< char * > metadataList_;
  vector< string > metadataList_;
  vector< int > metadataSizeList_;

  N_IO_CmdParse & commandLine_;
  string fileName_;
};

#endif //N_IO_DISTRIBUTIONTOOL_H

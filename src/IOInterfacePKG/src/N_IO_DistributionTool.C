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
// Filename       : $RCSfile: N_IO_DistributionTool.C,v $
//
// Purpose        : Defines the DistributionTool class.  Distribution tool
//                  buffers and distributes circuit blocks (and related data
//                  such as option blocks, metadata, etc) for/during parsing.
// Special Notes  :
//
// Creator        : Eric Rankin, SNL
//
// Creation Date  : 03/12/2003
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.88 $
//
// Revision Date  : $Date: 2014/02/24 23:49:22 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>


#include <N_IO_fwd.h>
#include <N_IO_DistributionTool.h>
#include <N_PDS_Comm.h>
#include <N_IO_PkgOptionsMgr.h>
#include <N_IO_CircuitBlock.h>
#include <N_IO_CmdParse.h>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Function      : DistributionTool::DistributionTool
// Purpose       : ctor
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
DistributionTool::DistributionTool( CircuitBlock * cktBlk,
                                              CmdParse & cp,
                                              N_PDS_Comm * pdsCommPtr = NULL )
: cktBlk_( cktBlk ),
  currProc_(0),
  commandLine_(cp)

#ifdef PMPDE
  ,
  usingMPDE_(false)
#endif

#ifdef Xyce_PARALLEL_MPI

  ,
  pdsCommPtr_( pdsCommPtr ),
  charBufferSize_(0),
  charBufferPos_(0),
  charBuffer_(0)

#endif

{ }


//-----------------------------------------------------------------------------
// Function      : DistributionTool::~DistributionTool
// Purpose       : dtor
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
DistributionTool::~DistributionTool()
{

}

//-----------------------------------------------------------------------------
// Function      : DistributionTool::registerPkgOptionsMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, 1437
// Creation Date : 10/21/08
//-----------------------------------------------------------------------------
bool DistributionTool::registerPkgOptionsMgr( PkgOptionsMgr *pkgOptPtr )
{
  pkgOptMgrPtr_ = pkgOptPtr;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DistributionTool::deviceCount
// Purpose       : set total device count and determine device/proc ratio
// Special Notes : a noop in serial
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void DistributionTool::deviceCount( int devices )
{

#ifdef Xyce_PARALLEL_MPI

  // determine how many devices to send to each far proc
  procDeviceCount_ = devices / numProcs_;

#ifdef PMPDE
  // save number of devices
  devices_ = devices;
#endif

#endif

}


//-----------------------------------------------------------------------------
// Function      : DistributionTool::circuitContext
// Purpose       : send circuit context to all procs
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool DistributionTool::circuitContext( CircuitContext * const
 circuitContexts )
{
  // store until explicit send
  circuitContexts_ = circuitContexts;

#ifdef Xyce_PARALLEL_MPI

  // resize buffer if necessary
  int tmpSize = circuitContexts_->packedByteCount() + sizeof( int );

  if ( tmpSize > charBufferSize_ )
  {
    charBufferSize_ = tmpSize;
  }

  // set flag to tx global data
  circuitContextReady_ = true;

  // try to broadcast
  return broadcastGlobalData();

#else

  // FIXME:  setting of circuit context here (on proc 0)
  // FIXME:  to be consistent, the disttool would do this right now
  //cktBlk_->receiveCircuitContext( *circuitContexts_ );
  return true;

#endif

}


#ifdef Xyce_PARALLEL_MPI

//-----------------------------------------------------------------------------
// Function      : DistributionTool::bcastCircuitContext
// Purpose       : send circuit context to all procs
// Special Notes : KRS, 12/14/07:  We're augmenting this to broadcast
//                 info in commandLine_ to all procs (info related to creating
//                 netlists that have resistors from "dangling" nodes to
//                 ground)
// Scope         : private
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void DistributionTool::bcastCircuitContext()
{
  int bsize = 0;
  int pos = 0;

  // pack circuit contexts
  circuitContexts_->pack( charBuffer_, charBufferSize_, pos, pdsCommPtr_ );

  bsize=pos;

  // broadcast
  pdsCommPtr_->bcast( &bsize, 1, 0 );
  pdsCommPtr_->bcast( charBuffer_, bsize, 0 );



  // transform booleans into integers and transmit them
  bool netlistcopy = commandLine_.getHangingResistor().getNetlistCopy();
  bool oneterm = commandLine_.getHangingResistor().getOneTerm();
  bool nodcpath = commandLine_.getHangingResistor().getNoDCPath();

  int nlcopyint = 0;
  int otint = 0;
  int nodcint = 0;

  if (netlistcopy)
    nlcopyint = 1;
  if (oneterm)
    otint = 1;
  if (nodcpath)
    nodcint = 1;

  pdsCommPtr_->bcast(&nlcopyint,1,0);
  pdsCommPtr_->bcast(&otint,1,0);
  pdsCommPtr_->bcast(&nodcint,1,0);

  std::string onetermres(commandLine_.getHangingResistor().getOneTermRes());
  std::string nodcpathres(commandLine_.getHangingResistor().getNoDCPathRes());
  int otreslength = onetermres.size();
  int nodcreslength = nodcpathres.size();

  pdsCommPtr_->bcast(&otreslength,1,0);
  pdsCommPtr_->bcast(&onetermres[0],otreslength,0);
  pdsCommPtr_->bcast(&nodcreslength,1,0);
  pdsCommPtr_->bcast(&nodcpathres[0],nodcreslength,0);

}


//-----------------------------------------------------------------------------
// Function      : DistributionTool::receiveCircuitContext
// Purpose       : receive circuit context from proc 0
// Special Notes : KRS, 12/14/07:  We're augmenting this to receive info in
//                 commandLine_ from proc 0 (info related to creating
//                 netlists that have resistors from "dangling" nodes to
//                 ground) and set those same attributes in commandLine_ on
//                 every other processor.
// Scope         : private
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool DistributionTool::receiveCircuitContext()
{
  int bsize = 0;
  int pos = 0;

  // receive contexts
  pdsCommPtr_->bcast( &bsize, 1, 0 );
  pdsCommPtr_->bcast( charBuffer_, bsize, 0 );

  // unpack circuit contexts; set top level parent to NULL
  circuitContexts_ = cktBlk_->getCircuitContextPtr();

  circuitContexts_->setParentContextPtr( NULL );
  circuitContexts_->unpack( charBuffer_, bsize, pos, pdsCommPtr_ );

  // parse circuit context on far end
  cktBlk_->receiveCircuitContext( *circuitContexts_ );


  // Initialize integer-ized booleans
  int nlcopyint = 0;
  int otint = 0;
  int nodcint = 0;

  // Receive integer-ized booleans
  pdsCommPtr_->bcast(&nlcopyint,1,0);
  pdsCommPtr_->bcast(&otint,1,0);
  pdsCommPtr_->bcast(&nodcint,1,0);

  //Set the appropriate booleans in commandLine_
  if (nlcopyint == 1)
    commandLine_.getHangingResistor().setNetlistCopy(true);
  if (otint == 1)
    commandLine_.getHangingResistor().setOneTerm(true);
  if (nodcint == 1)
    commandLine_.getHangingResistor().setNoDCPath(true);

  //Receive the strings
  int otreslength=0;
  pdsCommPtr_->bcast(&otreslength,1,0);
  std::string onetermres;
  onetermres.resize(otreslength);
  pdsCommPtr_->bcast(&onetermres[0],otreslength,0);
  commandLine_.getHangingResistor().setOneTermRes(onetermres);

  int nodcreslength=0;
  pdsCommPtr_->bcast(&nodcreslength,1,0);
  std::string nodcpathres;
  nodcpathres.resize(nodcreslength);
  pdsCommPtr_->bcast(&nodcpathres[0],nodcreslength,0);
  commandLine_.getHangingResistor().setNoDCPathRes(nodcpathres);

return true;
}

#endif


//-----------------------------------------------------------------------------
// Function      : DistributionTool::updateCircuitOptions
// Purpose       :
// Special Notes :
// Scope         :
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void DistributionTool::updateCircuitOptions(
 const std::list< Util::OptionBlock > & updatedOptionList )
{
  // clear old option block?
  options_.clear();

  // set the option block to the the new value and register
  options_ = updatedOptionList;
// DEBUG_ELR: do this later?
  registerCircuitOptions();

  // set flags to trigger distribution of update later
  circuitOptionsReady_ = false;

}

//-----------------------------------------------------------------------------
// Function      : DistributionTool::circuitOptions
// Purpose       : stage options for tx
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool DistributionTool::circuitOptions( const std::list< Util::OptionBlock > &
  options )
{
  // store until explicit send
  options_ = options;

#ifdef Xyce_PARALLEL_MPI

  // resize buffer if necessary
  int tmpSize = sizeof( int );
  std::list< Util::OptionBlock >::const_iterator it_obL = options_.begin();
  std::list< Util::OptionBlock >::const_iterator it_oeL = options_.end();

  for( ; it_obL != it_oeL; ++it_obL )
  {
    tmpSize += it_obL->packedByteCount();
  }

  if ( tmpSize > charBufferSize_ )
  {
    charBufferSize_ = tmpSize;
  }

  // set flag to tx global data
  circuitOptionsReady_ = true;

  // try to broadcast
  broadcastGlobalData();

#endif

  // register options on proc 0
  return registerCircuitOptions();
}


#ifdef Xyce_PARALLEL_MPI

//-----------------------------------------------------------------------------
// Function      : DistributionTool::bcastCircuitOptions
// Purpose       : send option blocks to all procs
// Special Notes :
// Scope         : private
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void DistributionTool::bcastCircuitOptions()
{
  int bsize = 0;
  int pos = 0;

  std::list< Util::OptionBlock >::iterator it_obL = options_.begin();
  std::list< Util::OptionBlock >::iterator it_oeL = options_.end();

  // pack options
  int count = options_.size();
  pdsCommPtr_->pack( &count, 1, charBuffer_, charBufferSize_, pos );
  for( ; it_obL != it_oeL; ++it_obL )
  {
    it_obL->pack( charBuffer_, charBufferSize_, pos, pdsCommPtr_ );
  }

  bsize = pos;

  // broadcast options
  pdsCommPtr_->bcast( &bsize, 1, 0 );
  pdsCommPtr_->bcast( charBuffer_, bsize, 0 );

}


//-----------------------------------------------------------------------------
// Function      : DistributionTool::receiveCircuitOptions
// Purpose       : receive option blocks from proc 0
// Special Notes :
// Scope         :
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool DistributionTool::receiveCircuitOptions()
{
  int bsize = 0;
  int pos = 0;
  int size;

  // receive from proc 0
  pdsCommPtr_->bcast( &bsize, 1, 0 );
  pdsCommPtr_->bcast( charBuffer_, bsize, 0 );

  // unpack options
  pdsCommPtr_->unpack( charBuffer_, bsize, pos, &size, 1 );
  for( int i = 0; i < size; ++i )
  {
    Util::OptionBlock anOptionBlock;
    anOptionBlock.unpack( charBuffer_, bsize, pos, pdsCommPtr_ );
    options_.push_back( anOptionBlock );
  }

  // register circuit options (far proc)
  return registerCircuitOptions();

}


#endif


//-----------------------------------------------------------------------------
// Function      : DistributionTool::registerCircuitOptions
// Purpose       : register options
// Special Notes :
// Scope         : private
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool DistributionTool::registerCircuitOptions()
{

  std::string netListFile("");
  if (commandLine_.getArgumentValue("netlist") != "")
  {
    netListFile = commandLine_.getArgumentValue("netlist");
  }

  // validate the ref count pointer we're going to use
  if( pkgOptMgrPtr_ == 0 )
  {
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL,
	    "DistributionTool::registerOptions() called with null PkgOptionsMgr\n");
  }

#if 1

  std::list<N_UTL_OptionBlock>::iterator iterOB = options_.begin();
  std::list<N_UTL_OptionBlock>::iterator endOB = options_.end();

  //register options with pkg option mgr (new options control technique)
  for( ; iterOB != endOB; ++iterOB )
  {

#ifdef PMPDE
    // look for .options mpde flag and set distribution mode
    if( !usingMPDE_ && ( "MPDE" == iterOB->name ) )
    {
      setupMPDE();
    }
#endif

    if( (iterOB->getName() == "GLOBAL") && (iterOB->getStatus() != Xyce::Util::PROCESSED_STATE) )
    {
      pkgOptMgrPtr_->submitOptions( *iterOB, netListFile );
      iterOB->setStatus(Xyce::Util::PROCESSED_STATE);
    }
  }

  for( iterOB = options_.begin(); iterOB != endOB; ++iterOB )
  {
    if( iterOB->getName() != "GLOBAL" && iterOB->getName() != "PARALLEL" &&
     iterOB->getName() != "PRINT" )
    {
      if(iterOB->getStatus() != Xyce::Util::PROCESSED_STATE)
      {
        pkgOptMgrPtr_->submitOptions( *iterOB, netListFile );
        iterOB->setStatus(Xyce::Util::PROCESSED_STATE);
      }
    }
  }

  for( iterOB = options_.begin(); iterOB != endOB; ++iterOB )
  {
    if( (iterOB->getName() == "PRINT") && (iterOB->getStatus() != Xyce::Util::PROCESSED_STATE))
    {
      pkgOptMgrPtr_->submitOptions( *iterOB, netListFile );
      iterOB->setStatus(Xyce::Util::PROCESSED_STATE);
    }
  }

#else

  // FIXME this was incorrect... see DistribMgr

  std::list< Util::OptionBlock >::const_iterator it_obL = options_.begin();
  std::list< Util::OptionBlock >::const_iterator it_oeL = options_.end();

  // register options with pkg option mgr (the new options control technique)
  for( ; it_obL != it_oeL; ++it_obL )
  {
    // original logic:  GLOBAL or (!GLOBAL and !PARALLEL and !PRINT) or PRINT
    if ( ( it_obL->name == "GLOBAL" ) || ( it_obL->name == "PRINT" ) ||
       ( it_obL->name != "PARALLEL" ) )
    {
      pkgOptMgrPtr_->submitOptions( *it_obL, netListFile );
    }
  }

#endif

  // FIXME this should be in some other class (CB) for consistency
  return true;
}

//-----------------------------------------------------------------------------
// Function      : DistributionTool::circuitDeviceLine
// Purpose       : Send a circuit device line to current proc
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool DistributionTool::circuitDeviceLine(std::vector< SpiceSeparatedFieldTool::StringToken > & deviceLine )
{

#ifdef Xyce_PARALLEL_MPI

  if( currProc_ != 0 )
  {
    int size = 0;

    // count line type
    size += sizeof(char);

    // count line size
    size += sizeof(int);

    // count device line
    int deviceLineSize = deviceLine.size();
    for (int i = 0; i < deviceLineSize; ++i)
    {
      size += deviceLine[i].packedByteCount();
    }

    // flush buffer as needed
    send(size);

    // pack line type; "d"evice line
    char lineType = 'd';
    pdsCommPtr_->pack( &lineType, 1, charBuffer_, charBufferSize_, charBufferPos_ );

    // pack device line size
    pdsCommPtr_->pack( &deviceLineSize, 1, charBuffer_, charBufferSize_, charBufferPos_ );

    // pack device line
    for (int i = 0; i < deviceLineSize; ++i)
    {
      deviceLine[ i ].pack( charBuffer_, charBufferSize_, charBufferPos_, pdsCommPtr_ );
    }

    // increment device line counter
    deviceLinesSent_++;

    // manage proc change if n/p reached
    if (deviceLinesSent_ >= procDeviceCount_)
    {
      int minus1 = -1;

      // flush buffer
      send();

#ifdef PMPDE

      if( !usingMPDE_ )
      {

#endif

        // allow this proc to exit receive loop and begin processing
        pdsCommPtr_->send( &minus1, 1, currProc_ );

#ifdef PMPDE

      }

      else
      {
        // broadcast to all procs to exit receive loop and begin processing
        pdsCommPtr_->bcast( &minus1, 1, 0 );

        // proc 0 processes last inconsequential lines (e.g. comments, .END);
        currProc_ = 0;
        return false;
      }

#endif

      // reset device line counter
      deviceLinesSent_ = 0;

      // move to next processor
      currProc_++;

      // switch to proc 0 if nearly complete
      if (currProc_ == numProcs_)
      {
        currProc_ = 0;
      }

      // otherwise, prepare currProc for next set of device lines
      else
      {
        int length = fileName_.size();
        lineType = 'f';

        // flush buffer as needed
        send(sizeof(char) + sizeof(int) + length);

        // pack the filename
        pdsCommPtr_->pack( &lineType, 1, charBuffer_, charBufferSize_, charBufferPos_ );
        pdsCommPtr_->pack( &length, 1, charBuffer_, charBufferSize_, charBufferPos_ );
        pdsCommPtr_->pack( fileName_.c_str(), length, charBuffer_, charBufferSize_, charBufferPos_ );

        // pack the subcircuit names
        int subcircuitNamesSize = subcircuitNames_.size();
        for (int i = 0; i < subcircuitNamesSize; ++i)
        {
          if( currProc_ != 0 )
          {
            packSubcircuitData(subcircuitNames_[i], subcircuitNodes_[i],
                               subcircuitPrefixes_[i], subcircuitParams_[i]);
          }
        }
      }
    }

#ifdef PMPDE

    // mpde mode tell proc 0 to process each line
    return (!usingMPDE_);

#endif

    return true;
  }

  else

#endif

    // let proc 0 parse line
    return false;
}

//-----------------------------------------------------------------------------
// Function      : DistributionTool::endDeviceLines()
// Purpose       : Make sure other processors have been told that device processing
//                 is ended.  This can be a problem with small device count circuits
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 03/06/06
//-----------------------------------------------------------------------------
void DistributionTool::endDeviceLines()
{

#ifdef Xyce_PARALLEL_MPI

  if (currProc_ > 0)
  {
    int minus1 = -1;

    // flush buffer
    send();

    // end parsing on current node
    pdsCommPtr_->send( &minus1, 1, currProc_ );
    ++currProc_;

    // end parsing on remaining nodes
    for ( ; currProc_ < numProcs_ ; ++currProc_)
    {
      pdsCommPtr_->send( &minus1, 1, currProc_ );
    }

    // complete any remaining parser action on node 0
    currProc_ = 0;
  }

#endif

}

//-----------------------------------------------------------------------------
// Function      : DistributionTool::circuitStart
// Purpose       : change current subcircuit context after a new subcircuit is
//               : started
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool DistributionTool::circuitStart( std::string const & subcircuitName,
                                          std::list<std::string> const & nodes,
					  std::string const & prefix,
                                          std::vector<N_DEV_Param> const & params )
{

#ifdef Xyce_PARALLEL_MPI

  if( currProc_ != 0 )
  {
    // save new context data on stacks
    subcircuitNames_.push_back( subcircuitName );
    subcircuitPrefixes_.push_back( prefix );
    subcircuitNodes_.push_back( nodes );
    subcircuitParams_.push_back( params );

    // pack data into buffer
    packSubcircuitData( subcircuitName, nodes, prefix, params );
  }

#endif

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DistributionTool::setFileName
// Purpose       : Change name of netlist file
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 03/10/06
//-----------------------------------------------------------------------------
void DistributionTool::setFileName(std::string const & fileNameIn)
{

#ifdef Xyce_PARALLEL_MPI

  fileName_ = fileNameIn;

  char lineType = 'f';
  int length = fileName_.size();

  // flush buffer as needed
  send(sizeof(char) + sizeof(int) + length);

  pdsCommPtr_->pack( &lineType, 1, charBuffer_, charBufferSize_, charBufferPos_ );
  pdsCommPtr_->pack( &length, 1, charBuffer_, charBufferSize_, charBufferPos_ );
  pdsCommPtr_->pack( fileName_.c_str(), length, charBuffer_, charBufferSize_, charBufferPos_ );

#endif

}

#ifdef Xyce_PARALLEL_MPI


//-----------------------------------------------------------------------------
// Function      : DistributionTool::packSubcircuitData
// Purpose       : pack subcircuit data into buffer; helper for circuitStart()
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool DistributionTool::packSubcircuitData( std::string const & subcircuitName,
                                                std::list<std::string> const & nodes,
                                                std::string const & prefix,
                                                std::vector<N_DEV_Param> const & params )
{
  int size = 0;

  // count line type, name size, and name chars
  size += sizeof(char) + sizeof(int) + subcircuitName.length();

  // count # of nodes,
  size += sizeof(int);

  // count node sizes and node chars
  std::list<std::string>::const_iterator it_sbL = nodes.begin();
  std::list<std::string>::const_iterator it_seL = nodes.end();
  for( ; it_sbL != it_seL; ++it_sbL  )
  {
    size += sizeof(int) + it_sbL->length();
  }

  // count prefix size and prefix chars
  size += sizeof(int) + prefix.length();

  // count # of params
  size += sizeof(int);

  // count params
  std::vector<N_DEV_Param>::const_iterator it_pbL = params.begin();
  std::vector<N_DEV_Param>::const_iterator it_peL = params.end();
  for( ; it_pbL != it_peL; ++it_pbL  )
  {
    size += it_pbL->packedByteCount();
  }

  // flush buffer as needed
  send(size);


  // pack line type; subcircuit "s"tart
  char lineType = 's';
  pdsCommPtr_->pack( &lineType, 1, charBuffer_, charBufferSize_, charBufferPos_ );

  // pack name size and name chars
  int length = subcircuitName.length();
  pdsCommPtr_->pack( &length, 1, charBuffer_, charBufferSize_, charBufferPos_ );
  pdsCommPtr_->pack( subcircuitName.c_str(), length, charBuffer_, charBufferSize_, charBufferPos_ );

  // pack # of nodes
  size = nodes.size();
  pdsCommPtr_->pack( &size, 1, charBuffer_, charBufferSize_, charBufferPos_ );

  // pack nodes name sizes and chars
  it_sbL = nodes.begin();
  for( ; it_sbL != it_seL; ++it_sbL  )
  {
    length = it_sbL->length();
    pdsCommPtr_->pack( &length, 1, charBuffer_, charBufferSize_, charBufferPos_ );
    pdsCommPtr_->pack( it_sbL->c_str(), length, charBuffer_, charBufferSize_, charBufferPos_ );
  }

  // pack prefix size and prefix chars
  length = prefix.length();
  pdsCommPtr_->pack( &length, 1, charBuffer_, charBufferSize_, charBufferPos_ );
  pdsCommPtr_->pack( prefix.c_str(), length, charBuffer_, charBufferSize_, charBufferPos_ );

  // pack # of params
  size = params.size();
  pdsCommPtr_->pack( &size, 1, charBuffer_, charBufferSize_, charBufferPos_ );

  // pack params
  it_pbL = params.begin();
  for( ; it_pbL != it_peL; ++it_pbL  )
  {
    it_pbL->pack( charBuffer_, charBufferSize_, charBufferPos_, pdsCommPtr_ );
  }

  return true;
}

#endif


//-----------------------------------------------------------------------------
// Function      : DistributionTool::circuitEnd
// Purpose       : change current subcircuit context to previous context
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool DistributionTool::circuitEnd()
{

#ifdef Xyce_PARALLEL_MPI

  if( currProc_ != 0 )
  {
    // remove latest context data on stacks
    subcircuitNames_.pop_back();
    subcircuitPrefixes_.pop_back();
    subcircuitNodes_.pop_back();
    subcircuitParams_.pop_back();

    // flush buffer as needed
    send(sizeof(char));

    // pack line type; subcircuit "e"nd
    char lineType = 'e';
    pdsCommPtr_->pack( &lineType, 1, charBuffer_, charBufferSize_, charBufferPos_ );
  }

#endif

  return true;
}


//-----------------------------------------------------------------------------
// Function      : DistributionTool::done
// Purpose       : clean up after reaching the end of the entire circuit
// Special Notes : far procs get cleaned up after switch to next proc
//               : while they are in receive loop, not here
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool DistributionTool::done()
{

#ifdef Xyce_PARALLEL_MPI

  int isUpdateRequired = ( circuitOptionsReady_ == false );

  // send update flag
  pdsCommPtr_->bcast( &isUpdateRequired, 1, 0 );

  if ( isUpdateRequired )
  {
    if ( pdsCommPtr_->procID() == 0 )
    {
      bcastCircuitOptions();
    }
    else
    {
      options_.clear();
      receiveCircuitOptions();
    }
  }

  if( charBufferSize_ > 0 )
  {
    delete [] charBuffer_;
  }

#endif
  return true;
}


#ifdef Xyce_PARALLEL_MPI

//-----------------------------------------------------------------------------
// Function      : DistributionTool::send
// Purpose       : send buffer from proc 0 and adjust buffer limits
// Special Notes : size == -1, and no args send(), will force a transmission
// Scope         :
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
void DistributionTool::send(int size)
{

  // transmit buffer if next set of object will fill it completely, or forced
  if ((charBufferPos_ + size >= charBufferSize_) || (size == -1))
  {

#ifdef Xyce_DEBUG_DISTRIBUTED_PARSER

    cerr << "node " << pdsCommPtr_->procID() << " packed "
         << charBufferPos_ << " bytes into " << charBufferSize_ << " byte buffer"
         << " [data size:  " << size << "]" << endl
         << "node " << currProc_ << " was sent "
         << deviceLinesSent_ << " of " << procDeviceCount_ << " devices "
         << (float)deviceLinesSent_ / procDeviceCount_ * 100.0
         << ( size == -1 ? "% complete  [flushed]" : "% complete"  ) << endl;

#endif

#ifdef PMPDE

    if( !usingMPDE_ )
    {

#endif

      // tx buffer size
      pdsCommPtr_->send( &charBufferPos_, 1, currProc_ );

      // tx buffer contents
      pdsCommPtr_->send( charBuffer_, charBufferPos_, currProc_ );

#ifdef PMPDE

    }

    // broadcast all communication for MPDE
    else
    {
      pdsCommPtr_->bcast( &charBufferPos_, 1, 0 );
      pdsCommPtr_->bcast( charBuffer_, charBufferPos_, 0 );
    }

#endif

    // reset counters
    charBufferPos_ = 0;

    // adjust buffer size if objects larger than existing buffer
    if (size > charBufferSize_)
    {

      // free memory
      // delete [] charBuffer_;

      // set new buffer size
      charBufferSize_ = size;

      // reserve memory; including offset
      charBuffer_ = new char[charBufferSize_ + sizeof(char) + sizeof(int)];

#ifdef Xyce_DEBUG_DISTRIBUTED_PARSER

      cerr << "node " << pdsCommPtr_->procID()
           << " resized buffer to " << charBufferSize_
           << " with offset " << sizeof(char) + sizeof(int) << endl;

#endif

    }
  }
}


//-----------------------------------------------------------------------------
// Function      : DistributionTool::receive
// Purpose       : receive all data from proc 0
// Special Notes :
// Scope         :
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool DistributionTool::receive()
{
  bool processingLine = true;
  int size, bsize, length, pos, i;
  char lineType;

  N_ERH_ErrorMgr::safeBarrier(pdsCommPtr_->comm());     // Slave procs call (3)

  // receive buffer size
  pdsCommPtr_->bcast( &charBufferSize_, 1, 0 );

  // check for halt request
  if (charBufferSize_ < 0)
  {
    return false;
  }

  // reserve memory; including offset
  charBuffer_ = new char[charBufferSize_ + sizeof(char) + sizeof(int)];

  // receive global data

  if (!receiveCircuitOptions())
  {
    std::string msg("Internal error registering/receiving circuit options");
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, msg );
  }
  else
  {
    N_ERH_ErrorMgr::safeBarrier(pdsCommPtr_->comm());     // Slave procs call (4)
  }

  receiveCircuitContext();

  pdsCommPtr_->barrier();

  // free memory
  // delete [] charBuffer_;

  std::vector<char *> bufs;
  std::vector<int> bufSize;
  char *currBuffer_;
  unsigned int ib;

#ifdef PMPDE
  if( !usingMPDE_ )
  {
#endif

    // receive the device lines first
    while (true)
    {
      pdsCommPtr_->recv( &bsize, 1, 0 );

      // check for halt request
      if (bsize < 0)
      {
        break;
      }

      currBuffer_ = new char[bsize];
      bufs.push_back(currBuffer_);
      bufSize.push_back(bsize);

      pdsCommPtr_->recv(currBuffer_, bsize, 0);
    }

#ifdef PMPDE
    // parse the device lines second
    processReceivedLines( bufs, bufSize );
  }

  else
  {

    // interleave broadcasts with parsing
    while (true)
    {
      pdsCommPtr_->bcast( &bsize, 1, 0 );

      // check for halt request
      if (bsize < 0)
      {
        break;
      }

      currBuffer_ = new char[bsize];
      bufs.push_back(currBuffer_);
      bufSize.push_back(bsize);

      pdsCommPtr_->bcast( currBuffer_, bsize, 0 );

      processReceivedLines( bufs, bufSize );
    }
  }
#endif


#ifndef PMPDE
  // process received device lines, contexts info, and exit calls
  for (ib=0 ; ib<bufs.size() ; ++ib)
  {
    currBuffer_ = bufs[ib];
    bsize = bufSize[ib];

    // get ready to process buffer
    pos = 0;
    processingLine = true;

    // break apart buffered lines and parse
    while( pos < bsize )
    {
      // get the linetype marker that indicates incoming line
      pdsCommPtr_->unpack( currBuffer_, bsize, pos, &lineType, 1 );

      // process line according to type
      switch( lineType )
      {
        case 'd': // process a device line
        {
          // unpack the device line
          pdsCommPtr_->unpack( currBuffer_, bsize, pos, &size, 1 );
          std::vector< SpiceSeparatedFieldTool::StringToken > deviceLine( size );

          for( i = 0; i < size; ++i )
          {
            deviceLine[i].unpack( currBuffer_, bsize, pos, pdsCommPtr_ );
          }

          // hand to circuit block for processing
          cktBlk_->handleDeviceLine( deviceLine );

          break;
        }

        case 's': // subcircuit start found
        {
          // get the subcircuit name
          pdsCommPtr_->unpack( currBuffer_, bsize, pos, &length, 1 );
          std::string subcircuitName(std::string( ( currBuffer_ + pos ), length ));
          pos += length;

          // get the nodes for this subcircuit call
          std::list<std::string> nodes;
          pdsCommPtr_->unpack( currBuffer_, bsize, pos, &size, 1 );
          for( i = 0; i < size; ++i )
          {
            pdsCommPtr_->unpack( currBuffer_, bsize, pos, &length, 1 );
            nodes.push_back( std::string( ( currBuffer_ + pos ), length ) );
            pos += length;
          }

          // get the prefix
          pdsCommPtr_->unpack( currBuffer_, bsize, pos, &length, 1 );
          std::string prefix(std::string( ( currBuffer_ + pos ), length ));
          pos += length;

          // get params
          std::vector<N_DEV_Param> params;
          pdsCommPtr_->unpack( currBuffer_, bsize, pos, &size, 1 );
          for( i = 0; i < size; ++i )
          {
            Device::Param param;
            param.unpack( currBuffer_, bsize, pos, pdsCommPtr_ );
            params.push_back( param );
          }

          // send to cktblk
          circuitContexts_->setContext( subcircuitName, prefix, nodes );
          circuitContexts_->resolve( params );

          break;
        }

        case 'e': // subcircuit end found
        {
          // adjust the circuit context pointer
          circuitContexts_->restorePreviousContext();

          break;
        }

        case 'f': // change netlist file name
        {
          pdsCommPtr_->unpack( currBuffer_, bsize, pos, &length, 1 );
          fileName_ = std::string( ( currBuffer_ + pos ), length );
          pos += length;
          cktBlk_->setFileName(fileName_);

          break;
        }

        default:  // something went wrong
        {
          std::string msg("Node ");
          msg += pdsCommPtr_->procID();
          msg += " received an invalid message type: \"";
          msg += lineType;
          msg += "\"\n";
          N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL_0, msg );
        }
      }
    }
    delete [] bufs[ib];
  }
#endif

  bufs.clear();
  bufSize.clear();

  checkNodeDevConflicts();

  return true;
}


#ifdef PMPDE
//-----------------------------------------------------------------------------
// Function      : DistributionTool::processReceivedLLines
// Purpose       : unpack and parse netlist lines
// Special Notes :
// Scope         : private
// Creator       : Eric Rankin
// Creation Date :
//-----------------------------------------------------------------------------
void DistributionTool::processReceivedLines( std::vector<char *> & bufs,
  std::vector<int> & bufSize )
{
  bool processingLine = true;
  int size, bsize, length, pos, i;
  char lineType;

  char *currBuffer_;

  // process received device lines, contexts info, and exit calls
  for (int ib=0 ; ib<bufs.size() ; ++ib)
  {
    currBuffer_ = bufs[ib];
    bsize = bufSize[ib];

    // get ready to process buffer
    pos = 0;
    processingLine = true;

    // break apart buffered lines and parse
    while( pos < bsize )
    {
      // get the linetype marker that indicates incoming line
      pdsCommPtr_->unpack( currBuffer_, bsize, pos, &lineType, 1 );

      // process line according to type
      switch( lineType )
      {
        case 'd': // process a device line
        {
          // unpack the device line
          pdsCommPtr_->unpack( currBuffer_, bsize, pos, &size, 1 );
          std::vector< SpiceSeparatedFieldTool::StringToken > deviceLine( size );

          for( i = 0; i < size; ++i )
          {
            deviceLine[i].unpack( currBuffer_, bsize, pos, pdsCommPtr_ );
          }

          // hand to circuit block for processing
          cktBlk_->handleDeviceLine( deviceLine );

          break;
        }

        case 's': // subcircuit start found
        {
          // get the subcircuit name
          pdsCommPtr_->unpack( currBuffer_, bsize, pos, &length, 1 );
          std::string subcircuitName(std::string( ( currBuffer_ + pos ), length ));
          pos += length;

          // get the nodes for this subcircuit call
          std::list<std::string> nodes;
          pdsCommPtr_->unpack( currBuffer_, bsize, pos, &size, 1 );
          for( i = 0; i < size; ++i )
          {
            pdsCommPtr_->unpack( currBuffer_, bsize, pos, &length, 1 );
            nodes.push_back( std::string( ( currBuffer_ + pos ), length ) );
            pos += length;
          }

          // get the prefix
          pdsCommPtr_->unpack( currBuffer_, bsize, pos, &length, 1 );
          std::string prefix(std::string( ( currBuffer_ + pos ), length ));
          pos += length;

          // get params
          std::vector<Device::Param> params;
          pdsCommPtr_->unpack( currBuffer_, bsize, pos, &size, 1 );
          for( i = 0; i < size; ++i )
          {
            Device::Param param;
            param.unpack( currBuffer_, bsize, pos, pdsCommPtr_ );
            params.push_back( param );
          }

          // send to cktblk
          circuitContexts_->setContext( subcircuitName, prefix, nodes );
          circuitContexts_->resolve( params );

          break;
        }

        case 'e': // subcircuit end found
        {
          // adjust the circuit context pointer
          circuitContexts_->restorePreviousContext();

          break;
        }

        case 'f': // change netlist file name
        {
          pdsCommPtr_->unpack( currBuffer_, bsize, pos, &length, 1 );
          fileName_ = std::string( ( currBuffer_ + pos ), length );
          pos += length;
          cktBlk_->setFileName(fileName_);

          break;
        }

        default:  // something went wrong
        {
          std::string msg("Node ");
          msg += pdsCommPtr_->procID();
          msg += " received an invalid message type: \"";
          msg += lineType;
          msg += "\"\n";
          N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL_0, msg );
        }
      }
    }

    delete bufs[ib];
  }

  bufs.clear();
  bufSize.clear();
}
#endif

#endif

//-----------------------------------------------------------------------------
// Function      : DistributionTool::checkNodeDevConflicts
// Purpose       : Check for name collisions between devices
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/20/06
//-----------------------------------------------------------------------------
void DistributionTool::checkNodeDevConflicts()
{
#ifdef Xyce_PARALLEL_MPI
  numProcs_ = pdsCommPtr_->numProc();

  if ( numProcs_ != 1 )
  {
    std::vector<std::string> deviceNames;
    std::vector<std::string>::iterator n_i, n_end;
    int byteCount, maxByteCount, length, pos, sendProc, recvProc, bsize;
    byteCount = cktBlk_->getDeviceNames(deviceNames);

    N_ERH_ErrorMgr::safeBarrier(pdsCommPtr_->comm());     // All procs call (5)

    pdsCommPtr_->maxAll(&byteCount, &maxByteCount, 1);
    std::vector<char> sendBuffer(byteCount);
    std::vector<char> recvBuffer(maxByteCount);
    n_end = deviceNames.end();
    pos = 0;
    for (n_i = deviceNames.begin() ; n_i != n_end ; ++n_i)
    {
      length = n_i->size();
      pdsCommPtr_->pack (&length, 1, &sendBuffer[0], byteCount, pos);
      pdsCommPtr_->pack (n_i->c_str(), length, &sendBuffer[0], byteCount, pos);
    }

    sendProc = pdsCommPtr_->procID()+1;
    recvProc = pdsCommPtr_->procID()-1;
    if (sendProc == numProcs_) sendProc = 0;
    if (recvProc == -1) recvProc = numProcs_-1;
    while (sendProc != pdsCommPtr_->procID())
    {
      pdsCommPtr_->iRecv( &bsize, 1, recvProc );
      pdsCommPtr_->send( &byteCount, 1, sendProc );
      pdsCommPtr_->waitAll();
      pdsCommPtr_->iRecv( &recvBuffer[0], bsize, recvProc );
      pdsCommPtr_->send( &sendBuffer[0], byteCount, sendProc );
      pdsCommPtr_->waitAll();
      deviceNames.resize(0);
      pos = 0;
      while (pos < bsize)
      {
        pdsCommPtr_->unpack(&recvBuffer[0], bsize, pos, &length, 1);
        deviceNames.push_back(std::string((&recvBuffer[0])+pos, length));
        pos += length;
      }

      cktBlk_->checkDeviceNames(deviceNames);
      N_ERH_ErrorMgr::safeBarrier(pdsCommPtr_->comm());          // All procs call (5+N)
      if (++sendProc == numProcs_) sendProc = 0;
      if (--recvProc == -1) recvProc = numProcs_-1;
    }
  }
#endif
}

//-----------------------------------------------------------------------------
// Function      : DistributionTool::start
// Purpose       : prepare proc 0 to distribute, others receive
// Special Notes :
// Scope         : public
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool DistributionTool::start()
{

#ifdef Xyce_PARALLEL_MPI

  // procID == 0 may not be true if we're running in a hierarchical parallel context
  // so check if numProc() == 1 too.
  if (pdsCommPtr_->procID() == 0)
  {
    numProcs_ = pdsCommPtr_->numProc();
    ( numProcs_ == 1) ? currProc_ = 0 : currProc_ = 1;
    deviceLinesSent_ = 0;
    charBufferSize_ = 250000;
    charBufferPos_ = 0;
    circuitContextReady_ = false;
    circuitOptionsReady_ = false;
  }

  // put far end procs in receive loop
  else
  {
    receive();
  }

#endif

  return true;
}


#ifdef Xyce_PARALLEL_MPI

//-----------------------------------------------------------------------------
// Function      : DistributionTool::broadcastGlobalData
// Purpose       : send bufsize, options, metatdata, and context to procs
// Special Notes :
// Scope         : private
// Creator       :
// Creation Date :
//-----------------------------------------------------------------------------
bool DistributionTool::broadcastGlobalData()
{
  // make sure data is ready to tx
  if( circuitContextReady_ && circuitOptionsReady_ && numProcs_ >= 1 )
  {

    N_ERH_ErrorMgr::safeBarrier(pdsCommPtr_->comm());    // Master proc calls (3)

    // reserve memory; including offset
    charBuffer_ = new char[charBufferSize_ + sizeof(char) + sizeof(int)];

    // tx buffer size
    pdsCommPtr_->bcast( &charBufferSize_, 1, 0 );

    // tx global data
    bcastCircuitOptions();
    N_ERH_ErrorMgr::safeBarrier(pdsCommPtr_->comm());    // Master proc calls (4)
    bcastCircuitContext();

    pdsCommPtr_->barrier();
  }

  return true;
}

#endif


#ifdef PMPDE
//-----------------------------------------------------------------------------
// Function      : DistributionTool::setupMPDE
// Purpose       : setup distribution mode for MPDE
// Special Notes :
// Scope         : private
// Creator       : Eric Rankin
// Creation Date :
//-----------------------------------------------------------------------------
void DistributionTool::setupMPDE()
{
  // flag mpde mode
  usingMPDE_ = true;

#ifdef Xyce_PARALLEL_MPI

  // raise count to send entire circuit to all procs
  procDeviceCount_ = devices_;

#endif

}
#endif

} // namespace IO
} // namespace Xyce

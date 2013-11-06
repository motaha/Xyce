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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_IO_RestartMgr.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 7/19/01
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.53.2.3 $
//
// Revision Date  : $Date: 2013/10/03 17:23:43 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


//------- Std Includes -----------

#include <N_UTL_Misc.h>

#include <sstream>

#include <fstream>

#ifdef HAVE_ALGORITHM
#include <algorithm>
#else
#ifdef HAVE_ALGO_H
#include <algo.h>
#else
error Must have either <algorithm> or <algo.h>!
#endif
#endif

#ifdef HAVE_CSTDIO
#include <cstdio>
#else
#include <stdio.h>
#endif

//------- Xyce Includes ----------

#include <N_IO_RestartMgr.h>
#include <N_IO_RestartNode.h>

#include <N_TOP_Topology.h>

#include <N_DEV_DeviceInterface.h>

#include <N_ANP_AnalysisInterface.h>

#include <N_PDS_Manager.h>
#include <N_PDS_Comm.h>

#include <N_UTL_OptionBlock.h>
#include <N_UTL_Functors.h>

#include <N_IO_PkgOptionsMgr.h>

#include <N_ERH_ErrorMgr.h>

#include <N_IO_CmdParse.h>

//------ Misc Includes --------

//-----------------------------------------------------------------------------
// Function      : N_IO_RestartMgr::factory
// Purpose       : singleton access
// Special Notes : static
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 7/19/01
//-----------------------------------------------------------------------------
N_IO_RestartMgr * N_IO_RestartMgr::factory(N_IO_CmdParse & cl)
{
  N_IO_RestartMgr * rmPtr = new N_IO_RestartMgr(cl);
  return rmPtr;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_RestartMgr::N_IO_RestartMgr
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 7/19/01
//-----------------------------------------------------------------------------
N_IO_RestartMgr::N_IO_RestartMgr(N_IO_CmdParse & cl)
: pdsMgrPtr_(0),
  topMgrPtr_(0),
  devIntPtr_(0),
  anaIntPtr_(0),
  restartFlag_(false),
  restartFileName_(""),
  restartJobName_(""),
  pack_(true),
  commandLine_(cl),
  initialSaveInterval_(0.0)
{
}

//-----------------------------------------------------------------------------
// Function      : N_IO_RestartMgr::registerPkgOptionsMgr
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, 1437
// Creation Date : 10/21/08
//-----------------------------------------------------------------------------
bool N_IO_RestartMgr::registerPkgOptionsMgr( RCP<N_IO_PkgOptionsMgr> pkgOptPtr )
{
  pkgOptMgrPtr_ = pkgOptPtr;
  string netListFile("");
  if (commandLine_.getArgumentValue("netlist") != "")
  {
    netListFile = commandLine_.getArgumentValue("netlist");
  }
  pkgOptMgrPtr_->submitRegistration( "RESTART", netListFile,
                                  new N_IO_RestartMgr_OptionsReg( this ) );
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_RestartMgr::registerRestartOptions
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 7/19/01
//-----------------------------------------------------------------------------
bool N_IO_RestartMgr::registerRestartOptions(const N_UTL_OptionBlock & OB)
{
  list<N_UTL_Param>::const_iterator iterPL = OB.getParams().begin();
  bool done = false;

  //For windoze/mingw, packed restart data blows chunks, so we always turn it off
#ifdef Xyce_RESTART_NOPACK
  pack_ = false;
#endif

      ostringstream ost;

  while( iterPL != OB.getParams().end() )
  {
    ExtendedString currTag ( iterPL->tag() );
    currTag.toUpper();

    if( currTag == "TIME" )
    {
      double t = iterPL->dVal();
      ++iterPL;
      if( iterPL == OB.getParams().end() )
      {
        N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL_0,
         ".OPTIONS Restart: Badly formed INITIAL_INTERVAL\n" );
      }
      double iv = iterPL->dVal();
      saveIntervalPairs_.push_back( pair<double,double>(t,iv) );
    }

    if( currTag == "PACK" )
    {
      pack_ = iterPL->iVal();

#ifdef Xyce_PARALLEL_MPI
      if (!pack_ && !(pdsMgrPtr_->getPDSComm()->isSerial()))
      {
        N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_FATAL_0,
                                "Restart Data Must Be Packed for Parallel Runs, Please remove 'pack=0' from .OPTIONS RESTART line" );
      }
#endif
    }

    if( currTag == "INITIAL_INTERVAL" )
    {
      initialSaveInterval_ = iterPL->dVal();
    }

    if( currTag == "FILE" )
    {
      restartFileName_ = iterPL->sVal();
      restartFlag_ = true;
    }

    if( currTag == "JOB" )
    {
      restartJobName_ = iterPL->sVal();
    }

    if( currTag == "START_TIME" )
    {
      double t = iterPL->dVal();
      ost << t;
      restartFlag_ = true;
    }

    // move on to next tag
    ++iterPL;

    // FIXME ill-formed .OPTIONS RESTART lines may get through
  }

  if( ( restartFlag_ ) && ( restartFileName_ == "" ) )
  {
    restartFileName_ = restartJobName_ + ost.str();
  }

#ifdef Xyce_RESTART_NOPACK
  if( pack_ )
  {
    pack_ = false;
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_WARNING, "Turning off restart data packing!!! since this is broken in Windows!!!\n" );
  }
#endif

#ifdef Xyce_PARALLEL_MPI
  if( !pack_  && !(pdsMgrPtr_->getPDSComm()->isSerial()) )
  {
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL_0, ".OPTIONS Restart: Parallel Restarts must be packed!\n" );
  }
#endif

#ifdef Xyce_DEBUG_RESTART
  N_PDS_Comm * comm = pdsMgrPtr_->getPDSComm();
  if( comm->procID() == 0 )
  {
    string netListFile("");
    if (commandLine_.getArgumentValue("netlist") != "")
    {
      netListFile = commandLine_.getArgumentValue("netlist");
    }

    cout << "RESTART OPTIONS SETUP: ";
    cout << netListFile << endl;
    cout << "-----------------------\n";
    cout << "restartFileName: " << restartFileName_ << endl;
    cout << "restartJobName: " << restartJobName_ << endl;
    cout << "isRestart: " << restartFlag_ << endl;
    cout << "initial Interval: " << initialSaveInterval_ << endl;
    cout << "pack data: " << pack_ << endl;
    for( int i = 0; i < saveIntervalPairs_.size(); ++i )
      cout << saveIntervalPairs_[i].first << " " << saveIntervalPairs_[i].second << endl;
    cout << "-----------------------\n";
  }
#endif

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_RestartMgr::getRestartIntervals
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 7/19/01
//-----------------------------------------------------------------------------
bool N_IO_RestartMgr::getRestartIntervals(double & initialInterval,
                                 vector< pair<double,double> > & intervalPairs)
{
  initialInterval = initialSaveInterval_;
  intervalPairs = saveIntervalPairs_;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_RestartMgr::dumpRestartData
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 7/19/01
//-----------------------------------------------------------------------------
bool N_IO_RestartMgr::dumpRestartData(const double & time)
{
  if( pdsMgrPtr_ == 0 || topMgrPtr_ == 0 || devIntPtr_ == 0 || anaIntPtr_ == 0 )
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL,
      " ERROR: Restart Manager cannot access a package manager\n" );

  N_PDS_Comm * comm = pdsMgrPtr_->getPDSComm();

  bool success = (comm!=0);

  int procID = comm->procID();
  int numProcs = comm->numProc();

  //Get Restart Node Data and setup buffers
  //---------------------------------------
  int dataSize = sizeof(int);
  vector<N_IO_RestartNode*> nodeVec;

  topMgrPtr_->getRestartNodes( nodeVec );

  int nodeCount = nodeVec.size();
  double nC = static_cast<double>(nodeCount);

  double gNC = 0;
  comm->sumAll( &nC, &gNC, 1 );

  int globalNodeCount = static_cast<int>(gNC);

  for( int i = 0; i < nodeCount; ++i )
    dataSize += nodeVec[i]->packedByteCount();

  int maxSize = 0;
  int bsize;
  char * buf;
#ifdef Xyce_PARALLEL_MPI
  double dS = dataSize;
  double mS;
  comm->maxAll( &dS, &mS, 1 );
  maxSize = static_cast<int>(mS);
  if( procID == 0 )
  {
    buf = new char[maxSize];
    bsize = maxSize;
  }
  else
#endif
  {
    buf = new char[dataSize];
    bsize = dataSize;
    maxSize = dataSize;
  }

  //Setup Proc 0 and store its nodes
  //--------------------------------
  // This may be a bit inefficient, having a ofstream on the stack
  // for all processes rather than making one dynamically for proc=0
  // however, gcc 4.0 optimized the destruction improperly when it's
  // done dynamiclly.
  ofstream * outStream = NULL;
  ofstream outStreamSt;
  string outName;

  int proc = 0;
  if( procID == 0 )
  {
    ostringstream ost;
    ost << restartJobName_ << time;
    outName.assign(ost.str());
    outStreamSt.open( outName.c_str() );
    outStream = &outStreamSt;

    if( !outStream->is_open() )
    {
      N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_WARNING_0, "Cannot Open CheckPoint File: " + outName + "\n" );
      return false;
    }

#ifdef Xyce_DEBUG_RESTART
    cout << "DUMPING RESTART: " << outName << endl;
    cout << "----------------\n";
    cout << "numProcs: " << numProcs << " maxSize: " << maxSize << endl;
    cout << "proc: " << proc << " dataSize: " << dataSize << endl;
    cout << "nodeCount: " << nodeCount << endl;
    cout << "pack: " << pack_ << endl;
#endif

    int packed = pack_?1:0;
    (*outStream) << numProcs << " " << maxSize << " " << packed << " " << globalNodeCount << " ";
  }


  if( pack_ )
  {
    //Pack Restart Nodes
    int pos = 0;
    comm->pack( &nodeCount, 1, buf, bsize, pos );
    for( int i = 0; i < nodeCount; ++i )
    {
      nodeVec[i]->pack( buf, bsize, pos, comm );
    }

    if( procID == 0 )
    {
      (*outStream) << procID << " " << dataSize << " ";
      outStream->write( buf, dataSize );
    }
  }
  else
  {
    (*outStream) << " " << nodeCount << " ";
    for( int i = 0; i < nodeCount; ++i )
    {
      nodeVec[i]->dump(*outStream);
    }
  }

#ifdef Xyce_PARALLEL_MPI
  //PARALLEL, get and store nodes from other procs
  //----------------------------------------------
  comm->barrier();

  int size;
  for( proc = 1; proc < numProcs; ++proc )
  {
    if( procID == 0 )
    {
      comm->recv( &size, 1, proc );
      comm->recv( buf, size, proc );

      (*outStream) << proc << " " << size << " ";
      outStream->write( buf, size );
    }
    else if( procID == proc )
    {
      comm->send( &dataSize, 1, 0 );
      comm->send( buf, dataSize, 0 );
    }

    comm->barrier();
  }
#endif

  //Add Time Int stuff from proc 0
  //------------------------------
  if( procID == 0 )
  {
    //dataSize = restartDataSize();
    dataSize = anaIntPtr_->restartDataSize( pack_ );
    if( dataSize > maxSize )
    {
      delete [] buf;
      buf = 0;
      maxSize = dataSize;
      buf = new char[maxSize];
      bsize = maxSize;
    }

#ifdef Xyce_DEBUG_RESTART
    cout << "before tia RESTART dump: maxSize = " << maxSize << endl;
#endif
    int pos = 0;
    success &= anaIntPtr_->dumpRestartData( buf, bsize, pos, comm, pack_ );

    (*outStream) << dataSize << " ";
    outStream->write( buf, dataSize );
  }

#ifdef Xyce_PARALLEL_MPI
  comm->barrier();
#endif

#ifdef Xyce_DEBUG_RESTART
  if (procID == 0) cout << "---------------------------\n";
#endif

  //Add Device Mgr stuff from proc 0
  //------------------------------
  if( procID == 0 )
  {
    //dataSize = restartDataSize();
    dataSize = devIntPtr_->restartDataSize( pack_ );
    if( dataSize > maxSize )
    {
      delete [] buf;
      buf = 0;
      maxSize = dataSize;
      buf = new char[maxSize];
      bsize = maxSize;
    }

#ifdef Xyce_DEBUG_RESTART
    cout << "before dev RESTART dump: maxSize = " << maxSize << endl;
#endif
    int pos = 0;
    success &= devIntPtr_->dumpRestartData( buf, bsize, pos, comm, pack_ );

    (*outStream) << dataSize << " ";
    outStream->write( buf, dataSize );
  }

#ifdef Xyce_PARALLEL_MPI
  comm->barrier();
#endif

#ifdef Xyce_DEBUG_RESTART
  if (procID == 0) cout << "---------------------------\n";
#endif

  if( (procID == 0) && outStreamSt.is_open() )
  {
    outStreamSt.close();
  }

  int nodeSize = nodeVec.size();
  for (int i = 0; i < nodeSize; ++i)
    delete nodeVec[i];

  if( bsize != 0 ) delete [] buf;
  bsize = 0;

  return success;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_RestartMgr::restoreRestartData
// Purpose       :
// Special Notes : version used with distributed parser
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 7/19/01
//-----------------------------------------------------------------------------
bool N_IO_RestartMgr::restoreRestartData()
{
  if( pdsMgrPtr_ == 0 || topMgrPtr_ == 0 || devIntPtr_ == 0 || anaIntPtr_ == 0 )
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL,
      " ERROR: Restart Manager cannot access a package manager\n" );

  N_PDS_Comm * comm = pdsMgrPtr_->getPDSComm();

  bool success = ( comm != 0 );

  int procID = comm->procID();
  int numProcs = comm->numProc();

  char * buf = 0;
  int dataSize;
  int pos = 0;
  int bsize = 0;
  int fProcID;

  ifstream * inStream;

  int nodeCount;
  int oldNumProcs, maxSize;

  //Table of restart nodes used by all procs
  vector<N_IO_RestartNode*> nodeTable;

  //1. Proc 0 reads in data pushing NumDevices/NumProcs to every processor
  if (procID == 0)
  {
    inStream = new ifstream( restartFileName_.c_str());

    if (!inStream->is_open())
      N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_FATAL_0,
                              "Cannot Open CheckPoint File: " +
                              restartFileName_ + " for Restart\n" );

    int packed;
    int TotNumDevs;
    (*inStream) >> oldNumProcs >> maxSize >> packed >> TotNumDevs;
    pack_ = packed;

    char * buf1 = 0;
    int    bsize1 = 0;
    int    pos1   = 0;

#ifdef Xyce_DEBUG_RESTART
    cout << "RESTORING RESTART DATA\n";
    cout << "------------------------\n";
    cout << "oldNumProcs: " << oldNumProcs << endl;
    cout << "maxSize: " << maxSize << endl;
    cout << "packed: " << pack_ << endl;
    cout << "TotNumDevs: " << TotNumDevs << endl;
#endif

    int DevsPerProc = TotNumDevs/numProcs + 1;

    N_IO_RestartNode * nodeP;
    int CurrProc = 1;
    for( int oldProc = 0; oldProc < oldNumProcs; ++oldProc )
    {
      if( pack_ )
      {
        (*inStream) >> fProcID >> dataSize;
        char dummy = 'x';
        inStream->read( &dummy, 1 );

#ifdef Xyce_DEBUG_RESTART
        cout << "fProcID: " << fProcID << endl;
        cout << "dataSize: " << dataSize << endl;
#endif

        bsize = dataSize;
        buf = new char[bsize];

        inStream->read( buf, dataSize );

        pos = 0;
        comm->unpack( buf, bsize, pos, &nodeCount, 1 );
      }
      else
        (*inStream) >> nodeCount;

#ifdef Xyce_DEBUG_RESTART
      cout << "nodeCount: " << nodeCount << endl;
#endif

      for( int i = 0; i < nodeCount; ++i )
      {
        nodeP = new N_IO_RestartNode();

        if( pack_ )
          nodeP->unpack( buf, bsize, pos, comm );
        else
          nodeP->restore(*inStream);

        nodeTable.push_back( nodeP );

        if( ((nodeTable.size()==DevsPerProc) || (i==nodeCount-1)&&(oldProc==oldNumProcs-1)) )
        {
          //Pack up data and send to CurrProc
          int Size = 0;
          int nTSize = nodeTable.size();
          Size += sizeof(int);
          for( int j = 0; j < nTSize; ++j )
            Size += nodeTable[j]->packedByteCount();

          if( Size > bsize1 )
          {
            if( bsize1 ) delete [] buf1;
            bsize1 = Size;
            buf1 = new char[bsize1];
          }

          pos1 = 0;
          comm->pack( &nTSize, 1, buf1, bsize1, pos1 );
          for( int j = 0; j < nTSize; ++j )
            nodeTable[j]->pack( buf1, bsize1, pos1, comm );

          if( CurrProc < numProcs )
          {
            comm->send( &Size, 1, CurrProc );
            comm->send( buf1, Size, CurrProc );
            ++CurrProc;
          }

          for_each( nodeTable.begin(), nodeTable.end(), DeletePtr<N_IO_RestartNode>() );
          nodeTable.clear();
        }
      }
    }

    if( bsize ) delete [] buf;
    buf = buf1;
    bsize = bsize1;
  }

#ifdef Xyce_PARALLEL_MPI
  //Everybody else receive their set of devices
  if( procID != 0 )
  {
    comm->recv( &bsize, 1, 0 );
    buf = new char[bsize];
    comm->recv( buf, bsize, 0 );
  }
#ifdef Xyce_DEBUG_RESTART
  cout << "Finished ReadIn of Restart Node Data on Proc: " << procID << endl;
#endif

  comm->barrier();
#endif

#ifdef Xyce_PARALLEL_MPI
  //Set Buffers to largest necessary size
  double sz = static_cast<double>(bsize);
  double msz;
  comm->maxAll( &sz, &msz, 1 );
  int bsize2 = static_cast<int>(msz);
  char * buf2 = new char[bsize2];

  memcpy( buf2, buf, bsize );

  delete [] buf;
  bsize = bsize2;
  buf = new char[bsize];

  swap( buf, buf2 );
#endif

  //Check devices in your buffer, then push around ring
  for( int i = 0; i < numProcs; ++i )
  {
    pos = 0;
    comm->unpack( buf, bsize, pos, &nodeCount, 1 );
    nodeTable.resize(nodeCount);

#ifdef Xyce_DEBUG_RESTART
    if( !procID ) cout << "Restart Ring Stage: " << i << " " << nodeCount << endl;
#endif

    for( int j = 0; j < nodeCount; ++j )
    {
      nodeTable[j] = new N_IO_RestartNode();
      nodeTable[j]->unpack( buf, bsize, pos, comm );
    }

    topMgrPtr_->restoreRestartNodes( nodeTable );

    for_each( nodeTable.begin(), nodeTable.end(), DeletePtr<N_IO_RestartNode>() );
    nodeTable.clear();

#ifdef Xyce_PARALLEL_MPI
    int sendProc = procID-1;
    if( sendProc == -1 ) sendProc = numProcs-1;
    int recvProc = procID+1;
    if( recvProc == numProcs ) recvProc = 0;

    int MaxSendSize = 1e6;

    int NumSends = 1;
    int SendSize = bsize;
    if( bsize > MaxSendSize )
    {
      NumSends = bsize/MaxSendSize;
      if( bsize != NumSends*MaxSendSize ) ++NumSends;
      SendSize = MaxSendSize;
    }

#ifdef Xyce_DEBUG_RESTART
    cout << "NumSends: " << NumSends << endl;
#endif

    for( int j = 0; j < NumSends; ++j )
    {
      int Offset = j*SendSize;

      if( j == (NumSends-1) ) SendSize = bsize - Offset;

      comm->iRecv( buf2+Offset, SendSize, sendProc );

      comm->rSend( buf+Offset, SendSize, recvProc );

      comm->waitAll();
    }

    swap( buf, buf2 );
#endif
  }

#ifdef Xyce_PARALLEL_MPI
  bsize2 = 0;
  delete [] buf2;
  buf2 = 0;
#endif

// read in Time data
  if( procID == 0 )
  {
    int proc;
    (*inStream) >> dataSize;
    char dummy = 'x';
    inStream->read( &dummy, 1 );

#ifdef Xyce_DEBUG_RESTART
    cout << "RESTORING TIME INT STUFF\n";
    cout << "dataSize: " << dataSize << endl;
    cout << "-----------------------\n";
#endif

    if( dataSize > maxSize )
    {
      delete [] buf;
      maxSize = dataSize;
      buf = new char[maxSize];
      bsize = maxSize;
    }

    inStream->read( buf, dataSize );
  }

#ifdef Xyce_PARALLEL_MPI
  comm->barrier();

  comm->bcast( &dataSize, 1, 0 );
  if( procID != 0 )
    if( dataSize > bsize )
    {
      delete [] buf;
      buf = new char[dataSize];
      bsize = dataSize;
    }

  comm->bcast( buf, dataSize, 0 );
#endif

  pos = 0;
  anaIntPtr_->restoreRestartData( buf, bsize, pos, comm, pack_ );

// read in Dev data
  if( procID == 0 )
  {
    int proc;
    (*inStream) >> dataSize;
    char dummy = 'x';
    inStream->read( &dummy, 1 );

#ifdef Xyce_DEBUG_RESTART
    cout << "RESTORING DEV MGR STUFF\n";
    cout << "dataSize: " << dataSize << endl;
    cout << "-----------------------\n";
#endif

    if( dataSize > maxSize )
    {
      delete [] buf;
      maxSize = dataSize;
      buf = new char[maxSize];
      bsize = maxSize;
    }

    inStream->read( buf, dataSize );

    delete inStream;
  }

#ifdef Xyce_PARALLEL_MPI
  comm->barrier();

  comm->bcast( &dataSize, 1, 0 );
  if( procID != 0 )
    if( dataSize > bsize )
    {
      delete [] buf;
      buf = new char[dataSize];
      bsize = dataSize;
    }

  comm->bcast( buf, dataSize, 0 );
#endif

  pos = 0;
  devIntPtr_->restoreRestartData( buf, bsize, pos, comm, pack_ );
  delete [] buf;


  return success;
}

//-----------------------------------------------------------------------------
// Function      : N_IO_RestartMgr::restartDataSize
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 7/19/01
//-----------------------------------------------------------------------------
int N_IO_RestartMgr::restartDataSize()
{
  if( pdsMgrPtr_ == 0 || topMgrPtr_ == 0 || devIntPtr_ == 0 || anaIntPtr_ == 0 )
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL,
      " ERROR: Restart Manager cannot access a package manager\n" );

  int devSize = devIntPtr_->restartDataSize( pack_ );
  int anaSize = anaIntPtr_->restartDataSize( pack_ );
  int size = anaSize+devSize;

  return size;
}

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
// Filename       : $RCSfile: N_PDS_GlobalAccessor.C,v $
//
// Purpose        : Simple migrator utility using Zoltan utilities
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 06/06/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.17.6.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:49 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// --------- Standard Includes ----------

// --------- Xyce Includes --------------

#include <N_ERH_ErrorMgr.h>

#include <N_LAS_MultiVector.h>

#include <N_PDS_GlobalAccessor.h>
#include <N_PDS_ParComm.h>
#include <N_PDS_MPIComm.h>
#include <N_PDS_ParMap.h>

#ifdef Xyce_PARALLEL_MPI
//  #include <GSComm_Plan.h>
//  #include <GSComm_Comm.h>
 #include <Epetra_MpiComm.h>
 #include <Epetra_MpiDistributor.h>
 #include <Epetra_SerialComm.h>
 #include <Epetra_SerialDistributor.h>
#else
 #include <Epetra_SerialDistributor.h>
#endif

// --------- Other Includes -------------

//-----------------------------------------------------------------------------
// Function      : N_PDS_GlobalAccessor::N_PDS_GlobalAccessor
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Compuational Sciences
// Creation Date : 02/06/01
//-----------------------------------------------------------------------------
N_PDS_GlobalAccessor::N_PDS_GlobalAccessor( N_PDS_Comm * comm )
 : pdsComm_(comm),
   arrayReceiveGIDs_(0),
   arrayReceiveLIDs_(0),
   arrayReceiveProcs_(0),
   arraySendGIDs_(0),
   arraySendLIDs_(0),
   arraySendProcs_(0),
   numReceiveObjs_(0),
   numSendObjs_(0),
   recvBuf_(0),
   recvBufSize_(0),
   sendBuf_(0),
   sendBufSize_(0)
                   
#ifdef Xyce_PARALLEL_MPI
//                  ,
//                  GSComm_(0),
//                  GSPlan_(0)
#endif
                  ,
                  distributor_(0)
{
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_GlobalAccessor::~N_PDS_GlobalAccessor
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Compuational Sciences
// Creation Date : 06/07/00
//-----------------------------------------------------------------------------
N_PDS_GlobalAccessor::~N_PDS_GlobalAccessor()
{

  if( arrayReceiveGIDs_ ) delete [] arrayReceiveGIDs_;
  if( arrayReceiveLIDs_ ) delete [] arrayReceiveLIDs_;
  if( arrayReceiveProcs_ ) delete [] arrayReceiveProcs_;
  if( recvBuf_ ) delete [] recvBuf_;

  if( arraySendGIDs_ ) delete [] arraySendGIDs_;
  if( arraySendLIDs_ ) delete [] arraySendLIDs_;
  if( arraySendProcs_ ) delete [] arraySendProcs_;
  if( sendBuf_ ) delete [] sendBuf_;

  if( distributor_ ) delete distributor_;

}

//-----------------------------------------------------------------------------
// Function      : N_PDS_GlobalAccessor::generateMigrationPlan
// Purpose       : generates migration plan (Comm_Obj) using Zoltan
// 		   utilities
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/07/00
//-----------------------------------------------------------------------------
void N_PDS_GlobalAccessor::generateMigrationPlan()
{
  numReceiveObjs_ = externGIDVector_.size();

  map<int,int> sortMap;

  // sort externGIDVector_ numerically by processor
  for( int i = 0; i < numReceiveObjs_; ++i )
    sortMap[ externGIDVector_[i].second ]++;

  int count = 0;
  int tmp;
  for( map<int,int>::iterator it_iiM = sortMap.begin();
       it_iiM != sortMap.end(); ++it_iiM )
  {
    tmp = count;
    count += it_iiM->second;
    it_iiM->second = tmp;
  }

  if( !arrayReceiveGIDs_ )
  {
    arrayReceiveGIDs_ = new int[ numReceiveObjs_ ];
    arrayReceiveLIDs_ = new int[ numReceiveObjs_ ];
    arrayReceiveProcs_ = new int[ numReceiveObjs_ ];
  }
  else
  {
    delete [] arraySendGIDs_;
    delete [] arraySendLIDs_;
    delete [] arraySendProcs_;

    delete [] sendBuf_;
    delete [] recvBuf_;
  }

  int loc, val1, val2;
  for( int i = 0; i < numReceiveObjs_; ++i )
  {
    val1 = externGIDVector_[i].first;
    val2 = externGIDVector_[i].second;
    loc = sortMap[ val2 ];
    arrayReceiveGIDs_[ loc ] = val1;
    arrayReceiveLIDs_[ loc ] = val1;
    arrayReceiveProcs_[ loc ] = val2;
    sortMap[ val2 ]++;
  }

#ifdef Xyce_DEBUG_PARALLEL
  cout << "N_PDS_GlobalAccessor::generateMigrationPlan:" << endl;
  cout << " setup numRecvObjs: " << numReceiveObjs_ << endl;
  cout << " setup numSendObjs: " << numSendObjs_ << endl;

  for( int i = 0; i < numReceiveObjs_; ++i )
    cout << "  " << arrayReceiveGIDs_[i] << " " << arrayReceiveProcs_[i]
	<< endl;
#endif

#ifdef Xyce_PARALLEL_MPI

  if( !distributor_ ) 
  {
    if( pdsComm_->isSerial() )
    {
      distributor_ = new Epetra_SerialDistributor(
                    *(dynamic_cast<Epetra_SerialComm*>(pdsComm_->petraComm())));
    }
    else
    {
      distributor_ = new Epetra_MpiDistributor(
                    *(dynamic_cast<Epetra_MpiComm*>(pdsComm_->petraComm())) );
      // can only call CreateFromRecvs if we're parallel with more than one proc
      distributor_->CreateFromRecvs( numReceiveObjs_, arrayReceiveGIDs_, arrayReceiveProcs_,
	                         true, numSendObjs_, arraySendGIDs_, arraySendProcs_ );
    }
  }

  sendGIDVector_.resize( numSendObjs_ );

  for( int i = 0; i < numSendObjs_; ++i )
    sendGIDVector_[i] = pair<int,int>( arraySendGIDs_[i],
			arraySendProcs_[i] );

  sendBufSize_ = numSendObjs_ * (sizeof(int)+sizeof(double));
  recvBufSize_ = numReceiveObjs_ * (sizeof(int)+sizeof(double));
  sendBuf_ = new char[ sendBufSize_ ];
  recvBuf_ = new char[ recvBufSize_ ];

#ifdef Xyce_DEBUG_PARALLEL
  cout << "Created Migration Plan: " << endl;
  cout << " numRecvObjs: " << numReceiveObjs_ << endl;
  for( int i = 0; i < numReceiveObjs_; ++i )
    cout << "  " << arrayReceiveGIDs_[i] << " " << arrayReceiveProcs_[i]
	<< endl;
  cout << " numSendObjs: " << numSendObjs_ << endl;
  for( int i = 0; i < numSendObjs_; ++i )
    cout << "  " << arraySendGIDs_[i] << " " << arraySendProcs_[i]
	<< endl;
#endif

#endif /* Xyce_PARALLEL_MPI */

}

//-----------------------------------------------------------------------------
// Function      : N_PDS_GlobalAccessor::migrateMultiVector
// Purpose       : migrates nonlocal parts of multivector based on
//                 migration plan (theZoltanCommObjPtr_)
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/07/00
//-----------------------------------------------------------------------------
void N_PDS_GlobalAccessor::migrateMultiVector( N_LAS_MultiVector * mVector )
{

#ifdef Xyce_DEBUG_PARALLEL
  cout << "migrating vector<<<<<<<<<<<<<<" << endl;
  cout << " " << numSendObjs_ << " " << numReceiveObjs_ << endl;
#endif

#ifdef Xyce_PARALLEL_MPI

  int pos = 0;
  int idx;
  double val;
  for( int i = 0; i < numSendObjs_; ++i )
  {
    val = mVector->getElementByGlobalIndex( sendGIDVector_[i].first );

    idx = sendGIDVector_[i].first;
    pdsComm_->pack( &idx, 1, sendBuf_, sendBufSize_, pos );
    pdsComm_->pack( &val, 1, sendBuf_, sendBufSize_, pos );
  }

  distributor_->Do( sendBuf_, sizeof(int)+sizeof(double), recvBufSize_, recvBuf_ );

  mVector->clearExternVectorMap();

  pos = 0;
  for( int i = 0; i < numReceiveObjs_; ++i )
  {
    pdsComm_->unpack( recvBuf_, numReceiveObjs_*(sizeof(int)+sizeof(double)), pos, &idx, 1 );
    pdsComm_->unpack( recvBuf_, numReceiveObjs_*(sizeof(int)+sizeof(double)), pos, &val, 1 );

    mVector->addElementToExternVectorMap( idx, val );
  }

#endif /* Xyce_PARALLEL_MPI */

}

//-----------------------------------------------------------------------------
// Function      : N_PDS_GlobalAccessor::migrateIntArray
// Purpose       : migrates nonlocal parts of integer array based on
//                 migration plan (theZoltanCommObjPtr_)
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/07/00
//-----------------------------------------------------------------------------
void N_PDS_GlobalAccessor::migrateIntArray( map<int,int> & sendMap,
	map<int,int> & recvMap )
{

#ifdef Xyce_PARALLEL_MPI

  map<int,int>::iterator it_i2M, end_i2M;

#ifdef Xyce_DEBUG_PARALLEL
  cout << "Send Map" << endl;
  it_i2M = sendMap.begin();
  end_i2M = sendMap.end();
  for( ; it_i2M != end_i2M; ++it_i2M )
    cout << "  " << it_i2M->first << " " << it_i2M->second << endl;
  cout << endl;
#endif

  int pos = 0;
  int idx;
  int val;
  for( int i = 0; i < numSendObjs_; ++i )
  {
    idx = sendGIDVector_[i].first;
    val = sendMap[ sendGIDVector_[i].first ];

#ifdef Xyce_DEBUG_PARALLEL
    cout << "send values: " << sendGIDVector_[i].first <<
	" " << val << endl;
#endif

    pdsComm_->pack( &idx, 1, sendBuf_, 2 * numSendObjs_ * sizeof( int ), pos );
    pdsComm_->pack( &val, 1, sendBuf_, 2 * numSendObjs_ * sizeof( int ), pos );
  }

  distributor_->Do( sendBuf_, 2*sizeof(int), recvBufSize_, recvBuf_ );

  recvMap.clear();

  pos = 0;
  for( int i = 0; i < numReceiveObjs_; ++i )
  {
    pdsComm_->unpack( recvBuf_, 2 * numReceiveObjs_ * sizeof( int ), pos, &idx, 1 );
    pdsComm_->unpack( recvBuf_, 2 * numReceiveObjs_ * sizeof( int ), pos, &val, 1 );

    recvMap[ idx ] = val;
  }

#ifdef Xyce_DEBUG_PARALLEL
  cout << "Recv Map" << endl;
  it_i2M = recvMap.begin();
  end_i2M = recvMap.end();
  for( ; it_i2M != end_i2M; ++it_i2M )
    cout << "  " << it_i2M->first << " " << it_i2M->second << endl;
  cout << endl;
#endif

#endif /* Xyce_PARALLEL_MPI */

}

//-----------------------------------------------------------------------------
// Function      : N_PDS_GlobalAccessor::migrateIntVecs
// Purpose       : migrates nonlocal parts of integer vectors
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
void N_PDS_GlobalAccessor::migrateIntVecs( map< int,vector<int> > & sendMap,
	map< int,vector<int> > & recvMap )
{

#ifdef Xyce_PARALLEL_MPI

  map< int,vector<int> >::iterator it_i2M, end_i2M;

  int maxSize = 0;
  int maxGlobalSize;

  double tmpVar1, tmpVar2;

  it_i2M = sendMap.begin();
  end_i2M = sendMap.end();
  for( ; it_i2M != end_i2M; ++it_i2M )
    if( maxSize < it_i2M->second.size() ) maxSize = it_i2M->second.size();

  tmpVar1 = maxSize;
  pdsComm_->maxAll( &tmpVar1, &tmpVar2, 1 );
  maxGlobalSize = tmpVar2;

  int pos = 0;
  int idx, val;

  map<int,int> recvSizeMap;

  for( int i = 0; i < numSendObjs_; ++i )
  {
    idx = sendGIDVector_[i].first;
    val = sendMap[ idx ].size();

    pdsComm_->pack( &idx, 1, sendBuf_, 2 * numSendObjs_ * sizeof( int ), pos );
    pdsComm_->pack( &val, 1, sendBuf_, 2 * numSendObjs_ * sizeof( int ), pos );
  }

  distributor_->Do( sendBuf_, 2*sizeof(int), recvBufSize_, recvBuf_ );

  recvSizeMap.clear();

  pos = 0;
  for( int i = 0; i < numReceiveObjs_; ++i )
  {
    pdsComm_->unpack( recvBuf_, 2 * numReceiveObjs_ * sizeof( int ), pos, &idx, 1 );
    pdsComm_->unpack( recvBuf_, 2 * numReceiveObjs_ * sizeof( int ), pos, &val, 1 );

    recvSizeMap[ idx ] = val;
    recvMap[ idx ] = vector<int>(val,0);
  }

  for( int loc = 0; loc < maxGlobalSize; ++loc )
  {
    pos = 0;

    for( int i = 0; i < numSendObjs_; ++i )
    {
      idx = sendGIDVector_[i].first;
      val = -1;
      if( loc < sendMap[ idx ].size() )
        val = (sendMap[ idx ])[loc];

      pdsComm_->pack( &idx, 1, sendBuf_, 2 * numSendObjs_ * sizeof( int ), pos );
      pdsComm_->pack( &val, 1, sendBuf_, 2 * numSendObjs_ * sizeof( int ), pos );
    }

    distributor_->Do( sendBuf_, 2*sizeof(int), recvBufSize_, recvBuf_ );

    pos = 0;
    for( int i = 0; i < numReceiveObjs_; ++i )
    {
      pdsComm_->unpack( recvBuf_, 2 * numReceiveObjs_ * sizeof( int ), pos, &idx, 1 );
      pdsComm_->unpack( recvBuf_, 2 * numReceiveObjs_ * sizeof( int ), pos, &val, 1 );

      if( loc < recvMap[ idx ].size() ) (recvMap[ idx ])[ loc ] = val;
    }

  }

#endif /* Xyce_PARALLEL_MPI */

}


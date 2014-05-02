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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_PDS_ParComm.C,v $
//
// Purpose        : Implementation file for the parallel communication
//                  class for Xyce.
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 06/26/01
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.21 $
//
// Revision Date  : $Date: 2014/02/24 23:49:25 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

#ifdef Xyce_PARALLEL_MPI
 #include <mpi.h>
#endif

// ----------   Xyce Includes   ----------

#include <N_PDS_ParComm.h>

#ifdef Xyce_PARALLEL_MPI
 #include <N_PDS_MPIComm.h>
#endif

#include <N_ERH_ErrorMgr.h>

// ----------  Other Includes   ----------

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParComm::N_PDS_ParComm
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
N_PDS_ParComm::N_PDS_ParComm( int iargs,
                              char * cargs[] )
  :
  petraComm_(0),
  petraCommOwned_(true)
#ifdef Xyce_PARALLEL_MPI
  ,
  mpiCommOwned_(true)
#endif
{
#ifdef Xyce_PARALLEL_MPI
  mpiComm_ = new N_PDS_MPIComm( iargs, cargs );

  petraComm_ = new Epetra_MpiComm( mpiComm_->comm() );
  petraCommOwned_ = true;
#else
  petraComm_ = new Epetra_SerialComm();
  petraCommOwned_ = true;
#endif

  isSerial_ = ( numProc() == 1 );
}

#ifdef Xyce_PARALLEL_MPI
//-----------------------------------------------------------------------------
// Function      : N_PDS_ParComm::N_PDS_ParComm
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 03/16/06
//-----------------------------------------------------------------------------
N_PDS_ParComm::N_PDS_ParComm( MPI_Comm comm )
  :
  petraCommOwned_(true),
  mpiCommOwned_(true)
{
  mpiComm_ = new N_PDS_MPIComm( comm );
  petraComm_ = new Epetra_MpiComm( mpiComm_->comm() );

  isSerial_ = ( numProc() == 1 );
}
#endif

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParComm::N_PDS_ParComm
// Purpose       : Copy constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
N_PDS_ParComm::N_PDS_ParComm(const N_PDS_ParComm &right)
  :
  isSerial_(right.isSerial_) ,
  petraComm_(right.petraComm_),
  petraCommOwned_(false)
#ifdef Xyce_PARALLEL_MPI
  ,
  mpiComm_(right.mpiComm_),
  mpiCommOwned_(false)
#endif
{
}

Xyce::Parallel::Machine N_PDS_ParComm::comm() const {
#ifdef Xyce_PARALLEL_MPI
      return mpiComm_->comm();
#else
      return 0;
#endif
    }
    
//-----------------------------------------------------------------------------
// Function      : N_PDS_ParComm::operator=
// Purpose       : Assignment operator.
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
N_PDS_ParComm & N_PDS_ParComm::operator=(const N_PDS_ParComm &right)
{
  isSerial_ = right.isSerial_;
  petraComm_  = right.petraComm_;
  petraCommOwned_ = false;
#ifdef Xyce_PARALLEL_MPI
  mpiComm_  = right.mpiComm_;
  mpiCommOwned_ = false;
#endif

  return *this;
}
//-----------------------------------------------------------------------------
// Function      : N_PDS_ParComm::~N_PDS_ParComm
// Purpose       : Destructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoeksta, SNL, Parallel Computational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
N_PDS_ParComm::~N_PDS_ParComm()
{
  if( petraCommOwned_ )
    if( petraComm_ ) delete petraComm_;

#ifdef Xyce_PARALLEL_MPI
  if( mpiCommOwned_ )
    if( mpiComm_ ) delete mpiComm_;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParComm::numProc
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/27/01
//-----------------------------------------------------------------------------
int N_PDS_ParComm::numProc() const
{
#ifdef Xyce_PARALLEL_MPI
  return mpiComm_->numProc();
#else
  return 1;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParComm::procID
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/27/01
//-----------------------------------------------------------------------------
int N_PDS_ParComm::procID() const
{
#ifdef Xyce_PARALLEL_MPI
  return mpiComm_->procID();
#else
  return 0;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParComm::scanSum
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_ParComm::scanSum( const double * vals, double * sums,
		const int & count ) const
{
  return ( petraComm_->ScanSum( const_cast<double *> (vals), sums,
		const_cast<int &> (count) ) == 0 );
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_Comm::sumAll
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_ParComm::sumAll( const double * vals, double * sums,
		const int & count ) const
{
  return ( petraComm_->SumAll( const_cast<double *> (vals), sums,
		const_cast<int &> (count) ) == 0 );
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParComm::maxAll
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_ParComm::maxAll( const double * vals, double * maxs,
		const int & count ) const
{
  return ( petraComm_->MaxAll( const_cast<double *> (vals), maxs,
		const_cast<int &> (count) ) == 0 );
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParComm::minAll
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/30/01
//-----------------------------------------------------------------------------
bool N_PDS_ParComm::minAll( const double * vals, double * mins,
		const int & count ) const
{
  return ( petraComm_->MinAll( const_cast<double *> (vals), mins,
		const_cast<int &> (count) ) == 0 );
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParComm::scanSum
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Dave Shirley, PSSI
// Creation Date : 06/18/05
//-----------------------------------------------------------------------------
bool N_PDS_ParComm::scanSum( const int * vals, int * sums,
		const int & count ) const
{
  return ( petraComm_->ScanSum( const_cast<int *> (vals), sums,
		const_cast<int &> (count) ) == 0 );
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_Comm::sumAll
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Dave Shirley, PSSI
// Creation Date : 06/18/05
//-----------------------------------------------------------------------------
bool N_PDS_ParComm::sumAll( const int * vals, int * sums,
		const int & count ) const
{
  return ( petraComm_->SumAll( const_cast<int *> (vals), sums,
		const_cast<int &> (count) ) == 0 );
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParComm::maxAll
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Dave Shirley, PSSI
// Creation Date : 06/18/05
//-----------------------------------------------------------------------------
bool N_PDS_ParComm::maxAll( const int * vals, int * maxs,
		const int & count ) const
{
  return ( petraComm_->MaxAll( const_cast<int *> (vals), maxs,
		const_cast<int &> (count) ) == 0 );
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParComm::minAll
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Dave Shirley, PSSI
// Creation Date : 06/18/05
//-----------------------------------------------------------------------------
bool N_PDS_ParComm::minAll( const int * vals, int * mins,
		const int & count ) const
{
  return ( petraComm_->MinAll( const_cast<int *> (vals), mins,
		const_cast<int &> (count) ) == 0 );
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParComm::bcast
// Purpose       : INTs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_ParComm::bcast( int * val, const int & count, const int & root )
									const
{
#ifdef Xyce_PARALLEL_MPI
  return mpiComm_->bcast( val, count, root );
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParComm::bcast
// Purpose       : LNGs
// Special Notes :
// Scope         : Public
// Creator       : Derek Barnes, SNL, Parallel Computational Sciences
// Creation Date : 08/01/01
//-----------------------------------------------------------------------------
bool N_PDS_ParComm::bcast( long * val, const int & count, const int & root )
									const
{
#ifdef Xyce_PARALLEL_MPI
  return mpiComm_->bcast( val, count, root );
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParComm::bcast
// Purpose       : CHARs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_ParComm::bcast( char * val, const int & count, const int & root )
									const
{
#ifdef Xyce_PARALLEL_MPI
  return mpiComm_->bcast( val, count, root );
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParComm::bcast
// Purpose       : DBLEs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_ParComm::bcast( double * val, const int & count, const int & root )
									const
{
#ifdef Xyce_PARALLEL_MPI
  return mpiComm_->bcast( val, count, root );
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParComm::send
// Purpose       : INTs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_ParComm::send( const int * val, const int & count, const int & dest)
									const
{
#ifdef Xyce_PARALLEL_MPI
  return mpiComm_->send( val, count, dest );
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParComm::send
// Purpose       : LNGSs
// Special Notes :
// Scope         : Public
// Creator       : Derek Barnes, SNL, Parallel Computational Sciences
// Creation Date : 08/01/01
//-----------------------------------------------------------------------------
bool N_PDS_ParComm::send( const long * val, const int & count, const int & dest)
									const
{
#ifdef Xyce_PARALLEL_MPI
  return mpiComm_->send( val, count, dest );
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParComm::send
// Purpose       : CHARs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_ParComm::send( const char * val, const int & count, const int & dest)
									const
{
#ifdef Xyce_PARALLEL_MPI
  return mpiComm_->send( val, count, dest );
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParComm::send
// Purpose       : DBLEs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_ParComm::send( const double * val, const int & count, const int & dest)									const
{
#ifdef Xyce_PARALLEL_MPI
  return mpiComm_->send( val, count, dest );
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParComm::recv
// Purpose       : INTs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_ParComm::recv( int * val, const int & count, const int & src ) const
{
#ifdef Xyce_PARALLEL_MPI
  return mpiComm_->recv( val, count, src );
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParComm::recv
// Purpose       : LNGSs
// Special Notes :
// Scope         : Public
// Creator       : Derek Barnes, SNL, Parallel Computational Sciences
// Creation Date : 08/01/01
//-----------------------------------------------------------------------------
bool N_PDS_ParComm::recv( long * val, const int & count, const int & src ) const
{
#ifdef Xyce_PARALLEL_MPI
  return mpiComm_->recv( val, count, src );
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParComm::recv
// Purpose       : CHARs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_ParComm::recv( char * val, const int & count, const int & src ) const
{
#ifdef Xyce_PARALLEL_MPI
  return mpiComm_->recv( val, count, src );
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParComm::recv
// Purpose       : DBLEs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_ParComm::recv( double * val, const int & count, const int & src ) const
{
#ifdef Xyce_PARALLEL_MPI
  return mpiComm_->recv( val, count, src );
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParComm::rSend
// Purpose       : INTs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_ParComm::rSend( const int * val, const int & count, const int & dest) const
{
#ifdef Xyce_PARALLEL_MPI
  return mpiComm_->rSend( val, count, dest );
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParComm::rSend
// Purpose       : LNGSs
// Special Notes :
// Scope         : Public
// Creator       : Derek Barnes, SNL, Parallel Computational Sciences
// Creation Date : 08/01/01
//-----------------------------------------------------------------------------
bool N_PDS_ParComm::rSend( const long * val, const int & count, const int & dest) const
{
#ifdef Xyce_PARALLEL_MPI
  return mpiComm_->rSend( val, count, dest );
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParComm::rSend
// Purpose       : CHARs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_ParComm::rSend( const char * val, const int & count, const int & dest) const
{
#ifdef Xyce_PARALLEL_MPI
  return mpiComm_->rSend( val, count, dest );
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParComm::rSend
// Purpose       : DBLEs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_ParComm::rSend( const double * val, const int & count, const int & dest) const
{
#ifdef Xyce_PARALLEL_MPI
  return mpiComm_->rSend( val, count, dest );
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParComm::iRecv
// Purpose       : INTs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_ParComm::iRecv( int * val, const int & count, const int & src ) const
{
#ifdef Xyce_PARALLEL_MPI
  return mpiComm_->iRecv( val, count, src );
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParComm::iRecv
// Purpose       : LNGSs
// Special Notes :
// Scope         : Public
// Creator       : Derek Barnes, SNL, Parallel Computational Sciences
// Creation Date : 08/01/01
//-----------------------------------------------------------------------------
bool N_PDS_ParComm::iRecv( long * val, const int & count, const int & src ) const
{
#ifdef Xyce_PARALLEL_MPI
  return mpiComm_->iRecv( val, count, src );
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParComm::iRecv
// Purpose       : CHARs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_ParComm::iRecv( char * val, const int & count, const int & src ) const
{
#ifdef Xyce_PARALLEL_MPI
  return mpiComm_->iRecv( val, count, src );
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParComm::iRecv
// Purpose       : DBLEs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_ParComm::iRecv( double * val, const int & count, const int & src ) const
{
#ifdef Xyce_PARALLEL_MPI
  return mpiComm_->iRecv( val, count, src );
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParComm::waitAll
// Purpose       : 
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_ParComm::waitAll()
{
#ifdef Xyce_PARALLEL_MPI
  return mpiComm_->waitAll();
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParComm::pack
// Purpose       : INTs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_ParComm::pack( const int * val, const int count, char * buf,
		const int size, int & pos ) const
{
#ifdef Xyce_PARALLEL_MPI
  return mpiComm_->pack( val, count, buf, size, pos );
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParComm::pack
// Purpose       : LNGs
// Special Notes :
// Scope         : Public
// Creator       : Derek Barnes, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_ParComm::pack( const long * val, const int count, char * buf,
		const int size, int & pos ) const
{
#ifdef Xyce_PARALLEL_MPI
  return mpiComm_->pack( val, count, buf, size, pos );
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParComm::pack
// Purpose       : CHARs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_ParComm::pack( const char * val, const int count, char * buf,
		const int size, int & pos ) const
{
#ifdef Xyce_PARALLEL_MPI
  return mpiComm_->pack( val, count, buf, size, pos );
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParComm::pack
// Purpose       : DBLEs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_ParComm::pack( const double * val, const int count, char * buf,
		const int size, int & pos ) const
{
#ifdef Xyce_PARALLEL_MPI
  return mpiComm_->pack( val, count, buf, size, pos );
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParComm::pack
// Purpose       : size_t's
// Special Notes :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 06/26/2013
//-----------------------------------------------------------------------------
bool N_PDS_ParComm::pack( const std::size_t * val, const int count, char * buf,
		const int size, int & pos ) const
{
#ifdef Xyce_PARALLEL_MPI
  return mpiComm_->pack( val, count, buf, size, pos );
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParComm::unpack
// Purpose       : INTs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_ParComm::unpack( const char * buf, const int size, int & pos,
		int * val, const int count ) const
{
#ifdef Xyce_PARALLEL_MPI
  return mpiComm_->unpack( buf, size, pos, val, count );
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParComm::unpack
// Purpose       : LNGs
// Special Notes :
// Scope         : Public
// Creator       : Derek Barnes, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_ParComm::unpack( const char * buf, const int size, int & pos,
		long * val, const int count ) const
{
#ifdef Xyce_PARALLEL_MPI
  return mpiComm_->unpack( buf, size, pos, val, count );
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParComm::unpack
// Purpose       : CHARs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_ParComm::unpack( const char * buf, const int size, int & pos,
		char * val, const int count ) const
{
#ifdef Xyce_PARALLEL_MPI
  return mpiComm_->unpack( buf, size, pos, val, count );
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParComm::unpack
// Purpose       : DBLEs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_ParComm::unpack( const char * buf, const int size, int & pos,
		double * val, const int count ) const
{
#ifdef Xyce_PARALLEL_MPI
  return mpiComm_->unpack( buf, size, pos, val, count );
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParComm::unpack
// Purpose       : size_t's
// Special Notes :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 06/26/2013
//-----------------------------------------------------------------------------
bool N_PDS_ParComm::unpack( const char * buf, const int size, int & pos,
		std::size_t * val, const int count ) const
{
#ifdef Xyce_PARALLEL_MPI
  return mpiComm_->unpack( buf, size, pos, val, count );
#else
  return false;
#endif
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParComm::barrier
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
void N_PDS_ParComm::barrier() const
{
#ifdef Xyce_PARALLEL_MPI
  mpiComm_->barrier();
#endif
}


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
// Filename       : $RCSfile: N_PDS_MPIComm.h,v $
//
// Purpose        : Specification file for the abstract parallel communication
//                  class for Xyce.  This class will contain parallel data and
//                  functions.
//
// Special Notes  : It assumes that all parallel communication is being handled
//                  through MPI.
//
// Creator        : Scott A. Hutchinson, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/08/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.20 $
//
// Revision Date  : $Date: 2014/02/24 23:49:25 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef Xyce_N_PDS_MPIComm_h
#define Xyce_N_PDS_MPIComm_h

#ifdef Xyce_PARALLEL_MPI

#include <N_ERH_ErrorMgr.h>

#include <mpi.h>
#include<list>

//-----------------------------------------------------------------------------
// Class         : N_PDS_MPIComm
// Purpose       : MPI specific comm object, hides MPI functionality
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Compuational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
class N_PDS_MPIComm
{

public:

  N_PDS_MPIComm( int iargs, char * cargs[] );
  N_PDS_MPIComm( const MPI_Comm mComm );

  //Destructor
  ~N_PDS_MPIComm();

private:
  //No copy or assignment
  N_PDS_MPIComm( const N_PDS_MPIComm & right );
  N_PDS_MPIComm & operator=(const N_PDS_MPIComm & right);

public:
  // Get my processor ID.
  int procID() const;

  // Get the total number of processors in this communicator.
  int numProc() const;

  MPI_Comm comm() const;

  // MPI_Bcast wrappers
  bool bcast( int * val, const int & count, const int & root ) const;
  bool bcast( char * val, const int & count, const int & root ) const;
  bool bcast( double * val, const int & count, const int & root ) const;
  bool bcast( long * val, const int & count, const int & root ) const;

  // MPI_Send wrappers
  bool send( const int * val, const int & count, const int & dest ) const;
  bool send( const char * val, const int & count, const int & dest ) const;
  bool send( const double * val, const int & count, const int & dest ) const;
  bool send( const long * val, const int & count, const int & dest ) const;

  // MPI_Irecv wrappers
  bool recv( int * val, const int & count, const int & src ) const;
  bool recv( char * val, const int & count, const int & src ) const;
  bool recv( double * val, const int & count, const int & src ) const;
  bool recv( long * val, const int & count, const int & src ) const;

  // MPI_Rsend wrappers
  bool rSend( const int * val, const int & count, const int & dest ) const;
  bool rSend( const char * val, const int & count, const int & dest ) const;
  bool rSend( const double * val, const int & count, const int & dest ) const;
  bool rSend( const long * val, const int & count, const int & dest ) const;

  // MPI_Irecv wrappers
  bool iRecv( int * val, const int & count, const int & src );
  bool iRecv( char * val, const int & count, const int & src );
  bool iRecv( double * val, const int & count, const int & src );
  bool iRecv( long * val, const int & count, const int & src );

  // MPI_Waitall wrappers
  bool waitAll();

  // MPI_Pack wrappers
  bool pack( const int * val, const int count, char * buf, const int size, int & pos ) const;
  bool pack( const char * val, const int count, char * buf, const int size, int & pos ) const;
  bool pack( const double * val, const int count, char * buf, const int size, int & pos ) const;
  bool pack( const long * val, const int count, char * buf, const int size, int & pos ) const;
  bool pack( const std::size_t * val, const int count, char * buf, const int size, int & pos ) const;

  // MPI_Unpack wrappers
  bool unpack( const char * buf, const int size, int & pos, int * val, const int count ) const;
  bool unpack( const char * buf, const int size, int & pos, char * val, const int count ) const;
  bool unpack( const char * buf, const int size, int & pos, double * val, const int count ) const;
  bool unpack( const char * buf, const int size, int & pos, long * val, const int count ) const;
  bool unpack( const char * buf, const int size, int & pos, std::size_t * val, const int count ) const;

  // Communicator Barrier function.
  // A no-op for a serial communicator.  For MPI, it causes each processor in
  // the communicator to wait until all processors have arrived.
  void barrier() const;

protected:

  // MPI communicator.
  MPI_Comm mpiComm_;
  bool mpiCommOwned_;

  mutable MPI_Status status_;
  std::list<MPI_Request> request_;
  
  // Note number of processors 
  int numProc_;

private:

  bool operator==(const N_PDS_Comm &right) const { return false; }
  bool operator!=(const N_PDS_Comm &right) const { return false; }

private:

};

//-----------------------------------------------------------------------------
// Function      : N_PDS_Comm::comm
// Purpose       : Returns MPI_Comm object
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 03/13/00
//-----------------------------------------------------------------------------
inline MPI_Comm N_PDS_MPIComm::comm() const
{
  return mpiComm_;
}

#endif

#endif


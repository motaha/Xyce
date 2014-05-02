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

//-----------------------------------------------------------------------------
// Filename       : $RCSfile: N_PDS_SerialParComm.h,v $
//
// Purpose        : Specification file for the serial communication
//                  class that must work within a parallel environment
//                  for Xyce.
//
// Special Notes  :
//
// Creator        : Richard Schiek, SNL, Parallel Computational Sciences
//
// Creation Date  : 08/30/06
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.11 $
//
// Revision Date  : $Date: 2014/02/24 23:49:25 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_PDS_SerialParComm_h
#define Xyce_N_PDS_SerialParComm_h

#include <N_PDS_Comm.h>

class N_PDS_ParComm;

//-----------------------------------------------------------------------------
// Class         : N_PDS_SerialParComm
// Purpose       : Serial communication class for Xyce.
// Special Notes :
// Creator       : Richard Schiek, SNL, Parallel Compuational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
class N_PDS_SerialParComm : public N_PDS_Comm
{

public:

  // Default constructor.
  N_PDS_SerialParComm( N_PDS_ParComm * aParComm );

  //Destructor
  ~N_PDS_SerialParComm();

  // Copy constructor - reuse the same Petra_Comm object.
  N_PDS_SerialParComm(const N_PDS_SerialParComm & right);

  // Assignment operator - reuse the same Petra_Comm object.
  N_PDS_SerialParComm & operator = (const N_PDS_SerialParComm & right);

  // Cloning
  N_PDS_Comm * clone() const { return new N_PDS_SerialParComm(* this); }

    Xyce::Parallel::Machine comm() const {
      return 0;
    }
    
  // Get my processor ID.
  int procID() const { return processorID_; }

  // Get the total number of processors in this communicator.
  int numProc() const { return 1; }

  // Is this a serial run?
  bool isSerial() const { return true; }

  // Get the last processor flag (by default true in serial)
  // NOTE:  This is used for generating augmented linear systems.
  bool isLastProc() const { return true; }

  // Wrappers for Petra_Comm functionality.
  bool scanSum(const double * vals, double * sums, const int & count) const;
  bool sumAll(const double * vals, double * sums, const int & count) const;
  bool maxAll(const double * vals, double * maxs, const int & count) const;
  bool minAll(const double * vals, double * mins, const int & count) const;

  bool scanSum(const int * vals, int * sums, const int & count) const;
  bool sumAll(const int * vals, int * sums, const int & count) const;
  bool maxAll(const int * vals, int * maxs, const int & count) const;
  bool minAll(const int * vals, int * mins, const int & count) const;

  // Wrapper for MPI broadcasts.
  bool bcast(int * val, const int & count, const int & root) const { return true; }
  bool bcast(char * val, const int & count, const int & root) const { return true; }
  bool bcast(double * val, const int & count, const int & root) const { return true; }
  bool bcast(long * val, const int & count, const int & root) const { return true; }

  // Wrappers for MPI Sends.
  bool send(const int * val, const int & count, const int & dest) const { return true; }
  bool send(const char * val, const int & count, const int & dest) const { return true; }
  bool send(const double * val, const int & count, const int & dest) const { return true; }
  bool send(const long * val, const int & count, const int & dest) const { return true; }

  // Wrappers for MPI Recvs.
  bool recv(int * val, const int & count, const int & src) const { return true; }
  bool recv(char * val, const int & count, const int & src) const { return true; }
  bool recv(double * val, const int & count, const int & src) const { return true; }
  bool recv(long * val, const int & count, const int & src) const { return true; }

  // Wrappers for MPI pack and unpack functionality.
  bool rSend( const int * val, const int & count, const int & dest ) const { return true; }
  bool rSend( const char * val, const int & count, const int & dest ) const { return true; }
  bool rSend( const double * val, const int & count, const int & dest ) const { return true; }
  bool rSend( const long * val, const int & count, const int & dest ) const { return true; }
                                                                                           
  // MPI_Irecv wrappers
  bool iRecv( int * val, const int & count, const int & src ) const { return true; }
  bool iRecv( char * val, const int & count, const int & src ) const { return true; }
  bool iRecv( double * val, const int & count, const int & src ) const { return true; }
  bool iRecv( long * val, const int & count, const int & src ) const { return true; }
                                                                                           
  // MPI_Waitall wrappers
  bool waitAll() { return true; }
                                                                                           
  // Wrappers for MPI pack and unpack functionality.
  bool pack(const int * val, const int count, char * buf, const int size, int & pos) const;
  bool pack(const char * val, const int count, char * buf, const int size, int & pos) const;
  bool pack(const double * val, const int count, char * buf, const int size, int & pos) const;
  bool pack(const long * val, const int count, char * buf, const int size, int & pos) const;
  bool pack(const std::size_t * val, const int count, char * buf, const int size, int & pos) const;

  bool unpack(const char * buf, const int size, int & pos, int * val, const int count) const;
  bool unpack(const char * buf, const int size, int & pos, char * val, const int count) const;
  bool unpack(const char * buf, const int size, int & pos, double * val, const int count) const;
  bool unpack(const char * buf, const int size, int & pos, long * val, const int count) const;
  bool unpack(const char * buf, const int size, int & pos, std::size_t * val, const int count) const;

  Epetra_Comm * petraComm() { return petraComm_; }

#ifdef Xyce_PARALLEL_MPI
  N_PDS_ParComm * parComm() { return backgroundN_PDS_ParComm_; }
#endif

  // Communicator Barrier function.
  // A no-op for a serial communicator.  For MPI, it causes each processor in
  // the communicator to wait until all processors have arrived.
  void barrier() const { }
  
#ifdef Xyce_PARALLEL_MPI
  // set the background MPI communicator.  This is used when a parallel 
  // Xyce build is being run on one processor.  Xyce must still call 
  // MPI_Init() but the rest of Xyce will not know that it's running in parallel
  // void setN_PDS_ParComm( N_PDS_ParComm * theN_PDS_ParComm );
#endif

private:

  // Serial-run flag.
  int processorID_;
  bool isSerial_;

  Epetra_Comm * petraComm_;
  bool petraCommOwned_;
  
#ifdef Xyce_PARALLEL_MPI
  N_PDS_ParComm * backgroundN_PDS_ParComm_;
#endif

};

#endif

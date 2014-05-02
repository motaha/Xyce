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
// Filename       : $RCSfile: N_PDS_MPIComm.C,v $
//
// Purpose        : Implementation file for the abstract parallel communication
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
// Revision Number: $Revision: 1.30 $
//
// Revision Date  : $Date: 2014/02/24 23:49:25 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


#ifdef Xyce_PARALLEL_MPI

// ---------- Standard Includes ----------
#include <N_UTL_Misc.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <string.h>
#include <iostream>
#include <mpi.h>

#ifdef Xyce_Dakota

#include <string>
#include <stdlib.h>
#include <iostream>
#include <sstream>

#endif
// ----------   Xyce Includes   ----------

#include <N_PDS_MPIComm.h>
#include <N_ERH_ErrorMgr.h>

// ----------  Other Includes   ----------

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::N_PDS_MPIComm
// Purpose       : Default constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
N_PDS_MPIComm::N_PDS_MPIComm( int iargs, char * cargs[] )
 : mpiCommOwned_(true)
{
  // Set the MPI communicator
  mpiComm_ = MPI_COMM_WORLD;

  // Error string
  const std::string errorMsg( "N_PDS_MPIComm::initMPI - MPI_Init failed.");
  
#ifdef HAVE_UNISTD_H
  //
  // In some implementations of MPICH using the ch_p4 communicator,
  // calling mpirun -np 1 Xyce , or just running Xyce (which implies -np 1)
  // causes MPI_Init() to reset the current working directory to the location
  // of the Xyce binary.  Here we record the current working directory
  // and reset it if it has changed by calling MPI_Init().  This requires
  // unistd.h so I'm enclosing this test in ifdefs so it doesn't break
  // other builds unnecessarily 
  //

  const int maxPathLen = 4096;
  char originalCWD[maxPathLen];
  char * getcwdRV1 = getcwd( originalCWD, maxPathLen );
  
  if ( MPI_SUCCESS != MPI_Init( &iargs, &cargs ) )
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, errorMsg);
  
  char currentCWD[maxPathLen];
  char * getcwdRV2 = getcwd( currentCWD, maxPathLen );
  
  if( (getcwdRV1 != NULL) && (getcwdRV2 != NULL) )
  {
    if( strcmp( originalCWD, currentCWD ) != 0 )
    {
      if( chdir( originalCWD ) != 0 )
      {
        // changing the directory back may have failed.  Issue a warning
        // and try and continue.
        const std::string cdErrorMsg1( "N_PDS_MPIComm::initMPI Resetting working directory failed.  Trying to continue.");
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_WARNING, cdErrorMsg1 );
      }
    }
  }
  else
  {
    const std::string cdErrorMsg2( "N_PDS_MPIComm::initMPI Couldn't get working directory.  Trying to continue.");
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_WARNING, cdErrorMsg2 );
  }
#endif 

#ifdef Xyce_Dakota
  // if -ppp <number> was specified, then split the communicator
  int procPerProblem = 0;
  for( int i=0; i<iargs; i++ )
  {
    std::string anArg( cargs[i] );
    if( anArg.compare( "-ppp" ) == 0 )
    {
      // next arg should be number of processors per problem
      i++;
      int nextI = i;
      if( nextI < iargs )
      { 
        std::string argValue( cargs[ nextI ] );
        std::stringstream iost;
        iost << argValue;
        iost >> procPerProblem;
      }
      else
      {
        // issue error that "-ppp " was not followed by a number
      }
      break;
    }
  }
  
  if( procPerProblem > 0 )
  {
    
    // split mpi communicator
    int processID = procID();
    MPI_Comm myNewMpiComm;
    Xyce::dout() << "N_PDS_MPIComm on procID, " << processID << ", procPerProblem = " << procPerProblem << std::endl;
    const std::string errorMsgForSplit( "N_PDS_MPIComm::initMPI - MPI_Comm_split failed.");
    if( MPI_SUCCESS != MPI_Comm_split( MPI_COMM_WORLD, processID, procPerProblem, &myNewMpiComm ) )
      N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, errorMsg);
    mpiComm_ = myNewMpiComm;  
  }
  numProc_ = numProc();

#endif


}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::N_PDS_MPIComm
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
N_PDS_MPIComm::N_PDS_MPIComm( const MPI_Comm mComm )
 : mpiComm_(mComm),
   mpiCommOwned_(false)
{
  numProc_ = numProc();
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::N_PDS_MPIComm
// Purpose       : Copy constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
N_PDS_MPIComm::N_PDS_MPIComm( const N_PDS_MPIComm & right )
 : mpiComm_(right.mpiComm_),
   mpiCommOwned_(false)
{
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::N_PDS_MPIComm
// Purpose       : Assignment
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
N_PDS_MPIComm & N_PDS_MPIComm::operator=( const N_PDS_MPIComm & right )
{
  mpiComm_ = right.mpiComm_;
  mpiCommOwned_ = false;

  return *this;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::~N_PDS_MPIComm
// Purpose       : Destructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
N_PDS_MPIComm::~N_PDS_MPIComm()
{
  if( mpiCommOwned_ )
    MPI_Finalize();
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::numProc
// Purpose       : Returns the number of processors in the current
//                 configuration.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 03/13/00
//-----------------------------------------------------------------------------
int N_PDS_MPIComm::numProc() const
{

  int size = 1;

  // Error string
  const std::string error_msg(
    "N_PDS_MPIComm::numProc -  MPI_Comm_size failed.");

  // Get the machine size.
  if (MPI_SUCCESS != MPI_Comm_size(mpiComm_, &size) )
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, error_msg);

  return size;

}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::procID
// Purpose       : Returns "this" processor's id.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 03/13/00
//-----------------------------------------------------------------------------
int N_PDS_MPIComm::procID() const
{
  int id = 0;

  // Error string
  const std::string error_msg(
    "N_PDS_Comm::getProcID -  MPI_Comm_rank failed.");

  // Get the machine size.
  if (MPI_SUCCESS != MPI_Comm_rank(mpiComm_, &id) )
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, error_msg);

  return id;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::barrier
// Purpose       : Communicator Barrier function.  A no-op for a serial
//                 communicator.  For MPI, it causes each processor in the
//                 communicator to wait until all processors have arrived.
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Computational Sciences
// Creation Date : 09/18/00
//-----------------------------------------------------------------------------
void N_PDS_MPIComm::barrier() const
{
  MPI_Barrier(mpiComm_);
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::bcast
// Purpose       : for INTs
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::bcast( int * val, const int & count, const int & root ) const
{
  MPI_Bcast( val, const_cast<int &> (count), MPI_INT,
		const_cast<int &> (root), mpiComm_ );
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::bcast
// Purpose       : for LNGs
// Special Notes :
// Scope         : Public
// Creator       : Derek Barnes, SNL, Parallel Computational Sciences
// Creation Date : 08/01/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::bcast( long * val, const int & count, const int & root ) const
{
  MPI_Bcast( val, const_cast<int &> (count), MPI_INT,
		const_cast<int &> (root), mpiComm_ );
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::bcast
// Purpose       : for CHARs
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::bcast( char * val, const int & count, const int & root ) const
{
  MPI_Bcast( val, const_cast<int &> (count), MPI_CHAR,
		const_cast<int &> (root), mpiComm_ );
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::bcast
// Purpose       : for DBLEs
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::bcast( double * val, const int & count, const int & root ) const
{
  MPI_Bcast( val, const_cast<int &> (count), MPI_INT,
		const_cast<int &> (root), mpiComm_ );
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::send
// Purpose       : for INTs
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::send( const int * val, const int & count, const int & dest ) const
{
  MPI_Send( const_cast<int *> (val), const_cast<int &> (count), MPI_INT,
		const_cast<int &> (dest), 0, mpiComm_ );
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::send
// Purpose       : for INTs
// Special Notes :
// Scope         : Public
// Creator       : Derek Barnes, SNL, Parallel Computational Sciences
// Creation Date : 08/01/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::send( const long * val, const int & count, const int & dest ) const
{
  MPI_Send( const_cast<long *> (val), const_cast<int &> (count), MPI_INT,
		const_cast<int &> (dest), 0, mpiComm_ );
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::send
// Purpose       : for CHARs
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::send( const char * val, const int & count, const int & dest ) const
{
  MPI_Send( const_cast<char *> (val), const_cast<int &> (count), MPI_CHAR,
		const_cast<int &> (dest), 0, mpiComm_ );
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::send
// Purpose       : for DBLEs
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::send( const double * val, const int & count, const int & dest ) const
{
  MPI_Send( const_cast<double *> (val), const_cast<int &> (count), MPI_DOUBLE,
		const_cast<int &> (dest), 0, mpiComm_ );
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::recv
// Purpose       : for INTs
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::recv( int * val, const int & count, const int & src ) const
{
  MPI_Recv( val, const_cast<int &> (count), MPI_INT, const_cast<int &> (src),
		0, mpiComm_, &status_ );
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::recv
// Purpose       : for LNGs
// Special Notes :
// Scope         : Public
// Creator       : Derek Barnes, SNL, Parallel Computational Sciences
// Creation Date : 08/01/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::recv( long * val, const int & count, const int & src )const
{
  MPI_Recv( val, const_cast<int &> (count), MPI_INT, const_cast<int &> (src),
		0, mpiComm_, &status_ );
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::recv
// Purpose       : for CHARs
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::recv( char * val, const int & count, const int & src ) const
{
  MPI_Recv( val, const_cast<int &> (count), MPI_CHAR, const_cast<int &> (src),
		0, mpiComm_, &status_ );
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::recv
// Purpose       : for DBLEs
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::recv( double * val, const int & count, const int & src ) const
{
  MPI_Recv( val, const_cast<int &> (count), MPI_DOUBLE, const_cast<int &> (src),
		0, mpiComm_, &status_ );
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::rSend
// Purpose       : for INTs
// Special Notes :
// Scope         : Public
// Creator       : Derek Barnes, SNL, Parallel Computational Sciences
// Creation Date : 08/01/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::rSend( const int * val, const int & count, const int & dest ) const
{
  MPI_Rsend( const_cast<int *> (val), const_cast<int &> (count), MPI_INT,
            const_cast<int &> (dest), 0, mpiComm_ );
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::rSend
// Purpose       : for LONGs
// Special Notes :
// Scope         : Public
// Creator       : Derek Barnes, SNL, Parallel Computational Sciences
// Creation Date : 08/01/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::rSend( const long * val, const int & count, const int & dest ) const
{
  MPI_Rsend( const_cast<long *> (val), const_cast<int &> (count), MPI_INT,
            const_cast<int &> (dest), 0, mpiComm_ );
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::rSend
// Purpose       : for CHARs
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::rSend( const char * val, const int & count, const int & dest ) const
{
  MPI_Rsend( const_cast<char *> (val), const_cast<int &> (count), MPI_CHAR,
		const_cast<int &> (dest), 0, mpiComm_ );
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::rSend
// Purpose       : for DBLEs
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::rSend( const double * val, const int & count, const int & dest ) const
{
  MPI_Rsend( const_cast<double *> (val), const_cast<int &> (count), MPI_DOUBLE,
		const_cast<int &> (dest), 0, mpiComm_ );
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::iRecv
// Purpose       : for INTs
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::iRecv( int * val, const int & count, const int & src )
{
  request_.push_front( MPI_Request() );
  MPI_Irecv( val, const_cast<int &> (count), MPI_INT, const_cast<int &> (src),
		0, mpiComm_, &request_.front() );
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::iRecv
// Purpose       : for LNGs
// Special Notes :
// Scope         : Public
// Creator       : Derek Barnes, SNL, Parallel Computational Sciences
// Creation Date : 08/01/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::iRecv( long * val, const int & count, const int & src )
{
  request_.push_front( MPI_Request() );
  MPI_Irecv( val, const_cast<int &> (count), MPI_INT, const_cast<int &> (src),
		0, mpiComm_, &request_.front() );
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::iRecv
// Purpose       : for CHARs
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::iRecv( char * val, const int & count, const int & src )
{
  request_.push_front( MPI_Request() );
  MPI_Irecv( val, const_cast<int &> (count), MPI_CHAR, const_cast<int &> (src),
		0, mpiComm_, &request_.front() );
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::Irecv
// Purpose       : for DBLEs
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::iRecv( double * val, const int & count, const int & src )
{
  request_.push_front( MPI_Request() );
  MPI_Irecv( val, const_cast<int &> (count), MPI_DOUBLE, const_cast<int &> (src),
		0, mpiComm_, &request_.front() );
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::waitAll
// Purpose       : 
// Special Notes :
// Scope         : Public
// Creator       : Dave Shirley, PSSI
// Creation Date : 02/09/06
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::waitAll()
{
  std::list<MPI_Request>::iterator r=request_.begin();
  std::list<MPI_Request>::iterator r_end = request_.end();

  for ( ; r!=r_end ; ++r)
  {
    MPI_Wait( &(*r), &status_ );
  }
  request_.clear();

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::pack
// Purpose       : for INTs
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::pack( const int * val, const int count, char * buf, const int size, int & pos ) const
{
#ifdef Xyce_DEBUG_PARALLEL
  MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
#endif
  int err=MPI_SUCCESS;
  err = MPI_Pack( const_cast<int *> (val), count, MPI_INT,
	   	  buf, size, &pos, mpiComm_ );
#ifdef Xyce_DEBUG_PARALLEL
  if (err != MPI_SUCCESS) {
    std::cerr << "Processor " << procID() << " encountered an error ("<< err << ") calling MPI_Pack on an MPI_INT" << std::endl;
    std::cerr << "count = " << count << ", pos = "<< pos << ", size = " << size << std::endl;
    char error_string[BUFSIZ];
    int length_of_error_string;

    std::cerr << "Packing array: ";
    for (int i=0; i<count; i++) 
      std::cerr << val[i];
    std::cerr << std::endl;
    MPI_Error_string(err, error_string, &length_of_error_string);
    fprintf(stderr, "%3d: %s\n", procID(), error_string);
    return false;
  }
#endif
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::pack
// Purpose       : for LNGs
// Special Notes :
// Scope         : Public
// Creator       : Derek Barnes, SNL, Parallel Computational Sciences
// Creation Date : 08/01/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::pack( const long * val, const int count, char * buf, const int size, int & pos ) const
{
#ifdef Xyce_DEBUG_PARALLEL
  MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
#endif
  int err=MPI_SUCCESS;
  err = MPI_Pack( const_cast<long *> (val), count, MPI_INT,
		buf, size, &pos, mpiComm_ );
#ifdef Xyce_DEBUG_PARALLEL
  if (err != MPI_SUCCESS) {
    std::cerr << "Processor " << procID() << " encountered an error ("<< err << ") calling MPI_Pack on an MPI_LONG_INT" << std::endl;
    std::cerr << "count = " << count << ", pos = "<< pos << ", size = " << size << std::endl;
    char error_string[BUFSIZ];
    int length_of_error_string;

    std::cerr << "Packing array: ";
    for (int i=0; i<count; i++) 
      std::cerr << val[i];
    std::cerr << std::endl;
    MPI_Error_string(err, error_string, &length_of_error_string);
    fprintf(stderr, "%3d: %s\n", procID(), error_string);
    return false;
  }
#endif
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::pack
// Purpose       : for CHARs
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::pack( const char * val, const int count, char * buf, const int size, int & pos ) const
{
#ifdef Xyce_DEBUG_PARALLEL
  MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
#endif
  int err=MPI_SUCCESS;
  err = MPI_Pack( const_cast<char *> (val), count, MPI_CHAR,
		buf, size, &pos, mpiComm_ );
#ifdef Xyce_DEBUG_PARALLEL
  if (err != MPI_SUCCESS) {
    std::cerr << "Processor " << procID() << " encountered an error (" << err << ") calling MPI_Pack on an MPI_CHAR" << std::endl;
    std::cerr << "count = " << count << ", pos = "<< pos << ", size = " << size << std::endl;
    char error_string[BUFSIZ];
    int length_of_error_string;

    std::cerr << "Packing array: ";
    for (int i=0; i<count; i++) 
      std::cerr << val[i];
    std::cerr << std::endl;
    MPI_Error_string(err, error_string, &length_of_error_string);
    for (int i=0; i<count; i++) 
      std::cerr << val[i];
    std::cerr << std::endl;
    fprintf(stderr, "%3d: %s\n", procID(), error_string);
    return false;
  }
#endif
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::pack
// Purpose       : for DBLEs
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::pack( const double * val, const int count, char * buf, const int size, int & pos ) const
{
#ifdef Xyce_DEBUG_PARALLEL
  MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
#endif
  int err=MPI_SUCCESS;
  err = MPI_Pack( const_cast<double *> (val), count, MPI_DOUBLE,
		buf, size, &pos, mpiComm_ );
#ifdef Xyce_DEBUG_PARALLEL
  if (err != MPI_SUCCESS) {
    std::cerr << "Processor " << procID() << " encountered an error (" << err << ") calling MPI_Pack on an MPI_DOUBLE" << std::endl;
    std::cerr << "count = " << count << ", pos = "<< pos << ", size = " << size << std::endl;
    char error_string[BUFSIZ];
    int length_of_error_string;

    std::cerr << "Packing array: ";
    for (int i=0; i<count; i++) 
      std::cerr << val[i];
    std::cerr << std::endl;
    MPI_Error_string(err, error_string, &length_of_error_string);
    fprintf(stderr, "%3d: %s\n", procID(), error_string);
    return false;
  }
#endif
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::pack
// Purpose       : for size_t's 
// Special Notes :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 06/26/2013
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::pack( const std::size_t * val, const int count, char * buf, const int size, int & pos ) const
{
#ifdef Xyce_DEBUG_PARALLEL
  MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
#endif
  int err=MPI_SUCCESS;
  err = MPI_Pack( const_cast<size_t*> (val), count, MPI_INT,
		buf, size, &pos, mpiComm_ );
#ifdef Xyce_DEBUG_PARALLEL
  if (err != MPI_SUCCESS) {
    std::cerr << "Processor " << procID() << " encountered an error ("<< err << ") calling MPI_Pack on an MPI_LONG_INT" << std::endl;
    std::cerr << "count = " << count << ", pos = "<< pos << ", size = " << size << std::endl;
    char error_string[BUFSIZ];
    int length_of_error_string;

    std::cerr << "Packing array: ";
    for (int i=0; i<count; i++) 
      std::cerr << val[i];
    std::cerr << std::endl;
    MPI_Error_string(err, error_string, &length_of_error_string);
    fprintf(stderr, "%3d: %s\n", procID(), error_string);
    return false;
  }
#endif
  return true;
}



//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::unpack
// Purpose       : for INTs
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::unpack( const char * buf, const int size, int & pos, int * val, const int count ) const
{
  MPI_Unpack( const_cast<char *> (buf), size, &pos, val,
		count, MPI_INT, mpiComm_ );
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::unpack
// Purpose       : for LNGs
// Special Notes :
// Scope         : Public
// Creator       : Derek Barnes, SNL, Parallel Computational Sciences
// Creation Date : 08/01/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::unpack( const char * buf, const int size, int & pos, long * val, const int count ) const
{
  MPI_Unpack( const_cast<char *> (buf), size, &pos, val,
		count, MPI_INT, mpiComm_ );
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::unpack
// Purpose       : for CHARs
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::unpack( const char * buf, const int size, int & pos, char * val, const int count ) const
{
  MPI_Unpack( const_cast<char *> (buf), size, &pos, val,
		count, MPI_CHAR, mpiComm_ );
  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::unpack
// Purpose       : for DBLEs
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::unpack( const char * buf, const int size, int & pos, double * val, const int count ) const
{
  MPI_Unpack( const_cast<char *> (buf), size, &pos, val,
		count, MPI_DOUBLE, mpiComm_ );
  return true;
}


//-----------------------------------------------------------------------------
// Function      : N_PDS_MPIComm::unpack
// Purpose       : for size_t's 
// Special Notes :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 06/26/2013
//-----------------------------------------------------------------------------
bool N_PDS_MPIComm::unpack( const char * buf, const int size, int & pos, std::size_t * val, const int count ) const
{
  MPI_Unpack( const_cast<char *> (buf), size, &pos, val,
		count, MPI_INT, mpiComm_ );
  return true;
}

#endif


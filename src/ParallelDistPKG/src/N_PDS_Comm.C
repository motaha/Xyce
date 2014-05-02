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
// Filename       : $RCSfile: N_PDS_Comm.C,v $
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
// Revision Number: $Revision: 1.18 $
//
// Revision Date  : $Date: 2014/02/24 23:49:25 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#ifdef Xyce_PARALLEL_MPI
#include <mpi.h>
#endif

#include <N_PDS_Comm.h>
#include <N_ERH_ErrorMgr.h>

#include <Epetra_Comm.h>

namespace Xyce {
namespace Parallel {

//-----------------------------------------------------------------------------
// Function      : Communicator::Communicator
// Purpose       : Default constructor
// Special Notes :
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 03/13/00
//-----------------------------------------------------------------------------
Communicator::Communicator()
  : commOwned_(true)
{
  // Initialize
  isSerial_ = true;

#ifdef Xyce_PARALLEL_MPI

  // Set the MPI communicator
  mpiComm_ = MPI_COMM_WORLD;

  // Set the serial flag based upon the machine size.
  int size = 1;

  MPI_Comm_size(mpiComm_, &size);
  isSerial_ = (size == 1);
#endif

  // Call the Petra constructor for the true Petra comm.
#ifdef Xyce_PARALLEL_MPI
  libComm_ = new Epetra_Comm(mpiComm_);
#else
  libComm_ = new Epetra_Comm();
#endif

}

//-----------------------------------------------------------------------------
// Function      : Communicator::~Communicator
// Purpose       : Destructor
// Special Notes : Virtual
// Scope         : Public
// Creator       : Scott A. Hutchinson, SNL, Parallel Computational Sciences
// Creation Date : 03/13/00
//-----------------------------------------------------------------------------
Communicator::~Communicator()
{
  if( commOwned_ )
    if( libComm_ != NULL && libComm_ != 0 ) delete libComm_;

#ifdef Xyce_PARALLEL_MPI
  mpiComm_ = MPI_COMM_NULL;
#endif
}

#ifdef Xyce_PARALLEL_MPI

//-----------------------------------------------------------------------------
// Function      : Communicator::scanSum
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool Communicator::scanSum( double * vals, double * sums, int count ) const
{
  return ( libComm_->ScanSum( vals, sums, count ) == 0 );
}

//-----------------------------------------------------------------------------
// Function      : Communicator::sumAll
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool Communicator::sumAll( double * vals, double * sums, int count ) const
{
  return ( libComm_->SumAll( vals, sums, count ) == 0 );
}

//-----------------------------------------------------------------------------
// Function      : Communicator::maxAll
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool Communicator::maxAll( double * vals, double * maxs, int count ) const
{
  return ( libComm_->MaxAll( vals, maxs, count ) == 0 );
}

//-----------------------------------------------------------------------------
// Function      : Communicator::minAll
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool Communicator::minAll( double * vals, double * mins, int count ) const
{
  return ( libComm_->MinAll( vals, maxs, count ) == 0 );
}
#endif

//-----------------------------------------------------------------------------
// Function      : Communicator::scanSum
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Dave Shirley, PSSI
// Creation Date : 06/18/05
//-----------------------------------------------------------------------------
bool Communicator::scanSum( int * vals, int * sums, int count ) const
{
  return ( libComm_->ScanSum( vals, sums, count ) == 0 );
}

//-----------------------------------------------------------------------------
// Function      : Communicator::sumAll
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Dave Shirley, PSSI
// Creation Date : 06/18/05
//-----------------------------------------------------------------------------
bool Communicator::sumAll( int * vals, int * sums, int count ) const
{
  return ( libComm_->SumAll( vals, sums, count ) == 0 );
}

//-----------------------------------------------------------------------------
// Function      : Communicator::maxAll
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Dave Shirley, PSSI
// Creation Date : 06/18/05
//-----------------------------------------------------------------------------
bool Communicator::maxAll( int * vals, int * maxs, int count ) const
{
  return ( libComm_->MaxAll( vals, maxs, count ) == 0 );
}

//-----------------------------------------------------------------------------
// Function      : Communicator::minAll
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Dave Shirley, PSSI
// Creation Date : 06/18/05
//-----------------------------------------------------------------------------
bool Communicator::minAll( int * vals, int * mins, int count ) const
{
  return ( libComm_->MinAll( vals, maxs, count ) == 0 );
}

} // namespace Parallel
} // namespace Xyce

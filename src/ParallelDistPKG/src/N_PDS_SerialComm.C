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
// Filename       : $RCSfile: N_PDS_SerialComm.C,v $
//
// Purpose        : Implementation file for the serial communication
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
// Revision Number: $Revision: 1.16 $
//
// Revision Date  : $Date: 2014/02/24 23:49:25 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------

#include <N_PDS_SerialComm.h>

#include <N_ERH_ErrorMgr.h>

// ----------  Other Includes   ----------

#include <Epetra_SerialComm.h>

//-----------------------------------------------------------------------------
// Function      : N_PDS_SerialComm::N_PDS_SerialComm
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
N_PDS_SerialComm::N_PDS_SerialComm()
  :
  isSerial_(true),
  petraComm_(0),
  petraCommOwned_(true)
{
  petraComm_ = new Epetra_SerialComm();
  petraCommOwned_ = true;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_SerialComm::N_PDS_SerialComm
// Purpose       : Copy constructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
N_PDS_SerialComm::N_PDS_SerialComm(const N_PDS_SerialComm &right)
  :
  isSerial_(right.isSerial_) ,
  petraComm_(right.petraComm_),
  petraCommOwned_(false)
{
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_SerialComm::~N_PDS_SerialComm
// Purpose       : Destructor
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoeksta, SNL, Parallel Computational Sciences
// Creation Date : 06/26/01
//-----------------------------------------------------------------------------
N_PDS_SerialComm::~N_PDS_SerialComm()
{
  if( petraCommOwned_ )
    if( petraComm_ ) delete petraComm_;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_Comm::scanSum
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_SerialComm::scanSum( const double * vals, double * sums,
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
bool N_PDS_SerialComm::sumAll( const double * vals, double * sums,
		const int & count ) const
{
  return ( petraComm_->SumAll( const_cast<double *> (vals), sums,
		const_cast<int &> (count) ) == 0 );
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_SerialComm::maxAll
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_SerialComm::maxAll( const double * vals, double * maxs,
		const int & count ) const
{
  return ( petraComm_->MaxAll( const_cast<double *> (vals), maxs,
		const_cast<int &> (count) ) == 0 );
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_SerialComm::minAll
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL, Parallel Computational Sciences
// Creation Date : 08/30/01
//-----------------------------------------------------------------------------
bool N_PDS_SerialComm::minAll( const double * vals, double * mins,
		const int & count ) const
{
  return ( petraComm_->MinAll( const_cast<double *> (vals), mins,
		const_cast<int &> (count) ) == 0 );
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_Comm::scanSum
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Dave Shirley, PSSI
// Creation Date : 06/18/05
//-----------------------------------------------------------------------------
bool N_PDS_SerialComm::scanSum( const int * vals, int * sums,
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
bool N_PDS_SerialComm::sumAll( const int * vals, int * sums,
		const int & count ) const
{
  return ( petraComm_->SumAll( const_cast<int *> (vals), sums,
		const_cast<int &> (count) ) == 0 );
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_SerialComm::maxAll
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Dave Shirley, PSSI
// Creation Date : 06/18/05
//-----------------------------------------------------------------------------
bool N_PDS_SerialComm::maxAll( const int * vals, int * maxs,
		const int & count ) const
{
  return ( petraComm_->MaxAll( const_cast<int *> (vals), maxs,
		const_cast<int &> (count) ) == 0 );
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_SerialComm::minAll
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Dave Shirley, PSSI
// Creation Date : 06/18/05
//-----------------------------------------------------------------------------
bool N_PDS_SerialComm::minAll( const int * vals, int * mins,
		const int & count ) const
{
  return ( petraComm_->MinAll( const_cast<int *> (vals), mins,
		const_cast<int &> (count) ) == 0 );
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_SerialComm::pack
// Purpose       : INTs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_SerialComm::pack( const int * val, const int count, char * buf,
		const int size, int & pos ) const
{
  int bC = count * sizeof(int);
  memcpy( (void*)(buf+pos), (void*)(val), bC );
  pos += bC;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_SerialComm::pack
// Purpose       : LNGs
// Special Notes :
// Scope         : Public
// Creator       : Derek Barnes, SNL, Parallel Computational Sciences
// Creation Date : 08/01/01
//-----------------------------------------------------------------------------
bool N_PDS_SerialComm::pack( const long * val, const int count, char * buf,
		const int size, int & pos ) const
{
  int bC = count * sizeof(long);
  memcpy( (void*)(buf+pos), (void*)(val), bC );
  pos += bC;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_SerialComm::pack
// Purpose       : CHARs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_SerialComm::pack( const char * val, const int count, char * buf,
		const int size, int & pos ) const
{
  int bC = count * sizeof(char);
  memcpy( (void*)(buf+pos), (void*)(val), bC );
  pos += bC;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_SerialComm::pack
// Purpose       : DBLEs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_SerialComm::pack( const double * val, const int count, char * buf,
		const int size, int & pos ) const
{
  int bC = count * sizeof(double);
  memcpy( (void*)(buf+pos), (void*)(val), bC );
  pos += bC;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_SerialComm::pack
// Purpose       : size_t's
// Special Notes :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 06/26/2013
//-----------------------------------------------------------------------------
bool N_PDS_SerialComm::pack( const std::size_t * val, const int count, char * buf,
		const int size, int & pos ) const
{
  int bC = count * sizeof(size_t);
  memcpy( (void*)(buf+pos), (void*)(val), bC );
  pos += bC;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_SerialComm::unpack
// Purpose       : INTs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_SerialComm::unpack( const char * buf, const int size, int & pos,
		int * val, const int count ) const
{
  int bC = count * sizeof(int);
  memcpy( (void*)(val), (void*)(buf+pos), bC );
  pos += bC;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_SerialComm::unpack
// Purpose       : LNGs
// Special Notes :
// Scope         : Public
// Creator       : Derek Barnes, SNL, Parallel Computational Sciences
// Creation Date : 08/01/01
//-----------------------------------------------------------------------------
bool N_PDS_SerialComm::unpack( const char * buf, const int size, int & pos,
		long * val, const int count ) const
{
  int bC = count * sizeof(long);
  memcpy( (void*)(val), (void*)(buf+pos), bC );
  pos += bC;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_SerialComm::unpack
// Purpose       : CHARs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_SerialComm::unpack( const char * buf, const int size, int & pos,
		char * val, const int count ) const
{
  int bC = count * sizeof(char);
  memcpy( (void*)(val), (void*)(buf+pos), bC );
  pos += bC;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_SerialComm::unpack
// Purpose       : DBLEs
// Special Notes :
// Scope         : Public
// Creator       : Robert J Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 02/26/01
//-----------------------------------------------------------------------------
bool N_PDS_SerialComm::unpack( const char * buf, const int size, int & pos,
		double * val, const int count ) const
{
  int bC = count * sizeof(double);
  memcpy( (void*)(val), (void*)(buf+pos), bC );
  pos += bC;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_SerialComm::unpack
// Purpose       : size_t's 
// Special Notes :
// Scope         : Public
// Creator       : Eric R. Keiter, SNL
// Creation Date : 06/26/2013
//-----------------------------------------------------------------------------
bool N_PDS_SerialComm::unpack( const char * buf, const int size, int & pos,
		std::size_t * val, const int count ) const
{
  int bC = count * sizeof(size_t);
  memcpy( (void*)(val), (void*)(buf+pos), bC );
  pos += bC;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_SerialComm::operator=
// Purpose       : Assignment operator.
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 06/27/01
//-----------------------------------------------------------------------------
N_PDS_SerialComm & N_PDS_SerialComm::operator=(const N_PDS_SerialComm &right)
{
  isSerial_ = right.isSerial_;
  petraComm_  = right.petraComm_;
  petraCommOwned_ = false;

  return *this;
}


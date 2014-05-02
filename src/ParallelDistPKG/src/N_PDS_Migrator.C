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
// Filename       : $RCSfile: N_PDS_Migrator.C,v $
//
// Purpose        : Static methods for migrating data
//
// Special Notes  :
//
// Creator        : Rob Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/20/01
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

// ----------   Xyce Includes   ----------

#include <N_PDS_Migrator.h>

#include <N_PDS_ParComm.h>
#include <N_PDS_MPIComm.h>
#include <N_PDS_CommFactory.h>

#ifdef Xyce_PARALLEL_MPI
#include <Epetra_MpiComm.h>
#include <Epetra_MpiDistributor.h>

#include <N_UTL_Packable.h>
#endif

#include <N_UTL_Functors.h>

// ----------   Other Includes   ----------

#include <algorithm>

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParDir::N_PDS_ParDir
// Purpose       : Constructor
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/01
//-----------------------------------------------------------------------------
N_PDS_Migrator::N_PDS_Migrator( N_PDS_Comm * comm)
 : pdsComm_(comm),
   commOwned_(false)
{
  if( !pdsComm_ )
  {
    pdsComm_ = N_PDS_CommFactory::create();
    commOwned_ = true;
  }
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_Migrator::~N_PDS_Migrator
// Purpose       : Destructor
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/01
//-----------------------------------------------------------------------------
N_PDS_Migrator::~N_PDS_Migrator()
{
  if( commOwned_ ) delete pdsComm_;
}

#ifdef Xyce_PARALLEL_MPI

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParDir::migratePackable
// Purpose       : Migrate packable objects
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/01
//-----------------------------------------------------------------------------
void N_PDS_Migrator::migratePackable( const std::vector<int> & procVec,
				      const std::vector<Packable *> & exportVec,
				      std::vector<Packable *> & importVec )
{
  int exportCount = procVec.size();

  Epetra_MpiDistributor distributor( *(dynamic_cast<Epetra_MpiComm*>(pdsComm_->petraComm())) );

  std::vector<int> assign( procVec );
  std::vector<Packable*> exportV( exportVec );

  if( !IsSorted( assign ) ) SortContainer2( assign, exportV );

  int max_size = 0;
  for( int i = 0; i < exportCount; ++i )
    max_size = Xycemax( max_size, exportV[i]->packedByteCount() );

  int importCount;
  distributor.CreateFromSends( exportCount, &(assign[0]), true, importCount );

  double d_max_size = max_size;
  double d_max_all;
  pdsComm_->maxAll( &d_max_size, &d_max_all, 1 );
  int max_all = (int)d_max_all;

  std::vector<char> exports( max_all * exportCount );
  char * imports = new char[ max_all * importCount ];

  int pos;
  for( int i = 0; i < exportCount; ++i )
  {
    pos = max_all * i;
    exportV[i]->pack( &(exports[0]), ( max_all * exportCount ), pos, pdsComm_ );
  }

  int importSize = max_all * importCount;
  distributor.Do( &(exports[0]), max_all, importSize, imports );

  importVec.resize( importCount );
  for( int i = 0; i < importCount; ++i )
  {
    pos = max_all * i;
    importVec[i] = exportVec[0]->instance();
    importVec[i]->unpack( &(imports[0]), ( max_all * importCount ), pos, pdsComm_ );
  }

  delete [] imports;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParDir::reverseMigratePackable
// Purpose       : Migrate packable objects
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/01
//-----------------------------------------------------------------------------
void N_PDS_Migrator::reverseMigratePackable( const std::vector<int> & procVec,
				      const std::vector<Packable *> & exportVec,
				      std::vector<Packable *> & importVec )
{
  int importCount = procVec.size();

  Epetra_MpiDistributor distributor(
                 *(dynamic_cast<Epetra_MpiComm*>(pdsComm_->petraComm())) );

  std::vector<int> assign( procVec );
  sort( assign.begin(), assign.end() );

  int exportCount;
  distributor.CreateFromSends( importCount, &(assign[0]), true, exportCount );

  if( exportCount != exportVec.size() )
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL,
     "Reverse Export Size Error! (Packable)\n" );

  int max_size = 0;
  for( int i = 0; i < exportCount; ++i )
    max_size = Xycemax( max_size, exportVec[i]->packedByteCount() );
  double d_max_size = max_size;
  double d_max_all;
  pdsComm_->maxAll( &d_max_size, &d_max_all, 1 );
  int max_all = (int)d_max_all;

  std::vector<char> exports( max_all * exportCount );
  char * imports = new char[ max_all * importCount ];

  int pos;
  for( int i = 0; i < exportCount; ++i )
  {
    pos = max_all * i;
    exportVec[i]->pack( &(exports[0]), ( max_all * exportCount ), pos, pdsComm_ );
  }

  int importSize = max_all * importCount;
  distributor.DoReverse( &(exports[0]), max_all, importSize, imports );

  importVec.resize( importCount );
  for( int i = 0; i < importCount; ++i )
  {
    pos = max_all * i;
    importVec[i] = exportVec[0]->instance();
    importVec[i]->unpack( &(imports[0]), ( max_all * importCount ), pos, pdsComm_ );
  }

  delete [] imports;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParDir::migrateString
// Purpose       : Migrate strings
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/01
//-----------------------------------------------------------------------------
void N_PDS_Migrator::migrateString( const std::vector<int> & procVec,
				    const std::vector<std::string> & exportVec,
				    std::vector<std::string> & importVec )
{
  int exportCount = procVec.size();

  Epetra_MpiDistributor distributor( *(dynamic_cast<Epetra_MpiComm*>(pdsComm_->petraComm())) );

  std::vector<int> assign( procVec );
  std::vector<std::string> exportV( exportVec );

  if( !IsSorted( assign ) ) SortContainer2( assign, exportV );

  int importCount;
  distributor.CreateFromSends( exportCount, &(assign[0]), true, importCount );

  int max_size = 0;
  for( int i = 0; i < exportCount; ++i )
    max_size = Xycemax( max_size, exportV[i].length() + sizeof(int) );
  double d_max_size = max_size;
  double d_max_all;
  pdsComm_->maxAll( &d_max_size, &d_max_all, 1 );
  int max_all = static_cast<int> (d_max_all);

  std::vector<char> exports( max_all * exportCount );
  char * imports = new char[ max_all * importCount ];

  int pos;
  int length;
  char * charBuf = new char[ max_all ];
  for( int i = 0; i < exportCount; ++i )
  {
    pos = max_all * i;
    length = exportV[i].length();
    pdsComm_->pack( &length, 1, &(exports[0]), max_all * exportCount, pos );
    exportVec[i].copy( charBuf, std::string::npos );
    pdsComm_->pack( charBuf, length, &(exports[0]), max_all * exportCount, pos );
  }

  int importSize = max_all * importCount;
  distributor.Do( &(exports[0]), max_all, importSize, imports );

  importVec.resize( importCount );
  for( int i = 0; i < importCount; ++i )
  {
    pos = max_all * i;
    pdsComm_->unpack( &(imports[0]), max_all * importCount, pos, &length, 1 );
    pdsComm_->unpack( &(imports[0]), max_all * importCount, pos, charBuf, length );
    importVec[i] = std::string( charBuf, length );
  }

  delete [] charBuf;

  delete [] imports;

}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParDir::migrateInt
// Purpose       : Migrate ints
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/01
//-----------------------------------------------------------------------------
void N_PDS_Migrator::migrateInt( const std::vector<int> & procVec,
			         const std::vector< std::pair<std::string,int> > & exportVec,
				 std::vector< std::pair<std::string,int> > & importVec )
{
  int exportCount = procVec.size();

  Epetra_MpiDistributor distributor( *(dynamic_cast<Epetra_MpiComm*>(pdsComm_->petraComm())) );

  std::vector<int> assign( procVec );
  std::vector< std::pair<std::string,int> > exportV( exportVec );

  if( !IsSorted( assign ) ) SortContainer2( assign, exportV );

  int importCount;
  distributor.CreateFromSends( exportCount, &(assign[0]), true, importCount );

  int max_size = 0;
  for( int i = 0; i < exportCount; ++i )
    max_size = Xycemax( max_size, exportV[i].first.length() + 2 * sizeof(int) );
  double d_max_size = max_size;
  double d_max_all;
  pdsComm_->maxAll( &d_max_size, &d_max_all, 1 );
  int max_all = static_cast<int> (d_max_all);

  std::vector<char> exports( max_all * exportCount );
  char * imports = new char[ max_all * importCount ];

  int pos, val;
  int length;
  char * charBuf = new char[ max_all ];
  for( int i = 0; i < exportCount; ++i )
  {
    pos = max_all * i;
    length = exportV[i].first.length();
    pdsComm_->pack( &length, 1, &(exports[0]), max_all * exportCount, pos );
    exportVec[i].first.copy( charBuf, std::string::npos );
    pdsComm_->pack( charBuf, length, &(exports[0]), max_all * exportCount, pos );
    val = exportV[i].second;
    pdsComm_->pack( &val, 1, &(exports[0]), max_all * exportCount, pos );
  }

  int importSize = max_all * importCount;
  distributor.Do( &(exports[0]), max_all, importSize, imports );

  importVec.resize( importCount );
  for( int i = 0; i < importCount; ++i )
  {
    pos = max_all * i;
    pdsComm_->unpack( &(imports[0]), max_all * importCount, pos, &length, 1 );
    pdsComm_->unpack( &(imports[0]), max_all * importCount, pos, charBuf, length );
    importVec[i].first = std::string( charBuf, length );
    pdsComm_->unpack( &(imports[0]), max_all * importCount, pos, &(importVec[i].second), 1 );
  }

  delete [] charBuf;
  delete [] imports;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParDir::reverseMigrateInt
// Purpose       : Migrate ints
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/01
//-----------------------------------------------------------------------------
void N_PDS_Migrator::reverseMigrateInt( const std::vector<int> & procVec,
			                const std::vector< std::pair<std::string,int> > & exportVec,
			                std::vector< std::pair<std::string,int> > & importVec )
{
  int pid = pdsComm_->procID();

  int importCount = procVec.size();

  Epetra_MpiDistributor distributor( *(dynamic_cast<Epetra_MpiComm*>(pdsComm_->petraComm())) );

  std::vector<int> assign( procVec );
  sort( assign.begin(), assign.end() );

  int exportCount;
  distributor.CreateFromSends( importCount, &(assign[0]), true, exportCount );

  if( exportCount != exportVec.size() )
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL,
     "Reverse Export Size Error! (Int)\n" );

  int max_size = 0;
  for( int i = 0; i < exportCount; ++i )
    max_size = Xycemax( max_size, exportVec[i].first.length() +
			2 * sizeof(int) );
  double d_max_size = max_size;
  double d_max_all;
  pdsComm_->maxAll( &d_max_size, &d_max_all, 1 );
  int max_all = static_cast<int>(d_max_all);

  std::vector<char> exports( max_all * exportCount );
  char * imports = new char[ max_all * importCount ];

  int pos, val, length;
  char * charBuf = new char[ max_all ];
  for( int i = 0; i < exportCount; ++i )
  {
    pos = max_all * i;
    length = exportVec[i].first.length();
    pdsComm_->pack( &length, 1, &(exports[0]), max_all * exportCount, pos );
    exportVec[i].first.copy( charBuf, std::string::npos );
    pdsComm_->pack( charBuf, length, &(exports[0]), max_all * exportCount, pos );
    val = exportVec[i].second;
    pdsComm_->pack( &val, 1, &(exports[0]), max_all * exportCount, pos );
  }

  int importSize = max_all * importCount;
  distributor.DoReverse( &(exports[0]), max_all, importSize, imports );

  importVec.resize( importCount );
  for( int i = 0; i < importCount; ++i )
  {
    pos = max_all * i;
    pdsComm_->unpack( &(imports[0]), max_all * importCount, pos, &length, 1 );
    pdsComm_->unpack( &(imports[0]), max_all * importCount, pos, charBuf, length );
    importVec[i].first = std::string( charBuf, length );
    pdsComm_->unpack( &(imports[0]), max_all * importCount, pos, &(importVec[i].second), 1 );
  }

  delete [] charBuf;

  delete [] imports;

}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParDir::migrateIntVec
// Purpose       : Migrate int vectors
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/01
//-----------------------------------------------------------------------------
void N_PDS_Migrator::migrateIntVec( const std::vector<int> & procVec,
	            const std::vector< std::pair< std::string,std::vector<int> > > & exportVec,
		    std::vector< std::pair< std::string,std::vector<int> > > & importVec )
{
  int exportCount = procVec.size();

  Epetra_MpiDistributor distributor( *(dynamic_cast<Epetra_MpiComm*>(pdsComm_->petraComm())) );

  std::vector<int> assign( procVec );
  std::vector< std::pair< std::string,std::vector<int> > > exportV( exportVec );

  if( !IsSorted( assign ) ) SortContainer2( assign, exportV );

  int importCount;
  distributor.CreateFromSends( exportCount, &(assign[0]), true, importCount );

  int max_size = 0;
  for( int i = 0; i < exportCount; ++i )
    max_size = Xycemax( max_size, exportVec[i].first.length() +
	(exportV[i].second.size() + 2) * sizeof(int) );
  double d_max_size = max_size;
  double d_max_all;
  pdsComm_->maxAll( &d_max_size, &d_max_all, 1 );
  int max_all = static_cast<int>(d_max_all);

  std::vector<char> exports( max_all * exportCount );
  char * imports = new char[ max_all * importCount ];

  int pos, length, val;
  char * charBuf = new char[ max_all ];
  for( int i = 0; i < exportCount; ++i )
  {
    pos = max_all * i;
    length = exportV[i].first.length();
    pdsComm_->pack( &length, 1, &(exports[0]), max_all * exportCount, pos );
    exportV[i].first.copy( charBuf, std::string::npos );
    pdsComm_->pack( charBuf, length, &(exports[0]), max_all * exportCount, pos );
    length = exportV[i].second.size();
    pdsComm_->pack( &length, 1, &(exports[0]), max_all * exportCount, pos );
    for( int j = 0; j < length; ++j )
    {
      val = exportV[i].second[j];
      pdsComm_->pack( &val, 1, &(exports[0]), max_all * exportCount, pos );
    }
  }

  int importSize = max_all * importCount;
  distributor.Do( &(exports[0]), max_all, importSize, imports );

  importVec.resize( importCount );
  for( int i = 0; i < importCount; ++i )
  {
    pos = max_all * i;
    pdsComm_->unpack( &(imports[0]), max_all * importCount, pos, &length, 1 );
    pdsComm_->unpack( &(imports[0]), max_all * importCount, pos, charBuf, length );
    importVec[i].first = std::string( charBuf, length );
    pdsComm_->unpack( &(imports[0]), max_all * importCount, pos, &length, 1 );
    importVec[i].second.resize( length );
    for( int j = 0; j < length; ++j )
    {
      pdsComm_->unpack( &(imports[0]), max_all * importCount, pos, &val, 1 );
      importVec[i].second[j] = val;
    }
  }

  delete [] charBuf;
  delete [] imports;
}

//-----------------------------------------------------------------------------
// Function      : N_PDS_ParDir::reverseMigrateIntVec
// Purpose       : Migrate int vectors
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 05/20/01
//-----------------------------------------------------------------------------
void N_PDS_Migrator::reverseMigrateIntVec( const std::vector<int> & procVec,
	            const std::vector< std::pair< std::string,std::vector<int> > > & exportVec,
		    std::vector< std::pair< std::string,std::vector<int> > > & importVec )
{
  int importCount = procVec.size();

  Epetra_MpiDistributor distributor( *(dynamic_cast<Epetra_MpiComm*>(pdsComm_->petraComm())) );

  std::vector<int> assign( procVec );
  sort( assign.begin(), assign.end() );

  int exportCount;
  distributor.CreateFromSends( importCount, &(assign[0]), true, exportCount );

  if( exportCount != exportVec.size() )
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL,
     "Reverse Export Size Error! (IntVec)\n" );

  int max_size = 0;
  for( int i = 0; i < exportCount; ++i )
    max_size = Xycemax( max_size, exportVec[i].first.length() +
		(exportVec[i].second.size() + 2) * sizeof(int) );
  double d_max_size = max_size;
  double d_max_all;
  pdsComm_->maxAll( &d_max_size, &d_max_all, 1 );
  int max_all = static_cast<int>(d_max_all);

  std::vector<char> exports( max_all * exportCount );
  char * imports = new char[ max_all * importCount ];

  int pos, length, val;
  char * charBuf = new char[ max_all ];
  for( int i = 0; i < exportCount; ++i )
  {
    pos = max_all * i;
    length = exportVec[i].first.length();
    pdsComm_->pack( &length, 1, &(exports[0]), max_all * exportCount, pos );
    exportVec[i].first.copy( charBuf, std::string::npos );
    pdsComm_->pack( charBuf, length, &(exports[0]), max_all * exportCount, pos );
    length = exportVec[i].second.size();
    pdsComm_->pack( &length, 1, &(exports[0]), max_all * exportCount, pos );
    for( int j = 0; j < length; ++j )
    {
      val = exportVec[i].second[j];
      pdsComm_->pack( &val, 1, &(exports[0]), max_all * exportCount, pos );
    }
  }

  int importSize = max_all * importCount;
  distributor.DoReverse(  &(exports[0]), max_all, importSize, imports );

  importVec.resize( importCount );
  for( int i = 0; i < importCount; ++i )
  {
    pos = max_all * i;
    pdsComm_->unpack( &(imports[0]), max_all * importCount, pos, &length, 1 );
    pdsComm_->unpack( &(imports[0]), max_all * importCount, pos, charBuf, length );
    importVec[i].first = std::string( charBuf, length );
    pdsComm_->unpack( &(imports[0]), max_all * importCount, pos, &length, 1 );
    importVec[i].second.resize( length );
    for( int j = 0; j < length; ++j )
    {
      pdsComm_->unpack( &(imports[0]), max_all * importCount, pos, &val, 1 );
      importVec[i].second[j] = val;
    }
  }

  delete [] charBuf;

  delete [] imports;
}

#endif


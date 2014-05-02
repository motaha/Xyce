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
// Filename       : $RCSfile: N_PDS_Directory.h,v $
//
// Purpose        : Generic templated distributed directory
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 08/07/03
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.14 $
//
// Revision Date  : $Date: 2014/02/24 23:49:25 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_PDS_Directory_h
#define Xyce_N_PDS_Directory_h

// ---------- Standard Includes ----------

#include <N_UTL_Misc.h>
#include <N_TOP_Misc.h> // To define NodeID for the topology package

#ifdef HAVE_CMATH
#include <cmath>
#else
#include <math.h>
#endif

#include <map>
#include <iostream>
#include <vector>

#include <Teuchos_RCP.hpp>
using Teuchos::RCP;
using Teuchos::rcp;

// ---------- Forward Declarations -------

// ----------   Xyce Includes   ----------

#include <N_ERH_ErrorMgr.h>

#include <N_UTL_Functors.h>

// ----------   Other Includes   ----------

namespace Xyce {
namespace Parallel {

//-----------------------------------------------------------------------------
// Class         : N_PDS_Directory
// Purpose       : Distributed directory for object lookup
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Compuational Sciences
// Creation Date : 05/18/01
//-----------------------------------------------------------------------------
template <typename KT, typename DT, class DH, class AC, class MG>
class Directory
{

public:

  typedef typename std::map< KT, RCP<DT> >         DataMap;
  typedef typename DataMap::iterator        DataMapIter;
  typedef typename DataMap::const_iterator  DataMapCIter;

  typedef typename std::multimap< KT, RCP<DT> >    DataRecvMap;
  typedef typename DataRecvMap::iterator        DataRecvMapIter;
  typedef typename DataRecvMap::const_iterator  DataRecvMapCIter;

  typedef typename std::vector<KT>          KeyList;
  typedef typename KeyList::iterator        KeyListIter;
  typedef typename KeyList::const_iterator  KeyListCIter;

  typedef typename std::vector<int>   ProcList;
  typedef typename ProcList::iterator ProcListIter;

  typedef typename std::pair<int,KT> ProcKeyPair;
  typedef typename std::vector<ProcKeyPair> ProcKeyList;
  typedef typename ProcKeyList::iterator ProcKeyListIter;

  typedef typename AC::iterator       ContainerIter;
  typedef typename AC::const_iterator ContainerCIter;

  // Constructors
  Directory( MG migrate,
	     DH distHash )
  : migrate_(migrate),
    distHash_(distHash)
  {}

  // Destructor
  ~Directory() {}

private:
  // No public copy construction, assignment, or equality operators
  Directory( const Directory & );

  Directory & operator=( const Directory & );

  bool operator==( const Directory & ) const;
  bool operator!=( const Directory & ) const;

public:

  // Add objects from directory.
  void addEntries( DataMap const & entries );

  // Remove objects from directory.
  void deleteEntries( KeyList & keys );

  // Get the items in the directory.
  void getEntries( KeyList & keys,
                   DataMap & entries );

  AC & container() { return container_; }
  ContainerIter & begin() { return container_.begin(); }
  ContainerIter & end() { return container_.end(); }

protected:

#ifdef Xyce_PARALLEL_MPI
  void pushKeys_( KeyList &, KeyList &, ProcList & );
  void pushData_( DataMap const &, DataRecvMap &, ProcList & );
#endif

  MG migrate_;
  DH distHash_;
  AC container_;

};

//-----------------------------------------------------------------------------
// Function      : Xyce::Parallel::Hash
// Purpose       : Basic hash class with impl. for std::string
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/20/03
//-----------------------------------------------------------------------------
template <typename T>
class Hash
{
 public:
  int operator()( const T & in ) { assert(0); return 0; }
  int size();
};

template <>
class Hash<std::string>
{
  int size_;

 public:

  Hash( int size )
  : size_( size )
  {}

  int operator()( const std::string & in )
  {
    int slen = in.length();
    int sum = 0;
    for( int i = 0; i < slen; ++i )
      sum += static_cast<int>( in[i] ); 

    return static_cast<int>( fmod( static_cast<double>( sum ), static_cast<double>(size_) ) );
  }

  int size() {return size_;}
};

template <>
class Hash<NodeID>
{
  int size_;

 public:

  Hash( int size )
  : size_( size )
  {}

  // Perform a hash on a pair<string,int>
  int operator()( const NodeID & in )
  {
    int slen = in.first.length();
    int sum = 0;
    for( int i = 0; i < slen; ++i )
      sum += static_cast<int>( (in.first)[i] ); 
    sum += in.second;

    return static_cast<int>( fmod( static_cast<double>( sum ), static_cast<double>(size_) ) );
  }

  int size() {return size_;}
};


//-----------------------------------------------------------------------------
// Function      : Xyce::Parallel::Directory::addEntries
// Purpose       : Add entries to directory
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/7/03
//-----------------------------------------------------------------------------
template <typename KT, typename DT, class DH, class AC, class MG>
void
Directory<KT,DT,DH,AC,MG>::
addEntries( DataMap const & entries )
{
#ifdef Xyce_PARALLEL_MPI

  DataRecvMap newEntries;
  ProcList procs;
  pushData_( entries, newEntries, procs );

  DataRecvMapCIter citDM = newEntries.begin();
  DataRecvMapCIter cendDM = newEntries.end();

#else

  DataMapCIter citDM  = entries.begin();
  DataMapCIter cendDM = entries.end();

#endif

  for( ; citDM != cendDM; ++citDM )
      container_.insert( *citDM );
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Parallel::Directory::deleteEntries
// Purpose       : Delete entries from directory
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/7/03
//-----------------------------------------------------------------------------
template <typename KT, typename DT, class DH, class AC, class MG>
void
Directory<KT,DT,DH,AC,MG>::
deleteEntries( KeyList & keys )
{
#ifdef Xyce_PARALLEL_MPI

  KeyList newKeys;
  ProcList procs;
  pushKeys_( keys, newKeys, procs );

  KeyListCIter citKL = newKeys.begin();
  KeyListCIter cendKL = newKeys.end();

#else

  KeyListCIter citKL  = keys.begin();
  KeyListCIter cendKL = keys.end();

#endif

  for( ; citKL != cendKL; ++citKL )
    container_.erase( *citKL );
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Parallel::Directory::getEntries
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/7/03
//-----------------------------------------------------------------------------
template <typename KT, typename DT, class DH, class AC, class MG>
void
Directory<KT,DT,DH,AC,MG>::
getEntries( KeyList & keys,
            DataMap & entries )
{
#ifdef Xyce_PARALLEL_MPI

  //Push Keys to owning processors
  KeyList newKeys;
  ProcList procs;
  pushKeys_( keys, newKeys, procs );
  //int numProcs = procs.size();

  KeyListCIter citKL  = newKeys.begin();
  KeyListCIter cendKL = newKeys.end();

  //Rvs migrate to move data from directory back to requesting procs
  DataMap newEntries;
  for( ; citKL != cendKL; ++citKL )
  {
    if( !container_.count( *citKL ) )
    {
      std::ostringstream os;
      os << "Data not in directory: " << *citKL << std::endl;
      N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL, os.str() );
    }
    newEntries[*citKL] = (container_.lower_bound( *citKL ))->second;
  }

  migrate_.rvs( procs, newKeys, newEntries, entries );

#else

  KeyListCIter citKL  = keys.begin();
  KeyListCIter cendKL = keys.end();
  for( ; citKL != cendKL; ++citKL )
  {
    if( !container_.count( *citKL ) )
    {
      std::ostringstream os;
      os << "Data not in directory: " << *citKL << std::endl;
      N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL, os.str() );
    }
    entries[*citKL] = (container_.lower_bound( *citKL ))->second;
  }

#endif
}

#ifdef Xyce_PARALLEL_MPI

//-----------------------------------------------------------------------------
// Function      : Xyce::Parallel::Directory::pushKeys_
// Purpose       : Pushes keys to owning processors based on hash function
// Special Notes : Forces keys to be ordered by proc number
// Scope         : Protected
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/7/03
//-----------------------------------------------------------------------------
template <typename KT, typename DT, class DH, class AC, class MG>
void
Directory<KT,DT,DH,AC,MG>::
pushKeys_( KeyList & sKeys,
           KeyList & rKeys,
	   ProcList & procs )
{

  KeyListCIter itKL  = sKeys.begin();
  KeyListCIter endKL = sKeys.end();

  procs.clear();
  for( ; itKL != endKL; ++itKL )
    procs.push_back( distHash_(*itKL) );

  if( !IsSorted( procs ) ) SortContainer2( procs, sKeys );
  
  migrate_( procs, sKeys, rKeys );
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Parallel::Directory::pushData_
// Purpose       : Pushes data to owning processors based on hash function
// Special Notes :
// Scope         : Protected
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/7/03
//-----------------------------------------------------------------------------
template <typename KT, typename DT, class DH, class AC, class MG>
void
Directory<KT,DT,DH,AC,MG>::
pushData_( DataMap const & sData,
           DataRecvMap & rData,
	   ProcList & procs )
{
  DataMapCIter itDM  = sData.begin();
  DataMapCIter endDM = sData.end();

  procs.clear();
  for( ; itDM != endDM; ++itDM )
    procs.push_back( distHash_(itDM->first) );

  migrate_( procs, sData, rData );
}

#endif

} //namespace Parallel
} //namespace Xyce

#endif

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
// Filename       : $RCSfile: N_UTL_Functors.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Robert J. Hoekstra
//
// Creation Date  : 7/16/01
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.11 $
//
// Revision Date  : $Date: 2014/02/24 23:49:28 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef _N_UTL_Functors_h
#define _N_UTL_Functors_h 1

// ---------- Standard Includes ----------

#include <functional>
#include <map>

// ----------   Xyce Includes   ----------

#include <N_UTL_Xyce.h>

//-----------------------------------------------------------------------------
// Class         : DeletePtr
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 7/16/01
//-----------------------------------------------------------------------------
template < typename T >
struct DeletePtr : public std::unary_function < const T *, void >
{
  void operator() (const T * ptr) const { delete ptr; }
};

//-----------------------------------------------------------------------------
// Class         : FirstOfPair
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 7/16/01
//-----------------------------------------------------------------------------
template < typename T, typename U >
struct FirstOfPair : public std::unary_function < const T &, const U & >
{
  const U & operator() (const T & ref) const { return ref.first; }
};

//-----------------------------------------------------------------------------
// Class         : SortContainer2
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/23/04
//-----------------------------------------------------------------------------
template < typename T, typename U >
void SortContainer2( T & firstContainer, U & secondContainer )
{
  typedef typename std::multimap< typename T::value_type, typename U::value_type> UTMultiMap;

  UTMultiMap SortMap;

  typename T::iterator iterT = firstContainer.begin();
  typename T::iterator endT = firstContainer.end();
  typename U::iterator iterU = secondContainer.begin();
  typename U::iterator endU = secondContainer.end();

  for( ; (iterT!=endT)||(iterU!=endU) ; ++iterT, ++iterU )
    SortMap.insert( typename UTMultiMap::value_type( *iterT, *iterU ) );

  firstContainer.clear();
  secondContainer.clear();

  typename UTMultiMap::iterator iterUTM = SortMap.begin();
  typename UTMultiMap::iterator endUTM = SortMap.end();

  for( ; iterUTM != endUTM; ++iterUTM )
  {
    firstContainer.push_back( iterUTM->first );
    secondContainer.push_back( iterUTM->second );
  }
}

//-----------------------------------------------------------------------------
// Class         : IsSorted
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 8/23/04
//-----------------------------------------------------------------------------
template < typename T >
bool IsSorted( T & container )
{
  if( container.size() < 2 ) return true;

  typename T::iterator iterT = container.begin();
  typename T::iterator endT = container.end();
  typename T::iterator iterTPlus = iterT;
  iterTPlus++;

  for( ; iterTPlus != endT; ++iterT, ++iterTPlus )
    if( !(*iterT<*iterTPlus) ) return false;

  return true;
}

//-----------------------------------------------------------------------------
// Class         : N_UTL_Less 
// Purpose       : Class implemented for use with RedStorm PGI compilers
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Org. 1437
// Creation Date : 3/27/06
//-----------------------------------------------------------------------------
template < typename T, typename ST >
class LessThan : public std::less<ST>
{
   public:
     bool operator() ( const T& val1, const ST& val2 ) const { return val1 < val2; }
};

#endif


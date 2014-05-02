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
// Filename       : $RCSfile: N_UTL_Graph.h,v $
//
// Purpose        : Simple undirected graph class
//
// Special Notes  : 
//
// Creator        : Robert J. Hoekstra, SNL, Electrical & MicroSystems
//
// Creation Date  : 8/10/06
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.23 $
//
// Revision Date  : $Date: 2014/02/24 23:49:28 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_UTL_Graph_H
#define Xyce_UTL_Graph_H

// ---------- Standard Includes ----------

#include <vector>
#include <map>
#include <queue>
#include <set>
#include <utility>
#include <algorithm>

#include <iostream>

// ----------   Xyce Includes   ----------

// ----------  Other Includes   ----------

//-----------------------------------------------------------------------------
// Class         : N_UTL_Graph_Comparator
// Purpose       : Simple class to provide a comparison for removing keys from the 
//                 N_UTL_Graph class
// Special Notes :
// Scope         : Public
// Creator       : Heidi K. Thornquist, SNL
//-----------------------------------------------------------------------------
template <typename Index>
class N_UTL_Graph_Comparator
{
  public:

  // Constructor 
  N_UTL_Graph_Comparator( const std::vector<Index> & inputIDs )
  : numIDs_(inputIDs.size()), inputIDs_(inputIDs)
  {}

  // Binary search
  int binary_search( const std::vector<Index> & array, Index id, int low, int high )
  {
    if (high < low)
      return -1;
    int mid = (low+high) / 2;
    if (array[mid] > id)
      return binary_search( array, id, low, mid-1 );
    else if (array[mid] < id)
      return binary_search( array, id, mid+1, high );
    else 
      return mid;
  }

  // Comparison operator/functor
  bool operator()( Index id )
  {
    bool ret = false;
    int idx = binary_search( inputIDs_, id, 0, numIDs_-1 );
    if (idx > -1)
      ret = true;

    return ret;
  }
 
  private:

  int numIDs_; 
  std::vector<Index> inputIDs_;

};

//-----------------------------------------------------------------------------
// Class         : N_UTL_Graph
// Purpose       : Simple undirected graph
// Special Notes :
// Scope         : Public
// Creator       : Robert J. Hoekstra, SNL
//-----------------------------------------------------------------------------
template <typename Key1Type, typename Key2Type, typename DataType>
class N_UTL_Graph
{
 public:

  typedef int Index;
  typedef typename std::map<Index,Key1Type> Key1Map;
  typedef typename std::map<Key1Type,Index> Index1Map;
  typedef typename std::map<Key1Type,DataType> Data1Map;
  typedef typename std::map<Index,Key2Type> Key2Map;
  typedef typename std::map<Key2Type,Index> Index2Map;
  typedef typename std::map<Key2Type,DataType> Data2Map;
  typedef typename std::map<Key1Type,Key2Type> Key12Map;
  typedef typename std::map<Key2Type,Key1Type> Key21Map;

  // Constructors
  N_UTL_Graph() : numRemovedNodes_(0) {}

  // Destructor
  virtual ~N_UTL_Graph() {}

 private:

  // Copy constructor (private).
  // Not supported for now
  N_UTL_Graph(const N_UTL_Graph&);

  // Assignment operator (private).
  // Not supported for now
  N_UTL_Graph& operator=(const N_UTL_Graph&);

  // Equality Operations (private).
  // Not supported for now
  int operator==(const N_UTL_Graph& right) const;
  int operator!=(const N_UTL_Graph& right) const;

 public:

  bool insertNode(const Key1Type& key1, const Key2Type& key2, const std::vector<Key1Type> adj, DataType& data);
  void chgKey2(const Key1Type& key1, const Key2Type& key2);

  //deleteNode(KeyType&);

  int numNodes() const;
  bool checkKey(const Key1Type& key) const;
  bool checkKey(const Key2Type& key) const;

  Key1Type getKey1(const Key2Type& key) const;
  Key2Type getKey2(const Key1Type& key) const;

  DataType& getData(const Key1Type& key);
  DataType& getData(const Key2Type& key);

  std::vector<DataType>& getData(const std::vector<Key1Type>& keys);
  std::vector<DataType>& getData(const std::vector<Key2Type>& keys);

  int numAdjacent(const Key1Type& key);
  int numAdjacent(const Key2Type& key);

  std::vector<Key1Type> getAdjacent(const Key1Type& key);
  std::vector<Key2Type> getAdjacent(const Key2Type& key);

  Index getIndex(const Key1Type & key) { return rvsKeys1_[ key ]; }
  std::vector<Index> & getAdjacentRow(const Index idx) { return adjacencyGraph_[idx]; }

  void addToAdjacent(const Key1Type& oldkey, const Key1Type& key, std::vector<Key1Type> & newAdjVec);
  void addToAdjacent(const Key2Type& oldkey, const Key2Type& key, std::vector<Key2Type> & newAdjVec);

  void replaceAdjacent(const Key1Type & oldKey1, const Key1Type & newKey1);
  void replaceAdjacent(const Key2Type & oldKey2, const Key2Type & newKey2);

  void removeKey(const Key1Type oldKey1);
  void removeKey(const Key2Type oldKey2);

  void removeKeys(const std::vector<Key1Type> & oldKeys1);
  void removeKeys(const std::vector<Key2Type> & oldKeys2);

  std::vector<Key1Type> getSingletons();

  int checkGraphState();

  // Give access to the data map in the event that the nodes need 
  // to be traversed without an ordered list (BFT) being generated.  
  const std::map<Key1Type,DataType>& getData1Map() { return data1_; }

  const std::vector<Key1Type>& getBFT();
  int generateBFT();
  int generateBFT(const Key1Type& key);
  Key1Type getCenter(double threshold, int maxTries);

  void print(std::ostream & ostr) const;

 private:

  int generateBFT_(const Index& start);

  std::vector< std::vector<Index> > adjacencyGraph_;

  // if we remove any nodes from the adjacencyGraph_, then we'll have empty
  // rows.  We could erase these, but then we would have to reindx the keymaps.
  // For now we'll just count the number of removals
  int numRemovedNodes_;

  Key1Map keys1_;
  Index1Map rvsKeys1_;
  Data1Map data1_;
  
  mutable Key12Map keys12_;
  mutable Key21Map keys21_;

  std::vector<Index> bft_;
  std::vector<Key1Type> bftKeys_;
};

//-----------------------------------------------------------------------------
// Function      : N_UTL_Graph::insertNode
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Robert J. Hoekstra, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
template <typename Key1Type, typename Key2Type, typename DataType>
inline bool N_UTL_Graph<Key1Type,Key2Type,DataType>::insertNode(
    const Key1Type& key1, const Key2Type& key2, 
    const std::vector<Key1Type> adj, DataType& data)
{
  if(checkKey(key1)) return false;
  std::vector<Index> ids;

  size_t size = adj.size();
  for(size_t i = 0; i < size; ++i)
    ids.push_back(rvsKeys1_[adj[i]]);
  adjacencyGraph_.push_back(ids);

  //add back edges
  int loc = adjacencyGraph_.size()-1;
  for(size_t i = 0; i < size; ++i)
    adjacencyGraph_[ rvsKeys1_[adj[i]] ].insert(adjacencyGraph_[ rvsKeys1_[adj[i]] ].begin(), loc);

  Index currIndex = adjacencyGraph_.size()-1;
  keys1_[currIndex] = key1;
  rvsKeys1_[key1] = currIndex;
  data1_[key1] = data;
  keys12_[key1] = key2;
  keys21_[key2] = key1;

  return true;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Graph::chgKey2
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Robert J. Hoekstra, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
template <typename Key1Type, typename Key2Type, typename DataType>
inline void N_UTL_Graph<Key1Type,Key2Type,DataType>::chgKey2
  (const Key1Type& key1, const Key2Type& key2)
{
  keys12_[key1] = key2;
  keys21_[key2] = key1; 
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Graph::numNodes
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Robert J. Hoekstra, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
template <typename Key1Type, typename Key2Type, typename DataType>
inline int N_UTL_Graph<Key1Type,Key2Type,DataType>::numNodes() const
{
  return adjacencyGraph_.size() - numRemovedNodes_;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Graph::checkKey
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Robert J. Hoekstra, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
template <typename Key1Type, typename Key2Type, typename DataType>
inline bool N_UTL_Graph<Key1Type,Key2Type,DataType>::checkKey
  (const Key1Type& key) const
{
  if(rvsKeys1_.count(key)) return true;
  return false;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Graph::checkKey
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Robert J. Hoekstra, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
template <typename Key1Type, typename Key2Type, typename DataType>
inline bool N_UTL_Graph<Key1Type,Key2Type,DataType>::checkKey
  (const Key2Type& key) const
{
  if(keys21_.count(key)) return true;
  return false;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Graph::getKey1
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Robert J. Hoekstra, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
template <typename Key1Type, typename Key2Type, typename DataType>
inline Key1Type N_UTL_Graph<Key1Type,Key2Type,DataType>::getKey1
  (const Key2Type& key) const
{
  return keys21_[key];
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Graph::getKey2
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Robert J. Hoekstra, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
template <typename Key1Type, typename Key2Type, typename DataType>
inline Key2Type N_UTL_Graph<Key1Type,Key2Type,DataType>::getKey2
  (const Key1Type& key) const
{
  return keys12_[key];
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Graph::getData
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Robert J. Hoekstra, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
template <typename Key1Type, typename Key2Type, typename DataType>
inline DataType& N_UTL_Graph<Key1Type,Key2Type,DataType>::getData
  (const Key1Type& key)
{
  return data1_[key];
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Graph::getData
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Robert J. Hoekstra, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
template <typename Key1Type, typename Key2Type, typename DataType>
inline DataType& N_UTL_Graph<Key1Type,Key2Type,DataType>::getData
  (const Key2Type& key)
{
  return getData(keys21_[key]);
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Graph::getData
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Robert J. Hoekstra, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
template <typename Key1Type, typename Key2Type, typename DataType>
inline std::vector<DataType>& N_UTL_Graph<Key1Type,Key2Type,DataType>::getData
  (const std::vector<Key1Type>& keys)
{
  std::vector<DataType> data;
  for(int i = 0; i < keys.size(); ++i)
    data.push_back(data1_[keys[i]]);

  return data;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Graph::getData
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Robert J. Hoekstra, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
template <typename Key1Type, typename Key2Type, typename DataType>
inline std::vector<DataType>& N_UTL_Graph<Key1Type,Key2Type,DataType>::getData
  (const std::vector<Key2Type>& keys)
{
  std::vector<Key1Type> keys1;
  for(int i = 0; i < keys.size(); ++i)
    keys1.push_back(keys21_[keys[i]]);

  return getData(keys1);
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Graph::numAdjacent
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Robert J. Hoekstra, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
template <typename Key1Type, typename Key2Type, typename DataType>
inline int N_UTL_Graph<Key1Type,Key2Type,DataType>::numAdjacent
  (const Key1Type& key)
{
  Index id = rvsKeys1_[key];
  return adjacencyGraph_[id].size();
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Graph::numAdjacent
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Robert J. Hoekstra, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
template <typename Key1Type, typename Key2Type, typename DataType>
inline int N_UTL_Graph<Key1Type,Key2Type,DataType>::numAdjacent
  (const Key2Type& key)
{
  return numAdjacent(keys21_[key]);
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Graph::getAdjacent
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Robert J. Hoekstra, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
template <typename Key1Type, typename Key2Type, typename DataType>
inline std::vector<Key1Type> N_UTL_Graph<Key1Type,Key2Type,DataType>::getAdjacent
  (const Key1Type& key)
{
  std::vector<Key1Type> adjKeys;
  Index id = rvsKeys1_[key];
  int size = adjacencyGraph_[id].size();
  for(int i = 0; i < size; ++i)
    adjKeys.push_back(keys1_[ adjacencyGraph_[id][i] ]);
  return adjKeys;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Graph::getAdjacent
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Robert J. Hoekstra, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
template <typename Key1Type, typename Key2Type, typename DataType>
inline std::vector<Key2Type> N_UTL_Graph<Key1Type,Key2Type,DataType>::getAdjacent
  (const Key2Type& key)
{
  std::vector<Key1Type> keys1 = getAdjacent(keys21_[key]);

  std::vector<Key2Type> keys2;
  for(unsigned int i = 0; i < keys1.size(); ++i)
    keys2.push_back(keys12_[keys1[i]]);
  return keys2;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Graph::addToAdjacent
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
template <typename Key1Type, typename Key2Type, typename DataType>
inline void N_UTL_Graph<Key1Type,Key2Type,DataType>::addToAdjacent
  (const Key1Type& oldkey, const Key1Type& key, std::vector<Key1Type> & newAdjVec)
{
  // originally I had this as setAdjacent() where it would clear and then 
  // write a new adjacency info to adjacencyGraph_.
  // A flaw with using this as setAdjacent is that it can change the order of 
  // edge nodes.  potentially reversing an element on another part of the circuit.
  // so rather than erasing and then setting a new adjacent list, we'll just add to what
  // is there
  
  // add in new values
  Index id = rvsKeys1_[key];
  Index oldId = rvsKeys1_[oldkey];
  int extraElementsSize = newAdjVec.size();
  for(int i=0; i< extraElementsSize; i++)
  {
    adjacencyGraph_[id].push_back(rvsKeys1_[newAdjVec[i]]);
    // and any new edges
    // need to do this by changing the edge's oldId to the id of the new key
    Index edgeIndex = rvsKeys1_[newAdjVec[i]];
    std::vector< Index >::iterator endLoc = adjacencyGraph_[ edgeIndex ].end();
    std::vector< Index >::iterator edgeLoc = find(adjacencyGraph_[ edgeIndex ].begin(), endLoc, oldId);
    if(edgeLoc == endLoc)
    {
      adjacencyGraph_[edgeIndex].push_back(id);
    }
    else
    {
      *edgeLoc = id;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Graph::addToAdjacent
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
template <typename Key1Type, typename Key2Type, typename DataType>
inline void N_UTL_Graph<Key1Type,Key2Type,DataType>::addToAdjacent
  (const Key2Type& oldkey, const Key2Type& key, std::vector<Key2Type> & newAdjVec)
{
  // originally I had this as setAdjacent() where it would clear and then 
  // write a new adjacency info to adjacencyGraph_.
  // A  flaw with using this as setAdjacent is that it can change the order of 
  // edge nodes.  potentially reversing an element on another part of the circuit.
  // so rather than erasing and then setting a new adjacent list, we'll just add to what
  // is there
  
  // add in new values
  Index id = rvsKeys1_[ keys21_[key]]; 
  Index oldId = rvsKeys1_[ keys21_[oldkey]];
  int extraElementsSize = newAdjVec.size();
  for(int i=0; i< extraElementsSize; i++)
  {
    adjacencyGraph_[id].push_back(rvsKeys1_[newAdjVec[i]]);
    // and any new edges
    // need to do this by changing the edge's oldId to the id of the new key
    Index edgeIndex = rvsKeys1_[newAdjVec[i]];
    std::vector< Index >::iterator endLoc = adjacencyGraph_[ edgeIndex ].end();
    std::vector< Index >::iterator edgeLoc = find(adjacencyGraph_[ edgeIndex ].begin(), endLoc, oldId);
    if(edgeLoc == endLoc)
    {
      adjacencyGraph_[edgeIndex].push_back(id);
    }
    else
    {
      *edgeLoc = id;
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Graph::replaceAdjacent
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
template <typename Key1Type, typename Key2Type, typename DataType>
inline void N_UTL_Graph<Key1Type,Key2Type,DataType>::replaceAdjacent
  (const Key1Type & oldKey1, const Key1Type & newKey1)
{
  // get id's for these keys
  Index oldId = rvsKeys1_[ oldKey1 ];
  Index newId = rvsKeys1_[ newKey1 ];
  
  // loop through adjancy graph and replace any references to oldKey with
  // newKey only if newKey doesn't already exist on that row
  int numAdjRows = adjacencyGraph_.size();
  for(int i=0; i<numAdjRows; ++i)
  {
    std::vector< Index >::iterator beginLoc = adjacencyGraph_[i].begin();
    std::vector< Index >::iterator endLoc = adjacencyGraph_[i].end();
    // look for the old key
    std::vector< Index >::iterator oldKeyLoc = find(beginLoc, endLoc, oldId);
    if(oldKeyLoc != endLoc)
    {
      // found old key, so overwrite it
      *oldKeyLoc = newId;
      // 
      // original logic removed the old one and added a new one
      // this could break ordering of the id's which would be bad
      // // found old key, so erase it
      // adjacencyGraph_[i].erase(oldKeyLoc);
      // // also check for new key before inserting it
      // std::vector< Index >::iterator newKeyLoc = find(beginLoc, endLoc, newId);
      // if(newKeyLoc == endLoc)
      // {
      //   // new key isn't there so add it in
      //   adjacencyGraph_[i].push_back(newId);
      // }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Graph::replaceAdjacent
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
template <typename Key1Type, typename Key2Type, typename DataType>
inline void N_UTL_Graph<Key1Type,Key2Type,DataType>::replaceAdjacent
  (const Key2Type & oldKey2, const Key2Type & newKey2)
{
  // get id's for these keys
  Index oldId = rvsKeys1_[ keys21_[oldKey2] ];
  Index newId = rvsKeys1_[ keys21_[newKey2] ];
  
  // loop through adjancy graph and replace any references to oldKey with
  // newKey only if newKey doesn't already exist on that row
  int numAdjRows = adjacencyGraph_.size();
  for(int i=0; i<numAdjRows; ++i)
  {
    std::vector< Index >::iterator beginLoc = adjacencyGraph_[i].begin();
    std::vector< Index >::iterator endLoc = adjacencyGraph_[i].end();
    // look for the old key
    std::vector< Index >::iterator oldKeyLoc = find(beginLoc, endLoc, oldId);
    if(oldKeyLoc != endLoc)
    {
      // found old key, so overwrite it
      *oldKeyLoc = newId;
      // 
      // original logic removed the old one and added a new one
      // this could break ordering of the id's which would be bad
      // // found old key, so erase it
      // adjacencyGraph_[i].erase(oldKeyLoc);
      // // also check for new key before inserting it
      // std::vector< Index >::iterator newKeyLoc = find(beginLoc, endLoc, newId);
      // if(newKeyLoc == endLoc)
      // {
      //   // new key isn't there so add it in
      //   adjacencyGraph_[i].push_back(newId);
      // }
    }
  }
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Graph::removeKey
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
template <typename Key1Type, typename Key2Type, typename DataType>
inline void N_UTL_Graph<Key1Type,Key2Type,DataType>::removeKey
  (const Key1Type oldKey1)
{
  // get keys and id
  Index id = rvsKeys1_[ oldKey1 ];
  Key2Type oldKey2 = keys12_[ oldKey1 ];
  
  // clear the row it uses
  // this leaves an empty row in place, but that's ok
  adjacencyGraph_[id].clear();
  
  // erase row that this key points to in adjacency graph
  // std::vector< std::vector< Index > >::iterator oldRowItr = adjacencyGraph_.begin() + id;
  // adjacencyGraph_.erase(oldRowItr);
  // if I erase the row it throws off the key maps.  So I need to fix up those 
  
  // search through adjacencyGraph_ and remove all references to this key
  // then remove it from the maps as well
  int numAdjRows = adjacencyGraph_.size();
  for(int i=0; i<numAdjRows; ++i)
  {
    if(!adjacencyGraph_[i].empty())
    {
      std::vector< Index >::iterator beginLoc = adjacencyGraph_[i].begin();
      std::vector< Index >::iterator endLoc = adjacencyGraph_[i].end();
      // look for the old key
      std::vector< Index >::iterator oldKeyLoc = find(beginLoc, endLoc, id);
      if(oldKeyLoc != endLoc)
      {
        // found old key, so erase it
        adjacencyGraph_[i].erase(oldKeyLoc);
        i--; // ugly.  we need to ensure we find all of the oldKeyLoc values.
      }
    }
  }
  
  // clean up maps
  keys1_.erase(id);         // Key1Map keys1_;
  rvsKeys1_.erase(oldKey1); // Index1Map rvsKeys1_;
  data1_.erase(oldKey1);    // Data1Map data1_;
  keys12_.erase(oldKey1);   // Key12Map keys12_;
  keys21_.erase(oldKey2);    // Key21Map keys21_;
  
  numRemovedNodes_++;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Graph::removeKey
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
template <typename Key1Type, typename Key2Type, typename DataType>
inline void N_UTL_Graph<Key1Type,Key2Type,DataType>::removeKey
  (const Key2Type oldKey2)
{
  // get keys and id
  Key1Type oldKey1 = keys21_[ oldKey2 ];
  Index id = rvsKeys1_[oldKey1];
  
  // clear the row it uses
  adjacencyGraph_[id].clear();
  
  // search through adjacencyGraph_ and remove all references to this key
  // then remove it from the maps as well
  int numAdjRows = adjacencyGraph_.size();
  for(int i=0; i<numAdjRows; ++i)
  {
    if(!adjacencyGraph_[i].empty())
    {
      std::vector< Index >::iterator beginLoc = adjacencyGraph_[i].begin();
      std::vector< Index >::iterator endLoc = adjacencyGraph_[i].end();
      // look for the old key
      std::vector< Index >::iterator oldKeyLoc = find(beginLoc, endLoc, id);
      if(oldKeyLoc != endLoc)
      {
        // found old key, so erase it
        adjacencyGraph_[i].erase(oldKeyLoc);
        i--;
      }
    }
  }
  
  // clean up maps
  keys1_.erase(id);         // Key1Map keys1_;
  rvsKeys1_.erase(oldKey1); // Index1Map rvsKeys1_;
  data1_.erase(oldKey1);    // Data1Map data1_;
  keys12_.erase(oldKey1);   // Key12Map keys12_;
  keys21_.erase(oldKey2);    // Key21Map keys21_;
  
  numRemovedNodes_++;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Graph::removeKeys
// Purpose       : 
// Special Notes : This is a collective removal from the adjacencyGraph for efficiency.
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
template <typename Key1Type, typename Key2Type, typename DataType>
inline void N_UTL_Graph<Key1Type,Key2Type,DataType>::removeKeys
  (const std::vector<Key1Type> & oldKeys1)
{
  int numKeysToRemove = oldKeys1.size();
  std::vector<Index> ids(numKeysToRemove);

  // Go through all the keys, retrieve their ids and clean up the maps
  for(int i=0; i<numKeysToRemove; ++i)
  {
    // get keys and id
    ids[i] = rvsKeys1_[ oldKeys1[i] ];
    Key2Type oldKey2 = keys12_[ oldKeys1[i] ];
 
    // clean up maps
    keys1_.erase(ids[i]);         // Key1Map keys1_;
    rvsKeys1_.erase(oldKeys1[i]); // Index1Map rvsKeys1_;
    data1_.erase(oldKeys1[i]);    // Data1Map data1_;
    keys12_.erase(oldKeys1[i]);   // Key12Map keys12_;
    keys21_.erase(oldKey2);    // Key21Map keys21_;
  
    // clear the row it uses
    // this leaves an empty row in place, but that's ok
    adjacencyGraph_[ids[i]].clear();
  }

  // Sort the list of ids to expedite searching
  std::sort( ids.begin(), ids.end() );

  // search through adjacencyGraph_ and remove all references to this key
  int numAdjRows = adjacencyGraph_.size();
  for(int i=0; i<numAdjRows; ++i)
  {
    if(!adjacencyGraph_[i].empty())
    {
      std::vector< Index >::iterator beginLoc = adjacencyGraph_[i].begin();
      std::vector< Index >::iterator endLoc = adjacencyGraph_[i].end();
      adjacencyGraph_[i].erase( remove_if(beginLoc, endLoc, N_UTL_Graph_Comparator<Index>(ids)), endLoc);
    }
  }

  numRemovedNodes_ += numKeysToRemove;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Graph::removeKeys
// Purpose       : 
// Special Notes : This is a collective removal from the adjacencyGraph for efficiency.
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
template <typename Key1Type, typename Key2Type, typename DataType>
inline void N_UTL_Graph<Key1Type,Key2Type,DataType>::removeKeys
  (const std::vector<Key2Type> & oldKeys2)
{
  int numKeysToRemove = oldKeys2.size();
  std::vector<Index> ids(numKeysToRemove);

  // Go through all the keys, retrieve their ids and clean up the maps
  for(int i=0; i<numKeysToRemove; ++i)
  {
    // get keys and id
    Key1Type oldKey1 = keys21_[ oldKeys2[i] ];
    ids[i] = rvsKeys1_[ oldKey1 ];

    // clean up maps
    keys1_.erase(ids[i]);       // Key1Map keys1_;
    rvsKeys1_.erase(oldKey1);   // Index1Map rvsKeys1_;
    data1_.erase(oldKey1);      // Data1Map data1_;
    keys12_.erase(oldKey1);     // Key12Map keys12_;
    keys21_.erase(oldKeys2[i]); // Key21Map keys21_;
  
    // clear the row it uses
    // this leaves an empty row in place, but that's ok
    adjacencyGraph_[ids[i]].clear();
  }

  // Sort the list of ids to expedite searching
  std::sort( ids.begin(), ids.end() );

  // search through adjacencyGraph_ and remove all references to this key
  int numAdjRows = adjacencyGraph_.size();
  for(int i=0; i<numAdjRows; ++i)
  {
    if(!adjacencyGraph_[i].empty())
    {
      std::vector< Index >::iterator beginLoc = adjacencyGraph_[i].begin();
      std::vector< Index >::iterator endLoc = adjacencyGraph_[i].end();
      adjacencyGraph_[i].erase( remove_if(beginLoc, endLoc, N_UTL_Graph_Comparator<Index>(ids)), endLoc);
    }
  }
  
  numRemovedNodes_ += numKeysToRemove;
}


//-----------------------------------------------------------------------------
// Function      : N_UTL_Graph::getSingletons
// Purpose       : 
// Special Notes : This method returns the viable nodes that are graph singletons.
// Scope         : public
// Creator       : Heidi Thornquist, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
template <typename Key1Type, typename Key2Type, typename DataType>
inline std::vector<Key1Type> N_UTL_Graph<Key1Type,Key2Type,DataType>::getSingletons()
{
  std::vector<Key1Type> singletonKeys;

  int numAdjRows = adjacencyGraph_.size();
  for (int i=0; i<numAdjRows; ++i)
  {
    if (adjacencyGraph_[i].empty() && keys1_.count(i)>0)
      singletonKeys.push_back(keys1_[i]);
  }
 
  return singletonKeys;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Graph::checkGraphState
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
template <typename Key1Type, typename Key2Type, typename DataType>
inline int N_UTL_Graph<Key1Type,Key2Type,DataType>::checkGraphState()
{
  // this is for debugging
  // it checks that the maps are consisitent (i.e. they all have the same
  // keys and indices).
  // then it traverses the adjacencyGraph_ and makes sure that it only 
  // points to items that are in the key maps
  
  int numAdjRows = adjacencyGraph_.size();
  for(int i=0; i<numAdjRows; ++i)
  {
    int numNodes = adjacencyGraph_[i].size();
    for(int j=0; j< numNodes; j++)
    {
      // check index in all maps
      Index testIndex = adjacencyGraph_[i][j];
      Key1Type key1val = keys1_[ testIndex ];
      if(keys1_.count(testIndex) == 0)
      {
        return 1;
      }
      if(rvsKeys1_.count(key1val) == 0)
      {
        return 2;
      }
      if(data1_.count(key1val) == 0)
      {
        return 3;
      }
      if(keys12_.count(key1val) == 0)
      {
        return 4;
      }
    }
  }
  return 0;
}
//
//#endif // supernode
 
//-----------------------------------------------------------------------------
// Function      : N_UTL_Graph::getBFT
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Robert J. Hoekstra, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
template <typename Key1Type, typename Key2Type, typename DataType>
inline const std::vector<Key1Type>& N_UTL_Graph<Key1Type,Key2Type,DataType>::getBFT()
{
  if(bft_.empty()) generateBFT();
  return bftKeys_;
}
    
//-----------------------------------------------------------------------------
// Function      : N_UTL_Graph::generateBFT
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Robert J. Hoekstra, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
template <typename Key1Type, typename Key2Type, typename DataType>
inline int N_UTL_Graph<Key1Type,Key2Type,DataType>::generateBFT()
{
#ifdef Xyce_GRAPH_DEBUG
  int state = checkGraphState();
  if(state != 0)
  {
    Xyce::dout() << "Graph is in inconsistent state!  checkGraphState returned: " << state << std::endl;
  }
#endif

  Index firstIndex = (*(keys1_.begin())).first;
  return generateBFT_(firstIndex);
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Graph::generateBFT
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Robert J. Hoekstra, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
template <typename Key1Type, typename Key2Type, typename DataType>
inline int N_UTL_Graph<Key1Type,Key2Type,DataType>::generateBFT
  (const Key1Type& key)
{
#ifdef Xyce_GRAPH_DEBUG
  int state = checkGraphState();
  if(state != 0)
  {
    Xyce::dout() << "Graph is in inconsistent state!  checkGraphState returned: " << state << std::endl;
  }
#endif

  Index id = rvsKeys1_[key];
  return generateBFT_(id);
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Graph::getCenter
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Robert J. Hoekstra, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
template <typename Key1Type, typename Key2Type, typename DataType>
inline Key1Type N_UTL_Graph<Key1Type,Key2Type,DataType>::getCenter(double threshold, int maxTries)
{
  int cutoff = int (numNodes()*threshold);
 
  int numTries = 1;
  for(typename Key1Map::const_iterator it_k1m = keys1_.begin(); it_k1m != keys1_.end(); ++it_k1m)
  {
    int lvl = generateBFT_((*it_k1m).first);
    if((lvl < cutoff) || (numTries == maxTries)) 
      return (*it_k1m).second;

    numTries++;
  }

  return (*(keys1_.begin())).second;
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Graph::print
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Robert J. Hoekstra, SNL
// Creation Date : 
//-----------------------------------------------------------------------------
template <typename Key1Type, typename Key2Type, typename DataType>
inline void N_UTL_Graph<Key1Type,Key2Type,DataType>::print(std::ostream & ostr) const
{
  ostr << "-------------------- Basic Graph ----------------------------\n";
  ostr << "Adjacency Graph\n";
  ostr << "---------------\n";
  for(size_t i = 0; i < adjacencyGraph_.size(); ++i)
  {
    ostr << "Node " << i << " : ";
    for(size_t j = 0; j < adjacencyGraph_[i].size(); ++j)
      ostr << " " << adjacencyGraph_[i][j];
    ostr << std::endl;
  }
  ostr << "---------------\n";
  ostr << "Key1Map\n";
  for(typename Key1Map::const_iterator it_k1m = keys1_.begin(); it_k1m != keys1_.end(); ++it_k1m)
    ostr << it_k1m->first << ":" << it_k1m->second << std::endl;
  ostr << "-------\n";
  ostr << "Index1Map\n";
  for(typename Index1Map::const_iterator it_i1m = rvsKeys1_.begin(); it_i1m != rvsKeys1_.end(); ++it_i1m)
    ostr << it_i1m->first << ":" << it_i1m->second << std::endl;
  ostr << "-------\n";
  ostr << "Data1Map\n";
  for(typename Data1Map::const_iterator it_d1m = data1_.begin(); it_d1m != data1_.end(); ++it_d1m)
    ostr << it_d1m->first << ":" << it_d1m->second << std::endl;
  ostr << "-------\n";
  ostr << "Key12Map\n";
  for(typename Key12Map::const_iterator it_k12m = keys12_.begin(); it_k12m != keys12_.end(); ++it_k12m)
    ostr << it_k12m->first << ":" << it_k12m->second << std::endl;
  ostr << "-------\n";
  ostr << "Key21Map\n";
  for(typename Key21Map::const_iterator it_k21m = keys21_.begin(); it_k21m != keys21_.end(); ++it_k21m)
    ostr << it_k21m->first << ":" << it_k21m->second << std::endl;
  ostr << "-------\n";
  ostr << "BFT\n";
  for(size_t i = 0; i < bft_.size(); ++i)
    ostr << bft_[i] << ":" << bftKeys_[i] << std::endl;
  ostr << "-------\n";
  ostr << "-------------------- Basic Graph END ------------------------\n";
}

//-----------------------------------------------------------------------------
// Function      : N_UTL_Graph::generateBFT_
// Purpose       : 
// Special Notes : 
// Scope         : private
// Creator       : Robert J. Hoekstra, SNL
//               : Heidi Thornquist, SNL (modified search for efficiency)
// Creation Date : 
//-----------------------------------------------------------------------------
template <typename Key1Type, typename Key2Type, typename DataType>
inline int N_UTL_Graph<Key1Type,Key2Type,DataType>::generateBFT_
  (const Index& start)
{
  bft_.clear();
  bftKeys_.clear();

  // Work queue
  std::queue< std::pair<Index, int> > idQueue;

  // Keep track of which IDs have been found (false=not found, true=found)
  std::vector<bool> foundIds(adjacencyGraph_.size(), false);  

  // Keep track of all the levels if the graph is a forest
  std::vector<int> levels;  
  
  Index localCopyOfStart = start;

  // Before we push back "start" need to verify that it's valid
  bool startIsValid = false;
  int numAdjacentRows = adjacencyGraph_.size();
  while(!startIsValid)
  {
    if(keys1_.count(localCopyOfStart) > 0)
    {
      startIsValid=true;
    }
    else
    {
      localCopyOfStart++;
      if(localCopyOfStart > (numAdjacentRows-1))
      {
        localCopyOfStart = 0;
      }
    }
  }
 
  // Initialize level, root, and work queue
  int level = 0;
  unsigned int root = 0;
  idQueue.push(std::make_pair(localCopyOfStart,level));
  bft_.push_back(localCopyOfStart);
  foundIds[localCopyOfStart] = true;

  while(!idQueue.empty())
  {
    Index currId = idQueue.front().first;
    level = idQueue.front().second;
    idQueue.pop();
    
    for(unsigned int i = 0; i < adjacencyGraph_[currId].size(); ++i)
    {
      Index adjId = adjacencyGraph_[currId][i];
      if(!foundIds[adjId])
      {
        idQueue.push(std::make_pair(adjId,level+1));
        bft_.push_back(adjId);
        foundIds[adjId] = true;
      }
    }

    // If adjacencyGraph_.size() isn't the full size of the problem, then there is a forest.
    // At this point, the level should be reset to 0, the current level should be stored in "levels",
    // and the search for the next viable root should start at "root".
    if(idQueue.empty() && ((int)bft_.size()!=numNodes()))
    {
      levels.push_back(level);    // Store level of completed tree.
      level = 0;                  // Initialize level since we are ordering a new tree.
     
      // Start from root, since we know everything before root has been seen 
      for( ; root < adjacencyGraph_.size(); ++root)
      {
        // we only count need to find Id's for places on the graph with greater than zero size
        if(adjacencyGraph_[root].size() > 0 && !foundIds[root])
        {
          idQueue.push(std::make_pair(root,level));
          bft_.push_back(root);
          foundIds[root] = true;
          break;
        }
      }
    }
  }

  int numIds = bft_.size();
  for(int i = 0; i < numIds; ++i)
    bftKeys_.push_back(keys1_[bft_[i]]);

  // The returned level is the maximum of all the trees, if there are multiple trees.
  if (levels.size() > 1)
  {
    std::vector<Index>::iterator max_level = std::max_element( levels.begin(), levels.end());
    level = *max_level;
  }
    
  return level;
}



#endif

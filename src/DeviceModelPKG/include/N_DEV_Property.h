//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2013  Sandia Corporation
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

/**
 * @file   N_DEV_Property.h
 * @author David G. Baur  KTech Corp.  Sandia National Laboratories 9143 <dgbaur@sandia.gov>
 * @date   Mon Apr 22 09:49:11 2013
 */
//-----------------------------------------------------------------------------
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.3.2.3 $
//
// Revision Date  : $Date: 2013/10/03 17:23:37 $
//
// Current Owner  : $Author: tvrusso $
//
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Property_h
#define Xyce_N_DEV_Property_h
#include <Xyce_config.h>

#ifdef HAVE_FUNCTIONAL
#include <functional>
#else
#ifdef HAVE_TR1_FUNCTIONAL
#include <tr1/functional>
#else
#error neither functional nor tr1/functional found.  Unable to proceed.
#endif
#endif
#ifdef HAVE_UNORDERED_MAP
#include <unordered_map>
#else
#ifdef HAVE_TR1_UNORDERED_MAP
#include <tr1/unordered_map>
#else
#error neither unordered_map nor tr1/unordered_map found.  Unable to proceed.
#endif
#endif

namespace Xyce {
namespace Device {


template<class C, class T>
class Entry;

/**
 * Property Entry Interface
 *
 *
 */
template<>
class Entry<void, void>
{
public:
/**
 * Returns the type_info of the data type being stored in the entry.
 *
 * @return
 */
  virtual const type_info &type() const = 0;
};

/**
 *
 *
 * @param memberValuePtr
 *
 * @return
 */
template<class T, class C>
class Entry : public Entry<void, void>
{
public:
  /**
   * Constructs a PropertyMap entry without tracking if value was given.
   *
   * @param memberValuePtr Member pointer containing value
   *
   * @return
   */
  Entry(T C::*memberValuePtr)
    : memberValuePtr_(memberValuePtr),
      memberGivenPtrSet_(false),
      memberGivenPtr_(0)
  {}

/**
 * Constructs a PropertyMap entry with tracking if value was given.
 *
 * @param memberValuePtr Member pointer containing vlaue
 * @param memberGivenPtr Member pointer containing given state
 *
 * @return
 */  Entry(T C::*memberValuePtr, bool C::*memberGivenPtr)
    : memberValuePtr_(memberValuePtr),
      memberGivenPtrSet_(true),
      memberGivenPtr_(memberGivenPtr)
  {}

/**
 * Returns type_info of this entry.
 *
 * @return type_info of this entry
 */
  const type_info &type() const {
    return typeid(T);
  }

/**
 * Returns the member pointer containing value.
 *
 * @return member pointer containing value.
 */
  const T C::*getMemberValuePtr() const {
    return memberValuePtr_;
  }

  T C::*getMemberValuePtr() {
    return memberValuePtr_;
  }

  /**
   * Returns the member pointer containing given state
   *
   * @return member pointer containing given state
   */
  const bool C::*getMemberGivenPtr() const {
    return memberGivenPtr_;
  }

  bool C::*getMemberGivenPtr() {
    return memberGivenPtr_;
  }

/**
 * Returns true if the member pointer given state is valid.
 *
 *
 * @return true if the member pointer given state is valid
 */
  bool getMemberGivenSet() const {
    return memberGivenPtrSet_;
  }

private:
  T C::*        memberValuePtr_;                //< Member pointer containing value
  bool C::*     memberGivenPtr_;                //< Member pointer containing given state
  bool          memberGivenPtrSet_;             //< True if memberGivenPtr_ is valid
};


template<class C>
class PropertyMap
{
public:
  typedef std::tr1::unordered_map<std::string, Entry<void, void> *> EntryMap;

  /**
   * Insert a new property into the property map.
   *
   * Property names are case sensitive and must be unique.
   *
   * @param name
   * @param member_pointer
   */
  template<class T>
  void insert(const std::string &name, T C::*member_pointer) {
    Entry<T, C> *entry = new Entry<T, C>(member_pointer);

    if (!entry_map.insert(typename EntryMap::value_type(name, entry)).second)
      throw std::runtime_error("Entry " + std::string(name) + " already exists in property map");
  }

  /**
   * Insert a new property into the property map and set *given_member_pointer to true when setValue() sets the value.
   *
   * Property names are case sensitive and must be unique.
   *
   * @param name
   * @param member_pointer
   * @param given_member_pointer
   */
  template<class T>
  void insert(const std::string &name, T C::*member_pointer, bool C::*given_member_pointer) {
    Entry<T, C> *entry = new Entry<T, C>(member_pointer, given_member_pointer);

    if (!entry_map.insert(typename EntryMap::value_type(name, entry)).second)
      throw std::runtime_error("Entry " + std::string(name) + " already exists in property map");
  }

  template<class T>
  Entry<T, C> *getEntry(const std::string &name) {
    typename EntryMap::const_iterator it = entry_map.find(name);
    if (it == entry_map.end())
      return 0;

    Entry<T, C> *entry = static_cast<Entry<T, C> *>((*it).second);
    if (entry->type() != typeid(T))
      throw std::runtime_error("Entry " + std::string(name) + " is of type " + entry->type().name() + ", type " + typeid(T).name() + " requested");

    return entry;
  }

  template<class T>
  const Entry<T, C> *getEntry(const std::string &name) const {
    typename EntryMap::const_iterator it = entry_map.find(name);
    if (it == entry_map.end())
      return 0;

    const Entry<T, C> *entry = static_cast<const Entry<T, C> *>((*it).second);
    if (entry->type() != typeid(T))
      throw std::runtime_error("Entry " + std::string(name) + " is of type " + entry->type().name() + ", type " + typeid(T).name() + " requested");

    return entry;
  }

private:
  EntryMap entry_map;
};

template<typename R, class P, class C>
const R &getValue(const P &property_map, const std::string &name, const C &instance) {
  const Entry<R, C> *entry = property_map.template getEntry<R>(name);
  if (!entry)
    throw std::runtime_error("Entry " + std::string(name) + " not found in " + instance.getName());

  const R C::*p = entry->getMemberValuePtr();

  return instance.*p;
}

template<typename T, class P, class C>
void setValue(P &property_map, const std::string &name, C &instance, T value) {
  Entry<T, C> *entry = property_map.template getEntry<T>(name);
  if (!entry)
    throw std::runtime_error("Entry " + std::string(name) + " not found in " + instance.getName());

  T C::*p = entry->getMemberValuePtr();
  instance.*p = value;

  if (entry->getMemberGivenSet()) {
    bool C::*g = entry->getMemberGivenPtr();
    instance.*g = true;
  }
}

template<typename R, class P, class C>
const bool &wasGiven(const P &property_map, const std::string &name, const C &instance) {
  Entry<R, C> *entry = property_map.template getEntry<R>(name);
  if (!entry)
    throw std::runtime_error("Entry " + std::string(name) + " not found in " + instance.getName());

  if (entry->getMemberGivenSet()) {
    bool C::*g = entry->getMemberGivenPtr();
    return instance.*g;
  }

  throw std::runtime_error("Entry " + std::string(name) + " is not tracking if value was given in " + instance.getName());
}

} // namespace Device
} // namespace Xyce

#endif // Xyce_N_DEV_Property_h

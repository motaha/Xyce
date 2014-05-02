#if 0
//-------------------------------------------------------------------------
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
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_UTL_Registry.C,v $
//
// Purpose        : 
//
// Special Notes  : 
//
// Creator        : David Baur
//
// Creation Date  : 3/20/2013
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.5 $
//
// Revision Date  : $Date: 2014/02/24 23:49:28 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------


#include <stdexcept>
#include <sstream>
#include <string>

#include <N_UTL_Registry.h>
#include <N_UTL_Demangle.h>

#ifdef XYCE_DLOPEN_ENABLED
#include <dlfcn.h>
#endif

namespace Xyce {
namespace Plugin {

Registry &
Registry::rootInstance()
{
  static Registry registry;

  return registry;
}



void
Registry::registerIt(
  const KeyPair &       key_pair,
  void *                func_ptr)
{
  // slibout.m(Slib::LOG_PLUGIN) << "Registering " << key_pair.second
  //                             << " of type " << demangle(key_pair.first->key())
  //                             << " at " << func_ptr << stk::diag::dendl;

  RegistryMap::const_iterator registry_entry = getRegistryMap().find(key_pair);
  if (registry_entry != getRegistryMap().end() && (*registry_entry).second != func_ptr) {
    std::ostringstream strout;
    strout << "Function with signature " << demangle((*registry_entry).first.first->key())
           << " and derived key '" << (*registry_entry).first.second
           << "' already registered to create function at address " << (*registry_entry).second;
    throw std::invalid_argument(strout.str());
  }
  getRegistryMap()[key_pair] = func_ptr;
}


void *
Registry::getPluginPtr(
  const KeyPair &      key_pair) const
{
  void *creator_function = getFuncPtr(key_pair);
  if (creator_function)
    return creator_function;
  else {
    std::ostringstream strout;

    strout << "User plugin creator function with base class '" << demangle(name_pair.first->name())
           << "' and derived class name '" << name_pair.second
           << "' not found in registry";
    throw std::invalid_argument(strout.str());
  }
}


void *
Registry::getFunctionPtr(
  const NamePair &      name_pair) const
{
  void *creator_function = getFuncPtr(name_pair);
  if (creator_function)
    return creator_function;
  else {
    std::ostringstream strout;
    strout << "User subroutine " << name_pair.second << "\n"
           << " with signature " << demangle(name_pair.first->name()) << "\n"
           << " not found in registry";
    throw std::invalid_argument(strout.str());
  }
}


Registry *
Registry::getFactoryPtr(
  const NamePair &      name_pair) const
{
  Registry *creator_function = (Registry *) getFuncPtr(name_pair);
  if (creator_function)
    return creator_function;
  else {
    std::ostringstream strout;
    strout << "Registry does not contain function with signature " << demangle(name_pair.first->name())
           << " and derived name '" << name_pair.second << "'";
    throw std::invalid_argument(strout.str());
  }
}


void *
Registry::getFuncPtr(
  const NamePair &      name_pair) const
{
  RegistryMap::const_iterator registry_entry = getRegistryMap().find(name_pair);
  return registry_entry == getRegistryMap().end() ? NULL : (*registry_entry).second;
}


std::vector<std::string>
Registry::getDerivedNames(
  const std::type_info &        type) const
{
  std::vector<std::string> derived_names;

  for (RegistryMap::const_iterator it = getRegistryMap().begin(); it != getRegistryMap().end(); ++it)
    if (*(*it).first.first == type)
      derived_names.push_back((*it).first.second);

  return derived_names;
}


typedef void (*dl_register_t)();

void
Registry::registerDL(
  const char *          so_path,
  const char *          function_name)
{
#ifdef XYCE_DLOPEN_ENABLED
  slibout.m(Slib::LOG_PLUGIN) << "Loading dynamic library " << so_path << stk::diag::dendl;
  void *dl = dlopen(so_path, RTLD_NOW);
  if (!dl){
    throw std::runtime_error(dlerror());
  }

  if (function_name) {
    std::string s = std::strlen(function_name) ? function_name : "dl_register";

    dl_register_t f = (dl_register_t) dlsym(dl, s.c_str());
    if (!f) {
      s = s + XYCE_FORTRAN_SUFFIX;

      f = (dl_register_t) dlsym(dl, s.c_str());
    }

    if (f) {
      slibout.m(Slib::LOG_PLUGIN) << "Executing dynamic library " << so_path << " function " << s << "()" << stk::diag::dendl;
      (*f)();
    }
    else {
      if (std::strlen(function_name)) {
        std::ostringstream str;
        str << "Registration function " << function_name << " not found in " << so_path;
        throw std::runtime_error(str.str().c_str());
      }
    }
  }

#else
  throw std::runtime_error("Dynamic linkage is not supported on this platform");
#endif
}


template <>
void *
Registry::getsym<void *>(
  const char *  sym)
{
#ifdef XYCE_DLOPEN_ENABLED
  void *s = NULL;
  void *dl = dlopen(NULL, RTLD_LAZY);
  if (dl) {
    s = dlsym(dl, sym);
    dlclose(dl);
  }

  return s;
#else
  return NULL;
#endif
}


std::ostream &
Registry::verbose_print(
  std::ostream &                os) const
{
  for (RegistryMap::const_iterator it = getRegistryMap().begin(); it != getRegistryMap().end(); ++it)
    os << (*it).first.second << " of type " << demangle((*it).first.first->name()) << " at " << (*it).second << std::endl;
  return os;
}

} // namespace Plugin
} // namespace Xyce

#endif

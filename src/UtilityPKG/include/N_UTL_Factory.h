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
// Filename       : $RCSfile: N_UTL_Factory.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Robert Hoekstra
//
// Creation Date  : 2/07/02
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.12 $
//
// Revision Date  : $Date: 2014/02/24 23:49:28 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef  _FACTORY_H
#define  _FACTORY_H

// ---------- Standard Includes ----------

#include <map>
#include <vector>
#include <string>

// ----------   Xyce Includes   ----------

// ----------   Fwd Declares    ----------

namespace XUTL
{

  typedef std::string defaultKeyType;

  //---------------------------------------------------------------------------
  // Class         : Factory
  // Purpose       : Templated Factory Implementation
  // Special Notes :
  // Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
  // Creation Date : 2/07/02
  //---------------------------------------------------------------------------
  template < class manufacturedType, typename keyType = defaultKeyType >
  class Factory
  {

    // Use auto_ptr if available.
    typedef manufacturedType * manufacturedTypePtr;
    typedef manufacturedTypePtr(* CreateFn) ();
    typedef std::map < keyType, CreateFn > FnRegistry;
    FnRegistry registry_;

    // Default constructor.
    Factory() { }

    // Copy constructor - not implemented.
    Factory(const Factory &);

    // Assignment operator - not implemented.
    Factory & operator = (const Factory &);

  public:

    // Singleton
    static Factory & instance()
    {
      static Factory bf;
      return bf;
    }

    bool RegCreateFn(const keyType &, CreateFn);

    manufacturedTypePtr Create(const keyType &);
  };

  // Factory Methods:

  template < class manufacturedType, typename keyType >
  bool Factory < manufacturedType,
                 keyType >::RegCreateFn(const keyType & className,
                                        CreateFn func)
  {
    registry_[className] = func;
    return (func != static_cast < CreateFn > (0));
  }

  template < class manufacturedType, typename keyType >
  typename Factory < manufacturedType,
                     keyType >::manufacturedTypePtr Factory < manufacturedType,
                                                              keyType >::Create(const keyType & className)
  {
    if (registry_.count(className))
      return registry_[className] ();
    return manufacturedTypePtr(0);
  }

  //---------------------------------------------------------------------------
  // Class         : RegisterInFactory
  // Purpose       : Templated Factory Helper
  // Special Notes :
  // Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
  // Creation Date : 2/07/02
  //---------------------------------------------------------------------------
  template < class ancestorType, class manufacturedType, typename keyType =
             defaultKeyType >
  class RegisterInFactory
  {

    typedef ancestorType * ancestorTypePtr;

  public:
    static ancestorTypePtr CreateInstance();

    RegisterInFactory(const keyType &);
  };

  // RegisterInFactory Methods:

  template < class ancestorType, class manufacturedType, typename keyType >
  typename RegisterInFactory < ancestorType, manufacturedType,
  	                       keyType >::ancestorTypePtr RegisterInFactory < ancestorType,
  	                                                                      manufacturedType, keyType >::CreateInstance()
  {
    return ancestorTypePtr(new manufacturedType);
  }

  template < class ancestorType, class manufacturedType, typename keyType >
  RegisterInFactory < ancestorType, manufacturedType,
                      keyType >::RegisterInFactory(const keyType & id)
  {
    Factory < ancestorType, keyType >::instance().RegCreateFn(id, CreateInstance);
  }

}

#endif

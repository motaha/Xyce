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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_Param.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric Keiter
//
// Creation Date  : 5/22/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.21.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:37 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef  N_DEV_PARAM_H
#define  N_DEV_PARAM_H

// ---------- Standard Includes ----------

#include <iosfwd>
#include <string>
#include <vector>
#include <map>

#include <N_DEV_fwd.h>
#include <N_UTL_NoCase.h>
#include <N_UTL_Param.h>

// ----------   Forward Declarations   ----------


namespace Xyce {
namespace Device {

class ParamData;

template <class T>
struct DataTypeTrait;

template<>
struct DataTypeTrait<std::string>
{
    enum {type = STR};
};

template<>
struct DataTypeTrait<double>
{
    enum {type = DBLE};
};

template<>
struct DataTypeTrait<int>
{
    enum {type = INT};
};

template<>
struct DataTypeTrait<long>
{
    enum {type = LNG};
};

template<>
struct DataTypeTrait<bool>
{
    enum {type = BOOL};
};

template<>
struct DataTypeTrait<std::vector<std::string> >
{
    enum {type = STR_VEC};
};

template<>
struct DataTypeTrait<std::vector<int> >
{
    enum {type = INT_VEC};
};

template<>
struct DataTypeTrait<std::vector<double> >
{
    enum {type = DBLE_VEC};
};

template<>
struct DataTypeTrait<CompositeMap>
{
    enum {type = COMPOSITE};
};

//-----------------------------------------------------------------------------
// Class         : N_DEV_Param
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/10/01
//-----------------------------------------------------------------------------
class Param : public N_UTL_Param
{
public:

  // Constructors
  Param();
  Param( const std::string & t, const string & v, const bool & g = false );
  Param( const std::string & t, const double & v, const bool & g = false );
  Param( const std::string & t, const int & v, const bool & g = false );
  Param( const std::string & t, const long & v, const bool & g = false );
  Param( const std::string & t, const bool & v, const bool & g = false );
  Param( const std::string & t, const char * v, const bool & g = false );
  Param( const std::string & t, const std::vector<std::string> & v, const bool & g = false );
  Param( const std::string & t, const vector<double> & v, const bool & g = false );
  Param( Param const& rhsParam );

  Param & operator=(Param const& rhsParam);

  // Destructor
  ~Param();

  // Methods to set reset value
  void setGiven( const bool & g );
  void setDefault( const bool & d );

  // Methods to get values
  const bool & given() const;
  const bool & default_val() const;

  Packable * instance() const;
  int packedByteCount() const;

  void pack( char * buf, int bsize, int & pos, N_PDS_Comm * comm ) const;
  void unpack( char * pB, int bsize, int & pos, N_PDS_Comm * comm );

  void print();

private:

  ParamData* data_;

};

typedef std::map<std::string, std::vector<Param>, LessNoCase> DeviceParamMap;

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::Param N_DEV_Param;

#endif



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
// Filename      : $RCSfile: N_DEV_Param.C,v $
//
// Purpose       :
//
// Special Notes :
//
// Creator       : Robert Hoekstra, SNL
//
// Creation Date : 5/15/01
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.19.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:38 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------
#include <iostream>

// ----------   Xyce Includes   ----------

#include <N_DEV_Param.h>
#include <N_PDS_Comm.h>
#include <N_ERH_ErrorMgr.h>

// ---------  Other Includes  -----------

// ---------  Helper Classes ------------

namespace Xyce {
namespace Device {

class ParamData
{
  public:
    ParamData( bool g = false ) :
      given_(g), default_(false)
    {}

    bool given_;
    bool default_;
};

//-----------------------------------------------------------------------------
// Function      : Param::Param
// Purpose       : Constructor
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 5/16/01
//-----------------------------------------------------------------------------
Param::Param() : N_UTL_Param()
 , data_( new ParamData )
{
}

//-----------------------------------------------------------------------------
// Function      : Param::Param
// Purpose       : Constructor
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 5/16/01
//-----------------------------------------------------------------------------
Param::Param( const string & t, const string & v,
                          const bool & g ) : N_UTL_Param( t, v )
 , data_( new ParamData(g) )
{
}

//-----------------------------------------------------------------------------
// Function      : Param::Param
// Purpose       : Constructor
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 5/16/01
//-----------------------------------------------------------------------------
Param::Param( const string & t, const double & v,
                          const bool & g ) : N_UTL_Param( t, v )
 , data_( new ParamData(g) )
{
}

//-----------------------------------------------------------------------------
// Function      : Param::Param
// Purpose       : Constructor
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 5/16/01
//-----------------------------------------------------------------------------
Param::Param( const string & t, const int & v,
                          const bool & g ) : N_UTL_Param( t, v )
 , data_( new ParamData(g) )
{
}

//-----------------------------------------------------------------------------
// Function      : Param::Param
// Purpose       : Constructor
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 03/14/06
//-----------------------------------------------------------------------------
Param::Param( const string & t, const long & v,
                          const bool & g ) : N_UTL_Param( t, v )
 , data_( new ParamData(g) )
{
}

//-----------------------------------------------------------------------------
// Function      : Param::Param
// Purpose       : Constructor
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 03/14/06
//-----------------------------------------------------------------------------
Param::Param( const string & t, const bool & v,
                          const bool & g ) : N_UTL_Param( t, v )
 , data_( new ParamData(g) )
{
}

//-----------------------------------------------------------------------------
// Function      : Param::Param
// Purpose       : Constructor
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 03/14/06
//-----------------------------------------------------------------------------
Param::Param( const string & t, const char * v,
                          const bool & g ) : N_UTL_Param( t, string(v) )
 , data_( new ParamData(g) )
{
}


//-----------------------------------------------------------------------------
// Function      : Param::Param
// Purpose       : Constructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Modeling, 1445
// Creation Date : 12/23/10
//-----------------------------------------------------------------------------
Param::Param( const string & t, const vector<string> & v, const bool & g  )
  : N_UTL_Param( t, v ) , data_( new ParamData(g) )
{
}

//-----------------------------------------------------------------------------
// Function      : Param::Param
// Purpose       : Constructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical Modeling, 1445
// Creation Date : 12/23/10
//-----------------------------------------------------------------------------
Param::Param( const string & t, const vector<double> & v, const bool & g )
  : N_UTL_Param( t, v ) , data_( new ParamData(g) )
{
}


//-----------------------------------------------------------------------------
// Function      : Param::Param
// Purpose       : Copy Constructor
// Special Notes :
// Scope         : public
// Creator       : Lon Waters, SNL
// Creation Date : 5/16/01
//-----------------------------------------------------------------------------
Param::Param( Param const& rhsParam )
  : N_UTL_Param( rhsParam )
 , data_( new ParamData(rhsParam.data_->given_) )
{
  data_->default_ = rhsParam.data_->default_;
}

//-----------------------------------------------------------------------------
// Function      : Param::operator=
// Purpose       : assignment operator
// Special Notes :
// Scope         :
// Creator       : Lon Waters, SNL
// Creation Date : 05/18/01
//-----------------------------------------------------------------------------
Param & Param::operator=(Param const& rhsParam)
{
  N_UTL_Param::operator=(rhsParam);
  data_->given_ = rhsParam.data_->given_;
  data_->default_ = rhsParam.data_->default_;

  return *this;
}

//-----------------------------------------------------------------------------
// Function      : Param::~Param
// Purpose       : Destructor
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL
// Creation Date : 5/16/01
//-----------------------------------------------------------------------------
Param::~Param()
{
  if( data_ != NULL ) delete data_;
  data_ = NULL;
}

//-----------------------------------------------------------------------------
// Function      : Param::setGiven
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Rob Hoekstra, SNL
// Creation Date : 05/18/01
//-----------------------------------------------------------------------------
void Param::setGiven( const bool & g )
{
  data_->given_ = g;
}

//-----------------------------------------------------------------------------
// Function      : Param::setDefault
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/11/05
//-----------------------------------------------------------------------------
void Param::setDefault( const bool & d )
{
  data_->default_ = d;
}

//-----------------------------------------------------------------------------
// Function      : Param::given
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Rob Hoekstra, SNL
// Creation Date : 05/18/01
//-----------------------------------------------------------------------------
const bool & Param::given() const
{
  return data_->given_;
}

//-----------------------------------------------------------------------------
// Function      : Param::default_val
// Purpose       :
// Special Notes :
// Scope         :
// Creator       : Dave Shirley, PSSI
// Creation Date : 01/11/05
//-----------------------------------------------------------------------------
const bool & Param::default_val() const
{
  return data_->default_;
}

//-----------------------------------------------------------------------------
// Function      : instance
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL
// Creation Date : 5/15/01
//-----------------------------------------------------------------------------
Packable * Param::instance() const
{
  return new Param();
}

//-----------------------------------------------------------------------------
// Function      : packedByteCount
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL
// Creation Date : 5/15/01
//-----------------------------------------------------------------------------
int Param::packedByteCount() const
{
  //N_UTL_Param info
  int byteCount = N_UTL_Param::packedByteCount();

  //given & default
  byteCount += sizeof(int);

  return byteCount;

}

//-----------------------------------------------------------------------------
// Function      : pack
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL
// Creation Date : 5/15/01
//-----------------------------------------------------------------------------
void Param::pack( char * buf, const int bsize, int & pos,
		N_PDS_Comm * comm ) const
{
#ifdef Xyce_COUNT_PACKED_BYTES
  int predictedPos = pos+packedByteCount();
#endif

  N_UTL_Param::pack( buf, bsize, pos, comm );

  //pack given_
  int dg = (data_->given_?1:0) + 2*(data_->default_?1:0);
  comm->pack( &dg, 1, buf, bsize, pos );

#ifdef Xyce_COUNT_PACKED_BYTES
  if (pos != predictedPos)
  {
    string msg = "Predicted pos does not match actual pos in Param::pack";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING, msg );
  }
#endif
}

//-----------------------------------------------------------------------------
// Function      : unpack
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Rob Hoekstra, SNL
// Creation Date : 5/15/01
//-----------------------------------------------------------------------------
void Param::unpack( char * pB, int bsize, int & pos, N_PDS_Comm * comm )
{

  N_UTL_Param::unpack( pB, bsize, pos, comm );

  //unpack given_
  int dg;
  comm->unpack( pB, bsize, pos, &dg, 1 );
  data_->given_ = ( dg%2 != 0 );
  data_->default_ = ( dg >= 2 );

}

//-----------------------------------------------------------------------------
// Function      : print
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Dave Shirley, PSSI
// Creation Date : 12/09/05
//-----------------------------------------------------------------------------
void Param::print()
{

  N_UTL_Param::print();
  cout << "Given: " << data_->given_ << "  Default: " << data_->default_ << endl;

}

} // namespace Device
} // namespace Xyce

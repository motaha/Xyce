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
// Filename       : $RCSfile: N_UTL_OptionBlock.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Rob Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 05/15/01
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.26 $
//
// Revision Date  : $Date: 2014/02/24 23:49:28 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>

#include <N_UTL_OptionBlock.h>
#include <N_PDS_Comm.h>
#include <N_ERH_Message.h>

namespace Xyce {
namespace Util {

//-----------------------------------------------------------------------------
// Function      : OptionBlock::OptionBlock
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/15/01
//-----------------------------------------------------------------------------
OptionBlock::OptionBlock(const std::string & n_str)
  : name_(n_str),
    netlistLocation_(),
    status_(NEW_STATE),
    params_()
{
}

//-----------------------------------------------------------------------------
// Function      : OptionBlock::OptionBlock
// Purpose       : copy constructor
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/15/01
//-----------------------------------------------------------------------------
OptionBlock::OptionBlock(const OptionBlock & right)
  : name_(right.name_),
    netlistLocation_(right.netlistLocation_),
    status_(right.status_),
    params_(right.params_)
{
}

//-----------------------------------------------------------------------------
// Function      : OptionBlock::~OptionBlock
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/15/01
//-----------------------------------------------------------------------------
OptionBlock::~OptionBlock()
{}


//-----------------------------------------------------------------------------
// Function      : OptionBlock::compareParamLists
// Purpose       : Compares contained Param lists
// Special Notes :
// Scope         : public
// Creator       : Rich Schiek, Electrical Systems Modeling
// Creation Date : 2/6/12
//-----------------------------------------------------------------------------
bool OptionBlock::compareParamLists( const OptionBlock & right ) const
{
  bool match = true;
  if( params_.size() == right.params_.size() )
  {
    // length of params list is the same, so they may match
    std::list< Param >::const_iterator existingListItr = params_.begin();
    std::list< Param >::const_iterator existingListEnd = params_.end();
    std::list< Param >::const_iterator newListItr = right.params_.begin();
    std::list< Param >::const_iterator newListEnd = right.params_.end();
    while( existingListItr != existingListEnd )
    {
      if( !(existingListItr->deepCompare(*newListItr)) )
      {
        // Param objects didn't match
        // thus these two lists are not the same
        match = false;
        break;
      }
      existingListItr++;
      newListItr++;
    }
  }
  else
  {
    match = false;
  }
  return match;
}


//-----------------------------------------------------------------------------
// Function      : OptionBlock::operator=
// Purpose       : "=" operator
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/15/01
//-----------------------------------------------------------------------------
OptionBlock & OptionBlock::operator=(const OptionBlock & right)
{
  name_   = right.name_;
  status_ = right.status_;
  params_ = right.params_;

  return *this;
}

//-----------------------------------------------------------------------------
// Function      : OptionBlock::operator==
// Purpose       : "==" operator
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/15/01
//-----------------------------------------------------------------------------
int OptionBlock::operator==(const OptionBlock & right) const
{
  return name_ == right.name_;
}

//-----------------------------------------------------------------------------
// Function      : OptionBlock::operator!=
// Purpose       : "!=" operator
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/15/01
//-----------------------------------------------------------------------------
int OptionBlock::operator!=(const OptionBlock & right) const
{
  return name_ != right.name_;
}

//-----------------------------------------------------------------------------
// Function      : OptionBlock::clear
// Purpose       : Empties out the block.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/15/01
//-----------------------------------------------------------------------------
void OptionBlock::clear()
{
  name_ = "";
  status_ = NEW_STATE;

  params_.clear();
}

//-----------------------------------------------------------------------------
// Function      : OptionBlock::instance
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/28/01
//-----------------------------------------------------------------------------
Packable * OptionBlock::instance() const
{
  return new OptionBlock();
}

//-----------------------------------------------------------------------------
// Function      : OptionBlock::packedByteCount
// Purpose       : Counts bytes needed to pack block.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/15/01
//-----------------------------------------------------------------------------
int OptionBlock::packedByteCount() const
{

  int byteCount = 0;

  int size, length, i;

  std::list<Param>::const_iterator it_tpL;

  //----- count name
  length = name_.length();
  byteCount += sizeof(int);
  byteCount += length;

  //----- count netlistLocation
  length = netlistLocation_.getPath().length();
  byteCount += sizeof(int);
  byteCount += length;
  byteCount += sizeof(int);

  //----- count status
  byteCount += sizeof(int);

  //----- count params
  size = params_.size();
  byteCount += sizeof(int);
  it_tpL = params_.begin();
  for (i = 0; i < size; ++i, ++it_tpL)
    byteCount += it_tpL->packedByteCount();

  return byteCount;

}

//-----------------------------------------------------------------------------
// Function      : OptionBlock::pack
// Purpose       : Packs OptionBlock into char buffer using MPI_PACK.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/15/01
//-----------------------------------------------------------------------------
void OptionBlock::pack(char * buf, int bsize, int & pos, N_PDS_Comm * comm) const
{

#ifdef Xyce_DEBUG_TOPOLOGY
  Xyce::dout() << "Packing OptionBlock: " << getName() << std::endl;
#endif

  int size, length;
  int i;
  std::list<Param>::const_iterator it_tpL;

  //----- pack name
  length = name_.length();
  comm->pack( &length, 1, buf, bsize, pos );
  comm->pack( name_.c_str(), length, buf, bsize, pos );

  //----- pack netlist location
  length = netlistLocation_.getPath().length();
  comm->pack( &length, 1, buf, bsize, pos );
  comm->pack( netlistLocation_.getPath().c_str(), length, buf, bsize, pos );
  int line_number = netlistLocation_.getLineNumber();
  comm->pack( &line_number, 1, buf, bsize, pos );

  //----- pack status
  comm->pack( &status_, 1, buf, bsize, pos );

  //----- pack params
  size = params_.size();
  comm->pack( &size, 1, buf, bsize, pos );
  for (i = 0, it_tpL = params_.begin(); i < size; ++i, ++it_tpL)
    it_tpL->pack( buf, bsize, pos, comm );

#ifdef Xyce_DEBUG_TOPOLOGY
  Xyce::dout() << "Packed " << pos << " of " << bsize << " bytes for OptionBlock: " <<
    getName() << std::endl;
#endif

}

//-----------------------------------------------------------------------------
// Function      : OptionBlock::unpack
// Purpose       : Unpacks OptionBlock from char buffer using MPI_UNPACK.
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/15/01
//-----------------------------------------------------------------------------
void OptionBlock::unpack(char * pB, int bsize,int & pos, N_PDS_Comm * comm)
{

  int size, length;
  int i;

  //----- unpack name
  comm->unpack( pB, bsize, pos, &length, 1 );
  name_ = std::string( (pB+pos), length);
  pos += length;

#ifdef Xyce_DEBUG_TOPOLOGY
  Xyce::dout() << "Unpacking OptionBlock: " << getName() << std::endl;
#endif

  //----- unpack netlist location
  comm->unpack( pB, bsize, pos, &length, 1 );
  netlistLocation_.setPath(std::string( (pB+pos), length));
  pos += length;
  int line_number;
  comm->unpack( pB, bsize, pos, &line_number, 1 );
  netlistLocation_.setLineNumber(line_number);

  //----- unpack status
  comm->unpack( pB, bsize, pos, &status_, 1 );

  //----- unpack params
  comm->unpack( pB, bsize, pos, &size, 1 );

  params_.clear();
  Param param;
  for( i = 0; i < size; ++i )
  {
    param.unpack( pB, bsize, pos, comm );
    params_.push_back( param );
  }

#ifdef Xyce_DEBUG_TOPOLOGY
  Xyce::dout() << "Unpacked " << pos << " of " << bsize << " bytes for OptionBlock: " <<
    getName() << std::endl;
#endif

}


//-----------------------------------------------------------------------------
// Function      : OptionBlock::tagExists
// Purpose       : Searches param list and returns true if tag exists
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystem Modeling
// Creation Date : 3/11/2009
//-----------------------------------------------------------------------------
bool OptionBlock::tagExists( const std::string iTag ) const
{
  bool retValue=false;

  std::list<Param>::const_iterator paramListIt = params_.begin();
  std::list<Param>::const_iterator paramListItEnd = params_.end();
  while( paramListIt != paramListItEnd )
  {
    if(paramListIt->tag() == iTag)
    {
      retValue = true;
      break;
    }
    paramListIt++;
  }

  return retValue;
}


//-----------------------------------------------------------------------------
// Function      : OptionBlock::getTagValueAsString
// Purpose       : Searches param list and returns tag's value as a string
// Special Notes : empty string returned if tag isn't found
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystem Modeling
// Creation Date : 3/11/2009
//-----------------------------------------------------------------------------
std::string OptionBlock::getTagValueAsString( std::string iTag ) const
{
  std::string retValue;

  std::list<Param>::const_iterator paramListIt = params_.begin();
  std::list<Param>::const_iterator paramListItEnd = params_.end();
  while( paramListIt != paramListItEnd )
  {
    if(paramListIt->tag() == iTag)
    {
      retValue = paramListIt->stringValue();
      break;
    }
    paramListIt++;
  }

  return retValue;
}

//-----------------------------------------------------------------------------
// Function      : OptionBlock::getTagValueAsDouble
// Purpose       : Searches param list and returns tag's value as a double
// Special Notes : zero returned if tag is not found
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystem Modeling
// Creation Date : 3/11/2009
//-----------------------------------------------------------------------------
double OptionBlock::getTagValueAsDouble( std::string iTag ) const
{
  double retValue;

  std::list<Param>::const_iterator paramListIt = params_.begin();
  std::list<Param>::const_iterator paramListItEnd = params_.end();
  while( paramListIt != paramListItEnd )
  {
    if(paramListIt->tag() == iTag)
    {
      retValue = paramListIt->getImmutableValue<double>();
      break;
    }
    paramListIt++;
  }

  return retValue;
}

//-----------------------------------------------------------------------------
// Function      : OptionBlock::getTagValueAsInt
// Purpose       : Searches param list and returns tag's value as a real
// Special Notes : zero it returned if tag is not found
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystem Modeling
// Creation Date : 3/11/2009
//-----------------------------------------------------------------------------
int OptionBlock::getTagValueAsInt( std::string iTag ) const
{
  int retValue;

  std::list<Param>::const_iterator paramListIt = params_.begin();
  std::list<Param>::const_iterator paramListItEnd = params_.end();
  while( paramListIt != paramListItEnd )
  {
    if(paramListIt->tag() == iTag)
    {
      retValue = paramListIt->getImmutableValue<int>();
      break;
    }
    paramListIt++;
  }

  return retValue;
}

//-----------------------------------------------------------------------------
// Function      : OptionBlock::operator<<
// Purpose       : "<<" operator
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/15/01
//-----------------------------------------------------------------------------
std::ostream & operator<<(std::ostream & os, const OptionBlock & mb)
{
  std::list<Param>::const_iterator it_pL, end_pL;

  os << "Option Block" << std::endl;
  os << " name:  " << mb.getName() << std::endl;
  os << " status: " << mb.getStatus() << std::endl;
  os << " Params" << std::endl;
  os << " -------------" << std::endl;

  it_pL=mb.getParams().begin();
  end_pL=mb.getParams().end();
  for ( ; it_pL != end_pL; ++it_pL)
  {
    os << *it_pL;
  }
  os << " -------------" << std::endl;
  os << std::endl;

  return os;
}

} // namespace Device
} // namespace Xyce

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

//-----------------------------------------------------------------------------
// Filename       : $RCSfile: N_UTL_OptionBlock.h,v $
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
// Revision Number: $Revision: 1.11.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:52 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_UTL_OptionBlock_h
#define Xyce_N_UTL_OptionBlock_h

// ---------- Standard Includes ----------
#include <string>
#include <list>
#include <iosfwd>

// ----------   Xyce Includes   ----------
#include <N_UTL_fwd.h>
#include <N_UTL_Param.h>
#include <N_UTL_Misc.h>
#include <N_UTL_Xyce.h>
#include <N_UTL_Packable.h>

namespace Xyce {
namespace Util {

// this enum defines the current status of an OptionBlock
// it's new when parsed and not yet sent to the package that needs it.
// it is PROCESSED after it has been communicated to a given package.
// it is MODIFIED if later parsing has changed it due to parameter or
// expression resolution.  MODIFIED blocks are re-sent to packages and
// then marked as PROCESSED.  Thus packages can know if an options block
// is being resent because of resolved item.
enum { NEW_STATE, MODIFIED_STATE, PROCESSED_STATE };

//-----------------------------------------------------------------------------
// Class         : N_DEV_OptionBlock
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/15/01
//-----------------------------------------------------------------------------
class OptionBlock : public Packable
{
  public:
    typedef std::list<N_UTL_Param> ParameterList;

    // Constructor
    OptionBlock(const std::string & n_str = std::string());

    // Copy constructor.
    OptionBlock(const OptionBlock & right);

    // Destructor
    ~OptionBlock();

    // compare contained lists
    // the equity operator just compares the "name" field.  This
    // will compare the N_UTL_Param lists as well.
    bool compareParamLists( const OptionBlock & right ) const;

    // Assignment operator.
    OptionBlock & operator = (const OptionBlock & right);

    // Equality operator.
    int operator == (const OptionBlock & right) const;

    // Non-equality operator.
    int operator != (const OptionBlock & right) const;

    // Empties out the block.
    void clear();

    //beginning of params
    ParameterList::iterator begin() {
      return params_.begin();
    }

    //end of params
    ParameterList::iterator end() {
      return params_.end();
    }

    // Packing functionality.
    Packable * instance() const;

    // Counts bytes needed to pack block.
    int packedByteCount() const;

    // Packs OptionBlock into char buffer using MPI_PACK.
    void pack(char * buf, int bsize, int & pos, N_PDS_Comm * comm) const;

    // Unpacks OptionBlock from char buffer using MPI_UNPACK.
    void unpack(char * pB, int bsize, int & pos, N_PDS_Comm * comm);

    // simpe get routines
    bool tagExists( const std::string iTag ) const;
    std::string getTagValueAsString( std::string iTag ) const;
    double getTagValueAsDouble( std::string iTag ) const;
    int getTagValueAsInt( std::string iTag ) const;

    const std::list<N_UTL_Param> &getParams() const {
      return params_;
    }

    std::list<N_UTL_Param> &getParams() {
      return params_;
    }

    const std::string &getName() const {
      return name_;
    }

    void setName(const std::string &name) {
      name_ = name;
    }

    int getStatus() const {
      return status_;
    }

    void setStatus(int status) {
      status_ = status;
    }

  private:

    // Not encapsulated, probably should be.
    std::string name_;
    int status_;  // state of block: NEW_STATE, MODIFIED_STATE, PROCESSED_STATE
    list < N_UTL_Param > params_;
};

std::ostream &operator<<(std::ostream &os, const OptionBlock &options);

} // namespace Device
} // namespace Xyce


typedef Xyce::Util::OptionBlock N_UTL_OptionBlock;

#endif

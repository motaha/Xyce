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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_UTL_NetlistLocation.h,v $
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
// Revision Number: $Revision: 1.2.2.1 $
//
// Revision Date  : $Date: 2014/02/26 19:14:55 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef Xyce_N_UTL_NetlistLocation_h
#define Xyce_N_UTL_NetlistLocation_h

#include <iosfwd>
#include <string>

namespace Xyce {

struct NetlistLocation
{
  NetlistLocation()
    : path_(),
      lineNumber_(0)
  {}

  NetlistLocation(const std::string &path, const int line_number)
    : path_(path),
      lineNumber_(line_number)
  {}

  NetlistLocation(const NetlistLocation &netlist_location)
    : path_(netlist_location.path_),
      lineNumber_(netlist_location.lineNumber_)
  {}

  NetlistLocation &operator=(const NetlistLocation &netlist_location) 
  {
    path_ = netlist_location.path_;
    lineNumber_ = netlist_location.lineNumber_;

    return *this;
  }

  NetlistLocation &setPath(const std::string &path) 
  {
    path_ = path;
    return *this;
  }

  const std::string &getPath() const 
  {
    return path_;
  }

  NetlistLocation &setLineNumber(int line_number) 
  {
    lineNumber_ = line_number;
    return *this;
  }

  int getLineNumber() const 
  {
    return lineNumber_;
  }

private:
  std::string   path_;
  int           lineNumber_;
};

std::ostream &operator<<(std::ostream &os, const NetlistLocation &netlist_locationc);

} // namespace Xyce

#endif // Xyce_N_UTL_NetlistLocation_h


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
// Filename       : $RCSfile: N_ERH_Message.C,v $
//
// Purpose        : 
//
// Special Notes  : 
//
// Creator        : David Baur
//
// Creation Date  : 1/7/2014
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision $
//
// Revision Date  : $Date: 2014/02/27 18:38:36 $
//
// Current Owner  : $Author $
//-------------------------------------------------------------------------


#include <list>
#include <string>
#include <sstream>
#include <utility>
#include <vector>

#include <N_ERH_Message.h>
#include <N_ERH_Messenger.h>

namespace Xyce {
namespace Report {

namespace {

//-----------------------------------------------------------------------------
// Function      : prefix
// Purpose       : Construct message parameters based on the report_mask and the
//                 message content.
// Special Notes :
// Scope         : Private
// Creator       : David Baur
// Creation Date : 1/7/2014
//-----------------------------------------------------------------------------
std::ostream &prefix(std::ostream &os, unsigned report_mask)
{
  
  if (report_mask & MSG_USER)
    os << "User ";

  if (report_mask & MSG_DEVEL)
    os << "Developer ";

  if ((report_mask & MSG_TYPE_MASK) == MSG_FATAL)
    os << "error";

  if ((report_mask & MSG_TYPE_MASK) == MSG_ERROR)
    os << "error";

  if ((report_mask & MSG_TYPE_MASK) == MSG_WARNING)
    os << "warning";

  if ((report_mask & MSG_TYPE_MASK) == MSG_INFORMATION)
    os << "info";

  if ((report_mask & MSG_TYPE_MASK) == MSG_DEBUG)
    os << "debug";

  return os;
}

}  // namespace <unnamed>

MessageCode
MessageCode::s_defaultMessageCode(100000000);

Message::Message(
  MessageType   message_type,
  MessageCode   message_code)
  : messageType_(message_type),
    messageCode_(message_code),
    oss_(),
    netlistLocation_(),
    functionName_(0)
{}

static const bool PARSABLE_LINE = false;

Message::~Message()
{
  std::ostringstream os;

  if (PARSABLE_LINE) {
    if (netlistLocation_.getLineNumber() > 0)
      os << netlistLocation_.getPath() << ":" << netlistLocation_.getLineNumber() << ": ";

    prefix(os, messageType_);
    os << ": ";
  }
  else {
    prefix(os, messageType_);

    if (netlistLocation_.getLineNumber() > 0)
      os << " in file " << netlistLocation_.getPath() << " at or near line " << netlistLocation_.getLineNumber() << "\n";
    else
      os << ": ";
  }

  if (functionName_)
    os << "function " << functionName_ << ":\n";

  os << oss_.str();

  Xyce::Report::report_message(os.str().c_str(), messageType_, messageCode_);
}

} // namespace Report
} // namespace Xyce

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
// Filename       : $RCSfile: N_ERH_Message.h,v $
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
// Revision Number: $Revision: 1.8.2.2 $
//
// Revision Date  : $Date: 2014/03/03 18:29:27 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------


#ifndef Xyce_N_ERH_Message_h
#define Xyce_N_ERH_Message_h

#include <string>
#include <sstream>
#include <cstddef>

#include <N_ERH_fwd.h>

#include <N_UTL_NetlistLocation.h>

namespace Xyce {
namespace Report {

/**
 * @brief Typedef <b>MessageId</b> defines a message identifier.
 *
 * Message identifiers must be consist from reference to reference, unique for each instance, and
 * yet consist within each instance across multiple processors.  To meet these criteria, the message
 * identifier is implemented as a static memory location.  It must be declared by the application
 * developer as static.  This results in the linker selecting an address for each instance, which
 * never changes from reference to reference, is unique for each instance and the same regardless of
 * executable (assuming the executable is mapped into the same memory location for each process).
 * In order to remove the pointer-ness of the static memory location, the address is cast to a
 * pointer difference type which is an integral type by subtracting the zero pointer from it.
 *
 */
typedef std::ptrdiff_t MessageId;

/**
 * @brief Enumeration <b>ThrottleGroup</b> lists defined throttling groups.
 *
 * When messages are throttled, the throttle count may be reset at varior points during an
 * application run.  Some throttles defined for the application, while other may be reset at each
 * time step or other interval.  This allows warnings to be repeated at each time step rather than
 * cut off.
 *
 */
enum ThrottleGroup {
  MSG_APPLICATION       = 0,
  MSG_TIME_STEP         = 1,
  MSG_SOLVER            = 2
};

/**
 * @brief Class <b>Throttle</b> describes the cutoff limits for a message throttle.
 *
 */
struct Throttle
{
  /**
   * Creates a new <b>Throttle</b> instance.
   *
   * @param cutoff		a <b>size_t</b> value to display before the message is no longer
   *                            displayed.
   *
   * @param group		an <b>int</b> value to identify the throttle group that this message
   *                            belongs to.
   *
   */
  Throttle(size_t cutoff, int group)
    : m_cutoff(cutoff),
      m_group(group),
      m_count(0)
  {}

  size_t        m_cutoff;                       ///< Maximum number to display
  int           m_group;                        ///< Throttle group of message
  size_t        m_count;                        ///< Number which have been displayed
};

/**
 * @brief Class <b>MessageCode</b> declares a message identifier and throttle characteristics for a
 * message.  THESE MUST BE DECLARED STATIC.
 *
 * All messages have an associated message code.  This message code is used to identify a message
 * for throttling and aggregation.
 *
 * Message identifiers must be consist from reference to reference, unique for each instance, and
 * yet consist within each instance across multiple processors.  To meet these criteria, the message
 * identifier is implemented as a static memory location.  It must be declared by the application
 * developer as static.  This results in the linker selecting an address for each instance, which
 * never changes from reference to reference, is unique for each instance and the same regardless of
 * executable (assuming the executable is mapped into the same memory location for each process).
 * In order to remove the pointer-ness of the static memory location, the address is cast to a
 * pointer difference type which is an integral type by subtracting the zero pointer from it.
 *
 */
struct MessageCode
{
  /**
   * Creates a new <b>MessageCode</b> instance.
   *
   * @param throttle_cutoff	a <b>size_t</b> value to display before the message is no longer
   *                            displayed.
   *
   * @param throttle_group	an <b>int</b> value to identify the throttle group that this message
   *                            belongs to.
   *
   */
  MessageCode(size_t throttle_cutoff = 5, int throttle_group = MSG_APPLICATION)
    : m_id(&m_id - (MessageId *) 0),
      m_throttle(throttle_cutoff, throttle_group)
  {}

  /**
   * Creates a new <b>MessageCode</b> instance.  Be particularly careful when using this
   * constructor.  The message_id value must be the same on all processors for deferred message
   * reporting to work properly.
   *
   * @param message_id		a <b>MessageId</b> value of the message id.  This value must be the
   *                            same for each message across all processors.
   *
   * @param throttle_cutoff	a <b>size_t</b> value to display before the message is no longer
   *                            displayed.
   *
   * @param throttle_group	an <b>int</b> value to identify the throttle group that this message
   *                            belongs to.
   *
   */
  MessageCode(MessageId message_id, size_t throttle_cutoff, int throttle_group)
    : m_id(message_id),
      m_throttle(throttle_cutoff, throttle_group)
  {}

  static MessageCode    s_defaultMessageCode;   ///< Default message code

  MessageId             m_id;                   ///< Message identifier
  Throttle              m_throttle;             ///< Throttle characteristics
};

typedef std::ostream &(*OStreamFunctionPtr)(std::ostream &);
typedef std::ios_base &(*IOSBaseFunctionPtr)(std::ios_base &);

/**
 * @brief Enumeration <b>MessageType</b> declares the global message types.
 *
 * Currently warning and doomed (error) message types are defined.  Additional type may be added
 * after MSG_DOOMED.  The MSG_SYMMETRIC bit indicates that the message was generated identically
 * across all processors.
 *
 */
enum MessageType {
  MSG_TYPE_MASK         = 0x000000FF,           ///< Mask of levels
  MSG_WARNING           =          0,           ///< Message is a warning
  MSG_ERROR             =          1,           ///< Message is an error, but processing can proceed
  MSG_FATAL             =          2,           ///< Message is a fatal error, seg fault imminent
  MSG_EXCEPTION         =          3,           ///< Message is an exception
  MSG_INFORMATION       =          4,           ///< Message is informational
  MSG_DEBUG             =          5,           ///< Message is debug

  MSG_USER              = 0x00000100,
  MSG_DEVEL             = 0x00000200,

  MSG_SYMMETRIC         = 0x80000000,           ///< Message is symmetrical, same on all processors
  MSG_DEFERRED          = 0x40000000,           ///< Message is deferred, forward to 0, aggregated
  MSG_UNUSED0           = 0x20000000,
  MSG_UNUSED1           = 0x10000000,

  MSG_TERMINATE         = 0x00010000,

  // user error types:
  USR_FATAL        = MSG_USER | MSG_FATAL | MSG_TERMINATE,
  USR_ERROR        = MSG_USER | MSG_ERROR,
  USR_WARNING      = MSG_USER | MSG_WARNING,
  USR_FATAL_0      = USR_FATAL | MSG_SYMMETRIC,
  USR_ERROR_0      = USR_ERROR | MSG_SYMMETRIC,
  USR_WARNING_0    = USR_WARNING | MSG_SYMMETRIC,

  // developer error types:
  DEV_FATAL        = MSG_DEVEL | MSG_FATAL | MSG_TERMINATE,
  DEV_WARNING      = MSG_DEVEL | MSG_WARNING,
  DEV_FATAL_0      = DEV_FATAL | MSG_SYMMETRIC,
  DEV_WARNING_0    = DEV_WARNING | MSG_SYMMETRIC,
};

/**
 * @brief Function <b>operator<<</b> writes the message type name to the output stream.  If the
 * symmetric bit is set, "parallel" is prefixed to the name.
 *
 * @param os		a <b>std::ostream</b> reference to the output stream.
 *
 * @param message_type	a <b>MessageType</b> const reference to write the name of.
 *
 * @return		a <b>std::ostream</b> reference to os
 */
std::ostream &operator<<(std::ostream &os, const MessageType &message_type);

/**
 * Message creates an error message to be submitted to the error manager
 * on object destruction.
 *
 * Operator<< (put-to) can be used to populate a std::ostringstream
 * contained in the message class to make assembly of a meaningful error
 * meessage painless.
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 
 * @date   Tue Sep 17 09:28:54 2013
 */
class Message
{
public:
  Message(MessageType message_type, MessageCode message_code = MessageCode::s_defaultMessageCode);

  virtual ~Message();

private:
  Message(const Message &right);
  Message &operator=(const Message &);

public:
  Message &at(const NetlistLocation &netlist_location) 
  {
    netlistLocation_ = netlist_location;

    return *this;
  }

  Message &at(const std::string &path, int line_number) 
  {
    netlistLocation_ = NetlistLocation(path, line_number);

    return *this;
  }

  Message &die()
  {
    messageType_ |= MSG_TERMINATE;
    return *this;
  }

  std::ostringstream &os()
  {
    return oss_;
  }

  template<class T>
  Message &operator<<(const T &t)
  {   
    os() << t;
    return *this;
  }

  Message &operator<<(OStreamFunctionPtr f)
  {
    f(os());
    return *this;
  }

  Message &operator<<(IOSBaseFunctionPtr f)
  {
    f(os());
    return *this;
  }

private:
  unsigned            messageType_;
  MessageCode         messageCode_;
  std::ostringstream  oss_;
  NetlistLocation     netlistLocation_;

protected:
  const char *        functionName_;
};

struct UserFatal : public Message
{
  UserFatal(MessageCode message_code = MessageCode::s_defaultMessageCode)
    : Message(USR_FATAL, message_code)
  {}
};

struct UserFatal0 : public Message
{
  UserFatal0(MessageCode message_code = MessageCode::s_defaultMessageCode)
    : Message(USR_FATAL_0, message_code)
  {}
};

struct UserError : public Message
{
  UserError(MessageCode message_code = MessageCode::s_defaultMessageCode)
    : Message(USR_ERROR, message_code)
  {}
};

struct UserError0 : public Message
{
  UserError0(MessageCode message_code = MessageCode::s_defaultMessageCode)
    : Message(USR_ERROR_0, message_code)
  {}
};

struct UserWarning : public Message
{
  UserWarning(MessageCode message_code = MessageCode::s_defaultMessageCode)
    : Message(USR_WARNING, message_code)
  {}
};

struct UserWarning0 : public Message
{
  UserWarning0(MessageCode message_code = MessageCode::s_defaultMessageCode)
    : Message(USR_WARNING_0, message_code)
  {}
};

struct DevelFatal : public Message
{
  DevelFatal(MessageCode message_code = MessageCode::s_defaultMessageCode)
    : Message(DEV_FATAL, message_code)
  {}

  Message &in(const char *function_name)
  {
    functionName_ = functionName_;

    return *this;
  }
};

struct DevelFatal0 : public Message
{
  DevelFatal0(MessageCode message_code = MessageCode::s_defaultMessageCode)
    : Message(DEV_FATAL_0, message_code)
  {}

  Message &in(const char *function_name)
  {
    functionName_ = functionName_;

    return *this;
  }
};

struct DevelWarning : public Message
{
  DevelWarning(MessageCode message_code = MessageCode::s_defaultMessageCode)
    : Message(DEV_WARNING, message_code)
  {}

  Message &in(const char *function_name)
  {
    functionName_ = functionName_;

    return *this;
  }
};

struct DevelWarning0 : public Message
{
  DevelWarning0(MessageCode message_code = MessageCode::s_defaultMessageCode)
    : Message(DEV_WARNING_0, message_code)
  {}

  Message &in(const char *function_name)
  {
    functionName_ = functionName_;

    return *this;
  }
};

} // namespace Report
} // namespace Xyce

#endif // Xyce_N_ERH_Message_h

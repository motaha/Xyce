/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <N_UTL_ReportHandler.h>

#include <iostream>
#include <stdexcept>

namespace Xyce {

namespace {

REH s_reportHandler = &default_report_handler;

}

void
report(
  const char *		message,
  unsigned              type)
{
    (*s_reportHandler)(message, type);
}


void
default_report_handler(
  const char *		message,
  unsigned              type)
{
  std::cout << "Message type " << type << ": " << message << std::endl;
}


REH
set_report_handler(
  REH		        reh)
{
  /* %TRACE[ON]% */  /* %TRACE% */
  if (!reh)
    throw std::runtime_error("Cannot set report handler to NULL");

  REH prev_reh = s_reportHandler;
  s_reportHandler = reh;

  return prev_reh;
}

} // namespace Xyce

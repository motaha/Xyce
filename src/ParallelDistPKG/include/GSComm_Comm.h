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
// Filename       : $RCSfile: GSComm_Comm.h,v $
//
// Purpose        : 
//
// Special Notes  : 
//
// Creator        : 
//
// Creation Date  : 
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.7 $
//
// Revision Date  : $Date: 2014/02/24 23:49:25 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------
#ifndef _GSCOMM_COMM_H_
#define _GSCOMM_COMM_H_

//
// This plan is constructed for gather/scatter operations
// OO version derived from Zoltan's LB_comm_create

#include "GSComm_Plan.h"

class GSComm_Comm
{

public:

  // Default constructor
  GSComm_Comm();

  // Destructor
  ~GSComm_Comm();

  bool Do(GSComm_Plan & comm_plan, const int & tag, char * export_objs,
          const int & obj_size, char * import_objs);

  bool Do_Posts(GSComm_Plan & comm_plan, const int & tag, char * export_objs,
                const int & obj_size, char * import_objs);

  bool Do_Waits(GSComm_Plan & comm_plan, const int & tag, char * export_objs,
                const int & obj_size, char * import_objs);

  bool DoReverse(GSComm_Plan & comm_plan, const int & tag, char * export_objs,
                 const int & obj_size, char * import_objs);

  bool DoReverse_Posts(GSComm_Plan & comm_plan, const int & tag,
                       char * export_objs, const int & obj_size,
                       char * import_objs);

  bool DoReverse_Waits(GSComm_Plan & comm_plan, const int & tag,
                       char * export_objs, const int & obj_size,
                       char * import_objs);

private:

  char * recv_array_;
  char * send_array_;

  GSComm_Plan * comm_plan_reverse_;

};


// GSComm_Comm constructor
inline GSComm_Comm::GSComm_Comm()
  :
  recv_array_(0),
  send_array_(0),
  comm_plan_reverse_(0)
{
}

// GSComm_Comm destructor
inline GSComm_Comm::~GSComm_Comm()
{
  if (send_array_ != 0) delete[] send_array_;

  if (comm_plan_reverse_ != 0) delete comm_plan_reverse_;
}

// GSComm_Comm Do method
inline bool GSComm_Comm::Do(GSComm_Plan & comm_plan, const int & tag,
                            char * export_objs, const int & obj_size,
                            char * import_objs)
{
  bool comm_flag = true;

  Do_Posts(comm_plan, tag, export_objs, obj_size, import_objs);
  Do_Waits(comm_plan, tag, export_objs, obj_size, import_objs);

  return comm_flag;
}

// GSComm_Comm DoReverse method
inline bool GSComm_Comm::DoReverse(GSComm_Plan & comm_plan, const int & tag,
                                   char * export_objs, const int & obj_size,
                                   char * import_objs)
{
  bool comm_flag = true;

  DoReverse_Posts(comm_plan, tag, export_objs, obj_size, import_objs);
  DoReverse_Waits(comm_plan, tag, export_objs, obj_size, import_objs);

  return comm_flag;
}

#endif /* _GSCOMM_COMM_H_ */

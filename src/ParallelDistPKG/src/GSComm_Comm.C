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
// Filename       : $RCSfile: GSComm_Comm.C,v $
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
// Revision Number: $Revision: 1.8 $
//
// Revision Date  : $Date: 2014/02/24 23:49:25 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------


/* DEBUG:  missing standard header */

#include <Xyce_config.h>


#include "GSComm_Comm.h"

//---------------------------------------------------------------------------
//Do_Posts Method
//---------------------------------------------------------------------------
bool GSComm_Comm::Do_Posts( GSComm_Plan & comm_plan,
				const int & tag,
				char * export_objs,
				const int & obj_size,
				char * import_objs )
{
  int i, j, k;

  int my_proc = 0;
  int self_recv_address = 0;

  MPI_Comm_rank( comm_plan.comm_, &my_proc );

  recv_array_ = import_objs;

  j = 0;
  k = 0;

  for( i = 0; i < comm_plan.nrecvs_; ++i )
  {
    if( comm_plan.procs_from_[i] != my_proc )
    {
      MPI_Irecv( &(recv_array_[j]), comm_plan.lengths_from_[i] * obj_size,
        MPI_CHAR, comm_plan.procs_from_[i], tag, comm_plan.comm_,
        &(comm_plan.request_[k]) );
      ++k;
    }
    else
      self_recv_address = j;

    j += comm_plan.lengths_from_[i] * obj_size;
  }

  MPI_Barrier( comm_plan.comm_ );

  int self_num, self_index;

  if (comm_plan.max_send_length_* obj_size>0) 
    send_array_ = new char[ comm_plan.max_send_length_ * obj_size ];

  j = 0;
  for( i = 0; i < comm_plan.nsends_; ++i )
  {
    if( comm_plan.procs_to_[i] != my_proc )
    {
      int offset = 0;
      for( k = 0; k < comm_plan.lengths_to_[i]; ++k )
      {
        memcpy( &(send_array_[offset]), 
          &(export_objs[comm_plan.indices_to_[j]*obj_size]), obj_size );
        ++j;
        offset += obj_size;
      }
      //   cout << "my_proc = " << my_proc << " length = " << comm_plan.lengths_to_[i] * obj_size 
      //   << " send to = " << comm_plan.procs_to_[i] << " tag = " << tag << endl;
      MPI_Rsend( send_array_, comm_plan.lengths_to_[i] * obj_size,
	MPI_CHAR, comm_plan.procs_to_[i], tag, comm_plan.comm_ );
    }
    else
    {
      self_num = i;
      self_index = j;
      j += comm_plan.lengths_to_[i];
    }
  }

  if( comm_plan.self_msg_ )
    for( k = 0; k < comm_plan.lengths_to_[self_num]; ++k )
    {
      memcpy( &(recv_array_[self_recv_address]),
         &(export_objs[comm_plan.indices_to_[self_index]*obj_size]),
          obj_size );
      ++self_index;
      self_recv_address += obj_size;
    }

  if (send_array_!=0) {
    delete [] send_array_;
    send_array_ = 0;
  }
  return true;
}

//---------------------------------------------------------------------------
//Do_Waits Method
//---------------------------------------------------------------------------
bool GSComm_Comm::Do_Waits( GSComm_Plan & comm_plan,
			    const int & tag,
			    char * export_objs,
		 	    const int & obj_size,
			    char * import_objs )
{
   if( comm_plan.nrecvs_ > 0 )
    MPI_Waitall( comm_plan.nrecvs_ - comm_plan.self_msg_, 
		comm_plan.request_, comm_plan.status_ );

  return true;
}

//---------------------------------------------------------------------------
//DoReverse_Posts Method
//---------------------------------------------------------------------------
bool GSComm_Comm::DoReverse_Posts( GSComm_Plan & comm_plan,
				       const int & tag,
				       char * export_objs,
				       const int & obj_size,
				       char * import_objs )
{
  int i;
  int my_proc = 0;

  MPI_Comm_rank( comm_plan.comm_, &my_proc );

  int total_send_length = 0;
  for( i = 0; i < comm_plan.nsends_; ++i )
    total_send_length += comm_plan.lengths_to_[i];

  int max_recv_length = 0;
  for( i = 0; i < comm_plan.nrecvs_; ++i )
    if( comm_plan.procs_from_[i] != my_proc )
      if( comm_plan.lengths_from_[i] > max_recv_length )
        max_recv_length = comm_plan.lengths_from_[i];

  comm_plan_reverse_ = new GSComm_Plan;

  comm_plan_reverse_->lengths_to_ = comm_plan.lengths_from_;
  comm_plan_reverse_->procs_to_ = comm_plan.procs_from_;
  comm_plan_reverse_->indices_to_ = comm_plan.indices_from_;
  comm_plan_reverse_->lengths_from_ = comm_plan.lengths_to_;
  comm_plan_reverse_->procs_from_ = comm_plan.procs_to_;
  comm_plan_reverse_->indices_from_ = comm_plan.indices_to_;
  comm_plan_reverse_->nrecvs_ = comm_plan.nsends_;
  comm_plan_reverse_->nsends_ = comm_plan.nrecvs_;
  comm_plan_reverse_->self_msg_ = comm_plan.self_msg_;
  comm_plan_reverse_->max_send_length_ = max_recv_length;
  comm_plan_reverse_->total_recv_length_ = total_send_length;
  comm_plan_reverse_->comm_ = comm_plan.comm_;

  comm_plan_reverse_->request_ = new MPI_Request[ comm_plan_reverse_->nrecvs_ ];
  comm_plan_reverse_->status_= new MPI_Status[ comm_plan_reverse_->nrecvs_ ];

  comm_plan_reverse_->no_delete_ = true;

  bool comm_flag = Do_Posts( *comm_plan_reverse_, tag, 
	export_objs, obj_size, import_objs );

  return comm_flag;
}

//---------------------------------------------------------------------------
//DoReverse_Waits Method
//---------------------------------------------------------------------------
bool GSComm_Comm::DoReverse_Waits( GSComm_Plan & comm_plan,
			           const int & tag,
			           char * export_objs,
		 	           const int & obj_size,
			           char * import_objs )
{
  if( comm_plan_reverse_ == 0 ) return false;

  bool comm_flag = Do_Waits( *comm_plan_reverse_, tag, export_objs, obj_size,
	import_objs );

  if (comm_plan_reverse_!=0) {
    delete comm_plan_reverse_;
    comm_plan_reverse_ = 0;
  }

  return comm_flag;
}


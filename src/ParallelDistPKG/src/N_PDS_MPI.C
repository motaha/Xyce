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
// Filename       : $RCSfile: N_PDS_MPI.C,v $
//
// Purpose        : 
//                  
//                  
//
// Special Notes  : 
//                  
//
// Creator        : David Baur
//
// Creation Date  : 
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.4.2.1 $
//
// Revision Date  : $Date: 2014/02/27 00:52:18 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifdef Xyce_PARALLEL_MPI

#include <sstream>

#include <N_PDS_MPI.h>

namespace Xyce {
namespace Parallel {

template struct Loc<int>;
template struct Loc<double>;
template struct Loc<float>;

MPI_Datatype
double_complex_type()
{
  static MPI_Datatype s_mpi_double_complex;
  static bool initialized = false;

  if (!initialized) {
    initialized = true;

    MPI_Type_contiguous(2, MPI_DOUBLE, &s_mpi_double_complex);
    MPI_Type_commit(&s_mpi_double_complex);
  }
  return s_mpi_double_complex;
}

MPI_Datatype
float_complex_type()
{
  static MPI_Datatype s_mpi_float_complex;
  static bool initialized = false;

  if (!initialized) {
    initialized = true;

    MPI_Type_contiguous(2, MPI_FLOAT, &s_mpi_float_complex);
    MPI_Type_commit(&s_mpi_float_complex);
  }
  return s_mpi_float_complex;
}


// #ifdef MPI_LONG_LONG_INT
// MPI_Datatype
// long_long_int_int_type()
// {
//   static MPI_Datatype s_mpi_long_long_int_int;
//   static bool initialized = false;

//   int B[] = {2, 1};
//   MPI_Aint D[] = {0, 8};
//   MPI_Datatype T[] = {MPI_LONG_LONG_INT, MPI_INT};
  
//   if (!initialized) {
//     initialized = true;

//     MPI_Type_struct(2, B, D, T, &s_mpi_long_long_int_int);
//     MPI_Type_commit(&s_mpi_long_long_int_int);
//   }
//   return s_mpi_long_long_int_int;
// }
// #endif


MPI_Datatype
double_double_int_type()
{
  static MPI_Datatype s_mpi_double_double_int;
  static bool initialized = false;

  int B[] = {2, 1};
  MPI_Aint D[] = {0, 16};
  MPI_Datatype T[] = {MPI_DOUBLE, MPI_INT};
  
  
  if (!initialized) {
    initialized = true;

    MPI_Type_struct(2, B, D, T, &s_mpi_double_double_int);
    MPI_Type_commit(&s_mpi_double_double_int);
  }
  return s_mpi_double_double_int;
}


namespace {

extern "C" {
  void
  mpi_double_complex_sum(
    void *		invec,
    void *		inoutvec,
    int *		len,
    MPI_Datatype *	datatype)
  {
    std::complex<double> *complex_in = static_cast<std::complex<double> *>(invec);
    std::complex<double> *complex_inout = static_cast<std::complex<double> *>(inoutvec);

    for (int i = 0; i < *len; ++i)
      complex_inout[i] += complex_in[i];
  }
} // extern "C"

} // namespace <unnamed>


MPI_Op
double_complex_sum_op()
{
  static MPI_Op s_mpi_double_complex_sum;
  static bool initialized = false;

  if (!initialized) {
    initialized = true;

    MPI_Op_create(mpi_double_complex_sum, true, &s_mpi_double_complex_sum);
  }
  return s_mpi_double_complex_sum;
}

namespace {

const Parallel::ReduceSet *s_currentReduceSet = 0;

extern "C" {
  typedef void (*ParallelReduceOp)
  (void * inv, void * outv, int *, MPI_Datatype *);
}

void
all_reduce(
  MPI_Comm		arg_comm,
  ParallelReduceOp	arg_op,
  void *		arg_in,
  void *		arg_out,
  unsigned		arg_len)
{
  MPI_Op mpi_op = MPI_OP_NULL ;

  MPI_Op_create(arg_op, 0, & mpi_op);

  // The SUN was buggy when combining an
  // MPI_Allreduce with a user defined operator,
  // use reduce/broadcast instead.

  const int result = MPI_Allreduce(arg_in,arg_out,arg_len,MPI_BYTE,mpi_op,arg_comm);

  MPI_Op_free(& mpi_op);

  if (MPI_SUCCESS != result) {
    std::ostringstream msg ;
    msg << "Xyce::MPI::all_reduce FAILED: MPI_Allreduce = " << result;
    throw std::runtime_error(msg.str());
  }
}

struct ReduceCheck : public ReduceInterface
{
  ReduceCheck()
  {}

  void setSize(unsigned size) {
    m_size = size;
  }

  virtual void size(void *&inbuf) const {
    unsigned *t = align_cast<unsigned>(inbuf);
    t += sizeof(unsigned);
    inbuf = t;
  }

  virtual void copyin(void *&inbuf) const {
    unsigned *t = align_cast<unsigned>(inbuf);
    *t++ = m_size;
    inbuf = t;
  }

  virtual void copyout(void *&outbuf) const {
    unsigned *t = align_cast<unsigned>(outbuf);

    unsigned size = *t++;
    if (m_size != size)
      throw std::runtime_error("size mismatch");

    outbuf = t;
  }

  virtual void op(void *&inbuf, void *&outbuf) const {
    unsigned *tin = align_cast<unsigned>(inbuf);
    unsigned *tout = align_cast<unsigned>(outbuf);

    *tout = std::min(*tout, *tin);

    inbuf = ++tin;
    outbuf = ++tout;
  }

private:
  unsigned	m_size;
};

} // namespace <unnamed>


ReduceSet::ReduceSet()
{
  add(new ReduceCheck);
}


ReduceSet::~ReduceSet()
{
  for (ReduceVector::const_iterator it = m_reduceVector.begin(); it != m_reduceVector.end(); ++it)
    delete (*it);
}


size_t
ReduceSet::size() const {
  void *buffer_end = 0;

  for (ReduceVector::const_iterator it = m_reduceVector.begin(); it != m_reduceVector.end(); ++it)
    (*it)->size(buffer_end);

  ReduceCheck *reduce_check = static_cast<ReduceCheck *>(m_reduceVector.front());
  reduce_check->setSize(reinterpret_cast<char *>(buffer_end) - (char *) 0);

  return reinterpret_cast<char *>(buffer_end) - (char *) 0;
}

void
ReduceSet::copyin(void * const buffer_in) const {
  void *inbuf = buffer_in;

  for (ReduceVector::const_iterator it = m_reduceVector.begin(); it != m_reduceVector.end(); ++it)
    (*it)->copyin(inbuf);
}

void
ReduceSet::copyout(void * const buffer_out) const {
  void *outbuf = buffer_out;

  for (ReduceVector::const_iterator it = m_reduceVector.begin(); it != m_reduceVector.end(); ++it)
    (*it)->copyout(outbuf);
}

void
ReduceSet::op(void * const buffer_in, void * const buffer_out) const {
  void *inbuf = buffer_in;
  void *outbuf = buffer_out;

  for (ReduceVector::const_iterator it = m_reduceVector.begin(); it != m_reduceVector.end(); ++it)
    (*it)->op(inbuf, outbuf);
}

void ReduceSet::void_op(void * inv, void * outv, int *, MPI_Datatype *) {
  s_currentReduceSet->op(inv, outv);
}


void
ReduceSet::add(
  ReduceInterface *	reduce_interface)
{
  m_reduceVector.push_back(reduce_interface);
}


void
AllReduce(
  MPI_Comm		comm,
  const ReduceSet &	reduce_set)
{
  size_t size = reduce_set.size();

  if (size) {
    char *input_buffer  = new char[size];
    char *output_buffer = new char[size];
    void *inbuf = (void *) input_buffer;
    void *outbuf = (void *) output_buffer;

    s_currentReduceSet = &reduce_set;

    ParallelReduceOp f = reinterpret_cast<ParallelReduceOp>(& ReduceSet::void_op);

    reduce_set.copyin(inbuf);
    all_reduce(comm, f, inbuf, outbuf, size);
    reduce_set.copyout(outbuf);
    delete [] output_buffer;
    delete [] input_buffer;
  }
}

void
AllWriteString(
  MPI_Comm              mpi_comm,
  std::ostream &        os,
  const std::string &   message)
{
  if (mpi_parallel_run(mpi_comm)) {
    const int root = 0;
    const unsigned size = Parallel::size(mpi_comm);
    const unsigned rank = Parallel::rank(mpi_comm);

    int result;

    // Gather the send counts on root processor
    int send_count = message.size();

    std::vector<int> recv_count(size, 0);
    int * const recv_count_ptr = &recv_count[0];

    result = MPI_Gather(& send_count, 1, MPI_INT,
                        recv_count_ptr, 1, MPI_INT,
                        root, mpi_comm);

    if (MPI_SUCCESS != result) {
      std::ostringstream message;
      message << "stk::all_write FAILED: MPI_Gather = " << result;
      throw std::runtime_error(message.str());
    }

    // Receive counts are only non-zero on the root processor:
    std::vector<int> recv_displ(size + 1, 0);
    for (unsigned i = 0; i < size; ++i) {
      recv_displ[i+1] = recv_displ[i] + recv_count[i];
    }

    const unsigned recv_size = (unsigned) recv_displ[ size ];
    std::vector<char> buffer(recv_size);

    {
      const char * const send_ptr = message.c_str();
      char * const recv_ptr = recv_size ? & buffer[0] : (char *) NULL;
      int * const recv_displ_ptr = & recv_displ[0];

      result = MPI_Gatherv((void*) send_ptr, send_count, MPI_CHAR,
                           recv_ptr, recv_count_ptr, recv_displ_ptr, MPI_CHAR,
                           root, mpi_comm);
    }

    if (MPI_SUCCESS != result) {
      std::ostringstream message;
      message << "stk::all_write FAILED: MPI_Gatherv = " << result;
      throw std::runtime_error(message.str());
    }

    if (root == (int) rank) {
      for (unsigned i = 0; i < size; ++i) {
        if (recv_count[i]) {
          char * const ptr = & buffer[ recv_displ[i] ];
          os.write(ptr, recv_count[i]);
          os << std::endl;
        }
      }
      os.flush();
    }
  }
  else
    os << message;
}

} // namespace Parallel
} // namespace Xyced

#endif // Xyce_PARALLEL_MPI

#ifndef XYCE_UTIL_Serial_h
#define XYCE_UTIL_Serial_h

#include <vector>

#include <N_PDS_fwd.h>
#include <N_PDS_ParallelMachine.h>

#ifndef Xyce_PARALLEL_MPI

static const int MPI_SUM = 0;
static const int MPI_MAX = 0;
static const int MPI_MIN = 0;

namespace Xyce {
namespace Parallel {

template<class T>
inline void
AllReduce(Machine comm, int op, T *src_dest, size_t size)
{}

template<class T>
inline void
AllReduce(Machine comm, int op, const T *source, T *dest, size_t size)
{
  std::copy(source, source + size, dest);
}

template<class T>
inline void
AllReduce(Machine comm, int op, std::vector<T> &src_dest)
{}

template<class T>
inline void
Broadcast(Machine comm, T *src_dest, size_t len, int root)
{}

inline void
AllWriteString(
  Machine               comm,
  std::ostream &        os,
  const std::string &   message)
{
  os << message;
}

inline void Barrier(Machine com)
{}

} // namespace Parallel
} // namespace Xyce

#endif // Xyce_PARALLEL_MPI

#endif // Xyce_UTIL_Serial_h

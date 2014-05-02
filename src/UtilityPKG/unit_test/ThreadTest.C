// g++ ThreadTest.C -I../include -I ../../ParallelDistPKG/include -lpthread

#include <Xyce_config.h>

#include <iosfwd>
#include <fstream>
#include <sstream>
#include <vector>

#include <N_UTL_PThread.h>
#include <N_UTL_LogStream.h>

typedef std::vector<Xyce::Util::xyce_pthread_t> ThreadVector;

struct ThreadData
{
  ThreadData(int r, int s)
    : rank(r),
      size(s)
  {}
  
  int           rank;
  int           size;
};


void *run_thread(void *arg)
{
  ThreadData *thread_data = static_cast<ThreadData *>(arg);
  int rank = thread_data->rank;
  int size = thread_data->size;

  Xyce::Util::xyce_pthread_t self = Xyce::Util::xyce_pthread_self();

  std::ostringstream path;
  path << "log." << rank << "." << size;

  std::ofstream log_os(path.str().c_str());
  log_os << "Logging started on thread rank " << rank << " of " << size << " on thread " << (void *)(self) << std::endl;
  for (int i = 0; i < 1000; ++i) {
    if (i%10 == 2)
      Xyce::addThreadStream(&log_os);

    if (i%10 == 9)
      Xyce::removeThreadStream(&log_os);

    Xyce::lout() << "This is step " << i << " on thread rank " << rank << " of " << size << " on thread " << self << std::endl;
  }

  Xyce::removeThreadStream(&log_os);
  log_os << "Logging complete on thread rank " << rank << " of " << size << " on thread " << self << std::endl;
  log_os.close();

  return 0;
}


int
main(int argc, char **argv)
{
  Xyce::initializeLogStreamByThread();

  int thread_count = 4;
  if (argc > 1) {
    std::istringstream is(argv[1]);
    is >> thread_count;
  }

  std::vector<Xyce::Util::xyce_pthread_t> thread_list(thread_count);

  for (ThreadVector::iterator it = thread_list.begin(); it != thread_list.end(); ++it) {
    ThreadData *thread_data = new ThreadData(it - thread_list.begin(), thread_list.size());

    Xyce::Util::xyce_pthread_create(&*it, 0, run_thread, thread_data);
  }
  
  for (ThreadVector::iterator it = thread_list.begin(); it != thread_list.end(); ++it)
    Xyce::Util::xyce_pthread_join(*it, 0);
}

// Hack since utilite really can't depend on ParallelDist and the test lives in utility
#ifdef Xyce_PARALLEL_MPI
namespace Xyce {
namespace Parallel {

void AllWriteString(MPI_Comm comm, std::ostream &os, const std::string &message) {
  os << message;
}

}
}
#endif // Xyce_PARALLEL_MPI

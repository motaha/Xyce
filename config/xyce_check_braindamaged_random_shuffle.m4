dnl this is a test for a serious kludge around some linuxen.
dnl random_shuffle is used in N_IO_DistribMgr.C, but at least with some 
dnl versions of the GNU STL and gnu libc, -ansi causes an inconsistent
dnl set of -D flags and causes N_IO_DistribMgr.C not to compile when 
dnl compiling for parallel.  Let's see what we can do about that.

AC_DEFUN([XYCE_CHECK_BRAINDAMAGED_RANDOM_SHUFFLE],
[AC_CACHE_CHECK(whether your gcc has a braindamaged STL algorithm random_shuffle,
ac_cv_xyce_braindamaged_random_shuffle,
[AC_LANG_PUSH(C++)
 AC_COMPILE_IFELSE(
[AC_LANG_PROGRAM([[
#ifdef HAVE_CSTDLIB
#include <cstdlib>
#else
#include <stdlib.h>
#endif

#include <vector>
#ifdef HAVE_ALGORITHM
#include <algorithm>
#else
#include <algo.h>
#endif
]],
[[

  std::vector<int> ipVec(5,0);

  std::random_shuffle(ipVec.begin(),ipVec.end());
]])], ac_cv_xyce_braindamaged_random_shuffle=no , ac_cv_xyce_braindamaged_random_shuffle=yes)
AC_LANG_POP(C++)
])
])

AC_DEFUN([XYCE_CHECK_IOTA],
[AC_CACHE_CHECK(whether your STL algorithm includes iota,
ac_cv_xyce_have_iota,
[AC_LANG_PUSH(C++)
 AC_COMPILE_IFELSE(
[AC_LANG_PROGRAM([[#ifdef HAVE_CSTDLIB
#include <cstdlib>
#else
#include <stdlib.h>
#endif

#include <vector>
#include <numeric>
#ifdef HAVE_ALGORITHM
#include <algorithm>
#else
#include <algo.h>
#endif
]],
[[
  std::vector<int> ipVec(5,0);

  std::iota(ipVec.begin(),ipVec.end(),1);
]])], ac_cv_xyce_have_iota=yes , ac_cv_xyce_have_iota=no)
AC_LANG_POP(C++)
])
if test "$ac_cv_xyce_have_iota" = yes; then
  AC_DEFINE(HAVE_IOTA,1,[define if the Standard Template Library algorithms include iota])
fi
])

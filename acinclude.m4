##########Does the stuff I wrote for Xyce work?
######################################################################
# These macros from the GNU Autoconf Macro Archive at 
# http://www.gnu.org/software/ac-archive/
######################################################################
AC_DEFUN([AC_CXX_NAMESPACES],
[AC_CACHE_CHECK(whether the compiler implements namespaces,
ac_cv_cxx_namespaces,
[AC_LANG_PUSH(C++)
 AC_COMPILE_IFELSE(
[AC_LANG_PROGRAM([[namespace Outer { namespace Inner { int i = 0; }}]],
                [[using namespace Outer::Inner; return i;]])],
 ac_cv_cxx_namespaces=yes, ac_cv_cxx_namespaces=no)
 AC_LANG_POP(C++)
])
if test "$ac_cv_cxx_namespaces" = yes; then
  AC_DEFINE(HAVE_NAMESPACES,,[define if the compiler implements namespaces])
fi
])

AC_DEFUN([AC_CXX_HAVE_NUMERIC_LIMITS],
[AC_CACHE_CHECK(whether the compiler has numeric_limits<T>,
ac_cv_cxx_have_numeric_limits,
[AC_REQUIRE([AC_CXX_NAMESPACES])
 AC_LANG_PUSH(C++)
 AC_COMPILE_IFELSE(
[AC_LANG_PROGRAM([[#include <limits>
]],[[double e = std::numeric_limits<double>::epsilon(); return 0;]])],
 ac_cv_cxx_have_numeric_limits=yes, ac_cv_cxx_have_numeric_limits=no)
 AC_LANG_POP(C++)
])
if test "$ac_cv_cxx_have_numeric_limits" = yes; then
  AC_DEFINE(HAVE_NUMERIC_LIMITS,1,[define if the compiler has numeric_limits<T>])
fi
])

########################################################################
# We use sprintf sometimes.  This requires stdio somehow.  Some C++ allow
# #include <cstdio>, where that is possible use it, otherwise
# #include <stdio.h>
########################################################################
AC_DEFUN([AC_CXX_HAVE_CSTDIO],
[AC_CACHE_CHECK(whether the compiler has cstdio,
ac_cv_cxx_have_cstdio,
[AC_REQUIRE([AC_CXX_NAMESPACES])
 AC_LANG_PUSH(C++)
 AC_COMPILE_IFELSE(
[AC_LANG_PROGRAM([[#include <cstdio>
]],[[char foobie[128];std::sprintf(foobie,"%d\n",1); return 0;]])],
 ac_cv_cxx_have_cstdio=yes, ac_cv_cxx_have_cstdio=no)
 AC_LANG_POP(C++)
])
if test "$ac_cv_cxx_have_cstdio" = yes; then
  AC_DEFINE(HAVE_CSTDIO,1,[define if the compiler has cstdio include])
fi
])

########################################################################
# Exception handling for floating point exceptions can be enabled on
# linux.  As this is a common development and deployment platform, such
# handling on this platform should prevent this type of error from
# creeping into the code.
########################################################################
AC_DEFUN([AC_CXX_HAVE_FENV],
[AC_CACHE_CHECK(whether the compiler has fenv.h,
ac_cv_cxx_have_fenv_h,
[AC_REQUIRE([AC_CXX_NAMESPACES])
 AC_LANG_PUSH(C++)
 AC_COMPILE_IFELSE(
[AC_LANG_PROGRAM([[#include <fenv.h>
]],[[std::feenableexcept(FE_DIVBYZERO); return 0;]])],
 ac_cv_cxx_have_fenv_h=yes, ac_cv_cxx_have_fenv_h=no)
 AC_LANG_POP(C++)
])
if test "$ac_cv_cxx_have_fenv_h" = yes; then
  AC_DEFINE(HAVE_LINUX_EXCEPTIONS, 1,[define if the compiler has fenv include])
fi
])

########################################################################
#STL Algorithms
########################################################################
#we need to know whether to use algortihm or algo.h
########################################################################
AC_DEFUN([AC_CXX_HAVE_ALGORITHM],
[AC_CACHE_CHECK(whether the compiler has algorithm,
ac_cv_cxx_have_algorithm,
[AC_REQUIRE([AC_CXX_NAMESPACES])
 AC_REQUIRE([AC_CXX_HAVE_STL])
 AC_LANG_PUSH(C++)
 AC_COMPILE_IFELSE(
[AC_LANG_PROGRAM([[#include <algorithm>
]],[[int baz; baz = std::min(1,2); return 0;]])],
 ac_cv_cxx_have_algorithm=yes, ac_cv_cxx_have_algorithm=no)
 AC_LANG_POP(C++)
])
if test "$ac_cv_cxx_have_algorithm" = yes; then
  AC_DEFINE(HAVE_ALGORITHM,1,[define if the compiler has algorithm])
fi
])
########################################################################
AC_DEFUN([AC_CXX_HAVE_ALGO_H],
[AC_CACHE_CHECK(whether the compiler has algo.h,
ac_cv_cxx_have_algo_h,
[AC_REQUIRE([AC_CXX_NAMESPACES])
 AC_REQUIRE([AC_CXX_HAVE_STL])
 AC_LANG_PUSH(C++)
 AC_COMPILE_IFELSE(
[AC_LANG_PROGRAM([[#include <algo.h>
]],[[int baz; baz = std::min(1,2); return 0;]])],
 ac_cv_cxx_have_algo_h=yes, ac_cv_cxx_have_algo_h=no)
 AC_LANG_POP(C++)
])
if test "$ac_cv_cxx_have_algo_h" = yes; then
  AC_DEFINE(HAVE_ALGO_H,1,[define if the compiler has algo.h])
fi
])

########################################################################
# Check if the STL implementation requires that pair.h be included
# Note that we presuppose that string exists... probably a problem
# In fact, this is a kludgy test... really what we test is if the pair
# works without pair.h included, then assume that including pair.h
# fixes it.  Think about this a little harder
########################################################################
AC_DEFUN([AC_CXX_NEED_PAIR_H],
[AC_CACHE_CHECK(whether use of stl strings and pairs requires including pair.h,
ac_cv_cxx_need_pair_h,
[ AC_REQUIRE([AC_CXX_NAMESPACES])
 AC_REQUIRE([AC_CXX_HAVE_STL])
 AC_LANG_PUSH(C++)
 AC_COMPILE_IFELSE(
[AC_LANG_PROGRAM([[ #include <string>
]],[[typedef std::pair<std::string,int> NodeID]])],
 ac_cv_cxx_need_pair_h=no, ac_cv_cxx_need_pair_h=yes)
 AC_LANG_POP(C++)
] )
if test "$ac_cv_cxx_need_pair_h" = yes; then
  AC_DEFINE(NEED_PAIR_H,1,[define if stl string header does not cause pair.h to be included])
fi
])

######################################################################
# we use M_PI somewhere, need to define it.  Check if math.h defines it
######################################################################
AC_DEFUN([AC_HEADER_MATH_H_HAS_M_PI],
[AC_CACHE_CHECK(whether math.h defines M_PI,
 ac_cv_header_math_h_has_m_pi,
[AC_LANG_PUSH(C)
 AC_COMPILE_IFELSE(
[AC_LANG_PROGRAM([[
#include <math.h>
]],
[[
#ifndef M_PI
error not defined
#endif
]])],ac_cv_header_math_h_has_m_pi=yes,ac_cv_header_math_h_has_m_pi=no)
AC_LANG_POP(C)
] )
if test "$ac_cv_header_math_h_has_m_pi" = yes; then
  AC_DEFINE(MATH_H_HAS_M_PI,1,[define if math.h defines M_PI])
fi
])

AC_DEFUN([AC_CXX_HAVE_STL],
[AC_CACHE_CHECK(whether the compiler supports Standard Template Library,
ac_cv_cxx_have_stl,
[AC_REQUIRE([AC_CXX_NAMESPACES])
 AC_LANG_PUSH(C++)
 AC_COMPILE_IFELSE(
[AC_LANG_PROGRAM([[#include <list>
#include <deque>
]],[[std::list<int> x; x.push_back(5);
std::list<int>::iterator iter = x.begin(); if (iter != x.end()) ++iter; return 0;]])],
 ac_cv_cxx_have_stl=yes, ac_cv_cxx_have_stl=no)
 AC_LANG_POP(C++)
])
if test "$ac_cv_cxx_have_stl" = yes; then
  AC_DEFINE(HAVE_STL,,[define if the compiler supports Standard Template Library])
fi
])


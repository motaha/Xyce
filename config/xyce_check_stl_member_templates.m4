AC_DEFUN([XYCE_CHECK_STL_MEMBER_TEMPLATES],
[AC_CACHE_CHECK(whether your STL supports member templates,
ac_cv_xyce_stl_member_templates,
[AC_LANG_PUSH(C++)
 AC_COMPILE_IFELSE(
[AC_LANG_PROGRAM([[#include <vector>
#include <list>
#ifdef HAVE_NAMESPACES
using namespace std;
#endif
]],
[[
	vector<int> v(2);
	v[0] = v[1] = 0;
	list<char> l( v.begin(), v.end() );
]])], ac_cv_xyce_stl_member_templates=yes , ac_cv_xyce_stl_member_templates=no)
AC_LANG_POP(C++)
])
if test "$ac_cv_xyce_stl_member_templates" = no; then
  AC_DEFINE(BAD_STL,1,[define if the STL cannot support member templates])
fi
])

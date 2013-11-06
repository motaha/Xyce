#The DEC compiler uses a braindead way of handling template instantiation.
# check if we're using it. 
AC_DEFUN([XYCE_TEST_USING_DEC_CXX],
[AC_CACHE_CHECK(whether your C++ compiler is the DEC compiler,
ac_cv_cxx_xyce_using_dec_cxx,
[AC_LANG_PUSH(C++)
 AC_COMPILE_IFELSE(
[AC_LANG_PROGRAM([[]],[[
#ifdef __DECCXX
#define A 1
#else
error not a dec compiler
#endif
int i = 1;
]])],ac_cv_cxx_xyce_using_dec_cxx=yes, ac_cv_cxx_xyce_using_dec_cxx=no)
AC_LANG_POP(C++)
])
if test "$ac_cv_cxx_xyce_using_dec_cxx" = yes; then
  USING_DEC_CXX=yes
fi
])

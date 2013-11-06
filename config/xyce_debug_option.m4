dnl   XYCE_DEBUG_OPTION(chaco,no,[chaco],Xyce_USE_CHACO,USE_CHACO)
dnl  will add -DXyce_USE_CHACO to CPPFLAGS and define USE_CHACO=yes if the
dnl  user 
AC_DEFUN([XYCE_DEBUG_OPTION],
[dnl
dnl get the arguments
m4_pushdef([xyce_name], [$1]) dnl
m4_pushdef([xyce_name_upcase], [m4_toupper(xyce_name)]) dnl
dnl
dnl default
m4_pushdef([xyce_default], [$2]) dnl
dnl
dnl symbol
m4_pushdef([xyce_symbol], m4_default($4, [Xyce_]xyce_name_upcase))
dnl
AC_ARG_ENABLE([$1], [AS_HELP_STRING([--enable-$1], [enable $3.])], [],
              [enable_[]xyce_name[]=xyce_default; enableval=$enable_[]xyce_name[]])
dnl AC_ARG_WITH($1,[  --with-$1                     DOES ABSOLUTELY NOTHING!.],[
dnl AC_MSG_ERROR([invalid option --with-$1.  Perhaps you meant --enable-$1.])])
dnl
if test "x$enable_[]xyce_name[]" != "xno"; then
   CPPFLAGS="-D[]xyce_symbol[] $CPPFLAGS"
fi
dnl
m4_ifval($5, [$5=$enable_[]xyce_name[]])
dnl
m4_popdef([xyce_name], [xyce_name_upcase], [xyce_default]) dnl
m4_popdef([xyce_symbol]) dnl
])

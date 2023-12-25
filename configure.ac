#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.71])
AC_INIT([simulWaterHoseWall],[1.0.2],[dominiqu-jc@users.sourceforge.net])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([simulWaterHoseWall/main.c])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC

# Checks for libraries.
AC_CHECK_LIB([m], [sqrt])

# Checks for header files.
AC_CHECK_HEADERS([inttypes.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_TYPE_UINT32_T

# Checks for library functions.
AC_CHECK_FUNCS([sqrt strrchr])

AC_CONFIG_FILES([Makefile
                 doc/Makefile
                 examples/Makefile
                 simulWaterHoseWall/Makefile])
AC_OUTPUT

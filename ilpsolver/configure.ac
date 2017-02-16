#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT(ilpsolver, 1.0, mingfu.shao@cs.cmu.edu)
AM_INIT_AUTOMAKE()
AC_CONFIG_SRCDIR([src/main.cc])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries.
AS_IF([test "s$GUROBI_HOME" = "s"], [AC_MSG_ERROR([*** Please set environment variable GUROBI_HOME ***])]) 
AS_IF([test "s$GRB_LICENSE_FILE" = "s"], [AC_MSG_ERROR([*** Please set environment variable GRB_LICENSE_FILE ***])]) 

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT

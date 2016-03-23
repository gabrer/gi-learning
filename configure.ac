dnl Process this file with autoconf to produce a configure script.

# autoconf minimum version
AC_PREREQ(2.59)


AC_INIT(GI-learning, 1.0, [pietro.cottone@gmail.com], [GI-learning])

#unique-file-in-source-dir is some file that is in the package's source directory; configure checks for this file's existence to make sure that the directory that it is told contains the source code in fact does. 
#Occasionally people accidentally specify the wrong directory with --srcdir; this is a safety check. See configure Invocation, for more information.
AC_CONFIG_SRCDIR([src/GI-learning.cpp])

# make config header
AC_CONFIG_HEADERS([config.h])

AM_MAINTAINER_MODE

# auxiliary files
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])

#obsolete
#AC_CANONICAL_SYSTEM

# init automake
AM_INIT_AUTOMAKE([1.11 subdir-objects -Wall -Werror])
AM_PROG_AR

#da Gabriele per .so
#LT_INIT

# "configure" creates Makefile
#AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES(Makefile src/Makefile)

#openMP
AC_OPENMP
AC_SUBST(AM_CXXFLAGS,"$OPENMP_CXXFLAGS -fopenmp")
#AX_OPENMP


## Checks for libraries.
#have_pthreads=no
#AC_SEARCH_LIBS([pthread_create], [pthread], [have_pthreads=yes])
## Checks for header files.
#AC_CHECK_HEADERS([stdlib.h])
#if test "x${have_pthreads}" = xyes; then
#AC_CHECK_HEADERS([pthread.h], [], [have_pthreads=no])
#fi
#if test "x${have_pthreads}" = xno; then
#AC_MSG_ERROR([
#------------------------------------------
#The pthread library and header file
#required to build jupiter. Stopping...
#Check 'config.log' for more information.
#------------------------------------------])
#fi

# LIBTOOL
LT_INIT
AC_DISABLE_SHARED

# find and probe C++ compiler and 
AC_PROG_CXX


# ------Check for boost_filesystem --------
BOOST_REQUIRE([1.54])
BOOST_FILESYSTEM
BOOST_SYSTEM
#AC_CHECK_LIB(boost_filesystem, main,, [
#	AC_MSG_ERROR([You need the C++ boost filesystem library.])])
#-------------------------------------------

# ----------Check for libxml++ -------------
#PKG_CHECK_MODULES(SOMEAPP, libxml++-2.6 >= 2.38.0)
#    #AC_SUBST(SOMEAPP_CFLAGS)
#    #AC_SUBST(SOMEAPP_LIBS)
#	
#	PKG_CHECK_MODULES(XML,
#   [libxml++-2.6 >= 2.38.0],
#  ,
#    [
#      [echo "SDL_MIXER not found via pkg_config, checking with AC_CHECK_LIB:"]
#      AC_CHECK_LIB(libxml++-2.6,
#	 main,
#         ,
#         [AC_MSG_ERROR([SDL_mixer not found!
#http://www.libsdl.org/projects/SDL_mixer])]
#      )
#      [echo "SDL_MIXER successfully located"]
#    ]
#)
#
#AC_DEFINE([HAVE_LIBSDL_MIXER],[1],[Define to 1 if you have the
#`SDL_mixer` library])
#CFLAGS="$CFLAGS $SDL_MIXER_CFLAGS"
#LIBS="$LIBS $SDL_MIXER_LIBS"
PKG_CHECK_MODULES([LIBXMLXX], [libxml++-2.6])
#------------------------------------------------
	
	




# default compiler options
: ${CXXFLAGS=""}



#make this data visible in Makefile.am
#AC_SUBST(DEPS_CFLAGS)
#AC_SUBST(DEPS_LIBS)

# create output
AC_OUTPUT


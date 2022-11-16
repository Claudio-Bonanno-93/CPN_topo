AC_INIT([CPN_topo], [1.0.0], [claudiobonanno93@gmail.com])
AM_INIT_AUTOMAKE(subdir-objects)
AC_PROG_CC
AC_CONFIG_SRCDIR([./include/macro.h])

# overwrite default values
if test x${ac_cv_env_CFLAGS_set} != xset; then
  CFLAGS=''
fi
if test x${ac_cv_env_CPPFLAGS_set} != xset; then
  CPPFLAGS=''
fi

#--------- Custom checks ----------

# check for openssl/md5.h
AC_CHECK_HEADERS([openssl/md5.h], [has_md5=yes], [AC_MSG_WARN([

|---------------------------------------------------------------------------------------------|
|                                                                                             |
|   ***WARNING: Header openssl/md5.h not found. No hash will be given to configurations!***   |
|                                                                                             |
|---------------------------------------------------------------------------------------------|

])])

# if openssl/md5.h is found, define ENABLE_MD5_HASH macro and link the library
if test x${has_md5} = xyes; then
  AC_DEFINE(ENABLE_MD5_HASH, , [Enable computation of the hash of the configuration])
fi

AM_CONDITIONAL([LINK_MD5_H], [test x${has_md5} = xyes])

#----------------------------------

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT

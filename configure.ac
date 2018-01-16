#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([libdecaf], [0.9.4], [krion@FreeBSD.org])
AC_CONFIG_SRCDIR([aux/like_eddsa.cxx])
# AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([1.9])
AC_PREFIX_DEFAULT([/usr/local])

COMPILE_ARGS="${ac_configure_args}"
AC_DEFINE_UNQUOTED(COMPILE_ARGS, "$COMPILE_ARGS")

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_LN_S
AC_PROG_MAKE_SET

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lc':
AC_CHECK_LIB([c], [main])
# FIXME: Replace `main' with a function in `-ldecaf':
AC_CHECK_LIB([decaf], [main])

dnl initial checks; tweakings to CFLAGS and LDFLAGS
dnl see final checks for tweakings to LIBS
if test "x$ac_cv_prog_gcc" = xyes ; then
  CFLAGS="-O2 ${CFLAGS}"
fi

dnl set debug level
AC_MSG_CHECKING([whether to enable debugging compiler options])
AC_ARG_ENABLE(debug,
        [  --enable-debug                   Enable debugging compiler options (default: no)],
        [
                if test x$enableval = x"yes" ; then
                        AC_MSG_RESULT(yes)
                        tmp_CFLAGS=`echo $CFLAGS | sed 's/O2/O0/g'`
                        CFLAGS="$tmp_CFLAGS"
                        CFLAGS="$CFLAGS -g -W -Wall"
                else
                        AC_MSG_RESULT(no)
                fi
        ],
        [
                AC_MSG_RESULT(no)
        ]
)

dnl set relax level
AC_MSG_CHECKING([whether to relax compiler optimizations])
AC_ARG_ENABLE(relax,
[  --enable-relax                   Relax compiler optimization (default: no)],
        [
                if test x$enableval = x"yes" ; then
                        AC_MSG_RESULT(yes)
                        tmp_CFLAGS=`echo $CFLAGS | sed 's/O2/O0/g'`
                        CFLAGS="$tmp_CFLAGS"
                else
                        AC_MSG_RESULT(no)
                fi
        ],
        [
                AC_MSG_RESULT(no)
        ]
)

dnl: shared object handling
AC_MSG_CHECKING([whether to disable linking against shared objects])
AC_ARG_ENABLE(so,
        [  --disable-so                     Disable linking against shared objects (default: no)],
        [
        if test x$enableval = x"yes" ; then
          AC_MSG_RESULT(no)
          AC_CHECK_FUNC(dlopen, [ USING_DLOPEN="yes" ], [])
          AC_CHECK_LIB(dl, dlopen,
                [ USING_DLOPEN="yes"
                  LIBS="${LIBS} -ldl" ],
                  [])
          if test x"$USING_DLOPEN" != x"yes"; then
                AC_MSG_ERROR(Unable to find dlopen(). Try with --disable-so)
          fi
        else
          AC_MSG_RESULT(yes)
          if test "x$ac_cv_prog_gcc" = xyes ; then
                LDFLAGS="-static ${LDFLAGS}"
          fi
        fi
        ],
        [
        AC_MSG_RESULT(no)
        AC_CHECK_FUNC(dlopen, [ USING_DLOPEN="yes" ], [])
        AC_CHECK_LIB(dl, dlopen,
                [ USING_DLOPEN="yes"
                  LIBS="${LIBS} -ldl" ],
                  [])
        if test x"$USING_DLOPEN" != x"yes"; then
                AC_MSG_ERROR(Unable to find dlopen(). Try with --disable-so)
        fi
        ]
)
dnl finish: shared object handling

dnl os specific flags
case "$host_os" in
  Sun*)
    AC_DEFINE(SOLARIS, 1)
    LIBS="-lresolv -lsocket -lnsl ${LIBS}"
    ;;
  *BSD)
    AC_DEFINE(BSD, 1)
    ;;
esac

dnl cpu specific flags
case "$host_cpu" in
  sun*)
  AC_DEFINE(CPU_sparc, 1)
  ;;
esac

AC_CHECK_PROG(MAKE, gmake, gmake)
if test x"$MAKE" = x""; then
  AC_CHECK_PROG(MAKE, make, make)
fi
AC_SUBST(MAKE)
AC_PROG_MAKE_SET

dnl
dnl some systems have __progname ; if this is the case and we play around argv
dnl we need to enable a minor hack to make things work nicely.
dnl
AC_MSG_CHECKING(for __progname)
AC_TRY_LINK([ extern char *__progname; ],
[ __progname = "test"; ],
[AC_MSG_RESULT(yes); AC_DEFINE(PROGNAME, 1)], [AC_MSG_RESULT(no)])

dnl
dnl Some checks to understand whether we need to instruct the linker for
dnl exporting collector symbols to dynamically loaded classifiers.
dnl
dnl OS'es with ELF executables using the GNU linker (Linux and recent *BSD,
dnl in rare cases Solaris) typically need '-Wl,-export-dynamic'; some SYSv4
dnl systems instead need '-Wl,-Bexport'; AIX 4.x wants -Wl,-bexpall,-brtl'.
dnl

AC_MSG_CHECKING(for extra flags needed to export symbols)
if test "x$ac_cv_prog_gcc" = xyes ; then
        save_ldflags="${LDFLAGS}"
        LDFLAGS="-Wl,--export-dynamic ${save_ldflags}"
        AC_TRY_LINK(,,
                [ AC_MSG_RESULT(--export-dynamic) ],
                [ LDFLAGS="-Wl,-Bexport ${save_ldflags}"
                  AC_TRY_LINK(,,
                        [ AC_MSG_RESULT(-Bexport) ],
                        [ AC_MSG_RESULT(none)
                          LDFLAGS="${save_ldflags}" ]
                  )
                ]
        )
else
        AC_MSG_RESULT(none)
fi

dnl
dnl some C compilers (ex. MipsPRO) don't know how to parse the 'inline' keyword
dnl
AC_MSG_CHECKING(for static inline)
AC_TRY_COMPILE([#include <stdio.h>]

static inline func()
{
}
,
[
  func();
],
[AC_MSG_RESULT(yes)], [AC_MSG_RESULT(no); AC_DEFINE(NOINLINE, 1)])

dnl
dnl Check for architecture endianess: big | little
dnl
dnl XXX: switch to manually define this feature
ac_cv_endianess="unknown"
if test x"$ac_cv_endianess" = x"unknown"; then
  AC_MSG_CHECKING(endianess)
  AC_TRY_RUN([main () {

  union
  {
    long l;
    char c[sizeof (long)];
  } u;
  u.l = 1;
  exit (u.c[sizeof (long) - 1] == 1);
  }], ac_cv_endianess="little", ac_cv_endianess="big", ac_cv_endianess="little")
  AC_MSG_RESULT($ac_cv_endianess)
fi

if test x"$ac_cv_endianess" = x"big"; then
  AC_DEFINE(IM_BIG_ENDIAN, 1)
fi
if test x"$ac_cv_endianess" = x"little"; then
  AC_DEFINE(IM_LITTLE_ENDIAN, 1)
fi

dnl
dnl Check for unaligned memory access; based entirely over
dnl AC_LBL_UNALIGNED_ACCESS
dnl
dnl XXX: switch to manually define this feature;
ac_cv_unaligned="unknown"
case "$host_cpu" in
  alpha*|arm*|hp*|mips*|sh*|sparc*|ia64|nv1)
    ac_cv_unaligned="fail"
    AC_MSG_CHECKING(unaligned accesses)
    AC_MSG_RESULT($ac_cv_unaligned)
    ;;
esac

if test x"$ac_cv_unaligned" = x"unknown"; then
  AC_MSG_CHECKING(unaligned accesses)
  cat > conftest.c << EOF
    #include <sys/types.h>
    #include <sys/wait.h>
    #include <stdio.h>
    unsigned char a[[5]] = { 1, 2, 3, 4, 5 };
    main () {
        unsigned int i;
        pid_t pid;
        int status;
        /* avoid "core dumped" message */
        pid = fork();
        if (pid <  0) exit(2);
        if (pid > 0) {
                /* parent */
                pid = waitpid(pid, &status, 0);
                if (pid < 0) exit(3);
                exit(!WIFEXITED(status));
        }
        /* child */
        i = *(unsigned int *)&a[[1]];
        printf("%d\n", i);
        exit(0);
    }
EOF
  ${CC-cc} -o conftest $CFLAGS $CPPFLAGS $LDFLAGS \
        conftest.c $LIBS >/dev/null 2>&1
  if test ! -x conftest ; then
    ac_cv_unaligned="fail"
  else
    ./conftest >conftest.out
    if test ! -s conftest.out ; then ac_cv_unaligned="fail"
    else ac_cv_unaligned="ok"
    fi
  fi
  rm -f conftest* core core.conftest
  AC_MSG_RESULT($ac_cv_unaligned)
fi

if test x"$ac_cv_unaligned" = x"fail"; then
  AC_DEFINE(NEED_ALIGN, 1)
fi

# Checks for header files.
AC_CHECK_HEADERS([fcntl.h inttypes.h limits.h stdint.h stdlib.h string.h sys/param.h sys/time.h unistd.h])

dnl Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_TYPE(u_int64_t, [AC_DEFINE(HAVE_U_INT64_T, 1)])
AC_CHECK_TYPE(u_int32_t, [AC_DEFINE(HAVE_U_INT32_T, 1)])
AC_CHECK_TYPE(u_int16_t, [AC_DEFINE(HAVE_U_INT16_T, 1)])
AC_CHECK_TYPE(u_int8_t, [AC_DEFINE(HAVE_U_INT8_T, 1)])
AC_CHECK_TYPE(uint64_t, [AC_DEFINE(HAVE_UINT64_T, 1)])
AC_CHECK_TYPE(uint32_t, [AC_DEFINE(HAVE_UINT32_T, 1)])
AC_CHECK_TYPE(uint16_t, [AC_DEFINE(HAVE_UINT16_T, 1)])
AC_CHECK_TYPE(uint8_t, [AC_DEFINE(HAVE_UINT8_T, 1)])

AC_MSG_CHECKING([whether to enable 64bit counters])
AC_ARG_ENABLE(64bit,
        [  --enable-64bit                   Enable 64bit counters (default: yes)],
        [
          if test x$enableval = x"yes" ; then
            AC_MSG_RESULT(yes)
            AC_DEFINE(HAVE_64BIT_COUNTERS, 1)
          else
            AC_MSG_RESULT(no)
          fi
        ],
        [
          AC_MSG_RESULT(yes)
          AC_DEFINE(HAVE_64BIT_COUNTERS, 1)
        ]
)

# Checks for library functions.
AC_TYPE_SIGNAL
AC_FUNC_MALLOC
AC_CHECK_FUNCS([gettimeofday memmove memset pow sqrt])

dnl final checks
dnl trivial solution to portability issue
CFLAGS="${CFLAGS} ${INCLUDES}"
INCLUDES=""

echo "
PLATFORM ..... : `uname -m`
OS ........... : `uname -rs` (`uname -n`)
COMPILER ..... : ${CC}
CFLAGS ....... : ${CFLAGS}
LIBS ......... : ${LIBS}
LDFLAGS ...... : ${LDFLAGS}
PLUGINS ...... : ${PLUGINS}

Now type 'make' to compile the source code.

Need for documentation and examples?
Start by reading the README.md file

Star, watch or contribute to the project on sourceforge:
https://sourceforge.net/projects/ed448goldilocks/

"

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

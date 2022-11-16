bin_PROGRAMS = cpn \
               multicanonic_cpn \
               cpn_conf_check

cpn_conf = ./include/cpn_cmplx_op.h \
           ./lib/cpn_cmplx_op.c \
           ./include/cpn_conf.h \
           ./lib/cpn_conf_def.c \
           ./lib/cpn_meas.c \
           ./lib/cpn_update.c

cpn_param = ./include/cpn_param.h \
            ./lib/cpn_param.c

geometry = ./include/geometry.h \
           ./lib/geometry.c

cpn_upd_multicanonic = ./include/cpn_multicanonic.h \
                       ./lib/cpn_update_multicanonic.c 

rng = ./include/rng.h \
      ./lib/rng.c

endianness = ./include/endianness.h \
             ./lib/endianness.c

cpn_SOURCES = ./src/cpn.c \
              ./include/macro.h \
              ${rng} \
              ${cpn_conf} \
              ${cpn_param} \
              ${geometry} \
              ${endianness}

multicanonic_cpn_SOURCES = ./src/multicanonic_cpn.c \
                           ./include/macro.h \
                           ${rng} \
                           ${cpn_conf} \
                           ${cpn_param} \
                           ${geometry} \
                           ${endianness} \
                           ${cpn_upd_multicanonic}

cpn_conf_check_SOURCES = ./src/cpn_conf_check.c \
                         ./include/macro.h \
                         ${rng} \
                         ${cpn_conf} \
                         ${cpn_param} \
                         ${geometry} \
                         ${endianness}

if LINK_MD5_H 
LDADD = -lm -lcrypto
else
LDADD = -lm
endif

AM_CFLAGS = -std=c99 -O3 -Wall -Wextra -pedantic -Wconversion -Werror

CLEANFILES = *.dat* *.dat_back* *.dat_replica*

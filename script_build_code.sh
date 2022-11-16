#!/bin/bash

# build code
aclocal
autoconf
autoheader
automake --add-missing --foreign

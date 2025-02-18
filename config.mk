OPENCILK_BIN?=
CC=$(OPENCILK_BIN)clang
CXX=$(OPENCILK_BIN)clang++
CILKFLAG=-fopencilk
OPT?=-O2
DBG?=-g3
MACHINE?=$(shell uname -m)
ARCH_x86_64=-march=x86-64-v3	# uname -m on Linux and Darwin
ARCH_amd64=-march=x86-64-v3	# uname -m on FreeBSD
ARCH_aarch64=
ARCH=$(ARCH_$(MACHINE))
TOOL=

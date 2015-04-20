#!/bin/bash
#################################################################################
#
# MAKEFILE WITH AUTOMATIC DEPENDENCIES for program prova.f90
#
#################################################################################
#
# Date of last update: 08/05/2013
#
#
# Author: Marco Giometto
#
# Environment: Ubuntu 12.04
#
#
# Description:
# 
# Some precompiler flags are there and used in the OUTPUT_SPC.f90, this is done
# so that anyone who uses the code can implement his output in specific standalone 
# modules without affecting the memory allocation at compile time of the code.
#
#
#################################################################################

#executable name
EXE = program.exe

#directory for the .o files
OBJDIR = ./src/.obj
MODDIR = ./src/.mod

FOBJ=./src/.obj/MAIN.o 

#directories for sources
SRCS = ./src/*.f90

#settings for DORA with INTEL compiler (ifort)
FC = ftn
FFLAGS = -O3

all: clean createdir program.exe

# preprocessor flags settings
COMPSTR = '$$(FC) -c -o $$@ $$(FFLAGS) $$<'

createdir:
	mkdir -p ./src/.obj ./src/.mod

program.exe: $(FOBJ)
	$(FC) -o $@ $(FFLAGS) $(LDFLAGS) $(FOBJ) $(LIBS)

./src/.obj/MAIN.o : ./src/MAIN.f90 
	$(FC) -c -o $@ $(FFLAGS) $<

#cleaning out everything
clean:
	rm -Rf ./src/.obj
	rm -Rf ./src/.mod
	rm -f program.exe

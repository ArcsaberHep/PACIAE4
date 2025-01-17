################################################################################
# Makefile is a part of the PACIAE event generator.
# Copyright (C) 2024 PACIAE Group.
# PACIAE is licensed under the GNU GPL v2 or later, see LICENCE for details.
# Open source: https://github.com/ArcsaberHep/PACIAE4
# Author: An-Ke Lei, January 2024 - November 2024.
#
# This is is the Makefile used to complie and build PACIAE program linked with
#  PYTHIA 8 on POSIX systems.
#
# Usage:
#  Type "make -j" to genarate PACIAE running program.
#  Type "make PACIAE -j" to genarate both PACIAE and averaging programs.
#  Type "make clean" to clean temporary files (.o, .log, .out etc.).
#
# Note that you need "make, gfortran and g++" installed on your system.
# For help using the make command please consult the local system documentation,
#  i.e. "man make" or "make --help".
#
#                                               By An-Ke at CCNU on 16/01/2024
#                                  Last updated by An-Ke at UiO  on 17/01/2025
################################################################################

################################################################################
# VARIABLES: Definition of the relevant variables and SHELL commands.
################################################################################

# Sets the shell.
SHELL=/usr/bin/env bash

# Copys "../pythia8311/Makefile.inc" to "./" if "Makefile.inc" does not exist.
ifeq ("$(wildcard ./Makefile.inc)","")
  ifneq ("$(wildcard ../pythia8311/Makefile.inc)","")
    CP_MAKE_INC:=$(shell cp ../pythia8311/Makefile.inc ./)
  endif
endif

# Includes the configuration from PYTHIA 8.
-include Makefile.inc

# Checks PYTHIA 8 distribution (use local version first, then installed one).
ifneq ("$(wildcard ../pythia8311/lib/libpythia8.*)","")
  PREFIX_LIB:=../pythia8311/lib
  PREFIX_INCLUDE:=../pythia8311/include
endif

# PACIAE distribution.
DIR_TOP:=$(shell pwd)
DIR_INC:=./include
DIR_LIB:=
# DIR_LIB:=./lib
DIR_SRC:=./src
DIR_BIN:=./bin
DIR_OBJ:=./obj
DIR_ETC:=./etc
DIR_DOC:=
# DIR_DOC:=./doc
DIR_LOG:=./log
DIR_SIM:=./sim

# Creates folders and move files.
LOCAL_MKDIR1:=$(shell mkdir -p $(DIR_INC) $(DIR_LIB) $(DIR_SRC))
LOCAL_MKDIR2:=$(shell mkdir -p $(DIR_BIN) $(DIR_OBJ) $(DIR_ETC))
LOCAL_MKDIR3:=$(shell mkdir -p $(DIR_DOC) $(DIR_LOG) $(DIR_SIM))

# Firstly, deltes /bin/xPaciae40.x, /bin/xRms_analysis.x
RM_BIN_X:=$(shell rm -f ./bin/xPaciae40.x ./bin/xRms_analysis.x)

# Moves "./usu.dat" to "./etc" if it exists.
ifneq ("$(wildcard ./usu.dat)","")
  MV_USU1:=$(shell mv -f ./usu.dat $(DIR_ETC)/)
endif
# Moves "./src/usu.dat" to "./etc" if it exists.
ifneq ("$(wildcard $(DIR_SRC)/usu.dat)","")
  MV_USU2:=$(shell mv -f $(DIR_SRC)/usu.dat $(DIR_ETC)/)
endif

# Moves "./pythia6_extra.cfg" to "./etc" if it exists.
ifneq ("$(wildcard ./pythia6_extra.cfg)","")
  MV_PY6CFG1:=$(shell mv -f ./pythia6_extra.cfg $(DIR_ETC)/)
endif
# Moves "./src/pythia6_extra.cfg" to "./etc" if it exists.
ifneq ("$(wildcard $(DIR_SRC)/pythia6_extra.cfg)","")
  MV_PY6CFG2:=$(shell mv -f $(DIR_SRC)/pythia6_extra.cfg $(DIR_ETC)/)
endif

# Move "./pythia8_extra.cfg" to "./etc" as if it exists.
ifneq ("$(wildcard ./pythia8_extra.cfg)","")
  MV_PY8CFG1:=$(shell mv -f ./pythia8_extra.cfg $(DIR_ETC)/)
endif
# Move "./src/pythia8_extra.cfg" to "./etc" as if it exists.
ifneq ("$(wildcard $(DIR_SRC)/pythia8_extra.cfg)","")
  MV_PY8CFG2:=$(shell mv -f $(DIR_SRC)/pythia8_extra.cfg $(DIR_ETC)/)
endif

# Moves source code to "./src" if they exist.
SRC_FILE:= ./*.f* ./*.F* ./*.c* ./*.C*
ifneq ("$(filter-out ./%.cfg, $(wildcard $(SRC_FILE)))","")
  MV_SRC:=$(shell mv -t $(DIR_SRC)/ $(SRC_FILE))
endif

# Moves header files to "./include" if they exist.
INC_FILE:= ./*.h* ./*.H*
ifneq ("$(wildcard $(INC_FILE))","")
  MV_INC:=$(shell mv -t $(DIR_INC)/ $(INC_FILE))
endif

# Moves log files to "./log" if they exist.
LOG_FILE:= ./*.log
ifneq ("$(wildcard $(LOG_FILE))","")
  MV_LOG:=$(shell mv -t $(DIR_LOG)/ $(LOG_FILE))
endif

# Cleans the log files in "./sim" folder.
SIM_LOG_CLEAN:=$(shell rm -f $(DIR_SIM)/*.log)

# Compilation flags. Two different compilers would be used: gfortran and g++.
# Part of them were set in "Makefile.inc".
# C++ compilation flags.
CXX_COMMON:=-I$(PREFIX_INCLUDE) -I$(DIR_INC) $(CXX_COMMON) $(GZIP_LIB) -w
PYTHIA8=$(PREFIX_LIB)/libpythia8$(LIB_SUFFIX)
# Fortran compilation flags.
GFO:=gfortran
#debug
#GFO_COMMON:=-W -Wall -Wshadow -g -fdefault-integer-8 -fbounds-check
GFO_COMMON:=-O2 -W -Wall -Wshadow -g -fdefault-integer-8 -w

# Link flags. Note linker is gfortran.
GFO_LINK:=$(GFO_COMMON)
GFO_LINK+=-I$(PREFIX_INCLUDE) -I$(DIR_INC)
GFO_LINK+=-L$(DIR_LIB) -Wl,-rpath,$(DIR_LIB)
GFO_LINK+=-L$(PREFIX_LIB) -Wl,-rpath,$(PREFIX_LIB) -lpythia8
GFO_LINK+=-ldl -lstdc++
#debug
# GFO_LINK+= -fbacktrace
# GFO_LINK+= -ffpe-trap=invalid,zero,overflow,underflow,denormal,inexact

################################################################################
# RULES: Definition of the rules used to build the PACIAE40.
################################################################################

# PACIAE40 program, i.e. the name of executable file.
TARGET:=xPaciae40.x
TARG_AVE:=xRms_analysis.x
OBJECTS:=$(DIR_OBJ)/Main_40.o $(DIR_OBJ)/Parini_40.o
OBJECTS+=$(DIR_OBJ)/Parcas_40.o $(DIR_OBJ)/Sfm_40.o
OBJECTS+=$(DIR_OBJ)/Coales_40.o $(DIR_OBJ)/Hadcas_40.o
OBJECTS+=$(DIR_OBJ)/Analy_40.o $(DIR_OBJ)/P_40.o
OBJECTS+=$(DIR_OBJ)/Pythia8_fort_interface.o $(DIR_OBJ)/Pythia8CppInterface.o
OBJECTS+=$(DIR_OBJ)/PaciaeUserHooks.o
# OBJECTS=$( wildcard $(DIR_OBJ)/*.o )
OBJ_AVE:=Rms_analysis.o

# Rules without physical targets (secondary expansion for specific rules).
#.SECONDEXPANSION:

.PHONY: all copy_x copy_usu copy_py6cfg copy_py8cfg

# All targets (no default behavior).
all: $(DIR_BIN)/$(TARGET) copy_x copy_usu copy_py6cfg copy_py8cfg
	$(info )
	$(info PACIAE40 program generation succeeded!)
	$(info )

# Copys executable file to "./sim".
copy_x: $(DIR_BIN)/$(TARGET)
	cp -f $(DIR_BIN)/$(TARGET) $(DIR_SIM)

# Copys "./etc/usu.dat" to "./sim" if it does not exist.
copy_usu: $(DIR_ETC)/usu.dat
	echo "no" | cp -i $(DIR_ETC)/usu.dat $(DIR_SIM)/usu.dat

# Copys "./etc/pythia6_extra.cfg" to "./sim" if it does not exist.
copy_py6cfg: $(DIR_ETC)/pythia6_extra.cfg
	echo "no" | cp -i $(DIR_ETC)/pythia6_extra.cfg $(DIR_SIM)/pythia6_extra.cfg

# Copys "./etc/pythia8_extra.cfg" to "./sim" if it does not exist.
copy_py8cfg: $(DIR_ETC)/pythia8_extra.cfg
	echo "no" | cp -i $(DIR_ETC)/pythia8_extra.cfg $(DIR_SIM)/pythia8_extra.cfg

.PHONY: PACIAE

PACIAE: $(DIR_BIN)/$(TARGET) $(DIR_BIN)/$(TARG_AVE)
	$(info )
	$(info PACIAE40 program generation succeeded!)
	$(info )

# PYTHIA library.
$(PYTHIA8):
	$(error Error: PYTHIA8 must be built, please run "./configure" and                 "make" in the top PYTHIA8 directory)

# Generate "xPaciae40.x" by linking "PACIAE (Fortran 77/90 + C++)" with
#  "Pythia8_fort_interface (Fortran 90/2003)", "Pythia8CppInterface (C++)",
#  "PYTHIA 6 (FORTRAN 77)" and "PYTHIA 8 (C++)" altogether.
$(DIR_BIN)/$(TARGET): $(OBJECTS)
	$(GFO) -o $@ $^ $(GFO_LINK)

# Generate "Main_40.o".
$(DIR_OBJ)/Main_40.o: $(DIR_SRC)/Main_40.f90
	$(GFO) -c -o $@ $^ $(GFO_COMMON)

# Generate "Parini_40.o".
$(DIR_OBJ)/Parini_40.o: $(DIR_SRC)/Parini_40.f90
	$(GFO) -c -o $@ $^ $(GFO_COMMON)

# Generate "Parcas_40.o".
$(DIR_OBJ)/Parcas_40.o: $(DIR_SRC)/Parcas_40.f90
	$(GFO) -c -o $@ $^ $(GFO_COMMON)

# Generate "Sfm_40.o".
$(DIR_OBJ)/Sfm_40.o: $(DIR_SRC)/Sfm_40.f90
	$(GFO) -c -o $@ $^ $(GFO_COMMON)

# Generate "Coales_40.o".
$(DIR_OBJ)/Coales_40.o: $(DIR_SRC)/Coales_40.f90
	$(GFO) -c -o $@ $^ $(GFO_COMMON)

# Generate "Hadcas_40.o".
$(DIR_OBJ)/Hadcas_40.o: $(DIR_SRC)/Hadcas_40.f90
	$(GFO) -c -o $@ $^ $(GFO_COMMON)

# Generate "Analy_40.o".
$(DIR_OBJ)/Analy_40.o: $(DIR_SRC)/Analy_40.f90
	$(GFO) -c -o $@ $^ $(GFO_COMMON)

# Generate "P_40.o".
$(DIR_OBJ)/P_40.o: $(DIR_SRC)/P_40.f
	$(GFO) -c -o $@ $^ $(GFO_COMMON)

# Generate "Pythia8_fort_interface.o".
$(DIR_OBJ)/Pythia8_fort_interface.o: $(DIR_SRC)/Pythia8_fort_interface.f90
	$(GFO) -c -o $@ $^ $(GFO_COMMON)

# Generate "Pythia8CppInterface.o".
$(DIR_OBJ)/Pythia8CppInterface.o: $(DIR_SRC)/Pythia8CppInterface.cpp
	$(CXX) -c -o $@ $^ $(CXX_COMMON)

# Generate "PaciaeUserHooks.o".
$(DIR_OBJ)/PaciaeUserHooks.o: $(DIR_SRC)/PaciaeUserHooks.cpp
	$(CXX) -c -o $@ $^ $(CXX_COMMON)

# Plugin libraries.
# lib%.so: %.cc $(PYTHIA)
# 	$(CXX) $< -o $@ -w $(CXX_COMMON) $(CXX_SHARED) -Wl,--no-as-needed

# Generate "xRms_analysis.x".
$(DIR_BIN)/$(TARG_AVE): $(DIR_OBJ)/$(OBJ_AVE)
	$(GFO) -o $@ $^ $(GFO_COMMON)

# Generate "Rms_analysis.o".
$(DIR_OBJ)/$(OBJ_AVE): $(DIR_SRC)/Rms_analysis.f90
	$(GFO) -c -o $@ $^ $(GFO_COMMON)

# Clean.
.PHONY: clean
clean:
	rm -f $(DIR_OBJ)/*.o; rm -f $(DIR_OBJ)/*.mod;
	rm -f $(DIR_BIN)/*.x; rm -f $(DIR_LOG)/*.log;
	rm -f $(DIR_SIM)/*.o; rm -f $(DIR_SIM)/*.mod;
	rm -f $(DIR_SIM)/*.x; rm -f $(DIR_SIM)/*.log;
	rm -f $(DIR_SIM)/*.out;
	rm -f `ls $(DIR_SIM)/*.dat|egrep -v $(DIR_SIM)/usu.dat`;
	rm -f `ls $(DIR_SIM)/*.cfg|egrep -v $(DIR_SIM)/pythia6_extra.cfg`;
	rm -f `ls $(DIR_SIM)/*.cfg|egrep -v $(DIR_SIM)/pythia8_extra.cfg`;
	rm -f *.o *.mod *.x *.log *.dat *.cfg;

# Clean "Rms_analysis.o" and "xRms_analysis.x".
.PHONY: clean_rms_analysis
clean_rms_analysis:
	rm -f $(DIR_OBJ)/Rms_analysis.o; rm -f $(DIR_BIN)/xRms_analysis.x

# Clean "Rms_analysis.o" and "xRms_analysis.x".
.PHONY: clean_sim
clean_sim:
	rm -f $(DIR_SIM)/*.o; rm -f $(DIR_SIM)/*.mod;
	rm -f $(DIR_SIM)/*.x; rm -f $(DIR_SIM)/*.log;
	rm -f $(DIR_SIM)/*.out;
	rm -f `ls $(DIR_SIM)/*.dat|egrep -v $(DIR_SIM)/usu.dat`;
	rm -f `ls $(DIR_SIM)/*.cfg|egrep -v $(DIR_SIM)/pythia6_extra.cfg`;
	rm -f `ls $(DIR_SIM)/*.cfg|egrep -v $(DIR_SIM)/pythia8_extra.cfg`;

# Clean "$(DIR_LOG)/*.log".
.PHONY: clean_log
clean_log:
	rm -f $(DIR_LOG)/*.log;

# Clean "bin/*.x".
.PHONY: clean_bin
clean_bin:
	rm -f $(DIR_BIN)/*.x;

# Clean "obj/*.o obj/*.mod".
.PHONY: clean_obj
clean_obj:
	rm -f $(DIR_OBJ)/*.o $(DIR_OBJ)/*.mod;

# Clean "sim directory".
.PHONY: clean_SIM_dir
clean_SIM_dir:
	rm -rI $(DIR_SIM);

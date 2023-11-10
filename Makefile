# **************************************************************
# *   Based on C++ Makefile Template by Arash Partow (2003)    *
# *                                                            *
# * URL: http://www.partow.net/programming/makefile/index.html *
# * http://www.opensource.org/licenses/MIT                     *
# *                                                            *
# **************************************************************

CXX      := -g++
CXXFLAGS := -pedantic-errors -Wall -Wextra -Werror

#################
# library flags #
#################

# GCC_LIB = location of /lib64 directory in GCC installation
BASE_LDFLAGS  := -L${GCC_LIB} -lstdc++ -lm -std=c++17
# OpenMPI_LIB = location of /lib directory in OpenMPI installation
OPENMPI_LDFLAGS := -L${OpenMPI_LIB} -lmpi -pthread -Wl,-rpath -Wl,${OpenMPI_LIB} -Wl,--enable-new-dtags
# fmt_LIB = location of /lib directory in fmt installation
FMT_LDFLAGS := -L${fmt_LIB} -lfmt
# LAMMPS_LIB location of /lib directory in LAMMPS installation
LAMMPS_LDFLAGS := ${LAMMPS_LIB}/liblammps_twistable_BD_OMP.so

LDFLAGS := ${BASE_LDFLAGS} ${OPENMPI_LDFLAGS} ${FMT_LDFLAGS} ${LAMMPS_LDFLAGS}

################
# header flags #
################

# GCC_INC = location of /include directory in GCC installation
BASE_INCLUDE  := -Iinclude/ -I${GCC_INC}
# OpenMPI_INC = location of /include directory in OpenMPI installation
OPENMPI_INCLUDE := -I${OpenMPI_INC}
# fmt_INC = location of /include directory in fmt installation
FMT_INCLUDE := -I${fmt_INC}
# location of header files in LAMMPS source code
LAMMPS_INCLUDE := -I/home/ben/Software/Suites/LAMMPS/lammps/src

INCLUDE := ${BASE_INCLUDE} ${OPENMPI_INCLUDE} ${FMT_INCLUDE} ${LAMMPS_INCLUDE}

#####################
# build information #
#####################

BUILD    := ./build
OBJ_DIR  := $(BUILD)/objects
APP_DIR  := $(BUILD)/apps

TARGET   := btree_chromo

SRC      := $(wildcard src/*.cpp)

OBJECTS  := $(SRC:%.cpp=$(OBJ_DIR)/%.o)
DEPENDENCIES \
         := $(OBJECTS:.o=.d)

all: build $(APP_DIR)/$(TARGET)

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -MMD -o $@

$(APP_DIR)/$(TARGET): $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -o $(APP_DIR)/$(TARGET) $^ $(LDFLAGS)

-include $(DEPENDENCIES)

.PHONY: all build clean debug release info

build:
	@mkdir -p $(APP_DIR)
	@mkdir -p $(OBJ_DIR)

debug: CXXFLAGS += -DDEBUG -g -ggdb3
debug: all

release: CXXFLAGS += -O2
release: all

clean:
	-@rm -rvf $(OBJ_DIR)/*
	-@rm -rvf $(APP_DIR)/*

info:
	@echo "[*] Application dir: ${APP_DIR}     "
	@echo "[*] Object dir:      ${OBJ_DIR}     "
	@echo "[*] Sources:         ${SRC}         "
	@echo "[*] Objects:         ${OBJECTS}     "
	@echo "[*] Dependencies:    ${DEPENDENCIES}"
	@echo "[*] Include:         ${INCLUDE}     "

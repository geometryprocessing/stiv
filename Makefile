
# Things for users to define
#
#include make.inc


IAGM_LIB := libiagm.a

CXX := mpic++
CXXFLAGS := -std=c++11 -Wall -Wextra -pipe -fPIC -O3 -fopenmp # -shared -no-inline-max-size -no-inline-max-total-size

EIGEN_INC := /opt/local/include/eigen3
INCS := -Iinc -Isrc -I$(EIGEN_INC)

#LIB_DIR := -L/opt/intel/compilers_and_libraries_2018.2.164/mac/compiler/lib
#LIBS := -lCGAL -lgmp -llapack -liomp5
#LIBS := -lCGAL -lgmp -llapack #-lgomp

VPATH = ./src
SRCS := $(wildcard ./src/*.cpp)

BUILD_DIR := build/
OBJS := $(addprefix $(BUILD_DIR), $(addsuffix .o, $(basename $(notdir $(SRCS)))))


all: $(IAGM_LIB)

$(IAGM_LIB): $(OBJS)
	mkdir -p lib
	/bin/rm -f lib/$(IAGM_LIB)
	ar ruv lib/$(IAGM_LIB) $(OBJS)
	#$(CXX) -dynamiclib -Wl,-install_name,./lib/$(IAGM_LIB) -o lib/$(IAGM_LIB) $(OBJS) $(LIBS) $(LIB_DIR)

$(BUILD_DIR)%.o: %.cpp
	$(shell mkdir -p $(BUILD_DIR))
	$(CXX) $(CXXFLAGS) $(INCS) -c -o $@ $<

depend:
	$(CXX) $(INCS) -MM $(SRCS) | sed 's/.*\.o/build\/&/' > .depend

ifeq (.depend, $(wildcard .depend))
  include .depend
endif

clean: 
	-rm -f lib/$(IAGM_LIB) $(BUILD_DIR)*.o .depend

rebuild: clean depend all


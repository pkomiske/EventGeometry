CXX = g++
SRCS = $(shell ls *.cpp)
CXXFLAGS = -O3 -Wall -std=c++14 -g -ffast-math

ifeq ($(shell uname), Darwin)
	CXXFLAGS += -Xpreprocessor -fopenmp
	LIBRARIES += -lomp
endif
ifeq ($(shell uname), Linux)
	CXXFLAGS += -fopenmp
	LDFLAGS += -fopenmp
endif

FASTJETINC = $(shell fastjet-config --cxxflags)
FASTJETLIB = $(shell fastjet-config --libs)
INCLUDES  += $(FASTJETINC) -I.. -I../Wasserstein -I../Wasserstein/examples/include
LIBRARIES += $(FASTJETLIB) -lz

CNPY = ../Wasserstein/examples/src/cnpy.o

# http://make.mad-scientist.net/papers/advanced-auto-dependency-generation/#combine
DEPDIR = .deps
COMPILE.cpp = $(CXX) $(CXXFLAGS) $(INCLUDES) -c $(<) -o $(@) -MT $(@) -MMD -MP -MF $(DEPDIR)/$*.d
%.o : %.cpp
%.o : %.cpp $(DEPDIR)/%.d | $(DEPDIR)
	$(COMPILE.cpp)

.PHONY: all clean
all: emd_example pairwise_emds_example

emd_example: emd_example.o $(CNPY)
	$(CXX) -o $@ $^ $(LIBRARIES) $(LDFLAGS)

pairwise_emds_example: pairwise_emds_example.o $(CNPY)
	$(CXX) -o $@ $^ $(LIBRARIES) $(LDFLAGS)

$(CNPY):
	make -C ../Wasserstein/examples src/cnpy.o

clean:
	rm -rfv *.o *_example $(DEPDIR) ../Wasserstein/examples/src/*.o

$(DEPDIR): ; @mkdir -p $@
DEPFILES := $(SRCS:%.cpp=$(DEPDIR)/%.d)
$(DEPFILES):

include $(wildcard $(DEPFILES))

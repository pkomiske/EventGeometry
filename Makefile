# If you are using this Makefile standalone and fastjet-config is not
# in your path, edit this line to specify the full path
FASTJETCONFIG=fastjet-config
PREFIX=$(shell $(FASTJETCONFIG) --prefix)
CXX=g++
CXXFLAGS+=-O3 -Wall -g -std=c++14 -fPIC -Xpreprocessor -fopenmp
install_script = $(SHELL) ../utils/install-sh
check_script = ../utils/check.sh

# global contrib-wide Makefile include may override some of the above
# variables (leading "-" means don't give an error if you can't find
# the file)
-include ../.Makefile.inc

#------------------------------------------------------------------------
# things that are specific to this contrib
NAME=EventGeometry
SRCS=
EXAMPLES=example
INSTALLED_HEADERS=EventGeometry.hh Wasserstein/wasserstein/CorrelationDimension.hh Wasserstein/wasserstein/EMD.hh
INSTALLED_INTERNAL_HEADERS=Wasserstein/wasserstein/internal/EMDUtils.hh Wasserstein/wasserstein/internal/Event.hh Wasserstein/wasserstein/internal/HistogramUtils.hh Wasserstein/wasserstein/internal/PairwiseDistance.hh Wasserstein/wasserstein/internal/NetworkSimplex.hh
FASTJET_SWIG_INTERFACE=/usr/local/opt/fastjet-3.3.4/pyinterface/fastjet.i
SWIG_NUMPY=true
SWIGINTERFACE=EventGeometry.i
SWIGFLAGS=-w509,511 -keyword -fastproxy -IWasserstein/wasserstein
#------------------------------------------------------------------------

CXXFLAGS+= $(shell $(FASTJETCONFIG) --cxxflags)
LDFLAGS += -lm $(shell $(FASTJETCONFIG) --libs)

OBJS  := $(SRCS:.cc=.o)
EXAMPLES_SRCS  = $(EXAMPLES:=.cc)

install_HEADER  = $(install_script) -c -m 644
install_LIB     = $(install_script) -c -m 644
install_DIR     = $(install_script) -d
install_DATA    = $(install_script) -c -m 644
install_PROGRAM = $(install_script) -c -s
install_SCRIPT  = $(install_script) -c

ifdef SWIGINTERFACE
SWIG_SRC=Py$(NAME).cc
SRCS+=$(SWIG_SRC)
endif
ifdef SWIG_NUMPY
SWIGINTERFACE+=Wasserstein/swig/numpy.i
SWIGFLAGS+=-DSWIG_NUMPY
endif
SWIGINTERFACE+=$(FASTJET_SWIG_INTERFACE)
SWIGFLAGS+=-DFASTJET_SWIG_INTERFACE=$(FASTJET_SWIG_INTERFACE)

.PHONY: clean distclean examples check install all swig

# http://make.mad-scientist.net/papers/advanced-auto-dependency-generation/#combine
DEPDIR = .deps
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.d
COMPILE.cc = $(CXX) $(DEPFLAGS) $(CXXFLAGS) -c

%.o : %.cc
%.o : %.cc $(DEPDIR)/%.d | $(DEPDIR)
	$(COMPILE.cc) $(OUTPUT_OPTION) $<

# compilation of the code (default target)
all: lib$(NAME).a

lib$(NAME).a: $(OBJS) 
	ar cru lib$(NAME).a $(OBJS)
	ranlib lib$(NAME).a

# building the examples
examples: $(EXAMPLES)

# the following construct makes it possible to automatically build
# each of the examples listed in $EXAMPLES
$(EXAMPLES): % : %.o
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

# check that everything went fine
check: examples
	@for prog in $(EXAMPLES); do\
	  $(check_script) $${prog} ../data/single-event.dat || exit 1; \
	done
	@echo "All tests successful"

# cleaning the directory
clean:
	rm -fv *~ *.o *.a *.so Py$(NAME).cc $(NAME).py
	rm -rfv __pycache__ build .deps

distclean: clean
	rm -f lib$(NAME).a $(EXAMPLES)

# install things in PREFIX/...
install:
	$(install_DIR) $(PREFIX)/include/fastjet/contrib
	$(install_DIR) $(PREFIX)/include/fastjet/contrib/internal
	for header in $(INSTALLED_HEADERS); do\
	  $(install_HEADER) $$header $(PREFIX)/include/fastjet/contrib/;\
	done
	for header in $(INSTALLED_INTERNAL_HEADERS); do\
	  $(install_HEADER) $$header $(PREFIX)/include/fastjet/contrib/internal;\
	done

depend:
	makedepend -Y --   -- $(SRCS) $(EXAMPLES_SRCS)

swig: $(SWIG_SRC)

$(SWIG_SRC): $(SWIGINTERFACE) $(INSTALLED_HEADERS) $(INSTALLED_INTERNAL_HEADERS)
	swig -c++ -python $(SWIGFLAGS) -I$(PREFIX)/include -o $@ $<

$(DEPDIR): ; @mkdir -p $@

DEPFILES := $(SRCS:%.cc=$(DEPDIR)/%.d)

$(DEPFILES):

include $(wildcard $(DEPFILES))

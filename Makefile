#-------------------------------------------------------------
# Makefile for cmc++ library
#-------------------------------------------------------------
#include ../root_dir.mk # must be included first
include ./make_options.mk
#-------------------------------------------------------------
# Source files
SRCS = scheduler/mpi_comm.cpp 
SRCS+= scheduler/cmdargs.cpp 
SRCS+= scheduler/inputparams.cpp 
SRCS+= scheduler/taskparams.cpp 
SRCS+= scheduler/worker.cpp 
SRCS+= scheduler/master_scheduler.cpp
SRCS+= scheduler/scheduler.cpp
#SRCS+= xml/pugixml.cpp 
SRCS+= expression/complex_expression.cpp
SRCS+= utils/utils.cpp 
SRCS+= mcdata/mcdata.cpp
SRCS+= mcdata/mc_observable.cpp
SRCS+= lattice/lattice.cpp
SRCS+= lattice/latticelibrary.cpp
SRCS+= model/quantum_op.cpp
SRCS+= model/strmatrix.cpp
SRCS+= model/model_term.cpp
SRCS+= model/model.cpp
SRCS+= model/modellibrary.cpp
SRCS+= diag/hamiltonian.cpp
SRCS+= diag/kspace.cpp
SRCS+= diag/kspace_kpath.cpp
SRCS+= diag/bandstruct.cpp
SRCS+= diag/wavefunction.cpp
SRCS+= diag/observables.cpp
SRCS+= diag/diag.cpp
SRCS+= main.cpp
VMC_SRCS = $(addprefix src/,$(SRCS))
#-------------------------------------------------------------
# Headers
HDRS=    scheduler/mpi_comm.h \
         scheduler/optionparser.h scheduler/cmdargs.h \
         scheduler/inputparams.h scheduler/worker.h scheduler/task.h \
         scheduler/scheduler.h \
         expression/complex_expression.h \
         utils/utils.h \
         utils/curve_fit.h \
	 mcdata/mcdata.h  \
	 mcdata/mc_observable.h  \
         lattice/constants.h lattice/lattice.h \
         model/model_params.h  model/quantum_op.h \
	 model/strmatrix.h \
	 model/model_term.h \
	 model/quantum_op.h \
	 model/model.h \
	 diag/matrix.h \
	 diag/kspace.h \
	 diag/hamiltonian.h \
	 diag/bandstruct.h \
	 diag/wavefunction.h \
	 diag/observables.h \
	 diag/diag.h \
	 vmc/bandstruct.h \
#-------------------------------------------------------------
VMC_HDRS = $(addprefix src/,$(HDRS))
MUPARSER_LIB = $(PROJECT_ROOT)/src/expression/muparserx/libmuparserx.a
#-------------------------------------------------------------
# Target
ifeq ($(MPI), HAVE_BOOST_MPI)
  TAGT=diag_mpi.x
else
  TAGT=diag.x
endif

# Put all auto generated stuff to this build dir.
ifeq ($(BUILD_DIR), $(CURDIR))
  $(error In-source build is not allowed, choose another build directory)
endif

# All .o files go to BULD_DIR
OBJS=$(patsubst %.cpp,$(BUILD_DIR)/%.o,$(VMC_SRCS))
# GCC/Clang will create these .d files containing dependencies.
DEPS=$(patsubst %.o,%.d,$(OBJS)) 

.PHONY: all
all: $(TAGT) #$(INCL_HDRS)

$(TAGT): $(OBJS) $(MUPARSER_LIB)
	$(CXX) -o $(TAGT) $(OBJS) $(LDFLAGS) $(LIBS) $(MUPARSER_LIB)

%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) -o $@ $<

# Include all .d files
-include $(DEPS)

$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	@echo "$(CXX) -c $(CXXFLAGS) -o $(@F) $(<F)"
	@$(CXX) -MMD -c $(CXXFLAGS) -o $@ $<

$(VMC_INCLDIR)/%.h: %.h 
	@mkdir -p $(@D)
	@echo "Copying $< to 'include'" 
	@cp -f $< $@

$(MUPARSER_LIB):
	@cd ./src/expression/muparserx/ && $(MAKE)

# installation
#prefix = ../install#/usr/local
#libdir = $(prefix)/lib
#includedir = $(prefix)/include/cmc++

.PHONY: install
install:	
	@echo "Already installed in $(VMC_LIBDIR) and $(VMC_INCLDIR)" 

.PHONY: clean
clean:	
	@echo "Removing temporary files in the build directory"
	@rm -f $(OBJS) $(DEPS) 
	@echo "Removing $(TAGT)"
	@rm -f $(TAGT) 

.PHONY: bclean
bclean:	
	@echo "Removing temporary files in the build directory"
	@rm -f $(OBJS) $(DEPS) 

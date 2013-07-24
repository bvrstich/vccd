###############################################################################
#
#  makefile template for the sources
#
###############################################################################

# -----------------------------------------------------------------------------
#   Sources for all modules
# -----------------------------------------------------------------------------
BINNAME = mps_btas
CPPSRC	= mps_btas.cpp\
           Tensor.cpp\
           MPS.cpp\
           Tools.cpp


OBJ	= $(CPPSRC:.cpp=.o)

# -----------------------------------------------------------------------------
#   These are the standard libraries, include paths and compiler settings
# -----------------------------------------------------------------------------

BRIGHT_ROOT= .

BLITZINC=`pkg-config --cflags blitz`
BLITZLIB=`pkg-config --libs blitz`

BOOSTINC=-I/usr/include/boost
BOOSTLIB=-L$(BOOSTDIR)/lib -lboost_serialization

BTASINC=-I/usr/include/btas
BTASLIB= /usr/lib

INCLUDE = ./include $(BOOSTINC) $(BLITZINC) $(BTASINC)

LIBS= -llapack -lblas $(BOOSTLIB) $(BLITZLIB) $(BTASLIB)/libbtas.a

CC	= gcc
CXX	= g++

# -----------------------------------------------------------------------------
#   Compiler & Linker flags
# -----------------------------------------------------------------------------
CFLAGS	= -I$(INCLUDE) -g -O3 -flto -fopenmp
LDFLAGS	= -g -O3 -flto -fopenmp

# =============================================================================
#   Targets & Rules
# =============================================================================
all:
	@echo
	@echo '  +++ Building $(BINNAME)...'
	@echo	
	$(MAKE) $(BRIGHT_ROOT)/$(BINNAME)
	@if test $?; then \
	   echo; echo '*************** FAILED! ***************' ; echo; \
	 else \
	   echo; echo '  +++ $(BINNAME) has been built successfully!'; \
	   echo; \
	 fi

# -----------------------------------------------------------------------------
#   The default way to compile all source modules
# -----------------------------------------------------------------------------
%.o:	%.for makefile
	@echo; echo "Compiling $(@:.o=.for) ..."
	$(FF) -c $(FFLAGS) $(SFLAGS) $(@:.o=.for) -o $@

%.o:	%.c makefile
	@echo; echo "Compiling $(@:.o=.c) ..."
	$(CC) -c $(CFLAGS) $(SFLAGS) $(@:.o=.c) -o $@

%.o:	%.cpp makefile
	@echo; echo "Compiling $(@:.o=.cpp) ..."
	$(CXX) -c $(CFLAGS) $(SFLAGS) $(DEFS) $(@:.o=.cpp) -o $@

# -----------------------------------------------------------------------------
#   Link everything together
# -----------------------------------------------------------------------------
$(BRIGHT_ROOT)/$(BINNAME):	makefile $(OBJ) 
	@echo; echo "Linker: creating $(BRIGHT_ROOT)/$(BINNAME) ..."
	$(CXX) $(LDFLAGS) $(SFLAGS) -o $(BRIGHT_ROOT)/$(BINNAME) $(OBJ) $(LIBS)

# -----------------------------------------------------------------------------
#   Create everything newly from scratch
# -----------------------------------------------------------------------------
new:	clean all

# -----------------------------------------------------------------------------
#   Clean up all object files
# -----------------------------------------------------------------------------
clean:
	@echo -n '  +++ Cleaning all object files ... '
	@echo -n $(OBJ)
	@rm -f $(OBJ)
	@echo 'Done.'

# -----------------------------------------------------------------------------
#   Make new documentation using doxygen
# -----------------------------------------------------------------------------
doc:
	@doxygen doc-config

# ====================== End of file 'makefile.in' ========================== #

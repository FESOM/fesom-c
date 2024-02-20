# Define compiler
FC = gfortran

# Define source files
SOURCES = fv_wrapper.F90 fv_d3.F90

# Define object files
OBJECTS = $(SOURCES:.F90=.o)

# Target to create, including a pre-compilation clean
all: clean $(OBJECTS)

# Rule to create object files
%.o: %.F90
	$(FC) -c $<

# Clean command to remove object and module files
clean:
	rm -f $(OBJECTS) *.mod

# Suffixes
.SUFFIXES: .F90 .o

# Dummy targets
.PHONY: all clean


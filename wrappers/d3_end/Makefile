# Define compiler
FC = ifort

# Define source files
SOURCES = fv_wrapper.F90 fv_d3.F90 fv_run_d3.F90

# Define object files
OBJECTS = $(SOURCES:.F90=.o)

# Name of the executable
EXECUTABLE = fv_run_d3

# Target to create, including a pre-compilation clean
all: clean $(EXECUTABLE)

# Rule to compile and link the main program into an executable
$(EXECUTABLE): $(OBJECTS)
	$(FC) -o $@ $(OBJECTS)

# Rule to create object files
%.o: %.F90
	$(FC) -c $<

# Clean command to remove object and module files
clean:
	rm -f $(OBJECTS) *.mod $(EXECUTABLE)

# Suffixes
.SUFFIXES: .F90 .o

# Dummy targets
.PHONY: all clean

# The compiler
FC = gfortran 
EXEC = barotrp
FCFLAGS = -g -fbounds-check
FCFLAGS += -I/usr/local/netcdf/include/

 
# libraries needed for linking, unused in the examples
LDFLAGS = -L/usr/local/netcdf/lib/ -lnetcdf

# List of executables to be built within the package	
OBJECTS =  module_constants.o module_io.o module_solve.o baro_trp.o
	
all: $(OBJECTS)
	$(FC) -o $(EXEC) $(OBJECTS) $(FCFLAGS) $(LDFLAGS)
		
%.o: %.f90
	$(FC) -c $< $(FCFLAGS) $(LDFLAGS)

clean:
	rm $(OBJECTS) $(EXEC) *.mod
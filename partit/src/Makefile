######################################################
# Makefile                
######################################################

######################################################
# Include architecture-specific definitions

include Makefile.in

###### Objects for Mesh Partitioning ################################################
# modules
MOD_INI =  fort_part.o  oce_modules.o  gen_modules_config.o gen_modules_partitioning.o gen_modules_rotate_grid.o  

OBJ_INI =  fvom_init.o \
           oce_local.o \
           gen_comm.o


# Name of executable

EXE_INI = fvom_ini.x

# TARGETS

default: run_ini

run_ini: CPP_DEFS+=-DFVOM_INIT
run_ini: cleanomod $(MOD_INI) $(OBJ_INI)
	@echo "======= Building FESOM paritioning program =========="
	$(LD)  $(OPT) -o $(EXE_INI)  $(MOD_INI) $(OBJ_INI) $(LIB_METIS)
	cp -pf $(EXE_INI) ../bin/.

.SUFFIXES: .c .F90 .o 

.c.o :
	$(CC) $(COPT) $(METIS_DEFS) $(CPP_DEFS) $(PARMS_DEFS) $(METIS_INC) $(PARMS_INC) $(MPI_INCLUDE) -c $*.c	

#.F90.o :
#	@echo $@
#	$(CPP) $(CPP_DEFS) $(CPP_SOL) $(PETSC_INC) $(PARMS_INC) $(PETSCCONF_INC) $(FOS_INC) $(MPI_INC) $(NCINC) $(LIBS_SLV) < $*.F90 > $*.f90
#	$(FC) $(OPT) $(CPP_SOL) $(CPP_DEFS) $(OASIS3_INC) $(FOS_INC)\
#	$(PARMS_INC) $(MPI_INCLUDE) $(NCINC)  -c $*.f90
.F90.o :
	@echo $@
	$(FC) $(MPI_INC) $(NCINC) $(LIBS_SLV) $(OPT) $(CPP_SOL) $(CPP_DEFS) $(OASIS3_INC) $(OMP_DEFS)\
	$(PARMS_INC) $(MPI_INCLUDE) $(NCINC)  -c $*.F90

clean :
	rm -f *.o *.mod *~ *.f90 $(EXE_INI)

cleanomod:
	rm -f *.o *.mod
	CPP_DEFS=; export CPP_DEFS

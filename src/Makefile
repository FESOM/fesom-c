#Makefile

include Makefile.in

.SUFFIXES: .F90 .o

LIBS = $(GOTMLIB)  $(NETCDFLIB) $(PNETCDFLIB)
INCS = $(GOTMINCS) $(NETCDFINC) $(PNETCDFINC)

OBJECTS= fv_var.o                  \
      gen_modules_partitioning_c.o \
      gen_modules_rotate_grid_c.o  \
      gen_halo_exchange_c.o        \
      fv_utilit.o               \
      oce_mesh_c.o                 \
      fv_vert_coord.o           \
      fv_sbc.o                  \
      fv_obc.o                  \
      fv_rivers.o               \
      fv_ic.o                   \
      fv_read_run.o             \
      fv_ini.o                  \
      fv_mesh_array.o           \
      fv_gradients.o            \
      fv_3D_vel.o               \
      fv_advection.o            \
      fv_average_dynamic.o      \
      fv_mixing.o               \
      fv_output.o               \
      fv_pressure_new.o         \
      fv_viscosity.o            \
      fv_tracer.o               \
      fv_restart.o               \
      fv_tide_potential.o       \
      fv_sbcmpi.o          \
      fv_ncoutput.o             \
      fv_ncoutputmpi.o          \
      fv_ice_setup_step.o \
      fv_sediment.o      \
      fv_sediment_adv.o      \
      fv_output_nc_TMP.o \
      fv_main.o  




INCS_INI = $(METISINC)
LIBS_INI = $(METISLIB)

OBJ_C_INI = metis_wrapper.o  
SRC_INI =  fv_var.F90 \
	   fv_utilit.F90 \
	   fv_read_init.F90 \
	   fv_mesh_array.F90 \
	   fv_distribute_mesh.F90 \
           fv_comm.F90 \
           fv_init.F90


#--------------------------------------------------

#$(EXE): $(OBJECTS)  $(GOTMDIR)/lib/libturbulence.a
$(EXE): $(OBJECTS) 
	 $(LINK) $(FCFLAGS) -o $(EXE) $(OBJECTS) $(LIBS)

$(EXE_INI): $(OBJ_C_INI) $(SRC_INI)
	 $(LINK) $(FCFLAGS) -o $(EXE_INI) $(SRC_INI) $(OBJ_C_INI) $(LIBS_INI)

#GOTM_turbulence: 
#	$(MAKE) -C $(GOTMDIR) all

all:
#	$(MAKE) GOTM_turbulence
	$(MAKE) $(EXE)

clean :
	rm -f *~ *.o *.mod $(EXE) $(EXE_INI)

allclean :
	$(MAKE) clean
#	$(MAKE) -C $(GOTMDIR) realclean


run : $(EXE)         
	./$(EXE)

#--------------------------------------------------

.F90.o: $(MODULES)
	$(FC)  $(FCFLAGS) $(INCS) -c $*.F90

.c.o:
	$(CC) $(CCFLAGS) $(INCS_INI) -c $*.c

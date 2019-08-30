# beps.make
# make file to compile beps program

# CFLAGS=-std=c99 -O3
 CFLAGS=-std=c99 -g

# first name a variable objects for all object files
objects = aerodynamic_conductance.o beps401b_glb_hr.o calc_temp_leaf.o \
  inter_prg.o evaporation_canopy.o evaporation_soil.o init_soil.o \
  lai_calc.o meteo_pack.o netRadiation.o photosyn_gs.o \
  rainfall.o readparam_l.o s_coszs.o plant_respir_l.o\
  sensible_heat.o snowpack.o soil_thermal_regime.o \
  soil_water_stress.o surface_temp.o soilresp.o\
  transpiration.o updatesoilmoisture.o vcmax_jmax_mod.o \
  readcontrolf_b.o readlonlat.o read_yrdata_c.o read_input2.o readcoef_l.o read_all_lai.o

# name a variable sources for all source files
sources = aerodynamic_conductance.c beps401b_glb_hr.c calc_temp_leaf.c \
  inter_prg.c evaporation_canopy.c evaporation_soil.c init_soil.c \
  lai_calc.c meteo_pack.c netRadiation.c photosyn_gs.c \
  rainfall.c readparam_l.c s_coszs.c plant_respir_l.c\
  sensible_heat.c snowpack.c soil_thermal_regime.c \
  soil_water_stress.c surface_temp.c soilresp.c\
  transpiration.c updatesoilmoisture.c vcmax_jmax_mod.c \
  readcontrolf_b.c readlonlat.c read_yrdata_c.c read_input2.c readcoef_l.c read_all_lai.c

# now give target as makebeps with objects as variable dependencies + command line

beps: $(objects) 
#	icc -g -lm -o bepsb $(objects)
	gcc $(CFLAGS) -o bepsb $(objects) -lm

# list the dependencies for object files - those header files which help build objects

# how to build all object files from all dependent source files

$(objects): $(sources)
#	icc -g -c $(sources)
	gcc $(CFLAGS) -c $(sources)

clean:
	rm $(objects) bepsb


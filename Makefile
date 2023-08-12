.SUFFIXES: .f .f90 .o

F90 = gfortran # compiile with gfortran

CMPLFLG = -c -O3 -fbacktrace

OBJS = Multipoles.o LongRangePotential.o General_Coordinates_Format.o Coordinate_Transformation.o Index_Searcher.o Approx_1_Sph2.o Approx_2_Sph2.o Approx_3_Sph2.o Approx_4_Sph2.o Approx_5_Sph2.o Approx_6_Sph2.o Approx_7_Sph2.o Approx_8_Sph2.o Induction_4_Sph2.o  Induction_5_Sph2.o Induction_6_Sph2.o Induction_7_Sph2.o Induction_8_Sph2.o HyperPolarizability_6_Sph2.o HyperPolarizability_7_Sph2.o Dispersion_6_Sph2.o Dispersion_7_Sph2.o Dispersion_8_Sph2.o T_lk.o T_l0.o T_ll.o T_saved.o
OBJS_Test =  test.o LongRangePotential.o General_Coordinates_Format.o Coordinate_Transformation.o Index_Searcher.o Approx_1_Sph2.o Approx_2_Sph2.o Approx_3_Sph2.o Approx_4_Sph2.o Approx_5_Sph2.o Approx_6_Sph2.o Approx_7_Sph2.o Approx_8_Sph2.o Induction_4_Sph2.o  Induction_5_Sph2.o Induction_6_Sph2.o Induction_7_Sph2.o Induction_8_Sph2.o HyperPolarizability_6_Sph2.o HyperPolarizability_7_Sph2.o Dispersion_6_Sph2.o Dispersion_7_Sph2.o Dispersion_8_Sph2.o T_lk.o T_l0.o T_ll.o T_saved.o

all : $(OBJS)
	$(F90) $(OBJS) -o ejec.x
#	rm *.o 
#	rm *.mod 

build-run-test : $(OBJS_Test)
	$(F90) $(OBJS_Test) -o ejec_test.x
	./ejec_test.x

clean:
	del *.exe
	del *.x
	del *.o
	del *.mod

$(OBJS) :
.f90.o:
	$(F90) $(CMPLFLG) $<
.f.o:
	$(F90) $(CMPLFLG) $<



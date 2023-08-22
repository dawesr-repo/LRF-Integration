.SUFFIXES: .f .f90 .o

F90 = gfortran # compiile with gfortran

CMPLFLG = -c -O3 -fbacktrace

OBJS = 		 T_saved.o General_Coordinates_Format.o Coordinate_Transformation.o Index_Searcher.o Approx_1_Sph2.o Approx_2_Sph2.o Approx_3_Sph2.o Approx_4_Sph2.o Approx_5_Sph2.o Approx_6_Sph2.o Approx_7_Sph2.o Approx_8_Sph2.o Induction_4_Sph2.o  Induction_5_Sph2.o Induction_6_Sph2.o Induction_7_Sph2.o Induction_8_Sph2.o HyperPolarizability_6_Sph2.o HyperPolarizability_7_Sph2.o Dispersion_6_Sph2.o Dispersion_7_Sph2.o Dispersion_8_Sph2.o T_lk.o T_l0.o T_ll.o LongRangePotential.o evaluateLRF.o min_example.o 
OBJS_Test =  T_saved.o General_Coordinates_Format.o Coordinate_Transformation.o Index_Searcher.o Approx_1_Sph2.o Approx_2_Sph2.o Approx_3_Sph2.o Approx_4_Sph2.o Approx_5_Sph2.o Approx_6_Sph2.o Approx_7_Sph2.o Approx_8_Sph2.o Induction_4_Sph2.o  Induction_5_Sph2.o Induction_6_Sph2.o Induction_7_Sph2.o Induction_8_Sph2.o HyperPolarizability_6_Sph2.o HyperPolarizability_7_Sph2.o Dispersion_6_Sph2.o Dispersion_7_Sph2.o Dispersion_8_Sph2.o T_lk.o T_l0.o T_ll.o LongRangePotential.o test.o 


OBJS_PROD =  production.o  min_example.o 
OBJS_PROD_Test =  production.o    test.o
OBJS_PROD_Build =  build_production.o 


all : $(OBJS)
	$(F90) $(OBJS) -o ejec.x
	./ejec.x


build-run-test : $(OBJS_Test)
	$(F90) $(OBJS_Test) -o ejec_test.x
	./ejec_test.x

build-prod : $(OBJS_PROD_Build)
	$(F90) $(OBJS_PROD_Build) -o ejec_build_prod.x
	./ejec_build_prod.x

run-prod : $(OBJS_PROD)
	$(F90) $(OBJS_PROD) -o ejec_prod.x
	./ejec_prod.x

run-prod-test : $(OBJS_PROD_Test)
	$(F90) $(OBJS_PROD_Test) -o ejec_test_prod.x
	./ejec_test_prod.x

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



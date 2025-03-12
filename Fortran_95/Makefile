.SUFFIXES: .f .f90 .o

F90 = gfortran # compile with gfortran

CMPLFLG = -c -O3 -fbacktrace


OBJS = 		 Fitting_Constant_v2.o Geometry_Constant_v2.o helper_functions.o min_example.o
OBJS_Test =  Fitting_Constant_v2.o Geometry_Constant_v2.o helper_functions.o Testing_v2.o test_v2.o


OBJS_PROD =       LRF.o  min_example.o 
OBJS_PROD_Test =  LRF.o  Testing_v2.o test_v2.o
OBJS_PROD_Build = build_production.o 



all : $(OBJS)
	$(F90) $(OBJS) -o ejec.x
	./ejec.x


test : $(OBJS_Test)
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



.SUFFIXES: .f .f90 .o

F90 = gfortran # compiile with gfortran

CMPLFLG = -c -O3 -fbacktrace

OBJS = 		 Geometry_Constant_v2.o FittingConstant.o helperFunc.o Multipole_Sph3.o Induction_Sph3.o Dispersion_Sph3.o min_example.o 
OBJS_Test =  Geometry_Constant_v2.o FittingConstant.o helperFunc.o Multipole_Sph3.o Induction_Sph3.o Dispersion_Sph3.o test_v2.o


OBJS_PROD =       LRF.o  min_example.o 
OBJS_PROD_Test =  LRF.o  test.o
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



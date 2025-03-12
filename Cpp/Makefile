.SUFFIXES: .cpp .h .o

CPP = g++ # compile with g++

CMPLFLG = -c


OBJS = 		 main.o
OBJS_Test =  Fitting_Constant_v2.o Geometry_Constant_v2.o helper_functions.o Testing_v2.o test_v2.o


OBJS_PROD =       LRF.o  min_example.o 
OBJS_PROD_Test =  LRF.o  Testing_v2.o test_v2.o
OBJS_PROD_Build = build_production.o 



all : $(OBJS)
	$(CPP) $(OBJS) -o ejec.x
	./ejec.x


test : $(OBJS_Test)
	$(CPP) $(OBJS_Test) -o ejec_test.x
	./ejec_test.x

build-prod : $(OBJS_PROD_Build)
	$(CPP) $(OBJS_PROD_Build) -o ejec_build_prod.x
	./ejec_build_prod.x

run-prod : $(OBJS_PROD)
	$(CPP) $(OBJS_PROD) -o ejec_prod.x
	./ejec_prod.x

run-prod-test : $(OBJS_PROD_Test)
	$(CPP) $(OBJS_PROD_Test) -o ejec_test_prod.x
	./ejec_test_prod.x

clean:
	del *.exe
	del *.x
	del *.o
	del *.mod

$(OBJS) :
.cpp.o:
	$(CPP) $(CMPLFLG) $<
.h.o:
	$(CPP) $(CMPLFLG) $<



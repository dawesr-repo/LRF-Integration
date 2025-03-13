.SUFFIXES: .cpp .h .o

CPP = g++ # compile with g++

CMPLFLG = -c


OBJS = 		 calculate_tensor.o interactions.o evaluate_lrf.o potential_energy_surface.o  main.o

all : $(OBJS)
	$(CPP) $(OBJS) -o ejec.x
	./ejec.x


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



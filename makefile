all:twogrid_POD BuildFEM.o

OBJS1=FEMheat.o src/BuildFEM.o
OBJS3=twogrid_POD.o src/BuildFEM.o

include ${PHG_MAKEFILE_INC}

twogrid_POD.o: /opt/local-MVAPICH2/phg-0.9.2/lib/libphg.so twogrid_POD.c fun.h
FEMheat.o: /opt/local-MVAPICH2/phg-0.9.2/lib/libphg.so FEMheat.c fun.h

FEMheat:/opt/local-MVAPICH2/phg-0.9.2/lib/libphg.so ${OBJS1}
	${LINKER} ${USER_LDFLAGS} ${LDFLAGS} -o $@ ${OBJS1}${USER_LIB} ${LIBS}

twogrid_POD:/opt/local-MVAPICH2/phg-0.9.2/lib/libphg.so ${OBJS3}
	${LINKER} ${USER_LDFLAGS} ${LDFLAGS} -o $@ ${OBJS3}${USER_LIB} ${LIBS}
clean:
	rm -f *.o *.text FEMheat twogrid_POD 


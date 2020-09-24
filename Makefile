EXECS=fox
MPICC?=mpicc

all: fox

fox: fox.c
        ${MPICC} -o fox fox.c


clean:
        rm -f fox

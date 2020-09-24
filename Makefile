EXECS=lab5
MPICC?=mpicc

all: lab5

ratetest: fox.c
        ${MPICC} -o fox fox.c


clean:
        rm -f lab5

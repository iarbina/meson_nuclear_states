# Makefile to run program Klein-Gordon Equation solver

FC=gfortran

SRC=src
OBJ=obj
BIN=bin

OBJECTS=inp_data.o spls3.o potentials.o Numerov.o nuc-modular.o

exeprog.exe: ${OBJECTS}
	${FC} -o ${BIN}/exeprog.exe $(addprefix ${OBJ}/, ${OBJECTS})

inp_data.o: ${SRC}/inp_data.f90
	${FC} -J${OBJ} -c ${SRC}/inp_data.f90 -o ${OBJ}/inp_data.o

spls3.o: ${SRC}/spls3.f
	${FC} -c ${SRC}/spls3.f -o ${OBJ}/spls3.o

potentials.o: ${SRC}/potentials.f90
	${FC} -J${OBJ} -c ${SRC}/potentials.f90 -o ${OBJ}/potentials.o

Numerov.o: ${SRC}/Numerov.f90
	${FC} -I${OBJ} -c ${SRC}/Numerov.f90 -o ${OBJ}/Numerov.o

nuc-modular.o: ${SRC}/nuc-modular.f90
	${FC} -I${OBJ} -c ${SRC}/nuc-modular.f90 -o ${OBJ}/nuc-modular.o

clean:
	rm -f ${BIN}/exeprog.exe ${OBJ}/*.mod ${OBJ}/*.o

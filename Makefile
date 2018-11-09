# Makefile to run program Klein-Gordon Equation solver

FC=gfortran

DEB=-cpp -DDEBUG=1 -DDEBUG_FILE=1000
SRC=src
OBJ=obj
BIN=bin

OBJECTS=ansi_colors.o inp_data.o spls.o spls3.o potentials.o Numerov.o nuc-modular.o

exeprog.exe: ${OBJECTS}
	${FC} ${DEB} -o ${BIN}/exeprog.exe $(addprefix ${OBJ}/, ${OBJECTS})

ansi_colors.o: ${SRC}/ansi_colors.f90
	${FC} ${DEB} -J${OBJ} -c ${SRC}/ansi_colors.f90 -o ${OBJ}/ansi_colors.o

inp_data.o: ${SRC}/inp_data.f90
	${FC} ${DEB} -J${OBJ} -c ${SRC}/inp_data.f90 -o ${OBJ}/inp_data.o

spls.o: ${SRC}/spls.f90
	${FC} ${DEB} -c ${SRC}/spls.f90 -o ${OBJ}/spls.o

spls3.o: ${SRC}/spls3.f
	${FC} ${DEB} -c ${SRC}/spls3.f -o ${OBJ}/spls3.o

potentials.o: ${SRC}/potentials.f90
	${FC} ${DEB} -J${OBJ} -c ${SRC}/potentials.f90 -o ${OBJ}/potentials.o

Numerov.o: ${SRC}/Numerov.f90
	${FC} ${DEB} -I${OBJ} -c ${SRC}/Numerov.f90 -o ${OBJ}/Numerov.o

nuc-modular.o: ${SRC}/nuc-modular.f90
	${FC} ${DEB} -I${OBJ} -c ${SRC}/nuc-modular.f90 -o ${OBJ}/nuc-modular.o

clean:
	rm -f ${BIN}/exeprog.exe ${OBJ}/*.mod ${OBJ}/*.o fort.1000

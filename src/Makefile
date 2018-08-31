# Makefile to run program Klein-Gordon Equation solver

FC=gfortran

OBJ=inp_data.o spls3.o potentials.o Numerov.o nuc-modular.o

exeprog: ${OBJ}
	${FC} -o exeprog inp_data.o spls3.o potentials.o Numerov.o nuc-modular.o

inp_data.o: inp_data.f90
	${FC} -c inp_data.f90

spls3.o: spls3.f
	${FC} -c spls3.f

potentials.o: potentials.f90
	${FC} -c potentials.f90

Numerov.o: Numerov.f90
	${FC} -c Numerov.f90

nuc-modular.o: nuc-modular.f90
	${FC} -c nuc-modular.f90

clean:
	rm -f exeprog *.mod *.o

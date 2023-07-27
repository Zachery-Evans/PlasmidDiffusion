#SRC=.
#TARGET=./target
#CFLAGS=-g -Wall -Wextra 
ARGS= ''

all:
	nvcc doublePlasmid.cu -o mc.run
	nvcc geometry.cu -o g.run
	nvcc uncertainty.cu -o u.run
	./mc.run $(ARGS)

default:
	main.cu

main.c:
	nvcc $(CFLAGS) $(SRC)/*.cu

test:
	main.cu
	$(TARGET)/sim.run $(ARGS)

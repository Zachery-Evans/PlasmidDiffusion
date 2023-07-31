#SRC=.
#TARGET=./target
#CFLAGS=-o $(TARGET)/sim.run
ARGS= ''

all:
	nvcc -lm -lcurand testing.cu -o mc.run
	./mc.run

default:
	main.cu

main.c:
	nvcc $(CFLAGS) $(SRC)/*.cu

test:
	main.cu
	$(TARGET)/sim.run $(ARGS)

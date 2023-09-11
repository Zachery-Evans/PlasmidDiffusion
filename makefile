#SRC=.
#TARGET=./target
#CFLAGS=-o $(TARGET)/sim.run
ARGS= ''

all:
	nvcc -lcurand -lcudadevrt -rdc=true -Xlinker=-rpath,$$CUDA_PATH/lib64 testing.cu -o mc.run
	./mc.run

default:
	main.cu

main.c:
	nvcc $(CFLAGS) $(SRC)/*.cu

test:
	main.cu
	$(TARGET)/sim.run $(ARGS)

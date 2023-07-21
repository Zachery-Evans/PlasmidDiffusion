#SRC=.
#TARGET=./target
#CFLAGS=-g -Wall -Wextra 
ARGS= ''

all:
	cpp -lm -O3 doublePlasmid.cpp -o mc.run
	./mc.run $(ARGS)
	
default:
	main.cpp

main.c:
	ccpp $(CFLAGS) $(SRC)/*.cpp

test:
	main.cpp
	$(TARGET)/sim.run $(ARGS)

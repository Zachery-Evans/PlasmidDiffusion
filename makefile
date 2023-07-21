#SRC=.
#TARGET=./target
#CFLAGS=-g -Wall -Wextra 
ARGS= ''

all:
	g++ -pthread -lm -O3 doublePlasmid.cpp -o mc.run
	./mc.run $(ARGS)

default:
	main.cpp

main.c:
	cpp $(CFLAGS) $(SRC)/*.cpp

test:
	main.cpp
	$(TARGET)/sim.run $(ARGS)

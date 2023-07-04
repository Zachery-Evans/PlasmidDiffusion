#SRC=.
#TARGET=./target
#CFLAGS=-g -Wall -Wextra 
ARGS= ''

all:
	cc -lm -O3 doublePlasmid.c -o mc.run
	cc -lm -O3 geometry.c -o geo.run
	cc -lm -O3 uncertainty.c -o unc.run
	./mc.run $(ARGS)

default:
	main.c

main.c:
	cc $(CFLAGS) $(SRC)/*.c

test:
	main.c
	$(TARGET)/sim.run $(ARGS)

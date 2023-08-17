#SRC=.
#TARGET=./target
#CFLAGS=-g -Wall -Wextra 
ARGS= ''

all:
	cc -lm -O3 radiusgyration.c -o mc.run
	./mc.run $(ARGS)

default:
	main.c

main.c:
	cc $(CFLAGS) $(SRC)/*.c

test:
	main.c
	$(TARGET)/sim.run $(ARGS)

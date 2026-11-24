CC      = mpicc
CFLAGS  = -Wall -O2
NP     ?= 4
ARGS   ?= 10 7 5 2.0

producto_matrices_1d: producto_matrices_1d.c
	$(CC) $(CFLAGS) $< -o $@

run: producto_matrices_1d
	mpirun -np $(NP) --oversubscribe ./producto_matrices_1d $(ARGS)

clean:
	rm -f producto_matrices_1d

.PHONY: run clean

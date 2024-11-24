all: db5242

db5242: db5242.c
	gcc -mavx2 -o db5242 db5242.c

clean:
	rm -f db5242
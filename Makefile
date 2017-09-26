# standard compile options for the c++ executable
CC = gcc
EXE = mosher
OBJS =  main.o myrand.o particles.o physics.o
CFLAGS = -O3 -Wall -g
LDFLAGS = -g
LIBS = -lm

# default super-target
all: $(EXE)

# the standard executable
$(EXE): $(OBJS)
	$(CC) $(LDFLAGS) $^ -o $@ $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c $<

main.o:      main.c myrand.h particles.h
myrand.o:    myrand.c    myrand.h
particles.o: particles.c particles.h
physics.o:   physics.c   physics.h

code.pdf: code.md
	pandoc $< -o $@

code.md: particles.h physics.h myrand.h particles.c physics.c main.c
	cp code_notes.md code.md
	ldoc $^ >> code.md

.PHONY: all clean mac-time mac-ctr

mac-time: mosher
	instruments -t "Time Profiler" -D ~/Desktop/mosher.trace ./mosher

mac-ctr: mosher
	instruments -t "Core Data" -D ~/Desktop/mosher.trace ./mosher

tidy:
	@find | egrep "#" | xargs rm -f
	@find | egrep "\~" | xargs rm -f
	@find | egrep ".txt" | xargs rm -f

clean:
	rm -f $(EXE) *.o

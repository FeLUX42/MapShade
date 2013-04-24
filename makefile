CC=gcc
CFLAGS=-c -Wall
LDFLAGS=-lm -lpng -Wall
SOURCES=shademap.c
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=daylightmap

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
	
.PHONY: clean

clean:
	rm -f *.o *~ $(EXECUTABLE)




#$(ODIR)/%.o: %.c $(DEPS)
#	$(CC) -c -o $@ $< $(CFLAGS)

#hellomake: $(OBJ)
#	gcc -o $@ $^ $(CFLAGS) $(LIBS)

#.PHONY: clean

#clean:
#	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 

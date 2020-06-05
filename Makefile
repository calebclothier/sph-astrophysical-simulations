CC         = g++
EXECUTABLE = main
CFLAGS     = -c -Wall -DUSING_OSX

SOURCES    = main.cpp PARTICLE_SYSTEM.cpp BH_TREE.cpp PARTICLE.cpp 
OBJECTS    = $(SOURCES:.cpp=.o)

$(EXECUTABLE): $(OBJECTS)
	$(CC) -framework GLUT -framework OpenGL -framework Cocoa $(OBJECTS) -o $(EXECUTABLE) -Wno-deprecated-declarations

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@ -Wno-deprecated-declarations

clean:
	rm -rf *.o $(EXECUTABLE)
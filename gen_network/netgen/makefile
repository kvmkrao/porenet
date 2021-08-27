CC 		= c++
CFLAGS        	= -c -O3
LDFLAGS 	= -lm
PROG		= netgen_linux_x86
SOURCES		= \
		node.cpp \
		rockElem.cpp \
		problem.cpp \
		main.cpp 

OBJECTS		= $(SOURCES:.cpp=.o)

all: $(SOURCES) $(PROG)
	
$(PROG): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
        
clean:
	rm -f netgen *.o 

        

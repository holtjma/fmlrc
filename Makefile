CC=g++
CFLAGS=-O3 -std=c++11 -Wall -pthread
LDFLAGS=-O3 -std=c++11 -Wall -pthread
SOURCES=alignment_util.cpp base_bwt.cpp bit_array.cpp csa_bwt.cpp file_iterators.cpp main.cpp rle_bwt.cpp string_util.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=fmlrc

CON_SOURCES=converter/converter_main.cpp
CON_OBJECTS=$(CON_SOURCES:.cpp=.o)
CON_EXEC=fmlrc-convert

all: $(EXECUTABLE) $(CON_EXEC)

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

$(CON_EXEC): $(CON_OBJECTS)
	$(CC) $(LDFLAGS) $(CON_OBJECTS) -o $@

.cpp.o:
	$(CC) -c $(CFLAGS) $< -o $@
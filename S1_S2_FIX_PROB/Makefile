HEADER1 = msp.hpp
HEADER2 = fun.hpp
TARGET = main

OBJS := $(patsubst %.cpp,%.o,$(wildcard *.cpp))
CC = mpic++
OPTION = -O3

$(TARGET): $(OBJS)
	$(CC) -o $(TARGET) $(OBJS) $(OPTION) -lm 

%.o: %.cpp $(HEADER1) $(HEADER2)
	$(CC) $(CFLAGS) -c $<

clean:
	rm -rf *.o

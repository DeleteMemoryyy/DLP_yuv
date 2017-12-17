CXX = g++
MAKE = make
CONVERT = Convert.cpp
EXE = convert

LIBS = -lGL `pkg-config --static --libs glfw3`
CXXFLAGS = -I UI_LIB/ -I UI_LIB/glfw/include/ `pkg-config --cflags glfw3`
CXXFLAGS += -Wall -Wformat -g -mmmx -msse2 -mavx -mavx2 -ffloat-store

all:$(CONVERT) Convert.h
	$(CXX) -o $(EXE) $(OBJS) $(CONVERT) $(CXXFLAGS) $(LIBS)

.cpp.o:
	$(CXX) $(CXXFLAGS) $(LIBS) -c -o $@ $<

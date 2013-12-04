#!/bin/make

TARGET = ./main$(EXEEXT)
#SRCS = main.cpp vector3.cpp screen.cpp color.cpp light.cpp photon.cpp sphere.cpp polygon3.cpp scene.cpp aabb3.cpp aabb3n.cpp
#SRCS = ulist.c hamming.c main.c
#SRCS = unit.c unit32.c hamming.c main.c
SRCS = main.cpp
CXX = g++
CXXFLAGS = -O3 -Wall

#OPENCVINC = `pkg-config --cflags opencv`
#OPENCVLIB = `pkg-config --libs opencv`

LDFLAGS  = -fopenmp

OBJS = $(SRCS:.cpp=.o)

.PHONY: all clean
#.SUFFIXES: .c .o
.SUFFIXES: .cpp .o


all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(OBJS) -o $(TARGET) $(LDFLAGS) $(OPENCVLIB)

.cpp.o:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OPENCVINC) -c $<

depend:
	$(CXX) -MM $(INCLUDE) $(CXXFLAGS) $(SRCS) > dependencies

clean:
	rm -rf $(OBJS) $(TARGET)


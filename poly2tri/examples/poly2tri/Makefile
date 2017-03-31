####### Compiler, tools and options

CC	=	gcc
CXX	=	g++
CFLAGS	=	-Wall -W -O2
CXXFLAGS=	-Wall -W -O2
INCPATH	=	
LINK	=	g++
LFLAGS	=	

TAR	=	tar -cf
GZIP	=	gzip -9f

####### Files

HEADERS =	defs.h \
		geometry.h \
		splay.h \
		parse.h
SOURCES =	geometry.cc \
		predicates.cc \
		main.cc
OBJECTS =	geometry.o \
		predicates.o \
		main.o
DIST	=	
TARGET	=	poly2tri

####### Implicit rules

.SUFFIXES: .cpp .cxx .cc .C .c

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $@ $<

.cxx.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $@ $<

.cc.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $@ $<

.C.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $@ $<

.c.o:
	$(CC) -c $(CFLAGS) $(INCPATH) -o $@ $<

####### Build rules


all: $(TARGET)

$(TARGET): $(OBJECTS) 
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJECTS)

dist:
	$(TAR) poly2tri.tar poly2tri.pro $(SOURCES) $(HEADERS) $(DIST)
	$(GZIP) poly2tri.tar

clean:
	-rm -f $(OBJECTS) $(TARGET)
	-rm -f *~ core

####### Compile
geometry.o: geometry.cc \
		geometry.h \
		defs.h \
		splay.h

predicates.o: predicates.cc

main.o: main.cc \
		geometry.h \
		defs.h \
		splay.h \
		parse.h


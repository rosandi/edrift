
SRC = buffer.cpp drift.cpp equ.cpp driftlib.cpp

LIB = libdrift
OBJ = $(SRC:.cpp=.o)
CXX=mpic++ -Wall
CXXFLAGS = -O -I../ -I/opt/ompi/include -I$(HOME)/include
LIBFLAGS = -L$(HOME)/lib
ARCHIVE = ar
ARCHFLAG = -rc
LINK = mpic++
LINKFLAGS =	-O
USRLIB =
SYSLIB =

all: 	$(OBJ)
	$(ARCHIVE) $(ARFLAGS) $(LIB).a $(OBJ)

shlib: $(OBJ)
	mpicxx -shared -o $(LIB).so $(LIBFLAGS) $(OBJ) -lm -lomdtool -lmuparser

%.o:%.cpp
	$(CXX) $(CXXFLAGS) -fPIC -c $<

clean:
	-rm -f *.o *~ $(LIB).a $(LIB).so

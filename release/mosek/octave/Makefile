PROJ=mosek

INCL = ..
LIBS = ..

CXXFLAGS=-Wall -Wshadow -std=c++98 -ggdb -O3 -fPIC
CXXFLAGS2=-Wall -std=c++98 -ggdb -O3 -fPIC

MKOCT=mkoctfile

all:
	@ echo "Starting ${PROJ} creation..."
	@ echo "Compiling project..."
	$(MKOCT) ${PROJ}.cpp -o __${PROJ}_qp__.oct -I$(INCL) -L$(LIBS) -lmosek64 && echo "Success"

clean:
	rm *.o *.so *.oct

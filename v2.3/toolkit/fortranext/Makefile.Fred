LD = g++
LDFLAGS += -L../fastnlotoolkit/.libs -Wl,-rpath,../fastnlotoolkit/.libs -lfastnlotoolkit
# C++ compiler settings
GXX = g++
CFLAGS += -I../fastnlotoolkit/include/fastnlotk/
# Fortran compiler settings
FC = gfortran
LDFLAGS += -lgfortran

all: fortran_example

fortran_example: fastNLOFI_cc.o fastNLOFI_f.o
	$(LD) $(LDFLAGS) -o fortran_example fastNLOFI_cc.o fastNLOFI_f.o

fastNLOFI_cc.o: fastNLOFI_cc.cc
	$(GXX) $(CFLAGS) -o fastNLOFI_cc.o -c fastNLOFI_cc.cc

fastNLOFI_f.o: fastNLOFI_f.f
	$(FC) -o fastNLOFI_f.o -c fastNLOFI_f.f

clean:
	rm -rf *.o fortran_example

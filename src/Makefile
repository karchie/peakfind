MEX=mex
CC=gcc
FC=gcc
FFLAGS=-fPIC -g --fixed-line-length-132 -fcray-pointer
CFLAGS=-std=c99 -fPIC -g
LDLIBS=-lblas -lm -lglfortran

.SUFFIXES: .mexa64

TARGETS= test_peakfind peakfind.so peakfind.mexa64

OBJECTS=find_peaks.o sphereblur.o pflog.o pferr.o \
	fftsol.o hsphere3d.o imgpad.o imgvalx.o matopr.o polfit.o

all: $(TARGETS)

clean:
	rm -f *.[oa] $(TARGETS)

test_peakfind: test_peakfind.c $(OBJECTS)
	$(CC) -g $(LDFLAGS) -o $@ $^ $(LDLIBS)

lib: peakfind.so

peakfind.so: $(OBJECTS)
	$(CC) -shared $^ -o peakfind.so $(LDLIBS)

find_peaks.o: find_peaks.c peakfind.h pf_util.h

sphereblur.o: sphereblur.c pf_util.h

pflog.o: pflog.c pf_util.h

pferr.o: pferr.c pf_util.h

%.mexa64 : %.c
	$(MEX) $(MEXFLAGS) $^ -o $@ $(LDFLAGS) $(LDLIBS)

peakfind.mexa64: peakfind.c $(OBJECTS)

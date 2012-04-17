MEX=mex
CC=gcc
FC=gcc
FFLAGS=-fPIC -g --fixed-line-length-132 -fcray-pointer
CFLAGS=-std=c99 -fPIC -g
LDLIBS=-lblas -lm -lgfortran

.SUFFIXES: .mex

TARGETS= test_peakfind peakfind.so \
	peakfind.mex sphereblur.mex pf_set_log.mex

OBJECTS=find_peaks.o blur.o pflog.o pferr.o \
	fftsol.o hsphere3d.o imgpad.o imgvalx.o matopr.o polfit.o

all: $(TARGETS)

clean:
	rm -f $(OBJECTS) $(TARGETS) *~

test_peakfind: test_peakfind.c $(OBJECTS)
	$(CC) -g $(LDFLAGS) -o $@ $^ $(LDLIBS)

lib: peakfind.so

peakfind.so: $(OBJECTS)
	$(CC) $(LDFLAGS) -shared $^ -o peakfind.so $(LDLIBS)

find_peaks.o: find_peaks.c peakfind.h pf_util.h

blur.o: blur.c pf_util.h

pflog.o: pflog.c pf_util.h

pferr.o: pferr.c pf_util.h

%.mex : %.c
	$(MEX) $(MEXFLAGS) $^ -o $* $(LDFLAGS) $(LDLIBS)

peakfind.mex: peakfind.c $(OBJECTS)

sphereblur.mex: sphereblur.c $(OBJECTS)

pf_set_log.mex: pf_set_log.c $(OBJECTS)
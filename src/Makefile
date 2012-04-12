CC=gcc
FC=gcc
FFLAGS=-fPIC -g --fixed-line-length-132 -fcray-pointer
CFLAGS=-fPIC -g
LDLIBS=-lblas -lm

TARGETS= test_peakfind peakfind.so

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

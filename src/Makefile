FC=gcc
FFLAGS=-fPIC -g --fixed-line-length-132 -fcray-pointer
CFLAGS=-fPIC -g
LDLIBS=-lblas -lm

all: test_peakfind peakfind.so

test_peakfind: test_peakfind.c peakfind.h peakfind.so
	$(CC) -g $(LDFLAGS) -o $@ $< peakfind.so $(LDLIBS)

lib: peakfind.so

peakfind.so: peakfind.o fftsol.o hsphere3d.o imgpad.o imgvalx.o \
 matopr.o polfit.o
	$(CC) -shared $^ -o peakfind.so $(LDLIBS)

peakfind.o: peakfind.c peakfind.h

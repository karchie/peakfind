Building this package is easy enough, once you have the tools.

On Mac OS X: 

Install gfortran: http://gcc.gnu.org/wiki/GFortranBinaries

You'll need to set your path so that the gfortran-savvy gcc comes
before the system installed gcc, so something like:

export PATH=/usr/local/gfortran/bin:$PATH

You'll also need Matlab (specifically, mex) in your path.

This still won't quite work when mex is doing the compiling. This
works for me:

LDFLAGS=-L/usr/local/gfortran/lib/gcc/x86_64-apple-darwin10/4.6.2 make

The makefile is still a little messed up on Mac -- the mex suffix is
wrong. This will be fixed someday.
# Compile this with:
# Note it doesn't really need g++ except that the library routines use
# C++ style comments and therefore need to be compiled with g++. To
# get the linker to work properly, this then has to be compiled with
# g++ as well!

g++ -o chisig chisig.c -I$HOME/libsrc/numerics -L$HOME/lib -lnumerics
g++ -o chitab chitab.c -I$HOME/libsrc/numerics -L$HOME/lib -lnumerics

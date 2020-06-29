atl = /usr/lib/x86_64-linux-gnu/atlas/
#atl = /data/prog/lib/atlas/build/lib/
cla = ~/clapack/CLAPACK-3.2.1/
inatl = -I /data/prog/lib/atlas/include

lib = ./library/
in = -I $(lib)
incla = -I $(cla)INCLUDE
#clapack = $(incla) $(inatl) $(cla)*.a  $(atl)libcblas.a $(atl)libatlas.a $(cla)F2CLIBS/libf2c.a -lm
clapack = $(incla) $(inatl) $(cla)*.a $(cla)F2CLIBS/libf2c.a -lm -lblas -llapack


all:
	g++ -O3 -o run-sa joint-state.cpp $(lib)analysis.cpp $(lib)model.cpp $(lib)observation.cpp $(lib)tensor.cpp $(lib)math_lib.cpp $(in) $(incla) $(inatl) $(clapack)

clean:
	rm *.o

#src = strategy-3
#all:
#	mpicxx -O3 -c $(src).cpp $(in)
#	mpicxx -O3 -c $(lib)analysis.cpp $(in)
#	mpicxx -O3 -c $(lib)model.cpp $(in)
#	mpicxx -O3 -c $(lib)observation.cpp $(in)
#	mpicxx -O3 -c $(lib)tensor.cpp $(in)
#	mpicxx -O3 -c $(lib)math_lib.cpp $(in) $(incla) $(inatl)
#	mpicxx -O3 -o run-SA-9 $(src).o analysis.o observation.o tensor.o math_lib.o model.o $(clapack)
#	rm *.o


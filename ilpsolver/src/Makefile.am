bin_PROGRAMS=ilpsolver

# GUROBI_HOME, BOOST_HOME required as environmental variables

ilpsolver_CPPFLAGS = -I/opt/gurobi/include
ilpsolver_LDFLAGS = -L/opt/gurobi/lib -Wl,-rpath,/opt/gurobi/lib
ilpsolver_LDADD = -lgurobi_c++ -lgurobi70

ilpsolver_SOURCES = ilpsolver.h ilpsolver.cc \
					config.h config.cc \
					main.cc

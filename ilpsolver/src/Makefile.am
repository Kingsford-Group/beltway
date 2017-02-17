bin_PROGRAMS=ilpsolver

# GUROBI_HOME, BOOST_HOME required as environmental variables

ilpsolver_CPPFLAGS = -I$(GUROBI_HOME)/include
ilpsolver_LDFLAGS = -L$(GUROBI_HOME)/lib
ilpsolver_LDADD = -lgurobi_c++ -lgurobi

ilpsolver_SOURCES = ilpsolver.h ilpsolver.cc \
					main.cc

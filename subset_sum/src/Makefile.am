bin_PROGRAMS=test_subset_sum subset_sum_ilp subset_sum_maxflow

BOOST_DIR=/opt/local/stow/boost-1.55.0-gcc

subset_sum_ilp_CPPFLAGS = -I/opt/local/stow/gurobi702/include
subset_sum_ilp_LDFLAGS = -L/opt/local/stow/gurobi702/lib -Wl,-rpath,/opt/local/stow/gurobi702/lib
subset_sum_ilp_LDADD = -lgurobi_c++ -lgurobi70

subset_sum_maxflow_CPPFLAGS = -I$(BOOST_DIR)/include
subset_sum_maxflow_LDFLAGS = -L$(BOOST_DIR)/lib -Wl,-rpath,$(BOOST_DIR)/lib
subset_sum_maxflow_LDADD = -lboost_graph

test_subset_sum_SOURCES = test_subset_sum.cc \
						  pointer.cc pointer.h \
						  subset_sum.cc subset_sum.h 

subset_sum_ilp_SOURCES = subset_sum_ilp_driver.cc \
						 subset_sum_ilp.cc subset_sum_ilp.h \
						 pointer.cc pointer.h \
						 subset_sum.cc subset_sum.h

subset_sum_maxflow_SOURCES = subset_sum_maxflow_driver.cc \
							 subset_sum_maxflow.cc subset_sum_maxflow.h \
							 pointer.cc pointer.h \
							 subset_sum.cc subset_sum.h \
							 mygraph.h mygraph.cc 

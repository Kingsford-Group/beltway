bin_PROGRAMS=test_subset_sum subset_sum_ilp subset_sum_maxflow

#if we use gurobi
subset_sum_ilp_CPPFLAGS = -I/opt/local/stow/gurobi702/include
subset_sum_ilp_LDFLAGS = -L/opt/local/stow/gurobi702/lib -Wl,-rpath,/opt/local/stow/gurobi702/lib
subset_sum_ilp_LDADD = -lgurobi_c++ -lgurobi70

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
							 subset_sum.cc subset_sum.h

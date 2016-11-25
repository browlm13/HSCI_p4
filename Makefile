###################################################################
#  Makefile for project 4
#
#  Daniel R. Reynolds
#  SMU Mathematics
#  Math 3316
#  31 October 2015
###################################################################

# compiler & flags
CXX = g++
CXXFLAGS = -O2 -std=c++11

# makefile targets
all : test_Gauss2.exe	test_int.exe	test_gen.exe	test_adapt.exe	test_carbon.exe

test_carbon.exe : more_general_composite_int.cpp adaptive_int.cpp test_carbon.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

test_gen.exe : test_gen.cpp more_general_composite_int.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

test_int.exe : test_int.cpp composite_int.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

test_adapt.exe : more_general_composite_int.cpp test_adapt.cpp adaptive_int.cpp 
	$(CXX) $(CXXFLAGS) $^ -o $@

test_Gauss2.exe : test_Gauss2.cpp composite_Gauss2.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

clean :
	\rm -f *.o *.txt

realclean : clean
	\rm -f *.exe *~


####### End of Makefile #######

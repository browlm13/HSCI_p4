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
all : test_Gauss2.exe	test_int.exe

test_int.exe : test_int.cpp composite_int.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

test_Gauss2.exe : test_Gauss2.cpp composite_Gauss2.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

clean :
	\rm -f *.o *.txt

realclean : clean
	\rm -f *.exe *~


####### End of Makefile #######

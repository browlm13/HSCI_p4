###################################################################
#  Makefile for project 4
#
#  Laurence Brown
#  SMU Mathematics
#  Math 3316
#  1 December 2016
###################################################################

# compiler & flags
CXX = g++
CXXFLAGS = -O2 -std=c++11

# makefile targets
all : test_Gauss2.exe	test_Gausx.exe

test_Gausx.exe : matrix.cpp test_Gaussx.cpp composite_Gaussx.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@


test_Gauss2.exe : test_Gauss2.cpp composite_Gauss2.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

clean :
	\rm -f *.o *.txt

realclean : clean
	\rm -f *.exe *~


####### End of Makefile #######

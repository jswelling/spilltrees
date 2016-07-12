
# CXXFLAGS= -I. -g -fPIC -DSWIGRUNTIME_DEBUG
CXXFLAGS= -I. -g -fPIC
LDFLAGS= -g
PYTHONINCCMD= python -c 'import sys,os; print os.path.join(sys.prefix,"include","python%s"%sys.version[:3])'


targets: spilltree_tester spilltree.py _spilltree.so 

spilltree_tester: spilltree_tester.o
	$(CXX) ${LDFLAGS} -o $@ spilltree_tester.o

spilltree_tester.o:: spilltree_tester.cc spilltree.h

swig: spilltree.h spilltree.i
	swig -python -c++ -o spilltree_wrap.cpp spilltree.i

spilltree.py: spilltree.h spilltree.i swig
	echo "swig has been run"

spilltree_wrap.cpp: spilltree.h spilltree.i swig
	echo "swig has been run"

spilltree_wrap.o: spilltree_wrap.cpp
	$(CXX) $(CXXFLAGS) -I`${PYTHONINCCMD}` -c spilltree_wrap.cpp

_spilltree.so: spilltree_wrap.o
	$(CXX) $(CXXFLAGS) -shared -o $@ spilltree_wrap.o

clean:
	rm -f *.o
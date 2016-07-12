%{
using namespace std;
#include <spilltree.h>
#ifdef DOUBLE
typedef double FloatType;
#else
typedef float FloatType;
#endif

#ifndef DIM
#define DIM 2
#endif
%}

#ifdef DOUBLE
typedef double FloatType;
#define PRECSTRING "double"
#else
typedef float FloatType;
#define PRECSTRING "float"
#endif
#ifndef DIM
#define DIM 2
#endif

%define DOCSTRING(mytype)
"This extension module provides support for the Hybrid Spill-Tree
nearest neighbor finding algorithm as described in 'An Investigation
of Practical Approximate Nearest Neighbor Algorithms' by Liu, Moore,
Gray and Yang.  The algorithm can be made exact (at some cost in
running time) by setting the constructor parameter tau equal to 0.0.
The source code is based on a C++ template which can provide trees
of arbitrary dimensionality.  This extension was configured to
use " ## mytype ## " precision internally."
%enddef

%module(docstring=DOCSTRING( PRECSTRING )) spilltree

%include "typemaps.i";
%feature("autodoc",1);

%exception {
   try {
      $action
   } catch (char* s) {
      fprintf(stderr,"Exception: %s\n",s);
      PyErr_SetString(PyExc_RuntimeError, s);
      return NULL;
   } catch (char const* s) {
      fprintf(stderr,"Exception: %s\n",s);
      PyErr_SetString(PyExc_RuntimeError, s);
      return NULL;
   }
}

%{
class _HookHolder: public SpTree<_HookHolder,DIM>::TreeEntry {
public:
  PyObject* hook;
  SpTree<_HookHolder,DIM>::Point p;
  _HookHolder();
  _HookHolder(PyObject* obj, FloatType loc[DIM]);
  virtual const SpTree<_HookHolder,DIM>::Point& getLoc() const  { return p; };
  PyObject* getHook() const { return hook; }
};

class JWrapNeighborInfo {
public:
  PyObject* obj;
  FloatType sepsqr;
  JWrapNeighborInfo() { obj= NULL; sepsqr= 0.0; }
  JWrapNeighborInfo( PyObject* obj_in, FloatType sepsqr_in )
  { obj= obj_in; sepsqr= sepsqr_in; }
};

class JWrapSpTree {
public: 
  JWrapSpTree( PyObject* input, 
	       const FloatType tau=0.3, const FloatType rho=0.7);
  ~JWrapSpTree();
  int getNDim();
  void dump(FILE* ofile);
  JWrapNeighborInfo findApproxNearest( const _HookHolder& p );
  void findApproxKNearest( JWrapNeighborInfo* neighbors_out, 
			   int* nFound_out,
			   const int k, const _HookHolder& p );
protected:
  SpTree<_HookHolder,DIM>* tree;
};

_HookHolder::_HookHolder()
{
  for (int i=0; i<DIM; i++) p.x[i]= 0.0;
  hook= NULL;
}

_HookHolder::_HookHolder( PyObject* obj, FloatType loc[DIM])
{
  for (int i=0; i<DIM; i++) p.x[i]= loc[i];
  hook= obj;
}

JWrapSpTree::JWrapSpTree( PyObject* input,
			  const FloatType tau, const FloatType rho )
{
  /* Check if it is a list */
  if (PyList_Check(input)) {
    int n = PyList_Size(input);
    if (n==0) throw "Can't build a spilltree with zero samples";
    _HookHolder** v = new _HookHolder*[n];
    for (int i = 0; i < n; i++) {
      void* hhPtr;
      int res1= SWIG_ConvertPtr(PyList_GetItem(input,i),
				&hhPtr, SWIGTYPE_p__HookHolder, 0 | 0);
      if (!SWIG_IsOK(res1)) {
	PyErr_SetString(PyExc_ValueError,"Input sequence element of inappropriate type");
	delete [] v;
	throw "Input sequence element of inappropriate type";
      }
      
      v[i]= reinterpret_cast<_HookHolder*>(hhPtr);
    }
    tree= new SpTree<_HookHolder,DIM>( v, n, tau, rho );
    delete [] v;
  } else {
    PyErr_SetString(PyExc_ValueError,"Expected a sequence");
    throw "Expected a sequence";
  }
}

JWrapSpTree::~JWrapSpTree()
{
  delete tree;
}

int JWrapSpTree::getNDim()
{
  return tree->getNDim();
}

void JWrapSpTree::dump(FILE* ofile)
{
  tree->dump(ofile);
}

JWrapNeighborInfo 
JWrapSpTree::findApproxNearest( const _HookHolder& p )
{
  SpTree<_HookHolder,DIM>::NeighborInfo* thisNbr= tree->findApproxNearest(p);
  JWrapNeighborInfo sni= JWrapNeighborInfo(thisNbr->nbr->getHook(), 
						     thisNbr->sepsqr);
  delete thisNbr;
  return sni;
}

void JWrapSpTree::findApproxKNearest(JWrapNeighborInfo* neighbors_out,
					  int* nFound_out,
					  const int k, const _HookHolder& p)
{
  SpTree<_HookHolder,DIM>::NeighborInfo* hooks_out=
    new SpTree<_HookHolder,DIM>::NeighborInfo[k];
  int nFound= tree->findApproxKNearest( hooks_out, p, k );
  for (int i=0; i<nFound; i++) 
    neighbors_out[i]= JWrapNeighborInfo(hooks_out[i].nbr->getHook(),
					     hooks_out[i].sepsqr);

  delete [] hooks_out;
  *nFound_out= nFound;
}

%}

%typemap(in) FloatType loc[DIM] (FloatType temp[DIM]) {
  int i;
  if (!PySequence_Check($input)) {
    PyErr_SetString(PyExc_ValueError,"Expected a sequence");
    return NULL;
  }
  if (PySequence_Length($input) != DIM) {
    PyErr_SetString(PyExc_ValueError,"Location array size mismatch");
    return NULL;
  }
  for (i = 0; i < DIM; i++) {
    PyObject *o = PySequence_GetItem($input,i);
    if (PyNumber_Check(o)) {
      temp[i] = (FloatType) PyFloat_AsDouble(o);
    } else {
      PyErr_SetString(PyExc_ValueError,"Sequence elements must be numbers");
      return NULL;
    }
  }
  $1 = temp;
}

class _HookHolder {
public:
  _HookHolder(PyObject* obj, FloatType loc[DIM]);
  PyObject* getHook();
};

%feature("pythonprepend") 
JWrapSpTree::JWrapSpTree(PyObject* input,
			 const FloatType tau=0.3, const FloatType rho=0.7 ) 
%{
        hookList= []
        for a in input:
            new_HookHolder= _HookHolder(a,a.getLoc())
            hookList.append(new_HookHolder)
        input = hookList
%}

%feature("pythonappend") 
JWrapSpTree::JWrapSpTree(PyObject* input,
			 const FloatType tau=0.3, const FloatType rho=0.7 ) 
%{
        self.hookList= hookList
%}

%feature("pythonappend")
JWrapSpTree::~JWrapSpTree()
%{
        self.hookList= None
%}

%typemap(in) (FILE* ofile) {
  if (PyFile_Check($input)) {
    $1= PyFile_AsFile($input);
  } else {
    PyErr_SetString(PyExc_ValueError,"Input must be a file");
    return NULL;
  }
}

%feature("pythonprepend") 
JWrapSpTree::findApproxNearest(const _HookHolder& p) 
%{
        p = _HookHolder(p, p.getLoc())
%}

%typemap(out) JWrapNeighborInfo {
  PyObject* o= $1.obj;
  FloatType sepsqr= $1.sepsqr;
  PyObject* tuple= PyTuple_New(2);
  PyTuple_SetItem(tuple,0,o);
  Py_INCREF(o); // because SetItem steals refs
  PyObject* o2= PyFloat_FromDouble(sepsqr);
  PyTuple_SetItem(tuple,1,o2);
  $result= tuple;
}

%feature("pythonprepend") 
JWrapSpTree::findApproxKNearest( JWrapNeighborInfo* neighbors_out,
				      int* nFound_out_p,
				      const int k, const _HookHolder& p ) 
%{
        p= _HookHolder(p,p.getLoc())
%}

%typemap(in, numinputs=1) (JWrapNeighborInfo* neighbors_out, 
			   int* nFound_out_p, const int k)

{
  // typemap(in, numinputs=1) on neighbors_out
  $3= (int)PyInt_AsLong($input);
  $2= new int;
  $1= new JWrapNeighborInfo[$3];
  // end typemap(in, numinputs=1) on neighbors_out
}

%typemap(in) const _HookHolder& p {
  // typemap(in) for _HookHolder& p
  void* hhPtr;
  int res1= SWIG_ConvertPtr($input, &hhPtr, SWIGTYPE_p__HookHolder, 0 | 0);
  if (!SWIG_IsOK(res1)) {
    SWIG_exception_fail(SWIG_ArgError(res1), "in '$symname', argument " "$argnum" " of type '$1_type'"); 
  }
  
  $1= reinterpret_cast<_HookHolder*>(hhPtr);
  // end typemap(in) for _HookHolder& p
}

%typemap(argout) (JWrapNeighborInfo* neighbors_out,
                  int* nFound_out_p, const int k)
{
  // typemap(argout) for neighbors_out
  JWrapNeighborInfo* neighbors_out= $1;
  
  PyObject* outList= PyList_New(*$2);

  for (int i=0; i<*$2; i++) {
    PyObject* tuple= PyTuple_New(2);
    Py_INCREF(neighbors_out[i].obj); // because SetItem steals refs
    PyTuple_SetItem(tuple,0,neighbors_out[i].obj);
    PyTuple_SetItem(tuple,1,PyFloat_FromDouble(neighbors_out[i].sepsqr));
    PyList_SetItem(outList,i,tuple);
  }
  // Clean up things allocated in the input typemap
  delete $2;
  delete [] $1;

  // Replace the original result with the list we've just built
  Py_XDECREF($result);
  $result= outList;
  // end typemap(argout) for neighbors_out
}

%rename (_NeighborInfo) JWrapNeighborInfo;
class JWrapNeighborInfo {
public:
  PyObject* obj;
  FloatType sepsqr;
  JWrapNeighborInfo( PyObject* obj_in, FloatType sepsqr_in )
  { obj= obj_in; sepsqr= sepsqr_in; }
};

%rename (SpTree) JWrapSpTree;
class JWrapSpTree {                              
public: 
  %feature(autodoc,"SpTree([ObjWithLoc, ObjWithLoc, ...], tau=0.3, rho=0.7) -> SpTree") JWrapSpTree;
  %feature(docstring,"""
The objects in the input list must be instances of classes having a
method:

  thisObjWithLoc.getLoc() -> [FloatType, FloatType... ]

where the coordinate list length matches the dimensionality.  The
parameters tau and rho have the meanings described in the paper
by Liu, Moore, Gray and Yang.  Set tau=0.0 for exact neighbors.
""") JWrapSpTree;
  JWrapSpTree( PyObject* input,
 	       const FloatType tau=0.3, const FloatType rho=0.7 );
  %feature("autodoc",1);
  ~JWrapSpTree();
  %feature(docstring,"Returns the dimensionality of the tree") getNDim;
  int getNDim();
  %feature(docstring,"Writes a summary of the tree structure to the given file") dump;
  void dump(FILE* ofile);
  %feature("autodoc","findApproxNearest(ObjWithLoc) -> (neighborObjWithLoc, sepsqr)") findApproxNearest;
  %feature("docstring","sepsqr is the square of the separation distance") findApproxNearest;
  JWrapNeighborInfo findApproxNearest( const _HookHolder& p );
  %feature("autodoc","findApproxKNearest(ObjWithLoc) -> [(obj,sepsqr),(obj,sepsqr),...]") findApproxKNearest;
  %feature("docstring","The entries in the output list are in order of increasing sepsqr") findApproxKNearest;
  void findApproxKNearest( JWrapNeighborInfo* neighbors_out, 
			   int* nFound_out_p,
			   const int k, const _HookHolder& p );
};



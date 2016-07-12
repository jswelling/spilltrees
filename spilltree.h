///////////////////////////////
// This software or portions of this software is Copyright © 2008, 
// Pittsburgh Supercomputing Center (PSC).
//
// Author: Joel Welling, PSC
// 
// Permission to use, copy, modify, and distribute this software and its 
// documentation without fee for personal  use or non-commercial use within 
// your organization is hereby granted, provided that the above copyright 
// notice is preserved in all copies and that the copyright and this 
// permission notice appear in supporting documentation.  PSC makes no 
// representations about the suitability of this software for any purpose.  
// It is provided "as is" without express or implied warranty.
////////////////////////////////

///////////////////////////////////////////////////////////////////////
// $Id: spilltree.h,v 1.15 2009-08-19 17:11:23 welling Exp $
///////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////
// This is an implementation of the Hybrid Spill-Tree nearest 
// neighbor finding algorithm as described in "An Investigation of
// Practical Approximate Nearest Neighbor Algorithms" by Liu, Moore, 
// Gray and Yang
///////////////////////////////////////////////////////////////////////

#include <math.h>
#include <assert.h>
#include <set>
#include <exception>
#include <algorithm>

//////////////
// Notes-
// -Expand the union to make the whole thing a bit smaller?  Or clearer?
//////////////

typedef float FloatType;

#define LEAF_LIMIT 5

template<class DerivedFmTreeEntry, int NDIM>
class SpTree {
 public:

  class Vector;

  class Point {
  public:
    FloatType x[NDIM];
    Point() { for (int i=0; i<NDIM; i++) x[i]= 0.0; }
    FloatType sep( Point& other );
    Vector operator-( const Point& other ) const {
      Vector result;
      for (int i=0; i<NDIM; i++) result.x[i]= this->x[i]-other.x[i];
      return result;
    }
    Point operator+( const Vector& other ) const {
      Point result;
      for (int i=0; i<NDIM; i++) result.x[i]= this->x[i] + other.x[i];
      return result;
    }
    FloatType sepSqr( const Point& other ) const {
      Vector v= *this - other;
      return v.dot(v);
    }
  };
  
  class Vector {
  public:
    FloatType x[NDIM];
    Vector() { for (int i=0; i<NDIM; i++) x[i]= 0.0; }
    Vector( const FloatType* v, const int n ) {
      if (n<NDIM) throw "not enough values to make a vector!";
      for (int i=0; i<NDIM; i++) x[i]=v[i];
    }
    Vector operator*( const FloatType scale ) const {
      Vector result;
      for (int i=0; i<NDIM; i++) result.x[i]= scale*this->x[i];
      return result;
    }
    FloatType dot( const Vector& other ) const {
      FloatType sum= this->x[0]*other.x[0];
      for (int i=1; i<NDIM; i++) sum += this->x[i]*other.x[i];
      return sum;
    }
  };

  // Classes to be inserted should derive from this interface
  class TreeEntry {
  public:
    virtual ~TreeEntry() {}
    static const int getNDim() { return NDIM; }
    virtual const Point& getLoc() const = 0;
  };

  class NeighborInfo {
  public:
    const DerivedFmTreeEntry* nbr;
    FloatType sepsqr;
    NeighborInfo() { sepsqr= 0.0; nbr= NULL; }
    NeighborInfo( const DerivedFmTreeEntry* nbr_in, FloatType sepsqr_in ) 
      { nbr= nbr_in; sepsqr= sepsqr_in; }
    bool operator<(const NeighborInfo& other) const
    { return (this->sepsqr < other.sepsqr); }
  };

 protected:
  Point A;
  Vector normVec;
  FloatType halfOverlap;
  int isLeaf;
  union {
    struct {
      long nLeafEntries;
      DerivedFmTreeEntry** leafEntries;
    } leafData;
    struct {
      SpTree* leftKid;
      SpTree* rightKid;
    } nodeData;
  } d;
  const DerivedFmTreeEntry* 
    findMostDistant(const DerivedFmTreeEntry& p1, 
		    DerivedFmTreeEntry const * const * candidates, 
		    int n)
  {
    if (n<=0) return NULL;
    
    const DerivedFmTreeEntry* farthest= candidates[0];
    FloatType dist= p1.getLoc().sepSqr(farthest->getLoc());
    for (int i=1; i<n; i++) {
      FloatType thisDist= p1.getLoc().sepSqr(candidates[i]->getLoc());
      if (thisDist>dist) {
	dist= thisDist;
	farthest= candidates[i];
      }
    }
    
    return farthest;
  }


  int selectIntoLeft(DerivedFmTreeEntry const * const * p, const int n, 
		     DerivedFmTreeEntry const ** kidVec)
  {
    int nKids= 0;
    for (int i=0; i<n; i++) {
      const DerivedFmTreeEntry* samp= p[i];
      Vector v= samp->getLoc()-A;
      FloatType proj= v.dot(normVec);
      if (proj <= halfOverlap) {
	kidVec[nKids++]= samp;
      }
    }
    return nKids;
  }

  int selectIntoRight(DerivedFmTreeEntry const * const * p, const int n, 
		      DerivedFmTreeEntry const ** kidVec)
  {
    int nKids= 0;
    for (int i=0; i<n; i++) {
      const DerivedFmTreeEntry* samp= p[i];
      Vector v= p[i]->getLoc()-A;
      FloatType proj= v.dot(normVec);
      if (proj >= -halfOverlap) kidVec[nKids++]= (DerivedFmTreeEntry*)samp;
    }
    return nKids;
  }

  int mergeSortUnique(NeighborInfo* out, const int nOut,
		      const NeighborInfo* buf1, const int nBuf1,
		      const NeighborInfo* buf2, const int nBuf2)
  {
    // sort in
    int i1= 0;
    int i2= 0;
    int iOut= 0;
    while (iOut<nOut) {
      if (i1<nBuf1) {
	if (i2<nBuf2) {
	  if (buf1[i1].nbr==buf2[i2].nbr)
	    i2++; // Skip; duplicate entry
	  else if (buf1[i1].sepsqr<=buf2[i2].sepsqr)
	    out[iOut++]= buf1[i1++];
	  else 
	    out[iOut++]= buf2[i2++];
	}
	else {
	  out[iOut++]= buf1[i1++];
	}
      }
      else {
	if (i2<nBuf2) {
	  out[iOut++]= buf2[i2++];
	}
	else {
	  break;
	}
      }
    }
    //   fprintf(stderr,"Merge sort! %d from %d %d\n",nOut, nBuf1, nBuf2);
    //   for (int i=0; i<nBuf1; i++) 
    //     fprintf(stderr,"buf1 %d: %g\n",i,buf1[i].sepsqr);
    //   for (int i=0; i<nBuf2; i++) 
    //     fprintf(stderr,"buf2 %d: %g\n",i,buf2[i].sepsqr);
    //   for (int i=0; i<nOut; i++) 
    //     fprintf(stderr,"out %d: %g\n",i,out[i].sepsqr);
    return iOut;
  }

  void printLeafData(FILE* ofile, int indent=0) const {
    for (int i=0; i<d.leafData.nLeafEntries; i++) { 
      Point p= d.leafData.leafEntries[i]->getLoc();
      fprintf(ofile,"%*s  ( ",indent,"");
      for (int j=0; j<NDIM; j++) {
	if (j==0) fprintf(ofile,"%g",p.x[j]);
	else fprintf(ofile,", %g",p.x[j]);
      }
      fprintf(ofile," )\n");
    }
  }

 public:

  SpTree() // produce an empty SpTree
    {
      // Creates an empty but internally consistent node
      halfOverlap= 0.0;
      isLeaf= 1; 
      d.leafData.nLeafEntries= 0;
      d.leafData.leafEntries= NULL;
    }

  SpTree(DerivedFmTreeEntry const * const * p, 
	 const int n, const FloatType tau_in= 0.3,
	 const FloatType rho_in= 0.7 )  throw (const char*) // build
    {
      if (n<=LEAF_LIMIT) {
	isLeaf= 1;
	halfOverlap= 0.0;
	d.leafData.nLeafEntries= n;
	d.leafData.leafEntries= new DerivedFmTreeEntry*[n];
	memcpy(d.leafData.leafEntries, p, n*sizeof(DerivedFmTreeEntry*));
      }
      else {
	isLeaf= 0;
	const DerivedFmTreeEntry* lc= findMostDistant(*p[0], p, n);
	if (lc==NULL) throw("Tried to build SpTree with a single point");
	const DerivedFmTreeEntry* rc= findMostDistant(*lc, p, n);
	if (lc==rc) {
	  // throw("Only one point or all points in same location");
	  // No choice but to put them all in one place, but this is probably
	  // a case of bad (degenerate) input data.
	  isLeaf= 1;
	  halfOverlap= 0.0;
	  d.leafData.nLeafEntries= n;
	  d.leafData.leafEntries= new DerivedFmTreeEntry*[n];
	  memcpy(d.leafData.leafEntries, p, n*sizeof(DerivedFmTreeEntry*));
	  fprintf(stderr,"Probably degenerate input data:\n");
	  printLeafData(stderr);
	}
	else {
	  Vector u= rc->getLoc()-lc->getLoc();
	  Vector halfVec= u*0.5;
	  FloatType normSqr= halfVec.dot(halfVec);
	  if (normSqr==0.0)
	    throw "Only one point, or all points are at same location";
	  normVec= halfVec*(1.0/sqrt(normSqr));
	  A= lc->getLoc() + halfVec;
	  halfOverlap= tau_in*halfVec.dot(normVec);
	  
	  DerivedFmTreeEntry const ** kidVec= new const DerivedFmTreeEntry*[n];
	  if (kidVec==NULL) {
	    throw("Out of memory allocating kidVec?");
	  }
	    
	  // We will build left and right kids on separate scans to
	  // reduce the amount of temporary memory that gets allocated.
	  // Start with the left kid.
	  //
	  // If the node is too imbalanced, switch to a metric tree
	  // by setting overlap to zero.  This is what makes this a
	  // hybrid Sp-Tree.
	  
	  int nKids= selectIntoLeft(p, n, kidVec);
	  if ((nKids>rho_in*n) || (n-nKids>rho_in*n)) {
	    halfOverlap= 0.0;
	    nKids= selectIntoLeft(p,n,kidVec);
	    // Test for a pathological case- points are far enough
	    // apart that a separation vector can be found, but when
	    // they are sorted floating point rounding error puts them
	    // all on the same side.  If so, they all go in the leaf 
	    // node. 
	    if (nKids==n) {
	      // throw("Points are too close together to sort");
	      isLeaf= 1;
	      d.leafData.nLeafEntries= n;
	      d.leafData.leafEntries= new DerivedFmTreeEntry*[n];
	      memcpy(d.leafData.leafEntries, p, n*sizeof(DerivedFmTreeEntry*));
	      fprintf(stderr,"Probably degenerate input data:\n");
	      printLeafData(stderr);
	      delete kidVec;
	      return;
	    }
	  }
	  d.nodeData.leftKid= 
	    new SpTree<DerivedFmTreeEntry,NDIM>( kidVec, nKids, tau_in, 
						 rho_in );
	  
	  // And now the right kid
	  nKids= selectIntoRight(p, n, kidVec);
	  d.nodeData.rightKid= 
	    new SpTree<DerivedFmTreeEntry,NDIM>( kidVec, nKids, tau_in, 
						 rho_in );
	  
	  delete kidVec;
	}
      }
    }
  
  virtual ~SpTree()
    {
      if (isLeaf) { 
	delete [] d.leafData.leafEntries; 
      }
      else {
	delete d.nodeData.rightKid;
	delete d.nodeData.leftKid;
      }
    }

  int getNDim() { return NDIM; }

  void dump(FILE* ofile, int indent=0, int deltaIndent=4)
  {
    if (isLeaf) {
      fprintf(ofile,"%*sleaf: %ld nodes\n",
	      indent,"",d.leafData.nLeafEntries);
      //printLeafData(ofile, indent);
    }
    else {
      fprintf(ofile,"%*sSplit, halfOverlap= %g\n",
	      indent,"",halfOverlap);
      d.nodeData.leftKid->dump(ofile,indent+deltaIndent);
      d.nodeData.rightKid->dump(ofile,indent+deltaIndent);
    }
  }

  NeighborInfo* findApproxNearest(const DerivedFmTreeEntry& p)
  {
    const Point& loc= p.getLoc();
    if (isLeaf) {
      assert( d.leafData.nLeafEntries!=0 );
      // fprintf(stderr,"Sorting through %d others in leaf\n",d.leafData.nLeafEntries);
      DerivedFmTreeEntry* result= d.leafData.leafEntries[0];
      FloatType minSepSqr= result->getLoc().sepSqr(loc);
      for (int i=1; i<d.leafData.nLeafEntries; i++) {
	DerivedFmTreeEntry* te= d.leafData.leafEntries[i];
	FloatType sep= te->getLoc().sepSqr(loc);
	if (sep<minSepSqr) {
	  result= te;
	  minSepSqr= sep;
	}
      }
      return new NeighborInfo(result,minSepSqr);
    }
    else {
      Vector v= loc-A;
      FloatType proj= v.dot(normVec);
      NeighborInfo* candidate= NULL;
      // fprintf(stderr,"proj= %g -> %s kid\n",proj,(proj<0.0)?"left":"right");
      if (proj<0.0) candidate= d.nodeData.leftKid->findApproxNearest(p);
      else candidate= d.nodeData.rightKid->findApproxNearest(p);
      if (halfOverlap!=0.0) {
	// This is a spilltree cell; do defeatist search
	return candidate;
      }
      else {
	// This is a metric cell; check other branch
	if (proj*proj<candidate->sepsqr) {
	  NeighborInfo* otherCandidate;
	  if (proj<0.0) 
	    otherCandidate= d.nodeData.rightKid->findApproxNearest(p);
	  else 
	    otherCandidate= d.nodeData.leftKid->findApproxNearest(p);
	  if (candidate->sepsqr<=otherCandidate->sepsqr) {
	    delete otherCandidate;
	    return candidate;
	  }
	  else {
	    delete candidate;
	    return otherCandidate;
	  }
	}
	else return candidate;
      }
    }
  }
  
  // next returns nFound; output array must contain at least k elements
  int findApproxKNearest(NeighborInfo* neighbors, 
			 const DerivedFmTreeEntry& p, const int k)
  {
    const Point& loc= p.getLoc();
    if (isLeaf) {
      assert( d.leafData.nLeafEntries!=0 );
      // fprintf(stderr,"Sorting through %d others in leaf\n",d.leafData.nLeafEntries);
      if (d.leafData.nLeafEntries<=k) {
	// Copy in *all* this leaf's entries
	//fprintf(stderr,"Leaf copy of %d for quota %d\n",d.leafData.nLeafEntries,k);
	for (int i=0; i<d.leafData.nLeafEntries; i++) {
	  DerivedFmTreeEntry* te= d.leafData.leafEntries[i];
	  neighbors[i].nbr= te;
	  neighbors[i].sepsqr= te->getLoc().sepSqr(loc);
	}
	// The rest of the algorithm expects the list to be sorted
	std::sort(neighbors,neighbors+d.leafData.nLeafEntries);
	return d.leafData.nLeafEntries;
      }
      else {
	// Copy in only the k closest.  We use a temporary array 
	// allocated on the stack.
	//fprintf(stderr,"Sorting %d into %d spaces\n",d.leafData.nLeafEntries,k);
	NeighborInfo niTable[LEAF_LIMIT];
	for (int i=0; i<d.leafData.nLeafEntries; i++) {
	  DerivedFmTreeEntry* te= d.leafData.leafEntries[i];
	  niTable[i].nbr= te;
	  niTable[i].sepsqr= te->getLoc().sepSqr(loc);
	}
	// Sorting uses the < operator for the NeighborInfo class
	std::sort(niTable,niTable+d.leafData.nLeafEntries);
	for (int i=0; i<k; i++) neighbors[i]= niTable[i];
	return k;
      }
    }
    else {
      Vector v= loc-A;
      FloatType proj= v.dot(normVec);
      SpTree* closeKid= NULL;
      SpTree* farKid= NULL;
      if (proj<0.0) {
	closeKid= d.nodeData.leftKid;
	farKid= d.nodeData.rightKid;
      }
      else {
	closeKid= d.nodeData.rightKid;
	farKid= d.nodeData.leftKid;
      }
      int closeCount= closeKid->findApproxKNearest(neighbors, p, k);
      if (closeCount==k && halfOverlap!=0.0) {
	// This is a SpTree cell; use approx method and have faith.
	//fprintf(stderr,"spilltree approximate accept of %d\n",k);
	return closeCount;
      }
      else {
	// Check max separation against sep to other kid.
	// If others are possible, make a new buffer, get k from that side,
	// and merge sort them.
	FloatType farthestSepSqr= neighbors[closeCount-1].sepsqr;
	if (closeCount<k || proj*proj<farthestSepSqr) {
	  NeighborInfo* buf1= new NeighborInfo[closeCount];
	  NeighborInfo* buf2= new NeighborInfo[k];
	  int farCount= farKid->findApproxKNearest(buf2, p, k);
	  
	  // move original values out of the way
	  for (int i=0; i<closeCount; i++) buf1[i]= neighbors[i];
	  
	  // Merge unique entries
	  int iOut= mergeSortUnique( neighbors, k, 
				     buf1, closeCount, buf2, farCount );
	  
	  // Clean up and return
	  delete [] buf1;
	  delete [] buf2;
	  return iOut;
	}
	else {
	  // We haven't got any more!
	  //fprintf(stderr,"Other branch out of range!\n");
	  return closeCount;
	}
      }
    }
  }
};

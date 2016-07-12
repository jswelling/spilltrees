////////////////////////////////
// This software or portions of this software is Copyright © 2008, 
// Pittsburgh Supercomputing Center (PSC).
// 
// Permission to use, copy, modify, and distribute this software and its 
// documentation without fee for personal  use or non-commercial use within 
// your organization is hereby granted, provided that the above copyright 
// notice is preserved in all copies and that the copyright and this 
// permission notice appear in supporting documentation.  PSC makes no 
// representations about the suitability of this software for any purpose.  
// It is provided "as is" without express or implied warranty.
////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <algorithm>

using namespace std;
#include <spilltree.h>


static char rcsid[] = "$Id: spilltree_tester.cc,v 1.7 2008-10-09 21:52:06 welling Exp $";

#define DIM 3

class TestPt: public SpTree<TestPt,DIM>::TreeEntry {
public:
  SpTree<TestPt,DIM>::Point p;
  int index;
  static int currentIndex;
  TestPt() { index= currentIndex++; }
  virtual const SpTree<TestPt,DIM>::Point& getLoc() const  { return p; };
  static TestPt* randomLocation();
};

int TestPt::currentIndex= 0;

TestPt* 
TestPt::randomLocation() 
{
  TestPt* pt= new TestPt();
  for (int i=0; i<DIM; i++) pt->p.x[i]= (FloatType)drand48();
  return pt;
};

static SpTree<TestPt,DIM>::NeighborInfo* 
findNbrSequential(const TestPt& loc, TestPt** v, int npts)
{
  TestPt* best= v[0];
  SpTree<TestPt,DIM>::Point locPt= loc.getLoc();
  SpTree<TestPt,DIM>::Point p= best->getLoc();
  FloatType bestSepSqr= locPt.sepSqr(p);
  for (int i=1; i<npts; i++) {
    p= v[i]->getLoc();
    FloatType sep= locPt.sepSqr(p);
    if (sep<bestSepSqr) {
      best= v[i];
      bestSepSqr= sep;
    }
  }
  // fprintf(stderr,"Sequential search: points %d and %d have sepSqr %g\n",
  //  loc.index, best->index, bestSepSqr);
  return new SpTree<TestPt,DIM>::NeighborInfo(best,bestSepSqr);
}

static int 
findKNbrSequential( SpTree<TestPt,DIM>::NeighborInfo* seqNbrTbl, 
		    int k, const TestPt& pt, TestPt** v, int n )
{
  // We'll just sort it; this method doesn't have to be fast.
  SpTree<TestPt,DIM>::NeighborInfo* allPts= 
    new SpTree<TestPt,DIM>::NeighborInfo[n];
  SpTree<TestPt,DIM>::Point loc= pt.getLoc();
  for (int i=0; i<n; i++) {
    allPts[i].nbr= v[i];
    allPts[i].sepsqr= v[i]->getLoc().sepSqr(loc);
  }
  std::sort(allPts, allPts+n);
  int nOut= (k<n) ? k : n;
  for (int i=0; i<nOut; i++) seqNbrTbl[i]= allPts[i];
  return nOut;
}

static void 
printUsageAndExit(const char* prog)
{
  fprintf(stderr,"Usage: %s [-d] [-n nNbrs] [-t nTests] NumPoints\n",prog);
  fprintf(stderr,"     -d means dump the tree\n");
  fprintf(stderr,"     nNbrs==0 (the default) means find only nearest\n");
  fprintf(stderr,"     nTests (default 10) is number of trials\n");
  exit(-1);
}

int 
main(int argc, char* argv[])
{
  int dumpFlag= 0;
  int nNbrs= 0;
  int nTests= 10;
  int n= 0;
  int iArg;

  for (iArg=1; iArg<argc; iArg++) {
    if (argv[iArg][0]=='-') {
      if (!strcmp(argv[iArg],"-d")) dumpFlag= 1;
      else if (!strcmp(argv[iArg],"-n")) {
	if (iArg==argc) printUsageAndExit(argv[0]);
	nNbrs= atoi(argv[++iArg]); 
      }
      else if (!strcmp(argv[iArg],"-t")) {
	if (iArg==argc) printUsageAndExit(argv[0]);
	nTests= atoi(argv[++iArg]); 
      }
      else printUsageAndExit(argv[0]);
    }
    else {
      if (n==0) n= atoi(argv[iArg]);
      else break;
    }
  }

  if (iArg!=argc) printUsageAndExit(argv[0]);

  printf("%d trials on %d points; finding %d neighbors\n",
	 nTests,n,nNbrs);

  TestPt** v= new TestPt*[n];
  for (int i=0; i<n; i++) v[i]= TestPt::randomLocation();
  SpTree<TestPt,DIM>* tree= new SpTree<TestPt,DIM>(v,n,0.0);
  if (dumpFlag) tree->dump(stdout);

  if (nNbrs==0) {

    // test findApproxNearest
    for (int i=0; i<nTests; i++) {
      TestPt* myLocPtr= TestPt::randomLocation();
      SpTree<TestPt,DIM>::NeighborInfo* 
	seqNbr= findNbrSequential( *myLocPtr, v, n );
      SpTree<TestPt,DIM>::NeighborInfo* 
	treeNbr= tree->findApproxNearest(*myLocPtr);
      printf("TestPt %d (sepsqr %g) or %d (sepsqr %g) is closest to sample point %d\n",
	     ((TestPt*)(seqNbr->nbr))->index, seqNbr->sepsqr,
	     ((TestPt*)(treeNbr->nbr))->index, treeNbr->sepsqr,
	     myLocPtr->index);
      delete treeNbr;
      delete seqNbr;
      delete myLocPtr;
    }

  }
  else {

    // Test findApproxKNearest
    for (int i=0; i<nTests; i++) {
      TestPt* myLocPtr= TestPt::randomLocation();
      SpTree<TestPt,DIM>::NeighborInfo* 
	seqNbrTbl= new SpTree<TestPt,DIM>::NeighborInfo[nNbrs];
      SpTree<TestPt,DIM>::NeighborInfo* 
	treeNbrTbl= new SpTree<TestPt,DIM>::NeighborInfo[nNbrs];
      int nSeqNbrs= findKNbrSequential( seqNbrTbl, nNbrs, *myLocPtr, v, n );
      int nTreeNbrs= tree->findApproxKNearest(treeNbrTbl, *myLocPtr, nNbrs);
      
      printf("###########################################\n");
      printf("Test point %d has ID %d\n",i,myLocPtr->index);
      
      printf("\tSeq.ID\tSeq.sepsqr\tTree.ID\tTree.sepsqr\n");
      int count= 0;
      for (int i=0; i<nNbrs; i++) {
	printf("%d\t",i);
	if (i<nSeqNbrs) printf("%d\t%g\t",
			       ((TestPt*)(seqNbrTbl[i].nbr))->index,
			       seqNbrTbl[i].sepsqr);
	else printf("   \t   \t");
	if (i<nTreeNbrs) printf("%d\t%g\n",
				((TestPt*)(treeNbrTbl[i].nbr))->index,
				treeNbrTbl[i].sepsqr);
	else printf("\n");
      }
      delete treeNbrTbl;
      delete seqNbrTbl;
      delete myLocPtr;
    }
  }

  delete tree;
  for (int i=0; i<n; i++) delete v[i];
  delete [] v;

}

// This code is part of the project "Ligra: A Lightweight Graph Processing
// Framework for Shared Memory", presented at Principles and Practice of 
// Parallel Programming, 2013.
// Copyright (c) 2013 Julian Shun and Guy Blelloch
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#ifndef LIGRA_H
#define LIGRA_H
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cstring>
#include <string>
#include <algorithm>
#include "parallel.h"
#include "gettime.h"
#include "utils.h"
#include "vertex.h"
#include "compressedVertex.h"
#include "vertexSubset.h"
#include "graph.h"
#include "IO.h"
#include "parseCommandLine.h"
#include "gettime.h"
#include "config.h"
#include "store.h"
using namespace std;

#include <src/object_store/RDMAObjectStore.h>
#ifdef SIRIUS
Store my_store;
#endif

#include "store.h"

//*****START FRAMEWORK*****

//options to edgeMap for different versions of dense edgeMap (default is DENSE)
enum options { DENSE, DENSE_FORWARD };

template <class vertex, class F>
bool* edgeMapDense(graph<vertex> GA, bool* vertexSubset, F &f, bool parallel = 0) {
  cout << "edgeMapDense" << endl;
  long numVertices = GA.n;
#ifndef SIRIUS
  vertex *G = GA.V;
#endif
  bool* next = newA(bool,numVertices);
  {parallel_for (long i=0; i<numVertices; i++) {
    next[i] = 0;
    if (f.cond(i)) {
#ifdef SIRIUS
        // I don't know if this thing changes the vertex or not
        // worst case it does
      //vertex* v = objStore[i + 1];
      vertex* v = reinterpret_cast<vertex*>(getStore(i + 1));
      v->decodeInNghBreakEarly(i, vertexSubset, f, next, parallel);
      //objStore.put(&v, sizeof(v), i + 1);
#else
      G[i].decodeInNghBreakEarly(i, vertexSubset, f, next, parallel);
#endif
    }
  }}
  return next;
}

template <class vertex, class F>
bool* edgeMapDenseForward(graph<vertex> GA, bool* vertexSubset, F &f) {
  cout << "edgeMapDenseForward" << endl;
  long numVertices = GA.n;
  vertex *G = GA.V;
  bool* next = newA(bool,numVertices);
  {parallel_for(long i=0;i<numVertices;i++) next[i] = 0;}
  {parallel_for (long i=0; i<numVertices; i++){
    if (vertexSubset[i]) {
      G[i].decodeOutNgh(i, vertexSubset, f, next);
    }
  }}
  return next;
}

template <class vertex, class F>
pair<long,uintE*> edgeMapSparse(vertex* frontierVertices, uintE* indices, 
        uintT* degrees, uintT m, F &f, 
        long remDups=0, uintE* flags=NULL) {
  uintT* offsets = degrees;
  long outEdgeCount = sequence::plusScan(offsets, degrees, m);
  uintE* outEdges = newA(uintE,outEdgeCount);
  {parallel_for (long i = 0; i < m; i++) {
      uintT v = indices[i], o = offsets[i];
      vertex vert = frontierVertices[i]; 
      vert.decodeOutNghSparse(v, o, f, outEdges);
    }}
  uintE* nextIndices = newA(uintE, outEdgeCount);
  if(remDups) remDuplicates(outEdges,flags,outEdgeCount,remDups);
  // Filter out the empty slots (marked with -1)
  long nextM = sequence::filter(outEdges,nextIndices,outEdgeCount,nonMaxF());
  free(outEdges);
  return pair<long,uintE*>(nextM, nextIndices);
}

// decides on sparse or dense base on number of nonzeros in the active vertices
template <class vertex, class F>
vertexSubset edgeMap(graph<vertex> GA, vertexSubset &V, F f, intT threshold = -1, 
		 char option=DENSE, bool remDups=false) {
  long numVertices = GA.n, numEdges = GA.m;
  if(threshold == -1) threshold = numEdges/20; //default threshold
#ifndef SIRIUS
  vertex *G = GA.V;
#endif
  long m = V.numNonzeros();
  if (numVertices != V.numRows()) {
    cout << "edgeMap: Sizes Don't match" << endl;
    abort();
  }
  // used to generate nonzero indices to get degrees
  uintT* degrees = newA(uintT, m);
  vertex* frontierVertices;
  V.toSparse();
  frontierVertices = newA(vertex,m);
  {parallel_for (long i=0; i < m; i++){
  //  cout << "Get from store i: " << i << endl;
#ifdef SIRIUS
    vertex v = *reinterpret_cast<vertex*>(getStore(V.s[i] + 1));
#else
    vertex v = G[V.s[i]];
#endif
    degrees[i] = v.getOutDegree();
    frontierVertices[i] = v;
    }}
  cout << "degrees and frontierVertices done" << endl;
  uintT outDegrees = sequence::plusReduce(degrees, m);
  if (outDegrees == 0) return vertexSubset(numVertices);
  if (m + outDegrees > threshold) { 
    cout << "edgeMapDense" << endl;
    V.toDense();
    free(degrees);
    free(frontierVertices);
    bool* R = (option == DENSE_FORWARD) ? 
      edgeMapDenseForward(GA,V.d,f) : 
      edgeMapDense(GA, V.d, f, option);
    vertexSubset v1 = vertexSubset(numVertices, R);
    return v1;
  } else { 
    cout << "edgeMapSparse" << endl;
    pair<long,uintE*> R = 
      remDups ? 
      edgeMapSparse(frontierVertices, V.s, degrees, V.numNonzeros(), f, 
		    numVertices, GA.flags) :
      edgeMapSparse(frontierVertices, V.s, degrees, V.numNonzeros(), f);
    free(degrees);
    free(frontierVertices);
    return vertexSubset(numVertices, R.first, R.second);
  }
}

template<class F>
vertexSubset  edgeMap2(graph<symmetricVertex> GA, 
        vertexSubset &V, F f, intT threshold = -1,
        char option=DENSE, bool remDups=false) {
  long numVertices = GA.n, numEdges = GA.m;
  if(threshold == -1) threshold = numEdges/20; //default threshold
#ifndef SIRIUS
  symmetricVertex *G = GA.V;
#endif
  long m = V.numNonzeros();
  if (numVertices != V.numRows()) {
    cout << "edgeMap: Sizes Don't match" << endl;
    abort();
  }
  // used to generate nonzero indices to get degrees
  uintT* degrees = newA(uintT, m);
  symmetricVertex* frontierVertices;
  V.toSparse();
  frontierVertices = newA(symmetricVertex,m);
  {parallel_for (long i=0; i < m; i++){
  //  cout << "Get from store i: " << i << endl;
#ifdef SIRIUS
    symmetricVertex v = *reinterpret_cast<symmetricVertex*>(getStore(V.s[i] + 1));
#else
    symmetricVertex v = G[V.s[i]];
#endif
    degrees[i] = v.getOutDegree();
    frontierVertices[i] = v;
    }}
  cout << "degrees and frontierVertices done" << endl;
  uintT outDegrees = sequence::plusReduce(degrees, m);
  if (outDegrees == 0) return vertexSubset(numVertices);
  if (m + outDegrees > threshold) { 
    cout << "edgeMapDense" << endl;
    V.toDense();
    free(degrees);
    free(frontierVertices);
    bool* R = (option == DENSE_FORWARD) ? 
      edgeMapDenseForward(GA,V.d,f) : 
      edgeMapDense(GA, V.d, f, option);
    vertexSubset v1 = vertexSubset(numVertices, R);
    return v1;
  } else { 
    cout << "edgeMapSparse" << endl;
    pair<long,uintE*> R = 
      remDups ? 
      edgeMapSparse(frontierVertices, V.s, degrees, V.numNonzeros(), f, 
		    numVertices, GA.flags) :
      edgeMapSparse(frontierVertices, V.s, degrees, V.numNonzeros(), f);
    free(degrees);
    free(frontierVertices);
    return vertexSubset(numVertices, R.first, R.second);
  }
}

//*****VERTEX FUNCTIONS*****

//Note: this is the optimized version of vertexMap which does not
//perform a filter
template <class F>
void vertexMap(vertexSubset V, F add) {
  long n = V.numRows(), m = V.numNonzeros();
  if(V.isDense) {
    {parallel_for(long i=0;i<n;i++)
	if(V.d[i]) add(i);}
  } else {
    {parallel_for(long i=0;i<m;i++)
	add(V.s[i]);}
  }
}

//Note: this is the version of vertexMap in which only a subset of the
//input vertexSubset is returned
template <class F>
vertexSubset vertexFilter(vertexSubset V, F filter) {
  long n = V.numRows(), m = V.numNonzeros();
  V.toDense();
  bool* d_out = newA(bool,n);
  {parallel_for(long i=0;i<n;i++) d_out[i] = 0;}
  {parallel_for(long i=0;i<n;i++)
      if(V.d[i]) d_out[i] = filter(i);}
  return vertexSubset(n,d_out);
}

//cond function that always returns true
inline bool cond_true (intT d) { return 1; }

template<class vertex>
void Compute(graph<vertex>&, commandLine);

int parallel_main(int argc, char* argv[]) {
#ifdef SIRIUS
    init_store();
#endif
  commandLine P(argc,argv," [-s] <inFile>");
  char* iFile = P.getArgument(0);
  bool symmetric = P.getOptionValue("-s");
  bool compressed = P.getOptionValue("-c");
  bool binary = P.getOptionValue("-b");
  long rounds = P.getOptionLongValue("-rounds",3);
  //if (compressed) {
  //  if (symmetric) {
  //      std::cout << "compressedSymmetricVertex" << std::endl;
  //    graph<compressedSymmetricVertex> G =
  //      readCompressedGraph<compressedSymmetricVertex>(iFile,symmetric); //symmetric graph
  //    Compute(G,P);
  //    for(int r=0;r<rounds;r++) {
  //      startTime();
  //      Compute(G,P);
  //      nextTime("Running time");
  //    }
  //    G.del();
  //  } else {
  //      std::cout << "compressedAsymmetricVertex" << std::endl;
  //    graph<compressedAsymmetricVertex> G =
  //      readCompressedGraph<compressedAsymmetricVertex>(iFile,symmetric); //asymmetric graph
  //    Compute(G,P);
  //    if(G.transposed) G.transpose();
  //    for(int r=0;r<rounds;r++) {
  //      startTime();
  //      Compute(G,P);
  //      nextTime("Running time");
  //      if(G.transposed) G.transpose();
  //    }
  //    G.del();
  //  }
  //} else {
  //  if (symmetric) {
        std::cout << "symmetricVertex" << std::endl;
      graph<symmetricVertex> G =
        readGraph<symmetricVertex>(iFile,compressed,symmetric,binary); //symmetric graph
      Compute(G,P);
      for(int r=0;r<rounds;r++) {
        startTime();
        Compute(G,P);
        nextTime("Running time");
      }
      G.del();
    //} else {
    //    std::cout << "asymmetricVertex" << std::endl;
    //  graph<asymmetricVertex> G =
    //    readGraph<asymmetricVertex>(iFile,compressed,symmetric,binary); //asymmetric graph
    //  Compute(G,P);
    //  if(G.transposed) G.transpose();
    //  for(int r=0;r<rounds;r++) {
    //    startTime();
    //    Compute(G,P);
    //    nextTime("Running time");
    //    if(G.transposed) G.transpose();
    //  }
    //  G.del();
  //  }
  //}
}
#endif

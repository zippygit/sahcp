// Copyright (c) 2014 Argonne National Laboratory. All rights reserved.

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <numeric>

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::pair;
using std::string;
using std::ifstream;
using std::istream;
using std::stringstream;

// Forward declarations
void loadInsertedCycle( vector<vector<int> >& cycles_, int n_, ifstream& stdoutFile_ );
void loadUniqueFoundCycles( vector<vector<int> >& cycles_, int n_, pair<int,int>& findings_, ifstream& outyFile_ );
void countOverlapsWithInsertedCycle( vector<vector<int> >& cycles_, int n_, vector<int>& counts_ );

template< class T >
inline std::ostream& operator<<( std::ostream& os_, const std::vector<T>& v_ )
{
  os_ << "[";
  for ( typename vector<T>::const_iterator ii = v_.begin(); ii != v_.end(); ++ii ) { os_ << " " << *ii; }
  os_ << " ]";
  return os_;
}

int main( int argc, char* argv[] )
{
  int N = atoi( argv[ 1 ] );
  string stdoutFile = argv[ 2 ];
  string outyFile = argv[ 3 ];

  ifstream stdoutIn( stdoutFile.c_str(), ifstream::in );
  ifstream outyIn( outyFile.c_str(), ifstream::in );

  int nDistincCycles = 1;

  vector<vector<int> > cycles;
  loadInsertedCycle( cycles, N, stdoutIn );

  pair<int,int> findings; // first: how many times inserted cycle found; second: how many other unique cycles found
  loadUniqueFoundCycles( cycles, N, findings, outyIn );
  cout << endl;
  cout << "Found " << findings.first << " unique HC ; found the inserted cycle " << findings.second << " times" << endl;

  vector<int> edgesAlsoInInsertedCycle( ( cycles.size() - 1 ), 0 ); // initialize counts to zero
  // cout << "DEBUG: cycles.size()=" << cycles.size() << " edgesAlsoInInsertedCycle.size()=" << edgesAlsoInInsertedCycle.size() << endl;
  countOverlapsWithInsertedCycle( cycles, N, edgesAlsoInInsertedCycle );
  cout << "Edges in unique cycles that are also in inserted cycle: " << edgesAlsoInInsertedCycle << endl;
  int totalEdgesAlsoInInsertedCycle = std::accumulate( edgesAlsoInInsertedCycle.begin(), edgesAlsoInInsertedCycle.end(), 0 );
  // cout << "Total edges from all unique cycles that are also in inserted cycle: " << totalEdgesAlsoInInsertedCycle << endl;
  int nTrials = cycles.size() - 1;
  cout << "Average edges from all unique cycles that are also in inserted cycle: " << totalEdgesAlsoInInsertedCycle/nTrials << endl; 
}

void loadInsertedCycle( vector<vector<int> >& cycles_, int n_, ifstream& stdoutIn_ )
{
  string line;
  string icString = "inserted cycle = [";
  while ( getline( stdoutIn_, line ) )
  {
    // cout << "DEBUG: in loadInsertedCycle(): line=" << line << endl;
    if ( line.find( icString ) != string::npos )
    {
      size_t startPos = line.find( "[" ) + 1;
      size_t endPos = line.find( "]" ) - 1;
      size_t len = endPos - startPos + 1;

      string insertedCycleValues = line.substr( startPos, len );
      // cout << "insertedCycleValues = " << insertedCycleValues << endl;
      stringstream ss( insertedCycleValues );
      vector<int> insertedCycle;
      int vertex;
      while ( ss >> vertex ) insertedCycle.push_back( vertex );
      vector<int>::iterator minLoc = std::min_element( insertedCycle.begin(), insertedCycle.end() );
      // cout << "insertedCycle = " << insertedCycle << endl;
      std::rotate( insertedCycle.begin(), minLoc, insertedCycle.end() );
      cout << "reordered insertedCycle = " << insertedCycle << endl;
      cycles_.push_back( insertedCycle );
      return;
    }
  }
  cerr << "Missing inserted cycle in stdout file" << endl;
  exit(1);
}

void loadUniqueFoundCycles( vector<vector<int> >& cycles_, int n_, pair<int,int>& findings_, 
                            ifstream& outyIn_ )
{
  int nInsertedCycleFound = 0;
  int nUniqueCyclesFound = 0;
  bool firstTimeFindingInsertedCycle = true;
  string line;
  string hcString = "Found HC = [ ";
  while ( getline( outyIn_, line ) )
  {
    // cout << "DEBUG: in loadUniqueFoundCycles(): line=" << line << endl;
    if ( line.find( hcString ) != string::npos )
    {
      size_t startPos = line.find( "[" ) + 1;
      size_t endPos = line.find( "]" ) - 1;
      size_t len = endPos - startPos + 1;

      string hcValues = line.substr( startPos, len );
      // cout << "hcValues = " << hcValues << endl;
      stringstream ss( hcValues );
      vector<int> hc;
      int vertex;
      while ( ss >> vertex ) hc.push_back( vertex );
      vector<int>::iterator minLoc = std::min_element( hc.begin(), hc.end() );
      cout << "DEBUG: hc = " << hc << endl;
      std::rotate( hc.begin(), minLoc, hc.end() );
      bool uniqueHC = true;
      for ( vector<vector<int> >::iterator it = cycles_.begin(); it != cycles_.end(); it++ ) {
        if ( hc == *it ) {
          uniqueHC = false;
          if ( it == cycles_.begin() ) { nInsertedCycleFound += 1; }
          if ( firstTimeFindingInsertedCycle ) {
            nUniqueCyclesFound++;
            firstTimeFindingInsertedCycle = false;
          }
          break;
        }
      }
      if ( uniqueHC ) {
        cout << "reordered found hc = " << hc << endl;
        cycles_.push_back( hc );
        nUniqueCyclesFound++;
      }
    }
  }
  // int nCyclesFound = cycles_.size() - 1;
  // findings_.second = nInsertedCycleFound;
  findings_.first = nUniqueCyclesFound;
  findings_.second = nInsertedCycleFound;
  return;
}

void countOverlapsWithInsertedCycle( vector<vector<int> >& cycles_, int n_, vector<int>& counts_ )
{
  vector<int>& insertedCycle = cycles_[0];
  int nUniqueCyclesNonInserted = cycles_.size() - 1;
  int node, nextNode, insertedNode, nextInsertedNode, insertedNodeIndex, nextInsertedNodeIndex;
  bool edgeInInsertedCycle = false;
  for ( int i = 1; i <= nUniqueCyclesNonInserted; i++ ) {
    vector<int>& cycle = cycles_[i];
    for ( int j = 0; j < n_; j++ ) {
      node = cycle[j];
      nextNode = cycle[ ( j + 1 ) % n_ ];
      // cout << "DEBUG: i=" << i << " j=" << j << " node=" << node << " nextNode=" << nextNode << endl;
      vector<int>::iterator placeInInsertedCycle = std::find( insertedCycle.begin(), insertedCycle.end(), node );
      if ( placeInInsertedCycle != insertedCycle.end() ) {
        insertedNodeIndex = placeInInsertedCycle - insertedCycle.begin();
        nextInsertedNodeIndex = ( insertedNodeIndex + 1 ) % n_;
        insertedNode = insertedCycle[ insertedNodeIndex ];
        nextInsertedNode = insertedCycle[ nextInsertedNodeIndex ];
        // cout << "DEBUG: insertedNodeIndex=" << insertedNodeIndex << " nextInsertedNodeIndex=" << nextInsertedNodeIndex
        //      << " insertedNode=" << insertedNode << " nextInsertedNode=" << nextInsertedNode << endl;
        if ( nextInsertedNode == nextNode ) {
          counts_[ i - 1 ]++;
          // cout << "DEBUG: counts_[i-1 = " << (i-1) << "] = " << counts_[ i - 1 ] << endl;
        }
      }
    }
  }
  return;
}

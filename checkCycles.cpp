// Copyright (c) 2014 Argonne National Laboratory. All rights reserved.
// cd /projects/catalyst/zippy/hcp/mpi/ ; g++ -std=c++0x -o checkCycles checkCycles.cpp

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <numeric>
#include <iomanip>
#include <regex>

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::pair;
using std::string;
using std::ifstream;
using std::istream;
using std::stringstream;
using std::regex;
using std::sregex_token_iterator;

typedef vector<pair<long long int, long long int> > CompactGraph_t;

// Forward declarations
void loadInsertedCycle( vector<vector<int> >& cycles_, int n_, ifstream& stdoutFile_ );
void loadUniqueFoundCycles( vector<vector<int> >& cycles_, int n_, pair<int,int>& findings_, ifstream& outyFile_ );
void countOverlapsWithInsertedCycle( vector<vector<int> >& cycles_, int n_, vector<int>& counts_ );
void loadRandomTour( vector<vector<int> >& cycles_, int n_, ifstream& stdoutIn_ );
void countOverlapsWithRandomTour( vector<vector<int> >& cycles_, int n_, vector<int>& counts_ );
void loadRandomEdges( CompactGraph_t& randomEdges_, int n_, ifstream& stdoutIn_ );
void countOverlapsWithRandomEdges( CompactGraph_t& randomEdges_, vector<vector<int> >& cycles_, int n_, vector<int>& counts_ );

template< class T >
inline std::ostream& operator<<( std::ostream& os_, const std::pair<T,T>& p_ )
{
  os_ << "( " << p_.first << ", " << p_.second << " )";
  return os_;
}

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

  countOverlapsWithInsertedCycle( cycles, N, edgesAlsoInInsertedCycle );
  cout << "Edges in unique cycles that are also in inserted cycle: " << edgesAlsoInInsertedCycle << endl;

  int nUnique = findings.first;
  int nUniqueNotInserted = findings.first;
  if ( findings.second > 0 ) { nUniqueNotInserted--; } // don't include inserted cycle in average calculation
  int totalEdgesAlsoInInsertedCycle;
  totalEdgesAlsoInInsertedCycle = std::accumulate( edgesAlsoInInsertedCycle.begin(), edgesAlsoInInsertedCycle.end(), 0 );
  double averageEdgesAlsoInInsertedCycle;
  if ( ( nUnique == 1 ) && ( findings.second >= 1 ) ) {
    averageEdgesAlsoInInsertedCycle = (double) N;
  } else {
    averageEdgesAlsoInInsertedCycle = totalEdgesAlsoInInsertedCycle/( (double) nUniqueNotInserted );
  }
  cout << "Average edges from all unique cycles that are also in inserted cycle: " << std::fixed << std::setprecision(2) << averageEdgesAlsoInInsertedCycle << endl;    

  int minEdgesAlsoInInsertedCycle = N;
  int maxEdgesAlsoInInsertedCycle = N;
  if ( cycles.size() != 1 ) {
    minEdgesAlsoInInsertedCycle = *std::min_element( edgesAlsoInInsertedCycle.begin(), edgesAlsoInInsertedCycle.end() );
    maxEdgesAlsoInInsertedCycle = *std::max_element( edgesAlsoInInsertedCycle.begin(), edgesAlsoInInsertedCycle.end() );
  }
  cout << "Range of edges from all unique cycles that are also in inserted cycle: [ " << minEdgesAlsoInInsertedCycle
       << " - " << maxEdgesAlsoInInsertedCycle << " ]" << endl;

  // Do the same calculation comparing against a random tour (I know, I know, source reuse...)
  stdoutIn.clear();    // Rewind, step 1
  stdoutIn.seekg( 0 ); // Rewind, step 2
  loadRandomTour( cycles, N, stdoutIn );
  vector<int> edgesAlsoInRandomTour( ( cycles.size() - 1 ), 0 ); // initialize counts to zero
  countOverlapsWithRandomTour( cycles, N, edgesAlsoInRandomTour );
  cout << "Edges in unique cycles that are also in random tour: " << edgesAlsoInRandomTour << endl;

  int totalEdgesAlsoInRandomTour;
  totalEdgesAlsoInRandomTour = std::accumulate( edgesAlsoInRandomTour.begin(), edgesAlsoInRandomTour.end(), 0 );
  double averageEdgesAlsoInRandomTour;
  if ( ( nUnique == 1 ) && ( findings.second >= 1 ) ) {
    averageEdgesAlsoInRandomTour = (double) N;
  } else {
    averageEdgesAlsoInRandomTour = totalEdgesAlsoInRandomTour/( (double) nUniqueNotInserted);
  }
  cout << "Average edges from all unique cycles that are also in random tour: " << std::fixed << std::setprecision(3) << averageEdgesAlsoInRandomTour << endl;    

  int minEdgesAlsoInRandomTour = N;
  int maxEdgesAlsoInRandomTour = N;
  if ( cycles.size() != 1 ) {
    minEdgesAlsoInRandomTour = *std::min_element( edgesAlsoInRandomTour.begin(), edgesAlsoInRandomTour.end() );
    maxEdgesAlsoInRandomTour = *std::max_element( edgesAlsoInRandomTour.begin(), edgesAlsoInRandomTour.end() );
  }
  cout << "Range of edges from all unique cycles that are also in random tour: [ " << minEdgesAlsoInRandomTour
       << " - " << maxEdgesAlsoInRandomTour << " ]" << endl;

  // Do similar calculation comparing against a random set of N_m +1 edges
  stdoutIn.clear();    // Rewind, step 1
  stdoutIn.seekg( 0 ); // Rewind, step 2

  //zippyNope vector<CompactGraph_t> randomEdges;
  CompactGraph_t randomEdges;

  loadRandomEdges( randomEdges, N, stdoutIn );
  vector<int> edgesAlsoInRandomEdges( ( cycles.size() - 1 ), 0 ); // initialize counts to zero
  countOverlapsWithRandomEdges( randomEdges, cycles, N, edgesAlsoInRandomEdges );
  cout << "Edges in unique cycles that are also in random edges: " << edgesAlsoInRandomEdges << endl;

  int totalEdgesAlsoInRandomEdges;
  totalEdgesAlsoInRandomEdges = std::accumulate( edgesAlsoInRandomEdges.begin(), edgesAlsoInRandomEdges.end(), 0 );
  double averageEdgesAlsoInRandomEdges;
  if ( ( nUnique == 1 ) && ( findings.second >= 1 ) ) {
    averageEdgesAlsoInRandomEdges = (double) N;
  } else {
    averageEdgesAlsoInRandomEdges = totalEdgesAlsoInRandomEdges/( (double) nUniqueNotInserted);
  }
  cout << "Average edges from all unique cycles that are also in random edges: " << std::fixed << std::setprecision(3) << averageEdgesAlsoInRandomEdges << endl;    

  int minEdgesAlsoInRandomEdges = N;
  int maxEdgesAlsoInRandomEdges = N;
  if ( cycles.size() != 1 ) {
    minEdgesAlsoInRandomEdges = *std::min_element( edgesAlsoInRandomEdges.begin(), edgesAlsoInRandomEdges.end() );
    maxEdgesAlsoInRandomEdges = *std::max_element( edgesAlsoInRandomEdges.begin(), edgesAlsoInRandomEdges.end() );
  }
  cout << "Range of edges from all unique cycles that are also in random edges: [ " << minEdgesAlsoInRandomEdges
       << " - " << maxEdgesAlsoInRandomEdges << " ]" << endl;

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

void loadRandomTour( vector<vector<int> >& cycles_, int n_, ifstream& stdoutIn_ )
{
  string line;
  string icString = "random tour = [";
  while ( getline( stdoutIn_, line ) )
  {
    if ( line.find( icString ) != string::npos )
    {
      size_t startPos = line.find( "[" ) + 1;
      size_t endPos = line.find( "]" ) - 1;
      size_t len = endPos - startPos + 1;
      string randomTourValues = line.substr( startPos, len );
      stringstream ss( randomTourValues );

      vector<int> randomTour;
      int vertex;
      while ( ss >> vertex ) randomTour.push_back( vertex );
      vector<int>::iterator minLoc = std::min_element( randomTour.begin(), randomTour.end() );
      std::rotate( randomTour.begin(), minLoc, randomTour.end() );
      cout << "reordered randomTour = " << randomTour << endl;
      cycles_[ 0 ] = randomTour;
      return;
    }
  }
  cerr << "Missing random tour in stdout file" << endl;
  exit(1);
}

void loadRandomEdges( CompactGraph_t& randomEdges_, int n_, ifstream& stdoutIn_ )
{
  string line;
  string icString = "random edges = [";
  while ( getline( stdoutIn_, line ) )
  {
    if ( line.find( icString ) != string::npos )
    {
      size_t startPos = line.find( "{" ) + 1;
      size_t endPos = line.find( "}" ) - 1;
      size_t len = endPos - startPos + 1;
      string randomEdgesValues = line.substr( startPos, len );
      stringstream ssRandomEdgeValues( randomEdgesValues );
      long long int startNode, endNode;
      while ( ssRandomEdgeValues >> startNode >> endNode ) {
        randomEdges_.push_back( pair<long long int, long long int>( startNode, endNode ) );
      }
      cout << "DEBUG: randomEdges_ = " << randomEdges_ << endl;
      return;
    }
  }
  cerr << "Missing random edges in stdout file" << endl;
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
      // cout << "DEBUG: hc = " << hc << endl;
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
        cout << "DEBUG: rotated found unique hc = " << hc << endl;
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

void countOverlapsWithRandomTour( vector<vector<int> >& cycles_, int n_, vector<int>& counts_ )
{
  vector<int>& insertedCycle = cycles_[0];
  int nUniqueCyclesNonInserted = cycles_.size() - 1;
  int node, nextNode, insertedNode, nextInsertedNode, insertedNodeIndex, nextInsertedNodeIndex;
  bool edgeInRandomTour = false;
  for ( int i = 1; i <= nUniqueCyclesNonInserted; i++ ) {
    vector<int>& cycle = cycles_[i];
    for ( int j = 0; j < n_; j++ ) {
      node = cycle[j];
      nextNode = cycle[ ( j + 1 ) % n_ ];
      vector<int>::iterator placeInRandomTour = std::find( insertedCycle.begin(), insertedCycle.end(), node );
      if ( placeInRandomTour != insertedCycle.end() ) {
        insertedNodeIndex = placeInRandomTour - insertedCycle.begin();
        nextInsertedNodeIndex = ( insertedNodeIndex + 1 ) % n_;
        insertedNode = insertedCycle[ insertedNodeIndex ];
        nextInsertedNode = insertedCycle[ nextInsertedNodeIndex ];
        if ( nextInsertedNode == nextNode ) {
          counts_[ i - 1 ]++;
        }
      }
    }
  }
  return;
}

void countOverlapsWithRandomEdges( CompactGraph_t& randomEdges_, vector<vector<int> >& cycles_, int n_, vector<int>& counts_ )
{
  typedef pair<long long int, long long int> edge_t;
  int nUniqueCyclesNonInserted = cycles_.size() - 1;
  int node, nextNode, insertedNode, nextInsertedNode, insertedNodeIndex, nextInsertedNodeIndex;
  for ( int i = 1; i <= nUniqueCyclesNonInserted; i++ ) {
    vector<int>& cycle = cycles_[i];
    for ( int j = 0; j < n_; j++ ) {
      node = cycle[j];
      nextNode = cycle[ ( j + 1 ) % n_ ];
      edge_t edge( node, nextNode );
      CompactGraph_t::iterator placeInRandomEdges = std::find( randomEdges_.begin(), randomEdges_.end(), edge );
      if ( placeInRandomEdges != randomEdges_.end() ) {
        counts_[ i - 1 ]++;
      }
    }
  }
  return;
}

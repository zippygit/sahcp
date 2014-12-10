// Copyright (c) 2014 Argonne National Laboratory. All rights reserved.

// cd /projects/catalyst/zippy/hcp/mpi/ ; mpic++11 -g -i8 -O3 -ffast-math -o xm2rchcp m2rchcp.cpp;

// MPI version with two levels of parallelism (trials and permutations)
// This one repeats the annealing schedule with same setup but different RN
// sequence when the moves and annealing start.

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <algorithm>
#include <math.h>
#include <random>
#include <mpi.h>
#include <iomanip>
#ifdef ZIPPY_USE_BOOST
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::pair;
using std::random_shuffle;
using std::string;
using std::ifstream;
using std::ofstream;
using std::ostream;
using std::stringstream;
using std::ostringstream;
using std::uniform_real_distribution;
using std::default_random_engine;
using std::minstd_rand;
using std::uniform_int_distribution;


template< class T >
inline std::ostream& operator<<( std::ostream& os_, const std::vector<T>& v_ )
{
  os_ << "[";
  for ( typename vector<T>::const_iterator ii = v_.begin(); ii != v_.end(); ++ii ) { os_ << " " << *ii; }
  os_ << " ]";
  return os_;
}

typedef vector<long long> Tour_t;

class HCPGraph
{
public:

  HCPGraph( long long N_ ) : N_m( N_ )
  {
    N2_m = N_m * N_m;
    graph_m = new double[ N_ * N_ ];
    for ( long long i = 0; i < ( N_ * N_ ); i++ ) { graph_m[ i ] = 1.0; }
  }
  virtual ~HCPGraph() { delete [] graph_m; }
  
  void addRandomEdges( long long nEdges_, long long seed_, double weight_ )
  {
    //trySomethingDifferentInDesperation std::seed_seq seq{1 + seed_, 2 + seed_, 3 + seed_, 4 + seed_, 5 + seed_, N_m + seed_ };
    //trySomethingDifferentInDesperation std::mt19937 gen(seq);
    //zippyMake64bit std::mt19937 gen(seed_);
    std::mt19937_64 gen(seed_);
    //zippyMake64bit std::uniform_int_distribution<> dis( 0, N_m - 1 );
    std::uniform_int_distribution<long long int> dis( 0, N_m - 1 );

    long long startNode, endNode;
    bool acceptedEdge;
    for ( long long edge = 0; edge < nEdges_; edge++ ) {
      acceptedEdge = false;
      while ( !acceptedEdge ) {
        startNode = dis( gen );
        do { endNode = dis( gen ); } while ( endNode == startNode );
        if ( ! hasEdge( startNode, endNode ) ) {
          addEdge( startNode, endNode, weight_ );
          acceptedEdge = true;
        }
      }
    }
  }

  void addNonRandomEdges( long long nEdges_, double weight_, long long myRank_ )
  {
    long long edgesAdded = 0;
    for ( long long stride = 1; ( ( stride < ( N_m - 1 ) ) && ( edgesAdded < nEdges_ ) ); stride++ ) {
      for ( long long startNode = 0; ( ( startNode <= ( N_m - 1 - stride ) ) && ( edgesAdded < nEdges_ ) ); startNode += ( stride + 1 ) ) {
        long long endNode = startNode + stride;
        if ( !hasEdge( startNode, endNode ) && ( edgesAdded < nEdges_ ) ) {
          addEdge( startNode, endNode, weight_ );
          edgesAdded += 1;
        }
      }
    }
    if ( myRank_ == 0 ) cout << "rank0 anre(): added " << edgesAdded << " edges; nEdges_="<< nEdges_ << endl;
  }

  void insertCycle( long long cycleSeed_, long long myRank_ )
  {
    Tour_t tour( N_m );
    for ( long long ii = 0; ii < N_m; ii++ ) {
      tour[ ii ] = ii;
    }
    long long seedForCycle = ( cycleSeed_ + 1 ) * 3;
    //zippyMake64bit std::mt19937 shuffleGen( cycleSeed_ );
    std::mt19937_64 shuffleGen( cycleSeed_ );
    shuffle( &tour[ 0 ], &tour[ N_m ], shuffleGen );

    for ( long long ii = 0; ii < N_m - 1; ii++ ) {
      addEdge( tour[ii], tour[ ii + 1 ], 0.0 );
    }
    addEdge( tour[ N_m - 1 ], tour[ 0 ], 0.0 );
    if ( myRank_ == 0 ) { cout << "debug rank " << myRank_ << ":: inserted cycle = " << tour; }
  }

  void writeNRandomEdges( long long seed_, std::ostream& outy_, long long myRank_ )
  {
    std::seed_seq seq{1 + seed_, 2 + seed_, 3 + seed_, 4 + seed_, 5 + seed_, N_m + seed_ };
    //zippyMake64bit std::mt19937 gen(seq);
    std::mt19937_64 gen(seq);
    //zippyMake64bit std::uniform_int_distribution<> dis( 0, N_m - 1 );
    std::uniform_int_distribution<long long int> dis( 0, N_m - 1 );

    // ??? how to do this?
    // long long startNode, endNode;
    // bool acceptedEdge;
    // for ( long long edge = 0; edge < N_m; edge++ ) {
    //   acceptedEdge = false;
    //   while ( !acceptedEdge ) {
    //     startNode = dis( gen );
    //     do { endNode = dis( gen ); } while ( endNode == startNode );
    //     if ( hasEdge( startNode, endNode ) ) {
    //       cout << "rank " << myRank_ << ":: N random edges - "
    //       acceptedEdge = true;
    //     }
    //   }
    // }
  }

  void insertFakeEdges( long long nFakeEdges_, double fakeEdgeWeight_, long long randomGraphSeed_ )
  {
    long long seedForFakes = ( randomGraphSeed_ + 1 ) * 2;
    addRandomEdges( nFakeEdges_, seedForFakes, fakeEdgeWeight_ );
  }
  
  inline bool hasEdge( long long node1_, long long node2_ )
  {
    return ( graph_m[ node1_ * N_m + node2_ ] != 1.0 );
  }

  inline double getWeight( long long node1_, long long node2_ )
  {
    return graph_m[ node1_ * N_m + node2_ ];
  }

  inline long long nNodes() { return N_m; }

  inline long long nEdges()
  {
    long long nE = 0;
    for ( long long ii = 0; ii < N_m; ii++ ) {
      for ( long long jj = 0; jj < N_m; jj++ ) {
        if ( hasEdge( ii, jj ) ) { nE += 1; }
      }
    }
    return nE;
  }

  inline void addEdge( long long node1_, long long node2_, double weight_ )
  {
    graph_m[ node1_ * N_m + node2_ ] = weight_;
  }

  string toString() {
    stringstream ss;
    ss << "edges: [";
    for ( long long ii = 0; ii < N_m; ii++ ) {
      for ( long long jj = 0; jj < N_m; jj++ ) {
        if ( graph_m[ ii * N_m + jj ] != 1.0 ) {
          ss << "(" << ii << ", " << jj << ", weight: " << graph_m[ ii * N_m + jj ] << "), ";
        }
      }
    }
    return ss.str();
  }

private:
  double* graph_m;
  long long N_m;
  long long N2_m;
  bool seeded_m = false;
};

class HCPParameters
{
public:
  HCPParameters( ) { }

  HCPParameters( string& inputFileName_ ) : fileName( inputFileName_ ) { parseInput(); }

  HCPParameters( string fileName_, long long nMultipliers_, double minMultiplier_, double maxMultiplier_, 
                 long long nTrials_, long long firstTrial_, double nFakeEdgesMultiplier_, double fakeEdgeWeight_, 
                 double tFactor_, double t0_, long long naSteps_, double kMultiplier_, long long maxPermMultiplier_, 
                 bool retry_, bool randomEdges_, bool insertRandomCycle_, long long randomEdgesSeed_ )
    : fileName( fileName_ ), nMultipliers( nMultipliers_ ), minMultiplier( minMultiplier_ ), 
      maxMultiplier( maxMultiplier_ ), nTrials( nTrials_ ), firstTrial( firstTrial_ ), 
      nFakeEdgesMultiplier( nFakeEdgesMultiplier_ ), fakeEdgeWeight( fakeEdgeWeight_ ), 
      tFactor( tFactor_ ), t0( t0_ ), naSteps( naSteps_ ), kMultiplier( kMultiplier_ ), 
      maxPermMultiplier( maxPermMultiplier_ ), retry( retry_ ), randomEdges( randomEdges_ ),
      insertRandomCycle( insertRandomCycle_ ), randomEdgesSeed( randomEdgesSeed_) { }

  void parseInput()
    {
      string key, value, line;
      std::ifstream inpStream( ( fileName ) );
      while ( getline( inpStream, key, '=' ) ) {
        key.erase( remove_if( key.begin(), key.end(), isspace ), key.end());
        stringstream ssKey( key );
        ssKey >> key;
        getline( inpStream, value );
        stringstream ssValue( value );
        if ( key == "nMultipliers" ) {
          ssValue >> nMultipliers;
        } else if ( key == "minMultiplier" ) {
          ssValue >> minMultiplier;
        } else if ( key == "maxMultiplier" ) {
          ssValue >> maxMultiplier;
        } else if ( key == "nTrials" ) {
          ssValue >> nTrials;
        } else if ( key == "firstTrial" ) {
          ssValue >> firstTrial;
        } else if ( key == "nFakeEdgesMultiplier" ) {
          ssValue >> nFakeEdgesMultiplier;
        } else if ( key == "fakeEdgeWeight" ) {
          ssValue >> fakeEdgeWeight;
        } else if ( key == "" ) {
          ssValue >> fakeEdgeWeight;
        } else if ( key == "fakeEdgeWeight" ) {
          ssValue >> fakeEdgeWeight;
        } else if ( key == "tFactor" ) {
          ssValue >> tFactor;
        } else if ( key == "t0" ) {
          ssValue >> t0;
        } else if ( key == "naSteps" ) {
          ssValue >> naSteps;
        } else if ( key == "kMultiplier" ) {
          ssValue >> kMultiplier;
        } else if ( key == "maxPermMultiplier" ) {
          ssValue >> maxPermMultiplier;
        } else if ( key == "retry" ) {
          ssValue >> retry;
          cout << "streamed ssValue into retry, now retry = " << retry << endl;
        } else if ( key == "randomEdges" ) {
          ssValue >> randomEdges;
        } else if ( key == "insertRandomCycle" ) {
          ssValue >> insertRandomCycle;
        } else if ( key == "randomEdgesSeed" ) {
          ssValue >> randomEdgesSeed;
        }
      }
    }

  string toString() {
    stringstream ss;
    ss << "fileName=" << fileName << " nMultipliers=" << nMultipliers
       << " minMultiplier=" << minMultiplier << " maxMultiplier=" << maxMultiplier
       << " nTrials=" << nTrials << " firstTrial=" << firstTrial
       << " nFakeEdgesMultiplier=" << nFakeEdgesMultiplier
       << " fakeEdgeWeight=" << fakeEdgeWeight << " tFactor=" << tFactor << " t0=" << t0
       << " naSteps=" << naSteps << " kMultiplier=" << kMultiplier
       << " maxPermMultiplier=" << maxPermMultiplier << " retry=" << retry
       << " randomEdges=" << randomEdges << " insertRandomCycle=" << insertRandomCycle
       << " randomEdgesSeed=" << randomEdgesSeed;
    return ss.str();
  }

  string fileName;
  long long nMultipliers;
  double minMultiplier;
  double maxMultiplier;
  long long nTrials = 1;
  long long firstTrial = 0;
  double nFakeEdgesMultiplier;
  double fakeEdgeWeight;
  double tFactor;
  double t0;
  long long naSteps;
  double kMultiplier;
  long long maxPermMultiplier;
  bool retry = false;
  bool randomEdges = true;
  bool insertRandomCycle = true;
  long long randomEdgesSeed = 0;
};

double costFunction( HCPGraph& g_, Tour_t& tour_, bool debug_ )
{
  long long N = g_.nNodes();
  long long Nm1 = N - 1;
  double cost = 0.0;
  for ( long long vert = 0; vert < Nm1; vert++ ) {
    cost += g_.getWeight( tour_[ vert ], tour_[ vert + 1 ] );
  }
  cost += g_.getWeight( tour_[ N - 1 ], tour_[ 0 ] );
  return cost;
}

void lkTransport( Tour_t& tour_, Tour_t& newTour_, long long L_,
                  uniform_int_distribution<long long int>& pickNodes_, minstd_rand& pickNodesGen_ )
{
  // Pick pair of locations (A,B) with B>A, a third location C, and insert the 
  // segment A->B after location C [Lin-Kernighan "transport" move from TSP]
  // L = len( tour_ )

  // Pick two distinct nodes at random:
  long long node1, node2;
  node1 = pickNodes_( pickNodesGen_ );
  do {
    node2 = pickNodes_( pickNodesGen_ );
  } while ( node2 == node1 );
  long long locAIndex = node1;
  long long locBIndex = node2;
  if ( locBIndex < locAIndex ) {
    long long tempIndex = locAIndex;
    locAIndex = locBIndex;
    locBIndex = tempIndex;
  }
  long long coinFlip = 0;
  while ( ( locBIndex - locAIndex ) >= ( L_ - 2 ) ) {
    if ( coinFlip == 0 ) {
      locAIndex += 1;
      coinFlip = 1;
    } else {
      locBIndex -= 1;
      coinFlip = 0;
    }
  }

  // Third location; moved segment will be inserted between this and the next location:
  long long locCIndex = pickNodes_( pickNodesGen_ );
  while ( ( locCIndex >= locAIndex ) and ( locCIndex <= locBIndex ) ) {
    locCIndex = pickNodes_( pickNodesGen_ );
  }
  if ( locCIndex > locBIndex ) {
    long long lSegment1 = ( locBIndex - locAIndex + 1 ); // length from A thru B
    long long lSegment2 = locCIndex - ( locBIndex + 1 ) + 1; // length from B+1 thru C : locCIndex - ( locBIndex + 1 ) + 1
    for ( long long n = 0; n < locAIndex; n++ ) { newTour_[ n ] = tour_[ n ]; }
    for ( long long n = locAIndex; n < ( locAIndex + lSegment2 ); n++ ) {
      newTour_[ n ] = tour_[ n + lSegment1 ];
    }
    for ( long long n = locAIndex + lSegment2; n < locAIndex + lSegment2 + lSegment1 ; n++ ) {
      newTour_[ n ] = tour_[ n - lSegment2 ];
    }
    for ( long long n = ( locCIndex + 1 ); n < L_; n++ ) {
      newTour_[ n ] = tour_[ n ];
    }
  } else { // locCIndex < locAIndex
    long long lSegment1 = locAIndex - 1 - ( locCIndex + 1 ) + 1; // length from C+1 through A-1
    long long lSegment2 = ( locBIndex - locAIndex + 1 ); // length from A thru B
    for ( long long n = 0; n <= locCIndex; n++ ) { newTour_[ n ] = tour_[ n ]; }
    for ( long long n = locCIndex + 1; n < ( locCIndex + 1 + lSegment2 ); n++ ) {
      newTour_[ n ] = tour_[ n + lSegment1 ];
    }
    for ( long long n = ( locCIndex + 1 + lSegment2 ); n < ( locCIndex + 1 + lSegment2 +lSegment1 ); n++ ) {
      newTour_[ n ] = tour_[ n - lSegment2 ];
    }
    for ( long long n = ( locBIndex + 1 ); n < L_; n++ ) {
      newTour_[ n ] = tour_[ n ];
    }
  }
}

void swap( Tour_t& tour_, Tour_t& newTour_, long long L_,
           uniform_int_distribution<long long int>& pickNodes_, minstd_rand& pickNodesGen_ )
{
  // Swap two locations in the tour:
  // Pick two (different) locations at random:
  long long node1, node2;
  node1 = pickNodes_( pickNodesGen_ );
  do {
    node2 = pickNodes_( pickNodesGen_ );
  } while ( node2 == node1 );
  // Swap the locations:
  newTour_[ node1 ] = tour_[ node2 ];
  newTour_[ node2 ] = tour_[ node1 ];
  return;
}

long long move( Tour_t& tour_, Tour_t& newTour_, long long permutation_, long long L_,
          uniform_int_distribution<long long int>& pickNodes_, minstd_rand& pickNodesGen_ )
{
  if ( ( permutation_ % 2 ) == 0 ) {
    lkTransport( tour_, newTour_, L_, pickNodes_, pickNodesGen_ );
    return 1;
  } else {
    swap( tour_, newTour_, L_, pickNodes_, pickNodesGen_ );
    return 2;
  }
  return 3;
}

void windUp( bool foundHC_, Tour_t& tour_, long long step_, double t_, long long permutation_, std::ostream& outy_, long long trial_,
             HCPGraph& g_ ) //zippy:added g_ for debugging
{
  if ( foundHC_ ) {
    outy_ << "trial " << trial_ << " step " << step_ << " t=" << t_ << " perm " << permutation_
          << " Found HC = " << tour_ << " DEBUG: costFunction(tour_)=" << costFunction( g_, tour_, false ) << endl;
    cout << "trial " << trial_ << " step " << step_ << " t=" << t_ << " perm " << permutation_
         << " Found HC = " << tour_ << " DEBUG: costFunction(tour_)=" << costFunction( g_, tour_, false ) << endl;
  } else {
    outy_ << "trial " << trial_ << " Found no HC." << endl;
  }
}

// fwd decl:
bool hcp( HCPGraph& g_, HCPParameters& config_, ostream& outy_, long long trial_, MPI_Comm& trialComm_ );

void sweep( long long N_, string& inputFileName_ )
{
  int mpiError = 0, nRanks = 0, myRank = 0;
  mpiError = MPI_Comm_size( MPI_COMM_WORLD, &nRanks );
  mpiError = MPI_Comm_rank( MPI_COMM_WORLD, &myRank );
  mpiError = MPI_Barrier( MPI_COMM_WORLD );

  time_t startPoint, startPointM, endPoint, endPointM;
  time( &startPoint );

  long long NLogN = (long long)( N_ * log( N_ ) );
  if ( myRank == 0 ) { cout << "N_ = " << N_ << " NLogN = " << NLogN << endl; }

  string outputFileName;
  long long maxPermMultiplier;
  long long maxPermutations;
  long long nMultipliers;
  double minMultiplier = 99.0;
  double maxMultiplier = 99.0;
  long long nTrials;
  long long firstTrial;
  long long nFakeEdges = 0;
  double nFakeEdgesMultiplier;
  double fakeEdgeWeight = 0.0;
  bool retry = 0;
  double tFactor = 99.0;
  double t0 = 99.0;
  long long naSteps;
  double kMultiplier = 99.0;
  bool randomEdges = 1;
  bool insertRandomCycle = 1;
  long long randomEdgesSeed = 0;

  ostream* outyPtr;

  // Read input parameters
  HCPParameters* configPtr;
  long long cStringLength;
  outputFileName = inputFileName_ + ".outy";
  outputFileName = "outy." + inputFileName_.substr( 7 ); // Assume input filename begins with "config."
  if ( myRank == 0 )
  {
    configPtr = new HCPParameters( inputFileName_ );
    HCPParameters& config = *configPtr;
    outyPtr = new ofstream( outputFileName, ofstream::out );
    ostream& outy = *outyPtr;
    outy << "N_ = " << N_ << " NLogN = " << NLogN << endl;
    outy << "config = " << config.toString() << endl;
    maxPermMultiplier = config.maxPermMultiplier;
    maxPermutations = (long long)( maxPermMultiplier * pow( N_, 2 ) );
    outy << "maxPermutations = " << maxPermutations << endl;

    nMultipliers = config.nMultipliers;
    minMultiplier = config.minMultiplier;
    maxMultiplier = config.maxMultiplier;
    nTrials = config.nTrials;
    firstTrial = config.firstTrial;
    nFakeEdges = 0;
    nFakeEdgesMultiplier = config.nFakeEdgesMultiplier;
    fakeEdgeWeight = config.fakeEdgeWeight;
    tFactor = config.tFactor;
    t0 = config.t0;
    naSteps = config.naSteps;
    kMultiplier = config.kMultiplier;
    retry = config.retry;
    randomEdges = config.randomEdges;
    insertRandomCycle = config.insertRandomCycle;
    randomEdgesSeed = config.randomEdgesSeed;

    MPI_Bcast( &nMultipliers, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD );
    MPI_Bcast( &minMultiplier, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &maxMultiplier, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &nTrials, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD );
    MPI_Bcast( &firstTrial, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD );
    MPI_Bcast( &nFakeEdgesMultiplier, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &fakeEdgeWeight, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    mpiError = MPI_Barrier( MPI_COMM_WORLD );
    MPI_Bcast( &t0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &tFactor, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &naSteps, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD );
    MPI_Bcast( &kMultiplier, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &maxPermMultiplier, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD );
    long long retryInt = 0;
    if ( retry ) { retryInt = 1; }
    MPI_Bcast( &retryInt, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD );
    long long randomEdgesInt = 1;
    if ( !randomEdges ) { randomEdgesInt = 0; }
    MPI_Bcast( &randomEdgesInt, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD );
    long long insertRandomCycleInt = 1;
    if ( !insertRandomCycle ) { insertRandomCycleInt = 0; }
    MPI_Bcast( &insertRandomCycleInt, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD );
    MPI_Bcast( &randomEdgesSeed, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD );
  }
  else
  {
    outyPtr = new ostringstream( "" );
    ostream& outy = *outyPtr;

    MPI_Bcast( &nMultipliers, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD );
    MPI_Bcast( &minMultiplier, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &maxMultiplier, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &nTrials, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD );
    MPI_Bcast( &firstTrial, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD );
    MPI_Bcast( &nFakeEdgesMultiplier, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &fakeEdgeWeight, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    mpiError = MPI_Barrier( MPI_COMM_WORLD );
    MPI_Bcast( &t0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &tFactor, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &naSteps, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD );
    MPI_Bcast( &kMultiplier, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( &maxPermMultiplier, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD );
    long long retryInt = 0;
    MPI_Bcast( &retryInt, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD );
    if ( retryInt == 1 ) { retry = true; }
    long long randomEdgesInt = 1;
    MPI_Bcast( &randomEdgesInt, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD );
    if ( randomEdgesInt == 0 ) { randomEdges = false; }
    long long insertRandomCycleInt = 1;
    MPI_Bcast( &insertRandomCycleInt, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD );
    if ( insertRandomCycleInt == 0 ) { insertRandomCycle = false; }
    MPI_Bcast( &randomEdgesSeed, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD );

    configPtr = new HCPParameters( inputFileName_, nMultipliers, minMultiplier, maxMultiplier, nTrials,
                                   firstTrial, nFakeEdgesMultiplier, fakeEdgeWeight, tFactor, t0, naSteps,
                                   kMultiplier, maxPermMultiplier, retry, randomEdges, insertRandomCycle,
                                   randomEdgesSeed );
  }

  mpiError = MPI_Barrier( MPI_COMM_WORLD );
  ostream& outy = *outyPtr;
  
  HCPParameters& config = *configPtr;
  fakeEdgeWeight = config.fakeEdgeWeight;

  if ( ( myRank == 0 ) || ( myRank == 1 ) ) {
    cout << "rank " << myRank << ":: config = " << config.toString() << endl;
  }
  mpiError = MPI_Barrier( MPI_COMM_WORLD );

  double deltaMultiplier = ( maxMultiplier - minMultiplier ) / ( 1.0 * nMultipliers );

  for ( long long m = 0; m < nMultipliers; m++ ) {
    time( &startPointM );
    double multiplier = minMultiplier + m * deltaMultiplier;
    long long M = (long long)( multiplier * NLogN );
    // if ( ( M % 2 ) == 0 ) { M++; } //zippy try always odd.
    // if ( M == 572 ) { M = 571; } // zippy hack
    //zippy hack:
    double dNLogN = N_ * log( N_ );
    double randRoads = multiplier * dNLogN;
    M = (long long)( randRoads );
    //zippy hack.
    // M++; //zippy hack
    // long long randomGraphSeed = multiplier;
    // long long randomGraphSeed = M;
    // long long randomGraphSeed = M + myRank; //hack debug
    //trydifferent20141126 long long randomGraphSeed = M + M * randomEdgesSeed;
    long long randomGraphSeed = randomEdgesSeed; //trydifferent20141126
    HCPGraph g( N_ );
    if ( randomEdges ) {
      g.addRandomEdges( M, randomGraphSeed, 0.0 );
    } else {
      g.addNonRandomEdges( M, 0.0, myRank );
    }
    if ( myRank == 0 ) {
      outy << "debug: before fakes/tour, g.nEdges() = " << g.nEdges() << endl;
    }
    randomGraphSeed *= 2;
    if ( insertRandomCycle ) {
      g.insertCycle( randomGraphSeed, myRank );
      if ( myRank == 0 ) {
        outy << "debug: after tour insert, g.nEdges() = " << g.nEdges() << endl;
      }
    } else {
      // This hack makes the checkCycle code still work when there's no actual cycle inserted:
      Tour_t bogusTour( N_ );
      for ( long long nn = 0; nn < N_; nn++ ) { bogusTour[ nn ] = 0; }
      if ( myRank == 0 ) { cout << "debug rank " << myRank << ":: inserted cycle = " << bogusTour; }
    }

    // long long nRandomEdgesSeed = randomGraphSeed + 222333444555;
    // g.writeNRandomEdges( nRandomEdgesSeed, outy, myRank );

    if ( nFakeEdgesMultiplier > 0.0 ) { nFakeEdges = (long long)( nFakeEdgesMultiplier * NLogN ); }
    if ( myRank == 0 ) {
      cout << "nFakeEdges =" << nFakeEdges << endl;
      outy << "nFakeEdges =" << nFakeEdges << endl;
    }
    randomGraphSeed *= 3;
    if ( nFakeEdges > 0 ) { g.insertFakeEdges( nFakeEdges, fakeEdgeWeight, randomGraphSeed ); }
    if ( myRank == 0 ) {
      outy << "debug: after fakes insert, g.nnEdges()=" << g.nEdges() << endl;
      cout << "M = " << M << " multiplier = " << std::fixed << std::setprecision(4) << multiplier << endl;
      outy << "M = " << M << " multiplier = " << multiplier << "g = " << g.toString() << endl;
    }
    bool foundHC = false;
    long long cyclesFound = 0;

    // Divide trials as evenly among MPI ranks, and run blocks of trials in
    // parallel. Assumed that things divide out evenly, so user beware.

    long long nRanksPerTrial = nRanks / nTrials;
    if ( nRanksPerTrial == 31 ) {
      cout << "nRanksPerTrial=" << nRanksPerTrial << " nRanks=" << nRanks << " nTrials=" << nTrials << endl;
    }
    MPI_Comm trialComm;
    long long color = myRank % nTrials;
    long long key = myRank / nTrials;
    MPI_Barrier( MPI_COMM_WORLD );
    mpiError = MPI_Comm_split( MPI_COMM_WORLD, color, key, &trialComm );
    int nRanksPerTrialTest = 0, myTrialRank = 0;
    mpiError = MPI_Comm_size( trialComm, &nRanksPerTrialTest );
    mpiError = MPI_Comm_rank( trialComm, &myTrialRank );
    if ( nRanksPerTrialTest != nRanksPerTrial ) {
      cerr << "Error in configuration. nRanksPerTrialTest=" << nRanksPerTrialTest
           << " but nRanksPerTrial=" << nRanksPerTrial << endl;
      MPI_Abort( MPI_COMM_WORLD, 3737 );
    }
    for ( long long rank = 0; rank < nRanks; rank++ ) {
      if ( myRank == rank ) {
        if ( ( myRank % 64 ) == 0 ) {
          cerr << "debug: myRank = " << myRank << " color = " << color << " key = " << key << " myTrialRank = " << myTrialRank
               << " nRanksPerTrial = " << nRanksPerTrial << endl;
        }
      }
      mpiError = MPI_Barrier( MPI_COMM_WORLD );
    }
    bool die = false;
    for ( long long rnk = 0; rnk < nRanks; rnk++ ) {
      if ( myRank == rnk ) {
        if ( myTrialRank != key ) {
          cerr << "Error in configuration. myTrialRank=" << myTrialRank << " but key=" << key << endl;
          die = true;
        }
      }
      MPI_Barrier( MPI_COMM_WORLD );
    }
    if ( die ) { MPI_Abort( MPI_COMM_WORLD, 3838 ); }

    long long myFirstTrial = color;
    long long myNTrials = 1;

    for ( long long trial = myFirstTrial; trial < ( myFirstTrial + myNTrials ); trial++ ) {
      foundHC = hcp( g, config, outy, trial, trialComm );
      if ( foundHC ) { cyclesFound += 1; }
      mpiError = MPI_Barrier( trialComm );
    }
    if ( myTrialRank == 0 ) {
      cout << "DEBUG RE-ATTEMPT: myRank=" << myRank << " back from hcp()" << endl;
    }
    mpiError = MPI_Barrier( MPI_COMM_WORLD );
    cout.flush();
    if ( myTrialRank == 0 ) {
      cout << "DEBUG RE-ATTEMPT: myRank=" << myRank << " past first global barrier" << endl;
    }
    cout.flush();

    // mpiError = MPI_Barrier( MPI_COMM_WORLD );
    // if ( myTrialRank == 0 ) {
    //   cout << "DEBUG RE-ATTEMPT: myRank=" << myRank << " past second global barrier" << endl;
    // }
    // cout.flush();

    // Sum up totals from all ranks:
    long long totalCyclesFoundInMyTrial = 0;
    mpiError = MPI_Allreduce( &cyclesFound, &totalCyclesFoundInMyTrial, 1, MPI_LONG_LONG, MPI_MAX, trialComm );
    if ( myTrialRank == 0 ) {
      cout << "DEBUG RE-ATTEMPT: myRank=" << myRank << " past totalCyclesFoundInMyTrial allreduce" << endl;
    }
    mpiError = MPI_Barrier( MPI_COMM_WORLD );
    cout.flush();
    if ( myTrialRank == 0 ) {
      cout << "DEBUG RE-ATTEMPT: myRank=" << myRank << " past third global barrier" << endl;
    }
    cout.flush();
    long long totalCyclesFound = 0;
    mpiError = MPI_Allreduce( &totalCyclesFoundInMyTrial, &totalCyclesFound, 1, MPI_LONG_LONG, MPI_SUM, 
                              MPI_COMM_WORLD );
    if ( myTrialRank == 0 ) {
      cout << "DEBUG RE-ATTEMPT: myRank=" << myRank << " past totalCyclesFound allreduce" << endl;
    }
    totalCyclesFound /= nRanksPerTrial;
    if ( myRank == 0 ) {
      cout << "rank 0:: cycles found = " << cyclesFound << endl;
      cout << "     ...totalCyclesFound =" << totalCyclesFound << endl;
    }

    mpiError = MPI_Barrier( MPI_COMM_WORLD );
    time( &endPointM );
    double cpuTimeSeconds = difftime( endPointM, startPointM );
    if ( myTrialRank == 0 ) {
      outy << "rank " << myRank << ":: M=" << M << " multiplier=" << multiplier << " cyclesFound = "
           << cyclesFound <<  " cpuTime =" <<  cpuTimeSeconds << endl;
    }
  }

  mpiError = MPI_Barrier( MPI_COMM_WORLD );
  time( &endPoint );
  double cpuTimeSeconds = difftime( endPoint, startPoint );
  outy.flush();
  cout.flush();
  if ( myRank == 0 ) {
    cout << "Total CPU time = "  << std::fixed << std::setprecision(1) << cpuTimeSeconds << endl;
    outy << "Total CPU time = "  << std::fixed << std::setprecision(1) << cpuTimeSeconds << endl;
    outy.flush();
    cout.flush();
  }

  // Append all non-0 ranks' output to outy file:
  if ( myRank == 0 ) {
    outy.flush();
    ((ofstream*)(outyPtr))->close();
  }
  mpiError = MPI_Barrier( MPI_COMM_WORLD );
  for ( long long rank = 1; rank < nRanks; rank++ ) {
    if ( rank == myRank ) {
      ofstream fileOuty( outputFileName, std::ios::out | std::ios::app );
      fileOuty << ((ostringstream*)(outyPtr))->str();
      fileOuty.flush();
      fileOuty.close();
    }
    mpiError = MPI_Barrier( MPI_COMM_WORLD );
  }

  if ( mpiError != 0 ) {
    cout << "rank " << myRank << ":: mpiError = " << mpiError << endl;
    cout.flush();
  }

  // delete [] configPtr;
  // delete [] outyPtr;
  mpiError = MPI_Barrier( MPI_COMM_WORLD );
  return;
}

bool hcp( HCPGraph& g_, HCPParameters& config_, ostream& outy_, long long trial_, MPI_Comm& trialComm_ )
{
  int mpiError = 0, nRanks = 0, myRank = 0;
  mpiError = MPI_Comm_size( MPI_COMM_WORLD, &nRanks );
  mpiError = MPI_Comm_rank( MPI_COMM_WORLD, &myRank );

  int nRanksPerTrial = 0, myPermRank = 0;
  mpiError = MPI_Comm_size( trialComm_, &nRanksPerTrial );
  mpiError = MPI_Comm_rank( trialComm_, &myPermRank );

  cout.flush();

  time_t startPoint;
  time( &startPoint );

  long long n = g_.nNodes();

  // Read input parameters
  long long naSteps = config_.naSteps; // Maximim number of annealing steps allowed
  double tFactor = config_.tFactor; // Annealing factor (multiply T by this each annealing step)
  double t0 = config_.t0; // Initial T value
  double maxPermMultiplier = config_.maxPermMultiplier;
  long long maxPermutations = (long long)( maxPermMultiplier * pow( n, 2 ) );
  double kMultiplier = config_.kMultiplier;
  double k = kMultiplier * maxPermutations; // Boltzmann-like constant in exponential
  bool retry = config_.retry;
  bool randomEdges = config_.randomEdges;
  bool insertRandomCycle = config_.insertRandomCycle;

  Tour_t tour( n );
  for ( long long nn = 0; nn < n; nn++ ) { tour[ nn ] = nn; }

  //wasButUnused20141011long long seedForShuffle = trial_ * 942;
  //trouble20141011 std::mt19937 shuffleGen( trial_ );
  long long seedForShuffle = trial_ * 99942 + 77742;
  //zippyMake64bit std::mt19937 shuffleGen( seedForShuffle );
  std::mt19937_64 shuffleGen( seedForShuffle );

  shuffle( &tour[ 0 ], &tour[ n ], shuffleGen );
  if ( myRank == 0 ) {
    cout << "debug rank " << myRank << ":: Trial " << trial_ << " initial tour=" << tour << endl;
  }

  // Try eliminating modulo and rand() in moves functions
//  minstd_rand pickNodesGen( trial_ + 3377 * myPermRank ); // ???RIGHT?
  minstd_rand pickNodesGen( trial_ + 3377 * myPermRank + myPermRank ); // try this
  uniform_int_distribution<long long int> pickNodes( 0, ( n - 1 ) );

  double currentCost = costFunction( g_, tour, false );
  if ( currentCost == 0.0 ) {
    if ( myPermRank == 0 ) {
      cout << "rank " << myRank << ":: currentCost == 0.0 upon entry to hcp()" << endl;
      windUp( true, tour, 0, t0, 0, outy_, trial_, g_ ); //zippy:added g_ for debugging
    }
    return true;
  }
  Tour_t newTour( n );
  newTour = tour;
  Tour_t copyOfTour( n );
  copyOfTour = tour;

  // C++11:
  double lowerBound = 0.0, upperBound = 1.0;
  uniform_real_distribution<double> unif( lowerBound, upperBound );
  //zippytrybetter default_random_engine re
  std::random_device rand_dev; // try this for seeding
  std::mt19937_64 re( rand_dev() );

  // Manage multiple ranks in this trial:
  long long myNPerms = maxPermutations / nRanksPerTrial;
  long long myFirstPerm = myPermRank * myNPerms;
  long long myLastPerm = myFirstPerm + myNPerms;
  // for ( long long serializedRank = 0; serializedRank < nRanks; serializedRank++ ) {
  //   if ( myRank == serializedRank ) {
  //     cout << "debug rank " << myRank << " myPermRank=" << myPermRank << ":: myFirstPerm=" << myFirstPerm
  //          << " myLastPerm=" << myLastPerm << "myFirstPerm+myLastPerm=" << (myFirstPerm+myLastPerm) << endl;
  //   }
  // }
 
  // To re-try failures:
  long long seed2, nAttempts = 1;
  long long seedForMetropolis = trial_; // equivalent to original
  if ( retry ) { nAttempts = 2; }
  for ( long long attempt = 0; attempt < nAttempts; attempt++ ) {
    tour = copyOfTour;
    newTour = tour;
    if ( attempt == 0 ) {
      re.seed( seedForMetropolis );
    } else {
      currentCost = costFunction( g_, tour, false );
      seed2 = trial_ * 637 + 2929;
      re.seed( seed2 );
      if ( myPermRank == 0 ) {
        cout << "rank " << myRank << ":: Re-attempting Trial " << trial_ << " with new seed=" << seed2 << endl;
      }
    }

    double t = t0;
    long long zeroDeltas = 0;
    long long acceptedBads = 0;

    mpich_struct_mpi_double_int inS[1], outS[1];
    inS[0].i = myPermRank;

    // declare out here for 2-way parallelism:
    double cost;
    long long permutation;
    long long moveType;
    double deltaCost = 0.0;
    //zippyIsThisABug? mpiError = MPI_Barrier( MPI_COMM_WORLD ); //zippydebug
    mpiError = MPI_Barrier( trialComm_ ); //zippydebugDifferent
    for ( long long step = 1; step <= naSteps; step++ ) {
      if ( myRank == 0 ) {
        cout << "debug rank " << myRank << ":: Trial " << trial_<< std::fixed << std::setprecision(1)
             << " currentCost=" << currentCost
             << " Attempt " << attempt << " Annealing step: " << step << endl;
        cout.unsetf(std::ios_base::floatfield);
      }
      if ( attempt == 1 ) {
        if ( myPermRank == 0 ) {
          cout << "re-attempt, debug rank " << myRank << ":: Trial " << trial_<< std::fixed << std::setprecision(1)
               << " currentCost=" << currentCost
               << " Attempt " << attempt << " Annealing step: " << step << endl;
          cout.unsetf(std::ios_base::floatfield);
        }
      }
      bool foundHC = false;
      for ( permutation = myFirstPerm; ( ( permutation < myLastPerm ) && !foundHC ); permutation++ ) {
        newTour = tour;
        moveType = move( tour, newTour, permutation, n, pickNodes, pickNodesGen );
        cost = costFunction( g_, newTour, false );
        if ( cost == 0.0 ) { 
          cout << "DEBUG: found cost=0 myRank=" << myRank << " trial_=" << trial_ << endl;
          foundHC = true;
          currentCost = cost;
          tour = newTour;
        } else {
          deltaCost = cost - currentCost;
          if ( deltaCost <= 0.0 ) { //zippydebug tried < 0.0
            if ( deltaCost == 0.0 ) {
              zeroDeltas += 1;
            }
            tour = newTour;
            currentCost = cost;
          } else {
            if ( cost > 0.0 ) { //addedFor2WayPar
              double diceRoll = unif( re );
              double threshold = exp( -fabs( deltaCost )/( kMultiplier * t ) );
              if ( diceRoll < threshold ) {
                acceptedBads += 1;
                tour = newTour;
                currentCost = cost;
              }
            }
          }
        }
      }
      mpiError = MPI_Barrier( trialComm_ );

      // Find the minimum-cost of all the ranks' tours and reset all to that:
      inS[0].d = currentCost;
      double costHere = currentCost;
      // zippy: try shaking it up sometimes:
      // double diceRoll2 = unif( re );
      // double threshold = 1.0e-2;
      // double diceRoll3 = 0.0;
      // mpiError = MPI_Allreduce( &diceRoll2, &diceRoll3, 1, MPI_DOUBLE, MPI_MIN, trialComm_ );
      // if ( (myRank==0) || (myRank==1) || (myRank==16) ) {
      //   cout << "debug: shake it up: rank " << myRank << " diceRoll3=" << diceRoll3 << endl;
      // }
      // if ( diceRoll3 < threshold ) {
      //   acceptedBads += 1;
      //   mpiError = MPI_Allreduce( inS, outS, 1, MPI_DOUBLE_INT, MPI_MAXLOC, trialComm_ );
      //   double cc = outS[0].d;
      //   cout << "debug: bad max accepted (shake it up): rank " << myRank << " currentCost=" << currentCost << " cc=" << cc
      //        << " diceRoll2=" << diceRoll2 << " diceRoll3=" << diceRoll3 << " myPermRank=" << myPermRank << " trial_="
      //        << trial_ << endl;
      // } else {
        mpiError = MPI_Allreduce( inS, outS, 1, MPI_DOUBLE_INT, MPI_MINLOC, trialComm_ );
      // }
      currentCost = outS[0].d;
      cost = currentCost;
      newTour = tour; //zippy: maybe this will fix the bug...
      mpiError = MPI_Bcast( &newTour[0], n, MPI_LONG_LONG, outS[0].i, trialComm_ );
      tour = newTour;
      if ( cost == 0.0 ) {
        long long totalAcceptedBads = 0;
        mpiError = MPI_Reduce( &acceptedBads, &totalAcceptedBads, 1, MPI_LONG_LONG, MPI_SUM, 0,
                               trialComm_ );
        if ( myPermRank == 0 ) {
          cout << "rank " << myRank << ":: Found HC; acceptedBads = " << totalAcceptedBads << endl;
          windUp( true, newTour, step, t, permutation, outy_, trial_, g_ ); //zippy:added g_ for debugging
        }
        if ( attempt == 1 ) {
          if ( myPermRank == 0 ) {
            cout << "rank " << myRank << ":: on re-attempt, about to return from hcp" << endl;
          }
        }
        return true;
      }

      t *= tFactor;  // Anneal, end of current SA step
    }

    if ( myPermRank == 0 ) {
      cout << "rank = " << myRank << ":: zeroDeltas = " << zeroDeltas << " acceptedBads = " << acceptedBads
           << endl;
    }
  } // to re-try failures.

  if ( myPermRank == 0 ) {
      cout << "rank " << myRank << ":: trial " << trial_ <<  " Found no HC <<  currentCost =" << std::fixed << std::setprecision(1) << currentCost << endl;
      outy_ << "rank " << myRank << ":: trial " << trial_ <<  " Found no HC <<  currentCost =" << std::fixed << std::setprecision(1) << currentCost << endl;
  }
  return false;
}

int main( int argc, char* argv[] )
{
  int mpiError = 0, nRanks = 0, myRank = 0;
  mpiError = MPI_Init( &argc, &argv );
  MPI_Comm_size( MPI_COMM_WORLD, &nRanks );
  MPI_Comm_rank( MPI_COMM_WORLD, &myRank );

  long long N = atoi( argv[ 1 ] );
  string inputFilename = argv[ 2 ];
  if ( myRank == 0 ) {
    cout << "rank " << myRank << ":: N=" << N << " inputFilename=" << inputFilename << endl;
  }
  sweep( N, inputFilename );

  mpiError = MPI_Barrier( MPI_COMM_WORLD );
  MPI_Finalize();

  return 0;
}

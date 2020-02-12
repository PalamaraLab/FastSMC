#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <bitset>
#include <string.h>
#include <stdlib.h>
#include <inttypes.h>
#include <time.h>
#include <sys/time.h>
#include <cstdint>
#include <stdint.h>
#include <unistd.h>
#include <cmath>
#include <math.h>

// Includes for seeding
#include <boost/unordered_map.hpp>

#include "Data.hpp"
#include "DecodingParams.hpp"
#include "DecodingQuantities.hpp"
#include "FileUtils.hpp"
#include "StringUtils.hpp"
#include "Timer.hpp"
#include "HMM.cpp"

using namespace std;

vector< int > hist_ctr;
vector< float > hist_win;
HMM hmm;

// Set the word/hash sizes here at compile time, must be able to cast from a ulong
float PAR_MIN_MATCH = 1;
bool PAR_DIAGNOSTICS = 0;
int PAR_GAP = 1;
int MAX_seeds = 0;
int GLOBAL_READ_WORDS = 0;
int GLOBAL_CURRENT_WORD = 0;
int GLOBAL_SKIPPED_WORDS = 0;
int GEN_THRESHOLD = 100;
bool PAR_BIN_OUT = 0;
bool PAR_HAPLOID = 1;

bool is_j_above_diag;
unsigned int windowSize;
unsigned int w_i;
unsigned int w_j;
int jobID, jobs;

vector< double > *all_markers;
unsigned long int num_ind;
unsigned long int num_ind_tot;

const int CONST_READ_AHEAD = 10;
const int WORD_SIZE = 64;
typedef uint64_t hash_size;

struct Marker {
  string id;
  unsigned long int pos;
  double cm;

  string print() {
    stringstream ss;
    ss << id << '\t' << pos << '\t' << cm << endl;
    return ss.str();
  }
};

class Individuals {
  unsigned int idnum;
public:
  bitset<WORD_SIZE> hap[ CONST_READ_AHEAD ];

  void clear( int w ) {
    hap[ w % CONST_READ_AHEAD ].reset();
  }
  void setMarker( int w , int bit ) {
    hap[ w % CONST_READ_AHEAD ].set( bit );
  }
  hash_size getWordHash( int w ) {
    return hap[ w % CONST_READ_AHEAD ].to_ulong();
  }
  string getWordString( int w ) {
    return hap[ w % CONST_READ_AHEAD ].to_string();
  }
  unsigned int getNum() { return idnum; }
  Individuals(unsigned int);
};
Individuals::Individuals( unsigned int iid ) {
  idnum = iid;
  for ( int w = 0; w < CONST_READ_AHEAD ; w++ ) clear( w );
}

vector< Individuals > all_ind;

// Convenience function to compute genetic distance between two words (start of w1 and end of w2)
double cmBetween( int w1 , int w2 ) {
  int end = WORD_SIZE * w2 + WORD_SIZE - 1;
  if ( end >= all_markers->size() ) end = all_markers->size() - 1;
  return 100*( (*all_markers)[end] - (*all_markers)[WORD_SIZE * w1]);
}

bool isHapInJob( unsigned int i ) {
  return ( (i >= (w_i-1)*windowSize && i < w_i*windowSize) || (i >= (w_j-1)*windowSize && i < w_j*windowSize) || ( jobs == jobID && i >= (w_j-1)*windowSize ) );
}

bool isSampleInJob( unsigned int i ) {
  return ( ( i >= (uint)((w_i-1)*windowSize)/2 && i < (uint)(w_i*windowSize)/2 ) || ( i >= (uint)((w_j-1)*windowSize)/2 && i < (uint)(w_j*windowSize)/2 ) || ( jobs == jobID && i >= (uint)((w_j-1)*windowSize)/2 ) );
}

struct Match {
  int interval[2] = {0,0};
  unsigned int gaps = 0;
  // pair : identifiers for the corresponding Individualss in all_ind
  void print( pair<unsigned int,unsigned int> p ) {
    double mlen = cmBetween(interval[0],interval[1]);
    if ( mlen >= PAR_MIN_MATCH ) {
      hmm.decodeFromGERMLINE(p.first, p.second, interval[0] * WORD_SIZE, interval[1] * WORD_SIZE + WORD_SIZE - 1);
    }
  }

  void extend( int w ) {
    if ( interval[1] < w ) interval[1] = w;
  }
  void addGap() {
    gaps++;
  }
  Match(int);
  Match(void);
};
Match::Match( int i ) { interval[0] = interval[1] = i; }
Match::Match() { interval[0] = interval[1] = 0; }

/* Object for storing extension between pairs of Individualss */
class ExtendHash {
  boost::unordered_map<unsigned long int, Match > extend_hash;
  unsigned long int num = 0;
  // Empty Match to insert into hash
  Match m;
  // Iterator for testing insertion
  std::pair<boost::unordered::iterator_detail::iterator<boost::unordered::detail::ptr_node<std::pair<const unsigned long int, Match > > >, bool> extend_ret;
public:
  ExtendHash(unsigned long int);
  // Compute pair of Individualss from location indicator
  pair<unsigned int,unsigned int> locationToPair( unsigned long int loc ) {
    pair<unsigned int,unsigned int> p;
    // round everyone down to the nearest haplotype
    if ( !PAR_HAPLOID ) {
      p.second = 2 * (loc % num);
      p.first = 2 * ((loc - p.second/2) / num);
    } else {
      p.second = loc % num;
      p.first = (loc - p.second) / num;
    }
    return p;
  }
  // Compute location from pair of Individualss
  unsigned long int pairToLocation( unsigned int i , unsigned int j ) {

    if ( !PAR_HAPLOID ) {
      // round everyone down to the nearest haplotype
      i = (i - (i % 2)) / 2;
      j = (j - (j % 2)) / 2;
    }
    unsigned long int loc = (i > j) ? j * num + i : i * num + j;
    return loc;
  }
  // Extend or add a given pair in the current hash
  // unsigned int i,j : identifiers for the two Individualss
  // int w : current word # to extend or add
  void extendPair( unsigned int i , unsigned int j , int w ) {
    m.interval[0] = GLOBAL_CURRENT_WORD;
    // Find/extend this location in the hash
    extend_ret = extend_hash.insert( pair< unsigned long int , Match >( pairToLocation(i,j) , m) );
    (extend_ret.first->second).extend(w);
  }

  // Remove all pairs that were not extended beyond w
  // int w : word # to remove prior to
  void clearPairsPriorTo(int w) {
    for ( auto it = extend_hash.begin() ; it != extend_hash.end() ; ) {
      if ( it->second.interval[1] < w ) {
        it->second.print( locationToPair(it->first) );
        it = extend_hash.erase( it );
      } else {
        if ( it->second.interval[1] < GLOBAL_CURRENT_WORD ) it->second.addGap();
        it++;
      }
    }
  }

  // Remove all pairs that were not extended beyond w
  // int w : word # to remove prior to
  void extendAllPairsTo(int w) {
    for ( auto it = extend_hash.begin() ; it != extend_hash.end() ; it++ ) it->second.interval[1] = w;
  }

  // Remove all pairs
  // int w : word # to remove prior to
  void clearAllPairs() {
    for ( auto it = extend_hash.begin() ; it != extend_hash.end() ; ) {
      it->second.print( locationToPair(it->first) );
      it = extend_hash.erase( it );
    }
  }
  int size() {
    return extend_hash.size();
  }
};
ExtendHash::ExtendHash(unsigned long int n) { num = n; }

/* Object for storing initial word seeds */
class SeedHash {
  boost::unordered_map<hash_size, vector<unsigned int> > seed_hash;
  // Empty vector to insert into the seed hash
  vector<unsigned int> vec;
  // Iterator for testing insertion of elements
  // std::pair<boost::unordered::iterator_detail::iterator<boost::unordered::detail::ptr_node<std::pair<const unsigned int, vector<unsigned int> > > >, bool> seed_ret;
public:
  void insertIndividuals( unsigned int i , hash_size word ) {
    auto seed_ret = seed_hash.insert( pair<hash_size, vector<unsigned int>> ( word , vec ) );
    (seed_ret.first->second).push_back( i );
  }
  void clear() {
    seed_hash.clear();
  }
  int size() {
    return seed_hash.size();
  }

  // Generate a new hash for this vector of Individualss
  unsigned long subHash( ExtendHash * e , vector<unsigned int> vec , int w ) {
    SeedHash cur_sh;
    // seed the next word from this subset of Individualss
    for ( int i = 0 ; i < vec.size() ; i++ ) cur_sh.insertIndividuals( vec[i] , all_ind[ vec[i] ].getWordHash( w ) );
    // recursion:
    // cerr << "\tsubhash seeds: " << w << " " << cur_sh.size() << endl;
    return cur_sh.extendAllPairs( e , w );
  }
  // Extend/save all pairs in the current hash
  // ExtendHash * e : Pointer to ExtendHash which will be called for each pair
  // returns : number of pairs evaluated
  unsigned long extendAllPairs( ExtendHash * e , int w ) {
    unsigned long tot_pairs = 0;
    for ( auto it = seed_hash.begin() ; it != seed_hash.end() ; ++it ) {

      // *** As long as the # of pairs is high, generate a sub-hash for the next word
      // *** Only store pairs of Individualss that have collision in a small hash
      // *** Extend only to the haplotypes that seeded here
      if ( MAX_seeds != 0 && it->second.size() > MAX_seeds && w + 1 < GLOBAL_READ_WORDS ) {
        // recursively generate a sub-hash
        // IMPORTANT: if we run out of buffered words then this seed does not get analyzed
        if ( w + 1 < GLOBAL_READ_WORDS ) tot_pairs += subHash( e , it->second , w + 1 );
        else GLOBAL_SKIPPED_WORDS++;
      } else {
        //tot_pairs += it->second.size() * (it->second.size() - 1) / 2;
        for ( int i = 0 ; i < it->second.size() ; i++ ) {
          for ( int ii = i+1 ; ii < it->second.size() ; ii++ ) {

            unsigned int ind_i = std::max(it->second[i], it->second[ii]);
            unsigned int ind_j = std::min(it->second[i], it->second[ii]);

            // for the last job only
            if ( jobID == jobs ) {
              if ( all_ind[ind_i].getNum() >= (w_i-1)*windowSize && all_ind[ind_j].getNum() >= (w_j-1)*windowSize ) {
                if ( all_ind[ind_j].getNum() <  (w_j-1)*windowSize + (all_ind[ind_i].getNum() - (w_i-1)*windowSize) ) {
                  e->extendPair( ind_j , ind_i , w );
                  tot_pairs++;
                }
              }
            }

              // for all other jobs
            else if ( (all_ind[ind_i].getNum() >= (w_i-1)*windowSize && all_ind[ind_i].getNum() < w_i*windowSize) && (all_ind[ind_j].getNum() >= (w_j-1)*windowSize && all_ind[ind_j].getNum() < w_j*windowSize) ) {
              if ( is_j_above_diag && all_ind[ind_j].getNum() <  (w_j-1)*windowSize + (all_ind[ind_i].getNum() - (w_i-1)*windowSize) ) {
                e->extendPair( ind_j , ind_i , w );
                tot_pairs++;
              } else if ( !is_j_above_diag && all_ind[ind_j].getNum() >= (w_j-1)*windowSize + (all_ind[ind_i].getNum() - (w_i-1)*windowSize) ) {
                e->extendPair( ind_j , ind_i , w );
                tot_pairs++;
              }
            }

          }
        }
      }
    }
    return tot_pairs;
  }
};

double get_cpu_time(){
  return (double)clock() / CLOCKS_PER_SEC;
}

inline bool fileExists(const std::string& name) {
  ifstream f(name.c_str());
  return f.good();
}

int main (int argc, char* argv[])
{
  srand(1234);

  DecodingParams params;

  double TIME_start = get_cpu_time();
  double TIME_prev = TIME_start;
  float PAR_MIN_MAF = 0;
  float PAR_skip = 0;

  string line, discard;
  bool opt_error = 0;
  int c;

  // parse input arguments
  if (!params.processCommandLineArgs(argc, argv)) {
    cerr << "Error processing command line; exiting." << endl;
    exit(1);
  }

  PAR_MIN_MATCH = params.min_m;
  PAR_MIN_MAF = params.min_maf;
  PAR_GAP = params.gap;
  PAR_skip = params.skip;
  MAX_seeds = params.max_seeds;

  // used for benchmarking
  Timer timer;
  // read decoding quantities from file
  std::string str(params.decodingQuantFile.c_str());
  DecodingQuantities decodingQuantities(params.decodingQuantFile.c_str());
  printf("\n*** Read decoding quantities in %.3f seconds. ***\n\n", timer.update_time());

  cout << "---------------------------" << endl;
  cout << "      	 READING DATA     	" << endl;
  cout << "---------------------------" << endl;
  cout << "CSFS samples: " << decodingQuantities.CSFSSamples << endl;
  unsigned int sequenceLength = Data::countHapLines(params.inFileRoot.c_str());

  //ifstream file_haps(params.inFileRoot + ".hap");
  FileUtils::AutoGzIfstream file_haps;
  if (fileExists(params.inFileRoot + ".hap.gz")) {
    file_haps.openOrExit(params.inFileRoot + ".hap.gz");
  } else if (fileExists(params.inFileRoot + ".hap")) {
    file_haps.openOrExit(params.inFileRoot + ".hap");
  } else if (fileExists(params.inFileRoot + ".haps.gz")) {
    file_haps.openOrExit(params.inFileRoot + ".haps.gz");
  } else if (fileExists(params.inFileRoot + ".haps")) {
    file_haps.openOrExit(params.inFileRoot + ".haps");
  } else {
    cerr << "ERROR. Could not find hap file in " + params.inFileRoot + ".hap.gz, " + params.inFileRoot + ".hap, " + params.inFileRoot + ".haps.gz, or " + params.inFileRoot + ".haps" << endl;
    exit(1);
  }

  //ifstream file_samp(params.inFileRoot + ".samples");
  FileUtils::AutoGzIfstream file_samp;
  if (fileExists(params.inFileRoot + ".samples")) {
    file_samp.openOrExit(params.inFileRoot + ".samples");
  } else if (fileExists(params.inFileRoot + ".sample")) {
    file_samp.openOrExit(params.inFileRoot + ".sample");
  } else {
    cerr << "ERROR. Could not find sample file in " + params.inFileRoot + ".sample or " + params.inFileRoot + ".samples" << endl;
    exit(1);
  }

  ifstream file_genm(params.map);
  if(!file_genm ) { cerr << params.map << " could not be opened" << endl; return -1; }

  string map_field[3];
  stringstream ss;

  // *** read genetic map
  vector< pair<unsigned long int,double> > genetic_map;
  unsigned int cur_g = 0;
  while(getline(file_genm,line)) {
    ss.clear(); ss.str(line);
    ss >> map_field[0] >> map_field[1] >> map_field[2];
    if ( map_field[0] == "position" || map_field[0] == "" ) continue;
    genetic_map.push_back( pair<unsigned long int,double>( stol(map_field[0]) , stod(map_field[2]) ) );
    if ( cur_g > 0 && (genetic_map[ cur_g ].first < genetic_map[ cur_g - 1 ].first || genetic_map[ cur_g ].second < genetic_map[ cur_g - 1 ].second) ) {
      cerr << "ERROR: genetic map not in sorted order at line\n" << line << endl;
      return -1;
    }
    cur_g++;
  }
  file_genm.close();
  if ( genetic_map.size() < 2 ) { cerr << "ERROR: genetic map must have at least two valid entries" << endl; return -1; }

  //cerr << "*** runtime : " << get_cpu_time() - TIME_start << "\t";
  cerr << genetic_map.size() << " genetic map entries read" << endl;
  printf("\n*** Read genetic map in %.3f seconds. ***\n\n", timer.update_time());

  Data data(params.inFileRoot.c_str(), sequenceLength, decodingQuantities.CSFSSamples, params.foldData, params.usingCSFS, params.jobInd, params.jobs, genetic_map);
  printf("\n*** Read haps in %.3f seconds. ***\n\n", timer.update_time());

  //freeing memory
  genetic_map.clear();
  vector< pair<unsigned long int,double> >().swap( genetic_map );

  cout << "---------------------------" << endl;
  cout << "   INFERRING IBD SEGMENTS  " << endl;
  cout << "---------------------------" << endl;

  hmm.buildHMM(data, decodingQuantities, params);
  hmm.decodeAll(params.jobs, params.jobInd);

  if (!params.GERMLINE) {
    exit(1);
  }

  jobs = params.jobs;
  jobID = params.jobInd;
  windowSize = data.windowSize;
  w_i = data.w_i;
  w_j = data.w_j;
  is_j_above_diag = data.is_j_above_diag;

  // *** read Individuals information
  unsigned int linectr = 0;
  while(getline(file_samp,line)) {
    vector <string> splitStr;
    istringstream iss(line); string buf; while (iss >> buf) splitStr.push_back(buf);

    // Skip first two lines (header) if present
    if ((splitStr[0] == "ID_1" && splitStr[1] == "ID_2" && splitStr[2] == "missing")
        || (splitStr[0] == "0" && splitStr[1] == "0" && splitStr[2] == "0")) {
      continue;
    }

    ss.clear(); ss.str( line );
    ss >> map_field[0] >> map_field[1];
    if ( isSampleInJob( linectr ) ) {
      if ( PAR_HAPLOID ) {
        all_ind.push_back( Individuals(2*linectr) );
        all_ind.push_back( Individuals(2*linectr+1) );
      } else {
        all_ind.push_back( Individuals(2*linectr) );
        all_ind.push_back( Individuals(2*linectr) );
      }
    }
    linectr++;
  }
  file_samp.close();
  num_ind_tot = data.sampleSize*2;
  num_ind = all_ind.size();

  cerr << num_ind / 2 << " sample identifiers read" << endl;

  Marker cur_marker;
  // track position through genetic map
  cur_g = 0;
  unsigned int snp_ctr;
  char al[2] , inp;
  string cur_al;

  // Storage for seeds
  SeedHash seeds;

  hash_size word[2];

  // Storage for extensions
  ExtendHash extend( PAR_HAPLOID ? num_ind : num_ind / 2 );

  // Hash individual words
  GLOBAL_READ_WORDS = 0;
  GLOBAL_CURRENT_WORD = 0;

  string marker_id;
  unsigned long int marker_pos;
  all_markers = &data.geneticPositions;
  while( 1 ) {
    snp_ctr = 0;
    while( getline(file_haps,line) )
    {
      // read the meta data
      ss.clear(); ss.str( line );
      ss >> map_field[0] >> marker_id >> marker_pos >> al[0] >> al[1];
      if(map_field[0] == "") continue;

      // restrict on MAF
      if ( PAR_MIN_MAF > 0 ) {
        int maf_ctr = 0;
        for( int i=0; i< num_ind_tot ; i++) {
          ss >> inp;
          if ( inp == '1' ) maf_ctr++;
        }
        float maf = (float) maf_ctr / num_ind_tot;
        if ( maf < PAR_MIN_MAF || maf > 1 - PAR_MIN_MAF ) continue;

        // re-load the data
        ss.clear(); ss.str( line );
        ss >> map_field[0] >> map_field[0] >> map_field[0] >> al[0] >> al[1];
      }

      // read haplotype
      unsigned int hap_ctr = 0;
      for( unsigned int i=0; i<num_ind_tot ; i++ )
      {
        ss >> inp;
        if ( isSampleInJob(i/2) ) {
          if ( inp == '1' ) all_ind[hap_ctr].setMarker( GLOBAL_READ_WORDS , snp_ctr );
          hap_ctr++;
        }
      }
      snp_ctr++;

      if ( snp_ctr % WORD_SIZE == 0 ) {
        if ( ++GLOBAL_READ_WORDS >= CONST_READ_AHEAD ) break;
        else cerr << "*** loading word buffer " << GLOBAL_READ_WORDS << " / " << CONST_READ_AHEAD << endl;
        snp_ctr = 0 ;
      }
    }

    // end if read all data
    if( GLOBAL_CURRENT_WORD >= GLOBAL_READ_WORDS ) {
      break;
    }

    for ( unsigned int i = 0 ; i < num_ind ; i++ ) {
      seeds.insertIndividuals( i , all_ind[i].getWordHash( GLOBAL_CURRENT_WORD ) );
    }

    GLOBAL_SKIPPED_WORDS = 0;
    int cur_seeds = seeds.size();
    unsigned long cur_pairs = 0;

    // skip low-complexity words
    if ( (float) cur_seeds / num_ind > PAR_skip ) {
      cur_pairs = seeds.extendAllPairs( &extend , GLOBAL_CURRENT_WORD );
      extend.clearPairsPriorTo( GLOBAL_CURRENT_WORD - PAR_GAP );
    } else {
      cerr << "low complexity word - " << cur_seeds << " - skipping" << endl;
      extend.extendAllPairsTo( GLOBAL_CURRENT_WORD );
    }

    seeds.clear();

    for( unsigned int i=0; i< num_ind ; i++) {
      all_ind[i].clear( GLOBAL_CURRENT_WORD );
    }
    GLOBAL_CURRENT_WORD++;
  }

  extend.clearAllPairs();
  file_haps.close();

  hmm.finishFromGERMLINE();
  cerr << "processed " << GLOBAL_CURRENT_WORD * WORD_SIZE << " / " << all_markers->size() << " SNPs" << endl;

  printf("\n*** Inference done in %.3f seconds. ***\n\n", timer.update_time());

  return 0;
}


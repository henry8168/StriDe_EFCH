#ifndef TYPEDEF_INTERVALMATCHMAP_H
#define TYPEDEF_INTERVALMATCHMAP_H

//#include "HashMap.h"
#include "tbb/concurrent_hash_map.h"
#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"
#include <string>
#include <vector>
#include "BWTAlgorithms.h"

//typedef google::sparse_hash_map<int64_t, std::vector<int> > IntervalMatchMap;

/*
// Structure that defines hashing and comparison operations for user's type.
struct MyHashCompare {
    static size_t hash( const std::string& x ) {
        size_t h = 0;
        for( const char* s = x.c_str(); *s; ++s )
            h = (h*17)^*s;
        return h;
    }
    //! True if strings are equal
    static bool equal( const std::string& x, const std::string& y ) {
        return x==y;
    }
};
*/
// A concurrent hash table that maps strings to ints.
typedef tbb::concurrent_hash_map<int64_t , std::vector<std::pair<int, unsigned short> > > IntervalMatchMap;  //int64_t代表interval的lower bound，後方的vector代表存入的long read names和該次的steps

class intervalPackage
{
    public:
    BWTInterval FBintervalRBWTx, FEintervalRBWTx;
    BWTInterval RBintervalBWTx, REintervalBWTx;
};

class itemsForOutput
{
    public:
        int lrid;
        std::string lrname;
        DNAString mergedread;
        intervalPackage FRBE[2];
        bool pass;
    
        itemsForOutput():pass(true)
        {
            
        }
};

#endif
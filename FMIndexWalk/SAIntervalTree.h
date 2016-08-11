//----------------------------------------------
// Copyright 2014 National Chung Cheng University
// Written by Yao-Ting Huang
// Released under the GPL
//-----------------------------------------------
//
//
// Re-written from Jared Simpson's StringThreaderNode and StringThreader class
// The search tree represents a traversal through implicit FM-index graph
//
#ifndef SAINTERVALTREE_H
#define SAINTERVALTREE_H

#include <list>
#include "BWT.h"
//#include "HashMap.h"
#include "typedef_IntervalMatchMap.h"
#include "BWTAlgorithms.h"
#include <omp.h>

// Typedefs
class SAIntervalNode;
typedef std::list<SAIntervalNode*> STNodePtrList;

// Object to hold the result of the threading process
struct SAIntervalNodeResult
{
    std::string thread;
	size_t SAICoverage;
    std::vector<int64_t> pathIntervalBWT[2];  //mod2
    std::vector<int64_t> pathIntervalRBWT[2];
    std::vector<unsigned short> order;
};
typedef std::vector<SAIntervalNodeResult> SAIntervalNodeResultVector;

// A node in the threading tree
class SAIntervalNode
{
    public:
        
        //
        // Functions
        //
        SAIntervalNode(const std::string* pQuery, SAIntervalNode* parent);
        ~SAIntervalNode();

        // Add a child node to this node with the given label
        // Returns a pointer to the created node
        SAIntervalNode* createChild(const std::string& label);

        // Extend the label of this node by l
        void extend(const std::string& ext);

        // Return a suffix of length l of the string represented by this node
        std::string getSuffix(size_t l) const;

        // Return the complete sequence of the string represented by the branch
        std::string getFullString() const;


        // Initialize or update the alignment data
        //void computeInitialAlignment(const std::string& initialLabel, int queryAlignmentEnd, int bandwidth);
        void computeInitial(const std::string& initialLabel);

        // Recursive function to print all the strings represented
        // by this node and all its children.
        void printAllStrings(const std::string& parent) const;


        size_t getKmerCount(){return m_kmerCount;};
        void addKmerCount(size_t kmercount){ m_kmerCount+=kmercount; };
        // size_t getKmerSize(){return m_kmersize;};
        // void setKmerSize(size_t kmersize){ m_kmersize=kmersize; };

        BWTInterval fwdInterval;
        BWTInterval rvcInterval;
        
        //Interval
        std::vector<int64_t> m_pathIntervalBWT[2]; //mod3
        std::vector<int64_t> m_pathIntervalRBWT[2];

    //private:

        //
        // Data
        //

        // The extension string from the parent
        std::string m_label;
        size_t m_kmerCount;
        // size_t m_kmersize;


        // The query string being threaded through the graph
        const std::string* m_pQuery;
        
        std::vector<unsigned short> m_order;//mod4

        // The parent node, can be NULL
        SAIntervalNode* m_pParent;
        STNodePtrList m_children;

};

class SAIntervalTree
{
    public:
        SAIntervalTree(const std::string* pQuery,
                       SAIntervalNodeResult* resultthis,
                       intervalPackage &intervalP,
                       IntervalMatchMap* intervalMatchMapBWT,
                       IntervalMatchMap* intervalMatchMapRBWT,
                       bool &uniqueCase,
                       size_t coverage,
                       size_t compressionLevel,
                       size_t minOverlap,
					   size_t maxOverlap,
                       size_t MaxLength,
                       size_t MaxLeaves,
                       //size_t meaninsertsize,
                       //size_t stddev,
					   BWTIndexSet indices,
                       std::string secondread,
                       size_t SA_threshold=3,
                       bool KmerMode=false);

        ~SAIntervalTree();

        //return the merged string
        //bool mergeTwoReads(StringVector & mergeReads);
        int mergeTwoReads(std::string &mergedseq);
		size_t getKmerCoverage(){return m_maxKmerCoverage;};
		size_t getMaxUsedLeaves(){return m_maxUsedLeaves;};
		bool isBubbleCollapsed(){return m_isBubbleCollapsed;}

        // Print all the strings represented by the tree
        void printAll();

    private:

        //
        // Functions
        //
        void extendLeaves();
        void attempToExtend(STNodePtrList &newLeaves);
        void refineSAInterval(size_t newKmerSize);
        std::vector<std::pair<std::string, BWTIntervalPair> > getFMIndexExtensions(SAIntervalNode* pNode);

        // Check if the leaves can be extended no further
        bool isTerminated(SAIntervalNodeResultVector& results);
        bool isTwoReadsOverlap(std::string & mergedseq);
        size_t calculateKmerCoverage (const std::string & seq , size_t kmerLength , const BWT* pBWT);
		bool replaceLowFreqKmer (std::string & seq , size_t kmerLength);

        void removeLeavesByRepeatKmer();

        //
        // Data
        //
        const std::string* m_pQuery;
        SAIntervalNodeResult* m_resultthis; ///mod1-1
        IntervalMatchMap *m_intervalMatchMapBWT;
        IntervalMatchMap *m_intervalMatchMapRBWT;
        IntervalMatchMap::accessor immAccessor;
        size_t m_coverage;
        size_t m_compressionLevel;
        size_t m_minOverlap;
		size_t m_maxOverlap;
        size_t m_MaxLength;
        size_t m_MaxLeaves;
        //size_t m_meaninsertsize;
        //size_t m_stddev;
		BWTIndexSet m_indices;
        std::string m_secondread;
        size_t m_min_SA_threshold;
        bool m_kmerMode;
        
        bool fwdReadsExisted[2]; //does beginning kmer or ending kmer exist in hash table? ([0] for begining, [1] for ending) //mod1-2
        bool rvcReadsExisted[2];
        int64_t intervalGap;
        size_t atmosttimess;
        BWTInterval beginningkmerBWT,beginningkmerRBWT,endingkmerBWT,endingkmerRBWT; //mod5
        unsigned short beginningkmerOrder, endingkmerOrder;
        std::vector<std::pair<int, unsigned short> > fwdBeginningIdOrder; //mod1-3
        std::vector<std::pair<int, unsigned short> > rvcBeginningIdOrder;
        std::vector<std::pair<int, unsigned short> > fwdEndingIdOrder;
        std::vector<std::pair<int, unsigned short> > rvcEndingIdOrder;
        int FBsize,FEsize,RBsize,REsize;
        int FBid, FEid, RBid, REid;
        unsigned short FBorder, FEorder, RBorder, REorder;
        unsigned short cur_order;
        unsigned short orderGap;

        SAIntervalNode* m_pRootNode;
        STNodePtrList m_leaves;
        DenseHashMap<std::string, size_t, StringHasher> m_KmerIndexMap;
        size_t m_currentLength;
		size_t m_currentKmerSize;
		size_t m_maxKmerCoverage;
		size_t m_maxUsedLeaves;
		bool m_isBubbleCollapsed;

        BWTInterval m_fwdTerminatedInterval;   //in rBWT 
        BWTInterval m_rvcTerminatedInterval;   //in BWT
};

#endif

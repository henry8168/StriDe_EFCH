///----------------------------------------------
// Copyright 2014 National Chung Cheng University
// Written by Yao-Ting Huang
// Released under the GPL
//-----------------------------------------------
//
// SAIntervalTree - Iteratively construct a
// string representing a walk through an assembly graph
// matching a query sequence.
//
// The assembly graph is abstractly represented as
// an FM-index.
//
#include "SAIntervalTree.h"
#include "BWTAlgorithms.h"
#include <vector>

//
// SAIntervalNode
//
SAIntervalNode::SAIntervalNode(const std::string* pQuery, SAIntervalNode* parent) : 
									   m_kmerCount(0), m_pQuery(pQuery), m_pParent(parent)
{

}

// Destructor, recurisvely delete the children of the node
SAIntervalNode::~SAIntervalNode()
{
    // Delete children
    for(STNodePtrList::iterator iter = m_children.begin(); iter != m_children.end(); ++iter)
        delete *iter;

}

// Return a suffix of length l of the path from the root to this node
std::string SAIntervalNode::getSuffix(size_t l) const
{
    size_t n = m_label.size();
    if(l <= n)
    {
        return m_label.substr(n - l, l);
    }
    else
    {
        assert(m_pParent != NULL);
        return m_pParent->getSuffix(l - n) + m_label;
    }
}

// Return the full string of the path from the root to this node
std::string SAIntervalNode::getFullString() const
{
    if(m_pParent == NULL)
        return m_label;
    else
        return m_pParent->getFullString() + m_label;
}

// Create a new child node with the given label. Returns a pointer to the new node.
SAIntervalNode* SAIntervalNode::createChild(const std::string& label)
{
    SAIntervalNode* pAdded = new SAIntervalNode(m_pQuery, this);
    m_children.push_back(pAdded);

    //assert(!m_alignmentColumns.empty());
    //pAdded->computeExtendedAlignment(label, m_alignmentColumns.back());
    pAdded->extend(label);

    return pAdded;
}

// Extend the label of this node
void SAIntervalNode::extend(const std::string& ext)
{
    assert(!ext.empty());
    //assert(!m_alignmentColumns.empty());
    m_label.append(ext);
}


void SAIntervalNode::computeInitial(const std::string& initialLabel)
{
    m_label = initialLabel;

}


// Print the string(s) represented by this node and its children
void SAIntervalNode::printAllStrings(const std::string& parent) const
{
    if(m_children.empty())
    {
        std::cout << ">\n" << parent + m_label << "\n";
    }
    else
    {
        for(STNodePtrList::const_iterator iter = m_children.begin(); iter != m_children.end(); ++iter)
            (*iter)->printAllStrings(parent + m_label);
    }
}

//
// Class: SAIntervalTree
SAIntervalTree::SAIntervalTree(const std::string* pQuery,
                               SAIntervalNodeResult* resultthis,
                               intervalPackage &intervalP,
                               IntervalMatchMap* intervalMatchMapBWT,
                               IntervalMatchMap* intervalMatchMapRBWT,
                               bool& uniqueCase,
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
                               size_t SA_threshold,
                               bool KmerMode) :
                               m_pQuery(pQuery), m_resultthis(resultthis), m_intervalMatchMapBWT(intervalMatchMapBWT), m_intervalMatchMapRBWT(intervalMatchMapRBWT), m_coverage(coverage),
                               m_compressionLevel(compressionLevel), m_minOverlap(minOverlap), m_maxOverlap(maxOverlap), m_MaxLength(MaxLength),
                               m_MaxLeaves(MaxLeaves), /*m_meaninsertsize(meaninsertsize), m_stddev(stddev),*/ m_indices(indices), 
                               m_secondread(secondread), m_min_SA_threshold(SA_threshold),
                               m_kmerMode(KmerMode), m_maxKmerCoverage(0), m_maxUsedLeaves(0), m_isBubbleCollapsed(false)
                               
{
    // Create the root node containing the seed string
    m_pRootNode = new SAIntervalNode(pQuery, NULL);
    m_pRootNode->computeInitial(*pQuery);   //store initial str of root
    m_leaves.push_back(m_pRootNode);
    m_currentLength=pQuery->length();
	m_currentKmerSize=m_minOverlap;
    
    if(m_compressionLevel)
    {
        fwdReadsExisted[0] = false;
        rvcReadsExisted[0] = false;
        fwdReadsExisted[1] = false;
        rvcReadsExisted[1] = false;
        /*
        fwdBeginningIdOrder.clear();
        rvcBeginningIdOrder.clear();
        fwdEndingIdOrder.clear();
        rvcEndingIdOrder.clear();
        */
        intervalGap = (int)(m_coverage*0.005*(11-m_compressionLevel));
        if(intervalGap<1) intervalGap=1;
        cur_order = 0;
        //orderGap = m_minOverlap;
        orderGap = 0;
        
        atmosttimess = m_coverage/intervalGap;
    }
    
	//beginning kmer is a suffix of first read
    //initialize the beginning kmer SA intervals with kmer length=m_minOverlap
    std::string beginningkmer=pQuery->substr(m_currentLength-m_minOverlap); //起始kmer。substr(a-b)代表從字串的(a-b)位置開始，一直到字尾為子字串; ///now
    m_pRootNode->fwdInterval=BWTAlgorithms::findInterval( m_indices.pRBWT, reverse(beginningkmer));
    m_pRootNode->rvcInterval=BWTAlgorithms::findInterval( m_indices.pBWT, reverseComplement(beginningkmer));
    
    if(m_compressionLevel) ///mod1-1
    {
     beginningkmerRBWT.lower=m_pRootNode->fwdInterval.lower;   //mod7-1
     beginningkmerRBWT.upper=m_pRootNode->fwdInterval.upper;
     beginningkmerBWT.lower=m_pRootNode->rvcInterval.lower;
     beginningkmerBWT.upper=m_pRootNode->rvcInterval.upper;
    
     beginningkmerOrder=0;
    
     intervalP.FBintervalRBWTx = beginningkmerRBWT;
     intervalP.RBintervalBWTx = beginningkmerBWT;
     
     //std::cout<<"iSize: "<<beginningkmerRBWT.upper - beginningkmerRBWT.lower<<std::endl;
     //std::cout<<"iSize: "<<beginningkmerBWT.upper - beginningkmerBWT.lower<<std::endl;
     
     size_t timess=0; //正股起始檢查
     for(int64_t i=(int64_t)(m_pRootNode->fwdInterval.lower); i<(int64_t)(m_pRootNode->fwdInterval.upper) && timess<atmosttimess; i+=intervalGap) //檢查正股起始kmer的interval是否包含於已存之interval
     {
        if((*m_intervalMatchMapRBWT).find(immAccessor,i))
        {
            fwdReadsExisted[0] = true;
            int sizeVector = immAccessor->second.size();
            for(int p=0;p<sizeVector;p++)
            {
                fwdBeginningIdOrder.push_back(immAccessor->second.at(p));
            }
        }
        immAccessor.release();
        timess++;
     }
     //if(fwdBeginningIdOrder.size()>10) std::cout<<fwdBeginningIdOrder.size()<<std::endl;
     
     timess=0; //反股起始檢查
     for(int64_t i=(int64_t)(m_pRootNode->rvcInterval.lower); i<(int64_t)(m_pRootNode->rvcInterval.upper) && timess<atmosttimess; i+=intervalGap) //檢查反股起始kmer的interval是否包含於已存之interval
     {
        if((*m_intervalMatchMapBWT).find(immAccessor,i))
        {
            rvcReadsExisted[0] = true;
            int sizeVector = immAccessor->second.size();
            for(int p=0;p<sizeVector;p++)
            {
                rvcBeginningIdOrder.push_back(immAccessor->second.at(p));
            }
        }
        immAccessor.release();
        timess++;
     }
     //if(rvcBeginningIdOrder.size()>10) std::cout<<rvcBeginningIdOrder.size()<<std::endl;
    }

	//ending kmer is a prefix of second read
    //initialize the ending SA intervals with kmer length=m_minOverlap
    std::string endingkmer=secondread.substr(0,m_minOverlap); //終點kmer，secondread.substr(a,b)代表抓secondread從位置a處往後b長度的子字串
    m_fwdTerminatedInterval=BWTAlgorithms::findInterval( m_indices.pRBWT, reverse(endingkmer));
    m_rvcTerminatedInterval=BWTAlgorithms::findInterval( m_indices.pBWT, reverseComplement(endingkmer));
    
    if(m_compressionLevel)///mod1-2
    {
     endingkmerRBWT.lower=m_fwdTerminatedInterval.lower; //mod7-2-1
     endingkmerRBWT.upper=m_fwdTerminatedInterval.upper;
     endingkmerBWT.lower=m_rvcTerminatedInterval.lower;
     endingkmerBWT.upper=m_rvcTerminatedInterval.upper;
    
     intervalP.FEintervalRBWTx = endingkmerRBWT;
     intervalP.REintervalBWTx = endingkmerBWT;
     
     //std::cout<<"iSize: "<<endingkmerRBWT.upper - endingkmerRBWT.lower<<std::endl;
     //std::cout<<"iSize: "<<endingkmerBWT.upper - endingkmerBWT.lower<<std::endl;
     
     size_t timess=0; //正股終點檢查
     for(int64_t i=(int64_t)(m_fwdTerminatedInterval.lower); i<(int64_t)(m_fwdTerminatedInterval.upper) && timess<atmosttimess; i+=intervalGap) //檢查正股終點kmer的interval是否包含於已存之interval
     {
        if((*m_intervalMatchMapRBWT).find(immAccessor,i))
        {
            fwdReadsExisted[1] = true;
            int sizeVector = immAccessor->second.size();
            for(int p=0;p<sizeVector;p++)
            {
                fwdEndingIdOrder.push_back(immAccessor->second.at(p));
            }
        }
        immAccessor.release();
        timess++;
     }
     //if(fwdEndingIdOrder.size()>10) std::cout<<fwdEndingIdOrder.size()<<std::endl;
     
     timess=0; //反股終點檢查
     for(int64_t i=(int64_t)(m_rvcTerminatedInterval.lower); i<(int64_t)(m_rvcTerminatedInterval.upper) && timess<atmosttimess; i+=intervalGap) //檢查反股終點kmer的interval是否包含於已存之interval
     {
        if((*m_intervalMatchMapBWT).find(immAccessor,i))
        {
            rvcReadsExisted[1] = true;
            int sizeVector = immAccessor->second.size();
            for(int p=0;p<sizeVector;p++)
            {
                rvcEndingIdOrder.push_back(immAccessor->second.at(p));
            }
        }
        immAccessor.release();
        timess++;
     }
     //if(rvcEndingIdOrder.size()>10) std::cout<<rvcEndingIdOrder.size()<<std::endl;
    }
    
    if(m_compressionLevel) ///mod2
    {
     
     if(fwdReadsExisted[0] && fwdReadsExisted[1] && uniqueCase) //正正
     {
        //std::cout<<"3.FF\n";
        FBsize = fwdBeginningIdOrder.size();
        FEsize = fwdEndingIdOrder.size();
        //std::cout<<m_pQueryId<<"'s fwdReadsExisted[0] and fwdReadsExisted[1] are true.\n";
        for(int i=0; i<FBsize && uniqueCase; i++)
            for(int j=0; j<FEsize && uniqueCase; j++)
            {
                FBid = fwdBeginningIdOrder.at(i).first;
                FEid = fwdEndingIdOrder.at(j).first;
                FBorder = fwdBeginningIdOrder.at(i).second;
                FEorder = fwdEndingIdOrder.at(j).second;
                //std::cout<<longReadNames_FwdBeginning.at(i)<<"\n";
                //std::cout<<longReadNames_FwdEnding.at(j)<<"\n\n";
                if((FBid == FEid) && (FBorder+orderGap)<=FEorder)
                {
                    uniqueCase=0;
                    //std::cout<<"["<<m_pQueryId<<"] was discarded by long read ["<<FBid<<"]"<<"FF\n";
                }
            }
     }
     
     if(rvcReadsExisted[0] && rvcReadsExisted[1] && uniqueCase) //反反
     {
        //std::cout<<"4.RR\n";
        RBsize = rvcBeginningIdOrder.size();
        REsize = rvcEndingIdOrder.size();
        //std::cout<<m_pQueryId<<"'s rvcReadsExisted[0] and rvcReadsExisted[1] are true.\n";
        for(int i=0; i<RBsize && uniqueCase; i++)
            for(int j=0; j<REsize && uniqueCase; j++)
            {
                RBid = rvcBeginningIdOrder.at(i).first;
                REid = rvcEndingIdOrder.at(j).first;
                RBorder = rvcBeginningIdOrder.at(i).second;
                REorder = rvcEndingIdOrder.at(j).second;
                //std::cout<<longReadNames_RvcBeginning.at(i)<<"\n";
                //std::cout<<longReadNames_RvcEnding.at(j)<<"\n\n";
                if((RBid == REid) && (RBorder+orderGap)<=REorder)
                {
                    uniqueCase=0;
                    //std::cout<<"["<<m_pQueryId<<"] was discarded by long read ["<<RBid<<"]"<<"RR\n";
                }
            }
     }
     
     if(fwdReadsExisted[0] && rvcReadsExisted[1] && uniqueCase) //正反
     {
        //std::cout<<"1.FR\n";
        FBsize = fwdBeginningIdOrder.size();
        REsize = rvcEndingIdOrder.size();
        //std::cout<<m_pQueryId<<"'s fwdReadsExisted[0] and fwdReadsExisted[1] are true.\n";
        for(int i=0; i<FBsize && uniqueCase; i++)
            for(int j=0; j<REsize && uniqueCase; j++)
            {
                FBid = fwdBeginningIdOrder.at(i).first;
                REid = rvcEndingIdOrder.at(j).first;
                FBorder = fwdBeginningIdOrder.at(i).second;
                REorder = rvcEndingIdOrder.at(j).second;
                //std::cout<<longReadNames_FwdBeginning.at(i)<<"\n";
                //std::cout<<longReadNames_FwdEnding.at(j)<<"\n\n";
                if((FBid == REid) && (FBorder+orderGap)<=REorder)
                {
                    uniqueCase=0;
                    //std::cout<<"["<<m_pQueryId<<"] was discarded by long read ["<<FBid<<"]"<<"FR\n";
                }
            }
     }
     
     if(rvcReadsExisted[0] && fwdReadsExisted[1] && uniqueCase) //反正
     {
        //std::cout<<"2.RF\n";
        RBsize = rvcBeginningIdOrder.size();
        FEsize = fwdEndingIdOrder.size();
        //std::cout<<m_pQueryId<<"'s rvcReadsExisted[0] and rvcReadsExisted[1] are true.\n";
        for(int i=0; i<RBsize && uniqueCase; i++)
            for(int j=0; j<FEsize && uniqueCase; j++)
            {
                RBid = rvcBeginningIdOrder.at(i).first;
                FEid = fwdEndingIdOrder.at(j).first;
                RBorder = rvcBeginningIdOrder.at(i).second;
                FEorder = fwdEndingIdOrder.at(j).second;
                //std::cout<<longReadNames_RvcBeginning.at(i)<<"\n";
                //std::cout<<longReadNames_RvcEnding.at(j)<<"\n\n";
                if((RBid == FEid) && (RBorder+orderGap)<=FEorder)
                {
                    uniqueCase=0;
                    //std::cout<<"["<<m_pQueryId<<"] was discarded by long read ["<<RBid<<"]"<<"RF\n";
                }
            }
     }
     
    }
    
	//std::cout << m_minOverlap << ":" << beginningkmer << ":" << endingkmer << "\n";
	//getchar();
}

//
SAIntervalTree::~SAIntervalTree()
{
    // Recursively destroy the tree
    delete m_pRootNode;
}

//On success return the length of merged string
int SAIntervalTree::mergeTwoReads(std::string &mergedseq)
{
    SAIntervalNodeResultVector results;
	
    if( isTwoReadsOverlap(mergedseq)) //若起點跟終點overlap,直接合併
    {
        //mod9
        (*m_resultthis).thread = mergedseq;
        
        (*m_resultthis).pathIntervalRBWT[0].push_back(beginningkmerRBWT.lower);
        (*m_resultthis).pathIntervalRBWT[1].push_back(beginningkmerRBWT.upper);
        (*m_resultthis).pathIntervalBWT[0].push_back(beginningkmerBWT.lower);
        (*m_resultthis).pathIntervalBWT[1].push_back(beginningkmerBWT.upper);
        (*m_resultthis).order.push_back(0);
        
        (*m_resultthis).pathIntervalRBWT[0].push_back(endingkmerRBWT.lower);
        (*m_resultthis).pathIntervalRBWT[1].push_back(endingkmerRBWT.upper);
        (*m_resultthis).pathIntervalBWT[0].push_back(endingkmerBWT.lower);
        (*m_resultthis).pathIntervalBWT[1].push_back(endingkmerBWT.upper);
        (*m_resultthis).order.push_back(1);
        
		return 1;
    }

	//BFS search from 1st to 2nd read via FM-index walk
    while(!m_leaves.empty() && m_leaves.size() <= m_MaxLeaves && m_currentLength <=m_MaxLength)
    {
        // ACGT-extend the leaf nodes via updating existing SA interval
        extendLeaves();
				
		if(m_leaves.size()>m_maxUsedLeaves)	 m_maxUsedLeaves=m_leaves.size();
		// std::cout << m_currentKmerSize << ":" << m_currentLength << ":" << m_leaves.size() << "\n";	

		//see if terminating string is reached
		if(isTerminated(results))
			break;
    }
		
	//find the path with maximum kmer coverage
	if( results.size()>0 )
	{
		//if multiple paths are bubbles collapsing all together at terminal
		if(results.size()==m_leaves.size()) 
		{
			// std::cout << m_maxUsedLeaves << "\t" << results.size() <<"\n" << results[0].thread << "\n";
			m_isBubbleCollapsed=true;
		}
		std::string tmpseq;
        int j=-1; //mod5-1
		for (size_t i = 0 ; i < results.size() ;i++)
		{	//bug fix: m_secondread may be shorter than m_minOverlap
			if(m_secondread.length()>m_minOverlap)
				tmpseq=results[i].thread+m_secondread.substr(m_minOverlap);
			else
				tmpseq=results[i].thread;
			
			size_t cov = calculateKmerCoverage (tmpseq, m_minOverlap, m_indices.pBWT); //計算這條long read路徑上每kmerSize/2段的intervalSize總和
			// size_t cov=results[i].SAICoverage;
			if (  cov > m_maxKmerCoverage )
			{
				mergedseq=tmpseq;
				m_maxKmerCoverage=cov;
                j=i; //mod5-2
			}
		}
		
        if(m_compressionLevel) //mod4
        {
            (*m_resultthis) = results.at(j);
            /*
            int sizeInterval = results[j].pathInterval[0].size();
            for(int k=0;k<sizeInterval;k++)
            {
                //int strlength = m_pQueryId.length();
                int64_t t_interval=results[j].pathInterval[0].at(k);
                (*m_intervalMatchMap).insert(immAccessor,t_interval);
                    immAccessor->second.push_back(m_pQueryId);
                immAccessor.release();
            }
            */
        }
        
        
		// std::cout << ">\n" << *m_pQuery << "\n>\n" << reverseComplement(m_secondread);
		// std::cout << mergedseq.length() << "\t" << m_maxKmerCoverage <<  "\t" << (double)m_maxKmerCoverage/mergedseq.length()  << "\n";
		return 1;
    }

	// if(m_leaves.size() > 512){
		// std::cout << m_leaves.size() << "\t" << m_currentLength<< "\t" << m_currentLength-m_pQuery->length()<< "\n";
		// printAll();
		// getchar();
	// }
	
    //Did not reach the terminal kmer
    if(m_leaves.empty())
        return -1;	//high error
    else if(m_currentLength>m_MaxLength)
        return -2;	//exceed search depth
    else if(m_leaves.size() > m_MaxLeaves)
        return -3;	//too much repeats
	else
		return -4;
}

// Print the string represented by every node
void SAIntervalTree::printAll()
{
    std::cout << "Print all: \n";
    m_pRootNode->printAllStrings("");
}

// Extend each leaf node

void SAIntervalTree::attempToExtend(STNodePtrList &newLeaves) //newLeaves 代表某層需要被篩選的葉子 
{
    for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
        std::vector< std::pair<std::string, BWTIntervalPair> > extensions;
        extensions = getFMIndexExtensions(*iter); //

        // Either extend the current node or branch it
        // If no extension, do nothing and this node
        // is no longer considered a leaf
        if(extensions.size() == 1)
        {
            // Single extension, do not branch
            (*iter)->extend(extensions.front().first);
            (*iter)->fwdInterval=extensions.front().second.interval[0]; //fwdInterval.lower & fwdInterval.upper
            (*iter)->rvcInterval=extensions.front().second.interval[1]; //rvcInterval.lower & rvcInterval.upper
			if((*iter)->fwdInterval.isValid()) (*iter)->addKmerCount( (*iter)->fwdInterval.size());
			if((*iter)->rvcInterval.isValid()) (*iter)->addKmerCount( (*iter)->rvcInterval.size());
			
			if(m_compressionLevel) ///mod3-1
            {
                (*iter)->m_pathIntervalRBWT[0].push_back((*iter)->fwdInterval.lower);//mod3-1
                (*iter)->m_pathIntervalRBWT[1].push_back((*iter)->fwdInterval.upper);
                (*iter)->m_pathIntervalBWT[0].push_back((*iter)->rvcInterval.lower);
                (*iter)->m_pathIntervalBWT[1].push_back((*iter)->rvcInterval.upper);
                (*iter)->m_order.push_back(cur_order);
            }
            
            newLeaves.push_back(*iter);
        }
        else if(extensions.size() > 1)
        {
            // Branch
            for(size_t i = 0; i < extensions.size(); ++i)
            {
                SAIntervalNode* pChildNode = (*iter)->createChild(extensions[i].first);
                pChildNode->fwdInterval=extensions[i].second.interval[0]; //fwdInterval.lower & fwdInterval.upper
                pChildNode->rvcInterval=extensions[i].second.interval[1]; //rvcInterval.lower & rvcInterval.upper
				//inherit accumulated kmerCount from parent
				pChildNode->addKmerCount( (*iter)->getKmerCount() );
				if(pChildNode->fwdInterval.isValid()) pChildNode->addKmerCount( pChildNode->fwdInterval.size());
				if(pChildNode->rvcInterval.isValid()) pChildNode->addKmerCount( pChildNode->rvcInterval.size());
                
                if(m_compressionLevel) ///mod3-2
                {
                    pChildNode->m_pathIntervalRBWT[0].push_back(pChildNode->fwdInterval.lower);
                    pChildNode->m_pathIntervalRBWT[1].push_back(pChildNode->fwdInterval.upper);
                    pChildNode->m_pathIntervalBWT[0].push_back(pChildNode->rvcInterval.lower);
                    pChildNode->m_pathIntervalBWT[1].push_back(pChildNode->rvcInterval.upper);
                    pChildNode->m_order.push_back(cur_order);
                }
                
                newLeaves.push_back(pChildNode);
            }
        }
    }	
}

void SAIntervalTree::extendLeaves()
{
    STNodePtrList newLeaves;
	
	//attempt to extend one base for each leave
    attempToExtend(newLeaves);
	
    //shrink the SAIntervals in case overlap is larger than read length, which lead to empty newLeaves
    if(!m_kmerMode  &&  newLeaves.empty() )
    {
        refineSAInterval(m_minOverlap);
        attempToExtend(newLeaves);
    }	
	
	//extension succeed
    if(!newLeaves.empty()){
		m_currentKmerSize++;
        m_currentLength++;
	}

    m_leaves.clear();
    m_leaves = newLeaves;

	if(!m_leaves.empty() && (m_kmerMode || m_currentKmerSize >= m_maxOverlap) )
		refineSAInterval(m_minOverlap);

}

// Check for leaves whose extension has terminated. If the leaf has
// terminated, the walked string and coverage is pushed to the result vector
bool SAIntervalTree::isTerminated(SAIntervalNodeResultVector& results)
{
	bool found = false;

    for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
        BWTInterval currfwd=(*iter)->fwdInterval;
        BWTInterval currrvc=(*iter)->rvcInterval;

        // assert(currfwd.isValid() || currrvc.isValid());

		//If terminating kmer is a substr, the current SA interval is a sub-interval of the terminating interval
        bool isFwdTerminated = currfwd.isValid() && currfwd.lower >= m_fwdTerminatedInterval.lower
                            && currfwd.upper <= m_fwdTerminatedInterval.upper;
														
        bool isRvcTerminated=currrvc.isValid() && currrvc.lower >= m_rvcTerminatedInterval.lower
                            && currrvc.upper <= m_rvcTerminatedInterval.upper;

        if(isFwdTerminated || isRvcTerminated)
        {
            std::string STNodeStr = (*iter)->getFullString();
            SAIntervalNodeResult STresult;
            STresult.thread=STNodeStr;
			STresult.SAICoverage=(*iter)->getKmerCount();
			
			endingkmerOrder=cur_order+1; //mod7-2-2
            
            STresult.pathIntervalBWT[0].push_back(beginningkmerBWT.lower); //mod8-1
            STresult.pathIntervalBWT[1].push_back(beginningkmerBWT.upper);
            STresult.pathIntervalRBWT[0].push_back(beginningkmerRBWT.lower);
            STresult.pathIntervalRBWT[1].push_back(beginningkmerRBWT.upper);
            STresult.order.push_back(beginningkmerOrder);
            
            STresult.pathIntervalBWT[0].push_back(endingkmerBWT.lower); //mod8-2
            STresult.pathIntervalBWT[1].push_back(endingkmerBWT.upper);
            STresult.pathIntervalRBWT[0].push_back(endingkmerRBWT.lower);
            STresult.pathIntervalRBWT[1].push_back(endingkmerRBWT.upper);
            STresult.order.push_back(endingkmerOrder);
            
            while((*iter)!=NULL) //mod6
            {
                int sizeIntervalBWT = (*iter)->m_pathIntervalBWT[0].size(); ///mod3-1
                for(int i=0;i<sizeIntervalBWT;i++)
                {
                    STresult.pathIntervalBWT[0].push_back((*iter)->m_pathIntervalBWT[0].at(i));
                    STresult.pathIntervalBWT[1].push_back((*iter)->m_pathIntervalBWT[1].at(i));
                    STresult.pathIntervalRBWT[0].push_back((*iter)->m_pathIntervalRBWT[0].at(i));
                    STresult.pathIntervalRBWT[1].push_back((*iter)->m_pathIntervalRBWT[1].at(i));
                    STresult.order.push_back((*iter)->m_order.at(i));
                }
                /*
                (*iter)->m_pathIntervalBWT[0].clear(); ///mod4-1-1
                (*iter)->m_pathIntervalBWT[1].clear();
                (*iter)->m_pathIntervalRBWT[0].clear(); ///mod4-1-2
                (*iter)->m_pathIntervalRBWT[1].clear();
                (*iter)->m_order.clear(); ///mod4-1-3
                */
                (*iter) = (*iter)->m_pParent;
            }

            //compute the merged pos right next to the kmer on 2nd read.
            results.push_back(STresult);
            found =  true;
            /*
            STresult.pathIntervalBWT[0].clear(); ///mod4-2
            STresult.pathIntervalBWT[1].clear();
            STresult.pathIntervalRBWT[0].clear();
            STresult.pathIntervalRBWT[1].clear();
            STresult.order.clear();
            */
        }
    }

    return found;
}

bool SAIntervalTree::isTwoReadsOverlap(std::string & mergedseq)
{
    //case 1: 1st read sense overlap to 2nd read at exact m_minOverlap bases
    if(BWTInterval::equal(m_pRootNode->fwdInterval, m_fwdTerminatedInterval))
    {
        mergedseq= (*m_pQuery)+m_secondread.substr(m_minOverlap);
        return true;
    }

    //case 2: 1st read sense overlap 2nd read
    std::string secondLeftKmer=m_secondread.substr(0,m_minOverlap);
	//assume overlap can't exceed 100 bp
    size_t pos=m_pQuery->find(secondLeftKmer, m_pQuery->length()>=200?m_pQuery->length()-200:0);
    if(pos!=std::string::npos)	
    {
		//make sure entire suffix of 1st read after pos matches the prefix of 2nd read
		if( m_pQuery->substr(pos) == m_secondread.substr(0, m_pQuery->length()-pos) )
		{
			mergedseq=m_pQuery->substr(0,pos)+m_secondread;
			return true;
		}
    }

    //case 3: 1st read antisense overlap with 2nd read, or 1st read is substr of 2nd read
	//This is rare case and we don't do this in m_kmerMode during island joint
	if(m_kmerMode) return false;
    std::string firstLeftKmer=m_pQuery->substr(0,m_minOverlap);
    pos=m_secondread.find(firstLeftKmer);
	//assume antisense overlap can't exceed 50bp due to rare cases
    if(pos!=std::string::npos && pos <=50)
    {
        //make sure entire suffix of 2nd read after pos matches the prefix of 1st read
		if( m_secondread.substr(pos) ==  m_pQuery->substr(0, m_secondread.length()-pos))
		{
			//return overlapped portion
			mergedseq=m_secondread.substr(pos);
			return true;
		}
    }

    return false;

}

//update SA intervals of each leaf, which corresponds to one-base extension
std::vector<std::pair<std::string, BWTIntervalPair> > SAIntervalTree::getFMIndexExtensions(SAIntervalNode* pNode)
{
    std::vector<std::pair<std::string, BWTIntervalPair> > out;

    for(int i = 1; i < BWT_ALPHABET::size; ++i) //i=A,C,G,T
    {
        char b = BWT_ALPHABET::getChar(i);

        //update forward Interval using extension b
        BWTInterval fwdProbe=pNode->fwdInterval;
        if(fwdProbe.isValid())
            BWTAlgorithms::updateInterval(fwdProbe, b, m_indices.pRBWT);

        //update reverse complement Interval using extension rcb
        BWTInterval rvcProbe=pNode->rvcInterval;
		char rcb=BWT_ALPHABET::getChar(5-i); //T,G,C,A
        if(rvcProbe.isValid())
            BWTAlgorithms::updateInterval(rvcProbe, rcb, m_indices.pBWT);

        size_t bcount = 0;
        if(fwdProbe.isValid())
            bcount += fwdProbe.size();
        if(rvcProbe.isValid())
            bcount += rvcProbe.size();
		
		//min freq at fwd and rvc bwt
        if(bcount >= m_min_SA_threshold) //超過門檻值,保留進 out vector
        {
			// if(bcount>50)
				// std::cout << m_currentKmerSize << ":" << bcount <<"\n";
            // extend to b
            std::string tmp;
            tmp.append(1,b);
            BWTIntervalPair bip;
            bip.interval[0]=fwdProbe;
            bip.interval[1]=rvcProbe;
            out.push_back(std::make_pair(tmp, bip));
        }
    }// end of ACGT
    cur_order++;
    
    return out;
}

size_t SAIntervalTree::calculateKmerCoverage (const std::string & seq , size_t kmerLength , const BWT* pBWT) //seq為被篩選的long reads群
{
	if (seq.length() < kmerLength) return 0;

	size_t cov = 0 ;
	for (size_t i=0; i<=seq.length()-kmerLength;i+=kmerLength/2)
		cov += BWTAlgorithms::countSequenceOccurrences(seq.substr(i,kmerLength) , pBWT ); //居然TM是這樣計算的
	
	return cov;
}

// replace each kmer with highest one at each locus
bool SAIntervalTree::replaceLowFreqKmer (std::string & seq , size_t kmerLength)
{
	bool changed = false;
	
	for (size_t i=0; i <=seq.length()-kmerLength; i++)
	{
		//Forward kmer should be computed reversely using pRBWT for forward extension
		BWTInterval fwdProbe=BWTAlgorithms::findInterval(m_indices.pRBWT, reverse(seq.substr(i, kmerLength-1)));
		BWTInterval rvcProbe=BWTAlgorithms::findInterval(m_indices.pBWT, reverseComplement(seq.substr(i, kmerLength-1)));
		
		size_t maxcov=0;
		for(int j = 1; j < BWT_ALPHABET::size; ++j) //j=A,C,G,T
		{
			char b = BWT_ALPHABET::getChar(j);

			//update forward Interval using extension b
			if(fwdProbe.isValid())
				BWTAlgorithms::updateInterval(fwdProbe, b, m_indices.pRBWT);

			//update reverse complement Interval using extension rcb
			char rcb=BWT_ALPHABET::getChar(5-i); //T,G,C,A
			if(rvcProbe.isValid())
				BWTAlgorithms::updateInterval(rvcProbe, rcb, m_indices.pBWT);

			size_t bcount = 0;
			if(fwdProbe.isValid())
				bcount += fwdProbe.size();
			if(rvcProbe.isValid())
				bcount += rvcProbe.size();

			if(bcount > maxcov) {
				maxcov = bcount;
				seq.replace(i+kmerLength-1, 1, 1, b);
				changed = true;
			}
		}
	}
	
	return changed;
}

// Refine SA intervals of each leave with a new kmer
void SAIntervalTree::refineSAInterval(size_t newKmerSize) //extension長到max_overlap後，縮回min_overlap長度
{
	assert(m_currentLength >= newKmerSize);

    for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
        // reset the SA intervals using original m_minOverlap
        std::string pkmer = (*iter)->getSuffix(newKmerSize);
		(*iter)->fwdInterval=BWTAlgorithms::findInterval(m_indices.pRBWT, reverse(pkmer));
		(*iter)->rvcInterval=BWTAlgorithms::findInterval(m_indices.pBWT, reverseComplement(pkmer));
    }

	m_currentKmerSize=newKmerSize;
}

/***Dead code***/

// Remove leaves with two or more same kmers
void SAIntervalTree::removeLeavesByRepeatKmer()
{
    STNodePtrList newLeaves;

    for(STNodePtrList::iterator iter = m_leaves.begin(); iter != m_leaves.end(); ++iter)
    {
        /*
        GAGGCAGTTGAGGCAGTTGAGGCAGTTGAGGCAGTTGAGGCAGTTGAGGCAGTTGAGGCAGTTGAGGCAGTTGAGGCAGTTGAGGCAGTTGAGGCAGTTG
        TGGATTCCAGATTGTTCGAGGAGAATTTGGTGGAGCTACGCGGGATCGAACCGCGGACCTCTTGCATGCCATGCAAGCGCTCTCCCAGCTGAGCTATAACCCCTTGGATTCCAGATTGTTCGAGGAGAATTTGGTGGAGCTACGCGGGATCGAACCGCGGACCTCTTGCATGCCATGCAAGCGCTCTCCCAGCTGAGCTATAACC
        GAGAGGGACTCGAACCCTCACACCCGGGGGGCACTAACACCTGAAGCTAGCGCGTCTACCAATTCCGCCACCTTCGCACATCGGGTTATCAGTCTGGATTTACATGCTGTCTGATAAAAGCATGGTGCGAAGAGAGGGACTCGAACCCTCACACCCGGGGGGCACTAACACCTGAAGCTAGCGCGTCTACCAATTCCGCCACCTTCGCACATCGGGTTATCAGTCTGGATTT
        GCATATCCATCCCACCAGCACATCGACCTATCGACTCATCAGTTCATCAGTTCATCAGTTCATCAGTTCATCAGTTCATCAGTTCATCAGTTCATCAGTTCATCAGTTCATCAGTTCATCAGTTCATCAGTTCATCAGTTCATC
        CATCGGCGTCAGCCTGCTGGGCTTCACCCATCAGGGCAACAAGTGGCTGTGGCAGCAGGCCAGGGCCGCTCTTCCCTCCCTCAAGGGGGAGCTGGTGGCGGGGGGGGGGGGGGGGGGGGGGGGGGCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
        */
        std::string STNodeStr = (*iter)->getFullString();
        std::string fwdrepeatunit = STNodeStr.substr(STNodeStr.size()-m_minOverlap);
        std::string revrepeatunit = reverseComplement(fwdrepeatunit);
        size_t index1=STNodeStr.find(fwdrepeatunit);
        size_t index2=STNodeStr.find(revrepeatunit);

        if(index1 == (STNodeStr.size()- m_minOverlap) && index2 == std::string::npos)
        {
            newLeaves.push_back(*iter);
        }
    }

    m_leaves=newLeaves;
}

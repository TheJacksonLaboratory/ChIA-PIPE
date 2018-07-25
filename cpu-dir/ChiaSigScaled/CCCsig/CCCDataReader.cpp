// ChiaSig: detection of significant interactions in ChIA-PET data using Non-central hypergeometric distribution (NCHG)
// This file is subject to the terms and conditions defined in file 'LICENSE', which is part of this source code package.
// Copyright (2014) Jonas Paulsen

#include <set>
#include <fstream>
#include <sstream>
#include <iostream>
#include <assert.h>
#include "CCCDataReader.h"
#include "Segment.h"

using namespace std;

CCCDataReader::CCCDataReader(string file, map<ChromosomeIndexType, CCCMatrixInteraction>& cm)
: contactMatrices(cm), filePath(file)
{
	// The constructor determines the unique segments, and nothing more.
	// These are needed for all downstream analyses.
	// Currently, this is only based on intra-data.
	
	ifstream inFile(filePath.c_str());
	if (!inFile.good()) {
		cerr << "Problem opening file " << filePath.c_str() << endl;
		assert(inFile.good());
	}
	string chrL, chrR;
	ChromosomeIndexType nChrL, nChrR;
	int startL, endL, startR, endR, nij;
	
	// Get information about (unique) Segments:
	ChromosomeIndexType nIndex = 0;
	map<string, ChromosomeIndexType >::iterator it_chromToIndex;
	
	map<ChromosomeIndexType, set<Segment> > segs;
	while (inFile >> chrL >> startL >> endL >> chrR >> startR >> endR >> nij) {
		assert(chrL == chrR); // Only support intra-data at the moment.
		
		it_chromToIndex = chromToIndex.find(chrL);
		if (it_chromToIndex==chromToIndex.end()) {
			chromToIndex[chrL] = nIndex;
			chromosomes.push_back(chrL);
			nChrL = nIndex;
			nIndex++;
		} else {
			nChrL = it_chromToIndex->second;
		}
		it_chromToIndex = chromToIndex.find(chrR);
		if (it_chromToIndex==chromToIndex.end()) {
			chromToIndex[chrR] = nIndex;
			chromosomes.push_back(chrR);
			nChrR = nIndex;
			nIndex++;
		} else {
			nChrR = it_chromToIndex->second;
		}
		
		Segment segL(nChrL, startL, endL);
		Segment segR(nChrR, startR, endR);
		
		segs[nChrL].insert(segL);
		segs[nChrR].insert(segR);
		
	}
	
	// Loop through chromosomes, and build empty matrices:
	set<Segment>::size_type totalAnchor = 0;
	Interaction interaction(0,0,1.0);
	for(map<ChromosomeIndexType, set<Segment> >::iterator it = segs.begin(); it != segs.end(); it++) {
		//CCCMatrixInteraction mat(it->second, it->second, interaction, true, chromosomes[it->first]);
		//contactMatrices[it->first] = mat;
		contactMatrices.insert(pair<ChromosomeIndexType, CCCMatrixInteraction >(it->first, CCCMatrixInteraction()));
		contactMatrices[it->first].LoadInteractions(it->second, it->second, interaction, true, chromosomes[it->first]);
		totalAnchor += it->second.size();
		fprintf(stderr, "[M::%s:%d] #anchor(%s)=%lu\n", __func__, __LINE__, chromosomes[it->first].c_str(), it->second.size());
	}
	inFile.close();
	fprintf(stderr, "[M::%s:%d] #anchor(genome)=%lu\n", __func__, __LINE__, totalAnchor);
}


void CCCDataReader::buildContactMatrices() {
	ifstream inFile(filePath.c_str());
	string chrL, chrR;
	int startL, endL, startR, endR, nij;
	assert(inFile.good());
	while (inFile >> chrL >> startL >> endL >> chrR >> startR >> endR >> nij) {
		assert(chrL == chrR); // Only support intra-data at the moment.
		ChromosomeIndexType nChrL = chromToIndex.at(chrL);
		Segment segL(nChrL, startL, endL);
		Segment segR(chromToIndex.at(chrR), startR, endR);
		contactMatrices[nChrL].setElement(segL, segR, nij);
	}
	inFile.close();
	
	vector<string>::size_type totalInteractions = 0;
	for(vector<string>::size_type i=0; i<chromosomes.size(); ++i) {
		totalInteractions += contactMatrices[i].getInteractionCount();
		fprintf(stderr, "[M::%s:%d] #interaction(%s)=%lu\n", __func__, __LINE__, chromosomes[i].c_str(), contactMatrices[i].getInteractionCount());
	}
	fprintf(stderr, "[M::%s:%d] #interaction(genome)=%lu\n", __func__, __LINE__, totalInteractions);
}

CCCMatrixInteraction& CCCDataReader::getContactMatrix(ChromosomeIndexType chr) {
	return contactMatrices[chr];
}

void CCCDataReader::getChromosomes(vector<string>& res) {
	res.clear(); res.reserve(chromosomes.size());
	for(ChromosomeIndexType i = 0; i<chromosomes.size(); ++i) {
		res.push_back(chromosomes[i]);
	}
}

bool CCCDataReader::isChromosomePresent(string& chr) const {
	for(ChromosomeIndexType i = 0; i<chromosomes.size(); ++i) {
		if (chromosomes[i]==chr) return true;
	}
	return false;
}

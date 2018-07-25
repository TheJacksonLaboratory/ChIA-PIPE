// ChiaSig: detection of significant interactions in ChIA-PET data using Non-central hypergeometric distribution (NCHG)
// This file is subject to the terms and conditions defined in file 'LICENSE', which is part of this source code package.
// Copyright (2014) Jonas Paulsen

#ifndef GUARD_CCCDataReader
#define GUARD_CCCDataReader

#include <string>
#include <map>
#include <vector>

#include "CCCMatrix.h"

class CCCDataReader {
public:
	CCCDataReader(std::string, std::map<ChromosomeIndexType, CCCMatrixInteraction>&);
	void buildContactMatrices();
	CCCMatrixInteraction &getContactMatrix(ChromosomeIndexType);
	void getChromosomes(std::vector<std::string>&);
	bool isChromosomePresent(std::string& chr) const;

private:
	std::string filePath;
	std::map<ChromosomeIndexType, CCCMatrixInteraction >& contactMatrices;
	std::vector<std::string> chromosomes;
	std::map<std::string, ChromosomeIndexType > chromToIndex;
};

#endif

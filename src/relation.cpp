#include "relation.h"
#include "MurmurHash3.h"
#include "OT.h"
#include "OEP.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include "PRNG.h"
#include "PSI.h"
#include "GC.h"

using namespace std;

Relation::Relation(const Relation& R)
{
	numAttrs = R.numAttrs;
	annots1 = R.annots1;
	annots2 = R.annots2;
	annotsShared = R.annotsShared;
	tuplesPublic = R.tuplesPublic;
	annotsBool = R.annotsBool;
	attrNames = new string[numAttrs];
	copy(R.attrNames, R.attrNames + numAttrs, attrNames);
	attrTypes = new DataType[numAttrs]();
	copy(R.attrTypes, R.attrTypes + numAttrs, attrTypes);
	tuples.reserve(R.tuples.size());
	for (auto row : R.tuples)
	{
		uint64_t* thisRow = new uint64_t[numAttrs];
		copy(row, row + numAttrs, thisRow);
		tuples.push_back(thisRow);
	}
}

Relation::Relation(string* attrNames, int numAttrs, bool annotsBool)
{
	this->attrNames = new string[numAttrs];
	copy(attrNames, attrNames + numAttrs, this->attrNames);
	this->attrTypes = new DataType[numAttrs]();
	this->numAttrs = numAttrs;
	this->annotsShared = false;
	this->tuplesPublic = false;
	this->annotsBool = annotsBool;
}

uint64_t hashString(string str)
{
	char out[16];
	MurmurHash3_x64_128(str.c_str(), (int)str.length(), 0, out);
	return *((uint64_t*)out);
}

void Relation::loadData(const char* filePath, string annotName)
{
	unordered_map<string, int> inverseAttrMap;
	for (int i = 0; i < numAttrs; i++)
		inverseAttrMap[attrNames[i]] = i;

	// open file filePath
	ifstream fin(filePath, ios::in);


	if (!fin.is_open())
		throw "Cannot open file!";

	string element;
	int numColumns;
	fin >> numColumns;
	fin >> ws;// remove the newline character

	int* indexMap = new int[numColumns];
	for (int i = 0; i < numColumns; i++)
	{
		getline(fin, element, '|');
		auto it = inverseAttrMap.find(element);
		if (it != inverseAttrMap.end())
			indexMap[i] = it->second;
		else if (element == annotName)
			indexMap[i] = -2;
		else
			indexMap[i] = -1;
	}
	fin >> ws;

	for (int i = 0; i < numColumns; i++)
	{
		getline(fin, element, '|');
		int index = indexMap[i];
		if (index >= 0)
		{
			if (element == "int")
				attrTypes[index] = DataType::INT;
			else if (element == "decimal")
				attrTypes[index] = DataType::DECIMAL;
			else if (element == "date")
				attrTypes[index] = DataType::DATE;
			else if (element == "string")
				attrTypes[index] = DataType::STRING;
			else if (element == "char")
				attrTypes[index] = DataType::CHAR;
			else
				throw "Unknown data type!";
		}
	}
	fin >> ws;

	while (!fin.eof())
	{
		uint64_t* row = new uint64_t[numAttrs]; // annotation not included
		int year, month, day;
		float f_value;
		for (int i = 0; i < numColumns; i++)
		{
			getline(fin, element, '|');
			int index = indexMap[i];
			if (index >= 0)
			{
				switch (attrTypes[index])
				{
				case DataType::INT:
					row[index] = stoi(element);
					break;
				case DataType::DECIMAL:
					f_value = stof(element);
					row[index] = (uint64_t)(100 * f_value);
					break;
				case DataType::DATE:
					if (sscanf(element.c_str(), "%4d-%2d-%2d", &year, &month, &day) != 3)
						throw "Read date error!";
					row[index] = year * 10000 + month * 100 + day;
					break;
				case DataType::STRING:
					row[index] = hashString(element);
					break;
				case DataType::CHAR:
					row[index] = element[0];
					break;
				}
			}
			else if (index == -2) // annotation
				annots1.push_back(stoi(element));

		}
		tuples.push_back(row);
		fin >> ws;
	}
	fin.close();
	annots2.resize(annots1.size());
	delete[] indexMap;
}

void Relation::reveal()
{
	if (!annotsShared)
		return;
	for (int i = 0; i < annots1.size(); i++)
	{
		if (annotsBool)
		{
			addComm(1);
			annots1[i] ^= annots2[i];
		}
		else
		{
			addComm(16);
			annots1[i] += annots2[i];
		}
		annots2[i] = 0;
	}
	annotsShared = false;
}

void Relation::print(bool showRealAnnot, bool showZeroAnnotTuples, int limit_size)
{
	for (int i = 0; i < numAttrs; i++)
		std::cout << attrNames[i] << '\t';
	std::cout << "annotation\n";
	uint64_t i_value, year, month, day;
	float f_value;
	for (int i = 0; i < tuples.size(); ++i)
	{
		uint16_t realAnnot = annots1[i] + annots2[i];
		if (realAnnot == 0 && !showZeroAnnotTuples)
			continue;
		if (limit_size-- <= 0)
			break;
		for (int j = 0; j < numAttrs; ++j)
		{
			switch (attrTypes[j])
			{
			case DataType::INT:
			case DataType::STRING:
				std::cout << (int)tuples[i][j] << '\t';
				break;
			case DataType::DATE:
				i_value = tuples[i][j];
				day = i_value % 100;
				i_value /= 100;
				month = i_value % 100;
				i_value /= 100;
				year = i_value;
				std::cout << year << '-' << setfill('0') << setw(2) << month << '-' << setfill('0') << setw(2) << day << '\t';
				break;
			case DataType::DECIMAL:
				f_value = (float)tuples[i][j] / 100;
				std::cout << setprecision(2) << fixed << f_value << '\t';
				break;
			case DataType::CHAR:
				std::cout << (char)tuples[i][j] << '\t';
				break;
			}
		}
		uint16_t out = showRealAnnot ? realAnnot : annots1[i];
		std::cout << out << '\n';
	}
	cout << '\n';
}

void Relation::permute(int* permutedIndices)
{
	int numRows = tuples.size();
	auto origTuples = tuples;
	for (int i = 0; i < numRows; i++)
		tuples[i] = origTuples[permutedIndices[i]];

	uint16_t* Z1 = new uint16_t[numRows];
	copy(annots1.begin(), annots1.end(), Z1);
	for (int i = 0; i < numRows; i++)
		annots1[i] = Z1[permutedIndices[i]];
	if (annotsShared)
	{
		uint16_t* Z2 = new uint16_t[numRows];
		oblivPermutation(permutedIndices, (uint16_t*)&annots2[0], numRows, Z1, Z2);
		for (int i = 0; i < numRows; i++)
		{
			annots1[i] += Z1[i];
			annots2[i] = Z2[i];
		}
		delete[] Z2;
	}
	delete[] Z1;

}

void Relation::sort()
{
	int numRows = (int)tuples.size();
	int* permutedIndices = new int[numRows];
	for (int i = 0; i < numRows; i++)
		permutedIndices[i] = i;
	std::sort(permutedIndices, permutedIndices + numRows, [&](int i, int j) {
		for (int n = 0; n < numAttrs; n++)
		{
			if (tuples[i][n] < tuples[j][n])
				return true;
			else if (tuples[i][n] > tuples[j][n])
				return false;
		}
		return false;
		});
	permute(permutedIndices);
	delete[] permutedIndices;
}

void Relation::project(string* projectAttrNames, int numProjectAttrs)
{
	unordered_map<string, int> invAttrMap;
	int* indexMap = new int[numProjectAttrs];
	DataType* projectAttrTypes = new DataType[numProjectAttrs];
	for (int i = 0; i < numAttrs; i++)
		invAttrMap[attrNames[i]] = i;
	for (int i = 0; i < numProjectAttrs; i++)
	{
		auto it = invAttrMap.find(projectAttrNames[i]);
		if (it != invAttrMap.end())
			indexMap[i] = it->second;
		else
			throw "Project attribute error!";
	}
	for (int i = 0; i < tuples.size(); i++)
	{
		uint64_t* row = tuples[i];
		uint64_t* projectRow = new uint64_t[numProjectAttrs];
		for (int i = 0; i < numProjectAttrs; i++)
			projectRow[i] = row[indexMap[i]];
		tuples[i] = projectRow;
		delete[] row;
	}
	for (int i = 0; i < numProjectAttrs; i++)
		projectAttrTypes[i] = attrTypes[indexMap[i]];
	numAttrs = numProjectAttrs;
	delete[] indexMap;
	delete[] attrNames;
	delete[] attrTypes;
	attrTypes = projectAttrTypes;
	attrNames = new string[numAttrs];
	copy(projectAttrNames, projectAttrNames + numAttrs, attrNames);
}

void Relation::annotProject(string* projectAttrNames, int numProjectAttrs)
{
	if (annotsShared)
	{
		oblivAnnotProject(projectAttrNames, numProjectAttrs);
		return;
	}
	for (int i = 0; i < tuples.size(); i++)
		if (annots1[i] != 0)
			annots1[i] = 1;
	project(projectAttrNames, numProjectAttrs);
	sort();
	for (int i = 0; i < tuples.size() - 1; i++)
	{
		bool sameTuple = true;
		for (int j = 0; j < numAttrs; j++)
		{
			if (tuples[i][j] != tuples[i + 1][j])
			{
				sameTuple = false;
				break;
			}
		}
		if (sameTuple)
		{
			annots1[i + 1] |= annots1[i];
			annots1[i] = 0;
			for (int j = 0; j < numAttrs; j++) // replace with dummy tuple
				tuples[i][j] = gRNG.nextUInt64();
		}
	}
	annotsBool = true;
}

void Relation::oblivAnnotProject(string* projectAttrNames, int numProjectAttrs)
{
	project(projectAttrNames, numProjectAttrs);
	sort();
	int numRows = tuples.size();

	// Since the differences of time and cost on OEP between boolean and arithmetic sharing are not large, 
	// we do not implement OEP on boolean sharing. Therefore, we first run OEP and then do equality test
	bool* b1 = new bool[numRows];
	bool* b2 = new bool[numRows];
	bool* c = new bool[numRows - 1](); // the indicators showing whether two neighbor tuples are DIFFERENT

	for (int i = 0; i < numRows; i++)
		if (annotsBool)
		{
			b1[i] = annots1[i] & 1;
			b2[i] = annots2[i] & 1;
		}
		else
			equalityTest(annots1[i], -annots2[i], 16, b1[i], b2[i]);

	for (int i = 0; i < numRows - 1; i++)
	{
		for (int j = 0; j < numAttrs; j++)
		{
			if (tuples[i][j] != tuples[i + 1][j])
			{
				c[i] = true;
				break;
			}
		}
	}
	using namespace GC;
	vector<Gate*> gates;
	vector<Gate*> outGates;
	auto gb2 = new InGate * [numRows]; // the input gates of b2
	for (int i = 0; i < numRows; i++)
	{
		gb2[i] = new InGate(b2[i], false);
		gates.push_back(gb2[i]);
	}
	auto gb = new Gate * [numRows];
	for (int i = 0; i < numRows; i++)
	{
		gb[i] = new GHXGate(b1[i], gb2[i]);
		gates.push_back(gb[i]);
	}
	// b_i' = b_i & c_i
	// b_{i+1}'=(b_i ^ c_i ^ b_i') & b_{i+1}'
	for (int i = 0; i < numRows - 1; i++)
	{
		Gate* gbi = new GHAGate(c[i], gb[i]);
		gates.push_back(gbi);
		outGates.push_back(gbi);
		Gate* gbi2 = new GHXGate(c[i], gb[i]);
		gates.push_back(gbi2);
		gbi2 = new XORGate(gbi2, gbi);
		gates.push_back(gbi2);
		gbi2 = new ANDGate(gbi2, gb[i + 1]);
		gates.push_back(gbi2);
		gb[i] = gbi;
		gb[i + 1] = gbi2;
	}
	outGates.push_back(gb[numRows - 1]);

	vector<bool> z1, z2;
	GC::GC(gates, outGates, z1, z2);

	for (int i = 0; i < numRows - 1; i++)
	{
		if (c[i]) // i and i+1 are different
			z1[i] = !z1[i];
		else // replace with dummy
			for (int j = 0; j < numAttrs; j++)
				tuples[i][j] = gRNG.nextUInt64();
	}
	if (c[numRows - 2])
		z1[numRows - 1] = !z1[numRows - 1];
	copy(z1.begin(), z1.end(), annots1.begin());
	copy(z2.begin(), z2.end(), annots2.begin());
	annotsBool = true;

	for (auto gate : gates)
		delete gate;
	delete[] b1;
	delete[] b2;
	delete[] c;

}

void Relation::agg(string* groupByAttrNames, int numGroupByAttrs)
{
	if (annotsShared && !tuplesPublic)
	{
		oblivAgg(groupByAttrNames, numGroupByAttrs);
		return;
	}
	size_t numRows = tuples.size();
	project(groupByAttrNames, numGroupByAttrs);
	sort();
	for (int i = 0; i < numRows - 1; i++)
	{
		bool sameTuple = true;
		for (int j = 0; j < numAttrs; j++)
		{
			if (tuples[i][j] != tuples[i + 1][j])
			{
				sameTuple = false;
				break;
			}
		}
		if (sameTuple)
		{
			annots1[i + 1] += annots1[i];
			annots1[i] = 0;
			annots2[i + 1] += annots2[i];
			annots2[i] = 0;
			for (int j = 0; j < numAttrs; j++) // replace with dummy tuple
				tuples[i][j] = gRNG.nextUInt64();
		}
	}
	sort();
	annotsBool = false;
}

void Relation::oblivAgg(string* groupByAttrNames, int numGroupByAttrs)
{
	size_t numRows = tuples.size();
	project(groupByAttrNames, numGroupByAttrs);
	sort();

	uint16_t msg0[2];
	uint16_t msg1[2];
	uint16_t msgb[2];
	for (int i = 0; i < numRows - 1; i++)
	{
		// a,b,c,d and messages belong to Bob
		uint16_t a = annots2[i], b = annots2[i + 1];
		uint16_t c = gRNG.nextUInt16();
		uint16_t d = gRNG.nextUInt16();
		msg0[0] = a - c;
		msg0[1] = b - d;
		msg1[0] = -c;
		msg1[1] = a + b - d;
		bool sameTuple = 1;
		for (int j = 0; j < numAttrs; j++)
		{
			if (tuples[i][j] != tuples[i + 1][j])
			{
				sameTuple = 0;
				break;
			}
		}
		OT::OT((uint8_t*)msg0, (uint8_t*)msg1, 16 / 8 * 2, sameTuple, (uint8_t*)msgb);
		if (sameTuple)
		{
			annots1[i + 1] += annots1[i] + msgb[1];
			annots1[i] = msgb[0];
			for (int j = 0; j < numAttrs; j++) // replace with dummy tuple
				tuples[i][j] = gRNG.nextUInt64();
		}
		else
		{
			annots1[i] += msgb[0];
			annots1[i + 1] += msgb[1];
		}
		annots2[i] = c;
		annots2[i + 1] = d;
	}
	annotsShared = true;
	annotsBool = false;
}

uint64_t hashMultiColumns(uint64_t* values, int numColumns)
{
	if (numColumns <= 0)
		throw "Error num columns to hash!";
	if (numColumns == 1)
		return values[0];
	uint64_t out[2];
	MurmurHash3_x64_128(values, numColumns * 8, 0, out);
	return out[0];
}

void Relation::oblivSemiJoin(Relation& child, const unordered_map<string, string>& joinAttrNameMap)
{
	Relation myProjectRelation(*this);
	string* joinAttrNames = new string[child.numAttrs];
	for (int i = 0; i < child.numAttrs; i++)
		joinAttrNames[i] = joinAttrNameMap.find(child.attrNames[i])->second;
	myProjectRelation.project(joinAttrNames, child.numAttrs);
	delete[] joinAttrNames;
	int myRowNum = (int)tuples.size();
	int numJoinAttrs = child.numAttrs;

	// Alice 
	unordered_map<uint64_t, vector<int>> rowIndexMap;
	for (int i = 0; i < myRowNum; i++)
		rowIndexMap[hashMultiColumns(myProjectRelation.tuples[i], numJoinAttrs)].push_back(i);

	vector<uint64_t> myHashValues;
	vector<vector<int>> firstPermutedIndices;
	myHashValues.reserve(myRowNum);
	firstPermutedIndices.reserve(myRowNum);
	for (auto mapPair : rowIndexMap)
	{
		myHashValues.push_back(mapPair.first);
		firstPermutedIndices.push_back(mapPair.second);
	}
	for (int i = 0; i < myRowNum - firstPermutedIndices.size(); i++)
		myHashValues.push_back(gRNG.nextUInt64()); // append dummy values
	vector<int> secondPermutedIndices;
	vector<uint16_t> childAnnots1;

	// Bob
	vector<uint64_t> childHashValues;
	childHashValues.reserve(child.tuples.size());
	for (auto row : child.tuples)
		childHashValues.push_back(hashMultiColumns(row, numJoinAttrs));
	vector<uint16_t> childAnnots2;

	// PSI
	PSI psi(myHashValues, childHashValues);
	if (child.annotsShared)
		psi.sendSharedPayloads(child.annots1, child.annots2, childAnnots1, childAnnots2);
	else
		psi.sendPayloads(child.annots1, childAnnots1, childAnnots2);
	psi.getAliceIndicesHashed(secondPermutedIndices);
	int bucketSize = (int)secondPermutedIndices.size();
	int* invTotalPermutedIndices = new int[myRowNum];
	for (int i = 0; i < bucketSize; i++)
	{
		int fPIIndex = secondPermutedIndices[i];
		if (fPIIndex < firstPermutedIndices.size()) // ignore dummy tuples
			for (auto index : firstPermutedIndices[fPIIndex])
				invTotalPermutedIndices[index] = i;
	}

	uint16_t* Z1 = new uint16_t[myRowNum];
	uint16_t* Z2 = new uint16_t[myRowNum];
	OEP(invTotalPermutedIndices, bucketSize, myRowNum, (uint16_t*)&childAnnots2[0], Z1, Z2);
	for (int i = 0; i < myRowNum; i++)
	{
		Z1[i] += childAnnots1[invTotalPermutedIndices[i]];
		if (child.annotsBool)
			oblivExtTransfer(0, 0, annots1[i], annots2[i], Z1[i] & 1, Z2[i] & 1, annots1[i], annots2[i]);
		else if (annotsBool)
			oblivExtTransfer(0, 0, Z1[i], Z2[i], annots1[i] & 1, annots2[i] & 1, annots1[i], annots2[i]);
		else
			shareMul(annots1[i], annots2[i], Z1[i], Z2[i], annots1[i], annots2[i]);
	}
	if (!child.annotsBool)
		annotsBool = false;
	delete[] invTotalPermutedIndices;
	delete[] Z1;
	delete[] Z2;
	annotsShared = true;
}

void Relation::oblivSemiJoin(Relation& child)
{
	unordered_map<string, string> joinAttrNameMap;
	for (int i = 0; i < child.numAttrs; i++)
		joinAttrNameMap[child.attrNames[i]] = child.attrNames[i];
	oblivSemiJoin(child, joinAttrNameMap);
}

void Relation::oblivSemiJoin2(Relation& child, const unordered_map<string, string> &joinAttrNameMap)
{
	int numJoinAttrs = joinAttrNameMap.size();
	string* joinAttrNamesInChild = new string[numJoinAttrs];
	int i = 0;
	for (auto it : joinAttrNameMap)
		joinAttrNamesInChild[i++] = it.first;
	Relation childAnnotProjectRelation(child);
	childAnnotProjectRelation.annotProject(joinAttrNamesInChild, numJoinAttrs);
	oblivSemiJoin(childAnnotProjectRelation, joinAttrNameMap);
}

void Relation::oblivSemiJoin2(Relation& child)
{
	unordered_map<string, string> joinAttrNameMap;
	for (int i = 0; i < child.numAttrs; i++)
		joinAttrNameMap[child.attrNames[i]] = child.attrNames[i];
	oblivSemiJoin2(child, joinAttrNameMap);
}

void Relation::revealNonZeroAnnotTuples()
{
	int numRows = (int)tuples.size();
	vector<uint64_t*> nonZeroAnnotTuples;
	vector<uint16_t> newAnnot1, newAnnot2;
	for (int i = 0; i < numRows; i++)
	{
		bool z1 = annots1[i] & 1, z2 = annots2[i] & 1;
		if (!annotsBool)
		{
			equalityTest(annots1[i], -annots2[i], 16, z1, z2);
			z1 = !z1;
		}
		addComm(1); // Reveals equality
		if (z1 ^ z2) // non-zero annot tuple
		{
			nonZeroAnnotTuples.push_back(tuples[i]);
			newAnnot1.push_back(annots1[i]);
			newAnnot2.push_back(annots2[i]);
			addComm(64 * numAttrs); // sends the tuple
		}
	}
	tuples = nonZeroAnnotTuples;
	annots1 = newAnnot1;
	annots2 = newAnnot2;
	tuplesPublic = true;
}

Relation Relation::join(Relation &lr, Relation& rr, const unordered_map<string, string> &joinAttrNameMap, vector<int>& AliceIndices, vector<int>& BobIndices)
{
	if (!rr.tuplesPublic)
		throw "Error: Join with non-public relation!";
	// We still need to store the corresponding indices map during the join, so as to reveal the annot at last

	// Step 1: define attribute names & types
	string* joinAttrNames = new string[lr.numAttrs + rr.numAttrs];
	int* attrIndex = new int[lr.numAttrs + rr.numAttrs];
	int joinNumAttrs = 0;
	for (int i = 0; i < lr.numAttrs; ++i) {
		attrIndex[joinNumAttrs] = i;
		joinAttrNames[joinNumAttrs++] = lr.attrNames[i];
	}
	for (int i = 0; i < rr.numAttrs; ++i) 
		if (joinAttrNameMap.find(rr.attrNames[i]) == joinAttrNameMap.end()) 
		{
			attrIndex[joinNumAttrs] = i;
			joinAttrNames[joinNumAttrs++] = rr.attrNames[i];
		}
	Relation joinRelation(joinAttrNames, joinNumAttrs, false);
	
	int joinAttrsId = 0;
	for (int i = 0; i < lr.numAttrs; ++i)
		joinRelation.attrTypes[joinAttrsId++] = lr.attrTypes[i];
	for (int i = 0; i < rr.numAttrs; ++i) {
		if (joinAttrNameMap.find(rr.attrNames[i]) == joinAttrNameMap.end())
			joinRelation.attrTypes[joinAttrsId++] = rr.attrTypes[i];
	}

	// Step 2: hash join two table
	Relation AliceProjectRelation(lr);
	Relation BobProjectRelation(rr);
	int projectSize = joinAttrNameMap.size(), attrNamesId = 0;
	string* AliceJoinAttrNames = new string[projectSize];
	string* BobJoinAttrNames = new string[projectSize];
	for (auto mapPair : joinAttrNameMap) {
		AliceJoinAttrNames[attrNamesId] = mapPair.second;
		BobJoinAttrNames[attrNamesId] = mapPair.first;
		attrNamesId++;
	}
	AliceProjectRelation.project(AliceJoinAttrNames, projectSize);
	BobProjectRelation.project(BobJoinAttrNames, projectSize);
	delete[] AliceJoinAttrNames;
	delete[] BobJoinAttrNames;

	unordered_map<string, int> AliceInvAttrMap, BobInvAttrMap;
	for (int i = 0; i < AliceProjectRelation.numAttrs; ++i) {
		AliceInvAttrMap[AliceProjectRelation.attrNames[i]] = i;
	}
	for (int i = 0; i < BobProjectRelation.numAttrs; ++i) {
		BobInvAttrMap[BobProjectRelation.attrNames[i]] = i;
	}
	int* AliceIndex = new int[projectSize];
	int* BobIndex = new int[projectSize];
	int indexid = 0;
	for (auto projectPair : joinAttrNameMap) {
		AliceIndex[indexid] = AliceInvAttrMap[projectPair.second];
		BobIndex[indexid] = BobInvAttrMap[projectPair.first];
		indexid++;
	}

	unordered_map<uint64_t, vector<int>> AliceRowIndexMap, BobRowIndexMap;
	for (int i = 0; i < AliceProjectRelation.tuples.size(); ++i)
		AliceRowIndexMap[hashMultiColumns(AliceProjectRelation.tuples[i], projectSize)].push_back(i);
	for (int i = 0; i < BobProjectRelation.tuples.size(); ++i)
		BobRowIndexMap[hashMultiColumns(BobProjectRelation.tuples[i], projectSize)].push_back(i);

	for (auto mapPair : AliceRowIndexMap) {
		if (BobRowIndexMap.find(mapPair.first) == BobRowIndexMap.end())
			continue;
		vector<int> AliceIds = mapPair.second;
		vector<int> BobIds = BobRowIndexMap[mapPair.first];
		for (int i = 0; i < AliceIds.size(); ++i) {
			for (int j = 0; j < BobIds.size(); ++j) {
				int ai = AliceIds[i], bi = BobIds[j];
				// check equality
				bool equality = true;
				for (int equid = 0; equid < projectSize; ++equid) {
					if (AliceProjectRelation.tuples[ai][AliceIndex[equid]] != BobProjectRelation.tuples[bi][BobIndex[equid]]) {
						equality = false;
						break;
					}
				}
				if (equality) {
					uint64_t* tuple = new uint64_t[joinRelation.numAttrs];
					for (int i = 0; i < lr.numAttrs; i++)
						tuple[i] = lr.tuples[ai][attrIndex[i]];
					for (int i = lr.numAttrs; i < joinRelation.numAttrs; i++)
						tuple[i] = rr.tuples[bi][attrIndex[i]];
					joinRelation.tuples.push_back(tuple);
					AliceIndices.push_back(ai);
					BobIndices.push_back(bi);
				}
			}
		}
	}
	joinRelation.annots1.resize(joinRelation.tuples.size());
	joinRelation.annots2.resize(joinRelation.tuples.size());
	joinRelation.tuplesPublic = true;
	delete[] joinAttrNames;
	delete[] attrIndex;
	return joinRelation;
}

// This works only for two table join. For more tables it is similar (TODO)
Relation Relation::oblivJoin(Relation &lr, Relation& rr, const unordered_map<string, string> &joinAttrNameMap) {
	// Reveal phase
	// This only works for the case revealing the results to both
	// For the case only one party should know the result, we just need an OT (TODO)
	Relation R1(lr);
	Relation R2(rr);
	R1.revealNonZeroAnnotTuples();
	R2.revealNonZeroAnnotTuples();

	// Join phase
	vector<int> AliceIndices, BobIndices;
	int originAliceSize = R1.tuples.size();
	int originBobSize = R2.tuples.size();
	Relation joinRelation = join(R1, R2, joinAttrNameMap, AliceIndices, BobIndices);
	joinRelation.annotsShared = true;
	int joinRowNum = joinRelation.tuples.size();

	// Alice send the size of |J*| to Bob
	addComm(32);

	// compute annotations phase
	uint16_t* z1 = new uint16_t[joinRowNum];
	uint16_t* z2 = new uint16_t[joinRowNum];
	OEP((int*)&AliceIndices[0], originAliceSize, joinRowNum, (uint16_t*)&R1.annots2[0], z1, z2);
	for (int i = 0; i < joinRowNum; ++i)
		z1[i] += R1.annots1[AliceIndices[i]];

	uint16_t* z3 = new uint16_t[joinRowNum];
	uint16_t* z4 = new uint16_t[joinRowNum];
	OEP((int*)&BobIndices[0], originBobSize, joinRowNum, (uint16_t*)&R2.annots2[0], z3, z4);
	for (int i = 0; i < joinRowNum; ++i) {
		z3[i] += R2.annots1[BobIndices[i]];
		if (R2.annotsBool)
			oblivExtTransfer(0, 0, z1[i], z2[i], z3[i] & 1, z4[i] & 1, joinRelation.annots1[i], joinRelation.annots2[i]);
		else if (R1.annotsBool)
			oblivExtTransfer(0, 0, z3[i], z4[i], z1[i] & 1, z2[i] & 1, joinRelation.annots1[i], joinRelation.annots2[i]);
		else
			shareMul(z1[i], z2[i], z3[i], z4[i], joinRelation.annots1[i], joinRelation.annots2[i]);
	}
	joinRelation.annotsBool = R1.annotsBool && R2.annotsBool;

	delete[] z1;
	delete[] z2;
	delete[] z3;
	delete[] z4;
	return joinRelation;
}

Relation::~Relation()
{
	delete[] attrNames;
	delete[] attrTypes;
	for (auto row : tuples)
		delete[] row;
}


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

using namespace std;

Relation::Relation(Relation& R)
{
	numAttrs = R.numAttrs;
	annot1 = R.annot1;
	annot2 = R.annot2;
	annotsShared = R.annotsShared;
	tuplesPublic = R.tuplesPublic;
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

Relation::Relation(string* attrNames, int numAttrs)
{
	this->attrNames = new string[numAttrs];
	copy(attrNames, attrNames + numAttrs, this->attrNames);
	this->attrTypes = new DataType[numAttrs]();
	this->numAttrs = numAttrs;
	this->annotsShared = false;
	this->tuplesPublic = false;
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
				annot1.push_back(stoi(element));

		}
		tuples.push_back(row);
		fin >> ws;
	}
	fin.close();
	annot2.resize(annot1.size());
	delete[] indexMap;
}

void Relation::reveal()
{
	for (int i = 0; i < annot1.size(); i++)
	{
		addComm(16);
		annot1[i] += annot2[i];
		annot2[i] = 0;
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
		uint16_t realAnnot = annot1[i] + annot2[i];
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
		uint16_t out = showRealAnnot ? realAnnot : annot1[i];
		std::cout << out << '\n';
	}
	cout << '\n';
}

// Note: this project operation does not elimiate duplicate tuples!
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

// This corresponds to the pi_1 operator, which eliminates duplicate tuples
// It sets the annotation of a tuple as 1 if at least one of its duplicates has non-zero annotation
void Relation::oblivProject(string* projectAttrNames, int numProjectAttrs)
{
	int numRows = (int)tuples.size();
	bool* b1 = new bool[numRows];
	bool* b2 = new bool[numRows];
	for (int i = 0; i < numRows; i++)
		equalityTest(annot1[i], -annot2[i], 16, b1[i], b2[i]);
	throw "Oblivious projection not implemented yet!";
	// Design the AND switch gate like oblivAgg
}

void Relation::permuteTuples(int* permutedIndices, bool permuteAnnot)
{
	auto origTuples = tuples;
	for (int i = 0; i < tuples.size(); i++)
		tuples[i] = origTuples[permutedIndices[i]];
	if (permuteAnnot)
	{
		auto origAnnot1 = annot1;
		for (int i = 0; i < tuples.size(); i++)
			annot1[i] = origAnnot1[permutedIndices[i]];
	}
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
	int* permutedIndices = new int[numRows];
	for (int i = 0; i < numRows; i++)
		permutedIndices[i] = i;
	sort(permutedIndices, permutedIndices + numRows, [&](int i, int j) {
		for (int n = 0; n < numAttrs; n++)
		{
			if (tuples[i][n] < tuples[j][n])
				return true;
			else if (tuples[i][n] > tuples[j][n])
				return false;
		}
		return false;
		});
	permuteTuples(permutedIndices);
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
			annot1[i + 1] += annot1[i];
			annot1[i] = 0;
			annot2[i + 1] += annot2[i];
			annot2[i] = 0;
			for (int j = 0; j < numAttrs; j++) // replace with dummy tuple
				tuples[i][j] = gRNG.nextUInt64();
		}
	}
}

void Relation::oblivAgg(string* groupByAttrNames, int numGroupByAttrs)
{
	size_t numRows = tuples.size();
	project(groupByAttrNames, numGroupByAttrs);
	int* permutedIndices = new int[numRows];
	for (int i = 0; i < numRows; i++)
		permutedIndices[i] = i;
	sort(permutedIndices, permutedIndices + numRows, [&](int i, int j) {
		for (int n = 0; n < numAttrs; n++)
		{
			if (tuples[i][n] < tuples[j][n])
				return true;
			else if (tuples[i][n] > tuples[j][n])
				return false;
		}
		return false;
		});
	permuteTuples(permutedIndices);
	uint16_t* Z1 = new uint16_t[numRows];
	uint16_t* Z2 = new uint16_t[numRows];
	oblivPermutation(permutedIndices, (uint16_t*)&annot2[0], (int)numRows, Z1, Z2);
	for (int i = 0; i < numRows; i++)
		Z1[i] += annot1[i];

	uint16_t msg0[2];
	uint16_t msg1[2];
	uint16_t msgb[2];
	// The following code works only for little-endian
	for (int i = 0; i < numRows - 1; i++)
	{
		// a,b,c,d and messages belong to Bob
		uint16_t a = Z2[i], b = Z2[i + 1];
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
			Z1[i + 1] += Z1[i] + msgb[1];
			Z1[i] = msgb[0];
			for (int j = 0; j < numAttrs; j++) // replace with dummy tuple
				tuples[i][j] = gRNG.nextUInt64();
		}
		else
		{
			Z1[i] += msgb[0];
			Z1[i + 1] += msgb[1];
		}
		Z2[i] = c;
		Z2[i + 1] = d;
	}
	copy(Z1, Z1 + numRows, annot1.begin());
	copy(Z2, Z2 + numRows, annot2.begin());
	delete[] permutedIndices;
	delete[] Z1;
	delete[] Z2;
	annotsShared = true;
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

void Relation::oblivSemiJoin(Relation& child, unordered_map<string, string> joinAttrNameMap)
{
	Relation myProjectRelation(*this);
	string* joinAttrNames = new string[child.numAttrs];
	for (int i = 0; i < child.numAttrs; i++)
		joinAttrNames[i] = joinAttrNameMap[child.attrNames[i]];
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
	vector<uint16_t> childAnnot1;

	// Bob
	vector<uint64_t> childHashValues;
	childHashValues.reserve(child.tuples.size());
	for (auto row : child.tuples)
		childHashValues.push_back(hashMultiColumns(row, numJoinAttrs));
	vector<uint16_t> childAnnot2;

	// PSI
	if (child.annotsShared)
		PSIwithSharedPayload(myHashValues, childHashValues, child.annot1, child.annot2, secondPermutedIndices, childAnnot1, childAnnot2);
	else
		PSIwithPayload(myHashValues, childHashValues, child.annot1, secondPermutedIndices, childAnnot1, childAnnot2);
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
	OEP(invTotalPermutedIndices, bucketSize, myRowNum, (uint16_t*)&childAnnot2[0], Z1, Z2);
	for (int i = 0; i < myRowNum; i++)
	{
		Z1[i] += childAnnot1[invTotalPermutedIndices[i]];
		shareMul(annot1[i], annot2[i], Z1[i], Z2[i], annot1[i], annot2[i]);
	}
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

void Relation::revealNonZeroAnnotTuples()
{
	int numRows = (int)tuples.size();
	vector<uint64_t*> nonZeroAnnotTuples;
	vector<uint16_t> newAnnot1, newAnnot2;
	for (int i = 0; numRows; i++)
	{
		bool z1 = 0, z2 = 0;
		equalityTest(annot1[i], -annot2[i], 16, z1, z2);
		addComm(1); // Reveals equality
		z1 ^= z2; 
		if (!z1) // non-zero annot tuple
		{
			nonZeroAnnotTuples.push_back(tuples[i]);
			newAnnot1.push_back(annot1[i]);
			newAnnot2.push_back(annot2[i]);
			addComm(64 * numAttrs); // sends the tuple
		}
	}
	tuples = nonZeroAnnotTuples;
	annot1 = newAnnot1;
	annot2 = newAnnot2;
	tuplesPublic = true;
}

void Relation::join(Relation& child)
{
	if (!child.tuplesPublic)
		throw "Error: Join with non-public relation!";
	// We still need to store the corresponding indices map during the join, so as to reveal the annot at last
	throw "Join operator not implemented yet!";

}

Relation::~Relation()
{
	delete[] attrNames;
	delete[] attrTypes;
	for (auto row : tuples)
		delete[] row;
}
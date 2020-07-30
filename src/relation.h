#pragma once
#include <vector>
#include <string>
#include <unordered_map>

// File Structure:
// The first line: an integer indicating number of columns in the file
// The second line: attribute names (and annotation names, e.g. q3_annotation, q8_annotation1)
// The third line: attribute types (options: int, decimal, date, string, char)
// The attribute types for annotations can be arbitrary, which do not work
// Other lines: data


class Relation {
	enum class DataType
	{
		INT,
		DECIMAL,
		DATE,
		STRING,
		CHAR
	};
public:
	Relation(const Relation& R);
	// Annotation name must not be included in attrNames
	Relation(std::string* attrNames, int numAttrs, bool annotsBool);
	void loadData(const char* filePath, std::string annotationName);
	void print(bool showRealAnnot = true, bool showZeroAnnotTuples = false, int limit_size = 100);
	void reveal(); // sends shares to reveal annotations
	void sort();
	void agg(std::string* groupByAttrNames, int numGroupByAttrs);
	void oblivAgg(std::string* groupByAttrNames, int numGroupByAttrs);

	// This corresponds to the pi_1 operator, which eliminates duplicate tuples (to zero-annotated dummy tuples)
	// It sets annotation of a tuple as 1 if at least one of its duplicates has non-zero annotation
	void annotProject(std::string* projectAttrNames, int numProjectAttrs);
	// The shared version of annotProject
	void oblivAnnotProject(std::string* projectAttrNames, int numProjectAttrs);

	void oblivSemiJoin(Relation& child);
	// Note: the keys of the map is attrNames of child relation, while the values are corresponding join attrNames of this relation
	void oblivSemiJoin(Relation& child, const std::unordered_map<std::string, std::string>& joinAttrNameMap);

	// The second type of semijoin that is used to remove dangling tuples 
	void oblivSemiJoin2(Relation& child);
	void oblivSemiJoin2(Relation& child, const std::unordered_map<std::string, std::string>& joinAttrNameMap);

	void revealNonZeroAnnotTuples();
	static Relation oblivJoin(Relation& lr, Relation& rr, const std::unordered_map<std::string, std::string>& joinAttrNameMap);
	~Relation();
private:
	std::string* attrNames;
	DataType* attrTypes;
	int numAttrs;
	std::vector< uint64_t* > tuples;
	std::vector< uint16_t > annots1;
	std::vector< uint16_t > annots2;
	bool annotsShared; // If this is true, the owner of this relation does not know the annots
	bool tuplesPublic; // If this is true, both parties know the tuples (but not the annots)
	bool annotsBool; // If this is true, then only the last bits of the annotations make sense, which indicate whether tuples are non-dummy or not

	// Note: this project operation does not elimiate duplicate tuples!
	void project(std::string* projectAttrNames, int numGroupByAttrs);
	void permute(int* permutedIndices);
	
	static Relation join(Relation& lr, Relation& rr, const std::unordered_map<std::string, std::string>& joinAttNameMap, std::vector<int>& AliceIndices, std::vector<int>& BobIndices);

};


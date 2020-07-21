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
    Relation(Relation& R);
    // Annotation name must not be included in attrNames
    Relation(std::string* attrNames, int numAttrs);
    void loadData(const char* filePath, std::string annotationName);
    void print(bool showRealAnnot = true, bool showZeroAnnotTuples = false, int limit_size = 100);
    void reveal();
    void agg(std::string* groupByAttrNames, int numGroupByAttrs);
    void oblivAgg(std::string* groupByAttrNames, int numGroupByAttrs);
    void oblivProject(std::string* projectAttrNames, int numProjectAttrs);
    void oblivSemiJoin(Relation& child );
    // Note: the keys of the map is attrNames of child relation, while the values are corresponding join attrNames of this relation
    void oblivSemiJoin(Relation &child, std::unordered_map<std::string,std::string> joinAttrNameMap);
    void revealNonZeroAnnotTuples();
    void join(Relation& child);
    ~Relation();
private:
    std::string* attrNames;
    DataType* attrTypes;
    int numAttrs;
    std::vector< uint64_t* > tuples;
    std::vector< uint16_t > annot1;
    std::vector< uint16_t > annot2;
    bool annotsShared; // If this is true, the owner of this relation does not know the annots
    bool tuplesPublic; // If this is true, both parties know the tuples (but not the annots)
    void project(std::string* projectAttrNames, int numGroupByAttrs);
    void permuteTuples(int* permutedIndices, bool permuteAnnot = true);
};


#include "src/relation.h"
#include "src/OEP.h"
#include "src/utils.h"
#include <string>
#include <ctime>
#include <iostream>

using namespace std;

#ifdef _WIN32
#include <direct.h>
#define chdir _chdir
#else
#include "unistd.h"
#endif

void q3(char* datapath)
{
    chdir(datapath);
    string annot = "q3_annot";
    string customerAttrs[] = { "c_custkey" };
    
    Relation customer(customerAttrs, 1);
    customer.loadData("customer.tbl", annot);

    string ordersAttrs[] = { "o_custkey", "o_orderkey", "o_orderdate", "o_shippriority" };
    Relation orders(ordersAttrs, 4);
    orders.loadData("orders.tbl", annot);

    string lineitemAttrs[] = { "l_orderkey" };
    Relation lineitem(lineitemAttrs, 1);
    lineitem.loadData("lineitem.tbl", annot);
    lineitem.agg(lineitemAttrs, 1);
    
    clock_t start, end;
    start = clock();
    unordered_map<string, string> C2O={{"c_custkey","o_custkey"}};
    orders.oblivSemiJoin(customer, C2O);
    unordered_map<string, string> L2O = { {"l_orderkey","o_orderkey"} };
    orders.oblivSemiJoin(lineitem, L2O);
    string groupByAttrs[] = { "o_orderkey", "o_orderdate", "o_shippriority" };
    orders.oblivAgg(groupByAttrs, 3);
    orders.reveal();
    end = clock();

    char comm_text[20];
    getComm(comm_text);
    cout << "Q3: " << 1.0 * (end - start) / CLOCKS_PER_SEC << " s, " << comm_text << endl;
    orders.print();
    resetComm();
}

int main()
{
    char datapath[] = "data";
    q3(datapath);
    return 0;
}
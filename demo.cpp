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

void q3()
{
    string annot = "q3_annot";
    string customerAttrs[] = { "c_custkey" };
    
    Relation customer(customerAttrs, sizeof(customerAttrs)/sizeof(string), true);
    customer.loadData("customer.tbl", annot);

    string ordersAttrs[] = { "o_custkey", "o_orderkey", "o_orderdate", "o_shippriority" };
    Relation orders(ordersAttrs, sizeof(ordersAttrs)/sizeof(string), true);
    orders.loadData("orders.tbl", annot);

    string lineitemAttrs[] = { "l_orderkey" };
    Relation lineitem(lineitemAttrs, sizeof(lineitemAttrs)/sizeof(string), false);
    lineitem.loadData("lineitem.tbl", annot);
    lineitem.agg(lineitemAttrs, sizeof(lineitemAttrs) / sizeof(string));

    clock_t start, end;
    start = clock();
    unordered_map<string, string> C2O={{"c_custkey","o_custkey"}};
    orders.oblivSemiJoin(customer, C2O);
    unordered_map<string, string> L2O = { {"l_orderkey","o_orderkey"} };
    orders.oblivSemiJoin(lineitem, L2O);
    string groupByAttrs[] = { "o_orderkey", "o_orderdate", "o_shippriority" };
    orders.oblivAgg(groupByAttrs, sizeof(groupByAttrs)/sizeof(string));
    orders.reveal();
    end = clock();

    char comm_text[20];
    getComm(comm_text);
    cout << "Q3: " << 1.0 * (end - start) / CLOCKS_PER_SEC << " s, " << comm_text << endl;
    orders.print();
    resetComm();
}

void q10()
{
    string annot = "q10_annot";

    string customerAttrs[] = { "c_custkey", "c_name", "c_nationkey" };
    Relation customer(customerAttrs, sizeof(customerAttrs)/sizeof(string), false);
    customer.loadData("customer.tbl", annot);

    string ordersAttrs[] = { "o_custkey", "o_orderkey" };
    Relation orders(ordersAttrs, sizeof(ordersAttrs)/sizeof(string), false);
    orders.loadData("orders.tbl", annot);

    string lineitemAttrs[] = { "l_orderkey" };
    Relation lineitem(lineitemAttrs, sizeof(lineitemAttrs) / sizeof(string), false);
    lineitem.loadData("lineitem.tbl", annot);
    lineitem.agg(lineitemAttrs, sizeof(lineitemAttrs) / sizeof(string));

    clock_t start, end;
    start = clock();
    unordered_map<string, string> L2O = { {"l_orderkey","o_orderkey"} };
    orders.oblivSemiJoin(lineitem, L2O);

    string ordersAgg[] = { "o_custkey" };
    orders.oblivAgg(ordersAgg, sizeof(ordersAgg) / sizeof(string));

    unordered_map<string, string> O2C = { {"o_custkey", "c_custkey"} };
    customer.oblivSemiJoin(orders, O2C);
    customer.reveal();
    end = clock();

    char comm_text[20];
    getComm(comm_text);
    cout << "Q10: " << 1.0 * (end - start) / CLOCKS_PER_SEC << " s, " << comm_text << endl;
    customer.print();
    resetComm();
}


void q8()
{
    string annot = "q8_annot";

    string customerAttrs[] = { "c_custkey" };
    Relation customer(customerAttrs, sizeof(customerAttrs) / sizeof(string), true);
    customer.loadData("customer.tbl", annot);

    string ordersAttrs[] = { "o_orderkey", "o_custkey", "o_year" };
    Relation orders1(ordersAttrs, sizeof(ordersAttrs) / sizeof(string), true);
    orders1.loadData("orders.tbl", annot);

    string lineitemAttrs[] = { "l_orderkey", "l_suppkey", "l_partkey" };
    Relation lineitem1(lineitemAttrs, sizeof(lineitemAttrs) / sizeof(string), false);
    lineitem1.loadData("lineitem.tbl", annot);
    lineitem1.agg(lineitemAttrs, sizeof(lineitemAttrs) / sizeof(string)); 

    string partAttrs[] = { "p_partkey"};
    Relation part(partAttrs, sizeof(partAttrs) / sizeof(string), true);
    part.loadData("part.tbl", annot);

    string supplierAttrs[] = { "s_suppkey" };
    Relation supplier1(supplierAttrs, sizeof(supplierAttrs) / sizeof(string), true);
    Relation supplier2(supplier1);
    supplier1.loadData("supplier.tbl", "q8_annot1");
    supplier2.loadData("supplier.tbl", "q8_annot2");

    clock_t start, end;
    start = clock();

    unordered_map<string, string> C2O = { {"c_custkey", "o_custkey"} };
    orders1.oblivSemiJoin(customer, C2O);

    unordered_map<string, string> P2L = { {"p_partkey","l_partkey"} };
    lineitem1.oblivSemiJoin(part, P2L);

    unordered_map<string, string> S2L = { {"s_suppkey","l_suppkey"} };
    string lineitemAgg[] = { "l_orderkey" };
    Relation lineitem2(lineitem1);
    lineitem1.oblivSemiJoin(supplier1, S2L);
    lineitem1.oblivAgg(lineitemAgg, sizeof(lineitemAgg) / sizeof(string));
    lineitem2.oblivSemiJoin(supplier2, S2L);
    lineitem2.oblivAgg(lineitemAgg, sizeof(lineitemAgg) / sizeof(string));

    string ordersAgg[] = { "o_year" };
    unordered_map<string, string> L2O = { {"l_orderkey", "o_orderkey"} };
    Relation orders2(orders1);
    orders1.oblivSemiJoin(lineitem1, L2O);
    orders1.oblivAgg(ordersAgg, sizeof(ordersAgg) / sizeof(string));
    orders1.reveal();
    orders2.oblivSemiJoin(lineitem2, L2O);
    orders2.oblivAgg(ordersAgg, sizeof(ordersAgg) / sizeof(string));
    orders2.reveal(); // Actually we should perform oblivious division before revealing, but we leave this to future work
    end = clock();

    char comm_text[20];
    getComm(comm_text);
    cout << "Q8: " << 1.0 * (end - start) / CLOCKS_PER_SEC << " s, " << comm_text << endl;
    orders1.print();
    orders2.print();
    resetComm();
}

void q18()
{
    string annot = "q18_annot";

    string customerAttrs[] = { "c_custkey", "c_name" };
    Relation customer(customerAttrs, sizeof(customerAttrs)/sizeof(string), false);
    customer.loadData("customer.tbl", annot);

    string ordersAttrs[] = { "o_custkey", "o_orderkey", "o_orderdate", "o_totalprice" };
    Relation orders(ordersAttrs, sizeof(ordersAttrs) / sizeof(string), false);
    orders.loadData("orders.tbl", annot);

    string lineitemAttrs[] = { "l_orderkey" };
    Relation lineitem(lineitemAttrs, sizeof(lineitemAttrs) / sizeof(string), false);
    lineitem.loadData("lineitem.tbl", annot);
    lineitem.agg(lineitemAttrs, sizeof(lineitemAttrs) / sizeof(string));

    clock_t start, end;
    start = clock();
    unordered_map<string, string> L2O = { {"l_orderkey","o_orderkey"} };
    orders.oblivSemiJoin(lineitem, L2O);

    unordered_map<string, string> C2O = { {"c_custkey","o_custkey"} };
    orders.oblivSemiJoin2(customer, C2O);

    unordered_map<string, string> O2C = { {"o_custkey", "c_custkey"} };
    customer.oblivSemiJoin2(orders, O2C);

    Relation result = Relation::oblivJoin(customer, orders, O2C);
    result.reveal();
    result.sort();
    end = clock();

    char comm_text[20];
    getComm(comm_text);
    cout << "Q18: " << 1.0 * (end - start) / CLOCKS_PER_SEC << " s, " << comm_text << endl;
    result.print(false, true); // check that there is no zero-annotated tuples in result
    resetComm();
}

int main()
{
    char datapath[] = "data";
    chdir(datapath);
    while (1)
    {
        cout << "Input query id [3/8/10/18]: ";
        int qid;
        cin >> qid;
        switch (qid)
        {
        case 3: q3(); break;
        case 8: q8(); break;
        case 10: q10(); break;
        case 18: q18(); break;
        default: cout << "Input Error!" << endl;
        }
    }
    return 0;
}
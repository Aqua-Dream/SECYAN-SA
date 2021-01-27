# Secure Yannakakis Stand-Alone Demo
A stand-alone implementation for Secure Yannakakis protocol. Visit [SECYAN](https://github.com/hkustDB/SECYAN/) for standard version.

## Requirements

1. [NTL](https://www.shoup.net/ntl/)
2. [Blake3](https://github.com/BLAKE3-team/BLAKE3)

## Instructions
To run the demo with cmake, you can run the following instructions through terminal.

```
mkdir build
cd build
cmake ../
make
cp demo ../
cd ..
./demo
```

Please note that the exetuable file `demo` must be at the top level directory.

## Data
We have provided 1MB data in the `data` folder. It also contains a 10MB data zip file. If you want to use your own data, please make sure the formats of your files are consistent with that we provided.

## Communication Cost

| Query   |  1 MB         | 10 MB         |
| -------- | -----------  | ------------  | 
| 3           | 7.10 MB     |79.53 MB    |
| 8           | 64.13 MB   |739.01 MB  |
| 10         | 6.01 MB     |67.92 MB    |
| 18         | 10.46 MB   |116.29 MB  |

## Example

For TPC-H Q3:
``` sql
select o_orderkey, o_orderdate, o_shippriority,
    sum(l_extendedprice * (1 - l_discount)) as revenue
from customer, orders, lineitem
where c_mktsegment = 'AUTOMOBILE'
    and c_custkey = o_custkey
    and l_orderkey = o_orderkey
    and o_orderdate < date '1995-03-13'
    and l_shipdate > date '1995-03-13'
group by o_orderkey, o_orderdate, o_shippriority;
```

``` bash
./demo
Input query id [3/8/10/18]: 3
Q3: 0.53 s, 7.10 MB
o_orderkey      o_orderdate     o_shippriority  annotation
1092    1995-03-04      0       8006
1830    1995-02-23      0       7164
2053    1995-02-07      0       12143
3110    1994-12-17      0       2937
3814    1995-02-22      0       12594
4134    1995-01-12      0       12116
4227    1995-02-24      0       8724
4550    1994-12-29      0       898
4707    1995-02-27      0       5718
4960    1995-02-26      0       11275
5312    1995-02-24      0       6176
```

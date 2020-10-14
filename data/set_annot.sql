-- avoid overflow as we are using 16-bit annotation
-- It requires more work to support annotations for any number of bits
UPDATE lineitem
SET l_extendedprice = l_extendedprice / 10;

UPDATE partsupp
SET ps_supplycost = ps_supplycost / 10;

ALTER TABLE orders
ADD o_year INT;

UPDATE orders 
SET o_year = extract(year from o_orderdate);

-- Q3
ALTER TABLE customer
ADD q3_annot INT DEFAULT 0;
ALTER TABLE lineitem
ADD q3_annot INT DEFAULT 0;
ALTER TABLE orders
ADD q3_annot INT DEFAULT 0;

UPDATE orders
SET q3_annot = 1 
WHERE o_orderdate < date '1995-03-13';

UPDATE lineitem
SET q3_annot = FLOOR(l_extendedprice * (1 - l_discount)+0.5) 
WHERE l_shipdate > date '1995-03-13';

UPDATE customer
SET q3_annot = 1
WHERE c_mktsegment = 'AUTOMOBILE' ;


ALTER TABLE customer
ADD q10_annot INT DEFAULT 1;
ALTER TABLE lineitem
ADD q10_annot INT DEFAULT 0;
ALTER TABLE orders
ADD q10_annot INT DEFAULT 0;

UPDATE lineitem
SET q10_annot =  FLOOR(l_extendedprice * (1 - l_discount)+0.5) 
WHERE l_returnflag = 'R';

UPDATE orders
SET q10_annot = 1 
WHERE o_orderdate >= date '1993-08-01' AND o_orderdate < date '1993-11-01';

-- Q18

ALTER TABLE customer
ADD q18_annot INT DEFAULT 1;
ALTER TABLE lineitem
ADD q18_annot INT DEFAULT 0;
ALTER TABLE orders
ADD q18_annot INT DEFAULT 1;

UPDATE lineitem
SET q18_annot=l_quantity
WHERE l_orderkey in (select l_orderkey from (
    (select l_orderkey
    from lineitem
    group by l_orderkey
    having sum(l_quantity) > 300) AS temp
) );

-- Q8
ALTER TABLE customer
ADD q8_annot INT DEFAULT 0;
ALTER TABLE lineitem
ADD q8_annot INT DEFAULT 0;
ALTER TABLE orders
ADD q8_annot INT DEFAULT 0;
ALTER TABLE part
ADD q8_annot INT DEFAULT 0;
ALTER TABLE supplier
ADD q8_annot1 INT DEFAULT 0;
ALTER TABLE supplier
ADD q8_annot2 INT DEFAULT 1;

UPDATE customer
SET q8_annot = 1
WHERE c_nationkey in (8,9,12,18,21);

UPDATE lineitem
SET q8_annot = FLOOR(l_extendedprice * (1 - l_discount)+0.5);

UPDATE part
SET q8_annot = 1
WHERE p_type = 'SMALL PLATED COPPER';

UPDATE orders
SET q8_annot = 1
WHERE o_orderdate between date '1995-01-01' and date '1996-12-31';

UPDATE supplier
SET q8_annot1 = 1
WHERE s_nationkey=8;

--Q9
ALTER TABLE lineitem
ADD q9_annot1 INT;
ALTER TABLE lineitem
ADD q9_annot2 INT;
ALTER TABLE orders
ADD q9_annot INT DEFAULT 1;
ALTER TABLE part
ADD q9_annot INT DEFAULT 0;
ALTER TABLE partsupp
ADD q9_annot1 INT DEFAULT 1;
ALTER TABLE partsupp
ADD q9_annot2 INT;

UPDATE part
SET q9_annot = 1
WHERE p_name like '%green%';

UPDATE lineitem
SET q9_annot1 = FLOOR(l_extendedprice * (1 - l_discount)+0.5);

UPDATE lineitem
SET q9_annot2 = FLOOR(l_quantity+0.5);

UPDATE partsupp
SET q9_annot2 = FLOOR(ps_supplycost+0.5);


ALTER TABLE supplier
ADD q9_annot0 INT DEFAULT 0;
UPDATE supplier
SET q9_annot0 = 1
WHERE s_nationkey=0;

ALTER TABLE supplier
ADD q9_annot1 INT DEFAULT 0;
UPDATE supplier
SET q9_annot1 = 1
WHERE s_nationkey=1;

ALTER TABLE supplier
ADD q9_annot2 INT DEFAULT 0;
UPDATE supplier
SET q9_annot2 = 1
WHERE s_nationkey=2;

ALTER TABLE supplier
ADD q9_annot3 INT DEFAULT 0;
UPDATE supplier
SET q9_annot3 = 1
WHERE s_nationkey=3;

ALTER TABLE supplier
ADD q9_annot4 INT DEFAULT 0;
UPDATE supplier
SET q9_annot4 = 1
WHERE s_nationkey=4;

ALTER TABLE supplier
ADD q9_annot5 INT DEFAULT 0;
UPDATE supplier
SET q9_annot5 = 1
WHERE s_nationkey=5;

ALTER TABLE supplier
ADD q9_annot6 INT DEFAULT 0;
UPDATE supplier
SET q9_annot6 = 1
WHERE s_nationkey=6;

ALTER TABLE supplier
ADD q9_annot7 INT DEFAULT 0;
UPDATE supplier
SET q9_annot7 = 1
WHERE s_nationkey=7;

ALTER TABLE supplier
ADD q9_annot8 INT DEFAULT 0;
UPDATE supplier
SET q9_annot8 = 1
WHERE s_nationkey=8;

ALTER TABLE supplier
ADD q9_annot9 INT DEFAULT 0;
UPDATE supplier
SET q9_annot9 = 1
WHERE s_nationkey=9;

ALTER TABLE supplier
ADD q9_annot10 INT DEFAULT 0;
UPDATE supplier
SET q9_annot10 = 1
WHERE s_nationkey=10;

ALTER TABLE supplier
ADD q9_annot11 INT DEFAULT 0;
UPDATE supplier
SET q9_annot11 = 1
WHERE s_nationkey=11;


ALTER TABLE supplier
ADD q9_annot12 INT DEFAULT 0;
UPDATE supplier
SET q9_annot12 = 1
WHERE s_nationkey=12;

ALTER TABLE supplier
ADD q9_annot13 INT DEFAULT 0;
UPDATE supplier
SET q9_annot13 = 1
WHERE s_nationkey=13;

ALTER TABLE supplier
ADD q9_annot14 INT DEFAULT 0;
UPDATE supplier
SET q9_annot14 = 1
WHERE s_nationkey=14;

ALTER TABLE supplier
ADD q9_annot15 INT DEFAULT 0;
UPDATE supplier
SET q9_annot15 = 1
WHERE s_nationkey=15;

ALTER TABLE supplier
ADD q9_annot16 INT DEFAULT 0;
UPDATE supplier
SET q9_annot16 = 1
WHERE s_nationkey=16;

ALTER TABLE supplier
ADD q9_annot17 INT DEFAULT 0;
UPDATE supplier
SET q9_annot17 = 1
WHERE s_nationkey=17;

ALTER TABLE supplier
ADD q9_annot18 INT DEFAULT 0;
UPDATE supplier
SET q9_annot18 = 1
WHERE s_nationkey=18;

ALTER TABLE supplier
ADD q9_annot19 INT DEFAULT 0;
UPDATE supplier
SET q9_annot19 = 1
WHERE s_nationkey=19;

ALTER TABLE supplier
ADD q9_annot20 INT DEFAULT 0;
UPDATE supplier
SET q9_annot20 = 1
WHERE s_nationkey=20;

ALTER TABLE supplier
ADD q9_annot21 INT DEFAULT 0;
UPDATE supplier
SET q9_annot21 = 1
WHERE s_nationkey=21;

ALTER TABLE supplier
ADD q9_annot22 INT DEFAULT 0;
UPDATE supplier
SET q9_annot22 = 1
WHERE s_nationkey=22;

ALTER TABLE supplier
ADD q9_annot23 INT DEFAULT 0;
UPDATE supplier
SET q9_annot23 = 1
WHERE s_nationkey=23;

ALTER TABLE supplier
ADD q9_annot24 INT DEFAULT 0;
UPDATE supplier
SET q9_annot24 = 1
WHERE s_nationkey=24;

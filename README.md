
# PTAB for Top-k Subgraph Matching

## Introduction
This algorithm is built on the [gStore system](https://github.com/pkumod/gStore), you can see it for more instruction on how to use this system. 
The implementations of algorithms can be found at `Database/Optimizer` and `Query/topk`.

## Build
After download the source codes, build the PTAB by:
```shell
make -j pre
make -j
```

## How to Use PTAB
### Prepare Dataset if Needed
If needed, you should clean your dataset and make sure every node only has one label (type). The codes in folder `PrepareDataset` helps for that.
You can build a database on the uncleaned dataset and output the needed records. Then rebuild the database on the cleaned datasets.

### Build the Database and HopIndex
Use `./bin/gbuild [db_name] [db_file_path]` to build the database. Then use `./bin/indexBuilder [db_name]` to build the hop index (the parameter $H$ is default to be $2$).

### Run Top-k Subgraph Matching Query
You can input the top-k subgraph matching query with the following format:
```sparql
select [returned variable 1] [returned variable 2] ...
where { pattern 1. pattern 2. ...}
order by [score_function] limit [k]
```
An example:
```sparql
select ?v0 ?v1 ?v2 ?N1 ?N2 where {
  ?upt <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://www4.wiwiss.fu-berlin.de/bizer/bsbm/v01/vocabulary/ProductType>.
  ?v0 <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> ?upt.
  ?v0 <http://purl.org/dc/elements/1.1/publisher> ?v1 .
  ?v0 <http://purl.org/dc/elements/1.1/date> ?v2 .
  ?v0 <http://www4.wiwiss.fu-berlin.de/bizer/bsbm/v01/vocabulary/productPropertyNumeric1> ?N1.
  ?v0 <http://www4.wiwiss.fu-berlin.de/bizer/bsbm/v01/vocabulary/productPropertyNumeric5> ?N2.
} order by (?N1 - ?N2) limit 5
```

And the shell command is:
```shell
./bin/gquery [db_name] [query_path] [result_path] [top-k_method]
```
The `top-k_method` specify the top-k execution method. The available methods are:
* PTAB, the method introduced in our paper.
* DP-B.
* kTPM.
* SAE.
* Take-all.
* Take2.
* Eager.
* Naive (sort after matching).
`PATB` is the default method.


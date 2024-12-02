# Maximal Similar-weight Biclique Enumeration for Large Bipartite Graphs
## Introduction

For the details, please refer to our ICDE'2025 paper
"Maximal Similar-weight Biclique Enumeration for Large Bipartite Graphs"

If you have any further questions, please feel free to contact us.

Please cite our paper, if you use our source code.

## Compile
Under the root directory of the project, execute the following commands to compile the source code.

```zsh
g++ -std=c++11 -O3 mswbe.cpp -o mswbe
g++ main.cpp BCE.cpp naive.cpp mswbe.cpp -o naive
```

## Input
The graph is edge-weighted.
Each graph starts with 'N1 N2 M' where N1 and N2 are the number of vertices in left and right and M is the number of edges.
An edge is formatted as 'VertexID VertexID Weight'. 
Note that we require that the vertex id is started from 0 and the range is [0, N1 + N2 - 1]. 
The following is an input sample. 

Example:

```zsh
5 5 23
4 9 2
3 7 5
4 6 3
0 5 2
1 6 3
0 8 4
2 5 3
1 9 1
2 8 3
4 5 4
3 9 1
3 6 3
0 7 3
2 7 1
1 5 2
1 8 5
4 7 2
3 5 4
3 8 2
0 9 2
2 9 3
1 7 4
2 6 2
```

We can execute a enumeration by specifying the following parameter.

$code_dir $data_dir$ ${dvalues[d]} ${avalues[a]} ${bvalues[b]} ${reductions[r]} ${orders[o]} ${methods[m]} 

```zsh

./mswbe /test.txt 1 4 4 1 1
./naive /test.txt 1 4 4 1 1
```

## Experiment Datasets
The real world datasets used in our paper can be downloaded [here](http://konect.cc/).



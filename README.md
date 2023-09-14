# gcs_hadamard
Construction of Hadamard matrices from GCS, and verifications of their existence.
## Constrution of Hadamard matrices
Run test_gen_hadamard_from_pairs.m\
n1 = [l1, l2, ...lk], n2 = [m1, m2, ..., mk].

## Verifications of existence of 4-phase GCS quads and GCS octets
Set two parameters in gcs.cpp: \
long long n = (long long) 1e13; // maximum number, depends on the memory capacity.\
int tab_max_exp = 36; // a tradeoff between time complexity and space complexity, here it consume 2^36 bytes = 64GB.\
Then compile the cpp files and run the executable file:\
1. cd gcs_verifications\
2. clang++ gcs.cpp gen_S4.cpp -std=c++20 -fopenmp -O3\
3. ./a.out\
The density of gcs quad is saved in density.txt.\

## ## Verifications of existence of 4-phase GCS quads and GCS octets
Set two parameters in gcs.cpp: \
int len_exp = 44; // maximum number 2^44, depends on the memory capacity, here it consumes about 2.6TB.\
int tab_max_exp = 36; // a tradeoff between time complexity and space complexity, here it consume 2^36 bytes = 64GB.\
Then compile the cpp files and run the executable file:\
1. cd hadamard_verifications\
2. clang++ gcs_hadamard.cpp gen_S4.cpp -std=c++20 -fopenmp -O3.\
3. ./a.out\


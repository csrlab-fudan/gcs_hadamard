# gcs_hadamard
Construction of Hadamard matrices from GCS pairs, and verifications of their existence.
## Constrution of Hadamard matrices
Run test_gen_hadamard_from_pairs.m\
n1 = [l1, l2, ...lk], n2 = [m1, m2, ..., mk].

## Verifications of the existence of 4-phase GCS quads and GCS octets
Set two parameters in gcs.cpp:\
long long n = (long long) 1e13; // maximum number, restricted by the memory capacity, here it consumes about 3TB.\
int tab_max_exp = 36; // a tradeoff between the time complexity and the space complexity, here it consume $2^{36}$ bytes = 64GB.\
Then compile the cpp files and run the executable file: (the g++ compiler is less efficient)
1. cd gcs_verifications
2. clang++ gcs.cpp gen_S4.cpp -std=c++20 -fopenmp -O3
3. ./a.out

The density of gcs quad is saved in density.txt.

## Verifications of the asymptotic existence of Hadamard matrices
Set two parameters in gcs_hadamard.cpp: \
int len_exp = 44; // maximum number $2^{44}$, restricted by the memory capacity, here it consumes about 2.6TB.\
int tab_max_exp = 36; // a tradeoff between the time complexity and the space complexity, here it consume $2^{36}$ bytes = 64GB.\
Then compile the cpp files and run the executable file:  (the g++ compiler is less efficient)
1. cd hadamard_verifications
2. clang++ gcs_hadamard.cpp gen_S4.cpp -std=c++20 -fopenmp -O3
3. ./a.out


#ifndef __GEN_S4_H__
#define __GEN_S4_H__
#include <unordered_set>
#include <set>
#include <math.h>
using namespace std;
class S4 {
    public:
        unordered_set<long long> set_large;
        bool *lookup_table;
        long long criteria;
        vector<long long> s4_vec;
        vector<long long>::iterator new_it;
        S4(int c) {
            criteria = c;
            long long tab_len = pow(2, criteria);
            lookup_table = (bool*)malloc(sizeof(bool)*(tab_len));
            lookup_table[0] = true;
            #pragma omp parallel for
            for (long long i=1; i<tab_len; i++) {
                lookup_table[i] = false;
            }
        }
        void update(vector<long long> &vec, long long len, bool *exist);
        bool check(long long n);
        void reset();
        void replace(vector<long long> &v);
};

void gen_S4(S4 &s4, vector<long long> &cp, long long n);

#endif

#include <stdio.h>
#include <vector>
#include <time.h>
#include <omp.h>
#include <math.h>
#include <algorithm>
#include "gen_S4.h"
using namespace std;


void gen_C(vector<long long> &C, long long n, int r) {
    for (int a=0; pow(2, a)<=n; a++) {
        long long ta = pow(2, a);
        for (int b=0; ta*pow(3, b)<=n; b++) {
            long long tb = ta*pow(3, b);
            for (int c=0; tb*pow(5, c)<=n; c++) {
                long long tc = tb*pow(5, c);
                for (int d=0; tc*pow(11, d)<=n; d++) {
                    long long td = tc*pow(11, d);
                    for (int e=0; td*pow(13, e)<=n; e++) {
                        int u = min(c+e, a);
                        if (b+c+d+e>a+u+r) continue;
                        C.push_back(td*pow(13, e));
                    }
                }
            }
        }
    }
    sort(C.begin(), C.end());
}

void gen_mx_sum(S4 &E, vector<long long> &C, long long n, int m) {
    vector<long long> vec(C.begin(), C.end());
    bool *exist = (bool*)malloc(sizeof(bool) * C.size());
    for (long long i=0; i<C.size(); i++) {
	exist[i] = true;
    }
    for (auto it1 = C.begin(); it1 != C.end(); it1++) {
        long long factor = n/m - (*it1);
        auto end = upper_bound(C.begin(), it1+1, factor);

        long long i = 0;
        for (auto it2 = C.begin(); it2 != end; it2++) {
            vec[i++] = m*(*it1+*it2);
        }
        E.update(vec, i, exist);
    }
    free(exist);
}

void gen_2x_square(S4 &s4, vector<long long>::iterator C_begin, vector<long long>::iterator C_end, long long n) {
    long long idx = s4.s4_vec.size();
    vector<long long> vec(distance(C_begin, C_end));
    bool* exist = (bool*)malloc(sizeof(bool) * vec.size());
    for (long long i=0; i<vec.size(); i++) {
	exist[i] = true;
    }
    for (auto it_C = C_begin; it_C != C_end; it_C++) {
        // all for parallel update
        long long factor = n/((*it_C)*2);
        auto end = upper_bound(C_begin, it_C+1, factor);

        long long i = 0;
        for (auto it_D = C_begin; it_D != end; it_D++) {
            vec[i++] = 2*(*it_C)*(*it_D);
        }
        s4.update(vec, i, exist);
    }
    free(exist);
}


template <typename it>
void gen_mx_product(S4 &s4, it C_begin, it C_end, it D_begin, it D_end, long long n, int m) {
    long long idx = s4.s4_vec.size();
    vector<long long> vec(distance(D_begin, D_end));
    bool* exist = (bool*)malloc(sizeof(bool) * vec.size());
    for (long long i=0; i<vec.size(); i++) {
	exist[i] = true;
    }
    for (auto it_C = C_begin; it_C != C_end; it_C++) {
        long long factor = n/((*it_C)*m);
        auto end = upper_bound(D_begin, D_end, factor);

        long long i = 0;
        for (auto it_D = D_begin; it_D != end; it_D++) {
            vec[i++] = m*(*it_C)*(*it_D);
        }
        s4.update(vec, i, exist);
    }
    s4.new_it = s4.s4_vec.begin()+idx;
    free(exist);
}

void gen_S4(S4 &s4, vector<long long> &cp, long long n) {
    long long num;
    // FILE* stream;
    //freopen("out.txt", "a", stdout);
    setvbuf(stdout, (char *)NULL, _IONBF, 0);
    printf("\n\nn = 2^%d", (int)log2(n));

// gen C and cp
    time_t start, end;
    time(&start);
    vector<long long> C;
    gen_C(C, n, 1);
    gen_C(cp, n, 2);
    time(&end);
    printf("\nC and cp: time=%ld; Len: %zd, %zd", end-start, C.size(), cp.size());

//  gen D
    time(&start);
    vector<long long> D;
    for (long long i=1; i<=77; i+=2) {
        D.push_back(i);
    }
    for (auto x : C) {
        num = 2*x+1;
        if (num>n || x<=38) continue;
        D.push_back(num);
    }
    sort(D.begin(), D.end());
    time(&end);
    printf("\nD: time=%ld; Len: %zd", end-start, D.size());

// 2(2CD U D U 2(cp+cp))^2
    time(&start);
    bool* exist = (bool*)malloc(sizeof(bool) * D.size());
    for (long long i=0; i<D.size(); i++) {
	exist[i] = true;
    }
    s4.update(D, D.size(), exist);

    gen_mx_product(s4, D.begin(), D.end(), C.begin(), C.end(), n, 2);
    sort(s4.s4_vec.begin(), s4.s4_vec.end());
    time(&end);
    printf("\n2CD: time=%ld; Len: %zd; Large_set size: %zd", end-start, s4.s4_vec.size(), s4.set_large.size());
    //printf("Min of large set: %lld", *min_element(s4.set_large.begin(), s4.set_large.end()));
    time(&start);
    gen_mx_sum(s4, cp, n, 2);
    time(&end);
    printf("\n2(cp+cp): time=%ld; Len: %zd; Large_set size: %zd", end-start, s4.s4_vec.size(), s4.set_large.size());
    time(&start);
    vector<long long> E;
    E = s4.s4_vec;
    sort(E.begin(), E.end());
    s4.reset();
    gen_2x_square(s4, E.begin(), E.end(), n);
    time(&end);
    printf("\n2(2CD U D U 2(cp+cp))^2: time=%ld; Len: %zd; Large_set size: %zd", end-start, s4.s4_vec.size(), s4.set_large.size());

// recursive update
    time(&start);
    vector<long long> F;
    s4.new_it = s4.s4_vec.begin();
    int i = 0;
    while (s4.new_it != s4.s4_vec.end()) {
        i++;
        auto F_size = F.size();
        F.insert(F.end(), s4.new_it, s4.s4_vec.end());
        sort(F.begin(), F.end());


        s4.replace(E);
        time(&end);
        printf("\nReplace with E%d: time=%ld", i - 1, end - start);

        gen_mx_product(s4, C.begin(), C.end(), F.begin() + F_size, F.end(), n, 2); // The second vector is ordered
        sort(s4.new_it, s4.s4_vec.end());
        vector<long long> vec(s4.new_it, s4.s4_vec.end());
        time(&end);
        printf("\nE%d: time=%ld; Len of Enew%d: %zd", i, end - start, i, vec.size());

        s4.replace(F);
        time(&end);
        printf("\nReplace with F%d: time=%ld", i - 1, end - start);
        if (vec.size() == 0) break;

	auto idx = s4.s4_vec.size();
        gen_2x_square(s4, vec.begin(), vec.end(), n);
        time(&end);
        printf("\nF%d from square: time=%ld; Len: %zd; Large_set size: %zd", i, end - start, s4.s4_vec.size(), s4.set_large.size());
        if (vec.size() < E.size()) gen_mx_product(s4, vec.begin(), vec.end(), E.begin(), E.end(), n, 2); // E is ordered
        else gen_mx_product(s4, E.begin(), E.end(), vec.begin(), vec.end(), n, 2);
	s4.new_it = s4.s4_vec.begin()+idx;

        time(&end);
        printf("\nF%d heterogenous: time=%ld; Len: %zd; Large_set size: %zd", i, end - start, s4.s4_vec.size(), s4.set_large.size());
        E.insert(E.end(), vec.begin(), vec.end());
        sort(E.begin(), E.end());
    }
    vector<long long>().swap(E);
    vector<long long>().swap(F);
    time(&end);
    printf("\nQuad: time=%ld; Len: %zd; Large_set size: %zd", end - start, s4.s4_vec.size(), s4.set_large.size());
    


// add D and C+C to s4
    time(&start);
    s4.update(D, D.size(), exist);
    free(exist);

// S4 X C 
    time(&start);
    vector<long long> vec = s4.s4_vec;
    sort(vec.begin(), vec.end());
    gen_mx_product(s4, C.begin() + 1, C.end(), vec.begin(), vec.end(), n, 1);
    vector<long long>().swap(vec);
    time(&end);
    printf("\nS4 x C, time=%ld; Len: %zd; Large_set size: %zd", end - start, s4.s4_vec.size(), s4.set_large.size());
    

//  CxC + CxC
    time(&start);
    gen_mx_sum(s4, cp, n, 1);
    time(&end);
    printf("\ncp + cp, time=%ld; Len: %zd; Large_set size: %zd", end-start, s4.s4_vec.size(), s4.set_large.size());
    sort((s4.s4_vec).begin(), (s4.s4_vec).end());
    return;
}


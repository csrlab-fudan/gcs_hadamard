#include <vector>
#include <iostream>
#include <fstream>
#include <time.h>
#include <omp.h>
#include "gen_S4.h"
#include <math.h>
#include <algorithm>
#include <stdio.h>
using namespace std;


void S4::update(vector<long long> &a, long long len, bool *exist) {
    auto s4_vec_size = s4_vec.size();
    #pragma omp parallel for num_threads(16)
    for (long long i = 0; i < len; i++) {
        if (!check(a[i])) exist[i] = false; // dirty write
        //else exist[i] = false;
    }
    for (auto i = 0; i < len; i++) {
        if (exist[i]) continue;
	    exist[i] = true; // recover
        s4_vec.push_back(a[i]);
        if (a[i] >> criteria > 0) set_large.insert(a[i]);
        else lookup_table[a[i]] = true;
    } 
}

bool S4::check(long long n) {
    if (n>>criteria > 0) {
        return set_large.contains(n);
    }
    else {
        return lookup_table[n];
    }    
}

void S4::reset() {
    set_large.clear();
    // set_large.reserve(pow(2, 31));
    lookup_table[0] = true;
    long long imax = pow(2, criteria);
    #pragma omp parallel for 
    for (long long i = 1; i < imax; i++) {
        lookup_table[i] = false;
    }
    s4_vec.clear();
}

void S4::replace(vector<long long> &v) {
    reset();
    s4_vec = v;
    new_it = s4_vec.begin();
    for (auto x : v) {
        if (x>>criteria > 0) set_large.insert(x);
        else lookup_table[x] = true;
    }  
}

int main() {
    setvbuf(stdout, (char *)NULL, _IONBF, 0);
    time_t start, end;
    time(&start);
    long long n = (long long)1e13;
    int tab_max_exp = 36; // maximum tab len exp
    S4 s4(tab_max_exp);
    
    vector<long long> cp;
    gen_S4(s4, cp, n);
// quad
    // FILE* stream;
    //freopen("out.txt", "a", stdout);
    setvbuf(stdout, (char *)NULL, _IONBF, 0);
    cout << endl << "Level 2: ";
    long long num = 0;
    vector<pair<long long, double>> density;
    for (long long j=1; j<=1e4; j++) {
        if (s4.check(j)) {
            num++;
            density.push_back(make_pair(j, double(num)/double(j)));
        }
	}

    long long steps = (n-1e4)/1e4;
    for (long long j = 1e4 + steps; j<=n; j += steps) {
         auto it = upper_bound(s4.s4_vec.begin(), s4.s4_vec.end(), j);
         auto dis = distance(s4.s4_vec.begin(), it);
         density.push_back(make_pair(j, double(dis)/double(j)));
    }
    
    ofstream f("density.txt", ios::out);
    for (auto [x, y] : density) f << x << " " << y << endl;
    f.close();

    cout << endl << "Level 3: ";
    int percent = 0;
    time(&start);
    long long eight_max = n;
    #pragma omp parallel for schedule(static, 8)
    for (long long j=1; j<n+1; j++) {
        if (j%(n/100)==0) {
            #pragma omp critical
            {time(&end);
            cout << "  " << ++percent << "%, time=" << end-start << ";";}
        }
        if (j>eight_max) continue;
        
        auto it = lower_bound(s4.s4_vec.begin(), s4.s4_vec.end(), j);
        if (it != s4.s4_vec.end() && *it == j) continue;
        std::reverse_iterator<decltype(it)> rit{it};
        rit++;
        bool found = false;
        for (; rit != s4.s4_vec.rend() && *rit >= j>>1 && !found; rit++) {
            found = s4.check(j-*rit);
        }

        if (!found) {
            #pragma omp critical
            eight_max = min(eight_max, j-1);
        }

    }

    cout << endl << "S8: Maximum Num: " << eight_max;

    free(s4.lookup_table);
    return 0;
}

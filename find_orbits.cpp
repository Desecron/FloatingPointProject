/* 
 * Hongyan Wang
 * Dec 20, 2018
 */

#include <iostream>
#include <fstream>
#include <bitset>
#include <climits>
#include <stdint.h>
#include <unordered_map>
#include <map>
#include <math.h>  
#include <vector>
#include <time.h>
#include <stdio.h>
#include <algorithm>
#include <iomanip>

using namespace std;

/* range of valid float numbers */
#define UPPERBOUNDARY 1065353215 // 0-01111110-11111111111111111111111
#define LOWERBOUNDARY 0
#define SMALLESTDENORM 1
#define LARGESTDENORM 8388607
#define PI 3.141592653589793

/* representation of single-precision float number*/
struct FloatNumber{
    uint32_t e;
    uint32_t m;
};

/* show binary bits */
template<typename T>
void showBinrep(const T& a)
{
    const char* beg = reinterpret_cast<const char*>(&a);
    const char* end = beg + sizeof(a);
    while(beg != end)
        cout << bitset<CHAR_BIT>(*beg++) << ' ';
    cout << '\n';
}

/* transfer int to float and keep all bits unchanged */
float intToFloat(const uint32_t a){
    return *((float*)&a);
}

/* transfer float to int and keep all bits unchanged */
uint32_t floatToInt(const float b){
    return *((uint32_t*)&b);
}

/* keep only values between 0 and 1 */
bool checkValid(const uint32_t num){
    float temp = intToFloat(num);
    if((temp >= 0) && (temp < 1.0)){
        return true;
    } else {
        return false;
    }
}

/* multiply by beta, module 1. Be careful about rounding */
float fTrans(const float num, const double beta1, const double beta2){
    float result = fmod(num * beta1 / beta2, 1.0);
    if (result == 1) {
        return 0.0;
    } else {
        return result;
    }
}

/* construct a float number from m and e */
float bitsToFloat(const FloatNumber f){
    uint32_t e = f.e;
    uint32_t m = f.m;
    uint32_t num = 0;
    uint32_t mask = e << 23;
    num = num | mask;
    mask = m;
    num = num | m;
    return intToFloat(num);
}

/* extract m and e from a float number */
FloatNumber floatToBits(const float num){
    uint32_t num_int = floatToInt(num);
    uint32_t e = (num_int >> 23);
    uint32_t mask = 1;
    mask = (mask << 23) - 1;
    uint32_t m = num_int & mask;
    FloatNumber result;
    result.e = e;
    result.m = m;
    return result;
}


/* Test functions */
double xPowFun(const double num, const int n){
    return pow(num,n);
}

/*
 * orbits: save all the orbits
 * transients: 
 * myflag: records which orbit each number belongs to
 * beta1, beta2: transform function parameter, beta = beta1/beta2
 */
void getOrbits(vector< vector<uint32_t> >& orbits, unordered_map<uint32_t, double>& transients, uint32_t *myflag, const double beta1, const double beta2) {

    /* flag means orbit number */
    uint32_t count, flag = 0; 

    /* Find all the orbits */
    for(count = LOWERBOUNDARY; count <= UPPERBOUNDARY; count++){

        /* skip the subnormal numbers */
        if (count >= SMALLESTDENORM && count <= LARGESTDENORM) {
            continue;
        }

        uint32_t startPoint = count;

        vector<uint32_t> currOrbit;

        /* it means we already visited this orbit */
        if(myflag[startPoint] > 0){ 
            continue;
        } 

        flag++;

        uint32_t currPoint = startPoint;

        while(myflag[currPoint] == 0){
            myflag[currPoint] = flag;

            FloatNumber currPointE = floatToBits(intToFloat(currPoint));
            int eNum = currPointE.e;
            transients[flag] = transients[flag] + pow(2, eNum - 150);

            const float currPointFloat = fTrans(intToFloat(currPoint), beta1, beta2);
            currPoint = floatToInt(currPointFloat);
            if (currPoint >= SMALLESTDENORM && currPoint <= LARGESTDENORM) {
                cout<<"Subnormal numbers appear!"<<endl;
            }
        }

        if(myflag[currPoint] == flag){ // We reach the same orbit

            /* The first element in the orbits vector represents the orbit number */
            currOrbit.push_back(flag);


            uint32_t tmpPoint = currPoint;
            currOrbit.push_back(currPoint);
            const float currPointFloat = fTrans(intToFloat(currPoint), beta1, beta2);
            currPoint = floatToInt(currPointFloat);
            
            while(currPoint != tmpPoint){
                currOrbit.push_back(currPoint);
                const float currPointFloat = fTrans(intToFloat(currPoint), beta1, beta2);
                currPoint = floatToInt(currPointFloat);
            }
            orbits.push_back(currOrbit);
        } else {
            // all the points belong to a different orbit
            uint32_t realFlag = myflag[currPoint];
            uint32_t tmpPoint = currPoint;
            currPoint = startPoint;
            while(currPoint != tmpPoint){
                myflag[currPoint] = realFlag;

                FloatNumber currPointE = floatToBits(intToFloat(currPoint));
                int eNum = currPointE.e;
                transients[flag] = transients[flag] - pow(2, eNum - 150);
                transients[realFlag] = transients[realFlag] + pow(2, eNum - 150);

                const float currPointFloat = fTrans(intToFloat(currPoint), beta1, beta2);
                currPoint = floatToInt(currPointFloat);
            }
            flag--;
        }
    }
}

static bool compare (const vector<uint32_t>& a, const vector<uint32_t>& b) {
    return a.size() < b.size();
}

int main()
{   
    auto t = clock();
    vector< vector<uint32_t> > orbits; // save different orbits
    unordered_map<uint32_t, double> transients; 

    uint32_t *myflag = new uint32_t[UPPERBOUNDARY - LOWERBOUNDARY + 1]; // record each integer belongs to which orbit
    memset(myflag, 0, sizeof(myflag));
    const double beta1 = 15.0;
    const double beta2 = 1.0;

    getOrbits(orbits, transients, myflag, beta1, beta2); 
    
    // sort the orbits by the orbits length
    sort(orbits.begin(), orbits.end(), compare);  

    double total_measures = 0.0;

    for (auto it : transients) {
        if (it.second == 0) {
            continue;
        }
        total_measures += it.second;
    }

    cout<<"The total measures is "<<total_measures<<endl;

    unordered_map<uint32_t,vector<double>> divided_transients;
    unordered_map<uint32_t, double> orbit_measures;
    for(int currPoint = LOWERBOUNDARY; currPoint <= UPPERBOUNDARY; currPoint++){
        // skip the subnormal numbers
        if (currPoint >= SMALLESTDENORM && currPoint <= LARGESTDENORM) {
            continue;
        }
        FloatNumber currPointE = floatToBits(intToFloat(currPoint));
        int eNum = currPointE.e;
        divided_transients[myflag[currPoint]].push_back(pow(2, eNum - 150));
    }
    for (int i = 0; i < orbits.size(); i++) {
        int orbit_number = orbits[i][0];
        sort(divided_transients[orbit_number].begin(),divided_transients[orbit_number].end());
        double tmp_measure_sum = 0.0;
        for (int k = 0; k < divided_transients[orbit_number].size(); k++) {
            tmp_measure_sum += divided_transients[orbit_number][k];
        }
        orbit_measures[orbit_number] = tmp_measure_sum;
    }

    /*Write to csv file*/
    ofstream myfile;
    myfile.open("beta15.csv");
    for (int q = 1; q <= 100; q++) {
        for(int i = 0; i < orbits.size(); i++) {
            double sum = 0.0;
            // notice j starts from 1 because the first element represents orbit number
            for (int j = 1; j < orbits[i].size(); j++) {
                sum += xPowFun((double)intToFloat(orbits[i][j]), q);
            }
            myfile<<q<<","<<(i+1)<<","<<orbits[i].size()-1<<","<<sum<<","<<orbit_measures[orbits[i][0]]/total_measures<<"\n";
        }
    }
    myfile.close();

    delete[] myflag;

    t = clock() - t;
    printf ("It took me %lu clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);

    return 0;
} 
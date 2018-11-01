#ifndef _REASON_
#define _REASON_


#define PI 3.141592653589

#include<map>
#include<set>
#include<cmath>
#include<deque>
using namespace std;


vector<double> getSigma2();
void getSampleData(deque<vector<double>>&dataList);
void getSample(vector<double>Sample);
vector<double>getParent(int x);
vector<int>getOmegRank(vector<double> parentNode);
vector<int> getSumSample(const vector<double>&tauDiscreteList, const vector<int>&tauLoopNumber);
vector<int> getOmegaMinusTao(const vector<int>&tauLoopNumber, const vector<int>&omegaNodeNumber);
vector<vector<double>>getProbabilityData(const vector<int>&tauSumNumberList);
map<vector<double>, int>getProbability(const vector<int>&tauSumNumberList, const vector<double>&tauDiscreteList, const vector<int>&tauLoopNumber, const vector<int>&omegaNodeNumber);
double xZero(map<double, double>&muBeitaDict, double sigma2);
void loopMain(const vector<int>&loopNumber, const vector<int>&omegaNodeNumber, const vector<double>&sampleList, const vector<double>&parentNode, map<vector<int>, double>&tauCList, int k);
double Discrete(double number, int col);
void loopHead(vector<int>omegaNodeNumber, vector<double>sampleList, vector<double>parentNode, int k);
void loopStart(int i);
void showPath();

#endif
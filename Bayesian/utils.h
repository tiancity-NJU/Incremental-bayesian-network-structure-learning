#ifndef _UTIL_
#define _UTIL_


#include <stdio.h>
#include <iostream>
#include <string>
#include<vector>
#include<deque>
#include<cmath>
#include<fstream>  
#include<stdlib.h>
#include<string>   
#include<stack>
#include<sstream> 
#include<assert.h>
#include<Windows.h>
using namespace std;




typedef __int64 int64_t;
typedef long long int64t;
int64_t GetTime();





extern vector<double>meanList;
extern vector<double>standardList;    


double S2D(const string str);
deque<vector<double>> loadtxt(const char* filename, int d);
vector<vector<double>> loadShowB(const char* filename, int d);
double Shape(vector<double> lineStr);
double Sign(double op);
void GetColData(const deque<vector<double>>&rawX, unsigned i, vector<double>&data);
void Scale(const deque<vector<double>>& rawX, deque<vector<double>>&X);
vector<vector<double> > Zero(int x, int y);
void RemoveCol(deque<vector<double>>&data, int d);
void RemoveCols(deque<vector<double>>&data, int dim1, int dim2);
void RemoveElem(const vector<vector<double>>& B, int i, int j, vector<double>&result);
double Norm(vector<double> sample);
void Dot_vv(const deque<vector<double>>&a, const vector<double>&b, vector<double>&result);
double Dot_av(const vector<double>&a, const vector<double>&b);
void GetCol_S(const deque<vector<double>>&data, int d, vector<double>&result);
void GetCol_M(const vector<vector<double>>&data, int d, vector<double>&result);
void Dot_vv(const deque<vector<double>>&a, const vector<double>&b, vector<double>&result, int col1, int col2); //将矩阵a的col1和col2列排除，然后和b进行矩阵乘法，结果放到result里面

/*推理的相关辅助函数*/
vector<vector<int>>Combinations(vector<int>list, int len);


#endif

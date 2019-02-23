#include"utils.h"
#include"SparseBN.h"
#include"Reason.h"
#include"batchSparseBN.h"
#include<map>
#include <iomanip>
using namespace std;

map<double, double> loadMap(const char* filename)
{
	ifstream inFile(filename, ios::in);
	if (!inFile.is_open())
	{
		cout << "Error opening file"; exit(1);
	}

	string lineStr;
	map<double, double> strMap;
	while (getline(inFile, lineStr))	//以行为单位读入lineStr
	{
		stringstream ss(lineStr);
		string key, value;
		getline(ss, key, ':');	//以，为分割读ss入str
		getline(ss, value);
		strMap[S2D(key)] = S2D(value);
	}
	return strMap;
}

int main(){
	
	//cout << "读取时间..." << ((double)end - (double)start)/1000<< endl;
	
	//SBN("E:\\fangtian\\S410039\\DataDim57-new.csv", 56, 2, 0.1);     //初始化，并且跳过csv的第一行和第一列（一个是特征名，一个是日期）

	/*    SBN(path,dim,penalty,ralationthreshold)       */

	string path;
	int dims=0;
	double penalty=1.0, relation=0.1;
	string param_txt;
	string sigma2_txt;
	string edge_txt;
	cout << "file path: ";
	cin >> path;
	cout << "dim: ";
	cin >> dims;
	cout << "penalty param:";
	cin >> penalty;
	cout << "relation threshold: ";
	cin >> relation;
	cout << "param matrix save path:(such as：E:\\param.txt)";
	cin >> param_txt;
	cout << "sigma2 save path:(such as：E:\\sigma2.txt)";
	cin >> sigma2_txt;
	cout << "edge save path:(such as: E:\\edge.txt)";
	cin >> edge_txt;
	
	cout << "loading data....." << endl;
	SBN((char*)path.data(), dims, penalty, relation);


	Train();
	cout << "*******************************computing sigma2*****************************************" << endl;
	vector<double>sigma2=getSigma2();

	cout << "param is      " << path.data() << "  dim is" << dim << " penalty " << penalty << " relation threshold  " << relation<<endl;

	cout << "**********************************DAG**********************************" << endl;
	printDAG();
	cout << endl << endl;
	cout << "************************************B*******************************" << endl;
	printB();
	cout << endl << endl;
	cout << "**********************************showB********************************" << endl;
	printShowB();

	int edge = 0;
	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			if (DAG[i][j] == 1)
				edge++;
		}
	}

	ofstream MatrixB;
	MatrixB.open(param_txt);


	ofstream Sigma;
	Sigma.open(sigma2_txt);

	ofstream Edge;
	Edge.open(edge_txt);

	for (int i = 0; i < dim - 1; i++)
	{
		for (int j = 0; j < dim-1; j++)
		{
			MatrixB << showB[i][j] << ",";
		}
		MatrixB << showB[i][dim-1];
			MatrixB << "\n";
	}

	for (int i = 0; i < sigma2.size()-1; i++)
	{
		Sigma << sigma2[i] << ",";
	}
	Sigma << sigma2[sigma2.size() - 1];

	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			if (DAG[i][j] == 1)			// print('存在边 ', i + 1, j + 1)
				Edge << "G(" << i + 1 << ")" << "(" << j + 1 << ")" << "= 1" << endl;
		}
	}
	MatrixB.close();
	Sigma.close();
	Edge.close();

	cout << "**********************************统计结果为*******************************" << endl;
	cout << "node: " << dim << "edge: " << edge << endl;
	cout << "键入任意数字终止......" << endl;
	int tmp;
	cin >> tmp;

	/**
	cout << "输入任意数字进入推断" << endl;
	cin >> qq;
	
	cout << "***********************************下面进入推断**************************************" << endl;
	
	getSample(rawX[10]);    //将原始数据的一项投入，得到正则化交付给sampleList
	cout << "*******************************成功获取待检测样本*************************************" << endl;
	getSampleData(X);      //  将正则化的训练样本加上待检测样本再加上两个最值，离散化得到DiscreX
	cout << "******************************************成功离散化**********************************" << endl;
	getSigma2();
	cout << "******************************************成功获取sigma*******************************" << endl;
	loopStart(0);

	cout << "**********************************异常路径*************************************" << endl;
	showPath();
	**/
	
}


#include"BatchSparseBN.h"
#include"SparseBN.h"
#include"utils.h"

using namespace std;


double batchLamda1;
double batchLamda2;
double batchRelationThreshold;


int batchNumber;   //批数

void Split(deque<vector<double>>&rawData, deque<vector<double>>&Data, int start, int end)    //将全样本分成对应位置的批次
{
	assert(start >= 0);
	assert(end < rawData.size());
	Data.clear();
	for (int i = start; i < end; i++)
	{
		Data.push_back(rawData[i]);
	}
}



void BatchSBN(char* filename, int d, double penalty, double relation, int batch)
{

	dim = d;
	batchNumber = batch;
	batchRelationThreshold = relation;
	rawX = loadtxt(filename, dim);

	Split(rawX, X, 0, rawX.size() / batch);    //取前十分之一批次

	Scale(rawX, X);   //正则化
	cout << "第一个样本是.............." << endl;
	for (int i = 0; i < dim; i++)
	{
		cout << X[0][i] << " ";
	}
	num = X.size();
	B = Zero(dim - 1, dim);     //(p-1)*p 的
	showB = Zero(dim - 1, dim);         //用来显示 DAG 用的 B，在原始B上淘汰关系阈值以下的，并进行相应剪枝，总之 DAG 和这个 showB 是对应的，后面的推理也是基于这个showB， 不直接对B剪枝是消除B在增量更新时候的误差累计
	batchLamda1 = penalty*(0.12*num);    			//penalty
	batchLamda2 = (((num - 1)*(num - 1)*dim) / batchLamda1 - batchLamda1)*1.2;

	DAG = Zero(dim, dim);     //DAG
	P = Zero(dim, dim);
}

vector<double> BatchBFS(int i){
	//得到DAG中i为起点的闭包
	for (int j = 0; j < dim; j++)
		P[i][j] = 0;
	queue<double> q;
	q.push(i);
	while (!q.empty())
	{
		double node = q.front();
		q.pop();
		if (P[i][node] == 0)
		{
			P[i][node] = 1;
			for (int j = 0; j < dim; j++)
			{
				if (DAG[node][j] == 1 && P[i][j] == 0)
					q.push(j);
			}
		}
	}
	return P[i];    //更新P矩阵一列，确保闭包
}


void BatchRefreshDAG(int i){                       //BCD 对50维每一维迭代的时候都刷新对应行的DAG，这样BFS才有效
	for (int j = 0; j < dim; j++)
		DAG[j][i] = 0;                              //对应列清空，重新刷新一遍    为保险，B不断增量，DAG修改时每次全部刷新

	for (int j = 0; j < (dim - 1); j++)
	{
		if (j >= i)
		{
			if (B[j][i] != 0)
				DAG[j + 1][i] = 1;

		}
		else
		{
			if (B[j][i] != 0)
				DAG[j][i] = 1;
		}
	}
}


vector<double> BatchGetBij(int i, int j)
{
	vector<double>result;
	result.clear();
	RemoveElem(B, i, j, result);
	return result;
}



void BatchGetXij(deque<vector<double>>&result, int dim1, int dim2)
{
	int k;
	if (dim1 < dim2){ k = dim1; dim1 = dim2; dim2 = k; }
	for (unsigned i = 0; i < result.size(); i++)
	{
		result[i].erase(result[i].begin() + dim1);
		result[i].erase(result[i].begin() + dim2);
	}
	return;
}


void BatchShooting(int i, int shootIter)
{
	vector<double> xi;

	GetCol_S(X, i, xi);


	for (int t = 0; t < shootIter; t++)     //最大迭代次数
	{
		
		vector<double>LastColB;
		GetCol_M(B, i, LastColB);
		double temp = Norm(LastColB);     //每次循环把上一次的Bi 存到temp中
		for (int j = 0; j < dim; j++)    //J 的每一趟循环，更新B（ (p-1)*p  ）的一列 ， 即 p-1 个元素， 特别的，当j>i的时候，索引要-1，因为占据了本来i的位置，才变成p-1维
		{

			if (j == i)
				continue;    //保证每次j的迭代运算 (p-1) 次

			//这是   n*(dim-2)*(dim-2)*1的慢速计算，现在用的Dot_vv进行了优化，又进行了30%的加速
			//deque<vector<double>> test(X);  getXij(test, i, j);     //   test 是X去除两列以后的结果      20ms
			//vector <double>xb;  Dot_vv(test, getBij(i, j), xb);      //50ms

			vector<double> tmp;

			vector<double>xb;  Dot_vv(X, BatchGetBij(i, j), xb, i, j);


			for (unsigned k = 0; k < xi.size(); k++)
				tmp.push_back(xi[k] - xb[k]);			 //(n*1) - (n*48).*(48*1)    得到的分子 tmp shape is （n*1）


			vector<double> xj;
			GetCol_S(X, j, xj);     					//xj  (n*1)
			double part1 = Dot_av(tmp, xj);

			double Denominator = Dot_av(xi, xi);

			if (Denominator == 0){             //  个别特征自乘全为0，确保分母不会变为0
				part1 = part1 / num;
			}
			else{
				part1 = part1 / Denominator;                //该part1 对应B的第一部分
			}


			//  part2不怎么耗时间
			double part2;
			if (Denominator == 0){
				part2 = (batchLamda1 + batchLamda2*P[i][j]) / num;
			}
			else{
				part2 = (batchLamda1 + batchLamda2*P[i][j]) / Denominator;
			}

			double sign = Sign(part1);

			double result;
			if (fabs(part1) - part2 < 0)
				result = 0;
			else
				result = (fabs(part1) - part2)*sign;
			// 在 shooting 训练里面先不管关系阈值等问题
			if (j < i)
				B[j][i] = result;
			else
				B[j - 1][i] = result;

		}
		GetCol_M(B, i, LastColB);
		double newtemp = Norm(LastColB);
		//std::cout << "对于shooting里的" << i << "次迭代，Bi变成了: " << newtemp << endl;
		if (fabs(newtemp - temp) < 0.001)      //如果差值小于0.01 则提前退出迭代，最多100次迭代
			break;
	}

}


void BatchBCD(int BCDIter)
{
	for (int t = 0; t < BCDIter; t++)
	{
		cout << "进入第"<<t<<"轮BCD迭代................." << endl;
		for (int i = 0; i < dim; i++)		 //反复迭代更新 Bi
		{
			vector<double> Pij = BatchBFS(i);
			cout << "第" << t << "轮BCD的第" << i << "维迭代............." << endl;
			BatchShooting(i, 100);
			BatchRefreshDAG(i);    //更新DAG 中 i  对应的那一列
		}
	}
}


void BatchCheckAndPrune()  // 对 B 剪枝
{
	for (int i = 0; i < (dim - 1); i++)
	{
		for (int j = 0; j < dim; j++)
			showB[i][j] = B[i][j];              // 先把原始的 B 移植到 showB 上来
	}


	for (int i = 0; i < (dim - 1); i++)
	{
		for (int j = i + 1; j < dim; j++)
		{
			if (fabs(showB[i][j]) > fabs(showB[j - 1][i]))
				showB[j - 1][i] = 0;
			else
				showB[i][j] = 0;
		}
	}


	for (int i = 0; i < (dim - 1); i++)  // 关系阈值判定
	{
		for (int j = 0; j < dim; j++)
			if (fabs(showB[i][j]) < batchRelationThreshold)
				showB[i][j] = 0;
	}

	for (int i = 0; i < dim; i++)		// 先清空  DAG
	{
		for (int j = 0; j < dim; j++)
			DAG[i][j] = 0;
	}

	for (int i = 0; i < (dim - 1); i++)		 // 最后再重新刷新一遍 DAG
	{
		for (int j = 0; j < dim; j++)
		{
			if (showB[i][j] != 0)
			{
				if (i >= j)
					DAG[i + 1][j] = 1;
				else
					DAG[i][j] = 1;
			}
		}
	}
}


void BatchTrain()
{
	time_t start, end;
	time(&start);
	BatchBCD(5);   //先对第一批进行训练
	
	for (int i = 1; i < batchNumber; i++)
	{
		cout << "*******************************开始第" << i + 1 << "批数据的训练*********************************" << endl;
		//Split(rawX, X, (int)(i*(int)rawX.size() / batchNumber), (int)((i + 1)*(int)rawX.size() / batchNumber));
		Split(rawX, X, (int)((i - 1)*num + (num / 2)),(int)(i*num + (num / 2)));
		
		Scale(rawX, X);

		BatchBCD(2);
		cout << "***********************************新一轮**************************" << endl;
	}
	
	cout << "****************************************最后对总体样本再次循环一次***********************************" << endl;
	X.clear();
	X = rawX;
	Scale(rawX, X);

	batchLamda1 = batchLamda1*batchNumber;
	num = X.size();
	batchLamda2 = (((num - 1)*(num - 1)*dim) / batchLamda1 - batchLamda1)*1.1;


	BatchBCD(10);
	time(&end);
	cout << "Batch训练时间..................................................." << difftime(end, start) << "秒" << endl;
	BatchCheckAndPrune();
	
}


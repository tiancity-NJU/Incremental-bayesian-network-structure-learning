#include"SparseBN.h"
#include"utils.h"
#include"Reason.h"
using namespace std;

map<int, vector<int>>path;    //路径字典

vector<vector<double>>DiscreX;     //对应python版本的dataListNew

vector<double>sigma2;     //从SpareBN中获取simga2


vector<double>sampleList; //待检测样本




void getSample(vector<double>Sample)      //将到来的待检测样本正则化,  离散化是全体的工作
{
	double temp;
	for (int i = 0; i < dim; i++)
	{
		temp = (Sample[i] - meanList[i]) / standardList[i];
		sampleList.push_back(temp);
	}
}

void getSampleData(deque<vector<double>>&dataList)    //均匀离散化，十等分
{
	vector<double>Max;  //存储各个维度的最值
	vector<double>Min;
	int i, j, x;
	double FeatMax, FeatMin;

	for (i = 0; i < dim; i++)
	{
		FeatMax = dataList[0][i];
		FeatMin = dataList[0][i];
		for (j = 0; j < dataList.size(); j++)
		{
			if (dataList[j][i]>FeatMax) FeatMax = dataList[j][i];
			if (dataList[j][i]<FeatMin) FeatMin = dataList[j][i];
		}
		if (sampleList[i]>FeatMax) FeatMax = sampleList[i];
		if (sampleList[i]<FeatMin) FeatMin = sampleList[i];

		Max.push_back(FeatMax);     //各维特征的最大值
		Min.push_back(FeatMin);     //各维特征的最小值
	}

	vector<double>DiscreSample;
	for (i = 0; i < dataList.size(); i++)   //对于每一个样本i
	{
		for (j = 0; j < dataList[i].size(); j++)   //对于每一个特征 j
		{
			for (x = 0; x < 10; x++)
			{
				if (Min[j] + (x + 1)*(Max[j] - Min[j]) / 10 >= dataList[i][j])
				{
					DiscreSample.push_back(x + 1);
					break;
				}
			}
		}
		DiscreX.push_back(DiscreSample);   //获得离散化的数据样本
		DiscreSample.clear();
	}

	DiscreSample.clear();
	
	for (i = 0; i < sampleList.size(); i++)   //对于每一个特征 j
	{
		for (x = 0; x < 10; x++)
		{
			
			if (Min[i] + (x + 1)*(Max[i] - Min[i]) / 10 >= sampleList[i])
			{
				DiscreSample.push_back(x + 1);
				break;
			}
		}
	}
	DiscreX.push_back(DiscreSample);   //获得离散化的数据样本

	
	sampleList = DiscreSample;
	// 再把各维最值作为最后两行加入其中    所以DiscreX是一个[num+1（待检测）+2（最值）]*dim的一个矩阵
	DiscreX.push_back(Max);
	DiscreX.push_back(Min);
	
	//cout << "正则化数据为***********************************************************" << endl;
	//for (int k = 0; k < DiscreX.size(); k++)
	//{
	//	for (int j = 0; j < DiscreX[k].size(); j++)
	//	{
	//		cout << DiscreX[k][j] << ' ';
	//	}
	//	cout << DiscreX[k].size();
	//	cout << endl;
	//}
	
}


void loopStart(int i)
{
	cout << "进入了loopStart....................." << i<<endl;

	vector<double>parentNode = getParent(i);    //获得showB矩阵中i对应列的扩展向量（50*1）

	vector<int>omegaNodeNumber = getOmegRank(parentNode);
	if (omegaNodeNumber.size() > 0){
		loopHead(omegaNodeNumber, sampleList, parentNode, i);
	}
	else{
		vector<int>p;
		p.push_back(-1);
		path[i] = p;       //没有父亲的话   置-1
		cout << "结点" << i << "无父节点" << endl;
	}
}


vector<double>getParent(int x)    //获取B矩阵特定列构成的向量,自己也要扩充，补成dim维
{

	vector<double>list1;
	for (int i = 0; i < dim; i++)
	{
		if (i < x) list1.push_back(showB[i][x]);
		else if (i == x) list1.push_back(0.0);
		else list1.push_back(showB[i - 1][x]);
	}

	for (int i = 0; i < list1.size(); i++)
	{
		cout << list1[i] << ' ';
	}
	cout << endl;
	return list1;
}


vector<int>getOmegRank(vector<double> parentNode)
{
	vector<int>parentNodeNumber;
	for (int i = 0; i < parentNode.size(); i++)
	{
		if (parentNode[i] != 0)
			parentNodeNumber.push_back(i);
	}
	return parentNodeNumber;
}

void loopHead(vector<int>omegaNodeNumber, vector<double>sampleList, vector<double>parentNode, int k)    //omegaNodeNumber是int型的父亲集合 sampleList是待检测样本矩阵，parentNode是showB中对应列的数据
{
	cout << endl;
	cout << "************************追踪第" << k << "维的影响父节点*****************************" << endl;
	if (path.count(k) != 0){
		cout << "---------------" << k << "维已经遍历-----------------" << endl;
		return;
	}
	vector<double>a;
	map<vector<int>, double>tauCList;    //具体不清楚
	for (int i = 0; i < omegaNodeNumber.size(); i++)
	{
		vector<vector<int>>combins = Combinations(omegaNodeNumber, i + 1);
		cout << "计算组合C" << omegaNodeNumber.size() << ',' << i + 1 << endl;
		for (int j = 0; j < combins.size(); j++)
		{
			for (int p = 0; p < combins[j].size(); p++) cout << combins[j][p] << ' ';
			cout << endl;
		
			loopMain(combins[j], omegaNodeNumber, sampleList, parentNode, tauCList, k);
			
		}
	}

	//  下面离散化
	double jDiscrete = Discrete(sampleList[k], k);
	double min = abs(tauCList.begin()->second - jDiscrete);



	vector<int> keyValue;

	for (map<vector<int>, double>::iterator it = tauCList.begin(); it != tauCList.end(); it++)   //寻找一个什么最小值
	{
		
		if (abs(it->second - jDiscrete) <= min)
		{
			min = abs(it->second - jDiscrete);
			keyValue = it->first;
		}
		
	}
	path[k] = keyValue;       
	
	for (int i = 0; i < keyValue.size(); i++)
	{
		if (keyValue[i] >= dim)
		{
			break;
		}
		loopStart(keyValue[i]);
	}
}


void loopMain(const vector<int>&loopNumber, const vector<int>&omegaNodeNumber, const vector<double>&sampleList, const vector<double>&parentNode, map<vector<int>, double>&tauCList, int k)
{
	vector<int>tauLoopNumber;    //  combins[j]对应omegaNodeNumber里的坐标
	vector<double>tauDiscreteList;   //  对应的tau离散化数值
	for (int i = 0; i < loopNumber.size(); i++)
	{
		tauLoopNumber.push_back(omegaNodeNumber[loopNumber[i]]);
		
	}


	cout << sampleList.size() << endl;
	cout << endl;
	for (int i = 0; i < tauLoopNumber.size(); i++)
	{
		tauDiscreteList.push_back(sampleList[tauLoopNumber[i]]);                   //  tauLoopNumber 太大了
	}

	vector<int>tauSumNumberList = getSumSample(tauDiscreteList, tauLoopNumber);



	vector<int>omegaMinusTao = getOmegaMinusTao(tauLoopNumber, omegaNodeNumber);

	
	map<vector<double>, int>omegaMinusTauDict = getProbability(tauSumNumberList, tauDiscreteList, tauLoopNumber, omegaNodeNumber);
	map<double, double>muBeitaDict;
	double sum = 0;
	double sum1 = 0;
	for (int i = 0; i < tauLoopNumber.size(); i++)
	{
		sum1 = parentNode[tauLoopNumber[i]] * tauDiscreteList[i];
	}


	for (map<vector<double>, int>::iterator it = omegaMinusTauDict.begin(); it != omegaMinusTauDict.end(); it++)
	{
		sum = sum1;
		for (int j = 0; j < (it->first).size(); j++){
			sum += (it->first)[j] * parentNode[omegaMinusTao[j]];
		}
	
		muBeitaDict[sum] = (double)((double)omegaMinusTauDict[it->first] / (double)tauSumNumberList.size());
	}
	
	double cvalue = xZero(muBeitaDict, sigma2[k]);

	tauCList[tauLoopNumber] = cvalue;
}


vector<int> getSumSample(const vector<double>&tauDiscreteList, const vector<int>&tauLoopNumber)
{
	
	int row = DiscreX.size() - 2;
	set<vector<double>>tauSumSet;               //  不同顺序的vector<double>不知道算不算不同的成员
	tauSumSet.insert(tauDiscreteList);
	vector<int>tauSumNumberList;
	for (int i = 0; i < row; i++)
	{
		vector<double>tempList;
		for (int j = 0; j < tauLoopNumber.size(); j++)
		{
			tempList.push_back(DiscreX[i][tauLoopNumber[j]]);
		}
		
		if (tauSumSet.count(tempList) != 0)
		{
			
			tauSumNumberList.push_back(i);
		}
		
	}

	return tauSumNumberList;

}



/*
vector<int> getSumSample(const vector<double>&tauDiscreteList, const vector<int>&tauLoopNumber)
{

	int row = DiscreX.size() - 2;
	for (int i = 0; i < tauDiscreteList.size(); i++)
	{
		cout << tauDiscreteList[i] << ' ';
	}
	cout << endl;
	vector<int>tauSumNumberList;
	for (int i = 0; i < row; i++)
	{
		vector<double>tempList;
		for (int j = 0; j < tauLoopNumber.size(); j++)
		{
			tempList.push_back(DiscreX[i][tauLoopNumber[j]]);
		}

		if (tempList==tauDiscreteList)
		{
			cout << "已经压入***************************" << endl;
			tauSumNumberList.push_back(i);
		}
		else { cout << "meiyou fdfsdg*******************" << endl; }
	}

	return tauSumNumberList;

}
*/
vector<int> getOmegaMinusTao(const vector<int>&tauLoopNumber, const vector<int>&omegaNodeNumber)     //没写好
{
	vector<int>omegaNodeNumberTmp;
	for (int i = 0; i < omegaNodeNumber.size(); i++)
	{
		omegaNodeNumberTmp.push_back(omegaNodeNumber[i]);
	}

	for (int j = 0; j < tauLoopNumber.size(); j++)
	{
		vector<int>::iterator it = find(omegaNodeNumberTmp.begin(), omegaNodeNumberTmp.end(), tauLoopNumber[j]);
		if (it != omegaNodeNumberTmp.end())              //  vector的查找  找到返回迭代器对应位置，没有的话返回end()
		{
			omegaNodeNumberTmp.erase(it);      //  tauLoopNumber在omegaNodeNumberTmp中的话则删除
		}
	}
	return omegaNodeNumberTmp;
}



double Evalf(map<double,double>&muBeitaDict,double sigma2,double x)
{
	
	double f = 0;
	for (map<double, double>::iterator it = muBeitaDict.begin(); it != muBeitaDict.end(); it++)
	{
		f = f - (x - float(it->first)) / ((sqrt(2 * PI*sigma2))*sqrt(sigma2))*exp(-pow(x - float(it->first), 2) / (2 * sigma2))*muBeitaDict[it->first];    //   it->first==key   it->second==value
	}
	return f;
}



double xZero(map<double, double>&muBeitaDict, double sigma2)
{
	
	double f = 0.0;


	vector<double>keyList;

	for (map<double, double>::iterator it = muBeitaDict.begin(); it != muBeitaDict.end(); it++)
	{
		
		keyList.push_back(it->first);
	}

	if (keyList.size() == 0) {
		 return 0.0;
	}

	double maxmiu = keyList[0];
	double minmiu = keyList[0];
	for (int k = 0; k < keyList.size(); k++){
		if (keyList[k] > maxmiu) maxmiu = keyList[k];
		if (keyList[k] < minmiu) minmiu = keyList[k];
	}

	double midmiu = 0.0;

	
	double high = Evalf(muBeitaDict, sigma2, maxmiu);
	double low = Evalf(muBeitaDict, sigma2, minmiu);

	


	int i = 0;
	while (high*low < 0)
	{
		midmiu = (maxmiu + minmiu) / 2;
		if (abs(midmiu - maxmiu) < 0.000000001)
			break;

		double middle = Evalf(muBeitaDict,sigma2,midmiu);
		if (low*middle < 0)
		{
			maxmiu = midmiu;
			high = Evalf(muBeitaDict, sigma2, maxmiu);
			i = i + 1;
		}
		else if (middle*high < 0)
		{
			minmiu = midmiu;
			low = Evalf(muBeitaDict, sigma2, minmiu);
			i = i + 1;
		}
	}

	return midmiu;
}


double Discrete(double number, int col)
{
	int row = DiscreX.size();
	for (int x = 0; x < 10; x++)
	{
		if (DiscreX[row - 1][col] + (x + 1)*(DiscreX[row - 2][col] - DiscreX[row - 1][col]) / 10 >= number)
		{
			number = x + 1;
			break;
		}
	}
	return number;
}

map<vector<double>, int>getProbability(const vector<int>&tauSumNumberList, const vector<double>&tauDiscreteList, const vector<int>&tauLoopNumber, const vector<int>&omegaNodeNumber)
{
	vector<int>omegaMinusTao = getOmegaMinusTao(tauLoopNumber, omegaNodeNumber);
	vector<vector<double>>probabilityData = getProbabilityData(tauSumNumberList);

	set<vector<double>>key;
	map<vector<double>, int>dictory;
	vector<vector<double>>tmplist1;
	for (int i = 0; i < probabilityData.size(); i++)
	{
		vector<double>tmplist2;
		for (int j = 0; j < omegaMinusTao.size(); j++)
		{
			tmplist2.push_back(probabilityData[i][j]);
		}
		tmplist1.push_back(tmplist2);

		if (key.count(tmplist2) == 0){     //not in key
			dictory[tmplist2] = 1;
			key.insert(tmplist2);
		}
		else
		{
			dictory[tmplist2] = dictory[tmplist2] + 1;
		}
	}
	
	return dictory;
}


vector<vector<double>>getProbabilityData(const vector<int>&tauSumNumberList)
{
	
	vector<vector<double>>probabilityData;
	for (int i = 0; i < tauSumNumberList.size(); i++)
	{
		probabilityData.push_back(DiscreX[tauSumNumberList[i]]);      // DiscreX是离散化后的样本+待检测样本+最大最小值，应该就是对应python的dataListNew
	}
	return probabilityData;
}


void showPath()     //打印path
{
	map<int, vector<int>>::iterator it = path.begin();
	for (; it != path.end(); it++)
	{
		cout << it->first << ':';
		for (int i = 0; i < it->second.size(); i++)
		{
			cout << it->second[i] << ' ';
		}
		cout << endl;
	}
}


vector<double> getSigma2()     //新建得到各个结点的Sigma2
{
	double temp = 0.0;
	for (int i = 0; i < dim; i++)
	{
		temp = Parameter(i);
		sigma2.push_back(temp);
	}

	cout << "the sigma of every node is "<< endl<<'[';
	for (int i = 0; i < sigma2.size(); i++)
	{
		cout << sigma2[i] << ' ';
	}
	cout << ']';

	return sigma2;
}














































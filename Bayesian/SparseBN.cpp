#include"SparseBN.h"
#include"utils.h"
#include"Reason.h"
using namespace std;

deque<vector<double>> rawX;
deque<vector<double>> X;
vector<vector<double>> B;
vector<vector<double>> showB;
vector<vector<double>> DAG;
vector<vector<double>> P;
double lamda1;
long double lamda2;
int dim;
int num;
double relationThreshold;


void SBN(char* filename, int d, double penalty, double relation){

	dim = d;
	relationThreshold = relation;
	rawX = loadtxt(filename, dim);
	
	cout << "standard......................" << endl;
	cout << "样本大小为...................." << rawX.size() << endl;
	X = rawX;
	Scale(rawX, X);   //standard   0 mean and 1 variance



	num = X.size();
	B = Zero(dim - 1, dim);     //(p-1)*p 的
	showB = Zero(dim - 1, dim);         //用来显示 DAG 用的 B，在原始B上淘汰关系阈值以下的，并进行相应剪枝，总之 DAG 和这个 showB 是对应的，后面的推理也是基于这个showB， 不直接对B剪枝是消除B在增量更新时候的误差累计
	lamda1 = penalty*(0.11*num);    			//penalty
	//lamda2 = (((num - 1)*(num - 1)*dim) / lamda1 - lamda1)*10;
	
	
	lamda2 = 1000 * num;
	
	
	cout << num << ' ' << dim << ' ' << lamda1 <<"得到的lamda2 "<<lamda2<<endl;
	//vector<double> lamdaset=[Linspace(penalty,(dim*(num-1)*(num-1))/lamda1-lamda1,10)];   //set  lamda2 from lamda1 to (n-1)2 *p/lamda1 - lamda1   （备用，暂时没用上，两个lamda作为后续参数调节）

	DAG = Zero(dim, dim);     //DAG
	P = Zero(dim, dim);
}


vector<double> BFS(int i){      //   速度已检验
	//得到DAG中i为起点的闭包 
     
	for (int j = 0; j < dim; j++)
		P[i][j] = 0;    //i对应的行清空，重新BFS生成
	queue<double> q;
	q.push(i);
	while (!q.empty())
	{
		double node = q.front();
		q.pop();
		if (P[i][node] == 0)
		{
			P[i][node] = 1;
			for (int adj = 0; adj < dim; adj++)
			{
				if (DAG[node][adj] == 1 && P[i][adj] == 0)
					q.push(adj);
			}
		}
	}
	
	
	return P[i];    //更新P矩阵一列，确保闭包
}


void refreshDAG(int i){                       //BCD 对50维每一维迭代的时候都刷新对应行的DAG，这样BFS才有效
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

vector<double> getBij(int i, int j)  // get Bi / j    dim is(p - 2) * 1
{
	vector<double>result;
	result.clear();
	RemoveElem(B, i, j,result);
	return result;
}



void getXij(deque<vector<double>>&result, int dim1, int dim2)      //我们最多只需要移除两列         这个速度比python快，对于五万数据  这个6ms python14ms
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


void Shooting(int i, int shootIter)
{
	vector<double> xi;
	
	GetCol_S(X, i,xi);

	
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
			
			

			vector<double>xb;  Dot_vv(X, getBij(i, j), xb, i, j);
			
			
			

			for (unsigned k = 0; k < xi.size(); k++)
				tmp.push_back(xi[k] - xb[k]);			 //(n*1) - (n*48).*(48*1)    得到的分子 tmp shape is （n*1）


			vector<double> xj;
			GetCol_S(X, j,xj);     					//xj  (n*1)
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
				part2 = (lamda1 + lamda2*P[i][j]) / num;
			}
			else{
				part2 = (lamda1 + lamda2*P[i][j]) / Denominator;
			}


			//cout << "两个part " << part1 << "    " << part2 <<"  此时的Pij  "<<P[i][j]<<endl;

			if (P[i][j] == 1 && abs(part1) > part2)  cout << "*************出现了避环惩罚项不够的情况！！！！输入任意数字继续运行********************" << endl;

			double sign = Sign(part1);

			//cout << "part1 part2 " << part1 << ' ' << part2 << endl;

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
		std::cout << "shooting iter " << i << "，Bi->: " << newtemp << endl;
		
		if (fabs(newtemp - temp) < 0.01)      //如果差值小于0.01 则提前退出迭代，最多100次迭代
			break;
	}

}


void BCD(int BCDIter)
{
	for (int t = 0; t < BCDIter; t++)
	{
		cout << "进入下一个BCD" << endl;
		cout << "进入第" << t << "轮BCD迭代................." << endl;
		for (int i = 0; i < dim; i++)		 //反复迭代更新 Bi
		{
			int64_t start, end;
			start = GetTime();
			vector<double> Pij = BFS(i);
			cout << "第" << t << "轮BCD的第" << i << "维迭代............." << endl;
			Shooting(i, 100);
			refreshDAG(i);    //更新DAG 中 i  对应的那一列

			int sss = 0;
			for (int k = 0; k < dim; k++)
			{
				if (DAG[k][i] == 1) sss++;
			}
			cout << "该列现在有非零值" << sss << "个" << endl;

			for (int k = 0; k < dim - 1;k++)
			{
				cout << B[k][i] << ',';
			}
			cout << endl << endl << endl;


			end = GetTime();
			cout << (end - start)/1000 << endl;
		}
	}
}




void CheckAndPrune()  // 对 B 剪枝，避免成环    对于一个比较合适的阈值设定，可以保证不会成环，该函数是以防万一           或者 删除关系阈值以下的边，减少稀疏度，但也一定程度削减了精度
{
	//vector<vector<double>> showB;
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
			if (fabs(showB[i][j]) < relationThreshold)
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



void Train()	 // 用 BCD 训练新的 B, 检查环和关系阈值剪枝，更新二维矩阵 DAG
{
	time_t start, end;
	time(&start);
	BCD(10);  // 训练过程中不对 B 和 DAG过多干预，仅仅维持同步，训练好了结构以后开始检查环并剪枝
	time(&end);
	cout << "BCD结束" << endl;
	CheckAndPrune();
	cout << "剪枝结束" << endl;
	cout << "训练时间..................................................." << difftime(end, start) << "秒" << endl;
}



void update(vector<double> sample)     //可增量化
{
	rawX.pop_front();
	rawX.push_back(sample);   // 将sample封装成一行，即将rawX第一行删了再补新样本到最后一行
	Scale(rawX,X);

	printf("try to update..........\n");

	BCD(1);				// 增量以后更新DAG
	CheckAndPrune();   //对新的 B 进行检查，并根据关系阈值等参数更新 B 和 DAG

	printf("update end...............\n");
}

/****************************************打印***************************************************/
void printDAG()  // 输出 DAG
{
	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			if (DAG[i][j] == 1)			// print('存在边 ', i + 1, j + 1)
				cout << "G[" << i+1 << "]" << "[" << j+1 << "]" << "= 1" << endl;
			//print('G[', i, ',', j, ']=1');
		}
	}
}

void printB()
{
	for (int i = 0; i < dim - 1; i++)
	{
		for (int j = 0; j < dim; j++)
			cout << B[i][j] << " , ";
		cout << endl;
	}
}

void printShowB()
{
	for (int i = 0; i < dim - 1; i++)
	{
		for (int j = 0; j < dim; j++)
		
			cout << showB[i][j] << " , ";
		cout << endl;
	}
}


double Expection(int x)   // 参数 x 的期望
{
	vector<double>Xx;
	GetCol_S(X, x,Xx);
	double total = 0;
	for (unsigned i = 0; i < Xx.size(); i++)
		total += Xx[i];  // 累加第二列
	return total / num;
}

double Cov(int x, int y)
{
	vector<double>mul;
	vector<double>Xx, Xy;
	GetCol_S(X, x,Xx);
	GetCol_S(X, y,Xy);
	for (unsigned i = 0; i < Xy.size(); i++)			//点乘x,y列
		mul.push_back(Xx[i] * Xy[i]);
	//mul=self.X[:,x]*self.X[:,y]
	double total = 0;
	for (unsigned i = 0; i < mul.size(); i++)		//mul求和
		total += mul[i];
	total = total / num;   //  计算 E[XY]
	double covXY = total - Expection(x)*Expection(y);
	return covXY;
}

double Parameter(int x){         //计算 P(X| u1 u2 u3 u4.......) 的线性高斯所对应的方差
	double temp = 0;
	int posi = 0, posj = 0;
	for (int i = 0; i < dim - 1; i++)
		for (int j = 0; j < dim - 1; j++){
			if (i >= x) posi = i + 1;
			else posi = i;

			if (j >= x) posj = j + 1;
			else posj = j;
			temp += showB[i][x] * showB[j][x] * Cov(posi, posj);
		}
	double sigma2 = Cov(x, x) - temp;
	return sigma2;
}


/*
int main(){
	
	
	vector<int>aa,bb;
	aa.push_back(3);
	aa.push_back(4);
	aa.push_back(7);

	vector<int>::iterator it = find(aa.begin(), aa.end(), 4);
	if (it != aa.end()){
		aa.erase(it);
		cout << "yeah" << endl;
		for (int k = 0; k < aa.size(); k++) cout << aa[k] << ' ';
	}
	else cout << "fuck" << endl;




	/*
	SBN("D:\\data8000.csv",62,1.0,0.02);
	train();

	printDAG();
	cout << "--------------------------------------------" << endl;
	printB();
	cout << "********************************************" << endl;
    printShowB();
	*/
	/*
	int num = 0;
	for (int i = 0; i < sbn.B.size(); i++)
		for (int j = 0; j < sbn.B[0].size(); j++)
		{
			if (sbn.B[i][j] != 0) num++;
		}
	cout << "Mat B nonzeros " << num << endl;

	num = 0;
	for (int i = 0; i < sbn.showB.size(); i++)
		for (int j = 0; j < sbn.showB[0].size(); j++)
		{
			if (sbn.showB[i][j] != 0) num++;
		}
	cout << "Mat showB nonzeros " << num << endl;

	num = 0;
	for (int i = 0; i < sbn.DAG.size(); i++)
		for (int j = 0; j < sbn.DAG[0].size(); j++)
		{
			if (sbn.DAG[i][j] != 0) num++;
		}
	cout << "Mat DAG nonzeros " << num << endl;
	*/
//	return 0;


// VIXIndex.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//	//---------------
		//计算VixIndex program
	//time 20220305
	//author zhangting 

	//	/----------------

#include <iostream>
#include <fstream>
#include<vector>
#include <map>
#include <string>
#include<sstream>
#include <stdio.h>
#include<math.h>

using namespace std;

typedef struct  OptionDataType {
	double strick;
	double put_bid;
	double put_ask;
	double  call_bid;
	double call_ask;
}OptionDataType;

typedef struct  DataType {
	double strick; //行权价格
	string type;  //买卖方向 call put
	double mix;     //Pk(i)
	double variance;//方差
}DataType;


#define None   -999999
double rates[2] = { 0.000305,0.000286 };
double minutesOfTotalYear = (60 * 24 * 365);
int Nt[2] = { 855 + 510 + 34560,   854 + 900 + 44640 };  //近月合约生育到期时间 （单位分钟）

double T1 = Nt[0] / minutesOfTotalYear;
double T2 = (Nt[1]) / minutesOfTotalYear;
double T[2] = { T1,T2 };  //NT/N365

long GetSpecialChaNum(const char *_str, char ch) {
	long _retNum = 0;
	while (*_str) {
		if (*_str == ch) {
			_retNum++;
		}
		_str++;
	}
	return _retNum;
}

void GetOptionData(std::vector<string> &filepathVec, std::map<int, std::vector<OptionDataType>> & optionDatas) {

	for (int i = 0; i < 2; i++) {
		ifstream fp(filepathVec[i]);
		string line;
		bool isJumbHead = true;
		while (getline(fp, line)) {
			if (isJumbHead) {
				isJumbHead = false;
				continue;
			}
			std::vector<double> rowData;
			string  number;
			istringstream readstr(line);
			int delimiterNum = GetSpecialChaNum(line.c_str(), ',');
			for (int j = 0; j < delimiterNum + 1; j++) {
				getline(readstr, number, ',');
				rowData.push_back(std::stod(number));
			}
			OptionDataType optiondatas;
			optiondatas.strick = rowData[0];
			optiondatas.put_bid = rowData[1];
			optiondatas.put_ask = rowData[2];
			optiondatas.call_bid = rowData[3];
			optiondatas.call_ask = rowData[4];
			optionDatas[i].push_back(optiondatas);
		}
	}
}

//获取近远期价格水平值
void GetNearAndFarFutureValue(std::map<int, std::vector<OptionDataType>> & optionDatasMap, std::vector<double> &resultVec)
{
	for (int i = 0; i < 2; i++) {
		double mindiff = None;
		double diff = None;
		double fstrike = None;//行权价格
		double fcall = None;
		double fput = None;
		for (auto &col : optionDatasMap[i]) {
			diff = std::abs(((col.put_ask + col.put_bid) / 2) - ((col.call_ask + col.call_bid) / 2));
			if (mindiff == None || diff < mindiff) {
				mindiff = diff;
				fstrike = col.strick;
				fcall = (col.put_ask + col.put_bid) / 2;
				fput = (col.call_ask + col.call_bid) / 2;
			}
		}
		double f = fstrike + exp(rates[0] * T[0]) *(fcall - fput);
		resultVec.push_back(f);
	}


}

void GetValueOfK0(std::vector<double> fValueVector, map<int, std::vector<DataType>> &selectOptVecs, map<int, std::vector<OptionDataType>> & optionDatas, std::vector<double> &resultVec) {
	for (int k = 0; k < 2; k++)
	{
		int i = 0;
		int k0i = 0;
		for (auto &col : optionDatas[k]) {
			if (col.strick < fValueVector[0]) { //找到接近F值得位置  确定好k0的值
				resultVec[k] = col.strick;
				k0i = i;
				i += 1;
			}
		}
		OptionDataType d = optionDatas[k][k0i];
		DataType data;
		data.strick = d.strick;
		data.type = "put /call average"; //put /call average
		data.mix = (((d.put_ask + d.put_bid) / 2) + ((d.call_ask + d.call_bid) / 2)) / 2;
		std::vector<double> selectOptVec;
		selectOptVecs[k].push_back(data);

		i = k0i - 1;  //向上移动一个位置
		bool b = true;
		double previousbid = None;
		while (b  && i >= 0) {
			OptionDataType d = optionDatas[k][i];  //获取上一行的数据
			if (d.call_bid > 0) {
				DataType data;
				data.strick = d.strick;
				data.type = "put "; //put /call average
				data.mix = (d.call_ask + d.call_bid) / 2;
				selectOptVecs[k].insert(selectOptVecs[k].begin(), data); //按照顺序插入
			}
			else {
				if (previousbid == 0) b = false;
			}
			previousbid = d.call_bid;
			i -= 1;
		}
		//=============
		i = k0i + 1;  //向下移动一个位置
		b = true;
		previousbid = None;
		while (b  && i < optionDatas[k].size()) {
			OptionDataType d = optionDatas[k][i];  //获取上一行的数据
			if (d.put_bid > 0) {
				DataType data;
				data.strick = d.strick;
				data.type = "call "; //put /call average
				data.mix = (d.put_ask + d.put_bid) / 2;
				selectOptVecs[k].push_back(data);
			}
			else {
				if (previousbid == 0) b = false;
			}
			previousbid = d.put_bid;
			i += 1;
		}
	}



}
//calculate volatility for both near-term and next-term and next-term
void GetVolatilityValue(map<int, std::vector<DataType>> &selectOptVecs)
{
	for (int k = 0; k < 2; k++)
	{
		int i = 0;
		for (auto &d : selectOptVecs[k]) {
			double deltak = 0;
			if (i == 0) {
				deltak = selectOptVecs[k][1].strick - selectOptVecs[k][0].strick;
			}
			else if (i == selectOptVecs[k].size() - 1) {
				deltak = selectOptVecs[k][i].strick - selectOptVecs[k][i - 1].strick;
			}
			else {
				deltak = (selectOptVecs[k][i + 1].strick - selectOptVecs[k][i - 1].strick) / 2;
			}
			double contributionbystrike = (deltak / (d.strick * d.strick)) *exp(rates[0] * T[0]) * d.mix;
			selectOptVecs[k][i].variance = contributionbystrike; //每个vector 元素都加一了
			i += 1;

		}
	}
}

double GetVarianceValue(map<int, std::vector<DataType>> &selectOptVecs, const std::vector<double> &F, const std::vector<double> &K0)
{
	if (selectOptVecs.size() != 2) {
		return -1;
	}
	double aggregatedcontributionbystrike1 = 0, aggregatedcontributionbystrike2 = 0;
	for (auto &d : selectOptVecs[0]) {
		aggregatedcontributionbystrike1 += d.variance;
	}
	aggregatedcontributionbystrike1 = (2 / T[0])*aggregatedcontributionbystrike1;

	for (auto &d : selectOptVecs[1]) {
		aggregatedcontributionbystrike2 += d.variance;
	}
	aggregatedcontributionbystrike2 = (2 / T[0])*aggregatedcontributionbystrike2;
	double sigmasquared1 = 0, sigmasquared2 = 0;
	sigmasquared1 = aggregatedcontributionbystrike1 - (1 / T[0])*(F[0] / K0[0] - 1)*(F[0] / K0[0] - 1);
	sigmasquared2 = aggregatedcontributionbystrike1 - (1 / T[1])*(F[1] / K0[1] - 1)*(F[1] / K0[1] - 1);

	//Step 6: Calculate the 30 - day weighted average of sigmasquared[0] and sigmasquared[1](p9)
	int N30 = 30 * 1440;
	int N365 = 365 * 1440;
	double VIX = 100 * sqrt((T[0] * sigmasquared1*(Nt[1] - N30) / (Nt[1] - Nt[0]) + T[1] * sigmasquared2*(N30 - Nt[0]) / (Nt[1] - Nt[0]))*N365 / N30);
	return VIX;

}
int main()
{

	printf("Nt1: %d  Nt2: %d \n", Nt[0], Nt[1]);
	printf("T1: %f  T2: %f \n", T[0], T[1]);

	std::map<int, std::vector<OptionDataType>>optonDatasMap;
	std::vector<string> filePathVec{ "near.csv","next.csv" };
	//Step 1 获取option近期 和后期数据
	GetOptionData(filePathVec, optonDatasMap);

	//Step2 获远期价格水平值 fVauleVec[0] fVauleVec[1]---->  F1  F2
	std::vector<double> F;

	GetNearAndFarFutureValue(optonDatasMap, F);
	if (F.size() == 2) {
		printf("F1: %f  F2: %f \n", F[0], F[1]);
	}
	//Step3 计算F这个数下面最近一档的行权价K0
	map<int, std::vector<DataType>> selectOptVecs;
	std::vector<double> K0;
	K0.resize(2);

	GetValueOfK0(F, selectOptVecs, optonDatasMap, K0);

	//#Step 4: Calculate volatility for both near - term and next - term options
	//求出了方差
	GetVolatilityValue(selectOptVecs);

	//Step5  求方差
	double VIX = GetVarianceValue(selectOptVecs, F, K0);
	std::cout << "VIX Value is :" << VIX << endl;

	exit(0);

}


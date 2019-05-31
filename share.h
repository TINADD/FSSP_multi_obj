#pragma once
#ifndef SHARE_H
#include<iostream>
#include<algorithm>
#include<vector>
#include<set>
#include<map>
#include<fstream>
#include<sstream>
#include<cmath>
#include<string>
#include<cstdlib>
#include<climits>
#include<ctime>
#include<cstdlib>
using namespace std;

class Known {
public:
	static map<string, string> pars_wl; //参数
	static int job_count, stage_count, machine_count;
	static vector<int> due_date;//工件的截止日期
	static vector<double> ep;//电价
	static vector<vector<int>> bpc, pc;
	static vector <vector <double>> bec;
	static int iter;  //当前迭代次数
	static unsigned long next; //随机种子


	static vector<pair<vector<bool>, vector<bool>>> hs_1;
	static vector<pair<vector<bool>, vector<bool>>> hs_2;

	static vector<vector<vector<map<int, int>>>> tabuTable; //基于动作tabuTable

	static double runtime; //程序运行时间
	static time_t st; //记录程序开始时间
	static time_t et;//记录程序结束时间
	static time_t arch_et;
	static int genes;
	static int pops;
	static int TS_iter;
	static int disturb_its;
	static double acp; //接受互不支配的概率
	static double tb_acp; //禁忌解存放互不支配解的概率
	static double w; //变异概率
	static double c1;
	static double c2;
	static double Pc;//遗传算法-交叉概率
	static double Pm;//遗传算法-变异概率
	static double itr_time;
	void read_file();
	void print_file();
	void write_head();
	int myrand(void);
	void mysrand(unsigned seed);
};
//--操作类
class Operation {
public:
	int jno;//curr_job_no
	int sno;//curr_stage_no
	int machine;//分配的机器
	int m_speed;//分配的速度--其实是选择的处理时间--m_speed = 对应的下标
	int po;//在当前机器上的位置
	int st;//操作的开始time
	int et;//结束time
	int process_time;//操作总的处理时间
	double energy_cost;//消耗的能源成本
	Operation();
	double update_energy_cost();//计算energy_cost
};

class TabuOp {
public:
	int jno = -1;//工件号
	int sno = -1;//工序号
	int mno = -1;//机器号
	int mv = -1;//机器速度下标
	int mpo = -1;//机器上的位置
};

class Idles {
public:
	//空闲时间段的开始-结束
	vector<pair<int, int>> idles; //idles[i].first = job[opjno].st,idles[i].second = job[opjno].et;
	vector<int> opjno;
	int jnum = 0;//操作数
	int ava_time = 0;//机器最早可用时间
};


#endif // !SHARE_H

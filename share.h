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
	static map<string, string> pars_wl; //����
	static int job_count, stage_count, machine_count;
	static vector<int> due_date;//�����Ľ�ֹ����
	static vector<double> ep;//���
	static vector<vector<int>> bpc, pc;
	static vector <vector <double>> bec;
	static int iter;  //��ǰ��������
	static unsigned long next; //�������


	static vector<pair<vector<bool>, vector<bool>>> hs_1;
	static vector<pair<vector<bool>, vector<bool>>> hs_2;

	static vector<vector<vector<map<int, int>>>> tabuTable; //���ڶ���tabuTable

	static double runtime; //��������ʱ��
	static time_t st; //��¼����ʼʱ��
	static time_t et;//��¼�������ʱ��
	static time_t arch_et;
	static int genes;
	static int pops;
	static int TS_iter;
	static int disturb_its;
	static double acp; //���ܻ���֧��ĸ���
	static double tb_acp; //���ɽ��Ż���֧���ĸ���
	static double w; //�������
	static double c1;
	static double c2;
	static double Pc;//�Ŵ��㷨-�������
	static double Pm;//�Ŵ��㷨-�������
	static double itr_time;
	void read_file();
	void print_file();
	void write_head();
	int myrand(void);
	void mysrand(unsigned seed);
};
//--������
class Operation {
public:
	int jno;//curr_job_no
	int sno;//curr_stage_no
	int machine;//����Ļ���
	int m_speed;//������ٶ�--��ʵ��ѡ��Ĵ���ʱ��--m_speed = ��Ӧ���±�
	int po;//�ڵ�ǰ�����ϵ�λ��
	int st;//�����Ŀ�ʼtime
	int et;//����time
	int process_time;//�����ܵĴ���ʱ��
	double energy_cost;//���ĵ���Դ�ɱ�
	Operation();
	double update_energy_cost();//����energy_cost
};

class TabuOp {
public:
	int jno = -1;//������
	int sno = -1;//�����
	int mno = -1;//������
	int mv = -1;//�����ٶ��±�
	int mpo = -1;//�����ϵ�λ��
};

class Idles {
public:
	//����ʱ��εĿ�ʼ-����
	vector<pair<int, int>> idles; //idles[i].first = job[opjno].st,idles[i].second = job[opjno].et;
	vector<int> opjno;
	int jnum = 0;//������
	int ava_time = 0;//�����������ʱ��
};


#endif // !SHARE_H

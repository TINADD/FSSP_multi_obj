#pragma once

#ifndef  NSGAII_H
#include "share.h"

/*start-������ҵ�������ͷ�ļ�*/
class Individual;
class Fac;
class Individual {
public:
	vector<vector<Operation>> sol;//�������
	vector<vector<vector<int>>> m_job;//���������Ӧ��--�����Ĵ�������ż���

	int tot_tardiness; //�ܵ��ӳ�ʱ��
	double tot_energy_cost;//�ܵ��ܺĳɱ�
	double cd;//ӵ��ϵ��

	Individual();
	void init_indi(int &rm);//�������һ������
	void _init_indi(int &rm); //�в��Ե�����һ������
	void update_sol(); //���㵱ǰ����Ŀ��ֵ
	void _update_sol();//�������ȼ���Ŀ��ֵ

	void update_tot_tardiness();//�ܵ��ӳ�ʱ��
	void update_tot_energy_cost();//�ܵ��ܺĳɱ�
	void update_tbb(TabuOp &tbo,int &pos); //���½��ɱ�

	void cal_machine_job();
	void print_sol(); //��ӡ�����
	void print_sol_();
	void LS_based_move(int &pos);//���ڶ����Ľ�������
	void local_search(int &pos); //�򵥾ֲ�����-first improved


	bool operator<(Individual &b) const;
	bool operator==(Individual &b) const;
	bool check();
};
class Fac
{
public:
	/*��Ա����*/
	static vector<Individual> populations;//��Ⱥ
	static vector<Individual> offSpring; //offSpring.size == pops; offSpring = GA(populations);	
	static vector<Individual> Rt;//Rt = offSpring + populations ;Rt.size == 2*pops;
	vector<Individual> archive; //������е���������populations��������Ž�
	/*��Ա����*/
	Fac();
	void start_up();//��������
	void start_down();//�������

	void init_pops();//���������Ⱥ
	void NSGAII(); //ʹ��NSGHAII�㷨

	vector<vector<int>> fast_non_dominated_sort();//��֧������

	void crowding_distance_assignment(vector<int> &non_domi);

	void GA(int &iter); //�Ŵ��㷨

	void find_opt();//����NSGAII+GA�㷨

	void writesol();// �ѽ�д�뵽�ļ���
	void push_into_arch(Individual &sol);

	void print_pops(vector<Individual> &popu);//���������Ⱥ
};



/*end-������ҵ�������*/

#endif // ! FSSP_H
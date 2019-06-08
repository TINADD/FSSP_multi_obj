#pragma once
#ifndef MOPSO_H
#define MOPSO_H
#include"share.h"

/*ȫ�ֱ���*/
class Solution;
class Operation;
class Factory;

// ��ǰ�����Ĳ�������
using ele = vector<Operation>;

//--���������
class Solution {
public:
	vector<vector<ele>> sol_1;//version_2

	int tot_tardiness;//Ŀ�꺯��1:���ӳ�
	double tot_energy_cost;//Ŀ�꺯��2:�ܵ��ܺĳɱ�
	double crowd_dis;//ӵ����
	Solution();
	int generate_sol_(int rm); //�������һ����
	int generate_sol(int rm); //�в��Ե�����һ����

	void update_sol();//�������Ȱ汾���ڼ����ʱ����ܻ���������ڹ����ϵ�˳��
	void update_sol_no();//���������Ȱ汾�������ڼ���Ŀ��ֵʱ���ı乤��˳��
	void calcu_tot_tardiness();//����Ŀ��ֵ1--���ӳ�
	void calcu_tot_energy_cost();//����Ŀ��ֵ2--���ܺ�

	void print_sol(); //��ӡ��Ľṹ

	bool operator<(const Solution &sol) const;
	bool operator==(const Solution &sol) const;
	bool check_sol(); //�жϽ��Ƿ���ȷ

};

static Solution gbest; //ȫ������

//--������
class Particle {
public:
	Solution  pbest; //�����ӵ���ʷ���Ž�
	Solution  sol;//��ǰ���ӽ������

	Particle();
	void print_(); //��ӡĿ��ֵ����ʱû��

	void tabu_search(int &index);
	void local_search();//�ֲ������ҵ��ȵ�ǰ���õ�-�ı��ڻ����ϵ�˳��
	void LS_bm(int &pos); //���ڶ����Ľ�������
	void LS_bm2(int &pos); //first improve ���ڶ����Ľ�������
	void update_tbb(TabuOp &tbo,int &pos); //����λ��pos�Ľ��ɱ�
	//bool is_in_tbb(int &pos); //�жϵ�ǰ�����Ƿ񱻽���

	void update_HT(Solution &sol, int index);//����hash��
	bool func_HT(Solution &sol, int index);//����hashֵ
	void update_position(int &in); //��ɢ
};

//--������
class Factory {
public:
	static vector<Particle> par_swarm;//����Ⱥ
	static vector<Solution> archive;//�浵��--������ʷ���������������н�

	Factory();
	void initial_sol();//���ɳ�ʼ��	

	void control();//���ڽ�Ľ������
	void control_based_move();//���ڶ����Ľ������
	void control_ls();//�򵥾ֲ��������

	bool is_into_archive(Solution &sol);
	void print_archive();
	void print_pops();//��ӡ����Ⱥ

	void random_par();//�������һ������

	void infile(); //д�ļ�
	void writesol();//�����Ž�д���ļ���

};
#endif
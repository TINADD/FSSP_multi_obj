#pragma once
#ifndef MOPSO_H
#define MOPSO_H
#include"share.h"

/*全局变量*/
class Solution;
class Operation;
class Factory;

// 当前机器的操作集合
using ele = vector<Operation>;

//--解决方案类
class Solution {
public:
	vector<vector<ele>> sol_1;//version_2

	int tot_tardiness;//目标函数1:总延迟
	double tot_energy_cost;//目标函数2:总的能耗成本
	int crowd_dis;//拥挤度

	Solution();
	int generate_sol_(int rm); //随机生成一个解
	int generate_sol(int rm); //有策略的生成一个解

	void update_sol();//主动调度版本（在计算的时候可能会调整操作在工件上的顺序）
	void update_sol_no();//非主动调度版本（即：在计算目标值时不改变工件顺序）
	void calcu_tot_tardiness();//计算目标值1--总延迟
	void calcu_tot_energy_cost();//计算目标值2--总能耗

	void print_sol(); //打印解的结构

	bool operator<(const Solution &sol) const;
	bool operator==(const Solution &sol) const;
	bool check_sol(); //判断解是否正确

};

static Solution gbest; //全局最优

//--粒子类
class Particle {
public:
	Solution  pbest; //该粒子的历史最优解
	Solution  sol;//当前粒子解决方案

	Particle();
	void print_(); //打印目标值，暂时没用

	void tabu_search(int &index);
	void local_search();//局部搜索找到比当前更好的-改变在机器上的顺序
	void LS_bm(int &pos); //基于动作的禁忌搜索
	void LS_bm2(int &pos); //first improve 基于动作的禁忌搜索
	void update_tbb(TabuOp &tbo,int &pos); //更新位置pos的禁忌表
	//bool is_in_tbb(int &pos); //判断当前动作是否被禁忌

	void update_HT(Solution &sol, int index);//更新hash表
	bool func_HT(Solution &sol, int index);//计算hash值
	void update_position(int &in); //发散
};

//--工厂类
class Factory {
public:
	static vector<Particle> par_swarm;//粒子群
	static vector<Solution> archive;//存档集--迭代历史上所有最优帕累托解

	Factory();
	void initial_sol();//生成初始解	

	void control();//基于解的禁忌入口
	void control_based_move();//基于动作的禁忌入口
	void control_ls();//简单局部搜索入口

	bool is_into_archive(Solution &sol);
	void print_archive();
	void print_pops();//打印粒子群

	void random_par();//随机生成一个粒子

	void infile(); //写文件
	void writesol();//将最优解写入文件里

};
#endif
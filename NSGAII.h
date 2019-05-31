#pragma once

#ifndef  NSGAII_H
#include "share.h"

/*start-柔性作业车间调度头文件*/
class Individual;
class Fac;
class Individual {
public:
	vector<vector<Operation>> sol;//解决方案
	vector<vector<vector<int>>> m_job;//解决方案对应的--机器的处理工件编号集合

	int tot_tardiness; //总的延迟时间
	double tot_energy_cost;//总的能耗成本
	double cd;//拥挤系数

	Individual();
	void init_indi(int &rm);//随机生成一个个体
	void _init_indi(int &rm); //有策略的生成一个个体
	void update_sol(); //计算当前个体目标值
	void _update_sol();//主动调度计算目标值

	void update_tot_tardiness();//总的延迟时间
	void update_tot_energy_cost();//总的能耗成本
	void update_tbb(TabuOp &tbo,int &pos); //更新禁忌表

	void cal_machine_job();
	void print_sol(); //打印个体解
	void print_sol_();
	void LS_based_move(int &pos);//基于动作的禁忌搜素
	void local_search(int &pos); //简单局部搜索-first improved


	bool operator<(Individual &b) const;
	bool operator==(Individual &b) const;
	bool check();
};
class Fac
{
public:
	/*成员变量*/
	static vector<Individual> populations;//种群
	static vector<Individual> offSpring; //offSpring.size == pops; offSpring = GA(populations);	
	static vector<Individual> Rt;//Rt = offSpring + populations ;Rt.size == 2*pops;
	vector<Individual> archive; //存放所有迭代过程中populations里面的最优解
	/*成员函数*/
	Fac();
	void start_up();//程序启动
	void start_down();//程序结束

	void init_pops();//随机生成种群
	void NSGAII(); //使用NSGHAII算法

	vector<vector<int>> fast_non_dominated_sort();//非支配排序

	void crowding_distance_assignment(vector<int> &non_domi);

	void GA(int &iter); //遗传算法

	void find_opt();//调用NSGAII+GA算法

	void writesol();// 把解写入到文件中
	void push_into_arch(Individual &sol);

	void print_pops(vector<Individual> &popu);//输出父代种群
};



/*end-柔性作业车间调度*/

#endif // ! FSSP_H
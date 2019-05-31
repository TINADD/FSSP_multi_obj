#include "MOPSO.h"

vector<Particle> Factory::par_swarm;//粒子群
vector<Solution> Factory::archive;//存档集--迭代历史上所有最优帕累托解
Known known = Known();

/*基于基本处理时间排序--stage one*/
bool cmp_bpt(pair<int, int> &m1, pair<int, int> &m2)
{
	//m1.first = jno
	return known.bpc[0][m1.first] < known.bpc[0][m2.first];
}


/*基于上一道工序的结束时间排序*/
bool cmp_et(pair<int, int> &m1, pair<int, int> &m2)
{
	//m1.second = 当前工件上一道工序的结束时间
	return m1.second < m2.second;
}


//待看
bool cmp_machine(pair<int, int> &m1, pair<int, int> &m2)
{
	return  m1.second < m2.second;
}

//解决方案类函数实现
Solution::Solution()
{
	sol_1 = vector<vector<ele>>(known.stage_count, vector<ele>(known.machine_count));
	tot_tardiness = 0;
	tot_energy_cost = 0.;
}

//生成初始解--基于策略
int Solution::generate_sol(int rm)
{
	vector<pair<int, int>> jobs(known.job_count);//jobs[i].first = 工件号；jobs[i].second = 上一道工序完成时间
	int mac = 0, speed = 0, v = 0;
	for (int i = 0; i < jobs.size(); ++i)
	{
		jobs[i].first = i;
		jobs[i].second = 0;
	}
	vector<Idles> midles(known.machine_count);

	sort(jobs.begin(), jobs.end(), cmp_bpt); //sort based basement processing time

	for (int j = 0; j < jobs.size(); ++j)
	{
		Operation tmp = Operation();
		tmp.sno = 0; //工序号
		tmp.jno = jobs[j].first; //工件号

		srand(rand()%100000 + rm);
		mac = rand() % known.machine_count;
		//cout << "m" << mac << '\t';
		srand(rand()%100000 + rm);
		speed = rand() % known.pc[0].size();//处理时间对应的下标
		//cout << speed << '\t';
		v = known.pc[0][speed];//处理时间

		tmp.machine = mac;
		tmp.po = midles[mac].jnum;// 在当前机器上的位置，从0开始算
		++midles[mac].jnum;
		tmp.m_speed = speed;
		tmp.st = midles[mac].ava_time;
		tmp.process_time = known.bpc[0][tmp.jno] + v;
		midles[mac].ava_time = jobs[j].second = tmp.et = (tmp.st + tmp.process_time);

		tmp.energy_cost = tmp.update_energy_cost();
		sol_1[0][mac].push_back(tmp);
	}

	for (int s = 1; s < known.stage_count; ++s)
	{
		midles = vector<Idles>(known.machine_count);
		sort(jobs.begin(), jobs.end(), cmp_et); //sort based pre_stage et
		//sort(jobs.begin(),jobs.end(),cmp_machine);
		for (int j = 0; j < jobs.size(); ++j)
		{
			Operation tmp = Operation();
			tmp.sno = s;
			tmp.jno = jobs[j].first;
			srand(rand()%1000000 + rm);
			mac = rand() % known.machine_count;
			//cout << "m" << mac << '\t';

			srand(rand()%1000000+rm);
			speed = rand() % known.pc[s].size();//处理时间对应的下标
			//cout << speed << '\t';

			v = known.pc[s][speed]; //处理时间

			tmp.machine = mac;
			tmp.po = midles[mac].jnum;
			++midles[mac].jnum;
			tmp.m_speed = speed;
			tmp.st = max(jobs[j].second, midles[mac].ava_time);
			tmp.process_time = known.bpc[s][tmp.jno] + v;
			midles[mac].ava_time = jobs[j].second = tmp.et = tmp.st + tmp.process_time;
			tmp.energy_cost = tmp.update_energy_cost();

			sol_1[s][mac].push_back(tmp);
		}
	}
	calcu_tot_tardiness();
	calcu_tot_energy_cost();
	return 1;
}

int Solution::generate_sol_(int rm)
{
	for (int s = 0; s < known.stage_count; ++s)
	{
		for (int j = 0; j < known.job_count; ++j)
		{
			Operation oper = Operation();
			oper.jno = j;
			oper.sno = s;
			//随机生成一台机器&&一个处理时间
			//srand(clock() + rand() % 1000 + rm);
			oper.machine = rand() % known.machine_count;
			//srand(clock() + rm);
			oper.m_speed = rand() % known.pc[s].size();
			sol_1[s][oper.machine].push_back(oper);
		}
	}
	update_sol();
	return 0;
}

void Solution::print_sol()
{
	//打印每台机器的处理时间段序列
	cout << "机器的处理时间段序列" << endl;
	for (int s = 0; s < known.stage_count; ++s)
	{
		int js = 0;
		for (int m = 0; m < known.machine_count; ++m)
		{
			cout << "[ ";
			for (int op = 0; op < sol_1[s][m].size(); ++op)
			{
				++js;
				cout << sol_1[s][m][op].st << "," << sol_1[s][m][op].et << " ";

			}
			cout << "]	";
		}
		cout << "== " << js << endl;
	}

	cout << "工件的时间段序列：" << endl;
	//打印每个工件的时间段序列
	vector<vector<pair<int, int>>> tj(known.job_count, vector<pair<int, int>>(known.stage_count));
	for (int s = 0; s < known.stage_count; ++s)
	{
		for (int m = 0; m < known.machine_count; ++m)
		{
			ele &e = sol_1[s][m];
			for (int op = 0; op < e.size(); ++op)
			{
				tj[e[op].jno][s].first = e[op].st;
				tj[e[op].jno][s].second = e[op].et;
			}
		}

	}
	for (int j = 0; j < known.job_count; ++j)
	{
		cout << "job " << j << ":  ";
		for (int s = 0; s < known.stage_count; ++s)
		{
			cout << tj[j][s].first << "," << tj[j][s].second << "  ";
		}
		cout << '\n';
	}

	cout << "目标值： " << tot_tardiness << "	" << tot_energy_cost << endl;
}


/*重载<运算符*/
bool Solution::operator<(const Solution &sol) const
{
	if ((tot_tardiness < sol.tot_tardiness) && (tot_energy_cost <= sol.tot_energy_cost) ||
		(tot_tardiness < sol.tot_tardiness) && (tot_energy_cost < sol.tot_energy_cost)) return true;
	return false;
}
/*重载==运算符*/
bool Solution::operator==(const Solution &sol) const
{
	if (((tot_tardiness < sol.tot_tardiness) && (tot_energy_cost > sol.tot_energy_cost))
		|| ((tot_tardiness > sol.tot_tardiness) && (tot_energy_cost < sol.tot_energy_cost))) return true;
	return false;
}

bool Solution::check_sol()
{
	bool flag = true;
	vector<int> pre_et(known.job_count, 0);
	for (int s = 0; s < known.stage_count; ++s)
	{
		set<int> ops; //判断是否每个工件都被处理了
		int js = 0; //判断当前工序所有机器上的操作总数是否等于job_count
		for (int m = 0; m < known.machine_count; ++m)
		{
			vector<int> seq; //判断当前机器是否在某个时刻只处理一个工件
			for (int op = 0; op < sol_1[s][m].size(); ++op)
			{
				++js;
				ops.insert(sol_1[s][m][op].jno);
				seq.push_back(sol_1[s][m][op].st);
				seq.push_back(sol_1[s][m][op].et);
				if (pre_et[sol_1[s][m][op].jno] > sol_1[s][m][op].st) return false;
				pre_et[sol_1[s][m][op].jno] = sol_1[s][m][op].et;
			}

			for (int si = 2; si < ((int)seq.size() - 1); si += 2)
			{
				if (seq[si] < seq[si - 1]) return false;
			}
		}
		if (js != known.job_count || ops.size() != known.job_count) return false;
	}
	return flag;
}

/*计算总延迟*/
void Solution::calcu_tot_tardiness()
{
	tot_tardiness = 0;
	for (int m = 0; m < known.machine_count; ++m)
	{
		for (int op = 0; op < sol_1[known.stage_count - 1][m].size(); ++op)
		{
			int job = sol_1[known.stage_count - 1][m][op].jno;
			int et = sol_1[known.stage_count - 1][m][op].et;
			tot_tardiness += max(0, et - known.due_date[job]);
		}
	}

}

/*计算总能耗*/
void Solution::calcu_tot_energy_cost()
{
	tot_energy_cost = 0.;
	for (int s = 0; s < known.stage_count; ++s)
	{
		for (int m = 0; m < known.machine_count; ++m)
		{
			for (int op = 0; op < sol_1[s][m].size(); ++op)
			{
				tot_energy_cost += sol_1[s][m][op].energy_cost;
			}
		}
	}
}

void Solution::update_sol()
{   //调整顺序的计算
	vector<int> jobs(known.job_count, 0);//jobs[i] = 上一道工序完成时间 
	double cmp1 = 0., cmp2 = 0.;
	//stage 1
	for (int s = 0; s < known.stage_count; ++s)
	{
		vector<Idles> midles(known.machine_count);//当前阶段每台机器的空闲时间段集合
		for (int m = 0; m < known.machine_count; ++m)
		{
			ele &ele_ = sol_1[s][m]; //当前机器的操作集合
			for (int op = 0; op < ele_.size(); ++op)
			{
				ele_[op].process_time = known.bpc[s][ele_[op].jno] + known.pc[s][ele_[op].m_speed];
				if (!midles[m].ava_time)
				{
					//是当前机器处理的第一个操作
					ele_[op].st = max(midles[m].ava_time, jobs[ele_[op].jno]);
					jobs[ele_[op].jno] = midles[m].ava_time = ele_[op].et = ele_[op].st + ele_[op].process_time;
					ele_[op].energy_cost = ele_[op].update_energy_cost();
				}
				else
				{
					pair<int, int> prep;
					bool flag = false;  //判断当前操作是否插入了前面空闲时间
					//没有可用空闲段
					ele_[op].st = max(midles[m].ava_time, jobs[ele_[op].jno]);
					ele_[op].et = ele_[op].st + ele_[op].process_time;
					cmp2 = ele_[op].energy_cost = ele_[op].update_energy_cost(); //正常顺序的能耗成本

					for (int i = 0; i < midles[m].idles.size(); ++i)
					{//遍历机器的空闲时间段

						if (max(midles[m].idles[i].first, jobs[ele_[op].jno]) + ele_[op].process_time <= midles[m].idles[i].second)
						{//可以插入，但还要判断插入会不会使能耗变小

							Operation tmpop = ele_[op];
							tmpop.st = max(midles[m].idles[i].first, jobs[ele_[op].jno]);
							tmpop.et = tmpop.st + ele_[op].process_time;
							cmp1 = tmpop.update_energy_cost();
							if (cmp1 <= cmp2)
							{
								//可以插入
								prep.first = midles[m].idles[i].first;
								prep.second = tmpop.st;
								ele_[op] = tmpop;
								ele_[op].energy_cost = cmp1;

								jobs[ele_[op].jno] = midles[m].idles[i].first = ele_[op].et;

								//改变当前操作在机器中的处理顺序
								ele_.insert(ele_.begin() + i + 1, ele_[op]);
								ele_.erase(ele_.begin() + op + 1);

								//增加空闲时间段
								midles[m].idles.insert(midles[m].idles.begin() + i, prep);
								flag = true;
								break;
							}
							else continue;
						}
					}
					if (!flag)
					{
						//没有可用空闲段
						prep.first = midles[m].ava_time;
						prep.second = ele_[op].st;

						jobs[ele_[op].jno] = midles[m].ava_time = ele_[op].et;
						midles[m].idles.push_back(prep);
					}
				}
			}
			for (int i = 0; i < ele_.size(); ++i)
			{
				ele_[i].po = i;
			}
		}
	}
	calcu_tot_tardiness();
	calcu_tot_energy_cost();
	return;
}

void Solution::update_sol_no()
{
	vector<int> jobs(known.job_count, 0);//jobs[i] = 上一道工序完成时间 
	for (int s = 0; s < known.stage_count; ++s)
	{
		vector<Idles> midles(known.machine_count);//当前阶段每台机器的空闲时间段集合
		for (int m = 0; m < known.machine_count; ++m)
		{
			ele &ele_ = sol_1[s][m]; //当前机器的操作集合
			for (int op = 0; op < ele_.size(); ++op)
			{
				ele_[op].po = op + 1;
				ele_[op].process_time = known.bpc[s][ele_[op].jno] + known.pc[s][ele_[op].m_speed];
				ele_[op].st = max(midles[m].ava_time, jobs[ele_[op].jno]);
				midles[m].ava_time = ele_[op].et = ele_[op].st + ele_[op].process_time;
				ele_[op].energy_cost = ele_[op].update_energy_cost();
				jobs[ele_[op].jno] = ele_[op].et;

			}
		}
	}
	calcu_tot_tardiness();
	calcu_tot_energy_cost();
	return;
}

//粒子类函数实现
Particle::Particle()
{
	sol = Solution();
}

void Particle::print_()
{
	cout << sol.tot_tardiness << "	 " << sol.tot_energy_cost << endl;

}

void Particle::update_HT(Solution &sol, int index)
{
	/*for (int s=0;s<known.stage_count;++s)
	{
		int step = known.machine_count*known.pc[s].size();
		for (int m=0;m<known.machine_count;++m)
		{
			ele &e = sol.sol_1[s][m];
			for (int op=0;op<e.size();++op)
			{
				int key = e[op].jno*step + m * known.pc[s].size() + e[op].m_speed;
				HT[s][key] = 1;
			}
		}
	}*/
	vector<pair<long long, long long>> stas(known.stage_count);
	for (int s = 0; s < known.stage_count; ++s)
	{
		vector<int> js(known.job_count);
		pair<long long, long long> v(pair<long long, long long>(0, 0));
		for (int m = 0; m < known.machine_count; ++m)
		{
			ele &e = sol.sol_1[s][m];
			for (int op = 0; op < e.size(); ++op)
			{
				js[e[op].jno] = (e[op].machine * 10 + e[op].m_speed) * 100 + e[op].po;
				v.first += ((e[op].jno + 1) * 2 - 1)*js[e[op].jno];
				v.second += ((e[op].jno + 1) * 2)*js[e[op].jno];
			}
		}
		stas[s] = v;
	}
	pair<long long, long long> mv(pair<long long, long long>(0, 0));
	for (int s = 0; s < stas.size(); ++s)
	{
		mv.first += (2 * s + 1)*stas[s].first;
		mv.second += 2 * (s + 1)*stas[s].second;
	}
	int pre = mv.first / pow(10, 7);
	int back = mv.first % 1000000;
	int pre2 = mv.second / pow(10, 6);
	int back2 = mv.second % 100000;
	known.hs_1[index].first[pre] = known.hs_1[index].second[back] = true;
	known.hs_2[index].first[pre2] = known.hs_2[index].second[back2] = true;
	return;
}

bool Particle::func_HT(Solution &sol, int index)
{
	bool flag = false;
	vector<pair<long long, long long>> stas(known.stage_count);
	for (int s = 0; s < known.stage_count; ++s)
	{
		vector<int> js(known.job_count);
		pair<long long, long long> v(pair<long long, long long>(0, 0));
		for (int m = 0; m < known.machine_count; ++m)
		{
			ele &e = sol.sol_1[s][m];
			for (int op = 0; op < e.size(); ++op)
			{
				js[e[op].jno] = (e[op].machine * 10 + e[op].m_speed) * 100 + e[op].po;
				v.first += ((e[op].jno + 1) * 2 - 1)*js[e[op].jno];
				v.second += ((e[op].jno + 1) * 2)*js[e[op].jno];
			}
		}
		stas[s] = v;
	}
	pair<long long, long long> mv(pair<long long, long long>(0, 0));
	for (int s = 0; s < stas.size(); ++s)
	{
		mv.first += (2 * s + 1)*stas[s].first;
		mv.second += 2 * (s + 1)*stas[s].second;
	}
	int pre = mv.first / pow(10, 7);
	int back = mv.first % 1000000;
	int pre2 = mv.second / pow(10, 6);
	int back2 = mv.second % 100000;
	flag = (known.hs_1[index].first[pre] && known.hs_1[index].second[back] && known.hs_2[index].first[pre2] && known.hs_2[index].second[back2]);
	return flag;
}

void Particle::update_position(int &in)//当前粒子在群体中的位置
{
	Factory fac = Factory();
	double w_ = 0.;
	//newma:新机器 newpt:新处理时间  newp:新的位置 
	int newma = 0, newpt = 0, newp = 0, newop = 0;

	int tm = 0, tpt = 0, tpo = 0;//中间变量
	bool flag = false;
	Operation tmp_op;
	vector<vector<ele>> temp_sol;

	srand(clock() + rand() % 10000 + in);
	w_ = (rand() % 10000)*0.0001;

	if (w_ < known.w)
	{
		//进行变异
		for (int s = 0; s < known.stage_count; ++s)
		{
			for (int m = 0; m < known.machine_count; ++m)
			{
				ele &e = sol.sol_1[s][m];
				for (int op = 0; op < e.size(); ++op)
				{
					if (rand() % 5 == 0) //有待调整
					{
						tm = e[op].machine;
						tpt = e[op].m_speed;
						tpo = e[op].po;
						//随机选择一个机器上的操作和当前操作交换
						srand(clock());
						newma = rand() % known.machine_count;
						if (sol.sol_1[s][newma].size() == 0)
						{
							e[op].machine = newma;
							e[op].po = 0;
							sol.sol_1[s][newma].push_back(e[op]);
							e.erase(e.begin() + op);
							--op;
							continue;
						}
						else
						{
							//在新的机器上选一个新的位置
							newop = rand() % sol.sol_1[s][newma].size();

							e[op].machine = newma;
							e[op].m_speed = sol.sol_1[s][newma][newop].m_speed;
							e[op].po = sol.sol_1[s][newma][newop].po;
							sol.sol_1[s][newma][newop].machine = tm;
							sol.sol_1[s][newma][newop].m_speed = tpt;
							sol.sol_1[s][newma][newop].po = tpo;
						}
					}
				}
			}
		}
	}
	sol.update_sol();

	flag = fac.is_into_archive(sol);
	if (flag) fac.archive.push_back(sol);

	srand(clock() + rand() % 10000 + in);
	w_ = (rand() % 10000)*0.0001;
	if (w_ < known.c1)
	{
		//随机选一个粒子的pbest进行交叉
		newp = rand() % known.pops;

		Solution &pbest = Factory::par_swarm[newp].pbest;

		vector<bool> croPB(known.job_count); //croPB[i] = true:选择pbest；croPB[i]=false:保留自己
		//存放交叉后的结果
		temp_sol = vector<vector<ele>>(known.stage_count, vector<ele>(known.machine_count));

		for (int s = 0; s < known.stage_count; ++s)
		{
			for (int j = 0; j < known.job_count; ++j)
			{
				//srand(clock() + rand() % 100);
				croPB[j] = rand() % 2;
			}
			for (int m = 0; m < known.machine_count; ++m)
			{
				for (int op1 = 0, op2 = 0; op1 < sol.sol_1[s][m].size() || op2 < pbest.sol_1[s][m].size(); ++op1, ++op2)
				{
					if (op1 < sol.sol_1[s][m].size() && !croPB[sol.sol_1[s][m][op1].jno]) temp_sol[s][m].push_back(sol.sol_1[s][m][op1]);
					if (op2 < pbest.sol_1[s][m].size() && croPB[pbest.sol_1[s][m][op2].jno]) temp_sol[s][m].push_back(pbest.sol_1[s][m][op2]);
				}
			}
		}

		sol.sol_1 = temp_sol;
		sol.update_sol();

		flag = fac.is_into_archive(sol);
		if (flag) fac.archive.push_back(sol);
	}

	srand(clock() + rand() % 10000 + in);
	w_ = (rand() % 10000)*0.0001;
	if (w_ < known.c2)
	{
		vector<double> ras(known.job_count);
		temp_sol = vector<vector<ele>>(known.stage_count, vector<ele>(known.machine_count));

		for (int s = 0; s < known.stage_count; ++s)
		{
			for (int i = 0; i < known.job_count; ++i)
			{
				srand(clock() + rand() % 1000);
				ras[i] = (rand() % 10000) *0.0001;
			}
			for (int m = 0; m < known.machine_count; ++m)
			{
				for (int op1 = 0, op2 = 0; op1 < sol.sol_1[s][m].size() || op2 < gbest.sol_1[s][m].size(); ++op1, ++op2)
				{
					if (op1 < sol.sol_1[s][m].size() && ras[sol.sol_1[s][m][op1].jno] < 0.3) temp_sol[s][m].push_back(sol.sol_1[s][m][op1]);
					if (op2 < gbest.sol_1[s][m].size() && ras[gbest.sol_1[s][m][op2].jno] >= 0.3) temp_sol[s][m].push_back(gbest.sol_1[s][m][op2]);
				}

			}
		}
		sol.sol_1 = temp_sol;
		sol.update_sol();
	}

}


//工厂类函数实现
Factory::Factory()
{

}


/*生成初始解*/
void Factory::initial_sol()
{
	while (par_swarm.size() < known.pops)
	{
		Particle tmp_par = Particle();
		tmp_par.sol.generate_sol(par_swarm.size());
		tmp_par.pbest = tmp_par.sol;
		if (!par_swarm.size()) gbest = tmp_par.sol;
		else {
			if (tmp_par.sol.tot_tardiness <= gbest.tot_tardiness)
			{
				gbest = tmp_par.sol;
			}
		}
		/*更新archive*/
		bool flag = is_into_archive(tmp_par.sol);
		if (flag) archive.push_back(tmp_par.sol);
		par_swarm.push_back(tmp_par);
	}
}

/*控制函数*/
void Factory::control()
{
	known.st = clock();
	known.et = clock();
	bool flag = false;
	known.iter = 0;
	cout << "开始读文件......" << endl;
	known.read_file();
	known.print_file();
	known.write_head();
	cout << "结束读文件......" << endl;
	cout << "开始生成初始解......" << endl;
	
	//初始化禁忌表，每个禁忌动作映射为一个值
	known.tabuTable = vector<vector<vector<map<int, int>>>>(known.pops, 
		vector<vector<map<int, int>>>(known.stage_count, vector<map<int, int>>(known.job_count)));

	initial_sol();
	//print_pops(); //只输出sol的目标值
	//cout << endl;
	cout << "结束生成初始解......" << endl;
	cout << "迭代开始："  << '\n';

	if (known.pars_wl["al_name"] == "_HPSO_") 
	{
		//基于动作禁忌搜索
		control_based_move();
	}
	else if (known.pars_wl["al_name"] == "_HPSOII_") 
	{
		//just 局部搜索
		control_ls();
	}
	known.et = clock();
	known.runtime = (double)(known.et - known.st) / CLOCKS_PER_SEC;
	//cout << "用时：" << known.runtime<<"S" << '\n';
	writesol(); //将archive写入文件中

}

void Factory::control_based_move()
{
	bool flag = false;
	while (((double)(known.et - known.st) / CLOCKS_PER_SEC) < known.itr_time)
	{
		int pi = 0;
		for (auto &p : par_swarm)
		{
			//发散
			p.update_position(pi);
			//cout << "变异后" << p.sol.tot_tardiness << ' ' << p.sol.tot_energy_cost << endl;
			bool flag = is_into_archive(p.sol);
			if (flag)
			{
				archive.push_back(p.sol);
				Known::arch_et = clock();
			}
			if (p.sol < gbest)
			{
				gbest = p.sol;
			}
			else if (p.sol.tot_tardiness < gbest.tot_tardiness)
			{
				gbest = p.sol;
			}
			else if (p.sol.tot_energy_cost < gbest.tot_energy_cost)
			{
				//以一定概率更新，（多大概率呢？暂定0.5）
				if ((rand() % 1000) * 0.001 < 0.0)
					gbest = p.sol;
			}
			if (p.sol < p.pbest)
			{
				p.pbest = p.sol;
			}
			else if (p.sol.tot_tardiness < p.pbest.tot_tardiness)
				p.pbest = p.sol;
			else if (p.sol.tot_energy_cost < p.pbest.tot_energy_cost)
			{
				if ((rand() % 100) * 0.01 < 0.0)
					p.pbest = p.sol;
			}
			//先向tot_atrdiness收敛
			for (int ite = 0; ite <known.TS_iter; ++ite)
			{
				p.LS_bm2(pi);
				//cout << "搜索后：" << p.sol.tot_tardiness << ' ' << p.sol.tot_energy_cost << endl;
				flag = is_into_archive(p.sol);
				if (flag)
				{
					archive.push_back(p.sol);
					Known::arch_et = clock();
				}
				if (p.sol < gbest)
				{
					gbest = p.sol;
				}
				else if (p.sol.tot_tardiness < gbest.tot_tardiness)
				{
					gbest = p.sol;
				}
				else if (p.sol.tot_energy_cost < gbest.tot_energy_cost)
				{
					//以一定概率更新，（多大概率呢？暂定0.5）
					if ((rand() % 1000) * 0.0001 < 0.0)
						gbest = p.sol;
				}
				if (p.sol < p.pbest)
				{
					p.pbest = p.sol;
				}
				else if (p.sol.tot_tardiness < p.pbest.tot_tardiness)
					p.pbest = p.sol;
				else if (p.sol.tot_energy_cost < p.pbest.tot_energy_cost)
				{
					if ((rand() % 1000) * 0.001 < 0.00)
						p.pbest = p.sol;
				}
			}
			++pi;
			
		}
		++known.iter;
		//par_swarm[0].sol = gbest;
		if (known.iter % known.disturb_its == 0)
		{ //对粒子群进行扰动
			sort(archive.begin(), archive.end(), [](Solution &s1, Solution &s2) {return s1.tot_tardiness < s2.tot_tardiness; });
			for (int po = 0,a=0; po < known.pops &&a<archive.size(); ++po,++a)
			{
				//int spo = rand() % archive.size();
				//cout << spo <<endl;
				par_swarm[po].sol = archive[a];
			}
			//gbest = archive[rand() % archive.size()];
		}
		known.et = clock();
		//cout << known.iter << "  " << gbest.tot_tardiness << " " << gbest.tot_energy_cost << endl;
	}
}

void Factory::control_ls()
{
	bool flag = false;
	while ((known.et - known.st) / CLOCKS_PER_SEC < Known::itr_time)
	{
		int pi = 0;
		for (auto &p : par_swarm)
		{
			//发散
			p.update_position(pi);
			bool flag = is_into_archive(p.sol);
			if (flag)
			{
				archive.push_back(p.sol);
				Known::arch_et = clock();
			}
			if (p.sol < gbest)
			{
				gbest = p.sol;
			}
			else if (p.sol.tot_tardiness < gbest.tot_tardiness)
			{
				gbest = p.sol;
			}
			if (p.sol.tot_tardiness <= p.pbest.tot_tardiness)
				p.pbest = p.sol;
			else if (p.sol.tot_energy_cost < p.pbest.tot_energy_cost)
			{
				if ((rand() % 1000) * 0.001 < 0.05)
					p.pbest = p.sol;
			}
			
			//先向tot_atrdiness收敛
			for (int ite = 0; ite < known.TS_iter; ++ite)
			{
				p.local_search();
				bool flag = is_into_archive(p.sol);
				if (flag)
				{
					archive.push_back(p.sol);
					Known::arch_et = clock();
				}
				if (p.sol < gbest)
				{
					gbest = p.sol;
				}
				else if (p.sol.tot_tardiness < gbest.tot_tardiness)
				{
					gbest = p.sol;
				}
				if (p.sol.tot_tardiness <= p.pbest.tot_tardiness)
					p.pbest = p.sol;
				else if (p.sol.tot_energy_cost < p.pbest.tot_energy_cost)
				{
					if ((rand() % 1000) * 0.001 < 0.05)
						p.pbest = p.sol;
				}
			}
			++pi;
		}
		++known.iter;
		gbest = archive[rand() % archive.size()];
		if (known.iter % stoi(known.pars_wl["disturb_its"]) == 0)
		{
			for (int p = 0; p < stoi(known.pars_wl["pops"]); ++p)
			{
				int spo = rand() % archive.size();
				par_swarm[p].sol = archive[spo];
			}
		}
		known.et = clock();
		cout << known.iter << "  " << gbest.tot_tardiness << " " << gbest.tot_energy_cost << endl;
		
	}
}


void Particle::tabu_search(int &index)
{
	/*
	将当前粒子的邻域分为禁忌邻域和非禁忌邻域，优先从非禁忌邻域中选取解，没有非禁忌的：以一定的概率从禁忌的选取否则变异一下当前解
	*/
	int nm = 0, npt = 0, npo = 0;//新机器、新处理时间、新位置
	Solution tabubest = Solution(), non_tabubest = Solution();
	tabubest.tot_energy_cost = non_tabubest.tot_energy_cost = INT_MAX;
	tabubest.tot_tardiness = non_tabubest.tot_tardiness = INT_MAX;
	vector<Solution> tabu, nontabu;
	vector<Solution> better;
	vector<Solution> nobetter;
	Solution backup;
	bool flag = false;
	Factory fac = Factory();
	//改变处理速度
	for (int s = 0; s < known.stage_count; ++s)
	{
		for (int m = 0; m < known.machine_count; ++m)
		{
			ele &e = sol.sol_1[s][m];
			for (int n = 0; n < e.size(); ++n)
			{
				for (int v = 0; v < known.pc[s].size(); ++v)
				{
					if (sol.sol_1[s][m][n].m_speed == v) continue;
					Solution tmp_sol_ = this->sol;
					tmp_sol_.sol_1[s][m][n].m_speed = v;
					tmp_sol_.update_sol();
					//tmp_sol_.update_sol_no();

					flag = fac.is_into_archive(tmp_sol_);
					if (flag) Factory::archive.push_back(tmp_sol_);


					if (!func_HT(tmp_sol_, index))
					{
						//非禁忌--待改进
						if (tmp_sol_ < non_tabubest)
						{

							/*for (int i=0;i<nontabu.size();++i)
							{
								if (tmp_sol_ < nontabu[i])
								{
									nontabu.erase(nontabu.begin()+i);
									--i;
								}
							}*/
							/*sol = tmp_sol_;
							return;*/
							nontabu.push_back(tmp_sol_);
							non_tabubest = tmp_sol_;
						}
						else if (tmp_sol_ == non_tabubest)
						{
							nontabu.push_back(tmp_sol_);
						}
					}
					else
					{
						//禁忌
						if (tmp_sol_ < tabubest)
						{
							/*for (int i = 0; i < tabu.size(); ++i)
							{
								if (tmp_sol_ < tabu[i])
								{
									tabu.erase(tabu.begin() + i);
									--i;
								}
							}*/
							tabu.push_back(tmp_sol_);
							tabubest = tmp_sol_;
						}
						else if (tmp_sol_ == tabubest)
						{
							tabu.push_back(tmp_sol_);
						}
					}
				}
			}
		}
	}

	//改变分配的机器*/
	for (int s = 0; s < known.stage_count; ++s)
	{

		for (int m = 0; m < known.machine_count; ++m)
		{
			ele &e = sol.sol_1[s][m];
			for (int op = 0; op < e.size(); ++op)
			{

				Solution tmp_sol_ = this->sol;
				Operation tmp_o = e[op];

				nm = rand() % known.machine_count;
				if (sol.sol_1[s][nm].size() == 0)
				{
					tmp_sol_.sol_1[s][m].erase(tmp_sol_.sol_1[s][m].begin() + op);
					tmp_sol_.sol_1[s][nm].push_back(tmp_o);
				}
				else
				{
					npo = rand() % sol.sol_1[s][nm].size();
					tmp_sol_.sol_1[s][m][op] = tmp_sol_.sol_1[s][nm][npo];
					tmp_sol_.sol_1[s][nm][npo] = tmp_o;
				}
				tmp_sol_.update_sol();
				//tmp_sol_.update_sol_no();

				flag = fac.is_into_archive(tmp_sol_);
				if (flag) Factory::archive.push_back(tmp_sol_);
				if (!func_HT(tmp_sol_, index))
				{
					//非禁忌
					if (tmp_sol_ < non_tabubest)
					{
						/*for (int i = 0; i < nontabu.size(); ++i)
						{
							if (tmp_sol_ < nontabu[i])
								{
									nontabu.erase(nontabu.begin() + i);
									--i;
								}
							}*/
							/*sol = tmp_sol_;
							return;*/
						nontabu.push_back(tmp_sol_);
						non_tabubest = tmp_sol_;
					}
					else if (tmp_sol_ == non_tabubest)
					{
						nontabu.push_back(tmp_sol_);
					}
				}
				else
				{
					//禁忌
					if (tmp_sol_ < tabubest)
					{
						/*for (int i = 0; i < tabu.size(); ++i)
						{
								if (tmp_sol_ < tabu[i])
								{
									tabu.erase(tabu.begin() + i);
									--i;
								}
							}*/
						tabu.push_back(tmp_sol_);
						tabubest = tmp_sol_;
					}
					else if (tmp_sol_ == tabubest)
					{
						tabu.push_back(tmp_sol_);
					}
				}


			}
		}
	}

	//改变机器的工件顺序
	for (int s = 0; s < known.stage_count; ++s)
	{
		for (int m = 0; m < known.machine_count; ++m)
		{
			int size_ = sol.sol_1[s][m].size();
			for (int i = 0; i < size_; ++i)
			{
				for (int j = i + 1; j < size_; ++j)
				{
					Solution tmp_sol_ = sol;
					Operation tmp_o = tmp_sol_.sol_1[s][m][i];
					tmp_sol_.sol_1[s][m][i] = tmp_sol_.sol_1[s][m][j];
					tmp_sol_.sol_1[s][m][j] = tmp_o;
					tmp_sol_.update_sol();
					//tmp_sol_.update_sol_no();

					flag = fac.is_into_archive(tmp_sol_);
					if (flag) Factory::archive.push_back(tmp_sol_);
					if (!func_HT(tmp_sol_, index))
					{
						//非禁忌
						if (tmp_sol_ < non_tabubest)
						{
							/*for (int i = 0; i < nontabu.size(); ++i)
							{
								if (tmp_sol_ < nontabu[i])
									{
										nontabu.erase(nontabu.begin() + i);
										--i;
									}
								}*/
								/*sol = tmp_sol_;
								 return;*/
							nontabu.push_back(tmp_sol_);
							non_tabubest = tmp_sol_;
						}
						else if (tmp_sol_ == non_tabubest)
						{
							nontabu.push_back(tmp_sol_);
						}
					}
					else
					{
						//禁忌
						if (tmp_sol_ < tabubest)
						{
							/*for (int i = 0; i < tabu.size(); ++i)
							{
									if (tmp_sol_ < tabu[i])
									{
										tabu.erase(tabu.begin() + i);
										--i;
									}
								}*/
							tabu.push_back(tmp_sol_);
							tabubest = tmp_sol_;
						}
						else if (tmp_sol_ == tabubest)
						{
							tabu.push_back(tmp_sol_);
						}
					}
				}
			}
		}
	}
	update_HT(this->sol, index);

	if (nontabu.size() > 0)
	{
		//srand(clock() + tabu.size() + nontabu.size());
		int po = rand() % nontabu.size();
		this->sol = nontabu[po];
	}
	else if (tabu.size() > 0)
	{
		//srand(clock() + tabu.size() + nontabu.size());
		int po = rand() % tabu.size();
		this->sol = tabu[po];
	}
	return;
}

void Particle::local_search()
{
	//改变处理速度
	bool flag = false;
	Factory fac = Factory();
	for (int s = 0; s < known.stage_count; ++s)
	{
		for (int m = 0; m < known.machine_count; ++m)
		{
			ele &e = sol.sol_1[s][m];
			for (int n = 0; n < e.size(); ++n)
			{
				for (int v = 0; v < known.pc[s].size(); ++v)
				{
					if (sol.sol_1[s][m][n].m_speed == v) continue;
					Solution tmp_sol_ = sol;
					tmp_sol_.sol_1[s][m][n].m_speed = v;

					tmp_sol_.update_sol();
					flag = fac.is_into_archive(tmp_sol_);
					if (flag)
					{
						Factory::archive.push_back(tmp_sol_);
						Known::arch_et = clock();
					}
					if (tmp_sol_ < sol)
					{
						sol = tmp_sol_;
						return;
					}
					else if (tmp_sol_ == sol)
					{
						srand(clock());
						if ((rand() % 10000) * 0.0001 < known.acp)
						{
							sol = tmp_sol_;
							return;
						}
					}
				}
			}
		}
	}
	//改变分配的机器*/
	for (int s = 0; s < known.stage_count; ++s)
	{
		for (int m = 0; m < known.machine_count; ++m)
		{
			ele &e = sol.sol_1[s][m];
			for (int op = 0; op < e.size(); ++op)
			{
				for (int m2 = 0; m2 < known.machine_count; ++m2)
				{
					Solution tmp_sol_ = sol;
					Operation tmp_o = e[op];
					if (tmp_sol_.sol_1[s][m2].size() == 0)
					{
						tmp_o.machine = m2;
						tmp_o.po = 0;
						tmp_sol_.sol_1[s][m].erase(tmp_sol_.sol_1[s][m].begin() + op);
						tmp_sol_.sol_1[s][m2].push_back(tmp_o);
					}
					else
					{
						int mop = rand() % tmp_sol_.sol_1[s][m2].size();
						tmp_sol_.sol_1[s][m][op] = tmp_sol_.sol_1[s][m2][mop];
						tmp_sol_.sol_1[s][m2][mop] = tmp_o;
					}
					tmp_sol_.update_sol();
					flag = fac.is_into_archive(tmp_sol_);
					if (flag)
					{
						Factory::archive.push_back(tmp_sol_);
						Known::arch_et = clock();
					}
					if (tmp_sol_ < sol)
					{
						sol = tmp_sol_;
						return;
					}
					else if (tmp_sol_ == sol)
					{
						srand(clock());
						if ((rand() % 10000) * 0.0001 < known.acp)
						{
							sol = tmp_sol_;
							return;
						}
					}
				}

			}
		}
	}
	//改变机器的工件顺序
	for (int s = 0; s < known.stage_count; ++s)
	{
		for (int m = 0; m < known.machine_count; ++m)
		{
			int size_ = sol.sol_1[s][m].size();
			for (int i = 0; i < size_; ++i)
			{
				for (int j = i + 1; j < size_; ++j)
				{
					Solution tmp_sol_ = sol;
					Operation tmp_o = tmp_sol_.sol_1[s][m][i];
					tmp_sol_.sol_1[s][m][i] = tmp_sol_.sol_1[s][m][j];
					tmp_sol_.sol_1[s][m][j] = tmp_o;
					tmp_sol_.update_sol();
					flag = fac.is_into_archive(tmp_sol_);
					if (flag)
					{
						Factory::archive.push_back(tmp_sol_);
						Known::arch_et = clock();
					}
					if (tmp_sol_ < sol)
					{
						sol = tmp_sol_;
						return;
					}
					else if (tmp_sol_ == sol)
					{
						srand(clock());
						if ((rand() % 10000) * 0.0001 < known.acp)
						{
							sol = tmp_sol_;
							return;
						}
					}
				}
			}
		}
	}
	return;
}

void Particle::LS_bm(int & pos)
{
	/*
	将当前粒子的邻域分为禁忌邻域和非禁忌邻域，优先从非禁忌邻域中选取解，
	没有非禁忌的：以一定的概率从禁忌的选取否则变异一下当前解
	每个待选解都对应有一个禁忌操作，用来更新禁忌表
	*/
	int nm = 0, npt = 0, npo = 0;//新机器、新处理时间、新位置
	Solution tabubest = Solution(), non_tabubest = Solution();
	tabubest.tot_energy_cost = non_tabubest.tot_energy_cost = INT_MAX;
	tabubest.tot_tardiness = non_tabubest.tot_tardiness = INT_MAX;
	vector<pair<Solution,TabuOp>> tabu, nontabu;

	Solution backup;
	bool flag = false;
	Factory fac = Factory();
	TabuOp to = TabuOp();
	int tov = 0;//操作to映射成的解
	//改变处理速度
	for (int s = 0; s < known.stage_count; ++s)
	{
		to.sno = s;
		for (int m = 0; m < known.machine_count; ++m)
		{
			ele &e = sol.sol_1[s][m];
			to.mno = m;
			for (int n = 0; n < e.size(); ++n)
			{
				to.jno = e[n].jno;
				to.mv = sol.sol_1[s][m][n].m_speed;
				to.mpo = sol.sol_1[s][m][n].po;
				for (int v = 0; v < known.pc[s].size(); ++v)
				{
					if (sol.sol_1[s][m][n].m_speed == v) continue;

					pair<Solution,TabuOp> tmp_sol_(this->sol,to);
					tmp_sol_.first.sol_1[s][m][n].m_speed = v;
					tmp_sol_.first.update_sol();
					//tmp_sol_.update_sol_no();

					flag = fac.is_into_archive(tmp_sol_.first);
					if (flag) Factory::archive.push_back(tmp_sol_.first);

					tov = ((to.mno * 10) + to.mv) * 100 + to.mpo;
					if (known.iter >= known.tabuTable[pos][s][to.jno][tov])
					{
						//非禁忌--待改进
						if (tmp_sol_.first < non_tabubest)
						{

					
							/*sol = tmp_sol_;
							return;*/
							nontabu.push_back(tmp_sol_);
							non_tabubest = tmp_sol_.first;
						}
						else if (tmp_sol_.first == non_tabubest)
						{
							nontabu.push_back(tmp_sol_);
						}
					}
					else
					{
						//禁忌
						if (tmp_sol_.first < tabubest)
						{
							/*for (int i = 0; i < tabu.size(); ++i)
							{
								if (tmp_sol_ < tabu[i])
								{
									tabu.erase(tabu.begin() + i);
									--i;
								}
							}*/
							tabu.push_back(tmp_sol_);
							tabubest = tmp_sol_.first;
						}
						else if (tmp_sol_.first == tabubest)
						{
							tabu.push_back(tmp_sol_);
						}
					}
				}
			}
		}
	}

	//改变分配的机器*/
	for (int s = 0; s < known.stage_count; ++s)
	{
		to.sno = s;
		for (int m = 0; m < known.machine_count; ++m)
		{
			ele &e = sol.sol_1[s][m];
			to.mno = m;
			for (int op = 0; op < e.size(); ++op)
			{
				to.jno = e[op].jno;
				to.mv = e[op].m_speed;
				to.mpo = e[op].po;

				pair<Solution,TabuOp> tmp_sol_(this->sol,to);
				Operation tmp_o = e[op];

				nm = rand() % known.machine_count;
				if (sol.sol_1[s][nm].size() == 0)
				{
					tmp_sol_.first.sol_1[s][m].erase(tmp_sol_.first.sol_1[s][m].begin() + op);
					tmp_sol_.first.sol_1[s][nm].push_back(tmp_o);
				}
				else
				{
					npo = rand() % sol.sol_1[s][nm].size();
					tmp_sol_.first.sol_1[s][m][op] = tmp_sol_.first.sol_1[s][nm][npo];
					tmp_sol_.first.sol_1[s][nm][npo] = tmp_o;
				}
				tmp_sol_.first.update_sol();
				//tmp_sol_.update_sol_no();

				flag = fac.is_into_archive(tmp_sol_.first);
				if (flag) Factory::archive.push_back(tmp_sol_.first);

				tov = ((to.mno * 10) + to.mv) * 100 + to.mpo;
				if (known.iter >= known.tabuTable[pos][s][to.jno][tov])
				{
					//非禁忌
					if (tmp_sol_.first < non_tabubest)
					{
						/*for (int i = 0; i < nontabu.size(); ++i)
						{
							if (tmp_sol_ < nontabu[i])
								{
									nontabu.erase(nontabu.begin() + i);
									--i;
								}
							}*/
							/*sol = tmp_sol_;
							return;*/
						nontabu.push_back(tmp_sol_);
						non_tabubest = tmp_sol_.first;
					}
					else if (tmp_sol_.first == non_tabubest)
					{
						nontabu.push_back(tmp_sol_);
					}
				}
				else
				{
					//禁忌
					if (tmp_sol_.first < tabubest)
					{
						/*for (int i = 0; i < tabu.size(); ++i)
						{
								if (tmp_sol_ < tabu[i])
								{
									tabu.erase(tabu.begin() + i);
									--i;
								}
							}*/
						tabu.push_back(tmp_sol_);
						tabubest = tmp_sol_.first;
					}
					else if (tmp_sol_.first == tabubest)
					{
						tabu.push_back(tmp_sol_);
					}
				}
			}
		}
	}

	//改变机器的工件顺序
	for (int s = 0; s < known.stage_count; ++s)
	{
		to.sno = s;
		for (int m = 0; m < known.machine_count; ++m)
		{
			to.mno = m;
			int size_ = sol.sol_1[s][m].size();
			for (int i = 0; i < size_; ++i)
			{
				to.jno = sol.sol_1[s][m][i].jno;
				to.mpo = sol.sol_1[s][m][i].machine;
				to.mv = sol.sol_1[s][m][i].m_speed;
				for (int j = i + 1; j < size_; ++j)
				{
					pair<Solution,TabuOp> tmp_sol_(sol,to);
					Operation tmp_o = tmp_sol_.first.sol_1[s][m][i];
					tmp_sol_.first.sol_1[s][m][i] = tmp_sol_.first.sol_1[s][m][j];
					tmp_sol_.first.sol_1[s][m][j] = tmp_o;
					tmp_sol_.first.update_sol();
					//tmp_sol_.update_sol_no();

					flag = fac.is_into_archive(tmp_sol_.first);
					if (flag) Factory::archive.push_back(tmp_sol_.first);
					tov = ((to.mno * 10) + to.mv) * 100 + to.mpo;
					if (known.iter >= known.tabuTable[pos][s][to.jno][tov])
					{
						//非禁忌
						if (tmp_sol_.first < non_tabubest)
						{
							/*for (int i = 0; i < nontabu.size(); ++i)
							{
								if (tmp_sol_ < nontabu[i])
									{
										nontabu.erase(nontabu.begin() + i);
										--i;
									}
								}*/
								/*sol = tmp_sol_;
								 return;*/
							nontabu.push_back(tmp_sol_);
							non_tabubest = tmp_sol_.first;
						}
						else if (tmp_sol_.first == non_tabubest)
						{
							nontabu.push_back(tmp_sol_);
						}
					}
					else
					{
						//禁忌
						if (tmp_sol_.first < tabubest)
						{
							/*for (int i = 0; i < tabu.size(); ++i)
							{
									if (tmp_sol_ < tabu[i])
									{
										tabu.erase(tabu.begin() + i);
										--i;
									}
								}*/
							tabu.push_back(tmp_sol_);
							tabubest = tmp_sol_.first;
						}
						else if (tmp_sol_.first == tabubest)
						{
							tabu.push_back(tmp_sol_);
						}
					}
				}
			}
		}
	}
	
	if (nontabu.size() > 0)
	{
		srand(clock() + tabu.size() + nontabu.size());
		int po = rand() % nontabu.size();
		to = nontabu[po].second;
		this->sol = nontabu[po].first;
		update_tbb(to, pos);
	}
	else if (tabu.size() > 0)
	{
		srand(clock() + tabu.size() + nontabu.size());
		int po = rand() % tabu.size();
		to = tabu[po].second;
		this->sol = tabu[po].first;
		update_tbb(to, pos);
	}
	return;
}

void Particle::LS_bm2(int & pos)
{
	
	/*
	将当前粒子的邻域分为禁忌邻域和非禁忌邻域，优先从非禁忌邻域中选取解，
	没有非禁忌的：以一定的概率从禁忌的选取否则变异一下当前解
	每个待选解都对应有一个禁忌操作，用来更新禁忌表
	*/
	int nm = 0, npt = 0, npo = 0;//新机器、新处理时间、新位置
	vector<pair<Solution, TabuOp>> tabu, nontabu;

	bool flag = false;
	Factory fac = Factory();

	TabuOp to = TabuOp(); //用来存放需要禁忌的动作  
	int tov = 0;//操作to映射成的解

	srand(clock() + rand() % 100 + pos);
	int ns = rand() % Known::stage_count;

	/*可以考虑随机选一个工序开始*/
	//改变处理速度
	for (int s = ns; s < known.stage_count; ++s)
	{
		to.sno = s;
		for (int m = 0; m < known.machine_count; ++m)
		{
			ele &e = sol.sol_1[s][m];
			to.mno = m;
			for (int n = 0; n < e.size(); ++n)
			{
				to.jno = e[n].jno;
				to.mv = sol.sol_1[s][m][n].m_speed;
				to.mpo = sol.sol_1[s][m][n].po;
				for (int v = 0; v < known.pc[s].size(); ++v)
				{
					if (sol.sol_1[s][m][n].m_speed == v) continue;

					pair<Solution, TabuOp> tmp_sol_(this->sol, to); 
					tmp_sol_.first.sol_1[s][m][n].m_speed = v;
					tmp_sol_.first.update_sol();

					flag = fac.is_into_archive(tmp_sol_.first);
					if (flag)
					{
						Factory::archive.push_back(tmp_sol_.first);
						Known::arch_et = clock();
					}

					tov = ((to.mno * 10) + to.mv) * 100 + to.mpo;
					if (known.iter >= known.tabuTable[pos][s][to.jno][tov])
					{
						//非禁忌--待改进
						if (tmp_sol_.first < sol)
						{
							sol = tmp_sol_.first;
							
							update_tbb(to,pos);
							return;
						}else if (tmp_sol_.first.tot_tardiness < sol.tot_tardiness)
						{
							if ((rand() % 10000) * 0.0001 < stod(Known::pars_wl["acp1"]))
							{
								sol = tmp_sol_.first;
								update_tbb(to, pos);
								return;
							}
						}
						else if (tmp_sol_.first.tot_energy_cost < sol.tot_energy_cost)
						{
							if ((rand() % 10000) * 0.0001 < Known::acp)
							{
								sol = tmp_sol_.first;
								update_tbb(to, pos);
								return;
							}
						}
					}
					else
					{
						//禁忌
						if (tmp_sol_.first < sol)
						{
							tabu.push_back(tmp_sol_);
						}else if (tmp_sol_.first.tot_tardiness < sol.tot_tardiness)
						{
							if ((rand() % 10000) * 0.0001 < stod(Known::pars_wl["tb_acp1"]))
								tabu.push_back(tmp_sol_);
						}
						else if (tmp_sol_.first.tot_energy_cost < sol.tot_energy_cost)
						{
							if ((rand() % 10000) * 0.0001 < Known::tb_acp)
								tabu.push_back(tmp_sol_);
						}
					}
				}
			}
		}
	}

	for (int s = 0; s < ns; ++s)
	{
		to.sno = s;
		for (int m = 0; m < known.machine_count; ++m)
		{
			ele &e = sol.sol_1[s][m];
			to.mno = m;
			for (int n = 0; n < e.size(); ++n)
			{
				to.jno = e[n].jno;
				to.mv = sol.sol_1[s][m][n].m_speed;
				to.mpo = sol.sol_1[s][m][n].po;
				for (int v = 0; v < known.pc[s].size(); ++v)
				{
					if (sol.sol_1[s][m][n].m_speed == v) continue;

					pair<Solution, TabuOp> tmp_sol_(this->sol, to);
					tmp_sol_.first.sol_1[s][m][n].m_speed = v;
					tmp_sol_.first.update_sol();

					flag = fac.is_into_archive(tmp_sol_.first);
					if (flag)
					{
						Factory::archive.push_back(tmp_sol_.first);
						Known::arch_et = clock();
					}

					tov = ((to.mno * 10) + to.mv) * 100 + to.mpo;
					if (known.iter >= known.tabuTable[pos][s][to.jno][tov])
					{
						//非禁忌--待改进
						if (tmp_sol_.first < sol)
						{
							sol = tmp_sol_.first;
							update_tbb(to, pos);
							return;
						}
						else if (tmp_sol_.first.tot_tardiness < sol.tot_tardiness)
						{
							if ((rand() % 10000) * 0.0001 < stod(Known::pars_wl["acp1"]))
							{
								sol = tmp_sol_.first;
								update_tbb(to, pos);
								return;
							}
						}
						else if (tmp_sol_.first.tot_energy_cost < sol.tot_energy_cost)
						{
							if ((rand() % 10000) * 0.0001 < Known::acp)
							{
								sol = tmp_sol_.first;
								update_tbb(to, pos);
								return;
							}
						}
					}
					else
					{
						//禁忌
						if (tmp_sol_.first < sol)
						{
							tabu.push_back(tmp_sol_);
						}
						else if (tmp_sol_.first.tot_tardiness < sol.tot_tardiness)
						{
							if ((rand() % 10000) * 0.0001 < stod(Known::pars_wl["tb_acp1"]))
								tabu.push_back(tmp_sol_);
						}
						else if (tmp_sol_.first.tot_energy_cost < sol.tot_energy_cost)
						{
							if ((rand() % 10000) * 0.0001 < Known::tb_acp)
								tabu.push_back(tmp_sol_);
						}
					}
				}
			}
		}
	}
	//改变分配的机器*/
	for (int s = ns; s < known.stage_count; ++s)
	{
		to.sno = s;
		for (int m = 0; m < known.machine_count; ++m)
		{
			ele &e = sol.sol_1[s][m];
			to.mno = m;
			for (int op = 0; op < e.size(); ++op)
			{
				to.jno = e[op].jno;
				to.mv = e[op].m_speed;
				to.mpo = e[op].po;

				for (int m2 = 0; m2 < known.machine_count; ++m2)
				{
					if (m2 == m) continue;
					pair<Solution, TabuOp> tmp_sol_(this->sol, to);
					Operation tmp_o = e[op];
					if (tmp_sol_.first.sol_1[s][m2].size() == 0)
					{
						tmp_o.machine = m2;
						tmp_o.po = 0;
						tmp_sol_.first.sol_1[s][m].erase(tmp_sol_.first.sol_1[s][m].begin() + op);
						tmp_sol_.first.sol_1[s][m2].push_back(tmp_o);
					}
					else
					{
						int mop = rand() % tmp_sol_.first.sol_1[s][m2].size();
						tmp_sol_.first.sol_1[s][m][op] = tmp_sol_.first.sol_1[s][m2][mop];
						tmp_sol_.first.sol_1[s][m2][mop] = tmp_o;
					}
					tmp_sol_.first.update_sol();

					flag = fac.is_into_archive(tmp_sol_.first);
					if (flag)
					{
						Factory::archive.push_back(tmp_sol_.first);
						Known::arch_et = clock();
					}
					tov = ((to.mno * 10) + to.mv) * 100 + to.mpo;
					if (known.iter >= known.tabuTable[pos][s][to.jno][tov])
					{
						//非禁忌--待改进
						if (tmp_sol_.first < sol)
						{
							sol = tmp_sol_.first;
							update_tbb(to, pos);
							return;
						}
						else if (tmp_sol_.first.tot_tardiness < sol.tot_tardiness)
						{
							if ((rand() % 10000) * 0.0001 < stod(Known::pars_wl["acp1"]))
							{
								sol = tmp_sol_.first;
								update_tbb(to, pos);
								return;
							}
						}
						else if (tmp_sol_.first.tot_energy_cost < sol.tot_energy_cost)
						{
							if ((rand() % 10000) * 0.0001 < Known::acp)
							{
								sol = tmp_sol_.first;
								update_tbb(to, pos);
								return;
							}
						}
					}
					else
					{
						//禁忌
						if (tmp_sol_.first < sol)
						{
							tabu.push_back(tmp_sol_);
						}
						else if (tmp_sol_.first.tot_tardiness < sol.tot_tardiness)
						{
							if ((rand() % 10000) * 0.0001 < stod(Known::pars_wl["tb_acp1"]))
								tabu.push_back(tmp_sol_);
						}
						else if (tmp_sol_.first.tot_energy_cost < sol.tot_energy_cost)
						{
							if ((rand() % 10000) * 0.0001 < Known::tb_acp)
								tabu.push_back(tmp_sol_);
						}
					}
				}
			}
		}
	}

	//改变分配的机器*/
	for (int s = 0; s < ns; ++s)
	{
		to.sno = s;
		for (int m = 0; m < known.machine_count; ++m)
		{
			ele &e = sol.sol_1[s][m];
			to.mno = m;
			for (int op = 0; op < e.size(); ++op)
			{
				to.jno = e[op].jno;
				to.mv = e[op].m_speed;
				to.mpo = e[op].po;

				for (int m2=0;m2<known.machine_count;++m2) 
				{
					if (m2 == m) continue;
					pair<Solution, TabuOp> tmp_sol_(this->sol, to);
					Operation tmp_o = e[op];
					if (tmp_sol_.first.sol_1[s][m2].size() == 0)
					{
						tmp_o.machine = m2;
						tmp_o.po = 0;
						tmp_sol_.first.sol_1[s][m].erase(tmp_sol_.first.sol_1[s][m].begin() + op);
						tmp_sol_.first.sol_1[s][m2].push_back(tmp_o);
					}
					else
					{
						int mop = rand() % tmp_sol_.first.sol_1[s][m2].size();
						tmp_sol_.first.sol_1[s][m][op] = tmp_sol_.first.sol_1[s][m2][mop];
						tmp_sol_.first.sol_1[s][m2][mop] = tmp_o;
					}
					tmp_sol_.first.update_sol();

					flag = fac.is_into_archive(tmp_sol_.first);
					if (flag)
					{
						Factory::archive.push_back(tmp_sol_.first);
						Known::arch_et = clock();
					}
					tov = ((to.mno * 10) + to.mv) * 100 + to.mpo;
					if (known.iter >= known.tabuTable[pos][s][to.jno][tov])
					{
						//非禁忌--待改进
						if (tmp_sol_.first < sol)
						{
							sol = tmp_sol_.first;
							update_tbb(to, pos);
							return;
						}
						else if (tmp_sol_.first.tot_tardiness < sol.tot_tardiness)
						{
							if ((rand() % 10000) * 0.0001 < stod(Known::pars_wl["acp1"]))
							{
								sol = tmp_sol_.first;
								update_tbb(to, pos);
								return;
							}
						}
						else if (tmp_sol_.first.tot_energy_cost < sol.tot_energy_cost)
						{
							if ((rand() % 10000) * 0.0001 < Known::acp)
							{
								sol = tmp_sol_.first;
								update_tbb(to, pos);
								return;
							}
						}
					}
					else
					{
						//禁忌
						if (tmp_sol_.first < sol)
						{
							tabu.push_back(tmp_sol_);
						}
						else if (tmp_sol_.first.tot_tardiness < sol.tot_tardiness)
						{
							if ((rand() % 10000) * 0.0001 < stod(Known::pars_wl["tb_acp1"]))
								tabu.push_back(tmp_sol_);
						}
						else if (tmp_sol_.first.tot_energy_cost < sol.tot_energy_cost)
						{
							if ((rand() % 10000) * 0.0001 < Known::tb_acp)
								tabu.push_back(tmp_sol_);
						}
					}
				}
			}
		}
	}

	//改变机器的工件顺序
	for (int s = ns; s < Known::stage_count; ++s)
	{
		to.sno = s;
		for (int m = 0; m < known.machine_count; ++m)
		{
			to.mno = m;
			int size_ = sol.sol_1[s][m].size();
			for (int i = 0; i < size_; ++i)
			{
				to.jno = sol.sol_1[s][m][i].jno;
				to.mpo = i;
				to.mv = sol.sol_1[s][m][i].m_speed;
				for (int j = i + 1; j < size_; ++j)
				{
					pair<Solution, TabuOp> tmp_sol_(sol,  to);
					Operation tmp_o = tmp_sol_.first.sol_1[s][m][i];
					tmp_sol_.first.sol_1[s][m][i] = tmp_sol_.first.sol_1[s][m][j];
					tmp_sol_.first.sol_1[s][m][j] = tmp_o;
					tmp_sol_.first.update_sol();

					flag = fac.is_into_archive(tmp_sol_.first);
					if (flag)
					{
						Factory::archive.push_back(tmp_sol_.first);
						Known::arch_et = clock();
					}
					tov = ((to.mno * 10) + to.mv) * 100 + to.mpo;
					if (known.iter >= known.tabuTable[pos][s][to.jno][tov])
					{
						//非禁忌
						//非禁忌--待改进
						if (tmp_sol_.first < sol)
						{
							sol = tmp_sol_.first;
							update_tbb(to, pos);
							return;
						}
						else if (tmp_sol_.first.tot_tardiness < sol.tot_tardiness)
						{
							if ((rand() % 10000) * 0.0001 < stod(Known::pars_wl["acp1"]))
							{
								sol = tmp_sol_.first;
								update_tbb(to, pos);
								return;
							}
						}
						else if (tmp_sol_.first.tot_energy_cost < sol.tot_energy_cost)
						{
							if ((rand() % 10000) * 0.0001 < Known::acp)
							{
								sol = tmp_sol_.first;
								update_tbb(to, pos);
								return;
							}
						}
					}
					else
					{
						//禁忌
						if (tmp_sol_.first < sol)
						{
							tabu.push_back(tmp_sol_);
						}
						else if (tmp_sol_.first.tot_tardiness < sol.tot_tardiness)
						{
							if ((rand() % 10000) * 0.0001 < stod(Known::pars_wl["tb_acp1"]))
								tabu.push_back(tmp_sol_);
						}
						else if (tmp_sol_.first.tot_energy_cost < sol.tot_energy_cost)
						{
							if ((rand() % 10000) * 0.0001 < Known::tb_acp)
								tabu.push_back(tmp_sol_);
						}
					}
				}
			}
		}
	}

	for (int s = 0; s < ns; ++s)
	{

		to.sno = s;
		for (int m = 0; m < known.machine_count; ++m)
		{
			to.mno = m;
			int size_ = sol.sol_1[s][m].size();
			for (int i = 0; i < size_; ++i)
			{
				to.jno = sol.sol_1[s][m][i].jno;
				to.mpo = i;
				to.mv = sol.sol_1[s][m][i].m_speed;
				for (int j = i + 1; j < size_; ++j)
				{
					pair<Solution, TabuOp> tmp_sol_(sol, to);
					Operation tmp_o = tmp_sol_.first.sol_1[s][m][i];
					tmp_sol_.first.sol_1[s][m][i] = tmp_sol_.first.sol_1[s][m][j];
					tmp_sol_.first.sol_1[s][m][j] = tmp_o;
					tmp_sol_.first.update_sol();

					flag = fac.is_into_archive(tmp_sol_.first);
					if (flag)
					{
						Factory::archive.push_back(tmp_sol_.first);
						Known::arch_et = clock();
					}
					tov = ((to.mno * 10) + to.mv) * 100 + to.mpo;
					if (known.iter >= known.tabuTable[pos][s][to.jno][tov])
					{
						//非禁忌
						//非禁忌--待改进
						if (tmp_sol_.first < sol)
						{
							sol = tmp_sol_.first;
							update_tbb(to, pos);
							return;
						}
						else if (tmp_sol_.first.tot_tardiness < sol.tot_tardiness)
						{
							if ((rand() % 10000) * 0.0001 < stod(Known::pars_wl["acp1"]))
							{
								sol = tmp_sol_.first;
								update_tbb(to, pos);
								return;
							}
						}
						else if (tmp_sol_.first.tot_energy_cost < sol.tot_energy_cost)
						{
							if ((rand() % 10000) * 0.0001 < Known::acp)
							{
								sol = tmp_sol_.first;
								update_tbb(to, pos);
								return;
							}
						}
					}
					else
					{
						//禁忌
						if (tmp_sol_.first < sol)
						{
							tabu.push_back(tmp_sol_);
						}
						else if (tmp_sol_.first.tot_tardiness < sol.tot_tardiness)
						{
							if ((rand() % 10000) * 0.0001 < stod(Known::pars_wl["tb_acp1"]))
								tabu.push_back(tmp_sol_);
						}
						else if (tmp_sol_.first.tot_energy_cost < sol.tot_energy_cost)
						{
							if ((rand() % 10000) * 0.0001 < Known::tb_acp)
								tabu.push_back(tmp_sol_);
						}
					}
				}
			}
		}
	}
	if (tabu.size() > 0)
	{
		srand(clock() + tabu.size());
		int po = rand() % tabu.size();
		to = tabu[po].second;
		this->sol = tabu[po].first;
		update_tbb(to, pos);
	}
	return;
}

void Particle::update_tbb(TabuOp & tbo, int & pos)
{
	int tbov = (tbo.mno * 10 + tbo.mv) * 100 + tbo.mpo;
	known.tabuTable[pos][tbo.sno][tbo.jno][tbov] = known.iter+rand()%10;
}



bool Factory::is_into_archive(Solution & sol)
{
	bool flag = false, flag_ = true;
	if (!archive.size())
	{
		flag = true;
	}
	else
	{
		for (int i = 0; i < archive.size(); ++i)
		{
			if (sol.tot_tardiness == archive[i].tot_tardiness
				&& sol.tot_energy_cost == archive[i].tot_energy_cost)
				return false;
			if (sol < archive[i])
			{
				flag = true;
				archive.erase(archive.begin() + i);
				--i;
			}
			else if (sol == archive[i])
			{
				flag = true;
			}
			else
			{
				flag_ = false;
			}
		}
	}
	return flag && flag_;
}


void Factory::print_archive()
{
	cout << "print archive:" << '\n';
	sort(archive.begin(), archive.end(), [](Solution &s1, Solution &s2) {return s1.tot_tardiness < s2.tot_tardiness; });
	for (auto &sol : archive)
	{
		cout << sol.tot_tardiness << "	" << sol.tot_energy_cost << '\n';
	}
}

void Factory::print_pops()
{
	for (int p = 0; p < stoi(known.pars_wl["pops"]); ++p)
	{
		cout << "第" << p << "个粒子:    ";
		par_swarm[p].print_();
	}
}

void Factory::random_par()
{
	while (par_swarm.size() < stoi(known.pars_wl["pops"]))
	{
		Particle tmp_par = Particle();
		int r = rand() % 3;
		if (r == 1)
			tmp_par.sol.generate_sol_(par_swarm.size());
		else
			tmp_par.sol.generate_sol(par_swarm.size());
		tmp_par.pbest = tmp_par.sol;
		gbest = gbest < tmp_par.sol ? gbest : tmp_par.sol;

		/*更新archive*/
		bool flag = is_into_archive(tmp_par.sol);
		if (flag) archive.push_back(tmp_par.sol);
		par_swarm.push_back(tmp_par);
	}
	/*srand(clock() + iter + rand()%10+1);
	Particle par = Particle();
	Operation oper = Operation();
	for (int s = 0; s < known.stage_count; ++s)
	{
		oper.sno = s;
		for (int o = 0; o < job_count; ++o)
		{
			oper.jno = o;
			//随机生成一个机器
			srand(clock() + iter + rand() % 20);
			oper.machine = rand() % machine_count;
			oper.m_speed = rand() % pt[s].size();
			par.sol.sol_1[s][oper.machine].push_back(oper);
		}
	}

	par.sol.update_sol();
	par.pbest = par.sol;
	gbest = gbest < par.pbest ? gbest : par.pbest;
	return par;*/
	return;
}

void Factory::infile()
{
	string iname = "gbest";
	iname += ".txt";
	ofstream onf(iname, ios::out);
	if (!onf)
	{
		cout << "error" << endl;
		exit(0);
	}
	onf << "*** display solution " << " ***, tt: " << gbest.tot_tardiness << "\ttec:" << gbest.tot_energy_cost << '\n';
	// for each machine to display solution
	for (int stage_i = 0; stage_i < known.stage_count; stage_i++)
	{
		onf << "stage " << stage_i + 1 << ":" << '\n';
		for (int mach_i = 0; mach_i < known.machine_count; mach_i++)
		{
			ele &e = gbest.sol_1[stage_i][mach_i];
			onf << "m " << mach_i + 1 << ", " << e.size() << ":\t";
			for (int oper_i = 0; oper_i < e.size(); oper_i++)
				onf << e[oper_i].jno + 1 << "," << e[oper_i].sno + 1 << "\t"
				<< (e[oper_i].process_time - known.pc[stage_i][e[oper_i].m_speed]) << "\t" << known.pc[stage_i][e[oper_i].m_speed] << "\t";
			onf << '\n';
		}
	}
}

void Factory::writesol()
{
	//将archive写入文件中
	double fin = (double)(known.arch_et - known.st) / CLOCKS_PER_SEC;
	ofstream ofile(Known::pars_wl["ofn"], ios::app); //以追加的方式写入
	ofile << Known::pars_wl["run_num"] << '\t' << Known::iter << '\t'
		<< Known::runtime << '\t' << archive.size() << endl;
	ofile << fin << endl;
	sort(archive.begin(), archive.end(), [](Solution &s1, Solution &s2) {return s1.tot_tardiness < s2.tot_tardiness; });
	for (int i = 0; i < archive.size(); ++i)
	{
		ofile << archive[i].tot_tardiness << "  " << archive[i].tot_energy_cost << ", ";
	}
	ofile << endl;
	ofile.close();
	 
	/*Known::pars_wl["ofn2"] += (Known::pars_wl["job_id"]+Known::pars_wl["al_name"] + Known::pars_wl["fn"]);
    ofile.open(Known::pars_wl["ofn2"], ios::app); //以追加的方式写入
	ofile << fin << "  ";
	for (int i = 0; i < archive.size(); ++i)
	{
		ofile << archive[i].tot_tardiness << "  " << archive[i].tot_energy_cost << ", ";
	}
	ofile << endl;
	ofile.close();*/
}

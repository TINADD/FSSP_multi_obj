// MO_NSGAII.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
#include "NSGAII.h"

Known kno = Known();
vector<Individual> Fac::populations = vector<Individual>();//种群
vector<Individual>  Fac::offSpring = vector<Individual>();//offSpring.size == pops; offSpring = GA(populations);	
vector<Individual> Fac::Rt = vector<Individual>();//Rt = offSpring + populations ;Rt.size == 2*pops;
//写文件
void Fac::writesol()
{
	double fin = (double)(Known::arch_et - Known::st) / CLOCKS_PER_SEC;
	Known::runtime = (double)(Known::et - Known::st) / CLOCKS_PER_SEC;	

	string filename = Known::pars_wl["ofn"];
	ofstream ofile(filename, ios::app); //以追加的方式写入
	ofile << Known::pars_wl["run_num"] << '\t' << Known::iter << '\t'
		<< Known::runtime << '\t' << populations.size() << endl;
	ofile << fin << endl;
	sort(populations.begin(), populations.end(), [](Individual &i1, Individual &i2) {return i1.tot_tardiness < i2.tot_tardiness; });
	for (int i = 0; i < populations.size(); ++i)
	{
		ofile << populations[i].tot_tardiness << "  " << populations[i].tot_energy_cost << ", ";
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

void Fac::push_into_arch(Individual &sol)
{
	bool flag = false, flag_ = true;
	if (!archive.size())
	{
		archive.push_back(sol);
		return;
	}
	else
	{
		for (int i = 0; i < archive.size(); ++i)
		{
			if (sol.tot_tardiness == archive[i].tot_tardiness
				&& sol.tot_energy_cost == archive[i].tot_energy_cost)
				return;
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
		if (flag && flag_)
		{
			archive.push_back(sol);
			Known::arch_et = clock();
		}
	}
	return;
}

/*打印种群*/
void Fac::print_pops(vector<Individual> &popu)
{
	for (int indi = 0; indi < popu.size(); ++indi)
	{
		cout << "个体" << indi << "	" << popu[indi].tot_tardiness << "	" << popu[indi].tot_energy_cost << endl;
	}
}

/*局部搜索*/
void Individual::LS_based_move(int &pos)
{
	/*
	将当前粒子的邻域分为禁忌邻域和非禁忌邻域，优先从非禁忌邻域中选取解，
	没有非禁忌的：以一定的概率从禁忌的选取否则变异一下当前解
	每个待选解都对应有一个禁忌操作，用来更新禁忌表
	*/
	int nm = 0, npt = 0, npo = 0;//新机器、新处理时间、新位置
	Individual tabubest = Individual(), non_tabubest = Individual();
	tabubest.tot_energy_cost = non_tabubest.tot_energy_cost = INT_MAX;
	tabubest.tot_tardiness = non_tabubest.tot_tardiness = INT_MAX;
	vector<pair<Individual, TabuOp>> tabu, nontabu;

	bool flag = false;
	Fac fac = Fac();
	TabuOp to = TabuOp();
	int tov = 0;//操作to映射成的解
	int jno1 = 0,jno2=0;

	srand(clock() + rand() % 100 + pos);
	int ns = rand()%Known::stage_count;//从当前工序开始寻找邻域,first improved
	//改变处理速度
	for (int s = ns; s < Known::stage_count; ++s)
	{
		to.sno = s;
		for (int m = 0; m < Known::machine_count; ++m)
		{
			vector<int> &e = m_job[s][m];
			to.mno = m;
			for (int n = 0; n < e.size(); ++n)
			{
				to.jno = e[n];
				to.mv = sol[s][to.jno].m_speed;
				to.mpo = sol[s][to.jno].po;

				for (int v = 0; v < Known::pc[s].size(); ++v)
				{
					if (to.mv == v) continue;

					pair<Individual, TabuOp> tmp_sol_(*this, to);
					tmp_sol_.first.sol[s][to.jno].m_speed = v;
					tmp_sol_.first.update_sol();
					//tmp_sol_.update_sol_no();

					tov = ((to.mno * 10) + to.mv) * 100 + to.mpo;
					if (Known::iter >= Known::tabuTable[pos][s][to.jno][tov])
					{
						//非禁忌--待改进
						if (tmp_sol_.first < non_tabubest)
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
	for (int s = 0; s < Known::stage_count; ++s)
	{
		to.sno = s;
		for (int m = 0; m < Known::machine_count; ++m)
		{
			vector<int> &e = m_job[s][m];
			to.mno = m;
			for (int op = 0; op < e.size(); ++op)
			{
				to.jno = e[op];
				to.mv = sol[s][to.jno].m_speed;
				to.mpo = sol[s][to.jno].po;

				pair<Individual, TabuOp> tmp_sol_(*this, to);
				jno1 = to.jno;
				nm = rand() % Known::machine_count;
				if (m_job[s][nm].size() == 0)
				{
					for (int o = to.mpo+1; o < tmp_sol_.first.m_job[s][m].size(); ++o)
					{
						int jno = tmp_sol_.first.m_job[s][m][o];
						tmp_sol_.first.sol[s][jno].po--;
					}
					tmp_sol_.first.m_job[s][m].erase(tmp_sol_.first.m_job[s][m].begin() + to.mpo);
					tmp_sol_.first.m_job[s][nm].push_back(to.jno);
					tmp_sol_.first.sol[s][to.jno].machine = nm;
					tmp_sol_.first.sol[s][to.jno].po = 0;
				}
				else
				{
					npo = rand() % m_job[s][nm].size();
					jno2 = m_job[s][nm][npo];
					//更新机器上的工件号
					swap(tmp_sol_.first.m_job[s][m][to.mpo], tmp_sol_.first.m_job[s][nm][npo]);
					//更新工件的位置
					swap(tmp_sol_.first.sol[s][jno1].po, tmp_sol_.first.sol[s][jno2].po);
					//更新工件的机器号
					tmp_sol_.first.sol[s][jno1].machine = nm;
					tmp_sol_.first.sol[s][jno2].machine = m;
				}
				tmp_sol_.first.update_sol();
				//tmp_sol_.update_sol_no();

				tov = ((to.mno * 10) + to.mv) * 100 + to.mpo;
				if (Known::iter >= Known::tabuTable[pos][s][to.jno][tov])
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
	for (int s = 0; s < Known::stage_count; ++s)
	{
		to.sno = s;
		for (int m = 0; m < Known::machine_count; ++m)
		{
			to.mno = m;
			for (int i = 0; i < m_job[s][m].size(); ++i)
			{
				jno1 = m_job[s][m][i];
				to.jno = jno1;
				to.mpo = m;
				to.mv = sol[s][jno1].m_speed;
				for (int j = i + 1; j < m_job[s][m].size(); ++j)
				{
					jno2 = m_job[s][m][j];
					pair<Individual, TabuOp> tmp_sol_(*this, to);
					Operation tmp_o = tmp_sol_.first.sol[s][m];
					tmp_sol_.first.sol[s][jno1].po = j;
					tmp_sol_.first.sol[s][jno2].po = i;
					swap(tmp_sol_.first.m_job[s][m][i], tmp_sol_.first.m_job[s][m][j]);
					tmp_sol_.first.update_sol();
					
					tov = ((to.mno * 10) + to.mv) * 100 + to.mpo;
					if (Known::iter >= Known::tabuTable[pos][s][to.jno][tov])
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
		*this = nontabu[po].first;
		update_tbb(nontabu[po].second, pos);
	}
	else if (tabu.size() > 0)
	{
		srand(clock() + tabu.size() + nontabu.size());
		int po = rand() % tabu.size();
		*this = tabu[po].first;
		update_tbb(tabu[po].second, pos);

	}
	return;
}

void Individual::local_search(int & pos)
{
	srand(clock() + rand() % 100 + pos);
	int ns = rand() % Known::stage_count;
	int jno = 0;
	int npo = 0;
	for (int s=ns;s<Known::stage_count;++s) 
	{
		for (int m=0;m<Known::machine_count;++m) 
		{
			vector<int> &e = m_job[s][m];
			for (int j=0;j<e.size();++j) 
			{
				jno = e[j];
				//改变工件处理时间
				for (int v=0;v<Known::pc[s].size();++v) 
				{
					if (v == sol[s][jno].m_speed) continue;
					Individual tmp = *this;
					tmp.sol[s][jno].m_speed = v;
					tmp._update_sol();
					if (tmp < *this) 
					{
						*this = tmp;
						return;
					}
					else if (tmp.tot_tardiness < this->tot_tardiness) 
					{
						if ((rand() % 10000) * 0.0001 < stod(Known::pars_wl["acp"]))
						{
							*this = tmp;
							return;
						}
					}
				}

				for (int m2=0;m2<Known::machine_count;++m2) 
				{
					if (m2 == m) continue;
					Individual tmp = *this;
					tmp.sol[s][jno].machine = m2;
					tmp.m_job[s][m].erase(tmp.m_job[s][m].begin()+j);
					if (tmp.m_job[s][m2].size() == 0) tmp.m_job[s][m2].push_back(jno);
					else
					{
						npo = rand() % tmp.m_job[s][m2].size();
						tmp.m_job[s][m2].insert(tmp.m_job[s][m2].begin()+npo,jno);
					}
					tmp._update_sol();
					if (tmp < *this)
					{
						*this = tmp;
						return;
					}
					else if (tmp.tot_tardiness < this->tot_tardiness)
					{
						if ((rand() % 10000) * 0.0001 < stod(Known::pars_wl["acp"]))
						{
							*this = tmp;
							return;
						}
					}
				
				}
				//改变工件顺序
				for (int jc = j+1;jc<e.size();++jc)
				{
					Individual tmp = *this;
					swap(tmp.m_job[s][m][j],tmp.m_job[s][m][jc]);
					tmp._update_sol();
					if (tmp < *this)
					{
						*this = tmp;
						return;
					}
					else if (tmp.tot_tardiness < this->tot_tardiness)
					{
						if ((rand() % 10000) * 0.0001 < stod(Known::pars_wl["acp"]))
						{
							*this = tmp;
							return;
						}
					}
				}
			}
		}
	}


	for (int s=0;s<ns;++s) 
	{
		for (int m = 0; m < Known::machine_count; ++m)
		{
			vector<int> &e = m_job[s][m];
			for (int j = 0; j < e.size(); ++j)
			{
				jno = e[j];
				//改变工件处理时间
				for (int v = 0; v < Known::pc[s].size(); ++v)
				{
					if (v == sol[s][jno].m_speed) continue;
					Individual tmp = *this;
					tmp.sol[s][jno].m_speed = v;
					tmp._update_sol();
					if (tmp < *this)
					{
						*this = tmp;
						return;
					}
					else if (tmp == *this && tmp.tot_tardiness < this->tot_tardiness)
					{
						if ((rand() % 10000) * 0.0001 < stod(Known::pars_wl["acp"]))
						{
							*this = tmp;
							return;
						}
					}
				}

				for (int m2 = 0; m2 < Known::machine_count; ++m2)
				{
					if (m2 == m) continue;
					Individual tmp = *this;
					tmp.sol[s][jno].machine = m2;
					tmp.m_job[s][m].erase(tmp.m_job[s][m].begin() + j);
					if (tmp.m_job[s][m2].size() == 0) tmp.m_job[s][m2].push_back(jno);
					else
					{
						npo = rand() % tmp.m_job[s][m2].size();
						tmp.m_job[s][m2].insert(tmp.m_job[s][m2].begin() + npo, jno);
					}
					tmp._update_sol();
					if (tmp < *this)
					{
						*this = tmp;
						return;
					}
					else if (tmp == *this && tmp.tot_tardiness < this->tot_tardiness)
					{
						if ((rand() % 10000) * 0.0001 < stod(Known::pars_wl["acp"]))
						{
							*this = tmp;
							return;
						}
					}

				}
				//改变工件顺序
				for (int jc = j + 1; jc < e.size(); ++jc)
				{
					Individual tmp = *this;
					swap(tmp.m_job[s][m][j], tmp.m_job[s][m][jc]);
					tmp._update_sol();
					if (tmp < *this)
					{
						*this = tmp;
						return;
					}
					else if (tmp == *this && tmp.tot_tardiness < this->tot_tardiness)
					{
						if ((rand() % 10000) * 0.0001 < stod(Known::pars_wl["acp"]))
						{
							*this = tmp;
							return;
						}
					}
				}
			}
		}
	}
	return;
}

Fac::Fac()
{
	for (int s = 0; s < Known::stage_count; ++s)
	{ //遍历每个阶段
		for (int j = 0; j < Known::job_count; ++j)
		{ //遍历当前阶段的每个工件
			//随机生成一个机器以及对应的速度
			srand(clock());
			int mm = rand() % Known::machine_count;
			//	cout << mm << "-";
				//srand(clock());
			int mv = Known::pc[s][rand() % Known::pc[s].size()];
			//	cout << mv<<" ";

		}
		//	cout << endl;
	}
}

/*启动函数*/
void Fac::start_up()
{
	cout << "--程序开始--" << endl;
	/*读取文件filename，初始化调度信息*/
	cout << "...读文件..." << endl;
	kno.read_file();
	kno.print_file();

	kno.write_head();

	/*随机构造初始解*/
	cout << "...构造初始解..." << endl;
	init_pops();

	print_pops(populations);
	cout << "启发式算法开始..." << endl;
	 
	/*核心算法运行*/
	find_opt();

	start_down();
}


/*随机生成种群*/
void Fac::init_pops()
{
	while (populations.size() != Known::pops)
	{
		Individual tmp_indi = Individual();
		int rm = populations.size();
		tmp_indi._init_indi(rm); //有策略的随机生成较优解
		tmp_indi._update_sol(); //主动调度
		//tmp_indi.print_sol_();
		push_into_arch(tmp_indi);
		populations.push_back(tmp_indi);
	}
}


void Fac::NSGAII()
{
	Rt = populations;
	for (int i = 0; i < offSpring.size(); ++i)
	{
		//push_into_arch(offSpring[i]);
		Rt.push_back(offSpring[i]);
	}
	/*sort(Rt.begin(), Rt.end(), [](Individual &i1, Individual &i2) {return i1.tot_tardiness < i2.tot_tardiness; });
	cout << "Rt:" << endl;
	for (int rt=0;rt<Rt.size();++rt) 
	{
		cout << "第" << rt << "粒子: " << Rt[rt].tot_tardiness << " " << Rt[rt].tot_energy_cost << endl;
	}*/
	//cout << "开始分层... " << endl;
	/*
		将合并后的解分层,non_domi[i][j]=个体在种群Rt中的下标
	*/
	vector<vector<int>> non_domi = fast_non_dominated_sort();

	/*	cout << "分层后的结果: (输出的是个体的下标)" << endl;
		for (int le = 0; le < non_domi.size(); ++le)
		{
			cout << "level " << le << ":	";
			for (const auto x : non_domi[le])
			{
				cout << x << "	";
			}
			cout << endl;
		}*/
	//cout << "分层结束" << endl;
	vector<Individual> new_p;//存储分层后的前n项
	for (int i = 0; i < non_domi.size(); ++i)
	{
		if (non_domi[i].size() + new_p.size() <= Known::pops)
		{
			for (auto indi : non_domi[i])
			{
				new_p.push_back(Rt[indi]);
			}
		}
		else if (new_p.size() < Known::pops)
		{
			/*计算non_domi[i]的拥挤距离*/
			//cout << "计算拥挤距离开始..." << endl;
			crowding_distance_assignment(non_domi[i]);
			//cout << "计算拥挤距离结束..." << endl;
			int count = 0;
			while (new_p.size() < Known::pops &&count < non_domi[i].size())
			{
				new_p.push_back(Rt[non_domi[i][count++]]);
			}

		}
	}
	populations = new_p;
}

/*拥挤距离计算*/
void Fac::crowding_distance_assignment(vector<int> &non_domi) {
	int l = non_domi.size();
	for (int i : non_domi)
	{
		Rt[i].cd = 0.;
	}
	//cout << "tt排序start: " << endl;
	sort(non_domi.begin(), non_domi.end(), [this](int i1, int i2) {return Rt[i1].tot_tardiness < Rt[i2].tot_tardiness; });
	Rt[non_domi[0]].cd = Rt[non_domi[l - 1]].cd = INT_MAX;
	for (int i = 1; i < l - 1; ++i)
	{
		Rt[non_domi[i]].cd += (double)(Rt[non_domi[i + 1]].tot_tardiness - Rt[non_domi[i - 1]].tot_tardiness) / (double)(Rt[non_domi[l - 1]].tot_tardiness - Rt[non_domi[0]].tot_tardiness);
	}
	//cout << "tt排序end..." << endl;
	//cout << "ect排序start..." << endl;
	sort(non_domi.begin(), non_domi.end(), [this](int i1, int i2) {return Rt[i1].tot_energy_cost < Rt[i2].tot_energy_cost; });
	Rt[non_domi[0]].cd = Rt[non_domi[l - 1]].cd = INT_MAX;
	for (int i = 1; i < l - 1; ++i)
	{
		Rt[non_domi[i]].cd += (double)(Rt[non_domi[i + 1]].tot_energy_cost - Rt[non_domi[i - 1]].tot_energy_cost) / (Rt[non_domi[l - 1]].tot_energy_cost - Rt[non_domi[0]].tot_energy_cost);
	}
	//	cout << "ect排序end..." << endl;
	sort(non_domi.begin(), non_domi.end(), [this](int i1, int i2) {return Rt[i1].cd > Rt[i2].cd; });
}


/*遗传Algorithm*/
void Fac::GA(int &iter)
{
	//遗传算法--生成子代
	offSpring.clear();
	int cross2 = 0;
	int jno = 0;
	//先交叉
	for (int p=0;p<populations.size();++p) 
	{
		//随机选一个个体跟当前粒子交叉
		srand(clock() + rand()%100 + p);
		cross2 = rand() % populations.size();
		while (cross2 == p)
		{
			cross2 = rand() % populations.size();
		}
		Individual tid = Individual();
		for (int s=0;s<Known::stage_count;++s) 
		{
			for (int j=0;j<Known::job_count;++j) 
			{
				if (rand() % 1000 * 0.001 < 0.5) 
				{//选取当前粒子的op
					tid.sol[s][j] = populations[p].sol[s][j];
					tid.sol[s][j].po = tid.m_job[s][tid.sol[s][j].machine].size();
					tid.m_job[s][tid.sol[s][j].machine].push_back(j);
				}else
				{//选择另一个个体的op
					tid.sol[s][j] = populations[cross2].sol[s][j];
					tid.sol[s][j].po = tid.m_job[s][tid.sol[s][j].machine].size();
					tid.m_job[s][tid.sol[s][j].machine].push_back(j);
				}
			}
		}
		offSpring.push_back(tid);
	}

	//输出交叉后的offSpring
	/*for (int p=0;p<offSpring.size();++p) 
	{
		cout << "第" << p << "个: " << endl;
		for (int s=0;s<Known::stage_count;++s) 
		{
			for (int j=0;j<Known::job_count;++j) 
			{
				 cout << "[ " << offSpring[p].sol[s][j].machine << ","<< offSpring[p].sol[s][j].m_speed << "," << offSpring[p].sol[s][j].po << " ]";
				
			}
			for (int m=0;m<Known::machine_count;++m) 
			{
				for (int op=0;op<offSpring[p].m_job[s][m].size();++op) 
				{
					cout << offSpring[p].m_job[s][m][op] << " ";
				}
			
			}
			cout << endl;
		}
	
	}*/
	//变异offSpring
	for (int p=0;p<offSpring.size();++p) 
	{
		for (int s=0;s<Known::stage_count;++s) 
		{
			for (int j=0;j<Known::job_count;++j) 
			{
				if (rand()%3 ==0)
				{
					//当前操作需要变异
					Operation &op = offSpring[p].sol[s][j];
					//随机选择一个机器上的操作和当前操作交换
					int newma = rand() % Known::machine_count;
					if (offSpring[p].m_job[s][newma].size() == 0)
					{
						for (int o = op.po+1; o < offSpring[p].m_job[s][op.machine].size(); ++o)
						{
							jno = offSpring[p].m_job[s][op.machine][o];
							offSpring[p].sol[s][jno].po--;
						}
						offSpring[p].m_job[s][op.machine].erase(offSpring[p].m_job[s][op.machine].begin() + op.po);
						offSpring[p].m_job[s][newma].push_back(op.jno);
						op.po = 0;
						continue;
					}
					else 
					{
						int newop = rand() % offSpring[p].m_job[s][newma].size();
						jno = offSpring[p].m_job[s][newma][newop];
						//更新两个操作的  位置/机器号/机器上的工件号
						//cout << op.machine << "  " << op.po << " " << offSpring[p].m_job[s][op.machine].size() << " " << 6 << endl;
						offSpring[p].m_job[s][op.machine][op.po] = jno;

						offSpring[p].m_job[s][newma][newop] = j;

						offSpring[p].sol[s][jno].po = op.po;

						offSpring[p].sol[s][j].po = newop;

						offSpring[p].sol[s][jno].machine = op.machine;

						op.machine = newma;
					}
				}
				
			}
		}
		offSpring[p].update_sol();
		push_into_arch(offSpring[p]);
	}
	return;
}

/*非支配排序*/
vector<vector<int>> Fac::fast_non_dominated_sort() {
	vector<vector<int>> ans;//分层后的结果--个体所在种群的下标
	/*论文思路*/
		/*
		Sp_np.first = i的支配集合
		Sp_np.second = 支配i的个体counter
		*/
	vector<pair<vector<int>, int>> Sp_np = vector<pair<vector<int>, int>>(Rt.size()); //Sp_np.size() == Rt.size;
	vector<int> tmp;
	for (int i = 0; i < Rt.size(); ++i) { //遍历Rt，分别计算每个个体的Sp&&np
		Sp_np[i].second = 0;
		for (int j = 0; j < Rt.size(); ++j) {
			if (j != i) {
				if (Rt[i] < Rt[j]) { /*i支配j*/
					Sp_np[i].first.push_back(j);
				}
				else if (Rt[j] < Rt[i]) {
					Sp_np[i].second++;
				}
			}
		}
		if (Sp_np[i].second == 0) tmp.push_back(i);//此时tmp=第0层的元素下标

	}
	/*cout << "level 0:  ";
	for (int i0 = 0; i0 < tmp.size(); ++i0) {
		cout << tmp[i0] << " ";
	}
	cout << '\n';*/
	ans.push_back(tmp);
	for (int i = 0; i < ans.size(); ++i)
	{
		tmp.clear();//将要存储i+1层的元素下标
		for (auto p : ans[i])
		{
			for (auto q : Sp_np[p].first)
			{
				--Sp_np[q].second;
				if (Sp_np[q].second == 0)
					tmp.push_back(q);
			}
		}
		//bug：由于ans[i]的元素支配集合里会有重复的，导致tmp里面会有重复解
		/*			cout << "level " << i + 1 <<":	";
					for (int ii = 0; ii < tmp.size();++ii) {
						cout << tmp[ii] << " ";
					}*/
		if (tmp.size() != 0)
		{
			ans.push_back(tmp);
		}
	}
	//ans的每层会有重复解吗？
	return ans;

}


/*核心算法*/
void Fac::find_opt()
{
	Known::iter = 0;
	Known::st = clock();
	Known::et = clock();
	while (((double)(Known::et - Known::st) / CLOCKS_PER_SEC) < Known::itr_time)
	{
		/*根据当前种群，运用遗传算法生成子代*/
		for (int p = 0; p < populations.size(); ++p)
		{
			//每个粒子进行多代局部搜索
			for (int i = 0; i < Known::TS_iter; ++i)
			{
				populations[p].local_search(p);
				push_into_arch(populations[p]);
			}
		}

		//cout << "GA" << endl;
		GA(Known::iter);
		//cout << "NSGAII" << endl;
		/*运用NSGAII生成新的种群*/
		NSGAII();
		++Known::iter;
		Known::et = clock();
		for (int i=0;i<populations.size();++i) 
		{ 
			//cout << populations[i].tot_tardiness << "  " << populations[i].tot_energy_cost << endl;
			push_into_arch(populations[i]);
		}
		//cout << "第" << Known::iter << "次迭代:" << endl;
		
		if (Known::iter % Known::disturb_its ==0) 
		{
			//添加扰动
			for (int i=0;i<archive.size();++i) 
			{
				srand(clock()+rand()%1000+i);
				int index = rand() % populations.size();
				populations[index] = archive[i];
			}
		}
	}
}

void Individual::init_indi(int &rm)
{//随机生成一个个体
	for (int s = 0; s < Known::stage_count; ++s)
	{ //遍历每个阶段
		for (int j = 0; j < Known::job_count; ++j)
		{ //遍历当前阶段的每个工件
			//随机生成一个机器以及对应的速度
			sol[s][j].jno = j;
			sol[s][j].sno = s;
			srand(clock() + rand() % 1000 + rm);
			sol[s][j].machine = rand() % Known::machine_count;
			srand(clock()+rm);
			sol[s][j].m_speed = rand() % Known::pc[s].size(); //下标
			sol[s][j].po = m_job[s][sol[s][j].machine].size();
			m_job[s][sol[s][j].machine].push_back(j);
		}
	}
	return;
}

void Individual::_init_indi(int &rm)
{
	//有策略的生成初始解
	vector<pair<int, int>> jobs(Known::job_count);//jobs[i].first = 工件号；jobs[i].second = 上一道工序完成时间
//	vector<int> pre_et(known.job_count,0);
	int mac = 0, speed = 0, v = 0;
	for (int i = 0; i < jobs.size(); ++i)
	{
		jobs[i].first = i;
		jobs[i].second = 0;
	}
	sort(jobs.begin(), jobs.end(), [](pair<int, int> &j1, pair<int, int> &j2) {return Known::bpc[0][j1.first] < Known::bpc[0][j2.first]; }); //sort based basement processing time
	vector<int> mavt(Known::machine_count,0); //记录每台机器最早可用时间

	for (int j = 0; j < jobs.size(); ++j)
	{
		Operation tmp = Operation();
		tmp.sno = 0;
		tmp.jno = jobs[j].first;

		srand(clock() + rand()%100+rm);
		mac = rand() % Known::machine_count;
		srand(clock() + rm + rand() % 100);
		speed = rand() % Known::pc[0].size();//下标
		v = Known::pc[0][speed];

		tmp.machine = mac;
		tmp.m_speed = speed;

		tmp.po = m_job[0][mac].size();
		m_job[0][mac].push_back(tmp.jno);

		tmp.st= mavt[mac];
		tmp.process_time = Known::bpc[0][tmp.jno] + v;
		mavt[mac] = jobs[j].second = tmp.et = (tmp.st + tmp.process_time);

		tmp.energy_cost = tmp.update_energy_cost();
		sol[0][tmp.jno]=tmp;
	}

	for (int s = 1; s < Known::stage_count; ++s)
	{
		mavt = vector<int>(Known::machine_count,0);
		sort(jobs.begin(), jobs.end(), [](pair<int, int> &j1, pair<int, int> &j2) {return j1.second <  j2.second; }); //sort based pre_stage et
		for (int j = 0; j < jobs.size(); ++j)
		{
			Operation tmp = Operation();
			tmp.sno = s;
			tmp.jno = jobs[j].first;

			srand(clock() + rand() % 100 + rm);
			mac = rand() % Known::machine_count;
			srand(clock() + rm +rand()%100);
			speed = rand() % Known::pc[s].size();//下标
			v = Known::pc[s][speed];

			tmp.machine = mac;
			tmp.m_speed = speed;

			tmp.po = m_job[s][mac].size();
			m_job[s][mac].push_back(tmp.jno);

			tmp.st = max(jobs[j].second,mavt[mac]);
			tmp.process_time = Known::bpc[s][tmp.jno] + v;
			mavt[mac] = jobs[j].second = tmp.et = (tmp.st + tmp.process_time);

			tmp.energy_cost = tmp.update_energy_cost();
			sol[s][tmp.jno] = tmp;
		}
	}
	update_tot_tardiness();
	update_tot_energy_cost();
	return;
}

Individual::Individual()
{
	sol = vector<vector<Operation>>(Known::stage_count, vector<Operation>(Known::job_count));
	m_job = vector<vector<vector<int>>>(Known::stage_count, vector<vector<int>>(Known::machine_count));
	tot_tardiness = 0;
	tot_energy_cost = cd = 0.;
}

void Individual::_update_sol()
{ //主动调度
	vector<int> jobs(Known::job_count, 0);//jobs[i] = 工件i上一道工序完成时间 
	double cmp1 = 0., cmp2 = 0.;
	int jno = 0;
	//stage 1
	for (int s = 0; s < Known::stage_count; ++s)
	{
		vector<Idles> midles(Known::machine_count);//当前阶段每台机器的空闲时间段集合
		//midles[i].opjno = 已经确定顺序的机器编号
		for (int m = 0; m < Known::machine_count; ++m)
		{
			vector<int> &ele_ = m_job[s][m]; //当前机器的操作集合
			for (int op = 0; op < ele_.size(); ++op)
			{
				jno = ele_[op];//当前编号
				sol[s][jno].process_time = Known::bpc[s][jno] + Known::pc[s][sol[s][jno].m_speed];
				if (!midles[m].ava_time)
				{
					//是当前机器处理的第一个操作
					sol[s][jno].st = max(midles[m].ava_time, jobs[jno]);
					jobs[jno] = midles[m].ava_time = sol[s][jno].et = sol[s][jno].st + sol[s][jno].process_time;
					sol[s][jno].energy_cost = sol[s][jno].update_energy_cost();
					midles[m].opjno.push_back(jno);
				}
				else
				{
					bool flag = false;  //判断当前操作是否插入了前面空闲时间
					//没有可用空闲段
					sol[s][jno].st = max(midles[m].ava_time, jobs[jno]);
					sol[s][jno].et = sol[s][jno].st + sol[s][jno].process_time;
					cmp2 = sol[s][jno].energy_cost = sol[s][jno].update_energy_cost(); //正常顺序的能耗成本

					for (int i = 1; i < midles[m].opjno.size(); ++i)
					{//遍历机器的空闲时间段
						int jno1 = midles[m].opjno[i-1];
						int jno2 = midles[m].opjno[i];
						if (max(sol[s][jno1].et, jobs[jno]) + sol[s][jno].process_time <= sol[s][jno2].st)
						{//可以插入，但还要判断插入会不会使能耗变小
							Operation tmpop = sol[s][jno];
							tmpop.st = max(sol[s][jno1].et, jobs[jno]);
							tmpop.et = tmpop.st + sol[s][jno].process_time;
							cmp1 = tmpop.update_energy_cost();
							if (cmp1 <= cmp2)
							{
								//可以插入
								midles[m].opjno.insert(midles[m].opjno.begin()+i,jno);
								sol[s][jno].st = tmpop.st;
								sol[s][jno].et = tmpop.et;
								sol[s][jno].energy_cost = cmp1;
								
								jobs[jno]  = sol[s][jno].et;
								flag = true;
								break;
							}
							else continue;
						}
					}
					if (!flag)
					{
						//没有可用空闲段
						jobs[jno] = midles[m].ava_time = sol[s][jno].et;
						midles[m].opjno.push_back(jno);
					}
				}
			}
			m_job[s][m] = midles[m].opjno;
			for (int i = 0; i < midles[m].opjno.size(); ++i)
			{
				jno = midles[m].opjno[i];
				sol[s][jno].po = i;
			}
		}
	}
	update_tot_tardiness();
	update_tot_energy_cost();
	return;
}

void Individual::update_sol()
{
	int jno = 0;
	vector<int> jpre_et(Known::job_count, 0); //工件上一道工序的结束时间
	for (int s=0;s<Known::stage_count;++s) 
	{
		//vector<int>	mavt(Known::machine_count,0); //机器最早可用时间
		for (int m=0;m<Known::machine_count;++m) 
		{
			int avt = 0;
			for (int op=0;op<m_job[s][m].size();++op) 
			{
				jno = m_job[s][m][op];
				sol[s][jno].st = max(jpre_et[jno],avt);
				sol[s][jno].process_time = Known::bpc[s][jno] + Known::pc[s][sol[s][jno].m_speed];
				jpre_et[jno] = avt = sol[s][jno].et = sol[s][jno].st + sol[s][jno].process_time;
				sol[s][jno].update_energy_cost();
			}
		}
	}
	update_tot_tardiness();
	update_tot_energy_cost();
	return ;
}

void Individual::update_tot_tardiness()
{
	tot_tardiness = 0;
	for (int j = 0; j < Known::job_count; ++j)
	{
		//每个工件
		tot_tardiness += max(0, (sol[Known::stage_count - 1][j].et - Known::due_date[j]));
	}
}


void Individual::update_tot_energy_cost()
{
	tot_energy_cost = 0.;
	for (int s = 0; s < Known::stage_count; ++s)
	{
		for (int j = 0; j < Known::job_count; ++j)
		{
			tot_energy_cost += sol[s][j].energy_cost;
		}
	}
}

void Individual::update_tbb(TabuOp & tbo, int & pos)
{
	int tbov = (tbo.mno * 10 + tbo.mv) * 100 + tbo.mpo;
	Known::tabuTable[pos][tbo.sno][tbo.jno][tbov] = Known::iter + rand() % 10;
}

bool Individual::operator<(Individual &b) const
{
	if (tot_tardiness <= b.tot_tardiness && tot_energy_cost < b.tot_energy_cost || tot_tardiness < b.tot_tardiness && tot_energy_cost <= b.tot_energy_cost) return true;
	return false;
}

bool Individual::operator==(Individual & b) const
{
	if (tot_tardiness < b.tot_tardiness && tot_energy_cost > b.tot_energy_cost || tot_tardiness > b.tot_tardiness && tot_energy_cost < b.tot_energy_cost) return true;
	return false;
}

bool Individual::check()
{
	bool flag = true;
	for (int s = 0; s < Known::stage_count; ++s)
	{
		vector<vector<int>> seqs(Known::machine_count);
		set<int> ops; //判断是否每个工件都被处理了
		for (int j = 0; j < Known::job_count; ++j)
		{
			ops.insert(sol[s][j].jno);
			seqs[sol[s][j].machine].push_back(sol[s][j].et);
			seqs[sol[s][j].machine].push_back(sol[s][j].et);
			if (sol[s][j].st < (s == 0 ? 0 : sol[s - 1][j].et)) return false;
		}
		for (const auto &ma : seqs)
		{
			for (int si = 2; si < ((int)ma.size() - 1); si += 2)
			{
				if (ma[si] < ma[si - 1])  return false;
			}
		}

		if (ops.size() != Known::job_count) return false;
	}
	return flag;
}

/*void Individual::assign(const Individual &indi, Factory &f) {
	for (int s = 0; s < Known::stage_count; ++s) {
		for (int j = 0; j < Known::job_count; ++j) {
			this->sol[s][j].m_v = indi.sol[s][j].m_v;
		}
	}
}*/

void Individual::cal_machine_job()
{
	for (int s = 0; s < Known::stage_count; ++s)
	{
		for (int j = 0; j < Known::job_count; ++j)
		{
			m_job[s][sol[s][j].machine].push_back(j);
		}
	}
}

void Individual::print_sol()
{
	for (int s = 0; s < Known::stage_count; ++s)
	{
		vector<vector<pair<int, int>>> ma(Known::machine_count);
		for (int j = 0; j < Known::job_count; ++j)
		{
			Operation &op = sol[s][j];
			pair<int, int> jv;
			jv.first = op.jno;
			jv.second = op.m_speed;
			ma[op.machine].push_back(jv);
		}

		for (auto &ele : ma)
		{
			cout << "[	";
			for (int i = 0; i < ele.size(); ++i)
			{
				cout << ele[i].first << ',' << ele[i].second << "	";
			}
			cout << "]";
		}
		cout << endl;
	}
}

void Individual::print_sol_()
{
	cout << "个体	" << tot_tardiness << "	" << tot_energy_cost << endl;
	for (int s = 0; s < Known::stage_count; ++s)
	{
		cout << "stage" << s + 1 << endl;
		for (int j = 0; j < Known::job_count; ++j)
		{
			//输出选择的 机器&&处理时间
			cout << sol[s][j].machine << "--" << sol[s][j].m_speed << "	";
		}
		cout << endl;
	}
}

void Fac::start_down()
{
	cout << "...最终结果..." << endl;
	/*cout << "final populations: " << endl;
	sort(populations.begin(), populations.end(), [](Individual &i1, Individual &i2) {return i1.tot_tardiness < i2.tot_tardiness; });
	for (int i = 0; i < populations.size(); ++i)
	{
		cout << i << ":		" << populations[i].tot_tardiness << "		" << populations[i].tot_energy_cost << endl;
	}*/
	//cout << "final archive: " << endl;
	sort(archive.begin(), archive.end(), [](Individual &i1, Individual &i2) {return i1.tot_tardiness < i2.tot_tardiness; });
	/*for (int i = 0; i < archive.size(); ++i)
	{
		cout << i << ":		" << archive[i].tot_tardiness << "		" << archive[i].tot_energy_cost << endl;
        //输出工件的st&&et
		Individual &indi = populations[i];
		for (int s=0;s<Known::stage_count;++s) 
		{
			for (int j=0;j<Known::job_count;++j) 
			{
				cout << indi.sol[s][j].st << "， " << indi.sol[s][j].et << "  ";
			}
			cout << endl;
		}
		//输出每台机器上工件的st,et
		for (int s=0;s<Known::stage_count;++s) 
		{
			for (int m=0;m<Known::machine_count;++m) 
			{
				vector<int> &e = indi.m_job[s][m];
				cout << "[ ";
				for (int o=0;o<e.size();++o) 
				{
					cout << indi.sol[s][e[o]].st << "， " << indi.sol[s][e[o]].et<<" ";
				}
				cout << " ]";
			}
			cout << endl;
		}
	}
	*/
	cout << "...开始写文件..." << endl;
	writesol();
	cout << "--程序结束--" << endl;
}


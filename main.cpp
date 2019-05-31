#include "MOPSO.h"
#include "NSGAII.h"

int main(int argc, char*  argv[])
{
	map<string, string> &wl = Known::pars_wl;
	//srand(time(0));
	//windows环境参数设置
	wl["genes"] = "300";
	wl["pops"] = "50"; //NSGAII默认100代
	wl["TS_iter"] = "5";
	wl["disturb_its"] = "10";
	wl["acp"] = "0.0";//接受互不支配解的概率
	wl["acp1"] = "0.0";
	wl["tb_acp"] = "0.0"; 
	wl["tb_acp1"] = "0.0";
	wl["w"] = "0.0"; //变异概率
	wl["c1"] = "0.0";//跟pbest的交叉概率
	wl["c2"] = "0.0";//跟gbest的交叉概率
	wl["al_name"] = "_HPSO_";//当前算法名称
	wl["ifn"] = "./Instances/100_10_8_4_10_3.txt";
	wl["ofn2"] = "./sols_no_info/";  //用于matlab分析的解
	wl["ofn"] = "./sols/";  //用于表格分析的解集合
	wl["fn"] = "100_10_8_4_10_3.txt";
	wl["run_num"] = "1"; //运行次数
	wl["Pc"] = "0.1";  //操作变异概率
	wl["Pm"] = "0.1";
	wl["itr_time"] = "3600"; //运行10min

	//linux环境下的参数赋值
	for (int i = 1; i < argc; i += 2)
	{
		if (i + 1 < argc)
			wl[argv[i]] = argv[i + 1];
		else cout << "main函数参数不成对！" << endl;
	}

	if (wl["al_name"] == "_HPSO_" || wl["al_name"] == "_HPSOI_" ||wl["al_name"] == "_HPSOII_") 
	{
		//基于解禁忌的PSO算法
		Factory fac = Factory();
		fac.control();
	}
	else if (wl["al_name"] == "_NSGAII_") 
	{
		//NSGAII算法
		Fac fac = Fac();
		fac.start_up();
	}
	

#ifdef _WIN32
	system("pause");
#endif
}
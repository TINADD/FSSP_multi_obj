#include "MOPSO.h"
#include "NSGAII.h"

int main(int argc, char*  argv[])
{
	map<string, string> &wl = Known::pars_wl;
	//srand(time(0));
	//windows������������
	wl["genes"] = "300";
	wl["pops"] = "1"; //NSGAIIĬ��100��
	wl["TS_iter"] = "50";
	wl["disturb_its"] = "1";
	wl["acp"] = "0.0";//���ܻ���֧���ĸ���
	wl["acp1"] = "0.0";
	wl["tb_acp"] = "0.0"; 
	wl["tb_acp1"] = "0.0";
	wl["w"] = "1.0"; //�������
	wl["c1"] = "0.5";//��pbest�Ľ������
	wl["c2"] = "0.8";//��gbest�Ľ������
	wl["al_name"] = "_HPSO_";//��ǰ�㷨����
	wl["ifn"] = "./Instances/100_10_8_4_10_3.txt";
	wl["ofn2"] = "./newsols/";  //����matlab�����Ľ�
	wl["ofn"] = "./Solutions/";  //���ڱ������Ľ⼯��
	wl["fn"] = "100_10_8_4_10_3.txt";
	wl["run_num"] = "1"; //���д���
	wl["Pc"] = "0.1";  //�����������
	wl["Pm"] = "0.1";
	wl["itr_time"] = "3600"; //����10min

	//linux�����µĲ�����ֵ
	for (int i = 1; i < argc; i += 2)
	{
		if (i + 1 < argc)
			wl[argv[i]] = argv[i + 1];
		else cout << "main�����������ɶԣ�" << endl;
	}

	if (wl["al_name"] == "_HPSO_" || wl["al_name"] == "_HPSOI_" ||wl["al_name"] == "_HPSOII_") 
	{
		//���ڽ���ɵ�PSO�㷨
		Factory fac = Factory();
		fac.control();
	}
	else if (wl["al_name"] == "_NSGAII_") 
	{
		//NSGAII�㷨
		Fac fac = Fac();
		fac.start_up();
	}
	

#ifdef _WIN32
	system("pause");
#endif
}
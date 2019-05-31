#include "share.h"
map<string, string> Known::pars_wl;

int Known::job_count = 0;
int Known::stage_count = 0;
int Known::machine_count = 0;
vector<int> Known::due_date;
vector<double> Known::ep = vector<double>(24);
vector<vector<int>> Known::bpc; //基本处理时间
vector<vector<int>> Known::pc; //可选处理时间
vector<vector<double>> Known::bec; //基本能耗
vector<pair<vector<bool>, vector<bool>>> Known::hs_1;
vector<pair<vector<bool>, vector<bool>>> Known::hs_2;
vector<vector<vector<map<int, int>>>> Known::tabuTable;
double Known::runtime;
double Known::itr_time;
unsigned long Known::next = 1;


time_t Known::st;
time_t Known::et;
time_t Known::arch_et;

int Known::iter = 0;
int Known::genes;
int Known::pops;
int Known::TS_iter;
int Known::disturb_its;
double Known::acp;
double Known::tb_acp;
double Known::w;
double Known::c1;
double Known::c2;
double Known::Pc;
double Known::Pm;

void Known::read_file()
{
	string s1, s2;
	ifstream infile(pars_wl["ifn"]);
	if (!infile)
	{
		std::cout << "Can't open the file!" << infile.is_open() << std::endl;
		return;
	}
	infile >> s1 >> s1 >> s1 >> job_count;
	infile >> s1 >> s1 >> s1 >> stage_count;
	infile >> s1 >> s1 >> s1 >> machine_count;
	bpc = vector<vector<int>>(stage_count, vector<int>(job_count));
	bec = vector<vector<double>>(stage_count, vector<double>(job_count));
	//电价
	for (int i = 0; i < 7; ++i) {
		ep[i] = 80.;
	}
	for (int i = 7; i < 15; ++i) {
		ep[i] = 160.;
	}
	for (int i = 15; i < 20; ++i) {
		ep[i] = 240.;
	}
	for (int i = 20; i < 22; ++i) {
		ep[i] = 160.;
	}
	for (int i = 22; i < 24; ++i) {
		ep[i] = 80.;

	}
	due_date = vector<int>(job_count);
	pc = vector<vector<int>>(stage_count);
	getline(infile, s1);
	getline(infile, s1);
	/*读入baseline processing time*/
	for (int i = 0; i < stage_count; ++i) {
		infile >> s1 >> s2;
		for (int j = 0; j < job_count; ++j) {
			infile >> bpc[i][j];
		}
	}
	/*读入能耗*/
	getline(infile, s1);
	getline(infile, s1);
	for (int i = 0; i < stage_count; ++i) {
		infile >> s1 >> s1;
		for (int j = 0; j < job_count; ++j) {
			infile >> bec[i][j];
		}
	}
	/*读入dute_date*/
	infile >> s1 >> s1;
	for (int i = 0; i < job_count; ++i) {
		infile >> due_date[i];
	}

	getline(infile, s1);
	getline(infile, s1);
	/*读入处理速度*/
	for (int i = 0; i < stage_count; ++i) {
		getline(infile, s1);
		std::stringstream ss(s1);
		ss >> s2 >> s2;
		int speed = 0;
		while (ss >> speed) {
			pc[i].push_back(speed);
		}

	}
	infile.close();
}

void Known::print_file()
{
	cout << "number of jobs:	" << job_count << endl;
	cout << "number of stages:	" << stage_count << endl;
	cout << "number of machines:	" << machine_count << endl;
	cout << "baseline processing time:	" << '\n';
	for (int s = 0; s < stage_count; ++s) {
		cout << "stage	" << s + 1 << ":	";
		for (int j = 0; j < job_count; ++j) {
			cout << "	" << bpc[s][j];
		}
		cout << '\n';
	}

	cout << "baseline energy consumption:" << '\n';
	for (int s = 0; s < stage_count; ++s) {
		cout << "stage	" << s + 1 << ":	";
		for (int j = 0; j < job_count; ++j) {
			cout << "	" << bec[s][j];
		}
		cout << '\n';
	}
	cout << "due date:	";
	for (int j = 0; j < job_count; ++j) {
		cout << "	" << due_date[j];
	}
	cout << '\n' << "processing speed:" << '\n';
	for (int s = 0; s < stage_count; ++s) {
		cout << "stage	" << s + 1 << ":	";
		for (int j = 0; j < pc[s].size(); ++j) {
			cout << "	" << pc[s][j];
		}
		cout << '\n';
	}
}

void Known::write_head()
{
	pars_wl["ofn"] += (pars_wl["job_id"] + pars_wl["al_name"] + pars_wl["fn"]);
	ofstream ofile(pars_wl["ofn"], ios::app); //以追加的方式写入
	for (auto pa : pars_wl)
	{
		ofile << pa.first << ":" << pa.second << '\t';
	}
	ofile << endl;
	ofile << "run_num   " << '\t' << "iter" << '\t' << "run_time" << '\t' << "num_archive" << endl;
	ofile.close();
	genes = stoi(pars_wl["genes"]);
	pops = stoi(pars_wl["pops"]);
	TS_iter = stoi(pars_wl["TS_iter"]);
	disturb_its = stoi(pars_wl["disturb_its"]);
	acp = stod(pars_wl["acp"]);
	w = stod(pars_wl["w"]);
	c1 = stod(pars_wl["c1"]);
	c2 = stod(pars_wl["c2"]);
	Pc = stod(pars_wl["Pc"]);
	Pm = stod(pars_wl["Pm"]);
	itr_time = stod(pars_wl["itr_time"]);
}

int Known::myrand(void)
{
	//next = next * 1103515245 + 12345;
	//next = next * 10 + 12345;
    //Known kno = Known();
	//kno.end();
	//next = rand() % 100;
	//cout << Known::interval << endl;
	//int seed = Known::interval*pow(10, 6);
	//srand(clock()+rand()%100+rand()%10);
	//int seed = ((rand() % 1000) * 10 + rand() % 10) * 100 + rand() % 100;
	//srand(time(0)+ rand()%100 + rand()%10);
	//srand(seed+clock());
	//return((unsigned)(next / 65536) % 32768);
	return rand();
}

void Known::mysrand(unsigned seed)
{
	next = seed;
}


Operation::Operation()
{
	jno = sno = machine = 0;
	m_speed = 0;
	st = et = process_time = 0;
	energy_cost = 0.;
}

/*计算当前操作的能耗*/
double Operation::update_energy_cost()
{
	/*计算当前操作能耗成本*/
	double v = Known::pc[sno][m_speed];
	double EP = (1. + 0.6*(v*v) / (Known::bpc[sno][jno] * Known::bpc[sno][jno]) - 1.4 * v / Known::bpc[sno][jno])
		*Known::bec[sno][jno] * Known::bpc[sno][jno] / (v + Known::bpc[sno][jno]);
	double final_ec = 0.;
	for (int i = st; i < et; ++i)
	{
		final_ec += EP * Known::ep[i % 24];//越界了 should %24
	}
	energy_cost = final_ec;
	return final_ec;
}

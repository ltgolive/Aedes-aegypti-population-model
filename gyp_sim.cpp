//gyp_sim.cpp was developed by Dr Penelope A. Hancock
//The development of gyp_sim.cpp is described in Hancock et al. 2016, "Predicting Wolbachia
//invasion dynamics in Aedes aegypti populations using models of density-dependent
//demographic traits", BMC Biology.
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <float.h>
#include <stdio.h>
#include <string>
#include <sstream>
#include <iomanip>

using namespace std;

typedef vector<int> Row_Int;
typedef vector<Row_Int> Matrix_Int;

typedef vector<double> Row_Double;
typedef vector<Row_Double> Matrix_Double;

//some of these vectors that are getting passed around could just be global...
void read_data(ofstream& cohort_means_out, ofstream& cohort_means_wolb_out, ofstream& cohort_stds_out, ofstream& cohort_stds_wolb_out, ofstream& mu_p_out,ofstream& mu_p_wolb_out, ofstream& L_out, ofstream& L_wolb_out, ofstream& A_out, ofstream& A_wolb_out, ofstream& freqA2_out, ofstream& lambda_out, ifstream& rel_siz_in);
double Run_Model(const vector<int> max_dt, ofstream& mu_p_out,  ofstream& freqA2_out, ofstream& lambda_out);

const int no_cohorts =1000;
const int no_pdays = 170;
const int no_sdays = 169;
const int maxtime = 2500;

vector<float> mean_dt(no_cohorts), mean_dt_wolb(no_cohorts), std_dt(no_cohorts), std_dt_wolb(no_cohorts), max_dt(no_cohorts), mu_p(maxtime), mu_p_wolb(maxtime), L(maxtime),  L_wolb(maxtime), A(maxtime), A_ovipos(maxtime), A_wolb(maxtime), A_ovipos_wolb(maxtime), L_data(no_pdays);

Matrix_Double P(maxtime, Row_Double(3)), P_wolb(maxtime, Row_Double(3));

vector<int> pupae(no_pdays),pdate(no_pdays),hatch(no_cohorts),hdate(no_cohorts),larvae_sampled(no_cohorts), pupae_sampled(no_pdays);

stringstream strstm;

int release_day;
double release_size, DI_cost, DI_L_cost;

int main(){

//Input and Output files

ofstream cohort_means_out,cohort_means_wolb_out, cohort_stds_out, cohort_stds_wolb_out, mu_p_out, mu_p_wolb_out, A_out, A_wolb_out, L_out, L_wolb_out, freqA2_out, lambda_out;

ifstream rel_siz_in;

read_data(cohort_means_out,cohort_means_wolb_out, cohort_stds_out, cohort_stds_wolb_out, mu_p_out, mu_p_wolb_out, L_out, L_wolb_out, A_out, A_wolb_out, freqA2_out, lambda_out, rel_siz_in);

//Assign Cohorts to blocks and maximum development time
vector<int> max_dt(no_cohorts);
for (int i=0; i<no_cohorts; i++){
	max_dt.at(i) = maxtime - hdate.at(i);
}

Run_Model(max_dt, mu_p_out, freqA2_out, lambda_out);

for (int j=0; j<no_cohorts; j++) {
	cohort_means_out << mean_dt.at(j) << " ";
	cohort_means_wolb_out << mean_dt_wolb.at(j) << " ";
	cohort_stds_out << std_dt.at(j) << " ";
	cohort_stds_wolb_out << std_dt_wolb.at(j) << " ";
}
cohort_means_out << endl;
cohort_means_wolb_out << endl;
cohort_stds_out << endl;
cohort_stds_wolb_out << endl;

for (int j=0; j<maxtime; j++) {
	mu_p_out << mu_p.at(j) << " ";
	mu_p_wolb_out << mu_p_wolb.at(j) << " ";
}
mu_p_out << endl;
mu_p_wolb_out << endl;

for (int j=0; j<maxtime; j++) {
	L_out << L.at(j) << " ";
	L_wolb_out << L_wolb.at(j) << " ";
	A_out << A.at(j) << " ";
	A_wolb_out << A_wolb.at(j) << " ";
}
L_out << endl;
A_out << endl;


return 0;

}

//---------------------------------function "read_data"---------------------------------------
void read_data(ofstream& cohort_means_out, ofstream& cohort_means_wolb_out, ofstream& cohort_stds_out,ofstream& cohort_stds_wolb_out, ofstream& mu_p_out, ofstream& mu_p_wolb_out, ofstream& L_out, ofstream& L_wolb_out, ofstream& A_out, ofstream& A_wolb_out, ofstream& freqA2_out, ofstream& lambda_out, ifstream& rel_siz_in){


ifstream hdate_data, hatch_data;
hdate_data.open("hdate3.txt");
for (int i=0; i<no_cohorts; i++){
hdate_data >> hdate.at(i);
}

//Open output files named in inits file
string cohort_means_file,cohort_means_wolb_file, cohort_stds_file, cohort_stds_wolb_file, mu_p_file, mu_p_wolb_file, L_file, L_wolb_file, A_file, A_wolb_file, rel_siz_file, freqA2_file, freqL_file, lambda_file;
cin >> cohort_means_file >> cohort_means_wolb_file >> cohort_stds_file >> cohort_stds_wolb_file  >> mu_p_file >> L_file >> A_file >> mu_p_wolb_file >> L_wolb_file>> A_wolb_file >>  freqA2_file >>lambda_file >> rel_siz_file;
cout <<   " cohort_means_file " << cohort_means_file<< " cohort_means_wolb_file " << cohort_means_wolb_file << " cohort_stds_file " << cohort_stds_file << " cohort_stds_wolb_file " << cohort_stds_wolb_file <<  " mu_p_file " << mu_p_file << " L_file " << L_file << " A_file "<<  A_file  << " mu_p_wolb_file " << mu_p_wolb_file << " L_wolb_file " << L_wolb_file << " A_wolb_file " << A_wolb_file << " release_siz_file " << rel_siz_file << " freqA2_file " << freqA2_file << " lambda_file " << lambda_file << endl;

strstm.clear();
strstm.str("");
strstm <<cohort_means_file;
string filename4=strstm.str();
cohort_means_out.open(filename4.c_str());

strstm.clear();
strstm.str("");
strstm <<cohort_stds_file;
string filename5=strstm.str();
cohort_stds_out.open(filename5.c_str());

strstm.clear();
strstm.str("");
strstm <<mu_p_file;
string filename10=strstm.str();
mu_p_out.open(filename10.c_str());

strstm.clear();
strstm.str("");
strstm <<L_file;
string filename11=strstm.str();
L_out.open(filename11.c_str());

strstm.clear();
strstm.str("");
strstm <<A_file;
string filename17=strstm.str();
A_out.open(filename17.c_str());

strstm.clear();
strstm.str("");
strstm <<mu_p_wolb_file;
string filename18=strstm.str();
mu_p_wolb_out.open(filename18.c_str());

strstm.clear();
strstm.str("");
strstm <<L_wolb_file;
string filename19=strstm.str();
L_wolb_out.open(filename19.c_str());

strstm.clear();
strstm.str("");
strstm <<A_wolb_file;
string filename20=strstm.str();
A_wolb_out.open(filename20.c_str());

strstm.clear();
strstm.str("");
strstm <<cohort_means_wolb_file;
string filename21=strstm.str();
cohort_means_wolb_out.open(filename21.c_str());

strstm.clear();
strstm.str("");
strstm <<cohort_stds_wolb_file;
string filename22=strstm.str();
cohort_stds_wolb_out.open(filename22.c_str());

strstm.clear();
strstm.str("");
strstm <<rel_siz_file;
string filename24=strstm.str();
rel_siz_in.open(filename24.c_str());

strstm.clear();
strstm.str("");
strstm <<freqA2_file;
string filename25=strstm.str();
freqA2_out.open(filename25.c_str());

strstm.clear();
strstm.str("");
strstm <<lambda_file;
string filename26=strstm.str();
lambda_out.open(filename26.c_str());

cin >> release_day;
rel_siz_in >> release_size;
cin >> DI_cost;
cin >> DI_L_cost;
cout << "DI_cost " << DI_cost << " (additional density-INdependent daily mortality experienced by adults in the field environment) " << endl;
cout << "DI_L_cost " << DI_L_cost << " (additional density-INdependent daily mortality experienced by larvae in the field environment) " << endl;

}

//---------------------------------function "Run_Model"----------------------------------------
double Run_Model(const vector<int> max_dt, ofstream& mu_p_out, ofstream& freqA2_out, ofstream& lambda_out){

vector<double>  non_emerg_prob(no_cohorts), non_emerg_prob_wolb(no_cohorts), L_avg_cohort(no_cohorts), hatch_sim(no_cohorts), hatch_sim_wolb(no_cohorts), A_wolb_imm(maxtime), A_ovipos_wolb_imm(maxtime), mn_dt_gam(no_cohorts), mn_dt_wolb_gam(no_cohorts), std_dt_gam(no_cohorts), std_dt_wolb_gam(no_cohorts),no_emerg_tot(no_cohorts), no_emerg_tot_wolb(no_cohorts), L_avg(no_cohorts), lambda(no_cohorts), lambda_wolb(no_cohorts);
vector<int> emerg_flag(no_cohorts);
//Make an array to store the number of larvae of each age at each time step
Matrix_Double L_cohort(maxtime, Row_Double(no_cohorts));
Matrix_Double L_cohort_wolb(maxtime, Row_Double(no_cohorts));
//Make an array to store the number of individuals from each cohort that emerge as pupae at each time 
Matrix_Double emerg_record(maxtime, Row_Double(no_cohorts)), emerg_record_wolb(maxtime, Row_Double(no_cohorts));
double tau_p, surv_L, L_avg1, L_cum2, L_cum, denom2, Dt_shp, Dt_shp_wolb, Dt_scl, Dt_scl_wolb, prob, prob_wolb, no_emerge,no_emerge_wolb, H_cum, av_growth, L_av1, L_av2, L_avg_f, sh, w, freqA, lambda_max, survA_imm, intc_TL, intc_TL_wolb, alpha_TL, alpha_TL_wolb, exp_TL, exp_TL_wolb, a_value, a_value_wolb, b_value, b_value_wolb, b_value2, b_value2_wolb, intc_f, alpha_f, survA, survA_wolb, lambda1, lambda_wolb1, freq_end, lambda_min, mean_max,std_max,std_min, freq_end2, freqL_end,mean_dt1,mean_lambda,mean_std_dt;
int index,first_hatch_date, max_cohort, L_max, time_lag, count, tlagh, tlagg,tlagi,tlagp, min_cohort, max_dt_int, release_end, release_day_init, dens_lag;

first_hatch_date=hdate.at(0);
tlagh = 5;//lag between oviposition and hatching
tlagg = 6;//minimum time between emergence and first oviposition
tlagi = 21;//time over which larval density is averaged
tlagp = 2;//time required for pupal development
intc_f=28.0;
alpha_f=-3.3;
intc_TL = 1.8;
intc_TL_wolb = 1.8;
alpha_TL = 0.536;
alpha_TL_wolb = 0.536;
exp_TL = 0.533;
exp_TL_wolb = 0.533;
a_value = 0.22;
a_value_wolb = 0.22;
b_value = 0.0168;
b_value_wolb = 0.0168;
b_value2 = 0.867;
b_value2_wolb = 0.867;
max_dt_int=100;
survA = 1-0.03-DI_cost;
lambda_max=14;
lambda_min=0.5;
mean_max = 60;
std_max = 40;
std_min = 1.0;
L_max = 4200;


//Wolbachia parameters
sh=0.99;
w=0.01;

release_day_init = release_day;
release_end = release_day + 89;

survA_imm=survA;
survA_wolb = survA;

//Initialise L_cohort
for (int i=0; i<maxtime; i++){
	for (int j=0; j< no_cohorts; j++){
		L_cohort[i][j]=0;
		L_cohort_wolb[i][j]=0;
	}
}

//Initialise L
for (int i=0; i<maxtime; i++) {L.at(i)=0; L_wolb.at(i)=0;}


//Initialise emerg_record array.  
for (int i=0; i<maxtime; i++){
	for (int j=0; j<no_cohorts; j++){
		emerg_record[i][j]=0;
		emerg_record_wolb[i][j]=0;
	}
}

//Initialise non_emerg_prob array.  
for (int i=0; i<no_cohorts; i++){
	non_emerg_prob.at(i)=1;
	non_emerg_prob_wolb.at(i)=1;
}

//Initialise emerg_flag
for (int i=0; i<no_cohorts; i++) emerg_flag.at(i)=0;

//Initialise L_avg_cohort & L_avg
for (int i=0; i<no_cohorts; i++) {L_avg_cohort.at(i) = 0; L_avg.at(i)=0;}

//Clear mu_p
for (int i=0; i<maxtime; i++) {mu_p.at(i) = 0; mu_p_wolb.at(i)=0;}

//Clear mn_dt_gam
for (int i=0; i<no_cohorts; i++) {mn_dt_gam.at(i) = 999; mn_dt_wolb_gam.at(i) = 999;}

//Clear std_dt
for (int i=0; i<no_cohorts; i++) {std_dt_gam.at(i) = 1; std_dt_wolb_gam.at(i)=1;}

for (int i=0; i<no_cohorts; i++) {
	mean_dt.at(i)=0;
	std_dt.at(i)=0;
	no_emerg_tot.at(i)=0;
	no_emerg_tot_wolb.at(i)=0;
}

//Initialize P and P_wolb
for (int i=0; i<maxtime; i++){
	for (int j=0; j<3; j++){
		P[i][j]=0;
		P_wolb[i][j]=0;
	}
}

//Initialize A
A.at(0) = 2;
for (int i=1; i<maxtime; i++) A.at(i)=0;

//Initialize A_wolb
A_wolb.at(0) = 0;
for (int i=1; i<maxtime; i++) A_wolb.at(i)=0;

//Initialize A_ovipos and A_ovipos_wolb
for (int i=0; i<maxtime; i++) {
	A_ovipos.at(i) = 0;
	A_ovipos_wolb.at(i) = 0;
	A_wolb_imm.at(i) = 0;
	A_ovipos_wolb_imm.at(i) = 0;
}

//Initialize hatch_sim and hatch_sim_wolb
for (int i=0; i<no_cohorts; i++) {
	hatch_sim.at(i) = 0;
	hatch_sim_wolb.at(i) = 0;
}

//Main time loop
for (int time1=1; time1<maxtime; time1++){

//-------------------LARVAE--------------------------------------------------

	surv_L = 0.95 - DI_L_cost;

	max_cohort=-1;
	for (int i=0; i<no_cohorts; i++) if (hdate.at(i)<=time1) max_cohort = i;
	min_cohort=no_cohorts-1;
	for (int i=no_cohorts-1; i>=0; i--) if (hdate.at(i)>=time1 - max_dt_int) min_cohort = i;		

	for (int cohort=min_cohort; cohort<=max_cohort; cohort++){

		L_cohort[time1][cohort] = L_cohort[time1-1][cohort] *  surv_L;
		L_cohort_wolb[time1][cohort] = L_cohort_wolb[time1-1][cohort] * surv_L;
		if (time1 == hdate.at(cohort)){
//One week larval density lagged by 6 days
			L_avg.at(cohort) = 0;

			dens_lag = time1-tlagh-tlagi-tlagg;
			if (dens_lag>=hdate.at(1)){
			//lags: between oviposition and hatching, time over which larval density is averaged, time to become gravid after pupation
				for (int j=time1-tlagh-tlagi-tlagg; j<time1-tlagh-tlagg; j++) L_avg.at(cohort)+=L.at(j) + L_wolb.at(j);
					L_avg.at(cohort) = L_avg.at(cohort)/tlagi; 
			}
			else if (dens_lag>=0) {
				if (time1-tlagh-tlagg<=hdate.at(1)) L_avg.at(cohort)=60;//if the end of the lag is before the first hatch date
				if (time1-tlagh-tlagg>hdate.at(1)) {
					for (int j=dens_lag; j<=hdate.at(1); j++) L_avg.at(cohort)+=60;
					for (int j=hdate.at(1)+1; j<=time1-tlagh-tlagg; j++) L_avg.at(cohort)+=L.at(j) + L_wolb.at(j);
					L_avg.at(cohort) = L_avg.at(cohort)/tlagi;
					
				} 
			}
			else L_avg.at(cohort) = 60; 

			L_avg_f = L_avg.at(cohort);

			lambda1 = intc_f + alpha_f*log(L_avg_f);

			if (lambda1<lambda_min) lambda1=lambda_min;
			if (lambda1>lambda_max) lambda1 = lambda_max;
			lambda_wolb1 = lambda1;

			lambda.at(cohort) = lambda1;
			lambda_wolb.at(cohort) = lambda_wolb1;

 			if (time1>tlagh+tlagg){
				hatch_sim.at(cohort) = lambda1*A_ovipos.at(time1-tlagh)+ w*lambda_wolb1*A_ovipos_wolb.at(time1-tlagh)+ w*lambda_max*A_ovipos_wolb_imm.at(time1-tlagh);
				hatch_sim_wolb.at(cohort) = lambda_wolb1*(1-w)*A_ovipos_wolb.at(time1-tlagh) + lambda_max*(1-w)*A_ovipos_wolb_imm.at(time1-tlagh);
			}

			L_cohort[time1][cohort] = hatch_sim.at(cohort);	
			L_cohort_wolb[time1][cohort] = hatch_sim_wolb.at(cohort);

		}//end of if time1==hdate.at(cohort) loop	

//-------------------PUPAE------------------------------------------------------

		if (time1>hdate.at(cohort)+4 && emerg_flag.at(cohort)==0){ //if the cohort could have begun emergence and less than 1 individual have emerged so far
		//check that this ordering of the cohorts and times is okay
 			L_avg1 = 0;
			L_cum2 = 0;
			denom2 = 0;
			for (int etime = hdate.at(cohort)+5; etime <=time1; etime++){//all emergence times "etime" past to present
				L_cum=0;
				H_cum=0;
				count=cohort;
				for (int k=hdate.at(cohort); k<etime; k++) {
					L_cum += (L.at(k) + L_wolb.at(k));//cumulative exposure experienced at day before emergence time
				}
				L_cum-=H_cum;
				L_cum2 += (emerg_record[etime][cohort] + emerg_record_wolb[etime][cohort])*L_cum/(etime-hdate.at(cohort));//accumulating daily average exposure * no. emerged
				denom2 += emerg_record[etime][cohort] + emerg_record_wolb[etime][cohort];//cumulative no. emerged
			}
			if (denom2>0) L_avg1 = L_cum2/denom2;
			else L_avg1 = L_cum/(time1 - hdate.at(cohort));

			mn_dt_gam.at(cohort) = intc_TL + alpha_TL * pow(L_avg1, exp_TL);
			mn_dt_wolb_gam.at(cohort) = mn_dt_gam.at(cohort);

			if (mn_dt_gam.at(cohort)<0) mn_dt_gam.at(cohort)=0;
			if (mn_dt_wolb_gam.at(cohort)<0) mn_dt_wolb_gam.at(cohort)=0;

			if (denom2>1) {
				emerg_flag.at(cohort)=1;
				L_avg_cohort.at(cohort) = L_avg1;
			}

			std_dt_gam.at(cohort) = a_value + b_value * pow(L_avg1, b_value2);
			std_dt_wolb_gam.at(cohort) = a_value_wolb + b_value_wolb * pow(L_avg1, b_value2_wolb);

			if (L_avg1>L_max) std_dt_gam.at(cohort) = std_max;
			if (L_avg1>L_max) std_dt_wolb_gam.at(cohort) = std_max;
			if (std_dt_gam.at(cohort) <std_min) std_dt_gam.at(cohort) = std_min;
			if (std_dt_wolb_gam.at(cohort) <std_min) std_dt_wolb_gam.at(cohort) = std_min;
			if (L_avg1>L_max) {mn_dt_gam.at(cohort) = mean_max; mn_dt_wolb_gam.at(cohort) = mean_max;std_dt_gam.at(cohort) = std_max;std_dt_wolb_gam.at(cohort) = std_max;}

		}//end of if (time1>hdate.at(cohort)+4 && emerg_flag.at(cohort)==0)
		
		if (time1>hdate.at(cohort)+4){
			Dt_shp = pow(mn_dt_gam.at(cohort)/std_dt_gam.at(cohort),2);
			Dt_scl = mn_dt_gam.at(cohort)/Dt_shp;			
			Dt_shp_wolb = pow(mn_dt_wolb_gam.at(cohort)/std_dt_wolb_gam.at(cohort),2);
			Dt_scl_wolb = mn_dt_wolb_gam.at(cohort)/Dt_shp_wolb;			
			time_lag = time1 - hdate.at(cohort);
			//if (Dt_shp>0 && Dt_scl>0) 
			if (time_lag>5 && mn_dt_gam.at(cohort)>0 && Dt_shp>0 && Dt_scl>0 && Dt_shp<200 && Dt_scl<200){				 
				prob = gsl_cdf_gamma_P(time_lag-5,Dt_shp,Dt_scl) - gsl_cdf_gamma_P(time_lag-5-1,Dt_shp,Dt_scl);
			}
			else prob=0;
			if (time_lag>5 && mn_dt_wolb_gam.at(cohort)>0 && Dt_shp_wolb>0 && Dt_scl_wolb>0 && Dt_shp_wolb<200 && Dt_scl_wolb<200){
				prob_wolb = gsl_cdf_gamma_P(time_lag-5,Dt_shp_wolb,Dt_scl_wolb) - gsl_cdf_gamma_P(time_lag-5-1,Dt_shp_wolb,Dt_scl_wolb);
			}

			else prob_wolb=0;

			if (mn_dt_gam.at(cohort)==0){
				Dt_shp = 9.0;
				Dt_scl = 0.2;
				prob = gsl_cdf_gamma_P(time_lag-4,Dt_shp,Dt_scl) - gsl_cdf_gamma_P(time_lag-4-1,Dt_shp,Dt_scl);
			   
			}
			if (mn_dt_wolb_gam.at(cohort)==0){
			 	Dt_shp_wolb = 9.0;
			 	Dt_scl_wolb = 0.2;
				prob_wolb = gsl_cdf_gamma_P(time_lag-4,Dt_shp_wolb,Dt_scl_wolb) - gsl_cdf_gamma_P(time_lag-4-1,Dt_shp_wolb,Dt_scl_wolb);
			}

			no_emerge = L_cohort[time1-1][cohort] * prob/non_emerg_prob.at(cohort);
			no_emerge_wolb = L_cohort_wolb[time1-1][cohort] * prob_wolb/non_emerg_prob_wolb.at(cohort);

			if (no_emerge > L_cohort[time1-1][cohort]){
				no_emerge=L_cohort[time1-1][cohort];
			 	prob=1;
			}

			if (no_emerge<0){no_emerge=0; prob=0;}
			if (prob<0) prob=0;

			if (no_emerge_wolb > L_cohort_wolb[time1-1][cohort]){
				no_emerge_wolb=L_cohort_wolb[time1-1][cohort];
				prob_wolb=1;
			}

			if (no_emerge_wolb<0){no_emerge_wolb=0; prob_wolb=0;}
			if (prob_wolb<0) prob_wolb=0;

			non_emerg_prob.at(cohort) *= (1-prob/non_emerg_prob.at(cohort));
			non_emerg_prob_wolb.at(cohort) *= (1-prob_wolb/non_emerg_prob_wolb.at(cohort));
			L_cohort[time1][cohort] -= no_emerge;
			L_cohort_wolb[time1][cohort] -= no_emerge_wolb;
			emerg_record[time1][cohort]=no_emerge;
			emerg_record_wolb[time1][cohort]=no_emerge_wolb;
			mean_dt.at(cohort)+=(time1-hdate.at(cohort))*no_emerge;
			mean_dt_wolb.at(cohort)+=(time1-hdate.at(cohort))*no_emerge_wolb;
			no_emerg_tot.at(cohort)+=no_emerge;
			no_emerg_tot_wolb.at(cohort)+=no_emerge_wolb;
			P[time1][0] += no_emerge;
			P_wolb[time1][0] += no_emerge_wolb;
		}//end of if (time1>hdate.at(cohort)+4) loop

		L.at(time1)+=L_cohort[time1][cohort];
		L_wolb.at(time1)+=L_cohort_wolb[time1][cohort];


	}//end of cohort loop


//-------------------ADULTS ----------------------------------------------------

	P[time1][1] = P[time1-1][0]*survA;
	P_wolb[time1][1] = P_wolb[time1-1][0]*survA_wolb; 

	A[time1] = A[time1-1]*survA;
	A[time1] += P[time1-1][1]*survA;

	A_wolb[time1] = A_wolb[time1-1]*survA_wolb;
	A_wolb[time1] += P_wolb[time1-1][1]*survA_wolb;


	
//account for wolbachia additions 
	if (time1>release_day_init) A_wolb_imm[time1] = A_wolb_imm[time1-1]*survA_imm;

	if (time1==release_day && time1<release_end){
		A_wolb_imm[time1]+=release_size;
		release_day+=7;
	}

	
	if (time1>=tlagg) {
		if (A[time1-tlagg]>0) freqA = (A_wolb[time1-tlagg] + A_wolb_imm[time1-tlagg])/(A[time1-tlagg] + A_wolb[time1-tlagg] + A_wolb_imm[time1-tlagg]);
		else freqA=0;

		A_ovipos[time1] = 0.5*(1-sh*freqA)*A[time1-tlagg+tlagp] * pow(survA,tlagg-tlagp-1);
		A_ovipos_wolb[time1] = 0.5*A_wolb[time1-tlagg+tlagp] * pow(survA_wolb,tlagg-tlagp-1);
		A_ovipos_wolb_imm[time1] = 0.5*A_wolb_imm[time1-tlagg+tlagp] * pow(survA_imm, tlagg-tlagp-1);
	}


}//end of time loop


//work out cohort development time means and standard deviations
for (int i=0; i<no_cohorts; i++) {
	if (no_emerg_tot.at(i)>0) mean_dt.at(i)/=no_emerg_tot.at(i);
	else mean_dt.at(i)=0;
	if (no_emerg_tot_wolb.at(i)>0) mean_dt_wolb.at(i)/=no_emerg_tot_wolb.at(i);
	else mean_dt_wolb.at(i)=0;
	for (int j=0; j<maxtime; j++) {
		if (j>hdate.at(i)+4){
			std_dt.at(i) +=(j-hdate.at(i)-mean_dt.at(i))*(j-hdate.at(i)-mean_dt.at(i))*emerg_record[j][i];
			std_dt_wolb.at(i) += (j-hdate.at(i)-mean_dt_wolb.at(i))*(j-hdate.at(i)-mean_dt_wolb.at(i))*emerg_record_wolb[j][i];
		}
   	}
	if (no_emerg_tot.at(i)>0) std_dt.at(i) = sqrt(std_dt.at(i)/no_emerg_tot.at(i));
	else std_dt.at(i)=0;
	if (no_emerg_tot_wolb.at(i)>0) std_dt_wolb.at(i) = sqrt(std_dt_wolb.at(i)/no_emerg_tot_wolb.at(i));
	else std_dt_wolb.at(i)=0;
}



for (int i=0; i<maxtime; i++) {
	for (int j=0; j<no_cohorts; j++) {
		mu_p.at(i) += emerg_record[i][j];
		mu_p_wolb.at(i) += emerg_record_wolb[i][j];
	}
}


for (int i=0; i<no_cohorts; i++) lambda_out << lambda.at(i) << endl;

for (int i=0; i<maxtime; i++) A_wolb.at(i) += A_wolb_imm.at(i);

freq_end2 = A_wolb.at(release_end)/(A_wolb.at(release_end) + A.at(release_end));
freqA2_out << freq_end2 << endl;

freqL_end = L_wolb.at(release_end)/(L_wolb.at(release_end) + L.at(release_end));

//Work out the mean fecundity and development times

mean_dt1=0;mean_lambda=0;mean_std_dt=0;
for (int i=300; i<800; i++){
	mean_dt1+=mean_dt.at(i);
	mean_std_dt+=std_dt.at(i);
	mean_lambda += lambda.at(i);
} 
cout << endl;
cout << "mean_dt " << mean_dt1/500 <<" (the average (over time) of the mean development time of (uninfected) larvae) " << endl;
cout << "mean_lambda " << mean_lambda/500 << " (the average (over time) of per-capita female fecundity) " << endl;
cout << "mean_std_dt " << mean_std_dt/500 << " (the average (over time) of the standard deviation of the development time of (uninfected) larvae " <<  endl;

return 0;			
	

}




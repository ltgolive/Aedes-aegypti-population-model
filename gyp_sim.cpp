// The original gyp_sim.cpp was developed by Dr Penelope A. Hancock
// The development of gyp_sim.cpp is described in Hancock et al. 2016, "Predicting Wolbachia invasion dynamics in Aedes aegypti populations using models of density-dependent demographic traits", BMC Biology.
// This version (3) has been streamlined, annotated and excludes all Wolbachia infected mosquitoes and related parameters to only represent the life-cycle of a wild type population.

// Import libraries
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
#include <iomanip>

using namespace std;

// Type definitions (data structures)
typedef vector<int> Row_Int;
typedef vector<Row_Int> Matrix_Int;

typedef vector<double> Row_Double;
typedef vector<Row_Double> Matrix_Double;

// Declare functions - could some of these vectors be globals instead?
void read_data(ofstream& cohort_means_out, ofstream& cohort_stds_out, ofstream& mu_p_m_out,ofstream& L_m_out, ofstream& A_m_out, ofstream& mu_p_f_out,ofstream& L_f_out, ofstream& A_f_out, ofstream& lambda_out);
double Run_Model(ofstream& mu_p_m_out, ofstream& mu_p_f_out, ofstream& lambda_out);

// Globals
const int no_cohorts = 1000; // number of mosquito cohorts
const int maxtime = 2500; // simulation length (days; ~7 years)

// Declare data containers
vector<float> mean_dt(no_cohorts), std_dt(no_cohorts), A_ovipos(maxtime);
vector<double> mu_p_m(maxtime), L_m(maxtime), A_m(maxtime), mu_p_f(maxtime), L_f(maxtime), A_f(maxtime);
Matrix_Double P_m(maxtime, Row_Double(3)), P_f(maxtime, Row_Double(3));
vector<int> hdate(no_cohorts);

// Declare variables
double DI_cost, DI_L_cost;

//------------------------------------------- MAIN function -------------------------------------------

int main(){

// Declare input and output files
ofstream cohort_means_out, cohort_stds_out, mu_p_m_out, L_m_out, A_m_out, mu_p_f_out, L_f_out, A_f_out, lambda_out;

// Calls function: read_data
read_data(cohort_means_out, cohort_stds_out, mu_p_m_out, L_m_out, A_m_out, mu_p_f_out, L_f_out, A_f_out, lambda_out);

// Calls function: Run_Model - calls in output file streams that need to be written to during the simulation rather than after
Run_Model(mu_p_m_out, mu_p_f_out, lambda_out);

// Write cohort development times to file
for (int j=0; j<no_cohorts; j++) {
	cohort_means_out << mean_dt.at(j) << " ";
	cohort_stds_out << std_dt.at(j) << " ";
}
cohort_means_out << endl;
cohort_stds_out << endl;

// Write daily pupal emergence to file
for (int j=0; j<maxtime; j++) {
	mu_p_m_out << mu_p_m.at(j) << " ";
	mu_p_f_out << mu_p_f.at(j) << " ";
}
mu_p_m_out << endl;
mu_p_f_out << endl;

// Write daily larval and adult counts to file
for (int j=0; j<maxtime; j++) {
	L_m_out << L_m.at(j) << " ";
	L_f_out << L_f.at(j) << " ";
	A_m_out << A_m.at(j) << " ";
	A_f_out << A_f.at(j) << " ";
}
L_m_out << endl; 
L_f_out << endl;
A_m_out << endl; 
A_f_out << endl;

return 0;

}

//------------------------------------------- READ DATA function -------------------------------------------

// Call in input and output files
void read_data(ofstream& cohort_means_out, ofstream& cohort_stds_out, ofstream& mu_p_m_out, ofstream& L_m_out, ofstream& A_m_out, ofstream& mu_p_f_out, ofstream& L_f_out, ofstream& A_f_out, ofstream& lambda_out){

// Read in hatch dates
ifstream hdate_data, hatch_data;
hdate_data.open("hdate3.txt");
for (int i=0; i<no_cohorts; i++){
hdate_data >> hdate.at(i);
}

// Open input and output files in inits file
string cohort_means_file, cohort_stds_file, mu_p_m_file, L_m_file, A_m_file, mu_p_f_file, L_f_file, A_f_file, lambda_file;
cin >> cohort_means_file >> cohort_stds_file >> mu_p_m_file >> L_m_file >> A_m_file >> mu_p_f_file >> L_f_file >> A_f_file >> lambda_file;

// Print file names to console
cout << " cohort_means_file " << cohort_means_file << " cohort_stds_file " << cohort_stds_file << " mu_p_m_file " << mu_p_m_file << " L_m_file " << L_m_file << " A_m_file " <<  A_m_file  << " mu_p_f_file " << mu_p_f_file << " L_f_file " << L_f_file << " A_f_file " <<  A_f_file  << " lambda_file " << lambda_file << endl;

// Output files
cohort_means_out.open(cohort_means_file); // mean development times per cohort
cohort_stds_out.open(cohort_stds_file); // standard deviations of development times per cohort
mu_p_m_out.open(mu_p_m_file); // number of pupae eclosing per day (males)
L_m_out.open(L_m_file); // total number of larvae per day (males)
A_m_out.open(A_m_file); // total number of adults per day (males)
mu_p_f_out.open(mu_p_f_file); 
L_f_out.open(L_f_file); 
A_f_out.open(A_f_file); 
lambda_out.open(lambda_file); // per-capita female fecundity per cohort

// Read in key parameter values in inits file
cin >> DI_cost;
cin >> DI_L_cost;
cout << "DI_cost " << DI_cost << " (additional density-INdependent daily mortality experienced by adults in the field environment) " << endl;
cout << "DI_L_cost " << DI_L_cost << " (additional density-INdependent daily mortality experienced by larvae in the field environment) " << endl;

}

//------------------------------------------- RUN MODEL function -------------------------------------------

// Call in input and output files
double Run_Model(ofstream& mu_p_m_out, ofstream& mu_p_f_out, ofstream& lambda_out){

// Create data containers (by cohort or day)
vector<double> non_emerg_prob(no_cohorts), L_avg_cohort(no_cohorts), L_avg(no_cohorts), mn_dt_gam(no_cohorts), std_dt_gam(no_cohorts), no_emerg_tot(no_cohorts), hatch_sim(no_cohorts), lambda(no_cohorts);
vector<int> emerg_flag(no_cohorts);

Matrix_Double L_cohort_m(maxtime, Row_Double(no_cohorts)); // number of larvae per day per cohort (2500 days x 1000 cohorts) (males)
Matrix_Double L_cohort_f(maxtime, Row_Double(no_cohorts)); 
Matrix_Double emerg_record_m(maxtime, Row_Double(no_cohorts)); // number of larvae emerging as pupae per day per cohort (2500 days x 1000 cohorts) (males)
Matrix_Double emerg_record_f(maxtime, Row_Double(no_cohorts)); 

// Declare temporary variables
double surv_L, L_avg1, L_cum2, L_cum, denom2, Dt_shp, Dt_scl, prob, no_emerge, no_emerge_m, no_emerge_f, L_avg_f, lambda_max, survA_imm, intc_TL, alpha_TL, exp_TL, a_value, b_value, b_value2, intc_f, alpha_f, survA, lambda1, lambda_min, mean_max, std_max, std_min, mean_dt1, mean_lambda, mean_std_dt, sex_ratio, emerg_record, L_cohort;
int first_hatch_date, max_cohort, L_max, time_lag, tlagh, tlagg, tlagi, tlagp, min_cohort, max_dt_int, dens_lag;

// Time lags (locals)
first_hatch_date=hdate.at(0);
tlagh = 5; // oviposition to egg hatching
tlagg = 6; // emerging as adults to being ready to lay eggs
tlagi = 21; // time over which larval density is averaged (to calculate effect of density)
tlagp = 2; // pupal development time

// Biological constants (locals)
intc_f = 28.0; // fecundity log-linear relationship - intercept (fitted from experimental data)
alpha_f = -3.3; // fecundity log-linear relationship - slope (fitted from experimental data)

intc_TL = 1.8; // mean larval development time gamma distribution - intercept (fitted from experimental data)
alpha_TL = 0.536; // mean larval development time gamma distribution - slope (fitted from experimental data)
exp_TL = 0.533; // mean larval development time gamma distribution - exponent (fitted from experimental data)

a_value = 0.22; // std dev larval development time gamma distribution - intercept (fitted from experimental data)
b_value = 0.0168; // std dev larval development time gamma distribution - slope (fitted from experimental data)
b_value2 = 0.867; // std dev larval development time gamma distribution - exponent (fitted from experimental data)

sex_ratio = 0.5; // sex ratio (offspring)

// Biological limits (locals)
max_dt_int = 100; // cohorts older than 100 days are ignored
survA = 1-0.03-DI_cost; // daily survival rate of adults - base mortality is 3% per day, plus any additional cost (DI_cost)
lambda_max = 14; // a female can lay at most 14 eggs per day
lambda_min = 0.5; // a female lays at least 0.5 eggs per day even under very crowded conditions
mean_max = 60; // if crowding is extreme, mean development time is capped at 60 days (no minimum?)
std_max = 40; // if crowding is extreme, std dev development time is capped at ±60 days from mean
std_min = 1.0; // if crowding is minimal, std dev development time is capped at ±1 day from mean
L_max = 4200; // the larval density above which the caps on development time kick in

survA_imm=survA; // survival rate for adults (same as above)

// Initialising cohorts
// Arrays (by day and cohort)
for (int i=0; i<maxtime; i++){
	for (int j=0; j<no_cohorts; j++){
		L_cohort_m[i][j] = 0; // number of larvae (males)
		L_cohort_f[i][j] = 0; 
		emerg_record_m[i][j] = 0; // number of new pupae (males)
		emerg_record_f[i][j] = 0; 
	}
}

// Arrays (by cohort)  
for (int i=0; i<no_cohorts; i++){
	non_emerg_prob.at(i) = 1; // cumulative probability that a larva has not yet emerged as a pupa (1 = 100%; decreases over time)
	emerg_flag.at(i) = 0; // flag to track if enough pupae have emerged to start calculating statistics
	L_avg_cohort.at(i) = 0; // running average larval density experienced (estimated)
	L_avg.at(i) = 0; // average larval density experienced (final)
	mn_dt_gam.at(i) = 999; // mean larval development time (used for gamma distribution)
	std_dt_gam.at(i) = 1; // standard deviation of larval development time (used for gamma distribution)
	mean_dt.at(i) = 0; // mean larval development time (final)
	std_dt.at(i) = 0; // standard deviation of larval development time (final)
	no_emerg_tot.at(i) = 0; // number of new larvae
	hatch_sim.at(i) = 0; // number of larvae hatching into each cohort (i.e. size of cohort at time of hatching)
}

// Arrays (by day) 
for (int i=0; i<maxtime; i++) {
	A_ovipos.at(i) = 0; // number of reproductively active females
	L_m.at(i) = 0; // number of larvae (males)
	L_f.at(i) = 0; 
	mu_p_m.at(i) = 0; // number of new adults (males)
	mu_p_f.at(i) = 0; 
	// With 3 columns, one for each pupal stage (0-2 days)
	for (int j=0; j<3; j++){ 
		P_m[i][j]=0; // number of pupae (males)
		P_f[i][j]=0; 
	}
}

// Arrays (by day - starts at 1 then 0 for all other days)
A_m.at(0) = 1;
A_f.at(0) = 1;
for (int i=1; i<maxtime; i++) {
	A_m.at(i)=0; // number of adults (males)
	A_f.at(i)=0; 
}


//------------------------------------- RUN MODEL function: Time Loop -------------------------------------

// The whole code repeats itself every day from day 1 to 2500
for (int time1=1; time1<maxtime; time1++){

	// Daily survival rate of larvae
	surv_L = 0.95 - DI_L_cost;

	// Track how many active cohorts the code needs to loop over
	max_cohort=-1; // find most recent cohort that has already hatched 
	for (int i=0; i<no_cohorts; i++) if (hdate.at(i)<=time1) max_cohort = i;
	min_cohort=no_cohorts-1; // find oldest cohort still worth tracking (ignore if older than 100 days [max_dt_int])
	for (int i=no_cohorts-1; i>=0; i--) if (hdate.at(i)>=time1 - max_dt_int) min_cohort = i;	

	//----------------------------------- RUN MODEL function: Cohort Loop ----------------------------------

	// The whole code repeats itself for every active cohort
	for (int cohort=min_cohort; cohort<=max_cohort; cohort++){

		//------------------------------------------------- Larvae -----------------------------------------

		// Calculate number of larvae alive today
		L_cohort_m[time1][cohort] = L_cohort_m[time1-1][cohort] * surv_L; // larvae alive yesterday * survival rate (males)
		L_cohort_f[time1][cohort] = L_cohort_f[time1-1][cohort] * surv_L; 
		
		// Check if today is a hatch date - if yes, calculate number of new larvae hatching today
		if (time1 == hdate.at(cohort)){

			// Calculate starting point of look-back window 
			dens_lag = time1-tlagh-tlagg-tlagi; // 32 days ago, when parent larvae started developing

			// Calculate average larval density experienced by this cohort (parents)
			L_avg.at(cohort) = 0; // initialise

			// Case 1: full window within simulation (i.e. start ≥ first hatch date)
			if (dens_lag>=hdate.at(1)){
				// Total larvae (21 days)/21 days = larval density experienced
				for (int j=time1-tlagh-tlagg-tlagi; j<time1-tlagh-tlagg; j++) L_avg.at(cohort) += L_m.at(j) + L_f.at(j);
					L_avg.at(cohort) = L_avg.at(cohort)/tlagi; 
			}

			// Case 2: window partially overlaps with simulation (i.e. start < first hatch date, but still after day 0)
			else if (dens_lag>=0) {
				// Case 2a: end < first hatch date (i.e. full window outside simulation) - I don't think this can happen, first hatch date is day 4
				if (time1-tlagh-tlagg<=hdate.at(1)) L_avg.at(cohort) = 60; // default average density (60)
				// Case 2b: end > first hatch date (i.e. window partially within simulation)
				if (time1-tlagh-tlagg>hdate.at(1)) {
					for (int j=dens_lag; j<=hdate.at(1); j++) L_avg.at(cohort) += 60; // section before first hatch date = default (60 larvae per day)
					for (int j=hdate.at(1)+1; j<=time1-tlagh-tlagg; j++) L_avg.at(cohort) += L_m.at(j) + L_f.at(j); // section after first hatch date = total larvae per day (real data)
					L_avg.at(cohort) = L_avg.at(cohort)/tlagi;
				} 
			}

			// Case 3: window starts before day 0
			else L_avg.at(cohort) = 60; // default average density (60)

			// Calculate per-capita female fecundity of mothers (lambda1 -average number of eggs laid per female per day)
			L_avg_f = L_avg.at(cohort); // average larval density experienced by mothers (as calculated above)

			lambda1 = intc_f + alpha_f*log(L_avg_f); // calculate fecundity using log-linear function (higher density experienced = lower fecundity)

			// Limit fecundity to biologically realistic limits
			if (lambda1<lambda_min) lambda1 = lambda_min; // female cannot lay less than 0.5 eggs per day
			if (lambda1>lambda_max) lambda1 = lambda_max; // female cannot lay more than 14 eggs per day

			// Store fecundity for this cohort
			lambda.at(cohort) = lambda1;

			// Calculate number of new larvae hatching today
 			if (time1>tlagh+tlagg){ // prevents hatching from happening until enough time has passed for first adults in simulation to lay eggs (11 days)
				hatch_sim.at(cohort) = lambda1*A_ovipos.at(time1-tlagh); // fecundity * number of adult females ready to lay eggs (at c - Th)
			}

			// Store newly hatched larvae as starting larval population for this new cohort
			L_cohort_m[time1][cohort] = sex_ratio*hatch_sim.at(cohort);	
			L_cohort_f[time1][cohort] = (1-sex_ratio)*hatch_sim.at(cohort);	

		} // End loop for newly hatching larvae	


		//------------------------------------------------- Pupae ------------------------------------------

		// Calculate mean and std dev of larval development time for each cohort (expected)
		if (time1>hdate.at(cohort)+4 && emerg_flag.at(cohort)==0){ // only run if ≥5 days since cohort hatched (min dev time) and <2 pupae have emerged
		
			// Initialise
			L_cum2 = 0; // weighted average
			denom2 = 0; // number of larvae emerged as pupae
			L_avg1 = 0; // average larval density

			// Calculate weighted average larval density experienced by each cohort
			for (int etime = hdate.at(cohort)+5; etime <=time1; etime++){ // for each emergence day since 5 days since cohort hatched

				// Initialise
				L_cum = 0; // cumulative larval density

				// Calculate cumulative number of larvae from hatch to emergence day
				for (int k=hdate.at(cohort); k<etime; k++) {
					L_cum += L_m.at(k) + L_f.at(k);
				}

				// Calculate total larvae emerging (sum by sex)
				emerg_record = emerg_record_m[etime][cohort] + emerg_record_f[etime][cohort];

				// Calculate weighted average
				L_cum2 += emerg_record * L_cum/(etime-hdate.at(cohort)); // number emerged * daily average
				
				// Store cumulative number of larvae emerged as pupae
				denom2 += emerg_record;
			}

			// Check when to use weighted average
			if (denom2>0) L_avg1 = L_cum2/denom2; // yes - larvae have started emerging
			else L_avg1 = L_cum/(time1 - hdate.at(cohort)); // no - no larvae have emerged yet (use unweighted)

			// Calculate expected mean of larval development time
			mn_dt_gam.at(cohort) = intc_TL + alpha_TL * pow(L_avg1, exp_TL); // calculate mean using power-law formula (higher density experienced = longer development time)
			if (mn_dt_gam.at(cohort)<0) mn_dt_gam.at(cohort)=0; // floor at 0 to avoid nonsensical values
			
			// Calculate expected std dev of larval development time
			std_dt_gam.at(cohort) = a_value + b_value * pow(L_avg1, b_value2); // calculate std dev using power-law formula (higher density experienced = more variable development time)
			
			// Limit mean and std dev to biologically realistic limits
			if (L_avg1>L_max) {mn_dt_gam.at(cohort) = mean_max; std_dt_gam.at(cohort) = std_max;} // if density exceeded (L_max = 4200) then mean and std dev set to max allowed values
			if (std_dt_gam.at(cohort) < std_min) std_dt_gam.at(cohort) = std_min; // if std dev is less than minimum (1) then set to minimum

			// Flag - if ≥2 pupae have emerged, stop updating mean and std dev
			if (denom2>1) {
				emerg_flag.at(cohort)=1;
				L_avg_cohort.at(cohort) = L_avg1;
			}

		} // End loop for calculating mean and std dev of larval development time (expected)
		
		// Calculate number of larvae emerging as pupae
		if (time1>hdate.at(cohort)+4){

			// Convert mean and std dev to scale and shape parameters
			Dt_shp = pow(mn_dt_gam.at(cohort)/std_dt_gam.at(cohort),2);
			Dt_scl = mn_dt_gam.at(cohort)/Dt_shp;						
			time_lag = time1 - hdate.at(cohort);

			// Calculate probability of emerging today
			if (time_lag>5 && mn_dt_gam.at(cohort)>0 && Dt_shp>0 && Dt_scl>0 && Dt_shp<200 && Dt_scl<200){				 
				prob = gsl_cdf_gamma_P(time_lag-5,Dt_shp,Dt_scl) - gsl_cdf_gamma_P(time_lag-5-1,Dt_shp,Dt_scl);
			}
			else prob=0;

			// Safety checks
			if (mn_dt_gam.at(cohort)==0){
				Dt_shp = 9.0;
				Dt_scl = 0.2;
				prob = gsl_cdf_gamma_P(time_lag-4,Dt_shp,Dt_scl) - gsl_cdf_gamma_P(time_lag-4-1,Dt_shp,Dt_scl);
			}

			// Calculate total number of larvae (sum by sex)
			L_cohort = L_cohort_m[time1-1][cohort] + L_cohort_f[time1-1][cohort];
			
			// Convert probability to numbers
			no_emerge = L_cohort * prob/non_emerg_prob.at(cohort); // total
			no_emerge_m = sex_ratio * no_emerge; // males
            no_emerge_f = (1-sex_ratio) * no_emerge; // females
			
			// Sanity checks
			if (no_emerge > L_cohort){
				no_emerge=L_cohort;
			 	prob=1;
			}

			if (no_emerge<0){no_emerge=0; prob=0;}
			if (prob<0) prob=0;

			// Updates based on calculations
			non_emerg_prob.at(cohort) *= (1-prob/non_emerg_prob.at(cohort)); // update probability of not emerging
	
			L_cohort_m[time1][cohort] -= no_emerge_m; // remove emerging larvae from larvae compartment (males)
            L_cohort_f[time1][cohort] -= no_emerge_f; 
			
			emerg_record_m[time1][cohort] = no_emerge_m; // store emerging larvae (males)
			emerg_record_f[time1][cohort] = no_emerge_f; 

			mean_dt.at(cohort)+=(time1-hdate.at(cohort))*no_emerge; // accumulate mean development time
			no_emerg_tot.at(cohort)+=no_emerge;

			P_m[time1][0] += no_emerge_m; // add emerging larvae to pupae compartment (males)
			P_f[time1][0] += no_emerge_f; 

		} // End loop for calculating number of larvae emerging as pupae

		// Store total number of larvae
		L_m.at(time1) += L_cohort_m[time1][cohort]; // add up total larvae from each cohort
		L_f.at(time1) += L_cohort_f[time1][cohort];

	} // End loop for calculations for each active cohort


	//--------------------------------------------- Adults --------------------------------------------

	// Calculate number of pupae alive today
	P_m[time1][1] = P_m[time1-1][0] * survA; // move into second pupal day
	P_f[time1][1] = P_f[time1-1][0] * survA;

	// Calculate number of adults alive today
	A_m[time1] = A_m[time1-1] * survA;
	A_f[time1] = A_f[time1-1] * survA;

	// Calculate number of pupae graduating into adults today
	A_m[time1] += P_m[time1-1][1] * survA;
	A_f[time1] += P_f[time1-1][1] * survA;

	// Calculate number of adult females ready to lay eggs
	if (time1>=tlagg) { // adults must be at least 6 days old to be ready to lay eggs
		A_ovipos[time1] = A_f[time1-tlagg+tlagp] * pow(survA,tlagg-tlagp-1); // females that became adults 4 days (eclosed 6 days ago) that are still alive
	}

} // End time loop

//--------------------------------------------- Model Outputs ------------------------------------------

// Calculate observed mean and standard deviation of larval development time for each cohort
for (int i=0; i<no_cohorts; i++) {
	if (no_emerg_tot.at(i)>0) mean_dt.at(i)/=no_emerg_tot.at(i); // mean
	else mean_dt.at(i)=0;
	for (int j=0; j<maxtime; j++) { // accumulate squared deviations
		if (j>hdate.at(i)+4){
			emerg_record = emerg_record_f[j][i] + emerg_record_m[j][i];
			std_dt.at(i) +=(j-hdate.at(i)-mean_dt.at(i))*(j-hdate.at(i)-mean_dt.at(i))*emerg_record;
		}
   	}
	if (no_emerg_tot.at(i)>0) std_dt.at(i) = sqrt(std_dt.at(i)/no_emerg_tot.at(i)); // std dev
	else std_dt.at(i)=0;
}

// Calculate total number of pupae that emerged each day
for (int i=0; i<maxtime; i++) {
	for (int j=0; j<no_cohorts; j++) {
		mu_p_m.at(i) += emerg_record_m[i][j];
        mu_p_f.at(i) += emerg_record_f[i][j];
	}
}

// Write per-capita female fecundity for each cohort to the output file
for (int i=0; i<no_cohorts; i++) lambda_out << lambda.at(i) << endl;

// Calculate summary statistics
mean_dt1=0;mean_lambda=0;mean_std_dt=0;
for (int i=300; i<800; i++){ // averaged across cohorts 300 to 800
	mean_dt1+=mean_dt.at(i); // mean
	mean_std_dt+=std_dt.at(i); // fecundity
	mean_lambda += lambda.at(i); // std dev
} 

// Print summary statistics
cout << endl;
cout << "mean_dt " << mean_dt1/500 <<" (the average (over time) of the mean development time of (wildtype) larvae) " << endl;
cout << "mean_lambda " << mean_lambda/500 << " (the average (over time) of per-capita female fecundity) " << endl;
cout << "mean_std_dt " << mean_std_dt/500 << " (the average (over time) of the standard deviation of the development time of (wildtype) larvae " <<  endl;

return 0;			

}




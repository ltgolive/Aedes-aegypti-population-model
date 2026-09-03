// The original gyp_sim.cpp was developed by Dr Penelope A. Hancock
// The development is described in Hancock et al. 2016, "Predicting Wolbachia invasion dynamics in Aedes aegypti populations using models of density-dependent demographic traits", BMC Biology.

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

// Declare functions
void read_data(ofstream& cohort_means_out, ofstream& cohort_stds_out, ofstream& lambda_out, ofstream& Lm_out, ofstream& Lf_out, ofstream& Pm_emerg_out, ofstream& Pf_emerg_out, ofstream& Am_out, ofstream& Af_out, ofstream& Lm_gm_out, ofstream& Lf_gm_out, ofstream& Pm_gm_emerg_out, ofstream& Am_gm_out, ofstream& Am_gm_release_out, ofstream& Am_freq_out, ofstream& Am_gm_freq_out, ofstream& Am_gm_release_freq_out);
double Run_Model(ofstream& lambda_out);

// Declare global variables
const int no_cohorts = 1000; // number of mosquito cohorts
const int maxtime = 2500; // simulation length (days; ~7 years)
int release_size, release_day, release_period;
double survA_wt, survA_rel_perc, survA_gm_perc;
double survL_wt, survL_gm_perc;
double lambda_intc, lambda_alpha;
double dt_mean_intc, dt_mean_alpha, dt_mean_exp;
double dt_std_intc, dt_std_alpha, dt_std_exp;

// Declare global data containers
vector<float> dt_mean(no_cohorts), dt_std(no_cohorts);
vector<double> Lm(maxtime), Pm_emerg(maxtime), Am(maxtime), Lf(maxtime), Pf_emerg(maxtime), Af(maxtime); // wildtype
Matrix_Double Pm(maxtime, Row_Double(2)), Pf(maxtime, Row_Double(2)); // wildtype
vector<double> Lm_gm(maxtime), Pm_gm_emerg(maxtime), Am_gm(maxtime), Am_gm_release(maxtime), Lf_gm(maxtime); // GM
Matrix_Double Pm_gm(maxtime, Row_Double(2)); // GM
vector<double> Am_freq(maxtime), Am_gm_freq(maxtime), Am_gm_release_freq(maxtime);
vector<double> Af_mated_wt(maxtime), Af_mated_gm(maxtime), Af_mated_rel(maxtime);
vector<int> hdate(no_cohorts);


//------------------------------------------- MAIN function -------------------------------------------

int main(){

// Declare input and output files
ofstream cohort_means_out, cohort_stds_out, lambda_out;
ofstream Lm_out, Lf_out, Pm_emerg_out, Pf_emerg_out, Am_out, Af_out;
ofstream Lm_gm_out, Lf_gm_out, Pm_gm_emerg_out, Am_gm_out, Am_gm_release_out;
ofstream Am_freq_out, Am_gm_freq_out, Am_gm_release_freq_out;

// Calls function: read_data
read_data(cohort_means_out, cohort_stds_out, lambda_out, Lm_out, Lf_out, Pm_emerg_out, Pf_emerg_out, Am_out, Af_out, Lm_gm_out, Lf_gm_out, Pm_gm_emerg_out, Am_gm_out, Am_gm_release_out, Am_freq_out, Am_gm_freq_out, Am_gm_release_freq_out);

// Calls function: Run_Model - calls in output file streams that need to be written to during the simulation rather than after
Run_Model(lambda_out);

// Write cohort development times to file
for (int j=0; j<no_cohorts; j++) {
	cohort_means_out << dt_mean.at(j) << " ";
	cohort_stds_out << dt_std.at(j) << " ";
}
cohort_means_out << endl;
cohort_stds_out << endl;

// Write daily larval count, pupal emergence and adult count to file
for (int j=0; j<maxtime; j++) {
	Lm_out << Lm.at(j) << " ";
	Lf_out << Lf.at(j) << " ";
	Lm_gm_out << Lm_gm.at(j) << " ";
    Lf_gm_out << Lf_gm.at(j) << " ";

	Pm_emerg_out << Pm_emerg.at(j) << " ";
	Pf_emerg_out << Pf_emerg.at(j) << " ";
	Pm_gm_emerg_out << Pm_gm_emerg.at(j) << " ";

	Am_out << Am.at(j) << " ";
	Af_out << Af.at(j) << " ";
    Am_gm_out << Am_gm.at(j) << " ";
    Am_gm_release_out << Am_gm_release.at(j) << " ";

	Am_freq_out << Am_freq.at(j) << " ";
	Am_gm_freq_out << Am_gm_freq.at(j) << " ";
	Am_gm_release_freq_out << Am_gm_release_freq.at(j) << " ";
}
Lm_out << endl; 
Lf_out << endl;
Lm_gm_out << endl;
Lf_gm_out << endl;

Pm_emerg_out << endl;
Pf_emerg_out << endl;
Pm_gm_emerg_out << endl;

Am_out << endl; 
Af_out << endl;
Am_gm_out << endl;
Am_gm_release_out << endl;

Am_freq_out << endl;
Am_gm_freq_out << endl;
Am_gm_release_freq_out << endl;

return 0;

}

//------------------------------------------- READ DATA function -------------------------------------------

// Call in input and output files
void read_data(ofstream& cohort_means_out, ofstream& cohort_stds_out, ofstream& lambda_out, ofstream& Lm_out, ofstream& Lf_out, ofstream& Pm_emerg_out, ofstream& Pf_emerg_out, ofstream& Am_out, ofstream& Af_out, ofstream& Lm_gm_out, ofstream& Lf_gm_out, ofstream& Pm_gm_emerg_out, ofstream& Am_gm_out, ofstream& Am_gm_release_out, ofstream& Am_freq_out, ofstream& Am_gm_freq_out, ofstream& Am_gm_release_freq_out){

// Read in hatch dates
ifstream hdate_data;
hdate_data.open("hdate3.txt");
for (int i=0; i<no_cohorts; i++){
hdate_data >> hdate.at(i);
}

// Open input and output files in inits file
string cohort_means_file, cohort_stds_file, lambda_file;
string Lm_file, Lf_file, Pm_emerg_file, Pf_emerg_file, Am_file, Af_file;
string Lm_gm_file, Lf_gm_file, Pm_gm_emerg_file, Am_gm_file, Am_gm_release_file;
string Am_freq_file, Am_gm_freq_file, Am_gm_release_freq_file;

cin >> cohort_means_file >> cohort_stds_file >> lambda_file >> Lm_file >> Lf_file >> Pm_emerg_file >> Pf_emerg_file >> Am_file >> Af_file >> Lm_gm_file >> Lf_gm_file >> Pm_gm_emerg_file >> Am_gm_file >> Am_gm_release_file >> Am_freq_file >> Am_gm_freq_file >> Am_gm_release_freq_file;

// Print file names to console
cout << " cohort_means_file " << cohort_means_file << " cohort_stds_file " << cohort_stds_file << " lambda_file " << lambda_file << endl;
cout << " Lm_file " << Lm_file << " Lf_file " << Lf_file << " Lm_gm_file " << Lm_gm_file << " Lf_gm_file " << Lf_gm_file << endl;
cout << " Pm_emerg_file " << Pm_emerg_file << " Pf_emerg_file " << Pf_emerg_file << " Pm_gm_emerg_file " << Pm_gm_emerg_file << endl;
cout << " Am_file " << Am_file << " Af_file " << Af_file << " Am_gm_file " << Am_gm_file << " Am_gm_release_file " << Am_gm_release_file << endl;
cout << " Am_freq_file " << Am_freq_file << " Am_gm_freq_file " << Am_gm_freq_file << " Am_gm_release_freq_file " << Am_gm_release_freq_file << endl;

// Output files
cohort_means_out.open(cohort_means_file); // mean development times per cohort
cohort_stds_out.open(cohort_stds_file); // standard deviations of development times per cohort
lambda_out.open(lambda_file); // per-capita female fecundity per cohort

Lm_out.open(Lm_file); // total number of larvae per day     
Lf_out.open(Lf_file);   

Pm_emerg_out.open(Pm_emerg_file); // total number of pupae emerging per day       
Pf_emerg_out.open(Pf_emerg_file); 

Am_out.open(Am_file); // total number of adults per day                  
Af_out.open(Af_file);   

Lm_gm_out.open(Lm_gm_file);              
Lf_gm_out.open(Lf_gm_file);  

Pm_gm_emerg_out.open(Pm_gm_emerg_file); 
  
Am_gm_out.open(Am_gm_file);                
Am_gm_release_out.open(Am_gm_release_file);

Am_freq_out.open(Am_freq_file);
Am_gm_freq_out.open(Am_gm_freq_file);
Am_gm_release_freq_out.open(Am_gm_release_freq_file);

// Read in key parameter values in inits file and print to check correctly read
cin >> survA_wt >> survA_gm_perc >> survA_rel_perc;
cin >> survL_wt >> survL_gm_perc;
cin >> lambda_intc >> lambda_alpha;
cin >> dt_mean_intc >> dt_mean_alpha >> dt_mean_exp;
cin >> dt_std_intc >> dt_std_alpha >> dt_std_exp;
cin >> release_size >> release_day >> release_period;

cout << "Adult daily survival: WT: " << survA_wt << " 2nd-gen: " << survA_gm_perc << " Released: " << survA_rel_perc << endl;
cout << "Larval daily survival: WT: " << survL_wt << " 2nd-gen: " << survL_gm_perc << endl;
cout << "Density-dependent parameters:" << endl;
cout << "Fecundity: " << "lambda_intc: " << lambda_intc << " lambda_alpha: " << lambda_alpha << endl;
cout << "Development time (mean): " << "dt_mean_intc: " << dt_mean_intc << " dt_mean_alpha: " << dt_mean_alpha << " dt_mean_exp: " << dt_mean_exp << endl;
cout << "Development time (std): " << "dt_std_intc: " << dt_std_intc << " dt_std_alpha: " << dt_std_alpha << " dt_std_exp: " << dt_std_exp << endl;
cout << "Release size: " << release_size << endl;
cout << "Release day: " << release_day << endl;
cout << "Release period: " << release_period << " weeks" <<endl;
}

//------------------------------------------- RUN MODEL function -------------------------------------------

// Call in input and output files
double Run_Model(ofstream& lambda_out){

// Time lags
int tlagh = 5; // oviposition to egg hatching
int tlagl = 21; // time over which larval density is averaged
int tlagp = 2; // pupal development time
int tlagg = 4; // mating to oviposition
int tlagm = 2; // emerging as adults and being ready to mate (males)
int tlagf = 1; // emerging as adults and being ready to mate (females)

// Mortality
double survA_rel = survA_wt * survA_rel_perc; 
double survA_gm = survA_wt * survA_gm_perc; 
double survL_gm = survL_wt * survL_gm_perc;

// Fecundity caps (log-linear relationship)
double lambda_max = 14; 
double lambda_min = 0.5;

// Pupal emergence caps (gamma distribution)
double dt_mean_max = 60; // if crowding is extreme, mean development time is capped at 60 days
double dt_std_max = 40; // if crowding is extreme, std dev development time is capped at ±60 days from mean
double dt_std_min = 1.0; // if crowding is minimal, std dev development time is capped at ±1 day from mean
int L_max = 4200; // the larval density above which the above caps on development time kick in

// Mating competitiveness
double comp_wt = 1.0;
double comp_gm = 1.0;
double comp_rel = 1.0;

// Other
double sex_ratio = 0.5; 
int cohort_age_max = 100; // cohorts older than 100 days are ignored

// Temporary variables
int max_cohort, min_cohort; // active cohorts

double L_avg_f, lambda1; // hatching larvae/fecundity
int start, end; 

double L_cum, L_cum2, denom, L_avg_gam1, P_emerg_cohort; // expected mean & std dev
int time_lag;

double dt_shp, dt_scl, prob, Pm_emerge_today, Pf_emerge_today, Pm_gm_emerge_today, Pf_gm_emerge_today; // pupal emergence


// Create local data containers
// By day and cohort (2500 days x 1000 cohorts)
// Wildtype
Matrix_Double Lm_cohort(maxtime, Row_Double(no_cohorts, 0.0)); // total number of larvae (male)
Matrix_Double Lf_cohort(maxtime, Row_Double(no_cohorts, 0.0)); 
Matrix_Double Pm_emerg_cohort(maxtime, Row_Double(no_cohorts, 0.0)); // number of new pupae emerging (male)
Matrix_Double Pf_emerg_cohort(maxtime, Row_Double(no_cohorts, 0.0)); 
// GM
Matrix_Double Lm_gm_cohort(maxtime, Row_Double(no_cohorts, 0.0)); 
Matrix_Double Lf_gm_cohort(maxtime, Row_Double(no_cohorts, 0.0)); 
Matrix_Double Pm_gm_emerg_cohort(maxtime, Row_Double(no_cohorts, 0.0)); 
Matrix_Double Pf_gm_emerg_cohort(maxtime, Row_Double(no_cohorts, 0.0));

// By cohort 
vector<double> non_emerg_prob(no_cohorts, 1.0); // cumulative probability that a larva has not yet emerged as a pupa (1 = 100%; decreases over time)
vector<int> emerg_flag(no_cohorts, 0); // flag to track if enough pupae have emerged to start calculating statistics
vector<double> L_avg_gam(no_cohorts, 0.0); // running average larval density experienced (estimated)
vector<double> L_avg(no_cohorts, 0.0); // average larval density experienced (final)
vector<double> dt_mean_gam(no_cohorts, 999.0); // mean larval development time (gamma distribution)
vector<double> dt_std_gam(no_cohorts, 1.0); // std dev of larval development time (gamma distribution)
vector<double> P_emerg(no_cohorts, 0.0); // total number of pupae that emerged
vector<double> P_gm_emerg(no_cohorts, 0.0);
vector<double> lambda(no_cohorts, 0.0); // per-capita female fecundity


// Initialise global data containers
// By cohort
for (int i=0; i<no_cohorts; i++){
	dt_mean.at(i) = 0.0; // mean larval development time (final)
	dt_std.at(i) = 0.0; // std dev of larval development time (final)
}

// By day
for (int i=0; i<maxtime; i++) {
	Lm.at(i) = 0.0; // number of larvae
	Lf.at(i) = 0.0;
	Lm_gm.at(i) = 0.0; 
	Lf_gm.at(i) = 0.0; 

	Pm_emerg.at(i) = 0.0; // number of new pupae
	Pf_emerg.at(i) = 0.0; 
	Pm_gm_emerg.at(i) = 0.0;

	Af_mated_wt.at(i)  = 0.0; // number of females ready to lay eggs by mating group
	Af_mated_gm.at(i)  = 0.0;
	Af_mated_rel.at(i) = 0.0;
	 
	Am_gm.at(i) = 0.0; // number of adults
	Am_gm_release.at(i) = 0.0;

	Am_freq.at(i) = 0.0; // frequency of males ready to mate
    Am_gm_freq.at(i) = 0.0;
    Am_gm_release_freq.at(i) = 0.0;
	
	// With 2 columns, one for each pupal stage (0-2 days)
	for (int j=0; j<2; j++){ 
		Pm[i][j] = 0.0; // number of pupae
		Pf[i][j] = 0.0; 
		Pm_gm[i][j] = 0.0; 
	}
}

// By day - starts at 1 then 0 for all other days (wildtype)
Am.at(0) = 1.0;
Af.at(0) = 1.0;
for (int i=1; i<maxtime; i++) {
	Am.at(i) = 0.0; // number of adults
	Af.at(i) = 0.0; 
}


//------------------------------------- RUN MODEL function: Time Loop -------------------------------------

// The whole code repeats itself every day from day 1 to 2500
for (int time1=1; time1<maxtime; time1++){

	// Track how many active cohorts the code needs to loop over
	max_cohort=-1; // find most recent cohort that has already hatched 
	for (int i=0; i<no_cohorts; i++) if (hdate.at(i)<=time1) max_cohort = i;
	min_cohort=no_cohorts-1; // find oldest cohort still worth tracking (ignore if older than 100 days [cohort_age_max])
	for (int i=no_cohorts-1; i>=0; i--) if (hdate.at(i)>=time1 - cohort_age_max) min_cohort = i;	


	//----------------------------------- RUN MODEL function: Cohort Loop ----------------------------------

	// The whole code repeats itself for every active cohort
	for (int cohort=min_cohort; cohort<=max_cohort; cohort++){

		//-------------------------------------------------- Eggs ------------------------------------------

		//------------------------------------------------- Larvae -----------------------------------------

		// Calculate number of yesterday's larvae still alive
		Lm_cohort[time1][cohort] = Lm_cohort[time1-1][cohort] * survL_wt; // larvae alive yesterday * survival rate (males)
		Lf_cohort[time1][cohort] = Lf_cohort[time1-1][cohort] * survL_wt; 

		Lm_gm_cohort[time1][cohort] = Lm_gm_cohort[time1-1][cohort] * survL_gm;
		Lf_gm_cohort[time1][cohort] = Lf_gm_cohort[time1-1][cohort] * survL_gm;
		
		// Check if today is a hatch date - if yes, calculate number of new larvae hatching today
		if (time1 == hdate.at(cohort)){

			// Calculate start/end point of look-back window 
			start = time1-tlagh-tlagg-tlagf-tlagp-tlagl; // 33 days ago, when female parent larvae started developing (note density-dependence only applied to mothers)
			end = time1-tlagh-tlagg-tlagf-tlagp; // 12 days ago, when female parent larvae emerged as pupae

			// Calculate average larval density experienced by this cohort (parents)
			L_avg.at(cohort) = 0; 

			// Case 1: full window within simulation (i.e. start ≥ first hatch date)
			if (start>=hdate.at(1)){
				// Total larvae (21 days)/21 days = larval density experienced
				for (int j=start; j<end; j++) L_avg.at(cohort) += Lm.at(j) + Lf.at(j) + Lm_gm.at(j) + Lf_gm.at(j);
					L_avg.at(cohort) = L_avg.at(cohort)/tlagl; 
			}

			// Case 2: window partially overlaps with simulation (i.e. start < first hatch date, but still after day 0)
			else if (start>=0) {
				// Case 2a: end < first hatch date (i.e. full window outside simulation) - I don't think this can happen, first hatch date is day 4
				if (end<=hdate.at(1)) L_avg.at(cohort) = 60; // default average density (60)
				// Case 2b: end > first hatch date (i.e. window partially within simulation)
				if (end>hdate.at(1)) {
					for (int j=start; j<=hdate.at(1); j++) L_avg.at(cohort) += 60; // section before first hatch date = default (60 larvae per day)
					for (int j=hdate.at(1)+1; j<=end; j++) L_avg.at(cohort) += Lm.at(j) + Lf.at(j) + Lm_gm.at(j) + Lf_gm.at(j); // section after first hatch date = total larvae per day (real data)
					L_avg.at(cohort) = L_avg.at(cohort)/tlagl;
				} 
			}

			// Case 3: window starts before day 0
			else L_avg.at(cohort) = 60; // default average density (60)

			// Calculate per-capita female fecundity of mothers 
			L_avg_f = L_avg.at(cohort);
			if (L_avg_f > 0.0){
				lambda1 = lambda_intc + lambda_alpha * log(L_avg_f); // calculate fecundity using log-linear function (higher density experienced = lower fecundity)

				// Limit fecundity to biologically realistic limits
				if (lambda1<lambda_min) lambda1 = lambda_min; // female cannot lay less than 0.5 eggs per day
				if (lambda1>lambda_max) lambda1 = lambda_max; // female cannot lay more than 14 eggs per day
			} else {
				lambda1 = lambda_max; // no crowding = maximum fecundity
			}

			// Store fecundity for this cohort
			lambda.at(cohort) = lambda1;

			// Calculate number of new larvae hatching today
 			if (time1>tlagh+tlagg+tlagm){ // prevents hatching from happening until enough time has passed for first adults in simulation to lay eggs
				
				Lm_cohort[time1][cohort] = sex_ratio 	 * lambda1 * Af_mated_wt.at(time1-tlagh-tlagg) * pow(survA_wt, tlagg);
				Lf_cohort[time1][cohort] = (1-sex_ratio) * lambda1 * Af_mated_wt.at(time1-tlagh-tlagg) * pow(survA_wt, tlagg);	

				Lm_gm_cohort[time1][cohort] = sex_ratio 	* lambda1 * Af_mated_rel.at(time1-tlagh-tlagg) * pow(survA_gm, tlagg);	
				Lf_gm_cohort[time1][cohort] = (1-sex_ratio) * lambda1 * Af_mated_rel.at(time1-tlagh-tlagg) * pow(survA_gm, tlagg);
				
			}

		} // End loop for newly hatching larvae	


		//------------------------------------------------- Pupae ------------------------------------------

		// Calculate mean and std dev of larval development time for each cohort (expected)
		if (time1>hdate.at(cohort)+4 && emerg_flag.at(cohort)==0){ // only run if ≥5 days since cohort hatched (min dev time) and <2 pupae have emerged
		
			L_cum2     = 0;
			denom 	   = 0;
			L_avg_gam1 = 0; 

			// Calculate weighted average larval density experienced by each cohort
			for (int etime = hdate.at(cohort)+5; etime <=time1; etime++){ // for each emergence day since 5 days since cohort hatched

				L_cum = 0;

				// Calculate cumulative number of larvae from hatch to emergence day
				for (int k=hdate.at(cohort); k<etime; k++) {
					L_cum += Lm.at(k) + Lf.at(k) + Lm_gm.at(k) + Lf_gm.at(k);
				}

				// Calculate total emerged pupae (by day and cohort)
				P_emerg_cohort = Pm_emerg_cohort[etime][cohort] + Pf_emerg_cohort[etime][cohort] + Pm_gm_emerg_cohort[etime][cohort];

				// Calculate weighted average
				L_cum2 += P_emerg_cohort * L_cum/(etime-hdate.at(cohort)); // number emerged * daily average
				
				// Store cumulative number of emerged pupae
				denom += P_emerg_cohort;
			}

			// Check when to use weighted average
			if (denom>0) L_avg_gam1 = L_cum2/denom; // yes - pupae have started emerging
			else L_avg_gam1 = L_cum/(time1 - hdate.at(cohort)); // no - no pupae have emerged yet (use unweighted)

			// Calculate expected mean of larval development time
			dt_mean_gam.at(cohort) = dt_mean_intc + dt_mean_alpha * pow(L_avg_gam1, dt_mean_exp); // calculate mean using power-law formula (higher density experienced = longer development time)
			if (dt_mean_gam.at(cohort)<0) dt_mean_gam.at(cohort)=0; // floor at 0 to avoid nonsensical values
			
			// Calculate expected std dev of larval development time
			dt_std_gam.at(cohort) = dt_std_intc + dt_std_alpha * pow(L_avg_gam1, dt_std_exp); // calculate std dev using power-law formula (higher density experienced = more variable development time)
			
			// Limit mean and std dev to biologically realistic limits
			if (L_avg_gam1>L_max) {dt_mean_gam.at(cohort) = dt_mean_max; dt_std_gam.at(cohort) = dt_std_max;} // if density exceeded (L_max = 4200) then mean and std dev set to max allowed values
			if (dt_std_gam.at(cohort) < dt_std_min) dt_std_gam.at(cohort) = dt_std_min; // if std dev is less than minimum (1) then set to minimum

			// Flag - if ≥2 pupae have emerged, stop updating mean and std dev
			if (denom>1) {
				emerg_flag.at(cohort)=1;
				L_avg_gam.at(cohort) = L_avg_gam1;
			}

		} // End loop for calculating mean and std dev of larval development time (expected)
		
		// Calculate number of larvae emerging as pupae
		if (time1>hdate.at(cohort)+4){

			// Convert mean and std dev to scale and shape parameters
			dt_shp = pow(dt_mean_gam.at(cohort)/dt_std_gam.at(cohort),2);
			dt_scl = dt_mean_gam.at(cohort)/dt_shp;						
			time_lag = time1 - hdate.at(cohort);

			// Calculate probability of emerging today
			if (time_lag>5 && dt_mean_gam.at(cohort)>0 && dt_shp>0 && dt_scl>0 && dt_shp<200 && dt_scl<200){				 
				prob = gsl_cdf_gamma_P(time_lag-5,dt_shp,dt_scl) - gsl_cdf_gamma_P(time_lag-5-1,dt_shp,dt_scl);
			}
			else prob=0;

			// Safety checks
			if (dt_mean_gam.at(cohort)==0){
				dt_shp = 9.0;
				dt_scl = 0.2;
				prob = gsl_cdf_gamma_P(time_lag-4,dt_shp,dt_scl) - gsl_cdf_gamma_P(time_lag-4-1,dt_shp,dt_scl);
			}
			
			// Convert probability to numbers
			Pm_emerge_today = Lm_cohort[time1-1][cohort] * prob/non_emerg_prob.at(cohort);
			Pf_emerge_today = Lf_cohort[time1-1][cohort] * prob/non_emerg_prob.at(cohort);
			Pm_gm_emerge_today = Lm_gm_cohort[time1-1][cohort] * prob / non_emerg_prob.at(cohort);
			Pf_gm_emerge_today = Lf_gm_cohort[time1-1][cohort] * prob / non_emerg_prob.at(cohort);
			
			// Sanity checks
			if (Pm_emerge_today > Lm_cohort[time1-1][cohort]){Pm_emerge_today = Lm_cohort[time1-1][cohort]; prob=1;} // prevents more pupae emerging than larvae exist
			if (Pf_emerge_today > Lf_cohort[time1-1][cohort]){Pf_emerge_today = Lf_cohort[time1-1][cohort]; prob=1;}
			if (Pm_gm_emerge_today > Lm_gm_cohort[time1-1][cohort]) {Pm_gm_emerge_today = Lm_gm_cohort[time1-1][cohort]; prob=1;}
			if (Pf_gm_emerge_today > Lf_gm_cohort[time1-1][cohort]) {Pf_gm_emerge_today = Lf_gm_cohort[time1-1][cohort]; prob=1;}
			
			if (Pm_emerge_today<0){Pm_emerge_today=0; prob=0;} // prevents negative emergence
			if (Pf_emerge_today<0){Pf_emerge_today=0; prob=0;}
			if (Pm_gm_emerge_today<0){Pm_gm_emerge_today=0; prob=0;}
			if (Pf_gm_emerge_today<0){Pf_gm_emerge_today=0; prob=0;}
			
			if (prob<0) prob=0; // further prevents negative probability of emergence

			// Updates based on calculations
			non_emerg_prob.at(cohort) *= (1-prob/non_emerg_prob.at(cohort)); // update probability of not emerging
	
			Lm_cohort[time1][cohort] -= Pm_emerge_today; // remove emerging larvae from larvae compartment
            Lf_cohort[time1][cohort] -= Pf_emerge_today;
			Lm_gm_cohort[time1][cohort] -= Pm_gm_emerge_today;
			Lf_gm_cohort[time1][cohort] -= Pf_gm_emerge_today;
			
			Pm_emerg_cohort[time1][cohort] = Pm_emerge_today; // store emerging pupae by day into each cohort
			Pf_emerg_cohort[time1][cohort] = Pf_emerge_today;
			Pm_gm_emerg_cohort[time1][cohort]  = Pm_gm_emerge_today;
			Pf_gm_emerg_cohort[time1][cohort] = Pf_gm_emerge_today; // stored to monitor how many larvae are removed but not added to pupal compartment below

			P_emerg.at(cohort) 	  += Pm_emerge_today + Pf_emerge_today; // store total number of WT pupae that emerged by cohort
			P_gm_emerg.at(cohort) += Pm_gm_emerge_today;

			Pm[time1][0] 	+= Pm_emerge_today; // add emerging pupae to pupae compartment (males)
			Pf[time1][0] 	+= Pf_emerge_today;
			Pm_gm[time1][0] += Pm_gm_emerge_today;

			dt_mean.at(cohort) += (time1-hdate.at(cohort)) * (Pm_emerge_today + Pf_emerge_today + Pm_gm_emerge_today); // accumulate mean development time

		} // End loop for calculating number of larvae emerging as pupae

		// Store total number of larvae
		Lm.at(time1) += Lm_cohort[time1][cohort]; // add up total larvae from each cohort
		Lf.at(time1) += Lf_cohort[time1][cohort];

		Lm_gm.at(time1) += Lm_gm_cohort[time1][cohort];
        Lf_gm.at(time1) += Lf_gm_cohort[time1][cohort];

	} // End loop for calculations for each active cohort


	//--------------------------------------------- Adults --------------------------------------------

	// Calculate number of pupae alive today
	Pm[time1][1] = Pm[time1-1][0] * survA_wt; // move into second pupal day
	Pf[time1][1] = Pf[time1-1][0] * survA_wt;

	Pm_gm[time1][1] = Pm_gm[time1-1][0] * survA_gm;

	// Calculate number of adults alive today
	Am[time1] = Am[time1-1] * survA_wt;
	Af[time1] = Af[time1-1] * survA_wt;

	Am_gm[time1] 		 = Am_gm[time1-1] 		  * survA_gm;
	Am_gm_release[time1] = Am_gm_release[time1-1] * survA_rel;

	// Calculate number of pupae graduating into adults today
	Am[time1] += Pm[time1-1][1] * survA_wt;
	Af[time1] += Pf[time1-1][1] * survA_wt;

	Am_gm[time1] += Pm_gm[time1-1][1] * survA_gm;

	// Release event
	for (int r = 0; r < release_period; r++){ // number of releases/weeks
    if (time1 == release_day + r * 7){ // check if today is a release day
        Am_gm_release.at(time1) += release_size; // same release size every time
    }
	}

	// Calculate observed male frequencies (different from mating frequencies due to time lag and mating competitiveness)
	double Am_total_obs = Am[time1] + Am_gm[time1] + Am_gm_release[time1];

	if (Am_total_obs > 0.0){
    	Am_freq[time1]             = Am[time1]             / Am_total_obs;
    	Am_gm_freq[time1]          = Am_gm[time1]          / Am_total_obs;
    	Am_gm_release_freq[time1]  = Am_gm_release[time1]  / Am_total_obs;
	} else {
    	Am_freq[time1]             = 0.0;
    	Am_gm_freq[time1]          = 0.0;
    	Am_gm_release_freq[time1]  = 0.0;
	}


	//--------------------------------------------- Mating -------------------------------------------- 


	// Calculate mating frequencies by male type
	double Am_freq_wt, Am_freq_gm, Am_freq_rel;
	if (time1>=tlagm){ // adult males must be at least 2 days old to be able to mate
    	
		// Number of adult males ready to mate (male adults alive 2 days ago still alive today)
		double Am_mate = comp_wt * (Am[time1-tlagm] * pow(survA_wt, tlagm));
    	double Am_gm_mate = comp_gm * (Am_gm[time1-tlagm] * pow(survA_gm, tlagm));
    	double Am_gm_release_mate = comp_rel * (Am_gm_release[time1-tlagm] * pow(survA_rel, tlagm));

		// Total males ready to mate
		double Am_total = Am_mate + Am_gm_mate + Am_gm_release_mate;

		// Calculate mating frequencies
		if (Am_total > 0.0){
    		Am_freq_wt  = Am_mate             / Am_total;
    		Am_freq_gm  = Am_gm_mate          / Am_total;
    		Am_freq_rel = Am_gm_release_mate  / Am_total;
		} else {
    		Am_freq_wt  = 0.0;
    		Am_freq_gm  = 0.0;
    		Am_freq_rel = 0.0;
		}
	}

	// Calculate number of ovipositing females alive today 
	Af_mated_wt[time1]  = Af_mated_wt[time1-1]  * survA_wt;
	Af_mated_gm[time1]  = Af_mated_gm[time1-1]  * survA_wt;
	Af_mated_rel[time1] = Af_mated_rel[time1-1] * survA_wt;

	// Calculate new adult females ready to lay eggs 
	if (time1>=tlagm+tlagg) { // adult males need to be ready to mate, then females need to be ready to lay eggs after mating
		double Af_mate = Af[time1-tlagf] * pow(survA_wt, tlagf); // females ready to mate (i.e. females alive 5 days ago that were still alive the next day)
		double Af_unmated = Af_mate - Af_mated_wt[time1] - Af_mated_gm[time1] - Af_mated_rel[time1]; // unmated females ready to mate 4 days ago 
		Af_unmated = max(0.0, Af_unmated);

		Af_mated_wt[time1] += Af_unmated * Am_freq_wt;
		Af_mated_gm[time1] += Af_unmated * Am_freq_gm;
		Af_mated_rel[time1] += Af_unmated * Am_freq_rel;
	}
	
} // End time loop

//--------------------------------------------- Model Outputs ------------------------------------------

// Calculate observed mean and standard deviation of larval development time for each cohort
for (int i=0; i<no_cohorts; i++) {
	if (P_emerg.at(i) + P_gm_emerg.at(i) > 0) dt_mean.at(i) /= (P_emerg.at(i) + P_gm_emerg.at(i)); // mean
	else dt_mean.at(i)=0;
	for (int j=0; j<maxtime; j++) { // accumulate squared deviations
		if (j>hdate.at(i)+4){
			P_emerg_cohort = Pf_emerg_cohort[j][i] + Pm_emerg_cohort[j][i] + Pm_gm_emerg_cohort[j][i];
			dt_std.at(i) +=(j-hdate.at(i)-dt_mean.at(i))*(j-hdate.at(i)-dt_mean.at(i))*P_emerg_cohort;
		}
   	}
	if (P_emerg.at(i) + P_gm_emerg.at(i) > 0) dt_std.at(i) = sqrt(dt_std.at(i) / (P_emerg.at(i) + P_gm_emerg.at(i))); // std dev
	else dt_std.at(i)=0;
}

// Calculate total number of pupae that emerged each day
for (int i=0; i<maxtime; i++) {
	for (int j=0; j<no_cohorts; j++) {
		Pm_emerg.at(i) += Pm_emerg_cohort[i][j];
        Pf_emerg.at(i) += Pf_emerg_cohort[i][j];

		Pm_gm_emerg.at(i) += Pm_gm_emerg_cohort[i][j];
	}
}

// Write per-capita female fecundity for each cohort to the output file
for (int i=0; i<no_cohorts; i++) lambda_out << lambda.at(i) << endl;

// Calculate summary statistics
double dt_mean_sum = 0, dt_std_sum = 0, lambda_sum = 0;
for (int i=300; i<800; i++){ // averaged across cohorts 300 to 800
	dt_mean_sum+=dt_mean.at(i); // mean
	dt_std_sum+=dt_std.at(i); // std dev
	lambda_sum += lambda.at(i); // lambda
} 

// Print summary statistics
cout << endl;
cout << "Average mean development time of larvae " << dt_mean_sum/500 << endl;
cout << "Average per-capita female fecundity " << lambda_sum/500 << endl;
cout << "Average standard deviation of the development time of larvae " << dt_std_sum/500 << endl;

return 0;			

}




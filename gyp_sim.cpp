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
void read_data(ofstream& cohort_means_out, ofstream& cohort_stds_out, ofstream& lambda_out, ofstream& Lm_out, ofstream& Lf_out, ofstream& Pm_emerg_out, ofstream& Pf_emerg_out, ofstream& Am_out, ofstream& Af_out, ofstream& Lm_gm_out, ofstream& Lf_gm_out, ofstream& Pm_gm_emerg_out, ofstream& Am_gm_out, ofstream& Am_gm_release_out);
double Run_Model(ofstream& lambda_out);

// Declare global variables
const int no_cohorts = 1000; // number of mosquito cohorts
const int maxtime = 2500; // simulation length (days; ~7 years)
double costA, costL;
int release_size, release_day;

// Declare global data containers
vector<float> dt_mean(no_cohorts), dt_std(no_cohorts), A_ovipos(maxtime);
vector<double> Pm_emerg(maxtime), Lm(maxtime), Am(maxtime), Pf_emerg(maxtime), Lf(maxtime), Af(maxtime); // wildtype
Matrix_Double Pm(maxtime, Row_Double(3)), Pf(maxtime, Row_Double(3)); // wildtype
vector<double> Pm_gm_emerg(maxtime), Lm_gm(maxtime), Am_gm(maxtime), Lf_gm(maxtime), Am_gm_release(maxtime); // GM
Matrix_Double Pm_gm(maxtime, Row_Double(3)); // GM
vector<int> hdate(no_cohorts);


//------------------------------------------- MAIN function -------------------------------------------

int main(){

// Declare input and output files
ofstream cohort_means_out, cohort_stds_out, lambda_out;
ofstream Lm_out, Lf_out, Pm_emerg_out, Pf_emerg_out, Am_out, Af_out;
ofstream Lm_gm_out, Lf_gm_out, Pm_gm_emerg_out, Am_gm_out, Am_gm_release_out;

// Calls function: read_data
read_data(cohort_means_out, cohort_stds_out, lambda_out, Lm_out, Lf_out, Pm_emerg_out, Pf_emerg_out, Am_out, Af_out, Lm_gm_out, Lf_gm_out, Pm_gm_emerg_out, Am_gm_out, Am_gm_release_out);

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

return 0;

}

//------------------------------------------- READ DATA function -------------------------------------------

// Call in input and output files
void read_data(ofstream& cohort_means_out, ofstream& cohort_stds_out, ofstream& lambda_out, ofstream& Lm_out, ofstream& Lf_out, ofstream& Pm_emerg_out, ofstream& Pf_emerg_out, ofstream& Am_out, ofstream& Af_out, ofstream& Lm_gm_out, ofstream& Lf_gm_out, ofstream& Pm_gm_emerg_out, ofstream& Am_gm_out, ofstream& Am_gm_release_out){

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

cin >> cohort_means_file >> cohort_stds_file >> lambda_file >> Lm_file >> Lf_file >> Pm_emerg_file >> Pf_emerg_file >> Am_file >> Af_file >> Lm_gm_file >> Lf_gm_file >> Pm_gm_emerg_file >> Am_gm_file >> Am_gm_release_file;

// Print file names to console
cout << " cohort_means_file " << cohort_means_file << " cohort_stds_file " << cohort_stds_file << " lambda_file " << lambda_file << endl;
cout << " Lm_file " << Lm_file << " Lf_file " << Lf_file << " Pm_emerg_file " << Pm_emerg_file << " Pf_emerg_file " << Pf_emerg_file << " Am_file " << Am_file << " Af_file " << Af_file << endl;
cout << " Lm_gm_file " << Lm_gm_file << " Lf_gm_file " << Lf_gm_file << " Pm_gm_emerg_file " << Pm_gm_emerg_file << " Am_gm_file " << Am_gm_file << " Am_gm_release_file " << Am_gm_release_file << endl;

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

// Read in key parameter values in inits file and print to check correctly read
cin >> costA;
cin >> costL;
cin >> release_size;
cin >> release_day;

cout << "costA " << costA << " (additional density-INdependent daily mortality experienced by adults in the field environment) " << endl;
cout << "costL " << costL << " (additional density-INdependent daily mortality experienced by larvae in the field environment) " << endl;
cout << "Release size: " << release_size << endl;
cout << "Release day: " << release_day << endl;
}

//------------------------------------------- RUN MODEL function -------------------------------------------

// Call in input and output files
double Run_Model(ofstream& lambda_out){

// Time lags
int tlagh = 5; // oviposition to egg hatching
int tlagl = 21; // time over which larval density is averaged
int tlagp = 2; // pupal development time
int tlagg = 4; // emerging as adults to being ready to lay eggs (females)
int tlagm = 2; // emerging as adults and being ready to mate (males)

// Mortality
double survA = 1-0.03-costA; // daily survival rate of adults - base mortality is 3% per day, plus any additional cost (costA)
double survL = 0.95-costL; // daily survival rate of larvae
double survLf_gm = 0.001; // daily survival rate of GM female larvae

// Fecundity (log-linear relationship)
double lambda_intc = 28.0; // intercept 
double lambda_alpha = -3.3; // slope 

double lambda_max = 14; // a female can lay at most 14 eggs per day
double lambda_min = 0.5; // a female lays at least 0.5 eggs per day even under very crowded conditions

// Pupal emergence (gamma distribution)
double dt_mean_intc = 1.8; // mean larval development time - intercept
double dt_mean_alpha = 0.536; // mean larval development time - slope
double dt_mean_exp = 0.533; // mean larval development time - exponent

double dt_std_intc = 0.22; // std dev larval development time - intercept
double dt_std_alpha = 0.0168; // std dev larval development time - slope
double dt_std_exp = 0.867; // std dev larval development time - exponent

double dt_mean_max = 60; // if crowding is extreme, mean development time is capped at 60 days
double dt_std_max = 40; // if crowding is extreme, std dev development time is capped at ±60 days from mean
double dt_std_min = 1.0; // if crowding is minimal, std dev development time is capped at ±1 day from mean

int L_max = 4200; // the larval density above which the above caps on development time kick in

// Other
double sex_ratio = 0.5; 
int cohort_age_max = 100; // cohorts older than 100 days are ignored

// Temporary variables
int max_cohort, min_cohort; // active cohorts

double Am_total, Am_gm_release_freq, Am_gm_freq, Am_freq; // mating equation
int t_mate;

double L_avg_f, lambda1; // hatching larvae/fecundity
int start, end; 

double L_cum, L_cum2, denom, L_avg_gam1, P_emerg_cohort; // expected mean & std dev
int time_lag;

double dt_shp, dt_scl, prob, L_cohort, P_emerge_today, Pm_gm_emerge_today; // pupal emergence


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

// By cohort 
vector<double> non_emerg_prob(no_cohorts, 1.0); // cumulative probability that a larva has not yet emerged as a pupa (1 = 100%; decreases over time)
vector<int> emerg_flag(no_cohorts, 0); // flag to track if enough pupae have emerged to start calculating statistics
vector<double> L_avg_gam(no_cohorts, 0.0); // running average larval density experienced (estimated)
vector<double> L_avg(no_cohorts, 0.0); // average larval density experienced (final)
vector<double> dt_mean_gam(no_cohorts, 999.0); // mean larval development time (gamma distribution)
vector<double> dt_std_gam(no_cohorts, 1.0); // std dev of larval development time (gamma distribution)
vector<double> P_emerg(no_cohorts, 0.0); // total number of pupae that emerged
vector<double> P_gm_emerg(no_cohorts, 0.0);
vector<double> L_hatch(no_cohorts, 0.0); // number of larvae hatching into each cohort (i.e. initial population)
vector<double> lambda(no_cohorts, 0.0); // per-capita female fecundity


// Initialise global data containers
// By cohort
for (int i=0; i<no_cohorts; i++){
	dt_mean.at(i) = 0.0; // mean larval development time (final)
	dt_std.at(i) = 0.0; // std dev of larval development time (final)
}

// By day
for (int i=0; i<maxtime; i++) {
	// Wildtype
	A_ovipos.at(i) = 0.0; // number of reproductively active females
	Lm.at(i) = 0.0; // number of larvae (males)
	Lf.at(i) = 0.0; 
	Pm_emerg.at(i) = 0.0; // number of new pupae (males)
	Pf_emerg.at(i) = 0.0; 
	// GM
	Am_gm.at(i) = 0.0;
	Am_gm_release.at(i) = 0.0;
	Lm_gm.at(i) = 0.0; 
	Lf_gm.at(i) = 0.0; 
	Pm_gm_emerg.at(i) = 0.0; 

	// With 3 columns, one for each pupal stage (0-2 days)
	for (int j=0; j<3; j++){ 
		// Wildtype
		Pm[i][j] = 0.0; // number of pupae (males)
		Pf[i][j] = 0.0; 
		// GM
		Pm_gm[i][j] = 0.0; // number of pupae (males)
	}
}

// By day - starts at 1 then 0 for all other days (wildtype)
Am.at(0) = 1.0;
Af.at(0) = 1.0;
for (int i=1; i<maxtime; i++) {
	Am.at(i) = 0.0; // number of adults (males)
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

	// Calculate male type frequency for mating equation
	t_mate = tlagh+tlagg; // time of mating (as soon as females emerge +1 day, since calculated based on yesterday's counts)

	Am_freq            = 1.0;
	Am_gm_freq         = 0.0;
	Am_gm_release_freq = 0.0;

	if (time1>tlagh+tlagg){ // wait until enough time has passed for second generation adults in simulation to emerge
		Am_total = Am.at(time1-t_mate) + Am_gm.at(time1-t_mate) + Am_gm_release.at(time1-t_mate);

		if (Am_total > 0.0){
			Am_gm_release_freq = Am_gm_release.at(time1-t_mate) / Am_total; // % released GM males
			Am_gm_freq = Am_gm.at(time1-t_mate) / Am_total; // % second generation GM males
			Am_freq = 1.0 - Am_gm_freq - Am_gm_release_freq; // % wildtype males
		}
	}


	//----------------------------------- RUN MODEL function: Cohort Loop ----------------------------------

	// The whole code repeats itself for every active cohort
	for (int cohort=min_cohort; cohort<=max_cohort; cohort++){

		//------------------------------------------------- Larvae -----------------------------------------

		// Calculate number of yesterday's larvae still alive
		Lm_cohort[time1][cohort] = Lm_cohort[time1-1][cohort] * survL; // larvae alive yesterday * survival rate (males)
		Lf_cohort[time1][cohort] = Lf_cohort[time1-1][cohort] * survL; 

		Lm_gm_cohort[time1][cohort] = Lm_gm_cohort[time1-1][cohort] * survL;
		Lf_gm_cohort[time1][cohort] = Lf_gm_cohort[time1-1][cohort] * survLf_gm;
		
		// Check if today is a hatch date - if yes, calculate number of new larvae hatching today
		if (time1 == hdate.at(cohort)){

			// Calculate start/end point of look-back window 
			start = time1-tlagh-tlagg-tlagp-tlagl; // 32 days ago, when parent larvae started developing
			end = time1-tlagh-tlagg-tlagp; // 11 days ago, when parent larvae emerged as pupae

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
 			if (time1>tlagh+tlagg){ // prevents hatching from happening until enough time has passed for first adults in simulation to lay eggs (9 days)
				L_hatch.at(cohort) = lambda1 * A_ovipos.at(time1-tlagh); // fecundity * number of adult females ready to mate
			}

			// Store newly hatched larvae as starting larval population for this new cohort
			Lm_cohort[time1][cohort] = sex_ratio * L_hatch.at(cohort) * Am_freq;	
			Lf_cohort[time1][cohort] = (1-sex_ratio) * L_hatch.at(cohort) * Am_freq;	

			Lm_gm_cohort[time1][cohort] = sex_ratio * L_hatch.at(cohort) * Am_gm_release_freq;	
			Lf_gm_cohort[time1][cohort] = (1-sex_ratio) * L_hatch.at(cohort) * Am_gm_release_freq;

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
			L_cohort = Lm_cohort[time1-1][cohort] + Lf_cohort[time1-1][cohort]; // total wildtype larvae

			P_emerge_today = L_cohort * prob/non_emerg_prob.at(cohort); // total
			Pm_gm_emerge_today = Lm_gm_cohort[time1-1][cohort] * prob / non_emerg_prob.at(cohort);
			
			// Sanity checks
			if (P_emerge_today > L_cohort){P_emerge_today = L_cohort; prob=1;} // prevents more pupae emerging than larvae exist
			if (Pm_gm_emerge_today > Lm_gm_cohort[time1-1][cohort]) {Pm_gm_emerge_today = Lm_gm_cohort[time1-1][cohort]; prob = 1;}
			
			if (P_emerge_today<0){P_emerge_today=0; prob=0;} // prevents negative emergence
			if (Pm_gm_emerge_today<0){Pm_gm_emerge_today=0; prob=0;}
			
			if (prob<0) prob=0; // further prevents negative probability of emergence

			// Updates based on calculations
			non_emerg_prob.at(cohort) *= (1-prob/non_emerg_prob.at(cohort)); // update probability of not emerging
	
			Lm_cohort[time1][cohort] -= sex_ratio * P_emerge_today; // remove emerging larvae from larvae compartment
            Lf_cohort[time1][cohort] -= (1-sex_ratio) * P_emerge_today;
			Lm_gm_cohort[time1][cohort] -= Pm_gm_emerge_today;
			
			Pm_emerg_cohort[time1][cohort] = sex_ratio * P_emerge_today; // store emerging pupae by day into each cohort
			Pf_emerg_cohort[time1][cohort] = (1-sex_ratio) * P_emerge_today;
			Pm_gm_emerg_cohort[time1][cohort]  = Pm_gm_emerge_today;

			P_emerg.at(cohort) += P_emerge_today; // store total number of pupae that emerged by cohort
			P_gm_emerg.at(cohort) += Pm_gm_emerge_today;

			Pm[time1][0] += sex_ratio * P_emerge_today; // add emerging pupae to pupae compartment (males)
			Pf[time1][0] += (1-sex_ratio) * P_emerge_today;
			Pm_gm[time1][0] += Pm_gm_emerge_today;

			dt_mean.at(cohort) += (time1-hdate.at(cohort)) * (P_emerge_today + Pm_gm_emerge_today); // accumulate mean development time

		} // End loop for calculating number of larvae emerging as pupae

		// Store total number of larvae
		Lm.at(time1) += Lm_cohort[time1][cohort]; // add up total larvae from each cohort
		Lf.at(time1) += Lf_cohort[time1][cohort];

		Lm_gm.at(time1) += Lm_gm_cohort[time1][cohort];
        Lf_gm.at(time1) += Lf_gm_cohort[time1][cohort];

	} // End loop for calculations for each active cohort


	//--------------------------------------------- Adults --------------------------------------------

	// Calculate number of pupae alive today
	Pm[time1][1] = Pm[time1-1][0] * survA; // move into second pupal day
	Pf[time1][1] = Pf[time1-1][0] * survA;

	Pm_gm[time1][1] = Pm_gm[time1-1][0] * survA;

	// Calculate number of adults alive today
	Am[time1] = Am[time1-1] * survA;
	Af[time1] = Af[time1-1] * survA;

	Am_gm[time1] = Am_gm[time1-1] * survA;
	Am_gm_release[time1] = Am_gm_release[time1-1] * survA;

	// Calculate number of pupae graduating into adults today
	Am[time1] += Pm[time1-1][1] * survA;
	Af[time1] += Pf[time1-1][1] * survA;

	Am_gm[time1] += Pm_gm[time1-1][1] * survA;

	// Release event (only one release currently)
	if (time1 == release_day){
    	Am_gm_release.at(time1) += release_size;
	}

	// Calculate number of adult females ready to lay eggs
	if (time1>=tlagg) { // adult females must be at least 4 days old to be ready to lay eggs
		A_ovipos[time1] = Af[time1-tlagg] * pow(survA,tlagg-1); // females that emerged as adults 4 days ago that are still alive (adjusted survival ^3 days)
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




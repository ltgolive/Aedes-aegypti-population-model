//gyp_sim.cpp was developed by Dr Penelope A. Hancock
//The development of gyp_sim.cpp is described in Hancock et al. 2016, "Predicting Wolbachia
//invasion dynamics in Aedes aegypti populations using models of density-dependent
//demographic traits", BMC Biology.

// This section imports libraries needed (like in R) -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#include <iostream> // reading/writing files
#include <fstream> // reading/writing files
#include <cmath> // math functions
#include <algorithm>
#include <vector>
#include <gsl/gsl_rng.h> // statistical distributions
#include <gsl/gsl_randist.h> // statistical distributions
#include <gsl/gsl_cdf.h> // statistical distributions
#include <float.h>
#include <stdio.h>
#include <string>
#include <sstream>
#include <iomanip>

using namespace std; // standard library names (avoid writing std:: before a common library name every time (e.g. std::string can just be string))
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Type definitions (data structures) and globals

// TYPE DEFINITIONS
// Type definitions are basically schortcuts/renames for data structures (like vectors, or matrix_double - this is a list of lists/table)
typedef vector<int> Row_Int; // you can write "Row_Int row" instead of "std::vector<int> row" for a row of integers
typedef vector<Row_Int> Matrix_Int; // uses the above Row_Int, so creates a matrix/table of integers

typedef vector<double> Row_Double; // same as above but using double (with decimals) instead of integers (round numbers)
typedef vector<Row_Double> Matrix_Double;
//NOTE: first line could be replaced in modern C++ with --> using Row_Int = vector<int> (same with other lines)

// They are essentially a heads-up to the compiler saying "these functions exist (or will exist) and here's what they look like — I'll give you the full details later."
// In C++, before you can use a function, the compiler needs to have seen at least its declaration.
// void read_data(...) — declares a function called read_data that returns nothing (void) for output (and 1 input) files (i.e. they are currently blank)
// double Run_Model(...) — declares a function called Run_Model that returns a double (0 for now) for 3 output files and a vector.
// NOTE: I don't think max_dt is ever used?
//some of these vectors that are getting passed around could just be global...
void read_data(ofstream& cohort_means_out, ofstream& cohort_means_wolb_out, ofstream& cohort_stds_out, ofstream& cohort_stds_wolb_out, ofstream& mu_p_out,ofstream& mu_p_wolb_out, ofstream& L_out, ofstream& L_wolb_out, ofstream& A_out, ofstream& A_wolb_out, ofstream& freqA2_out, ofstream& lambda_out, ifstream& rel_siz_in);
double Run_Model(const vector<int> max_dt, ofstream& mu_p_out,  ofstream& freqA2_out, ofstream& lambda_out);

// GLOBALS
// These are constant/fixed integers (i.e. globals)
const int no_cohorts = 1000; // number of mosquito cohorts (egg batches) being tracked. Hmmm... not sure I like specifying this
const int no_pdays = 170; // "number of pupal days". It sizes arrays L_data, pupae, pdate and pupae_sampled, none of which are then used. Remove?
const int no_sdays = 169; // number of sample days - not used? remove?
const int maxtime = 2500; // simulation length (days; ~7 years)

// Global variables (data containers)
// Declares a long list (float/decimal vectors) for the mean and std dev time for each of the 1000 cohorts, and other variables over the max time - I think L_data isn't used
// e.g. one mean dev time per cohort; number of pupae emerging per day (mu_p)
vector<float> mean_dt(no_cohorts), mean_dt_wolb(no_cohorts), std_dt(no_cohorts), std_dt_wolb(no_cohorts), max_dt(no_cohorts), mu_p(maxtime), mu_p_wolb(maxtime), L(maxtime),  L_wolb(maxtime), A(maxtime), A_ovipos(maxtime), A_wolb(maxtime), A_ovipos_wolb(maxtime), L_data(no_pdays);

// Declares two matrices; P and P-wolb. Tracks pupae at different stages of development.
// One row per day (i.e. 2500 total) - for each row, 3 columns: new pupae, day-old puape and a third, unused column
Matrix_Double P(maxtime, Row_Double(3)), P_wolb(maxtime, Row_Double(3));

// Declares vectors (nothing is used other than hdate...)
// e.g. hdate: hatch date for each of the 1000 cohorts
vector<int> pupae(no_pdays),pdate(no_pdays),hatch(no_cohorts),hdate(no_cohorts),larvae_sampled(no_cohorts), pupae_sampled(no_pdays);

// A tool for building and manipulating text strings; reads/writes to a string in memory
stringstream strstm;

int release_day;
double release_size, DI_cost, DI_L_cost;

//--------------------------------- MAIN function ---------------------------------------
// Coordinates the overall flow of the program: 
// Setting up input and output files, assigning cohort properties, running the model, and then writing all the results to the output files.

int main(){

// Declaring input and output files
ofstream cohort_means_out,cohort_means_wolb_out, cohort_stds_out, cohort_stds_wolb_out, mu_p_out, mu_p_wolb_out, A_out, A_wolb_out, L_out, L_wolb_out, freqA2_out, lambda_out;
ifstream rel_siz_in;

// Calls read_data function
read_data(cohort_means_out,cohort_means_wolb_out, cohort_stds_out, cohort_stds_wolb_out, mu_p_out, mu_p_wolb_out, L_out, L_wolb_out, A_out, A_wolb_out, freqA2_out, lambda_out, rel_siz_in);

// Assign maximum development times to each cohort - appears to be redundant
// i.e. the number of days remaining in the simulation from when that cohort hatches
// So a cohort that hatches late in the simulation has fewer days left and therefore a shorter maximum possible development time.
vector<int> max_dt(no_cohorts);
for (int i=0; i<no_cohorts; i++){
	max_dt.at(i) = maxtime - hdate.at(i);
}

// Calls Run_Model function
// It passes in the maximum development times and three output file streams that need to be written to during the simulation rather than after (after this line/function)
Run_Model(max_dt, mu_p_out, freqA2_out, lambda_out);

// Writing cohort development time statistics to file
// This loops through all 1000 cohorts and writes four statistics for each one to their respective output files - each value is separated by a space
// These are the OBSERVED statistics calculated from the full emergence record
for (int j=0; j<no_cohorts; j++) {
	cohort_means_out << mean_dt.at(j) << " "; // mean development time of uninfected larvae
	cohort_means_wolb_out << mean_dt_wolb.at(j) << " "; // mean development time of infected larvae
	cohort_stds_out << std_dt.at(j) << " "; // standard deviation of development time of uninfected larvae
	cohort_stds_wolb_out << std_dt_wolb.at(j) << " "; // standard deviation of development time of infected larvae
}
// New line (basically like typing "enter" on keyboard) 
cohort_means_out << endl; 
cohort_means_wolb_out << endl; 
cohort_stds_out << endl; 
cohort_stds_wolb_out << endl; 

// Writing daily pupal emergence to file
// This loops through all 2500 days and writes the total number of uninfected and infected pupae that eclosed on each day - each value is separated by a space
for (int j=0; j<maxtime; j++) {
	mu_p_out << mu_p.at(j) << " ";
	mu_p_wolb_out << mu_p_wolb.at(j) << " ";
}
mu_p_out << endl;
mu_p_wolb_out << endl;

// Writing daily larval and adult counts to file
// This loops through all 2500 days and writes four population counts for each day - each value is separated by a space
for (int j=0; j<maxtime; j++) {
	L_out << L.at(j) << " "; // Total uninfected larvae present
	L_wolb_out << L_wolb.at(j) << " "; // Total infected larvae present
	A_out << A.at(j) << " "; // Total uninfected adults present
	A_wolb_out << A_wolb.at(j) << " "; // Total infected adults present
}
L_out << endl;
A_out << endl;


return 0;

}

//---------------------------------function "read_data"---------------------------------------
// This function handles all the input and output file management for the model
// It reads in the hatch dates, opens all the output files where results will be saved, and reads in the key parameter values that control the simulation.

void read_data(ofstream& cohort_means_out, ofstream& cohort_means_wolb_out, ofstream& cohort_stds_out,ofstream& cohort_stds_wolb_out, ofstream& mu_p_out, ofstream& mu_p_wolb_out, ofstream& L_out, ofstream& L_wolb_out, ofstream& A_out, ofstream& A_wolb_out, ofstream& freqA2_out, ofstream& lambda_out, ifstream& rel_siz_in){

// Opening and reading the hatch date file of the 1000 cohorts
// Stores data in the hdate array
ifstream hdate_data, hatch_data;
hdate_data.open("hdate3.txt");
for (int i=0; i<no_cohorts; i++){
hdate_data >> hdate.at(i);
}

// Open and read in the names of all the output files from the initialisation file (gyp_sim_inits.txt)
string cohort_means_file,cohort_means_wolb_file, cohort_stds_file, cohort_stds_wolb_file, mu_p_file, mu_p_wolb_file, L_file, L_wolb_file, A_file, A_wolb_file, rel_siz_file, freqA2_file, freqL_file, lambda_file;
cin >> cohort_means_file >> cohort_means_wolb_file >> cohort_stds_file >> cohort_stds_wolb_file  >> mu_p_file >> L_file >> A_file >> mu_p_wolb_file >> L_wolb_file>> A_wolb_file >>  freqA2_file >>lambda_file >> rel_siz_file;
// Print file names to console
cout <<   " cohort_means_file " << cohort_means_file<< " cohort_means_wolb_file " << cohort_means_wolb_file << " cohort_stds_file " << cohort_stds_file << " cohort_stds_wolb_file " << cohort_stds_wolb_file <<  " mu_p_file " << mu_p_file << " L_file " << L_file << " A_file "<<  A_file  << " mu_p_wolb_file " << mu_p_wolb_file << " L_wolb_file " << L_wolb_file << " A_wolb_file " << A_wolb_file << " release_siz_file " << rel_siz_file << " freqA2_file " << freqA2_file << " lambda_file " << lambda_file << endl;

// Opening output files - maybe this can be looped/shortened
// 1: Clears the string stream (strstm.clear() and strstm.str("")) — resetting it so it doesn't carry over content from the previous file name
// 2: Loads the file name into the string stream (strstm << cohort_means_file)
// 3: Converts it to a string (strstm.str())
// 4: Opens the file for writing (cohort_means_out.open(filename4.c_str()))

// Mean development times of uninfected larvae per cohort
strstm.clear();
strstm.str("");
strstm <<cohort_means_file;
string filename4=strstm.str();
cohort_means_out.open(filename4.c_str());

// Standard deviations of development times of uninfected larvae per cohort
strstm.clear();
strstm.str("");
strstm <<cohort_stds_file;
string filename5=strstm.str();
cohort_stds_out.open(filename5.c_str());

// Number of uninfected pupae eclosing each day
strstm.clear();
strstm.str("");
strstm <<mu_p_file;
string filename10=strstm.str();
mu_p_out.open(filename10.c_str());

// Number of uninfected larvae present each day
strstm.clear();
strstm.str("");
strstm <<L_file;
string filename11=strstm.str();
L_out.open(filename11.c_str());

// Number of uninfected adults present each day
strstm.clear();
strstm.str("");
strstm <<A_file;
string filename17=strstm.str();
A_out.open(filename17.c_str());

// Number of infected pupae eclosing each day
strstm.clear();
strstm.str("");
strstm <<mu_p_wolb_file;
string filename18=strstm.str();
mu_p_wolb_out.open(filename18.c_str());

// Number of infected larvae present each day
strstm.clear();
strstm.str("");
strstm <<L_wolb_file;
string filename19=strstm.str();
L_wolb_out.open(filename19.c_str());

// Number of infected adults present each day
strstm.clear();
strstm.str("");
strstm <<A_wolb_file;
string filename20=strstm.str();
A_wolb_out.open(filename20.c_str());

// Mean development times of infected larvae per cohort
strstm.clear();
strstm.str("");
strstm <<cohort_means_wolb_file;
string filename21=strstm.str();
cohort_means_wolb_out.open(filename21.c_str());

// Standard deviations of development times of infected larvae per cohort
strstm.clear();
strstm.str("");
strstm <<cohort_stds_wolb_file;
string filename22=strstm.str();
cohort_stds_wolb_out.open(filename22.c_str());

// INPUT FILE: release size
strstm.clear();
strstm.str("");
strstm <<rel_siz_file;
string filename24=strstm.str();
rel_siz_in.open(filename24.c_str());

// Wolbachia frequency on the final day of release
strstm.clear();
strstm.str("");
strstm <<freqA2_file;
string filename25=strstm.str();
freqA2_out.open(filename25.c_str());

// Per-capita female fecundity for each cohort
strstm.clear();
strstm.str("");
strstm <<lambda_file;
string filename26=strstm.str();
lambda_out.open(filename26.c_str());

// Reading key parameter values (not files per se)
cin >> release_day; // the day on which the first Wolbachia release occurs (read from the initialisation file)
rel_siz_in >> release_size; // the number of infected mosquitoes released each week (read from the separate release size file)
cin >> DI_cost; // the additional density-independent daily mortality experienced by adults in the field
cin >> DI_L_cost; // the additional density-independent daily mortality experienced by larvae in the field

// Printing mortality parameter values to console
cout << "DI_cost " << DI_cost << " (additional density-INdependent daily mortality experienced by adults in the field environment) " << endl;
cout << "DI_L_cost " << DI_L_cost << " (additional density-INdependent daily mortality experienced by larvae in the field environment) " << endl;

}

//---------------------------------function "Run_Model"----------------------------------------
// This is where all variables are set up before we run the model

// This basically says "I am starting the Run_Model function" which starts from here "{" until the end of the file "}".
// I think the things in brackets here are the main outputs after running the model
double Run_Model(const vector<int> max_dt, ofstream& mu_p_out, ofstream& freqA2_out, ofstream& lambda_out){

// This section creates all the data containers (lists and tables) that will be used within the function
// Now we're creating the containers that we told it earlier we would create and what format they would have
// "non_emerg_prob" is the probability per cohort that a larva has not yet emerged as pupa (1 [100%] and decreases over time)
// "emerg_flag" a flag to track if enough pupae have emerged per cohort to start calculating stats
vector<double>  non_emerg_prob(no_cohorts), non_emerg_prob_wolb(no_cohorts), L_avg_cohort(no_cohorts), hatch_sim(no_cohorts), hatch_sim_wolb(no_cohorts), A_wolb_imm(maxtime), A_ovipos_wolb_imm(maxtime), mn_dt_gam(no_cohorts), mn_dt_wolb_gam(no_cohorts), std_dt_gam(no_cohorts), std_dt_wolb_gam(no_cohorts),no_emerg_tot(no_cohorts), no_emerg_tot_wolb(no_cohorts), L_avg(no_cohorts), lambda(no_cohorts), lambda_wolb(no_cohorts);
vector<int> emerg_flag(no_cohorts);
//Make an array to store the number of larvae of each age at each time step
// Large tables (2500 days x 1000 cohorts) storing how many larvae are alive each day
Matrix_Double L_cohort(maxtime, Row_Double(no_cohorts));
Matrix_Double L_cohort_wolb(maxtime, Row_Double(no_cohorts));
//Make an array to store the number of individuals from each cohort that emerge as pupae at each time 
// Similar to above, stores how many larvae emerged as pupae each day
Matrix_Double emerg_record(maxtime, Row_Double(no_cohorts)), emerg_record_wolb(maxtime, Row_Double(no_cohorts));
// These declare many individual variables — single numbers rather than lists (double and integer formats) — that will be used as temporary working values during the simulation. 
// e.g. surv_L is daily larval survival rate, and prob will hold the probability of a larva emerging on a given day.
double tau_p, surv_L, L_avg1, L_cum2, L_cum, denom2, Dt_shp, Dt_shp_wolb, Dt_scl, Dt_scl_wolb, prob, prob_wolb, no_emerge,no_emerge_wolb, H_cum, av_growth, L_av1, L_av2, L_avg_f, sh, w, freqA, lambda_max, survA_imm, intc_TL, intc_TL_wolb, alpha_TL, alpha_TL_wolb, exp_TL, exp_TL_wolb, a_value, a_value_wolb, b_value, b_value_wolb, b_value2, b_value2_wolb, intc_f, alpha_f, survA, survA_wolb, lambda1, lambda_wolb1, freq_end, lambda_min, mean_max,std_max,std_min, freq_end2, freqL_end,mean_dt1,mean_lambda,mean_std_dt;
int index,first_hatch_date, max_cohort, L_max, time_lag, count, tlagh, tlagg,tlagi,tlagp, min_cohort, max_dt_int, release_end, release_day_init, dens_lag;

// Time lags
first_hatch_date=hdate.at(0);
tlagh = 5; // lag between oviposition and hatching
tlagg = 6; // minimum time between emergence as adults and first oviposition (first?)
tlagi = 21; // time over which larval density is averaged (to calculate effect of density)
tlagp = 2; // time required for pupal development into adults

// Biological constants/parameters
// The _wolb versions of each are set to identical values, meaning Wolbachia infection is assumed not to affect these traits in this version of the model.
intc_f=28.0; // eggs laid per female per day; fecundity decreases with density (so it does incorporate carrying capacity?)
alpha_f=-3.3; // ""
intc_TL = 1.8; // how long larvae take to develop into pupae
intc_TL_wolb = 1.8; // ""
alpha_TL = 0.536; // ""
alpha_TL_wolb = 0.536; // ""
exp_TL = 0.533; // ""
exp_TL_wolb = 0.533; // ""
a_value = 0.22; // variability in development time
a_value_wolb = 0.22; // ""
b_value = 0.0168; // ""
b_value_wolb = 0.0168; // ""
b_value2 = 0.867; // ""
b_value2_wolb = 0.867; // ""

// Biological limits (upper and lower)
max_dt_int=100; //  cohorts older than 100 days are ignored to save computation time
survA = 1-0.03-DI_cost; // the daily survival rate of adults. Base mortality is 3% per day, plus any additional cost from Wolbachia (DI_cost)
lambda_max=14; // a female can lay at most 14 eggs per day
lambda_min=0.5; // a female lays at least 0.5 eggs per day even under very crowded conditions
mean_max = 60; // if crowding is extreme, development time is capped at 60 days
std_max = 40; 
std_min = 1.0;
L_max = 4200; // the larval density above which the caps on development time kick in


//Wolbachia-specific parameters
sh=0.99; // cytoplasmic incompatibility strength. When an uninfected female mates with a Wolbachia-infected male, 99% of her eggs fail to survive
w=0.01; // the rate at which Wolbachia-infected mothers accidentally produce uninfected offspring (1%)

release_day_init = release_day;
release_end = release_day + 89; // infected mosquitoes are released weekly for 89 days

survA_imm=survA; // survival rate (set above) is the same for infected and uninfected mosquitoes (why create survA_imm and survA_wolb?)
survA_wolb = survA;

// Initialising cohorts - setting everything to zero
// This long block simply sets all the tables and lists created earlier to their starting values before the simulation begins
// Larval, pupal and adult counts
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

//Initialise non_emerg_prob array - starts at 1 (100%) i.e. no larvae have emerged yet
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
// The whole code repeats itself every day from day 1 to 2500
// On each day, it works through every cohort (a batch of eggs laid on the same day).
for (int time1=1; time1<maxtime; time1++){

//-------------------LARVAE--------------------------------------------------
// In summary: each day, the code survives existing larvae, hatches any new cohorts (with a biologically informed hatch number based on adult density and fecundity), and totals everything up (total happens at the end of the pupae section).

// This calculates the daily survival rate of larvae. 
// It starts from a base of 95% survival, then subtracts any extra mortality specified for the field environment (DI_L_cost). 
// So if you set DI_L_cost = 0.1, then larvae have a 85% chance of surviving each day.
	surv_L = 0.95 - DI_L_cost;

	// These 4 lines track how many active cohorts the code needs to loop over (i.e. min and max cohorts)
	// This finds the most recent cohort that has already hatched by today (time1). 
	// It loops through all cohorts and keeps updating max_cohort to the last one whose hatch date is on or before today. 
	// Starting at -1 means "none yet". 
	max_cohort=-1;
	for (int i=0; i<no_cohorts; i++) if (hdate.at(i)<=time1) max_cohort = i;
	// This finds the oldest cohort still worth tracking — any cohort that hatched within the last 100 days (max_dt_int = 100). 
	// Cohorts older than that are assumed to have all developed into adults already and can be ignored.
	min_cohort=no_cohorts-1;
	for (int i=no_cohorts-1; i>=0; i--) if (hdate.at(i)>=time1 - max_dt_int) min_cohort = i;		

	// Now the code loops through every active cohort (from oldest still alive to newest). 
	// For each cohort, it takes the number of larvae alive yesterday and multiplies by the survival rate to get today's survivors — done separately for uninfected (L_cohort) and Wolbachia-infected (L_cohort_wolb) larvae.
	for (int cohort=min_cohort; cohort<=max_cohort; cohort++){

		L_cohort[time1][cohort] = L_cohort[time1-1][cohort] * surv_L;
		L_cohort_wolb[time1][cohort] = L_cohort_wolb[time1-1][cohort] * surv_L;

		// This checks: is today the hatch date for this cohort? 
		// If yes, the following block runs to set up the new cohort of freshly hatched larvae.
		if (time1 == hdate.at(cohort)){

// This block figures out the average larval density that the mothers of this cohort experienced when they were larvae themselves.
// This affects how many eggs they lay (fecundity).

// One week larval density lagged by 6 days
			L_avg.at(cohort) = 0; // average larval density experienced by the mothers of this cohort when they were larvae

			// dens_lag calculates the starting point of this look-back window, accounting for three delays:
			// tlagh = 5: days between egg-laying and hatching
			// tlagg = 6: days between emerging from pupa and first egg-laying
			// tlagi = 21: the window of days over which larval density is averaged (this is the larval development period - The idea is that a female's fecundity is shaped by the larval density she was exposed to throughout her entire larval life, not just on a single day. Averaging over 21 days is a way of capturing that cumulative experience.)
			// The code has to look back in time by a combined lag of 5 + 6 + 21 = 32 days, and average over a 21-day window
			dens_lag = time1-tlagh-tlagi-tlagg;

			// These 3 if/else branches handle edge cases 
			// Case 1: all data exists within simulation. If the start of the 21-day window falls on or after the first hatch date, there is real simulation data for the entire window. 
			// The code simply sums up the total larvae (infected + uninfected) across all 21 days and divides by 21 to get the average. This is the normal case.
			if (dens_lag>=hdate.at(1)){
			//lags: between oviposition and hatching, time over which larval density is averaged, time to become gravid after pupation
				for (int j=time1-tlagh-tlagi-tlagg; j<time1-tlagh-tlagg; j++) L_avg.at(cohort)+=L.at(j) + L_wolb.at(j);
					L_avg.at(cohort) = L_avg.at(cohort)/tlagi; 
			}
			// Case 2: the window partially overlaps with the simulation. This triggers when dens_lag is zero or positive, but still falls before the first hatch date — meaning the window exists within the simulation's time frame, but partially before any larvae existed.
			else if (dens_lag>=0) {
				// Case 2a: If even the end of the 21-day window is before the first hatch date, then there's no real data at all within the window. The entire average is set to the default value of 60.
				if (time1-tlagh-tlagg<=hdate.at(1)) L_avg.at(cohort)=60; // if the end of the lag is before the first hatch date
				// Case 2b: If the window is partly before and partly after the first hatch date then the code fills the early (pre-hatch) portion with the default value of 60 per day, and uses real simulation data for the later portion. 
				// The two halves are then combined and divided by 21 to get a blended average.
				if (time1-tlagh-tlagg>hdate.at(1)) {
					for (int j=dens_lag; j<=hdate.at(1); j++) L_avg.at(cohort)+=60;
					for (int j=hdate.at(1)+1; j<=time1-tlagh-tlagg; j++) L_avg.at(cohort)+=L.at(j) + L_wolb.at(j);
					L_avg.at(cohort) = L_avg.at(cohort)/tlagi;
					
				} 
			}
			// Case 3: The window starts before day zero. If dens_lag is negative — meaning the look-back window starts before the simulation even began — the code doesn't attempt any calculation at all. 
			// It simply assigns the default value of 60 as the average larval density.
			else L_avg.at(cohort) = 60; 

			// This calculates per-capita female fecundity (how many eggs each female lays per day on average). 
			// It uses a log-linear formula: fecundity decreases as larval density increases (density dependence). 
			// The result is then clamped between a minimum (0.5) and maximum (14) to keep it biologically realistic.
			L_avg_f = L_avg.at(cohort);

			lambda1 = intc_f + alpha_f*log(L_avg_f);

			if (lambda1<lambda_min) lambda1=lambda_min;
			if (lambda1>lambda_max) lambda1 = lambda_max;
			lambda_wolb1 = lambda1;

			lambda.at(cohort) = lambda1;
			lambda_wolb.at(cohort) = lambda_wolb1;

			// This calculates how many larvae hatch today from eggs laid 5 days ago (time1-tlagh). There are two totals:
			// Uninfected larvae (hatch_sim): come from uninfected mothers, plus a small fraction (w = 0.01, the "leakage" parameter) from Wolbachia-infected mothers who occasionally produce uninfected offspring.
			// Infected larvae (hatch_sim_wolb): come from Wolbachia-infected mothers (both naturally infected A_ovipos_wolb and recently released A_ovipos_wolb_imm), minus that same small leakage fraction.
 			if (time1>tlagh+tlagg){
				hatch_sim.at(cohort) = lambda1*A_ovipos.at(time1-tlagh)+ w*lambda_wolb1*A_ovipos_wolb.at(time1-tlagh)+ w*lambda_max*A_ovipos_wolb_imm.at(time1-tlagh);
				hatch_sim_wolb.at(cohort) = lambda_wolb1*(1-w)*A_ovipos_wolb.at(time1-tlagh) + lambda_max*(1-w)*A_ovipos_wolb_imm.at(time1-tlagh);
			}

			// The freshly calculated hatch numbers are stored as the starting larval population for this new cohort today. 
			// This overwrites the survival calculation from earlier, because today is hatch day — there were no larvae from this cohort yesterday.
			L_cohort[time1][cohort] = hatch_sim.at(cohort);	
			L_cohort_wolb[time1][cohort] = hatch_sim_wolb.at(cohort);

		}//end of if time1==hdate.at(cohort) loop	

//-------------------PUPAE------------------------------------------------------

// For each active cohort, this section works out how many larvae transition into pupae today. 
// It does this using a statistical distribution (a gamma distribution) to spread emergence across time — because not all larvae in a cohort develop at exactly the same rate. 
// The spread of that distribution depends on the larval density the cohort experienced.

// This block of code updates the model's estimate of the mean and spread of development times for each cohort, based on how larval density has actually played out so far.
// The first line ensures this code only runs if 1) at least 5 days have passed since this cohort hatched (larvae need a minimum time before any can pupate), and
// 2) emerg_flag is still 0, meaning fewer than 2 individuals from this cohort have emerged so far
		if (time1>hdate.at(cohort)+4 && emerg_flag.at(cohort)==0){ // if the cohort could have begun emergence and less than 1 individual have emerged so far
		// check that this ordering of the cohorts and times is okay
 			L_avg1 = 0;
			L_cum2 = 0;
			denom2 = 0;

			// Calculating average larval density experienced by this cohort (for each possible emergence day (etime) from day 5 onwards)
			// 1) Calculates the cumulative larval density experienced from hatch up to that emergence day (L_cum)
			// 2) Converts that to a daily average by dividing by the number of days lived (etime - hdate)
			// 3) Weights that average by how many individuals actually emerged on that day (emerg_record)
			// i.e. weighted average density — individuals that emerged later (having experienced more days of density) contribute more to the average, and it's weighted by how many actually emerged at each time. This accounts for the fact that faster-developing larvae experienced fewer days of crowding.
			for (int etime = hdate.at(cohort)+5; etime <=time1; etime++){ // all emergence times "etime" past to present
				L_cum=0;
				H_cum=0;
				count=cohort;
				for (int k=hdate.at(cohort); k<etime; k++) {
					L_cum += (L.at(k) + L_wolb.at(k)); // cumulative exposure experienced at day before emergence time
				}
				L_cum-=H_cum;
				L_cum2 += (emerg_record[etime][cohort] + emerg_record_wolb[etime][cohort])*L_cum/(etime-hdate.at(cohort));//accumulating daily average exposure * no. emerged
				denom2 += emerg_record[etime][cohort] + emerg_record_wolb[etime][cohort];//cumulative no. emerged
			}

			// If any larvae have already emerged, use the weighted average just calculated. 
			// If none have emerged yet (nothing to weight by), fall back to a simple unweighted average of cumulative density divided by days so far.
			if (denom2>0) L_avg1 = L_cum2/denom2;
			else L_avg1 = L_cum/(time1 - hdate.at(cohort));

			// Calculates the expected mean development time for this cohort using the above average larval density and a power-law formula (with parameters intc_TL = 1.8, alpha_TL = 0.536, exp_TL = 0.533). 
			// Higher density → longer development time. The result is floored at zero to avoid nonsensical negative values. 
			// Infected and uninfected larvae get the same mean development time here.
			mn_dt_gam.at(cohort) = intc_TL + alpha_TL * pow(L_avg1, exp_TL);
			mn_dt_wolb_gam.at(cohort) = mn_dt_gam.at(cohort);

			if (mn_dt_gam.at(cohort)<0) mn_dt_gam.at(cohort)=0;
			if (mn_dt_wolb_gam.at(cohort)<0) mn_dt_wolb_gam.at(cohort)=0;

			// Setting the emerg_flag: once more than 1 individual has emerged from this cohort, the flag is set to 1. 
			// This means the density estimate is now considered stable enough to lock in — the code won't keep updating it on future days. This is why the whole block above was wrapped in emerg_flag==0.
			// Stopping the mean and standard deviation from updating once emergence starts helps avoid feedback loops so the shape of the gamma distribution doesn't keep changing once it's already underway. It basically says "we've seen enough, commit to this estimate"
			// I would move this to the end of this loop i.e. after the checks for mean & std dev
			if (denom2>1) {
				emerg_flag.at(cohort)=1;
				L_avg_cohort.at(cohort) = L_avg1;
			}

			// Calculating the standard deviation (spread) of development time - because not all larvae within one cohort emerge on the same day. Higher density -> more variable development/higher spread
			std_dt_gam.at(cohort) = a_value + b_value * pow(L_avg1, b_value2);
			std_dt_wolb_gam.at(cohort) = a_value_wolb + b_value_wolb * pow(L_avg1, b_value2_wolb);

			// Caps/checks for the mean & standard deviation of development time
			if (L_avg1>L_max) std_dt_gam.at(cohort) = std_max; // if density exceeded (L_max = 4200) then mean and std dev set to max allowed values
			if (L_avg1>L_max) std_dt_wolb_gam.at(cohort) = std_max; // same as above but for infected
			if (std_dt_gam.at(cohort) <std_min) std_dt_gam.at(cohort) = std_min; // minimum values - prevents collapsing to zero
			if (std_dt_wolb_gam.at(cohort) <std_min) std_dt_wolb_gam.at(cohort) = std_min;
			if (L_avg1>L_max) {mn_dt_gam.at(cohort) = mean_max; mn_dt_wolb_gam.at(cohort) = mean_max;std_dt_gam.at(cohort) = std_max;std_dt_wolb_gam.at(cohort) = std_max;}

		}//end of if (time1>hdate.at(cohort)+4 && emerg_flag.at(cohort)==0)
		
		// Converting mean and standard deviation to scale and shape parameters for gamma distribution
		if (time1>hdate.at(cohort)+4){
			Dt_shp = pow(mn_dt_gam.at(cohort)/std_dt_gam.at(cohort),2); // convert std dev to shape
			Dt_scl = mn_dt_gam.at(cohort)/Dt_shp; // convert mean to scale			
			Dt_shp_wolb = pow(mn_dt_wolb_gam.at(cohort)/std_dt_wolb_gam.at(cohort),2); // convert std dev to shape (infected)
			Dt_scl_wolb = mn_dt_wolb_gam.at(cohort)/Dt_shp_wolb; // convert mean to scale (infected)			
			time_lag = time1 - hdate.at(cohort);

			// Calculating probability of emerging today (time_lag = days since this cohort hatched) based on gamma distribution
			// It does this by taking the probability of emerging by today and subtracting the probability of having emerged by yesterday
			// The -5 offset reflects the minimum development time of 5 days
			//if (Dt_shp>0 && Dt_scl>0) 
			if (time_lag>5 && mn_dt_gam.at(cohort)>0 && Dt_shp>0 && Dt_scl>0 && Dt_shp<200 && Dt_scl<200){				 
				prob = gsl_cdf_gamma_P(time_lag-5,Dt_shp,Dt_scl) - gsl_cdf_gamma_P(time_lag-5-1,Dt_shp,Dt_scl); // uninfected
			}
			else prob=0;
			if (time_lag>5 && mn_dt_wolb_gam.at(cohort)>0 && Dt_shp_wolb>0 && Dt_scl_wolb>0 && Dt_shp_wolb<200 && Dt_scl_wolb<200){
				prob_wolb = gsl_cdf_gamma_P(time_lag-5,Dt_shp_wolb,Dt_scl_wolb) - gsl_cdf_gamma_P(time_lag-5-1,Dt_shp_wolb,Dt_scl_wolb); // infected
			}

			else prob_wolb=0;

			// Safety checks
			// If scale and shape exceed max of 200 or are zero, probability is set to zero to avoid mathematical errors
			if (mn_dt_gam.at(cohort)==0){ // uninfected
				Dt_shp = 9.0;
				Dt_scl = 0.2;
				prob = gsl_cdf_gamma_P(time_lag-4,Dt_shp,Dt_scl) - gsl_cdf_gamma_P(time_lag-4-1,Dt_shp,Dt_scl);
			   
			}
			if (mn_dt_wolb_gam.at(cohort)==0){ // infected
			 	Dt_shp_wolb = 9.0;
			 	Dt_scl_wolb = 0.2;
				prob_wolb = gsl_cdf_gamma_P(time_lag-4,Dt_shp_wolb,Dt_scl_wolb) - gsl_cdf_gamma_P(time_lag-4-1,Dt_shp_wolb,Dt_scl_wolb);
			}

			// Converting probability to actual numbers and updating larval counts
			// The number emerging today is calculated as: yesterday's larval count × the probability of emerging today, divided by non_emerg_prob 
			// non_emerg_prob is a running tracker of the probability of not yet having emerged, used to correctly condition the probability on larvae that are still present. 
			no_emerge = L_cohort[time1-1][cohort] * prob/non_emerg_prob.at(cohort);
			no_emerge_wolb = L_cohort_wolb[time1-1][cohort] * prob_wolb/non_emerg_prob_wolb.at(cohort);

			// Sanity checks
			// Ensure no_emerge doesn't exceed the number of larvae actually present or go negative.
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

			// Some more updates based on above calculations
			non_emerg_prob.at(cohort) *= (1-prob/non_emerg_prob.at(cohort));
			non_emerg_prob_wolb.at(cohort) *= (1-prob_wolb/non_emerg_prob_wolb.at(cohort));
			L_cohort[time1][cohort] -= no_emerge; // emerging larvae are removed from the larval pool
			L_cohort_wolb[time1][cohort] -= no_emerge_wolb;
			emerg_record[time1][cohort]=no_emerge; // emerging larvae number is recorded in emerg_record for future reference
			emerg_record_wolb[time1][cohort]=no_emerge_wolb;
			mean_dt.at(cohort)+=(time1-hdate.at(cohort))*no_emerge;
			mean_dt_wolb.at(cohort)+=(time1-hdate.at(cohort))*no_emerge_wolb;
			no_emerg_tot.at(cohort)+=no_emerge;
			no_emerg_tot_wolb.at(cohort)+=no_emerge_wolb;
			P[time1][0] += no_emerge; // emerging larvae are added to the pupal pool (P[time1][0])
			P_wolb[time1][0] += no_emerge_wolb;
		}//end of if (time1>hdate.at(cohort)+4) loop

		L.at(time1)+=L_cohort[time1][cohort];
		L_wolb.at(time1)+=L_cohort_wolb[time1][cohort];


	}//end of cohort loop


//-------------------ADULTS ----------------------------------------------------
// Surviving adults and pupae advancing
// This section updates the adult population each day, accounting for survival, pupae emerging into adults, and Wolbachia-infected mosquito releases.
// It also calculates how many adults are ready to lay eggs.
// Note pupal stage is only 2 days long.

	// Yesterday's newly emerged pupae (P[time1-1][0]) survive and move into their second pupal day (P[time1][1]).
	P[time1][1] = P[time1-1][0]*survA;
	P_wolb[time1][1] = P_wolb[time1-1][0]*survA_wolb; 

	A[time1] = A[time1-1]*survA; // Yesterday's adults are multiplied by the daily adult survival rate to give today's surviving adults (uninfected).
	A[time1] += P[time1-1][1]*survA; // Yesterday's second-day pupae (P[time1-1][1]) survive and graduate into adults, being added to today's adult count (uninfected).

	A_wolb[time1] = A_wolb[time1-1]*survA_wolb; // Same as above but for infected.
	A_wolb[time1] += P_wolb[time1-1][1]*survA_wolb; // Same as above but for infected.


	
    // Weekly releases of Wolbachia-infected mosquitoes into the population.

	// Whatever infected immigrants were alive yesterday, a proportion of them survive to today.
	// This line runs every day after the first release day, not only on release days.
	if (time1>release_day_init) A_wolb_imm[time1] = A_wolb_imm[time1-1]*survA_imm;

	// On each release day, release_size infected adults are added to A_wolb_imm — a separate tracker for recently released immigrants, distinct from the naturally produced infected adults. 
	if (time1==release_day && time1<release_end){ // release schedule
		A_wolb_imm[time1]+=release_size;
		release_day+=7; // After each release, release_day is advanced by 7 so the next release happens one week later. Releases stop once release_end is reached (90 days after the first release).
	}

	// Calculates how many adults are ready to lay eggs today, done separately for uninfected, naturally infected, and released infected adults.
	// Only females lay eggs, so counts are multiplied by 0.5
	// There is a 6-day lag (tlagg) between emerging from the pupa and being ready to lay eggs
	// Cytoplasmic incompatibility is applied to uninfected females via (1 - sh*freqA) — when the Wolbachia frequency (freqA) is high, more uninfected females mate with infected males and produce inviable eggs, reducing their effective egg-laying. sh = 0.99 means near-complete incompatibility
	// The adults from tlagg days ago are survival-adjusted forward to today using pow(survA, tlagg-tlagp-1), accounting for the fact that some will have died in the intervening days
	if (time1>=tlagg) {
		if (A[time1-tlagg]>0) freqA = (A_wolb[time1-tlagg] + A_wolb_imm[time1-tlagg])/(A[time1-tlagg] + A_wolb[time1-tlagg] + A_wolb_imm[time1-tlagg]);
		else freqA=0;

		A_ovipos[time1] = 0.5*(1-sh*freqA)*A[time1-tlagg+tlagp] * pow(survA,tlagg-tlagp-1); // uninfected
		A_ovipos_wolb[time1] = 0.5*A_wolb[time1-tlagg+tlagp] * pow(survA_wolb,tlagg-tlagp-1); // infected
		A_ovipos_wolb_imm[time1] = 0.5*A_wolb_imm[time1-tlagg+tlagp] * pow(survA_imm, tlagg-tlagp-1); // released - I thought only males were released
	}


	// Summary: Each day the adults section: carries yesterday's adults forward with survival, graduates 2-day-old pupae into adults, adds any scheduled Wolbachia releases, 
	// and calculates how many females of each type are ready to lay eggs today — factoring in maturation time, sex ratio, and cytoplasmic incompatibility.

}//end of main time loop

//------------------------- OUTPUTS----------------------------------------------------------


// Calculating the actual cohort mean and standard deviation of development time for each cohort, based on the full emergence record.
// These are the observed development time statistics from the simulation, as opposed to the predicted ones calculated during the time loop from the density formula.
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


// Summing daily pupal emergence across cohorts
// During the time loop, emerg_record tracked how many larvae from each cohort emerged as pupae on each day. 
// This section simply sums across all cohorts to give the total number of uninfected and infected pupae emerging on each day of the simulation — stored in mu_p and mu_p_wolb.
for (int i=0; i<maxtime; i++) {
	for (int j=0; j<no_cohorts; j++) {
		mu_p.at(i) += emerg_record[i][j];
		mu_p_wolb.at(i) += emerg_record_wolb[i][j];
	}
}

// Writing the per-capita female fecundity for each cohort to the output file
// Each cohort's fecundity was calculated during the time loop based on the larval density experienced by that cohort's mothers — this line just saves those values
for (int i=0; i<no_cohorts; i++) lambda_out << lambda.at(i) << endl;

// Calculating Wolbachia frequency at end of release
// Released immigrant infected adults (A_wolb_imm) are merged into the main infected adult count (A_wolb) — up to this point they had been tracked separately
for (int i=0; i<maxtime; i++) A_wolb.at(i) += A_wolb_imm.at(i);

// Wolbachia frequency on the final day of the release period - proportion of infected adults out of all adults
freq_end2 = A_wolb.at(release_end)/(A_wolb.at(release_end) + A.at(release_end));
freqA2_out << freq_end2 << endl; // saved to output file

freqL_end = L_wolb.at(release_end)/(L_wolb.at(release_end) + L.at(release_end));

// Calculating and printing summary statistics - mean fecundity and development times - averaged across cohorts 300 to 800
// This time window is chosen to represent the stable middle period of the simulation, avoiding the initialisation period at the start and the release period at the end
mean_dt1=0;mean_lambda=0;mean_std_dt=0;
for (int i=300; i<800; i++){
	mean_dt1+=mean_dt.at(i);
	mean_std_dt+=std_dt.at(i);
	mean_lambda += lambda.at(i);
} 
cout << endl;
cout << "mean_dt " << mean_dt1/500 <<" (the average (over time) of the mean development time of (uninfected) larvae) " << endl; // mean, uninfected larvae
cout << "mean_lambda " << mean_lambda/500 << " (the average (over time) of per-capita female fecundity) " << endl; // female fecunity
cout << "mean_std_dt " << mean_std_dt/500 << " (the average (over time) of the standard deviation of the development time of (uninfected) larvae " <<  endl; // standard deviation

return 0;			

// Summary: Once the main time loop finishes, this section tidies everything up 
// computing the true observed development time statistics from the full emergence record, summing pupal emergence across cohorts, 
// merging the immigrant and resident infected populations, calculating the final Wolbachia invasion frequency, and printing summary statistics to the console.

}




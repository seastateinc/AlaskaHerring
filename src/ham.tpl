// -------------------------------------------------------------------------- //
//               Age-structured model for Alaska herring stocks               //
//                                                                            //
//                               VERSION 0.1                                  //
//                                Jan  2015                                   //
//                                                                            //
//                                 AUTHORS                                    //
//                              Sherri Dressel                                //                                 
//                        sherri.dressel@alaska.gov                           //
//                               Sara Miller                                  //
//                          sara.miller@alaska.gov                            //
//                               Kray Van Kirk                                //
//                          kray.vankirk@alaska.gov                           //
//                                                                            //
//                   Built on code developed by Peter Hulson                  //
//                            pete.hulson@noaa.gov                            //
//                                                                            //
//                           Layout and references                            //
//                              Steven Martell                                //
//                          martell.steve@gmail.com                           //
//                                                                            //
// CONVENTIONS: Formatting conventions are based on                           //
//              The Elements of C++ Style (Misfeldt et al. 2004)              //
//                                                                            //
//                                                                            //
// NAMING CONVENTIONS:                                                        //
//             Macros       -> UPPERCASE                                      //
//             Constants    -> UpperCamelCase                                 //
//             Functions    -> lowerCamelCase                                 //
//             Variables    -> lowercase                                      //
//                                                                            //
//                                                                            //
//                                                                            //
// -------------------------------------------------------------------------- //
//-- CHANGE LOG:                                                             -//
//--  Jan 2015 - revision of legacy code::                                   -//
//--               :variable naming conventions                              -//
//--               :intra-annual calendar                                    -//
//--               :standardization of units across stocks                   -//
//--               :modification for potential code distribution             -//
//--  May 2015 - added code and scripts to a git-hub repo.                   -//
//--                                                                         -//
// -------------------------------------------------------------------------- //

DATA_SECTION

// |---------------------------------------------------------------------------|
// | CHECK FOR OPTIONAL COMMAND LINE ARGUMENTS & SET FLAGS
// |---------------------------------------------------------------------------|
// | b_simulation_flag	-> flag for running in simulation mode
// | rseed			-> random number seed for simulation 
	int b_simulation_flag;
	int rseed;
	LOCAL_CALCS
		int on = 0;
		b_simulation_flag = 0;
		if (ad_comm::argc > 1)
		{
			int on = 0;
			rseed  = 0;
			if ( (on=option_match(ad_comm::argc,ad_comm::argv,"-sim")) > -1 )
			{
				b_simulation_flag = 1;
				rseed = atoi(ad_comm::argv[on+1]);
			}
		}
	END_CALCS
// |---------------------------------------------------------------------------|


// |---------------------------------------------------------------------------|
// | STRINGS FOR INPUT FILES                                                   |
// |---------------------------------------------------------------------------|
// |- The files below are listed in the 'model.dat' file;
// |   nothing else should be in the 'model.dat' file
// |- These files should be named by the stock they are modeling;
// |   example: "sitka.dat", "seymour.dat"
// |- DO NOT use two word names such as "seymour cove.dat"
// |
// | DEBUG_FLAG             : Boolean Flag used for manual debugging
// | DataFile               : data to condition the assessment model    
// | ControlFile            : controls for years, phases, and block options 
	init_int DEBUG_FLAG;
	init_adstring DataFile;      
  init_adstring ControlFile;


// |---------------------------------------------------------------------------|
// | READ CONTENTS OF DATA FILE
// |---------------------------------------------------------------------------|
	!! ad_comm::change_datafile_name(DataFile);

// |---------------------------------------------------------------------------|
// | Model dimensions.
// |---------------------------------------------------------------------------|
// | dat_syr			-> first year of data
// | day_nyr			-> last year of data
// | mod_syr 			-> first year of model
// | mod_nyr			-> last year of model
// | sage					-> first age class
// | nage 				-> plus group age class
	init_int dat_syr;
	init_int dat_nyr;
	init_int mod_syr;
	init_int mod_nyr;
	init_int sage;
	init_int nage;
	vector age(sage,nage);
	!! age.fill_seqadd(sage,1);

// |---------------------------------------------------------------------------|
// | Time series data.  Catch in short tons. Comps in proportions.
// |---------------------------------------------------------------------------|
// | data_catch		-> Catch: colnames(Year,Catch,log.se)
// | data_sp_waa 	-> Spawn Weight-at-age (Year,weight-at-age)
// | data_cm_waa	-> Commercial catch weight-at-age (Year, weight-at-age)
// | data_cm_comp	-> Commercial catch composition (Year, age proporitons)
// | data_sp_comp -> Spawn Sample age composition (Year, age proportions)
// | data_egg_dep -> Egg deposition survey (Year, Index, log.se)
	init_matrix   data_catch(dat_syr,dat_nyr,1,3);
	init_matrix  data_sp_waa(dat_syr,dat_nyr,sage-1,nage);
	init_matrix  data_cm_waa(dat_syr,dat_nyr,sage-1,nage);
	init_matrix data_cm_comp(dat_syr,dat_nyr,sage-1,nage);
	init_matrix data_sp_comp(dat_syr,dat_nyr,sage-1,nage);
	init_matrix data_egg_dep(dat_syr,dat_nyr,1,3);
	init_matrix data_mileday(dat_syr,dat_nyr,1,3);

// |---------------------------------------------------------------------------|
// | END OF DATA FILE
// |---------------------------------------------------------------------------|
	init_int eof; 
	!! if(eof != 999){cout<<"Error reading data file, aborting."<<endl; exit(1);}



// |---------------------------------------------------------------------------|
// | READ CONTENTS OF CONTROL FILE
// |---------------------------------------------------------------------------|
	!! ad_comm::change_datafile_name(ControlFile);


// |-------------------------------------------------------------------------|
// | DESIGN MATRIX FOR PARAMETER CONTROLS                                    |
// |-------------------------------------------------------------------------|
// | - theta_DM -> theta is a vector of estimated parameters.
  int n_theta;
  !! n_theta = 5;
  init_matrix theta_DM(1,n_theta,1,4);
  vector    theta_ival(1,n_theta);
  vector      theta_lb(1,n_theta);
  vector      theta_ub(1,n_theta);
  ivector    theta_phz(1,n_theta);
	!! theta_lb  = column(theta_DM,2);
	!! theta_ub  = column(theta_DM,3);
	!! theta_phz = ivector(column(theta_DM,4));
	

// |---------------------------------------------------------------------------|
// | Controls for time-varying maturity
// |---------------------------------------------------------------------------|
	init_int nMatBlocks;
	init_ivector nMatBlockYear(1,nMatBlocks);

// |---------------------------------------------------------------------------|
// | Controls for time-varying maturity
// |---------------------------------------------------------------------------|
	init_int nMortBlocks;
	init_matrix cntrl_NaturalMortality(1,nMortBlocks,1,5);
	ivector             nMortBlockYear(1,nMortBlocks);
	ivector                   mort_phz(1,nMortBlocks);
	vector                     mort_lb(1,nMortBlocks);
	vector                     mort_ub(1,nMortBlocks);
	vector                   mort_ival(1,nMortBlocks);
	!! nMortBlockYear = ivector(column(cntrl_NaturalMortality,1));
	!! mort_ival      = column(cntrl_NaturalMortality,2);
	!! mort_lb        = column(cntrl_NaturalMortality,3);
	!! mort_ub        = column(cntrl_NaturalMortality,4);
	!! mort_phz       = ivector(column(cntrl_NaturalMortality,5));
	
	
	
INITIALIZATION_SECTION
  theta theta_ival;
  log_m mort_ival;
  
  


PARAMETER_SECTION

// |---------------------------------------------------------------------------|
// | POPULATION PARAMETERS
// |---------------------------------------------------------------------------|
// | - theta(1) -> log natural mortality
// | - theta(2) -> log initial average age-3 recruitment for ages 4-9+ in dat_styr
// | - theta(3) -> log average age-3 recruitment from dat_styr to dat_endyr
// | - theta(4) -> log of unfished recruitment.
// | - theta(5) -> log of recruitment compensation (reck > 1.0)
  init_bounded_number_vector theta(1,n_theta,theta_lb,theta_ub,theta_phz);
  number log_natural_mortality;
  number log_rinit;
  number log_rbar;
  number log_ro;
  number log_reck;

// |---------------------------------------------------------------------------|
// | MATURITY PARAMETERS
// |---------------------------------------------------------------------------|
// | - mat_params[1] -> Age at 50% maturity
// | - mat_params[2] -> Slope at 50% maturity
	init_bounded_vector_vector mat_params(1,nMatBlocks,1,2,0,10,2);
	matrix mat(mod_syr,mod_nyr,sage,nage);

// |---------------------------------------------------------------------------|
// | NATURAL MORTALITY PARAMETERS
// |---------------------------------------------------------------------------|
// | - mat_params[1] -> Age at 50% maturity
// | - mat_params[2] -> Slope at 50% maturity
	init_bounded_number_vector log_m(1,nMortBlocks,mort_lb,mort_ub,mort_phz);
	matrix Mij(mod_syr,mod_nyr,sage,nage);



	objective_function_value f;

PROCEDURE_SECTION

// |---------------------------------------------------------------------------|
// | RUN STOCK ASSEAAMENT MODEL ROUTINES
// |---------------------------------------------------------------------------|
// | PSUEDOCODE:
// | - initialize model parameters.
// | - initialize Maturity Schedule information.
// | - get natural mortality schedules.
// |    - get fisheries selectivity schedules.
// | - initialize State variables
// | - update State variables
// |---------------------------------------------------------------------------|
  initializeModelParameters();
  if(DEBUG_FLAG) cout<<"--> Ok after initializeModelParameters <--"<<endl;

	initializeMaturitySchedules();
  if(DEBUG_FLAG) cout<<"--> Ok after initializeAgeSchedules    <--"<<endl;

  getNaturalMortality();
  if(DEBUG_FLAG) cout<<"--> Ok after getNaturalMortality       <--"<<endl;
  

// |---------------------------------------------------------------------------|



	f = sum(square(theta));

FUNCTION void initializeModelParameters()
  log_natural_mortality = theta(1);
  log_rinit             = theta(2);
  log_rbar              = theta(3);
  log_ro                = theta(4);
  log_reck              = theta(5);

FUNCTION void initializeMaturitySchedules() 
	int iyr = mod_syr;
	mat.initialize();
	for(int h = 1; h <= nMatBlocks; h++) {
		dvariable mat_a = mat_params(h,1);
		dvariable mat_b = mat_params(h,2);

		// fill maturity array using logistic function
		do{
			mat(iyr++) = plogis(age,mat_a,mat_b);
		} while(iyr <= nMatBlockYear(h));	
	}
	//COUT(mat);
	
FUNCTION void getNaturalMortality()
	int iyr = mod_syr;
	Mij.initialize();
	for(int h = 1; h <= nMortBlocks; h++){
		dvariable mi = exp(log_m(h));

		// fill mortality array by block
		do{
			Mij(iyr++) = mi;
		} while(iyr <= nMortBlockYear(h));
	}	
	//COUT(Mij);


GLOBALS_SECTION
	#include <admodel.h>
	#include <string.h>
	#include <time.h>



  #undef REPORT
  #define REPORT(object) report << #object "\n" << setw(8) \
  << setprecision(4) << setfixed() << object << endl;

  #undef COUT
  #define COUT(object) cout << #object "\n" << setw(6) \
  << setprecision(3) << setfixed() << object << endl;

  template<typename T>
  dvar_vector plogis(const dvector x, T location, T scale)
  {
  	return(1.0 / (1.0 + exp(-(x-location)/scale)));
  }


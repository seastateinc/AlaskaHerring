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
	init_int dat_eof; 
	!! if(dat_eof != 999){cout<<"Error reading data file, aborting."<<endl; exit(1);}



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
  !! theta_ival = column(theta_DM,1);
	!! theta_lb  = column(theta_DM,2);
	!! theta_ub  = column(theta_DM,3);
	!! theta_phz = ivector(column(theta_DM,4));
	

// |---------------------------------------------------------------------------|
// | Controls for time-varying maturity
// |---------------------------------------------------------------------------|
	init_int mat_phz;
	init_int nMatBlocks;
	init_ivector nMatBlockYear(1,nMatBlocks);

// |---------------------------------------------------------------------------|
// | Controls for natural mortality rate deviations in each block.
// |---------------------------------------------------------------------------|
	init_int mort_dev_phz;
	init_int nMortBlocks;
	init_ivector nMortBlockYear(1,nMortBlocks);
	

// |---------------------------------------------------------------------------|
// | Controls for selectivity parameters
// |---------------------------------------------------------------------------|
// | - nSlxCols 		» number of columns in selectivity design matrix
// | - nSlxBlks 			» number of selectivity blocks/patterns 
// | - selex_cont			» matrix of controls to be read in from control file.
	int nSlxCols;
	!! nSlxCols = 9;
	init_int nSlxBlks;
	init_matrix selex_cont(1,nSlxBlks,1,nSlxCols);
	ivector       nSelType(1,nSlxBlks);
	ivector       nslx_phz(1,nSlxBlks);
	ivector      nslx_rows(1,nSlxBlks);
	ivector      nslx_cols(1,nSlxBlks);
	ivector      	nslx_syr(1,nSlxBlks);
	ivector      	nslx_nyr(1,nSlxBlks);


	LOCAL_CALCS
		nSelType = ivector(column(selex_cont,2));
		nslx_phz = ivector(column(selex_cont,7));
		nslx_syr = ivector(column(selex_cont,8));
		nslx_nyr = ivector(column(selex_cont,9));

		// determine dimensions for log_slx_pars ragged object.
		for(int h = 1; h <= nSlxBlks; h++){
			nslx_rows(h) = 1;
			switch(nSelType(h)){
				case 1: // logistic 2-parameters
					nslx_cols = int(2);
				break;
			}
		}
	END_CALCS

// |---------------------------------------------------------------------------|
// | END OF Control FILE
// |---------------------------------------------------------------------------|
	init_int ctl_eof;
	LOCAL_CALCS
		if(ctl_eof != 999){
			cout<<"Error reading control file, aborting."<<ctl_eof<<endl; 
			exit(1);
		}
	END_CALCS

INITIALIZATION_SECTION
	theta theta_ival;
  


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
  init_bounded_dev_vector log_rinit_devs(sage+1,nage,-5.0,5.0,2);
  init_bounded_dev_vector log_rbar_devs(mod_syr,mod_nyr,-5.0,5.0,2);


// |---------------------------------------------------------------------------|
// | MATURITY PARAMETERS
// |---------------------------------------------------------------------------|
// | - mat_params[1] -> Age at 50% maturity
// | - mat_params[2] -> Slope at 50% maturity
	init_bounded_matrix mat_params(1,nMatBlocks,1,2,0,10,mat_phz);
	matrix mat(mod_syr,mod_nyr,sage,nage);

// |---------------------------------------------------------------------------|
// | NATURAL MORTALITY PARAMETERS
// |---------------------------------------------------------------------------|
// | - log_m_dev 		-> deviations in natural mortality for each block.
// | - Mij					-> Array for natural mortality rate by year and age.
	init_bounded_dev_vector log_m_dev(1,nMortBlocks,-5.0,5.0,mort_dev_phz);
	matrix Mij(mod_syr,mod_nyr,sage,nage);

// |---------------------------------------------------------------------------|
// | SELECTIVITY PARAMETERS
// |---------------------------------------------------------------------------|
// | - log_slx_pars	» parameters for selectivity models (ragged object).
	init_bounded_matrix_vector log_slx_pars(1,nSlxBlks,1,nslx_rows,1,nslx_cols,-25,25,nslx_phz);
	matrix log_slx(mod_syr,mod_nyr,sage,nage);
	LOCAL_CALCS
		if( ! global_parfile ){
			for(int h = 1; h <= nSlxBlks; h++){
				switch(nSelType(h)){
					case 1: //logistic
						log_slx_pars(h,1,1) = log(selex_cont(h,3));
						log_slx_pars(h,1,2) = log(selex_cont(h,4));
					break; 
				}
			}	
		}

	END_CALCS


// |---------------------------------------------------------------------------|
// | VECTORS
// |---------------------------------------------------------------------------|
// | - ssb 			» spawning stock biomass at the time of spawning.
	vector ssb(mod_syr,mod_nyr);	

// |---------------------------------------------------------------------------|
// | MATRIXES
// |---------------------------------------------------------------------------|
// | - Nij 			» numbers-at-age N(syr,nyr,sage,nage)
// | - Nij 			» mature numbers-at-age O(syr,nyr,sage,nage)
// | - Pij 			» numbers-at-age P(syr,nyr,sage,nage) post harvest.
// | - Sij 			» selectivity-at-age 
// | - Qij 			» vulnerable proportions-at-age
// | - Cij    	» predicted catch-at-age in numbers.
	matrix Nij(mod_syr,mod_nyr+1,sage,nage);
	matrix Oij(mod_syr,mod_nyr+1,sage,nage);
	matrix Pij(mod_syr,mod_nyr+1,sage,nage);
	matrix Sij(mod_syr,mod_nyr+1,sage,nage);
	matrix Qij(mod_syr,mod_nyr+1,sage,nage);
	matrix Cij(mod_syr,mod_nyr+1,sage,nage);

	objective_function_value f;



PROCEDURE_SECTION
	
// |---------------------------------------------------------------------------|
// | RUN STOCK ASSEAAMENT MODEL ROUTINES
// |---------------------------------------------------------------------------|
// | PSUEDOCODE:
// | - initialize model parameters.
// | - initialize Maturity Schedule information.
// | - get natural mortality schedules.
// | - get fisheries selectivity schedules.
// | - initialize State variables
// | - update State variables
// | 		- calculate spawning stock biomass
// | 		- calculate age-composition residuals
// |---------------------------------------------------------------------------|
  
  initializeModelParameters();
  if(DEBUG_FLAG) cout<<"--> Ok after initializeModelParameters      <--"<<endl;

	initializeMaturitySchedules();
  if(DEBUG_FLAG) cout<<"--> Ok after initializeMaturitySchedules    <--"<<endl;

  calcNaturalMortality();
  if(DEBUG_FLAG) cout<<"--> Ok after calcNaturalMortality           <--"<<endl;
  
  calcSelectivity();
  if(DEBUG_FLAG) cout<<"--> Ok after calcSelectivity                <--"<<endl;
  
  initializeStateVariables();
  if(DEBUG_FLAG) cout<<"--> Ok after initializeStateVariables       <--"<<endl;
  
  updateStateVariables();
  if(DEBUG_FLAG) cout<<"--> Ok after updateStateVariables           <--"<<endl;
	
	calcSpawningStockRecruitment();
  if(DEBUG_FLAG) cout<<"--> Ok after calcSpawningStockRecruitment   <--"<<endl;

  calcAgeCompResiduals();
  if(DEBUG_FLAG) cout<<"--> Ok after calcAgeCompResiduals           <--"<<endl;

  exit(1);
// |---------------------------------------------------------------------------|



	f = sum(square(theta));

FUNCTION void initializeModelParameters()

  log_natural_mortality = theta(1);
  log_rinit             = theta(2);
  log_rbar              = theta(3);
  log_ro                = theta(4);
  log_reck              = theta(5);
  //COUT(log_natural_mortality);


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
	

FUNCTION void calcNaturalMortality()
	
	int iyr = mod_syr;
	Mij.initialize();
	for(int h = 1; h <= nMortBlocks; h++){
		dvariable mi = exp(log_natural_mortality + log_m_dev(h));

		// fill mortality array by block
		do{
			//cout<<iyr<<"\t"<<theta<<endl;
			Mij(iyr++) = mi;
		} while(iyr <= nMortBlockYear(h));
	}		


FUNCTION void calcSelectivity()
	/**
		- Loop over each of the selectivity block/pattern
		- Determine which selectivity type is being used.
		- get parameters from log_slx_pars
		- calculate the age-dependent selectivity pattern
		- fill selectivty array for that block.
		- selectivity is scaled to have a mean = 1 across all ages.
	*/
	dvariable p1,p2;
	dvar_vector slx(sage,nage);
	log_slx.initialize();
	
	for(int h = 1; h <= nSlxBlks; h++){

		switch(nSelType(h)){
			case 1: //logistic
				p1  = mfexp(log_slx_pars(h,1,1));
				p2  = mfexp(log_slx_pars(h,1,2));
				slx = plogis(age,p1,p2);
			break;
		}
		

		for(int i = nslx_syr(h); i <= nslx_nyr(h); i++){
			log_slx(i) = log(slx) - log(mean(slx));
		}
	}
	Sij.sub(mod_syr,mod_nyr) = mfexp(log_slx);


FUNCTION void initializeStateVariables()
	/**
		- Set initial values for numbers-at-age matrix in first year
		  and sage recruits for all years.
		*/
	Nij.initialize();

	// initialize first row of numbers-at-age matrix
	dvar_vector lx(sage,nage);
	for(int j = sage; j <= nage; j++){

		lx(j) = exp(-Mij(mod_syr,j)*(j-sage));
		if( j==nage ) lx(j) /= (1.0-exp(-Mij(mod_syr,j)));

		if( j > sage ){
			Nij(mod_syr)(j) = mfexp(log_rinit + log_rinit_devs(j)) * lx(j);			
		}
	} 


	// iniitialize first columb of numbers-at-age matrix
	for(int i = mod_syr; i <= mod_nyr; i++){
		Nij(i,sage) = exp(log_rbar + log_rbar_devs(i));
	}
	//COUT(lx);
	//COUT(Nij);


FUNCTION void updateStateVariables()
	/**
		- Update the numbers-at-age conditional on the catch-at-age.
		- Assume a pulse fishery.
		- step 1 » calculate a vector of vulnerable-numbers-at-age
		- step 2 » calculate vulnerable proportions-at-age.
		- step 3 » calc average weight of catch (wbar) conditional on Qij.
		- step 4 » calc catch-at-age | catch in biomass Cij = Ct/wbar * Qij.
		- step 5 » update numbers-at-age (using a very dangerous difference eqn.)
		*/

		Qij.initialize();
		Cij.initialize();
		Pij.initialize();
		dvariable wbar;		// average weight of the catch.
		dvar_vector vj(sage,nage);
		dvar_vector pj(sage,nage);
		dvar_vector sj(sage,nage);
		

		for(int i = mod_syr; i <= mod_nyr; i++){

			// step 1.
			vj = elem_prod(Nij(i),Sij(i));

			// step 2.
			Qij(i) = vj / sum(vj);

			// step 3.
			dvector wa = data_cm_waa(i)(sage,nage);
			wbar = wa * Qij(i);

			// step 4.
			Cij(i) = data_catch(i,2) / wbar * Qij(i);

			// step 5.
			sj = mfexp(-Mij(i));
			Pij(i) = Nij(i) - Cij(i); // should use posfun here
			Nij(i+1)(sage+1,nage) =++ elem_prod(Pij(i)(sage,nage-1),sj(sage,nage-1));
			Nij(i+1)(nage) += Pij(i,nage) * sj(nage);
		}
		// COUT(Nij)
		// cross check... Looks good.
		// COUT(Cij(mod_syr) * data_cm_waa(mod_syr)(sage,nage));


FUNCTION void calcSpawningStockRecruitment()
	/**
		- The functional form of the stock recruitment model follows that of a 
			Ricker model, where R = so * SSB * exp(-beta * SSB).  The two parameters
			so and beta where previously estimated as free parameters in the old
			herring model.  Herein this fucntion I derive so and beta from the 
			leading parameters Ro and reck; Ro is the unfished sage recruits, and reck
			is the recruitment compensation parameter, or the relative improvement in
			juvenile survival rates as the spawning stock SSB tends to 0.  Simply a 
			multiple of the replacement line Ro/Bo.

			At issue here is time varying maturity and time-varying natural mortality.
			When either of these two variables are assumed to change over time, then
			the underlying stock recruitment relationship will also change. This 
			results in a non-stationary distribution.  For the purposes of this 
			assessment model, I use the average mortality and maturity schedules to
			derive the spawning boimass per recruit, which is ultimately used in 
			deriving the parameters for the stock recruitment relationship.
		*/
	for(int i = mod_syr; i <= mod_nyr; i++){
		Oij(i) = elem_prod(mat(i),Nij(i));
		ssb(i) = (Oij(i) - Cij(i)) * data_sp_waa(i)(sage,nage);
	}

	// average natural mortality
	dvar_vector mbar(sage,nage);
	int n = Mij.rowmax() - Mij.rowmin() + 1;
	mbar  = colsum(Mij)/n;

	// average maturity
	dvar_vector mat_bar(sage,nage);
	mat_bar = colsum(mat)/n;

	COUT(mat(mod_syr));
	COUT(mat_bar)
	exit(1);

	// Ricker stock-recruitment function 
	// so = reck/phiE; where reck > 1.0
	// beta = log(reck)/(ro * phiE)
	dvariable reck = mfexp(log_reck);




FUNCTION void calcAgeCompResiduals()
	/**
		- Commercial catch-age comp residuals
		- Spawning survey catch-age comp residuals.
		*/

		dvar_matrix pred_cm_comp(mod_syr,mod_nyr,sage,nage);
		dvar_matrix resd_cm_comp(mod_syr,mod_nyr,sage,nage);
		dvar_matrix pred_sp_comp(mod_syr,mod_nyr,sage,nage);
		dvar_matrix resd_sp_comp(mod_syr,mod_nyr,sage,nage);

		resd_cm_comp.initialize();
		resd_sp_comp.initialize();
		for(int i = mod_syr; i <= mod_nyr; i++){
			
			// commercial age-comp prediction 
			pred_cm_comp(i) = Qij(i);
			if( data_cm_comp(i,sage) >= 0 ){
				resd_cm_comp(i) = data_cm_comp(i)(sage,nage) - pred_cm_comp(i);
			}

			// spawning age-comp prediction
			pred_sp_comp(i) = Oij(i) / sum(Oij(i));
			if( data_sp_comp(i,sage) >= 0 ){
				resd_sp_comp(i) = data_sp_comp(i)(sage,nage) - pred_sp_comp(i);
			}
		}

		//COUT(resd_sp_comp);

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


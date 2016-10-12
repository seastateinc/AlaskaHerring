	#include <admodel.h>
	#include <string.h>
	#include <time.h>
	#undef EOF
	#define EOF 999
	#undef TINY
	#define TINY  1.0e-10
	#undef REPORT
	#define REPORT(object) report << #object "\n" << setw(8) \
	<< setprecision(4) << setfixed() << object << endl;
	#undef COUT
	#define COUT(object) cout << #object "\n" << setw(6) \
	<< setprecision(3) << setfixed() << object << endl;
	template<typename T>
	dvar_vector plogis(const dvector x, T location, T scale)
	{
		return(1.0 / (1.0 + mfexp(-(x-location)/scale)));
	}
	dvar_vector plogis(const dvector x, dvariable location, dvariable scale)
	{
		return(1.0 / (1.0 + mfexp(-(x-location)/scale)));
	}
	dvar_vector plogis(const dvector x, prevariable location, prevariable scale)
	{
		return(1.0 / (1.0 + mfexp(-(x-location)/scale)));
	}
	template<typename T>
	dvar_vector plogis95(const dvector x, T a50, T a95)
	{
		return(1.0 / (1.0 + mfexp(-log(19)*(x-a50)/(a95-a50))));
	}
	dvariable dmvlogistic(const dmatrix o, const dvar_matrix& p,dvar_matrix& nu, double& tau2,const double minp)
	{ //returns the negative loglikelihood using the MLE for the variance
	/*
		This is a modified version of the dmvlogistic negative log likelihood
		where proportions at age less than minp are pooled into the consecutive 
		age-classes.  See last paragraph in Appendix A of Richards, Schnute and
		Olsen 1997. 
		
		NB minp must be greater than 0, otherwise this algorithm returns an
		error if one of the observed proportions is zero.
		
		-1) first count the number of observations for each year > minp
		-2) normalized observed and predicted age-proportions
		-3) loop over ages, and check if observed proportion is < minp
				-if yes, then add observed proprtion to bin k
				-if no then add observed proportion to bin k and increment
				 bin k if k is currently less than the number of bins.
		-4) do the same grouping for the predicted proportions.
		-5) use ivector iiage to correctly assign residuals into nu
		
		FEB 8, 2011.  Fixed a bug in the variance calculation & 
		likelihood scaling that was discovered at the 2011 Hake
		assessment STAR panel in Seattle.
	*/
	RETURN_ARRAYS_INCREMENT();
	int i,j,k,n;
	int age_counts=0;
	int a = o.colmin();
	int A=o.colmax();
	double tiny=0.001/(A-a+1);
	int t=o.rowmin();
	int T=o.rowmax();
	dvariable tau_hat2;   //mle of the variance
	//dvar_matrix nu(t,T,a,A);
	nu.initialize();
	
	for(i=t; i<=T; i++)
	{ 
		n=0;
		dvector oo = o(i)/sum(o(i));
		dvar_vector pp = p(i)/sum(p(i));
		
		//count # of observations greater than minp (2% is a reasonable number)
		for(j=a;j<=A;j++)
			if(oo(j) > minp)n++;
		
		ivector iiage(1,n);
		dvector o1(1,n); o1.initialize();
		dvar_vector p1(1,n); p1.initialize();
		k=1;
		for(j=a;j<=A;j++)
		{
			if(oo(j)<=minp)
			{
				o1(k)+=oo(j);
				p1(k)+=pp(j);
			}
			else
			{
				o1(k)+=oo(j);
				p1(k)+=pp(j);
				if(k<=n)iiage(k)=j;   //ivector for the grouped residuals
				if(k<n) k++;
			}
		}
		
		//assign residuals to nu based on iiage index
		dvar_vector t1 = log(o1)-log(p1) - mean(log(o1)-log(p1));
		
		for(j=1;j<=n;j++)
			nu(i)(iiage(j))=t1(j);
		
		age_counts += n-1;
	}
	//Depricated  Wrong Variance & likelihood calculation.
	//tau_hat2 = 1./(age_counts*(T-t+1))*norm2(nu);
	//dvariable nloglike =(age_counts*(T-t+1))*log(tau_hat2);
	
	//Feb 8, 2011  Fixed variance & likelihood
	tau_hat2 = 1./(age_counts)*norm2(nu);
	dvariable nloglike =(age_counts)*log(tau_hat2);
	tau2=value(tau_hat2); //mle of the variance 
	RETURN_ARRAYS_DECREMENT();
	return(nloglike);
	}
	double dicValue;
	double dicNoPar;
 	void function_minimizer::mcmc_eval(void)
 	{
 		// |---------------------------------------------------------------------------|
 		// | Added DIC calculation.  Martell, Jan 29, 2013                             |
 		// |---------------------------------------------------------------------------|
 		// | DIC = pd + dbar
 		// | pd  = dbar - dtheta  (Effective number of parameters)
 		// | dbar   = expectation of the likelihood function (average f)
 		// | dtheta = expectation of the parameter sample (average y) 
 	  gradient_structure::set_NO_DERIVATIVES();
 	  initial_params::current_phase=initial_params::max_number_phases;
 	  uistream * pifs_psave = NULL;
 	#if defined(USE_LAPLACE)
 	#endif
 	#if defined(USE_LAPLACE)
 	    initial_params::set_active_random_effects();
 	    int nvar1=initial_params::nvarcalc(); 
 	#else
 	  int nvar1=initial_params::nvarcalc(); // get the number of active parameters
 	#endif
 	  int nvar;
	  
 	  pifs_psave= new
 	    uistream((char*)(ad_comm::adprogram_name + adstring(".psv")));
 	  if (!pifs_psave || !(*pifs_psave))
 	  {
 	    cerr << "Error opening file "
 	            << (char*)(ad_comm::adprogram_name + adstring(".psv"))
 	       << endl;
 	    if (pifs_psave)
 	    {
 	      delete pifs_psave;
 	      pifs_psave=NULL;
 	      return;
 	    }
 	  }
 	  else
 	  {     
 	    (*pifs_psave) >> nvar;
 	    if (nvar!=nvar1)
 	    {
 	      cout << "Incorrect value for nvar in file "
 	           << "should be " << nvar1 << " but read " << nvar << endl;
 	      if (pifs_psave)
 	      {
 	        delete pifs_psave;
 	        pifs_psave=NULL;
 	      }
 	      return;
 	    }
 	  }
  
 	  int nsamp = 0;
 	  double sumll = 0;
 	  independent_variables y(1,nvar);
 	  independent_variables sumy(1,nvar);
 	  do
 	  {
 	    if (pifs_psave->eof())
 	    {
 	      break;
 	    }
 	    else
 	    {
 	      (*pifs_psave) >> y;
 	      sumy = sumy + y;
 	      if (pifs_psave->eof())
 	      {
 	      	double dbar = sumll/nsamp;
 	      	int ii=1;
 	      	y = sumy/nsamp;
 	      	initial_params::restore_all_values(y,ii);
 	        initial_params::xinit(y);   
 	        double dtheta = 2.0 * get_monte_carlo_value(nvar,y);
 	        double pd     = dbar - dtheta;
 	        double dic    = pd + dbar;
 	        dicValue      = dic;
 	        dicNoPar      = pd;
 	        cout<<"Number of posterior samples    = "<<nsamp    <<endl;
 	        cout<<"Expectation of log-likelihood  = "<<dbar     <<endl;
 	        cout<<"Expectation of theta           = "<<dtheta   <<endl;
 	        cout<<"Number of estimated parameters = "<<nvar1    <<endl;
 		    	cout<<"Effective number of parameters = "<<dicNoPar <<endl;
 		    	cout<<"DIC                            = "<<dicValue <<endl;
 	        break;
 	      }
 	      int ii=1;
 	      initial_params::restore_all_values(y,ii);
 	      initial_params::xinit(y);   
 	      double ll = 2.0 * get_monte_carlo_value(nvar,y);
 	      sumll    += ll;
 	      nsamp++;
 	      // cout<<sumy(1,3)/nsamp<<" "<<get_monte_carlo_value(nvar,y)<<endl;
 	    }
 	  }
 	  while(1);
 	  if (pifs_psave)
 	  {
 	    delete pifs_psave;
 	    pifs_psave=NULL;
 	  }
 	  return;
 	}
#include <admodel.h>
#include <contrib.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <gdbprintlib.cpp>

#include <ham.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
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
			retro_yrs = 0;
			if ( (on=option_match(ad_comm::argc,ad_comm::argv,"-retro")) > -1 )
			{
				retro_yrs = atoi(ad_comm::argv[on+1]);
				cout<<"|—————————————————————————————————————————————————|\n";
				cout<<"| Implementing Retrospective analysis             |\n";
				cout<<"|—————————————————————————————————————————————————|\n";
				cout<<"| Number of retrospective years = "<<retro_yrs<<endl;
			}
		}
		
  DEBUG_FLAG.allocate("DEBUG_FLAG");
  DataFile.allocate("DataFile");
  ControlFile.allocate("ControlFile");
 ad_comm::change_datafile_name(DataFile);
  dat_syr.allocate("dat_syr");
  dat_nyr.allocate("dat_nyr");
  mod_syr.allocate("mod_syr");
  mod_nyr.allocate("mod_nyr");
  sage.allocate("sage");
  nage.allocate("nage");
  rec_syr = mod_syr + sage;
  age.allocate(sage,nage);
 age.fill_seqadd(sage,1);
  nFecBlocks.allocate("nFecBlocks");
  nFecBlockYears.allocate(1,nFecBlocks,"nFecBlockYears");
  fec_slope.allocate(1,nFecBlocks,"fec_slope");
  fec_inter.allocate(1,nFecBlocks,"fec_inter");
  data_ct_raw.allocate(dat_syr,dat_nyr,1,3,"data_ct_raw");
  data_sp_waa.allocate(dat_syr,dat_nyr,sage-1,nage,"data_sp_waa");
  data_cm_waa.allocate(dat_syr,dat_nyr,sage-1,nage,"data_cm_waa");
  data_cm_comp.allocate(dat_syr,dat_nyr,sage-1,nage,"data_cm_comp");
  data_sp_comp.allocate(dat_syr,dat_nyr,sage-1,nage,"data_sp_comp");
  data_egg_dep.allocate(dat_syr,dat_nyr,1,3,"data_egg_dep");
  data_mileday.allocate(dat_syr,dat_nyr,1,3,"data_mileday");
  avg_sp_waa.allocate(sage,nage);
		int n = data_sp_waa.rowmax() - data_sp_waa.rowmin() + 1;
		avg_sp_waa = colsum(data_sp_waa)(sage,nage) / n;
  Eij.allocate(mod_syr,mod_nyr,sage,nage);
		int iyr = mod_syr;
		
		for(int h = 1; h <= nFecBlocks; h++){
			do{
				Eij(iyr) = 1.e-6 *
								(data_sp_waa(iyr)(sage,nage) * fec_slope(h) - fec_inter(h));
				iyr ++;
			}while(iyr <= nFecBlockYears(h));
		}
  dat_eof.allocate("dat_eof");
 if(dat_eof != 999){cout<<"Error reading data file, aborting."<<endl; exit(1);}
 ad_comm::change_datafile_name(ControlFile);
 n_theta = 6;
  theta_DM.allocate(1,n_theta,1,7,"theta_DM");
  theta_ival.allocate(1,n_theta);
  theta_lb.allocate(1,n_theta);
  theta_ub.allocate(1,n_theta);
  theta_phz.allocate(1,n_theta);
  theta_iprior.allocate(1,n_theta);
  theta_p1.allocate(1,n_theta);
  theta_p2.allocate(1,n_theta);
 theta_ival = column(theta_DM,1);
 theta_lb  = column(theta_DM,2);
 theta_ub  = column(theta_DM,3);
 theta_phz = ivector(column(theta_DM,4));
 theta_iprior = ivector(column(theta_DM,5));
 theta_p1 = column(theta_DM,6);
 theta_p2 = column(theta_DM,7);
  nMatBlocks.allocate("nMatBlocks");
  maturity_cont.allocate(1,4,1,nMatBlocks,"maturity_cont");
  mat_a50.allocate();
  mat_a95.allocate();
  mat_phz.allocate(1,nMatBlocks);
  nMatBlockYear.allocate(1,nMatBlocks);
		mat_a50 = maturity_cont(1);
		mat_a95 = maturity_cont(2);
		mat_phz = ivector(maturity_cont(3));
		nMatBlockYear = ivector(maturity_cont(4));
  mort_type.allocate("mort_type");
  mort_dev_phz.allocate("mort_dev_phz");
  nMortBlocks.allocate("nMortBlocks");
  nMortBlockYear.allocate(1,nMortBlocks,"nMortBlockYear");
 nSlxCols = 9;
  nSlxBlks.allocate("nSlxBlks");
  selex_cont.allocate(1,nSlxBlks,1,nSlxCols,"selex_cont");
  nSelType.allocate(1,nSlxBlks);
  nslx_phz.allocate(1,nSlxBlks);
  nslx_rows.allocate(1,nSlxBlks);
  nslx_cols.allocate(1,nSlxBlks);
  nslx_syr.allocate(1,nSlxBlks);
  nslx_nyr.allocate(1,nSlxBlks);
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
  nMiscCont.allocate("nMiscCont");
  dMiscCont.allocate(1,nMiscCont,"dMiscCont");
  data_catch.allocate(dat_syr,dat_nyr,1,3);
		data_catch = data_ct_raw;
		for( int i = dat_syr; i <= dat_nyr; i++ ) {
			data_catch(i,2) = dMiscCont(1) * data_ct_raw(i,2);
		}
		if(dMiscCont(2)) cout<<"Condition model on Ft"<<endl;
  ctl_eof.allocate("ctl_eof");
		if(ctl_eof != 999){
			cout<<"Error reading control file, aborting."<<ctl_eof<<endl; 
			exit(1);
		}
 mod_nyr = mod_nyr - retro_yrs;
 nf = 0;
}

void model_parameters::initializationfunction(void)
{
  theta.set_initial_value(theta_ival);
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  theta.allocate(1,n_theta,theta_lb,theta_ub,theta_phz,"theta");
  log_natural_mortality.allocate("log_natural_mortality");
  #ifndef NO_AD_INITIALIZE
  log_natural_mortality.initialize();
  #endif
  log_rinit.allocate("log_rinit");
  #ifndef NO_AD_INITIALIZE
  log_rinit.initialize();
  #endif
  log_rbar.allocate("log_rbar");
  #ifndef NO_AD_INITIALIZE
  log_rbar.initialize();
  #endif
  log_ro.allocate("log_ro");
  #ifndef NO_AD_INITIALIZE
  log_ro.initialize();
  #endif
  log_reck.allocate("log_reck");
  #ifndef NO_AD_INITIALIZE
  log_reck.initialize();
  #endif
  log_sigma_r.allocate("log_sigma_r");
  #ifndef NO_AD_INITIALIZE
  log_sigma_r.initialize();
  #endif
  log_rinit_devs.allocate(sage+1,nage,-15.0,15.0,2,"log_rinit_devs");
  log_rbar_devs.allocate(mod_syr,mod_nyr+1,-15.0,15.0,2,"log_rbar_devs");
  mat_params.allocate(1,nMatBlocks,1,2,0,100,mat_phz,"mat_params");
  mat.allocate(mod_syr,mod_nyr,sage,nage,"mat");
  #ifndef NO_AD_INITIALIZE
    mat.initialize();
  #endif
		cout<<"Good to here"<<endl;
		cout<<mat_params(1)<<endl;
		if( !global_parfile ) {
			for(int h = 1; h <= nMatBlocks; h++){
				mat_params(h,1) = mat_a50(h);
				mat_params(h,2) = mat_a95(h);				
			}
		}
		cout<<mat_params(1)<<endl;
		
  log_m_devs.allocate(1,nMortBlocks,-15.0,15.0,mort_dev_phz,"log_m_devs");
  Mij.allocate(mod_syr,mod_nyr,sage,nage,"Mij");
  #ifndef NO_AD_INITIALIZE
    Mij.initialize();
  #endif
  log_slx_pars.allocate(1,nSlxBlks,1,nslx_rows,1,nslx_cols,-25,25,nslx_phz,"log_slx_pars");
  log_slx.allocate(mod_syr,mod_nyr,sage,nage,"log_slx");
  #ifndef NO_AD_INITIALIZE
    log_slx.initialize();
  #endif
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
 int phz; phz = dMiscCont(2)==0?-1:1;
  log_ft_pars.allocate(mod_syr,mod_nyr,-30.,3.0,phz,"log_ft_pars");
		if(b_simulation_flag) log_ft_pars = log(0.1);
  ro.allocate("ro");
  #ifndef NO_AD_INITIALIZE
  ro.initialize();
  #endif
  reck.allocate("reck");
  #ifndef NO_AD_INITIALIZE
  reck.initialize();
  #endif
  so.allocate("so");
  #ifndef NO_AD_INITIALIZE
  so.initialize();
  #endif
  beta.allocate("beta");
  #ifndef NO_AD_INITIALIZE
  beta.initialize();
  #endif
  fore_sb.allocate("fore_sb");
  #ifndef NO_AD_INITIALIZE
  fore_sb.initialize();
  #endif
  fore_vb.allocate("fore_vb");
  #ifndef NO_AD_INITIALIZE
  fore_vb.initialize();
  #endif
  ghl.allocate("ghl");
  #ifndef NO_AD_INITIALIZE
  ghl.initialize();
  #endif
  ssb.allocate(mod_syr,mod_nyr,"ssb");
  #ifndef NO_AD_INITIALIZE
    ssb.initialize();
  #endif
  recruits.allocate(rec_syr,mod_nyr+1,"recruits");
  #ifndef NO_AD_INITIALIZE
    recruits.initialize();
  #endif
  spawners.allocate(rec_syr,mod_nyr+1,"spawners");
  #ifndef NO_AD_INITIALIZE
    spawners.initialize();
  #endif
  resd_rec.allocate(rec_syr,mod_nyr+1,"resd_rec");
  #ifndef NO_AD_INITIALIZE
    resd_rec.initialize();
  #endif
  pred_egg_dep.allocate(mod_syr,mod_nyr,"pred_egg_dep");
  #ifndef NO_AD_INITIALIZE
    pred_egg_dep.initialize();
  #endif
  resd_egg_dep.allocate(mod_syr,mod_nyr,"resd_egg_dep");
  #ifndef NO_AD_INITIALIZE
    resd_egg_dep.initialize();
  #endif
  pred_mileday.allocate(mod_syr,mod_nyr,"pred_mileday");
  #ifndef NO_AD_INITIALIZE
    pred_mileday.initialize();
  #endif
  resd_mileday.allocate(mod_syr,mod_nyr,"resd_mileday");
  #ifndef NO_AD_INITIALIZE
    resd_mileday.initialize();
  #endif
  pred_catch.allocate(mod_syr,mod_nyr,"pred_catch");
  #ifndef NO_AD_INITIALIZE
    pred_catch.initialize();
  #endif
  resd_catch.allocate(mod_syr,mod_nyr,"resd_catch");
  #ifndef NO_AD_INITIALIZE
    resd_catch.initialize();
  #endif
  Nij.allocate(mod_syr,mod_nyr+1,sage,nage,"Nij");
  #ifndef NO_AD_INITIALIZE
    Nij.initialize();
  #endif
  Oij.allocate(mod_syr,mod_nyr+1,sage,nage,"Oij");
  #ifndef NO_AD_INITIALIZE
    Oij.initialize();
  #endif
  Pij.allocate(mod_syr,mod_nyr+1,sage,nage,"Pij");
  #ifndef NO_AD_INITIALIZE
    Pij.initialize();
  #endif
  Sij.allocate(mod_syr,mod_nyr+1,sage,nage,"Sij");
  #ifndef NO_AD_INITIALIZE
    Sij.initialize();
  #endif
  Qij.allocate(mod_syr,mod_nyr+1,sage,nage,"Qij");
  #ifndef NO_AD_INITIALIZE
    Qij.initialize();
  #endif
  Cij.allocate(mod_syr,mod_nyr+1,sage,nage,"Cij");
  #ifndef NO_AD_INITIALIZE
    Cij.initialize();
  #endif
  Fij.allocate(mod_syr,mod_nyr+1,sage,nage,"Fij");
  #ifndef NO_AD_INITIALIZE
    Fij.initialize();
  #endif
  pred_cm_comp.allocate(mod_syr,mod_nyr,sage,nage,"pred_cm_comp");
  #ifndef NO_AD_INITIALIZE
    pred_cm_comp.initialize();
  #endif
  resd_cm_comp.allocate(mod_syr,mod_nyr,sage,nage,"resd_cm_comp");
  #ifndef NO_AD_INITIALIZE
    resd_cm_comp.initialize();
  #endif
  pred_sp_comp.allocate(mod_syr,mod_nyr,sage,nage,"pred_sp_comp");
  #ifndef NO_AD_INITIALIZE
    pred_sp_comp.initialize();
  #endif
  resd_sp_comp.allocate(mod_syr,mod_nyr,sage,nage,"resd_sp_comp");
  #ifndef NO_AD_INITIALIZE
    resd_sp_comp.initialize();
  #endif
  f.allocate("f");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  fpen.allocate("fpen");
  #ifndef NO_AD_INITIALIZE
  fpen.initialize();
  #endif
  sd_terminal_ssb.allocate("sd_terminal_ssb");
  sd_forecast_ssb.allocate("sd_forecast_ssb");
  sd_projected_ssb.allocate("sd_projected_ssb");
}

void model_parameters::preliminary_calculations(void)
{

#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
	/* 
	 * SIMULATION MODEL SWITCH
	 */
	if( b_simulation_flag && rseed > 0) {
		cout<<"|--------------------------|"<<endl;
		cout<<"| RUNNING SIMULATION MODEL |"<<endl;
		cout<<"|--------------------------|"<<endl;
		
		char type;
		do
		{
				cout<<"| Continue? [y]es or [n]o "<<endl;
				cin >> type;
		}
		while( !cin.fail() && type!='y' && type!='n' );
		if( type =='y' ){
			runSimulationModel(rseed);
		} else {
			exit(1);
		}
	} else if ( b_simulation_flag && rseed < 0){
		runSimulationModel(rseed);
	}
}

void model_parameters::userfunction(void)
{
  f =0.0;
	initializeModelParameters();
	if(DEBUG_FLAG) cout<<"--> Ok after initializeModelParameters      <--"<<endl;
	initializeMaturitySchedules();
	if(DEBUG_FLAG) cout<<"--> Ok after initializeMaturitySchedules    <--"<<endl;
	calcNaturalMortality();
	if(DEBUG_FLAG) cout<<"--> Ok after calcNaturalMortality           <--"<<endl;
	calcSelectivity();
	if(DEBUG_FLAG) cout<<"--> Ok after calcSelectivity                <--"<<endl;
	if( dMiscCont(2) ) {
		calcFishingMortalitiy();
		if(DEBUG_FLAG) cout<<"--> Ok after calcFishingMortalitiy          <--"<<endl;  
	}
	initializeStateVariables();
	if(DEBUG_FLAG) cout<<"--> Ok after initializeStateVariables       <--"<<endl;
	updateStateVariables();
	if(DEBUG_FLAG) cout<<"--> Ok after updateStateVariables           <--"<<endl;
	calcSpawningStockRecruitment();
	if(DEBUG_FLAG) cout<<"--> Ok after calcSpawningStockRecruitment   <--"<<endl;
	calcAgeCompResiduals();
	if(DEBUG_FLAG) cout<<"--> Ok after calcAgeCompResiduals           <--"<<endl;
	calcEggSurveyResiduals();
	if(DEBUG_FLAG) cout<<"--> Ok after calcEggSurveyResiduals         <--"<<endl;
	calcMiledaySurveyResiduals();
	if(DEBUG_FLAG) cout<<"--> Ok after calcMiledaySurveyResiduals     <--"<<endl;
	//if( dMiscCont(2) ) {
		calcCatchResiduals();
		if(DEBUG_FLAG) cout<<"--> Ok after calcCatchResiduals           <--"<<endl; 
	// }
	calcObjectiveFunction(); nf++;
	if(DEBUG_FLAG) cout<<"--> Ok after calcObjectiveFunction          <--"<<endl;
	sd_terminal_ssb = ssb(mod_nyr);
	if( last_phase() ) {
		runForecast();
	}
	if(mceval_phase()) {
		writePosteriorSamples();
	}
}

void model_parameters::runForecast()
{
	/** 
	Conduct a 1-year ahead forcasts based on SR
	PSUEDOCODE:
		-1 declare variables for forcasting.
			* recruitment, spawning biomass, numbers-at-age
			* selectivity curve
		-2 Calculate F-at-age conditional on harvest rule.
			* harvest_rate = 20% || set TAC option
		-3 Update state variables from pyr=mod_nyr+1 to pyr+2
			* recruitment based on ssb(pyr-sage)
		-4 Compute GHLs given threshold and target harvest rate
			* user specifies threshold and harvest rate in control file.
	**/
	int nyr = mod_nyr;
	int pyr = nyr+1;
	dvariable fore_rt;	// sage recruits
	dvar_vector fore_nj(sage,nage);	//numbers-at-age
	dvar_vector fore_cj(sage,nage); //catch-at-age
	fore_rt = so * ssb(pyr-sage) * exp(-beta*ssb(pyr-sage));
	fore_nj = Nij(pyr); fore_nj(sage) = fore_rt;
	fore_vb = fore_nj * elem_prod(Sij(nyr),data_cm_waa(nyr)(sage,nage));
	fore_sb = fore_nj * elem_prod(mat(nyr),data_sp_waa(nyr)(sage,nage));
	sd_forecast_ssb = fore_sb;
	// GHL for pyr
	double ssb_threshold = dMiscCont(3);
	double target_hr = dMiscCont(4);
	dvariable hr = (2.0 + 8.0*fore_sb / dMiscCont(5))/100.0;
	if( hr > target_hr) {
		hr = target_hr;
	} else if( fore_sb < ssb_threshold ) {
		hr = 0.0;
	}
	ghl = hr * fore_sb;
	//cout<<"harvest rate = "<<hr<<" GHL = "<<ghl<<endl;
	// update state variables to pyr+1 so you can predict
	// the effect of the 2016 fishery on the 2017 spawning stock.
	// predicted catch-at-age
	dvar_vector pa = elem_prod(fore_nj,Sij(nyr));
	pa /= sum(pa);
	dvariable wbar = pa * data_cm_waa(nyr)(sage,nage);
	fore_cj = ghl/wbar * pa;
	fore_nj = elem_prod(fore_nj - fore_cj,mfexp(-Mij(nyr)));
	fore_nj(sage) = so * ssb(pyr+1-sage) * exp(-beta*ssb(pyr-sage));
	fore_sb = fore_nj * elem_prod(mat(nyr),data_sp_waa(nyr)(sage,nage));
	sd_projected_ssb = fore_sb;
}

void model_parameters::writePosteriorSamples()
{
	/**
	- This function is only envoked when the -mceval
		command line option is implemented.
	*/
	if(nf==1){
		ofstream ofs("ssb.ps");
	}
	ofstream ofs("ssb.ps",ios::app);
	ofs<<ssb<<endl;
}

void model_parameters::runSimulationModel(const int& rseed)
{
	/*
		PSUEDOCODE:
		1) initialize model parameters based on pin file, or control file.
		2) initialize maturity schedules for all block years.
		3) initialize natural mortality schedules for all block years.
		4) calculate selectivity parameters.
		5) generate random normal deviates for process/observation errors.
		6) initialize state variables.
		7) update state variables conditioned on the observed catch data.
		8) calculate age-composition residuals and over-write input data in memory.
		9) calculate egg survey residuals and simulate fake survey.
	*/
	if(global_parfile) {
		cout<<"\nUsing pin file for simulation parameter values.\n"<<endl;
	}
	// 1) initialize model parameters based on pin file, or control file.
	initializeModelParameters();
	// 2) initialize maturity schedules for all block years.
	initializeMaturitySchedules();
	// 3) initialize natural mortality schedules for all block years.
	calcNaturalMortality();
	// 4) calculate selectivity parameters.
	calcSelectivity();
	// 4b) calculate fishing mortality.
	calcFishingMortalitiy();
	// 5) generate random normal deviates for process errors:
	//    - log_m_devs
	//    - log_rinit_devs
	//    - log_rbar_devs
	random_number_generator rng(rseed);
	dvector epsilon_m_devs(1,nMortBlocks);
	dvector epsilon_rbar_devs(mod_syr,mod_nyr+1);
	dvector epsilon_rinit_devs(sage+1,nage);
	double sigma_m_devs = 0.1;
	double sigma_rbar_devs = 0.4;
	double sigma_rinit_devs = 0.4;
	epsilon_m_devs.fill_randn(rng);
	epsilon_rbar_devs.fill_randn(rng);
	epsilon_rinit_devs.fill_randn(rng);
	log_m_devs = dvar_vector(epsilon_m_devs * sigma_m_devs - 0.5 * square(sigma_m_devs));
	log_rbar_devs = dvar_vector(epsilon_rbar_devs * sigma_rbar_devs - 0.5 * square(sigma_rbar_devs));
	log_rinit_devs = dvar_vector(epsilon_rinit_devs * sigma_rinit_devs - 0.5 * square(sigma_rinit_devs));
	// Not sure if the following should be done. It should produce a less
	// biased MLE of the average recruitment, but uncertainty is biased downwards.
	// ensure random deviates satisfy sum dev = 0 constraint.
	// log_m_devs -= mean(log_m_devs);
	// log_rbar_devs -= mean(log_rbar_devs);
	// log_rinit_devs -= mean(log_rinit_devs);
	// 6) initialize state variables.
	initializeStateVariables();
	// 7) update state variables conditioned on the observed catch data.
	updateStateVariables();
	// 7a) calcSpawningRecruitment
	calcSpawningStockRecruitment();
	// 8) calculate age-composition residuals and over-write input data in memory.
	//    - Note the age-comps are sampled from a multivariate logistic dist.
	calcAgeCompResiduals();
	for(int i = mod_syr; i <= mod_nyr; i++) {
		dvector t1 = value(pred_cm_comp(i));
		dvector t2 = value(pred_sp_comp(i));
		// cout<<t2<<endl;
		data_cm_comp(i)(sage,nage) = rmvlogistic(t1,0.30,rseed + i);
		data_sp_comp(i)(sage,nage) = rmvlogistic(t2,0.30,rseed - i);
	}
	// 9) calculate egg survey residuals and simulate fake survey.
	calcEggSurveyResiduals();
	dvector epsilon_egg_dep(mod_syr,mod_nyr);
	epsilon_egg_dep.fill_randn(rng);
	for(int i = mod_syr; i <= mod_nyr; i++) {
		if( data_egg_dep(i,2) > 0 ) {
			double sd = data_egg_dep(i,3);
			data_egg_dep(i,2) = value(pred_egg_dep(i)) 
													* exp(epsilon_egg_dep(i) * sd - 0.5 * square(sd));  
		}
	}
	// 10) calculate mile_days
	calcMiledaySurveyResiduals();
	dvector epsilon_mileday(mod_syr,mod_nyr);
	epsilon_mileday.fill_randn(rng);
	for(int i = mod_syr; i <= mod_nyr; i++) {
		if( data_mileday(i,2) > 0 ) {
			double sd = data_mileday(i,3);
			data_mileday(i,2) = value(pred_mileday(i))
													* exp(epsilon_mileday(i) * sd - 0.5 * square(sd));
		}
	}
}

void model_parameters::initializeModelParameters()
{
	fpen = 0;
	log_natural_mortality = theta(1);
	log_rinit             = theta(2);
	log_rbar              = theta(3);
	log_ro                = theta(4);
	log_reck              = theta(5);
	log_sigma_r						= theta(6);
	// COUT(theta);
}

void model_parameters::initializeMaturitySchedules()
{
	int iyr = mod_syr;
	int jyr;
	mat.initialize();
	for(int h = 1; h <= nMatBlocks; h++) {
		dvariable mat_a50 = mat_params(h,1);
		dvariable mat_a95 = mat_params(h,2);
		jyr = h != nMatBlocks ? nMatBlockYear(h) : nMatBlockYear(h)-retro_yrs;
		// if( h != nMatBlocks ){
		// 	jyr = nMatBlockYear(h);
		// } else {
		// 	jyr = nMatBlockYear(h) - retro_yrs;
		// }
		// fill maturity array using logistic function
		do{
			mat(iyr++) = plogis95(age,mat_a50,mat_a95);
		} while(iyr <= jyr);  
	}
}

void model_parameters::calcNaturalMortality()
{
	int iyr = mod_syr;
	int jyr;
	Mij.initialize();
	switch(mort_type) {
		case 1:	// constant M within block.
			for(int h = 1; h <= nMortBlocks; h++){
				dvariable mi = mfexp(log_natural_mortality + log_m_devs(h));
				jyr = h != nMortBlocks?nMortBlockYear(h):nMortBlockYear(h)-retro_yrs;
				// fill mortality array by block
				do{
					Mij(iyr++) = mi;
				} while(iyr <= jyr);
			}
		break;
		case 2: // cubic spline  
			dvector iiyr = dvector((nMortBlockYear - mod_syr)/(mod_nyr-mod_syr));
			dvector jjyr(mod_syr,mod_nyr);
			jjyr.fill_seqadd(0,double(1.0/(mod_nyr-mod_syr)));
			dvar_vector mi = log_natural_mortality + log_m_devs;
			vcubic_spline_function cubic_spline_m(iiyr,mi);
			dvar_vector mtmp = cubic_spline_m(jjyr);
			for(int i=mod_syr; i <= mod_nyr; i++)
			{
				Mij(i) = mfexp(mtmp(i));
			}
		break;
	}
}

void model_parameters::calcSelectivity()
{
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
				slx = plogis(age,p1,p2) + TINY;
			break;
		}
		int jyr = h != nSlxBlks ? nslx_nyr(h):nslx_nyr(h)-retro_yrs;
		for(int i = nslx_syr(h); i <= jyr; i++){
			log_slx(i) = log(slx) - log(mean(slx));
		}
	}
	Sij.sub(mod_syr,mod_nyr) = mfexp(log_slx);
}

void model_parameters::calcFishingMortalitiy()
{
	/**
		- Calculate Fishing mortality, and then Zij = Mij + Fij
		*/
	for(int i = mod_syr; i <= mod_nyr; i++) {
		Fij(i) = exp(log_ft_pars(i)) * Sij(i);
	}
}

void model_parameters::initializeStateVariables()
{
	/**
		- Set initial values for numbers-at-age matrix in first year
			and sage recruits for all years.
		*/
	Nij.initialize();
	// initialize first row of numbers-at-age matrix
	// lx is a vector of survivorship (probability of surviving to age j)
	dvar_vector lx(sage,nage);
	for(int j = sage; j <= nage; j++){
		lx(j) = exp(-Mij(mod_syr,j)*(j-sage));
		if( j==nage ) lx(j) /= (1.0-exp(-Mij(mod_syr,j)));
		if( j > sage ){
			Nij(mod_syr)(j) = mfexp(log_rinit + log_rinit_devs(j)) * lx(j);     
		}
	} 
	// iniitialize first column of numbers-at-age matrix
	for(int i = mod_syr; i <= mod_nyr + 1; i++){
		Nij(i,sage) = mfexp(log_rbar + log_rbar_devs(i));
	}
	//COUT(lx);
	//COUT(Nij);
}

void model_parameters::updateStateVariables()
{
	/**
	- Update the numbers-at-age conditional on the catch-at-age.
	- Assume a pulse fishery.
	- step 1 calculate a vector of vulnerable-numbers-at-age
	- step 2 calculate vulnerable proportions-at-age.
	- step 3 calc average weight of catch (wbar) conditional on Qij.
	- step 4 calc catch-at-age | catch in biomass Cij = Ct/wbar * Qij.
	- step 5 condition on Ft or else condition on observed catch.
	- step 6 update numbers-at-age (using a very dangerous difference eqn.)
	*/
	Qij.initialize();
	Cij.initialize();
	Pij.initialize();
	dvariable wbar;   // average weight of the catch.
	dvar_vector vj(sage,nage);
	//dvar_vector pj(sage,nage);
	dvar_vector sj(sage,nage);
	for(int i = mod_syr; i <= mod_nyr; i++){
		// step 1.
		vj = elem_prod(Nij(i),Sij(i));
		// step 2.
		Qij(i) = vj / sum(vj);
		// ADF&G's approach.
		if( !dMiscCont(2) ) {
			// step 3.
			dvector wa = data_cm_waa(i)(sage,nage);
			wbar = wa * Qij(i);
			// step 4.
			Cij(i) = data_catch(i,2) / wbar * Qij(i);
			// should use posfun here
			Pij(i) = posfun(Nij(i) - Cij(i),0.01,fpen); 
			// step 6. update numbers at age    
			sj = mfexp(-Mij(i));
			Nij(i+1)(sage+1,nage) =++ elem_prod(Pij(i)(sage,nage-1),
																					sj(sage,nage-1));
			Nij(i+1)(nage) += Pij(i,nage) * sj(nage);       
		} 
		// step 5.
		// Condition on Ft
		else { 
			// add flexibility here for Popes approx, or different seasons               
			Pij(i) = elem_prod( Nij(i), exp(-Fij(i)) );
			Cij(i) = elem_prod( Nij(i), 1.-exp(-Fij(i)));
			// the following assumes fishery is at the start of the year.
			dvar_vector zi = Mij(i) + Fij(i);
			Nij(i+1)(sage+1,nage) = ++ elem_prod(Nij(i)(sage,nage-1),
																					 mfexp(-zi(sage,nage-1)));
			Nij(i+1)(nage) += Nij(i,nage) * mfexp(-zi(nage));
		}
	}
	// cross check... Looks good.
	// COUT(Cij(mod_syr) * data_cm_waa(mod_syr)(sage,nage));
}

void model_parameters::calcSpawningStockRecruitment()
{
	/**
	The functional form of the stock recruitment model follows that of a 
	Ricker model, where R = so * SSB * exp(-beta * SSB).  The two parameters
	so and beta where previously estimated as free parameters in the old
	herring model.  Herein this fucntion I derive so and beta from the 
	leading parameters Ro and reck; Ro is the unfished sage recruits, and reck
	is the recruitment compensation parameter, or the relative improvement in
	juvenile survival rates as the spawning stock SSB tends to 0. Simply a 
	multiple of the replacement line Ro/Bo.
	At issue here is time varying maturity and time-varying natural mortality.
	When either of these two variables are assumed to change over time, then
	the underlying stock recruitment relationship will also change. This 
	results in a non-stationary distribution.  For the purposes of this 
	assessment model, I use the average mortality and maturity schedules to
	derive the spawning boimass per recruit, which is ultimately used in 
	deriving the parameters for the stock recruitment relationship.
		*/
	/*
		Spoke to Sherri about this. Agreed to change the equation to prevent
	*/
	for(int i = mod_syr; i <= mod_nyr; i++){
		//Oij(i) = elem_prod(mat(i),Nij(i));
		//ssb(i) = (Oij(i) - Cij(i)) * data_sp_waa(i)(sage,nage);
		Oij(i) = elem_prod(mat(i),Nij(i)-Cij(i));
		ssb(i) = Oij(i) * data_sp_waa(i)(sage,nage);
	}
	// average natural mortality
	dvar_vector mbar(sage,nage);
	int n = Mij.rowmax() - Mij.rowmin() + 1;
	mbar  = colsum(Mij)/n;
	// average maturity
	dvar_vector mat_bar(sage,nage);
	mat_bar = colsum(mat)/n;
	// unfished spawning biomass per recruit
	dvar_vector lx(sage,nage);
	lx(sage) = 1.0;
	for(int j = sage + 1; j <= nage; j++){
		lx(j) = lx(j-1) * mfexp(-mbar(j-1));
		if(j == nage){
			lx(j) /= 1.0 - mfexp(-mbar(j));
		}
	}
	dvariable phie = lx * elem_prod(avg_sp_waa,mat_bar);
	// Ricker stock-recruitment function 
	// so = reck/phiE; where reck > 1.0
	// beta = log(reck)/(ro * phiE)
	// Beverton Holt use:
	// beta = (reck - 1.0)/(ro *phiE)
	ro   = mfexp(log_ro);
	reck = mfexp(log_reck) + 1.0;
	so   = reck/phie;
	beta = log(reck) / (ro * phie);
	spawners = ssb(mod_syr,mod_nyr-sage+1).shift(rec_syr);
	recruits = elem_prod( so*spawners , mfexp(-beta*spawners) );
	resd_rec = log(column(Nij,sage)(rec_syr,mod_nyr+1)+TINY) 
						 - log(recruits+TINY);
}

void model_parameters::calcAgeCompResiduals()
{
	/**
	- Commercial catch-age comp residuals
	- Spawning survey catch-age comp residuals.
	*/
	resd_cm_comp.initialize();
	resd_sp_comp.initialize();
	for(int i = mod_syr; i <= mod_nyr; i++){
		// commercial age-comp prediction 
		pred_cm_comp(i) = Qij(i);
		if( data_cm_comp(i,sage) >= 0 ){
			resd_cm_comp(i) = data_cm_comp(i)(sage,nage) - pred_cm_comp(i);
		}
		// spawning age-comp prediction
		pred_sp_comp(i) = (Oij(i)+TINY) / sum(Oij(i)+TINY);
		if( data_sp_comp(i,sage) >= 0 ){
			resd_sp_comp(i) = data_sp_comp(i)(sage,nage) - pred_sp_comp(i);
		}
	}
	//COUT(resd_sp_comp);
}

void model_parameters::calcEggSurveyResiduals()
{
	/**
	- Observed egg data is in trillions of eggs
	- Predicted eggs is the mature female numbers-at-age multiplied 
		by the fecundity-at-age, which comes from a regession of 
		fecundity = slope * obs_sp_waa - intercept
	- Note Eij is the Fecundity-at-age j in year i.
	- 
	*/
	resd_egg_dep.initialize();
	for(int i = mod_syr; i <= mod_nyr; i++){
		pred_egg_dep(i) = (0.5 * Oij(i)) * Eij(i);
		if(data_egg_dep(i,2) > 0){
			resd_egg_dep(i) = log(data_egg_dep(i,2)) - log(pred_egg_dep(i));
		}
	}
}

void model_parameters::calcMiledaySurveyResiduals()
{
	/**
	- Assumed index from aerial survey is a relative abundance
	index. The slope of the regression ln(SSB) = q * ln(MileMilt) + 0
	is computed on the fly in case of missing data. 
	- See Walters and Ludwig 1994 for more details.
	*/
	resd_mileday.initialize();
	pred_mileday.initialize();
	int n = 1;
	dvar_vector zt(mod_syr,mod_nyr); zt.initialize();
	dvariable zbar = 0;
	for(int i = mod_syr; i <= mod_nyr; i++){
		if(data_mileday(i,2) > 0){
			zt(i) = log(data_mileday(i,2)) - log(ssb(i));
			zbar  = zt(i)/n + zbar *(n-1)/n;
			n++ ;
		}   
	}
	pred_mileday = ssb * exp(zbar);
	resd_mileday = zt - zbar;
	// COUT(resd_mileday);
}

void model_parameters::calcCatchResiduals()
{
	/**
	- Catch residuals assuming a lognormal error structure.
	*/
	pred_catch.initialize();
	resd_catch.initialize();
	for(int i = mod_syr; i <= mod_nyr; i++) {
		if(data_catch(i,2) > 0) {
			pred_catch(i) = Cij(i) * data_cm_waa(i)(sage,nage);
			resd_catch(i) = log(data_catch(i,2)) - log(pred_catch(i));
		}
	}
}

void model_parameters::calcObjectiveFunction()
{
	/**
	- THIS FUNCTION IS ORGANIZED AS FOLLOWS:
		1. Negative logliklihoods (nll)
		2. Penalized loglikelihoods (penll)
	*/
	// 1. Negative loglikelihoods
	dvar_vector nll(1,7);
	nll.initialize();
	// Mulitvariate logistic likelihood for composition data.
	double sp_tau2;
	double minp = 0.00;
	dmatrix d_sp_comp = trans(trans(data_sp_comp).sub(sage,nage)).sub(mod_syr,mod_nyr);
	nll(1) = dmvlogistic(d_sp_comp,pred_sp_comp,resd_sp_comp,sp_tau2,minp);
	// Mulitvariate logistic likelihood for composition data.
	double cm_tau2;
	dmatrix d_cm_comp  = trans(trans(data_cm_comp).sub(sage,nage)).sub(mod_syr,mod_nyr);
	nll(2) = dmvlogistic(d_cm_comp,pred_cm_comp,resd_cm_comp,cm_tau2,minp);
	// Negative loglikelihood for egg deposition data
	dvector std_egg_dep = TINY + column(data_egg_dep,3)(mod_syr,mod_nyr);
	nll(3) = dnorm(resd_egg_dep,std_egg_dep);
	// Negative loglikelihood for milt mile day
	dvector std_mileday = TINY + column(data_mileday,3)(mod_syr,mod_nyr);
	nll(4) = dnorm(resd_mileday,std_mileday);
	// Negative loglikelihood for stock-recruitment data
	dvariable std_rec = sqrt(1.0/exp(log_sigma_r));
	nll(5) = dnorm(resd_rec,std_rec);
	// Negative loglikelihood for catch data
	dvector std_catch = column(data_catch,3)(mod_syr,mod_nyr);
	nll(6) = dnorm(resd_catch,std_catch);
	// 2. Penalized logliklelihoods
	dvar_vector penll(1,4);
	penll.initialize();
	// | Terminal phase of estimation.
	dvariable log_fbar = mean(log_ft_pars);
	if( last_phase() ){
		penll(1) = dnorm(log_rinit_devs,0.0,5.0);
		penll(2) = dnorm(log_rbar_devs,0.0,5.0);
		penll(3) = dnorm(log_fbar,0.2,2.0);
	} else {
		penll(1) = dnorm(log_rinit_devs,0.0,1.0);
		penll(2) = dnorm(log_rbar_devs,0.0,1.0);
		penll(3) = dnorm(log_fbar,0.2,0.07);      
	}
	//if( dMiscCont(2) ) nll(7) = norm2(resd_catch);
	//// cout<<(data_sp_comp)<<endl;
	//f  = sum(nll) + 1000.0 * fpen;
	//if( dMiscCont(2) ){
	//  if( !last_phase() ) {
	//    f += 100. * square(mean(log_ft_pars)-log(0.2));
	//
	//  } else {
	//    f += 100. * square(mean(log_ft_pars)-log(0.2));
	//  }
	//}
	f = sum(nll) + sum(penll) + sum(calcPriors()) + fpen;
	if(DEBUG_FLAG){
		COUT(f);
		COUT(nll);
		COUT(penll);
		COUT(calcPriors());
		if(fpen > 0 ){COUT(fpen);}
	}  
	//FUNCTION dvariable dnorm(const dvar_vector residual, const dvector std)
	//RETURN_ARRAYS_INCREMENT();
	//double pi=3.141593;
	//int n=size_count(residual);
	//dvector var=elem_prod(std,std);
	//dvar_vector SS=elem_prod(residual,residual);
	//RETURN_ARRAYS_DECREMENT();
	//return(0.5*n*log(2.*pi)+sum(log(std))+sum(elem_div(SS,2.*var)));
}

dvar_vector model_parameters::calcPriors()
{
	// |---------------------------------------------------------------------------------|
	// | PRIORS FOR LEADING PARAMETERS p(theta)
	// |---------------------------------------------------------------------------------|
	// | - theta_prior is a switch to determine which statistical distribution to use.
	// |
	dvariable ptmp; 
	dvar_vector priors(1,n_theta);
	priors.initialize();
	for(int i=1;i<=n_theta;i++)
	{
		ptmp = 0;
		// for(j=1;j<=ipar_vector(i);j++)
		// {
			if( active(theta(i)) )
			{
				switch(theta_iprior(i))
				{
				case 1:   //normal
					ptmp += dnorm(theta(i),theta_p1(i),theta_p2(i));
					break;
				case 2:   //lognormal CHANGED RF found an error in dlnorm prior. rev 116
					ptmp += dlnorm(theta(i),theta_p1(i),theta_p2(i));
					break;
				case 3:   //beta distribution (0-1 scale)
					double lb,ub;
					lb=theta_lb(i);
					ub=theta_ub(i);
					ptmp += dbeta((theta(i)-lb)/(ub-lb),theta_p1(i),theta_p2(i));
					break;
				case 4:   //gamma distribution
					ptmp += dgamma(theta(i),theta_p1(i),theta_p2(i));
					break;
				default:  //uniform density
					if(theta_lb(i) > theta_ub(i)){
						ptmp += log(1./(theta_lb(i)-theta_ub(i)));
					} else {
						ptmp += log(1./(theta_ub(i)-theta_lb(i)));
					}
					break;
				}
			}
		// }
		priors(i) = ptmp; 
	}
	return(priors);
}

void model_parameters::report(const dvector& gradients)
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
	REPORT(dat_syr);
	REPORT(dat_nyr);
	REPORT(mod_syr);
	REPORT(mod_nyr);
	REPORT(sage);
	REPORT(nage);
	REPORT(nFecBlocks);
	REPORT(nFecBlockYears);
	REPORT(fec_slope);
	REPORT(fec_inter);
	REPORT(data_ct_raw);
	REPORT(data_sp_waa);
	REPORT(data_cm_waa);
	REPORT(data_cm_comp);
	REPORT(data_sp_comp);
	REPORT(data_egg_dep);
	REPORT(data_mileday);
	REPORT(EOF)
	ivector year(mod_syr,mod_nyr);
	year.fill_seqadd(mod_syr,1);
	ivector years(mod_syr,mod_nyr+1);
	years.fill_seqadd(mod_syr,1);
	ivector rec_years(rec_syr,mod_nyr+1);
	rec_years.fill_seqadd(rec_syr,1);
	ivector iage = ivector(age);
	REPORT(iage);
	REPORT(year);
	REPORT(years);
	REPORT(rec_years);
	REPORT(data_catch);
	REPORT(pred_catch);
	REPORT(ssb);
	REPORT(spawners);
	REPORT(recruits);
	REPORT(pred_egg_dep);
	REPORT(Nij);
	REPORT(Oij);
	REPORT(Pij);
	REPORT(Cij);
	REPORT(Sij);
	REPORT(Qij);
	REPORT(Mij);
	REPORT(resd_egg_dep);
	REPORT(resd_rec);
	REPORT(resd_catch);
	REPORT(resd_sp_comp);
	REPORT(resd_cm_comp);
	REPORT(fore_sb);
	REPORT(fore_vb);
	REPORT(ghl);
	REPORT(theta_ival);
}

void model_parameters::final_calcs()
{
	//system("cp ham.rep run1.rep");
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::set_runtime(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}

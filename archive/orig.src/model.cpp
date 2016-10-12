	#include <admodel.h>
	#include <string.h>
	#include <time.h>
        adstring model_name;
        adstring data_file;
	time_t start,finish;
	long hour,minute,second;
	double elapsed_time;
	adstring BaseFileName;
	adstring ReportFileName;
	adstring NewFileName;
	adstring stripExtension(adstring fileName)
	{
		/*
		This function strips the file extension
		from the fileName argument and returns
		the file name without the extension.
		*/
		const int length = fileName.size();
		for (int i=length; i>=0; --i)
		{
			if (fileName(i)=='.')
			{
				return fileName(1,i-1);
			}
		}
		return fileName;
	}
#include <admodel.h>
#include <contrib.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <model.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  pad_evalout = new ofstream("evalout.prj");;
  DataFile.allocate("DataFile");
  ControlFile.allocate("ControlFile");
  Graphics.allocate("Graphics");
 BaseFileName = stripExtension(DataFile);  
 ReportFileName = BaseFileName + adstring(".rep");
 cout<<"You are modeling the "<<BaseFileName<<" stock of herring"<<endl;
 cout<<""<<endl;
 time(&start);
 ad_comm::change_datafile_name(ControlFile);
  nages.allocate("nages");
  dat_styr.allocate("dat_styr");
  dat_endyr.allocate("dat_endyr");
  mod_styr.allocate("mod_styr");
  mod_endyr.allocate("mod_endyr");
  Year.allocate(mod_styr,mod_endyr+1);
  dyrs=dat_endyr-dat_styr+1;
  myrs=mod_endyr-mod_styr+1;
  md_offset=mod_styr-dat_styr;
  Year.fill_seqadd(mod_styr,1);   
  Thresh.allocate("Thresh");
  Thresh_denom.allocate("Thresh_denom");
  Threshold.allocate(1,myrs,"Threshold");
  fw_a_a.allocate(1,nages,"fw_a_a");
cout <<myrs<<endl; 
  F_Bk.allocate("F_Bk");
  F_Bk_Yrs.allocate(1,F_Bk+1,"F_Bk_Yrs");
  F_slope.allocate(1,F_Bk,"F_slope");
  F_inter.allocate(1,F_Bk,"F_inter");
  F_Bk_Idx.allocate(1,F_Bk+1);
  S_Bk.allocate("S_Bk");
  S_Bk_Yrs.allocate(1,S_Bk+1,"S_Bk_Yrs");
  S_Bk_Idx.allocate(1,S_Bk+1);
  mat_Bk.allocate("mat_Bk");
  mat_Bk_Yrs.allocate(1,mat_Bk+1,"mat_Bk_Yrs");
  mat_Bk_Idx.allocate(1,mat_Bk+1);
  gs_Bk.allocate("gs_Bk");
  gs_Bk_Yrs.allocate(1,gs_Bk+1,"gs_Bk_Yrs");
  gs_Bk_Idx.allocate(1,gs_Bk+1);
  ph_Int.allocate("ph_Int");
  ph_S.allocate("ph_S");
  ph_mat_a.allocate("ph_mat_a");
  ph_gs_a.allocate("ph_gs_a");
  ph_gs_b.allocate("ph_gs_b");
  ph_mat_b.allocate("ph_mat_b");
  ph_Rec.allocate("ph_Rec");
  ph_Ric.allocate("ph_Ric");
  ph_md.allocate("ph_md");
  lC.allocate("lC");
  lS.allocate("lS");
  lR.allocate("lR");
  lE.allocate(1,dyrs,"lE");
  lM.allocate(1,dyrs,"lM");
  eof1.allocate("eof1");
    if(eof1==42) cout << BaseFileName<<".ctl has been read correctly!"<<endl;
    else 
    {    
         cout <<"|----------------------------------------------------------------------|"<<endl;   
         cout <<"|      Red alert! Captain to bridge! The .ctl file is compromised!     |"<<endl;
         cout <<"|----------------------------------------------------------------------|"<<endl; 
	 cout <<"|      Last integer read is "<<eof1<<", but the file *should* end with 42      |"<<endl;
         cout <<"| Please check the .ctl file for errors and make sure the above calls  |"<<endl;
         cout <<"|              are matched exactly by the file's contents              |"<<endl;
         cout <<"|----------------------------------------------------------------------|"<<endl;
    exit(1); 
    }
   // | The lines below populate the indices for fecundity, mortality,
   // | selectivity and maturity
   for (int i=1;i<=F_Bk+1;i++)
     {
       F_Bk_Idx(i)=F_Bk_Yrs(i)-mod_styr+1;
     } 
   for (int i=1;i<=S_Bk+1;i++)
     {
       S_Bk_Idx(i)=S_Bk_Yrs(i)-mod_styr+1;
     }
   for (int i=1;i<=mat_Bk+1;i++)
     {
       mat_Bk_Idx(i)=mat_Bk_Yrs(i)-mod_styr+1;
     }
   for (int i=1;i<=gs_Bk+1;i++)
     {
       gs_Bk_Idx(i)=gs_Bk_Yrs(i)-mod_styr+1;
     }
 ad_comm::change_datafile_name(DataFile);
  tcb.allocate(1,dyrs,"tcb");
  obs_sp_waa.allocate(1,dyrs,1,nages,"obs_sp_waa");
  obs_c_waa.allocate(1,dyrs,1,nages,"obs_c_waa");
  obs_c_comp.allocate(1,dyrs,1,nages,"obs_c_comp");
  obs_sp_comp.allocate(1,dyrs,1,nages,"obs_sp_comp");
  tot_obs_egg.allocate(1,dyrs,"tot_obs_egg");
  mile_days.allocate(1,dyrs,"mile_days");
  tcbm.allocate(1,dyrs);
  eof2.allocate("eof2");
    if(eof2==42) cout << BaseFileName<<".dat has been read correctly!"<<endl;
    else 
    {       
         cout <<"|----------------------------------------------------------------------|"<<endl;   
         cout <<"|   ** Red alert! Captain to bridge! The .dat file is compromised! **  |"<<endl;
         cout <<"|----------------------------------------------------------------------|"<<endl; 
	 cout <<"|      Last integer read is "<<eof2<<", but the file *should* end with 42      |"<<endl;
         cout <<"| Please check the .dat file for errors and make sure the above calls  |"<<endl;
         cout <<"|              are matched exactly by the file's contents              |"<<endl;
         cout <<"|----------------------------------------------------------------------|"<<endl;
    exit(1); 
    }
    // |This line converts annual catch biomass from short tons to metric tons
    for (int i=1;i<=dyrs;i++)
    {
       tcbm(i)=0.90718*tcb(i);
    } 
    // |This line extracts the last line of spawning weight-at-age for forecasting
    //fw_a_a = obs_sp_waa(1,dyrs);
   
 ad_comm::change_datafile_name(Graphics);
  mod_yrs.allocate(1,myrs,"mod_yrs");
  yminusthree.allocate(1,myrs,"yminusthree");
  yminustwo.allocate(1,myrs,"yminustwo");
  yminusone.allocate(1,myrs,"yminusone");
  yminusthreeFOR.allocate("yminusthreeFOR");
  yminustwoFOR.allocate("yminustwoFOR");
  yminusoneFOR.allocate("yminusoneFOR");
  eof3.allocate("eof3");
    if(eof3==42) cout << BaseFileName<<"_graphics.ctl has been read correctly!"<<endl;
    else 
    {       
         cout <<"|----------------------------------------------------------------------|"<<endl;   
         cout <<"|   Red alert! Captain to bridge! The graphics file is compromised!    |"<<endl;
         cout <<"|----------------------------------------------------------------------|"<<endl; 
	 cout <<"|      Last integer read is "<<eof3<<", but the file *should* end with 42      |"<<endl;
         cout <<"|  Please check the graphics file for errors and make sure the above   |"<<endl;
         cout <<"|          calls are matched exactly by the file's contents            |"<<endl;
         cout <<"|----------------------------------------------------------------------|"<<endl;
    exit(1); 
    }
   
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  init_age_3.allocate(1,myrs,0,2500,ph_Rec,"init_age_3");
  init_pop.allocate(1,nages-1,0,500,ph_Int,"init_pop");
  log_alpha.allocate(-10,0,ph_Ric,"log_alpha");
  log_beta.allocate(ph_Ric,"log_beta");
  alpha.allocate("alpha");
  #ifndef NO_AD_INITIALIZE
  alpha.initialize();
  #endif
  beta.allocate("beta");
  #ifndef NO_AD_INITIALIZE
  beta.initialize();
  #endif
  M.allocate(1,S_Bk,.0001,1,ph_S,"M");
  Sur.allocate(1,myrs,1,nages,"Sur");
  #ifndef NO_AD_INITIALIZE
    Sur.initialize();
  #endif
  S.allocate(1,S_Bk,"S");
  #ifndef NO_AD_INITIALIZE
    S.initialize();
  #endif
  S_for.allocate("S_for");
  #ifndef NO_AD_INITIALIZE
  S_for.initialize();
  #endif
  mat_a.allocate(1,mat_Bk,1,10,ph_mat_a,"mat_a");
  mat_b.allocate(1,mat_Bk,0,5,ph_mat_b,"mat_b");
  gs_a.allocate(1,gs_Bk,2,10,ph_gs_a,"gs_a");
  gs_b.allocate(1,gs_Bk,0,5,ph_gs_b,"gs_b");
  md_c.allocate(200,600,ph_Int,"md_c");
  GS.allocate(1,myrs,1,nages,"GS");
  #ifndef NO_AD_INITIALIZE
    GS.initialize();
  #endif
  GS_Sc.allocate(1,myrs,1,nages,"GS_Sc");
  #ifndef NO_AD_INITIALIZE
    GS_Sc.initialize();
  #endif
  Mat.allocate(1,myrs,1,nages,"Mat");
  #ifndef NO_AD_INITIALIZE
    Mat.initialize();
  #endif
  int1.allocate(1,myrs,"int1");
  #ifndef NO_AD_INITIALIZE
    int1.initialize();
  #endif
  int2.allocate(1,myrs,1,nages,"int2");
  #ifndef NO_AD_INITIALIZE
    int2.initialize();
  #endif
  int3.allocate(1,myrs,"int3");
  #ifndef NO_AD_INITIALIZE
    int3.initialize();
  #endif
  mat_for.allocate(1,nages,"mat_for");
  #ifndef NO_AD_INITIALIZE
    mat_for.initialize();
  #endif
  naa.allocate(1,myrs,1,nages,"naa");
  #ifndef NO_AD_INITIALIZE
    naa.initialize();
  #endif
  sel_naa.allocate(1,myrs,1,nages,"sel_naa");
  #ifndef NO_AD_INITIALIZE
    sel_naa.initialize();
  #endif
  sel_naa_prop.allocate(1,myrs,1,nages,"sel_naa_prop");
  #ifndef NO_AD_INITIALIZE
    sel_naa_prop.initialize();
  #endif
  est_c_naa.allocate(1,myrs,1,nages,"est_c_naa");
  #ifndef NO_AD_INITIALIZE
    est_c_naa.initialize();
  #endif
  est_sp_naa.allocate(1,myrs,1,nages,"est_sp_naa");
  #ifndef NO_AD_INITIALIZE
    est_sp_naa.initialize();
  #endif
  est_sp_comp.allocate(1,myrs,1,nages,"est_sp_comp");
  #ifndef NO_AD_INITIALIZE
    est_sp_comp.initialize();
  #endif
  est_sp_baa.allocate(1,myrs,1,nages,"est_sp_baa");
  #ifndef NO_AD_INITIALIZE
    est_sp_baa.initialize();
  #endif
  est_mat_naa.allocate(1,myrs,1,nages,"est_mat_naa");
  #ifndef NO_AD_INITIALIZE
    est_mat_naa.initialize();
  #endif
  est_mat_baa.allocate(1,myrs,1,nages,"est_mat_baa");
  #ifndef NO_AD_INITIALIZE
    est_mat_baa.initialize();
  #endif
  post_naa.allocate(1,myrs,1,nages,"post_naa");
  #ifndef NO_AD_INITIALIZE
    post_naa.initialize();
  #endif
  est_egg_naa.allocate(1,myrs,1,nages,"est_egg_naa");
  #ifndef NO_AD_INITIALIZE
    est_egg_naa.initialize();
  #endif
  tot_sel_N.allocate(1,myrs,"tot_sel_N");
  #ifndef NO_AD_INITIALIZE
    tot_sel_N.initialize();
  #endif
  tot_sp_B.allocate(1,myrs,"tot_sp_B");
  #ifndef NO_AD_INITIALIZE
    tot_sp_B.initialize();
  #endif
  tot_mat_N.allocate(1,myrs,"tot_mat_N");
  #ifndef NO_AD_INITIALIZE
    tot_mat_N.initialize();
  #endif
  tot_sp_N.allocate(1,myrs,"tot_sp_N");
  tot_mat_B.allocate(1,myrs,"tot_mat_B");
  tot_post_N.allocate(1,myrs,"tot_post_N");
  tot_est_egg.allocate(1,myrs,"tot_est_egg");
  N.allocate(1,myrs,"N");
  SR.allocate(1,myrs-3,"SR");
  #ifndef NO_AD_INITIALIZE
    SR.initialize();
  #endif
  M_D.allocate(1,myrs,"M_D");
  #ifndef NO_AD_INITIALIZE
    M_D.initialize();
  #endif
  for_naa.allocate(1,nages,"for_naa");
  #ifndef NO_AD_INITIALIZE
    for_naa.initialize();
  #endif
  for_mat_naa.allocate(1,nages,"for_mat_naa");
  #ifndef NO_AD_INITIALIZE
    for_mat_naa.initialize();
  #endif
  for_mat_baa.allocate(1,nages,"for_mat_baa");
  #ifndef NO_AD_INITIALIZE
    for_mat_baa.initialize();
  #endif
  for_mat_prop.allocate(1,nages,"for_mat_prop");
  #ifndef NO_AD_INITIALIZE
    for_mat_prop.initialize();
  #endif
  for_mat_b_prop.allocate(1,nages,"for_mat_b_prop");
  #ifndef NO_AD_INITIALIZE
    for_mat_b_prop.initialize();
  #endif
  for_mat_B.allocate("for_mat_B");
  #ifndef NO_AD_INITIALIZE
  for_mat_B.initialize();
  #endif
  for_mat_B_st.allocate("for_mat_B_st");
  #ifndef NO_AD_INITIALIZE
  for_mat_B_st.initialize();
  #endif
  for_tot_mat_N.allocate("for_tot_mat_N");
  #ifndef NO_AD_INITIALIZE
  for_tot_mat_N.initialize();
  #endif
  HR.allocate("HR");
  #ifndef NO_AD_INITIALIZE
  HR.initialize();
  #endif
  HR_p.allocate("HR_p");
  #ifndef NO_AD_INITIALIZE
  HR_p.initialize();
  #endif
  GHL.allocate("GHL");
  #ifndef NO_AD_INITIALIZE
  GHL.initialize();
  #endif
  FIGDATA.allocate(1,myrs,1,42,"FIGDATA");
  #ifndef NO_AD_INITIALIZE
    FIGDATA.initialize();
  #endif
  FIGDATAAGE.allocate(1,nages,1,3,"FIGDATAAGE");
  #ifndef NO_AD_INITIALIZE
    FIGDATAAGE.initialize();
  #endif
  res_c_comp.allocate(1,myrs,1,nages,"res_c_comp");
  #ifndef NO_AD_INITIALIZE
    res_c_comp.initialize();
  #endif
  res_sp_comp.allocate(1,myrs,1,nages,"res_sp_comp");
  #ifndef NO_AD_INITIALIZE
    res_sp_comp.initialize();
  #endif
  res_tot_egg.allocate(1,myrs,"res_tot_egg");
  #ifndef NO_AD_INITIALIZE
    res_tot_egg.initialize();
  #endif
  res_SR.allocate(1,myrs-3,"res_SR");
  #ifndef NO_AD_INITIALIZE
    res_SR.initialize();
  #endif
  SSQC.allocate("SSQC");
  #ifndef NO_AD_INITIALIZE
  SSQC.initialize();
  #endif
  SSQSp.allocate("SSQSp");
  #ifndef NO_AD_INITIALIZE
  SSQSp.initialize();
  #endif
  wSSQE.allocate(1,myrs,"wSSQE");
  #ifndef NO_AD_INITIALIZE
    wSSQE.initialize();
  #endif
  WSSQE.allocate("WSSQE");
  #ifndef NO_AD_INITIALIZE
  WSSQE.initialize();
  #endif
  SSQR.allocate("SSQR");
  #ifndef NO_AD_INITIALIZE
  SSQR.initialize();
  #endif
  M_DR.allocate(1,myrs,"M_DR");
  #ifndef NO_AD_INITIALIZE
    M_DR.initialize();
  #endif
  wSSQM.allocate(1,myrs,"wSSQM");
  #ifndef NO_AD_INITIALIZE
    wSSQM.initialize();
  #endif
  WSSQM.allocate("WSSQM");
  #ifndef NO_AD_INITIALIZE
  WSSQM.initialize();
  #endif
  f.allocate("f");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  n_obs_c_comp.allocate(1,myrs,1,nages,"n_obs_c_comp");
  #ifndef NO_AD_INITIALIZE
    n_obs_c_comp.initialize();
  #endif
  n_obs_sp_comp.allocate(1,myrs,1,nages,"n_obs_sp_comp");
  #ifndef NO_AD_INITIALIZE
    n_obs_sp_comp.initialize();
  #endif
  n_tot_obs_egg.allocate(1,myrs,"n_tot_obs_egg");
  #ifndef NO_AD_INITIALIZE
    n_tot_obs_egg.initialize();
  #endif
  n_d.allocate(1,myrs+3,"n_d");
  #ifndef NO_AD_INITIALIZE
    n_d.initialize();
  #endif
  w_d.allocate(1,myrs+3,"w_d");
  #ifndef NO_AD_INITIALIZE
    w_d.initialize();
  #endif
  lnL_d.allocate(1,myrs+3,"lnL_d");
  #ifndef NO_AD_INITIALIZE
    lnL_d.initialize();
  #endif
  n_C.allocate("n_C");
  #ifndef NO_AD_INITIALIZE
  n_C.initialize();
  #endif
  n_S.allocate("n_S");
  #ifndef NO_AD_INITIALIZE
  n_S.initialize();
  #endif
  n_R.allocate("n_R");
  #ifndef NO_AD_INITIALIZE
  n_R.initialize();
  #endif
  n.allocate("n");
  #ifndef NO_AD_INITIALIZE
  n.initialize();
  #endif
  sig_1.allocate("sig_1");
  #ifndef NO_AD_INITIALIZE
  sig_1.initialize();
  #endif
  lnL.allocate("lnL");
  #ifndef NO_AD_INITIALIZE
  lnL.initialize();
  #endif
  AIC.allocate("AIC");
  #ifndef NO_AD_INITIALIZE
  AIC.initialize();
  #endif
  AICc.allocate("AICc");
  #ifndef NO_AD_INITIALIZE
  AICc.initialize();
  #endif
  p.allocate("p");
  #ifndef NO_AD_INITIALIZE
  p.initialize();
  #endif
  n_My.allocate(1,myrs,"n_My");
  #ifndef NO_AD_INITIALIZE
    n_My.initialize();
  #endif
  n_M.allocate("n_M");
  #ifndef NO_AD_INITIALIZE
  n_M.initialize();
  #endif
  sig_d.allocate(1,myrs+3,"sig_d");
  #ifndef NO_AD_INITIALIZE
    sig_d.initialize();
  #endif
  SpawnBio.allocate(1,myrs,"SpawnBio");
  SpawnBioFor.allocate("SpawnBioFor");
  GHLsd.allocate("GHLsd");
  AICcsd.allocate("AICcsd");
  //Convert observed spawner-at-age to proportions
  //for (int i=1; i<=dyrs; i++)
  //{
   //obs_sp_comp(i)/=sum(obs_sp_comp(i)+0.00001);
  //}
}

void model_parameters::preliminary_calculations(void)
{

#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
  // |---------------------------------------------------------------------------------|
  // | Needed? .PIN file better option?
  // |---------------------------------------------------------------------------------|
  log_beta=-10;                       
  mat_b=1;
}

void model_parameters::userfunction(void)
{
  f =0.0;
  ofstream& evalout= *pad_evalout;
  get_parameters();
  Time_Loop();
  get_residuals();
  evaluate_the_objective_function();
  if(last_phase())
    {
      get_forecast();
    } 
  if(sd_phase())
    {
      compute_AICc();
     get_FIGDATA();
     output_FIGDATA();
     get_FIGDATAAGE();
     output_FIGDATAAGE();
      get_sdreports();
      get_report();
    }
  if(mceval_phase())
    {
     evalout<<for_mat_B_st<<" "
          <<init_age_3<<" "
          <<init_pop<<" "<<endl;
    }
}

void model_parameters::get_parameters(void)
{
  ofstream& evalout= *pad_evalout;
  /**
    * Defines matrices for survival, maturity, and selectivity
    * based on the values entered in the Control File for:
    * S_Bk, mat_BK, and gs_BK
  **/
  Sur.initialize();
  Mat.initialize();
  GS.initialize();
  GS_Sc.initialize();
  int1.initialize();
  S_for.initialize();
  mat_for.initialize();
  //Survival
  for (int t=1;t<=S_Bk;t++)
  {
    S(t)=exp(-M(t));
      for (int i=S_Bk_Idx(t);i<=S_Bk_Idx(t+1);i++)
      {
         for (int j=1;j<=nages;j++)
         {
           Sur(i,j)=exp(-M(t));
         }
      }
  }
  S_for=Sur(myrs,1);
  //Maturity
  for (int t=1;t<=mat_Bk;t++)
  {
    for (int i=mat_Bk_Idx(t);i<=mat_Bk_Idx(t+1);i++)
       {
         for (int j=1;j<=nages;j++)
          {
            Mat(i,j)=1/(1+exp(-1.0*mat_b(t)*((j+2)-mat_a(t))));
          }
       }
   }
   mat_for = Mat(myrs);
  //Gear selectivity-at-age (scaled)
  for (int t=1;t<=gs_Bk;t++)
    {
      for (int i=gs_Bk_Idx(t);i<=gs_Bk_Idx(t+1);i++)
        {
          for (int j=1;j<=nages;j++)
            {
              GS(i,j)=1/(1+exp(-1.0*gs_b(t)*((j+2)-gs_a(t))));
            }
        }
    }
  for (int i=1;i<=myrs;i++)
    {
      int1(i)=max(GS(i));   // This is not differentiab
        for (int j=1;j<=nages;j++)
        {
          GS_Sc(i,j)=GS(i,j)/int1(i);
        }
     }
}

void model_parameters::Time_Loop(void)
{
  ofstream& evalout= *pad_evalout;
  naa.initialize();
  N.initialize();
  sel_naa.initialize();
  tot_sel_N.initialize();
  sel_naa_prop.initialize();
  est_c_naa.initialize();
  est_sp_naa.initialize();
  tot_sp_N.initialize();
  est_sp_comp.initialize();
  est_sp_baa.initialize();
  tot_sp_B.initialize();
  est_mat_baa.initialize();
  tot_mat_B.initialize();
  post_naa.initialize(); 
  tot_post_N.initialize();
  est_egg_naa.initialize();
  tot_est_egg.initialize();
  SR.initialize();
  int2.initialize();
  int3.initialize();
  M_D.initialize();
  /**
    *Pre-Fishery & Catch: This loop is the basis for the model. To perform
    *calculations we first set up the first year, get catch, then proceed to
    *calculate for remaining years
  **/
  //----------------------------------------------------------------------------
  // YEAR ONE
  //----------------------------------------------------------------------------
  
  for(int i=1;i<=1;i++)
   {
        naa(i,1)=init_age_3(i);
     for(int j=2;j<=nages;j++)
       {  //recruitment vector - year 1
         naa(i,j)=init_pop(j-1);              //initial population - year 1
       }
     for(int j=1;j<=nages;j++)
       {
         sel_naa(i,j)=naa(i,j)*GS_Sc(i,j);    //numbers-at-age vulnerable to gear (make sure this statement is separate from above)
       }
   
     tot_sel_N=rowsum(sel_naa);               //total numbers vulnerable to gear
   
    for(int j=1;j<=nages;j++)
     {
       sel_naa_prop(i,j)=sel_naa(i,j)/tot_sel_N(i);          //proportion-at-age: gear
       int2(i,j)=sel_naa_prop(i,j)*obs_c_waa(i+md_offset,j); //proporton * wt_g
     }
     int3=rowsum(int2);
    for(int j=1;j<=nages;j++)
      {
        est_c_naa(i,j)=tcbm(i+md_offset)*sel_naa_prop(i,j)/int3(i); //est. caa as % from tcbm
      }
    for(int j=1;j<=nages;j++)
      {
        post_naa(i,j)=naa(i,j)-est_c_naa(i,j); //numbers - catch
      }
  } 
  //----------------------------------------------------------------------------
  // END FIRST YEAR LOOP
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  // ALL OTHER YEARS
  //----------------------------------------------------------------------------
  for(int i=2;i<=myrs;i++)
    {
    
    for(int j=2;j<=nages;j++)
      {   
        naa(i,1)=init_age_3(i);                 //recruitment vector, years > 1     
        naa(i,j)=post_naa(i-1,j-1)*Sur(i,j-1);  //naa: (numbers-catch)*survival
      }
    for(int j=nages;j<=nages;j++)
      {
        naa(i,j)=post_naa(i-1,j-1)*Sur(i,j-1)+post_naa(i-1,j)*Sur(i,j); //+ class, naa
      }
    for(int j=1;j<=nages;j++)
      {
        sel_naa(i,j)=naa(i,j)*GS_Sc(i,j); //naa vulnerable to gear
      }
    
     tot_sel_N=rowsum(sel_naa);
    for(int j=1;j<=nages;j++)
      {
        sel_naa_prop(i,j)=sel_naa(i,j)/tot_sel_N(i);           //proportion-at-age: gear
        int2(i,j)=sel_naa_prop(i,j)*obs_c_waa(i+md_offset,j);  //proporton * wt_kg
      }
      int3=rowsum(int2);
    for(int j=1;j<=nages;j++)
      {
        est_c_naa(i,j)=tcbm(i+md_offset)*sel_naa_prop(i,j)/int3(i); //est. caa as % from tcb
      }
    for(int j=1;j<=nages;j++)
      {
        post_naa(i,j)=naa(i,j)-est_c_naa(i,j); //numbers-catch
      }
      
  }  //End other years Loop
    N=rowsum(naa);                  // total abundance (millions)
    tot_post_N = rowsum(post_naa);  // total abundance [numbers - catch] (millions)
  
  /**
    *spawning population quantities
    *Recall that "spawning population" is distinct from "mature population"
    *because spawning = mature - fishery
    *If there is NO fishery, spawning population = mature population
  **/
  
  //----------------------------------------------------------------------------
  // SPAWNERS & MATURE
  //----------------------------------------------------------------------------
  for (int i=1;i<=myrs;i++)
    {
      
      for(int j=1;j<=nages;j++)
        {
          est_sp_naa(i,j)=(Mat(i,j)*naa(i,j))-est_c_naa(i,j);          //spawning numbers-at-age (differs for stocks) For Craig should be Mat *(naa-est_c_na)a
        }}
    tot_sp_N=rowsum(est_sp_naa);                                       //total spawning numbers
  
  for (int i=1;i<=myrs;i++)
    {
      
      for(int j=1;j<=nages;j++)
        {
          est_sp_baa(i,j)=obs_sp_waa(i+md_offset,j)*est_sp_naa(i,j);   //spawning biomass-at-age
        }}
    tot_sp_B=rowsum(est_sp_baa);                                       //total spawning biomass
  for (int i=1;i<=myrs;i++)
    {
      
      for(int j=1;j<=nages;j++)
        {
          est_mat_naa(i,j)=(Mat(i,j)*naa(i,j));                        //mature-at-age
        }}
    tot_mat_N=rowsum(est_mat_naa);                                     //total mature numbers
      
  for (int i=1;i<=myrs;i++)
    {
      
      for(int j=1;j<=nages;j++)
        {
          est_mat_baa(i,j)=est_c_naa(i,j)*obs_c_waa(i+md_offset,j)+est_sp_naa(i,j)*obs_sp_waa(i+md_offset,j); //mature biomass at age
        }}
    tot_mat_B=rowsum(est_mat_baa);                                     //total mature biomass
  for (int i=1;i<=myrs;i++)
    {
   for(int j=1;j<=nages;j++)
        {
          est_sp_comp(i,j)=est_sp_naa(i,j)/tot_sp_N(i);                //spawning % at age
        }}
  //----------------------------------------------------------------------------
  // EGG DEPOSITION
  //----------------------------------------------------------------------------
  for (int t=1;t<=F_Bk;t++)
    {
      for (int i=F_Bk_Idx(t);i<=F_Bk_Idx(t+1);i++)
        {
          for (int j=1;j<=nages;j++)
            {
              est_egg_naa(i,j)=0.5*est_sp_naa(i,j)*(F_slope(t)*obs_sp_waa(i+md_offset,j)-F_inter(t))*0.000001;
            }
        cout<<est_egg_naa(i)<<endl;
        }
    }
    exit(1);
  tot_est_egg=rowsum(est_egg_naa);
  //----------------------------------------------------------------------------
  // SR FUNCTION
  //----------------------------------------------------------------------------
  alpha=exp(log_alpha);
  beta=exp(log_beta);
  for(int i=1;i<=myrs-3;i++)
    {
      SR(i)=alpha*tot_sp_B(i)*exp(-1.0*beta*tot_sp_B(i));
    }
 
  //Mile-days of milt
    for (int i=1;i<=myrs;i++){
    M_D(i)=0.5*tot_sp_B(i)/md_c;}
}

void model_parameters::get_residuals(void)
{
  ofstream& evalout= *pad_evalout;
  res_c_comp.initialize();
  res_sp_comp.initialize();
  res_tot_egg.initialize();
  res_SR.initialize();
  M_DR.initialize();
  //----------------------------------------------------------------------------
  // CATCH AGE COMPOSITION
  //----------------------------------------------------------------------------
  for (int i=1;i<=myrs;i++) 
    {
      for (int j=1;j<=nages;j++)
        {
          if (obs_c_comp(i+md_offset,1)<0)
            {
              res_c_comp(i,j)=0;
            }
          else
            {
              res_c_comp(i,j)=obs_c_comp(i+md_offset,j)-sel_naa_prop(i,j);
            }
        }
    }
  //----------------------------------------------------------------------------
  // SPAWNING AGE COMPOSITION
  //----------------------------------------------------------------------------
  for (int i=1;i<=myrs;i++)
    {
     for (int j=1;j<=nages;j++)
       {
         if (obs_sp_comp(i+md_offset,j)<0)
           {
             res_sp_comp(i,j)=0;
           }
        else
          {
            res_sp_comp(i,j)=obs_sp_comp(i+md_offset,j)-est_sp_comp(i,j);
          }
       }
    }
  //----------------------------------------------------------------------------
  // EGG DEPOSITION
  //----------------------------------------------------------------------------
  for (int i=1;i<=myrs;i++)
    {
      if (tot_obs_egg(i+md_offset)<0){res_tot_egg(i)=0;}
      else{res_tot_egg(i)=log(tot_obs_egg(i+md_offset))-log(tot_est_egg(i));
      wSSQE(i)=lE(i+md_offset)*square(res_tot_egg(i));
    }}
  //----------------------------------------------------------------------------
  // SPAWNER-RECRUIT FUNCTION
  //----------------------------------------------------------------------------
  for (int i=1;i<=myrs-3;i++)
   {
     res_SR(i)=log(naa(i+3,1))-log(SR(i));
   }
  //----------------------------------------------------------------------------
  // MILE DAYS OF MILT
  //----------------------------------------------------------------------------
   for (int i=1;i<=myrs;i++)
    {
      if (mile_days(i+md_offset)<0)
        {
          M_DR(i)=0;
        }
      else
        {
          M_DR(i)=log(mile_days(i+md_offset))-log(M_D(i));
        }
      wSSQM(i)=lM(i+md_offset)*square(M_DR(i));
    }
  WSSQM=sum(wSSQM);
  SSQC=norm2(res_c_comp);
  SSQSp=norm2(res_sp_comp);
  WSSQE=sum(wSSQE);
  SSQR=norm2(res_SR);
}

void model_parameters::evaluate_the_objective_function(void)
{
  ofstream& evalout= *pad_evalout;
  f=lC*SSQC+lS*SSQSp+WSSQE+lR*SSQR;
}

void model_parameters::get_forecast(void)
{
  ofstream& evalout= *pad_evalout;
  for_naa(1)=alpha*tot_sp_B(myrs-2)*exp(-1.0*beta*tot_sp_B(myrs-2)); //forecast age 3 numbers
  
  for (int j=2;j<=nages-1;j++)
    {
      for_naa(j)=post_naa(myrs,j-1)*S_for;                           //forecast naa, ages 4 - 7
    }
  for (int j=nages;j<=nages;j++)
    {
      for_naa(j)=post_naa(myrs,j-1)*S_for+post_naa(myrs,j)*S_for;    //forecast naa, age 8
    }
  for (int j=1;j<=nages;j++)
    {
      for_mat_naa(j)=for_naa(j)*Mat(myrs,j);    //forecast mature numbers at age
    }
  for_tot_mat_N=sum(for_mat_naa);
  for (int j=1;j<=nages;j++)
    {
      for_mat_baa(j)=for_mat_naa(j)*fw_a_a(j);  //forecast mature biomass at age
    }
  for_mat_B=sum(for_mat_baa);
  for (int j=1;j<=nages;j++)
    {
      for_mat_prop(j)=for_mat_naa(j)/for_tot_mat_N;  //forecast % mature at age
    }
  for (int j=1;j<=nages;j++)
    {
      for_mat_b_prop(j) = for_mat_baa(j)/for_mat_B;  //forecast % mature at age biomass
    }
  
  //GHL CALCS
  for_mat_B_st=for_mat_B/0.90718;                   //RETURN TO SHORT TONS
   HR_p=(2+8*for_mat_B_st/Thresh_denom)/100;        //Region-specific harvest rate (make sure all are in tons) CHECK!!!
  if(HR_p>0.2){HR=0.2;}
   else{HR=HR_p;}
  if(for_mat_B_st<Thresh){HR=0;}
  if(HR==0.2){GHL=0.2*for_mat_B_st;}
  else{GHL=(2+8*for_mat_B_st/Thresh_denom)/100*for_mat_B_st;}             //GHL calculation differs for Sitka
   
                                                  
}

void model_parameters::get_FIGDATA(void)
{
  ofstream& evalout= *pad_evalout;
  FIGDATA.initialize();
  for (int i=1;i<=myrs;i++){
  for (int j=1;j<=1;j++){FIGDATA(i,j)=yminusthree(i);} //mature biomass (tons) (Figure 3)
  for (int j=2;j<=2;j++){FIGDATA(i,j)=yminustwo(i);} //mature biomass (tons) (Figure 3)
  for (int j=3;j<=3;j++){FIGDATA(i,j)=yminusone(i);} //mature biomass  (tons)(Figure 3)
  for (int j=4;j<=4;j++){FIGDATA(i,j)=tot_obs_egg(i+md_offset);}//observed total yearly egg deposition in trillions (Figure 1, Figure 2, Figure 3)
  for (int j=5;j<=5;j++){FIGDATA(i,j)=tot_est_egg(i);}//estimated total yearly egg deposition (Figure 1)
  for (int j=6;j<=6;j++){FIGDATA(i,j)=tot_mat_B(i);}// total mature biomass (metric tonnes) (Figure 2, Figure 3, Figure 6, Figure 14)
  for (int j=7;j<=7;j++){FIGDATA(i,j)=tcb(i+md_offset);}////observed catch biomass (tons)(Figure 2, Figure 3)
  for (int j=8;j<=8;j++){FIGDATA(i,j)=Threshold(i);}}//threshold vector (tons)(Figure 2, Figure 3)
  for (int i=1;i<=myrs-3;i++){
  for (int j=9;j<=9;j++){FIGDATA(i,j)=res_SR(i);}}//spawner-recruit residuals (Figure 4)
  for (int i=1;i<=myrs;i++){
  for (int j=10;j<=10;j++){FIGDATA(i,j)=res_tot_egg(i);}//egg deposition residuals (Figure 4)
  for (int j=11;j<=11;j++){FIGDATA(i,j)=init_age_3(i);}//initial age 3 (millions of fish) (Figure 5; Figure 10)
  for (int j=12;j<=12;j++){FIGDATA(i,j)=tot_sp_B(i);}// total spawning biomass (metric tons) (Figure 6; Figure 10, Figure 14)  
  for (int j=13;j<=13;j++){FIGDATA(i,j)=tot_sp_N(i);} // total spawning abundance [millions](Figure 6)
  for (int j=14;j<=14;j++){FIGDATA(i,j)=tot_mat_N(i);} // total mature abundance[millions] (Figure 6)
  for (int j=15;j<=15;j++){FIGDATA(i,j)=tot_post_N(i);}// total population [mature+ immature] - catch [millions](Figure 6)    
  for (int j=16;j<=16;j++){FIGDATA(i,j)=N(i);}// total abundance (mature + immature [millions] (Figure 6)
  for (int j=17;j<=22;j++){FIGDATA(i,j)=sel_naa_prop(i,j-16);}//proportion of N selected by gear (estimated) (Figure 8)
  for (int j=23;j<=28;j++){FIGDATA(i,j)=obs_c_comp(i+md_offset,j-22);}//observed catch compostion (Figure 8)
  for (int j=29;j<=34;j++){FIGDATA(i,j)=est_sp_comp(i,j-28);}//estimated spawner (cast net) age composition (Figure 9)
  for (int j=35;j<=40;j++){FIGDATA(i,j)=obs_sp_comp(i+md_offset,j-34);}}//observed spawner (cast net) age composition (Figure 9)
  for (int i=1;i<=myrs-3;i++){
  for (int j=41;j<=41;j++){FIGDATA(i,j)=SR(i);}}//(Figure 10)
  for (int i=1;i<=myrs;i++){
  for (int j=42;j<=42;j++){FIGDATA(i,j)=mod_yrs(i);}}//model years
}

void model_parameters::output_FIGDATA(void)
{
  ofstream& evalout= *pad_evalout;
 ofstream figdata("FIGDATA.dat");
 figdata<<"yminusthree yminustwo yminusone tot_obs_egg tot_est_egg tot_mat_B tcb Threshold res_SR res_tot_egg init_age_3 tot_sp_B tot_sp_N tot_mat_N tot_post_N N sel_naa_prop3 sel_naa_prop4 sel_naa_prop5 sel_naa_prop6 sel_naa_prop7 sel_naa_prop8 obs_c_comp3 obs_c_comp4 obs_c_comp5 obs_c_comp6 obs_c_comp7 obs_c_comp8 est_sp_comp3 est_sp_comp4 est_sp_comp5 est_sp_comp6 est_sp_comp7 est_sp_comp8 obs_sp_comp3 obs_sp_comp4 obs_sp_comp5 obs_sp_comp6 obs_sp_comp7 obs_sp_comp8 SR Year"<<endl;
 figdata<<FIGDATA<<endl;
}

void model_parameters::get_FIGDATAAGE(void)
{
  ofstream& evalout= *pad_evalout;
  FIGDATAAGE.initialize();
  for (int i=1;i<=nages;i++){
 //Mature biomass at age forecast
  for (int j=1;j<=1;j++){FIGDATAAGE(i,j)=for_mat_baa(i);}//forecasted mature biomass at age (metric tons) (Figure 11)
 //Mature numbers at age forecast (%)
  for (int j=2;j<=2;j++){FIGDATAAGE(i,j)=for_mat_prop(i);}//projected % mature #s at age (Figure 12)
 //forecast weight at age
  for (int j=3;j<=3;j++){FIGDATAAGE(i,j)=fw_a_a(i);}}//forecasted weight-at-age (Figure 13) 
  
}

void model_parameters::output_FIGDATAAGE(void)
{
  ofstream& evalout= *pad_evalout;
 ofstream figdataage("FIGDATAAGE.dat");
 figdataage<<"for_mat_baa for_mat_prop fw_a_a"<<endl;
 figdataage<<FIGDATAAGE<<endl;  
}

void model_parameters::compute_AICc(void)
{
  ofstream& evalout= *pad_evalout;
  n_tot_obs_egg.initialize();
  n_obs_c_comp.initialize();
  n_C.initialize();
  n_obs_sp_comp.initialize();
  n_S.initialize();
  n_R.initialize();
  n_d.initialize();
  n.initialize();
  sig_1.initialize();
  //sig_d.initialize();
  w_d.initialize();
  lnL_d.initialize();
  lnL.initialize();
  AIC.initialize();
  AICc.initialize();
  p.initialize();
 // n_My.initialize();
 // n_M.initialize();
  //Compute sample sizes
  
  //Egg Deposition
  for (int i=1;i<=myrs;i++)
    {
      n_tot_obs_egg(i)=1;
    }
  //SR function
  n_R=myrs-3;
  //Catch age comp
  for (int i=1;i<=myrs;i++)
    {
      for (int j=1;j<=nages;j++)
        {
          if (obs_c_comp(i+md_offset,j)<0)
            {
              n_obs_c_comp(i,j)=0;
            }
          else
            {
              n_obs_c_comp(i,j)=1;
            }
        }
    }
  n_C=sum(n_obs_c_comp);
  //Spawning age comp
  for (int i=1;i<=myrs;i++)
    {
      for (int j=1;j<=nages;j++)
        {
          if (obs_sp_comp(i+md_offset,1)<0)
            {
              n_obs_sp_comp(i,j)=0;
            }
          else
            {
              n_obs_sp_comp(i,j)=1;
            }
        }
     }
  n_S=sum(n_obs_sp_comp);
  //Mile-days
 // for (int i=1;i<=myrs;i++){
  //  if(mile_days(i+md_offset)>0){n_My(i)=1;}
   // else{n_My(i)=0;}}
 // n_M=sum(n_My);
  //Set up sample size vector
  n_d(1)=n_C;
  n_d(2)=n_S;
  n_d(3)=n_R;
 // n_d(3)=n_M;
  for (int i=1;i<=myrs;i++)
   {
     n_d(i+3)=n_tot_obs_egg(i);
   }
  n=sum(n_d);
  //Set up weighting vector
  w_d(1)=lC;
  w_d(2)=lS;
 // w_d(3)=lM;
  w_d(3)=lR;
  for (int i=1;i<=myrs;i++){
  w_d(i+3)=lE(i+md_offset);}
  //Set up sigma vector
  sig_1=f/n;
  //sig_d(1)=sig_1;
  //for (int i=1;i<=myrs+3;i++){
  //sig_d(i)=sig_1/w_d(i+1);}
  //Calculate log likelihood
  for (int i=1;i<=myrs+3;i++){ 
  if(w_d(i)>0){
   lnL_d(i)=-n_d(i)/2*(log(2*3.141593*sig_1/w_d(i))+1);}
  else{
   lnL_d(i)=0;}}
  lnL=sum(lnL_d);
  //Compute AICc
  p=initial_params::nvarcalc();
  AIC=-2*lnL+2*p;
  AICc=AIC+2*p*(p+1)/(n-p-1);
}

void model_parameters::get_sdreports(void)
{
  ofstream& evalout= *pad_evalout;
  for (int i=1;i<=myrs;i++){
  SpawnBio(i)=tot_sp_B(i);}
  SpawnBioFor=for_mat_B_st;
  AICcsd=AICc;
  GHLsd=GHL;
}

void model_parameters::get_report(void)
{
  ofstream& evalout= *pad_evalout;
       ofstream Report("Report.csv");
    int vsize = Year.size();
    Report << "MODEL RESULTS" <<endl;
 
    Report<<"Objective function value"<<","<<","<<","<<f<<endl;
    Report<<"AICc"<<","<<","<<","<<AICc<<endl;
    Report<<"  "<<endl;
    Report<<"Dataset Components:"<<endl;
    Report<<"Catch SSQ"<<","<<","<<","<<SSQC<<endl;
    Report<<"Spawning SSQ"<<","<<","<<","<<SSQSp<<endl;
    Report<<"Egg dep. SSQ"<<","<<","<<","<<WSSQE<<endl;
    Report<<"Ricker SSQ"<<","<<","<<","<<SSQR<<endl;
    Report<<"  "<<endl;
    Report<<"  "<<endl;
    Report<<"FORECAST INFORMATION:"<<endl;
    Report<<"Mature biomass forecast (short tons)"<<","<<","<<","<<SpawnBioFor<<endl;
    Report<<"Threshold (short tons)"<<","<<","<<","<<Thresh<<endl;
    Report<<"GHL (short tons)"<<","<<","<<","<<GHL<<endl;
    Report<<"Harvest Rate"<<","<<","<<","<<HR<<endl;
    Report<<"  "<<endl;
    Report<<"  "<<endl;
    Report<<","<<","<<","<<","<<"Age 3"<<","<<"Age 4"<<","<<"Age 5"<<","<<"Age 6"<<","<<"Age 7"<<","<<"Age 8+"<<endl;
    Report<<"Mature and immature numbers-at-age forecast (millions)"<<","<<","<<","<<","<<for_naa[1]<<","<<for_naa[2]<<","<<for_naa[3]<<","<<for_naa[4]<<","<<for_naa[5]<<","<<for_naa[6]<<endl;
    Report<<"Mature numbers-at-age forecast (millions)"<<","<<","<<","<<","<<for_mat_naa[1]<<","<<for_mat_naa[2]<<","<<for_mat_naa[3]<<","<<for_mat_naa[4]<<","<<for_mat_naa[5]<<","<<for_mat_naa[6]<<endl;
    Report<<"Proportion mature numbers-at-age forecast"<<","<<","<<","<<","<<for_mat_prop[1]<<","<<for_mat_prop[2]<<","<<for_mat_prop[3]<<","<<for_mat_prop[4]<<","<<for_mat_prop[5]<<","<<for_mat_prop[6]<<endl;
    Report<<"Mature biomass-at-age forecast (tons)"<<","<<","<<","<<","<<for_mat_baa[1]/0.90718<<","<<for_mat_baa[2]/0.90718<<","<<for_mat_baa[3]/0.90718<<","<<for_mat_baa[4]/0.90718<<","<<for_mat_baa[5]/0.90718<<","<<for_mat_baa[6]/0.90718<<endl;
    Report<<"Proportion mature biomass-at-age forecast"<<","<<","<<","<<","<<for_mat_b_prop[1]<<","<<for_mat_b_prop[2]<<","<<for_mat_b_prop[3]<<","<<for_mat_b_prop[4]<<","<<for_mat_b_prop[5]<<","<<for_mat_b_prop[6]<<endl;
    Report<<"Weight-at-age used in forecast (g)"<<","<<","<<","<<","<<fw_a_a[1]<<","<<fw_a_a[2]<<","<<fw_a_a[3]<<","<<fw_a_a[4]<<","<<fw_a_a[5]<<","<<fw_a_a[6]<<endl;
    Report<<"Maturity used for forecast"<<","<<","<<","<<","<<mat_for[1]<<","<<mat_for[2]<<","<<mat_for[3]<<","<<mat_for[4]<<","<<mat_for[5]<<","<<mat_for[6]<<endl;
    Report<<"Survival forecast"<<","<<","<<","<<","<<S_for<<","<<S_for<<","<<S_for<<","<<S_for<<","<<S_for<<","<<S_for<<endl;
    Report<<"  "<<endl;
    Report<<"  "<<endl;
   Report<<"Year"<<","<<"Mature biomass (tons)"<<","<<"Spawning biomass (tons)"<<","<<"Catch (tons)"<<","<<"Mature and immature age-3 abundance (recruitment in millions)"<<endl;
   for(int n; n<=vsize-2; n++)
   Report<<Year[n+mod_styr]<<","<<tot_mat_B[n+1]/0.90718<<","<<tot_sp_B[n+1]/0.90718<<","<<tcb[n+1+md_offset]<<","<<init_age_3[n+1]<<endl;
   Report<<"  "<<endl;
   Report<<"Mature and immature numbers-at-age (millions)"<<endl;
   Report<<"Year"<<","<<"Age 3"<<","<<"Age 4"<<","<<"Age 5"<<","<<"Age 6"<<","<<"Age 7"<<","<<"Age 8+"<<endl;
   for(int n; n<=vsize-2; n++)
   Report<<Year[n+mod_styr]<<","<<naa(n+1,1)<<","<<naa(n+1,2)<<","<<naa(n+1,3)<<","<<naa(n+1,4)<<","<<naa(n+1,5)<<","<<naa(n+1,6)<<endl;
   Report<<"  "<<endl;
   Report<<"Mature numbers-at-age (millions)"<<endl;
   Report<<"Year"<<","<<"Age 3"<<","<<"Age 4"<<","<<"Age 5"<<","<<"Age 6"<<","<<"Age 7"<<","<<"Age 8+"<<endl;
   for(int n; n<=vsize-2; n++)
   Report<<Year[n+mod_styr]<<","<<est_mat_naa(n+1,1)<<","<<est_mat_naa(n+1,2)<<","<<est_mat_naa(n+1,3)<<","<<est_mat_naa(n+1,4)<<","<<est_mat_naa(n+1,5)<<","<<est_mat_naa(n+1,6)<<endl;
   Report<<"  "<<endl;
   Report<<"Spawning numbers-at-age (millions)"<<endl;
   Report<<"Year"<<","<<"Age 3"<<","<<"Age 4"<<","<<"Age 5"<<","<<"Age 6"<<","<<"Age 7"<<","<<"Age 8+"<<endl;
   for(int n; n<=vsize-2; n++)
   Report<<Year[n+mod_styr]<<","<<est_sp_naa(n+1,1)<<","<<est_sp_naa(n+1,2)<<","<<est_sp_naa(n+1,3)<<","<<est_sp_naa(n+1,4)<<","<<est_sp_naa(n+1,5)<<","<<est_sp_naa(n+1,6)<<endl;
   Report<<"  "<<endl;
   Report<<"Estimated catch-at-age (millions)"<<endl;
   Report<<"Year"<<","<<"Age 3"<<","<<"Age 4"<<","<<"Age 5"<<","<<"Age 6"<<","<<"Age 7"<<","<<"Age 8+"<<endl;
   for(int n; n<=vsize-2; n++)
   Report<<Year[n+mod_styr]<<","<<est_c_naa(n+1,1)<<","<<est_c_naa(n+1,2)<<","<<est_c_naa(n+1,3)<<","<<est_c_naa(n+1,4)<<","<<est_c_naa(n+1,5)<<","<<est_c_naa(n+1,6)<<endl;
   Report<<"  "<<endl;
   Report<<"Estimated proportion of each age class that is caught (e.g. number of age-3 fish caught divided by the number of age-3 mature at age)"<<endl;
   Report<<"Year"<<","<<"Age 3"<<","<<"Age 4"<<","<<"Age 5"<<","<<"Age 6"<<","<<"Age 7"<<","<<"Age 8+"<<endl;
   for(int n; n<=vsize-2; n++)
   Report<<Year[n+mod_styr]<<","<<est_c_naa(n+1,1)/est_sp_naa(n+1,1)<<","<<est_c_naa(n+1,2)/est_sp_naa(n+1,2)<<","<<est_c_naa(n+1,3)/est_sp_naa(n+1,3)<<","<<est_c_naa(n+1,4)/est_sp_naa(n+1,4)<<","<<est_c_naa(n+1,5)/est_sp_naa(n+1,5)<<","<<est_c_naa(n+1,6)/est_sp_naa(n+1,6)<<endl;
   Report<<"  "<<endl;
   Report<<"Mature biomass-at-age"<<endl;
   Report<<"Year"<<","<<"Age 3"<<","<<"Age 4"<<","<<"Age 5"<<","<<"Age 6"<<","<<"Age 7"<<","<<"Age 8+"<<endl;
   for(int n; n<=vsize-2; n++)
   Report<<Year[n+mod_styr]<<","<<est_mat_baa(n+1,1)/0.90718<<","<<est_mat_baa(n+1,2)/0.90718<<","<<est_mat_baa(n+1,3)/0.90718<<","<<est_mat_baa(n+1,4)/0.90718<<","<<est_mat_baa(n+1,5)/0.90718<<","<<est_mat_baa(n+1,6)/0.90718<<endl;
   Report<<"  "<<endl;
   Report<<"Spawning biomass-at-age (tons)"<<endl;
   Report<<"Year"<<","<<"Age 3"<<","<<"Age 4"<<","<<"Age 5"<<","<<"Age 6"<<","<<"Age 7"<<","<<"Age 8+"<<endl;
   for(int n; n<=vsize-2; n++)
   Report<<Year[n+mod_styr]<<","<<est_sp_baa(n+1,1)/0.90718<<","<<est_sp_baa(n+1,2)/0.90718<<","<<est_sp_baa(n+1,3)/0.90718<<","<<est_sp_baa(n+1,4)/0.90718<<","<<est_sp_baa(n+1,5)/0.90718<<","<<est_sp_baa(n+1,6)/0.90718<<endl;
   Report<<"  "<<endl;
   Report<<"Survival"<<endl;
   Report<<"Year"<<","<<"Age 3"<<","<<"Age 4"<<","<<"Age 5"<<","<<"Age 6"<<","<<"Age 7"<<","<<"Age 8+"<<endl;
   for(int n; n<=vsize-2; n++)
   Report<<Year[n+mod_styr]<<","<<Sur(n+1,1)<<","<<Sur(n+1,2)<<","<<Sur(n+1,3)<<","<<Sur(n+1,4)<<","<<Sur(n+1,5)<<","<<Sur(n+1,6)<<endl;
   Report<<"  "<<endl;
   Report<<"Maturity"<<endl;
   Report<<"Year"<<","<<"Age 3"<<","<<"Age 4"<<","<<"Age 5"<<","<<"Age 6"<<","<<"Age 7"<<","<<"Age 8+"<<endl;
   for(int n; n<=vsize-2; n++)
   Report<<Year[n+mod_styr]<<","<<Mat(n+1,1)<<","<<Mat(n+1,2)<<","<<Mat(n+1,3)<<","<<Mat(n+1,4)<<","<<Mat(n+1,5)<<","<<Mat(n+1,6)<<endl;
   Report<<"  "<<endl;
   Report<<"Gear Selectivity"<<endl;
   Report<<"Year"<<","<<"Age 3"<<","<<"Age 4"<<","<<"Age 5"<<","<<"Age 6"<<","<<"Age 7"<<","<<"Age 8+"<<endl;
   for(int n; n<=vsize-2; n++)
   Report<<Year[n+mod_styr]<<","<<GS_Sc(n+1,1)<<","<<GS_Sc(n+1,2)<<","<<GS_Sc(n+1,3)<<","<<GS_Sc(n+1,4)<<","<<GS_Sc(n+1,5)<<","<<GS_Sc(n+1,6)<<endl;
   Report<<"  "<<endl;
   Report<<"Observed spawning age-composition (cast net)"<<endl;
   Report<<"Year"<<","<<"Age 3"<<","<<"Age 4"<<","<<"Age 5"<<","<<"Age 6"<<","<<"Age 7"<<","<<"Age 8+"<<endl;
   for(int n; n<=vsize-2; n++)
   Report<<Year[n+mod_styr]<<","<<obs_sp_comp(n+10,1)<<","<<obs_sp_comp(n+10,2)<<","<<obs_sp_comp(n+10,3)<<","<<obs_sp_comp(n+10,4)<<","<<obs_sp_comp(n+10,5)<<","<<obs_sp_comp(n+10,6)<<endl;
   Report<<"  "<<endl;
   Report<<"Observed commercial seine age-composition (spring seine)"<<endl;
   Report<<"Year"<<","<<"Age 3"<<","<<"Age 4"<<","<<"Age 5"<<","<<"Age 6"<<","<<"Age 7"<<","<<"Age 8+"<<endl;
   for(int n; n<=vsize-2; n++)
   Report<<Year[n+mod_styr]<<","<<obs_c_comp(n+10,1)<<","<<obs_c_comp(n+10,2)<<","<<obs_c_comp(n+10,3)<<","<<obs_c_comp(n+10,4)<<","<<obs_c_comp(n+10,5)<<","<<obs_c_comp(n+10,6)<<endl;
   Report<<"  "<<endl;
   Report.close();
 
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
   //The reports below have to included in the output or the R figures will not work properly
  //Note: you can't name the same variable two different names or else the reptoRlist function in program R will give an error
  report<<"Scaled_Gear_Selectivity"<<endl;//needs to be in report for figures (Figure 7)
  report<<GS_Sc<<endl;
  report<<"Maturity"<<endl;//needs to be in report for figures (Figure 7)
  report<<Mat<<endl;
  report<<"Survival"<<endl;//needs to be in report for figures (Figure 7)
  report<<Sur<<endl;
  report<<"GHL"<<endl;//needs to be in report for figures (Figure 14)
  report<<GHL<<endl;
  report<<"yminusthreeFOR"<<endl;//needs to be in report for figures-make sure in tons (Figure 3)
  report<<yminusthreeFOR<<endl;
  report<<"yminustwoFOR"<<endl;//needs to be in report for figures-make sure in tons (Figure 3)
  report<<yminustwoFOR<<endl;
  report<<"yminusoneFOR"<<endl;//needs to be in report for figures-make sure in tons (Figure 3)
  report<<yminusoneFOR<<endl;
  report<<"for_mat_B_st"<<endl;//needs to be in report for figures-make sure in tons (Figure 2, Figure 3, Figure 14)
  report<<for_mat_B_st<<endl;
  report<<"Thresh"<<endl;//current threshold-make sure in tons//needs to be in report for figures (Figure 14)
  report<<Thresh<<endl;
  report<<"res_sp_comp"<<endl;//needs to be in report for figures (Figure 15)
  report<<res_sp_comp<<endl;
  report<<"est_sp_comp"<<endl;//needs to be in report for figures (Figure 16)
  report<<est_sp_comp<<endl;
  report<<"res_c_comp"<<endl;//needs to be in report for figures (Figure 17)
  report<<res_c_comp<<endl;
  report<<""<<endl;
  report<<"MODEL_RESULTS:"<<endl;
  report<<"Objective_function_value"<<endl;
  report<<f<<endl;
  report<<"Dataset_components:"<<endl;
  report<<"Catch SSQ"<<endl;
  report<<SSQC<<endl;
  report<<"Spawning_SSQ"<<endl;
  report<<SSQSp<<endl;
  report<<"Egg_dep_SSQ"<<endl;
  report<<WSSQE<<endl;
  report<<"Mile_Days SSQ"<<endl;;
  report<<WSSQM<<endl;
  report<<"Ricker_SSQ"<<endl;
  report<<SSQR<<endl;
  report<<"FORECAST INFORMATION:"<<endl;
  report<<"Harvest Rate"<<endl;
  report<<HR<<endl;
  report<<"Numbers-at-age forecast"<<endl;
  report<<for_naa<<endl;
  report<<"Mature biomass-at-age forecast"<<endl;
  report<<for_mat_baa<<endl;
  report<<"Mature biomass at age"<<endl;
  report<<est_mat_baa<<endl;
  report<<"spawning biomass at age"<<endl;
  report<<est_sp_baa<<endl;
}

void model_parameters::set_runtime(void)
{
  dvector temp1("{5000 5000 5000 5000}");
  maximum_function_evaluations.allocate(temp1.indexmin(),temp1.indexmax());
  maximum_function_evaluations=temp1;
  dvector temp("{0.0001}");
  convergence_criteria.allocate(temp.indexmin(),temp.indexmax());
  convergence_criteria=temp;
}

void model_parameters::final_calcs()
{
	//  Print run time statistics to the screen.
	time(&finish);
	elapsed_time=difftime(finish,start);
	hour=long(elapsed_time)/3600;
	minute=long(elapsed_time)%3600/60;
	second=(long(elapsed_time)%3600)%60;
	cout<<endl<<endl<<"*******************************************"<<endl;
	cout<<"--Start time: "<<ctime(&start)<<endl;
	cout<<"--Finish time: "<<ctime(&finish)<<endl;
	cout<<"--Runtime: ";
	cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
        cout<<""<<endl;
	cout<<"--Objective function value: "<<f<<endl;
        cout<<""<<endl;
	cout<<"--AIC: "<<AICc<<endl;
        cout<<""<<endl;
	cout<<"--Maximum gradient component: "<<objective_function_value::gmax<<endl;
        cout<<""<<endl;
	cout<<"*******************************************"<<endl;
        cout<<""<<endl;
        cout<< "O frabjous day!"<<endl;
        cout<< "The sheep frolic!"<<endl;
        cout<<""<<endl;
        cout<<"        ...moo..."<<endl;
        cout<<"            | "<<endl;
        cout<<"            | "<<endl;
        cout<<"            | "<<endl;
        cout<<"             _.%%%%%%%%%%%%%             " <<endl;
        cout<<"            //-_%%%%%%%%%%%%%            " <<endl;
        cout<<"           (_ %\\%%%%%%%%%%%%%%~             "<<endl;
        cout<<"               %%%%%%%%%%%%%%             "<<endl;
        cout<<"                 %%%%%*%%%%              "<<endl;
        cout<<"            ,,,,,,||,,,,||,,,,,         "<<endl;
        cout<<""<<endl;
		cout<<"            ▄▄▄▄▄▄▄▄▄▄▄▄▄ "<<endl;
		cout<<"         ▄▀▀═════════════▀▀▄ "<<endl;
		cout<<"        █═══════════════════█ "<<endl;
		cout<<"       █═════════════════════█ "<<endl;
		cout<<"      █═══▄▄▄▄▄▄▄═══▄▄▄▄▄▄▄═══█ "<<endl;
		cout<<"     █═══█████████═█████████═══█ "<<endl;
		cout<<"     █══██▀    ▀█████▀    ▀██══█ "<<endl;
		cout<<"    ██████  █▀█  ███  █▀█  ██████ "<<endl;
		cout<<"    ██████  ▀▀▀  ███  ▀▀▀  ██████ "<<endl;
		cout<<"     █══▀█▄    ▄██ ██▄    ▄█▀══█ "<<endl;
		cout<<"     █════▀█████▀   ▀█████▀════█ "<<endl;
		cout<<"     █═════════════════════════█ "<<endl;
		cout<<"     █═════════════════════════█ "<<endl;
		cout<<"     █═══════▀▄▄▄▄▄▄▄▄▄▀═══════█ "<<endl;
		cout<<"     █═════════════════════════█ "<<endl;
		cout<<"    ▐▓▓▌═════════════════════▐▓▓▌ "<<endl;
		cout<<"    ▐▐▓▓▌▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▐▓▓▌▌ "<<endl;
		cout<<"    █══▐▓▄▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▄▓▌══█ "<<endl;
		cout<<"   █══▌═▐▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▌═▐══█ "<<endl;
		cout<<"   █══█═▐▓▓▓▓▓▓▄▄▄▄▄▄▄▓▓▓▓▓▓▌═█══█ "<<endl;
		cout<<"   █══█═▐▓▓▓▓▓▓▐██▀██▌▓▓▓▓▓▓▌═█══█ "<<endl;
		cout<<"   █══█═▐▓▓▓▓▓▓▓▀▀▀▀▀▓▓▓▓▓▓▓▌═█══█ "<<endl;
		cout<<"   █══█▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓█══█ "<<endl;
		cout<<"  ▄█══█▐▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▌█══█▄ "<<endl;
		cout<<"  █████▐▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓█ █████ "<<endl;
		cout<<"  ██████▐▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▌███████ "<<endl;
		cout<<"   ▀█▀█  ▐▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▌   █▀█▀ "<<endl;
		cout<<"          ▐▓▓▓▓▓▓▌▐▓▓▓▓▓▓▌ "<<endl;
		cout<<"           ▐▓▓▓▓▌  ▐▓▓▓▓▌ "<<endl;
		cout<<"          ▄████▀    ▀████▄ "<<endl;
		cout<<"          ▀▀▀▀        ▀▀▀▀ "<<endl;
	
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{
  delete pad_evalout;
  pad_evalout = NULL;
}

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
  arrmblsize=5000000;
  gradient_structure::set_MAX_NVAR_OFFSET(5000);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(800000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(800000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);
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

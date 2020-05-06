//We need to translate this code to TMB. 
//version 1.3 is denoted as v13; 
//Revision
//(1) in implementing random walk for fishing mortality
//(2) in treatment of cpue data as the stochastic index data, where It = q*Bt*exp(error_t)
//(3) considering age-specific M instead of just M  
//
//data file name: mackereldata.DAT
//Author: Saang-Yoon Hyun, Kyuhan Kim, Jinwoo Gim, Minkyu Bahk and Doyul Kim as of April 21; 

DATA_SECTION
  init_int nages; 		// 6
  init_int nlengths;  	// 42 classes
  init_int nyrs; 		// 22 years

  init_matrix yieldcpue(1,nyrs,1,3);
  vector years(1,nyrs); 
  !!years=column(yieldcpue,1); 
  
  init_int indexMinyrLD;  //2000 ==> index of 6; 
 
  int nyrsLD;       //nyears of the length data 
  !!nyrsLD=nyrs-5;  //i.e., 18 years of 2000-2017; //<= as of 4/27
   
  init_vector x(1,nlengths); //discrete lengths;
  vector L(1,nlengths);      //discrete lengths;
  !!L=x; 
  
  init_matrix lengthfq(1,nyrsLD,1,nlengths);  //length frequency by year and length class                                                              
                             
  vector yield(1,nyrs); // in MT
  !!yield=column(yieldcpue,2);
  vector CPUE(1,nyrs);  // in MT/haul
  !!CPUE=column(yieldcpue,3);
  vector Effort(1,nyrs); 
  !!Effort=elem_div(yield,CPUE); // vector(1,23) in hauls
  
  init_vector par_logistic(1,2); //For the maturation rate by length class;
  number b0;
  number b1;
  !!b0=par_logistic(1); 
  !!b1=par_logistic(2); 
  vector maturation(1,nlengths); 
  
  init_vector musig2_r(1,2);
  number mu_r;
  number sig2_r;
  !!mu_r=musig2_r(1);
  !!sig2_r=musig2_r(2);
  
  init_vector abWL(1,2); 
  number aWL; 
  number bWL;  
  !!aWL=abWL(1);
  !!bWL=abWL(2);
  
  vector Wt(1,nlengths);
  
  init_number powerML;  //power (=beta*b1) in M = par*L^power
                        //ad hoc value = 3.425 * 0.305, where 3.425 from JW, and 0.305 from Lorenzen (1996);
                         
  //init_vector M(1,nlengths); // length-specific M 
  init_number lambda1; // lambda in the objective function
  init_number lambda2; // lambda in the objective function
  
  init_number logLinf;
  number Linf;
  !!Linf=mfexp(logLinf);
  
  init_int choiceSR;  //choice of a stock-recruitment model;  
  
  int r;
  !!r=1; //recruitment is defined as the pop size at age 1; 
  int ncohorts;
  !!ncohorts=nyrs;
  
  ivector SamSize(1,nyrsLD);   //annual sample size for the length data;    
  //init_number sig2_yield;
  
  !!cout<<"choiceSR: "<<choiceSR<<endl; 

  
PARAMETER_SECTION
  init_number logF1(1);                        //log(F in the first time, i.e., first year);
  init_bounded_vector F_devs(1,nyrs-1,lowFdev,upFedv,1);         //random walk deviations for annual fishing mortality in years after 1st year; 
  init_bounded_vector logRec(1,indexMinyrLD,10.0,28.0,2); //age-1 N in 1996 ~ 2000; 
  init_bounded_number logkappa(-3.0,-0.80,3);
  init_bounded_number logL50(3.0,4.0,4);     //selectivity 
  init_bounded_number gamma(0.05,1.5,5);     //selectivity
  init_bounded_number sigmaL(0.0,3.0,6);     //the uncertainty in the L_{a+1} equation
  init_number log_parML(7);     //the parameter in M = par * L^(3.425*0.305);
  
  //Under B-H, Recruits_{t+1} = Spawners_t/(a+b*Spawners_t); 
  //Under Ricker, Recruits_{t+1} = a*Spawners_t*exp(-b*Spawners_t); 
  init_bounded_number logaSR(0.0,10.0,2);    //parameter a in the stock-recruitment model 
  init_bounded_number logbSR(-30.0,0.0,2);     //parameter b in the stock-recruitment model 
                                            
  vector M(1,nyrs);  
  vector F_t(1,nyrs);      //instantaneous fishing mortality by year 
                          //previous notation: F_yr(1,nyrs); 
  vector logFt(1,nyrs);   //log(Ft);
  
  number q;
  vector temelems(indexMinyrLD,nyrs); 
  
  vector Recruits(1,nyrs+1); 
  number kappa; 
  number L50; 
  number aSR; 
  number bSR; 
  
  vector Sel(1,nlengths); 
  number Rho; 
  number kkk; //for the cumulative purpose
  matrix f(1,nages,1,nlengths);  //length frequency as pmf                          
  3darray pp(1,nages,1,nlengths,1,nlengths); //pp(1 To Ages,1 To L,1 To x); 
                                             // L -- (growth) --> x
  vector Mu(1,nlengths); //differ by length class
  vector SS(1,nages);    //assumed to be constant over length classes
 
 
  matrix F_tx(1,nyrs,1,nlengths);  //F_{year,xind}; 
  matrix Z(1,nyrs,1,nlengths);
  matrix ExpZ(1,nyrs,1,nlengths);  //Exp(-Z); 
  
  3darray NL(1,nages,1,ncohorts+1,1,nlengths); 
  matrix N(1,nages,1,nyrs+1);  //the same index as VB         
 
  3darray SpawnersL(1,nages,1,ncohorts+1,1,nlengths);
  matrix Spawners(1,nages,1,nyrs+1);
  
  number SumP;
  vector p(1,nlengths); 
  vector p_plus(1,nlengths); //for the last age class
  
  number CNum;
  number CWt;
  vector TCatch(1,nyrs);
  matrix Catch(1,nyrs,1,nlengths);  
  vector Yieldhat(1,nyrs); 
  vector Pop(1,nyrs); 
  3darray ENx(1,nages,1,ncohorts,1,nlengths); 
  vector EN(1,ncohorts); 
  vector B(1,ncohorts);
  vector EB(1,ncohorts); 
  
  matrix hatlengthfq(1,nyrsLD,1,nlengths);    //predicted length-frequency
 
  number logmult; //log(multinomial); 
  vector elem_obj2(indexMinyrLD,nyrs);  //elements in part 2 of the objective function;
  
  number lognormal;  //log(normal); 
  number maxloglike; // for calculation of AIC
  number aic; 
  //number sig2_yield;
  
  objective_function_value obj;
  
PRELIMINARY_CALCS_SECTION
  //The maturation rate by length class;
  for(int i=1;i<=nlengths;i++)
      maturation(i)=mfexp(b0+b1*x(i))/(1+mfexp(b0+b1*x(i)));
  
  //Calculate the weight relationship
  Wt=aWL*pow(x,bWL)/1000; //the division of 1000 is to convert gram to kg 
  
  //annual sample size for the length data
  for(int i=1;i<=nyrsLD;i++) 
     SamSize(i)=sum(lengthfq(i));
 
  cout<<"maturation: "<<maturation<<endl;  
  cout<<"lengthfq(1,nyrsLD,1,nlengths): "<<endl; 
  cout<<lengthfq<<endl;
  cout<<"SamSize: "<<SamSize<<endl; 
  
PROCEDURE_SECTION  
  kappa=mfexp(logkappa); 
  L50=mfexp(logL50);
  parML=mfexp(log_parML);   
  aSR=mfexp(logaSR); 
  bSR=mfexp(logbSR);   
  
  //instantaneous natural mortality by length class
  for(int xind=1; xind<=nlengths; xind++) 
      M(xind)=parML*(pow(x(xind),powerML); 
  
  //instantaneous fishing mortality by time, i.e., year;
  log_F(1)=log_F1;                 
  F_t(1)=mfexp(log_F(1));
  for(int t=2; t<=nyrs; t++)  {
     log_F(t)=log_F(t-1)+F_devs(t-1); 
     F_t(t)=mfexp(log_F(t)); 
  };    
  
  //
  for(int i=1;i<=indexMinyrLD;i++)
      Recruits(i)=mfexp(logRec(i));
         
  obj=0.0;
  
  //Calculate the selectivity for each length class; 
  Sel=1.0/(1+mfexp(-1.0*gamma*(x-L50))); 
 
  //Calculate fishing mortality for each length using F = q*Effort*Sel;
  //Sum_over x_F_{yr,x}  =/= F_{yr}; "=/=" "is not equal to"
  for(int t=1;t<=nyrs;t++) 
    for(int xind=1;xind<=nlengths;xind++)  { 
       F_tx(t,xind)=F_t(t)*Sel(xind); //F_tx; //FM in the VB code;
       Z(t,xind)=M(xind)+F_tx(t,xind); 
       ExpZ(t,xind)=mfexp(-1.0*Z(t,xind)); 
    };       
    
  //Calculating the length frequency dsn of the recruits;
  //recruitment is at one year of age;
  //f(1,x) is the same for all cohorts so this can be in the initial calculations;
  Rho=mfexp(-1.0*kappa); //Brody coefficient
  SS(1)=sig2_r; 
  kkk=0.0;

  for(int xind=1;xind<=nlengths;xind++) {
     f(1,xind)=mfexp(-1.0*square(x(xind)-mu_r)/(2.0*SS(1))); //Eq. 1 in Quinn et al. (1998) in the AK Sea Grant College Program
     kkk=kkk+f(1,xind); //similar to the cumulative normal prob.
  };
   
  for(int xind=1;xind<=nlengths;xind++) {
     f(1,xind)=f(1,xind)/kkk;  //normalize; //Eq. 1 in Quinn et al (1998) in the AK Sea Grant College Program
     Mu(xind)=Linf-(Linf-x(xind))*Rho; 
  };
 
  for(int a=2;a<=nages;a++)  {
      SS(a)=square(sigmaL)*(1.0-pow(Rho,(2.0*a-2.0*r)))/(1.0-square(Rho))+(pow(Rho,(2.0*a-2.0*r)))*sig2_r; 
             //this SS is for Shrimp; //see Cohen and Fishman (1980)
  };
    
  for(int a=1;a<=nages;a++) {
     for(int xind=1;xind<=nlengths;xind++)  {
        kkk=0.0; 
        for(int Lind=1;Lind<=nlengths;Lind++)  {
            pp(a,Lind,xind)=mfexp(-1.0/(2.0*SS(a))*square(L(Lind)-Mu(xind)));
                                     //Eq. 13 in Quinn et al. (1998); //It is expressed as f(L|x).
            kkk=kkk+pp(a,Lind,xind); //similar to the cumulative normal prob.
        };
        for(int Lind=1;Lind<=nlengths;Lind++)
            pp(a,Lind,xind)=pp(a,Lind,xind)/kkk;  //Eq. 13 in Quinn et al. (1998); //It is expressed as f(L|x).
     };       
  };
     
  //cout<<"pp: "<<endl; 
  //cout<<pp<<endl; 
      
  //Start of cohort loop
  int m;     //m: time index
  for(int c=1;c<=ncohorts+1;c++)  {   //c: cohort index
     int a=1; 
     m=c;          //m: time index
     
     N(a,m)=0.0;  
     Spawners(a,m)=0.0; 
     if(c<=indexMinyrLD) {
        for(int xind=1;xind<=nlengths;xind++)  {
           NL(a,m,xind)=Recruits(m)*f(a,xind);   //<***************************************
           N(a,m)=N(a,m)+NL(a,m,xind);   
        
           SpawnersL(a,m,xind)=NL(a,m,xind)*maturation(xind); 
           Spawners(a,m)=Spawners(a,m)+SpawnersL(a,m,xind); 
        };
     }
     else   {
	    if(choiceSR==1) 
	       N(a,m)=sum(column(Spawners,m-1))/(aSR+bSR*sum(column(Spawners,m-1)));  //B-H   
	    else if(choiceSR==2) 
	       N(a,m)=aSR*sum(column(Spawners,m-1))*mfexp(-1.0*bSR*sum(column(Spawners,m-1)));  //Ricker   
	     
	    for(int xind=1;xind<=nlengths;xind++)  {
	       NL(a,m,xind) = N(a,m)*f(a,xind);
	    
     	   SpawnersL(a,m,xind)=NL(a,m,xind)*maturation(xind); 
           Spawners(a,m)=Spawners(a,m)+SpawnersL(a,m,xind);   
        }; 
     };    
 
     for(int a=2;a<=nages;a++)   {  //note age starts at 2, not 1.
        m=a+c-1;    //m: cohort index

        if(m<=nyrs+1)   {
	       SumP=0.0;  
           for(int Lind=1;Lind<=nlengths;Lind++)  {
               p(Lind)=0;    
               for(int xind=1;xind<=nlengths;xind++)
                   p(Lind)=p(Lind)+f(a-1,xind)*ExpZ(m-1,xind)*pp(a,Lind,xind);  //Eq. 14 in the AK Sea Grant's paper 
                                                                                //written by Quinn, Turnbull, and Fu
               SumP=SumP+p(Lind);
           };  
          
           if(a!=nages) {
              N(a,m)=0.0;  
              Spawners(a,m)=0.0; 
              for(int Lind=1;Lind<=nlengths;Lind++)  {
                   f(a,Lind)=p(Lind)/SumP;  //normalize;
                   
                   NL(a,m,Lind)=N(a-1,m-1)*p(Lind);
                              //cf. NL(a,m,Lind)=N(m-1,a-1)*f(a-1,Lind);
                   N(a,m)=N(a,m)+NL(a,m,Lind);
                   
                   SpawnersL(a,m,Lind)=NL(a,m,Lind)*maturation(Lind); 
                   Spawners(a,m)=Spawners(a,m)+SpawnersL(a,m,Lind); 
              };
           }
           else if(a==nages) {
               for(int Lind=1;Lind<=nlengths;Lind++)  {
                   p_plus(Lind)=0.0;    
                   for(int xind=1;xind<=nlengths;xind++)
                        p_plus(Lind)=p_plus(Lind)+f(a,xind)*ExpZ(m-1,xind)*pp(a,Lind,xind);  // f(a,xind) is f(6,xind)
               };   
	           
	           N(a,m)=0.0;  
	           Spawners(a,m)=0.0; 
               for(int Lind=1;Lind<=nlengths;Lind++)  {
                   NL(a,m,Lind)=N(a-1,m-1)*p(Lind)+N(a,m-1)*p_plus(Lind);   //new p_plus
                              //NL(a,m,Lind)=N(a-1,m-1)*f(a-1,Lind)+N(a,m-1)*f(a,Lind); 
                   N(a,m)=N(a,m)+NL(a,m,Lind);
                    
                   SpawnersL(a,m,Lind)=NL(a,m,Lind)*maturation(Lind); 
                   Spawners(a,m)=Spawners(a,m)+SpawnersL(a,m,Lind); 
               };     
               f(a)=NL(a,m)/N(a,m);  //normalize; f(6) is a vector    
	       };     
       };
    };
  };  //cohort index ends here
    
  for(int i=indexMinyrLD+1;i<=nyrs+1;i++)
        Recruits(i)=N(1,i);
        
  //cout<<"NL(a,m,xind): "<<endl;
  //cout<<NL<<endl;              
  //cout<<"N(a,m): "<<endl;      
  //cout<<N<<endl;     
  //cout<<"transpose of N(a,m) ==> N(m,a): "<<endl; 
  //cout<<trans(N)<<endl; 
  //cout<<"f(a,Lind): "<<endl; 
  //cout<<f<<endl;          
  
  cout<<"N(a,m)"<<endl;
  cout<<N<<endl; 
  cout<<"Spawners(a,m)"<<endl;   
  cout<<Spawners<<endl; 
    
  cout<<"aSR, bSR: "<<aSR<<" "<<bSR<<endl; 
    
  //================  //note m starts at indexMinyrLD;
  for(int m=indexMinyrLD;m<=nyrs;m++) {  //note m starts at indexMinyrLD;
	  TCatch(m)=0.0;
	  Catch(m)=0.0;       //Catch(1,nyrs,1,nlengths);  
	  Yieldhat(m)=0.0; 
	  Pop(m)=0.0;      
	  EN(m)=0.0;       
	  B(m)=0.0;        
	  EB(m)=0.0;      
	   
	  for(int a=1;a<=nages;a++)   {
	      for(int xind=1;xind<=nlengths;xind++) {
		      CNum=NL(a,m,xind)*(F_tx(m,xind)/Z(m,xind))*(1-ExpZ(m,xind));
		      CWt=CNum*Wt(xind);  //in kg
		      Catch(m,xind)=Catch(m,xind)+CNum;
		      TCatch(m)=TCatch(m)+CNum; 
		      Yieldhat(m)=Yieldhat(m)+CWt; 
		     
		      Pop(m)=Pop(m)+NL(a,m,xind);   //population; 
		      ENx(a,m,xind)=NL(a,m,xind)*Sel(xind); //Exploitable population; 
		      EN(m)=EN(m)+ENx(a,m,xind);            //Exploitable population; 
		      B(m)=B(m)+NL(a,m,xind)*Wt(xind);      //Biomass;  //in kg
		      EB(m)=EB(m)+ENx(a,m,xind)*Wt(xind);   //expoitable biomass; 
	      };
      };  
  };  //m ends here;  
 
  cout<<"Yieldhat: "<<Yieldhat<<endl; 
  
  //The expected length-frequency 
  for(int m=indexMinyrLD;m<=nyrs;m++)  //note m starts at indexMinyrLD, 5;
     for(int xind=1;xind<=nlengths;xind++)    
	    hatlengthfq(m-5,xind)=(Catch(m,xind)/sum(Catch(m)))*SamSize(m-5);  //<== as of 4/27
	 
  // cout<<"hatlengthfq(nyrsLD, lengths): "<<endl;    
  // cout<<hatlengthfq<<endl; 
  
  //part 1 of the objective funcion: multinomial
  logmult=0.0;
  for(int i=indexMinyrLD;i<=nyrs;i++) {  //indexMinyrLD = 6;      
     logmult+=gammln(SamSize(i-5)+1);    //<== as of 4/27
     for(int xind=1;xind<=nlengths;xind++) {     
        logmult+=-1.0*gammln((lengthfq(i-5,xind)+1))+lengthfq((i-5),xind)*log(Catch(i,xind)/sum(Catch(i)));      //<== as of 4/27        
     };
  };     
  
  obj+=lambda1*(-1.0*logmult); //lamda1*(the negative multinomial likelihood);
  
  maxloglike=0.0; 
  maxloglike+=logmult; 
  
  //part 2 of the objective function: lognormal for yield data
  for(int m=indexMinyrLD;m<=nyrs;m++) 
     elem_obj2(m)=square( log(yield(m))-log(Yieldhat(m)/1000)+sig2_yield/(2.0*square(Yieldhat(m)/1000))  );  //in MT 
                                 //the above is based on the delta method
  obj+=lambda2*(0.5*(nyrs-indexMinyrLD+1)*log(sig2_yield)+sum(elem_obj2)/(2.0*sig2_yield));  //lamda2*(the negative normal loglikelihood);

  cout<<"elem_obj2: "<<endl; 
  cout<<elem_obj2<<endl; 
  
  lognormal1=0.0; 
  lognormal1=(-0.5*(nyrs-indexMinyrLD+1)*log(2*M_PI)-0.5*(nyrs-indexMinyrLD+1)*log(sig2_yield)-sum(elem_obj2)/(2.0*sig2_yield));
  
  maxloglike+=lognormal1;  
  
  
  //part 3 of the objective function: lognormal for cpue data
  //log(CPUE) ~ normal (log(q*Bt), sig2_cpue); 
  for(int m=indexMinyrLD;m<=nyrs;m++) 
      temelems(m) = log(CPUE(m))-log(B(m));
  
  q=mfexp(sum(temelmens)/(nyrs-indexMindyrLD+1));
          
  for(int m=indexMinyrLD;m<=nyrs;m++) 
      elem_obj3(m)=square(log(CPUE(m))-log(q)-log(B(m)));      
    
  sig2_logCPUE=sum(elem_obj3)/(nyrs-indexMinyrLD+1); 
  
  obj+=lambda3*(0.5*(nyrs-indexMinyrLD+1)*log(sig2_logCPUE)+sum(elem_obj3)/(2.0*sig2_logCPUE));  //lamda2*(the negative normal loglikelihood);

  //cout<<"elem_obj3: "<<endl; 
  //cout<<elem_obj3<<endl; 
  
  lognormal2=0.0; 
  lognormal2=(-0.5*(nyrs-indexMinyrLD+1)*log(2*M_PI)-0.5*(nyrs-indexMinyrLD+1)*log(sig2_logCPUE)+sum(elem_obj2)/(2.0*sig2_logCPUE));  
 
  maxloglike+=lognormal2;  

  /////////////////////////
  aic=-2.0*maxloglike+2.0*(nyrs+5);   
  
 
  
REPORT_SECTION
  report<<"Yr Recruits Fyr Yield Yieldhat"<<endl; 
  for(int i=1;i<=nyrs;i++)  
     report<<years(i)<<" "<<Recruits(i)<<" "<<F_yr(i)<<" "<<yield(i)<<" "<<Yieldhat(i)/1000<<endl; 
     
  report<<"yrs Pop B.in.kg ExploitB"<<endl;   
  for(int i=1;i<=nyrs;i++)     
     report<<years(i)<<" "<<Pop(i)<<" "<<B(i)<<" "<<EB(i)<<endl;
     
  report<<"yrs N"<<endl; 
  for(int i=1;i<=nyrs;i++)
     report<<years(i)<<" "<<trans(N)(i)<<endl; 
  report<<"2018"<<" "<<trans(N)(nyrs+1)<<endl; 
  
  report<<"yrs Spawners"<<endl; 
  for(int i=1;i<=nyrs;i++)
     report<<years(i)<<" "<<trans(Spawners)(i)<<endl; 
  report<<"2018"<<" "<<trans(Spawners)(nyrs+1)<<endl; 
  
  cout<<"Hi"<<endl; 
  //exit(13);
     
  report<<"lengthfq: "<<lengthfq<<endl; 
  report<<"Catchhat: "<<Catch<<endl; 
  report<<"max.grad: "<<objective_function_value::gmax<<endl; 
  report<<"Sample size in the length frequency: "<<SamSize<<endl;

  report<<"x(1,nlengths): "<<x<<endl;   
  report<<"Mu(1,nlengths): "<<Mu<<endl;

  report<<"SS(1,nages): "<<SS<<endl;
  report<<"M, sigmaL, and aic: "<<M<<" "<<sigmaL<<" "<<aic<<endl; 
  
  report<<"sig2_yield: "<<sig2_yield<<endl; 
  
  //exit(13); 
  
RUNTIME_SECTION
  maximum_function_evaluations 100,150,300,10000
  convergence_criteria .01,.0001,1e-7

TOP_OF_MAIN_SECTION
  gradient_structure::set_MAX_NVAR_OFFSET(1000);  //maximum number of depdendent variables of 400 exceeded 
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(1000);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(100000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(1000000);
  arrmblsize=900000;
  
GLOBALS_SECTION
  #include <admodel.h>
  #include <math.h>
  #include <stdio.h>
  #include <stddef.h>
  #include <stdlib.h>



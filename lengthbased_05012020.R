library(TMB);

# Convert the code in ADMB to TMB
# Target stock: the Korean mackerel, i.e., understand the data

sbmack <- "
 //A size-based (sb) model for the Korean mackerel (mack) stock assessment;
 //Author: Saang-Yoon Hyun, KH Kim, JW Gim, MK Bahk, and Doyul Kim as of April 29, 2020; 
 
 #include <TMB.hpp>
 
 // pass missing values
 template<class Type>
 bool isNA(Type x){
    return R_IsNA(asDouble(x));
 }
  
 // square
 template<class Type>
 Type square(Type i) {
    return i*i; 
 }
 
 //objective function
 template<class Type>
 Type objective_function<Type>::operator() () {
 //Data
 DATA_INTEGER(nages);    //number of imaginary age classes; 
 DATA_INTEGER(nlengths); //number of length classes;
 DATA_INTEGER(nyrs);     //number of years of data on yield and cpue 
 DATA_MATRIX(yieldcpue); //23 rows (23 years) and 3 columns (yr, yield, cpue);
                         //year 1995 through year 2017
 DATA_INTEGER(indexMinyrLD); //index of 6
                             //1995+(6-1)= 2000, the first year of length data
 DATA_VECTOR(x);          //a vector that has the midpoint of each length class;
 DATA_MATRIX(lengthfrq); //length frequency data (18 years x 42 classes); 
                         //Do NOT use IMATRIX;  //Do found it; 
                         
 DATA_VECTOR(par_logistic); //two parameters in the logistic relation of sexual maturation rate against body length; 
 DATA_VECTOR(musig2_r);   //the mean and the variance of lengths at age 1;
 DATA_VECTOR(abWL);       //alpha and beta in the weight-length relationship; 
 DATA_VECTOR(M);          //Instantaneous natural mortality by length class; 
                          //We could allow M to vary by length class 
                          //this is NM(x) in Quinn's code;
 DATA_VECTOR(lambdas);    //two weights allocated to the corresponding likelihood //a little different from ADMB TPL code 
 DATA_SCALAR(Linf);       //the possible maximum of the body length; 
 DATA_SCALAR(neff);       //effective sample size in the multinomial likelihood for the length data; //a little different from ADMB TPL code
 DATA_INTEGER(choiceSR);  //choice of a stock-recruitment model (1 = B_H or 2 = Ricker);  
 DATA_SCALAR(sig2_logYield); //the variance of log(yield);   
 
 
 //Parameters
 PARAMETER_VECTOR(logRec); //log(Rec) whose size is the same as indexMinyrLD (6)
 PARAMETER(logkappa);      //
 PARAMETER(logL50);        //gear selectivity parameter;
 PARAMETER(loggamma);      //gear selectivity parameter;
 PARAMETER(logsig_L);      //log(sig_L), where L is the L_{a+1} equation; 

 PARAMETER(log_F1);        //log(F in the first time, i.e., first year);
 PARAMETER_VECTOR(F_devs); //random walk deviations for annual fishing mortality in the following year;
 PARAMETER(logq_index);    //this q is from index (or cpue) = q*Bimoass;
 
 
 PARAMETER(logaSR);        //parameter a in the stock-recruitment model 
 PARAMETER(logbSR);        //parameter b in the stock-recruitment model
 
 //Derived quantities
 vector<Type> L=x;              //Length after growth; 
 Type nyrsLD=(nyrs-Type(5));    //the number of years of length data;  //i.e., 18 years of 2000-2017; 
 Type kappa=exp(logkappa); 
 Type Rho=exp(Type(-1.0)*kappa);//Brody coefficient;
 Type L50=exp(logL50);          //gear selectivity
 Type gamma=exp(loggamma);      //gear selectivity
 Type sig_L=exp(logsig_L);      //log(sig_L), where L is the L_{a+1} equation; 

 Type r = Type(1);              //recruit is defined as the pop size at age 1; 
 Type ncohorts=nyrs; 
 Type q_index=exp(logq_index);  //this q is from index (or cpue) = q*Bimoass;
 
 vector<Type> Rec(nyrs+1);      //recruitment
 Type aSR=exp(logaSR); 
 Type bSR=exp(logbSR); 
 
 vector<Type> Wt(nlengths);          //body weight by length 
 vector<Type> Sel(nlengths);         //gear selectivity by length class; 
 vector<Type> Maturation(nlengths);  //maturation rate by length class; 

 matrix<Type> f(nages,nlengths);     //length frequency as pmf //f(a,x,sex) as 3d_array in Quinn's code; 
 array<Type> pp(nlengths,nlengths,nages); //pp(L,x,a); //at the level of each age,
                                                        //where nlengths (row) x nlengths (columns) at each level;
 vector<Type> Mu(nlengths);   //diffe by length
 vector<Type> SS(nages);      //differ by age
 
 vector<Type> F(nyrs);             //instantaneous fishing mortality by year; 
 vector<Type> log_F(nyrs);           //log(F);  
 matrix<Type> F_tx(nyrs,nlengths); //instantaneous fishing mortality;
 matrix<Type> Z_tx(nyrs,nlengths); 
 matrix<Type> ExpZ_tx(nyrs,nlengths); 
 
 array<Type> NL(ncohorts+1,nlengths,nages); //at the level of each age
 matrix<Type> N(nages,nyrs); 
 
 array<Type> SpawnersL(ncohorts+1,nlengths,nages); //at the level of each age
 matrix<Type> Spawners(nages,nyrs+1);   
 
 vector<Type> p(nlengths); 
 vector<Type> p_plus(nlengths);  //the last age class; 
 
 vector<Type> TCatch(nyrs); 
 matrix<Type> Catch(nyrs,nlengths); 
 vector<Type> Yieldhat(nyrs);
 vector<Type> Pop(nyrs); 
 array<Type> ENx(ncohorts,nlenghts,nages); 
 vector<Type> EN(ncohorts);
 vector<Type> B(ncohorts); 
 vector<Type> EB(ncohorts); 
 
 matrix<Type> hatlengthfq(nyrsLD, nlengths);   //predicted length frequency
 vector<Type> SamSize(nyrsLD); 
 
 vector<Type> elem_obj2(indexMinyrLD,nyrs);    //elements in the normal likelihood for log(Yield) 
 
 vector<Type> nll(2);      //elements of the objective funtion, which is the negative loglikelihood;
 
 //instantaneous fishing mortality by time, i.e., year;
 log_F(0)=log_F1;                 //don't forget zero index in TMB
 F(0)=exp(log_F(0));
 for(int t=1; t<nyrs; t++)  {
    log_F(t)=log_F(t-1)+F_devs(t-1); 
    F(t)=exp(log_F(t)); 
 };    
 
 //Weight, gear selectivity and maturation rate by length
 for(int xind=0; xind<nlengths; xind++)  {  //don't forget zero index in TMB
     Wt(xind)=abWL(0)*pow(x(xind),abWL(1))/Type(1000);   //the division of 1000 is to convert gram to kg;
     Sel(xind)=Type(1.0)/(Type(1.0)+exp(-gamma*(x(xind)-L50)) );  //gear selectivity 
     Maturation(xind)=exp(par_logistic(0)+par_logistic(1)*x(xind))/(1+exp(par_logistic(0)+par_logistic(1)*x(xind))  ) //maturation rate 
     
     //mortality and survival rate by time and length
     for(int t=0; t<nyrs; t++)  { //t is year; //it is m, month in Quinn's code; 
           F_tx(t,xind)=Sel(xind)*F(t); 
           Z_tx(t,xind)=M(xind)+F_tx(t,xind);
           ExpZ_tx(t,xind)=exp(Type(-1.0)*Z_tx(t,xind)); //survival; 
     };
 };
 
 //LVB body growth; 
 SS(0)=sig2_r;          //don't forget zero index in TMB 
 Type kkk=Type(0.0); 
 
 for(int xind=0; xind<nlengths; xind++)  {
    f(0,xind)=exp(Type(-1.0)*square(x(xind)-mu_r)/(Type(2.0)*SS(0)) );  
    kkk=kkk+f(0,xind); 
 };
 
 for(int xind=0; xind<nlengths; xind++)  {
    f(0,xind)=f(0,xind)/kkk;       //normalize;  //this f(a,x) is to be used later (see the code below);
    Mu(xind)=Linf-(Linf-x(xind))*Rho; 
 };
 
 for(int a=1; a<nages; a++)  { //be careful the index in TMB; //i.e., a=1 means the 2nd element; 
    //this SS is from Cohen and Fishman (1980); //it was used for the shrimp in the Quinn's paper;
    SS(a)=square(sigmaL)*(1.0-pow(Rho,(Type(2.0)*a-Type(2.0)*r)))/(Type(1.0)-square(Rho))+(pow(Rho,(Type(2.0)*a-Type(2.0)*r)))*sig2_r; 
 };
 
 for(int a=0; a<nages; a++){
   for(int xind=0; xind<nlengths; xind++) {
      kkk=Type(0.0); 
      for(int Lind=0; Lind<nlengths; Lind++)  {
          pp(Lind,xind,a)=exp(Type(-1.0)*square(L(Lind)-Mu(xind))/(Type(2.0)*SS(a)) );  // f(L|x) in Quinn et al. (1998); 
          kkk=kkk+pp(Lind,xind,a);
      };
      
      for(int Lind=0; Lind<nlengths; Lind++)
         pp(Lind,xind,a)=pp(Lind,xind,a)/kkk;  //normalize f(L|x);
   };
 };
 
 //recruitment
 for(int i=0;i<indexMinyrLD;i++)
      Rec(i)=exp(logRec(i)); 
 
 //Start of cohort loop
 int m; 
 int a;
 for(int coh=0; coh<ncohorts; coh++) {  //cohort 
    a=Type(0); 
    m=coh; 
        
    N(a,m)=Type(0.0); //JW suggests this; //see its corresponding code in TPL file, lengthSA_m.TPL.
    Spawners(a,m)=Type(0.0); 
    
    if(coh<=indexMinyrLD) {    
        for(int xind=0; xind<nlengths; xind++) {
           NL(m,xind,a)=Rec(m)*f(a,xind); 
           N(a,m)=N(a,m)+NL(m,xind,a); 
           
           SpawnersL(m,xind,a)=NL(m,xind,a)*maturation(xind); 
           Spawners(a,m)=Spawners(a,m)+SpawnersL(m,xind,a); 
        };
    }
    else {
      if(choiceSR==1) 
	       N(a,m)=sum(column(Spawners,m-1))/(aSR+bSR*sum(column(Spawners,m-1)));  //B-H   
	    else if(choiceSR==2) 
	       N(a,m)=aSR*sum(column(Spawners,m-1))*mfexp(Type(-1.0)*bSR*sum(column(Spawners,m-1)));  //Ricker   
	     
	    for(int xind=1;xind<=nlengths;xind++)  {
	       NL(m,xind,a) = N(a,m)*f(a,xind);
	    
     	   SpawnersL(m,xind,a)=NL(m,xind,a)*maturation(xind); 
         Spawners(a,m)=Spawners(a,m)+SpawnersL(m,xind,a);   
        }; 
    };
  
    for(int a=1;a<nages;a++) { //a=1 means age of 2;  
        m=a+coh-1; 
        
        if(m <= nyrs+1)  {
            Type SumP=Type(0.0);
            for(int Lind=0;Lind<nlengths;Lind++)  {
                p(Lind)=0.0;
                for(int xind=0;xind<nlengths;xind++)
                    p(Lind)=p(Lind)+f(a-1,xind)*ExpZ(m-1,xind)*pp(Lind,xind,a);  //Eq. 14 in the AK Sea Grant's paper 
                                                                                //written by Quinn, Turnbull, and Fu
                SumP=SumP+p(Lind); 
            };
        
            if(a!=nages)  {
                
                N(a,m)=Type(0.0);
                Spawners(a,m)=Type(0.0); 
                
                for(int Lind=1;Lind<=nlengths;Lind++)  {
                    f(a,Lind)=p(Lind)/SumP;  //normalize;
                   
                    NL(a,m,Lind)=N(a-1,m-1)*p(Lind);
                              //cf. NL(a,m,Lind)=N(m-1,a-1)*f(a-1,Lind);
                    N(a,m)=N(a,m)+NL(m,Lind,a);
                   
                    SpawnersL(m,Lind,a)=NL(m,Lind,a)*maturation(Lind); 
                    Spawners(a,m)=Spawners(a,m)+SpawnersL(m,Lind,a); 
                };
            }
            
            else if(a==nages) {
                for(int Lind=1;Lind<=nlengths;Lind++)  {
                   p_plus(Lind)=Type(0.0);    
                   for(int xind=1;xind<=nlengths;xind++)
                        p_plus(Lind)=p_plus(Lind)+f(a,xind)*ExpZ(m-1,xind)*pp(Lind,xind,a);  // f(a,xind) is f(6,xind)
                };   
	           
	              N(a,m)=0.0;  
	              Spawners(a,m)=0.0; 
              
                for(int Lind=1;Lind<=nlengths;Lind++)  {
                   NL(m,Lind,a)=N(a-1,m-1)*p(Lind)+N(a,m-1)*p_plus(Lind);   //new p_plus
                              //NL(m,Lind,a)=N(a-1,m-1)*f(a-1,Lind)+N(a,m-1)*f(a,Lind); 
                   N(a,m)=N(a,m)+NL(m,Lind,a);
                    
                   SpawnersL(m,Lind,a)=NL(m,Lind,a)*maturation(Lind); 
                   Spawners(a,m)=Spawners(a,m)+SpawnersL(m,Lind,a); 
                };     
                f(a)=NL(a,m)/N(a,m);  //normalize; f(6) is a vector    
            };
        };  //m ends here 
    }; //a ends here       
 };  //this brace corresponds to cohort;
 
 for(int i=indexMinyrLD;i<nyrs+1;i++)
      Rec(i)=N(1,i);
        
 //================  //note m starts at indexMinyrLD;
 for(int m=indexMinyrLD;m<=nyrs;m++) {  //note m starts at indexMinyrLD;
	  TCatch(m)=Type(0.0);
	  Catch(m)=Type(0.0);       //Catch(1,nyrs,1,nlengths);  
	  Yieldhat(m)=Type(0.0); 
	  Pop(m)=Type(0.0);      
	  EN(m)=Type(0.0);       
	  B(m)=Type(0.0);        
	  EB(m)=Type(0.0);      
	   
	  for(int a=1;a<=nages;a++)   {
	      for(int xind=1;xind<=nlengths;xind++) {
		        CNum=NL(m,xind,a)*(F_tx(m,xind)/Z(m,xind))*(1-ExpZ(m,xind));
		        CWt=CNum*Wt(xind);  //in kg
		        Catch(m,xind)=Catch(m,xind)+CNum;
		        TCatch(m)=TCatch(m)+CNum; 
		        Yieldhat(m)=Yieldhat(m)+CWt; 
		     
		        Pop(m)=Pop(m)+NL(m,xind,a);   //population; 
		        ENx(m,xind,a)=NL(m,xind,a)*Sel(xind); //Exploitable population; 
		        EN(m)=EN(m)+ENx(m,xind,a);            //Exploitable population; 
		        B(m)=B(m)+NL(m,xind,a)*Wt(xind);      //Biomass;  //in kg
		        EB(m)=EB(m)+ENx(m,xind,a)*Wt(xind);   //expoitable biomass; 
	      };
     };  
  };  //m ends here;  
 
  //The expected length-frequency 
  SamSize=lengthfq.rowwise().sum();    //Do suggestion works;
  for(int m=indexMinyrLD;m<=nyrs;m++)  //note m starts at indexMinyrLD, 6;
     for(int xind=1;xind<=nlengths;xind++)    
	      hatlengthfq(m-5,xind)=(Catch(m,xind)/sum(Catch(m)))*SamSize(m-5);  //<== as of 4/27
	 
	//objective function
  nll.setZero(); 
	
  //part 1 of the objective funcion: multinomial
  Type logmult=Type(0.0);
  for(int i=indexMinyrLD-1;i<nyrs;i++) {  //indexMinyrLD = 6;     
     logmult+=lgamma(SamSize(i-5)+1);    //<== as of 4/27
     for(int xind=1;xind<=nlengths;xind++) {     
        logmult+=Type(-1.0)*lgamma((lengthfq(i-5,xind)+1))+lengthfq((i-5),xind)*log(Catch(i,xind)/sum(Catch(i)));      //<== as of 4/27        
     };
  };     
 
  nll(0) +=lambdas(0)*(Type(-1.0)*logmult); 
  
  Type maxloglike=Type(0.0); 
  maxloglike+=logmult; 
 
  //part 2 of the objective function: lognormal
  for(int m=indexMinyrLD-1;m<=nyrs-1;m++) 
      nll(1) +=lambdas(1)*Type(-1.0)*dnorm(log(yield(m)), log(Yieldhat(m)/Type(1000))+sig2_yield/(2.0*square(Yieldhat(m)/1000)), sqrt(sig2_yield), true);
             //elem_obj2(m)=square( log(yield(m))-log(Yieldhat(m)/1000)+sig2_yield/(2.0*square(Yieldhat(m)/1000))  );  //in MT 
             //the above is based on the delta method
             //obj+=lambda2*(0.5*(nyrs-indexMinyrLD+1)*log(sig2_yield)+sum(elem_obj2)/(2.0*sig2_yield));  //lamda2*(the negative normal loglikelihood);

  Type lognormal=Type(0.0); 
  lognormal=(-0.5*(nyrs-indexMinyrLD+1)*log(2*M_PI)-0.5*(nyrs-indexMinyrLD+1)*log(sig2_yield)-sum(elem_obj2)/(2.0*sig2_yield));
  
  maxloglike+=lognormal;  
  
  //Type aic=Type(-2.0)*maxloglike+Type(2.0)*(the number of free parameters);   
  
  Type jnll=nll.sum(); 
  

}  "
 
#compile
write(sbmack, file="sbmack.cpp"); 
compile("sbmack.cpp");
dyn.load(dynlib("sbmack")); 

 
 


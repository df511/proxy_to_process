/*********************************************************************************
***** This macro assigns both dimensions of the 2021 canadian           	 *****
***** deprivation index to any administrative file with a postal        	 *****
***** code variable and, optionally, a census subdivision variable.     	 *****
***** Four versions, canadian, regional, metropolitan and by zone,     	  	 *****
***** are assigned. A simplified version of the PCCF+, a               	  	 *****
***** Statistics Canada product, is used to link the dissemination      	 *****
***** area with the postal code. More specifically, three               	 *****
***** conversion tables were used:                                      	 *****
*****  - WCF from the PCCF+                                             	 *****
*****  - PCCFDUPS from the PCCF+                                        	 *****
*****  - PCCFUNIQ from the PCCF+                                        	 *****
*****                                                                   	 *****
***** An additional file (MUNIC2021) creates a correspondance of        	 *****
***** municipal codes in Canada from 1991 to 2021.			            	 *****
*****                                                                   	 *****
***** Macro parameters:                                                 	 *****
*****  - IN: SAS input file on which the user wants to add the          	 *****
*****        2021 material and social deprivation index (4 versions)    	 *****
*****  - PCODE: six-character postal code variable name                 	 *****
*****  - MUNIC: census subdivison variable name. By default,            	 *****
*****           0 = the index is assigned with the postal code only     	 *****
*****  - OUT: SAS output file which will be identical to the                 *****
*****         input file, but enhanced with the deprivation index variables  *****  
*****                                                                        *****
***** Authors: Denis Hamel, INSPQ                                    		 *****			
*****          Philippe Gamache, INSPQ                               		 *****			
*****          Jérémie Sylvain-Morneau, INSPQ                           	 *****
***** Last modified: November 15 2023     		                       		 *****	
***** PCCF+ version : 8A (december 2022)                                	 *****       
 *********************************************************************************/



/****************************************************************************************/
/* 1. establish the path to the folder where the 4 datafiles are                        */
%let repertory = D:\Projets;



/****************************************************************************************/
/* 2. define the four parameters (by default, munic = 0)     						    */
%let in = input_file;
%let pcode = postalcode;
%let munic = 0;
%let out = output_file;



/****************************************************************************************/
/* 3. Start of the macro. Do not change anything from here.						    	*/
%macro assignation2021(in=, pcode=, munic=, out=);

filename fcp "&repertoire\fcp.txt";
filename pccfuniq "&repertoire\fccpindiceuniq.txt";
filename pccfdups "&repertoire\fccpindicedouble.txt";
filename munic "&repertoire\munic2021.txt";

data wcf;
  infile fcp;
  input munic 1-12 pcode $ 13-18 propmunic 19-23 .2 propcp 5.2 first 29 
        QuintMat 30 QuintSoc 31 QuintMatRC 32 QuintSocRC 33 QuintMatZone 34 QuintSocZone 35 QuintMatRMR 36 QuintSocRMR 37 
        zone 38 PR $ 39-40 css 42 RC $ 43-57 rmr 58-60;
run;

data pccfuniq;
  infile pccfuniq;
  input munic 1-12 pcode $ 13-18 propmunic 19-23 .2 propcp 5.2 first 29 
        QuintMat 30 QuintSoc 31 QuintMatRC 32 QuintSocRC 33 QuintMatZone 34 QuintSocZone 35 QuintMatRMR 36 QuintSocRMR 37 
        zone 38 PR $ 39-40 css 42 RC $ 43-57 rmr 58-60;
run;

data pccfdups;
  infile pccfdups;
  input munic 1-12 pcode $ 13-18 propmunic 19-23 .2 propcp 5.2 first 29 
        QuintMat 30 QuintSoc 31 QuintMatRC 32 QuintSocRC 33 QuintMatZone 34 QuintSocZone 35 QuintMatRMR 36 QuintSocRMR 37 
        zone 38 PR $ 39-40 css 42 RC $ 43-57 rmr 58-60;
run;

data munic;
  infile munic;
  input munic 1-7 mun2021 8-14;
run;

data wcf;
  set wcf;
  nuniq = _n_;
run;

data pccfuniq;
  set pccfuniq;
  nuniq = _n_;
run;

data pccfdups;
  set pccfdups;
  nuniq = _n_;
run;


data fichin;
  %if &munic ne 0 %then %do; 
    length munic 8.; 
  %end;
  set &in;
  numero=_n_;
  nba=ranuni(766);
  pcode=compress(&pcode);
  %if &munic ne 0 %then 
     %do; 
       		if substr(pcode,1,1) in ("A")                 then munic=1000000+mod(&munic,1000000);
       else if substr(pcode,1,1) in ("B")                 then munic=1200000+mod(&munic,1200000);
       else if substr(pcode,1,1) in ("C")                 then munic=1100000+mod(&munic,1100000);
       else if substr(pcode,1,1) in ("E")                 then munic=1300000+mod(&munic,1300000);
       else if substr(pcode,1,1) in ("G","H","J")         then munic=2400000+mod(&munic,2400000);
       else if substr(pcode,1,1) in ("K","L","M","N","P") then munic=3500000+mod(&munic,3500000);
       else if substr(pcode,1,1) in ("R")                 then munic=4600000+mod(&munic,4600000);
       else if substr(pcode,1,1) in ("S")                 then munic=4700000+mod(&munic,4700000);
       else if substr(pcode,1,1) in ("T")                 then munic=4800000+mod(&munic,4800000);
       else if substr(pcode,1,1) in ("V")                 then munic=5900000+mod(&munic,5900000);
       else if substr(pcode,1,1) in ("X")                 then munic=6100000+mod(&munic,6100000);
       else if substr(pcode,1,1) in ("Y")                 then munic=6000000+mod(&munic,6000000);
     %end;
run;

%presence(var=quintmat);

%if &presence ne 0 %then %do;
data fichin;
  set fichin;
  where zone in (.,0);
run;
%end;

data fichin;
  set fichin;
  keep numero nba pcode %if &munic ne 0 %then %do; munic %end;;
run;


/*********************************************************************************************************/
/* Deprivation index assignment with the postal code and the census subdivision variable (if munic ne 0) */

%if &munic ne 0 %then %do;

proc sort data=munic;by munic;run;
proc sort data=fichin;by munic;run;

data fichin;
  merge munic(in=in2) fichin(in=in1);
  by munic;
  if in1*in2 then munic=mun2021;
  else if in2 then delete;
run;

/***** With PCCFUNIQ *****/
/*** Canada ***/
proc sql;
  create table test as
  select fichin.pcode, fichin.munic, fichin.numero, fichin.nba, pccfuniq.first, 
         pccfuniq.QuintMat, pccfuniq.QuintSoc, pccfuniq.quintmatrc, pccfuniq.quintsocrc, 
         pccfuniq.quintmatzone, pccfuniq.quintsoczone, pccfuniq.quintmatrmr, pccfuniq.quintsocrmr, 
         pccfuniq.propmunic, pccfuniq.nuniq, pccfuniq.zone, pccfuniq.pr, pccfuniq.css, pccfuniq.RC, pccfuniq.rmr
  from fichin, pccfuniq
  where fichin.pcode=pccfuniq.pcode and fichin.munic=pccfuniq.munic
  order by numero,nuniq;
quit;

proc sort data=test;by pcode munic numero;run;

/*** Choice of one deprivation index if there's more than one index by postal code ***/
data test1;
  set test;
  by pcode munic numero;
  retain cumul find direct;
  if first.numero then do; 
    find=0; 
    cumul=0; 
    if last.numero then direct=1; 
    else direct=0;
  end;
  if find = 0 then do;
    cumul=cumul+(propmunic/100);
	if nba < cumul then do; 
      find=1; output; 
    end;
  end;
run; 

proc sort data=test1;by numero;run;
proc sort data=fichin;by numero;run;

data appmunic1 nonapp1;
  merge fichin(in=in1) test1(in=in2c);
  by numero;
  if in1*in2c then output appmunic1;
  else if in1 then output nonapp1;
run;

/***** With WCF *****/
/*** Canada ***/
proc sql;
  create table test as
  select fichin.pcode, fichin.munic, fichin.numero, fichin.nba, wcf.first, 
         wcf.quintmat, wcf.quintsoc, wcf.quintmatrc, wcf.quintsocrc, 
         wcf.quintmatzone, wcf.quintsoczone, wcf.quintmatrmr, wcf.quintsocrmr, 
         wcf.pr, wcf.propmunic, wcf.nuniq, wcf.zone, wcf.css, wcf.RC, wcf.rmr
  from fichin, wcf
  where fichin.pcode=wcf.pcode and fichin.munic=wcf.munic
  order by numero,nuniq;
quit;

proc sort data=test;by pcode munic numero;run;

/*** Choice of one deprivation index if there's more than one index by postal code ***/
data test1;
  set test;
  by pcode munic numero;
  retain cumul find direct;
  if first.numero then do; 
    find=0; 
    cumul=0; 
    if last.numero then direct=1; 
    else direct=0;
  end;
  if find = 0 then do;
    cumul=cumul+(propmunic/100);
	if nba < cumul then do; 
      find=1; output; 
    end;
  end;
run; 

proc sort data=test1;by numero;run;

proc sort data=nonapp1;by numero;run; 

data appmunic2 nonapp2;
  merge nonapp1(in=in1) test1(in=in2c);
  by numero;
  if in1*in2c then output appmunic2;
  else if in1 then output nonapp2;
run;

/***** With PCCFDUPS *****/
/*** Canada ***/
proc sql;
  create table test as
  select fichin.pcode, fichin.munic, fichin.numero, fichin.nba, pccfdups.first, 
         pccfdups.quintmat, pccfdups.quintsoc, pccfdups.quintmatrc, pccfdups.quintsocrc, 
         pccfdups.quintmatzone, pccfdups.quintsoczone, pccfdups.quintmatrmr, pccfdups.quintsocrmr, 
         pccfdups.propmunic, pccfdups.nuniq, pccfdups.css, pccfdups.zone, pccfdups.pr, pccfdups.RC, pccfdups.rmr 
  from fichin, pccfdups
  where fichin.pcode=pccfdups.pcode and fichin.munic=pccfdups.munic
  order by numero,nuniq;
quit;

proc sort data=test;by pcode munic numero;run;

/*** Choice of one deprivation index if there's more than one index by postal code ***/
data test1;
  set test;
  by pcode munic numero;
  retain cumul find direct;
  if first.numero then do; 
    find=0; 
    cumul=0; 
    if last.numero then direct=1; 
    else direct=0;
  end;
  if find = 0 then do;
    cumul=cumul+(propmunic/100);
	if nba < cumul then do; 
      find=1; output; 
    end;
  end;
run; 

proc sort data=test1;by numero;run;
proc sort data=nonapp2;by numero;run; 

data appmunic3 fichin;
  merge nonapp2(in=in1) test1(in=in2c);
  by numero;
  if in1*in2c then output appmunic3;
  else if in1 then output fichin;
run;

%end;


/**********************************************************************/
/* Deprivation index assignment with the postal code only (munic = 0) */

/***** With PCCFUNIQ *****/
/*** Canada ***/
proc sql;
  create table test as
  select fichin.pcode, fichin.numero, fichin.nba, pccfuniq.first, 
         pccfuniq.quintmat, pccfuniq.quintsoc, pccfuniq.quintmatrc, pccfuniq.quintsocrc, 
         pccfuniq.quintmatzone, pccfuniq.quintsoczone, pccfuniq.quintmatrmr, pccfuniq.quintsocrmr,
         pccfuniq.propcp, pccfuniq.nuniq, pccfuniq.zone, pccfuniq.pr, pccfuniq.css, pccfuniq.RC, pccfuniq.rmr
  from fichin, pccfuniq
  where fichin.pcode=pccfuniq.pcode
  order by numero,nuniq;
quit;

proc sort data=test;by pcode numero;run;

/*** Choice of one deprivation index if there's more than one index by postal code ***/
data test1;
  set test;
  by pcode numero;
  retain cumul find direct;
  if first.numero then do; 
    find=0; 
    cumul=0; 
    if last.numero then direct=1; 
    else direct=0;
  end;
  if find = 0 then do;
    cumul=cumul+(propcp/100);
	if nba < cumul then do; 
      find=1; output; 
    end;
  end;
run; 

proc sort data=test1;by numero;run;
proc sort data=fichin;by numero;run; 

data appcp1 nonapp3;
  merge fichin(in=in1) test1(in=in2c);
  by numero;
  if in1*in2c then output appcp1;
  else if in1 then output nonapp3;
run;

/***** With WCF *****/
/*** Canada ***/
proc sql;
  create table test as
  select fichin.pcode, fichin.numero, fichin.nba, wcf.first, 
         wcf.quintmat, wcf.quintsoc, wcf.quintmatrc, wcf.quintsocrc, 
         wcf.quintmatzone, wcf.quintsoczone, wcf.quintmatrmr, wcf.quintsocrmr, 
         wcf.pr, wcf.propcp, wcf.nuniq, wcf.zone, wcf.css, wcf.RC, wcf.rmr
  from fichin, wcf
  where fichin.pcode=wcf.pcode
  order by numero,nuniq;
quit;

proc sort data=test;by pcode numero;run;

/*** Choice of one deprivation index if there's more than one index by postal code ***/
data test1;
  set test;
  by pcode numero;
  retain cumul find direct;
  if first.numero then do; 
    find=0; 
    cumul=0; 
    if last.numero then direct=1; 
    else direct=0;
  end;
  if find = 0 then do;
    cumul=cumul+(propcp/100);
	if nba < cumul then do; 
      find=1; output; 
    end;
  end;
run; 

proc sort data=test1;by numero;run;
proc sort data=nonapp3;by numero;run; 

data appcp2 nonapp4;
  merge nonapp3(in=in1) test1(in=in2c);
  by numero;
  if in1*in2c then output appcp2;
  else if in1 then output nonapp4;
run;

/***** With PCCFDUPS *****/
/*** Canada ***/
proc sql;
  create table test as
  select fichin.pcode, fichin.numero, fichin.nba, pccfdups.first, 
         pccfdups.quintmat, pccfdups.quintsoc, pccfdups.quintmatrc, pccfdups.quintsocrc, 
         pccfdups.quintmatzone, pccfdups.quintsoczone, pccfdups.quintmatrmr, pccfdups.quintsocrmr, 
         pccfdups.propcp, pccfdups.nuniq, pccfdups.zone, pccfdups.pr, pccfdups.css, pccfdups.RC, pccfdups.rmr
  from fichin, pccfdups
  where fichin.pcode=pccfdups.pcode
  order by numero,nuniq;
quit;

proc sort data=test;by pcode numero;run;

/*** Choice of one deprivation index if there's more than one index by postal code ***/
data test1;
  set test;
  by pcode numero;
  retain cumul find direct;
  if first.numero then do; 
    find=0; 
    cumul=0; 
    if last.numero then direct=1; 
    else direct=0;
  end;
  if find = 0 then do;
    cumul=cumul+(propcp/100);
	if nba < cumul then do; 
      find=1; output; 
    end;
  end;
run; 

proc sort data=test1;by numero;run;
proc sort data=nonapp4;by numero;run; 

data appcp3 nonapp5;
  merge nonapp4(in=in1) test1(in=in2c);
  by numero;
  if in1*in2c then output appcp3;
  else if in1 then output nonapp5;
run;


data tout;
  set %if &munic ne 0 %then %do; appmunic1(in=in1) appmunic2(in=in2) appmunic3(in=in3) %end; 
                                 appcp1(in=in4) appcp2(in=in5) appcp3(in=in6) nonapp5(in=in7);
  if in7 then direct=0;
  if in5 then do; 
    direct=1;
  end;

  %if &munic ne 0 %then %do; 
     if in1 then base=1;
     else if in2 then do; 
       direct=1;base=2;
     end;
	 else if in3 then base=3;
   %end;
   if in4 then base=1;
   else if in5 then base=2;
   else if in6 then base=3;

   if in7=1 then do;
     quintmat = 0;quintsoc = 0;quintmatrc=0;quintsocrc=0;quintmatzone=0;quintsoczone=0;quintmatrmr=0;quintsocrmr=0;zone = 0;css=0; rmr=0;
   end;
   keep numero base quintmat quintsoc quintmatrc quintsocrc quintmatzone quintsoczone quintmatrmr quintsocrmr zone CSS pr RC rmr;
run;

proc sort;by numero;run; 

data fichin;
  set &in;
  numero=_n_;
run;

data tout;
  retain pr CSS QuintMat QuintSoc 
         RC QuintMatRC QuintSocRC 
         zone QuintMatZone QuintSocZone 
         rmr QuintMatRMR QuintSocRMR;
  set tout;
run;



/* English translation */
data tout;

	retain pr CSS QuintMat QuintSoc 
           CR QuintMatRC QuintSocRC 
           zone QuintMatZone QuintSocZone 
           rmr QuintMatRMR QuintSocRMR;

	set tout;

	length CR $30.;
		 if RC = "BC"          then CR = "British Columbia";
	else if RC = "Prairies"    then CR = "Prairies";
	else if RC = "Ontario"     then CR = "Ontario";
	else if RC = "Québec"      then CR = "Quebec";
	else if RC = "Atlantique"  then CR = "Atlantic Canada";
	else                            CR = "";

	rename CSS = SAC
           QuintMatRC = QuintMatCR
		   QuintSocRC = QuintSocCR
		   RMR = CMA
		   QuintMatRMR = QuintMatCMA
		   QuintSocRMR = QuintSocCMA;

	drop RC;
run;


data &out;
  merge fichin tout;
  by numero;
  DATE = '20231115';
  *format groupe fgroupe. base fbase.;
  drop numero;
run;


proc datasets library=work;
  delete corr2 wcf2 pccfuniq2 pccfdups2 test test1 test2 test3 fichin tout appmunic1 appmunic2 appmunic3 appmunic4 appcp1 appcp2 appcp3 appcp4 nonapp1 nonapp2 
         nonapp3 nonapp4 nonapp5 nonapp6 nonapp7 pccfdup2; 
quit;
%mend assignation2021;


%macro presence(var=);
%global presence;
%let dsid=%sysfunc(open(fichin,i));
    %if &dsid>0 %then
       %do;
	     %let presence=%sysfunc(varnum(&dsid,&var));
		 /*** the macrovariable presence = 0 if the variable is not in the file ****/
	   %end;
%put &presence;
%let rc = %sysfunc(close(&dsid));
%mend presence;



/***** Example *****/
%assignation2021(in=&in,pcode=&pcode,munic=&munic,out=&out);



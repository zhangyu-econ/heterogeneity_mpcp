clear all
set matsize 800

global DIR /home/yu/Desktop/yu_cex

global CEXDIR $DIR/data
global PSIDDIR $DIR/bpp2008_6811
global PSIDORGDIR $DIR/bpp2008_original
cd $PSIDDir

set more off
*use $PSIDDIR/data3_cexlong_altiv_yrlinear_fstmp, clear  	//imputed using CEX 72-11
use $PSIDDIR/data3_cexlong_altiv_yrlinear_fstmp_7297, clear 	//imputed using CEX 72-98

**********************************************************************
**********SECTION 1:                                       ***********
**********HOUSE KEEPING **********************************************
**********************************************************************

order person year age agew empst house food fout fstmp bppndur ndur y ftaxsim educexpense hlth* ///
	gas-trans autoinsurance ///
	childcare homemaintn furnishing clothing tripsvacation recreation proptax ///
	state smsa fsize sex ncars race empst self marit educw kids fchg newhd split

tsset person year
gen t=1
egen count = total(t), by(person)
drop t

*gsort -count person year

/*bpp*/
sort year intid
merge year intid using $PSIDDIR/bpp_orig_sample
drop if _merge==2
drop _merge 
*drop if age<30|age>65
drop if age<25|age>65
*drop if year<1980|year>1996
drop if seosample|immsample|latsample
*drop if !bpp
drop if sex == 2 // female heads
drop if race >1  // not enough data to estimate non-white age profiles
replace lumpsum  = 0 if year<=1981
replace trpaidout= 0 if trpaidout==.
replace ftaxsim = ftax if year<1977
replace ftax = ftaxsim if year>1990

gen income = y-asset-(ftaxsim)*((y-asset)/y)   //NO ASSET INCOME BEFORE 1974
*gen income = y-asset -(ftax)*((y-asset)/y)

gen logy = log(income)-lpcpi
*gen logy = log(income-trpaidout+lumpsum)
gen logc = log(bppndur)
*bppndur already real values
*gen logc = log(food+fout+fstmp)-lpfood

recode empst (1=1) (2 3 = 2) (4 = 3) (0 5 6 = .)
save temp, replace

**********************************************************************
**********SECTION 2:                                      ************
**********REMOVE AGE EFFECT AND YEAR EFFECT IN RAW INCOME ************
**************LATER REMOVE COHORT EFFECT BY FIRST DIFFERENCING *******
**********************************************************************

cap log close
log using $DIR/stata05112014_bootstrap_psid1968_1997_remove_age_year_7297, replace text
/*generate age profiles of income*/

//METHOD 1: AGE AND YEAR DUMMIES (AND DEMOGRAPHIC DUMMIES)
/*
use temp, clear

qui reg logy i.year##(i.educ i.race /*i.empst*/ i.region) /*i.age#i.educ */ i.yb  i.fsize i.kids ///
	extra bigcity kidsout
*reg logy i.year
qui predict uy if e(sample), res
qui reg logc i.year##(i.educ i.race /*i.empst*/ i.region) /*i.age#i.educ */ i.yb  i.fsize i.kids ///
	extra bigcity kidsout
*reg logc i.year
qui predict uc if e(sample), res
*/

//METHOD 2: HP FILTER

/*
sort educ race age
collapse logy logc, by(educ race age)
gen raceduc = race*10+educ
tsset raceduc age
tsfilter hp clogy clogc = logy logc, smooth(100) trend(hlogy hlogc)
twoway (connected logc age if educ==3 & race==1) ///
		(connected logc age if educ==2 & race==1) ///
		(connected hlogc age if educ==3 & race==1) ///
		(connected hlogc age if educ==2 & race==1), name(hp_cons)
twoway (connected logy age if educ==3 & race==1) ///
		(connected logy age if educ==2 & race==1) ///
		(connected hlogy age if educ==3 & race==1) ///
		(connected hlogy age if educ==2 & race==1), name(hp_income)
gen dhlogy = D.hlogy
gen dlogy = D.logy
twoway (connected dlogy age if educ==3 & race==1) ///
		(connected dlogy age if educ==2 & race==1) ///
		(connected dhlogy age if educ==3 & race==1) ///
		(connected dhlogy age if educ==2 & race==1), name(hp_dy)
drop raceduc logy logc clogy clogc
sort educ race age
save profile, replace

use temp, clear
sort educ race age
merge educ race age using profile
erase profile.dta
tab _merge
gen lhy = logy - hlogy
qui reg lhy i.year
qui predict uy if e(sample), res
qui reg logc i.year
qui predict uc if e(sample), res
*/

//METHOD 3: FRACTIONAL POLYNOMIALS 

/*
u temp, clear
set more off
xi i.year
fracpoly, degree(3): regress logy age _Iyear* if educ==3
fracpred fplogy, for(age)
fracpoly, degree(3): regress logy age _Iyear* if educ==2
fracpred temp2, for(age)
fracpoly, degree(3): regress logy age _Iyear* if educ==1
fracpred temp1, for(age)
replace fplogy = temp2 if educ==2
replace fplogy = temp1 if educ==1
drop temp1 temp2

fracpoly, degree(3): regress logy age _Iyear* if educ==3&empst==1
fracpred fplogye, for(age)
fracpoly, degree(3): regress logy age _Iyear* if educ==2&empst==1
fracpred temp2, for(age)
fracpoly, degree(3): regress logy age _Iyear* if educ==1&empst==1
fracpred temp1, for(age)
replace fplogye = temp2 if educ==2
replace fplogye = temp1 if educ==1
drop temp2 temp1

fracpoly, degree(3): regress logc age _Iyear* if educ==3
fracpred fplogc, for(age)
fracpoly, degree(3): regress logc age _Iyear* if educ==2
fracpred temp2, for(age)
fracpoly, degree(3): regress logc age _Iyear* if educ==1
fracpred temp1, for(age)
replace fplogc = temp2 if educ==2
replace fplogc = temp1 if educ==1
drop temp1 temp2

fracpoly, degree(3): regress logc age _Iyear* if educ==3&empst==1
fracpred fplogce, for(age)
fracpoly, degree(3): regress logc age _Iyear* if educ==2&empst==1
fracpred temp2, for(age)
fracpoly, degree(3): regress logc age _Iyear* if educ==1&empst==1
fracpred temp1, for(age)
replace fplogce = temp2 if educ==2
replace fplogce = temp1 if educ==1
drop temp2 temp1

collapse fplogy fplogye fplogc fplogce, by(educ year age)
sort educ year age
save profile_fp, replace

use temp, clear

sort educ year age
merge educ year age using profile_fp
*erase profile_hp.dta
tab _merge
drop _merge

*gen uy = logy - fplogy
gen uy = logy - fplogye

qui reg logc i.year
qui predict uc if e(sample), res
*gen uc = logc - fplogc
*/

//METHOD 4: PARTIALLY LINEAR MODEL (Villaverde-Krueger 2004)

//PARLIN: estimate y = m(z) + Xb + e (semiparametrics; epanechnikov kernel)

cap mata mata drop parlin()
mata 
function parlin(yvar,nonparvar,xvars,touse,nonpar,fitted)
{	
  st_view(y   =0,.,yvar,touse) //dependent var: can be anything
  st_view(z   =0,.,nonparvar,touse) //nonpar argument: integer variable
  st_view(X   =0,.,xvars,touse) //controls: can be anything
  
  minZ = min(z)
  maxZ = max(z)
  nZ   = max(z)-min(z)+1
  
  nObs = rows(y)
  S = J(nZ,nObs,.)
  h = 5
  
  for(i=1; i<=nZ; i++) {
  for(j=1; j<=nObs; j++) {
  u = min((abs(minZ+i-1-z[j])/h,1))
  S[i,j]=(.75/h)*(1-u^2) 
  }
  }
  
  yTildeGrid = (S*y):/(S*J(nObs,1,1))
  XTildeGrid = (S*X):/(S*J(nObs,1,1))
  
  yResid = J(nObs,1,.)
  XResid = J(nObs, cols(X),.)
  for(j=1; j<=nObs; j++) {
    yResid[j] = y[j] - yTildeGrid[z[j]-minZ+1]
    XResid[j,.]  = X[j,.]  - XTildeGrid[z[j]-minZ+1,.]
  }
  
  bHat = invsym(XResid'*XResid)*XResid'*yResid
  yZGrid = (S*(y - X*bHat)):/(S*J(nObs,1,1))
  
  yZ = J(nObs,1,.)
  for(j=1; j<=nObs; j++) {
    yZ[j] = yZGrid[z[j]-minZ+1]
  }
  
  yHat = yZ + X*bHat
  
  st_addvar("double",nonpar)
  st_addvar("double",fitted)
  st_store(.,nonpar,touse,yZ)
  st_store(.,fitted,touse,yHat)
}
mata mosave parlin(),replace	
end

use temp, clear
//logy: create selectvar and controls
gen touse =!missing(logy*age*year)
tab year if logy <., gen(yrd)
drop yrd1
*mata: parlin("logy","age","yrd*","touse","logyAge","logyHat")

//empst-nonresticed, three education groups: logy
replace touse = !missing(logy*age*year)&educ==3
mata: parlin("logy","age","yrd*","touse","logyAge","logyHat")
replace touse = !missing(logy*age*year)&educ==2
mata: parlin("logy","age","yrd*","touse","tempA2","tempH2")
replace touse = !missing(logy*age*year)&educ==1
mata: parlin("logy","age","yrd*","touse","tempA1","tempH1")
replace logyAge = tempA2 if educ==2
replace logyAge = tempA1 if educ==1
replace logyHat = tempH2 if educ==2
replace logyHat = tempH1 if educ==1
drop tempA* tempH*

//empst-resticed, three education groups: logy
replace touse = !missing(logy*age*year)&educ==3&empst==1
mata: parlin("logy","age","yrd*","touse","logyEAge","logyEHat")
replace touse = !missing(logy*age*year)&educ==2&empst==1
mata: parlin("logy","age","yrd*","touse","tempA2","tempH2")
replace touse = !missing(logy*age*year)&educ==1&empst==1
mata: parlin("logy","age","yrd*","touse","tempA1","tempH1")
replace logyEAge = tempA2 if educ==2
replace logyEAge = tempA1 if educ==1
replace logyEHat = tempH2 if educ==2
replace logyEHat = tempH1 if educ==1
drop tempA* tempH*

drop yrd*
tab year if logc <., gen(yrd)
drop yrd1

//empst-nonresticed, three education groups: logc
replace touse = !missing(logc*age*year)&educ==3
mata: parlin("logc","age","yrd*","touse","logcAge","logcHat")
replace touse = !missing(logc*age*year)&educ==2
mata: parlin("logc","age","yrd*","touse","tempA2","tempH2")
replace touse = !missing(logc*age*year)&educ==1
mata: parlin("logc","age","yrd*","touse","tempA1","tempH1")
replace logcAge = tempA2 if educ==2
replace logcAge = tempA1 if educ==1
replace logcHat = tempH2 if educ==2
replace logcHat = tempH1 if educ==1
drop tempA* tempH*

//empst-resticed, three education groups: logc
replace touse = !missing(logc*age*year)&educ==3&empst==1
mata: parlin("logc","age","yrd*","touse","logcEAge","logcEHat")
replace touse = !missing(logc*age*year)&educ==2&empst==1
mata: parlin("logc","age","yrd*","touse","tempA2","tempH2")
replace touse = !missing(logc*age*year)&educ==1&empst==1
mata: parlin("logc","age","yrd*","touse","tempA1","tempH1")
replace logcEAge = tempA2 if educ==2
replace logcEAge = tempA1 if educ==1
replace logcEHat = tempH2 if educ==2
replace logcEHat = tempH1 if educ==1
drop tempA* tempH*

drop yrd*

gen uy = logy - logyHat
*gen uy = logy - logyEHat

qui reg logc i.year
qui predict uc if e(sample), res
*gen uc = logc - logcEHat
*gen uc = logc - logcHat

//OPTION A1: ADJUST LOG CONSUMPTION BY EQUIVALENCE SCALE (CANNOT USE WITH METHOD 1)


**********************************************************************
**********SECTION 3:                                         *********
***********HETEROGENEITY REGRESSIONS**********************************
**********************************************************************


************PREPARING FOR AGE HETEROGENEITY IN CONS INSUR************

*drop if empst>1
qui tsset person year
qui gen duc = D.uc
qui gen duy = D.uy
gen ivperm = F2.duy+F.duy+duy+L.duy
ivreg duc (duy=ivperm), nocons
gen duycol = duy*(educ==3)
gen ivpermcol = ivperm*(educ==3)
ivreg duc (duy duycol=ivperm ivpermcol), nocons
*ivreg duc (duy duycol=ivperm ivpermcol) if year<1985, nocons
*ivreg duc (duy duycol=ivperm ivpermcol) if year>1985, nocons

cap drop *rage
gen rphi_gt_rage = .
gen rphi_le_rage = .
gen rphi_at_rage = .
gen sphi_gt_rage = .
gen sphi_le_rage = .
gen sphi_at_rage = .
gen tphi_gt_le_rage = .
gen rage = .

*************Age profile of consumption insurance******************
set more off

forvalues i=27/62 {
	qui local j = `i'+1
	qui cap drop  duyolderthan`i' ivpermolderthan`i' duyyoungerthan`j' ivpermyoungerthan`j'
	qui gen duyolderthan`i' = duy*(age>`i'&!missing(age))
	qui gen ivpermolderthan`i' = ivperm*(age>`i'&!missing(age))
	qui gen duyyoungerthan`j' = duy*(age<=`i'&!missing(age))
	qui gen ivpermyoungerthan`j' = ivperm*(age<=`i'&!missing(age))
	ivreg duc (duyyoungerthan`j' duyolderthan`i' = ivpermyoungerthan`j' ivpermolderthan`i'), nocons
	qui replace rage = `i' in `i'
	qui replace rphi_gt_rage = _b[duyolderthan`i'] in `i'
	qui replace sphi_gt_rage = _se[duyolderthan`i'] in `i'
	qui replace rphi_le_rage = _b[duyyoungerthan`j'] in `i'
	qui replace sphi_le_rage = _se[duyyoungerthan`j'] in `i'
	test duyyoungerthan`j'= duyolderthan`i'
	qui replace tphi_gt_le_rage = r(p) in `i'
	ivreg duc (duy = ivperm) if (age>=(`i'-2) &age<=(`i'+2) &!missing(age)), nocons
	qui replace rphi_at_rage = _b[duy] in `i'
	qui replace sphi_at_rage = _se[duy] in `i'
	qui cap drop  duyolderthan`i' ivpermolderthan`i' duyyoungerthan`j' ivpermyoungerthan`j'
}

graph drop _all
twoway 	(connected rphi_gt_rage rage if rage<60) ///
	(connected rphi_le_rage rage if rage<60) ///
	(connected rphi_at_rage rage if rage<62), name(age)

*graph drop age
	
gen duyagegt52age 	= duy*(age>52)*(age-53)
gen duyagegt52age2 	= duy*(age>52)*(age-53)^2
gen duyagelt30age       = duy*(age<30)*(age-30)

gen ivpermagegt52age 	= ivperm*(age>52)*(age-53)
gen ivpermagegt52age2 	= ivperm*(age>52)*(age-53)^2
gen ivpermagelt30age    = ivperm*(age<30)*(age-30)

****AGE REGRESSION: PARAMETRIC*****
ivreg2 duc (duy duyagegt52age = ivperm ivpermagegt52age) ///
	, nocons cluster(person)

ivreg2 duc (duy duyagegt52age* = ivperm ivpermagegt52age*) ///
	, nocons cluster(person)
test duyagegt52age duyagegt52age2

****AGE REGRESSION: PARAMETRIC, NON-RETIRED*****
ivreg2 duc (duy duyagegt52age = ivperm ivpermagegt52age) if empst<3 ///
	, nocons cluster(person)

ivreg2 duc (duy duyagegt52age* = ivperm ivpermagegt52age*) if empst<3 ///
	, nocons cluster(person)
test duyagegt52age duyagegt52age2

****AGE REGRESSION: PARAMETRIC, YOUTH EFFECT*****
ivreg2 duc (duy duyagegt52age duyagelt30age = ivperm ivpermagegt52age ivpermagelt30age) ///
	, nocons cluster(person)
	
*****PREPARING FOR WEALTH, INCOME, ETC****************

cap drop price
gen     price=exp(lpcpi)

*****Interest rate series ******
gen		rate =4.22	if year ==1967
replace rate =5.66	if year ==1968
replace rate =8.20	if year ==1969
replace rate =7.18	if year ==1970
replace rate =4.66	if year ==1971
replace rate =4.43	if year ==1972
replace rate =8.73	if year ==1973
replace rate =10.50	if year ==1974
replace rate =5.82	if year ==1975
replace rate =5.05	if year ==1976
replace rate =5.54	if year ==1977
replace rate =7.93	if year ==1978
replace rate =11.19	if year ==1979
replace rate =13.36	if year ==1980
replace rate =16.38	if year ==1981
replace rate =12.26	if year ==1982
replace rate =9.09	if year ==1983
replace rate =10.23	if year ==1984
replace rate =8.10	if year ==1985
replace rate =6.81	if year ==1986
replace rate =6.66	if year ==1987
replace rate =7.57	if year ==1988
replace rate =9.22	if year ==1989
replace rate =8.10	if year ==1990
replace rate =5.69	if year ==1991
replace rate =3.52	if year ==1992
replace rate =3.02	if year ==1993
replace rate =4.20	if year ==1994
replace rate =5.84	if year ==1995
replace rate =5.30	if year ==1996
replace rate =5.46	if year ==1997

gen tot_ass=((asset/rate)/price)+(house/price)
gen liq_ass = asset/rate/price
gen ayratio=tot_ass/(income/price)
gen liqratio = liq_ass/(income/price)
gen homeval = house/price
gen hvalratio = homeval/(income/price)

*****************Liq Assets and Consumption Insurance***************
cap drop *qliq *_liq_*
gen rphi_gt_qliq = .
gen rphi_le_qliq = .
gen sphi_gt_qliq = .
gen sphi_le_qliq = .
gen tphi_gt_le_qliq = .
gen rqliq = .

set more off
foreach i of numlist 10 20 30 40 50 60 70 80 90{
	qui cap drop  *_liq_*`i'
	qui egen pctile_liq_`i' = pctile(liq_ass), p(`i') by(year)
	qui gen duy_liq_le_`i' = duy*(liq_ass<=pctile_liq_`i'&!missing(liq_ass))
	qui gen ivperm_liq_le_`i' = ivperm*(liq_ass<=pctile_liq_`i'&!missing(liq_ass))
	qui gen duy_liq_gt_`i' = duy*(liq_ass>pctile_liq_`i'&!missing(liq_ass))
	qui gen ivperm_liq_gt_`i' = ivperm*(liq_ass>pctile_liq_`i'&!missing(liq_ass))
	bootstrap, reps(200): ivreg duc (duy_liq_le_`i' duy_liq_gt_`i' = ivperm_liq_le_`i' ivperm_liq_gt_`i'), nocons
	qui replace rqliq = `i' in `i'
	qui replace rphi_gt_qliq = _b[duy_liq_gt_`i'] in `i'
	qui replace sphi_gt_qliq = _se[duy_liq_gt_`i'] in `i'
	qui replace rphi_le_qliq = _b[duy_liq_le_`i'] in `i'
	qui replace sphi_le_qliq = _se[duy_liq_le_`i'] in `i'
	test duy_liq_le_`i'= duy_liq_gt_`i'
	qui replace tphi_gt_le_qliq = r(p) in `i'
	qui cap drop  *_liq_*`i'
}
	

cap drop meanliq
egen meanliq = mean(liq_ass), by(year)
qui gen duy_liq_le_mean = duy*(liq_ass<=meanliq&!missing(liq_ass))
qui gen duy_liq_gt_mean = duy*(liq_ass>meanliq&!missing(liq_ass))
qui gen ivperm_liq_le_mean = ivperm*(liq_ass<=meanliq&!missing(liq_ass))
qui gen ivperm_liq_gt_mean = ivperm*(liq_ass>meanliq&!missing(liq_ass))
bootstrap, reps(200): ivreg duc (duy_liq_le_mean duy_liq_gt_mean = ivperm_liq_le_mean ivperm_liq_gt_mean), nocons
test duy_liq_le_mean = duy_liq_gt_mean


***************Liq/Y ratio and Consumption Insurance****************

set more off
cap drop *qlratio *_lratio_*
gen rphi_gt_qlratio = .
gen rphi_le_qlratio = .
gen sphi_gt_qlratio = .
gen sphi_le_qlratio = .
gen tphi_gt_le_qlratio = .
gen rqlratio = .

set more off
foreach i of numlist 10 20 30 40 50 60 70 80 90{
  qui cap drop  *_lratio_*`i'
  qui egen pctile_lratio_`i' = pctile(liqratio), p(`i') by(year)
  qui gen duy_lratio_le_`i' = duy*(liqratio<=pctile_lratio_`i'&!missing(liqratio))
  qui gen ivperm_lratio_le_`i' = ivperm*(liqratio<=pctile_lratio_`i'&!missing(liqratio))
  qui gen duy_lratio_gt_`i' = duy*(liqratio>pctile_lratio_`i'&!missing(liqratio))
  qui gen ivperm_lratio_gt_`i' = ivperm*(liqratio>pctile_lratio_`i'&!missing(liqratio))
  bootstrap, reps(200): ivreg duc (duy_lratio_le_`i' duy_lratio_gt_`i' = ivperm_lratio_le_`i' ivperm_lratio_gt_`i'), nocons
  qui replace rqlratio = `i' in `i'
  qui replace rphi_gt_qlratio = _b[duy_lratio_gt_`i'] in `i'
  qui replace sphi_gt_qlratio = _se[duy_lratio_gt_`i'] in `i'
  qui replace rphi_le_qlratio = _b[duy_lratio_le_`i'] in `i'
  qui replace sphi_le_qlratio = _se[duy_lratio_le_`i'] in `i'
  test duy_lratio_le_`i'= duy_lratio_gt_`i'
  qui replace tphi_gt_le_qlratio = r(p) in `i'
  qui cap drop  *_lratio_*`i'
}

cap drop meanliqy
egen meanliqy = mean(liqratio), by(year)
qui gen duy_lratio_le_mean = duy*(liqratio<=meanliqy&!missing(liqratio))
qui gen duy_lratio_gt_mean = duy*(liqratio>meanliqy&!missing(liqratio))
qui gen ivperm_lratio_le_mean = ivperm*(liqratio<=meanliqy&!missing(liqratio))
qui gen ivperm_lratio_gt_mean = ivperm*(liqratio>meanliqy&!missing(liqratio))
bootstrap, reps(200): ivreg duc (duy_lratio_le_mean duy_lratio_gt_mean = ivperm_lratio_le_mean ivperm_lratio_gt_mean), nocons
test duy_lratio_le_mean = duy_lratio_gt_mean


***************Home Values and Consumption Insurance*************

set more off
cap drop *qhomeval *_homeval_*
gen rphi_gt_qhomeval = .
gen rphi_le_qhomeval = .
gen sphi_gt_qhomeval = .
gen sphi_le_qhomeval = .
gen tphi_gt_le_qhomeval = .
gen rqhomeval = .

set more off
foreach i of numlist 10 20 30 40 50 60 70 80 90{
  qui cap drop  *_homeval_*`i'
  qui egen pctile_homeval_`i' = pctile(homeval), p(`i') by(year)
  qui gen duy_homeval_le_`i' = duy*(homeval<=pctile_homeval_`i'&!missing(homeval))
  qui gen ivperm_homeval_le_`i' = ivperm*(homeval<=pctile_homeval_`i'&!missing(homeval))
  qui gen duy_homeval_gt_`i' = duy*(homeval>pctile_homeval_`i'&!missing(homeval))
  qui gen ivperm_homeval_gt_`i' = ivperm*(homeval>pctile_homeval_`i'&!missing(homeval))
  bootstrap, reps(200): ivreg duc (duy_homeval_le_`i' duy_homeval_gt_`i' = ivperm_homeval_le_`i' ivperm_homeval_gt_`i'), nocons
  qui replace rqhomeval = `i' in `i'
  qui replace rphi_gt_qhomeval = _b[duy_homeval_gt_`i'] in `i'
  qui replace sphi_gt_qhomeval = _se[duy_homeval_gt_`i'] in `i'
  qui replace rphi_le_qhomeval = _b[duy_homeval_le_`i'] in `i'
  qui replace sphi_le_qhomeval = _se[duy_homeval_le_`i'] in `i'
  test duy_homeval_le_`i'= duy_homeval_gt_`i'
  qui replace tphi_gt_le_qhomeval = r(p) in `i'
  qui cap drop  *_homeval_*`i'
}

cap drop meanhomeval
egen meanhomeval = mean(homeval), by(year)
qui gen duy_homeval_le_mean = duy*(homeval<=meanhomeval&!missing(homeval))
qui gen duy_homeval_gt_mean = duy*(homeval>meanhomeval&!missing(homeval))
qui gen ivperm_homeval_le_mean = ivperm*(homeval<=meanhomeval&!missing(homeval))
qui gen ivperm_homeval_gt_mean = ivperm*(homeval>meanhomeval&!missing(homeval))
bootstrap, reps(200): ivreg duc (duy_homeval_le_mean duy_homeval_gt_mean = ivperm_homeval_le_mean ivperm_homeval_gt_mean), nocons
test duy_homeval_le_mean = duy_homeval_gt_mean


*****************HomeVal-to-income Ratio and Consumption Insurance********

set more off
cap drop *qhvalratio 
gen rphi_gt_qhvalratio = .
gen rphi_le_qhvalratio = .
gen sphi_gt_qhvalratio = .
gen sphi_le_qhvalratio = .
gen tphi_gt_le_qhvalratio = .
gen rqhvalratio = .

set more off
foreach i of numlist 10 20 30 40 50 60 70 80 90{
  qui cap drop  *_hvalratio_*`i'
  qui egen pctile_hvalratio_`i' = pctile(hvalratio), p(`i') by(year)
  qui gen duy_hvalratio_le_`i' = duy*(hvalratio<=pctile_hvalratio_`i'&!missing(hvalratio))
  qui gen ivperm_hvalratio_le_`i' = ivperm*(hvalratio<=pctile_hvalratio_`i'&!missing(hvalratio))
  qui gen duy_hvalratio_gt_`i' = duy*(hvalratio>pctile_hvalratio_`i'&!missing(hvalratio))
  qui gen ivperm_hvalratio_gt_`i' = ivperm*(hvalratio>pctile_hvalratio_`i'&!missing(hvalratio))
  bootstrap, reps(200): ivreg duc (duy_hvalratio_le_`i' duy_hvalratio_gt_`i' = ivperm_hvalratio_le_`i' ivperm_hvalratio_gt_`i'), nocons
  qui replace rqhvalratio = `i' in `i'
  qui replace rphi_gt_qhvalratio = _b[duy_hvalratio_gt_`i'] in `i'
  qui replace sphi_gt_qhvalratio = _se[duy_hvalratio_gt_`i'] in `i'
  qui replace rphi_le_qhvalratio = _b[duy_hvalratio_le_`i'] in `i'
  qui replace sphi_le_qhvalratio = _se[duy_hvalratio_le_`i'] in `i'
  test duy_hvalratio_le_`i'= duy_hvalratio_gt_`i'
  qui replace tphi_gt_le_qhvalratio = r(p) in `i'
  qui cap drop  *_hvalratio_*`i'
}

cap drop meanhy
egen meanhy = mean(hvalratio),by(year)
qui gen duy_hvalratio_le_mean = duy*(hvalratio<=meanhy&!missing(hvalratio))
qui gen duy_hvalratio_gt_mean = duy*(hvalratio>meanhy&!missing(hvalratio))
qui gen ivperm_hvalratio_le_mean = ivperm*(hvalratio<=meanhy&!missing(hvalratio))
qui gen ivperm_hvalratio_gt_mean = ivperm*(hvalratio>meanhy&!missing(hvalratio))
bootstrap, reps(200): ivreg duc (duy_hvalratio_le_mean duy_hvalratio_gt_mean = ivperm_hvalratio_le_mean ivperm_hvalratio_gt_mean), nocons
test duy_hvalratio_le_mean = duy_hvalratio_gt_mean


***************Total Assets and Consumption Insurance*************

set more off
cap drop *qtotass *_totass_*
gen rphi_gt_qtotass = .
gen rphi_le_qtotass = .
gen sphi_gt_qtotass = .
gen sphi_le_qtotass = .
gen tphi_gt_le_qtotass = .
gen rqtotass = .

set more off
foreach i of numlist 10 20 30 40 50 60 70 80 90{
  qui cap drop  *_totass_*`i'
  qui egen pctile_totass_`i' = pctile(tot_ass), p(`i') by(year)
  qui gen duy_totass_le_`i' = duy*(tot_ass<=pctile_totass_`i'&!missing(tot_ass))
  qui gen ivperm_totass_le_`i' = ivperm*(tot_ass<=pctile_totass_`i'&!missing(tot_ass))
  qui gen duy_totass_gt_`i' = duy*(tot_ass>pctile_totass_`i'&!missing(tot_ass))
  qui gen ivperm_totass_gt_`i' = ivperm*(tot_ass>pctile_totass_`i'&!missing(tot_ass))
  bootstrap, reps(200): ivreg duc (duy_totass_le_`i' duy_totass_gt_`i' = ivperm_totass_le_`i' ivperm_totass_gt_`i'), nocons
  qui replace rqtotass = `i' in `i'
  qui replace rphi_gt_qtotass = _b[duy_totass_gt_`i'] in `i'
  qui replace sphi_gt_qtotass = _se[duy_totass_gt_`i'] in `i'
  qui replace rphi_le_qtotass = _b[duy_totass_le_`i'] in `i'
  qui replace sphi_le_qtotass = _se[duy_totass_le_`i'] in `i'
  test duy_totass_le_`i'= duy_totass_gt_`i'
  qui replace tphi_gt_le_qtotass = r(p) in `i'
  qui cap drop  *_totass_*`i'
}

cap drop meanass
egen meanass = mean(tot_ass), by(year)
qui gen duy_totass_le_mean = duy*(tot_ass<=meanass&!missing(tot_ass))
qui gen duy_totass_gt_mean = duy*(tot_ass>meanass&!missing(tot_ass))
qui gen ivperm_totass_le_mean = ivperm*(tot_ass<=meanass&!missing(tot_ass))
qui gen ivperm_totass_gt_mean = ivperm*(tot_ass>meanass&!missing(tot_ass))
bootstrap, reps(200): ivreg duc (duy_totass_le_mean duy_totass_gt_mean = ivperm_totass_le_mean ivperm_totass_gt_mean), nocons
test duy_totass_le_mean = duy_totass_gt_mean


*****************Wealth-to-income Ratio and Consumption Insurance********

set more off
cap drop *qayratio 
gen rphi_gt_qayratio = .
gen rphi_le_qayratio = .
gen sphi_gt_qayratio = .
gen sphi_le_qayratio = .
gen tphi_gt_le_qayratio = .
gen rqayratio = .

set more off
foreach i of numlist 10 20 30 40 50 60 70 80 90{
  qui cap drop  *_ayratio_*`i'
  qui egen pctile_ayratio_`i' = pctile(ayratio), p(`i') by(year)
  qui gen duy_ayratio_le_`i' = duy*(ayratio<=pctile_ayratio_`i'&!missing(ayratio))
  qui gen ivperm_ayratio_le_`i' = ivperm*(ayratio<=pctile_ayratio_`i'&!missing(ayratio))
  qui gen duy_ayratio_gt_`i' = duy*(ayratio>pctile_ayratio_`i'&!missing(ayratio))
  qui gen ivperm_ayratio_gt_`i' = ivperm*(ayratio>pctile_ayratio_`i'&!missing(ayratio))
  bootstrap, reps(200): ivreg duc (duy_ayratio_le_`i' duy_ayratio_gt_`i' = ivperm_ayratio_le_`i' ivperm_ayratio_gt_`i'), nocons
  qui replace rqayratio = `i' in `i'
  qui replace rphi_gt_qayratio = _b[duy_ayratio_gt_`i'] in `i'
  qui replace sphi_gt_qayratio = _se[duy_ayratio_gt_`i'] in `i'
  qui replace rphi_le_qayratio = _b[duy_ayratio_le_`i'] in `i'
  qui replace sphi_le_qayratio = _se[duy_ayratio_le_`i'] in `i'
  test duy_ayratio_le_`i'= duy_ayratio_gt_`i'
  qui replace tphi_gt_le_qayratio = r(p) in `i'
  qui cap drop  *_ayratio_*`i'
}

cap drop meanay
egen meanay = mean(ayratio),by(year)
qui gen duy_ayratio_le_mean = duy*(ayratio<=meanay&!missing(ayratio))
qui gen duy_ayratio_gt_mean = duy*(ayratio>meanay&!missing(ayratio))
qui gen ivperm_ayratio_le_mean = ivperm*(ayratio<=meanay&!missing(ayratio))
qui gen ivperm_ayratio_gt_mean = ivperm*(ayratio>meanay&!missing(ayratio))
bootstrap, reps(200): ivreg duc (duy_ayratio_le_mean duy_ayratio_gt_mean = ivperm_ayratio_le_mean ivperm_ayratio_gt_mean), nocons
test duy_ayratio_le_mean = duy_ayratio_gt_mean


***************** Residual income and consumption insurance *************

set more off
cap drop *quy 
gen rphi_gt_quy = .
gen rphi_le_quy = .
gen sphi_gt_quy = .
gen sphi_le_quy = .
gen tphi_gt_le_quy = .
gen rquy = .

set more off
foreach i of numlist 10 20 30 40 50 60 70 80 90{
  qui cap drop  *_uy_*`i'
  qui egen pctile_uy_`i' = pctile(uy), p(`i') by(year)
  qui gen duy_uy_le_`i' = duy*(uy<=pctile_uy_`i')
  qui gen ivperm_uy_le_`i' = ivperm*(uy<=pctile_uy_`i')
  qui gen duy_uy_gt_`i' = duy*(uy>pctile_uy_`i')
  qui gen ivperm_uy_gt_`i' = ivperm*(uy>pctile_uy_`i')
  bootstrap, reps(200): ivreg duc (duy_uy_le_`i' duy_uy_gt_`i' = ivperm_uy_le_`i' ivperm_uy_gt_`i'), nocons
  qui replace rquy = `i' in `i'
  qui replace rphi_gt_quy = _b[duy_uy_gt_`i'] in `i'
  qui replace sphi_gt_quy = _se[duy_uy_gt_`i'] in `i'
  qui replace rphi_le_quy = _b[duy_uy_le_`i'] in `i'
  qui replace sphi_le_quy = _se[duy_uy_le_`i'] in `i'
  test duy_uy_le_`i'= duy_uy_gt_`i'
  qui replace tphi_gt_le_quy = r(p) in `i'
  qui cap drop  *_uy_*`i'
}

cap drop meanuy
egen meanuy = mean(uy), by(year)
qui gen duy_uy_le_mean = duy*(uy<=meanuy)
qui gen duy_uy_gt_mean = duy*(uy>meanuy)
qui gen ivperm_uy_le_mean = ivperm*(uy<=meanuy)
qui gen ivperm_uy_gt_mean = ivperm*(uy>meanuy)
bootstrap, reps(200): ivreg duc (duy_uy_le_mean duy_uy_gt_mean = ivperm_uy_le_mean ivperm_uy_gt_mean), nocons
test duy_uy_le_mean = duy_uy_gt_mean


******************** Raw Income and Consumption insurance ****************

set more off
cap drop *qy 
gen rphi_gt_qy = .
gen rphi_le_qy = .
gen sphi_gt_qy = .
gen sphi_le_qy = .
gen tphi_gt_le_qy = .
gen rqy = .

set more off
foreach i of numlist 10 20 30 40 50 60 70 80 90{
  qui cap drop  *_y_*`i'
  qui egen pctile_y_`i' = pctile(y), p(`i') by(year)
  qui gen duy_y_le_`i' = duy*(y<=pctile_y_`i')
  qui gen ivperm_y_le_`i' = ivperm*(y<=pctile_y_`i')
  qui gen duy_y_gt_`i' = duy*(y>pctile_y_`i')
  qui gen ivperm_y_gt_`i' = ivperm*(y>pctile_y_`i')
  bootstrap, reps(200): ivreg duc (duy_y_le_`i' duy_y_gt_`i' = ivperm_y_le_`i' ivperm_y_gt_`i'), nocons
  qui replace rqy = `i' in `i'
  qui replace rphi_gt_qy = _b[duy_y_gt_`i'] in `i'
  qui replace sphi_gt_qy = _se[duy_y_gt_`i'] in `i'
  qui replace rphi_le_qy = _b[duy_y_le_`i'] in `i'
  qui replace sphi_le_qy = _se[duy_y_le_`i'] in `i'
  test duy_y_le_`i'= duy_y_gt_`i'
  qui replace tphi_gt_le_qy = r(p) in `i'
  qui cap drop  *_y_*`i'
}

cap drop meany
egen meany = mean(y), by(year)
qui gen duy_y_le_mean = duy*(y<=meany)
qui gen duy_y_gt_mean = duy*(y>meany)
qui gen ivperm_y_le_mean = ivperm*(y<=meany)
qui gen ivperm_y_gt_mean = ivperm*(y>meany)
bootstrap, reps(200): ivreg duc (duy_y_le_mean duy_y_gt_mean = ivperm_y_le_mean ivperm_y_gt_mean), nocons
test duy_y_le_mean = duy_y_gt_mean


******************** IVPERM and Consumption insurance ****************

set more off
cap drop *qivperm *_ivperm_*
gen rphi_gt_qivperm = .
gen rphi_le_qivperm = .
gen sphi_gt_qivperm = .
gen sphi_le_qivperm = .
gen tphi_gt_le_qivperm = .
gen rqivperm = .

set more off
foreach i of numlist 10 20 30 40 50 60 70 80 90{
  qui cap drop  *_ivperm_*`i'
  qui egen pctile_ivperm_`i' = pctile(ivperm), p(`i') by(year)
  qui gen duy_ivperm_le_`i' = duy*(ivperm<=pctile_ivperm_`i')
  qui gen ivperm_ivperm_le_`i' = ivperm*(ivperm<=pctile_ivperm_`i')
  qui gen duy_ivperm_gt_`i' = duy*(ivperm>pctile_ivperm_`i')
  qui gen ivperm_ivperm_gt_`i' = ivperm*(ivperm>pctile_ivperm_`i')
  bootstrap, reps(200): ivreg duc (duy_ivperm_le_`i' duy_ivperm_gt_`i' = ivperm_ivperm_le_`i' ivperm_ivperm_gt_`i'), nocons
  qui replace rqivperm = `i' in `i'
  qui replace rphi_gt_qivperm = _b[duy_ivperm_gt_`i'] in `i'
  qui replace sphi_gt_qivperm = _se[duy_ivperm_gt_`i'] in `i'
  qui replace rphi_le_qivperm = _b[duy_ivperm_le_`i'] in `i'
  qui replace sphi_le_qivperm = _se[duy_ivperm_le_`i'] in `i'
  test duy_ivperm_le_`i'= duy_ivperm_gt_`i'
  qui replace tphi_gt_le_qivperm = r(p) in `i'
  qui cap drop  *_ivperm_*`i'
}

cap drop meanivperm
egen meanivperm = mean(ivperm), by(year)
qui gen duy_ivperm_le_mean = duy*(ivperm<=meanivperm)
qui gen duy_ivperm_gt_mean = duy*(ivperm>meanivperm)
qui gen ivperm_ivperm_le_mean = ivperm*(ivperm<=meanivperm)
qui gen ivperm_ivperm_gt_mean = ivperm*(ivperm>meanivperm)
bootstrap, reps(200): ivreg duc (duy_ivperm_le_mean duy_ivperm_gt_mean = ivperm_ivperm_le_mean ivperm_ivperm_gt_mean), nocons
test duy_ivperm_le_mean = duy_ivperm_gt_mean

save D:\stata04042014_bootstrap_psid_1968_1997_remove_age_year, replace

******************** Shock Magnitude and Consumption insurance ****************

set more off
cap drop *qivperm *_ivperm_*
gen rphi_gt_qivperm = .
gen rphi_le_qivperm = .
gen sphi_gt_qivperm = .
gen sphi_le_qivperm = .
gen tphi_gt_le_qivperm = .
gen rqivperm = .

set more off
foreach i of numlist 10 20 30 40 50 60 70 80 90{
  qui cap drop  *_ivperm_*`i'
  qui egen pctile_ivperm_`i' = pctile(ivperm), p(`i') by(year)
  qui gen duy_ivperm_le_`i' = duy*(ivperm<=pctile_ivperm_`i')
  qui gen ivperm_ivperm_le_`i' = ivperm*(ivperm<=pctile_ivperm_`i')
  qui gen duy_ivperm_gt_`i' = duy*(ivperm>pctile_ivperm_`i')
  qui gen ivperm_ivperm_gt_`i' = ivperm*(ivperm>pctile_ivperm_`i')
  bootstrap, reps(200): ivreg duc (duy_ivperm_le_`i' duy_ivperm_gt_`i' = ivperm_ivperm_le_`i' ivperm_ivperm_gt_`i'), nocons
  qui replace rqivperm = `i' in `i'
  qui replace rphi_gt_qivperm = _b[duy_ivperm_gt_`i'] in `i'
  qui replace sphi_gt_qivperm = _se[duy_ivperm_gt_`i'] in `i'
  qui replace rphi_le_qivperm = _b[duy_ivperm_le_`i'] in `i'
  qui replace sphi_le_qivperm = _se[duy_ivperm_le_`i'] in `i'
  test duy_ivperm_le_`i'= duy_ivperm_gt_`i'
  qui replace tphi_gt_le_qivperm = r(p) in `i'
  qui cap drop  *_ivperm_*`i'
}

cap drop meanivperm
egen meanivperm = mean(ivperm), by(year)
qui gen duy_ivperm_le_mean = duy*(ivperm<=meanivperm)
qui gen duy_ivperm_gt_mean = duy*(ivperm>meanivperm)
qui gen ivperm_ivperm_le_mean = ivperm*(ivperm<=meanivperm)
qui gen ivperm_ivperm_gt_mean = ivperm*(ivperm>meanivperm)
bootstrap, reps(200): ivreg duc (duy_ivperm_le_mean duy_ivperm_gt_mean = ivperm_ivperm_le_mean ivperm_ivperm_gt_mean), nocons
test duy_ivperm_le_mean = duy_ivperm_gt_mean

save $DIR/stata05112014_bootstrap_psid1968_1997_remove_age_year_7297, replace

******************** End Heterogeneity ***************************


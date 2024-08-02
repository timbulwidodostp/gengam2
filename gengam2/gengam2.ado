*! version 8.0  27apr2002
program define gengam2, eclass byable(recall)
	/* Author: Anirban Basu */
	version 8.0

	syntax [varlist(default=empty)] [if] [in] /*
		*/ [, CLuster(string) CMD Level(integer $S_level) /*
		*/ Robust  FRailty(string) ANCillary(string) ANC2(string) /*
		*/ TIme TR noHR noSHow SCore(passthru) *]
	if _by() {
		_byoptnotallowed score() `"`score'"'
	}
		
capture drop _h* _c*
local chk1: word count `varlist'
if `chk1'==0 {
	di in red "Cannot create schoenfeld residuals for intercept-only model" _n
	di in red "Use -streg- for intercept only model"
			exit 198
}
	

	local passthru `score'
	streg `varlist' `w' `if' `in', cluster(`cluster') dist(gamma) `robust' /*
				*/ frailty(`frailty') `time' `tr' `nohr' `noshow' /*
				*/ `passthru'  ancillary(`ancillary') anc2(`anc2') `options'

quietly{
	predict _xb, xb
	replace _xb=_xb-_b[_cons]
	
	gen  __sig=exp(_b[ln_sig:_cons])
	ereturn local anclist "`ancillary'"
	local wrds : word count `e(anclist)'	
      	tokenize "`e(anclist)'"
 		local i=1
		while `i' <=`wrds' {
		replace __sig = __sig*exp(_b[ln_sig:``i'']*``i'')
		local i=`i' +1
		}
	quietly summ __sig
	global s=_result(3)

	gen __kap=_b[kappa:_cons]
	ereturn local anc2list "`anc2'"
	local wrdk : word count `e(anc2list)'
      	tokenize "`e(anc2list)'"
 		local i=1
		while `i' <=`wrdk' {
		replace __kap = __kap + (_b[ln_sig:``i'']*``i'')
		local i=`i' +1
		}	
	quietly summ __kap
	global k=_result(3)

	local _try=""
	ereturn local xlist "`varlist'"
	local wrd : word count `e(xlist)'
      	tokenize "`e(xlist)'"
 		local i=1
		while `i' <=`wrd' {
		  tempvar expxb mult id sumexp summul exb
		  sort _t
		  gen double `expxb'= exp(-_xb/($s*$k))
		  gen double `mult' = `expxb'*``i''		
		  gen `id'=_N -_n
		  sort `id'
		  gen double `sumexp'=sum(`expxb')
		  gen double `summul' = sum(`mult')
		  sort _t
		  gen double `exb'=`summul'/`sumexp'		  
		  gen _h`i'=``i''-`exb'		  
		  local  _try = "`_try'" + " _h`i'"		  
		local i = `i' +1
	}
	
	ereturn local vl_sch "`_try'"


/* Incase ancillary and anc2 options are parameterized then
   it is needed to get var-cov mattix in terms of ln_sig and kappa */

   mat _V0=e(V)
   local wrd1=`wrd' +1	
   mat __A11=I(`wrd1')
   local anctot=`wrds' + `wrdk' + 2
   mat __A12=J(`wrd1', `anctot',0)		
   mat __A21=J(2, `wrd1',0)	
   local wrds1 = `wrds' +1	
   mat __A22=J(2, `wrds1',0)   	
   local wrdk1 = `wrdk' +1
   mat __A23=J(2, `wrdk1',0)  

      	tokenize "`e(anclist)'"
 		local i=1
		while `i' <=`wrds' {
		tempvar ds
		gen  double `ds' = ``i''
		quietly summ `ds'		
		scalar __s`i'=_result(3)
		tempvar drop ds
		local i = `i' +1
		}
			
		scalar __s`wrds1'=1			
		
	
			
	local i=1
	while `i' <= `wrds1' {
		mat __A22[1, `i'] = __s`i'
	local i = `i' +1
	}

      	tokenize "`e(anc2)'"
 		local i=1
		while `i' <=`wrdk' {
		tempvar ds
		gen  double `ds' = ``i''
		quietly summ `ds'
		scalar __k`i'=_result(3)
		tempvar drop ds
		local i = `i' +1
		}
	
		scalar __k`wrdk1'=1		
	
	
		
	local i=1
	while `i' <= `wrdk1' {
		mat __A23[2, `i'] = __k`i'
	local i = `i' +1
	}

	mat __A1 = __A11, __A12
	mat __A2 = __A21, __A22, __A23
	mat __A0 = __A1 \ __A2

	mat ___V = __A0*_V0*__A0'

global rr=rowsof(___V)	
mat ___V11=_V0[1..`wrd1', 1..`wrd1']
mat ___V12 = ___V[1..`wrd1', `wrd1'+1..$rr]
mat ___V21 = ___V[`wrd1'+1..$rr, 1..`wrd1']
mat ___V22 = ___V[`wrd1'+1..$rr, `wrd1'+1..$rr]

mat colnames ___V12 = ln_sig:_cons kappa:_cons
mat rownames ___V21 = ln_sig:_cons kappa:_cons
mat colnames ___V22 = ln_sig:_cons kappa:_cons
mat rownames ___V22 = ln_sig:_cons kappa:_cons

mat ___V1 = ___V11, ___V12 
mat ___V2 = ___V21, ___V22
mat _V1 = ___V1 \ ___V2


/* Calculate Var-Cov Matrix for beta/(sigma*kappa) */
	
	mat _B1=e(b)
	mat _V2=J(`wrd', `wrd', 0)
	local i=1
	while `i' <= `wrd' {
		local j=1
		while `j' <= `wrd' {
		mat _V2[`i', `j'] = _V1[`i', `j']
		local j = `j' +1
		}
	local i = `i' +1
	}
	
		
	mat _V3=J(`wrd', 2, 0)
	local i=1
	while `i' <= `wrd' {
		local j=1
		while `j' <= 2 {
		local k= `wrd'+ 1 + `j'
		mat _V3[`i', `j'] = _V1[`i', `k']
		local j = `j' +1
		}
	local i = `i' +1
	}


	mat _V4=J(2, `wrd', 0)
	local i=1
	while `i' <= 2 {
		local k= `wrd'+ 1 + `i'
		local j=1
		while `j' <= `wrd' {		
		mat _V4[`i', `j'] = _V1[`k', `j']
		local j = `j' +1
		}
	local i = `i' +1
	}

	mat _V5=J(2, 2, 0)
	local i=1
	while `i' <= 2 {
		local k= `wrd'+ 1 + `i'
		local j=1
		while `j' <= 2 {	
		local l= `wrd'+ 1 + `j'	
		mat _V5[`i', `j'] = _V1[`k', `l']
		local j = `j' +1
		}
	local i = `i' +1
	}
	
	mat _V11=_V2, _V3
	mat _V21=_V4, _V5
	mat _Var=_V11 \ _V21
	
	local db = - 1/($s*$k)
	
	mat _R1=J(`wrd', 1, `db')
	mat _R11=diag(_R1)
	
	mat _R21=J(2, `wrd', 0)
		local j=1
		while `j' <= `wrd' {	
		mat _R21[1, `j'] = _B1[1,`j']/($s*$k)
		local j = `j' +1
		}
		local j=1
		while `j' <= `wrd' {	
		mat _R21[2, `j'] = _B1[1,`j']/($s*($k^2))
		local j = `j' +1
		}

	mat _R=_R11\_R21
	mat _V=_R'*_Var*_R
	*mat _VI = syminv(_V)
	
/* Calc Scaled Schoenfield Res */
	local _try=""
	local i=1
		while `i' <= `wrd' {
		gen _c`i'=0
		local j=1
		while `j' <= `wrd' {
		  replace _c`i' = _c`i' + _h`j'*_V[`j',`i']		  
		  local j= `j' +1
		}
		local  _try = "`_try'" + " _c`i'"	
		local i = `i' +1
	}

	local i=1
	while `i' <= `wrd'{
		replace _c`i' = -(_B1[1, `i']/($s*$k)) + (_N*_c`i')
	local i=`i' +1
	}

ereturn local vl_ssc "`_try'"

local i=1 
tokenize "`e(xlist)'"
while `i' <=`chk1'{
label var _h`i' "``i''"
label var _c`i' "``i''"
local i = `i' +1
}


/*
local i=1 
global rname " "
global eqname ""
tokenize "`e(xlist)'"
while `i' <=`wrd'{

 global rname = "$rname"  +  "``i'' "
 global eqname = "$eqname" + "_t "
local i= `i' +1
}

mat rownames _V1=$rname _cons _cons _cons
mat colnames _V1=$rname _cons _cons _cons
mat roweq _V1 = $eqname _t ln_sig kappa
mat coleq _V1 = $eqname _t ln_sig kappa
*/


local wrd3 = `wrd' + 3
mat __eb=e(b)
mat __b1=__eb[1, 1..`wrd1']

mat __b2=J(1,2,0)
mat __b2[1,1]=ln($s)
mat __b2[1,2]=$k
mat colnames __b2 = ln_sig:_cons kappa:_cons
mat __b = __b1, __b2
mat rownames __b=y1




/* ln_sig for testing */

tokenize "`e(anclist)'"
local i=1
global lnseq ="[ln_sig]_cons"
while `i' <=`wrds' {
	global lnseq = "$lnseq" + " + [ln_sig]``i''"
	local i= `i' +1
}



/* kappa for testing */

tokenize "`e(anc2)'"
local i=1
global lnkeq ="[kappa]_cons"
while `i' <=`wrdk' {
	global lnkeq = "$lnkeq" + " + [kappa]``i''"
	local i= `i' +1
}




tempname b V oldest
mat `b' = __b
mat `V' = _V1




tempname origb origV
mat `origb'=e(b)
mat `origV'=e(V)


local _nobs=e(N)
local _df = e(df_m)



nobreak {
_estimates hold `oldest'
ereturn post `b' `V', o(`_nobs') 
ereturn local cmd "gamma"
ereturn local cmd2 "streg"


/* Test for Gamma */
testnl exp([ln_sig]_cons)=[kappa]_cons
local _Gamma=r(chi2)

 /* Test for Log N */
test [kappa]_cons=0
local _Lnr= r(chi2)

 /* Test for Weibull */
test [kappa]_cons=1
local _Wbl=r(chi2)

/* Test for Exp */
testnl ([ln_sig]_cons=0) ([kappa]_cons=1)
local _Expt=r(chi2)

_estimates unhold `oldest'
}


ereturn scalar xN=`wrd'
ereturn matrix bred __b
ereturn matrix Vred _V1
ereturn scalar sigma=$s
ereturn scalar kappa=$k




	noi di in gr _n _n _n "Tests for identifying distributions" _n	
	noi di in smcl in gr "{hline 33}{c TT}{hline 51}"
	noi di in smcl in gr /*
		*/"Distributions                    {c |}     chi2       df       Prob>chi2"	
	noi di in smcl in gr "{hline 33}{c |}{hline 51}"



	noi di in smcl in gr /*
		*/"Std. Gamma   (kappa = sigma)     {c |} "   in ye _col(37) %9.2f `_Gamma' 	_col(50) %5.0f 1  _col(64) %5.4f chiprob(1,`_Gamma')	
	noi di in smcl in gr /*
		*/"Log Normal   (kappa = 0)         {c |} "   in ye _col(37) %9.2f `_Lnr' 	_col(50) %5.0f 1   _col(64) %5.4f chiprob(1,`_Lnr')	
	noi di in smcl in gr /*
		*/"Weibull      (kappa = 1)         {c |} "   in ye _col(37) %9.2f `_Wbl' 	_col(50) %5.0f 1    _col(64) %5.4f chiprob(1,`_Wbl')	
	noi di in smcl in gr /*
		*/"Exponential  (kappa = sigma = 1) {c |} "   in ye _col(37) %9.2f `_Expt' 	_col(50) %5.0f 2   _col(64) %5.4f chiprob(2,`_Expt')	

	noi di in smcl in gr "{hline 33}{c BT}{hline 51}"

noi di in gr _n _n in ye "New Variables for Schoenfield Residuals created: _h*, _c*" _n


tempvar drop _all

drop __sig __kap
mat drop  _V2 _V3 _V4 _V5 _V11 _V21 _Var  _R1 _R11 _R21  _B1 ___V11 ___V12 /*
*/ ___V21 ___V22 ___V1 ___V2 __b1 __b2
drop _xb

}
end


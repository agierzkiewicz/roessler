/* ==========================================================
 * Periodic orbits in the Roessler system
 * Anna Gierzkiewicz and Piotr Zgliczyński
 * ==========================================================
 * Program: 02-Roessler_a47.cpp
 * Date: 15 XII 2020
 * ==========================================================
 * The program contains the proofs related to the Roessler
 * system with a = 4.7, declared in Section 5.
 * ==========================================================*/ 

#include <iostream>
using namespace std;

#include "capd/capdlib.h"
using namespace capd;
using namespace capd::alglib;
using namespace capd::matrixAlgorithms;

#include "utils.h"

// Main function

int main()
{
  	cout.precision(16);
	try
	{
		///===================== definition of the Roessler system with a=4.7 =====================	
		cout << boolalpha;  	
		cout << "===========================================================" << endl;
		cout << "|| Roessler system, a = 4.7" << endl;
		system3d roessler47(interval(4.7));						// The Roessler system with a=4.7 defined
		
		///===================== variables used in Procedure 1:  =====================
		IVector pCloseTo5periodic({-3.860531150992989, 0.0375827119820099});		// starting point for INM 
		IVector containsPeriodicPoint;												// will contain periodic point, also in Procs. 3, 4
		
		///===================== variables used in Procedure 2:  =====================		
		IVector p2({-6.858260447127058, 0.03505366666527084});	
		IMatrix M({{-1., 0.000842495}, {-0.000842495, -1.}});
		
		vector<HSet2D> L(3);														// L[0], L[2] lie on the section, we do not use L[1]
			L[0] = HSet2D(M * IVector({-1.96783, 0.}) + p2 , M , DVector({1.02,2e-4}));			
			L[2] = HSet2D(M * IVector({0.454186, 0.}) + p2 , M , DVector({0.476895,2e-4}));
		
		///===================== variables used in Procedure 3:  =====================		
		IVector pCloseTo2periodic({interval(-4.883924258743838, -4.883924258742273),interval(0.03663128109720402, 0.03663128109729561)});	// starting point for INM 
		
		///===================== variables used in Procedure 4:  =====================		
		IVector pCloseTo4periodic({interval(-6.332251180191445, -6.3322511801862),interval(0.03544579954748642, 0.03544579954785883)});		// starting point for INM 
		
		
		///===================== variables used in Procedure 5:  =====================		
		HSet2D A (IVector ({-6.19384, 0.0356629}) , IMatrix ({{-1., 0.000777754}, {-0.000777754, -1.}}) , DVector({2.66856,4e-4}));			// Attractor's container
		
		///===================== variables used in Procedure 6:  =====================		
		IVector S (2);				
		//IVector S ({interval(-7.186544568881281, -7.165195528898398),interval(0.0344908202443027, 0.03522742411001666)});					// will contain a possible 3-periodic point
		
		///===================== choice of Procedure:  =====================
		int whichCase;
		
		do{
		cout << "===========================================================" << endl;
		cout << "|| Choose the number of procedure You wish to perform:	  " << endl;
		cout << "|| 1. Proof of a 5-periodic orbit's existence by INM	  " << endl;
		cout << "|| 2. Proof of the chain horizontal covering relations	  " << endl;
		cout << "|| 	N0 =Pc=> N2 =Pc=> N2 =Pc^3=> N0					  " << endl;
		cout << "|| 3. Proof of a 2-periodic orbit's existence by INM	  " << endl;
		cout << "|| 4. Proof of a 4-periodic orbit's existence by INM	  " << endl;
		cout << "|| 5. Proof of forward invariance of A					  " << endl;
		cout << "|| 6. Proof of a 3-periodic orbit's non-existence by INM " << endl;
		cout << "|| other: exit	" << endl;
		cout << "===========================================================" << endl;
		cout << "Procedure " ;
		cin >> whichCase;
		cout << endl;
		
		switch (whichCase)
			{
				///===================== Procedure 1:  =====================
				case 1:
				{
					containsPeriodicPoint=roessler47.anyStationaryPoint(pCloseTo5periodic,5); 		// now it contains a 5-periodic point
					cout << "-----------------------------------------------------------" << endl;		
					cout << "5-periodic orbit for Poincare map of the Roessler system with a = 4.7 contained in the set: " << endl;
					cout << "     p1 = " << containsPeriodicPoint << endl;
			
					for(int i=0;i<5;++i)															// prints all orbit
					{
						containsPeriodicPoint=roessler47.P2d(containsPeriodicPoint);
						cout << "P^" << i+1 << "(p1) = " << containsPeriodicPoint << endl;
					};			
					cout << "-----------------------------------------------------------" << endl << endl;			
					break;
				}
				
				///===================== Procedure 2:  =====================
				case 2:
				{					
					cout << "N0 covers N2? ... " << roessler47.covers2D(L[0],L[2],60) << endl;							// if L[0] =P=> L[2]
																																// divided horizontally into 60 pieces
					cout << "N2 covers N2? ... " << roessler47.covers2D(L[2],L[2],180) << endl;							// if L[2] =P=> L[2]
																																// divided horizontally into 180 pieces
					cout << "3rd iterate of N2 covers N0? ... " << roessler47.covers2D(L[2],L[0],500,3,3) << endl;		// if L[2] =P^3=> L[0]
					cout << "----------------------------------------" << endl << endl;											// divided into 500x3 pieces	
					break;
				}			
			
				///===================== Procedure 3:  =====================
				case 3:
				{
					containsPeriodicPoint = roessler47.anyStationaryPoint(pCloseTo2periodic,2);		// now it contains a 2-periodic point		
					cout << "---------------------------------------------------" << endl;
					cout << "2-periodic orbit:" << endl;
					cout << "---------------------------------------------------" << endl;
					cout << "p1 = " << containsPeriodicPoint << endl;
					for(int i=0;i<2;++i)															// prints all orbit
					{
						containsPeriodicPoint=roessler47.P2d(containsPeriodicPoint);
						cout << "P^" << i+1 << "(p1) = " << containsPeriodicPoint << endl;
					};	
					cout << "---------------------------------------------------" << endl << endl;	
					break;
				}	
			
				///===================== Procedure 4:  =====================
				case 4:
				{		
					containsPeriodicPoint = roessler47.anyStationaryPoint(pCloseTo4periodic,4);		// now it contains a 4-periodic point
					cout << "---------------------------------------------------" << endl;
					cout << "4-periodic orbit:" << endl;
					cout << "---------------------------------------------------" << endl;
					cout << "p1 = " << containsPeriodicPoint << endl;
					for(int i=0;i<4;++i)															// prints all orbit
					{
						containsPeriodicPoint=roessler47.P2d(containsPeriodicPoint);
						cout << "P^" << i+1 << "(p1) = " << containsPeriodicPoint << endl;
					};	
					cout << "----------------------------------------" << endl << endl;
					break;
				}
				
				///===================== Procedure 5:  =====================
			
				case 5:
				{		
					cout << "Is A mapped inside A? ... " << roessler47.inside(A,A,450) << endl;		// check if P(A) < A divided horizontally into 450 pieces
					cout << "----------------------------------------" << endl << endl;
					break;
				}
			
				///===================== Procedure 6:  =====================
				case 6:
				{		
					cout << "We detect the area which is not mapped away by the 3rd iteration of P: " << endl;	
					S = roessler47.whatIsNotMappedOutside(A,500,10,3);
					cout << "The 3-periodic point of the attractor may lie only in the box S = " << S << "." << endl;
					cout << "----------------------------------------" << endl;		
					cout << "Check if S contains only one stationary point for P^3: " << endl;
					cout << "Only one 3-periodic point in S? ... " << roessler47.newtonDivided(S,3,3,3) << endl;		
					cout << "----------------------------------------" << endl;		
					cout << "Check if S contains only one stationary point for P: " << endl;
					cout << "Only one 1-periodic point in S? ... " << roessler47.newtonDivided(S,1,1) << endl;	
					cout << "----------------------------------------" << endl << endl;
					break;
				}
				
				default:
					cout << "Goodbye!" << endl;
					break;
				}
			}
			while (whichCase>=1 and whichCase <=6);
		
	}
	catch(exception& e)
  	{
    		cout << "\n\nException caught: "<< e.what() << endl;
  	}
  return 0;
} 

/* OUTPUT:
===========================================================
|| Roessler system, a = 4.7
===========================================================
|| Choose the number of procedure You wish to perform:	  
|| 1. Proof of a 5-periodic orbit's existence by INM	  
|| 2. Proof of the chain horizontal covering relations	  
|| 	N0 =Pc=> N2 =Pc=> N2 =Pc^3=> N0					  
|| 3. Proof of a 2-periodic orbit's existence by INM	  
|| 4. Proof of a 4-periodic orbit's existence by INM	  
|| 5. Proof of forward invariance of A					  
|| 6. Proof of a 3-periodic orbit's non-existence by INM 
|| other: exit	
===========================================================
Procedure 1

A stationary point for P^5 in = {[-3.860531150996197, -3.8605311509898],[0.03758271198192245, 0.03758271198209733]} proven.
-----------------------------------------------------------
5-periodic orbit for Poincare map of the Roessler system with a = 4.7 contained in the set: 
     p1 = {[-3.860531150996197, -3.8605311509898],[0.03758271198192245, 0.03758271198209733]}
P^1(p1) = {[-6.8191913250046, -6.819191324993896],[0.03508216704659804, 0.03508216704660764]}
P^2(p1) = {[-7.831961535752061, -7.831961535733791],[0.0343731813277388, 0.03437318132799605]}
P^3(p1) = {[-5.748327298250972, -5.748327298207539],[0.03590379426606567, 0.03590379426663539]}
P^4(p1) = {[-8.736464582732129, -8.736464582710001],[0.03378750071674676, 0.03378750071723056]}
P^5(p1) = {[-3.860531151013839, -3.860531150972143],[0.03758271198159423, 0.03758271198242558]}
-----------------------------------------------------------

===========================================================
|| Choose the number of procedure You wish to perform:	  
|| 1. Proof of a 5-periodic orbit's existence by INM	  
|| 2. Proof of the chain horizontal covering relations	  
|| 	N0 =Pc=> N2 =Pc=> N2 =Pc^3=> N0					  
|| 3. Proof of a 2-periodic orbit's existence by INM	  
|| 4. Proof of a 4-periodic orbit's existence by INM	  
|| 5. Proof of forward invariance of A					  
|| 6. Proof of a 3-periodic orbit's non-existence by INM 
|| other: exit	
===========================================================
Procedure 2

N0 covers N2? ... true
N2 covers N2? ... true
3rd iterate of N2 covers N0? ... true
----------------------------------------

===========================================================
|| Choose the number of procedure You wish to perform:	  
|| 1. Proof of a 5-periodic orbit's existence by INM	  
|| 2. Proof of the chain horizontal covering relations	  
|| 	N0 =Pc=> N2 =Pc=> N2 =Pc^3=> N0					  
|| 3. Proof of a 2-periodic orbit's existence by INM	  
|| 4. Proof of a 4-periodic orbit's existence by INM	  
|| 5. Proof of forward invariance of A					  
|| 6. Proof of a 3-periodic orbit's non-existence by INM 
|| other: exit	
===========================================================
Procedure 3

A stationary point for P^2 in = {[-4.883924258743846, -4.883924258742264],[0.03663128109720364, 0.03663128109729597]} proven.
---------------------------------------------------
2-periodic orbit:
---------------------------------------------------
p1 = {[-4.883924258743846, -4.883924258742264],[0.03663128109720364, 0.03663128109729597]}
P^1(p1) = {[-8.220951552235825, -8.220951552233453],[0.03411620084268562, 0.03411620084269375]}
P^2(p1) = {[-4.883924258746204, -4.883924258739907],[0.03663128109720341, 0.03663128109729621]}
---------------------------------------------------

===========================================================
|| Choose the number of procedure You wish to perform:	  
|| 1. Proof of a 5-periodic orbit's existence by INM	  
|| 2. Proof of the chain horizontal covering relations	  
|| 	N0 =Pc=> N2 =Pc=> N2 =Pc^3=> N0					  
|| 3. Proof of a 2-periodic orbit's existence by INM	  
|| 4. Proof of a 4-periodic orbit's existence by INM	  
|| 5. Proof of forward invariance of A					  
|| 6. Proof of a 3-periodic orbit's non-existence by INM 
|| other: exit	
===========================================================
Procedure 4

A stationary point for P^4 in = {[-6.332251180191433, -6.332251180186209],[0.03544579954748502, 0.03544579954786024]} proven.
---------------------------------------------------
4-periodic orbit:
---------------------------------------------------
p1 = {[-6.332251180191433, -6.332251180186209],[0.03544579954748502, 0.03544579954786024]}
P^1(p1) = {[-8.470343654125783, -8.47034365411986],[0.03395554856430392, 0.03395554856440434]}
P^2(p1) = {[-4.362667259915144, -4.3626672599017],[0.0371024558505368, 0.03710245585075974]}
P^3(p1) = {[-7.571669540362099, -7.571669540341765],[0.0345497016342597, 0.03454970163427861]}
P^4(p1) = {[-6.33225118021206, -6.332251180165581],[0.03544579954737354, 0.03544579954797172]}
----------------------------------------

===========================================================
|| Choose the number of procedure You wish to perform:	  
|| 1. Proof of a 5-periodic orbit's existence by INM	  
|| 2. Proof of the chain horizontal covering relations	  
|| 	N0 =Pc=> N2 =Pc=> N2 =Pc^3=> N0					  
|| 3. Proof of a 2-periodic orbit's existence by INM	  
|| 4. Proof of a 4-periodic orbit's existence by INM	  
|| 5. Proof of forward invariance of A					  
|| 6. Proof of a 3-periodic orbit's non-existence by INM 
|| other: exit	
===========================================================
Procedure 5

Is A mapped inside A? ... true
----------------------------------------

===========================================================
|| Choose the number of procedure You wish to perform:	  
|| 1. Proof of a 5-periodic orbit's existence by INM	  
|| 2. Proof of the chain horizontal covering relations	  
|| 	N0 =Pc=> N2 =Pc=> N2 =Pc^3=> N0					  
|| 3. Proof of a 2-periodic orbit's existence by INM	  
|| 4. Proof of a 4-periodic orbit's existence by INM	  
|| 5. Proof of forward invariance of A					  
|| 6. Proof of a 3-periodic orbit's non-existence by INM 
|| other: exit	
===========================================================
Procedure 6

We detect the area which is not mapped away by the 3rd iteration of P: 
The 3-periodic point of the attractor may lie only in the box S = {[-7.186544631101602, -7.165195528898398],[0.0344908202443027, 0.03530742411001666]}.
----------------------------------------
Check if S contains only one stationary point for P^3: 
Only one 3-periodic point in S? ... The box {[-7.179428263700533, -7.172311896299465],[0.03476302153287402, 0.03503522282144534]} contains a unique stationary point for P^3
 in = {[-7.175759361589162, -7.175485411855703],[0.03479757277291359, 0.0348476222075226]}
true
----------------------------------------
Check if S contains only one stationary point for P: 
Only one 1-periodic point in S? ... The box {[-7.186544631101602, -7.165195528898398],[0.0344908202443027, 0.03530742411001666]} contains a unique stationary point for P^1
 in = {[-7.175676999516465, -7.17562056752148],[0.03482379625256042, 0.03482750993787596]}
true
----------------------------------------

===========================================================
|| Choose the number of procedure You wish to perform:	  
|| 1. Proof of a 5-periodic orbit's existence by INM	  
|| 2. Proof of the chain horizontal covering relations	  
|| 	N0 =Pc=> N2 =Pc=> N2 =Pc^3=> N0					  
|| 3. Proof of a 2-periodic orbit's existence by INM	  
|| 4. Proof of a 4-periodic orbit's existence by INM	  
|| 5. Proof of forward invariance of A					  
|| 6. Proof of a 3-periodic orbit's non-existence by INM 
|| other: exit	
===========================================================
Procedure 0

Goodbye!

real	0m53,169s
user	0m44,965s
sys	0m0,012s
* 
*/

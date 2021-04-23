/* ==========================================================
 * Periodic orbits in Roessler system
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
		IVector pCloseTo5periodic({-3.885277116885549,0.03755839144423033});		// starting point for INM 
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
					cout << "N0 covers N2? ... " << roessler47.covers2D(L[0],L[2],30) << endl;							// if L[0] =P=> L[2]
																																// divided horizontally into 30 pieces
					cout << "N2 covers N2? ... " << roessler47.covers2D(L[2],L[2],90) << endl;							// if L[2] =P=> L[2]
																																// divided horizontally into 90 pieces
					cout << "3rd iterate of N2 covers N0? ... " << roessler47.covers2D(L[2],L[0],410,1,3) << endl;		// if L[2] =P^3=> L[0]
					cout << "----------------------------------------" << endl << endl;											// divided horizontally into 410 pieces	
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

A stationary point for P^5 in = {[-3.885277116910041, -3.885277116888829],[0.03755839144432487, 0.03755839144485094]} proven.
-----------------------------------------------------------
5-periodic orbit for Poincare map of the Roessler system with a = 4.7 contained in the set: 
     p1 = {[-3.885277116910041, -3.885277116888829],[0.03755839144432487, 0.03755839144485094]}
P^1(p1) = {[-6.858260447161843, -6.858260447127058],[0.03505366666527084, 0.03505366666529888]}
P^2(p1) = {[-7.766631245400215, -7.766631245340529],[0.0344171339279908, 0.03441713392882635]}
P^3(p1) = {[-5.895584354629448, -5.895584354490923],[0.03578591178680491, 0.03578591178862134]}
P^4(p1) = {[-8.722396020084886, -8.722396020038246],[0.0337962993633148, 0.03379629936509001]}
P^5(p1) = {[-3.885277116944484, -3.885277116854386],[0.03755839144371548, 0.03755839144546031]}
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

A stationary point for P^2 in = {[-4.88392425874385, -4.883924258742262],[0.03663128109720355, 0.03663128109729606]} proven.
---------------------------------------------------
2-periodic orbit:
---------------------------------------------------
p1 = {[-4.88392425874385, -4.883924258742262],[0.03663128109720355, 0.03663128109729606]}
P^1(p1) = {[-8.220951552235828, -8.22095155223345],[0.03411620084268561, 0.03411620084269377]}
P^2(p1) = {[-4.883924258746212, -4.883924258739899],[0.03663128109720328, 0.03663128109729633]}
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
P^1(p1) = {[-8.470343654125783, -8.470343654119862],[0.03395554856430393, 0.03395554856440434]}
P^2(p1) = {[-4.36266725991514, -4.3626672599017],[0.03710245585053683, 0.03710245585075971]}
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
The 3-periodic point of the attractor may lie only in the box S = {[-7.186544568881281, -7.165195528898398],[0.0344908202443027, 0.03522742411001666]}.
----------------------------------------
Check if S contains only one stationary point for P^3: 
The box {[-7.17942822222032, -7.172311875559359],[0.03473635486620735, 0.034981889488112]} contains a unique stationary point for P^3
 in = {[-7.175736618877388, -7.175516652987057],[0.03480642820042148, 0.03484077666061459]}
Only one 3-periodic point in S? ... true
----------------------------------------
Check if S contains only one stationary point for P: 
The box {[-7.186544568881281, -7.165195528898398],[0.0344908202443027, 0.03522742411001666]} contains a unique stationary point for P^1
 in = {[-7.175673679228632, -7.17562445696159],[0.03482423076005173, 0.03482706657491325]}
Only one 1-periodic point in S? ... true
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
Procedure 7

Goodbye!
50.16user 0.00system 0:58.89elapsed 85%CPU (0avgtext+0avgdata 5608maxresident)k
0inputs+0outputs (0major+270minor)pagefaults 0swaps


------------------
(program exited with code: 0)
Press return to continue

*/

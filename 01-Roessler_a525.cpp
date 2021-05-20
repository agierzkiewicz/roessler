/* ==========================================================
 * Periodic orbits in the Roessler system
 * Anna Gierzkiewicz and Piotr Zgliczyński
 * ==========================================================
 * Program: 01-Roessler_a525.cpp
 * Date: 15 XII 2020
 * ==========================================================
 * The program contains the proofs related to the Roessler
 * system with a = 5.25, declared in Section 4.
 * ==========================================================*/ 

#include <iostream>
using namespace std;

#include "capd/capdlib.h"
using namespace capd;
using namespace capd::alglib;
using namespace capd::vectalg;
using namespace capd::matrixAlgorithms;

#include "utils.h"

// Main function
int main()
{
  	cout.precision(16);
	try
	{		
		///===================== definition of the Roessler system with a=5.25 =====================		
		cout << boolalpha;  
		cout << "===========================================================" << endl;
		cout << "|| Roessler system, a = 5.25" << endl;
		system3d roessler525(interval(5.25));						// The Roessler system with a=5.25 defined
		
		///===================== variables used in Procedure 1:  =====================
		IVector pCloseTo3periodic({-3.46641, 0.0346316});			// starting point for INM 
		IVector contains3PeriodicPoint(pCloseTo3periodic);			// will contain 3-periodic point
		
		///===================== variables used in Procedure 2:  =====================
		IVector p2({-6.264007533274157, 0.03265435884602701});	
		IMatrix M({{-1., 0.000706767}, {-0.000706767, -1.}});
						
		vector<HSet2D> L(2);										// L[0], L[1] lie on the section
		L[0] = HSet2D (M*IVector ({-1.23094, 0.}) + p2, M, DVector({1.41278,7e-4}));	
		L[1] = HSet2D (M*IVector ({1.84699, 0.}) + p2, M, DVector({1.55949,7e-4}));
		
		///===================== choice of Procedure:  =====================
		int whichCase;
		
		do{
		cout << "===========================================================" << endl;
		cout << "|| Choose the number of procedure You wish to perform:" << endl;
		cout << "|| 1. Proof of a 3-periodic orbit's existence by INM" << endl;
		cout << "|| 2. Proof of the chain horizontal covering relations" << endl;
		cout << "|| 	N0 =Pc=> N1 =Pc=> N1 =Pc=> N0" << endl;
		cout << "|| other: exit" << endl;
		cout << "===========================================================" << endl;
		cout << "Procedure " ;
		cin >> whichCase;
		cout << endl;
		
		switch (whichCase)
			{
				///===================== Procedure 1:  =====================
				case 1:
				{
					contains3PeriodicPoint = roessler525.anyStationaryPoint(pCloseTo3periodic,3);				// now p3Periodic contains a 3-periodic point
					cout << "-----------------------------------------------------------" << endl;
					cout << "3-periodic orbit for Poincare map of the Roessler system with a = 5.25 contained in the set: " << endl;
					cout << "     p1 = " << contains3PeriodicPoint << endl;	
					for(int i=0;i<3;++i)																		// prints all orbit
						{
						contains3PeriodicPoint=roessler525.P2d(contains3PeriodicPoint);
						cout << "P^" << i+1 << "(p1) = " << contains3PeriodicPoint << endl;
						};			
					cout << "-----------------------------------------------------------" << endl << endl;	
					break;
				}
				
				///===================== Procedure 2:  =====================
				case 2:		
				{
					cout << "N0 covers N1? ... " << roessler525.covers2D(L[0],L[1],80) << endl;			// check if L[0] =P=> L[1]
																												// divided horizontally into 80 pieces
					cout << "N1 covers N1? ... " << roessler525.covers2D(L[1],L[1],140) << endl;		// if L[1] =P=> L[1]
																												// divided horizontally into 140 pieces
					cout << "N1 covers N0? ... " << roessler525.covers2D(L[1],L[0],100,2) << endl;		// if L[1] =P=> L[0]
					cout << "----------------------------------------" << endl << endl;									// divided into 100x2 pieces
					break;
				}
				
				///===================== Exit:  =====================
				default:
					cout << "Goodbye!" << endl;
					break;
				}
			}
			while (whichCase==1 or whichCase ==2);							// repeat the choice
			
	}
	catch(exception& e)
  	{
    		cout << "\n\nException caught: "<< e.what() << endl;
  	}
  return 0;
} 

/* FULL OUTPUT:
===========================================================
|| Roessler system, a = 5.25
===========================================================
|| Choose the number of procedure You wish to perform:
|| 1. Proof of a 3-periodic orbit's existence by INM
|| 2. Proof of the chain horizontal covering relations
|| 	N0 =Pc=> N1 =Pc=> N1 =Pc=> N0
|| other: exit
===========================================================
Procedure 1

A stationary point for P^3 in = {[-3.466415205012921, -3.466415205008745],[0.0346316054764013, 0.03463160547651117]} proven.
-----------------------------------------------------------
3-periodic orbit for Poincare map of the Roessler system with a = 5.25 contained in the set: 
     p1 = {[-3.466415205012921, -3.466415205008745],[0.0346316054764013, 0.03463160547651117]}
P^1(p1) = {[-6.264007533282313, -6.264007533274737],[0.03265435884611445, 0.03265435884612054]}
P^2(p1) = {[-9.748889918093207, -9.748889918088979],[0.03075287338063169, 0.03075287338070213]}
P^3(p1) = {[-3.466415205015202, -3.46641520500646],[0.03463160547637927, 0.0346316054765332]}
-----------------------------------------------------------

===========================================================
|| Choose the number of procedure You wish to perform:
|| 1. Proof of a 3-periodic orbit's existence by INM
|| 2. Proof of the chain horizontal covering relations
|| 	N0 =Pc=> N1 =Pc=> N1 =Pc=> N0
|| other: exit
===========================================================
Procedure 2

N0 covers N1? ... true
N1 covers N1? ... true
N1 covers N0? ... true
----------------------------------------

===========================================================
|| Choose the number of procedure You wish to perform:
|| 1. Proof of a 3-periodic orbit's existence by INM
|| 2. Proof of the chain horizontal covering relations
|| 	N0 =Pc=> N1 =Pc=> N1 =Pc=> N0
|| other: exit
===========================================================
Procedure 0

Goodbye!

real	0m4,368s
user	0m1,186s
sys	0m0,004s


*/

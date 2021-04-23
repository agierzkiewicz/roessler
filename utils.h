//////////////////////////////////////////////////////////////////////////////
///
///  @file utils.h
///  
///  @author ag  @date   Mar 5, 2020
//////////////////////////////////////////////////////////////////////////////


#ifndef _EXAMPLES_PROJECTSTARTER_UTILS_H_
#define _EXAMPLES_PROJECTSTARTER_UTILS_H_

#include "capd/intervals/lib.h"
#include "capd/capdlib.h"
#include <iostream>
#include "capd/covrel/HSet2D.h"
#include "capd/dynsys/DynSysMap.h"
#include "capd/covrel/HSetMD.h"
using namespace std;
using namespace capd;			

typedef capd::covrel::HSet2D<DMatrix, IMatrix> HSet2D;
typedef capd::dynsys::DynSysMap<IMap> DynSysMap;
typedef capd::covrel::HSetMD<DMatrix,IMatrix> HSet;
typedef capd::covrel::GridSet<IMatrix> GridSet;


//////////////////////////////////////////////////////////////////////////////
///	Some auxiliary functions
//////////////////////////////////////////////////////////////////////////////

// cuts the first coordinate of an interval vector
IVector cut(IVector x);	
DVector cut(DVector x);	
IMatrix cut(IMatrix M);	
			
//prepends an interval 2-vector (or 2x2 matrix) to 3-vector with zero
IVector expand(IVector x);	
DVector expand(DVector x);
IMatrix expand(IMatrix M);	

/////////////////////////////////////////////////////////////////////////
/// SecMap class
/////////////////////////////////////////////////////////////////////////

// 2D Poincare Map image of an interval vector on the x=0 section 
class SecMap	
{
public:
	IPoincareMap *P;
	SecMap(IPoincareMap &Q) {P=&Q;}	
	SecMap(){};						
	
	IVector image(const IVector &x,IMatrix &DP, const int iter = 1) const;
	IVector image(const IVector &x, const int iter = 1) const;
	IMatrix derivative(const IVector &x, const int iter = 1) const;
	IVector operator()(const IVector &x,IMatrix &DP, const int iter = 1) const {return image(x,DP,iter);}
	IVector operator()(const IVector &x, const int iter = 1) const {return image(x,iter);}
	IMatrix operator[](const IVector &x) const {return derivative(x);}
	
};

//////////////////////////////////////////////////////////////////////////////
///	system3d structure
//////////////////////////////////////////////////////////////////////////////

// stores the Roessler system
struct system3d
{
int order;
IMap vf;
ICoordinateSection section;
IOdeSolver solver;
IPoincareMap P;
SecMap P2d;
int iteration;

system3d(const interval &aa = interval(5.7)) :									// constructor, default parameter a=5.7
	order(20),
	vf("par:a,b; var:x,y,z; fun: -y - z, x + b*y, b + z*(x - a);"),					// vector field
	section(3,0),																	// x=0 section
	solver(vf,order),
	P(solver,section,poincare::MinusPlus),
	P2d(P)
	{
		//interval aa = interval(5.7);													// Roessler typical parameters
		interval bb = interval(0.2);
		vf.setParameter("a",aa);													
		vf.setParameter("b",bb); 
	}

IVector anyStationaryPoint(const IVector &x0, int iter =1);
bool covers2D(const HSet2D &hset1, const HSet2D &hset2, int howManyPiecesH=1, int howManyPiecesV=1, int iteration =1);
bool inside(const HSet2D &hset1, const HSet2D &hset2, int howManyPiecesH=1, int howManyPiecesV=1, int iteration = 1);
IVector whatIsNotMappedOutside(const HSet2D &hset1, int howManyPiecesH=1, int howManyPiecesV=1, int iteration = 1);
bool newtonDivided(const IVector &box, int howManyPiecesH, int howManyPiecesV, int iteration = 1);

};

#endif //_EXAMPLES_PROJECTSTARTER_UTILS_H_



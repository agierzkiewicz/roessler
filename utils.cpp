//////////////////////////////////////////////////////////////////////////////
///
///  @file utils.cpp
///  
///  @author ag  @date   Mar 5, 2020
//////////////////////////////////////////////////////////////////////////////

// Standard libraries and the CAPD library
#include <iostream>
#include "utils.h"
#include "capd/covrel/HSet2D.h"
#include "capd/dynsys/DynSysMap.h"
#include "capd/covrel/HSetMD.h"
#include "capd/intervals/lib.h"
#include "capd/capdlib.h"
//#include <stdexcept>
//#include <sstream>
//#include "capd/vectalg/iobject.hpp"
//#include "capd/dynset/C0DoubletonSet.h"
//#include "capd/geomset/CenteredDoubletonSet.hpp"
//#include "capd/matrixAlgorithms/floatMatrixAlgorithms.hpp"
using namespace capd;
using namespace capd::alglib;
using namespace capd::matrixAlgorithms;
using namespace std;

typedef capd::covrel::HSet2D<DMatrix, IMatrix> HSet2D;
typedef capd::dynsys::DynSysMap<IMap> DynSysMap;
typedef capd::covrel::HSetMD<DMatrix,IMatrix> HSet;
typedef capd::covrel::GridSet<IMatrix> GridSet;

//////////////////////////////////////////////////////////////////////////////
///	Some auxiliary functions
//////////////////////////////////////////////////////////////////////////////

// cut
// cuts the first coordinate of an interval 3-vector
IVector cut(IVector x)						
{
	if(x.dimension()==2) return x;
	IVector y(2);
	y[0]=x[1];
	y[1]=x[2];
	return y;
}

DVector cut(DVector x)						
{
	if(x.dimension()==2) return x;
	DVector y(2);
	y[0]=x[1];
	y[1]=x[2];
	return y;
}

// cuts an interval 3x3-matrix to 2x2
IMatrix cut(IMatrix M)						
{
	if(M.numberOfColumns()==2 and M.numberOfRows()==2) return M;
	IMatrix y(2,2);
	y[0][0]=M[1][1];
	y[0][1]=M[1][2];
	y[1][0]=M[2][1];
	y[1][1]=M[2][2];
	return y;
}

//expand
//prepends an interval 2-vector to 3-vector with zero 
IVector expand(IVector x)					
{	
	if(x.dimension()==3) return x;										
	IVector y({0.,0.,0.});
	y[1]=x[0];
	y[2]=x[1];
	return y;
}
DVector expand(DVector x)					
{	
	if(x.dimension()==3) return x;										
	DVector y({0.,0.,0.});
	y[1]=x[0];
	y[2]=x[1];
	return y;
}

//expands an interval 2x2 matrix to 3x3 with identity
IMatrix expand(IMatrix M)					
{		
	if(M.numberOfColumns()==3 and M.numberOfRows()==3) return M;									
	IMatrix y(3,3);
	y[1][1]=M[0][0];
	y[1][2]=M[0][1];
	y[2][1]=M[1][0];
	y[2][2]=M[1][1];
	y[0][1]=y[1][0]=y[2][0]=y[0][2]=interval(0.);
	y[0][0]=interval(1.);
	return y;
}

/////////////////////////////////////////////////////////////////////////
/// SecMap class
/////////////////////////////////////////////////////////////////////////

// image of the iteration iter, stores the derivative in DP
IVector SecMap::image(const IVector &x,IMatrix &DP, const int iter) const  
{
	IVector X({0.,0.,0.});
	X[1]=x[0];
	X[2]=x[1];
	
	C1Rect2Set Q(X);
	IMatrix DP3(3,3);
	
	IVector Y=(*P)(Q,DP3,iter);
	DP3=P->computeDP(Y,DP3,iter);
	
	IVector y(2);
	y[0]=Y[1];
	y[1]=Y[2];
	
	for(int i=0;i<2;++i)
	{
		for(int j=0;j<2;++j) DP[i][j]=DP3[i+1][j+1];
	}
	return y;
}

// image of the iteration iter only, no derivative calculation (faster)
IVector SecMap::image(const IVector &x, const int iter) const		
{
	IVector X({0.,0.,0.});
	X[1]=x[0];
	X[2]=x[1];
	
	C0Rect2Set Q(X);
	
	IVector Y=(*P)(Q,iter);
	
	IVector y(2);
	y[0]=Y[1];
	y[1]=Y[2];
	
	return y;
}

// derivative only
IMatrix SecMap::derivative(const IVector &x, const int iter) const		
{
	IMatrix DP(2,2);
	image(x,DP,iter);
	return DP;
}


//////////////////////////////////////////////////////////////////////////////
///	system3d structure methods
//////////////////////////////////////////////////////////////////////////////

// looks for a stationary point close to the starting point x0
IVector system3d::anyStationaryPoint(const IVector &x0, int iter)
{	  	
	IVector x = x0;
	IVector Px, xPrev, N(2);			
	IMatrix Id = IMatrix({{1.,0.},{0.,1.}});	// Identity 2x2 matrix	
	IMatrix DP(2,2);							
	int a=1;	
        try{
			do
			{
				x=midVector(x);					// Preliminary iteration to get a 'stable' point
				xPrev = x;
				Px=P2d(x,DP,iter);
				x=x-gauss(DP-Id,Px-x);
				++a;
			} while(not(subset(xPrev,x)) && a<50);

			/// if x is a good candidate, then we try a proper interval Newton method:
			if(a<50)
			{		
				const double delta = pow(10.0,-9);			// the size of the box X around x
				IVector X(2);
				X[0]=x[0]+delta*interval(-1.0,1.0);
				X[1]=x[1]+delta*interval(-1.0,1.0);	
				P2d(X,DP,iter);								// calculate DP = [dP^iter(X)]
				N = x-gauss(DP-Id,P2d(x,iter)-x);				// interval Newton
						
				if(subsetInterior(N,X))						// is N contained in the primary neighbourhood?;
				{
					cout << "A stationary point for P^" << iter << " in = " << N << " proven." << endl;
				}
			}
			else cout << "For the starting point " << x0 << " no success." << endl;
			} catch(exception& e){
              cout << "For the starting point " << x0 << " exception caught: " << e.what() << endl;
            }

return N;
}	

// checks if hset1 covers hset2 horizontally
bool system3d::covers2D(const HSet2D &hset1, const HSet2D &hset2, int howManyPiecesH, int howManyPiecesV, int iteration)
{	
	GridSet gridSet(2);																				// will contain the 2D grids
	HSet2D hset2Base(IVector({0.,0.}),IMatrix::Identity(2),hset2.get_r());
		
	IMatrix expandedGridSetCoordinateSystem = expand(hset1.coordinateSystem());						// hset1 coordinate system expanded to 3D
	IMatrix expandedHSet2InvCoordinateSystem = expand(hset2.invCoordinateSystem());					// hset2 inverse coordinate system expanded to 3D
	IVector expandedHSet2Center = expand(hset2.center());											// center of hset2 expanded to 3D
	IVector expandedGridSetBox;																		// expanded boxes of 2D grids
	
	C0Rect2Set Set3d(expand(hset1.center()));														
	IVector Pset2d(2);  
	interval rt(0.);
	
	bool liesAcross = 1;																			// check if hset1 is mapped across hset2	
		hset1.gridSet(gridSet,howManyPiecesH,howManyPiecesV);										// grid of the whole hset1
		expandedGridSetBox = expand(gridSet.box());
	
		for(auto i = gridSet.begin();i!=gridSet.end();++i)
		{
			Set3d = C0Rect2Set(expand(*i),expandedGridSetCoordinateSystem,expandedGridSetBox);		// tiny parallelogram in 3D
			for(int j=0;j<iteration-1;++j)															// iteration-1 of P		
				P(Set3d);
			Pset2d = cut(P(Set3d,expandedHSet2Center,expandedHSet2InvCoordinateSystem,rt));			// the last iterate converted to straight coordinates in 2D
			liesAcross = liesAcross && hset2Base.across(Pset2d);									// does the image lie between horizontal edges of hset2 or outside?
		}		
		if(!liesAcross) return false;																// does not lie across => no covering
	
	bool leftEdgeToLeft = 1;																		// check if left edge of hset1 is mapped to the left of hset2
		hset1.gridLeftEdge(gridSet,howManyPiecesV);	  												// grid of left edge of hset1
		//IMatrix expandedGridSetCoordinateSystem = expand(gridSet.coordinateSystem());
		expandedGridSetBox = expand(gridSet.box());
		
		for(auto i = gridSet.begin();i!=gridSet.end();++i)
		{
			Set3d = C0Rect2Set(expand(*i),expandedGridSetCoordinateSystem,expandedGridSetBox);
			for(int j=0;j<iteration-1;++j)
				P(Set3d);
			Pset2d = cut(P(Set3d,expandedHSet2Center,expandedHSet2InvCoordinateSystem,rt));
			leftEdgeToLeft = leftEdgeToLeft && hset2Base.onLeft(Pset2d);							// does the image lie to the left of hset2?
		}
	
	bool rightEdgeToRight = 1;																		// check if right edge of hset1 is mapped to the right of hset2
		hset1.gridRightEdge(gridSet,howManyPiecesV);	  											// grid of right edge of hset1	  
		expandedGridSetBox = expand(gridSet.box());
	
		for(auto i = gridSet.begin();i!=gridSet.end();++i)
		{
			Set3d = C0Rect2Set(expand(*i),expandedGridSetCoordinateSystem,expandedGridSetBox);
			for(int j=0;j<iteration-1;++j)
				P(Set3d);
			Pset2d = cut(P(Set3d,expandedHSet2Center,expandedHSet2InvCoordinateSystem,rt));
			rightEdgeToRight = rightEdgeToRight && hset2Base.onRight(Pset2d);						// does the image lie to the right of hset2?
		}
	
	bool leftEdgeToRight = 1;																		// check if left edge of hset1 is mapped to the right of hset2
		hset1.gridLeftEdge(gridSet,howManyPiecesV);	
		expandedGridSetBox = expand(gridSet.box());
		
		for(auto i = gridSet.begin();i!=gridSet.end();++i)
		{
			Set3d = C0Rect2Set(expand(*i),expandedGridSetCoordinateSystem,expandedGridSetBox);
			for(int j=0;j<iteration-1;++j)
				P(Set3d);
			Pset2d = cut(P(Set3d,expandedHSet2Center,expandedHSet2InvCoordinateSystem,rt));
			leftEdgeToRight = leftEdgeToRight && hset2Base.onRight(Pset2d);
		}
	
	bool rightEdgeToLeft = 1;																		// check if right edge of hset1 is mapped to the left of hset2
		hset1.gridRightEdge(gridSet,howManyPiecesV);	  
		expandedGridSetBox = expand(gridSet.box());
	
		for(auto i = gridSet.begin();i!=gridSet.end();++i)
		{
			Set3d = C0Rect2Set(expand(*i),expandedGridSetCoordinateSystem,expandedGridSetBox);
			for(int j=0;j<iteration-1;++j)
				P(Set3d);
			Pset2d = cut(P(Set3d,expandedHSet2Center,expandedHSet2InvCoordinateSystem,rt));
			rightEdgeToLeft = rightEdgeToLeft && hset2Base.onLeft(Pset2d);
		}
		
	return (liesAcross and ((leftEdgeToLeft and rightEdgeToRight) or (leftEdgeToRight and rightEdgeToLeft)));	// 2D covering condition
}

// checks if the image of hset1 lies inside hset2
bool system3d::inside(const HSet2D &hset1, const HSet2D &hset2, int howManyPiecesH, int howManyPiecesV, int iteration) 
{
	GridSet gridSet(2);																		// stores grid
	hset1.gridSet(gridSet,howManyPiecesH,howManyPiecesV);		
	HSet2D hset2Base(IVector({0.,0.}),IMatrix::Identity(2),hset2.get_r());
  
	C0Rect2Set Set3d(expand(hset1.center()));
	IVector Pset2d(2);  
	interval rt(0.);
	IMatrix expandedGridSetCoordinateSystem = expand(gridSet.coordinateSystem());			// expanded to 3d, where P is defined
	IVector expandedGridSetBox = expand(gridSet.box());
	IMatrix expandedHSet2InvCoordinateSystem = expand(hset2.invCoordinateSystem());
	IVector expandedHSet2Center = expand(hset2.center());
	
	bool liesInside = 1;
	for(auto i = gridSet.begin();i!=gridSet.end();++i)
	{
		Set3d = C0Rect2Set(expand(*i),expandedGridSetCoordinateSystem,expandedGridSetBox);
		for(int j=0;j<iteration-1;++j)
			P(Set3d);
		Pset2d = cut(P(Set3d,expandedHSet2Center,expandedHSet2InvCoordinateSystem,rt));		// the image of the small box
		liesInside=liesInside && hset2Base.inside(Pset2d);									// do all lie inside?
		
		if(!liesInside) 
			{
				cout << "Does not lie inside, for this tiny box: " << Pset2d << " does not lie inside the box " << hset2Base << endl;
				return 0;
			}
	}
    return liesInside;
}

// checks what part of the image of hset1 does not lie outside itself
IVector system3d::whatIsNotMappedOutside(const HSet2D &hset1, int howManyPiecesH, int howManyPiecesV, int iteration) 
{
	//Grid
	GridSet gridSet(2);																		// will contain the 2D grid
	hset1.gridSet(gridSet,howManyPiecesH,howManyPiecesV);		
	DVector r({(hset1.get_r())[0]/howManyPiecesH,(hset1.get_r())[1]/howManyPiecesV});
	HSet2D hset1TinyBase(IVector({0.,0.}),IMatrix::Identity(2),r);
  
	C0Rect2Set Set3d(expand(hset1.center()));
	IVector Pset2d, problematic(2);  
	interval rt(0.);
	IMatrix expandedGridSetCoordinateSystem = expand(gridSet.coordinateSystem());			// expanded to 3d, where P is defined
	IVector expandedGridSetBox = expand(gridSet.box());
	IMatrix expandedHSet1InvCoordinateSystem = expand(hset1.invCoordinateSystem());
	IVector expandedHSet1Center = expand(hset1.center());

	
	int a = 1;
	for(auto i = gridSet.begin();i!=gridSet.end();++i)
	{		
		hset1TinyBase = HSet2D(hset1.invCoordinateSystem()*((*i)-hset1.center()),IMatrix::Identity(2),r);
		Set3d = C0Rect2Set(expand(*i),expandedGridSetCoordinateSystem,expandedGridSetBox);
		for(int j=0;j<iteration-1;++j)
			P(Set3d);
		Pset2d = cut(P(Set3d,expandedHSet1Center,expandedHSet1InvCoordinateSystem,rt));		//the image of the small box
		
		if(!hset1TinyBase.outside(Pset2d) and a==1) 										// if mapped outside for first time
			{
				problematic = hset1TinyBase.box()+hset1TinyBase.center();					// stores in problematic
				++a;
			}
		else if(!hset1TinyBase.outside(Pset2d))												// if mapped outside next time
			{
				problematic = intervalHull(problematic, hset1TinyBase.box()+hset1TinyBase.center());
				++a;																		// joins to problematic
			};
	};
  
    return hset1.coordinateSystem()*problematic + hset1.center();							// problematic in affine coordinate system
}

// checks if box contains a unique stationary point via INM
bool system3d::newtonDivided(const IVector &box, int howManyPiecesH, int howManyPiecesV, int iteration) 
{	
	IVector X, Px, N, x(2);  
	double dx = (box[0].rightBound()-box[0].leftBound())/howManyPiecesH;
	double dy = (box[1].rightBound()-box[1].leftBound())/howManyPiecesV;
	
	bool newtonDetermines = 1;
	int a = 0;										// no of periodic points
	IMatrix m1, DP(2,2);
	IMatrix Id = IMatrix({{1.,0.},{0.,1.}});		// Identity 2x2 matrix
	
	for(int i = 0; i < howManyPiecesH ; ++i)
	{	
		for(int j = 0; j < howManyPiecesV ; ++j)
		{
			X = IVector({interval((box[0]).leftBound()+i*dx,(box[0]).leftBound()+(i+1)*dx),interval((box[1]).leftBound()+j*dy,(box[1]).leftBound()+(j+1)*dy)});	
			x= midVector(X);	
			P2d(X,DP,iteration);
			Px= P2d(x,iteration);
			N=x-gauss(DP-Id,Px-x);					// interval Newton
						
			if(subsetInterior(N,X))					// if N is contained in the primary neighbourhood
			{
				cout << "The box " << X << " contains a unique stationary point for P^" << iteration << endl;
				cout << " in = " << N << endl;					
				++a;								// count stationary points
			}
			else if(intersectionIsEmpty(N,X))
			{
				//cout << "no stationary points here! " << endl;
			}
			else
			{	
				//cout << "we know nothing here!" << endl;
				newtonDetermines = 0;
			};
		};		
	}  
    return (newtonDetermines and (a==1));
}




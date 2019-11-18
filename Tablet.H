#ifndef __LEVODOPA_TABLET__
#define __LEVODOPA_TABLET__

#include "fvCFD.H"
namespace Foam { 
class Tablet  {
	
	/*======================================
		     	public
	======================================*/
	public:

	Tablet(Time& iRunTime,
	       fvMesh& iStomach,    
	       fvMesh& iSmallIntestine,
	       dictionary& tabletDict,
	       dimensionedScalar   istomachVolume,
	       volScalarField& iIntestineV, // intestinal volume
	       volVectorField& iIntestineU); 

	// Set the dosage in the stomach and the intestine. 
	void setInitialConditions(volScalarField& stomach,volScalarField& intestine); 

	// get the stomach source. 
	volScalarField& getStomachSource();
	volScalarField& getIntestineSource();


	/*======================================
		     	private
	======================================*/
	private:


	Time  & __RunTime; 
	fvMesh& __Stomach; 
	fvMesh& __SmallIntestine;

	dimensionedScalar  __LocationDistance;
	label   __LocationLabel; // -1 = stomach, ID = Small Intestine Cell. 

	const volScalarField& __IntestineV; // intestinal volume
	const volVectorField& __IntestineU; 
	volScalarField __IntestineSource; 
	volScalarField __StomachSource; 
 
	dimensionedScalar&  __stomachVolume;
	dimensionedScalar   __tabletSize; // in kg. 

	dimensionedScalar   __Stomach_Burst; 
	dimensionedScalar   __Stomach_ResidenceTime; 

	dimensionedScalar   __StartErode; // The time the tablet starts the erosion. 

	dimensionedScalar   __ErosionRate; 

};
}
#endif 

#include "Tablet.H"


Foam::Tablet::Tablet(Time& iRunTime,
	       fvMesh& iStomach, 
	       fvMesh& iSmallIntestine,
	       dictionary& tabletDict,
	       dimensionedScalar   istomachVolume,
	       volScalarField& iIntestineV, // intestinal volume
	       volVectorField& iIntestineU) :
		__RunTime(iRunTime),
		__Stomach(iStomach), 
		__SmallIntestine(iSmallIntestine),
		__IntestineV(iIntestineV),
		__IntestineU(iIntestineU),
		__IntestineSource(
			IOobject
			(
			    "TabletIntestine",    // dictionary name
			    __RunTime.constant(),     // dict is found in "constant"
			    __SmallIntestine,                   // registry for the dict
			    IOobject::NO_READ,    // must exist, otherwise failure
			    IOobject::NO_WRITE      // dict is only read by the solver
			),
			__SmallIntestine,
		 	dimensionedScalar("zero0",dimMass/dimTime/dimVolume,scalar(0))
		), 
		__StomachSource(
			IOobject
			(
			    "TabletStomach",    // dictionary name
			    __RunTime.constant(),     // dict is found in "constant"
			    __Stomach,                   // registry for the dict
			    IOobject::NO_READ,    // must exist, otherwise failure
			    IOobject::NO_WRITE      // dict is only read by the solver
			),
			__Stomach,
		 	dimensionedScalar("zero1",dimMass/dimTime/dimVolume,scalar(0))
		),
		__stomachVolume(istomachVolume)
	{

		__LocationLabel    	= tabletDict.lookupOrDefault("TabletLocation",-1); // start in the stomach. 
		__tabletSize.dimensions().reset(dimMass);
		__Stomach_Burst.dimensions().reset(dimMass);
		__Stomach_ResidenceTime.dimensions().reset(dimTime);
		__StartErode.dimensions().reset(dimTime);
		__LocationDistance.dimensions().reset(dimLength);
		__ErosionRate.dimensions().reset(dimMass/dimTime/dimMass);
		__Stomach_ErosionRate.dimensions().reset(dimMass/dimTime/dimMass);

		__tabletSize  		= dimensionedScalar(tabletDict.lookup("tabletSize")); 
		__Stomach_Burst 	= dimensionedScalar(tabletDict.subDict("Stomach").lookup("Burst") );
		__Stomach_ResidenceTime = dimensionedScalar(tabletDict.subDict("Stomach").lookup("Residence") );
		__StartErode		= dimensionedScalar(tabletDict.lookup("ErosionStart")); 
		__ErosionRate		= dimensionedScalar(tabletDict.lookup("ErosionRate")); 

		__Stomach_ErosionRate   = dimensionedScalar(tabletDict.subDict("Stomach").lookupOrDefault("ErosionRate",__ErosionRate)); 
		

	}




// Set the dosage in the stomach and the intestine. 
void Foam::Tablet::setInitialConditions(volScalarField& stomach,volScalarField& intestine) {
	stomach[0] += (__Stomach_Burst/__stomachVolume).value();
	__tabletSize -= __Stomach_Burst;

}

// get the stomach source. 
volScalarField& Foam::Tablet::getStomachSource() {
	__StomachSource[0] = 0;
	
	if (__LocationLabel == -1) { 

		if (__RunTime.time() > __StartErode) { 
			dimensionedScalar eroded = __tabletSize*__Stomach_ErosionRate*__RunTime.deltaT();
			__StomachSource[0] = (eroded/__stomachVolume).value();
		

			__tabletSize -= eroded;
			if (__tabletSize.value() < 0) {
					__tabletSize.value() = 0;
			}
		}
	
		if (__RunTime.time() > __Stomach_ResidenceTime) { 
			__LocationLabel = 0; 
			__LocationDistance = dimensionedScalar("LocationDistance",dimLength,scalar(0)); 
		}
	}	
	return __StomachSource;
}

volScalarField& Foam::Tablet::getIntestineSource() { 
	forAll(__IntestineSource,celli) { 
			__IntestineSource[celli] = 0;
	}

	if (__LocationLabel > -1) { 
		if (__RunTime.time() > __StartErode) { 

			dimensionedScalar eroded = __tabletSize*__ErosionRate*__RunTime.deltaT();
			__IntestineSource[__LocationLabel]  += (eroded/__IntestineV[__LocationLabel] ).value();


			__tabletSize -= eroded;
			if (__tabletSize.value() < 0) {
				__tabletSize.value() = 0;
			}
		}

		// advance the tablet.
		__LocationDistance.value() += (__IntestineU[__LocationLabel].component(0)*__RunTime.deltaT()).value();
		// check if passed the boundary, and if you did update the label (and make 

		if (__LocationDistance.value() > __SmallIntestine.Cf()[__LocationLabel+1].component(0)) { 
			__LocationLabel++;
		}

		if (__LocationDistance.value() > __SmallIntestine.bounds().span().component(0)) { 
			__LocationLabel = -2;
		}
		Info << " ------- " << __IntestineSource[__LocationLabel] << " -- " << __LocationLabel<< " ---- " << __tabletSize << endl;

	} //.. if. 	


	return __IntestineSource;
}


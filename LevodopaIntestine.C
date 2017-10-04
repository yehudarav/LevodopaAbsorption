/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
	LevedopaIntestine

Description
 
=======

	A basic solver for stomach-small intestine processes of levodopa absorption. 


	The intestinal concentration is A/(CrossSection*dh) 

	

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "radiationModel.H"
#include "turbulentTransportModel.H"

#include "fvIOoptionList.H"
#include "pimpleControl.H"
#include "interpolation.H"
#include "Random.H"
#include "meshSearch.H"

#include <sstream>
#include <fstream>



//#include "stdlib.h"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    Info << " \n\nIntestinalLevodopa  :  0.0.1 (b)" << endl; 
    Info << " -------------------------- " << endl; 

     const scalar pi = 3.1428;
	
    // --------------------------------------------------- Stomach ---------------------------------------

    // Create the Mesh
    Foam::fvMesh Stomach 
    (
        Foam::IOobject
        (
            "Stomach",
            "0",
            runTime,
            Foam::IOobject::MUST_READ
        )
    );

    // Read the constants 
    IOdictionary StomachPropertiesDict
    (
        IOobject
        (
            "StomachProperties",    // dictionary name
            runTime.constant(),     // dict is found in "constant"
            Stomach,                   // registry for the dict
            IOobject::MUST_READ,    // must exist, otherwise failure
            IOobject::NO_WRITE      // dict is only read by the solver
        )
    ); 
    dimensionedScalar stomachemptying(StomachPropertiesDict.lookup("stomachemptying"));    
    dimensionedScalar stomachvolume(StomachPropertiesDict.lookup("stomachvolume"));    

    // Create the variables. 
    volScalarField LevodopaStomach 
    (
        IOobject
        (
            "LevodopaStomach",
            runTime.timeName(),
            Stomach,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        Stomach
    );
   
    volScalarField TotalLevodopaStomachEmptying(IOobject("TotalLevodopaStomachEmptying",runTime.timeName(),Stomach,IOobject::NO_READ,IOobject::NO_WRITE),
					   Stomach,dimensionedScalar("zero",dimMass,scalar(0)) );

    volScalarField LevodopaStomachEmptying(IOobject("LevodopaStomachEmptying",runTime.timeName(),Stomach,IOobject::NO_READ,IOobject::NO_WRITE),
					   Stomach,dimensionedScalar("zero",dimMass/dimTime,scalar(0)) );

    // --------------------------------------------------- Small Intestine ---------------------------------------

    // Create the Mesh
    Foam::fvMesh SmallIntestine 
    (
        Foam::IOobject
        (
            "SmallIntestine",
            "0",
            runTime,
            Foam::IOobject::MUST_READ
        )
    );


    // Read the constants 
    IOdictionary SmallIntestinePropertiesDict
    (
        IOobject
        (
            "SmallIntestineProperties",    // dictionary name
            runTime.constant(),     // dict is found in "constant"
            SmallIntestine,                   // registry for the dict
            IOobject::MUST_READ,    // must exist, otherwise failure
            IOobject::NO_WRITE      // dict is only read by the solver
        )
    ); 

    // Read the constants 
    IOdictionary SmallIntestineLevodopaPropertiesDict
    (
        IOobject
        (
            "SmallIntestineLevodopaProperties",    // dictionary name
            runTime.constant(),     // dict is found in "constant"
            SmallIntestine,                   // registry for the dict
            IOobject::MUST_READ,    // must exist, otherwise failure
            IOobject::NO_WRITE      // dict is only read by the solver
        )
    ); 

    
    // Create the variables. 
    volScalarField LevodopaSmallIntestine 
    (
        IOobject
        (
            "LevodopaSmallIntestine",
            runTime.timeName(),
            SmallIntestine,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        SmallIntestine
    );

   
    dimensionedVector SI_U(SmallIntestinePropertiesDict.lookup("SI_U"));
    dimensionedScalar SI_Radius(SmallIntestinePropertiesDict.lookup("SI_Radius"));	
    dimensionedScalar SI_L(SmallIntestinePropertiesDict.lookup("SI_Length"));
    dimensionedScalar SI_Duodenum(SmallIntestinePropertiesDict.lookup("SI_DuodenumLength"));
    dimensionedScalar SI_Jejunum(SmallIntestinePropertiesDict.lookup("SI_JejunumLength"));

    dimensionedScalar SIemptying(SmallIntestinePropertiesDict.lookup("Emptying"));    
	
    volVectorField SmallIntestineVelocity 
    (
        IOobject
        (
            "SmallIntestineVelocity",
            runTime.timeName(),
            SmallIntestine,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        SmallIntestine,
	SI_U
    );

	

    volTensorField SmallIntestineDispersion 
    (
        IOobject
        (
            "SmallIntestineDispersion",
            runTime.timeName(),
            SmallIntestine,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        SmallIntestine,
	dimensionedTensor(SmallIntestinePropertiesDict.lookup("SI_D"))
    );

   // Setting the last cell to be 0, 
   // the intestine should be emptied (velocity = velocity of the rest of the intestine) according to the intestinal emptying parameter. 
   //SmallIntestineVelocity[SmallIntestine.C().size()-1].component(0) = 0;
   //SmallIntestineVelocity[SmallIntestine.C().size()-2].component(0) = 0;
   // need to set the boundary condition to be 0. not the center of the cell. 

    surfaceScalarField SI_phi(
        IOobject
        (
            "SmallIntestineVelocity",
            runTime.timeName(),
            SmallIntestine,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	linearInterpolate(SmallIntestineVelocity) & SmallIntestine.Sf()
    );

    label IleoCecalpatchID 		= SmallIntestine.boundaryMesh().findPatchID("IleoCecal"); 
    const polyPatch& IleoCecalcPatch 	= SmallIntestine.boundaryMesh()[IleoCecalpatchID]; 


	forAll(IleoCecalcPatch, facei) 
	{ 
		SI_phi.boundaryField()[IleoCecalpatchID][facei] =0; 
	} 


    label PyloruspatchID 		= SmallIntestine.boundaryMesh().findPatchID("Pylorus"); 
    const polyPatch& PylorusPatch 	= SmallIntestine.boundaryMesh()[PyloruspatchID]; 


	forAll(PylorusPatch, facei) 
	{ 
		SI_phi.boundaryField()[PyloruspatchID][facei] =0; 
	} 


    // Amplification (without the microvilli) per unit area of intestine. 
    volScalarField SI_SurfaceArea(IOobject("SI_SurfaceArea",runTime.timeName(),SmallIntestine,IOobject::NO_READ,IOobject::NO_WRITE),
					   SmallIntestine,dimensionedScalar("one",dimLength,scalar(1)) );
    

    volScalarField SI_dh(IOobject("SI_dH",runTime.timeName(),SmallIntestine,IOobject::NO_READ,IOobject::NO_WRITE),
					   SmallIntestine,dimensionedScalar("one",dimLength,scalar(1)) );


    volScalarField SI_FoldingAmplification(IOobject("SI_SurfaceArea",runTime.timeName(),SmallIntestine,IOobject::NO_READ,IOobject::NO_WRITE),
					   SmallIntestine,dimensionedScalar("one",dimless,scalar(1)) );


    volScalarField SI_V = pi*SI_Radius*SI_Radius*SI_dh;


    // folding = 3 from 10 to 50. 
    // decrease linearly from 3 to 1 in the range 50->282. 
    const volScalarField& X = SmallIntestine.C().component(0); 
    forAll(SI_SurfaceArea,cellID) { 
	// 		Krecking folding contribution. 
			scalar folding_amplification;

			if (X[cellID] < SI_Duodenum.value()/2) { 
				folding_amplification = 1;
			} else if (X[cellID] > SI_Duodenum.value()/2 && X[cellID] < SI_Jejunum.value()/2) { 
				folding_amplification = 3;
			} else { 
				folding_amplification = 3 - 2*(X[cellID] - SI_Jejunum.value()/2) / (SI_L.value()-SI_Jejunum.value()/2); 
			}

			SI_FoldingAmplification[cellID]=folding_amplification;

	// 		Villi contribution. 
			scalar villiradius = 50*1e-6; // 50 micron.
			scalar villiheight = (800 - 300*X[cellID]/SI_L.value())*1e-6; // 800 to 500 micron.  
			
			scalar villi_tip_area  = 2*pi*villiradius*villiradius; // area of a half sphere. 
			scalar villi_body_area = 2*pi*villiradius*villiheight;
			
			scalar villi_density   = 25e6; // 25 mm^-2 = 25e6m**2. 
			scalar villi_amplification = villi_density*(villi_tip_area+villi_body_area);

			SI_SurfaceArea[cellID] 		= 2*pi*SI_Radius.value()*folding_amplification*villi_amplification;
			SI_dh[cellID]  = SmallIntestine.V()[cellID];
    }

	

    volScalarField SI_SurfaceArea_Microvilli(IOobject("SI_SurfaceArea",runTime.timeName(),SmallIntestine,IOobject::NO_READ,IOobject::NO_WRITE),
					   SmallIntestine,dimensionedScalar("one",dimless,scalar(25)) );   

    dimensionedScalar SI_CrossSection(SmallIntestinePropertiesDict.lookup("SI_CrossSection"));

    // Levodopa properties. 
    dimensionedScalar Levodopa_LNAA_Vmax(SmallIntestineLevodopaPropertiesDict.lookup("LNAA_Vmax"));
    dimensionedScalar Levodopa_LNAA_Km(SmallIntestineLevodopaPropertiesDict.lookup("LNAA_Km"));

    volScalarField SI_Emptying(IOobject("SI_Emptying",runTime.timeName(),SmallIntestine,IOobject::NO_READ,IOobject::AUTO_WRITE),
					     SmallIntestine,dimensionedScalar("one",dimMass/dimTime/dimVolume,scalar(0)) );

    volScalarField SI_LevodopaAbsorption(IOobject("SI_LevodopaAbsorption",runTime.timeName(),SmallIntestine,IOobject::NO_READ,IOobject::AUTO_WRITE),
					     SmallIntestine,dimensionedScalar("one",dimMass/dimTime,scalar(0)) );


    volScalarField SI_LevodopaAbsorption_Coeff(IOobject("SI_LevodopaAbsorption",runTime.timeName(),SmallIntestine,IOobject::NO_READ,IOobject::NO_WRITE),
					     SmallIntestine,dimensionedScalar("one",dimless/dimTime,scalar(0)) );

  
    volScalarField SI_LevodopaTotalAbsorption(IOobject("SI_LevodopaTotalAbsorption",runTime.timeName(),SmallIntestine,IOobject::NO_READ,IOobject::NO_WRITE),
					     SmallIntestine,dimensionedScalar("one",dimMass,scalar(0)) );


    // --------------------------------------------------- Epithelium ---------------------------------------
    // Create the Mesh
    Foam::fvMesh Epithelium
    (
        Foam::IOobject
        (
            "Epithelium",
            "0",
            runTime,
            Foam::IOobject::MUST_READ
        )
    );

    // Read the constants 
    IOdictionary EpitheliumLevodopaPropertiesDict
    (
        IOobject
        (
            "EpitheliumLevodopaProperties",    // dictionary name
            runTime.constant(),     // dict is found in "constant"
            Epithelium,                   // registry for the dict
            IOobject::MUST_READ,    // must exist, otherwise failure
            IOobject::NO_WRITE      // dict is only read by the solver
        )
    ); 

    
    // Create the variables. 
    volScalarField LevodopaEpithelium 
    (
        IOobject
        (
            "LevodopaEpithelium",
            runTime.timeName(),
            Epithelium,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        Epithelium
    );

    dimensionedScalar Epithelium_width(SmallIntestinePropertiesDict.lookup("SI_Epithelium_H"));

    volScalarField Epithelium_CrossSection(IOobject("Epithelium_Volume",runTime.timeName(),Epithelium,IOobject::NO_READ,IOobject::NO_WRITE),
				     Epithelium,
				     dimensionedScalar("one",dimArea,scalar(1)) );

    volScalarField Epithelium_Volume(IOobject("Epithelium_Volume",runTime.timeName(),Epithelium,IOobject::NO_READ,IOobject::NO_WRITE),
				     Epithelium,
				     dimensionedScalar("one",dimVolume,scalar(1)) );


    volScalarField Epithelium_SurfaceArea(IOobject("Epithelium_SurfaceArea",runTime.timeName(),Epithelium,IOobject::NO_READ,IOobject::NO_WRITE),
				     Epithelium,
				     dimensionedScalar("one",dimLength,scalar(1)) );


    forAll(Epithelium_Volume,CellID) { 
		Epithelium_CrossSection[CellID]  = 2*pi*SI_Radius.value()*Epithelium_width.value()*SI_FoldingAmplification[CellID];
		Epithelium_Volume[CellID] 	 = Epithelium_CrossSection[CellID]*SI_dh[CellID];
		Epithelium_SurfaceArea[CellID]	 = SI_SurfaceArea[CellID]; // transfer mesh. 
    }






    volScalarField Epithelium_LevodopaAbsorption_From_SI(IOobject("Epithelium_LevodopaAbsorption_From_SI",runTime.timeName(),Epithelium,IOobject::NO_READ,IOobject::NO_WRITE),
					     		 Epithelium,
							 dimensionedScalar("one",dimMass/(dimVolume*dimTime),scalar(0)) );

    volScalarField Epithelium_TotalLevodopaAbsorption_From_SI(IOobject("Epithelium_TotalLevodopaAbsorption_From_SI",runTime.timeName(),Epithelium,IOobject::NO_READ,IOobject::NO_WRITE),
					     		 Epithelium,
							 dimensionedScalar("one",dimMass,scalar(0)) );

    dimensionedScalar Levodopa_COMT_Vmax(EpitheliumLevodopaPropertiesDict.lookup("COMT_Vmax"));
    dimensionedScalar Levodopa_COMT_Km  (EpitheliumLevodopaPropertiesDict.lookup("COMT_Km"  ));

    volScalarField Epithelium_COMT_TotalMetabolism(IOobject("Epithelium_COMT_Metabolism",runTime.timeName(),Epithelium,IOobject::NO_READ,IOobject::NO_WRITE),
					     Epithelium,dimensionedScalar("one",dimMass,scalar(0)) );


    volScalarField Epithelium_Flux_To_Body(IOobject("Epithelium_Flux_To_Body",runTime.timeName(),Epithelium,IOobject::NO_READ,IOobject::NO_WRITE),
					     Epithelium,dimensionedScalar("one",dimMass/dimTime,scalar(0)) );

	
    volScalarField Epithelium_TotalFlux_To_Body(IOobject("Epithelium_Flux_To_Body",runTime.timeName(),Epithelium,IOobject::NO_READ,IOobject::NO_WRITE),
					     Epithelium,dimensionedScalar("one",dimMass,scalar(0)) );



    // --------------------------------------------------- Body ---------------------------------------
    // Create the Mesh
    Foam::fvMesh Body
    (
        Foam::IOobject
        (
            "Body",
            "0",
            runTime,
            Foam::IOobject::MUST_READ
        )
    );

    // Read the constants 
    IOdictionary BodyLevodopaPropertiesDict
    (
        IOobject
        (
            "BodyLevodopaProperties",    // dictionary name
            runTime.constant(),     // dict is found in "constant"
            Body,                   // registry for the dict
            IOobject::MUST_READ,    // must exist, otherwise failure
            IOobject::NO_WRITE      // dict is only read by the solver
        )
    ); 

    
    // Create the variables. 
    volScalarField LevodopaBody
    (
        IOobject
        (
            "LevodopaBody",
            runTime.timeName(),
            Body,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        Body
    );

    dimensionedScalar bodyelimination(BodyLevodopaPropertiesDict.lookup("elimination"));    
    dimensionedScalar bodyvolume(BodyLevodopaPropertiesDict.lookup("volume"));    

    volScalarField Body_Flux_To_Body(IOobject("Body_Flux_To_Body",runTime.timeName(),Body,IOobject::NO_READ,IOobject::NO_WRITE),
					     Body,dimensionedScalar("one",dimMass/(dimVolume*dimTime),scalar(0)) );

    volScalarField TotalElimination(IOobject("TotalElimination",runTime.timeName(),Body,IOobject::NO_READ,IOobject::NO_WRITE),
					     Body,dimensionedScalar("one",dimMass,scalar(0)) );


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


    Info<< "\nStarting time loop\n" << endl;


    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        //#include "readTimeControls.H"
        //#include "CourantNo.H"   
        //#include "setDeltaT.H"
	
	{ 	// Stomach
		fvScalarMatrix LevodopaStomachEqn (
			fvm::ddt(LevodopaStomach) == fvm::Sp(-stomachemptying,LevodopaStomach)			
		); 	
		
		LevodopaStomachEqn.solve();
		LevodopaStomachEmptying  = fvc::Sp(stomachemptying,LevodopaStomach)*stomachvolume; 
		TotalLevodopaStomachEmptying += LevodopaStomachEmptying*runTime.deltaT();
		
		Info << endl << "Stomach" << endl;
		Info << "--------" << endl; 
		Info << "In Stomach " << (LevodopaStomach*stomachvolume)->internalField()[0] << endl;
		Info << "Emptied " << TotalLevodopaStomachEmptying.internalField()[0] << endl;
		Info << "In Stomach + Emptied = " << (LevodopaStomach*stomachvolume)->internalField()[0] << " + " 
						      << TotalLevodopaStomachEmptying.internalField()[0] << " = " 
						      << (LevodopaStomach*stomachvolume +TotalLevodopaStomachEmptying)->internalField()[0]
						      << endl;


	} //.. stomach 


	{ // Small intestine

		// Build the Sp coefficient. 
		// remember that SI_SurfaceArea is for unit length (m**2/m=m). and so we have to multiply it to get the total amount absorbed. 
		//			         m          *      []                 *      kg/s/m**2   /      (kg/m3)                              /	   m**2	       = kg*m2/(m2*s*kg) = 1/s. 
		//													     =  m**2/s               /	   m**2	       = 1/s.               
		SI_LevodopaAbsorption_Coeff = SI_SurfaceArea*SI_SurfaceArea_Microvilli*Levodopa_LNAA_Vmax/(Levodopa_LNAA_Km + LevodopaSmallIntestine)/SI_CrossSection;
	
		// SI_LevodopaAbsorption_Coeff*LevodopaSmallIntestine = g/(s*m**3) that is the amount absorbed per unit length. 
		//     --- divide by pi*R^2 to get the concentration per cell. 
		// 

		// Intestinal phase
	
		volScalarField TimeStepAbs = SI_LevodopaAbsorption_Coeff*LevodopaSmallIntestine;
		
		// Empty to the first cell. 
		SI_Emptying[0]  = LevodopaStomachEmptying[0]/SI_V[0];

		fvScalarMatrix LevodopaIntestineEqn (
			fvm::ddt(LevodopaSmallIntestine)  
			+ fvm::div(SI_phi,LevodopaSmallIntestine)  
			- fvm::laplacian(SmallIntestineDispersion,LevodopaSmallIntestine) 
					== 
			SI_Emptying
			-TimeStepAbs
		      
		); 	

		LevodopaIntestineEqn.solve();

		SI_LevodopaAbsorption       = TimeStepAbs*SI_V; // g/(s*m**3)*m**3 = g/s
		SI_LevodopaTotalAbsorption += SI_LevodopaAbsorption*runTime.deltaT(); // g/s *s = g. 


		Info << endl << "Small Intestine" << endl;
		Info <<         "---------------" << endl;
		Info << "Intestinal Levodopa " << sum(LevodopaSmallIntestine*SI_V) << endl;
		Info << "Absorption " << sum(SI_LevodopaTotalAbsorption) << endl;

		Info << "Total " << sum(SI_LevodopaTotalAbsorption)+sum(LevodopaSmallIntestine*SI_V) << endl;


		if (runTime.time() > SIemptying) { 

			forAll(IleoCecalcPatch, facei) 
			{ 
				SI_phi = linearInterpolate(SmallIntestineVelocity) & SmallIntestine.Sf();
			} 

		}

		

	} //.. small intestine. 

	

	{ // Epithelium

		// Build the Sp coefficient. 
		//			            		        m           *      kg/s/m2     /	  m2		 /              (kg/m3)                    * kg/m3             = kg/m3/s. 
		volScalarField Epithelium_Basal_Absorption = (Epithelium_SurfaceArea*Levodopa_LNAA_Vmax/Epithelium_CrossSection/(Levodopa_LNAA_Km + LevodopaEpithelium))*LevodopaEpithelium;


		// Should be g/(s*m**3). 
		forAll(SI_SurfaceArea,cellID) { 
			Epithelium_LevodopaAbsorption_From_SI[cellID] = SI_LevodopaAbsorption[cellID] /Epithelium_Volume[cellID];
		}		

		//                                                (kg/m3)/s   *  (kg/m3)         /(kg/m3) 			= kg/m3/s
 		volScalarField Epithelium_COMT_Metabolism = Levodopa_COMT_Vmax*LevodopaEpithelium/(Levodopa_COMT_Km + LevodopaEpithelium);

		fvScalarMatrix EpitheliumEqn (
			fvm::ddt(LevodopaEpithelium) 
					== 
			 Epithelium_LevodopaAbsorption_From_SI
			-Epithelium_COMT_Metabolism
			-Epithelium_Basal_Absorption
		); 	

		EpitheliumEqn.solve();
		Epithelium_TotalLevodopaAbsorption_From_SI 	+= Epithelium_LevodopaAbsorption_From_SI*Epithelium_Volume*runTime.deltaT();
		Epithelium_COMT_TotalMetabolism 		+= Epithelium_COMT_Metabolism*Epithelium_Volume*runTime.deltaT();
    		Epithelium_Flux_To_Body		 		 = Epithelium_Basal_Absorption*Epithelium_Volume;
		Epithelium_TotalFlux_To_Body			+= Epithelium_Flux_To_Body*runTime.deltaT();



		Info << endl << "   Epithelium  " << endl;
		Info <<         "---------------" << endl;
		Info << "Epithelium Levodopa " << sum(LevodopaEpithelium*Epithelium_Volume) << endl;
		Info << "Intesintal Absorption " << sum(Epithelium_TotalLevodopaAbsorption_From_SI) << endl;
		Info << "COMT metabolism  " << sum(Epithelium_COMT_TotalMetabolism) << endl;
		Info << "Basal Absorption " << sum(Epithelium_TotalFlux_To_Body) << endl;

	} // .. Epithelium 
  
   

	{ // BODY 

		Body_Flux_To_Body[0] = sum(Epithelium_Flux_To_Body).value()/bodyvolume.value();

 		volScalarField Body_Levodoa_Elimination = bodyelimination*LevodopaBody;

		fvScalarMatrix BodyEqn (
			fvm::ddt(LevodopaBody) 
					== 
			Body_Flux_To_Body
			-Body_Levodoa_Elimination

		); 	

		BodyEqn.solve();

		TotalElimination += Body_Levodoa_Elimination*bodyvolume*runTime.deltaT();

		Info << endl << "   Body  " << endl;
		Info <<         "---------------" << endl;
		Info << "Body Levodopa    " << sum(LevodopaBody*bodyvolume) << endl;
		Info << "Body Elimination " << sum(TotalElimination) << endl;
	
	}



	runTime.write();

	Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
			<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
			<< nl << endl;
            
    }

    Info<< "End\n" << endl;

    return 0;
}


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


	We know that the COMT has negligible effect on the bioavailability. 
	Therefore, we don't include it in this version of the code. 
	

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
//#include "singlePhaseTransportModel.H"
//#include "turbulentTransportModel.H"

#include "pimpleControl.H"
#include "interpolation.H"
#include "Random.H"
#include "meshSearch.H"

#include <sstream>
#include <fstream>

#include "Tablet.H"


//#include "stdlib.h"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    Info << " \n\nIntestinalLevodopa  :  0.1.0" << endl; 
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
    dimensionedScalar currentstomachemptying(stomachemptying);

    List<dimensionedScalar*> lagStartList; 
    List<dimensionedScalar*> lagEndList; 

    if (StomachPropertiesDict.found("lagStartList")) { 

		scalarList lagStartScalarList(StomachPropertiesDict.lookup("lagStartList"));
		scalarList lagEndWordList(StomachPropertiesDict.lookup("lagEndList"));

		forAll(lagStartScalarList, i) 
		{


			lagStartList.append(new dimensionedScalar("x",dimTime,lagStartScalarList[i])); 
			lagEndList.append  (new dimensionedScalar("x",dimTime,lagEndWordList[i])); 
		}

    } else { 

	    

	    dimensionedScalar lagstart(StomachPropertiesDict.lookup("lagStart"));    
	    dimensionedScalar lagend(StomachPropertiesDict.lookup("lagEnd"));    

	    lagStartList.append(new dimensionedScalar(lagstart));
	    lagEndList.append(  new dimensionedScalar(lagend)); 
    }


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
   
    volScalarField TotalLevodopaStomachEmptying(
					IOobject("Stomach_TotalLevodopaEmptying",
					runTime.timeName(),
					Stomach,
					IOobject::NO_READ,IOobject::AUTO_WRITE),
					Stomach,
					dimensionedScalar("zero",dimMass,scalar(0)) );

    volScalarField LevodopaStomachEmptying(IOobject("LevodopaStomachEmptying",runTime.timeName(),Stomach,IOobject::NO_READ,IOobject::AUTO_WRITE),
					   Stomach,dimensionedScalar("zero",dimMass/dimTime,scalar(0)) );

    volScalarField StomachTabletErosion(IOobject("TabletStomach_Erosion",runTime.timeName(),Stomach,IOobject::NO_READ,IOobject::NO_WRITE),
					   Stomach,dimensionedScalar("zero",dimMass/dimTime/dimVolume,scalar(0)) );

    volScalarField Total_StomachTabletErosion(IOobject("Stomach_TotalTabletErosion",runTime.timeName(),Stomach,IOobject::NO_READ,IOobject::AUTO_WRITE),
					   Stomach,dimensionedScalar("zero",dimMass,scalar(0)) );

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
	
    const bool absorptionWindow = SmallIntestinePropertiesDict.lookupOrDefault<bool>("AbsorptionWindow", false);

    if (absorptionWindow) { 
		Info << " ---------------> Using the absorption window " << endl;
    } else { 
		Info << " ---------------> NO absorption window " << endl;
    } 

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

	//Info << SmallIntestineVelocity << endl;	

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
   SmallIntestineVelocity[SmallIntestine.C().size()-1].component(0) = 0;
   

    surfaceScalarField SI_phi(
        IOobject
        (
            "SmallIntestinePhi",
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
		SI_phi.boundaryFieldRef()[IleoCecalpatchID][facei] =0; 
	} 


    label PyloruspatchID 		= SmallIntestine.boundaryMesh().findPatchID("Pylorus"); 
    const polyPatch& PylorusPatch 	= SmallIntestine.boundaryMesh()[PyloruspatchID]; 


	forAll(PylorusPatch, facei) 
	{ 
		SI_phi.boundaryFieldRef()[PyloruspatchID][facei] =0; 
	} 


    // Amplification (without the microvilli) per unit area of intestine. 
    volScalarField SI_SurfaceArea(IOobject("SI_SurfaceArea",runTime.timeName(),SmallIntestine,IOobject::NO_READ,IOobject::AUTO_WRITE),
					   SmallIntestine,dimensionedScalar("one",dimLength,scalar(1)) );
    

    volScalarField SI_dh(IOobject("SI_dH",runTime.timeName(),SmallIntestine,IOobject::NO_READ,IOobject::NO_WRITE),
					   SmallIntestine,dimensionedScalar("one",dimLength,scalar(1)) );


    volScalarField SI_FoldingAmplification(IOobject("SI_FoldingAmplification",runTime.timeName(),SmallIntestine,IOobject::NO_READ,IOobject::AUTO_WRITE),
					   SmallIntestine,dimensionedScalar("one",dimless,scalar(1)) );


    volScalarField SI_V(IOobject("SI_V",runTime.timeName(),SmallIntestine,IOobject::NO_READ,IOobject::AUTO_WRITE),pi*SI_Radius*SI_Radius*SI_dh);


    // folding = 3 from 10 to 50. 
    // decrease linearly from 3 to 1 in the range 50->282. 

    forAll(SI_SurfaceArea,cellID) { 
	// 		Krecking folding contribution. 
			scalar folding_amplification;
			const volScalarField& X = SmallIntestine.C().component(0); 

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

			bool setSI = true; 

			if (absorptionWindow && X[cellID] > 0.94) { 
				setSI = false; 
			}
		
			if (setSI) { 
				SI_SurfaceArea[cellID] 		= 2*pi*SI_Radius.value()*folding_amplification*villi_amplification;
			} else { 
				SI_SurfaceArea[cellID]          = 0;
			}

			SI_dh[cellID]  = SmallIntestine.V()[cellID];
    }

	//Info << SI_SurfaceArea;

    volScalarField SI_SurfaceArea_Microvilli(IOobject("SI_SurfaceArea",runTime.timeName(),SmallIntestine,IOobject::NO_READ,IOobject::NO_WRITE),
					   SmallIntestine,dimensionedScalar("one",dimless,scalar(25)) );   

    volScalarField SI_TabletErosion(IOobject("SI_TabletErosion",runTime.timeName(),SmallIntestine,IOobject::NO_READ,IOobject::NO_WRITE),
					   SmallIntestine,dimensionedScalar("zero",dimMass/dimTime/dimVolume,scalar(0)) );

    volScalarField SI_TotalTabletErosion(IOobject("SI_TotalTabletErosion",runTime.timeName(),SmallIntestine,IOobject::NO_READ,IOobject::AUTO_WRITE),
					   SmallIntestine,dimensionedScalar("zero",dimMass,scalar(0)) );

    dimensionedScalar SI_CrossSection(SmallIntestinePropertiesDict.lookup("SI_CrossSection"));

    // Levodopa properties. 
    dimensionedScalar Levodopa_LNAA_Vmax(SmallIntestineLevodopaPropertiesDict.lookup("LNAA_Vmax"));
    dimensionedScalar Levodopa_LNAA_Km(SmallIntestineLevodopaPropertiesDict.lookup("LNAA_Km"));

    volScalarField SI_Emptying(IOobject("SI_Emptying",runTime.timeName(),SmallIntestine,IOobject::NO_READ,IOobject::AUTO_WRITE),
					     SmallIntestine,dimensionedScalar("one",dimMass/dimTime/dimVolume,scalar(0)) );

    volScalarField SI_LevodopaAbsorption(IOobject("SI_LevodopaAbsorption",runTime.timeName(),SmallIntestine,IOobject::NO_READ,IOobject::AUTO_WRITE),
					     SmallIntestine,dimensionedScalar("one",dimMass/dimTime,scalar(0)) );


    volScalarField SI_LevodopaAbsorption_Coeff(IOobject("SI_LevodopaAbsorption_coeff",runTime.timeName(),SmallIntestine,IOobject::NO_READ,IOobject::NO_WRITE),
					     SmallIntestine,dimensionedScalar("one",dimless/dimTime,scalar(0)) );

  
    volScalarField SI_LevodopaTotalAbsorption(IOobject("SI_LevodopaTotalAbsorption",runTime.timeName(),SmallIntestine,IOobject::NO_READ,IOobject::AUTO_WRITE),
					     SmallIntestine,dimensionedScalar("one",dimMass,scalar(0)) );


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

    // Create the variables. 
    volScalarField LevodopaBodyC2
    (
        IOobject
        (
            "LevodopaBodyC2",
            runTime.timeName(),
            Body,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        Body
    );

    dimensionedScalar k10(BodyLevodopaPropertiesDict.lookup("k10"));
    dimensionedScalar k12(BodyLevodopaPropertiesDict.lookup("k12"));    
    dimensionedScalar k21(BodyLevodopaPropertiesDict.lookup("k21"));        
    dimensionedScalar V0(BodyLevodopaPropertiesDict.lookup("volume"));    

    volScalarField Body_Flux_To_Body(IOobject("Body_Flux_To_Body",runTime.timeName(),Body,IOobject::NO_READ,IOobject::NO_WRITE),
					     Body,dimensionedScalar("one",dimMass/(dimVolume*dimTime),scalar(0)) );

    volScalarField TotalBodyAbsorption(IOobject("TotalAbsorption",runTime.timeName(),Body,IOobject::NO_READ,IOobject::NO_WRITE),
					     Body,dimensionedScalar("one",dimMass,scalar(0)) );


    volScalarField TotalElimination(IOobject("TotalElimination",runTime.timeName(),Body,IOobject::NO_READ,IOobject::NO_WRITE),
					     Body,dimensionedScalar("one",dimMass,scalar(0)) );

    volScalarField C1_To_C2(IOobject("C1_To_C2",runTime.timeName(),Body,IOobject::NO_READ,IOobject::NO_WRITE),
					     Body,dimensionedScalar("one",dimMass,scalar(0)) );
    volScalarField C2_To_C1(IOobject("C2_To_C1",runTime.timeName(),Body,IOobject::NO_READ,IOobject::NO_WRITE),
					     Body,dimensionedScalar("one",dimMass,scalar(0)) );


    // --------------------------------------------------- Tablets ---------------------------------------
    IOdictionary tabletDict
    (
        IOobject
        (
            "Tablets",    // dictionary name
            runTime.constant(),     // dict is found in "constant"
            Stomach,                   // registry for the dict
            IOobject::MUST_READ,    // must exist, otherwise failure
            IOobject::NO_WRITE      // dict is only read by the solver
        )
    );
    wordList tabletNameList
	    (
		tabletDict.lookup("tablets")
	    );
    
    PtrList<Tablet> tabletList(tabletNameList.size());

    forAll(tabletNameList, tabletid)
    {
	tabletList.set(tabletid,new Tablet(runTime,
					   Stomach, 
					   SmallIntestine,
					   tabletDict.subDict(tabletNameList[tabletid]),
					   stomachvolume,
					   SI_V,
					   SmallIntestineVelocity)); 
    } 

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


    Info<< "\nStarting time loop\n" << endl;
    forAll(tabletList, tabletid)
    {
	tabletList[tabletid].setInitialConditions(LevodopaStomach,LevodopaSmallIntestine);
    }

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        //#include "readTimeControls.H"
        //#include "CourantNo.H"   
        //#include "setDeltaT.H"
	
	{ 	// Stomach
		currentstomachemptying = stomachemptying;
		forAll(lagStartList,i) { 
			if ((runTime.time() > *(lagStartList[i])) && (runTime.time() < *(lagEndList[i]))) { 
				currentstomachemptying.value() = 0;
			}
		}

		StomachTabletErosion.ref()[0] = 0;
		forAll(tabletList, tabletid)
		{
			StomachTabletErosion += tabletList[tabletid].getStomachSource();
		}


		Total_StomachTabletErosion += StomachTabletErosion*runTime.deltaT()*stomachvolume;

		fvScalarMatrix LevodopaStomachEqn (
			fvm::ddt(LevodopaStomach) == fvm::Sp(-currentstomachemptying,LevodopaStomach) + StomachTabletErosion
		); 	
	
		LevodopaStomachEqn.solve();
		LevodopaStomachEmptying  = fvc::Sp(currentstomachemptying,LevodopaStomach)*stomachvolume; 
		TotalLevodopaStomachEmptying += LevodopaStomachEmptying*runTime.deltaT();
		
		Info << endl << "Stomach" << endl;
		Info << "--------" << endl; 
		Info << "In Stomach " << (LevodopaStomach*stomachvolume)->ref()[0] << endl;
		Info << "Emptied " << TotalLevodopaStomachEmptying()[0] << endl;
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

		forAll(SI_TabletErosion,cellid) { 
			SI_TabletErosion.ref()[cellid] = 0;
		}

		forAll(tabletList, tabletid)
		{
			SI_TabletErosion += tabletList[tabletid].getIntestineSource();
		}

		SI_TotalTabletErosion+= SI_TabletErosion*SI_V*runTime.deltaT();

		fvScalarMatrix LevodopaIntestineEqn (
			fvm::ddt(LevodopaSmallIntestine)  
			+ fvm::div(SI_phi,LevodopaSmallIntestine)  
			- fvm::laplacian(SmallIntestineDispersion,LevodopaSmallIntestine) 
					== 
			SI_Emptying
			+ SI_TabletErosion 
			- TimeStepAbs
		      
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

			SmallIntestineVelocity[SmallIntestine.C().size()-1] = SI_U.value();

		}

		

	} //.. small intestine. 

	{ // BODY 

		Body_Flux_To_Body[0] = sum(SI_LevodopaAbsorption).value()/V0.value();
 		volScalarField Body_Levodoa_C2ToC1 = k21*LevodopaBodyC2;
		fvScalarMatrix BodyEqn (
			fvm::ddt(LevodopaBody) 
					== 
			Body_Flux_To_Body
			-fvm::Sp(k12,LevodopaBody)
			- fvc::Sp(k12,LevodopaBody)
			+ Body_Levodoa_C2ToC1
			-fvm::Sp(k10,LevodopaBody) //Body_Levodoa_Elimination


		); 	
 		volScalarField Body_Levodoa_Elimination = fvc::Sp(k10,LevodopaBody);
 		volScalarField Body_Levodoa_C1ToC2 = fvc::Sp(k12,LevodopaBody);

		fvScalarMatrix BodyEqncC2 (
			fvm::ddt(LevodopaBodyC2) 
					== 
			  fvc::Sp(k12,LevodopaBody)
			- Body_Levodoa_C2ToC1
		); 	


		BodyEqn.solve();
		BodyEqncC2.solve();


		TotalBodyAbsorption += Body_Flux_To_Body*V0*runTime.deltaT();
		TotalElimination += Body_Levodoa_Elimination*V0*runTime.deltaT();

		Info << endl << "   Body  " << endl;
		Info <<         "---------------" << endl;
		Info << "Body Levodopa    " << sum(LevodopaBody*V0) << endl;
		Info << "Body Absorption  " << sum(TotalBodyAbsorption) << endl;
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


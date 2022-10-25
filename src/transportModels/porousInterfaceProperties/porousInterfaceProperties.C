/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

// #include "dimensionSets.H"
#include "porousInterfaceProperties.H"
#include "alphaContactAngleFvPatchScalarField.H"
#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcAverage.H"
#include "fvCFD.H"


// * * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * //

const Foam::scalar Foam::interfaceProperties::convertToRad =
    Foam::constant::mathematical::pi/180.0;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Correction for the boundary condition on the unit normal nHat on
// walls to produce the correct contact angle.

// The dynamic contact angle is calculated from the component of the
// velocity on the direction of the interface, parallel to the wall.
void Foam::interfaceProperties::correctContactAngle
(
    surfaceVectorField::Boundary& nHatb,
    const surfaceVectorField::Boundary& gradAlphaf
) const
{
    const fvMesh& mesh = alpha1_.mesh();
    const volScalarField::Boundary& abf = alpha1_.boundaryField();

    const fvBoundaryMesh& boundary = mesh.boundary();

    forAll(boundary, patchi)
    {
        if (isA<alphaContactAngleFvPatchScalarField>(abf[patchi]))
        {
            alphaContactAngleFvPatchScalarField& acap =
                const_cast<alphaContactAngleFvPatchScalarField&>
                (
                    refCast<const alphaContactAngleFvPatchScalarField>
                    (
                        abf[patchi]
                    )
                );

            fvsPatchVectorField& nHatp = nHatb[patchi];
            const scalarField theta
            (
                convertToRad*acap.theta(U_.boundaryField()[patchi], nHatp)
            );

            const vectorField nf
            (
                boundary[patchi].nf()
            );

            // Reset nHatp to correspond to the contact angle

            const scalarField a12(nHatp & nf); 

            const scalarField b1(cos(theta));

            scalarField b2(nHatp.size());
            forAll(b2, facei)
            {
                b2[facei] = cos(acos(a12[facei]) - theta[facei]);
            }

            const scalarField det(1.0 - a12*a12);

            scalarField a((b1 - a12*b2)/det);
            scalarField b((b2 - a12*b1)/det);

            nHatp = a*nf + b*nHatp;
            nHatp /= (mag(nHatp) + deltaN_.value());

            acap.gradient() = (nf & nHatp)*mag(gradAlphaf[patchi]);
            acap.evaluate();
        }
    }
}


// The dynamic contact angle is calculated from the component of the
// velocity on the direction of the interface, parallel to the porosity wall
void Foam::interfaceProperties::correctContactAngleDB
(
    surfaceVectorField& nHatb, 
    surfaceVectorField& gradAlphaf
) const
{
    // Reading Values From Main Simulation
    volScalarField Solid_ = alpha1_.mesh().lookupObject<volScalarField>("Solid"); 
    scalar SmootherCount = 0;

    // Defining the Solid Surface Interface
    const surfaceScalarField SolidSurface(fvc::interpolate(Solid_));
    //forAll(SolidSurface, facei) //setting all face B values to 1 at solid interface and zero otherwise
    //        {
    //            if (SolidSurface[facei]>0) //>0 ==1
    //            {SolidSurface[facei]=1;}
    //        }
  
    //Calculating surface gradient at the face and smoothing
    const volVectorField gradSolid(fvc::grad(Solid_, "nHatSolid"));
    surfaceVectorField gradSolidf(fvc::interpolate(gradSolid));
    surfaceVectorField nHatSolidfv(gradSolidf/(mag(gradSolidf) + deltaN_));

    //the previous smoothing scheme of the porous boundary  
    if (porousBoundaryNormSmoothingMethod_ == 0)
    {
        while (SmootherCount <= porousBoundaryNormSmoothingCycles_)
        {
        if (SmootherCount < porousBoundaryNormSmoothingCycles_)
            {gradSolidf = SolidSurface*fvc::interpolate(fvc::average(gradSolidf));}
        else 
            {gradSolidf = fvc::interpolate(fvc::average(gradSolidf));}
        SmootherCount +=1;  
        }
        nHatSolidfv = gradSolidf/(mag(gradSolidf) + deltaN_);
    }

    // improved smoothing scheme of the porous boundary  
    if (porousBoundaryNormSmoothingMethod_ == 1)
    {
        volVectorField nHatSolidv(gradSolid/(mag(gradSolid) + deltaN_));
        volScalarField nHatSolidvMag0(mag(nHatSolidv));
        volScalarField nHatSolidvMag(max(nHatSolidvMag0,deltaN_.value()));

        while (SmootherCount <= porousBoundaryNormSmoothingCycles_)
        {
            nHatSolidv = nHatSolidvMag0*fvc::average(fvc::interpolate(nHatSolidv*nHatSolidvMag))/fvc::average(fvc::interpolate(nHatSolidvMag));
            nHatSolidv = nHatSolidv/(mag(nHatSolidv)+deltaN_.value());
            SmootherCount +=1;  
        }
        nHatSolidfv = fvc::interpolate(nHatSolidv);
    }

    //Converting Input Contact Angle to Radians
    scalar theta0Rad_ = convertToRad*theta0_;
 
	    //dot product between both normals = cos(thetaI)
            const surfaceScalarField a12(nHatb & nHatSolidfv);
	
	    //cos(theta0)
            const scalar b1(cos(theta0Rad_));

	    //cos(thetaDiff)
            surfaceScalarField b2(cos(acos(a12) - theta0Rad_));

	    //1-cos^2(thetaI)
            const surfaceScalarField det(1.0 - a12*a12);

	    //cos(theta0) - cos(thetaI)*cos(thetaDiff)/ ( 1-cos^2(thetaI) )
            surfaceScalarField a((b1 - a12*b2)/det);

	    //cos(thetaDiff) - cos(thetaI)*cos(theta0)/ ( 1-cos^2(thetaI) )
            surfaceScalarField b((b2 - a12*b1)/det);

	    // correction of the interface and re-normalizing
	    nHatb = a*nHatSolidfv + b*nHatb;
            nHatb /= (mag(nHatb) + deltaN_.value());
}

/////Calculating the curvature
void Foam::interfaceProperties::calculateK()
{
    
    const fvMesh& mesh = alpha1_.mesh();
    const surfaceVectorField& Sf = mesh.Sf();
    volScalarField Solid_ = alpha1_.mesh().lookupObject<volScalarField>("Solid"); 
    surfaceScalarField Fluidf_ = alpha1_.mesh().lookupObject<surfaceScalarField>("Fluidf");  

    volScalarField talpha2(1.0-alpha1_);
    dimensionedScalar small("small", dimless, 1e-10);
    //alpha in fluid region is extrapolated to the near porous interface cell
    //if alpha in solid =0 (gas)
    if (alphaInPorousRegion_ == 0)  
    {
        alpha1_ave = (fvc::average(fvc::interpolate(max(talpha2,small),"Fluid_c2f")*Fluidf_)/fvc::average(Fluidf_));
        alpha1_ave = (1.0-Solid_)*talpha2 + Solid_*alpha1_ave;
        alpha1_ave = 1.0 - alpha1_ave;
    }
   //if alpha in solid =1 (water)
    if (alphaInPorousRegion_ == 1)  
    {
        alpha1_ave = (fvc::average(fvc::interpolate(max(alpha1_,small),"Fluid_c2f")*Fluidf_)/fvc::average(Fluidf_));
        alpha1_ave = (1.0-Solid_)*alpha1_ + Solid_*alpha1_ave;
    }

    // clip alpha field to exclude the localized unphysical variations in alpha close to 0 or 1 in the interface vicinity
    volScalarField alpha1s ( min(max((1.0+2.0*alphaThresholdK_)*(alpha1_ave-0.5)+0.5,0.0),1.0) );  
    // smooth the alpha field in the fluid region
    for (int i=0;i<nSK_;i++)
    {
        alpha1s = Solid_*alpha1_ave + (1.0-Solid_)*(cSK_ * fvc::average(fvc::interpolate(alpha1s*(1.0-Solid_))) \
            /fvc::average(fvc::interpolate(max(1.0-Solid_,small))) + (1.0 - cSK_) * alpha1s);
    }
     //gradient interpolation
    const volVectorField gradAlpha(fvc::grad(alpha1s, "nHat")); 
    surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha)); 
    surfaceVectorField nHatfv_g (gradAlphaf/(mag(gradAlphaf) + deltaN_)); 
    
    //normal interpolation
    nI_ = gradAlpha/(Foam::mag(gradAlpha) + deltaN_);
    surfaceVectorField nHatfv_n(fvc::interpolate(nI_));   

    //hybrid interpolation
    surfaceVectorField nHatfv_ave (nHatfv_n*hybridNlgCoeff_+nHatfv_g*(1.0-hybridNlgCoeff_));

    if (activatePorousContactAngle_==1)
    {
        correctContactAngleDB(nHatfv_ave,gradAlphaf);
    }
    correctContactAngle(nHatfv_ave.boundaryFieldRef(), gradAlphaf.boundaryField());

    // Face unit interface normal flux
    nHatf_ = nHatfv_ave & Sf; 
    
    // Complex expression for curvature.
    // Correction is formally zero but numerically non-zero.
    //  KCorr=0 simple, KCorr=1 Complex (default simple)
    K_ = -fvc::div(nHatf_) + (nI_ & fvc::grad(nHatfv_ave) & nI_)*KCorr_;

    K_.correctBoundaryConditions();

    //the curvation correction based on the SSF scheme
    if (nSC_>0)
    {
        volScalarField alpha1c_ = Foam::min(1.0, Foam::max(alpha1_ave, 0.0))();  
      
        volScalarField w = Foam::sqrt(alpha1c_*(1.0 - alpha1c_) + 1e-6)();
        volScalarField factor = (2.0 * Foam::sqrt(alpha1c_*(1.0 - alpha1c_)))();
        volScalarField KStar = (fvc::average(fvc::interpolate(K_*w*(1.0-Solid_)))/fvc::average(fvc::interpolate(w*max(1.0-Solid_,small))))();
        volScalarField Ks = (factor * K_ + (1.0 - factor) * KStar)();

        for (int i=1;i<nSC_;i++)
        {
            KStar = fvc::average(fvc::interpolate(Ks*w*(1.0-Solid_)))/fvc::average(fvc::interpolate(w*max(1.0-Solid_,small)));
            Ks = factor * K_ + (1.0 - factor) * KStar;
        }

        Kf_ = fvc::interpolate(w*Ks*(1.0-Solid_))/fvc::interpolate(w*max(1.0-Solid_,small));
    }

}


void Foam::interfaceProperties::calculateFc() 
{
    const fvMesh& mesh = alpha1_.mesh();
    surfaceScalarField Solidf = alpha1_.mesh().lookupObject<surfaceScalarField>("Solidf");  

    // clip alpha field to exclude the localized unphysical variations in alpha close to 0 or 1 in the interface vicinity
    alpha1_pc = min(max((1.0+2.0*alphaThresholdFc_)*(alpha1_ave-0.5)+0.5,0.0),1.0); 
    //sharp the alpha field based on the SSF scheme
    alpha1_pc = 1.0/(1.0-cPc_)*(min( max(alpha1_pc,cPc_/2.0), (1.0-cPc_/2.0) ) - cPc_/2.0);

    alpha1_pc.correctBoundaryConditions();

    deltasf_ = fvc::snGrad(alpha1_pc);
    stf_ = sigmaKSSF()*deltasf_;
    fc_ = fvc::reconstruct(stf_*(1-Solidf)*mesh.magSf());
}
 
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceProperties::interfaceProperties
(
    const volScalarField& alpha1,
    const volVectorField& U,
    const IOdictionary& dict
)
:
    transportPropertiesDict_(dict),
    cAlpha_
    (
        readScalar
        (
            alpha1.mesh().solverDict(alpha1.name()).lookup("cAlpha")
        )
    ),

    theta0_ 
    (
           transportPropertiesDict_.lookupOrDefault("theta0",90)
    ),

    cSK_
    (
        readScalar
        (
            alpha1.mesh().solverDict(alpha1.name()).lookup("cSK")
        )
    ),

    nSK_
    (
        readScalar
        (
            alpha1.mesh().solverDict(alpha1.name()).lookup("nSK")
        )
    ),

    cPc_
    (
       readScalar
       (
           alpha1.mesh().solverDict(alpha1.name()).lookup("cPc")
       )
    ),

    nSC_
    (
        readScalar
        (
            alpha1.mesh().solverDict(alpha1.name()).lookup("nSC")
        )
    ),

    KCorr_
    (
        readScalar
        (
            alpha1.mesh().solverDict(alpha1.name()).lookup("KCorr")
        )
    ),

    porousBoundaryNormSmoothingCycles_
    (
        alpha1.mesh().solverDict(alpha1.name()).lookupOrDefault("porousBoundaryNormSmoothingCycles",1)
    ),

    porousBoundaryNormSmoothingMethod_
    (
        alpha1.mesh().solverDict(alpha1.name()).lookupOrDefault("porousBoundaryNormSmoothingMethod",1)
    ),

    alphaInPorousRegion_
    (
        alpha1.mesh().solverDict(alpha1.name()).lookupOrDefault("alphaInPorousRegion",1)
    ),

    hybridNlgCoeff_
    (
        alpha1.mesh().solverDict(alpha1.name()).lookupOrDefault("hybridNlgCoeff",1.0)
    ),

    alphaThresholdK_
   (
       alpha1.mesh().solverDict(alpha1.name()).lookupOrDefault("alphaThresholdK",0.0)
    ),

    alphaThresholdFc_
    (
        alpha1.mesh().solverDict(alpha1.name()).lookupOrDefault("alphaThresholdFc",0.001)
    ),

    activatePorousContactAngle_
    (
        transportPropertiesDict_.lookupOrDefault("activatePorousContactAngle",0)
    ),

    sigmaPtr_(surfaceTensionModel::New(dict, alpha1.mesh())),

    deltaN_
    (
        "deltaN",
        1e-8/pow(average(alpha1.mesh().V()), 1.0/3.0)
    ), 

    alpha1_(alpha1),

    U_(U),

    nHatf_
    (
        IOobject
        (
            "nHatf",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar("nHatf", dimArea, 0.0)
    ),

    //gas-fluid interface normal vector at cell center
    nI_
    (
        IOobject
        (
            "nI",
            alpha1_.time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1_.mesh(),
        dimensionedVector(dimless, vector(0.0,0.0,0.0)) //(0.0 0.0 0.0)
    ),

    K_
    (
        IOobject
        (
            "interfaceProperties_K",
            alpha1_.time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1_.mesh(),
        dimensionedScalar("K", dimless/dimLength, 0.0)
    ),

    Kf_
    (
        IOobject
        (
            "interfaceProperties:Kf",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar("Kf", dimless/dimLength, 0.0)
    ),

    alpha1_ave
    (
        IOobject
        (
            "alpha.wetting.forKlg",
            alpha1_.time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1_
    ),

    alpha1_pc
    (
        IOobject
        (
            "alpha1_pc",
            alpha1_.time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        alpha1_
    ),

    fc_
    (
        IOobject
        (
            "fc",
            alpha1_.time().timeName(),
            alpha1_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alpha1_.mesh(),
        dimensionedVector(Foam::dimPressure/dimLength, vector(0.0,0.0,0.0)) //(0.0 0.0 0.0)
    ),
    
    deltasf_
    (
        IOobject
        (
            "deltasf",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar( dimless/dimLength, 0.0)
    ),

    stf_
    (
        IOobject
        (
            "stf",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar( dimMass/(dimLength*dimLength*dimTime*dimTime), 0.0)
    )
{
    calculateK();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
Foam::tmp<Foam::surfaceScalarField>
Foam::interfaceProperties::sigmaKSSF() const
{
    if (nSC_>0)
        {return fvc::interpolate( sigmaPtr_->sigma() )*Kf_;}    // smooth K_
    else
        {return fvc::interpolate( (sigmaPtr_->sigma())*K_);}  //don't smooth K_ 
}

Foam::tmp<Foam::surfaceScalarField>
Foam::interfaceProperties::surfaceTensionForce() const
{
   return stf_; 
}


Foam::tmp<Foam::surfaceScalarField>
Foam::interfaceProperties::deltasf() const
{
   return deltasf_; 
}


Foam::tmp<Foam::volScalarField>
Foam::interfaceProperties::Klg() const
{
    return K_;
}

Foam::tmp<Foam::volVectorField>
Foam::interfaceProperties::nlg() const
{
    return nI_;
}

Foam::tmp<Foam::volVectorField>
Foam::interfaceProperties::fc() const
{
    return fc_;
}


Foam::tmp<Foam::volScalarField>
Foam::interfaceProperties::alphapc() const
{
    return alpha1_pc;
}

Foam::tmp<Foam::volScalarField>
Foam::interfaceProperties::nearInterface() const
{
    return pos(alpha1_ - 0.01)*pos(0.99 - alpha1_);
}


void Foam::interfaceProperties::correct()
{
    calculateK();
    calculateFc();
}


bool Foam::interfaceProperties::read()
{
    alpha1_.mesh().solverDict(alpha1_.name()).lookup("cAlpha") >> cAlpha_;
    sigmaPtr_->readDict(transportPropertiesDict_);

    return true;
}


// ************************************************************************* //

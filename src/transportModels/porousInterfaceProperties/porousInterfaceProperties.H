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

Class
    Foam::interfaceProperties

Description
    Contains the interface properties.

    Properties to aid interFoam:
    -# Correct the alpha boundary condition for dynamic contact angle.
    -# Calculate interface curvature.

SourceFiles
    interfaceProperties.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceProperties_H
#define interfaceProperties_H

#include "IOdictionary.H"
#include "scalar.H"
#include "surfaceTensionModel.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class interfaceProperties Declaration
\*---------------------------------------------------------------------------*/

class interfaceProperties
{
    // Private data

        //- Keep a reference to the transportProperties dictionary
        const dictionary& transportPropertiesDict_;

        //- Compression coefficient
        scalar cAlpha_;

	//- Contact Angle
	scalar theta0_;

        //- smoothing coefficient in the SSF scheme
        scalar cSK_;
        //-number of smoothing cycle in the SSF scheme
        label nSK_;
        //- Sharp force coefficient in the SSF scheme (put 0.98-0.99 for static problems, 0.4-0.5 for dynamic)
        scalar cPc_;
        //-number of smoothing cycle for curvature K_ in the SSF scheme
        label nSC_;
        //-curvature K_ calculation model: KCorr=0 simple, KCorr=1 Complex (default simple)
        scalar KCorr_;
        //-number of smoothing cycle for the porous boundary normal
        scalar porousBoundaryNormSmoothingCycles_;
        //- model index of the smoothing scheme for the porous boundary normal
        int porousBoundaryNormSmoothingMethod_;
        //- indicate the initial phase in the pore space of the porous region
        //0: gas, 1: fluid
        int alphaInPorousRegion_;
        //- coefficient of the hybrid scheme for the gas-fluid interface normal
        scalar hybridNlgCoeff_;
        //- coefficient of the alpha clipping scheme for the curvation computation
        scalar alphaThresholdK_;
        //- coefficient of the alpha clipping scheme for the surface tension computation
        scalar alphaThresholdFc_;

        //- Activate Contact Angle
        scalar activatePorousContactAngle_;

        //- Surface tension
        autoPtr<surfaceTensionModel> sigmaPtr_;

        //- Stabilisation for normalisation of the interface normal
        const dimensionedScalar deltaN_;
	    const volScalarField& alpha1_;
        const volVectorField& U_;
        surfaceScalarField nHatf_;
        volVectorField nI_;
        volScalarField K_;
        surfaceScalarField Kf_;
        volScalarField alpha1_ave;
        volScalarField  alpha1_pc;
        volVectorField fc_;

        surfaceScalarField deltasf_;
        surfaceScalarField stf_;

    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        interfaceProperties(const interfaceProperties&);
        void operator=(const interfaceProperties&);

        //- Correction for the boundary condition on the unit normal nHat on
        //  walls to produce the correct contact dynamic angle
        //  calculated from the component of U parallel to the wall
        void correctContactAngle
        (
            surfaceVectorField::Boundary& nHat,
            const surfaceVectorField::Boundary& gradAlphaf
        ) const;

	//- Correction for the boundary condition on the unit normal SolidnHat on
        //  porous walls to produce the correct contact dynamic angle
        //  calculated from the component of U parallel to the wall
        void correctContactAngleDB
        (
            surfaceVectorField& nHat,
	    surfaceVectorField& gradAlphaf
        ) const;

        //- Re-calculate the interface curvature
        void calculateK();
        void calculateFc();


public:

    //- Conversion factor for degrees into radians
    static const scalar convertToRad;


    // Constructors

        //- Construct from volume fraction field gamma and IOdictionary
        interfaceProperties
        (
	    const volScalarField& alpha1,
	    const volVectorField& U,
            const IOdictionary&
        );


    // Member Functions

        scalar cAlpha() const
        {
            return cAlpha_;
        }

        const dimensionedScalar& deltaN() const
        {
            return deltaN_;
        }

        const surfaceScalarField& nHatf() const
        {
            return nHatf_;
        }

        tmp<surfaceScalarField> sigmaKSSF() const;

        tmp<surfaceScalarField> surfaceTensionForce() const;

        tmp<surfaceScalarField> deltasf() const;

        tmp<volScalarField> Klg() const;

        tmp<volVectorField> nlg() const;

        tmp<volVectorField> fc() const;

        tmp<volScalarField> alphapc() const;

        tmp<volScalarField> nearInterface() const;

        void correct();

        //- Read transportProperties dictionary
        bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

// Defenition of Solid Indicator
volScalarField Solid
(
    IOobject
    (
        "Solid",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar ("Solid", dimensionSet(0,0,0,0,0,0,0), 0) 
);

//Labeling the places that have porous media
forAll(eps,celli) 
{
	if(eps[celli]>=0.99) 
	{	
	Solid[celli]=0; 
	}
	else
	{
	Solid[celli]=1; 
	}
}

forAll(mesh.boundary(), patchi) 
{
	forAll(Solid.boundaryField()[patchi],facei)
	{
		if(eps.boundaryField()[patchi][facei]>=0.99) 
		{Solid.boundaryFieldRef()[patchi][facei]=0;}
		else
		{Solid.boundaryFieldRef()[patchi][facei]=1;}
	}
}

surfaceScalarField Solidf
(
    IOobject
    (
        "Solidf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ, 
        IOobject::AUTO_WRITE
    ),
    fvc::interpolate(Solid)
);

//choose the interface surface tension model
// 0: vof model for the fluid region, 1: porous model for the porous region
int mode_interface_fc(transportProperties.lookupOrDefault("Mifc", 0));

//Eliminating intermediate values of Solid Indicator Function
forAll(Solidf,facei)
{
	if(mode_interface_fc==0)
    {
        if(Solidf[facei]<1) // if <1 ==0.  
        {	
            //when calculting the surface tension
            //porous interface is treated as the fluid phase 
            Solidf[facei]=0; 
        }
    }
    else
    {
        if(Solidf[facei]>0) // if <1 ==0.  
        {	
            //when calculting the surface tension
            //porous interface is treated as the porous phase 
            Solidf[facei]=1; 
        } 
    }
}

// Defenition of Solid Indicator
surfaceScalarField Fluidf
(
    IOobject
    (
        "Fluidf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ, //READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    max(1.0-fvc::interpolate(Solid),SMALL)
);
Fluidf = max(1.0-Solidf,SMALL);

// Defenition of Porous Indicator
volScalarField Porous
(
    IOobject
    (
        "Porous",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar ("Porous", dimensionSet(0,0,0,0,0,0,0), 0) 
);

forAll(eps,celli) 
{
	if(eps[celli]>0.01 && eps[celli]<0.99) 
	{	
	Porous[celli]=1; 
	}
	else
	{
	Porous[celli]=0; 
	}
}

    Info<< "Reading g\n" << endl;
    uniformDimensionedVectorField g
    (
        IOobject
        (
            "g",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    Info<< "Reading transportProperties\n" << endl;
    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar Sc(transportProperties.lookup("Sc"));

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Prepared for compressibleTwoPhaseMixtureThermo
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimDensity
    );

    Info<< "Constructing compressibleTwoPhaseMixtureThermo\n" << endl;
    compressibleTwoPhaseMixtureThermo thermo(mesh);

    volScalarField& YLiq = thermo.YLiq();
    volScalarField& YGas = thermo.YGas();
    volScalarField& YbarLiq = thermo.YbarLiq();
    volScalarField& YbarGas = thermo.YbarGas();

    Info<< "Reading thermophysical properties\n" << endl;
    volScalarField& p = thermo.p();
    volScalarField& T = thermo.T();
    volScalarField& rhoLiq = thermo.thermoLiq().rho();
    volScalarField& rhoGas = thermo.thermoGas().rho();
    const volScalarField& psiLiq = thermo.thermoLiq().psi();
    const volScalarField& psiGas = thermo.thermoGas().psi();
    volScalarField& psi = const_cast<volScalarField&>(thermo.psi());

    Info<< "Calculating mass flux field phi\n" << endl;
    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(rho * U) & mesh.Sf()
    );

    Info<< "Creating turbulence model\n" << endl;
    autoPtr<compressible::RASModel> turbulence
    (
        compressible::RASModel::New(rho, U, phi, thermo)
    );

    Info<< "Creating field dpdt\n" << endl;
    volScalarField dpdt
    (
        IOobject
        (
            "dpdt",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("dpdt", p.dimensions()/dimTime, 0)
    );
    
    volScalarField EkChange
    (
	IOobject
	(
	    "EkChange",
	    runTime.timeName(),
	    mesh
	),
	mesh,
	dimensionedScalar("dpdt", p.dimensions()/dimTime, 0)
    );
    
    Info<< "Creating field gradp\n" << endl;
    volVectorField gradp
    (
        IOobject
        (
            "gradp",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("gradp", p.dimensions()/dimLength, 0) * pTraits<vector>::zero
    );
    
    volVectorField limGradP
    (
	"limGradP",
	gradp
    );
    
    Info<< "Creating field kinetic energy K\n" << endl;
    volScalarField K("K", 0.5*magSqr(U));
    
    DynamicList<label> badQualityCells(0);
    
    surfaceScalarField pos
    (
	IOobject
	(
	    "pos",
	    runTime.timeName(),
	    mesh
	),
	mesh,
	dimensionedScalar("pos", dimless, 1.0)
    );
    
    surfaceScalarField neg
    (
	IOobject
	(
	    "neg",
	    runTime.timeName(),
	    mesh
	),
	mesh,
	dimensionedScalar("neg", dimless, -1.0)
    );
    
    volScalarField& rhoLiq0 = thermo.thermoLiq().rhoZero();
    
    Info << "max/min rhoLiq0" << max(rhoLiq0).value() << "/" << min(rhoLiq0).value() << endl;
    
    volScalarField alphaSqrRhoLiq0
    (
	"alphaSqrRhoLiq0",
	YLiq*YbarLiq*rhoLiq0
    );
    
    surfaceScalarField alphaSqrRhoLiq0_pos
    (
	"alphaSqrRhoLiq0_pos",
	fvc::interpolate(alphaSqrRhoLiq0, pos, "reconstruct(p)")
    );
    
    surfaceScalarField alphaSqrRhoLiq0_neg
    (
	"alphaSqrRhoLiq0_neg",
	fvc::interpolate(alphaSqrRhoLiq0, neg, "reconstruct(p)")
    );

    Info<< "Declare and initialise the cumulative continuity error\n" << endl;
    scalar cumulativeContErr = 0;
    
    surfaceScalarField liqDiffusionFlux
    (
	"liqDiffusionFlux",
	fvc::interpolate(turbulence->muEff()) * ( (fvc::snGrad(YLiq)) * mesh.magSf()) * 0.0
    );
    
    surfaceScalarField gasDiffusionFlux
    (
	"gasDiffusionFlux",
	liqDiffusionFlux
    );

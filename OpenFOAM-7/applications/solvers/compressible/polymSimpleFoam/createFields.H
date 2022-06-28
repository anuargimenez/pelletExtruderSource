Info<< "Reading thermophysical properties\n" << endl;

autoPtr<fluidThermo> pThermo //create the object "thermo"
(
    fluidThermo::New(mesh)
);
fluidThermo& thermo = pThermo();

thermo.validate(args.executable(), "h", "e");

volScalarField& p = thermo.p();

volScalarField Cp = thermo.Cp(); // To be able to plot it in post processing
volScalarField kappa = thermo.kappa(); // To be able to plot it in post processing

volScalarField rho //create scalar field rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

Info<< "Reading field U\n" << endl; //create vector field U
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

#include "compressibleCreatePhi.H" //creates and initializes the relative face-flux field phi

pressureControl pressureControl(p, rho, simple.dict());

mesh.setFluxRequired(p.name());

Info<< "Creating turbulence model\n" << endl; 
autoPtr<compressible::turbulenceModel> turbulence //turbulence model is created  --> Foam::compressible::turbulenceModel  steckt in turbulentFluidThermoModel.H
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

dimensionedScalar initialMass = fvc::domainIntegrate(rho);

#include "createMRF.H"
#include "createFvOptions.H"

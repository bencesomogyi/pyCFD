    Info << "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    const volScalarField::GeometricBoundaryField& pbf = p.boundaryField();
    wordList pbt = pbf.types();

    volScalarField pp
    (
        IOobject
        (

            "pp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", p.dimensions(), 0.0),
        pbt
    );


    // Set pp boundary values
    forAll(pp.boundaryField(), patchi)
    {
   	if (isType<fixedValueFvPatchScalarField>(pp.boundaryField()[patchi]))
        {
           fixedValueFvPatchScalarField& ppbound = refCast<fixedValueFvPatchScalarField>(pp.boundaryField()[patchi]);

           ppbound == scalarField(ppbound.size(),0.0);
        }

    }


    Info << "Reading field U\n" << endl;
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


#   include "createPhi.H"

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("SIMPLE"), pRefCell, pRefValue);


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

    dimensionedScalar nu
    (
        transportProperties.lookup("nu")
    );



//SIMPLE correction 

volScalarField DU = 1.0/UEqn.A();


volVectorField gradP  = fvc::grad(p);
surfaceScalarField DUf("DUf",linearInterpolate(DU));
surfaceVectorField gradp_avg_f = linearInterpolate(gradP);
surfaceScalarField gradp_f = fvc::snGrad(p);//snGrad is a scalar!

surfaceVectorField U_avg_f = linearInterpolate(U);

// Rhie-Chow interplation
phi = (U_avg_f & mesh.Sf()) - DUf*( (gradp_f*mesh.magSf()) - (gradp_avg_f&mesh.Sf()));

pp = scalar(0.0)*pp; // re-init p correction
pp.correctBoundaryConditions(); // re-apply bcs

fvScalarMatrix ppEqn
(
    - fvm::laplacian(DUf, pp, "laplacian(pDiff,pp)") // tell the solver which scheme to use
    + fvc::div(phi)
);

ppEqn.setReference(pRefCell, pRefValue);

ppEqn.solve();

phi += ppEqn.flux();

scalar URF = mesh.equationRelaxationFactor("pp");
p += URF*pp;
p.correctBoundaryConditions();

U -= fvc::grad(pp)*DU;
U.correctBoundaryConditions();



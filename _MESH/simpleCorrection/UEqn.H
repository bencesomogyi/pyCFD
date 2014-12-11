    // Solve the Momentum equation

    fvVectorMatrix UEqn
    (
            fvm::div(phi, U)
          - fvm::laplacian(nu, U)
    );

    UEqn.relax();

    solve
    (
        UEqn == -fvc::grad(p)
    );




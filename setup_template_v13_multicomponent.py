import os

def create_base_template():
    base_dir = os.path.join("templates", "base_case_multicomp")
    dirs = [
        os.path.join(base_dir, "0"),
        os.path.join(base_dir, "constant"),
        os.path.join(base_dir, "constant", "triSurface"),
        os.path.join(base_dir, "system")
    ]

    for d in dirs:
        os.makedirs(d, exist_ok=True)

    print(f"Created directories in {base_dir}")

    # --- 1. PHYSICS ---
    mom_transport = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v13                                   |
|   \\\\  /    A nd           | Website:  www.openfoam.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile { version 2.0; format ascii; class dictionary; object momentumTransport; }

simulationType RAS;

RAS { model kOmegaSST; turbulence on; printCoeffs on; }
"""

    thermophysical_transport = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v13                                   |
|   \\\\  /    A nd           | Website:  www.openfoam.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile { version 2.0; format ascii; class dictionary; object thermophysicalTransport; }

RAS
{
    model   unityLewisEddyDiffusivity;
    Prt     0.85;
    Sct     0.85;
}
"""

    combustion_properties = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v13                                   |
|   \\\\  /    A nd           | Website:  www.openfoam.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile { version 2.0; format ascii; class dictionary; object combustionProperties; }

combustionModel none;
"""

    g_file = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v13                                   |
|   \\\\  /    A nd           | Website:  www.openfoam.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile { version 2.0; format ascii; class uniformDimensionedVectorField; object g; }

dimensions [0 1 -2 0 0 0 0];
value (0 0 -9.81);
"""

    thermo_props = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v13                                   |
|   \\\\  /    A nd           | Website:  www.openfoam.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile { version 2.0; format ascii; class dictionary; object thermophysicalProperties; }

thermoType
{
    type heRhoThermo;
    mixture multiComponentMixture;
    transport const;
    thermo janaf;
    equationOfState perfectGas;
    specie specie;
    energy sensibleEnthalpy;
}

defaultSpecie Air;
species ( Air CH4 );

    Air { specie { molWeight 28.96; } thermodynamics { Tlow 200; Thigh 6000; Tcommon 1000; highCpCoeffs (3.0879271 0.001245971 -4.237188e-07 6.747207e-11 -3.97077e-15 -995.2627 5.959609); lowCpCoeffs (3.568396 -0.0006787294 1.5537e-06 -3.29937e-12 -4.686e-12 -995.2627 3.67482); } transport { mu 1.8e-05; Pr 0.7; } }
    CH4 { specie { molWeight 16.04; } thermodynamics { Tlow 200; Thigh 6000; Tcommon 1000; highCpCoeffs (0.074851495 0.0133909467 -5.73285809e-06 1.22292535e-09 -1.01815251e-13 -9468.34459 18.437318); lowCpCoeffs (5.14987613 -0.0136709788 4.91800599e-05 -4.84743026e-08 1.66693956e-11 -10246.6476 -4.64130376); } transport { mu 1.8e-05; Pr 0.7; } }
"""

    # --- 2. MESH & CONTROL ---
    surface_features = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v13                                   |
|   \\\\  /    A nd           | Website:  www.openfoam.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile { version 2.0; format ascii; class dictionary; object surfaceFeaturesDict; }

surfaces ( "room.stl" );
extractionMethod extractFromSurface;
includedAngle 150;
subsetFeatures { nonManifoldEdges no; openEdges yes; }
"""

    control_dict = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v13                                   |
|   \\\\  /    A nd           | Website:  www.openfoam.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile { version 2.0; format ascii; class dictionary; object controlDict; }

application foamRun;
solver multicomponentFluid;
startFrom startTime;
startTime 0;
stopAt endTime;
endTime 600;
deltaT 0.0001;
writeControl runTime;
writeInterval 50;
purgeWrite 0;
writeFormat ascii;
writePrecision 6;
writeCompression off;
timeFormat general;
timePrecision 6;
runTimeModifiable true;
adjustTimeStep yes;
maxCo 1.0;
maxDeltaT 0.05;

functions
{
    outletFlow
    {
        type surfaceFieldValue;
        libs ("libfieldFunctionObjects.so");
        writeControl timeStep;
        writeInterval 1;
        log true;
        writeFields true;
        surfaceFormat vtk;
        regionType patches;
        patches ( placeholder_outlet_patch );
        operation sum;
        fields ( phi );
    }
}
"""

    block_mesh = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v13                                   |
|   \\\\  /    A nd           | Website:  www.openfoam.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile { version 2.0; format ascii; class dictionary; object blockMeshDict; }

convertToMeters 1;

vertices
(
    (-0.5 -0.5 -0.5) ( 4.5 -0.5 -0.5) ( 4.5 4.5 -0.5) (-0.5 4.5 -0.5)
    (-0.5 -0.5 4.5) ( 4.5 -0.5 4.5) ( 4.5 4.5 4.5) (-0.5 4.5 4.5)
);

blocks ( hex (0 1 2 3 4 5 6 7) (30 30 30) simpleGrading (1 1 1) );

edges ();

boundary
(
    bounding_box
    {
        type patch;
        faces
        (
            (0 1 5 4) (1 2 6 5) (2 3 7 6) (3 0 4 7) (0 3 2 1) (4 5 6 7)
        );
    }
);

mergePatchPairs ();
"""

    snappy_hex_mesh = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v13                                   |
|   \\\\  /    A nd           | Website:  www.openfoam.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile { version 2.0; format ascii; class dictionary; object snappyHexMeshDict; }

castellatedMesh true; snap true; addLayers false;

geometry { room.stl { type triSurfaceMesh; file "room.stl"; name room; } }

castellatedMeshControls
{
    maxLocalCells 100000; maxGlobalCells 2000000; minRefinementCells 10;
    maxLoadUnbalance 0.10; nCellsBetweenLevels 3; allowFreeStandingZoneFaces true;
    features ( { file "room.eMesh"; level 1; } );
    refinementSurfaces { room { level (2 2); patchInfo { type wall; } } }
    resolveFeatureAngle 30; insidePoint (2.0 2.0 2.0);
    refinementRegions {}
}

snapControls { nSmoothPatch 3; tolerance 2.0; nSolveIter 30; nRelaxIter 5; nFeatureSnapIter 10; implicitFeatureSnap false; explicitFeatureSnap true; multiRegionFeatureSnap false; }

addLayersControls { relativeSizes true; layers {} expansionRatio 1.0; finalLayerThickness 0.3; minThickness 0.1; nGrow 0; featureAngle 60; nRelaxIter 3; nSmoothSurfaceNormals 1; nSmoothNormals 3; nSmoothThickness 10; maxFaceThicknessRatio 0.5; maxThicknessToMedialRatio 0.3; minMedianAxisAngle 90; nBufferCellsNoExtrude 0; nLayerIter 50; }

meshQualityControls { #include "meshQualityDict" }

mergeTolerance 1e-6;
"""

    mesh_quality = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v13                                   |
|   \\\\  /    A nd           | Website:  www.openfoam.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile { version 2.0; format ascii; class dictionary; object meshQualityDict; }

maxNonOrtho 65; maxBoundarySkewness 20; maxInternalSkewness 4; maxConcave 80;
minVol 1e-13; minTetQuality 1e-30; minArea -1; minTwist 0.05; minDeterminant 0.001;
minFaceWeight 0.05; minVolRatio 0.01; minTriangleTwist -1; nSmoothScale 4; errorReduction 0.75;
"""

    fv_options = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v13                                   |
|   \\\\  /    A nd           | Website:  www.openfoam.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile { version 2.0; format ascii; class dictionary; object fvOptions; }

leakSource
{
    type semiImplicitSource;
    active true;
    selectionMode points;
    points ( (placeholder_x placeholder_y placeholder_z) );
    volumeMode absolute;
    injectionRateSuSp
    {
        placeholder_species_name (placeholder_rate_kgs 0);
    }
}
"""

    fv_schemes = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v13                                   |
|   \\\\  /    A nd           | Website:  www.openfoam.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile { version 2.0; format ascii; class dictionary; object fvSchemes; }

ddtSchemes { default Euler; }
gradSchemes { default Gauss linear; grad(p_rgh) Gauss linear; grad(U) cellLimited Gauss linear 1; }
divSchemes
{
    default Gauss limitedLinear 1;
    div(phi,U) Gauss linearUpwind grad(U);
    div(phi,k) Gauss upwind;
    div(phi,epsilon) Gauss upwind;
    div(phi,omega) Gauss upwind;
    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
wallDist { method meshWave; }
"""

    fv_solution = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v13                                   |
|   \\\\  /    A nd           | Website:  www.openfoam.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile { version 2.0; format ascii; class dictionary; object fvSolution; }

solvers
{
    "rho.*"
    {
        solver diagonal;
    }
    "p"
    {
        solver GAMG; tolerance 1e-6; relTol 0.05; smoother DICGaussSeidel;
        agglomerator faceAreaPair; mergeLevels 1; cacheAgglomeration true;
    }
    "pFinal"
    {
        solver GAMG; tolerance 1e-7; relTol 0; smoother DICGaussSeidel;
        agglomerator faceAreaPair; mergeLevels 1; cacheAgglomeration true;
    }
    "p_rgh"
    {
        solver GAMG; tolerance 1e-6; relTol 0.05; smoother DICGaussSeidel;
        agglomerator faceAreaPair; mergeLevels 1; cacheAgglomeration true;
    }
    "p_rghFinal"
    {
        solver GAMG; tolerance 1e-7; relTol 0; smoother DICGaussSeidel;
        agglomerator faceAreaPair; mergeLevels 1; cacheAgglomeration true;
    }
    "(U|h|k|epsilon|omega|K|Yi|Air|N2|O2|CO2|H2O|H2|CH4|C2H6|C3H8|C4H10|C2H4|C2H2|NH3|CO|H2S|SO2|Cl2|He|Ar|C6H6|C7H8|C8H10|C6H14|C5H12|CH3OH|NO2|NO|HCN).*"
    {
        solver PBiCGStab; preconditioner DILU; tolerance 1e-6; relTol 0.1;
    }
}

PIMPLE
{
    momentumPredictor yes;
    nOuterCorrectors 3;
    nCorrectors 1;
    nNonOrthogonalCorrectors 0;
    
    residualControl
    {
        U 1e-4;
        p_rgh 1e-3;
        h 1e-4;
        "Yi.*" 1e-4;
        "(CH4|H2|C3H8|C2H6|C4H10|NH3|Cl2|H2S|CO|CO2|Air|He|O2|N2|Ar|SO2|NO2|NO|HCN).*" 1e-4;
    }
}
"""

    # --- 3. INITIAL CONDITIONS ---
    u_file = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v13                                   |
|   \\\\  /    A nd           | Website:  www.openfoam.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile { version 2.0; format ascii; class volVectorField; object U; }

dimensions [0 1 -1 0 0 0 0];
internalField uniform (0 0 0);

boundaryField
{
    ".*" { type noSlip; }
    bounding_box { type inletOutlet; inletValue uniform (0 0 0); value uniform (0 0 0); }
    placeholder_inlet_patch { type fixedValue; value uniform (1 0 0); }
    placeholder_outlet_patch { type inletOutlet; inletValue uniform (0 0 0); value uniform (0 0 0); }
}
"""

    p_rgh_file = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v13                                   |
|   \\\\  /    A nd           | Website:  www.openfoam.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile { version 2.0; format ascii; class volScalarField; object p_rgh; }

dimensions [1 -1 -2 0 0 0 0];
internalField uniform 101325;

boundaryField
{
    ".*" { type fixedFluxPressure; value uniform 101325; }
    bounding_box { type fixedValue; value uniform 101325; }
    placeholder_outlet_patch { type fixedValue; value uniform 101325; }
}
"""

    p_file = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v13                                   |
|   \\\\  /    A nd           | Website:  www.openfoam.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile { version 2.0; format ascii; class volScalarField; object p; }

dimensions [1 -1 -2 0 0 0 0];
internalField uniform 101325;

boundaryField
{
    ".*" { type calculated; value uniform 101325; }
    bounding_box { type calculated; value uniform 101325; }
}
"""

    t_file = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v13                                   |
|   \\\\  /    A nd           | Website:  www.openfoam.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile { version 2.0; format ascii; class volScalarField; object T; }

dimensions [0 0 0 1 0 0 0];
internalField uniform 300;

boundaryField
{
    ".*" { type zeroGradient; }
    bounding_box { type inletOutlet; inletValue uniform 300; value uniform 300; }
    placeholder_inlet_patch { type fixedValue; value uniform 300; }
}
"""

    k_file = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v13                                   |
|   \\\\  /    A nd           | Website:  www.openfoam.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile { version 2.0; format ascii; class volScalarField; object k; }

dimensions [0 2 -2 0 0 0 0];
internalField uniform 0.01;

boundaryField
{
    ".*" { type kqRWallFunction; value uniform 0.01; }
    bounding_box { type inletOutlet; inletValue uniform 0.01; value uniform 0.01; }
    placeholder_inlet_patch { type fixedValue; value uniform 0.01; }
    placeholder_outlet_patch { type inletOutlet; inletValue uniform 0.01; value uniform 0.01; }
}
"""

    eps_file = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v13                                   |
|   \\\\  /    A nd           | Website:  www.openfoam.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile { version 2.0; format ascii; class volScalarField; object epsilon; }

dimensions [0 2 -3 0 0 0 0];
internalField uniform 0.01;

boundaryField
{
    ".*" { type epsilonWallFunction; value uniform 0.01; }
    bounding_box { type inletOutlet; inletValue uniform 0.01; value uniform 0.01; }
    placeholder_inlet_patch { type fixedValue; value uniform 0.01; }
    placeholder_outlet_patch { type inletOutlet; inletValue uniform 0.01; value uniform 0.01; }
}
"""

    omega_file = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v13                                   |
|   \\\\  /    A nd           | Website:  www.openfoam.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile { version 2.0; format ascii; class volScalarField; object omega; }

dimensions [0 0 -1 0 0 0 0];
internalField uniform 0.01;

boundaryField
{
    ".*" { type omegaWallFunction; value uniform 0.01; }
    bounding_box { type inletOutlet; inletValue uniform 0.01; value uniform 0.01; }
    placeholder_inlet_patch { type fixedValue; value uniform 0.01; }
    placeholder_outlet_patch { type inletOutlet; inletValue uniform 0.01; value uniform 0.01; }
}
"""

    nut_file = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v13                                   |
|   \\\\  /    A nd           | Website:  www.openfoam.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile { version 2.0; format ascii; class volScalarField; object nut; }

dimensions [0 2 -1 0 0 0 0];
internalField uniform 0;

boundaryField
{
    ".*" { type nutkWallFunction; value uniform 0; }
    bounding_box { type calculated; value uniform 0; }
    placeholder_inlet_patch { type calculated; value uniform 0; }
    placeholder_outlet_patch { type calculated; value uniform 0; }
}
"""

    alphat_file = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v13                                   |
|   \\\\  /    A nd           | Website:  www.openfoam.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile { version 2.0; format ascii; class volScalarField; object alphat; }

dimensions [1 -1 -1 0 0 0 0];
internalField uniform 0;

boundaryField
{
    ".*" { type calculated; value uniform 0; }
    bounding_box { type calculated; value uniform 0; }
    placeholder_inlet_patch { type calculated; value uniform 0; }
    placeholder_outlet_patch { type calculated; value uniform 0; }
}
"""

    files = {
        os.path.join(base_dir, "system", "blockMeshDict"): block_mesh,
        os.path.join(base_dir, "system", "snappyHexMeshDict"): snappy_hex_mesh,
        os.path.join(base_dir, "system", "surfaceFeaturesDict"): surface_features,
        os.path.join(base_dir, "system", "meshQualityDict"): mesh_quality,
        os.path.join(base_dir, "system", "controlDict"): control_dict,
        os.path.join(base_dir, "system", "fvOptions"): fv_options,
        os.path.join(base_dir, "system", "fvSchemes"): fv_schemes,
        os.path.join(base_dir, "system", "fvSolution"): fv_solution,
        os.path.join(base_dir, "constant", "thermophysicalProperties"): thermo_props,
        os.path.join(base_dir, "constant", "thermophysicalTransport"): thermophysical_transport,
        os.path.join(base_dir, "constant", "combustionProperties"): combustion_properties,
        os.path.join(base_dir, "constant", "g"): g_file,
        os.path.join(base_dir, "constant", "momentumTransport"): mom_transport,
        os.path.join(base_dir, "0", "U"): u_file,
        os.path.join(base_dir, "0", "p_rgh"): p_rgh_file,
        os.path.join(base_dir, "0", "p"): p_file,
        os.path.join(base_dir, "0", "T"): t_file,
        os.path.join(base_dir, "0", "k"): k_file,
        os.path.join(base_dir, "0", "epsilon"): eps_file,
        os.path.join(base_dir, "0", "omega"): omega_file,
        os.path.join(base_dir, "0", "nut"): nut_file,
        os.path.join(base_dir, "0", "alphat"): alphat_file,
    }

    for path, content in files.items():
        with open(path, "w") as f:
            f.write(content)
        print(f"Written: {path}")

    print("✅ Template updated: Base files generated.")

if __name__ == "__main__":
    create_base_template()
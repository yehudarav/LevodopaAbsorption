    singlePhaseTransportModel laminarTransport(U, phi);

    // Thermal expansion coefficient [1/K] 1/Theta0
    dimensionedScalar beta(laminarTransport.lookup("beta"));

    // Laminar Prandtl number
    dimensionedScalar Pr(laminarTransport.lookup("Pr"));

    // Turbulent Prandtl number
    dimensionedScalar Prt(laminarTransport.lookup("Prt"));

	// Read the fcoriolis. if not present set to 0. 
	dimensionedVector zerocoriolis("zerovec1",dimensionSet(0, 0, -1, 0, 0, 0, 0),vector(0,0,0));
	dimensionedVector fCoriolis(laminarTransport.lookupOrDefault("fCoriolis",zerocoriolis));
	

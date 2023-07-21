# magnetism-app
A Shiny app to demonstrate the Ising model of magnetism

The Ising model is a simple model of how a lattice of magnetic particles, or "spins", can form complex ordering via a simple interaction.
In this model, spins only interact with their nearest neighbors, and prefer to either align or anti-align with neighboring spins for ferromagnetic and antiferromagnetic materials, respectively. The model allows spins to slip if it is energetically favorable for them to do so, or if they can get enough random energy based on the temperature of the simulated material. This randomness allows spins to flip into a higher energy state with a probability given by exp(-1/T), where T is the temperature measured in units of the interaction energy between neighboring spins.

This app lets the user explore the Ising model by creating a lattice of spins (either ferromagnetic or antiferromagnetic), setting a temperature, and an external magnetic field. The app simulates how the lattice would evolve from a random initial state. The user can either evolve the lattice one step at a time (a "step" allows all the spins to have a chance to flip once), multiple steps at a time, or to run continuously so that the lattice will constantly evolve. The current state of the lattice is shown at the top of the main panel, and a time series plot of the lattice's total magnetic moment (the sum of all the spins, with +1 for spin up and -1 for spin down) is shown at the bottom of the main panel.

The user is invited to play around with the temperature, the type of interaction, the external field, and the lattice size and see how complex magnetic order that spans the whole lattice (such as domains and boundaries) can form even in this simple, local model!

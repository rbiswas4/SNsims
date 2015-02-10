Background:
-----------
We already have an object that represents a broad category of SN, viz SNCosmo.Model . A specific SN with specified properties (ie. model parameters) can be instantiated by the use of the sncosmo.Model, and sncosmo.Model.set functions. 

We would like to use SN generated using this method as a building block towards simulating a SN survey. The discussion here is about supporting different 
usecases for simulation using such building blocks, ie. we assume that individual SN can be simulated using sncosmo.model .

Different types of simulations:
-------------------------------
There are two types of details in simulations to consider:
- details with regards to model of SN
- details with regards to modelling the observations, and therefore outputs.

A list of possible types of SN Model issues in simulations we may want to support:

1. Volumetric SN simulation: Basic inputs: number of SN per comoving unit volume, per rest frame time as a function of redsfhift, redshift dependent distributions of model parameters, observational parameters. 
2.  Based on a galaxy catalog: A set of galaxies with properties at different redshifts. Rules about SN forming in galaxies, and/or a rate of SN as in volumetric  simulations.
3. Given a set of SN (perhapse simulated along the lines of (1.)), associate galaxies to the SN using a set of rules.
4. Simulate SNe at specific redshift(s).
5. Simulate SNe without noise realization.
6. Only simulate "vanilla" supernovae, i.e. an average SN Ia with no Galactic or galactic dust extinction.

Other use cases would be interesting to learn about. It is not entirely clear that having both (2.) and (3.) are essential. Two applications that may requre having each of these are 
- (For 3): Imagine you already have a set of SNIa from some survey, and you want a distribution of galaxies consistent with the SNIa you already have and a set of rules.
- (For 2): If you want a generative model, after which you would like to change survey strategy/selection obtaining different sets of SN in each case, 3 is a better option.

What kind of output would be desirable?:
---------------------------------------
The minimal output is a set of summary parameters (metadata) that describe the characteristics of the SN, along with the simulation truths of the supernova, as well as the photometry data (MJD, Filter, Flux, ZeroPoint, Magsys)

In a simulation, we want 
- potentially multiple models of SN. eg. Type Ia and different core-collapse  models. More generally, we may want to use even multiple models of Type Ia (eg. prompt and tardy components.)
- Enable easy draws of parameter of sncosmo.model 
- Assert that models be drawn with certain frequencies. In practice, we need multiple ways of doing this. 
        - The most frequently used method will be setting the probabilities through a rate. eg.  for a volumetric simulation of a survey of time length T, we will simulate V(z) * rate * (1 +z) SNIa with t0 drawn randomly for times of the simulation. 
        - A second method to think of is two species of SN, which together match a rate, but we would like to assert that ratio of their numbers is a:b 

What implementations should the abstract class provide?
-------------------------------------------------------

- Convert rates to an expected number of SN in redshift and time bin and angular bin
- Provide functions to sample a ra, dec bin area uniformly
-  

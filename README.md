
# Contact Process Ising Model
The Contact Process Ising Model (CPIM) picks up on the Contact Process (CP)

	https://github.com/jekeymer/Contact-Process/wiki 

as a model of 2D spatial growth and mix it to the Ising Model (IM), 
a model for the magnetic interaction.

The CPIM stand for a Lattice model mixing the Contact Process (CP)

	https://github.com/jekeymer/Contact-Process

and the Ising Model (IM)

	https://en.wikipedia.org/wiki/Ising_model

together in a combined model. Even thoughh Metropolis and Gillespie are not compatible.

The basic idea is to merge both models in order to model the spatial biology of an Ising-like 
genetic network expressed in a quasi-2D colony of bacteria
growing on a surface of solid agar. 

The genetic network consist of a bi-stable genetic system of 
reporter genes (RFP vs GFP) coupled by auto-inducer molecules (HSLs)

By fair we mean here using the Gillespie algorithm

	https://en.wikipedia.org/wiki/Gillespie_algorithm

in the Monte Carlo method updating the Lattice. 

This is a fundamental requirement as there are Lattice states with two possible reactions.
Of these only differentiation is commensurable as it does represent a rate. Metropolis flips are not.

Lattice states are: 
	0: vacancy, representing empty sites suitable for colonization by its neightboors (plotted in black) 
	2: undifferentiated state (plotted in white) 
	-1: spin down, representing a site expressing one reporter gene (plotted in green) 
	+1: spin up, representing a site expressing the other reporter (plotted in magenta)

Occupied sites (states: -1,2,+1) can have spin flips or differentiation reactions 
together with particle death reactions. 

So we use only Gillespie (version Partial Gillespie) for the differentiation reaction function 
	update_lattice_2
For the spin flisp, we only do it if the cell survides the time step. Thus no Gillespie.

To compare what happens otherwise, we use Gillespie for flips and see that the Metropolis Algorithm is
not compatible with Gillespie as it is not suppose to represent rates. 
We see then (if ran) it affects the contact process death rate in an artificial fashion.
For this test we made the function:
	update_lattice_1

In the Gillespie algorithm we use 2 random numbers. 
The first one is used to decide which one of the competing reactions might take place. 
Then, a second random number is used to decide if the chosen reaction, indeed take place (or not). 
This way, the Poisson process is respected and probabilities are well defined.


To compile type:
	  make

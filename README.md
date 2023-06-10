
# Contact Process Ising Model
The Contact Process Ising Model (CPIM) picks up on the Contact Process (CP)

	https://github.com/jekeymer/Contact-Process/wiki 

as a model of 2D spatial growth and mix it to the Ising Model (IM), 
a model for the magnetic interaction.

The letters CPIM stand for Contact Process (CP)

	https://github.com/jekeymer/Contact-Process

together with the Ising Model (IM)

	https://en.wikipedia.org/wiki/Ising_model


Even though Metropolis and Gillespie are not fully compatible, here we experiment with the idea of using both together. Two versions are tested: Full Gillespie and Partial Gillespie.

The basic idea is to merge both models in order to model the spatial biology of an Ising-like (synthetic) genetic network expressed in a quasi-2D colony of bacteria
growing on a surface of solid agar. 

The genetic network consist of a bi-stable genetic system (switch) of two
reporter genes (RFP vs GFP) which are spatially coupled by diffusive auto-inducer molecules (HSLs) which act as quoromones (hormones)as they are know from the biology of quorum sensing.

	https://en.wikipedia.org/wiki/Quorum_sensing


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

Occupied sites can be in states: {-1,2,+1}. Transformations of these states can ocurr due to spin flip or differentiation reactions. Birth and death process lead to colonization together with particle death reactions which are incorporated in the Contact process. Thus the full system of states of the process {0, -1, +1, 2}.

METROPOLIS vs GILLESPIE

So we use only Gillespie (version Partial Gillespie) for the differentiation reaction function. This is implemented in the function.
	
	update_lattice_2

For the spin flips, we only ran the Metropolis algorithm only if the site visited survives the CP Monte Carlo step. Thus no Gillespie is used here!

To compare what happens otherwise (using Gillespie), we also implement a version of the update lattice code using the Gillespie algorith, for flips and see how the Metropolis Algorithm is not compatible with Gillespie. 

	https://en.wikipedia.org/wiki/Metropolis%E2%80%93Hastings_algorithm

This is due to the fact that Metropolis is not suppose to represent true rates. It is only meant to draw configurations consistent with the ensamble average defined by the Isin Model's statistical mechanics.

We see that if Gillespie is ran along with the Metropolis algorithm, the flipping operator affects the rates of the Contact Process in an artificial fashion. To see this the interested hacker can test by herself. For this, we made the function:

	update_lattice_1

which can be ran by changing the following line at the top of the code

#define GILLESPIE_OPTION 2 // 1 : full ; 2: partial

GILLESPIE

In the Gillespie algorithm we use 2 random numbers. The first one is used to decide which one of the competing reactions might take place. Then, a second random number is used to decide if the chosen reaction, indeed take place (or not). This way, the Poisson process is respected and probabilities are well defined.

COMPILE

To compile just type:
	  make

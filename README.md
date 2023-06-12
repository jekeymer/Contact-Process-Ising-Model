
# Contact Process Ising Model
The Contact Process Ising Model (CPIM) picks up on the Contact Process (CP)

	https://github.com/jekeymer/Contact-Process/wiki 

as a model of 2D spatial growth and mix it to the Ising Model (IM), 
a model for the magnetic interaction.

The letters CPIM stand for Contact Process (CP)

	https://github.com/jekeymer/Contact-Process

together with the Ising Model (IM)

	https://en.wikipedia.org/wiki/Ising_model


The basic idea is to merge both models in order to model the spatial biology of an Ising-like (synthetic) genetic network expressed in a quasi-2D colony of bacteria
growing on a surface of solid agar. 

Notice that as cells are macroscopic object, we separate two different time-scales and organizational (spatial) scales. 

(i) one is the scale of the full cell, which undergoes birth, death, colonization process. We model this by the CP. 
(ii) another scale is the scale of chemical reactions withing living cells. These reactions are only considered possible as long as the cell is alive. That is it survives the Monte Carlo step of the CP. Thus, spin flips of the IM & differentiation (gentic induction) only take place once a cell survive. 

The genetic network we are interested in modelling is made of a bi-stable genetic switch consisting of two reporter genes (RFP vs GFP) which are spatially coupled by diffusive auto-inducer molecules (HSLs) which act as quoromones (hormones) as they are know from the biology of quorum sensing.

	https://en.wikipedia.org/wiki/Quorum_sensing


Lattice states of the model are: 

	0: vacancy, representing empty sites suitable for colonization by its neightboors (plotted in black) 

	2: undifferentiated state (plotted in white) 

	-1: spin down, representing a site expressing one reporter gene (plotted in green) 

	+1: spin up, representing a site expressing the other reporter (plotted in magenta)

Occupied sites can be in states: {-1,2,+1}. Transformations of these states can ocurr due to spin flip or differentiation reactions. Birth and death process lead to colonization together with particle death reactions which are incorporated in the Contact process. Thus the full system of states of the process is {0, -1, +1, 2}.


For the spin flips, we only ran the Metropolis algorithm only if the site visited survives the CP Monte Carlo step. 

	https://en.wikipedia.org/wiki/Metropolis%E2%80%93Hastings_algorithm

Notice that Metropolis is not suppose to represent true rates. It is only meant to draw configurations consistent with the ensamble average defined by the Isin Model's statistical mechanics. As the temporal scale of gene induction is assumed to be bigger than birth death process we mantain the Metropolis algorithm to trigger spin flips presenting genetic induction.

COMPILE

To compile the code just type:
	  make

or use gcc and the Gtk configuration tool by typing:

	 gcc CPIM.c mt64.c -lm -o CPIM `pkg-config --cflags gtk+-3.0` `pkg-config --libs gtk+-3.0`



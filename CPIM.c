// C code to expand the Contact Process (CP) code to incorporate 
// a basic Ising Model into the occupancy states of the CP

#include <stdlib.h>
#include <gtk/gtk.h> /* GUI, Gtk library */
#include "mt64.h"    /* Pseudo-random number generation MT library (64 bit) */
#include <math.h>    /* Math to transform random n from continuous to discrete */
#include <time.h>    /* Used to seed pseudo-random number generator */
#include <stdio.h>

/* Lattice Size */
#define X_SIZE 256
#define Y_SIZE 256

/* Defaulfs */
#define SAMPLE_RATE 100
// default birth/colonization rate/probability and scale ranges
#define BETA  0.003
#define BETA_STEP 0.00001
#define BETA_MIN 0.000000
#define BETA_MAX 0.005
// default mortality/extinction rate/probability and scale ranges
#define DELTA  0.0001
#define DELTA_STEP  0.00001
#define DELTA_MIN 0.000000
#define DELTA_MAX 0.005
// default differentiation rate/probability and scale ranges
#define ALPHA  0.1
#define ALPHA_STEP 0.01
#define ALPHA_MIN  0.00
#define ALPHA_MAX  1.0
// strength of the coupling in positive terms (J =  -1*COUPLING kBT units)
// we should have 1/2 if we do not want to douple count pairs
// therefore
// use a positive number!
#define COUPLING (1)
// default Temperature and scale ranges
#define TEMPERATURE 2.269
#define TEMPERATURE_STEP 0.0000001
#define TEMPERATURE_MIN  0.0000001
#define TEMPERATURE_MAX  15 
// default interaction radius
#define RADIUS 1
// default initial condition chosen
#define INIT 1



/* Structure with the simulation data */
struct simulation
  {
  int lattice_configuration[X_SIZE][Y_SIZE]; /* Store latice configuration */
  gint run;                   /* Time handler tag */
  gboolean running;           /* Are we running? */
  int init_option;            /* Choice of initial condition*/
  gboolean initialized;       /* Have we been initialized? */
  int generation_time;        /* Generations simulated */
  int Ising_neighboorhood;    /* Ising Neighboorhood: r=1 (NN) vs r=2 (NNN)*/
  int occupancy;              /* Lattice occupancy */
  int vacancy;                /* Lattice vacancy*/
  int up;                     /* Number of spins in the up   (+1) state */
  int down;                   /* Number of spins in the down (-1) state */
  int display_rate;           /* Display rate: to paint the lattice*/
  double birth_rate;          /* Contact Process' birth */
  double death_rate;          /* Contact Process' death */
  double differentiation_rate;/* Differentiation into spin state */
  double T;                   /* Ising's temperature */
  double J;                   /* Ising's coupling: ferro (-kB) or anti-ferro (+kB) */
  double lamda_rate;          /* Contact-Ising Monte Carlo biass*/
} s ;        // instance s of the structure to hold the simulation


// Gdk Pixel Buffer functions Implemented at the end of document.
void put_pixel(GdkPixbuf *pixbuf, int x, int y, guchar red, guchar green, guchar blue, guchar alpha);
static void paint_a_background (gpointer data);
static void paint_lattice (gpointer data);


// Other Funtions

double local_energy (int x, int y)
	{
  // Energy of site at coordinate (x,y)
  double energy; //in kB*T units
	int up = 0;   
	int down = 0; 
	if (s.Ising_neighboorhood == 1)  // Nearest Neighboorhood (NN) has 4 sites
    	{ 
      // we check the South (S) neighboor (#1)
      if (s.lattice_configuration[x][(int)((Y_SIZE + y+1)%Y_SIZE)] == 1){up++;}
    	 else if (s.lattice_configuration[x][(int)((Y_SIZE + y+1)%Y_SIZE)] == -1){down++;}
    	// we check the North (N) neighboor (#2)
      if (s.lattice_configuration[x][(int)((Y_SIZE + y-1)%Y_SIZE)] == 1){up++;}
    	 else if (s.lattice_configuration[x][(int)((Y_SIZE + y-1)%Y_SIZE)] == -1){down++;}
    	// we check the West (W) neighboor  (#3)
      if (s.lattice_configuration[(int)((X_SIZE + x-1)%X_SIZE)][y] == 1){up++;}
    	 else if (s.lattice_configuration[(int)((X_SIZE + x-1)%X_SIZE)][y] == -1){down++;}
    	// we check the East (E) neighboor  (#4)
      if (s.lattice_configuration[(int)((X_SIZE + x+1)%X_SIZE)][y] == 1){up++;}
    	 else if (s.lattice_configuration[(int)((X_SIZE + x+1)%X_SIZE)][y] == -1){down++;}
    	}
    else if (s.Ising_neighboorhood == 2) //Next Nearest Neighboorhood (NNN) has 12 sites
        { 
        // we check the South (S) neighboor (#1)
        if (s.lattice_configuration[x][(int)((Y_SIZE + y+1)%Y_SIZE)] == 1){up++;}
    	   else if (s.lattice_configuration[x][(int)((Y_SIZE + y+1)%Y_SIZE)] == -1){down++;}
    	  // we check the North (N) neighboor (#2)
        if (s.lattice_configuration[x][(int)((Y_SIZE + y-1)%Y_SIZE)] == 1){up++;}
    	   else if (s.lattice_configuration[x][(int)((Y_SIZE + y-1)%Y_SIZE)] == -1){down++;}
    	  // we check the West (W) neighboor (#3)
        if (s.lattice_configuration[(int)((X_SIZE + x-1)%X_SIZE)][y] == 1){up++;}
    	   else if (s.lattice_configuration[(int)((X_SIZE + x-1)%X_SIZE)][y] == -1){down++;}
    	  // we check the East (E) neighboor (#4)
        if (s.lattice_configuration[(int)((X_SIZE + x+1)%X_SIZE)][y] == 1){up++;}
    	   else if (s.lattice_configuration[(int)((X_SIZE + x+1)%X_SIZE)][y] == -1){down++;}
        // we check the South-South (SS) neighboor (#5)
        if (s.lattice_configuration[x][(int)((Y_SIZE + y+2)%Y_SIZE)] == 1){up++;}
    	   else if (s.lattice_configuration[x][(int)((Y_SIZE + y+2)%Y_SIZE)] == -1){down++;}
    	  // we check the North-North (NN) neighboor (#6)
        if (s.lattice_configuration[x][(int)((Y_SIZE + y-2)%Y_SIZE)] == 1){up++;}
    	   else if (s.lattice_configuration[x][(int)((Y_SIZE + y-2)%Y_SIZE)] == -1){down++;}
    	  // we check the West-West (WW) neighboor   (#7)
        if (s.lattice_configuration[(int)((X_SIZE + x-2)%X_SIZE)][y] == 1){up++;}
    	   else if (s.lattice_configuration[(int)((X_SIZE + x-2)%X_SIZE)][y] == -1){down++;}
    	  // we chack the East-East (EE) neighboor   (#8)
        if (s.lattice_configuration[(int)((X_SIZE + x+2)%X_SIZE)][y] == 1){up++;}
    	   else if (s.lattice_configuration[(int)((X_SIZE + x+2)%X_SIZE)][y] == -1){down++;}
    	  // we check the South-West (SW) neighboor  (#9)
        if (s.lattice_configuration[(int)((X_SIZE + x-1)%X_SIZE)][(int)((Y_SIZE + y+1)%Y_SIZE)] == 1){up++;}
    	   else if (s.lattice_configuration[(int)((X_SIZE + x-1)%X_SIZE)][(int)((Y_SIZE + y+1)%Y_SIZE)] == -1){down++;}
    	  // we check the North-East (NE) neighboor  (#10)
        if (s.lattice_configuration[(int)((X_SIZE + x+1)%X_SIZE)][(int)((Y_SIZE + y-1)%Y_SIZE)] == 1){up++;}
    	   else if (s.lattice_configuration[(int)((X_SIZE + x+1)%X_SIZE)][(int)((Y_SIZE + y-1)%Y_SIZE)] == -1){down++;}
    	  // we check the North-West (NW) neighboor  (#11)
        if (s.lattice_configuration[(int)((X_SIZE + x-1)%X_SIZE)][(int)((Y_SIZE + y-1)%Y_SIZE)] == 1){up++;}
    	   else if (s.lattice_configuration[(int)((X_SIZE + x-1)%X_SIZE)][(int)((Y_SIZE + y-1)%Y_SIZE)] == -1){down++;}
       	// we check the South-East (SE) neighboor  (#12)
        if (s.lattice_configuration[(int)((X_SIZE + x+1)%X_SIZE)][(int)((Y_SIZE + y+1)%Y_SIZE)] == 1){up++;}
    	   else if (s.lattice_configuration[(int)((X_SIZE + x+1)%X_SIZE)][(int)((Y_SIZE + y+1)%Y_SIZE)] == -1){down++;}
    	  }
	energy =  s.J * (double) (s.lattice_configuration[x][y] * (up-down));
  return energy;
	}



/* Update function */
int update_lattice (gpointer data)
  {
  // int random_neighbor;
  int random_neighbor_state, random_neighbor;
  double random_spin;
  // Energies
  double spin_energy, spin_energy_diff;
  // Probability of reactions
  double transition_probability;
  int random_x_coor, random_y_coor;
  // For the Contact Process we always consider NN interactions
  for (int site = 0; site < (int) (Y_SIZE*X_SIZE); site++)
    {
    /* Pick a random focal site */
    random_x_coor = (int) floor (genrand64_real1 ()* X_SIZE);
    random_y_coor = (int) floor (genrand64_real1 ()* Y_SIZE);
    switch (s.lattice_configuration[random_x_coor][random_y_coor])
      {
      case 0: /* Site is empty */
      /* Chose a random neighbor from the num_neighbors posible ones */
      random_neighbor = (int) floor (genrand64_real3()* 4);
			switch(random_neighbor)
					{
					case 0: // South
							random_neighbor_state = s.lattice_configuration[random_x_coor][(int) ((Y_SIZE + random_y_coor-1)%Y_SIZE)]; 
							break;
					case 1: // North
							random_neighbor_state =	s.lattice_configuration[random_x_coor][(int) ((Y_SIZE + random_y_coor+1)%Y_SIZE)]; 
							break;
					case 2: // East
							random_neighbor_state =	s.lattice_configuration[(int) ((X_SIZE + random_x_coor-1)%X_SIZE)][random_y_coor];
							break;
					case 3: // West
							random_neighbor_state =	s.lattice_configuration[(int) ((X_SIZE + random_x_coor+1)%X_SIZE)][random_y_coor];
							break;
					}
        /* If its random neighbor is occupied: put a copy at the focal site
           with probability brith_rate * dt */
        if (genrand64_real2 () < s.birth_rate)
           {
           switch(random_neighbor_state)
             {
              case 2: 
                s.lattice_configuration[random_x_coor][random_y_coor] = 2;
                s.occupancy ++; s.vacancy --;
               break;
              case 1: 
                s.lattice_configuration[random_x_coor][random_y_coor] = 1;
                s.occupancy ++; s.vacancy --;
                s.up ++;
               break;
              case -1:
                s.lattice_configuration[random_x_coor][random_y_coor] = -1;
                s.occupancy ++;s.vacancy --;
                s.down ++;
               break; 
              case 0:
                s.lattice_configuration[random_x_coor][random_y_coor] = 0;
                break;
             }
          }
        break; /* break case 0 */
      case 2: /* Focal point is in the occupied, undifferentiated state */
        // First we check if the site survives
        // No need for Gillespie as cells are macroscopic compare to its 
        // inner components which can undertake reactions only if the cell
        // indeed exists
        if (genrand64_real2 () < s.death_rate)
                       {
                        s.lattice_configuration[random_x_coor][random_y_coor] = 0;
                        s.occupancy --; s.vacancy ++;
                       }
             else if (genrand64_real2 () < s.differentiation_rate)
                      {
                       /* Set an occupied site in the middle of the lattice */
                       random_spin = (int) ((genrand64_int64 () % 2) * 2) - 1;
                      if (random_spin == 1)
                          {
                           s.lattice_configuration[random_x_coor][random_y_coor] = random_spin;
                           s.up ++;
                           }
                       else if (random_spin == -1)
                           {
                           s.lattice_configuration[random_x_coor][random_y_coor] = random_spin;
                           s.down ++;
                           }
                      }
        break;
      case 1: /* Focal point is in the up (+1) state */
        // We skip Gillespie because of separation of scales
        spin_energy = local_energy (random_x_coor, random_y_coor);
        spin_energy_diff = -(2) * spin_energy;
        transition_probability = exp (-spin_energy_diff/s.T);
        if (genrand64_real2 () < s.death_rate)
                        {
                        s.lattice_configuration[random_x_coor][random_y_coor] = 0;
                        s.occupancy --; s.vacancy ++;
                        s.up --;
                        }
                else if (spin_energy_diff < 0 ||
                                              genrand64_real2 () < transition_probability)
                        {
                        s.lattice_configuration[random_x_coor][random_y_coor] = -1;
                        s.up --;
                        s.down ++;
                        }
        break;
      case -1: /* Focal point is in the down (-1) state */
        // We skip Gillespie because of separation of scales
        spin_energy = local_energy (random_x_coor, random_y_coor);
        spin_energy_diff = -(2) * spin_energy;
        transition_probability = exp (-spin_energy_diff/s.T);
        if (genrand64_real2 () < s.death_rate)
                        {
                        s.lattice_configuration[random_x_coor][random_y_coor] = 0;
                        s.occupancy --; s.vacancy ++;
                        s.down --;
                        }
                else if (spin_energy_diff < 0 ||
                                              genrand64_real2 () < transition_probability)
                        {
                        s.lattice_configuration[random_x_coor][random_y_coor] = 1;
                        s.up ++;
                        s.down --;
                        }
        break;
      }
    }
  s.generation_time ++;
  if(s.generation_time%s.display_rate == 0)
      {
      paint_lattice (data);
      g_print ("Gen: %d \t Vacancy: %f \t Occupancy: %f \t Up: %f \t Down: %f\n",
           s.generation_time, (double) s.vacancy / (double) (Y_SIZE*X_SIZE), (double) s.occupancy/(double) (Y_SIZE*X_SIZE), (double) s.up/(double) (s.occupancy), (double) s.down/(double) s.occupancy);
     }
  return 0;
}



/* Time handler to connect update function to the gtk loop */
gboolean time_handler (gpointer data)
   { 
    update_lattice (data);
    return TRUE;
    }


/* Callback to initialize lattice*/
static void init_lattice (GtkWidget *widget, gpointer data)
  {
  int random_spin;
  int x,y;
  /* Fill the lattice with 0s (unoccupied state) */
  for (x = 0; x < X_SIZE; x++)
    {
    for (y = 0; y < Y_SIZE; y++)
      {
      s.lattice_configuration[x][y]= 0;
      }
    }
  s.occupancy = 0;
  s.up =  0;
  s.down = 0;
  s.vacancy = (int) X_SIZE*Y_SIZE;
  switch(s.init_option)
    {
      case 1:
            /* Set an occupied site in the middle of the lattice */
            random_spin = (int) ((genrand64_int64 () % 2) * 2) - 1;
            if (random_spin == 1)
              {
              s.lattice_configuration[(int) X_SIZE/2][(int) Y_SIZE/2] = random_spin;
              s.up ++; s.vacancy--; s.occupancy++;
              }
            else if (random_spin == -1)
              {
              s.lattice_configuration[(int) X_SIZE/2][(int) Y_SIZE/2] = random_spin;
              s.down ++; s.vacancy --; s.occupancy++;
              }
            break;
      case 2:
            /* Set an undifferentiated site in the middle of the lattice*/
              s.lattice_configuration[(int) X_SIZE/2][(int) Y_SIZE/2] = 2;
              s.vacancy--; s.occupancy++;

            break;
      case 3:
            // Set a small (r=2) cluster with undifferentiated sites in the middle of the lattice
           for (x = (int) X_SIZE/2 - 2 ; x < (int) X_SIZE/2 + 2; x++)
                                for (y = (int) X_SIZE/2 - 2; y < (int) X_SIZE/2 + 2; y++)
                                        {
                                        s.lattice_configuration[x][y]=2;
                                        s.occupancy ++; s.vacancy --;
                                        }
                        break;
            break;
      case 4:
            for (x = (int) X_SIZE/2 - 2 ; x < (int) X_SIZE/2 + 2; x++)
                                for (y = (int) X_SIZE/2 - 2; y < (int) X_SIZE/2 + 2; y++)
                                    {
                                      random_spin = (int) ((genrand64_int64 () % 2) * 2) - 1;
                                      if (random_spin == 1)
                                        {
                                        s.lattice_configuration[x][y] = random_spin;
                                        s.up ++; s.vacancy--; s.occupancy++;
                                        }
                                        else if (random_spin == -1)
                                           {
                                           s.lattice_configuration[x][y] = random_spin;
                                           s.down ++; s.vacancy --; s.occupancy++;
                                           }
                                    }
         
            break;
      case 5:
            // Se a lattice fully occupied with undufferenciated particels
            for (x = 0; x < (int) X_SIZE; x++)
               for (y = 0; y < (int) Y_SIZE; y++)
                    {
                    s.lattice_configuration[x][y]=2;
                    s.occupancy ++; s.vacancy --;
                    }
            break;
    }
   s.initialized = TRUE;
   s.generation_time = 0;
   paint_lattice (data);
   g_print ("Lattice initialized\n");
  }

// Stop simulation control
static void stop_simulation (gpointer data)
  {
  if (s.running)
    {
    g_source_remove (s.run);
    s.running = FALSE;
    g_print ("Simulation stopped\n");
    }
  }


/* Callback to launch dialog with info (github's wiki) */
static void on_button_show_about(GtkWidget *widget, gpointer data)
  {
   GdkPixbuf *pixbuf = gdk_pixbuf_new_from_file("X-Institute_logo_small.tif", NULL);
   GtkWidget *dialog = gtk_about_dialog_new();
   gtk_about_dialog_set_program_name (GTK_ABOUT_DIALOG(dialog),
                                    "Contact Process Ising Model App");
   gtk_about_dialog_set_version(GTK_ABOUT_DIALOG(dialog), "version 0.1, 2023");
   gtk_about_dialog_set_copyright(GTK_ABOUT_DIALOG(dialog),"Open Source Code");
   gtk_about_dialog_set_comments(GTK_ABOUT_DIALOG(dialog),
     "The Contact Process Ising Model (CPIM).");
   gtk_about_dialog_set_website(GTK_ABOUT_DIALOG(dialog),
     "https://github.com/jekeymer/Contact-Process-Ising-Model/wiki");
   gtk_about_dialog_set_logo(GTK_ABOUT_DIALOG(dialog), pixbuf);
   g_object_unref(pixbuf), pixbuf = NULL;
   gtk_dialog_run(GTK_DIALOG (dialog));
   gtk_widget_destroy(dialog);
  }


/* Callback to start simulation */
static void on_button_start_simulation (GtkWidget *button, gpointer data)
  {
  if(!s.running && s.initialized)
    {
    s.run = g_idle_add ((GSourceFunc) time_handler, GTK_IMAGE (data));
    s.running = TRUE;
    g_print ("Simulation started\n");
    }
  }



/* Callback to stop simulation */
static void on_button_stop_simulation (GtkWidget *button, gpointer data)
  {
  stop_simulation (data);
  }


/* Callback to change Initial conditions -- dirty */
/* get_active() method is cleaner as I could use only one handler */
// Init 1
static void on_radio_initial_condition_1 (GtkWidget *button, gpointer data)
  {
  char *id_radio = (char*)data;g_print("%s\n", id_radio);
  s.init_option = 1;
  }
// Init 2
static void on_radio_initial_condition_2 (GtkWidget *button, gpointer data)
  {
  char *id_radio = (char*)data;g_print("%s\n", id_radio);
  s.init_option = 2;
  }
// Init 3
static void on_radio_initial_condition_3 (GtkWidget *button, gpointer data)
  {
  char *id_radio = (char*)data;g_print("%s\n", id_radio);
  s.init_option = 3;
  }
// Init 4
static void on_radio_initial_condition_4 (GtkWidget *button, gpointer data)
  {
  char *id_radio = (char*)data;g_print("%s\n", id_radio);
  s.init_option = 4;
  }
// Init 5
static void on_radio_initial_condition_5 (GtkWidget *button, gpointer data)
  {
  char *id_radio = (char*)data;g_print("%s\n", id_radio);
  s.init_option = 5;
  }

/* Callback to change Ising NN (r=1) vs NNN (r=2) conditions -- dirty */
/* get_active() methosh is cleaner as I could use only one handler */
// NN; r = 1
static void on_radio_NN (GtkWidget *button, gpointer data)
  {
  char *id_radio = (char*)data;g_print("%s\n", id_radio);
  s.Ising_neighboorhood = 1;
  }
// NNN; r =2
static void on_radio_NNN (GtkWidget *button, gpointer data)
  {
  char *id_radio = (char*)data;g_print("%s\n", id_radio);
  s.Ising_neighboorhood = 2;
  }

/* Callback to change Ising J = -kB (ferro) vs J = +kB (anti-ferro) -- dirty */
/* get_active() methosh is cleaner as I could use only one handler */
// Ferro:       J = -kB
static void on_radio_ferro (GtkWidget *button, gpointer data)
  {
    char *id_radio = (char*)data;g_print("%s\n", id_radio);
    s.J = -1 * (float) COUPLING;
    }
// anti Ferro:  J = +kB
static void on_radio_anti_ferro(GtkWidget *button, gpointer data)
  {
    char *id_radio = (char*)data;g_print("%s\n", id_radio);
    s.J =  1 * (float) COUPLING;
    }


/*  Callback to respond Gtk scale slide move event */
static void birth_rate_scale_moved (GtkRange *range, gpointer user_data)
  {
  gdouble pos = gtk_range_get_value (range);
  s.birth_rate = (float) pos;
  }


/*  Callback to respond Gtk scale slide move event */
static void differenciation_rate_scale_moved (GtkRange *range, gpointer user_data)
  {
  gdouble pos = gtk_range_get_value (range);
  s.differentiation_rate = (float) pos;;
  }


/*  Callback to respond Gtk scale slide move event */
static void death_rate_scale_moved (GtkRange *range, gpointer user_data)
  {
  gdouble pos = gtk_range_get_value (range);
  s.death_rate = (float) pos;
  }


/*  Callback to respond Gtk scale slide move event */
static void temperature_scale_moved (GtkRange *range, gpointer user_data)
  {
  gdouble pos = gtk_range_get_value (range);
  s.T = (float) pos;
  }


/*  Callback to respond Gtk scale slide move event */
static void display_rate_scale_moved (GtkRange *range, gpointer user_data)
  {
  gdouble pos = gtk_range_get_value (range);
  s.display_rate = (int) pos;
  }


static void initialize_simulation(void)
{
  /* Initialize Mersenne Twister algorithm for random number genration */
  unsigned int seed = (unsigned int) time (NULL);
  init_genrand64 (seed);

  /* Set default parameters of the simulation */
  //initial condition option
   s.init_option = (int) INIT;

  // Contact Process
  s.birth_rate = (double) BETA;
  s.death_rate = (double) DELTA;

  // Cell differenciation
  s.differentiation_rate = (double) ALPHA;

  // Ising Model
  // interaction radius
  s.Ising_neighboorhood = (int) RADIUS;
  // Temperature
  s.T = (double) TEMPERATURE;
  // Spin coupling
  s.J = -1 * (double) COUPLING;
  /* Set simulation flags */
  s.running = FALSE;
  s.initialized = FALSE;
  // Display rate to paint the lattice
  s.display_rate = (int) SAMPLE_RATE;
}


/* Activate function */
static void activate (GtkApplication *app, gpointer user_data)
  {
  initialize_simulation();
  // This function should only contain Gtk stuff
  /* General Gtk widgets for the Window packing */
  GtkWidget *window, *grid, *image_lattice, *label, *frame, *notebook, *box, *scale, *radio, *separator;
  GdkPixbuf *pixbuf;

  // Parameters Section
  /* Create a Gtk Notebook to hold pages of parameters */
  notebook = gtk_notebook_new ();
  // Pre-Pack all Radio buttons in a box to be framed
  // and then packed inyo the grid
  grid = gtk_grid_new();

  // Pre-pack: Coupling Neighborhood & Coupling Strengh Radio Buttons
  //
  // We make a Frame to back the Neighborhood size
  frame =  gtk_frame_new ("Neighborhood size");
  gtk_frame_set_label_align (GTK_FRAME(frame),0,0);
  // We make a box to hold stuff together
  box = gtk_box_new (GTK_ORIENTATION_VERTICAL, 0);
  // We make a radio button for the NN (r=1) choice
  radio = gtk_radio_button_new_with_label(NULL, "Nearest Neighboors (NN); r=1 ");
  g_signal_connect(GTK_TOGGLE_BUTTON(radio), "pressed", G_CALLBACK(on_radio_NN), (gpointer)"NN interaction selected");
  gtk_box_pack_start(GTK_BOX(box), radio, TRUE, TRUE, 0);
  // We make a radio button for the NNN (r=2) choice
  radio = gtk_radio_button_new_with_label_from_widget(GTK_RADIO_BUTTON(radio), "Next Nearest Neighboors (NNN); r=2 ");
  g_signal_connect(GTK_TOGGLE_BUTTON(radio), "pressed", G_CALLBACK(on_radio_NNN), (gpointer)"NNN interaction selected");
  gtk_box_pack_start(GTK_BOX(box), radio, TRUE, TRUE, 0);
  // Add the packing box to a Neighboorhood size Frame
  gtk_container_add (GTK_CONTAINER (frame), box);
  gtk_grid_attach (GTK_GRID (grid), frame, 0, 0, 1, 1);

  // Add a verical separator to the parameter grid for order
  separator = gtk_separator_new (GTK_ORIENTATION_VERTICAL);
  gtk_grid_attach (GTK_GRID (grid), separator, 1, 0, 1, 2);
  //
  // We make another Frame for the Coupling Strength
  //(ferro vs anti-ferro magnetic)
  frame =  gtk_frame_new ("Coupling strength");
  gtk_frame_set_label_align (GTK_FRAME(frame),1,0);
  // Again we make a box to hold stuff together
  box = gtk_box_new (GTK_ORIENTATION_VERTICAL, 0);
  // We make a radio button for the ferro magnetic case
  radio = gtk_radio_button_new_with_label(NULL, "  J = - 1 * kB");
  g_signal_connect(GTK_TOGGLE_BUTTON(radio), "pressed", G_CALLBACK(on_radio_ferro), (gpointer)"Ferro-Magnetic (J < 0) interaction selected");
  gtk_box_pack_start(GTK_BOX(box), radio, TRUE, TRUE, 0);
  // We make a radio button for the anti-ferro magnetic case
  radio = gtk_radio_button_new_with_label_from_widget(GTK_RADIO_BUTTON(radio), "  J = +1 * kB");
  g_signal_connect(GTK_TOGGLE_BUTTON(radio), "pressed", G_CALLBACK(on_radio_anti_ferro), (gpointer)"Anti-Ferro-Magnetic (J > 0) interaction selected");
  gtk_box_pack_start(GTK_BOX(box), radio, TRUE, TRUE, 0);
  // Add the packing box to the Coupling strength Frame
  gtk_container_add (GTK_CONTAINER (frame), box);
  gtk_grid_attach (GTK_GRID (grid), frame, 2, 0, 1, 2);

  // Now we are done with pre-packing the radio choices in Frames


  // Contact  Process (CP) Box
  // to control the parameters of the Contact Process
  // we mnake a box
  box = gtk_box_new (GTK_ORIENTATION_VERTICAL, 0);
  // BIRTH
  // scale bar to set birth rate
  scale = gtk_scale_new_with_range(GTK_ORIENTATION_HORIZONTAL,(gdouble) BETA_MIN, (gdouble) BETA_MAX, (gdouble) BETA_STEP);
  // we set it to its default value
  gtk_range_set_value (GTK_RANGE(scale), (gfloat) BETA);
  g_signal_connect (scale, "value-changed", G_CALLBACK (birth_rate_scale_moved), NULL);
  // we pack it in a Frame
  frame = gtk_frame_new("Birth rate");
  gtk_container_add (GTK_CONTAINER (frame), scale);
  // we add that Frame to the CP box
  gtk_container_add (GTK_CONTAINER (box), frame);
  // DEATH
  // scale bar to set death rate
  scale = gtk_scale_new_with_range(GTK_ORIENTATION_HORIZONTAL, (gdouble) DELTA_MIN, (gdouble) DELTA_MAX, (gdouble) DELTA_STEP);
  // we set it to its default value
  gtk_range_set_value (GTK_RANGE(scale), (gfloat) DELTA);
  g_signal_connect (scale, "value-changed", G_CALLBACK (death_rate_scale_moved), NULL);
  // we pack it in a Frame
  frame = gtk_frame_new ("Death rate");
  gtk_container_add (GTK_CONTAINER (frame), scale);
  // we add that Frame to the CP box
  gtk_container_add (GTK_CONTAINER (box), frame);
  // ALPHA
  // scale bar to set the cell differenciation rate
  scale = gtk_scale_new_with_range(GTK_ORIENTATION_HORIZONTAL,(gdouble) ALPHA_MIN, (gdouble) ALPHA_MAX, (gdouble) ALPHA_STEP);
  // we set it to its default value
  gtk_range_set_value (GTK_RANGE(scale), (gfloat) ALPHA);
  g_signal_connect (scale, "value-changed", G_CALLBACK (differenciation_rate_scale_moved), NULL);
  // we pack it in a Frame
  frame = gtk_frame_new ("Differenciation rate");
  gtk_container_add (GTK_CONTAINER (frame), scale );
  // we add that Frame to the CP box
  gtk_container_add (GTK_CONTAINER (box), frame);

  /* Make a Contact Process label and put it with its box in the Notebook*/
  label = gtk_label_new ("Contact Process");
  gtk_notebook_append_page (GTK_NOTEBOOK (notebook), box, label);


  // Ising Model (IM) Box
  // to control the parameter of the Ising model
  box = gtk_box_new (GTK_ORIENTATION_VERTICAL, 0);
  // TEMPERATURE
  // make a scale bar to set Temperature
  scale = gtk_scale_new_with_range(GTK_ORIENTATION_HORIZONTAL, (gdouble) TEMPERATURE_MIN, (gdouble) TEMPERATURE_MAX, (gdouble) TEMPERATURE_STEP);
  // we set it to its default value
  gtk_range_set_value (GTK_RANGE(scale), (gfloat) TEMPERATURE);
  g_signal_connect (scale, "value-changed", G_CALLBACK (temperature_scale_moved), NULL);
  frame =  gtk_frame_new ("Temperature");
  gtk_container_add (GTK_CONTAINER (frame), scale);
  gtk_container_add (GTK_CONTAINER (box), frame);

  // Add a verical separator to the parameter grid for order
  separator = gtk_separator_new (GTK_ORIENTATION_HORIZONTAL);
  gtk_container_add (GTK_CONTAINER (box), separator);

  // add the grid contasining the pre-packed radio buttons to control spin coupling
  frame =  gtk_frame_new ("Spin interaction");
  gtk_container_add (GTK_CONTAINER (frame), grid);
  gtk_container_add (GTK_CONTAINER (box), frame);

  /* make a IM label and put it with its box in the Notebook*/
  label = gtk_label_new ("Ising Model");
  gtk_notebook_append_page (GTK_NOTEBOOK (notebook), box, label);



  // Initial Conditions (IC) Box
  // make a box to hold a bunch of  radio buttons reppresenting different initial configurations
  box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
  // make buttons and pack them up
  // -1-
  radio = gtk_radio_button_new_with_label(NULL, "Single differenciated site (spin up or down)");
  g_signal_connect(GTK_TOGGLE_BUTTON(radio), "pressed", G_CALLBACK(on_radio_initial_condition_1), (gpointer)"option 1 selected");
  gtk_box_pack_start(GTK_BOX(box), radio, TRUE, TRUE, 0);
  // -2-
  radio = gtk_radio_button_new_with_label_from_widget(GTK_RADIO_BUTTON(radio), "Single un-differenciated site");
  g_signal_connect(GTK_TOGGLE_BUTTON(radio), "pressed", G_CALLBACK(on_radio_initial_condition_2), (gpointer)"option 2 selected");
  gtk_box_pack_start(GTK_BOX(box), radio, TRUE, TRUE, 0);
  // -3-
  radio = gtk_radio_button_new_with_label_from_widget(GTK_RADIO_BUTTON(radio), "Single un-differenciated cluster");
  g_signal_connect(GTK_TOGGLE_BUTTON(radio), "pressed", G_CALLBACK(on_radio_initial_condition_3) , (gpointer)"option 3 selected");
  gtk_box_pack_start(GTK_BOX(box), radio, TRUE, TRUE, 0);
  // -4-
  radio = gtk_radio_button_new_with_label_from_widget(GTK_RADIO_BUTTON(radio), "Single differenciated cluster (spin up or down)");
  g_signal_connect(GTK_TOGGLE_BUTTON(radio), "pressed", G_CALLBACK(on_radio_initial_condition_4), (gpointer)"option 4 selected");
  gtk_box_pack_start(GTK_BOX(box), radio, TRUE, TRUE, 0);
  // -5-
  radio = gtk_radio_button_new_with_label_from_widget(GTK_RADIO_BUTTON(radio), "Fully occupied un-differenciated lattice");
  g_signal_connect(GTK_TOGGLE_BUTTON(radio), "pressed", G_CALLBACK(on_radio_initial_condition_5), (gpointer)"option 4 selected");
  gtk_box_pack_start(GTK_BOX(box), radio, TRUE, TRUE, 0);
  // make IC label for Initial Conditions page and put it in the Notebook
  label = gtk_label_new ("Init Lattice");
  frame =  gtk_frame_new ("Configuration at genesis  (t = 0)");
  gtk_frame_set_label_align (GTK_FRAME(frame),0,1);
  /* --- Add the page with the frame and label --- */
  gtk_notebook_append_page (GTK_NOTEBOOK (notebook), frame, label);
  // Add IC box to frame to be put on the frame of the third page of the Notebook
  gtk_container_add (GTK_CONTAINER (frame), box);


  // Parameters section final touch
  /* We make a final frame for the whole Notebook */
  GtkWidget *parameters_frame;
  parameters_frame =  gtk_frame_new ("Parameters");
  gtk_frame_set_label_align (GTK_FRAME(parameters_frame),0,1);

  /*  add the notebook to the recently made frame*/
  gtk_container_add (GTK_CONTAINER (parameters_frame), notebook);



  // MAIN WINDOW
  /* Create a new window, and set its title */
  window = gtk_application_window_new (app);
  gtk_window_set_title (GTK_WINDOW (window), "Contact Process Ising Model");
  gtk_window_set_resizable (GTK_WINDOW (window), FALSE);

  /* Use Gtk grid to pack our widgets in the Main App Window */
  grid = gtk_grid_new ();


  /* Add the frame to the Application window's main grid*/
  // Pack it into the window main grid
  gtk_grid_attach (GTK_GRID (grid), parameters_frame, 0, 0, 5, 3);



  // PIX BUFFER
  /* Pixel buffer @ start up and default canvas display */
  pixbuf = gdk_pixbuf_new (GDK_COLORSPACE_RGB, 0, 8, X_SIZE, Y_SIZE);
  image_lattice = gtk_image_new_from_pixbuf (pixbuf);
  paint_a_background (image_lattice);
  // We place the image on row 7 of our grid spanning 5 columns
  gtk_grid_attach (GTK_GRID (grid), image_lattice, 0, 7, 5, 1); 
  /* Position (0,3) spanning 5 col and 1 row */


  // Simulation CONTROLS
  GtkWidget *button, *ctrl_frame, *button_box;

  button_box = gtk_button_box_new (GTK_ORIENTATION_HORIZONTAL);
  ctrl_frame =  gtk_frame_new ("Simulation Control");
  // Initialize Lattice
  button = gtk_button_new_with_label ("Init");
  g_signal_connect (button, "clicked", G_CALLBACK (init_lattice), GTK_IMAGE (image_lattice));
  gtk_container_add (GTK_CONTAINER (button_box), button);
  // Start
  button = gtk_button_new_with_label ("Start");
  g_signal_connect (button, "clicked", G_CALLBACK (on_button_start_simulation), GTK_IMAGE (image_lattice));
  gtk_container_add (GTK_CONTAINER (button_box), button);
  // Stop
  button = gtk_button_new_with_label ("Stop");
  g_signal_connect (button, "clicked", G_CALLBACK (on_button_stop_simulation), NULL);
  gtk_container_add (GTK_CONTAINER (button_box), button);
  // About (to be fair to the user)
  button = gtk_button_new_with_label ("About");
  g_signal_connect (button, "clicked", G_CALLBACK (on_button_show_about), NULL);
  gtk_container_add (GTK_CONTAINER (button_box), button);
  // Quit
  button = gtk_button_new_with_label ("Quit");
  g_signal_connect (button, "clicked", G_CALLBACK (on_button_stop_simulation), NULL);
  g_signal_connect_swapped (button, "clicked", G_CALLBACK (gtk_widget_destroy), window);
  gtk_container_add (GTK_CONTAINER (button_box), button);

  // We now add all the buttons in their specialized to a frame
  gtk_container_add (GTK_CONTAINER (ctrl_frame), button_box);
  // We finally place the frame on row 8 of our grid spanning 5 columns
  gtk_grid_attach (GTK_GRID (grid), ctrl_frame, 0, 8, 5, 1);


  // FRAME_SAMPLE_RATE
  // scale bar to set display rate
  scale = gtk_scale_new_with_range(GTK_ORIENTATION_HORIZONTAL,1,100,10);
  gtk_range_set_value (GTK_RANGE(scale), (gfloat) SAMPLE_RATE);
  g_signal_connect (scale, "value-changed", G_CALLBACK (display_rate_scale_moved), NULL);
  // we pack it in a Frame
  frame = gtk_frame_new ("Display rate");
  gtk_container_add (GTK_CONTAINER (frame), scale );
  // we add that Frame to the CP box
  // gtk_container_add (GTK_CONTAINER (box), frame);
  gtk_grid_attach (GTK_GRID (grid), frame, 0, 9, 5, 1);


  /* Pack the main grid into the window */
  gtk_container_add (GTK_CONTAINER (window), grid);
  /* Show the window and all widgets in it */
  gtk_widget_show_all (window);
}




/* Implementation of put pixel function. Code retrieved from:
   https://developer.gnome.org/gdk-pixbuf/stable/gdk-pixbuf-The-GdkPixbuf-Structure.html */
void put_pixel (GdkPixbuf *pixbuf, int x, int y,
                guchar red, guchar green, guchar blue, guchar alpha)
  {
  guchar *pixels, *p;
  int rowstride, numchannels;
  numchannels = gdk_pixbuf_get_n_channels(pixbuf);
  rowstride = gdk_pixbuf_get_rowstride(pixbuf);
  pixels = gdk_pixbuf_get_pixels(pixbuf);
  p = pixels + y * rowstride + x * numchannels;
  p[0] = red;	p[1] = green; p[2] = blue; p[3] = alpha;
  return;
  }



/* Creates a pixel buffer and paints an image to display as default canvas */
static void paint_a_background (gpointer data)
  {
  GdkPixbuf *p;
  p = gdk_pixbuf_new (GDK_COLORSPACE_RGB, 0, 8, X_SIZE, Y_SIZE);
  /* Paint a background canvas for start up image */
  int x, y;
  for (x = 0; x < X_SIZE; x++)
    	{
        for (y = 0; y < Y_SIZE; y++)
            {
            put_pixel (p, (int) x, (int) y,
                       (guchar) x, (guchar) y, (guchar) x, 255);
            }
      }
  gtk_image_set_from_pixbuf (GTK_IMAGE (data), GDK_PIXBUF (p));
  g_object_unref (p);
  }


/* Function that paints the pixel buffer with the simulation data   */
// The states are:
//  0: vacant
// -1: spin down
// +1: spin up
//  2: undifferenciated
static void paint_lattice (gpointer data)
  {
  // we make a Gdk pixbuffer to paint configurations
  GdkPixbuf *p;
  p = gdk_pixbuf_new (GDK_COLORSPACE_RGB, 0, 8, X_SIZE, Y_SIZE);
  /* Paint lattice configuration to a pixel buffer */
  int x, y;
  for (x = 0; x < X_SIZE; x++)
    {
    for (y = 0; y < Y_SIZE; y++)
      {
      switch (s.lattice_configuration[x][y])
        {
        case 0:	/* Empty (vacant) site  (black) */
          put_pixel (p, (int) x, (int) y, 
                     (guchar) 0, (guchar) 0, (guchar) 0, 255);
          break;
        case -1:	/* Spin down (occupied) site (green) */
          put_pixel (p, (int) x, (int) y, 
                       (guchar) 0, (guchar) 255, (guchar) 0, 255);
          break;
        case 1:	  /* Spin  up  (occupied) site (magenta) */
          put_pixel (p, (int) x, (int) y, 
                     (guchar) 255, (guchar) 0, (guchar) 255, 255);
          break;
        case 2:	/*  Un-differentiated (occupied) site (white) */
          put_pixel (p, (int) x, (int) y, 
                     (guchar) 255, (guchar) 255, (guchar) 255, 255);
          break;
        }
      }
    }
  gtk_image_set_from_pixbuf (GTK_IMAGE (data), GDK_PIXBUF (p));
  g_object_unref (p);
  }


/* Main function spanning a Gtk Application object */
int main (int argc, char **argv)
  {
  GtkApplication *app;
  int status;
  app = gtk_application_new ("keymer.lab.contact_process_ising_model", G_APPLICATION_FLAGS_NONE);
  g_signal_connect (app, "activate", G_CALLBACK (activate), NULL);
  status = g_application_run (G_APPLICATION (app), argc, argv);
  g_object_unref (app);
  return status;
  }

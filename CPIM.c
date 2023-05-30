/* Should write a description of what this code does */

/* Libraries */
//#include <windows.h>
#include <stdlib.h>
#include <gtk/gtk.h> /* GUI, Gtk library */
#include "mt64.h" /* Pseudo-random number generation MT library (64 bit) */
#include <math.h> /* Math to transform random n from continuous to discrete */
#include <time.h> /* Used to seed pseudo-random number generator */
#include <stdio.h>


/* Constants */
#define X_SIZE 200
#define Y_SIZE 200


/* Structure with the simulation data
   i.e. a class-like object that allow us to keep track of data */
struct simulation
  {
  int lattice_configuration[X_SIZE][Y_SIZE]; /* Store latice configuration */
  gint run;                   /* Time handler tag */
  gboolean running;           /* Are we running? */
  gboolean initialized;       /* Have we been initialized? */
  int generation_time;        /* Generations simulated */
  int influence_radius;       /* Distance (in number of sites) of Ising-like 
                                 interactions */
  double occupancy;           /* Lattice occupancy */
  double up;                  /* Keeps track of spins in the up (+1) state */
  double down;                /* Keeps track of spins in the down (-1) state */
  double birth_rate;          /* Contact Process' birth rate */
  double death_rate;          /* Contact Process' death rate */
  double differentiation_rate;/* Rate at which occupied sites get a 
                                 differentiated state */
  double temperature;         /* Ising model's temperature (T) */
  double coupling;            /* Ising model's coupling (J) parameter */
  double magnetic_field;      /* Ising model's magnetic field (B) */
  } s;                        /* Instance s of struct */


static void stop_simulation (gpointer data)
  {
  if (s.running)
    {
    g_source_remove (s.run);
    s.running = FALSE;
    g_print ("Simulation stopped\n");
    }
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


/* Function that creates a new pixel buffer and paints an image to display as
   default canvas */
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


/* Function that paints the pixel buffer with the simulation data
   */
static void paint_lattice (gpointer data)
  {
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
        case 0:	/* Empty site. Paint white */
          put_pixel (p, (int) x, (int) y, 
                     (guchar) 255, (guchar) 255, (guchar) 255, 255);
          break;
        // case 2:	/* Occupied, un-differentiated site. Paint black */
        //   put_pixel (p, (int) x, (int) y, 
        //              (guchar) 48, (guchar) 48, (guchar) 48, 255);
        //   break;
        case 1:	/* Occupied, +1 site. Paint red */
          put_pixel (p, (int) x, (int) y, 
                     (guchar) 249, (guchar) 237, (guchar) 105, 255);
          break;
        case -1:	/* Occupied, -1 site. Paint green */
          put_pixel (p, (int) x, (int) y, 
                     (guchar) 106, (guchar) 44, (guchar) 112, 255);
          break;
        }
      }
    }
  gtk_image_set_from_pixbuf (GTK_IMAGE (data), GDK_PIXBUF (p));
  g_object_unref (p);
  }


/* Function to get the closest neighbors (separated by r sites) of a given site 
   in the lattice */
void get_closest_neighbors(int x, int y, int r, int* neighbors)
  {
  int n = 0;
  for (int i = -r; i <= r; i++) 
    {
    for (int j = -r; j <= r; j++) 
      {
      if (i == 0 && j == 0) {continue;}
			else if (abs(i) + abs(j) > r) {continue;}
      neighbors[n++] = s.lattice_configuration[(int) (X_SIZE + x + i) % X_SIZE]
                                              [(int) (Y_SIZE + y + j) % Y_SIZE];
      }
    }
  }


/* Function used to compute the energy value of the site located at (x, y) */
double compute_energy (int x, int y)
  {
  double energy;
  int spin;
  double neighborhood_configuration = 0;
  int num_neighbors = pow(s.influence_radius, 2) + pow(s.influence_radius + 1, 2);
  int neighborhood[num_neighbors];

  spin = s.lattice_configuration[x][y];

  get_closest_neighbors (x, y, s.influence_radius, neighborhood);
  
  for (int n = 0; n < num_neighbors; n++) 
    {
    // /* If neighbor has no spin, go to the next one */
    // if (neighborhood[n] == 2) {continue;} 
    neighborhood_configuration += neighborhood[n];
    }
  energy = spin*(s.coupling * neighborhood_configuration - s.magnetic_field);
  return energy;
  }




/* Update function which simulates the stochastic process and updates the
   configuration of the lattice */
int update_lattice (gpointer data)
  {
  // int random_neighbor;
  int random_neighbor_state;
  // int random_spin;
  double spin_energy, spin_energy_diff;
  double transition_probability;
  int random_x_coor, random_y_coor;

  int contact_process_radius = 1;
  int num_neighbors = pow(contact_process_radius, 2) + pow(contact_process_radius + 1, 2);
  int neighbors[num_neighbors];

  for (int site = 0; site < (int) (Y_SIZE*X_SIZE); site++)
    {
    /* Pick a random focal site */
    random_x_coor = (int) floor (genrand64_real3 ()* X_SIZE);
    random_y_coor = (int) floor (genrand64_real3 ()* Y_SIZE);
    switch (s.lattice_configuration[random_x_coor][random_y_coor])
      {
      case 0: /* Site is empty */
        /* Chose a random neighbor from the num_neighbors posible ones */
        get_closest_neighbors (random_x_coor, random_y_coor, 
                               contact_process_radius, neighbors);
        random_neighbor_state = neighbors[(int) floor (genrand64_real3 () 
                                          * (num_neighbors))];

        /* If its random neighbor is occupied: put a copy at the focal site 
           with probability brith_rate * dt */
        if (genrand64_real2 () < s.birth_rate)
          {
          if (random_neighbor_state == 2)
            {
            s.lattice_configuration[random_x_coor][random_y_coor] = 2;
            s.occupancy ++; 
            }
          else if (random_neighbor_state == 1)
            {
            s.lattice_configuration[random_x_coor][random_y_coor] = 1; 
            s.occupancy ++; 
            s.up ++;
            }
          else if (random_neighbor_state == -1)
            {
            s.lattice_configuration[random_x_coor][random_y_coor] = -1;
            s.occupancy ++;
            s.down ++;
            }
          }
        break; /* break case 0 */
      // case 2: /* Focal point is in the occupied, undifferentiated state */
      //   if (genrand64_real2 () < s.death_rate)
      //     {
      //     s.lattice_configuration[random_x_coor][random_y_coor] = 0;
      //     s.occupancy --;
      //     }
      //   else if (genrand64_real2 () < s.differentiation_rate)
      //     {
      //     /* Set an occupied site in the middle of the lattice */
      //     random_spin = ((genrand64_int64 () % 2) * 2) - 2;
      //     if (random_spin == 1)
      //       {
      //       s.lattice_configuration[random_x_coor][random_y_coor] = random_spin;
      //       s.up ++;
      //       }
      //     else if (random_spin == -1)
      //       {
      //       s.lattice_configuration[random_x_coor][random_y_coor] = random_spin;
      //       s.down ++;
      //       }
      //     }
      //   break;
      case 1: /* Focal point is in the up (+1) state */
        if (genrand64_real2 () < s.death_rate)
          {
          s.lattice_configuration[random_x_coor][random_y_coor] = 0;
          s.occupancy --;
          s.up --;
          }
        else
          {
          spin_energy = compute_energy (random_x_coor, random_y_coor);
          spin_energy_diff = -2 * spin_energy;
          transition_probability = exp (-spin_energy_diff/s.temperature);
          if (spin_energy_diff < 0 || 
              genrand64_real2 () < transition_probability)
            {
            s.lattice_configuration[random_x_coor][random_y_coor] = -1;
            s.up --;
            s.down ++;
            }
          }
        break;
      case -1: /* Focal point is in the down (-1) state */
        if (genrand64_real2 () < s.death_rate)
          {
          s.lattice_configuration[random_x_coor][random_y_coor]= 0;
          s.occupancy --;
          s.down --;
          }
        else
          {
          spin_energy = compute_energy (random_x_coor, random_y_coor);
          spin_energy_diff = -2 * spin_energy;
          transition_probability = exp (-spin_energy_diff/s.temperature);
          if (spin_energy_diff < 0 || 
              genrand64_real2 () < transition_probability)
            {
            s.lattice_configuration[random_x_coor][random_y_coor] = 1;
            s.down --;
            s.up ++;
            }
          }
        break;
      }
    }
  s.generation_time ++;
  paint_lattice (data);
  g_print ("Gen: %d \t Occupancy: %f \t Up: %f \t Down: %f\n", 
           s.generation_time, s.occupancy/(X_SIZE*Y_SIZE), 
           s.up/s.occupancy, s.down/s.occupancy);
  /* This is a simple occupancy check to avoid keep running the simulation
    when there's no particle left on the lattice */
  if (s.occupancy <= 0) {stop_simulation (data);}
  return 0;
  }


/* Time handler that connects the update function to the gtk loop for its
   scomputation */
gboolean time_handler (gpointer data)
  {
  update_lattice (data);
  return TRUE;
  }


/* Callback to initialize the lattice button click */
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
  s.occupancy = 0.0;
  s.up = 0;
  s.down = 0;

  // /* Set a site in the moddle of the lattice occupied (un-differntiated)  */
  // s.lattice_configuration[(int) X_SIZE/2][(int) Y_SIZE/2] = 2;



  /* Set an occupied site in the middle of the lattice */
  random_spin = (int) ((genrand64_int64 () % 2) * 2) - 1;
  if (random_spin == 1)
    {
    s.lattice_configuration[(int) X_SIZE/2][(int) Y_SIZE/2] = random_spin;
    s.up ++;
    }
  else if (random_spin == -1)
    {
    s.lattice_configuration[(int) X_SIZE/2][(int) Y_SIZE/2] = random_spin;
    s.down ++;
    }
  s.occupancy ++;

  s.initialized = TRUE;
  s.generation_time = 0;
  paint_lattice (data);
  g_print ("Lattice initialized\n");
  }



// HERE GOES THE ABOUT DIALOG BOX For info at a website: lab wiki on the contact process
static void show_about(GtkWidget *widget, gpointer data)
        {
        GdkPixbuf *pixbuf = gdk_pixbuf_new_from_file("kimero_LAB_transparent.tiff", NULL);
        GtkWidget *dialog = gtk_about_dialog_new();
        gtk_about_dialog_set_program_name (GTK_ABOUT_DIALOG(dialog),
                                    "Contact Process Ising Model App");
        gtk_about_dialog_set_version(GTK_ABOUT_DIALOG(dialog), "version 0.1, 2023");
        gtk_about_dialog_set_copyright(GTK_ABOUT_DIALOG(dialog),"Open Source Code");
        gtk_about_dialog_set_comments(GTK_ABOUT_DIALOG(dialog),
     "The Contact process Ising model.");
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


/*  Callback to respond Gtk scale slide move event */
static void influence_radius_scale_moved (GtkRange *range, gpointer user_data)
  {
  GtkWidget *label = user_data;
  gdouble pos = gtk_range_get_value (range);
  s.influence_radius = (float) pos;
  gchar *str = g_strdup_printf ("radius = %d", (int) pos);
  gtk_label_set_text (GTK_LABEL (label), str);
  g_free (str);
  }


/*  Callback to respond Gtk scale slide move event */
static void birth_rate_scale_moved (GtkRange *range, gpointer user_data)
  {
  GtkWidget *label = user_data;
  gdouble pos = gtk_range_get_value (range);
  s.birth_rate = (float) pos;
  gchar *str = g_strdup_printf ("beta = %.2f", pos);
  gtk_label_set_text (GTK_LABEL (label), str);
  g_free (str);
  }


/*  Callback to respond Gtk scale slide move event */
static void death_rate_scale_moved (GtkRange *range, gpointer user_data)
  {
  GtkWidget *label = user_data;
  gdouble pos = gtk_range_get_value (range);
  s.death_rate = (float) pos;
  gchar *str = g_strdup_printf ("delta = %.2f", pos);
  gtk_label_set_text (GTK_LABEL (label), str);
  g_free (str);
  }


/*  Callback to respond Gtk scale slide move event */
static void temperature_scale_moved (GtkRange *range, gpointer user_data)
  {
  GtkWidget *label = user_data;
  gdouble pos = gtk_range_get_value (range);
  s.temperature = (float) pos;
  gchar *str = g_strdup_printf ("temperature = %.2f", pos);
  gtk_label_set_text (GTK_LABEL (label), str);
  g_free (str);
  }


/*  Callback to respond Gtk scale slide move event */
static void coupling_scale_moved (GtkRange *range, gpointer user_data)
  {
  GtkWidget *label = user_data;
  gdouble pos = gtk_range_get_value (range);
  s.coupling = (float) pos;
  gchar *str = g_strdup_printf ("coupling = %.2f", pos);
  gtk_label_set_text (GTK_LABEL (label), str);
  g_free (str);
  }


/*  Callback to respond Gtk scale slide move event */
static void magnetic_field_scale_moved (GtkRange *range, gpointer user_data)
  {
  GtkWidget *label = user_data;
  gdouble pos = gtk_range_get_value (range);
  s.magnetic_field = (float) pos;
  gchar *str = g_strdup_printf ("magnetic field = %.2f", pos);
  gtk_label_set_text (GTK_LABEL (label), str);
  g_free (str);
  }





/* Activate function with all widget creation and initialization. Ideally this
should be handled by a GtkBuilder object */
static void activate (GtkApplication *app, gpointer user_data)
  {
  /* Declare a bunch of Gtk widgets for the GUI */
  GtkWidget *window, *grid, *button, *image_lattice;
  GdkPixbuf *pixbuf;

  /* To control parameters of the process */
  GtkWidget *influence_radius_scale, *influence_radius_label;
  GtkWidget *birth_rate_scale, *birth_rate_label;
  GtkWidget *death_rate_scale, *death_rate_label;
  GtkWidget *temperature_scale, *temperature_label;
  GtkWidget *coupling_scale, *coupling_label;
  GtkWidget *magnetic_field_scale, *magnetic_field_label;

  /* Initialize Mersenne Twister algorithm for random number genration */
  unsigned int seed = (unsigned int) time (NULL);
  init_genrand64 (seed);

  /* Set default parameters of the simulation */
  s.influence_radius = 1;
  s.birth_rate = 1.00;
  s.death_rate = 0.00;
  // s.differentiation_rate = 1.0;
  s.temperature = 0.001;
  s.coupling = -1.00;
  s.magnetic_field = 0.00;

  /* Set simulation flags */
  s.running = FALSE;
  s.initialized = FALSE;

  /* Create a new window, and set its title */
  window = gtk_application_window_new (app);
  gtk_window_set_title (GTK_WINDOW (window), "Contact Process Ising Model");
  gtk_window_set_resizable (GTK_WINDOW (window), FALSE);

  /* Here we make a grid that is going pack our widgets */
  grid = gtk_grid_new ();
  /* Pack the grid into the window */
  gtk_container_add (GTK_CONTAINER (window), grid);

  /* Vecinity scale slide bar */
  influence_radius_scale = gtk_scale_new_with_range (GTK_ORIENTATION_HORIZONTAL, 1, 3, 1);
  influence_radius_label = gtk_label_new ("radius");
  g_signal_connect (influence_radius_scale, "value-changed", G_CALLBACK (influence_radius_scale_moved),  influence_radius_label);
  gtk_grid_attach (GTK_GRID (grid), influence_radius_scale, 0, 0, 1, 1); /* Position (0,0) spanning 1 col and 1 row */
  gtk_grid_attach (GTK_GRID (grid),  influence_radius_label, 1, 0, 1, 1); /* Position (1,0) spanning 1 col and 1 row */

  /* Birth rate scale slide bar */
  birth_rate_scale = gtk_scale_new_with_range (GTK_ORIENTATION_HORIZONTAL, 0.0, 1.0, 0.01);
  birth_rate_label = gtk_label_new ("beta");
  gtk_range_set_value (GTK_RANGE (birth_rate_scale), 1.0);
  g_signal_connect (birth_rate_scale, "value-changed", G_CALLBACK (birth_rate_scale_moved), birth_rate_label);
  gtk_grid_attach (GTK_GRID (grid), birth_rate_scale, 0, 1, 1, 1); /* Position (0,1) spanning 1 col and 1 row */
  gtk_grid_attach (GTK_GRID (grid), birth_rate_label, 1, 1, 1, 1); /* Position (1,1) spanning 1 col and 1 row */

  /* Death rate scale slide bar */
  death_rate_scale = gtk_scale_new_with_range (GTK_ORIENTATION_HORIZONTAL, 0.0, 1.0, 0.01);
  death_rate_label = gtk_label_new ("delta"); //LABEL to be shown J
  g_signal_connect (death_rate_scale, "value-changed", G_CALLBACK (death_rate_scale_moved), death_rate_label);
  gtk_grid_attach (GTK_GRID (grid), death_rate_scale, 2, 1, 1, 1); /* Position (2,1) spanning 1 col and 1 row */
  gtk_grid_attach (GTK_GRID (grid), death_rate_label, 3, 1, 1, 1); /* Position (3,1) spanning 1 col and 1 row */

  /* Temperature (T) scale slide bar */
  temperature_scale = gtk_scale_new_with_range (GTK_ORIENTATION_HORIZONTAL, 0.001, 50, 0.001);
  temperature_label = gtk_label_new ("temperature"); //LABEL to be shown T
  g_signal_connect (temperature_scale, "value-changed", G_CALLBACK (temperature_scale_moved), temperature_label);
  gtk_grid_attach (GTK_GRID (grid), temperature_scale, 0, 2, 1, 1); /* Position (0,2) spanning 1 col and 1 row */
  gtk_grid_attach (GTK_GRID (grid), temperature_label, 1, 2, 1, 1); /* Position (1,2) spanning 1 col and 1 row */

  /* Coupling (J) scale slide bar */
  coupling_scale = gtk_scale_new_with_range (GTK_ORIENTATION_HORIZONTAL, -1.0, 1.0, 0.01);
  coupling_label = gtk_label_new ("coupling"); //LABEL to be shown J
  g_signal_connect (coupling_scale,"value-changed", G_CALLBACK (coupling_scale_moved), coupling_label);
  gtk_grid_attach (GTK_GRID (grid), coupling_scale, 2, 2, 1, 1); /* Position (2,2) spanning 1 col and 1 row */
  gtk_grid_attach (GTK_GRID (grid), coupling_label, 3, 2, 1, 1); /* Position (3,2) spanning 1 col and 1 row */

  /* Magnetic field (B) scale slide bar */
  magnetic_field_scale = gtk_scale_new_with_range (GTK_ORIENTATION_HORIZONTAL, -1.0, 1.0, 0.01);
  magnetic_field_label = gtk_label_new ("magnetic field"); //LABEL to be shown J
  gtk_range_set_value (GTK_RANGE (magnetic_field_scale), 0.0);
  g_signal_connect (magnetic_field_scale, "value-changed", G_CALLBACK (magnetic_field_scale_moved), magnetic_field_label);
  gtk_grid_attach (GTK_GRID (grid), magnetic_field_scale, 0, 3, 1, 1); /* Position (2,3) spanning 1 col and 1 row */
  gtk_grid_attach (GTK_GRID (grid), magnetic_field_label, 1, 3, 1, 1); /* Position (3,3) spanning 1 col and 1 row */
  

  /* Pixel buffer @ start up and default canvas display */
  pixbuf = gdk_pixbuf_new (GDK_COLORSPACE_RGB, 0, 8, X_SIZE, Y_SIZE);
  image_lattice = gtk_image_new_from_pixbuf (pixbuf);
  paint_a_background (image_lattice);
  gtk_grid_attach (GTK_GRID (grid), image_lattice, 0, 4, 5, 1); /* Position (0,3) spanning 5 col and 1 row */

  /* ----------------------------  INIT BUTTON  ----------------------------- */
  button = gtk_button_new_with_label ("Initialize");
  g_signal_connect (button, "clicked", G_CALLBACK (init_lattice), GTK_IMAGE (image_lattice));
  gtk_grid_attach (GTK_GRID (grid), button, 0, 5, 1, 1); /* Position (0,4) spanning 1 col and 1 row */
  /* ----------------------------  START BUTTON  ---------------------------- */
  button = gtk_button_new_with_label ("Start");
  g_signal_connect (button, "clicked", G_CALLBACK (on_button_start_simulation), GTK_IMAGE (image_lattice));
  gtk_grid_attach (GTK_GRID(grid), button, 1, 5, 1, 1); /* Position (1,4) spanning 1 col and 1 row */
  /* ----------------------------  STOP BUTTON  ----------------------------- */
  button = gtk_button_new_with_label ("Stop");
  g_signal_connect (button, "clicked", G_CALLBACK (on_button_stop_simulation), NULL);
  gtk_grid_attach (GTK_GRID(grid), button, 2, 5, 1, 1); /* Position (2,4) spanning 1 col and 1 row */

// -----   ABOUT ? BUTTON  ----
	button = gtk_button_new_with_label ("?");
	g_signal_connect (button, "clicked", G_CALLBACK(show_about), GTK_WINDOW(window));
	gtk_grid_attach (GTK_GRID (grid), button, 3, 5, 1, 1); // position (3,3) spanning 1 col and 1 raw

  /* ----------------------------  QUIT BUTTON  ----------------------------- */
  button = gtk_button_new_with_label ("Quit");
  g_signal_connect_swapped (button, "clicked", G_CALLBACK (gtk_widget_destroy), window);
  gtk_grid_attach (GTK_GRID (grid), button, 4, 5, 1, 1); /* Position (4,4) spanning 1 col and 1 row */

  // Show the window and all widgets
  gtk_widget_show_all (window);
  }


/* main function */
int main (int argc, char **argv)
  {
  GtkApplication *app;
  int status;
  app = gtk_application_new ("keymer.lab.contact_process", G_APPLICATION_FLAGS_NONE);
  g_signal_connect (app, "activate", G_CALLBACK (activate), NULL);
  status = g_application_run (G_APPLICATION (app), argc, argv);
  g_object_unref (app);
  return status;
  }

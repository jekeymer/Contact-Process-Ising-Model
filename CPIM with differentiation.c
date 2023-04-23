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
#define X_SIZE 250
#define Y_SIZE 250

/* Structure with the simulation data
   i.e. a class-like object that allow us to keep track of data */
struct simulation
  {
  int lattice_configuration[X_SIZE][Y_SIZE]; /* Store latice configuration */
  gint run;                   /* Time handler tag */
  gboolean running;           /* Are we running? */
  gboolean initialized;       /* Have we been initialized? */
  int generation_time;        /* Generations simulated */
  int influence_radius;
  double occupancy;           /* Lattice occupancy */
  double up;
  double down;
  double birth_rate;
  double death_rate;
  double differentiation_rate;
  double temperature;         /* Ising model's temperature (T) */
  double coupling;            /* Ising model's coupling (J) parameter */
  } s;                        /* Instance s of struct */


/* Declare put_pixel function to access individual pixel data on a pixel buffer.
   Implemented at the end of document. */
void put_pixel (GdkPixbuf *pixbuf, int x, int y,
                guchar red, guchar green, guchar blue, guchar alpha);


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


/* Function that keeps track of each site's configuration
   */
static void count (gpointer data)
  {
  GdkPixbuf *p;
  p = gdk_pixbuf_new (GDK_COLORSPACE_RGB, 0, 8, X_SIZE, Y_SIZE);
  int x, y;
  int RP = 0;
  int GP = 0;
  int BP = 0;
  int WP = 0;
  for (x = 0; x < X_SIZE; x++)
  {
  for (y = 0; y < Y_SIZE; y++)
  {
  switch (s.lattice_configuration[x][y])
  {
  case 0:	/* Empty (white )*/
    				WP++;
  put_pixel (p, (int) x, (int) y, (guchar) 255, (guchar) 255, (guchar) 255, 255);
  break;
  case 1:	/* Used*/
    				BP++;
  put_pixel (p, (int) x, (int) y, (guchar) 48, (guchar) 48, (guchar) 48, 255);
  break;
  case 2:	//+1: Rojo
    				RP++;
  put_pixel (p, (int) x, (int) y, (guchar) 255, (guchar) 227, (guchar) 138, 255); 
  break;
  case 3:	//-1: Verde
    				GP++;
  put_pixel (p, (int) x, (int) y, (guchar) 149, (guchar) 255, (guchar) 211, 255);
  break;
  default:
  put_pixel (p, (int) x, (int) y, (guchar) 0, (guchar) 0, (guchar) 0, 255);
  break;
  }
  }
  }
  FILE *fp;
  fp = fopen ("Ising.txt", "a+");		//Crea/Abre el archivo y escribe en él
  //fprintf(fp, "Gen:\tOccupancy:\tCoup:\tJ:\tWithe:\tBlack:\tRed:\tGreen:\tX*Y:\tW+B+R+G:\n");
  //fprintf(fp, "%u\t%f\t%.1f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", s.generation_time, (s.occupancy/(X_SIZE * Y_SIZE)),s.temperature, s.J, WP, BP, RP, GP,(X_SIZE*Y_SIZE), (WP+BP+RP+GP)); 			//Imprime en el archivo
  fclose (fp);
    if (p != NULL)
        {
        gchar *str = g_strdup_printf ("T_%.2f_J_%.2f_.png", s.temperature, s.coupling);
        cairo_surface_t *surf = cairo_image_surface_create (CAIRO_FORMAT_RGB24, X_SIZE, Y_SIZE);
        cairo_t *cr = cairo_create (surf);
        gdk_cairo_set_source_pixbuf (cr, p, 0, 0);
        cairo_paint (cr);
        cairo_surface_write_to_png (surf, str);
        g_print ("Screenshot saved as T_%.2f_J_%.2f.png\n", s.temperature, s.coupling);
        }
    else
        {
        g_print ("Unable to get the screenshot.\n");
        }
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
          put_pixel (p, (int) x, (int) y, (guchar) 255, (guchar) 255, (guchar) 255, 255);
          break;
        case 1:	/* Occupied site. Paint black */
          put_pixel (p, (int) x, (int) y, (guchar) 48, (guchar) 48, (guchar) 48, 255);
          break;
        case 2:	/* Occupied, +1 site. Paint red */
          put_pixel (p, (int) x, (int) y, (guchar) 249, (guchar) 237, (guchar) 105, 255);
          break;
        case 3:	/* Occupied, -1 site. Paint green */
          put_pixel (p, (int) x, (int) y, (guchar) 106, (guchar) 44, (guchar) 112, 255);
          break;
        }
      }
    }
  gtk_image_set_from_pixbuf (GTK_IMAGE (data), GDK_PIXBUF (p));
  g_object_unref (p);
  }


/* Function to get the closest neighbors (separated by r sites) of a given site in the lattice */
void get_closest_neighbors(int x, int y, int r, int* neighbors)
  {
  int n = 0;
  for (int i = -r; i <= r; i++) 
    {
    for (int j = -r; j <= r; j++) 
      {
      if (i == 0 && j == 0) {continue;}
      neighbors[n++] = s.lattice_configuration[(int) (X_SIZE + x + i) % X_SIZE][(int) (Y_SIZE + y + j) % Y_SIZE];
      }
    }
  }


/* Function used to compute the energy value of the site located at (x, y) */
int compute_energy (int x, int y)
  {
  int energy;
  int spin;
  int neighborhood_configuration = 0;
  int num_neighbors = 4*s.influence_radius*(s.influence_radius + 1);
  int neighbors_arr[num_neighbors];
  if      (s.lattice_configuration[x][y] == 2){spin = 1;}
	else if (s.lattice_configuration[x][y] == 3){spin =-1;}

  get_closest_neighbors (x, y, s.influence_radius, neighbors_arr);
  int n = 0;
  for (int n = 0; n < num_neighbors; n++) 
    {
    if (neighbors_arr[n] == 2) {neighborhood_configuration += 1;}
    else if (neighbors_arr[n] == 3) {neighborhood_configuration -= 1;}
    }
  energy = s.coupling * spin * neighborhood_configuration;
  return energy;
  }


/* Update function which simulates the stochastic process and updates the
   configuration of the lattice */
int update_lattice (gpointer data)
  {
  // int random_neighbor;
  int random_neighbor_state;
  int random_spin_value;
  int spin_energy, spin_energy_diff;
  double transition_probability;
  int random_x_coor, random_y_coor;

  int num_neighbors = 4*s.influence_radius*(s.influence_radius + 1);
  int neighbors[num_neighbors];
  int sites = 0;
  for (sites; sites < (int) (Y_SIZE*X_SIZE); sites++)
    {
    /* Pick a random focal site */
    random_x_coor = (int) floor (genrand64_real3 ()*X_SIZE);
    random_y_coor = (int) floor (genrand64_real3 ()*Y_SIZE);
    switch (s.lattice_configuration[random_x_coor][random_y_coor])
      {
      case 0: /* Site is empty */
      /* Chose a random neighbor from the 4 posible ones */
        get_closest_neighbors (random_x_coor, random_y_coor, s.influence_radius, neighbors);
        random_neighbor_state = neighbors[(int) floor (genrand64_real3 () * num_neighbors)];
        /* If its random neighbor is occupied: put a 1 copy at the focal site with probability brith_rate * dt */
        if (genrand64_real2 () < s.birth_rate)
          {
          if (random_neighbor_state == 1)
            {
            s.lattice_configuration[random_x_coor][random_y_coor] = 1; s.occupancy ++;
            }
          else if (random_neighbor_state == 2)
            {
            s.lattice_configuration[random_x_coor][random_y_coor] = 2; s.occupancy ++; s.up ++;
            }
          else if (random_neighbor_state == 3)
            {
            s.lattice_configuration[random_x_coor][random_y_coor] = 3; s.occupancy ++; s.down ++;
            }
          }
        break; /* break case 0 */
      case 1: /* If a particle is present at the focal site, it can die with probability dead_rate* dt */
        if (genrand64_real2 () < s.death_rate)
          {
          s.lattice_configuration[random_x_coor][random_y_coor]= 0;
          s.occupancy --;
          }
        else
          {
          if (genrand64_real2 () < s.differentiation_rate)
            {
            random_spin_value = (int) floor (genrand64_real3 ()*2);
            if (random_spin_value == 0)
              {
              s.lattice_configuration[random_x_coor][random_y_coor]= 2; s.up ++;
              }
            else if (random_spin_value == 1)
              {
              s.lattice_configuration[random_x_coor][random_y_coor]= 3; s.down ++;
              }
            }
          }
        break; /* break case 1 */

      case 2: /* Focal point is in state up (+1) */
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
            s.lattice_configuration[random_x_coor][random_y_coor] = 3;
            s.up --;
            s.down ++;
            }
          }
        break;
      case 3: /* Focal point is in state down (-1) */
        if (genrand64_real2 () < s.death_rate)
          {
          s.lattice_configuration[random_x_coor][random_y_coor]= 0;
          s.occupancy --;
          }
        else
          {
          spin_energy = compute_energy (random_x_coor, random_y_coor);
          spin_energy_diff = -2 * spin_energy;
          transition_probability = exp (-spin_energy_diff/s.temperature);
          if (spin_energy_diff < 0 || 
              genrand64_real2 () < transition_probability)
            {
            s.lattice_configuration[random_x_coor][random_y_coor] = 2;
            s.down --;
            s.up ++;
            }
          }
        break;
      }
    }
  s.generation_time ++;
  paint_lattice (data);
  //g_print ("\tGeneration:\t%u\t\tOccupancy:\t%f\n", s.generation_time, (s.occupancy/(X_SIZE * Y_SIZE)));
  g_print ("\tGen:\t%u\n", s.generation_time);
  // if (s.generation_time == 200) 			//Si el tiempo generacional alcanza un cierto número, la simulación se detiene
  //   {
  //   stop_simulation (data);
  //   }
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
  /* Set an occupied site in the middle of the lattice */
  s.lattice_configuration[X_SIZE/2][Y_SIZE/2] = 1;
  s.occupancy ++;

  s.initialized = TRUE;
  s.generation_time = 0;
  paint_lattice (data);
  g_print ("Lattice initialized\n");
  }


/* Callback to start simulation */
static void start_simulation (GtkWidget *button, gpointer data)
  {
  if(!s.running && s.initialized)
    {
    s.run = g_idle_add ((GSourceFunc) time_handler, GTK_IMAGE (data));
    s.running = TRUE;
    g_print ("Simulation started\n");
    }
  }


/* Callback to stop simulation */
static void stop_simulation (GtkWidget *button, gpointer data)
  {
  if (s.running)
    {
    g_source_remove (s.run);
    s.running = FALSE;
    g_print ("Simulation Stopped\n");
    }
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
  gchar *str = g_strdup_printf ("coupling = %d", (int) pos);
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

  /* Initialize Mersenne Twister algorithm for random number genration */
  unsigned int seed = (unsigned int) time (NULL);
  init_genrand64 (seed);

  /* Set default parameters of the simulation */
  s.influence_radius = 1;
  s.birth_rate = 1.00;
  s.death_rate = 0.00;
  s.differentiation_rate = 0.10;
  s.temperature = 0.00000001;
  s.coupling = -1.00;

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
  influence_radius_scale = gtk_scale_new_with_range (GTK_ORIENTATION_HORIZONTAL, 1, 2, 1);
  influence_radius_label = gtk_label_new ("radius");
  g_signal_connect (influence_radius_scale, "value-changed", G_CALLBACK (influence_radius_scale_moved),  influence_radius_label);
  gtk_grid_attach (GTK_GRID (grid), influence_radius_scale, 0, 0, 1, 1); // position (0,0) spanning 2 col and 1 row
  gtk_grid_attach (GTK_GRID (grid),  influence_radius_label, 1, 0, 1, 1); // position (1,0) spanning 3 col and 1 row

  /* Birth rate scale slide bar */
  birth_rate_scale = gtk_scale_new_with_range (GTK_ORIENTATION_HORIZONTAL, 0.0, 1.0, 0.01);
  birth_rate_label = gtk_label_new ("beta");
  g_signal_connect (birth_rate_scale, "value-changed", G_CALLBACK (birth_rate_scale_moved), birth_rate_label);
  gtk_grid_attach (GTK_GRID (grid), birth_rate_scale, 0, 1, 1, 1); // position (0,0) spanning 2 col and 1 row
  gtk_grid_attach (GTK_GRID (grid), birth_rate_label, 1, 1, 1, 1); // position (1,0) spanning 3 col and 1 row

  /* Death rate scale slide bar */
  death_rate_scale = gtk_scale_new_with_range (GTK_ORIENTATION_HORIZONTAL, 0.0, 1.0, 0.01);
  death_rate_label = gtk_label_new ("delta"); //LABEL to be shown J
  g_signal_connect (death_rate_scale, "value-changed", G_CALLBACK (death_rate_scale_moved), death_rate_label);
  gtk_grid_attach (GTK_GRID (grid), death_rate_scale, 2, 1, 1, 1); // position (3,0) spanning 2 col and 1 row
  gtk_grid_attach (GTK_GRID (grid), death_rate_label, 3, 1, 1, 1); // position (4,0) spanning 3 col and 1 row

  /* Temperature (T) scale slide bar */
  temperature_scale = gtk_scale_new_with_range (GTK_ORIENTATION_HORIZONTAL, 0.001, 10, 0.001);
  temperature_label = gtk_label_new ("temperature"); //LABEL to be shown T
  g_signal_connect (temperature_scale, "value-changed", G_CALLBACK (temperature_scale_moved), temperature_label);
  gtk_grid_attach (GTK_GRID (grid), temperature_scale, 0, 2, 1, 1); // position (0,0) spanning 2 col and 1 row
  gtk_grid_attach (GTK_GRID (grid), temperature_label, 1, 2, 1, 1); // position (1,0) spanning 3 col and 1 row

  /* Coupling (J) scale slide bar */
  coupling_scale = gtk_scale_new_with_range (GTK_ORIENTATION_HORIZONTAL, -1, 1, 1);
  coupling_label = gtk_label_new ("coupling"); //LABEL to be shown J
  g_signal_connect (coupling_scale,"value-changed", G_CALLBACK (coupling_scale_moved), coupling_label);
  gtk_grid_attach (GTK_GRID (grid), coupling_scale, 2, 2, 1, 1); // position (3,0) spanning 2 col and 1 row
  gtk_grid_attach (GTK_GRID (grid), coupling_label, 3, 2, 1, 1); // position (4,0) spanning 3 col and 1 row

  /* Pixel buffer @ start up and default canvas display */
  pixbuf = gdk_pixbuf_new (GDK_COLORSPACE_RGB, 0, 8, X_SIZE, Y_SIZE);
  image_lattice = gtk_image_new_from_pixbuf (pixbuf);
  paint_a_background (image_lattice);
  gtk_grid_attach (GTK_GRID (grid), image_lattice, 0, 3, 5, 1); // position (0,1) spanning 5 col and 1 row)

  /* ----------------------------  INIT BUTTON  ----------------------------- */
  button = gtk_button_new_with_label ("Initialize");
  g_signal_connect (button, "clicked", G_CALLBACK (init_lattice), GTK_IMAGE (image_lattice));
  gtk_grid_attach (GTK_GRID (grid), button, 0, 4, 1, 1); // position (0,3) spanning 1 col and 1 row)
  /* ----------------------------  START BUTTON  ---------------------------- */
  button = gtk_button_new_with_label ("Start");
  g_signal_connect (button, "clicked", G_CALLBACK (start_simulation), GTK_IMAGE (image_lattice));
  gtk_grid_attach (GTK_GRID(grid), button, 1, 4, 1, 1); // position (1,3) spanning 1 col and 1 row)
  /* ----------------------------  STOP BUTTON  ----------------------------- */
  button = gtk_button_new_with_label ("Stop");
  g_signal_connect (button, "clicked", G_CALLBACK(stop_simulation), NULL);
  g_signal_connect (button, "clicked", G_CALLBACK(count), NULL);
  gtk_grid_attach (GTK_GRID(grid), button, 2, 4, 1, 1); // position (2,3) spanning 1 col and 1 row)
  /* ----------------------------  QUIT BUTTON  ----------------------------- */
  button = gtk_button_new_with_label ("Quit");
  g_signal_connect_swapped (button, "clicked", G_CALLBACK (gtk_widget_destroy), window);
  gtk_grid_attach (GTK_GRID (grid), button, 4, 4, 1, 1); // position (4,3) spanning 1 col and 1 row)

  // Show the window and all widgets
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


/* main function */
int main (int argc, char **argv)
  {
  GtkApplication *app;
  int status;
  app = gtk_application_new ("fflab.cpim_01", G_APPLICATION_DEFAULT_FLAGS);
  g_signal_connect (app, "activate", G_CALLBACK (activate), NULL);
  status = g_application_run (G_APPLICATION (app), argc, argv);
  g_object_unref (app);
  return status;
  }

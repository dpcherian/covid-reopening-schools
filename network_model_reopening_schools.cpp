#include <iostream>
#include <stdio.h>
#include <random>
#include <algorithm>
#include <math.h>
#include <time.h>
#include <array>
#include <stdlib.h>
#include <string.h>

using namespace std;

/* Random Number Generators **********************************/

std::random_device rd;
std::mt19937 gen(rd());
// std::mt19937 gen(-1);

double uniform()        // Generate a uniform random number between [0,1)
{
    std::uniform_real_distribution<> d(0.0,1.0);
    return d(gen);
}

int randint(int b)      // Generate a uniform random number in [0, b)
{
    std::uniform_int_distribution<> d(0,b-1);
    return d(gen);
}

void shuffle(int start, int end, int *array)
{
    if (end-start > 1)
    {
        for (int i = start; i < end; i++)
        {
          int j = i + rand() / (RAND_MAX / (end - i) + 1);
          int t = array[j];
          array[j] = array[i];
          array[i] = t;
        }
    }
}

void shuffle(int start, int end, vector<int> array)
{
    if (end-start > 1)
    {
        for (int i = start; i < end; i++)
        {
          int j = i + rand() / (RAND_MAX / (end - i) + 1);
          int t = array[j];
          array[j] = array[i];
          array[i] = t;
        }
    }
}


/*****************************************************************/


/* Global variables **********************************************/

/*
Below are certain parameters used by the code. They are divided into three sections:
1. User defined parameters: users can change these as they see fit. These decide the
   population size, the initial conditions of the code, contact parameters, vaccination
   parameters and so on.

2. Constants used by the code: These are set constants that the code requires to run.
   Changing any of these values will either break the code or make it run unpredictably.
   Do not edit them unless you know what you're doing.

3. Variables set during runtime: These are variables that will be reset every run, and
   are initialised globally for ease-of-access. There is no need to edit them.
*/

// *****************************************************************//
// ****************** 1. USER DEFINED VARIABLES ********************//
// *****************************************************************//

double dt = 0.1;
char output_folder[100] = "tmp";
bool print_log = true;

const int mult = 1;              // Multiplier to scale up problem

const int n_pop = mult*20316;    // Total population
const int n_loc = mult*6623;     // Total number of locations
const int n_net = mult*123;      // Total number of networks
const int n_hospitals = mult*1;  // Set the number of hospitals.

const int start_schools = n_hospitals; // Where school locations start (DON'T change this!)
const int n_schools = mult*1;          // Set the number of schools (1)

const int classrooms = 100;            // Number of classrooms per school

double initial_asymptomatic_fraction = 1.00/100; // Initial fraction of asymptomatic
double initial_recovered_fraction    = 30.0/100; // Initial fraction of recovered
double initial_vaccinated_fraction   = 20.0/100; // Initial fraction of vaccinated


double Cpars[] = {1,  1,   1,    1,   0.1}; // Contact parameters for transmission
               // SA  SP  S-MI  S-SI  S-H

/* Vaccination parameters ********************/

double dvr = 0.0; // Daily Vaccination Rate (units: percentage of total population per day)
double infection_reduction       = 1.0 - 0.441;
double gamma_fractional_increase = 1.0 - 0.598;
double transmission_reduction    = 0.4;

int vaccination_strategy = 0;  // Can be 0 (unset), ascending by age (1), or descending by age (-1)

/* School parameters *************************/

double unlockSchoolsOn = 100;	// units: days
bool start_with_schools_locked = true;          // Don't allow people to go to schools. Default: true (start schools locked)
bool shuffle_teachers = false;                  // Allow teachers move between classrooms. Default: false (don't let them move)


// *****************************************************************//
// ****************** 2. CONSTANTS USED BY CODE ********************//
// *****************************************************************//

const char states[][3] = {"S","E", "A","P", "MI","SI","R","H", "D"};
                        // 0   1    2    3   4    5    6   7    8
const int n_states = sizeof(states)/sizeof(states[0]); // Number of states

const int S = 0;
const int E = 1;
const int A = 2;
const int P = 3;
const int MI= 4;
const int SI= 5;
const int R = 6;
const int H = 7;
const int D = 8;

const int RAT = 0;
const int PCR = 1;

const int person_attr = 6;  // Number of attributes a person can be assigned in the pop array (state, age, "something", home location, work location, current location)

const int s = 0;   // Index of "state" attribute

const int a = 1;   // Index of "age" attribute
const int x = 2;   // Index of "some factor" attribute

const int h = 3;   // Index of "home loc" attribute
const int w = 4;   // Index of "work loc" attribute
const int c = 5;   // Index of "current loc" attribute

const int max_age     = 9;
const int min_adult_age=2;

const int ascending = 1;       // Possible vaccination strategy (ascending)
const int descending=-1;       // Possible vaccination strategy (descending)

const int op_width = 1 + n_states + 2;                // Width of output array: <time>(1) <Number of states> <Vaccines Administered> <Background Seropositivity>
const int age_op_width = 1 + (n_states+2)*max_age;    // Width of age-wise output array <time>(1) (<Number of states> <Vaccines administered> <Background Seropositivity> ) x max_age


// *****************************************************************//
// ************** 3. VARIABLES SET DURING RUN-TIME *****************//
// *****************************************************************//

double rate[max_age][n_states][n_states] = {};
int n_asym = 0;      // To store: initial number of asymptomatics
int n_rec  = 0;      // To store: initial number of recovered

/* Population arrays ******************************/

int pop[n_pop][person_attr] = {};

int n[n_states];
int n_per_location[n_loc][n_states];
int n_per_age[max_age][n_states];

bool is_adult[n_pop]     = {};


/* Location arrays *******************************/

vector<vector<int>> people_linked_to(n_loc);
int len_ppl_linked_to[n_loc] = {};

vector<vector<vector<int>>> n_per_room(n_loc);
bool is_school[n_loc] = {};
int n_rooms[n_loc]    = {};

int current_room[n_pop]  = {};
int assigned_room[n_pop] = {};

/* Vaccination arrays ***************************/

int ramp_up_time = 90+14;

bool is_vaccinated[n_pop] = {};
int vaccinated_on[n_pop] = {};
int vaccs_per_state[n_states] = {};

int vaccines_administered = 0;
int agewise_vaccines_administered[max_age] = {};


/************************************************/


void readcsv(char *filename) // read in the synthetic population from csv file
{
int num_to_prevaccinate = (int)(initial_vaccinated_fraction*n_pop);

FILE* fp = fopen(filename, "rt");
if (fp==NULL) {printf("No file with name %s!\n",filename);exit(2);}

char *buffer,*tofree,*value;

int row   = -1;    // Initial row is -1, to ignore the headers
int column = 0;

while (row<n_pop)
	{
	tofree=buffer=(char *)malloc(1024*sizeof(char));

	if (fgets(buffer,1024, fp)==NULL) {printf("End of file reached, only %d rows read\n",row);}

	column = 0;

	if (row == -1) {row++;continue;} // because you need to discard the first row

	// Get each column as a string, convert to number for required columns
	while ((value=strsep(&buffer, ","))!=NULL)
		{
      pop[row][column] = atoi(value);
      if(column==a){pop[row][column] = (int)pop[row][column]/10;} // If the column is the age, divide it by 10 to get to the decade

      if(column==s){                                              // If the column is the state, set the initial state
        double r = uniform();                                     // Random number to decide initial state of agents
        if(r<initial_asymptomatic_fraction){pop[row][column] = A; n_asym++;}
        else if(r>initial_asymptomatic_fraction && r<initial_asymptomatic_fraction+initial_recovered_fraction){ pop[row][column]=R; n_rec++;}
        else{pop[row][column] = S;}
      }
		column++;
		}

    // Prevaccinate 80% of the 60+ individuals
    is_vaccinated[row] = (num_to_prevaccinate>0 && pop[row][a] >= 6 && uniform()<0.8) ? true : false;
    vaccinated_on[row] = is_vaccinated[row] ? -1000 : 1000;
    if(is_vaccinated[row]){
      num_to_prevaccinate--;
      vaccs_per_state[pop[row][s]]++;
      agewise_vaccines_administered[pop[row][a]]++;
      vaccines_administered++;
    }

    int state       = pop[row][s];
    int age         = pop[row][a];
    int work_loc    = pop[row][w];
    int current_loc = pop[row][c];



    n_per_location[current_loc][state] += 1;  // Increment number of disease-state in current location
    n_per_age[age][state]              += 1;  // Increment number of disease state per age

    if(age>=min_adult_age){is_adult[row] = true;}

    assigned_room[row] = randint(n_rooms[work_loc]);
    current_room[row]  = 0;

    n_per_room[current_loc][current_room[row]][state] += 1; // Increase number of current state in current room of current location

	row++;

	free(tofree); // I don't fully understand why I need to keep freeing buffer and reallocating it, and why it needs to be done like this, but anything else doesn't seem to work!
	}

  //*********** Prevaccinate randomly if any doses are left ********//

  while(num_to_prevaccinate>0){                            // While "pre"vaccines are still available
    int p = randint(n_pop);                                // Choose a person at random
    if(!is_vaccinated[p]){                                 // If they have not been vaccinated
      is_vaccinated[p] = is_adult[p] ? true : false;       // and are adults, vaccinate them
      vaccinated_on[p] = is_vaccinated[p] ? -1000 : 1000;
      if(is_vaccinated[p]){
        num_to_prevaccinate--;
        vaccs_per_state[pop[p][s]]++;
        agewise_vaccines_administered[pop[p][a]]++;
        vaccines_administered++;
      }
    }
  }

  //************ Done prevaccinating *********************************//

fclose(fp);

//***** CREATE ARRAY OF PEOPLE-LINKED-TO LOCATION ***************************************************//

for(int i=0;i<n_loc;i++){
  people_linked_to[i] = vector<int>();          // Empty vector for each location

  int counter = 0;                              // Counter to count people in location
  for(int j=0; j<n_pop;j++){
      if(pop[j][h] == i || pop[j][w] == i){     // If the person's home or work location is i
          people_linked_to[i].push_back(j);     // Add them to this location in people_linked_to
          counter++;
      }
  }

  len_ppl_linked_to[i] = counter;               // Array containing the total number of people linked to each location
}

//***** DONE CREATING ARRAY OF PEOPLE-LINKED-TO LOCATION ********************************************//


//***** CREATE n ARRAY *************************************//

  n[S] = n_pop-n_asym-n_rec;
  n[E] = 0;
  n[A] = n_asym;
  n[P] = 0;
  n[MI]= 0;
  n[SI]= 0;
  n[R] = n_rec;
  n[H] = 0;
  n[D] = 0;

//********* DONE CREATING n ARRAY **************************//

}


void createPopulation(char *filename){

  // Reset arrays and some variables, and load population from csv file
  n_asym = 0;
  n_rec  = 0;
  vaccines_administered = 0;

  for(int age=0;age<max_age;age++){agewise_vaccines_administered[age]=0;}

  // *************************************************//
  // *************** LOCATION SETUP ******************//
  // *************************************************//

  // Assign rooms to each location, and initialise the n_per_room vector.
  // Reset all classroom numbers

  for(int i=0; i<n_loc;i++){

    // Reset n_per_location
    for(int j=0;j<n_states;j++){n_per_location[i][j]=0;}

    // Decide if a location is a school or not
    if(i<start_schools+n_schools && i>=start_schools){
      is_school[i] = true;
    }
    else{ is_school[i] = false; }

    // Assigning quantities per room:
    if(!is_school[i]){ n_rooms[i] = 1;} else{ n_rooms[i] = classrooms;}  // If it's not a school, assign only 1 room, else assign `classrooms` rooms.

    n_per_room[i] = vector<vector<int>>(n_rooms[i]);
    for(int room=0; room<n_rooms[i];room++){          // For each room in locations
      n_per_room[i][room] = vector<int>(n_states,0);  // assign a vector which -- for each room -- contains the number of individuals per state
    }
  }


  // *************************************************//
  // ****************** AGE SETUP ********************//
  // *************************************************//

  // Reset n_per_age
  for(int i=0;i<max_age;i++){for(int j=0; j<n_states;j++){n_per_age[i][j]=0;}}


  // *************************************************//
  // **************** PERSON SETUP *******************//
  // *************************************************//

  readcsv(filename);

  // *************************************************//
  // **************** PRINT POPULATION ***************//
  // *************************************************//

  // for(int i=0; i<n_pop;i++){
  //   for(int j=0;j<person_attr;j++){
  //     printf("%i, ", pop[i][j]);
  //   }
  //   printf("%s",is_vaccinated[i]?"1":"0");
  //   printf("\n");
  // }
  // exit(2);

}


char outputFilename[1000];
void writetofile(int output[][op_width], int age_output[][age_op_width], int tf, double time_taken, int details[9],int iter){
  // FOR REFERENCE: int details[] = {locations_moved, hcw_recovered, hcw};

  // Write output to a file

  errno = 0;      // Variable to store error number in case file open does not work

  double rn1 = uniform();
  double rn2 = uniform();

  sprintf(outputFilename,"./%s/Total_IR_%g_IV_%g_VR_%g_%s_UnlockSchoolsOn_%g_%lf%lf-%i.txt",output_folder,initial_recovered_fraction*100, initial_vaccinated_fraction*100, dvr,vaccination_strategy==ascending ? "Ascending" : vaccination_strategy == descending ? "Descending" : "Unset",unlockSchoolsOn,rn1,rn2,iter);
  FILE *fpt=(FILE *)fopen(outputFilename,"wt"); // Open the file to print output

  if(fpt){
    if(print_log){
      fprintf(fpt,"###### TEST LOG ####################\n");
      fprintf(fpt,"# Time taken               : %.2f s\n",time_taken);
      fprintf(fpt,"# Locations Moved in total : %d\n", details[0]);
      fprintf(fpt,"# Total recovered HCW      : %d\n", details[1]);
      fprintf(fpt,"# Total HCW                : %d\n", details[2]);

      fprintf(fpt,"# Vaccination strategy     : %s\n", vaccination_strategy==ascending ? "Ascending" : vaccination_strategy == descending ? "Descending" : "Unset");
      fprintf(fpt,"# Total Vaccines Given     : %d\n", vaccines_administered);
      fprintf(fpt,"# Schools unlocked on day  : %g\n", unlockSchoolsOn);

      fprintf(fpt,"# Rate Array: \n");
      for(int i=0; i<n_states;i++){fprintf(fpt,"# ");for(int j=0; j<n_states;j++){fprintf(fpt,"%5g ", rate[0][i][j]);}fprintf(fpt,"\n");} // NOTE! : Change to correct rates
      fprintf(fpt,"###### END LOG #####################\n");
      fprintf(fpt,"#\n");
      fprintf(fpt, "#%s %s %s %s %s %s %s %s %s %s %s %s\n","Day","nS","nE","nA","nP","nMI","nSI","nR","nH","nD","nV","nBS");
    }

    for(int i=0;i<tf;i++){
      for(int j=0; j<op_width;j++){
        fprintf(fpt, "%i ", output[i][j]);
      }
      fprintf(fpt, "\n");
    }

    fclose(fpt);  // Close the file

  }
  else{printf("File error, %d\n",errno);}

  errno = 0;      // Variable to store error number in case file open does not work
  char ageOutputFilename[1000];

  sprintf(ageOutputFilename,"./%s/AgeStratified_IR_%g_IV_%g_VR_%g_%s_UnlockSchoolsOn_%g_%lf%lf-%i.txt",output_folder,initial_recovered_fraction*100, initial_vaccinated_fraction*100, dvr,vaccination_strategy==ascending ? "Ascending" : vaccination_strategy == descending ? "Descending" : "Unset",unlockSchoolsOn,rn1,rn2,iter);
  FILE *fpt1=(FILE *)fopen(ageOutputFilename,"wt"); // Open the file to print output

  if(fpt){
    if(print_log){
      fprintf(fpt1,"###### TEST LOG ####################\n");
      fprintf(fpt1,"# Time taken               : %.2f s\n",time_taken);
      fprintf(fpt1,"# Main file                : %s \n",outputFilename);

      fprintf(fpt1,"# Locations Moved in total : %d\n", details[7]);
      fprintf(fpt1,"# Total recovered HCW      : %d\n", details[8]);
      fprintf(fpt1,"# Total HCW                : %d\n", details[9]);

      fprintf(fpt1,"# Vaccination strategy     : %s\n", vaccination_strategy==ascending ? "Ascending" : vaccination_strategy == descending ? "Descending" : "Unset");
      fprintf(fpt1,"# Total Vaccines Given     : %d\n", vaccines_administered);
      fprintf(fpt1,"# Schools unlocked on day  : %g\n", unlockSchoolsOn);

      fprintf(fpt1,"# Rate Array: \n");
      for(int i=0; i<n_states;i++){fprintf(fpt1,"# ");for(int j=0; j<n_states;j++){fprintf(fpt1,"%5g ", rate[0][i][j]);}fprintf(fpt1,"\n");} // NOTE! : Change to correct rates
      fprintf(fpt1,"###### END LOG #####################\n");
      fprintf(fpt1,"#\n");

      fprintf(fpt1, "#Day ");
      for(int j=0;j<max_age;j++){
        for(int i=0;i<n_states;i++){
          fprintf(fpt1, "%s%i ",states[i],(j+1)*10);
        }
        fprintf(fpt1, "%s%i ","V",(j+1)*10);
        fprintf(fpt1, "%s%i ","BS",(j+1)*10);
      }
      fprintf(fpt1,"\n");
    }

    for(int i=0;i<tf;i++){
      for(int j=0; j<(n_states+2)*max_age + 1;j++){
        fprintf(fpt1, "%i ", age_output[i][j]);
      }
      fprintf(fpt1, "\n");
    }

    fclose(fpt1);  // Close the file

  }
  else{printf("File error, %d\n",errno);}

}

// do one run of the epidemic for tf days
void run(int tf, double start_vacc, double dvr, double unlockSchoolsOn, int iter){

  clock_t start, end;                                     // Measuring how long the function takes to run
  double cpu_time_used;
  start = clock();

  int locations_moved    = 0;

  double r[n_states][n_states] = {};                      // Array to store rates per event

  double alphaH = 1 - Cpars[4];                           // Quantity by which V is reduced per hospitalised individual
  double alphaT = 1 - transmission_reduction;             // Effective reduction in transmission per person due to vaccination

  double t = 0.0;
  int day  = 0;
  bool midday_move_completed = false;

  bool first_move_done =false;
  bool second_move_done=false;
  bool third_move_done =false;

  bool lock_schools = start_with_schools_locked;           //Reset schools to start off locked

  int output[tf+1][op_width];         // Output array, to be printed to file
  int age_output[tf+1][age_op_width]; // Output array per age


  // Computing Background Seropositivity **************************************//

  double some_protection = 0;
  double agewise_some_protection[max_age] = {};
  double agewise_populations[max_age] = {};

  for(int p=0;p<n_pop;p++){
    agewise_populations[pop[p][a]]++;
    if(pop[p][s]==S && !is_vaccinated[p]){some_protection++; agewise_some_protection[pop[p][a]]++;}
  }

  double bs = n_pop - some_protection;
  double agewise_bs[max_age] = {};

  for(int age=0;age<max_age;age++){
    agewise_bs[age] = agewise_populations[age] - agewise_some_protection[age];
  }

  // First line of output
  output[day][0] = day; for(int q=0;q<n_states;q++){output[day][q+1]=n[q];}
  output[day][n_states+1] = vaccines_administered;
  output[day][n_states+2] = bs;

  age_output[day][0] = day; int counter = 1; for(int r=0;r<max_age;r++){for(int q=0;q<n_states;q++){age_output[day][counter]=n_per_age[r][q];counter++;} age_output[day][counter]=agewise_vaccines_administered[r];counter++;
  age_output[day][counter]=agewise_bs[r];counter++;}

  double exit_rate[n_states] = {};

  while(t<tf){

    if(t>unlockSchoolsOn && lock_schools==true){  // Check if schools need to be unlocked
      lock_schools=false;
    }


    //******************** Moving people around deterministically (HOME TO WORK) ********************//

    if(midday_move_completed == false && t - day >  0.5){
      midday_move_completed = true;

      for(int i=0;i<n_pop;i++){

        if(pop[i][s]!=H && pop[i][s]!=D){
          locations_moved++;
          int home_loc = pop[i][h];
          int work_loc = pop[i][w];
          if(pop[i][c]==home_loc)if(lock_schools==false || ( lock_schools==true && !is_school[work_loc])){ pop[i][c] = work_loc; n_per_location[home_loc][pop[i][s]]--; n_per_location[work_loc][pop[i][s]]++; n_per_room[home_loc][current_room[i]][pop[i][s]]--; current_room[i]=assigned_room[i]; n_per_room[work_loc][assigned_room[i]][pop[i][s]]++; }
        }
      }
    }

    //******************** Move teachers between classrooms in schools ********************//

    if(lock_schools==false && ((first_move_done==false && t-day>0.5 && t-day<0.5+0.17) || (second_move_done==false && t-day>0.5+0.17 && t-day<0.5+0.33) || (third_move_done==false && t-day>0.5+0.33 && t-day<0.5+0.5)) ){

      for(int i=start_schools; i<start_schools+n_schools; i++){ // For every school location

        int num_students_at_school = 0;

        for(int j=0; j<len_ppl_linked_to[i];j++){   // Loop over all people in this school
          int person = people_linked_to[i][j];
          int person_state = pop[person][s];
          if(pop[person][c] == i){                  // If the person is currently in the school, add them to a list.
            num_students_at_school++;

            if(first_move_done==false){             // If they're just entering the school, allow for possibility to move to a random room
              int set_classroom = assigned_room[person];
              if(set_classroom != current_room[person]){ printf("Set classroom is %i current room is %i\n",set_classroom,current_room[person] );}
              current_room[person] = set_classroom; // Set the person to their assigned classroom.
            }
            else if(is_adult[person] && shuffle_teachers){
              int from_room = current_room[person];  // Get current room of teacher in school i
              int to_room = randint( n_rooms[i] );   // Choose random room from the number of rooms in this school

              current_room[ person ] = to_room;

              n_per_room[i][from_room][person_state]--;
              n_per_room[i][to_room][person_state]++;
            }

            if(pop[person][s] == H){printf("ERROR! Hospitalised in School!\n"); exit(12);}
          }
        }

        int temp_tot = 0;

        for(int classroom=0; classroom<n_rooms[i]; classroom++){
          for(int state = 0; state<n_states; state++){
              temp_tot += n_per_room[i][classroom][state];
          }
        }


        if(temp_tot!=num_students_at_school){printf("Error! Sum of n_per_classroom %4i, num students %4i \n",temp_tot, num_students_at_school ); exit(2);}

      }// End loop over schools

      if(first_move_done==false){first_move_done=true; }
      else if(second_move_done==false){second_move_done=true;}
      else if(third_move_done==false){third_move_done=true;}
      else{printf("Error, all three moves are done, but movement is still occurring.\n");}

    }
	  
    // now loop over all people in all rooms and update their infection state if needed
    for(int i=0; i<n_loc; i++){

      for(int room=0; room<n_rooms[i];room++){

        int N = 0;

        int counter = 0;

        vector<int> ind;

        for(int state = 0; state<n_states;state++){N += n_per_room[i][room][state];}                 // Total number currently in this room

        for(int p=0;p<len_ppl_linked_to[i];p++){                                                     // For all people linked to this school
          int person = people_linked_to[i][p];
          if(pop[person][c]==i && current_room[person] == room){ ind.push_back(person); counter++;}  // If the person is in this school and is currently in this room, add them to ind
        }

        int sum=0;for(int state=0;state<n_states;state++){sum+=n_per_room[i][room][state];}

        if(counter!=N){printf("Error in location %i, N = %i, sum %i, counter = %i at time %lf\n",i,N,sum,counter,t );}

        if(N==0){continue;}

        int vacc_by_state_in_room[n_states]={}; // Number of vaccinated individuals in this room
        int tot_vacc = 0;

        for(int j = 0; j<N;j++){      // Go over people in this room
            int p = ind[j];           // index of person
            if(is_vaccinated[p]==true){vacc_by_state_in_room[pop[p][s]]++; tot_vacc++; }
        }

        shuffle(0,N,ind);      // Shuffle list of people currently in locations

        for(int j=0;j<N;j++){  // Loop over the people currently in the location

          int age = pop[ind[j]][a];

          double new_gamma = rate[age][E][A]/(rate[age][E][A]+rate[age][E][P]);    // Default vals: If the person isn't a vaccinated subject, f is just gamma
          double lambda_E = rate[age][E][A]+rate[age][E][P];
          double beta_multiplier = 1.0;                                            // Default beta-multiplier (reduction) = 1 unless vaccinated

          if(is_vaccinated[ind[j]]){
            int vday   = vaccinated_on[ind[j]];
            double this_gamma = rate[age][E][A]/(rate[age][E][A]+rate[age][E][P]);

            beta_multiplier = std::min(std::max(1 + (infection_reduction - 1)*(t - vday)/ramp_up_time, 1 + (infection_reduction - 1)), 1.0);
            new_gamma = std::min(std::max(this_gamma + (1 - (gamma_fractional_increase/infection_reduction)*(1 - this_gamma) - this_gamma)*(t - vday)/ramp_up_time, this_gamma), this_gamma + (1 - (gamma_fractional_increase/infection_reduction)*(1 - this_gamma) - this_gamma));
          }

          int people_here = 0; for(int state=0;state<n_states;state++){people_here += n_per_room[i][room][state];}

          double V = people_here - alphaH*n_per_room[i][room][H];   // Spatial damping parameter (adjusted by alphaH for hospitals)

          // S->E
          r[S][E] = beta_multiplier * rate[age][S][E] * 1/V * (Cpars[0]*(n_per_room[i][room][A]  - vacc_by_state_in_room[A]*alphaT) +    // SA
                                                               Cpars[1]*(n_per_room[i][room][P]  - vacc_by_state_in_room[P]*alphaT) +    // SP
                                                               Cpars[2]*(n_per_room[i][room][MI] - vacc_by_state_in_room[MI]*alphaT) +   // S-MI
                                                               Cpars[3]*(n_per_room[i][room][SI] - vacc_by_state_in_room[SI]*alphaT) +   // S-SI
                                                               Cpars[4]*n_per_room[i][room][H]);                                                                          // SH
          r[E][A]  =  new_gamma * lambda_E;      // E->A
          r[E][P]  = (1 - new_gamma) * lambda_E; // E->P
          r[A][R]  = rate[age][A][R];            // A->R
          r[P][MI] = rate[age][P][MI];           // P->MI
          r[P][SI] = rate[age][P][SI];           // P->SI
          r[MI][R] = rate[age][MI][R];           // MI->R
          r[SI][R] = rate[age][SI][R];           // SI->R
          r[SI][H] = rate[age][SI][H];           // SI->H
          r[H][R]  = rate[age][H][R];            // H->R
          r[H][D]  = rate[age][H][D];            // H->D

          exit_rate[S] = r[S][E];
          exit_rate[E] = r[E][A] + r[E][P];
          exit_rate[A] = r[A][R];
          exit_rate[P] = r[P][MI] + r[P][SI];
          exit_rate[MI]= r[MI][R];
          exit_rate[SI]= r[SI][R] + r[SI][H];
          exit_rate[R] = 0;
          exit_rate[H] = r[H][R]+r[H][D];
          exit_rate[D] = 0;

          if(r[S][E] <0 || r[E][A]<0 || r[E][P]<0){printf("Negative rates %lf %lf %lf\n", r[S][E], r[E][A],r[E][P]);}

          if(r[S][E] <0 || r[E][A]<0 || r[E][P]<0){
            double a = beta_multiplier;
            double b = rate[age][S][E];
            double c = 1/V;
            double d =  (Cpars[0]*(n_per_room[i][room][A]  - vacc_by_state_in_room[A]*alphaT) +    // SA
                         Cpars[1]*(n_per_room[i][room][P]  - vacc_by_state_in_room[P]*alphaT) +    // SP
                         Cpars[2]*(n_per_room[i][room][MI] - vacc_by_state_in_room[MI]*alphaT) +   // S-MI
                         Cpars[3]*(n_per_room[i][room][SI] - vacc_by_state_in_room[SI]*alphaT) +   // S-SI
                         Cpars[4]*n_per_room[i][room][H]);
            printf("Location %i Beta Multiplier %lf RateSE %lf 1/V %lf I/N %lf \n",i,a,b,c,d );
          }

          int from = pop[ind[j]][s];
          if(uniform()<exit_rate[from]*dt){                             // If the person is selected to change infection state
            double p = uniform();
            double temp=0;

            for(int to=0;to<n_states;to++){                             // Loop over possible "to" states

              temp += (r[from][to]/exit_rate[from]);

              if(p<temp){                                               // If such a transition must occur,
                pop[ind[j]][s] = to;                                    // Send this person to the "to" state.
                n_per_location[i][from]--; n_per_location[i][to]++;
                n[from]--; n[to]++;                                     // Change the values of n[from] and n[to]
                n_per_age[age][from]--;n_per_age[age][to]++;
                n_per_room[i][room][from]--; n_per_room[i][room][to]++;

                vacc_by_state_in_room[from]--;vacc_by_state_in_room[to]++;

                if(to==H){
                  // Move them to the hospital

                  int hosp = randint(n_hospitals);                          // Choose a hospital at random
                  while(hosp==pop[ind[j]][w]){hosp = randint(n_hospitals);} // Make sure a HCW is not hospitalised in the same hospital
                  int hosp_room = randint(n_rooms[hosp]);                   // Random room in the hospital to send person to

                  people_linked_to[hosp].push_back(ind[j]); len_ppl_linked_to[hosp]++; // Link individual permanently to that hospital

                  n_per_location[ pop[ind[j]][c] ][H]--;      // Decrement number of "from" in current location
                  n_per_room[ pop[ind[j]][c] ][room][H]--;
                  pop[ind[j]][c] = hosp;                      // Send them to a random hospitals
                  current_room[ind[j]] = hosp_room;
                  n_per_location[ pop[ind[j]][c] ][H]++;      // Increment number of "to" in current location
                  n_per_room[pop[ind[j]][c]][hosp_room][H]++; // Send them to a random room in the hospital.
                }

                else if(from==H && (to == R || to == D)){     // If they exit from H, move them home

                  int home = pop[ind[j]][h];
                  int home_room = randint(n_rooms[home]);

                   // Send recovered or dead who were hospitalised home
                  n_per_location[ pop[ind[j]][c] ][to]--;     // Decrement number of to (=R,D) in current locations
                  n_per_room[ pop[ind[j]][c] ][room][to]--;   // Decrement number of people in this room
                  pop[ind[j]][c] = home;                      // Send them home
                  current_room[ind[j]] = home_room;
                  n_per_location[ pop[ind[j]][c] ][to]++;     // Increment number of to (=R,D) in current locations
                  n_per_room[ pop[ind[j]][c] ][ home_room ][to]++; // Send them to a random room in their house
                }

                break; // Exit the "to" loop, move to next person
              } // End if condition (if p<temp)
            } // End loop over "to" states
          } // End if condition (infection state == from)
        } // End loop over people in location

      /************ END CHANGE STATE OF POP ***********/

      }

    } // End of location loop.

    /**************************** THE END OF DAY(S) ****************************/

    if(t>=day+1){

          //************ VACCINATION ******************//

          if(n[R]>=start_vacc/100 * n_pop){          // Start vaccination only after some fraction of recovered is attained

            int vaccs_available = dvr/100 * n_pop;   // Total number of available vaccines per day

            vector<int> vacc_targets;                // List of all target ages
            int len_vacc_targets=0;                  // Number of target ages

            if(vaccination_strategy == descending){
              if(t<0){
                vacc_targets.push_back(6);
                vacc_targets.push_back(7);
                vacc_targets.push_back(8);
                vacc_targets.push_back(9);

                len_vacc_targets = 4;
              }
              else if(t>0 && t<=30){
                vacc_targets.push_back(4);
                vacc_targets.push_back(5);
                vacc_targets.push_back(6);
                vacc_targets.push_back(7);
                vacc_targets.push_back(8);
                vacc_targets.push_back(9);

                len_vacc_targets = 6;
              }
              else{
                vacc_targets.push_back(2);
                vacc_targets.push_back(3);
                vacc_targets.push_back(4);
                vacc_targets.push_back(5);
                vacc_targets.push_back(6);
                vacc_targets.push_back(7);
                vacc_targets.push_back(8);
                vacc_targets.push_back(9);

                len_vacc_targets = 8;
              }
            }
            else if(vaccination_strategy == ascending){
              if(t<=30){
                vacc_targets.push_back(2);
                vacc_targets.push_back(3);

                len_vacc_targets = 2;
              }
              else if(t>30 && t<=60){
                vacc_targets.push_back(2);
                vacc_targets.push_back(3);
                vacc_targets.push_back(4);
                vacc_targets.push_back(5);

                len_vacc_targets = 4;
              }
              else{
                vacc_targets.push_back(2);
                vacc_targets.push_back(3);
                vacc_targets.push_back(4);
                vacc_targets.push_back(5);
                vacc_targets.push_back(6);
                vacc_targets.push_back(7);
                vacc_targets.push_back(8);
                vacc_targets.push_back(9);

                len_vacc_targets = 8;
              }
            }
            else{
              printf("Error! Vaccination strategy wrong or not set!\n");
              exit(2);
            }

            vector<int> list_of_vacc;         // List of targeted vaccinated individuals
            vector<int> list_of_others;       // List of random individuals for any remaining vaccines

            int n_to_vacc = 0;                // Number of targeted vaccinations
            int n_others  = 0;                // Number of other individuals

            int backup[n_states]={};

            for(int i=0;i<n_pop;i++){
              bool added = false;
              for(int j=0;j<len_vacc_targets;j++){
                int v_target = vacc_targets[j];
                if(pop[i][a]==v_target && pop[i][s]!=MI && pop[i][s]!=SI && pop[i][s]!=H && pop[i][s]!=D && !is_vaccinated[i]){
                  backup[pop[i][s]]++;
                  list_of_vacc.push_back(i); n_to_vacc++; added = true; break;
                }
              }
              if(added == false && pop[i][s]!=MI && pop[i][s]!=SI && pop[i][s]!=H && pop[i][s]!=D && !is_vaccinated[i]){
                list_of_others.push_back(i); n_others++;
              }
            }

            int vaccs_done_today = std::min(n_to_vacc,vaccs_available);
            shuffle(0,n_to_vacc,list_of_vacc);

            for(int i=0;i<vaccs_done_today;i++){    // Vaccinate chosen individuals
              int vi = list_of_vacc[i];
              is_vaccinated[vi]=true;
              vaccinated_on[vi] = day;
              vaccs_available--;
              vaccs_per_state[pop[vi][s]]++;

              vaccines_administered++;
              agewise_vaccines_administered[pop[vi][a]]++;
            }

            int random_vaccs_done_today=0;

            // If any vaccines are still available, and there are still people to vaccinate, vaccinate them randomly

            if(n_others>0 && vaccs_available>0){
               random_vaccs_done_today = std::min(n_others,vaccs_available);

              shuffle(0,n_others,list_of_others);

              for(int i=0;i<random_vaccs_done_today;i++){
                int vi = list_of_others[i];
                is_vaccinated[vi]=true;
                vaccinated_on[vi] = day;
                vaccs_available--;
                vaccs_per_state[pop[vi][s]]++;

                vaccines_administered++;
                agewise_vaccines_administered[pop[vi][a]]++;
              }
            }
          }


      //******************** Moving people around deterministically (WORK TO HOME) ********************//

      for(int i=0;i<n_pop;i++){
        if(pop[i][s]!=H && pop[i][s]!=D){

          locations_moved++;

          int home_loc = pop[i][h];
          int work_loc = pop[i][w];
          if(pop[i][c]==work_loc){ pop[i][c] = home_loc; n_per_location[work_loc][pop[i][s]]--; n_per_location[home_loc][pop[i][s]]++; n_per_room[work_loc][current_room[i]][pop[i][s]]--; current_room[i] = randint(n_rooms[home_loc]); n_per_room[home_loc][current_room[i]][pop[i][s]]++; }
        }
      }

      // Reset movement flags
      midday_move_completed = false;
      first_move_done       = false;
      second_move_done      = false;
      third_move_done       = false;

      //***********************************************************************************************//

      day++; // Increment the day

      // Computing Background Seropositivity
      double some_protection = 0;
      double agewise_some_protection[max_age] = {};
      double agewise_populations[max_age] = {};

      for(int p=0;p<n_pop;p++){
        agewise_populations[pop[p][a]]++;
        if(pop[p][s]==S && !is_vaccinated[p]){some_protection++; agewise_some_protection[pop[p][a]]++;}
      }

      double bs = n_pop - some_protection;
      double agewise_bs[max_age] = {};

      for(int age=0;age<max_age;age++){
        agewise_bs[age] = agewise_populations[age] - agewise_some_protection[age];
      }

      //***************************** Write the output to an array *****************************//

      output[day][0] = day; for(int q=0;q<n_states;q++){output[day][q+1]=n[q];}
      output[day][n_states+1] = vaccines_administered;
      output[day][n_states+2] = bs;

      age_output[day][0] = day; int counter = 1; for(int r=0;r<max_age;r++){for(int q=0;q<n_states;q++){age_output[day][counter]=n_per_age[r][q];counter++;} age_output[day][counter] = agewise_vaccines_administered[r]; counter++;
      age_output[day][counter] = agewise_bs[r]; counter++;}

      //*****************************************************************************************//

    }

    /********************* END OF END OF DAY(S) ********************************/

    t+=dt;
  }  // End of While Loop above.

  //*** Checking fraction of HCW who contracted the disease during the pandemic ***//
  int hcw_recovered = 0;
  int hcw = 0;
  for(int i=0;i<n_pop;i++){               // Find all HCW, and mark those that have recovered or died.
    if(pop[i][w]<n_hospitals){
      hcw++;
      if(pop[i][s]==R || pop[i][s]==D){hcw_recovered++;}
    }
  }

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

  int details[] = {locations_moved, hcw_recovered, hcw};
  writetofile(output, age_output, tf, cpu_time_used, details, iter);


}// End of run function.



int main() {

  char filename[1000]="synthetic_population.csv";

  /************ PARAMETERS ************/

  double lambda_S = 0.7;
  double lambda_E = 1/4.5;
  double lambda_A = 1.0/8;
  double lambda_P = 1/1.1;
  double lambda_MI = 1.0/8;
  double lambda_SI = 1/1.5;
  double lambda_H  = 1/18.1;


  /************* AGE FACTOR *************/

  double lambda_S_array[] = {1.0*lambda_S, 1.0*lambda_S, 1.0*lambda_S, 1.0*lambda_S, 1.0*lambda_S, 1.0*lambda_S, 1.0*lambda_S, 1.0*lambda_S, 1.0*lambda_S, 1.0*lambda_S};

  double gamma_array[] = {1-0.5, 1-0.55,1- 0.6, 1-0.65,1-0.7, 1-0.75,1-0.8, 1-0.85,1- 0.9, 1-0.9};

  double delta_array[] = {1-0.0005,  1-0.00165,  1-0.00720,  1-0.02080,  1-0.03430,  1-0.07650,  1-0.13280,  1-0.20655,  1-0.24570,  1-0.24570};

  double sigma_array[] = {0.00002, 0.00002, 0.0001, 0.00032, 0.00098, 0.00265, 0.00766, 0.02439, 0.08292, 0.16190};


  double lambda_E_array[]  = {lambda_E , lambda_E  , lambda_E , lambda_E ,  lambda_E  , lambda_E  , lambda_E, lambda_E  , lambda_E  , lambda_E};
  double lambda_A_array[]  = {lambda_A , lambda_A  , lambda_A , lambda_A ,  lambda_A  , lambda_A  , lambda_A, lambda_A  , lambda_A  , lambda_A};
  double lambda_P_array[]  = {lambda_P , lambda_P  , lambda_P , lambda_P ,  lambda_P  , lambda_P  , lambda_P, lambda_P  , lambda_P  , lambda_P};

  double lambda_MI_array[] = {lambda_MI, lambda_MI , lambda_MI , lambda_MI  , lambda_MI , lambda_MI , lambda_MI , lambda_MI , lambda_MI , lambda_MI};
  double lambda_SI_array[] = {lambda_SI, lambda_SI , lambda_SI , lambda_SI  , lambda_SI , lambda_SI , lambda_SI , lambda_SI , lambda_SI , lambda_SI};

  double lambda_H_array[]  = {lambda_H , lambda_H , lambda_H , lambda_H , lambda_H , lambda_H , lambda_H , lambda_H , lambda_H , lambda_H};


  /************* RATE ARRAY *************/

  for(int i=0;i<max_age;i++){
    rate[i][S][E] = lambda_S_array[i];                       // S -> E
    rate[i][E][A] = gamma_array[i]*lambda_E_array[i];        // E -> A
    rate[i][E][P] = (1 - gamma_array[i])*lambda_E_array[i];  // E -> P
    rate[i][A][R] = lambda_A_array[i];                       // A -> R
    rate[i][P][MI]= delta_array[i]*lambda_P_array[i];        // P -> MI
    rate[i][P][SI]= (1-delta_array[i])*lambda_P_array[i];    // P -> SI
    rate[i][MI][R]= lambda_MI_array[i];                      // MI -> R
    rate[i][SI][R]= 0;                                       // SI -> R (Nobody recovers directly from SI)
    rate[i][SI][H]= lambda_SI_array[i];                      // SI -> H
    rate[i][H][D] = sigma_array[i]*lambda_H_array[i];        // H -> R
    rate[i][H][R] = (1-sigma_array[i])*lambda_H_array[i];    // H -> R
  }


  /******* SIMULATION PARAMETERS *******/

  int tf = 200;                            // Total simulation run time in days


  //********************* Monte Carlo run with the above specifications *********************//

  int mc_runs = 1;

  double unlockSchoolsOnArray[] = {0,10,20,30,40,50,60,70,80,90,100};
  int n_unlocks = sizeof(unlockSchoolsOnArray)/sizeof(unlockSchoolsOnArray[0]);

  double irf[] = {30.0/100.0, 50.0/100.0};
  int n_irf    = sizeof(irf)/sizeof(irf[0]);

  for(int i=0;i<mc_runs;i++){
    dvr = 0.2;                                      // Daily vaccination rate (in percentage of total population per day)
    vaccination_strategy=descending;                // First vaccinate 40+ for a month, then 20+

    for(int r=0;r<n_irf;r++){
      initial_recovered_fraction = irf[r];          // Initial fraction of recovered

      for(int j=0;j<n_unlocks;j++){
        unlockSchoolsOn = unlockSchoolsOnArray[j];  // Unlock schools on day

        int start_vacc = 0.0;                       // Start vaccination immediately

        createPopulation(filename);                 // Create a random population with a fixed infection seed (default: 1%) and recovered fraction (default: 30%), set by initial_asymptomatic_fraction
                                                    // Resets the pop array, n_per_location array, and people_linked_to array.

        run(tf,
            start_vacc,
            dvr,
            unlockSchoolsOn,
            i);
       }
    }

  }

 }

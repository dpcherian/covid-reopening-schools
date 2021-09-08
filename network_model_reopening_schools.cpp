//Including age stratification and schools
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

int poisson(double mu)  // Generate a Poisson Random Number with mean mu
{
    std::poisson_distribution<> d(mu);
    return d(gen);
}

double uniform()        // Generate a uniform random number between [0,1)
{
    std::uniform_real_distribution<> d(0.0,1.0);
    return d(gen);
}

int randint(int a,int b) // Generate a uniform random integer in [a, b)
{
    std::uniform_int_distribution<> d(a,b-1);
    return d(gen);
}


int randint(int b) // Generate a uniform random number in [0, b)
{
    std::uniform_int_distribution<> d(0,b-1);
    return d(gen);
}

/*****************************************************************/


/* Global variables *******************************************/

double dt = 0.1;
char output_folder[100] = "data";

const int mult = 1;              // Multiplier to scale up problem (1: n_pop = 1 x 10,000, 2: n_pop = 20,000, etc.)

const int n_pop = mult*20316;    // Total population
const int n_loc = mult*6623;     // Total number of locations
const int n_net = mult*123;      // Total number of networks
const int n_overlap = n_net+1;   // Start position of home locations + 1
const int n_hospitals = mult*1;  // Set the number of hospitals.

const int start_schools = n_hospitals; // Where school locations start
const int n_schools = mult*1;   // Set the number of schools (10)

const int classrooms = 100;       // Number of classrooms per school

int n_asym = 0;//mult*100;            // Initial number of asymptomatics
int n_rec  = 0;//mult*3000;            // Initial number of asymptomatics

double initial_asymptomatic_fraction = 1.00/100;
double initial_recovered_fraction    = 30.0/100;
double initial_vaccinated_fraction   = 20.0/100;

const int maxppl = 10000;        // Maximum number of people who can be linked to a location (to avoid 2GB limit on arrays) [Not needed anymore]

char states[][3] = {"S","E", "A","P", "MI","SI","R","H", "D"};
                  // 0   1    2    3   4    5    6   7    8

const int S = 0;
const int E = 1;
const int A = 2;
const int P = 3;
const int MI= 4;
const int SI= 5;
const int R = 6;
const int H = 7;
const int D = 8;

const int s = 0;   // Index of "state" attribute

const int a = 1;   // Index of "age" attribute
const int x = 2;   // Index of "some factor" attribute


const int h = 3;   // Index of "home loc" attribute
const int w = 4;   // Index of "work loc" attribute
const int c = 5;   // Index of "current loc" attribute

const int n_states = sizeof(states)/sizeof(states[0]); // Number of states

const int person_attr = 6;     // Number of attributes in the pop array (currently: state, age, "something", home location, work location, currently location)
const int max_age     = 9;

double Cpars[] = {1,  1,   1,    1,   0.1,  0.1}; // Contact parameters for transmission
               // SA  SP  S-MI  S-SI  S-H   S-Q

/* Initial populations ***************************/

int n[n_states];// {n_pop-n_asym, n_asym, 0, 0,  0,  0, 0 };
                     // [S,        A,    P, MI, SI, R, H ]
 int n_per_age[max_age][n_states];//   // NEW, for age-stratification


/*************************************************/

int total_loc_confined_time = 14;
int total_isolation_time    = 14;
int if_positive_test_after  = 14;

int days_bw_hq_tests = 0;
int days_bw_lq_tests = 0;

double rate[max_age][n_states][n_states] = {};
// double age_factor[max_age][n_states] = {};

int pop[n_pop][person_attr] = { };
int n_per_location[n_loc][n_states];

vector<vector<int>> people_linked_to(n_loc);
int len_ppl_linked_to[n_loc] = {};

vector<vector<vector<int>>> n_per_room(n_loc);
bool is_school[n_loc] = {};
int n_rooms[n_loc] = {};


// New vaccination arrays
bool is_vaccinated[n_pop] = {};
int vaccinated_on[n_pop] = {};
int vaccs_per_state[n_states] = {}; // NEW: Remove quick?
int vaccines_administered = 0;
int agewise_vaccines_administered[max_age] = {};

int current_room[n_pop] = {};
int assigned_room[n_pop]={};
bool is_adult[n_pop]={};

// New variables for vaccinated //

double infection_reduction       = 1.0 - 0.441;
double gamma_fractional_increase = 1.0 - 0.598;
double transmission_reduction    = 0.4;

int vaccination_strategy = 0;
int ascending = 1;
int descending=-1;
// int vaccines_administered = 0;

double unlockSchoolsAt = 100;
double dvr = 0.0;

//*****************************//


bool quarantine_confined = true;                   // Change to false if you don't want to confine individuals who've tested positive
bool lock_schools = true;                          // NEW: Don't allow people to go to schools. Default: true (start schools locked)
bool shuffle_teachers = false;                     // NEW: Variable that decides whether teachers move between classrooms. Default: false (don't let them move.)

const int op_width = 1 + n_states + 3 + 2 + 1 + 1 + 1;     // Width of output array: <time>(1) <Number of states> <Testing details>(3) <Home quarantine details>(3)
const int age_op_width = 1 + (n_states+2)*max_age;         // Width of age-wise output array <time>(1) (<Number of states> <Vaccines administered> <Background Seropositivity> ) x max_age
int positives = 0;
/*****************************************************************/

/* Other functions ***********************************************/

int sum(int array[],int len){

  int sum = 0;

  for(int i =0; i<len;i++){
    sum += array[i];
  }

  return sum;
}

int sum(bool array[],int len){

  int sum = 0;

  for(int i =0; i<len;i++){
    sum += array[i];
  }

  return sum;
}


double sum(double array[],int len){

  double sum = 0;

  for(int i =0; i<len;i++){
    sum += array[i];
  }

  return sum;
}

int sum(int array[][n_states], const int len, int axis){

  int sum = 0;

  for(int i =0; i<len;i++){
    sum += array[i][axis];
  }

  return sum;
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


void readcsv(char *filename)
{
int num_to_prevaccinate = (int)(initial_vaccinated_fraction*n_pop);
// printf("Number to prevacc at start %i\n",num_to_prevaccinate );

FILE* fp = fopen(filename, "rt");
if (fp==NULL) {printf("no such file\n");exit(0);}

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

      if(column==s){                                              // Setting initial states
        double r = uniform();
        if(r<initial_asymptomatic_fraction){pop[row][column] = A; n_asym++;}
        else if(r>initial_asymptomatic_fraction && r<initial_asymptomatic_fraction+initial_recovered_fraction){ pop[row][column]=R; n_rec++;}
        else{pop[row][column] = S;}
      }

		column++;
		}

    is_vaccinated[row] = (num_to_prevaccinate>0 && pop[row][a] >= 6 && uniform()<0.8) ? true : false;
    vaccinated_on[row] = is_vaccinated[row] ? -1000 : 1000;
    if(is_vaccinated[row]){
      num_to_prevaccinate--;
      vaccs_per_state[pop[row][s]]++;
      agewise_vaccines_administered[pop[row][a]]++;
      vaccines_administered++;
    }

    n_per_location[pop[row][c]][pop[row][s]] += 1;  // Increment number of disease-state in current location
    n_per_age     [pop[row][a]][pop[row][s]] += 1;  // Increment number of disease state per age

    if(pop[row][a]>1){is_adult[row] = true;}

    int work_loc    = pop[row][w];
    int current_loc = pop[row][c];

    assigned_room[row] = randint(n_rooms[work_loc]);
    current_room[row]  = 0;
    // printf("loc %i rooms %i assigned %i\n",work_loc, n_rooms[work_loc], assigned_room[row] );

    n_per_room[current_loc][current_room[row]][pop[row][s]] += 1; // Increase number of current state in current room of current location

	row++;

	free(tofree); // I don't fully understand why I need to keep freeing buffer and reallocating it, and why it needs to be done like this, but anything else doesn't seem to work!
	}

  //*********** Prevaccinate randomly if any doses are left ********//
  // printf("Num to prevacc before random %i\n", num_to_prevaccinate);
  while(num_to_prevaccinate>0){
    int p = randint(n_pop);
    if(!is_vaccinated[p]){
      is_vaccinated[p] = is_adult[p] ? true : false;
      vaccinated_on[p] = is_vaccinated[p] ? -1000 : 1000;
      if(is_vaccinated[p]){
        num_to_prevaccinate--;
        vaccs_per_state[pop[p][s]]++;
        agewise_vaccines_administered[pop[p][a]]++;
        vaccines_administered++;
      }
    }
  }
  // printf("Num to prevacc after random %i, vaccs administered %i\n", num_to_prevaccinate,vaccines_administered);
  //************ Done prevaccinating *********************************//

fclose(fp);

//***** CREATE PEOPLE-LINKED-TO ARRAY ***************************************************//

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

//***** DONE CREATING PEOPLE-LINKED-TO ARRAY ********************************************//



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
  // Create a people_linked_to list of all people associated with a location (either as work or home).
  // Reset all classroom numbers

  for(int i=0; i<n_loc;i++){

    // Reset n_per_location
    for(int j=0;j<n_states;j++){n_per_location[i][j]=0;}

    // Deciding if a location is a school or not
    if(i<start_schools+n_schools && i>=start_schools){
      is_school[i] = true;
    }
    else{ is_school[i] = false; }

    // Assigning quantities per room:

    if(!is_school[i]){ n_rooms[i] = 1;} else{ n_rooms[i] = classrooms;}  // If it's not a school, assign only 1 room, else assign `classroom` rooms.

    n_per_room[i] = vector<vector<int>>(n_rooms[i]);
    for(int room=0; room<n_rooms[i];room++){        // For each room in locations
      n_per_room[i][room] = vector<int>(n_states,0);// assign a vector which -- for each room -- contains the number of individuals per state
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

  // Todo: Read from CSV here

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
void writetofile(int output[][op_width], int age_output[][age_op_width], int tf, double Tpars[2][4], double begin_at, double test_frac, double time_taken, int details[9],int iter){
  // FOR REFERENCE: int details[] = {quarantine_confined, lock_homes, quarantine_when_sample_taken, lq_tests_conducted, hq_tests_conducted, tests_conducted, results_declared, locations_moved, hcw_recovered};
  // Write output to a FILE
  errno = 0;      // Variable to store error number in case file open does not work

  // sprintf(outputFilename,"./%s/output/Targeted_BeginAt_%g_DTR_%g_RAT_%g_%g_PCR_%g_%g_%lf%lf-%i.txt",names[iter],begin_at,test_frac,Tpars[0][0],Tpars[0][3],Tpars[1][0],Tpars[1][3],uniform(),uniform(),iter);
  double rn1 = uniform();
  double rn2 = uniform();

  // printf("./%s/VR_%g_%s_UnlockSchoolsAt_%g_%lf%lf-%i.txt",output_folder,dvr,vaccination_strategy==ascending ? "Ascending" : vaccination_strategy == descending ? "Descending" : "Unset",unlockSchoolsAt,rn1,rn2,iter);

  sprintf(outputFilename,"./%s/Total_IR_%g_IV_%g_VR_%g_%s_UnlockSchoolsAt_%g_%lf%lf-%i.txt",output_folder,initial_recovered_fraction*100, initial_vaccinated_fraction*100, dvr,vaccination_strategy==ascending ? "Ascending" : vaccination_strategy == descending ? "Descending" : "Unset",unlockSchoolsAt,rn1,rn2,iter);
  FILE *fpt=(FILE *)fopen(outputFilename,"wt"); // Open the file to print output

  if(fpt){

    fprintf(fpt,"###### TEST LOG ####################\n");
    fprintf(fpt,"# Time taken               : %.2f s\n",time_taken);

    fprintf(fpt,"# Test Parameters: \n");
    fprintf(fpt,"# %.2f %.2f %.2f %.2f\n", Tpars[0][0], Tpars[0][1], Tpars[0][2], Tpars[0][3]);
    fprintf(fpt,"# %.2f %.2f %.2f %.2f\n", Tpars[1][0], Tpars[1][1], Tpars[1][2], Tpars[1][3]);

    fprintf(fpt,"# Confined Less Infective? : %s\n", details[0] ? "True": "False");
    fprintf(fpt,"# Homes Quarantined?       : %s\n", details[1] ? "True": "False");
    fprintf(fpt,"# Quarantine when sampled? : %s\n", details[2] ? "True": "False");
    fprintf(fpt,"# Testing Started When     : %.2f%% recovered\n", begin_at);
    fprintf(fpt,"# Fraction Tested Daily    : %.2f%%\n", test_frac);

    fprintf(fpt,"# LQ Tests Done in total   : %d\n", details[3]);
    fprintf(fpt,"# HQ Tests Done in total   : %d\n", details[4]);
    fprintf(fpt,"# All Tests Done in total  : %d\n", details[5]);
    fprintf(fpt,"# Results Given in total   : %d\n", details[6]);
    fprintf(fpt,"# Locations Moved in total : %d\n", details[7]);
    fprintf(fpt,"# Total recovered HCW      : %d\n", details[8]);
    fprintf(fpt,"# Total HCW                : %d\n", details[9]);

    fprintf(fpt,"# Vaccination strategy     : %s\n", vaccination_strategy==ascending ? "Ascending" : vaccination_strategy == descending ? "Descending" : "Unset");
    fprintf(fpt,"# Total Vaccines Given     : %d\n", vaccines_administered);
    fprintf(fpt,"# Schools unlocked on day  : %g\n", unlockSchoolsAt);

    fprintf(fpt,"# Rate Array: \n");
    for(int i=0; i<n_states;i++){fprintf(fpt,"# ");for(int j=0; j<n_states;j++){fprintf(fpt,"%5g ", rate[0][i][j]);}fprintf(fpt,"\n");}// NOTE! : Change to correct rates
    fprintf(fpt,"###### END LOG #####################\n");
    fprintf(fpt,"#\n");
    fprintf(fpt, "#%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n","Day","nS","nE","nA","nP","nMI","nSI","nR","nH","nD","nV","nBS","PCR_conducted_today", "RAT_conducted_today", "Tests_remaining_today", "Agents_currently_confined", "Quarantines_removed_today","Locations_in_quarantine_today" );

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

  // sprintf(ageOutputFilename,"./%s/output/AgeStratified_BeginAt_%g_DTR_%g_RAT_%g_%g_PCR_%g_%g_%lf%lf-%i.txt",names[iter],begin_at,test_frac,Tpars[0][0],Tpars[0][3],Tpars[1][0],Tpars[1][3],uniform(),uniform(),iter);
  sprintf(ageOutputFilename,"./%s/AgeStratified_IR_%g_IV_%g_VR_%g_%s_UnlockSchoolsAt_%g_%lf%lf-%i.txt",output_folder,initial_recovered_fraction*100, initial_vaccinated_fraction*100, dvr,vaccination_strategy==ascending ? "Ascending" : vaccination_strategy == descending ? "Descending" : "Unset",unlockSchoolsAt,rn1,rn2,iter);
  FILE *fpt1=(FILE *)fopen(ageOutputFilename,"wt"); // Open the file to print output

  if(fpt){

    fprintf(fpt1,"###### TEST LOG ####################\n");
    fprintf(fpt1,"# Time taken               : %.2f s\n",time_taken);
    fprintf(fpt1,"# Main file                : %s \n",outputFilename);

    fprintf(fpt1,"# Test Parameters: \n");
    fprintf(fpt1,"# %.2f %.2f %.2f %.2f\n", Tpars[0][0], Tpars[0][1], Tpars[0][2], Tpars[0][3]);
    fprintf(fpt1,"# %.2f %.2f %.2f %.2f\n", Tpars[1][0], Tpars[1][1], Tpars[1][2], Tpars[1][3]);

    fprintf(fpt1,"# Confined Less Infective? : %s\n", details[0] ? "True": "False");
    fprintf(fpt1,"# Homes Quarantined?       : %s\n", details[1] ? "True": "False");
    fprintf(fpt1,"# Quarantine when sampled? : %s\n", details[2] ? "True": "False");
    fprintf(fpt1,"# Testing Started When     : %.2f%% recovered\n", begin_at);
    fprintf(fpt1,"# Fraction Tested Daily    : %.2f%%\n", test_frac);

    fprintf(fpt1,"# LQ Tests Done in total   : %d\n", details[3]);
    fprintf(fpt1,"# HQ Tests Done in total   : %d\n", details[4]);
    fprintf(fpt1,"# All Tests Done in total  : %d\n", details[5]);
    fprintf(fpt1,"# Results Given in total   : %d\n", details[6]);
    fprintf(fpt1,"# Locations Moved in total : %d\n", details[7]);
    fprintf(fpt1,"# Total recovered HCW      : %d\n", details[8]);
    fprintf(fpt1,"# Total HCW                : %d\n", details[9]);

    fprintf(fpt,"# Vaccination strategy     : %s\n", vaccination_strategy==ascending ? "Ascending" : vaccination_strategy == descending ? "Descending" : "Unset");
    fprintf(fpt,"# Total Vaccines Given     : %d\n", vaccines_administered);
    fprintf(fpt,"# Schools unlocked on day  : %g\n", unlockSchoolsAt);

    fprintf(fpt1,"# Rate Array: \n");
    for(int i=0; i<n_states;i++){fprintf(fpt1,"# ");for(int j=0; j<n_states;j++){fprintf(fpt1,"%5g ", rate[0][i][j]);}fprintf(fpt1,"\n");}// NOTE! : Change to correct rates
    fprintf(fpt1,"###### END LOG #####################\n");
    fprintf(fpt1,"#\n");
    // fprintf(fpt, "# %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n","Day","nS","nA","nP","nMI","nSI","nR","nH","PCR_conducted_today", "RAT_conducted_today", "Tests_remaining_today", "Agents_currently_confined", "Quarantines_removed_today","Locations_in_quarantine_today" );
    fprintf(fpt1, "#Day ");
    for(int j=0;j<max_age;j++){
      for(int i=0;i<n_states;i++){
        fprintf(fpt1, "%s%i ",states[i],(j+1)*10);
      }
      fprintf(fpt1, "%s%i ","V",(j+1)*10);
      fprintf(fpt1, "%s%i ","BS",(j+1)*10);
    }
    fprintf(fpt1,"\n");

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


void Targeted_Run(double Tpars[][4], int tf, bool lock_homes, bool quarantine_when_sample_taken, double begin_at, double test_frac, double start_vacc, double dvr, double unlockSchoolsAt, int iter){

  clock_t start, end; // Measuring how long the function takes to run
  double cpu_time_used;
  start = clock();

  bool is_confined[n_pop]  = {};
  bool being_tested[n_pop] = {};
  bool loc_confined[n_loc] = {};

  // New school array
  // bool is_school[n_loc] = {};
  // for(int i=start_schools;i<start_schools+n_schools;i++){is_school[i]=true;}

  //* Added to keep track of retesting RAT negative symptomatics **//
  int last_test_type[n_pop]={};
  std::fill_n(last_test_type,n_pop,-1000);  // Set all values of last_test_type to -1000
  int last_test_result[n_pop]={};
  std::fill_n(last_test_result,n_pop,-1000);  // Set all values of last_test_result to -1000
  /*****************************************************************/

  int test_result[n_pop]   = {};
  int next_test_date[n_pop]= {};

  int result_declared_date[n_pop]={};
  std::fill_n(result_declared_date,n_pop,-1000);  // Set all values of result_declared_date to -1000

  int loc_confined_time[n_loc] = {};
  std::fill_n(loc_confined_time,n_loc,-1000);     // Set all values of loc_confined_time to -1000

  int person_isolated_time[n_pop] = {};
  std::fill_n(person_isolated_time, n_pop,-1000);

  //   HQ=>PCR     LQ=>RAT
  int hq_tests_conducted = 0;
  int lq_tests_conducted = 0;
  int tests_conducted    = 0;
  int results_declared   = 0;
  int locations_moved    = 0;

  double r[n_states][n_states] = {};                      // Array to store rates per event

  double alphaH = 1 - Cpars[4];                           // Quantity by which V is reduced per hospitalised individual
  double alphaQ = 1 - Cpars[5];

  double alphaT = 1 - transmission_reduction;             // NEW: effective reduction in number of people in location due to vaccination

  int tests_available_daily = test_frac/100 * n_pop;
  int tests_remaining_today = test_frac/100 * n_pop;

  int lq_tests_daily = Tpars[0][3]*tests_available_daily;
  int lq_tests_today = Tpars[0][3]*tests_available_daily;

  int hq_tests_daily = Tpars[1][3]*tests_available_daily;
  int hq_tests_today = Tpars[1][3]*tests_available_daily;

  double lq_sens = Tpars[0][0];
  double lq_spec = Tpars[0][1];
  double lq_delay= Tpars[0][2];

  double hq_sens = Tpars[1][0];
  double hq_spec = Tpars[1][1];
  double hq_delay= Tpars[1][2];

  double t = 0.0;
  int day  = 0;
  bool midday_move_completed = false;

  bool first_move_done =false;
  bool second_move_done=false;
  bool third_move_done =false;

  bool lock_schools = true; //NEW: Start with schools always locked!

  int output[tf+1][op_width];   // Output array, to be printed to file
  int age_output[tf+1][age_op_width]; // NEW: to print out output per age


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
  output[day][n_states+3]=hq_tests_conducted;output[day][n_states+4]=lq_tests_conducted;output[day][n_states+5]=tests_remaining_today;output[day][n_states+6]=0;output[day][n_states+7]=0;output[day][n_states+8]=0;

  age_output[day][0] = day; int counter = 1; for(int r=0;r<max_age;r++){for(int q=0;q<n_states;q++){age_output[day][counter]=n_per_age[r][q];counter++;} age_output[day][counter]=agewise_vaccines_administered[r];counter++;
  age_output[day][counter]=agewise_bs[r];counter++;}

  double exit_rate[n_states] = {};

  while(t<tf){

    if(t>unlockSchoolsAt && lock_schools==true){  // Check if schools need to be unlocked
      lock_schools=false;
      printf("Schools unlocked at t=%.2f\n",t);
    }


    // Moving people around deterministically (HOME TO WORK)
    if(midday_move_completed == false && t - day >  0.5){
      midday_move_completed = true;

      for(int i=0;i<n_pop;i++){

        if(pop[i][s]!=H && pop[i][s]!=D && is_confined[i] == false && loc_confined[pop[i][c]]==false){
          locations_moved++;
          int home_loc = pop[i][h];
          int work_loc = pop[i][w];
          if(pop[i][c]==home_loc)if(lock_schools==false || ( lock_schools==true && !is_school[work_loc])){ pop[i][c] = work_loc; n_per_location[home_loc][pop[i][s]]--; n_per_location[work_loc][pop[i][s]]++; n_per_room[home_loc][current_room[i]][pop[i][s]]--; current_room[i]=assigned_room[i]; n_per_room[work_loc][assigned_room[i]][pop[i][s]]++; }
        }
      }
    }

    if(lock_schools==false && ((first_move_done==false && t-day>0.5 && t-day<0.5+0.17) || (second_move_done==false && t-day>0.5+0.17 && t-day<0.5+0.33) || (third_move_done==false && t-day>0.5+0.33 && t-day<0.5+0.5)) ){
      // Move them within the school classrooms

      for(int i=start_schools; i<start_schools+n_schools; i++){ // For every school location
        // int school = i-start_schools;

        // //******** Print out initial populations ********//
        // if(first_move_done==false){
        //   printf("Start\n");
        //   for(int j=0;j<n_rooms[i]; j++){
        //     for(int k=0; k<n_states;k++){
        //       printf("%i ",n_per_room[i][j][k] );
        //     }
        //     printf("\n");
        //   }
        //   printf("End\n");
        // }
        // // **************************************************//

        int num_students_at_school = 0;

        for(int j=0; j<len_ppl_linked_to[i];j++){   // Loop over all people in this school
          int person = people_linked_to[i][j];
          int person_state = pop[person][s];
          if(pop[person][c] == i){                  // If the person is currently in the school, add them to a list.
            num_students_at_school++;

            if(first_move_done==false){                   // If they're just entering the school
              // int random_classroom = randint(n_rooms[i]); // choose random classroom from all possible rooms in this location
              int set_classroom = assigned_room[person];
              if(set_classroom != current_room[person]){printf("Set classroom is %i current room is %i\n",set_classroom,current_room[person] );}
              current_room[person] = set_classroom; // Set the person to their assigned classroom.
              // n_per_room[i][set_classroom][person_state]++;
            }
            else if(is_adult[person] && shuffle_teachers){
              int from_room = current_room[person];  // Get current room of teacher in school i
              int to_room = randint( n_rooms[i] ); // Choose random room from the number of rooms in this school

              current_room[ person ] = to_room;

              n_per_room[i][from_room][person_state]--;
              n_per_room[i][to_room][person_state]++;
            }

            if(pop[person][s] == H){printf("ERROR! Hospitalised in School!\n"); exit(12);}
          }
        }

        int temp_tot = 0;
        // int old_tot = 0;

        for(int classroom=0; classroom<n_rooms[i]; classroom++){
          for(int state = 0; state<n_states; state++){
              temp_tot += n_per_room[i][classroom][state];
          }
          // printf("Temp Total for classroom %i = %i\n",classroom,  temp_tot-old_tot);
          // old_tot=temp_tot;
        }
      //   if(i == start_schools){
      //   //********************** Print out classrooms **********************//
      //   printf("School %i Students %i Temp_tot %i Rooms %i \n",i,num_students_at_school, temp_tot, n_rooms[i] );
      //   for(int j=0;j<n_rooms[i];j++){
      //     for(int s=S;s<n_states;s++){printf("%i ",n_per_room[i][j][s] );}
      //     printf("\n");
      //   }
      //   printf("***********************************\n");
      //   // exit(11);
      //   //*****************************************************************//
      // }

        if(temp_tot!=num_students_at_school){printf("Error! Sum of n_per_classroom %4i, num students %4i \n",temp_tot, num_students_at_school ); exit(2);}

      }// End loop over schools

      // if(first_move_done==false){first_move_done=true; printf("Day %i  First Move Done\n",day);}
      // else if(second_move_done==false){second_move_done=true;printf("Day %i Second Move Done\n",day);}
      // else if(third_move_done==false){third_move_done=true;printf("Day %i  Third Move Done\n",day);}
      // else{printf("Error, all three moves are done, but movement is still occurring.\n");}
      if(first_move_done==false){first_move_done=true; }
      else if(second_move_done==false){second_move_done=true;}
      else if(third_move_done==false){third_move_done=true;}
      else{printf("Error, all three moves are done, but movement is still occurring.\n");}

    }

    for(int i=0; i<n_loc; i++){
      // printf("Location %4i Rooms %4i\n",i,n_rooms[i]);
      for(int room=0; room<n_rooms[i];room++){

        int N = 0;

        int counter = 0;

        vector<int> ind;

        for(int state = 0; state<n_states;state++){N += n_per_room[i][room][state];} // Total number currently in this room

        int conf_by_state_in_loc[n_states]={};

        for(int p=0;p<len_ppl_linked_to[i];p++){ // For all people linked to this school
          int person = people_linked_to[i][p];

          if(pop[person][c]==i && current_room[person] == room){ ind.push_back(person); counter++;}  // If the person is in this school and is currently in this room, add them to ind
          else if(pop[person][c]==i && is_confined[person]){ conf_by_state_in_loc[pop[person][s]]++;}
        }

        int sum=0;for(int qq=0;qq<n_states;qq++){sum+=n_per_room[i][room][qq];}

        if(counter!=N){printf("Error in location %i, N = %i, sum %i, counter = %i at time %lf\n",i,N,sum,counter,t );}

        if(N==0){continue;}

        int conf_by_state_in_room[n_states]={};
        int tot_conf = 0;

        int vacc_by_state_in_room[n_states]={}; // New: To reduce transmission of vaccinated individuals
        int tot_vacc = 0;

        for(int j = 0; j<N;j++){      // Go over people in this room
            int p = ind[j];           // index of person
            if(quarantine_confined == true && is_confined[p] == true) {conf_by_state_in_room[pop[p][s]]++; tot_conf++; } // If the person is confined, increment the conf_by_state_in_loc of their state
            else if(is_vaccinated[p]==true){vacc_by_state_in_room[pop[p][s]]++; tot_vacc++; }
        }


        if(loc_confined[i]==true && conf_by_state_in_loc[S]+conf_by_state_in_loc[E]+conf_by_state_in_loc[A]+conf_by_state_in_loc[P]+conf_by_state_in_loc[MI]+conf_by_state_in_loc[SI]+conf_by_state_in_loc[R]+conf_by_state_in_loc[H] == 0){loc_confined[i]=false;loc_confined_time[i]=-1000;} // Unlock the house if there are no confined people

        shuffle(0,N,ind); // Shuffle list of people currently in locations

        for(int j=0;j<N;j++){                    // Loop over the people currently in the location

          int age = pop[ind[j]][a];

          // NEW : ADDING THE VACCINE MULTIPLIER:
          double f = rate[age][E][A]/(rate[age][E][A]+rate[age][E][P]);            // default vals: If the person isn't a vaccinated test subject, f is just gamma (beta not reduced)
          double lambda_E = rate[age][E][A]+rate[age][E][P];
          double beta_multiplier = 1.0;                                            // default beta-multipler (reduction) = 1 unless vaccinated

          if(is_vaccinated[ind[j]]){

            int vday   = vaccinated_on[ind[j]];

            double this_gamma = rate[age][E][A]/(rate[age][E][A]+rate[age][E][P]);

            //      if(test_subjects[tindex][2] == 1){ f = std::min(this_gamma + (0.8*(1 - this_gamma)/42)*(t - vday),this_gamma + (0.8*(1 - this_gamma)/42)*28);} // first vaccination didn't work perfectly
            // else if(test_subjects[tindex][2] == 2){ f = std::min(this_gamma + (0.8*(1 - this_gamma)/42)*(t - vday),this_gamma + (0.8*(1 - this_gamma)/42)*42);}      // second vaccination worked perfectly
            int fd = 90+14;//28;
            // beta_multiplier = infection_reduction;
            beta_multiplier = std::min(std::max(1 + (infection_reduction - 1)*(t - vday)/fd, 1 + (infection_reduction - 1)), 1.0);

            f = std::min(std::max(this_gamma + (1 - (gamma_fractional_increase/infection_reduction)*(1 - this_gamma) - this_gamma)*(t - vday)/fd, this_gamma), this_gamma + (1 - (gamma_fractional_increase/infection_reduction)*(1 - this_gamma) - this_gamma));


            // f = std::min(this_gamma + (gamma_fractional_increase*(1 - this_gamma)/fd)*(t - vday),this_gamma + (gamma_fractional_increase*(1 - this_gamma)/fd)*fd); // the multipler changes from gamma to this.

            //****** To print out gamma for a single vaccinated individual ******//
            // if(ind[j]==31){printf("Vaccinated person %i age %i gamma %lf Vaccinated on %i Today %i \n",ind[j],age,f,vaccinated_on[ind[j]], day);}
            // if(ind[j]==31){printf("%lf %lf \n",t,f);}
            // **************************************************************** //
          }


          int newN = n_per_room[i][room][S]+n_per_room[i][room][E]+n_per_room[i][room][A]+n_per_room[i][room][P]+n_per_room[i][room][MI]+n_per_room[i][room][SI]+n_per_room[i][room][R]+n_per_room[i][room][H]+n_per_room[i][room][D];

          double V = newN - alphaH*n_per_room[i][room][H];   // Spatial damping parameter (adjusted by alphaH for hospitals)

          // S->E
          r[S][E] = beta_multiplier * rate[age][S][E] * 1/V * (Cpars[0]*(n_per_room[i][room][A] - conf_by_state_in_room[A]*alphaQ - vacc_by_state_in_room[A]*alphaT) +    // SA
                                                               Cpars[1]*(n_per_room[i][room][P] - conf_by_state_in_room[P]*alphaQ - vacc_by_state_in_room[P]*alphaT) +    // SP
                                                               Cpars[2]*(n_per_room[i][room][MI]-conf_by_state_in_room[MI]*alphaQ - vacc_by_state_in_room[MI]*alphaT) +   // S-MI
                                                               Cpars[3]*(n_per_room[i][room][SI]-conf_by_state_in_room[SI]*alphaQ - vacc_by_state_in_room[SI]*alphaT) +   // S-SI
                                                               Cpars[4]*n_per_room[i][room][H]);                                                                          // SH
          // E->A
          r[E][A] =  f * lambda_E;

          // E->P
          r[E][P] = (1 - f) * lambda_E;

          // A->R
          r[A][R] = rate[age][A][R];
          // P->MI
          r[P][MI] = rate[age][P][MI];
          // P->SI
          r[P][SI] = rate[age][P][SI];
          // MI->R
          r[MI][R] = rate[age][MI][R];
          // SI->R
          r[SI][R] = rate[age][SI][R];
          // SI->H
          r[SI][H] = rate[age][SI][H];
          // H->R
          r[H][R] = rate[age][H][R];
          // H->D
          r[H][D] = rate[age][H][D];

          // printf("S->E %lf\n",r[S][E] );
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
            double d =  (Cpars[0]*(n_per_room[i][room][A] - conf_by_state_in_room[A]*alphaQ - vacc_by_state_in_room[A]*alphaT) +    // SA
                         Cpars[1]*(n_per_room[i][room][P] - conf_by_state_in_room[P]*alphaQ - vacc_by_state_in_room[P]*alphaT) +    // SP
                         Cpars[2]*(n_per_room[i][room][MI]-conf_by_state_in_room[MI]*alphaQ - vacc_by_state_in_room[MI]*alphaT) +   // S-MI
                         Cpars[3]*(n_per_room[i][room][SI]-conf_by_state_in_room[SI]*alphaQ - vacc_by_state_in_room[SI]*alphaT) +   // S-SI
                         Cpars[4]*n_per_room[i][room][H]);
            printf("Location %i Beta Multiplier %lf RateSE %lf 1/V %lf I/N %lf \n",i,a,b,c,d );
          }

          int from = pop[ind[j]][s];
          if(uniform()<exit_rate[from]*dt){            // If the person is selected to move
            double p = uniform();
            double temp=0;

            for(int to=0;to<n_states;to++){    // Loop over possible "to" states


              temp += (r[from][to]/exit_rate[from]);
              if(p<temp){                      // If such a transition must occur,
                // printf("From %s to %s\n",states[from], states[to] );
                pop[ind[j]][s] = to;           // Send this person to the "to" state.

                n_per_location[i][from]--; n_per_location[i][to]++;
                n[from]--; n[to]++;            // Change the values of n[from] and n[to]
                n_per_age[age][from]--;n_per_age[age][to]++;
                n_per_room[i][room][from]--; n_per_room[i][room][to]++;

                if(is_confined[ind[j]]==true){
                    conf_by_state_in_loc[from]--;      // change the number of confined in different states accordingly
                    conf_by_state_in_loc[to]++;
                    if(conf_by_state_in_loc[S]<0 || conf_by_state_in_loc[E]<0 || conf_by_state_in_loc[A]<0||conf_by_state_in_loc[P]<0||conf_by_state_in_loc[MI]<0||conf_by_state_in_loc[SI]<0||conf_by_state_in_loc[R]<0||conf_by_state_in_loc[H]<0){printf("Negative confs\n");}
                }
                vacc_by_state_in_room[from]--;vacc_by_state_in_room[to]++;

                if(to==H){
                  // Move them to the hospital

                  int hosp = randint(n_hospitals);                          // Choose a hospital at random
                  while(hosp==pop[ind[j]][w]){hosp = randint(n_hospitals);} // Make sure a HCW is not hospitalised in the same hospital
                  int hosp_room = randint(n_rooms[hosp]);                   // Random room in the hospital to send person to

                  people_linked_to[hosp].push_back(ind[j]); len_ppl_linked_to[hosp]++; // Link individual permanently to that hospital

                  if(is_confined[ind[j]]==true){
                    is_confined[ind[j]] = false;           // Remove confinement
                    person_isolated_time[ind[j]] = -1000;  // Added new: reset time of isolation.
                    conf_by_state_in_loc[H]--;             // Reduce number of confined hospitalised in this location if the person was confined (since they've moved to a hospital)
                    conf_by_state_in_room[H]--;
                  }


                  n_per_location[ pop[ind[j]][c] ][H]--; // Decrement number of "from" in current location
                  n_per_room[ pop[ind[j]][c] ][room][H]--;
                  pop[ind[j]][c] = hosp;                 // Send them to a random hospitals
                  current_room[ind[j]] = hosp_room;
                  n_per_location[ pop[ind[j]][c] ][H]++; // Increment number of "to" in current location
                  n_per_room[pop[ind[j]][c]][hosp_room][H]++; // Send them to a random room in the hospital.


                }
                else if(from==H && (to == R || to == D)){
                  // Remove confinement and move them home
                  int home = pop[ind[j]][h];
                  int home_room = randint(n_rooms[home]);

                   is_confined[ind[j]] = false; // Remove confinement

                   // Send recovered who were hospitalised home
                  n_per_location[ pop[ind[j]][c] ][R]--; // Decrement number of k (=R) in current locations
                  n_per_room[ pop[ind[j]][c] ][room][R]--; // Decrement number of people in this room
                  pop[ind[j]][c] = home;       // Send them home
                  current_room[ind[j]] = home_room;
                  n_per_location[ pop[ind[j]][c] ][R]++; // Increment number of k (=R) in current locations
                  n_per_room[ pop[ind[j]][c] ][ home_room ][R]++; // Send them to a random room in their house
                }

                break; // Exit the "to" loop, move to next person
              } // End if condition (if p<temp)
            } // End loop over "to" states
          } // End if condition (infection state == from)
        } // End loop over people in location

      /************ END CHANGE STATE OF POP ***********/

      }

      // if(N==0){continue;}
    } // End of location loop.

    /**************************** THE END OF DAY(S) ****************************/

    if(t>=day+1){

            /************* TARGETED TESTING ******************/

            if(n[R]>=begin_at/100 * n_pop){

              // Find all symptomatics

              int list_of_sym[n_pop]={};
              int n_sym = 0;

              int sym_rat_neg[n_pop]={};            // Symptomatics who tested negative on an RAT test that was declared greater than or equal to 7 days ago
              int n_srn = 0;

              int list_of_remaining[n_pop]={};       // Array to store remaining people for random testing
              int n_remaining  = 0;                  // Number of such people

              int not_eligible_for_testing = 0;      // Number of people not eligible for testing today

              for(int i=0; i<n_pop;i++){

                /******************** Find list of targets and remaining people **********************/

                if((pop[i][s]==MI||pop[i][s]==SI) && day>=next_test_date[i] && being_tested[i]==false && is_confined[i]==false){
                  list_of_sym[n_sym] = i; n_sym++;

                  // Retesting symptomatics who received a negative RAT test a week ago
                  if(day>=result_declared_date[i]+7 && last_test_type[i]==0 && last_test_result[i]==-1){sym_rat_neg[n_srn] = i; n_srn++;} // Make separate list of symptomatics who
                                                                                                                                          // last tested negative on a RAT, more than 7 days ago
                }
                else if(pop[i][s]!=H && pop[i][s]!=D && being_tested[i]==false && is_confined[i]==false){
                  list_of_remaining[n_remaining] = i; n_remaining++;
                }
                else{ not_eligible_for_testing++;}

                /****************************** Done finding lists **********************************/
              }

              list_of_sym[n_sym] = -1;                 // Set -1 to mark the end of this array
              list_of_remaining[n_remaining] = -1;     // Set -1 to mark end of this array

              sym_rat_neg[n_srn] = -1;                 // Set -1 to mark end of this array


              if(n_sym>0){
                int targeted_tests_done_today =   std::min(tests_remaining_today,n_sym);  // Tests done this dt is the minimum of the tests available and
                                                                                           // the people to be tested. As the day progresses, the tests
                                                                                           // available drops lower until it's 0, and no testing happens.

                shuffle(0, n_sym, list_of_sym); // Shuffle list of people to be tested.

                // for(int lo=targeted_tests_done_today;lo<n_sym;lo++){ list_of_remaining[n_remaining] = list_of_sym[lo];n_remaining++; } // Add those who couldn't be tested to the n_remaining
                                                                                                                                        // (this is a little pointess: if there are enough tests_conducted
                                                                                                                                        // targeted individuals will always be tested, until no tests
                                                                                                                                        // remain. Meaning if these people couldn't be tested in targeted testing
                                                                                                                                        // no random testing is going to happen! But anyway....)

                // list_of_remaining[n_remaining] = -1;  // Reset position of -1 to mark new end of this array

                  /*** REARRANGE list_of_sym SO THAT sym_rat_neg COME FIRST! *****/
                  int counter = 0;

                  for(int j=0;j<n_srn;j++){                 // For every sym_rat_neg
                    for(int k=0;k<n_sym;k++){               // Go over the remaining symptomatics
                      if(list_of_sym[k] == sym_rat_neg[j]){ // When you find this symptomatic in that list
                        int temp = list_of_sym[counter];       // Swap the lowest person on the list eligible
                        list_of_sym[counter] = list_of_sym[k];
                        list_of_sym[k] = temp;

                        counter++; break;                      // Increase the counter and break this loop
                      }
                    }
                  }


                for(int j=0; j<targeted_tests_done_today;j++){

                  int si = list_of_sym[j];  // Individual to test

                  if(day >= next_test_date[si] && being_tested[si]==false && tests_remaining_today>0){
                    // If so, perform a tests
                    being_tested[si] = true; tests_conducted++; tests_remaining_today--;

                    if(quarantine_when_sample_taken==true){   // Quarantine as soon as sample is taken
                      // Move them Home
                      n_per_location[pop[si][c]][pop[si][s]]--;
                      n_per_room[pop[si][c]][current_room[si]][s]--;
                      pop[si][c] = pop[si][h];
                      current_room[si] = randint(n_rooms[pop[si][c]]);
                      n_per_location[pop[si][c]][pop[si][s]]++;
                      n_per_room[pop[si][c]][current_room[si]][s]--;


                      is_confined[si] = quarantine_confined;
                      person_isolated_time[si] = day;             // Added new: time of isolation.

                      loc_confined[pop[si][h]] = lock_homes;      // Lock home depending on variable `lock_homes`.
                      loc_confined_time[pop[si][h]] = day;
                    }

                    // Targeted testing using HQ tests unless there are none
                    int test_type = 1;

                    if(hq_tests_today<=0){test_type = 0;} // If no HQ tests available, give them LQ (one of the two is guaranteed, since tests_remaining_today>0)

                    last_test_type[si] = test_type;       // NEW: Keeping track of last test type

                    if(test_type==0){
                      // Do a low quality test

                      lq_tests_conducted++; lq_tests_today--;

                      next_test_date[si]      = day + days_bw_lq_tests;   // Next day to be considered for a test
                      result_declared_date[si]= day + lq_delay;           // Days to wait before result result_declared_date

                      if(pop[si][s]>S && pop[si][s]<R){                   // If the person isn't susceptible or recovered or hospitalised
                        if(uniform()<lq_sens){test_result[si]=1;}         // and the test comes back positive, set their test_result
                        else{test_result[si]=-1;}                         // otherwise it's negative
                      }
                      else if(pop[si][s]!=H){                             // On the other hand, if they are S or R (not testing Hospitalised)
                        if(uniform()>lq_spec){test_result[si]=1;}         // and the test comes back false positive, set their test_result
                        else{test_result[si]=-1;}                         // otherwise it's negative
                      }
                    }

                    else if(test_type==1){
                      // Do a high quality test

                      hq_tests_conducted++; hq_tests_today--;

                      next_test_date[si]      = day + days_bw_hq_tests;   // Next day to be considered for a test
                      result_declared_date[si]= day + hq_delay;           // Days to wait before result result_declared_date

                      if(pop[si][s]>S && pop[si][s]<R){                  // If the person isn't susceptible or recovered or hospitalised
                        if(uniform()<hq_sens){test_result[si]=1;}         // and the test comes back positive, set their test_result
                        else{test_result[si]=-1;}                         // otherwise it's negative
                      }
                      else if(pop[si][s]!=H){                             // On the other hand, if they are S or R (not testing Hospitalised)
                        if(uniform()>hq_spec){test_result[si]=1;}         // and the test comes back false positive, set their test_result
                        else{test_result[si]=-1;}                         // otherwise it's negative
                      }
                    }
                  }
                }
              }

              /*************RANDOMLY TESTING POPULATION*************/

              // Test remaining people in population (in the array list_of_remaining) randomly

              int remaining_tests_done_today =   std::min(tests_remaining_today,n_remaining);  // Remaining tests done is the minimum of the tests available and
                                                                                               // available drops lower until it's 0, and no testing happens.
                                                                                               // the people to be tested. As the day progresses, the tests

              if(remaining_tests_done_today>0){               // If there are any tests remaining

                  shuffle(0,n_remaining,list_of_remaining);  // Shuffle list of remaining people

                  // Test the first "remaining_tests_done_today" people

                  for(int j=0; j<remaining_tests_done_today;j++){

                    int ri = list_of_remaining[j];           // Individual to test

                    // ** THE PART BELOW EXCEPT FOR SELECTING TESTS IS ESSENTIALLY THE SAME AS FOR TARGETED TESTING (with si -> ri)

                    being_tested[ri] = true; tests_conducted++; tests_remaining_today--;

                    if(quarantine_when_sample_taken==true){     // Quarantine as soon as sample is taken
                      // Move them Home
                      n_per_location[pop[ri][c]][pop[ri][s]]--;
                      n_per_room[pop[ri][c]][current_room[ri]][s]--;
                      pop[ri][c] = pop[ri][h];
                      current_room[ri] = randint(n_rooms[pop[ri][c]]);
                      n_per_location[pop[ri][c]][pop[ri][s]]++;
                      n_per_room[pop[ri][c]][current_room[ri]][s]--;

                      is_confined[ri] = quarantine_confined;
                      person_isolated_time[ri] = day;             // Added new: time of isolation.

                      loc_confined[pop[ri][h]] = lock_homes;      // Lock home depending on variable `lock_homes`.
                      loc_confined_time[pop[ri][h]] = day;

                    }

                    int test_type = randint(2);           // Returns either 0 and 1 with equal probability.

                    last_test_type[ri] = test_type;       // Keeping track of last test type

                    if(test_type==0 && lq_tests_today<=0){test_type=1;}      // If no LQ tests, give them HQ
                    else if(test_type==1 && hq_tests_today<=0){test_type=0;} // and vice versa

                    if(test_type==0){
                      // Do a low quality test

                      lq_tests_conducted++; lq_tests_today--;

                      next_test_date[ri]      = day + days_bw_lq_tests;   // Next day to be considered for a test
                      result_declared_date[ri]= day + lq_delay;           // Days to wait before result result_declared_date

                      if(pop[ri][s]>S && pop[ri][s]<R){                  // If the person isn't susceptible or recovered or hospitalised
                        if(uniform()<lq_sens){test_result[ri]=1;}         // and the test comes back positive, set their test_result
                        else{test_result[ri]=-1;}                         // otherwise it's negative
                      }
                      else if(pop[ri][s]!=H){                             // On the other hand, if they are S or R (not testing Hospitalised)
                        if(uniform()>lq_spec){test_result[ri]=1;}         // and the test comes back false positive, set their test_result
                        else{test_result[ri]=-1;}                         // otherwise it's negative
                      }
                    }

                    else if(test_type==1){
                      // Do a high quality test

                      hq_tests_conducted++; hq_tests_today--;

                      next_test_date[ri]      = day + days_bw_hq_tests;   // Next day to be considered for a test
                      result_declared_date[ri]= day + hq_delay;           // Days to wait before result result_declared_date

                      if(pop[ri][s]>S && pop[ri][s]<R){                  // If the person isn't susceptible or recovered or hospitalised
                        if(uniform()<hq_sens){test_result[ri]=1;}         // and the test comes back positive, set their test_result
                        else{test_result[ri]=-1;}                         // otherwise it's negative
                      }
                      else if(pop[ri][s]!=H){                             // On the other hand, if they are S or R (not testing Hospitalised)
                        if(uniform()>hq_spec){test_result[ri]=1;}         // and the test comes back false positive, set their test_result
                        else{test_result[ri]=-1;}                         // otherwise it's negative
                      }
                    }
                    // END OF RANDOM TESTING
                  }
              }

              /***********************DONE TESTING!***********************/

              for(int i=0;i<n_pop;i++){

                  /*********************DECLARING RESULTS********************/

                  if(day == result_declared_date[i] && being_tested[i]==true){

                    // First set them to no longer being tested
                    being_tested[i] = false;
                    results_declared++;

                    // Declare results_declared

                    if(test_result[i]==1 && pop[i][s] != H){  // If the result is positive, and the person hasn't already moved to a hospital
                      positives++;
                      is_confined[i] = quarantine_confined; // Confine them if quarantine_confined==true
                      person_isolated_time[i] = day;        // Time of isolation.

                      last_test_result[i] = 1;              // Keep track of retesting RAT negative symptomatics

                      // Move them home
                      n_per_location[pop[i][c]][pop[i][s]]--;              // Decremement number in current location
                      n_per_room[pop[i][c]][current_room[i]][pop[i][s]]--; // Decremement number in current room
                      if(pop[i][c] == pop[i][w]){ pop[i][c] = pop[i][h]; current_room[i] = randint(n_rooms[pop[i][h]]);}  // If they're at work, send them home to a random room.
                      n_per_location[pop[i][c]][pop[i][s]]++;              // Incremement number in current location.
                      n_per_room[pop[i][c]][current_room[i]][pop[i][s]]++; // Increment number in current room

                      loc_confined[pop[i][h]] = lock_homes;                // Reconfine their homes
                      loc_confined_time[pop[i][h]] = day;

                      next_test_date[i] = day + if_positive_test_after; // TO make compatible with applet. (CHECK!) NEW: changed from += 14 to today + 14 days.
                    }
                    else if(test_result[i] == -1 && pop[i][s] != H){       // If they're negative and have not been hospitalised (CHANGED!)
                      is_confined[i] = false;  // Remove confinement
                      person_isolated_time[i] = -1000;        // Reset time of isolation.

                      last_test_result[i] = -1;              // NEW: Added to keep track of retesting RAT negative symptomatics

                    }

                    if(pop[i][s]==H || pop[i][s]==D){is_confined[i]=false; person_isolated_time[i]=-1000;last_test_result[i]=0;} // If they've been hospitalised, do the same as for negative results,
                                                                                                                 // only set their last test result to 0.(CHECK!!!!)

                    test_result[i] = 0;        // Reset test result to 0.
                  }

                  /***********************DONE DECLARING!**********************/

                  /******** Remove confinement if 14 days have passed ********/

                  if(is_confined[i]==true && day >= person_isolated_time[i]+total_isolation_time){ // total_isolation_time = 14 days
                    is_confined[i]=false;
                    person_isolated_time[i] = -1000;
                  }
              }
            }

          //************ VACCINATION ******************//

          if(n[R]>=start_vacc/100 * n_pop){ // Start vaccination when number of recovered > start_vacc%.
            // printf("t = %.2f, vaccines administered = %d\n",t,vaccines_administered);
            int vaccs_available = dvr/100 * n_pop;
            // printf("Vaccines available %i\n", vaccs_available);

            vector<int> vacc_targets;
            int len_vacc_targets=0;

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
              printf("Error! Vaccination strategy wrong or not set! \n");
              exit(2);
            }

            vector<int> list_of_vacc;
            vector<int> list_of_others;

            int n_to_vacc = 0;
            int n_others = 0;

            int backup[n_states]={};

            for(int i=0;i<n_pop;i++){
              bool added = false;
              for(int j=0;j<len_vacc_targets;j++){
                int v_target = vacc_targets[j];
                if(pop[i][a]==v_target && pop[i][s]!=MI && pop[i][s]!=SI && pop[i][s]!=H && pop[i][s]!=D && !is_vaccinated[i] && !being_tested[i]){
                  backup[pop[i][s]]++;
                  list_of_vacc.push_back(i); n_to_vacc++; added = true; break;
                }
              }
              if(added == false && pop[i][s]!=MI && pop[i][s]!=SI && pop[i][s]!=H && pop[i][s]!=D && !is_vaccinated[i] && !being_tested[i]){
                list_of_others.push_back(i); n_others++;
              }
            }


            int vaccs_done_today = std::min(n_to_vacc,vaccs_available);
            shuffle(0,n_to_vacc,list_of_vacc);

            for(int i=0;i<vaccs_done_today;i++){
              int vi = list_of_vacc[i];
              is_vaccinated[vi]=true;
              vaccinated_on[vi] = day;
              vaccs_available--;
              vaccs_per_state[pop[vi][s]]++;

              vaccines_administered++;
              agewise_vaccines_administered[pop[vi][a]]++;
            }
            int random_vaccs_done_today=0;
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

            // printf("Day %3i Targeted vaccines %3i Random vaccines %3i Total remaining %i\n",day,vaccs_done_today, random_vaccs_done_today,vaccs_available );
            // printf("%i %i %i %i ",day,vaccs_done_today, random_vaccs_done_today,vaccs_available );
            // for(int i=0;i<n_states;i++){printf("%i ",vaccs_per_state[i] );}printf("\n");
            // printf("States of targets "); for(int i=0;i<n_states;i++){printf("%i ", backup[i]);}printf("\n");

          }


            // Lock or unlock homes
            int unlockedtoday = 0;
            for(int i=0; i<n_loc;i++){    // This can probably be added to the first location loop, if needed.

              if(loc_confined[i]==true && day-loc_confined_time[i] >= total_loc_confined_time){
                // for(int j = 0; j<n_pop; j++){if(pop[j][c]==i && is_confined[j]==true){printf("Problem in location %i\n", i);}}
                loc_confined[i]      = false; // Remove confinement
                loc_confined_time[i] = - 1000;
                unlockedtoday++;
              }

            }


      // Moving people around deterministically (WORK TO HOME)

      for(int i=0;i<n_pop;i++){
        if(pop[i][s]!=H && pop[i][s]!=D && is_confined[i] == false && loc_confined[pop[i][c]]==false){

          locations_moved++;

          int home_loc = pop[i][h];
          int work_loc = pop[i][w];
          if(pop[i][c]==work_loc){ pop[i][c] = home_loc; n_per_location[work_loc][pop[i][s]]--; n_per_location[home_loc][pop[i][s]]++; n_per_room[work_loc][current_room[i]][pop[i][s]]--; current_room[i] = randint(n_rooms[home_loc]); n_per_room[home_loc][current_room[i]][pop[i][s]]++; }
        }
      }

      midday_move_completed = false;
      first_move_done       = false;
      second_move_done      = false;
      third_move_done       = false;

      // Count locations currently quarantined
      int locked = 0;
      for(int hq=0;hq<n_loc;hq++){if(loc_confined_time[hq]>=0 && lock_homes==true){locked++;}}

      // Count people currently confined
      int currently_confined = 0;
      for(int hq=0;hq<n_pop;hq++){if(is_confined[hq]==true){currently_confined++;}}

      // Increment the day, write the output to an array, and reset the number of tests  //
      day++;

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

      output[day][0] = day; for(int q=0;q<n_states;q++){output[day][q+1]=n[q];}
      output[day][n_states+1] = vaccines_administered;
      output[day][n_states+2] = bs;
      output[day][n_states+3]=hq_tests_conducted;output[day][n_states+4]=lq_tests_conducted;output[day][n_states+5]=tests_remaining_today;output[day][n_states+6]=currently_confined;output[day][n_states+7]=unlockedtoday;output[day][n_states+8]=locked;
      tests_remaining_today = tests_available_daily;
      lq_tests_today = lq_tests_daily;
      hq_tests_today = hq_tests_daily;

      age_output[day][0] = day; int counter = 1; for(int r=0;r<max_age;r++){for(int q=0;q<n_states;q++){age_output[day][counter]=n_per_age[r][q];counter++;} age_output[day][counter] = agewise_vaccines_administered[r]; counter++;
      age_output[day][counter] = agewise_bs[r]; counter++;}

    }

    /********************* END OF END OF DAY(S) ********************************/

    t+=dt;
  }  // End of While Loop above.

  // Checking fraction of HCW who contracted the disease during the pandemic //
  int hcw_recovered = 0;
  int hcw = 0;
  for(int i=0;i<n_pop;i++){               // Find all HCW, and mark those that have recovered.
    if(pop[i][w]<n_hospitals){
      hcw++;
      if(pop[i][s]==R || pop[i][s]==D){hcw_recovered++;}
    }
  }

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

  int details[] = {quarantine_confined, lock_homes, quarantine_when_sample_taken, lq_tests_conducted, hq_tests_conducted, tests_conducted, results_declared, locations_moved, hcw_recovered,hcw};
  writetofile(output, age_output, tf, Tpars, begin_at, test_frac,cpu_time_used,details,iter);


}// End of TargetedTesting function.



int main() {

  char filename[1000]="synthetic_population.csv";
  // for(int tmp=0;tmp<10;tmp++){printf("%s\n","Blank");}
  // createPopulation(filename);

  /************ PARAMETERS ************/

  double lambda_S = 0.7;  // Changed from 0.35
  // double gamma    = 0.5;  // Fraction going from S->A

  double lambda_E = 1/4.5;

  double lambda_A = 1.0/8;//0.143;
  double lambda_P = 1/1.1;//0.5;
  // double delta    = 0.85;  // Fraction going from P->MI

  double lambda_MI = 1.0/8;//0.1;
  double lambda_SI = 1/1.5;//0.5;
  // double sigma     = 0.8;

  double lambda_H  = 1/18.1;//0.1;

  /**************************************/


  /************* AGE FACTOR *************/
  // double lambda_S_array[] = {0.34*lambda_S, 0.67*lambda_S, 1.0*lambda_S, 1.0*lambda_S, 1.0*lambda_S, 1.0*lambda_S, 1.0*lambda_S, 1.24*lambda_S, 1.47*lambda_S, 1.47*lambda_S};
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

  /**************************************/



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

  /**************************************/


  /********* TESTING PARAMETERS *********/

  bool lock_homes=false;                   // Change to true to quarantine homes
  bool quarantine_when_sample_taken=false; // Change to true to quarantine homes when sample is taken
  double begin_at  = 100;                  // Start testing when this percentage of the population has recovered
  double test_frac = 0.0;                  // Daily testing rate (in percentage)

  double rat_delay = 0;                    // RAT test delay
  double pcr_delay = 0;                    // PCR test delay

  double rat_fraction = 0;                 // RAT fraction in mixture
  double pcr_fraction = 1-rat_fraction;    // PCR fraction in mixture

  double Tpars[2][4] = {{0.5, 0.98, rat_delay, rat_fraction},   // {{RAT Sensitivity, RAT Specificity, RAT Test Delay, RAT Fraction in Mixture},
                        {1.0,  1.0, pcr_delay, pcr_fraction}};  //  {PCR Sensitivity, PCR Specificity, PCR Test Delay, PCR Fraction in Mixture}}

  /**********************************/

  int tf = 200;                            // Total simulation run time in days



  //** Monte Carlo run with the above specifications **//

  int mc_runs = 1;

  double unlockSchoolsAtArray[] = {0,10,20,30,40,50,60,70,80,90,100};
  int n_unlocks = sizeof(unlockSchoolsAtArray)/sizeof(unlockSchoolsAtArray[0]);

  double irf[] = {30.0/100.0, 50.0/100.0};
  int n_irf    = sizeof(irf)/sizeof(irf[0]);

  for(int i=0;i<mc_runs;i++){
    dvr = 0.2;
    vaccination_strategy=descending;

    for(int r=0;r<n_irf;r++){
      initial_recovered_fraction = irf[r];

      for(int j=0;j<n_unlocks;j++){
        unlockSchoolsAt = unlockSchoolsAtArray[j];

        createPopulation(filename);              // Create a random population with a fixed infection seed (default: 1%) and recovered fraction (default: 30%), set by initial_asymptomatic_fraction
                                                 // Resets the pop array, n_per_location array, and people_linked_to array.

        Targeted_Run(Tpars,
                     tf,
                     lock_homes,
                     quarantine_when_sample_taken,
                     begin_at,
                     test_frac,
                     0.0, //<-- Start vacc at
                     dvr,
                     unlockSchoolsAt,//<-- unlockSchoolsAt
                     0);
       }
    }


  }



 }

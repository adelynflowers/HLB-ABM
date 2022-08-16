#include "../headers/bioABM.h"
using namespace std;
namespace bioABM {

/*******************************************************************************************
 * MODEL PARAMETERS
 * Parameters and globals used in the model. These are
 * largely taken from Lee et. al (2015) and are available 
 * to change if future work demands it.
 * *****************************************************************************************/
#pragma region
int maxFlushAge = 30; // Maximum Flush Age
int flushEmerging = 20; // Flush Shoots Emerging
int eggAdultTransition_25C = 17;
int eggAdultTransition_28C = 14;
int eggAdultTransition = 17; //Egg to Adult Transition
int durationYoungFlush_25C = 13;
int durationYoungFlush_28C = 16;
int durationYoungFlush = 13; //Duration of young flush
float proportionMigrating = 0.4; //Proportion of Migrating Adults
float withinRowP = 0.95; //Within-row probability
float betweenRowP = 0.5; //Between-row probability
int eggDuration_25C = 4;
int eggDuration_28C = 3;
int eggDuration = 4; //Egg Duration
int nymphDuration_25C = 13;
int nymphDuration_28C = 11;
int nymphDuration = 13; //Nymph Duration
int shootCapacity = 40; //Maximum shoots a flush can support
int shootEggCapacity = 40; // Daily flush shoot egg capacity
int eggsPerFemaleAdult = 10; //Eggs per Female Adult
float transmissionFlushNymph = 0.083; //Transmission from flushes to nymphs
float transmissionAdultFlush = 0.3; //Transmission from adults to flushes
int latentPeriod = 15; //Latent Period, T=15 (see Supplement 1)
int asymptomaticLength = 548; //Days infected before becoming symptomatic, ~1.5 years
float eggSurvivalP_25C = 0.8614;
float eggSurvivalP_28C = 0.8343;
float eggSurvivalP = 0.8614; // Egg/nymph survival probability
float adultSurvivalP_25C = 0.9847;
float adultSurvivalP_28C = 0.9659;
float adultSurvivalP = 0.9847; //Adult survival probability OG: 0.9847
int nymphMinAgeToInfect = eggAdultTransition; // Minimum age to transmit HLB to flush, paper says adult??
int nymphMinAgeToBeInfected_25C = 9;
int nymphMinAgeToBeInfected_28C = 8;
int nymphMinAgeToBeInfected = 9; //Minimum age to receive HLB from flush
int modelDuration = 366 * 2; //Duration of model in days
int springFlushStart = 80;
int springFlushEnd = 140;
int summerFlushStart = 180;
int summerFlushEnd = 210;
int fallFlushStart = 250;
int fallFlushEnd = 265;
string invasionDays;
queue<int> invasionDays_q;
string invasionModalities;
queue<int> invasionModalities_q;
int carryingCapacity = 40000;
int invasionModality = 1;
int invasionGrove = 0;
double borderCrossingP = 0.01;
bool outputFlag = false;
float initialInfectedPortion = 0.18;
int initialNumPsyllids = 300;

// Lattice dimensions
// paper size: 69(nR) x 157 (rL)
const int rowLength = 69;
const int numRows = 157;
const int hBorders = 0;
const int vBorders = 0;


#pragma endregion

/**********************************************************************************************
 * GLOABLS
 * These are implementation related globals that support
 * the model as it runs
 * ********************************************************************************************/
#pragma region
// SQL and simple helper func
int experimentID;
void setExperimentID(int id) {
    experimentID = id;
}

bool isFlushingPeriod = false;
bool modelStarted = false;
int modelDay = -1;
ofstream csvFile;
string csvName;

boost::random::lagged_fibonacci607 rng(std::random_device{}());
boost::random::uniform_01<> gen;

//Migration probabilities moving clockwise, {North, East, South, West}
vector<double> migrationProbabilities = {betweenRowP/2, withinRowP/2, betweenRowP/2, withinRowP};
//Used to add to the current coordinate of migrating psyllid
vector<coord> migrationDiffs  = {coord(1,0), coord(0,1), coord(-1,0), coord(0,-1)};

#pragma endregion

/***********************************************************************************************
 * CLASSES AND STRUCTS
 * These are the core classes and structs that make up the model agents
 * *********************************************************************************************/
#pragma region

/********************************
 * FLUSH SHOOT
 * A single flush shoot, is part
 * of a greater flush patch
 * ******************************/
struct FlushShoot {
    int age = 0;
    bool infected = false;
    bool symptomatic = false;
    int numEggs = 0;
    int daysAsymptomatic = 0;
    bool alive = true;
    bool bark = false;
    int numPsyllids = 0;
    int numInfectedPsyllids = 0;
};

/*******************************
 * FLUSH PATCH
 * The core grid unit, holds flush
 * shoots. The term flush patch
 * is used interchangebly with
 * tree in our discussions/paper
 * *****************************/
struct FlushPatch {
    int age = 0;
    bool alive = true;
    int oldInfectedShoots = 0;
    int oldUninfectedShoots = 0;
    int numPsyllids_male = 0;
    int numPsyllids_female = 0;
    int numInfectedPsyllids_male = 0;
    int numInfectedPsyllids_female = 0;
    bool symptomatic = false;
    int daysInfected = 0;
    bool infectious = false;
    array<int,17> numNymphs;
    array<int,17> numInfectedNymphs;
    array<int, 30> numShoots;
    array<int, 30> numInfectedShoots;

    //Constructor
    FlushPatch() {      
        for (int i = 0; i < 17; i++) {
            numNymphs[i] = 0;
            numInfectedNymphs[i] = 0;
        }
        for (int i = 0; i < 30; i++) {
            numShoots[i] = 0;
            numInfectedShoots[i] = 0;
        }
    }

    //Used on tree death
    void kill() {
        oldInfectedShoots = 0;
        oldUninfectedShoots = 0;
        numShoots.fill(0);
        numInfectedShoots.fill(0);
        clearPsyllids();
        alive = false;
        symptomatic = false;
    }

    //Debugging helper
    bool validate() {
        for (int i = 0; i < 17; i++) {
            if (numNymphs[i] < 0 || numInfectedNymphs[i] < 0) {
                return false;
            }
        }
        for (int i = 0; i < 30; i++) {
            if (numShoots[i] < 0 || numInfectedShoots[i] < 0) {
                return false;
            }
        }
        return true;
    }

    //Age getter
    int getAge() {
        return age;
    }

    //Psyllid getter
    int getNumPsyllids() {
        int numPsyllids = 0;
        numPsyllids += numPsyllids_male;
        numPsyllids += numPsyllids_female;
        numPsyllids += numInfectedPsyllids_male;
        numPsyllids += numInfectedPsyllids_female;
        //numPsyllids += accumulate(numNymphs.begin(), numNymphs.end(), 0);
        //numPsyllids += accumulate(numInfectedNymphs.begin(), numInfectedNymphs.end(), 0);
        return numPsyllids;
    }

    //Calculate HLB severity, which is measured as the proportion of infected flushes
    double getHLBSeverity() {
        int uninfected = accumulate(numShoots.begin(), numShoots.end(), 0);
        int infected = accumulate(numInfectedShoots.begin(), numInfectedShoots.end(), 0);
        double hlbNum = (double)infected + (double)oldInfectedShoots;
        double hlbDenom = (double)uninfected + (double)oldUninfectedShoots + (double)hlbNum;
       
        if (hlbDenom == 0) {
            return 0;
        }
        else {
            assert((hlbNum / hlbDenom) >= 0 && (hlbNum / hlbDenom) <= 1);
            double severity = hlbNum / hlbDenom;
            return severity;
        }

    }

    //Remove all psyllids
    void clearPsyllids() {
        numPsyllids_female = 0;
        numPsyllids_male = 0;
        numInfectedPsyllids_female = 0;
        numInfectedPsyllids_male = 0;
        numNymphs.fill(0);
        numInfectedNymphs.fill(0);
        return;
    }

    //Wrapper used somewhere, unsure why this is here
    int getTotalPsyllids() {
        return getNumPsyllids();
    }

    //Helper to add a psyllid based on its characteristics
    void placePsyllid(bool female, bool adult, bool infected, int age = -1) {
        if (adult && female && infected) {
            numInfectedPsyllids_female++;
        }
        else if (adult && female && !infected) {
            numPsyllids_female++;
        }
        else if (adult && !female && infected) {
            numInfectedPsyllids_male++;
        }
        else if (adult && !female && !infected) {
            numPsyllids_male++;
        }
        else if (!adult && !infected) {
            numNymphs[age]++;
        }
        else if (!adult && infected) {
            numInfectedNymphs[age]++;
        }
    }
};

//Used in determining the number of neighbors a cell has
enum PositionType {MIDDLE, EDGE, CORNER};

typedef boost::tuple<int, int> coord;
typedef array<FlushPatch, rowLength> FlushRow;
//11x25 lattice that represents flush patches
array<FlushRow, numRows> lattice; 

#pragma endregion

/***********************************************************************************************
 * HELPERS
 * These are helper functions that help keep the code
 * somewhat readable. Usually in place so that the behavioral model can call these to get
 * values
 * *********************************************************************************************/
#pragma region

/****************************************************
* MULTINOM
* Generates a pull from a multinomial distribution with
* k components, N size, and vector of probabilities p.
* Results are stored in res parameter, which must be of
* length k. Based on description of conditional
* distribution method in Davis (1993)
******************************************************/
void multinom(int k, int N, double* p, int* res) {
    for (int i = 0; i < k; i++) {
        int xSum = 0;
        double prob_denom_alt = 1;
        for (int j = 0; j < i; j++) {
            prob_denom_alt -= p[j];
            xSum += res[j];
        }
        int trials = N - xSum;
        double binom_prob = p[i] / prob_denom_alt;
        boost::binomial_distribution<int> puller(trials, binom_prob);
        res[i] = puller(rng);
    }
}


/*****************************************
 * INT RAND
 * Used for generating random integers in 
 * [min,max]
 * ***************************************/
int intRand(const int& min, const int& max) {
    static thread_local std::mt19937 generator(std::random_device{}());
    std::uniform_int_distribution<int> distribution(min, max);
    return distribution(generator);
}

/****************************************
 * DOUBLE RAND
 * Used for generating random doubles in
 * [min, max)
 * **************************************/
double doubleRand(double min, double max) {
    static thread_local std::mt19937 generator(std::random_device{}());
    std::uniform_real_distribution<double> distribution(min, max);
    //return gen(rng);
    return distribution(generator);
}

/**************************************************************
 * Discrete Probability Match
 * Takes in a discrete vector of probabilities, and returns the
 * corresponding index based on a uniform pull 
 * and those probabilities
 * ***********************************************************/
int discreteProbabilityMatch(vector<double> probabilities) {
    double pull;
    pull = doubleRand(0, 1);
    double cumSum = 0;
    int resultIdx = -1;
    for (int i = 0; i < probabilities.size(); i++) {
        if (pull >= cumSum && pull <= (cumSum + probabilities[i])) {
            resultIdx = i;
            break;
        }
        else {
            cumSum += probabilities[i];
        }
    }
    if (resultIdx < 0) {
        if (pull == 1) {
            resultIdx = probabilities.size() - 1;
        }
        else {
            cout << "Well there's ya problem\n";
            assert(resultIdx >= 0);
        }

    }
    return resultIdx;
}
/*****************************************************
 * STRING SPLIT
 * Homemade string split function ripped from stack 
 * overflow
 * ***************************************************/
vector<string> split(string s, string delimiter) {
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    string token;
    vector<string> res;

    while ((pos_end = s.find (delimiter, pos_start)) != string::npos) {
        token = s.substr (pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back (token);
    }

    res.push_back (s.substr (pos_start));
    return res;
}

/*********************
 * Grid Size Helpers
 *********************/
int getNumRows() {
    return numRows;
}

int getRowLength() {
    return rowLength;
}
/***********************
* Get Flushing Bounds
************************/
int getFallStart() {
    return fallFlushStart;
}
int getSpringStart() {
    return springFlushStart;
}
int getSummerStart() {
    return summerFlushStart;
}
/*************************************************************
* Get Age at
* Returns flush patch age at coordinates
*************************************************************/
int getAgeAt(int i, int j, int differential) {
    return (lattice[i][j].getAge() + differential) / 365;
}

bool isTreeAlive(int i, int j) {
    if (isValidCoordinate(coord(i,j))) {
        return lattice[i][j].alive;
    }
    else {
        return false;
    }
}

/***************************************
 * Is Symptomatic
 * Returns true if flush patch at 
 * coordinates is symptomatic
 * *************************************/
bool isSymptomatic(int i,int j) {
    return lattice[i][j].symptomatic;
}
/************************************************************
* Get Severity At
* Returns flush patch severity at coordinates
*************************************************************/
double getSeverityAt(int i, int j) {
    return lattice[i][j].getHLBSeverity();
}

/**************************************************************
* Get Psyllids At
* Returns number of psyllids at coordinates
**************************************************************/
double getPsyllidsAt(int i, int j) {
    if (i > numRows || j > rowLength || i < 0 || j < 0) {
        return 0;
    }
    else {
        return lattice[i][j].getTotalPsyllids();
    }
}


/**************************************************************
* Get Model Day
* For econ wrapper
***************************************************************/
int getModelDay() {
    return modelDay;
}

/***************************************************************
* Get Model Duration
* For econ wrapper
***************************************************************/
int getModelDuration() {
    return modelDuration;
}

/********************************
 * Write CSV batch
 * Writes the current model state
 * to the opened csv file
 * *******************************/
void write_csv_batch() {
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < rowLength; j++) {
            int numPsyllids = lattice[i][j].numPsyllids_male + 
                lattice[i][j].numPsyllids_female + 
                accumulate(lattice[i][j].numNymphs.begin(), lattice[i][j].numNymphs.end(), 0);
            int numInfectedPsyllids = lattice[i][j].numInfectedPsyllids_male 
                + lattice[i][j].numInfectedPsyllids_female +
                accumulate(lattice[i][j].numInfectedNymphs.begin(), lattice[i][j].numInfectedNymphs.end(), 0);
            double severity = lattice[i][j].getHLBSeverity();
            if (severity > 1) {
                cout << "WARNING: HLB ERROR AT (" << i << ", " << j << ")\n";
            }
            csvFile << modelDay << "," << i << "," << j << "," << numPsyllids << "," 
                << numInfectedPsyllids << "," << std::fixed 
                << setprecision(5) << severity << "," << experimentID <<  "," << lattice[i][j].alive << "," << lattice[i][j].symptomatic << "\n";
        }
    }
}

/*********************************
 * Count psyllids
 * Gets a total model psyllid count
 * *******************************/
int countPsyllids() {
    int num = 0;
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < rowLength; j++) {
            int cnt = lattice[i][j].getTotalPsyllids();
            if (cnt > 0) {
                cout << "Found " << cnt << " at " << i << "," << j <<  " on day " << modelDay << endl;
                num += cnt;
            }
        }
    }
    return num;
}

/********************************
 * Validate lattice
 * Validates each cell of the 
 * lattice
 * ******************************/
bool validateLattice() {
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < rowLength; j++) {
            if (!lattice[i][j].validate()) {
                cout << "Failed validation at (" << i << ", " << j << ")\n";
                return false;
            }
        }
    }
    return true;
}

/*****************************
 * Finish run
 * Closes the output file
 * ***************************/
void finishRun() {
    csvFile.close();
}

#pragma endregion

/*************************************************************************************************
 * MODEL FUNCTIONALITY
 * These are functions that the model uses to run, but are not activities
 * ***********************************************************************************************/
#pragma region
/*************************************************************
 * Update Seasonal Parameters
 * Updates parameters that vary based on operative temperature
 * ***********************************************************/
void updateSeasonalParameters(int relativeT) {
    if ( (relativeT >= springFlushStart && relativeT <= springFlushEnd) ||
         (relativeT >= fallFlushStart && relativeT <= fallFlushEnd) ) {
             eggAdultTransition = eggAdultTransition_25C;
             durationYoungFlush = durationYoungFlush_25C;
             eggDuration = eggDuration_25C;
             nymphDuration = nymphDuration_25C;
             eggSurvivalP = eggSurvivalP_25C;
             adultSurvivalP = adultSurvivalP_25C;
         }
    else if (relativeT >= summerFlushStart && relativeT <= summerFlushEnd) {
        eggAdultTransition = eggAdultTransition_28C;
        durationYoungFlush = durationYoungFlush_28C;
        eggDuration = eggDuration_28C;
        nymphDuration = nymphDuration_28C;
        eggSurvivalP = eggSurvivalP_28C;
        adultSurvivalP = adultSurvivalP_28C;
    }
}

/**************************************************************
 * Is Flushing Period
 * ************************************************************/
void setFlushingPeriod(int period) {
    int t = period % 365;
    if ( (t >= springFlushStart && t <= springFlushEnd) ||
         (t >= summerFlushStart && t <= summerFlushEnd) || 
         (t >= fallFlushStart   && t <= fallFlushEnd)) {
        isFlushingPeriod = true;
        updateSeasonalParameters(t);
    }
    else {
        isFlushingPeriod = false;
    }
}


/**************************************************************
 * crossesBorder
 * Determines if movement between the two cells is a border 
 * crossing
 * ************************************************************/
bool crossesBorder(coord a, coord b) {
    int aQuadrant_i = a.get<0>() / (numRows / (hBorders+1));
    int aQuadrant_j = a.get<1>() / (rowLength / (vBorders+1));
    int bQuadrant_i = b.get<0>() / (numRows / (hBorders+1));
    int bQuadrant_j = b.get<1>() / (rowLength / (vBorders+1));
    if (aQuadrant_i != bQuadrant_i || aQuadrant_j != bQuadrant_j) {
        return true;
    }
    else {
        return false;
    }
}

/***************************************************************
 * Is Valid Coordinate
 * ************************************************************/
bool isValidCoordinate(coord pos) {
    int row = pos.get<0>();
    int col = pos.get<1>();
    if (row >= 0 && row < numRows && col >= 0 && col < rowLength) {
        if (lattice[row][col].alive) {
            return true;
        }
        else {
            return false;
        }
    }
    else {
        return false;
    }
}

/**************************************************************
 * Determine position type
 * Determines if a coordinate is a middle, edge, or corner
 * ************************************************************/
PositionType determinePositionType(coord pos) {
    int row = pos.get<0>();
    int col = pos.get<1>();
    array<coord, 4> neighbors = { coord(row,col+1), coord(row, col-1), coord(row+1, col), coord(row-1, col)};
    int numNeighbors = 0;
    for (int i = 0; i < neighbors.size(); i ++) {
        if (isValidCoordinate(neighbors[i])) {
            numNeighbors += 1;
        }
    }
    if (numNeighbors < 2 || numNeighbors > 4) {
        assert(numNeighbors >= 2 && numNeighbors <= 4);
    }
    switch(numNeighbors) {
        case 2:
            return PositionType::CORNER;
        case 3:
            return PositionType::EDGE;
        case 4:
            return PositionType::MIDDLE;
        default:
            assert(false); //ERROR
            return PositionType::MIDDLE;
    }
}

/*************************************************
 * Initialize Lattice
 * ***********************************************/
void initializeLattice() {
    for (int i = 0; i < numRows; i++) {
        FlushRow row;
        for (int j = 0; j < rowLength; j++) {
            FlushPatch patch;
            row[j] = patch;
        }
        lattice[i] = row;
    }
}

/***************************************************
* getGroveBounds
* Returns the row and column upper/lower bounds
* based on a two digit identifier, indices are
* 0: row bounds
* 1: column bounds
* Bounds are lower inclusive, upper exclusive
****************************************************/
vector<coord> getGroveBounds(int identifier) {
    int rowID = identifier / 10;
    int colID = identifier % 10;
    
    //First row bounds
    int gamma_r = numRows / (hBorders + 1);
    int rowLB = gamma_r * rowID;
    int rowUB = gamma_r * (rowID + 1);
    coord rowBounds(rowLB, rowUB);

    //Same for column
    int gamma_c = rowLength / (vBorders + 1);
    int colLB = gamma_c * colID;
    int colUB = gamma_c * (colID + 1);
    coord colBounds(colLB, colUB);

    //Package and return
    vector<coord> bounds;
    bounds.push_back(rowBounds);
    bounds.push_back(colBounds);
    return bounds;
}

#pragma endregion

/**************************************************************************************************
 * MANAGEMENT FUNCTIONS
 * These are functions associated with ACP/HLB management strategies
 * ************************************************************************************************/
#pragma region
/*****************************************************************
* rogueTreeAt
******************************************************************/
void rogueTreeAt(int i, int j) {
    if (isValidCoordinate(coord(i,j))) {
        lattice[i][j].kill();
    }
    return;
}

/****************************************************************
* Spray trees
******************************************************************/
void sprayTrees(double efficacy, vector<coord> locations) {
    int psyllidsRemoved = 0;
    //cout << "Spraying on day: " << modelDay << " with efficacy " << efficacy << endl;
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < rowLength; j++) {
            //int r = locations[i].get<0>();
            //int c = locations[i].get<1>();
        
            //int beforePsyllids = lattice[r][c].getTotalPsyllids();
            lattice[i][j].numPsyllids_male = ceil((1.0 - efficacy) * (double)lattice[i][j].numPsyllids_male);
            lattice[i][j].numPsyllids_female = ceil((1.0 - efficacy) * (double)lattice[i][j].numPsyllids_female);
            lattice[i][j].numInfectedPsyllids_male = ceil((1.0 - efficacy) * (double)lattice[i][j].numInfectedPsyllids_male);
            lattice[i][j].numInfectedPsyllids_female = ceil((1.0 - efficacy) * (double)lattice[i][j].numInfectedPsyllids_female);
            for (int k = 0; k < 17; k++) {
                lattice[i][j].numNymphs[k] = ceil((1.0 - efficacy) * (double)lattice[i][j].numNymphs[k]);
                lattice[i][j].numInfectedNymphs[k] = ceil((1.0 - efficacy) * (double)lattice[i][j].numInfectedNymphs[k]);
            }
        }
    }
        /*int startingMales = lattice[r][c].numPsyllids_male;
        for (int j = 0; j < startingMales; j++) {
            if (doubleRand(0, 1) <= efficacy) {
                lattice[r][c].numPsyllids_male--;
            }
        }
        int startingFemales = lattice[r][c].numPsyllids_female;
        for (int j = 0; j < startingFemales; j++) {
            if (doubleRand(0, 1) <= efficacy) {
                lattice[r][c].numPsyllids_female--;
            }
        }
        int startingMales_i = lattice[r][c].numInfectedPsyllids_male;
        for (int j = 0; j < startingMales_i; j++) {
            if (doubleRand(0, 1) <= efficacy) {
                lattice[r][c].numInfectedPsyllids_male--;
            }
        }
        int startingFemales_i = lattice[r][c].numInfectedPsyllids_female;
        for (int j = 0; j < lattice[r][c].numInfectedPsyllids_female; j++) {
            if (doubleRand(0, 1) <= efficacy) {
                lattice[r][c].numInfectedPsyllids_female--;
            }
        }
        for (int j = 0; j < 17; j++) {
            int startingNymphs = lattice[r][c].numNymphs[j];
            for (int k = 0; k < startingNymphs; k++) {
                if (doubleRand(0, 1) <= efficacy) {
                    lattice[r][c].numNymphs[j]--;
                }
            }
            int startingNymphs_i = lattice[r][c].numInfectedNymphs[j];
            for (int k = 0; k < startingNymphs_i; k++) {
                if (doubleRand(0, 1) <= efficacy) {
                    lattice[r][c].numInfectedNymphs[j]--;
                }
            }
            //lattice[r][c].numNymphs[j] = ceil((1.0 - efficacy) * (double)lattice[r][c].numNymphs[j]);
            //lattice[r][c].numInfectedNymphs[j] = ceil((1.0 - efficacy) * (double)lattice[r][c].numInfectedNymphs[j]);
        }
        int afterPsyllids = lattice[r][c].getTotalPsyllids();
        psyllidsRemoved += beforePsyllids - afterPsyllids;*/
    //}
    //cout << "Spray removed " << psyllidsRemoved << endl;
}

#pragma endregion

/***************************************************************************************************
 * PSYLLID INVASION
 * These are functions that are used for psyllid invasion events
 * ************************************************************************************************/
#pragma region
/***************************************************
* uniformPsyllidDistribution
* Distributes a number of uninfected psyllids evenly 
* amongst cells not in the specified vector
* *************************************************/
void uniformPsyllidDistribution(double percent, int numPsyllids, vector<coord> occupied, int groveID) {
    //Determine number of trees to infect
    int rowsPerGrove = numRows / (hBorders + 1);
    int rowLPerGrove = rowLength / (vBorders + 1);
    int totalTrees = rowsPerGrove * rowLPerGrove;
    int availableTrees = totalTrees - occupied.size();
    int treesToInfect = floor((double)availableTrees * percent);

    //Set up coordinates available
    vector<coord> groveTrees = getGroveBounds(groveID);
    int psyllidsPerTree = floor((double)numPsyllids / (double)treesToInfect);

    //For using std::find
    
    vector<coord> cells;
    for (int i = 0; i < treesToInfect; i++) {
        //Choose a tree
        coord c;
        c = coord(intRand(groveTrees[0].get<0>(), groveTrees[0].get<1>() - 1), intRand(groveTrees[1].get<0>(), groveTrees[1].get<1>() - 1));
        //Reroll until it's not one of the off-limit trees or one we already picked
        while (find(occupied.begin(), occupied.end(), c) != occupied.end() ||
               find(cells.begin(), cells.end(), c) != cells.end()) {
            c = coord(intRand(groveTrees[0].get<0>(), groveTrees[0].get<1>() - 1), intRand(groveTrees[1].get<0>(), groveTrees[1].get<1>() - 1));
        }
        cells.push_back(c);
  
    }
    for (int i = 0; i < cells.size(); i++) {
        for (int j = 0; j < psyllidsPerTree; j++) {
            bool infected = doubleRand(0, 1) < initialInfectedPortion;
            bool female = doubleRand(0, 1) < 0.5;
            lattice[cells[i].get<0>()][cells[i].get<1>()].placePsyllid(female, true, infected);
        }
    }
}


/***************************************************
* invasionModality1
* 200 psyllids, evenly distributed 
* amongst 4 trees in the SW corner
****************************************************/
vector<coord> invasionModality1(int groveID) {
    int psyllidsPerTree = 50;
    vector<coord> bounds = getGroveBounds(groveID);
    vector<coord> cells;

    //Collect the 4 trees in the SW corner
    cells.push_back(coord(bounds[0].get<1>() - 1, bounds[1].get<0>()));
    cells.push_back(coord(bounds[0].get<1>() - 2, bounds[1].get<0>()));
    cells.push_back(coord(bounds[0].get<1>() - 1, bounds[1].get<0>() + 1));
    cells.push_back(coord(bounds[0].get<1>() - 2, bounds[1].get<0>() + 1));

    //Place psyllids
    for (int i = 0; i < cells.size(); i++) {
        for (int j = 0; j < psyllidsPerTree; j++) {
            bool infected = doubleRand(0,1) < initialInfectedPortion;
            bool female = doubleRand(0,1) < 0.5;
            lattice[cells[i].get<0>()][cells[i].get<1>()].placePsyllid(female, true, infected);
        }
    }
    //cout << "Psyllids placed: " << numPsyllids << endl;
    return cells;

}

/***************************************************************
* invasionModality2
*   Modality 1 AND 35% of randomly selected trees from
*   remaining patches are occupied by 200 uninfected 
*   ACP distributed evenly
****************************************************************/
void invasionModality2(int groveID) {
    vector<coord> cells = invasionModality1(groveID);
    uniformPsyllidDistribution(0.35, 200, cells, groveID);
}

/****************************************************************
* invasionModality3
* 200 psyllids, placed on 25% of trees on the southern edge
* and 100% trees on the eastern edge
****************************************************************/
vector<coord> invasionModality3(int groveID) {
    vector<coord> cells;
    //Get bounds
    vector<coord> bounds = getGroveBounds(groveID);

    //Determine number of psyllids per tree
    int southernTrees = bounds[1].get<1>() - bounds[1].get<0>();
    southernTrees = floor((double)southernTrees * 0.25);
    int easternTrees = bounds[0].get<1>() - bounds[0].get<0>();
    int totalTrees = southernTrees + easternTrees;
    int psyllidsPerTree = floor(200 / (double)totalTrees);

    //Add eastern trees
    for (int i = bounds[0].get<0>(); i < bounds[0].get<1>(); i++) {
        coord c(i, bounds[1].get<1>() - 1);
        cells.push_back(c);
    }

    //Add southern trees
    for (int i = 0; i < southernTrees; i++) {
        coord c(bounds[0].get<1>() - 1, intRand(bounds[1].get<0>(), bounds[1].get<1>() - 2));
        while (find(cells.begin(), cells.end(), c) != cells.end()) {
            c = coord(bounds[0].get<1>() - 1, intRand(bounds[1].get<0>(), bounds[1].get<1>() - 2));
        }
        cells.push_back(c);
    }

    //Distribute psyllids
    for (int i = 0; i < cells.size(); i++) {
        for (int j = 0; j < psyllidsPerTree; j++) {
            bool infected = doubleRand(0, 1) < initialInfectedPortion;
            bool female = doubleRand(0, 1) < 0.5;
            lattice[cells[i].get<0>()][cells[i].get<1>()].placePsyllid(female, true, infected);
        }
    }
    return cells;
}

/***************************************************************
* invasionModality4
* Modality 3 AND 35% of randomly selected trees from 
* remaining patches are occupied by 200 uninfected
* ACP distributed evenly
* **************************************************************/
void invasionModality4(int groveID) {
    vector<coord> cells = invasionModality3(groveID);
    uniformPsyllidDistribution(0.35, 200, cells, groveID);
}

/****************************************************************
* invasionModality5
* 9 trees distributed around the center of the grove
* are occupied by 198 psyllids
* ***************************************************************/
vector<coord> invasionModality5(int groveID) {
    //Collect the 9 middle cells
    vector<coord> cells;
    vector<coord> bounds = getGroveBounds(groveID);
    int center_i = ceil(  ( (double)bounds[0].get<1>() - 1 - (double)bounds[0].get<0>() ) / 2) + bounds[0].get<0>();
    int center_j = ceil(  ( (double)bounds[1].get<1>() - 1 - (double)bounds[1].get<0>() ) / 2) + bounds[1].get<0>();
    for (int i = center_i - 1; i <= center_i + 1; i++) {
        for (int j = center_j - 1; j <= center_j + 1; j++) {
            coord c(i, j);
            cells.push_back(c);
        }
    }
    //Distribute psyllids
    int psyllidsPerTree = 22;
    for (int i = 0; i < cells.size(); i++) {
        for (int j = 0; j < psyllidsPerTree; j++) {
            bool infected = doubleRand(0, 1) < initialInfectedPortion;
            bool female = doubleRand(0, 1) < 0.5;
            lattice[cells[i].get<0>()][cells[i].get<1>()].placePsyllid(female, true, infected);
        }
    }

    return cells;

}

/***************************************************************
* invasionModality6
* Modality 5 AND 35% of randomly selected trees from
* remaining patches are occupied by 200 uninfected 
* ACP distributed evenly
****************************************************************/
void invasionModality6(int groveID) {
    vector<coord> cells = invasionModality5(groveID);
    uniformPsyllidDistribution(0.35, 200, cells, groveID);
}


/***************************************************************
 * Place Initial Psyllids
 * INVASION SPATIAL MODALITIES:
 *      1: 200 psyllids, evenly distributed 
 *          amongst 4 trees in the SW corner
 *      2: Modality 1 AND 35% of randomly selected trees from
 *         remaining patches are occupied by 200 uninfected 
           ACP distributed evenly
 *      3: 200 psyllids, placed on 25% of trees on the southern edge
 *         and 100% of trees on the eastern edge
 *      4: Modality 3 AND 35% of randomly selected trees from 
 *         remaining patches are occupied by 200 uninfected
 *         ACP distributed evenly
 *      5: 10 trees distributed around the center of the grove
 *         are occupied by 200 psyllids
 *      6: Modality 5 AND 35% of randomly selected trees from
 *         remaining patches are occupied by 200 uninfected 
 *         ACP distributed evenly
 * ************************************************************/
void placeInitialPsyllids(int invasionModality, int groveID) {
    switch (invasionModality) {
        case 1:
            invasionModality1(groveID);
            break;
        case 2:
            invasionModality2(groveID);
            break;
        case 3:
            invasionModality3(groveID);
            break;
        case 4:
            invasionModality4(groveID);
            break;
        case 5:
            invasionModality5(groveID);
            break;
        case 6:
            invasionModality6(groveID);
            break;
        default:
            invasionModality1(groveID);
    }
    return;
}
#pragma endregion

/****************************************************************************************************
 * INITIALIZATION
 * These are functions related to initializing the model
 * **************************************************************************************************/
#pragma region
/************************************************
 * Initialize Model
 * Called from econ
 * **********************************************/
void initializeModel() {
    initializeLattice();
    return;
}

/************************************
 * Initialize CSV
 * Open the csv file for writing
 * **********************************/
void initializeCSV() {
    try {
        csvFile.open(csvName);
    }
    catch (exception e) {
        cout << e.what() << endl;
    }
    csvFile << "t,i,j,numPsyllids,numInfectedPsyllids,hlbSeverity, experiment_id,alive,symptomatic\n";
}

/*****************************************
 * Parse parameter file
 * Imports parameters from a json config
 * ***************************************/
void parseParameterFile(string fileName) {
    try {
        ifstream is(fileName);
        cereal::JSONInputArchive archive(is);
        archive(maxFlushAge, flushEmerging, eggAdultTransition,
            durationYoungFlush, proportionMigrating, withinRowP,
            betweenRowP, eggDuration, nymphDuration, shootCapacity,
            shootEggCapacity, eggsPerFemaleAdult, transmissionFlushNymph,
            transmissionAdultFlush, latentPeriod, eggSurvivalP, adultSurvivalP,
            nymphMinAgeToInfect, nymphMinAgeToBeInfected,
            modelDuration, csvName, initialInfectedPortion, initialNumPsyllids,
            invasionDays, carryingCapacity, borderCrossingP, springFlushStart,
            springFlushEnd, summerFlushStart, summerFlushEnd, fallFlushStart,
            fallFlushEnd, invasionModalities, invasionGrove, outputFlag);

        vector<string> invasionDays_str = split(invasionDays, ",");
        for (int i = 0; i < invasionDays_str.size(); i++) {
            invasionDays_q.push(stoi(invasionDays_str[i]));
        }
        vector<string> invasionModalities_str = split(invasionModalities, ",");
        for (int i = 0; i < invasionModalities_str.size(); i++) {
            invasionModalities_q.push(stoi(invasionModalities_str[i]));
        }
    }
    catch (exception e) {
        cout << "ERROR WITH BIO JSON" << e.what() << endl;
        exit(-1);
    }
}
#pragma endregion


/****************************************************************************************************
 * ACTIVITIES 
 * These are the core functions of the model, associated with activities
 * ***************************************************************************************************/
#pragma region

/************************************************
 * Birth New Flush
 * Emergence of new flush:
 * - For each flush patch, birth 20 shoots
 * **********************************************/
void birthNewFlush() {
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < rowLength; j++) {
            //Dont birth flush on dead trees
            if (!lattice[i][j].alive) {
                continue;
            }
            double birthInfectChance = lattice[i][j].getHLBSeverity();
            //assert(birthInfectChance < 1);
            if (birthInfectChance > 0) {
                for (int k = 0; k < flushEmerging; k++) {
                    if (doubleRand(0, 1) <= birthInfectChance) {
                        lattice[i][j].numInfectedShoots[0]++;
                    }
                    else {
                        lattice[i][j].numShoots[0]++;
                    }
                }
            }
            else {
                lattice[i][j].numShoots[0] += flushEmerging;
            }
        }
    }
}

/*************************************************
 * Migration
 * Psyllids migrate between flush patches:
 * - Calculate probabilities for rook contiguous
 *   locations based on current model state
 * - For each subpopulation of psyllid 
 *   (uninfected,infected,male,female),
 *   use these probabilities with a multinomial 
 *   distribution to determine the migration of
 *   psyllids
 * ***********************************************/
void migration() {
    //ROUGH FIX: copied lattice for storing differentials
    //INDEX AS: mtrx[y*numRows + x]
    int* diffMales = new int[numRows*rowLength];
    int* diffMales_i = new int[numRows*rowLength];
    int* diffFemales = new int[numRows*rowLength];
    int* diffFemales_i = new int[numRows*rowLength];
    std::fill(diffMales, diffMales + (numRows*rowLength), 0);
    std::fill(diffFemales, diffFemales + (numRows*rowLength), 0);
    std::fill(diffMales_i, diffMales_i + (numRows*rowLength), 0);
    std::fill(diffFemales_i, diffFemales_i + (numRows*rowLength), 0);


    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < rowLength; j++) {
            //1. Calculate probabilities 
            double leftProb, rightProb, upProb, downProb, noneProb;
            int total = lattice[i][j].getTotalPsyllids();

            //1a. Determine left right probabilities
            double inRowProb = proportionMigrating * withinRowP;
            if (isValidCoordinate(coord(i, j + 1)) && isValidCoordinate(coord(i, j - 1))) {
                leftProb = rightProb = inRowProb * 0.5;
            }
            else if (!isValidCoordinate(coord(i, j + 1)) && isValidCoordinate(coord(i, j - 1))) {
                leftProb = inRowProb;
                rightProb = 0;
            }
            else if (isValidCoordinate(coord(i, j + 1)) && !isValidCoordinate(coord(i, j - 1))) {
                rightProb = inRowProb;
                leftProb = 0;
            }
            else {
                leftProb = rightProb = 0;
            }
            
            //1b. Determine up-down probabilities
            double betweenRowProb = proportionMigrating * betweenRowP;
            if (isValidCoordinate(coord(i + 1, j)) && isValidCoordinate(coord(i - 1, j))) {
                upProb = downProb = betweenRowProb * 0.5;
            }
            else if (!isValidCoordinate(coord(i + 1, j)) && isValidCoordinate(coord(i - 1, j))) {
                downProb = betweenRowProb;
                upProb = 0;
            }
            else if (isValidCoordinate(coord(i + 1, j)) && !isValidCoordinate(coord(i - 1, j))) {
                upProb = betweenRowProb;
                downProb = 0;
            }
            else {
                upProb = downProb = 0;
            }
            
            //1c. Apply border crossing modifiers
            if (upProb != 0 && crossesBorder(coord(i, j), coord(i + 1, j))) {
                upProb = upProb * borderCrossingP;
            }
            if (downProb != 0 && crossesBorder(coord(i, j), coord(i - 1, j))) {
                downProb = downProb * borderCrossingP;
            }
            if (leftProb != 0 && crossesBorder(coord(i, j), coord(i, j - 1))) {
                leftProb = leftProb * borderCrossingP;
            }
            if (rightProb != 0 && crossesBorder(coord(i, j), coord(i, j + 1))) {
                rightProb = rightProb * borderCrossingP;
            }

            //1d. calculate remainder
            noneProb = 1 - leftProb - rightProb - upProb - downProb;
            //2. Pull from multinomial for each population
            int multiPull[5] = {};
            double mProbs[5] = { leftProb, rightProb, upProb, downProb, noneProb };
            
            //2a. Uninfected male
            if (lattice[i][j].numPsyllids_male > 0) {
                std::fill(multiPull, multiPull + 5, 0);
                multinom(5, lattice[i][j].numPsyllids_male, mProbs, multiPull);
                assert(accumulate(multiPull, multiPull + 5, 0) == lattice[i][j].numPsyllids_male);
                if (leftProb > 0) {
                    diffMales[numRows*(j-1) + i] += multiPull[0];
                }
                if (rightProb > 0) {
                    //lattice[i][j + 1].numPsyllids_male += multiPull[1];
                    diffMales[numRows*(j+1) + i] += multiPull[1];
                }
                if (upProb > 0) {
                    //lattice[i + 1][j].numPsyllids_male += multiPull[2];
                    diffMales[numRows*j + i+1] += multiPull[2];
                }
                if (downProb > 0) {
                    //lattice[i - 1][j].numPsyllids_male += multiPull[3];
                    diffMales[numRows*j + i-1] += multiPull[3];
                }
                //lattice[i][j].numPsyllids_male -= (lattice_old[i][j].numPsyllids_male - multiPull[4]);
                diffMales[numRows*j + i] -= (lattice[i][j].numPsyllids_male - multiPull[4]);
            }

            //2b. Uninfected female
            if (lattice[i][j].numPsyllids_female > 0) {
                std::fill(multiPull, multiPull + 5, 0);
                multinom(5, lattice[i][j].numPsyllids_female, mProbs, multiPull);
                assert(accumulate(multiPull, multiPull + 5, 0) == lattice[i][j].numPsyllids_female);
                if (leftProb > 0) {
                    //lattice[i][j - 1].numPsyllids_female += multiPull[0];
                    diffFemales[numRows*(j-1) + i] += multiPull[0];
                }
                if (rightProb > 0) {
                    //lattice[i][j + 1].numPsyllids_female += multiPull[1];
                    diffFemales[numRows*(j+1) + i] += multiPull[1];
                }
                if (upProb > 0) {
                    //lattice[i + 1][j].numPsyllids_female += multiPull[2];
                    diffFemales[numRows*j + i+1] += multiPull[2];
                }
                if (downProb > 0) {
                    //lattice[i - 1][j].numPsyllids_female += multiPull[3];
                    diffFemales[numRows*j + i-1] += multiPull[3];
                }
                //lattice[i][j].numPsyllids_female -= (lattice_old[i][j].numPsyllids_female - multiPull[4]);
                diffFemales[numRows*j + i] -= (lattice[i][j].numPsyllids_female - multiPull[4]);
            }
            //2c. Infected male
            if (lattice[i][j].numInfectedPsyllids_male > 0) {
                std::fill(multiPull, multiPull + 5, 0);
                multinom(5, lattice[i][j].numInfectedPsyllids_male, mProbs, multiPull);
                assert(accumulate(multiPull, multiPull + 5, 0) == lattice[i][j].numInfectedPsyllids_male);
                if (leftProb > 0) {
                    //lattice[i][j - 1].numInfectedPsyllids_male += multiPull[0];
                    diffMales_i[numRows*(j-1) + i] += multiPull[0];
                }
                if (rightProb > 0) {
                    //lattice[i][j + 1].numInfectedPsyllids_male += multiPull[1];
                    diffMales_i[numRows*(j+1) + i] += multiPull[1];
                }
                if (upProb > 0) {
                    //lattice[i + 1][j].numInfectedPsyllids_male += multiPull[2];
                    diffMales_i[numRows*j + i+1] += multiPull[2];
                }
                if (downProb > 0) {
                    //lattice[i - 1][j].numInfectedPsyllids_male += multiPull[3];
                    diffMales_i[numRows*j + i-1] += multiPull[3];
                }
                //lattice[i][j].numInfectedPsyllids_male -= (lattice_old[i][j].numInfectedPsyllids_male - multiPull[4]);
                diffMales_i[numRows*j + i] -= (lattice[i][j].numInfectedPsyllids_male - multiPull[4]);
            }

            //2d. Infected female
            if (lattice[i][j].numInfectedPsyllids_female > 0) {
                std::fill(multiPull, multiPull + 5, 0);
                multinom(5, lattice[i][j].numInfectedPsyllids_female, mProbs, multiPull);
                assert(accumulate(multiPull, multiPull + 5, 0) == lattice[i][j].numInfectedPsyllids_female);
                if (leftProb > 0) {
                    //lattice[i][j - 1].numInfectedPsyllids_female += multiPull[0];
                    diffFemales_i[numRows*(j-1) + i] += multiPull[0];
                }
                if (rightProb > 0) {
                    //lattice[i][j + 1].numInfectedPsyllids_female += multiPull[1];
                    diffFemales_i[numRows*(j+1) + i] += multiPull[1];
                }
                if (upProb > 0) {
                    //lattice[i + 1][j].numInfectedPsyllids_female += multiPull[2];
                    diffFemales_i[numRows*j + i+1] += multiPull[2];
                }
                if (downProb > 0) {
                    //lattice[i - 1][j].numInfectedPsyllids_female += multiPull[3];
                    diffFemales_i[numRows*j + i-1] += multiPull[3];
                }
                //lattice[i][j].numInfectedPsyllids_female -= (lattice_old[i][j].numInfectedPsyllids_female - multiPull[4]);
                diffFemales_i[numRows*j + i] -= (lattice[i][j].numInfectedPsyllids_female - multiPull[4]);
            }
        }
    }

    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < rowLength; j++) {
            lattice[i][j].numPsyllids_male += diffMales[numRows*j + i];
            lattice[i][j].numPsyllids_female += diffFemales[numRows*j + i];
            lattice[i][j].numInfectedPsyllids_male += diffMales_i[numRows*j + i];
            lattice[i][j].numInfectedPsyllids_female += diffFemales_i[numRows*j + i];
        }
    }
    delete diffMales;
    delete diffFemales;
    delete diffMales_i;
    delete diffFemales_i;
    
}

/**************************************************
 * Psyllid aging
 * Psyllids age and become adults:
 * - Calculate nymph/adult survival probabilities
 *   based on current carrying capacity
 * - Age each subpopulation of adult psyllids
 * - Age each subpopulation of nymphs
 * ************************************************/
void psyllidAging() {
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < rowLength; j++) {
            //Calculate modifier
            double modifier = 1;
            if (lattice[i][j].getTotalPsyllids() > carryingCapacity) {
                modifier = carryingCapacity / (double)lattice[i][j].getTotalPsyllids();
            }
            double nymphSurvivalChance = eggSurvivalP * modifier;
            double adultSurvivalChance = adultSurvivalP * modifier;
            int totalPsyllids = lattice[i][j].getTotalPsyllids();

            
            //Age adults
            
            //UI Male
            int ui_male = lattice[i][j].numPsyllids_male;
            if (ui_male > 0) {
                boost::random::binomial_distribution<int> pAdult(ui_male, adultSurvivalChance);
                lattice[i][j].numPsyllids_male = pAdult(rng);
            }
            //I Male
            int i_male = lattice[i][j].numInfectedPsyllids_male;
            if (i_male > 0) {
                boost::random::binomial_distribution<int> pAdult(i_male, adultSurvivalChance);
                lattice[i][j].numInfectedPsyllids_male = pAdult(rng);
            }
            
            //UI Female
            int ui_female = lattice[i][j].numPsyllids_female;
            if (ui_female > 0) {
                boost::random::binomial_distribution<int> pAdult(ui_female, adultSurvivalChance);
                lattice[i][j].numPsyllids_female = pAdult(rng);
            }
            //I Female
            int i_female = lattice[i][j].numInfectedPsyllids_female;
            if (i_female > 0) {
                boost::random::binomial_distribution<int> pAdult(i_female, adultSurvivalChance);
                lattice[i][j].numInfectedPsyllids_female = pAdult(rng);
            }
            //Age nymphs
            for (int k = 16; k >= 0; k--) {
                //Max age nymphs
                if (k == 16) {
                    // UNINFECTED
                    //How many age up
                    int startingNymphs_ui = lattice[i][j].numNymphs[k];
                    if (startingNymphs_ui > 0) {
                        boost::random::binomial_distribution<int> pAge_ui(startingNymphs_ui, nymphSurvivalChance);
                        int numNymphsAging_ui = pAge_ui(rng);
                        //How many are female
                        boost::random::binomial_distribution<int> pSex_ui(numNymphsAging_ui, 0.5);
                        int numFemaleNymphs_ui = pSex_ui(rng);
                        lattice[i][j].numPsyllids_female += numFemaleNymphs_ui;
                        lattice[i][j].numPsyllids_male += numNymphsAging_ui - numFemaleNymphs_ui;
                    }
                    // INFECTED
                    //How many age up
                    int startingNymphs_i = lattice[i][j].numInfectedNymphs[k];
                    if (startingNymphs_i > 0) {
                        boost::random::binomial_distribution<int> pAge_i(startingNymphs_i, nymphSurvivalChance);
                        int numNymphsAging_i = pAge_i(rng);
                        //How many are female
                        boost::random::binomial_distribution<int> pSex_i(numNymphsAging_i, 0.5);
                        int numFemaleNymphs_i = pSex_i(rng);
                        lattice[i][j].numInfectedPsyllids_female += numFemaleNymphs_i;
                        lattice[i][j].numInfectedPsyllids_male += numNymphsAging_i - numFemaleNymphs_i;
                    }
                }



                if (k == 0) {
                    lattice[i][j].numNymphs[k] = 0;
                    lattice[i][j].numInfectedNymphs[k] = 0;
                }
                else {
                    //UNINFECTED
                    lattice[i][j].numNymphs[k] = 0;
                    if (lattice[i][j].numNymphs[k - 1] > 0) {
                        boost::random::binomial_distribution<int> pAge_ui(lattice[i][j].numNymphs[k - 1], nymphSurvivalChance);
                        lattice[i][j].numNymphs[k] = pAge_ui(rng);
                    }
                    //INFECTED
                    lattice[i][j].numInfectedNymphs[k] = 0;
                    if (lattice[i][j].numInfectedNymphs[k - 1] > 0) {
                        boost::random::binomial_distribution<int> pAge_i(lattice[i][j].numInfectedNymphs[k - 1], nymphSurvivalChance);
                        lattice[i][j].numInfectedNymphs[k] = pAge_i(rng);
                    }
                }
            }


        }
    }
}

/*****************************************************
 * Egg management
 * Egg laying:
 * - Viable mothers lay eggs based on model state
 * ****************************************************/
void eggManagement() {
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < rowLength; j++) {
            int numMothers = lattice[i][j].numPsyllids_female + lattice[i][j].numInfectedPsyllids_female;
            int numViableShoots = accumulate(lattice[i][j].numShoots.begin(), lattice[i][j].numShoots.end(), 0);
            numViableShoots += accumulate(lattice[i][j].numInfectedShoots.begin(), lattice[i][j].numInfectedShoots.end(), 0);
            int totalNymphs = min(numMothers * eggsPerFemaleAdult, numViableShoots * shootEggCapacity);
            lattice[i][j].numNymphs[0] += totalNymphs;
        }
    }
}

/**************************************************
 * Disease transmission
 * Transmission of HLB between ACP and flush: 
 * - Transmit HLB from flush -> psyllids
 * - Transmit HLB from psyllids -> flush
 * ************************************************/
void diseaseTransmission() {
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < rowLength; j++) {
            //1. Transmit from flush to psyllids
            int numUninfectedShoots = accumulate(lattice[i][j].numShoots.begin(), lattice[i][j].numShoots.end(), 0);
            int numInfectedShoots = accumulate(lattice[i][j].numInfectedShoots.begin(), lattice[i][j].numInfectedShoots.end(), 0);
            if (lattice[i][j].infectious) {
                double infectedproportion = 0;
                if (numUninfectedShoots > 0) {
                    infectedproportion = (double)numInfectedShoots / ((double)numInfectedShoots + (double)numUninfectedShoots);
                }
                for (int k = nymphMinAgeToBeInfected; k < 17; k++) {
                    int ui_nymphs = lattice[i][j].numNymphs[k];
                    ui_nymphs = floor(ui_nymphs * infectedproportion);
                    int nymphsInfected = floor(ui_nymphs * transmissionFlushNymph);
                    lattice[i][j].numNymphs[k] -= nymphsInfected;
                    lattice[i][j].numInfectedNymphs[k] += nymphsInfected;
                }
            }

            //2. Transmit from psyllids to flush
            if (numUninfectedShoots != 0) {
                int totalShoots = numUninfectedShoots + numInfectedShoots;
                double probabilities[18] = {};
                double finalProb = 1;
                for (int k = 0; k < 17; k++) {
                    probabilities[k] = (double)lattice[i][j].numShoots[k] / (double)totalShoots * (double)transmissionAdultFlush;
                    finalProb -= probabilities[k];
                }
                probabilities[17] = finalProb;
                int deltaInfected[18] = {};
                int totalInfectedPsyllids = lattice[i][j].numInfectedPsyllids_female + lattice[i][j].numInfectedPsyllids_male;
                multinom(18, totalInfectedPsyllids, probabilities, deltaInfected);
                for (int k = 0; k < 17; k++) {
                    lattice[i][j].numShoots[k] -= min(deltaInfected[k], lattice[i][j].numShoots[k]);
                    lattice[i][j].numInfectedShoots[k] += min(deltaInfected[k], lattice[i][j].numShoots[k]);
                }
            }
        }
    }
}

/***************************************************
 * Age Flush
 * Flush aging:
 * - Flush age, those reaching maximum age become
 *   hardened
 * *************************************************/
void ageFlush() {
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < rowLength; j++) {
            for (int k = 29; k >= 0; k--) {
                //Oldest shoots graduate to old shoots
                if (k == 29) {
                    lattice[i][j].oldUninfectedShoots += lattice[i][j].numShoots[k];
                    lattice[i][j].oldInfectedShoots += lattice[i][j].numInfectedShoots[k];
                }

                // All ages shift forward 1
                if (k == 0) {
                    lattice[i][j].numShoots[k] = 0;
                    lattice[i][j].numInfectedShoots[k] = 0;
                }
                else {
                    lattice[i][j].numShoots[k] = lattice[i][j].numShoots[k - 1];
                    lattice[i][j].numInfectedShoots[k] = lattice[i][j].numInfectedShoots[k - 1];
                }
            }

            //Infection management
            
            if (lattice[i][j].getHLBSeverity() > 0) {
                if (lattice[i][j].daysInfected <= asymptomaticLength) {
                    lattice[i][j].daysInfected++;
                    if (lattice[i][j].daysInfected == latentPeriod) {
                        lattice[i][j].infectious = true;
                    }
                    if (lattice[i][j].daysInfected == asymptomaticLength) {
                        lattice[i][j].symptomatic = true;
                    }
                }
            }
        }
    }
}

#pragma endregion


/******************************************
 * Advance biological model
 * The wrapper function for all activities,
 * called by the econ layer
 * ****************************************/
void advanceBiologicalModel() {
    //First time function is called
    if (!modelStarted) {
        initializeModel();
        if (outputFlag) {
            initializeCSV();
        }
        modelDay = 0;
        modelStarted = true;
    }
    //Invasion day activities
    if (!invasionDays_q.empty() && modelDay == invasionDays_q.front()) {
        placeInitialPsyllids(invasionModalities_q.front(), invasionGrove);
        invasionDays_q.pop();
        invasionModalities_q.pop();
    }

    //Keep track of flushing periods
    setFlushingPeriod(modelDay);

    //Main activities
    if (isFlushingPeriod) {
       birthNewFlush();
    }

    migration();
    eggManagement();

    diseaseTransmission();

    ageFlush();

    psyllidAging();

    if (outputFlag) {
        write_csv_batch();
    }
    modelDay++;
}






} //end namespace
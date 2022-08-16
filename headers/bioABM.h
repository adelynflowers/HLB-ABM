#ifndef BIO_ABM
#define BIO_ABM
/********************************
 * INCLUDES
 * ******************************/
#include <iostream>
#include <vector>
#include <queue>
#include <array>
#include <boost/random.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/assign.hpp>
#include <memory>
#include <list>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <string>
#include <algorithm>
#include <cereal/archives/json.hpp>
#include <cereal/types/memory.hpp>
#include <ctime>
#include <chrono>
#include <atomic>
#include <random>
#include <numeric>
#include <stdlib.h>
#include <iostream>


namespace bioABM {
    int getModelDay();
    int getModelDuration();

    int getNumRows();
    int getRowLength();
    
    // Returns age at coordinates, with an optional age differential in days
    int getAgeAt(int i, int j, int differential = 0);

    //returns an array of flush starting days (spring, summer, fall)
    int getFallStart();
    int getSummerStart();
    int getSpringStart();
    int countPsyllids();
    //Returns severity at coordinates
    double getSeverityAt(int i, int j);

    //Returns true if tree is symptomatic
    bool isSymptomatic(int i, int j);

    //Returns psyllids at coordinates
    double getPsyllidsAt(int i, int j);

    typedef boost::tuple<int, int> coord;
    //Functions to interface w/ from economic model
    void rogueTreeAt(int i, int j);

    void sprayTrees(double, std::vector<coord>);

    bool isValidCoordinate(coord);

    void parseParameterFile(std::string);

    void advanceBiologicalModel();
    
    //void advanceBiologicalModel_parallel();

    bool isTreeAlive(int i, int j);

    void finishRun();

    void bioTestSuite();
    
    void setExperimentID(int);
}

















#endif
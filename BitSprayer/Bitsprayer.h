#ifndef MEDICALPREDICTOR_BITSPRAYER_H
#define MEDICALPREDICTOR_BITSPRAYER_H

#include <cmath>
#include <list>
#include <sstream>
#include <vector>
#include <iostream>
#include <time.h>
#include <fstream>

using namespace std;

class Bitsprayer {
  public:
    Bitsprayer();           //creates an unallocated bitspray
    Bitsprayer(int states);      //create a bitspray with buffer S states
    Bitsprayer(Bitsprayer &other);  //copy constructor
    ~Bitsprayer();                //destructor

    int create(int states);
    int randomize();
    int copy(Bitsprayer &other);
    int print();
    int print(ostream &aus);
    static int destroy();

    int twoPtCrossover(Bitsprayer &other);
    int mutate(int numMuts);

    vector<int> getBitsVec(int len);

  private:
    int initInput;
    int numStates;
    int initState;
    int curState;
    vector<int> buf;
    vector<vector<int> > transitions;
    vector<vector<vector<int> > > responses;
};

#endif //MEDICALPREDICTOR_BITSPRAYER_H

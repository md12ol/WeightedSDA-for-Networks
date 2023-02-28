#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

class Graph {
public:
    Graph();
    explicit Graph(int nn);

    int fill(vector<int> &weights, bool diag);
    vector<int> SIR(double alpha, int p0);
    void print(ostream &out); 
    vector<int> weightHist();

protected:
    static int infect(int nin, double alpha);
    int numNodes;
    int numEdges;
    int totWeight;
    int maxWeight;
    vector<vector<int> > adjM;
};


#endif //THADSNBIT_GRAPH_H

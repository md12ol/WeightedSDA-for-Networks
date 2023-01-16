#include "Bitsprayer.h"

#define VERBOSEBIT false
#define MAXVAL 6

Bitsprayer::Bitsprayer() {

    initInput = -1;
    numStates = -1;
    initState = -1;
    curState = -1;
//    transitions = {{}};
//    responses = {{{}}};
    if (VERBOSEBIT) cout << "Bitsprayer Made" << endl;
}

Bitsprayer::Bitsprayer(int states) {

    initInput = -1;
    numStates = states;
    initState = -1;
    curState = -1;
    transitions.reserve(states);
    for (vector<int> v: transitions){
        v.reserve(MAXVAL);
    }
    responses.reserve(states);
    for (vector<vector<int> > v: responses){
        v.reserve(MAXVAL);
    }

//    transitions = {{}};
//    responses = {{{}}};
    if (VERBOSEBIT) cout << "Bitsprayer Made w " << states << " states" << endl;
}

Bitsprayer::Bitsprayer(Bitsprayer &other) {
    copy(other);
}

Bitsprayer::~Bitsprayer() {
    destroy();
}

int Bitsprayer::create(int states) {
    initInput = (int)lrand48() % MAXVAL;
    numStates = states;
    initState = 0;
    curState = -1;

    vector<int> oneState;
    for (int s = 0; s < numStates; ++s) {
        oneState.clear();
        for (int i = 0; i < MAXVAL; ++i) {
            oneState.push_back((int) lrand48() % numStates);
        }
        transitions.push_back(oneState);
    }

    vector<int> oneResponse;
    vector <vector<int>> oneStateResps;
    int respSize;
    for (int s = 0; s < numStates; ++s) {
        oneStateResps.clear();
        for (int t = 0; t < MAXVAL; ++t) {
            oneResponse.clear();
            respSize = (int)lrand48() % 2 + 1;
            for (int i = 0; i < respSize; ++i) {
                oneResponse.push_back((int)lrand48() % MAXVAL);
            }
            oneStateResps.push_back(oneResponse);
        }
        responses.push_back(oneStateResps);
    }
    if (VERBOSEBIT) cout << "Bitsprayer Created w " << states << " states" << endl;
    return 0;
}

int Bitsprayer::randomize() {
    initInput = (int) lrand48() % 2;
    // TODO: Change start state?
//    initState = (int) lrand48()%numStates;
    curState = -1;

    vector<int> oneResponse;
    int respSize;
    for (int s = 0; s < numStates; ++s) {
        transitions.at(s).at(0) = (int)lrand48() % numStates;
        transitions.at(s).at(1) = (int)lrand48() % numStates;
        for (int t = 0; t < 2; ++t) {
            oneResponse.clear();
            respSize = (int)lrand48() % 2 + 1;
            for (int i = 0; i < respSize; ++i) {
                oneResponse.push_back((int)lrand48() % MAXVAL);
            }
            responses.at(s).at(t) = oneResponse;
        }
    }
    return 0;
}

int Bitsprayer::copy(Bitsprayer &other) {
    initInput = other.initInput;
    numStates = other.numStates;
    initInput = other.initState;
    curState = other.curState;

    transitions = other.transitions;
    responses = other.responses;
    if (VERBOSEBIT) cout << "Bitsprayer Copied" << endl;
    return 0;
}

int Bitsprayer::print() {
    print(cout);
    return 0;
}

int Bitsprayer::print(ostream &aus) {
    aus << initState << " <- " << initInput << endl;
    for (int s = 0; s < numStates; ++s) {
        for (int t = 0; t < MAXVAL; ++t) {
            aus << s << " -> " << transitions.at(s).at(t) << " [";
            for (int v: responses.at(s).at(t)) {
                aus << " " << v;
            }
            aus << " ]" << endl;
        }
    }
    if (VERBOSEBIT) cout << "Bitsprayer Printed" << endl;
    return 0;
}

int Bitsprayer::destroy(){
    return 0;
}

int Bitsprayer::twoPtCrossover(Bitsprayer &other) {
    int cp1, cp2;
    int swapInt;
    vector<int> swapVec;

    if (numStates != other.numStates){
        return 1;
    }

    do {
        cp1 = (int) lrand48() % numStates;
        cp2 = (int) lrand48() % numStates;
        if (cp1 > cp2) {
            swapInt = cp1;
            cp1 = cp2;
            cp2 = swapInt;
        }
    } while (cp1==cp2);

    if (cp1 == 0) {
        swapInt = initInput;
        initInput = other.initInput;
        other.initInput = swapInt;
    }

    for (int s = cp1; s < cp2; s++) {
        swapVec = transitions.at(s);
        transitions.at(s) = other.transitions.at(s);
        other.transitions.at(s) = swapVec;
        swapVec = responses.at(s).at(0);
        responses.at(s).at(0) = other.responses.at(s).at(0);
        other.responses.at(s).at(0) = swapVec;
        swapVec = responses.at(s).at(1);
        responses.at(s).at(1) = other.responses.at(s).at(1);
        other.responses.at(s).at(1) = swapVec;
    }
    return 0;
}

int Bitsprayer::mutate(int numMuts){
    int mutPt;
    vector<int> oneResponse;
    int respSize;

    for (int m = 0; m < numMuts; ++m) {
        mutPt = (int) lrand48() % (2 * numStates + 1);

        if (mutPt == 0){
            initInput = 1-initInput;
            return 0;
        }
        mutPt = (mutPt - 1)/2;
        int transNum = (int) lrand48() % MAXVAL;
        if ((int)lrand48() % 2 == 0){ // Mutate transition
            transitions.at(mutPt).at(transNum) = (int) lrand48() % numStates;
        } else { // Mutate response
            oneResponse.clear();
            respSize = (int) lrand48() % 2 + 1;
            for (int i = 0; i < respSize; ++i) {
                oneResponse.push_back((int) lrand48() % MAXVAL);
            }
            responses.at(mutPt).at(transNum) = oneResponse;
        }
    }
    return 0;
}

vector<int> Bitsprayer::getBitsVec(int len){
    vector<int> rtn;
    int nextBit;
    curState = initState;
    buf.clear();
    rtn.push_back(initInput);
    buf.push_back(initInput);

    while (rtn.size() < len){
        nextBit = buf.front();
        buf.erase(buf.begin());
        for(int i: responses.at(curState).at(nextBit)){
            buf.push_back(i);
            rtn.push_back(i);
        }
        curState = transitions.at(curState).at(nextBit);
    }
    return rtn;
}

int printVec(vector<int> vec){
    for (int i:vec) {
        cout<<i<<" ";
    }
    cout<<endl;
    return 0;
}

//int main() {
//    Bitsprayer b(5);
//    b.create(5);
//    b.randomize();
//    b.print();
//
////    Bitsprayer b2(5);
////    b2.create(5);
////    b2.randomize();
////    b2.print();
////
////    b.twoPtCrossover(b2);
////    b.print();
////    b2.print();
//
//    printVec(b.getBitsVec(20));
//
//    if (VERBOSEBIT) cout << "DONE!" << endl;
//    return 0;
//}
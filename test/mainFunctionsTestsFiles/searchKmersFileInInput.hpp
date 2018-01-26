#include <iostream>
#include <fstream>
#include <thread>
#include <chrono>
#include <string>
#include <vector>

#include "../functionTestsFiles/utilitiesTests.hpp"

void searchKmers( vector<string>* kmers, string searchFileName, struct KmersStats *ks);
int main(int argc, char const *argv[]);
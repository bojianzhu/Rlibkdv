#include "alg_visual.h"
#include <Rcpp.h>
using namespace Rcpp;


//' @title kernel density visualization in C++
//' @name kdvCpp
//' @param args arguments for kdv
//' @return the kdv result
//' @keywords internal
// [[Rcpp::export]]
std::string kdvCpp(CharacterVector args){
    int argc = args.size() + 1;
    char **argv = new char*[argc];
    // argv[0] = "kdv";
    for (int i = 0; i < args.size(); i++) {
        argv[i + 1] = args[i];
    }
    alg_visual algorithm;
    algorithm.load_datasets_CSV(argv);
	std::string result = algorithm.compute(argc, argv);
    delete[] argv;
	algorithm.clear_basic_memory();
    return result;

}

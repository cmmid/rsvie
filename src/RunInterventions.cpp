#include <Rcpp.h>

using namespace Rcpp;
// [[Rcpp::plugins("cpp14")]]
// [[Rcpp::depends(BH)]]    

#include "RunInterventions.h"


RCPP_MODULE(RunInterventionsModule) {
    class_<RunInterventions>( "RunInterventions" )
    .constructor<double, double, NumericVector>()
    .field( "p_mat", &RunInterventions::p_mat )
    .field( "pVHR", &RunInterventions::pVHR )
    .field( "pHR", &RunInterventions::pHR )
    .field( "pLR", &RunInterventions::pLR )
    
    .field( "cnt_matrix_p", &RunInterventions::cnt_matrix_p )
    .field( "cnt_matrix_p_h", &RunInterventions::cnt_matrix_p_h )
    .field( "pwp_p", &RunInterventions::pwp_p )
    .field( "pwn_p", &RunInterventions::pwn_p )
    .field( "nwp_p", &RunInterventions::nwp_p )
    .field( "nwn_p", &RunInterventions::nwn_p )

    .field( "cnt_matrix_c", &RunInterventions::cnt_matrix_c )
    .field( "cnt_matrix_c_h", &RunInterventions::cnt_matrix_c_h )
    .field( "pwp_c", &RunInterventions::pwp_c )
    .field( "pwn_c", &RunInterventions::pwn_c )
    .field( "nwp_c", &RunInterventions::nwp_c )
    .field( "nwn_c", &RunInterventions::nwn_c )
    
    .field( "u_p", &RunInterventions::u_p )
    
    .field( "run_start", &RunInterventions::run_start )
    .field( "run_burn", &RunInterventions::run_burn )
    .field( "run_full", &RunInterventions::run_full )

    .method( "Sample", &RunInterventions::Sample )
    .method( "StatesValues", &RunInterventions::StatesValues )


    ;
}

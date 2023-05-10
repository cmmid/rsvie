//
//  RunIntervention.h
//  
//
//  Created by David Hodgson on 05/07/2021.
//

#ifndef RunIntervention_h
#define RunIntervention_h

#include <Rcpp.h>
#include <RcppEigen.h>
#include <random>
#include <boost/random.hpp>
#include <boost/math/distributions.hpp>
#include "ascent/Ascent.h"

//#include "manipContacts.h"

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins("cpp14")]]
using namespace Rcpp;
using namespace std;
using namespace Eigen;
using namespace asc;

// Stuff for random number generation
std::random_device dev;
std::mt19937 engine(dev());
typedef boost::mt19937 PRNG_s;
PRNG_s rng(engine());

class RunInterventions {
public:
    int A;

    vector<double > p_mat;
    vector<double > u_p;

    MatrixXd cnt_matrix_p;
    MatrixXd cnt_matrix_p_h;
    MatrixXd pwp_p;
    MatrixXd pwn_p;
    MatrixXd nwp_p;
    MatrixXd nwn_p;
    
    MatrixXd cnt_matrix_c;
    MatrixXd cnt_matrix_c_h;
    MatrixXd pwp_c;
    MatrixXd pwn_c;
    MatrixXd nwp_c;
    MatrixXd nwn_c;


    int rg;
    int sg;
    int ag;
    
    List Sample(List vac_calendar, List vac_dose, double cov_c, List vac_info, VectorXd posteriors);
    MatrixXd StatesValues(List vac_calendar, List vac_dose, double cov_c, List vac_info, VectorXd posteriors);

    /**********************************************************************/
    /** CwX Contacts for persons who are cocooned ****/
    /**********************************************************************/
    
    vector<vector<double > > get_cwn(double prop_c, char s)
    {
        vector<vector<double > > cwn_e(A,vector<double>(A, 0));
        if (s == 'p')
        {
            for (int i = 0; i < 12; i++)        //infants with
                for (int j = 0; j < A; j++)        //infants
                    cwn_e[i][j] = (nwn_p(i, j));
            
            for (int i = 0; i < 12; i++)
                for (int j = 0; j < 12; j++)        //infants
                    cwn_e[i][j] = (nwn_p(i, j))*(1-prop_c); //contacts made with cocooned infants
            
            for (int i = 18; i < 21; i++)        //infants with
                for (int j = 0; j < A; j++)        //infants
                    cwn_e[i][j] = (pwn_p(i, j));
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 0; j < 12; j++) // all
                    cwn_e[i][j] = (cnt_matrix_p(i, j) - cnt_matrix_p_h(i, j))*(1-prop_c);

        }
        else if (s == 'c')
        {
            for (int i = 0; i < 12; i++)        //infants with
                for (int j = 0; j < A; j++)        //infants
                    cwn_e[i][j] = (nwn_c(i, j));
            
            for (int i = 0; i < 12; i++)
                for (int j = 0; j < 12; j++)        //infants
                    cwn_e[i][j] = (nwn_c(i, j))*(1-prop_c); //contacts made with cocooned infants
            
            for (int i = 18; i < 21; i++)        //infants with
                for (int j = 0; j < A; j++)        //infants
                    cwn_e[i][j] = (pwn_c(i, j));
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 0; j < 12; j++) // all
                    cwn_e[i][j] = (cnt_matrix_c(i, j) - cnt_matrix_c_h(i, j))*(1-prop_c);
        }
        else
            cout << "Error cont" << endl;
        
        return cwn_e;
    }
    
    vector<vector<double > > get_cwp(double prop_c, char s)
    {
        vector<vector<double > > cwp_e(A,vector<double>(A,0));
        
        if (s == 'p')
        {
            for (int i = 0; i < 12; i++) //infants with/
                for (int j = 18; j < 21; j++) //cba
                    cwp_e[i][j] = (cnt_matrix_p(i, j) - cnt_matrix_p_h(i, j))*(p_mat[j])*(1-prop_c);
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 18; j < 21; j++) // infants
                    cwp_e[i][j] = (pwp_p(i, j))*(1-prop_c); // parental w. infant plus non-parental with infant who are cocooned
    
        }
        else if (s == 'c')
        {
            for (int i = 0; i < 12; i++) //infants with/
                for (int j = 18; j < 21; j++) //cba
                    cwp_e[i][j] = (cnt_matrix_c(i, j) - cnt_matrix_c_h(i, j))*(p_mat[j])*(1-prop_c);
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 18; j < 21; j++) // infants
                    cwp_e[i][j] = (pwp_c(i, j))*(1-prop_c); // parental w. infant plus non-parental with infant who are cocooned
        }
        else
            cout << "Error cont" << endl;
        
        return cwp_e;
    }
    
    vector<vector<double > > get_cwc(double prop_c, char s)
    {
        vector<vector<double > > cwc_e(A,vector<double>(A,0));

        if (s == 'p')
        {
            for (int i = 0; i < 12; i++)
                for (int j = 0; j < 12; j++)        //infants
                    cwc_e[i][j] = (nwn_p(i, j))*(prop_c); //contacts made with cocooned infants
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 0; j < 12; j++) // all
                    cwc_e[i][j] = cnt_matrix_p_h(i, j) + (cnt_matrix_p(i, j) - cnt_matrix_p_h(i, j))*(prop_c);
            
            for (int i = 0; i < 12; i++) //infants with/
                for (int j = 18; j < 21; j++) //cba
                    cwc_e[i][j] = cnt_matrix_p_h(i, j)/2.0 + (cnt_matrix_p(i, j) - cnt_matrix_p_h(i, j))*(p_mat[j])*(prop_c);
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 18; j < 21; j++) // infants
                    cwc_e[i][j] = (pwp_p(i, j))*(prop_c); // parental w. infant plus non-parental with infant who are cocooned
    
        }
        else if (s == 'c')
        {
            for (int i = 0; i < 12; i++)
                for (int j = 0; j < 12; j++)        //infants
                    cwc_e[i][j] = (nwn_c(i, j))*(prop_c); //contacts made with cocooned infants
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 0; j < 12; j++) // all
                    cwc_e[i][j] = cnt_matrix_c_h(i, j) + (cnt_matrix_c(i, j) - cnt_matrix_c_h(i, j))*(prop_c);
            
            for (int i = 0; i < 12; i++) //infants with/
                for (int j = 18; j < 21; j++) //cba
                    cwc_e[i][j] = cnt_matrix_c_h(i, j)/2.0 + (cnt_matrix_c(i, j) - cnt_matrix_c_h(i, j))*(p_mat[j])*(prop_c);
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 18; j < 21; j++) // infants
                    cwc_e[i][j] = (pwp_c(i, j))*(prop_c); // parental w. infant plus non-parental with infant who are cocooned

        }
        else
            cout << "Error cont" << endl;
        return cwc_e;
    }
    
    /**********************************************************************/
    /** PwX Contacts for persons who are mothers but not cocooned ****/
    /**********************************************************************/
    
    
    vector<vector<double > > get_pwn(double prop_c, char s)
    {
        vector<vector<double > > pwn_e(A,vector<double>(A,0));
        
        if (s == 'p')
        {
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 0; j < A; j++) // all
                    pwn_e[i][j] = pwn_p(i, j); // all contacts
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 0; j < 12; j++) // all
                    pwn_e[i][j] = cnt_matrix_p_h(i, j) + (cnt_matrix_p(i, j) - cnt_matrix_p_h(i, j))*(1-prop_c);
        }
        else if (s == 'c')
        {
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 0; j < A; j++) // all
                    pwn_e[i][j] = pwn_c(i, j); // all contacts
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 0; j < 12; j++) // all
                    pwn_e[i][j] = cnt_matrix_c_h(i, j) + (cnt_matrix_c(i, j) - cnt_matrix_c_h(i, j))*(1-prop_c);
        }
        else
            cout << "Error cont" << endl;
        
        return pwn_e;
    }
    
    vector<vector<double > > get_pwp(double prop_c, char s)
    {
        vector<vector<double > > pwp_e(A,vector<double>(A,0));
        
        if (s == 'p')
        {
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 18; j < 21; j++) // all
                    pwp_e[i][j] = pwp_p(i, j)*(1-prop_c); // all contacts
        }
        else if (s == 'c')
        {
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 18; j < 21; j++) // all
                    pwp_e[i][j] = pwp_c(i, j)*(1-prop_c); // all contacts
        }
        else
            cout << "Error cont" << endl;
        
        return pwp_e;
    }
    
    vector<vector<double > > get_pwc(double prop_c, char s)
    {
        // Non cocooned parent -> cocooned infant (family and outside family)
        vector<vector<double > > pwc_e(A,vector<double>(A,0));
        
        if (s == 'p')
        {
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 0; j < 12; j++) // all
                    pwc_e[i][j] = (cnt_matrix_p(i, j) - cnt_matrix_p_h(i, j))*(prop_c);
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 18; j < 21; j++) // infants
                    pwc_e[i][j] = pwn_p(i, j)*prop_c; //
        }
        else if (s == 'c')
        {
            // Parent -> parents (all)
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 0; j < 12; j++) // all
                    pwc_e[i][j] = (cnt_matrix_c(i, j) - cnt_matrix_c_h(i, j))*(prop_c);
            
            for (int i = 18; i < 21; i++) // cba with.
                for (int j = 18; j < 21; j++) // infants
                    pwc_e[i][j] = pwn_c(i, j)*prop_c; //
        }
        else
            cout << "Error cont" << endl;
        
        return pwc_e;
    }
    
    
 /**********************************************************************/
/** NwX Contacts for persons who are neither cocooned nor mothers ****/
 /**********************************************************************/
    
    vector<vector<double > > get_nwn(double prop_c, char s)
    {
        // cnt_matrix_p cnt_matrix_p_h
        vector<vector<double > > nwn_e(A,vector<double>(A,0));
        if (s == 'p')
        {
            for (int i = 0; i < A; i++) // population with.
                for (int j = 0; j < A; j++) //     infants
                    nwn_e[i][j] = nwn_p(i, j); //number of contacts made by population with non-cocooned prop_ci
            
            for (int i = 0; i < A; i++) // population with.
                for (int j = 0; j < 12; j++) //     infants
                    nwn_e[i][j] = nwn_p(i, j)*(1-prop_c); //number
        }
        else if (s == 'c')
        {
            for (int i = 0; i < A; i++) // population with.
                for (int j = 0; j < A; j++) //     infants
                    nwn_e[i][j] = nwn_c(i, j); //number of contacts made by population with non-cocooned prop_ci
            
            for (int i = 0; i < A; i++) // population with.
                for (int j = 0; j < 12; j++) //     infants
                    nwn_e[i][j] = nwn_c(i, j)*(1-prop_c); //numbe

        }
        else
            cout << "Error cont" << endl;

        return nwn_e;
    }
    
    vector<vector<double > > get_nwp(double prop_c, char s)
    {
        // Population with non-cocooned parents)
        vector<vector<double > > nwp_e(A,vector<double>(A,0));
        if (s == 'p')
        {
            for (int i = 0; i < 12; i++) // infants
                for (int j = 18; j < 21; j++) // persons of cba
                    nwp_e[i][j] = (cnt_matrix_p_h(i, j)/2.0) + (cnt_matrix_p(i, j) - cnt_matrix_p_h(i, j))*p_mat[j]*(1-prop_c);
            
            for (int i = 13; i < 25; i++) // infants
                for (int j = 18; j < 21; j++) // persons of cba
                    nwp_e[i][j] = cnt_matrix_p(i, j)*p_mat[j]*(1-prop_c);
        }
        else if (s == 'c')
        {
            for (int i = 0; i < 12; i++) // infants
                for (int j = 18; j < 21; j++) // persons of cba
                    nwp_e[i][j] = (cnt_matrix_c_h(i, j)/2.0) + (cnt_matrix_c(i, j) - cnt_matrix_c_h(i, j))*p_mat[j]*(1-prop_c);
            
            for (int i = 13; i < 25; i++) // infants
                for (int j = 18; j < 21; j++) // persons of cba
                    nwp_e[i][j] = cnt_matrix_c(i, j)*p_mat[j]*(1-prop_c);
        }
        else
            cout << "Error cont" << endl;
        
        return nwp_e;
    }
    
    vector<vector<double > > get_nwc(double prop_c, char s)
    {
        
        vector<vector<double > > nwc_e(A,vector<double>(A,0));
        if (s == 'p')
        {
            for (int i = 0; i < A; i++) // population with.
                for (int j = 0; j < 12; j++) //     infants
                    nwc_e[i][j] = nwn_p(i, j)*(prop_c); //numbe
            
            for (int i = 0; i < 12; i++) // infants
                for (int j = 18; j < 21; j++) // persons of cba
                    nwc_e[i][j] = (cnt_matrix_p(i, j) - cnt_matrix_p_h(i, j))*p_mat[j]*(prop_c);
            
            for (int i = 13; i < 25; i++) // infants
                for (int j = 18; j < 21; j++) // persons of cba
                    nwc_e[i][j] = cnt_matrix_p(i, j)*p_mat[j]*(prop_c);
        }
        else if (s == 'c')
        {
            for (int i = 0; i < A; i++) // population with.
                for (int j = 0; j < 12; j++) //     infants
                    nwc_e[i][j] = nwn_c(i, j)*(prop_c); //numbe
            
            for (int i = 0; i < 12; i++) // infants
                for (int j = 18; j < 21; j++) // persons of cba
                    nwc_e[i][j] = (cnt_matrix_c(i, j) - cnt_matrix_c_h(i, j))*p_mat[j]*(prop_c);
            
            for (int i = 13; i < 25; i++) // infants
                for (int j = 18; j < 21; j++) // persons of cba
                    nwc_e[i][j] = cnt_matrix_c(i, j)*p_mat[j]*(prop_c);
        }
        else
            cout << "Error cont" << endl;
        
        return nwc_e;
    }
    
    vector< double >  populationPerAgeGroup, eta, modelIncidencePerTime, pA, ep_t;
    double dailyBirthRate, totPopulation;
    NumericVector ageStratification;
    
    int dayNoAfterBurn, weekNo, monthNo;
    double valueLogLikelihood;
    
    // Defined later but of size A
    double currentODETime;
    double run_start, run_burn, run_full, dt;

    RunInterventions(double dailyBirthRate_t, double totPopulation_t, NumericVector ageStratification_t): dailyBirthRate(dailyBirthRate_t), totPopulation(totPopulation_t), ageStratification(ageStratification_t){
           
            A = ageStratification.size();
            eta.push_back(0);
            for (int i = 0; i < A-1; i++){
                populationPerAgeGroup.push_back(dailyBirthRate*365*(ageStratification[i+1]-ageStratification[i]));
                eta.push_back( 1.0/(365.0*(ageStratification[i+1] - ageStratification[i])) );
                modelIncidencePerTime.push_back(0);
            }
            modelIncidencePerTime.push_back(0);
            populationPerAgeGroup.push_back(totPopulation - (dailyBirthRate*365)*ageStratification[A-1]);
            eta.push_back(dailyBirthRate/(totPopulation - (dailyBirthRate*365)*ageStratification[A-1]) );

            dt = 1;
            currentODETime = 0;
            dayNoAfterBurn = 0;
            valueLogLikelihood = 0;
            
            rg = 45;
            sg = rg * 3;
            ag = sg * 6 + 23;
    }
    
    vector<double > pVHR;
    vector<double > pHR;
    vector<double > pLR;
    
    // Related to parameter values
    NumericVector parameterValues;
    
    void ParameterValuesforODE(VectorXd currentParamValues){

        this->ep_t.clear();
        this->pA.clear();
        
        NumericVector parameterValuesTemp(25);
        
         parameterValuesTemp["xi"] = currentParamValues(0);
         parameterValuesTemp["si"] = currentParamValues(1);
         parameterValuesTemp["ga0"] = currentParamValues(2);
         parameterValuesTemp["g1"] = currentParamValues(3);
         parameterValuesTemp["g2"] = currentParamValues(4);
         parameterValuesTemp["om"] = currentParamValues(5);
         parameterValuesTemp["pA1"] = currentParamValues(6);
         parameterValuesTemp["pA2"] = currentParamValues(7);
         parameterValuesTemp["pA3"] = currentParamValues(8);
         parameterValuesTemp["pA4"] = currentParamValues(9);
         parameterValuesTemp["alpha_i"] =currentParamValues(10);
         parameterValuesTemp["d1"] = currentParamValues(11);
         parameterValuesTemp["d2"] = currentParamValues(12);
         parameterValuesTemp["d3"] = currentParamValues(13);
         parameterValuesTemp["phi"] = currentParamValues(14);
         parameterValuesTemp["qp"] = currentParamValues(15);
         parameterValuesTemp["qc"] = currentParamValues(16);
         parameterValuesTemp["b1"] = currentParamValues(17);
         parameterValuesTemp["psi"] = currentParamValues(18);
         parameterValuesTemp["c5ep1"] = currentParamValues(19);
         parameterValuesTemp["c5ep2"] = currentParamValues(20);
         parameterValuesTemp["ep5"] = currentParamValues(21);
         parameterValuesTemp["ep6"] = currentParamValues(22);
         parameterValuesTemp["I1"] = currentParamValues(23);
         parameterValuesTemp["I2"] = currentParamValues(24);
        this->parameterValues = parameterValuesTemp;
        
        // Define ep_t;
        for (int a = 0; a < 16; a++)
            this->ep_t.push_back(exp(parameterValues["c5ep1"] + a*parameterValues["c5ep2"]));
        for (int a = 16; a < 23; a++)
            this->ep_t.push_back(parameterValues["ep5"]);
        for (int a = 23; a < 25; a++)
            this->ep_t.push_back(parameterValues["ep6"]);
        // define pA;
        for (int a = 0; a < 12; a++)
            this->pA.push_back(parameterValues["pA1"]);
        for (int a = 12; a < 16; a++)
            this->pA.push_back(parameterValues["pA2"]);
        for (int a = 16; a < 18; a++)
            this->pA.push_back(parameterValues["pA3"]);
        for (int a = 18; a < 25; a++)
            this->pA.push_back(parameterValues["pA4"]);
    }
    
    vector<double > initial_M()
    {
        double xi = 1.0/this->parameterValues["xi"];
        boost::math::exponential_distribution <> exp ( xi );
        
        vector<double > init_con;
        for (int i = 0; i< this->A - 1; i++)
        {
            double init_con_temp = (cdf(exp, 365*ageStratification[i+1])-cdf(exp, 365*ageStratification[i]))/((365*ageStratification[i+1]-365*ageStratification[i])*xi);
            init_con.push_back(init_con_temp*populationPerAgeGroup[i]);
        }
        init_con.push_back( (cdf(exp, 365*90)-cdf(exp, 365*ageStratification[A-1]))/((365*90-365*ageStratification[A-1])*xi)*populationPerAgeGroup[A-1]);
        return init_con;
    }

    double poisson_cdf(double l, double a, double x){
      if( l == 0.0 || a == 0.0){
        boost::math::poisson_distribution<> p(0.000001); return cdf(p,x);
      }
      else{
        boost::math::poisson_distribution<> p(l*a); return cdf(p,x);
      }
    }

    vector< double >  initialProportionExposure(double l, double a1, double a2){
        vector< double >  prop(this->A);
        prop[0] = abs(poisson_cdf(l,a2,0)-poisson_cdf(l,a1,0))/((a2-a1)*l);
        prop[1] = abs(poisson_cdf(l,a2,1)-poisson_cdf(l,a1,1))/((a2-a1)*l);
        prop[2] = abs(poisson_cdf(l,a2,2)-poisson_cdf(l,a1,2))/((a2-a1)*l);
        prop[3] = 1 - (prop[2]+prop[1]+prop[0]);
        return prop;
    }

    vector<double > generateInitialStates(double cov_c)
    {
        vector<double > initialStates;
        vector<double > populationMatPro = initial_M();
        
        double s_prop, r_prop;
        double pI1; double pI2; double pI3; double pI4;
        
        double I1 = this->parameterValues["I1"];
        double I2 = this->parameterValues["I2"];
        double I3 = 0.5;
        
        double si = 1.0/this->parameterValues["si"];
        double g0 = 1.0/this->parameterValues["ga0"];
        double g1 = 1.0/(this->parameterValues["ga0"]*this->parameterValues["g1"]);
        double g2 = 1.0/(this->parameterValues["ga0"]*this->parameterValues["g1"]*this->parameterValues["g2"]);
        double d1 = this->parameterValues["d1"];
        double d2 = this->parameterValues["d1"]*this->parameterValues["d2"];
        double d3 = this->parameterValues["d1"]*this->parameterValues["d2"]*this->parameterValues["d3"];

        double a1, a2;
        
        for(int a = 0; a < A; a++)
        {
            if (a < A-1){
                a1 = this->ageStratification(a); a2 =  this->ageStratification(a+1);
            }
            else{
                a1 =  this->ageStratification(a); a2 = 90;
            }
            
            vector<double > propEachExposureGroup = initialProportionExposure(I3, a1, a2);
            pI1 = propEachExposureGroup[0]; pI2 = propEachExposureGroup[1]; pI3 = propEachExposureGroup[2]; pI4 = propEachExposureGroup[3];
            /** Get ratio of number of number of each person in each compoartment group **/
            double age_size = populationPerAgeGroup[a] - populationMatPro[a];
            
            for (int s = 0; s < 6; s++)
            {
                if (a < 12)
                {
                    if (s == 0 || s == 3)
                        s_prop = (0.0);
                    else if (s == 1 || s == 4)
                        s_prop = (cov_c);
                    else if (s == 2 || s == 5)
                        s_prop = (1-cov_c);
                    else
                        cout << "OOPS" << '\n';
                }
                else
                {
                    if (s == 0 || s == 3)
                        s_prop = p_mat[a]*(1-cov_c);
                    else if (s == 1 || s == 4)
                        s_prop = p_mat[a]*(cov_c);
                    else if (s == 2 || s == 5)
                        s_prop = 1-p_mat[a];
                    
                    else
                        cout << "OOPS" << '\n';
                }
                for (int r = 0; r < 3; r++)
                {
                    if (r == 0)
                        r_prop = pVHR[a];
                    else if (r == 1)
                        r_prop = pHR[a];
                    else if (r == 2)
                        r_prop = pLR[a];
                    else
                        cout << "ERROR" << endl;
                    
                    initialStates.push_back(r_prop*s_prop*populationMatPro[a]); // Number  in M group
                    
                    initialStates.push_back(r_prop*s_prop*pI1*age_size*(1.0 - I1)*(1.0-I2));      //Sus
                    initialStates.push_back(r_prop*s_prop*pI1*age_size*I1*si/(si+g0));       //Exp
                    initialStates.push_back(r_prop*s_prop*pI1*age_size*I1*g0/(si+g0)*pA[a]);
                    initialStates.push_back(r_prop*s_prop*pI1*age_size*I1*g0/(si+g0)*(1-pA[a]));   //Inf S
                    initialStates.push_back(r_prop*s_prop*pI1*age_size*(1.0 - I1)*I2);      //Rec
                    
                    initialStates.push_back(r_prop*s_prop*pI2*age_size*(1.0 - d1*I1)*(1.0-I2));      //Sus
                    initialStates.push_back(r_prop*s_prop*pI2*age_size*d1*I1*si/(si+g1));       //Exp
                    initialStates.push_back(r_prop*s_prop*pI2*age_size*d1*I1*g1/(si+g1)*pA[a]);
                    initialStates.push_back(r_prop*s_prop*pI2*age_size*d1*I1*g1/(si+g1)*(1-pA[a]));   //Inf S
                    initialStates.push_back(r_prop*s_prop*pI2*age_size*(1.0 - d1*I1)*I2);      //Rec
                    
                    initialStates.push_back(r_prop*s_prop*pI3*age_size*(1.0 - d2*I1)*(1.0-I2));      //S/Sus
                    initialStates.push_back(r_prop*s_prop*pI3*age_size*d2*I1*si/(si+g2));       //Exp
                    initialStates.push_back(r_prop*s_prop*pI3*age_size*d2*I1*g2/(si+g2)*pA[a]);
                    initialStates.push_back(r_prop*s_prop*pI3*age_size*d2*I1*g2/(si+g2)*(1-pA[a]));   //Inf S
                    initialStates.push_back(r_prop*s_prop*pI3*age_size*(1.0 - d2*I1)*I2);      //ec
                    
                    initialStates.push_back(r_prop*s_prop*pI4*age_size*(1.0 - d3*I1)*(1.0-I2));        //Sus
                    initialStates.push_back(r_prop*s_prop*pI4*age_size*d3*I1*si/(si+g2));       //Exp
                    initialStates.push_back(r_prop*s_prop*pI4*age_size*d3*I1*g2/(si+g2)*pA[a]);
                    initialStates.push_back(r_prop*s_prop*pI4*age_size*d3*I1*g2/(si+g2)*(1-pA[a]));   //Inf S
                    initialStates.push_back(r_prop*s_prop*pI4*age_size*(1.0 - d3*I1)*I2);     //Rec
                    
                    for (int j = 0; j < (this->rg - 21); j++) {
                        initialStates.push_back(0);     //Rec
                    }
                }
            }
            for (int j = 0; j < 23; j++)
                initialStates.push_back(0.0);      //inf exposure 1
        }
        return initialStates;
    }
    
    void getWeeklyIncidence(vector<double> &x0, NumericMatrix &sampleWeeklyIncidence, NumericMatrix &no_doses, bool epFlag)
    {
        if (this->dayNoAfterBurn % 365 == 0) {
            this->dayNoAfterBurn = 0;
            for (int a = 0; a < this->A; a++) {
                x0[this->ag*a + this->sg*6 + 1] = 0;
                x0[this->ag*a + this->sg*6 + 2] = 0;
                x0[this->ag*a + this->sg*6 + 3] = 0;
                x0[this->ag*a + this->sg*6 + 4] = 0;
                for (int j = 0; j < 9; j++) {
                    x0[this->ag*a + this->sg*6 + 8 + j] = 0.0; //Incidence at t_d = 0;
                }
            }
        }
      //if (this->dayNoAfterBurn == 0){
      //  for (int a = 0; a < this->A; a++) {
       //     x0[this->ag*a + this->sg*6 + 1] = 0;
        //    x0[this->ag*a + this->sg*6 + 2] = 0;
         //   x0[this->ag*a + this->sg*6 + 3] = 0;
          //  x0[this->ag*a + this->sg*6 + 4] = 0;
       //     for (int j = 0; j < 9; j++) {
       //         x0[this->ag*a + this->sg*6 + 8 + j] = 0.0; //Incidence at t_d = 0;
       //     }
      //  }
     // }
        if (this->dayNoAfterBurn%7 == 0 && this->dayNoAfterBurn > 0)
        {
            for (int a = 0; a < this->A; a++) {
                for (int j = 0; j < 9; j++)
                {
                    if (epFlag)
                      sampleWeeklyIncidence(this->weekNo, 9*a + j) = x0[this->ag*a + this->sg*6 + 8 + j]*this->ep_t[a];
                    else
                      sampleWeeklyIncidence(this->weekNo, 9*a + j) = x0[this->ag*a + this->sg*6 + 8 + j]; //Incidence at t_d = 7;
                  
                    x0[this->ag*a + this->sg*6 + 8 + j] = 0.0;
                }
            }
            for (int a = 0; a < this->A; a++) {
                no_doses(this->weekNo, 0) += (x0[this->ag*a + this->sg*6 + 1]); // pal
                no_doses(this->weekNo, 1) += (x0[this->ag*a + this->sg*6 + 2]); // mab
                no_doses(this->weekNo, 2) += (x0[this->ag*a + this->sg*6 + 3]); // lav
                no_doses(this->weekNo, 3) += (x0[this->ag*a + this->sg*6 + 4]); // mat
                
                x0[this->ag*a + this->sg*6 + 1] = 0;
                x0[this->ag*a + this->sg*6 + 2] = 0;
                x0[this->ag*a + this->sg*6 + 3] = 0;
                x0[this->ag*a + this->sg*6 + 4] = 0;
            }
            this->weekNo++;
        }

        this->dayNoAfterBurn++;
    }
};


class ODE_desc
{
    RunInterventions* finELL;
  
public:
    double M, S0, S1, S2, S3, SN, E0, E1, E2, E3, A0, A1, A2, A3, I0, I1, I2, I3, R0, R1, R2, R3, M_, S0_, S1_, S2_, S3_, E0_, E1_, E2_, E3_, A0_, A1_, A2_, A3_, I0_, I1_, I2_, I3_, R0_, R1_, R2_, R3_, N;
    double Mv, S0v, S1v, S2v, S3v, E0v, E1v, E2v, E3v, A0v, A1v, A2v, A3v, I0v, I1v, I2v, I3v, R0v, R1v, R2v, R3v, Mv_, S0v_, S1v_, S2v_, S3v_, E0v_, E1v_, E2v_, E3v_, A0v_, A1v_, A2v_, A3v_, I0v_, I1v_, I2v_, I3v_, R0v_, R1v_, R2v_, R3v_;
    double xi, si, ga0, ga1, ga2, ga3, d1, d2, d3, a1, a2, a3, alpha_i, rho, om, b1, qp, qc, psi, phi, beta;
    int A;
    double p_vul;
    VectorXd I_n;     VectorXd I_c;     VectorXd I_p;
    VectorXd I_n_v;   VectorXd I_c_v;   VectorXd I_p_v;
    
    VectorXd N_tot;
    VectorXd N_tot_n;         VectorXd N_tot_c;         VectorXd N_tot_p;
    VectorXd N_tot_n_v;       VectorXd N_tot_c_v;       VectorXd N_tot_p_v;
    double N_tot_n_v_inv, N_tot_p_v_inv, N_tot_c_v_inv, N_tot_n_inv, N_tot_p_inv, N_tot_c_inv;
    
    VectorXd prop_n;  VectorXd prop_c;  VectorXd prop_p;
    VectorXd prop_nv; VectorXd prop_cv; VectorXd prop_pv;
    VectorXd prop_empty;
    
    VectorXd PS;
    
    double I_temp_p, I_temp_c, I_temp_n, I_temp_p_v, I_temp_c_v, I_temp_n_v;
    double x_tot, x_tot_1, x_tot_2, x_tot_3;
    double xi_b;
    
    vector<double > pVHR;
    vector<double > pHR;
    vector<double > pLR;
    vector<double > p_mat;

    double dailyBirthRate;
    NumericMatrix contactMatrixPhy, contactMatrixCon;
    vector< double >  populationPerAgeGroup, eta, pA;
    vector<double > VHR_g = {0.653012, 0.291842, 1., 1., 1., 0.0321774, 1., 1., 0,0,0,0,0,0,0,0,0,0,0};
    vector<double > LR_g = {1.00126, 1.00168, 1., 1., 1., 1.00067, 1., 1., 1.00002, 1, 1, 1, 1, 1, 1, 1, 1};
    
    MatrixXd cal_pal, cal_mAB_VHR, cal_mAB_HR, cal_mAB_LR, cal_LAV_HR, cal_LAV_LR, cal_mat_LR;
    MatrixXd cal_mAB_VHR_step_out, cal_mAB_VHR_step_in, cal_mAB_HR_step_out, cal_mAB_HR_step_in, cal_mAB_LR_step_out, cal_mAB_LR_step_in;


    MatrixXd cal_pal_t, cal_mAB_VHR_t, cal_mAB_HR_t, cal_mAB_LR_t, cal_LAV_HR_t, cal_LAV_LR_t, cal_mat_LR_t;
    MatrixXd cal_mAB_VHR_step_out_t, cal_mAB_VHR_step_in_t, cal_mAB_HR_step_out_t, cal_mAB_HR_step_in_t, cal_mAB_LR_step_out_t, cal_mAB_LR_step_in_t;

    MatrixXd cal_pal_dose, cal_mAB_VHR_dose, cal_mAB_HR_dose, cal_mAB_LR_dose, cal_LAV_HR_dose, cal_LAV_LR_dose, vac_cal_dose;
    MatrixXd cal_mAB_VHR_step_out_dose, cal_mAB_VHR_step_in_dose, cal_mAB_HR_step_out_dose, cal_mAB_HR_step_in_dose, cal_mAB_LR_step_out_dose, cal_mAB_LR_step_in_dose;

    List vac_calendar, vac_info, vac_info_vhr, vac_dose;
    double cov_c, eff_wane, eff_wane_lav, eff_wane_mat, eff_wane_vhr, eff_wane_mab, eff_wane_step, xi_boost;
    bool direct;

    double u18p;
    double u19p;
    double u20p;
    int rg, sg, ag;
    
  ODE_desc(RunInterventions* finELL_t, List vac_calendar_t, List vac_dose_t, List vac_info_t, double cov_c_t): finELL(finELL_t), vac_calendar(vac_calendar_t), vac_dose(vac_dose_t), vac_info(vac_info_t),  cov_c(cov_c_t){

        xi = 1.0/finELL->parameterValues["xi"];
        si = 1.0/finELL->parameterValues["si"];
        ga0 = 1.0/(finELL->parameterValues["ga0"]);
        ga1 = 1.0/(finELL->parameterValues["ga0"]*finELL->parameterValues["g1"]);
        ga2 = 1.0/(finELL->parameterValues["ga0"]*finELL->parameterValues["g1"]*finELL->parameterValues["g2"]);
        ga3 = ga2;
        om = 1.0/finELL->parameterValues["om"];

        rho = 1.0;
        alpha_i = finELL->parameterValues["alpha_i"];
        d1 = finELL->parameterValues["d1"];
        d2 = finELL->parameterValues["d1"]*finELL->parameterValues["d2"];
        d3 = finELL->parameterValues["d1"]*finELL->parameterValues["d2"]*finELL->parameterValues["d3"];
        a1 = 1.0, a2 = 1.0, a3 = 1.0;

        phi = finELL->parameterValues["phi"];
        qp = finELL->parameterValues["qp"];
        qc = finELL->parameterValues["qc"];
        b1 = finELL->parameterValues["b1"];
        psi = finELL->parameterValues["psi"];

        A = finELL->A;
        eta = finELL->eta;
      
        u18p = finELL->u_p[0];
        u19p = finELL->u_p[1];
        u20p = finELL->u_p[2];

        populationPerAgeGroup = finELL->populationPerAgeGroup;
        pVHR = finELL->pVHR;
        pHR = finELL->pHR;
        pLR = finELL->pLR;
        p_mat = finELL->p_mat;

        dailyBirthRate = finELL->dailyBirthRate;
        pA = finELL->pA;
      
        I_n = VectorXd::Zero(A);
        I_c = VectorXd::Zero(A);
        I_p = VectorXd::Zero(A);
        I_n_v = VectorXd::Zero(A);
        I_c_v = VectorXd::Zero(A);
        I_p_v = VectorXd::Zero(A);
      
        N_tot = VectorXd::Zero(A);
        N_tot_n = VectorXd::Zero(A);
        N_tot_c = VectorXd::Zero(A);
        N_tot_p = VectorXd::Zero(A);
        N_tot_n_v = VectorXd::Zero(A);
        N_tot_c_v = VectorXd::Zero(A);
        N_tot_p_v = VectorXd::Zero(A);
      
        prop_n = VectorXd::Zero(A);   prop_c = VectorXd::Zero(A);   prop_p = VectorXd::Zero(A);
        prop_nv = VectorXd::Zero(A);  prop_cv = VectorXd::Zero(A);  prop_pv = VectorXd::Zero(A);
        prop_empty = VectorXd::Zero(A);
      
        PS = VectorXd::Zero(finELL->rg);

        // This needs sorting out somehow
         cal_mAB_VHR_t = vac_calendar["mAB_VHR"];
         cal_mAB_LR_t = vac_calendar["mAB_LR"];
         cal_LAV_LR_t = vac_calendar["LAV_LR"];
         cal_mat_LR_t = vac_calendar["mat_LR"];
      
         cal_mAB_VHR_step_out_t = vac_calendar["mAB_VHR_custom_out"];
         cal_mAB_LR_step_out_t = vac_calendar["mAB_LR_custom_out"];

        cal_mAB_VHR_dose = vac_dose["mAB_VHR"];
        cal_mAB_LR_dose = vac_dose["mAB_LR"];
        cal_LAV_LR_dose = vac_dose["LAV_LR"];
        vac_cal_dose = vac_dose["mat_LR"];
      
        cal_mAB_VHR_step_out_dose = vac_dose["mAB_VHR_custom_out"];
        cal_mAB_LR_step_out_dose = vac_dose["mAB_LR_custom_out"];

        cal_mAB_HR_step_out_dose = MatrixXd::Zero(365, A);
        cal_mAB_HR_dose = MatrixXd::Zero(365, A);
        cal_LAV_HR_dose = MatrixXd::Zero(365, A);

        direct = vac_info["direct"];

        // Sort out efficacy stuff here too
        eff_wane_lav = vac_info["lav_mass"];
        eff_wane_mat = vac_info["mat_mass"];
        eff_wane_mab = vac_info["mab_mass"];
        eff_wane_vhr = vac_info["mab_vhr"];
        eff_wane_step = 0;

        xi_boost = 1;

        rg = finELL->rg;
        sg = finELL->sg;
        ag = finELL->ag;
   };
  
    
  
  void operator() (  vector< double >  &x , vector< double >  &dxdt , const double  t, MatrixXd  &collect_protect )
    {
      /*/////////////*/
      /* CALENDARS /*/
      /*/////////////*/
      MatrixXd vac_cal_vhr; MatrixXd vac_cal;
      if (t < 2 * 365)
      {
          cal_pal = MatrixXd::Zero(365, A);
          vac_cal = MatrixXd::Zero(365, A);
          cal_mAB_VHR = MatrixXd::Zero(365, A);
          cal_mAB_HR = MatrixXd::Zero(365, A);
          cal_mAB_LR = MatrixXd::Zero(365, A);
          cal_LAV_HR = MatrixXd::Zero(365, A);
          cal_LAV_LR = MatrixXd::Zero(365, A);

        /* THE STEP FUNCTION IMMUNITY */
        cal_mAB_VHR_step_out = MatrixXd::Zero(365, A);
        cal_mAB_HR_step_out = MatrixXd::Zero(365, A);
        cal_mAB_LR_step_out = MatrixXd::Zero(365, A);


      }
      else{
            cal_mAB_VHR = cal_mAB_VHR_t;
            cal_mAB_LR = cal_mAB_LR_t;
            cal_LAV_LR = cal_LAV_LR_t;          
            vac_cal = cal_mat_LR_t;


       /* THE STEP FUNCTION IMMUNITY */
        cal_mAB_VHR_step_out = cal_mAB_VHR_step_out_t;
        cal_mAB_LR_step_out = cal_mAB_LR_step_out_t;
      }

      /*/////////////*/
      /* POPULATIONS /*/
      /*/////////////*/
      //Find sum of every group
      //   Rcpp::Rcout << "Info on popoulations" << std::endl;

      double tot_coc = 0, tot_coc_v = 0, tot_inf = 0, tot_inf_v = 0;
      
      for(int a = 0; a < A; a++)
      {
      /*    double r_prop = 0;
          double tot_temp_n = 0, tot_temp_c = 0, tot_temp_p = 0, tot_temp_nv = 0, tot_temp_cv = 0, tot_temp_pv = 0;
          double tot_temp_vhr = 0, tot_temp_vhrv = 0, tot_temp_hr = 0, tot_temp_hrv = 0, tot_temp_lr = 0, tot_temp_lrv = 0;

          for (int j = 0; j < 24; j++){
              for (int r = 0; r < 3; r++){
                  tot_temp_p +=  x[a*ag+0*sg+r*rg+j];        tot_temp_c += x[a*ag+1*sg+r*rg+j];            tot_temp_n += x[a*ag+2*sg+r*rg+j];
                  tot_temp_pv += x[a*ag+3*sg+r*rg+j];       tot_temp_cv += x[a*ag+4*sg+r*rg+j];           tot_temp_nv += x[a*ag+5*sg+r*rg+j];
              }
          }
          
          for (int j = 0; j < 24; j++){
              for (int s = 0; s < 3; s++){
                  tot_temp_vhr +=  x[a*ag+s*sg+0*rg+j];        tot_temp_hr += x[a*ag+s*sg+1*rg+j];            tot_temp_lr += x[a*ag+s*sg+2*rg+j];
                  tot_temp_vhrv += x[a*ag+(s+3)*sg+0*rg+j];       tot_temp_hrv += x[a*ag+(s+3)*sg+1*rg+j];           tot_temp_lrv += x[a*ag+(s+3)*sg+2*rg+j];
              }
          }
          
          if (a < 12){
              tot_coc += tot_temp_c;                  tot_coc_v += tot_temp_cv;
              tot_inf += tot_temp_c + tot_temp_n;     tot_inf_v += tot_temp_cv + tot_temp_nv;
          }
        // for (int t = 0; t < 365; t++) {
        //  tot_temp_n = square_cal_in[t, a] - square_cal_out[t, a]
        //  }
        //*/
       //   N_tot_n[a] = tot_temp_n;    N_tot_c[a] = 0;    N_tot_p[a] = tot_temp_p;
       //   N_tot_n_v[a] = tot_temp_nv; N_tot_c_v[a] = 0; N_tot_p_v[a] = tot_temp_pv;
            if (a < 12) {
                N_tot_n[a] = populationPerAgeGroup[a] *  (1 - cov_c);    N_tot_c[a] = populationPerAgeGroup[a] * cov_c; N_tot_p[a] = 0;
                N_tot_n_v[a] = populationPerAgeGroup[a] *  (1 - cov_c); N_tot_c_v[a] = populationPerAgeGroup[a] * cov_c; N_tot_p_v[a] = 0;
            } else {
                N_tot_n[a] = populationPerAgeGroup[a] * (1 - p_mat[a]);    N_tot_c[a] = populationPerAgeGroup[a] * p_mat[a] * cov_c;    N_tot_p[a] = populationPerAgeGroup[a] * p_mat[a] * (1 - cov_c);
                N_tot_n_v[a] = populationPerAgeGroup[a] *  (1 - p_mat[a]); N_tot_c_v[a] = populationPerAgeGroup[a] * p_mat[a] * cov_c; N_tot_p_v[a] = populationPerAgeGroup[a] * p_mat[a] * (1 - cov_c);               
            }
      }
      
      /*/////////////*/
      /* 2. Force of Infection /*/
      /*/////////////*/
     //  Rcpp::Rcout << "Info on FOI" << std::endl;
 
      double phi_c = cov_c;
      double phi_c_v = cov_c;
      
      //  2.1 Get contact matrices
      vector<vector< double > > cnt_matrix_cwc_p = finELL->get_cwc(phi_c,'p');
      vector<vector< double > > cnt_matrix_cwp_p = finELL->get_cwp(phi_c,'p');
      vector<vector< double > > cnt_matrix_cwn_p = finELL->get_cwn(phi_c,'p');
      vector<vector< double > > cnt_matrix_pwc_p = finELL->get_pwc(phi_c,'p');
      vector<vector< double > > cnt_matrix_pwp_p = finELL->get_pwp(phi_c,'p');
      vector<vector< double > > cnt_matrix_pwn_p = finELL->get_pwn(phi_c,'p');
      vector<vector< double > > cnt_matrix_nwc_p = finELL->get_nwc(phi_c,'p');
      vector<vector< double > > cnt_matrix_nwp_p = finELL->get_nwp(phi_c,'p');
      vector<vector< double > > cnt_matrix_nwn_p = finELL->get_nwn(phi_c,'p');
      
      vector<vector< double > > cnt_matrix_cwc_c = finELL->get_cwc(phi_c,'c');
      vector<vector< double > > cnt_matrix_cwp_c = finELL->get_cwp(phi_c,'c');
      vector<vector< double > > cnt_matrix_cwn_c = finELL->get_cwn(phi_c,'c');
      vector<vector< double > > cnt_matrix_pwc_c = finELL->get_pwc(phi_c,'c');
      vector<vector< double > > cnt_matrix_pwp_c = finELL->get_pwp(phi_c,'c');
      vector<vector< double > > cnt_matrix_pwn_c = finELL->get_pwn(phi_c,'c');
      vector<vector< double > > cnt_matrix_nwc_c = finELL->get_nwc(phi_c,'c');
      vector<vector< double > > cnt_matrix_nwp_c = finELL->get_nwp(phi_c,'c');
      vector<vector< double > > cnt_matrix_nwn_c = finELL->get_nwn(phi_c,'c');

      vector<vector< double > > cnt_matrix_cwc_p_v = finELL->get_cwc(phi_c_v,'p');
      vector<vector< double > > cnt_matrix_cwp_p_v = finELL->get_cwp(phi_c_v,'p');
      vector<vector< double > > cnt_matrix_cwn_p_v = finELL->get_cwn(phi_c_v,'p');
      vector<vector< double > > cnt_matrix_pwc_p_v = finELL->get_pwc(phi_c_v,'p');
      vector<vector< double > > cnt_matrix_pwp_p_v = finELL->get_pwp(phi_c_v,'p');
      vector<vector< double > > cnt_matrix_pwn_p_v = finELL->get_pwn(phi_c_v,'p');
      vector<vector< double > > cnt_matrix_nwc_p_v = finELL->get_nwc(phi_c_v,'p');
      vector<vector< double > > cnt_matrix_nwp_p_v = finELL->get_nwp(phi_c_v,'p');
      vector<vector< double > > cnt_matrix_nwn_p_v = finELL->get_nwn(phi_c_v,'p');
      
      vector<vector< double > > cnt_matrix_cwc_c_v = finELL->get_cwc(phi_c_v,'c');
      vector<vector< double > > cnt_matrix_cwp_c_v = finELL->get_cwp(phi_c_v,'c');
      vector<vector< double > > cnt_matrix_cwn_c_v = finELL->get_cwn(phi_c_v,'c');
      vector<vector< double > > cnt_matrix_pwc_c_v = finELL->get_pwc(phi_c_v,'c');
      vector<vector< double > > cnt_matrix_pwp_c_v = finELL->get_pwp(phi_c_v,'c');
      vector<vector< double > > cnt_matrix_pwn_c_v = finELL->get_pwn(phi_c_v,'c');
      vector<vector< double > > cnt_matrix_nwc_c_v = finELL->get_nwc(phi_c_v,'c');
      vector<vector< double > > cnt_matrix_nwp_c_v = finELL->get_nwp(phi_c_v,'c');
      vector<vector< double > > cnt_matrix_nwn_c_v = finELL->get_nwn(phi_c_v,'c');
      
      //  2.2 Seasonal forcing
      //beta = (1 + b1*cos((t/365.0-phi)*2*PI));
      double t1 = (int)t%365;
      beta = (1 + b1*(1 + exp(-((t1/365.0 - phi))*((t1/365.0 - phi))/(2*psi*psi))));

      /*/////////////*/
      /*  3. DYNAMIC MATERNAL PROTECTION /*/
      /*/////////////*/

     // num_vec pos_wcb_ag = {18,21};
      double sum_wcb = 0;
      double sum_wcb_v = 0;
      double CB2_temp, CB2_temp_v;
      CB2_temp = CB2_temp_v = 0;
      for(int a = 18; a < 21 ; a++)
      {
          for (int r = 0; r < 3; r++)
          {
              CB2_temp += (x[a*ag+2*sg+r*rg+1] + x[a*ag+2*sg+r*rg+6] + x[a*ag+2*sg+r*rg+11] + x[a*ag+2*sg+r*rg+16] + x[a*ag+2*sg+r*rg+2] + x[a*ag+2*sg+r*rg+7] + x[a*ag+2*sg+r*rg+12] + x[a*ag+2*sg+r*rg+17]);
              CB2_temp_v += (x[a*ag+5*sg+r*rg+1] + x[a*ag+5*sg+r*rg+6] + x[a*ag+5*sg+r*rg+11] + x[a*ag+5*sg+r*rg+16] + x[a*ag+5*sg+r*rg+2] + x[a*ag+5*sg+r*rg+7] + x[a*ag+5*sg+r*rg+12] + x[a*ag+5*sg+r*rg+17]);
          }
      }
      sum_wcb = CB2_temp/((double)N_tot_n[18]+(double)N_tot_n[19]+(double)N_tot_n[20]);
      sum_wcb_v = (CB2_temp_v)/((double)N_tot_n[18]+(double)N_tot_n[19]+(double)N_tot_n[20]);
      
      double tot_C = N_tot_c[18] + N_tot_c[19] + N_tot_c[20];
      double tot_C_inv;
      if (tot_C < 1.0)
          tot_C_inv = 0;
      else
          tot_C_inv = 1.0/tot_C;
      /*/////////////*/
      /* 4. ODES /*/
      /*/////////////*/
      //4.1 Age groups
      double protectpal;
      double protectmabs;
      double protectLAV;
      double protectmat;
      double step_in, step_out, step_out_dose, step_in_S0, step_in_S1;
      double loss_mat, loss_lav_S0, loss_lav_S1, loss_lav_S2, loss_lav_S3;

     //  Rcpp::Rcout << "ODEs" << std::endl;

      for(int a = 0; a < A; a++)
      {
          double r_prop = 0;
          protectpal = protectmabs = protectLAV = protectmat = 0; step_out_dose = 0, step_in = 0,  step_out = 0; step_in_S0 = 0,  step_in_S1 = 0;
          I_temp_c = 0.0;    I_temp_p = 0.0; I_temp_n = 0.0;
          for (int k = 0; k < A; k++)
          {
              if (N_tot_n[k] < 0.1){N_tot_n_inv = 0;}else{N_tot_n_inv=1.0/N_tot_n[k];}
              if (N_tot_c[k] < 0.1){N_tot_c_inv = 0;}else{N_tot_c_inv=1.0/N_tot_c[k];}
              if (N_tot_p[k] < 0.1){N_tot_p_inv = 0;}else{N_tot_p_inv=1.0/N_tot_p[k];}
              
              for (int r = 0; r < 3; r++)
              {
                  I_temp_p += (x[k*ag+0*sg+r*rg+3]*alpha_i+x[k*ag+0*sg+r*rg+4]+a1*(x[k*ag+0*sg+r*rg+8]*alpha_i+x[k*ag+0*sg+r*rg+9])+a2*(x[k*ag+0*sg+r*rg+13]*alpha_i+x[k*ag+0*sg+r*rg+14])+a3*(x[k*ag+0*sg+r*rg+18]*alpha_i+x[k*ag+0*sg+r*rg+19]))*(qp*(cnt_matrix_pwp_p[a][k]+qc*cnt_matrix_pwp_c[a][k]))*N_tot_p_inv +(x[k*ag+1*sg+r*rg+3]*alpha_i+x[k*ag+1*sg+r*rg+4]+a1*(x[k*ag+1*sg+r*rg+8]*alpha_i+x[k*ag+1*sg+r*rg+9])+a2*(x[k*ag+1*sg+r*rg+13]*alpha_i+x[k*ag+1*sg+r*rg+14])+a3*(x[k*ag+1*sg+r*rg+18]*alpha_i+x[k*ag+1*sg+r*rg+19]))*(qp*(cnt_matrix_pwc_p[a][k]+qc*cnt_matrix_pwc_c[a][k]))*N_tot_c_inv +(x[k*ag+2*sg+r*rg+3]*alpha_i+x[k*ag+2*sg+r*rg+4]+a1*(x[k*ag+2*sg+r*rg+8]*alpha_i+x[k*ag+2*sg+r*rg+9])+a2*(x[k*ag+2*sg+r*rg+13]*alpha_i+x[k*ag+2*sg+r*rg+14])+a3*(x[k*ag+2*sg+r*rg+18]*alpha_i+x[k*ag+2*sg+r*rg+19]))*(qp*(cnt_matrix_pwn_p[a][k]+qc*cnt_matrix_pwn_c[a][k]))*N_tot_n_inv;
                  
                  I_temp_c += (x[k*ag+0*sg+r*rg+3]*alpha_i+x[k*ag+0*sg+r*rg+4]+a1*(x[k*ag+0*sg+r*rg+8]*alpha_i+x[k*ag+0*sg+r*rg+9])+a2*(x[k*ag+0*sg+r*rg+13]*alpha_i+x[k*ag+0*sg+r*rg+14])+a3*(x[k*ag+0*sg+r*rg+18]*alpha_i+x[k*ag+0*sg+r*rg+19]))*(qp*(cnt_matrix_cwp_p[a][k]+qc*cnt_matrix_cwp_c[a][k]))*N_tot_p_inv +(x[k*ag+1*sg+r*rg+3]*alpha_i+x[k*ag+1*sg+r*rg+4]+a1*(x[k*ag+1*sg+r*rg+8]*alpha_i+x[k*ag+1*sg+r*rg+9])+a2*(x[k*ag+1*sg+r*rg+13]*alpha_i+x[k*ag+1*sg+r*rg+14])+a3*(x[k*ag+1*sg+r*rg+18]*alpha_i+x[k*ag+1*sg+r*rg+19]))*(qp*(cnt_matrix_cwc_p[a][k]+qc*cnt_matrix_cwc_c[a][k]))*N_tot_c_inv +(x[k*ag+2*sg+r*rg+3]*alpha_i+x[k*ag+2*sg+r*rg+4]+a1*(x[k*ag+2*sg+r*rg+8]*alpha_i+x[k*ag+2*sg+r*rg+9])+a2*(x[k*ag+2*sg+r*rg+13]*alpha_i+x[k*ag+2*sg+r*rg+14])+a3*(x[k*ag+2*sg+r*rg+18]*alpha_i+x[k*ag+2*sg+r*rg+19]))*(qp*(cnt_matrix_cwn_p[a][k]+qc*cnt_matrix_cwn_c[a][k]))*N_tot_n_inv;
                  
                  I_temp_n += (x[k*ag+0*sg+r*rg+3]*alpha_i+x[k*ag+0*sg+r*rg+4]+a1*(x[k*ag+0*sg+r*rg+8]*alpha_i+x[k*ag+0*sg+r*rg+9])+a2*(x[k*ag+0*sg+r*rg+13]*alpha_i+x[k*ag+0*sg+r*rg+14])+a3*(x[k*ag+0*sg+r*rg+18]*alpha_i+x[k*ag+0*sg+r*rg+19]))*(qp*(cnt_matrix_nwp_p[a][k]+qc*cnt_matrix_nwp_c[a][k]))*N_tot_p_inv +(x[k*ag+1*sg+r*rg+3]*alpha_i+x[k*ag+1*sg+r*rg+4]+a1*(x[k*ag+1*sg+r*rg+8]*alpha_i+x[k*ag+1*sg+r*rg+9])+a2*(x[k*ag+1*sg+r*rg+13]*alpha_i+x[k*ag+1*sg+r*rg+14])+a3*(x[k*ag+1*sg+r*rg+18]*alpha_i+x[k*ag+1*sg+r*rg+19]))*(qp*(cnt_matrix_nwc_p[a][k]+qc*cnt_matrix_nwc_c[a][k]))*N_tot_c_inv +(x[k*ag+2*sg+r*rg+3]*alpha_i+x[k*ag+2*sg+r*rg+4]+a1*(x[k*ag+2*sg+r*rg+8]*alpha_i+x[k*ag+2*sg+r*rg+9])+a2*(x[k*ag+2*sg+r*rg+13]*alpha_i+x[k*ag+2*sg+r*rg+14])+a3*(x[k*ag+2*sg+r*rg+18]*alpha_i+x[k*ag+2*sg+r*rg+19]))*(qp*(cnt_matrix_nwn_p[a][k]+qc*cnt_matrix_nwn_c[a][k]))*N_tot_n_inv;
              }

          }

          
          I_temp_c_v = 0.0;  I_temp_p_v = 0.0; I_temp_n_v = 0.0;
          for (int k = 0; k < A ; k++)
          {
              if (N_tot_n_v[k] < 1){N_tot_n_v_inv = 0;}else{N_tot_n_v_inv=1.0/N_tot_n_v[k];}
              if (N_tot_c_v[k] < 1){N_tot_c_v_inv = 0;}else{N_tot_c_v_inv=1.0/N_tot_c_v[k];}
              if (N_tot_p_v[k] < 1){N_tot_p_v_inv = 0;}else{N_tot_p_v_inv=1.0/N_tot_p_v[k];}
              
              for (int r = 0; r < 3; r++)
              {
                  I_temp_p_v += (x[k*ag+3*sg+r*rg+3]*alpha_i+x[k*ag+3*sg+r*rg+4]+a1*(x[k*ag+3*sg+r*rg+8]*alpha_i+x[k*ag+3*sg+r*rg+9])+a2*(x[k*ag+3*sg+r*rg+13]*alpha_i+x[k*ag+3*sg+r*rg+14])+a3*(x[k*ag+3*sg+r*rg+18]*alpha_i+x[k*ag+3*sg+r*rg+19]))*(qp*(cnt_matrix_pwp_p_v[a][k]+qc*cnt_matrix_pwp_c_v[a][k]))*N_tot_p_v_inv + (x[k*ag+4*sg+r*rg+3]*alpha_i+x[k*ag+4*sg+r*rg+4]+a1*(x[k*ag+4*sg+r*rg+8]*alpha_i+x[k*ag+4*sg+r*rg+9])+a2*(x[k*ag+4*sg+r*rg+13]*alpha_i+x[k*ag+4*sg+r*rg+14])+a3*(x[k*ag+4*sg+r*rg+18]*alpha_i+x[k*ag+4*sg+r*rg+19]))*(qp*(cnt_matrix_pwc_p_v[a][k]+qc*cnt_matrix_pwc_c_v[a][k]))*N_tot_c_v_inv  + (x[k*ag+5*sg+r*rg+3]*alpha_i+x[k*ag+5*sg+r*rg+4]+a1*(x[k*ag+5*sg+r*rg+8]*alpha_i+x[k*ag+5*sg+r*rg+9])+a2*(x[k*ag+5*sg+r*rg+13]*alpha_i+x[k*ag+5*sg+r*rg+14])+a3*(x[k*ag+5*sg+r*rg+18]*alpha_i+x[k*ag+5*sg+r*rg+19]))*(qp*(cnt_matrix_pwn_p_v[a][k]+qc*cnt_matrix_pwn_c_v[a][k]))*N_tot_n_v_inv;
                  
                  I_temp_c_v += (x[k*ag+3*sg+r*rg+3]*alpha_i+x[k*ag+3*sg+r*rg+4]+a1*(x[k*ag+3*sg+r*rg+8]*alpha_i+x[k*ag+3*sg+r*rg+9])+a2*(x[k*ag+3*sg+r*rg+13]*alpha_i+x[k*ag+3*sg+r*rg+14])+a3*(x[k*ag+3*sg+r*rg+18]*alpha_i+x[k*ag+3*sg+r*rg+19]))*(qp*(cnt_matrix_cwp_p_v[a][k]+qc*cnt_matrix_cwp_c_v[a][k]))*N_tot_p_v_inv +(x[k*ag+4*sg+r*rg+3]*alpha_i+x[k*ag+4*sg+r*rg+4]+a1*(x[k*ag+4*sg+r*rg+8]*alpha_i+x[k*ag+4*sg+r*rg+9])+a2*(x[k*ag+4*sg+r*rg+13]*alpha_i+x[k*ag+4*sg+r*rg+14])+a3*(x[k*ag+4*sg+r*rg+18]*alpha_i+x[k*ag+4*sg+r*rg+19]))*(qp*(cnt_matrix_cwc_p_v[a][k]+qc*cnt_matrix_cwc_c_v[a][k]))*N_tot_c_v_inv +(x[k*ag+5*sg+r*rg+3]*alpha_i+x[k*ag+5*sg+r*rg+4]+a1*(x[k*ag+5*sg+r*rg+8]*alpha_i+x[k*ag+5*sg+r*rg+9])+a2*(x[k*ag+5*sg+r*rg+13]*alpha_i+x[k*ag+5*sg+r*rg+14])+a3*(x[k*ag+5*sg+r*rg+18]*alpha_i+x[k*ag+5*sg+r*rg+19]))*(qp*(cnt_matrix_cwn_p_v[a][k]+qc*cnt_matrix_cwn_c_v[a][k]))*N_tot_n_v_inv;
                  
                  I_temp_n_v += (x[k*ag+3*sg+r*rg+3]*alpha_i+x[k*ag+3*sg+r*rg+4]+a1*(x[k*ag+3*sg+r*rg+8]*alpha_i+x[k*ag+3*sg+r*rg+9])+a2*(x[k*ag+3*sg+r*rg+13]*alpha_i+x[k*ag+3*sg+r*rg+14])+a3*(x[k*ag+3*sg+r*rg+18]*alpha_i+x[k*ag+3*sg+r*rg+19]))*(qp*(cnt_matrix_nwp_p_v[a][k]+qc*cnt_matrix_nwp_c_v[a][k]))*N_tot_p_v_inv + (x[k*ag+4*sg+r*rg+3]*alpha_i+x[k*ag+4*sg+r*rg+4]+a1*(x[k*ag+4*sg+r*rg+8]*alpha_i+x[k*ag+4*sg+r*rg+9])+a2*(x[k*ag+4*sg+r*rg+13]*alpha_i+x[k*ag+4*sg+r*rg+14])+a3*(x[k*ag+4*sg+r*rg+18]*alpha_i+x[k*ag+4*sg+r*rg+19]))*(qp*(cnt_matrix_nwc_p_v[a][k]+qc*cnt_matrix_nwc_c_v[a][k]))*N_tot_c_v_inv + (x[k*ag+5*sg+r*rg+3]*alpha_i+x[k*ag+5*sg+r*rg+4]+a1*(x[k*ag+5*sg+r*rg+8]*alpha_i+x[k*ag+5*sg+r*rg+9])+a2*(x[k*ag+5*sg+r*rg+13]*alpha_i+x[k*ag+5*sg+r*rg+14])+a3*(x[k*ag+5*sg+r*rg+18]*alpha_i+x[k*ag+5*sg+r*rg+19]))*(qp*(cnt_matrix_nwn_p_v[a][k]+qc*cnt_matrix_nwn_c_v[a][k]))*N_tot_n_v_inv;
              }
          }

          if (direct)
          {
              I_temp_p_v = I_temp_p; I_temp_c_v = I_temp_c; I_temp_n_v = I_temp_n;
          }
          
          int pj = max(ag*(a-1),0);
          int cj = max(ag*(a),0);
          double muBp, muBc, muBn, muBpv, muBcv, muBnv, mu_mat, cl, Icp, Icc, up, kpc, kpd, pro, pro_v_in, pro_v_out, rp;
          double p_vulp, p_vulc, p_vuln, p_vulpv, p_vulcv, p_vulnv;
          double xi_bp, xi_bc, xi_bn, xi_bpv, xi_bcv, xi_bnv;
          double ej1  = eta[a+1];
          double ej   = eta[a];
          double lossP, lossMS0, lossMS1, lossMS2;
          
          // Birth rate into each social group
          // Rcpp::Rcout << "Birth rate into each social group: " << a << std::endl;

          if (a == 0)
          {
              p_vulp = 0; p_vulc = sum_wcb; p_vuln = sum_wcb; p_vulpv = 0; p_vulcv = sum_wcb_v; p_vulnv = sum_wcb_v;
              
              muBp = 0;    muBc = dailyBirthRate*(cov_c);   muBn = dailyBirthRate*(1-cov_c);
              muBpv = 0;  muBcv = dailyBirthRate*(cov_c);  muBnv = dailyBirthRate*(1-cov_c);
              xi_bp = 1;    xi_bc = 1;   xi_bn = 1;
              xi_bpv = 1;   xi_bcv = xi_boost;  xi_bnv = 1;
          }
          else
          {
              p_vulp = 0; p_vulc = 0; p_vuln = 0; p_vulpv = 0; p_vulcv = 0; p_vulnv = 0;
              muBp = 0; muBc = 0; muBn = 0;
              muBpv = 0; muBcv = 0; muBnv = 0;
              xi_bp = 1;    xi_bc = 1;   xi_bn = 1;
              xi_bpv = 1;   xi_bcv = 1;  xi_bnv = 1;
          }
          

          for(int s = 0; s < 6; s++)
          {
              kpd = kpc = pro = pro_v_in = pro_v_out = 0;
              cl = 0;     Icp = 1;    Icc = 0;    up = 0;
              
              double u, In, mu, vac_c_o, vac_c_i;
              int kp;
              if (a < 12)
              {
                  if (s == 0)
                  {   kp = 0*sg;       u = 0;             In = I_temp_p;   mu = muBp;  mu_mat = 0; p_vul = p_vulp; xi_b = xi_bp;
                  }
                  else if (s == 1)
                  {   kp = 1*sg;       u = 1;       In = I_temp_c;   mu = muBc;  mu_mat = 0; p_vul = p_vulc; xi_b = xi_bc;
                  }
                  else if (s == 2)
                  {   kp = 2*sg;       u = 1;     In = I_temp_n;   mu = muBn;  mu_mat = 0; p_vul = p_vuln; xi_b = xi_bn;
                  }
                  else if (s == 3)
                  {   kp = 3*sg;       u = 0;             In = I_temp_p_v; mu = muBpv; mu_mat = 0; p_vul = p_vulpv; xi_b = xi_bpv;
                  }
                  else if (s == 4)
                  {   kp = 4*sg;       u = 1;       In = I_temp_c_v; mu = muBcv*(1-vac_cal(t1,a)); mu_mat = muBcv*vac_cal(t1,a); p_vul = p_vulcv;  xi_b = xi_bcv;
                  }
                  else
                  {   kp = 5*sg;       u = 1;     In = I_temp_n_v; mu = muBnv; mu_mat = 0; p_vul = p_vulnv;  xi_b = xi_bnv;
                  }
              }
              else
              {
                  if (s == 0)
                  {   kp = 0*sg;       u = p_mat[a]*(1-cov_c);        In = I_temp_p;   mu = mu_mat = 0;  p_vul = p_vulp;
                  }
                  else if (s == 1)
                  {   kp = 1*sg;       u = p_mat[a]*(cov_c);          In = I_temp_c;   mu = mu_mat = 0;  p_vul = p_vulc;
                  }
                  else if (s == 2)
                  {   kp = 2*sg;       u = (1-p_mat[a]);              In = I_temp_n;   mu = mu_mat = 0;  p_vul = p_vuln;
                  }
                  else if (s == 3)
                  {   kp = 3*sg;       u = p_mat[a]*(1-cov_c);        In = I_temp_p_v; mu = mu_mat = 0; p_vul = p_vulpv;
                  }
                  else if (s == 4)
                  {   kp = 4*sg;       u = p_mat[a]*(cov_c);          In = I_temp_c_v; mu = mu_mat = 0; p_vul = p_vulcv;
                  }
                  else
                  {   kp = 5*sg;       u = (1-p_mat[a]);             In = I_temp_n_v; mu = mu_mat = 0; p_vul = p_vulnv;
                  }
              }
              
                for (int r = 0; r < 3; r++)
                {
                    // Sort out maternal vaccination protection in infants
                    if ((a == 0) & (s == 4)){   
                        kp = 4*sg; 
                        if (r == 0)  {
                            mu = muBcv; mu_mat = 0;
                        } else {
                            mu = muBcv * (1-vac_cal(t1, 0)); mu_mat = muBcv * vac_cal(t1, 0);
                        }
                    } else if ((a == 0) & (s == 5)) {
                        kp = 5*sg; 
                        if (r == 0)  {
                            mu = muBnv; mu_mat = 0;
                        } else {
                            mu = muBnv; mu_mat = 0;
                        }
                    } 
                  //   Rcpp::Rcout << "Risk groups parameters: " << r << std::endl;

                  double PST = 0;
                  step_out = 0; step_in = 0; step_in_S0 = 0;  step_in_S1 = 0; 
                  x_tot = x_tot_1 = x_tot_2 = 0; x_tot_3 = 0;
                  if (r==0){
                      rp = pVHR[a];}
                  else if (r==1){
                      rp = pHR[a];}
                  else if (r==2){
                      rp = pLR[a];}
                  else
                      cout << "ERROR" << endl;
                  
                  if (r == 2)
                  {
                   //   Rcpp::Rcout << "Place where r == 2." << std::endl;
                      if (a==18 || a==19 || a==20)
                      {
                          if (a==18)
                              up = u18p;
                          else if (a==19)
                              up = u19p;
                          else if (a==20)
                              up = u20p;
                          else
                              cout << "ERROR" << endl;
                      }
                  }
              
                  if (s < 3)
                  {
                    //  Rcpp::Rcout << "Place where s < 3." << std::endl;
                      if (a < 12)
                      {
                          // AGING first 12 months, all risk and social groups are pooled
                            for (int i = 0; i < 21; i++)
                                PS[i] = x[pj+s*sg+0*rg+i] + x[pj+s*sg+1*rg+i] + x[pj+s*sg+2*rg+i];
                                
                            for (int i = 21; i < rg; i++)
                                PS[i] = x[pj+s*sg+0*rg+i] + x[pj+s*sg+1*rg+i] + x[pj+s*sg+2*rg+i];
                      }
                      else
                      {
                          // AGING, >1 years, all risk and social groups are pooled
                          for (int i = 0; i < rg; i++)
                                PS[i] = (x[pj+0*sg+0*rg+i]+x[pj+0*sg+1*rg+i]+x[pj+0*sg+2*rg+i] + x[pj+1*sg+0*rg+i]+x[pj+1*sg+1*rg+i]+x[pj+1*sg+2*rg+i] + x[pj+2*sg+0*rg+i]+x[pj+2*sg+1*rg+i]+x[pj+2*sg+2*rg+i]);
                      }
                  }
                  else
                  {
                        if (a < 12)
                        {
                            // AGING first 12 months, all risk and social groups are pooled
                            for (int i = 0; i < 21; i++)
                                PS[i] = x[pj+s*sg+0*rg+i] + x[pj+s*sg+1*rg+i] + x[pj+s*sg+2*rg+i]; //x[pj+3*sg+0*rg+i] + x[pj+3*sg+1*rg+i] + x[pj+3*sg+2*rg+i] + x[pj+4*sg+0*rg+i] + x[pj+4*sg+1*rg+i] + x[pj+4*sg+2*rg+i] + x[pj+5*sg+0*rg+i] + x[pj+5*sg+1*rg+i] + x[pj+5*sg+2*rg+i];
                                
                            // AGING, all social groups are pooled
                            for (int i = 21; i < rg; i++) 
                                PS[i] = x[pj+s*sg+0*rg+i] + x[pj+s*sg+1*rg+i] + x[pj+s*sg+2*rg+i]; 
                      }
                      else
                      {
                          // AGING, >1 years, all risk and social groups are pooled
                            for (int i = 0; i < rg; i++)
                                PS[i] = x[pj+3*sg+0*rg+i] + x[pj+3*sg+1*rg+i] + x[pj+3*sg+2*rg+i] + x[pj+4*sg+0*rg+i] + x[pj+4*sg+1*rg+i] + x[pj+4*sg+2*rg+i] + x[pj+5*sg+0*rg+i] + x[pj+5*sg+1*rg+i] + x[pj+5*sg+2*rg+i];
                      }
                  }
                  
                    double cpmu, cpo, cMmu, cMo, cLo, cMa, cpmu_dose, cpo_dose, cMmu_dose, cMo_dose, cLo_dose, cMa_dose;

                    int p = a*ag + s*sg + r*rg;
                    int q = a*ag + kpc + r*rg;
                    int o = 0;
                    if (s > 3)
                        o = a*ag + (s-3)*sg + r*rg;
                  
                    if (s < 3)
                    {
                        step_out = 0; step_in = 0; step_in_S0 = 0; step_in_S1 = 0;
                        cpmu = 0; cpmu_dose = 0;
                        cpo = 0; cpo_dose = 0;
                        cMmu = 0; cMmu_dose = 0;
                        cMo = 0; cMo_dose = 0;
                        cLo = 0; cLo_dose = 0;
                        cMa = 0; cMa_dose = 0;
                        SN = 1;
                        lossMS2 = lossMS1 = lossMS0 = lossP = 0;
                        loss_mat = loss_lav_S0 = loss_lav_S1 = loss_lav_S2 = loss_lav_S3 = lossMS0 = lossMS1 = lossMS2 = 0;
                    }
                    else
                    {
                        // Get the waning rates for each of the categories 
                        if (r == 0) {
                            loss_mat = loss_lav_S0 = loss_lav_S1 = loss_lav_S2 = loss_lav_S3 = 0;
                            lossMS0 = x[p + 26]*(eff_wane_vhr);
                            lossMS1 = x[p + 29]*(eff_wane_vhr);
                            lossMS2 = x[p + 32]*(eff_wane_vhr); 
                        } else if (r == 2){
                            loss_mat = x[p + 21]*(eff_wane_mat);
                            lossMS0 = x[p + 26]*(eff_wane_mab);
                            lossMS1 = x[p + 29]*(eff_wane_mab);
                            lossMS2 = x[p + 32]*(eff_wane_mab); 
                            loss_lav_S0 = x[p + 35]*(eff_wane_lav); 
                            loss_lav_S1 = x[p + 38]*(eff_wane_lav); 
                            loss_lav_S2 = x[p + 41]*(eff_wane_lav); 
                            loss_lav_S3 = x[p + 44]*(eff_wane_lav); 

                        } else if (r == 1) {
                            loss_mat = loss_lav_S0 = loss_lav_S1 = loss_lav_S2 = loss_lav_S3 = lossMS0 = lossMS1 = lossMS2 = 0;
                        }

                        step_out = 0; step_in = 0; step_in_S0 = 0; step_in_S1 = 0;
                        cpmu = 0; cpmu_dose = 0;
                        cpo = 0; cpo_dose = 0;
                        cMmu = 0; cMmu_dose = 0;
                        cMo = 0; cMo_dose = 0;
                        cLo = 0; cLo_dose = 0;
                        cMa = 0; cMa_dose = 0;

                        if (a == 0){
                            M = x[a*ag + (s-3)*sg + r*rg + 0]; 
                            S0 = x[a*ag + (s-3)*sg + r*rg + 1]; S1 = x[a*ag + (s-3)*sg + r*rg + 6]; S2 = x[a*ag + (s-3)*sg + r*rg + 11]; S3 = x[a*ag + (s-3)*sg + r*rg + 16];

                            if (r == 0){
                                if( s == 5) {
                                    step_in_S0 = 0; // collect_protect(t, 0);
                                    step_in_S1 = 0; //collect_protect(t, 1);
                                    step_out_dose += cal_mAB_VHR_step_out_dose(t1, 0);
                                }
                                step_out = cal_mAB_VHR_step_out(t1, 0);
                                
                                cMo = cal_mAB_VHR(t1, 0);
                                cMo_dose = cal_mAB_VHR_dose(t1, 0);
                            } else if (r == 1){
                                if( s == 5) {
                                    step_in_S0 = 0; //collect_protect(t, 2);
                                    step_in_S1 = 0;  //collect_protect(t, 3);
                                    step_out_dose += cal_mAB_HR_step_out_dose(t1, 0);
                                }
                                step_out = cal_mAB_HR_step_out(t1, 0);
                                
                                cMo = cal_mAB_HR(t1, 0);
                                cMo_dose = cal_mAB_HR_dose(t1, 0);
                            } else {
                                if( s == 5) {
                                    step_in_S0 = collect_protect(t, 4) + collect_protect(t, 0) + collect_protect(t, 2);
                                    step_in_S1 = collect_protect(t, 5) + collect_protect(t, 1) + collect_protect(t, 3);
                                    step_out_dose += cal_mAB_LR_step_out_dose(t1, 0);
                                }
                                step_out = cal_mAB_LR_step_out(t1, 0);
                                cMo = cal_mAB_LR(t1, 0);
                                cMo_dose = cal_mAB_LR_dose(t1, 0);
                            }
                        }
                        else if (a > 0) {

                            M = x[a*ag + (s-3)*sg + r*rg + 0]; 

                            if (r == 0){
                                if( s == 5) {
                                    step_in_S0 = 0; // collect_protect(t, 6 * (a) + 0);
                                    step_in_S1 = 0; // collect_protect(t, 6 * (a) + 1);
                                    step_out_dose += cal_mAB_VHR_step_out_dose(t1, a);
                                }
                                step_out = cal_mAB_VHR_step_out(t1, a);
                                
                                cMo = cal_mAB_VHR(t1, a);
                                cMo_dose = cal_mAB_VHR_dose(t1, a);
                                S0 = x[a*ag + (s-3)*sg + r*rg + 1]; S1 = x[a*ag + (s-3)*sg + r*rg + 6]; S2 = x[a*ag + (s-3)*sg + r*rg + 11]; S3 = x[a*ag + (s-3)*sg + r*rg + 16];
                                SN = 1;
                            }
                            else if (r == 1){
                                if( s == 5) {
                                    step_in_S0 = 0; // collect_protect(t, 6 * (a) + 2);
                                    step_in_S1 = 0; // collect_protect(t, 6 * (a) + 3);
                                    step_out_dose += cal_mAB_HR_step_out_dose(t1, a);
                                }
                                step_out = cal_mAB_HR_step_out(t1, a);
                                
                                cMo = cal_mAB_HR(t1, a); // Rcpp::Rcout << "A" << endl;
                                cLo = cal_LAV_HR(t1, a);// Rcpp::Rcout << "B" << endl;
                                cMo_dose = cal_mAB_HR_dose(t1, a);// Rcpp::Rcout << "C" << endl;
                                cLo_dose = cal_LAV_HR_dose(t1, a);// Rcpp::Rcout << "D" << endl;
                                S0 = x[a*ag + (s-3)*sg + r*rg + 1]; S1 = x[a*ag + (s-3)*sg + r*rg + 6]; S2 = x[a*ag + (s-3)*sg + r*rg + 11]; S3 = x[a*ag + (s-3)*sg + r*rg + 16]; //Rcpp::Rcout << "E" << endl;
                                SN = 1;// Rcpp::Rcout << "F" << endl;
                            }
                            else
                            {
                                if( s == 5) {
                                    step_in_S0 = collect_protect(t, 6 * (a) + 4) + collect_protect(t, 6 * (a) + 0) + collect_protect(t, 6 * (a) + 2);
                                    step_in_S1 = collect_protect(t, 6 * (a) + 5) + collect_protect(t, 6 * (a) + 1) + collect_protect(t, 6 * (a) + 3);
                                    step_out_dose += cal_mAB_LR_step_out_dose(t1, a);
                                }
                                    step_out = cal_mAB_LR_step_out(t1, a);  
                                    cMo = cal_mAB_LR(t1, a);
                                    cMo_dose = cal_mAB_LR_dose(t1, a);
                                if (s == 4)
                                {
                                    cLo = vac_cal(t1, a);
                                    cLo_dose = vac_cal_dose(t1, a);
                                    S0 = x[a*ag + (s-3)*sg + r*rg + 1]; S1 = x[a*ag + (s-3)*sg + r*rg + 6]; S2 = x[a*ag + (s-3)*sg + r*rg + 11]; S3 = x[a*ag + (s-3)*sg + r*rg + 16];
                                }
                                else
                                {
                                    SN = 1;
                                    cLo = cal_LAV_LR(t1, a);
                                    cLo_dose = cal_LAV_LR_dose(t1, a);
                                    S0 = x[a*ag + (s-3)*sg + r*rg + 1]; S1 = x[a*ag + (s-3)*sg + r*rg + 6]; S2 = x[a*ag + (s-3)*sg + r*rg + 11]; S3 = x[a*ag + (s-3)*sg + r*rg + 16];
                              }
                          }
                      }
                      else
                      {
                          
                      }
                  }
                  
                  for (int i = 0; i < 21; i++)
                  {
                      x_tot += x[o+i];
                      PST += PS[i];
                      if (PST < 1)
                          PST = 1;
                  }
                  
              //    Rcpp::Rcout << "keeping track of doses below" << std::endl;

                  protectpal += cpo_dose*x_tot;
                  protectmabs += cMo_dose*x_tot;
                  protectLAV += cLo_dose*x_tot;
                  protectmat += cMa_dose*x_tot;
                  
                  if (s == 5) {
                    if (r == 0) {
                        collect_protect(t + eff_wane_step, 6 * (a + 5) + 0) += x_tot_1 * step_out;
                        collect_protect(t + eff_wane_step, 6 * (a + 5) + 1) += x_tot_2 * step_out;
                    } else if (r == 1) {
                        collect_protect(t + eff_wane_step, 6 * (a + 5) + 2) += x_tot_1 * step_out;
                        collect_protect(t + eff_wane_step, 6 * (a + 5) + 3) += x_tot_2 * step_out;
                    } else if (r == 2) {
                        collect_protect(t + eff_wane_step, 6 * (a + 5) + 4) += x_tot_1 * step_out;
                        collect_protect(t + eff_wane_step, 6 * (a + 5) + 5) += x_tot_2 * step_out;
                    }
                  }

                    VectorXd mab_trans = VectorXd::Zero(12);
                    VectorXd lav_trans = VectorXd::Zero(4);

                    // Check for truncation values for Monoclonal antibodies
                    for (int i = 0; i < 12; i++)
                    {
                        double state = x[p+i];
                        double prot = x[o+i]*(cMo);
                        if ((state - prot) < 0) {
                            mab_trans[i] =  x[p+i]*(cMo);
                        } else {
                            mab_trans[i] = prot;
                        }
                    }
                    for (int i = 0; i < 2; i++)
                        x_tot_1 += mab_trans[i];
                  
                    for (int i = 2; i < 7; i++)
                        x_tot_2 += mab_trans[i];

                    for (int i = 7; i < 12; i++)
                        x_tot_3 += mab_trans[i];

                    // Check for truncation values for LAV antibodies
                    vector<double > pos = {1, 6, 11, 16};
                    for (int i = 0; i < 4; i++) {
                        double state = x[p+pos[i]];
                        double prot = x[o+pos[i]]*(cLo);
                        if ((state - prot) < 0) {
                            lav_trans[i] =  x[p+pos[i]]*(cLo);
                        } else {
                            lav_trans[i] = prot;
                        }
                    }

                  // Need to add a truncation mechanism to ensure that x[p+i] is always greater than 0

                  dxdt[p+0] = (1.0-p_vul)*mu*rp - x[p+0]*xi*xi_b - (x[p+0])*ej1 + PS[0]*ej*rp*u - mab_trans[0]; 
                  
                  dxdt[p+1] = p_vul*mu*rp  + x[p+0]*xi*xi_b + lossMS0 - x[p+1]*In*beta - (x[p+1])*ej1 + PS[1]*ej*rp*u - mab_trans[1] - x[o+1]*step_out - lav_trans[0] + step_in_S0 + x[p+23]*(eff_wane_mat);
                  dxdt[p+2] = x[p+1]*In*beta                 - x[p+2]*si          - (x[p+2])*ej1 + PS[2]*ej*rp*u -mab_trans[2] - x[p+2]*step_out;
                  dxdt[p+3] = x[p+2]*si*pA[a]                - x[p+3]*ga0*rho     -  (x[p+3])*ej1 + PS[3]*ej*rp*u - mab_trans[3] - x[p+3]*step_out;
                  dxdt[p+4] = x[p+2]*si*(1.0-pA[a])          - x[p+4]*ga0         - (x[p+4])*ej1 + PS[4]*ej*rp*u - mab_trans[4] - x[p+4]*step_out;
                  dxdt[p+5] = x[p+4]*ga0 + x[p+3]*ga0*rho    - x[p+5]*om          - (x[p+5])*ej1 + PS[5]*ej*rp*u - mab_trans[5] - x[p+5]*step_out;
                  
                  dxdt[p+6] = x[p+5]*om                      - d1*x[p+6]*In*beta + lossMS1  - (x[p+6])*ej1 + PS[6]*ej*rp*u - mab_trans[6] + loss_lav_S0 - lav_trans[1] + step_in_S1;
                  dxdt[p+7] = d1*x[p+6]*In*beta              - x[p+7]*si          - (x[p+7])*ej1 + PS[7]*ej*rp*u - mab_trans[7];
                  dxdt[p+8] = x[p+7]*si*pA[a]                - x[p+8]*ga1*rho     - (x[p+8])*ej1 + PS[8]*ej*rp*u - mab_trans[8];
                  dxdt[p+9] = x[p+7]*si*(1.0-pA[a])          - x[p+9]*ga1         - (x[p+9])*ej1 + PS[9]*ej*rp*u - mab_trans[9];
                  dxdt[p+10] = x[p+9]*ga1 + x[p+8]*ga1*rho   - x[p+10]*om         - (x[p+10])*ej1 + PS[10]*ej*rp*u - mab_trans[10];
                  
                  dxdt[p+11] = x[p+10]*om                    - d2*x[p+11]*In*beta + lossMS2 - (x[p+11])*ej1 + PS[11]*ej*rp*u - mab_trans[11] + loss_lav_S1 - lav_trans[2];
                  dxdt[p+12] = d2*x[p+11]*In*beta            - x[p+12]*si         - (x[p+12])*ej1 + PS[12]*ej*rp*u;
                  dxdt[p+13] = x[p+12]*si*pA[a]              - x[p+13]*ga2*rho    - (x[p+13])*ej1 + PS[13]*ej*rp*u;
                  dxdt[p+14] = x[p+12]*si*(1.0-pA[a])        - x[p+14]*ga2        - (x[p+14])*ej1 + PS[14]*ej*rp*u;
                  dxdt[p+15] = x[p+14]*ga2 + x[p+13]*ga2*rho - x[p+15]*om         - (x[p+15])*ej1 + PS[15]*ej*rp*u;
                  
                  dxdt[p+16] = x[p+15]*om + x[p+20]*om       - d3*x[p+16]*In*beta - (x[p+16])*ej1 + PS[16]*ej*rp*u  + loss_lav_S2 + loss_lav_S3 - lav_trans[3];
                  dxdt[p+17] = d3*x[p+16]*In*beta            - x[p+17]*si         - (x[p+17])*ej1 + PS[17]*ej*rp*u;
                  dxdt[p+18] = x[p+17]*si*pA[a]              - x[p+18]*ga3*rho    - (x[p+18])*ej1 + PS[18]*ej*rp*u;
                  dxdt[p+19] = x[p+17]*si*(1.0-pA[a])        - x[p+19]*ga3        - (x[p+19])*ej1 + PS[19]*ej*rp*u;
                  dxdt[p+20] = x[p+19]*ga3 + x[p+18]*ga3*rho - x[p+20]*om         - (x[p+20])*ej1 + PS[20]*ej*rp*u;
                  
                  // Maternal vaccination
                  dxdt[p+21] = mu_mat*rp - x[p+21]*eff_wane_mat - x[p+21]*ej1 + PS[21]*ej*u*rp;  
                  dxdt[p+22] = x[p+21]*eff_wane_mat - x[p+22]*eff_wane_mat - x[p+22]*ej1 + PS[22]*ej*u*rp;  
                  dxdt[p+23] = x[p+22]*eff_wane_mat - x[p+23]*eff_wane_mat - x[p+23]*ej1 + PS[23]*ej*u*rp;  

                  // Monoclonal protection
                  dxdt[p+24] = x_tot_1 - x[p+24]*eff_wane_mab  - x[p+24]*ej1 + PS[24]*ej*u*rp;  // Monoclonal protection for exposure 1 (group 1, erlang-3 distributed)
                  dxdt[p+25] = x[p+24]*eff_wane_mab - x[p+25]*eff_wane_mab  - x[p+25]*ej1 + PS[25]*ej*u*rp;  // Monoclonal protection for exposure 1 (group 2, erlang-3 distributed)
                  dxdt[p+26] = x[p+25]*eff_wane_mab - lossMS0 - x[p+26]*ej1 + PS[26]*ej*u*rp;   // Monoclonal protection for exposure 1 (group 3, erlang-3 distributed)

                  dxdt[p+27] = x_tot_2 - x[p+27]*eff_wane_mab - x[p+27]*ej1 + PS[27]*ej*u*rp; // Monoclonal protection for exposure 2 (group 1, erlang distributed)
                  dxdt[p+28] = x[p+27]*eff_wane_mab - x[p+28]*eff_wane_mab - x[p+28]*ej1 + PS[28]*ej*u*rp;  // Monoclonal protection for exposure 1 (group 1, erlang distributed)
                  dxdt[p+29] = x[p+28]*eff_wane_mab - lossMS1 - x[p+29]*ej1 + PS[29]*ej*u*rp; // Monoclonal protection for exposure 2 (group 2, erlang distributed)

                  dxdt[p+30] = x_tot_3 - x[p+30]*eff_wane_mab - x[p+30]*ej1 + PS[30]*ej*u*rp; // Monoclonal protection for exposure 3 (group 1, erlang distributed)
                  dxdt[p+31] = x[p+30]*eff_wane_mab - x[p+31]*eff_wane_mab - x[p+31]*ej1 + PS[31]*ej*u*rp;  // Monoclonal protection for exposure 1 (group 1, erlang distributed)
                  dxdt[p+32] = x[p+31]*eff_wane_mab - lossMS2 - x[p+32]*ej1 + PS[32]*ej*u*rp;  // Monoclonal protection for exposure 3 (group 2, erlang distributed)
      
                   // LAV protection
                  dxdt[p+33] = lav_trans[0] - x[p+33]*eff_wane_lav - x[p+33]*ej1 + PS[33]*ej*u*rp; // LAV protection for exposure 1 (group 1, erlang-3 distributed)
                  dxdt[p+34] = x[p+33]*eff_wane_lav - x[p+34]*eff_wane_mab - x[p+34]*ej1 + PS[34]*ej*u*rp;  // Monoclonal protection for exposure 1 (group 2, erlang-3 distributed)
                  dxdt[p+35] = x[p+34]*eff_wane_lav - loss_lav_S0 - x[p+35]*ej1 + PS[35]*ej*u*rp; // LAV protection for exposure 1 (group 3, erlang-3 distributed)

                  dxdt[p+36] = lav_trans[1] - x[p+36]*eff_wane_lav - x[p+36]*ej1 + PS[36]*ej*u*rp; // LAV protection for exposure 2 (group 1, erlang-3 distributed)
                  dxdt[p+37] = x[p+36]*eff_wane_lav - x[p+37]*eff_wane_mab - x[p+37]*ej1 + PS[37]*ej*u*rp;  // Monoclonal protection for exposure 1 (group 2, erlang-3 distributed)
                  dxdt[p+38] = x[p+37]*eff_wane_lav - loss_lav_S1 - x[p+38]*ej1 + PS[38]*ej*u*rp; // LAV protection for exposure 2 (group 3, erlang-3 distributed)

                  dxdt[p+39] = lav_trans[2] - x[p+39]*eff_wane_lav - x[p+39]*ej1 + PS[39]*ej*u*rp; // LAV protection for exposure 3 (group 1, erlang-3 distributed)
                  dxdt[p+40] = x[p+39]*eff_wane_lav - x[p+40]*eff_wane_mab - x[p+40]*ej1 + PS[40]*ej*u*rp;  // Monoclonal protection for exposure 1 (group 2, erlang-3 distributed)
                  dxdt[p+41] = x[p+40]*eff_wane_lav - loss_lav_S2 - x[p+41]*ej1 + PS[41]*ej*u*rp; // LAV protection for exposure 3 (group 3, erlang-3 distributed)

                  dxdt[p+42] = lav_trans[3] - x[p+42]*eff_wane_lav - x[p+42]*ej1 + PS[42]*ej*u*rp; // LAV protection for exposure 4 (group 1, erlang-3 distributed)
                  dxdt[p+43] = x[p+42]*eff_wane_lav - x[p+43]*eff_wane_mab - x[p+43]*ej1 + PS[43]*ej*u*rp;  // Monoclonal protection for exposure 1 (group 2, erlang-3 distributed)
                  dxdt[p+44] = x[p+43]*eff_wane_lav - loss_lav_S3 - x[p+44]*ej1 + PS[44]*ej*u*rp; // LAV protection for exposure 4 (group 3, erlang-3 distributed)
                }
            }
          // Vaccine groups
          dxdt[a*ag + 6*sg + 0] = si*(x[cj+3*sg+0*rg+2] + x[cj+4*sg+0*rg+2] + x[cj+5*sg+0*rg+2] + x[cj+3*sg+1*rg+2] + x[cj+4*sg+1*rg+2] + x[cj+5*sg+1*rg+2] + x[cj+3*sg+2*rg+2] + x[cj+4*sg+2*rg+2] + x[cj+5*sg+2*rg+2]);
          dxdt[a*ag + 6*sg + 1] = protectpal;
          dxdt[a*ag + 6*sg + 2] = protectmabs + step_out_dose;
         // cout << protectmabs << endl;
          /*dxdt[a*ag + 6*sg + 2] = protectmabs + vac_cal(t1, a)*(x[a*ag+4*sg+2*rg+1]+x[a*ag+4*sg+2*rg+6]+x[a*ag+4*sg+2*rg+11]+x[a*ag+4*sg+2*rg+16])
          + (cal_LAV_HR(t1, a))*((x[a*ag+0*sg+1*rg+1]+x[a*ag+0*sg+1*rg+6]+x[a*ag+0*sg+1*rg+11]+x[a*ag+0*sg+1*rg+16]) + (x[a*ag+1*sg+1*rg+1]+x[a*ag+1*sg+1*rg+6]+x[a*ag+1*sg+1*rg+11]+x[a*ag+1*sg+1*rg+16]))
          + (cal_LAV_LR(t1, a))*((x[a*ag+0*sg+2*rg+1]+x[a*ag+0*sg+2*rg+6]+x[a*ag+0*sg+2*rg+11]+x[a*ag+0*sg+2*rg+16]) + (x[a*ag+1*sg+2*rg+1]+x[a*ag+1*sg+2*rg+6]+x[a*ag+1*sg+2*rg+11]+x[a*ag+1*sg+2*rg+16]));
          */
          dxdt[a*ag + 6*sg + 3] = protectLAV;
          
          dxdt[a*ag + 6*sg + 4] = protectmat;
          dxdt[a*ag + 6*sg + 5] = 0;
          dxdt[a*ag + 6*sg + 6] = 0;
          dxdt[a*ag + 6*sg + 7] = 0;
          
          // Monitoring Parent
          dxdt[a*ag + 6*sg + 8] =  si*(x[cj+3*sg+0*rg+2]+x[cj+3*sg+0*rg+7]+x[cj+3*sg+0*rg+12]+x[cj+3*sg+0*rg+17]); //VHR
          dxdt[a*ag + 6*sg + 9] =  si*(x[cj+3*sg+1*rg+2]+x[cj+3*sg+1*rg+7]+x[cj+3*sg+1*rg+12]+x[cj+3*sg+1*rg+17]); //HR
          dxdt[a*ag + 6*sg + 10] = si*(x[cj+3*sg+2*rg+2]+x[cj+3*sg+2*rg+7]+x[cj+3*sg+2*rg+12]+x[cj+3*sg+2*rg+17]); //LR
          
          // Monitoring Cocoon
          dxdt[a*ag + 6*sg + 11] = si*(x[cj+4*sg+0*rg+2]+x[cj+4*sg+0*rg+7]+x[cj+4*sg+0*rg+12]+x[cj+4*sg+0*rg+17]); //VHR
          dxdt[a*ag + 6*sg + 12] = si*(x[cj+4*sg+1*rg+2]+x[cj+4*sg+1*rg+7]+x[cj+4*sg+1*rg+12]+x[cj+4*sg+1*rg+17]); //HR
          dxdt[a*ag + 6*sg + 13] = si*(x[cj+4*sg+2*rg+2]+x[cj+4*sg+2*rg+7]+x[cj+4*sg+2*rg+12]+x[cj+4*sg+2*rg+17]); //LR
          
          // Monitoring Neither
          dxdt[a*ag + 6*sg + 14] = si*(x[cj+5*sg+0*rg+2]+x[cj+5*sg+0*rg+7]+x[cj+5*sg+0*rg+12]+x[cj+5*sg+0*rg+17]); //VHR
          dxdt[a*ag + 6*sg + 15] = si*(x[cj+5*sg+1*rg+2]+x[cj+5*sg+1*rg+7]+x[cj+5*sg+1*rg+12]+x[cj+5*sg+1*rg+17]); //HR
          dxdt[a*ag + 6*sg + 16] = si*(x[cj+5*sg+2*rg+2]+x[cj+5*sg+2*rg+7]+x[cj+5*sg+2*rg+12]+x[cj+5*sg+2*rg+17]); //LR
          
          x[a*ag + 6*sg + 17] = N_tot_n[a];
          x[a*ag + 6*sg + 18] = N_tot_c[a];
          x[a*ag + 6*sg + 19] = N_tot_p[a];
          
          x[a*ag + 6*sg + 20] = N_tot_n_v[a];
          x[a*ag + 6*sg + 21] = N_tot_c_v[a];
          x[a*ag + 6*sg + 22] = N_tot_p_v[a];
      }
    }
};

void check_stability(vector<double>& x0)  {
    for (int i = 0; i < x0.size(); i ++) {
        if (x0[i] < 0) {
            Rcpp::Rcout << "i: " << i << ". x0: " << x0[i] << std::endl;
            x0[i] = 0;
        }
    }
};



List RunInterventions::Sample(List vac_calendar, List vac_dose, double cov_c, List vac_info, VectorXd posteriors)
{
    // Restart values
    this->currentODETime = this->run_start;
    VectorXd inc_tot(this->A);
    inc_tot = VectorXd::Zero(this->A);
    this->weekNo = 0;
    this->dayNoAfterBurn = 0;

    VectorXd currentParamValues = posteriors;
    ParameterValuesforODE(currentParamValues);
    MatrixXd collect_protect = MatrixXd::Zero(365 * 12 + 365, 30 * 6);
     // Set up and Run ODE solver 
    EulerTalt<state_t> integrator;
    SystemTalt<state_t, system_t> System;
    asc::Recorder recorder;
    vector< double > x0 = generateInitialStates(cov_c);
    ODE_desc ODE_desc_inst(this, vac_calendar, vac_dose, vac_info, cov_c);

    NumericMatrix sampleWeeklyIncidence(521, this->A * 9);
    NumericMatrix no_doses(521, 4);

    while (this->currentODETime < (this->run_full + this->run_burn)){

        integrator(ODE_desc_inst, x0, this->currentODETime, this->dt, collect_protect);
        if (this->currentODETime > this->run_burn){

            getWeeklyIncidence(x0, sampleWeeklyIncidence, no_doses, FALSE);
          
        }
    }
    x0.clear();

    return Rcpp::List::create(Rcpp::Named("inci") = sampleWeeklyIncidence,
                              Rcpp::Named("doses") = no_doses);
}

MatrixXd RunInterventions::StatesValues(List vac_calendar, List vac_dose, double cov_c, List vac_info, VectorXd posteriors)
{
  // Restart values
  this->currentODETime = this->run_start;
  VectorXd inc_tot(this->A);
  inc_tot = VectorXd::Zero(this->A);
  this->weekNo = 0;
  this->dayNoAfterBurn = 0;
    
  VectorXd currentParamValues = posteriors;
  ParameterValuesforODE(currentParamValues);
  MatrixXd collect_protect = MatrixXd::Zero(365 * 12 + 365, 30 * 6);
  // Set up and Run ODE solver
    EulerTalt<state_t> integrator;
    SystemTalt<state_t, system_t> System;
    asc::Recorder recorder;
  vector< double > x0 = generateInitialStates(cov_c);
  int x0_N = x0.size();
  MatrixXd xmat = MatrixXd::Zero( 365 * 12 + 365, x0_N);
  ODE_desc ODE_desc_inst(this, vac_calendar, vac_dose, vac_info, cov_c);

    NumericMatrix sampleWeeklyIncidence(521, this->A*9);
    NumericMatrix no_doses(521, 4);


    while (this->currentODETime < (this->run_full + this->run_burn)){
        for (int j = 0; j < x0_N; j++){
            xmat(this->currentODETime, j) = x0[j];
        }
        integrator(ODE_desc_inst, x0, this->currentODETime, this->dt, collect_protect);

      if (this->currentODETime > this->run_burn){

        getWeeklyIncidence(x0, sampleWeeklyIncidence, no_doses, FALSE);
          
      }
    }

    return xmat;
}


#endif /* RunIntervention_h */

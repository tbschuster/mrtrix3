#include "command.h"
#include "dwi/sdeconv/response.h"
#include "dwi/shells.h"


using namespace MR;
using namespace App;

void usage () {

    AUTHOR = "Joe Bloggs";
    SYNOPSIS = "Estimate response function coefficients for a given b-value";
    DESCRIPTION
        + "What the command does"
        + "multiple paragraphs can be added like this";

    ARGUMENTS
        + Argument ("arg", "argument description").type_text();
}



// Create a class to represent the response
//    public attributes of the class:
//      1. Method to retrieve the coeffs at a given b-value
//      2. Method to query the max harmonic order (that it stores)

class Response{
    using BValueList = decltype(std::declval<const Eigen::MatrixXd>().col(0));

    public:
        Eigen::VectorXd retrieve_coeffs(const BValueList bval, const Eigen::Vector3d& dir)
        {
            // this method retrieves the SH coefficients for a given b-value
            Eigen::MatrixXd M, Y;
            Eigen::VectorXd amplitudes, b;

            // coef_vec = S - M
            // M[i] = Y[i]


            // std::distance(x.begin(), std::min_element(x.begin(), x.end()));
        }


        int max_harmonic()b
        {
            
        }
};

void run()
{

}
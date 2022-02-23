#include "command.h"
#include "dwi/sdeconv/response.h"
// #include "dwi/shells.h"


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



void run()
{
    DWI::SDeconv::Response response(argument[0]);

    VAR(response.lmax());
    VAR(response.coeffs(3500));
}
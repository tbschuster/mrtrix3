#include "command.h"
// #include "dwi/sdeconv/response.h"

using namespace MR;
using namespace App;

void usage () {

    AUTHOR = "Joe Bloggs";
    SYNOPSIS = "example minimal command";
    DESCRIPTION
        + "What the command does"
        + "multiple paragraphs can be added like this";

    ARGUMENTS
        + Argument ("arg", "argument description").type_text();
}


void run()
{

}
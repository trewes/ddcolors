#ifndef DDCOLORS_SCIP_INTERFACE_H
#define DDCOLORS_SCIP_INTERFACE_H

#include "DecisionDiagram.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/pub_misc.h"


double compute_flow_solution_SCIP(DecisionDiagram &dd, Model model = IP, int coloring_upper_bound = -1);

double
validate_flow_solution_SCIP(DecisionDiagram &dd, double flow_value, Model model = IP, int coloring_upper_bound = -1);



static const int npoints = 10;
static const unsigned int randseed = 42;

static
SCIP_RETCODE setupProblem(
        SCIP*                 scip,               /**< SCIP data structure */
        SCIP_RANDNUMGEN*      randnumgen          /**< random number generator */
)
{
    SCIP_VAR* a;
    SCIP_VAR* b;
    SCIP_VAR* r;

    char name[SCIP_MAXSTRLEN];
    int i;

    /* create empty problem */
    SCIP_CALL( SCIPcreateProbBasic(scip, "circle") );

    /* create variables and add to problem */
    SCIP_CALL( SCIPcreateVarBasic(scip, &a, "a", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
    SCIP_CALL( SCIPcreateVarBasic(scip, &b, "b", -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
    SCIP_CALL( SCIPcreateVarBasic(scip, &r, "r", 0.0, SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );

    SCIP_CALL( SCIPaddVar(scip, a) );
    SCIP_CALL( SCIPaddVar(scip, b) );
    SCIP_CALL( SCIPaddVar(scip, r) );

    /* create soc constraints, add to problem, and forget about them */
    for( i = 0; i < npoints; ++i )
    {
        SCIP_CONS* cons;
        SCIP_VAR* ab[2];
        SCIP_Real xy[2];

        ab[0] = a;
        ab[1] = b;
        xy[0] = -SCIPrandomGetReal(randnumgen, 1.0, 10.0);
        xy[1] = -SCIPrandomGetReal(randnumgen, 1.0, 10.0);

        (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "point%d", i);
        SCIP_CALL( SCIPcreateConsBasicSOC(scip, &cons, name, 2, ab, NULL, xy, 0.0, r, 1.0, 0.0) );
        SCIP_CALL( SCIPaddCons(scip, cons) );
        SCIP_CALL( SCIPreleaseCons(scip, &cons) );
    }

    /* release variables */
    SCIP_CALL( SCIPreleaseVar(scip, &a) );
    SCIP_CALL( SCIPreleaseVar(scip, &b) );
    SCIP_CALL( SCIPreleaseVar(scip, &r) );

    return SCIP_OKAY;
}

static
SCIP_RETCODE runCircle(void)
{
    SCIP* scip;
    SCIP_RANDNUMGEN* randnumgen;

    SCIP_CALL( SCIPcreate(&scip) );
    SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

    SCIPinfoMessage(scip, NULL, "\n");
    SCIPinfoMessage(scip, NULL, "*********************************************\n");
    SCIPinfoMessage(scip, NULL, "* Running Smallest Enclosing Circle Problem *\n");
    SCIPinfoMessage(scip, NULL, "*********************************************\n");
    SCIPinfoMessage(scip, NULL, "\n");

    /* create random number generator */
    SCIP_CALL( SCIPcreateRandom(scip, &randnumgen, randseed, TRUE) );

    SCIP_CALL( setupProblem(scip, randnumgen) );

    SCIPinfoMessage(scip, NULL, "Original problem:\n");
    SCIP_CALL( SCIPprintOrigProblem(scip, NULL, "cip", FALSE) );

    SCIPinfoMessage(scip, NULL, "\nSolving...\n");
    SCIP_CALL( SCIPsolve(scip) );

    if( SCIPgetNSols(scip) > 0 )
    {
        SCIPinfoMessage(scip, NULL, "\nSolution:\n");
        SCIP_CALL( SCIPprintSol(scip, SCIPgetBestSol(scip), NULL, FALSE) );
    }

    /* free random number generator */
    SCIPfreeRandom(scip, &randnumgen);

    SCIP_CALL( SCIPfree(&scip) );

    return SCIP_OKAY;
}


#endif //DDCOLORS_SCIP_INTERFACE_H

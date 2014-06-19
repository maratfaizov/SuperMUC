#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <papi.h>↲
#include "xread.h"
#include "xwrite.h"

#define NUM_EVENTS 4
/** Peak Performance of one core of superMUC fat node 9600 MFlops*/
#define PEAK_MFLOPS 9600

int error();

int main( int argc, char *argv[] ) {
    if ( argc < 4 ) {
            printf( "Usage: %s <format> <input file> <output prefix>\n", argv[0] );
            return EXIT_FAILURE;
    }
    long_long time_init, time_compute, time_output, time_total;
    int nodeCnt;
    int **points;
    int **elems;

    char *type_flag = argv[1];
    char *file_in = argv[2];
    char *file_out = argv[3];

    int status = 0;
    int isBinary = 0;

    if (strcmp(type_flag, "text") == 0) {
            isBinary = 0;
    } else if (strcmp(type_flag, "bin") == 0) {
            isBinary = 1;
    } else        {
            printf( "Error: <format> must be either 'text' or 'bin'", argv[0] );
            return EXIT_FAILURE;
    }

    /** internal cells start and end index*/
    int nintci, nintcf;
    /** external cells start and end index. The external cells are only ghost cells. They are accessed only through internal cells*/
    int nextci, nextcf;
    /** link cell-to-cell array. Stores topology information*/
    int **lcc;
    /** red-black colouring of the cells*/
    int *nboard;
    /** boundary coefficients for each volume cell */
    double *bs, *be, *bn, *bw, *bl, *bh, *bp, *su;

    /** Event Set #1 for TCM and TCA */
    int retval, EventSet = PAPI_NULL;
    long_long values[NUM_EVENTS];

    /** Event Set #2 for OPS and time */
    int EventSet2 = PAPI_NULL;
    long_long values2[NUM_EVENTS], start_cycles, end_cycles, start_usec, end_usec;

    /** Initialization of the PAPI library */
    retval = PAPI_library_init( PAPI_VER_CURRENT );
    if( retval != PAPI_VER_CURRENT ) {
        fprintf( stderr, "PAPI library initialization error!\n" );
        exit( 1 );
    }

    /** Create a file for counters */
    char filename[255];
    strcpy( filename, "pstats_" );
    strcat( filename, argv[3] );
    strcat( filename, ".dat" );
    FILE *doo;
    doo = fopen( filename, "w" );

    /** Initialization of vtk parameters */
    char output_file_su_vtk[255];
    strcpy( output_file_su_vtk, argv[3] );
    strcat( output_file_su_vtk, "_SU.vtk" );
    char output_file_var_vtk[255];
    strcpy( output_file_var_vtk, argv[3] );
    strcat( output_file_var_vtk, "_VAR.vtk" );
    char output_file_cgup_vtk[255];
    strcpy( output_file_cgup_vtk, argv[3] );
    strcat( output_file_cgup_vtk, "_CGUP.vtk" );

    /* For event sets counting (two sets)*/

    /* Create the Event Set #1 */
    if (PAPI_create_eventset(&EventSet) != PAPI_OK) {
        error(1);
    }

    /* Create the Event Set #2*/
    if (PAPI_create_eventset(&EventSet2) != PAPI_OK) {
        error(1);
    }

    /* Add L2 Total Cache misses to our Event Set */
    if( PAPI_add_event( EventSet, PAPI_L2_TCM ) != PAPI_OK ) {
        error(1);
    }

    /* Add L2 Total Cache accesses to our Event Set */
    if( PAPI_add_event( EventSet, PAPI_L2_TCA ) != PAPI_OK ) {
        error(1);
    }

    /* Add L3 Total Cache misses to our Event Set */
    if( PAPI_add_event( EventSet, PAPI_L3_TCM ) != PAPI_OK ) {
        error(1);
    }

    /* Add L3 Total Cache accesses to our Event Set */
    if( PAPI_add_event( EventSet, PAPI_L3_TCA ) != PAPI_OK ) {
        error(1);
    }

    /* Add Floating Point Operations to our Event Set */
    if( PAPI_add_event( EventSet2, PAPI_FP_OPS ) != PAPI_OK ) {
        error(1);
    }

    int step = 1;
    while( step < 3 ) {
        if ( step == 1 ) {
        /* Start counting events in the Event Set #1 */
        if( PAPI_start( EventSet ) != PAPI_OK ) {
            error( 1 );
            }
        } else if( step == 2 ) {
            /** Stop counting events in the Event Set #1*/
            if( PAPI_stop( EventSet, values ) != PAPI_OK ) {
                error( 1 );
            }

            /* Start counting events in the Event Set #2 */
            if( PAPI_start( EventSet2 ) != PAPI_OK ) {
                printf("3222 OK\n");
                error( 1 );
            }

            /* Gets the starting time in microseconds */
            start_usec = PAPI_get_real_usec();
        }

        /* initialization  */
        // read-in the input file

        int f_status = 0;
        if (isBinary) {
            int f_status = read_binary( file_in, &nintci, &nintcf, &nextci, &nextcf, &lcc,
                                       &bs, &be, &bn, &bw, &bl, &bh, &bp, &su, &nboard );
        } else {
            int f_status = read_formatted( file_in, &nintci, &nintcf, &nextci, &nextcf, &lcc,
                                        &bs, &be, &bn, &bw, &bl, &bh, &bp, &su, &nboard );
        }
        if ( f_status != 0 ) {
            printf( "failed to initialize data!\n" );
            return EXIT_FAILURE;
        }

        // allocate arrays used in gccg
        int nomax = 3;
        /** the reference residual*/
        double resref = 0.0;
        /** the ratio between the reference and the current residual*/
        double ratio;

        /** array storing residuals */
        double* resvec = ( double * ) calloc( sizeof( double ), (nintcf + 1) );
        /** the variation vector -> keeps the result in the end */
        double* var = ( double * ) calloc( sizeof( double ), (nextcf + 1) );

        /** the computation vectors */
        double* direc1 = ( double * ) calloc( sizeof( double ), (nextcf + 1) );
        double* direc2 = ( double * ) calloc( sizeof( double ), (nextcf + 1) );

        /** additional vectors */
        double* cgup = ( double * ) calloc( sizeof( double ), (nextcf + 1) );
        double* oc = ( double * ) calloc( sizeof( double ), (nintcf + 1) );
        double* cnorm = ( double * ) calloc( sizeof( double ), (nintcf + 1) );
        double* adxor1 = ( double * ) calloc( sizeof( double ), (nintcf + 1) );
        double* adxor2 = ( double * ) calloc( sizeof( double ), (nintcf + 1) );
        double* dxor1 = ( double * ) calloc( sizeof( double ), (nintcf + 1) );
        double* dxor2 = ( double * ) calloc( sizeof( double ), (nintcf + 1) );

        // initialize the reference residual
        for ( int nc = nintci; nc <= nintcf; nc++ ) {
            resvec[nc] = su[nc];
            resref = resref + resvec[nc] * resvec[nc];
        }
        resref = sqrt(resref);
        if ( resref < 1.0e-15 ) {
            printf("i/o - error: residue sum less than 1.e-15 - %lf\n", resref);
            return EXIT_FAILURE;
        }

        // initialize the arrays
        for ( int nc = 0; nc <= 10; nc++ ) {
            oc[nc] = 0.0;
            cnorm[nc] = 1.0;
        }

        for ( int nc = nintci; nc <= nintcf; nc++ ) {
            cgup[nc] = 0.0;
            var[nc] = 0.0;
        }

        for ( int nc = nextci; nc <= nextcf; nc++ ) {
            var[nc] = 0.0;
            cgup[nc] = 0.0;
            direc1[nc] = 0.0;
            bs[nc] = 0.0;
            be[nc] = 0.0;
            bn[nc] = 0.0;
            bw[nc] = 0.0;
            bl[nc] = 0.0;
            bh[nc] = 0.0;
        }

        for ( int nc = nintci; nc <= nintcf; nc++ ) {
            cgup[nc] = 1.0 / bp[nc];
        }

        int if1 = 0;
        int if2 = 0;
        int iter = 1;
        int nor = 1;
        int nor1 = nor - 1;
        /* finished initalization */

        if ( step == 1 ) {
            /* Read the counting events in the Event Set #1 */
            if( PAPI_read( EventSet, values ) != PAPI_OK ) {
            error(1);
            }
            fprintf( doo, "INPUT PAPI_L2_TCM %lld\nINPUT PAPI_L2_TCA %lld\n", values[0], values[1]);
            fprintf( doo, "INPUT PAPI_L3_TCM %lld\nINPUT PAPI_L3_TCA %lld\n", values[2], values[3]);
            fprintf( doo, "INPUT L2MissRate %.4f%\n", ( double )values[0] / ( double )values[1] * 100 );
            fprintf( doo, "INPUT L3MissRate %.4f%\n", ( double )values[2] / ( double )values[3] * 100 );

            /* Reset the counting events in the Event Set */
            if (PAPI_reset( EventSet) != PAPI_OK) {
                error(1);
            }
        } else if ( step == 2 ) {
            /* Gets the ending time in microseconds */
            end_usec = PAPI_get_real_usec();

            time_init = end_usec-start_usec;

            /* Read the counting events in the Event Set #2 */
            if( PAPI_read( EventSet2, values2 ) != PAPI_OK ) {
                error(1);
            }

            fprintf(doo, "INPUT ExecutionTime %lld\n", time_init );
            fprintf(doo, "INPUT PAPI_FP_OPS %lld\n", values2[0] );
            fprintf(doo, "INPUT Mflops %.4f\n", ( double )values2[0] / ( double )time_init );
            fprintf(doo, "INPUT Utilization %.4f%\n", ( double )values2[0] / ( double )(time_init * PEAK_MFLOPS) *100 );

            /* Reset the counting events in the Event Set */
            if (PAPI_reset( EventSet2) != PAPI_OK) {
                error(1);
            }
        }

        /* Gets the starting time in microseconds again */
        start_usec = PAPI_get_real_usec();

        /* start computation loop */
        while ( iter < 10000 ) {
            /* start phase 1 */

            // update the old values of direc
            for ( int nc = nintci; nc <= nintcf; nc++ ) {
                    direc1[nc] = direc1[nc] + resvec[nc] * cgup[nc];
            }

            // compute new guess (approximation) for direc
            for ( int nc = nintci; nc <= nintcf; nc++ ) {
                    direc2[nc] = bp[nc] * direc1[nc] - bs[nc] * direc1[lcc[0][nc]]
                                     - bw[nc] * direc1[lcc[3][nc]] - bl[nc] * direc1[lcc[4][nc]]
                                     - bn[nc] * direc1[lcc[2][nc]] - be[nc] * direc1[lcc[1][nc]]
                                     - bh[nc] * direc1[lcc[5][nc]];
            } /* end phase 1 */

            /*  start phase 2 */
            // execute normalization steps
            double oc1, oc2, occ;
           if ( nor1 == 1 ) {
                    oc1 = 0;
                    occ = 0;
                    for ( int nc = nintci; nc <= nintcf; nc++ ) {
                            occ = occ + adxor1[nc] * direc2[nc];
                    }
                    oc1 = occ / cnorm[1];
                    for ( int nc = nintci; nc <= nintcf; nc++ ) {
                            direc2[nc] = direc2[nc] - oc1 * adxor1[nc];
                            direc1[nc] = direc1[nc] - oc1 * dxor1[nc];
                    }
                    if1++;
            } else if ( nor1 == 2 ) {
                oc1 = 0;
                occ = 0;
                    for ( int nc = nintci; nc <= nintcf; nc++ ) {
                            occ = occ + adxor1[nc] * direc2[nc];
                        oc1 = occ / cnorm[1];
                        oc2 = 0;
                        occ = 0;
                        }

                        for ( int nc = nintci; nc <= nintcf; nc++ ) {
                        occ = occ + adxor2[nc] * direc2[nc];

                        oc2 = occ / cnorm[2];
                        }
                        for ( int nc = nintci; nc <= nintcf; nc++ ) {
                                direc2[nc] = direc2[nc] - oc1 * adxor1[nc] - oc2 * adxor2[nc];
                                direc1[nc] = direc1[nc] - oc1 * dxor1[nc] - oc2 * dxor2[nc];
                        }

                        if2++;
                }

                cnorm[nor] = 0;
                double omega = 0;

                // compute the new residual
                for ( int nc = nintci; nc <= nintcf; nc++ ) {
                        cnorm[nor] = cnorm[nor] + direc2[nc] * direc2[nc];
                        omega = omega + resvec[nc] * direc2[nc];
                }
                omega = omega / cnorm[nor];
                omega = omega / cnorm[nor];

                double resnew = 0.0;
                for ( int nc = nintci; nc <= nintcf; nc++ ) {
                        var[nc] = var[nc] + omega * direc1[nc];
                        resvec[nc] = resvec[nc] - omega * direc2[nc];
                        resnew = resnew + resvec[nc] * resvec[nc];
                }
                resnew = sqrt(resnew);
                ratio = resnew / resref;

                // exit on no improvements of residual
                if ( ratio <= 1.0e-10 ) {
                        break;
                }
                iter++;

                // prepare additional arrays for the next iteration step
                if ( nor == nomax ) {
                    nor = 1;
                } else {
                    if ( nor == 1 ) {
                        for ( int nc = nintci; nc <= nintcf; nc++ ) {
                            dxor1[nc] = direc1[nc];
                            adxor1[nc] = direc2[nc];
                        }
                    } else if ( nor == 2 ) {
                        for ( int nc = nintci; nc <= nintcf; nc++ ) {
                            dxor2[nc] = direc1[nc];
                            adxor2[nc] = direc2[nc];
                        }
                    }
                    nor++;
                }
                nor1 = nor - 1;
            }/* end phase 2 */

            /* finished computation loop */

            if ( step == 1 ){
            /* Read the counting events in the Event Set #1 */
                if( PAPI_read( EventSet, values ) != PAPI_OK ) {
                    error(1);
                }

                fprintf( doo, "COMP PAPI_L2_TCM %lld\nCOMP PAPI_L2_TCA %lld\n", values[0], values[1]);
                fprintf( doo, "COMP PAPI_L3_TCM %lld\nCOMP PAPI_L3_TCA %lld\n", values[2], values[3]);

                fprintf( doo, "COMP L2MissRate %.4f%\n", ( double )values[0] / ( double )values[1] * 100);
                fprintf( doo, "COMP L3MissRate %.4f%\n", ( double )values[2] / ( double )values[3] * 100 );

                /* Reset the counting events in the Event Set */
                if (PAPI_reset( EventSet) != PAPI_OK) {
                    error(1);
                }
            } else if( step == 2 ) {
                /* Gets the ending time in microseconds */
                end_usec = PAPI_get_real_usec();

                time_compute = end_usec-start_usec;

                /* Read the counting events in the Event Set #2 */
                if( PAPI_read( EventSet2, values2 ) != PAPI_OK ) {
                    error(1);
                }

                fprintf(doo, "СOMP ExecutionTime %lld\n", time_compute);
                fprintf(doo, "СOMP PAPI_FP_OPS %lld\n", values2[0]);
                fprintf(doo, "СOMP Mflops %.4f\n", ( double )values2[0] / ( double )time_compute);
                fprintf(doo, "COMP Utilization %.4f%\n", ( double)values2[0] / ( double) (time_compute * PEAK_MFLOPS) * 100);

                /* Reset the counting events in the Event Set */
                if (PAPI_reset( EventSet2) != PAPI_OK) {
                    error(1);
                }
            }

            /* Gets the starting time in microseconds */
            start_usec = PAPI_get_real_usec();

            /* write to vtk file */
            if( vol2mesh( nintci, nintcf, lcc, &nodeCnt, &points, &elems ) != 0 ) {
                printf( "failed to convert to unstructured mesh topology and geometry!\n" );
                return EXIT_FAILURE;
            }
            if( write_result_vtk( output_file_su_vtk, nintci, nintcf, nodeCnt, points, elems, su ) != 0 ) {
                printf( "error when trying to write to file %s\n", output_file_su_vtk );
            }
            if( write_result_vtk( output_file_var_vtk, nintci, nintcf,
                                  nodeCnt, points, elems, var ) != 0 ) {
                printf( "error when trying to write to file %s\n", output_file_var_vtk );
            }
            if( write_result_vtk( output_file_cgup_vtk, nintci, nintcf,
                                  nodeCnt, points, elems, cgup ) != 0 ) {
                printf( "error when trying to write to file %s\n", output_file_cgup_vtk );
            }
            if( step == 2 ) {
                /* write output file  */
                if ( write_result(file_in, file_out, nintci, nintcf, var, iter, ratio) != 0 ){
                    printf( "error when trying to write to file %s\n", file_out );
                }
            }

            if ( step == 1 ) {
            /* Read the counting events in the Event Set #1 */
                if( PAPI_read( EventSet, values ) != PAPI_OK ) {
                    error(1);
                }

                fprintf( doo, "OUTPUT PAPI_L2_TCM %lld\nOUTPUT PAPI_L2_TCA %lld\n", values[0], values[1]);
                fprintf( doo, "OUTPUT PAPI_L3_TCM %lld\nOUTPUT PAPI_L3_TCA %lld\n", values[2], values[3]);
                fprintf( doo, "OUTPUT L2MissRate %.4f%\n", ( double)values[0] / ( double )values[1] * 100);
                fprintf( doo, "OUTPUT L3MissRate %.4f%\n", ( double)values[2] / ( double )values[3] * 100);

                /* Reset the counting events in the Event Set */
                if (PAPI_reset( EventSet) != PAPI_OK) {
                    error(1);
                }
            } else  if( step == 2 ) {
                /* Gets the ending time in microseconds */
                end_usec = PAPI_get_real_usec();

                time_output = end_usec-start_usec;

                time_total = time_init+time_compute+time_output;

                /* Read the counting events in the Event Set #2 */
                if( PAPI_read( EventSet2, values2 ) != PAPI_OK ) {
                    error(1);
                }

                fprintf(doo, "OUTPUT ExecutionTime %lld\n", time_output, values2[0]);
                fprintf(doo, "OUTPUT PAPI_FP_OPS %lld\n", values2[0]);
                fprintf(doo, "OUTPUT Mflops %.4f\n", ( double )values2[0] / ( double )time_output);
                fprintf(doo, "OUTPUT Utilization %.4f%\n", (double)values2[0]/(time_output*PEAK_MFLOPS)*100);

                fprintf(doo, "TOTAL ExecutionTime %lld\n", time_total);
                fprintf(doo, "INTPUT Input/TotalTime %.4f%\n", ( double ) time_init/ time_total *100);
                fprintf(doo, "COMP Comput/TotalTime %.4f%\n", ( double ) time_compute/ time_total * 100);
                fprintf(doo, "OUTPUT Output/TotalTime %.4f%\n", ( double ) time_output/ time_total * 100);

                /* Reset the counting events in the Event Set */
                if (PAPI_reset( EventSet2) != PAPI_OK) {
                    error(1);
                }
            }

            /* Free all the dynamically allocated memory */
            free( direc2 );
            free( direc1 );
            free( dxor2 );
            free( dxor1 );
            free( adxor2 );
            free( adxor1 );
            free( cnorm );
            free( oc );
            free( var );
            free( cgup );
            free( resvec );
            free( su );
            free( bp );
            free( bh );
            free( bl );
            free( bw );
            free( bn );
            free( be );
            free( bs );

            step++;
        }

        printf( "Simulation completed successfully!\n" );

        return EXIT_SUCCESS;
    }

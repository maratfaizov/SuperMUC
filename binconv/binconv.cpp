#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string.h>

int main( int argc, char *argv[] ){
    if ( argc < 2 ){
        printf( "No proper input and output filename given" );
        return -1;
    }

    char *input_filename = argv[1];
    char *output_filename = argv[2];

    FILE *f_in = fopen(input_filename, "r");
    FILE *f_out = fopen(output_filename, "wb");

    int i;
    int NINTCI;
    int NINTCF;
    int iBuffer;
    int iBuffer_LCC[6];
    double dBuffer[8];

    fscanf( f_in, "%d", &NINTCI )
    fwrite( &NINTCI, sizeof(int), 1, f_out );

    fscanf( f_in, "%d", &NINTCF);
    fwrite( &NINTCF, sizeof(int), 1, f_out );

    fscanf( f_in, "%d", &iBuffer);
    fwrite( &iBuffer, sizeof(int), 1, f_out );

    fscanf( f_in, "%d", &iBuffer);
    fwrite( &iBuffer, sizeof(int), 1, f_out );


        for (i = NINTCI; i <= NINTCF; i++){
                fscanf( f_in, "%d", &iBuffer_LCC[0] );
                fscanf( f_in, "%d", &iBuffer_LCC[1] );
                fscanf( f_in, "%d", &iBuffer_LCC[2] );
                fscanf( f_in, "%d", &iBuffer_LCC[3] );
                fscanf( f_in, "%d", &iBuffer_LCC[4] );
                fscanf( f_in, "%d", &iBuffer_LCC[5] );

                fwrite( &iBuffer_LCC, sizeof(int), sizeof(iBuffer_LCC) / sizeof(int), f_out );
        }

        for (i = NINTCI; i <= NINTCF; i++){
                fscanf( f_in, "%lf", &dBuffer[0] );
                fscanf( f_in, "%lf", &dBuffer[1] );
                fscanf( f_in, "%lf", &dBuffer[2] );
                fscanf( f_in, "%lf", &dBuffer[3] );
                fscanf( f_in, "%lf", &dBuffer[4] );
                fscanf( f_in, "%lf", &dBuffer[5] );
                fscanf( f_in, "%lf", &dBuffer[6] );
                fscanf( f_in, "%lf", &dBuffer[7] );

                fwrite( &dBuffer, sizeof(double), sizeof(dBuffer) / sizeof(double), f_out );
        }

        for (i = NINTCI; i <= NINTCF; i++){
                fscanf( f_in, "%d", &iBuffer);
                fwrite( &iBuffer, sizeof(int), 1, f_out );
        }

        fclose(f_in);
        fclose(f_out);
        printf("Binary file succesfully created\n");
}

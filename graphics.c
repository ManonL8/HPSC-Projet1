#include <stdio.h>
#include <stdlib.h>
#include <math.h>


/* Function parameters:
 - argv[0]: ./graphics
 - argv[1]: .dat
 - argv[2]: .txt

*/

// Function to find cross product of two vector array.
void crossProduct(float *vect_A, float *vect_B, float *cross_P)

{
    cross_P[0] = vect_A[1]*vect_B[2] - vect_A[2]*vect_B[1];
    cross_P[1] = vect_A[2]*vect_B[0] - vect_A[0]*vect_B[2];
    cross_P[2] = vect_A[0]*vect_B[1] - vect_A[1]*vect_B[0];
}

// Function to find scalar product of two vectors
float ScalarProduct(float *vect_A, float *vect_B) 
{ 
    float scalarProduct; 
    for (int i = 0; i < 3; i++) 
       scalarProduct += vect_A[i]*vect_B[i]; 
    return scalarProduct; 
}

// Function to find the norm of a vector
float norm(float *vect)
{
    float norm = 0;

    for (int i = 0 ; i < 3 ; i++)
    {
      norm += vect[i]*vect[i];
    }
    norm = sqrt(norm);
    return norm;
}

// Function to compute a vector from two points
float *vect(float *point_A, float *point_B, float *vect)
{
    vect[0] = point_B[0] - point_A[0];
    vect[1] = point_B[1] - point_A[1];
    vect[2] = point_B[2] - point_A[2];
}

/*----------------MAIN-----------------*/

int main(int argc, char **argv){
    if(argc<3)
        printf("\nEnter A Valid Number Of Arguments: 3");
    else if(argc == 3)
    {

        /*--------Store the triangle vertex positions in an array--------*/

        unsigned int nb_triangles = 0;

        char *model_file = argv[1];
        char *param_file = argv[2];
        FILE *model = fopen(model_file, "r");
        if(model == NULL)
        {
          printf("Error: can't open file to read\n");
          exit(0);
        }
        /* Seek to the beginning of the file */
        fseek(model, 0, SEEK_SET);

        /* Read and display data */
        int count = fread(&nb_triangles, sizeof(unsigned), 1, model); // Count the number of triangles (first data in the .data file)

        // Printing data to check validity
        printf("Data read from file: %d \n", nb_triangles);
        printf("Elements read: %d\n", count);

        float vertices[nb_triangles][9]; // Matrix containing the triangles vertices

        for (int i = 0; i < nb_triangles; i++) // Loop on the rows, 1 row of size 9 per triangle
        {
            for (int j = 0; j < 9; j++)
            {
                fread(&vertices[i][j], sizeof(float), 1, model);
                printf("%f, ", vertices[i][j]);
            }
            printf("\n");
        }
        fclose(model);

        /*------------Store the T shading parameters------------*/

        float lightBeam[3];
        lightBeam[0] = lightBeam[1] = lightBeam[2] = -1/sqrt(3);
        float vect_AB[3], vect_AC[3];
        float point_A[3], point_B[3], point_C[3];
        float normal[3];
        float s[nb_triangles];
        float vectorProduct[3];

        for (int i = 0; i < nb_triangles; i++)
        {
            point_A[0] = vertices[i][0];
            point_A[1] = vertices[i][1];
            point_A[2] = vertices[i][2];

            //printf("%lf, %lf, %lf\n", point_A[0], point_A[1], point_A[2]);

            point_B[0] = vertices[i][3];
            point_B[1] = vertices[i][4];
            point_B[2] = vertices[i][5];

            point_C[0] = vertices[i][6];
            point_C[1] = vertices[i][7];
            point_C[2] = vertices[i][8];

            vect(point_A, point_B, vect_AB);
            vect(point_A, point_B, vect_AC);

            crossProduct(vect_AB, vect_AC, vectorProduct);
            normal[0] = vectorProduct[0]/(norm(vect_AB)*norm(vect_AC));
            normal[1] = vectorProduct[1]/(norm(vect_AB)*norm(vect_AC));
            normal[2] = vectorProduct[2]/(norm(vect_AB)*norm(vect_AC));

            s[i] = ScalarProduct(normal, lightBeam); 

            if(s[i] <= 0)
            {
                s[i] = 0;
            }
            printf("%f\n", s[i]);

        }

    }
    return 0;
}

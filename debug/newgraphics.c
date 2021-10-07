#include <stdio.h>
#include <stdlib.h>
#include <math.h>


typedef struct Pixel{
  unsigned char r;
  unsigned char g;
  unsigned char b;
  float d;
}PIX;

// Function to find cross product of two vector array.
void crossProduct(float *vect_A, float *vect_B, float *cross_P);
void crossProduct(float *vect_A, float *vect_B, float *cross_P){
    cross_P[0] = (vect_A[1]*vect_B[2]) - (vect_A[2]*vect_B[1]);
    cross_P[1] = (vect_A[2]*vect_B[0]) - (vect_A[0])*(vect_B[2]);
    cross_P[2] = (vect_A[0]*vect_B[1]) - (vect_A[1])*(vect_B[0]);
}

// Function to find scalar product of two vectors
float scalarProduct(float *vect_A, float *vect_B);
float scalarProduct(float *vect_A, float *vect_B){
    float scalarProduct = 0.0;
    for(size_t i = 0; i < 3; i++){
      scalarProduct += (vect_A[i])*(vect_B[i]);
    }

    return scalarProduct;
}

// Function to find the norm of a vector
float norm(float *vect);
float norm(float *vect){
    float norm = 0.0;
    for(size_t i = 0; i < 3; i++){
      norm += (vect[i])*(vect[i]);
    }
    norm = sqrt(norm);

    return norm;
}

// Function to compute a vector from two points
void vect(float* point_A, float* point_B, float *vect_AB);
void vect(float* point_A, float* point_B, float *vect_AB){
    vect_AB[0] = point_B[0] - point_A[0];
    vect_AB[1] = point_B[1] - point_A[1];
    vect_AB[2] = point_B[2] - point_A[2];
}
// Function to compute the product between a matrix and a vector
void matriceMult(size_t m, size_t n, float **mat1, float *mat2, float *mult);
void matriceMult(size_t m, size_t n, float **mat1, float *mat2, float *mult){
    for(size_t i = 0; i < m; i++){
        mult[i] = 0;
        for(size_t j = 0; j < n; j++){
            mult[i] += mat1[i][j] * mat2[j];
        }
    }
}

void changeOfCoordinates(float **view_mat, float **projection_mat, float* vertex);
void changeOfCoordinates(float **view_mat, float **projection_mat, float* vertex){
  float result[3], coordinates[4];
  coordinates[0] = vertex[0];
  coordinates[1] = vertex[1];
  coordinates[2] = vertex[2];
  coordinates[3] = 1;

  matriceMult(3, 4, view_mat, coordinates, result);
  vertex[0] = coordinates[0] = result[0];
  vertex[1] = coordinates[1] = result[1];
  vertex[2] = coordinates[2] = result[2];

  projection_mat[0][0] = projection_mat[0][0]/(vertex[2]);
  projection_mat[1][1] = projection_mat[1][1]/(vertex[2]);
  projection_mat[2][2] = projection_mat[2][2]/(vertex[2]);
  projection_mat[2][3] = projection_mat[2][3]/(vertex[2]);

  matriceMult(3, 4, projection_mat, coordinates, result);
  vertex[0] = result[0];
  vertex[1] = result[1];
  vertex[2] = result[2];
}

void freeMatrix(PIX **matrix, unsigned n);
void freeMatrix(PIX **matrix, unsigned n){
  for(unsigned i = 0; i < n; i++){
    free(matrix[i]);
  }
  free(matrix);
}

int max(int a,int b);
int max(int a,int b){
  if(a >= b)
    return a;
  else
    return b;
}

int min(int a,int b);
int min(int a,int b){
  if(a <= b)
    return a;
  else
    return b;
}
/*----------------MAIN-----------------*/
/* Function parameters:
 - argv[0]: ./graphics
 - argv[1]: .dat
 - argv[2]: .txt
*/

int main(int argc, char **argv){
    if(argc<3)
        printf("\nEnter A Valid Number Of Arguments: 3");
    else if(argc == 3){

        /*--------Fetch the elements from the model format .dat file--------*/

        unsigned nb_triangles = 0;

        char *model = argv[1];
        char *param = argv[2];
        FILE *model_file = fopen(model, "rb");
        if(model_file == NULL){
          printf("Error: can't open file to read\n");
          exit(0);
        }
        /* Seek to the beginning of the file */
        fseek(model_file, 0, SEEK_SET);

        /* Read and display data */
        unsigned int count = fread(&nb_triangles, sizeof(unsigned int), 1, model_file); // Fetch the number of triangles (first data in the .data file)

        // Printing data to check validity
        printf("Data read from file: %d \n", nb_triangles);
        printf("Elements read: %d\n", count);

        // array
        float *vertices = malloc(9*nb_triangles*sizeof(float));
        if(!vertices){
          printf("Error: vertices is a null pointer\n");
          exit(0);
        }

        //fseek(model_file, 1, SEEK_SET);

        count = 0;
        count = fread(vertices, sizeof(float), 9*nb_triangles, model_file);

        for(size_t i = 0; i < 9*nb_triangles; i++){
          printf("%f\n",vertices[i]);
        }

        if(count != 9*nb_triangles){
          printf("Error: fread() function has not read all the elements of the file\n");
          exit(0);
        }

        fclose(model_file);


        /*--------Fetch the scene parameters from the .txt file--------*/

        FILE *param_file = fopen(param, "r");

        unsigned char object_color[3];
        unsigned char background_color[3];
        float camera_coordinates[3];
        unsigned int hw[2];
        float theta_y;
        float nf[2];

        if(!param_file){
          printf("Error: can't open the param file to read\n");
          exit(0);
        }
        /* Seek to the beginning of the file */
        /*
        fseek(param_file, 0, SEEK_SET);

        /* Read and display data */

        for(size_t i = 0; i < 3; i++){
            fscanf(param_file, "%hhd", &object_color[i]);
        }

        for(size_t i = 0; i < 3; i++){
            fscanf(param_file, "%hhd", &background_color[i]);
        }

        for(size_t i = 0; i < 3; i++){
            fscanf(param_file, "%f", &camera_coordinates[i]);
        }

        for(size_t i = 0; i < 2; i++){
            fscanf(param_file, "%d", &hw[i]);
        }

        fscanf(param_file, "%f", &theta_y);

        for(size_t i = 0; i < 2; i++){
            fscanf(param_file, "%f", &nf[i]);
        }

        fclose(param_file);

        /*------------Store the T shading parameters------------*/

        float lightBeam[3];
        float s[nb_triangles];

        lightBeam[0] = -1.0/sqrt(3);
        lightBeam[1] = -1.0/sqrt(3);
        lightBeam[2] = -1.0/sqrt(3);

        /*-------Change of coordinates matrices (4)-------*/

        float vect_OM[3];

        POINT point_M;
        vect_OM[0] = point_M[0] = camera_coordinates[0];
        vect_OM[1] = point_M[1] = camera_coordinates[1];
        vect_OM[2] = point_M[2] = camera_coordinates[2];

        float d = sqrt(point_M[0]*point_M[0] + point_M[1]*point_M[1] + point_M[2]*point_M[2]);
        float d_g = sqrt(point_M[0]*point_M[0] + point_M[2]*point_M[2]);

        float e_u[3];

        float e_v[3]

        float e_w[3]
        float view_mat[3][4];

        e_u[0] = view_mat[0][0] = (point_M[2])/d_g;
        e_u[1] = view_mat[0][1] = 0;
        e_u[2] = view_mat[0][2] = -(point_M[0])/d_g;

        e_v[0] = view_mat[1][0] = -(point_M[0] * point_M[1])/(d_g*d);
        e_v[1] = view_mat[1][1] = (point_M[2] * point_M[2] + point_M[0] * point_M[0])/(d_g*d);
        e_v[2] = view_mat[1][2] = -(point_M[2] * point_M[1])/(d_g*d);

        e_w[0] = view_mat[2][0] = (point_M[0])/d;
        e_w[1] = view_mat[2][1] = (point_M[1])/d;
        e_w[2] = view_mat[2][2] = (point_M[2])/d;

        view_mat[0][3] = -scalarProduct(vect_OM, e_u);
        view_mat[1][3] = -scalarProduct(vect_OM, e_v);
        view_mat[2][3] = -scalarProduct(vect_OM, e_w);

        /*-------Change of coordinates matrices (5)-------*/

        float r = hw[1]/hw[0];
        float S_y = 1/tan(theta_y/2);
        float S_x = S_y/r;
        float projection_mat[3][4] = malloc(3*sizeof(float*));


        float k = (nf[0] - nf[1])/(2*nf[0]*nf[1]);
        for(size_t i = 0; i < 3; i++){
            for(size_t j = 0; j <4; j++){
              if (i == 0 && j == 0)
                projection_mat[i][j] = k*S_x;
              else if (i == 1 && j == 1)
                projection_mat[i][j] = k*S_y;
              else if (i == 2 && j == 2)
                projection_mat[i][j] = k*(nf[0] + nf[1])/(nf[0] - nf[1]);
              else if (i == 2 && j == 3)
                projection_mat[i][j] = k*(-1);
              else
                projection_mat[i][j] = 0;
            }
        }


        float vectorProduct[3];
        //printf("AVANT LA BOUCLE\n");
        for(size_t i = 0; i < nb_triangles; i++){
          /*-------Shading parameter and normal-------*/

          vect(vertices[i][0], vertices[i][1], vect_AB);
          vect(vertices[i][0], vertices[i][2], vect_AC);

          crossProduct(vect_AB, vect_AC, vectorProduct);

          norm_vectorProduct = norm(vectorProduct);
          normal[0] = (vectorProduct[0])/norm_vectorProduct);
          normal[1] = (vectorProduct[1])/norm_vectorProduct);
          normal[2] = (vectorProduct[2])/norm_vectorProduct);
          s[i] = scalarProduct(normal, lightBeam);

          if(s[i] < 0)
              s[i] = -s[i];
          else
            s[i] = 0;
          //printf("%f\n", s[i]);

          /*-------Change of coordinates formula-------*/
          changeOfCoordinates(view_mat, projection_mat, vertices[i + 0*nb_triangles]);
          changeOfCoordinates(view_mat, projection_mat, vertices[i + 1*nb_triangles]);
          changeOfCoordinates(view_mat, projection_mat, vertices[i + 2*nb_triangles]);
        }

        /*--------------------------RASTERIZATION-----------------------------*/

        PIX **p = malloc(hw[0]*sizeof(PIX *));
        if(!p){
          printf("Error: allocation error for the p array");
          exit(0);
        }
        for(size_t i = 0; i < hw[0]; i++){
          p[i] = malloc(hw[1]*sizeof(PIX));
          if(!p[i]){
            printf("Error: allocation error for the p array");
            exit(0);
          }
        }

        for(size_t i = 0; i < hw[0]; i++){
          for(size_t j = 0; j < hw[1]; j++){
            p[i][j].r = background_color[0];
            p[i][j].g = background_color[1];
            p[i][j].b = background_color[2];
            p[i][j].d = 1;
          }
        }

        // First we loop over the nb_triangles
        // Ici faire un vecteur temporaire qui stocke les points du current triangle


        float curr_A[3], curr_B[3], curr_C[3], vec_AB[3], vec_AC[3], vec_BC[3], vec_BA[3], vec_P[3];
        float scal_prod1, scal_prod2, w_A, w_B, w_C, d_P;
        float max_x, min_x, max_y, min_y;

        for(size_t k = 0; k < nb_triangles; k++){
          for(size_t l = 0; l < 3; l++){
            curr_A[l] = vertices[l + 9*k];
            curr_B[l] = vertices[3 + l + 9*k];
            curr_C[l] = vertices[6 + l + 9*k];
            printf("vertices abs: %f\n", vertices[l + 9*k]);
            printf("vertices ord: %f\n", vertices[l + 9*k]);
            printf("vertices alt: %f\n", vertices[l + 9*k]);
          }

          vect(curr_A, curr_B, vec_AB);
          vect(curr_A, curr_C, vec_AC);
          vect(curr_B, curr_C, vec_BC);
          vect(curr_B, curr_A, vec_BA);
          // On parcourt les pixels et pour chaque pixel on vérifie s'il est dans
          // le current triangle, optimisation ?

          // Ici l'idée ça serait d'avoir une estimation du contour du triangle
          // et de regarder seulement les pixels qui seraient à l'intérieur de
          // cette estimation et non pas de regarder chacun des pixels

          min_x = min(curr_A[0], curr_B[0]);
          min_x = min(min_x, curr_C[0]);

          max_x = max(curr_A[0], curr_B[0]);
          max_x = max(max_x, curr_C[0]);

          min_y = min(curr_A[1], curr_B[1]);
          min_y = min(min_x, curr_C[1]);

          max_y = max(curr_A[1], curr_B[1]);
          max_y = max(max_y, curr_C[1]);

          printf("minx : %f\n", min_x);
          printf("maxx : %f\n", max_x);
          printf("miny : %f\n", min_y);
          printf("maxy : %f\n", max_y);
          size_t j_min_x = ((min_x + 1)*hw[1])/2 - 0.5;
          size_t j_max_x = ((max_x + 1)*hw[1])/2 - 0.5;
          size_t i_min_y = ((min_y - 1)*hw[0])/(-2) - 0.5;
          size_t i_max_y = ((max_y - 1)*hw[0])/(-2) - 0.5;

          float pix_point[3];

          for(size_t i = i_min_y; i <= i_max_y ; i++){
            for(size_t j = j_min_x; j <= j_max_x; j++){
              // changer les coordonnées en xs, ys
              pix_point[0] = 2*(j + 0.5)/hw[1] - 1;
              pix_point[1] = -2*(i + 0.5)/hw[0] + 1;
              pix_point[2] = 0;
              vect(curr_A, pix_point, vec_P);
              crossProduct(vec_AB, vec_P, vec_AB);
              crossProduct(vec_P, vec_AC, vec_AC);
              scal_prod1 = scalarProduct(vec_AB, vec_AC);

              vect(curr_B, pix_point, vec_P);
              crossProduct(vec_BC, vec_P, vec_BC);
              crossProduct(vec_P, vec_BA, vec_BA);
              scal_prod2 = scalarProduct(vec_BC, vec_BA);

              if(scal_prod1 >= 0 && scal_prod2 >= 0){
                w_A = ((curr_B[0] - curr_A[0])*(pix_point[1] - curr_A[1])
                - (curr_B[1] - curr_A[1])*(pix_point[0] - curr_A[0]))
                /((curr_B[0] - curr_A[0])*(curr_B[1] - curr_A[1]) -
                (curr_B[1] - curr_A[1])*(curr_C[0] - curr_A[0]));
                //wA = ((xB − xA)(yP − yA) − (yB − yA)(xP − xA))/((xB − xA)(yC − yA) − (yB − yA)(xC − xA));

                w_B = ((curr_A[0] - curr_C[0])*(pix_point[1] - curr_C[1]) -
                (curr_A[1] - curr_C[1])*(pix_point[0] - curr_C[0]))
                /((curr_A[0] - curr_C[0])*(curr_B[1] - curr_C[1]) -
                (curr_A[1] - curr_B[1])*(curr_B[0] - curr_C[0]));
                //wB = ((xA − xC)(yP − yC) − (yA − yC)(xP − xC))/((xA − xC)(yB − yC) − (yA − yC)(xB − xC));

                w_C = 1 - w_A - w_B;

                d_P = w_A * curr_A[2] + w_B * curr_B[2] + w_C * curr_C[2];

                if(d_P  < p[i][j].d){
                  p[i][j].d = d_P;
                  p[i][j].r = s[k]*object_color[0];
                  p[i][j].g = s[k]*object_color[1];
                  p[i][j].b = s[k]*object_color[2];
                }
              }
            }
          }
        }

        /*-------------------Writing the results in .ppm file----------------*/
       FILE *Image_2D = fopen("Image2D.ppm", "wb"); /* b - binary mode */
       fprintf(Image_2D, "P6\n%d %d\n255\n", hw[1], hw[0]);

        for (size_t i = 0; i < hw[0]; ++i){
          for (size_t j = 0; j < hw[1]; ++j){
            static unsigned char color[3];
            color[0] = p[i][j].r;
            color[1] = p[i][j].g;
            color[2] = p[i][j].b;
            fwrite(color, 1, 3, Image_2D);
          }
        }


        fclose(Image_2D);

        // Ici je me disais si le fonction freeMatrix pouvait pas prendre void **matrix
        // en argument ça serait plus simple
        // Freeing the matrix of pixels
        for(size_t i = 0; i < hw[0]; i++){
          free(p[i]);
        }
        free(p);
    }

    return 0;
}

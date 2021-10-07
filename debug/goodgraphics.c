#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
The structure point is defined by :
 - its coordinates in the world frame
 - its coordinates in the view frame
*/

typedef struct Point{
  float abs;
  float ord;
  float alt;
}POINT;
okok
 /*
 The structure vector is defined by :
  - its coordinates in the world frame
  - its coordinates in the view frame
  - its norm
*/

typedef struct Vector{
  float abs;
  float ord;
  float alt;
  float norm;
}VEC;

typedef struct Pixel{
  unsigned char r;
  unsigned char g;
  unsigned char b;
  float d;
}PIX;

// Function to find cross product of two vector array.
void crossProduct(VEC *vect_A, VEC *vect_B, VEC *cross_P);
void crossProduct(VEC *vect_A, VEC *vect_B, VEC *cross_P){
    cross_P->abs = (vect_A->abs)*(vect_B->alt) - (vect_A->alt)*(vect_B->ord);
    cross_P->ord = (vect_A->alt)*(vect_B->abs) - (vect_A->abs)*(vect_B->alt);
    cross_P->alt = (vect_A->abs)*(vect_B->ord) - (vect_A->ord)*(vect_B->abs);
}

// Function to find scalar product of two vectors
float scalarProduct(VEC *vect_A, VEC *vect_B);
float scalarProduct(VEC *vect_A, VEC *vect_B){
    float scalarProduct = 0.0;

    scalarProduct += (vect_A->abs)*(vect_B->abs);
    scalarProduct += (vect_A->ord)*(vect_B->ord);
    scalarProduct += (vect_A->alt)*(vect_B->alt);

    return scalarProduct;
}

// Function to find the norm of a vector
float norm(VEC *vect);
float norm(VEC *vect){
    float norm = 0.0;

    norm += (vect->abs)*(vect->abs);
    norm += (vect->ord)*(vect->ord);
    norm += (vect->alt)*(vect->alt);
    norm = sqrt(norm);

    return norm;
}

// Function to compute a vector from two points
void vect(POINT point_A, POINT point_B, VEC *vect_AB);
void vect(POINT point_A, POINT point_B, VEC *vect_AB){
    vect_AB->abs = point_B.abs - point_A.abs;
    vect_AB->ord = point_B.ord - point_A.ord;
    vect_AB->alt = point_B.alt - point_A.alt;
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

void changeOfCoordinates(float **view_mat, float **projection_mat, POINT vertex);
void changeOfCoordinates(float **view_mat, float **projection_mat, POINT vertex){
  float result[3], coordinates[4];
  coordinates[0] = vertex.abs;
  coordinates[1] = vertex.ord;
  coordinates[2] = vertex.alt;
  coordinates[3] = 1;

  matriceMult(3, 4, view_mat, coordinates, result);
  vertex.abs = coordinates[0] = result[0];
  vertex.ord = coordinates[1] = result[1];
  vertex.alt = coordinates[2] = result[2];

  projection_mat[0][0] = projection_mat[0][0]/(vertex.alt);
  projection_mat[1][1] = projection_mat[1][1]/(vertex.alt);
  projection_mat[2][2] = projection_mat[2][2]/(vertex.alt);
  projection_mat[2][3] = projection_mat[2][3]/(vertex.alt);

  matriceMult(3, 4, projection_mat, coordinates, result);
  vertex.abs = result[0];
  vertex.ord = result[1];
  vertex.alt = result[2];
}

void freeMatrix(float **matrix, unsigned n);
void freeMatrix(float **matrix, unsigned n){
  for(unsigned i = 0; i < n; i++){
    free(matrix[i]);
  }
  free(matrix);
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
        //fseek(model_file, 0, SEEK_SET);

        /* Read and display data */
        unsigned int count = fread(&nb_triangles, sizeof(unsigned int), 1, model_file); // Fetch the number of triangles (first data in the .data file)

        // Printing data to check validity
        printf("Data read from file: %d \n", nb_triangles);
        printf("Elements read: %d\n", count);

        // Matrix whose elements are POINT structures containing the triangles
        // vertices
        POINT **vertices = malloc(nb_triangles*sizeof(POINT*));
        if(vertices == NULL){
          printf("Error: vertices is a null pointer\n");
          exit(0);
        }

        for(size_t i = 0; i < nb_triangles; i++){
          vertices[i] = malloc(3*sizeof(POINT));
          if(vertices[i] == NULL){
            printf("Error: vertices is a null pointer\n");
            exit(0);
          }
        }

        //fseek(model_file, 1, SEEK_SET);

        count = 0;
        for(size_t i = 0; i < nb_triangles; i++){
          for(size_t j = 0; j < 3; j++){
            count += fread(&(vertices[i][j].abs), sizeof(float), 1, model_file);
            count += fread(&(vertices[i][j].ord), sizeof(float), 1, model_file);
            count += fread(&(vertices[i][j].alt), sizeof(float), 1, model_file);
            printf("vertices abs: %f\n", vertices[i][j].abs);
            printf("vertices ord: %f\n", vertices[i][j].ord);
            printf("vertices alt: %f\n", vertices[i][j].alt);

          }
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

        if(param_file == NULL){
          printf("Error: can't open file to read\n");
          exit(0);
        }
        /* Seek to the beginning of the file */
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

        VEC *lightBeam = malloc(sizeof(VEC));
        if(lightBeam == NULL){
          printf("Error: Problem of allocation");
          exit(0);
        }

        lightBeam->abs = -1.0/sqrt(3);
        lightBeam->ord = -1.0/sqrt(3);
        lightBeam->alt = -1.0/sqrt(3);

        VEC *vect_AB = malloc(sizeof(VEC));
        if(vect_AB == NULL){
          printf("Error: Problem of allocation");
          exit(0);
        }

        VEC *vect_AC = malloc(sizeof(VEC));
        if(vect_AC == NULL){
          printf("Error: Problem of allocation");
          exit(0);
        }

        VEC *normal = malloc(sizeof(VEC));
        if(vect_AC == NULL){
          printf("Error: Problem of allocation");
          exit(0);
        }

        VEC *vectorProduct = malloc(sizeof(VEC));
        if(vectorProduct == NULL){
          printf("Error: Problem of allocation");
          exit(0);
        }

        float s[nb_triangles];

        /*-------Change of coordinates matrices (4)-------*/

        VEC *vect_OM = malloc(sizeof(VEC));
        if(vect_OM == NULL){
          printf("Error: Problem of allocation");
          exit(0);
        }

        POINT point_M;
        vect_OM->abs = point_M.abs = camera_coordinates[0];
        vect_OM->ord = point_M.ord = camera_coordinates[1];
        vect_OM->alt = point_M.alt = camera_coordinates[2];

        float d = sqrt(point_M.abs*point_M.abs + point_M.ord*point_M.ord + point_M.alt*point_M.alt);
        float d_g = sqrt(point_M.abs*point_M.abs + point_M.alt*point_M.alt);

        VEC *e_u = malloc(sizeof(VEC));
        if(e_u == NULL){
          printf("Error: Problem of allocation");
          exit(0);
        }

        VEC *e_v = malloc(sizeof(VEC));
        if(e_v == NULL){
          printf("Error: Problem of allocation");
          exit(0);
        }

        VEC *e_w = malloc(sizeof(VEC));
        if(e_w == NULL){
          printf("Error: Problem of allocation");
          exit(0);
        }
        float **view_mat = malloc(3*sizeof(float*));
        if(view_mat == NULL){
          printf("Error: mat1 is a null pointer\n");
          exit(0);
        }

        for(size_t i = 0; i < 3; i++){
          view_mat[i] = malloc(4*sizeof(float));
          if(view_mat[i] == NULL){
            printf("Error: mat1 is a null pointer\n");
            exit(0);
          }
        }

        e_u->abs = view_mat[0][0] = (point_M.alt)/d_g;
        e_u->ord = view_mat[0][1] = 0;
        e_u->alt = view_mat[0][2] = -(point_M.abs)/d_g;

        e_v->abs = view_mat[1][0] = -(point_M.abs * point_M.ord)/(d_g*d);
        e_v->ord = view_mat[1][1] = (point_M.alt * point_M.alt + point_M.abs * point_M.abs)/(d_g*d);
        e_v->alt = view_mat[1][2] = -(point_M.alt * point_M.ord)/(d_g*d);

        e_w->abs = view_mat[2][0] = (point_M.abs)/d;
        e_w->ord = view_mat[2][1] = (point_M.ord)/d;
        e_w->alt = view_mat[2][2] = (point_M.alt)/d;

        view_mat[0][3] = -scalarProduct(vect_OM, e_u);
        view_mat[1][3] = -scalarProduct(vect_OM, e_v);
        view_mat[2][3] = -scalarProduct(vect_OM, e_w);

        /*-------Change of coordinates matrices (5)-------*/

        float r = hw[1]/hw[0];
        float S_y = 1/tan(theta_y/2);
        float S_x = S_y/r;
        float **projection_mat = malloc(3*sizeof(float*));
        if(projection_mat == NULL){
          printf("Error: projection_mat is a null pointer\n");
          exit(0);
        }

        for(size_t i = 0; i < 3; i++){
          projection_mat[i] = malloc(4*sizeof(float));
          if(projection_mat[i] == NULL){
            printf("Error: projection_mat[%ld] is a null pointer\n", i);
            exit(0);
          }
        }
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


        //printf("AVANT LA BOUCLE\n");
        for(size_t i = 0; i < nb_triangles; i++){
          /*-------Shading parameter and normal-------*/

          vect(vertices[i][0], vertices[i][1], vect_AB);
          vect(vertices[i][0], vertices[i][2], vect_AC);

          crossProduct(vect_AB, vect_AC, vectorProduct);

          vectorProduct->norm = norm(vectorProduct);
          normal->abs = (vectorProduct->abs)/(vectorProduct->norm);
          normal->ord = (vectorProduct->ord)/(vectorProduct->norm);
          normal->alt = (vectorProduct->alt)/(vectorProduct->norm);
          s[i] = scalarProduct(normal, lightBeam);

          if(s[i] < 0)
              s[i] = -s[i];
          else
            s[i] = 0;
          printf("%f\n", s[i]);

          /*-------Change of coordinates formula-------*/
          changeOfCoordinates(view_mat, projection_mat, vertices[i][0]);
          changeOfCoordinates(view_mat, projection_mat, vertices[i][1]);
          changeOfCoordinates(view_mat, projection_mat, vertices[i][2]);
        }
        free(lightBeam);
        free(vect_AB);
        free(vect_AC);
        free(normal);
        free(vectorProduct);
        free(vect_OM);
        free(e_u);
        free(e_v);
        free(e_w);
        freeMatrix(view_mat, 3);
        freeMatrix(projection_mat, 3);

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

        POINT *curr = malloc(3*sizeof(POINT));
        if(!curr){
          printf("Error: allocation error for curr_triangle");
          exit(0);
        }

        VEC *vec_AB = malloc(sizeof(VEC));
        if(vec_AB == NULL){
          printf("Error: Problem of allocation");
          exit(0);
        }
        VEC *vec_AC = malloc(sizeof(VEC));
        if(vec_AB == NULL){
          printf("Error: Problem of allocation");
          exit(0);
        }
        VEC *vec_BC = malloc(sizeof(VEC));
        if(vec_BC == NULL){
          printf("Error: Problem of allocation");
          exit(0);
        }
        VEC *vec_BA = malloc(sizeof(VEC));
        if(vec_BA == NULL){
          printf("Error: Problem of allocation");
          exit(0);
        }
        VEC *vec_P = malloc(sizeof(VEC));
        if(vec_P == NULL){
          printf("Error: Problem of allocation");
          exit(0);
        }

        POINT pix_point;
        float scal_prod1, scal_prod2, w_A, w_B, w_C, d_P;
        float max_x, min_x, max_y, min_y;

        for(size_t k = 0; k < nb_triangles; k++){
          for(size_t l = 0; l < 3; l++){
            curr[l].abs = vertices[k][l].abs;
            curr[l].ord = vertices[k][l].ord;
            curr[l].alt = vertices[k][l].alt;
            printf("vertices abs: %f\n", vertices[k][l].abs);
            printf("vertices ord: %f\n", vertices[k][l].ord);
            printf("vertices alt: %f\n", vertices[k][l].alt);
          }

          vect(curr[0], curr[1], vec_AB);
          vect(curr[0], curr[2], vec_AC);
          vect(curr[1], curr[2], vec_BC);
          vect(curr[1], curr[0], vec_BA);
          // On parcourt les pixels et pour chaque pixel on vérifie s'il est dans
          // le current triangle, optimisation ?

          // Ici l'idée ça serait d'avoir une estimation du contour du triangle
          // et de regarder seulement les pixels qui seraient à l'intérieur de
          // cette estimation et non pas de regarder chacun des pixels

          min_x = min_y = 10000000000;
          max_x = max_y = -1000000000;
          for(size_t m = 0; m < 3; m++){
            if(curr[m].abs < min_x)
            {
              min_x = curr[m].abs;
              printf("%f\n", curr[m].abs);
            }

            if(curr[m].abs > max_x)
              max_x = curr[m].abs;
            if(curr[m].ord < min_y)
              min_y = curr[m].ord;
            if(curr[m].ord > max_y)
              max_y = curr[m].ord;
          }

          printf("minx : %f\n", min_x);
          printf("maxx : %f\n", max_x);
          printf("miny : %f\n", min_y);
          printf("maxy : %f\n", max_y);
          size_t j_min_x = ((min_x + 1)*hw[1])/2 - 0.5;
          size_t j_max_x = ((max_x + 1)*hw[1])/2 - 0.5;
          size_t i_min_y = ((min_y - 1)*hw[0])/(-2) - 0.5;
          size_t i_max_y = ((max_y - 1)*hw[0])/(-2) - 0.5;

          for(size_t i = i_min_y; i <= i_max_y ; i++){
            for(size_t j = j_min_x; j <= j_max_x; j++){
              // changer les coordonnées en xs, ys
              pix_point.abs = 2*(j + 0.5)/hw[1] - 1;
              pix_point.ord = -2*(i + 0.5)/hw[0] + 1;
              pix_point.alt = 0;
              vect(curr[0], pix_point, vec_P);
              crossProduct(vec_AB, vec_P, vec_AB);
              crossProduct(vec_P, vec_AC, vec_AC);
              scal_prod1 = scalarProduct(vec_AB, vec_AC);

              vect(curr[1], pix_point, vec_P);
              crossProduct(vec_BC, vec_P, vec_BC);
              crossProduct(vec_P, vec_BA, vec_BA);
              scal_prod2 = scalarProduct(vec_BC, vec_BA);

              if(scal_prod1 >= 0 && scal_prod2 >= 0){
                w_A = ((curr[1].abs - curr[0].abs)*(pix_point.ord - curr[0].ord)
                - (curr[1].ord - curr[0].ord)*(pix_point.abs - curr[0].abs))
                /((curr[1].abs - curr[0].abs)*(curr[2].ord - curr[0].ord) -
                (curr[1].ord - curr[0].ord)*(curr[2].abs - curr[0].abs));
                //wA = ((xB − xA)(yP − yA) − (yB − yA)(xP − xA))/((xB − xA)(yC − yA) − (yB − yA)(xC − xA));

                w_B = ((curr[0].abs - curr[2].abs)*(pix_point.ord - curr[2].ord) -
                (curr[0].ord - curr[2].ord)*(pix_point.abs - curr[2].abs))
                /((curr[0].abs - curr[2].abs)*(curr[1].ord - curr[2].ord) -
                (curr[0].ord - curr[2].ord)*(curr[1].abs - curr[2].abs));
                //wB = ((xA − xC)(yP − yC) − (yA − yC)(xP − xC))/((xA − xC)(yB − yC) − (yA − yC)(xB − xC));

                w_C = 1 - w_A - w_B;

                d_P = w_A * curr[0].alt + w_B * curr[1].alt + w_C * curr[2].alt;

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

        for (size_t i = 0; i < hw[0]; i++){
          for (size_t j = 0; j < hw[1]; j++){
            static unsigned char color[3];
            color[0] = p[i][j].r;
            color[1] = p[i][j].g;
            color[2] = p[i][j].b;
            fwrite(color, 1, 3, Image_2D);
          }
        }


        fclose(Image_2D);

        for(size_t i = 0; i < nb_triangles; i++){
          free(vertices[i]);
        }
        free(vertices);
        free(curr);
        free(vec_AB);
        free(vec_AC);
        free(vec_BC);
        free(vec_BA);
        free(vec_P);
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

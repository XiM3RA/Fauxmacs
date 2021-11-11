#include <stdio.h>
#include <stdlib.h>

int main()
{
    int N;
    double delta_t, t_end;
    char *filename = "run-params.txt";
    FILE *fp = fopen(filename, "r");

    if (fp == NULL)
    {
        printf("Error: could not open file %s", filename);
        return 1;
    }

    fscanf (fp, "%d", &N);
    printf("%d\n", N);
    fscanf (fp, "%lf", &delta_t);
    printf ("%lf\n", delta_t);

    fclose(fp);

    printf("%d %lf %lf\n", N, delta_t, t_end);
    return 0;
}
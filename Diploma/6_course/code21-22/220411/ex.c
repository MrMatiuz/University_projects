#include <stdlib.h>
#include <stdio.h>

#define min(a,b)  (((a) < (b)) ? (a) : (b))

// int main(void)
// {
// 	double a[] = {9, 9, 9, 0, 0, 0, 9, 9, 0, 0, 9, 9};
// 	double b[] = {5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5};
// 	double c[sizeof(a) / sizeof(a[0])];
// 	int *intersection = NULL;
// 	int j=0;
// 	int needed_intersection = 1;

// 	for(int i=0; i < sizeof(a) / sizeof(a[0]); i++)
// 	{
// 		c[i] = (a[i] < b[i]);
// 		printf("%f ", a[i]);
// 		printf("%f ", b[i]);
// 		printf("%f\n", c[i]);
// 	}

// 	for(int i = 0; i < sizeof(a) / sizeof(a[0]) - 1; i++)
// 	{
// 		if (c[i+1] - c[i] == 1) 
// 		{
// 			intersection = (int*)realloc(intersection, (j+1)*sizeof(int));
// 			intersection[j] = i;
// 			j++;
// 		}
// 	}
// 	needed_intersection = intersection[min(needed_intersection, sizeof(intersection) / sizeof(intersection[0])) - 1];


// 	for(int i = 0; i < sizeof(intersection) / sizeof(intersection[0]); i++)
// 		printf("%d ", intersection[i]);
// 	printf("\n");

// 	printf("%d ", needed_intersection);

// 	return 0;
// }


double* get_inttersection(double a[], double b[], int arr_len, int needed_intersection_num)
{
	double c[arr_len];
	double * res = malloc(sizeof(double)*arr_len);
	int *intersection = NULL;
	int j=0;

	for(int i=0; i < arr_len; i++)
	{
		c[i] = (a[i] < b[i]);
		printf("%f ", a[i]);
		printf("%f ", b[i]);
		printf("%f\n", c[i]);
	}

	for(int i = 0; i < arr_len - 1; i++)
	{
		if (c[i+1] - c[i] == 1) 
		{
			intersection = (int*)realloc(intersection, (j+1)*sizeof(int));
			intersection[j] = i;
			j++;
		}
	}

	if (j != 0)
	{
		needed_intersection_num = min(needed_intersection_num, j);
		printf("%d %d %d", sizeof(intersection), sizeof(intersection[0]), j);

		for(int i = 0; i < j; i++)
			printf("%c ", intersection[i]);
		printf("\n");

		for(int i = 0; i < arr_len; i++)
			res[i] = (i > intersection[needed_intersection_num - 1]) ? a[i] : b[i];

		// for(int i = 0; i < arr_len; i++)
		// 	printf("%f\n", res[i]);
	}
	else
		res = a;

	return res;
}

int main(void)
{

	// double a[] = {9, 9, 9, 0, 0, 0, 9, 9, 0, 0, 9, 9};
	// double b[] = {5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5};

	double a[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	double b[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

	int needed_intersection_num = 1;
	int arr_len = sizeof(a) / sizeof(a[0]);
	double *r;

	r = get_inttersection(a, b, arr_len, needed_intersection_num);

	int l = 12;
	for(int i = 0; i < l; i++)
		printf("%f\n", r[i]);
	return 0;


}

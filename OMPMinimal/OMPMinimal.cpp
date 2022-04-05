#include <iostream>
#include "omp.h"


//Шаги
#define hStart 0.01
#define hEnd 0.000025

//Область
#define top 5
#define right 5
#define left -5
#define bottom -5


using namespace std;


double func(double x, double y) {
	return 9.2 * x - 0.6 * y + exp(0.81 * x * x + 0.19 * y * y);
	//9.2*x-0.6*y+exp(0.81*x^2+0.19*y^2)
}

int main()
{
	int threads = 4;
	double start = omp_get_wtime();
	double** minimals = new double*[threads];
	for (int i = 0; i < threads; i++)minimals[i] = new double[3];
	int xIterations = (right - left) / hStart;
	bool flag = false;

#pragma omp parallel for num_threads(threads) firstprivate(flag)
		for (int i = 0; i < xIterations; i++)
		{
			double x = left + i * hStart;
			for (double j = bottom; j < top; j += hStart)
			{
				if (!flag) {
					minimals[omp_get_thread_num()][1] = x;
					minimals[omp_get_thread_num()][2] = j;
					minimals[omp_get_thread_num()][0] = func(x, j);
					flag = true;
				}
				else if (func(x, j) < minimals[omp_get_thread_num()][0]) {
					minimals[omp_get_thread_num()][1] = x;
					minimals[omp_get_thread_num()][2] = j;
					minimals[omp_get_thread_num()][0] = func(x, j);
				}
			}
		}

		

		double* globalMin = new double[3];
		globalMin[0] = minimals[0][0];
		globalMin[1] = minimals[0][1];
		globalMin[2] = minimals[0][2];
		for (int i = 0; i < threads; i++) {
			if (minimals[i][0] < globalMin[0]) {
				globalMin[0] = minimals[i][0];
				globalMin[1] = minimals[i][1];
				globalMin[2] = minimals[i][2];
			}
		}
		
		

		flag = false;
		xIterations = (hStart*2)/hEnd;

#pragma omp parallel for num_threads(threads) firstprivate(flag)
		for (int i = 0; i < xIterations; i++)
		{
			double x = globalMin[1] - hStart + i * hEnd;
			for (double j = bottom; j < top; j += hEnd)
			{
				if (!flag) {
					minimals[omp_get_thread_num()][1] = x;
					minimals[omp_get_thread_num()][2] = j;
					minimals[omp_get_thread_num()][0] = func(x, j);
					flag = true;
				}
				else if (func(x, j) < minimals[omp_get_thread_num()][0]) {
					minimals[omp_get_thread_num()][1] = x;
					minimals[omp_get_thread_num()][2] = j;
					minimals[omp_get_thread_num()][0] = func(x, j);
				}
			}
		}


		globalMin[0] = minimals[0][0];
		globalMin[1] = minimals[0][1];
		globalMin[2] = minimals[0][2];
		for (int i = 0; i < threads; i++) {
			if (minimals[i][0] < globalMin[0]) {
				globalMin[0] = minimals[i][0];
				globalMin[1] = minimals[i][1];
				globalMin[2] = minimals[i][2];
			}
		}

		cout << "Min point x: " << globalMin[1] << " y: " << globalMin[2] << " F: " << globalMin[0] << endl;
		cout << "Time: " << omp_get_wtime() - start << endl;


		system("pause");
}

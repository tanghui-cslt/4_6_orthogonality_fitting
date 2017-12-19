#include <iostream>
using namespace std;
#define Error 1e-3

void  scanf_data(double ** &x_fx, int &m, int &n);		//输入函数
void init_data(double **&p, double *&a, double *&b, double *&d, int m, int n);
void calc_express(double *a, double *b, double *d, double **express, int times, int n);
double  solution(double ** x_fx, int m, int n);
double calc_error(double **x_fx, double *d, double **p, int m, int n);

int main()
{
	double **x_fx = nullptr;					// x,fx的二维数组
	int n = 1;									// 方程组的次数
	int m = 0;									// 观测数据的个数
	scanf_data(x_fx, m, n);
	while (solution(x_fx, m, n) > Error)
	{
		n++;

	}
	cout << "第" << n << "次满足误差" << endl;

	getchar();
	getchar();
	return 0;
}
void init_data(double **&p, double *&a, double *&b, double *&d, int m, int n)
{
	p = new double *[n];
	for (int i = 0; i < n; i++)
	{
		p[i] = new double[m];
	}
	a = new double[n];
	b = new double[n];
	d = new double[n];
}
void  scanf_data(double ** &x_fx, int &m, int &n)
{

	cout << "********请输入数值的个数:**********\n";
	cin >> m;

	cout << "********请输入方程的次数:**********\n";
	cin >> n;

	x_fx = new double *[2];

	cout << "********请先输入" << m << "个x,然后输入" << m << "个f(x):**********\n";


	for (int i = 0; i < 2; i++)
	{
		x_fx[i] = new double[m];
		for (int j = 0; j < m; j++)
		{
			cin >> x_fx[i][j];
		}
	}
}

double fun(double x[], double y[], double pi[], int m)
{
	double sum = 0;
	for (int i = 0; i < m; i++)
	{
		sum += x[i] * y[i] * pi[i];
	}
	return sum;
}

double fun(double  y[], double pi[], int m)//重载函数
{
	double *x = new double[m];
	for (int i = 0; i < m; i++)
	{
		x[i] = 1;
	}
	return fun(x, y, pi, m);
}


double solution(double ** x_fx, int m, int n)
{
	double *a = NULL;
	double *b = NULL;
	double **p = nullptr;
	double *d = nullptr;

	init_data(p, a, b, d, m, n + 1);

	for (int i = 0; i < m; i++)
	{
		p[0][i] = 1;
	}

	d[0] = fun(x_fx[1], p[0], m) / fun(p[0], p[0], m);
	a[1] = fun(x_fx[0], p[0], p[0], m) / fun(p[0], p[0], m);
	for (int i = 0; i < m; i++)
	{
		p[1][i] = (x_fx[0][i] - a[1]);
	}

	//利用 a ,b和 p来求解展开公式

	d[1] = fun(x_fx[1], p[1], m) / fun(p[1], p[1], m);

	for (int i = 1; i < n; i++)
	{
		a[i + 1] = fun(x_fx[0], p[i], p[i], m) / fun(p[i], p[i], m);
		b[i] = fun(p[i], p[i], m) / fun(p[i - 1], p[i - 1], m);

		for (int j = 0; j < m; j++)
		{
			p[i + 1][j] = (x_fx[0][j] - a[i + 1])*p[i][j] - b[i] * p[i - 1][j];
		}
		d[i + 1] = fun(x_fx[1], p[i + 1], m) / fun(p[i + 1], p[i + 1], m);
	}

	//cout << "*****a*****\n";
	//for (int i = 0; i < n + 1; i++)
	//{
	//	cout << a[i] << " ";
	//}
	//cout << endl;
	/*cout << "*****d*****\n";
	for (int i = 0; i < n + 1; i++)
	{
		cout << d[i] << " ";
	}
	cout << endl;

	cout << "\n**p**\n";
	for (int i = 0; i < n + 1; i++)
	{
		for (int j = 0; j < m; j++)
		{
			cout << p[i][j] << " ";
		}
		cout << endl;
	}*/

	double * *express = new double *[n + 1];
	for (int i = 0; i < n + 1; i++)
	{
		express[i] = new double[n + 1];
		for (int j = 0; j < n + 1; j++)
		{
			express[i][j] = 0;
		}
	}

	int times = 0;
	calc_express(a, b, d, express, times, n);

	return calc_error(x_fx, d, p, m, n);

}
void calc_express(double *a, double *b, double *d, double **express, int times, int n)//递归展开多项式
{
	if (times == 0)
	{
		express[0][0] = 1;
	}
	else if (times == 1)
	{
		express[1][1] = express[0][0];
		express[1][0] = (-1)*a[1] * express[0][0];
	}
	else
	{
		for (int i = 0; i <= times; i++)
		{
			if (i > 0) //x*p[times-1]
				express[times][i] += express[times - 1][i - 1];

			if (i < times)
				express[times][i] += (-1)*a[times] * express[times - 1][i];

			if (i < times - 1)
				express[times][i] += (-1)*b[times - 1] * express[times - 2][i];
		}
	}
	if (times == n)
	{

	
		cout << "\n\n!!!!-----------阶数为" << n << "时----------!!!!!\n";
		double *ans = new double[n + 1];
		/*cout << "矩阵Q的表达式为：\n";
		for (int i = 0; i <= times; i++)
		{
			for (int j = 0; j <= times; j++)
			{
				cout << express[i][j] << " ";
			}
			cout << endl;
		}*/

		for (int i = 0; i <= times; i++)
		{
			ans[i] = 0;
			for (int j = 0; j <= times; j++)
			{
				ans[i] += d[j] * express[j][i];

			}

		}
		cout << "\n----拟合多项式为----\n";
		for (int i = 0; i <= times; i++)
		{
			if (i == 0)
			{
				cout << "\t"<<ans[i] << "+";
			}
			else if (i == times)
			{
				cout << ans[i];
				for (int j = 0; j < i; j++)
				{
					if (j == 0)
					{
						cout << "x";
					}
					else
					{
						cout << "*x";
					}
				}
			}
			else
			{
				cout << ans[i];
				for (int j = 0; j < i; j++)
				{
					if (j == 0)
					{
						cout << "x";
					}
					else
					{
						cout << "*x";
					}
				}
				cout << "+";
			}
		}
		cout << endl;

		return;
	}

	calc_express(a, b, d, express, times + 1, n);
}

double calc_error(double **x_fx, double *d, double **p, int m, int n)
{
	double error = 0;

	for (int i = 0; i < m; i++)
	{
		error += pow(x_fx[1][i], 2);
	}

	double temp_error = 0;
	for (int i = 0; i < n + 1; i++)
	{
		temp_error += d[i] * d[i] * fun(p[i], p[i], m);
	}
	cout << "\n误差为 " << error - temp_error << endl;

	return error - temp_error;
}
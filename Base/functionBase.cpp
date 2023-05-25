#include<iostream>
#include<fstream>
#include"functionBase.h"
#include"const.h"

Constant constants;

void changemiu(double x)
{
	constants.changemiu(x);
}

void changevN(double x)
{
	constants.changevN(x);
}

void changechi(double x)
{
	constants.changechi(x);
}

void changen(double x)
{
	constants.changen(x);
}

void changeconst(double miu, double vN, double chi, double n)
{
	changemiu(miu);
	changevN(vN);
	changechi(chi);
	changen(n);
}

double coth(double x)
{
	double ans = 1 / tanh(x);
	return ans;
}

double ln(double x)
{
	return log(x);
}

double L(double beta)
{
	double ans = 0;
	if (beta != 0)
	{
		ans = coth(beta) - 1 / beta;
	}
	return ans;
}

double I3(double beta)
{
	double n = constants.getn();
	double ans = pow(sqrt(n) * L(beta), 3);
	return ans;
}

double J1(double beta)
{
	double n = constants.getn();
	double ans = 3.0 / 2.0 * (n * L(beta) * L(beta) - 1);
	return ans;
}

double J2(double beta)
{
	double n = constants.getn();
	double ans = sqrt(n) * L(beta);
	ans = 3 * ((ans * ans - 1) / 2) * ((ans * ans - 1) / 2);
	return ans;
}

double J3(double beta)
{
	double n = constants.getn();
	double ans = sqrt(n) * L(beta);
	ans = ans * ans - 1;
	ans = ans / 2;
	ans = pow(ans, 3);
	return ans;
}

double E0(double beta)
{
	double ans = J1(beta) / 3;
	return ans;
}

double A0(double beta)
{
	double vN = constants.getvN();
	double chi = constants.getchi();
	double mu = constants.getmiu();
	double x = I3(beta);
	double ans = -1 / vN * ((x - 1.0) * log(x / (x - 1.0)) + chi / x) - mu / vN;
	return ans;
}

double A1(double beta)
{
	double vN = constants.getvN();
	double chi = constants.getchi();
	double mu = constants.getmiu();
	double x = I3(beta);
	double ans = -mu / vN - log(x / (x - 1.0)) / vN + 1 / vN / x + chi / vN / (x * x);
	return ans;
}

double A2(double beta)
{
	double vN = constants.getvN();
	double chi = constants.getchi();
	double mu = constants.getmiu();
	double x = I3(beta);
	double ans = (-2.0 * chi * x + 2.0 * chi + x) / (x - 1.0) / vN / (x * x * x);
	ans = ans / 2;
	return ans;
}

double A3(double beta)
{
	double vN = constants.getvN();
	double chi = constants.getchi();
	double mu = constants.getmiu();
	double x = I3(beta);
	double ans = ((6.0 * chi - 3.0) * x * x + (-12.0 * chi + 2.0) * x + 6.0 * chi) / pow(x - 1.0, 2.0) / vN / (x * x * x * x);
	ans = ans / 6;
	return ans;
}

double a0(double beta)
{
	double A_0 = A0(beta);
	double A_1 = A1(beta);
	double A_2 = A2(beta);
	double A_3 = A3(beta);
	double I30 = I3(beta);
	double ans = A_0 - A_1 * I30 + A_2 * pow(I30, 2) - A_3 * pow(I30, 3);
	return ans;
}

double a1(double beta)
{
	double A_1 = A1(beta);
	double A_2 = A2(beta);
	double A_3 = A3(beta);
	double I30 = I3(beta);
	double ans = A_1 - 2 * A_2 * I30 + 3 * A_3 * pow(I30, 2);
	return ans;
}

double a2(double beta)
{
	double A_2 = A2(beta);
	double A_3 = A3(beta);
	double I30 = I3(beta);
	double ans = A_2 - 3 * A_3 * I30;
	return ans;
}

double a3(double beta)
{
	return A3(beta);
}

double H1(double beta)
{
	double n = constants.getn();
	double ans = 3.0 * n * (1 / tanh(beta) - 1 / beta) * (1.0 - 1 / pow(tanh(beta), 2.0) + 1 / (beta * beta));
	return ans;
}

double H2(double beta)
{
	double n = constants.getn();
	double t2 = 1 / tanh(beta);
	double t3 = t2 * t2;
	double t4 = n * t3;
	double t6 = beta * beta;
	double t7 = 1 / t6;
	double t10 = t3 * t3;
	double t15 = t6 * t6;
	double t19 = n * t2;
	double t24 = 1 / beta;
	double t31 = 3.0 * n - 12.0 * t4 + 6.0 * n * t7 + 9.0 * n * t10 - 6.0 * t4 * t7 + 9.0 * n / t15 - 6.0 * t19 / beta / t6 + 6.0 * t19 * t24 - 6.0 * n * t2 * t3 * t24;
	return t31;
}

double H3(double beta)
{
	double n = constants.getn();
	double t1 = 1 / tanh(beta);
	double t2 = n * t1;
	double t4 = t1 * t1;
	double t6 = n * t1 * t4;
	double t8 = beta * beta;
	double t10 = 1 / beta / t8;
	double t13 = t4 * t4;
	double t20 = 1 / t8;
	double t25 = t8 * t8;
	double t33 = 1 / beta;
	double t34 = n * t33;
	double t41 = -24.0 * t2 + 60.0 * t6 - 18.0 * n * t10 - 36.0 * n * t1 * t13 + 18.0 * n * t4 * t10 - 18.0 * t2 * t20 + 18.0 * t6 * t20 - 36.0 * n / beta / t25 + 18.0 * t2 / t25 + 6.0 * t34 - 24.0 * t34 * t4 + 18.0 * n * t13 * t33;
	return t41;
}

double c1(double beta)
{
	double h1 = H1(beta);
	double ans = 1.0 / h1;
	return ans;
}

double c2(double beta)
{
	double h1 = H1(beta);
	double h2 = H2(beta);
	double h3 = H3(beta);
	double ans = -h2 / pow(h1, 3);
	return ans;
}

double c3(double beta)
{
	double h1 = H1(beta);
	double h2 = H2(beta);
	double h3 = H3(beta);
	double ans = (3 * pow(h2, 2) - h1 * h3) / pow(h1, 5);
	return ans;
}

double b1(double beta)
{
	double n = constants.getn();
	double ans = 1 / tanh(beta) + 1 / beta + beta - beta / pow(tanh(beta), 2.0) - cosh(beta) / sinh(beta);
	ans = n * ans;
	return ans;
}

double b2(double beta)
{
	double n = constants.getn();
	double ans = 1.0 - 2.0 / pow(tanh(beta), 2.0) - 1 / (beta * beta) - 2.0 / tanh(beta) * beta + 2.0 * beta / pow(tanh(beta), 3.0) + pow(cosh(beta), 2.0) / pow(sinh(beta), 2.0);
	ans = n * ans;
	return ans;
}

double b3(double beta)
{
	double n = constants.getn();
	double ans = -6.0 / tanh(beta) + 6.0 / pow(tanh(beta), 3.0) + 2.0 / (beta * beta * beta) - 2.0 * beta + 8.0 * beta / pow(tanh(beta), 2.0) - 6.0 * beta / pow(tanh(beta), 4.0) - 2.0 * pow(cosh(beta), 3.0) / pow(sinh(beta), 3.0) + 2.0 * cosh(beta) / sinh(beta);
	ans = n * ans;
	return ans;
}

double B1(double beta)
{
	double ans = b1(beta) * c1(beta);
	return ans;
}

double B2(double beta)
{
	double ans = b2(beta) * c1(beta) * c1(beta) + b1(beta) * c2(beta);
	ans = ans / 2;
	return ans;
}

double B3(double beta)
{
	double ans = b3(beta) * pow(c1(beta), 3) + 3 * b2(beta) * c1(beta) * c2(beta) + b1(beta) * c3(beta);
	ans = ans / 6;
	return ans;
}

double alpha(double beta)
{
	double a_1 = a1(beta);
	double a_2 = a2(beta);
	double a_3 = a3(beta);

	double B_1 = B1(beta);
	double B_2 = B2(beta);
	double B_3 = B3(beta);

	double J10 = J1(beta);

	double ans = 0;
	ans = a_1 + 2 * a_2 + 3 * a_3 + B_1 - 2 * B_2 * J10 + 3 * B_3 * pow(J10, 2);
	return ans;
}

double lambda(double beta)
{
	double a_1 = a1(beta);
	double a_2 = a2(beta);
	double a_3 = a3(beta);

	double B_1 = B1(beta);
	double B_2 = B2(beta);
	double B_3 = B3(beta);

	double J10 = J1(beta);

	double ans = 0;
	ans = 9 * a_3 + 4 * a_2 + a_1 + 2 * B_2 - 6 * B_3 * J10;
	return ans;
}

double G(double beta)
{
	double a_1 = a1(beta);
	double a_2 = a2(beta);
	double a_3 = a3(beta);

	double ans = 0;
	ans = -3 * a_3 - 2 * a_2 - a_1;
	return ans;
}

double l(double beta)
{
	double a_1 = a1(beta);
	double a_2 = a2(beta);
	double a_3 = a3(beta);

	double B_3 = B3(beta);

	double ans = 0;
	ans = 9.0 / 2.0 * a_3 - 0.5 * a_1 + 3 * B_3;
	return ans;
}

double m(double beta)
{
	double a_1 = a1(beta);
	double a_3 = a3(beta);

	double ans = -3 * a_3 + a_1;
	return ans;
}

double n(double beta)
{

	double a_1 = a1(beta);
	double a_2 = a2(beta);
	double a_3 = a3(beta);

	double ans = 0;
	ans = 12 * a_3 + 8 * a_2 + 4 * a_1;
	return ans;
}

double Sigma(double beta)
{
	double sigma = 0;
	double ALPHA = alpha(beta);
	double LAMBDA = lambda(beta);
	double g = G(beta);
	double L = l(beta);
	double N = n(beta);
	double e0 = E0(beta);
	sigma = ALPHA + (-ALPHA + 3 * LAMBDA + 2 * g) * e0 + (1.5 * ALPHA - 3 * LAMBDA - 2 * g + 9 * L + N) * pow(e0, 2);
	return sigma;
}

double Solve(double SIGMA)
{
	double ans = 0;
	double beta = 0;
	double delta = 0.000001;
	double I30 = I3(beta);
	double sigma = 0;
	bool T = true;
	while (I30 < 1)
	{
		beta = beta + 0.01;
		I30 = I3(beta);
	}
	sigma = Sigma(beta) - SIGMA;
	while (T)
	{
		if (sigma > 0)
		{
			beta = beta - delta;
			sigma = Sigma(beta) - SIGMA;
			T = sigma > 0;
		}
		else
		{
			beta = beta + delta;
			sigma = Sigma(beta) - SIGMA;
			T = sigma <= 0;
		}
	}
	ans = beta;
	return ans;
}

double W(double beta)
{
	double n = constants.getn();
	double miu = constants.getmiu();
	double vN = constants.getvN();
	double chi = constants.getchi();
	double ans = 0;
	double ans_1 = n * (beta * L(beta) + ln(beta / sinh(beta)));
	double I = I3(beta);
	double ans_2 = -1 / vN * ((I - 1) * ln(I / (I - 1)) + chi / I) - miu / vN * (I - 1);
	if (I == 1)
	{
		ans = 0;
	}
	else
	{
		ans = ans_1 + ans_2;
	}
	return ans;
}

double w(double beta)
{
	double ans = 0;
	ans = alpha(beta);
	return ans;
}
//Ã»±àÍê¡£

double solve()
{
	double X = 0;
	double delta = 0.0001;
	double beta = ln(1 + X);
	double betabefore = 0;
	double betaafter = 0;
	double ans = 0;
	bool T = true;

	while (T)
	{
		T = I3(beta) < 1;
		X = X + delta;
		beta = ln(1 + X);
	}

	T = true;
	while (T)
	{
		X = X + delta;
		beta = ln(1 + X);
		betabefore = ln(1 + X - delta);
		betaafter = ln(1 + X + delta);
		T = (W(betaafter) - W(beta)) * (W(beta) - W(betabefore)) > 0;
		//std::cout << I3(beta) << "   " << (W(betaafter) - W(beta)) * (W(beta) - W(betabefore)) << std::endl;
	}
	std::cout << "Finishied" << std::endl;
	//system("pause");
	//std::cout << std::endl << std::endl;
	return beta;
}

void output()
{
	using namespace std;
	ofstream outfile;
	ofstream outconstants;
	int i = 0;
	double x = 0;
	double miu = -pow(10, x);
	double delta = 0.01;
	double ans = 0;
	double ans1 = 0;
	double ans2 = 0;
	double ans3 = 0;

	miu = -delta;
	//outfile.open("data//mu.csv");
	outconstants.open("data//const.csv");
	while (miu > -2)
	{
		changemiu(miu);
		ans = Solve(0);
		//outfile << -x << "," << log10(I3(ans1) - 1) << "," << log10(I3(ans2) - 1) << "," << log10(I3(ans3) - 1) << endl;
		outconstants << miu << "," << lambda(ans) << "," << G(ans) << "," << l(ans) << "," << m(ans) << "," << n(ans) << endl;
		//cout << -x << " Finished " << Sigma(ans1) << "  " << Sigma(ans2) << "  " << Sigma(ans3) << endl;
		cout << miu << " Finished " << Sigma(ans) << endl;
		x = x - delta;
		miu = miu - delta;
	}
	//outfile.close();
	outconstants.close();
}

void outLambdaAboutChi()
{
	using namespace std;
	ofstream outfile;
	double delta = 0.01;
	double chi = 0.1;

	double miu[5] = { -0.4,-0.6,-0.8,-1.0,-1.2 };
	double ans[5];

	int i = 0;

	outfile.open("data//lambda.csv");
	while (chi < 0.5)
	{
		changechi(chi);
		outfile << chi << ",";
		i = 0;
		while (i < 5)
		{
			changemiu(miu[i]);
			ans[i] = Solve(0);
			outfile << lambda(ans[i]) << ",";
			i = i++;
		}
		outfile << endl;
		cout << "chi=" << chi << " Finished" << endl;
		chi = chi + delta;
	}
	outfile.close();
}

void GaboutChi()
{
	using namespace std;
	ofstream outfile;
	double delta = 0.01;
	double chi = 0.1;

	double miu[5] = { -0.4,-0.6,-0.8,-1.0,-1.2 };
	double ans[5];

	int i = 0;

	outfile.open("data//G.csv");
	while (chi < 0.5)
	{
		changechi(chi);
		outfile << chi << ",";
		i = 0;
		while (i < 5)
		{
			changemiu(miu[i]);
			ans[i] = Solve(0);
			outfile << G(ans[i]) << ",";
			i = i++;
		}
		outfile << endl;
		cout << "chi=" << chi << " Finished" << endl;
		chi = chi + delta;
	}
	outfile.close();
}

void lAboutChi()
{
	using namespace std;
	ofstream outfile;
	double delta = 0.01;
	double chi = 0.1;

	double miu[5] = { -0.4,-0.6,-0.8,-1.0,-1.2 };
	double ans[5];

	int i = 0;

	outfile.open("data//l.csv");
	while (chi < 0.5)
	{
		changechi(chi);
		outfile << chi << ",";
		i = 0;
		while (i < 5)
		{
			changemiu(miu[i]);
			ans[i] = Solve(0);
			outfile << l(ans[i]) << ",";
			i = i++;
		}
		outfile << endl;
		cout << "chi=" << chi << " Finished" << endl;
		chi = chi + delta;
	}
	outfile.close();
}

void mAboutChi()
{
	using namespace std;
	ofstream outfile;
	double delta = 0.01;
	double chi = 0.1;

	double miu[5] = { -0.4,-0.6,-0.8,-1.0,-1.2 };
	double ans[5];

	int i = 0;

	outfile.open("data//m.csv");
	while (chi < 0.5)
	{
		changechi(chi);
		outfile << chi << ",";
		i = 0;
		while (i < 5)
		{
			changemiu(miu[i]);
			ans[i] = Solve(0);
			outfile << m(ans[i]) << ",";
			i = i++;
		}
		outfile << endl;
		cout << "chi=" << chi << " Finished" << endl;
		chi = chi + delta;
	}
	outfile.close();
}

void nAboutChi()
{
	using namespace std;
	ofstream outfile;
	double delta = 0.01;
	double chi = 0.1;

	double miu[5] = { -0.4,-0.6,-0.8,-1.0,-1.2 };
	double ans[5];

	int i = 0;

	outfile.open("data//n.csv");
	while (chi < 0.5)
	{
		changechi(chi);
		outfile << chi << ",";
		i = 0;
		while (i < 5)
		{
			changemiu(miu[i]);
			ans[i] = Solve(0);
			outfile << n(ans[i]) << ",";
			i = i++;
		}
		outfile << endl;
		cout << "chi=" << chi << " Finished" << endl;
		chi = chi + delta;
	}
	outfile.close();
}

void any()
{
	using namespace std;
	ofstream outfile;

	double numberOfLinks = 100;
	double ans = 0;

	outfile.open("data//number.csv");
	while (numberOfLinks < 10000)
	{
		changen(numberOfLinks);
		numberOfLinks = numberOfLinks + 10;
		ans = Solve(0);
		outfile << numberOfLinks << "," << I3(ans) << endl;
		cout << numberOfLinks << "  Down" << endl;
	}
	outfile.close();
}

void outLambdaAboutVn()
{
	using namespace std;
	ofstream outfile;
	double delta = 0.001;
	double vN = 0.01;

	double miu[5] = { -0.4,-0.6,-0.8,-1.0,-1.2 };
	double ans[5];

	int i = 0;

	outfile.open("data//lambdaAboutVn.csv");
	while (vN <= 0.1)
	{
		changevN(vN);
		outfile << vN << ",";
		i = 0;
		while (i < 5)
		{
			changemiu(miu[i]);
			ans[i] = Solve(0);
			if (i != 4)
			{
				outfile << lambda(ans[i]) << ",";
			}
			else
			{
				outfile << lambda(ans[i]);
			}
			i = i++;
		}
		outfile << endl;
		cout << "vN= " << vN << " Finished" << endl;
		vN = vN + delta;
	}
	outfile.close();
}

void GAboutVn()
{
	using namespace std;
	ofstream outfile;
	double delta = 0.001;
	double vN = 0.01;

	double miu[5] = { -0.4,-0.6,-0.8,-1.0,-1.2 };
	double ans[5];

	int i = 0;

	outfile.open("data//GAboutVn.csv");
	while (vN <= 0.1)
	{
		changevN(vN);
		outfile << vN << ",";
		i = 0;
		while (i < 5)
		{
			changemiu(miu[i]);
			ans[i] = Solve(0);
			if (i != 4)
			{
				outfile << G(ans[i]) << ",";
			}
			else
			{
				outfile << G(ans[i]);
			}
			i = i++;
		}
		outfile << endl;
		cout << "vN= " << vN << " Finished" << endl;
		vN = vN + delta;
	}
	outfile.close();
}

void lAboutVn()
{
	using namespace std;
	ofstream outfile;
	double delta = 0.001;
	double vN = 0.01;

	double miu[5] = { -0.4,-0.6,-0.8,-1.0,-1.2 };
	double ans[5];

	int i = 0;

	outfile.open("data//lAboutVn.csv");
	while (vN <= 0.1)
	{
		changevN(vN);
		outfile << vN << ",";
		i = 0;
		while (i < 5)
		{
			changemiu(miu[i]);
			ans[i] = Solve(0);
			if (i != 4)
			{
				outfile << l(ans[i]) << ",";
			}
			else
			{
				outfile << l(ans[i]);
			}
			i = i++;
		}
		outfile << endl;
		cout << "vN= " << vN << " Finished" << endl;
		vN = vN + delta;
	}
	outfile.close();
}

void mAboutVn()
{
	using namespace std;
	ofstream outfile;
	double delta = 0.001;
	double vN = 0.01;

	double miu[5] = { -0.4,-0.6,-0.8,-1.0,-1.2 };
	double ans[5];

	int i = 0;

	outfile.open("data//mAboutVn.csv");
	while (vN <= 0.1)
	{
		changevN(vN);
		outfile << vN << ",";
		i = 0;
		while (i < 5)
		{
			changemiu(miu[i]);
			ans[i] = Solve(0);
			if (i != 4)
			{
				outfile << m(ans[i]) << ",";
			}
			else
			{
				outfile << m(ans[i]);
			}
			i = i++;
		}
		outfile << endl;
		cout << "vN= " << vN << " Finished" << endl;
		vN = vN + delta;
	}
	outfile.close();
}

void nAboutVn()
{
	using namespace std;
	ofstream outfile;
	double delta = 0.001;
	double vN = 0.01;

	double miu[5] = { -0.4,-0.6,-0.8,-1.0,-1.2 };
	double ans[5];

	int i = 0;

	outfile.open("data//nAboutVn.csv");
	while (vN <= 0.1)
	{
		changevN(vN);
		outfile << vN << ",";
		i = 0;
		while (i < 5)
		{
			changemiu(miu[i]);
			ans[i] = Solve(0);
			if (i != 4)
			{
				outfile << n(ans[i]) << ",";
			}
			else
			{
				outfile << n(ans[i]);
			}
			i = i++;
		}
		outfile << endl;
		cout << "vN= " << vN << " Finished" << endl;
		vN = vN + delta;
	}
	outfile.close();
}

void outWaboutMiuAndI3()
{
	using namespace std;

	ofstream outfile;

	double delta = 0.01;
	double miu = -3;
	double deltabeta = 0.0001;
	double beta = 0;
	double I = I3(beta);
	double w = 0;
	bool T = I > 1 && I < 10;

	outfile.open("data//W.csv");
	while (miu < -0.001)
	{
		changemiu(miu);
		beta = deltabeta;
		I = I3(beta);
		while (I < 11)
		{
			T = I > 1 && I < 11;
			if (T)
			{
				w = W(beta);
				outfile << miu << "," << I << "," << w << endl;
				//cout << "miu= " << miu << " V/V0=" << I << " W=" << w << " Finished" << endl;
			}
			beta = beta + deltabeta;
			I = I3(beta);
		}
		//cout << miu << endl;
		miu = miu + delta;
	}
	outfile.close();
}

void outVaboutMiu()
{
	using namespace std;

	ofstream outfile;
	double miu = 0;
	double ans1 = 0;
	double ans2 = 0;
	double ans3 = 0;
	double delta = 0.1;
	double X = 0.5;
	int i = 0;
	
	outfile.open("data/V.csv");
	while (i < 60)
	{
		miu = -pow(10, X);
		changemiu(miu);
		changevN(0.001);
		ans1 = Solve(0);
		ans1 = I3(ans1);
		ans1 = log10(ans1 - 1);

		changevN(0.01);
		ans2 = Solve(0);
		ans2 = I3(ans2);
		ans2 = log10(ans2 - 1);

		changevN(0.1);
		ans3 = Solve(0);
		ans3 = I3(ans3);
		ans3 = log10(ans3 - 1);

		outfile << -miu << "," << -X << "," << ans1 << "," << ans2 << "," << ans3 << endl;
		cout << "X=" << X << endl;
		X = X - delta;
		i = i++;
	}
	outfile.close();
}

double E(double beta)
{
	double LAMBDA = lambda(beta);
	double g = G(beta);
	double ans = g * (3 * LAMBDA + 2 * g) / (LAMBDA + g);
	return ans;
}

void outEaboutMiu()
{
	using namespace std;
	ofstream outfile;
	double miu = -3.0;
	double delta = 0.01;
	double beta;
	int i = 0;

	outfile.open("data/EaboutMiu.csv");
	while (i < 300)
	{
		changemiu(miu);
		beta = Solve(0);
		outfile << miu << "," << E(beta) << endl;
		miu = miu + delta;
		i = i++;
	}
	outfile.close();
}

double NU(double beta)
{
	double ans = 0;
	ans = lambda(beta) / 2;
	ans = ans / (lambda(beta) + G(beta));
	return ans;
}

double constB1(double beta)
{
	double g = G(beta);
	double LAMBDA = lambda(beta);
	double L = l(beta);
	double M = m(beta);
	double N = n(beta);
	double nu = NU(beta);

	double ans1 = LAMBDA * (1 - 2 * nu);
	double ans2 = 3 * g - g * nu * nu;
	double ans3 = 2 * M * (1 - 2 * nu * nu - nu);
	double ans4 = N * (nu + nu * nu);

	double ans = ans1 + ans2 + ans3 + ans4;

	return ans;
}

double constB2(double beta)
{
	double g = G(beta);
	double LAMBDA = lambda(beta);
	double L = l(beta);
	double M = m(beta);
	double N = n(beta);
	double nu = NU(beta);

	double ans1 = -LAMBDA * (0.5 + nu * nu);
	double ans2 = -g * nu * nu;
	double ans3 = -L * pow((1 - 2 * nu), 2);
	double ans4 = -2 * M * nu * (1 + nu) + N * nu;

	double ans = ans1 + ans2 + ans3 + ans4;

	return ans;
}

double DELTA(double beta)
{
	double g = G(beta);
	double LAMBDA = lambda(beta);
	double b1 = constB1(beta);
	double b2 = constB2(beta);

	double ans = g * b2 - (LAMBDA + g) * b1;
	ans = ans / (g * (3 * LAMBDA + 2 * g));

	return ans;
}

void DELTAaboutMiu()
{
	using namespace std;
	ofstream outfile;
	double miu = -3.0;
	double delta = 0.01;
	double beta;
	int i = 0;

	outfile.open("data/DELTAaboutMiu.csv");
	while (i < 300)
	{
		changemiu(miu);
		beta = Solve(0);
		outfile << miu << "," << DELTA(beta) << endl;
		miu = miu + delta;
		i = i++;
	}
	outfile.close();
}

void DELTAaboutChi()
{
	using namespace std;
	ofstream outfile;
	double miu[10] = { -1.0,-0.9,-0.8,-0.7,-0.6 ,-0.5,-0.4,-0.3,-0.2,-0.1 };
	double chi = 0.1;
	double delta = 0.01;
	double beta = 0;
	double ans[10];
	int i = 0;
	int j = 0;

	outfile.open("data/delta-chi.csv");
	while (i < 40)
	{
		changechi(chi);
		outfile << chi;
		j = 0;
		while (j < 10)
		{
			changemiu(miu[j]);
			ans[j] = Solve(0);
			ans[j] = DELTA(ans[j]);
			outfile << "," << ans[j];
			j = j++;
		}
		outfile << endl;
		chi = chi + delta;
		i = i++;
	}
	outfile.close();
}

void DELTAaboutVn()
{
	using namespace std;
	ofstream outfile;
	double miu[10] = { -1.0,-0.9,-0.8,-0.7,-0.6 ,-0.5,-0.4,-0.3,-0.2,-0.1 };
	double vN = 0.01;
	double delta = 0.001;
	double beta = 0;
	double ans[10];
	int i = 0;
	int j = 0;

	outfile.open("data/delta-vn.csv");
	while (i < 90)
	{
		changevN(vN);
		outfile << vN;
		j = 0;
		while (j < 10)
		{
			changemiu(miu[j]);
			ans[j] = Solve(0);
			ans[j] = DELTA(ans[j]);
			outfile << "," << ans[j];
			j = j++;
		}
		outfile << endl;
		vN = vN + delta;
		i = i++;
	}
	outfile.close();
}

void strainL()
{
	using namespace std;
	ofstream outfile;
	double ans = Solve(0);
	double delta = DELTA(ans);
	double e = E(ans);
	double T = 0; 
	double deltaT = 1;
	double strain = T / e + delta * pow(T / e, 2);

	outfile.open("data/T.csv");
	while (T < 1000)
	{
		strain = T / e + delta * pow(T / e, 2);
		outfile << T << "," << strain << "," << T / e << endl;
		T = T + deltaT;
	}
	outfile.close();
}
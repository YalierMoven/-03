#pragma once

#include<iostream>

class Constant
{
public:

	Constant();
	Constant(double, double, double, double);

	double getmiu();
	double getvN();
	double getchi();
	double getn();

	void changemiu(double);
	void changevN(double);
	void changechi(double);
	void changen(double);

private:

	double miu;//Protential miu_wave = μ/kT
	double vN;//vN放在一起
	double chi;//χ
	double n;//每条链上的链数

};
#include<iostream>
#include"const.h"

Constant::Constant()
{
	miu = -1;//Protential miu_wave = μ/kT
	vN = 0.001;//vN放在一起
	chi = 0.4;//χ
	n = 1000;//每条链上的链数
}

Constant::Constant(double MIU, double Vn, double CHI, double N)
{
	miu = MIU;
	vN = Vn;
	chi = CHI;
	n = N;
}

double Constant::getmiu()
{
	return miu;
}

double Constant::getvN()
{
	return vN;
}

double Constant::getchi()
{
	return chi;
}

double Constant::getn()
{
	return n;
}


void Constant::changemiu(double x)
{
	miu = x;
}

void Constant::changevN(double x)
{
	vN = x;
}

void Constant::changechi(double x)
{
	chi = x;
}

void Constant::changen(double x)
{
	n = x;
}
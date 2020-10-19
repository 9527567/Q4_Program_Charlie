/*
		���ߣ�Charlie_Xia
		���ܣ��򵥵�����Ԫ����

*/

#include <iostream>
#include <vector>
#include "Elemstiffness.hpp"
#include "readData.hpp"
#include <cmath>
#include <numeric>
#include<math.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdio.h>
#include<stdlib.h>
using namespace std;


void ex1()
{
	vector<vector<double>> nodes;			// ����ڵ���������
	nodes = readData_1();					// ��ȡ�ڵ�����
	vector<vector<double>> loads;			// �����غ�����
	loads = readData_2();					// 
	vector<vector<int>> elements;		// ���嵥Ԫ�Ľڵ���
	elements = readData_3();				// 
	vector<vector<double>> bounds;			// ����߽�����
	bounds = readData_4();					// 

	//		���ݽڵ�����������նȾ���
	//		��ʼ������նȾ���

	vector<double> b(2 * nodes.size(),0);
	vector<vector<double>> K(2 * nodes.size(),b);

	// ���b vector ���һ����ڴ�
	vector<double>().swap(b);

	// ��ʼ����Ԫ�նȾ���
	for (int j = 0; j < 8; j++)
	{
		b.push_back(0);
	}
	vector<vector<double>>Ke(8,b);
	// ���b vector ���һ����ڴ�
	vector<double>().swap(b);
	
	for (int i = 0; i < sizeof(elements); i++)				//k����Ԫ��k�ε�Ԫ�նȾ���
	{
		// ��Ԫ�նȾ���
		Ke = ElementStiffness(30000000, 0.3, 1, nodes[elements[i][0] - 1][0], nodes[elements[i][0] - 1][1], 
			nodes[elements[i][1] - 1][0], nodes[elements[i][1] - 1][1], nodes[elements[i][2] - 1][0], 
			nodes[elements[i][2] - 1][1], nodes[elements[i][3] - 1][0], nodes[elements[i][3] - 1][1]);

		//���ϵ�����նȾ�����
		for (int j = 0; j < 4; j++)
		{
			for (int p = 0; p < 4; p++)
			{
				K[2 *  elements[i][j] - 2][2 * elements[i][p] - 2] = Ke[2 * j][2 * p] + K[2 * elements[i][j] - 2][2 * elements[i][p] - 2];
				K[2 * elements[i][j] - 2][2 * elements[i][p] - 2 + 1] = Ke[2 * j][2 * p + 1] + K[2 * elements[i][j] - 2][2 * elements[i][p] - 2 + 1];
				K[2 * elements[i][j] - 2 + 1][2 * elements[i][p] - 2] = Ke[2 * j + 1][2 * p] + K[2 * elements[i][j] - 2 + 1][2 * elements[i][p] - 2];
				K[2 * elements[i][j] - 2 + 1][2 * elements[i][p] - 2 + 1] = Ke[2 * j + 1][2 * p + 1] + K[2 * elements[i][j] - 2 + 1][2 * elements[i][p] - 2 + 1];
			}
		}
	}
	// �鿴�նȾ���ʵ��������С�뵱ǰʵ�ʰ�����Ԫ�ظ���
	cout << K.capacity() << " " <<K.size();
	
	/*
	//����߽����������öԽ�Ԫ�ظ�1���Զ�����նȾ���K���غ���������
	for (int i = 0; i < sizeof(bounds); i++)
	{
		if (bounds[i][0] == 0)
			0;//������
		else
		{
			K[2 * i][2 * i] = 1;
			loads[i][0] = 0;
			for (int j = 0; j < 2 * i; j++)
				K[2 * i][j] = 0;
			for (int j = 2 * i + 1; j < 2 * sizeof(bounds); j++)
				K[2 * i][j] = 0;

			for (int j = 0; j < 2 * i; j++)
				K[j][2*i] = 0;
			for (int j = 2 * i + 1; j < 2 * sizeof(bounds); j++)
				K[j][2*i] = 0;
		}
		if (bounds[i][1] == 0)
			0;//������
		else	
		{
			K[2 * i + 1][2 * i + 1] = 1;
			loads[i][1] = 0;
			for (int j = 0; j < 2 * i + 1; j++)
				K[2 * i + 1][j] = 0;
			for (int j = 2 * i + 2; j < 2 * sizeof(bounds); j++)
				K[2 * i + 1][j] = 0;

			for (int j = 0; j < 2 * i + 1; j++)
				K[j][2 * i + 1] = 0;
			for (int j = 2 * i + 2; j < 2 * sizeof(bounds); j++)
				K[j][2 * i + 1] = 0;
		}

	}



	//�ⷽ��
	
	vector<int> F;
	
	for (int i = 0; i < 2 * sizeof(loads)-1; i=i+2)  //�Ѷ�ά����loadsת��Ϊһά����F���������
	{
		F[i] = loads[i / 2][0];
		F[i+1]= loads[i / 2][1];
	}
	
	double x[2 * sizeof(nodes)][2 * sizeof(nodes)];
	for (int i = 0; i < 2 * sizeof(nodes); i++) 
		x[0][i]=0;
	
	double e;
	e = 0.02;
	void cau();
	{

		for (int k = 1; k <= 10; k++)
		{
			for (int i = 0; i < 2 * sizeof(nodes); i++)
			{
				x[k][i] = 1.0 / K[i][i];
				double re = F[i];
				for (int j = 0; j < 2 * sizeof(nodes); j++)
				{
					if (j != i)
						re -= K[i][j] * x[k - 1][j];
				}
				x[k][i] *= re;
			}
			for (int i = 0; i < 2 * sizeof(nodes); i++)
			{
				cout << x[k][i] << " ";
				cout << endl;
			}
			bool judge = true;
			for (int i = 0; i < 2 * sizeof(nodes); i++)
				if (fabs(x[k - 1][i] - x[k][i]) > e) {
					judge = false; break;
				}
			if (judge == true) return;

		}
	}
	*/


}








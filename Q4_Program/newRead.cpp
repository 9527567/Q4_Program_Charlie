#pragma once
#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<string>
#include<stdio.h>
#include "newRead.hpp"

using namespace std;


//��һ������ȡ�ڵ�����
struct dataX_r node()
{
    dataX_r c;
    //���Ȱ������ַ������뵽�ַ���A��
    vector<string> A;
    string a;
    ifstream infile;
    infile.open(".\\mesh.k", ios::in);
    if (!infile)
    {
        cout << "fail to open file" << endl;
    }
    int i = 0;
    while (getline(infile, a))     //k�ļ��е������ַ����������ַ������Ѿ�������A
    {
        A.push_back(a);
    }



    //��ȡ�ڵ�����
    //�Ȱѽڵ�������ַ�������
    vector<string> tempnodes;
    for (int i = 0; i < A.size(); i++)
    {
        if (A[i] == "*NODE")
        {
            for (int j = i + 1; j < A.size(); j++)
            {
                if (A[j] == "*ELEMENT_SHELL")
                    break;
                tempnodes.push_back(A[j]);

            }
        }
    }
  
   //��ʼ��nodes
    vector<double> b(2, 0);
    for (int i = 0; i < tempnodes.size(); i++)
    {
        c.nodes.push_back(b);
    }
    
    //�ѽڵ�������ַ����������
    for (int i =0; i < tempnodes.size(); i++)
    {
       string temp = tempnodes[i];
        string x = temp.substr(9, 15);
        string y=temp.substr(25, 15);
        c.nodes[i][0] = stod(x);
        c.nodes[i][1] = stod(y);
    }

    //��ȡÿ����Ԫ�Ľڵ����
    //ȡÿ����Ԫ�Ľڵ����
    vector<string> tempelements;
    for (int i = 0; i < A.size(); i++)
    {
        if (A[i] == "*ELEMENT_SHELL")
        {
            for (int j = i + 1; j < A.size(); j++)
            {
                if (A[j] == "*LOAD_NODE_POINT")
                    break;
                tempelements.push_back(A[j]);

            }
        }
    }
   

    //��ʼ��elements
    vector<double> b1(4, 0);
    for (int i = 0; i < tempelements.size(); i++)
    {
        c.elements.push_back(b1);
    }

    //�ѵ�Ԫ�ڵ���ŵ��ַ����������
    for (int i = 0; i < tempelements.size(); i++)
    {
        string temp = tempelements[i];
        string n1 = temp.substr(21, 3);
        string n2 = temp.substr(29, 3);
        string n3 = temp.substr(37, 3);
        string n4 = temp.substr(45, 3);
        c.elements[i][0] = stod(n1);
        c.elements[i][1] = stod(n2);
        c.elements[i][2] = stod(n3);
        c.elements[i][3] = stod(n4);
        /*cout << c.elements[i][0] << "  " << c.elements[i][1] << "  " << c.elements[i][2] << "  " << c.elements[i][3] << endl;*/
    }

   
    //��ȡʩ�ӵ�����

    vector<string> temploads;
    for (int i = 0; i < A.size(); i++)
    {
        if (A[i] == "*LOAD_NODE_POINT")
        {
            for (int j = i + 1; j < A.size(); j++)
            {
                if (A[j] == "*BOUNDARY_SPC_NODE")
                    break;
                temploads.push_back(A[j]);

            }
        }
    }

    //��ʼ��loads
    vector<double> b2(2, 0);
    for (int i = 0; i < c.nodes.size(); i++)
    {
        c.loads.push_back(b2);
    }

    //��loads�ַ����������
    for (int i = 2; i < temploads.size(); i++)
    {
        string temp = temploads[i];
        string n = temp.substr(7, 3);
        string x = temp.substr(36, 4);
        string y = temp.substr(49, 1);
        int num = stoi(n);
        c.loads[num-1][0] = stod(x);
        c.loads[num-1][1] = stod(y);
       
    }

    /*for (int i = 0; i < c.loads.size(); i++)
    {
        cout << c.loads[i][0] << "  " << c.loads[i][1] << endl;
    }*/
    


    //��ȡ�߽�����

    vector<string> tempbounds;
    for (int i = 0; i < A.size(); i++)
    {
        if (A[i] == "*BOUNDARY_SPC_NODE")
        {
            for (int j = i + 1; j < A.size(); j++)
            {
                if (A[j] == "*END")
                    break;
                tempbounds.push_back(A[j]);

            }
        }
    }


    //��ʼ��bounds
    vector<double> b3(2, 0);
    for (int i = 0; i < c.nodes.size(); i++)
    {
        c.bounds.push_back(b3);
    }

 
    //��bounds�ַ����������
    for (int i = 2; i < tempbounds.size(); i++)
    {
        string temp = tempbounds[i];
        string n = temp.substr(7, 3);
        string x = temp.substr(29, 1);
        string y = temp.substr(39, 1);
        int num = stoi(n);
        if (x!=" ")
        c.bounds[num - 1][0] = stod(x);
        if (y!=" ")
        c.bounds[num - 1][1] = stod(y);

    }

   

    return c;
}
















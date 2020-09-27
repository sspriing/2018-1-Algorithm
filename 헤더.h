#include <iostream>
#include <string>
#include <fstream>
#include <time.h>
#include <ctime>
#include <algorithm>
#include <vector>
#include "WELL1024a.h"

using namespace std;


/*BWT �˰��򿡼� ���� Ŭ����*/
class BWTmatrix {
	int frontNum = 0;
	int endNum = 0;
	int order;
public:
	string front;
	string end;
	string all;
	BWTmatrix() { front = ""; end = ""; all = ""; }
	void setOrder(int order) {
		this->order = order;
	}
	void setFront(char front) {
		this->front += front;
	}
	void setEnd(char end) {
		this->end += end;
	}
	string getFront() {
		return front;
	}
	string getEnd() {
		return end;
	}
	void setFrontNum(int frontNum) {
		this->frontNum = frontNum;
	}
	void setEndNum(int endNum) {
		this->endNum = endNum;
	}
	int getOrder() {
		return order;
	}
	int getFrontNum() {
		return frontNum;
	}
	int getEndNum() {
		return endNum;
	}
};

/*�ؽ� �ε��� �˰��򿡼� ���� Ŭ����*/
class node {
public:
	int address = -1;
	node* next = NULL;
	node() { address = -1; next = NULL; }
};


bool operator < (const BWTmatrix &a, const BWTmatrix & b) {
	return a.all < b.all;
}

void makeMygenome(string myGenome);   //myGenome�� ����
void makeRefgenome(string myGenome, string refGenome);   //������ refGenome�� ����
void makeShort(string myGenome, string* shorts);   //short reads�� ����


/*�������� ������� myGenome�� �����Ѵ�.*/
void bruteFindShorts(string refGenome, string recoverGenome, string* shorts);   


/*BWT �˰����� ����Ͽ� myGenome�� �����Ѵ�*/
void BWTrecover(string input, string recoverGenome, string* shorts);   
void makeMatrix(string input);	//BWT �˰��򿡼� ����� �迭 ����
void sortMatrix();	//makeMatrix���� ���� �迭�� ���ĺ������� ����
void checkDup();	//�չ��ڿ� �޹��ڿ��� �ߺ��Ǵ� ���ڵ鿡 ���ڸ� �ű�
void findInput(string* shorts, string recoverGenome);	//������ �迭���� short reads���� ã�´�
int findC(char c, int num);		//�� ���ڿ����� ã���� �ϴ� ������ �ε����� ã�Ƽ� ��ȯ�Ѵ�.

/*�ؽ� �ε��� �˰����� ����Ͽ� myGenome�� �����Ѵ�*/
void hashRecover(string input,string recoverGenome, string* shorts);
void makeHash(string* shorts);	//short reads���� �ؽ��Ͽ� �迭�� ����
void findHash(string input, string recoverGenome, string* shorts);	//reference Genome���� short reads���� ã�Ƽ� myGenome ����

void cmpGenome(string myGenome, string reGenome);	//�˰������ ��Ȯ���� ���ϴ� �Լ�

#define SHORTNUM 400
#define SHORTLONG 30   //5~30
#define MISMATCH 3
#define MAXDNA 2000//10000~

#define MATRIXLONG 30	//BWT ���� sorting �Ҷ� �ʿ��� ���ڿ��� ����


node* hashTable = new node[pow(4, SHORTLONG/3)];

vector<BWTmatrix> output(MAXDNA + 2);
int dupNum[4][2] = { 0, };   //A,C,G,T ��
#include <iostream>
#include <string>
#include <fstream>
#include <time.h>
#include <ctime>
#include <algorithm>
#include <vector>
#include "WELL1024a.h"

using namespace std;


/*BWT 알고리즘에서 사용될 클래스*/
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

/*해시 인덱싱 알고리즘에서 사용될 클래스*/
class node {
public:
	int address = -1;
	node* next = NULL;
	node() { address = -1; next = NULL; }
};


bool operator < (const BWTmatrix &a, const BWTmatrix & b) {
	return a.all < b.all;
}

void makeMygenome(string myGenome);   //myGenome을 생성
void makeRefgenome(string myGenome, string refGenome);   //참고할 refGenome을 생성
void makeShort(string myGenome, string* shorts);   //short reads를 생성


/*직선적인 방법으로 myGenome을 복구한다.*/
void bruteFindShorts(string refGenome, string recoverGenome, string* shorts);   


/*BWT 알고리즘을 사용하여 myGenome을 복구한다*/
void BWTrecover(string input, string recoverGenome, string* shorts);   
void makeMatrix(string input);	//BWT 알고리즘에서 사용할 배열 생성
void sortMatrix();	//makeMatrix에서 만든 배열을 알파벳순으로 정렬
void checkDup();	//앞문자와 뒷문자에서 중복되는 문자들에 숫자를 매김
void findInput(string* shorts, string recoverGenome);	//생성한 배열에서 short reads들을 찾는다
int findC(char c, int num);		//앞 문자열에서 찾고자 하는 문자의 인덱스를 찾아서 반환한다.

/*해시 인덱시 알고리즘을 사용하여 myGenome을 복구한다*/
void hashRecover(string input,string recoverGenome, string* shorts);
void makeHash(string* shorts);	//short reads들을 해싱하여 배열로 저장
void findHash(string input, string recoverGenome, string* shorts);	//reference Genome에서 short reads들을 찾아서 myGenome 복구

void cmpGenome(string myGenome, string reGenome);	//알고리즘들의 정확도를 구하는 함수

#define SHORTNUM 400
#define SHORTLONG 30   //5~30
#define MISMATCH 3
#define MAXDNA 2000//10000~

#define MATRIXLONG 30	//BWT 에서 sorting 할때 필요한 문자열의 길이


node* hashTable = new node[pow(4, SHORTLONG/3)];

vector<BWTmatrix> output(MAXDNA + 2);
int dupNum[4][2] = { 0, };   //A,C,G,T 순
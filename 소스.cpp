#define _CRT_SECURE_NO_WARNINGS
#include "���.h"

using namespace std;

int main() {
	string myGenome = "myGenome.txt", refGenome = "refGenome.txt", trivialGenome = "trivialRecoverGenome.txt", BWTGenome = "BWTRecoverGenome.txt", hashGenome="hashRecoverGenome.txt";
	string* shorts = new string[MAXDNA];
	float trivialTime, BWTTime, hashTime;

	srand(time(NULL));
	clock_t start;
	ifstream input;

	//������ ����
	makeMygenome(myGenome);
	makeRefgenome(myGenome, refGenome);
	makeShort(myGenome, shorts);

	cout << "DNA ��: " << MAXDNA << "   short read ����: " << SHORTLONG << "   short read ����: " << SHORTNUM <<"	mismatch ����: "<<MISMATCH<< endl;
	time_t curTime;
	struct tm *curr_tm;

	curTime = time(NULL);
	curr_tm = localtime(&curTime);
	cout << curr_tm->tm_hour << "�� " << curr_tm->tm_min << "�� " << curr_tm->tm_sec << "��" << endl << endl;



	start = clock();
	bruteFindShorts(refGenome, trivialGenome, shorts);
	trivialTime = (float)(clock() - start) / CLOCKS_PER_SEC;
	cout << "����ð�: " << trivialTime << endl << "-----------------------" << endl;

	start = clock();
	BWTrecover(refGenome, BWTGenome, shorts);
	BWTTime = (float)(clock() - start) / CLOCKS_PER_SEC;
	cout << "����ð�: " << BWTTime << "(sec)" << endl << "-----------------------" << endl;

	start = clock();
	hashRecover(refGenome, hashGenome, shorts);
	hashTime = (float)(clock() - start) / CLOCKS_PER_SEC;
	cout << "����ð�: " << hashTime << endl << "-----------------------" << endl;

	cout << endl << "*****<<trivial Algorithm>>*****" << endl;
	cmpGenome(myGenome, trivialGenome); cout << "   ����ð�: " << trivialTime << endl;

	cout << endl << "*****<<BWT Algorithm>>*****" << endl;
	cmpGenome(myGenome, BWTGenome);  cout << "   ����ð�: " << BWTTime << endl;

	cout << endl << "*****<<HASH INDEXING Algorithm>>*****" << endl;
	cmpGenome(myGenome, hashGenome);  cout << "   ����ð�: " << hashTime << endl;
}


void makeMygenome(string myGenome) {
	char dna[4] = { 'A','C','G','T' };
	ofstream outFile(myGenome);
	for (int i = 0; i < MAXDNA; i++) {
		outFile << dna[rand() % 4];
	}
	outFile.close();
}

void makeRefgenome(string myGenome, string refGenome) {
	char dna[4] = { 'A','C','G','T' };
	char c;
	int flag;
	ofstream outFile(refGenome);
	ifstream inFile(myGenome);

	for (int i = 0; i < MAXDNA; i++) {
		flag = rand() % 100;   //�����ڸ� �ٲ��� ������ �����ϴ� �÷���
		inFile.get(c);         //�ѱ��ڸ� �ϴ� �о� ���´�
		if (flag < 5)         //Ȯ���� 5%�� �����ߴ�. flag���� 5���� ������ ������ ���� �ٲ۴�.
			outFile << dna[rand() % 4];
		else
			outFile << c;
	}

	outFile << '$';

	outFile.close();
}

void makeShort(string myGenome, string* shorts) {
	ifstream inFile(myGenome);
	int index = 0, i = 0, j, k = 0;
	unsigned int init[32];
	srand((unsigned int)time(NULL));

	//WELL Random �� �ʱ�ȭ �ϱ� ����, C ǥ�� rand() �Լ��� �̿��Ͽ� init ���� �����մϴ�
	for (int i = 0; i<32; i++) {
		init[i] = rand() << 16 | rand();
	}

	InitWELLRNG1024a(init); // WELL Random �ʱ�ȭ

	int flag;
	char c;

	for (i = 0; i < SHORTNUM; i++) {   //shortread�� ������ �ε����� 0���� (�ؽ�Ʈ������ ���� - shortread�� ����) ���̿��� �������� �����Ѵ�.
		index =(int)((double)WELLRNG1024a()*(MAXDNA+1-SHORTLONG));
		inFile.seekg(index);
		for (int j = 0; j < SHORTLONG; j++) {
			inFile.get(c);      //�տ��� ������ �ε������� �����ϴ� ���ڿ��� shortread�� �����Ѵ�.
			shorts[i] += c;
		}
	}
	cout << "short reads �����Ϸ�.." << endl;
}


void bruteFindShorts(string refGenome, string recoverGenome, string* shorts) {
	int i, j, k, mis = 0;
	char c;
	bool match = false;
	string semiShort;
	ifstream inFile(refGenome);
	ofstream outFile(recoverGenome);

	cout << endl << "*****<< " << refGenome << "  trivial search>>*****" << endl << endl;
	for (i = 0; i < SHORTNUM; i++) {
		if (i % 10 == 0) {
			cout << (double)i / (double)SHORTNUM * 100.0 << "% �Ϸ�.." << endl;
		}
		while (!inFile.eof()) {
			k = inFile.tellg();
			mis = 0;
			semiShort = "";
			for (j = 0; j < SHORTLONG; j++) {
				inFile.get(c);
				semiShort += c;
				if (c != shorts[i][j]) {
					if (mis >= MISMATCH) {   //�ٸ� �����ڰ� MISMATCH �̻��̸� Ʋ�������� �����Ѵ�.
						match = false;
						inFile.seekg(-j, ios::cur);   //Ŀ���� j��ŭ �̵����� �ٽ� ���ڿ��� ���Ѵ�.
						break;
					}
					else mis++;
				}
				else match = true;
			}
			if (match == true) {
				outFile.seekp(k);
				outFile << shorts[i];
			}
			if (inFile.fail())break;
		}
		inFile.clear();
		inFile.seekg(0);
	}
	inFile.close();
	outFile.close();
}

void makeMatrix(string input) {
	clock_t start;
	//start = clock();
	ofstream a;
	ifstream ref;
	int i = 0, j;
	char c;
	ref.open(input);

	for (i = 0; i <= MAXDNA; i++) {
		ref.get(c);
		output[i].setFront(c);
		output[i].setOrder(i);
	}

	ref.clear();

	for (i = 0; i <= MAXDNA; i++) { //���� ����
		output[i].setEnd(output[(MAXDNA + i) % (MAXDNA + 1)].getFront()[0]);
	}

	for (i = 0; i <= MAXDNA; i++) {   //sort �ϱ� ���� �ʿ��� ���ڿ� ����
		for (j = 1; j < MATRIXLONG; j++) {
			output[i].all += output[(i + j) % (MAXDNA + 1)].end;
		}
	}

}

void sortMatrix() {   //���ĺ������� ���ڿ� ����
	sort(output.begin(), output.end());
}

void checkDup() {
	int i; char c1, c2;

	for (i = 1; i <= MAXDNA + 1; i++) {
		c1 = output[i].getFront()[0];
		c2 = output[i].getEnd()[0];
		switch (c1)
		{
		case'A':
			dupNum[0][0]++;
			output[i].setFrontNum(dupNum[0][0]);
			break;
		case'C':
			dupNum[1][0]++;
			output[i].setFrontNum(dupNum[1][0]);
			break;
		case'G':
			dupNum[2][0]++;
			output[i].setFrontNum(dupNum[2][0]);
			break;
		case'T':
			dupNum[3][0]++;
			output[i].setFrontNum(dupNum[3][0]);
			break;
		default:
			break;
		}

		switch (c2) {
		case'A':
			dupNum[0][1]++;
			output[i].setEndNum(dupNum[0][1]);
			break;
		case'C':
			dupNum[1][1]++;
			output[i].setEndNum(dupNum[1][1]);
			break;
		case'G':
			dupNum[2][1]++;
			output[i].setEndNum(dupNum[2][1]);
			break;
		case'T':
			dupNum[3][1]++;
			output[i].setEndNum(dupNum[3][1]);
			break;
		default:
			break;
		}
	}

}

void findInput(string* shorts, string recoverGenome) {
	ofstream outFile(recoverGenome);
	int sindex = 0, eindex = MAXDNA + 2, k = 3, num = MAXDNA, m, findNum = 0, mismatch = 0;
	double percent;
	char c, n;
	string tmp, repair;
	bool find = false;


	for (int i = 0; i < SHORTNUM; i++) {
		int* finds = new int[100];
		c = shorts[i][SHORTLONG - 1];

		if (i % 10 == 0) {
			cout << (double)i / (double)SHORTNUM * 100.0 << "% �Ϸ�.." << endl;
		}

		switch (c) {   //ã�ƾ��� ������ �� �� ���ڰ� ���� �˾Ƴ���.
		case'A':
			sindex = 1; eindex = dupNum[0][0] + 1;
			break;
		case'C':
			sindex = dupNum[0][0] + 1; eindex = dupNum[0][0] + dupNum[1][0] + 1;
			break;
		case'G':
			sindex = dupNum[0][0] + dupNum[1][0] + 1; eindex = dupNum[0][0] + dupNum[1][0] + dupNum[2][0] + 1;
			break;
		case'T':
			sindex = dupNum[0][0] + dupNum[1][0] + dupNum[2][0] + 1; eindex = MAXDNA + 2;
			break;
		default:
			break;
		}

		for (int j = sindex; j < eindex; j++) {      //�ش� �Ǵ� �������� �˻�
			m = j;
			mismatch = 0;
			for (k = 2; k <= SHORTLONG; k++) {
				n = shorts[i][SHORTLONG - k];   //n�� ���� ã�ƾ� �ϴ� ����.

				find = false;

				if (output[m].getEnd()[0] == n) {//������ ������ ã�ƾ��ϴ� ���ڿ� ��ġ�ϸ� ��� �˻� ����(������ ���� ���ڸ� �Ӹ����� ã�´�.)
					m = findC(output[m].getEnd()[0], output[m].getEndNum());   // m = ������ ��ġ�ϴ� �Ӹ��� �ε���																			  
					find = true;
				}
				else if (mismatch<MISMATCH&&output[m].getEnd() != "$") { //������ ������ ã�ƾ��ϴ� ���ڿ� ��ġ�ϸ� ��� �˻� ����(������ ���� ���ڸ� �Ӹ����� ã�´�.)
					m = findC(output[m].getEnd()[0], output[m].getEndNum());   // m = ������ ��ġ�ϴ� �Ӹ��� �ε���
					mismatch++;
					find = true;
				}

				if (find == false) {  break; }
			}
			if (find == true) {   //�� �������� ��� ���� 
				finds[findNum] = output[j].getOrder() - SHORTLONG + 2;
				outFile.seekp(finds[findNum] - 1);
				outFile << shorts[i];
				findNum++;
			}
		}

		sort(finds, finds + findNum);
		delete[] finds;
		findNum = 0;
	}

}

int findC(char c, int num) {
	int sindex, eindex, i;

	switch (c) {   
	case'A':
		sindex = 2;
		break;
	case'C':
		sindex = dupNum[0][0] + 2; 
		break;
	case'G':
		sindex = dupNum[0][0] + dupNum[1][0] + 2; 
		break;
	case'T':
		sindex = dupNum[0][0] + dupNum[1][0] + dupNum[2][0] + 2;
		break;
	default:
		cout << "error.. " << c << num << endl;
		return -1;
		break;
	}

	i = sindex + num - 1;

	return i;
}

void BWTrecover(string input, string recoverGenome, string* shorts) {
	//BWT ����

	cout << endl << "*****<< " << input << "  BWT recover>>*****" << endl << endl;

	cout << input << " �� ���� ������ �迭 ����.." << endl << endl;
	makeMatrix(input);

	cout << endl << "���ĺ������� ����.." << endl << endl;
	sortMatrix();

	cout << endl << "���� ��ũ üũ.." << endl << endl;
	checkDup();

	cout << endl << "MY GENOME ����.." << endl << endl;
	findInput(shorts, recoverGenome);
}

void hashRecover(string input,string recoverGenome, string* shorts) {
	cout << endl << "*****<< " << input << "  Hash Indexing search>>*****" << endl << endl;
	cout << "�ؽ� ���̺� ����..."<<endl;
	makeHash(shorts);
	cout << "MY GENOME ����.." << endl;
	findHash(input, recoverGenome, shorts);
}

void makeHash(string* shorts) {
	int i, j;
	long long int hashValue = 0;
	for (i = 0; i < SHORTNUM; i++) {
		hashValue = 0;
		for (j = 0; j < SHORTLONG / 3; j++) {	//SHORT������ 1/3 ������ �ؽ����� ��� �Ѵ�.
			switch (shorts[i][j]) {		//�ؽ� �� ���
			case'A':
				hashValue = 0 + hashValue*4;
				break;
			case'C':
				hashValue = 1 + hashValue * 4;
				break;
			case'G':
				hashValue = 2 + hashValue * 4;
				break;
			case'T':
				hashValue = 3 + hashValue * 4;
				break;
			default:
				break;
			}
		}
		if (hashTable[hashValue].address == -1)  
			hashTable[hashValue].address = i;
		else {
			node* newNode = new node;
			node* tmp = &hashTable[hashValue];
			while (tmp->next != NULL) {
				tmp = tmp->next;
			}
			newNode->address = i;
			tmp->next = newNode;
		}
	}
}

void findHash(string input, string recoverGenome, string* shorts) {
	ifstream original(input);
	ofstream recover(recoverGenome);
	long long int hashValue = 0;
	int mis = 0;
	bool flag;
	char c;

	string tmp;

	for (int i = 0; i < MAXDNA; i++) {
		hashValue = 0;
		mis = 0;
		flag = false;
		original.seekg(i, ios::beg);

		if (i % 10000 == 0) {
			cout << (double)i / (double)MAXDNA * 100.0 << "% �Ϸ�.." << endl;
		}

		for (int j = 0; j < SHORTLONG / 3; j++) {
			original.get(c);
			tmp += c;
			switch (c) {
			case'A':
				hashValue = 0 + hashValue * 4;
				break;
			case'C':
				hashValue = 1 + hashValue * 4;
				break;
			case'G':
				hashValue = 2 + hashValue * 4;
				break;
			case'T':
				hashValue = 3 + hashValue * 4;
				break;
			default:
				break;
			}
		}

		if (hashTable[hashValue].address != -1) {	//�˰� �ִ� �ؽ� �� �϶� shorts ���� ã�´�.
			original.seekg(-SHORTLONG / 3, ios::cur);
			for (int k = 0; k < SHORTLONG; k++) {	
				original.get(c);
				if (mis > MISMATCH)break;
				if (shorts[hashTable[hashValue].address][k] != c) mis++;
				else flag = true;
			}

			if (flag == true) {
				recover.seekp(i);
				recover << shorts[hashTable[hashValue].address];
			}
		}

		
	}
	
}

void cmpGenome(string myGenome, string reGenome) {
	ifstream original(myGenome);
	ifstream recover(reGenome);

	char c1, c2;
	int correct = 0;
	double percent;

	while (!original.eof() && !recover.eof()) {
		original.get(c1);
		recover.get(c2);
		if (c1 == c2)correct++;
	}

	percent = (double)correct / (double)MAXDNA * 100.0;
	cout << "correct: " << correct << "   ";
	cout << "��Ȯ��: " << percent << "%";
}
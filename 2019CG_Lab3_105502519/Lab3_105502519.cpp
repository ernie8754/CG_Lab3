#include <gl/GLUT.H>
#include <iostream>
#include <algorithm>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <string>
#include <fstream>
#include "Lab3_105502519.h"
#define PI 3.14159265f

std::fstream file;
float Ar;
bool HaveSquare, HaveTriangle;
int H, W;
float vxl, vxr, vyb, vyt;
float Ex, Ey, Ez, COIx, COIy, COIz, Tilt, Hither, Yon, Hav;
bool noback=false;

struct matrix
{
	std::vector<std::vector<float>> M;
	matrix multi(matrix M2) {
		matrix M3;
		M3.setSize(M.size(), M2.M[0].size());
		for (size_t i = 0; i < M.size(); i++) {
			for (size_t j = 0; j < M2.M[0].size(); j++) {
				M3.M[i][j] = 0;
				for (size_t k = 0; k < M[0].size(); k++) {
					M3.M[i][j] = M3.M[i][j] + M[i][k] * M2.M[k][j];
				}
			}
		}
		return M3;
	}
	matrix cross(matrix M2) {
		matrix M3;
		M3.M = {
			{M[0][1] * M2.M[0][2] - M[0][2] * M2.M[0][1], 
			-(M[0][0] * M2.M[0][2] - M[0][2] * M2.M[0][0]), 
			M[0][0] * M2.M[0][1] - M[0][1] * M2.M[0][0]}
		};
		return M3;
	}
	void PD() {
		for (size_t i = 0; i < M[0].size(); i++) {
			for (size_t j = 0; j < M.size(); j++) {
				M[j][i] /= M[M.size() - 1][i];
				
			}
			//std::cout << M[M.size() - 2][i];
		}
	}
	void general() {
		float l = sqrt(M[0][0] * M[0][0] + M[0][1] * M[0][1] + M[0][2] * M[0][2]);
		M[0][0] /= l;
		M[0][1] /= l;
		M[0][2] /= l;
	}
	void setSize(size_t row, size_t col) {
		M.resize(row);
		for (size_t i = 0; i < M.size(); i++) {
			M[i].resize(col);
		}
	}
};
matrix TM;
matrix EM;
matrix PM;
void showMI(matrix M) {
	for (size_t i = 0; i < M.M.size(); i++) {
		std::cout << "\n[ ";
		for (size_t j = 0; j < M.M[0].size(); j++) {
			std::cout << M.M[i][j] << ", ";
		}
		std::cout << "]\n";
	}
}
struct obj {
	matrix vertices;
	std::vector<std::vector<int>> lineInfo;
};
struct line {
	float x1, y1, z1, w1, x2, y2, z2, w2;
	line(float x1, float y1, float z1, float w1, float x2, float y2, float z2, float w2) {
		this->x1 = x1;
		this->y1 = y1;
		this->z1 = z1;
		this->w1 = w1;
		this->x2 = x2;
		this->y2 = y2;
		this->z2 = z2;
		this->w2 = w2;
	}
	void show() {
		if(x1 > 1 || x2 > 1 || x1<-1||x2<-1||
			y1 > 1 || y2 > 1 || y1 < -1 || y2 < -1 ||
			z1 > 1 || z2 > 1 || z1 < -1 || z2 < -1 )
		std::cout << "\n[ "<< x1 << ", " << x2 
			<< "\n" << y1 << ", " << y2
			<< "\n" << z1 << ", " << z2
			<< "\n" << w1 << ", " << w2 << "]\n";
	}
	matrix vec() {
		matrix M;
		M.M = { {x2 - x1, y2 - y1, z2 - z1} };
		return M;
	}
};
struct plane {
	std::vector<line> line;
};
std::vector<matrix> squares;
std::vector<matrix> triangles;
std::vector<obj> objs;
std::vector<line> lines;
std::vector<line> Showlines;
//matrix square;
//matrix triangle;
class Point
{
public:
	int x;
	int y;
	Point(int x, int y) {
		this->x = x;
		this->y = y;
	}
};
std::vector<Point> AllPoints;
std::vector<Point> ShowPoint;
matrix eyeMatrix(float Ex, float Ey, float Ez, float COIx, float COIy, float COIz, float Tilt) {
	matrix TiltM;
	matrix GRM;
	matrix GRM_V1;
	matrix GRM_V2;
	matrix GRM_Vz;
	matrix EyeTM;
	matrix MirrorM;
	float ZL = sqrt((COIx - Ex)*(COIx - Ex) + (COIy - Ey)*(COIy - Ey) + (COIz - Ez)*(COIz - Ez));
	GRM_Vz.M = { {(COIx-Ex) / ZL, (COIy - Ey) / ZL, (COIz - Ez) / ZL} };
	//showM(GRM_Vz);
	GRM_V1.M = { {0,1,0} };
	GRM_V1 = GRM_V1.cross(GRM_Vz);
	GRM_V1.general();
	GRM_V2 = GRM_Vz.cross(GRM_V1);
	GRM_V2.general();
	TiltM.M = { 
		{cos(Tilt*PI / 180.0f), sin(Tilt*PI / 180.0f), 0, 0},
		{-sin(Tilt*PI / 180.0f), cos(Tilt*PI / 180.0f), 0, 0},
		{0, 0, 1, 0},
		{0, 0, 0, 1}
	};
	
	GRM.M={
		{GRM_V1.M[0][0], GRM_V1.M[0][1], GRM_V1.M[0][2], 0},
		{GRM_V2.M[0][0], GRM_V2.M[0][1], GRM_V2.M[0][2], 0},
		{GRM_Vz.M[0][0], GRM_Vz.M[0][1], GRM_Vz.M[0][2], 0},
		{0, 0, 0, 1}
	};
	//showM(GRM);
	MirrorM.M = {
		{-1, 0, 0, 0},
		{0, 1, 0, 0},
		{0, 0, 1, 0},
		{0, 0, 0, 1},
	};
	EyeTM.M = {
		{1, 0, 0, -Ex},
		{0, 1, 0, -Ey},
		{0, 0, 1, -Ez},
		{0, 0, 0, 1}
	};
	/*showMI(GRM);
	showMI(MirrorM);
	showMI(TiltM);*/
	return TiltM.multi(MirrorM.multi(GRM.multi(EyeTM)));
}
matrix projectMatrix(float Ar, float H, float Y, float theta){
	matrix PM;
	PM.M= {
		{1, 0, 0, 0},
		{0, Ar, 0, 0},
		{0, 0, Y*tan(theta*PI / 180.0f) / (Y - H), H*Y*tan(theta*PI / 180.0f) / (H - Y)},
		{0, 0, tan(theta*PI / 180.0f), 0}
	};
	//showMI(PM);
	return PM;
}
void writeDot(int x, int y) {
	glColor3f(0.0f, 0.0f, 0.0f);
	glBegin(GL_POINTS);
	glVertex2f(x, y);
	glEnd();
	glFlush();
	Point temp = Point(x, y);
	AllPoints.push_back(temp);
}
void writeLine(int x1, int y1, int x2, int y2) {
	writeDot(x1, y1);
	writeDot(x2, y2);
	bool swap = false;
	bool isNegX = false;
	bool isNegY = false;
	//std::cout << "x1:" << x1 << std::endl << "y1:" << y1 << std::endl << "x2:" << x2 << std::endl << "y2:" << y2 << std::endl;
	if (abs(x2 - x1) < abs(y2 - y1)) {
		swap = true;
		//std::cout << "swap";
	}
	if (x2 < x1) {
		isNegX = true;
		//std::cout << "NegX";
	}
	if (y2 < y1) {
		isNegY = true;
		//std::cout << "NegY";
	}
	if (isNegX) {
		x2 = -x2;
		x1 = -x1;
	}
	if (isNegY) {
		y2 = -y2;
		y1 = -y1;
	}
	if (swap) {
		std::swap(x1, y1);
		std::swap(x2, y2);
	}
	int x = x1, y = y1;
	int a = abs(y2 - y1);
	int b = -abs(x2 - x1);
	int d = (2 * a + b) / 2;
	//std::cout << "a:" << a << std::endl << "b:" << b << std::endl << "d:" << d << std::endl;
	//std::cout << "----------------" << std::endl;
	while (x < x2)
	{
		//std::cout << d << std::endl;
		if (d <= 0) {
			x++;
			d += a;
		}
		else {
			x++;
			y++;
			d += (a + b);
		}
		int newX = x, newY = y;
		if (swap) {
			std::swap(newX, newY);
		}
		if (isNegX) {
			newX = -newX;
		}
		if (isNegY) {
			newY = -newY;
		}
		//std::cout << "x:" << newX << std::endl << "y:" << newY << std::endl;
		writeDot(newX, newY);
	}
	//std::cout << "----------------" << std::endl;
}
void reset() {
	TM.M = { {1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1} };
}
void clearData() {
	HaveSquare = false, HaveTriangle = false;
	squares.clear();
	triangles.clear();
}
void Translation(float x, float y, float z) {
	matrix Translate;
	Translate.M = { {1,0,0,x},{0,1,0,y},{0,0,1,z},{0,0,0,1} };
	TM = Translate.multi(TM);
}
matrix Translation(float x, float y, float z, matrix Matrix) {
	matrix Translate;
	Translate.M = { {1,0,0,x},{0,1,0,y},{0,0,1,z},{0,0,0,1} };
	Matrix = Translate.multi(Matrix);
	return Matrix;
}
void Scaling(float x, float y, float z) {
	matrix scale;
	scale.M = { {x,0,0,0},{0,y,0,0},{0,0,z,0},{0,0,0,1} };
	TM = scale.multi(TM);
}
matrix Scaling(float x, float y, float z, matrix Matrix) {
	matrix scale;
	scale.M = { {x,0,0,0},{0,y,0,0},{0,0,z,0},{0,0,0,1} };
	Matrix = scale.multi(Matrix);
	return Matrix;
}
void Rotating(float theta, int dir) {
	matrix Rotate;
	switch (dir)
	{
	case 0:
		
		Rotate.M = { {cos(theta*PI / 180.0f),-sin(theta*PI / 180.0f),0,0},{sin(theta*PI / 180.0f),cos(theta*PI / 180.0f),0,0},{0,0,1,0},{0,0,0,1} };
		TM = Rotate.multi(TM);
		break;
	case 1:
		Rotate.M = { {1,0,0,0},{0,cos(theta*PI / 180.0f),-sin(theta*PI / 180.0f),0},{0,sin(theta*PI / 180.0f),cos(theta*PI / 180.0f),0},{0,0,0,1} };
		TM = Rotate.multi(TM);
		break;
	case 2:
		Rotate.M = { {cos(theta*PI / 180.0f),0,sin(theta*PI / 180.0f),0},{0,1,0,0},{-sin(theta*PI / 180.0f),0,cos(theta*PI / 180.0f),0},{0,0,0,1} };
		TM = Rotate.multi(TM);
		break;
	default:
		break;
	}
	//matrix Rotate;
	//Rotate.M = { {cos(theta*PI / 180.0f),-sin(theta*PI / 180.0f),0},{sin(theta*PI / 180.0f),cos(theta*PI / 180.0f),0},{0,0,1} };
	//float x = TM.M[0][2];
	//float y = TM.M[1][2];
	//Translation(-x, -y);
	//TM = Rotate.multi(TM);
	//Translation(x, y);
}
matrix Rotating(float theta, int dir, matrix Matrix) {
	matrix Rotate;
	switch (dir)
	{
	case 0:
		Rotate.M = { {cos(theta*PI / 180.0f),-sin(theta*PI / 180.0f),0,0},{sin(theta*PI / 180.0f),cos(theta*PI / 180.0f),0,0},{0,0,1,0},{0,0,0,1} };
		Matrix = Rotate.multi(Matrix);
		break;
	case 1:
		Rotate.M = { {1,0,0,0},{0,cos(theta*PI / 180.0f),-sin(theta*PI / 180.0f),0},{0,sin(theta*PI / 180.0f),cos(theta*PI / 180.0f),0},{0,0,0,1} };
		Matrix = Rotate.multi(Matrix);
		break;
	case 2:
		Rotate.M = { {cos(theta*PI / 180.0f),0,sin(theta*PI / 180.0f),0},{0,1,0,0},{-sin(theta*PI / 180.0f),0,cos(theta*PI / 180.0f),0},{0,0,0,1} };
		Matrix = Rotate.multi(Matrix);
		break;
	default:
		break;
	}
	//matrix Rotate;
	//Rotate.M = { {cos(theta*PI / 180.0f),-sin(theta*PI / 180.0f),0},{sin(theta*PI / 180.0f),cos(theta*PI / 180.0f),0},{0,0,1} };
	//float x = Matrix.M[0][2];
	//float y = Matrix.M[1][2];
	////Matrix = Translation(-x, -y, Matrix);
	//Matrix = Rotate.multi(Matrix);
	////Matrix = Translation(x, y, Matrix);
	return Matrix;
}
void RemoveBack(obj Obj) {
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {

		}
	}
}
bool clipLine(line *l, float t, float m, float c[12], bool type) {
	float x = l->x1 + t * (l->x2 - l->x1);
	float y = l->y1 + t * (l->y2 - l->y1);
	float z = l->z1 + t * (l->z2 - l->z1);
	float w = l->w1 + t * (l->w2 - l->w1);
	if (type) {
		l->x1 = x;
		l->y1 = y;
		l->z1 = z;
		l->w1 = w;
		//std::cout << "jj" << w + x << "kkl";
		//l->show();
	}
	else {
			
		l->x2 = x;
		l->y2 = y;
		l->z2 = z;
		l->w2 = w;
		//l->show();
	}
	if (
			w - x >= 0 &&
			w + x >= 0 &&
			w - y >= 0 &&
			w + y >= 0 &&
			z >= 0 &&
			w - z >= 0
			)
		return true;
	return false;
}
void updateC(float c[12], line *l) {
	c[0] = l->w1 + l->x1;
	c[1] = l->w1 - l->x1;
	c[2] = l->w2 + l->x2;
	c[3] = l->w2 - l->x2;
	c[4] = l->w1 + l->y1;
	c[5] = l->w1 - l->y1;
	c[6] = l->w2 + l->y2;
	c[7] = l->w2 - l->y2;
	c[8] = l->z1;
	c[9] = l->w1 - l->z1;
	c[10] = l->z2;
	c[11] = l->w2 - l->z2;
}
bool ClipCheck(line *l) {
	float c[12];
	updateC(c, l);
	if (c[0] < 0  && c[2] < 0) {
		return false;
	}
	else if (c[1] < 0 && c[3] < 0) {
		return false;
	}
	else if (c[4] < 0 && c[6] < 0) {
		return false;
	}
	else if (c[5] < 0 && c[7] < 0) {
		return false;
	}
	else if (c[8] < 0 && c[10] < 0) {
		return false;
	}
	else if (c[9] < 0 && c[11] < 0) {
		return false;
	}
	else if (c[0] >= 0 && c[1] >= 0 && c[2] >= 0 && c[3] >= 0 &&
		c[4] >= 0 && c[5] >= 0 && c[6] >= 0 && c[7] >= 0 &&
		c[8] >= 0 && c[9] >= 0 && c[10] >= 0 && c[11] >= 0)
	{
		//l->show();
		return true;
	}
	else {
		bool clipC=false;
		if (c[0] < 0) {
			
			//l->show();
			//std::cout << "jjjj\n";
			if (clipLine(l, c[0] / (c[0] - c[2]), c[2] - c[0], c, 1)) {
				
				//l->show();
				clipC = true;
			}
			updateC(c, l);
		}
		if (c[2] < 0) {
			//l->show();
			if (clipLine(l, c[0] / (c[0] - c[2]), c[2] - c[0], c, 0)) {
				
				clipC = true;
			}
			updateC(c, l);
		}
		if (c[1] < 0) {
			if (clipLine(l, c[1] / (c[1] - c[3]), c[3] - c[1], c, 1)) {
				//l->show();
				
				clipC = true;
			}
			updateC(c, l);
		}
		if (c[3] < 0) {
			if (clipLine(l, c[1] / (c[1] - c[3]), c[3] - c[1], c, 0)) {
				//l->show();
				
				clipC = true;
			}
			updateC(c, l);
		}
		if (c[4] < 0) {
			if (clipLine(l, c[4] / (c[4] - c[6]), c[6] - c[4], c, 1)) {
				//l->show();
				
				clipC = true;
			}
			updateC(c, l);
		}
		if (c[6] < 0) {
			if (clipLine(l, c[4] / (c[4] - c[6]), c[6] - c[4], c, 0)) {
				//l->show();
				
				clipC = true;
			}
			updateC(c, l);
		}
		if (c[5] < 0) {
			if (clipLine(l, c[5] / (c[5] - c[7]), c[7] - c[5], c, 1)) {
				//l->show();
				
				clipC = true;
			}
			updateC(c, l);
		}
		if (c[7] < 0) {
			if (clipLine(l, c[5] / (c[5] - c[7]), c[7] - c[5], c, 0)) {
				//l->show();
				
				clipC = true;
			}
			updateC(c, l);
		}
		if (c[8] < 0) {
			if (clipLine(l, c[8] / (c[8] - c[10]), c[10] - c[8], c, 1)) {
				//l->show();
				
				clipC = true;
			}
			updateC(c, l);
		}
		if (c[10] < 0) {
			if (clipLine(l, c[8] / (c[8] - c[10]), c[10] - c[8], c, 0)) {
				//l->show();
				
				clipC = true;
			}
			updateC(c, l);
		}
		if (c[9] < 0) {
			if (clipLine(l, c[9] / (c[9] - c[11]), c[11] - c[9], c, 1)) {
				//l->show();
				
				clipC = true;
			}
			updateC(c, l);
		}
		if (c[11] < 0) {
			if (clipLine(l, c[9] / (c[9] - c[11]), c[11] - c[9], c, 0)) {
				//l->show();
				
				clipC = true;
			}
			updateC(c, l);
		}
		return clipC;
	}
}
void clipping() {
	for (size_t i = 0; i < objs.size(); i++) {
		//std::cout << "plane: " << objs[i].lineInfo.size() << "\n";
		Showlines.reserve(objs[i].lineInfo.size() * 3);
		matrix tempM;
		tempM = PM.multi(EM.multi(objs[i].vertices));
		tempM.PD();
		for (size_t j = 0; j < objs[i].lineInfo.size(); j++) {
			plane TempPlane;
			TempPlane.line.reserve(objs[i].lineInfo[j].size());
			//showMI(tempM);
			//std::cout <<"line: " <<objs[i].lineInfo[j].size() << "\n";
			for (size_t k = 1; k <= objs[i].lineInfo[j].size(); k++) {
				line temp = line(
					tempM.M[0][objs[i].lineInfo[j][k % objs[i].lineInfo[j].size()] - 1],
					tempM.M[1][objs[i].lineInfo[j][k % objs[i].lineInfo[j].size()] - 1],
					tempM.M[2][objs[i].lineInfo[j][k % objs[i].lineInfo[j].size()] - 1],
					tempM.M[3][objs[i].lineInfo[j][k % objs[i].lineInfo[j].size()] - 1],
					tempM.M[0][objs[i].lineInfo[j][k - 1] - 1],
					tempM.M[1][objs[i].lineInfo[j][k - 1] - 1],
					tempM.M[2][objs[i].lineInfo[j][k - 1] - 1],
					tempM.M[3][objs[i].lineInfo[j][k - 1] - 1]);
				TempPlane.line.emplace_back(temp);
				//temp.show();
				//lines.push_back(temp);
				
				//std::cout << "kk\n";

			}
			if (TempPlane.line.size() > 1) {
				matrix NV=TempPlane.line[0].vec().cross(TempPlane.line[1].vec());
				if (NV.M[0][2] < 0 || !noback) {
					for (size_t k = 0; k < TempPlane.line.size(); k++) {
						if (ClipCheck(&TempPlane.line[k])) {
							Showlines.emplace_back(TempPlane.line[k]);
						}
					}
				}
			}
			
		}
	}
}
//void Square() {
//	matrix square;
//	square.M = { {1,-1,-1,1},{1,1,-1,-1},{1,1,1,1} };
//	square = TM.multi(square);
//	squares.push_back(square);
//	HaveSquare = true;
//}
//void Triangle() {
//	matrix triangle;
//	triangle.M = { {0,-1,1},{1,-1,-1},{1,1,1} };
//	HaveTriangle = true;
//	triangle = TM.multi(triangle);
//	triangles.push_back(triangle);
//}
void ViewPort(float vxl, float vxr, float vyb, float vyt) {
	writeLine(round(vxl), round(vyb), round(vxl), round(vyt));
	//writeLine(10, 10, 550, 550);
	writeLine(round(vxr), round(vyb), round(vxr), round(vyt));
	writeLine(round(vxl), round(vyb), round(vxr), round(vyb));
	writeLine(round(vxl), round(vyt), round(vxr), round(vyt));
	//std::cout <<"\n"<< Showlines.size();
	//std::cout << "\n" << lines.size();
	for (size_t i = 0; i < Showlines.size(); i++) {
		matrix showM;
		showM.M = { 
			{Showlines[i].x1, Showlines[i].x2},
			{Showlines[i].y1, Showlines[i].y2},
			{Showlines[i].z1, Showlines[i].z2},
			{Showlines[i].w1, Showlines[i].w2}
		};
		//std::cout << "point:" << showM.M[0][0] << " " << showM.M[1][0] << std::endl;
		//std::cout << "point2:" << showM.M[0][1] << " " << showM.M[1][1] << std::endl;
		showM = Translation(1, 1, 0, showM);
		//showMI(showM);
		showM = Scaling((vxr - vxl)/(2), (vyt - vyb)/(2), 1, showM);
		showM = Translation(vxl, vyb, 0, showM);
		writeLine(round(showM.M[0][0]), round(showM.M[1][0]), round(showM.M[0][1]), round(showM.M[1][1]));
		/*std::cout << "point:" << round(showM.M[0][0]) << " " << round(showM.M[1][0]) << std::endl;
		std::cout << "point2:" << round(showM.M[0][1]) << " " << round(showM.M[1][1]) << std::endl;*/
	}
}
void Initial(void)//初始化函数 
{
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);//白色背景，前3个是RGB，最后是Alpha值，用来控制透明，1.0表示完全不透明
	glMatrixMode(GL_PROJECTION);//OpenGL按照三维方式来处理图像，所以需要一个投影变换将三维图形投影到显示器的二维空间中
	gluOrtho2D(0.0, W, 0.0, H);//指定使用正投影将一个x坐标在0~200，y坐标0~150范围内的矩形坐标区域投影到显示器窗口
	glClear(GL_COLOR_BUFFER_BIT);
}
void ShowDis() {
	std::string input;
	if (!file)
		exit(0);
	file >> input;
	std::cout << input;
	if (input == "end") {
		exit(0);
	}
	else if (input == "#") {
		std::getline(file, input);
		std::cout << input;
	}
	else if (input == "scale") {
		float x, y, z;
		file >> x >> y >> z;
		Scaling(x, y, z);
	}
	else if (input == "rotate") {
		float x, y, z;
		file >> x >> y >> z;
		Rotating(z, 0);
		Rotating(x, 1);
		Rotating(y, 2);
	}
	else if (input == "translate") {
		float x, y, z;
		file >> x >> y >> z;
		Translation(x, y, z);
	}
	else if (input == "reset") {
		reset();
	}
	/*else if (input == "square") {
		Square();
	}
	else if (input == "triangle") {
		Triangle();
	}*/
	//else if (input == "clearScreen") {
	//	//glFlush();
	//	/*for (size_t i = 0; i < ShowPoint.size(); i++) {
	//		glColor3f(1.0f, 1.0f, 1.0f);
	//		glBegin(GL_POINTS);
	//		glVertex2d(ShowPoint[i].x, ShowPoint[i].y);
	//		glEnd();
	//		glFlush();
	//	}*/
		/*glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT);
		glFlush();
		ShowPoint.clear();
		AllPoints.clear();*/
		//}
		/*else if (input == "clearData") {
			clearData();
		}*/
		/*else if (input == "view") {
			float wxl, wxr, wyb, wyt, vxl, vxr, vyb, vyt;
			file >> wxl >> wxr >> wyb >> wyt >> vxl >> vxr >> vyb >> vyt;
			ViewPort(wxl, wxr, wyb, wyt, vxl, vxr, vyb, vyt);
		}*/
	else if (input == "nobackfaces") {
		noback = true;
	}
	else if (input == "object") {
		std::string fname;
		file >> fname;
		//fname = 
		//std::cout << fname[0];
		std::fstream asc;
		asc.open(fname, std::ios::in);
		if (!asc)
			return;
		int num1, num2;
		asc >> num1 >> num2;
		obj temp;
		temp.vertices.setSize(4, num1);
		//std::cout << num1<< num2;
		for (int i = 0; i < num1; i++) {
			//std::cout << i << "lll" << "\n";
			asc >> temp.vertices.M[0][i] >> temp.vertices.M[1][i] >> temp.vertices.M[2][i];
			temp.vertices.M[3][i] = 1;
		}
		//std::cout << fname;
		temp.lineInfo.resize(num2);
		//std::cout << fname;
		for (int i = 0; i < num2; i++) {
			int Pnum;
			asc >> Pnum;
			temp.lineInfo[i].resize(Pnum);
			for (int j = 0; j < Pnum; j++) {
				int PID;
				asc >> PID;
				temp.lineInfo[i][j] = PID;
				//std::cout << i << "1" << "\n";
			}
			//std::cout << i << "22" << "\n";
		}
		//std::cout << "22" << "\n";
		temp.vertices = TM.multi(temp.vertices);
		//temp.vertices = PM.multi(temp.vertices);
		objs.push_back(temp);
		//std::cout << "3" << "\n";
		asc.close();
	}
	else if (input == "observer") {
		
		file >> Ex >> Ey >> Ez >>
			COIx >> COIy >> COIz >>
			Tilt >>
			Hither >> Yon >> Hav;
		//std::cout << "\n" << Hav;
		EM = eyeMatrix(Ex, Ey, Ez, COIx, COIy, COIz, Tilt);
		//showM(PM);
	}
	else if (input == "viewport") {
		file >> vxl >> vxr >> vyb >> vyt;
		PM = projectMatrix((vxr - vxl) / (vyt - vyb), Hither, Yon, Hav);
		vxl = W * (vxl + 1) / 2;
		//std::cout << "\n" << vxl;
		vxr = W * (vxr + 1) / 2;
		//std::cout <<"\n"<< vxr;
		vyb = H * (vyb + 1) / 2;
		//std::cout << "\n" << vyb;
		vyt = H * (vyt + 1) / 2;
		//std::cout << "\n" << vyt;
		
	}
	else if (input == "display") {
		lines.clear();
		Showlines.clear();
		glFlush();
		glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT);
		glFlush();
		ShowPoint.clear();
		AllPoints.clear();
		clipping();
		ViewPort(vxl, vxr, vyb, vyt);
		//std::cout << lines.size();
		system("pause");
	}
	std::cout << std::endl;
}
void myDisplay(void)//显示回调函数
{
	while (true) {
		std::string input;
		if (!file)
			exit(0);
		file >> input;
		std::cout << input;
		if (input == "end") {
			exit(0);
		}
		else if (input == "#") {
			std::getline(file, input);
			std::cout << input;
		}
		else if (input == "scale") {
			float x, y, z;
			file >> x >> y >> z;
			Scaling(x, y, z);
		}
		else if (input == "rotate") {
			float x, y, z;
			file >> x >> y >> z;
			Rotating(z, 0);
			Rotating(x, 1);
			Rotating(y, 2);
		}
		else if (input == "translate") {
			float x, y, z;
			file >> x >> y >> z;
			Translation(x, y, z);
		}
		else if (input == "reset") {
			reset();
		}
		/*else if (input == "square") {
			Square();
		}
		else if (input == "triangle") {
			Triangle();
		}*/
		//else if (input == "clearScreen") {
		//	//glFlush();
		//	/*for (size_t i = 0; i < ShowPoint.size(); i++) {
		//		glColor3f(1.0f, 1.0f, 1.0f);
		//		glBegin(GL_POINTS);
		//		glVertex2d(ShowPoint[i].x, ShowPoint[i].y);
		//		glEnd();
		//		glFlush();
		//	}*/
			/*glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
			glClear(GL_COLOR_BUFFER_BIT);
			glFlush();
			ShowPoint.clear();
			AllPoints.clear();*/
			//}
			/*else if (input == "clearData") {
				clearData();
			}*/
			/*else if (input == "view") {
				float wxl, wxr, wyb, wyt, vxl, vxr, vyb, vyt;
				file >> wxl >> wxr >> wyb >> wyt >> vxl >> vxr >> vyb >> vyt;
				ViewPort(wxl, wxr, wyb, wyt, vxl, vxr, vyb, vyt);
			}*/
		else if (input == "nobackfaces") {
			noback = true;
		}
		else if (input == "object") {
			std::string fname;
			file >> fname;
			//fname = 
			//std::cout << fname[0];
			std::fstream asc;
			asc.open(fname, std::ios::in);
			if (!asc)
				return;
			int num1, num2;
			asc >> num1 >> num2;
			obj temp;
			temp.vertices.setSize(4, num1);
			//std::cout << num1<< num2;
			for (int i = 0; i < num1; i++) {
				//std::cout << i << "lll" << "\n";
				asc >> temp.vertices.M[0][i] >> temp.vertices.M[1][i] >> temp.vertices.M[2][i];
				temp.vertices.M[3][i] = 1;
			}
			//std::cout << fname;
			temp.lineInfo.resize(num2);
			//std::cout << fname;
			for (int i = 0; i < num2; i++) {
				int Pnum;
				asc >> Pnum;
				temp.lineInfo[i].resize(Pnum);
				for (int j = 0; j < Pnum; j++) {
					int PID;
					asc >> PID;
					temp.lineInfo[i][j] = PID;
					//std::cout << i << "1" << "\n";
				}
				//std::cout << i << "22" << "\n";
			}
			//std::cout << "22" << "\n";
			temp.vertices = TM.multi(temp.vertices);
			//temp.vertices = PM.multi(temp.vertices);
			objs.push_back(temp);
			//std::cout << "3" << "\n";
			asc.close();
		}
		else if (input == "observer") {
			float Ex, Ey, Ez, COIx, COIy, COIz, Tilt, Hither, Yon, Hav;
			file >> Ex >> Ey >> Ez >>
				COIx >> COIy >> COIz >>
				Tilt >>
				Hither >> Yon >> Hav;
			//std::cout << "\n" << Hav;
			EM = eyeMatrix(Ex, Ey, Ez, COIx, COIy, COIz, Tilt);
			PM = projectMatrix(Ar, Hither, Yon, Hav);
			//showM(PM);
		}
		else if (input == "viewport") {
			file >> vxl >> vxr >> vyb >> vyt;
			vxl = W * (vxl + 1) / 2;
			//std::cout << "\n" << vxl;
			vxr = W * (vxr + 1) / 2;
			//std::cout <<"\n"<< vxr;
			vyb = H * (vyb + 1) / 2;
			//std::cout << "\n" << vyb;
			vyt = H * (vyt + 1) / 2;
			//std::cout << "\n" << vyt;
		}
		else if (input == "display") {
			lines.clear();
			Showlines.clear();
			glFlush();
			glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
			glClear(GL_COLOR_BUFFER_BIT);
			glFlush();
			ShowPoint.clear();
			AllPoints.clear();
			clipping();
			ViewPort(vxl, vxr, vyb, vyt);
			//std::cout << lines.size();
			system("pause");
		}
		std::cout << std::endl;
	}
}
int main(int argc, char * argv[])//这是使用glut库函数进行窗口管理
{
	system("pause");
	reset();
	clearData();
	if (argv[1])
		file.open(argv[1], std::ios::in);
	glutInit(&argc, argv);//使用glut库需要进行初始化
	glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE);//设定窗口显示模式，颜色模型和缓存，这里是RGB颜色模型和单缓存
	glutInitWindowPosition(200, 200);//设定窗口的初始位置，屏幕左上角为原点，单位为像素
	file >> W >> H;
	//W *= 2;
	//H *= 2;
	Ar = W * 1.0f / H;
	glutInitWindowSize(W, H);//设定窗口的大小
	glutCreateWindow("My Progect 2");//创建一个窗口，参数是窗口标题名
	//glutDisplayFunc(&myDisplay);//将myDisplay指定为当前窗口的显示内容函数
	glutIdleFunc(&ShowDis);
	Initial();
	//glutMouseFunc(&mouseButton);
	//glutKeyboardFunc(&keyboard);
	//glutMotionFunc(&mouseDrag);
	glutMainLoop();//使窗口框架运行起来，使显示回调函数开始工作


	return 0;
}
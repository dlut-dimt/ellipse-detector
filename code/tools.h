/*
This code is intended for academic use only.
You are free to use and modify the code, at your own risk.

If you use this code, or find it useful, please refer to the paper:


The comments in the code refer to the abovementioned paper.
If you need further details about the code or the algorithm, please contact me at:

lianbosong@foxmail.com

last update: 
*/


#pragma once
#include "stdafx.h"
#include "CNEllipseDetector.h"

#define V2SP Point2f p3,Point2f p2,Point2f p1,Point2f p4,Point2f p5,Point2f p6

void MultiImage_OneWin(const std::string& MultiShow_WinName, const vector<Mat>& SrcImg_V, CvSize SubPlot, CvSize ImgMax_Size);
Point2f lineCrossPoint(Point2f l1p1,Point2f l1p2,Point2f l2p1,Point2f l2p2 );
void point2Mat(Point2f p1,Point2f p2,float mat[2][2]);
float value4SixPoints( V2SP );
void PyrDown(string picName);
Mat matResize(Mat src,double scale);
void SaveEllipses(const string& fileName, const vector<Ellipse>& ellipses);
// 14pr
void SaveEllipses(const string& workingDir, const string& imgName, const vector<Ellipse>& ellipses /*, const vector<double>& times*/);
void LoadGT(vector<Ellipse>& gt, const string& sGtFileName, bool bIsAngleInRadians = true);
// useless 
bool LoadTest(vector<Ellipse>& ellipses, const string& sTestFileName, vector<double>& times, bool bIsAngleInRadians = true);

bool TestOverlap(const Mat1b& gt, const Mat1b& test, float th);
int Count(const vector<bool> v);
float Evaluate(const vector<Ellipse>& ellGT, const vector<Ellipse>& ellTest, const float th_score, const Mat3b& img);
// end 14pr

//saalt
void salt(cv::Mat& image, int n);
//salt


//绘制选出的弧段
/**
输入：一个象限的弧段容器
输出：显示弧段，一个象限的一个颜色。颜色随机
*/
void showEdge(vector<vector<Point>> points_,Mat& picture);
//绘制选出的弧段

//file operation
int writeFile(string fileName_cpp,vector<string> vsContent);
int readFile(string fileName_cpp);
int readFileByChar(string fileName_split);
void Trim(string &str);

/******分离特定格式的数据******/
//C++中没有Split()这个方法，需要自定义函数分离数据，而C#和Java中有这个方法
vector<string> getStr(string str);
/******分离特定格式的数据******/
/**
* path:目录
* files：用于保存文件名的vector
* r：是否需要遍历子目录
*/
void listDir(string real_dir,vector<string>& files,bool r = false);
//file operation

//90 82 49 255
//#define V2SP Point2f p3,Point2f p2,Point2f p1,Point2f p6,Point2f p5,Point2f p4
//113 77 69 335
//#define V2SP Point2f p1,Point2f p5,Point2f p3,Point2f p6,Point2f p2,Point2f p4
//116 77 74 337
//#define V2SP Point2f p1,Point2f p5,Point2f p3,Point2f p4,Point2f p2,Point2f p6
//123 82 84 316
//#define V2SP Point2f p3,Point2f p2,Point2f p1,Point2f p4,Point2f p5,Point2f p6

//93 80 48 246
//#define V2SP Point2f p1,Point2f p2,Point2f p3,Point2f p4,Point2f p5,Point2f p6
//92 75 46 240
//#define V2SP Point2f p1,Point2f p2,Point2f p5,Point2f p6,Point2f p3,Point2f p4
//99 82 61 306
//#define V2SP Point2f p1,Point2f p2,Point2f p3,Point2f p6,Point2f p4,Point2f p5
//99 82 61 306
//#define V2SP Point2f p1,Point2f p2,Point2f p4,Point2f p5,Point2f p3,Point2f p6
//100 82 61 306
//#define V2SP Point2f p1,Point2f p2,Point2f p3,Point2f p5,Point2f p4,Point2f p6
//99 82 61 306
//#define V2SP Point2f p1,Point2f p2,Point2f p4,Point2f p6,Point2f p3,Point2f p5

//100 82 60 301
//#define V2SP Point2f p1,Point2f p3,Point2f p2,Point2f p4,Point2f p5,Point2f p6
//116 77 76 341
//#define V2SP Point2f p1,Point2f p3,Point2f p5,Point2f p6,Point2f p2,Point2f p4
//98 82 59 300
//#define V2SP Point2f p1,Point2f p3,Point2f p2,Point2f p6,Point2f p4,Point2f p5
//100 87 55 292
//#define V2SP Point2f p1,Point2f p3,Point2f p4,Point2f p5,Point2f p2,Point2f p6
//125 82 59 300
//#define V2SP Point2f p1,Point2f p3,Point2f p2,Point2f p5,Point2f p4,Point2f p6
//102 87 55 292
//#define V2SP Point2f p1,Point2f p3,Point2f p4,Point2f p6,Point2f p2,Point2f p5

//100 82 60 301
//#define V2SP Point2f p1,Point2f p4,Point2f p2,Point2f p3,Point2f p5,Point2f p6
//116 77 76 341
//#define V2SP Point2f p1,Point2f p4,Point2f p5,Point2f p6,Point2f p2,Point2f p3
//98 82 59 300
//#define V2SP Point2f p1,Point2f p4,Point2f p2,Point2f p6,Point2f p3,Point2f p5
//100 87 55 292
//#define V2SP Point2f p1,Point2f p4,Point2f p3,Point2f p5,Point2f p2,Point2f p6
//125 82 59 300
//#define V2SP Point2f p1,Point2f p4,Point2f p2,Point2f p5,Point2f p3,Point2f p6
//102 87 55 292
//#define V2SP Point2f p1,Point2f p4,Point2f p3,Point2f p6,Point2f p2,Point2f p5

//98 87 53 291
//#define V2SP Point2f p1,Point2f p5,Point2f p2,Point2f p3,Point2f p4,Point2f p6
//139 77 69 335
//#define V2SP Point2f p1,Point2f p5,Point2f p4,Point2f p6,Point2f p2,Point2f p3
//98 82 61 302
//#define V2SP Point2f p1,Point2f p5,Point2f p2,Point2f p6,Point2f p3,Point2f p4
//100 77 74 337
//#define V2SP Point2f p1,Point2f p5,Point2f p3,Point2f p4,Point2f p2,Point2f p6
//125 87 53 291
//#define V2SP Point2f p1,Point2f p5,Point2f p2,Point2f p4,Point2f p3,Point2f p6
//102 77 69 335
//#define V2SP Point2f p1,Point2f p5,Point2f p3,Point2f p6,Point2f p2,Point2f p4

//98 87 53 291
//#define V2SP Point2f p1,Point2f p6,Point2f p2,Point2f p3,Point2f p4,Point2f p5
//139 77 69 335
//#define V2SP Point2f p1,Point2f p6,Point2f p4,Point2f p5,Point2f p2,Point2f p3
//98 82 61 302
//#define V2SP Point2f p1,Point2f p6,Point2f p2,Point2f p5,Point2f p3,Point2f p4
//100 77 74 337
//#define V2SP Point2f p1,Point2f p6,Point2f p3,Point2f p4,Point2f p2,Point2f p5
//125 87 53 291
//#define V2SP Point2f p1,Point2f p6,Point2f p2,Point2f p4,Point2f p3,Point2f p5
//102 77 69 335
//#define V2SP Point2f p1,Point2f p6,Point2f p3,Point2f p5,Point2f p2,Point2f p4

/*
This code is intended for academic use only.
You are free to use and modify the code, at your own risk.

If you use this code, or find it useful, please refer to the paper:


The comments in the code refer to the abovementioned paper.
If you need further details about the code or the algorithm, please contact me at:

lianbosong@foxmail.com

last update: 
*/

#include "tools.h"

Point2f lineCrossPoint(Point2f l1p1,Point2f l1p2,Point2f l2p1,Point2f l2p2 )
{
	Point2f crossPoint;
	float k1,k2,b1,b2;
	if (l1p1.x==l1p2.x&&l2p1.x==l2p2.x){
		crossPoint=Point2f(0,0);//无效点
		return crossPoint;
	}
	if (l1p1.x==l1p2.x)
	{
		crossPoint.x=l1p1.x;
		k2=(l2p2.y-l2p1.y)/(l2p2.x-l2p1.x);
		b2=l2p1.y-k2*l2p1.x;
		crossPoint.y=k2*crossPoint.x+b2;
		return crossPoint;
	}
	if (l2p1.x==l2p2.x)
	{
		crossPoint.x=l2p1.x;
		k2=(l1p2.y-l1p1.y)/(l1p2.x-l1p1.x);
		b2=l1p1.y-k2*l1p1.x;
		crossPoint.y=k2*crossPoint.x+b2;
		return crossPoint;
	}

	k1=(l1p2.y-l1p1.y)/(l1p2.x-l1p1.x);
	k2=(l2p2.y-l2p1.y)/(l2p2.x-l2p1.x);
	b1=l1p1.y-k1*l1p1.x;
	b2=l2p1.y-k2*l2p1.x;
	if (k1==k2)
	{
		crossPoint=Point2f(0,0);//无效点
	}
	else
	{
		crossPoint.x=(b2-b1)/(k1-k2);
		crossPoint.y=k1*crossPoint.x+b1;
	}
	return crossPoint;
}

void point2Mat(Point2f p1,Point2f p2,float mat[2][2])
{
	mat[0][0]=p1.x;
	mat[0][1]=p1.y;
	mat[1][0]=p2.x;
	mat[1][1]=p2.y;
}
float value4SixPoints( V2SP )
{
	float result=1;
	Mat A,B,C;
	float matB[2][2],matC[2][2];
	Point2f v,w,u;
	v=lineCrossPoint(p1,p2,p3,p4);
	w=lineCrossPoint(p5,p6,p3,p4);
	u=lineCrossPoint(p5,p6,p1,p2);

	point2Mat(u,v,matB);
	point2Mat(p1,p2,matC);
	B=Mat(2,2,CV_32F,matB);
	C=Mat(2,2,CV_32F,matC);
	A=C*B.inv();
	/*cout<<"u:\t"<<u<<endl;
	cout<<"v:\t"<<v<<endl;
	cout<<"B:\t"<<B<<endl;
	cout<<A<<endl;*/
	result*=A.at<float>(0,0)*A.at<float>(1,0)/(A.at<float>(0,1)*A.at<float>(1,1));

	point2Mat(p3,p4,matC);
	point2Mat(v,w,matB);
	B=Mat(2,2,CV_32F,matB);
	C=Mat(2,2,CV_32F,matC);
	A=C*B.inv();
	result*=A.at<float>(0,0)*A.at<float>(1,0)/(A.at<float>(0,1)*A.at<float>(1,1));

	point2Mat(p5,p6,matC);
	point2Mat(w,u,matB);
	B=Mat(2,2,CV_32F,matB);
	C=Mat(2,2,CV_32F,matC);
	A=C*B.inv();
	result*=A.at<float>(0,0)*A.at<float>(1,0)/(A.at<float>(0,1)*A.at<float>(1,1));
	return result;
}

void MultiImage_OneWin(const std::string& MultiShow_WinName, const vector<Mat>& SrcImg_V, CvSize SubPlot, CvSize ImgMax_Size)  
{  
	//Reference : http://blog.csdn.net/yangyangyang20092010/article/details/21740373  

	//************* Usage *************//  
	//vector<Mat> imgs(4);  
	//imgs[0] = imread("F:\\SA2014.jpg");  
	//imgs[1] = imread("F:\\SA2014.jpg");  
	//imgs[2] = imread("F:\\SA2014.jpg");  
	//imgs[3] = imread("F:\\SA2014.jpg");  
	//MultiImage_OneWin("T", imgs, cvSize(2, 2), cvSize(400, 280));  

	//Window's image  
	Mat Disp_Img;  
	//Width of source image  
	CvSize Img_OrigSize = cvSize(SrcImg_V[0].cols, SrcImg_V[0].rows);  
	//******************** Set the width for displayed image ********************//  
	//Width vs height ratio of source image  
	float WH_Ratio_Orig = Img_OrigSize.width/(float)Img_OrigSize.height;  
	CvSize ImgDisp_Size = cvSize(100, 100);  
	if(Img_OrigSize.width > ImgMax_Size.width)  
		ImgDisp_Size = cvSize(ImgMax_Size.width, (int)(ImgMax_Size.width/WH_Ratio_Orig));  
	else if(Img_OrigSize.height > ImgMax_Size.height)  
		ImgDisp_Size = cvSize((int)(ImgMax_Size.height*WH_Ratio_Orig), ImgMax_Size.height);  
	else  
		ImgDisp_Size = cvSize(Img_OrigSize.width, Img_OrigSize.height);  
	//******************** Check Image numbers with Subplot layout ********************//  
	int Img_Num = (int)SrcImg_V.size();  
	if(Img_Num > SubPlot.width * SubPlot.height)  
	{  
		cout<<"Your SubPlot Setting is too small !"<<endl;  
		exit(0);  
	}  
	//******************** Blank setting ********************//  
	CvSize DispBlank_Edge = cvSize(80, 60);  
	CvSize DispBlank_Gap  = cvSize(15, 15);  
	//******************** Size for Window ********************//  
	Disp_Img.create(Size(ImgDisp_Size.width*SubPlot.width + DispBlank_Edge.width + (SubPlot.width - 1)*DispBlank_Gap.width,   
		ImgDisp_Size.height*SubPlot.height + DispBlank_Edge.height + (SubPlot.height - 1)*DispBlank_Gap.height), CV_8UC3);  
	Disp_Img.setTo(0);//Background  
	//Left top position for each image  
	int EdgeBlank_X = (Disp_Img.cols - (ImgDisp_Size.width*SubPlot.width + (SubPlot.width - 1)*DispBlank_Gap.width))/2;  
	int EdgeBlank_Y = (Disp_Img.rows - (ImgDisp_Size.height*SubPlot.height + (SubPlot.height - 1)*DispBlank_Gap.height))/2;  
	CvPoint LT_BasePos = cvPoint(EdgeBlank_X, EdgeBlank_Y);  
	CvPoint LT_Pos = LT_BasePos;  

	//Display all images  
	for (int i=0; i < Img_Num; i++)  
	{  
		//Obtain the left top position  
		if ((i%SubPlot.width == 0) && (LT_Pos.x != LT_BasePos.x))  
		{  
			LT_Pos.x = LT_BasePos.x;  
			LT_Pos.y += (DispBlank_Gap.height + ImgDisp_Size.height);  
		}  
		//Writting each to Window's Image  
		Mat imgROI = Disp_Img(Rect(LT_Pos.x, LT_Pos.y, ImgDisp_Size.width, ImgDisp_Size.height));  
		resize(SrcImg_V[i], imgROI, Size(ImgDisp_Size.width, ImgDisp_Size.height));  

		LT_Pos.x += (DispBlank_Gap.width + ImgDisp_Size.width);  
	}  

	//Get the screen size of computer  #include<windows.h> #include <WinUser.h>
	int Scree_W = 1366;//GetSystemMetrics(SM_CXSCREEN);  
	int Scree_H = 768;//GetSystemMetrics(SM_CYSCREEN);  
	//cout<<Scree_W<<"\t"<<Scree_H<<endl;
	cvNamedWindow(MultiShow_WinName.c_str(), CV_WINDOW_NORMAL);
	cvMoveWindow(MultiShow_WinName.c_str(),(Scree_W - Disp_Img.cols)/2 ,(Scree_H - Disp_Img.rows)/2);//Centralize the window  
	IplImage tmp_Disp_Img= IplImage(Disp_Img);
	cvShowImage(MultiShow_WinName.c_str(), &tmp_Disp_Img);  
	cvWaitKey(0);  
	cvDestroyWindow(MultiShow_WinName.c_str());  
}  


void PyrDown(string picName)
{
	Mat img1=imread(picName);
	Mat img2;
	Size sz;
	//金字塔向下或者向上采样操作 ，基本不改变图像长宽比率
	//pyrDown(img1,img2,sz,BORDER_DEFAULT);
	pyrUp(img1,img2,sz,BORDER_DEFAULT);
	pyrUp(img2,img2,sz,BORDER_DEFAULT);
	namedWindow("WindowOrg");
	namedWindow("WindowNew");
	imshow("WindowOrg",img1);
	imshow("WindowNew",img2);

	waitKey(10000);
}
Mat matResize(Mat src,double scale){
	Mat img2;
	bool showtimeandpic=false;
	if(!showtimeandpic){
		Size dsize = Size(int(src.cols*scale),int(src.rows*scale));
		img2 = Mat(dsize,CV_32S);
		resize(src, img2,dsize,CV_INTER_CUBIC);
	}
	else{
		clock_t start_time=clock();
		{
			Size dsize = Size(int(src.cols*scale),int(src.rows*scale));
			img2 = Mat(dsize,CV_32S);
			resize(src, img2,dsize,CV_INTER_CUBIC);
		}
		clock_t end_time=clock();
		cout<< "Running time is: "<<static_cast<double>(end_time-start_time)/CLOCKS_PER_SEC*1000<<"ms"<<endl;//输出运行时间

		//CV_INTER_NN - 最近邻差值,
		//CV_INTER_LINEAR -  双线性差值 (缺省使用) 
		//CV_INTER_AREA -  使用象素关系重采样。当图像缩小时候，该方法
		//可以避免波纹出现。当图像放大时，类似于  CV_INTER_NN  方法.. 
		//CV_INTER_CUBIC -  立方差值. 
		namedWindow("WindowOrg",CV_WINDOW_AUTOSIZE);
		namedWindow("WindowNew",CV_WINDOW_AUTOSIZE);
		imshow("WindowOrg",src);
		imshow("WindowNew",img2);

		waitKey(1000);
	}
	return img2;
}

//绘制选出的弧段
void showEdge(vector<vector<Point>> points_,Mat& picture)
{
	srand( (unsigned)time( NULL ));
	int radius=1;
	Point center;
	
	int sEdge=points_.size();
	Point prev_point;
	Point current_point;
	for (int iEdge=0;iEdge<sEdge;iEdge++){
		int r=rand()%256;
		int g=rand()%256;
		int b=rand()%256;
		Scalar color=Scalar(b,g,r);
		vector<Point> Edge=points_.at(iEdge);
		int sPoints=Edge.size();
		for(int iPoint=0;iPoint<sPoints-1;iPoint++){
			center=Edge.at(iPoint);
			//参数为：承载的图像、圆心、半径、颜色、粗细、线型  
			//circle(picture,center,radius,color); 
			prev_point=Edge.at(iPoint);
			current_point=Edge.at(iPoint+1);
			//Mat to IplImage to cvArr
			IplImage ipl_img = picture;
			cvLine(&ipl_img, prev_point, current_point, color, 1, CV_AA);
		}
	}
}
// file operation
int writeFile(string fileName_cpp,vector<string> vsContent){
	string line="";
	vector<string> data;
	vector<string> data_split;
	ofstream out(fileName_cpp);
	if(!out)
	{
		cout<<"读写文件失败"<<endl;
		return -1;
	}
	for(vector<string>::iterator i=vsContent.begin();i<vsContent.end();i++){
		out<<*i<<endl;
	}
	out.close();
	return 1;
}

int readFile(string fileName_cpp){
	string line="";
	vector<string> data;
	ifstream in(fileName_cpp);
	if(!in)
	{
		cout<<"读写文件失败"<<endl;
		return -1;
	}
	while(getline(in,line))
	{
		data.push_back(line);     //读取文件每一行数据，并放到“容器”里面
	}
	in.close();
	/******遍历data里面的内容******/
	for(unsigned int i=0;i<data.size();i++)
	{
		cout<<data.at(i)<<endl;
	}
	return 0;
	/******遍历data里面的内容******/
}
int readFileByChar(string fileName_split){

	string line="";
	vector<string> data;
	vector<string> data_split;
	ifstream in_split(fileName_split);
	if(!in_split)
	{
		cout<<"读写文件失败"<<endl;
		return -1;
	}
	while(getline(in_split,line))
	{
		data_split.push_back(line);     //读取文件每一行数据，并放到“容器”里面
	}
	in_split.close();
	/******读文件******/
	/******提取split.txt文件里面的数据******/

	/******遍历data_split里面的内容(数据分离)******/
	for(unsigned int i=0;i<data_split.size();i++)
	{
		cout<<"--------------------"<<endl;
		for(unsigned int j=0;j<getStr(data_split.at(i)).size();j++)
		{
			cout<<getStr(data_split.at(i)).at(j)<<endl;
		}
	}
	/******遍历data_split里面的内容(数据分离)******/
	return 0;
}
void Trim(string &str)
{
	int s=str.find_first_not_of(" \t\n");
	int e=str.find_last_not_of(" \t\n");
	str=str.substr(s,e-s+1);
}
/******分离特定格式的数据******/
//C++中没有Split()这个方法，需要自定义函数分离数据，而C#和Java中有这个方法
vector<string> getStr(string str)
{
	int j=0;
	string a[100];
	vector<string> v_a;
	//Split()
	for(unsigned int i=0;i<str.size();i++)
	{
		if((str[i]!=',')&&str[i]!='\0')
		{
			a[j]+=str[i];
		}
		else j++;
	}

	for(int k=0;k<j+1;k++)
	{
		v_a.push_back(a[k]);
	}
	return v_a;
}
/******分离特定格式的数据******/

/**
* path:目录
* files：用于保存文件名的vector
* r：是否需要遍历子目录
*/
void listDir(string real_dir,vector<string>& files,bool r){
	DIR *pDir;
	struct dirent *ent;
	string childpath;
	string absolutepath;
	pDir = opendir(real_dir.c_str());
	while ((ent = readdir(pDir)) != NULL){
		if (ent->d_type & DT_DIR){
			if (strcmp(ent->d_name, ".") == 0 || strcmp(ent->d_name, "..") == 0){
				continue;
			}
			if(r){ //如果需要遍历子目录
				childpath=real_dir+ent->d_name;
				listDir(childpath,files);
			}
		}
		else{
			absolutepath= real_dir+ent->d_name;
			files.push_back(ent->d_name);//文件名
		}
	}
	sort(files.begin(),files.end());//排序
}
void SaveEllipses(const string& fileName, const vector<Ellipse>& ellipses){
	unsigned n = ellipses.size();
	vector<string> resultString;
	stringstream resultsitem;
	// Save number of ellipses
	resultsitem << n ;
	resultString.push_back(resultsitem.str());
	// Save ellipses
	for (unsigned i = 0; i<n; ++i)
	{
		const Ellipse& e = ellipses[i];
		resultsitem.str("");
		resultsitem << e._xc << "\t" << e._yc << "\t" 
			<< e._a << "\t" << e._b << "\t" 
			<< e._rad << "\t" << e._score;
		resultString.push_back(resultsitem.str());
	}
	writeFile(fileName,resultString);
	for (int i=0;i<resultString.size();i++){
		cout<<resultString[i]<<endl;
	}
}

// 14pr
// Should be checked
void SaveEllipses(const string& workingDir, const string& imgName, const vector<Ellipse>& ellipses /*, const vector<double>& times*/) 
{
	string path(workingDir + "/" + imgName + ".txt");
	ofstream out(path, ofstream::out | ofstream::trunc);
	if (!out.good())
	{
		cout << "Error saving: " << path << endl;
		return;
	}
	unsigned n = ellipses.size();

	// Save number of ellipses
	out << n << "\n";
	// Save ellipses
	for (unsigned i = 0; i<n; ++i)
	{
		const Ellipse& e = ellipses[i];
		out << e._xc << "\t" << e._yc << "\t" << e._a << "\t" << e._b << "\t" << e._rad << "\t" << e._score << "\n";
	}
	out.close();
}

void LoadGT(vector<Ellipse>& gt, const string& sGtFileName, bool bIsAngleInRadians)
{
	ifstream in(sGtFileName);
	if (!in.good())
	{
		cout << "Error opening: " << sGtFileName << endl;
		return;
	}

	unsigned n;
	in >> n;

	gt.clear();
	gt.reserve(n);

	while (in.good() && n--)
	{
		Ellipse e;
		in >> e._xc >> e._yc >> e._a >> e._b >> e._rad;

		if (!bIsAngleInRadians)
		{
			// convert to radians
			e._rad = float(e._rad * CV_PI / 180.0);
		}

		if (e._a < e._b)
		{
			float temp = e._a;
			e._a = e._b;
			e._b = temp;

			e._rad = e._rad + float(0.5*CV_PI);
		}

		e._rad = fmod(float(e._rad + 2.f*CV_PI), float(CV_PI));
		e._score = 1.f;
		gt.push_back(e);
	}
	in.close();
}

// Should be checked
bool LoadTest(vector<Ellipse>& ellipses, const string& sTestFileName, vector<double>& times, bool bIsAngleInRadians)
{
	ifstream in(sTestFileName);
	if (!in.good())
	{
		cout << "Error opening: " << sTestFileName << endl;
		return false;
	}

	times.resize(6);
	in >> times[0] >> times[1] >> times[2] >> times[3] >> times[4] >> times[5];

	unsigned n;
	in >> n;

	ellipses.clear();

	if (n == 0) return true;

	ellipses.reserve(n);

	while (in.good() && n--)
	{
		Ellipse e;
		in >> e._xc >> e._yc >> e._a >> e._b >> e._rad >> e._score;

		if (!bIsAngleInRadians)
		{
			e._rad = e._rad * float(CV_PI / 180.0);
		}

		e._rad = fmod(float(e._rad + 2.0*CV_PI), float(CV_PI));

		if ((e._a > 0) && (e._b > 0) && (e._rad >= 0))
		{
			ellipses.push_back(e);
		}
	}
	in.close();

	// Sort ellipses by decreasing score
	sort(ellipses.begin(), ellipses.end());

	return true;
}
// Should be checked !!!!!
//  TestOverlap
float Evaluate(const vector<Ellipse>& ellGT, const vector<Ellipse>& ellTest, const float th_score, const Mat3b& img)
{
	float threshold_overlap = 0.8f;
	//float threshold = 0.95f;

	unsigned sz_gt = ellGT.size();
	unsigned size_test = ellTest.size();

	unsigned sz_test = unsigned(min(1000, int(size_test)));

	vector<Mat1b> gts(sz_gt);
	vector<Mat1b> tests(sz_test);
	//绘制每个目标椭圆
	for (unsigned i = 0; i<sz_gt; ++i)
	{
		const Ellipse& e = ellGT[i];

		Mat1b tmp(img.rows, img.cols, uchar(0));
		ellipse(tmp, Point((int)e._xc, (int)e._yc), Size((int)e._a, (int)e._b), e._rad * 180.0 / CV_PI, 0.0, 360.0, Scalar(255), -1);
		gts[i] = tmp;
	}
	//绘制检测的椭圆
	for (unsigned i = 0; i<sz_test; ++i)
	{
		const Ellipse& e = ellTest[i];

		Mat1b tmp(img.rows, img.cols, uchar(0));
		ellipse(tmp, Point((int)e._xc, (int)e._yc), Size((int)e._a, (int)e._b), e._rad * 180.0 / CV_PI, 0.0, 360.0, Scalar(255), -1);
		tests[i] = tmp;
	}

	Mat1b overlap(sz_gt, sz_test, uchar(0));
	for (int r = 0; r < overlap.rows; ++r)
	{
		for (int c = 0; c < overlap.cols; ++c)
		{
			//重叠区域占真个区域的比例 与比上并大于阈值 为255
			overlap(r, c) = TestOverlap(gts[r], tests[c], threshold_overlap) ? uchar(255) : uchar(0);
		}
	}

	int counter = 0;

	vector<bool> vec_gt(sz_gt, false);
	//矩阵每行有一个就代表找到
	for (unsigned int i = 0; i < sz_test; ++i)
	{
		//const Ellipse& e = ellTest[i];
		for (unsigned int j = 0; j < sz_gt; ++j)
		{
			if (vec_gt[j]) { continue; }

			bool bTest = overlap(j, i) != 0;

			if (bTest)
			{
				vec_gt[j] = true;
				break;
			}
		}
	}

	int tp = Count(vec_gt);
	int fn = int(sz_gt) - tp;
	int fp = size_test - tp; // !!!!

	float pr(0.f);
	float re(0.f);
	float fmeasure(0.f);

	if (tp == 0)
	{
		if (fp == 0)
		{
			pr = 1.f;
			re = 0.f;
			fmeasure = (2.f * pr * re) / (pr + re);
		}
		else
		{
			pr = 0.f;
			re = 0.f;
			fmeasure = 0.f;
		}
	}
	else
	{
		pr = float(tp) / float(tp + fp);
		re = float(tp) / float(tp + fn);
		fmeasure = (2.f * pr * re) / (pr + re);
	}

	return fmeasure;
}

bool TestOverlap(const Mat1b& gt, const Mat1b& test, float th)
{
	float fAND = float(countNonZero(gt & test));
	float fOR = float(countNonZero(gt | test));
	float fsim = fAND / fOR;

	return (fsim >= th);
}

int Count(const vector<bool> v)
{
	int counter = 0;
	for (unsigned i = 0; i<v.size(); ++i)
	{
		if (v[i]) { ++counter; }
	}
	return counter;
}

void salt(cv::Mat& image, int n){
	for(int k=0; k<n; k++){
		int i = rand()%image.cols;
		int j = rand()%image.rows;

		if(image.channels() == 1){
			image.at<uchar>(j,i) = 255;
		}else{
			image.at<cv::Vec3b>(j,i)[0] = 255;
			image.at<cv::Vec3b>(j,i)[1] = 255;
			image.at<cv::Vec3b>(j,i)[2] = 255;
		}
	}
}
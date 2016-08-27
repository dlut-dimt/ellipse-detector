/*
This code is intended for academic use only.
You are free to use and modify the code, at your own risk.

If you use this code, or find it useful, please refer to the paper:


The comments in the code refer to the abovementioned paper.
If you need further details about the code or the algorithm, please contact me at:

lianbosong@foxmail.com

last update: 
*/

#include "stdafx.h"
#include "tools.h"
#include "CNEllipseDetector.h"
using namespace std;
using namespace cv;

string SWORKINGDIR="EllipseDataset/";
string DBNAME="Dataset#2";//"good2";//"/PrasadImages-DatasetPrasad";//"/RandomImages-Dataset#1";//"/BRhoChange";//
string TESTIMGNAME="027_0003.jpg";//"163_0037.jpg";//"003_0144.jpg";//"193_0094.jpg";//"175_0026.jpg";//"im3949.jpg";//"177_0074.jpg";//"";//"140_0038.jpg";//"bike_0068.jpg";//
int MethodId=1;
//椭圆检测参数设置
float	fThScoreScore = 0.7f;	//0.8
float	fMinReliability	= 0.5f;	// Const parameters to discard bad ellipses 0.4
float	fTaoCenters = 0.05f;//0.05 	
int		ThLength=16;//16
float	MinOrientedRectSide=3.0f;

float	scale=1;

/*函数声明*/
float showT(string sWorkingDir,string imagename, CNEllipseDetector cned,vector<Ellipse> ellsCned,float fThScoreScore,bool showpic=false);
vector<double> OnImage(string sWorkingDir,string imagename,float fThScoreScore,float fMinReliability,float fTaoCenters,bool showpic=false);
vector<double> OnImage_salt(string sWorkingDir,string imagename,int saltrate,float fThScoreScore,float fMinReliability,float fTaoCenters,bool showpic=false);
vector<double> database(string dirName,float fThScoreScore,float fMinReliability,float fTaoCenters);
void OnVideo();
void SetParameter(map<char,string> Parameter);
void ReadParameter(int argc,char** argv);
/*函数实现*/
void showTime(vector<double> times){
	cout << "--------------------------------" << endl;
	cout << "Execution Time: " << endl;
	cout << "Edge Detection: \t" << times[0] << endl;
	cout << "Pre processing: \t" << times[1] << endl;
	cout << "Grouping:       \t" << times[2] << endl;
	cout << "Estimation:     \t" << times[3] << endl;
	cout << "Validation:     \t" << times[4] << endl;
	cout << "Clustering:     \t" << times[5] << endl;
	cout << "--------------------------------" << endl;
	cout << "Total:	         \t" << accumulate(times.begin(),times.begin()+6,0.0) << endl;
	cout << "F-Measure:      \t" << times[6] << endl;
	cout << "--------------------------------" << endl;
	if(times.size()==9){
		cout << "countsOfFindEllipse \t"<<times[7]<<endl;
		cout << "countsOfGetFastCenter \t"<<times[8]<<endl;
	}
}
float showT(string sWorkingDir,string imagename, CNEllipseDetector cned,vector<Ellipse> ellsCned,float fThScoreScore, bool showpic){
	vector<Ellipse> gt;
	LoadGT(gt, sWorkingDir + "/gt/" + "gt_" + imagename + ".txt", false);
	string filename = sWorkingDir + "/images/" + imagename;
	Mat3b image = imread(filename);
	Mat3b resultImage = image.clone();

	// Draw GT ellipses
	for (unsigned i = 0; i < gt.size(); ++i)
	{
		Ellipse& e = gt[i];
		Scalar color(0, 0, 255);
		ellipse(resultImage, Point(cvRound(e._xc), cvRound(e._yc)), Size(cvRound(e._a), cvRound(e._b)), e._rad*180.0 / CV_PI, 0.0, 360.0, color, 3);
	}

	cned.DrawDetectedEllipses(resultImage, ellsCned);

	Mat3b res = image.clone();
	float fmeasure = Evaluate(gt, ellsCned, fThScoreScore, res);

	if (showpic){
		cout << "F-Measure: " << fmeasure << endl;
		imshow("Cned", resultImage);
		cvSaveImage("result/resultImage.jpg",&IplImage(resultImage));
	}
	return fmeasure;
}
void OnVideo()
{
	VideoCapture cap(0);
	if(!cap.isOpened()) return;

	CNEllipseDetector cned;

	Mat1b gray;
	while(true)
	{	
		Mat3b image;
		cap >> image;
		cvtColor(image, gray, CV_BGR2GRAY);	

		vector<Ellipse> ellipses;

		//Find Ellipses		
		cned.Detect(gray, ellipses);
		cned.DrawDetectedEllipses(image,ellipses);
		imshow("Output", image);

			
		if(waitKey(10) >= 0) break;
	}
}

vector<double> OnImage(string filename,float fThScoreScore,float fMinReliability,float fTaoCenters,bool showpic)
{
	vector<double> times;
	CNEllipseDetector cned;
	// Read image
	Mat3b image = imread(filename);
	
	if(!image.data){
		cout<<filename<<" not exist"<<endl;
		return times;
	}
	Size sz = image.size();

	// Convert to grayscale
	Mat1b gray;
	cvtColor(image, gray, CV_BGR2GRAY);

	// Parameters Settings (Sect. 4.2)
	int		iThLength = ThLength;//过滤太短的弧
	float	fThObb = MinOrientedRectSide;
	float	fThPos = 1.0f;
	//float	fTaoCenters = 0.018f;//0.05
	int 	iNs = 16;//弦数
	float	fMaxCenterDistance = sqrt(float(sz.width*sz.width + sz.height*sz.height)) * fTaoCenters;

	//float	fThScoreScore = 0.5f;//0.8	
	//fTaoCenters = 0.05f;
	// Other constant parameters settings.
	// Gaussian filter parameters, in pre-processing
	Size	szPreProcessingGaussKernelSize	= Size(5,5);
	double	dPreProcessingGaussSigma		= 1.0;

	float	fDistanceToEllipseContour		= 0.1f;	// (Sect. 3.3.1 - Validation)
	//float	fMinReliability					= 0.4f;	// Const parameters to discard bad ellipses 0.5

	// Initialize Detector with selected parameters
	cned.SetParameters	(	szPreProcessingGaussKernelSize,	
		dPreProcessingGaussSigma,		
		fThPos,
		fMaxCenterDistance,
		iThLength,
		fThObb,
		fDistanceToEllipseContour,		
		fThScoreScore,
		fMinReliability,		
		iNs
		);
	// Detect 
	vector<Ellipse> ellsCned;
	Mat1b gray_clone=gray.clone();
	cned.Detect(gray_clone, ellsCned);
	times = cned.GetTimes();

	Mat3b resultImage = image.clone();
	cned.DrawDetectedEllipses(resultImage, ellsCned);
	imshow("Cned", resultImage);
	if(showpic){
		_mkdir("result");
		cvSaveImage("result/resultImage.jpg",&IplImage(resultImage));
		SaveEllipses("result/result.txt", ellsCned);
	}
	double fmeasure = 0;//showT(sWorkingDir,imagename, cned,ellsCned,0.8f,showpic);
	times.push_back(fmeasure);
	times.push_back(cned.countsOfFindEllipse);
	times.push_back(cned.countsOfGetFastCenter);
	return times;
}

void SetParameter(map<char,string> Parameter){
	/**
	-N image Name
	-D DataSet Name
	-S The threshold of ellipse score
	-R The threshold of Reliability
	-C The threshold of CenterDis
	-M The method id
	-P
	*/
	map<char,string>::iterator it;
	for(it=Parameter.begin();it!=Parameter.end();++it){
		switch(it->first){
			case 'N': TESTIMGNAME=it->second;break;
			case 'D': DBNAME=it->second;break;
			case 'S': fThScoreScore=atof(it->second.c_str());break;
			case 'C': fTaoCenters=atof(it->second.c_str());break;
			case 'R': fMinReliability=atof(it->second.c_str());break;
			case 'M': MethodId=atof(it->second.c_str());break;
			case 'P': SWORKINGDIR=it->second.c_str();break;
		}
	}
}
void ReadParameter(int argc,char** argv){
    if(argc<=1||argv[1][0]!='-'){
		if(argc>1)cout<<argv[1]<<endl;
		return;
    }
    map<char,string> Parameter;
    int targc=argc;
    int cur=1,pre=1;
	string sPa;
	char cP=argv[pre][1];
    for(int i=2;i<argc;i++){
        string stmp=argv[i];
        if(0==stmp.find("-")){
            cur=i;
			cP=argv[pre][1];
			Trim(sPa);
            Parameter.insert(pair<char,string>(cP,sPa));
			pre=cur;sPa="";
        }else{
			sPa=sPa+" "+argv[i];
		}
    }
	cP=argv[pre][1];
	Trim(sPa);
    Parameter.insert(pair<char,string>(cP,sPa));
	map<char,string>::iterator it;
    for(it=Parameter.begin();it!=Parameter.end();++it)
        cout<<"key: "<<it->first <<" value: "<<it->second<<endl;
	SetParameter(Parameter);
	
}


vector<double> OnImage(string sWorkingDir,string imagename,float fThScoreScore,float fMinReliability,float fTaoCenters,bool showpic)
{
	CNEllipseDetector cned;

	string filename = sWorkingDir + "/images/" + imagename;
	// Read image
	Mat3b image = imread(filename);
	if(!image.data){
		cout<<filename<<" not exist"<<endl;
		return vector<double>(0);
	}
	Size sz = image.size();

	// Convert to grayscale
	Mat1b gray;
	cvtColor(image, gray, CV_BGR2GRAY);

	// Parameters Settings (Sect. 4.2)
	int		iThLength = ThLength;//过滤太短的弧
	float	fThObb = MinOrientedRectSide;
	float	fThPos = 1.0f;
	//float	fTaoCenters = 0.018f;//0.05
	int 	iNs = 16;//弦数
	float	fMaxCenterDistance = sqrt(float(sz.width*sz.width + sz.height*sz.height)) * fTaoCenters;

	//float	fThScoreScore = 0.5f;//0.8	
	//fTaoCenters = 0.05f;
	// Other constant parameters settings.
	// Gaussian filter parameters, in pre-processing
	Size	szPreProcessingGaussKernelSize	= Size(5,5);
	double	dPreProcessingGaussSigma		= 1.0;

	float	fDistanceToEllipseContour		= 0.1f;	// (Sect. 3.3.1 - Validation)
	//float	fMinReliability					= 0.4f;	// Const parameters to discard bad ellipses 0.5

	// Initialize Detector with selected parameters
	cned.SetParameters	(	szPreProcessingGaussKernelSize,	
		dPreProcessingGaussSigma,		
		fThPos,
		fMaxCenterDistance,
		iThLength,
		fThObb,
		fDistanceToEllipseContour,		
		fThScoreScore,
		fMinReliability,		
		iNs
		);
	// Detect 
	vector<Ellipse> ellsCned;
	Mat1b gray_clone=gray.clone();
	cned.Detect(gray_clone, ellsCned);
	vector<double> times = cned.GetTimes();
	if(showpic){
		_mkdir("result");
		SaveEllipses("result/result.txt", ellsCned);
	}
	double fmeasure = showT(sWorkingDir,imagename, cned,ellsCned,0.8f,showpic);
	times.push_back(fmeasure);
	times.push_back(cned.countsOfFindEllipse);
	times.push_back(cned.countsOfGetFastCenter);
	return times;
}
vector<double> database(string sWorkingDir,float fThScoreScore,float fMinReliability,float fTaoCenters){
	vector<string> resultString;
	resultString.push_back("picName,Edge Detection,Pre processing,Grouping,Estimation,Validation,Clustering,WholeTime,F-measure,countsOfFindEllipse,countsOfGetFastCenter");
	vector<double> allTimes(9,0.0);
	string dirName=sWorkingDir+"/images/";
	string imageName;
	double countsOfFindEllipse=0;
	double countsOfGetFastCenter=0;
	vector<string> files;
	listDir(dirName,files);

	imageName=files.at(1);
	vector<double> results=OnImage(sWorkingDir,imageName,fThScoreScore,fMinReliability,fTaoCenters);
	double fmeasure = results[6];
	vector<double> times = results;
	double wholetime=accumulate(times.begin(),times.begin()+6,0.0);
	fmeasure=0.0;
	for(unsigned int i=0;i<files.size();i++)
	{
		stringstream resultsitem;
		imageName=files.at(i);
		results=OnImage(sWorkingDir,imageName,fThScoreScore,fMinReliability,fTaoCenters);
		fmeasure += results[6];
		countsOfFindEllipse+=results[7];
		countsOfGetFastCenter+=results[8];
		times = results;
		wholetime=accumulate(times.begin(),times.begin()+6,0.0);
		for (unsigned int j=0;j<allTimes.size();j++){
			allTimes[j]+=times[j];
		}
		resultsitem<<imageName<<","<<times[0]<<","<<times[1]<<","<<times[2]<<","<<times[3]<<","<<times[4]<<","<<times[5]<<","<<wholetime<<","<<times[6]<<","<<times[7]<<","<<times[8];
		resultString.push_back(resultsitem.str());
	}
	for (unsigned int i=0;i<allTimes.size();i++){
		times[i]=allTimes[i]/(1.0*files.size());
	}
	wholetime=accumulate(times.begin(),times.begin()+6,0.0);
	fmeasure=fmeasure/(1.0*files.size());
	cout << "--------------------------------" << endl;
	cout << "DataSet:        \t"<<sWorkingDir<<endl;
	cout << "Execution Time" << endl;
	cout << "Edge Detection: \t" << times[0] << endl;
	cout << "Pre processing: \t" << times[1] << endl;
	cout << "Grouping:       \t" << times[2] << endl;
	cout << "Estimation:     \t" << times[3] << endl;
	cout << "Validation:     \t" << times[4] << endl;
	cout << "Clustering:     \t" << times[5] << endl;
	cout << "--------------------------------" << endl;
	cout << "Total:	         \t" << accumulate(times.begin(),times.begin()+6,0.0) << endl;
	cout << "F-Measure: " << fmeasure << endl;
	cout << "--------------------------------" << endl;
	cout << "countsOfFindEllipse \t"<<times[7]<<endl;
	cout << "countsOfGetFastCenter \t"<<times[8]<<endl;
	stringstream resultsitem;
	resultsitem<<"average"<<","<<times[0]<<","<<times[1]<<","<<times[2]<<","<<times[3]<<","<<times[4]<<","<<times[5]<<","<<wholetime<<","<<times[6]<<","<<times[7]<<","<<times[8];
	resultString.push_back(resultsitem.str());
	writeFile("dataset.csv",resultString);
	return times;
}

int main_OnePic()
{
	string sWorkingDir = SWORKINGDIR+DBNAME;
	string imagename = TESTIMGNAME;

	vector<double> results=OnImage(sWorkingDir,imagename,fThScoreScore,fMinReliability,fTaoCenters,true);
	if(!results.empty()){
		showTime(results);//显示运行信息
		waitKey(0);
		cvDestroyWindow("Cned");
	}
	
	return 0;
}
int main_OneDB()
{
	string sWorkingDir = SWORKINGDIR+DBNAME;
	string imagename = TESTIMGNAME;
	cout<<DBNAME<<"\t"<<fThScoreScore<<" "<<fMinReliability<<" "<<fTaoCenters<<endl;
	database(sWorkingDir,fThScoreScore,fMinReliability,fTaoCenters);
	return 0;
}
int main_allDB(int argc, char** argv)
{
	string sWorkingDirPath=SWORKINGDIR;
	string sWorkingDirName[3] = {"PrasadImages-DatasetPrasad","RandomImages-Dataset#1","Dataset#2"};
	float afThScoreScore[6]={0.3f,0.4f,0.5f,0.6f,0.7f,0.8f};//0-5
	float afMinReliability[6]={0.3f,0.4f,0.5f,0.6f,0.7f,0.8f};//0-5
	float afTaoCenters[6]={0.01f,0.02f,0.03f,0.04f,0.05f,0.6f};//0-5
	int iDir=0;
	int iScore=0;
	int iReliability=0;
	int iTaoCenter=0;
	string sWorkingDir;

	float	fThScoreScore = 0.5f;	//0.8
	float	fMinReliability	= 0.50f;	// Const parameters to discard bad ellipses 0.4
	float	fTaoCenters = 0.05f;//0.05 
	string oldname="dataset.csv";
	string newname="";

	vector<double> times;
	char oneTimeName[100];
	
	char cnewDir[100];
	sprintf(cnewDir,"normal/");
	string newDir=cnewDir;
	_mkdir(cnewDir);
	vector<string> resultString;
	resultString.push_back("iDir,fThScoreScore,fMinReliability,fTaoCenters,Edge Detection,Pre processing,Grouping,Estimation,Validation,Clustering,WholeTime,F-measure,countsOfFindEllipse,countsOfGetFastCenter");

	for(iDir=0;iDir<3;iDir++){
		sWorkingDir=sWorkingDirPath+sWorkingDirName[iDir];
		for(iScore=0;iScore<6;iScore++){
			fThScoreScore = afThScoreScore[iScore];
			for(iReliability=0;iReliability<6;iReliability++){
				system("cls");
				fMinReliability = afMinReliability[iReliability];
				for(iTaoCenter=0;iTaoCenter<6;iTaoCenter++){
					fTaoCenters = afTaoCenters[iTaoCenter];
					times=database(sWorkingDir,fThScoreScore,fMinReliability,fTaoCenters);
					double wholetime=accumulate(times.begin(),times.begin()+6,0.0);
					sprintf(oneTimeName,"our_%s_%4.03f_%4.03f_%4.03f.csv",sWorkingDirName[iDir].c_str(),fThScoreScore,fMinReliability,fTaoCenters);
					newname=newDir+oneTimeName;
					rename( oldname.c_str(),newname.c_str());
					cout<<newname<<endl;
					stringstream resultsitem;

					resultsitem<<iDir<<","<<fThScoreScore<<","<<fMinReliability<<","<<fTaoCenters<<","
						<<times[0]<<","<<times[1]<<","<<times[2]<<","<<times[3]<<","<<times[4]<<","<<times[5]<<","
						<<wholetime<<","<<times[6]<<","<<times[7]<<","<<times[8];
					resultString.push_back(resultsitem.str());
				}
			}
		}
	}
	
	writeFile(newDir+"our_allaverage.csv",resultString);
	return 0;
}
//salt 
vector<double> OnImage_salt(string sWorkingDir,string imagename,int saltrate,float fThScoreScore,float fMinReliability,float fTaoCenters,bool showpic)
{
	CNEllipseDetector cned;

	string filename = sWorkingDir + "/images/" + imagename;
	// Read image
	Mat3b image = imread(filename);
	Size sz = image.size();
	int n=sz.width*sz.height;
	
	// Convert to grayscale
	//medianBlur(image.clone(),image,3);

	// gray
	Mat1b gray;
	cvtColor(image, gray, CV_BGR2GRAY);
	//salt 添加
	salt(gray, n*saltrate/100);
	medianBlur(gray,gray,3);
	// Parameters Settings (Sect. 4.2)
	int		iThLength = 16;//过滤太短的弧
	float	fThObb = 3.0f;
	float	fThPos = 1.0f;
	//float	fTaoCenters = 0.018f;//0.05
	int 	iNs = 16;//弦数
	float	fMaxCenterDistance = sqrt(float(sz.width*sz.width + sz.height*sz.height)) * fTaoCenters;

	//float	fThScoreScore = 0.5f;//0.8	
	//fTaoCenters = 0.05f;
	// Other constant parameters settings.
	// Gaussian filter parameters, in pre-processing
	Size	szPreProcessingGaussKernelSize	= Size(5,5);
	double	dPreProcessingGaussSigma		= 1.0;

	float	fDistanceToEllipseContour		= 0.1f;	// (Sect. 3.3.1 - Validation)
	//float	fMinReliability					= 0.4f;	// Const parameters to discard bad ellipses 0.5

	// Initialize Detector with selected parameters
	cned.SetParameters	(	szPreProcessingGaussKernelSize,	
		dPreProcessingGaussSigma,		
		fThPos,
		fMaxCenterDistance,
		iThLength,
		fThObb,
		fDistanceToEllipseContour,		
		fThScoreScore,
		fMinReliability,		
		iNs
		);
	// Detect
	vector<Ellipse> ellsCned;
	Mat1b gray_clone=gray.clone();
	cned.Detect(gray_clone, ellsCned);
	vector<double> times = cned.GetTimes();
	double fmeasure = showT(sWorkingDir,imagename, cned,ellsCned,0.8f,showpic);
	times.push_back(fmeasure);
	times.push_back(cned.countsOfFindEllipse);
	times.push_back(cned.countsOfGetFastCenter);
	return times;
}
vector<double> database_salt(string sWorkingDir,int saltrate,float fThScoreScore,float fMinReliability,float fTaoCenters){
	vector<string> resultString;
	resultString.push_back("picName,Edge Detection,Pre processing,Grouping,Estimation,Validation,Clustering,WholeTime,F-measure,countsOfFindEllipse,countsOfGetFastCenter");
	vector<double> allTimes(9,0.0);
	string dirName=sWorkingDir+"/images/";
	string imageName;
	double countsOfFindEllipse=0;
	double countsOfGetFastCenter=0;
	vector<string> files;
	listDir(dirName,files);

	imageName=files.at(1);
	vector<double> results=OnImage(sWorkingDir,imageName,fThScoreScore,fMinReliability,fTaoCenters);
	double fmeasure = results[6];
	vector<double> times = results;
	double wholetime=accumulate(times.begin(),times.begin()+6,0.0);
	fmeasure=0.0;
	for(unsigned int i=0;i<files.size();i++)
	{
		stringstream resultsitem;
		imageName=files.at(i);
		results=OnImage_salt(sWorkingDir,imageName,saltrate,fThScoreScore,fMinReliability,fTaoCenters);
		fmeasure += results[6];
		countsOfFindEllipse+=results[7];
		countsOfGetFastCenter+=results[8];
		times = results;
		wholetime=accumulate(times.begin(),times.begin()+6,0.0);

		for (unsigned int j=0;j<allTimes.size();j++){
			allTimes[j]+=times[j];
		}
		resultsitem<<imageName<<","<<times[0]<<","<<times[1]<<","<<times[2]<<","<<times[3]<<","<<times[4]<<","<<times[5]<<","<<wholetime<<","<<times[6]<<","<<times[7]<<","<<times[8];
		resultString.push_back(resultsitem.str());
		//if((i+1)%10==0){cout<<endl;}//10个数据一个回车
	}
	for (unsigned int i=0;i<allTimes.size();i++){
		times[i]=allTimes[i]/(1.0*files.size());
	}
	wholetime=accumulate(times.begin(),times.begin()+6,0.0);
	fmeasure=fmeasure/(1.0*files.size());
	cout << "--------------------------------" << endl;
	cout << "DataSet:        \t"<<sWorkingDir<<endl;
	cout << "Execution Time" << endl;
	cout << "Edge Detection: \t" << times[0] << endl;
	cout << "Pre processing: \t" << times[1] << endl;
	cout << "Grouping:       \t" << times[2] << endl;
	cout << "Estimation:     \t" << times[3] << endl;
	cout << "Validation:     \t" << times[4] << endl;
	cout << "Clustering:     \t" << times[5] << endl;
	cout << "--------------------------------" << endl;
	cout << "Total:	         \t" << accumulate(times.begin(),times.begin()+6,0.0) << endl;
	cout << "F-Measure: " << fmeasure << endl;
	cout << "--------------------------------" << endl;
	cout << "countsOfFindEllipse \t"<<times[7]<<endl;
	cout << "countsOfGetFastCenter \t"<<times[8]<<endl;
	stringstream resultsitem;
	resultsitem<<"average"<<","<<times[0]<<","<<times[1]<<","<<times[2]<<","<<times[3]<<","<<times[4]<<","<<times[5]<<","<<wholetime<<","<<times[6]<<","<<times[7]<<","<<times[8];
	resultString.push_back(resultsitem.str());
	writeFile("dataset.csv",resultString);
	return times;
}
int main_salt(int argc, char** argv)
{//salt
	string sWorkingDirPath=SWORKINGDIR;
	string sWorkingDirName[3] = {"PrasadImages-DatasetPrasad","RandomImages-Dataset#1","Dataset#2"};
	//string sWorkingDirName[3] = {"good2"};
	float afThScoreScore[6]={0.3f,0.4f,0.5f,0.6f,0.7f,0.8f};//0-5
	float afMinReliability[6]={0.3f,0.4f,0.5f,0.6f,0.7f,0.8f};//0-5
	float afTaoCenters[6]={0.01f,0.02f,0.03f,0.04f,0.05f,0.6f};//0-5
	int saltrates[6]={3,6,9,12,15,18};//0-5
	int iDir=0;
	int iScore=0;
	int iReliability=0;
	int iTaoCenter=0;
	int isaltrate=0;
	string sWorkingDir;
	sWorkingDir=sWorkingDirPath+sWorkingDirName[iDir];


	float	fThScoreScore = 0.5f;	//0.8
	float	fMinReliability	= 0.50f;	// Const parameters to discard bad ellipses 0.4
	float	fTaoCenters = 0.05f;//0.05 
	string oldname="dataset.csv";
	string newname="";

	vector<double> times;
	char oneTimeName[100];
	char partDataName[50];
	for(isaltrate=0;isaltrate<6;isaltrate++){
		int saltrate=saltrates[isaltrate];
		char cnewDir[20];
		sprintf(cnewDir,"%02d/",saltrate);
		string newDir=cnewDir;
		_mkdir(cnewDir);
		vector<string> resultString;
		resultString.push_back("iDir,fThScoreScore,fMinReliability,fTaoCenters,Edge Detection,Pre processing,Grouping,Estimation,Validation,Clustering,Totaltimes,F-measure,countsOfFindEllipse,countsOfGetFastCenter");
		for(iDir=0;iDir<3;iDir++){
			sWorkingDir=sWorkingDirPath+sWorkingDirName[iDir];
			for(iScore=0;iScore<6;iScore++){
				fThScoreScore = afThScoreScore[iScore];
				for(iReliability=0;iReliability<6;iReliability++){
					fMinReliability = afMinReliability[iReliability]; 
					for(iTaoCenter=0;iTaoCenter<6;iTaoCenter++){
						fTaoCenters = afTaoCenters[iTaoCenter];
						times=database_salt(sWorkingDir,saltrate,fThScoreScore,fMinReliability,fTaoCenters);
						double wholetime=accumulate(times.begin(),times.begin()+6,0.0);
						sprintf(oneTimeName,"our_%s_%02d_%4.03f_%4.03f_%4.03f.csv",sWorkingDirName[iDir].c_str(),saltrate,fThScoreScore,fMinReliability,fTaoCenters);
						newname=newDir+oneTimeName;
						rename( oldname.c_str(),newname.c_str());
						cout<<newname<<endl;
						stringstream resultsitem;
						resultsitem<<iDir<<","<<fThScoreScore<<","<<fMinReliability<<","<<fTaoCenters<<","
							<<times[0]<<","<<times[1]<<","<<times[2]<<","<<times[3]<<","<<times[4]<<","<<times[5]<<","
							<<wholetime<<","<<times[6]<<","<<times[7]<<","<<times[8];
						resultString.push_back(resultsitem.str());
					}
				}
			}
		}
		sprintf(partDataName,"%02d/our_%02d_allaverage.csv",saltrate,saltrate);
		writeFile(partDataName,resultString);
	}
	return 0;
}

vector<double> detect4s(string sWorkingDir, string imageName){
	vector<double>  Counts;
	string wholeImageName=sWorkingDir+"/images/"+imageName;

	bool showpic=false;

	//show edges
	Mat1b gray;
	cvtColor(imread(wholeImageName), gray, CV_BGR2GRAY);
	CNEllipseDetector cned;
	int		iThLength = ThLength;//过滤太短的弧
	float	fThPos = 1.0f;
	float	fThObb = MinOrientedRectSide;

	float	fMaxCenterDistance = 10;
	int 	iNs = 16;//弦数

	Size	szPreProcessingGaussKernelSize	= Size(5,5);
	double	dPreProcessingGaussSigma		= 1.0;

	float	fDistanceToEllipseContour		= 0.1f;	// (Sect. 3.3.1 - Validation)

	// Initialize Detector with selected parameters
	cned.SetParameters	(	szPreProcessingGaussKernelSize,	
		dPreProcessingGaussSigma,		
		fThPos,
		fMaxCenterDistance,
		iThLength,
		fThObb,
		fDistanceToEllipseContour,		
		fThScoreScore,
		fMinReliability,		
		iNs
		);

	Mat1b gray_clone=gray.clone();
	int EdgeNumber=cned.showEdgeInPic(gray_clone,showpic);
	Counts.push_back(EdgeNumber);

	vector<double> results;
	// 用于防止缓存影响速度
	myselect1=false;myselect2=false;myselect3=false;
	results=OnImage(sWorkingDir,imageName,fThScoreScore,fMinReliability,fTaoCenters);
	double fmeasure = results[6];
	double wholetime=accumulate(results.begin(),results.begin()+6,0.0);
	// 4个
	myselect1=true;myselect2=false;myselect3=false;
	results=OnImage(sWorkingDir,imageName,fThScoreScore,fMinReliability,fTaoCenters);
	fmeasure = results[6];
	wholetime=accumulate(results.begin(),results.begin()+6,0.0);
	Counts.push_back(results[7]);
	Counts.push_back(results[8]);
	Counts.push_back(wholetime);
	Counts.push_back(fmeasure);

	myselect1=true;myselect2=true;myselect3=false;
	results=OnImage(sWorkingDir,imageName,fThScoreScore,fMinReliability,fTaoCenters);
	fmeasure = results[6];
	wholetime=accumulate(results.begin(),results.begin()+6,0.0);
	Counts.push_back(results[7]);
	Counts.push_back(results[8]);
	Counts.push_back(wholetime);
	Counts.push_back(fmeasure);

	myselect1=true;myselect2=true;myselect3=true;
	results=OnImage(sWorkingDir,imageName,fThScoreScore,fMinReliability,fTaoCenters);
	fmeasure = results[6];
	wholetime=accumulate(results.begin(),results.begin()+6,0.0);
	Counts.push_back(results[7]);
	Counts.push_back(results[8]);
	Counts.push_back(wholetime);
	Counts.push_back(fmeasure);

	myselect1=false;myselect2=false;myselect3=false;
	results=OnImage(sWorkingDir,imageName,fThScoreScore,fMinReliability,fTaoCenters);
	fmeasure = results[6];
	wholetime=accumulate(results.begin(),results.begin()+6,0.0);
	Counts.push_back(results[7]);
	Counts.push_back(results[8]);
	Counts.push_back(wholetime);
	Counts.push_back(fmeasure);

	return Counts;
}
vector<double> database_4s(string sWorkingDir){
	vector<string> resultString;
	resultString.push_back("ImageName,EdgeNumber,cOFE1,cOGFC1,select100Time,F-m1,cOFE2,cOGFC2,select110Time,F-m2,cOFE3,cOGFC3,select111Time,F-m3,cOFE,cOGFC,select000Time,F-m");

	string dirName=sWorkingDir+"/images/";
	vector<string> files;
	listDir(dirName,files);
	string imageName;

	imageName=files.at(1);
	vector<double> Counts=detect4s(sWorkingDir,imageName);
	vector<double> AllCounts(17,0.0);
	for(unsigned int i=0;i<files.size();i++)
	{
		imageName=files.at(i);
		Counts=detect4s(sWorkingDir,imageName);
		stringstream resultsitem;
		resultsitem<<imageName<<","<<Counts[0]
			<<","<<Counts[1]<<","<<Counts[2]<<","<<Counts[3]<<","<<Counts[4]
			<<","<<Counts[5]<<","<<Counts[6]<<","<<Counts[7]<<","<<Counts[8]
			<<","<<Counts[9]<<","<<Counts[10]<<","<<Counts[11]<<","<<Counts[12]
			<<","<<Counts[13]<<","<<Counts[14]<<","<<Counts[15]<<","<<Counts[16];
		resultString.push_back(resultsitem.str());
		
		for (unsigned int j=0;j<AllCounts.size();j++){
			AllCounts[j]+=Counts[j];
		}
	}
	for (unsigned int j=0;j<AllCounts.size();j++){
		Counts[j]=AllCounts[j]/(files.size());
	}
	cout << "--------------------------------" << endl;
	cout << "--------------------------------" << endl;
	cout << "DataSet:        \t"<<sWorkingDir<<endl;
	cout<<"Edge Number:\t"<<Counts[0]<<endl;
	cout<<"tTCNl:\t"<<tTCNl<<endl;
	cout<<"tCNC:\t"<<tCNC<<endl;
	cout << "--------------------------------" << endl;
	cout<<"countsOfFindEllipse:\t"<<Counts[1]<<endl;
	cout<<"countsOfGetFastCenter\t"<<Counts[2]<<endl<<endl;
	cout<<"exec_time1:\t\t"<<Counts[3]<<endl;
	cout<<"F-M 1:\t"<<Counts[4]<<endl;
	cout << "--------------------------------" << endl;
	cout<<"countsOfFindEllipse:\t"<<Counts[5]<<endl;
	cout<<"countsOfGetFastCenter\t"<<Counts[6]<<endl<<endl;
	cout<<"exec_time2:\t\t"<<Counts[7]<<endl;
	cout<<"F-M 2:\t"<<Counts[8]<<endl;
	cout << "--------------------------------" << endl;
	cout<<"countsOfFindEllipse:\t"<<Counts[9]<<endl;
	cout<<"countsOfGetFastCenter\t"<<Counts[10]<<endl<<endl;
	cout<<"exec_time3:\t\t"<<Counts[11]<<endl;
	cout<<"F-M 3:\t"<<Counts[12]<<endl;
	cout << "--------------------------------" << endl;
	cout<<"countsOfFindEllipse:\t"<<Counts[13]<<endl;
	cout<<"countsOfGetFastCenter\t"<<Counts[14]<<endl<<endl;
	cout<<"exec_time:\t\t"<<Counts[15]<<endl;
	cout<<"F-M :\t"<<Counts[16]<<endl;
	cout << "--------------------------------" << endl;
	cout<<"rate:\t"<<Counts[3]/Counts[15]<<"\t"<<Counts[7]/Counts[15]
	<<"\t"<<Counts[11]/Counts[15]<<endl;
	stringstream resultsitem;
	resultsitem<<"average"<<","<<Counts[0]
			<<","<<Counts[1]<<","<<Counts[2]<<","<<Counts[3]<<","<<Counts[4]
			<<","<<Counts[5]<<","<<Counts[6]<<","<<Counts[7]<<","<<Counts[8]
			<<","<<Counts[9]<<","<<Counts[10]<<","<<Counts[11]<<","<<Counts[12]
			<<","<<Counts[13]<<","<<Counts[14]<<","<<Counts[15]<<","<<Counts[16];
	resultString.push_back(resultsitem.str());
	writeFile("dataset.csv",resultString);
	return Counts;
}
int main_TCN()
{// 遍历tTCNl
	string sWorkingDirPath=SWORKINGDIR;
	string sWorkingDirName[3] = {"PrasadImages-DatasetPrasad","RandomImages-Dataset#1","good2"};
	
	//float tTCN[6]={0,1.0f/64,2.0f/64,4.0f/64,6.0f/64,8.0f/64};
	float tTCN[6]={0,1,2,3,4,5};
	vector<double>  Counts;

	int iDir=0;
	string sWorkingDir;
	sWorkingDir=sWorkingDirPath+sWorkingDirName[iDir];

	float	fThScoreScore = 0.5f;	//0.8
	float	fMinReliability	= 0.50f;	// Const parameters to discard bad ellipses 0.4
	float	fTaoCenters = 0.05f;//0.05 
	string oldname="dataset.csv";
	string newname="";

	vector<double> times;
	
	char oneTimeName[100];
	vector<string> resultString;
	resultString.push_back("iDir,tTCNl,EdgeNumber,select100Time,cOFE1,cOGFC1,F-m1,select110Time,cOFE2,cOGFC2,F-m2,select111Time,cOFE3,cOGFC3,F-m3,select000Time,cOFE,cOGFC,F-m");
	MinOrientedRectSide=0;
	for(iDir=0;iDir<3;iDir++){
		for(int itTCN=0;itTCN<6;itTCN++){
			tTCNl=tTCN[itTCN];
			MinOrientedRectSide=tTCNl;
			cout<<tTCNl<<endl;
			sWorkingDir=sWorkingDirPath+sWorkingDirName[iDir];
			Counts=database_4s(sWorkingDir);
			sprintf(oneTimeName,"EdgeCompare_%s_%4.03f.csv",sWorkingDirName[iDir].c_str(),tTCNl);
			newname=oneTimeName;
			rename( oldname.c_str(),newname.c_str());
			stringstream resultsitem;
			resultsitem<<iDir<<","<<tTCNl<<","<<Counts[0]
				<<","<<Counts[1]<<","<<Counts[2]<<","<<Counts[3]<<","<<Counts[4]
				<<","<<Counts[5]<<","<<Counts[6]<<","<<Counts[7]<<","<<Counts[8]
				<<","<<Counts[9]<<","<<Counts[10]<<","<<Counts[11]<<","<<Counts[12]
				<<","<<Counts[13]<<","<<Counts[14]<<","<<Counts[15]<<","<<Counts[16];
			resultString.push_back(resultsitem.str());
		}
		
	}
	
	writeFile("our_allEdgeCompareaverage.csv",resultString);
	return 0;
}
int main_CNC(){
	string sWorkingDirPath=SWORKINGDIR;
	string sWorkingDirName[3] = {"PrasadImages-DatasetPrasad","RandomImages-Dataset#1","Dataset#2"};
	
	float tCNCs[21]={0,0.1f,0.2f,0.3f,0.4f,0.5f,0.6f,1,2,3,4,5,10,20,30,40,50,100,200,500,1000};
	vector<double>  Counts;

	int iDir=0;
	string sWorkingDir;
	sWorkingDir=sWorkingDirPath+sWorkingDirName[iDir];
	string oldname="dataset.csv";
	string newname="";

	vector<double> times;
	
	char oneTimeName[100];
	vector<string> resultString;
	resultString.push_back("iDir,tTCNl,EdgeNumber,cOFE1,cOGFC1,select100Time,F-m1,cOFE2,cOGFC2,select110Time,F-m2,cOFE3,cOGFC3,select111Time,F-m3,cOFE,cOGFC,select000Time,F-m");
	MinOrientedRectSide=0;
	for(iDir=0;iDir<3;iDir++){
		for(int itCNCs=0;itCNCs<21;itCNCs++){
			tCNC=tCNCs[itCNCs];
			cout<<tCNC<<endl;
			sWorkingDir=sWorkingDirPath+sWorkingDirName[iDir];
			Counts=database_4s(sWorkingDir);
			sprintf(oneTimeName,"EdgeCompare_CNC_%s_%4.03f.csv",sWorkingDirName[iDir].c_str(),tCNC);
			newname=oneTimeName;
			rename( oldname.c_str(),newname.c_str());
			stringstream resultsitem;
			resultsitem<<iDir<<","<<tCNC<<","<<Counts[0]
				<<","<<Counts[1]<<","<<Counts[2]<<","<<Counts[3]<<","<<Counts[4]
				<<","<<Counts[5]<<","<<Counts[6]<<","<<Counts[7]<<","<<Counts[8]
				<<","<<Counts[9]<<","<<Counts[10]<<","<<Counts[11]<<","<<Counts[12]
				<<","<<Counts[13]<<","<<Counts[14]<<","<<Counts[15]<<","<<Counts[16];
			resultString.push_back(resultsitem.str());
		}
		
	}
	
	writeFile("our_allEdgeCompareCNCaverage.csv",resultString);
	return 0;
}

int main(int argc, char** argv)
{
	tCNC=0.2f;
	fThScoreScore = 0.6f;	//0.8
	fMinReliability	= 0.4f;	// Const parameters to discard bad ellipses 0.4
	fTaoCenters = 0.04f;//0.05 	
	ThLength=16;//16
	MinOrientedRectSide=3.0f;
	if(argc==2){
		string filename= argv[1];
		if(NULL!=strstr(filename.c_str(),".jpg")||NULL!=strstr(filename.c_str(),".bmp")
			||NULL!=strstr(filename.c_str(),".JPG")||NULL!=strstr(filename.c_str(),".png")){	
			cout<<filename<<endl;
			vector<double> results=OnImage(filename,fThScoreScore,fMinReliability, fTaoCenters, true);
			if(!results.empty()){
				showTime(results);//显示运行信息
				waitKey(0);	cvDestroyWindow("Cned");
				MethodId=0;
				system("pause");
			}
		}
	}
	ReadParameter(argc,argv);
	switch(MethodId){
	case 1:main_OnePic();//单张图片
		break;
	case 2:main_OneDB();//单个数据库
		break;
	case 3:main_allDB(argc, argv);//论文中三个数据库，不同参数的实验
		break;
	case 4:main_salt(argc, argv);
		break;
	case 5:main_TCN();
		break;
	case 6:main_CNC();
		break;
	case 7:OnVideo();
		break;
	case 8:{
		main_allDB(argc, argv);

		tCNC=0.2f;ThLength=16;tTCNl=3;
		fThScoreScore = 0.7f;fMinReliability = 0.5f;fTaoCenters = 0.05f;//0.05 	
		MinOrientedRectSide=3.0f;
		main_salt(argc, argv);

		tCNC=0.2f;ThLength=16;tTCNl=3;
		fThScoreScore = 0.7f;fMinReliability = 0.5f;fTaoCenters = 0.05f;//0.05 	
		MinOrientedRectSide=3.0f;
		main_TCN();

		tCNC=0.2f;ThLength=16;tTCNl=3;
		fThScoreScore = 0.7f;fMinReliability = 0.5f;fTaoCenters = 0.05f;//0.05 	
		MinOrientedRectSide=3.0f;
		main_CNC();
		}break;
	default: break;
	}
	system("pause");
	//main_normalDB(argc, argv);//论文中三个数据库，不同参数的实验
	//OnVideo();//摄像头
	return 0;	   
}

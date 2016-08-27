// stdafx.h : 标准系统包含文件的包含文件，
// 或是经常使用但不常更改的
// 特定于项目的包含文件
//
#pragma once
#include <opencv2/opencv.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/core/types_c.h>

#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/features2d/features2d.hpp>

#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

#include <numeric>
#include <unordered_map>
#include <vector>

#include <time.h>
#include <iostream>
#include <fstream>
#include <math.h>

#include <direct.h>

#include <sys/stat.h> 

#ifdef LINUX
#include<dirent.h>
#include </usr/include/x86_64-linux-gnu/sys/io.h>
#else 
#include <io.h>
#include "dirent.h"
#endif

using namespace std;using namespace cv;

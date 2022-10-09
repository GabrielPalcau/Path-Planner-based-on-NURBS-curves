// ******************************
// 
// Obstacle detection algorithm (C++)
// 
// @ Copyright 2022, Gabriel Lucian Palcau
//
// ******************************

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <iostream>
#include <fstream>

using namespace cv;
using namespace std;

const string pathFile[] = {"C:/Users/palca/Desktop/FileProgrammiTesi/Obstacle_Position.txt"};
const string path = "C:/Users/palca/Desktop/FileProgrammiTesi/Ostacoli1.bmp";

int heightInitial = 0;
int widthInitial = 0;
const int heightTarget = 600;
const int widthTarget = 1000;

void getContours(const Mat& imgCanny, const Mat& img, const Mat& img_copy);

int main()
{
    Mat img, imgCanny, imgDil, imgMod, imgGrid;
    img = imread(path, IMREAD_COLOR);

    // Check for failure
    if (img.empty())
    {
        cout << "Could not open or find the image" << endl;
        cin.get(); //wait for any key press
        return -1;
    }
    resize(img, img, Size(600,400), INTER_LINEAR);
    heightInitial = img.rows;
    widthInitial = img.cols;


    imshow("Original Image", img);

    
    img.copyTo(imgGrid);
    cvtColor(img, imgMod, COLOR_BGR2GRAY);

    GaussianBlur(imgMod, imgMod, Size(5, 5), 25, 25);
    Canny(imgMod, imgCanny, 50, 160);

    Mat kernel = getStructuringElement(MORPH_RECT, Size(3, 3));
    dilate(imgCanny, imgDil, kernel);

    getContours(imgDil, img, imgGrid);

    imwrite("C:/Users/palca/Desktop/FileImmaginiTesi/immagineOBS.png", img);

    imshow("Detected Obstacles", img);
    imshow("GaussianBlur Image", imgMod);
    imshow("Canny Detection", imgCanny);
    imshow("DilateCanny", imgDil);
    imshow("GridMap 0/1 Obstacle/Non-Obstacle", imgGrid);
    waitKey(0);

    return 0;
}

void getContours(const Mat& imgCanny, const Mat& img, const Mat& img_copy) {

    vector<vector<Point>> contours;
    vector<Vec4i> hierarchy;

    findContours(imgCanny, contours, hierarchy, RETR_EXTERNAL, CHAIN_APPROX_NONE);
    drawContours(img, contours, -1, Scalar(0, 255, 0), 2);

    ofstream File;
    File.open(pathFile[0], ios::out);

    //Fill contours
    vector< vector<Point> > hull(contours.size());

    for (int i = 0; i < contours.size(); i++)
        convexHull(Mat(contours[i]), hull[i], false);

    fillPoly(img_copy, hull, Scalar(0, 0, 0));

    //Center computation
    vector<Point> centers;
    for (int i = 0; i < contours.size(); i++) {
        Moments M = moments(contours[i]);
        Point center(M.m10 / M.m00, M.m01 / M.m00);
        centers.push_back(center);
        circle(img, centers[i], 1, CV_RGB(255, 0, 0), 5);
    }


    //Max distance from center
    float max_dist = 0.0;
    float max_dist1 = 0.0;
    float dist = 0.0;
    float dist1 = 0.0;

    vector<Point> perimeter;
    for (int i = 0; i < contours.size(); i++) {
        perimeter = contours[i];
        int x = 0;
        int y = 0;
        max_dist = 0.0;

        for (int j = 0; j < size(perimeter); j++) {
            x = perimeter[j].x;
            y = perimeter[j].y;

            dist = sqrt(pow(x - centers[i].x, 2) + pow(y - centers[i].y, 2));
            dist1 = sqrt(pow(x * widthTarget / widthInitial - centers[i].x * widthTarget / widthInitial, 2) + pow(y * heightTarget / heightInitial - centers[i].y * heightTarget / heightInitial, 2));
            if (dist > max_dist) {
                max_dist = dist;
                max_dist1 = dist1;
            }
        }

        circle(img, centers[i], max_dist, CV_RGB(0, 0, 255), 2);
        File << (int)ceil(centers[i].x * widthTarget / widthInitial) << " " << (int)ceil((heightInitial - centers[i].y) * heightTarget / heightInitial) << " " << (int)ceil(max_dist1) << endl;

    }

    File.close();

}


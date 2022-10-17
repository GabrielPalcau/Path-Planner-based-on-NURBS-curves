// ******************************
// 
// NURBS path planner algorithm (C++)
// 
// @ Copyright 2022, Gabriel Lucian Palcau
//
// ******************************

#include <iostream>
#include "generic_curve.hpp"
#include "curve.hpp"
#include "circular_arc.hpp"
#include <vector>
#include "sisl.h"
#include "straight_line.hpp"
#include <fstream>

#ifdef _WIN32
#include <GL/glut.h>
#elif defined(_APPLE_)
#endif

#define CURVE_EVALUATIONS 500
using namespace std;

bool flag3d = false;
static int wayPointsFlag = 0;
int cntFile = 0;

// ------------------Global Variable ------------------

shared_ptr<GenericCurve> genericCurve = nullptr;
vector<Eigen::Vector3d> ctrlPoints{};

vector<shared_ptr<CircularArc>> obs;
vector<shared_ptr<CircularArc>> obs_external;
vector<Eigen::Vector3d> obs_points{};

vector<Eigen::Vector3d> pathRobot{};
vector<Eigen::Vector3d> wayPoints{};

vector<Eigen::Vector3d> pathRobotGlobal{};

int RobotPosition = 0;
int initSizePulse = 5;
int statePulse = 1;

static Eigen::Vector3d START{0,0,0};
static Eigen::Vector3d GOAL{0,0,0};

// -- OpenGL drawing prototype functions --
void DisplayCtrPoints();
void init2d();
void DrawSinglePoints(const vector<Eigen::Vector3d>& points);
void DrawObs(const shared_ptr<CircularArc>& ob, const bool flag);
void Draw2dNURBS();
void DrawpathRobotGlobal();

// -- Prototype Function for Manimulation of the SISL curve --

void UpdateNurbs(int& degree, vector<double>& knots, const vector<double>& weights);
Eigen::Vector4d firstIntersection(vector<Eigen::Vector3d>& intersections, const shared_ptr<GenericCurve>& nurbs, const bool flag);
int FindObstacleFromPoint(const shared_ptr<CircularArc>& pc1, const Eigen::Vector3d& point);
int UpdateParameterOrder(vector<double>& paramOrder);
double UpdateParameterDistance(const Eigen::Vector3d& Point, const shared_ptr<GenericCurve>& nurbs);
Eigen::Vector3d FindCenterMedian(const shared_ptr<CircularArc>& C1, const shared_ptr<CircularArc>& C2);
vector<Eigen::Vector3d> ExternalIntersection(const int pos_first, const vector<int>& VecCenterIntersection, const Eigen::Vector3d& Last, const Eigen::Vector3d& LastL, const int cnt, const int Median, const bool Median1);
void printNurbsInformation(const vector<double>& knots);
void pathRobotCalculation();
int EvaluatePointInObstacle(const Eigen::Vector3d& point);
void ExternalRadiusOpt(vector<int>& radius);
void CalculateWayPoints(vector<Eigen::Vector3d>& wayPointsClosest);
void timer(int);
void timer1(int);
void timer2(int);
void display();
int pathgeneration();





// --------------- MAIN ------------------

int main(int argc, char** argv) {

    //Acquiring Start and Goal position
    cout << "******* Insert Start and Goal position *******" << endl;
    cout << "Insert Start position, x and y: " << endl;
    cin >> START.x();
    cin >> START.y();
    cout << "Insert Goal position, x and y: " << endl;
    cin >> GOAL.x();
    cin >> GOAL.y();

    //Trajectory printed on OpenGL
    int x = pathgeneration();
    if (x == -10){
        cerr << "Error, too complex environment";
        return -1;
    }

    glutInit(&argc, argv);
    glutInitWindowSize(1000, 600);
    glutInitWindowPosition(175, 40);
    glutCreateWindow("Trajectory");
    init2d();

    glutDisplayFunc(display);
    glutTimerFunc(0,timer,0);
    glutTimerFunc(0,timer1,0);
    glutTimerFunc(0,timer2,0);
    glDisable(GL_TEXTURE_2D);
    glutMainLoop();

    return 0;
}

// ------ OpenGL drawing functions -------

void DisplayCtrPoints() 
{
    glPointSize(5.0);
    glColor3f(0.0f, 0.0f, 1.0f);

    glBegin(GL_POINTS);
    for (const auto& coord : ctrlPoints)
        glVertex2f((GLfloat)coord[0], (GLfloat)coord[1]);
    glEnd();
}

void init2d()
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0.0, 1000, 0.0, 600);
    glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
    glClear(GL_COLOR_BUFFER_BIT);
}

void DrawSinglePoints(const vector<Eigen::Vector3d>& points)
{
    glPointSize(5.0);
    glColor3f(0.0f, 0.0f, 0.0f);

    glBegin(GL_POINTS);
    for (const auto& i : points)
        glVertex2f((float)i[0], (float)i[1]);
    glEnd();
}

void DrawObs(const shared_ptr<CircularArc>& ob, const bool flag)
{
    int left = 0;
    auto* discrete = new double[3 * CURVE_EVALUATIONS];

    for (int j = 0; j < CURVE_EVALUATIONS; j++)
    {
        double t = ob->StartParameter_s() + (ob->EndParameter_s() - ob->StartParameter_s()) * j / (CURVE_EVALUATIONS - 1.0);
        int stat;
        s1227(ob->CurvePtr(), 0, t, &left, discrete + 3 * j, &stat);

        if (stat != 0)
            cerr << "s1227 returned status: " << stat << endl;
    }

    if (flag3d) {
        glLineWidth(2.0);
        glBegin(GL_LINE_STRIP);
        for (int j = 0; j < CURVE_EVALUATIONS; j++) {
            glColor3f(0.0f, 0.5f, 1.0f);
            glVertex3f((float)discrete[0 + (3 * j)] * 10 / 1000, (float)discrete[1 + (3 * j)] * 6 / 600,
                (float)discrete[2 + (3 * j)] * 3 / 100);
        }
        glEnd();
        glFlush();
    }
    else {
        glLineWidth(1.0);
        glBegin(GL_LINE_STRIP);
        for (int j = 0; j < CURVE_EVALUATIONS; j++) {
            glColor3f(0.0f, 0.5f, 1.0f);

            if (flag)
                glVertex2f((float)ob->CentrePoint().x(), (float)ob->CentrePoint().y());
            glVertex2f((float)discrete[0 + (3 * j)], (float)discrete[1 + (3 * j)]);
        }
        glEnd();
    }
    delete[] discrete;
}

void Draw2dNURBS()
{
    glLineWidth(2.0);
    glBegin(GL_LINE_STRIP);
    glColor3f(1.0f, 0.0f, 0.0f);
    for (const auto& j : pathRobot)
        glVertex2f((float)j.x(), (float)j.y());

    glEnd();
}

void DrawpathRobotGlobal()
{
    glLineWidth(2.0);
    glColor3f(0.0f, 0.0f, 1.0f);
    glBegin(GL_LINE_STRIP);
    for (const auto& ele : pathRobotGlobal)
        glVertex2f((GLfloat)ele.x(), (GLfloat)ele.y());
    
    glEnd();
}

// ----- Function for Manimulation of the SISL curve -----

void UpdateNurbs(int& degree, vector<double>& knots, const vector<double>& weights)
{
    unsigned int nCtrl = ctrlPoints.size();
    if (ctrlPoints.size() == 3)
        degree = 2;
    else if (ctrlPoints.size() >= 4)
        degree = 3;

    unsigned int n_knots = degree + 1 + nCtrl; //Mathematical Formula for NURBS

    int value = 0;
    knots.clear();
    for (int i = 0; i < n_knots; i++) {
        if (i < degree + 1) {
            knots.push_back(value);
        }
        else {
            if ((n_knots - i) > degree + 1) {
                value++;
                knots.push_back(value);
            }
            else {
                value++;
                knots.push_back(value);
                value--;
            }
        }
    }
}

Eigen::Vector4d firstIntersection(vector<Eigen::Vector3d>& intersections, const shared_ptr<GenericCurve>& nurbs, const bool flag)
{
    double min = 9999;
    double eps = 5;
    Eigen::Vector3d temp;
    Eigen::Vector4d temp1 = { 0,0,0,0 };
    shared_ptr<StraightLine> ptr_Line = nullptr;

    for (int j = 0; j < 10 * CURVE_EVALUATIONS; j++) {

        double t = nurbs->StartParameter_s() + (nurbs->EndParameter_s() - nurbs->StartParameter_s()) * j / (10 * CURVE_EVALUATIONS - 1.0);
        nurbs->FromAbsSislToPos(t, temp);

        if (!flag) {
            for (const auto& ele : intersections) {
                ptr_Line = make_shared<StraightLine>(temp, ele, 3, 2);
                if (ptr_Line->Length() <= eps) {
                    auto elef = ele;
                    auto it = remove(intersections.begin(), intersections.end(), ele);
                    intersections.erase(intersections.end());
                    return { elef.x(), elef.y(), elef.z(), t };
                }
            }
        }
        else {

            tuple<double, double> closest = nurbs->FindClosestPoint(intersections[0]);
            auto abscissa_m1 = (double)(get<0>(closest));
            auto ClosestPoint = nurbs->At(double(abscissa_m1));
            ptr_Line = make_shared<StraightLine>(ClosestPoint, temp, 3, 2);

            if (ptr_Line->Length() <= min) {
                min = ptr_Line->Length();
                temp1 = { ClosestPoint.x(),ClosestPoint.y(),ClosestPoint.z(),t };
            }
        }
    }
    return temp1;
}

int FindObstacleFromPoint(const shared_ptr<CircularArc>& pc1, const Eigen::Vector3d& point)
{
    bool inside = false;
    Eigen::Vector3d temp = { point.x(),point.y(),point.z() };
    shared_ptr<StraightLine> line = nullptr;
    line = make_shared<StraightLine>(pc1->CentrePoint(), temp, 3, 2);
    shared_ptr<StraightLine> line1 = nullptr;
    line1 = make_shared<StraightLine>(pc1->CentrePoint(), pc1->StartPoint(), 3, 2);

    if (line->Length() <= (line1->Length() + 1))
        inside = true;
    return inside;
}

int UpdateParameterOrder(vector<double>& paramOrder)
{
    auto itP = paramOrder.begin();
    bool flagP = false;
    int index = -1;
    double x = -1;

    if (paramOrder.size() == 1)
        index = 1;
    else {
        x = paramOrder.back();
        paramOrder.pop_back();
        for (int p = 0; p < paramOrder.size(); p++) {
            if (x <= paramOrder.at(p)) {
                flagP = true;
                index = p + 1;
                paramOrder.insert(itP + p, x);
                break;
            }
        }
        if (!flagP) {
            paramOrder.push_back(x);
            index = (int)paramOrder.size();
        }
    }
    return index;
}

double UpdateParameterDistance(const Eigen::Vector3d& Point, const shared_ptr<GenericCurve>& nurbs)
{
    float eps = 1;
    Eigen::Vector3d temp;
    Eigen::Vector3d cPoint = Point;
    shared_ptr<StraightLine> Line = nullptr;

    tuple<double, double> closest = nurbs->FindClosestPoint(cPoint);
    auto abscissa_m1 = (double)(get<0>(closest));
    auto ClosestPoint = nurbs->At(double(abscissa_m1));

    for (int j = 0; j < 5 * CURVE_EVALUATIONS; j++) {

        double t = nurbs->StartParameter_s() + (nurbs->EndParameter_s() - nurbs->StartParameter_s()) * j / (5 * CURVE_EVALUATIONS - 1.0);
        nurbs->FromAbsSislToPos(t, temp);
        Line = make_shared<StraightLine>(ClosestPoint, temp, 3, 2);

        if (Line->Length() <= eps)
            return t;
    }
    return -1;
}

Eigen::Vector3d FindCenterMedian(const shared_ptr<CircularArc>& C1, const shared_ptr<CircularArc>& C2)
{
    shared_ptr<StraightLine> ptr_Line = nullptr;
    ptr_Line = make_shared<StraightLine>(C1->CentrePoint(), C2->CentrePoint(), 3, 2);

    vector<Eigen::Vector3d> temp{};
    vector<Eigen::Vector3d> intersections{};

    temp = ptr_Line->Intersection(C1);
    intersections.push_back(temp[0]);

    temp = ptr_Line->Intersection(C2);
    intersections.push_back(temp[0]);

    Eigen::Vector3d Medium{};
    Medium.x() = (intersections[0].x() + intersections[1].x()) / 2;
    Medium.y() = (intersections[0].y() + intersections[1].y()) / 2;
    Medium.z() = 0;
    return Medium;
}

vector<Eigen::Vector3d> ExternalIntersection(const int pos_first, const vector<int>& VecCenterIntersection, const Eigen::Vector3d& Last, const Eigen::Vector3d& LastL, const int cnt, const int Median, const bool Median1)
{
    const vector<shared_ptr<CircularArc>>& circle = obs_external;
    vector<Eigen::Vector3d> Intersections{};

    if (cnt == 2 && !Median1) {
        for (const int& ele : VecCenterIntersection) {
            auto IntersectionsTemp = circle.at(pos_first)->Intersection(circle.at(ele));
            for (const auto& ele1 : IntersectionsTemp)
                Intersections.push_back(ele1);
        }
    }
    else if (cnt == 1 && Median1) {

        auto IntersectionsTemp = circle.at(pos_first)->Intersection(circle.at(Median));
        for (const auto& ele1 : IntersectionsTemp)
            Intersections.push_back(ele1);
        IntersectionsTemp = circle.at(VecCenterIntersection.at(0))->Intersection(circle.at(VecCenterIntersection.at(1)));
        for (const auto& ele1 : IntersectionsTemp)
            Intersections.push_back(ele1);
    }
    else {
        Intersections = circle[pos_first]->Intersection(circle[Median]);
    }

    shared_ptr<StraightLine> Line = nullptr;
    double min = 9999;
    double max = 9999;
    Eigen::Vector3d MinPoint;
    Eigen::Vector3d MaxPoint;

    for (const auto& ele : Intersections) {
        Line = make_shared<StraightLine>(Last, ele, 3, 2);
        if (Line->Length() <= min) {
            min = Line->Length();
            MinPoint = ele;
        }

        Line = make_shared<StraightLine>(LastL, ele, 3, 2);
        if (Line->Length() <= max) {
            max = Line->Length();
            MaxPoint = ele;
        }
    }

    vector<Eigen::Vector3d> TwoPoints{};
    TwoPoints.push_back(MinPoint);
    TwoPoints.push_back(MaxPoint);
    return TwoPoints;
}

void printNurbsInformation(const vector<double>& knots)
{
    cout << "Generated NURBS:" << endl;
    cout << "----- Degree: " << genericCurve->Degree() << endl;
    cout << "----- Control Points: " << endl;
    cout << "----- ----- [ ";
    for (const auto& ele : ctrlPoints)
        cout << ele.x() << " ";
    cout << endl << "              ";
    for (const auto& ele : ctrlPoints)
        cout << ele.y() << " ";
    cout << "]" << endl;
    cout << "----- Knot Vector: " << endl;
    cout << "----- ----- [ ";
    for (const auto& ele : knots)
        cout << ele << " ";
    cout << "]" << endl;
}

void pathRobotCalculation()
{
    pathRobot.clear();

    int left = 0;
    auto* discrete = new double[3 * 10 * CURVE_EVALUATIONS];
    double len = 0.0;
    double lenTemp = 0.0;
    shared_ptr<StraightLine> line = nullptr;
    Eigen::Vector3d Apoint;
    Eigen::Vector3d Bpoint;

    for (int j = 0; j < 10 * CURVE_EVALUATIONS; j++) {
        double t = genericCurve->StartParameter_s() +
            (genericCurve->EndParameter_s() - genericCurve->StartParameter_s()) * j / (10 * CURVE_EVALUATIONS - 1.0);

        int stat;
        s1227(genericCurve->CurvePtr(), 0, t, &left, discrete + 3 * j, &stat);
        if (stat != 0)
            cerr << "s1227 returned status: " << stat << endl;

        if (j == 1) {
            Bpoint = { discrete[0 + (3 * j)], discrete[1 + (3 * j)], 0.0 };
            Apoint = { discrete[0 + (3 * (j - 1))], discrete[1 + (3 * (j - 1))], 0.0 };
            pathRobot.push_back(Apoint);
            pathRobot.push_back(Bpoint);
            len = 1;
        }

        if (j > 1)
        {
            Bpoint = { discrete[0 + (3 * j)], discrete[1 + (3 * j)], 0.0 };
            line = make_shared<StraightLine>(pathRobot.back(), Bpoint, 3, 2);
            lenTemp = line->Length();
            if (lenTemp >= (len - len / 8)) {
                pathRobot.push_back(Bpoint);
            }
        }
    }
}

int EvaluatePointInObstacle(const Eigen::Vector3d& point)
{
    int ObsToPoint = 0;
    for (int z = 0; z < obs.size(); z++) {
        if (FindObstacleFromPoint(obs.at(z), point)) {
            ObsToPoint = z;
            break;
        }
    }
}

void ExternalRadiusOpt(vector<int>& radius)
{
    vector<Eigen::Vector3d> intersections{};
    int times = 0;
    Eigen::Vector3d axis = { 0,0,1 };
    bool limit = false;

    for (int i = 0; i < obs_external.size(); i++) {
        times = 0;
        limit = false;
        for (const auto& ele : obs) {
            intersections = obs_external.at(i)->Intersection(ele);
            if (!intersections.empty())
                times++;
        }
        if (times > 2) {
            radius.at(i) = radius.at(i) - 1;
            if (radius.at(i) <= (obs.at(i)->StartPoint().x() - obs.at(i)->CentrePoint().x()) + 3) {
                radius.at(i) = (obs.at(i)->StartPoint().x() - obs.at(i)->CentrePoint().x()) + 5;
                limit = true;
            }
            obs_external.at(i) = make_shared<CircularArc>(2 * 3.1415926, axis, obs_points.at(i) + Eigen::Vector3d{ (double)radius.at(i), 0.0, 0.0 }, obs_points.at(i), 3, 3);
            if (!limit)
                i--;
        }
    }
}

void CalculateWayPoints(vector<Eigen::Vector3d>& wayPointsClosest)
{
    shared_ptr<StraightLine> line = nullptr;
    wayPoints.clear();

    for (auto& j : wayPointsClosest) {
        for (int i = 0; i < pathRobot.size(); i++) {
            line = make_shared<StraightLine>(j, pathRobot.at(i), 3, 2);
            if (line->Length() < 10) {
                wayPoints.push_back(pathRobot.at(i + 30));
                break;
            }
        }
    }
}

void timer(int)
{
    glutPostRedisplay();
    glutTimerFunc(1000 / 200, timer, 0);

    if (wayPointsFlag == 0) {
        for (const auto& ele : wayPoints) {
            if (ele == pathRobot[RobotPosition]) {
                wayPointsFlag = 1;
                pathRobotGlobal.push_back(pathRobot.at(RobotPosition));
            }
        }
    }

    if (wayPointsFlag == 0) {
        pathRobotGlobal.push_back(pathRobot.at(RobotPosition));
        RobotPosition++;
    }

    if (RobotPosition >= pathRobot.size())
        RobotPosition = (int)pathRobot.size() - 1;
}

void timer1(int)
{
    glutPostRedisplay();
    glutTimerFunc(1000 / 10, timer1, 0);
    switch (statePulse)
    {
    case 1:
        if (initSizePulse <= 10)
            initSizePulse++;
        else
            statePulse = 2;
        break;
    case 2:
        if (initSizePulse >= 5)
            initSizePulse--;
        else
            statePulse = 1;
        break;
    default:
        cerr << "Error in switching state Pulse of the robot position" << endl;
    }
}

void timer2(int)
{
    glutTimerFunc(3000, timer2, 0);

    if (wayPointsFlag == 1)
    {
        int x = pathgeneration();
        wayPointsFlag = 0;
        if (x == -10) {
            cerr << "Error too complex environment";
        }
    }
}

void display()
{
    glClear(GL_COLOR_BUFFER_BIT);

    for (const auto& ob : obs) {
        DrawObs(ob, true);
    }
    for (const auto& ob : obs_external) {
        DrawObs(ob, false);
    }

    Draw2dNURBS();
    DisplayCtrPoints();
    DrawSinglePoints(wayPoints);
    DrawpathRobotGlobal();

    glPointSize((float)initSizePulse);
    glColor3f(1.0f, 0.0f, 0.0f);
    glBegin(GL_POINTS);
    glVertex2f((GLfloat)pathRobot[RobotPosition].x(), (GLfloat)pathRobot[RobotPosition].y());
    glEnd();

    glutSwapBuffers();
}

int pathgeneration()
{

    cntFile++;

    // Defining the path of the obstacles txt
    string pathObstacleTxT;
    try {
        if (cntFile >= 1)
            pathObstacleTxT = "C:/Users/palca/Desktop/FileProgrammiTesi/Obstacle_Position.txt";
        else if (cntFile > 1)
            pathObstacleTxT = "C:/Users/palca/Desktop/FileProgrammiTesi/Obstacle_Position1.txt";
    }
    catch (...) {
        cerr << "Error while trying to read the path where the obstacle position are written !" << endl;
        return -1;
    }

    //Pre-Processing for saving the obstacle position from file (control Points)
    static vector<int> obs_radius{};
    static vector<Eigen::Vector3d> ctrlPoints_All;

    if (cntFile == 1)
        ctrlPoints_All.emplace_back(START);

    bool flagDynamicPath = false;
    bool flagSameObs = false;
    int x, y, r;
    Eigen::Vector3d ctrlTemp{};
    stringstream ss;
    string line;
    ifstream dataIN;

    dataIN.open(pathObstacleTxT);
    if (dataIN.is_open())
    {
        while (getline(dataIN, line))
        {
            flagSameObs = false;
            ss.str(line);
            ss >> x;
            ss >> y;
            ss >> r;

            if (r > 5) {

                if (cntFile == 1) {
                    ctrlPoints_All.emplace_back(x, y, 0);
                    obs_radius.push_back(r);
                }
                else {
                    for (const auto& ele : ctrlPoints_All) {
                        if ((x == (int)ele.x()) && (y == (int)ele.y())) {
                            flagSameObs = true;
                        }
                    }
                    if (!flagSameObs) {
                        flagDynamicPath = true;
                        ctrlPoints_All.insert(ctrlPoints_All.begin() + 1, Eigen::Vector3d{ (double)x,(double)y,0 });
                        obs_radius.insert(obs_radius.begin() + 1, r);
                    }
                }
            }
            ss.clear();
        }
        dataIN.close();
    }
    else {
        cerr << "Unable to open file txt where the position of the obstacles are saved!" << endl;
        return -2;
    }

    if (cntFile == 1)
        ctrlPoints_All.emplace_back(GOAL);

    // If new obstacles are not detected, return to the original path
    if (!flagDynamicPath && cntFile > 1) {
        wayPoints.erase(wayPoints.begin());
        return 0;
    }

    // If new obstacles are found, go to dynamic path
    if (wayPointsFlag == 1 && cntFile > 1)
    {
        ctrlPoints_All.at(0) = wayPoints.at(0);
        wayPoints.erase(wayPoints.begin());
        RobotPosition = 0;

        obs.clear();
        obs_external.clear();
        obs_points.clear();
    }

    const int nAll = (int)ctrlPoints_All.size();
    const int nObs = nAll - 2;
    auto itAll = ctrlPoints_All.begin();

    // Defining the FISRT parameters for the NURBS, from start to goal, it will be replaced later
    const int dimension = 3;
    int degree = 1;

    ctrlPoints.clear();
    ctrlPoints.push_back(*itAll);
    ctrlPoints.push_back(*(ctrlPoints_All.end() - 1));
    auto it = ctrlPoints.begin();

    vector<double> weights{};
    weights.push_back(1.0);
    weights.push_back(1.0);
    auto it1 = weights.begin();

    vector<double> knots{};
    const vector<double> coeff{};

    // ---- Defining only the obstacles control points ----

    for (int i = 1; i < ctrlPoints_All.size() - 1; i++)
        obs_points.emplace_back(ctrlPoints_All.at(i));

    // ---- Defining the SISL curve for obstacle and external_obstacle circle ----

    const int safetyDim = 15;
    vector<int> obs_ext_radius{};

    for (int i = 0; i < nObs; i++)
        obs_ext_radius.push_back(obs_radius.at(i) + safetyDim);

    const int dimensionObs = 3;
    const int orderObs = 3;
    const Eigen::Vector3d axis = { 0.0, 0.0, 1.0 };

    for (int i = 0; i < nObs; i++) {
        obs.push_back(make_shared<CircularArc>(2 * 3.1415926, axis, obs_points.at(i) + Eigen::Vector3d{ (double)obs_radius.at(i), 0.0, 0.0 }, obs_points.at(i), dimensionObs, orderObs));
        obs_external.push_back(make_shared<CircularArc>(2 * 3.1415926, axis, obs_points.at(i) + Eigen::Vector3d{ (double)obs_ext_radius.at(i), 0.0, 0.0 }, obs_points.at(i), dimensionObs, orderObs));
    }

    // Function for External radius optimization, maximum three obstacles intersections
    ExternalRadiusOpt(obs_ext_radius);


    // ----Definition of variables ----

    bool flagEnd = false;

    vector<Eigen::Vector3d> intersections_All = {};
    vector<Eigen::Vector3d> intersections{};
    vector<int> position{};
    vector<double> paramOrder;
    vector<int> visited;
    vector<Eigen::Vector3d> freeControlObs;

    Eigen::Vector3d ClosestPoint = {};
    Eigen::Vector3d second{};
    Eigen::Vector3d first{};

    int c = 0;
    int pos_first = -1;
    int previous = -2;
    int index = -1;
    int VectorAdder = 0;

    shared_ptr<StraightLine> ptr_straightLine = nullptr;
    vector<Eigen::Vector3d> wayPointsClosest{};

    // Cicle for the NURBS path generation

    for (int k = 0; k < 20; k++) {

        cout << endl << endl;
        cout << "Evaluating intersection of obstacle: iteration n." << k << endl;

        UpdateNurbs(degree, knots, weights);
        it = ctrlPoints.begin();
        it1 = weights.begin();

        genericCurve = make_shared<GenericCurve>(degree, knots, ctrlPoints, weights, coeff, dimension, degree + 1);

        intersections_All.clear();
        position.clear();
        pos_first = -1;

        //---- Calculating all the NURBS points intersection ----
        for (int j = 0; j < obs_points.size(); j++) {

            intersections.clear();
            intersections = genericCurve->Intersection(obs[j]);

            if (!intersections.empty())
                position.push_back(j);

            for (const auto& ele : intersections)
                intersections_All.push_back(ele);
        }

        if (!intersections_All.empty()) {

            //Find the first intersection, the obstacle, and the exit point from the obstacle
            Eigen::Vector4d temp1 = firstIntersection(intersections_All, genericCurve, false);
            first = { temp1.x(),temp1.y(), temp1.z() };
            paramOrder.push_back(temp1.w());
            auto itP = paramOrder.begin();
            index = UpdateParameterOrder(paramOrder);

            for (int w = 0; w < obs_points.size(); w++) {
                if (FindObstacleFromPoint(obs.at(w), first)) {
                    pos_first = w;
                    break;
                }
            }

            temp1 = firstIntersection(intersections_All, genericCurve, false);
            second = { temp1.x(),temp1.y(), temp1.z() };

            if (!FindObstacleFromPoint(obs.at(pos_first), second)) {
                second = {};
                second = first;
            }

            // Evaluate if the first obstacle (external circle) intersects others obstacles
            auto itV = find(visited.begin(), visited.end(), pos_first);
            vector<int> VecCenterIntersection{};
            for (int q = 0; q < obs_points.size(); q++) {
                if (pos_first != q) {
                    auto CenterIntersection = obs_external[pos_first]->Intersection(obs[q]);
                    if (!CenterIntersection.empty())
                        VecCenterIntersection.push_back(q);
                }
            }

            if (itV != visited.end()) {

                // --- We already visited this obstacle ---
                cout << "Obstacle already visited.." << endl;

                ptr_straightLine = nullptr;

                // Calculation of the point to be added on the external radius for avoiding the obstacle
                double x_abs = abs(first.x() - obs_points[pos_first].x());
                double y_abs = abs(first.y() - obs_points[pos_first].y());
                double kk = obs_ext_radius[pos_first] / min(x_abs, y_abs);

                if (first.y() <= obs_points[pos_first].y()) {
                    if (first.x() <= obs_points[pos_first].x()) {
                        ptr_straightLine = make_shared<StraightLine>(obs_points[pos_first], obs_points[pos_first] + Eigen::Vector3d{ -kk * x_abs, -kk * y_abs, 0.0 }, 3, 2);
                    }
                    else {
                        ptr_straightLine = make_shared<StraightLine>(obs_points[pos_first], obs_points[pos_first] + Eigen::Vector3d{ +kk * x_abs, -kk * y_abs, 0.0 }, 3, 2);
                    }
                }
                else {
                    if (first.x() <= obs_points[pos_first].x()) {
                        ptr_straightLine = make_shared<StraightLine>(obs_points[pos_first], obs_points[pos_first] + Eigen::Vector3d{ -kk * x_abs,kk * y_abs,0.0 }, 3, 2);
                    }
                    else {
                        ptr_straightLine = make_shared<StraightLine>(obs_points[pos_first], obs_points[pos_first] + Eigen::Vector3d{ kk * x_abs,kk * y_abs,0.0 }, 3, 2);
                    }
                }

                intersections.clear();
                intersections = ptr_straightLine->Intersection(obs_external[pos_first]);

                ctrlPoints.insert(it + index, intersections.at(0));
                weights.insert(it1 + index, 1.4);
                UpdateNurbs(degree, knots, weights);
                genericCurve = std::make_shared<GenericCurve>(degree, knots, ctrlPoints, weights, coeff, dimension, degree + 1);
                for (int s = 0; s < paramOrder.size(); s++)
                    paramOrder[s] = UpdateParameterDistance(ctrlPoints[s + 1], genericCurve);

                wayPointsClosest.pop_back();
                tuple<double, double> closest = genericCurve->FindClosestPoint(ctrlPoints[ctrlPoints.size() - 2]);
                auto abscissa_m1 = (double)(get<0>(closest));
                ClosestPoint = genericCurve->At(double(abscissa_m1));
                wayPointsClosest.push_back(ClosestPoint);

                printNurbsInformation(knots);
                continue;

            }
            else {

                //---- First time we visit this obstacle ----
                cout << "Obstacle visited for the first time.. " << endl;

                visited.push_back(pos_first);

                ctrlPoints.insert(it + index, first);
                weights.insert(it1 + index, 0.8);
                UpdateNurbs(degree, knots, weights);
                genericCurve = std::make_shared<GenericCurve>(degree, knots, ctrlPoints, weights, coeff, dimension, degree + 1);
                for (int s = 0; s < paramOrder.size(); s++)
                    paramOrder[s] = UpdateParameterDistance(ctrlPoints[s + 1], genericCurve);

                it = ctrlPoints.begin();
                it1 = weights.begin();

                // Calculate external point to be addend to bring the NURBS away from the obstacle
                tuple<double, double> closest = genericCurve->FindClosestPoint(obs_points[pos_first]);
                auto abscissa_m1 = (double)(get<0>(closest));
                ClosestPoint = genericCurve->At(double(abscissa_m1));
                double x_abs = abs(ClosestPoint.x() - obs_points[pos_first].x());
                double y_abs = abs(ClosestPoint.y() - obs_points[pos_first].y());
                double kk = obs_ext_radius[pos_first] / min(x_abs, y_abs);
                kk = kk + 1;

                Eigen::Vector3d finalPoint1;

                if (x_abs < 0.001 || y_abs < 0.001) {
                    if (x_abs < 0.001) {
                        if (ClosestPoint.y() <= obs_points[pos_first].y())
                            finalPoint1 = obs_points[pos_first] + Eigen::Vector3d{ 0.0, -(double)obs_ext_radius.at(pos_first) - 10, 0.0 };
                        else
                            finalPoint1 = obs_points[pos_first] + Eigen::Vector3d{ 0.0, +(double)obs_ext_radius.at(pos_first) + 10, 0.0 };
                    }
                    if (y_abs < 0.001) {
                        if (ClosestPoint.x() <= obs_points[pos_first].x())
                            finalPoint1 =
                            obs_points[pos_first] + Eigen::Vector3d{ -(double)obs_ext_radius.at(pos_first) - 10, 0.0, 0.0 };
                        else
                            finalPoint1 =
                            obs_points[pos_first] + Eigen::Vector3d{ +(double)obs_ext_radius.at(pos_first) + 10, 0.0, 0.0 };
                    }
                }
                else {
                    if (ClosestPoint.y() <= obs_points[pos_first].y()) {
                        if (ClosestPoint.x() <= obs_points[pos_first].x()) {
                            finalPoint1 = obs_points[pos_first] + Eigen::Vector3d{ -kk * x_abs, -kk * y_abs, 0.0 };
                        }
                        else {
                            finalPoint1 = obs_points[pos_first] + Eigen::Vector3d{ +kk * x_abs, -kk * y_abs, 0.0 };
                        }
                    }
                    else {
                        if (ClosestPoint.x() <= obs_points[pos_first].x()) {
                            finalPoint1 = obs_points[pos_first] + Eigen::Vector3d{ -kk * x_abs, +kk * y_abs, 0.0 };
                        }
                        else {
                            finalPoint1 = obs_points[pos_first] + Eigen::Vector3d{ +kk * x_abs, +kk * y_abs, 0.0 };
                        }
                    }
                }

                ptr_straightLine = nullptr;
                ptr_straightLine = make_shared<StraightLine>(obs_points[pos_first], finalPoint1, 3, 2);
                freeControlObs = ptr_straightLine->Intersection(obs_external[pos_first]);

                // Evaluating in the freeControlObs is inside another obstacle
                if (EvaluatePointInObstacle(freeControlObs.at(0)) == 0) {

                    intersections.clear();
                    intersections.push_back(freeControlObs[0]);
                    temp1 = firstIntersection(intersections, genericCurve, true);
                    paramOrder.push_back(temp1.w());
                    intersections.clear();

                    index = UpdateParameterOrder(paramOrder);
                    ctrlPoints.insert(it + index, freeControlObs[0]);
                    weights.insert(it1 + index, 2);
                    UpdateNurbs(degree, knots, weights);
                    genericCurve = std::make_shared<GenericCurve>(degree, knots, ctrlPoints, weights, coeff, dimension, degree + 1);
                    for (int s = 0; s < paramOrder.size(); s++)
                        paramOrder[s] = UpdateParameterDistance(ctrlPoints[s + 1], genericCurve);

                    it = ctrlPoints.begin();
                    it1 = weights.begin();
                }

                // Evaluating if we have 2 or 3 near obstacles
                Eigen::Vector3d freeMediumPoint = {};
                int Median = -1;
                int cnt = 0;
                bool flagFreeMediumPoint = false;
                bool Median1 = false;

                if (!VecCenterIntersection.empty()) {
                    ptr_straightLine = nullptr;
                    double minParMedian = 999;
                    Median = -1;
                    for (const auto& ele : VecCenterIntersection) {
                        ptr_straightLine = make_shared<StraightLine>(obs_points[pos_first], obs_points[ele], 3, 2);
                        intersections.clear();
                        intersections = genericCurve->Intersection(ptr_straightLine);
                        if (!intersections.empty()) {
                            cnt++;
                            flagFreeMediumPoint = true;
                            temp1 = firstIntersection(intersections, genericCurve, false);
                            if (temp1.w() <= minParMedian) {
                                Median = ele;
                                minParMedian = temp1.w();
                            }
                        }
                    }

                    if (VecCenterIntersection.size() > 1) {
                        ptr_straightLine = make_shared<StraightLine>(obs_points[VecCenterIntersection.at(0)],
                            obs_points[VecCenterIntersection.at(1)], 3, 2);
                        intersections.clear();
                        intersections = genericCurve->Intersection(ptr_straightLine);
                        temp1 = firstIntersection(intersections, genericCurve, false);
                        intersections.emplace_back(temp1.x(), temp1.y(), temp1.z());
                        auto intersections1to2 = obs.at(VecCenterIntersection.at(0))->Intersection(obs_external.at(VecCenterIntersection.at(1)));
                        auto intersections2to1 = obs_external.at(VecCenterIntersection.at(1))->Intersection(obs.at(VecCenterIntersection.at(0)));

                        if (!intersections.empty() && (temp1.w() >= minParMedian) && (!intersections1to2.empty() || !intersections2to1.empty()))
                            Median1 = true;
                    }
                }

                // Evaluating what control point need to be added in the case of 2 or 3 near obstacles
                if (flagFreeMediumPoint) {

                    freeMediumPoint = FindCenterMedian(obs[pos_first], obs[Median]);

                    if (EvaluatePointInObstacle(freeMediumPoint) != 0)
                        return -10;

                    Eigen::Vector3d Last{};
                    if (EvaluatePointInObstacle(freeControlObs[0]) == 0)
                        Last = ctrlPoints[ctrlPoints.size() - 3];
                    else
                        Last = ctrlPoints[ctrlPoints.size() - 2];

                    vector<Eigen::Vector3d> TwoPointsExternal{};
                    TwoPointsExternal = ExternalIntersection(pos_first, VecCenterIntersection, Last, ctrlPoints[ctrlPoints.size() - 1], cnt, Median, Median1);

                    vector<Eigen::Vector3d> TempVector{};
                    TempVector.push_back(TwoPointsExternal[0]);

                    shared_ptr<StraightLine> LineTemp = nullptr;
                    LineTemp = make_shared<StraightLine>(TwoPointsExternal[0], ctrlPoints[ctrlPoints.size() - 2], 3, 2);

                    if (EvaluatePointInObstacle(TwoPointsExternal[0]) == 0) {
                        if (LineTemp->Length() < 15) {
                            ctrlPoints.at(index) = TwoPointsExternal[0];
                            weights.at(index) = 6;
                            UpdateNurbs(degree, knots, weights);
                            genericCurve = std::make_shared<GenericCurve>(degree, knots, ctrlPoints, weights, coeff,
                                dimension, degree + 1);
                            it = ctrlPoints.begin();
                            it1 = weights.begin();
                        }
                        else {
                            auto temp4 = firstIntersection(TempVector, genericCurve, true);
                            paramOrder.push_back(temp4.w());

                            index = UpdateParameterOrder(paramOrder);
                            ctrlPoints.insert(it + index, TwoPointsExternal[0]);
                            weights.insert(it1 + index, 10);
                            UpdateNurbs(degree, knots, weights);
                            genericCurve = std::make_shared<GenericCurve>(degree, knots, ctrlPoints, weights, coeff,
                                dimension, degree + 1);
                            it = ctrlPoints.begin();
                            it1 = weights.begin();
                        }
                    }

                    for (int s = 0; s < paramOrder.size(); s++)
                        paramOrder[s] = UpdateParameterDistance(ctrlPoints[s + 1], genericCurve);

                    TempVector.clear();
                    TempVector.push_back(freeMediumPoint);
                    auto temp4 = firstIntersection(TempVector, genericCurve, true);
                    paramOrder.push_back(temp4.w());

                    index = UpdateParameterOrder(paramOrder);
                    ctrlPoints.insert(it + index, freeMediumPoint);
                    weights.insert(it1 + index, 6);
                    UpdateNurbs(degree, knots, weights);
                    genericCurve = std::make_shared<GenericCurve>(degree, knots, ctrlPoints, weights, coeff, dimension, degree + 1);
                    it = ctrlPoints.begin();
                    it1 = weights.begin();

                    for (int s = 0; s < paramOrder.size(); s++)
                        paramOrder[s] = UpdateParameterDistance(ctrlPoints[s + 1], genericCurve);

                    if ((VecCenterIntersection.size() == 2) && (Median != -1)) {
                        auto itR = remove(VecCenterIntersection.begin(), VecCenterIntersection.end(), Median);
                        VecCenterIntersection.pop_back();
                    }

                    // Evaluating if we have 3 near obstacle
                    if (cnt == 2 || Median1) {

                        Eigen::Vector3d Baricenter;
                        Baricenter.x() = (obs_points[pos_first].x() + obs_points[Median].x() + obs_points[VecCenterIntersection[0]].x()) / 3;
                        Baricenter.y() = (obs_points[pos_first].y() + obs_points[Median].y() + obs_points[VecCenterIntersection[0]].y()) / 3;
                        Baricenter.z() = 0;

                        TempVector.clear();
                        TempVector.push_back(Baricenter);
                        temp4 = firstIntersection(TempVector, genericCurve, true);
                        paramOrder.push_back(temp4.w());

                        index = UpdateParameterOrder(paramOrder);
                        ctrlPoints.insert(it + index, Baricenter);
                        weights.insert(it1 + index, 5);
                        UpdateNurbs(degree, knots, weights);
                        genericCurve = std::make_shared<GenericCurve>(degree, knots, ctrlPoints, weights, coeff, dimension, degree + 1);
                        it = ctrlPoints.begin();
                        it1 = weights.begin();

                        for (int s = 0; s < paramOrder.size(); s++)
                            paramOrder[s] = UpdateParameterDistance(ctrlPoints[s + 1], genericCurve);

                        //Evaluation path between near obstacles
                        int from = -1;
                        if (Median1)
                            from = Median;
                        if (cnt == 2)
                            from = pos_first;

                        freeMediumPoint = FindCenterMedian(obs[from], obs[VecCenterIntersection[0]]);
                        TempVector.clear();
                        TempVector.push_back(freeMediumPoint);
                        temp4 = firstIntersection(TempVector, genericCurve, true);
                        paramOrder.push_back(temp4.w());

                        index = UpdateParameterOrder(paramOrder);
                        ctrlPoints.insert(it + index, freeMediumPoint);
                        weights.insert(it1 + index, 5);
                        UpdateNurbs(degree, knots, weights);
                        genericCurve = std::make_shared<GenericCurve>(degree, knots, ctrlPoints, weights, coeff, dimension, degree + 1);
                        it = ctrlPoints.begin();
                        it1 = weights.begin();

                        for (int s = 0; s < paramOrder.size(); s++)
                            paramOrder[s] = UpdateParameterDistance(ctrlPoints[s + 1], genericCurve);

                        TempVector.clear();

                        if (EvaluatePointInObstacle(TwoPointsExternal[1]) == 0) {

                            TempVector.push_back(TwoPointsExternal[1]);
                            temp4 = firstIntersection(TempVector, genericCurve, true);
                            paramOrder.push_back(temp4.w());

                            index = UpdateParameterOrder(paramOrder);
                            ctrlPoints.insert(it + index, TwoPointsExternal[1]);
                            weights.insert(it1 + index, 5);
                            UpdateNurbs(degree, knots, weights);
                            genericCurve = std::make_shared<GenericCurve>(degree, knots, ctrlPoints, weights, coeff,
                                dimension, degree + 1);
                            it = ctrlPoints.begin();
                            it1 = weights.begin();

                            for (int s = 0; s < paramOrder.size(); s++)
                                paramOrder[s] = UpdateParameterDistance(ctrlPoints[s + 1], genericCurve);
                        }

                    }
                    else {

                        // We have only two near obstacles
                        TempVector.clear();
                        TempVector.push_back(TwoPointsExternal[1]);
                        temp4 = firstIntersection(TempVector, genericCurve, true);
                        paramOrder.push_back(temp4.w());

                        index = UpdateParameterOrder(paramOrder);
                        ctrlPoints.insert(it + index, TwoPointsExternal[1]);
                        weights.insert(it1 + index, 4);
                        UpdateNurbs(degree, knots, weights);
                        genericCurve = std::make_shared<GenericCurve>(degree, knots, ctrlPoints, weights, coeff, dimension, degree + 1);
                        it = ctrlPoints.begin();
                        it1 = weights.begin();

                        for (int s = 0; s < paramOrder.size(); s++)
                            paramOrder[s] = UpdateParameterDistance(ctrlPoints[s + 1], genericCurve);
                    }
                }
            }

            tuple<double, double> closest = genericCurve->FindClosestPoint(ctrlPoints[ctrlPoints.size() - 2]);
            auto abscissa_m1 = (double)(get<0>(closest));
            ClosestPoint = genericCurve->At(double(abscissa_m1));
            wayPointsClosest.push_back(ClosestPoint);

            printNurbsInformation(knots);

        }
        else
            break;
    }
    pathRobotCalculation();
    CalculateWayPoints(wayPointsClosest);
}
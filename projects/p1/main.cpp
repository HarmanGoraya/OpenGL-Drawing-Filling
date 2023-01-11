#include "util.hpp"
#include "comm.hpp"

// force imgui to use GLAD
// #define IMGUI_IMPL_OPENGL_LOADER_GLAD
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl2.h>

#include <iostream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <cstddef>
#include <iostream>
#include <fstream>
#include <iomanip>
#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif



std::vector<Polygon> polygons;
int pid = 0; // current polygon
static bool show_gui = false;
std::string inputfile;

const int LeftBit = 0x1;
const int RightBit = 0x2;
const int BottomBit = 0x4;
const int TopBit = 0x8;

struct Edges{
    float miny;
    float maxy;
    float inverse_slope;
    float xval;
};

std::vector<Edges> AllEdges;
std::vector<Edges> ActiveEdges;


typedef enum { Left, Right, Bottom, Top } Boundary;

bool myfunction(const Edges &edges1, const Edges &edges2){
    return edges1.miny < edges2.miny;
}

bool sortbyX(const Edges &edges1, const Edges &edges2){
    return edges1.xval < edges2.xval;
}

void UpdateCenterofMass(){
    for (auto& p : polygons) {
        float cx = 0.f, cy = 0.f;
        for (auto& pt : p.points) {
            cx += pt.x;
            cy += pt.y;
        }
        cx /= (float)p.points.size();
        cy /= (float)p.points.size();
        p.cx = cx;
        p.cy = cy;
    }
}


static bool
CapturedByGUI()
{
    ImGuiIO& io = ImGui::GetIO();
    return (io.WantCaptureMouse);
}


void
KeyCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
        // take a screen shot
        ScreenShot(window, "p1.jpg");
        // close window
        glfwSetWindowShouldClose(window, GLFW_TRUE);
    }
    else if(key == GLFW_KEY_G && action == GLFW_PRESS){
        show_gui = !show_gui;

    };
  
}

void
CursorPositionCallback(GLFWwindow* window, double xpos, double ypos)
{
    if (!CapturedByGUI()) {
        int left_state  = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
        int right_state = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT);
        // left click
        if (left_state == GLFW_PRESS) {

        }
        else {
        }

        // right click
        if (right_state == GLFW_PRESS) {
        }
        else {
        }
    }
}

void lineDDA(void *window,float x0, float y0, float xEnd, float yEnd) {
    int dx = xEnd - x0;
    int dy = yEnd - y0;
    int steps, k;

    float xIncrement, yIncrement;
    float x = x0;
    float y = y0;

    if (abs(dx) > abs(dy)) {
        steps = abs(dx);
    } else {
        steps = abs(dy);
    }
    xIncrement = float(dx) / float(steps);
    yIncrement = float(dy) / float(steps);



    MakePix(window,round(x),round(y));
    for (k = 0; k < steps; k++) {
        x += xIncrement;
        y += yIncrement;
        MakePix(window,round(x),round(y));

    }
} 

void Bresenham(GLFWwindow* window,int x1, int y1, int x2, int y2){ 

    int x = x1;
    int y = y1;
    int interchange = 0;

    int dx = (x2 - x1);
    int dy = (y2 - y1);

    int S1,S2;

    if(dx < 0){
        S1 = -1;
    }
    else if(dx == 0){
        S1 = 0;
    }
    else{
        S1 = 1;
    }

    if(dy < 0){
        S2 = -1;
    }
    else if(dy == 0){
        S2 = 0;
    }
    else{
        S2 = 1;
    }

    //Special Cases

    if(dy == 0){ // Horizontal Line
        MakePix(window,x,y);
        for(int i = 0;i <abs(dx);i++){
            x += S1;
            MakePix(window,x,y);
        }
    }
    else if(dx == 0){ // Vertical Line
        MakePix(window,x,y);
        for(int i = 0;i<abs(dy);i++){
            y += S2;
            MakePix(window,x,y);
        }
    }

    else if(dx == dy){ // Same Slope Lines
        MakePix(window,x,y);
        for(int i = 0;i<abs(dx);i++){
            x += S1;
            y += S2;
            MakePix(window,x,y);
        }
    }

    // Main Brensenham Algorithim

    else{
        if(abs(dy) > abs(dx)){ // Make switch for slope > 1
            int temp = abs(dx);
            dx = abs(dy);
            dy = temp;
            interchange = 1;
        }

        int E = 2*abs(dy) - abs(dx);
        int A = 2 * abs(dy);
        int B = 2*abs(dy) - 2*abs(dx);

        MakePix(window,x,y); // Plot first point

        for(int i = 0;i<abs(dx);i++){
            if(E<0){
                if(interchange == 1){
                    y = y + S2;
                    E = E + A;
                }
                else {
                    x = x + S1;
                    E = E + A;
                }
            }
            else{
                y = y + S2;
                x = x + S1;
                E = E + B;
            }
            MakePix(window,x,y);
        }

    }


}

void writefile(const std::string filename){
    char buff[FILENAME_MAX]; //create string buffer to hold path
    GetCurrentDir( buff, FILENAME_MAX );
    std::string cwd(buff);
    std::string file = cwd + "/" + filename;
    std::ofstream myfile(inputfile);
   // myfile.open ("/Users/harman/Downloads/ecs175-demo.v2/projects/p1/sample.txt");
    myfile<<polygons.size();
    myfile<<"\n";
    for(int i =0;i<polygons.size();i++){
        myfile<<"\n";
        myfile<<polygons[i].points.size();
        myfile<<"\n";
        for(int j = 0;j<polygons[i].points.size();j++){
            myfile<<std::fixed<<std::setprecision(1)<<float(polygons[i].points.at(j).x);
            myfile<<" ";
            myfile<<std::fixed<<std::setprecision(1)<<float(polygons[i].points.at(j).y);
            myfile<<"\n";
        }
    }

    myfile.close();
}

int CodeInside(int code){ 
    return !code;
}
int reject(int code1, int code2){ 
    return (code1 & code2);
}
int accept(int code1, int code2){ 
    return !(code1 | code2);
}

unsigned char encode(Point point,float xmin,float xmax, float ymin, float ymax){ 
    
    unsigned char code = 0x00;

    if (point.x < xmin) code = code | LeftBit;
    if (point.x > xmax) code = code | RightBit;
    if (point.y < ymin) code = code | BottomBit;
    if (point.y > ymax) code = code | TopBit;
    return (code);
}

void CohenSuther(Point p1, Point p2, float xmin, float xmax, float ymin, float ymax){ 

    unsigned char code1, code2;
    bool done, drawline = false;
    float m;

    while(!done){
        code1 = encode(p1,xmin, xmax, ymin, ymax);
        code2 = encode(p1,xmin, xmax, ymin, ymax);
        if(accept(code1,code2)){ // Both points are wihtin the clip polygon
            done = true;
            drawline = true;
        }
        else{
            if(reject(code1,code2)){ // both points are outside the clip polygon
                done = true;
            }
            else{
                if(CodeInside(code1)){ // If first point is inside the clip box
                    Point temp;  // Swap Points; Make P1 be the point outside the box
                    temp = p1;
                    p1 = p2;
                    p2 = temp;

                    unsigned char tmp; // Also have to swap the codes as well because we swapped the points
                    tmp = code1;
                    code1 = code2;
                    code2 = tmp;
                }
                if (p2.x != p1.x) {
                    m = (p2.y - p1.y) / (p2.x - p1.x); // Find slope in order to find intersection point
                }
                if (code1 & LeftBit) {
                    p1.y += (xmin - p1.x) * m;
                    p1.x = xmin;
                }
                else if (code1 & RightBit) {
                    p1.y += (xmax - p1.x) * m;
                    p1.x = xmax;
                }
                else if (code1 & BottomBit) { /* Need to update p1.x for nonvertical lines only. */
                    if (p2.x != p1.x) {
                        p1.x += (ymin - p1.y) / m;
                        p1.y = ymin;
                    }
                }
                else{
                    if (code1 & TopBit) {
                        if (p2.x != p1.x) {
                            p1.x += (ymax - p1.y) / m;
                            p1.y = ymax;
                        }
                    }
                }


            }
        }



    }
}

bool inside(float x, float y,Boundary b, float xmin, float ymin, float xmax, float ymax){ 

    switch(b){

        case 0: if (x < xmin) return false; break;
        case 1: if (x > xmax) return false; break;
        case 2: if (y < ymin) return false; break;
        case 3: if (y > ymax) return false; break;
    }
    return true;

}

Point intersect(Point p1, Point p2, Boundary winEdge, float xmin, float xmax, float ymin, float ymax){ 

    Point iPt;
    float m;

    if(p1.x != p2.x){
        m = (p1.y - p2.y) / (p1.x - p2.x);
    }

    switch (winEdge) {
        case Left: iPt.x = xmin; iPt.y = p2.y + (xmin - p2.x) * m; break;
        case Right: iPt.x = xmax; iPt.y = p2.y + (xmax - p2.x) * m; break;
        case Bottom: iPt.y = ymin; if (p1.x != p2.x) iPt.x = p2.x + (ymin - p2.y) / m; else iPt.x = p2.x; break;
        case Top: iPt.y = ymax; if (p1.x != p2.x) iPt.x = p2.x + (ymax - p2.y) / m; else iPt.x = p2.x; break;
    }

    return iPt;
}

std::vector<Point> SutherLand(const std::vector<Point>&Intial,std::vector<Point>ClipPoints,float xmin,float xmax,float ymin, float ymax){

    std::vector<Point> Output; // End Result Vector
    std::vector<Point>Temp; // Used to store Clipped Vertices from each iteration

    Point PolyEdge1;
    Point PolyEdge2;

    for(int i = 0;i<4;i++){ // 4 becuase clipping window is a square so four edges
        Boundary b;

        if(i != 0){
            Output = Temp;
            Temp.clear();
        }
        else{
            Output = Intial;
        }


        switch(i){ // Determine ClipEdge
            case 0: b = Left;break;
            case 1: b = Right;break;
            case 2: b = Bottom;break;
            case 3: b = Top;break;
            default:break;
        }

        for(int j = 0;j<Output.size();j++){

            if(j != Output.size() - 1){ // Not Last Edge
                PolyEdge1 = Output[j];
                PolyEdge2 = Output[j+1];

            }
            else{ // Last Edge
                PolyEdge1 = Output[j];
                PolyEdge2 = Output[0];
            }


            if(inside(PolyEdge1.x,PolyEdge1.y,b,xmin,ymin,xmax,ymax) && // Both points are inside push second point
                inside(PolyEdge2.x,PolyEdge2.y,b,xmin,ymin,xmax,ymax)){

                    Temp.push_back(PolyEdge2);
                }


            else if(inside(PolyEdge1.x,PolyEdge1.y,b,xmin,ymin,xmax,ymax) and // First point is inside and second is outside
                    !inside(PolyEdge2.x,PolyEdge2.y,b,xmin,ymin,xmax,ymax)){
                Temp.push_back(intersect(PolyEdge1,PolyEdge2,b,xmin,xmax,ymin,ymax));

            }

            else if(!inside(PolyEdge1.x,PolyEdge1.y,b,xmin,ymin,xmax,ymax) and // Both points are outside
                    !inside(PolyEdge2.x,PolyEdge2.y,b,xmin,ymin,xmax,ymax)){

            }

            else if(!inside(PolyEdge1.x,PolyEdge1.y,b,xmin,ymin,xmax,ymax) and // Second point is inside first is outside
                    inside(PolyEdge2.x,PolyEdge2.y,b,xmin,ymin,xmax,ymax)){

                Temp.push_back(intersect(PolyEdge1,PolyEdge2,b,xmin,xmax,ymin,ymax));
                Temp.push_back(PolyEdge2);

            }


        }
    }

    Output = Temp;
    return Output;
}

void RemoveEdges(int scanline){
    std::vector<int>indices;
    for(int i = 0;i<ActiveEdges.size();i++){
        if(int(ceil(ActiveEdges[i].maxy)) == scanline){
            indices.push_back(i);

        }
    }

    for(int j = 0;j<indices.size();j++){
        if(j == 0){
            ActiveEdges.erase(ActiveEdges.begin() + indices[j]);
        }
        else{
            ActiveEdges.erase(ActiveEdges.begin() + (indices[j] - j));
        }
    }
}

int Scan(void* window, int ymax){
    int scanline;

    if(AllEdges.empty()){
        return -1;
    }


    ActiveEdges.push_back(AllEdges[0]);


    AllEdges.erase(AllEdges.begin());

    scanline = int(ceil(ActiveEdges[0].miny));

    for(scanline;scanline < ymax;scanline++){
        for(auto edge:AllEdges){
            if(int(ceil(edge.miny)) == scanline && isinf(edge.inverse_slope) == 0){
                ActiveEdges.push_back(edge);
            }
        }
        std::sort(ActiveEdges.begin(),ActiveEdges.end(),sortbyX);

        RemoveEdges(scanline);
        for(int i = 0;i<ActiveEdges.size();i+=2){
            for(int j = int(ceil(ActiveEdges[i].xval));j<int(ceil(ActiveEdges[i+1].xval));j++){
                MakePix(window,j,scanline);
            }
            ActiveEdges[i].xval += ActiveEdges[i].inverse_slope;
            ActiveEdges[i+1].xval += ActiveEdges[i+1].inverse_slope;
        }

    }

    AllEdges.clear();
    ActiveEdges.clear();

    return 0;
}


void
DrawCall(GLFWwindow* window, int algo)
{
    float ratio;
    int   width, height;
    glfwGetFramebufferSize(window, &width, &height);
    ratio = width / (float)height;
    //glViewport(250, 0, width, height);
    //glClear(GL_COLOR_BUFFER_BIT);
    //glMatrixMode(GL_PROJECTION);
    //glLoadIdentity();
    //glOrtho(-ratio, ratio, -1.f, 1.f, 1.f, -1.f);

    for(int j = 0;j<polygons.size();j++){
        for(int i = 0;i<polygons[j].clipped.size();i++){
            if(i == polygons[j].clipped.size() - 1){
                Edges temp;
                temp.miny = std::min(polygons[j].clipped.at(i).y,polygons[j].clipped.at(0).y);
                temp.maxy = std::max(polygons[j].clipped.at(i).y,polygons[j].clipped.at(0).y);
                polygons[j].clipped.at(i).y > polygons[j].clipped.at(0).y ? temp.xval = polygons[j].clipped.at(0).x :
                        temp.xval = polygons[j].clipped.at(i).x;
                temp.inverse_slope = 1 / ((polygons[j].clipped.at(i).y - polygons[j].clipped.at(0).y) /
                        (polygons[j].clipped.at(i).x - polygons[j].clipped.at(0).x));

                if(isinf(temp.inverse_slope) == 0){
                    AllEdges.push_back(temp);
                }

                if(algo == 0){
                    lineDDA(window,polygons[j].clipped.at(i).x,polygons[j].clipped.at(i).y,
                            polygons[j].clipped.at(0).x,polygons[j].clipped.at(0).y);
                }
                else{
                    Bresenham(window,polygons[j].clipped.at(i).x,polygons[j].clipped.at(i).y,
                            polygons[j].clipped.at(0).x,polygons[j].clipped.at(0).y);
                }

            }
            else{

                Edges temp;
                temp.miny = std::min(polygons[j].clipped.at(i).y,polygons[j].clipped.at(i + 1).y);
                temp.maxy = std::max(polygons[j].clipped.at(i).y,polygons[j].clipped.at(i + 1).y);
                polygons[j].clipped.at(i).y > polygons[j].clipped.at(i+1).y ? temp.xval = polygons[j].clipped.at(i+1).x :
                        temp.xval = polygons[j].clipped.at(i).x;
                temp.inverse_slope = 1 / ((polygons[j].clipped.at(i+1).y - polygons[j].clipped.at(i).y) /
                                                                       (polygons[j].clipped.at(i+1).x - polygons[j].clipped.at(i).x));
                if(isinf(temp.inverse_slope) == 0){
                    AllEdges.push_back(temp);
                }

                if(algo == 0){
                    lineDDA(window,polygons[j].clipped.at(i).x,polygons[j].clipped.at(i).y,
                            polygons[j].clipped.at(i+1).x,polygons[j].clipped.at(i+1).y);
                }
                else{
                    Bresenham(window,polygons[j].clipped.at(i).x,polygons[j].clipped.at(i).y,
                            polygons[j].clipped.at(i+1).x,polygons[j].clipped.at(i+1).y);
                }

            }

        }
    }

    std::sort(AllEdges.begin(), AllEdges.end(),myfunction);

    
}

void transform(float tx, float ty){

    for(auto & point : polygons[pid].points){
        point.x = point.x + tx;
        point.y = point.y + ty;

    }

    UpdateCenterofMass();
} 

void scale(float sx, float sy){ // Book 7-1

    for(int i = 0;i<polygons[pid].points.size();i++){
        polygons[pid].points.at(i).x = polygons[pid].points.at(i).x * sx + polygons[pid].cx * (1 - sx);
        polygons[pid].points.at(i).y = polygons[pid].points.at(i).y * sy + polygons[pid].cy * (1 - sy);

    }

    UpdateCenterofMass();
}

void rotate(double r){

    r = r * M_PI / 180; // Convert to Radians

    for(int i = 0;i<polygons[pid].points.size();i++){
        float temp = polygons[pid].cx + (polygons[pid].points.at(i).x - polygons[pid].cx)*cos(r) // Do not want to update x as it is used in y calculation
                                        - (polygons[pid].points.at(i).y - polygons[pid].cy)*sin(r);

        polygons[pid].points.at(i).y = polygons[pid].cy + (polygons[pid].points.at(i).x - polygons[pid].cx)*sin(r)
                                       + (polygons[pid].points.at(i).y - polygons[pid].cy)*cos(r);

        polygons[pid].points.at(i).x = temp;

    }

    UpdateCenterofMass(); // Make sure to update center of mass for new points


} 

void ImGui_Setup(GLFWwindow* window){

    ImGui::Begin("Options");

    ImGui::Text("Line Drawing Algorithim");
    static int Algorithim = -1;
    const char* items[] = {"DDA", "Bresenheim"};
    ImGui::ListBox("",&Algorithim,items,2,2);
    static bool fill = false;
    ImGui::SameLine();
    ImGui::Checkbox("Fill",&fill);


    static float xmin = 0.0;
    static float ymin = 0.0;
    static float xmax = 400.0;
    static float ymax = 400.0;

    static float inputs[] = {0.0,400.0,0.0,400.0};



    ImGui::Text("xmin");
    ImGui::SameLine(65);
    ImGui::Text("xmax");
    ImGui::SameLine(130);
    ImGui::Text("ymin");
    ImGui::SameLine(195);
    ImGui::Text("ymax");

    ImGui::InputFloat4("vals",inputs);

    xmin = inputs[0];
    ymin = inputs[2];
    xmax = inputs[1];
    ymax = inputs[3];

    if(xmin > xmax){
        xmax = 400.0;
    }
    if(ymin > ymax){
        ymax = 400.0;
    }

    std::vector<Point>ClipPoints;
    ClipPoints.push_back({xmin,ymax}); // Left Clipper
    ClipPoints.push_back({xmin,ymin});

    ClipPoints.push_back({xmax,ymin}); // Right Clipper
    ClipPoints.push_back({xmax,ymax});

    ClipPoints.push_back({xmin,ymin}); // Bottom Clipper
    ClipPoints.push_back({xmax,ymin});

    ClipPoints.push_back({xmax,ymax}); // Top Clipper
    ClipPoints.push_back({xmax,ymax});

    if(Algorithim != -1){

        for(int i = 0;i<polygons.size();i++){
            polygons[i].clipped = SutherLand(polygons[i].points,ClipPoints,xmin,xmax,ymin,ymax);
        }

        DrawCall(window,Algorithim);
        if(fill){
            Scan(window,int(ceil(ymax)));
        }

    }

    ImGui::Text("Select Polygon to manipulate");
    ImGui::SliderInt("PolygonID",&pid,0,polygons.size()-1);

    ImGui::Text("Input Translation Vector");

    static float tranx = 0.0;
    static float trany = 0.0;

    ImGui::InputFloat("X",&tranx);
    ImGui::InputFloat("Y",&trany);

    static bool check = false;
    ImGui::Checkbox("Apply Translation",&check);

    if(check){
        transform(tranx,trany);
        check = false;
    }

    ImGui::Text("Scale Polygon");
    static float scalept = 0.0;

    ImGui::InputFloat("Scale",&scalept);
    static bool check1 = false;

    ImGui::Checkbox("Scale Polygon",&check1);

    if(check1){
        scale(scalept,scalept);
        check1 = false;
    }

    static float rotation = 0.0;
    static bool check5 = false;
    ImGui::Checkbox("Rotate",&check5);
    ImGui::SliderFloat("Rotation Angle", &rotation,-2.0,2.0);
    if(check5){
        rotate(rotation);
    }

    char buffer[30];

    std::string buf;
    static bool check4 = false;
    static int counter = 0;
    ImGui::Checkbox("Make File",&check4);

    if(check4){
        counter++;
        buf = "Output" + std::to_string(counter) + ".txt";
        writefile(buf);
        check4 = false;
    }

    ImGui::End();


}


int
main(const int argc, const char** argv)
{


  if (argc < 2)
    throw std::runtime_error("missing input file");

  inputfile = argv[1];

  
  ReadFile(argv[1]);

  // Compute the Center of Mass
  for (auto& p : polygons) {        
    float cx = 0.f, cy = 0.f;
    for (auto& pt : p.points) {
      cx += pt.x;
      cy += pt.y;
    }
    cx /= (float)p.points.size();
    cy /= (float)p.points.size();
    p.cx = cx;
    p.cy = cy;
  }

    for(int i = 0;i<polygons.size();i++){
        for(int j = 0;j<polygons[i].points.size();j++){
            polygons[i].points.at(j).x = polygons[i].points.at(j).x * 30;
            polygons[i].points.at(j).y = polygons[i].points.at(j).y * 30;
        }
    }

    UpdateCenterofMass();


    int width = 400, height = 400;

    // Initialize GLFW
    glfwSetErrorCallback(ErrorCallback);
    if (!glfwInit()) {
        exit(EXIT_FAILURE);
    }

    // Provide Window Hint
    glfwWindowHint(GLFW_RESIZABLE, GL_TRUE);

    // OpenGL Setup
    GLFWwindow* window = NULL;
    window = glfwCreateWindow(width, height, "ECS 175 Renderer", NULL, NULL);
    if (!window) {
        glfwTerminate();
        throw std::runtime_error("Failed to create GLFW window");
    }

    // Ready
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    // Load GLAD symbols
    int err = gladLoadGLLoader((GLADloadproc)glfwGetProcAddress) == 0;
    if (err) {
        throw std::runtime_error("Failed to initialize OpenGL loader!");
    }

    // Callback
    glfwSetKeyCallback(window, KeyCallback);
    glfwSetCursorPosCallback(window, CursorPositionCallback);

    // Execute
    const char* glsl_version = "#version 110";

    {
        // Setup Dear ImGui context
        IMGUI_CHECKVERSION();
        ImGui::CreateContext();

        // Setup Dear ImGui style
        ImGui::StyleColorsDark(); // or ImGui::StyleColorsClassic();

        // Initialize Dear ImGui
        ImGui_ImplGlfw_InitForOpenGL(window, true);
        ImGui_ImplOpenGL2_Init();
    }

    while (!glfwWindowShouldClose(window)) {

        glClear(GL_COLOR_BUFFER_BIT);
        glFlush();

        ImGui_ImplOpenGL2_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        ImGui_Setup(window);

        ImGui::Render();
        ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());


        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    ImGui_ImplOpenGL2_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();


    // Exit
    glfwDestroyWindow(window);
    glfwTerminate();

  return EXIT_SUCCESS;
}

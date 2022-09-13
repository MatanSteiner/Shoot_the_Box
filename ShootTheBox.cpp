#include <iostream>
#include <vector>
#include <array>
#include <chrono>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <Windows.h>

#define PI 3.1415926//53589793238


const short screenwidth = 1920*0.7, screenhight = 1080*0.7, screenmargin=2;



void GLAPIENTRY MessageCallback(GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, const GLchar* message, const void* userParam)
{

    std::cout << std::endl << "GL Message:" << std::endl << "Type: " << type << ", ID: " << id << ", Severity: " << severity << std::endl;
    std::cout << "Message:" << std::endl << message << std::endl;
}



class VertexBuffer
{
    unsigned int m_RendererID, m_size, m_type;

public:
   
    VertexBuffer(const void* data, unsigned int size, unsigned int type)
        :m_size(size), m_type(type)
    {
        glGenBuffers(1, &m_RendererID);
        glBindBuffer(GL_ARRAY_BUFFER, m_RendererID);
        glBufferData(GL_ARRAY_BUFFER, size, data, type);
    }

    VertexBuffer(unsigned int size, unsigned int type)
        :m_size(size), m_type(type)
    {
        glGenBuffers(1, &m_RendererID);
        glBindBuffer(GL_ARRAY_BUFFER, m_RendererID);
    }

    VertexBuffer()
        :m_size(sizeof(float)), m_type(GL_STATIC_DRAW)
    {
        glGenBuffers(1, &m_RendererID);
        glBindBuffer(GL_ARRAY_BUFFER, m_RendererID);
    }

    ~VertexBuffer()
    {
        glDeleteBuffers(1, &m_RendererID);
    }


    void Bind()
    {
        glBindBuffer(GL_ARRAY_BUFFER, m_RendererID);

    }
    void Unbind()
    {
        glBindBuffer(GL_ARRAY_BUFFER,0);

    }

    void newData(const void* data)
    {
        Bind();
        glBufferData(GL_ARRAY_BUFFER, m_size, data, m_type);
        
    }

    void newData (const void* data, unsigned int size, unsigned int type)
    {
        Bind();
        m_size = size;
        m_type = type;
        glBufferData(GL_ARRAY_BUFFER, m_size, data, m_type);
       
    }

};

class IndexBuffer
{

    unsigned int m_RendererID;
    unsigned int m_count, m_type;

public:

    IndexBuffer(const unsigned int* data, unsigned int count, unsigned int type)
        :m_count(count), m_type(type)
    {
        glGenBuffers(1, &m_RendererID);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_RendererID);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, count * sizeof(unsigned int), data, type);
    }

    IndexBuffer(unsigned int count, unsigned int type)
        :m_count(count), m_type(type)
    {
        glGenBuffers(1, &m_RendererID);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_RendererID);
    }

    IndexBuffer()
        :m_count(1), m_type(GL_STATIC_DRAW)
    {
        glGenBuffers(1, &m_RendererID);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_RendererID);
    }


    ~IndexBuffer()
    {
        glDeleteBuffers(1, &m_RendererID);
    }

    void Bind()
    {
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_RendererID);

    }

    void Unbind()
    {
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    }

    void newData(const unsigned int* data)
    {
        Bind();
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, m_count * sizeof(unsigned int), data, m_type);
    }

    void newData(const unsigned int* data, unsigned int count, unsigned int type)
    {
        m_count = count;
        m_type = type;
        Bind();
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, count * sizeof(unsigned int), data, type);
    }

};


struct ShaderProgramSource
{
    std::string VertexSource;
    std::string FragmentSource;
};



 class Shader
 {
     int colourlocation;
 public:

     unsigned int shader;

     Shader(const std::array<float,4>& colour)
     {
         ShaderProgramSource source = parseShader("Shaders/");

         shader = createShader(source.VertexSource, source.FragmentSource);
         glUseProgram(shader);

         colourlocation = glGetUniformLocation(shader, "u_color");
         glUniform4f(colourlocation,colour[0], colour[1], colour[2], colour[3]);
     }

     Shader()
     {
         ShaderProgramSource source = parseShader("Shaders/");

         shader = createShader(source.VertexSource, source.FragmentSource);
         glUseProgram(shader);

         colourlocation = glGetUniformLocation(shader, "u_color");
         glUniform4f(colourlocation, 1.0,1.0, 1.0, 1.0);
     }
     ~Shader()
     {
         glDeleteProgram(shader);

     }

     ShaderProgramSource parseShader(const std::string& filepath)
     {
         std::ifstream stream(filepath + "Vertex.shader");
         std::string line;
         std::stringstream ss[2];
         for (int i = 0; i < 2; i++) {
             while (getline(stream, line))
                 ss[i] << line << '\n';
             stream.close();
             stream.open(filepath + "Fragment.shader");
         }
         stream.close();

         return { ss[0].str(), ss[1].str() };
     }

     unsigned int compileShader(unsigned int type, const std::string& source) {

         unsigned int id = glCreateShader(type);
         const char* src = source.c_str();
         glShaderSource(id, 1, &src, nullptr);
         glCompileShader(id);

         // troubleshooting:
         int result;
         glGetShaderiv(id, GL_COMPILE_STATUS, &result);
         if (!result) {

             int length;
             glGetShaderiv(id, GL_INFO_LOG_LENGTH, &length);

             //char message[length] <- This doesn't work, so I'm allocating memory on the stack as a char pointer:
             char* message = (char*)alloca(length * sizeof(char));
             glGetShaderInfoLog(id, length, &length, message);

             std::cout << "Failed to compile " << (type == GL_VERTEX_SHADER ? "vertex" : "fragment") << " shader" << std::endl << message << std::endl;
             glDeleteShader(id);
             return 0;
         }


         return id;
     }


     unsigned int createShader(const std::string& vertexShader, const std::string& fragmentShader) {

         unsigned int program = glCreateProgram();
         unsigned int vs = compileShader(GL_VERTEX_SHADER, vertexShader);
         unsigned int fs = compileShader(GL_FRAGMENT_SHADER, fragmentShader);

         glAttachShader(program, vs);
         glAttachShader(program, fs);
         glLinkProgram(program);
         glValidateProgram(program);

         glDeleteShader(vs);
         glDeleteShader(fs);

         return program;

     }

     void colour(const std::array<float, 4>& colour)
     {
         glUniform4f(colourlocation, colour[0], colour[1], colour[2], colour[3]);
     }


 };
 


 struct coord {
     float x, y;
     
     coord (float X, float Y)
         : x(X), y(Y)
     {}

     coord ()
         : x(0), y(0)
     {}

     coord screenNorm() const { 
         return { x / (screenwidth / 2), y / (screenhight / 2) }; 
     }

     float vecSize() const {
         return sqrt(x * x + y * y); 
     }

     coord unitVec() const {
         return { x / vecSize(), y / vecSize()};
     }

     coord operator+ (const coord &other)  const {
         return { x + other.x, y + other.y };
     }

     coord operator- (const coord &other) const {
         return { x - other.x, y - other.y };
     }

     float operator* (const coord &other) const {
         return { x * other.x + y * other.y };
     }

     coord operator* (const float &cnst) const {
         return { x * cnst , y * cnst };
     }

 };

 

 class Box
 {

     float m_alpha, m_r;
     unsigned int m_drawmode = GL_TRIANGLES, m_dynDraw = GL_DYNAMIC_DRAW;

 public:


     unsigned int indices[6] = {
       2,3,0,
       0,1,2
     };
     coord r0, c;
     unsigned int RendererID;
     float omega, m, I, theta0 = 0.0;
     coord v;
     float normPositions[8];
     std::array<coord, 4> positions;
     VertexBuffer vb;
     IndexBuffer ib;
     std::array<float, 4> colour = { 1.0, 1.0, 1.0, 1.0 };
     Shader* shader;


     Box(coord r0, float h, float b, float m, float omega0, coord v0, Shader* shader)
         :r0(r0), omega(omega0), v(v0), m(m), shader(shader)
     {
         m_r = sqrt(h * h + b * b) / 2;
         m_alpha = 2 * atan(h / b);
         I = m * m_r * m_r / 3;
         setPositions(0.0);
         initGL();
     }

     Box(coord r0, float h, float b, float m, Shader* shader)
         :r0(r0), omega(0.0), v({ 0,0 }), m(m), shader(shader)
     {

         if (m <= 0)
             m_dynDraw = GL_STATIC_DRAW;

         if (m < 0)
         {
             m_drawmode = GL_LINE_LOOP;
             m = 0;
         }

         m_r = sqrt(h * h + b * b) / 2;
         m_alpha = 2 * atan(h / b);
         I = m * m_r * m_r / 3;
         setPositions(0.0);
         initGL();
     }



     void initGL()
     {

         vb.newData(normPositions, sizeof(normPositions), m_dynDraw);

         ib.newData(indices, 6, m_dynDraw);
     }

     ~Box()
     {

     }

     void  setPositions(const float& t)
     {

         float theta = -m_alpha / 2 + theta0;

         c = r0 + v * t;
         for (int i = 0; i < 8; i += 2)
         {
             positions[i / 2].x = m_r * cos(omega * t + theta) + c.x;
             positions[i / 2].y = m_r * sin(omega * t + theta) + c.y;
             normPositions[i] = positions[i / 2].x / (screenwidth / 2);
             normPositions[i + 1] = positions[i / 2].y / (screenhight / 2);
             (i % 4) ? theta += (PI - m_alpha) : theta += m_alpha;
         }

     }

     void resize(coord r0new, float h, float b)
     {
         if (m != 0)
             return;
         m_r = sqrt(h * h + b * b) / 2;
         m_alpha = 2 * atan(h / b);
         r0 = r0new;
         setPositions(0.0);
     }

     void timeReset(float t)
     {
         r0 = r0 + v * t;
         theta0 = theta0 + omega * t;

     }

     void draw()
     {
         shader->colour(colour);
         vb.newData(normPositions);
         glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(float) * 2, 0);
         ib.Bind();
         glDrawElements(m_drawmode, 6, GL_UNSIGNED_INT, nullptr);

     }

     void unbind()
     {
         vb.Unbind();
         ib.Unbind();

     }

 };


 class Circle
 {


     int m_triMesh;

 public:

     float r;
     float m, I;
     coord r0;
     unsigned int* indices_p;
     std::vector<unsigned int> indices;
     unsigned int RendererID;
     coord v;
     coord c;
     float* normPositions_p;
     std::vector<float> positions, normPositions;
     VertexBuffer vb;
     IndexBuffer ib;
     std::array<float, 4> colour = {1.0, 1.0, 1.0, 1.0};
     Shader* shader;


     Circle(coord r0, float r, int triMesh, float m, coord v0, Shader* shader)
         :r0(r0), r(r), v(v0), m_triMesh(triMesh), m(m), shader(shader)
     {
         initIndicis();
         setPositions(0.0);
         initGL();
     }

     Circle(coord r0, float r, int triMesh, float m, Shader* shader)
         :r0(r0), r(r), m_triMesh(triMesh), m(m), shader(shader)
     {
         v.x = 0.0;
         v.y = 0.0;
         initIndicis();
         setPositions(0.0);
         initGL();
     }

     void initIndicis()
     {

         I = m * r * r / 2;

         indices.resize(m_triMesh * 3);

         for (int i = 0; i < m_triMesh; i++)
         {
             indices[3 * i] = 0;
             indices[3 * i + 1] = i + 1;
             indices[3 * i + 2] = i + 2;
         }
         indices[indices.size() - 1] = 1;
         indices_p = &indices[0];

     }

     void initGL()
     {
         vb.newData(normPositions_p, sizeof(float) * (m_triMesh * 2 + 2), GL_DYNAMIC_DRAW);

         ib.newData(indices_p, m_triMesh * 3, GL_DYNAMIC_DRAW);
     }

     ~Circle()
     {

     }

     void  setPositions(const float& t)
     {
         positions.resize(m_triMesh * 2 + 2);
         normPositions.resize(m_triMesh * 2 + 2);

         positions[0] = r0.x + v.x * t;
         positions[1] = r0.y + v.y * t;
         normPositions[0] = positions[0] / (screenwidth / 2);
         normPositions[1] = positions[1] / (screenhight / 2);

         const float deltatheta = 2 * PI / m_triMesh;
         float theta = 0.0;

         for (int i = 2; i < m_triMesh * 2 + 2; i += 2)
         {
             positions[i] = r * cos(theta) + r0.x + v.x * t;
             positions[i + 1] = r * sin(theta) + r0.y + v.y * t;
             normPositions[i] = positions[i] / (screenwidth / 2);
             normPositions[i + 1] = positions[i + 1] / (screenhight / 2);
             theta += deltatheta;
         }

         normPositions_p = &normPositions[0];
         c.x = positions[0];
         c.y = positions[1];

     }

     void draw()
     {
         shader->colour(colour);
         vb.newData(normPositions_p);
         glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(float) * 2, 0);
         ib.Bind();
         glDrawElements(GL_TRIANGLES, m_triMesh * 3, GL_UNSIGNED_INT, nullptr);

     }

     void unbind()
     {
         vb.Unbind();
         ib.Unbind();

     }

 };


 class Line
 {

 public:

     coord positions[2];
     float normPositions[4];
     unsigned int indices[2] = { 0,1 };
     int length;
     float beta = 1E09;
     coord r0, r1;
     unsigned int RendererID;
     VertexBuffer vb;
     IndexBuffer ib;
     std::array<float, 4> colour = { 1.0, 1.0, 1.0, 1.0 };
     Shader* shader;


     Line(coord r0, int length, float beta, Shader* shader)
         :r0(r0), length(length), beta(beta), shader(shader)
     {
         initGL();
         setPositions();
     }

     Line(coord r0, coord r1, Shader* shader)
         :r0(r0), r1(r1), shader(shader)
     {
         initGL();
         setPositions();
     }

     void initGL()
     {
         vb.newData(normPositions, sizeof(normPositions), GL_DYNAMIC_DRAW);

         ib.newData(indices, 2, GL_DYNAMIC_DRAW);

     }

     void  setPositions()
     {

         positions[0] = r0;
         if (beta == 1E09)
            positions[1] = r1;
         else
             positions[1] = { r0.x + length * cos(beta), r0.y + length * sin(beta)};

         for (int i = 0; i <= 2; i += 2)
         {
             normPositions[i] = positions[i / 2].x / (screenwidth / 2);
             normPositions[i + 1] = positions[i / 2].y / (screenhight / 2);
         }

     }

     void draw()
     {
         shader->colour(colour);
         vb.newData(normPositions);
         glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(float) * 2, 0);
         ib.Bind();
         glDrawElements(GL_LINES, 2, GL_UNSIGNED_INT, nullptr);

     }

     void unbind()
     {
         vb.Unbind();
         ib.Unbind();

     }

 };



 bool contact(const std::array<coord, 4>& bv/*box vertices*/, const coord& c/*circle centre*/, const float& r, coord* contactSurface, coord* contactNorm)
 {
     for (char i = 0; i < 4; i++)
     {

         coord v = bv[i] - c;
         
         int j;

         (i == 0) ? (j = 3) : (j = i - 1);

         coord u = bv[i] - bv[j];
         coord u_unit = u.unitVec();

         coord impact = bv[j] + (u_unit * (u.vecSize() - v * u_unit));

         // distance between bv[i] and c < r and the coordinates are within the box's vertices
         
         if ((v * v) - ((u_unit * v) * (u_unit * v)) <= r * r &&
             impact.x >= min(bv[j].x, bv[i].x) && impact.x <= max(bv[j].x, bv[i].x) &&
             impact.y >= min(bv[j].y, bv[i].y) && impact.y <= max(bv[j].y, bv[i].y))
         {
             if (contactNorm) 
             {
                 *contactSurface = u_unit;
                 int k;
                 (j == 0) ? (k = 3) : (k = j - 1);
                 *contactNorm = (bv[j] - bv[k]).unitVec();
             }
             return 1;
         }
     }
     return 0;

 }
 short hitormiss(const std::array<coord, 4>& bv, const coord* gate)
 {
         char n = 0;
     for (char i = 0; i < 4; i++)
     {
         if (bv[i].x <= -screenwidth / 2 + screenmargin ||
             bv[i].y <= -screenhight / 2 + screenmargin || bv[i].y >= screenhight / 2 - screenmargin)
             return -1;

         if (bv[i].x >= screenwidth / 2 - screenmargin)
             if (bv[i].y > gate[1].y || bv[i].y < gate[0].y && bv[i].x <= screenwidth / 2 + screenmargin)
                 return -1;
             else
                 n++;
     }
     if (contact(bv, gate[0], 1, nullptr, nullptr))
         return -1;
     if (contact(bv, gate[1], 1, nullptr, nullptr))
         return -1;

     if (n == 4)
         return 1;
     return 0;

 }

 short bulletout(coord c)
 {
     if (c.x > screenwidth / 2 || c.y > screenhight / 2 || c.y < -screenhight / 2)
         return -1;
     return 0;
 }
 
 void momentumTransferNoFric(Box& box, Circle& circle, const coord& contactSurface, const coord& contactNorm)
 {
     float v_box_norm=box.v*contactNorm, v_box_II=box.v*contactSurface, v_circle_norm=circle.v*contactNorm, v_circle_II = circle.v * contactSurface;
     float d = (box.c - circle.c) * contactSurface;
     
     float A = circle.m * d * v_circle_norm / box.I + box.omega;
     float B = circle.m * v_circle_norm / box.m + v_box_norm;

     float a = circle.m + circle.m * circle.m * d * d / box.I + circle.m * circle.m / box.m;
     float b = -2 * circle.m * (A * d + B);
     float c = box.I * A * A + box.m * B * B - box.I * box.omega * box.omega - circle.m * v_circle_norm * v_circle_norm - box.m * v_box_norm * v_box_norm;

     float v_circle_new = ( - b + sqrt(b * b - 4 * a * c)) / 2 / a;
     float v_box_new = (circle.m / box.m) * (v_circle_norm - v_circle_new) + v_box_norm;
     box.omega = box.omega + (circle.m * d / box.I)*(v_circle_norm - v_circle_new);
     circle.v = (contactNorm * v_circle_new) + (contactSurface * v_circle_II);
     box.v= (contactNorm * v_box_new) + (contactSurface * v_box_II);
     circle.r0 = circle.c;
     box.r0 = box.c;

 }

 void momentumTransfer(Box& box, Circle& circle,  coord &contactSurface,  coord &contactNorm, const float &mu)
 {
     
     
     float u1_norm = box.v * contactNorm, u1_par = box.v * contactSurface, v1_norm = circle.v * contactNorm, v1_par = circle.v * contactSurface;
     float d_par = (box.c-circle.c) * contactSurface, d_norm = ((box.c - circle.c) * contactNorm);
     //float mu = 0.8;
     float mct =  -mu * v1_par / circle.v.vecSize();

     float A1 = -d_par + d_norm * mct;
     float A2 = A1 * v1_norm + box.I * box.omega / circle.m;

     float B = v1_par - mct * v1_norm ;

     float C1 = circle.m * mct / box.m;
     float C2 = C1 * v1_norm + u1_par;
     float C3 = circle.m * v1_norm/box.m + u1_norm;

     float a = circle.m * circle.m / box.I * A1 * A1 + circle.m * (1 + mct * mct) + box.m * C1 * C1 + circle.m * circle.m / box.m;
     float b = 2 * circle.m * (( - circle.m / box.I) * A1 * A2 + B * mct - C3) - 2 * box.m * C2 * C1;
     float c = -box.I * box.omega * box.omega - circle.m * (v1_norm * v1_norm + v1_par * v1_par) - box.m * (u1_norm * u1_norm + u1_par * u1_par)
         + circle.m * (circle.m * A2 * A2 / box.I + B * B) + box.m * (C2 * C2 + C3 * C3);

     float v2_norm = (-b + sqrt(b * b - 4 * a * c)) / 2 / a;
     float v2_par = mct * (v2_norm - v1_norm) + v1_par;

     float u2_norm = (circle.m / box.m) * (v1_norm - v2_norm) + u1_norm;
     float u2_par = mct * (u2_norm - u1_norm) + u1_par;

     box.omega = box.omega + (circle.m / box.I) * ((v1_norm - v2_norm)*d_par+(v2_par-v1_par)*d_norm);


     circle.v = (contactNorm * v2_norm) + (contactSurface * v2_par);
     box.v = (contactNorm * u2_norm) + (contactSurface * u2_par);
     circle.r0 = circle.c;
     box.r0 = box.c;

 }

 bool getKey(float& v, float& beta, coord& r0, coord& r0_2, const float vmax)
 {

     if (GetAsyncKeyState(VK_UP) && r0.y < screenhight / 2) {
         r0.y += 2;
         r0_2.y = r0.y;
         return 0;
     }
     if (GetAsyncKeyState(VK_DOWN) && r0.y> -screenhight/2) {
         r0.y -= 2;
         r0_2.y = r0.y;
         return 0;
     }
     if (GetAsyncKeyState(VK_LEFT) && v>0) {
        v-=5;
         return 0;
     }
     if (GetAsyncKeyState(VK_RIGHT) && v<=vmax) {
         v+=5;
         return 0;
     }
     if ((GetAsyncKeyState('s') || GetAsyncKeyState('S')) && beta>-PI/2) {
         beta -= PI / 360;
         return 0;
     }
     if ((GetAsyncKeyState('w') || GetAsyncKeyState('W')) && beta<PI/2) {
         beta += PI / 360;
         return 0;
     }
     
     if (GetAsyncKeyState(VK_SPACE) || GetAsyncKeyState(VK_RETURN))
         return 1;

     return 0;
     
 }
 struct LevelVars
 {
     const float box_b[4] = {0.09, 0.07, 0.09, 0.09};
     const float box_h[4]={0.09, 0.07, 0.2 , 0.09};
     const float mu[4]={0.8, 0.8, 0.8, 0.8};
     const float m_cirle[4] = { 1, 1, 0.5, 1 };
     const float gate1[4]={0, 0.4, 0.1, 0.1};
     const float gate2[4]={0.3, 0.2, -0.1, -0.1};
     const float box_omega0[4]={0.0, 0.0, 0.0, 0.5};
     const coord box_v0[4] = { {0.0, 0.0}, {0.0, 0.0} , {0.0, 0.0}, {0.0, -100.0} };
     const coord box_r0[4] = { {0.0,0.0}, {0.0, 0.0} , {0.0, 0.0}, {0.0, 0.4} };

 };



int main()
{
    /* Initialize the library */
    if (!glfwInit())
        return -1;

    
    /* Create a windowed mode window and its OpenGL context */
    GLFWwindow* window = glfwCreateWindow(screenwidth, screenhight, "Shoot the Box", /*glfwGetPrimaryMonitor()*/NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        return -1;
    }

    /* Make the window's context current */
    glfwMakeContextCurrent(window);

    if (glewInit() != GLEW_OK)
        return -1;
    std::cout << glGetString(GL_VERSION) << std::endl;

    //glEnable(GL_DEBUG_OUTPUT);
    //glDebugMessageCallback(MessageCallback, 0);


    //glfwSwapInterval(1);

    coord r0 ((float) - screenwidth * 0.93 / 2, 0);


    LevelVars lv;
  
    Shader shader;
    glEnableVertexAttribArray(0);

    char level = 0;
    const char levels = sizeof(lv.mu)/sizeof(lv.mu[0]);

    while (level < levels)
    {
        Box box(lv.box_r0[level]*screenhight, lv.box_h[level] * (float)screenhight, lv.box_b[level] * (float)screenhight, 1,lv.box_omega0[level],lv.box_v0[level], & shader);
        box.colour = { 0.2f, 0.3f, 0.8f, 1.0f };

        Circle circle(r0, 8, 20, lv.m_cirle[level], {0.0,0.0}, &shader);
        circle.colour = { 0.2f, 0.8f, 0.2f, 1.0f };

        Line aimline(r0, screenwidth * 0.25, 0.0, &shader);
        aimline.colour = { 0.8f, 0.2f, 0.2f, 1.0f };

        coord contactSurface, contactNorm;

        float vsize = 200.0, vmax = 500.0;

        bool shoot = 0;


        float mu = lv.mu[level];

        coord gate[2] = { {(float)screenwidth / 2 - (float)screenmargin, (float)screenhight * lv.gate1[level]}, {(float)screenwidth / 2 - (float)screenmargin , (float)screenhight * lv.gate2[level]}};

        if (gate[0].y > gate[1].y)
        {
            coord tmp = gate[0];
            gate[0] = gate[1];
            gate[1] = tmp;
        }

        float v0boxdims[4] = { -(float)screenwidth * 0.75 * 0.5, (float)screenhight * 0.95 * 0.5 , (float)screenhight * 0.04, (float)screenwidth * 0.2 };

        Box v0outerbox({ v0boxdims[0], v0boxdims[1] }, v0boxdims[2], v0boxdims[3] + 6, -1.0, &shader);
        Box v0innerbox({ v0boxdims[0] - (v0boxdims[3] / 2) * (1 - vsize / vmax), v0boxdims[1] }, v0boxdims[2] * 0.85, v0boxdims[3] * vsize / vmax, 0.0, &shader);
        Box borderbox({ 0.0,0.0 }, (float)screenhight - (float)screenmargin * 2, (float)screenwidth - (float)screenmargin * 2, -1.0, &shader);
        Line gateline(gate[0], gate[1], &shader);
        v0innerbox.colour = { 0.8f, 0.3f, 0.8f, 1.0f };
        v0outerbox.colour = { 0.8f, 0.8f, 0.8f, 1.0f };
        borderbox.colour = { 1.0f, 0.1f, 0.1f, 1.0f };

        short hom = 0;

        auto t0 = std::chrono::high_resolution_clock::now();
        /* Loop until the user closes the window */
        while (/*!glfwWindowShouldClose(window) &&*/ !shoot && !hom)
        {

            glClear(GL_COLOR_BUFFER_BIT);

            box.draw();


            circle.draw();


            aimline.draw();
            v0innerbox.draw();
            v0outerbox.draw();
            borderbox.draw();
            gateline.draw();

            shoot = getKey(vsize, aimline.beta, circle.r0, aimline.r0, vmax);

            circle.setPositions(0.0);

            aimline.setPositions();

            v0innerbox.resize({ v0boxdims[0] - (v0boxdims[3] / 2) * (1 - vsize / vmax), v0boxdims[1] }, v0boxdims[2] * 0.85, (v0boxdims[3]) * vsize / vmax);

            
            std::chrono::duration<float> t = std::chrono::high_resolution_clock::now() - t0;
            box.setPositions(t.count());
            hom = hitormiss(box.positions, gate);
           
            /* Swap front and back buffers */
            glfwSwapBuffers(window);

            /* Poll for and process events */
            glfwPollEvents();

            if (GetAsyncKeyState(VK_ESCAPE) || glfwWindowShouldClose(window))
                glfwTerminate();
        }


        circle.v = { vsize * cos(aimline.beta), vsize * sin(aimline.beta) };

        if (box.v.vecSize())
        {
            std::chrono::duration<float> t = std::chrono::high_resolution_clock::now() - t0;
            box.timeReset(t.count());
        }
        t0 = std::chrono::high_resolution_clock::now();


        bool madecontact = 0;

        while (!glfwWindowShouldClose(window) && !hom)
        {

            glClear(GL_COLOR_BUFFER_BIT);

            box.draw();
            circle.draw();

            borderbox.draw();
            gateline.draw();


            std::chrono::duration<float> t = std::chrono::high_resolution_clock::now() - t0;

            box.setPositions(t.count());
            circle.setPositions(t.count());

            if (contact(box.positions, circle.c, circle.r, &contactSurface, &contactNorm))
            {
                auto t_tmp = std::chrono::high_resolution_clock::now();
                box.timeReset(t.count());
                t0 = t_tmp;
                momentumTransfer(box, circle, contactSurface, contactNorm, mu);
                madecontact = 1;

            }
            hom = hitormiss(box.positions, gate);

            if (!madecontact && !hom)                
                hom = bulletout(circle.c);

            /* Swap front and back buffers */
            glfwSwapBuffers(window);

            /* Poll for and process events */
            glfwPollEvents();

            if (GetAsyncKeyState(VK_ESCAPE))
                glfwTerminate();
        }

        if (hom == 1)
            level++;
    }
    glfwTerminate(); 
}
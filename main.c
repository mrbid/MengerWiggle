/*
    James William Fletcher (github.com/mrbid)
        December 2022

    L45: #define FUN
    - Comment this out and any deviation from the
    intended maxfps of 144 will cause the simulation
    speed to change. Going over the maximum frame
    rate the computer can render will slow down
    the simulation but upto that point it will be
    speeding up. Reducing the framerate from 144
    will also slow down the simulation.
    - Don't comment it out and the simulation speed
    should be stable at all framerates, even auto-
    correcting at high framerates. But less fun.

    Most people don't detect their delta time step
    they just use the current frame delta.

    #$%*(*(*&%^&*()&^%%$#$%^&*()(&*^%$#$$%^&^%$#%^&*^%$#%^&*()(*&^%$#@$%^&*()*&^%$#%^&*(
        original git: https://github.com/mrbid/MengerCube
        !!! DECEMBER UPDATE
            View and Normal matrix elements are being randomly expanded and contracted.
    }{:{:|}{(*&^%$@#$%^}|{:<>?<>^&*<>?@#$%^&{}|:"<><>?{}|:%^&*)_@#:"{@#$>:$%^&()@}|>_@#$
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <sys/file.h>
#include <stdint.h>
#include <unistd.h>

#pragma GCC diagnostic ignored "-Wunused-result"

#define uint GLushort
#define sint GLshort
#define f32 GLfloat

#include "inc/gl.h"
#define GLFW_INCLUDE_NONE
#include "inc/glfw3.h"

#ifndef __x86_64__
    #define NOSSE
#endif

#define SEIR_RAND       // uncommenting this define will enable the MMX random (it's a marginally slower)
//#define REGULAR_PHONG   // or Blinn-Phong by default
#define FUN             // uncomment this for stable simulation speed at different frame rates

#include "inc/esAux3.h"
#include "inc/res.h"
#include "ncube.h"

//*************************************
// globals
//*************************************
GLFWwindow* window;
uint winw = 1024;
uint winh = 768;
double t = 0;   // time
f32 dt = 0;     // delta time
double fc = 0;  // frame count
double lfct = 0;// last frame count time
f32 aspect;
double x,y,lx,ly;
double rww, ww, rwh, wh, ww2, wh2;
double uw, uh, uw2, uh2; // normalised pixel dpi
double maxfps = 144.0;

// render state id's
GLint projection_id;
GLint modelview_id;
GLint normalmat_id = -1;
GLint position_id;
GLint lightpos_id;
GLint solidcolor_id;
GLint color_id;
GLint opacity_id;
GLint normal_id; // 

// render state matrices
mat projection;
mat view;

// models
ESModel mdlMenger;

// camera vars
#define FAR_DISTANCE 333.f
uint focus_cursor = 0;
double sens = 0.001f;
f32 xrot = 0.f;
f32 yrot = d2PI; // face on until [1]
f32 zoom = -14.0f; // -6.0f / -26.0f

// sim vars
vec lightpos = {0.f, 0.f, 0.f};
f32 r=0.f,g=0.f,b=0.f;

//*************************************
// utility functions
//*************************************
void timestamp(char* ts)
{
    const time_t tt = time(0);
    strftime(ts, 16, "%H:%M:%S", localtime(&tt));
}
float urandf()
{
    static const float RECIP_FLOAT_UINT64_MAX = 1.f/(float)UINT64_MAX;
    int f = open("/dev/urandom", O_RDONLY | O_CLOEXEC);
    uint64_t s = 0;
    read(f, &s, sizeof(uint64_t));
    close(f);
    return ((float)s) * RECIP_FLOAT_UINT64_MAX;
}
float urandfc()
{
    static const float RECIP_FLOAT_UINT64_MAX = 2.f/(float)UINT64_MAX;
    int f = open("/dev/urandom", O_RDONLY | O_CLOEXEC);
    uint64_t s = 0;
    read(f, &s, sizeof(uint64_t));
    close(f);
    return (((float)s) * RECIP_FLOAT_UINT64_MAX)-1.f;
}
float clamp(float f, float min, float max)
{
    if(f > max){return max;}
    else if(f < min){return min;}
    return f;
}
void stepTitle(float speed)
{
    static uint m = 0;
    static uint p = 0;
    static double lt = 0.0;

    static char m1[32] = {0};
    static uint m1s = 0;
    if(m1s == 0 && m == 0)
    {
        sprintf(m1, "Fancy a wiggle?");
        m1s = strlen(m1);
    }
    else if(m1s == 0 && m == 1)
    {
        sprintf(m1, "Current speed %.2f", speed);
        m1s = strlen(m1);
    }
    
    if(t > lt)
    {
        if(p == 0)
        {
            glfwSetWindowTitle(window, "L3 Menger Cube");
            lt = t+6.0;
            p++;
            return;
        }
        else if(p > 0 && p < m1s)
        {
            char t[32] = {0};
            for(uint i = 0; i < p; i++)
                t[i] = m1[i];
            glfwSetWindowTitle(window, t);
            p++;
        }
        else if(p == m1s)
        {
            glfwSetWindowTitle(window, m1);
            lt  = t+6.0;
            p   = 0;
            m   = 1 - m;
            m1s = 0;
            return;
        }
        lt = t+0.09+(urandf()*0.04);
    }
}

//*************************************
// update & render
//*************************************
void main_loop(uint dotick)
{
//*************************************
// camera
//*************************************
    if(focus_cursor == 1)
    {
        glfwGetCursorPos(window, &x, &y);
        xrot += (ww2-x)*sens;
        yrot += (wh2-y)*sens;
        glfwSetCursorPos(window, ww2, wh2);
    }

    mIdent(&view);
    mTranslate(&view, 0.f, 0.f, zoom);
    mRotate(&view, yrot, 1.f, 0.f, 0.f);
    mRotate(&view, xrot, 0.f, 0.f, 1.f);
    
    static f32 ss = 0.08f;
    if(focus_cursor == 0 && dotick == 1)
    {
#ifdef FUN
        // this is not stable at different framerates
        static f32 tft = 0.f;
        tft += dt;
        yrot += sinf(tft*0.001f)*-ss;
        ss += dt*0.000001f;
#else
        // this is stable at different framerates
        static f32 tft = -1.3f;
        tft += dt*ss;
        yrot = sinf(tft)*100.f;
        ss += dt*0.001f;
#endif
        xrot += dt*0.01f;
        r += urandfc()*dt*1.6f;
        g += urandfc()*dt*1.6f;
        b += urandfc()*dt*1.6f;
        r = clamp(r, -1.f, 1.f);
        g = clamp(g, -1.f, 1.f);
        b = clamp(b, -1.f, 1.f);
        glUniform3f(color_id, r, g, b);
        const f32 ft = tft*0.5f;
        glUniform3f(lightpos_id, sinf(ft) * 10.0f, cosf(ft) * 10.0f, sinf(ft) * 10.0f);
        stepTitle(ss);
    }

//*************************************
// render
//*************************************
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    static int st = 0;
    if(st == 0){st = time(0);}

    int ts = st+(int)t;
    float frac = t-floorf(t);
    if(frac > 0.5f){frac -= 1.f; frac = fabsf(frac);}
    srand(ts);

    const float ws = esRandFloat(0.f, 3.f); //randf(&ts)*6.f;
    const uint mode = esRand(0,1);
    const uint iter = esRand(0,16);

    static int lp = 0.f;
    if(ts != lp)
    {
        printf(":: %u %u %.2f\n", mode, iter, ws);
        lp = ts;
    }

    for(uint i = 0; i < iter; i++)
    {
        if(mode == 0)
        {
            const uint r = esRand(0,3), c = esRand(0,3);
            view.m[r][c] += view.m[r][c]*frac;
        }
        else
            view.m[esRand(0,3)][esRand(0,3)] += (esRandFloat(-1.f, 1.f)*ws)*frac;
    }

    glUniformMatrix4fv(modelview_id, 1, GL_FALSE, (GLfloat*) &view.m[0][0]);
    if(normalmat_id != -1)
    {
        mat inverted, normalmat;
        mInvert(&inverted.m[0][0], &view.m[0][0]);
        mTranspose(&normalmat, &inverted);

        for(uint i = 0; i < iter; i++)
        {
            if(mode == 0)
            {
                const uint r = esRand(0,3), c = esRand(0,3);
                normalmat.m[r][c] += normalmat.m[r][c]*frac;
            }
            else
                normalmat.m[esRand(0,3)][esRand(0,3)] += (esRandFloat(-1.f, 1.f)*ws)*frac;
        }
        
        glUniformMatrix4fv(normalmat_id, 1, GL_FALSE, (GLfloat*) &normalmat.m[0][0]);
    }
    glDrawElements(GL_TRIANGLES, ncube_numind, GL_UNSIGNED_INT, 0);

    glfwSwapBuffers(window);
}

//*************************************
// Input Handelling
//*************************************
static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if(action == GLFW_PRESS)
    {
        if(key == GLFW_KEY_F)
        {
            if(t-lfct > 2.0)
            {
                char strts[16];
                timestamp(&strts[0]);
                const double nfps = fc/(t-lfct);
                printf("[%s] FPS: %g\n", strts, nfps);
                maxfps = nfps;
                dt = 1.0f / (float)maxfps;
                lfct = t;
                fc = 0;
            }
        }
        else if(key == GLFW_KEY_Z)
        {
            shadeLambert1(&position_id, &projection_id, &modelview_id, &lightpos_id, &normal_id, &color_id, &opacity_id);
            glUniformMatrix4fv(projection_id, 1, GL_FALSE, (GLfloat*) &projection.m[0][0]);
            glUniform3f(lightpos_id, lightpos.x, lightpos.y, lightpos.z);
            glUniform1f(opacity_id, 1.0f);
            glUniform3f(color_id, r, g, b);
            normalmat_id = -1;
        }
        else if(key == GLFW_KEY_X)
        {
            shadePhong1(&position_id, &projection_id, &modelview_id, &normalmat_id, &lightpos_id, &normal_id, &color_id, &opacity_id);
            glUniformMatrix4fv(projection_id, 1, GL_FALSE, (GLfloat*) &projection.m[0][0]);
            glUniform3f(lightpos_id, lightpos.x, lightpos.y, lightpos.z);
            glUniform1f(opacity_id, 1.0f);
            glUniform3f(color_id, r, g, b);
        }
        else if(key == GLFW_KEY_A)
            glDisable(GL_BLEND);
        else if(key == GLFW_KEY_S)
            glEnable(GL_BLEND);
        else if(key == GLFW_KEY_M)
        {
            int ts = time(0);
            srand(ts);
            projection.m[esRand(0,3)][esRand(0,3)] += randfc(&ts)*0.3f;
            glUniformMatrix4fv(projection_id, 1, GL_FALSE, (GLfloat*) &projection.m[0][0]);
        }
        else if(key == GLFW_KEY_N)
        {
            mIdent(&projection);
            mPerspective(&projection, 60.0f, aspect, 0.01f, FAR_DISTANCE); 
            glUniformMatrix4fv(projection_id, 1, GL_FALSE, (GLfloat*) &projection.m[0][0]);
        }
    }
}

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    if(yoffset == -1)
        zoom += 0.06f * zoom;
    else if(yoffset == 1)
        zoom -= 0.06f * zoom;
    
    if(zoom > 0.f){zoom = 0.f;}
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    if(action == GLFW_PRESS)
    {
        if(button == GLFW_MOUSE_BUTTON_LEFT)
        {
            focus_cursor = 1 - focus_cursor;
            if(focus_cursor == 0)
                glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
            else
                glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_HIDDEN);
            glfwSetCursorPos(window, ww2, wh2);
            glfwGetCursorPos(window, &ww2, &wh2);
        }
        else if(button == GLFW_MOUSE_BUTTON_RIGHT)
        {
            r = urandfc(), g = urandfc(), b = urandfc();
            glUniform3f(color_id, r, g, b);
        }
    }
}

void window_size_callback(GLFWwindow* window, int width, int height)
{
    winw = width;
    winh = height;

    glViewport(0, 0, winw, winh);
    aspect = (f32)winw / (f32)winh;
    ww = winw;
    wh = winh;
    rww = 1/ww;
    rwh = 1/wh;
    ww2 = ww/2;
    wh2 = wh/2;
    uw = (double)aspect / ww;
    uh = 1 / wh;
    uw2 = (double)aspect / ww2;
    uh2 = 1 / wh2;

    mIdent(&projection);
    mPerspective(&projection, 60.0f, aspect, 0.01f, FAR_DISTANCE); 
}

//*************************************
// Process Entry Point
//*************************************
int main(int argc, char** argv)
{
    // allow custom msaa level
    int msaa = 16;
    if(argc >= 2){msaa = atoi(argv[1]);}

    // allow framerate cap
    if(argc >= 3){maxfps = atof(argv[2]);}

    // help
    printf("----\n");
    printf("L3 Menger Cube\n");
    printf("----\n");
    printf("James William Fletcher (github.com/mrbid)\n");
    printf("----\n");
    printf("Argv(2): msaa, maxfps\n");
    printf("e.g; ./uc 16 60\n");
    printf("----\n");
    printf("Left Click = Focus toggle camera control\n");
    printf("Right Click = Random Colour\n");
    printf("F = FPS to console.\n");
    printf("A = Opaque.\n");
    printf("S = Transparent.\n");
    printf("Z = Lambertian Shading.\n");
    printf("X = Phong Shading.\n");
    printf("----\n");

    // init glfw
    if(!glfwInit()){printf("glfwInit() failed.\n"); exit(EXIT_FAILURE);}
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
    glfwWindowHint(GLFW_SAMPLES, msaa);
    window = glfwCreateWindow(winw, winh, "L3 Menger Cube", NULL, NULL);
    if(!window)
    {
        printf("glfwCreateWindow() failed.\n");
        glfwTerminate();
        exit(EXIT_FAILURE);
    }
    const GLFWvidmode* desktop = glfwGetVideoMode(glfwGetPrimaryMonitor());
    glfwSetWindowPos(window, (desktop->width/2)-(winw/2), (desktop->height/2)-(winh/2)); // center window on desktop
    glfwSetWindowSizeCallback(window, window_size_callback);
    glfwSetKeyCallback(window, key_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetScrollCallback(window, scroll_callback);
    glfwMakeContextCurrent(window);
    gladLoadGL(glfwGetProcAddress);
    glfwSwapInterval(0); // 0 for immediate updates, 1 for updates synchronized with the vertical retrace, -1 for adaptive vsync

    // set icon
    glfwSetWindowIcon(window, 1, &(GLFWimage){16, 16, (unsigned char*)&icon_image.pixel_data});

//*************************************
// projection
//*************************************

    window_size_callback(window, winw, winh);

//*************************************
// bind vertex and index buffers
//*************************************

    // ***** BIND MENGER *****
    esBind(GL_ARRAY_BUFFER, &mdlMenger.vid, ncube_vertices, sizeof(ncube_vertices), GL_STATIC_DRAW);
    esBind(GL_ARRAY_BUFFER, &mdlMenger.nid, ncube_normals, sizeof(ncube_normals), GL_STATIC_DRAW);
    esBind(GL_ELEMENT_ARRAY_BUFFER, &mdlMenger.iid, ncube_indices, sizeof(ncube_indices), GL_STATIC_DRAW);

//*************************************
// compile & link shader programs
//*************************************

    //makeAllShaders();
    makeLambert1();
    makePhong1();

//*************************************
// configure render options
//*************************************

    // standard stuff
    glBlendFunc(GL_SRC_ALPHA, GL_ONE);
    glEnable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.13f, 0.13f, 0.13f, 0.0f);

    // setup shader
    shadePhong1(&position_id, &projection_id, &modelview_id, &normalmat_id, &lightpos_id, &normal_id, &color_id, &opacity_id);
    glUniformMatrix4fv(projection_id, 1, GL_FALSE, (GLfloat*) &projection.m[0][0]);
    glUniform3f(lightpos_id, lightpos.x, lightpos.y, lightpos.z);
    glUniform1f(opacity_id, 0.5f);
    
    // bind menger to render
    r = urandf(), g = urandf(), b = urandf();
    glUniform3f(color_id, r, g, b);

    glBindBuffer(GL_ARRAY_BUFFER, mdlMenger.vid);
    glVertexAttribPointer(position_id, 3, GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(position_id);

    glBindBuffer(GL_ARRAY_BUFFER, mdlMenger.nid);
    glVertexAttribPointer(normal_id, 3, GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(normal_id);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mdlMenger.iid);

//*************************************
// execute update / render loop
//*************************************

    // init
    t = glfwGetTime();
    lfct = t;
    dt = 1.0f / (float)maxfps; // fixed timestep delta-time

#ifndef FUN
    glfwSetWindowTitle(window, "Detecting frame rate...");
    yrot = sinf(-1.3f)*100.f; // [1]
    time_t ac = time(0) + 1;
    uint fct = 0;
#else
    const uint fct = 1;
#endif
    
    // fps accurate event loop
    useconds_t wait_interval = 1000000 / maxfps; // fixed timestep
    if(wait_interval == 0){wait_interval = 100;} // limited to 10,000 FPS maximum
    useconds_t wait = wait_interval;
    while(!glfwWindowShouldClose(window))
    {
        usleep(wait);
        t = glfwGetTime();

#ifndef FUN
        // auto correct max fps
        if(time(0) > ac)
        {
            const double nfps = fc/(t-lfct);
            if(fabs(nfps - maxfps) > 6.f)
            {
                char strts[16];
                timestamp(&strts[0]);
                printf("[%s] maxfps auto corrected from %.2f to %.2f.\n", strts, maxfps, nfps);
            }
            maxfps = nfps;
            dt = 1.0f / (float)maxfps;
            ac = time(0) + 6;
            fct = 1;
        }
#endif
        
        // don't tick our internal state until we know we have a decent delta-time [dt]
        glfwPollEvents();
        main_loop(fct);

        // accurate fps
        wait = wait_interval - (useconds_t)((glfwGetTime() - t) * 1000000.0);
        if(wait > wait_interval)
            wait = wait_interval;
        //printf("%u: %u - %u\n", wait_interval, wait, (useconds_t)((glfwGetTime() - t) * 1000000.0));
        
        fc++;
    }

    // done
    glfwDestroyWindow(window);
    glfwTerminate();
    exit(EXIT_SUCCESS);
    return 0;
}


# CSCI 5611 Project 2

## Kailash Kalyanasundaram

## Project Description

For project 2, I created a cloth simulation similar to the one in class. The simulation environment is meant to be space with the a fire blanket (cloth) colliding against earth. In this simulation, I created a 3D cloth and a sphere that has been textured as the earth as seen in the videos. The following features have been implemented: Multiple Ropes with the 3d Cloth Simulation seen in the videos (At 4 seconds in the first video and 10 and 33 inn the second video), and 3D Simulation & Rendering as seen in the videos using the camera class and Vec3 class from class (rest of the video during the camera change) for a total of 90 points. 

### Project Difficulties

I attempted the Continuum Fluid Simulation but I couldn't get the simulation to work properly as it was really hard to render the live fluid for me personally. Additionally, I found it difficult to get the cloth textured like fire, whereas the sphere was fairly easier. Overall, I think it is super cool to see the actual cloth interacting with the earth and it was a great project!

## Simulation Videos

[Video 1](https://drive.google.com/file/d/1ncRnkxaumhR1FAEz3kXsNFzf-8AP2-GU/preview)

[Video 2](https://drive.google.com/file/d/1rzIvqh_Ylra2qHRR_5lTI-_yt4P6S4Cp/preview)

## Code

### ClothSimulation.pde

```
//CSCI 5611 Project 2
// Kailash Kalyanasundaram

//Simulation Parameters
int nx = 20;
int ny = 20;
float dt = 0.035;
float ks = 40; 
float kd = 25; 
float l0 = 10;
Vec3 gravity = new Vec3(0, 10, 0);
Camera camera;
PImage earth; 
PShape globe;

//Initial positions and velocities of masses

Vec3 p[][] = new Vec3[nx][ny];
Vec3 v[][] = new Vec3[nx][ny];
Vec3 vn[][] = new Vec3[nx][ny];
Vec3 acc[][] = new Vec3[nx][ny];
Vec3 spherePos = new Vec3(9.5 *nx - 10, 40, 90 - 10 * ny);
Vec3 sphereVel = new Vec3(0,0,0);
float sphereRadius = 20;
float speed = 10;

void update(float dt){
  
  vn = v;
  //(new velocity buffer) 
  Vec3 velBuffer = new Vec3(0,0,0);

//Update vels. before pos.
  sphereVel.setToLength(speed);
  spherePos.add(velBuffer.times(dt));

  
//for i in range(nx-1): #horiz. 
//for j in range(ny):
//e = p[i+1,j] - p[i,j]
//l = np.sqrt(e.dot(e))
//e = e/l #normalize
//v1 = e.dot(v[i,j])
//v2 = e.dot(v[i+1,j])
//f = -ks*(l0-l)-kd*(v1-v2)
//vn[i,j] += f*e*dt
//vn[i+1,j] -= f*e*dt

  for (int i = 0; i < nx-1; i++) {
    for (int j = 0; j < ny; j++) {
      Vec3 e = p[i+1][j].minus(p[i][j]);
      float l = e.length();
      e.normalize();
      float v1 = dot(e, v[i][j]);
      float v2 = dot(e, v[i+1][j]);
      float f = -ks*(l0-l)-kd*(v1-v2);
      vn[i][j].add(e.times(f*dt));
      vn[i+1][j].subtract(e.times(f*dt));
    }
  }
  //for i in range(nx):  #vertical     
//for j in range(ny-1):
//e = p[i,j+1] - p[i,j]
//l = np.sqrt(e.dot(e))
//e = e/l #normalize
//v1 = e.dot(v[i,j])
//v2 = e.dot(v[i,j+1])
//f = -ks*(l0-l)-kd*(v1-v2)
//vn[i,j] += f*e*dt
//vn[i,j+1] -= f*e*dt

  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny-1; j++) {
      Vec3 e = p[i][j+1].minus(p[i][j]);
      float l = e.length();
      e.normalize();
      float v1 = dot(e, v[i][j]);
      float v2 = dot(e, v[i][j+1]);
      float f = -ks*(l0-l)-kd*(v1-v2);
      vn[i][j].add(e.times(f*dt));
      vn[i][j+1].subtract(e.times(f*dt));
    }
  }
  //vn += [0,-.1,0,] #gravity
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      vn[i][j].add(gravity.times(dt));
    }
  }
  //vn[0,:] = 0      #fix top row
  for(int i = 0; i < ny; i++){
     vn[0][i] = new Vec3(0,0,0); 
    }
  //v = vn #update vel.
  v = vn;
  //p += v*dt #update pos.
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      p[i][j].add(v[i][j].times(dt));
    }
  }
//  for i in range(nx):
//for j in range(ny):
//d = SpherePos.distTo(p[i,j])
//if d < sphereR+.09:
//n = -1*(SpherePos - p[i,j]) #sphere normal
//n.normalize(); n = [n[0],n[1],n[2]]
//bounce = np.multiply(np.dot(v[i,j],n),n)
//v[i,j] -= 1.5*bounce
//p[i,j] += np.multiply(.1 + sphereR - d, n) #move out
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      float d = spherePos.distanceTo(p[i][j]);
      if (d < sphereRadius + 0.09) {
        Vec3 n = (spherePos.minus(p[i][j])).times(-1);
        n.normalize();
        Vec3 bounce = n.times((dot(v[i][j],n)));
        v[i][j].subtract(bounce.times(1.5));
        p[i][j].add(n.times(0.1+sphereRadius-d));
      }
    }
  }
  
}

void keyPressed()
{
  if (key == ' ')
    paused = !paused;
  camera.HandleKeyPressed();
}

void keyReleased()
{
  camera.HandleKeyReleased();
}


//Create Window
String windowTitle = "Swinging Cloth";
void setup() {
   for(int i = 0; i < nx; i++){
     for(int j = 0; j < ny; j++){
       p[i][j] = new Vec3(10 * i, 0, 0 - 10 * j);
       v[i][j] = new Vec3(0, 0, 0);
       vn[i][j] = new Vec3(0,0,0);
       acc[i][j] = new Vec3(0,0,0);
     }
   }
  size(1024, 768, P3D);
  camera = new Camera();
  String http = "http://";
  earth = loadImage( http + "previewcf.turbosquid.com/Preview/2014/08/01__15_41_30/Earth.JPG5a55ca7f-1d7c-41d7-b161-80501e00d095Larger.jpg");
  globe = createShape(SPHERE, sphereRadius); 
  globe.setTexture(earth);
  frameRate(30);
  surface.setTitle(windowTitle);
}

void followMouse() {
  Vec3 mouse = new Vec3(mouseX,mouseY,0);
  Vec3 diff = mouse.minus(spherePos);
  spherePos = diff;
  sphereVel.add(diff.times(dt));
}

//Draw the scene: one sphere per mass, one line connecting each pair
boolean paused = true;
void draw() {
  background(0,0,0);
    if (!paused)
       update(dt);
    noStroke();
    lights();
    camera.Update(1.0/frameRate);
    pushMatrix();
    fill(255, 0, 0);
    translate(spherePos.x, spherePos.y, spherePos.z);
    //Received from https://forum.processing.org/two/discussion/13500/applying-a-texture-to-a-sphere
    shape(globe);
    //sphere(sphereRadius);
    if (key == 'm') {
       followMouse();
    }
    popMatrix();
    fill(0,0,0);
    //for(int i = 0; i < nx; i++){
    //  for(int j = 0; j < ny; j++){
    //    pushMatrix();
    //    fill(0,0,255);
    //    translate(p[i][j].x, p[i][j].y, p[i][j].z);
    //    sphere(2);
    //    popMatrix();
    //  }
    //}
    stroke(#FF4500);
    for(int i = 0; i < nx; i++){
      for(int j = 0; j < ny - 1; j++){
        line(p[i][j].x, p[i][j].y, p[i][j].z, p[i][j+1].x, p[i][j+1].y, p[i][j+1].z);
       
      }
    }
    for(int i = 0; i < nx - 1; i++){
      for(int j = 0; j < ny; j++){
        line(p[i][j].x, p[i][j].y, p[i][j].z, p[i+1][j].x, p[i+1][j].y, p[i+1][j].z);
      }
    }
}


///////////////////
// Vec3D Library
///////////////////

public class Vec3 {
  public float x, y, z;
  
  public Vec3(float x, float y, float z){
    this.x = x;
    this.y = y;
    this.z = z;
  }
  
  public String toString(){
    return "(" + x+ ", " + y + ", " + z +")";
  }
  
  public float length(){
    return sqrt(x*x+y*y+z*z);
  }
  
  public float lengthSqr(){
    return x*x+y*y+z*z;
  }
  
  public Vec3 plus(Vec3 rhs){
    return new Vec3(x+rhs.x, y+rhs.y,z+rhs.z);
  }
  
  public void add(Vec3 rhs){
    x += rhs.x;
    y += rhs.y;
    z += rhs.z;
  }
  
  public Vec3 minus(Vec3 rhs){
    return new Vec3(x-rhs.x, y-rhs.y, z-rhs.z);
  }
  
  public void subtract(Vec3 rhs){
    x -= rhs.x;
    y -= rhs.y;
    z -= rhs.z;
  }
  
  public Vec3 times(float rhs){
    return new Vec3(x*rhs, y*rhs, z*rhs);
  }
  
  public void mul(float rhs){
    x *= rhs;
    y *= rhs;
    z *= rhs;
  }
  
  public void normalize(){
    float magnitude = sqrt(x*x + y*y + z*z);
    x /= magnitude;
    y /= magnitude;
    z /= magnitude;
  }
  
  public Vec3 normalized(){
    float magnitude = sqrt(x*x + y*y + z*z);
    return new Vec3(x/magnitude, y/magnitude,z/magnitude);
  }
  
  public void clampToLength(float maxL){
    float magnitude = sqrt(x*x + y*y + z*z);
    if (magnitude > maxL){
      x *= maxL/magnitude;
      y *= maxL/magnitude;
      z *= maxL/magnitude;
    }
  }
  
  public void setToLength(float newL){
    float magnitude = sqrt(x*x + y*y + z*z);
    x *= newL/magnitude;
    y *= newL/magnitude;
    z *= newL/magnitude;
  }
  
  public float distanceTo(Vec3 rhs){
    float dx = rhs.x - x;
    float dy = rhs.y - y;
    float dz = rhs.z - z;
    return sqrt(dx*dx + dy*dy + dz*dz);
  }
  
}

Vec3 interpolate(Vec3 a, Vec3 b, float t){
  return a.plus((b.minus(a)).times(t));
}

float interpolate(float a, float b, float t){
  return a + ((b-a)*t);
}

float dot(Vec3 a, Vec3 b){
  return a.x*b.x + a.y*b.y + a.z*b.z;
}

Vec3 projAB(Vec3 a, Vec3 b){
  return b.times(a.x*b.x + a.y*b.y + a.z*b.z);
}
```
### Camera.pde

```
// Created for CSCI 5611 by Liam Tyler

// WASD keys move the camera relative to its current orientation
// Arrow keys rotate the camera's orientation
// Holding shift boosts the move speed

class Camera
{
  Camera()
  {
    position      = new PVector( 0, 0, 0 ); // initial position
    theta         = 0; // rotation around Y axis. Starts with forward direction as ( 0, 0, -1 )
    phi           = 0; // rotation around X axis. Starts with up direction as ( 0, 1, 0 )
    moveSpeed     = 50;
    turnSpeed     = 1.57; // radians/sec
    boostSpeed    = 10;  // extra speed boost for when you press shift
    
    // dont need to change these
    shiftPressed = false;
    negativeMovement = new PVector( 0, 0, 0 );
    positiveMovement = new PVector( 0, 0, 0 );
    negativeTurn     = new PVector( 0, 0 ); // .x for theta, .y for phi
    positiveTurn     = new PVector( 0, 0 );
    fovy             = PI / 4;
    aspectRatio      = width / (float) height;
    nearPlane        = 0.1;
    farPlane         = 10000;
  }
  
  void Update(float dt)
  {
    theta += turnSpeed * ( negativeTurn.x + positiveTurn.x)*dt;
    
    // cap the rotation about the X axis to be less than 90 degrees to avoid gimble lock
    float maxAngleInRadians = 85 * PI / 180;
    phi = min( maxAngleInRadians, max( -maxAngleInRadians, phi + turnSpeed * ( negativeTurn.y + positiveTurn.y ) * dt ) );
    
    // re-orienting the angles to match the wikipedia formulas: https://en.wikipedia.org/wiki/Spherical_coordinate_system
    // except that their theta and phi are named opposite
    float t = theta + PI / 2;
    float p = phi + PI / 2;
    PVector forwardDir = new PVector( sin( p ) * cos( t ),   cos( p ),   -sin( p ) * sin ( t ) );
    PVector upDir      = new PVector( sin( phi ) * cos( t ), cos( phi ), -sin( t ) * sin( phi ) );
    PVector rightDir   = new PVector( cos( theta ), 0, -sin( theta ) );
    PVector velocity   = new PVector( negativeMovement.x + positiveMovement.x, negativeMovement.y + positiveMovement.y, negativeMovement.z + positiveMovement.z );
    position.add( PVector.mult( forwardDir, moveSpeed * velocity.z * dt ) );
    position.add( PVector.mult( upDir,      moveSpeed * velocity.y * dt ) );
    position.add( PVector.mult( rightDir,   moveSpeed * velocity.x * dt ) );
    
    aspectRatio = width / (float) height;
    perspective( fovy, aspectRatio, nearPlane, farPlane );
    camera( position.x, position.y, position.z,
            position.x + forwardDir.x, position.y + forwardDir.y, position.z + forwardDir.z,
            upDir.x, upDir.y, upDir.z );
  }
  
  // only need to change if you want difrent keys for the controls
  void HandleKeyPressed()
  {
    if ( key == 'w' || key == 'W' ) positiveMovement.z = 1;
    if ( key == 's' || key == 'S' ) negativeMovement.z = -1;
    if ( key == 'a' || key == 'A' ) negativeMovement.x = -1;
    if ( key == 'd' || key == 'D' ) positiveMovement.x = 1;
    if ( key == 'q' || key == 'Q' ) positiveMovement.y = 1;
    if ( key == 'e' || key == 'E' ) negativeMovement.y = -1;
    
    if ( key == 'r' || key == 'R' ){
      Camera defaults = new Camera();
      position = defaults.position;
      theta = defaults.theta;
      phi = defaults.phi;
    }
    
    if ( keyCode == LEFT )  negativeTurn.x = 1;
    if ( keyCode == RIGHT ) positiveTurn.x = -0.5;
    if ( keyCode == UP )    positiveTurn.y = 0.5;
    if ( keyCode == DOWN )  negativeTurn.y = -1;
    
    if ( keyCode == SHIFT ) shiftPressed = true; 
    if (shiftPressed){
      positiveMovement.mult(boostSpeed);
      negativeMovement.mult(boostSpeed);
    }
    
  }
  
  // only need to change if you want difrent keys for the controls
  void HandleKeyReleased()
  {
    if ( key == 'w' || key == 'W' ) positiveMovement.z = 0;
    if ( key == 'q' || key == 'Q' ) positiveMovement.y = 0;
    if ( key == 'd' || key == 'D' ) positiveMovement.x = 0;
    if ( key == 'a' || key == 'A' ) negativeMovement.x = 0;
    if ( key == 's' || key == 'S' ) negativeMovement.z = 0;
    if ( key == 'e' || key == 'E' ) negativeMovement.y = 0;
    
    if ( keyCode == LEFT  ) negativeTurn.x = 0;
    if ( keyCode == RIGHT ) positiveTurn.x = 0;
    if ( keyCode == UP    ) positiveTurn.y = 0;
    if ( keyCode == DOWN  ) negativeTurn.y = 0;
    
    if ( keyCode == SHIFT ){
      shiftPressed = false;
      positiveMovement.mult(1.0/boostSpeed);
      negativeMovement.mult(1.0/boostSpeed);
    }
  }
  
  // only necessary to change if you want different start position, orientation, or speeds
  PVector position;
  float theta;
  float phi;
  float moveSpeed;
  float turnSpeed;
  float boostSpeed;
  
  // probably don't need / want to change any of the below variables
  float fovy;
  float aspectRatio;
  float nearPlane;
  float farPlane;  
  PVector negativeMovement;
  PVector positiveMovement;
  PVector negativeTurn;
  PVector positiveTurn;
  boolean shiftPressed;
};

// ----------- Example using Camera class -------------------- //
//Camera camera;

//void keyReleased()
//{
//  camera.HandleKeyReleased();
//}
```



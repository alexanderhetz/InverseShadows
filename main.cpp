/* Trying to produce objects casting certain shadows */

#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <time.h>
#include <bits/stdc++.h>

using namespace std;

const double PI=M_PI;
const int PIXELS=201;
const double BOXSIDE=10;
const double MINIMUMRADIUS=BOXSIDE/(2*PIXELS);
const double PN=.0001; // Regularization parameter. Cost of adding a sphere
const double A=.05; // Displacement amplitude for coordinates and radii
const double B=.3; // Probability of taking away a sphere
const double C=.7; // Probability of adding a sphere
const double D=.25; // Probability of adding or deleting spheres again
const int N=4;//4; // Pairs of mothers
const int M=32;//16; // Pairs of daughters
const int NIMAGES=6;

    double UniRand(){
        return (double)rand()/RAND_MAX;
    }

    double RandomRadius(){
        return BOXSIDE*((double)rand()/RAND_MAX)/50+MINIMUMRADIUS;
    }

    double RandomInBox(){
        return BOXSIDE*((double)rand()/RAND_MAX)-BOXSIDE/2;
    }

class Sphere{
public:
    double center[3];
    double radius;

    Sphere(){
    center[0]=0;
    center[1]=0;
    center[2]=0;
    radius=0;
    }

    Sphere(double c[3],double r){
        center[0]=c[0];
        center[1]=c[1];
        center[2]=c[2];
        radius=abs(r);
    }

    Sphere(double x,double y,double z,double r){
        center[0]=x;
        center[1]=y;
        center[2]=z;
        radius=abs(r);
    }

    void PrintSphere(){
    cout<<"{{"<<center[0]<<","<<center[1]<<","<<center[2]<<"},"<<radius<<"}";
    }

    void PrintSphereToFile(ofstream &newfile){
    newfile<<"{{"<<center[0]<<","<<center[1]<<","<<center[2]<<"},"<<radius<<"}";
    }
};

class Circle{
public:
    double center[2];
    double radius;

    Circle(){
    center[0]=0;
    center[1]=0;
    radius=0;
    }

    Circle(double c[2],double r){
        center[0]=c[0];
        center[1]=c[1];
        radius=abs(r);
    }

    Circle(double x,double y,double r){
        center[0]=x;
        center[1]=y;
        radius=abs(r);
    }

    Circle(Sphere s, double th, double ph){ //Proyection of a sphere in the angles th and ph. th is the angle of the radial vector with respect to the z-axis
                                            //ph is the angle of the radial vector with respect to the x-axis. The proyection plane is determined by the radial
                                            // vector and by a vector in the x-y plane.
        center[0]=s.center[0]*sin(ph)-s.center[1]*cos(ph);
        center[1]=s.center[0]*sin(th)*cos(ph)+s.center[1]*sin(th)*sin(ph)+s.center[2]*cos(th);
        radius=s.radius;
    }

    void PrintCircle(){
        cout<<"{{"<<center[0]<<","<<center[1]<<"}, "<<radius<<"}";
    }
};

class Image{ // Square binary images represented by booleans
public:
    bool canvas[PIXELS][PIXELS];

    Image(){
    for(int i=0;i<PIXELS;i++){
        for(int j=0;j<PIXELS;j++){
            canvas[i][j]=false;
        }
    }
    }

    void DrawCircle(Circle cir){ //Bresenham’s circle drawing algorithm modified to fill the circle
    int r=round(cir.radius*PIXELS/BOXSIDE);
    int x0=round((cir.center[0]+BOXSIDE/2)*PIXELS/BOXSIDE);
    int y0=round((cir.center[1]+BOXSIDE/2)*PIXELS/BOXSIDE);

    int f=1-r;
    int ddfx=0;
    int ddfy=-2*r;
    int x=0;
    int y=r;

    for(int i=-r;i<=r;i++){
        if(x0>=0 && x0<=PIXELS-1 && y0+i>=0 && y0+i<=PIXELS-1){canvas[x0][y0+i]=true;}
    }
    if(x0+r>=0 && x0+r<=PIXELS-1 && y0>=0 && y0<=PIXELS-1){canvas[x0+r][y0]=true;}
    if(x0-r>=0 && x0-r<=PIXELS-1 && y0>=0 && y0<=PIXELS-1){canvas[x0-r][y0]=true;}

    while(x<y){
        if(f>0){
            y--;
            ddfy+=2;
            f+=ddfy;
        }
        x++;
        ddfx+=2;
        f+=ddfx+1;

        for(int i=-y;i<=y;i++){
            if(x0+x>=0 && x0+x<=PIXELS-1 && y0+i>=0 && y0+i<=PIXELS-1){canvas[x0+x][y0+i]=true;}
            if(x0-x>=0 && x0-x<=PIXELS-1 && y0+i>=0 && y0+i<=PIXELS-1){canvas[x0-x][y0+i]=true;}
        }
        for(int i=-x;i<=x;i++){
            if(x0+y>=0 && x0+y<=PIXELS-1 && y0+i>=0 && y0+i<=PIXELS-1){canvas[x0+y][y0+i]=true;}
            if(x0-y>=0 && x0-y<=PIXELS-1 && y0+i>=0 && y0+i<=PIXELS-1){canvas[x0-y][y0+i]=true;}
        }
    }

    }

    void ImageToText(string title){ //Prints a text representation of the image, only for viewing purposes
    ofstream newfile;
    newfile.open (title);

        for(int i=0;i<PIXELS;i++){
        for(int j=0;j<PIXELS;j++){
            newfile << canvas[i][j] << " ";
        }
            newfile << "\n";
        }

    newfile.close();
    }

    void StringToImage(string title){ //Loads an image from a file formed by consecutive ones and zeroes
        ifstream file;
        file.open (title,ios::in);
        string line;
        getline(file,line);
        file.close();

        char charim[PIXELS][PIXELS];
        for(int i=0;i<PIXELS;i++){
        for(int j=0;j<PIXELS;j++){
            charim[i][j]=line[PIXELS*i+j];
        }
        }
        for(int i=0;i<PIXELS;i++){
        for(int j=0;j<PIXELS;j++){
            if(charim[i][j]=='1'){
            canvas[i][j]=true;
            }else{
            canvas[i][j]=false;}
        }
    }
    }

};

int SquaredImageDistance(Image im1,Image im2){
    int d=0;
    for(int i=0;i<PIXELS;i++){
        for(int j=0;j<PIXELS;j++){
            if(im1.canvas[i][j]!=im2.canvas[i][j]){
                d++;
            }
        }

        }
    return d;
}

struct Node{
public:
    Sphere sphere;
    Node* nextsphere;
};

class Model{ //Set of spheres represented as a linked list of spheres and other features

public:
    Node* head;
    int length;
    double fitness;

    bool operator < (const Model& m) const{ //Overloading of the < operator
        return (fitness < m.fitness);
    }

    Model(){
    head=NULL;
    length=0;
    fitness=NAN;
    }

    void AddSphere(Sphere s){
    Node* newnode=new Node;
    newnode->sphere=s;
    newnode->nextsphere=NULL;

    if(head==NULL){
        head=newnode;
    }else{
        newnode->nextsphere=head;
        head=newnode;
    }

    length++;
    }

    void AddSphereToLast(Sphere s){
    Node* temp=head;
    if(head==NULL){
        AddSphere(s);
    }else{while(temp->nextsphere!=NULL){
        temp=temp->nextsphere;
    }
    Node* newnode=new Node;
    newnode->sphere=s;
    newnode->nextsphere=NULL;
    temp->nextsphere=newnode;
    length++;
    }
    }

    int Length(){
    return length;
    }

    double Fitness(){
        return fitness;
    }

    Node* Head(){
        return head;
    }

    Sphere NthSphere(int n){

    Node* temp=head;
    for(int i=1;i<n;i++){
        temp=temp->nextsphere;
    }
    return temp->sphere;
    }

    void DeleteNthSphere(int n){

    if(n==1){
        delete head;
        head=head->nextsphere;
    }else{
        Node* temp=head;
        for(int i=1;i<n-1;i++){
            temp=temp->nextsphere;
        }
        Node* temppt=temp->nextsphere->nextsphere;
        delete temp->nextsphere;
        temp->nextsphere=temppt;
    }
    length--;
    }

    void DeleteAllSpheres(){ //This to avoid memory leaks
    Node* temp;
    while (head != NULL)
    {
    temp = head->nextsphere;
    delete head;
    head = temp;
    }
    length=0;
    }

    void PrintModelToTable(){
    Node* temp=head;
    cout<<"{";
    if(head!=NULL){
        while(temp->nextsphere!=NULL){
            temp->sphere.PrintSphere();
            cout<<",";
            temp=temp->nextsphere;
    }
    temp->sphere.PrintSphere();
    }
    cout<<"}";
    }

    void PrintModelToTableToFile(string title){
    Node* temp=head;
    ofstream newfile;
    newfile.open (title);
    newfile<<"{";
    if(head!=NULL){
        while(temp->nextsphere!=NULL){
            temp->sphere.PrintSphereToFile(newfile);
            newfile<<",";
            temp=temp->nextsphere;
    }
    temp->sphere.PrintSphereToFile(newfile);
    }
    newfile<<"}";
    newfile.close();
    }

    Image ProjectModelToImage(double th, double ph){
        Image im;
        Node* temp=head;
        while(temp!=NULL){
            im.DrawCircle(Circle(temp->sphere,th,ph));
            temp=temp->nextsphere;
        }
        return im;
    }

    void CalculateFitness(Image images[NIMAGES],double angles[NIMAGES][2]){
    Image modelims[NIMAGES];
    for(int i=0;i<NIMAGES;i++){
        modelims[i]=ProjectModelToImage(angles[i][0],angles[i][1]);
    }
    double fit=0;
    for(int i=0;i<NIMAGES;i++){
        fit+=SquaredImageDistance(modelims[i],images[i]);
    }
    fit=sqrt(fit)+length*PN;
    fitness=fit;
    }

    void GenerateRandomModel(int n){
        double r;
        for(int i=1;i<=n;i++){
            r=RandomRadius();
            Sphere s=Sphere(RandomInBox()*(BOXSIDE-r)/BOXSIDE,RandomInBox()*(BOXSIDE-r)/BOXSIDE,RandomInBox()*(BOXSIDE-r)/BOXSIDE,r);
            AddSphere(s);
        }
    }

};

class Madres{
public:
    Model generation[2*N*M+2];

    Madres(){
        for(int i=0;i<2*N*M+2;i++){
            generation[i]=Model();
        }
    }

    void GenerateRandomGeneration(int n){
        for(int i=0;i<2*N*M+2;i++){
        generation[i].GenerateRandomModel(n);
        }
    }

    void ImportGeneration(string title){
        Madres m=Madres();
        ifstream file;
        file.open (title,ios::in);
        string line;

        double c[3];
        double r;
        int l=0;

        for(int i=0;i<2*N*M+2;i++){
        getline(file,line,' ');
        m.generation[i].fitness=stod(line);

        getline(file,line,' ');
        l=stoi(line);

            for(int j=0;j<l;j++){
                for(int k=0;k<3;k++){
                    getline(file,line,' ');
                    c[k]=stod(line);
                }
                getline(file,line,' ');
                r=stod(line);

                m.generation[i].AddSphereToLast(Sphere(c,r));
            }

        }
        file.close();
        for(int i=0;i<2*N*M+2;i++){
            generation[i]=m.generation[i];
        }
    }

    void ExportGeneration(string title){
        ofstream newfile;
        newfile.open (title);
        int l=0;

        for(int i=0;i<2*N*M+2;i++){
        l=generation[i].Length();
        newfile<<generation[i].Fitness()<<" ";
        newfile<<l<<" ";

        Node* temp=generation[i].Head();
            for(int j=0;j<l;j++){
                for(int k=0;k<3;k++){
                    newfile<<temp->sphere.center[k]<<" ";
                }
                newfile<<temp->sphere.radius<<" ";
                temp=temp->nextsphere;
            }

        }
        newfile.close();
    }

    void GenerationFitness(Image images[NIMAGES],double angles[NIMAGES][2]){
        for(int i=0;i<2*N*M+2;i++){
            generation[i].CalculateFitness(images,angles);
        }
    }

    void SortByFitness(){
        sort(generation,generation+2*N*M+2);
    }

    void Hijas(Image images[NIMAGES],double angles[NIMAGES][2]){
        Madres m=Madres();

        m.generation[2*N*M]=generation[0]; //Keep two best individuals from the previous generation
        m.generation[2*N*M+1]=generation[1];

        Node* temp1;
        Node* temp2;
        for(int i=0;i<2*M;i++){ //Randomly take spheres from each mother to assemble a daughter
            for(int j=0;j<N;j++){
                temp1=generation[j].head;
                temp2=generation[j+N].head;
                while(temp1!=NULL&&temp2!=NULL){
                    if(UniRand()>=.5){
                        m.generation[N*i+j].AddSphereToLast(temp1->sphere);
                    }else{
                        m.generation[N*i+j].AddSphereToLast(temp2->sphere);
                    }
                temp1=temp1->nextsphere;
                temp2=temp2->nextsphere;
                }
                if(temp1==NULL){ //Only if the longest mother has better fitness her extra spheres are added
                     if(generation[j].fitness>=generation[j+N].fitness){
                        while(temp2!=NULL){
                            m.generation[N*i+j].AddSphereToLast(temp2->sphere);
                            temp2=temp2->nextsphere;
                        }
                     }
                }else if(temp2==NULL){
                    if(generation[j+N].fitness>=generation[j].fitness){
                        while(temp1!=NULL){
                            m.generation[N*i+j].AddSphereToLast(temp1->sphere);
                            temp1=temp1->nextsphere;
                        }
                    }
                }
            }
        }

        if(UniRand()>=.5){
            for(int i=0;i<2*N*M;i++){ //Randomly delete a sphere with a certain probability
                if(B>=UniRand()){
                    m.generation[i].DeleteNthSphere(rand()%m.generation[i].length+1);
                }
            }
            double r;
            for(int i=0;i<2*N*M;i++){ //Randomly add a new random sphere with a certain probability
                if(C>=UniRand()){
                    r=RandomRadius();
                    m.generation[i].AddSphereToLast(Sphere(RandomInBox()*(BOXSIDE-r)/BOXSIDE,RandomInBox()*(BOXSIDE-r)/BOXSIDE,RandomInBox()*(BOXSIDE-r)/BOXSIDE,r));
                }
            }

            if(D>=UniRand()){ //Randomly add and delete more spheres, but with less probability
                for(int i=0;i<2*N*M;i++){ //Randomly delete a sphere with a certain probability
                    if(B>=UniRand()){
                        m.generation[i].DeleteNthSphere(rand()%m.generation[i].length+1);
                    }
                }
                for(int i=0;i<2*N*M;i++){ //Randomly add a new random sphere with a certain probability
                    if(C>=UniRand()){
                    r=RandomRadius();
                    m.generation[i].AddSphereToLast(Sphere(RandomInBox()*(BOXSIDE-r)/BOXSIDE,RandomInBox()*(BOXSIDE-r)/BOXSIDE,RandomInBox()*(BOXSIDE-r)/BOXSIDE,r));
                    }
                }
            }

        }else{
            Node* temp;
            double newcoor;
            for(int i=0;i<2*N*M;i++){ //Randomly vary the positions and radii of spheres with probability inversely proportional to number of spheres
                temp=m.generation[i].head;
                while(temp!=NULL){
                    for(int k=0;k<3;k++){
                        newcoor=temp->sphere.center[k]+A*(2*UniRand()-1);
                        if((1.0/generation[0].length)>=UniRand()&&abs(newcoor)+temp->sphere.radius<=(BOXSIDE/2)){
                            temp->sphere.center[k]=newcoor;
                        }
                    }
                    newcoor=temp->sphere.radius+A*(2*UniRand()-1);
                    if((1.0/generation[0].length)>=UniRand()&&newcoor>=MINIMUMRADIUS){
                        temp->sphere.radius=newcoor;
                    }
                    temp=temp->nextsphere;
                }
            }
        }

        m.GenerationFitness(images,angles);
        m.SortByFitness();

        for(int i=2;i<2*N*M+2;i++){ // This to avoid memory leaks. Since we copied the first two "generation" at the beginning, we can't delete them now
            generation[i].DeleteAllSpheres();
        }

        for(int i=0;i<2*N*M+2;i++){
            generation[i]=m.generation[i];
        }
    }

};

int main()
{
    srand(time(0));
    time_t be,en;
    time(&be);

    //double angles[NIMAGES][2]={{0,0},{0,PI/2},{PI/2,0}};
    //double angles[NIMAGES][2]={{0,0},{0,2*PI/3},{0,4*PI/3}};
    //double angles[NIMAGES][2]={{PI/2,0},{PI/2+acos(-1.0/3),0},{PI/2+acos(-1.0/3),2*PI/3},{PI/2+acos(-1.0/3),4*PI/3}};
    //double angles[NIMAGES][2]={{0, 0}, {0, PI/4}, {0, PI/2}, {0, 3*PI/4}};
    //double angles[NIMAGES][2]={{0,0},{0,PI/2}};
    //double angles[NIMAGES][2]={{0,0},{0,PI/3},{0,2*PI/3}};
    //double angles[NIMAGES][2]={{0,0},{0,PI/5},{0,2*PI/5},{0,3*PI/5},{0,4*PI/5}};
    //double angles[NIMAGES][2] ={{0,0},{PI/4,PI/3},{0,2*PI/3},{PI/4,PI},{0,4*PI/3},{PI/4,5*PI/3}};
    double angles[NIMAGES][2] = {{PI/2 - atan(4.0/(3 + sqrt(5))),
    0}, {PI/2 - atan(pow(1 + sqrt(5),2)/2),
    PI/3}, {PI/2 - atan(4.0/(3 + sqrt(5))),
    2*PI/3}, {PI/2 - atan(pow(1 + sqrt(5),2)/2),
    PI}, {PI/2 - atan(4.0/(3 + sqrt(5))),
    4*PI/3}, {PI/2 - atan(pow(1 + sqrt(5),2)/2), 5*PI/3}};

    Image images[NIMAGES];
    images[0].StringToImage("c0.txt");
    images[1].StringToImage("c1.txt");
    images[2].StringToImage("c2.txt");
    images[3].StringToImage("c3.txt");
    images[4].StringToImage("c4.txt");
    images[5].StringToImage("c5.txt");

    Madres m;
    //m.GenerateRandomGeneration(64);
    //m.GenerationFitness(images,angles);
    //m.SortByFitness();
    m.ImportGeneration("generation.txt");
    for(int i=1;i<=4*1024;i++){
    m.Hijas(images,angles);
    time(&en);
    cout<<i<<": "<<m.generation[0].fitness<<" | "<<m.generation[0].length<<" | "<<100*pow(m.generation[0].fitness-PN*m.generation[0].length,2)/(NIMAGES*pow(PIXELS,2))<<" | "<<en-be<<"s"<<endl;
    }
    m.ExportGeneration("generation.txt");
    m.generation[0].PrintModelToTableToFile("model.txt");

    return 0;
}

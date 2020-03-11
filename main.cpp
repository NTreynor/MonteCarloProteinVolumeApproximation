

// Created by Nicholas Treynor 2/24

#include <iostream>
#include <stdio.h>      /* printf */
#include <math.h>       /* sqrt */
#include <fstream>
#include <sstream>
#include <string>
#include <ctime>
#include <cstdlib>

using namespace std;

class Point {
public:
    float getX() {
        return x;
    }

    float getY() {
        return y;
    }

    float getZ() {
        return z;
    }

    Point(float initX, float initY, float initZ) {
        x = initX;
        y = initY;
        z = initZ;
    }

private:
    float x;
    float y;
    float z;




};

class Ball {
public:
    float getX() {
        return x;
    }

    float getY() {
        return y;
    }

    float getZ() {
        return z;
    }

    float maxX() {
        return x + radius;
    }

    float maxY() {
        return y + radius;
    }

    float maxZ() {
        return z + radius;
    }

    float minX() {
        return x - radius;
    }

    float minY() {
        return y - radius;
    }

    float minZ() {
        return z - radius;
    }

    float getRadius() {
        return radius;
    }

    bool contains(Point point) {
        float pX = point.getX();
        float pY = point.getY();
        float pZ = point.getZ();
// Distance = d = ((x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2)^½
        if ((std::sqrt(((pX - this->x) * (pX - this->x)) + ((pY - this->y) * (pY - this->y)) + ((pZ - this->z) * (pZ - this->z)))) <= this->radius){ // You can remove the SQRT
            return true;
        } else {
            return false;
        }
        return false;
    }

    float distance(Point point) {
        float pX = point.getX();
        float pY = point.getY();
        float pZ = point.getZ();
// Distance = d = ((x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2)^½
        return (std::sqrt(((pX - this->x) * (pX - this->x)) + ((pY - this->y) * (pY - this->y)) + ((pZ - this->z) * (pZ - this->z))));
    }

    Ball(float initX, float initY, float initZ, float initRadius) {
        x = initX;
        y = initY;
        z = initZ;
        radius = initRadius;
    }

private:
    float x;
    float y;
    float z;
    float radius;
};

int main(int argc, char** argv) {


    // For testing purposes
/*
    Point p = Point(0.0, 0.0, 0.0);
    Ball b = Ball(0.0, 0.0, 0.0, 1);
    std::cout << "Ball b contains point p? " << boolalpha << b.contains(p) << std::endl;
    std::cout << "Distance: " << b.distance(p) << std::endl;
    std::cout << "Radius: " << b.getRadius() << std::endl;

    Point p2 = Point(0.0, 0.0, -0.0);
    Ball b2 = Ball(0.0, -0.0, 0.0, 0);
    std::cout << "Ball b2 contains point p2? " << boolalpha << b2.contains(p2) << std::endl;
    std::cout << "Distance: " << b2.distance(p2) << std::endl;
    std::cout << "Radius: " << b2.getRadius() << std::endl;

    Point p3 = Point(0.3, 0.5, 0.3);
    Ball b3 = Ball(0.0, 0.0, 0.0, 1);
    std::cout << "Ball b3 contains point p3? " << boolalpha << b3.contains(p3) << std::endl;
    std::cout << "Distance: " << b3.distance(p3) << std::endl;
    std::cout << "Radius: " << b3.getRadius() << std::endl;

    Point p4 = Point(5.0, 7.0, 2.0);
    Ball b4 = Ball(0.0, 2.0, 14.0, 9);
    std::cout << "Ball b4 contains point p4? " << boolalpha << b4.contains(p4) << std::endl;
    std::cout << "Distance: " << b4.distance(p4) << std::endl;
    std::cout << "Radius: " << b4.getRadius() << std::endl;
*/
    // Testing Completed.
///*
    FILE *textPtr = fopen(argv[1], "r");
    // std::cout << argv[1] << std::endl; // Print it to examine contents of the text file.

    int numInputs = 0;
    float tempX = 0;
    float tempY = 0;
    float tempZ = 0;
    float tempRadius = 0;

    float maxX;
    float maxY;
    float maxZ;
    float minX;
    float minY;
    float minZ;

    fscanf(textPtr, "        %d\n", &numInputs); // Parsing first line of text file

    cout << numInputs << std::endl; // For testing purposes

    Ball **array = new Ball *[numInputs]; // Initialize the array of spheres

    fscanf(textPtr, "  %f     %f     %f     %f\n", &tempX, &tempY, &tempZ, &tempRadius);
    array[0] = new Ball(tempX, tempY, tempZ, tempRadius);
    cout << "Ball " << 0 << " radius: " << array[0]->getRadius() << endl;
    cout << "Ball " << 0 << " X, Y, Z: " << array[0]->getX() << ", " << array[0]->getY() << ", " << array[0]->getZ()
         << endl;

    // After reading in the first of the objects, init. the min/max so we can begin to determine the bounds of the box "R" in the first pass.
    maxX = array[0]->maxX();
    maxY = array[0]->maxY();
    maxZ = array[0]->maxZ();
    minX = array[0]->minX();
    minY = array[0]->minY();
    minZ = array[0]->minZ();

    for (int i = 1; i < numInputs; i += 1) {
        fscanf(textPtr, "  %f     %f     %f     %f\n", &tempX, &tempY, &tempZ, &tempRadius);
        array[i] = new Ball(tempX, tempY, tempZ, tempRadius);
        cout << "Ball " << i << " radius: " << array[i]->getRadius() << endl;
        cout << "Ball " << i << " X, Y, Z: " << array[i]->getX() << ", " << array[i]->getY() << ", " << array[i]->getZ()
             << endl;

        // Values have been read into the new object. Now check them against current min/max and update accordingly

        if (array[i]->maxX() > maxX) {
            maxX = (array[i]->maxX());
            cout << "Global Max X updated: " << maxX << endl;
        }
        if (array[i]->maxY() > maxY) {
            maxY = (array[i]->maxY());
            cout << "Global Max Y updated: " << maxY << endl;
        }
        if (array[i]->maxZ() > maxZ) {
            maxZ = (array[i]->maxZ());
            cout << "Global Max Z updated: " << maxZ << endl;
        }

        // And now the min X/Y/Z

        if (array[i]->minX() < minX) {
            minX = (array[i]->minX());
            cout << "Global Min X updated: " << minX << endl;
        }
        if (array[i]->minY() < minY) {
            minY = (array[i]->minY());
            cout << "Global Min Y updated: " << minY << endl;
        }
        if (array[i]->minZ() < minZ) {
            minZ = (array[i]->minZ());
            cout << "Global Min Z updated: " << minZ << endl;
        }
    }

    cout << "Bounds of box: " << endl;
    cout << "X: " << minX << " to " << maxX << endl;
    cout << "Y: " << minY << " to " << maxY << endl;
    cout << "Z: " << minZ << " to " << maxZ << endl;

    float xLength = maxX - minX;
    float yLength = maxY - minY;
    float zLength = maxZ - minZ;

    float volumeOfBox = xLength * yLength * zLength; // Volume of box can now be calculated.

    // Seed random number generation.
    srand (static_cast <unsigned> (time(0)));

    // Begin random number generation.



    float randX = minX + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(maxX-minX)));
    float randY = minY + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(maxY-minY)));
    float randZ = minZ + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(maxZ-minZ)));

    cout << "Testing Random Number generation:" << endl;
    cout << "Random X value in range: " << randX << endl;
    cout << "Random Y value in range: " << randY << endl;
    cout << "Random Z value in range: " << randZ << endl;

    randX = minX + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(maxX-minX)));
    randY = minY + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(maxY-minY)));
    randZ = minZ + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(maxZ-minZ)));

    cout << "Testing Random Number generation:" << endl;
    cout << "Random X value in range: " << randX << endl;
    cout << "Random Y value in range: " << randY << endl;
    cout << "Random Z value in range: " << randZ << endl;



    int pointsToTest = stoi(argv[2]); // Probably will become a command line argument.

    Point **pointStorage = new Point *[pointsToTest];
    int Sum = 0;

for (int i = 0; i < pointsToTest; i++) {
    randX = minX + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(maxX-minX)));
    randY = minY + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(maxY-minY)));
    randZ = minZ + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(maxZ-minZ)));
    pointStorage[i] = new Point(randX, randY, randZ); // Initializing set of points to test
}


    float volumeOfProtein;
    float proportionInside;
    float proportionInsideSquared;
    float standardDeviation;

    float *volStorage = new float [pointsToTest];
    float *StdDevStorage = new float [pointsToTest];


    for (int i = 0; i < pointsToTest; i++){
        for (int j = 0; j < numInputs; j++){
            if (array[j]->contains(*pointStorage[i])){
                Sum += 1;
                cout << "Point located inside protein structure. Total: " << Sum << " / " << i+1 << endl;



                volumeOfProtein = volumeOfBox * Sum / (i+1);
                proportionInside = volumeOfProtein / volumeOfBox;
                proportionInsideSquared = proportionInside * proportionInside;
                standardDeviation = volumeOfBox * std::sqrt((proportionInside - proportionInsideSquared)/(i+1));

                cout << "Volume of protein calculated so far: " << volumeOfProtein << endl;
                cout << "Standard Deviation: " << standardDeviation << endl;

                volStorage[i] = volumeOfProtein;
                StdDevStorage[i] = standardDeviation;

                break;
            }
            if (j == numInputs - 1){
                Sum += 0;
                cout << "Point NOT located inside protein structure. Total: " << Sum << " / " << i+1 << endl;

                volumeOfProtein = volumeOfBox * Sum / (i+1);
                proportionInside = volumeOfProtein / volumeOfBox;
                proportionInsideSquared = proportionInside * proportionInside;
                standardDeviation = volumeOfBox * std::sqrt((proportionInside - proportionInsideSquared)/(i+1));

                cout << "Volume of protein calculated so far: " << volumeOfProtein << endl;
                cout << "Standard Deviation: " << standardDeviation << endl;

                volStorage[i] = volumeOfProtein;
                StdDevStorage[i] = standardDeviation;

                break;
            }
        }

    }

    volumeOfProtein = volumeOfBox * Sum / pointsToTest;
    proportionInside = volumeOfProtein / volumeOfBox;
    proportionInsideSquared = proportionInside * proportionInside;
    standardDeviation = volumeOfBox * std::sqrt((proportionInside - proportionInsideSquared)/pointsToTest);

    cout << endl;
    cout << "Completed testing. Total: " << Sum << " / " << pointsToTest << endl;
    cout << "Volume of protein: " << volumeOfProtein << endl;
    cout << "Volume of box: " << volumeOfBox << endl;
    cout << "Standard Deviation: " << standardDeviation << endl;

    float error = (100 - (100 * volumeOfProtein / 35490.34));
    float absError = abs(100 - (100 * volumeOfProtein / 35490.34));
    cout << "Error %: " << error << endl;


    // At this point, I am printing values for graphing purposes:

    cout  << 2 << " " << volStorage[1] << " " << StdDevStorage[1] << endl;
    cout  << 10 << " " << volStorage[9] << " " << StdDevStorage[9] << endl;
    cout  << 50 << " " << volStorage[49] << " " << StdDevStorage[49] << endl;
    cout  << 100 << " " << volStorage[99] << " " << StdDevStorage[99] << endl;
    cout  << 150 << " " << volStorage[149] << " " << StdDevStorage[149] << endl;
    for (int i = 250; i < pointsToTest; i += (pointsToTest/200)){
        cout  << i << " " << volStorage[i] << " " << StdDevStorage[i] << endl;
        //cout << "Points tested / Volume / StdDev: " << i << " " << volStorage[i] << " " << StdDevStorage[i] << endl;
    }



    // And now separating them to paste into excel

    cout  << 2 << endl;
    cout  << 10 << endl;
    cout  << 50 << endl;
    cout  << 100 << endl;
    cout  << 150 << endl;
    for (int i = 250; i < pointsToTest; i += (pointsToTest/200)){
        cout  << i << endl;
        //cout << "Points tested / Volume / StdDev: " << i << " " << volStorage[i] << " " << StdDevStorage[i] << endl;
    }

    cout   << " " << volStorage[1]  << endl;
    cout   << " " << volStorage[9] << endl;
    cout   << " " << volStorage[49] << endl;
    cout   << " " << volStorage[99] << endl;
    cout   << " " << volStorage[149] << endl;
    for (int i = 250; i < pointsToTest; i += (pointsToTest/200)){
        cout  << volStorage[i]  << endl;
        //cout << "Points tested / Volume / StdDev: " << i << " " << volStorage[i] << " " << StdDevStorage[i] << endl;
    }

    cout  << StdDevStorage[1] << endl;
    cout  <<  StdDevStorage[9] << endl;
    cout  <<  StdDevStorage[49] << endl;
    cout  <<  StdDevStorage[99] << endl;
    cout  <<  StdDevStorage[149] << endl;
    for (int i = 250; i < pointsToTest; i += (pointsToTest/200)){
        cout  <<  StdDevStorage[i] << endl;
        //cout << "Points tested / Volume / StdDev: " << i << " " << volStorage[i] << " " << StdDevStorage[i] << endl;
    }
    return 0;


    // For visualization, we want a table with the number of points tested, our guess at volume, standard deviation, and an asymptotic line with
    // Use matlab for visualizations.
}







/*
 * Basic methodology:
 *
 * int sum;
 *
    Point points [10000]; // Randomly generate these
    Ball proteinStructure [numInput]; // Read these in from the file

    for (Point x : points){
        for (Ball z : proteinStructure){
            if z.contains(x){
                sum +=1;
                break;
            }
        }
        sum += 0;
    }

*/








    /*
     * Proposed program structure:
     *
     * Program takes file name as command line arguments
     *
     * File consists of lines of a number indicating the number of amino acids, followed by X,Y,Z coordinates & a radius
     *
     * First step is to read in the number of amino acids that need to be kept track of
     *
     * Following this, initialize an array of "amino acid" objects of the correct size to store these values
     * Parse the file, reading in protein position and radius as we go into the various protein objects we construct
     * While doing this, keep track of min & max X,Y,Z positions to configure our zone. (Account for radius of objects as well whilst doing this)
     *
     * Once all have been read in, we can begin the problem solving.
     * Note min&max X,Y,Z values such that we can designate a rectangular prism that bounds the volume of the protein.
     *
     * Initialize a set of points at randomly distributed locations within the rectangular prism, which we know
     * the volume "V" to be (maxX-minX) * (maxY-minY) * (maxZ-minZ)
     *
     * For each point, determine if it is located inside the protein or not, by checking to see if the distance from the point generated
     * to the center of any of the "amino acid" balls is less than the radius of aforementioned ball.
     * Each time this is done, update a sum pointsContained which is equal to the number of points thus far contained in the circle, and a sum
     * pointsTested which is the number of points tested so far.
     *
     * Estimated volume can be calculated on the fly as pointsContained/pointsTested * V.
     * Visualization can be updated as this is done.
     *
     * Need to ask professor about displaying standard deviation as the prompt is rather mystifying.
     *
     * Suggest starting number of points should be 10000.
     */


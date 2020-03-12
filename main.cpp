

// Created by Nicholas Treynor 2/24

#include <iostream>
#include <cstdio>      /* printf */
#include <cmath>       /* sqrt */
#include <string>
#include <ctime>
#include <cstdlib>

using namespace std;


float generateRandomNumberInRange(float min, double max){
    float randX = min + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(max-min)));
    return randX;
}


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
        if ((((pX - this->x) * (pX - this->x)) + ((pY - this->y) * (pY - this->y)) + ((pZ - this->z) * (pZ - this->z))) < (this->radius*this->radius)){
            return true;
        } else {
            return false;
        }
    }



    Point generateRandomPointOnSphere(){

        // From Q(u,a) = (R*sqrt(1-(u*u))*cos(a), R*sqrt(1-(u*u))*sin(a), Ru
        // Where R is the radius, a is a randomly generated number from 0 to 2*pi, and u is a randomly generated number from [-1 to 1].

        float a = generateRandomNumberInRange(0.0, 2*M_PI); // Generate a
        float u = generateRandomNumberInRange(-1.0, 1); // Generate u

        // Now we calculate the position of the point:

        float pointX = this->radius*sqrt(1-(u*u))*cos(a);
        float pointY = this->radius*sqrt(1-(u*u))*sin(a);
        float pointZ = this->radius*u;

        // Make sure that the X,Y,Z points are adjusted by the centerpoint of the sphere!

        pointX += getX();
        pointY += getY();
        pointZ += getZ();

        Point NewPoint = Point(pointX, pointY, pointZ);

        return NewPoint;
    }

    float surfaceArea(){
        float R = this->radius;
        float surfaceArea = (4)*(M_PI)*(R*R); // 4 * pi * r^2
        return surfaceArea;
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

    // Command line arguments should be fed in in the following form:
    // String fileName Int PointsToTest
    // Like so:
    // Proteins.txt 50000


    // For testing purposes

    Point p = Point(0.0, 0.0, 0.0);
    Ball b = Ball(0.0, 0.0, 0.0, 1);
    std::cout << "Ball b contains point p? " << boolalpha << b.contains(p) << std::endl;
    std::cout << "Distance: " << b.distance(p) << std::endl;
    std::cout << "Radius: " << b.getRadius() << std::endl;

    Point bP = b.generateRandomPointOnSphere();
    float distance = b.distance(bP);
    cout << "distance from point bP on surface of b to the center of b: " << distance << endl;
    bool isAPointOnTheSurfaceContained = b.contains(bP);
    cout << "Is bP contained inside b? " << boolalpha << isAPointOnTheSurfaceContained << endl;
    cout << "Surface Area of b: ";
    cout << b.surfaceArea() << endl;
 ;/*
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

    Ball **ballStorage = new Ball *[numInputs]; // Initialize the array of spheres

    // Scan in the initial set of values from the first line
    fscanf(textPtr, "  %f     %f     %f     %f\n", &tempX, &tempY, &tempZ, &tempRadius);
    ballStorage[0] = new Ball(tempX, tempY, tempZ, tempRadius);
    cout << "Ball " << 0 << " radius: " << ballStorage[0]->getRadius() << endl;
    cout << "Ball " << 0 << " X, Y, Z: " << ballStorage[0]->getX() << ", " << ballStorage[0]->getY() << ", " << ballStorage[0]->getZ()
         << endl;

    // After reading in the first of the objects, init. the min/max so we can begin to determine the bounds of the box "R" in the first pass.
    maxX = ballStorage[0]->maxX();
    maxY = ballStorage[0]->maxY();
    maxZ = ballStorage[0]->maxZ();
    minX = ballStorage[0]->minX();
    minY = ballStorage[0]->minY();
    minZ = ballStorage[0]->minZ();

    // Now we begin looping through the remainder of the lines in the text file.
    for (int i = 1; i < numInputs; i += 1) {
        fscanf(textPtr, "  %f     %f     %f     %f\n", &tempX, &tempY, &tempZ, &tempRadius);
        ballStorage[i] = new Ball(tempX, tempY, tempZ, tempRadius);
        cout << "Ball " << i << " radius: " << ballStorage[i]->getRadius() << endl;
        cout << "Ball " << i << " X, Y, Z: " << ballStorage[i]->getX() << ", " << ballStorage[i]->getY() << ", " << ballStorage[i]->getZ()
             << endl;

        // Values have been read into the new object. Now check them against current min/max and update accordingly

        if (ballStorage[i]->maxX() > maxX) {
            maxX = (ballStorage[i]->maxX());
            cout << "Global Max X updated: " << maxX << endl;
        }
        if (ballStorage[i]->maxY() > maxY) {
            maxY = (ballStorage[i]->maxY());
            cout << "Global Max Y updated: " << maxY << endl;
        }
        if (ballStorage[i]->maxZ() > maxZ) {
            maxZ = (ballStorage[i]->maxZ());
            cout << "Global Max Z updated: " << maxZ << endl;
        }

        // And now the min X/Y/Z

        if (ballStorage[i]->minX() < minX) {
            minX = (ballStorage[i]->minX());
            cout << "Global Min X updated: " << minX << endl;
        }
        if (ballStorage[i]->minY() < minY) {
            minY = (ballStorage[i]->minY());
            cout << "Global Min Y updated: " << minY << endl;
        }
        if (ballStorage[i]->minZ() < minZ) {
            minZ = (ballStorage[i]->minZ());
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



    int pointsToTest = stoi(argv[2]); // Reading in the number of points to test from the input

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
            if (ballStorage[j]->contains(*pointStorage[i])){
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
    // float absError = abs(100 - (100 * volumeOfProtein / 35490.34));
    cout << "Error %: " << error << endl;

/*
    // At this point, I am printing values for graphing purposes:

    cout  << 2 << " " << volStorage[1] << " " << StdDevStorage[1] << endl;
    cout  << 10 << " " << volStorage[9] << " " << StdDevStorage[9] << endl;
    cout  << 50 << " " << volStorage[49] << " " << StdDevStorage[49] << endl;
    cout  << 100 << " " << volStorage[99] << " " << StdDevStorage[99] << endl;
    cout  << 150 << " " << volStorage[149] << " " << StdDevStorage[149] << endl;
    for (int i = 250; i < pointsToTest; i += (pointsToTest/200)){
        cout  << i << " " << volStorage[i-1] << " " << StdDevStorage[i-1] << endl;
        //cout << "Points tested / Volume / StdDev: " << i << " " << volStorage[i] << " " << StdDevStorage[i] << endl;
    }



    // And now separating them to paste into excel

    cout  << 2 << endl;
    cout  << 10 << endl;
    cout  << 50 << endl;
    cout  << 100 << endl;
    cout  << 150 << endl;
    for (int i = 250; i <= pointsToTest; i += (pointsToTest/200)){
        cout  << i << endl;
        //cout << "Points tested / Volume / StdDev: " << i << " " << volStorage[i] << " " << StdDevStorage[i] << endl;
    }

    cout   << " " << volStorage[1]  << endl;
    cout   << " " << volStorage[9] << endl;
    cout   << " " << volStorage[49] << endl;
    cout   << " " << volStorage[99] << endl;
    cout   << " " << volStorage[149] << endl;
    for (int i = 250; i <= pointsToTest; i += (pointsToTest/200)){
        cout  << volStorage[i-1]  << endl;
        //cout << "Points tested / Volume / StdDev: " << i << " " << volStorage[i] << " " << StdDevStorage[i] << endl;
    }

    cout  << StdDevStorage[1] << endl;
    cout  <<  StdDevStorage[9] << endl;
    cout  <<  StdDevStorage[49] << endl;
    cout  <<  StdDevStorage[99] << endl;
    cout  <<  StdDevStorage[149] << endl;
    for (int i = 250; i <= pointsToTest; i += (pointsToTest/200)){
        cout  <<  StdDevStorage[i-1] << endl;
        //cout << "Points tested / Volume / StdDev: " << i << " " << volStorage[i] << " " << StdDevStorage[i] << endl;
    }
    */





    // At this point, volume has been calculated. Now for surface area.

    /*
     * General strategy:
     *
     * Given a requested number of points to generate on the surface, we will calculate the total surface area of all spheres,
     * then calculate the required density of points on the surface to reach the desired amount of total test points. Given this,
     * we can then begin distributing points across the surfaces of the spheres, and calculate surface area in a manner very similar to that in our methodology for volume.
     */



    int pointsToTestOnSurface = stoi(argv[3]);
    Point **surfacePointStorage = new Point *[pointsToTestOnSurface];



    float totalSurfaceArea = 0;
    for (int i = 0; i < numInputs; i++){
        totalSurfaceArea += ballStorage[i]->surfaceArea();  // Sum total surface area
    }

    float surfaceDensityOfPoints = totalSurfaceArea / pointsToTestOnSurface;



    int pointsGeneratedSoFar = 0;
    for (int i = 0; i < numInputs; i++){ // For each ball we have in our model...
        int remainingPointsToGenerateForThisBall = ((ballStorage[i]->surfaceArea()) / (surfaceDensityOfPoints)); // Figure out how many points need to be generated
        while (remainingPointsToGenerateForThisBall > 0){                                                        // Then create them, and store them in an array.
            remainingPointsToGenerateForThisBall--;
            Point tempPoint = ballStorage[i]->generateRandomPointOnSphere();
            surfacePointStorage[pointsGeneratedSoFar] = new Point(tempPoint.getX(), tempPoint.getY(), tempPoint.getZ());
            pointsGeneratedSoFar++;
        }
    }



    // Points have now been generated. Now we begin the difficult part, testing these points.

    int pointsToTestSurfaceArea = pointsGeneratedSoFar;

    float surfaceAreaOfProtein;
    float proportionOnSurface;
    float proportionOnSurfaceSquared;
    float standardDeviationSurfaceArea;

    float *surfaceAreaStorage = new float [pointsToTestSurfaceArea];
    float *StdDevStorageSurfaceArea = new float [pointsToTestSurfaceArea];

    int surfaceSum = 0;


    for (int i = 0; i < pointsToTestSurfaceArea; i++){                          // For each point we generated on the surface...
        for (int j = 0; j < numInputs; j++){                                    // Loop across all balls we have read in...
            // cout << "Distance between ball and point: " << ballStorage[j]->distance(*surfacePointStorage[i]) << endl;
            // cout << "Radius of ball: " << ballStorage[j]->getRadius() << endl;
            if (ballStorage[j]->contains(*surfacePointStorage[i])){             // If a point is contained inside (IE, not on the surface of the protein) ...
                surfaceSum += 0;                                                // Do not increment the counter of how many are contained inside. Break from this loop and move to the next point, but update data.
                cout << "Point located inside protein structure. Total: " << surfaceSum << " / " << i+1 << endl;



                surfaceAreaOfProtein = totalSurfaceArea * surfaceSum / (i+1);
                proportionOnSurface = surfaceAreaOfProtein / totalSurfaceArea;
                proportionOnSurfaceSquared = proportionOnSurface * proportionOnSurface;
                standardDeviationSurfaceArea = totalSurfaceArea * std::sqrt((proportionOnSurface - proportionOnSurfaceSquared)/(i+1));

                cout << "Surface area of protein calculated so far: " << surfaceAreaOfProtein << endl;
                cout << "Standard Deviation Of Surface Area: " << standardDeviationSurfaceArea << endl;

                surfaceAreaStorage[i] = surfaceAreaOfProtein;
                StdDevStorageSurfaceArea[i] = standardDeviationSurfaceArea;

                break;
            }
            if (j == numInputs - 1){ // If all points have been tested and none have been contained...
                surfaceSum += 1;                                                    // If NOT contained, we know this point must lay on the surface, and so we increment the count
                cout << "Point located ON THE SURFACE of the protein structure. Total: " << surfaceSum << " / " << i+1 << endl;

                surfaceAreaOfProtein = totalSurfaceArea * surfaceSum / (i+1);
                proportionOnSurface = surfaceAreaOfProtein / totalSurfaceArea;
                proportionOnSurfaceSquared = proportionOnSurface * proportionOnSurface;
                standardDeviationSurfaceArea = totalSurfaceArea * std::sqrt((proportionOnSurface - proportionOnSurfaceSquared)/(i+1));

                cout << "Surface area of protein calculated so far: " << surfaceAreaOfProtein << endl;
                cout << "Standard Deviation Of Surface Area: " << standardDeviationSurfaceArea << endl;

                surfaceAreaStorage[i] = surfaceAreaOfProtein;
                StdDevStorageSurfaceArea[i] = standardDeviationSurfaceArea;

                break;
            }
        }

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


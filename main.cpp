

// Created by Nicholas Treynor 2/24

#include <iostream>
#include <cstdio>      /* printf */
#include <cmath>       /* sqrt */
#include <string>
#include <ctime>
#include <cstdlib>

using namespace std;


double generateRandomNumberInRange(double min, double max){
    double randX = min + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(max-min)));
    return randX;
}


class Point { // An object class that we will use to represent the points we will generate within our union of balls model of a protein structure.
public:
    double getX() {
        return x;
    }

    double getY() {
        return y;
    }

    double getZ() {
        return z;
    }

    Point(double initX, double initY, double initZ) { // Constructor for the object
        x = initX;
        y = initY;
        z = initZ;
    }

private:
    double x;
    double y;
    double z;




};

class Ball { // An object class that we will use to represent balls within our union of balls model of a protein.
public:
    double getX() {
        return x;
    }

    double getY() {
        return y;
    }

    double getZ() {
        return z;
    }


    // The following min/max x/y/z functions modify the center of the point by the radius in various directions in order to determine the true maximum and minimum X,Y, and Z coordinates of the sphere.
    // This will later be utilized to determine the global max and min X,Y, and Z points that will bound our box "R".
    double maxX() {
        return x + radius;
    }

    double maxY() {
        return y + radius;
    }

    double maxZ() {
        return z + radius;
    }

    double minX() {
        return x - radius;
    }

    double minY() {
        return y - radius;
    }

    double minZ() {
        return z - radius;
    }

    double getRadius() {
        return radius;
    }

    bool contains(Point point) {
        double pX = point.getX();
        double pY = point.getY();
        double pZ = point.getZ();
// Distance = d = ((x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2)^½
// We will use a slight modification to my first method, which was comparing distance to the radius of the sphere, and instead compare distance squared against radius squared, as it is easier to compute,
// as we do not need to utilize the Sqrt() function.
        if ((((pX - this->x) * (pX - this->x)) + ((pY - this->y) * (pY - this->y)) + ((pZ - this->z) * (pZ - this->z))) < (this->radius*this->radius)){
            return true;
        } else {
            return false;
        }
    }

    Point generateRandomPointOnSphere(){

        // This is Cory Simon's method for generating points uniformly over a unit circle, modified slightly to generate points just beyond the surface of the sphere this function is called on.
        // Writeup can be found here: http://corysimon.github.io/articles/uniformdistn-on-sphere/

        double theta = 2 * M_PI * generateRandomNumberInRange(0.0, 1);
        double phi = acos(1 - 2 * generateRandomNumberInRange(0.0, 1));

        // Calculate the position of the point on surface RELATIVE to the center point:

        double x = sin(phi) * cos(theta) * getRadius() * 1.000000000001; // Multiplying by 1.000000000001 allows us to ensure rounding errors do not cause the point to be counted as "contained" within the ball this point
        double y = sin(phi) * sin(theta) * getRadius() * 1.000000000001; // is generated on the surface of, without impacting in any significant way the accuracy of our program.
        double z = cos(phi) * getRadius() * 1.000000000001;




        double pointX = x;
        double pointY = y;
        double pointZ = z;

        // Now we make sure that the X,Y,Z points are adjusted by the center point of the sphere!

        pointX += getX();
        pointY += getY();
        pointZ += getZ();

        Point NewPoint = Point(pointX, pointY, pointZ);

        return NewPoint;
    }

    double surfaceArea(){
        double R = this->radius;
        double surfaceArea = (4)*(M_PI)*(R*R); // 4 * pi * r^2
        return surfaceArea;
    }

    double distance(Point point) {
        double pX = point.getX();
        double pY = point.getY();
        double pZ = point.getZ();
     // Distance formula: d = ((x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2)^½
        return (std::sqrt(((pX - this->x) * (pX - this->x)) + ((pY - this->y) * (pY - this->y)) + ((pZ - this->z) * (pZ - this->z))));
    }

    Ball(double initX, double initY, double initZ, double initRadius) { // Constructor for the object
        x = initX;
        y = initY;
        z = initZ;
        radius = initRadius;
    }

private:
    double x;
    double y;
    double z;
    double radius;
};



int main(int argc, char** argv) {

    // Command line arguments should be fed in in the following form:
    // String fileName Int PointsToTest Int PointsToTestOnSurface
    // Like so:
    // Proteins.txt 5000 10000

    // The above example would test 5000 points in an attempt to determine volume, and 10,000 points distributed across the surface of balls in order to determine surface area.

    FILE *textPtr = fopen(argv[1], "r"); // Open up our file from input to parse.

    int pointsToTest = stoi(argv[2]); // Reading in the number of points to test from the input (for volume)
    int pointsToTestOnSurface = stoi(argv[3]); // Reading in the number of points to test from the input (for surface area)

    // Here we initialize our containers, to be updated shortly.

    int numInputs = 0;
    double tempX = 0;
    double tempY = 0;
    double tempZ = 0;
    double tempRadius = 0;

    double maxX;
    double maxY;
    double maxZ;
    double minX;
    double minY;
    double minZ;

    fscanf(textPtr, "        %d\n", &numInputs); // Parsing first line of text file

    Ball **ballStorage = new Ball *[numInputs]; // Initialize the array of spheres using our command line inputs.

    // Scan in the initial set of values from the second line
    fscanf(textPtr, "  %lf     %lf     %lf     %lf\n", &tempX, &tempY, &tempZ, &tempRadius);
    ballStorage[0] = new Ball(tempX, tempY, tempZ, tempRadius);

    // After reading in the first of the objects, init. the min/max so we can begin to determine the bounds of the box "R" in the first pass.
    maxX = ballStorage[0]->maxX();
    maxY = ballStorage[0]->maxY();
    maxZ = ballStorage[0]->maxZ();
    minX = ballStorage[0]->minX();
    minY = ballStorage[0]->minY();
    minZ = ballStorage[0]->minZ();

    // Now we begin looping through the remainder of the lines in the text file.
    for (int i = 1; i < numInputs; i += 1) {
        fscanf(textPtr, "  %lf     %lf     %lf     %lf\n", &tempX, &tempY, &tempZ, &tempRadius);
        ballStorage[i] = new Ball(tempX, tempY, tempZ, tempRadius);

        // Values have been read into the new object. Now check them against current min/max and update accordingly
        if (ballStorage[i]->maxX() > maxX) {
            maxX = (ballStorage[i]->maxX()); // Global minimum X has been updated
        }
        if (ballStorage[i]->maxY() > maxY) {
            maxY = (ballStorage[i]->maxY()); // Global minimum Y has been updated
        }
        if (ballStorage[i]->maxZ() > maxZ) {
            maxZ = (ballStorage[i]->maxZ()); // Global minimum Z has been updated
        }

        // And now the min X/Y/Z

        if (ballStorage[i]->minX() < minX) {
            minX = (ballStorage[i]->minX()); // Global minimum X has been updated
        }
        if (ballStorage[i]->minY() < minY) {
            minY = (ballStorage[i]->minY()); // Global minimum Y has been updated
        }
        if (ballStorage[i]->minZ() < minZ) {
            minZ = (ballStorage[i]->minZ()); // Global minimum Z has been updated
        }
    }

// Now that we have the global min and max, we can determine the dimensions of our box R.
    double xLength = maxX - minX;
    double yLength = maxY - minY;
    double zLength = maxZ - minZ;

    double volumeOfBox = xLength * yLength * zLength; // Volume of box can now be calculated. This is our box "R".

    Point **pointStorage = new Point *[pointsToTest]; // Generating an array to store our points as we create them.

    // Seed random number generation, and generate containers for random values
    srand (static_cast <unsigned> (time(0)));

    double randX;
    double randY;
    double randZ;

    // Now we will fill our box R with randomly generated points, storing the points in an array
for (int i = 0; i < pointsToTest; i++) {
    randX = generateRandomNumberInRange(minX, maxX);
    randY = generateRandomNumberInRange(minY, maxY);
    randZ = generateRandomNumberInRange(minZ, maxZ);
    pointStorage[i] = new Point(randX, randY, randZ); // Initializing set of points to test
}


    int Sum = 0;                                             // Initialize our counter for the number of points contained within the protein structure.

    for (int i = 0; i < pointsToTest; i++){                  // For each point we generated in the box R....
        for (int j = 0; j < numInputs; j++){                 // loop across all balls we have read in...
            if (ballStorage[j]->contains(*pointStorage[i])){ // if a point is contained inside one of these balls, we...
                Sum += 1;                                    // increment the counter of how many are contained inside. Break from this loop and move to the next point.
                break;
            }
            if (j == numInputs - 1){                         // If the point has been checked against every ball and none have contained the point...
                Sum += 0;                                    //  we know this point must lay outside the protein structure, and so we do not increment the count, instead moving on.
                break;
            }
        }
    }

    // Now we can calculate our volume and standard deviation, using the sum we found earlier.
    // Proportion inside is equivalent to < f^2 >, and proportionInsideSquared is equivalent to < f >^2 in the prompt.
    // It's worth noting that < f^2 > and < f > are identical, as f(x) = f(x)^2 when f(x) can only ever be 1 or 0. (1 = 1 * 1, 0 = 0 * 0)
    double volumeOfProtein = volumeOfBox * Sum / pointsToTest;
    double proportionInside = volumeOfProtein / volumeOfBox;
    double proportionInsideSquared = proportionInside * proportionInside;
    double standardDeviation = volumeOfBox * std::sqrt((proportionInside - proportionInsideSquared)/pointsToTest);

    cout << endl;
    // cout << "Completed testing. Total: " << Sum << " / " << pointsToTest << endl;
    cout << "Volume of protein: " << volumeOfProtein << " cubic angstrom " << endl;
    // cout << "Volume of box: " << volumeOfBox << endl;
    cout << "Standard deviation of volume: " << standardDeviation << " cubic angstrom " << endl;


    // At this point, volume has been calculated. Now for surface area.

    /*
     * General strategy:
     *
     * Given a requested number of points to generate on the surface, we will calculate the total surface area of all spheres,
     * then calculate the required density of points on the surface to reach the desired amount of total test points. Given this,
     * we can then begin distributing points across the surfaces of the spheres, and calculate surface area in a manner very similar to that in our methodology for volume.
     */



    Point **surfacePointStorage = new Point *[pointsToTestOnSurface];

    double totalSurfaceArea = 0;
    for (int i = 0; i < numInputs; i++){
        totalSurfaceArea += ballStorage[i]->surfaceArea();  // Sum total surface area
    }

    double surfaceDensityOfPoints = totalSurfaceArea / pointsToTestOnSurface;    // Calculate the desired density of points on the surface of balls.

    int pointsGeneratedSoFar = 0;
    for (int i = 0; i < numInputs; i++){                                                                         // For each ball we have in our model...
        int remainingPointsToGenerateForThisBall = ((ballStorage[i]->surfaceArea()) / (surfaceDensityOfPoints)); // Figure out how many points need to be generated, utilizing surface area of a ball and desired density of points
        while (remainingPointsToGenerateForThisBall > 0){                                                        // Then create them, and store them in an array.
            remainingPointsToGenerateForThisBall--;
            Point tempPoint = ballStorage[i]->generateRandomPointOnSphere();                                     // Utilizing our method for random point generation
            surfacePointStorage[pointsGeneratedSoFar] = new Point(tempPoint.getX(), tempPoint.getY(), tempPoint.getZ());
            pointsGeneratedSoFar++;
        }
    }



    // Points have now been generated. Now we begin the time consuming part, testing these points.

    int pointsToTestSurfaceArea = pointsGeneratedSoFar;


    int surfaceSum = 0; // Initialize the sum we will use to count the number of points not contained within the other balls (That lay on the true surface of the protein)


    for (int i = 0; i < pointsToTestSurfaceArea; i++){                          // For each point we generated on the surface...
        for (int j = 0; j < numInputs; j++){                                    // Loop across all balls we have read in...

            if (ballStorage[j]->contains(*surfacePointStorage[i])){             // If a point is contained inside (IE, not on the surface of the protein) ...
                surfaceSum += 0;                                                // Do not increment the counter of how many are contained inside. Break from this loop and move to the next point, but update data.
                break;
            }
            if (j == numInputs - 1){                                            // If all points have been tested and none have been contained...
                surfaceSum += 1;                                                //  we know this point must lay on the surface, and so we increment the count of points that lay on the true surface of the protein
                break;
            }
        }
    }


    // Now we perform our final calculations and output our results after.
    // ProportionOnSurface is equivalent to < f^2 >, and proportionOnSurfaceSquared is equivalent to < f >^2 in the prompt.
    // It's worth noting that < f^2 > and < f > are identical, as f(x) = f(x)^2 when f(x) can only ever be 1 or 0. (1 = 1 * 1, 0 = 0 * 0)

    double surfaceAreaOfProtein = totalSurfaceArea * surfaceSum / (pointsToTestSurfaceArea);
    double proportionOnSurface = surfaceAreaOfProtein / totalSurfaceArea;
    double proportionOnSurfaceSquared = proportionOnSurface * proportionOnSurface;
    double standardDeviationSurfaceArea = totalSurfaceArea * std::sqrt((proportionOnSurface - proportionOnSurfaceSquared)/(pointsToTestSurfaceArea));

    cout << "Surface area of protein: " << surfaceAreaOfProtein << " square angstrom " << endl;
    cout << "Standard deviation of surface area: " << standardDeviationSurfaceArea << " square angstrom " << endl;

    return 0; // Concluded the program. End.

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

    After this has been done, we can calculate the proportion of points that fall inside the union of balls from "sum" and
    "numInput", and multiply that times the bounds of the box we generated the points within

*/








    /*
     * Proposed program structure:
     *
     * Program takes file name as command line arguments, in addition to a number of desired points to test
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


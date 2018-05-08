#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <limits>

#define DIMENSION 2
#define INPUTFILE "dbscan_input.txt"

enum state {visited, unvisited};
enum group {inCluster, noise, undecided};

struct point {
    double *coord;
    bool isCorePoint;
    state st;
    group gr;
};

struct edge;

struct cell {
    std::vector<edge*> edges;
    std::vector<point*> points;
    double startCoord[DIMENSION];
    bool hasEnoughPoints;
    bool coreCell;
    int pointCounter = 0;
    int id;                 //union find
    int parent;
};

struct edge {
    cell *firstCell, *secondCell;
};

struct level {
    double levelGridSize;
    std::vector<cell*> cells;
};

struct cluster {
    std::vector<cell*> cells;
};

//TODO upravit na citanie bodov z buniek
void printActualCluster(cluster &actualCluster) {
    for(int i=0; i<(int)actualCluster.cells.size() ; i++) {
        for(int j=0 ; j<actualCluster.cells[i]->points.size() ; j++){
            for(int k=0 ; k<DIMENSION ; k++)
                std::cout << actualCluster.cells[i]->points[j]->coord[k] << " ";

            if(actualCluster.cells[i]->points[j]->st == 0)
                std::cout << "visited, ";
            else
                std::cout << "unvisited, ";

            if(actualCluster.cells[i]->points[j]->gr == 0)
                std::cout << "inCluster";
            else if(actualCluster.cells[i]->points[j]->gr == 1)
                std::cout << "noise";
            else
                std::cout << "undecided";
            std::cout << "\n";
        }
    }
    getchar();
}

void printAllClusters(std::vector<cluster*> clusters) {
    for(int i=0 ; i<clusters.size() ; i++) {
        for(int j=0 ; j<clusters[i]->cells.size() ; j++){
            int k = clusters[i]->cells[j]->points.size();
            std::cout << i+1 << ". zhluk:";
            for(int l=0 ; l<k ; l++) {
                std::cout << "\n\t";
                for(int m=0 ; m<DIMENSION ; m++)
                std::cout << clusters[i]->cells[j]->points[l]->coord[m] << " ";
            }
            std::cout << "\n";
        }
    }
}

void printNoise(std::vector<point> &input) {

    int count = 0;

    std::cout << "\nNoise:\n";
    for(int i=0 ; i<(int)input.size() ; i++) {
        if(input[i].gr == noise) {
            count++;
            for(int j=0 ; j<DIMENSION ; j++)
                std::cout << input[i].coord[j] << " ";
            std::cout << "\n";
        }
    }
    if(count == 0)
        std::cout << "no noise\n";
}

void getInput(std::vector<point> &input, double *bordersMax, double *bordersMin){

    int index = 0;
    std::ifstream file;

    file.open(INPUTFILE);
    for(int i = 0 ; i<DIMENSION ; i++){
        bordersMax[i] = 0;
        bordersMin[i] = std::numeric_limits<double>::max();
    }

    if (file.is_open()) {
        while (!file.eof()) {
            input.push_back(point());
            input[index].coord = new double[DIMENSION];

            for(int i=0 ; i<DIMENSION ; i++) {
                file >> input[index].coord[i];
                if(input[index].coord[i] > bordersMax[i])
                    bordersMax[i] = input[index].coord[i];
                if(input[index].coord[i] < bordersMin[i])
                    bordersMin[i] = input[index].coord[i];
            }
            input[index].st = unvisited;
            input[index].gr = undecided;
            input[index].isCorePoint = false;

            index++;
        }
    }
    file.close();
}

void tagNeighbours(point actual, std::vector<point> &input, std::vector<point*> &neighbours, double maxEps){

    double dist;
    int i, j;

    for(i=0 ; i<(int)input.size() ; i++) {
        dist = 0;

        for (j=0 ; j<DIMENSION ; j++)
            dist += pow(actual.coord[j] - input[i].coord[j], 2);

        dist = sqrt(dist);

        if(dist < maxEps) {
            neighbours.push_back(&input[i]);
        }
    }
}
/*
cluster expand(std::vector<point*> &neighbours, cluster &actualCluster, double maxEps, int minPoints, std::vector<point> &input){

    int sum;
    int neighNum = neighbours.size();
    std::vector<point*> expandNeighbours;

    for(int i=0 ; i<neighNum ; i++){
        if(neighbours[i]->st == unvisited) {
            neighbours[i]->st = visited;
            neighbours[i]->gr = inCluster;
            actualCluster.points.push_back(neighbours[i]);

            expandNeighbours.clear();
            tagNeighbours(*neighbours[i], input, expandNeighbours, maxEps);
            sum = expandNeighbours.size();

            if(sum >= minPoints) {
                neighbours.insert(neighbours.end(), expandNeighbours.begin(), expandNeighbours.end());
                neighNum = neighbours.size();
            }
        }
    }

    return actualCluster;
}
*/
void recursiveCellInit(int actualDim, level *levelToInit, double gridSize, double *bordersMax, double *bordersMin, double *otherAxis){

    cell *initCell;

    if(actualDim == 0)
        return;

    if(actualDim == 1){
        for(double i = bordersMin[0] ; i<bordersMax[0] + gridSize ; i+=gridSize){
            initCell = new cell;
            initCell->hasEnoughPoints = false;
            initCell->coreCell = false;
            initCell->startCoord[0] = i;
            for(int j = DIMENSION-1 ; j>0 ; j--)
                initCell->startCoord[j] = otherAxis[j];
            levelToInit->cells.push_back(initCell);
        }
        return;
    }

    while(otherAxis[actualDim-1] < bordersMax[actualDim-1] + gridSize){
        recursiveCellInit(actualDim-1, levelToInit, gridSize, bordersMax, bordersMin, otherAxis);
        otherAxis[actualDim-1]+= gridSize;
    }
    otherAxis[actualDim-1] = bordersMin[actualDim-1];
}

void initLevel(level *levelToInit, double gridSize, double *bordersMax, double *bordersMin){

    //double rho = 0;
    double axis[DIMENSION];

    for(int i = 0 ; i<DIMENSION ; i++)
        axis[i] = bordersMin[i];

    //chcem n-rozmernu kocku (zlozenu z buniek prvej urovne)
    //nad celym priestorom v ktorom sa vyskytuju body
    for(double i = bordersMin[DIMENSION-1] ; i<bordersMax[DIMENSION-1] + gridSize ; i+=gridSize){
        axis[DIMENSION-1] = i;
        recursiveCellInit(DIMENSION-1, levelToInit, gridSize, bordersMax, bordersMin, axis);
    }
}

void assignPoints(std::vector<point> &input, level *specificLevel, double gridSize){

    for(int i=0 ; i<input.size() ; i++)
        for(int j = 0 ; j<specificLevel->cells.size() ; j++)
            if(input[i].coord[0] > specificLevel->cells[j]->startCoord[0]
                && input[i].coord[0] <= specificLevel->cells[j]->startCoord[0] + gridSize
                && input[i].coord[1] > specificLevel->cells[j]->startCoord[1]
                && input[i].coord[1] <= specificLevel->cells[j]->startCoord[1] + gridSize){
                specificLevel->cells[j]->pointCounter++;
                specificLevel->cells[j]->points.push_back(&input[i]);
                }
            else
                continue;

}

level *evaluateCells(level *specificLevel, int minPts){

    int counter = specificLevel->cells.size();
    level *nonEmpty = new level;
    int idCounter = 0;
    nonEmpty->levelGridSize = specificLevel->levelGridSize;

    for(int i = 0 ; i<counter ; i++)
        if(specificLevel->cells[i]->pointCounter > 0){
            specificLevel->cells[i]->hasEnoughPoints = true;
            specificLevel->cells[i]->id = idCounter;
            specificLevel->cells[i]->parent = idCounter++;
            nonEmpty->cells.push_back(specificLevel->cells[i]);
        }
    //free(specificLevel->cells.clear());
    //free(specificLevel);

    return nonEmpty;
}

void checkCoreClusters(level *specificLevel, int minPts, double maxEps, std::vector<point> &input){

    int counter = specificLevel->cells.size();
    int sum;
    std::vector<point*> neighbours;

    for(int i = 0 ; i<counter ; i++){
        if(specificLevel->cells[i]->pointCounter >= minPts)
            specificLevel->cells[i]->coreCell = true;
        else{
            for(int j = 0 ; j<specificLevel->cells[i]->pointCounter ; j++){
                tagNeighbours(*(specificLevel->cells[i]->points[j]), input, neighbours, maxEps);
                sum = neighbours.size();
                if(sum >= minPts){
                    specificLevel->cells[i]->points[j]->isCorePoint = true;
                    specificLevel->cells[i]->coreCell = true;
                }
                neighbours.clear();
            }
        }
    }
}

double getPointDistance(point *p1, point *p2) {

    double distance = 0;

    for (int i = 0 ; i<DIMENSION ; i++)
        distance += pow(p1->coord[i] - p2->coord[i], 2);
    distance = sqrt(distance);

    return distance;
}


double getCellDistance(cell *firstCell, cell *secondCell, double maxEps){

    double cellDistance;

    for(int i = 0 ; i<firstCell->pointCounter ; i++)
        for(int j = 0 ; j<secondCell->pointCounter ; j++){
            cellDistance = getPointDistance(firstCell->points[i], secondCell->points[j]);
            if(cellDistance <= maxEps)
                return cellDistance;
        }

    return maxEps + 1;
}

void addEdges(level *specificLevel, double maxEps){

    int counter = specificLevel->cells.size();
    double distance;

    for(int i = 0 ; i<counter ; i++){
        for(int j = i+1 ; j<counter ; j++){
            distance = getCellDistance(specificLevel->cells[i], specificLevel->cells[j], maxEps);
            if(distance < maxEps){
                edge *validEdge;
                validEdge = new edge();
                validEdge->firstCell = specificLevel->cells[i];
                validEdge->secondCell = specificLevel->cells[j];
                specificLevel->cells[i]->edges.push_back(validEdge);
            }
        }
    }
}

void unionCells(edge *clusterConnect){

    int root1, root2;
    root1 = clusterConnect->firstCell->parent;
    root2 = clusterConnect->secondCell->parent;

    if(root1 >= 0 && root2 >= 0)
        clusterConnect->secondCell->parent = root1;
}

void connectClusters(level *specificLevel){

    for(int i = 0 ; i<specificLevel->cells.size() ; i++){
        for( int j = 0 ; j<specificLevel->cells[i]->edges.size() ; j++){
            unionCells(specificLevel->cells[i]->edges[j]);
        }
    }
}

void fillClusters(level *specificLevel, std::vector<cluster> clusters){
    std::vector<int> clusterIndexes[specificLevel->cells.size()];   //TODO premenuj
    int clusterCounter = 0;
    cluster *addCluster;
    for(int i = 0 ; i<specificLevel->cells.size() ; i++){
        clusterIndexes[specificLevel->cells[i]->parent].push_back(i);
    }

    //vypis indexov buniek
    for(int i = 0 ; i<specificLevel->cells.size() ; i++){
        if(clusterIndexes[i].size() > 0){
            for(int j = 0 ; j<clusterIndexes[i].size() ; j++)
                std::cout << clusterIndexes[i][j] << " ";
            std::cout << "\n";
        }
    }


    addCluster = new cluster();
    clusters.push_back(*addCluster);
std::cout << clusters.size() << "\n";
return;
    for(int i = 0 ; i<specificLevel->cells.size() ; i++){
        if(clusterIndexes[i].size() > 0){
            addCluster = new cluster();
            clusters.push_back(*addCluster);
            //clusterCounter++;
            //for(int j = 0 ; j<clusterIndexes[i].size() ; j++){
            //    *clusters[clusterCounter-1].cells.push_back(specificLevel->cells[clusterIndexes[i][j]);
            //}
        }
    }
    std::cout << clusters.size() << "\n";
}

int main(int argc, char* argv[]) {
    std::vector<point> input;
    std::vector<cluster*> clusters;
    std::vector<point*> neighbours;
    std::vector<level*> all_levels;

    int minPts, numOfNeigh, clustIndex = 0;
    double maxEps, rho, gridSize;
    double *coordBordersMax = new double[DIMENSION];
    double *coordBordersMin = new double[DIMENSION];

    getInput(input, coordBordersMax, coordBordersMin);

    std::cout << "Zvolte minPts: ";
    std::cin >> minPts;
    std::cout << "\n" << "Zvolte maxEps: ";
    std::cin >> maxEps;
    //std::cout << "\n" << "Zvolte velkost rho: ";
    //std::cin >> rho;

    gridSize = maxEps / sqrt(2);
    level firstLevel;

    firstLevel.levelGridSize = gridSize;
    all_levels.push_back(&firstLevel);

    initLevel(all_levels[0], gridSize, coordBordersMax, coordBordersMin);

    //for(int i = 0 ; i<all_levels[0]->cells.size() ; i++){
    //    std::cout << all_levels[0]->cells[i]->startCoord[0] << " " << all_levels[0]->cells[i]->startCoord[1] << " " << all_levels[0]->cells[i]->startCoord[2] << "\n";
    //}
    //std::cout << "\n" << all_levels[0]->cells.size() << "\n\n";

    assignPoints(input, all_levels[0], gridSize);
    all_levels[0] = evaluateCells(all_levels[0], minPts);

    checkCoreClusters(all_levels[0], minPts, maxEps, input); //NeCore nehladat tagNeigh pre cely input, len v susedoch do eps
    addEdges(all_levels[0], maxEps); //TODO pozerat len susedne bunky

    std::cout << all_levels[0]->cells.size() << "\n\n";

    //for(int j = 0 ; j<all_levels[0]->cells.size() ; j++)
    //    std::cout << all_levels[0]->cells[j]->parent << "\n";

    connectClusters(all_levels[0]);

    //for(int j = 0 ; j<all_levels[0]->cells.size() ; j++)
    //    std::cout << all_levels[0]->cells[j]->parent << "\n";

    //fillClusters(all_levels[0], clusters);

    //miesto fillClustera
    std::vector<int> clusterIndexes[all_levels[0]->cells.size()];   //TODO premenuj
    int clusterCounter = 0;
    cluster *addCluster;
    for(int i = 0 ; i<all_levels[0]->cells.size() ; i++){
        clusterIndexes[all_levels[0]->cells[i]->parent].push_back(i);
    }

    for(int i = 0 ; i<all_levels[0]->cells.size() ; i++){
        if(clusterIndexes[i].size() > 0){
            for(int j = 0 ; j<clusterIndexes[i].size() ; j++)
                std::cout << clusterIndexes[i][j] << " ";
            std::cout << "\n";
        }
    }


    for(int i = 0 ; i<all_levels[0]->cells.size() ; i++){
        if(clusterIndexes[i].size() > 0){
            addCluster = new cluster();

            clusters.push_back(addCluster);
            clusterCounter++;
            for(int j = 0 ; j<clusterIndexes[i].size() ; j++){
                (*(clusters[clusterCounter-1])).cells.push_back(all_levels[0]->cells[clusterIndexes[i][j]]);
            }
        }
    }
    std::cout << "\n" << clusters.size() << "\n\n";


    //printAllClusters(clusters);
    printNoise(input);
    std::cout << "pocet zhlukov " << clusters.size();

    return 0;
}

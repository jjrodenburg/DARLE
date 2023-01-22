#include <iostream>
#include "main.h"
#include "LabPixel.h"
#include "Node.h"
#include "CLIInputParser.h"
#include <bitset>

#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/opencv.hpp>
#include <SDL.h>
#include <fstream>

using std::cout;
using namespace cv;
using std::vector;
using std::map;
using std::pair;
using std::get;

template<typename T>
std::vector<T> flatten(std::vector<std::vector<T>> const &vec);

void Draw(vector<vector<vector<float>>> matrix);
std::tuple<float, float, float> RGBToLAB(int r, int g, int b);
std::tuple<float, float, float> labToRGB(float l, float a, float b);
std::vector<std::string> DecodeFile(std::string fileName);
std::vector<std::string> DecodeFileBinary(std::string fileName);

vector<vector<vector<float>>> DecodeMatrix(std::vector<std::string> decodedFile);
Mat ReadImage(std::string filepath);
std::vector<std::vector<std::string>> CreateMatrix(Mat image, int userEDelta);
std::string WriteMatrixToFile(std::string originalFilePath, std::vector<std::vector<std::string>> matrix);
std::string WriteMatrixToFile_NoComma_SingleSeparator(std::string originalFilePath, std::vector<std::vector<std::string>> matrix);

int main(int argc, char *argv[]) {
    CLIInputParser cliOptions = CLIInputParser(argc, argv);

    std::string filePath = "";
    int edelta = 1;
    bool compress = false;
    bool render = false;

    if (cliOptions.cmdOptionExists("--f")) { // --f is input for filepath
        filePath = cliOptions.getCmdOption("--f");
    } else {
        cout << "please provide file";
        return 1;
    }

    if (cliOptions.cmdOptionExists("-c")) { // -c is for compress flag
        compress = true;
    }

    if (cliOptions.cmdOptionExists("-r")) { // -r is for render flag
        render = true;
    }

    if (cliOptions.cmdOptionExists("-e")) { // -r is for render flag
        edelta = stoi(cliOptions.getCmdOption("-e"));
    }

    if (compress) {
        Mat image = ReadImage(filePath);
        cout << "prepped file\n";

        std::vector<std::vector<std::string>> matrix = CreateMatrix(image, edelta);
        cout << "created matrix\n";

        if (cliOptions.cmdOptionExists("-x")) {
            filePath = WriteMatrixToFile_NoComma_SingleSeparator(filePath, matrix);
        } else {
            filePath = WriteMatrixToFile(filePath, matrix);
        }
        cout << "wrote matrix to file\n";
    }

    if (render) {
        cout << "loading matrix\n";
        std::vector<std::string> decodedFile;
        if (cliOptions.cmdOptionExists("-x")) {
            decodedFile = DecodeFile(filePath);
        } else {
            decodedFile = DecodeFileBinary(filePath);
        }

        cout << "loaded matrix\n";
        std::vector<std::vector<std::vector<float>>> decodedMatrix = DecodeMatrix(decodedFile);

        cout << "decoded matrix\n";

        Draw(decodedMatrix);
    }

    cout << "done!";
    return 0;
}


Mat ReadImage(std::string filepath) {
    return imread(filepath, IMREAD_UNCHANGED);
}

std::vector<std::vector<std::string>> CreateMatrix(Mat image, int userEDelta) {
    std::vector<std::vector<std::string>> matrix = std::vector<std::vector<std::string>>();

    AreaNode currentArea;
    currentArea.Location = -1;

    std::map<std::string, std::string> dictionary = std::map<std::string, std::string>();

    int location = 0;
    for(int y=-1;y<image.rows;y++)
    {
        std::vector<std::string> columnNodes;

        if (y == -1) {
            columnNodes.push_back(std::to_string(image.size().width) + "," + std::to_string(image.size().height));
            matrix.push_back(columnNodes);
            continue;
        }

        for(int x=0;x<image.cols;x++)
        {
            // get pixel
            Vec3b & color = image.at<Vec3b>(Point(x,y));


            std::tuple<float, float, float> labColors = RGBToLAB(color[2], color[1], color [0]); // todo figure out difference between stack and my implementation
            CIEDE2000::LAB labPixel = CIEDE2000::LAB();
            labPixel.l = get<0>(labColors);
            labPixel.a = get<1>(labColors);
            labPixel.b = get<2>(labColors);

            LabPixel memoryPixel = LabPixel(labPixel, y * x);

            double eDelta = CIEDE2000::CIEDE2000(currentArea.LabEncoding, labPixel);

            // redundant pixel per edelta (indistinguishable by human eye from previous pixel)
            if (eDelta <= userEDelta && currentArea.Location != -1) {
                location++;
                continue;
            }

            AreaNode area = AreaNode(memoryPixel.location, 0, memoryPixel.labEncoding);

            currentArea = area;

            std::string dictKey = std::to_string(labPixel.l) + "," + std::to_string(labPixel.a) + "," + std::to_string(labPixel.b);
            if (dictionary.contains(dictKey)) {
                columnNodes.push_back(std::to_string(location) + ',' + dictionary.at(dictKey));
            } else {
                columnNodes.push_back(std::to_string(location) + ',' + dictKey);
                dictionary[dictKey] = std::to_string(dictionary.size());
            }
            location++;
        }

        if (!columnNodes.empty()) matrix.push_back(columnNodes);
    }

    std::vector<std::string> dictionaryNodes;
    dictionaryNodes.push_back("*");
    for (auto const& [key, val] : dictionary)
    {
        dictionaryNodes.push_back(key + '#');
    }
    dictionaryNodes.push_back("*");

    matrix.insert(matrix.begin()+1, 1, dictionaryNodes);
    return matrix;
}

std::string WriteMatrixToFile(std::string originalFilePath, std::vector<std::vector<std::string>> matrix) {
    std::string filePath = originalFilePath.substr(0, originalFilePath.find_last_of('.')) + ".Δrle";
    std::ofstream file (filePath);

    if (file) {
        std::bitset<8> binary;

        for (int y = 0; y < matrix.size(); y++) {
            for (int x = 0; x < matrix[y].size(); x++) {
                for (std::size_t i = 0; i < matrix[y][x].size(); i++)
                {
                    file << std::bitset<8>((matrix[y][x]).c_str()[i]);
                }

                file << std::bitset<8>('|');
            }
            file << std::bitset<8>('\n');
        }

        file.close();
    }

    return filePath;
}

std::string WriteMatrixToFile_NoComma_SingleSeparator(std::string originalFilePath, std::vector<std::vector<std::string>> matrix) {
    std::string filePath = originalFilePath.substr(0, originalFilePath.find_last_of('.')) + ".Δrle";
    std::ofstream file (filePath);

    if (file) {
        for (int y = 0; y < matrix.size(); y++) {
            for (int x = 0; x < matrix[y].size(); x++) {
                file << matrix[y][x] << "|";
            }
            file << std::endl;
        }

        file.close();
    }

    return filePath;
}


void WriteMatrixToFile_NoComma_NoSeparator(std::vector<std::vector<std::string>> matrix) {
    std::ofstream file ("/Users/jortrodenburg/Desktop/School/Dissertation/code/dissertation/algorithm/code/test3.Δrle");
    if (file) {
        for (int y = 0; y < matrix.size(); y++) {
            for (int x = 0; x < matrix[y].size(); x++) {
                file << matrix[y][x];
            }
            file << std::endl;
        }

        file.close();
    }
}


void WriteMatrixToFile_TwoSeparatorsAndComma(std::vector<std::vector<std::string>> matrix) {
    std::ofstream file ("/Users/jortrodenburg/Desktop/School/Dissertation/code/dissertation/algorithm/code/sqaures.Δrle");
    if (file) {
        for (int y = 0; y < matrix.size(); y++) {
            for (int x = 0; x < matrix[y].size(); x++) {
                file << "[" << matrix[y][x] << "]";

                // do not add a comma if we are at the end of the file.
                if (y != matrix.size()-1 && x != matrix[y].size() - 1) {
                    file << ",";
                }
            }
            file << std::endl;
        }

        file.close();
    }
}

std::vector<std::string> DecodeFile(std::string fileName) {
    std::vector<std::string> rows;

    std::ifstream file(fileName);
    std::string line;

    while(std::getline(file,line, '\n'))
    {
        std::stringstream lineStream(line);
        std::string element;

        while(std::getline(lineStream,element,'\n'))
            rows.push_back(element);
    }

    return rows;
}


std::vector<std::string> DecodeFileBinary(std::string fileName) {
    std::vector<std::string> rows;

    std::ifstream file(fileName, std::ofstream::binary);
    std::ostringstream os;
    os << file.rdbuf();

    std::string s = os.str();

    std::string output = "";
    std::stringstream sstream(s);
    while (sstream.good())
    {
        std::bitset<8> bits;
        sstream >> bits;
        output += char(bits.to_ulong());
    }

    std::stringstream lineStream(output);
    std::string element;
    while (std::getline(lineStream, element, '\n')) {
        rows.push_back(std::move(element));
    }


    return rows;
}


vector<vector<vector<float>>> DecodeMatrix(vector<std::string> decodedFile) {
    vector<vector<vector<float>>> matrix = vector<vector<vector<float>>>();
    std::map<int, std::vector<float>> dictionary = std::map<int, std::vector<float>>();
    bool populateDictionary = false;

    vector<float> curElement = vector<float>();

    for (int row = 0; row < decodedFile.size(); row++) {
        vector<vector<float>> curRow = vector<vector<float>>();

        std::string columnElement = "";

        for (int column = 0; column <= decodedFile[row].size(); column++) {
            char currentChar = decodedFile[row][column];
            if (populateDictionary) {
                if (currentChar == '#') continue;
            }
            if (column == decodedFile[row].size()) {
                matrix.push_back(curRow);
                break;
            }

            if (currentChar == ',') {
                curElement.push_back(stof(columnElement));
                columnElement = "";
                continue;
            }

            if (currentChar == '|') {
                if (columnElement == "*") {
                    populateDictionary ^= true;
                    columnElement = "";
                    continue;
                }

                curElement.push_back(stof(columnElement));

                if (populateDictionary) {
                    dictionary[dictionary.size()] = curElement;
                } else {
                    curRow.push_back(curElement);
                }

                curElement = vector<float>();
                columnElement = "";

                continue;
            }

            columnElement += currentChar;
        }
    }

    auto flat = flatten(matrix);
    std::vector<int>::size_type size = flat.size();
    for (std::vector<int>::size_type i = 1; i <= size; i++) {
        if (flat[i].size() == 2) { // this means we need to insert a dictionary value
            vector<float> filledFromDictionary = vector<float>();
            filledFromDictionary.push_back(flat[i][0]);

            int dictKey = flat[i][1];
            filledFromDictionary.insert(filledFromDictionary.end(), std::make_move_iterator(dictionary[dictKey].begin()), std::make_move_iterator(dictionary[dictKey].end()));
            flat[i] = filledFromDictionary;
        }
    }



    auto flat2 = flatten(flat);

    vector<vector<float>> filled = vector<vector<float>>();
    vector<float> sizing = vector<float>();

    int width = flat2[0];
    int height = flat2[1];
    sizing.push_back(width);
    sizing.push_back(height);

    filled.push_back(sizing);

    int previousLocation = 0;
    int elementCount = 0;
    vector<float> currentElement = vector<float>();

    size = flat2.size();

    // start at 2 to skip sizing info
    for (std::vector<int>::size_type i = 2; i <= size; i++) {
        if (elementCount == 4) {
            filled.push_back(currentElement);
            currentElement = vector<float>();

            previousLocation = flat2[i - elementCount];
            elementCount = 0;

            if (i == size) break; // last row added to matrix
        }

        if (elementCount == 0) {
            int curLocation = flat2[i];
            if (curLocation != 0 && curLocation != previousLocation + 1) {
                while (curLocation != previousLocation + 1) {
                    auto lastPaintedElement = filled[filled.size()-1];
                    filled.push_back(lastPaintedElement);
                    ++previousLocation;
                }
            }
        } else {
            currentElement.push_back(flat2[i]);
        }

        elementCount++;
    }

    matrix = vector<vector<vector<float>>>();

    int curX = 0;
    int curY = 1;

    vector<vector<float>> row = vector<vector<float>>();
    row.push_back(filled[0]);
    matrix.push_back(row);

    row = vector<vector<float>>();
    for (int i = 1; i < filled.size(); i++) {
        if (curX > width - 1) {
            curX = 0;
            ++curY;
            matrix.push_back(row);
            row = vector<vector<float>>();
        }

        row.push_back(filled[i]);
        ++curX;
    }

    matrix.push_back(row);

    auto lastElement = matrix[matrix.size() - 1][matrix[matrix.size() -1].size()-1];

    if (matrix.size() <= height) {
        row = vector<vector<float>>();

        while (row.size() < width) {
            row.push_back(lastElement);
        }

        while (matrix.size() < height) {
            matrix.push_back(row);
        }
    }

    for (int row = 1; row < matrix.size(); row++) {
        while (matrix[row].size() < width) {
            matrix[row].push_back(lastElement);
        }
    }

    return matrix;
}


template<typename T>
std::vector<T> flatten(std::vector<std::vector<T>> const &vec)
{
    std::vector<T> flattened;
    for (auto const &v: vec) {
        flattened.insert(flattened.end(), v.begin(), v.end());
    }
    return flattened;
}

void Draw(vector<vector<vector<float>>> matrix) {
    SDL_Event event;
    SDL_Renderer *renderer;
    SDL_Window *window;

    SDL_Init(SDL_INIT_VIDEO);

    int width = matrix[0][0][0];
    int height = matrix[0][0][1];

    SDL_CreateWindowAndRenderer(width, height, SDL_WINDOW_RESIZABLE, &window, &renderer);

    SDL_RenderClear(renderer);

    for (int matrixY = 0; matrixY < matrix.size(); matrixY++) {
        if (matrixY == 0) continue; // size info is at matrixY 0

        for (int matrixX = 0; matrixX < matrix[matrixY].size(); matrixX++) {
            std::tuple<float, float, float> rgb = labToRGB(matrix[matrixY][matrixX][0], matrix[matrixY][matrixX][1], matrix[matrixY][matrixX][2]);
            SDL_SetRenderDrawColor(renderer, get<0>(rgb), get<1>(rgb), get<2>(rgb), 255);
            SDL_RenderDrawPoint(renderer, matrixX, matrixY-1);
        }
    }

    SDL_RenderPresent(renderer);
    cout << "rendered";
    while (1) {
        if (SDL_PollEvent(&event) && event.type == SDL_QUIT)
            break;
    }
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
}

// formula from // using https://web.archive.org/web/20111111080001/http://www.easyrgb.com/index.php?X=MATH&H=01#tex1
std::tuple<float, float, float> RGBToLAB(int r, int g, int b) {
    float varR = r / 255.0;
    float varG = g / 255.0;
    float varB = b / 255.0;

    if (varR > 0.04045) {
        varR = pow(((varR + 0.055) / 1.055), 2.4);
    } else {
        varR /= 12.92;
    }

    if (varG > 0.04045) {
        varG = pow(((varG + 0.055) / 1.055), 2.4);
    } else {
        varG /= 12.92;
    }

    if (varB > 0.04045) {
        varB = pow(((varB + 0.055) / 1.055), 2.4);
    } else {
        varB /= 12.92;
    }

    varR *= 100.;
    varG *= 100.;
    varB *= 100.;

    // XYZ to LAB
    float varX = (varR * 0.4124 + varG * 0.3576 + varB * 0.1805) / 95.047;
    float varY = (varR * 0.2126 + varG * 0.7152 + varB * 0.0722) / 100.;
    float varZ = (varR * 0.0193 + varG * 0.1192 + varB * 0.9505) / 108.883;

    if (varX > 0.008856) {
        varX = pow(varX, 1./3.);
    } else {
        varX = (7.787 * varX) + (16. / 116.);
    }

    if (varY > 0.008856) {
        varY = pow(varY,  1./3.);
    } else {
        varY = (7.787 * varY) + (16. / 116.);
    }

    if (varZ > 0.008856) {
        varZ = pow(varZ,  1./3);
    } else {
        varZ = (7.787 * varZ) + (16. / 116.);
    }

    return std::tuple<float, float, float>((116. * varY) - 16., 500. * (varX - varY), 200. * (varY - varZ));
}

// labToRGB converts a LAB value to RGB values. Returned is a tuple<int r, int g, int b>
std::tuple<float, float, float> labToRGB(float l, float a, float b) {
    // first convert LAB to XYZ
    float varY = (l + 16.) / 116.;
    float varX = a / 500. + varY;
    float varZ = varY - b / 200.;

    if (pow(varY, 3) > 0.008856) {
        varY = pow(varY, 3);
    } else {
        varY = (varY - 16. / 116.) / 7.787;
    }

    if (pow(varX, 3) > 0.008856) {
        varX = pow(varX, 3);
    } else {
        varX = (varX - 16. / 116.) / 7.787;
    }

    if (pow(varZ, 3) > 0.008856) {
        varZ = pow(varZ, 3);
    } else {
        varZ = (varZ - 16. / 116.) / 7.787;
    }

    // Convert XYZ to RGB
    varX = (95.047 * varX)/ 100.;
    varY = (100. * varY) / 100.;
    varZ = (108.883 * varZ) / 100.;

    float varR = varX * 3.2406 + varY * -1.5372 + varZ * (-0.4986);
    float varG = varX * (-0.9689) + varY * 1.8758 + varZ * 0.0415;
    float varB = varX * 0.0557 + varY * (-0.2040) + varZ * 1.0570;

    if (varR > 0.0031308) {
        varR = 1.055 * (pow(varR, (1 / 2.4))) - 0.055;
    } else {
        varR *= 12.92;
    }

    if (varG > 0.0031308) {
        varG = 1.055 * (pow(varG, (1 / 2.4))) - 0.055;
    } else {
        varG *= 12.92;
    }

    if (varB > 0.0031308) {
        varB = 1.055 * (pow(varB, (1 / 2.4))) - 0.055;
    } else {
        varB *= 12.92;
    }

    return std::tuple<float,float,float>(varR * 255., varG * 255., varB * 255.);
}


/*
 * CIEDE2000.cpp
 * Part of http://github.com/gfiumara/CIEDE2000 by Gregory Fiumara.
 * See LICENSE for details.
 */
#include <cmath>

/*******************************************************************************
 * Conversions.
 ******************************************************************************/

constexpr double
CIEDE2000::deg2Rad(
        const double deg)
{
    return (deg * (M_PI / 180.0));
}

constexpr double
CIEDE2000::rad2Deg(
        const double rad)
{
    return ((180.0 / M_PI) * rad);
}

double
CIEDE2000::CIEDE2000(
        const LAB &lab1,
        const LAB &lab2)
{
    /*
     * "For these and all other numerical/graphical 􏰀delta E00 values
     * reported in this article, we set the parametric weighting factors
     * to unity(i.e., k_L = k_C = k_H = 1.0)." (Page 27).
     */
    const double k_L = 1.0, k_C = 1.0, k_H = 1.0;
    const double deg360InRad = CIEDE2000::deg2Rad(360.0);
    const double deg180InRad = CIEDE2000::deg2Rad(180.0);
    const double pow25To7 = 6103515625.0; /* pow(25, 7) */

    /*
     * Step 1
     */
    /* Equation 2 */
    double C1 = sqrt((lab1.a * lab1.a) + (lab1.b * lab1.b));
    double C2 = sqrt((lab2.a * lab2.a) + (lab2.b * lab2.b));
    /* Equation 3 */
    double barC = (C1 + C2) / 2.0;
    /* Equation 4 */
    double G = 0.5 * (1 - sqrt(pow(barC, 7) / (pow(barC, 7) + pow25To7)));
    /* Equation 5 */
    double a1Prime = (1.0 + G) * lab1.a;
    double a2Prime = (1.0 + G) * lab2.a;
    /* Equation 6 */
    double CPrime1 = sqrt((a1Prime * a1Prime) + (lab1.b * lab1.b));
    double CPrime2 = sqrt((a2Prime * a2Prime) + (lab2.b * lab2.b));
    /* Equation 7 */
    double hPrime1;
    if (lab1.b == 0 && a1Prime == 0)
        hPrime1 = 0.0;
    else {
        hPrime1 = atan2(lab1.b, a1Prime);
        /*
         * This must be converted to a hue angle in degrees between 0
         * and 360 by addition of 2􏰏 to negative hue angles.
         */
        if (hPrime1 < 0)
            hPrime1 += deg360InRad;
    }
    double hPrime2;
    if (lab2.b == 0 && a2Prime == 0)
        hPrime2 = 0.0;
    else {
        hPrime2 = atan2(lab2.b, a2Prime);
        /*
         * This must be converted to a hue angle in degrees between 0
         * and 360 by addition of 2􏰏 to negative hue angles.
         */
        if (hPrime2 < 0)
            hPrime2 += deg360InRad;
    }

    /*
     * Step 2
     */
    /* Equation 8 */
    double deltaLPrime = lab2.l - lab1.l;
    /* Equation 9 */
    double deltaCPrime = CPrime2 - CPrime1;
    /* Equation 10 */
    double deltahPrime;
    double CPrimeProduct = CPrime1 * CPrime2;
    if (CPrimeProduct == 0)
        deltahPrime = 0;
    else {
        /* Avoid the fabs() call */
        deltahPrime = hPrime2 - hPrime1;
        if (deltahPrime < -deg180InRad)
            deltahPrime += deg360InRad;
        else if (deltahPrime > deg180InRad)
            deltahPrime -= deg360InRad;
    }
    /* Equation 11 */
    double deltaHPrime = 2.0 * sqrt(CPrimeProduct) *
                         sin(deltahPrime / 2.0);

    /*
     * Step 3
     */
    /* Equation 12 */
    double barLPrime = (lab1.l + lab2.l) / 2.0;
    /* Equation 13 */
    double barCPrime = (CPrime1 + CPrime2) / 2.0;
    /* Equation 14 */
    double barhPrime, hPrimeSum = hPrime1 + hPrime2;
    if (CPrime1 * CPrime2 == 0) {
        barhPrime = hPrimeSum;
    } else {
        if (fabs(hPrime1 - hPrime2) <= deg180InRad)
            barhPrime = hPrimeSum / 2.0;
        else {
            if (hPrimeSum < deg360InRad)
                barhPrime = (hPrimeSum + deg360InRad) / 2.0;
            else
                barhPrime = (hPrimeSum - deg360InRad) / 2.0;
        }
    }
    /* Equation 15 */
    double T = 1.0 - (0.17 * cos(barhPrime - CIEDE2000::deg2Rad(30.0))) +
               (0.24 * cos(2.0 * barhPrime)) +
               (0.32 * cos((3.0 * barhPrime) + CIEDE2000::deg2Rad(6.0))) -
               (0.20 * cos((4.0 * barhPrime) - CIEDE2000::deg2Rad(63.0)));
    /* Equation 16 */
    double deltaTheta = CIEDE2000::deg2Rad(30.0) *
                        exp(-pow((barhPrime - deg2Rad(275.0)) / deg2Rad(25.0), 2.0));
    /* Equation 17 */
    double R_C = 2.0 * sqrt(pow(barCPrime, 7.0) /
                            (pow(barCPrime, 7.0) + pow25To7));
    /* Equation 18 */
    double S_L = 1 + ((0.015 * pow(barLPrime - 50.0, 2.0)) /
                      sqrt(20 + pow(barLPrime - 50.0, 2.0)));
    /* Equation 19 */
    double S_C = 1 + (0.045 * barCPrime);
    /* Equation 20 */
    double S_H = 1 + (0.015 * barCPrime * T);
    /* Equation 21 */
    double R_T = (-sin(2.0 * deltaTheta)) * R_C;

    /* Equation 22 */
    double deltaE = sqrt(
            pow(deltaLPrime / (k_L * S_L), 2.0) +
            pow(deltaCPrime / (k_C * S_C), 2.0) +
            pow(deltaHPrime / (k_H * S_H), 2.0) +
            (R_T * (deltaCPrime / (k_C * S_C)) * (deltaHPrime / (k_H * S_H))));

    return (deltaE);
}

/*******************************************************************************
 * Operators.
 ******************************************************************************/

std::ostream&
operator<<(
        std::ostream &s,
        const CIEDE2000::LAB &labColor)
{
    return (s << "CIELAB(" << labColor.l << "," << labColor.a << "," <<
              labColor.b << ")");
}

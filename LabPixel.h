//
// Created by Jort Rodenburg on 12/4/22.
//

#ifndef ALGORITHM_LABPIXEL_H
#define ALGORITHM_LABPIXEL_H


#include "main.h"

class LabPixel {
private:


public:
    CIEDE2000::LAB labEncoding;
    int location;

    LabPixel(CIEDE2000::LAB colors, int location) : labEncoding(colors), location(location){};
};


struct LabPixelCompare {
    bool operator()(const LabPixel &right, const LabPixel &left) const{
        return right.location > left.location;
    }
};

#endif //ALGORITHM_LABPIXEL_H

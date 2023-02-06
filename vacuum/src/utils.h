#ifndef TIDYUP_UTILS_H
#define TIDYUP_UTILS_H

#include <vector>

#include <iostream>
#include <sstream>
#include <fstream>

bool parse_pg_sample_line(std::string& line){ // returns true if is sample pg line
    std::stringstream *line_stream = new std::stringstream(line);
    std::string col;

    // make sure it's CO
    std::getline(*line_stream, col, '\t');
    if(std::strcmp(col.c_str(),"@CO")!=0){
        delete line_stream;
        return false;
    }

    // check if ID == SAMPLE
    std::getline(*line_stream,col,':');
    if(std::strcmp(col.c_str(),"SAMPLE")!=0){
        delete line_stream;
        return false;
    }

    std::getline(*line_stream,col,'\n');
    line = col;
    delete line_stream;
    return true;
}

#endif //TIDYUP_UTILS_H

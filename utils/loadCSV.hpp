/*
*   load the correspondences from the CSV file
*   by R. Falque
*   21/11/2018
*/

#ifndef LOADCSV_H
#define LOADCSV_H

#include <fstream>      // std::ifstream
#include <vector>
#include <Eigen/Core>

inline int load_csv (const std::string & path, Eigen::MatrixXd &sources, Eigen::MatrixXd &targets) {
    std::ifstream indata;
    indata.open(path);
    std::string line;
    std::vector<double> values;
    uint rows = 0;
    while (std::getline(indata, line)) {
        std::stringstream lineStream(line);
        std::string cell;
        while (std::getline(lineStream, cell, ',')) {
            values.push_back(std::stod(cell));
        }
        ++rows;
    }
    Eigen::Map<const Eigen::Matrix<typename Eigen::MatrixXd::Scalar, 
                                     Eigen::MatrixXd::RowsAtCompileTime, 
                                     Eigen::MatrixXd::ColsAtCompileTime, 
                                     Eigen::RowMajor>> temp(values.data(), rows, values.size()/rows);

    //std::cout << temp << std::endl;
    //std::cout << temp.leftCols(3) << std::endl;
    //std::cout << temp.rightCols(3) << std::endl;

    sources = temp.leftCols(3);
    targets = temp.rightCols(3);

    return 0;
};

#endif

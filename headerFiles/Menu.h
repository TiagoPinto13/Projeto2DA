//
// Created by Diogo Vieira on 01/05/2024.
//

#ifndef PROJ2DA_MENU_H
#define PROJ2DA_MENU_H

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>

#include "Data.h"

class Menu {
public:
    Menu();
    void drawMenu();
    void drawTop();
    void drawBottom();
    void drawBacktracking(std::string vertex_id);
    void drawTriangular(std::string vertex_id);
    void drawCluster(std::string vertex_id);
    void drawMST(std::string vertex_id);
    void drawApproximationAnalysis(std::string vertex_id);

private:
    Data data_;
};

#endif //PROJ2DA_MENU_H
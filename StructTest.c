#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAXSIZE 10000
/* 本程序是采用对应的struct 来进行方程组的求解*/

/*创建结构，来处理对应的矩阵*/

struct Matrix{
    double data[MAXSIZE];
    int row;
    int col;
};




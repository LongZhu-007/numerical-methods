#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define MAXSIZE 10000

double *Gauss_Jordan(double *Arr, int row, int col, double *b, int cl);
double *Upper_Solver(double *Arr, int row, int col, double *b, int cl);
/* 定义一个结构体，用来储存矩阵的对应属性*/

struct Matrix{
    double data[MAXSIZE];
    int row;
    int col;
};

// 下面来创建函数进行矩阵的基本操作
double *Mat_Add(const struct Matrix *mat1, const struct Matrix *mat2){
    // 首先来判断只有同型矩阵才可以进行相加减
    if (mat1->row!=mat2->row || mat1->col!=mat2->col){
        printf("Sorry! These two matrix you input are not valid for Mat_Add.\n");
        exit(-1);
    }
    int i, size;
    double *res;
    size = mat1->row * mat1->col;
    res = (double *)calloc(size, sizeof(double));

    for (i=0; i<size; i++){
        *(res+i) = mat1->data[i] + mat2->data[i];
    }
    // 计算结束之后，来创建一个结构体，并返回
    return res;
}

// 下面来创建函数进行矩阵的相减
double *Mat_Sub(const struct Matrix *mat1, const struct Matrix *mat2){
    if (mat1->row!=mat2->row || mat1->col!=mat2->col){
        printf("Sorry! These two matrix you input are not valid for Mat_Add.\n");
        exit(-1);
    }
    int i, size;
    double *res;
    size = mat1->row * mat1->col;
    res = (double *)calloc(size, sizeof(double));

    for (i=0; i<size; i++){
        *(res+i) = mat1->data[i] + mat2->data[i];
    }
    // 计算结束之后，来创建一个结构体，并返回
    return res;
}

// 下面进行矩阵的转置
struct *Mat_Transpose(const struct Matrix *mat1){
    // 获取传入参数的row,col
    // 创建一个要返回的结构体
    struct Matrix *a;
    int n, m;
    for (n=0; n<mat1->row; n++){
        for(m=0; m<mat1->col; m++){
            a->data[m*mat1->row+n] = mat1->data[n*mat1->row+m];
        }
    }
    a->row = mat1->col;
    a->col = mat1->row;
    return a;
}

// 下面来创建一个方法，来实现列主元方程组的求解
double *Gauss_Jordan(double *Arr, int row, int col, double *b, int cl){
    // 首先进行列主元的寻找。
    if (col != cl){
        printf("Sorry! The arguments named Arr and b has not same column.\n");
        exit(-1);
    }
    int i, j, d, max_index;
    double temp;
    for (i=0; i<row; i++){
        max_index = i;
        for (j=i+1; j<row; j++){
            if (abs(Arr[j*row+i]) > abs(Arr[max_index*row+i])){
                max_index = j;
            }
        }
        // 循环结束后，得到对应的主元所在的行
        // 进行交换
        for (j=i; j<col; j++){
            temp = Arr[max_index*row+j];
            Arr[max_index*row+j] = Arr[i*row+j];
            Arr[i*row+j] = temp;
        }
        // 交换对应的b中的元素
        temp = b[max_index];
        b[max_index] = b[i];
        b[i] = temp;

        // 交换完毕之后，可以来进行消元计算了
        for (j=i+1; j<row; j++){
            temp = *(Arr+j*row+i) / *(Arr+i*row+i);
            for (d=i; d<col; d++){
                *(Arr+j*row+d) -= temp * *(Arr+i*row+d);
            }
            // 循环结束之后，同样需要消去对应的b数组中的元素
            *(b+j) -= *(b+i) * temp;
        }
    }
    // 循环结束后，便可以来进行求解计算了
    return Upper_Solver(Arr, row, col, b, cl);
}

// 下面来创建一个函数，进行上三角形函数的求解
double *Upper_Solver(double *Arr, int row, int col, double *b, int cl){
    int i, j;
    double *res;
    double sum_value;
    res = (double *)calloc(row, sizeof(double));
    for (i=row-1; i>=0; i--){
        sum_value = 0;
        for (j=i+1; j<col; j++){
            sum_value += *(Arr+i*row+j) * *(res+j);
        }
        // 循环结束之后，后面的求和计算便结束了
        *(res+i) = (*(b+i) - sum_value) / *(Arr+i*row+i);
    }
    return res;
}

// 下面来创建一个函数，来进行矩阵的三角分解
double *Matrix_LU(double *Arr, int row, int col, double *b, int cl){
    // 首先来判断当前的传入的两个参数是否满足方程组求解的条件
    // 要求系数矩阵Arr的列 col 必须与b数组的列数相等 cl
    if (col != cl){
        printf("Sorry! The arguments named Arr or b must be the same column.\n");
    }
    int i, j, d;
    double sum_value;
    /* LU分解的基本要求是先求解U矩阵中的行，再去求解L矩阵中的列*/
    for (i=0; i<row; i++){
        for (j=i; j<col; j++){
            // 这里的j是进行搜索的列索引
            // 下面需要来创建一个当前i行上面的行的求和
            sum_value = 0;
            for(d=0; d<i; d++){ // 这里的d对应的是U中的列，j是U中对应的行
                sum_value += *(Arr+i*row+j) * *(Arr+j*row+d);
            }
            // 循环结束之后，对于u中的中间求和的过程便计算完成了
            *(Arr+i*row+j) -= sum_value;
        }
        // 待U中的对应的一行计算完成之后，下面进行L中的列计算
        // 对于L中的中间变量-->求和是需要进行当前列的前面的列的求和
        // 那么这里的i的遍历将会变为矩阵中的列索引
        for (d=i+1; d<row; d++){
            // 这里的d是对应的L中的行
            // 需要来创建一个中间索引变量j, j是L中的当前寻找的列索引
            sum_value = 0;
            for(j=0; j<i; j++){
                sum_value += *(Arr+d*row+j) * *(Arr+j*row+i);
            }
            *(Arr+d*row+i) = (*(Arr+d*row+i)/ *(Arr+i*row+i));
        }
    }
    // 首先来进行下三角型方程组求解
    double *y;
    y = Down_Solver(Arr, row, col, b, cl);
    double *x;
    x = Upper_Solver(Arr, row, col, y, cl);
    free(y);
    return x;
}

// 我们在这里来创建一个函数，来计算下三角形系数方程组
double *Down_Solver(double *Arr, int row, int col, double *b, int cl){
    if (col != cl){
        printf("Sorry! These two mats you input should be the same size.\n");
        exit(-1);
    }
    int n, m;
    double *res;
    double sum_value;
    res = (double *)calloc(row, sizeof(double));
    for (n=0; n<row; n++){
        // 因为下三角方程组对应的res的值是从第一个元素开始计算的
        sum_value = 0;
        for (m=0; m<n; m++){
            sum_value += *(Arr+n*row+m) * *(res+m);
        }
        // 当循环结束之后，对应的解便计算出来了
        *(res+n) = (*(b+n) - sum_value);
    }
    return res;
}

// 下面来创建一个函数，来实现列主元LU分解
double *LU(double *Arr, int row, int col, double *b, int cl){
    if (col != cl){
        printf("Sorry! These mats you input are invalid.\n");
    }
    /*
        选主元三角分解法的主要步骤

    */
    int n, m,d, max_index;
    double sum_value, temp;
    // 因为是列主元消去法，这里我们首先要做的是处理对应的列
    for (n=0; n<row; n++){
        // 然后遍历当前主行所对应的列上的下面的所有的行
        max_index = n;
        for (m=n; m<row; m++){
            sum_value = 0;
            for (d=0; d<n; d++){
                sum_value += *(Arr+m*row+d) * *(Arr+d*row+n);
            }
            // 然后计算出对应的A_m_n上
            *(Arr+m*row+n) -= sum_value;
            if (abs(*(Arr+m*row+n))>abs(*(Arr+max_index*row+n))){
                max_index = m;
            }
        }
        // 然后当循环结束之后，便得到了主元所在的行
        // 下面进行数据的交换
        if (max_index != n){
            for (d=0; d<col; d++){
                temp = *(Arr+max_index*row+d);
                *(Arr+max_index*row+d) = *(Arr+n*row+d);
                *(Arr+n*row+d) = temp;
            }
            // 交换对应的b上的数据
            temp = *(b+max_index);
            *(b+max_index) = *(b+n);
            *(b+n) = temp;
        }
        // 交换完毕之后，下面进行l_i_r的计算
        // 且首先来计算出l_i_r数据, d是对应的行
        for (d=n+1; d<row; d++){
            *(Arr+d*row+n) /= *(Arr+n*row+n);
        }
        // 然后下面进行对应U中行的计算
        // 当转向行计算的时候，n在U中表示的是U当前遍历到的行
        // 下面该去寻找指定的列元素，下面的m是对应的列索引
        for (m=n+1; m<col; m++){
            // 遍历当前n行上面的所有的行
            sum_value = 0;
            for (d=0; d<n; d++){
                sum_value += *(Arr+n*row+d) * *(Arr+d*row+m);
            }
            // 然后循环结束之后，来计算得到对应的U_n_m
            *(Arr+n*row+m) -= sum_value;
        }
    }
    double *y;
    y = Down_Solver(Arr, row, col, b, cl);
    double *x;
    x = Upper_Solver(Arr, row, col, y, cl);
    free(y);
    return x;
}
// 我们在这里来创建一个方法，来通过Gauss_Jordan求解矩阵的逆
void Mat_Inv(double *data, int row){
    // 首先来创建一个内存空间，来存储我们所指定的B矩阵
    double *B;
    B = (double *)calloc(row, sizeof(double));
    int r, d, cl, max_index, count;
    double temp, effc;
    for (r=0; r<row; r++){
       // 找到最大主元所在的行
       memset(B, 0, row);
       *(B+r) = 1;
       max_index = r;
       for (d=r+1; d<row; d++){
            if (abs(*(data+d*row+r)) > abs(*(data+max_index*row+r))){
                max_index = r;
            }
       }
       // 循环结束后，便得到了最大的主元所在的行
       if (max_index !=r){
            // 如果对应的max_index 不等于 r, 那么我们来交换
            for (d=r; d<row; d++){
                temp = *(data+r*row+d);
                *(data+r*row+d) = *(data+max_index*row+d);
                *(data+max_index*row+d) = temp;
            }
            // 循环结束之后，那么我们便需要来交换对应的B中的elema
            temp = *(B+max_index);
            *(B+max_index) = *(B+r);
            *(B+r) = temp;
            count += 1;

       }
       // 这里首先判断当前主元是否为0
       if (abs(*(data+r*row+r))<1e-5){
            printf("Sorry! 您传入的矩阵不是满秩矩阵，无法求逆!\n");
            exit(-1);
       }

       // 然后下面我们需要来进行消元处理了, 消元必须从第一行开始执行
       for (d=0; d<row; d++){
            if (d==r){
                effc = *(data+d*row+r);

            }else{
                effc = *(data+d*row+r) / *(data+r*row+r);
            }
            // 然后下面我们需要来处理对应的行上的元素
            // 来处理对应的列上的元素
            for (cl=r; cl<row; cl++){
                *(data+d*row+cl) -= *(data+r*row+cl) * effc;
            }
            // 循环结束之后，来处理对应的B中的数据
            *(B+d) -= *(B+r)* effc;
       }
       // 循环结束之后，我们需要来将B中的数据写入到mat中的r列
       for (d=0; d<row; d++){
            *(data+d*row+r) = *(B+d);
       }
    }
    // 当循环结束之后，来判断对应的count
    if (count % 2 == 1){
        int start, end;
        start = 0;
        end = row-1;
        while (end > start){
            // 然后来进行对应的行遍历
            for (r=0; r<row; r++){
                temp = *(data+r*row+start);
                *(data+r*row+start) = *(data+r*row+end);
                *(data+r*row+end) = temp;
            }
            start += 1;
            end -= 1;
        }
    }
    // 不应该有返回，因为是直接修改的原矩阵
}

// 创建一个函数，来打印输出对应的1D数组
void print_1DArray(const double *Arr, int row){
    int n;
    for (n=0; n<row; n++){
        printf("%lf\t", *(Arr+n));
    }
    printf("\n");
}
void print_2DArray(const double *Arr, int row, int col){
    int n, m;
    for (n=0; n<row; n++){
        for (m=0; m<col; m++){
            printf("%lf\t", *(Arr+n*row+m));
        }
        printf("\n");
    }
}

int main()
{
    printf("Hello world!\n");
    double Arr[] = {1, 1, 1, 0, 4, -1, 0, 0, -2};
    double b[] = {6, 5, -6};
    printf("打印输出您所输入的数组和矩阵:\n");
    print_2DArray(Arr, 3, 3);
    printf("打印输出方程组的右端项:\n");
    print_1DArray(b, 3);
    double *res;
    res = Gauss_Jordan(Arr, 3, 3, b, 3);
    printf("采用高斯消去法来计算出对应的结果:\n");
    print_1DArray(res, 3);
    free(res);
    return 0;
}




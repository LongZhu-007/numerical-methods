#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define MAXSIZE 10000

double *Gauss_Jordan(double *Arr, int row, int col, double *b, int cl);
double *Upper_Solver(double *Arr, int row, int col, double *b, int cl);
/* ����һ���ṹ�壬�����������Ķ�Ӧ����*/

struct Matrix{
    double data[MAXSIZE];
    int row;
    int col;
};

// �����������������о���Ļ�������
double *Mat_Add(const struct Matrix *mat1, const struct Matrix *mat2){
    // �������ж�ֻ��ͬ�;���ſ��Խ�����Ӽ�
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
    // �������֮��������һ���ṹ�壬������
    return res;
}

// �����������������о�������
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
    // �������֮��������һ���ṹ�壬������
    return res;
}

// ������о����ת��
struct *Mat_Transpose(const struct Matrix *mat1){
    // ��ȡ���������row,col
    // ����һ��Ҫ���صĽṹ��
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

// ����������һ����������ʵ������Ԫ����������
double *Gauss_Jordan(double *Arr, int row, int col, double *b, int cl){
    // ���Ƚ�������Ԫ��Ѱ�ҡ�
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
        // ѭ�������󣬵õ���Ӧ����Ԫ���ڵ���
        // ���н���
        for (j=i; j<col; j++){
            temp = Arr[max_index*row+j];
            Arr[max_index*row+j] = Arr[i*row+j];
            Arr[i*row+j] = temp;
        }
        // ������Ӧ��b�е�Ԫ��
        temp = b[max_index];
        b[max_index] = b[i];
        b[i] = temp;

        // �������֮�󣬿�����������Ԫ������
        for (j=i+1; j<row; j++){
            temp = *(Arr+j*row+i) / *(Arr+i*row+i);
            for (d=i; d<col; d++){
                *(Arr+j*row+d) -= temp * *(Arr+i*row+d);
            }
            // ѭ������֮��ͬ����Ҫ��ȥ��Ӧ��b�����е�Ԫ��
            *(b+j) -= *(b+i) * temp;
        }
    }
    // ѭ�������󣬱������������������
    return Upper_Solver(Arr, row, col, b, cl);
}

// ����������һ�������������������κ��������
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
        // ѭ������֮�󣬺������ͼ���������
        *(res+i) = (*(b+i) - sum_value) / *(Arr+i*row+i);
    }
    return res;
}

// ����������һ�������������о�������Ƿֽ�
double *Matrix_LU(double *Arr, int row, int col, double *b, int cl){
    // �������жϵ�ǰ�Ĵ�������������Ƿ����㷽������������
    // Ҫ��ϵ������Arr���� col ������b������������ cl
    if (col != cl){
        printf("Sorry! The arguments named Arr or b must be the same column.\n");
    }
    int i, j, d;
    double sum_value;
    /* LU�ֽ�Ļ���Ҫ���������U�����е��У���ȥ���L�����е���*/
    for (i=0; i<row; i++){
        for (j=i; j<col; j++){
            // �����j�ǽ���������������
            // ������Ҫ������һ����ǰi��������е����
            sum_value = 0;
            for(d=0; d<i; d++){ // �����d��Ӧ����U�е��У�j��U�ж�Ӧ����
                sum_value += *(Arr+i*row+j) * *(Arr+j*row+d);
            }
            // ѭ������֮�󣬶���u�е��м���͵Ĺ��̱���������
            *(Arr+i*row+j) -= sum_value;
        }
        // ��U�еĶ�Ӧ��һ�м������֮���������L�е��м���
        // ����L�е��м����-->�������Ҫ���е�ǰ�е�ǰ����е����
        // ��ô�����i�ı��������Ϊ�����е�������
        for (d=i+1; d<row; d++){
            // �����d�Ƕ�Ӧ��L�е���
            // ��Ҫ������һ���м���������j, j��L�еĵ�ǰѰ�ҵ�������
            sum_value = 0;
            for(j=0; j<i; j++){
                sum_value += *(Arr+d*row+j) * *(Arr+j*row+i);
            }
            *(Arr+d*row+i) = (*(Arr+d*row+i)/ *(Arr+i*row+i));
        }
    }
    // �����������������ͷ��������
    double *y;
    y = Down_Solver(Arr, row, col, b, cl);
    double *x;
    x = Upper_Solver(Arr, row, col, y, cl);
    free(y);
    return x;
}

// ����������������һ����������������������ϵ��������
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
        // ��Ϊ�����Ƿ������Ӧ��res��ֵ�Ǵӵ�һ��Ԫ�ؿ�ʼ�����
        sum_value = 0;
        for (m=0; m<n; m++){
            sum_value += *(Arr+n*row+m) * *(res+m);
        }
        // ��ѭ������֮�󣬶�Ӧ�Ľ����������
        *(res+n) = (*(b+n) - sum_value);
    }
    return res;
}

// ����������һ����������ʵ������ԪLU�ֽ�
double *LU(double *Arr, int row, int col, double *b, int cl){
    if (col != cl){
        printf("Sorry! These mats you input are invalid.\n");
    }
    /*
        ѡ��Ԫ���Ƿֽⷨ����Ҫ����

    */
    int n, m,d, max_index;
    double sum_value, temp;
    // ��Ϊ������Ԫ��ȥ����������������Ҫ�����Ǵ����Ӧ����
    for (n=0; n<row; n++){
        // Ȼ�������ǰ��������Ӧ�����ϵ���������е���
        max_index = n;
        for (m=n; m<row; m++){
            sum_value = 0;
            for (d=0; d<n; d++){
                sum_value += *(Arr+m*row+d) * *(Arr+d*row+n);
            }
            // Ȼ��������Ӧ��A_m_n��
            *(Arr+m*row+n) -= sum_value;
            if (abs(*(Arr+m*row+n))>abs(*(Arr+max_index*row+n))){
                max_index = m;
            }
        }
        // Ȼ��ѭ������֮�󣬱�õ�����Ԫ���ڵ���
        // ����������ݵĽ���
        if (max_index != n){
            for (d=0; d<col; d++){
                temp = *(Arr+max_index*row+d);
                *(Arr+max_index*row+d) = *(Arr+n*row+d);
                *(Arr+n*row+d) = temp;
            }
            // ������Ӧ��b�ϵ�����
            temp = *(b+max_index);
            *(b+max_index) = *(b+n);
            *(b+n) = temp;
        }
        // �������֮���������l_i_r�ļ���
        // �������������l_i_r����, d�Ƕ�Ӧ����
        for (d=n+1; d<row; d++){
            *(Arr+d*row+n) /= *(Arr+n*row+n);
        }
        // Ȼ��������ж�ӦU���еļ���
        // ��ת���м����ʱ��n��U�б�ʾ����U��ǰ����������
        // �����ȥѰ��ָ������Ԫ�أ������m�Ƕ�Ӧ��������
        for (m=n+1; m<col; m++){
            // ������ǰn����������е���
            sum_value = 0;
            for (d=0; d<n; d++){
                sum_value += *(Arr+n*row+d) * *(Arr+d*row+m);
            }
            // Ȼ��ѭ������֮��������õ���Ӧ��U_n_m
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
// ����������������һ����������ͨ��Gauss_Jordan���������
void Mat_Inv(double *data, int row){
    // ����������һ���ڴ�ռ䣬���洢������ָ����B����
    double *B;
    B = (double *)calloc(row, sizeof(double));
    int r, d, cl, max_index, count;
    double temp, effc;
    for (r=0; r<row; r++){
       // �ҵ������Ԫ���ڵ���
       memset(B, 0, row);
       *(B+r) = 1;
       max_index = r;
       for (d=r+1; d<row; d++){
            if (abs(*(data+d*row+r)) > abs(*(data+max_index*row+r))){
                max_index = r;
            }
       }
       // ѭ�������󣬱�õ���������Ԫ���ڵ���
       if (max_index !=r){
            // �����Ӧ��max_index ������ r, ��ô����������
            for (d=r; d<row; d++){
                temp = *(data+r*row+d);
                *(data+r*row+d) = *(data+max_index*row+d);
                *(data+max_index*row+d) = temp;
            }
            // ѭ������֮����ô���Ǳ���Ҫ��������Ӧ��B�е�elema
            temp = *(B+max_index);
            *(B+max_index) = *(B+r);
            *(B+r) = temp;
            count += 1;

       }
       // ���������жϵ�ǰ��Ԫ�Ƿ�Ϊ0
       if (abs(*(data+r*row+r))<1e-5){
            printf("Sorry! ������ľ��������Ⱦ����޷�����!\n");
            exit(-1);
       }

       // Ȼ������������Ҫ��������Ԫ������, ��Ԫ����ӵ�һ�п�ʼִ��
       for (d=0; d<row; d++){
            if (d==r){
                effc = *(data+d*row+r);

            }else{
                effc = *(data+d*row+r) / *(data+r*row+r);
            }
            // Ȼ������������Ҫ�������Ӧ�����ϵ�Ԫ��
            // �������Ӧ�����ϵ�Ԫ��
            for (cl=r; cl<row; cl++){
                *(data+d*row+cl) -= *(data+r*row+cl) * effc;
            }
            // ѭ������֮���������Ӧ��B�е�����
            *(B+d) -= *(B+r)* effc;
       }
       // ѭ������֮��������Ҫ����B�е�����д�뵽mat�е�r��
       for (d=0; d<row; d++){
            *(data+d*row+r) = *(B+d);
       }
    }
    // ��ѭ������֮�����ж϶�Ӧ��count
    if (count % 2 == 1){
        int start, end;
        start = 0;
        end = row-1;
        while (end > start){
            // Ȼ�������ж�Ӧ���б���
            for (r=0; r<row; r++){
                temp = *(data+r*row+start);
                *(data+r*row+start) = *(data+r*row+end);
                *(data+r*row+end) = temp;
            }
            start += 1;
            end -= 1;
        }
    }
    // ��Ӧ���з��أ���Ϊ��ֱ���޸ĵ�ԭ����
}

// ����һ������������ӡ�����Ӧ��1D����
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
    printf("��ӡ����������������;���:\n");
    print_2DArray(Arr, 3, 3);
    printf("��ӡ�����������Ҷ���:\n");
    print_1DArray(b, 3);
    double *res;
    res = Gauss_Jordan(Arr, 3, 3, b, 3);
    printf("���ø�˹��ȥ�����������Ӧ�Ľ��:\n");
    print_1DArray(res, 3);
    free(res);
    return 0;
}




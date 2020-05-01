using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MatrixOperation
{
    public class General
    {
        //矩阵A与矩阵B相乘
        public static double [,] Multiply_AxB(double [,] matrixA, double[,] matrixB, int rowCountOfA, int colCountOfB, int CommonNum)
        {
            //矩阵A的行数
            //矩阵A的列数与B的行数
            //矩阵B的列数
            int l = rowCountOfA;
            int m = CommonNum;
            int n = colCountOfB;

            //矩阵A
            //矩阵B
            //A与B相乘的结果
            double[,] A = matrixA;
            double[,] B = matrixB;
            double[,] AxB = new double[l, n];

            for(int i = 0; i < l; i++)
            {
                for(int j = 0; j < n; j++)
                {
                    double[] theRow = new double[m];
                    double[] theCol = new double[m];
                    double sum = 0;
                    for(int k = 0; k < m; k++)
                    {
                        theRow[k] = A[i, k];
                        theCol[k] = B[k, j];
                        sum = sum + theRow[k] * theCol[k];
                    }

                    AxB[i, j] = sum;

                }
            }

            //返还结果值
            return AxB;
        }
        //矩阵A与矩阵B相乘，结果再与矩阵C相乘
        public static double[,] Multiply_AxBxC(double[,] matrixA, double[,] matrixB, double[,] matrixC, int rowCountOfA, int colCountOfC, int CommonNum1st, int CommonNum2nd)
        {
            //矩阵A的行数
            //矩阵A的列数与B的行数
            //矩阵B的列数与C的行数
            //矩阵C的列数
            int l = rowCountOfA;
            int m = CommonNum1st;
            int n = CommonNum2nd;
            int o = colCountOfC;

            //矩阵A
            //矩阵B
            //矩阵C
            //A与B相乘的结果
            double[,] A = matrixA;
            double[,] B = matrixB;
            double[,] AxB = new double[l, n];
            double[,] C = matrixC;
            double[,] AxBxC = new double[l, o];
            
            //将矩阵A与矩阵B相乘
            //得到临时结果矩阵AxB
            for (int i = 0; i < l; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    double[] theRow = new double[m];
                    double[] theCol = new double[m];
                    double sum = 0;
                    for (int k = 0; k < m; k++)
                    {
                        theRow[k] = A[i, k];
                        theCol[k] = B[k, j];
                        sum = sum + theRow[k] * theCol[k];
                    }

                    AxB[i, j] = sum;

                }
            }

            //将临时结果矩阵AxB与矩阵C相乘
            //得到最终的结果矩阵AxBxC
            for (int i = 0; i < l; i++)
            {
                for(int j = 0; j < o; j++)
                {
                    double[] theRow = new double[n];
                    double[] theCol = new double[n];
                    double sum = 0;
                    for (int k = 0; k < n; k++)
                    {
                        theRow[k] = AxB[i, k];
                        theCol[k] = C[k, j];
                        sum = sum + theRow[k] * theCol[k];
                    }

                    AxBxC[i, j] = sum;
                }
            }

            //返还结果值
            return AxBxC;
        }
        //矩阵转置
        public static double[,] Transpose(double[,] matrixA, int rowCountOfA, int colCountOfA)
        {
            double[,] A = matrixA;
            int rowCountOfA_transpose = colCountOfA;
            int colCountOfA_transpose = rowCountOfA;
            int rCOAt = rowCountOfA_transpose;
            int cCOAt = colCountOfA_transpose;
            double[,] A_transpose = new double[rCOAt, cCOAt];

            for(int i = 0; i < rCOAt; i++)
            {
                for(int j = 0; j < cCOAt; j++)
                {
                    A_transpose[i, j] = A[j, i];
                }
            }

            //返回结果值
            return A_transpose;

        }
        //矩阵求逆
        public static double[,] Inverse(double[,] matrixB, int orderNum)
        {
            //判断是否满秩
            bool IsFullRank = true;
            //n为阶级
            int n = orderNum;


            //####赋值####
            //矩阵B
            //矩阵B的逆矩阵
            //单位矩阵E和矩阵B组成的矩阵
            double[,] B_normal = matrixB;
            double[,] B_inverse = new double[n, n];
            double[,] EandB_normal = new double[n, 2 * n];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    if (i == j)
                        EandB_normal[i, j] = 1;
                    else
                        EandB_normal[i, j] = 0;

                }
                for (int k = n; k < 2 * n; k++)
                {
                    EandB_normal[i, k] = B_normal[i, k - n];
                }
            }



            //####计算####
            //中间变量数组，用于临时盛装值
            double[] rowHaveZero = new double[2 * n];
            //EB矩阵右边的n*n变为上三角矩阵
            for (int j = n; j < 2 * n; j++)
            {
                int firstRowN = j - n;
                int lastRowN = n;
                int colCount = 2 * n;
                //把EB中索引为j的列的值化为1
                for (int i = firstRowN; i < lastRowN; i++)
                {
                    //如果EBijNum值为0，就把0所在的行与此刻最后一行调换位置
                    //并且循环变量i的终止值减去1,直到EBijNum值不为0
                    //最多调换到0所在的行的下一行
                    double EBijNum = EandB_normal[i, j];
                    while (EBijNum == 0 && lastRowN > i + 1)
                    {
                        for (int k = 0; k < colCount; k++)
                        {
                            rowHaveZero[k] = EandB_normal[i, k];
                            EandB_normal[i, k] = EandB_normal[lastRowN - 1, k];
                            EandB_normal[lastRowN - 1, k] = rowHaveZero[k];
                        }
                        lastRowN -= 1;
                        EBijNum = EandB_normal[i, j];
                    }
                    //如果while循环是由第二个判断跳出
                    //即EBijNum始终为0
                    if (EBijNum == 0)
                    {
                        //循环变量i的终止值再减去1，然后跳出循环
                        lastRowN -= 1;
                        break;
                    }
                    //如果为负数，该行变号
                    if (EBijNum < 0)
                    {
                        for (int k = 0; k < colCount; k++)
                        {
                            EandB_normal[i, k] = (-1) * EandB_normal[i, k];
                        }
                        EBijNum = EandB_normal[i, j];
                    }
                    //将该值变为1，则其余值都除以EBijNum
                    for (int k = 0; k < colCount; k++)
                    {
                        EandB_normal[i, k] = EandB_normal[i, k] / EBijNum;
                    }
                }

                //自n列起，每列只保留一个1，呈阶梯状
                int secondRowN = firstRowN + 1;
                for (int i = secondRowN; i < lastRowN; i++)
                {
                    for (int k = 0; k < colCount; k++)
                    {
                        EandB_normal[i, k] = EandB_normal[i, k] - EandB_normal[firstRowN, k];
                    }
                }

                if (lastRowN == firstRowN)
                {
                    //矩阵不满秩
                    IsFullRank = false;
                    break;
                }
            }
            //不满秩，结束运算
            if (!IsFullRank)
            {
                double[,] error = new double[n, n];
                for (int i = 0; i < n; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        error[i, j] = 0;
                    }
                }
                //返还值均为0的矩阵
                return error;

            }

            //将上三角矩阵变为单位矩阵
            for (int j = 2 * n - 1; j > n; j--)
            {
                //firstRowN为参考行
                //secondRowN为运算行
                int firstRowN = j - n;
                int secondRowN = firstRowN - 1;
                int colCount = j + 1;
                //从最后一列起，每列只保留一个1，其余减为0
                for (int i = secondRowN; i > -1; i--)
                {
                    double EBijNum = EandB_normal[i, j];
                    for (int k = 0; k < colCount; k++)
                    {
                        EandB_normal[i, k] = EandB_normal[i, k] - EandB_normal[firstRowN, k] * EBijNum;
                    }
                }

            }



            //####提取逆矩阵####
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    B_inverse[i, j] = EandB_normal[i, j];
                }
            }

            //返还结果值
            return B_inverse;
        }
    }
}

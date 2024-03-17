# GeneralTrnansmissionCombining
# Description
The file GeneralTrnansmissionCombining.m contains a function that generates the real part of both the transmission and combining matrices for a $N_t\times N_r$ multi-hop MIMO channel.
The function outpus 
1. transmitMatrixMIMO which is a symbolic matrix that holds half of the transmission matrix. The full trnanmission matrix is
$\begin{bmatrix}X\\\Y\end{bmatrix}$
$transmitMatrix=\begin{bmatrix} transmitMatrixMIMOReal \\\ conj(transmitMatrixMIMOReal) \end{bmatrix}$
3. combMatReal which is a symbolc matrix that holds the half of the combining matrix. The full combining matrix is combMat=[combMatReal;conj(combMatReal)].

# Output Example
1. N_t=N_r=2
transmitMatrixMIMOReal = 
[x1,  0;
x2,  0;
x3,  0;
x4,  0;
 0, x1;
 0, x2;
 0, x3
 0, x4]
combMatReal = 
[y_1_1 + y_4_2 + y_6_1 + y_7_2;
y_2_1 - y_3_2 - y_5_1 + y_8_2;
y_2_2 + y_3_1 - y_5_2 - y_8_1;
y_4_1 - y_1_2 - y_6_2 + y_7_1]

2. N_t=2, N_r=3
transmitMatrixMIMOReal = 
[x1,  0;
x2,  0;
x3,  0;
x4,  0;
x5,  0
x6,  0;
x7,  0;
x8,  0;
 0, x1;
 0, x2;
 0, x3;
 0, x4;
 0, x5;
 0, x6;
 0, x7;
 0, x8]
combMatReal = 
[y_1_1 + y_4_3 + y_6_2 + y_10_2 + y_11_3 + y_15_1;
 y_2_1 - y_3_3 - y_5_2 - y_9_2 + y_12_3 - y_16_1;
 y_2_3 + y_3_1 - y_8_2 - y_9_3 - y_12_2 + y_13_1;
y_4_1 - y_1_3 + y_7_2 - y_10_3 + y_11_2 - y_14_1;
y_2_2 + y_5_1 + y_8_3 - y_11_1 - y_14_2 + y_15_3;
y_6_1 - y_1_2 - y_7_3 + y_12_1 + y_13_2 + y_16_3;
 y_6_3 - y_4_2 + y_7_1 - y_9_1 - y_13_3 + y_16_2;
y_3_2 - y_5_3 + y_8_1 + y_10_1 - y_14_3 - y_15_2]

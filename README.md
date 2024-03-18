# GeneralTrnansmissionCombining
# Description
The Matlab file GeneralTrnansmissionCombining.m contains a function that generates the real part of both the transmission and combining matrices for a $N_t\times N_r$ multi-hop MIMO channel.

The inputs to the function are
1. $N_t$ = number of transmit antennas (assuming the same number over all hops)
2. $N_r$ = number of receive antennas (assuming the same number over all hops)
3. A flag called **plotDebugFlag** that plots interim results such as the rate $1/2$ OSTBC of the $(N_t\times N_r)\times 1$ MISO channel for which the transmission-combining scheme is equivalent to.

The function outputs 
1. **transmitMatrixMIMO** which is a symbolic matrix that holds half of the transmission matrix. We note that the transmission matrix is composed of $\rho(N_t\times N_r)$ information symbols where $\rho(n)$ is the Hurwitz-Radon number (of a rate $1/2$ OSTBC for MISO with $N_t\times N_r$ transmit antennas). The full transmission matrix is

$$ transmitMatrix = \begin{bmatrix} 
   transmitMatrixMIMOReal \\
   conj(transmitMatrixMIMOReal)
   \end{bmatrix} $$

2. **combMatReal** is a symbolic matrix holding half of the combining matrix. Entry $y_{i,j}$ means taking the output of antenna $j$ that was received at time $i$. The full combining matrix is composed of stacking combMatReal with an identical matrix in which $y_{i,j^*}=y_{i,j+\rho(N_t\times N_r)}$

Note that the function does not provide the power normalization required to meet any power constraints.
  
# Output Examples
1. N_t=N_r=2

$$ transmitMatrixMIMOReal = \begin{bmatrix}
x_1 &  0 \\
x_2 &  0 \\
x_3 &  0 \\
x_4 &  0 \\
0 &  x_1 \\
0 & x_2 \\
0 & x_3 \\
0 & x_4 \end{bmatrix} $$

$$ combMatReal = \begin{bmatrix} 
y_{1,1} + y_{4,2} + y_{6,1} + y_{7,2} \\
y_{2,1} - y_{3,2} - y_{5,1} + y_{8,2} \\
y_{2,2} + y_{3,1} - y_{5,2} - y_{8,1} \\
y_{4,1} - y_{1,2} - y_{6,2} + y_{7,1} \end{bmatrix} $$

2. N_t=2, N_r=3

$$ transmitMatrixMIMOReal = \begin{bmatrix} 
x_1 & 0 \\
x_2 & 0 \\
x_3 & 0 \\
x_4 & 0 \\
x_5 & 0 \\
x_6 & 0 \\
x_7 & 0 \\
x_8 & 0 \\
0 & x_1 \\
0 & x_2 \\
0 & x_3 \\
0 & x_4 \\
0 & x_5 \\
0 & x_6 \\
0 & x_7 \\
0 & x_8 \end{bmatrix} $$

$$ combMatReal = \begin{bmatrix} 
y_{1,1} + y_{4,3} + y_{6,2} + y_{10,2} + y_{11,3} + y_{15,1} \\
y_{2,1} - y_{3,3} - y_{5,2} - y_{9,2} + y_{12,3} - y_{16,1} \\
y_{2,3} + y_{3,1} - y_{8,2} - y_{9,3} - y_{12,2} + y_{13,1} \\
y_{4,1} - y_{1,3} + y_{7,2} - y_{10,3} + y_{11,2} - y_{14,1} \\
y_{2,2} + y_{5,1} + y_{8,3} - y_{11,1} - y_{14,2} + y_{15,3} \\
y_{6,1} - y_{1,2} - y_{7,3} + y_{12,1} + y_{13,2} + y_{16,3} \\
y_{6,3} - y_{4,2} + y_{7,1} - y_{9,1} - y_{13,3} + y_{16,2} \\
y_{3,2} - y_{5,3} + y_{8,1} + y_{10,1} - y_{14,3} - y_{15,2} 
\end{bmatrix} $$

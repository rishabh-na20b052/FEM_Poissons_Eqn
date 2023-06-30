# FEM_Poissons_Eqn
Rishabh Upadhyay, NA20B052, IITM

I tried to approximate solution  of pdes using FEM in C language.
We will consider the following boundary value problem for the Poisson equation :
![image](https://github.com/rishabh-na20b052/FEM_Poissons_Eqn/assets/87112835/66e4dfcb-976e-401c-aad4-78fff419c6fb)

2D Poisson Equation is given by :
![image](https://github.com/rishabh-na20b052/FEM_Poissons_Eqn/assets/87112835/fd407ff1-b56e-48ad-8f6f-722535978d90)

I have used triangular finite elements for approximation :
![image](https://github.com/rishabh-na20b052/FEM_Poissons_Eqn/assets/87112835/bd5581a7-4e1c-497b-9704-3d857bd2d0d3)

We will approach the problem by forming Global system matrix and right-hand side vector :
![image](https://github.com/rishabh-na20b052/FEM_Poissons_Eqn/assets/87112835/14bf295a-a683-4987-8e7b-14e92c7f699f)

#include <stdio.h>
#include <math.h>

#define N_NODES 8  // Number of nodes
#define N_ELEMENTS 7  // Number of elements

// Function to assemble the global system matrix and right-hand side vector
void assembleSystem(double A[][N_NODES], double b[], double x[], double y[], int elements[][3])
{
    // Element stiffness matrix and load vector
    double Ke[3][3], fe[3];
    
    // Initialize the global matrices and vectors
    for (int i = 0; i < N_NODES; i++)
    {
        b[i] = 0.0;
        for (int j = 0; j < N_NODES; j++)
        {
            A[i][j] = 0.0;
        }
    }
    
    // Assemble the global system
    for (int e = 0; e < N_ELEMENTS; e++)
    {
        // Get the element nodes
        int n1 = elements[e][0];
        int n2 = elements[e][1];
        int n3 = elements[e][2];
        
        // Compute the element Jacobian
        double J = (x[n2] - x[n1]) * (y[n3] - y[n1]) - (y[n2] - y[n1]) * (x[n3] - x[n1]);
        
        // Compute the element inverse Jacobian
        double invJ = 1.0 / J;
        
        // Compute the element stiffness matrix
        Ke[0][0] = invJ * (pow(y[n2] - y[n3], 2) + pow(x[n3] - x[n2], 2));
        Ke[0][1] = invJ * (y[n2] - y[n3]) * (y[n3] - y[n1]) + invJ * (x[n3] - x[n2]) * (x[n1] - x[n3]);
        Ke[0][2] = invJ * (y[n3] - y[n1]) * (y[n1] - y[n2]) + invJ * (x[n1] - x[n3]) * (x[n2] - x[n1]);
        Ke[1][1] = invJ * (pow(y[n3] - y[n1], 2) + pow(x[n1] - x[n3], 2));
        Ke[1][2] = invJ * (y[n1] - y[n2]) * (y[n2] - y[n3]) + invJ * (x[n2] - x[n1]) * (x[n3] - x[n2]);
        Ke[2][2] = invJ * (pow(y[n1] - y[n2], 2) + pow(x[n2] - x[n1], 2));
        
        Ke[1][0] = Ke[0][1];
        Ke[2][0] = Ke[0][2];
        Ke[2][1] = Ke[1][2];
        
        // Compute the element load vector
        fe[0] = J / 6.0;
        fe[1] = J / 6.0;
        fe[2] = J / 6.0;
        
        // Assemble the element contributions into the global system
        A[n1][n1] += Ke[0][0];
        A[n1][n2] += Ke[0][1];
        A[n1][n3] += Ke[0][2];
        A[n2][n1] += Ke[1][0];
        A[n2][n2] += Ke[1][1];
        A[n2][n3] += Ke[1][2];
        A[n3][n1] += Ke[2][0];
        A[n3][n2] += Ke[2][1];
        A[n3][n3] += Ke[2][2];
        
        b[n1] += fe[0];
        b[n2] += fe[1];
        b[n3] += fe[2];
    }
}

// Function to solve the system of equations
void solveSystem(double A[][N_NODES], double b[], double x[])
{
    // Forward elimination (Gaussian elimination without pivoting)
    for (int k = 0; k < N_NODES - 1; k++)
    {
        for (int i = k + 1; i < N_NODES; i++)
        {
            double factor = A[i][k] / A[k][k];
            for (int j = k; j < N_NODES; j++)
            {
                A[i][j] -= factor * A[k][j];
            }
            b[i] -= factor * b[k];
        }
    }
    
    // Back substitution
    x[N_NODES - 1] = b[N_NODES - 1] / A[N_NODES - 1][N_NODES - 1];
    for (int i = N_NODES - 2; i >= 0; i--)
    {
        double sum = 0.0;
        for (int j = i + 1; j < N_NODES; j++)
        {
            sum += A[i][j] * x[j];
        }
        x[i] = (b[i] - sum) / A[i][i];
    }
}

// Function to evaluate the solution at a given point (x, y)
double evaluateSolution(double x[], double y[], double u[], double pointX, double pointY, int elements[][3])
{
    // Find the element that contains the point (x, y)
    int elementIndex = -1;
    for (int e = 0; e < N_ELEMENTS; e++)
    {
        double minX = fmin(fmin(x[elements[e][0]], x[elements[e][1]]), x[elements[e][2]]);
        double maxX = fmax(fmax(x[elements[e][0]], x[elements[e][1]]), x[elements[e][2]]);
        double minY = fmin(fmin(y[elements[e][0]], y[elements[e][1]]), y[elements[e][2]]);
        double maxY = fmax(fmax(y[elements[e][0]], y[elements[e][1]]), y[elements[e][2]]);
        
        if (pointX >= minX && pointX <= maxX && pointY >= minY && pointY <= maxY)
        {
            elementIndex = e;
            break;
        }
    }
    
    if (elementIndex == -1)
    {
        printf("Error: Point (%lf, %lf) is outside the domain.\n", pointX, pointY);
        return 0.0;
    }
    
    // Evaluate the solution at the point (x, y)
    double uPoint = 0.0;
    int n1 = elements[elementIndex][0];
    int n2 = elements[elementIndex][1];
    int n3 = elements[elementIndex][2];
    
    uPoint = u[n1] * (1 - pointX - pointY) + u[n2] * pointX + u[n3] * pointY;
    
    return uPoint;
}

int main()
{
    // Node coordinates
    double x[N_NODES] = {0.0, 0.5, 1.0, 0.0, 0.5, 1.0, 0.0, 0.5};
    double y[N_NODES] = {0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 1.0, 1.0};
    
    // Element-node connectivity
    int elements[N_ELEMENTS][3] = {{0, 1, 3}, {1, 4, 3}, {1, 2, 4}, {2,4,5},{3, 4, 7},{3, 6, 7},{4, 5, 7}};
    
    // Global system matrix and right-hand side vector
    double A[N_NODES][N_NODES];
    double b[N_NODES];
    
    // Solution vector
    double u[N_NODES];
    
    // Assemble the global system
    assembleSystem(A, b, x, y, elements);
    
    // Solve the system
    solveSystem(A, b, u);
    
    // Evaluate the solution at a point (0.75, 0.75) for example
    double pointX = 0.4;
    double pointY = 0.6;
    double uPoint = evaluateSolution(x, y, u, pointX, pointY, elements);
    printf("The solution at point (%lf, %lf) is %lf.\n", pointX, pointY, uPoint);

    return 0;
}

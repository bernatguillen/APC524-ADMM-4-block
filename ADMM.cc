#include <Eigen/Dense>

using namespace Eigen;

int ADMM(const double **C, const double **Aeq, const double *beq, const double **Ain, const double *bin,
         const int n, const int neq, const int nin) {
/* Function to run ADMM.  Parameters:
C   - n by n coefficient matrix
Aeq - neq by n^2 matrix encoding the equality constraints
beq - neq by 1 matrix such that Aeq*vec(X) = beq
Ain - nin by n^2 matrix encoding the inequality constraints
bin - nin by 1 matrix such that Ain*vec(X) <= bin
n   - size of coefficient matrix
neq - number of equality constraints
nin - number of inequality constraints
*/

}

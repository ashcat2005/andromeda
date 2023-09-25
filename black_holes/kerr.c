#include<math.h>


void kerrBH(double* q, double* result, double M, double a) {
    // Auxiliar Functions
    double Sigma = q[1]*q[1] + (a * cos(q[2]))*(a * cos(q[2]));
    double Delta = q[1]*q[1] - 2*M*q[1] + a*a;

    double W = -q[4]*(q[1]*q[1] + a*a) - a*q[7];
    double partXi = q[1]*q[1] + (q[7] + a*q[4])*(q[7] + a*q[4]) + a*a *(1 + q[4]*q[4])*cos(q[2])*cos(q[2]) + (q[7]*cos(q[2])/sin(q[2]))*(q[7]*cos(q[2])/sin(q[2]));
    double Xi = W*W - Delta*partXi;

    double dXidE = 2*W*(q[1]*q[1] + a*a) + 2.*a*Delta*(q[7] + a*q[4]*sin(q[2])*sin(q[2]));
    double dXidL = -2*a*W - 2*a*q[4]*Delta - 2*q[7]*Delta/(sin(q[2])*sin(q[2]));

    double dXidr = -4*q[1]*q[4]*W - 2*(q[1] - M)*partXi - 2*q[1]*Delta;

    double dAdr = (q[1] - M)/Sigma - (q[1]*Delta)/(Sigma*Sigma);
    double dBdr = -q[1]/Sigma/Sigma;
    double dCdr = dXidr/(2*Delta*Sigma) - (Xi*(q[1]-M))/(Sigma*Delta*Delta) - q[1]*Xi/(Delta*Sigma*Sigma);

    double auxth = (a*a) * cos(q[2])*sin(q[2]);

    double dAdth = Delta*auxth/(Sigma*Sigma);
    double dBdth = auxth/(Sigma*Sigma);
    double dCdth = ((1+q[4]*q[4])*auxth + q[7]*q[7] * cos(q[2])/(sin(q[2])*sin(q[2])*sin(q[2])))/Sigma + (Xi/(Delta*Sigma*Sigma))*auxth;

    // Geodesics differential equations 
    result[0] = dXidE/(2.*Delta*Sigma);
    result[1] = (Delta/Sigma)*q[5];
    result[2] = q[6]/Sigma;
    result[3] = - dXidL/(2.*Delta*Sigma);

    result[4] = 0.;
    result[5] = -dAdr*q[5]*q[5] - dBdr*q[6]*q[6] + dCdr;
    result[6] = -dAdth*q[5]*q[5] - dBdth*q[6]*q[6] + dCdth;
    result[7] = 0.;
}

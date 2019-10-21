#ifndef PRACTICA_INVARIANTS_H
#define PRACTICA_INVARIANTS_H

#endif //PRACTICA_INVARIANTS_H

double Eta(double g, double h, double k, double l, double m, double n) {
    return 4*Pow(g,2)*Pow(h,2) - 4*Pow(g,3)*k + 32*Pow(h,3)*l - 36*g*h*k*l - 27*Pow(k,2)*Pow(l,2) -
           16*g*Pow(h,2)*m + 24*Pow(g,2)*k*m + 72*h*k*l*m + 16*Pow(h,2)*Pow(m,2) - 48*g*k*Pow(m,2) +
           32*k*Pow(m,3) - 4*Pow(g,2)*h*n - 48*Pow(h,2)*l*n + 18*g*k*l*n + 16*g*h*m*n - 36*k*l*m*n - 16*h*Pow(m,2)*n +
           Pow(g,2)*Pow(n,2) + 24*h*l*Pow(n,2) - 4*g*m*Pow(n,2) + 4*Pow(m,2)*Pow(n,2) - 4*l*Pow(n,3);
}

double kappa(double h, double l, double g, double k, double m, double n) {
    return 16*(8*Pow(h,3)*l - 8*g*h*k*l - 9*Pow(k,2)*Pow(l,2) - 8*g*Pow(h,2)*m + 8*Pow(g,2)*k*m + 28*h*k*l*m -
               20*g*k*Pow(m,2) + 8*k*Pow(m,3) - 20*Pow(h,2)*l*n + 10*g*k*l*n + 12*g*h*m*n - 8*k*l*m*n - 8*h*Pow(m,2)*n -
               Pow(g,2)*Pow(n,2) + 8*h*l*Pow(n,2));
}

double D(double a, double b, double c, double d, double e, double f, double g, double h, double k, double l, double m, double n) {
    return 48*(-4*Pow(d,2)*Pow(e,2)*Pow(f,4)*Pow(g,2)*Pow(h,2) + 8*c*d*e*Pow(f,5)*Pow(g,2)*Pow(h,2) - 4*Pow(c,2)*Pow(f,6)*Pow(g,2)*Pow(h,2) -
               16*b*d*Pow(f,5)*Pow(g,3)*Pow(h,2) + 16*a*Pow(f,6)*Pow(g,3)*Pow(h,2) + 16*Pow(d,2)*Pow(e,3)*Pow(f,3)*g*Pow(h,3) -
               32*c*d*Pow(e,2)*Pow(f,4)*g*Pow(h,3) + 16*Pow(c,2)*e*Pow(f,5)*g*Pow(h,3) + 80*b*d*e*Pow(f,4)*Pow(g,2)*Pow(h,3) +
               16*b*c*Pow(f,5)*Pow(g,2)*Pow(h,3) - 96*a*e*Pow(f,5)*Pow(g,2)*Pow(h,3) - 16*Pow(d,2)*Pow(e,4)*Pow(f,2)*Pow(h,4) +
               32*c*d*Pow(e,3)*Pow(f,3)*Pow(h,4) - 16*Pow(c,2)*Pow(e,2)*Pow(f,4)*Pow(h,4) - 128*b*d*Pow(e,2)*Pow(f,3)*g*Pow(h,4) - 64*b*c*e*Pow(f,4)*g*Pow(h,4) +
               192*a*Pow(e,2)*Pow(f,4)*g*Pow(h,4) - 16*Pow(b,2)*Pow(f,4)*Pow(g,2)*Pow(h,4) + 64*b*d*Pow(e,3)*Pow(f,2)*Pow(h,5) +
               64*b*c*Pow(e,2)*Pow(f,3)*Pow(h,5) - 128*a*Pow(e,3)*Pow(f,3)*Pow(h,5) + 64*Pow(b,2)*e*Pow(f,3)*g*Pow(h,5) -
               64*Pow(b,2)*Pow(e,2)*Pow(f,2)*Pow(h,6) + 4*Pow(d,2)*Pow(e,2)*Pow(f,4)*Pow(g,3)*k - 8*c*d*e*Pow(f,5)*Pow(g,3)*k +
               4*Pow(c,2)*Pow(f,6)*Pow(g,3)*k + 16*b*d*Pow(f,5)*Pow(g,4)*k - 16*a*Pow(f,6)*Pow(g,4)*k - 16*Pow(d,2)*Pow(e,3)*Pow(f,3)*Pow(g,2)*h*k +
               32*c*d*Pow(e,2)*Pow(f,4)*Pow(g,2)*h*k - 16*Pow(c,2)*e*Pow(f,5)*Pow(g,2)*h*k - 80*b*d*e*Pow(f,4)*Pow(g,3)*h*k - 16*b*c*Pow(f,5)*Pow(g,3)*h*k +
               96*a*e*Pow(f,5)*Pow(g,3)*h*k + 8*Pow(d,2)*Pow(e,4)*Pow(f,2)*g*Pow(h,2)*k - 16*c*d*Pow(e,3)*Pow(f,3)*g*Pow(h,2)*k +
               8*Pow(c,2)*Pow(e,2)*Pow(f,4)*g*Pow(h,2)*k + 96*b*d*Pow(e,2)*Pow(f,3)*Pow(g,2)*Pow(h,2)*k + 48*b*c*e*Pow(f,4)*Pow(g,2)*Pow(h,2)*k -
               144*a*Pow(e,2)*Pow(f,4)*Pow(g,2)*Pow(h,2)*k + 32*Pow(b,2)*Pow(f,4)*Pow(g,3)*Pow(h,2)*k + 16*Pow(d,2)*Pow(e,5)*f*Pow(h,3)*k -
               32*c*d*Pow(e,4)*Pow(f,2)*Pow(h,3)*k + 16*Pow(c,2)*Pow(e,3)*Pow(f,3)*Pow(h,3)*k + 32*b*d*Pow(e,3)*Pow(f,2)*g*Pow(h,3)*k +
               32*b*c*Pow(e,2)*Pow(f,3)*g*Pow(h,3)*k - 64*a*Pow(e,3)*Pow(f,3)*g*Pow(h,3)*k - 128*Pow(b,2)*e*Pow(f,3)*Pow(g,2)*Pow(h,3)*k -
               64*b*d*Pow(e,4)*f*Pow(h,4)*k - 128*b*c*Pow(e,3)*Pow(f,2)*Pow(h,4)*k + 192*a*Pow(e,4)*Pow(f,2)*Pow(h,4)*k +
               96*Pow(b,2)*Pow(e,2)*Pow(f,2)*g*Pow(h,4)*k + 64*Pow(b,2)*Pow(e,3)*f*Pow(h,5)*k + 8*Pow(d,2)*Pow(e,4)*Pow(f,2)*Pow(g,2)*Pow(k,2) -
               16*c*d*Pow(e,3)*Pow(f,3)*Pow(g,2)*Pow(k,2) + 8*Pow(c,2)*Pow(e,2)*Pow(f,4)*Pow(g,2)*Pow(k,2) + 32*b*d*Pow(e,2)*Pow(f,3)*Pow(g,3)*Pow(k,2) +
               16*b*c*e*Pow(f,4)*Pow(g,3)*Pow(k,2) - 48*a*Pow(e,2)*Pow(f,4)*Pow(g,3)*Pow(k,2) - 16*Pow(b,2)*Pow(f,4)*Pow(g,4)*Pow(k,2) -
               16*Pow(d,2)*Pow(e,5)*f*g*h*Pow(k,2) + 32*c*d*Pow(e,4)*Pow(f,2)*g*h*Pow(k,2) - 16*Pow(c,2)*Pow(e,3)*Pow(f,3)*g*h*Pow(k,2) -
               96*b*d*Pow(e,3)*Pow(f,2)*Pow(g,2)*h*Pow(k,2) - 96*b*c*Pow(e,2)*Pow(f,3)*Pow(g,2)*h*Pow(k,2) + 192*a*Pow(e,3)*Pow(f,3)*Pow(g,2)*h*Pow(k,2) +
               64*Pow(b,2)*e*Pow(f,3)*Pow(g,3)*h*Pow(k,2) - 4*Pow(d,2)*Pow(e,6)*Pow(h,2)*Pow(k,2) + 8*c*d*Pow(e,5)*f*Pow(h,2)*Pow(k,2) -
               4*Pow(c,2)*Pow(e,4)*Pow(f,2)*Pow(h,2)*Pow(k,2) + 48*b*d*Pow(e,4)*f*g*Pow(h,2)*Pow(k,2) + 96*b*c*Pow(e,3)*Pow(f,2)*g*Pow(h,2)*Pow(k,2) -
               144*a*Pow(e,4)*Pow(f,2)*g*Pow(h,2)*Pow(k,2) + 16*b*d*Pow(e,5)*Pow(h,3)*Pow(k,2) + 80*b*c*Pow(e,4)*f*Pow(h,3)*Pow(k,2) -
               96*a*Pow(e,5)*f*Pow(h,3)*Pow(k,2) - 128*Pow(b,2)*Pow(e,3)*f*g*Pow(h,3)*Pow(k,2) - 16*Pow(b,2)*Pow(e,4)*Pow(h,4)*Pow(k,2) +
               4*Pow(d,2)*Pow(e,6)*g*Pow(k,3) - 8*c*d*Pow(e,5)*f*g*Pow(k,3) + 4*Pow(c,2)*Pow(e,4)*Pow(f,2)*g*Pow(k,3) + 16*b*d*Pow(e,4)*f*Pow(g,2)*Pow(k,3) +
               32*b*c*Pow(e,3)*Pow(f,2)*Pow(g,2)*Pow(k,3) - 48*a*Pow(e,4)*Pow(f,2)*Pow(g,2)*Pow(k,3) - 32*Pow(b,2)*Pow(e,2)*Pow(f,2)*Pow(g,3)*Pow(k,3) -
               16*b*d*Pow(e,5)*g*h*Pow(k,3) - 80*b*c*Pow(e,4)*f*g*h*Pow(k,3) + 96*a*Pow(e,5)*f*g*h*Pow(k,3) + 64*Pow(b,2)*Pow(e,3)*f*Pow(g,2)*h*Pow(k,3) -
               16*b*c*Pow(e,5)*Pow(h,2)*Pow(k,3) + 16*a*Pow(e,6)*Pow(h,2)*Pow(k,3) + 32*Pow(b,2)*Pow(e,4)*g*Pow(h,2)*Pow(k,3) + 16*b*c*Pow(e,5)*g*Pow(k,4) -
               16*a*Pow(e,6)*g*Pow(k,4) - 16*Pow(b,2)*Pow(e,4)*Pow(g,2)*Pow(k,4) + 8*Pow(d,3)*Pow(e,2)*Pow(f,3)*g*Pow(h,2)*l -
               16*c*Pow(d,2)*e*Pow(f,4)*g*Pow(h,2)*l + 8*Pow(c,2)*d*Pow(f,5)*g*Pow(h,2)*l + 48*b*Pow(d,2)*Pow(f,4)*Pow(g,2)*Pow(h,2)*l -
               48*a*d*Pow(f,5)*Pow(g,2)*Pow(h,2)*l + 16*Pow(d,3)*Pow(e,3)*Pow(f,2)*Pow(h,3)*l - 64*c*Pow(d,2)*Pow(e,2)*Pow(f,3)*Pow(h,3)*l +
               80*Pow(c,2)*d*e*Pow(f,4)*Pow(h,3)*l - 32*Pow(c,3)*Pow(f,5)*Pow(h,3)*l - 16*b*Pow(d,2)*e*Pow(f,3)*g*Pow(h,3)*l - 176*b*c*d*Pow(f,4)*g*Pow(h,3)*l +
               48*a*d*e*Pow(f,4)*g*Pow(h,3)*l + 144*a*c*Pow(f,5)*g*Pow(h,3)*l + 32*b*Pow(d,2)*Pow(e,2)*Pow(f,2)*Pow(h,4)*l - 32*b*c*d*e*Pow(f,3)*Pow(h,4)*l +
               96*a*d*Pow(e,2)*Pow(f,3)*Pow(h,4)*l + 192*b*Pow(c,2)*Pow(f,4)*Pow(h,4)*l - 288*a*c*e*Pow(f,4)*Pow(h,4)*l + 320*Pow(b,2)*d*Pow(f,3)*g*Pow(h,4)*l -
               288*a*b*Pow(f,4)*g*Pow(h,4)*l - 256*Pow(b,2)*d*e*Pow(f,2)*Pow(h,5)*l - 384*Pow(b,2)*c*Pow(f,3)*Pow(h,5)*l + 576*a*b*e*Pow(f,3)*Pow(h,5)*l +
               256*Pow(b,3)*Pow(f,2)*Pow(h,6)*l - 12*Pow(d,3)*Pow(e,2)*Pow(f,3)*Pow(g,2)*k*l + 24*c*Pow(d,2)*e*Pow(f,4)*Pow(g,2)*k*l -
               12*Pow(c,2)*d*Pow(f,5)*Pow(g,2)*k*l - 64*b*Pow(d,2)*Pow(f,4)*Pow(g,3)*k*l + 64*a*d*Pow(f,5)*Pow(g,3)*k*l - 4*Pow(d,3)*Pow(e,3)*Pow(f,2)*g*h*k*l +
               44*c*Pow(d,2)*Pow(e,2)*Pow(f,3)*g*h*k*l - 76*Pow(c,2)*d*e*Pow(f,4)*g*h*k*l + 36*Pow(c,3)*Pow(f,5)*g*h*k*l + 80*b*Pow(d,2)*e*Pow(f,3)*Pow(g,2)*h*k*l +
               208*b*c*d*Pow(f,4)*Pow(g,2)*h*k*l - 128*a*d*e*Pow(f,4)*Pow(g,2)*h*k*l - 160*a*c*Pow(f,5)*Pow(g,2)*h*k*l - 32*Pow(d,3)*Pow(e,4)*f*Pow(h,2)*k*l +
               88*c*Pow(d,2)*Pow(e,3)*Pow(f,2)*Pow(h,2)*k*l - 80*Pow(c,2)*d*Pow(e,2)*Pow(f,3)*Pow(h,2)*k*l + 24*Pow(c,3)*e*Pow(f,4)*Pow(h,2)*k*l -
               216*b*Pow(d,2)*Pow(e,2)*Pow(f,2)*g*Pow(h,2)*k*l + 160*b*c*d*e*Pow(f,3)*g*Pow(h,2)*k*l + 80*a*d*Pow(e,2)*Pow(f,3)*g*Pow(h,2)*k*l -
               232*b*Pow(c,2)*Pow(f,4)*g*Pow(h,2)*k*l + 208*a*c*e*Pow(f,4)*g*Pow(h,2)*k*l - 512*Pow(b,2)*d*Pow(f,3)*Pow(g,2)*Pow(h,2)*k*l +
               416*a*b*Pow(f,4)*Pow(g,2)*Pow(h,2)*k*l + 32*b*Pow(d,2)*Pow(e,3)*f*Pow(h,3)*k*l + 208*b*c*d*Pow(e,2)*Pow(f,2)*Pow(h,3)*k*l -
               304*a*d*Pow(e,3)*Pow(f,2)*Pow(h,3)*k*l - 304*b*Pow(c,2)*e*Pow(f,3)*Pow(h,3)*k*l + 368*a*c*Pow(e,2)*Pow(f,3)*Pow(h,3)*k*l +
               128*Pow(b,2)*d*e*Pow(f,2)*g*Pow(h,3)*k*l + 640*Pow(b,2)*c*Pow(f,3)*g*Pow(h,3)*k*l - 512*a*b*e*Pow(f,3)*g*Pow(h,3)*k*l +
               192*Pow(b,2)*d*Pow(e,2)*f*Pow(h,4)*k*l + 640*Pow(b,2)*c*e*Pow(f,2)*Pow(h,4)*k*l - 928*a*b*Pow(e,2)*Pow(f,2)*Pow(h,4)*k*l -
               640*Pow(b,3)*Pow(f,2)*g*Pow(h,4)*k*l - 256*Pow(b,3)*e*f*Pow(h,5)*k*l + 20*Pow(d,3)*Pow(e,4)*f*g*Pow(k,2)*l -
               76*c*Pow(d,2)*Pow(e,3)*Pow(f,2)*g*Pow(k,2)*l + 92*Pow(c,2)*d*Pow(e,2)*Pow(f,3)*g*Pow(k,2)*l - 36*Pow(c,3)*e*Pow(f,4)*g*Pow(k,2)*l +
               88*b*Pow(d,2)*Pow(e,2)*Pow(f,2)*Pow(g,2)*Pow(k,2)*l - 256*b*c*d*e*Pow(f,3)*Pow(g,2)*Pow(k,2)*l - 16*a*d*Pow(e,2)*Pow(f,3)*Pow(g,2)*Pow(k,2)*l +
               24*b*Pow(c,2)*Pow(f,4)*Pow(g,2)*Pow(k,2)*l + 160*a*c*e*Pow(f,4)*Pow(g,2)*Pow(k,2)*l + 192*Pow(b,2)*d*Pow(f,3)*Pow(g,3)*Pow(k,2)*l -
               128*a*b*Pow(f,4)*Pow(g,3)*Pow(k,2)*l + 4*Pow(d,3)*Pow(e,5)*h*Pow(k,2)*l + 4*c*Pow(d,2)*Pow(e,4)*f*h*Pow(k,2)*l -
               20*Pow(c,2)*d*Pow(e,3)*Pow(f,2)*h*Pow(k,2)*l + 12*Pow(c,3)*Pow(e,2)*Pow(f,3)*h*Pow(k,2)*l + 32*b*Pow(d,2)*Pow(e,3)*f*g*h*Pow(k,2)*l -
               16*b*c*d*Pow(e,2)*Pow(f,2)*g*h*Pow(k,2)*l + 144*a*d*Pow(e,3)*Pow(f,2)*g*h*Pow(k,2)*l + 368*b*Pow(c,2)*e*Pow(f,3)*g*h*Pow(k,2)*l -
               528*a*c*Pow(e,2)*Pow(f,3)*g*h*Pow(k,2)*l + 128*Pow(b,2)*d*e*Pow(f,2)*Pow(g,2)*h*Pow(k,2)*l - 256*Pow(b,2)*c*Pow(f,3)*Pow(g,2)*h*Pow(k,2)*l -
               64*a*b*e*Pow(f,3)*Pow(g,2)*h*Pow(k,2)*l + 24*b*Pow(d,2)*Pow(e,4)*Pow(h,2)*Pow(k,2)*l - 256*b*c*d*Pow(e,3)*f*Pow(h,2)*Pow(k,2)*l +
               160*a*d*Pow(e,4)*f*Pow(h,2)*Pow(k,2)*l + 88*b*Pow(c,2)*Pow(e,2)*Pow(f,2)*Pow(h,2)*Pow(k,2)*l - 16*a*c*Pow(e,3)*Pow(f,2)*Pow(h,2)*Pow(k,2)*l -
               896*Pow(b,2)*c*e*Pow(f,2)*g*Pow(h,2)*Pow(k,2)*l + 896*a*b*Pow(e,2)*Pow(f,2)*g*Pow(h,2)*Pow(k,2)*l +
               512*Pow(b,3)*Pow(f,2)*Pow(g,2)*Pow(h,2)*Pow(k,2)*l - 128*Pow(b,2)*d*Pow(e,3)*Pow(h,3)*Pow(k,2)*l -
               384*Pow(b,2)*c*Pow(e,2)*f*Pow(h,3)*Pow(k,2)*l + 640*a*b*Pow(e,3)*f*Pow(h,3)*Pow(k,2)*l + 512*Pow(b,3)*e*f*g*Pow(h,3)*Pow(k,2)*l +
               128*Pow(b,3)*Pow(e,2)*Pow(h,4)*Pow(k,2)*l - 4*c*Pow(d,2)*Pow(e,5)*Pow(k,3)*l + 8*Pow(c,2)*d*Pow(e,4)*f*Pow(k,3)*l -
               4*Pow(c,3)*Pow(e,3)*Pow(f,2)*Pow(k,3)*l - 40*b*Pow(d,2)*Pow(e,4)*g*Pow(k,3)*l + 128*b*c*d*Pow(e,3)*f*g*Pow(k,3)*l - 80*a*d*Pow(e,4)*f*g*Pow(k,3)*l -
               184*b*Pow(c,2)*Pow(e,2)*Pow(f,2)*g*Pow(k,3)*l + 176*a*c*Pow(e,3)*Pow(f,2)*g*Pow(k,3)*l - 192*Pow(b,2)*d*Pow(e,2)*f*Pow(g,2)*Pow(k,3)*l +
               256*Pow(b,2)*c*e*Pow(f,2)*Pow(g,2)*Pow(k,3)*l + 32*a*b*Pow(e,2)*Pow(f,2)*Pow(g,2)*Pow(k,3)*l - 128*Pow(b,3)*Pow(f,2)*Pow(g,3)*Pow(k,3)*l +
               32*b*c*d*Pow(e,4)*h*Pow(k,3)*l - 16*a*d*Pow(e,5)*h*Pow(k,3)*l + 64*b*Pow(c,2)*Pow(e,3)*f*h*Pow(k,3)*l - 80*a*c*Pow(e,4)*f*h*Pow(k,3)*l +
               128*Pow(b,2)*d*Pow(e,3)*g*h*Pow(k,3)*l + 384*Pow(b,2)*c*Pow(e,2)*f*g*h*Pow(k,3)*l - 640*a*b*Pow(e,3)*f*g*h*Pow(k,3)*l -
               256*Pow(b,3)*e*f*Pow(g,2)*h*Pow(k,3)*l + 128*Pow(b,2)*c*Pow(e,3)*Pow(h,2)*Pow(k,3)*l - 160*a*b*Pow(e,4)*Pow(h,2)*Pow(k,3)*l -
               256*Pow(b,3)*Pow(e,2)*g*Pow(h,2)*Pow(k,3)*l - 16*b*Pow(c,2)*Pow(e,4)*Pow(k,4)*l + 16*a*c*Pow(e,5)*Pow(k,4)*l -
               128*Pow(b,2)*c*Pow(e,3)*g*Pow(k,4)*l + 160*a*b*Pow(e,4)*g*Pow(k,4)*l + 128*Pow(b,3)*Pow(e,2)*Pow(g,2)*Pow(k,4)*l -
               4*Pow(d,4)*Pow(e,2)*Pow(f,2)*Pow(h,2)*Pow(l,2) + 8*c*Pow(d,3)*e*Pow(f,3)*Pow(h,2)*Pow(l,2) - 4*Pow(c,2)*Pow(d,2)*Pow(f,4)*Pow(h,2)*Pow(l,2) -
               48*b*Pow(d,3)*Pow(f,3)*g*Pow(h,2)*Pow(l,2) + 48*a*Pow(d,2)*Pow(f,4)*g*Pow(h,2)*Pow(l,2) - 64*b*Pow(d,3)*e*Pow(f,2)*Pow(h,3)*Pow(l,2) +
               160*b*c*Pow(d,2)*Pow(f,3)*Pow(h,3)*Pow(l,2) + 48*a*Pow(d,2)*e*Pow(f,3)*Pow(h,3)*Pow(l,2) - 144*a*c*d*Pow(f,4)*Pow(h,3)*Pow(l,2) +
               128*Pow(b,2)*Pow(d,2)*Pow(f,2)*Pow(h,4)*Pow(l,2) - 576*a*b*d*Pow(f,3)*Pow(h,4)*Pow(l,2) + 432*Pow(a,2)*Pow(f,4)*Pow(h,4)*Pow(l,2) +
               12*Pow(d,4)*Pow(e,2)*Pow(f,2)*g*k*Pow(l,2) - 24*c*Pow(d,3)*e*Pow(f,3)*g*k*Pow(l,2) + 12*Pow(c,2)*Pow(d,2)*Pow(f,4)*g*k*Pow(l,2) +
               96*b*Pow(d,3)*Pow(f,3)*Pow(g,2)*k*Pow(l,2) - 96*a*Pow(d,2)*Pow(f,4)*Pow(g,2)*k*Pow(l,2) + 20*Pow(d,4)*Pow(e,3)*f*h*k*Pow(l,2) -
               76*c*Pow(d,3)*Pow(e,2)*Pow(f,2)*h*k*Pow(l,2) + 92*Pow(c,2)*Pow(d,2)*e*Pow(f,3)*h*k*Pow(l,2) - 36*Pow(c,3)*d*Pow(f,4)*h*k*Pow(l,2) +
               80*b*Pow(d,3)*e*Pow(f,2)*g*h*k*Pow(l,2) - 368*b*c*Pow(d,2)*Pow(f,3)*g*h*k*Pow(l,2) - 32*a*Pow(d,2)*e*Pow(f,3)*g*h*k*Pow(l,2) +
               320*a*c*d*Pow(f,4)*g*h*k*Pow(l,2) + 96*b*Pow(d,3)*Pow(e,2)*f*Pow(h,2)*k*Pow(l,2) - 160*b*c*Pow(d,2)*e*Pow(f,2)*Pow(h,2)*k*Pow(l,2) +
               88*a*Pow(d,2)*Pow(e,2)*Pow(f,2)*Pow(h,2)*k*Pow(l,2) + 208*b*Pow(c,2)*d*Pow(f,3)*Pow(h,2)*k*Pow(l,2) - 256*a*c*d*e*Pow(f,3)*Pow(h,2)*k*Pow(l,2) +
               24*a*Pow(c,2)*Pow(f,4)*Pow(h,2)*k*Pow(l,2) + 352*Pow(b,2)*Pow(d,2)*Pow(f,2)*g*Pow(h,2)*k*Pow(l,2) + 320*a*b*d*Pow(f,3)*g*Pow(h,2)*k*Pow(l,2) -
               576*Pow(a,2)*Pow(f,4)*g*Pow(h,2)*k*Pow(l,2) - 384*Pow(b,2)*Pow(d,2)*e*f*Pow(h,3)*k*Pow(l,2) - 832*Pow(b,2)*c*d*Pow(f,2)*Pow(h,3)*k*Pow(l,2) +
               1472*a*b*d*e*Pow(f,2)*Pow(h,3)*k*Pow(l,2) + 192*a*b*c*Pow(f,3)*Pow(h,3)*k*Pow(l,2) - 576*Pow(a,2)*e*Pow(f,3)*Pow(h,3)*k*Pow(l,2) +
               256*Pow(b,3)*d*f*Pow(h,4)*k*Pow(l,2) + 384*a*Pow(b,2)*Pow(f,2)*Pow(h,4)*k*Pow(l,2) - Pow(d,4)*Pow(e,4)*Pow(k,2)*Pow(l,2) -
               16*c*Pow(d,3)*Pow(e,3)*f*Pow(k,2)*Pow(l,2) + 62*Pow(c,2)*Pow(d,2)*Pow(e,2)*Pow(f,2)*Pow(k,2)*Pow(l,2) -
               72*Pow(c,3)*d*e*Pow(f,3)*Pow(k,2)*Pow(l,2) + 27*Pow(c,4)*Pow(f,4)*Pow(k,2)*Pow(l,2) - 128*b*Pow(d,3)*Pow(e,2)*f*g*Pow(k,2)*Pow(l,2) +
               176*b*c*Pow(d,2)*e*Pow(f,2)*g*Pow(k,2)*Pow(l,2) + 32*a*Pow(d,2)*Pow(e,2)*Pow(f,2)*g*Pow(k,2)*Pow(l,2) +
               96*b*Pow(c,2)*d*Pow(f,3)*g*Pow(k,2)*Pow(l,2) - 32*a*c*d*e*Pow(f,3)*g*Pow(k,2)*Pow(l,2) - 144*a*Pow(c,2)*Pow(f,4)*g*Pow(k,2)*Pow(l,2) -
               352*Pow(b,2)*Pow(d,2)*Pow(f,2)*Pow(g,2)*Pow(k,2)*Pow(l,2) + 128*a*b*d*Pow(f,3)*Pow(g,2)*Pow(k,2)*Pow(l,2) +
               128*Pow(a,2)*Pow(f,4)*Pow(g,2)*Pow(k,2)*Pow(l,2) - 32*b*Pow(d,3)*Pow(e,3)*h*Pow(k,2)*Pow(l,2) +
               160*b*c*Pow(d,2)*Pow(e,2)*f*h*Pow(k,2)*Pow(l,2) - 96*a*Pow(d,2)*Pow(e,3)*f*h*Pow(k,2)*Pow(l,2) - 176*b*Pow(c,2)*d*e*Pow(f,2)*h*Pow(k,2)*Pow(l,2) +
               48*a*c*d*Pow(e,2)*Pow(f,2)*h*Pow(k,2)*Pow(l,2) - 144*b*Pow(c,3)*Pow(f,3)*h*Pow(k,2)*Pow(l,2) + 240*a*Pow(c,2)*e*Pow(f,3)*h*Pow(k,2)*Pow(l,2) +
               128*Pow(b,2)*Pow(d,2)*e*f*g*h*Pow(k,2)*Pow(l,2) + 576*Pow(b,2)*c*d*Pow(f,2)*g*h*Pow(k,2)*Pow(l,2) - 1088*a*b*d*e*Pow(f,2)*g*h*Pow(k,2)*Pow(l,2) -
               64*a*b*c*Pow(f,3)*g*h*Pow(k,2)*Pow(l,2) + 640*Pow(a,2)*e*Pow(f,3)*g*h*Pow(k,2)*Pow(l,2) + 896*Pow(b,2)*c*d*e*f*Pow(h,2)*Pow(k,2)*Pow(l,2) -
               896*a*b*d*Pow(e,2)*f*Pow(h,2)*Pow(k,2)*Pow(l,2) + 224*Pow(b,2)*Pow(c,2)*Pow(f,2)*Pow(h,2)*Pow(k,2)*Pow(l,2) -
               448*a*b*c*e*Pow(f,2)*Pow(h,2)*Pow(k,2)*Pow(l,2) + 224*Pow(a,2)*Pow(e,2)*Pow(f,2)*Pow(h,2)*Pow(k,2)*Pow(l,2) -
               768*Pow(b,3)*d*f*g*Pow(h,2)*Pow(k,2)*Pow(l,2) - 256*a*Pow(b,2)*Pow(f,2)*g*Pow(h,2)*Pow(k,2)*Pow(l,2) +
               256*Pow(b,3)*d*e*Pow(h,3)*Pow(k,2)*Pow(l,2) + 256*Pow(b,3)*c*f*Pow(h,3)*Pow(k,2)*Pow(l,2) - 1024*a*Pow(b,2)*e*f*Pow(h,3)*Pow(k,2)*Pow(l,2) -
               256*Pow(b,4)*Pow(h,4)*Pow(k,2)*Pow(l,2) + 32*b*c*Pow(d,2)*Pow(e,3)*Pow(k,3)*Pow(l,2) + 8*a*Pow(d,2)*Pow(e,4)*Pow(k,3)*Pow(l,2) -
               128*b*Pow(c,2)*d*Pow(e,2)*f*Pow(k,3)*Pow(l,2) + 64*a*c*d*Pow(e,3)*f*Pow(k,3)*Pow(l,2) + 144*b*Pow(c,3)*e*Pow(f,2)*Pow(k,3)*Pow(l,2) -
               120*a*Pow(c,2)*Pow(e,2)*Pow(f,2)*Pow(k,3)*Pow(l,2) + 128*Pow(b,2)*Pow(d,2)*Pow(e,2)*g*Pow(k,3)*Pow(l,2) -
               384*Pow(b,2)*c*d*e*f*g*Pow(k,3)*Pow(l,2) + 512*a*b*d*Pow(e,2)*f*g*Pow(k,3)*Pow(l,2) - 96*Pow(b,2)*Pow(c,2)*Pow(f,2)*g*Pow(k,3)*Pow(l,2) +
               64*a*b*c*e*Pow(f,2)*g*Pow(k,3)*Pow(l,2) - 320*Pow(a,2)*Pow(e,2)*Pow(f,2)*g*Pow(k,3)*Pow(l,2) + 512*Pow(b,3)*d*f*Pow(g,2)*Pow(k,3)*Pow(l,2) -
               128*a*Pow(b,2)*Pow(f,2)*Pow(g,2)*Pow(k,3)*Pow(l,2) - 256*Pow(b,2)*c*d*Pow(e,2)*h*Pow(k,3)*Pow(l,2) + 128*a*b*d*Pow(e,3)*h*Pow(k,3)*Pow(l,2) -
               256*Pow(b,2)*Pow(c,2)*e*f*h*Pow(k,3)*Pow(l,2) + 384*a*b*c*Pow(e,2)*f*h*Pow(k,3)*Pow(l,2) + 64*Pow(a,2)*Pow(e,3)*f*h*Pow(k,3)*Pow(l,2) -
               256*Pow(b,3)*d*e*g*h*Pow(k,3)*Pow(l,2) - 256*Pow(b,3)*c*f*g*h*Pow(k,3)*Pow(l,2) + 1024*a*Pow(b,2)*e*f*g*h*Pow(k,3)*Pow(l,2) -
               256*Pow(b,3)*c*e*Pow(h,2)*Pow(k,3)*Pow(l,2) + 512*a*Pow(b,2)*Pow(e,2)*Pow(h,2)*Pow(k,3)*Pow(l,2) + 512*Pow(b,4)*g*Pow(h,2)*Pow(k,3)*Pow(l,2) +
               128*Pow(b,2)*Pow(c,2)*Pow(e,2)*Pow(k,4)*Pow(l,2) - 128*a*b*c*Pow(e,3)*Pow(k,4)*Pow(l,2) - 16*Pow(a,2)*Pow(e,4)*Pow(k,4)*Pow(l,2) +
               256*Pow(b,3)*c*e*g*Pow(k,4)*Pow(l,2) - 512*a*Pow(b,2)*Pow(e,2)*g*Pow(k,4)*Pow(l,2) - 256*Pow(b,4)*Pow(g,2)*Pow(k,4)*Pow(l,2) +
               16*b*Pow(d,4)*Pow(f,2)*Pow(h,2)*Pow(l,3) - 16*a*Pow(d,3)*Pow(f,3)*Pow(h,2)*Pow(l,3) - 4*Pow(d,5)*Pow(e,2)*f*k*Pow(l,3) +
               8*c*Pow(d,4)*e*Pow(f,2)*k*Pow(l,3) - 4*Pow(c,2)*Pow(d,3)*Pow(f,3)*k*Pow(l,3) - 64*b*Pow(d,4)*Pow(f,2)*g*k*Pow(l,3) +
               64*a*Pow(d,3)*Pow(f,3)*g*k*Pow(l,3) - 80*b*Pow(d,4)*e*f*h*k*Pow(l,3) + 176*b*c*Pow(d,3)*Pow(f,2)*h*k*Pow(l,3) +
               64*a*Pow(d,3)*e*Pow(f,2)*h*k*Pow(l,3) - 160*a*c*Pow(d,2)*Pow(f,3)*h*k*Pow(l,3) + 128*Pow(b,2)*Pow(d,3)*f*Pow(h,2)*k*Pow(l,3) -
               736*a*b*Pow(d,2)*Pow(f,2)*Pow(h,2)*k*Pow(l,3) + 576*Pow(a,2)*d*Pow(f,3)*Pow(h,2)*k*Pow(l,3) + 8*b*Pow(d,4)*Pow(e,2)*Pow(k,2)*Pow(l,3) +
               64*b*c*Pow(d,3)*e*f*Pow(k,2)*Pow(l,3) + 32*a*Pow(d,3)*Pow(e,2)*f*Pow(k,2)*Pow(l,3) - 120*b*Pow(c,2)*Pow(d,2)*Pow(f,2)*Pow(k,2)*Pow(l,3) -
               128*a*c*Pow(d,2)*e*Pow(f,2)*Pow(k,2)*Pow(l,3) + 144*a*Pow(c,2)*d*Pow(f,3)*Pow(k,2)*Pow(l,3) + 192*Pow(b,2)*Pow(d,3)*f*g*Pow(k,2)*Pow(l,3) +
               128*a*b*Pow(d,2)*Pow(f,2)*g*Pow(k,2)*Pow(l,3) - 256*Pow(a,2)*d*Pow(f,3)*g*Pow(k,2)*Pow(l,3) + 64*Pow(b,2)*Pow(d,3)*e*h*Pow(k,2)*Pow(l,3) -
               704*Pow(b,2)*c*Pow(d,2)*f*h*Pow(k,2)*Pow(l,3) + 384*a*b*Pow(d,2)*e*f*h*Pow(k,2)*Pow(l,3) + 832*a*b*c*d*Pow(f,2)*h*Pow(k,2)*Pow(l,3) -
               256*Pow(a,2)*d*e*Pow(f,2)*h*Pow(k,2)*Pow(l,3) - 384*Pow(a,2)*c*Pow(f,3)*h*Pow(k,2)*Pow(l,3) - 128*Pow(b,3)*Pow(d,2)*Pow(h,2)*Pow(k,2)*Pow(l,3) +
               1024*a*Pow(b,2)*d*f*Pow(h,2)*Pow(k,2)*Pow(l,3) - 384*Pow(a,2)*b*Pow(f,2)*Pow(h,2)*Pow(k,2)*Pow(l,3) -
               64*Pow(b,2)*c*Pow(d,2)*e*Pow(k,3)*Pow(l,3) - 64*a*b*Pow(d,2)*Pow(e,2)*Pow(k,3)*Pow(l,3) + 384*Pow(b,2)*Pow(c,2)*d*f*Pow(k,3)*Pow(l,3) -
               256*a*b*c*d*e*f*Pow(k,3)*Pow(l,3) - 64*Pow(a,2)*d*Pow(e,2)*f*Pow(k,3)*Pow(l,3) - 288*a*b*Pow(c,2)*Pow(f,2)*Pow(k,3)*Pow(l,3) +
               384*Pow(a,2)*c*e*Pow(f,2)*Pow(k,3)*Pow(l,3) - 128*Pow(b,3)*Pow(d,2)*g*Pow(k,3)*Pow(l,3) - 768*a*Pow(b,2)*d*f*g*Pow(k,3)*Pow(l,3) +
               512*Pow(a,2)*b*Pow(f,2)*g*Pow(k,3)*Pow(l,3) + 512*Pow(b,3)*c*d*h*Pow(k,3)*Pow(l,3) - 256*a*Pow(b,2)*d*e*h*Pow(k,3)*Pow(l,3) -
               256*a*Pow(b,2)*c*f*h*Pow(k,3)*Pow(l,3) - 256*Pow(a,2)*b*e*f*h*Pow(k,3)*Pow(l,3) - 512*a*Pow(b,3)*Pow(h,2)*Pow(k,3)*Pow(l,3) -
               256*Pow(b,3)*Pow(c,2)*Pow(k,4)*Pow(l,3) + 256*a*Pow(b,2)*c*e*Pow(k,4)*Pow(l,3) + 128*Pow(a,2)*b*Pow(e,2)*Pow(k,4)*Pow(l,3) +
               512*a*Pow(b,3)*g*Pow(k,4)*Pow(l,3) + 16*b*Pow(d,5)*f*k*Pow(l,4) - 16*a*Pow(d,4)*Pow(f,2)*k*Pow(l,4) - 16*Pow(b,2)*Pow(d,4)*Pow(k,2)*Pow(l,4) -
               128*a*b*Pow(d,3)*f*Pow(k,2)*Pow(l,4) + 128*Pow(a,2)*Pow(d,2)*Pow(f,2)*Pow(k,2)*Pow(l,4) + 128*a*Pow(b,2)*Pow(d,2)*Pow(k,3)*Pow(l,4) +
               256*Pow(a,2)*b*d*f*Pow(k,3)*Pow(l,4) - 256*Pow(a,3)*Pow(f,2)*Pow(k,3)*Pow(l,4) - 256*Pow(a,2)*Pow(b,2)*Pow(k,4)*Pow(l,4) +
               8*Pow(d,3)*Pow(e,2)*Pow(f,3)*Pow(g,2)*h*m - 16*c*Pow(d,2)*e*Pow(f,4)*Pow(g,2)*h*m + 8*Pow(c,2)*d*Pow(f,5)*Pow(g,2)*h*m +
               32*b*Pow(d,2)*Pow(f,4)*Pow(g,3)*h*m - 32*a*d*Pow(f,5)*Pow(g,3)*h*m - 64*Pow(d,3)*Pow(e,3)*Pow(f,2)*g*Pow(h,2)*m +
               144*c*Pow(d,2)*Pow(e,2)*Pow(f,3)*g*Pow(h,2)*m - 96*Pow(c,2)*d*e*Pow(f,4)*g*Pow(h,2)*m + 16*Pow(c,3)*Pow(f,5)*g*Pow(h,2)*m -
               304*b*Pow(d,2)*e*Pow(f,3)*Pow(g,2)*Pow(h,2)*m + 16*b*c*d*Pow(f,4)*Pow(g,2)*Pow(h,2)*m + 352*a*d*e*Pow(f,4)*Pow(g,2)*Pow(h,2)*m -
               64*a*c*Pow(f,5)*Pow(g,2)*Pow(h,2)*m + 32*Pow(d,3)*Pow(e,4)*f*Pow(h,3)*m - 32*c*Pow(d,2)*Pow(e,3)*Pow(f,2)*Pow(h,3)*m -
               32*Pow(c,2)*d*Pow(e,2)*Pow(f,3)*Pow(h,3)*m + 32*Pow(c,3)*e*Pow(f,4)*Pow(h,3)*m + 352*b*Pow(d,2)*Pow(e,2)*Pow(f,2)*g*Pow(h,3)*m +
               544*b*c*d*e*Pow(f,3)*g*Pow(h,3)*m - 736*a*d*Pow(e,2)*Pow(f,3)*g*Pow(h,3)*m - 128*b*Pow(c,2)*Pow(f,4)*g*Pow(h,3)*m - 32*a*c*e*Pow(f,4)*g*Pow(h,3)*m -
               256*Pow(b,2)*d*Pow(f,3)*Pow(g,2)*Pow(h,3)*m + 320*a*b*Pow(f,4)*Pow(g,2)*Pow(h,3)*m - 128*b*Pow(d,2)*Pow(e,3)*f*Pow(h,4)*m -
               384*b*c*d*Pow(e,2)*Pow(f,2)*Pow(h,4)*m + 320*a*d*Pow(e,3)*Pow(f,2)*Pow(h,4)*m - 128*b*Pow(c,2)*e*Pow(f,3)*Pow(h,4)*m +
               320*a*c*Pow(e,2)*Pow(f,3)*Pow(h,4)*m + 64*Pow(b,2)*d*e*Pow(f,2)*g*Pow(h,4)*m + 320*Pow(b,2)*c*Pow(f,3)*g*Pow(h,4)*m -
               704*a*b*e*Pow(f,3)*g*Pow(h,4)*m + 128*Pow(b,2)*d*Pow(e,2)*f*Pow(h,5)*m + 128*Pow(b,2)*c*e*Pow(f,2)*Pow(h,5)*m +
               128*a*b*Pow(e,2)*Pow(f,2)*Pow(h,5)*m - 256*Pow(b,3)*Pow(f,2)*g*Pow(h,5)*m + 40*Pow(d,3)*Pow(e,3)*Pow(f,2)*Pow(g,2)*k*m -
               104*c*Pow(d,2)*Pow(e,2)*Pow(f,3)*Pow(g,2)*k*m + 88*Pow(c,2)*d*e*Pow(f,4)*Pow(g,2)*k*m - 24*Pow(c,3)*Pow(f,5)*Pow(g,2)*k*m +
               176*b*Pow(d,2)*e*Pow(f,3)*Pow(g,3)*k*m - 80*b*c*d*Pow(f,4)*Pow(g,3)*k*m - 192*a*d*e*Pow(f,4)*Pow(g,3)*k*m + 96*a*c*Pow(f,5)*Pow(g,3)*k*m -
               8*Pow(d,3)*Pow(e,4)*f*g*h*k*m + 8*c*Pow(d,2)*Pow(e,3)*Pow(f,2)*g*h*k*m + 8*Pow(c,2)*d*Pow(e,2)*Pow(f,3)*g*h*k*m - 8*Pow(c,3)*e*Pow(f,4)*g*h*k*m -
               160*b*Pow(d,2)*Pow(e,2)*Pow(f,2)*Pow(g,2)*h*k*m - 288*b*c*d*e*Pow(f,3)*Pow(g,2)*h*k*m + 416*a*d*Pow(e,2)*Pow(f,3)*Pow(g,2)*h*k*m +
               160*b*Pow(c,2)*Pow(f,4)*Pow(g,2)*h*k*m - 128*a*c*e*Pow(f,4)*Pow(g,2)*h*k*m + 256*Pow(b,2)*d*Pow(f,3)*Pow(g,3)*h*k*m -
               320*a*b*Pow(f,4)*Pow(g,3)*h*k*m + 16*Pow(d,3)*Pow(e,5)*Pow(h,2)*k*m - 96*c*Pow(d,2)*Pow(e,4)*f*Pow(h,2)*k*m +
               144*Pow(c,2)*d*Pow(e,3)*Pow(f,2)*Pow(h,2)*k*m - 64*Pow(c,3)*Pow(e,2)*Pow(f,3)*Pow(h,2)*k*m + 208*b*Pow(d,2)*Pow(e,3)*f*g*Pow(h,2)*k*m -
               608*b*c*d*Pow(e,2)*Pow(f,2)*g*Pow(h,2)*k*m + 96*a*d*Pow(e,3)*Pow(f,2)*g*Pow(h,2)*k*m + 208*b*Pow(c,2)*e*Pow(f,3)*g*Pow(h,2)*k*m +
               96*a*c*Pow(e,2)*Pow(f,3)*g*Pow(h,2)*k*m + 640*Pow(b,2)*d*e*Pow(f,2)*Pow(g,2)*Pow(h,2)*k*m - 384*Pow(b,2)*c*Pow(f,3)*Pow(g,2)*Pow(h,2)*k*m +
               128*a*b*e*Pow(f,3)*Pow(g,2)*Pow(h,2)*k*m - 128*b*Pow(d,2)*Pow(e,4)*Pow(h,3)*k*m + 544*b*c*d*Pow(e,3)*f*Pow(h,3)*k*m - 32*a*d*Pow(e,4)*f*Pow(h,3)*k*m +
               352*b*Pow(c,2)*Pow(e,2)*Pow(f,2)*Pow(h,3)*k*m - 736*a*c*Pow(e,3)*Pow(f,2)*Pow(h,3)*k*m - 1024*Pow(b,2)*d*Pow(e,2)*f*g*Pow(h,3)*k*m -
               1024*Pow(b,2)*c*e*Pow(f,2)*g*Pow(h,3)*k*m + 1664*a*b*Pow(e,2)*Pow(f,2)*g*Pow(h,3)*k*m + 512*Pow(b,3)*Pow(f,2)*Pow(g,2)*Pow(h,3)*k*m +
               320*Pow(b,2)*d*Pow(e,3)*Pow(h,4)*k*m + 64*Pow(b,2)*c*Pow(e,2)*f*Pow(h,4)*k*m - 704*a*b*Pow(e,3)*f*Pow(h,4)*k*m + 512*Pow(b,3)*e*f*g*Pow(h,4)*k*m -
               256*Pow(b,3)*Pow(e,2)*Pow(h,5)*k*m - 24*Pow(d,3)*Pow(e,5)*g*Pow(k,2)*m + 88*c*Pow(d,2)*Pow(e,4)*f*g*Pow(k,2)*m -
               104*Pow(c,2)*d*Pow(e,3)*Pow(f,2)*g*Pow(k,2)*m + 40*Pow(c,3)*Pow(e,2)*Pow(f,3)*g*Pow(k,2)*m - 208*b*Pow(d,2)*Pow(e,3)*f*Pow(g,2)*Pow(k,2)*m +
               608*b*c*d*Pow(e,2)*Pow(f,2)*Pow(g,2)*Pow(k,2)*m - 96*a*d*Pow(e,3)*Pow(f,2)*Pow(g,2)*Pow(k,2)*m - 208*b*Pow(c,2)*e*Pow(f,3)*Pow(g,2)*Pow(k,2)*m -
               96*a*c*Pow(e,2)*Pow(f,3)*Pow(g,2)*Pow(k,2)*m - 704*Pow(b,2)*d*e*Pow(f,2)*Pow(g,3)*Pow(k,2)*m + 64*Pow(b,2)*c*Pow(f,3)*Pow(g,3)*Pow(k,2)*m +
               576*a*b*e*Pow(f,3)*Pow(g,3)*Pow(k,2)*m + 8*c*Pow(d,2)*Pow(e,5)*h*Pow(k,2)*m - 16*Pow(c,2)*d*Pow(e,4)*f*h*Pow(k,2)*m +
               8*Pow(c,3)*Pow(e,3)*Pow(f,2)*h*Pow(k,2)*m + 160*b*Pow(d,2)*Pow(e,4)*g*h*Pow(k,2)*m - 288*b*c*d*Pow(e,3)*f*g*h*Pow(k,2)*m -
               128*a*d*Pow(e,4)*f*g*h*Pow(k,2)*m - 160*b*Pow(c,2)*Pow(e,2)*Pow(f,2)*g*h*Pow(k,2)*m + 416*a*c*Pow(e,3)*Pow(f,2)*g*h*Pow(k,2)*m +
               896*Pow(b,2)*d*Pow(e,2)*f*Pow(g,2)*h*Pow(k,2)*m + 896*Pow(b,2)*c*e*Pow(f,2)*Pow(g,2)*h*Pow(k,2)*m -
               1792*a*b*Pow(e,2)*Pow(f,2)*Pow(g,2)*h*Pow(k,2)*m - 256*Pow(b,3)*Pow(f,2)*Pow(g,3)*h*Pow(k,2)*m + 16*b*c*d*Pow(e,4)*Pow(h,2)*Pow(k,2)*m -
               64*a*d*Pow(e,5)*Pow(h,2)*Pow(k,2)*m - 304*b*Pow(c,2)*Pow(e,3)*f*Pow(h,2)*Pow(k,2)*m + 352*a*c*Pow(e,4)*f*Pow(h,2)*Pow(k,2)*m -
               384*Pow(b,2)*d*Pow(e,3)*g*Pow(h,2)*Pow(k,2)*m + 640*Pow(b,2)*c*Pow(e,2)*f*g*Pow(h,2)*Pow(k,2)*m + 128*a*b*Pow(e,3)*f*g*Pow(h,2)*Pow(k,2)*m -
               1024*Pow(b,3)*e*f*Pow(g,2)*Pow(h,2)*Pow(k,2)*m - 256*Pow(b,2)*c*Pow(e,3)*Pow(h,3)*Pow(k,2)*m + 320*a*b*Pow(e,4)*Pow(h,3)*Pow(k,2)*m +
               512*Pow(b,3)*Pow(e,2)*g*Pow(h,3)*Pow(k,2)*m - 80*b*c*d*Pow(e,4)*g*Pow(k,3)*m + 96*a*d*Pow(e,5)*g*Pow(k,3)*m +
               176*b*Pow(c,2)*Pow(e,3)*f*g*Pow(k,3)*m - 192*a*c*Pow(e,4)*f*g*Pow(k,3)*m + 64*Pow(b,2)*d*Pow(e,3)*Pow(g,2)*Pow(k,3)*m -
               704*Pow(b,2)*c*Pow(e,2)*f*Pow(g,2)*Pow(k,3)*m + 576*a*b*Pow(e,3)*f*Pow(g,2)*Pow(k,3)*m + 512*Pow(b,3)*e*f*Pow(g,3)*Pow(k,3)*m +
               32*b*Pow(c,2)*Pow(e,4)*h*Pow(k,3)*m - 32*a*c*Pow(e,5)*h*Pow(k,3)*m + 256*Pow(b,2)*c*Pow(e,3)*g*h*Pow(k,3)*m - 320*a*b*Pow(e,4)*g*h*Pow(k,3)*m -
               256*Pow(b,3)*Pow(e,2)*Pow(g,2)*h*Pow(k,3)*m - 16*Pow(d,4)*Pow(e,2)*Pow(f,2)*g*h*l*m + 32*c*Pow(d,3)*e*Pow(f,3)*g*h*l*m -
               16*Pow(c,2)*Pow(d,2)*Pow(f,4)*g*h*l*m - 96*b*Pow(d,3)*Pow(f,3)*Pow(g,2)*h*l*m + 96*a*Pow(d,2)*Pow(f,4)*Pow(g,2)*h*l*m -
               32*Pow(d,4)*Pow(e,3)*f*Pow(h,2)*l*m + 144*c*Pow(d,3)*Pow(e,2)*Pow(f,2)*Pow(h,2)*l*m - 192*Pow(c,2)*Pow(d,2)*e*Pow(f,3)*Pow(h,2)*l*m +
               80*Pow(c,3)*d*Pow(f,4)*Pow(h,2)*l*m + 176*b*Pow(d,3)*e*Pow(f,2)*g*Pow(h,2)*l*m + 400*b*c*Pow(d,2)*Pow(f,3)*g*Pow(h,2)*l*m -
               272*a*Pow(d,2)*e*Pow(f,3)*g*Pow(h,2)*l*m - 304*a*c*d*Pow(f,4)*g*Pow(h,2)*l*m - 64*b*Pow(d,3)*Pow(e,2)*f*Pow(h,3)*l*m +
               32*b*c*Pow(d,2)*e*Pow(f,2)*Pow(h,3)*l*m - 320*a*Pow(d,2)*Pow(e,2)*Pow(f,2)*Pow(h,3)*l*m - 736*b*Pow(c,2)*d*Pow(f,3)*Pow(h,3)*l*m +
               992*a*c*d*e*Pow(f,3)*Pow(h,3)*l*m + 96*a*Pow(c,2)*Pow(f,4)*Pow(h,3)*l*m - 1216*Pow(b,2)*Pow(d,2)*Pow(f,2)*g*Pow(h,3)*l*m +
               1664*a*b*d*Pow(f,3)*g*Pow(h,3)*l*m - 576*Pow(a,2)*Pow(f,4)*g*Pow(h,3)*l*m + 512*Pow(b,2)*Pow(d,2)*e*f*Pow(h,4)*l*m +
               1408*Pow(b,2)*c*d*Pow(f,2)*Pow(h,4)*l*m - 1216*a*b*d*e*Pow(f,2)*Pow(h,4)*l*m + 192*a*b*c*Pow(f,3)*Pow(h,4)*l*m -
               576*Pow(a,2)*e*Pow(f,3)*Pow(h,4)*l*m - 512*Pow(b,3)*d*f*Pow(h,5)*l*m - 768*a*Pow(b,2)*Pow(f,2)*Pow(h,5)*l*m - 44*Pow(d,4)*Pow(e,3)*f*g*k*l*m +
               100*c*Pow(d,3)*Pow(e,2)*Pow(f,2)*g*k*l*m - 68*Pow(c,2)*Pow(d,2)*e*Pow(f,3)*g*k*l*m + 12*Pow(c,3)*d*Pow(f,4)*g*k*l*m -
               368*b*Pow(d,3)*e*Pow(f,2)*Pow(g,2)*k*l*m + 80*b*c*Pow(d,2)*Pow(f,3)*Pow(g,2)*k*l*m + 416*a*Pow(d,2)*e*Pow(f,3)*Pow(g,2)*k*l*m -
               128*a*c*d*Pow(f,4)*Pow(g,2)*k*l*m - 16*Pow(d,4)*Pow(e,4)*h*k*l*m + 136*c*Pow(d,3)*Pow(e,3)*f*h*k*l*m -
               296*Pow(c,2)*Pow(d,2)*Pow(e,2)*Pow(f,2)*h*k*l*m + 248*Pow(c,3)*d*e*Pow(f,3)*h*k*l*m - 72*Pow(c,4)*Pow(f,4)*h*k*l*m +
               160*b*Pow(d,3)*Pow(e,2)*f*g*h*k*l*m + 480*b*c*Pow(d,2)*e*Pow(f,2)*g*h*k*l*m - 208*a*Pow(d,2)*Pow(e,2)*Pow(f,2)*g*h*k*l*m -
               64*b*Pow(c,2)*d*Pow(f,3)*g*h*k*l*m - 576*a*c*d*e*Pow(f,3)*g*h*k*l*m + 208*a*Pow(c,2)*Pow(f,4)*g*h*k*l*m +
               704*Pow(b,2)*Pow(d,2)*Pow(f,2)*Pow(g,2)*h*k*l*m - 1152*a*b*d*Pow(f,3)*Pow(g,2)*h*k*l*m + 640*Pow(a,2)*Pow(f,4)*Pow(g,2)*h*k*l*m +
               32*b*Pow(d,3)*Pow(e,3)*Pow(h,2)*k*l*m - 608*b*c*Pow(d,2)*Pow(e,2)*f*Pow(h,2)*k*l*m + 208*a*Pow(d,2)*Pow(e,3)*f*Pow(h,2)*k*l*m +
               400*b*Pow(c,2)*d*e*Pow(f,2)*Pow(h,2)*k*l*m + 400*a*c*d*Pow(e,2)*Pow(f,2)*Pow(h,2)*k*l*m + 368*b*Pow(c,3)*Pow(f,3)*Pow(h,2)*k*l*m -
               800*a*Pow(c,2)*e*Pow(f,3)*Pow(h,2)*k*l*m - 64*Pow(b,2)*Pow(d,2)*e*f*g*Pow(h,2)*k*l*m - 512*Pow(b,2)*c*d*Pow(f,2)*g*Pow(h,2)*k*l*m -
               1024*a*b*d*e*Pow(f,2)*g*Pow(h,2)*k*l*m - 640*a*b*c*Pow(f,3)*g*Pow(h,2)*k*l*m + 1472*Pow(a,2)*e*Pow(f,3)*g*Pow(h,2)*k*l*m +
               384*Pow(b,2)*Pow(d,2)*Pow(e,2)*Pow(h,3)*k*l*m - 1152*Pow(b,2)*c*d*e*f*Pow(h,3)*k*l*m + 640*a*b*d*Pow(e,2)*f*Pow(h,3)*k*l*m -
               64*Pow(b,2)*Pow(c,2)*Pow(f,2)*Pow(h,3)*k*l*m - 256*a*b*c*e*Pow(f,2)*Pow(h,3)*k*l*m + 832*Pow(a,2)*Pow(e,2)*Pow(f,2)*Pow(h,3)*k*l*m +
               2048*Pow(b,3)*d*f*g*Pow(h,3)*k*l*m - 512*a*Pow(b,2)*Pow(f,2)*g*Pow(h,3)*k*l*m - 1280*Pow(b,3)*d*e*Pow(h,4)*k*l*m - 1280*Pow(b,3)*c*f*Pow(h,4)*k*l*m +
               3328*a*Pow(b,2)*e*f*Pow(h,4)*k*l*m + 1024*Pow(b,4)*Pow(h,5)*k*l*m + 20*c*Pow(d,3)*Pow(e,4)*Pow(k,2)*l*m -
               76*Pow(c,2)*Pow(d,2)*Pow(e,3)*f*Pow(k,2)*l*m + 92*Pow(c,3)*d*Pow(e,2)*Pow(f,2)*Pow(k,2)*l*m - 36*Pow(c,4)*e*Pow(f,3)*Pow(k,2)*l*m +
               160*b*Pow(d,3)*Pow(e,3)*g*Pow(k,2)*l*m - 128*b*c*Pow(d,2)*Pow(e,2)*f*g*Pow(k,2)*l*m + 32*a*Pow(d,2)*Pow(e,3)*f*g*Pow(k,2)*l*m -
               464*b*Pow(c,2)*d*e*Pow(f,2)*g*Pow(k,2)*l*m - 16*a*c*d*Pow(e,2)*Pow(f,2)*g*Pow(k,2)*l*m + 48*b*Pow(c,3)*Pow(f,3)*g*Pow(k,2)*l*m +
               368*a*Pow(c,2)*e*Pow(f,3)*g*Pow(k,2)*l*m + 576*Pow(b,2)*Pow(d,2)*e*f*Pow(g,2)*Pow(k,2)*l*m + 128*Pow(b,2)*c*d*Pow(f,2)*Pow(g,2)*Pow(k,2)*l*m +
               704*a*b*d*e*Pow(f,2)*Pow(g,2)*Pow(k,2)*l*m - 64*a*b*c*Pow(f,3)*Pow(g,2)*Pow(k,2)*l*m - 1152*Pow(a,2)*e*Pow(f,3)*Pow(g,2)*Pow(k,2)*l*m -
               256*b*c*Pow(d,2)*Pow(e,3)*h*Pow(k,2)*l*m + 48*a*Pow(d,2)*Pow(e,4)*h*Pow(k,2)*l*m + 800*b*Pow(c,2)*d*Pow(e,2)*f*h*Pow(k,2)*l*m -
               288*a*c*d*Pow(e,3)*f*h*Pow(k,2)*l*m - 256*b*Pow(c,3)*e*Pow(f,2)*h*Pow(k,2)*l*m - 48*a*Pow(c,2)*Pow(e,2)*Pow(f,2)*h*Pow(k,2)*l*m -
               896*Pow(b,2)*Pow(d,2)*Pow(e,2)*g*h*Pow(k,2)*l*m - 896*Pow(b,2)*c*d*e*f*g*h*Pow(k,2)*l*m + 896*a*b*d*Pow(e,2)*f*g*h*Pow(k,2)*l*m -
               448*Pow(b,2)*Pow(c,2)*Pow(f,2)*g*h*Pow(k,2)*l*m + 1792*a*b*c*e*Pow(f,2)*g*h*Pow(k,2)*l*m - 448*Pow(a,2)*Pow(e,2)*Pow(f,2)*g*h*Pow(k,2)*l*m -
               1536*Pow(b,3)*d*f*Pow(g,2)*h*Pow(k,2)*l*m + 1280*a*Pow(b,2)*Pow(f,2)*Pow(g,2)*h*Pow(k,2)*l*m + 640*Pow(b,2)*c*d*Pow(e,2)*Pow(h,2)*Pow(k,2)*l*m +
               128*a*b*d*Pow(e,3)*Pow(h,2)*Pow(k,2)*l*m + 192*Pow(b,2)*Pow(c,2)*e*f*Pow(h,2)*Pow(k,2)*l*m - 512*a*b*c*Pow(e,2)*f*Pow(h,2)*Pow(k,2)*l*m -
               832*Pow(a,2)*Pow(e,3)*f*Pow(h,2)*Pow(k,2)*l*m + 1536*Pow(b,3)*d*e*g*Pow(h,2)*Pow(k,2)*l*m + 1536*Pow(b,3)*c*f*g*Pow(h,2)*Pow(k,2)*l*m -
               2560*a*Pow(b,2)*e*f*g*Pow(h,2)*Pow(k,2)*l*m + 1024*Pow(b,3)*c*e*Pow(h,3)*Pow(k,2)*l*m - 2048*a*Pow(b,2)*Pow(e,2)*Pow(h,3)*Pow(k,2)*l*m -
               2048*Pow(b,4)*g*Pow(h,3)*Pow(k,2)*l*m + 64*b*Pow(c,2)*d*Pow(e,3)*Pow(k,3)*l*m - 80*a*c*d*Pow(e,4)*Pow(k,3)*l*m -
               160*b*Pow(c,3)*Pow(e,2)*f*Pow(k,3)*l*m + 176*a*Pow(c,2)*Pow(e,3)*f*Pow(k,3)*l*m + 384*Pow(b,2)*c*d*Pow(e,2)*g*Pow(k,3)*l*m -
               640*a*b*d*Pow(e,3)*g*Pow(k,3)*l*m + 832*Pow(b,2)*Pow(c,2)*e*f*g*Pow(k,3)*l*m - 1024*a*b*c*Pow(e,2)*f*g*Pow(k,3)*l*m +
               576*Pow(a,2)*Pow(e,3)*f*g*Pow(k,3)*l*m - 256*Pow(b,3)*d*e*Pow(g,2)*Pow(k,3)*l*m - 256*Pow(b,3)*c*f*Pow(g,2)*Pow(k,3)*l*m -
               768*a*Pow(b,2)*e*f*Pow(g,2)*Pow(k,3)*l*m - 512*Pow(b,2)*Pow(c,2)*Pow(e,2)*h*Pow(k,3)*l*m + 512*a*b*c*Pow(e,3)*h*Pow(k,3)*l*m +
               64*Pow(a,2)*Pow(e,4)*h*Pow(k,3)*l*m - 1024*Pow(b,3)*c*e*g*h*Pow(k,3)*l*m + 2048*a*Pow(b,2)*Pow(e,2)*g*h*Pow(k,3)*l*m +
               1024*Pow(b,4)*Pow(g,2)*h*Pow(k,3)*l*m + 8*Pow(d,5)*Pow(e,2)*f*h*Pow(l,2)*m - 16*c*Pow(d,4)*e*Pow(f,2)*h*Pow(l,2)*m +
               8*Pow(c,2)*Pow(d,3)*Pow(f,3)*h*Pow(l,2)*m + 96*b*Pow(d,4)*Pow(f,2)*g*h*Pow(l,2)*m - 96*a*Pow(d,3)*Pow(f,3)*g*h*Pow(l,2)*m +
               128*b*Pow(d,4)*e*f*Pow(h,2)*Pow(l,2)*m - 416*b*c*Pow(d,3)*Pow(f,2)*Pow(h,2)*Pow(l,2)*m - 80*a*Pow(d,3)*e*Pow(f,2)*Pow(h,2)*Pow(l,2)*m +
               368*a*c*Pow(d,2)*Pow(f,3)*Pow(h,2)*Pow(l,2)*m - 256*Pow(b,2)*Pow(d,3)*f*Pow(h,3)*Pow(l,2)*m + 1472*a*b*Pow(d,2)*Pow(f,2)*Pow(h,3)*Pow(l,2)*m -
               1152*Pow(a,2)*d*Pow(f,3)*Pow(h,3)*Pow(l,2)*m + 4*Pow(d,5)*Pow(e,3)*k*Pow(l,2)*m + 4*c*Pow(d,4)*Pow(e,2)*f*k*Pow(l,2)*m -
               20*Pow(c,2)*Pow(d,3)*e*Pow(f,2)*k*Pow(l,2)*m + 12*Pow(c,3)*Pow(d,2)*Pow(f,3)*k*Pow(l,2)*m + 208*b*Pow(d,4)*e*f*g*k*Pow(l,2)*m +
               80*b*c*Pow(d,3)*Pow(f,2)*g*k*Pow(l,2)*m - 256*a*Pow(d,3)*e*Pow(f,2)*g*k*Pow(l,2)*m - 32*a*c*Pow(d,2)*Pow(f,3)*g*k*Pow(l,2)*m +
               48*b*Pow(d,4)*Pow(e,2)*h*k*Pow(l,2)*m - 288*b*c*Pow(d,3)*e*f*h*k*Pow(l,2)*m - 256*a*Pow(d,3)*Pow(e,2)*f*h*k*Pow(l,2)*m -
               48*b*Pow(c,2)*Pow(d,2)*Pow(f,2)*h*k*Pow(l,2)*m + 800*a*c*Pow(d,2)*e*Pow(f,2)*h*k*Pow(l,2)*m - 256*a*Pow(c,2)*d*Pow(f,3)*h*k*Pow(l,2)*m -
               1024*Pow(b,2)*Pow(d,3)*f*g*h*k*Pow(l,2)*m + 960*a*b*Pow(d,2)*Pow(f,2)*g*h*k*Pow(l,2)*m - 128*Pow(a,2)*d*Pow(f,3)*g*h*k*Pow(l,2)*m -
               384*Pow(b,2)*Pow(d,3)*e*Pow(h,2)*k*Pow(l,2)*m + 2432*Pow(b,2)*c*Pow(d,2)*f*Pow(h,2)*k*Pow(l,2)*m - 64*a*b*Pow(d,2)*e*f*Pow(h,2)*k*Pow(l,2)*m -
               1856*a*b*c*d*Pow(f,2)*Pow(h,2)*k*Pow(l,2)*m - 704*Pow(a,2)*d*e*Pow(f,2)*Pow(h,2)*k*Pow(l,2)*m + 960*Pow(a,2)*c*Pow(f,3)*Pow(h,2)*k*Pow(l,2)*m +
               512*Pow(b,3)*Pow(d,2)*Pow(h,3)*k*Pow(l,2)*m - 4096*a*Pow(b,2)*d*f*Pow(h,3)*k*Pow(l,2)*m + 1536*Pow(a,2)*b*Pow(f,2)*Pow(h,3)*k*Pow(l,2)*m -
               96*b*c*Pow(d,3)*Pow(e,2)*Pow(k,2)*Pow(l,2)*m - 32*a*Pow(d,3)*Pow(e,3)*Pow(k,2)*Pow(l,2)*m + 48*b*Pow(c,2)*Pow(d,2)*e*f*Pow(k,2)*Pow(l,2)*m +
               160*a*c*Pow(d,2)*Pow(e,2)*f*Pow(k,2)*Pow(l,2)*m + 240*b*Pow(c,3)*d*Pow(f,2)*Pow(k,2)*Pow(l,2)*m -
               176*a*Pow(c,2)*d*e*Pow(f,2)*Pow(k,2)*Pow(l,2)*m - 144*a*Pow(c,3)*Pow(f,3)*Pow(k,2)*Pow(l,2)*m - 256*Pow(b,2)*Pow(d,3)*e*g*Pow(k,2)*Pow(l,2)*m +
               128*Pow(b,2)*c*Pow(d,2)*f*g*Pow(k,2)*Pow(l,2)*m - 640*a*b*Pow(d,2)*e*f*g*Pow(k,2)*Pow(l,2)*m - 1088*a*b*c*d*Pow(f,2)*g*Pow(k,2)*Pow(l,2)*m +
               1024*Pow(a,2)*d*e*Pow(f,2)*g*Pow(k,2)*Pow(l,2)*m + 640*Pow(a,2)*c*Pow(f,3)*g*Pow(k,2)*Pow(l,2)*m +
               896*Pow(b,2)*c*Pow(d,2)*e*h*Pow(k,2)*Pow(l,2)*m - 896*Pow(b,2)*Pow(c,2)*d*f*h*Pow(k,2)*Pow(l,2)*m - 896*a*b*c*d*e*f*h*Pow(k,2)*Pow(l,2)*m +
               896*Pow(a,2)*d*Pow(e,2)*f*h*Pow(k,2)*Pow(l,2)*m + 896*a*b*Pow(c,2)*Pow(f,2)*h*Pow(k,2)*Pow(l,2)*m -
               896*Pow(a,2)*c*e*Pow(f,2)*h*Pow(k,2)*Pow(l,2)*m + 1024*Pow(b,3)*Pow(d,2)*g*h*Pow(k,2)*Pow(l,2)*m + 2560*a*Pow(b,2)*d*f*g*h*Pow(k,2)*Pow(l,2)*m -
               2304*Pow(a,2)*b*Pow(f,2)*g*h*Pow(k,2)*Pow(l,2)*m - 2816*Pow(b,3)*c*d*Pow(h,2)*Pow(k,2)*Pow(l,2)*m +
               512*a*Pow(b,2)*d*e*Pow(h,2)*Pow(k,2)*Pow(l,2)*m + 512*a*Pow(b,2)*c*f*Pow(h,2)*Pow(k,2)*Pow(l,2)*m +
               2304*Pow(a,2)*b*e*f*Pow(h,2)*Pow(k,2)*Pow(l,2)*m + 3072*a*Pow(b,3)*Pow(h,3)*Pow(k,2)*Pow(l,2)*m - 256*Pow(b,2)*Pow(c,2)*d*e*Pow(k,3)*Pow(l,2)*m +
               384*a*b*c*d*Pow(e,2)*Pow(k,3)*Pow(l,2)*m + 64*Pow(a,2)*d*Pow(e,3)*Pow(k,3)*Pow(l,2)*m - 384*Pow(b,2)*Pow(c,3)*f*Pow(k,3)*Pow(l,2)*m +
               832*a*b*Pow(c,2)*e*f*Pow(k,3)*Pow(l,2)*m - 704*Pow(a,2)*c*Pow(e,2)*f*Pow(k,3)*Pow(l,2)*m - 256*Pow(b,3)*c*d*g*Pow(k,3)*Pow(l,2)*m +
               1024*a*Pow(b,2)*d*e*g*Pow(k,3)*Pow(l,2)*m + 1024*a*Pow(b,2)*c*f*g*Pow(k,3)*Pow(l,2)*m - 768*Pow(a,2)*b*e*f*g*Pow(k,3)*Pow(l,2)*m +
               1536*Pow(b,3)*Pow(c,2)*h*Pow(k,3)*Pow(l,2)*m - 1536*a*Pow(b,2)*c*e*h*Pow(k,3)*Pow(l,2)*m - 768*Pow(a,2)*b*Pow(e,2)*h*Pow(k,3)*Pow(l,2)*m -
               3072*a*Pow(b,3)*g*h*Pow(k,3)*Pow(l,2)*m - 32*b*Pow(d,5)*f*h*Pow(l,3)*m + 32*a*Pow(d,4)*Pow(f,2)*h*Pow(l,3)*m - 16*b*Pow(d,5)*e*k*Pow(l,3)*m -
               80*b*c*Pow(d,4)*f*k*Pow(l,3)*m + 32*a*Pow(d,4)*e*f*k*Pow(l,3)*m + 64*a*c*Pow(d,3)*Pow(f,2)*k*Pow(l,3)*m + 64*Pow(b,2)*Pow(d,4)*h*k*Pow(l,3)*m +
               512*a*b*Pow(d,3)*f*h*k*Pow(l,3)*m - 512*Pow(a,2)*Pow(d,2)*Pow(f,2)*h*k*Pow(l,3)*m + 64*Pow(b,2)*c*Pow(d,3)*Pow(k,2)*Pow(l,3)*m +
               128*a*b*Pow(d,3)*e*Pow(k,2)*Pow(l,3)*m + 384*a*b*c*Pow(d,2)*f*Pow(k,2)*Pow(l,3)*m - 256*Pow(a,2)*Pow(d,2)*e*f*Pow(k,2)*Pow(l,3)*m -
               256*Pow(a,2)*c*d*Pow(f,2)*Pow(k,2)*Pow(l,3)*m - 768*a*Pow(b,2)*Pow(d,2)*h*Pow(k,2)*Pow(l,3)*m - 1536*Pow(a,2)*b*d*f*h*Pow(k,2)*Pow(l,3)*m +
               1536*Pow(a,3)*Pow(f,2)*h*Pow(k,2)*Pow(l,3)*m - 256*a*Pow(b,2)*c*d*Pow(k,3)*Pow(l,3)*m - 256*Pow(a,2)*b*d*e*Pow(k,3)*Pow(l,3)*m -
               256*Pow(a,2)*b*c*f*Pow(k,3)*Pow(l,3)*m + 512*Pow(a,3)*e*f*Pow(k,3)*Pow(l,3)*m + 2048*Pow(a,2)*Pow(b,2)*h*Pow(k,3)*Pow(l,3)*m -
               4*Pow(d,4)*Pow(e,2)*Pow(f,2)*Pow(g,2)*Pow(m,2) + 8*c*Pow(d,3)*e*Pow(f,3)*Pow(g,2)*Pow(m,2) - 4*Pow(c,2)*Pow(d,2)*Pow(f,4)*Pow(g,2)*Pow(m,2) -
               16*b*Pow(d,3)*Pow(f,3)*Pow(g,3)*Pow(m,2) + 16*a*Pow(d,2)*Pow(f,4)*Pow(g,3)*Pow(m,2) + 80*Pow(d,4)*Pow(e,3)*f*g*h*Pow(m,2) -
               192*c*Pow(d,3)*Pow(e,2)*Pow(f,2)*g*h*Pow(m,2) + 144*Pow(c,2)*Pow(d,2)*e*Pow(f,3)*g*h*Pow(m,2) - 32*Pow(c,3)*d*Pow(f,4)*g*h*Pow(m,2) +
               368*b*Pow(d,3)*e*Pow(f,2)*Pow(g,2)*h*Pow(m,2) - 80*b*c*Pow(d,2)*Pow(f,3)*Pow(g,2)*h*Pow(m,2) - 416*a*Pow(d,2)*e*Pow(f,3)*Pow(g,2)*h*Pow(m,2) +
               128*a*c*d*Pow(f,4)*Pow(g,2)*h*Pow(m,2) - 16*Pow(d,4)*Pow(e,4)*Pow(h,2)*Pow(m,2) - 32*c*Pow(d,3)*Pow(e,3)*f*Pow(h,2)*Pow(m,2) +
               96*Pow(c,2)*Pow(d,2)*Pow(e,2)*Pow(f,2)*Pow(h,2)*Pow(m,2) - 32*Pow(c,3)*d*e*Pow(f,3)*Pow(h,2)*Pow(m,2) -
               16*Pow(c,4)*Pow(f,4)*Pow(h,2)*Pow(m,2) - 320*b*Pow(d,3)*Pow(e,2)*f*g*Pow(h,2)*Pow(m,2) - 1184*b*c*Pow(d,2)*e*Pow(f,2)*g*Pow(h,2)*Pow(m,2) +
               1088*a*Pow(d,2)*Pow(e,2)*Pow(f,2)*g*Pow(h,2)*Pow(m,2) + 352*b*Pow(c,2)*d*Pow(f,3)*g*Pow(h,2)*Pow(m,2) + 32*a*c*d*e*Pow(f,3)*g*Pow(h,2)*Pow(m,2) +
               32*a*Pow(c,2)*Pow(f,4)*g*Pow(h,2)*Pow(m,2) + 992*Pow(b,2)*Pow(d,2)*Pow(f,2)*Pow(g,2)*Pow(h,2)*Pow(m,2) -
               1216*a*b*d*Pow(f,3)*Pow(g,2)*Pow(h,2)*Pow(m,2) + 128*Pow(a,2)*Pow(f,4)*Pow(g,2)*Pow(h,2)*Pow(m,2) + 64*b*Pow(d,3)*Pow(e,3)*Pow(h,3)*Pow(m,2) +
               576*b*c*Pow(d,2)*Pow(e,2)*f*Pow(h,3)*Pow(m,2) - 256*a*Pow(d,2)*Pow(e,3)*f*Pow(h,3)*Pow(m,2) + 576*b*Pow(c,2)*d*e*Pow(f,2)*Pow(h,3)*Pow(m,2) -
               768*a*c*d*Pow(e,2)*Pow(f,2)*Pow(h,3)*Pow(m,2) + 64*b*Pow(c,3)*Pow(f,3)*Pow(h,3)*Pow(m,2) - 256*a*Pow(c,2)*e*Pow(f,3)*Pow(h,3)*Pow(m,2) -
               320*Pow(b,2)*Pow(d,2)*e*f*g*Pow(h,3)*Pow(m,2) - 1216*Pow(b,2)*c*d*Pow(f,2)*g*Pow(h,3)*Pow(m,2) + 1600*a*b*d*e*Pow(f,2)*g*Pow(h,3)*Pow(m,2) -
               64*a*b*c*Pow(f,3)*g*Pow(h,3)*Pow(m,2) + 640*Pow(a,2)*e*Pow(f,3)*g*Pow(h,3)*Pow(m,2) - 64*Pow(b,2)*Pow(d,2)*Pow(e,2)*Pow(h,4)*Pow(m,2) -
               256*Pow(b,2)*c*d*e*f*Pow(h,4)*Pow(m,2) - 256*a*b*d*Pow(e,2)*f*Pow(h,4)*Pow(m,2) - 64*Pow(b,2)*Pow(c,2)*Pow(f,2)*Pow(h,4)*Pow(m,2) -
               256*a*b*c*e*Pow(f,2)*Pow(h,4)*Pow(m,2) - 64*Pow(a,2)*Pow(e,2)*Pow(f,2)*Pow(h,4)*Pow(m,2) + 512*Pow(b,3)*d*f*g*Pow(h,4)*Pow(m,2) +
               768*a*Pow(b,2)*Pow(f,2)*g*Pow(h,4)*Pow(m,2) + 48*Pow(d,4)*Pow(e,4)*g*k*Pow(m,2) - 184*c*Pow(d,3)*Pow(e,3)*f*g*k*Pow(m,2) +
               272*Pow(c,2)*Pow(d,2)*Pow(e,2)*Pow(f,2)*g*k*Pow(m,2) - 184*Pow(c,3)*d*e*Pow(f,3)*g*k*Pow(m,2) + 48*Pow(c,4)*Pow(f,4)*g*k*Pow(m,2) +
               208*b*Pow(d,3)*Pow(e,2)*f*Pow(g,2)*k*Pow(m,2) - 48*b*c*Pow(d,2)*e*Pow(f,2)*Pow(g,2)*k*Pow(m,2) -
               416*a*Pow(d,2)*Pow(e,2)*Pow(f,2)*Pow(g,2)*k*Pow(m,2) - 16*b*Pow(c,2)*d*Pow(f,3)*Pow(g,2)*k*Pow(m,2) + 416*a*c*d*e*Pow(f,3)*Pow(g,2)*k*Pow(m,2) -
               144*a*Pow(c,2)*Pow(f,4)*Pow(g,2)*k*Pow(m,2) - 480*Pow(b,2)*Pow(d,2)*Pow(f,2)*Pow(g,3)*k*Pow(m,2) + 704*a*b*d*Pow(f,3)*Pow(g,3)*k*Pow(m,2) -
               192*Pow(a,2)*Pow(f,4)*Pow(g,3)*k*Pow(m,2) - 32*c*Pow(d,3)*Pow(e,4)*h*k*Pow(m,2) + 144*Pow(c,2)*Pow(d,2)*Pow(e,3)*f*h*k*Pow(m,2) -
               192*Pow(c,3)*d*Pow(e,2)*Pow(f,2)*h*k*Pow(m,2) + 80*Pow(c,4)*e*Pow(f,3)*h*k*Pow(m,2) - 304*b*Pow(d,3)*Pow(e,3)*g*h*k*Pow(m,2) +
               400*b*c*Pow(d,2)*Pow(e,2)*f*g*h*k*Pow(m,2) + 96*a*Pow(d,2)*Pow(e,3)*f*g*h*k*Pow(m,2) + 400*b*Pow(c,2)*d*e*Pow(f,2)*g*h*k*Pow(m,2) -
               384*a*c*d*Pow(e,2)*Pow(f,2)*g*h*k*Pow(m,2) - 304*b*Pow(c,3)*Pow(f,3)*g*h*k*Pow(m,2) + 96*a*Pow(c,2)*e*Pow(f,3)*g*h*k*Pow(m,2) -
               704*Pow(b,2)*Pow(d,2)*e*f*Pow(g,2)*h*k*Pow(m,2) + 192*Pow(b,2)*c*d*Pow(f,2)*Pow(g,2)*h*k*Pow(m,2) - 64*a*b*d*e*Pow(f,2)*Pow(g,2)*h*k*Pow(m,2) +
               576*a*b*c*Pow(f,3)*Pow(g,2)*h*k*Pow(m,2) - 384*Pow(a,2)*e*Pow(f,3)*Pow(g,2)*h*k*Pow(m,2) + 352*b*c*Pow(d,2)*Pow(e,3)*Pow(h,2)*k*Pow(m,2) +
               32*a*Pow(d,2)*Pow(e,4)*Pow(h,2)*k*Pow(m,2) - 1184*b*Pow(c,2)*d*Pow(e,2)*f*Pow(h,2)*k*Pow(m,2) + 32*a*c*d*Pow(e,3)*f*Pow(h,2)*k*Pow(m,2) -
               320*b*Pow(c,3)*e*Pow(f,2)*Pow(h,2)*k*Pow(m,2) + 1088*a*Pow(c,2)*Pow(e,2)*Pow(f,2)*Pow(h,2)*k*Pow(m,2) +
               352*Pow(b,2)*Pow(d,2)*Pow(e,2)*g*Pow(h,2)*k*Pow(m,2) + 3200*Pow(b,2)*c*d*e*f*g*Pow(h,2)*k*Pow(m,2) - 832*a*b*d*Pow(e,2)*f*g*Pow(h,2)*k*Pow(m,2) +
               352*Pow(b,2)*Pow(c,2)*Pow(f,2)*g*Pow(h,2)*k*Pow(m,2) - 832*a*b*c*e*Pow(f,2)*g*Pow(h,2)*k*Pow(m,2) -
               1664*Pow(a,2)*Pow(e,2)*Pow(f,2)*g*Pow(h,2)*k*Pow(m,2) - 2048*Pow(b,3)*d*f*Pow(g,2)*Pow(h,2)*k*Pow(m,2) +
               512*a*Pow(b,2)*Pow(f,2)*Pow(g,2)*Pow(h,2)*k*Pow(m,2) - 1216*Pow(b,2)*c*d*Pow(e,2)*Pow(h,3)*k*Pow(m,2) - 64*a*b*d*Pow(e,3)*Pow(h,3)*k*Pow(m,2) -
               320*Pow(b,2)*Pow(c,2)*e*f*Pow(h,3)*k*Pow(m,2) + 1600*a*b*c*Pow(e,2)*f*Pow(h,3)*k*Pow(m,2) + 640*Pow(a,2)*Pow(e,3)*f*Pow(h,3)*k*Pow(m,2) +
               1024*Pow(b,3)*d*e*g*Pow(h,3)*k*Pow(m,2) + 1024*Pow(b,3)*c*f*g*Pow(h,3)*k*Pow(m,2) - 4096*a*Pow(b,2)*e*f*g*Pow(h,3)*k*Pow(m,2) +
               512*Pow(b,3)*c*e*Pow(h,4)*k*Pow(m,2) + 768*a*Pow(b,2)*Pow(e,2)*Pow(h,4)*k*Pow(m,2) - 1024*Pow(b,4)*g*Pow(h,4)*k*Pow(m,2) -
               4*Pow(c,2)*Pow(d,2)*Pow(e,4)*Pow(k,2)*Pow(m,2) + 8*Pow(c,3)*d*Pow(e,3)*f*Pow(k,2)*Pow(m,2) - 4*Pow(c,4)*Pow(e,2)*Pow(f,2)*Pow(k,2)*Pow(m,2) -
               16*b*c*Pow(d,2)*Pow(e,3)*g*Pow(k,2)*Pow(m,2) - 144*a*Pow(d,2)*Pow(e,4)*g*Pow(k,2)*Pow(m,2) - 48*b*Pow(c,2)*d*Pow(e,2)*f*g*Pow(k,2)*Pow(m,2) +
               416*a*c*d*Pow(e,3)*f*g*Pow(k,2)*Pow(m,2) + 208*b*Pow(c,3)*e*Pow(f,2)*g*Pow(k,2)*Pow(m,2) - 416*a*Pow(c,2)*Pow(e,2)*Pow(f,2)*g*Pow(k,2)*Pow(m,2) +
               224*Pow(b,2)*Pow(d,2)*Pow(e,2)*Pow(g,2)*Pow(k,2)*Pow(m,2) - 896*Pow(b,2)*c*d*e*f*Pow(g,2)*Pow(k,2)*Pow(m,2) -
               448*a*b*d*Pow(e,2)*f*Pow(g,2)*Pow(k,2)*Pow(m,2) + 224*Pow(b,2)*Pow(c,2)*Pow(f,2)*Pow(g,2)*Pow(k,2)*Pow(m,2) -
               448*a*b*c*e*Pow(f,2)*Pow(g,2)*Pow(k,2)*Pow(m,2) + 1344*Pow(a,2)*Pow(e,2)*Pow(f,2)*Pow(g,2)*Pow(k,2)*Pow(m,2) +
               1536*Pow(b,3)*d*f*Pow(g,3)*Pow(k,2)*Pow(m,2) - 1280*a*Pow(b,2)*Pow(f,2)*Pow(g,3)*Pow(k,2)*Pow(m,2) -
               80*b*Pow(c,2)*d*Pow(e,3)*h*Pow(k,2)*Pow(m,2) + 128*a*c*d*Pow(e,4)*h*Pow(k,2)*Pow(m,2) + 368*b*Pow(c,3)*Pow(e,2)*f*h*Pow(k,2)*Pow(m,2) -
               416*a*Pow(c,2)*Pow(e,3)*f*h*Pow(k,2)*Pow(m,2) + 192*Pow(b,2)*c*d*Pow(e,2)*g*h*Pow(k,2)*Pow(m,2) + 576*a*b*d*Pow(e,3)*g*h*Pow(k,2)*Pow(m,2) -
               704*Pow(b,2)*Pow(c,2)*e*f*g*h*Pow(k,2)*Pow(m,2) - 64*a*b*c*Pow(e,2)*f*g*h*Pow(k,2)*Pow(m,2) - 384*Pow(a,2)*Pow(e,3)*f*g*h*Pow(k,2)*Pow(m,2) -
               1024*Pow(b,3)*d*e*Pow(g,2)*h*Pow(k,2)*Pow(m,2) - 1024*Pow(b,3)*c*f*Pow(g,2)*h*Pow(k,2)*Pow(m,2) +
               4096*a*Pow(b,2)*e*f*Pow(g,2)*h*Pow(k,2)*Pow(m,2) + 992*Pow(b,2)*Pow(c,2)*Pow(e,2)*Pow(h,2)*Pow(k,2)*Pow(m,2) -
               1216*a*b*c*Pow(e,3)*Pow(h,2)*Pow(k,2)*Pow(m,2) + 128*Pow(a,2)*Pow(e,4)*Pow(h,2)*Pow(k,2)*Pow(m,2) -
               2048*Pow(b,3)*c*e*g*Pow(h,2)*Pow(k,2)*Pow(m,2) + 512*a*Pow(b,2)*Pow(e,2)*g*Pow(h,2)*Pow(k,2)*Pow(m,2) +
               2048*Pow(b,4)*Pow(g,2)*Pow(h,2)*Pow(k,2)*Pow(m,2) - 16*b*Pow(c,3)*Pow(e,3)*Pow(k,3)*Pow(m,2) + 16*a*Pow(c,2)*Pow(e,4)*Pow(k,3)*Pow(m,2) -
               480*Pow(b,2)*Pow(c,2)*Pow(e,2)*g*Pow(k,3)*Pow(m,2) + 704*a*b*c*Pow(e,3)*g*Pow(k,3)*Pow(m,2) - 192*Pow(a,2)*Pow(e,4)*g*Pow(k,3)*Pow(m,2) +
               1536*Pow(b,3)*c*e*Pow(g,2)*Pow(k,3)*Pow(m,2) - 1280*a*Pow(b,2)*Pow(e,2)*Pow(g,2)*Pow(k,3)*Pow(m,2) - 1024*Pow(b,4)*Pow(g,3)*Pow(k,3)*Pow(m,2) +
               8*Pow(d,5)*Pow(e,2)*f*g*l*Pow(m,2) - 16*c*Pow(d,4)*e*Pow(f,2)*g*l*Pow(m,2) + 8*Pow(c,2)*Pow(d,3)*Pow(f,3)*g*l*Pow(m,2) +
               48*b*Pow(d,4)*Pow(f,2)*Pow(g,2)*l*Pow(m,2) - 48*a*Pow(d,3)*Pow(f,3)*Pow(g,2)*l*Pow(m,2) + 16*Pow(d,5)*Pow(e,3)*h*l*Pow(m,2) -
               96*c*Pow(d,4)*Pow(e,2)*f*h*l*Pow(m,2) + 144*Pow(c,2)*Pow(d,3)*e*Pow(f,2)*h*l*Pow(m,2) - 64*Pow(c,3)*Pow(d,2)*Pow(f,3)*h*l*Pow(m,2) -
               304*b*Pow(d,4)*e*f*g*h*l*Pow(m,2) - 272*b*c*Pow(d,3)*Pow(f,2)*g*h*l*Pow(m,2) + 400*a*Pow(d,3)*e*Pow(f,2)*g*h*l*Pow(m,2) +
               176*a*c*Pow(d,2)*Pow(f,3)*g*h*l*Pow(m,2) + 32*b*Pow(d,4)*Pow(e,2)*Pow(h,2)*l*Pow(m,2) + 32*b*c*Pow(d,3)*e*f*Pow(h,2)*l*Pow(m,2) +
               352*a*Pow(d,3)*Pow(e,2)*f*Pow(h,2)*l*Pow(m,2) + 1088*b*Pow(c,2)*Pow(d,2)*Pow(f,2)*Pow(h,2)*l*Pow(m,2) -
               1184*a*c*Pow(d,2)*e*Pow(f,2)*Pow(h,2)*l*Pow(m,2) - 320*a*Pow(c,2)*d*Pow(f,3)*Pow(h,2)*l*Pow(m,2) +
               1472*Pow(b,2)*Pow(d,3)*f*g*Pow(h,2)*l*Pow(m,2) - 2752*a*b*Pow(d,2)*Pow(f,2)*g*Pow(h,2)*l*Pow(m,2) +
               1472*Pow(a,2)*d*Pow(f,3)*g*Pow(h,2)*l*Pow(m,2) - 256*Pow(b,2)*Pow(d,3)*e*Pow(h,3)*l*Pow(m,2) - 1664*Pow(b,2)*c*Pow(d,2)*f*Pow(h,3)*l*Pow(m,2) +
               704*a*b*Pow(d,2)*e*f*Pow(h,3)*l*Pow(m,2) - 1088*a*b*c*d*Pow(f,2)*Pow(h,3)*l*Pow(m,2) + 1472*Pow(a,2)*d*e*Pow(f,2)*Pow(h,3)*l*Pow(m,2) +
               192*Pow(a,2)*c*Pow(f,3)*Pow(h,3)*l*Pow(m,2) + 256*Pow(b,3)*Pow(d,2)*Pow(h,4)*l*Pow(m,2) + 1536*a*Pow(b,2)*d*f*Pow(h,4)*l*Pow(m,2) +
               768*Pow(a,2)*b*Pow(f,2)*Pow(h,4)*l*Pow(m,2) - 32*c*Pow(d,4)*Pow(e,3)*k*l*Pow(m,2) + 88*Pow(c,2)*Pow(d,3)*Pow(e,2)*f*k*l*Pow(m,2) -
               80*Pow(c,3)*Pow(d,2)*e*Pow(f,2)*k*l*Pow(m,2) + 24*Pow(c,4)*d*Pow(f,3)*k*l*Pow(m,2) - 232*b*Pow(d,4)*Pow(e,2)*g*k*l*Pow(m,2) -
               64*b*c*Pow(d,3)*e*f*g*k*l*Pow(m,2) + 416*a*Pow(d,3)*Pow(e,2)*f*g*k*l*Pow(m,2) + 8*b*Pow(c,2)*Pow(d,2)*Pow(f,2)*g*k*l*Pow(m,2) -
               208*a*c*Pow(d,2)*e*Pow(f,2)*g*k*l*Pow(m,2) + 80*a*Pow(c,2)*d*Pow(f,3)*g*k*l*Pow(m,2) + 384*Pow(b,2)*Pow(d,3)*f*Pow(g,2)*k*l*Pow(m,2) -
               416*a*b*Pow(d,2)*Pow(f,2)*Pow(g,2)*k*l*Pow(m,2) - 64*Pow(a,2)*d*Pow(f,3)*Pow(g,2)*k*l*Pow(m,2) + 208*b*c*Pow(d,3)*Pow(e,2)*h*k*l*Pow(m,2) +
               32*a*Pow(d,3)*Pow(e,3)*h*k*l*Pow(m,2) + 400*b*Pow(c,2)*Pow(d,2)*e*f*h*k*l*Pow(m,2) - 608*a*c*Pow(d,2)*Pow(e,2)*f*h*k*l*Pow(m,2) -
               800*b*Pow(c,3)*d*Pow(f,2)*h*k*l*Pow(m,2) + 400*a*Pow(c,2)*d*e*Pow(f,2)*h*k*l*Pow(m,2) + 368*a*Pow(c,3)*Pow(f,3)*h*k*l*Pow(m,2) +
               1472*Pow(b,2)*Pow(d,3)*e*g*h*k*l*Pow(m,2) - 1856*Pow(b,2)*c*Pow(d,2)*f*g*h*k*l*Pow(m,2) - 1024*a*b*Pow(d,2)*e*f*g*h*k*l*Pow(m,2) +
               4352*a*b*c*d*Pow(f,2)*g*h*k*l*Pow(m,2) - 64*Pow(a,2)*d*e*Pow(f,2)*g*h*k*l*Pow(m,2) - 2112*Pow(a,2)*c*Pow(f,3)*g*h*k*l*Pow(m,2) -
               640*Pow(b,2)*c*Pow(d,2)*e*Pow(h,2)*k*l*Pow(m,2) - 864*a*b*Pow(d,2)*Pow(e,2)*Pow(h,2)*k*l*Pow(m,2) -
               192*Pow(b,2)*Pow(c,2)*d*f*Pow(h,2)*k*l*Pow(m,2) + 1920*a*b*c*d*e*f*Pow(h,2)*k*l*Pow(m,2) - 640*Pow(a,2)*d*Pow(e,2)*f*Pow(h,2)*k*l*Pow(m,2) +
               32*a*b*Pow(c,2)*Pow(f,2)*Pow(h,2)*k*l*Pow(m,2) - 192*Pow(a,2)*c*e*Pow(f,2)*Pow(h,2)*k*l*Pow(m,2) -
               2944*Pow(b,3)*Pow(d,2)*g*Pow(h,2)*k*l*Pow(m,2) + 3840*a*Pow(b,2)*d*f*g*Pow(h,2)*k*l*Pow(m,2) - 1664*Pow(a,2)*b*Pow(f,2)*g*Pow(h,2)*k*l*Pow(m,2) +
               2816*Pow(b,3)*c*d*Pow(h,3)*k*l*Pow(m,2) + 1280*a*Pow(b,2)*d*e*Pow(h,3)*k*l*Pow(m,2) + 1280*a*Pow(b,2)*c*f*Pow(h,3)*k*l*Pow(m,2) -
               5888*Pow(a,2)*b*e*f*Pow(h,3)*k*l*Pow(m,2) - 4096*a*Pow(b,3)*Pow(h,4)*k*l*Pow(m,2) + 88*b*Pow(c,2)*Pow(d,2)*Pow(e,2)*Pow(k,2)*l*Pow(m,2) +
               96*a*c*Pow(d,2)*Pow(e,3)*Pow(k,2)*l*Pow(m,2) - 256*b*Pow(c,3)*d*e*f*Pow(k,2)*l*Pow(m,2) - 160*a*Pow(c,2)*d*Pow(e,2)*f*Pow(k,2)*l*Pow(m,2) +
               24*b*Pow(c,4)*Pow(f,2)*Pow(k,2)*l*Pow(m,2) + 208*a*Pow(c,3)*e*Pow(f,2)*Pow(k,2)*l*Pow(m,2) - 448*Pow(b,2)*c*Pow(d,2)*e*g*Pow(k,2)*l*Pow(m,2) +
               896*a*b*Pow(d,2)*Pow(e,2)*g*Pow(k,2)*l*Pow(m,2) + 896*Pow(b,2)*Pow(c,2)*d*f*g*Pow(k,2)*l*Pow(m,2) + 1792*a*b*c*d*e*f*g*Pow(k,2)*l*Pow(m,2) -
               1792*Pow(a,2)*d*Pow(e,2)*f*g*Pow(k,2)*l*Pow(m,2) - 896*a*b*Pow(c,2)*Pow(f,2)*g*Pow(k,2)*l*Pow(m,2) -
               448*Pow(a,2)*c*e*Pow(f,2)*g*Pow(k,2)*l*Pow(m,2) - 384*Pow(b,3)*Pow(d,2)*Pow(g,2)*Pow(k,2)*l*Pow(m,2) -
               2304*a*Pow(b,2)*d*f*Pow(g,2)*Pow(k,2)*l*Pow(m,2) + 2432*Pow(a,2)*b*Pow(f,2)*Pow(g,2)*Pow(k,2)*l*Pow(m,2) -
               704*Pow(b,2)*Pow(c,2)*d*e*h*Pow(k,2)*l*Pow(m,2) - 64*a*b*c*d*Pow(e,2)*h*Pow(k,2)*l*Pow(m,2) - 384*Pow(a,2)*d*Pow(e,3)*h*Pow(k,2)*l*Pow(m,2) +
               960*Pow(b,2)*Pow(c,3)*f*h*Pow(k,2)*l*Pow(m,2) - 1856*a*b*Pow(c,2)*e*f*h*Pow(k,2)*l*Pow(m,2) + 2432*Pow(a,2)*c*Pow(e,2)*f*h*Pow(k,2)*l*Pow(m,2) +
               3328*Pow(b,3)*c*d*g*h*Pow(k,2)*l*Pow(m,2) - 4352*a*Pow(b,2)*d*e*g*h*Pow(k,2)*l*Pow(m,2) - 4352*a*Pow(b,2)*c*f*g*h*Pow(k,2)*l*Pow(m,2) +
               2816*Pow(a,2)*b*e*f*g*h*Pow(k,2)*l*Pow(m,2) - 1920*Pow(b,3)*Pow(c,2)*Pow(h,2)*Pow(k,2)*l*Pow(m,2) +
               2816*a*Pow(b,2)*c*e*Pow(h,2)*Pow(k,2)*l*Pow(m,2) + 1408*Pow(a,2)*b*Pow(e,2)*Pow(h,2)*Pow(k,2)*l*Pow(m,2) +
               2048*a*Pow(b,3)*g*Pow(h,2)*Pow(k,2)*l*Pow(m,2) + 576*Pow(b,2)*Pow(c,3)*e*Pow(k,3)*l*Pow(m,2) - 736*a*b*Pow(c,2)*Pow(e,2)*Pow(k,3)*l*Pow(m,2) +
               128*Pow(a,2)*c*Pow(e,3)*Pow(k,3)*l*Pow(m,2) - 1152*Pow(b,3)*Pow(c,2)*g*Pow(k,3)*l*Pow(m,2) + 256*a*Pow(b,2)*c*e*g*Pow(k,3)*l*Pow(m,2) +
               128*Pow(a,2)*b*Pow(e,2)*g*Pow(k,3)*l*Pow(m,2) + 2048*a*Pow(b,3)*Pow(g,2)*Pow(k,3)*l*Pow(m,2) - 4*Pow(d,6)*Pow(e,2)*Pow(l,2)*Pow(m,2) +
               8*c*Pow(d,5)*e*f*Pow(l,2)*Pow(m,2) - 4*Pow(c,2)*Pow(d,4)*Pow(f,2)*Pow(l,2)*Pow(m,2) - 48*b*Pow(d,5)*f*g*Pow(l,2)*Pow(m,2) +
               48*a*Pow(d,4)*Pow(f,2)*g*Pow(l,2)*Pow(m,2) - 64*b*Pow(d,5)*e*h*Pow(l,2)*Pow(m,2) + 352*b*c*Pow(d,4)*f*h*Pow(l,2)*Pow(m,2) +
               16*a*Pow(d,4)*e*f*h*Pow(l,2)*Pow(m,2) - 304*a*c*Pow(d,3)*Pow(f,2)*h*Pow(l,2)*Pow(m,2) + 128*Pow(b,2)*Pow(d,4)*Pow(h,2)*Pow(l,2)*Pow(m,2) -
               1216*a*b*Pow(d,3)*f*Pow(h,2)*Pow(l,2)*Pow(m,2) + 992*Pow(a,2)*Pow(d,2)*Pow(f,2)*Pow(h,2)*Pow(l,2)*Pow(m,2) +
               160*b*c*Pow(d,4)*e*k*Pow(l,2)*Pow(m,2) + 24*a*Pow(d,4)*Pow(e,2)*k*Pow(l,2)*Pow(m,2) - 16*b*Pow(c,2)*Pow(d,3)*f*k*Pow(l,2)*Pow(m,2) -
               256*a*c*Pow(d,3)*e*f*k*Pow(l,2)*Pow(m,2) + 88*a*Pow(c,2)*Pow(d,2)*Pow(f,2)*k*Pow(l,2)*Pow(m,2) + 96*Pow(b,2)*Pow(d,4)*g*k*Pow(l,2)*Pow(m,2) -
               128*a*b*Pow(d,3)*f*g*k*Pow(l,2)*Pow(m,2) + 128*Pow(a,2)*Pow(d,2)*Pow(f,2)*g*k*Pow(l,2)*Pow(m,2) - 832*Pow(b,2)*c*Pow(d,3)*h*k*Pow(l,2)*Pow(m,2) +
               128*a*b*Pow(d,3)*e*h*k*Pow(l,2)*Pow(m,2) - 512*a*b*c*Pow(d,2)*f*h*k*Pow(l,2)*Pow(m,2) + 640*Pow(a,2)*Pow(d,2)*e*f*h*k*Pow(l,2)*Pow(m,2) +
               192*Pow(a,2)*c*d*Pow(f,2)*h*k*Pow(l,2)*Pow(m,2) + 1408*a*Pow(b,2)*Pow(d,2)*Pow(h,2)*k*Pow(l,2)*Pow(m,2) +
               2816*Pow(a,2)*b*d*f*Pow(h,2)*k*Pow(l,2)*Pow(m,2) - 1920*Pow(a,3)*Pow(f,2)*Pow(h,2)*k*Pow(l,2)*Pow(m,2) +
               224*Pow(b,2)*Pow(c,2)*Pow(d,2)*Pow(k,2)*Pow(l,2)*Pow(m,2) - 896*a*b*c*Pow(d,2)*e*Pow(k,2)*Pow(l,2)*Pow(m,2) -
               448*a*b*Pow(c,2)*d*f*Pow(k,2)*Pow(l,2)*Pow(m,2) + 896*Pow(a,2)*c*d*e*f*Pow(k,2)*Pow(l,2)*Pow(m,2) +
               224*Pow(a,2)*Pow(c,2)*Pow(f,2)*Pow(k,2)*Pow(l,2)*Pow(m,2) - 256*a*Pow(b,2)*Pow(d,2)*g*Pow(k,2)*Pow(l,2)*Pow(m,2) +
               1280*Pow(a,2)*b*d*f*g*Pow(k,2)*Pow(l,2)*Pow(m,2) - 1280*Pow(a,3)*Pow(f,2)*g*Pow(k,2)*Pow(l,2)*Pow(m,2) +
               2304*a*Pow(b,2)*c*d*h*Pow(k,2)*Pow(l,2)*Pow(m,2) + 512*Pow(a,2)*b*d*e*h*Pow(k,2)*Pow(l,2)*Pow(m,2) +
               512*Pow(a,2)*b*c*f*h*Pow(k,2)*Pow(l,2)*Pow(m,2) - 2816*Pow(a,3)*e*f*h*Pow(k,2)*Pow(l,2)*Pow(m,2) -
               5632*Pow(a,2)*Pow(b,2)*Pow(h,2)*Pow(k,2)*Pow(l,2)*Pow(m,2) - 384*a*Pow(b,2)*Pow(c,2)*Pow(k,3)*Pow(l,2)*Pow(m,2) +
               1024*Pow(a,2)*b*c*e*Pow(k,3)*Pow(l,2)*Pow(m,2) - 128*Pow(a,3)*Pow(e,2)*Pow(k,3)*Pow(l,2)*Pow(m,2) -
               512*Pow(a,2)*Pow(b,2)*g*Pow(k,3)*Pow(l,2)*Pow(m,2) + 16*b*Pow(d,6)*Pow(l,3)*Pow(m,2) - 16*a*Pow(d,5)*f*Pow(l,3)*Pow(m,2) -
               160*a*b*Pow(d,4)*k*Pow(l,3)*Pow(m,2) + 128*Pow(a,2)*Pow(d,3)*f*k*Pow(l,3)*Pow(m,2) + 512*Pow(a,2)*b*Pow(d,2)*Pow(k,2)*Pow(l,3)*Pow(m,2) -
               256*Pow(a,3)*d*f*Pow(k,2)*Pow(l,3)*Pow(m,2) - 512*Pow(a,3)*b*Pow(k,3)*Pow(l,3)*Pow(m,2) - 32*Pow(d,5)*Pow(e,3)*g*Pow(m,3) +
               80*c*Pow(d,4)*Pow(e,2)*f*g*Pow(m,3) - 64*Pow(c,2)*Pow(d,3)*e*Pow(f,2)*g*Pow(m,3) + 16*Pow(c,3)*Pow(d,2)*Pow(f,3)*g*Pow(m,3) -
               144*b*Pow(d,4)*e*f*Pow(g,2)*Pow(m,3) + 48*b*c*Pow(d,3)*Pow(f,2)*Pow(g,2)*Pow(m,3) + 160*a*Pow(d,3)*e*Pow(f,2)*Pow(g,2)*Pow(m,3) -
               64*a*c*Pow(d,2)*Pow(f,3)*Pow(g,2)*Pow(m,3) + 32*c*Pow(d,4)*Pow(e,3)*h*Pow(m,3) - 32*Pow(c,2)*Pow(d,3)*Pow(e,2)*f*h*Pow(m,3) -
               32*Pow(c,3)*Pow(d,2)*e*Pow(f,2)*h*Pow(m,3) + 32*Pow(c,4)*d*Pow(f,3)*h*Pow(m,3) + 96*b*Pow(d,4)*Pow(e,2)*g*h*Pow(m,3) +
               992*b*c*Pow(d,3)*e*f*g*h*Pow(m,3) - 736*a*Pow(d,3)*Pow(e,2)*f*g*h*Pow(m,3) - 320*b*Pow(c,2)*Pow(d,2)*Pow(f,2)*g*h*Pow(m,3) +
               32*a*c*Pow(d,2)*e*Pow(f,2)*g*h*Pow(m,3) - 64*a*Pow(c,2)*d*Pow(f,3)*g*h*Pow(m,3) - 1152*Pow(b,2)*Pow(d,3)*f*Pow(g,2)*h*Pow(m,3) +
               1472*a*b*Pow(d,2)*Pow(f,2)*Pow(g,2)*h*Pow(m,3) - 256*Pow(a,2)*d*Pow(f,3)*Pow(g,2)*h*Pow(m,3) - 256*b*c*Pow(d,3)*Pow(e,2)*Pow(h,2)*Pow(m,3) +
               64*a*Pow(d,3)*Pow(e,3)*Pow(h,2)*Pow(m,3) - 768*b*Pow(c,2)*Pow(d,2)*e*f*Pow(h,2)*Pow(m,3) + 576*a*c*Pow(d,2)*Pow(e,2)*f*Pow(h,2)*Pow(m,3) -
               256*b*Pow(c,3)*d*Pow(f,2)*Pow(h,2)*Pow(m,3) + 576*a*Pow(c,2)*d*e*Pow(f,2)*Pow(h,2)*Pow(m,3) + 64*a*Pow(c,3)*Pow(f,3)*Pow(h,2)*Pow(m,3) +
               192*Pow(b,2)*Pow(d,3)*e*g*Pow(h,2)*Pow(m,3) + 1472*Pow(b,2)*c*Pow(d,2)*f*g*Pow(h,2)*Pow(m,3) - 1088*a*b*Pow(d,2)*e*f*g*Pow(h,2)*Pow(m,3) +
               704*a*b*c*d*Pow(f,2)*g*Pow(h,2)*Pow(m,3) - 1664*Pow(a,2)*d*e*Pow(f,2)*g*Pow(h,2)*Pow(m,3) - 256*Pow(a,2)*c*Pow(f,3)*g*Pow(h,2)*Pow(m,3) +
               128*Pow(b,2)*c*Pow(d,2)*e*Pow(h,3)*Pow(m,3) + 128*a*b*Pow(d,2)*Pow(e,2)*Pow(h,3)*Pow(m,3) + 128*Pow(b,2)*Pow(c,2)*d*f*Pow(h,3)*Pow(m,3) +
               512*a*b*c*d*e*f*Pow(h,3)*Pow(m,3) + 128*Pow(a,2)*d*Pow(e,2)*f*Pow(h,3)*Pow(m,3) + 128*a*b*Pow(c,2)*Pow(f,2)*Pow(h,3)*Pow(m,3) +
               128*Pow(a,2)*c*e*Pow(f,2)*Pow(h,3)*Pow(m,3) - 256*Pow(b,3)*Pow(d,2)*g*Pow(h,3)*Pow(m,3) - 1536*a*Pow(b,2)*d*f*g*Pow(h,3)*Pow(m,3) -
               768*Pow(a,2)*b*Pow(f,2)*g*Pow(h,3)*Pow(m,3) + 16*Pow(c,2)*Pow(d,3)*Pow(e,3)*k*Pow(m,3) - 64*Pow(c,3)*Pow(d,2)*Pow(e,2)*f*k*Pow(m,3) +
               80*Pow(c,4)*d*e*Pow(f,2)*k*Pow(m,3) - 32*Pow(c,5)*Pow(f,3)*k*Pow(m,3) + 368*b*c*Pow(d,3)*Pow(e,2)*g*k*Pow(m,3) -
               64*a*Pow(d,3)*Pow(e,3)*g*k*Pow(m,3) - 800*b*Pow(c,2)*Pow(d,2)*e*f*g*k*Pow(m,3) + 96*a*c*Pow(d,2)*Pow(e,2)*f*g*k*Pow(m,3) +
               368*b*Pow(c,3)*d*Pow(f,2)*g*k*Pow(m,3) + 96*a*Pow(c,2)*d*e*Pow(f,2)*g*k*Pow(m,3) - 64*a*Pow(c,3)*Pow(f,3)*g*k*Pow(m,3) -
               576*Pow(b,2)*Pow(d,3)*e*Pow(g,2)*k*Pow(m,3) + 960*Pow(b,2)*c*Pow(d,2)*f*Pow(g,2)*k*Pow(m,3) + 1472*a*b*Pow(d,2)*e*f*Pow(g,2)*k*Pow(m,3) -
               2112*a*b*c*d*Pow(f,2)*Pow(g,2)*k*Pow(m,3) - 384*Pow(a,2)*d*e*Pow(f,2)*Pow(g,2)*k*Pow(m,3) + 768*Pow(a,2)*c*Pow(f,3)*Pow(g,2)*k*Pow(m,3) -
               320*b*Pow(c,2)*Pow(d,2)*Pow(e,2)*h*k*Pow(m,3) - 64*a*c*Pow(d,2)*Pow(e,3)*h*k*Pow(m,3) + 992*b*Pow(c,3)*d*e*f*h*k*Pow(m,3) +
               32*a*Pow(c,2)*d*Pow(e,2)*f*h*k*Pow(m,3) + 96*b*Pow(c,4)*Pow(f,2)*h*k*Pow(m,3) - 736*a*Pow(c,3)*e*Pow(f,2)*h*k*Pow(m,3) -
               1024*Pow(b,2)*c*Pow(d,2)*e*g*h*k*Pow(m,3) + 320*a*b*Pow(d,2)*Pow(e,2)*g*h*k*Pow(m,3) - 1024*Pow(b,2)*Pow(c,2)*d*f*g*h*k*Pow(m,3) -
               2304*a*b*c*d*e*f*g*h*k*Pow(m,3) + 1664*Pow(a,2)*d*Pow(e,2)*f*g*h*k*Pow(m,3) + 320*a*b*Pow(c,2)*Pow(f,2)*g*h*k*Pow(m,3) +
               1664*Pow(a,2)*c*e*Pow(f,2)*g*h*k*Pow(m,3) + 2304*Pow(b,3)*Pow(d,2)*Pow(g,2)*h*k*Pow(m,3) - 512*a*Pow(b,2)*d*f*Pow(g,2)*h*k*Pow(m,3) -
               256*Pow(a,2)*b*Pow(f,2)*Pow(g,2)*h*k*Pow(m,3) + 1472*Pow(b,2)*Pow(c,2)*d*e*Pow(h,2)*k*Pow(m,3) + 704*a*b*c*d*Pow(e,2)*Pow(h,2)*k*Pow(m,3) -
               256*Pow(a,2)*d*Pow(e,3)*Pow(h,2)*k*Pow(m,3) + 192*Pow(b,2)*Pow(c,3)*f*Pow(h,2)*k*Pow(m,3) - 1088*a*b*Pow(c,2)*e*f*Pow(h,2)*k*Pow(m,3) -
               1664*Pow(a,2)*c*Pow(e,2)*f*Pow(h,2)*k*Pow(m,3) - 2560*Pow(b,3)*c*d*g*Pow(h,2)*k*Pow(m,3) - 512*a*Pow(b,2)*d*e*g*Pow(h,2)*k*Pow(m,3) -
               512*a*Pow(b,2)*c*f*g*Pow(h,2)*k*Pow(m,3) + 6656*Pow(a,2)*b*e*f*g*Pow(h,2)*k*Pow(m,3) - 256*Pow(b,3)*Pow(c,2)*Pow(h,3)*k*Pow(m,3) -
               1536*a*Pow(b,2)*c*e*Pow(h,3)*k*Pow(m,3) - 768*Pow(a,2)*b*Pow(e,2)*Pow(h,3)*k*Pow(m,3) + 4096*a*Pow(b,3)*g*Pow(h,3)*k*Pow(m,3) +
               48*b*Pow(c,3)*d*Pow(e,2)*Pow(k,2)*Pow(m,3) - 64*a*Pow(c,2)*d*Pow(e,3)*Pow(k,2)*Pow(m,3) - 144*b*Pow(c,4)*e*f*Pow(k,2)*Pow(m,3) +
               160*a*Pow(c,3)*Pow(e,2)*f*Pow(k,2)*Pow(m,3) + 960*Pow(b,2)*Pow(c,2)*d*e*g*Pow(k,2)*Pow(m,3) - 2112*a*b*c*d*Pow(e,2)*g*Pow(k,2)*Pow(m,3) +
               768*Pow(a,2)*d*Pow(e,3)*g*Pow(k,2)*Pow(m,3) - 576*Pow(b,2)*Pow(c,3)*f*g*Pow(k,2)*Pow(m,3) + 1472*a*b*Pow(c,2)*e*f*g*Pow(k,2)*Pow(m,3) -
               384*Pow(a,2)*c*Pow(e,2)*f*g*Pow(k,2)*Pow(m,3) - 1536*Pow(b,3)*c*d*Pow(g,2)*Pow(k,2)*Pow(m,3) + 2560*a*Pow(b,2)*d*e*Pow(g,2)*Pow(k,2)*Pow(m,3) +
               2560*a*Pow(b,2)*c*f*Pow(g,2)*Pow(k,2)*Pow(m,3) - 4608*Pow(a,2)*b*e*f*Pow(g,2)*Pow(k,2)*Pow(m,3) - 1152*Pow(b,2)*Pow(c,3)*e*h*Pow(k,2)*Pow(m,3) +
               1472*a*b*Pow(c,2)*Pow(e,2)*h*Pow(k,2)*Pow(m,3) - 256*Pow(a,2)*c*Pow(e,3)*h*Pow(k,2)*Pow(m,3) + 2304*Pow(b,3)*Pow(c,2)*g*h*Pow(k,2)*Pow(m,3) -
               512*a*Pow(b,2)*c*e*g*h*Pow(k,2)*Pow(m,3) - 256*Pow(a,2)*b*Pow(e,2)*g*h*Pow(k,2)*Pow(m,3) - 4096*a*Pow(b,3)*Pow(g,2)*h*Pow(k,2)*Pow(m,3) +
               16*c*Pow(d,5)*Pow(e,2)*l*Pow(m,3) - 32*Pow(c,2)*Pow(d,4)*e*f*l*Pow(m,3) + 16*Pow(c,3)*Pow(d,3)*Pow(f,2)*l*Pow(m,3) +
               144*b*Pow(d,5)*e*g*l*Pow(m,3) + 48*b*c*Pow(d,4)*f*g*l*Pow(m,3) - 176*a*Pow(d,4)*e*f*g*l*Pow(m,3) - 16*a*c*Pow(d,3)*Pow(f,2)*g*l*Pow(m,3) -
               32*b*c*Pow(d,4)*e*h*l*Pow(m,3) - 128*a*Pow(d,4)*Pow(e,2)*h*l*Pow(m,3) - 736*b*Pow(c,2)*Pow(d,3)*f*h*l*Pow(m,3) + 544*a*c*Pow(d,3)*e*f*h*l*Pow(m,3) +
               352*a*Pow(c,2)*Pow(d,2)*Pow(f,2)*h*l*Pow(m,3) - 576*Pow(b,2)*Pow(d,4)*g*h*l*Pow(m,3) + 1664*a*b*Pow(d,3)*f*g*h*l*Pow(m,3) -
               1216*Pow(a,2)*Pow(d,2)*Pow(f,2)*g*h*l*Pow(m,3) + 640*Pow(b,2)*c*Pow(d,3)*Pow(h,2)*l*Pow(m,3) - 64*a*b*Pow(d,3)*e*Pow(h,2)*l*Pow(m,3) +
               1600*a*b*c*Pow(d,2)*f*Pow(h,2)*l*Pow(m,3) - 1216*Pow(a,2)*Pow(d,2)*e*f*Pow(h,2)*l*Pow(m,3) - 320*Pow(a,2)*c*d*Pow(f,2)*Pow(h,2)*l*Pow(m,3) -
               768*a*Pow(b,2)*Pow(d,2)*Pow(h,3)*l*Pow(m,3) - 1536*Pow(a,2)*b*d*f*Pow(h,3)*l*Pow(m,3) - 256*Pow(a,3)*Pow(f,2)*Pow(h,3)*l*Pow(m,3) -
               304*b*Pow(c,2)*Pow(d,3)*e*k*l*Pow(m,3) + 32*a*c*Pow(d,3)*Pow(e,2)*k*l*Pow(m,3) + 368*b*Pow(c,3)*Pow(d,2)*f*k*l*Pow(m,3) +
               208*a*Pow(c,2)*Pow(d,2)*e*f*k*l*Pow(m,3) - 304*a*Pow(c,3)*d*Pow(f,2)*k*l*Pow(m,3) + 192*Pow(b,2)*c*Pow(d,3)*g*k*l*Pow(m,3) -
               512*a*b*Pow(d,3)*e*g*k*l*Pow(m,3) - 640*a*b*c*Pow(d,2)*f*g*k*l*Pow(m,3) + 128*Pow(a,2)*Pow(d,2)*e*f*g*k*l*Pow(m,3) +
               576*Pow(a,2)*c*d*Pow(f,2)*g*k*l*Pow(m,3) + 832*Pow(b,2)*Pow(c,2)*Pow(d,2)*h*k*l*Pow(m,3) + 640*a*b*c*Pow(d,2)*e*h*k*l*Pow(m,3) +
               384*Pow(a,2)*Pow(d,2)*Pow(e,2)*h*k*l*Pow(m,3) - 256*a*b*Pow(c,2)*d*f*h*k*l*Pow(m,3) - 1152*Pow(a,2)*c*d*e*f*h*k*l*Pow(m,3) -
               64*Pow(a,2)*Pow(c,2)*Pow(f,2)*h*k*l*Pow(m,3) + 1280*a*Pow(b,2)*Pow(d,2)*g*h*k*l*Pow(m,3) - 4608*Pow(a,2)*b*d*f*g*h*k*l*Pow(m,3) +
               2816*Pow(a,3)*Pow(f,2)*g*h*k*l*Pow(m,3) - 5888*a*Pow(b,2)*c*d*Pow(h,2)*k*l*Pow(m,3) + 1280*Pow(a,2)*b*d*e*Pow(h,2)*k*l*Pow(m,3) +
               1280*Pow(a,2)*b*c*f*Pow(h,2)*k*l*Pow(m,3) + 2816*Pow(a,3)*e*f*Pow(h,2)*k*l*Pow(m,3) + 6144*Pow(a,2)*Pow(b,2)*Pow(h,3)*k*l*Pow(m,3) -
               576*Pow(b,2)*Pow(c,3)*d*Pow(k,2)*l*Pow(m,3) + 1472*a*b*Pow(c,2)*d*e*Pow(k,2)*l*Pow(m,3) - 384*Pow(a,2)*c*d*Pow(e,2)*Pow(k,2)*l*Pow(m,3) +
               192*a*b*Pow(c,3)*f*Pow(k,2)*l*Pow(m,3) - 832*Pow(a,2)*Pow(c,2)*e*f*Pow(k,2)*l*Pow(m,3) - 256*a*Pow(b,2)*c*d*g*Pow(k,2)*l*Pow(m,3) -
               256*Pow(a,2)*b*d*e*g*Pow(k,2)*l*Pow(m,3) - 256*Pow(a,2)*b*c*f*g*Pow(k,2)*l*Pow(m,3) + 2304*Pow(a,3)*e*f*g*Pow(k,2)*l*Pow(m,3) +
               1536*a*Pow(b,2)*Pow(c,2)*h*Pow(k,2)*l*Pow(m,3) - 4096*Pow(a,2)*b*c*e*h*Pow(k,2)*l*Pow(m,3) + 512*Pow(a,3)*Pow(e,2)*h*Pow(k,2)*l*Pow(m,3) +
               2048*Pow(a,2)*Pow(b,2)*g*h*Pow(k,2)*l*Pow(m,3) - 96*b*c*Pow(d,5)*Pow(l,2)*Pow(m,3) + 16*a*Pow(d,5)*e*Pow(l,2)*Pow(m,3) +
               80*a*c*Pow(d,4)*f*Pow(l,2)*Pow(m,3) + 320*a*b*Pow(d,4)*h*Pow(l,2)*Pow(m,3) - 256*Pow(a,2)*Pow(d,3)*f*h*Pow(l,2)*Pow(m,3) +
               640*a*b*c*Pow(d,3)*k*Pow(l,2)*Pow(m,3) - 128*Pow(a,2)*Pow(d,3)*e*k*Pow(l,2)*Pow(m,3) - 384*Pow(a,2)*c*Pow(d,2)*f*k*Pow(l,2)*Pow(m,3) -
               2048*Pow(a,2)*b*Pow(d,2)*h*k*Pow(l,2)*Pow(m,3) + 1024*Pow(a,3)*d*f*h*k*Pow(l,2)*Pow(m,3) - 1024*Pow(a,2)*b*c*d*Pow(k,2)*Pow(l,2)*Pow(m,3) +
               256*Pow(a,3)*d*e*Pow(k,2)*Pow(l,2)*Pow(m,3) + 256*Pow(a,3)*c*f*Pow(k,2)*Pow(l,2)*Pow(m,3) + 3072*Pow(a,3)*b*h*Pow(k,2)*Pow(l,2)*Pow(m,3) -
               16*Pow(c,2)*Pow(d,4)*Pow(e,2)*Pow(m,4) + 32*Pow(c,3)*Pow(d,3)*e*f*Pow(m,4) - 16*Pow(c,4)*Pow(d,2)*Pow(f,2)*Pow(m,4) -
               288*b*c*Pow(d,4)*e*g*Pow(m,4) + 192*a*Pow(d,4)*Pow(e,2)*g*Pow(m,4) + 96*b*Pow(c,2)*Pow(d,3)*f*g*Pow(m,4) - 32*a*c*Pow(d,3)*e*f*g*Pow(m,4) +
               32*a*Pow(c,2)*Pow(d,2)*Pow(f,2)*g*Pow(m,4) + 432*Pow(b,2)*Pow(d,4)*Pow(g,2)*Pow(m,4) - 576*a*b*Pow(d,3)*f*Pow(g,2)*Pow(m,4) +
               128*Pow(a,2)*Pow(d,2)*Pow(f,2)*Pow(g,2)*Pow(m,4) + 320*b*Pow(c,2)*Pow(d,3)*e*h*Pow(m,4) - 128*a*c*Pow(d,3)*Pow(e,2)*h*Pow(m,4) +
               320*b*Pow(c,3)*Pow(d,2)*f*h*Pow(m,4) - 384*a*Pow(c,2)*Pow(d,2)*e*f*h*Pow(m,4) - 128*a*Pow(c,3)*d*Pow(f,2)*h*Pow(m,4) -
               576*Pow(b,2)*c*Pow(d,3)*g*h*Pow(m,4) + 192*a*b*Pow(d,3)*e*g*h*Pow(m,4) - 1216*a*b*c*Pow(d,2)*f*g*h*Pow(m,4) +
               1408*Pow(a,2)*Pow(d,2)*e*f*g*h*Pow(m,4) + 512*Pow(a,2)*c*d*Pow(f,2)*g*h*Pow(m,4) - 64*Pow(b,2)*Pow(c,2)*Pow(d,2)*Pow(h,2)*Pow(m,4) -
               256*a*b*c*Pow(d,2)*e*Pow(h,2)*Pow(m,4) - 64*Pow(a,2)*Pow(d,2)*Pow(e,2)*Pow(h,2)*Pow(m,4) - 256*a*b*Pow(c,2)*d*f*Pow(h,2)*Pow(m,4) -
               256*Pow(a,2)*c*d*e*f*Pow(h,2)*Pow(m,4) - 64*Pow(a,2)*Pow(c,2)*Pow(f,2)*Pow(h,2)*Pow(m,4) + 768*a*Pow(b,2)*Pow(d,2)*g*Pow(h,2)*Pow(m,4) +
               1536*Pow(a,2)*b*d*f*g*Pow(h,2)*Pow(m,4) + 256*Pow(a,3)*Pow(f,2)*g*Pow(h,2)*Pow(m,4) + 96*b*Pow(c,3)*Pow(d,2)*e*k*Pow(m,4) +
               32*a*Pow(c,2)*Pow(d,2)*Pow(e,2)*k*Pow(m,4) - 288*b*Pow(c,4)*d*f*k*Pow(m,4) - 32*a*Pow(c,3)*d*e*f*k*Pow(m,4) + 192*a*Pow(c,4)*Pow(f,2)*k*Pow(m,4) +
               96*Pow(b,2)*Pow(c,2)*Pow(d,2)*g*k*Pow(m,4) + 832*a*b*c*Pow(d,2)*e*g*k*Pow(m,4) - 576*Pow(a,2)*Pow(d,2)*Pow(e,2)*g*k*Pow(m,4) +
               832*a*b*Pow(c,2)*d*f*g*k*Pow(m,4) - 512*Pow(a,2)*c*d*e*f*g*k*Pow(m,4) - 576*Pow(a,2)*Pow(c,2)*Pow(f,2)*g*k*Pow(m,4) -
               2304*a*Pow(b,2)*Pow(d,2)*Pow(g,2)*k*Pow(m,4) + 2560*Pow(a,2)*b*d*f*Pow(g,2)*k*Pow(m,4) - 768*Pow(a,3)*Pow(f,2)*Pow(g,2)*k*Pow(m,4) -
               576*Pow(b,2)*Pow(c,3)*d*h*k*Pow(m,4) - 1216*a*b*Pow(c,2)*d*e*h*k*Pow(m,4) + 512*Pow(a,2)*c*d*Pow(e,2)*h*k*Pow(m,4) + 192*a*b*Pow(c,3)*f*h*k*Pow(m,4) +
               1408*Pow(a,2)*Pow(c,2)*e*f*h*k*Pow(m,4) + 5120*a*Pow(b,2)*c*d*g*h*k*Pow(m,4) - 2048*Pow(a,2)*b*d*e*g*h*k*Pow(m,4) - 2048*Pow(a,2)*b*c*f*g*h*k*Pow(m,4) -
               3072*Pow(a,3)*e*f*g*h*k*Pow(m,4) + 768*a*Pow(b,2)*Pow(c,2)*Pow(h,2)*k*Pow(m,4) + 1536*Pow(a,2)*b*c*e*Pow(h,2)*k*Pow(m,4) +
               256*Pow(a,3)*Pow(e,2)*Pow(h,2)*k*Pow(m,4) - 6144*Pow(a,2)*Pow(b,2)*g*Pow(h,2)*k*Pow(m,4) + 432*Pow(b,2)*Pow(c,4)*Pow(k,2)*Pow(m,4) -
               576*a*b*Pow(c,3)*e*Pow(k,2)*Pow(m,4) + 128*Pow(a,2)*Pow(c,2)*Pow(e,2)*Pow(k,2)*Pow(m,4) - 2304*a*Pow(b,2)*Pow(c,2)*g*Pow(k,2)*Pow(m,4) +
               2560*Pow(a,2)*b*c*e*g*Pow(k,2)*Pow(m,4) - 768*Pow(a,3)*Pow(e,2)*g*Pow(k,2)*Pow(m,4) + 2048*Pow(a,2)*Pow(b,2)*Pow(g,2)*Pow(k,2)*Pow(m,4) +
               192*b*Pow(c,2)*Pow(d,4)*l*Pow(m,4) - 64*a*c*Pow(d,4)*e*l*Pow(m,4) - 128*a*Pow(c,2)*Pow(d,3)*f*l*Pow(m,4) - 288*a*b*Pow(d,4)*g*l*Pow(m,4) +
               320*Pow(a,2)*Pow(d,3)*f*g*l*Pow(m,4) - 704*a*b*c*Pow(d,3)*h*l*Pow(m,4) + 320*Pow(a,2)*Pow(d,3)*e*h*l*Pow(m,4) +
               64*Pow(a,2)*c*Pow(d,2)*f*h*l*Pow(m,4) + 768*Pow(a,2)*b*Pow(d,2)*Pow(h,2)*l*Pow(m,4) + 512*Pow(a,3)*d*f*Pow(h,2)*l*Pow(m,4) -
               928*a*b*Pow(c,2)*Pow(d,2)*k*l*Pow(m,4) + 192*Pow(a,2)*c*Pow(d,2)*e*k*l*Pow(m,4) + 640*Pow(a,2)*Pow(c,2)*d*f*k*l*Pow(m,4) +
               1664*Pow(a,2)*b*Pow(d,2)*g*k*l*Pow(m,4) - 1280*Pow(a,3)*d*f*g*k*l*Pow(m,4) + 3328*Pow(a,2)*b*c*d*h*k*l*Pow(m,4) - 1280*Pow(a,3)*d*e*h*k*l*Pow(m,4) -
               1280*Pow(a,3)*c*f*h*k*l*Pow(m,4) - 4096*Pow(a,3)*b*Pow(h,2)*k*l*Pow(m,4) + 384*Pow(a,2)*b*Pow(c,2)*Pow(k,2)*l*Pow(m,4) +
               256*Pow(a,3)*c*e*Pow(k,2)*l*Pow(m,4) - 2048*Pow(a,3)*b*g*Pow(k,2)*l*Pow(m,4) - 16*Pow(a,2)*Pow(d,4)*Pow(l,2)*Pow(m,4) +
               128*Pow(a,3)*Pow(d,2)*k*Pow(l,2)*Pow(m,4) - 256*Pow(a,4)*Pow(k,2)*Pow(l,2)*Pow(m,4) - 128*b*Pow(c,3)*Pow(d,3)*Pow(m,5) +
               64*a*Pow(c,2)*Pow(d,3)*e*Pow(m,5) + 64*a*Pow(c,3)*Pow(d,2)*f*Pow(m,5) + 576*a*b*c*Pow(d,3)*g*Pow(m,5) - 384*Pow(a,2)*Pow(d,3)*e*g*Pow(m,5) -
               256*Pow(a,2)*c*Pow(d,2)*f*g*Pow(m,5) + 128*a*b*Pow(c,2)*Pow(d,2)*h*Pow(m,5) + 128*Pow(a,2)*c*Pow(d,2)*e*h*Pow(m,5) +
               128*Pow(a,2)*Pow(c,2)*d*f*h*Pow(m,5) - 768*Pow(a,2)*b*Pow(d,2)*g*h*Pow(m,5) - 512*Pow(a,3)*d*f*g*h*Pow(m,5) + 576*a*b*Pow(c,3)*d*k*Pow(m,5) -
               256*Pow(a,2)*Pow(c,2)*d*e*k*Pow(m,5) - 384*Pow(a,2)*Pow(c,3)*f*k*Pow(m,5) - 2560*Pow(a,2)*b*c*d*g*k*Pow(m,5) + 1536*Pow(a,3)*d*e*g*k*Pow(m,5) +
               1536*Pow(a,3)*c*f*g*k*Pow(m,5) - 768*Pow(a,2)*b*Pow(c,2)*h*k*Pow(m,5) - 512*Pow(a,3)*c*e*h*k*Pow(m,5) + 4096*Pow(a,3)*b*g*h*k*Pow(m,5) +
               64*Pow(a,2)*c*Pow(d,3)*l*Pow(m,5) - 256*Pow(a,3)*Pow(d,2)*h*l*Pow(m,5) - 256*Pow(a,3)*c*d*k*l*Pow(m,5) + 1024*Pow(a,4)*h*k*l*Pow(m,5) -
               64*Pow(a,2)*Pow(c,2)*Pow(d,2)*Pow(m,6) + 256*Pow(a,3)*Pow(d,2)*g*Pow(m,6) + 256*Pow(a,3)*Pow(c,2)*k*Pow(m,6) - 1024*Pow(a,4)*g*k*Pow(m,6) -
               4*Pow(d,3)*Pow(e,2)*Pow(f,3)*Pow(g,3)*n + 8*c*Pow(d,2)*e*Pow(f,4)*Pow(g,3)*n - 4*Pow(c,2)*d*Pow(f,5)*Pow(g,3)*n -
               16*b*Pow(d,2)*Pow(f,4)*Pow(g,4)*n + 16*a*d*Pow(f,5)*Pow(g,4)*n + 12*Pow(d,3)*Pow(e,3)*Pow(f,2)*Pow(g,2)*h*n -
               20*c*Pow(d,2)*Pow(e,2)*Pow(f,3)*Pow(g,2)*h*n + 4*Pow(c,2)*d*e*Pow(f,4)*Pow(g,2)*h*n + 4*Pow(c,3)*Pow(f,5)*Pow(g,2)*h*n +
               64*b*Pow(d,2)*e*Pow(f,3)*Pow(g,3)*h*n + 32*b*c*d*Pow(f,4)*Pow(g,3)*h*n - 80*a*d*e*Pow(f,4)*Pow(g,3)*h*n - 16*a*c*Pow(f,5)*Pow(g,3)*h*n +
               24*Pow(d,3)*Pow(e,4)*f*g*Pow(h,2)*n - 80*c*Pow(d,2)*Pow(e,3)*Pow(f,2)*g*Pow(h,2)*n + 88*Pow(c,2)*d*Pow(e,2)*Pow(f,3)*g*Pow(h,2)*n -
               32*Pow(c,3)*e*Pow(f,4)*g*Pow(h,2)*n + 88*b*Pow(d,2)*Pow(e,2)*Pow(f,2)*Pow(g,2)*Pow(h,2)*n - 256*b*c*d*e*Pow(f,3)*Pow(g,2)*Pow(h,2)*n -
               16*a*d*Pow(e,2)*Pow(f,3)*Pow(g,2)*Pow(h,2)*n + 24*b*Pow(c,2)*Pow(f,4)*Pow(g,2)*Pow(h,2)*n + 160*a*c*e*Pow(f,4)*Pow(g,2)*Pow(h,2)*n +
               128*Pow(b,2)*d*Pow(f,3)*Pow(g,3)*Pow(h,2)*n - 160*a*b*Pow(f,4)*Pow(g,3)*Pow(h,2)*n - 32*Pow(d,3)*Pow(e,5)*Pow(h,3)*n +
               80*c*Pow(d,2)*Pow(e,4)*f*Pow(h,3)*n - 64*Pow(c,2)*d*Pow(e,3)*Pow(f,2)*Pow(h,3)*n + 16*Pow(c,3)*Pow(e,2)*Pow(f,3)*Pow(h,3)*n -
               304*b*Pow(d,2)*Pow(e,3)*f*g*Pow(h,3)*n + 208*b*c*d*Pow(e,2)*Pow(f,2)*g*Pow(h,3)*n + 368*a*d*Pow(e,3)*Pow(f,2)*g*Pow(h,3)*n +
               32*b*Pow(c,2)*e*Pow(f,3)*g*Pow(h,3)*n - 304*a*c*Pow(e,2)*Pow(f,3)*g*Pow(h,3)*n - 384*Pow(b,2)*d*e*Pow(f,2)*Pow(g,2)*Pow(h,3)*n -
               128*Pow(b,2)*c*Pow(f,3)*Pow(g,2)*Pow(h,3)*n + 640*a*b*e*Pow(f,3)*Pow(g,2)*Pow(h,3)*n + 192*b*Pow(d,2)*Pow(e,4)*Pow(h,4)*n -
               32*b*c*d*Pow(e,3)*f*Pow(h,4)*n - 288*a*d*Pow(e,4)*f*Pow(h,4)*n + 32*b*Pow(c,2)*Pow(e,2)*Pow(f,2)*Pow(h,4)*n + 96*a*c*Pow(e,3)*Pow(f,2)*Pow(h,4)*n +
               640*Pow(b,2)*d*Pow(e,2)*f*g*Pow(h,4)*n + 192*Pow(b,2)*c*e*Pow(f,2)*g*Pow(h,4)*n - 928*a*b*Pow(e,2)*Pow(f,2)*g*Pow(h,4)*n +
               128*Pow(b,3)*Pow(f,2)*Pow(g,2)*Pow(h,4)*n - 384*Pow(b,2)*d*Pow(e,3)*Pow(h,5)*n - 256*Pow(b,2)*c*Pow(e,2)*f*Pow(h,5)*n +
               576*a*b*Pow(e,3)*f*Pow(h,5)*n - 256*Pow(b,3)*e*f*g*Pow(h,5)*n + 256*Pow(b,3)*Pow(e,2)*Pow(h,6)*n - 36*Pow(d,3)*Pow(e,4)*f*Pow(g,2)*k*n +
               92*c*Pow(d,2)*Pow(e,3)*Pow(f,2)*Pow(g,2)*k*n - 76*Pow(c,2)*d*Pow(e,2)*Pow(f,3)*Pow(g,2)*k*n + 20*Pow(c,3)*e*Pow(f,4)*Pow(g,2)*k*n -
               184*b*Pow(d,2)*Pow(e,2)*Pow(f,2)*Pow(g,3)*k*n + 128*b*c*d*e*Pow(f,3)*Pow(g,3)*k*n + 176*a*d*Pow(e,2)*Pow(f,3)*Pow(g,3)*k*n -
               40*b*Pow(c,2)*Pow(f,4)*Pow(g,3)*k*n - 80*a*c*e*Pow(f,4)*Pow(g,3)*k*n - 128*Pow(b,2)*d*Pow(f,3)*Pow(g,4)*k*n + 160*a*b*Pow(f,4)*Pow(g,4)*k*n +
               36*Pow(d,3)*Pow(e,5)*g*h*k*n - 76*c*Pow(d,2)*Pow(e,4)*f*g*h*k*n + 44*Pow(c,2)*d*Pow(e,3)*Pow(f,2)*g*h*k*n - 4*Pow(c,3)*Pow(e,2)*Pow(f,3)*g*h*k*n +
               368*b*Pow(d,2)*Pow(e,3)*f*Pow(g,2)*h*k*n - 16*b*c*d*Pow(e,2)*Pow(f,2)*Pow(g,2)*h*k*n - 528*a*d*Pow(e,3)*Pow(f,2)*Pow(g,2)*h*k*n +
               32*b*Pow(c,2)*e*Pow(f,3)*Pow(g,2)*h*k*n + 144*a*c*Pow(e,2)*Pow(f,3)*Pow(g,2)*h*k*n + 384*Pow(b,2)*d*e*Pow(f,2)*Pow(g,3)*h*k*n +
               128*Pow(b,2)*c*Pow(f,3)*Pow(g,3)*h*k*n - 640*a*b*e*Pow(f,3)*Pow(g,3)*h*k*n + 8*c*Pow(d,2)*Pow(e,5)*Pow(h,2)*k*n -
               16*Pow(c,2)*d*Pow(e,4)*f*Pow(h,2)*k*n + 8*Pow(c,3)*Pow(e,3)*Pow(f,2)*Pow(h,2)*k*n - 232*b*Pow(d,2)*Pow(e,4)*g*Pow(h,2)*k*n +
               160*b*c*d*Pow(e,3)*f*g*Pow(h,2)*k*n + 208*a*d*Pow(e,4)*f*g*Pow(h,2)*k*n - 216*b*Pow(c,2)*Pow(e,2)*Pow(f,2)*g*Pow(h,2)*k*n +
               80*a*c*Pow(e,3)*Pow(f,2)*g*Pow(h,2)*k*n - 896*Pow(b,2)*d*Pow(e,2)*f*Pow(g,2)*Pow(h,2)*k*n + 896*a*b*Pow(e,2)*Pow(f,2)*Pow(g,2)*Pow(h,2)*k*n -
               256*Pow(b,3)*Pow(f,2)*Pow(g,3)*Pow(h,2)*k*n - 176*b*c*d*Pow(e,4)*Pow(h,3)*k*n + 144*a*d*Pow(e,5)*Pow(h,3)*k*n -
               16*b*Pow(c,2)*Pow(e,3)*f*Pow(h,3)*k*n + 48*a*c*Pow(e,4)*f*Pow(h,3)*k*n + 640*Pow(b,2)*d*Pow(e,3)*g*Pow(h,3)*k*n +
               128*Pow(b,2)*c*Pow(e,2)*f*g*Pow(h,3)*k*n - 512*a*b*Pow(e,3)*f*g*Pow(h,3)*k*n + 512*Pow(b,3)*e*f*Pow(g,2)*Pow(h,3)*k*n +
               320*Pow(b,2)*c*Pow(e,3)*Pow(h,4)*k*n - 288*a*b*Pow(e,4)*Pow(h,4)*k*n - 640*Pow(b,3)*Pow(e,2)*g*Pow(h,4)*k*n - 12*c*Pow(d,2)*Pow(e,5)*g*Pow(k,2)*n +
               24*Pow(c,2)*d*Pow(e,4)*f*g*Pow(k,2)*n - 12*Pow(c,3)*Pow(e,3)*Pow(f,2)*g*Pow(k,2)*n + 24*b*Pow(d,2)*Pow(e,4)*Pow(g,2)*Pow(k,2)*n -
               256*b*c*d*Pow(e,3)*f*Pow(g,2)*Pow(k,2)*n + 160*a*d*Pow(e,4)*f*Pow(g,2)*Pow(k,2)*n + 88*b*Pow(c,2)*Pow(e,2)*Pow(f,2)*Pow(g,2)*Pow(k,2)*n -
               16*a*c*Pow(e,3)*Pow(f,2)*Pow(g,2)*Pow(k,2)*n + 256*Pow(b,2)*d*Pow(e,2)*f*Pow(g,3)*Pow(k,2)*n - 192*Pow(b,2)*c*e*Pow(f,2)*Pow(g,3)*Pow(k,2)*n +
               32*a*b*Pow(e,2)*Pow(f,2)*Pow(g,3)*Pow(k,2)*n + 128*Pow(b,3)*Pow(f,2)*Pow(g,4)*Pow(k,2)*n + 208*b*c*d*Pow(e,4)*g*h*Pow(k,2)*n -
               160*a*d*Pow(e,5)*g*h*Pow(k,2)*n + 80*b*Pow(c,2)*Pow(e,3)*f*g*h*Pow(k,2)*n - 128*a*c*Pow(e,4)*f*g*h*Pow(k,2)*n -
               256*Pow(b,2)*d*Pow(e,3)*Pow(g,2)*h*Pow(k,2)*n + 128*Pow(b,2)*c*Pow(e,2)*f*Pow(g,2)*h*Pow(k,2)*n - 64*a*b*Pow(e,3)*f*Pow(g,2)*h*Pow(k,2)*n -
               256*Pow(b,3)*e*f*Pow(g,3)*h*Pow(k,2)*n + 48*b*Pow(c,2)*Pow(e,4)*Pow(h,2)*Pow(k,2)*n - 48*a*c*Pow(e,5)*Pow(h,2)*Pow(k,2)*n -
               512*Pow(b,2)*c*Pow(e,3)*g*Pow(h,2)*Pow(k,2)*n + 416*a*b*Pow(e,4)*g*Pow(h,2)*Pow(k,2)*n + 512*Pow(b,3)*Pow(e,2)*Pow(g,2)*Pow(h,2)*Pow(k,2)*n -
               64*b*Pow(c,2)*Pow(e,4)*g*Pow(k,3)*n + 64*a*c*Pow(e,5)*g*Pow(k,3)*n + 192*Pow(b,2)*c*Pow(e,3)*Pow(g,2)*Pow(k,3)*n -
               128*a*b*Pow(e,4)*Pow(g,2)*Pow(k,3)*n - 128*Pow(b,3)*Pow(e,2)*Pow(g,3)*Pow(k,3)*n + 12*Pow(d,4)*Pow(e,2)*Pow(f,2)*Pow(g,2)*l*n -
               24*c*Pow(d,3)*e*Pow(f,3)*Pow(g,2)*l*n + 12*Pow(c,2)*Pow(d,2)*Pow(f,4)*Pow(g,2)*l*n + 64*b*Pow(d,3)*Pow(f,3)*Pow(g,3)*l*n -
               64*a*Pow(d,2)*Pow(f,4)*Pow(g,3)*l*n + 12*Pow(d,4)*Pow(e,3)*f*g*h*l*n - 68*c*Pow(d,3)*Pow(e,2)*Pow(f,2)*g*h*l*n +
               100*Pow(c,2)*Pow(d,2)*e*Pow(f,3)*g*h*l*n - 44*Pow(c,3)*d*Pow(f,4)*g*h*l*n - 32*b*Pow(d,3)*e*Pow(f,2)*Pow(g,2)*h*l*n -
               256*b*c*Pow(d,2)*Pow(f,3)*Pow(g,2)*h*l*n + 80*a*Pow(d,2)*e*Pow(f,3)*Pow(g,2)*h*l*n + 208*a*c*d*Pow(f,4)*Pow(g,2)*h*l*n +
               48*Pow(d,4)*Pow(e,4)*Pow(h,2)*l*n - 184*c*Pow(d,3)*Pow(e,3)*f*Pow(h,2)*l*n + 272*Pow(c,2)*Pow(d,2)*Pow(e,2)*Pow(f,2)*Pow(h,2)*l*n -
               184*Pow(c,3)*d*e*Pow(f,3)*Pow(h,2)*l*n + 48*Pow(c,4)*Pow(f,4)*Pow(h,2)*l*n + 80*b*Pow(d,3)*Pow(e,2)*f*g*Pow(h,2)*l*n -
               208*b*c*Pow(d,2)*e*Pow(f,2)*g*Pow(h,2)*l*n + 8*a*Pow(d,2)*Pow(e,2)*Pow(f,2)*g*Pow(h,2)*l*n + 416*b*Pow(c,2)*d*Pow(f,3)*g*Pow(h,2)*l*n -
               64*a*c*d*e*Pow(f,3)*g*Pow(h,2)*l*n - 232*a*Pow(c,2)*Pow(f,4)*g*Pow(h,2)*l*n + 128*Pow(b,2)*Pow(d,2)*Pow(f,2)*Pow(g,2)*Pow(h,2)*l*n -
               128*a*b*d*Pow(f,3)*Pow(g,2)*Pow(h,2)*l*n + 96*Pow(a,2)*Pow(f,4)*Pow(g,2)*Pow(h,2)*l*n - 64*b*Pow(d,3)*Pow(e,3)*Pow(h,3)*l*n +
               96*b*c*Pow(d,2)*Pow(e,2)*f*Pow(h,3)*l*n + 368*a*Pow(d,2)*Pow(e,3)*f*Pow(h,3)*l*n + 96*b*Pow(c,2)*d*e*Pow(f,2)*Pow(h,3)*l*n -
               800*a*c*d*Pow(e,2)*Pow(f,2)*Pow(h,3)*l*n - 64*b*Pow(c,3)*Pow(f,3)*Pow(h,3)*l*n + 368*a*Pow(c,2)*e*Pow(f,3)*Pow(h,3)*l*n +
               576*Pow(b,2)*Pow(d,2)*e*f*g*Pow(h,3)*l*n + 128*Pow(b,2)*c*d*Pow(f,2)*g*Pow(h,3)*l*n - 640*a*b*d*e*Pow(f,2)*g*Pow(h,3)*l*n -
               512*a*b*c*Pow(f,3)*g*Pow(h,3)*l*n + 192*Pow(a,2)*e*Pow(f,3)*g*Pow(h,3)*l*n - 576*Pow(b,2)*Pow(d,2)*Pow(e,2)*Pow(h,4)*l*n -
               512*Pow(b,2)*c*d*e*f*Pow(h,4)*l*n + 832*a*b*d*Pow(e,2)*f*Pow(h,4)*l*n - 576*Pow(b,2)*Pow(c,2)*Pow(f,2)*Pow(h,4)*l*n +
               832*a*b*c*e*Pow(f,2)*Pow(h,4)*l*n + 96*Pow(a,2)*Pow(e,2)*Pow(f,2)*Pow(h,4)*l*n - 1280*Pow(b,3)*d*f*g*Pow(h,4)*l*n +
               1664*a*Pow(b,2)*Pow(f,2)*g*Pow(h,4)*l*n + 1536*Pow(b,3)*d*e*Pow(h,5)*l*n + 1536*Pow(b,3)*c*f*Pow(h,5)*l*n - 2560*a*Pow(b,2)*e*f*Pow(h,5)*l*n -
               1024*Pow(b,4)*Pow(h,6)*l*n - 18*Pow(d,4)*Pow(e,4)*g*k*l*n + 104*c*Pow(d,3)*Pow(e,3)*f*g*k*l*n - 172*Pow(c,2)*Pow(d,2)*Pow(e,2)*Pow(f,2)*g*k*l*n +
               104*Pow(c,3)*d*e*Pow(f,3)*g*k*l*n - 18*Pow(c,4)*Pow(f,4)*g*k*l*n + 80*b*Pow(d,3)*Pow(e,2)*f*Pow(g,2)*k*l*n +
               240*b*c*Pow(d,2)*e*Pow(f,2)*Pow(g,2)*k*l*n - 104*a*Pow(d,2)*Pow(e,2)*Pow(f,2)*Pow(g,2)*k*l*n - 32*b*Pow(c,2)*d*Pow(f,3)*Pow(g,2)*k*l*n -
               288*a*c*d*e*Pow(f,3)*Pow(g,2)*k*l*n + 104*a*Pow(c,2)*Pow(f,4)*Pow(g,2)*k*l*n + 128*Pow(b,2)*Pow(d,2)*Pow(f,2)*Pow(g,3)*k*l*n -
               128*a*b*d*Pow(f,3)*Pow(g,3)*k*l*n - 128*Pow(a,2)*Pow(f,4)*Pow(g,3)*k*l*n - 44*c*Pow(d,3)*Pow(e,4)*h*k*l*n + 100*Pow(c,2)*Pow(d,2)*Pow(e,3)*f*h*k*l*n -
               68*Pow(c,3)*d*Pow(e,2)*Pow(f,2)*h*k*l*n + 12*Pow(c,4)*e*Pow(f,3)*h*k*l*n - 128*b*Pow(d,3)*Pow(e,3)*g*h*k*l*n - 256*b*c*Pow(d,2)*Pow(e,2)*f*g*h*k*l*n -
               160*a*Pow(d,2)*Pow(e,3)*f*g*h*k*l*n - 256*b*Pow(c,2)*d*e*Pow(f,2)*g*h*k*l*n + 1088*a*c*d*Pow(e,2)*Pow(f,2)*g*h*k*l*n - 128*b*Pow(c,3)*Pow(f,3)*g*h*k*l*n -
               160*a*Pow(c,2)*e*Pow(f,3)*g*h*k*l*n - 1088*Pow(b,2)*Pow(d,2)*e*f*Pow(g,2)*h*k*l*n - 640*Pow(b,2)*c*d*Pow(f,2)*Pow(g,2)*h*k*l*n +
               1408*a*b*d*e*Pow(f,2)*Pow(g,2)*h*k*l*n + 768*a*b*c*Pow(f,3)*Pow(g,2)*h*k*l*n - 64*Pow(a,2)*e*Pow(f,3)*Pow(g,2)*h*k*l*n +
               416*b*c*Pow(d,2)*Pow(e,3)*Pow(h,2)*k*l*n - 232*a*Pow(d,2)*Pow(e,4)*Pow(h,2)*k*l*n - 208*b*Pow(c,2)*d*Pow(e,2)*f*Pow(h,2)*k*l*n -
               64*a*c*d*Pow(e,3)*f*Pow(h,2)*k*l*n + 80*b*Pow(c,3)*e*Pow(f,2)*Pow(h,2)*k*l*n + 8*a*Pow(c,2)*Pow(e,2)*Pow(f,2)*Pow(h,2)*k*l*n +
               896*Pow(b,2)*Pow(d,2)*Pow(e,2)*g*Pow(h,2)*k*l*n + 896*Pow(b,2)*c*d*e*f*g*Pow(h,2)*k*l*n - 896*a*b*d*Pow(e,2)*f*g*Pow(h,2)*k*l*n +
               896*Pow(b,2)*Pow(c,2)*Pow(f,2)*g*Pow(h,2)*k*l*n - 896*a*b*c*e*Pow(f,2)*g*Pow(h,2)*k*l*n - 896*Pow(a,2)*Pow(e,2)*Pow(f,2)*g*Pow(h,2)*k*l*n +
               2048*Pow(b,3)*d*f*Pow(g,2)*Pow(h,2)*k*l*n - 2304*a*Pow(b,2)*Pow(f,2)*Pow(g,2)*Pow(h,2)*k*l*n + 128*Pow(b,2)*c*d*Pow(e,2)*Pow(h,3)*k*l*n -
               512*a*b*d*Pow(e,3)*Pow(h,3)*k*l*n + 576*Pow(b,2)*Pow(c,2)*e*f*Pow(h,3)*k*l*n - 640*a*b*c*Pow(e,2)*f*Pow(h,3)*k*l*n +
               192*Pow(a,2)*Pow(e,3)*f*Pow(h,3)*k*l*n - 2560*Pow(b,3)*d*e*g*Pow(h,3)*k*l*n - 2560*Pow(b,3)*c*f*g*Pow(h,3)*k*l*n +
               3072*a*Pow(b,2)*e*f*g*Pow(h,3)*k*l*n - 1280*Pow(b,3)*c*e*Pow(h,4)*k*l*n + 1664*a*Pow(b,2)*Pow(e,2)*Pow(h,4)*k*l*n + 2560*Pow(b,4)*g*Pow(h,4)*k*l*n +
               12*Pow(c,2)*Pow(d,2)*Pow(e,4)*Pow(k,2)*l*n - 24*Pow(c,3)*d*Pow(e,3)*f*Pow(k,2)*l*n + 12*Pow(c,4)*Pow(e,2)*Pow(f,2)*Pow(k,2)*l*n -
               32*b*c*Pow(d,2)*Pow(e,3)*g*Pow(k,2)*l*n + 104*a*Pow(d,2)*Pow(e,4)*g*Pow(k,2)*l*n + 240*b*Pow(c,2)*d*Pow(e,2)*f*g*Pow(k,2)*l*n -
               288*a*c*d*Pow(e,3)*f*g*Pow(k,2)*l*n + 80*b*Pow(c,3)*e*Pow(f,2)*g*Pow(k,2)*l*n - 104*a*Pow(c,2)*Pow(e,2)*Pow(f,2)*g*Pow(k,2)*l*n -
               64*Pow(b,2)*Pow(d,2)*Pow(e,2)*Pow(g,2)*Pow(k,2)*l*n + 640*Pow(b,2)*c*d*e*f*Pow(g,2)*Pow(k,2)*l*n - 704*a*b*d*Pow(e,2)*f*Pow(g,2)*Pow(k,2)*l*n -
               64*Pow(b,2)*Pow(c,2)*Pow(f,2)*Pow(g,2)*Pow(k,2)*l*n - 704*a*b*c*e*Pow(f,2)*Pow(g,2)*Pow(k,2)*l*n +
               608*Pow(a,2)*Pow(e,2)*Pow(f,2)*Pow(g,2)*Pow(k,2)*l*n - 768*Pow(b,3)*d*f*Pow(g,3)*Pow(k,2)*l*n + 640*a*Pow(b,2)*Pow(f,2)*Pow(g,3)*Pow(k,2)*l*n -
               256*b*Pow(c,2)*d*Pow(e,3)*h*Pow(k,2)*l*n + 208*a*c*d*Pow(e,4)*h*Pow(k,2)*l*n - 32*b*Pow(c,3)*Pow(e,2)*f*h*Pow(k,2)*l*n +
               80*a*Pow(c,2)*Pow(e,3)*f*h*Pow(k,2)*l*n - 640*Pow(b,2)*c*d*Pow(e,2)*g*h*Pow(k,2)*l*n + 768*a*b*d*Pow(e,3)*g*h*Pow(k,2)*l*n -
               1088*Pow(b,2)*Pow(c,2)*e*f*g*h*Pow(k,2)*l*n + 1408*a*b*c*Pow(e,2)*f*g*h*Pow(k,2)*l*n - 64*Pow(a,2)*Pow(e,3)*f*g*h*Pow(k,2)*l*n +
               1024*Pow(b,3)*d*e*Pow(g,2)*h*Pow(k,2)*l*n + 1024*Pow(b,3)*c*f*Pow(g,2)*h*Pow(k,2)*l*n - 512*a*Pow(b,2)*e*f*Pow(g,2)*h*Pow(k,2)*l*n +
               128*Pow(b,2)*Pow(c,2)*Pow(e,2)*Pow(h,2)*Pow(k,2)*l*n - 128*a*b*c*Pow(e,3)*Pow(h,2)*Pow(k,2)*l*n + 96*Pow(a,2)*Pow(e,4)*Pow(h,2)*Pow(k,2)*l*n +
               2048*Pow(b,3)*c*e*g*Pow(h,2)*Pow(k,2)*l*n - 2304*a*Pow(b,2)*Pow(e,2)*g*Pow(h,2)*Pow(k,2)*l*n - 2048*Pow(b,4)*Pow(g,2)*Pow(h,2)*Pow(k,2)*l*n +
               64*b*Pow(c,3)*Pow(e,3)*Pow(k,3)*l*n - 64*a*Pow(c,2)*Pow(e,4)*Pow(k,3)*l*n + 128*Pow(b,2)*Pow(c,2)*Pow(e,2)*g*Pow(k,3)*l*n -
               128*a*b*c*Pow(e,3)*g*Pow(k,3)*l*n - 128*Pow(a,2)*Pow(e,4)*g*Pow(k,3)*l*n - 768*Pow(b,3)*c*e*Pow(g,2)*Pow(k,3)*l*n +
               640*a*Pow(b,2)*Pow(e,2)*Pow(g,2)*Pow(k,3)*l*n + 512*Pow(b,4)*Pow(g,3)*Pow(k,3)*l*n - 12*Pow(d,5)*Pow(e,2)*f*g*Pow(l,2)*n +
               24*c*Pow(d,4)*e*Pow(f,2)*g*Pow(l,2)*n - 12*Pow(c,2)*Pow(d,3)*Pow(f,3)*g*Pow(l,2)*n - 96*b*Pow(d,4)*Pow(f,2)*Pow(g,2)*Pow(l,2)*n +
               96*a*Pow(d,3)*Pow(f,3)*Pow(g,2)*Pow(l,2)*n - 24*Pow(d,5)*Pow(e,3)*h*Pow(l,2)*n + 88*c*Pow(d,4)*Pow(e,2)*f*h*Pow(l,2)*n -
               104*Pow(c,2)*Pow(d,3)*e*Pow(f,2)*h*Pow(l,2)*n + 40*Pow(c,3)*Pow(d,2)*Pow(f,3)*h*Pow(l,2)*n - 128*b*Pow(d,4)*e*f*g*h*Pow(l,2)*n +
               416*b*c*Pow(d,3)*Pow(f,2)*g*h*Pow(l,2)*n + 80*a*Pow(d,3)*e*Pow(f,2)*g*h*Pow(l,2)*n - 368*a*c*Pow(d,2)*Pow(f,3)*g*h*Pow(l,2)*n -
               144*b*Pow(d,4)*Pow(e,2)*Pow(h,2)*Pow(l,2)*n + 416*b*c*Pow(d,3)*e*f*Pow(h,2)*Pow(l,2)*n - 16*a*Pow(d,3)*Pow(e,2)*f*Pow(h,2)*Pow(l,2)*n -
               416*b*Pow(c,2)*Pow(d,2)*Pow(f,2)*Pow(h,2)*Pow(l,2)*n - 48*a*c*Pow(d,2)*e*Pow(f,2)*Pow(h,2)*Pow(l,2)*n +
               208*a*Pow(c,2)*d*Pow(f,3)*Pow(h,2)*Pow(l,2)*n - 64*Pow(b,2)*Pow(d,3)*f*g*Pow(h,2)*Pow(l,2)*n - 416*a*b*Pow(d,2)*Pow(f,2)*g*Pow(h,2)*Pow(l,2)*n +
               384*Pow(a,2)*d*Pow(f,3)*g*Pow(h,2)*Pow(l,2)*n + 768*Pow(b,2)*Pow(d,3)*e*Pow(h,3)*Pow(l,2)*n - 384*Pow(b,2)*c*Pow(d,2)*f*Pow(h,3)*Pow(l,2)*n -
               2112*a*b*Pow(d,2)*e*f*Pow(h,3)*Pow(l,2)*n + 1472*a*b*c*d*Pow(f,2)*Pow(h,3)*Pow(l,2)*n + 960*Pow(a,2)*d*e*Pow(f,2)*Pow(h,3)*Pow(l,2)*n -
               576*Pow(a,2)*c*Pow(f,3)*Pow(h,3)*Pow(l,2)*n - 768*Pow(b,3)*Pow(d,2)*Pow(h,4)*Pow(l,2)*n + 2560*a*Pow(b,2)*d*f*Pow(h,4)*Pow(l,2)*n -
               2304*Pow(a,2)*b*Pow(f,2)*Pow(h,4)*Pow(l,2)*n + 20*c*Pow(d,4)*Pow(e,3)*k*Pow(l,2)*n - 76*Pow(c,2)*Pow(d,3)*Pow(e,2)*f*k*Pow(l,2)*n +
               92*Pow(c,3)*Pow(d,2)*e*Pow(f,2)*k*Pow(l,2)*n - 36*Pow(c,4)*d*Pow(f,3)*k*Pow(l,2)*n + 104*b*Pow(d,4)*Pow(e,2)*g*k*Pow(l,2)*n -
               288*b*c*Pow(d,3)*e*f*g*k*Pow(l,2)*n - 32*a*Pow(d,3)*Pow(e,2)*f*g*k*Pow(l,2)*n - 104*b*Pow(c,2)*Pow(d,2)*Pow(f,2)*g*k*Pow(l,2)*n +
               240*a*c*Pow(d,2)*e*Pow(f,2)*g*k*Pow(l,2)*n + 80*a*Pow(c,2)*d*Pow(f,3)*g*k*Pow(l,2)*n + 128*Pow(b,2)*Pow(d,3)*f*Pow(g,2)*k*Pow(l,2)*n -
               64*a*b*Pow(d,2)*Pow(f,2)*Pow(g,2)*k*Pow(l,2)*n + 128*Pow(a,2)*d*Pow(f,3)*Pow(g,2)*k*Pow(l,2)*n + 32*b*c*Pow(d,3)*Pow(e,2)*h*k*Pow(l,2)*n +
               160*a*Pow(d,3)*Pow(e,3)*h*k*Pow(l,2)*n - 16*b*Pow(c,2)*Pow(d,2)*e*f*h*k*Pow(l,2)*n - 128*a*c*Pow(d,2)*Pow(e,2)*f*h*k*Pow(l,2)*n +
               368*b*Pow(c,3)*d*Pow(f,2)*h*k*Pow(l,2)*n - 464*a*Pow(c,2)*d*e*Pow(f,2)*h*k*Pow(l,2)*n + 48*a*Pow(c,3)*Pow(f,3)*h*k*Pow(l,2)*n -
               64*Pow(b,2)*Pow(d,3)*e*g*h*k*Pow(l,2)*n + 704*Pow(b,2)*c*Pow(d,2)*f*g*h*k*Pow(l,2)*n + 1408*a*b*Pow(d,2)*e*f*g*h*k*Pow(l,2)*n -
               1280*a*b*c*d*Pow(f,2)*g*h*k*Pow(l,2)*n - 1088*Pow(a,2)*d*e*Pow(f,2)*g*h*k*Pow(l,2)*n - 64*Pow(a,2)*c*Pow(f,3)*g*h*k*Pow(l,2)*n -
               1792*Pow(b,2)*c*Pow(d,2)*e*Pow(h,2)*k*Pow(l,2)*n + 896*a*b*Pow(d,2)*Pow(e,2)*Pow(h,2)*k*Pow(l,2)*n -
               448*Pow(b,2)*Pow(c,2)*d*f*Pow(h,2)*k*Pow(l,2)*n + 1792*a*b*c*d*e*f*Pow(h,2)*k*Pow(l,2)*n - 448*Pow(a,2)*d*Pow(e,2)*f*Pow(h,2)*k*Pow(l,2)*n -
               896*a*b*Pow(c,2)*Pow(f,2)*Pow(h,2)*k*Pow(l,2)*n + 896*Pow(a,2)*c*e*Pow(f,2)*Pow(h,2)*k*Pow(l,2)*n +
               128*Pow(b,3)*Pow(d,2)*g*Pow(h,2)*k*Pow(l,2)*n - 2816*a*Pow(b,2)*d*f*g*Pow(h,2)*k*Pow(l,2)*n + 3968*Pow(a,2)*b*Pow(f,2)*g*Pow(h,2)*k*Pow(l,2)*n +
               2304*Pow(b,3)*c*d*Pow(h,3)*k*Pow(l,2)*n - 256*a*Pow(b,2)*d*e*Pow(h,3)*k*Pow(l,2)*n - 256*a*Pow(b,2)*c*f*Pow(h,3)*k*Pow(l,2)*n -
               256*Pow(a,2)*b*e*f*Pow(h,3)*k*Pow(l,2)*n - 2048*a*Pow(b,3)*Pow(h,4)*k*Pow(l,2)*n + 32*b*Pow(c,2)*Pow(d,2)*Pow(e,2)*Pow(k,2)*Pow(l,2)*n -
               128*a*c*Pow(d,2)*Pow(e,3)*Pow(k,2)*Pow(l,2)*n - 32*b*Pow(c,3)*d*e*f*Pow(k,2)*Pow(l,2)*n + 176*a*Pow(c,2)*d*Pow(e,2)*f*Pow(k,2)*Pow(l,2)*n -
               144*b*Pow(c,4)*Pow(f,2)*Pow(k,2)*Pow(l,2)*n + 96*a*Pow(c,3)*e*Pow(f,2)*Pow(k,2)*Pow(l,2)*n + 320*Pow(b,2)*c*Pow(d,2)*e*g*Pow(k,2)*Pow(l,2)*n -
               576*a*b*Pow(d,2)*Pow(e,2)*g*Pow(k,2)*Pow(l,2)*n - 576*Pow(b,2)*Pow(c,2)*d*f*g*Pow(k,2)*Pow(l,2)*n + 384*a*b*c*d*e*f*g*Pow(k,2)*Pow(l,2)*n +
               320*Pow(a,2)*d*Pow(e,2)*f*g*Pow(k,2)*Pow(l,2)*n + 992*a*b*Pow(c,2)*Pow(f,2)*g*Pow(k,2)*Pow(l,2)*n -
               576*Pow(a,2)*c*e*Pow(f,2)*g*Pow(k,2)*Pow(l,2)*n - 128*Pow(b,3)*Pow(d,2)*Pow(g,2)*Pow(k,2)*Pow(l,2)*n +
               1024*a*Pow(b,2)*d*f*Pow(g,2)*Pow(k,2)*Pow(l,2)*n - 1280*Pow(a,2)*b*Pow(f,2)*Pow(g,2)*Pow(k,2)*Pow(l,2)*n +
               1024*Pow(b,2)*Pow(c,2)*d*e*h*Pow(k,2)*Pow(l,2)*n - 640*a*b*c*d*Pow(e,2)*h*Pow(k,2)*Pow(l,2)*n - 256*Pow(a,2)*d*Pow(e,3)*h*Pow(k,2)*Pow(l,2)*n +
               640*Pow(b,2)*Pow(c,3)*f*h*Pow(k,2)*Pow(l,2)*n - 1088*a*b*Pow(c,2)*e*f*h*Pow(k,2)*Pow(l,2)*n + 128*Pow(a,2)*c*Pow(e,2)*f*h*Pow(k,2)*Pow(l,2)*n -
               768*Pow(b,3)*c*d*g*h*Pow(k,2)*Pow(l,2)*n - 512*a*Pow(b,2)*d*e*g*h*Pow(k,2)*Pow(l,2)*n - 512*a*Pow(b,2)*c*f*g*h*Pow(k,2)*Pow(l,2)*n -
               512*Pow(a,2)*b*e*f*g*h*Pow(k,2)*Pow(l,2)*n - 1280*Pow(b,3)*Pow(c,2)*Pow(h,2)*Pow(k,2)*Pow(l,2)*n +
               1280*a*Pow(b,2)*c*e*Pow(h,2)*Pow(k,2)*Pow(l,2)*n - 256*Pow(a,2)*b*Pow(e,2)*Pow(h,2)*Pow(k,2)*Pow(l,2)*n +
               2560*a*Pow(b,3)*g*Pow(h,2)*Pow(k,2)*Pow(l,2)*n - 256*Pow(b,2)*Pow(c,3)*e*Pow(k,3)*Pow(l,2)*n + 128*a*b*Pow(c,2)*Pow(e,2)*Pow(k,3)*Pow(l,2)*n +
               192*Pow(a,2)*c*Pow(e,3)*Pow(k,3)*Pow(l,2)*n + 512*Pow(b,3)*Pow(c,2)*g*Pow(k,3)*Pow(l,2)*n - 512*a*Pow(b,2)*c*e*g*Pow(k,3)*Pow(l,2)*n +
               640*Pow(a,2)*b*Pow(e,2)*g*Pow(k,3)*Pow(l,2)*n - 512*a*Pow(b,3)*Pow(g,2)*Pow(k,3)*Pow(l,2)*n + 4*Pow(d,6)*Pow(e,2)*Pow(l,3)*n -
               8*c*Pow(d,5)*e*f*Pow(l,3)*n + 4*Pow(c,2)*Pow(d,4)*Pow(f,2)*Pow(l,3)*n + 64*b*Pow(d,5)*f*g*Pow(l,3)*n - 64*a*Pow(d,4)*Pow(f,2)*g*Pow(l,3)*n +
               96*b*Pow(d,5)*e*h*Pow(l,3)*n - 192*b*c*Pow(d,4)*f*h*Pow(l,3)*n - 80*a*Pow(d,4)*e*f*h*Pow(l,3)*n + 176*a*c*Pow(d,3)*Pow(f,2)*h*Pow(l,3)*n -
               192*Pow(b,2)*Pow(d,4)*Pow(h,2)*Pow(l,3)*n + 704*a*b*Pow(d,3)*f*Pow(h,2)*Pow(l,3)*n - 480*Pow(a,2)*Pow(d,2)*Pow(f,2)*Pow(h,2)*Pow(l,3)*n -
               80*b*c*Pow(d,4)*e*k*Pow(l,3)*n - 40*a*Pow(d,4)*Pow(e,2)*k*Pow(l,3)*n + 176*b*Pow(c,2)*Pow(d,3)*f*k*Pow(l,3)*n + 128*a*c*Pow(d,3)*e*f*k*Pow(l,3)*n -
               184*a*Pow(c,2)*Pow(d,2)*Pow(f,2)*k*Pow(l,3)*n - 128*Pow(b,2)*Pow(d,4)*g*k*Pow(l,3)*n - 128*a*b*Pow(d,3)*f*g*k*Pow(l,3)*n +
               128*Pow(a,2)*Pow(d,2)*Pow(f,2)*g*k*Pow(l,3)*n + 576*Pow(b,2)*c*Pow(d,3)*h*k*Pow(l,3)*n - 640*a*b*Pow(d,3)*e*h*k*Pow(l,3)*n -
               1024*a*b*c*Pow(d,2)*f*h*k*Pow(l,3)*n + 384*Pow(a,2)*Pow(d,2)*e*f*h*k*Pow(l,3)*n + 832*Pow(a,2)*c*d*Pow(f,2)*h*k*Pow(l,3)*n +
               128*a*Pow(b,2)*Pow(d,2)*Pow(h,2)*k*Pow(l,3)*n + 256*Pow(a,2)*b*d*f*Pow(h,2)*k*Pow(l,3)*n - 1152*Pow(a,3)*Pow(f,2)*Pow(h,2)*k*Pow(l,3)*n -
               320*Pow(b,2)*Pow(c,2)*Pow(d,2)*Pow(k,2)*Pow(l,3)*n + 512*a*b*c*Pow(d,2)*e*Pow(k,2)*Pow(l,3)*n +
               128*Pow(a,2)*Pow(d,2)*Pow(e,2)*Pow(k,2)*Pow(l,3)*n + 64*a*b*Pow(c,2)*d*f*Pow(k,2)*Pow(l,3)*n - 384*Pow(a,2)*c*d*e*f*Pow(k,2)*Pow(l,3)*n -
               96*Pow(a,2)*Pow(c,2)*Pow(f,2)*Pow(k,2)*Pow(l,3)*n + 640*a*Pow(b,2)*Pow(d,2)*g*Pow(k,2)*Pow(l,3)*n - 512*Pow(a,2)*b*d*f*g*Pow(k,2)*Pow(l,3)*n +
               512*Pow(a,3)*Pow(f,2)*g*Pow(k,2)*Pow(l,3)*n - 768*a*Pow(b,2)*c*d*h*Pow(k,2)*Pow(l,3)*n + 1024*Pow(a,2)*b*d*e*h*Pow(k,2)*Pow(l,3)*n +
               1024*Pow(a,2)*b*c*f*h*Pow(k,2)*Pow(l,3)*n - 256*Pow(a,3)*e*f*h*Pow(k,2)*Pow(l,3)*n - 512*Pow(a,2)*Pow(b,2)*Pow(h,2)*Pow(k,2)*Pow(l,3)*n +
               512*a*Pow(b,2)*Pow(c,2)*Pow(k,3)*Pow(l,3)*n - 768*Pow(a,2)*b*c*e*Pow(k,3)*Pow(l,3)*n - 128*Pow(a,3)*Pow(e,2)*Pow(k,3)*Pow(l,3)*n -
               512*Pow(a,2)*Pow(b,2)*g*Pow(k,3)*Pow(l,3)*n - 16*b*Pow(d,6)*Pow(l,4)*n + 16*a*Pow(d,5)*f*Pow(l,4)*n + 160*a*b*Pow(d,4)*k*Pow(l,4)*n -
               128*Pow(a,2)*Pow(d,3)*f*k*Pow(l,4)*n - 512*Pow(a,2)*b*Pow(d,2)*Pow(k,2)*Pow(l,4)*n + 256*Pow(a,3)*d*f*Pow(k,2)*Pow(l,4)*n +
               512*Pow(a,3)*b*Pow(k,3)*Pow(l,4)*n - 36*Pow(d,4)*Pow(e,3)*f*Pow(g,2)*m*n + 92*c*Pow(d,3)*Pow(e,2)*Pow(f,2)*Pow(g,2)*m*n -
               76*Pow(c,2)*Pow(d,2)*e*Pow(f,3)*Pow(g,2)*m*n + 20*Pow(c,3)*d*Pow(f,4)*Pow(g,2)*m*n - 160*b*Pow(d,3)*e*Pow(f,2)*Pow(g,3)*m*n +
               64*b*c*Pow(d,2)*Pow(f,3)*Pow(g,3)*m*n + 176*a*Pow(d,2)*e*Pow(f,3)*Pow(g,3)*m*n - 80*a*c*d*Pow(f,4)*Pow(g,3)*m*n - 72*Pow(d,4)*Pow(e,4)*g*h*m*n +
               248*c*Pow(d,3)*Pow(e,3)*f*g*h*m*n - 296*Pow(c,2)*Pow(d,2)*Pow(e,2)*Pow(f,2)*g*h*m*n + 136*Pow(c,3)*d*e*Pow(f,3)*g*h*m*n -
               16*Pow(c,4)*Pow(f,4)*g*h*m*n - 256*b*Pow(d,3)*Pow(e,2)*f*Pow(g,2)*h*m*n + 800*b*c*Pow(d,2)*e*Pow(f,2)*Pow(g,2)*h*m*n -
               48*a*Pow(d,2)*Pow(e,2)*Pow(f,2)*Pow(g,2)*h*m*n - 256*b*Pow(c,2)*d*Pow(f,3)*Pow(g,2)*h*m*n - 288*a*c*d*e*Pow(f,3)*Pow(g,2)*h*m*n +
               48*a*Pow(c,2)*Pow(f,4)*Pow(g,2)*h*m*n - 512*Pow(b,2)*Pow(d,2)*Pow(f,2)*Pow(g,3)*h*m*n + 512*a*b*d*Pow(f,3)*Pow(g,3)*h*m*n +
               64*Pow(a,2)*Pow(f,4)*Pow(g,3)*h*m*n + 80*c*Pow(d,3)*Pow(e,4)*Pow(h,2)*m*n - 192*Pow(c,2)*Pow(d,2)*Pow(e,3)*f*Pow(h,2)*m*n +
               144*Pow(c,3)*d*Pow(e,2)*Pow(f,2)*Pow(h,2)*m*n - 32*Pow(c,4)*e*Pow(f,3)*Pow(h,2)*m*n + 368*b*Pow(d,3)*Pow(e,3)*g*Pow(h,2)*m*n +
               400*b*c*Pow(d,2)*Pow(e,2)*f*g*Pow(h,2)*m*n - 800*a*Pow(d,2)*Pow(e,3)*f*g*Pow(h,2)*m*n - 608*b*Pow(c,2)*d*e*Pow(f,2)*g*Pow(h,2)*m*n +
               400*a*c*d*Pow(e,2)*Pow(f,2)*g*Pow(h,2)*m*n + 32*b*Pow(c,3)*Pow(f,3)*g*Pow(h,2)*m*n + 208*a*Pow(c,2)*e*Pow(f,3)*g*Pow(h,2)*m*n +
               192*Pow(b,2)*Pow(d,2)*e*f*Pow(g,2)*Pow(h,2)*m*n + 640*Pow(b,2)*c*d*Pow(f,2)*Pow(g,2)*Pow(h,2)*m*n - 512*a*b*d*e*Pow(f,2)*Pow(g,2)*Pow(h,2)*m*n +
               128*a*b*c*Pow(f,3)*Pow(g,2)*Pow(h,2)*m*n - 832*Pow(a,2)*e*Pow(f,3)*Pow(g,2)*Pow(h,2)*m*n - 736*b*c*Pow(d,2)*Pow(e,3)*Pow(h,3)*m*n +
               96*a*Pow(d,2)*Pow(e,4)*Pow(h,3)*m*n + 32*b*Pow(c,2)*d*Pow(e,2)*f*Pow(h,3)*m*n + 992*a*c*d*Pow(e,3)*f*Pow(h,3)*m*n -
               64*b*Pow(c,3)*e*Pow(f,2)*Pow(h,3)*m*n - 320*a*Pow(c,2)*Pow(e,2)*Pow(f,2)*Pow(h,3)*m*n - 64*Pow(b,2)*Pow(d,2)*Pow(e,2)*g*Pow(h,3)*m*n -
               1152*Pow(b,2)*c*d*e*f*g*Pow(h,3)*m*n - 256*a*b*d*Pow(e,2)*f*g*Pow(h,3)*m*n + 384*Pow(b,2)*Pow(c,2)*Pow(f,2)*g*Pow(h,3)*m*n +
               640*a*b*c*e*Pow(f,2)*g*Pow(h,3)*m*n + 832*Pow(a,2)*Pow(e,2)*Pow(f,2)*g*Pow(h,3)*m*n + 1024*Pow(b,3)*d*f*Pow(g,2)*Pow(h,3)*m*n -
               2048*a*Pow(b,2)*Pow(f,2)*Pow(g,2)*Pow(h,3)*m*n + 1408*Pow(b,2)*c*d*Pow(e,2)*Pow(h,4)*m*n + 192*a*b*d*Pow(e,3)*Pow(h,4)*m*n +
               512*Pow(b,2)*Pow(c,2)*e*f*Pow(h,4)*m*n - 1216*a*b*c*Pow(e,2)*f*Pow(h,4)*m*n - 576*Pow(a,2)*Pow(e,3)*f*Pow(h,4)*m*n -
               1280*Pow(b,3)*d*e*g*Pow(h,4)*m*n - 1280*Pow(b,3)*c*f*g*Pow(h,4)*m*n + 3328*a*Pow(b,2)*e*f*g*Pow(h,4)*m*n - 512*Pow(b,3)*c*e*Pow(h,5)*m*n -
               768*a*Pow(b,2)*Pow(e,2)*Pow(h,5)*m*n + 1024*Pow(b,4)*g*Pow(h,5)*m*n + 12*c*Pow(d,3)*Pow(e,4)*g*k*m*n - 68*Pow(c,2)*Pow(d,2)*Pow(e,3)*f*g*k*m*n +
               100*Pow(c,3)*d*Pow(e,2)*Pow(f,2)*g*k*m*n - 44*Pow(c,4)*e*Pow(f,3)*g*k*m*n + 48*b*Pow(d,3)*Pow(e,3)*Pow(g,2)*k*m*n -
               464*b*c*Pow(d,2)*Pow(e,2)*f*Pow(g,2)*k*m*n + 368*a*Pow(d,2)*Pow(e,3)*f*Pow(g,2)*k*m*n - 128*b*Pow(c,2)*d*e*Pow(f,2)*Pow(g,2)*k*m*n -
               16*a*c*d*Pow(e,2)*Pow(f,2)*Pow(g,2)*k*m*n + 160*b*Pow(c,3)*Pow(f,3)*Pow(g,2)*k*m*n + 32*a*Pow(c,2)*e*Pow(f,3)*Pow(g,2)*k*m*n +
               832*Pow(b,2)*Pow(d,2)*e*f*Pow(g,3)*k*m*n + 384*Pow(b,2)*c*d*Pow(f,2)*Pow(g,3)*k*m*n - 1024*a*b*d*e*Pow(f,2)*Pow(g,3)*k*m*n -
               640*a*b*c*Pow(f,3)*Pow(g,3)*k*m*n + 576*Pow(a,2)*e*Pow(f,3)*Pow(g,3)*k*m*n - 16*Pow(c,2)*Pow(d,2)*Pow(e,4)*h*k*m*n +
               32*Pow(c,3)*d*Pow(e,3)*f*h*k*m*n - 16*Pow(c,4)*Pow(e,2)*Pow(f,2)*h*k*m*n - 64*b*c*Pow(d,2)*Pow(e,3)*g*h*k*m*n + 208*a*Pow(d,2)*Pow(e,4)*g*h*k*m*n +
               480*b*Pow(c,2)*d*Pow(e,2)*f*g*h*k*m*n - 576*a*c*d*Pow(e,3)*f*g*h*k*m*n + 160*b*Pow(c,3)*e*Pow(f,2)*g*h*k*m*n -
               208*a*Pow(c,2)*Pow(e,2)*Pow(f,2)*g*h*k*m*n - 448*Pow(b,2)*Pow(d,2)*Pow(e,2)*Pow(g,2)*h*k*m*n - 896*Pow(b,2)*c*d*e*f*Pow(g,2)*h*k*m*n +
               1792*a*b*d*Pow(e,2)*f*Pow(g,2)*h*k*m*n - 896*Pow(b,2)*Pow(c,2)*Pow(f,2)*Pow(g,2)*h*k*m*n + 896*a*b*c*e*Pow(f,2)*Pow(g,2)*h*k*m*n -
               448*Pow(a,2)*Pow(e,2)*Pow(f,2)*Pow(g,2)*h*k*m*n - 1024*Pow(b,3)*d*f*Pow(g,3)*h*k*m*n + 2048*a*Pow(b,2)*Pow(f,2)*Pow(g,3)*h*k*m*n +
               400*b*Pow(c,2)*d*Pow(e,3)*Pow(h,2)*k*m*n - 304*a*c*d*Pow(e,4)*Pow(h,2)*k*m*n + 176*b*Pow(c,3)*Pow(e,2)*f*Pow(h,2)*k*m*n -
               272*a*Pow(c,2)*Pow(e,3)*f*Pow(h,2)*k*m*n - 512*Pow(b,2)*c*d*Pow(e,2)*g*Pow(h,2)*k*m*n - 640*a*b*d*Pow(e,3)*g*Pow(h,2)*k*m*n -
               64*Pow(b,2)*Pow(c,2)*e*f*g*Pow(h,2)*k*m*n - 1024*a*b*c*Pow(e,2)*f*g*Pow(h,2)*k*m*n + 1472*Pow(a,2)*Pow(e,3)*f*g*Pow(h,2)*k*m*n +
               1536*Pow(b,3)*d*e*Pow(g,2)*Pow(h,2)*k*m*n + 1536*Pow(b,3)*c*f*Pow(g,2)*Pow(h,2)*k*m*n - 2560*a*Pow(b,2)*e*f*Pow(g,2)*Pow(h,2)*k*m*n -
               1216*Pow(b,2)*Pow(c,2)*Pow(e,2)*Pow(h,3)*k*m*n + 1664*a*b*c*Pow(e,3)*Pow(h,3)*k*m*n - 576*Pow(a,2)*Pow(e,4)*Pow(h,3)*k*m*n +
               2048*Pow(b,3)*c*e*g*Pow(h,3)*k*m*n - 512*a*Pow(b,2)*Pow(e,2)*g*Pow(h,3)*k*m*n - 2048*Pow(b,4)*Pow(g,2)*Pow(h,3)*k*m*n +
               80*b*Pow(c,2)*d*Pow(e,3)*g*Pow(k,2)*m*n - 128*a*c*d*Pow(e,4)*g*Pow(k,2)*m*n - 368*b*Pow(c,3)*Pow(e,2)*f*g*Pow(k,2)*m*n +
               416*a*Pow(c,2)*Pow(e,3)*f*g*Pow(k,2)*m*n + 128*Pow(b,2)*c*d*Pow(e,2)*Pow(g,2)*Pow(k,2)*m*n - 64*a*b*d*Pow(e,3)*Pow(g,2)*Pow(k,2)*m*n +
               576*Pow(b,2)*Pow(c,2)*e*f*Pow(g,2)*Pow(k,2)*m*n + 704*a*b*c*Pow(e,2)*f*Pow(g,2)*Pow(k,2)*m*n - 1152*Pow(a,2)*Pow(e,3)*f*Pow(g,2)*Pow(k,2)*m*n -
               256*Pow(b,3)*d*e*Pow(g,3)*Pow(k,2)*m*n - 256*Pow(b,3)*c*f*Pow(g,3)*Pow(k,2)*m*n - 768*a*Pow(b,2)*e*f*Pow(g,3)*Pow(k,2)*m*n -
               96*b*Pow(c,3)*Pow(e,3)*h*Pow(k,2)*m*n + 96*a*Pow(c,2)*Pow(e,4)*h*Pow(k,2)*m*n + 704*Pow(b,2)*Pow(c,2)*Pow(e,2)*g*h*Pow(k,2)*m*n -
               1152*a*b*c*Pow(e,3)*g*h*Pow(k,2)*m*n + 640*Pow(a,2)*Pow(e,4)*g*h*Pow(k,2)*m*n - 1536*Pow(b,3)*c*e*Pow(g,2)*h*Pow(k,2)*m*n +
               1280*a*Pow(b,2)*Pow(e,2)*Pow(g,2)*h*Pow(k,2)*m*n + 1024*Pow(b,4)*Pow(g,3)*h*Pow(k,2)*m*n + 36*Pow(d,5)*Pow(e,3)*g*l*m*n -
               76*c*Pow(d,4)*Pow(e,2)*f*g*l*m*n + 44*Pow(c,2)*Pow(d,3)*e*Pow(f,2)*g*l*m*n - 4*Pow(c,3)*Pow(d,2)*Pow(f,3)*g*l*m*n +
               320*b*Pow(d,4)*e*f*Pow(g,2)*l*m*n - 32*b*c*Pow(d,3)*Pow(f,2)*Pow(g,2)*l*m*n - 368*a*Pow(d,3)*e*Pow(f,2)*Pow(g,2)*l*m*n +
               80*a*c*Pow(d,2)*Pow(f,3)*Pow(g,2)*l*m*n - 8*c*Pow(d,4)*Pow(e,3)*h*l*m*n + 8*Pow(c,2)*Pow(d,3)*Pow(e,2)*f*h*l*m*n +
               8*Pow(c,3)*Pow(d,2)*e*Pow(f,2)*h*l*m*n - 8*Pow(c,4)*d*Pow(f,3)*h*l*m*n + 208*b*Pow(d,4)*Pow(e,2)*g*h*l*m*n - 576*b*c*Pow(d,3)*e*f*g*h*l*m*n -
               64*a*Pow(d,3)*Pow(e,2)*f*g*h*l*m*n - 208*b*Pow(c,2)*Pow(d,2)*Pow(f,2)*g*h*l*m*n + 480*a*c*Pow(d,2)*e*Pow(f,2)*g*h*l*m*n +
               160*a*Pow(c,2)*d*Pow(f,3)*g*h*l*m*n - 128*Pow(b,2)*Pow(d,3)*f*Pow(g,2)*h*l*m*n + 960*a*b*Pow(d,2)*Pow(f,2)*Pow(g,2)*h*l*m*n -
               1024*Pow(a,2)*d*Pow(f,3)*Pow(g,2)*h*l*m*n + 96*b*c*Pow(d,3)*Pow(e,2)*Pow(h,2)*l*m*n - 304*a*Pow(d,3)*Pow(e,3)*Pow(h,2)*l*m*n -
               384*b*Pow(c,2)*Pow(d,2)*e*f*Pow(h,2)*l*m*n + 400*a*c*Pow(d,2)*Pow(e,2)*f*Pow(h,2)*l*m*n + 96*b*Pow(c,3)*d*Pow(f,2)*Pow(h,2)*l*m*n +
               400*a*Pow(c,2)*d*e*Pow(f,2)*Pow(h,2)*l*m*n - 304*a*Pow(c,3)*Pow(f,3)*Pow(h,2)*l*m*n - 2112*Pow(b,2)*Pow(d,3)*e*g*Pow(h,2)*l*m*n -
               64*Pow(b,2)*c*Pow(d,2)*f*g*Pow(h,2)*l*m*n + 4352*a*b*Pow(d,2)*e*f*g*Pow(h,2)*l*m*n - 1024*a*b*c*d*Pow(f,2)*g*Pow(h,2)*l*m*n -
               1856*Pow(a,2)*d*e*Pow(f,2)*g*Pow(h,2)*l*m*n + 1472*Pow(a,2)*c*Pow(f,3)*g*Pow(h,2)*l*m*n + 1664*Pow(b,2)*c*Pow(d,2)*e*Pow(h,3)*l*m*n +
               320*a*b*Pow(d,2)*Pow(e,2)*Pow(h,3)*l*m*n + 1664*Pow(b,2)*Pow(c,2)*d*f*Pow(h,3)*l*m*n - 2304*a*b*c*d*e*f*Pow(h,3)*l*m*n -
               1024*Pow(a,2)*d*Pow(e,2)*f*Pow(h,3)*l*m*n + 320*a*b*Pow(c,2)*Pow(f,2)*Pow(h,3)*l*m*n - 1024*Pow(a,2)*c*e*Pow(f,2)*Pow(h,3)*l*m*n +
               2816*Pow(b,3)*Pow(d,2)*g*Pow(h,3)*l*m*n - 4608*a*Pow(b,2)*d*f*g*Pow(h,3)*l*m*n + 1280*Pow(a,2)*b*Pow(f,2)*g*Pow(h,3)*l*m*n -
               3072*Pow(b,3)*c*d*Pow(h,4)*l*m*n - 2048*a*Pow(b,2)*d*e*Pow(h,4)*l*m*n - 2048*a*Pow(b,2)*c*f*Pow(h,4)*l*m*n + 5120*Pow(a,2)*b*e*f*Pow(h,4)*l*m*n +
               4096*a*Pow(b,3)*Pow(h,5)*l*m*n - 4*Pow(c,2)*Pow(d,3)*Pow(e,3)*k*l*m*n + 44*Pow(c,3)*Pow(d,2)*Pow(e,2)*f*k*l*m*n - 76*Pow(c,4)*d*e*Pow(f,2)*k*l*m*n +
               36*Pow(c,5)*Pow(f,3)*k*l*m*n - 160*b*c*Pow(d,3)*Pow(e,2)*g*k*l*m*n - 128*a*Pow(d,3)*Pow(e,3)*g*k*l*m*n + 1088*b*Pow(c,2)*Pow(d,2)*e*f*g*k*l*m*n -
               256*a*c*Pow(d,2)*Pow(e,2)*f*g*k*l*m*n - 160*b*Pow(c,3)*d*Pow(f,2)*g*k*l*m*n - 256*a*Pow(c,2)*d*e*Pow(f,2)*g*k*l*m*n - 128*a*Pow(c,3)*Pow(f,3)*g*k*l*m*n -
               64*Pow(b,2)*Pow(d,3)*e*Pow(g,2)*k*l*m*n - 1088*Pow(b,2)*c*Pow(d,2)*f*Pow(g,2)*k*l*m*n - 1280*a*b*Pow(d,2)*e*f*Pow(g,2)*k*l*m*n +
               1408*a*b*c*d*Pow(f,2)*Pow(g,2)*k*l*m*n + 704*Pow(a,2)*d*e*Pow(f,2)*Pow(g,2)*k*l*m*n - 64*Pow(a,2)*c*Pow(f,3)*Pow(g,2)*k*l*m*n -
               208*b*Pow(c,2)*Pow(d,2)*Pow(e,2)*h*k*l*m*n + 160*a*c*Pow(d,2)*Pow(e,3)*h*k*l*m*n - 576*b*Pow(c,3)*d*e*f*h*k*l*m*n +
               480*a*Pow(c,2)*d*Pow(e,2)*f*h*k*l*m*n + 208*b*Pow(c,4)*Pow(f,2)*h*k*l*m*n - 64*a*Pow(c,3)*e*Pow(f,2)*h*k*l*m*n + 1792*Pow(b,2)*c*Pow(d,2)*e*g*h*k*l*m*n -
               896*a*b*Pow(d,2)*Pow(e,2)*g*h*k*l*m*n + 1792*Pow(b,2)*Pow(c,2)*d*f*g*h*k*l*m*n - 5376*a*b*c*d*e*f*g*h*k*l*m*n + 1792*Pow(a,2)*d*Pow(e,2)*f*g*h*k*l*m*n -
               896*a*b*Pow(c,2)*Pow(f,2)*g*h*k*l*m*n + 1792*Pow(a,2)*c*e*Pow(f,2)*g*h*k*l*m*n + 256*Pow(b,3)*Pow(d,2)*Pow(g,2)*h*k*l*m*n +
               1536*a*Pow(b,2)*d*f*Pow(g,2)*h*k*l*m*n - 2816*Pow(a,2)*b*Pow(f,2)*Pow(g,2)*h*k*l*m*n - 64*Pow(b,2)*Pow(c,2)*d*e*Pow(h,2)*k*l*m*n -
               1024*a*b*c*d*Pow(e,2)*Pow(h,2)*k*l*m*n + 1472*Pow(a,2)*d*Pow(e,3)*Pow(h,2)*k*l*m*n - 2112*Pow(b,2)*Pow(c,3)*f*Pow(h,2)*k*l*m*n +
               4352*a*b*Pow(c,2)*e*f*Pow(h,2)*k*l*m*n - 1856*Pow(a,2)*c*Pow(e,2)*f*Pow(h,2)*k*l*m*n - 4096*Pow(b,3)*c*d*g*Pow(h,2)*k*l*m*n +
               5632*a*Pow(b,2)*d*e*g*Pow(h,2)*k*l*m*n + 5632*a*Pow(b,2)*c*f*g*Pow(h,2)*k*l*m*n - 5120*Pow(a,2)*b*e*f*g*Pow(h,2)*k*l*m*n +
               2816*Pow(b,3)*Pow(c,2)*Pow(h,3)*k*l*m*n - 4608*a*Pow(b,2)*c*e*Pow(h,3)*k*l*m*n + 1280*Pow(a,2)*b*Pow(e,2)*Pow(h,3)*k*l*m*n -
               2048*a*Pow(b,3)*g*Pow(h,3)*k*l*m*n - 32*b*Pow(c,3)*d*Pow(e,2)*Pow(k,2)*l*m*n + 80*a*Pow(c,2)*d*Pow(e,3)*Pow(k,2)*l*m*n +
               320*b*Pow(c,4)*e*f*Pow(k,2)*l*m*n - 368*a*Pow(c,3)*Pow(e,2)*f*Pow(k,2)*l*m*n - 1088*Pow(b,2)*Pow(c,2)*d*e*g*Pow(k,2)*l*m*n +
               1408*a*b*c*d*Pow(e,2)*g*Pow(k,2)*l*m*n - 64*Pow(a,2)*d*Pow(e,3)*g*Pow(k,2)*l*m*n - 64*Pow(b,2)*Pow(c,3)*f*g*Pow(k,2)*l*m*n -
               1280*a*b*Pow(c,2)*e*f*g*Pow(k,2)*l*m*n + 704*Pow(a,2)*c*Pow(e,2)*f*g*Pow(k,2)*l*m*n + 1024*Pow(b,3)*c*d*Pow(g,2)*Pow(k,2)*l*m*n -
               512*a*Pow(b,2)*d*e*Pow(g,2)*Pow(k,2)*l*m*n - 512*a*Pow(b,2)*c*f*Pow(g,2)*Pow(k,2)*l*m*n + 3072*Pow(a,2)*b*e*f*Pow(g,2)*Pow(k,2)*l*m*n -
               128*Pow(b,2)*Pow(c,3)*e*h*Pow(k,2)*l*m*n + 960*a*b*Pow(c,2)*Pow(e,2)*h*Pow(k,2)*l*m*n - 1024*Pow(a,2)*c*Pow(e,3)*h*Pow(k,2)*l*m*n +
               256*Pow(b,3)*Pow(c,2)*g*h*Pow(k,2)*l*m*n + 1536*a*Pow(b,2)*c*e*g*h*Pow(k,2)*l*m*n - 2816*Pow(a,2)*b*Pow(e,2)*g*h*Pow(k,2)*l*m*n -
               2048*a*Pow(b,3)*Pow(g,2)*h*Pow(k,2)*l*m*n - 16*c*Pow(d,5)*Pow(e,2)*Pow(l,2)*m*n + 32*Pow(c,2)*Pow(d,4)*e*f*Pow(l,2)*m*n -
               16*Pow(c,3)*Pow(d,3)*Pow(f,2)*Pow(l,2)*m*n - 160*b*Pow(d,5)*e*g*Pow(l,2)*m*n - 128*b*c*Pow(d,4)*f*g*Pow(l,2)*m*n + 208*a*Pow(d,4)*e*f*g*Pow(l,2)*m*n +
               80*a*c*Pow(d,3)*Pow(f,2)*g*Pow(l,2)*m*n - 128*b*c*Pow(d,4)*e*h*Pow(l,2)*m*n + 160*a*Pow(d,4)*Pow(e,2)*h*Pow(l,2)*m*n +
               416*b*Pow(c,2)*Pow(d,3)*f*h*Pow(l,2)*m*n - 288*a*c*Pow(d,3)*e*f*h*Pow(l,2)*m*n - 160*a*Pow(c,2)*Pow(d,2)*Pow(f,2)*h*Pow(l,2)*m*n +
               640*Pow(b,2)*Pow(d,4)*g*h*Pow(l,2)*m*n - 1152*a*b*Pow(d,3)*f*g*h*Pow(l,2)*m*n + 704*Pow(a,2)*Pow(d,2)*Pow(f,2)*g*h*Pow(l,2)*m*n -
               384*Pow(b,2)*c*Pow(d,3)*Pow(h,2)*Pow(l,2)*m*n + 576*a*b*Pow(d,3)*e*Pow(h,2)*Pow(l,2)*m*n - 64*a*b*c*Pow(d,2)*f*Pow(h,2)*Pow(l,2)*m*n +
               192*Pow(a,2)*Pow(d,2)*e*f*Pow(h,2)*Pow(l,2)*m*n - 704*Pow(a,2)*c*d*Pow(f,2)*Pow(h,2)*Pow(l,2)*m*n -
               256*a*Pow(b,2)*Pow(d,2)*Pow(h,3)*Pow(l,2)*m*n - 512*Pow(a,2)*b*d*f*Pow(h,3)*Pow(l,2)*m*n + 2304*Pow(a,3)*Pow(f,2)*Pow(h,3)*Pow(l,2)*m*n +
               144*b*Pow(c,2)*Pow(d,3)*e*k*Pow(l,2)*m*n + 32*a*c*Pow(d,3)*Pow(e,2)*k*Pow(l,2)*m*n - 528*b*Pow(c,3)*Pow(d,2)*f*k*Pow(l,2)*m*n -
               16*a*Pow(c,2)*Pow(d,2)*e*f*k*Pow(l,2)*m*n + 368*a*Pow(c,3)*d*Pow(f,2)*k*Pow(l,2)*m*n - 64*Pow(b,2)*c*Pow(d,3)*g*k*Pow(l,2)*m*n +
               768*a*b*Pow(d,3)*e*g*k*Pow(l,2)*m*n + 1408*a*b*c*Pow(d,2)*f*g*k*Pow(l,2)*m*n - 640*Pow(a,2)*Pow(d,2)*e*f*g*k*Pow(l,2)*m*n -
               1088*Pow(a,2)*c*d*Pow(f,2)*g*k*Pow(l,2)*m*n - 448*Pow(b,2)*Pow(c,2)*Pow(d,2)*h*k*Pow(l,2)*m*n + 896*a*b*c*Pow(d,2)*e*h*k*Pow(l,2)*m*n -
               896*Pow(a,2)*Pow(d,2)*Pow(e,2)*h*k*Pow(l,2)*m*n + 1792*a*b*Pow(c,2)*d*f*h*k*Pow(l,2)*m*n - 896*Pow(a,2)*c*d*e*f*h*k*Pow(l,2)*m*n -
               448*Pow(a,2)*Pow(c,2)*Pow(f,2)*h*k*Pow(l,2)*m*n - 2816*a*Pow(b,2)*Pow(d,2)*g*h*k*Pow(l,2)*m*n + 1536*Pow(a,2)*b*d*f*g*h*k*Pow(l,2)*m*n +
               256*Pow(a,3)*Pow(f,2)*g*h*k*Pow(l,2)*m*n + 2816*a*Pow(b,2)*c*d*Pow(h,2)*k*Pow(l,2)*m*n - 4352*Pow(a,2)*b*d*e*Pow(h,2)*k*Pow(l,2)*m*n -
               4352*Pow(a,2)*b*c*f*Pow(h,2)*k*Pow(l,2)*m*n + 3328*Pow(a,3)*e*f*Pow(h,2)*k*Pow(l,2)*m*n + 2048*Pow(a,2)*Pow(b,2)*Pow(h,3)*k*Pow(l,2)*m*n +
               640*Pow(b,2)*Pow(c,3)*d*Pow(k,2)*Pow(l,2)*m*n - 1088*a*b*Pow(c,2)*d*e*Pow(k,2)*Pow(l,2)*m*n + 128*Pow(a,2)*c*d*Pow(e,2)*Pow(k,2)*Pow(l,2)*m*n -
               64*a*b*Pow(c,3)*f*Pow(k,2)*Pow(l,2)*m*n + 576*Pow(a,2)*Pow(c,2)*e*f*Pow(k,2)*Pow(l,2)*m*n - 512*a*Pow(b,2)*c*d*g*Pow(k,2)*Pow(l,2)*m*n -
               512*Pow(a,2)*b*d*e*g*Pow(k,2)*Pow(l,2)*m*n - 512*Pow(a,2)*b*c*f*g*Pow(k,2)*Pow(l,2)*m*n - 768*Pow(a,3)*e*f*g*Pow(k,2)*Pow(l,2)*m*n -
               2304*a*Pow(b,2)*Pow(c,2)*h*Pow(k,2)*Pow(l,2)*m*n + 2560*Pow(a,2)*b*c*e*h*Pow(k,2)*Pow(l,2)*m*n + 1024*Pow(a,3)*Pow(e,2)*h*Pow(k,2)*Pow(l,2)*m*n +
               4096*Pow(a,2)*Pow(b,2)*g*h*Pow(k,2)*Pow(l,2)*m*n + 96*b*c*Pow(d,5)*Pow(l,3)*m*n - 16*a*Pow(d,5)*e*Pow(l,3)*m*n - 80*a*c*Pow(d,4)*f*Pow(l,3)*m*n -
               320*a*b*Pow(d,4)*h*Pow(l,3)*m*n + 256*Pow(a,2)*Pow(d,3)*f*h*Pow(l,3)*m*n - 640*a*b*c*Pow(d,3)*k*Pow(l,3)*m*n +
               128*Pow(a,2)*Pow(d,3)*e*k*Pow(l,3)*m*n + 384*Pow(a,2)*c*Pow(d,2)*f*k*Pow(l,3)*m*n + 2048*Pow(a,2)*b*Pow(d,2)*h*k*Pow(l,3)*m*n -
               1024*Pow(a,3)*d*f*h*k*Pow(l,3)*m*n + 1024*Pow(a,2)*b*c*d*Pow(k,2)*Pow(l,3)*m*n - 256*Pow(a,3)*d*e*Pow(k,2)*Pow(l,3)*m*n -
               256*Pow(a,3)*c*f*Pow(k,2)*Pow(l,3)*m*n - 3072*Pow(a,3)*b*h*Pow(k,2)*Pow(l,3)*m*n + 24*c*Pow(d,4)*Pow(e,3)*g*Pow(m,2)*n -
               80*Pow(c,2)*Pow(d,3)*Pow(e,2)*f*g*Pow(m,2)*n + 88*Pow(c,3)*Pow(d,2)*e*Pow(f,2)*g*Pow(m,2)*n - 32*Pow(c,4)*d*Pow(f,3)*g*Pow(m,2)*n +
               24*b*Pow(d,4)*Pow(e,2)*Pow(g,2)*Pow(m,2)*n - 256*b*c*Pow(d,3)*e*f*Pow(g,2)*Pow(m,2)*n + 208*a*Pow(d,3)*Pow(e,2)*f*Pow(g,2)*Pow(m,2)*n +
               88*b*Pow(c,2)*Pow(d,2)*Pow(f,2)*Pow(g,2)*Pow(m,2)*n - 160*a*c*Pow(d,2)*e*Pow(f,2)*Pow(g,2)*Pow(m,2)*n +
               96*a*Pow(c,2)*d*Pow(f,3)*Pow(g,2)*Pow(m,2)*n + 576*Pow(b,2)*Pow(d,3)*f*Pow(g,3)*Pow(m,2)*n - 736*a*b*Pow(d,2)*Pow(f,2)*Pow(g,3)*Pow(m,2)*n +
               128*Pow(a,2)*d*Pow(f,3)*Pow(g,3)*Pow(m,2)*n - 64*Pow(c,2)*Pow(d,3)*Pow(e,3)*h*Pow(m,2)*n + 144*Pow(c,3)*Pow(d,2)*Pow(e,2)*f*h*Pow(m,2)*n -
               96*Pow(c,4)*d*e*Pow(f,2)*h*Pow(m,2)*n + 16*Pow(c,5)*Pow(f,3)*h*Pow(m,2)*n - 800*b*c*Pow(d,3)*Pow(e,2)*g*h*Pow(m,2)*n +
               368*a*Pow(d,3)*Pow(e,3)*g*h*Pow(m,2)*n + 400*b*Pow(c,2)*Pow(d,2)*e*f*g*h*Pow(m,2)*n + 400*a*c*Pow(d,2)*Pow(e,2)*f*g*h*Pow(m,2)*n +
               208*b*Pow(c,3)*d*Pow(f,2)*g*h*Pow(m,2)*n - 608*a*Pow(c,2)*d*e*Pow(f,2)*g*h*Pow(m,2)*n + 32*a*Pow(c,3)*Pow(f,3)*g*h*Pow(m,2)*n +
               960*Pow(b,2)*Pow(d,3)*e*Pow(g,2)*h*Pow(m,2)*n - 704*Pow(b,2)*c*Pow(d,2)*f*Pow(g,2)*h*Pow(m,2)*n - 1856*a*b*Pow(d,2)*e*f*Pow(g,2)*h*Pow(m,2)*n -
               64*a*b*c*d*Pow(f,2)*Pow(g,2)*h*Pow(m,2)*n + 2432*Pow(a,2)*d*e*Pow(f,2)*Pow(g,2)*h*Pow(m,2)*n - 384*Pow(a,2)*c*Pow(f,3)*Pow(g,2)*h*Pow(m,2)*n +
               1088*b*Pow(c,2)*Pow(d,2)*Pow(e,2)*Pow(h,2)*Pow(m,2)*n - 320*a*c*Pow(d,2)*Pow(e,3)*Pow(h,2)*Pow(m,2)*n + 32*b*Pow(c,3)*d*e*f*Pow(h,2)*Pow(m,2)*n -
               1184*a*Pow(c,2)*d*Pow(e,2)*f*Pow(h,2)*Pow(m,2)*n + 32*b*Pow(c,4)*Pow(f,2)*Pow(h,2)*Pow(m,2)*n + 352*a*Pow(c,3)*e*Pow(f,2)*Pow(h,2)*Pow(m,2)*n -
               192*Pow(b,2)*c*Pow(d,2)*e*g*Pow(h,2)*Pow(m,2)*n + 32*a*b*Pow(d,2)*Pow(e,2)*g*Pow(h,2)*Pow(m,2)*n -
               640*Pow(b,2)*Pow(c,2)*d*f*g*Pow(h,2)*Pow(m,2)*n + 1920*a*b*c*d*e*f*g*Pow(h,2)*Pow(m,2)*n - 192*Pow(a,2)*d*Pow(e,2)*f*g*Pow(h,2)*Pow(m,2)*n -
               864*a*b*Pow(c,2)*Pow(f,2)*g*Pow(h,2)*Pow(m,2)*n - 640*Pow(a,2)*c*e*Pow(f,2)*g*Pow(h,2)*Pow(m,2)*n -
               1920*Pow(b,3)*Pow(d,2)*Pow(g,2)*Pow(h,2)*Pow(m,2)*n + 2816*a*Pow(b,2)*d*f*Pow(g,2)*Pow(h,2)*Pow(m,2)*n +
               1408*Pow(a,2)*b*Pow(f,2)*Pow(g,2)*Pow(h,2)*Pow(m,2)*n - 1664*Pow(b,2)*Pow(c,2)*d*e*Pow(h,3)*Pow(m,2)*n -
               1088*a*b*c*d*Pow(e,2)*Pow(h,3)*Pow(m,2)*n + 192*Pow(a,2)*d*Pow(e,3)*Pow(h,3)*Pow(m,2)*n - 256*Pow(b,2)*Pow(c,3)*f*Pow(h,3)*Pow(m,2)*n +
               704*a*b*Pow(c,2)*e*f*Pow(h,3)*Pow(m,2)*n + 1472*Pow(a,2)*c*Pow(e,2)*f*Pow(h,3)*Pow(m,2)*n + 2816*Pow(b,3)*c*d*g*Pow(h,3)*Pow(m,2)*n +
               1280*a*Pow(b,2)*d*e*g*Pow(h,3)*Pow(m,2)*n + 1280*a*Pow(b,2)*c*f*g*Pow(h,3)*Pow(m,2)*n - 5888*Pow(a,2)*b*e*f*g*Pow(h,3)*Pow(m,2)*n +
               256*Pow(b,3)*Pow(c,2)*Pow(h,4)*Pow(m,2)*n + 1536*a*Pow(b,2)*c*e*Pow(h,4)*Pow(m,2)*n + 768*Pow(a,2)*b*Pow(e,2)*Pow(h,4)*Pow(m,2)*n -
               4096*a*Pow(b,3)*g*Pow(h,4)*Pow(m,2)*n + 8*Pow(c,3)*Pow(d,2)*Pow(e,3)*k*Pow(m,2)*n - 16*Pow(c,4)*d*Pow(e,2)*f*k*Pow(m,2)*n +
               8*Pow(c,5)*e*Pow(f,2)*k*Pow(m,2)*n + 8*b*Pow(c,2)*Pow(d,2)*Pow(e,2)*g*k*Pow(m,2)*n + 80*a*c*Pow(d,2)*Pow(e,3)*g*k*Pow(m,2)*n -
               64*b*Pow(c,3)*d*e*f*g*k*Pow(m,2)*n - 208*a*Pow(c,2)*d*Pow(e,2)*f*g*k*Pow(m,2)*n - 232*b*Pow(c,4)*Pow(f,2)*g*k*Pow(m,2)*n +
               416*a*Pow(c,3)*e*Pow(f,2)*g*k*Pow(m,2)*n + 896*Pow(b,2)*c*Pow(d,2)*e*Pow(g,2)*k*Pow(m,2)*n - 896*a*b*Pow(d,2)*Pow(e,2)*Pow(g,2)*k*Pow(m,2)*n -
               448*Pow(b,2)*Pow(c,2)*d*f*Pow(g,2)*k*Pow(m,2)*n + 1792*a*b*c*d*e*f*Pow(g,2)*k*Pow(m,2)*n - 448*Pow(a,2)*d*Pow(e,2)*f*Pow(g,2)*k*Pow(m,2)*n +
               896*a*b*Pow(c,2)*Pow(f,2)*Pow(g,2)*k*Pow(m,2)*n - 1792*Pow(a,2)*c*e*Pow(f,2)*Pow(g,2)*k*Pow(m,2)*n -
               1152*Pow(b,3)*Pow(d,2)*Pow(g,3)*k*Pow(m,2)*n + 256*a*Pow(b,2)*d*f*Pow(g,3)*k*Pow(m,2)*n + 128*Pow(a,2)*b*Pow(f,2)*Pow(g,3)*k*Pow(m,2)*n -
               272*b*Pow(c,3)*d*Pow(e,2)*h*k*Pow(m,2)*n + 176*a*Pow(c,2)*d*Pow(e,3)*h*k*Pow(m,2)*n - 304*b*Pow(c,4)*e*f*h*k*Pow(m,2)*n +
               400*a*Pow(c,3)*Pow(e,2)*f*h*k*Pow(m,2)*n - 1856*Pow(b,2)*Pow(c,2)*d*e*g*h*k*Pow(m,2)*n + 4352*a*b*c*d*Pow(e,2)*g*h*k*Pow(m,2)*n -
               2112*Pow(a,2)*d*Pow(e,3)*g*h*k*Pow(m,2)*n + 1472*Pow(b,2)*Pow(c,3)*f*g*h*k*Pow(m,2)*n - 1024*a*b*Pow(c,2)*e*f*g*h*k*Pow(m,2)*n -
               64*Pow(a,2)*c*Pow(e,2)*f*g*h*k*Pow(m,2)*n + 3328*Pow(b,3)*c*d*Pow(g,2)*h*k*Pow(m,2)*n - 4352*a*Pow(b,2)*d*e*Pow(g,2)*h*k*Pow(m,2)*n -
               4352*a*Pow(b,2)*c*f*Pow(g,2)*h*k*Pow(m,2)*n + 2816*Pow(a,2)*b*e*f*Pow(g,2)*h*k*Pow(m,2)*n + 1472*Pow(b,2)*Pow(c,3)*e*Pow(h,2)*k*Pow(m,2)*n -
               2752*a*b*Pow(c,2)*Pow(e,2)*Pow(h,2)*k*Pow(m,2)*n + 1472*Pow(a,2)*c*Pow(e,3)*Pow(h,2)*k*Pow(m,2)*n -
               2944*Pow(b,3)*Pow(c,2)*g*Pow(h,2)*k*Pow(m,2)*n + 3840*a*Pow(b,2)*c*e*g*Pow(h,2)*k*Pow(m,2)*n - 1664*Pow(a,2)*b*Pow(e,2)*g*Pow(h,2)*k*Pow(m,2)*n +
               2048*a*Pow(b,3)*Pow(g,2)*Pow(h,2)*k*Pow(m,2)*n + 48*b*Pow(c,4)*Pow(e,2)*Pow(k,2)*Pow(m,2)*n - 48*a*Pow(c,3)*Pow(e,3)*Pow(k,2)*Pow(m,2)*n +
               384*Pow(b,2)*Pow(c,3)*e*g*Pow(k,2)*Pow(m,2)*n - 416*a*b*Pow(c,2)*Pow(e,2)*g*Pow(k,2)*Pow(m,2)*n - 64*Pow(a,2)*c*Pow(e,3)*g*Pow(k,2)*Pow(m,2)*n -
               384*Pow(b,3)*Pow(c,2)*Pow(g,2)*Pow(k,2)*Pow(m,2)*n - 2304*a*Pow(b,2)*c*e*Pow(g,2)*Pow(k,2)*Pow(m,2)*n +
               2432*Pow(a,2)*b*Pow(e,2)*Pow(g,2)*Pow(k,2)*Pow(m,2)*n + 2048*a*Pow(b,3)*Pow(g,3)*Pow(k,2)*Pow(m,2)*n +
               8*Pow(c,2)*Pow(d,4)*Pow(e,2)*l*Pow(m,2)*n - 16*Pow(c,3)*Pow(d,3)*e*f*l*Pow(m,2)*n + 8*Pow(c,4)*Pow(d,2)*Pow(f,2)*l*Pow(m,2)*n +
               208*b*c*Pow(d,4)*e*g*l*Pow(m,2)*n - 232*a*Pow(d,4)*Pow(e,2)*g*l*Pow(m,2)*n + 80*b*Pow(c,2)*Pow(d,3)*f*g*l*Pow(m,2)*n +
               160*a*c*Pow(d,3)*e*f*g*l*Pow(m,2)*n - 216*a*Pow(c,2)*Pow(d,2)*Pow(f,2)*g*l*Pow(m,2)*n - 576*Pow(b,2)*Pow(d,4)*Pow(g,2)*l*Pow(m,2)*n +
               320*a*b*Pow(d,3)*f*Pow(g,2)*l*Pow(m,2)*n + 352*Pow(a,2)*Pow(d,2)*Pow(f,2)*Pow(g,2)*l*Pow(m,2)*n + 96*b*Pow(c,2)*Pow(d,3)*e*h*l*Pow(m,2)*n +
               208*a*c*Pow(d,3)*Pow(e,2)*h*l*Pow(m,2)*n + 96*b*Pow(c,3)*Pow(d,2)*f*h*l*Pow(m,2)*n - 608*a*Pow(c,2)*Pow(d,2)*e*f*h*l*Pow(m,2)*n +
               208*a*Pow(c,3)*d*Pow(f,2)*h*l*Pow(m,2)*n + 1472*Pow(b,2)*c*Pow(d,3)*g*h*l*Pow(m,2)*n - 640*a*b*Pow(d,3)*e*g*h*l*Pow(m,2)*n -
               1024*a*b*c*Pow(d,2)*f*g*h*l*Pow(m,2)*n - 512*Pow(a,2)*Pow(d,2)*e*f*g*h*l*Pow(m,2)*n - 64*Pow(a,2)*c*d*Pow(f,2)*g*h*l*Pow(m,2)*n -
               1664*Pow(b,2)*Pow(c,2)*Pow(d,2)*Pow(h,2)*l*Pow(m,2)*n - 832*a*b*c*Pow(d,2)*e*Pow(h,2)*l*Pow(m,2)*n +
               352*Pow(a,2)*Pow(d,2)*Pow(e,2)*Pow(h,2)*l*Pow(m,2)*n - 832*a*b*Pow(c,2)*d*f*Pow(h,2)*l*Pow(m,2)*n + 3200*Pow(a,2)*c*d*e*f*Pow(h,2)*l*Pow(m,2)*n +
               352*Pow(a,2)*Pow(c,2)*Pow(f,2)*Pow(h,2)*l*Pow(m,2)*n - 1664*a*Pow(b,2)*Pow(d,2)*g*Pow(h,2)*l*Pow(m,2)*n +
               3840*Pow(a,2)*b*d*f*g*Pow(h,2)*l*Pow(m,2)*n - 2944*Pow(a,3)*Pow(f,2)*g*Pow(h,2)*l*Pow(m,2)*n + 6656*a*Pow(b,2)*c*d*Pow(h,3)*l*Pow(m,2)*n -
               512*Pow(a,2)*b*d*e*Pow(h,3)*l*Pow(m,2)*n - 512*Pow(a,2)*b*c*f*Pow(h,3)*l*Pow(m,2)*n - 2560*Pow(a,3)*e*f*Pow(h,3)*l*Pow(m,2)*n -
               6144*Pow(a,2)*Pow(b,2)*Pow(h,4)*l*Pow(m,2)*n + 80*b*Pow(c,3)*Pow(d,2)*e*k*l*Pow(m,2)*n - 216*a*Pow(c,2)*Pow(d,2)*Pow(e,2)*k*l*Pow(m,2)*n +
               208*b*Pow(c,4)*d*f*k*l*Pow(m,2)*n + 160*a*Pow(c,3)*d*e*f*k*l*Pow(m,2)*n - 232*a*Pow(c,4)*Pow(f,2)*k*l*Pow(m,2)*n -
               896*Pow(b,2)*Pow(c,2)*Pow(d,2)*g*k*l*Pow(m,2)*n - 896*a*b*c*Pow(d,2)*e*g*k*l*Pow(m,2)*n + 896*Pow(a,2)*Pow(d,2)*Pow(e,2)*g*k*l*Pow(m,2)*n -
               896*a*b*Pow(c,2)*d*f*g*k*l*Pow(m,2)*n + 896*Pow(a,2)*c*d*e*f*g*k*l*Pow(m,2)*n + 896*Pow(a,2)*Pow(c,2)*Pow(f,2)*g*k*l*Pow(m,2)*n +
               3968*a*Pow(b,2)*Pow(d,2)*Pow(g,2)*k*l*Pow(m,2)*n - 2816*Pow(a,2)*b*d*f*Pow(g,2)*k*l*Pow(m,2)*n + 128*Pow(a,3)*Pow(f,2)*Pow(g,2)*k*l*Pow(m,2)*n +
               1472*Pow(b,2)*Pow(c,3)*d*h*k*l*Pow(m,2)*n - 1024*a*b*Pow(c,2)*d*e*h*k*l*Pow(m,2)*n - 64*Pow(a,2)*c*d*Pow(e,2)*h*k*l*Pow(m,2)*n -
               640*a*b*Pow(c,3)*f*h*k*l*Pow(m,2)*n - 512*Pow(a,2)*Pow(c,2)*e*f*h*k*l*Pow(m,2)*n - 5120*a*Pow(b,2)*c*d*g*h*k*l*Pow(m,2)*n +
               5632*Pow(a,2)*b*d*e*g*h*k*l*Pow(m,2)*n + 5632*Pow(a,2)*b*c*f*g*h*k*l*Pow(m,2)*n - 4096*Pow(a,3)*e*f*g*h*k*l*Pow(m,2)*n -
               1664*a*Pow(b,2)*Pow(c,2)*Pow(h,2)*k*l*Pow(m,2)*n + 3840*Pow(a,2)*b*c*e*Pow(h,2)*k*l*Pow(m,2)*n - 2944*Pow(a,3)*Pow(e,2)*Pow(h,2)*k*l*Pow(m,2)*n -
               1024*Pow(a,2)*Pow(b,2)*g*Pow(h,2)*k*l*Pow(m,2)*n - 576*Pow(b,2)*Pow(c,4)*Pow(k,2)*l*Pow(m,2)*n + 320*a*b*Pow(c,3)*e*Pow(k,2)*l*Pow(m,2)*n +
               352*Pow(a,2)*Pow(c,2)*Pow(e,2)*Pow(k,2)*l*Pow(m,2)*n + 3968*a*Pow(b,2)*Pow(c,2)*g*Pow(k,2)*l*Pow(m,2)*n -
               2816*Pow(a,2)*b*c*e*g*Pow(k,2)*l*Pow(m,2)*n + 128*Pow(a,3)*Pow(e,2)*g*Pow(k,2)*l*Pow(m,2)*n -
               5120*Pow(a,2)*Pow(b,2)*Pow(g,2)*Pow(k,2)*l*Pow(m,2)*n - 144*b*Pow(c,2)*Pow(d,4)*Pow(l,2)*Pow(m,2)*n + 48*a*c*Pow(d,4)*e*Pow(l,2)*Pow(m,2)*n +
               96*a*Pow(c,2)*Pow(d,3)*f*Pow(l,2)*Pow(m,2)*n + 416*a*b*Pow(d,4)*g*Pow(l,2)*Pow(m,2)*n - 512*Pow(a,2)*Pow(d,3)*f*g*Pow(l,2)*Pow(m,2)*n +
               128*a*b*c*Pow(d,3)*h*Pow(l,2)*Pow(m,2)*n - 384*Pow(a,2)*Pow(d,3)*e*h*Pow(l,2)*Pow(m,2)*n + 640*Pow(a,2)*c*Pow(d,2)*f*h*Pow(l,2)*Pow(m,2)*n +
               512*Pow(a,2)*b*Pow(d,2)*Pow(h,2)*Pow(l,2)*Pow(m,2)*n - 2048*Pow(a,3)*d*f*Pow(h,2)*Pow(l,2)*Pow(m,2)*n +
               896*a*b*Pow(c,2)*Pow(d,2)*k*Pow(l,2)*Pow(m,2)*n - 896*Pow(a,2)*Pow(c,2)*d*f*k*Pow(l,2)*Pow(m,2)*n -
               2304*Pow(a,2)*b*Pow(d,2)*g*k*Pow(l,2)*Pow(m,2)*n + 2048*Pow(a,3)*d*f*g*k*Pow(l,2)*Pow(m,2)*n - 2560*Pow(a,2)*b*c*d*h*k*Pow(l,2)*Pow(m,2)*n +
               1536*Pow(a,3)*d*e*h*k*Pow(l,2)*Pow(m,2)*n + 1536*Pow(a,3)*c*f*h*k*Pow(l,2)*Pow(m,2)*n + 2048*Pow(a,3)*b*Pow(h,2)*k*Pow(l,2)*Pow(m,2)*n -
               256*Pow(a,2)*b*Pow(c,2)*Pow(k,2)*Pow(l,2)*Pow(m,2)*n - 768*Pow(a,3)*c*e*Pow(k,2)*Pow(l,2)*Pow(m,2)*n +
               2560*Pow(a,3)*b*g*Pow(k,2)*Pow(l,2)*Pow(m,2)*n + 32*Pow(a,2)*Pow(d,4)*Pow(l,3)*Pow(m,2)*n - 256*Pow(a,3)*Pow(d,2)*k*Pow(l,3)*Pow(m,2)*n +
               512*Pow(a,4)*Pow(k,2)*Pow(l,3)*Pow(m,2)*n + 16*Pow(c,3)*Pow(d,3)*Pow(e,2)*Pow(m,3)*n - 32*Pow(c,4)*Pow(d,2)*e*f*Pow(m,3)*n +
               16*Pow(c,5)*d*Pow(f,2)*Pow(m,3)*n + 368*b*Pow(c,2)*Pow(d,3)*e*g*Pow(m,3)*n - 304*a*c*Pow(d,3)*Pow(e,2)*g*Pow(m,3)*n -
               304*b*Pow(c,3)*Pow(d,2)*f*g*Pow(m,3)*n + 208*a*Pow(c,2)*Pow(d,2)*e*f*g*Pow(m,3)*n + 32*a*Pow(c,3)*d*Pow(f,2)*g*Pow(m,3)*n -
               576*Pow(b,2)*c*Pow(d,3)*Pow(g,2)*Pow(m,3)*n + 192*a*b*Pow(d,3)*e*Pow(g,2)*Pow(m,3)*n + 1472*a*b*c*Pow(d,2)*f*Pow(g,2)*Pow(m,3)*n -
               832*Pow(a,2)*Pow(d,2)*e*f*Pow(g,2)*Pow(m,3)*n - 384*Pow(a,2)*c*d*Pow(f,2)*Pow(g,2)*Pow(m,3)*n - 736*b*Pow(c,3)*Pow(d,2)*e*h*Pow(m,3)*n +
               352*a*Pow(c,2)*Pow(d,2)*Pow(e,2)*h*Pow(m,3)*n - 32*b*Pow(c,4)*d*f*h*Pow(m,3)*n + 544*a*Pow(c,3)*d*e*f*h*Pow(m,3)*n -
               128*a*Pow(c,4)*Pow(f,2)*h*Pow(m,3)*n + 832*Pow(b,2)*Pow(c,2)*Pow(d,2)*g*h*Pow(m,3)*n - 256*a*b*c*Pow(d,2)*e*g*h*Pow(m,3)*n -
               64*Pow(a,2)*Pow(d,2)*Pow(e,2)*g*h*Pow(m,3)*n + 640*a*b*Pow(c,2)*d*f*g*h*Pow(m,3)*n - 1152*Pow(a,2)*c*d*e*f*g*h*Pow(m,3)*n +
               384*Pow(a,2)*Pow(c,2)*Pow(f,2)*g*h*Pow(m,3)*n + 1536*a*Pow(b,2)*Pow(d,2)*Pow(g,2)*h*Pow(m,3)*n - 4096*Pow(a,2)*b*d*f*Pow(g,2)*h*Pow(m,3)*n +
               512*Pow(a,3)*Pow(f,2)*Pow(g,2)*h*Pow(m,3)*n + 640*Pow(b,2)*Pow(c,3)*d*Pow(h,2)*Pow(m,3)*n + 1600*a*b*Pow(c,2)*d*e*Pow(h,2)*Pow(m,3)*n -
               320*Pow(a,2)*c*d*Pow(e,2)*Pow(h,2)*Pow(m,3)*n - 64*a*b*Pow(c,3)*f*Pow(h,2)*Pow(m,3)*n - 1216*Pow(a,2)*Pow(c,2)*e*f*Pow(h,2)*Pow(m,3)*n -
               5888*a*Pow(b,2)*c*d*g*Pow(h,2)*Pow(m,3)*n + 1280*Pow(a,2)*b*d*e*g*Pow(h,2)*Pow(m,3)*n + 1280*Pow(a,2)*b*c*f*g*Pow(h,2)*Pow(m,3)*n +
               2816*Pow(a,3)*e*f*g*Pow(h,2)*Pow(m,3)*n - 768*a*Pow(b,2)*Pow(c,2)*Pow(h,3)*Pow(m,3)*n - 1536*Pow(a,2)*b*c*e*Pow(h,3)*Pow(m,3)*n -
               256*Pow(a,3)*Pow(e,2)*Pow(h,3)*Pow(m,3)*n + 6144*Pow(a,2)*Pow(b,2)*g*Pow(h,3)*Pow(m,3)*n + 48*b*Pow(c,4)*d*e*k*Pow(m,3)*n -
               16*a*Pow(c,3)*d*Pow(e,2)*k*Pow(m,3)*n + 144*b*Pow(c,5)*f*k*Pow(m,3)*n - 176*a*Pow(c,4)*e*f*k*Pow(m,3)*n + 192*Pow(b,2)*Pow(c,3)*d*g*k*Pow(m,3)*n -
               640*a*b*Pow(c,2)*d*e*g*k*Pow(m,3)*n + 576*Pow(a,2)*c*d*Pow(e,2)*g*k*Pow(m,3)*n - 512*a*b*Pow(c,3)*f*g*k*Pow(m,3)*n +
               128*Pow(a,2)*Pow(c,2)*e*f*g*k*Pow(m,3)*n - 256*a*Pow(b,2)*c*d*Pow(g,2)*k*Pow(m,3)*n - 256*Pow(a,2)*b*d*e*Pow(g,2)*k*Pow(m,3)*n -
               256*Pow(a,2)*b*c*f*Pow(g,2)*k*Pow(m,3)*n + 2304*Pow(a,3)*e*f*Pow(g,2)*k*Pow(m,3)*n - 576*Pow(b,2)*Pow(c,4)*h*k*Pow(m,3)*n +
               1664*a*b*Pow(c,3)*e*h*k*Pow(m,3)*n - 1216*Pow(a,2)*Pow(c,2)*Pow(e,2)*h*k*Pow(m,3)*n + 1280*a*Pow(b,2)*Pow(c,2)*g*h*k*Pow(m,3)*n -
               4608*Pow(a,2)*b*c*e*g*h*k*Pow(m,3)*n + 2816*Pow(a,3)*Pow(e,2)*g*h*k*Pow(m,3)*n + 2048*Pow(a,2)*Pow(b,2)*Pow(g,2)*h*k*Pow(m,3)*n -
               64*b*Pow(c,3)*Pow(d,3)*l*Pow(m,3)*n + 32*a*Pow(c,2)*Pow(d,3)*e*l*Pow(m,3)*n + 32*a*Pow(c,3)*Pow(d,2)*f*l*Pow(m,3)*n -
               512*a*b*c*Pow(d,3)*g*l*Pow(m,3)*n + 640*Pow(a,2)*Pow(d,3)*e*g*l*Pow(m,3)*n + 128*Pow(a,2)*c*Pow(d,2)*f*g*l*Pow(m,3)*n +
               1664*a*b*Pow(c,2)*Pow(d,2)*h*l*Pow(m,3)*n - 1024*Pow(a,2)*c*Pow(d,2)*e*h*l*Pow(m,3)*n - 1024*Pow(a,2)*Pow(c,2)*d*f*h*l*Pow(m,3)*n -
               512*Pow(a,2)*b*Pow(d,2)*g*h*l*Pow(m,3)*n + 2048*Pow(a,3)*d*f*g*h*l*Pow(m,3)*n - 4096*Pow(a,2)*b*c*d*Pow(h,2)*l*Pow(m,3)*n +
               1024*Pow(a,3)*d*e*Pow(h,2)*l*Pow(m,3)*n + 1024*Pow(a,3)*c*f*Pow(h,2)*l*Pow(m,3)*n + 4096*Pow(a,3)*b*Pow(h,3)*l*Pow(m,3)*n -
               512*a*b*Pow(c,3)*d*k*l*Pow(m,3)*n + 128*Pow(a,2)*Pow(c,2)*d*e*k*l*Pow(m,3)*n + 640*Pow(a,2)*Pow(c,3)*f*k*l*Pow(m,3)*n +
               3072*Pow(a,2)*b*c*d*g*k*l*Pow(m,3)*n - 2560*Pow(a,3)*d*e*g*k*l*Pow(m,3)*n - 2560*Pow(a,3)*c*f*g*k*l*Pow(m,3)*n -
               512*Pow(a,2)*b*Pow(c,2)*h*k*l*Pow(m,3)*n + 2048*Pow(a,3)*c*e*h*k*l*Pow(m,3)*n - 2048*Pow(a,3)*b*g*h*k*l*Pow(m,3)*n -
               128*Pow(a,2)*c*Pow(d,3)*Pow(l,2)*Pow(m,3)*n + 512*Pow(a,3)*Pow(d,2)*h*Pow(l,2)*Pow(m,3)*n + 512*Pow(a,3)*c*d*k*Pow(l,2)*Pow(m,3)*n -
               2048*Pow(a,4)*h*k*Pow(l,2)*Pow(m,3)*n + 192*b*Pow(c,4)*Pow(d,2)*Pow(m,4)*n - 128*a*Pow(c,3)*Pow(d,2)*e*Pow(m,4)*n - 64*a*Pow(c,4)*d*f*Pow(m,4)*n -
               928*a*b*Pow(c,2)*Pow(d,2)*g*Pow(m,4)*n + 640*Pow(a,2)*c*Pow(d,2)*e*g*Pow(m,4)*n + 192*Pow(a,2)*Pow(c,2)*d*f*g*Pow(m,4)*n +
               384*Pow(a,2)*b*Pow(d,2)*Pow(g,2)*Pow(m,4)*n + 256*Pow(a,3)*d*f*Pow(g,2)*Pow(m,4)*n - 704*a*b*Pow(c,3)*d*h*Pow(m,4)*n +
               64*Pow(a,2)*Pow(c,2)*d*e*h*Pow(m,4)*n + 320*Pow(a,2)*Pow(c,3)*f*h*Pow(m,4)*n + 3328*Pow(a,2)*b*c*d*g*h*Pow(m,4)*n - 1280*Pow(a,3)*d*e*g*h*Pow(m,4)*n -
               1280*Pow(a,3)*c*f*g*h*Pow(m,4)*n + 768*Pow(a,2)*b*Pow(c,2)*Pow(h,2)*Pow(m,4)*n + 512*Pow(a,3)*c*e*Pow(h,2)*Pow(m,4)*n -
               4096*Pow(a,3)*b*g*Pow(h,2)*Pow(m,4)*n - 288*a*b*Pow(c,4)*k*Pow(m,4)*n + 320*Pow(a,2)*Pow(c,3)*e*k*Pow(m,4)*n +
               1664*Pow(a,2)*b*Pow(c,2)*g*k*Pow(m,4)*n - 1280*Pow(a,3)*c*e*g*k*Pow(m,4)*n - 2048*Pow(a,3)*b*Pow(g,2)*k*Pow(m,4)*n +
               96*Pow(a,2)*Pow(c,2)*Pow(d,2)*l*Pow(m,4)*n - 640*Pow(a,3)*Pow(d,2)*g*l*Pow(m,4)*n + 512*Pow(a,3)*c*d*h*l*Pow(m,4)*n -
               1024*Pow(a,4)*Pow(h,2)*l*Pow(m,4)*n - 640*Pow(a,3)*Pow(c,2)*k*l*Pow(m,4)*n + 2560*Pow(a,4)*g*k*l*Pow(m,4)*n + 64*Pow(a,2)*Pow(c,3)*d*Pow(m,5)*n -
               256*Pow(a,3)*c*d*g*Pow(m,5)*n - 256*Pow(a,3)*Pow(c,2)*h*Pow(m,5)*n + 1024*Pow(a,4)*g*h*Pow(m,5)*n + 27*Pow(d,4)*Pow(e,4)*Pow(g,2)*Pow(n,2) -
               72*c*Pow(d,3)*Pow(e,3)*f*Pow(g,2)*Pow(n,2) + 62*Pow(c,2)*Pow(d,2)*Pow(e,2)*Pow(f,2)*Pow(g,2)*Pow(n,2) -
               16*Pow(c,3)*d*e*Pow(f,3)*Pow(g,2)*Pow(n,2) - Pow(c,4)*Pow(f,4)*Pow(g,2)*Pow(n,2) + 144*b*Pow(d,3)*Pow(e,2)*f*Pow(g,3)*Pow(n,2) -
               128*b*c*Pow(d,2)*e*Pow(f,2)*Pow(g,3)*Pow(n,2) - 120*a*Pow(d,2)*Pow(e,2)*Pow(f,2)*Pow(g,3)*Pow(n,2) +
               32*b*Pow(c,2)*d*Pow(f,3)*Pow(g,3)*Pow(n,2) + 64*a*c*d*e*Pow(f,3)*Pow(g,3)*Pow(n,2) + 8*a*Pow(c,2)*Pow(f,4)*Pow(g,3)*Pow(n,2) +
               128*Pow(b,2)*Pow(d,2)*Pow(f,2)*Pow(g,4)*Pow(n,2) - 128*a*b*d*Pow(f,3)*Pow(g,4)*Pow(n,2) - 16*Pow(a,2)*Pow(f,4)*Pow(g,4)*Pow(n,2) -
               36*c*Pow(d,3)*Pow(e,4)*g*h*Pow(n,2) + 92*Pow(c,2)*Pow(d,2)*Pow(e,3)*f*g*h*Pow(n,2) - 76*Pow(c,3)*d*Pow(e,2)*Pow(f,2)*g*h*Pow(n,2) +
               20*Pow(c,4)*e*Pow(f,3)*g*h*Pow(n,2) - 144*b*Pow(d,3)*Pow(e,3)*Pow(g,2)*h*Pow(n,2) - 176*b*c*Pow(d,2)*Pow(e,2)*f*Pow(g,2)*h*Pow(n,2) +
               240*a*Pow(d,2)*Pow(e,3)*f*Pow(g,2)*h*Pow(n,2) + 160*b*Pow(c,2)*d*e*Pow(f,2)*Pow(g,2)*h*Pow(n,2) + 48*a*c*d*Pow(e,2)*Pow(f,2)*Pow(g,2)*h*Pow(n,2) -
               32*b*Pow(c,3)*Pow(f,3)*Pow(g,2)*h*Pow(n,2) - 96*a*Pow(c,2)*e*Pow(f,3)*Pow(g,2)*h*Pow(n,2) - 256*Pow(b,2)*Pow(d,2)*e*f*Pow(g,3)*h*Pow(n,2) -
               256*Pow(b,2)*c*d*Pow(f,2)*Pow(g,3)*h*Pow(n,2) + 384*a*b*d*e*Pow(f,2)*Pow(g,3)*h*Pow(n,2) + 128*a*b*c*Pow(f,3)*Pow(g,3)*h*Pow(n,2) +
               64*Pow(a,2)*e*Pow(f,3)*Pow(g,3)*h*Pow(n,2) - 4*Pow(c,2)*Pow(d,2)*Pow(e,4)*Pow(h,2)*Pow(n,2) + 8*Pow(c,3)*d*Pow(e,3)*f*Pow(h,2)*Pow(n,2) -
               4*Pow(c,4)*Pow(e,2)*Pow(f,2)*Pow(h,2)*Pow(n,2) + 208*b*c*Pow(d,2)*Pow(e,3)*g*Pow(h,2)*Pow(n,2) + 24*a*Pow(d,2)*Pow(e,4)*g*Pow(h,2)*Pow(n,2) -
               160*b*Pow(c,2)*d*Pow(e,2)*f*g*Pow(h,2)*Pow(n,2) - 256*a*c*d*Pow(e,3)*f*g*Pow(h,2)*Pow(n,2) + 96*b*Pow(c,3)*e*Pow(f,2)*g*Pow(h,2)*Pow(n,2) +
               88*a*Pow(c,2)*Pow(e,2)*Pow(f,2)*g*Pow(h,2)*Pow(n,2) + 224*Pow(b,2)*Pow(d,2)*Pow(e,2)*Pow(g,2)*Pow(h,2)*Pow(n,2) +
               896*Pow(b,2)*c*d*e*f*Pow(g,2)*Pow(h,2)*Pow(n,2) - 448*a*b*d*Pow(e,2)*f*Pow(g,2)*Pow(h,2)*Pow(n,2) -
               896*a*b*c*e*Pow(f,2)*Pow(g,2)*Pow(h,2)*Pow(n,2) + 224*Pow(a,2)*Pow(e,2)*Pow(f,2)*Pow(g,2)*Pow(h,2)*Pow(n,2) -
               256*Pow(b,3)*d*f*Pow(g,3)*Pow(h,2)*Pow(n,2) + 512*a*Pow(b,2)*Pow(f,2)*Pow(g,3)*Pow(h,2)*Pow(n,2) + 160*b*Pow(c,2)*d*Pow(e,3)*Pow(h,3)*Pow(n,2) -
               144*a*c*d*Pow(e,4)*Pow(h,3)*Pow(n,2) - 64*b*Pow(c,3)*Pow(e,2)*f*Pow(h,3)*Pow(n,2) + 48*a*Pow(c,2)*Pow(e,3)*f*Pow(h,3)*Pow(n,2) -
               832*Pow(b,2)*c*d*Pow(e,2)*g*Pow(h,3)*Pow(n,2) + 192*a*b*d*Pow(e,3)*g*Pow(h,3)*Pow(n,2) - 384*Pow(b,2)*Pow(c,2)*e*f*g*Pow(h,3)*Pow(n,2) +
               1472*a*b*c*Pow(e,2)*f*g*Pow(h,3)*Pow(n,2) - 576*Pow(a,2)*Pow(e,3)*f*g*Pow(h,3)*Pow(n,2) + 256*Pow(b,3)*d*e*Pow(g,2)*Pow(h,3)*Pow(n,2) +
               256*Pow(b,3)*c*f*Pow(g,2)*Pow(h,3)*Pow(n,2) - 1024*a*Pow(b,2)*e*f*Pow(g,2)*Pow(h,3)*Pow(n,2) +
               128*Pow(b,2)*Pow(c,2)*Pow(e,2)*Pow(h,4)*Pow(n,2) - 576*a*b*c*Pow(e,3)*Pow(h,4)*Pow(n,2) + 432*Pow(a,2)*Pow(e,4)*Pow(h,4)*Pow(n,2) +
               256*Pow(b,3)*c*e*g*Pow(h,4)*Pow(n,2) + 384*a*Pow(b,2)*Pow(e,2)*g*Pow(h,4)*Pow(n,2) - 256*Pow(b,4)*Pow(g,2)*Pow(h,4)*Pow(n,2) +
               12*Pow(c,2)*Pow(d,2)*Pow(e,4)*g*k*Pow(n,2) - 24*Pow(c,3)*d*Pow(e,3)*f*g*k*Pow(n,2) + 12*Pow(c,4)*Pow(e,2)*Pow(f,2)*g*k*Pow(n,2) +
               96*b*c*Pow(d,2)*Pow(e,3)*Pow(g,2)*k*Pow(n,2) - 144*a*Pow(d,2)*Pow(e,4)*Pow(g,2)*k*Pow(n,2) + 176*b*Pow(c,2)*d*Pow(e,2)*f*Pow(g,2)*k*Pow(n,2) -
               32*a*c*d*Pow(e,3)*f*Pow(g,2)*k*Pow(n,2) - 128*b*Pow(c,3)*e*Pow(f,2)*Pow(g,2)*k*Pow(n,2) + 32*a*Pow(c,2)*Pow(e,2)*Pow(f,2)*Pow(g,2)*k*Pow(n,2) -
               96*Pow(b,2)*Pow(d,2)*Pow(e,2)*Pow(g,3)*k*Pow(n,2) - 384*Pow(b,2)*c*d*e*f*Pow(g,3)*k*Pow(n,2) + 64*a*b*d*Pow(e,2)*f*Pow(g,3)*k*Pow(n,2) +
               128*Pow(b,2)*Pow(c,2)*Pow(f,2)*Pow(g,3)*k*Pow(n,2) + 512*a*b*c*e*Pow(f,2)*Pow(g,3)*k*Pow(n,2) -
               320*Pow(a,2)*Pow(e,2)*Pow(f,2)*Pow(g,3)*k*Pow(n,2) + 256*Pow(b,3)*d*f*Pow(g,4)*k*Pow(n,2) - 512*a*Pow(b,2)*Pow(f,2)*Pow(g,4)*k*Pow(n,2) -
               368*b*Pow(c,2)*d*Pow(e,3)*g*h*k*Pow(n,2) + 320*a*c*d*Pow(e,4)*g*h*k*Pow(n,2) + 80*b*Pow(c,3)*Pow(e,2)*f*g*h*k*Pow(n,2) -
               32*a*Pow(c,2)*Pow(e,3)*f*g*h*k*Pow(n,2) + 576*Pow(b,2)*c*d*Pow(e,2)*Pow(g,2)*h*k*Pow(n,2) - 64*a*b*d*Pow(e,3)*Pow(g,2)*h*k*Pow(n,2) +
               128*Pow(b,2)*Pow(c,2)*e*f*Pow(g,2)*h*k*Pow(n,2) - 1088*a*b*c*Pow(e,2)*f*Pow(g,2)*h*k*Pow(n,2) + 640*Pow(a,2)*Pow(e,3)*f*Pow(g,2)*h*k*Pow(n,2) -
               256*Pow(b,3)*d*e*Pow(g,3)*h*k*Pow(n,2) - 256*Pow(b,3)*c*f*Pow(g,3)*h*k*Pow(n,2) + 1024*a*Pow(b,2)*e*f*Pow(g,3)*h*k*Pow(n,2) -
               48*b*Pow(c,3)*Pow(e,3)*Pow(h,2)*k*Pow(n,2) + 48*a*Pow(c,2)*Pow(e,4)*Pow(h,2)*k*Pow(n,2) + 352*Pow(b,2)*Pow(c,2)*Pow(e,2)*g*Pow(h,2)*k*Pow(n,2) +
               320*a*b*c*Pow(e,3)*g*Pow(h,2)*k*Pow(n,2) - 576*Pow(a,2)*Pow(e,4)*g*Pow(h,2)*k*Pow(n,2) - 768*Pow(b,3)*c*e*Pow(g,2)*Pow(h,2)*k*Pow(n,2) -
               256*a*Pow(b,2)*Pow(e,2)*Pow(g,2)*Pow(h,2)*k*Pow(n,2) + 512*Pow(b,4)*Pow(g,3)*Pow(h,2)*k*Pow(n,2) + 96*b*Pow(c,3)*Pow(e,3)*g*Pow(k,2)*Pow(n,2) -
               96*a*Pow(c,2)*Pow(e,4)*g*Pow(k,2)*Pow(n,2) - 352*Pow(b,2)*Pow(c,2)*Pow(e,2)*Pow(g,2)*Pow(k,2)*Pow(n,2) +
               128*a*b*c*Pow(e,3)*Pow(g,2)*Pow(k,2)*Pow(n,2) + 128*Pow(a,2)*Pow(e,4)*Pow(g,2)*Pow(k,2)*Pow(n,2) + 512*Pow(b,3)*c*e*Pow(g,3)*Pow(k,2)*Pow(n,2) -
               128*a*Pow(b,2)*Pow(e,2)*Pow(g,3)*Pow(k,2)*Pow(n,2) - 256*Pow(b,4)*Pow(g,4)*Pow(k,2)*Pow(n,2) - 36*c*Pow(d,4)*Pow(e,3)*g*l*Pow(n,2) +
               92*Pow(c,2)*Pow(d,3)*Pow(e,2)*f*g*l*Pow(n,2) - 76*Pow(c,3)*Pow(d,2)*e*Pow(f,2)*g*l*Pow(n,2) + 20*Pow(c,4)*d*Pow(f,3)*g*l*Pow(n,2) -
               144*b*Pow(d,4)*Pow(e,2)*Pow(g,2)*l*Pow(n,2) - 32*b*c*Pow(d,3)*e*f*Pow(g,2)*l*Pow(n,2) + 96*a*Pow(d,3)*Pow(e,2)*f*Pow(g,2)*l*Pow(n,2) +
               32*b*Pow(c,2)*Pow(d,2)*Pow(f,2)*Pow(g,2)*l*Pow(n,2) + 176*a*c*Pow(d,2)*e*Pow(f,2)*Pow(g,2)*l*Pow(n,2) -
               128*a*Pow(c,2)*d*Pow(f,3)*Pow(g,2)*l*Pow(n,2) - 256*Pow(b,2)*Pow(d,3)*f*Pow(g,3)*l*Pow(n,2) + 128*a*b*Pow(d,2)*Pow(f,2)*Pow(g,3)*l*Pow(n,2) +
               192*Pow(a,2)*d*Pow(f,3)*Pow(g,3)*l*Pow(n,2) + 40*Pow(c,2)*Pow(d,3)*Pow(e,3)*h*l*Pow(n,2) - 104*Pow(c,3)*Pow(d,2)*Pow(e,2)*f*h*l*Pow(n,2) +
               88*Pow(c,4)*d*e*Pow(f,2)*h*l*Pow(n,2) - 24*Pow(c,5)*Pow(f,3)*h*l*Pow(n,2) + 368*b*c*Pow(d,3)*Pow(e,2)*g*h*l*Pow(n,2) +
               48*a*Pow(d,3)*Pow(e,3)*g*h*l*Pow(n,2) - 16*b*Pow(c,2)*Pow(d,2)*e*f*g*h*l*Pow(n,2) - 464*a*c*Pow(d,2)*Pow(e,2)*f*g*h*l*Pow(n,2) +
               32*b*Pow(c,3)*d*Pow(f,2)*g*h*l*Pow(n,2) - 128*a*Pow(c,2)*d*e*Pow(f,2)*g*h*l*Pow(n,2) + 160*a*Pow(c,3)*Pow(f,3)*g*h*l*Pow(n,2) +
               640*Pow(b,2)*Pow(d,3)*e*Pow(g,2)*h*l*Pow(n,2) + 1024*Pow(b,2)*c*Pow(d,2)*f*Pow(g,2)*h*l*Pow(n,2) - 1088*a*b*Pow(d,2)*e*f*Pow(g,2)*h*l*Pow(n,2) -
               640*a*b*c*d*Pow(f,2)*Pow(g,2)*h*l*Pow(n,2) + 128*Pow(a,2)*d*e*Pow(f,2)*Pow(g,2)*h*l*Pow(n,2) - 256*Pow(a,2)*c*Pow(f,3)*Pow(g,2)*h*l*Pow(n,2) -
               416*b*Pow(c,2)*Pow(d,2)*Pow(e,2)*Pow(h,2)*l*Pow(n,2) + 208*a*c*Pow(d,2)*Pow(e,3)*Pow(h,2)*l*Pow(n,2) + 416*b*Pow(c,3)*d*e*f*Pow(h,2)*l*Pow(n,2) -
               48*a*Pow(c,2)*d*Pow(e,2)*f*Pow(h,2)*l*Pow(n,2) - 144*b*Pow(c,4)*Pow(f,2)*Pow(h,2)*l*Pow(n,2) - 16*a*Pow(c,3)*e*Pow(f,2)*Pow(h,2)*l*Pow(n,2) -
               448*Pow(b,2)*c*Pow(d,2)*e*g*Pow(h,2)*l*Pow(n,2) - 896*a*b*Pow(d,2)*Pow(e,2)*g*Pow(h,2)*l*Pow(n,2) -
               1792*Pow(b,2)*Pow(c,2)*d*f*g*Pow(h,2)*l*Pow(n,2) + 1792*a*b*c*d*e*f*g*Pow(h,2)*l*Pow(n,2) + 896*Pow(a,2)*d*Pow(e,2)*f*g*Pow(h,2)*l*Pow(n,2) +
               896*a*b*Pow(c,2)*Pow(f,2)*g*Pow(h,2)*l*Pow(n,2) - 448*Pow(a,2)*c*e*Pow(f,2)*g*Pow(h,2)*l*Pow(n,2) -
               1280*Pow(b,3)*Pow(d,2)*Pow(g,2)*Pow(h,2)*l*Pow(n,2) + 1280*a*Pow(b,2)*d*f*Pow(g,2)*Pow(h,2)*l*Pow(n,2) -
               256*Pow(a,2)*b*Pow(f,2)*Pow(g,2)*Pow(h,2)*l*Pow(n,2) - 384*Pow(b,2)*Pow(c,2)*d*e*Pow(h,3)*l*Pow(n,2) + 1472*a*b*c*d*Pow(e,2)*Pow(h,3)*l*Pow(n,2) -
               576*Pow(a,2)*d*Pow(e,3)*Pow(h,3)*l*Pow(n,2) + 768*Pow(b,2)*Pow(c,3)*f*Pow(h,3)*l*Pow(n,2) - 2112*a*b*Pow(c,2)*e*f*Pow(h,3)*l*Pow(n,2) +
               960*Pow(a,2)*c*Pow(e,2)*f*Pow(h,3)*l*Pow(n,2) + 2304*Pow(b,3)*c*d*g*Pow(h,3)*l*Pow(n,2) - 256*a*Pow(b,2)*d*e*g*Pow(h,3)*l*Pow(n,2) -
               256*a*Pow(b,2)*c*f*g*Pow(h,3)*l*Pow(n,2) - 256*Pow(a,2)*b*e*f*g*Pow(h,3)*l*Pow(n,2) - 768*Pow(b,3)*Pow(c,2)*Pow(h,4)*l*Pow(n,2) +
               2560*a*Pow(b,2)*c*e*Pow(h,4)*l*Pow(n,2) - 2304*Pow(a,2)*b*Pow(e,2)*Pow(h,4)*l*Pow(n,2) - 2048*a*Pow(b,3)*g*Pow(h,4)*l*Pow(n,2) -
               12*Pow(c,3)*Pow(d,2)*Pow(e,3)*k*l*Pow(n,2) + 24*Pow(c,4)*d*Pow(e,2)*f*k*l*Pow(n,2) - 12*Pow(c,5)*e*Pow(f,2)*k*l*Pow(n,2) -
               104*b*Pow(c,2)*Pow(d,2)*Pow(e,2)*g*k*l*Pow(n,2) + 80*a*c*Pow(d,2)*Pow(e,3)*g*k*l*Pow(n,2) - 288*b*Pow(c,3)*d*e*f*g*k*l*Pow(n,2) +
               240*a*Pow(c,2)*d*Pow(e,2)*f*g*k*l*Pow(n,2) + 104*b*Pow(c,4)*Pow(f,2)*g*k*l*Pow(n,2) - 32*a*Pow(c,3)*e*Pow(f,2)*g*k*l*Pow(n,2) -
               576*Pow(b,2)*c*Pow(d,2)*e*Pow(g,2)*k*l*Pow(n,2) + 992*a*b*Pow(d,2)*Pow(e,2)*Pow(g,2)*k*l*Pow(n,2) +
               320*Pow(b,2)*Pow(c,2)*d*f*Pow(g,2)*k*l*Pow(n,2) + 384*a*b*c*d*e*f*Pow(g,2)*k*l*Pow(n,2) - 576*Pow(a,2)*d*Pow(e,2)*f*Pow(g,2)*k*l*Pow(n,2) -
               576*a*b*Pow(c,2)*Pow(f,2)*Pow(g,2)*k*l*Pow(n,2) + 320*Pow(a,2)*c*e*Pow(f,2)*Pow(g,2)*k*l*Pow(n,2) + 512*Pow(b,3)*Pow(d,2)*Pow(g,3)*k*l*Pow(n,2) -
               512*a*Pow(b,2)*d*f*Pow(g,3)*k*l*Pow(n,2) + 640*Pow(a,2)*b*Pow(f,2)*Pow(g,3)*k*l*Pow(n,2) + 416*b*Pow(c,3)*d*Pow(e,2)*h*k*l*Pow(n,2) -
               368*a*Pow(c,2)*d*Pow(e,3)*h*k*l*Pow(n,2) - 128*b*Pow(c,4)*e*f*h*k*l*Pow(n,2) + 80*a*Pow(c,3)*Pow(e,2)*f*h*k*l*Pow(n,2) +
               704*Pow(b,2)*Pow(c,2)*d*e*g*h*k*l*Pow(n,2) - 1280*a*b*c*d*Pow(e,2)*g*h*k*l*Pow(n,2) - 64*Pow(a,2)*d*Pow(e,3)*g*h*k*l*Pow(n,2) -
               64*Pow(b,2)*Pow(c,3)*f*g*h*k*l*Pow(n,2) + 1408*a*b*Pow(c,2)*e*f*g*h*k*l*Pow(n,2) - 1088*Pow(a,2)*c*Pow(e,2)*f*g*h*k*l*Pow(n,2) -
               768*Pow(b,3)*c*d*Pow(g,2)*h*k*l*Pow(n,2) - 512*a*Pow(b,2)*d*e*Pow(g,2)*h*k*l*Pow(n,2) - 512*a*Pow(b,2)*c*f*Pow(g,2)*h*k*l*Pow(n,2) -
               512*Pow(a,2)*b*e*f*Pow(g,2)*h*k*l*Pow(n,2) - 64*Pow(b,2)*Pow(c,3)*e*Pow(h,2)*k*l*Pow(n,2) - 416*a*b*Pow(c,2)*Pow(e,2)*Pow(h,2)*k*l*Pow(n,2) +
               384*Pow(a,2)*c*Pow(e,3)*Pow(h,2)*k*l*Pow(n,2) + 128*Pow(b,3)*Pow(c,2)*g*Pow(h,2)*k*l*Pow(n,2) - 2816*a*Pow(b,2)*c*e*g*Pow(h,2)*k*l*Pow(n,2) +
               3968*Pow(a,2)*b*Pow(e,2)*g*Pow(h,2)*k*l*Pow(n,2) + 2560*a*Pow(b,3)*Pow(g,2)*Pow(h,2)*k*l*Pow(n,2) - 96*b*Pow(c,4)*Pow(e,2)*Pow(k,2)*l*Pow(n,2) +
               96*a*Pow(c,3)*Pow(e,3)*Pow(k,2)*l*Pow(n,2) + 128*Pow(b,2)*Pow(c,3)*e*g*Pow(k,2)*l*Pow(n,2) - 64*a*b*Pow(c,2)*Pow(e,2)*g*Pow(k,2)*l*Pow(n,2) +
               128*Pow(a,2)*c*Pow(e,3)*g*Pow(k,2)*l*Pow(n,2) - 128*Pow(b,3)*Pow(c,2)*Pow(g,2)*Pow(k,2)*l*Pow(n,2) +
               1024*a*Pow(b,2)*c*e*Pow(g,2)*Pow(k,2)*l*Pow(n,2) - 1280*Pow(a,2)*b*Pow(e,2)*Pow(g,2)*Pow(k,2)*l*Pow(n,2) -
               512*a*Pow(b,3)*Pow(g,3)*Pow(k,2)*l*Pow(n,2) + 8*Pow(c,2)*Pow(d,4)*Pow(e,2)*Pow(l,2)*Pow(n,2) - 16*Pow(c,3)*Pow(d,3)*e*f*Pow(l,2)*Pow(n,2) +
               8*Pow(c,4)*Pow(d,2)*Pow(f,2)*Pow(l,2)*Pow(n,2) + 160*b*c*Pow(d,4)*e*g*Pow(l,2)*Pow(n,2) + 24*a*Pow(d,4)*Pow(e,2)*g*Pow(l,2)*Pow(n,2) -
               16*b*Pow(c,2)*Pow(d,3)*f*g*Pow(l,2)*Pow(n,2) - 256*a*c*Pow(d,3)*e*f*g*Pow(l,2)*Pow(n,2) + 88*a*Pow(c,2)*Pow(d,2)*Pow(f,2)*g*Pow(l,2)*Pow(n,2) +
               128*Pow(b,2)*Pow(d,4)*Pow(g,2)*Pow(l,2)*Pow(n,2) + 128*a*b*Pow(d,3)*f*Pow(g,2)*Pow(l,2)*Pow(n,2) -
               352*Pow(a,2)*Pow(d,2)*Pow(f,2)*Pow(g,2)*Pow(l,2)*Pow(n,2) - 96*b*Pow(c,2)*Pow(d,3)*e*h*Pow(l,2)*Pow(n,2) -
               208*a*c*Pow(d,3)*Pow(e,2)*h*Pow(l,2)*Pow(n,2) - 96*b*Pow(c,3)*Pow(d,2)*f*h*Pow(l,2)*Pow(n,2) + 608*a*Pow(c,2)*Pow(d,2)*e*f*h*Pow(l,2)*Pow(n,2) -
               208*a*Pow(c,3)*d*Pow(f,2)*h*Pow(l,2)*Pow(n,2) - 1152*Pow(b,2)*c*Pow(d,3)*g*h*Pow(l,2)*Pow(n,2) - 64*a*b*Pow(d,3)*e*g*h*Pow(l,2)*Pow(n,2) +
               704*a*b*c*Pow(d,2)*f*g*h*Pow(l,2)*Pow(n,2) + 128*Pow(a,2)*Pow(d,2)*e*f*g*h*Pow(l,2)*Pow(n,2) + 576*Pow(a,2)*c*d*Pow(f,2)*g*h*Pow(l,2)*Pow(n,2) +
               1344*Pow(b,2)*Pow(c,2)*Pow(d,2)*Pow(h,2)*Pow(l,2)*Pow(n,2) - 448*a*b*c*Pow(d,2)*e*Pow(h,2)*Pow(l,2)*Pow(n,2) +
               224*Pow(a,2)*Pow(d,2)*Pow(e,2)*Pow(h,2)*Pow(l,2)*Pow(n,2) - 448*a*b*Pow(c,2)*d*f*Pow(h,2)*Pow(l,2)*Pow(n,2) -
               896*Pow(a,2)*c*d*e*f*Pow(h,2)*Pow(l,2)*Pow(n,2) + 224*Pow(a,2)*Pow(c,2)*Pow(f,2)*Pow(h,2)*Pow(l,2)*Pow(n,2) +
               2432*a*Pow(b,2)*Pow(d,2)*g*Pow(h,2)*Pow(l,2)*Pow(n,2) - 2304*Pow(a,2)*b*d*f*g*Pow(h,2)*Pow(l,2)*Pow(n,2) -
               384*Pow(a,3)*Pow(f,2)*g*Pow(h,2)*Pow(l,2)*Pow(n,2) - 4608*a*Pow(b,2)*c*d*Pow(h,3)*Pow(l,2)*Pow(n,2) +
               2560*Pow(a,2)*b*d*e*Pow(h,3)*Pow(l,2)*Pow(n,2) + 2560*Pow(a,2)*b*c*f*Pow(h,3)*Pow(l,2)*Pow(n,2) - 1536*Pow(a,3)*e*f*Pow(h,3)*Pow(l,2)*Pow(n,2) +
               2048*Pow(a,2)*Pow(b,2)*Pow(h,4)*Pow(l,2)*Pow(n,2) - 16*b*Pow(c,3)*Pow(d,2)*e*k*Pow(l,2)*Pow(n,2) +
               88*a*Pow(c,2)*Pow(d,2)*Pow(e,2)*k*Pow(l,2)*Pow(n,2) + 160*b*Pow(c,4)*d*f*k*Pow(l,2)*Pow(n,2) - 256*a*Pow(c,3)*d*e*f*k*Pow(l,2)*Pow(n,2) +
               24*a*Pow(c,4)*Pow(f,2)*k*Pow(l,2)*Pow(n,2) + 608*Pow(b,2)*Pow(c,2)*Pow(d,2)*g*k*Pow(l,2)*Pow(n,2) - 704*a*b*c*Pow(d,2)*e*g*k*Pow(l,2)*Pow(n,2) -
               64*Pow(a,2)*Pow(d,2)*Pow(e,2)*g*k*Pow(l,2)*Pow(n,2) - 704*a*b*Pow(c,2)*d*f*g*k*Pow(l,2)*Pow(n,2) + 640*Pow(a,2)*c*d*e*f*g*k*Pow(l,2)*Pow(n,2) -
               64*Pow(a,2)*Pow(c,2)*Pow(f,2)*g*k*Pow(l,2)*Pow(n,2) - 1280*a*Pow(b,2)*Pow(d,2)*Pow(g,2)*k*Pow(l,2)*Pow(n,2) +
               1024*Pow(a,2)*b*d*f*Pow(g,2)*k*Pow(l,2)*Pow(n,2) - 128*Pow(a,3)*Pow(f,2)*Pow(g,2)*k*Pow(l,2)*Pow(n,2) -
               1152*Pow(b,2)*Pow(c,3)*d*h*k*Pow(l,2)*Pow(n,2) + 704*a*b*Pow(c,2)*d*e*h*k*Pow(l,2)*Pow(n,2) + 576*Pow(a,2)*c*d*Pow(e,2)*h*k*Pow(l,2)*Pow(n,2) -
               64*a*b*Pow(c,3)*f*h*k*Pow(l,2)*Pow(n,2) + 128*Pow(a,2)*Pow(c,2)*e*f*h*k*Pow(l,2)*Pow(n,2) + 3072*a*Pow(b,2)*c*d*g*h*k*Pow(l,2)*Pow(n,2) -
               512*Pow(a,2)*b*d*e*g*h*k*Pow(l,2)*Pow(n,2) - 512*Pow(a,2)*b*c*f*g*h*k*Pow(l,2)*Pow(n,2) + 1024*Pow(a,3)*e*f*g*h*k*Pow(l,2)*Pow(n,2) +
               2432*a*Pow(b,2)*Pow(c,2)*Pow(h,2)*k*Pow(l,2)*Pow(n,2) - 2304*Pow(a,2)*b*c*e*Pow(h,2)*k*Pow(l,2)*Pow(n,2) -
               384*Pow(a,3)*Pow(e,2)*Pow(h,2)*k*Pow(l,2)*Pow(n,2) - 5120*Pow(a,2)*Pow(b,2)*g*Pow(h,2)*k*Pow(l,2)*Pow(n,2) +
               128*Pow(b,2)*Pow(c,4)*Pow(k,2)*Pow(l,2)*Pow(n,2) + 128*a*b*Pow(c,3)*e*Pow(k,2)*Pow(l,2)*Pow(n,2) -
               352*Pow(a,2)*Pow(c,2)*Pow(e,2)*Pow(k,2)*Pow(l,2)*Pow(n,2) - 1280*a*Pow(b,2)*Pow(c,2)*g*Pow(k,2)*Pow(l,2)*Pow(n,2) +
               1024*Pow(a,2)*b*c*e*g*Pow(k,2)*Pow(l,2)*Pow(n,2) - 128*Pow(a,3)*Pow(e,2)*g*Pow(k,2)*Pow(l,2)*Pow(n,2) +
               1536*Pow(a,2)*Pow(b,2)*Pow(g,2)*Pow(k,2)*Pow(l,2)*Pow(n,2) - 48*b*Pow(c,2)*Pow(d,4)*Pow(l,3)*Pow(n,2) + 16*a*c*Pow(d,4)*e*Pow(l,3)*Pow(n,2) +
               32*a*Pow(c,2)*Pow(d,3)*f*Pow(l,3)*Pow(n,2) - 128*a*b*Pow(d,4)*g*Pow(l,3)*Pow(n,2) + 192*Pow(a,2)*Pow(d,3)*f*g*Pow(l,3)*Pow(n,2) +
               576*a*b*c*Pow(d,3)*h*Pow(l,3)*Pow(n,2) + 64*Pow(a,2)*Pow(d,3)*e*h*Pow(l,3)*Pow(n,2) - 704*Pow(a,2)*c*Pow(d,2)*f*h*Pow(l,3)*Pow(n,2) -
               1280*Pow(a,2)*b*Pow(d,2)*Pow(h,2)*Pow(l,3)*Pow(n,2) + 1536*Pow(a,3)*d*f*Pow(h,2)*Pow(l,3)*Pow(n,2) +
               32*a*b*Pow(c,2)*Pow(d,2)*k*Pow(l,3)*Pow(n,2) - 192*Pow(a,2)*c*Pow(d,2)*e*k*Pow(l,3)*Pow(n,2) + 256*Pow(a,2)*Pow(c,2)*d*f*k*Pow(l,3)*Pow(n,2) +
               640*Pow(a,2)*b*Pow(d,2)*g*k*Pow(l,3)*Pow(n,2) - 768*Pow(a,3)*d*f*g*k*Pow(l,3)*Pow(n,2) - 768*Pow(a,2)*b*c*d*h*k*Pow(l,3)*Pow(n,2) -
               256*Pow(a,3)*d*e*h*k*Pow(l,3)*Pow(n,2) - 256*Pow(a,3)*c*f*h*k*Pow(l,3)*Pow(n,2) + 2048*Pow(a,3)*b*Pow(h,2)*k*Pow(l,3)*Pow(n,2) -
               128*Pow(a,2)*b*Pow(c,2)*Pow(k,2)*Pow(l,3)*Pow(n,2) + 512*Pow(a,3)*c*e*Pow(k,2)*Pow(l,3)*Pow(n,2) - 512*Pow(a,3)*b*g*Pow(k,2)*Pow(l,3)*Pow(n,2) -
               16*Pow(a,2)*Pow(d,4)*Pow(l,4)*Pow(n,2) + 128*Pow(a,3)*Pow(d,2)*k*Pow(l,4)*Pow(n,2) - 256*Pow(a,4)*Pow(k,2)*Pow(l,4)*Pow(n,2) +
               12*Pow(c,2)*Pow(d,3)*Pow(e,3)*g*m*Pow(n,2) - 20*Pow(c,3)*Pow(d,2)*Pow(e,2)*f*g*m*Pow(n,2) + 4*Pow(c,4)*d*e*Pow(f,2)*g*m*Pow(n,2) +
               4*Pow(c,5)*Pow(f,3)*g*m*Pow(n,2) + 240*b*c*Pow(d,3)*Pow(e,2)*Pow(g,2)*m*Pow(n,2) - 144*a*Pow(d,3)*Pow(e,3)*Pow(g,2)*m*Pow(n,2) +
               48*b*Pow(c,2)*Pow(d,2)*e*f*Pow(g,2)*m*Pow(n,2) - 176*a*c*Pow(d,2)*Pow(e,2)*f*Pow(g,2)*m*Pow(n,2) - 96*b*Pow(c,3)*d*Pow(f,2)*Pow(g,2)*m*Pow(n,2) +
               160*a*Pow(c,2)*d*e*Pow(f,2)*Pow(g,2)*m*Pow(n,2) - 32*a*Pow(c,3)*Pow(f,3)*Pow(g,2)*m*Pow(n,2) - 384*Pow(b,2)*Pow(d,3)*e*Pow(g,3)*m*Pow(n,2) -
               256*Pow(b,2)*c*Pow(d,2)*f*Pow(g,3)*m*Pow(n,2) + 832*a*b*Pow(d,2)*e*f*Pow(g,3)*m*Pow(n,2) + 384*a*b*c*d*Pow(f,2)*Pow(g,3)*m*Pow(n,2) -
               704*Pow(a,2)*d*e*Pow(f,2)*Pow(g,3)*m*Pow(n,2) + 64*Pow(a,2)*c*Pow(f,3)*Pow(g,3)*m*Pow(n,2) + 8*Pow(c,3)*Pow(d,2)*Pow(e,3)*h*m*Pow(n,2) -
               16*Pow(c,4)*d*Pow(e,2)*f*h*m*Pow(n,2) + 8*Pow(c,5)*e*Pow(f,2)*h*m*Pow(n,2) - 48*b*Pow(c,2)*Pow(d,2)*Pow(e,2)*g*h*m*Pow(n,2) -
               256*a*c*Pow(d,2)*Pow(e,3)*g*h*m*Pow(n,2) - 288*b*Pow(c,3)*d*e*f*g*h*m*Pow(n,2) + 800*a*Pow(c,2)*d*Pow(e,2)*f*g*h*m*Pow(n,2) +
               48*b*Pow(c,4)*Pow(f,2)*g*h*m*Pow(n,2) - 256*a*Pow(c,3)*e*Pow(f,2)*g*h*m*Pow(n,2) - 896*Pow(b,2)*c*Pow(d,2)*e*Pow(g,2)*h*m*Pow(n,2) +
               896*a*b*Pow(d,2)*Pow(e,2)*Pow(g,2)*h*m*Pow(n,2) + 896*Pow(b,2)*Pow(c,2)*d*f*Pow(g,2)*h*m*Pow(n,2) - 896*a*b*c*d*e*f*Pow(g,2)*h*m*Pow(n,2) -
               896*Pow(a,2)*d*Pow(e,2)*f*Pow(g,2)*h*m*Pow(n,2) + 896*Pow(a,2)*c*e*Pow(f,2)*Pow(g,2)*h*m*Pow(n,2) + 1536*Pow(b,3)*Pow(d,2)*Pow(g,3)*h*m*Pow(n,2) -
               1536*a*Pow(b,2)*d*f*Pow(g,3)*h*m*Pow(n,2) - 768*Pow(a,2)*b*Pow(f,2)*Pow(g,3)*h*m*Pow(n,2) - 416*b*Pow(c,3)*d*Pow(e,2)*Pow(h,2)*m*Pow(n,2) +
               368*a*Pow(c,2)*d*Pow(e,3)*Pow(h,2)*m*Pow(n,2) + 128*b*Pow(c,4)*e*f*Pow(h,2)*m*Pow(n,2) - 80*a*Pow(c,3)*Pow(e,2)*f*Pow(h,2)*m*Pow(n,2) +
               2432*Pow(b,2)*Pow(c,2)*d*e*g*Pow(h,2)*m*Pow(n,2) - 1856*a*b*c*d*Pow(e,2)*g*Pow(h,2)*m*Pow(n,2) + 960*Pow(a,2)*d*Pow(e,3)*g*Pow(h,2)*m*Pow(n,2) -
               384*Pow(b,2)*Pow(c,3)*f*g*Pow(h,2)*m*Pow(n,2) - 64*a*b*Pow(c,2)*e*f*g*Pow(h,2)*m*Pow(n,2) - 704*Pow(a,2)*c*Pow(e,2)*f*g*Pow(h,2)*m*Pow(n,2) -
               2816*Pow(b,3)*c*d*Pow(g,2)*Pow(h,2)*m*Pow(n,2) + 512*a*Pow(b,2)*d*e*Pow(g,2)*Pow(h,2)*m*Pow(n,2) +
               512*a*Pow(b,2)*c*f*Pow(g,2)*Pow(h,2)*m*Pow(n,2) + 2304*Pow(a,2)*b*e*f*Pow(g,2)*Pow(h,2)*m*Pow(n,2) - 256*Pow(b,2)*Pow(c,3)*e*Pow(h,3)*m*Pow(n,2) +
               1472*a*b*Pow(c,2)*Pow(e,2)*Pow(h,3)*m*Pow(n,2) - 1152*Pow(a,2)*c*Pow(e,3)*Pow(h,3)*m*Pow(n,2) + 512*Pow(b,3)*Pow(c,2)*g*Pow(h,3)*m*Pow(n,2) -
               4096*a*Pow(b,2)*c*e*g*Pow(h,3)*m*Pow(n,2) + 1536*Pow(a,2)*b*Pow(e,2)*g*Pow(h,3)*m*Pow(n,2) + 3072*a*Pow(b,3)*Pow(g,2)*Pow(h,3)*m*Pow(n,2) +
               80*b*Pow(c,3)*d*Pow(e,2)*g*k*m*Pow(n,2) - 32*a*Pow(c,2)*d*Pow(e,3)*g*k*m*Pow(n,2) + 208*b*Pow(c,4)*e*f*g*k*m*Pow(n,2) -
               256*a*Pow(c,3)*Pow(e,2)*f*g*k*m*Pow(n,2) + 128*Pow(b,2)*Pow(c,2)*d*e*Pow(g,2)*k*m*Pow(n,2) - 1088*a*b*c*d*Pow(e,2)*Pow(g,2)*k*m*Pow(n,2) +
               640*Pow(a,2)*d*Pow(e,3)*Pow(g,2)*k*m*Pow(n,2) - 256*Pow(b,2)*Pow(c,3)*f*Pow(g,2)*k*m*Pow(n,2) - 640*a*b*Pow(c,2)*e*f*Pow(g,2)*k*m*Pow(n,2) +
               1024*Pow(a,2)*c*Pow(e,2)*f*Pow(g,2)*k*m*Pow(n,2) - 256*Pow(b,3)*c*d*Pow(g,3)*k*m*Pow(n,2) + 1024*a*Pow(b,2)*d*e*Pow(g,3)*k*m*Pow(n,2) +
               1024*a*Pow(b,2)*c*f*Pow(g,3)*k*m*Pow(n,2) - 768*Pow(a,2)*b*e*f*Pow(g,3)*k*m*Pow(n,2) + 96*b*Pow(c,4)*Pow(e,2)*h*k*m*Pow(n,2) -
               96*a*Pow(c,3)*Pow(e,3)*h*k*m*Pow(n,2) - 1024*Pow(b,2)*Pow(c,3)*e*g*h*k*m*Pow(n,2) + 960*a*b*Pow(c,2)*Pow(e,2)*g*h*k*m*Pow(n,2) -
               128*Pow(a,2)*c*Pow(e,3)*g*h*k*m*Pow(n,2) + 1024*Pow(b,3)*Pow(c,2)*Pow(g,2)*h*k*m*Pow(n,2) + 2560*a*Pow(b,2)*c*e*Pow(g,2)*h*k*m*Pow(n,2) -
               2304*Pow(a,2)*b*Pow(e,2)*Pow(g,2)*h*k*m*Pow(n,2) - 3072*a*Pow(b,3)*Pow(g,3)*h*k*m*Pow(n,2) - 16*Pow(c,3)*Pow(d,3)*Pow(e,2)*l*m*Pow(n,2) +
               32*Pow(c,4)*Pow(d,2)*e*f*l*m*Pow(n,2) - 16*Pow(c,5)*d*Pow(f,2)*l*m*Pow(n,2) - 528*b*Pow(c,2)*Pow(d,3)*e*g*l*m*Pow(n,2) +
               368*a*c*Pow(d,3)*Pow(e,2)*g*l*m*Pow(n,2) + 144*b*Pow(c,3)*Pow(d,2)*f*g*l*m*Pow(n,2) - 16*a*Pow(c,2)*Pow(d,2)*e*f*g*l*m*Pow(n,2) +
               32*a*Pow(c,3)*d*Pow(f,2)*g*l*m*Pow(n,2) + 640*Pow(b,2)*c*Pow(d,3)*Pow(g,2)*l*m*Pow(n,2) - 64*a*b*Pow(d,3)*e*Pow(g,2)*l*m*Pow(n,2) -
               1088*a*b*c*Pow(d,2)*f*Pow(g,2)*l*m*Pow(n,2) + 576*Pow(a,2)*Pow(d,2)*e*f*Pow(g,2)*l*m*Pow(n,2) + 128*Pow(a,2)*c*d*Pow(f,2)*Pow(g,2)*l*m*Pow(n,2) +
               416*b*Pow(c,3)*Pow(d,2)*e*h*l*m*Pow(n,2) - 160*a*Pow(c,2)*Pow(d,2)*Pow(e,2)*h*l*m*Pow(n,2) - 128*b*Pow(c,4)*d*f*h*l*m*Pow(n,2) -
               288*a*Pow(c,3)*d*e*f*h*l*m*Pow(n,2) + 160*a*Pow(c,4)*Pow(f,2)*h*l*m*Pow(n,2) - 448*Pow(b,2)*Pow(c,2)*Pow(d,2)*g*h*l*m*Pow(n,2) +
               1792*a*b*c*Pow(d,2)*e*g*h*l*m*Pow(n,2) - 448*Pow(a,2)*Pow(d,2)*Pow(e,2)*g*h*l*m*Pow(n,2) + 896*a*b*Pow(c,2)*d*f*g*h*l*m*Pow(n,2) -
               896*Pow(a,2)*c*d*e*f*g*h*l*m*Pow(n,2) - 896*Pow(a,2)*Pow(c,2)*Pow(f,2)*g*h*l*m*Pow(n,2) - 2304*a*Pow(b,2)*Pow(d,2)*Pow(g,2)*h*l*m*Pow(n,2) +
               2560*Pow(a,2)*b*d*f*Pow(g,2)*h*l*m*Pow(n,2) + 1024*Pow(a,3)*Pow(f,2)*Pow(g,2)*h*l*m*Pow(n,2) - 384*Pow(b,2)*Pow(c,3)*d*Pow(h,2)*l*m*Pow(n,2) -
               64*a*b*Pow(c,2)*d*e*Pow(h,2)*l*m*Pow(n,2) - 704*Pow(a,2)*c*d*Pow(e,2)*Pow(h,2)*l*m*Pow(n,2) + 576*a*b*Pow(c,3)*f*Pow(h,2)*l*m*Pow(n,2) +
               192*Pow(a,2)*Pow(c,2)*e*f*Pow(h,2)*l*m*Pow(n,2) + 2816*a*Pow(b,2)*c*d*g*Pow(h,2)*l*m*Pow(n,2) - 4352*Pow(a,2)*b*d*e*g*Pow(h,2)*l*m*Pow(n,2) -
               4352*Pow(a,2)*b*c*f*g*Pow(h,2)*l*m*Pow(n,2) + 3328*Pow(a,3)*e*f*g*Pow(h,2)*l*m*Pow(n,2) - 256*a*Pow(b,2)*Pow(c,2)*Pow(h,3)*l*m*Pow(n,2) -
               512*Pow(a,2)*b*c*e*Pow(h,3)*l*m*Pow(n,2) + 2304*Pow(a,3)*Pow(e,2)*Pow(h,3)*l*m*Pow(n,2) + 2048*Pow(a,2)*Pow(b,2)*g*Pow(h,3)*l*m*Pow(n,2) -
               128*b*Pow(c,4)*d*e*k*l*m*Pow(n,2) + 80*a*Pow(c,3)*d*Pow(e,2)*k*l*m*Pow(n,2) - 160*b*Pow(c,5)*f*k*l*m*Pow(n,2) + 208*a*Pow(c,4)*e*f*k*l*m*Pow(n,2) -
               64*Pow(b,2)*Pow(c,3)*d*g*k*l*m*Pow(n,2) + 1408*a*b*Pow(c,2)*d*e*g*k*l*m*Pow(n,2) - 1088*Pow(a,2)*c*d*Pow(e,2)*g*k*l*m*Pow(n,2) +
               768*a*b*Pow(c,3)*f*g*k*l*m*Pow(n,2) - 640*Pow(a,2)*Pow(c,2)*e*f*g*k*l*m*Pow(n,2) - 512*a*Pow(b,2)*c*d*Pow(g,2)*k*l*m*Pow(n,2) -
               512*Pow(a,2)*b*d*e*Pow(g,2)*k*l*m*Pow(n,2) - 512*Pow(a,2)*b*c*f*Pow(g,2)*k*l*m*Pow(n,2) - 768*Pow(a,3)*e*f*Pow(g,2)*k*l*m*Pow(n,2) +
               640*Pow(b,2)*Pow(c,4)*h*k*l*m*Pow(n,2) - 1152*a*b*Pow(c,3)*e*h*k*l*m*Pow(n,2) + 704*Pow(a,2)*Pow(c,2)*Pow(e,2)*h*k*l*m*Pow(n,2) -
               2816*a*Pow(b,2)*Pow(c,2)*g*h*k*l*m*Pow(n,2) + 1536*Pow(a,2)*b*c*e*g*h*k*l*m*Pow(n,2) + 256*Pow(a,3)*Pow(e,2)*g*h*k*l*m*Pow(n,2) +
               4096*Pow(a,2)*Pow(b,2)*Pow(g,2)*h*k*l*m*Pow(n,2) + 192*b*Pow(c,3)*Pow(d,3)*Pow(l,2)*m*Pow(n,2) - 96*a*Pow(c,2)*Pow(d,3)*e*Pow(l,2)*m*Pow(n,2) -
               96*a*Pow(c,3)*Pow(d,2)*f*Pow(l,2)*m*Pow(n,2) - 64*a*b*c*Pow(d,3)*g*Pow(l,2)*m*Pow(n,2) - 256*Pow(a,2)*Pow(d,3)*e*g*Pow(l,2)*m*Pow(n,2) +
               128*Pow(a,2)*c*Pow(d,2)*f*g*Pow(l,2)*m*Pow(n,2) - 1792*a*b*Pow(c,2)*Pow(d,2)*h*Pow(l,2)*m*Pow(n,2) +
               896*Pow(a,2)*c*Pow(d,2)*e*h*Pow(l,2)*m*Pow(n,2) + 896*Pow(a,2)*Pow(c,2)*d*f*h*Pow(l,2)*m*Pow(n,2) +
               1280*Pow(a,2)*b*Pow(d,2)*g*h*Pow(l,2)*m*Pow(n,2) - 1536*Pow(a,3)*d*f*g*h*Pow(l,2)*m*Pow(n,2) + 4096*Pow(a,2)*b*c*d*Pow(h,2)*Pow(l,2)*m*Pow(n,2) -
               1024*Pow(a,3)*d*e*Pow(h,2)*Pow(l,2)*m*Pow(n,2) - 1024*Pow(a,3)*c*f*Pow(h,2)*Pow(l,2)*m*Pow(n,2) - 4096*Pow(a,3)*b*Pow(h,3)*Pow(l,2)*m*Pow(n,2) -
               64*a*b*Pow(c,3)*d*k*Pow(l,2)*m*Pow(n,2) + 128*Pow(a,2)*Pow(c,2)*d*e*k*Pow(l,2)*m*Pow(n,2) - 256*Pow(a,2)*Pow(c,3)*f*k*Pow(l,2)*m*Pow(n,2) -
               512*Pow(a,2)*b*c*d*g*k*Pow(l,2)*m*Pow(n,2) + 1024*Pow(a,3)*d*e*g*k*Pow(l,2)*m*Pow(n,2) + 1024*Pow(a,3)*c*f*g*k*Pow(l,2)*m*Pow(n,2) +
               1280*Pow(a,2)*b*Pow(c,2)*h*k*Pow(l,2)*m*Pow(n,2) - 1536*Pow(a,3)*c*e*h*k*Pow(l,2)*m*Pow(n,2) - 2048*Pow(a,3)*b*g*h*k*Pow(l,2)*m*Pow(n,2) +
               64*Pow(a,2)*c*Pow(d,3)*Pow(l,3)*m*Pow(n,2) - 256*Pow(a,3)*Pow(d,2)*h*Pow(l,3)*m*Pow(n,2) - 256*Pow(a,3)*c*d*k*Pow(l,3)*m*Pow(n,2) +
               1024*Pow(a,4)*h*k*Pow(l,3)*m*Pow(n,2) - 4*Pow(c,4)*Pow(d,2)*Pow(e,2)*Pow(m,2)*Pow(n,2) + 8*Pow(c,5)*d*e*f*Pow(m,2)*Pow(n,2) -
               4*Pow(c,6)*Pow(f,2)*Pow(m,2)*Pow(n,2) - 16*b*Pow(c,3)*Pow(d,2)*e*g*Pow(m,2)*Pow(n,2) + 88*a*Pow(c,2)*Pow(d,2)*Pow(e,2)*g*Pow(m,2)*Pow(n,2) +
               160*b*Pow(c,4)*d*f*g*Pow(m,2)*Pow(n,2) - 256*a*Pow(c,3)*d*e*f*g*Pow(m,2)*Pow(n,2) + 24*a*Pow(c,4)*Pow(f,2)*g*Pow(m,2)*Pow(n,2) +
               224*Pow(b,2)*Pow(c,2)*Pow(d,2)*Pow(g,2)*Pow(m,2)*Pow(n,2) - 448*a*b*c*Pow(d,2)*e*Pow(g,2)*Pow(m,2)*Pow(n,2) +
               224*Pow(a,2)*Pow(d,2)*Pow(e,2)*Pow(g,2)*Pow(m,2)*Pow(n,2) - 896*a*b*Pow(c,2)*d*f*Pow(g,2)*Pow(m,2)*Pow(n,2) +
               896*Pow(a,2)*c*d*e*f*Pow(g,2)*Pow(m,2)*Pow(n,2) - 384*a*Pow(b,2)*Pow(d,2)*Pow(g,3)*Pow(m,2)*Pow(n,2) +
               1024*Pow(a,2)*b*d*f*Pow(g,3)*Pow(m,2)*Pow(n,2) - 128*Pow(a,3)*Pow(f,2)*Pow(g,3)*Pow(m,2)*Pow(n,2) + 352*b*Pow(c,4)*d*e*h*Pow(m,2)*Pow(n,2) -
               304*a*Pow(c,3)*d*Pow(e,2)*h*Pow(m,2)*Pow(n,2) - 64*b*Pow(c,5)*f*h*Pow(m,2)*Pow(n,2) + 16*a*Pow(c,4)*e*f*h*Pow(m,2)*Pow(n,2) -
               832*Pow(b,2)*Pow(c,3)*d*g*h*Pow(m,2)*Pow(n,2) - 512*a*b*Pow(c,2)*d*e*g*h*Pow(m,2)*Pow(n,2) + 192*Pow(a,2)*c*d*Pow(e,2)*g*h*Pow(m,2)*Pow(n,2) +
               128*a*b*Pow(c,3)*f*g*h*Pow(m,2)*Pow(n,2) + 640*Pow(a,2)*Pow(c,2)*e*f*g*h*Pow(m,2)*Pow(n,2) + 2304*a*Pow(b,2)*c*d*Pow(g,2)*h*Pow(m,2)*Pow(n,2) +
               512*Pow(a,2)*b*d*e*Pow(g,2)*h*Pow(m,2)*Pow(n,2) + 512*Pow(a,2)*b*c*f*Pow(g,2)*h*Pow(m,2)*Pow(n,2) -
               2816*Pow(a,3)*e*f*Pow(g,2)*h*Pow(m,2)*Pow(n,2) + 128*Pow(b,2)*Pow(c,4)*Pow(h,2)*Pow(m,2)*Pow(n,2) -
               1216*a*b*Pow(c,3)*e*Pow(h,2)*Pow(m,2)*Pow(n,2) + 992*Pow(a,2)*Pow(c,2)*Pow(e,2)*Pow(h,2)*Pow(m,2)*Pow(n,2) +
               1408*a*Pow(b,2)*Pow(c,2)*g*Pow(h,2)*Pow(m,2)*Pow(n,2) + 2816*Pow(a,2)*b*c*e*g*Pow(h,2)*Pow(m,2)*Pow(n,2) -
               1920*Pow(a,3)*Pow(e,2)*g*Pow(h,2)*Pow(m,2)*Pow(n,2) - 5632*Pow(a,2)*Pow(b,2)*Pow(g,2)*Pow(h,2)*Pow(m,2)*Pow(n,2) -
               48*b*Pow(c,5)*e*k*Pow(m,2)*Pow(n,2) + 48*a*Pow(c,4)*Pow(e,2)*k*Pow(m,2)*Pow(n,2) + 96*Pow(b,2)*Pow(c,4)*g*k*Pow(m,2)*Pow(n,2) -
               128*a*b*Pow(c,3)*e*g*k*Pow(m,2)*Pow(n,2) + 128*Pow(a,2)*Pow(c,2)*Pow(e,2)*g*k*Pow(m,2)*Pow(n,2) -
               256*a*Pow(b,2)*Pow(c,2)*Pow(g,2)*k*Pow(m,2)*Pow(n,2) + 1280*Pow(a,2)*b*c*e*Pow(g,2)*k*Pow(m,2)*Pow(n,2) -
               1280*Pow(a,3)*Pow(e,2)*Pow(g,2)*k*Pow(m,2)*Pow(n,2) - 512*Pow(a,2)*Pow(b,2)*Pow(g,3)*k*Pow(m,2)*Pow(n,2) -
               144*b*Pow(c,4)*Pow(d,2)*l*Pow(m,2)*Pow(n,2) + 96*a*Pow(c,3)*Pow(d,2)*e*l*Pow(m,2)*Pow(n,2) + 48*a*Pow(c,4)*d*f*l*Pow(m,2)*Pow(n,2) +
               896*a*b*Pow(c,2)*Pow(d,2)*g*l*Pow(m,2)*Pow(n,2) - 896*Pow(a,2)*c*Pow(d,2)*e*g*l*Pow(m,2)*Pow(n,2) -
               256*Pow(a,2)*b*Pow(d,2)*Pow(g,2)*l*Pow(m,2)*Pow(n,2) - 768*Pow(a,3)*d*f*Pow(g,2)*l*Pow(m,2)*Pow(n,2) + 128*a*b*Pow(c,3)*d*h*l*Pow(m,2)*Pow(n,2) +
               640*Pow(a,2)*Pow(c,2)*d*e*h*l*Pow(m,2)*Pow(n,2) - 384*Pow(a,2)*Pow(c,3)*f*h*l*Pow(m,2)*Pow(n,2) - 2560*Pow(a,2)*b*c*d*g*h*l*Pow(m,2)*Pow(n,2) +
               1536*Pow(a,3)*d*e*g*h*l*Pow(m,2)*Pow(n,2) + 1536*Pow(a,3)*c*f*g*h*l*Pow(m,2)*Pow(n,2) + 512*Pow(a,2)*b*Pow(c,2)*Pow(h,2)*l*Pow(m,2)*Pow(n,2) -
               2048*Pow(a,3)*c*e*Pow(h,2)*l*Pow(m,2)*Pow(n,2) + 2048*Pow(a,3)*b*g*Pow(h,2)*l*Pow(m,2)*Pow(n,2) + 416*a*b*Pow(c,4)*k*l*Pow(m,2)*Pow(n,2) -
               512*Pow(a,2)*Pow(c,3)*e*k*l*Pow(m,2)*Pow(n,2) - 2304*Pow(a,2)*b*Pow(c,2)*g*k*l*Pow(m,2)*Pow(n,2) + 2048*Pow(a,3)*c*e*g*k*l*Pow(m,2)*Pow(n,2) +
               2560*Pow(a,3)*b*Pow(g,2)*k*l*Pow(m,2)*Pow(n,2) + 512*Pow(a,3)*Pow(d,2)*g*Pow(l,2)*Pow(m,2)*Pow(n,2) -
               1024*Pow(a,3)*c*d*h*Pow(l,2)*Pow(m,2)*Pow(n,2) + 2048*Pow(a,4)*Pow(h,2)*Pow(l,2)*Pow(m,2)*Pow(n,2) +
               512*Pow(a,3)*Pow(c,2)*k*Pow(l,2)*Pow(m,2)*Pow(n,2) - 2048*Pow(a,4)*g*k*Pow(l,2)*Pow(m,2)*Pow(n,2) - 96*b*Pow(c,5)*d*Pow(m,3)*Pow(n,2) +
               80*a*Pow(c,4)*d*e*Pow(m,3)*Pow(n,2) + 16*a*Pow(c,5)*f*Pow(m,3)*Pow(n,2) + 640*a*b*Pow(c,3)*d*g*Pow(m,3)*Pow(n,2) -
               384*Pow(a,2)*Pow(c,2)*d*e*g*Pow(m,3)*Pow(n,2) - 128*Pow(a,2)*Pow(c,3)*f*g*Pow(m,3)*Pow(n,2) - 1024*Pow(a,2)*b*c*d*Pow(g,2)*Pow(m,3)*Pow(n,2) +
               256*Pow(a,3)*d*e*Pow(g,2)*Pow(m,3)*Pow(n,2) + 256*Pow(a,3)*c*f*Pow(g,2)*Pow(m,3)*Pow(n,2) + 320*a*b*Pow(c,4)*h*Pow(m,3)*Pow(n,2) -
               256*Pow(a,2)*Pow(c,3)*e*h*Pow(m,3)*Pow(n,2) - 2048*Pow(a,2)*b*Pow(c,2)*g*h*Pow(m,3)*Pow(n,2) + 1024*Pow(a,3)*c*e*g*h*Pow(m,3)*Pow(n,2) +
               3072*Pow(a,3)*b*Pow(g,2)*h*Pow(m,3)*Pow(n,2) - 128*Pow(a,2)*Pow(c,3)*d*l*Pow(m,3)*Pow(n,2) + 512*Pow(a,3)*c*d*g*l*Pow(m,3)*Pow(n,2) +
               512*Pow(a,3)*Pow(c,2)*h*l*Pow(m,3)*Pow(n,2) - 2048*Pow(a,4)*g*h*l*Pow(m,3)*Pow(n,2) - 16*Pow(a,2)*Pow(c,4)*Pow(m,4)*Pow(n,2) +
               128*Pow(a,3)*Pow(c,2)*g*Pow(m,4)*Pow(n,2) - 256*Pow(a,4)*Pow(g,2)*Pow(m,4)*Pow(n,2) - 4*Pow(c,3)*Pow(d,2)*Pow(e,3)*g*Pow(n,3) +
               8*Pow(c,4)*d*Pow(e,2)*f*g*Pow(n,3) - 4*Pow(c,5)*e*Pow(f,2)*g*Pow(n,3) - 120*b*Pow(c,2)*Pow(d,2)*Pow(e,2)*Pow(g,2)*Pow(n,3) +
               144*a*c*Pow(d,2)*Pow(e,3)*Pow(g,2)*Pow(n,3) + 64*b*Pow(c,3)*d*e*f*Pow(g,2)*Pow(n,3) - 128*a*Pow(c,2)*d*Pow(e,2)*f*Pow(g,2)*Pow(n,3) +
               8*b*Pow(c,4)*Pow(f,2)*Pow(g,2)*Pow(n,3) + 32*a*Pow(c,3)*e*Pow(f,2)*Pow(g,2)*Pow(n,3) + 384*Pow(b,2)*c*Pow(d,2)*e*Pow(g,3)*Pow(n,3) -
               288*a*b*Pow(d,2)*Pow(e,2)*Pow(g,3)*Pow(n,3) - 64*Pow(b,2)*Pow(c,2)*d*f*Pow(g,3)*Pow(n,3) - 256*a*b*c*d*e*f*Pow(g,3)*Pow(n,3) +
               384*Pow(a,2)*d*Pow(e,2)*f*Pow(g,3)*Pow(n,3) - 64*a*b*Pow(c,2)*Pow(f,2)*Pow(g,3)*Pow(n,3) - 64*Pow(a,2)*c*e*Pow(f,2)*Pow(g,3)*Pow(n,3) -
               256*Pow(b,3)*Pow(d,2)*Pow(g,4)*Pow(n,3) + 256*a*Pow(b,2)*d*f*Pow(g,4)*Pow(n,3) + 128*Pow(a,2)*b*Pow(f,2)*Pow(g,4)*Pow(n,3) +
               176*b*Pow(c,3)*d*Pow(e,2)*g*h*Pow(n,3) - 160*a*Pow(c,2)*d*Pow(e,3)*g*h*Pow(n,3) - 80*b*Pow(c,4)*e*f*g*h*Pow(n,3) +
               64*a*Pow(c,3)*Pow(e,2)*f*g*h*Pow(n,3) - 704*Pow(b,2)*Pow(c,2)*d*e*Pow(g,2)*h*Pow(n,3) + 832*a*b*c*d*Pow(e,2)*Pow(g,2)*h*Pow(n,3) -
               384*Pow(a,2)*d*Pow(e,3)*Pow(g,2)*h*Pow(n,3) + 64*Pow(b,2)*Pow(c,3)*f*Pow(g,2)*h*Pow(n,3) + 384*a*b*Pow(c,2)*e*f*Pow(g,2)*h*Pow(n,3) -
               256*Pow(a,2)*c*Pow(e,2)*f*Pow(g,2)*h*Pow(n,3) + 512*Pow(b,3)*c*d*Pow(g,3)*h*Pow(n,3) - 256*a*Pow(b,2)*d*e*Pow(g,3)*h*Pow(n,3) -
               256*a*Pow(b,2)*c*f*Pow(g,3)*h*Pow(n,3) - 256*Pow(a,2)*b*e*f*Pow(g,3)*h*Pow(n,3) + 16*b*Pow(c,4)*Pow(e,2)*Pow(h,2)*Pow(n,3) -
               16*a*Pow(c,3)*Pow(e,3)*Pow(h,2)*Pow(n,3) + 128*Pow(b,2)*Pow(c,3)*e*g*Pow(h,2)*Pow(n,3) - 736*a*b*Pow(c,2)*Pow(e,2)*g*Pow(h,2)*Pow(n,3) +
               576*Pow(a,2)*c*Pow(e,3)*g*Pow(h,2)*Pow(n,3) - 128*Pow(b,3)*Pow(c,2)*Pow(g,2)*Pow(h,2)*Pow(n,3) +
               1024*a*Pow(b,2)*c*e*Pow(g,2)*Pow(h,2)*Pow(n,3) - 384*Pow(a,2)*b*Pow(e,2)*Pow(g,2)*Pow(h,2)*Pow(n,3) -
               512*a*Pow(b,3)*Pow(g,3)*Pow(h,2)*Pow(n,3) - 64*b*Pow(c,4)*Pow(e,2)*g*k*Pow(n,3) + 64*a*Pow(c,3)*Pow(e,3)*g*k*Pow(n,3) +
               192*Pow(b,2)*Pow(c,3)*e*Pow(g,2)*k*Pow(n,3) + 128*a*b*Pow(c,2)*Pow(e,2)*Pow(g,2)*k*Pow(n,3) - 256*Pow(a,2)*c*Pow(e,3)*Pow(g,2)*k*Pow(n,3) -
               128*Pow(b,3)*Pow(c,2)*Pow(g,3)*k*Pow(n,3) - 768*a*Pow(b,2)*c*e*Pow(g,3)*k*Pow(n,3) + 512*Pow(a,2)*b*Pow(e,2)*Pow(g,3)*k*Pow(n,3) +
               512*a*Pow(b,3)*Pow(g,4)*k*Pow(n,3) + 4*Pow(c,4)*Pow(d,2)*Pow(e,2)*l*Pow(n,3) - 8*Pow(c,5)*d*e*f*l*Pow(n,3) + 4*Pow(c,6)*Pow(f,2)*l*Pow(n,3) +
               176*b*Pow(c,3)*Pow(d,2)*e*g*l*Pow(n,3) - 184*a*Pow(c,2)*Pow(d,2)*Pow(e,2)*g*l*Pow(n,3) - 80*b*Pow(c,4)*d*f*g*l*Pow(n,3) +
               128*a*Pow(c,3)*d*e*f*g*l*Pow(n,3) - 40*a*Pow(c,4)*Pow(f,2)*g*l*Pow(n,3) - 320*Pow(b,2)*Pow(c,2)*Pow(d,2)*Pow(g,2)*l*Pow(n,3) +
               64*a*b*c*Pow(d,2)*e*Pow(g,2)*l*Pow(n,3) - 96*Pow(a,2)*Pow(d,2)*Pow(e,2)*Pow(g,2)*l*Pow(n,3) + 512*a*b*Pow(c,2)*d*f*Pow(g,2)*l*Pow(n,3) -
               384*Pow(a,2)*c*d*e*f*Pow(g,2)*l*Pow(n,3) + 128*Pow(a,2)*Pow(c,2)*Pow(f,2)*Pow(g,2)*l*Pow(n,3) + 512*a*Pow(b,2)*Pow(d,2)*Pow(g,3)*l*Pow(n,3) -
               768*Pow(a,2)*b*d*f*Pow(g,3)*l*Pow(n,3) - 128*Pow(a,3)*Pow(f,2)*Pow(g,3)*l*Pow(n,3) - 192*b*Pow(c,4)*d*e*h*l*Pow(n,3) +
               176*a*Pow(c,3)*d*Pow(e,2)*h*l*Pow(n,3) + 96*b*Pow(c,5)*f*h*l*Pow(n,3) - 80*a*Pow(c,4)*e*f*h*l*Pow(n,3) + 576*Pow(b,2)*Pow(c,3)*d*g*h*l*Pow(n,3) -
               1024*a*b*Pow(c,2)*d*e*g*h*l*Pow(n,3) + 832*Pow(a,2)*c*d*Pow(e,2)*g*h*l*Pow(n,3) - 640*a*b*Pow(c,3)*f*g*h*l*Pow(n,3) +
               384*Pow(a,2)*Pow(c,2)*e*f*g*h*l*Pow(n,3) - 768*a*Pow(b,2)*c*d*Pow(g,2)*h*l*Pow(n,3) + 1024*Pow(a,2)*b*d*e*Pow(g,2)*h*l*Pow(n,3) +
               1024*Pow(a,2)*b*c*f*Pow(g,2)*h*l*Pow(n,3) - 256*Pow(a,3)*e*f*Pow(g,2)*h*l*Pow(n,3) - 192*Pow(b,2)*Pow(c,4)*Pow(h,2)*l*Pow(n,3) +
               704*a*b*Pow(c,3)*e*Pow(h,2)*l*Pow(n,3) - 480*Pow(a,2)*Pow(c,2)*Pow(e,2)*Pow(h,2)*l*Pow(n,3) + 128*a*Pow(b,2)*Pow(c,2)*g*Pow(h,2)*l*Pow(n,3) +
               256*Pow(a,2)*b*c*e*g*Pow(h,2)*l*Pow(n,3) - 1152*Pow(a,3)*Pow(e,2)*g*Pow(h,2)*l*Pow(n,3) - 512*Pow(a,2)*Pow(b,2)*Pow(g,2)*Pow(h,2)*l*Pow(n,3) +
               64*b*Pow(c,5)*e*k*l*Pow(n,3) - 64*a*Pow(c,4)*Pow(e,2)*k*l*Pow(n,3) - 128*Pow(b,2)*Pow(c,4)*g*k*l*Pow(n,3) - 128*a*b*Pow(c,3)*e*g*k*l*Pow(n,3) +
               128*Pow(a,2)*Pow(c,2)*Pow(e,2)*g*k*l*Pow(n,3) + 640*a*Pow(b,2)*Pow(c,2)*Pow(g,2)*k*l*Pow(n,3) - 512*Pow(a,2)*b*c*e*Pow(g,2)*k*l*Pow(n,3) +
               512*Pow(a,3)*Pow(e,2)*Pow(g,2)*k*l*Pow(n,3) - 512*Pow(a,2)*Pow(b,2)*Pow(g,3)*k*l*Pow(n,3) - 48*b*Pow(c,4)*Pow(d,2)*Pow(l,2)*Pow(n,3) +
               32*a*Pow(c,3)*Pow(d,2)*e*Pow(l,2)*Pow(n,3) + 16*a*Pow(c,4)*d*f*Pow(l,2)*Pow(n,3) + 32*a*b*Pow(c,2)*Pow(d,2)*g*Pow(l,2)*Pow(n,3) +
               256*Pow(a,2)*c*Pow(d,2)*e*g*Pow(l,2)*Pow(n,3) - 192*Pow(a,2)*Pow(c,2)*d*f*g*Pow(l,2)*Pow(n,3) -
               128*Pow(a,2)*b*Pow(d,2)*Pow(g,2)*Pow(l,2)*Pow(n,3) + 512*Pow(a,3)*d*f*Pow(g,2)*Pow(l,2)*Pow(n,3) + 576*a*b*Pow(c,3)*d*h*Pow(l,2)*Pow(n,3) -
               704*Pow(a,2)*Pow(c,2)*d*e*h*Pow(l,2)*Pow(n,3) + 64*Pow(a,2)*Pow(c,3)*f*h*Pow(l,2)*Pow(n,3) - 768*Pow(a,2)*b*c*d*g*h*Pow(l,2)*Pow(n,3) -
               256*Pow(a,3)*d*e*g*h*Pow(l,2)*Pow(n,3) - 256*Pow(a,3)*c*f*g*h*Pow(l,2)*Pow(n,3) - 1280*Pow(a,2)*b*Pow(c,2)*Pow(h,2)*Pow(l,2)*Pow(n,3) +
               1536*Pow(a,3)*c*e*Pow(h,2)*Pow(l,2)*Pow(n,3) + 2048*Pow(a,3)*b*g*Pow(h,2)*Pow(l,2)*Pow(n,3) - 128*a*b*Pow(c,4)*k*Pow(l,2)*Pow(n,3) +
               192*Pow(a,2)*Pow(c,3)*e*k*Pow(l,2)*Pow(n,3) + 640*Pow(a,2)*b*Pow(c,2)*g*k*Pow(l,2)*Pow(n,3) - 768*Pow(a,3)*c*e*g*k*Pow(l,2)*Pow(n,3) -
               512*Pow(a,3)*b*Pow(g,2)*k*Pow(l,2)*Pow(n,3) - 32*Pow(a,2)*Pow(c,2)*Pow(d,2)*Pow(l,3)*Pow(n,3) - 128*Pow(a,3)*Pow(d,2)*g*Pow(l,3)*Pow(n,3) +
               512*Pow(a,3)*c*d*h*Pow(l,3)*Pow(n,3) - 1024*Pow(a,4)*Pow(h,2)*Pow(l,3)*Pow(n,3) - 128*Pow(a,3)*Pow(c,2)*k*Pow(l,3)*Pow(n,3) +
               512*Pow(a,4)*g*k*Pow(l,3)*Pow(n,3) - 80*b*Pow(c,4)*d*e*g*m*Pow(n,3) + 64*a*Pow(c,3)*d*Pow(e,2)*g*m*Pow(n,3) - 16*b*Pow(c,5)*f*g*m*Pow(n,3) +
               32*a*Pow(c,4)*e*f*g*m*Pow(n,3) + 64*Pow(b,2)*Pow(c,3)*d*Pow(g,2)*m*Pow(n,3) + 384*a*b*Pow(c,2)*d*e*Pow(g,2)*m*Pow(n,3) -
               256*Pow(a,2)*c*d*Pow(e,2)*Pow(g,2)*m*Pow(n,3) + 128*a*b*Pow(c,3)*f*Pow(g,2)*m*Pow(n,3) - 256*Pow(a,2)*Pow(c,2)*e*f*Pow(g,2)*m*Pow(n,3) -
               256*a*Pow(b,2)*c*d*Pow(g,3)*m*Pow(n,3) - 256*Pow(a,2)*b*d*e*Pow(g,3)*m*Pow(n,3) - 256*Pow(a,2)*b*c*f*Pow(g,3)*m*Pow(n,3) +
               512*Pow(a,3)*e*f*Pow(g,3)*m*Pow(n,3) - 32*b*Pow(c,5)*e*h*m*Pow(n,3) + 32*a*Pow(c,4)*Pow(e,2)*h*m*Pow(n,3) + 64*Pow(b,2)*Pow(c,4)*g*h*m*Pow(n,3) +
               512*a*b*Pow(c,3)*e*g*h*m*Pow(n,3) - 512*Pow(a,2)*Pow(c,2)*Pow(e,2)*g*h*m*Pow(n,3) - 768*a*Pow(b,2)*Pow(c,2)*Pow(g,2)*h*m*Pow(n,3) -
               1536*Pow(a,2)*b*c*e*Pow(g,2)*h*m*Pow(n,3) + 1536*Pow(a,3)*Pow(e,2)*Pow(g,2)*h*m*Pow(n,3) + 2048*Pow(a,2)*Pow(b,2)*Pow(g,3)*h*m*Pow(n,3) +
               96*b*Pow(c,5)*d*l*m*Pow(n,3) - 80*a*Pow(c,4)*d*e*l*m*Pow(n,3) - 16*a*Pow(c,5)*f*l*m*Pow(n,3) - 640*a*b*Pow(c,3)*d*g*l*m*Pow(n,3) +
               384*Pow(a,2)*Pow(c,2)*d*e*g*l*m*Pow(n,3) + 128*Pow(a,2)*Pow(c,3)*f*g*l*m*Pow(n,3) + 1024*Pow(a,2)*b*c*d*Pow(g,2)*l*m*Pow(n,3) -
               256*Pow(a,3)*d*e*Pow(g,2)*l*m*Pow(n,3) - 256*Pow(a,3)*c*f*Pow(g,2)*l*m*Pow(n,3) - 320*a*b*Pow(c,4)*h*l*m*Pow(n,3) +
               256*Pow(a,2)*Pow(c,3)*e*h*l*m*Pow(n,3) + 2048*Pow(a,2)*b*Pow(c,2)*g*h*l*m*Pow(n,3) - 1024*Pow(a,3)*c*e*g*h*l*m*Pow(n,3) -
               3072*Pow(a,3)*b*Pow(g,2)*h*l*m*Pow(n,3) + 64*Pow(a,2)*Pow(c,3)*d*Pow(l,2)*m*Pow(n,3) - 256*Pow(a,3)*c*d*g*Pow(l,2)*m*Pow(n,3) -
               256*Pow(a,3)*Pow(c,2)*h*Pow(l,2)*m*Pow(n,3) + 1024*Pow(a,4)*g*h*Pow(l,2)*m*Pow(n,3) + 16*b*Pow(c,6)*Pow(m,2)*Pow(n,3) -
               16*a*Pow(c,5)*e*Pow(m,2)*Pow(n,3) - 160*a*b*Pow(c,4)*g*Pow(m,2)*Pow(n,3) + 128*Pow(a,2)*Pow(c,3)*e*g*Pow(m,2)*Pow(n,3) +
               512*Pow(a,2)*b*Pow(c,2)*Pow(g,2)*Pow(m,2)*Pow(n,3) - 256*Pow(a,3)*c*e*Pow(g,2)*Pow(m,2)*Pow(n,3) - 512*Pow(a,3)*b*Pow(g,3)*Pow(m,2)*Pow(n,3) +
               32*Pow(a,2)*Pow(c,4)*l*Pow(m,2)*Pow(n,3) - 256*Pow(a,3)*Pow(c,2)*g*l*Pow(m,2)*Pow(n,3) + 512*Pow(a,4)*Pow(g,2)*l*Pow(m,2)*Pow(n,3) +
               16*b*Pow(c,5)*e*g*Pow(n,4) - 16*a*Pow(c,4)*Pow(e,2)*g*Pow(n,4) - 16*Pow(b,2)*Pow(c,4)*Pow(g,2)*Pow(n,4) - 128*a*b*Pow(c,3)*e*Pow(g,2)*Pow(n,4) +
               128*Pow(a,2)*Pow(c,2)*Pow(e,2)*Pow(g,2)*Pow(n,4) + 128*a*Pow(b,2)*Pow(c,2)*Pow(g,3)*Pow(n,4) + 256*Pow(a,2)*b*c*e*Pow(g,3)*Pow(n,4) -
               256*Pow(a,3)*Pow(e,2)*Pow(g,3)*Pow(n,4) - 256*Pow(a,2)*Pow(b,2)*Pow(g,4)*Pow(n,4) - 16*b*Pow(c,6)*l*Pow(n,4) + 16*a*Pow(c,5)*e*l*Pow(n,4) +
               160*a*b*Pow(c,4)*g*l*Pow(n,4) - 128*Pow(a,2)*Pow(c,3)*e*g*l*Pow(n,4) - 512*Pow(a,2)*b*Pow(c,2)*Pow(g,2)*l*Pow(n,4) +
               256*Pow(a,3)*c*e*Pow(g,2)*l*Pow(n,4) + 512*Pow(a,3)*b*Pow(g,3)*l*Pow(n,4) - 16*Pow(a,2)*Pow(c,4)*Pow(l,2)*Pow(n,4) +
               128*Pow(a,3)*Pow(c,2)*g*Pow(l,2)*Pow(n,4) - 256*Pow(a,4)*Pow(g,2)*Pow(l,2)*Pow(n,4));
}
#ifndef GAMMA_H_
#define GAMMA_H_

#include <math.h>
#include <limits.h>
#include <cfloat>

int DiscreteGamma(double* freqK,double* rK,double alfa,double beta,int K,bool median);
double LnGamma(double alpha);
double IncompleteGamma(double x,double alpha,double ln_gamma_alpha);
double PointNormal(double p);
double PointChi2(double prob,double v);

#define PointGamma(prob,alpha,beta)  PointChi2(prob,2.0*(alpha))/(2.0*(beta))

//=======================================================================================================
//DP-PIP
double PointChi2(double prob,double v){
	/*********************************************************/
	/* Inverse CDFs */
	/*********************************************************/
	/* returns z so that Prob{x<z}=prob where x is Chi2 distributed with df=v
   returns -1 if in error.   0.000002<prob<0.999998
   RATNEST FORTRAN by
       Best DJ & Roberts DE (1975) The percentage points of the
       Chi2 distribution.  Applied Statistics 24: 385-388.  (AS91)
   Converted into C by Ziheng Yang, Oct. 1993.
*/
   double aa=.6931471805, p=prob, g;
   double xx, c, ch, a=0,q=0,p1=0,p2=0,t=0,x=0,b=0,s1,s2,s3,s4,s5,s6;
   double e=.5e-6;

   if (p<.000002 || p>.999998 || v<=0) return ((double)-1);

   g = (double)LnGamma(v/2);
   xx=v/2;   c=xx-1;
   if (v >= -1.24*log(p)) goto l1;

   ch=pow((p*xx*exp(g+xx*aa)), 1/xx);
   if (ch-e<0) return (ch);
   goto l4;
l1:
   if (v>.32) goto l3;
   ch=0.4;   a=log(1-p);
l2:
   q=ch;  p1=1+ch*(4.67+ch);  p2=ch*(6.73+ch*(6.66+ch));
   t=-0.5+(4.67+2*ch)/p1 - (6.73+ch*(13.32+3*ch))/p2;
   ch-=(1-exp(a+g+.5*ch+c*aa)*p2/p1)/t;
   if (fabs(q/ch-1)-.01 <= 0) goto l4;
   else                       goto l2;

l3:
   x=(double)PointNormal (p);
   p1=0.222222/v;   ch=v*pow((x*sqrt(p1)+1-p1), 3.0);
   if (ch>2.2*v+6)  ch=-2*(log(1-p)-c*log(.5*ch)+g);
l4:
   q=ch;   p1=.5*ch;
   if ((t=(double)IncompleteGamma (p1, xx, g))<0) {
      printf("\nerr IncompleteGamma");
      return ((double)-1.);
   }
   p2=p-t;
   t=p2*exp(xx*aa+g+p1-c*log(ch));
   b=t/ch;  a=0.5*t-b*c;

   s1=(210+a*(140+a*(105+a*(84+a*(70+60*a))))) / 420;
   s2=(420+a*(735+a*(966+a*(1141+1278*a))))/2520;
   s3=(210+a*(462+a*(707+932*a)))/2520;
   s4=(252+a*(672+1182*a)+c*(294+a*(889+1740*a)))/5040;
   s5=(84+264*a+c*(175+606*a))/2520;
   s6=(120+c*(346+127*c))/5040;
   ch+=t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
   if (fabs(q/ch-1) > e) goto l4;

   return (double)(ch);
}
//=======================================================================================================
//DP-PIP
double PointNormal(double p){
/* returns z so that Prob{x<z}=prob where x ~ N(0,1) and (1e-12)<prob<1-(1e-12)
   returns (-9999) if in error
   Odeh RE & Evans JO (1974) The percentage points of the normal distribution.
   Applied Statistics 22: 96-97 (AS70)

   Newer methods:
     Wichura MJ (1988) Algorithm AS 241: the percentage points of the
       normal distribution.  37: 477-484.
     Beasley JD & Springer SG  (1977).  Algorithm AS 111: the percentage
       points of the normal distribution.  26: 118-121.

*/
/*    double a0=-.322232431088, a1=-1, a2=-.342242088547, a3=-.0204231210245; */
/*    double a4=-.453642210148e-4, b0=.0993484626060, b1=.588581570495; */
/*    double b2=.531103462366, b3=.103537752850, b4=.0038560700634; */
/*    double y, z=0, p=prob, p1; */

/*    p1 = (p<0.5 ? p : 1-p); */
/*    if (p1<1e-20) return (-INFINITY); */
/* /\*    if (p1<1e-20) return (-999.); *\/ */

/*    y = sqrt ((double)LOG(1/(p1*p1))); */
/*    z = y + ((((y*a4+a3)*y+a2)*y+a1)*y+a0) / ((((y*b4+b3)*y+b2)*y+b1)*y+b0); */
/*    return (p<0.5 ? -z : z); */

  static double zero = 0.0, one = 1.0, half = 0.5;
  static double split1 = 0.425, split2 = 5.0;
  static double const1 = 0.180625, const2 = 1.6;

  /* coefficients for p close to 0.5 */
  static double a[8] = {
    3.3871328727963666080e0,
    1.3314166789178437745e+2,
    1.9715909503065514427e+3,
    1.3731693765509461125e+4,
    4.5921953931549871457e+4,
    6.7265770927008700853e+4,
    3.3430575583588128105e+4,
    2.5090809287301226727e+3
  };
  static double b[8] = {
    0.0,
    4.2313330701600911252e+1,
    6.8718700749205790830e+2,
    5.3941960214247511077e+3,
    2.1213794301586595867e+4,
    3.9307895800092710610e+4,
    2.8729085735721942674e+4,
    5.2264952788528545610e+3
  };

  /* hash sum ab    55.8831928806149014439 */
  /* coefficients for p not close to 0, 0.5 or 1. */
  static double c[8] = {
    1.42343711074968357734e0,
    4.63033784615654529590e0,
    5.76949722146069140550e0,
    3.64784832476320460504e0,
    1.27045825245236838258e0,
    2.41780725177450611770e-1,
    2.27238449892691845833e-2,
    7.74545014278341407640e-4
  };
  static double d[8] = {
    0.0,
    2.05319162663775882187e0,
    1.67638483018380384940e0,
    6.89767334985100004550e-1,
    1.48103976427480074590e-1,
    1.51986665636164571966e-2,
    5.47593808499534494600e-4,
    1.05075007164441684324e-9
  };

  /* hash sum cd    49.33206503301610289036 */
  /* coefficients for p near 0 or 1. */
  static double e[8] = {
    6.65790464350110377720e0,
    5.46378491116411436990e0,
    1.78482653991729133580e0,
    2.96560571828504891230e-1,
    2.65321895265761230930e-2,
    1.24266094738807843860e-3,
    2.71155556874348757815e-5,
    2.01033439929228813265e-7
  };
  static double f[8] = {
    0.0,
    5.99832206555887937690e-1,
    1.36929880922735805310e-1,
    1.48753612908506148525e-2,
    7.86869131145613259100e-4,
    1.84631831751005468180e-5,
    1.42151175831644588870e-7,
    2.04426310338993978564e-15
  };

  /* hash sum ef    47.52583317549289671629 */
  double q, r, ret;

  q = p - half;
  if (fabs(q) <= split1) {
    r = const1 - q * q;
    ret = q * (((((((a[7] * r + a[6]) * r + a[5]) * r + a[4]) * r + a[3])
		 * r + a[2]) * r + a[1]) * r + a[0]) /
      (((((((b[7] * r + b[6]) * r + b[5]) * r + b[4]) * r + b[3])
	 * r + b[2]) * r + b[1]) * r + one);

    return (double)ret;
  }
  /* else */

  if (q < zero)
    r = p;
  else
    r = one - p;

  if (r <= zero)
    return (double)zero;

  r = sqrt(-log(r));
  if (r <= split2) {
    r -= const2;
    ret = (((((((c[7] * r + c[6]) * r + c[5]) * r + c[4]) * r + c[3])
	     * r + c[2]) * r + c[1]) * r + c[0]) /
      (((((((d[7] * r + d[6]) * r + d[5]) * r + d[4]) * r + d[3])
	 * r + d[2]) * r + d[1]) * r + one);
  }
  else {
    r -= split2;
    ret = (((((((e[7] * r + e[6]) * r + e[5]) * r + e[4]) * r + e[3])
	     * r + e[2]) * r + e[1]) * r + e[0]) /
      (((((((f[7] * r + f[6]) * r + f[5]) * r + f[4]) * r + f[3])
	 * r + f[2]) * r + f[1]) * r + one);
  }

  if (q < zero)
    ret = -ret;

  return (double)ret;
}
//=======================================================================================================
//DP-PIP
double IncompleteGamma(double x,double alpha,double ln_gamma_alpha){
/* returns the incomplete gamma ratio I(x,alpha) where x is the upper
	   limit of the integration and alpha is the shape parameter.
   returns (-1) if in error
   ln_gamma_alpha = ln(Gamma(alpha)), is almost redundant.
   (1) series expansion     if (alpha>x || x<=1)
   (2) continued fraction   otherwise
   RATNEST FORTRAN by
   Bhattacharjee GP (1970) The incomplete gamma integral.  Applied Statistics,
   19: 285-287 (AS32)
*/
   int i;
   double p=alpha, g=ln_gamma_alpha;
   double accurate=1e-8, overflow=1e30;
   double factor, gin=0, rn=0, a=0,b=0,an=0,dif=0, term=0, pn[6];

   if (fabs(x) < DBL_MIN) return ((double).0);
   if (x<0 || p<=0)        return ((double)-1);

   factor=exp(p*log(x)-x-g);
   if (x>1 && x>=p) goto l30;
   /* (1) series expansion */
   gin=1;  term=1;  rn=p;
 l20:
   rn++;
   term*=x/rn;   gin+=term;

   if (term > accurate) goto l20;
   gin*=factor/p;
   goto l50;
 l30:
   /* (2) continued fraction */
   a=1-p;   b=a+x+1;  term=0;
   pn[0]=1;  pn[1]=x;  pn[2]=x+1;  pn[3]=x*b;
   gin=pn[2]/pn[3];
 l32:
   a++;  b+=2;  term++;   an=a*term;
   for (i=0; i<2; i++) pn[i+4]=b*pn[i+2]-an*pn[i];
   if (fabs(pn[5]) < .0) goto l35;
   rn=pn[4]/pn[5];   dif=fabs(gin-rn);
   if (dif>accurate) goto l34;
   if (dif<=accurate*rn) goto l42;
 l34:
   gin=rn;
 l35:
   for (i=0; i<4; i++) pn[i]=pn[i+2];
   if (fabs(pn[4]) < overflow) goto l32;
   for (i=0; i<4; i++) pn[i]/=overflow;
   goto l32;
 l42:
   gin=1-factor*gin;

 l50:
   return (double)(gin);
}
//=======================================================================================================
//DP-PIP
double LnGamma(double alpha){
/* returns ln(gamma(alpha)) for alpha>0, accurate to 10 decimal places.
   Stirling's formula is used for the central polynomial part of the procedure.
   Pike MC & Hill ID (1966) Algorithm 291: Logarithm of the gamma function.
   Communications of the Association for Computing Machinery, 9:684
*/
   double x=alpha, f=0, z;
   if (x<7) {
      f=1;  z=x-1;
      while (++z<7)  f*=z;
      x=z;   f=-log(f);
   }
   z = 1/(x*x);
   return (double)(f + (x-0.5)*log(x) - x + .918938533204673
		   + (((-.000595238095238*z+.000793650793651)*z-.002777777777778)*z
		      +.083333333333333)/x);
}
//=======================================================================================================
//DP-PIP
int DiscreteGamma(double* freqK,double* rK,double alfa,double beta,int K,bool median){

	/* discretization of gamma distribution with equal proportions in each category */

	int i;
	double gap05=1.0/(2.0*(double)K);
	double t;
	double factor=alfa/beta*(double)K;
	double lnga1;

	if(K==1){
		freqK[0] = 1.0;
		rK[0] = 1.0;
		return 0;
	}

	if(median){

		for (i=0; i<K; i++){
			rK[i]=PointGamma((i*2.0+1)*gap05, alfa, beta);
		}
		for (i=0,t=0; i<K; i++){
			t+=rK[i];
		}
		for (i=0; i<K; i++){
			rK[i]*=factor/t;
		}

	}else{
		lnga1=LnGamma(alfa+1);
		for (i=0; i<K-1; i++){
			freqK[i]=PointGamma((i+1.0)/K, alfa, beta);
		}
		for (i=0; i<K-1; i++){
			freqK[i]=IncompleteGamma(freqK[i]*beta, alfa+1, lnga1);
		}
		rK[0] = freqK[0]*factor;
		rK[K-1] = (1-freqK[K-2])*factor;
		for (i=1; i<K-1; i++){
			rK[i] = (freqK[i]-freqK[i-1])*factor;
		}
	}

	for (i=0; i<K; i++){
		freqK[i]=1.0/K;
	}

	return (0);
}
//=======================================================================================================

#endif /* GAMMA_H_ */

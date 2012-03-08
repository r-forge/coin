#include <stdlib.h>
#include <stdio.h>
#include <iostream.h>
#include <math.h>
double *score;
 
void mult(double *p1, long lengtep1, double *p2, long &lengtep2, long xi)
{
    
    for (int j=0 ; j < lengtep1*2 ; j=j+2) 
    {
            
        p2[j] = p1[j] ;
        p2[j+1] = p1[j+1] + score[xi];  
    }
    lengtep2= lengtep1;
}

void mergesort(double *p, long lengtep, long templen, long be)
{
    double *copiep = new(double[be]);
    long t1=0; 
    long t2=0;
        
    for (long i=0 ; i< lengtep*2; i=i+1) copiep[i]=p[i];
    
    for (long j=0 ; j< lengtep*2; j=j+2)
    { 
    if (t1<=2*templen-2 && t2<=2*(lengtep-templen)-2)
         {
        if (copiep[t1+1] < copiep[2*templen+t2+1])
        {
            p[j]=copiep[t1];
            p[j+1]=copiep[t1+1];
            t1=t1+2;
        }
        else
        {
            p[j]=copiep[2*templen+t2];
            p[j+1]=copiep[2*templen+1+t2];
            t2=t2+2;
        }
         }
     else
         {
            if (t1>2*templen-2)
            {
                p[j]=copiep[2*templen+t2];
            p[j+1]=copiep[2*templen+1+t2];
            t2=t2+2; 
            }
            else
            {   
                p[j]=copiep[t1];
            p[j+1]=copiep[t1+1];
            t1=t1+2;
            }
         }          
        } 
}

// void compute()

void plus(double *p, long &lengtep, double *p2, long lengtep2)
{
    long elep=0;
        long k=0;
        long test=1;
    
    for (long i=0 ; i<lengtep2*2; i=i+2)
    {
    test=1;
        for (long j=elep ; j<lengtep*2 && test==1; j=j+2)
        {
            if (p2[i+1]-0.000001<=p[j+1] && p[j+1]<=p2[i+1]+0.000001) 
            {
                p[j] += p2[i];
                test=0;
                elep=i;             
            }
        } 
        if (test==1) 
            {
                p[2*lengtep + k] = p2[i];
                    p[2*lengtep + k+1]= p2[i+1];
                    k=k+2;
            }
        }
        lengtep=lengtep+k/2;
}

void printP(double *p,long lengte)
{
    for(long i=0 ; i< lengte*2; i=i+2)
    {
        cout <<p[i]<<"*X["<<p[i+1]<<"]";
        if (i<lengte*2-2) cout << "+ ";
        else cout << "\n";
    }   
    cout << endl;

}

void vulpoly(double *p, long &lengtep, int a, int b, long help3)
{       
    long h3=help3; 
        double *p3=new(double[h3]); 
        long lengtep3=0;
    long templen=0;
    
    for(int i=a; i<=b;i=i+1)
    {
            
        mult(p,lengtep,p3,lengtep3,i);          
        templen=lengtep; 
        plus(p,lengtep,p3,lengtep3);        
        mergesort(p,lengtep,templen,h3);        
//      printP(p,lengtep);  
    }
    
}

double aantalklein(double *p1, long &lengtep1, double *p2, long &lengtep2, double ob)
{
    double tot = 0;
    long tempel=0;
    long test=1;
    
    for (long i=0 ; i<2*lengtep1 ; i=i+2)
    {
        test=1;
        for (long j=tempel ; j<2*lengtep2 && test==1 ; j=j+2)
        {
            if (p1[i+1]+p2[2*lengtep2-1-j] <= ob)
            {
            tot = tot + p1[i]*p2[2*lengtep2-j-2];
            tempel=j;
            test=0;
            }
        }
    }
    return tot;
}

void cumulcoef(double *p, long &lengtep)
{
     double coef=0;
     for(long i=0 ; i<2*lengtep ; i=i+2)
        {
        p[i]=p[i]+coef;
        coef=p[i];
        }
}

double PHI_inverse(double ALPH)
{
     // Inverse cdf of a N(0,1)
     // PHI(x) := P(X<=x), for X ~ N(0,1);
     // then return PHI^(-1)(x);

     if (ALPH==0.5) return 0;
     else
     {
    const double C0=2.515517;
    const double C1=0.802853;
    const double C2=0.010328;
    const double D1=1.432788;
    const double D2=0.189269;
    const double D3=0.001308;
    double T,Z,zz,zz1;
        T= sqrt(log(1./(ALPH*ALPH)));
        Z = C0    + C1*T + C2*T*T;
        Z = Z/(1. + D1*T + D2*T*T + D3*T*T*T);
        zz =-T+Z;
        zz1= floor(1000*zz+0.5)/1000;       
    return(zz1);
     }
}

void main()
{
     double *poly1;
     double *poly2; 
     int n;
     double obs;
     long help1;
     long help2;
    
     double tot;
     double totprob=1;
     double prob=1; 

     cout << "Give sample size   : " ; cin >> n;
     cout << "Give observed value : " ; cin >> obs;
     score=new(double[n+1]); 
for (int s=1 ; s<=(n+1)/2; s++)
        {
        score[2*s-2]=0.5+(PHI_inverse(double(s)/(n+1)))/2;
        score[2*s-1] = 1-score[2*s-2];
        }
    for (int s1=0 ; s1<=n-1; s1++) cout << score[s1] <<" ";
    
   
    if (n<=30)
        {
    help1= long(pow(2,(n+1)/2+1));    // 
        help2=long(pow(2,n/2+1));
        }
    else
        {
        help1= 300000;    // 
        help2= 300000;
        }
    
     poly1 = new(double[help1]);         
     poly2 = new(double[help2]); 
             

     poly1[0]=1; poly1[1]=0;
     poly2[0]=1; poly2[1]=0;    
     long lengtep1 = 1;
     long lengtep2 = 1;
     
     vulpoly(poly2,lengtep2,(n-1)/2+1,n-1,help2);
     vulpoly(poly1,lengtep1,0,(n-1)/2,help1);
     
//     printP(poly2,lengtep2);
     
    
     cumulcoef(poly2,lengtep2);
 //    printP(poly2,lengtep2);
     tot=aantalklein(poly1,lengtep1,poly2,lengtep2,obs);
    

     totprob=pow(2,n);
     prob=tot/totprob;
     cout<<"Pr(SRN <="<< obs <<")= "<<prob<<"\n";
          
} 
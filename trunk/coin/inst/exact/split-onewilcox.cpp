#include <stdlib.h>
#include <stdio.h>
#include <iostream.h>

double *rang;
 
void mult(double *p1, int lengtep1, double *p2, int &lengtep2, int xi)
{
    
    for (int j=0 ; j < lengtep1*2 ; j=j+2) 
    {
            
        p2[j] = p1[j] ;
        p2[j+1] = p1[j+1] + rang[xi];   
    }
    lengtep2= lengtep1;
}

void mergesort(double *p, int lengtep, int templen, int be)
{
    double *copiep = new(double[be]);
    int t1=0; 
    int t2=0;
        
    for (int i=0 ; i< lengtep*2; i=i+1) copiep[i]=p[i];
    
    for (int j=0 ; j< lengtep*2; j=j+2)
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

void plus(double *p, int &lengtep, double *p2, int lengtep2)
{
    int elep=0;
        int k=0;
        int test=1;
    
    for (int i=0 ; i<lengtep2*2; i=i+2)
    {
    test=1;
        for (int j=elep ; j<lengtep*2 && test==1; j=j+2)
        {
            if (p[j+1]==p2[i+1]) 
            {
                p[j] += p2[i];
                test=0;
                elep=j;  // of i?           
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

void printP(double *p,int lengte)
{
    for(int i=0 ; i< lengte*2; i=i+2)
    {
        cout <<p[i]<<"*X["<<p[i+1]<<"]";
        if (i<lengte*2-2) cout << "+ ";
        else cout << "\n";
    }   
    cout << endl;

}

void vulpoly(double *p, int &lengtep, int a, int b)
{
       
        double *p3=new(double[b*(b+1)+2]); 
        int lengtep3=0;
    int templen=0;
    
    for(int i=a; i<=b;i=i+1)
    {
//          cout<<lengtep<< "\n";
        mult(p,lengtep,p3,lengtep3,i);      
        templen=lengtep; 
        plus(p,lengtep,p3,lengtep3);
        mergesort(p,lengtep,templen,b*(b+1)+2);
    }
//  cout<<lengtep<< "\n";
}

double aantalklein(double *p1, int &lengtep1, double *p2, int &lengtep2, double ob)
{
    double tot = 0;
    int tempel=0;
    int test=1;
    
    for (int i=0 ; i<2*lengtep1 ; i=i+2)
    {
        test=1;
        for (int j=tempel ; j<2*lengtep2 && test==1 ; j=j+2)
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

void cumulcoef(double *p, int &lengtep)
{
     double coef=0;
     for(int i=0 ; i<2*lengtep ; i=i+2)
        {
        p[i]=p[i]+coef;
        coef=p[i];
        }
}

void main()
{
     double *poly1;
     double *poly2; 
     int n;
     int obs;
     int help1;
     int help2;
    
     double tot;
     double totprob=1;
     double prob=1;
       


     cout << "Give sample size   : " ; cin >> n;
     cout << "Give observed value : " ; cin >> obs;
 
    help1= n/2*(n/2+1)+2;
    help2= n*(n+1)+2;
//     help1=5000;
//     help2=5000;
     poly1 = new(double[help1]);         
     poly2 = new(double[help2]);         

     rang = new(double[n+1]);
     for (int i=0 ; i<=n ; i++) rang[i]=i;
     
     poly1[0]=1; poly1[1]=0;
     poly2[0]=1; poly2[1]=0;    
     int lengtep1 = 1;
     int lengtep2 = 1;
     
     vulpoly(poly2,lengtep2,n/2+1,n);
     vulpoly(poly1,lengtep1,1,n/2);
    
//     printP(poly2,lengtep2);
     
    
     cumulcoef(poly2,lengtep2);
     tot=aantalklein(poly1,lengtep1,poly2,lengtep2,obs);
    
     for (int i=1 ; i<=n ; i++) totprob=totprob*2;
     prob=tot/totprob;
     cout<<tot<<" "<<prob<<" "<<totprob<< "\n";
          
} 
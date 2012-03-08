#include <iostream.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

struct celW
{
    long length;
    double *c;
    double *x;
};

celW **W1;
celW **W2;
double *rs;

double binomi(int m, int n)
{ 
        double bin=1;
        double bin1=1;
        double bin2=1;
        
        for(int i=1; i<=n;i++) bin1 = bin1*(m+1-i);
        for(int j=1; j<=n;j++) bin2= bin2*j;
        return bin = bin1/bin2;
}

void reserveerW1(int aa, int bb)
{
    long res = 0;
//  for (int h=1; h<=bb;h++) res +=aa+2-h;
//      res=100;
//      cout << res << "\n";
    W1 = new(celW*[aa+1]);

    for (int i=0 ; i<=aa ; i++)
        W1[i] = new(celW[bb+1]);
    for (int i=0 ; i<=aa ; i++)
        {
        for (int j=i ; j<=bb ; j++)
            {
            res=long(binomi(j,i));
            W1[i][j].c=new(double[res]);
            W1[i][j].x=new(double[res]);
//          cout << i <<" "<< j <<" "<< res << "\n";
            }
        }
}

void initW1(int aa,int bb)
{
    for (int i=1 ; i<=aa ; i++)
    for (int j=0 ; j<=bb ; j++)
        {
        W1[i][j].length =0;     
        }
    for (int j=0 ; j<=bb ; j++)
        {
        W1[0][j].length=1;
        W1[0][j].c[0]=1;
        W1[0][j].x[0]=0;
        }
    
}

void reserveerW2(int aa, int bb)
{
    long res = 0;
//  for (int h=1; h<=bb;h++) res +=aa+2-h;
//      res=100;
//      cout << res << "\n";
    W2 = new(celW*[aa+1]);

    for (int i=0 ; i<=aa ; i++)
        W2[i] = new(celW[bb+1]);
    for (int i=0 ; i<=aa ; i++)
        {
        for (int j=i ; j<=bb ; j++)
            {
            res=long(binomi(j,i));
            W2[i][j].c=new(double[res]);
            W2[i][j].x=new(double[res]);
//          cout << i <<" "<< j <<" "<< res << "\n";
            }
        }
}

void initW2(int aa,int bb)
{
    for (int i=1 ; i<=aa ; i++)
    for (int j=0 ; j<=bb ; j++)
        {
        W2[i][j].length =0;
        
        }
    for (int j=0 ; j<=bb ; j++)
        {
        W2[0][j].length=1;
        W2[0][j].c[0]=1;
        W2[0][j].x[0]=0;
        }
    
}

void printW(celW **W, int aa, int bb)
{
    if (W[aa][bb].length >0)
    {
                
        for (long j=0 ; j<W[aa][bb].length ; j++)
        {       
            cout << W[aa][bb].c[j]<<"*X["<<W[aa][bb].x[j]<<"] ";

            if (j < W[aa][bb].length-1) cout << "+ ";else cout << "\n";
        }
    }
    else cout << "Cel leeg..\n";
}

void mult(celW &tem, int a, int b, int rank)
{
        for (int j=0 ; j < tem.length ; j++)
        tem.x[j] = tem.x[j] + rs[rank];

}

void plus(celW **W, celW &tempie, int a, int b)
{
    int elep=0;
        int k=0;
        int test=1;
    
    for (int i=0; i<W[a][b-1].length; i=i+1)
    {
    test=1;
        for (int j=elep; j<tempie.length && test==1; j=j+1)
        {
            if (tempie.x[j]-0.000001<=W[a][b-1].x[i]
                            && W[a][b-1].x[i]<=tempie.x[j]+0.000001) 
            
            {
                tempie.c[j] += W[a][b-1].c[i];
                test=0;
                elep=j;             
            }
        } 
        if (test==1) 
            {
                tempie.c[tempie.length + k] = W[a][b-1].c[i];
                    tempie.x[tempie.length + k] = W[a][b-1].x[i];
                    k=k+1;
            }
        }
        tempie.length += k;

}

void mergesort(celW temptw, long tijd)
{
    celW copiep;
    copiep.c = new(double[temptw.length]);
    copiep.x = new(double[temptw.length]);
    int t1=0; 
    int t2=0;
        
    for (int i=0 ; i<temptw.length ; i=i+1) 
        {
        copiep.c[i]=temptw.c[i];
        copiep.x[i]=temptw.x[i];
        }
    
    for (int j=0 ; j<temptw.length; j=j+1)
    { 
    if (t1<=tijd-1 && t2<=temptw.length-tijd-1)
         {
        if (copiep.x[t1] < copiep.x[tijd+t2])
        {
            temptw.x[j]=copiep.x[t1];
            temptw.c[j]=copiep.c[t1];
            t1=t1+1;
        }
        else
        {
            temptw.x[j]=copiep.x[tijd+t2];
            temptw.c[j]=copiep.c[tijd+t2];
            t2=t2+1;
        }
         }
     else
         {
            if (t1>tijd-1)
            {
                temptw.x[j]=copiep.x[tijd+t2];
            temptw.c[j]=copiep.c[tijd+t2];
            t2=t2+1; 
            }
            else
            {   
                temptw.x[j]=copiep.x[t1];
            temptw.c[j]=copiep.c[t1];
            t1=t1+1;  
            }
         }          
        } 
}

void vulcel(celW **W, int i1, int j1, int r)
{
    // vul W[i][j];
    // gebruik W[i-1][j-1]*X^r ++ W[i1][j-1];
    
        long tijd;
    celW temp2;

    temp2.c = new(double[W[i1-1][j1-1].length+W[i1][j1-1].length]);
    temp2.x = new(double[W[i1-1][j1-1].length+W[i1][j1-1].length]);
    temp2.length = W[i1-1][j1-1].length;
        for(int j=0;j<temp2.length;j++)
            {
            temp2.c[j]=W[i1-1][j1-1].c[j];
            temp2.x[j]=W[i1-1][j1-1].x[j];
            }
    if (i1==j1)
        {       
            mult(temp2,i1-1,j1-1,r);            
            }
        else
            {           
            mult(temp2,i1-1,j1-1,r);                        
                tijd = temp2.length;                                
                plus(W,temp2,i1,j1);                            
                mergesort(temp2,tijd);                              
            }
    
    W[i1][j1].length = temp2.length;

    for (int j2=0 ; j2<temp2.length ; j2++)
    {
        W[i1][j1].c[j2] = temp2.c[j2];
        W[i1][j1].x[j2] = temp2.x[j2];
    }          
    delete (temp2.c);
    delete (temp2.x);

}

void spiegelW(celW **W,int ce,int bep, int start)     
{   
    double totsum=0;
    for (int r=0 ; r<bep ; r++) totsum += rs[start+r];
    
        long len=W[bep-ce][bep].length;
        
        for (int h=0 ; h< len ;h++)
            {
            W[ce][bep].length=W[bep-ce][bep].length;
            W[ce][bep].c[len-1-h]=W[bep-ce][bep].c[h];
        W[ce][bep].x[len-1-h]=totsum-W[bep-ce][bep].x[h];
//      cout << ce <<" " << bep <<" "<< W[ce][bep].x[len-1-h] <<"\n";
        }
}

void maakW(celW **W,int aa,int bb, int start)
{
    long i,j;
    int rank;
        int hulp;
    for(j=1 ; j<=bb; j++)  //verander naar 0!!
    {
        if (j<aa) hulp =j; 
        else hulp =aa;
        for (i=1 ; i<=hulp; i++)
        {   
            if (i<=j/2 || j==1)
            {   
                rank = start+j;
            vulcel(W,i,j,rank-1);
            }
            else
            {
            spiegelW(W,i,j,start);
            }                               
        }
    }
}

void cumulcoef(celW **W, int i1,int j1)
{
     double coef=0;
     for(int i=0 ; i<W[i1][j1].length ; i=i+1)
        {
        W[i1][j1].c[i]+=coef;
        coef=W[i1][j1].c[i];
        }
}

double aantalklein(int c, int b, double ob)
{
    double tot = 0;
    
    int test=1;
    int be=b/2;
    int bp=(b+1)/2;
    
    for (int h=0 ; h<=c ; h++)
        {
        int tempel=0;
        long le=W2[c-h][bp].length;
        for (int i=0 ; i<W1[h][be].length ; i=i+1)
            {
            test=1;
            for (int j=tempel ; j<le && test==1 ; j=j+1)
                {
                if (W1[h][be].x[i]+W2[c-h][bp].x[le-j-1]<=ob)
                {               
                tot = tot + W1[h][be].c[i]*W2[c-h][bp].c[le-j-1];
                tempel=j;
                test=0;
                }
                }
            }
        }
    return tot;
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
    int a,b,c,d;
    double tot,bino,prob;
    double ob;  
    cout << "Give smallest sample size : "; cin >> c;
    cout << "Give largest sample size: "; cin >> d;
    cout << "Give observed value : "; cin >> ob;
    b=c+d;
    rs = new(double[b+1]);
    for (int s=1 ; s<=(b+1)/2; s++)
        {
        rs[2*s-2]=PHI_inverse(double(s)/(b+1));
        rs[2*s-1] = -rs[2*s-2];
        }
//  for (int s1=0 ; s1<=b-1; s1++) cout << rs[s1] <<" ";
    cout << "\n";
        reserveerW1(c,(b+1)/2);
    initW1(c,(b+1)/2);
        reserveerW2(c,(b+1)/2);
    initW2(c,(b+1)/2);
//  printW(W1,0,0);
    maakW(W1,c,b/2,0);  
//        printW(W1,c,b/2); 
//        printW(W2,0,0);       
        maakW(W2,c,(b+1)/2,b/2);
//        cout << W[(a+1)/4(a+1)/2,(b+1)/2][(b+1)/2].length << "\n";
//  printW(W2,c,(b+1)/2);
    for (int u=0 ; u <= c ; u++) cumulcoef(W2,u,(b+1)/2);
    tot=aantalklein(c,b,ob);
    bino=binomi(b,c);
    prob= tot/bino; 
    cout <<"Pr(W <="<< ob <<")="<< prob << "\n";    
//  printW(W2,c,(b+1)/2);
}



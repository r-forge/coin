#include <iostream.h>
#include <stdlib.h>
#include <stdio.h>

struct celW
{
    long length;
    double *c;
    int *x;
};

celW **W1;
celW **W2;

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
            res=i*(j-i)+2;
            W1[i][j].c=new(double[res]);
            W1[i][j].x=new(int[res]);
            
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
            res=i*(j-i)+2;
            W2[i][j].c=new(double[res]);
            W2[i][j].x=new(int[res]);
            
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
        tem.x[j] = tem.x[j] + rank;

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
            if (W[a][b-1].x[i]==tempie.x[j]) 
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
    copiep.x = new(int[temptw.length]);
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
    temp2.x = new(int[W[i1-1][j1-1].length+W[i1][j1-1].length]);
    temp2.length = W[i1-1][j1-1].length;
        for(int j=0;j<=temp2.length;j++)
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
//  W[i1][j1].c = new(double[temp2.length]);
//  W[i1][j1].x = new(int[temp2.length]);

    for (long j=0 ; j<temp2.length ; j++)
    {
        W[i1][j1].c[j] = temp2.c[j];
        W[i1][j1].x[j] = temp2.x[j];
    }   

    delete (temp2.c);
    delete (temp2.x);

}

void spiegelW(celW **W,int ce,int bep, int start)     
{   
    int totsum=0;
    for (int r=1 ; r<=bep ; r++) totsum += start+r;
    
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
    int i,j;
    int rank;
        int hulp;
    for(j=1 ; j<=bb; j++)  //verander naar 0!!
    {
        if (j<aa) hulp =j; 
        else hulp =aa;
        for (i=1 ; i<=hulp; i++ )
        {   
            if (i<=j/2 || j==1)
            {   
                rank = start+j;
            vulcel(W,i,j,rank);
            }
            else
            {
            spiegelW(W,i,j,start);
            }
//      cout << i << " " << j <<" ";
//          printW(W,i,j);
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
                if (W1[h][be].x[i]+W2[c-h][bp].x[le-j-1] <=ob)
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

double binom(int m, int n)
{ 
        double bin=1;
        double bin1=1;
        double bin2=1;
        
        for(int i=1; i<=n;i++) bin1 = bin1*(m+1-i);
        for(int j=1; j<=n;j++) bin2= bin2*j;
        return bin = bin1/bin2;
}       

 
  
void main()
{
        
    int a,b,c,d,be,bp;
    double tot,bino,prob;
    double ob;
    cout << "Give smallest sample size : "; cin >> c;
    cout << "Give largest sample size: "; cin >> d;
    cout << "Give observed value : "; cin >> ob;
    b=c+d;
    be=b/2;
    bp=(b+1)/2;
        reserveerW1(c,bp);
    initW1(c,bp);
        reserveerW2(c,bp);
    initW2(c,bp);
//  printW(W1,0,0);
        maakW(W1,c,be,0);           
//        printW(W1,c,be); 
//        printW(W2,0,0);       
        maakW(W2,c,bp,be);
        
//        cout << W[(a+1)/4(a+1)/2,(b+1)/2][(b+1)/2].length << "\n";
//  printW(W2,c,bp);
    for (int u=0 ; u <= c ; u++) cumulcoef(W2,u,bp);
    tot=aantalklein(c,b,ob);
    bino=binom(b,c);
    prob= tot/bino;
    cout <<"Pr(W <="<< ob <<")="<< prob << "\n";
//  printW(W2,c,bp);
}



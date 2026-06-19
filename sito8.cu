/** 
@file=sito8.cu
*/
 
#include <omp.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
 
#define BUFSIZE 2048
#define BUFSIZE1 32768
#define BUFSIZE2 32*66564
#define NMAX 20
 
__global__ void test(char * BUFFOR1,int len,int print_if,int limit) {
  int tid = threadIdx.x+(blockIdx.x*blockDim.x);
    
  /*for (int j = 0;j<len;j++){
      BUFFOR[j]=BUFFOR1[j+tid*len]
    }*/
  //char * BUFFOR=BUFOR1

  if (tid>=limit){
    printf("%d return",tid);
    return;
  } 
  char *BUFFOR = BUFFOR1 + tid * len;
 
  //printf("%d,%s\n",tid,BUFFOR);
  int i,j,k,k3,k4,L,L1,z;
  double /* lambda, */ eps,g,h,ma,mn,norm,s,t,u,w;
  int cond;
  double d[NMAX+1], e[NMAX+1], e2[NMAX+1], Lb[NMAX+1];
  double x[NMAX+1];
  double a[NMAX*(NMAX-1)/2 + NMAX + 1];
  int n;
  int bit, poz, poz2;
  bit = 32;
  poz = 1;
  poz2 = 1;
  n = BUFFOR[0] - 63;
  a[0] = 0.0;
  for (i = 0; i < n; i++)
   for (j = 0; j<=i; j++)
    {
	  if (i==j) {a[poz2++] = 0; }
	  else {
        if (bit == 0) { bit = 32;  poz++; }
        if ((BUFFOR[poz] - 63) & bit)
                { a[poz2++] = 1; }
               else
                { a[poz2++] = 0; }
        bit = bit >> 1;
	  }
    }
  //dopelniony below
  /*for (i = 0; i < n; i++)
   for (j = 0; j<=i; j++)
    {
	  if (i==j) {a[poz2++] = 0; }
	  else {
        if (bit == 0) { bit = 32;  poz++; }
        if ((BUFFOR[poz] - 63) & bit)
                { a[poz2++] = 0; }
               else
                { a[poz2++] = 1; }
        bit = bit >> 1;
	  }
    }*/
 
  int k1 = 1;
  int k2 = n;
  // if ((1<=k1) && (k1<=k2) && (k2<=n))
   {
    i = 0;
    for (L=1;L<=n;L++) { i += L; d[L] = a[i]; /* printf("%Lf ",a[i]); */ }
 
    for (L=n;L>=2;L--)
     {
      i--; j = i; h = a[j]; s = 0;
      for (k=L-2;k>=1;k--) { i--; g = a[i]; s += g*g; }
      i--;
      if (s == 0) { e[L] = h; e2[L] = h*h; a[j] = 0.0; }
       else
        {
          s += h*h; e2[L] = s; g = sqrtf(s); if (h>=0.0) g=-g;
          e[L] = g;
          s = 1.0 / (s-h*g);
          a[j] = h - g; h = 0.0; L1 = L - 1; k3 = 1;
          for (j=1;j<=L1;j++)
           {
             k4 = k3; g = 0;
             for (k=1;k<=L1;k++) { g +=a[k4]*a[i+k]; 
                                   if (k<j)  z = 1; else z = k;
                                   k4 += z; }
             k3 += j; g *= s; e[j] = g; h += a[i+j]*g;
           }
          h *= 0.5*s; k3 = 1;
          for (j=1;j<=L1;j++)
           {
             s = a[i+j]; g = e[j]-h*s; e[j] = g;
             for (k=1;k<=j;k++) { a[k3] += -s*e[k]-a[i+k]*g; k3++; }
           }
        }
      h = d[L]; d[L] = a[i+L]; a[i+L] = h;
     }
    h = d[1]; d[1] = a[1]; a[1] = h; e[1] = 0.0; e2[1] = 0.0; s = d[n];
    t = fabsf(e[n]); mn = s - t; ma = s + t;
    for (i=n-1;i>=1;i--)
     {
      u = fabsf(e[i]); h = t + u; t = u; s = d[i]; u = s - h;
      if (u < mn) mn = u;
      u = s + h;
      if (u > ma) ma = u;
     }
    for (i=1;i<=n;i++) { Lb[i] = mn; x[i] = ma; }
    norm = fabsf(mn); s = fabsf(ma);
    if (s>norm) norm = s;
    w = ma; /* lambda = norm; */ eps = 7.28e-17*norm;
    for (k=k2;k>=k1;k--)
     {
      /* eps = 7.28e-17*norm; */ s = mn; i = k;
      do {cond = 0; g = Lb[i];
         if (s < g) s = g; else { i--; if (i>=k1) cond = 1; }
      } while (cond);
      g = x[k];
      if (w>g) w = g;
      while (w-s>2.91e-16*(fabsf(s)+fabsf(w))+eps)
       {
         if (floorf(w+10e-5)<s-10e-5) return;  // przedział nie zawiera liczby całkowitej
         L1 = 0; g = 1.0; t = 0.5*(s+w);
         for (i=1;i<=n;i++)
          {
            if (g!=0)  g = e2[i] / g; else g = fabsf(6.87e15*e[i]);
            g = d[i]-t-g;
            if (g<0) L1++;
          }
         if (L1<k1) { s = t; Lb[k1] = s; }
          else
           { if (L1<k)
               {
                 s = t; Lb[L1+1] = s;
                 if (x[L1]>t) x[L1] = t;
               }
              else w = t;
           }
      } // while
      u = 0.5*(s+w); x[k] = u;
	  if  (!(( ceilf(u) - u  < 10e-5 ) || ( u - floorf(u) < 10e-5 ))) { return; };
    }
  }   
  //for (i=0;i<=n;i++) printf("%f ",x[i]);
 
  //printf("GPU:[%s]\n",BUFFOR);
  if (print_if>0)
  printf("%s\n",BUFFOR);
}
 
 
/** 
 
*/
int main(int argc, char *argv[])
{
 
  char BUFFOR[BUFSIZE];
  char BUFFOR1[BUFSIZE1];
  char BUFFOR2[BUFSIZE2];
 
  char * cuda_bufor;
  char * cuda_bufor1;
  char * cuda_bufor2;
  int print_if = 1; 
  int active = 1;
  if (argc>1) {print_if=strtol(argv[1],NULL,10);}  
  if (argc>2) {active=strtol(argv[2],NULL,10);}  
 
  int i = 0,j=0,len;
  cudaMalloc((void**)&cuda_bufor,BUFSIZE1);
  cudaMalloc((void**)&cuda_bufor1,BUFSIZE);
  cudaMalloc((void**)&cuda_bufor2,BUFSIZE2);
  int k=0,iter=0,now=1;
  double start, fin,full_time1=0.0,full_time2=0.0,full_time3=0.0,full_time4=0.0,full_time5=0.0;
  while (fgets(BUFFOR,BUFSIZE-1,stdin)) {
    printf("%d\n",i);
    len = strlen(BUFFOR);
    for (int j1 = 0;j1<len-1;j1++){
      BUFFOR1[i*len+j1]=BUFFOR[j1];
      BUFFOR2[j*len+j1]=BUFFOR[j1];
    }
    BUFFOR1[i*len+len-1]='\0';
    BUFFOR2[j*len+len-1]='\0';
    i+=1;
    j+=1;
     // if (eigensymmatrix(BUFFOR)) 
    //printf("Main:%s",BUFFOR);
	  //BUFFOR[glen]='\0';
    if (active==1){
    start = omp_get_wtime();
    cudaMemcpy(cuda_bufor1,BUFFOR,BUFSIZE,cudaMemcpyHostToDevice);
    test<<<1,1>>>(cuda_bufor1,len,0,1);
    cudaDeviceSynchronize();	
    fin = omp_get_wtime();
    full_time1+=fin-start;
    }
    if (i==1024){
      start = omp_get_wtime();
      cudaMemcpy(cuda_bufor,BUFFOR1,BUFSIZE1,cudaMemcpyHostToDevice);
      test<<<1,1024>>>(cuda_bufor,len,1,1024);
      cudaDeviceSynchronize();	
      fin = omp_get_wtime();
      full_time2+=fin-start;
      

      start = omp_get_wtime();
      cudaMemcpy(cuda_bufor,BUFFOR1,BUFSIZE1,cudaMemcpyHostToDevice);
      test<<<1024,1>>>(cuda_bufor,len,0,1024);
      cudaDeviceSynchronize();	
      fin = omp_get_wtime();
      full_time4+=fin-start;

      start = omp_get_wtime();
      cudaMemcpy(cuda_bufor,BUFFOR1,BUFSIZE1,cudaMemcpyHostToDevice);
      test<<<32,32>>>(cuda_bufor,len,0,1024);
      cudaDeviceSynchronize();	
      fin = omp_get_wtime();
      full_time5+=fin-start;
      i=0;
    }
    if (j==66564){
      start = omp_get_wtime();
      cudaMemcpy(cuda_bufor2,BUFFOR2,BUFSIZE2,cudaMemcpyHostToDevice);
      test<<<256,256>>>(cuda_bufor2,len,print_if,66564);
      cudaDeviceSynchronize();	
      fin = omp_get_wtime();
      full_time3+=fin-start;
      j=0;
      k+=1;


      iter+=1;
      if (iter>=now)
      {
      printf("dla 66564 * %d\n",iter);

      if (active==1)
      printf("czas dla pojedynczych watkow = %f \n",full_time1 );
      printf("czas dla 1024 watkow = %f \n",full_time2 );
      printf("czas dla 1024 blokow z jednym watkiem = %f \n",full_time4 );
      printf("czas dla 32 blokow z 32 watkami = %f \n",full_time5 );
      printf("czas dla 256 blokow z 256 watkami = %f \n",full_time3 );
      now=now*2;
      }
    }
 
  } // while
  start = omp_get_wtime();
  cudaMemcpy(cuda_bufor,BUFFOR1,BUFSIZE1,cudaMemcpyHostToDevice);
  test<<<1,i>>>(cuda_bufor,len,0,i);
  cudaDeviceSynchronize();	
  fin = omp_get_wtime();
  full_time2+=fin-start;
  
  start = omp_get_wtime();
  cudaMemcpy(cuda_bufor,BUFFOR1,BUFSIZE1,cudaMemcpyHostToDevice);
  test<<<i,1>>>(cuda_bufor,len,0,i);
  cudaDeviceSynchronize();	
  fin = omp_get_wtime();
  full_time4+=fin-start;

  start = omp_get_wtime();
  cudaMemcpy(cuda_bufor,BUFFOR1,BUFSIZE1,cudaMemcpyHostToDevice);
  test<<<32,32>>>(cuda_bufor,len,0,i);
  cudaDeviceSynchronize();	
  fin = omp_get_wtime();
  full_time5+=fin-start;

  start = omp_get_wtime();
  cudaMemcpy(cuda_bufor2,BUFFOR2,BUFSIZE2,cudaMemcpyHostToDevice);
  test<<<256,256>>>(cuda_bufor2,len,print_if,j);
  cudaDeviceSynchronize();	
  fin = omp_get_wtime();
  full_time3+=fin-start;
  j=0;
  k+=1;


  if (active==1)
  printf("czas dla pojedynczych watkow = %f \n",full_time1 );
  printf("czas dla 1024 watkow = %f \n",full_time2 );
  printf("czas dla 1024 blokow z jednym watkiem = %f \n",full_time4 );
  printf("czas dla 32 blokow z 32 watkami = %f \n",full_time5 );
  printf("czas dla 256 blokow z 256 watkami = %f \n",full_time3 );
  cudaFree(cuda_bufor);
  cudaFree(cuda_bufor1);
  return EXIT_SUCCESS;
}
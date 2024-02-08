#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define row 32
#define col  64
#define itmax 100000
int main()
{
// FILE *fp;
FILE *fp2;
fp2 = fopen("outhex.txt","w");
fclose(fp2);

int stp,m,n;
int dep;
int diff;
float ETS;
int evnt3n;
int ecnt3n=0;
int proc3n[2]={0,0};
int fevnt;
int fdcnt=0; 
int occ[row][col];
int monomer(int occ[row][col],int ndispx[6],int ndispy[6]);
float time(float total,int process);
float edge(int x,int y,int occ[row][col]);
float probability[4][4]={{0.014,0.014,0.014,0.000},      //boron  //type 1
                         {4.610,2.220,0.014,0.014},
                         {8.170,4.530,2.450,2.140},
                         {0.000,7.280,4.150,2.560}};

float probability2[4][4]={{0.018,0.018,0.018,0.000},      //nitrogen //type 2
                          {4.610,2.150,0.018,0.018},
                          {5.410,4.120,2.810,1.860},
                          {0.000,3.850,3.100,2.470}};
   
                                      /* disp vect (NN): type 1 [0:2] & 2 [3:5] atom site*/
int ndispx[6]={ 0, 0, 1,-1, 0, 0};
int ndispy[6]={ 1,-1,-1, 1, 1,-1}; 
                                        /* disp vect (NNN) both type 1/2 atom site*/
int dispx[6]={-1,0,-1,+1,0,+1};
int dispy[6]={+2,+2,0,0,-2,-2};
                                        // saving NNN loop information in different processes 
int Xcorr[5][5][3000];
int Type[5][5][3000];
int Ycorr[5][5][3000];
float energy[5][5][3000];
int Xcorrs[5][5][3000];
int Ycorrs[5][5][3000];
int tsx,tsy;
int idiff;
int tSwing=0;
float total;
int process;
int breakit=0;
int desorp=0;

//  neib[i][j]={};   //occupancy and   matrics
for (int i=0;i<row;i++){for (int j=0;j<col;j++){ occ[i][j]=0;}}
 dep=0;
 diff=0;
 float xyz;
         // the main loop
         // add a output at end of this loop
 for (int i=1;i<=itmax;i++){
   printf("\n\n new iteration number : %d",i-1);
      process=0;
      xyz= monomer(occ,ndispx,ndispy);

      float Pf = 1 / (1+ (xyz / (row*col))*2*pow (10,5)) ;
      // float Pf =0.1;
      printf("\nthis is monomer %f and probability of deposition %f",xyz,Pf);
      float r1;
      r1=(float) random()/RAND_MAX;
      // printf("\nrandom number is : %f and  r_ad is : %f",r1,r_ad);
      if (r1>Pf){process=1;}
      else{process=2;}                //p=1 diffusion
      if ((dep-desorp)==615){break;}

      // printf("\nselected process if p=1 is diffision : %d",p); 

    if(process==1){                          //diffusion loop
   int cntdiff=0;
   total=0;
   int swings[2]={};
   int swing[2]={};
   float glass[5][5][2]={0};
   for (int x=0;x<row;x++){
   for (int y=0;y<col;y++){   
      if(occ[x][y]!=0){ 
            // printf("\niteration of : %d",x+y);
            int sxs,sys;
            int m=0,EN=0;                  // NN site around the diffusion site
            stp = occ[x][y];
            for (int j=0;j<3;j++){
                  sxs=(x+row+ndispx[3*(stp-1)+j])%row;
                  sys=(y+col+ndispy[3*(stp-1)+j])%col;  
                  if (occ[sxs][sys]!=0){
                  swings[0]=sxs;
                  swings[1]=sys;
                  m++;
                  }}
            if (m==0){
               float rate;
               int c2=glass[4][3][1];
                  Type[4][3][c2]=0;
                  Xcorrs[4][3][c2]=x;
                  Ycorrs[4][3][c2]=y;
                  Xcorr[4][3][c2]=x;
                  Ycorr[4][3][c2]=y;
                  if(stp==2){rate=pow(10,12)*exp(-0.44*7.252859);}
                  else{rate=pow(10,12)*exp(-0.78*7.252859);}
                  total = total + rate;
                  glass[4][3][0]=glass[4][3][0]+ rate;
                  glass[4][3][1]++;

            }            

            for (int j=0;j<6;j++){               //empty NNN sites i,j 

                  int sx,sy,site;
                  sx=(x+row+dispx[j])%row;
                  sy=(y+col+dispy[j])%col; 
                  if(occ[sx][sy]==0){
                  int n=0;
                  float rate=0;
                  int tsi,tsj;
                  for (int k=0;k<3;k++){ 
                     tsi=(sx+row+ndispx[3*(stp-1)+k])%row;
                     tsj=(sy+col+ndispy[3*(stp-1)+k])%col;
                     if (occ[tsi][tsj]!=0){
                        swing[0]=tsi;
                        swing[1]=tsj;
                        n++;
                        };}

                     if (m==1 && n==1){

                     if (swings[0]==swing[0] && swings[1]==swing[1] ){
                        int c1=glass[4][4][1];
                        Type[4][4][c1]=stp;
                        Xcorrs[4][4][c1]=x;
                        Ycorrs[4][4][c1]=y;
                        Xcorr[4][4][c1]=sx;
                        Ycorr[4][4][c1]=sy;
                        if(stp==1){rate=pow(10,12)*exp(-0.32*7.252859);}
                        else{rate=pow(10,12)*exp(-0.39*7.252859);}
                        total = total + rate;
                        glass[4][4][0]=glass[4][4][0]+ rate;
                        glass[4][4][1]++;
                     }}
                     
                     if (stp==1){
                     rate =pow(10,12)* exp(-probability[m][n]*7.252859);   //pow(10,12)*
                     // printf("\n the stp is %d so rate for boran",stp); 
                     }
                     else {
                     rate =pow(10,12)*exp(-probability2[m][n]*7.252859);}

                     glass[m][n][0]=glass[m][n][0]+rate;
                     int c=glass[m][n][1];
                     Type[m][n][c]=stp;
                     Xcorrs[m][n][c]=x;
                     Ycorrs[m][n][c]=y;
                     Xcorr[m][n][c]=sx;
                     Ycorr[m][n][c]=sy;
                     total = total + rate;
                     glass[m][n][1]++;

                     
                     ;}}


                  // edge and kink addition
                  edge(x,y,occ);
                  }}}

         printf("\nthis is 00 %f,%f",glass[0][0][1],glass[0][0][0]);
         printf("\tthis is 01 %f,%f",glass[0][1][1],glass[0][1][0]);
         printf("\nthis is 02 %f,%f",glass[0][2][1],glass[0][2][0]);
         printf("\tthis is 10  %f,%f",glass[1][0][1],glass[1][0][0]);
         printf("\nthis is 11 %f,%f",glass[1][1][1],glass[1][1][0]);
         printf("\tthis is 12  %f,%f",glass[1][2][1],glass[1][2][0]);
         printf("\nthis is 13 %f,%f",glass[1][3][1],glass[1][3][0]);
         printf("\tthis is 44  %f,%f",glass[4][4][1],glass[4][4][0]);
         printf("\n\nthis is 43  %f,%f",glass[4][3][1],glass[4][3][0]);






               if (total>0){
                  float r2,Renergy;
                  r2=(float) random()/RAND_MAX;
                  Renergy=r2*total;
                  int event[2]={},breakt=0;
                  float eenergy=0;

                  for (int p=0;p<5;p++){
                  for (int q=0;q<5;q++){
                     eenergy = eenergy + glass[p][q][0];
                        if (breakt==0 && eenergy>Renergy){
                           event[0]=p;event[1]=q;
                           breakt=1;
                        if(p==4 && q==3 ){desorp++;}
                        if(p==4 && q==4 )tSwing++;}

                           }}
                  printf("\nChosen one is %d,%d",event[0],event[1]);
                  float r3;
                  int l;
                  r3=(float) random()/RAND_MAX;
                  l=r3*glass[event[0]][event[1]][1];
                  printf("\ntotal %f and chosen one %d and random %f",glass[event[0]][event[1]][1],l,r3) ;
                  printf("\natom type %d shift from (%d,%d) to (%d,%d)",Type [event[0]][event[1]][l],Xcorrs[event[0]][event[1]][l],Ycorrs[event[0]][event[1]][l],Xcorr[event[0]][event[1]][l],Ycorr[event[0]][event[1]][l]);

                  occ[ Xcorrs[event[0]][event[1]][l] ] [Ycorrs[event[0]][event[1]][l] ] = 0;
                  occ[ Xcorr[event[0]][event[1]][l] ] [Ycorr[event[0]][event[1]][l] ] = Type [event[0]][event[1]][l] ;  
                  diff++;   
               }}
//                                           deposition loop
    if(process==2){
      float r;
      int si,sj,site;
      r=(float) random()/RAND_MAX;
      site=row*col*r;
      si=site%row;
      sj=site/row;
      stp=(site/row)%2+1;
      int xcorr[10]={};
      int ycorr[10]={};

      fevnt=0;
      if(occ[si][sj]==0){
         occ[si][sj]=stp;
         dep++; 
      }
      else{
         int cnt3n=0;
         for (int j=0;j<6;j++){
            int sxd,syd;
             sxd=(si+row+dispx[j])%row;
             syd=(sj+col+dispy[j])%col;   
             if(occ[sxd][syd]==0){
                xcorr[cnt3n]=sxd;           
                ycorr[cnt3n]=syd;          
                cnt3n++;}}
         if (cnt3n > 0){
            r=(float) random()/RAND_MAX;
            int m = r*cnt3n;
            printf("the random %f total emply %d and choosen one %d (%d,%d) to (%d,%d)",r,cnt3n,m,si,sj,xcorr[m],ycorr[m]);
            occ[xcorr[m]][ycorr[m]]=stp;
            dep++;
         }
         else{
            fevnt++;
         }
    }}
//                                        Printing output in File
      int flip=0;

      // ///////////////////////// for square lattice 

      // fp = fopen("outsq.txt","a");
      // if (fp == NULL)
      // {printf("Error opening the file" ); return -1;}

      // fprintf(fp,"ITEM: TIMESTEP \n %d\t%d\n",i,i);
      // fprintf(fp,"ITEM: NUMBER OF ATOMS \n %d\n",row*col);
      // fprintf(fp,"ITEM: BOX BOUNDS\n 0 \t%d\n 0 \t%d\n-0.5\t0.5\n",row,col);
      // fprintf(fp,"ITEM: ATOMS type x y\n");  

      // for (int a=0;a<row;a++){
      //     for (int b=0;b<col;b++){
      //       fprintf(fp,"%d %d %d\n",occ[a][b],a,b);
      // }}
      //    fclose(fp);

      // ///////////////////// for hexagonal layout

      fp2 = fopen("outhex.txt","a");
      if (fp2 == NULL)
      {printf("Error opening the file" ); return -1;}
      
      fprintf(fp2,"ITEM: TIMESTEP \n %d\t%d\n",i,i);
      fprintf(fp2,"ITEM: NUMBER OF ATOMS \n %d\n",row*col);
      fprintf(fp2,"ITEM: BOX BOUNDS\n 0 \t%d\n 0 \t%d\n-0.5\t0.5\n",row,col);
      fprintf(fp2,"ITEM: ATOMS type x y\n");  

      for (int a=0;a<row;a++){
          for (int b=0;b<col;b++){
            if(b%2==0 && flip==0){flip=flip+1;}
            else if(b%2==0 && flip==1){flip=flip-1;}
            // printf("%d",flip);
            if(flip==1){fprintf(fp2,"%d %d %d\n",occ[a][b],2*a+1,2*b);}
            else{fprintf(fp2,"%d %d %d\n",occ[a][b],2*a,2*b);}
      }}
      fclose(fp2);
      // time update
      ETS = ETS + time(total,process);
      printf("\nEstimated time of simulation %8.4e ns",ETS);
      }
      printf("\nno.failed process=%d\n number of swing = %d\n total deposition are : %d\n total diffusion are : %d\nEstimated time of simulation in nanoseconds %2.2e ns \n total desorption are : %d\n",fevnt,tSwing,dep,diff,ETS,desorp); 
}   //end of main loop


int monomer(int occ[row][col],int ndispx[6],int ndispy[6]){
   int m=0;
   for (int x=0;x<row;x++){
   for (int y=0;y<col;y++){
   int n=0;   
   if(occ[x][y]!=0){ 
      int tsi,tsj,stp;
      stp = occ[x][y];

      for (int k=0;k<3;k++){
         tsi=(x+row+ndispx[3*(stp-1)+k])%row;
         tsj=(y+col+ndispy[3*(stp-1)+k])%col;
         if(occ[tsi][tsj]==0){n++;}}

      if(n==3){m++;}
      
      }}}
   return m;
}
float time(float total, int process){
   float rate=0;
   rate = total + pow(10,12)* exp(-0.014*7.252859);
   float r,timestep;
   r=(float) random()/RAND_MAX;
   timestep= (- pow(10,9)*log(r))/rate;
   // printf("\n%f,%f,%f,%8.4e",total,rate,r,timestep);
   return timestep;
}


// int kink(int site, int stp, int occ[row][col]) {
//    int sy,sx,cy,cx;
//    int itr=0, ip;
//    int start=0,smax=5,pcnt;
//    int drow=(stp-1)*2;
//    int dcol;   
//                                     /* disp vect: (x/y) KINK[0:4],type 1 [0:1][0:4] & 2 [2:3][0:4]*/
//    int kdisp[4][18]={{-1,-1,-1, 0, 0, 0, 0,-1,-1,-1, 0, 1, 1, 1, 0},
//                      { 1, 2, 3, 2, 1,-1,-2,-1, 0, 1, 1, 0,-1,-2,-1},
//                      { 0,-1,-1,-1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0},
//                      {-1, 0, 1, 2, 1, 1, 2, 1, 0,-1,-1,-2,-3,-2,-1}};
//    int patt=0;
    
//    sy=site/col;
//    sx=site%col;


//                                               /*  Search for KINK(0)? */
//        pcnt=0;  
//        for (int m=0; m<3 ;m++){
//        for (ip=start;ip<start+smax;ip++){
//                                                 /*coordinates relative to the empty site */
//             cx=sx+kdisp[drow][ip + m*5];
//             cy=sy+kdisp[drow+1][ip + m*5];
      
//             if (cx<0) cx=cx+col;
//             if (cy<0) cy=cy+row; 
        
//             if (cx>=col) cx=cx-col;
//             if (cy>=row) cy=cy-row;   
//                                                 /*matrix notation, (row,col) for occ[][]*/ 
//             if (occ[cy][cx]!=0){ 
//                pcnt++;
//             }
//        }
//        if (pcnt==5){patt++;}
//        }
//    return patt;
// }


float edge(int x, int y, int occ[row][col]) {
   int cy,cx;
   int itr=0, ip;
   int start=0,smax=6,pcnt;
   int drow=(occ[x][y]-1)*2;    //(0,2)
   int dcol;

   int edisp[4][18] = {   {-1,-1,-2,-2,-2,-1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0},
                          { 1, 0, 1, 2, 3, 2, 1, 2, 3, 2, 1, 0,-1,-2,-3,-4,-3,-2},
                          { 0,-1,-1,-1, 0, 0, 0, 0, 0,-1,-1,-1, 1, 1, 2, 2, 2, 1},
                          { 1, 2, 3, 4, 3, 2,-1,-2,-3,-2,-1, 0,-1, 0,-1,-2,-3,-2}  };


   int patt=0;
    
                                              /*  Search for edge(1)? */
       pcnt=0;
       for (int m = 0; m<3 ;m++){
       for (ip=0; ip<smax; ip++){
                                                /*coordinates relative to the empty site */
            cx=x+row+edisp[drow][ip + m*6]%row;              //0 or 2
            cy=y+col+edisp[drow+1][ip + m*6]%col;            //1 or 3
    
            if (occ[cy][cx]!=0){pcnt++;}
       }
     if (pcnt==6){patt+=1;}
   }
   return patt;
}

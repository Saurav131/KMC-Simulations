#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#define row 70
#define col 80
#define temp 1073
#define iterations pow(10,4)
#define feed_w  2.5*pow(10,-4)
#define feed_se  2.5*pow(10,-4)
#define Boltzman_constant 1.380649*pow(10,-23)
#define ev_to_joule 1.602176634*pow(10,-19)
#define print_step pow(10,1)

// funcions used
int NN(int x, int y, int occ[row][col]);
int diffusion(int x ,int y , double local_event[10][2], int occ[row][col],float beta_ev);
int deposition(int occ[row][col],int type);
float samay(float total);

// Arrays used
int NNcorrx[3]={0};
int NNcorry[3]={0};
int NNtype[3]={0};
int swing[2];
 // first square bracket is for all the different events and second bracket is for storing energy for local_event, 
 // total possible local_event in the local_event class.
int point[10][10000][4]={0};

int main(){

MPI_Init(NULL, NULL);
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);


clock_t Seconds;
Seconds = clock();

FILE *fp;
fp = fopen("outhex.txt","w");
fclose(fp);
int m=0;
int occ[row][col]={0};
                                //this is silicon atom 
float beta_ev = (ev_to_joule/(Boltzman_constant*temp));  // value of 1/kbt in units of ev^-1
float Dep_W = 1 / (feed_w * row * col);
float Dep_Se = 1 / (feed_se * row * col);
float Edge_Dep_W = 1 /  (exp(0.55*beta_ev) * feed_w * row * col);
float Edge_Dep_Se = 1 / (exp(0.55*beta_ev) * feed_se * row * col);
float Kink_Dep_W = 1 / (exp(0.78*beta_ev) * feed_w * row * col);
float Kink_Dep_Se = 1 /(exp(0.78*beta_ev) * feed_se * row * col);
     
                                          //seed
if(world_rank==0){
    occ[row/2][col/2]=1;
    for (int ir=1;ir<50;ir++){
        int count=0;
    for(int k=0;k<row;k++){
    for(int j=0;j<col;j++){
    if( occ[k][j]==(ir%2+1) && count<(3*ir) && ((pow((k-(row/2)),2) + pow((j-(col/2)),2)) <145)) {
        m=NN(k,j,occ);
        if(m==0){
            for(int i=0;i<3;i++){
                if(occ[NNcorrx[i]][NNcorry[i]]==0){
                count++;
                occ[NNcorrx[i]][NNcorry[i]]=NNtype[i];}
                }}
        if(m==1 && count<(3*ir)){
            for(int i=0;i<2;i++){
                if(occ[NNcorrx[i]][NNcorry[i]]==0){
                count++;
                occ[NNcorrx[i]][NNcorry[i]]=NNtype[i];
                }}}
        if(m==2 && count<(3*ir)){
            for(int i=0;i<1;i++){
            if(occ[NNcorrx[i]][NNcorry[i]]==0){
                count++;
                occ[NNcorrx[i]][NNcorry[i]]=NNtype[i];
            }}}
    }}}}
printf("\n\n ------------------- SEED ----------------------");

}
    
    
MPI_Barrier(MPI_COMM_WORLD);
MPI_Bcast(occ, row * col, MPI_INT, 0, MPI_COMM_WORLD);


                  
                        // end of loop to create seed
                        // ////////////////Kmc Loop

float ETS=0;
int surfacedep=0, edgedep=1, kinkdep=1 ,type;


for(int itr=0;itr<iterations;itr++){
    double total=0;
    double local_event[10][2] = {0};
    double global_event[10][2] = {0};
    // printf("world rank :%d and size :%d ",world_rank,world_size);

// building the local_event catalogue


    for(int i=0;i<row;i++){
    for(int j=0;j<col;j++){
        float Rate;
        if(occ[i][j]!=0){   
            if ((i * j) % world_size == world_rank) {
                // printf("Rank %d handling i=%d, j=%d\n", world_rank, i, j);
                diffusion(i,j,local_event,occ,beta_ev);
    }}}}
MPI_Barrier(MPI_COMM_WORLD);
MPI_Reduce(local_event, global_event, 20, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);



if (world_rank ==0 ){
printf("\n iteration =%d , ETS= %8.4e",itr,ETS);
    for(int i=0;i<5;i++){
  //      printf("the local_event class %d, rate : %f , number : %f\n",i,global_event[i][0],global_event[i][1]);
        total=total+global_event[i][0];}

                                        //time Update
    ETS = ETS + samay(total);                        
   // printf("the total :%f  and ETS :%f   and timestep :%f",total,ETS,samay(total));                                    // selecting one local_event

    // surface, ring or edge diffusion



    if(ETS > surfacedep*Dep_Se ){
       // printf("\n\n surface deposition");
        type=0;
        deposition(occ,type);
        surfacedep++;}

    else if(ETS > edgedep*Edge_Dep_Se){
     //   printf("\n\n edge deposition");  
        type=1;
        deposition(occ,type);
        edgedep++;}

    else if (ETS > kinkdep*Kink_Dep_Se){
        printf("\n\n kink deposition");
        type=2;
        deposition(occ,type);
        kinkdep++;}

   // printf("\n");

    if(total>0){
        float r2,Renergy;
        int n,breakit=0;
        r2=(float) random()/RAND_MAX;
        Renergy=r2*total;
        float eenergy=0;
       // printf("random rate : %8.4f\n",Renergy);
        for(int p=0;p<10;p++){
            eenergy=eenergy+global_event[p][0];
            if(eenergy>Renergy && breakit==0){
                n=p;
                breakit=1;}}

        float r3;
        int l;
        r3=(float) random()/RAND_MAX;
        l=global_event[n][1];
      //  printf("%d\t",l);
        l=r3*l;
       // printf("%d\t",l);

     //   printf("\nthe point of diffusion is (%d,%d) and the local_event class is %d, and local_event no %d  \n" ,point[n][l][0],
      //  point[n][l][1], n ,l );
        occ[point[n][l][0]][point[n][l][1]]=0;
        occ[point[n][l][2]][point[n][l][3]]=(point[n][l][3]%2)+1;   
    }


    int i2;
    i2 = itr/print_step;
    i2 = i2*print_step;
    if(itr==1||i2==itr){
        int flip=0;
                                        // for hexagonal layout
        fp = fopen("outhex.txt","a");
        if (fp == NULL){printf("Error opening the file" ); return -1;}
        fprintf(fp,"ITEM: TIMESTEP \n %d\t%d\n",itr,itr);
        fprintf(fp,"ITEM: NUMBER OF ATOMS \n %d\n",row*col);
        fprintf(fp,"ITEM: BOX BOUNDS\n 0 \t%d\n 0 \t%d\n-0.5\t0.5\n",row,col);
        fprintf(fp,"ITEM: ATOMS type x y\n");
        for (int a=0;a<row;a++){
            for (int b=0;b<col;b++){
                if(b%2==0 && flip==0){flip=flip+1;}
                else if(b%2==0 && flip==1){flip=flip-1;}
                // printf("%d",flip);
                if(flip==1){fprintf(fp,"%d %d %d\n",occ[a][b],2*a,2*b);}
                else{fprintf(fp,"%d %d %d\n",occ[a][b],2*a+1,2*b);}
      }}
      fclose(fp);
    }
}
MPI_Bcast(occ, row * col, MPI_INT, 0, MPI_COMM_WORLD);
}                                      //iteration loop

                                       // Printing output
if(world_rank==0){
    for(int i=1;i<row;i++){
        printf("\n");
        for (int j=1;j<col;j++){
        printf("%d",occ[i][j]);
        };
    };
    // the simulation summary 
    printf("\n total surface kink and edge depostions are %d, %d, %d respectively",surfacedep,kinkdep,edgedep);
    // printf("%f \n\n",beta_ev);
    printf("\n%f\t%f\t%f\t%f\t%f\t%f\n",Dep_Se,Dep_W,Edge_Dep_Se,Edge_Dep_W,Kink_Dep_Se,Kink_Dep_W);
    printf("\n %f",feed_se);
    Seconds = clock()-Seconds;
    printf("\nthe total time taken in seconds are : %f\n",((float)Seconds)/CLOCKS_PER_SEC );
}
MPI_Finalize();
}
                                   //main loop ends

int deposition(int occ[row][col],int type){
    int rate,a=0,z=0,zz,x=0;
    float r;
    for(int i=0;i<row;i++){
    for(int j=0;j<col;j++){
        if(occ[i][j]==0 && NN(i,j,occ)==type){
            z++;
            }}}
    
    r=(float)random()/RAND_MAX;
    zz=z*r;
        
    for(int i=0;i<row && x==0 ;i++){
    for(int j=0;j<col && x==0 ;j++){
        
        if(occ[i][j]==0 && NN(i,j,occ)==type){a++;}
        
        if(zz==a){
            x=1;
            occ[i][j]=(j%2)+1;}
            }}

    return 0;
    }
// Few additional fuction are removed which are required to run these simulation. Contact me on Sauravmittal131@gmail.com
        s

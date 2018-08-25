
//传统方法计算适应值用支距离和层接近度//计算了评价次数和成功率

//目标路径路径2：等边三角形 1-2-3-5-6-7

#include <string.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#define SIZE 50//控制种群规模

#define SIZE1 128 //控制测试数据范围大小0-7000
#define NUM  100//实验次数，多次试验取平均值

#define CSIZE  33 //!!!!!!!!!!!!!!!!!!!代表3个数 需要的二进制存储位数
#define G 100  //控制进化代数
#define CRP  0.9//交叉算子
#define MUP  0.3//变异算子

#define PL 6//目标路径的节点数

#define N 3 // N  代表数组长度  year month day 3个数
//#define PN  3 //PN代表目标路径个数@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//#define CN  28 //CN 代表路径编码长度 CN=N-1  + N-2 + …… + 1

//#define FITNESS 0.19999//记录适应值大于FITNESS的就输出，当找到5个时，就退出程序
//#define DIF  1.0  //要求最少可以有相同位数的百分比

FILE *fp ;
#include "stdlib.h"
#include "time.h"

#define LEN  10//目标路径长度

//-----------------------------------------------------------------------------------
int s[PL];//统计目标路径的4个节点被执行的次数,在每一代种群执行之前清0，每次种群运行结束就统计适应值

//-----------------------------------------------------------------------------------


struct individual
{   char chrom[CSIZE+1];//染色体编码
	int p[N];    //测试数据
	int tpath[LEN];//穿越路径
	int branch[PL];//存放分支距离
	double fitness;//适应值函数
};
//==========================================================================
void selection(struct  individual population[])//
{
	int i,j,index;
	double p,sum=0.0;
	double cfitness[SIZE];  //cumulative fitness value
	struct individual newpopulation[SIZE];
	float max;
	//calculate relative fitness
	max=population[0].fitness;
	for(i=0;i<SIZE;i++)
	{
		if(max<population[i].fitness) max=population[i].fitness;
		sum+=population[i].fitness;
		//printf("population[%d].fitness=%f\n",i,population[i].fitness);//--------------------------------------------
	}
	for(i=0;i<SIZE;i++)
	{
		cfitness[i]=population[i].fitness/sum;
	}

	//calculate cumulation fitness
	for(i=1;i<SIZE;i++)
	{
		cfitness[i]=cfitness[i-1]+cfitness[i];
	}

	//selection operation
	for(i=0;i<SIZE;i++)
	{
		p=rand()%1000/1000.0;
		index=0;
		while(p>cfitness[index])
		{
			index++;
		}

		for(j=0;j<N;j++)
		newpopulation[i].p[j]=population[index].p[j];
	}
	for(i=0;i<SIZE;i++)
	{
		for(j=0;j<N;j++)
		population[i].p[j]=newpopulation[i].p[j];
	}
}
//==================================================================================

//编码函数
void code(int n[],char c[])
{
	int t,l,i,j;
	for(l=0;l<CSIZE;)
	{
		i=l/(CSIZE/N);t=n[i];j=0;//j为计数器，确保每个数据对应7位2进制数
	    while(t!=0)
		{
			if(t%2==1) {c[l++]='1';j++;}
			else  {c[l++]='0';j++; }
			t/=2;
		}
		while(j<(CSIZE/N))
		{ c[l++]='0';j++; }

	}

	c[l]='\0';
}
//================================================================================
//解码函数
void uncode(int n[],char c[])
{
	int t,l,i,j,m;
	for(l=0;l<CSIZE;)
	{
		i=l/(CSIZE/N);t=0;m=1;j=0;//j为计数器，确保每个数据对应7位2进制数
	    while(j<(CSIZE/N))
		{
			t=t+(c[l++]-'0')*m;
			m*=2;j++;
		}
	    n[i]=t%(SIZE1+1);
	}
}
//===========================================================================
//单点交叉函数
void cross(struct individual population[])
{
	int i,j,cn;//cn为交叉位置c1,c2为参与交叉的个体
	//int c1,c2;
	int index[SIZE],point,temp;
	float cm;//cm为交叉概率
	char t;
	////////////////////////////////
	//printf("====================\n");
    //  make a pair of individual randomly
	for(i=0;i<SIZE;i++)  index[i]=i;
	for(i=0;i<SIZE;i++)
	{
		point=rand()%(SIZE-i);
		temp=index[i];
		index[i]=index[point+i];
		index[point+i]=temp;
	}


    //one point crossover operation

	for(i=0;i<SIZE-1;i+=2)
	{
		cm=(float)(rand()%1000)/1000.0;//随机生成交叉算子

	    if(cm<=CRP)
		{
            cn=rand()%CSIZE;//printf("----------交叉点是%d\n",cn);
		    for(j=cn;j<CSIZE;j++)
			{
			  t=population[index[i]].chrom[j];
			  population[index[i]].chrom[j]=population[index[i+1]].chrom[j];
			  population[index[i+1]].chrom[j]=t;
			}
		}
	}
	return;
}
//============================================================================
//变异操作
void mutation(struct individual population[SIZE])
{
	int n,mn;double p; //mn为变异点
    for(n=0;n<SIZE;n++)
	{
		p=rand()%1000/1000.0; //随机生成变异算子
		if(p<MUP)
		{
			mn=rand()%CSIZE;//随机产生变异位
			population[n].chrom[mn]=((population[n].chrom[mn]=='1')?'0':'1');
		}
	}
}
//============================================================
int mmax(int a,int b)
{
  if(a>b)return a;
    else return b;
}
//============================================================================
int mmin(int a,int b)
{
  if(a<b)return a;
    else return b;
}
//==========================================================
/*原有的
//============================================================

int triangle(int a,int b,int c,int path[], int branch[],char s1[])
{
  int type,t;
  int i=0;//------------
  int x[14]={0};
for(t=0;t<PL;t++)
	  branch[t]=0;//分之距离初始化
  if(a>=b)
  {// 1
	    path[i++]=1;//---------------
		s[0]++;//--------------------路径节点1
	  t=a;a=b;b=t;
	   x[0]=1;
  }
  else
  {x[1]=1;branch[0]=a-b;}//分支距离1越大越好
  if(a>=c)
  {// 2
      path[i++]=2;//------------------
	  s[1]++;//-------------------------路径节点2
	  t=a;a=c;c=t;
	  x[2]=1;
  }
  else
	 { x[3]=1;branch[1]=a-c;}//分支距离2越大越好

  if(b>=c)
  {// 3
     path[i++]=3;//------------------
	 s[2]++;//----------------------路径节点3
     t=b;b=c;c=t;
	 x[4]=1;
  }
  else
  { x[5]=1;branch[2]=b-c;}//分支距离3越大越好

  if(a+b<=c)
  { // 4
     path[i++]=4;//------------------
     type=4;
	 branch[3]=a+b-c;//分支距离4越大越好
    x[6]=1;
  }
  else
  {   //5
	  x[7]=1;
	  path[i++]=5;//------------------
	 s[3]++;//----------------------路径节点4
	  if(a==b)
	  { //6
		  x[8]=1;
	     path[i++]=6;//------------------
	    s[4]++;//----------------------路径节点5

		  if(b==c)
		  {
			  x[9]=1;
			  type=1;// 7
  		      path[i++]=7;//------------------
			  s[5]++;//----------------------------路径节点6

              branch[5]=1;   //分支距离6
		  }
		  else
		  {
			   x[10]=1;
			  type=2;// 8
			  path[i++]=8;//------------------

		  }
	  }
	  else
	  {	 branch[4]=a-b;   //分支距离5越大越好

        x[11]=1;
		  if(b==c)
		  {
			  type=2;// 9
			  path[i++]=9;//------------------
			  x[12]=1;
		  }

		  else
		  {
			  type=3;// 10
			  path[i++]=10;//------------------
			  x[13]=1;
		  }

	  }
  }

  //将执行路径结果变成字符串存入数组中
for(i=0;i<14;i++)
{ //printf("%d",x[i]);
	if(x[i]==2) break;
	//if(x[i]==2) s[i]='2';
       if(x[i]==1) s1[i]='1';
        else  s1[i]='0';
}
s1[i]='\0';
//printf("\n");
//for(i=0;i<14;i++)
//printf("%c",s1[i]);
//printf("\n");
  return type;
}
 */
int triangle(int a,int b,int c,int path[], char s1[])
{
  int type,t;
  int i=0;//------------
  int x[14]={0};
//for(t=0;t<PL;t++)
	//  branch[t]=0;//分之距离初始化

  if(a>=b)
  {// 1
	    path[i++]=1;//---------------
		//s[0]++;//--------------------路径节点1
	  t=a;a=b;b=t;
	   x[0]=1;
  }
  else
  {x[1]=1;
  //branch[0]=a-b;
  }//分支距离1越大越好
  if(a>=c)
  {// 2
      path[i++]=2;//------------------
	 // s[1]++;//-------------------------路径节点2
	  t=a;a=c;c=t;
	  x[2]=1;
  }
  else
	 { x[3]=1;
  //branch[1]=a-c;
  }//分支距离2越大越好

  if(b>=c)
  {// 3
     path[i++]=3;//------------------
	// s[2]++;//----------------------路径节点3
     t=b;b=c;c=t;
	 x[4]=1;
  }
  else
  { x[5]=1;
  //branch[2]=b-c;
  }//分支距离3越大越好

  if(a+b<=c)
  { // 4
     path[i++]=4;//------------------
     type=4;
	// branch[3]=a+b-c;//分支距离4越大越好
    x[6]=1;
  }
  else
  {   //5
	  x[7]=1;
	  path[i++]=5;//------------------
	 //s[3]++;//----------------------路径节点4
	  if(a==b)
	  { //6
		  x[8]=1;
	     path[i++]=6;//------------------
	   // s[4]++;//----------------------路径节点5

		  if(b==c)
		  {
			  x[10]=1;
			  type=1;// 7
  		      path[i++]=7;//------------------
			 // s[5]++;//----------------------------路径节点6

             // branch[5]=1;   //分支距离6
		  }
		  else
		  {
			  x[11]=1;
			  type=2;// 8
			  path[i++]=8;//------------------

		  }
	  }
	  else
	  {	// branch[4]=a-b;   //分支距离5越大越好

        x[9]=1;
		  if(b==c)
		  {
			  type=2;// 9
			  path[i++]=9;//------------------
			  x[12]=1;
		  }

		  else
		  {
			  type=3;// 10
			  path[i++]=10;//------------------
			  x[13]=1;
		  }

	  }
  }
 /*printf("ljshi:");
  for(i=0;i<LEN;i++)
  {
	  if(path[i]==0)break;
	  printf("%d",path[i]);
  }
 printf("\n");

 */
  //将执行路径结果变成字符串存入数组中
for(i=0;i<14;i++)
{ //printf("%d",x[i]);
	if(x[i]==2) break;
	//if(x[i]==2) s[i]='2';
       if(x[i]==1) s1[i]='1';
        else  s1[i]='0';
}
s1[i]='\0';
//printf("\n");
//for(i=0;i<14;i++)
//printf("%c",s1[i]);
//printf("\n");
  return type;
}

//============================================================================
/*float fitness(int p[],int target[],int length,int branch[])//length路径长度，branch分支距离
{
  int i,j=0;
  for(i=0;i<length;i++)
         if(p[i]==target[i])j++;
		 else break;
 // printf("branch[%d]=%d\n",i,branch[i]);
  return j+1.0/(branch[i]+1);

}
*/
//============================================================================
float fitness0(int p[],int target[])//只计算节点穿越情况==本文新方法
{
	int i,j;
	float f=0;
/*
	for(j=0;j<PL;j++)
	   printf("%d次",s[j]);
	printf("\n");
*/
	for(i=0;i<LEN;i++)
	{
		if(p[i]==0)  break;
		for(j=0;j<PL;j++)//判断每目标路径中一个节点是否在穿越路径中出现
			if(p[i]==target[j])
			{   if(s[j]==0)  printf("****************!!!!!!!!!????????????\n");
				f+=1.0/s[j];
			}
	}
   // printf("f=%.3f  ",f);
	return f+1;
}
//============================================================================
float fitness(int p[],int target[],int branch[])//新方法权重*（传统方法 层接近度+分支距离）
{
	int i,j,cjjd,c,fzjl;
	float f=0,fq=0;

	//for(j=0;j<PL;j++)
	//   printf("%d次",s[j]);
	//printf("\n");

	//先计算权重
	for(i=0;i<LEN;i++)
	{
		if(p[i]==0)  break;
		for(j=0;j<PL;j++)//判断每目标路径中一个节点是否在穿越路径中出现
			if(p[i]==target[j])
			{   if(s[j]==0)  printf("****************!!!!!!!!!????????????\n");
				fq+=1.0/s[j];
			}
	}
  //    printf("fq=%.3f  ",fq);




	cjjd=0;fzjl=0;
	for(i=0;i<PL;i++)//统计从前向后相同节点个数
		if(p[i]!=target[i])break;
             else cjjd++;
	for(i=0;i<PL;i++)//判断每目标路径中一个节点是否在穿越路径中出现
	{  c=0;//标记是否穿越目标路径的第i个节点
       for(j=0;j<LEN;j++)
	   {
		   if(p[j]==0)  break;
		   if(p[j]==target[i]){c=1;break;}
	   }
	   if(c!=1)
	   {
		   //printf("branch[%d]=%d ",i,branch[i]);
		   fzjl+=branch[i];
	   }
	}
	fzjl=4*SIZE1+fzjl;//将分支距离处理成正的值
 //    printf("cjjd=%d ,fzjl=%d ",cjjd,fzjl);
 	f=fq*(1.0*cjjd/PL+1-pow(1.001,-fzjl));
     //	f=1.0*cjjd/PL+1-pow(1.001,-fzjl);
//	f=1.0*cjjd/PL ;
 //  printf("f=%.3f  ",f);
	return f ;
}
//============================================================================
main()
{
  struct
	  {
       	int dd;
		//int lj[LEN];
        //float ff;
	    int data[N];
	  }find[NUM];
       //find存放找到的数据，dd代表找到的代数、lj代表数据穿越的的路径、ff代表适应值；
          //data[N]代表测试数据，NUM代表试验次数

  clock_t start, finish; //=================time
  double duration; //=================time
  int count=0;
  struct individual population[SIZE];
 // int rr[NUM];  //存放每次实验找到最优解的进化代数
    char t1[SIZE][6]; /* 存放一次执行路径 *///@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  int t[LEN]={1,2,3,5,6,7};//目标路径
  int i,j,k;
  int ffind[NUM];//标记是否找到目标路径
  int y,ty=0; //统计实验次数
  char sz[SIZE][14]={0};//用以存放数据穿越的节点情况数组，0，1
  int sum1[14]={0};//用以存放所有数据运行程序后程序中每个节点被穿越情况=程序均衡性
  int del[SIZE][14]={0};//用以存放除某个数据外其他数据运行程序后程序中每个节点被穿越情况
 // int [SIZE]={0};//
  float balance=0,sbalance[SIZE]={0};//程序均衡性
  float eff[SIZE]={0};
 //影响因子
  float ba[7]={0};
  float ba1[SIZE][7]={0};
	float cgl=0;
	double pjcs=0;

  printf("=========================\n");

  if((fp=fopen(" 最大代数60000 传统方法作为适应值 20120512  三角形分类1-2048 种群200 　 路径2  5次_1.txt","w"))==NULL)   //打开一个文本文件用于存放生成的测试数据
        {
      	  printf("not pen file! fp\n");
      	  exit(0);
        }
   /* if((fpd=fopen("初始种群　20111103  三角形分类1-1024种群200 　 路径2  15次_2.txt","r"))==NULL)   //打开一个文本文件用于存放生成的测试数据
        {
      	  printf("not pen file! fpd\n");
      	  exit(0);
        }*/
srand((unsigned)time(NULL));//先使用随机数”种子”初始化
  start = clock(); //=================time
 //输出目标路径
  printf("目标路径：");
  fprintf(fp,"目标路径：");
  for(i=0;i<LEN;i++)
  { printf("%3d",t[i]);fprintf(fp,"%3d",t[i]);}
    printf("\n");
    fprintf(fp,"\n");

  for(y=0;y<NUM;y++)//一共进行NUM次试验
  {
      count=0; ffind[y]=0;
   //随机产生初始种群＋初始化个体穿越路径
     for(k=0;k<SIZE;k++)/* 随机生成SIZE组数据*/
		  {
			 for(i=0;i<N;i++)
			 { // fscanf(fpd,"%d",&population[k].p[i]);
			  // printf("%d  ",population[k].p[i]);
              population[k].p[i]=rand()%SIZE1+1;
               //fprintf(fpd,"%5d",population[k].p[i]);

			 //  fprintf(fpd,"%5d",population[k].p[i]);
			 }
                     // fprintf(fpd,"\n");
		  }//end of/* 随机生成SIZE组数据*/
      while(count<G)  //while-1 ***********
	  {
		  count++;  //进化代数加1
		  //----------------------------------------------------------------------------------------
         // for(k=0;k<PL;k++)
			//  s[k]=0;      //给s数组赋初值
		  //------------------------------------------------------------------------------------------
         balance=0;
          for(k=0;k<SIZE;k++)//for-1======
		  {   //将数据带入ymd程序执行

			 sbalance[k]=0;
			 eff[k]=0;
			 for(i=0;i<14;i++)
			 {del[k][i]=0;sum1[i]=0;}
			  for(i=0;i<7;i++)
              {ba[i]=0;ba1[k][i]=0;}


			  for(i=0;i<LEN;i++)
				  population[k].tpath[i]=0; //先给tpath赋初值

              triangle(population[k].p[0],population[k].p[1], population[k].p[2],population[k].tpath,t1[k]);
             //printf("---  执行路径是： %s  ----\n",t1[k]);

			  for(i=0;i<14;i++)
			  { sz[k][i]=t1[k][i];
//			   printf("%c ",sz[k][i]);
//               fprintf(fp,"%c  ",sz[k][i]);
			  }

//                printf("\n");
                //fprintf(fp,"\n");

              //判断是否找到目标路径　　
              for(i=0;i<LEN;i++)
			      if(t[i]!=population[k].tpath[i])break;
              if(i==LEN)//找到目标路径
			  {
				  ffind[y]=1;
//				 printf("y=%d,count=%d  \n",y,count);

				  find[y].data[0]=population[k].p[0];
				  find[y].data[1]=population[k].p[1];
				  find[y].data[2]=population[k].p[2];
				  break;
			  }
		  }// end of for -1=========

		  /*	for(j=0;j<PL;j++)
                 	   printf("%d次",s[j]);
            	printf("\n");
				*/

		  if(ffind[y]==1) break;
		   //计算程序均衡性及节点均衡性
              for(i=0;i<14;i++)
			  {      for(k=0;k<SIZE;k++)
						if(sz[k][i]!='0')
						sum1[i]=sum1[i]+1;//穿越各节点的流量
					  //printf("%d,",sum1[i]);
			  }
			  //printf("\n");
			  ba[0]=(sum1[0]+sum1[1]==0)?0:fabs(sum1[0]-sum1[1])/(sum1[0]+sum1[1]);
              ba[1]=(sum1[2]+sum1[3]==0)?0:fabs(sum1[2]-sum1[3])/(sum1[2]+sum1[3]);
			  ba[2]=(sum1[4]+sum1[5]==0)?0:fabs(sum1[4]-sum1[5])/(sum1[4]+sum1[5]);
			  ba[3]=(sum1[6]+sum1[7]==0)?0:fabs(sum1[6]-sum1[7])/(sum1[6]+sum1[7]);
			  ba[4]=(sum1[8]+sum1[11]==0)?0:fabs(sum1[8]-sum1[11])/(sum1[8]+sum1[11]);
			  ba[5]=(sum1[9]+sum1[10]==0)?0:fabs(sum1[9]-sum1[10])/(sum1[9]+sum1[10]);
			  ba[6]=(sum1[12]+sum1[13]==0)?0:fabs(sum1[12]-sum1[13])/(sum1[12]+sum1[13]);
			  balance=(ba[0]+ba[1]+ba[2]+ba[3]+ba[4]+ba[5]+ba[6])/7;
			  /*
		     for(k=0;k<7;k++)
              printf("%f,",ba[k]);
		      */
			 // printf("balance=%f,",balance);

//求删除各数据穿越情况后的均衡性
			  //printf("del:\n");
          for(j=0;j<SIZE;j++)
		  {
			  for(i=0;i<14;i++)
			  {      for(k=0;k<SIZE;k++)
						if(sz[k][i]!='0'&&j!=k)
						del[j][i]=del[j][i]+1;//穿越各节点的流量
					  //printf("%d,",del[j][i]);
				}
        //printf("\n");

		  }
		  //计算删除各数据的影响
		   for(j=0;j<SIZE;j++)
		   {  ba1[j][0]=(del[j][0]+del[j][1]==0)?0:fabs(del[j][0]-del[j][1])/(del[j][0]+del[j][1]);
              ba1[j][1]=(del[j][2]+del[j][3]==0)?0:fabs(del[j][2]-del[j][3])/(del[j][2]+del[j][3]);
			  ba1[j][2]=(del[j][4]+del[j][5]==0)?0:fabs(del[j][4]-del[j][5])/(del[j][4]+del[j][5]);
			  ba1[j][3]=(del[j][6]+del[j][7]==0)?0:fabs(del[j][6]-del[j][7])/(del[j][6]+del[j][7]);
			  ba1[j][4]=(del[j][8]+del[j][11]==0)?0:fabs(del[j][8]-del[j][11])/(del[j][8]+del[j][11]);
			  ba1[j][5]=(del[j][9]+del[j][10]==0)?0:fabs(del[j][9]-del[j][10])/(del[j][9]+del[j][10]);
			  ba1[j][6]=(del[j][12]+del[j][13]==0)?0:fabs(del[j][12]-del[j][13])/(del[j][12]+del[j][13]);
			  sbalance[j]=(ba1[j][0]+ba1[j][1]+ba1[j][2]+ba1[j][3]+ba1[j][4]+ba1[j][5]+ba1[j][6])/7;
		   }
          for(j=0;j<SIZE;j++)
		  {
		  //for(k=0;k<7;k++)
             // printf("%f,",ba1[j][k]);
		      // printf("sbalance[j]=%f,",sbalance[j]);

		  if(balance>=sbalance[j])
			   eff[j]=0.0001;
		  else
			  eff[j]=fabs(balance-sbalance[j]);
		  // printf("eff=%f,",eff[j]);

		  }
			  //printf("\n");

          //计算适应值
          for(k=0;k<SIZE;k++)
		  {
			  population[k].fitness=eff[k];

		  }
          selection(population);//选择
          //编码
          for(k=0;k<SIZE;k++)
                code(population[k].p,population[k].chrom); //编码
	       //对上一代种群进行交叉变异生成下一代种群
           cross(population);//单点交叉
           mutation(population);//变异操作
        	//解码
           for(k=0;k<SIZE;k++)
                 uncode(population[k].p,population[k].chrom); //解码
	  }//end of while-1 ********
	  ty+=count;
      find[y].dd=count;
  }// end of for(y=0;y<NUM;y++)//一共进行NUM次试验

  //输出时间
  finish = clock(); //=================time
  duration = (double)(finish - start) / CLOCKS_PER_SEC; //===================time
  printf( "average time:%f seconds\n", duration/NUM ); //===================time
  fprintf(fp, "average time:%f seconds\n", duration/NUM ); //===================time
  //-----------------------------------------------------------------
  //输出结果信息
  printf("&&&&&&&&&&&&&  %3d次试验， 平均运行代数为 %f  &&&&&&&&&&\n",NUM,1.0*ty/NUM);
  fprintf(fp,"&&&&&&&&&&&&&   %3d次试验， 平均运行代数为   %f  &&&&&&&&&&\n",NUM,1.0*ty/NUM);
  printf("&&&&&&&&&&&&&  每次找到目标的进化代数如下  &&&&&&&&&&\n");
  fprintf(fp,"&&&&&&&&&&&&&  每次找到目标的进化代数如下  &&&&&&&&&&\n");
  //printf("pjcs=%f\n",pjcs);
  for(i=0;i<NUM;i++)
  { pjcs=pjcs+find[i].dd*SIZE;
  //printf("find[i].dd=%d\n",find[i].dd);
  //printf("pjcs=%f\n",pjcs);
	  if(ffind[i]==1)
	  {  cgl=cgl+1;

		  printf("======================================\n\n========第%d次 %d代 =============\n",i+1,find[i].dd);
         fprintf(fp,"=========================================\n\n========第%d次 %d代 ==============\n",i+1,find[i].dd);
		 printf("----------------------数据是-----------------------\n");
         fprintf(fp,"------------------------数据是-- ---------------------------\n");
		  for(j=0;j<N;j++)
		  {
               printf("%6d",find[i].data[j]);
               fprintf(fp,"%6d",find[i].data[j]);

		  }
	      printf("\n");fprintf(fp,"\n");
	  }//end of if(ffind[i]==1)
	  else
	  {
		  printf("======================================\n\n========第%d次实验,  没找到测试数据\n",i+1);
		  fprintf(fp,"======================================\n\n========第%d次实验,  没找到测试数据\n",i+1);
	  }

  }
   printf("平均评价次数是：%f\n",pjcs/NUM);
  printf("成功率是：%f",cgl/NUM);
  fprintf(fp,"=======================\n种群规模为%d,测试数据范围大小%d,进化代数%d,实验次数%d\n", SIZE,SIZE1,G,NUM);
  fclose(fp);
//fclose(fpd);
  getch();

}


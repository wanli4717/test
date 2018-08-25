
//��ͳ����������Ӧֵ��֧����Ͳ�ӽ���//���������۴����ͳɹ���

//Ŀ��·��·��2���ȱ������� 1-2-3-5-6-7

#include <string.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#define SIZE 50//������Ⱥ��ģ

#define SIZE1 128 //���Ʋ������ݷ�Χ��С0-7000
#define NUM  100//ʵ��������������ȡƽ��ֵ

#define CSIZE  33 //!!!!!!!!!!!!!!!!!!!����3���� ��Ҫ�Ķ����ƴ洢λ��
#define G 100  //���ƽ�������
#define CRP  0.9//��������
#define MUP  0.3//��������

#define PL 6//Ŀ��·���Ľڵ���

#define N 3 // N  �������鳤��  year month day 3����
//#define PN  3 //PN����Ŀ��·������@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//#define CN  28 //CN ����·�����볤�� CN=N-1  + N-2 + ���� + 1

//#define FITNESS 0.19999//��¼��Ӧֵ����FITNESS�ľ���������ҵ�5��ʱ�����˳�����
//#define DIF  1.0  //Ҫ�����ٿ�������ͬλ���İٷֱ�

FILE *fp ;
#include "stdlib.h"
#include "time.h"

#define LEN  10//Ŀ��·������

//-----------------------------------------------------------------------------------
int s[PL];//ͳ��Ŀ��·����4���ڵ㱻ִ�еĴ���,��ÿһ����Ⱥִ��֮ǰ��0��ÿ����Ⱥ���н�����ͳ����Ӧֵ

//-----------------------------------------------------------------------------------


struct individual
{   char chrom[CSIZE+1];//Ⱦɫ�����
	int p[N];    //��������
	int tpath[LEN];//��Խ·��
	int branch[PL];//��ŷ�֧����
	double fitness;//��Ӧֵ����
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

//���뺯��
void code(int n[],char c[])
{
	int t,l,i,j;
	for(l=0;l<CSIZE;)
	{
		i=l/(CSIZE/N);t=n[i];j=0;//jΪ��������ȷ��ÿ�����ݶ�Ӧ7λ2������
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
//���뺯��
void uncode(int n[],char c[])
{
	int t,l,i,j,m;
	for(l=0;l<CSIZE;)
	{
		i=l/(CSIZE/N);t=0;m=1;j=0;//jΪ��������ȷ��ÿ�����ݶ�Ӧ7λ2������
	    while(j<(CSIZE/N))
		{
			t=t+(c[l++]-'0')*m;
			m*=2;j++;
		}
	    n[i]=t%(SIZE1+1);
	}
}
//===========================================================================
//���㽻�溯��
void cross(struct individual population[])
{
	int i,j,cn;//cnΪ����λ��c1,c2Ϊ���뽻��ĸ���
	//int c1,c2;
	int index[SIZE],point,temp;
	float cm;//cmΪ�������
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
		cm=(float)(rand()%1000)/1000.0;//������ɽ�������

	    if(cm<=CRP)
		{
            cn=rand()%CSIZE;//printf("----------�������%d\n",cn);
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
//�������
void mutation(struct individual population[SIZE])
{
	int n,mn;double p; //mnΪ�����
    for(n=0;n<SIZE;n++)
	{
		p=rand()%1000/1000.0; //������ɱ�������
		if(p<MUP)
		{
			mn=rand()%CSIZE;//�����������λ
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
/*ԭ�е�
//============================================================

int triangle(int a,int b,int c,int path[], int branch[],char s1[])
{
  int type,t;
  int i=0;//------------
  int x[14]={0};
for(t=0;t<PL;t++)
	  branch[t]=0;//��֮�����ʼ��
  if(a>=b)
  {// 1
	    path[i++]=1;//---------------
		s[0]++;//--------------------·���ڵ�1
	  t=a;a=b;b=t;
	   x[0]=1;
  }
  else
  {x[1]=1;branch[0]=a-b;}//��֧����1Խ��Խ��
  if(a>=c)
  {// 2
      path[i++]=2;//------------------
	  s[1]++;//-------------------------·���ڵ�2
	  t=a;a=c;c=t;
	  x[2]=1;
  }
  else
	 { x[3]=1;branch[1]=a-c;}//��֧����2Խ��Խ��

  if(b>=c)
  {// 3
     path[i++]=3;//------------------
	 s[2]++;//----------------------·���ڵ�3
     t=b;b=c;c=t;
	 x[4]=1;
  }
  else
  { x[5]=1;branch[2]=b-c;}//��֧����3Խ��Խ��

  if(a+b<=c)
  { // 4
     path[i++]=4;//------------------
     type=4;
	 branch[3]=a+b-c;//��֧����4Խ��Խ��
    x[6]=1;
  }
  else
  {   //5
	  x[7]=1;
	  path[i++]=5;//------------------
	 s[3]++;//----------------------·���ڵ�4
	  if(a==b)
	  { //6
		  x[8]=1;
	     path[i++]=6;//------------------
	    s[4]++;//----------------------·���ڵ�5

		  if(b==c)
		  {
			  x[9]=1;
			  type=1;// 7
  		      path[i++]=7;//------------------
			  s[5]++;//----------------------------·���ڵ�6

              branch[5]=1;   //��֧����6
		  }
		  else
		  {
			   x[10]=1;
			  type=2;// 8
			  path[i++]=8;//------------------

		  }
	  }
	  else
	  {	 branch[4]=a-b;   //��֧����5Խ��Խ��

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

  //��ִ��·���������ַ�������������
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
	//  branch[t]=0;//��֮�����ʼ��

  if(a>=b)
  {// 1
	    path[i++]=1;//---------------
		//s[0]++;//--------------------·���ڵ�1
	  t=a;a=b;b=t;
	   x[0]=1;
  }
  else
  {x[1]=1;
  //branch[0]=a-b;
  }//��֧����1Խ��Խ��
  if(a>=c)
  {// 2
      path[i++]=2;//------------------
	 // s[1]++;//-------------------------·���ڵ�2
	  t=a;a=c;c=t;
	  x[2]=1;
  }
  else
	 { x[3]=1;
  //branch[1]=a-c;
  }//��֧����2Խ��Խ��

  if(b>=c)
  {// 3
     path[i++]=3;//------------------
	// s[2]++;//----------------------·���ڵ�3
     t=b;b=c;c=t;
	 x[4]=1;
  }
  else
  { x[5]=1;
  //branch[2]=b-c;
  }//��֧����3Խ��Խ��

  if(a+b<=c)
  { // 4
     path[i++]=4;//------------------
     type=4;
	// branch[3]=a+b-c;//��֧����4Խ��Խ��
    x[6]=1;
  }
  else
  {   //5
	  x[7]=1;
	  path[i++]=5;//------------------
	 //s[3]++;//----------------------·���ڵ�4
	  if(a==b)
	  { //6
		  x[8]=1;
	     path[i++]=6;//------------------
	   // s[4]++;//----------------------·���ڵ�5

		  if(b==c)
		  {
			  x[10]=1;
			  type=1;// 7
  		      path[i++]=7;//------------------
			 // s[5]++;//----------------------------·���ڵ�6

             // branch[5]=1;   //��֧����6
		  }
		  else
		  {
			  x[11]=1;
			  type=2;// 8
			  path[i++]=8;//------------------

		  }
	  }
	  else
	  {	// branch[4]=a-b;   //��֧����5Խ��Խ��

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
  //��ִ��·���������ַ�������������
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
/*float fitness(int p[],int target[],int length,int branch[])//length·�����ȣ�branch��֧����
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
float fitness0(int p[],int target[])//ֻ����ڵ㴩Խ���==�����·���
{
	int i,j;
	float f=0;
/*
	for(j=0;j<PL;j++)
	   printf("%d��",s[j]);
	printf("\n");
*/
	for(i=0;i<LEN;i++)
	{
		if(p[i]==0)  break;
		for(j=0;j<PL;j++)//�ж�ÿĿ��·����һ���ڵ��Ƿ��ڴ�Խ·���г���
			if(p[i]==target[j])
			{   if(s[j]==0)  printf("****************!!!!!!!!!????????????\n");
				f+=1.0/s[j];
			}
	}
   // printf("f=%.3f  ",f);
	return f+1;
}
//============================================================================
float fitness(int p[],int target[],int branch[])//�·���Ȩ��*����ͳ���� ��ӽ���+��֧���룩
{
	int i,j,cjjd,c,fzjl;
	float f=0,fq=0;

	//for(j=0;j<PL;j++)
	//   printf("%d��",s[j]);
	//printf("\n");

	//�ȼ���Ȩ��
	for(i=0;i<LEN;i++)
	{
		if(p[i]==0)  break;
		for(j=0;j<PL;j++)//�ж�ÿĿ��·����һ���ڵ��Ƿ��ڴ�Խ·���г���
			if(p[i]==target[j])
			{   if(s[j]==0)  printf("****************!!!!!!!!!????????????\n");
				fq+=1.0/s[j];
			}
	}
  //    printf("fq=%.3f  ",fq);




	cjjd=0;fzjl=0;
	for(i=0;i<PL;i++)//ͳ�ƴ�ǰ�����ͬ�ڵ����
		if(p[i]!=target[i])break;
             else cjjd++;
	for(i=0;i<PL;i++)//�ж�ÿĿ��·����һ���ڵ��Ƿ��ڴ�Խ·���г���
	{  c=0;//����Ƿ�ԽĿ��·���ĵ�i���ڵ�
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
	fzjl=4*SIZE1+fzjl;//����֧���봦�������ֵ
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
       //find����ҵ������ݣ�dd�����ҵ��Ĵ�����lj�������ݴ�Խ�ĵ�·����ff������Ӧֵ��
          //data[N]����������ݣ�NUM�����������

  clock_t start, finish; //=================time
  double duration; //=================time
  int count=0;
  struct individual population[SIZE];
 // int rr[NUM];  //���ÿ��ʵ���ҵ����Ž�Ľ�������
    char t1[SIZE][6]; /* ���һ��ִ��·�� *///@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  int t[LEN]={1,2,3,5,6,7};//Ŀ��·��
  int i,j,k;
  int ffind[NUM];//����Ƿ��ҵ�Ŀ��·��
  int y,ty=0; //ͳ��ʵ�����
  char sz[SIZE][14]={0};//���Դ�����ݴ�Խ�Ľڵ�������飬0��1
  int sum1[14]={0};//���Դ�������������г���������ÿ���ڵ㱻��Խ���=���������
  int del[SIZE][14]={0};//���Դ�ų�ĳ�������������������г���������ÿ���ڵ㱻��Խ���
 // int [SIZE]={0};//
  float balance=0,sbalance[SIZE]={0};//���������
  float eff[SIZE]={0};
 //Ӱ������
  float ba[7]={0};
  float ba1[SIZE][7]={0};
	float cgl=0;
	double pjcs=0;

  printf("=========================\n");

  if((fp=fopen(" ������60000 ��ͳ������Ϊ��Ӧֵ 20120512  �����η���1-2048 ��Ⱥ200 �� ·��2  5��_1.txt","w"))==NULL)   //��һ���ı��ļ����ڴ�����ɵĲ�������
        {
      	  printf("not pen file! fp\n");
      	  exit(0);
        }
   /* if((fpd=fopen("��ʼ��Ⱥ��20111103  �����η���1-1024��Ⱥ200 �� ·��2  15��_2.txt","r"))==NULL)   //��һ���ı��ļ����ڴ�����ɵĲ�������
        {
      	  printf("not pen file! fpd\n");
      	  exit(0);
        }*/
srand((unsigned)time(NULL));//��ʹ������������ӡ���ʼ��
  start = clock(); //=================time
 //���Ŀ��·��
  printf("Ŀ��·����");
  fprintf(fp,"Ŀ��·����");
  for(i=0;i<LEN;i++)
  { printf("%3d",t[i]);fprintf(fp,"%3d",t[i]);}
    printf("\n");
    fprintf(fp,"\n");

  for(y=0;y<NUM;y++)//һ������NUM������
  {
      count=0; ffind[y]=0;
   //���������ʼ��Ⱥ����ʼ�����崩Խ·��
     for(k=0;k<SIZE;k++)/* �������SIZE������*/
		  {
			 for(i=0;i<N;i++)
			 { // fscanf(fpd,"%d",&population[k].p[i]);
			  // printf("%d  ",population[k].p[i]);
              population[k].p[i]=rand()%SIZE1+1;
               //fprintf(fpd,"%5d",population[k].p[i]);

			 //  fprintf(fpd,"%5d",population[k].p[i]);
			 }
                     // fprintf(fpd,"\n");
		  }//end of/* �������SIZE������*/
      while(count<G)  //while-1 ***********
	  {
		  count++;  //����������1
		  //----------------------------------------------------------------------------------------
         // for(k=0;k<PL;k++)
			//  s[k]=0;      //��s���鸳��ֵ
		  //------------------------------------------------------------------------------------------
         balance=0;
          for(k=0;k<SIZE;k++)//for-1======
		  {   //�����ݴ���ymd����ִ��

			 sbalance[k]=0;
			 eff[k]=0;
			 for(i=0;i<14;i++)
			 {del[k][i]=0;sum1[i]=0;}
			  for(i=0;i<7;i++)
              {ba[i]=0;ba1[k][i]=0;}


			  for(i=0;i<LEN;i++)
				  population[k].tpath[i]=0; //�ȸ�tpath����ֵ

              triangle(population[k].p[0],population[k].p[1], population[k].p[2],population[k].tpath,t1[k]);
             //printf("---  ִ��·���ǣ� %s  ----\n",t1[k]);

			  for(i=0;i<14;i++)
			  { sz[k][i]=t1[k][i];
//			   printf("%c ",sz[k][i]);
//               fprintf(fp,"%c  ",sz[k][i]);
			  }

//                printf("\n");
                //fprintf(fp,"\n");

              //�ж��Ƿ��ҵ�Ŀ��·������
              for(i=0;i<LEN;i++)
			      if(t[i]!=population[k].tpath[i])break;
              if(i==LEN)//�ҵ�Ŀ��·��
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
                 	   printf("%d��",s[j]);
            	printf("\n");
				*/

		  if(ffind[y]==1) break;
		   //�����������Լ��ڵ������
              for(i=0;i<14;i++)
			  {      for(k=0;k<SIZE;k++)
						if(sz[k][i]!='0')
						sum1[i]=sum1[i]+1;//��Խ���ڵ������
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

//��ɾ�������ݴ�Խ�����ľ�����
			  //printf("del:\n");
          for(j=0;j<SIZE;j++)
		  {
			  for(i=0;i<14;i++)
			  {      for(k=0;k<SIZE;k++)
						if(sz[k][i]!='0'&&j!=k)
						del[j][i]=del[j][i]+1;//��Խ���ڵ������
					  //printf("%d,",del[j][i]);
				}
        //printf("\n");

		  }
		  //����ɾ�������ݵ�Ӱ��
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

          //������Ӧֵ
          for(k=0;k<SIZE;k++)
		  {
			  population[k].fitness=eff[k];

		  }
          selection(population);//ѡ��
          //����
          for(k=0;k<SIZE;k++)
                code(population[k].p,population[k].chrom); //����
	       //����һ����Ⱥ���н������������һ����Ⱥ
           cross(population);//���㽻��
           mutation(population);//�������
        	//����
           for(k=0;k<SIZE;k++)
                 uncode(population[k].p,population[k].chrom); //����
	  }//end of while-1 ********
	  ty+=count;
      find[y].dd=count;
  }// end of for(y=0;y<NUM;y++)//һ������NUM������

  //���ʱ��
  finish = clock(); //=================time
  duration = (double)(finish - start) / CLOCKS_PER_SEC; //===================time
  printf( "average time:%f seconds\n", duration/NUM ); //===================time
  fprintf(fp, "average time:%f seconds\n", duration/NUM ); //===================time
  //-----------------------------------------------------------------
  //��������Ϣ
  printf("&&&&&&&&&&&&&  %3d�����飬 ƽ�����д���Ϊ %f  &&&&&&&&&&\n",NUM,1.0*ty/NUM);
  fprintf(fp,"&&&&&&&&&&&&&   %3d�����飬 ƽ�����д���Ϊ   %f  &&&&&&&&&&\n",NUM,1.0*ty/NUM);
  printf("&&&&&&&&&&&&&  ÿ���ҵ�Ŀ��Ľ�����������  &&&&&&&&&&\n");
  fprintf(fp,"&&&&&&&&&&&&&  ÿ���ҵ�Ŀ��Ľ�����������  &&&&&&&&&&\n");
  //printf("pjcs=%f\n",pjcs);
  for(i=0;i<NUM;i++)
  { pjcs=pjcs+find[i].dd*SIZE;
  //printf("find[i].dd=%d\n",find[i].dd);
  //printf("pjcs=%f\n",pjcs);
	  if(ffind[i]==1)
	  {  cgl=cgl+1;

		  printf("======================================\n\n========��%d�� %d�� =============\n",i+1,find[i].dd);
         fprintf(fp,"=========================================\n\n========��%d�� %d�� ==============\n",i+1,find[i].dd);
		 printf("----------------------������-----------------------\n");
         fprintf(fp,"------------------------������-- ---------------------------\n");
		  for(j=0;j<N;j++)
		  {
               printf("%6d",find[i].data[j]);
               fprintf(fp,"%6d",find[i].data[j]);

		  }
	      printf("\n");fprintf(fp,"\n");
	  }//end of if(ffind[i]==1)
	  else
	  {
		  printf("======================================\n\n========��%d��ʵ��,  û�ҵ���������\n",i+1);
		  fprintf(fp,"======================================\n\n========��%d��ʵ��,  û�ҵ���������\n",i+1);
	  }

  }
   printf("ƽ�����۴����ǣ�%f\n",pjcs/NUM);
  printf("�ɹ����ǣ�%f",cgl/NUM);
  fprintf(fp,"=======================\n��Ⱥ��ģΪ%d,�������ݷ�Χ��С%d,��������%d,ʵ�����%d\n", SIZE,SIZE1,G,NUM);
  fclose(fp);
//fclose(fpd);
  getch();

}


/*
align the qulity score
input: read.fq--there should be a qulity score
output: RL code qulity score ---- 8 records or 16 records
// output:  align result or statistic result,such as starray[][]-----should include the file "statisticChar.c"
score range: 33--126
*/
#include<stdlib.h>
#include<stdio.h>

//#define DEBUG
#define RECORDNUM  8
#define CLUSTER_RATE  0.4


#define BITSET(a,b) ( a|=(1<<b) )  //set 1
#define BITCLK(a,b) ( (a)&=~(1<<(b)) )  //set 0
#define BITTEST(a,b) ( (a)&(1<<(b)) )  //test 1 or 0
//public variate---------------------
int ReadLength, threshold_value, sameCount[RECORDNUM][RECORDNUM]={0}, a[RECORDNUM]={0}, b[RECORDNUM]={0}, **tranWord;
unsigned int *c, c_size;
//func-------------------------------
int CodeQuality(char *buf, char *str, int n) //return the length of buf
{
	//buf store RLcode, str is input seq, n is RL size
	char prech, ch;
	int  i, count=1, len=0;  //len means the length of buf
	n--;
	if( (prech=*str)<33 )
		return len;
	for(i=1; (ch=*(str+i))>32; i++)
	{
		if(prech==ch)
			count++;		
		else
		{
			if(count>n) //not >=	256
			{	
				if( count<256 )
					len+=sprintf(buf+len,"%c%c",-prech,count);
				else //up to 256
				{
					do
					{
						len+=sprintf(buf+len,"%c%c",-prech,255);
						count-=255;
					}while(count>255);
					if(count>n)
						len+=sprintf(buf+len,"%c%c",-prech,count);
					else
						while( count-- >0 )
							*(buf+len++)=prech;
				}
			}
			else			
				while( count-- >0 )
					*(buf+len++)=prech;
					//sprintf(buf+len,"%c",prech+OFFSET);
			
			count=1;
			prech=ch;
		}
	}
	if(count>n) //not >=			
		len+=sprintf(buf+len,"%c%c",-prech,count);			
	else			
		while( count-- >0 )
			*(buf+len++)=prech;
	
	//*(buf+len++)=10; //endline
	return len;
}

int Code2Quality(char *buf, char *str, int n) //mv code
{
	char prech, ch, pprech=0;
	int  i, count=1, len=0, countch=1;  //len means the length of buf
	n--;
	if( (prech=*str)<33 )
		return len;
	for(i=1; (ch=*(str+i))>32; )
	{
		if(prech==ch)
			count++;		
		else
		{	
			if(ch==pprech)
			{
				if(count<33)  //if count>32 the do next if
				{					
					i++;
					while( ch==*(str+i++) ) countch++;
					if(countch<256)
						len+=sprintf(buf+len,"%c%c",count,countch);
					else //countch>=256
					{
						len+=sprintf(buf+len,"%c%c",count,255);
						countch-=255;
						
						while(countch>256)
						{
							len+=sprintf(buf+len,"%c%c",-pprech,255);
							countch-=255;
						}
						if(countch>n)
						len+=sprintf(buf+len,"%c%c",-pprech,countch);
						else
							while( count-- >0 )
								*(buf+len++)=pprech;
					}
					pprech=0; countch=1;
					continue;
				}				
			}
			if(count>n) //not >=	256
			{	
				if( count<256 )
					len+=sprintf(buf+len,"%c%c",-prech,count);
				else //up to 256
				{
					do
					{
						len+=sprintf(buf+len,"%c%c",-prech,255);
						count-=255;
					}while(count>255);
					if(count>n)
						len+=sprintf(buf+len,"%c%c",-prech,count);
					else
						while( count-- >0 )
							*(buf+len++)=prech;
				}
				pprech=prech; //record last code char
			}
			else
			{
				while( count-- >0 )
					*(buf+len++)=prech;
					//sprintf(buf+len,"%c",prech+OFFSET);
				pprech=0; //fresh pprech
			}
			count=1;
			prech=ch;
		}
		i++;
	}
	if(count>n) //not >=			
		len+=sprintf(buf+len,"%c%c",-prech,count);			
	else			
		while( count-- >0 )
			*(buf+len++)=prech;
	
	//*(buf+len++)=10; //endline
	return len;
}
int main(int argc, char*argv[])
{
	int i,j, p, cluRecord=0, cluPara=1, readCount=0, codeCount=0, codeLen=0, temp[RECORDNUM]={0},t;
	
	FILE *readFile, *outFile;
	char str[4096]={0}, strtemp[4096]={0}, *outputbuf, *codebuf;

#ifdef DEBUG
	readFile=fopen("SRR1063349.fastq","rb");
	outFile =fopen("out.txt","wb");
#else
	if(argc<2)
	{ printf("errer parameter !\ncmd read.fastq outfile\n"); return -1; }

	readFile=fopen(argv[1],"rb");
	outFile =fopen(argv[2],"wb");
#endif

	if(!readFile || !outFile)
	{
		printf("error opening file !\n");
		return -1;
	}
	/*
	//init
	cluPara=1;
	for(i=0; i<RECORDNUM; i++)
	{
		str[i]=(char*)malloc(4096*sizeof(char));
		cluPara*=2;
	}
	outputbuf=(char*)malloc(4096*RECORDNUM*sizeof(char)); //enough size
	codebuf  =(char*)malloc(4096*RECORDNUM*sizeof(char)); //enough size
	tranWord=(int**)malloc(cluPara*sizeof(int *));
	for(i=0; i<cluPara; i++)
	{
		*(tranWord+i)=(int*)calloc(RECORDNUM+1,sizeof(int));
		p=1; t=0;
		for(j=0; j<RECORDNUM; j++)
		{
			if( BITTEST(i,j) )
			{
				tranWord[i][p]=j;
				p++;
			}
			else 
				temp[t++]=j;
		}
		tranWord[i][0]=p-1;    //the number of 1 in i
		for(j=0; j<t; j++)
			tranWord[i][j+p]=temp[j];
	}

	fgets(strtemp,4096,readFile);
	fgets(strtemp,4096,readFile);
	for(i=0; strtemp[i]>64; i++)
		;
	ReadLength=i; //just for fixed length
	threshold_value=ReadLength*CLUSTER_RATE;
	c_size=i/32+1;
	c=(unsigned int*)calloc(c_size,sizeof(unsigned int)); //bigger than the length cmp result
	fgets(strtemp,4096,readFile);
	j=0;
	while( fgets(str[j++],4096,readFile)!=NULL ) //qulity from fastq format
	{
		readCount++;
		if(j%RECORDNUM==0)
		{
			cluRecord=Qcluster(str,RECORDNUM,outputbuf); //-1 is a useless length
			//-----test--------
			printf("word: %d---%d\n%s\n",cluRecord,strlen(outputbuf),outputbuf);
			//LinkSeq(str,cluRecord);
			fprintf(outFile,"%c",cluRecord);  //write flag--char-8bit or short int--16bit
			codeLen=CodeQuality(codebuf,outputbuf,3); //input the two char *
			fwrite(codebuf,sizeof(char),codeLen,outFile);
			codeCount+=codeLen; //size of the codefile
			j=0;
		}
		fgets(strtemp,4096,readFile);
		fgets(strtemp,4096,readFile);
		fgets(strtemp,4096,readFile);
	}
	if(j>0)  //last lines <RECORDNUM
	{
		fprintf(outFile,"%c",0); //write flag 8bit or===16bit
		p=0;
		while( (*(outputbuf+p)=*(str[0]+p))>32 )
				p++; //in the end, p is the length of Q seq
		// *(outputbuf+p)=0;
		for(i=1; i<j; i++) //the last score seq
		{
			strncpy(outputbuf+p*i,str[i],p);			
		}
		*(outputbuf+j*p)=0;
		codeLen=CodeQuality(codebuf,outputbuf,3);
		fwrite(codebuf,sizeof(char),codeLen,outFile); //write code
		codeCount+=codeLen;
	}

	free(outputbuf); free(strtemp); free(c); //free the memory
	for(i=0; i<RECORDNUM; i++)
		free(str[i]);

	//output information
	
	*/
	codebuf  =(char*)malloc(4096*sizeof(char)); //enough size
	codeCount=0;
	fgets(strtemp,4096,readFile);
	fgets(strtemp,4096,readFile);
	fgets(strtemp,4096,readFile);
	while( fgets(str,4096,readFile)!=NULL )
	{	
		readCount++;			
		codeLen=Code2Quality(codebuf, str, 3); 
		codeLen=fwrite(codebuf,1,codeLen,outFile);		
		codeCount+=codeLen;

		fgets(strtemp,4096,readFile);
		fgets(strtemp,4096,readFile);
		fgets(strtemp,4096,readFile);
	}
	fclose(readFile);  fclose(outFile);
	printf("code file size: %d\nreadCount: %d\n",codeCount,readCount);
	return 0;
}

#include <stdio.h>
#include <stdlib.h>

extern void complement(const char* str, const int len, char* rev_str);

int main(int argc, char** argv)
{
  	char* original;
  	char* complement_str;
  	FILE * pFile;
   	char* infile;
	if (infile == NULL)
	{
		fprintf(stderr,"no mem\n");
		return 1;
	}
   	if (argc > 1)
		infile = argv[1];
	else 
	{
		fprintf(stderr,"expected input file as argument 1\n");
	}

   	pFile = fopen (infile , "r");
	
   	if (pFile == NULL) perror ("Error opening file");
   	else 
	{
		// determine file size
   		fseek(pFile, 0L, SEEK_END);
   		int filesize = ftell(pFile);
	   	fseek(pFile, 0L, SEEK_SET);
	   	//fprintf(stdout,"filesize: %i\n", filesize); 

		original = (char*) malloc (filesize*sizeof(char));
		complement_str = (char*) malloc (filesize*sizeof(char));

     		//fgets (original , filesize , pFile);
     		fgets (original , filesize , pFile);
		complement(original, filesize,  complement_str);
     		//puts (original);
     		puts (complement_str);
     		fclose (pFile);

		free(original);
		free(complement_str);
   	}


   	return 0;
}




void complement(const char* str, const int len, char* rev_str) 
{
	int i;
	for (i=0; i<len-1; i++)
	{
		switch (str[i])
		{
			case 'A': rev_str[i] = 'T'; break;
			case 'a': rev_str[i] = 't'; break;
			case 'C': rev_str[i] = 'G'; break;
			case 'c': rev_str[i] = 'g'; break;
			case 'G': rev_str[i] = 'C'; break;
			case 'g': rev_str[i] = 'c'; break;
			case 'T': rev_str[i] = 'A'; break;
			case 't': rev_str[i] = 'a'; break;
			case 'n': rev_str[i] = 'n'; break;
			case 'N': rev_str[i] = 'N'; break;
			default : rev_str[i] = str[i]; break; 
		}
	
	}	
   	return;
}

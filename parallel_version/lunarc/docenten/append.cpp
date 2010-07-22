#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdarg.h>
#include <gutil.h>

void append(int nfile,char* filename[],int nprocess) {

	// Extracts and appends output from process run directories to main directory
	// This will only work for: text files with or without a single header line,
	// labels line, or leading blank line
	// Other files you will have to write your own software and extract from each
	// individual run directory
	
	xtring text,outfile,infile;
	int i,j,pos;
	FILE* in,*out;
	bool firstline;
	
	for (i=0;i<nfile;i++) {
		outfile=filename[i];
		out=fopen(outfile,"wt");
		if (out) {
			for (j=0;j<nprocess;j++) {
				infile.printf("run%d/%s",j+1,filename[i]);
				in=fopen(infile,"rt");
				if (!in) {
					printf("append: could not open %s for input, aborting for %s\n",(char*)infile,(char*)outfile);
					goto abort;
				}
				printf("Reading data from %s\n",(char*)infile);
				firstline=true;
				while (!feof(in)) {
					readfor(in,"A#",&text);
					if (firstline && text!="") {
						pos=text.findnotoneof(" .-\t");
						if (pos)
							if (text[pos]>='0' && text[pos]<='9') firstline=false;
					}
					if (firstline && j) readfor(in,"A#",&text);
					firstline=false;
					if (text!="")
						fprintf(out,"%s\n",(char*)text);
				}
				fclose(in);
			}
		}
		else {
			printf("append: warning - could not open %s for output\n",(char*)outfile);
		}
		
		fclose(out);
		printf("Written %s\n",(char*)outfile);
abort:
	}
}

int main(int argc,char* argv[]) {

	// Usage: append <node-count> <file1> <file2> ... <filen>

	xtring numstr;
	int nprocess;
	
	if (argc<2) {
		printf("Usage: %s <node-count> <file1> <file2> ... <filen>\n",argv[0]);
		exit(99);
	}
	
	numstr=argv[1];
	if (!numstr.isnum()) {
		printf("Usage: %s <node-count> <file1> <file2> ... <filen>\n",argv[0]);
		exit(99);
	}
	
	nprocess=numstr.num();
	
	append(argc-2,&argv[2],nprocess);
	
	return 0;
}

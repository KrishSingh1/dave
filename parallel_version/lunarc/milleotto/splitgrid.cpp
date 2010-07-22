#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdarg.h>
#include <gutil.h>

void splitgrid(xtring gridfile,int nprocess) {

	// Opens the specified gridlist file and splits it into <nprocess> roughly equal
	// segments, placing one in each run directory
	
	int gridno;
	int ngrid,p,fac,count;
	xtring outfile,line,dummy;
	FILE* out;
	
	FILE* in=fopen(gridfile,"rt");
	if (!in) {
		printf("splitgrid: could not open %s for input\n",(char*)gridfile);
		exit(99);
	}
		
	// Work out number of lines in gridfile (should be one for each lon/lat pair)
	
	ngrid=0;
	while (!feof(in)) {
		dummy="";
		readfor(in,"A#",&dummy);
		if (!feof(in) || dummy!="") ngrid++;
	}
	
	fac=ngrid/nprocess;
	if (ngrid%nprocess) fac++;
	
	rewind(in);
	
	// Get filename part of gridfile (which might be a path)
	
	p=gridfile.find('/');
	while (p!=-1) {
		gridfile=gridfile.mid(p+1);
		p=gridfile.find('/');
	}
	
	for (gridno=0;gridno<nprocess;gridno++) {

		outfile.printf("./run%d/%s",gridno+1,(char*)gridfile);
		out=fopen(outfile,"wt");
		if (!out) {
			printf("splitgrid: could not open %s for output\n",(char*)outfile);
			exit(99);
		}

		count=0;

		while (!feof(in) && ((count%fac) || count==0 ||  gridno==nprocess-1)) {
			line="";
			readfor(in,"A#",&line);
			if (!feof(in) || line!="") fprintf(out,"%s\n",(char*)line);
			count++;
		}

		fclose(out);
		
		printf("Created %s\n",(char*)outfile);
	}
	fclose(in);
}

int main(int argc,char* argv[]) {

	xtring gridfile,numstr;
	int nprocess;
	
	if (argc!=3) {
		printf("Usage: %s <gridlist-file> <node-count>\n",argv[0]);
		exit(99);
	}
	
	gridfile=argv[1];
	numstr=argv[2];
	if (!numstr.isnum()) {
		printf("Usage: %s <gridlist-file> <node-count>\n",argv[0]);
		exit(99);
	}
	
	nprocess=numstr.num();
	
	splitgrid(gridfile,nprocess);
	
	return 0;
}

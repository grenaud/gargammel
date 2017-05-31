/*
 * fasta2fastas
 * Date: Jul-07-2016 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>
#include <gzstream.h>

#include "utils.h"

using namespace std;

int main (int argc, char *argv[]) {


    const string usage=string("\n This program takes a fasta file a chromosome\n")+
        " and produces 2 files for each IUPAC code\n\n"+
        " "+string(argv[0])+" [in fasta] [out fasta prefix] "+"\n\n"+
	""; 

    if( (argc!= 3) ||
        (argc== 2 && string(argv[1]) == "-h") ||
        (argc== 2 && string(argv[1]) == "-help") ||
        (argc== 2 && string(argv[1]) == "--help") ){
        cout<<usage<<endl;
        return 1;
    }

   string line;
   igzstream filein;
   ogzstream fileout1;
   ogzstream fileout2;

   string filenamein  = string(argv[1]);
   string filenameout = string(argv[2]);
   string filenameout1 = filenameout+"1.fa.gz";
   string filenameout2 = filenameout+"2.fa.gz";

   fileout1.open(filenameout1.c_str(),ios::out);
   fileout2.open(filenameout2.c_str(),ios::out);

   filein.open(filenamein.c_str(), ios::in);
   
   if (filein.good()){
       while ( getline (filein,line)){
	   if(line[0] == '>'){
	       fileout1<<line<<endl;
	       fileout2<<line<<endl;	       
	   }else{
	       string l1="";
	       string l2="";

	       for(unsigned int i=0;i<line.size();i++){
		   char c=line[i];

		   if(isValidDNA(toupper(c))){
		       l1+=c;
		       l2+=c;
		   }else{
		       if(c ==    'R'){
			   char c1;
			   char c2;
			   
			   c1='A';
			   c2='G';

			   if(randomBool()){
			       l1+=c1;
			       l2+=c2;
			   }else{
			       l1+=c2;
			       l2+=c1;		      
			   }

		       }else{
		       if(c ==    'Y'){
			   char c1;
			   char c2;
			   
			   c1='C';
			   c2='T';

			   if(randomBool()){
			       l1+=c1;
			       l2+=c2;
			   }else{
			       l1+=c2;
			       l2+=c1;		      
			   }

		       }else{			  			       
		       if(c ==    'S'){

			   char c1;
			   char c2;
			   
			   c1='C';
			   c2='G';

			   if(randomBool()){
			       l1+=c1;
			       l2+=c2;
			   }else{
			       l1+=c2;
			       l2+=c1;		      
			   }
				   
		       }else{	       		   
		       if(c ==    'W'){

			   char c1;
			   char c2;
			   
			   c1='A';
			   c2='T';

			   if(randomBool()){
			       l1+=c1;
			       l2+=c2;
			   }else{
			       l1+=c2;
			       l2+=c1;		      
			   }
			   
		       }else{
		       if(c ==    'K'){
			   
			   char c1;
			   char c2;
			   
			   c1='G';
			   c2='T';

			   if(randomBool()){
			       l1+=c1;
			       l2+=c2;
			   }else{
			       l1+=c2;
			       l2+=c1;		      
			   }

		       }else{
		       if(c ==    'M'){
			   char c1;
			   char c2;
			   
			   c1='A';
			   c2='C';

			   if(randomBool()){
			       l1+=c1;
			       l2+=c2;
			   }else{
			       l1+=c2;
			       l2+=c1;		      
			   }			   
		       }else{

			   //bizarre corner case in vcf-consensus
			   if( i<=(line.size()-5)  && 
			       (line[i+0] == '<')  &&
			       (line[i+1] == 'D')  &&
			       (line[i+2] == 'E')  &&
			       (line[i+3] == 'L')  &&
			       (line[i+4] == '>')  ){
			       i = i+4;
			       l1+="NNNNN";
			       l2+="NNNNN";
			       continue;
			   }else{
			       cerr<<"Unknown base "<<c<<endl;
			       return 1;
			   }
		       }//M
		       }//K
		       }//S
		       }//W
		       }//Y
		       }//R
		   
		   }//iupac
	       }//each chr
	       fileout1<<l1<<endl;
	       fileout2<<l2<<endl;	       

	   }

       }
       filein.close();
   }else{
       cerr << "Unable to open fasta file "<<filenamein<<endl;
       return 1;
   }

   fileout1.close();
   fileout2.close();

    return 0;
}


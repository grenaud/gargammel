/*
 * fragSim
 * Date: Jun-06-2016 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>
#include <gzstream.h>
#include <random>


#include <api/SamHeader.h>
#include <api/BamMultiReader.h>
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <api/BamAux.h>
#include "PutProgramInHeader.h"

#include "FastQParser.h"
#include "libgab.h"



using namespace std;


typedef struct { 
    double s[12];
} substitutionRates;



typedef struct { 
    unsigned int s[4];
} baseCountInt;



int main (int argc, char *argv[]) {

    string file5p="/dev/stdout";
    string file3p="/dev/stdout";
    bool singleStr       =false;
    bool doubleStr       =false;
    bool singAnddoubleStr=false;


    vector<substitutionRates> sub5p;
    vector<substitutionRates> sub3p;

    string damagepatternFile5p;
    string damagepatternFile3p;
    
    bool hFormat=false;

    const string usage=string("\t"+string(argv[0])+
                              " [options]  [5p.dat] [3p.dat]"+"\n\n"+
			      " This program reads 2 substitution matrix files from damage-patterns (https://bitbucket.org/ustenzel/damage-patterns/) and produces .prof files\n\n"+
			      
			      "\t\t"+"-5p\t[output file]\tOutput profile for the 5' end (Default: "+stringify(file5p)+")\n"+
			      "\t\t"+"-3p\t[output file]\tOutput profile for the 3' end (Default: "+stringify(file3p)+")\n"+

			      "\n\n\tYou can specify either one of the two:\n"+
			      "\t\t"+"-single\t\t\tReport the deamination profile of a single strand library  (Default: "+booleanAsString( singleStr )+")\n"+
			      "\t\t"+"-double\t\t\tReport  the deamination profile of a double strand library  (Default: "+booleanAsString( doubleStr )+")\n"+
			      "\n\tor specify this option:\n"+
			      "\t\t"+"-both\t\t\tReport both C->T and G->A regardless of stand  (Default: "+booleanAsString( singAnddoubleStr )+")\n"+

			      "\t\t"+"-h\t\t\tMore human readible output (Default: "+booleanAsString(hFormat)+")\n"+


                              "");
                              
    if( (argc== 1) ||
        (argc== 2 && string(argv[1]) == "-h") ||
        (argc== 2 && string(argv[1]) == "-help") ||
        (argc== 2 && string(argv[1]) == "--help") ){
        cout<<"Usage:"<<endl;
        cout<<""<<endl;
        cout<<usage<<endl;
        return 1;
    }

    for(int i=1;i<(argc-1);i++){ //all but the last 3 args

        if(string(argv[i]) == "-h"  ){
            hFormat=true;
            continue;
        }


        if(string(argv[i]) == "-5p" ){
            file5p = string(argv[i+1]);
            i++;
            continue;
        }

        if(string(argv[i]) == "-3p" ){
            file3p = string(argv[i+1]);
            i++;
            continue;
        }

        if(string(argv[i]) == "-both" ){
            //allStr           = false;
            singleStr        = false;
            doubleStr        = false;
            singAnddoubleStr = true;
            continue;
        }


        if(string(argv[i]) == "-single" ){
            //allStr    = false;
            singleStr        = true;
            doubleStr        = false;
	    singAnddoubleStr = false;
            continue;
        }

        if(string(argv[i]) == "-double" ){
            //allStr    = false;
            singleStr        = false;
            doubleStr        = true;
	    singAnddoubleStr = false;
            continue;
        }
    }


    damagepatternFile5p          = string(argv[argc-2]);
    damagepatternFile3p          = string(argv[argc-1]);
    
    ofstream file5pFP;

    if(file5p == "/dev/stdout"){
        file5pFP.open(file5p.c_str(), ofstream::out | ofstream::app);
    }else{
        file5pFP.open(file5p.c_str());
    }

    
    if (!file5pFP.is_open()){
        cerr << "Unable to write to 5p file "<<file5p<<endl;
        exit(1);
    }

    
    ofstream file3pFP;
    if(file3p == "/dev/stdout"){
        file3pFP.open(file3p.c_str(), ofstream::out | ofstream::app);
    }else{
        file3pFP.open(file3p.c_str());
    }

    if (!file3pFP.is_open()){
        cerr << "Unable to write to 3p file "<<file3p<<endl;
        exit(1);
    }




    string line;
    igzstream misincorpFileSt5p;
    igzstream misincorpFileSt3p;

    /////////////////////
    //    5p end
    /////////////////////
    
    misincorpFileSt5p.open(damagepatternFile5p.c_str(), ios::in);


    if (misincorpFileSt5p.good()){
		
	//1 line for header	
	if ( !getline (misincorpFileSt5p,line)){
	    cerr << "Parsing error, cannot get header in substitution file "<<damagepatternFile5p<<endl;
	    return 1;
	}

	vector<string> fieldsHeader = allTokens(line,'\t');

	if(fieldsHeader.size() != 13     ||
	   fieldsHeader[1]     != "A>C"  ||
	   fieldsHeader[12]     != "T>G" ){		
	    cerr << "Parsing error, wrong header for substitution file "<<damagepatternFile5p<<" should be: A>C A>G A>T C>A C>G C>T G>A G>C G>T T>A T>C T>G"<<endl;
	    return 1;
	 }

	 while ( getline (misincorpFileSt5p,line)){
	     vector<string> fields = allTokens(line,'\t');

	     substitutionRates sub;

	     for(int i=0;i<12;i++){
		 string s = fields[1+i];
		 vector<string> fs = allTokens(s,' ');
		 sub.s[i] = destringify<double>( fs[0] );
		 //cerr<<i<<" fs="<<fs[0]<<"# s="<<sub.s[i]<<endl;
		 //cerr<<i<<" "<<fs[0]<<endl;
	     }
	     
	     sub5p.push_back( sub );
	 } //end for getline for each line
		
    }else{//if not good
	cerr << "Unable to open substitution file "<<damagepatternFile5p<<endl;
	return 1;
    }
    misincorpFileSt5p.close();




    /////////////////////
    //    3p end
    /////////////////////
    
    misincorpFileSt3p.open(damagepatternFile3p.c_str(), ios::in);


    if (misincorpFileSt3p.good()){
		
	//1 line for header	
	if ( !getline (misincorpFileSt3p,line)){
	    cerr << "Parsing error, cannot get header in substitution file "<<damagepatternFile3p<<endl;
	    return 1;
	}

	vector<string> fieldsHeader = allTokens(line,'\t');

	if(fieldsHeader.size() != 13     ||
	   fieldsHeader[1]     != "A>C"  ||
	   fieldsHeader[12]     != "T>G" ){		
	    cerr << "Parsing error, wrong header for substitution file "<<damagepatternFile3p<<" should be: A>C A>G A>T C>A C>G C>T G>A G>C G>T T>A T>C T>G"<<endl;
	    return 1;
	 }

	 while ( getline (misincorpFileSt3p,line)){
	     vector<string> fields = allTokens(line,'\t');

	     substitutionRates sub;

	     for(int i=0;i<12;i++){
		 string s = fields[1+i];
		 vector<string> fs = allTokens(s,' ');
		 sub.s[i] = destringify<double>( fs[0] );		 
	     }
	     
	     sub3p.push_back( sub );
	 } //end for getline for each line
		
    }else{//if not good
	cerr << "Unable to open substitution file "<<damagepatternFile3p<<endl;
	return 1;
    }
    misincorpFileSt3p.close();

    //flip because the last line is the one next to the 3' end
    reverse(sub3p.begin(),sub3p.end());



    if(sub5p.size() != sub3p.size()){
	cerr << "Error in substitution file "<<damagepatternFile5p<<" the number of defined bases for the 5' and 3' ends differ"<<endl;
	return 1;
    }
    if(hFormat)
        file5pFP<<"pos\t";

    file5pFP<<"A>C\tA>G\tA>T\tC>A\tC>G\tC>T\tG>A\tG>C\tG>T\tT>A\tT>C\tT>G"<<endl;

    //compute actual substitutions 
    for(unsigned int i=0;i<sub5p.size();i++){
        if(hFormat)
            file5pFP<<printIntAsWhitePaddedString(i,int(log10(sub5p.size()))+1)<<"\t";
        
	vector<double> toprint;
	for(unsigned int k=0;k<12;k++){
	    toprint.push_back(sub5p[i].s[k]);
	    //cerr<<i<<" "<<k<<" "<<sub5p[i].s[k]<<endl;
	}
	
	file5pFP<<vectorToString(toprint,"\t")<<endl;

    }


    if(hFormat)
        file3pFP<<"pos\t";

    file3pFP<<"A>C\tA>G\tA>T\tC>A\tC>G\tC>T\tG>A\tG>C\tG>T\tT>A\tT>C\tT>G"<<endl;
    for(unsigned int i=0;i<sub3p.size();i++){
        if(hFormat){	    
	    file3pFP<<""<<printIntAsWhitePaddedString(i,int(log10(sub3p.size()))+1)<<"\t";      
        }

    
	vector<double> toprint;
	for(unsigned int k=0;k<12;k++){
	    toprint.push_back(sub3p[i].s[k]);
	    //cerr<<i<<" "<<k<<" "<<sub3p[i].s[k]<<endl;
	}
	
	file3pFP<<vectorToString(toprint,"\t")<<endl;
    }

    file5pFP.close();
    file3pFP.close();

    return 0;
}


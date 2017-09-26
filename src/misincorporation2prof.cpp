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
#include "utils.h"



using namespace std;


typedef struct { 
    double s[12];
} substitutionRates;

typedef struct { 
    unsigned int  s[12];
} substitutionRatesInt;


typedef struct { 
    unsigned int s[4];
} baseCountInt;



int main (int argc, char *argv[]) {



    vector<substitutionRates> sub5p;
    vector<substitutionRates> sub3p;

    string mapdamageFile;
    
    const string usage=string("\t"+string(argv[0])+
                              " [options]  [mis.txt]"+"\n\n"+
			      " This program reads a misincorporation file and produces .prof files"+

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


    mapdamageFile          = string(argv[argc-1]);


    // if(matrixFileSpecified && matrixSpecified){
    // 	cerr << "Do not specify both -matrix and -matrixfile "<<endl;
    // 	return 1;
    // }


    vector<substitutionRatesInt> subCount5p;
    vector<substitutionRatesInt> subCount3p;
    
    vector<baseCountInt> baseCount5p;
    vector<baseCountInt> baseCount3p;
    bool first5p=true;
    bool first3p=true;
    string line;
    igzstream misincorpFileSt;
    
    misincorpFileSt.open(mapdamageFile.c_str(), ios::in);


    if (misincorpFileSt.good()){
		
	//3 lines for header
	for(int i=0;i<3;i++){
	    if ( !getline (misincorpFileSt,line)){
		cerr << "Parsing error for misincorporation file "<<mapdamageFile<<endl;
		return 1;
	    }
	}
		
	if ( !getline (misincorpFileSt,line)){
	    cerr << "Parsing error for misincorporation file "<<mapdamageFile<<endl;
	    return 1;
	}

	int indexHeaderToTrueIndex [12];

	vector<string> fieldsHeader = allTokens(line,'\t');
	for(int i=0;i<12;i++){		    
		    
	    if( (fieldsHeader[9+i].size() != 3) || (fieldsHeader[9+i][1] != '>') ){
		cerr << "Parsing error for misincorporation file "<<mapdamageFile<<" wrong header field "<<fieldsHeader[9+i]<<endl;
		return 1;
	    }

	    int trueIndex=-1;
	    if(fieldsHeader[9+i][0]     == 'A'){
		if(fieldsHeader[9+i][2] == 'C')
		    trueIndex=0;
		if(fieldsHeader[9+i][2] == 'G')
		    trueIndex=1; 
		if(fieldsHeader[9+i][2] == 'T')
		    trueIndex=2; 
	    }else{
		if(fieldsHeader[9+i][0]     == 'C'){
		    if(fieldsHeader[9+i][2] == 'A')
			trueIndex=3; 
		    if(fieldsHeader[9+i][2] == 'G')
			trueIndex=4;
		    if(fieldsHeader[9+i][2] == 'T')
			trueIndex=5;
		}else{
		    if(fieldsHeader[9+i][0]     == 'G'){
			if(fieldsHeader[9+i][2] == 'A')
			    trueIndex=6;
			if(fieldsHeader[9+i][2] == 'C')
			    trueIndex=7; 
			if(fieldsHeader[9+i][2] == 'T')
			    trueIndex=8;
		    }else{
			if(fieldsHeader[9+i][0]     == 'T'){
			    if(fieldsHeader[9+i][2] == 'A')
				trueIndex=9; 
			    if(fieldsHeader[9+i][2] == 'C')
				trueIndex=10; 
			    if(fieldsHeader[9+i][2] == 'G')
				trueIndex=11; 
			}else{
			    cerr << "Parsing error for misincorporation file "<<mapdamageFile<<endl; return 1;
			}
		    }
		}
	    }
	    if(trueIndex==-1){
		cerr << "Parsing error for misincorporation file "<<mapdamageFile<<endl; return 1;
	    }
	    indexHeaderToTrueIndex[i] = trueIndex;
	    //cout<<i<<"\t"<<trueIndex<<"\t"<<fieldsHeader[9+i]<<endl;
	}

		
	while ( getline (misincorpFileSt,line)){
	    vector<string> fields = allTokens(line,'\t');

	    //3prime end
	    if(fields[1] == "3p"){
		if( first3p && (fields[2] == "-") ){
		    first3p=false;
		}
			
		if(first3p){//first iteration
		    substitutionRatesInt sub;
		    baseCountInt         bsc;
		    for(int i=0;i<4;i++)
			bsc.s[i] = destringify<unsigned int>(fields[4+i]);

		    for(int i=0;i<12;i++){
			// cout<<i<<" "<<indexHeaderToTrueIndex[i]<<endl;
			sub.s[ indexHeaderToTrueIndex[i] ] = destringify<unsigned int>(fields[9+i]);
		    }
		    subCount3p.push_back( sub);
		    baseCount3p.push_back(bsc);
		}else{
		    int indexToUse = destringify<int>(fields[3])-1;
			    
		    for(int i=0;i<4;i++)
			baseCount3p[indexToUse].s[i] += destringify<unsigned int>(fields[4+i]);
			    
		    for(int i=0;i<12;i++)
			subCount3p[ indexToUse].s[ indexHeaderToTrueIndex[i] ] += destringify<unsigned int>(fields[9+i]);
		}
		continue;
	    }



	    //5prime end
	    if(fields[1] == "5p"){
		if( first5p && (fields[2] == "-") ){
		    first5p=false;
		}
			
		if(first5p){//first iteration
		    substitutionRatesInt sub;
		    baseCountInt         bsc;
		    for(int i=0;i<4;i++)
			bsc.s[i] = destringify<unsigned int>(fields[4+i]);

		    for(int i=0;i<12;i++)
			sub.s[ indexHeaderToTrueIndex[i] ] = destringify<unsigned int>(fields[9+i]);
			
		    subCount5p.push_back( sub);
		    baseCount5p.push_back(bsc);
		}else{
		    int indexToUse = destringify<int>(fields[3])-1;
			    
		    for(int i=0;i<4;i++)
			baseCount5p[ indexToUse].s[i] += destringify<unsigned int>(fields[4+i]);

		    for(int i=0;i<12;i++)
			subCount5p[indexToUse].s[ indexHeaderToTrueIndex[i] ]   += destringify<unsigned int>(fields[9+i]);
		}
		continue;
	    }

	    //error if not 5p or 3p
	    cerr << "line from misincorporation.txt does not have a 5p or 3p in the 2nd column: "<<line<<endl;
	    return 1;

		    
		    
	}//end for getline for each line
		
    }else{//if not good
	cerr << "Unable to open misincorporation file "<<mapdamageFile<<endl;
	return 1;
    }

    if(subCount5p.size() != subCount3p.size()){
	cerr << "Error in misincorporation file "<<mapdamageFile<<" the number of defined bases for the 5' and 3' ends differ"<<endl;
	return 1;
    }

    cout<<"pos\tA>C\tA>G\tA>T\tC>A\tC>G\tC>T\tG>A\tG>C\tG>T\tT>A\tT>C\tT>G"<<endl;

    //compute actual substitutions 
    for(unsigned int i=0;i<subCount5p.size();i++){
	substitutionRates toadd;
	for(int b=0;b<4;b++){				    
	    for(int a=0;a<3;a++){		
		int idx = b*3+a;
		toadd.s[idx] = double(subCount5p[i].s[idx]) / double(baseCount5p[i].s[b]);
		// cout<<i<<" "<<idx<<" "<<subCount5p[i].s[idx]<<" "<<baseCount5p[i].s[b]<<endl;
	    }		    
	}

	sub5p.push_back(toadd);

	cout<<i;
	for(unsigned int k=0;k<12;k++)
	    cout<<"\t"<<sub5p[i].s[k];
	cout<<endl;

    }

    cout<<"pos\tA>C\tA>G\tA>T\tC>A\tC>G\tC>T\tG>A\tG>C\tG>T\tT>A\tT>C\tT>G"<<endl;
    for(unsigned int i=0;i<subCount3p.size();i++){
	substitutionRates toadd;
	for(int b=0;b<4;b++){				    
	    for(int a=0;a<3;a++){		
		int idx = b*3+a;
		toadd.s[idx] = double(subCount3p[i].s[idx]) / double(baseCount3p[i].s[b]);
	    }		    
	}
	sub3p.push_back(toadd);


	cout<<i;
	for(unsigned int k=0;k<12;k++)
	    cout<<"\t"<<sub3p[i].s[k];
	cout<<endl;

    }

    return 0;
}


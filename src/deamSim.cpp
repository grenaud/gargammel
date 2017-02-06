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

//#define DEBUG

using namespace std;

const uint32_t flagSingleReads =  4; // 00000100
const uint32_t flagFirstPair   = 77; // 01001101
const uint32_t flagSecondPair  =141; // 10001101


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
    string outFastagz                ;
    bool   outFastagzb         =false;

    string outBAM                    ;
    bool   outBAMb             =false;
    bool produceUnCompressedBAM=false;
    bool putDeamInName         =false;

    bool verbose               =false;
    //Deamination parameters
    bool useBriggs             =false;
    double vBrgs = 0.0;
    double lBrgs = 0.0;
    double dBrgs = 0.0;
    double sBrgs = 0.0;


    bool useBriggsss             =false;
    double sBrgs5p = 0.0;
    double sBrgs3p = 0.0;
    double lBrgs5p = 0.0;
    double lBrgs3p = 0.0;




    default_random_engine generator;
    geometric_distribution<int> overhang;

    geometric_distribution<int> overhanggeo5p;
    geometric_distribution<int> overhanggeo3p;

    vector<substitutionRates> sub5p;
    vector<substitutionRates> sub3p;

    string matrixFile;
    bool   matrixFileSpecified=false;

    string mapdamageFile;
    bool   mapdamageFileSpecified=false;
    string mapdamageProtocol;

    string matrix;
    bool   matrixSpecified    = false;
    bool   singleDeam         = false;
    bool   doubleDeam         = false;
    
    const string usage=string("\t"+string(argv[0])+
                              " [options]  [fasta or BAM file]"+"\n\n"+
			      " This program reads a fasta (default) or BAM file containing aDNA fragments and\n"+
			      " adds deamination according to a certain model file\n"+
			      " some model files are found the models/ directory\n"+
			      " if the input is fasta, the output will be fasta as well\n"+
			      "\n"+

			      " I/O Options:\n"+
			      "\t\t"+"-b\t"+"[BAM out]"   +"\t\t\t"+"Read BAM and write output as a BAM (default: fasta)"+"\n"+
			      "\t\t"+"-u\t"               +"\t\t\t\t"+"Produce uncompressed BAM (good for unix pipe)"+"\n"+

			      "\t\t"+"-o\t"+"[fasta out]" +"\t\t\t"+"Write fasta output as a zipped fasta"+"\n"+
			      "\t\t"+"-name"              +"\t\t\t\t\t"+"Put a tag in the read name with deam bases (Default "+booleanAsString(putDeamInName)+")"+"\n"+
			      "\t\t"+"-v\t"+""            +"\t\t\t\t"+"verbose mode"+"\n"+

			      "\n"
			      +" Mandatory deamination options:\n"+ 
			      "\tSpecify either:\n"+
			      
                              "\t\t"+"-mapdamage  [mis.txt] [protocol]" +"\t"+"Read the miscorporation file [mis.txt]"+"\n"+
                              "\t\t"+"                               " +"\t\t"+"produced by mapDamage"+"\n"+
                              "\t\t"+"                               " +"\t\t"+"[protocol] can be either \"single\" or \"double\" (without quotes)\n"+
                              "\t\t"+"                               " +"\t\t"+"Single strand will have C->T damage on both ends"+"\n"+
                              "\t\t"+"                               " +"\t\t"+"Double strand will have and C->T at the 5' end and G->A damage at the 3' end"+"\n"+
			      "\n"+

                              "\t\t"+"-matfile  [matrix file prefix]" +"\t\t"+"Read the matrix file of substitutions instead of the default "+"\n"+
                              "\t\t"+"                               " +"\t\t"+"Provide the prefix only, both files must end with"+"\n"+
                              "\t\t"+"                               " +"\t\t"+"5.dat and 3.dat"+"\n"+
			      "\n"+

                              "\t\t"+"-mat      [ancient DNA matrix]" +"\t\t"+"For default matrices, use either \"single\" or \"double\" (without quotes)"+"\n"+
                              "\t\t"+"                               " +"\t\t"+"Single strand will have C->T damage on both ends"+"\n"+
                              "\t\t"+"                               " +"\t\t"+"Double strand will have and C->T at the 5p end and G->A damage at the 3p end"+"\n"+

			      "\t\t"+"-damage     [v,l,d,s]" +"\t\t\t"+"For the Briggs et al. 2007 model"+"\n"+
                              "\t\t"+"                               " +"\t\t"+"The parameters must be comma-separated e.g.: -damage 0.03,0.4,0.01,0.7"+"\n"+
                              "\t\t"+"                               " +"\t\t"+"\tv: nick frequency"+"\n"+
                              "\t\t"+"                               " +"\t\t"+"\tl: length of overhanging ends (geometric parameter)"+"\n"+
                              "\t\t"+"                               " +"\t\t"+"\td: prob. of deamination of Cs in double-stranded parts"+"\n"+
                              "\t\t"+"                               " +"\t\t"+"\ts: prob. of deamination of Cs in single-stranded parts"+"\n"+

			      // "\t\t"+"-damagess     [d,s5,s3,l5,l3]" +"\t\t\t"+"Single-strand model"+"\n"+
                              // "\t\t"+"                               " +"\t\t"+"The parameters must be comma-separated e.g.: -damagess 0.01,0.5,0.7,0.6,0.7"+"\n"+
                              // "\t\t"+"                               " +"\t\t"+"\td: prob. of deamination of Cs in double-stranded parts"+"\n"+
                              // "\t\t"+"                               " +"\t\t"+"\ts5: prob. of deamination of Cs in single-stranded parts at the 5p end"+"\n"+
                              // "\t\t"+"                               " +"\t\t"+"\ts3: prob. of deamination of Cs in single-stranded parts at the 3p end"+"\n"+
                              // "\t\t"+"                               " +"\t\t"+"\tl5: length of overhanging ends at the 5p end (geometric parameter)"+"\n"+
                              // "\t\t"+"                               " +"\t\t"+"\tl3: length of overhanging ends at the 3p end (geometric parameter)"+"\n"+

                              // "\t\t"+"                               " +"\t\t"+"\tv: nick frequency"+"\n"+

                              // "\t\t"+"                               " +"\t\t"+"\td: prob. of deamination of Cs in double-stranded parts"+"\n"+

			      // dBrgs   = destringify<double>(temps[0]);
			      // sBrgs5p = destringify<double>(temps[1]);
			      // sBrgs3p = destringify<double>(temps[2]);
			      // lBrgs5p = destringify<double>(temps[3]);
			      // lBrgs3p = destringify<double>(temps[4]);


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

#ifdef DEBUG
    string  substiString[]={"A>C","A>G","A>T","C>A","C>G","C>T","G>A","G>C","G>T","T>A","T>C","T>G"};
#endif

    for(int i=1;i<(argc-1);i++){ //all but the last 3 args
	

	if(string(argv[i]) == "-mat" ){
	    matrix          = string(argv[i+1]);

	    if(matrix == "single"){
		singleDeam=true;
	    }else{
		if(matrix == "double" ){
		    doubleDeam=true;
		}else{
		    cerr << "Specify either \"single\" or \"double\" (without quotes), you entered "<<matrix<<endl;
		    return 1;
		}
	    }

	    matrixSpecified = true;
	    i++;
	    continue;
	}

	if(string(argv[i]) == "-matfile" ){
	    matrixFile          = string(argv[i+1]);
	    matrixFileSpecified = true;
	    i++;
	    continue;
	}

	if(string(argv[i]) == "-mapdamage" ){
	    mapdamageFile          = string(argv[i+1]);
	    mapdamageProtocol      = string(argv[i+2]);
	    mapdamageFileSpecified = true;

	    if(mapdamageProtocol != "single" &&
	       mapdamageProtocol != "double" ){
		cerr << "The protocol specified should be 'single' or 'double' without the quotes, got: "<<mapdamageProtocol<<endl;
		return 1;		
	    }
	    i++;
	    i++;
	    continue;
	}

	if(string(argv[i]) == "-damage" ){
	    string parametersB          = string(argv[i+1]);
	    useBriggs                   = true;
	    vector<string> temps=allTokens(parametersB,',');
	    if(temps.size()!=4){
		cerr << "Specify 4 comma-separated values for the Briggs model, you entered:"<<parametersB<<endl;
		return 1;
	    }
	    vBrgs=destringify<double>(temps[0]);
	    lBrgs=destringify<double>(temps[1]);
	    dBrgs=destringify<double>(temps[2]);
	    sBrgs=destringify<double>(temps[3]);	    
	    overhang = geometric_distribution<int>(lBrgs);
	    if( (vBrgs<0) || (vBrgs>1) ){ cerr << "Nick frequency must be between 0 and 1, entered:"<<vBrgs<<endl; return 1; }
	    if( (lBrgs<0) || (lBrgs>1) ){ cerr << "Geometric parameter for overhang must be between 0 and 1, entered:"<<lBrgs<<endl; return 1; }
	    if( (dBrgs<0) || (dBrgs>1) ){ cerr << "double-strand deamination prob. must be between 0 and 1, entered:"<<dBrgs<<endl; return 1; }
	    if( (sBrgs<0) || (sBrgs>1) ){ cerr << "single-strand deamination prob. must be between 0 and 1, entered:"<<sBrgs<<endl; return 1; }

	    i++;

	    continue;
	}



	if(string(argv[i]) == "-damagess" ){
	    string parametersB          = string(argv[i+1]);
	    useBriggsss                   = true;
	    vector<string> temps=allTokens(parametersB,',');
	    if(temps.size()!=5){
		cerr << "Specify 4 comma-separated values for the Briggs model, you entered:"<<parametersB<<endl;
		return 1;
	    }
	    //vBrgs=destringify<double>(temps[0]);

	    dBrgs   = destringify<double>(temps[0]);
	    sBrgs5p = destringify<double>(temps[1]);
	    sBrgs3p = destringify<double>(temps[2]);
	    lBrgs5p = destringify<double>(temps[3]);
	    lBrgs3p = destringify<double>(temps[4]);
	    
	    overhanggeo5p = geometric_distribution<int>(lBrgs5p);
	    overhanggeo3p = geometric_distribution<int>(lBrgs3p);

	    //if( (vBrgs<0) || (vBrgs>1) ){ cerr << "Nick frequency must be between 0 and 1, entered:"<<vBrgs<<endl; return 1; }
	    if( (lBrgs5p<0) || (lBrgs5p>1) ){ cerr << "Geometric parameter for overhang must be between 0 and 1, entered:"<<lBrgs5p<<endl; return 1; }
	    if( (lBrgs3p<0) || (lBrgs3p>1) ){ cerr << "Geometric parameter for overhang must be between 0 and 1, entered:"<<lBrgs3p<<endl; return 1; }

	    if( (dBrgs<0) || (dBrgs>1) ){ cerr << "double-strand deamination prob. must be between 0 and 1, entered:"<<dBrgs<<endl; return 1; }
	    if( (sBrgs5p<0) || (sBrgs5p>1) ){ cerr << "single-strand deamination prob. must be between 0 and 1, entered:"<<sBrgs5p<<endl; return 1; }
	    if( (sBrgs3p<0) || (sBrgs3p>1) ){ cerr << "single-strand deamination prob. must be between 0 and 1, entered:"<<sBrgs3p<<endl; return 1; }

	    i++;

	    continue;
	}

	if(string(argv[i]) == "-name" ){
	    putDeamInName = true;
	    continue;
	}

        if(string(argv[i]) == "-b" ){
            outBAM  = string(argv[i+1]);
            i++;
            outBAMb = true;
            continue;
        }

        if(string(argv[i]) == "-u"  ){
            produceUnCompressedBAM=true;
            continue;
        }

        if(string(argv[i]) == "-v"  ){
            verbose=true;
            continue;
        }

        if(string(argv[i]) == "-o" ){
            outFastagz  = string(argv[i+1]);
            i++;
            outFastagzb = true;
            continue;
        }

        cerr<<"Error: unknown option "<<string(argv[i])<<endl;
        return 1;
    }

    if(useBriggs   && (matrixSpecified || mapdamageFileSpecified)){
        cerr<<"Error: cannot specify both a file and Briggs parameter model"<<endl;
        return 1;
    }

    if(useBriggsss && (matrixSpecified || mapdamageFileSpecified)){
        cerr<<"Error: cannot specify both a file and Briggs parameter model"<<endl;
        return 1;
    }

    if(useBriggs && useBriggsss){
        cerr<<"Error: cannot specify both a file and Briggs parameter model"<<endl;
        return 1;
    }

    if(!useBriggs           && 
       !useBriggsss         && 
       !matrixSpecified     && 
       !matrixFileSpecified &&
       !mapdamageFileSpecified){
    	cerr << "Please specify the matrix to use or use the Briggs parameter model"<<endl;
    	return 1;
    }


    if(outFastagzb && outBAMb){ 
        cerr<<"Error: cannot specify both -o and -b"<<endl;
        return 1;
    }
    string inputFile = string(argv[argc-1]);


    // if(matrixFileSpecified && matrixSpecified){
    // 	cerr << "Do not specify both -matrix and -matrixfile "<<endl;
    // 	return 1;
    // }

    if(!useBriggs && !useBriggsss){    
	if(mapdamageFileSpecified){//we specified a mapDamage file

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
#ifdef DEBUG
		cout<<i<<"\t";
		for(unsigned int k=0;k<12;k++)
		    cout<<sub5p[i].s[k]<<" ";
		cout<<endl;
#endif
	    }

	    for(unsigned int i=0;i<subCount3p.size();i++){
		substitutionRates toadd;
		for(int b=0;b<4;b++){				    
		    for(int a=0;a<3;a++){		
			int idx = b*3+a;
			toadd.s[idx] = double(subCount3p[i].s[idx]) / double(baseCount3p[i].s[b]);
		    }		    
		}
		sub3p.push_back(toadd);

#ifdef DEBUG
		cout<<i<<"\t";
		for(unsigned int k=0;k<12;k++)
		    cout<<sub3p[i].s[k]<<" ";
		cout<<endl;
#endif

	    }
	    // return 1;
	    matrixSpecified=true; //fake a .dat file
	    if(mapdamageProtocol == "single" ){
		singleDeam=true;
	    }
	    
	    if(mapdamageProtocol == "double" ){
		doubleDeam=true;
	    }
	    
	}else{//else, we specify a .dat matrix
	    string deam5File;
	    string deam3File;
	
	    if(matrixFileSpecified){
		deam5File = matrixFile+"5.dat";
		deam3File = matrixFile+"3.dat";	
	    }else{   
		deam5File = getCWD(argv[0])+ "matrices/"+matrix+"-5.dat";
		deam3File = getCWD(argv[0])+ "matrices/"+matrix+"-3.dat";	
	    }
	
	    string line;

	    ///////////////////////////////////////
	    //READ file with 5p deamination rates /
	    ///////////////////////////////////////
	    igzstream deam5pFileSt;

	    deam5pFileSt.open(deam5File.c_str(), ios::in);


	    if (deam5pFileSt.good()){
		//header
		if ( !getline (deam5pFileSt,line)){
		    cerr << "Parsing error for 5p deamination file "<<deam5File<<endl;
		    return 1;
		}

		while ( getline (deam5pFileSt,line)){
		    vector<string> fields = allTokens(line,'\t');
		    vector<string> firstFields;
		    for(unsigned int i=1;i<fields.size();i++){
			vector<string> fields2 = allTokens(fields[i],' ');
			firstFields.push_back(fields2[0]);
		    }
	    
		    if(firstFields.size()!=12){
			cerr << "line from deamination does not have 12 fields "<<line<<" "<<firstFields.size()<<endl;
			return 1;
		    }
	    
		    substitutionRates toadd;
		    for(unsigned int i=0;i<firstFields.size();i++){		
			toadd.s[i]=destringify<double>(firstFields[i]);
#ifdef DEBUG
			cout<<sub5p.size()<<"\t"<<substiString[i]<<"\t"<<toadd.s[i]<<"\t"<<firstFields[i]<<endl;
#endif
		    }	


		    for(int b=0;b<4;b++){		
			double sum =0.0;
			for(int a=0;a<3;a++){		
			    int idx = b*3+a;
			    sum+=toadd.s[idx];
			}
			if(sum>1.0){
			    cerr<<"Problem with line "<<line<<" in file "<<deam5File<<" the sum of the substitution probabilities exceeds 1"<<endl;
			    return 1;
			}
		    }

		    sub5p.push_back(toadd);
		    //return 1;
		}
	
		deam5pFileSt.close();
	    }else{
		cerr << "Unable to open 5p deamination file "<<deam5File<<endl;
		return 1;
	    }


	    ///////////////////////////////////////
	    //READ file with 3p deamination rates /
	    ///////////////////////////////////////
	    igzstream deam3pFileSt;

	    deam3pFileSt.open(deam3File.c_str(), ios::in);


	    if (deam3pFileSt.good()){
		//header
		if ( !getline (deam3pFileSt,line)){
		    cerr << "Parsing error for 3p deamination file "<<deam3File<<endl;
		    return 1;
		}

		while ( getline (deam3pFileSt,line)){
		    vector<string> fields = allTokens(line,'\t');
		    vector<string> firstFields;
		    for(unsigned int i=1;i<fields.size();i++){
			vector<string> fields2 = allTokens(fields[i],' ');
			firstFields.push_back(fields2[0]);
		    }
	    
		    if(firstFields.size()!=12){
			cerr << "line from deamination does not have 12 fields "<<line<<" "<<firstFields.size()<<endl;
			return 1;
		    }
	    
		    substitutionRates toadd;
		    for(unsigned int i=0;i<firstFields.size();i++){		
			toadd.s[i]=destringify<double>(firstFields[i]);

#ifdef DEBUG
			cout<<sub3p.size()<<"\t"<<substiString[i]<<"\t"<<toadd.s[i]<<"\t"<<firstFields[i]<<endl;
#endif

			//cout<<toadd.s[i]<<"\t"<<firstFields[i]<<endl;
		    }	

		    for(int b=0;b<4;b++){		
			double sum =0.0;
			for(int a=0;a<3;a++){		
			    int idx = b*3+a;
			    sum+=toadd.s[idx];
			}
			if(sum>1.0){
			    cerr<<"Problem with line "<<line<<" in file "<<deam3File<<" the sum of the substitution probabilities exceeds 1"<<endl;
			    return 1;
			}
		    }
		

		    sub3p.push_back(toadd);
		    //return 1;
		}
	
		deam3pFileSt.close();
	    }else{
		cerr << "Unable to open 3p deamination file "<<deam3File<<endl;
		return 1;
	    }

	    // for(unsigned int i=0;i<sub3p.size();i++){		
	    // 	cout<<i<<"\t"<<arrayToString(sub3p[i].s,12,"\t")<<endl;
	    // }
	    reverse(sub3p.begin(),sub3p.end());
	
	    // for(unsigned int i=0;i<sub3p.size();i++){		
	    // 	cout<<i<<"\t"<<arrayToString(sub3p[i].s,12,"\t")<<endl;
	    // }
	}
    }else{
	//no file to read
    }

    //Fastq
    FastQParser * fp=NULL;
    FastQObj * fo;
    ogzstream outFastagzfp;

    //BAM
    BamReader reader;
    BamWriter writer;

    
    if(outBAMb){


	if ( reader.Open(inputFile) ) {
	    SamHeader  myHeader=reader.GetHeader();
	    SamProgram sp;
   
	    string pID          = "deamSim";   
	    string pName        = "deamSim";   
	    string pCommandLine = "";
	    for(int i=0;i<(argc);i++){
		pCommandLine += (string(argv[i])+" ");
	    }

	    putProgramInHeader(&myHeader,pID,pName,pCommandLine,returnGitHubVersion(string(argv[0]),".."));
	    if(produceUnCompressedBAM)  
		writer.SetCompressionMode(BamWriter::Uncompressed); 
	    
	    if( !writer.Open(outBAM,myHeader,reader.GetReferenceData() ) ) {
		cerr << "Could not open output BAM file  "<<outBAM << endl;
		return 1;   
	    }

	}else{
	    cerr << "Could not open input BAM file  "<<inputFile<< endl;
	    return 1;   
	}
    }else{
	fp = new FastQParser(inputFile,true);
	if(outFastagzb){

	    outFastagzfp.open(outFastagz.c_str(), ios::out);
	    if(!outFastagzfp.good()){       cerr<<"Cannot write to file "<<outFastagz<<endl; return 1; }
	    	   
	}
    }

    unsigned int   f       = 0;

    
    while(true){
	string def;
	string seq;

	BamAlignment al;

	if(outBAMb){
	    if(!reader.GetNextAlignment(al))
		break;
	    def = al.Name;
	    seq = al.QueryBases;
	}else{
	    if(!fp->hasData())
		break;
	    fo  = fp->getData();
	    def = *(fo->getID());
	    seq = *(fo->getSeq());
	}
	transform(seq.begin(), seq.end(),seq.begin(), ::toupper);

	bool deaminated=false;
	vector<int> deamPos;
	vector<int> deamPos3p;

	//cout<<dBrgs<<endl;
	if(f%100000 == 0 && f!=0){
	    cerr<<"Produced "<<thousandSeparator(f)<<endl;
	}
	f++;

	/////////////////////////
	//  ADDING DEAMINATION //
	/////////////////////////

	if(useBriggsss){
	    
	    int overhang5p=overhanggeo5p(generator);
	    int overhang3p=overhanggeo3p(generator);

	    if( (overhang5p+overhang3p)>=int(seq.size())){//all single strand

		//just apply sBrgs over the entire fragment
		for(int i=0;i<int(seq.size());i++){

		    //both in 5p and 3p overhang, pick closest
		    if( ((i+1)<=overhang5p) && ((int(seq.size())-i)<=overhang3p) ){
			if(i<int(double(seq.size())*0.5)){
			    if(seq[i] == 'C'){
				if( randomProb() < sBrgs5p ){
				    deamPos.push_back(i+1);
				    seq[i] = 'T';
				    deaminated=true;
				}
			    }
			}else{
			    if(seq[i] == 'C'){
				if( randomProb() < sBrgs3p ){
				    deamPos.push_back(i+1);
				    seq[i] = 'T';
				    deaminated=true;
				}
			    }
			}
			
		    }else{//just in one of them
			//just in 3p overhand
			if( !((i+1)<=overhang5p) &&  ((int(seq.size())-i)<=overhang3p) ){
			    if(seq[i] == 'C'){
				if( randomProb() < sBrgs3p ){
				    deamPos.push_back(i+1);
				    seq[i] = 'T';
				    deaminated=true;
				}
			    }

			}else{

			    //just 5' overhang
			    if(  ((i+1)<=overhang5p) && !((int(seq.size())-i)<=overhang3p) ){

				if(seq[i] == 'C'){
				    if( randomProb() < sBrgs5p ){
					deamPos.push_back(i+1);
					seq[i] = 'T';
					deaminated=true;
				    }
				}

			    }else{
				
				cerr<<"Internal error in full overhang, contact developers "<<al.Name<<endl; 
				return 1; 

			    }
			}		       
		    }//closes just in one of them
		}//end for each base
		

	    }else{
	
		for(int i=0;i<int(seq.size());i++){



			

		    //in 5p overhang
		    if((i+1)<=overhang5p){
#ifdef DEBUG
			cout<<"5p"<<endl;
#endif
			if(seq[i] == 'C'){
			    if( randomProb() < sBrgs5p ){
				deamPos.push_back(i+1);
				seq[i] = 'T';
				deaminated=true;
			    }
			}
		    }else{

			//in 3p overhang
			if((int(seq.size())-i)<=overhang3p){
#ifdef DEBUG
			    cout<<"3p"<<endl;
#endif
			    if(seq[i] == 'C'){
				if( randomProb() < sBrgs3p ){
				    deamPos.push_back(i+1);
				    seq[i] = 'T';
				    deaminated=true;
				}
			    }
			}else{ //in double strand
#ifdef DEBUG
			    cout<<"dd"<<endl;
#endif
			    if(seq[i] == 'C'){
				if( randomProb() < dBrgs ){
				    deamPos.push_back(i+1);
				    seq[i] = 'T';
				    deaminated=true;
				}
			    }
			}
		    }


		

		}//end loop each base

		if(verbose){
		    cerr<<"Read seq:"<<def<<" l="<<seq.size()<<" overhangs:"<<overhang5p<<"--"<<overhang3p;
		    		    
		    if(deaminated)
			cerr<<" added deamination "<<vectorToString(deamPos)<<endl;
		    else
			cerr<<" no deamination "<<endl;
		}
	
	    }//overhangs do not meet

	    

	}else{//not briggs single strand

	    if(useBriggs){
		int overhang5p=overhang(generator);
		int overhang3p=overhang(generator);
	    
		bool placedNick = false;
		int indexNick   = -1;
		if( (overhang5p+overhang3p)>=int(seq.size())){//all single strand

		    //just apply sBrgs over the entire fragment
		    for(int i=0;i<int(seq.size());i++){
			if(seq[i] == 'C'){
			    if( randomProb() < sBrgs ){
				deamPos.push_back(i+1);
				seq[i] = 'T';
				deaminated=true;
			    }
			}
		    }

		}else{
	
		    //Placing a nick (maybe)
		    for(int i=0;i<int(seq.size());i++){
#ifdef DEBUG
			cout<<i<<"\t"<<placedNick<<endl;
#endif
			if(!placedNick)
			    if(randomProb() < vBrgs){ //nick is present
				placedNick=true;
				indexNick =i;			
			    }


			if(placedNick){
			
			    //in 5p overhang
			    if((i+1)<=overhang5p){
#ifdef DEBUG
				cout<<"5pn"<<endl;
#endif
				if(seq[i] == 'C'){
				    if( randomProb() < sBrgs ){
					deamPos.push_back(i+1);
					seq[i] = 'T';
					deaminated=true;
				    }
				}
			    }else{

				//in 3p overhang
				if((int(seq.size())-i)<=overhang3p){
#ifdef DEBUG
				    cout<<"3pn"<<endl;
#endif
				    if(seq[i] == 'G'){
					if( randomProb() < sBrgs ){
					    deamPos.push_back(i+1);
					    seq[i] = 'A';
					    deaminated=true;
					}
				    }
				}else{ //in double strand
#ifdef DEBUG
				    cout<<"ddn"<<endl;
#endif						    
				    if(seq[i] == 'G'){
					if( randomProb() < dBrgs ){
					    deamPos.push_back(i+1);
					    seq[i] = 'A';
					    deaminated=true;
					}
				    }
				}
			    }


			}else{//No nick
			

			    //in 5p overhang
			    if((i+1)<=overhang5p){
#ifdef DEBUG
				cout<<"5p"<<endl;
#endif
				if(seq[i] == 'C'){
				    if( randomProb() < sBrgs ){
					deamPos.push_back(i+1);
					seq[i] = 'T';
					deaminated=true;
				    }
				}
			    }else{

				//in 3p overhang
				if((int(seq.size())-i)<=overhang3p){
#ifdef DEBUG
				    cout<<"3p"<<endl;
#endif
				    if(seq[i] == 'G'){
					if( randomProb() < sBrgs ){
					    deamPos.push_back(i+1);
					    seq[i] = 'A';
					    deaminated=true;
					}
				    }
				}else{ //in double strand
#ifdef DEBUG
				    cout<<"dd"<<endl;
#endif
				    if(seq[i] == 'C'){
					if( randomProb() < dBrgs ){
					    deamPos.push_back(i+1);
					    seq[i] = 'T';
					    deaminated=true;
					}
				    }
				}
			    }


		
			}//end no nick

		    }//end loop each base

		    if(verbose){
			cerr<<"Read seq:"<<def<<" l="<<seq.size()<<" overhangs:"<<overhang5p<<"--"<<overhang3p;
			if(placedNick)
			    cerr<<" nick: "<<indexNick;
		    
			if(deaminated)
			    cerr<<" added deamination "<<vectorToString(deamPos)<<endl;
			else
			    cerr<<" no deamination "<<endl;
		    }
	
		}//overhangs do not meet

	    }else{//no briggs

		if(matrixFileSpecified){
		
		    for(int i=0;i<int(seq.size());i++){
			if(!isResolvedDNA(seq[i]))
			    continue;
			double sumProb [5];
			double _sumProb[5];

			sumProb[0]  = 0.0;
			_sumProb[0] = 0.0;
		     

			//5p
			if(i<(int(seq.size())/2)){
			    int b = baseResolved2int(seq[i]);
			    double sum=0.0;
			    for(int a=0;a<4;a++){
				if(a==b)
				    continue;
				int idx;
				if(a<b)
				    idx = b*3+a;
				else
				    idx = b*3+(a-1);			     
				sum += sub5p[i].s[idx];
				_sumProb[a+1] = sub5p[i].s[idx];
			    }
			    _sumProb[b+1] = 1.0 - sum;

			    double s=0;
			    for(int a=0;a<5;a++){			     
				s += _sumProb[a];
				sumProb[a]  = s;
			    }

#ifdef DEBUG
			    for(int a=0;a<5;a++){
				cout<<"5p\t"<<i<<"\t"<<(int(seq.size())-i-1)<<"\t"<<"XACGT"[a]<<"\t"<<_sumProb[a]<<endl;
			    }
#endif
			 

			 
			    double p = randomProb();
			    //bool f=false;
			    for(int a=0;a<4;a++){
				if( (sumProb[a]    <= p ) 
				    &&
				    (sumProb[a+1]  >= p ) ){				 

				    seq[i] = "ACGT"[a];
				    //f=true;
				    if(a!=b){
					deamPos.push_back(i+1);
					deaminated=true;
				    }

				    break;
				}
			    }

			    // if(b==1 && 
			    //    i==0){
			    //     for(int a=0;a<5;a++){
			    // 	 cout<<"test\t5p\t"<<i<<"\t"<<(int(seq.size())-i-1)<<"\t"<<"XACGT"[a]<<"\t"<<sumProb[a]<<"\t"<<p<<"\t"<<vectorToString(deamPos)<<endl;
			    //     }
			    //     cout<<"deam\t5p\t"<<i<<"\t"<<(int(seq.size())-i-1)<<"\t"<<seq[i]<<"\t"<<p<<"\t"<<vectorToString(deamPos)<<endl;			     
			    // }


			    // if(b==1){
			    //     for(int a=0;a<5;a++){
			    // 	 cout<<"5p\t"<<i<<"\t"<<(int(seq.size())-i-1)<<"\t"<<"ACGT"[b]<<"\t"<<seq[i]<<endl;
			    //     }
			    // }
			 			 
			}
			//3p
			else{

			    int b = baseResolved2int(seq[i]);
			    double sum=0.0;
			    for(int a=0;a<4;a++){
				if(a==b)
				    continue;
				int idx;
				if(a<b)
				    idx = b*3+a;
				else
				    idx = b*3+(a-1);			     
				sum += sub3p[ int(seq.size())-i-1 ].s[idx];
				_sumProb[a+1] = sub3p[ int(seq.size())-i-1 ].s[idx];
			    }
			    _sumProb[b+1] = 1.0 - sum;


			    double s=0;
			    for(int a=0;a<5;a++){			     
				s += _sumProb[a];
				sumProb[a]  = s;
			    }
			 
#ifdef DEBUG
			    for(int a=0;a<5;a++){
				cout<<"3p\t"<<i<<"\t"<<(int(seq.size())-i-1)<<"\t3p\t"<<"XACGT"[a]<<"\t"<<sumProb[a]<<endl;
			    }
#endif


			 
			    double p = randomProb();
			    //bool f=false;
			    for(int a=0;a<4;a++){
				if( (sumProb[a]    <= p ) 
				    &&
				    (sumProb[a+1]  >= p ) ){				 
				    seq[i] = "ACGT"[a];
				    if(a!=b){
					deamPos.push_back( -(int(seq.size())-i ) );
					deaminated=true;
				    }
				    //f=true;
				    break;
				}
			    }

			    // if(b==1 && 
			    //    i==(int(seq.size())-1) ){

			    //     for(int a=0;a<5;a++){
			    // 	 cout<<"test\t3p\t"<<i<<"\t"<<(int(seq.size())-i-1)<<"\t"<<"XACGT"[a]<<"\t"<<sumProb[a]<<"\t"<<p<<"\t"<<vectorToString(deamPos)<<endl;
			    //     }

			    //     cout<<"deam\t3p\t"<<i<<"\t"<<(int(seq.size())-i-1)<<"\t"<<seq[i]<<"\t"<<p<<"\t"<<vectorToString(deamPos)<<endl; 	     
			    // }


			}
		    }

		}else{ //no matrix specified

		    //MATRIX FILE
		    if(singleDeam){ //Single strand will have C->T damage on both ends
	    
			//5p
			for(int i=0;
			    i<int(seq.size());
			    i++){
			    if(seq[i] == 'C'){
				if( randomProb() < sub5p[i].s[5] ){
				    deamPos.push_back(i+1);
				    seq[i] = 'T';
				    deaminated=true;
				}
			    }
			}

			//3p
			for(int i=0;
			    i<int(seq.size());
			    i++){
			    if(seq[i] == 'C'){
				if( randomProb() < sub3p[ int(seq.size())-i-1 ].s[5] ){		
				    deamPos.push_back( -(int(seq.size())-i));
				    seq[i] = 'T';
				    deaminated=true;
				}
			    }
			}



		    }

		    if(doubleDeam){ //Double strand will have and C->T at the 5p end and G->A damage at the 3p end
	    
			//5p
			for(int i=0;
			    i<int(seq.size());
			    i++){
			    if(seq[i] == 'C'){
				if( randomProb() < sub5p[i].s[5] ){
				    deamPos.push_back(i+1);
				    seq[i] = 'T';
				    deaminated=true;
				}
			    }
			}

			//3p
			for(int i=0;
			    i<int(seq.size());
			    i++){
			    if(seq[i] == 'G'){
				if( randomProb() < sub3p[ int(seq.size())-i-1 ].s[6] ){			
				    deamPos.push_back( -(int(seq.size())-i));
				    seq[i] = 'A';
				    deaminated=true;
				}

			    }
			}

		    }
		}//closes else no matrix specified

		if(verbose){
		    cerr<<"Read seq:"<<def;
		    if(deaminated)
			cerr<<" added deamination "<<vectorToString(deamPos)<<endl;
		    else
			cerr<<" no deamination "<<endl;
		}

	    }//closes else no briggs
	}//closes else no briggsss

	if(outBAMb){
	    al.QueryBases  = seq;
	    if(putDeamInName && deaminated){
		if(!al.AddTag("XD", "Z",vectorToString(deamPos)) ) {
		    cerr<<"Internal error, cannot add tag to BAM alignment "<<al.Name<<endl; 
		    return 1; 
		}	
	    }
	    writer.SaveAlignment(al);

	}else{

	    if(outFastagzb){
		if(putDeamInName && deaminated){
		    outFastagzfp<<*(fo->getID())<<"_DEAM:"<<vectorToString(deamPos)<<endl<<seq<<endl;
		}else{
		    outFastagzfp<<*(fo->getID())<<endl<<seq<<endl;
		}
	    }else{
		if(putDeamInName && deaminated){
		    cout<<*(fo->getID())<<"_DEAM:"<<vectorToString(deamPos)<<endl<<seq<<endl;
		}else{
		    cout<<*(fo->getID())<<endl<<seq<<endl;
		}
	    }
	}
    }

    if(outFastagzb){
	outFastagzfp.close();       
    }
    
    if(outBAMb){
        reader.Close();
        writer.Close();
    }

    cerr<<"Program "<<argv[0]<<" terminated succesfully, wrote "<<f<<" sequences"<<endl;

    return 0;
}


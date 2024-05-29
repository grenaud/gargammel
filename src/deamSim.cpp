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
    bool lastRowMatrix         =true;
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

    vector<substitutionRates> sub5pWiMeth;
    vector<substitutionRates> sub3pWiMeth;
    vector<substitutionRates> sub5pNoMeth;
    vector<substitutionRates> sub3pNoMeth;

    string profFile;
    bool   profFileSpecified=false;

    string matrixFile;
    bool   matrixFileSpecified=false;

    string matrixFileNoMeth;
    bool   matrixFileSpecifiedNoMeth  =false;
    string matrixFileWiMeth;
    bool   matrixFileSpecifiedWithMeth=false;

    string mapdamageFile;
    bool   mapdamageFileSpecified=false;
    string mapdamageProtocol;

    string matrix;
    bool   matrixSpecified    = false;
    bool   singleDeam         = false;
    bool   doubleDeam         = false;
    
    bool seedSpecified=false;
    unsigned int seed=0;
    
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
			      "\t\t"+"-name"              +"\t\t\t\t\t"+"Put a tag in the read name with deam bases (default "+booleanAsString(putDeamInName)+")"+"\n"+
			      "\t\t"+"-v\t"+""            +"\t\t\t\t"+"verbose mode"+"\n"+
			      "\t\t"+"-last\t"+""            +"\t\t\t\t"+"If matfile is used, do not use the substitution rates of the"+"\n"+
			      "\t\t"+"\t"+""             +"\t\t\t\t"+"last row over the rest of the molecule (default: no data = use last row)"+"\n"+
			      " Generic options:\n"+
			      "\t\t"+"--seed\t"+"[int]"       +"\t\t\t\t"+"Use [seed] as seed for the random number generator (default random seed each execution)"+"\n"+

			      
			      "\n"
			      +" Mandatory deamination options:\n"+ 
			      "\tSpecify either:\n"+
			      
                              "\t\t"+"-mapdamage  [mis.txt] [protocol]" +"\t"+"Read the miscorporation file [mis.txt]"+"\n"+
                              "\t\t"+"                               " +"\t\t"+"produced by mapDamage"+"\n"+
                              "\t\t"+"                               " +"\t\t"+"[protocol] can be either \"single\" or \"double\" (without quotes)\n"+
                              "\t\t"+"                               " +"\t\t"+"Single strand will have C->T damage on both ends"+"\n"+
                              "\t\t"+"                               " +"\t\t"+"Double strand will have and C->T at the 5' end and G->A damage at the 3' end"+"\n"+
			      "\n"+

                              "\t\t"+"-profile  [prof file prefix]" +"\t\t"+"Read the profile of substitutions instead of the default"+"\n"+
                              "\t\t"+"                               " +"\t\t"+"Provide the prefix only, both files must end with"+"\n"+
                              "\t\t"+"                               " +"\t\t"+"5p.prof and 3p.prof"+"\n"+
			      "\n"+

                              "\t\t"+"-matfile  [matrix file prefix]" +"\t\t"+"Read the matrix file of substitutions instead of the default"+"\n"+
                              "\t\t"+"                               " +"\t\t"+"Provide the prefix only, both files must end with"+"\n"+
                              "\t\t"+"                               " +"\t\t"+"5.dat and 3.dat"+"\n"+
			      "\n"+

                              "\t\t"+"-matfilenonmeth  [matrix file prefix]" +"\t"+"Read the matrix file of substitutions for non-methylated Cs"+"\n"+
                              "\t\t"+"                               " +"\t\t"+"Provide the prefix only, both files must end with"+"\n"+
                              "\t\t"+"                               " +"\t\t"+"5.dat and 3.dat"+"\n"+
			      "\n"+

                              "\t\t"+"-matfilemeth  [matrix file prefix]" +"\t"+"Read the matrix file of substitutions for methylated Cs"+"\n"+
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

	if(string(argv[i]) == "--seed" ){
	    seed=destringify<unsigned int>(argv[i+1]);
	    i++;
	    seedSpecified=true;
	    continue;
	}
	
	if(string(argv[i]) == "-profile" ){
	    profFile          = string(argv[i+1]);
	    profFileSpecified = true;
	    i++;
	    continue;
	}

	if(string(argv[i]) == "-matfile" ){
	    matrixFile          = string(argv[i+1]);
	    matrixFileSpecified = true;
	    i++;
	    continue;
	}

	if(string(argv[i]) == "-matfilenonmeth" ){
	    matrixFileNoMeth          = string(argv[i+1]);
	    matrixFileSpecifiedNoMeth = true;
	    i++;
	    continue;
	}

	if(string(argv[i]) == "-matfilemeth" ){
	    matrixFileWiMeth            = string(argv[i+1]);
	    matrixFileSpecifiedWithMeth = true;
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

        if(string(argv[i]) == "-last"  ){
            lastRowMatrix=false;
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

    if(!useBriggs                     && 
       !useBriggsss                   && 
       !matrixSpecified               && 
       !matrixFileSpecifiedWithMeth  &&
       !matrixFileSpecifiedNoMeth     &&
       !matrixFileSpecified           &&
       !mapdamageFileSpecified &&
       !profFileSpecified){
    	cerr << "Please specify the matrix to use or use the Briggs parameter model"<<endl;
    	return 1;
    }

    if(profFileSpecified && (matrixFileSpecified          |
			     matrixFileSpecifiedWithMeth  ||
			     matrixFileSpecifiedNoMeth  )  ){
    	cerr << "Please specify damage profile or substitution matrix"<<endl;
    	return 1;	
    }
    
    if(matrixFileSpecified          &&
       matrixFileSpecifiedWithMeth  &&
       matrixFileSpecifiedNoMeth    ){
    	cerr << "Please specify matrix with/without methylation but not the overall matrix"<<endl;
    	return 1;	
    }

    if( 
       (!matrixFileSpecifiedWithMeth  &&
	 matrixFileSpecifiedNoMeth    )
	||
       ( matrixFileSpecifiedWithMeth  &&
	 !matrixFileSpecifiedNoMeth    )
	){
    	cerr << "Please specify both the matrix with/without methylation, not just one"<<endl;
    	return 1;	
    }

    if(outFastagzb && outBAMb){ 
        cerr<<"Error: cannot specify both -o and -b"<<endl;
        return 1;
    }
    string inputFile = string(argv[argc-1]);


    if(seedSpecified){
	srand(   seed );
	srand48( seed );
	generator.seed(seed);
    }else{
    }
    

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
		cerr<<i<<"\t";
		for(unsigned int k=0;k<12;k++)
		    cerr<<sub5p[i].s[k]<<" ";
		cerr<<endl;
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
		cerr<<i<<"\t";
		for(unsigned int k=0;k<12;k++)
		    cerr<<sub3p[i].s[k]<<" ";
		cerr<<endl;
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
	 
	    //end if mapdamageFileSpecified   
	}else{//else, we specify a .dat matrix

	    //cerr<<"test1"<<endl;
	    if(matrixFileSpecifiedWithMeth || matrixFileSpecifiedNoMeth){
		//With methylation
		
		string deam5FileWiMeth;
		string deam3FileWiMeth;
		
		deam5FileWiMeth = matrixFileWiMeth+"5.dat";
		deam3FileWiMeth = matrixFileWiMeth+"3.dat";	
		string line;
#ifdef DEBUG
		cerr<<"Files deamination with methylation"<<endl;
		cerr<<deam5FileWiMeth<<endl;
		cerr<<deam3FileWiMeth<<endl;
#endif
		////////////////////////////////////////////////////////
		//READ file with 5p deamination rates with methylation/
		///////////////////////////////////////////////////////

		igzstream deam5pFileStWiMeth;

		deam5pFileStWiMeth.open(deam5FileWiMeth.c_str(), ios::in);


		if (deam5pFileStWiMeth.good()){
		    //header
		    if ( !getline (deam5pFileStWiMeth,line)){
			cerr << "Parsing error for 5p deamination file "<<deam5FileWiMeth<<endl;
			return 1;
		    }

		    while ( getline (deam5pFileStWiMeth,line)){

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
			    cerr<<sub5pWiMeth.size()<<"\t"<<substiString[i]<<"\t"<<toadd.s[i]<<"\t"<<firstFields[i]<<endl;
#endif
			}	


			for(int b=0;b<4;b++){		
			    double sum =0.0;
			    for(int a=0;a<3;a++){		
				int idx = b*3+a;
				sum+=toadd.s[idx];
			    }
			    if(sum>1.0){
				cerr<<"Problem with line "<<line<<" in file "<<deam5FileWiMeth<<" the sum of the substitution probabilities exceeds 1"<<endl;
				return 1;
			    }
			}

			sub5pWiMeth.push_back(toadd);
			//return 1;
		    }
	
		    deam5pFileStWiMeth.close();
		}else{
		    cerr << "Unable to open 5p deamination file "<<deam5FileWiMeth<<endl;
		    return 1;
		}


		////////////////////////////////////////////////////////
		//READ file with 3p deamination rates with methylation/
		///////////////////////////////////////////////////////
		igzstream deam3pFileStWiMeth;

		deam3pFileStWiMeth.open(deam3FileWiMeth.c_str(), ios::in);


		if (deam3pFileStWiMeth.good()){
		    //header
		    if ( !getline (deam3pFileStWiMeth,line)){
			cerr << "Parsing error for 3p deamination file "<<deam3FileWiMeth<<endl;
			return 1;
		    }

		    while ( getline (deam3pFileStWiMeth,line)){
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
			    cerr<<sub3pWiMeth.size()<<"\t"<<substiString[i]<<"\t"<<toadd.s[i]<<"\t"<<firstFields[i]<<endl;
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
				cerr<<"Problem with line "<<line<<" in file "<<deam3FileWiMeth<<" the sum of the substitution probabilities exceeds 1"<<endl;
				return 1;
			    }
			}
		

			sub3pWiMeth.push_back(toadd);
			//return 1;
		    }
	
		    deam3pFileStWiMeth.close();
		}else{
		    cerr << "Unable to open 3p deamination file "<<deam3FileWiMeth<<endl;
		    return 1;
		}

		reverse(sub3pWiMeth.begin(),sub3pWiMeth.end());



		//No methylation

		string deam5FileNoMeth;
		string deam3FileNoMeth;
		
		deam5FileNoMeth = matrixFileNoMeth+"5.dat";
		deam3FileNoMeth = matrixFileNoMeth+"3.dat";	

#ifdef DEBUG
		cerr<<"Files deamination without methylation"<<endl;
		cerr<<deam5FileNoMeth<<endl;
		cerr<<deam3FileNoMeth<<endl;
#endif

		//////////////////////////////////////////////////////////
		//READ file with 5p deamination rates without methylation/
		//////////////////////////////////////////////////////////
		igzstream deam5pFileStNoMeth;

		deam5pFileStNoMeth.open(deam5FileNoMeth.c_str(), ios::in);


		if (deam5pFileStNoMeth.good()){
		    //header
		    if ( !getline (deam5pFileStNoMeth,line)){
			cerr << "Parsing error for 5p deamination file "<<deam5FileNoMeth<<endl;
			return 1;
		    }

		    while ( getline (deam5pFileStNoMeth,line)){
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
			    cerr<<sub5pNoMeth.size()<<"\t"<<substiString[i]<<"\t"<<toadd.s[i]<<"\t"<<firstFields[i]<<endl;
#endif
			}	


			for(int b=0;b<4;b++){		
			    double sum =0.0;
			    for(int a=0;a<3;a++){		
				int idx = b*3+a;
				sum+=toadd.s[idx];
			    }
			    if(sum>1.0){
				cerr<<"Problem with line "<<line<<" in file "<<deam5FileNoMeth<<" the sum of the substitution probabilities exceeds 1"<<endl;
				return 1;
			    }
			}

			sub5pNoMeth.push_back(toadd);
			//return 1;
		    }
	
		    deam5pFileStNoMeth.close();
		}else{
		    cerr << "Unable to open 5p deamination file "<<deam5FileNoMeth<<endl;
		    return 1;
		}

		//////////////////////////////////////////////////////////
		//READ file with 3p deamination rates without methylation/
		//////////////////////////////////////////////////////////
		igzstream deam3pFileStNoMeth;

		deam3pFileStNoMeth.open(deam3FileNoMeth.c_str(), ios::in);


		if (deam3pFileStNoMeth.good()){
		    //header
		    if ( !getline (deam3pFileStNoMeth,line)){
			cerr << "Parsing error for 3p deamination file "<<deam3FileNoMeth<<endl;
			return 1;
		    }

		    while ( getline (deam3pFileStNoMeth,line)){
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
			    cerr<<sub3pNoMeth.size()<<"\t"<<substiString[i]<<"\t"<<toadd.s[i]<<"\t"<<firstFields[i]<<endl;
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
				cerr<<"Problem with line "<<line<<" in file "<<deam3FileNoMeth<<" the sum of the substitution probabilities exceeds 1"<<endl;
				return 1;
			    }
			}
		

			sub3pNoMeth.push_back(toadd);
			//return 1;
		    }
	
		    deam3pFileStNoMeth.close();
		}else{
		    cerr << "Unable to open 3p deamination file "<<deam3FileNoMeth<<endl;
		    return 1;
		}

		reverse(sub3pNoMeth.begin(),sub3pNoMeth.end());
		

		//end 	    if(matrixFileSpecifiedWithMeth || matrixFileSpecifiedNoMeth){
	    }//specified methylation
	    else{
		if(profFileSpecified){
		    string deam5File;
		    string deam3File;
	
		    deam5File = profFile+"5p.prof";
		    deam3File = profFile+"3p.prof";	
		
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
			    for(unsigned int i=0;i<fields.size();i++){
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
				cerr<<sub5p.size()<<"\t"<<substiString[i]<<"\t"<<toadd.s[i]<<"\t"<<firstFields[i]<<endl;
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
			    for(unsigned int i=0;i<fields.size();i++){
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
				cerr<<sub3p.size()<<"\t"<<substiString[i]<<"\t"<<toadd.s[i]<<"\t"<<firstFields[i]<<endl;
#endif

				//cerr<<toadd.s[i]<<"\t"<<firstFields[i]<<endl;
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
		    // //cerr<<"test2"<<endl;
		    // for(unsigned int i=0;i<sub3p.size();i++){		
		    //     cerr<<i<<"\t"<<arrayToString(sub3p[i].s,12,"\t")<<endl;
		    // }
		    //reverse(sub3p.begin(),sub3p.end());
	
		    // for(unsigned int i=0;i<sub3p.size();i++){		
		    //     cerr<<i<<"\t"<<arrayToString(sub3p[i].s,12,"\t")<<endl;
		    // }
		}else{

		    
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
				cerr<<sub5p.size()<<"\t"<<substiString[i]<<"\t"<<toadd.s[i]<<"\t"<<firstFields[i]<<endl;
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
				cerr<<sub3p.size()<<"\t"<<substiString[i]<<"\t"<<toadd.s[i]<<"\t"<<firstFields[i]<<endl;
#endif

				//cerr<<toadd.s[i]<<"\t"<<firstFields[i]<<endl;
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
		    // //cerr<<"test2"<<endl;
		    // for(unsigned int i=0;i<sub3p.size();i++){		
		    //     cerr<<i<<"\t"<<arrayToString(sub3p[i].s,12,"\t")<<endl;
		    // }
		    reverse(sub3p.begin(),sub3p.end());
	
		    // for(unsigned int i=0;i<sub3p.size();i++){		
		    //     cerr<<i<<"\t"<<arrayToString(sub3p[i].s,12,"\t")<<endl;
		    // }
		}
	    }
	}

	//end if(!useBriggs && !useBriggsss){    
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

	if(!(matrixFileSpecifiedWithMeth || matrixFileSpecifiedNoMeth)){	
	    transform(seq.begin(), seq.end(),seq.begin(), ::toupper);
	}

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
#ifdef DEBUG
		cerr<<"Briggs single-strand"<<endl;
#endif

	    int overhang5p=overhanggeo5p(generator);
	    int overhang3p=overhanggeo3p(generator);

	    if( (overhang5p+overhang3p)>=int(seq.size())){//all single strand

		//just apply sBrgs over the entire fragment
		for(int i=0;i<int(seq.size());i++){

		    //both in 5p and 3p overhang, pick closest
		    if( ((i+1)<=overhang5p) && ((int(seq.size())-i)<=overhang3p) ){
			if(i<int(double(seq.size())*0.5)){
			    if(seq[i] == 'C'){
				if( randomProb(!seedSpecified) < sBrgs5p ){
				    deamPos.push_back(i+1);
				    seq[i] = 'T';
				    deaminated=true;
				}
			    }
			}else{
			    if(seq[i] == 'C'){
				if( randomProb(!seedSpecified) < sBrgs3p ){
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
				if( randomProb(!seedSpecified) < sBrgs3p ){
				    deamPos.push_back(i+1);
				    seq[i] = 'T';
				    deaminated=true;
				}
			    }

			}else{

			    //just 5' overhang
			    if(  ((i+1)<=overhang5p) && !((int(seq.size())-i)<=overhang3p) ){

				if(seq[i] == 'C'){
				    if( randomProb(!seedSpecified) < sBrgs5p ){
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
			cerr<<"5p"<<endl;
#endif
			if(seq[i] == 'C'){
			    if( randomProb(!seedSpecified) < sBrgs5p ){
				deamPos.push_back(i+1);
				seq[i] = 'T';
				deaminated=true;
			    }
			}
		    }else{

			//in 3p overhang
			if((int(seq.size())-i)<=overhang3p){
#ifdef DEBUG
			    cerr<<"3p"<<endl;
#endif
			    if(seq[i] == 'C'){
				if( randomProb(!seedSpecified) < sBrgs3p ){
				    deamPos.push_back(i+1);
				    seq[i] = 'T';
				    deaminated=true;
				}
			    }
			}else{ //in double strand
#ifdef DEBUG
			    cerr<<"dd"<<endl;
#endif
			    if(seq[i] == 'C'){
				if( randomProb(!seedSpecified) < dBrgs ){
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

#ifdef DEBUG
		cerr<<"Briggs double-strand"<<endl;
#endif

	    if(useBriggs){
		int overhang5p=overhang(generator);
		int overhang3p=overhang(generator);
	    
		bool placedNick = false;
		int indexNick   = -1;
		if( (overhang5p+overhang3p)>=int(seq.size())){//all single strand

		    //just apply sBrgs over the entire fragment
		    for(int i=0;i<int(seq.size());i++){
			if(seq[i] == 'C'){
			    if( randomProb(!seedSpecified) < sBrgs ){
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
			cerr<<i<<"\t"<<placedNick<<endl;
#endif
			if(!placedNick)
			    if(randomProb(!seedSpecified) < vBrgs){ //nick is present
				placedNick=true;
				indexNick =i;			
			    }


			if(placedNick){
			
			    //in 5p overhang
			    if((i+1)<=overhang5p){
#ifdef DEBUG
				cerr<<"5pn"<<endl;
#endif
				if(seq[i] == 'C'){
				    if( randomProb(!seedSpecified) < sBrgs ){
					deamPos.push_back(i+1);
					seq[i] = 'T';
					deaminated=true;
				    }
				}
			    }else{

				//in 3p overhang
				if((int(seq.size())-i)<=overhang3p){
#ifdef DEBUG
				    cerr<<"3pn"<<endl;
#endif
				    if(seq[i] == 'G'){
					if( randomProb(!seedSpecified) < sBrgs ){
					    deamPos.push_back(i+1);
					    seq[i] = 'A';
					    deaminated=true;
					}
				    }
				}else{ //in double strand
#ifdef DEBUG
				    cerr<<"ddn"<<endl;
#endif						    
				    if(seq[i] == 'G'){
					if( randomProb(!seedSpecified) < dBrgs ){
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
				cerr<<"5p"<<endl;
#endif
				if(seq[i] == 'C'){
				    if( randomProb(!seedSpecified) < sBrgs ){
					deamPos.push_back(i+1);
					seq[i] = 'T';
					deaminated=true;
				    }
				}
			    }else{

				//in 3p overhang
				if((int(seq.size())-i)<=overhang3p){
#ifdef DEBUG
				    cerr<<"3p"<<endl;
#endif
				    if(seq[i] == 'G'){
					if( randomProb(!seedSpecified) < sBrgs ){
					    deamPos.push_back(i+1);
					    seq[i] = 'A';
					    deaminated=true;
					}
				    }
				}else{ //in double strand
#ifdef DEBUG
				    cerr<<"dd"<<endl;
#endif
				    if(seq[i] == 'C'){
					if( randomProb(!seedSpecified) < dBrgs ){
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

#ifdef DEBUG
		cerr<<"profFileSpecified         "<<profFileSpecified<<endl;
		cerr<<"matrixFileSpecified         "<<matrixFileSpecified<<endl;
		cerr<<"matrixFileSpecifiedNoMeth   "<<matrixFileSpecifiedNoMeth<<endl;
		cerr<<"matrixFileSpecifiedWithMeth "<<matrixFileSpecifiedWithMeth<<endl;

#endif

		if(profFileSpecified           ||
		   matrixFileSpecified         ||
		   matrixFileSpecifiedNoMeth   ||
		   matrixFileSpecifiedWithMeth ){
		
		    if(matrixFileSpecifiedNoMeth || matrixFileSpecifiedWithMeth ){
			//Methylation is taken into account

#ifdef DEBUG
			cerr<<"Methylation lastRowMatrix="<<lastRowMatrix<<endl;
#endif

			for(int i=0;i<int(seq.size());i++){
			    char originalBase = toupper(seq[i]);
			    bool isMethylated = (originalBase != seq[i]);
			    //cerr<<originalBase<<"\t"<<seq[i]<<endl;
			    vector<substitutionRates> * sub5pToUse;
			    vector<substitutionRates> * sub3pToUse;

			    if(isMethylated){
				sub5pToUse = &sub5pWiMeth ;
				sub3pToUse = &sub3pWiMeth ;
			    }else{
				sub5pToUse = &sub5pNoMeth ;
				sub3pToUse = &sub3pNoMeth ;
			    }
			    
			    if(!isResolvedDNA(originalBase))
				continue;

			    double sumProb [5];
			    double _sumProb[5];

			    sumProb[0]  = 0.0;
			    _sumProb[0] = 0.0;
		     

			    //5p
			    if(i<(int(seq.size())/2)){
				int b = baseResolved2int(originalBase);
				double sum=0.0;
				
				int iToUseSub5 = i;
				if(iToUseSub5>=int(sub5pToUse->size())){
				    if(lastRowMatrix){
					iToUseSub5 = int(sub5pToUse->size())-1;
				    }else{
					continue;
				    }
				}

				for(int a=0;a<4;a++){
				    if(a==b)
					continue;
				    int idx;
				    if(a<b)
					idx = b*3+a;
				    else
					idx = b*3+(a-1);			     
				    sum           += sub5pToUse->at(iToUseSub5).s[idx];
				    _sumProb[a+1]  = sub5pToUse->at(iToUseSub5).s[idx];
				}
				_sumProb[b+1] = 1.0 - sum;

				double s=0;
				for(int a=0;a<5;a++){			     
				    s += _sumProb[a];
				    sumProb[a]  = s;
				}

#ifdef DEBUG
				for(int a=0;a<5;a++){
				    cerr<<"M5p\t"<<i<<"\t"<<(int(seq.size())-i-1)<<"\t"<<iToUseSub5<<"\t"<<"XACGT"[a]<<"\t"<<_sumProb[a]<<"\t"<<isMethylated<<endl;
				}
#endif
			 

			 
				double p = randomProb(!seedSpecified);
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

			 			 
			    }
			    //3p
			    else{

				int b = baseResolved2int(originalBase);
				double sum=0.0;

				int iToUseSub3 = int(seq.size())-i-1;
				if(iToUseSub3>=int(sub3pToUse->size())){
				    if(lastRowMatrix){
					iToUseSub3 = int(sub3pToUse->size())-1;
				    }else{
					continue;
				    }
				}


				for(int a=0;a<4;a++){
				    if(a==b)
					continue;
				    int idx;
				    if(a<b)
					idx = b*3+a;
				    else
					idx = b*3+(a-1);			     
				    sum           += sub3pToUse->at( iToUseSub3 ).s[idx];
				    _sumProb[a+1]  = sub3pToUse->at( iToUseSub3 ).s[idx];
				}
				_sumProb[b+1] = 1.0 - sum;


				double s=0;
				for(int a=0;a<5;a++){			     
				    s += _sumProb[a];
				    sumProb[a]  = s;
				}
			 
#ifdef DEBUG
				for(int a=0;a<5;a++){
				    cerr<<"M3p\t"<<i<<"\t"<<(int(seq.size())-i-1)<<"\t"<<iToUseSub3<<"\t3p\t"<<"XACGT"[a]<<"\t"<<sumProb[a]<<"\t"<<isMethylated<<endl;
				}
#endif


			 
				double p = randomProb(!seedSpecified);
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



			    }
			}//end for
		    }

		    //else, we have a simple matrix no methylation
		    else{
#ifdef DEBUG
 			cerr<<"unmethylated"<<endl;
#endif
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
				
				
				int iToUseSub5 = i;
				if(iToUseSub5>=int(sub5p.size())){
				    if(lastRowMatrix){
					iToUseSub5 = int(sub5p.size())-1;
				    }else{
					continue;
				    }
				}


				for(int a=0;a<4;a++){
				    if(a==b)
					continue;
				    int idx;
				    if(a<b)
					idx = b*3+a;
				    else
					idx = b*3+(a-1);			     
				    sum += sub5p[iToUseSub5].s[idx];
				    _sumProb[a+1] = sub5p[iToUseSub5].s[idx];
				}
				_sumProb[b+1] = 1.0 - sum;

				double s=0;
				for(int a=0;a<5;a++){			     
				    s += _sumProb[a];
				    sumProb[a]  = s;
				}

#ifdef DEBUG
				for(int a=0;a<5;a++){
				    cerr<<"5p\t"<<i<<"\t"<<(int(seq.size())-i-1)<<"\t"<<iToUseSub5<<"\t"<<"XACGT"[a]<<"\t"<<_sumProb[a]<<endl;
				}
#endif
			 

			 
				double p = randomProb(!seedSpecified);
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

			 			 
			    }
			    //3p
			    else{

				int b = baseResolved2int(seq[i]);
				double sum=0.0;

				int iToUseSub3 = int(seq.size())-i-1;
				if(iToUseSub3>=int(sub3p.size())){
				    if(lastRowMatrix){
					iToUseSub3 = int(sub3p.size())-1;
				    }else{
					continue;
				    }
				}


				for(int a=0;a<4;a++){
				    if(a==b)
					continue;
				    int idx;
				    if(a<b)
					idx = b*3+a;
				    else
					idx = b*3+(a-1);			     
				    sum += sub3p[ iToUseSub3 ].s[idx];
				    _sumProb[a+1] = sub3p[ iToUseSub3 ].s[idx];
				}
				_sumProb[b+1] = 1.0 - sum;


				double s=0;
				for(int a=0;a<5;a++){			     
				    s += _sumProb[a];
				    sumProb[a]  = s;
				}
			 
#ifdef DEBUG
				for(int a=0;a<5;a++){
				    cerr<<"3p\t"<<i<<"\t"<<(int(seq.size())-i-1)<<"\t"<<(iToUseSub3)<<"\t3p\t"<<"XACGT"[a]<<"\t"<<sumProb[a]<<endl;
				}
#endif


			 
				double p = randomProb(!seedSpecified);
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



			    }
			}//end for
		    }//end if matrixFileSpecified

		    //no matrix specified
		}else{ 

#ifdef DEBUG
		    cerr<<"matrix file, singleDeam="<<singleDeam<<" doubleDeam="<<doubleDeam<<endl;
#endif

		    //MATRIX FILE
		    if(singleDeam){ //Single strand will have C->T damage on both ends
	    
			//5p
			for(int i=0;
			    i<int(seq.size());
			    i++){
			    if(seq[i] == 'C'){
				if( randomProb(!seedSpecified) < sub5p[i].s[5] ){
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
				if( randomProb(!seedSpecified) < sub3p[ int(seq.size())-i-1 ].s[5] ){		
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
				if( randomProb(!seedSpecified) < sub5p[i].s[5] ){
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
				if( randomProb(!seedSpecified) < sub3p[ int(seq.size())-i-1 ].s[6] ){			
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


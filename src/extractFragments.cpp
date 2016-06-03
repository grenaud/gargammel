/*
 * extractFragments
 * Date: Jun-02-2016 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>
#include <gzstream.h>
#include <random>

//#define DEBUG

#include "utils.h"

using namespace std;

typedef struct{
    double s[4];
} subrates;



bool isResolvedDNAstring(const string & s){
    for(unsigned int i=0;i<s.size();i++){
	if(!isResolvedDNA(s[i]))
	    return false;
    }
    return true;
}

int main (int argc, char *argv[]) {


    bool specifiedLoc    = false;
    bool specifiedScale  = false;
    bool specifiedLength = false;

    //log normal
    double location           = -1;
    double scale              = -1;
    default_random_engine generator;
    lognormal_distribution<double> distribution;    

    int   sizeFragments       = 20;
    unsigned int   nFragments = 10;
    string compFile           = "";
    bool  compFileSpecified   = false;
    string line;
    int    distFromEnd        = 1;
    
    vector<subrates> sub5pPlus;
    vector<subrates> sub5pMinus;
    vector<subrates> sub3pPlus;
    vector<subrates> sub3pMinus;

    string fileSizeFrag;
    bool   fileSizeFragB=false;
    int minimumFragSize=   0;
    int maximumFragSize=1000;


    const string usage=string("\nThis program takes a fasta file representing a chromosome and generates\n")+
	"aDNA fragments according to a certain distribution\n\n"+
	string(argv[0])+" [options]  [chromosome fasta] "+"\n\n"+
	
	"\t\t"+"-n\t"+"[number]" +"\t\t\t"+"Generate [number] fragments (default: "+stringify(nFragments)+")"+"\n"+
	"\n"+
	"\t\t"+"--comp\t"+"[file]"+"\t\t\t\t"+"Base composition for the fragments (default none)"+"\n"+
	"\t\t"+"--dist\t"+"[file]"+"\t\t\t\t"+"Distance from ends to consider  (default: "+stringify(distFromEnd)+")"+"\n"+
	"\t\t"+"\t"+""+"\t\t\t\t"+"if this is not specified, the base composition"+"\n"+
	"\t\t"+"\t"+""+"\t\t\t\t"+"will only reflect the chromosome file used"+"\n"+
	"\n"+
	"\tFragment size: \n"+
	"\t\t"+"-m\t"+"[length]" +"\t\t\t"+"Minimum fragments length  (default: "+stringify(minimumFragSize)+")"+"\n"+
	"\t\t"+"-M\t"+"[length]" +"\t\t\t"+"Maximum fragments length  (default: "+stringify(maximumFragSize)+")"+"\n"+

	"\tFragment size distribution: specify either one of the 3 possible options\n"+
	"\t\t"+"-l\t"+"[length]" +"\t\t\t"+"Generate fragments of fixed length  (default: "+stringify(sizeFragments)+")"+"\n"+
	"\t\t"+"-s\t"+"[file]"   +"\t\t\t\t"+"Open file with size distribution"+"\n"+
	"\n\t\tLength options:\n"+                                
	"\t\t\t"+"--loc\t"+"[file]"   +  "\t\t\t"+"Location for lognormal distribution (default none)"+"\n"+
	"\t\t\t"+"--scale\t"+"[file]"   +"\t\t\t"+"Scale for lognormal distribution      (default none)"+"\n"+
	
	
	// "\n\tOutput options:\n"+  
	// "\t\t"+"-seq  [fasta file]" +"\t\t"+"Output fasta file (default: stdout)"+"\n"+
	// "\t\t"+"-log  [log file]" +"\t\t"+"Output log  (default: stderr)"+"\n"+
	// "\t\t"+"-name [name]" +"\t\t\t"  +"Name  (default "+nameMT+") "+"\n"+
	// "\t\t"+"-qual [minimum quality]" +"\t\t"  +"Filter bases with quality less than this  (default "+stringify(minQ
	
	"";



    if( (argc== 1) ||
        (argc== 2 && string(argv[1]) == "-h") ||
        (argc== 2 && string(argv[1]) == "-help") ||
        (argc== 2 && string(argv[1]) == "--help") ){
        cout<<usage<<endl;
        return 1;
    }

    for(int i=1;i<(argc-1);i++){ //all but the last 3 args

        if(string(argv[i]) == "--comp" ){
	    compFile          = string(argv[i+1]);
            i++; 
	    compFileSpecified = true;
	    // cout<<"cf "<<compFile<<endl;
            continue;
        }

        if(string(argv[i]) == "--dist" ){
	    distFromEnd=destringify<int>(argv[i+1]);
            i++; 
            continue;
        }

	if(string(argv[i]) == "--loc" ){
            location=destringify<double>(argv[i+1]);
            i++;
            specifiedLoc=true;
            continue;
        }

        if(string(argv[i]) == "--scale" ){
            scale=destringify<double>(argv[i+1]);
            i++;
            specifiedScale=true;
            continue;
        }

        if(string(argv[i]) == "-m" ){
            minimumFragSize = destringify<int>(argv[i+1]);
            i++;
            continue;
        }

        if(string(argv[i]) == "-M" ){
            maximumFragSize = destringify<int>(argv[i+1]);
            i++;
            continue;
        }
	
	if(string(argv[i]) == "-s" ){
            fileSizeFrag=      string(argv[i+1]);
            i++;
	    fileSizeFragB=true;
            continue;
        }

        if(string(argv[i]) == "-l" ){
            sizeFragments=destringify<int>(argv[i+1]);
            i++;
	    specifiedLength = true;
            continue;
        }

	if(string(argv[i]) == "-n" ){
            nFragments=destringify<int>(argv[i+1]);
            i++;
            continue;
        }

        cerr<<"Error: unknown option "<<string(argv[i])<<endl;
        return 1;

    }

    //cerr<<"done"<<endl;


    if(specifiedLoc && !specifiedScale){
        cerr<<"Error: cannot specify --loc without --scale"<<endl;
        return 1;
    }

    if(specifiedScale && !specifiedLoc){
        cerr<<"Error: cannot specify --scale without --loc"<<endl;
        return 1;
    }

    if(specifiedScale && specifiedLoc){
	if(specifiedLength){ 
	    cerr<<"Error: cannot specify --scale and --loc with -l"<<endl;
	    return 1;
	}

	if(fileSizeFragB){ 
	    cerr<<"Error: cannot specify --scale and --loc with -s"<<endl;
	    return 1;
	}

	distribution = lognormal_distribution<double>(location,scale);


    }

    // for(int i=0;i<10000;i++){
    // 	double number = distribution(generator);
    // 	cout<<number<<endl;
    // }
    // return 1;

    if(specifiedLength && fileSizeFragB){ 
	cerr<<"Error: cannot specify -l and -s"<<endl;
	return 1;
    }

    vector<int> sizeFragList;
    if(fileSizeFragB){ 

	igzstream sizeFragfd;

	sizeFragfd.open(fileSizeFrag.c_str(), ios::in);

	if (sizeFragfd.good()){

	    while ( getline (sizeFragfd,line)){
		sizeFragList.push_back( destringify<int>( line ) );
	    }
	    sizeFragfd.close();

	}else{
	    cerr << "Unable to open size distribution file "<<fileSizeFrag<<endl;
	    return 1;
	}
    }
    
    // for(unsigned int i=0;i<sizeFragList.size();i++){
    // 	cout<<sizeFragList[i]<<endl;
    // }
    // return 1;

	




    int sub5pPlusi  = -1000;
    int sub5pMinusi = -1000;
    int sub3pPlusi  = -1000;
    int sub3pMinusi = -1000;

    if(compFileSpecified){

	igzstream compFilefd;

	compFilefd.open(compFile.c_str(), ios::in);

	if (compFilefd.good()){

	    while ( getline (compFilefd,line)){
		//cerr<<line<<endl;
		if(strBeginsWith(line,"#") ||
		   strBeginsWith(line,"Chr") )
		    continue;
		
		vector<string> fields = allTokens( line , '\t');
		if(fields.size()>7){
		    int d     = destringify<int>(fields[3]);
		    double aC = destringify<int>(fields[4]);
		    double cC = destringify<int>(fields[5]);
		    double gC = destringify<int>(fields[6]);
		    double tC = destringify<int>(fields[7]);
		    double bC = destringify<int>(fields[8]);

		    // cerr<<"d "<<d<<endl;

		    //cout<<line<<endl;
		    
		    if(abs(d)<=distFromEnd){
#ifdef DEBUG
			cerr<<line<<endl;
#endif
			if(fields[1] == "3p" && fields[2] == "+" ){
			    if(sub3pPlusi  == -1000){ sub3pPlusi  = d; }else{
				if(d<sub3pPlusi){
				    cerr<<"Parse error, wrong order of distance to fragment end  for line "<<line<<endl;
				    return 1;
				}				    
			    }

			    subrates toadd;
			    toadd.s[0] = aC/bC;
			    toadd.s[1] = cC/bC;
			    toadd.s[2] = gC/bC;
			    toadd.s[3] = tC/bC;
			    sub3pPlus.push_back(toadd);
			}

			if(fields[1] == "3p" && fields[2] == "-" ){
			    
			    if(sub3pMinusi  == -1000){ sub3pMinusi  = d; }else{
				if(d<sub3pMinusi){
				    cerr<<"Parse error, wrong order of distance to fragment end  for line "<<line<<endl;
				    return 1;
				}				    
			    }
			    
			    subrates toadd;
			    toadd.s[0] = aC/bC;
			    toadd.s[1] = cC/bC;
			    toadd.s[2] = gC/bC;
			    toadd.s[3] = tC/bC;
			    sub3pMinus.push_back(toadd);
			}

			if(fields[1] == "5p" && fields[2] == "+" ){
			    if(sub5pPlusi  == -1000){ sub5pPlusi  = d; }else{
				if(d<sub5pPlusi){
				    cerr<<"Parse error, wrong order of distance to fragment end  for line "<<line<<endl;
				    return 1;
				}				    
			    }

			    subrates toadd;
			    toadd.s[0] = aC/bC;
			    toadd.s[1] = cC/bC;
			    toadd.s[2] = gC/bC;
			    toadd.s[3] = tC/bC;
			    sub5pPlus.push_back(toadd);
			}

			if(fields[1] == "5p" && fields[2] == "-" ){
			    if(sub5pMinusi  == -1000){ sub5pMinusi  = d; }else{
				if(d<sub5pMinusi){
				    cerr<<"Parse error, wrong order of distance to fragment end  for line "<<line<<endl;
				    return 1;
				}				    
			    }
			    
			    subrates toadd;
			    toadd.s[0] = aC/bC;
			    toadd.s[1] = cC/bC;
			    toadd.s[2] = gC/bC;
			    toadd.s[3] = tC/bC;
			    sub5pMinus.push_back(toadd);
			}


		    }
		}
		//chrall+=line;

	    }
	    compFilefd.close();
	}else{
	    cerr << "Unable to open composition file "<<compFile<<endl;
	    return 1;
	}

    }else{//else no comp file
	distFromEnd=0;
    }


    // string bamfiletopen = string(argv[argc-1]);//bam file
    string fastaFile    = string(argv[argc-1]);//fasta file


#ifdef DEBUG    
    cerr<<"----------------"<<endl;
    cerr<<"5p +"<<endl;
    for(int i=0;i<2*distFromEnd;i++){
	cerr<<i<<"\t"<<sub5pPlus[i].s[0]<<"\t"<<sub5pPlus[i].s[1]<<"\t"<<sub5pPlus[i].s[2]<<"\t"<<sub5pPlus[i].s[3]<<endl;
    }
    cerr<<"----------------"<<endl;
    cerr<<"5p -"<<endl;
    for(int i=0;i<2*distFromEnd;i++){
	cerr<<i<<"\t"<<sub5pMinus[i].s[0]<<"\t"<<sub5pMinus[i].s[1]<<"\t"<<sub5pMinus[i].s[2]<<"\t"<<sub5pMinus[i].s[3]<<endl;
    }
    cerr<<"----------------"<<endl;
   cerr<<"3p +"<<endl;
    for(int i=0;i<2*distFromEnd;i++){
	cerr<<i<<"\t"<<sub3pPlus[i].s[0]<<"\t"<<sub3pPlus[i].s[1]<<"\t"<<sub3pPlus[i].s[2]<<"\t"<<sub3pPlus[i].s[3]<<endl;
    }
    cerr<<"----------------"<<endl;
    cerr<<"3p -"<<endl;
    for(int i=0;i<2*distFromEnd;i++){
	cerr<<i<<"\t"<<sub3pMinus[i].s[0]<<"\t"<<sub3pMinus[i].s[1]<<"\t"<<sub3pMinus[i].s[2]<<"\t"<<sub3pMinus[i].s[3]<<endl;
    }
    cerr<<"----------------"<<endl;
#endif

    //return 1;
    //TODO mmap
    igzstream myFile;
    string defline;
    myFile.open(fastaFile.c_str(), ios::in);
    string chrall="";
    getline (myFile,line);//header
    defline = line;
    if (myFile.good()){
	while ( getline (myFile,line)){
	    //cout<<line<<endl;
	    //chrall+=line;
	    chrall =line;
	}
	myFile.close();
    }else{
	cerr << "Unable to open file "<<fastaFile<<endl;
	return 1;
    }


    unsigned int   f       = 0;

    while(f<nFragments){

	int length      = 0;

	if(fileSizeFragB  ){ length=sizeFragList[ randomInt(0,int(sizeFragList.size())) ]; 	}
	if(specifiedLength){ length=sizeFragments;                                      	}
	if(specifiedScale ){ length=int(distribution(generator));                              	}
	


	int idx         = randomInt(distFromEnd,int(chrall.size())-length-distFromEnd);
	string temp     = chrall.substr(idx-distFromEnd,length+2*distFromEnd);
	bool plusStrand = randomBool();
	string preFrag  = "";
	string posFrag  = "";
	string frag     = "";

	if(!isResolvedDNAstring(temp))
	    continue;
	// cout<<length<<endl;
	// cout<<distFromEnd<<endl;
	//cerr<<"t:"<<temp<<endl;

	string deflineToPrint = defline+":";
	if(!plusStrand){
	    //cout<<"-"<<endl;	    

	    temp    = reverseComplement(temp);
	    preFrag = temp.substr(                 0 , distFromEnd);
	    frag    = temp.substr(       distFromEnd ,      length);
	    posFrag = temp.substr(length+distFromEnd , distFromEnd);

	    deflineToPrint = deflineToPrint+"-:"+stringify(idx)+":"+stringify(idx+length)+":"+stringify(length);
	}else{

	    //cout<<"+"<<endl;
	    preFrag = temp.substr(                 0 , distFromEnd);
	    frag    = temp.substr(       distFromEnd ,      length);
	    posFrag = temp.substr(length+distFromEnd , distFromEnd);

	    deflineToPrint = deflineToPrint+"+:"+stringify(idx)+":"+stringify(idx+length)+":"+stringify(length);
	}

	// cerr<<"---------"<<endl;
	// cerr<<temp<<endl;
	// cerr<<preFrag<<endl;
	// cerr<<frag<<endl;
	// cerr<<posFrag<<endl;


	if(!compFileSpecified){ //no composition file, just produce the sequences	    
	    cout<<frag<<endl;
	    f++;
	}else{
	    double probAcc = 1.0;
	    if(plusStrand){

		//5p before frag
		for(int i=0;i<distFromEnd;i++){
		    probAcc *= sub5pPlus[            i].s[baseResolved2int( preFrag[i]                   )];
		}

		//5p in frag
		for(int i=0;i<distFromEnd;i++){
		    probAcc *= sub5pPlus[i+distFromEnd].s[baseResolved2int(    frag[i]                   )];
		}

		//3p in frag
		for(int i=0;i<distFromEnd;i++){
		    probAcc *= sub3pPlus[            i].s[baseResolved2int(    frag[length-distFromEnd+i])];
		}

		//3p after frag
		for(int i=0;i<distFromEnd;i++){
		    probAcc *= sub3pPlus[i+distFromEnd].s[baseResolved2int( posFrag[i]                   )];
		}

		if(randomProb()<probAcc){
		    cout<<deflineToPrint<<"\n"<<frag<<endl;
		    f++;
		 }else{
		     //discard
		}
	    }else{

		//5p before frag
		for(int i=0;i<distFromEnd;i++){
		    probAcc *= sub5pMinus[            i].s[baseResolved2int( preFrag[i]                   )];
		}

		//5p in frag
		for(int i=0;i<distFromEnd;i++){
		    probAcc *= sub5pMinus[i+distFromEnd].s[baseResolved2int(    frag[i]                   )];
		}

		//3p in frag
		for(int i=0;i<distFromEnd;i++){
		    probAcc *= sub3pMinus[            i].s[baseResolved2int(    frag[length-distFromEnd+i])];
		}

		//3p after frag
		for(int i=0;i<distFromEnd;i++){
		    probAcc *= sub3pMinus[i+distFromEnd].s[baseResolved2int( posFrag[i]                   )];
		}

		if(randomProb()<probAcc){
		    cout<<deflineToPrint<<"\n"<<frag<<endl;
		    f++;
		 }else{
		     //discard
		}


	    }
	    
	}

    }



    return 0;
}


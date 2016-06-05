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
#include <cstring>

#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
//#define DEBUG

#include "utils.h"

using namespace std;

typedef struct{
    double s[4];
} subrates;


typedef struct{
    string name;
    uint64_t startIndexChr;
    uint64_t endIndexChr;
    uint64_t length;
} chrinfo;


/* it is the structure created by samtools faidx */
typedef struct faidx1_t{
    int32_t line_len, line_blen;
    int64_t len;
    uint64_t offset;
}faidx1_t,*FaidxPtr;


uint64_t  randGenomicCoord(const uint64_t & genomeLength){     

     while(true){
	 uint64_t toReturn = 
	     (((uint64_t) rand() <<  0) & 0x000000000000FFFFull) | 
	     (((uint64_t) rand() << 16) & 0x00000000FFFF0000ull) | 
	     (((uint64_t) rand() << 32) & 0x0000FFFF00000000ull) |
	     (((uint64_t) rand() << 48) & 0xFFFF000000000000ull);
	 //cout<<toReturn<<endl;
	 //should work if genome size is less than 2^64, 
	 //creates a slight but negligeable bias towards the begininning of the genome
	 return toReturn%genomeLength;
     }
 }


/**
 * wrapper for a mmap, a fileid and some faidx indexes
 */
class IndexedGenome{
private:
    /* used to get the size of the file */
    struct stat buf;
    /* genome fasta file file descriptor */
    int fd;
    /* the mmap (memory mapped) pointer */
    char *mapptr;
    /** reads an fill a string */
    bool readline(gzFile in,string& line){
	if(gzeof(in)) return false;
	line.clear();
	int c=-1;
	while((c=gzgetc(in))!=EOF && c!='\n') 
	    line+=(char)c;
	return true;
    }
    
public:
    /* maps a chromosome to the samtools faidx index */
    map<string,faidx1_t> name2index;

    /** constructor 
     * @param fasta: the path to the genomic fasta file indexed with samtools faidx
     */
    IndexedGenome(const char* fasta):fd(-1),mapptr(NULL){
	string faidx(fasta);
	//cout<<fasta<<endl;
	string line;
	faidx+=".fai";
	/* open *.fai file */
	//cout<<faidx<<endl;
	ifstream in(faidx.c_str(),ios::in);
	if(!in.is_open()){
	    cerr << "Cannot open faidx: " << faidx << endl;
	    exit(EXIT_FAILURE);
	}

	/* read indexes in fai file */
	while(getline(in,line,'\n')){
	    if(line.empty()) 
		continue;
	    const char* p=line.c_str();
	    char* tab=(char*)strchr(p,'\t');
	    if(tab==NULL) 
		continue;
	    string chrom(p,tab-p);
	    ++tab;
	    faidx1_t index;
	    if(sscanf(tab,"%ld\t%ld\t%d\t%d",&index.len, &index.offset, &index.line_blen,&index.line_len)!=4){
		cerr << "Cannot read index in "<< line << endl;
		exit(EXIT_FAILURE);
	    }
	    /* insert in the map(chrom,faidx) */
	    name2index.insert(make_pair(chrom,index));
	}
	/* close index file */
	in.close();
	
	/* get the whole size of the fasta file */
	if(stat(fasta, &buf)!=0){
	    perror("Cannot stat");
	    exit(EXIT_FAILURE);
	}
	
	/* open the fasta file */
	fd = open(fasta, O_RDONLY);
	if (fd == -1){
	    perror("Error opening file for reading");
	    exit(EXIT_FAILURE);
	}
	/* open a memory mapped file associated to this fasta file descriptor */
	mapptr = (char*)mmap(0, buf.st_size, PROT_READ, MAP_SHARED, fd, 0);
	if (mapptr == MAP_FAILED){
	    close(fd);
	    perror("Error mmapping the file");
	    exit(EXIT_FAILURE);
	}
    } //end Constructor

    /* destructor */
    ~IndexedGenome(){
	/* close memory mapped map */
	if(mapptr!=NULL && munmap(mapptr,buf.st_size) == -1){
	    perror("Error un-mmapping the file");
	}

	/* dispose fasta file descriptor */
	if(fd!=-1) 
	    close(fd);
    }

    
    string fetchSeq(const FaidxPtr faidx,int64_t index, int lengthFragment){
	//cout<<"fetchSeq"<<endl;
	int64_t index2=index;
	string strToPrint="";
	    
	for(int j=0;j<lengthFragment;j++){ //for each char
	    long pos= faidx->offset +
		index2 / faidx->line_blen * faidx->line_len +
		index2 % faidx->line_blen
		;
	    //cerr<<char(toupper(mapptr[pos]));
	    strToPrint+=char(toupper(mapptr[pos]));
	    index2++;
	}
	       
	return strToPrint;   	
    }

};

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
    bool noRev=false;

    const string usage=string("\nThis program takes a fasta file representing a chromosome and generates\n")+
	"aDNA fragments according to a certain distribution\n\n"+
	string(argv[0])+" [options]  [chromosome fasta] "+"\n\n"+
	
	"\t\t"+"-n\t"+"[number]" +"\t\t\t"+"Generate [number] fragments (default: "+stringify(nFragments)+")"+"\n"+
	"\n"+
	"\t\t"+"--comp\t"+"[file]"+"\t\t\t\t"+"Base composition for the fragments (default none)"+"\n"+
	"\t\t"+"--dist\t"+"[file]"+"\t\t\t\t"+"Distance from ends to consider  (default: "+stringify(distFromEnd)+")"+"\n"+
	"\t\t"+"\t"+""+"\t\t\t\t"+"if this is not specified, the base composition"+"\n"+
	"\t\t"+"\t"+""+"\t\t\t\t"+"will only reflect the chromosome file used"+"\n"+
	"\t\t"+"--norev\t"+""+"\t\t\t\t"+"Do not reverse complement (default: rev. comp half of seqs.)"+"\n"+

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

	if(string(argv[i]) == "--norev" ){
	    noRev=true;
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
    // if(compFileSpecified){
    // 	distFromEnd=0;
    // }
    //cerr<<distFromEnd<<endl;

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



    if(specifiedLength && fileSizeFragB){ 
	cerr<<"Error: cannot specify -l and -s"<<endl;
	return 1;
    }

    if(!specifiedScale && !specifiedLength  && !fileSizeFragB){
	cerr<<"Error: must specify -l, -s or either the log-normal parameters for the fragment size."<<endl;
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

    timeval time;
    gettimeofday(&time, NULL);
    srand(  long((time.tv_sec * 1000) + (time.tv_usec / 1000)) );





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
    IndexedGenome* genome=new IndexedGenome(fastaFile.c_str());
    cerr<<"Mapped "<<fastaFile<<" into memory"<<endl;

    uint64_t genomeLength=0;
    vector<chrinfo> chrFound;

    typedef map<string,faidx1_t>::iterator it_type;
    for(it_type iterator = genome->name2index.begin(); iterator != genome->name2index.end(); iterator++) {
	//cout<<iterator->first<<"\t"<<iterator->second.len<<"\t"<<genomeLength<<endl;
	
	chrinfo toadd;

	toadd.name          = iterator->first;
	toadd.startIndexChr = genomeLength+1;
	toadd.length        = iterator->second.len;
	toadd.endIndexChr   = genomeLength+iterator->second.len;
	genomeLength       += iterator->second.len;

	chrFound.push_back(toadd);
    }

    


    // //return 1;
    // igzstream myFile;
    // string defline;
    // myFile.open(fastaFile.c_str(), ios::in);
    // string chrall="";
    // getline (myFile,line);//header
    // defline = line;
    // if (myFile.good()){
    // 	while ( getline (myFile,line)){
    // 	    //cout<<line<<endl;
    // 	    //chrall+=line;
    // 	    chrall =line;
    // 	}
    // 	myFile.close();
    // }else{
    // 	cerr << "Unable to open file "<<fastaFile<<endl;
    // 	return 1;
    // }


    unsigned int   f       = 0;

    while(f<nFragments){

	int length      = 0;

	if(fileSizeFragB  ){ length=sizeFragList[ randomInt(0,int(sizeFragList.size())) ]; 	}
	if(specifiedLength){ length=sizeFragments;                                      	}
	if(specifiedScale ){ length=int(distribution(generator));                              	}
	

	//TODO generate chr and coord
	uint64_t idx;//         = randomInt(distFromEnd,int(chrall.size())-length-distFromEnd);
	uint64_t coord;
	bool found=false;
	//string temp     = chrall.substr(idx-distFromEnd,length+2*distFromEnd);
	string temp="";
	string deflineToPrint;// = defline+":";
	while(!found){
	    coord =randGenomicCoord(genomeLength);

	    for(unsigned int i=0;i<chrFound.size();i++){	    
		// cout<<"------"<<i<<"------"<<endl;
		// cout<<chrFound[i].name<<endl;
		// cout<<chrFound[i].startIndexChr<<endl;
		// cout<<chrFound[i].endIndexChr<<endl;
		
		if( (chrFound[i].startIndexChr+distFromEnd) <= coord && coord <= (chrFound[i].endIndexChr-length-distFromEnd+1)){
		    found=true;
		    idx = coord-chrFound[i].startIndexChr;
		    faidx1_t faidxForName=genome->name2index[ chrFound[i].name ];
		    temp=genome->fetchSeq(&faidxForName,coord-chrFound[i].startIndexChr-distFromEnd, length+2*distFromEnd);
		    //cout<<">coord "<<coord<<"\tname\t"<<chrFound[i].name<<"\tid=\t"<<idx<<"\tidx+l\t"<<idx+length<<endl;
		    deflineToPrint = ">"+chrFound[i].name+":";
		    //cout<<temp<<endl;
		    // toReturn.setChrName(       chrFound[i].name);
		    // toReturn.setStartCoord( coord-chrFound[i].startIndexChr);
		    // toReturn.setEndCoord(   coord-chrFound[i].startIndexChr+bpToExtract);
		    break;
		}		
	    }
	    if(found) break;
	    //cout<<"Not found coord="<<coord<<endl;
	}
	//return 1;

	bool plusStrand;
	if(noRev){
	    plusStrand = true;
	}else{
	    plusStrand = randomBool();
	}
	string preFrag  = "";
	string posFrag  = "";
	string frag     = "";

	if(!isResolvedDNAstring(temp))
	    continue;
	// cout<<length<<endl;
	// cout<<distFromEnd<<endl;
	//cerr<<"t:"<<temp<<endl;
	

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
	    cout<<deflineToPrint<<"\n"<<frag<<endl;
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


/*
 * fragSim
 * Date: Jun-02-2016 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>
#include <gzstream.h>
#include <inttypes.h>
#include <stdio.h>
#include <numeric>

#include <random>
#include <cstring>

#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <api/SamHeader.h>
#include <api/BamMultiReader.h>
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <api/BamAux.h>
#include "PutProgramInHeader.h"

#include "libgab.h"
#include "FastQParser.h"

//#define DEBUG


using namespace std;
using namespace BamTools;

const uint32_t flagSingleReads =  4; // 00000100


typedef struct{
    double s[4];
} subrates;

typedef struct{
    int    l;    //fragment length
    double cProb;//cumulative probability 
} freqFragLength;


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

inline double maxSubrates(const subrates & sub){     
    double maxsub = sub.s[0];
    if(sub.s[1]>maxsub) 
	maxsub = sub.s[1];
    if(sub.s[2]>maxsub) 
	maxsub = sub.s[2];
    if(sub.s[3]>maxsub) 
	maxsub = sub.s[3];
    return maxsub;
}

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


template <typename T>
T pdfNorm(T x, T m, T s){
    static const T inv_sqrt_2pi = 0.3989422804014327;
    T a = (x - m) / s;

    return inv_sqrt_2pi / s * exp(-T(0.5) * a * a);
}

double meanNorm(const vector<double> & v){
    double sumVector  = accumulate(v.begin(), v.end(), 0.0);
    return (sumVector / double(v.size()));
}


double stdevNorm(const vector<double> & v,const double meanV){
    vector<double> dTemp(v.size());
    transform(v.begin(), v.end(), dTemp.begin(),bind2nd(minus<double>(), meanV));
    double sumSq = inner_product(dTemp.begin(), dTemp.end(), dTemp.begin(), 0.0);
    return (sqrt(sumSq / v.size()));

}

//called only if uniqTags is used
void checkUniqueness(string * deflineToPrint,
		     map<string, unsigned int> * fragmentID2Count,
		     const bool tagb ,
		     const string & tag){

    if (  fragmentID2Count->find(*deflineToPrint) == fragmentID2Count->end() ) {
	fragmentID2Count->insert( pair<string,unsigned int>( *deflineToPrint, 1) );
	*deflineToPrint+="_1";
    }else{
	*deflineToPrint+="_"+stringify(++fragmentID2Count->at(*deflineToPrint));
    }
    
    if(tagb ){ //the tagb and not uniqTags was taken care of earlier
	*deflineToPrint=*deflineToPrint+"_"+tag;	    
    }
    
}

inline int selectFragmentLength(const bool fileSizeFragB,
				const bool fileSizeFragBFreq,
				const bool specifiedLength,
				const bool specifiedScale,
				const vector<int> & sizeFragList,
				const vector<freqFragLength> & fragLengthFreq,
				const int & sizeFragments,
				default_random_engine & generator,
				lognormal_distribution<double> & distribution){
    int length=0;
    //From file
    if(fileSizeFragB          ){ 
	int randIndex = randomInt(0,int(sizeFragList.size())-1);
	length=sizeFragList[ randIndex ]; 	

    }
    else
	//from frequency
	if(fileSizeFragBFreq      ){ 
		
	    double p          =  randomProb();
	    double cProbLower =  0.0;
	    bool foundL        =  false;
	    unsigned int indexFoundL=-1;
	    for(unsigned int i=0;i<fragLengthFreq.size();i++){
		if(cProbLower <=  p 
		   &&
		   p          <=  fragLengthFreq[i].cProb){						
		    indexFoundL = i;
		    foundL       = true;
		    cProbLower  = fragLengthFreq[i].cProb;
		    break;
		}
	    }

	    if(foundL){
		length=fragLengthFreq[ indexFoundL ].l; 	
	    }else{
		length=fragLengthFreq[ randomInt(0, int(fragLengthFreq.size())-1 ) ].l; 	
	    }
	}	
	else
	    //fixed length
	    if(specifiedLength    ){ length=sizeFragments;                                      	}
	    else
		//scale and loc
		if(specifiedScale ){ length=int(distribution(generator));                              	}

    return length;
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

    bool uppercase;

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
    IndexedGenome(const char* fasta,const char* fastai,bool _uppercase):fd(-1),mapptr(NULL),uppercase(_uppercase){
	//string faidx(fasta);
	//faidx+=".fai";
	
	string faidx(fastai);
	
	//cout<<fasta<<endl;
	string line;

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

	    vector<string> tokens=allTokens(line,'\t');
	    if(tokens.size()!=5){
		cerr<<"Parse error in "<<line<<endl;
		exit(EXIT_FAILURE);
	    }

	    // tokens[0] name
	    index.len       = destringify<int64_t>( tokens[1]);
	    index.offset    = destringify<uint64_t>(tokens[2]);
	    index.line_blen = destringify<int32_t>( tokens[3]);
	    index.line_len  = destringify<int32_t>( tokens[4]);

	    // if(sscanf(tab,"%"PRIu64"\t%"PRId64"\t%d\t%d",&index.len, &index.offset, &index.line_blen,&index.line_len)!=4){
	    // 	cerr << "Cannot read index in "<< line << endl;
	    // 	exit(EXIT_FAILURE);
	    // }
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
	//cerr<<"fetchSeq "<< uppercase<<endl;
	int64_t index2=index;
	string strToPrint="";
	    
	if(uppercase){
	    for(int j=0;j<lengthFragment;j++){ //for each char	    
		long pos= faidx->offset +
		    index2 / faidx->line_blen * faidx->line_len +
		    index2 % faidx->line_blen;
		strToPrint+=char(toupper(mapptr[pos]));
		index2++;
	    }
	}else{
	    for(int j=0;j<lengthFragment;j++){ //for each char	    
		long pos= faidx->offset +
		    index2 / faidx->line_blen * faidx->line_len +
		    index2 % faidx->line_blen;
		strToPrint+=char(        mapptr[pos] );
		index2++;
	    }
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
    double gcBias             = 0;
    bool   gcBiasB            = false;


    
    vector<subrates> sub5pPlus;
    vector<subrates> sub5pMinus;
    vector<subrates> sub3pPlus;
    vector<subrates> sub3pMinus;

    string fileSizeFrag;
    bool   fileSizeFragB=false;
    string fileSizeFragFreq;
    bool   fileSizeFragBFreq=false;

    string tmpDir = "/tmp/";
    
    int minimumFragSize=   0;
    int maximumFragSize=1000;
    bool noRev=false;

    string outFastagz       ;
    bool   outFastagzb=false;
    string outBAM           ;
    bool   outBAMb    =false;
    bool   tagb       =false;
    bool   uppercase  =true;
    bool   circ       =false;
    string circName   ="";
    uint64_t circOffset=0;

    
    string   tag              = "";
    bool uniqTags=false;
    bool fastqMode=false;
    
    bool produceUnCompressedBAM=false;
    map<string, unsigned int> fragmentID2Count;
    
    const string usage=string("\n This program takes a fasta file representing a chromosome and generates\n")+
	" aDNA fragments according to a certain distribution\n\n"+
	" "+string(argv[0])+" [options]  [chromosome fasta] "+"\n\n"+

	
	"\t\t"+"-n\t"+"[number]" +"\t\t\t"+"Generate [number] fragments (default: "+stringify(nFragments)+")"+"\n"+
	"\n"+
	"\t\t"+"--comp\t"+"[file]"+"\t\t\t\t"+"Base composition for the fragments (default none)"+"\n"+
	"\t\t"+"--dist\t"+"[file]"+"\t\t\t\t"+"Distance from ends to consider for base composition (default: "+stringify(distFromEnd)+")"+"\n"+
	"\t\t"+"\t"+""+"\t\t\t\t"+"if this is not specified, the base composition"+"\n"+
	"\t\t"+"\t"+""+"\t\t\t\t"+"will only reflect the chromosome file used"+"\n"+
	"\t\t"+"--norev\t"+""+"\t\t\t\t"+"Do not reverse complement (default: rev. comp half of seqs.)"+"\n"+
	"\t\t"+"--case\t"+""+"\t\t\t\t"+"Do not set the sequence to upper-case (default: uppercase the seqs.)"+"\n"+
	"\t\t"+"--fq\t"+""+"\t\t\t\t"+"The chromosome file is a fastq file. Suitable if you have pre-selected the fragments and\n"+
	"\t\t"+"    \t"+""+"\t\t\t\t"+"just want to trim them according to certain fragment lengths\n"+
	"\t\t"+"    \t"+""+"\t\t\t\t"+"please use either fastq.gz or fq.gz (Default:  "+booleanAsString(fastqMode)+")"+"\n"+
	"\t\t"+"    \t"+""+"\t\t\t\t"+"Please note that this mode will write fastq to STDOUT\n"+
	"\t\t"+"--circ\t"+"[REF NAME]"+"\t\t\t\t"+"Assume [REF NAME] is circular"+"\n"+

	"\n"+
	"\tOutput options\n"+
	"\t\t"+"-b\t"+"[bam out]"   +"\t\t\t"+"Write output as a BAM file (default: fasta in STDOUT)"+"\n"+
	"\t\t"+"-o\t"+"[fasta out]" +"\t\t\t"+"Write output as a zipped fasta (default: fasta in STDOUT)"+"\n"+
	"\t\t"+"-u"                 +"\t\t\t\t\t"+"Produce uncompressed BAM (good for UNIX pipe)"+"\n"+
	"\t\t"+"-tag\t" +"[tag]\t\t\t\t"+"Append this string to deflines or BAM tags (Default:  "+booleanAsString(tagb)+")"+"\n"+
	"\t\t"+"-tmp\t" +"[tmp dir]\t\t\t"+"Use this directory as the temporary dir for zipped files (default:  "+tmpDir+")"+"\n"+
	"\t\t"+"-uniq\t" +"\t\t\t\t"+"Make sure that the fragment names are unique by appending a suffix (default:  "+booleanAsString(uniqTags)+")"+"\n"+
	"\t\t"+"\t" +"\t\t\t\t"+"note: this might not be practical for large datasets"+"\n"+

	"\n"+
	"\tFragment size: \n"+
	"\t\t"+"-m\t"+"[length]" +"\t\t\t"+"Minimum fragments length < (default: "+stringify(minimumFragSize)+")"+"\n"+
	"\t\t"+"-M\t"+"[length]" +"\t\t\t"+"Maximum fragments length > (default: "+stringify(maximumFragSize)+")"+"\n"+
	"\n"+
	"\tFragment size distribution: specify either one of the 4 possible options\n"+
	"\t\t"+"-l\t"+"[length]" +"\t\t\t"+"Generate fragments of fixed length  (default: "+stringify(sizeFragments)+")"+"\n"+
	"\t\t"+"-s\t"+"[file]"   +"\t\t\t\t"+"Open file with size distribution (one fragment length per line)"+"\n"+
	"\t\t"+"-f\t"+"[file]"   +"\t\t\t\t"+"Open file with size frequency in the following format:"+"\n"+
	"\t\t"+"\t"   +""        +"\t\t\t\t\t"+"length[TAB]freq\tex:"+"\n"+
	"\t\t"+"\t"   +""        +"\t\t\t\t\t"+"40\t0.0525"+"\n"+
	"\t\t"+"\t"   +""        +"\t\t\t\t\t"+"41\t0.0491"+"\n"+
	"\t\t"+"\t"   +""        +"\t\t\t\t\t"+"..."+"\n"+


	"\n\t\tLength options:\n"+                                
	"\t\t\t"+"--loc\t"+"[file]"   +  "\t\t\t"+"Location for lognormal distribution (default none)"+"\n"+
	"\t\t\t"+"--scale\t"+"[file]"   +"\t\t\t"+"Scale for lognormal distribution      (default none)"+"\n"+
	"\n\t\tGC bias options:\n"+
	"\t\t"+"-gc\t"+"[gc bias]" +"\t\t\t"+"Use GC bias factor  (default: "+stringify(gcBias)+")"+"\n"+	
	
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

	if(string(argv[i]) == "-uniq" ){
	    uniqTags = true;
	    i++; 
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

	if(string(argv[i]) == "-o" ){
	    outFastagz  = string(argv[i+1]);
	    i++;
	    outFastagzb = true;
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

	if(string(argv[i]) == "--fq"  ){
	    fastqMode=true;
	    continue;
	}

	if(string(argv[i]) == "--case"  ){
	    uppercase=false;
	    continue;
	}

	if(string(argv[i]) == "--circ"  ){
	    circ=true;
	    circName = string(argv[i+1]);
	    i++;
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

	if(string(argv[i]) == "-f" ){
	    fileSizeFragFreq=      string(argv[i+1]);
	    i++;
	    fileSizeFragBFreq=true;
	    continue;
	}



	if(string(argv[i]) == "-gc" ){
	    gcBias  = destringify<double>(argv[i+1]);
	    i++;
	    gcBiasB = true;
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

	if(string(argv[i]) == "-tag" ){
	    tagb=true;
	    tag       = string(argv[i+1]);
	    i++;
	    continue;
	}

	if(string(argv[i]) == "-tmp" ){
	    tmpDir    = string(argv[i+1]);
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

    if(specifiedLength){
	if(sizeFragments>maximumFragSize){
	    cerr<<"Error: specified size "<<sizeFragments<<" is longer than the maximum allowed fragment size "<<maximumFragSize<<", use option -M to change this."<<endl;
	    return 1;
	}
    }

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

	if(fileSizeFragBFreq){ 
	    cerr<<"Error: cannot specify --scale and --loc with -f"<<endl;
	    return 1;
	}

	distribution = lognormal_distribution<double>(location,scale);

    }

    if(outFastagzb && outBAMb){ 
	cerr<<"Error: cannot specify both -o and -b"<<endl;
	return 1;
    }



    if(fileSizeFragB && fileSizeFragBFreq){ 
	cerr<<"Error: cannot specify -s with -f"<<endl;
	return 1;
    }


    if(specifiedLength && fileSizeFragB){ 
	cerr<<"Error: cannot specify -l and -s"<<endl;
	return 1;
    }

    if(specifiedLength && fileSizeFragBFreq){ 
	cerr<<"Error: cannot specify -l and -f"<<endl;
	return 1;
    }


    if(!specifiedScale && !specifiedLength  && !fileSizeFragB && !fileSizeFragBFreq){
	//cerr<<"Error: must specify -l, -s, -f  or either the log-normal parameters for the fragment size."<<endl;
	//return 1;
	//artificially setting specifiedLength to true
	specifiedLength = true;
    }

    // Use fragment size list
    vector<int> sizeFragList;
    if(fileSizeFragB){ 

	igzstream sizeFragfd;

	sizeFragfd.open(fileSizeFrag.c_str(), ios::in);

	if (sizeFragfd.good()){

	    while ( getline (sizeFragfd,line)){
		vector<string> fields = allTokens( line , '\t');
		if(fields.size() != 1){
		    cerr << "The following line "<<line<<" in file "<<fileSizeFrag<<" does not have 1 field"<<endl;
		    return 1;		    
		}

		int t = destringify<int>( line );
		if(t>=1)
		    sizeFragList.push_back( t );
	    }
	    sizeFragfd.close();

	}else{
	    cerr << "Unable to open size distribution file "<<fileSizeFrag<<endl;
	    return 1;
	}
    }


    // Use fragment size frequencies
    vector<freqFragLength> fragLengthFreq;

    if(fileSizeFragBFreq){ 

	igzstream sizeFragfd;

	sizeFragfd.open(fileSizeFragFreq.c_str(), ios::in);
	double cumulProb=0.0;
	if (sizeFragfd.good()){
	    while ( getline (sizeFragfd,line)){
		vector<string> fields = allTokens( line , '\t');
		if(fields.size() != 2){
		    cerr << "The following line "<<line<<" in file "<<fileSizeFragFreq<<" does not have 2 fields"<<endl;
		    return 1;		    
		}

		cumulProb    += destringify<double>( fields[1] ); 
		freqFragLength toadd;
		toadd.l       = destringify<int>(    fields[0] ); 
		toadd.cProb   = cumulProb;
		fragLengthFreq.push_back(toadd);
	    }
	    sizeFragfd.close();

	}else{
	    cerr << "Unable to open size frequency file "<<fileSizeFragFreq<<endl;
	    return 1;
	}

	if( (cumulProb<0.99) ||
	    (cumulProb>1.01) ){
	    cerr<<"Problem in file: "<<fileSizeFragFreq<<", the sum of frequencies does not sum to 1, sum="<<cumulProb<<endl;
	    return 1;
	}
	    
    }

    

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
	unsigned int aCtotal=0;
	unsigned int cCtotal=0;
	unsigned int gCtotal=0;
	unsigned int tCtotal=0;
	unsigned int bCtotal=0;
	subrates background;

	if (compFilefd.good()){
#ifdef DEBUG    
	    cerr<<"composition"<<endl<<"end\tstrand\tpos\tf(A)\tf(C)\tf(G)\tf(T)"<<endl;
#endif

	    while ( getline (compFilefd,line)){
		//cerr<<line<<endl;

		if(strBeginsWith(line,"#") ||
		   strBeginsWith(line,"Chr") )
		    continue;
		//unsigned int
		

		vector<string> fields = allTokens( line , '\t');
		if(fields.size()>7){
		    int d  = destringify<int>(fields[3]);
		    unsigned int aC = destringify<unsigned int>(fields[4]);
		    unsigned int cC = destringify<unsigned int>(fields[5]);
		    unsigned int gC = destringify<unsigned int>(fields[6]);
		    unsigned int tC = destringify<unsigned int>(fields[7]);
		    unsigned int bC = destringify<unsigned int>(fields[8]);

		    // cerr<<"d "<<d<<endl;

		    //cout<<line<<endl;
		    
		    if(abs(d)<=distFromEnd){

#ifdef DEBUG
			//cerr<<line<<endl;
#endif

			if(fields[1] == "3p" && fields[2] == "+" ){
			    if(sub3pPlusi  == -1000){ sub3pPlusi  = d; }else{
				if(d<sub3pPlusi){
				    cerr<<"Parse error, wrong order of distance to fragment end  for line "<<line<<endl;
				    return 1;
				}				    
			    }

			    subrates toadd;
			    toadd.s[0] = double(aC)/double(bC);
			    toadd.s[1] = double(cC)/double(bC);
			    toadd.s[2] = double(gC)/double(bC);
			    toadd.s[3] = double(tC)/double(bC);
			    sub3pPlus.push_back(toadd);
#ifdef DEBUG    
			    cerr<<"3p\t+\t"<<d<<"\t"<<toadd.s[0]<<"\t"<<toadd.s[1]<<"\t"<<toadd.s[2]<<"\t"<<toadd.s[3]<<endl;
#endif
			}

			if(fields[1] == "3p" && fields[2] == "-" ){
			    
			    if(sub3pMinusi  == -1000){ sub3pMinusi  = d; }else{
				if(d<sub3pMinusi){
				    cerr<<"Parse error, wrong order of distance to fragment end  for line "<<line<<endl;
				    return 1;
				}				    
			    }
			    
			    subrates toadd;
			    toadd.s[0] = double(aC)/double(bC);
			    toadd.s[1] = double(cC)/double(bC);
			    toadd.s[2] = double(gC)/double(bC);
			    toadd.s[3] = double(tC)/double(bC);
			    sub3pMinus.push_back(toadd);
#ifdef DEBUG    
			    cerr<<"3p\t-\t"<<d<<"\t"<<toadd.s[0]<<"\t"<<toadd.s[1]<<"\t"<<toadd.s[2]<<"\t"<<toadd.s[3]<<endl;
#endif

			}

			if(fields[1] == "5p" && fields[2] == "+" ){
			    if(sub5pPlusi  == -1000){ sub5pPlusi  = d; }else{
				if(d<sub5pPlusi){
				    cerr<<"Parse error, wrong order of distance to fragment end  for line "<<line<<endl;
				    return 1;
				}				    
			    }

			    subrates toadd;
			    toadd.s[0] = double(aC)/double(bC);
			    toadd.s[1] = double(cC)/double(bC);
			    toadd.s[2] = double(gC)/double(bC);
			    toadd.s[3] = double(tC)/double(bC);
			    sub5pPlus.push_back(toadd);
#ifdef DEBUG    
			    cerr<<"5p\t+\t"<<d<<"\t"<<toadd.s[0]<<"\t"<<toadd.s[1]<<"\t"<<toadd.s[2]<<"\t"<<toadd.s[3]<<endl;
#endif
			}

			if(fields[1] == "5p" && fields[2] == "-" ){
			    if(sub5pMinusi  == -1000){ sub5pMinusi  = d; }else{
				if(d<sub5pMinusi){
				    cerr<<"Parse error, wrong order of distance to fragment end  for line "<<line<<endl;
				    return 1;
				}				    
			    }
			    
			    subrates toadd;
			    toadd.s[0] = double(aC)/double(bC);
			    toadd.s[1] = double(cC)/double(bC);
			    toadd.s[2] = double(gC)/double(bC);
			    toadd.s[3] = double(tC)/double(bC);
			    sub5pMinus.push_back(toadd);
#ifdef DEBUG    
			    cerr<<"5p\t-\t"<<d<<"\t"<<toadd.s[0]<<"\t"<<toadd.s[1]<<"\t"<<toadd.s[2]<<"\t"<<toadd.s[3]<<endl;
#endif
			}


		    }else{
			//background
			
			aCtotal+=aC;
			cCtotal+=cC;
			gCtotal+=gC;
			tCtotal+=tC;
			bCtotal+=bC;

		    }
		}//if 7 fields
		//chrall+=line;
		

	    }
	    compFilefd.close();

	    background.s[0] = double(aCtotal)/double(bCtotal);
	    background.s[1] = double(cCtotal)/double(bCtotal);
	    background.s[2] = double(gCtotal)/double(bCtotal);
	    background.s[3] = double(tCtotal)/double(bCtotal);

#ifdef DEBUG    
	    cerr<<"all:\t*\t*\t"<<background.s[0]<<"\t"<<background.s[1]<<"\t"<<background.s[2]<<"\t"<<background.s[3]<<endl;
#endif

	    //correction
	    for(int i=0;i<2*distFromEnd;i++){
		double aT;
		double cT;
		double gT;
		double tT;

		aT = sub5pPlus[i].s[0]-background.s[0];
		cT = sub5pPlus[i].s[1]-background.s[1];
		gT = sub5pPlus[i].s[2]-background.s[2];
		tT = sub5pPlus[i].s[3]-background.s[3];
#ifdef DEBUG    
		cerr<<"5p\t+f\t"<<i<<"\t"<<aT<<"\t"<<cT<<"\t"<<gT<<"\t"<<tT<<"\t"<<(aT+cT+gT+tT)<<endl;
#endif

		sub5pPlus[i].s[0] = 0.25+aT;
		sub5pPlus[i].s[1] = 0.25+cT;
		sub5pPlus[i].s[2] = 0.25+gT;
		sub5pPlus[i].s[3] = 0.25+tT;
#ifdef DEBUG    
		cerr<<"5p\t+a\t"<<i<<"\t"<<sub5pPlus[i].s[0]<<"\t"<<sub5pPlus[i].s[1]<<"\t"<<sub5pPlus[i].s[2]<<"\t"<<sub5pPlus[i].s[3]<<"\t"<<(sub5pPlus[i].s[0]+sub5pPlus[i].s[1]+sub5pPlus[i].s[2]+sub5pPlus[i].s[3])<<endl;
#endif
		aT = sub5pMinus[i].s[0]-background.s[0];
		cT = sub5pMinus[i].s[1]-background.s[1];
		gT = sub5pMinus[i].s[2]-background.s[2];
		tT = sub5pMinus[i].s[3]-background.s[3];
#ifdef DEBUG    
		cerr<<"5p\t-f\t"<<i<<"\t"<<aT<<"\t"<<cT<<"\t"<<gT<<"\t"<<tT<<"\t"<<(aT+cT+gT+tT)<<endl;
#endif
		sub5pMinus[i].s[0] = 0.25+aT;
		sub5pMinus[i].s[1] = 0.25+cT;
		sub5pMinus[i].s[2] = 0.25+gT;
		sub5pMinus[i].s[3] = 0.25+tT;
#ifdef DEBUG    
		cerr<<"5p\t-a\t"<<i<<"\t"<<sub5pMinus[i].s[0]<<"\t"<<sub5pMinus[i].s[1]<<"\t"<<sub5pMinus[i].s[2]<<"\t"<<sub5pMinus[i].s[3]<<"\t"<<(sub5pMinus[i].s[0]+sub5pMinus[i].s[1]+sub5pMinus[i].s[2]+sub5pMinus[i].s[3])<<endl;
#endif
		aT = sub3pPlus[i].s[0]-background.s[0];
		cT = sub3pPlus[i].s[1]-background.s[1];
		gT = sub3pPlus[i].s[2]-background.s[2];
		tT = sub3pPlus[i].s[3]-background.s[3];
#ifdef DEBUG    
		cerr<<"3p\t+f\t"<<i<<"\t"<<aT<<"\t"<<cT<<"\t"<<gT<<"\t"<<tT<<"\t"<<(aT+cT+gT+tT)<<endl;
#endif
		sub3pPlus[i].s[0] = 0.25+aT;
		sub3pPlus[i].s[1] = 0.25+cT;
		sub3pPlus[i].s[2] = 0.25+gT;
		sub3pPlus[i].s[3] = 0.25+tT;
#ifdef DEBUG    
		cerr<<"3p\t+a\t"<<i<<"\t"<<sub3pPlus[i].s[0]<<"\t"<<sub3pPlus[i].s[1]<<"\t"<<sub3pPlus[i].s[2]<<"\t"<<sub3pPlus[i].s[3]<<"\t"<<(sub3pPlus[i].s[0]+sub3pPlus[i].s[1]+sub3pPlus[i].s[2]+sub3pPlus[i].s[3])<<endl;
#endif
		aT = sub3pMinus[i].s[0]-background.s[0];
		cT = sub3pMinus[i].s[1]-background.s[1];
		gT = sub3pMinus[i].s[2]-background.s[2];
		tT = sub3pMinus[i].s[3]-background.s[3];
#ifdef DEBUG    
		cerr<<"3p\t-f\t"<<i<<"\t"<<aT<<"\t"<<cT<<"\t"<<gT<<"\t"<<tT<<"\t"<<(aT+cT+gT+tT)<<endl;
#endif

		sub3pMinus[i].s[0] = 0.25+aT;
		sub3pMinus[i].s[1] = 0.25+cT;
		sub3pMinus[i].s[2] = 0.25+gT;
		sub3pMinus[i].s[3] = 0.25+tT;
#ifdef DEBUG    
		cerr<<"3p\t-a\t"<<i<<"\t"<<sub3pMinus[i].s[0]<<"\t"<<sub3pMinus[i].s[1]<<"\t"<<sub3pMinus[i].s[2]<<"\t"<<sub3pMinus[i].s[3]<<"\t"<<(sub3pMinus[i].s[0]+sub3pMinus[i].s[1]+sub3pMinus[i].s[2]+sub3pMinus[i].s[3])<<endl;
#endif
	    }
	}else{
	    cerr << "Unable to open composition file "<<compFile<<endl;
	    return 1;
	}

    }else{//else no comp file
	distFromEnd=0;
    }


    // string bamfiletopen = string(argv[argc-1]);//bam file

    string inFile       = string(argv[argc-1]);//fasta file
    string fastaFile;//    = string(argv[argc-1]);//fasta file
    string fastaFileFai;//    = string(argv[argc-1]);//fasta file
    
    double maxProbP = 1.0;
    double maxProbM = 1.0;


    for(int i=0;i<2*distFromEnd;i++){
	maxProbP=maxProbP*maxSubrates(sub5pPlus[i] );
    }
    for(int i=0;i<2*distFromEnd;i++){
	maxProbM=maxProbM*maxSubrates(sub5pMinus[i]);	
    }

    for(int i=0;i<2*distFromEnd;i++){
	maxProbP=maxProbP*maxSubrates(sub3pPlus[i] );		
    }
    
    for(int i=0;i<2*distFromEnd;i++){
	maxProbM=maxProbM*maxSubrates(sub3pMinus[i]);	
    }


#ifdef DEBUG    
    cerr<<"Maxprob.+\t"<<maxProbM<<"\tMaxprob.-\t"<<maxProbP<<endl;
    cerr<<"----------------"<<endl;
    cerr<<"5p + misincorporation"<<endl;
    for(int i=0;i<2*distFromEnd;i++){
	cerr<<i<<"\t"<<sub5pPlus[i].s[0]<<"\t"<<sub5pPlus[i].s[1]<<"\t"<<sub5pPlus[i].s[2]<<"\t"<<sub5pPlus[i].s[3]<<"\t"<<(sub5pPlus[i].s[0]+sub5pPlus[i].s[1]+sub5pPlus[i].s[2]+sub5pPlus[i].s[3])<<endl;
    }
    cerr<<"----------------"<<endl;
    cerr<<"5p - misincorporation"<<endl;
    for(int i=0;i<2*distFromEnd;i++){
	cerr<<i<<"\t"<<sub5pMinus[i].s[0]<<"\t"<<sub5pMinus[i].s[1]<<"\t"<<sub5pMinus[i].s[2]<<"\t"<<sub5pMinus[i].s[3]<<"\t"<<(sub5pMinus[i].s[0]+sub5pMinus[i].s[1]+sub5pMinus[i].s[2]+sub5pMinus[i].s[3])<<endl;
    }
    cerr<<"----------------"<<endl;
    cerr<<"3p + misincorporation"<<endl;
    for(int i=0;i<2*distFromEnd;i++){
	cerr<<i<<"\t"<<sub3pPlus[i].s[0]<<"\t"<<sub3pPlus[i].s[1]<<"\t"<<sub3pPlus[i].s[2]<<"\t"<<sub3pPlus[i].s[3]<<"\t"<<(sub3pPlus[i].s[0]+sub3pPlus[i].s[1]+sub3pPlus[i].s[2]+sub3pPlus[i].s[3])<<endl;
    }
    cerr<<"----------------"<<endl;
    cerr<<"3p - misincorporation"<<endl;
    for(int i=0;i<2*distFromEnd;i++){
	cerr<<i<<"\t"<<sub3pMinus[i].s[0]<<"\t"<<sub3pMinus[i].s[1]<<"\t"<<sub3pMinus[i].s[2]<<"\t"<<sub3pMinus[i].s[3]<<"\t"<<(sub3pMinus[i].s[0]+sub3pMinus[i].s[1]+sub3pMinus[i].s[2]+sub3pMinus[i].s[3])<<endl;
    }
    cerr<<"----------------"<<endl;
#endif

    //return 1;
    bool isGzipped=false;
    string tempfile_ = tmpDir+"/tmpfileXXXXXX";
    IndexedGenome* genome=NULL;

    FastQParser * fp=NULL;
    FastQObj * fo;

    uint64_t genomeLength=0;
    vector<chrinfo> chrFound;

    
    
    typedef map<string,faidx1_t>::iterator it_type;
    
    if(fastqMode){

	if(outBAMb){
	    cerr<<"Error: cannot specify output BAM for fastq mode."<<endl;
	    return 1;
	}
	if(outFastagzb){
	    cerr<<"Error: cannot specify output zipped fastq for fastq mode."<<endl;
	    return 1;
	}
	if(gcBiasB){
	    cerr<<"Error: cannot specify GC bias for fastq mode."<<endl;
	    return 1;
	}

	fp = new FastQParser(inFile,false);

    }else{
	if(strEndsWith(inFile,".gz")){
	

	    char tmpname[ tempfile_.size()];
	    strcpy(tmpname,tempfile_.c_str());
	    int fdtemp = mkstemp(tmpname);
	    if(fdtemp == -1){
		cerr<<"Cannot create temp file using pattern: "<<tempfile_<<", either choose a different temp dir or unzip the file"<<endl;
		return 1;
	    }
	    tempfile_ = string(tmpname);

	    cerr<<"File "<<inFile<<" is zipped, trying to write to temp. file: "<<tempfile_<<endl;

	    fastaFile    = tmpname;
	    fastaFileFai = inFile+".fai";	

	    igzstream fastagzfd;
	    string linetmp;

	    ofstream myFiletmp;	
	    myFiletmp.open(tmpname);

	    fastagzfd.open(inFile.c_str(), ios::in);

	    if (fastagzfd.good()){
		while ( getline (fastagzfd,linetmp)){
		    myFiletmp << linetmp << endl;
		}
	    }
	    myFiletmp.close();

	    cerr<<"unzipping is done"<<endl;

	    //fastaFile
	    isGzipped=true;
	}else{
	    fastaFile   =inFile;
	    fastaFileFai=inFile+".fai";	
	}
	genome=new IndexedGenome(fastaFile.c_str(),fastaFileFai.c_str(),uppercase);
	cerr<<"Mapped "<<fastaFile<<" into memory"<<endl;

	
	for(it_type iterator = genome->name2index.begin(); iterator != genome->name2index.end(); iterator++) {
#ifdef DEBUG    
	    cerr<<iterator->first<<"\t"<<iterator->second.len<<"\t"<<genomeLength<<endl;
#endif

	    chrinfo toadd;

	    toadd.name          = iterator->first;
	    toadd.startIndexChr = genomeLength+1;
	    toadd.length        = iterator->second.len;
	    toadd.endIndexChr   = genomeLength+iterator->second.len;
	    genomeLength       += iterator->second.len;

	    chrFound.push_back(toadd);
	}

    }//if not fastq
    

    //if circ was selected, try to check if the chr exists
    if(circ){
	faidx1_t faidxForName=genome->name2index[ circName ];			
	if(faidxForName.offset == 0){
	    cerr<<"ERROR: Cannot find chromosome "<<circName<<" which is specified via --circ"<<endl; 
	    return 1;
	}
	circOffset=faidxForName.offset;
    }
    

#ifdef DEBUG    
    cerr<<"genomeLength "<<genomeLength<<endl;
#endif

    double theoMeanGC=0.4;
    double theoStdvGC=0.1;
    double maxNormGC = 1;

    if(gcBiasB){
	cerr<<"Computing average GC content"<<endl;


	unsigned int f       = 0;

	vector<double> gcVal;

	while(f<1000000){
	    int length      = 50;
	    unsigned int gcCount = 0;
	    unsigned int atCount = 0;
	    if(length<distFromEnd)//if the length of the fragment is lesser than the distance from end, creating a slight bias against short fragments but acceptable one if distFromEnd is small enough
		continue;

	
	    //	    uint64_t idx;//         = randomInt(distFromEnd,int(chrall.size())-length-distFromEnd);
	    uint64_t coord;
	    bool found=false;
	    //string temp     = chrall.substr(idx-distFromEnd,length+2*distFromEnd);
	    string temp="";

	    while(!found){
		coord =randGenomicCoord(genomeLength);
	    

		for(unsigned int i=0;i<chrFound.size();i++){     
		    if( (chrFound[i].startIndexChr+distFromEnd) <= coord 
			&& 
			coord <= (chrFound[i].endIndexChr-length-distFromEnd+1)){
			found=true;
			//idx = coord-chrFound[i].startIndexChr;
			faidx1_t faidxForName=genome->name2index[ chrFound[i].name ];			
			temp=genome->fetchSeq(&faidxForName,coord-chrFound[i].startIndexChr-distFromEnd, length+2*distFromEnd);			
			break;
		    }
		}

		if(found){		
		    break;
		}
	    }

	    bool plusStrand;
	    if(noRev){
		plusStrand = true;
	    }else{
		plusStrand = randomBool();
	    }

	    if(!isResolvedDNAstring(temp))
		continue;
	    
	    if(!plusStrand){
		temp    = reverseComplement(temp);
	    }else{	   
	    }
	    
	    for(unsigned int i=0;i<temp.size();i++){
		//cout<<temp[i]<<endl;
		if(isResolvedDNA(temp[i])){		    
		    if(temp[i] == 'G' ||
		       temp[i] == 'C' ){
			gcCount++;
		    }else{
			atCount++;
		    }
		}	    
	    }
	    //cout<<(double(gcCount)/double(gcCount+atCount))<<endl;
	    gcVal.push_back( double(gcCount)/double(gcCount+atCount) );
	    f++;	
	}

	theoMeanGC  = meanNorm( gcVal                            );
	theoStdvGC  = stdevNorm(gcVal,      theoMeanGC           );
	maxNormGC   = pdfNorm(  theoMeanGC, theoMeanGC,theoStdvGC);

	cerr<<"GC content: mean:"<<theoMeanGC<<" stdev:"<<theoStdvGC<<endl;
    }

#ifdef DEBUG    
    cerr<<"done computing mean GC "<<endl;
#endif



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
    ogzstream outFastagzfp;

    if(outFastagzb){
	outFastagzfp.open(outFastagz.c_str(), ios::out);
	if(!outFastagzfp.good()){
	    cerr<<"Cannot write to file "<<outFastagz<<endl; 
	    return 1; 
	}
    }

    BamWriter writer;
    if(outBAMb){
	SamHeader header;
	string pID          = "fragSim";   
	string pName        = "fragSim";   
	string pCommandLine = "";
	for(int i=0;i<(argc);i++){
	    pCommandLine += (string(argv[i])+" ");
	}
	putProgramInHeader(&header,pID,pName,pCommandLine,returnGitHubVersion(string(argv[0]),".."));
	RefVector references;

	if(produceUnCompressedBAM) 
	    writer.SetCompressionMode(BamWriter::Uncompressed);
	
	if ( !writer.Open(outBAM,header,references) ) {
	    cerr<<"Cannot open write to BAM output file "<<outBAM<<endl;
	    return 1;
	}

    }

    unsigned int   f       = 0;
  
    if(fastqMode){
	while(true){

	    string def;
	    string seq;
	    string qal;
	    if(!fp->hasData())
		break;
	    fo  = fp->getData();
	    def = *(fo->getID());
	    seq = *(fo->getSeq());
	    qal = *(fo->getQual());

	    int length      = selectFragmentLength(fileSizeFragB,
						   fileSizeFragBFreq,
						   specifiedLength,
						   specifiedScale,
						   sizeFragList,
						   fragLengthFreq,
						   sizeFragments,
						   generator,
						   distribution);

	    cout<<def<<endl<<seq.substr(0,length)<<endl<<"+"<<endl<<qal.substr(0,length)<<endl;
	    f++;
	}
    }else{
     
	while(f<nFragments){


	    int length      = selectFragmentLength(fileSizeFragB,
						   fileSizeFragBFreq,
						   specifiedLength,
						   specifiedScale,
						   sizeFragList,
						   fragLengthFreq,
						   sizeFragments,
						   generator,
						   distribution);
	
	    if(length<distFromEnd)//if the length of the fragment is lesser than the distance from end, creating a slight bias against short fragments but acceptable one if distFromEnd is small enough
		continue;
	    if(length<minimumFragSize)
		continue;

	    if(length>maximumFragSize)
		continue;

	
	    uint64_t idx;//         = randomInt(distFromEnd,int(chrall.size())-length-distFromEnd);
	    uint64_t coord;
	    bool found=false;
	    //string temp     = chrall.substr(idx-distFromEnd,length+2*distFromEnd);
	    string temp="";
	    string deflineToPrint;// = defline+":";
	    while(!found){
		coord =randGenomicCoord(genomeLength);
	    
#ifdef DEBUG    
		cerr<<"coord "<<coord<<endl;
#endif

		for(unsigned int i=0;i<chrFound.size();i++){ 
		    
		    if( (chrFound[i].startIndexChr+distFromEnd) <= coord 
			&& 
			//(coord <= (chrFound[i].endIndexChr-length-distFromEnd+1))
			(coord <= (chrFound[i].endIndexChr-distFromEnd+1))
			){
			idx = coord-chrFound[i].startIndexChr;
			faidx1_t faidxForName=genome->name2index[ chrFound[i].name ];

			if(circ){//we consider some references to be circular
			    if( circOffset==faidxForName.offset){//it is the circular chromosome
				found=true;
				if((coord <= (chrFound[i].endIndexChr-length-distFromEnd+1))){//distance to end is reasonable, no need to circularize
				    
				}else{//need to circularize
				    string temp1=genome->fetchSeq(&faidxForName,
								  coord-chrFound[i].startIndexChr-distFromEnd,
								  chrFound[i].endIndexChr-(coord-distFromEnd));

				    int l2= (length+2*distFromEnd-(chrFound[i].endIndexChr-(coord-distFromEnd)));
				    
				    string temp2=genome->fetchSeq(&faidxForName,
				     				  0,
				     				  l2);				    
				    temp=temp1+temp2;

				    deflineToPrint = chrFound[i].name+":";					  
				    break;				    
				}
			    }else{ // not the circular chromosome
				if((coord <= (chrFound[i].endIndexChr-length-distFromEnd+1))){//distance to end is reasonable
				    found=true;
				}else{//too close to the end
				    break;
				}
			    }

			}else{
			    found=true;
			}
			
#ifdef DEBUG    
			cerr<<"found idx="<<idx<<" name="<<chrFound[i].name<<endl;
#endif

			temp=genome->fetchSeq(&faidxForName,coord-chrFound[i].startIndexChr-distFromEnd, length+2*distFromEnd);
			deflineToPrint = chrFound[i].name+":";
			break;
		    }
		}//end for each chr


		if(found){		
		    break;
		}

#ifdef DEBUG    
		cerr<<"coord not found "<<coord<<endl;
#endif

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
#ifdef DEBUG    	    
	    cerr<<"t:"<<temp<<"#\tl="<<length<<endl;	
#endif

	    if(!plusStrand){
#ifdef DEBUG    	    
		cerr<<"-"<<endl;	    
#endif
		temp    = reverseComplement(temp);
		preFrag = temp.substr(                 0 , distFromEnd);
		frag    = temp.substr(       distFromEnd ,      length);
		posFrag = temp.substr(length+distFromEnd , distFromEnd);

		deflineToPrint = deflineToPrint+"-:"+stringify(idx)+":"+stringify(idx+length)+":"+stringify(length);
	    }else{
#ifdef DEBUG    
		cerr<<"+"<<endl;
#endif
		preFrag = temp.substr(                 0 , distFromEnd);
		frag    = temp.substr(       distFromEnd ,      length);
		posFrag = temp.substr(length+distFromEnd , distFromEnd);

		deflineToPrint = deflineToPrint+"+:"+stringify(idx)+":"+stringify(idx+length)+":"+stringify(length);
	    }

#ifdef DEBUG    	    
	    cerr<<"t:"<<temp<<"#"<<endl;	
#endif

	    if(tagb && !uniqTags){
		deflineToPrint=deflineToPrint+tag;	    
	    }

	    if(gcBiasB){

		unsigned int gcCount=0;
		unsigned int atCount=0;

		for(unsigned int i=0;i<temp.size();i++){		
		    if(isResolvedDNA(temp[i])){		    
			if(temp[i] == 'G' ||
			   temp[i] == 'C' ){
			    gcCount++;
			}else{
			    atCount++;
			}
		    }	    
		}
		double gcContent =  double(gcCount)/double(gcCount+atCount);
		double pSurv     = 1/(1+exp(gcBias*(gcContent-theoMeanGC)));
		pSurv            = pSurv * (pdfNorm(  gcContent, theoMeanGC, theoStdvGC)/maxNormGC);

		double rP = randomProb();
		bool killed = (rP>pSurv);
		//	    cout<<gcContent<<"\t"<<pSurv<<"\t"<<rP<<"\t"<<killed<<endl;
		if(killed){//killed do not increase f
		    continue;
		}
	    }

	    // cerr<<"---------"<<endl;
	    // cerr<<temp<<endl;
	    // cerr<<preFrag<<endl;
	    // cerr<<frag<<endl;
	    // cerr<<posFrag<<endl;


	    if(!compFileSpecified){ //no composition file, just produce the sequences	    
#ifdef DEBUG    
		cerr<<"no comp file"<<endl;
#endif

		if(gcBiasB){


		    unsigned int gcCount=0;
		    unsigned int atCount=0;

		    for(unsigned int i=0;i<temp.size();i++){		
			if(isResolvedDNA(temp[i])){		    
			    if(temp[i] == 'G' ||
			       temp[i] == 'C' ){
				gcCount++;
			    }else{
				atCount++;
			    }
			}	    
		    }
		    double gcContent =  double(gcCount)/double(gcCount+atCount);
		    double pSurv     = 1/(1+exp(gcBias*(gcContent-theoMeanGC)));
		    pSurv            = pSurv * (pdfNorm(  gcContent, theoMeanGC, theoStdvGC)/maxNormGC);

		    double rP = randomProb();
		    bool killed = (rP>pSurv);
		    //	    cout<<gcContent<<"\t"<<pSurv<<"\t"<<rP<<"\t"<<killed<<endl;
		    if(killed){//killed do not increase f
			continue;
		    }

		
		}


		if(uniqTags){
		    checkUniqueness(&deflineToPrint,
				    &fragmentID2Count,
				    tagb ,
				    tag);
		}
	
		if(outFastagzb){
		    outFastagzfp<<">"<<deflineToPrint<<"\n"<<frag<<endl;
		}else{
		    if(outBAMb){
		    
			BamAlignment al;
			al.Name  = deflineToPrint;
			al.SetIsMapped (false);
			al.SetIsPaired(false);
			al.AlignmentFlag =  flagSingleReads;
			//     if(!al.AddTag("XI", "Z",string(index1S[i])) ) {cerr<<"Internal error, cannot add tag"<<endl; return 1; } 
			al.QueryBases  = frag;
			al.Qualities   = string(frag.size(), '!');

			writer.SaveAlignment(al);

		    }else{
			cout<<">"<<deflineToPrint<<"\n"<<frag<<endl;
		    }
		}


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

		    probAcc = probAcc/maxProbP;

#ifdef DEBUG    
		    cerr<<"+ strand comp file p[acc]="<<probAcc<<"\t"<<distFromEnd<<endl;
#endif

		    if(randomProb()<probAcc){

#ifdef DEBUG    
			cerr<<"accepted "<<probAcc<<endl;
#endif


			if(uniqTags){
			    checkUniqueness(&deflineToPrint,
					    &fragmentID2Count,
					    tagb ,
					    tag);
			}
	
		
			if(outFastagzb){
			    outFastagzfp<<">"<<deflineToPrint<<"\n"<<frag<<endl;
			}else{
			    if(outBAMb){

				BamAlignment al;
				al.Name  = deflineToPrint;
				al.SetIsMapped (false);
				al.SetIsPaired(false);
				al.AlignmentFlag =  flagSingleReads;
				//     if(!al.AddTag("XI", "Z",string(index1S[i])) ) {cerr<<"Internal error, cannot add tag"<<endl; return 1; } 
				al.QueryBases  = frag;
				al.Qualities   = string(frag.size(), '!');

				writer.SaveAlignment(al);

			    }else{
				cout<<">"<<deflineToPrint<<"\n"<<frag<<endl;
			    }
			}
		    
			f++;
		    }else{
			//discard

#ifdef DEBUG    
			cerr<<"rejected "<<probAcc<<endl;
#endif
			continue;
		    }

		}else{//minus strand

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

		    probAcc = probAcc/maxProbM;

#ifdef DEBUG    
		    cerr<<"- strand comp file p[acc]="<<probAcc<<"\t"<<distFromEnd<<endl;
#endif

		    if(randomProb()<probAcc){
			//cout<<deflineToPrint<<"\n"<<frag<<endl;
#ifdef DEBUG    
			cerr<<"accepted "<<probAcc<<endl;
#endif

			if(uniqTags){
			    checkUniqueness(&deflineToPrint,
					    &fragmentID2Count,
					    tagb ,
					    tag);
			}

			if(outFastagzb){
			    outFastagzfp<<">"<<deflineToPrint<<"\n"<<frag<<endl;
			}else{
			    if(outBAMb){
			    
				BamAlignment al;
				al.Name  = deflineToPrint;
				al.SetIsMapped (false);
				al.SetIsPaired(false);
				al.AlignmentFlag =  flagSingleReads;
				//     if(!al.AddTag("XI", "Z",string(index1S[i])) ) {cerr<<"Internal error, cannot add tag"<<endl; return 1; } 
				al.QueryBases  = frag;
				al.Qualities   = string(frag.size(), '!');

				writer.SaveAlignment(al);

			    }else{
				cout<<">"<<deflineToPrint<<"\n"<<frag<<endl;
			    }
			}

			f++;
		    }else{
			//discard
#ifdef DEBUG    
			cerr<<"rejected "<<probAcc<<endl;
#endif
			continue;
		    }


		}//minus strand
	    
	    }

	    if( (f-1)%100000 == 0 && (f-1)!=0){
		cerr<<"Produced "<<thousandSeparator(f-1)<<" out of "<< thousandSeparator(nFragments) <<" sequences"<<endl;
	    }

	}
    }

    if(outFastagzb){
	outFastagzfp.close();       
    }

    if(outBAMb){
	writer.Close();
    }


    if(isGzipped){
	cerr<<"removing temp file:"<<tempfile_<<endl;	   
	int rmrt=remove(tempfile_.c_str());
	if(rmrt == -1){
	    cerr<<"ERROR: Cannot remove temp file:"<<tempfile_<<endl;
	    return 1;
	}
    }
    cerr<<"Program "<<argv[0]<<" terminated succesfully, wrote "<<f<<" sequences"<<endl;
    return 0;
}


#include <iostream>
#include <fstream>
#include <gzstream.h>


#include <api/SamHeader.h>
#include <api/BamMultiReader.h>
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <api/BamAux.h>
#include "PutProgramInHeader.h"

#include "FastQParser.h"
#include "utils.h"

using namespace std;

const uint32_t flagSingleReads =  4; // 00000100
const uint32_t flagFirstPair   = 77; // 01001101
const uint32_t flagSecondPair  =141; // 10001101

typedef struct { 
    double s[12];
 } substitutionRates;

int main (int argc, char *argv[]) {

    string options_adapter_F_BAM="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG";
    string options_adapter_S_BAM="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTT";

    int desiredLength           =100;
    bool addInName              =false;
    // string fwdrfName            ="/dev/stdout";
    // string revrfName            ="/dev/null";

    string outBAM                    ;
    bool   outBAMb              =false;
    bool   outBAMp              =false;

    bool   outFastagzb          =false      ;
    string outFastagzfwd                    ;
    string outFastagzrev        ="/dev/null";

    string outArt              ="/dev/null";
    bool   outArts             =false;
    bool   outArtp             =false;
    bool   outArtb             =false;

    bool produceUnCompressedBAM =false;
    bool     tagb             =false;
    string   tag              = "";

    const string usage=string("\t"+string(argv[0])+
                              " [options]  [BAM/fasta file]"+"\n\n"+
			      " This program reads a fasta file containing aDNA fragments and\n"+
			      " splits them into two records, one containing the forward read\n"+
			      " and the second containing the reverse read (-fr,-rr)\n"+
			      " or into a single for single-end mode (-fr)"+
			      
			      "\n\tOptions:\n"+ 
			      "\n\t\tI/O options:\n"+ 
			      
			      "\t\t"+"-arts\t[out]" +"\t\t"+"Output single-end reads as ART (unzipped fasta) (Default: "+outArt+")"+"\n"+
			      "\t\t"+"-artp\t[out]" +"\t\t"+"Output reads as ART (unzipped fasta) (Default: "+outArt+")"+"\n"+
			      "\t\t"+""            +"\t\t\t"+"with wrap-around for paired-end mode"+"\n"+
		
			      "\t\t"+"-fr\t[out fwdr]" +"\t"+"Output forward read as zipped fasta (Default: "+outFastagzfwd+")"+"\n"+
			      "\t\t"+"-rr\t[out rwdr]" +"\t"+"Output reverse read as zipped fasta (Default: "+outFastagzrev+")"+"\n"+
			      "\t\t"+"-bs\t"+"[BAM out]"   +"\t"+"Read BAM and write output as a single-end BAM (Default: fasta)"+"\n"+
			      "\t\t"+"-bp\t"+"[BAM out]"   +"\t"+"Read BAM and write output as a single-end BAM (Default: fasta)"+"\n"+
			      "\t\t"+"-u"                 +"\t\t\t"+"Produce uncompressed BAM (good for unix pipe)"+"\n"+
			      "\n"+
			      
                              "\t"+"-f\t[seq]" +"\t\t\t"+"Adapter that is observed after the forward read (Default:  "+options_adapter_F_BAM.substr(0,10)+"...)"+"\n"+
                              "\t"+"-s\t[seq]" +"\t\t\t"+"Adapter that is observed after the reverse read (Default:  "+options_adapter_S_BAM.substr(0,10)+"...)"+"\n"+
                              "\t"+"-l\t[length]" +"\t\t"+"Desired read length  (Default:  "+stringify(desiredLength)+")"+"\n"+
			      "\t"+"-name" +"\t\t\t\t"+"Append BAM tags or to deflines if adapters are added (Default:  "+booleanAsString(tagb)+")"+"\n"+
			      "\t"+"-tag\t" +"[tag]\t\t\t"+"Append this string to deflines or BAM tags (Default:  "+booleanAsString(tagb)+")"+"\n"+

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

    for(int i=1;i<(argc-1);i++){ 
	
	if(string(argv[i]) == "-u"  ){
	    produceUnCompressedBAM=true;
	    continue;
	}

        if(string(argv[i]) == "-arts" ){
            outArt  = string(argv[i+1]);
            i++;
            outArtb = true;
            outArts = true;
            continue;
        }

        if(string(argv[i]) == "-artp" ){
            outArt  = string(argv[i+1]);
            i++;
            outArtb = true;
            outArtp = true;
            continue;
        }

        if(string(argv[i]) == "-bs" ){
            outBAM  = string(argv[i+1]);
            i++;
            outBAMb = true;
            outBAMp = false;
            continue;
        }

        if(string(argv[i]) == "-bp" ){
            outBAM  = string(argv[i+1]);
            i++;
            outBAMb = true;
            outBAMp = true;
            continue;
        }

	if(string(argv[i]) == "-l" ){
	    desiredLength          = destringify<int>(argv[i+1]);
	    i++;
	    continue;
	}

	if(string(argv[i]) == "-f" ){
	    options_adapter_F_BAM          = string(argv[i+1]);
	    i++;
	    continue;
	}

	if(string(argv[i]) == "-s" ){
	    options_adapter_S_BAM          = string(argv[i+1]);
	    i++;
	    continue;
	}

	if(string(argv[i]) == "-fr" ){
	    outFastagzfwd       = string(argv[i+1]);
	    i++;
	    outFastagzb         = true;
	    continue;
	}

	if(string(argv[i]) == "-rr" ){
	    outFastagzrev       = string(argv[i+1]);
	    i++;
	    outFastagzb         = true;
	    continue;
	}

	if(string(argv[i]) == "-name" ){
	    addInName=true;
	    continue;
	}

	if(string(argv[i]) == "-tag" ){
	    tagb=true;
	    tag       = string(argv[i+1]);
	    i++;
	    continue;
	}

        cerr<<"Error: unknown option "<<string(argv[i])<<endl;
        return 1;
    }

    if(outFastagzb && outFastagzfwd.empty()){
        cerr<<"Error: need to specify the output forward read "<<endl;
        return 1;
    }

    if(outFastagzb && outBAMb){
        cerr<<"Error: cannot produce both fastq and BAM"<<endl;
        return 1;
    }

    if(outFastagzb && outArtb){
        cerr<<"Error: cannot produce both fastq and ART output"<<endl;
        return 1;
    }

    if(outBAMb && outArtb){
        cerr<<"Error: cannot produce both BAM and ART output"<<endl;
        return 1;
    }

    string infileName = string(argv[argc-1]);
    
    //Fastq
    FastQParser * fp=NULL;
    FastQObj * fo;
    ogzstream outFastaFgzfp;
    ogzstream outFastaRgzfp;

    ofstream   outARTfp;


    //BAM
    BamReader reader;
    BamWriter writer;



    //INPUT OUTPUT
    if(outBAMb){

	if ( reader.Open(infileName) ) {
	    SamHeader  myHeader=reader.GetHeader();
	    SamProgram sp;
   
	    string pID          = "adptSim";   
	    string pName        = "adptSim";   
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
	    cerr << "Could not open input BAM file  "<<infileName<< endl;
	    return 1;   
	}

    }else{
	fp = new FastQParser(infileName,true);

	if( outFastagzb){	    
	    outFastaFgzfp.open(outFastagzfwd.c_str(), ios::out);
	    if(!outFastaFgzfp.good()){       cerr<<"Cannot write to file "<<outFastagzfwd<<endl; return 1; }

	    outFastaRgzfp.open(outFastagzrev.c_str(), ios::out);
	    if(!outFastaRgzfp.good()){       cerr<<"Cannot write to file "<<outFastagzrev<<endl; return 1; }
	}else{
	    if( outArtb){	    		
		outARTfp.open(outArt.c_str(), ios::out);
		if(!outARTfp.good()){ cerr<<"Cannot write to file "<<outArt<<endl; return 1; }		    
	    }
	    //nothing
	}
    }

    unsigned int   f       = 0;
    while(true){

	string namef ;
	string namer ;
	string seqf  ;
	string seqr  ;

        BamAlignment alf;
        BamAlignment alr;
	BamAlignment al;

	if(outBAMb){
	    if(!reader.GetNextAlignment(al))
		break;
	    namef = al.Name;
	    seqf  = al.QueryBases;
	    namer = namef;
	}else{
	    if(!fp->hasData())
		break;
	    fo    = fp->getData();
	    namef = *(fo->getID());
	    seqf  = *(fo->getSeq());	    
	    namer = namef+"_r";
	}
	transform(seqf.begin(), seqf.end(),seqf.begin(), ::toupper);
	seqr  = reverseComplement(seqf);
	//cerr<<"seqr "<<seqf<<endl;



	bool addedAdapter = false;
	bool addedRandom  = false;

	//forward
	if(int(seqf.size()) < desiredLength){
	    seqf=seqf+options_adapter_F_BAM;
	    seqf=seqf.substr(0,desiredLength);
	    addedAdapter=true;
	}else{
	    seqf=seqf.substr(0,desiredLength);	    
	}
	
	//if even with the adapter does not add up to length, add random seq
	if(int(seqf.size()) < desiredLength){
	    seqf=seqf+randomDNASeq(desiredLength-int(seqf.size()));
	    addedRandom=true;
	}

	//rev
	if(int(seqr.size()) < desiredLength){
	    seqr=seqr+options_adapter_S_BAM;
	    seqr=seqr.substr(0,desiredLength);
	}else{
	    seqr=seqr.substr(0,desiredLength);
	}
	
	if(int(seqr.size()) < desiredLength){
	    seqr=seqr+randomDNASeq(desiredLength-int(seqr.size()));
	}




	
	if(outBAMb){
            alf.QueryBases  = seqf;
	    alf.Qualities   = string(seqf.size(), '!');
            alr.QueryBases  = seqr;
	    alr.Qualities   = string(seqr.size(), '!');

	    alf.Name        = namef;
	    alr.Name        = namer;
	    
	    if(!outBAMp){
		alf.AlignmentFlag =  flagSingleReads;
	    }else{
		alf.AlignmentFlag =  flagFirstPair  ;
		alr.AlignmentFlag =  flagSecondPair ;
	    }

            if(addInName || tagb){
		string toaddtag="";
		if(addInName){
		    if(addedAdapter)
			toaddtag=toaddtag+"AD";	    
		    if(addedRandom)
			toaddtag=toaddtag+",RD";
		}

		if(tagb){
		    toaddtag=toaddtag+tag;	    
		}
		
                if(!alf.AddTag("XA", "Z",toaddtag) ) {
                    cerr<<"Internal error, cannot add tag to BAM alignment "<<alf.Name<<endl; 
                    return 1; 
                }

                if(!alr.AddTag("XA", "Z",toaddtag) ) {
                    cerr<<"Internal error, cannot add tag to BAM alignment "<<alr.Name<<endl; 
                    return 1; 
                }

            }
	    //add data and tags
	    
            writer.SaveAlignment(alf);
	    if(outBAMp)
		writer.SaveAlignment(alr);

	}else{//fasta
	    if(addInName  || tagb ){
		if(addInName){

		    if(addedAdapter)
			namef=namef+"_adpt";	    
		    if(addedRandom)
			namef=namef+"_rand";

		    if(addedAdapter)
			namer=namer+"_adpt";	    
		    if(addedRandom)
			namer=namer+"_rand";
		}

		if(tagb){
		    namef=namef+tag;	    
		    namer=namer+tag;	    
		}
		
	    }

	    if(outFastagzb){
		outFastaFgzfp<<namef<<endl<<seqf<<endl;
		outFastaRgzfp<<namer<<endl<<seqr<<endl;
	    }else{
		if(outArtb){
		    if(outArts)
			outARTfp<<namef<<endl<<seqf<<endl;
		    if(outArtp){
			seqr  = reverseComplement(seqr);
			outARTfp<<namef<<endl<<seqf<<seqr<<endl;
		    }
		}else{
		    cout<<namef<<endl<<seqf<<endl<<namer<<endl<<seqr<<endl;
		}
		
	    }
	}

	if(f%100000 == 0 && f!=0){
	    cerr<<"Produced "<<f<<" sequences"<<endl;
	}
	f++;
	
       
    }


    if(outFastagzb){
        outFastaFgzfp.close();
        outFastaRgzfp.close();
    }
    
    if(outBAMb){
        reader.Close();
        writer.Close();
    }

    if( outArtb){	    
	outARTfp.close();
    }

    cerr<<"Program "<<argv[0]<<" terminated succesfully, wrote "<<f<<" sequences"<<endl;



    return 0;
}


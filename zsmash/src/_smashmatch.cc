
#include <cstdio>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/tokenizer.hpp>
#include <boost/token_functions.hpp>
#include <exception>
#include <boost/timer/timer.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_histogram.h>

#include "semantic.h"
#include "config.h"

#define ZERO_ 1e-15

using namespace boost::program_options;

//-----------------------------------------------------

const string VERSION="v0.1 Copyright Ishanu Chattopadhyay 2016";
const string EMPTY_ARG_MESSAGE="Exiting. Type -h or --help for usage";

//-----------------------------------------------------

void err_(string s)
{
  cout << s << endl;
  exit(0);
};

//-----------------------------------------------------
ostream& operator << (ostream &out, map<unsigned int,unsigned int>& s)
{
  for(map<unsigned int,unsigned int>::iterator itr=s.begin();
      itr!=s.end();
      ++itr)
    out << itr->first << " " << itr->second << endl;
  return out;
}
//-----------------------------------------------------
ostream& operator << (ostream &out, vector<unsigned int>& s)
{
  for(vector<unsigned int>::iterator itr=s.begin();
      itr!=s.end();
      ++itr)
    out << *itr<< endl;
  return out;
}

//-----------------------------------------------------

PFSA read_mc(string filename, string TYPE__)
{
  connx aut;
  pitilde pit, Xpit;
  map < state, map < symbol, double > > pitilde_map;

  CONFIG modfile(filename);
  modfile.set_map<state,symbol,state>(aut,"CONNX");
  modfile.set_map<state,symbol,double>(pitilde_map,"PITILDE");

  PFSA *G;
  state c_state=0;
  for ( map < state, map < symbol,
	  double > >::iterator iti=pitilde_map.begin();
	iti != pitilde_map.end();
	iti++)
    {
      vector <double> vec_tmp;
      for (map <symbol, double>::iterator itj=iti->second.begin();
	   itj != iti->second.end();
	   itj++)
	vec_tmp.push_back(itj->second);
      pit[c_state++] = vec_tmp;
    }

  if(TYPE__=="XPFSA")
    {
      pitilde pit2;
      for (state  j=0; j < (int)aut.size();++j)
	pit2[j]= vector <double> (aut[0].size(),
				  1.0/(aut[0].size()+0.0));
      G = new PFSA(pit2,aut);
    }
  else
    G = new PFSA (pit,aut);
  G->set_Xpit(pit);

  return *G;
};

//-----------------------------------------------------
//-----------------------------------------------------

class lib_select_
{
protected:
  vector<matrix_dbl> H_;
  unsigned int ndata_;
  unsigned int nlib_;
  
  matrix_dbl LIB_;
  map<unsigned int,unsigned int> class_;

  size_t Nbins;
  
  map<unsigned int,matrix_dbl> cross_coeff;
  matrix_dbl class_prob_;
  vector<unsigned int> class_dec_;

public:

  //allocate Nbins
  lib_select_();
  lib_select_(vector<matrix_dbl>);


  void get_histogram(unsigned int nbins);
  void process_coeff_matrix();
  
  virtual matrix_dbl& get();
  vector<unsigned int>& class__();
  matrix_dbl& prob__();

};
//-----------------------------------------------------

lib_select_::lib_select_(vector<matrix_dbl> H)
{
  H_=H;
  nlib_=H.size();
  if(!H.empty())
    ndata_=H.begin()->size();
  else
    ndata_=0;
};

//-----------------------------------------------------
matrix_dbl& lib_select_::get()
{
  for(unsigned int i=0;i<H_.size();++i)
    for(matrix_dbl::iterator itr=H_[i].begin();
	itr!=H_[i].end();++itr)
      {
	double mn=10000000;
	for(map<unsigned int,double>::iterator itr_
	      =itr->second.begin();itr_!=itr->second.end();
	    ++itr_)
	  if(itr_->second<mn)
	    mn=itr_->second;
	      
	LIB_[itr->first][i]=mn;	  
      }

  //get histogram 
  get_histogram(10);
  process_coeff_matrix();

  return LIB_;
};
//-----------------------------------------------------
void lib_select_::process_coeff_matrix()
{
  for (map<unsigned int,matrix_dbl>::const_iterator itr
	 =cross_coeff.begin();
       itr!=cross_coeff.end();
       ++itr)
    for(matrix_dbl::const_iterator itr1
	  =itr->second.begin();
	itr1!=itr->second.end();
	++itr1)
      {
	class_prob_[itr->first][itr1->first]=1.0;
	for(map<unsigned int,double>::const_iterator
	      itr2=itr1->second.begin();
	    itr2!=itr1->second.end();
	    ++itr2)
	  if(itr1->first!=itr2->first) 
	    class_prob_[itr->first][itr1->first]*=itr2->second;
      }

  //normalization
  for(matrix_dbl::iterator itr
	=class_prob_.begin();
      itr!=class_prob_.end();
      ++itr)
    {
      double S=0.0;
      for(map<unsigned int,double>::iterator itr1
	    =itr->second.begin();
	  itr1!=itr->second.end();
	  ++itr1)
	S+=itr1->second;
      if(S>0.0)
	for(map<unsigned int,double>::iterator itr1
	      =itr->second.begin();
	    itr1!=itr->second.end();
	    ++itr1)
	  itr1->second/=S;	  
    }
  return;
}

//-----------------------------------------------------

void lib_select_::get_histogram(unsigned int nbins)
{
  Nbins=nbins;
  for (unsigned int samplei_=0;samplei_<ndata_;++samplei_)
    {
      map<unsigned int,gsl_histogram *> histograms_;
      map<unsigned int,gsl_histogram_pdf *> histograms_pdf;
      for (unsigned int libj_=0;libj_<nlib_;++libj_)
	{
	  double minval=0,maxval=1.0;
	  
	  histograms_[libj_] =  gsl_histogram_alloc (Nbins);
	  gsl_histogram_set_ranges_uniform (histograms_[libj_],
					    minval, maxval);
	  
	  for(map<unsigned int,double>::iterator itr
		=H_[libj_][samplei_].begin();
	      itr!=H_[libj_][samplei_].end();
	      ++itr)
	    gsl_histogram_increment (histograms_[libj_],
				     itr->second);
	}

      map<unsigned int, vector<double> > cum_prob_, rng_;
      for (map<unsigned int,
	     gsl_histogram *>::iterator itr
	     =histograms_.begin();
	   itr!=histograms_.end();
	   ++itr)
	{
	  //gsl_histogram_fprintf (stdout, itr->second, "%g", "%g");
	  histograms_pdf[itr->first]=gsl_histogram_pdf_alloc (Nbins);
	  gsl_histogram_pdf_init (histograms_pdf[itr->first],itr->second);
	  cum_prob_[itr->first]
	    =vector <double> (histograms_pdf[itr->first]->sum,
			      histograms_pdf[itr->first]->sum+Nbins);
	  /*
	  rng_[histograms_pdf[itr->first]]
	    =vector <double> (histograms_pdf[itr->first]->range,
			      histograms_pdf[itr->first]->range+Nbins);
	  */
	}

      for (map<unsigned int,
	     gsl_histogram *>::iterator itr
	     =histograms_.begin();
	   itr!=histograms_.end();
	   ++itr)
	for (map<unsigned int,
	       gsl_histogram *>::iterator itr1
	       =histograms_.begin();
	     itr1!=histograms_.end();
	     ++itr1)
	  {
	    double S=ZERO_, last_val=0.0;
	    for(size_t index=0;index<=Nbins;index++)
	      {
		S+=
		  cum_prob_[itr->first][index]
		  *(cum_prob_[itr1->first][index]
		    -last_val);
		last_val=cum_prob_[itr1->first][index];
	      }
	    cross_coeff[samplei_][itr->first][itr1->first]=S;
	  }

      //cout << "CROSS_COEFF ----------- " << endl;
      //cout << cross_coeff[samplei_] << endl;
      
      //cleanup
      for (map<unsigned int,
	     gsl_histogram *>::iterator itr
	     =histograms_.begin();
	   itr!=histograms_.end();
	   ++itr)
	{
	  gsl_histogram_free (itr->second);
	  gsl_histogram_pdf_free(histograms_pdf[itr->first]);
	}
    }
  return;
  
}

//-----------------------------------------------------
vector<unsigned int>& lib_select_::class__()
{
  class_dec_=vector<unsigned int>(ndata_,0);
  for(matrix_dbl::iterator itr=class_prob_.begin();
      itr!=class_prob_.end();
      ++itr)
    {
      double maxval=-1;
      unsigned int dec_=0;
      for(map<unsigned int,double>::iterator itr1=itr->second.begin();
	  itr1!=itr->second.end();
	  ++itr1)
	if(itr1->second > maxval)
	  {
	    maxval=itr1->second;
	    dec_=itr1->first;
	  }
      class_dec_[itr->first]=dec_;
    }

  return class_dec_;
}

//-----------------------------------------------------
matrix_dbl& lib_select_::prob__()
{  
  return class_prob_;
}
//-----------------------------------------------------
//-----------------------------------------------------
class lib_performance_: public lib_select_
{

public:

  lib_performance_();
  lib_performance_(vector<matrix_dbl> A): lib_select_(A){};

  matrix_dbl& get();
};
//-----------------------------------------------------
matrix_dbl& lib_performance_::get()
{
  return LIB_;
}
//-----------------------------------------------------

vector<option> ignore_numbers(vector<string>& args)
{
  vector<option> result;
  int pos = 0;
  while(!args.empty())
    {
      const auto& arg = args[0];
      double num;
      if(boost::conversion::try_lexical_convert(arg, num))
	{
	  result.push_back(option());
	  option& opt = result.back();

	  opt.position_key = pos++;
	  opt.value.push_back(arg);
	  opt.original_tokens.push_back(arg);

	  args.erase(args.begin());
	}
      else
	break;
    }

  return result;
}

void missing_alphabet_check(vector < Symbolic_string_ >& Svec)
{
  map<unsigned int,unsigned int> symset_;
  unsigned int MAX_ALPHABET_SIZE=0;
  for(unsigned int ni=0;ni<Svec.size();++ni)
    {
      set <symbol> symset;
      symbol_list_ thisstream=Svec[ni].get_symbol_list();
      for(unsigned int nj=0;nj<thisstream.size();++nj)
	symset.insert(thisstream[nj]);

      symset_[ni]=symset.size();
      if(symset.size()>MAX_ALPHABET_SIZE)
	MAX_ALPHABET_SIZE=symset.size();
    }
  for(map<unsigned int,unsigned int>::iterator itr=symset_.begin();
      itr!=symset_.end();
      ++itr)
    if(itr->second < MAX_ALPHABET_SIZE)
      {
	symbol_list_ upd_data=Svec[itr->first].get_symbol_list();
	Svec[itr->first]=Symbolic_string_(upd_data,MAX_ALPHABET_SIZE);
      }
  return;
};

//-----------------------------------------------------
//-----------------------------------------------------

int main(int argc, char *argv[])
{
  unsigned int len=100000;
  int NUM_EACH=100;
  unsigned int MODEL_RUN=100;
  unsigned int EVALUATION_DEPTH=8;
  vector <double> partition,partition_L;
  bool VERBOSE_=false;
  bool TIMER=true;
  bool DETAILED_OUTPUT=true;
  bool SAVE_PROB_=true;
  bool SAVE_CLASS_=true;
  
  string configfile="config.cfg", 
    datafile="data.txt",
    outfile="out",runfile="";
  vector<string> datafile_L,datafile_S;
  
  string DATA_DIR="across";
  string DATA_TYPE="continuous";
  string DATA_DIR_L="";
  string DATA_TYPE_L="";

  vector <bool> LIB_TYPE;
  
  options_description infor( "Program information");
  infor.add_options()
    ("help,h", "print help message.")
    ("version,V", "print version number");

  options_description usg( "Usage");
  usg.add_options()
    ("configfile,c",value<string>(),
     "config file [default: config.cfg]")
    ("datafile,f",value<vector< string> >()->multitoken(), "data file")
    ("libdatafile,F",value<vector<string> >()->multitoken(),
     "library data files")
    ("datalen,x",value< unsigned int >(),
     "data length max for model generated sequence")
    ("numeach,n",value< int >(),
     "number of repeats")
    ("datatype,T",value<string>(), "data type: continous or symbolic")
    ("datadir,D",value<string>(), "data direction: row or column")
    ("partition,P",value<vector<double> >()->multitoken(), "partition")
    ("libtype,L",value<vector<bool> >()->multitoken(),
     "specification if models or sequences are provided as library elements. True-> sequence, false -> file with model names. [true]")
    ("modelrun,m",value<unsigned int>(),
     "No. of model runs generated for comparison if models are used")
    ("timer,t",value< bool >(), "display timer [1 (true)] ");

  options_description outputopt( "Output options");
  outputopt.add_options()
    ("class,C",value<bool>(),
     "output major class [default: true]")
    ("prob,X",value<bool>(),
     "output class probabilities [default: true]")
    ("outfile,o",value<string>(),
     "output file prefix [default: out]")
    ("detailed_output,d",value< bool >(),
     "produce detailed output [true]")
    ("runfile,R",value< string >(),
     "filename to save  stream from machine run [default: - ]")
    ("verbose,v",value< bool >(), "verbosity [default: false] ");

  string tmpstr_="";
  options_description desc( (tmpstr_+"Library Match Via Data Smashing"+"\n"+VERSION).c_str());
  desc.add(infor).add(usg).add(outputopt);

  positional_options_description p;
  variables_map vm;

  if (argc == 1)
    err_(EMPTY_ARG_MESSAGE);
  try
    {
      store(command_line_parser(argc, argv)
	    .extra_style_parser(&ignore_numbers)
	    .options(desc)
	    .run(), vm);
      notify(vm);
    } 
  catch (std::exception &e)
    {
      cout << endl << e.what() 
	   << endl << desc << endl;
      return 1;
    }
  if (vm.count("help"))
    {
      cout << desc;
      err_("");
    }
  if (vm.count("version"))
    err_(VERSION);

  if (vm.count("timer"))
    TIMER=vm["timer"].as<bool>();
  if (vm.count("class"))
    SAVE_CLASS_=vm["class"].as<bool>();
  if (vm.count("prob"))
    SAVE_PROB_=vm["prob"].as<bool>();
  if (vm.count("detailed_output"))
    DETAILED_OUTPUT=vm["detailed_output"].as<bool>();
  if (vm.count("datafile"))
    datafile_S=vm["datafile"].as<vector<string> >();
  datafile=datafile_S[0];
  
  if (vm.count("libdatafile"))
    datafile_L=vm["libdatafile"].as<vector<string> >();

  LIB_TYPE=vector<bool> (datafile_L.size(),true);
  vector<bool> LIB_TYPE_tmp;
  
  if (vm.count("configfile"))
    configfile=vm["configfile"].as<string>();
  if (vm.count("outfile"))
    outfile=vm["outfile"].as<string>();
  if (vm.count("numeach"))
    NUM_EACH=vm["numeach"].as< int>();
  if (vm.count("runfile"))
    runfile=vm["runfile"].as<string>();
  if (vm.count("libtype"))
    LIB_TYPE_tmp=vm["libtype"].as<vector<bool> >();

  for(unsigned int i=0;i<LIB_TYPE_tmp.size();++i)
    if(i<LIB_TYPE.size())
      LIB_TYPE[i]=LIB_TYPE_tmp[i];
  
  if (vm.count("modelrun"))
    MODEL_RUN=vm["modelrun"].as<unsigned int>();

  CONFIG cfg(configfile);
  cfg.set(outfile,"OUTFILE_MAP");
  cfg.set(len,"DATA_LENGTH");
  cfg.set(NUM_EACH,"NUM_EACH");
  cfg.set(VERBOSE_,"VERBOSE");
  cfg.set(EVALUATION_DEPTH,"EVALUATION_DEPTH");
  cfg.set_vector<double>(partition,"PARTITION");
  cfg.set_vector<double>(partition_L,"PARTITION_L");
  cfg.set(DATA_TYPE,"DATA_TYPE");
  cfg.set(DATA_DIR,"DATA_DIR");
  cfg.set(DATA_TYPE_L,"DATA_TYPE_L");
  cfg.set(DATA_DIR_L,"DATA_DIR_L");

  if (vm.count("datalen"))
    len=vm["datalen"].as<unsigned int>();
  if (vm.count("datatype"))
    DATA_TYPE=vm["datatype"].as<string>();
  if (vm.count("datadir"))
    DATA_DIR=vm["datadir"].as<string>();
  if (vm.count("partition"))
    partition=vm["partition"].as<vector<double> >();
  if (vm.count("verbose"))
    VERBOSE_=vm["verbose"].as<bool>();

  if (DATA_DIR=="row")
    DATA_DIR="across";

  if (DATA_DIR=="column")
    DATA_DIR="up";

  if ((DATA_DIR!="up") && (DATA_DIR!="across"))
    err_("UNKNOWN DATA DIRECTION : row or column");

  if ((DATA_TYPE!="continuous") && (DATA_TYPE!="symbolic"))
    err_("UNKNOWN DATA TYPE : continuous or symbolic");

  if(DATA_TYPE_L=="")
    DATA_TYPE_L=DATA_TYPE;
  if(DATA_DIR_L=="")
    DATA_DIR_L=DATA_DIR;
  if(partition_L.empty() && !(partition.empty()))
    partition_L=partition;

//-----------------------------------------------------
//-----------------------------------------------------
  vector< vector<Symbolic_string_> > Svec_L;
  for(unsigned int i=0;i<datafile_L.size();++i)
    if(LIB_TYPE[i])
      {
	data_reader *R_L;
	if (DATA_TYPE_L=="continuous")
	  R_L = new data_reader(datafile_L[i],DATA_DIR_L,
				partition_L,len);
	else
	  R_L = new data_reader(datafile_L[i],DATA_DIR_L,len);

	vector < Symbolic_string_ > Svec__;
	Svec__=R_L->getsymbolic_string_vector();
	missing_alphabet_check(Svec__);
	Svec_L.push_back(Svec__);
      }
    else
      {
	vector<string> pfsafilenames;
	CONFIG pfsanames(datafile_L[i]);
	pfsanames.set_column_vector<string>(pfsafilenames,"LIBFILES");

	/*
	ifstream IN(datafile_L[i].c_str());
	string line;
	while (getline(IN,line))
	  {
	    stringstream ss(line);
	    string tmpstr;
	    while(ss>>tmpstr)
	      pfsafilenames.push_back(tmpstr);
	  }
	IN.close();
	*/
	
	vector<Symbolic_string_> svec_tmp;
	for(vector<string>::iterator itr=pfsafilenames.begin();
	    itr!=pfsafilenames.end();
	    ++itr)
	  for(unsigned int j=0;j<MODEL_RUN;++j)
	    svec_tmp.push_back(read_mc(*itr,"PFSA")
			       .gen_data(len));

	if(runfile=="")
	  runfile=outfile+"__tmp__";
	
	ofstream OUT(runfile.c_str());
	for(vector<Symbolic_string_>::iterator itr_=svec_tmp.begin();
	    itr_!=svec_tmp.end();++itr_)
	  OUT<<*itr_<<endl;
	OUT.close();

	data_reader *R_L1;
	if (DATA_TYPE=="continuous")
	  R_L1 = new data_reader(runfile,DATA_DIR_L,partition,len);
	else
	  R_L1 = new data_reader(runfile,DATA_DIR_L,len);

	vector < Symbolic_string_ > Svec__;
	Svec__=R_L1->getsymbolic_string_vector();
	missing_alphabet_check(Svec__);
	Svec_L.push_back(Svec__);

	//Svec_L.push_back(R_L1->getsymbolic_string_vector());	
	//	Svec_L.push_back(svec_tmp);
      }
//-----------------------------------------------------
//-----------------------------------------------------

  vector < Symbolic_string_ > Svec,Svec1;
  data_reader *R;
  if (DATA_TYPE=="continuous")
    R = new data_reader(datafile,DATA_DIR,partition,len);
  else
    R = new data_reader(datafile,DATA_DIR,len);
  Svec = R->getsymbolic_string_vector();

  missing_alphabet_check(Svec);

  vector<matrix_dbl> H,H1;
  for(unsigned int i=0;i<datafile_L.size();++i)
    if (TIMER)
      {
	timer::auto_cpu_timer t;
	Set_symbolic_string_ DSS(Svec,
				 Svec_L[i],NUM_EACH,
				 EVALUATION_DEPTH);
	H.push_back(DSS.distance_matrix()) ;
	if(DETAILED_OUTPUT)
	  {
	    ofstream outM((outfile+"_"+datafile_L[i]).c_str());
	    outM << H.back();
	    outM.close();
	  }
      }
    else
      {
	Set_symbolic_string_ DSS(Svec,
				 Svec_L[i],NUM_EACH,
				 EVALUATION_DEPTH);
	H.push_back(DSS.distance_matrix()) ;
	if(DETAILED_OUTPUT)
	  {
	    ofstream outM((outfile+"_"+datafile_L[i]).c_str());
	    outM << H.back();
	    outM.close();
	  }
      }


  lib_select_ L(H);
  L.get();

  if (SAVE_PROB_)
    {
      ofstream OUT((outfile+"_prob").c_str());
      OUT << L.prob__();
      OUT.close();
    }

  if (SAVE_CLASS_)
    {
      ofstream OUT((outfile+"_class").c_str());
      OUT << L.class__();
      OUT.close();
    }

  if(VERBOSE_)
    {
      cout << "CLASS" << endl;
      cout << L.class__();
      cout << "CLASS PROB" << endl;
      cout << L.prob__();
    }


  return 0;
}

//-----------------------------------------------------
//-----------------------------------------------------

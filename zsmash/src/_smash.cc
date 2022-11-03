#include <stdlib.h>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/tokenizer.hpp>
#include <boost/token_functions.hpp>
#include <exception>
#include <boost/timer/timer.hpp>


#include "semantic.h"
#include "config.h"


using namespace std;
using namespace boost::program_options;

const string VERSION="\nSMASH v1.31415 \nCopyright Ishanu Chattopadhyay 2017 UChicago";
const string EMPTY_ARG_MESSAGE="Exiting. Type -h or --help for usage";


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


matrix_dbl _smash(int argc, char *argv[], vector<vector<uint>> &data_in)
{
  int NUM_EACH=100;
  unsigned int len=1000,BEG=1,END=150,NAMEWIDTH=3;
  unsigned int COLUMNNUM=0;
  unsigned int EVALUATION_DEPTH=8;
  vector <double> partition;
  bool VERBOSE_=false, ONLY_SAE=false, ONLY_PAST=false;
  int HIST=1;
  bool TIMER=false;
  string configfile="", 
    datafile="",outfile="";
  
  string DATA_DIR="up";
  string DATA_TYPE="continuous";
  unsigned int iO=0,iE=0,jO=0,jE=0;


  options_description infor( "Program information");
  infor.add_options()
    ("help,h", "print help message.")
    ("version,V", "print version number");

  options_description usg( "Usage");
  usg.add_options()
    ("datafile,f",value<string>(), "datafile []")
    ("configfile,c",value<string>(), "config file  []")
    ("datadir,D",value<string>(), "data dir [column]")
    ("datatype,T",value<string>(), "data type [continuous]")
    ("seqlen,L",value<unsigned int>(), "max sequence length [100000000]")
    ("numeach,n",value<unsigned int>(), "number of reruns [100]")
    ("partition,p",value<vector<double>>()->multitoken(), "partition vector []")
    ("timer,t",value< bool >(), "display timer [1 (true)] ");
  
  options_description outputopt( "Output options");
  outputopt.add_options()
    ("outfile,o",value<string>(), "result file [default: H.dst]")
    ("verbose,v",value< bool >(), "verbosity level [default: 1] ");

  options_description desc( "\n### SMASH ### (ishanu chattopadhyay 2017)");
  desc.add(infor).add(usg).add(outputopt);

  positional_options_description p;
  variables_map vm;
  if (argc == 1)
    {
      cout << EMPTY_ARG_MESSAGE << endl;
      return {};
    }
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
      return {};
    }
  if (vm.count("help"))
    {
      cout << desc << endl;
      return {};
    }
  if (vm.count("version"))
    {
      cout << VERSION << endl; 
      return {};
    }

  if (vm.count("timer"))
    TIMER=vm["timer"].as<bool>();
  if (vm.count("datafile"))
    datafile=vm["datafile"].as<string>();
  if (vm.count("configfile"))
    configfile=vm["configfile"].as<string>();
  if (vm.count("outfile"))
    outfile=vm["outfile"].as<string>();
  if (vm.count("datadir"))
    DATA_DIR=vm["datadir"].as<string>();
  if (vm.count("datatype"))
    DATA_TYPE=vm["datatype"].as<string>();
  if (vm.count("seqlen"))
    len=vm["seqlen"].as<unsigned int>();
  if (vm.count("numeach"))
    NUM_EACH=vm["numeach"].as<unsigned int>();
  if (vm.count("partition"))
    partition=vm["partition"].as<vector<double> >();
  if (vm.count("verbose"))
    VERBOSE_=vm["verbose"].as<bool>();

  if(DATA_DIR=="row")
    DATA_DIR="across";
  if(DATA_DIR=="column")
    DATA_DIR="up";

  /** Read Config file*/
  if(configfile!="")
    {
      CONFIG cfg(configfile);
      cfg.set(BEG,"BEG");
      cfg.set(END,"END");
      cfg.set(NAMEWIDTH,"NAMEWIDTH");
      cfg.set(COLUMNNUM,"COLUMNNUM");
      cfg.set(outfile,"OUTFILE_MAP");
      cfg.set(len,"DATA_LENGTH");
      cfg.set(NUM_EACH,"NUM_EACH");
      cfg.set(VERBOSE_,"VERBOSE");
      cfg.set(ONLY_SAE,"ONLY_SAE");
      cfg.set(ONLY_PAST,"ONLY_PAST");
      cfg.set(HIST,"ONLY_PAST_HISTORY");
      cfg.set(EVALUATION_DEPTH,"EVALUATION_DEPTH");
      cfg.set_vector<double>(partition,"PARTITION");
      cfg.set(DATA_TYPE,"DATA_TYPE");
      cfg.set(DATA_DIR,"DATA_DIR");
      cfg.set(iO,"IO");
      cfg.set(iE,"IE");
      cfg.set(jO,"JO");
      cfg.set(jE,"JE");
    }

  if (VERBOSE_)
    {
      cout << "datafile:" << datafile << " outfile: "
	   << outfile << " configfile: "
	   <<  configfile << " BEG: "
	   << BEG<< endl;;
    }
  vector < Symbolic_string_ > Svec;
  /* if (datafile=="") */
  /*   get_continuous_DataMatrix(BEG,END,NAMEWIDTH, */
			      /* COLUMNNUM,partition, */
			      /* len,Svec); //SUEL3 */
  /* else */
    /* { */
      data_reader *R;
      /* if (DATA_TYPE=="continuous") */
	/* R = new data_reader(datafile,DATA_DIR,partition,len); */
      if (data_in.size() > 0)
    R = new data_reader(data_in,DATA_DIR,len);
      else if (datafile != "")
	R = new data_reader(datafile,DATA_DIR,len);
      else
    cout << "ERROR: No data passed to binary." << endl;
      Svec = R->getsymbolic_string_vector();
    

      missing_alphabet_check(Svec);
      
      if (VERBOSE_)
	{
	  cout << "reading " << DATA_TYPE << " data from file.. streams read: "
	       << Svec.size() << endl;
	}
    /* } */

  Set_symbolic_string_ *DSS;
  if (TIMER)
    {
      boost::timer::auto_cpu_timer t;
      if(iO!=iE)
	DSS = new Set_symbolic_string_ (Svec,NUM_EACH,EVALUATION_DEPTH,iO,iE,jO,jE);
      else
	DSS = new Set_symbolic_string_ (Svec,NUM_EACH,EVALUATION_DEPTH);
    }
  else
    {
      if(iO!=iE)
	DSS = new Set_symbolic_string_ (Svec,NUM_EACH,EVALUATION_DEPTH,iO,iE,jO,jE);
      else
	DSS = new Set_symbolic_string_ (Svec,NUM_EACH,EVALUATION_DEPTH);
    }
  matrix_dbl M1= DSS->distance_matrix() ;

  if (outfile != "")
    {
      ofstream outM(outfile.c_str());
      outM << M1;
      outM.close();
    }
  
  return M1;
}


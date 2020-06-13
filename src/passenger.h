#include <cstdio>
#include <getopt.h>
#include <cstdlib>
#include <cstring>

struct parameters {
  char* passenger_f;
  std::string urate;
};

struct parameters* interface(struct parameters* param,int argc, char *argv[]);
void delete_param(struct parameters* param);
void usage(void);

const char* program_name;

struct parameters* interface(struct parameters* param, int argc, char *argv[]){

  program_name = argv[0];
  int c;     // the next argument
  int help = 0;

  if (argc < 2) {
    usage();
    exit(0);
  }

  param = new struct parameters;
  param->passenger_f = new char;
 
  const struct option long_options[] ={
    {"passenger",1,0, 'p'},
    {"urate",1,0,'u'},
    {"help",0,0,'h'},
    {0, 0, 0, 0}
  };

  while (1) {

    int option_index = 0;
    c = getopt_long_only (argc, argv,"hp:u:",long_options, &option_index);

    if (c == -1) {
      break;
    }

    switch(c) {
    case 0:
      break;
    case 'p':
      param->passenger_f = optarg;
      break;
    case 'u':
      param->urate = std::string(optarg);
      break;
    case 'h':
      help = 1;
      break;
    case '?':
      help = 1;
      break;
    default:
      help = 1;
      break;
    }
  }

  if(help) {
    usage();
    delete_param(param);
    exit(0);
  }

  return param;
}

void usage()
{
  fprintf(stdout, "\n, Copyright (C) 2020 Sun, Ruping <ruping@umn.edu>\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "Usage: %s options [cell_id_file.gz] urate \n\n", program_name);
  fprintf(stdout, "-h --help    print the help message\n");
  fprintf(stdout, "-p --passenger  <filename>     passenger mutation file, gzipped\n");
  fprintf(stdout, "-u --urate      float          passenger mutation rate\n"); 
  fprintf(stdout, "\n");
}

void delete_param(struct parameters* param)
{
  delete(param->passenger_f);
  delete(param);
}

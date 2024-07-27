/*=============================================================================
# Filename: indexBuild.cpp
1. ./indexBuild                                        print the help message
2. ./indexBuild --help                                 simplified as -h, equal to 1
3. ./indexBuild db_folder                              load query from given path fro given database
=============================================================================*/

#include "../Database/Database.h"
#include "../Util/Util.h"

using namespace std;

void
help()
{
  printf("# Filename: indexBuild.cpp\n"
         "1. ./indexBuild                                        print the help message\n"
         "2. ./indexBuild --help                                 simplified as -h, equal to 1\n"
         "3. ./indexBuild db_folder                              load database from given path for given database");
}

int main(int argc, char * argv[])
{
  //chdir(dirname(argv[0]));
//#ifdef DEBUG
  Util util;
//#endif

  if (argc == 1 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0)
  {
    help();
    return 0;
  }
  cout << "index building..." << endl;
  if (argc < 2)
  {
    cout << "error: lack of DB_store to be queried" << endl;
    return 0;
  }
  {
    cout << "argc: " << argc << "\t";
    cout << "DB_store:" << argv[1] << "\t";
    cout << endl;
  }

  string db_folder = string(argv[1]);
  int len = db_folder.length();
  if(db_folder.substr(len-3, 3) == ".db")
  {
    cout<<"your database can not end with .db"<<endl;
    return -1;
  }

  Database _db(db_folder);
  _db.load(false,false,false);
  cout<<" finish load "<<endl;
  _db.BuildHopIndex();
  cout << "finish BuildIndex" << endl;
  return 0;
}

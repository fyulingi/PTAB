
#include "../Database/Database.h"
#include "../Util/Util.h"

using namespace std;

const string QUERY_PREFIX = "./Experiments/";
const string RESULT_PREFIX = "./ExperimentResults/Various_k_one_query/";
extern long UTIL_current_query_time;
extern vector<double> test_scores;
extern vector<vector<unsigned int>> result_ids;
extern bool k_value_feed;
extern int benchmark_k;
extern bool any_k;
extern int any_k_step;
extern vector<long> any_k_time;
extern std::vector<int> fq_num_vec;
extern std::vector<int> fr_num_vec;
extern std::vector<int> ow_num_vec;
extern bool iterator_count_enable;
extern int FQ_NUM;
extern int FR_NUM;
extern int OW_NUM;
extern RootStrategy global_root_strategy;
extern string root_var_name;

string TopKStrategyString(TopKStrategy method) {
    if (method == TopKStrategy::DP_B)
        return "DP_B";
    if (method == TopKStrategy::DP_P)
        return "DP_P";
    if (method == TopKStrategy::K_TPM)
        return "kTPM";
    if (method == TopKStrategy::DFS_OLD)
        return "DFS_OLD";
    if (method == TopKStrategy::PTAB)
        return "PTAB";
    if (method == TopKStrategy::TAKE2)
        return "Take2";
    if (method == TopKStrategy::TAKE_ALL)
        return "TAKE_ALL";
    if (method == TopKStrategy::EAGER)
        return "Eager";
    if (method == TopKStrategy::SAE)
        return "SAE";
    if (method == TopKStrategy::RankAfterMatching)
        return "Naive";
    return "Unknown";
}

vector<string> db_names{string("bsbm"), string("DBpedia3.5.1"), string("Yago")};
vector<string> query_folders{string(QUERY_PREFIX + "bsbmQueries"),
                             string(QUERY_PREFIX + "DBpedia3.5.1Queries"),
                             string(QUERY_PREFIX + "YagoQueries")};

vector<TopKStrategy> method_vec{TopKStrategy::PTAB,
                                TopKStrategy::TAKE2,
                                TopKStrategy::K_TPM,
                                TopKStrategy::EAGER,
                                TopKStrategy::DP_B};
//                                TopKStrategy::RankAfterMatching};

vector<vector<string>> query_names{{"1", "2", "3", "4", "5", "6"},
                                   {"1", "2", "3", "4", "5", "6"},
                                   {"1", "2", "3", "4", "5", "6"}};
constexpr int avg_times = 3;

void TestQueryWithK(Database &_db, const string &db_name, const string &query_folder, const string &query_name,
                    ofstream &result_file_stream, int step, int start, int end, bool warm_up) {
  string query_file = query_folder + "/q" + query_name + ".sql";
  auto query = Util::getQueryFromFile(query_file.c_str());
  if (query.empty()) {
    cout << "query file: " << query_file << " not exists." << endl;
    return;
  }

  ResultSet _rs;
  FILE *ofp = stdout;
  string msg;

  any_k_step = 1;

  for (auto method: method_vec) {
    string method_name = TopKStrategyString(method);
    cout << "method name: " << method_name << endl;

    // warm up
    any_k = false;
    k_value_feed = false;
    int ret = _db.query(query, _rs, ofp, true, false, nullptr, method);

    k_value_feed = true;
    any_k = true;

    benchmark_k = end;
    any_k_time.reserve(end);

    vector<double> times_add(end);
    for (unsigned int t = 0; t < avg_times; ++t) {
      any_k_time.clear();
      ret = _db.query(query, _rs, ofp, true, false, nullptr, method);
      for (unsigned index = start-1; index < end; index += step) times_add[index] += any_k_time[index];
    }

//    any_k_time.clear();
//    ret = _db.query(query, _rs, ofp, true, false, nullptr, method);
    for (unsigned index = start-1; index < end; index += step) result_file_stream << times_add[index]/(avg_times * 1.0) << (index == end-1 ? "\n" : ",");
//    for (auto t : any_k_time) result_file_stream << " " << t;
  }
}

int main(int argc, char *argv[]) {
  k_value_feed = false;
  Util util;
  cout << "Test Various K Value..." << endl;
  if (!util.dir_exist(RESULT_PREFIX))
    util.create_dir(RESULT_PREFIX);

  int start = 1;
  int end = 101;
  int step = 20;

  if (argc == 4) {
    start = atoi(argv[1]);
    end = atoi(argv[2]);
    step = atoi(argv[3]);
  }

  ofstream result_file_stream(RESULT_PREFIX + "various_k.txt");

  for (int i = 1; i < db_names.size(); ++i) {
    string db_name = db_names[i];
    cout << "DB name: " << db_name << endl;
    result_file_stream << db_name << endl;
    Database _db(db_name);
    _db.load();
    for (unsigned int q_i = 0; q_i < query_names[i].size(); q_i++) {
      cout << "query: " << query_folders[i] << "  " << query_names[i][q_i] << "...." << endl;
      TestQueryWithK(_db, db_names[i], query_folders[i], query_names[i][q_i], result_file_stream, step, start, end, true);
      result_file_stream << endl;
    }
    result_file_stream << endl;
  }
  return 0;
}

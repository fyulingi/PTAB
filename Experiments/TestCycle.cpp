
#include "../Database/Database.h"
#include "../Util/Util.h"

using namespace std;

const string QUERY_PREFIX = "./Experiments/";
const string RESULT_PREFIX = "./ExperimentResults/Cycle_time/";
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

/**
 * ref from: https://stackoverflow.com/questions/2342162/stdstring-formatting-like-sprintf
 */

template<typename ... Args>
std::string string_format(const std::string &format, Args ... args) {
  int size_s = std::snprintf(nullptr, 0, format.c_str(), args ...) + 1; // Extra space for '\0'
  if (size_s <= 0) { throw std::runtime_error("Error during formatting."); }
  auto size = static_cast<size_t>( size_s );
  auto buf = std::unique_ptr<char[]>(new char[size]);
  std::snprintf(buf.get(), size, format.c_str(), args ...);
  return std::string(buf.get(), buf.get() + size - 1); // We don't want the '\0' inside
}

string TopKStrategyString(TopKStrategy method) {
  if (method == TopKStrategy::DP_B)
    return "DP_B";
  if (method == TopKStrategy::DP_P)
    return "DP_P";
  if (method == TopKStrategy::K_TPM)
    return "K_TPM";
  if (method == TopKStrategy::DFS_OLD)
    return "DFS_OLD";
  if (method == TopKStrategy::PTAB)
    return "PTAB";
  if (method == TopKStrategy::TAKE2)
    return "TAKE2";
  if (method == TopKStrategy::TAKE_ALL)
    return "TAKE_ALL";
  if (method == TopKStrategy::EAGER)
    return "EAGER";
  if (method == TopKStrategy::SAE)
    return "SAE";
  if (method == TopKStrategy::RankAfterMatching)
    return "Naive";
  return "Unknown";
}

string RootStrategyString(RootStrategy method) {
  if (method == RootStrategy::Random)
    return "Random";
  if (method == RootStrategy::MinCand)
    return "MinCand";
  if (method == RootStrategy::MinDepth)
    return "MinDepth";
  if (method == RootStrategy::Ours)
    return "Ours";
  return "Unknown";
}

vector<string> db_names{string("DBpedia3.5.1")};
vector<string> query_folders{string(QUERY_PREFIX + "cycle_query")};

vector<RootStrategy> StrategyVec{
    RootStrategy::Random,
    RootStrategy::MinDepth,
    RootStrategy::MinCand,
    RootStrategy::Ours};

vector<vector<string>> query_names{{"1", "2", "3", "4", "5"}};

vector<TopKStrategy> cyclic_method_vec{TopKStrategy::PTAB,
                                       TopKStrategy::RankAfterMatching};

constexpr int avg_times = 3;

void TestQuery(Database &_db, ofstream &result_file_stream, const string &db_name, const string &query_folder, const string &query_name, bool warm_up) {
  result_file_stream << "q" + query_name << endl;
  string query_file = query_folder + "/q" + query_name + ".sql";
  auto query = Util::getQueryFromFile(query_file.c_str());
  if (query.empty()) {
    cout << "query file: " << query_file << " not exists." << endl;
    return;
  }

  ResultSet _rs;
  FILE *ofp = stdout;
  string msg;

  for (unsigned index = 0; index < cyclic_method_vec.size(); ++index) {
    result_file_stream << TopKStrategyString(cyclic_method_vec[index]) << (index == cyclic_method_vec.size()-1 ? "\n" : ",");
  }
  if (warm_up) { // warm up
    for (auto method: cyclic_method_vec) {
      cout << TopKStrategyString(method) << endl;
      int ret = _db.query(query, _rs, ofp, true, false, nullptr, method);
    }
  }

  for (unsigned index = 0; index < cyclic_method_vec.size(); ++index) {
    auto method = cyclic_method_vec[index];
    string method_name = TopKStrategyString(method);
    cout << "method name: " << method_name << endl;
    long total_time = 0;
    for (unsigned int t = 0; t < avg_times; ++t) {
      int ret = _db.query(query, _rs, ofp, true, false, nullptr, method);
      total_time += UTIL_current_query_time;
    }
    result_file_stream << string_format((index == cyclic_method_vec.size()-1 ? "%lf\n" : "%lf,"), (total_time * 1.0) / (avg_times * 1.0));
  }
}

int main(int argc, char *argv[]) {
  k_value_feed = false;
  Util util;
  cout << "Test Cycle..." << endl;
  if (!util.dir_exist(RESULT_PREFIX))
    util.create_dir(RESULT_PREFIX);

  ofstream result_file_stream(RESULT_PREFIX + "cycle_time.txt");
  for (int i = 0; i < db_names.size(); ++i) {
    string db_name = db_names[i];
    cout << "DB name: " << db_name << endl;
    result_file_stream << db_name << endl;
    Database _db(db_name);
    _db.load();
    for (unsigned int q_i = 0; q_i < query_names[i].size(); q_i++) {
      cout << "query: " << query_folders[i] << "  " << query_names[i][q_i] << "...." << endl;
      TestQuery(_db, result_file_stream, db_names[i], query_folders[i], query_names[i][q_i], true);
    }
  }

  return 0;
}
